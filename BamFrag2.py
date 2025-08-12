#!/usr/bin/python3 -tt

#Version 0.1 first production version
#Version 0.2 added support for maxlen, ShortFr is now computed by division on the number of  reads with length between size_cutoff_high and maxlen                                              
#Version 0.3 with four different cutoffs                    
# Version 0.4. new possibility to export only the n first bases of each read
# Version 0.5 computes also fraction of fragments with length below 150, ala Pedener
# Version 2, exports instead the absolute frequencies (counts) of fragment lengths between the outer cutoffs. Skipped the --endmotif option
# Version 2.1, exports also endmotifs counts to a dedicated file


import pysam
import sys
from argparse import ArgumentParser
import numpy
from scipy import stats
from itertools import product


########################################################
# MAIN FUNCTION
########################################################


def main():
    ########################
    #Parameters to be input and variables to initatilize
    ########################
    parser=ArgumentParser(description="Software to extract fragment lengths and genomic coordinates. Can also filter BAM file according to size cutoff.")
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file. Stdin if -.", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="output BAM file. Stdin if -.", required=False)
    parser.add_argument("--fragmentfrequencyfile", action="store", dest="fragmentfile", help="output file with relative fragment length frequencies", required=False)
    parser.add_argument("--statsfile", action="store", dest="statsfile", help="output descriptive statistics for fragment size filtration.", required=False)
    parser.add_argument("--short-cutoff-high", action="store", type=int,default=158,dest="shortcutoff_high", help="Cutoff for selection of reads for output BAM file. Only shorter fragments are kept. Both high and low cutoff need to be specified to work.", required=False)
    parser.add_argument("--short-cutoff-low", action="store", type=int,default=113,dest="shortcutoff_low", help="Cutoff for selection of reads for output BAM file. Only longer or identical length fragments are kept. Both high and low cutoff need to be specified to work. 0 means no low cutoff.", required=False)
    parser.add_argument("--long-cutoff-high", action="store", type=int,default=189,dest="longcutoff_high", help="Cutoff for the max length of reads to be used in computing ShortFr.", required=False)
    parser.add_argument("--long-cutoff-low", action="store", type=int,default=171,dest="longcutoff_low", help="Cutoff for the max length of reads to be used in computing ShortFr.", required=False)
    parser.add_argument("--endmotiffile", action="store", dest="endmotiffilename", help="Output counts of all possible 4-base endmotifs of length n to this file", required=False)


    
    o = parser.parse_args()

    #Require that both high and low cutoff are specified
#    if (o.cutoff_high and  not o.cutoff_low) or (not o.cutoff_high and o.cutoff.low):
 #       print("Both high and low fragment length cutoff need to be specified")
  #      sys.exit()

    # OPEN INPUT AND OUTPUT FILES
    # If the user want to read from stdin (Unix piping) instead, use -
    # Pysam handles that automatically
    
    inBam = pysam.AlignmentFile( o.infile, "rb" ,threads=2) # Open the input BAM file

    if o.outfile:
        if o.outfile == '-':
            outBam = sys.stdout
        else:
            outBam = pysam.AlignmentFile(o.outfile,"wb",template=inBam)

  
        #If the user want to write to stdout instead of a file
  
          

    if o.statsfile:
        statsfile=open(o.statsfile,'w')  # Output text file for fragment info
    
    totalcounter = 0 # To count the totalt number of seqs
    shorter150counter = 0 # To count reads shorter than 150 bases
    shortcounter = 0  # To count the number of seqs shorter than cutoff
    longcounter = 0 #To count the number of long reads, between size_cutoff_high and maxlen
    sizes = []   # List of fragment sizes for statistics

    if o.endmotiffilename:
        motifs = [''.join(m) for m in product('ACTG',repeat=4)] #itertools.product
        #generate all possible endmotifs
        motifdict = {m:0 for m in motifs}  # Generate empty dict for counting motifs

    # PROCESS INPUT FILE
    
    seqfetcher = inBam.fetch( until_eof = True )   #Open inBam for reading

    for seq in seqfetcher:   # Loop through the sequences
        totalcounter +=1
        fragmentlength = seq.query_length    #Reads the actual length of the read (adapters has been removed previously)
        sizes.append(fragmentlength)
        
        if (fragmentlength<150):
            shorter150counter +=1

        if (fragmentlength>=o.shortcutoff_low) and (fragmentlength<=o.shortcutoff_high):  #If fragment length is between the cutoffs, inclusive
            shortcounter +=1
            if o.outfile:
                outBam.write(seq)           #Write output fragments to bam, if specified
        elif (fragmentlength>=o.longcutoff_low) and (fragmentlength<=o.longcutoff_high):
            longcounter +=1

        startpos = seq.reference_start  #Fetch the start and stop positions in the 
        endpos = seq.reference_end      #reference genome, of the aligned part
        chrom = seq.reference_name        #Get the chromosome number

        if o.endmotiffilename:
            motif = seq.query_sequence[0:4]  # Extract the endmotif from the 5' end
            if motif in motifdict: # Don't count motifs with Ns
                motifdict[motif] += 1  # Count the endmotif
        
        # End for loop    

    #COMPUTE AND OUTPUT STATISTICS
    if o.statsfile:
            meansize = numpy.mean(sizes)
            mediansize = numpy.median(sizes)
            modesize = stats.mode(sizes)[0][0]
            minsize = numpy.min(sizes)
            maxsize = numpy.max(sizes)
            shortfraction = shortcounter/longcounter
            shorter150fraction=shorter150counter/totalcounter
            statsfile.write("Number\tMean\tMedian\tMode\tMin\tMax\tShortFr\tShorter150\n")
            statsfile.write(str(totalcounter)+"\t")
            statsfile.write(str.format("{0:.1f}",meansize)+ "\t")
            statsfile.write(str(mediansize)+ "\t")
            statsfile.write(str(modesize)+ "\t")
            statsfile.write(str(minsize) + "\t")
            statsfile.write(str(maxsize)+ "\t")
            statsfile.write(str.format("{0:.3f}",shortfraction)+ "\t")
            statsfile.write(str.format("{0:.3f}",shorter150fraction)+ "\n")
            statsfile.close()

    if o.fragmentfile:
            fragmentfile = open(o.fragmentfile,'w')  # Output text file for fragment frequencies
            sizelist = [l for l in range(50,251)]
            fragmentfrequencies = [sizes.count(l) for l in sizelist]
            
            for l in sizelist:  #Print header line with 
                fragmentfile.write(str(l)+"\t")
            fragmentfile.write("\n")

            for fr in fragmentfrequencies:
                fragmentfile.write(str(fr) + "\t")  #All frequencies beetween the limits
            fragmentfile.write("\n")
            fragmentfile.close()

    if o.endmotiffilename:
         motiffile = open(o.endmotiffilename,'w')  # Output text file for fragment frequencies
         for k in motifdict.keys():
             motiffile.write(str(k)+'\t')
         motiffile.write('\n')    
         for v in motifdict.values():
             motiffile.write(str(v)+'\t')
         motiffile.write('\n')
         motiffile.close()

            
                                       
    #CLOSE OTHER FILES
    inBam.close()
    if o.outfile:
          outBam.close()
    if o.fragmentfile:
        fragmentfile.close()

if __name__ == "__main__":
    main()
