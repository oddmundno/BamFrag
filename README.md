*****************************************************

BAMFRAG script README

Script that extracts fragment lengths, genomic coordinates and end motifs from BAM files

*****************************************************

The script was constructed for the publication Lapin et al. 2025...... The text below adds some additional method information that explains the usage of the script in the publication, for those that would like to reproduce on their own data. 


BIOINFORMATIC PROCESSING OF cfDNA FRAGMENTATION DATA

Fragmentation analysis was performed on cfDNA sequencing data from a previous study of the same patients as described above (1, 2). cfDNA sequencing reads that were not completely sequenced were first removed from aligned BAM files using sambamba version 0.8.0 (call ’sambamba view -F “[ZA]!=null” -f bam <infile.bam> <outfile.bam>). Adapter sequences were removed and multiple reads of the same cfDNA collapsed to consensus sequences using locally developed python scripts TaqXtractor and SSCScreator as previously described (2) with a minimum read number per consensus sequence of one (call ’HYTEC-pipeline.sh -i <infile.bam> -o <outfile.bam> -t ’–endtrimming 0’ -c ’–rep_filt 5 –minmem 1’ ”). Offtarget sequencing reads were extracted using samtools version 1.8 with the call “samtools view -hb -L <offtarget bedfile>”. The length frequencies and 3’ end motif sequences of consensus cfDNA molecules were extracted from aligned BAM files by locally developed python script BamFrag version 2.1 (call ’BamFrag2.py –infile <infile.bam> –outfile <outfile.bam> –short-cutoff-low 112 –short-cutoff-high 155 –long-cutoff-low 170 –long-cutoff-high 181 –statsfile <fragment statistics text file> –fragmentfrequencyfile <fragment frequency text file> –endmotiffile <end motif frequency text file>. Relative fragment length frequencies and 4-base end motifs were calculated by dividing the read numbers on the total number of off-target reads from each sample. 


CLUSTERING AND HEATMAP ANALYSIS

Fragment length and end motif data in the samples were first normalized against the mean frequencies in the samples from healthy volunteers and transformed by log2 transformation (denoted log2 ratio vs healthy). Patient and control cfDNA samples were clustered according to DMRs, fragment length frequencies between 100 and 200 base pairs of length, and 4-base 5’ end motifs by unsupervised hierarchical clustering within-group (patient/control), using the hclust function in R with a eucledian distance measure. A complete linkage agglomeration method was used for methylation data, and ward.D method was used for fragment length  and end motif frequencies. The order of DMR columns and end motif columns in the heatmaps were also computed by hierarchical clustering, using the same methods. 

MACHINE LEARNING-BASED MODELING OF ctDNA CONCENTRATION
Machine learning-based regression analysis was performed using the caret package for R, version 6.0-94. DMR, fragment length frequencies (between 50 and 200 bp) and end motif frequencies from 35 baseline patient samples (one with known high levels of ctDNA, but no point mutation in the employed sequencing panel was excluded) and 10 randomly selected healthy volunteers were centered and scaled using the “center” and “scale” methods of the preprocess function in the caret package. Known ctDNA concentrations, based on previous mutation-based measurements (1, 2), were modeled by general linear regression using the glmnet wrapper in caret, with the leave-one-out resampling and cross-validation method. Glmnet features elastic net (combining lasso and ridge regulization) built-in regularization. An example trainControl command was ’trcontroller = trainControl(method=”loocv”,number=10, repeats=3,savePredictions=’final’)’ and a training command: ’model = train(ctDNA ~ .,data=feature.dataframe,trControl=trcontroller,method=’glmnet’)’. Other machine learning models were also used in the early stages of the analysis but turned out to produce similar results. Separate models for the methylation, fragment length and end motif data were first built and then two ensemble models, either combining the three models in a linear combined model or entering all explanatory variables in a single training. The methylation model was strongly dominating both ensemble models. All models were validated in a separate sample set obtained during follow-up of the same patients and known to contain ctDNA by mutation analysis. Scatterplot of estimated and measured ctDNA concentrations in training and validations data sets are shown in Supplemental Figure 1. Root mean square error (RMSE) in training and validations sets are provided in supplemental Figure S2. 


REFERENCES

1.	Tjensvoll K, Lapin M, Gilje B, Garresori H, Oltedal S, Forthun RB, et al. Novel hybridization- and tag-based error-corrected method for sensitive ctDNA mutation detection using ion semiconductor sequencing. Sci Rep. 2022;12(1):5816.
2.	Lapin M, Edland KH, Tjensvoll K, Oltedal S, Austdal M, Garresori H, et al. Comprehensive ctDNA Measurements Improve Prediction of Clinical Outcomes and Enable Dynamic Tracking of Disease Progression in Advanced Pancreatic Cancer. Clin Cancer Res. 2023;29(7):1267-78.

