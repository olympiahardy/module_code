
Here we analysed RNA-sequencing read data from two Drosophila tissue samples. AssignmentScript.sh is the script submitted 
to the HPC completed in non-discovery mode where the reference transcriptome was a GTF file for the Drosophila genome. 

The file FinalLoopDiscov.sh is a similar script however this time we allowed for discovery by creating our own 
sample-specific transcriptome as a reference. 

The file AssignmentRScript.R is taking the generated .csv files from the non-discovery script output and carrying out 
DE analysis using DESeq2 and creating an MA plot to observe LFC shrinkage between the two tissue samples.

PP_script.sh is a script submitted to the HPC when we did our Pathogen Polyomics module investigating the genetic basis
of amphotericin B resistant Leismania mexicana.
