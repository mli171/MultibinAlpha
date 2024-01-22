library(readr)

rm(list = ls())

##################################################################

# Process each ".tsv" file, keep the necessary columns, save each one to a ".tsv.gz" file.

inputfileNames = dir("/Volumes/Work/TCR_data/emerson-2017-natgen", full.names = TRUE)
outputfileNames = dir("/Volumes/Work/TCR_data/emerson-2017-natgen", full.names = FALSE)
outputfileNames = unlist(strsplit(outputfileNames, split = ".tsv"))

# selected variables from original data file
VarNames_orig = c("sample_name", "species", "locus", "product_subtype", 
                  "sample_tags", "rearrangement", "amino_acid", "frame_type", 
                  "rearrangement_type", "seq_reads", "frequency", 
                  "productive_frequency", "v_family", "v_gene", 
                  "d_family", "d_gene", "j_family", "j_gene", 
                  "v_family_ties", "v_gene_ties", "d_family_ties", 
                  "d_gene_ties", "j_family_ties", "j_gene_ties")
# actual name used in later files (to match current data files)
tempVarNames <- c("sample_name", "species", "locus", "product_subtype", 
                  "sample_tags", "rearrangement", "amino_acid", "frame_type", 
                  "rearrangement_type", "reads", "frequency", 
                  "productive_frequency", "v_family", "v_gene", 
                  "d_family", "d_gene", "j_family", "j_gene", 
                  "v_family_ties", "v_gene_ties", "d_family_ties", 
                  "d_gene_ties", "j_family_ties", "j_gene_ties")

for(iii in 1:length(inputfileNames)){
  tmpdata = read.table(inputfileNames[iii], header = TRUE, sep = "\t", 
                       colClasses = "character")
  tmptargetdata = tmpdata[, VarNames_orig]
  colnames(tmptargetdata) = tempVarNames
  # save to ".tsv" file and zip the file for storage save
  write_tsv(tmptargetdata, file.path(dir="/Volumes/Work/TCR_data/TCR_zipfiles/", 
                                     paste0(outputfileNames[iii],".tsv.gz")))
  cat("Sample:", outputfileNames[iii], "is completed!\n")
}

##################################################################

# ### check the first sample
# fileNames = dir("/Volumes/Work/TCR_data/TCR_zipfiles/", full.names = TRUE)
# tempDataRaw = read.table(fileNames[1], header = TRUE, sep = "\t", colClasses = "character")
# dim(tempDataRaw)
