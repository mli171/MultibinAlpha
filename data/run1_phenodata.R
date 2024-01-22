
rm(list = ls())

##################################################################

#Process each ".gz" file, keep the necessary columns, save each one to a ".RData" file.

fileNames <- dir("TCR_zipfiles", full.names = TRUE)
fileNames <- fileNames[grep("gz", fileNames)]
sampleID <- unlist(lapply(strsplit(fileNames, "/"), function(x)x[length(x)]))
sampleID <- unlist(lapply(strsplit(sampleID, "[.]"), function(x)x[1]))
nSample <- length(sampleID)

tempVarNames <- c("rearrangement", "amino_acid", "frame_type",
                  "rearrangement_type", "reads", "frequency",
                  "productive_frequency", "v_family", "v_gene",
                  "d_family", "d_gene", "j_family", "j_gene",
                  "v_family_ties", "v_gene_ties", "d_family_ties",
                  "d_gene_ties", "j_family_ties", "j_gene_ties")

for(i in 1:nSample){
  tempDataRaw <- read.table(fileNames[i], header = TRUE, sep = "\t", colClasses = "character")
  tempDataRaw <- tempDataRaw[, tempVarNames]
  tempDataRaw$reads <- as.integer(tempDataRaw$reads)
  tempDataRaw$frequency <- as.numeric(tempDataRaw$frequency)
  tempDataRaw$productive_frequency <- as.numeric(tempDataRaw$productive_frequency)
  for(j in 1:ncol(tempDataRaw))if(is.character(tempDataRaw[, j]))tempDataRaw[tempDataRaw[, j] == "", j] <- NA
  save(tempDataRaw, file = paste0("TCR_RData/", sampleID[i], ".RData"))
  cat(i, sampleID[i], "\n")
}




##################################################################





##################################################################

#Read the first row of each file, extract the phenotype information, then generate "dataPheno".

tagVec <- rep(NA, nSample); names(tagVec) <- sampleID
for(i in 1:nSample){
  tempDataRaw <- read.table(fileNames[i], header = TRUE, sep = "\t", colClasses = "character", nrow = 1)
  tagVec[i] <- tempDataRaw$sample_tags
  cat(i, sampleID[i], "\n")
}

tagList <- list(); length(tagList) <- nSample; names(tagList) <- sampleID
for(i in 1:nSample){
  tempTag <- tagVec[i]
  tempVec <- unlist(strsplit(tempTag, ","))
  tempStrsplit <- strsplit(tempVec, ":")
  if(any(unlist(lapply(tempStrsplit, length)) != 2))stop("Something wrong!")
  tempVec <- unlist(lapply(tempStrsplit, function(x)x[2]))
  names(tempVec) <- unlist(lapply(tempStrsplit, function(x)x[1]))
  tagList[[i]] <- tempVec
}
#str(tagList)

tempNames <- names(table(unlist(lapply(tagList, names))))
tagMat <- matrix(NA, nSample, length(tempNames), dimnames = list(sampleID, tempNames))
for(i in 1:nSample){
  for(j in 1:length(tempNames)){
    tempChar <- tagList[[i]]
    tempChar <- tempChar[names(tempChar) == tempNames[j]]
    if(length(tempChar) == 0) next
    if(length(tempChar) > 1) tempChar <- paste(tempChar, collapse = ",")
    tagMat[i, j] <- tempChar
  }
}

dataPheno <- data.frame(tagMat, stringsAsFactors = FALSE, check.names = FALSE)
dataPheno$Age <- as.integer(gsub(" Years", "", dataPheno$Age))
#table(dataPheno$'Inferred CMV status', dataPheno$'Inferred CMV status (cross-validation)', useNA = "a")
#table(dataPheno$'Inferred CMV status', dataPheno$'Inferred CMV status (trained on Cohort 1)', useNA = "a")
#table(dataPheno$'Inferred CMV status (trained on Cohort 1)', dataPheno$Cohort, useNA = "a")
dataPheno$'Inferred CMV status'[dataPheno$Cohort == "Cohort 02"] <- dataPheno$'Inferred CMV status (trained on Cohort 1)'[dataPheno$Cohort == "Cohort 02"]
#table(dataPheno$'Inferred CMV status', useNA = "a")
#str(dataPheno)

#save(list = c("sampleID", "nSample", "dataPheno"), file = "dataPheno")
# load("dataPheno")

##################################################################




