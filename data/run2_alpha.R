
rm(list = ls())
load("dataPheno")

#Only the 666 "HIP" samples are used

sampleID <- sampleID[grep("HIP", sampleID)]
nSample <- length(sampleID)
dataPheno <- dataPheno[sampleID, ]

##################################################################

#Function that run the rarefaction and calculate the two alpha diversities (Shannon index and observed unique sequence)

funRarefaction <- function(x, depth = seq(1e4, 1e5, 1e4), nSimu = 10){
  nX <- length(x)
  nD <- length(depth)
  s <- sum(x)
  tempVec <- rep(1:nX, x)
  alphaMat1 <- alphaMat2 <- matrix(NA, nSimu, nD, dimnames = list(paste0("Simu", 1:nSimu), depth))
  for(i in 1:nD){
    tempDepth <- depth[i]
    if(tempDepth > s)next
    for(j in 1:nSimu){
      # rarefraction
      tempSample <- sample(tempVec, tempDepth)
      tempTable <- table(tempSample)
      if(sum(tempTable) != tempDepth)stop("something wrong!")
      # alpha diversity calculation
      tempP <- tempTable / tempDepth
      alphaMat1[j, i] <- -sum(tempP * log(tempP)) # Shannon index
      alphaMat2[j, i] <- length(tempP)            # Observed unique sequence
      cat(i, j, alphaMat1[j, i], alphaMat2[j, i], "\n")
    }
  }
  #	metrics = c("shannon", "observed")
  return(list(Shannon = alphaMat1, Observed = alphaMat2))
}

##################################################################

#Go over the ".RData" files and extract the number of each unique sequence and save them to "readsList".

fileNames <- dir("TCR_RData", full.names = TRUE); fileNames <- fileNames[grep("P", fileNames)]
readsList <- list(); length(readsList) <- nSample; names(readsList) <- sampleID
for(i in 1:length(fileNames)){
  load(fileNames[i])
  readsList[[i]] <- tempDataRaw$reads
  cat(i, fileNames[i], "\n")
}
#str(readsList)
#save(readsList, file = "readsList")
load("readsList")

totalReads <- unlist(lapply(readsList, sum))
#hist(totalReads, seq(0, 5e7, 5e5))

##################################################################

#Run rarefaction at different levels, calculate the alpha diversities, save to "dataAlpha".
#(It may take hours. To reduce time: 1. decrease nSimu; 2. shorten the number of depth levels; 3. try do it parallelly for each file.)

depth <- c(seq(1e5, 1e6, 1e5), seq(2e6, 1e7, 1e6))
nDepth <- length(depth)
nSimu <- 10

dataAlpha1 <- dataAlpha2 <- matrix(NA, nSample, nDepth, dimnames = list(sampleID, depth))
for(i in 1:length(fileNames)){
  load(fileNames[i])
  set.seed(i) #for replicability
  tempAlphaList <- funRarefaction(tempDataRaw$reads, depth = depth, nSimu = nSimu)
  tempList <- lapply(tempAlphaList, colMeans)
  dataAlpha1[i, ] <- tempList[[1]]
  dataAlpha2[i, ] <- tempList[[2]]
  cat(i, fileNames[i], dataAlpha1[i, ], dataAlpha2[i, ], "\n")
}

#save(list = c("totalReads", "dataAlpha1", "dataAlpha2"), file = "dataAlpha")
# load("dataAlpha")
