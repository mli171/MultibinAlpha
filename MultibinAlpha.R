rm(list = ls())

load("TCR_rarefy/readsList")
load("TCR_rarefy/dataPheno")
load("TCR_rarefy/dataAlpha")

sampleID <- sampleID[grep("P", sampleID)]
nSample <- length(sampleID)
dataPheno <- dataPheno[sampleID, ]

totalReads <- unlist(lapply(readsList, sum))

summary(dataPheno[,"Age"])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    1.00   31.00   41.00   39.54   49.00   74.00     111 
summary(as.factor(dataPheno[,"Biological Sex"]))
# Female   Male   NA's 
#    297    345     24 
summary(as.factor(dataPheno[,"inferred CMV status"]))
# Inferred CMV - Inferred CMV + Inferred CMV+-           NA's 
#            364            275              1             26 
summary(as.factor(dataPheno[,"Virus Diseases"]))
# Cytomegalovirus - Cytomegalovirus +              NA's 
#               352               289                25 

pdf("FigureS1.pdf", width = 20, height = 10)
par(mfrow = c(1, 2), cex.main = 3, cex.lab = 2, cex.axis = 2, mar = c(6, 6, 6, 1))
hist(totalReads[totalReads < 2e7] / 1e6, seq(0, 20, 0.5), col = "lightblue", 
     main = "Distribution of Total Reads", xlab = expression(paste("Total Reads (", ''%*%10^6, ")")))
hist(log10(totalReads), seq(4, 8, 0.1), col = "lightblue", 
     main = "Distribution of Total Reads (log scale)", xlab = "log (Total Reads)")
dev.off()



# Overall rarefying with various rarefaction levels ranging from 1e5 to 1e7

depth <- c(seq(1e5, 1e6, 1e5), seq(2e6, 1e7, 1e6))

pdf("FigureS3.pdf", width = 10, height = 15)
par(mfcol = c(3, 2), cex.main = 3, cex.axis = 2, cex.lab = 2, mar = c(4.6, 5.3, 4.6, 1))

plot(NA, xaxt = "n", xlim = c(0, 10), ylim = c(0, 6), type = "l", 
     xlab = expression(paste("Total Reads (", ''%*%10^6, ")")), 
     ylab = expression(paste("The number of unique sequences (", ''%*%10^5, ")")))
axis(1, depth[c(1, 10:19)] / 1e6, depth[c(1, 10:19)] / 1e6)
title("A",  adj = 0, line = 0.8, cex.main = 1.8)
for(i in 1:nrow(dataAlpha2))lines(depth / 1e6, dataAlpha2[i, ] / 1e5)

boxplot(dataAlpha2 / 1e5, at = depth / 1e5, outline = FALSE, ylim = c(0, 6),
        xlab = expression(paste("Total Reads (", ''%*%10^6, ")")), ylab = expression(paste("The number of unique sequences (", ''%*%10^5, ")")), names = NA)
title("B",  adj = 0, line = 1.4, cex.main = 2)
axis(1, depth[c(1, 10:19)] / 1e5, depth[c(1, 10:19)] / 1e6)
mtext(colSums(!is.na(dataAlpha2[, c(1, 10:19)])), 3, at = depth[c(1, 10:19)] / 1e5, font = 2)
mtext("Number of Samples Available", 3, line = 1.3, font = 4)
lines(lowess(rep(depth / 1e5, each = nSample)[!is.na(dataAlpha2)], dataAlpha2[!is.na(dataAlpha2)] / 1e5), col = rgb(1, 0, 0, 0.5), lwd = 4)

tempX <- totalReads
tempY <- dataAlpha2[, "rarefaction_1e+06"]
tempIdx0 <- which(totalReads < 1e7 & !is.na(tempY))
tempLowess0 <- lowess(tempX[tempIdx0], tempY[tempIdx0])
tempTest0 <- cor.test(tempX[tempIdx0], tempY[tempIdx0], method = "spearman")

plot(NA, xaxt = "n", yaxt = "n", 
     xlim = c(0, 10)*1e6, ylim = c(0, 6)*1e5,
     xlab = expression(paste("Total Reads (", ''%*%10^6, ")")),
     ylab = expression(paste("The number of unique sequences (", ''%*%10^5, ")")))
title("C",  adj = 0, line = 0.8, cex.main = 1.8)
points(x=tempX[tempIdx0], y=tempY[tempIdx0])
axis(1, depth[10:19], depth[10:19]/1e6)
axis(2, (0:6)*1e5, 0:6)
lines(tempLowess0, col = "red", lwd = 2)
tempPVec <- tempTest0$p.value
legend("bottomright", paste0("P = ", signif(tempPVec, 3)), cex = 1.5, bty = "n", text.col = "red")

plot(NA, xaxt = "n", yaxt = "n", 
     xlim = c(0, 10), ylim = c(5, 13), 
     xlab = expression(paste("Total Reads (", ''%*%10^6, ")")), 
     ylab = "Shannon index")
title("D",  adj = 0, line = 0.8, cex.main = 1.8)
axis(1, depth[c(1, 10:19)] / 1e6, depth[c(1, 10:19)] / 1e6)
axis(2, 5:13, 5:13)
for(i in 1:nrow(dataAlpha1))lines(depth / 1e6, dataAlpha1[i, ])

boxplot(dataAlpha1, at = depth / 1e5, outline = FALSE, ylim = c(5, 13),
        xlab = expression(paste("Total Reads (", ''%*%10^6, ")")), ylab = "Shannon index", names = NA)
title("E",  adj = 0, line = 1.4, cex.main = 1.8)
axis(1, depth[c(1, 10:19)] / 1e5, depth[c(1, 10:19)] / 1e6)
axis(2, 5:13, 5:13)
mtext(colSums(!is.na(dataAlpha1[, c(1, 10:19)])), 3, at = depth[c(1, 10:19)] / 1e5, font = 2)
mtext("Number of Samples Available", 3, line = 1.3, font = 4)
lines(lowess(rep(depth / 1e5, each = nSample)[!is.na(dataAlpha1)], dataAlpha1[!is.na(dataAlpha1)]), col = rgb(1, 0, 0, 0.5), lwd = 4)

tempX <- totalReads
tempY <- dataAlpha1[, "rarefaction_1e+06"]
tempIdx0 <- which(totalReads < 1e7 & !is.na(tempY))
tempLowess0 <- lowess(tempX[tempIdx0], tempY[tempIdx0])
tempTest0 <- cor.test(tempX[tempIdx0], tempY[tempIdx0], method = "spearman")

plot(NA, xaxt = "n", yaxt = "n",
     xlim = c(0, 10)*1e6, ylim = c(5, 13), 
     xlab = expression(paste("Total Reads (", ''%*%10^6, ")")),
     ylab = expression(paste("Shannon index")))
title("F",  adj = 0, line = 0.8, cex.main = 1.8)
points(x=tempX[tempIdx0], y=tempY[tempIdx0])
axis(1, depth[10:19], depth[10:19]/1e6)
axis(2, 5:13, 5:13)
lines(tempLowess0, col = "red", lwd = 2)
tempPVec <- tempTest0$p.value
legend("bottomright", paste0("P = ", signif(tempPVec, 3)), cex = 1.5, bty = "n", text.col = "red")

dev.off()





rm(list = ls())

library(gplots) #?plot.lowess

load("TCR_rarefy/readsList")
load("TCR_rarefy/dataPheno")
load("TCR_rarefy/dataAlpha")

sampleID  = sampleID[grep("P", sampleID)]
nSample   = length(sampleID)
dataPheno = dataPheno[sampleID, ]

ci = c(1e6, 2e6, 4e6, 8e6, 1e7)
mycolors = c("black", "red", "orange", "green", "blue", "pink", "cyan")

myID = function(tempX, tempY, ci){
  seg = c(min(tempX), ci, Inf)
  tempIdClc = vector("list", length(ci)+1)
  tempIdClc[[1]] = which(!is.na(tempY))
  for(i in 1:(length(ci)+1)){
    tempIdClc[[i+1]] = which(tempX >= seg[i] & tempX < seg[i+1] & !is.na(tempY))
  }
  return(tempIdClc)
}


pdf("TCR_rarefy/2BinsDetermine/Figure4.pdf", width = 17, height = 17)
par(mfrow = c(2, 2), cex.main = 2, cex.lab = 2, cex.axis = 2, mar = c(6, 6, 6, 1))

##########

myXaxis = c(0, 1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 1e7)

tempX = totalReads
tempY = dataAlpha2[, 1] #1e5
tempIdClc = myID(tempX, tempY, ci)

plot(tempX[tempIdClc[[1]]], tempY[tempIdClc[[1]]], col = "gray30", cex=1.5,
     xaxt = "n", yaxt = "n", xlim=c(0,10)*1e6, ylim = c(0, 10.5)*1e4,
     xlab = expression(paste("Library Size (", ''%*%10^6, ")")), 
     ylab = expression(paste("Unique sequence counts (", ''%*%10^4, ")")))
title("A: Unique sequence counts (rarefying to 1e5)",  adj = 0, line = 0.8)
axis(1, myXaxis, myXaxis/ 1e6)
axis(2, (0:10)*1e4, 0:10)
for(i in 1:length(ci)){abline(v = ci[i], lty = 2, col = "gray", lwd=3.5)}

tempPVec = rep(NA, length(tempIdClc))
for(i in 1:length(tempIdClc)){
  tempIdxx = tempIdClc[[i]]
  lines(lowess(tempX[tempIdxx], tempY[tempIdxx]), col=mycolors[i], lwd=3.5)
  tempcortest = cor.test(tempX[tempIdxx], tempY[tempIdxx], method = "spearman")
  if(tempcortest$p.value < 0.001){
    tempPVec[i] = "P < 0.001"
  }else{
    tempPVec[i] = paste0("P = ", format(round(tempcortest$p.value, 3), nsmall=3))
  }
}
legend("topleft", tempPVec, cex=2.1, bty = "n", text.col=mycolors, ncol=4)

##########

tempX = totalReads
tempY = dataAlpha2[, 1] #1e5
for(i in 3:length(tempIdClc)){
  tempIdxx = tempIdClc[[i]]
  tempY[tempIdxx] = dataAlpha2[tempIdxx, paste0("rarefaction_", ci[i-2])]
}

tempIdClc = myID(tempX, tempY, ci)

plot(tempX[tempIdClc[[1]]], tempY[tempIdClc[[1]]], col = "gray30", cex=1.5,
     xaxt = "n", yaxt = "n", xlim=c(0,10)*1e6, ylim = c(0, 10.5)*1e5,
     xlab = expression(paste("Library Size (", ''%*%10^6, ")")), 
     ylab = expression(paste("Unique sequence counts (", ''%*%10^5, ")")))
title("C: Unique sequence counts (rarefying within each bin)",  adj = 0, line = 0.8)
axis(1, myXaxis, myXaxis/ 1e6)
axis(2, (0:10)*1e5, 0:10)
for(i in 1:length(ci)){abline(v = ci[i], lty = 2, col = "gray", lwd=3.5)}

tempPVec = rep(NA, length(tempIdClc))
for(i in 1:length(tempIdClc)){
  tempIdxx = tempIdClc[[i]]
  lines(lowess(tempX[tempIdxx], tempY[tempIdxx]), col=mycolors[i], lwd=3.5)
  tempcortest = cor.test(tempX[tempIdxx], tempY[tempIdxx], method = "spearman")
  if(tempcortest$p.value < 0.001){
    tempPVec[i] = "P < 0.001"
  }else{
    tempPVec[i] = paste0("P = ", format(round(tempcortest$p.value, 3), nsmall=3))
  }
}
legend("topleft", tempPVec, cex=2.1, bty = "n", text.col=mycolors, ncol=4)

##########

tempX = totalReads
tempY = dataAlpha1[, 1] #1e5
tempIdClc = myID(tempX, tempY, ci)

plot(tempX[tempIdClc[[1]]]/ 1e6, tempY[tempIdClc[[1]]], col = "gray30", cex=1.5,
     xaxt = "n", yaxt = "n", xlim=c(0,10), ylim = c(5, 14.3),
     xlab = expression(paste("Library Size (", ''%*%10^6, ")")),  
     ylab = "Shannon index")
title("B: Shannon index (rarefying to 1e5)",  adj = 0, line = 0.8)
axis(1, myXaxis/1e6, myXaxis/ 1e6)
axis(2, 5:14, 5:14)
for(i in 1:length(ci)){abline(v = ci[i]/1e6, lty = 2, col = "gray", lwd=3.5)}
# abline(v = 1e6 / 1e6, lty = 2, col = "gray")

tempPVec = rep(NA, length(tempIdClc))
for(i in 1:length(tempIdClc)){
  tempIdxx = tempIdClc[[i]]
  lines(lowess(tempX[tempIdxx]/ 1e6, tempY[tempIdxx]), col=mycolors[i], lwd=3.5)
  tempcortest = cor.test(tempX[tempIdxx], tempY[tempIdxx], method = "spearman")
  if(tempcortest$p.value < 0.001){
    tempPVec[i] = "P < 0.001"
  }else{
    tempPVec[i] = paste0("P = ", format(round(tempcortest$p.value, 3), nsmall=3))
  }
}
legend("topleft", tempPVec, cex=2.1, bty = "n", text.col=mycolors, ncol=4)

##########

tempX = totalReads
tempY = dataAlpha1[, 1] #1e5
for(i in 3:length(tempIdClc)){
  tempIdxx = tempIdClc[[i]]
  tempY[tempIdxx] = dataAlpha1[tempIdxx, paste0("rarefaction_", ci[i-2])]
}

tempIdClc = myID(tempX, tempY, ci)

plot(tempX[tempIdClc[[1]]] / 1e6, tempY[tempIdClc[[1]]], col = "gray30", cex=1.5,
     xaxt = "n", yaxt = "n", xlim=c(0,10), ylim = c(5, 14.3),
     xlab = expression(paste("Library Size (", ''%*%10^6, ")")), 
     ylab = "Shannon index")
title("D: Shannon index (rarefying within each bin)",  adj = 0, line = 0.8)
axis(1, myXaxis/1e6, myXaxis/ 1e6)
axis(2, 5:14, 5:14)
for(i in 1:length(ci)){abline(v = ci[i]/1e6, lty = 2, col = "gray", lwd=3.5)}

tempPVec = rep(NA, length(tempIdClc))
for(i in 1:length(tempIdClc)){
  tempIdxx = tempIdClc[[i]]
  lines(lowess(tempX[tempIdxx]/ 1e6, tempY[tempIdxx]), col=mycolors[i], lwd=3.5)
  tempcortest = cor.test(tempX[tempIdxx], tempY[tempIdxx], method = "spearman")
  if(tempcortest$p.value < 0.001){
    tempPVec[i] = "P < 0.001"
  }else{
    tempPVec[i] = paste0("P = ", format(round(tempcortest$p.value, 3), nsmall=3))
  }
}
legend("topleft", tempPVec, cex=2.1, bty = "n", text.col=mycolors, ncol=4)
dev.off()









rm(list = ls())



red = c(071, 231, 247, 055)
green = c(083, 098, 170, 103)
blue = c(120, 084, 088, 149)
colors = rgb(red, green, blue, m=255)

library(gplots) #?plot.lowess


load("TCR_rarefy/readsList")
load("TCR_rarefy/dataPheno")
load("TCR_rarefy/dataAlpha")

sampleID <- sampleID[grep("P", sampleID)]
nSample <- length(sampleID)
dataPheno <- dataPheno[sampleID, ]


ci = c(1e6, 2e6, 4e6, 8e6, 1e7)

# Shannon index
myshannon = unlist(lapply(readsList, function(x){p <- x / sum(x); return(-sum(p * log(p)))}))
tmpshannon = dataAlpha1[, "rarefaction_1e+06"]
myshannon = cbind(myshannon, tmpshannon)
tmpshannon = dataAlpha1[, "rarefaction_1e+05"]
myshannon = cbind(myshannon, tmpshannon)
for(i in 1:length(ci)){
  tmpshannon = dataAlpha1[, paste0("rarefaction_", ci[i])]
  myshannon = cbind(myshannon, tmpshannon)
}
colnames(myshannon) = c("shannon", "shannon_overall", paste0("shannon_bin", 1:(length(ci)+1)))
dataPheno = cbind(dataPheno, myshannon)

# observed unique sequences
myobsvec = unlist(lapply(readsList, length))
tmpobsvec = dataAlpha2[, "rarefaction_1e+06"]
myobsvec = cbind(myobsvec, tmpobsvec)
tmpobsvec = dataAlpha2[, "rarefaction_1e+05"]
myobsvec = cbind(myobsvec, tmpobsvec)
for(i in 1:length(ci)){
  tmpobsvec = dataAlpha2[, paste0("rarefaction_", ci[i])]
  myobsvec = cbind(myobsvec, tmpobsvec)
}
colnames(myobsvec) = c("obsvec", "obsvec_overall", paste0("obsvec_bin", 1:(length(ci)+1)))
dataPheno = cbind(dataPheno, myobsvec)

# clinical covariates
dataPheno$yage        = dataPheno$Age
dataPheno$ysex        = (dataPheno$'Biological Sex' == "Female") + 0L
dataPheno$yinfercmv   = (dataPheno$'inferred CMV status' == "Inferred CMV +") + 0L
dataPheno$ycmv        = (dataPheno$'Virus Diseases' == "Cytomegalovirus +") + 0L
dataPheno$ytotalreads = totalReads




pdf("TCR_rarefy/5AppSection/Figure5.pdf", width = 23, height = 18)
par(mfrow=c(3,4), cex.main = 3, cex.axis = 2.6, cex.lab = 2.6, mar = c(6, 6, 6, 1))

plot(x=dataPheno$yage, y=dataPheno$obsvec_overall, xaxt="n", yaxt="n", ylim=c(0,4)*1e5,
     xlab = "", ylab = expression(paste("Unique sequences counts (", ''%*%10^5, ")")),
     main = "Age")
lines(lowess(dataPheno$obsvec_overall~dataPheno$yage), col=colors[2], lwd=4)
tmpp = cor.test(dataPheno$obsvec_overall, dataPheno$yage, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 40, y = 0.05*1e5, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, seq(0, 80, by=10), seq(0, 80, by=10))
axis(2, seq(0, 7e5, by=1e5), seq(0, 7e5, by=1e5)/1e5)


boxplot(obsvec_overall ~ ysex, dataPheno, xaxt="n", yaxt="n", col=colors[3:4], ylim=c(0,4)*1e5,
        xlab = "", ylab = expression(paste("Unique clonotypes counts (", ''%*%10^5, ")")),
        main = "Gender")
tmpp = cor.test(dataPheno$obsvec_overall, dataPheno$ysex, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 1.5, y = 0.05*1e5, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, 1:2, c("Female", "Male"))
axis(2, seq(0, 7e5, by=1e5), seq(0, 7e5, by=1e5)/1e5)




boxplot(obsvec_overall ~ yinfercmv, dataPheno, xaxt="n", yaxt="n", col=colors[3:4], ylim=c(0,4)*1e5,
        xlab = "", ylab = expression(paste("Unique clonotypes counts (", ''%*%10^5, ")")),
        main = "Inferred CMV")
tmpp = cor.test(dataPheno$obsvec_overall, dataPheno$yinfercmv, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 1.5, y = 0.05*1e5, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, 1:2, c("Inferred CMV-", "Inferred CMV+"))
axis(2, seq(0, 7e5, by=1e5), seq(0, 7e5, by=1e5)/1e5)


boxplot(obsvec_overall ~ ycmv, dataPheno, xaxt="n", yaxt="n", col=colors[3:4], ylim=c(0,4)*1e5,
        xlab = "", ylab = expression(paste("Unique clonotypes counts (", ''%*%10^5, ")")),
        main = "CMV")
tmpp = cor.test(dataPheno$obsvec_overall, dataPheno$ycmv, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 1.5, y = 0.05*1e5, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, 1:2, c("CMV -", "CMV +"))
axis(2, seq(0, 7e5, by=1e5), seq(0, 7e5, by=1e5)/1e5)




plot(x=dataPheno$yage, y=dataPheno$shannon_overall, ylim=c(5, 13), xaxt="n", yaxt="n", col=colors[1],
     xlab = "", ylab = "Shannon index",
     main = "Age")
lines(lowess(dataPheno$shannon_overall~dataPheno$yage), col=colors[2], lwd=4)
tmpp = cor.test(dataPheno$shannon_overall, dataPheno$yage, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 40, y = 5.1, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, seq(0, 80, by=10), seq(0, 80, by=10))
axis(2, 5:13, 5:13)

boxplot(shannon_overall ~ ysex, dataPheno, ylim=c(5, 13), xaxt="n", yaxt="n", col=colors[3:4],
        xlab = "", ylab = "Shannon index",
        main = "Gender")
tmpp = cor.test(dataPheno$shannon_overall, dataPheno$ysex, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 1.5, y = 5.1, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, 1:2, c("Female", "Male"))
axis(2, 5:13, 5:13)

boxplot(shannon_overall ~ yinfercmv, dataPheno, ylim=c(5, 13), xaxt="n", yaxt="n", col=colors[3:4],
        xlab = "", ylab = "Shannon index",
        main = "Inferred CMV")
tmpp = cor.test(dataPheno$shannon_overall, dataPheno$yinfercmv, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 1.5, y = 5.1, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, 1:2, c("Inferred CMV-", "Inferred CMV+"))
axis(2, 5:13, 5:13)

boxplot(shannon_overall ~ ycmv, dataPheno, ylim=c(5, 13), xaxt="n", yaxt="n", col=colors[3:4],
        xlab="", ylab = "Shannon index",
        main = "CMV")
tmpp = cor.test(dataPheno$shannon_overall, dataPheno$ycmv, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 1.5, y = 5.1, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, 1:2, c("CMV-", "CMV+"))
axis(2, 5:13, 5:13)


plot(x=dataPheno$yage, y=log10(dataPheno$ytotalreads), xaxt="n", yaxt="n", col=colors[1],
     xlab = "", ylab = expression('log'[10]('library size')), 
     main = "Age")
lines(lowess(log10(dataPheno$ytotalreads)~dataPheno$yage), col=colors[2], lwd=4)
tmpp = cor.test(dataPheno$ytotalreads, dataPheno$Age, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 40, y = 4.45, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, seq(0, 80, by=10), seq(0, 80, by=10))
axis(2, seq(4.5, 8, by=0.5), seq(4.5, 8, by=0.5))

boxplot(log10(ytotalreads) ~ ysex, dataPheno, xaxt="n", yaxt="n", col=colors[3:4],
        xlab = "", ylab = expression('log'[10]('library size')),
        main = "Gender")
tmpp = cor.test(dataPheno$ytotalreads, dataPheno$ysex, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 1.5, y = 4.45, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, 1:2, c("Female", "Male"))
axis(2, seq(4.5, 8, by=0.5), seq(4.5, 8, by=0.5))



boxplot(log10(ytotalreads) ~ yinfercmv, dataPheno, xaxt="n", yaxt="n", col=colors[3:4],
        xlab = "", ylab = expression('log'[10]('library size')), 
        main = "Inferred CMV")
tmpp = cor.test(dataPheno$ytotalreads, dataPheno$yinfercmv, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 1.5, y = 4.45, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, 1:2, c("Inferred CMV-", "Inferred CMV+"))
axis(2, seq(4.5, 8, by=0.5), seq(4.5, 8, by=0.5))

boxplot(log10(ytotalreads) ~ ycmv, dataPheno, xaxt="n", yaxt="n", col=colors[3:4],
        xlab = "", ylab = expression('log'[10]('library size')),
        main = "CMV")
tmpp = cor.test(dataPheno$ytotalreads, dataPheno$ycmv, method = "spearman")$p.value
if(tmpp < 0.0001){
  tmpp = formatC(tmpp, format="e", digits = 2)
}else{
  tmpp = round(tmpp, 4)
}
text(x = 1.5, y = 4.45, paste0("P-value = ", tmpp), col=colors[2], cex = 3)
axis(1, 1:2, c("CMV-", "CMV+"))
axis(2, seq(4.5, 8, by=0.5), seq(4.5, 8, by=0.5))

dev.off()







binFactor <- cut(totalReads, c(1e5, ci, Inf))
levels(binFactor) <- c(1e5, ci)

tempIdxSubsample  = which(!is.na(binFactor))
tempIdxOverall    = which(!is.na(binFactor) & totalReads >= 1e6)

seg = c(1e5, ci, Inf)
tempIdxBinClc = vector("list", length(ci)+1)
for(i in 1:(length(ci)+1)){
  tempIdxBinClc[[i]] = which(!is.na(binFactor) & totalReads >= seg[i] & totalReads < seg[i+1])
}

ApplyRes_obsvec = as.data.frame(matrix(0, nrow=6+length(ci)+1, ncol=6))
colnames(ApplyRes_obsvec) = c("samplesize", "age", "sex", "infercmv", "cmv", "totalreads")
rownames(ApplyRes_obsvec) = c("Norarefy", "Overall", "LOESS", paste0("Bin", 1:(length(ci)+1)), 
                              "Multi-bin-Equal", "Multi-bin-SSW", "Multi-bin-IVW")
ApplyRes_shannon = ApplyRes_obsvec

##################################################################
################################################################## No rarefaction
##################################################################

# Unique sequences counts
age_obsvec_m1        = summary(lm(obsvec ~ yage, data = dataPheno))
sex_obsvec_m1        = summary(lm(obsvec ~ ysex, data = dataPheno))
infercmv_obsvec_m1   = summary(lm(obsvec ~ yinfercmv, data = dataPheno))
cmv_obsvec_m1        = summary(lm(obsvec ~ ycmv, data = dataPheno))
totalreads_obsvec_m1 = summary(lm(obsvec ~ ytotalreads, data = dataPheno))

ApplyRes_obsvec["Norarefy", "samplesize"] = sum(1*(!is.na(dataPheno$obsvec)))
ApplyRes_obsvec["Norarefy", "age"] = paste0(round(age_obsvec_m1$coefficients[2,1]/(10^4),3), " (",
                                            formatC(age_obsvec_m1$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["Norarefy", "sex"] = paste0(round(sex_obsvec_m1$coefficients[2,1]/(10^4),3), " (",
                                            formatC(sex_obsvec_m1$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["Norarefy", "infercmv"] = paste0(round(infercmv_obsvec_m1$coefficients[2,1]/(10^4),3), " (",
                                                 formatC(infercmv_obsvec_m1$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["Norarefy", "cmv"] = paste0(round(cmv_obsvec_m1$coefficients[2,1]/(10^4),3), " (",
                                            formatC(cmv_obsvec_m1$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["Norarefy", "totalreads"] = paste0(round(totalreads_obsvec_m1$coefficients[2,1]/(10^-2),3), " (",
                                                   formatC(totalreads_obsvec_m1$coefficients[2,4], format = "e", digits = 2), ")")

# Shannon index 
age_shannon_m1        = summary(lm(shannon ~ yage, data = dataPheno))
sex_shannon_m1        = summary(lm(shannon ~ ysex, data = dataPheno))
infercmv_shannon_m1   = summary(lm(shannon ~ yinfercmv, data = dataPheno))
cmv_shannon_m1        = summary(lm(shannon ~ ycmv, data = dataPheno))
totalreads_shannon_m1 = summary(lm(shannon ~ ytotalreads, data = dataPheno))

ApplyRes_shannon["Norarefy", "samplesize"] = sum(1*(!is.na(dataPheno$shannon)))
ApplyRes_shannon["Norarefy", "age"] = paste0(round(age_shannon_m1$coefficients[2,1],3), " (",
                                             formatC(age_shannon_m1$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["Norarefy", "sex"] = paste0(round(sex_shannon_m1$coefficients[2,1],3), " (",
                                             formatC(sex_shannon_m1$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["Norarefy", "infercmv"] = paste0(round(infercmv_shannon_m1$coefficients[2,1],3), " (",
                                                  formatC(infercmv_shannon_m1$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["Norarefy", "cmv"] = paste0(round(cmv_shannon_m1$coefficients[2,1],3), " (",
                                             formatC(cmv_shannon_m1$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["Norarefy", "totalreads"] = paste0(round(totalreads_shannon_m1$coefficients[2,1]/(10^-7),3), " (",
                                                    formatC(totalreads_shannon_m1$coefficients[2,4], format = "e", digits = 2), ")")

##################################################################
################################################################## Overall rarefying
##################################################################
# discard samples with librarysize < 1e6 & subsample sequences to 1e6

# Unique sequences counts
age_obsvec_m3        = summary(lm(obsvec_overall ~ yage, data = dataPheno[tempIdxOverall, ]))
sex_obsvec_m3        = summary(lm(obsvec_overall ~ ysex, data = dataPheno[tempIdxOverall, ]))
infercmv_obsvec_m3   = summary(lm(obsvec_overall ~ yinfercmv, data = dataPheno[tempIdxOverall, ]))
cmv_obsvec_m3        = summary(lm(obsvec_overall ~ ycmv, data = dataPheno[tempIdxOverall, ]))
totalreads_obsvec_m3 = summary(lm(obsvec_overall ~ ytotalreads, data = dataPheno[tempIdxOverall, ]))

ApplyRes_obsvec["Overall","samplesize"] = sum(1*(!is.na(dataPheno$obsvec_overall[tempIdxOverall])))
ApplyRes_obsvec["Overall", "age"] = paste0(round(age_obsvec_m3$coefficients[2,1]/(10^4),3), " (",
                                           formatC(age_obsvec_m3$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["Overall", "sex"] = paste0(round(sex_obsvec_m3$coefficients[2,1]/(10^4),3), " (",
                                           formatC(sex_obsvec_m3$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["Overall", "infercmv"] = paste0(round(infercmv_obsvec_m3$coefficients[2,1]/(10^4),3), " (",
                                                formatC(infercmv_obsvec_m3$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["Overall", "cmv"] = paste0(round(cmv_obsvec_m3$coefficients[2,1]/(10^4),3), " (",
                                           formatC(cmv_obsvec_m3$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["Overall", "totalreads"] = paste0(round(totalreads_obsvec_m3$coefficients[2,1]/(10^-2),3), " (",
                                                  formatC(totalreads_obsvec_m3$coefficients[2,4], format = "e", digits = 2), ")")

# Shannon index 
age_shannon_m3        = summary(lm(shannon_overall ~ yage, data = dataPheno[tempIdxOverall, ]))
sex_shannon_m3        = summary(lm(shannon_overall ~ ysex, data = dataPheno[tempIdxOverall, ]))
infercmv_shannon_m3   = summary(lm(shannon_overall ~ yinfercmv, data = dataPheno[tempIdxOverall, ]))
cmv_shannon_m3        = summary(lm(shannon_overall ~ ycmv, data = dataPheno[tempIdxOverall, ]))
totalreads_shannon_m3 = summary(lm(shannon_overall ~ ytotalreads, data = dataPheno[tempIdxOverall, ]))

ApplyRes_shannon["Overall","samplesize"] = sum(1*(!is.na(dataPheno$shannon_overall[tempIdxOverall])))
ApplyRes_shannon["Overall", "age"] = paste0(round(age_shannon_m3$coefficients[2,1],3), " (",
                                            formatC(age_shannon_m3$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["Overall", "sex"] = paste0(round(sex_shannon_m3$coefficients[2,1],3), " (",
                                            formatC(sex_shannon_m3$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["Overall", "infercmv"] = paste0(round(infercmv_shannon_m3$coefficients[2,1],3), " (",
                                                 formatC(infercmv_shannon_m3$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["Overall", "cmv"] = paste0(round(cmv_shannon_m3$coefficients[2,1],3), " (",
                                            formatC(cmv_shannon_m3$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["Overall", "totalreads"] = paste0(round(totalreads_shannon_m3$coefficients[2,1]/(10^-7),3), " (",
                                                   formatC(totalreads_shannon_m3$coefficients[2,4], format = "e", digits = 2), ")")

##################################################################
################################################################## LOWESS fitting with span=0.75
##################################################################

myspan = 0.5

# Unique sequences counts
temploessOUS = loess(obsvec ~ ytotalreads, data = dataPheno, span = myspan)
yOUSloess = temploessOUS$residuals

age_obsvec_m8        = summary(lm(yOUSloess ~ yage, data = dataPheno))
sex_obsvec_m8        = summary(lm(yOUSloess ~ ysex, data = dataPheno))
infercmv_obsvec_m8   = summary(lm(yOUSloess ~ yinfercmv, data = dataPheno))
cmv_obsvec_m8        = summary(lm(yOUSloess ~ ycmv, data = dataPheno))
totalreads_obsvec_m8 = summary(lm(yOUSloess ~ ytotalreads, data = dataPheno))

ApplyRes_obsvec["LOESS","samplesize"] = sum(1*(!is.na(yOUSloess)))
ApplyRes_obsvec["LOESS", "age"] = paste0(round(age_obsvec_m8$coefficients[2,1]/(10^4),3), " (",
                                          formatC(age_obsvec_m8$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["LOESS", "sex"] = paste0(round(sex_obsvec_m8$coefficients[2,1]/(10^4),3), " (",
                                          formatC(sex_obsvec_m8$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["LOESS", "infercmv"] = paste0(round(infercmv_obsvec_m8$coefficients[2,1]/(10^4),3), " (",
                                               formatC(infercmv_obsvec_m8$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["LOESS", "cmv"] = paste0(round(cmv_obsvec_m8$coefficients[2,1]/(10^4),3), " (",
                                          formatC(cmv_obsvec_m8$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_obsvec["LOESS", "totalreads"] = paste0(round(totalreads_obsvec_m8$coefficients[2,1]/(10^-2),3), " (",
                                                 formatC(totalreads_obsvec_m8$coefficients[2,4], format = "e", digits = 2), ")")


# Shannon index 

temploessSI = loess(shannon ~ ytotalreads, data = dataPheno, span = myspan)
ySIloess = temploessSI$residuals

age_shannon_m8        = summary(lm(ySIloess ~ yage, data = dataPheno))
sex_shannon_m8        = summary(lm(ySIloess ~ ysex, data = dataPheno))
infercmv_shannon_m8   = summary(lm(ySIloess ~ yinfercmv, data = dataPheno))
cmv_shannon_m8        = summary(lm(ySIloess ~ ycmv, data = dataPheno))
totalreads_shannon_m8 = summary(lm(ySIloess ~ ytotalreads, data = dataPheno))

ApplyRes_shannon["LOESS","samplesize"] = sum(1*(!is.na(ySIloess)))
ApplyRes_shannon["LOESS", "age"] = paste0(round(age_shannon_m8$coefficients[2,1],3), " (",
                                           formatC(age_shannon_m8$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["LOESS", "sex"] = paste0(round(sex_shannon_m8$coefficients[2,1],3), " (",
                                           formatC(sex_shannon_m8$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["LOESS", "infercmv"] = paste0(round(infercmv_shannon_m8$coefficients[2,1],3), " (",
                                                formatC(infercmv_shannon_m8$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["LOESS", "cmv"] = paste0(round(cmv_shannon_m8$coefficients[2,1],3), " (",
                                           formatC(cmv_shannon_m8$coefficients[2,4], format = "e", digits = 2), ")")
ApplyRes_shannon["LOESS", "totalreads"] = paste0(round(totalreads_shannon_m8$coefficients[2,1]/(10^-7),3), " (",
                                                  formatC(totalreads_shannon_m8$coefficients[2,4], format = "e", digits = 2), ")")

##################################################################
################################################################## Multi-bin
##################################################################

age_obsvec = age_shannon = vector("list", length(ci)+1)
sex_obsvec = sex_shannon = vector("list", length(ci)+1)
infercmv_obsvec = infercmv_shannon = vector("list", length(ci)+1)
cmv_obsvec = cmv_shannon = vector("list", length(ci)+1)
totalreads_obsvec = totalreads_shannon = vector("list", length(ci)+1)
averagetotalreads = rep(0, length(ci)+1)

for(i in 1:(length(ci)+1)){
  
  tempIdxBin = tempIdxBinClc[[i]]
  
  averagetotalreads[i] = mean(dataPheno[tempIdxBin, "ytotalreads"])
  
  # Unique sequences counts
  age_obsvec_m        = summary(lm(as.formula(paste0("obsvec_bin", i, " ~ yage")), data = dataPheno[tempIdxBin, ]))
  sex_obsvec_m        = summary(lm(as.formula(paste0("obsvec_bin", i, " ~ ysex")), data = dataPheno[tempIdxBin, ]))
  infercmv_obsvec_m   = summary(lm(as.formula(paste0("obsvec_bin", i, " ~ yinfercmv")), data = dataPheno[tempIdxBin, ]))
  cmv_obsvec_m        = summary(lm(as.formula(paste0("obsvec_bin", i, " ~ ycmv")), data = dataPheno[tempIdxBin, ]))
  totalreads_obsvec_m = summary(lm(as.formula(paste0("obsvec_bin", i, " ~ ytotalreads")), data = dataPheno[tempIdxBin, ]))
  
  age_obsvec[[i]] = age_obsvec_m
  sex_obsvec[[i]] = sex_obsvec_m
  infercmv_obsvec[[i]] = infercmv_obsvec_m
  cmv_obsvec[[i]] = cmv_obsvec_m
  totalreads_obsvec[[i]] = totalreads_obsvec_m
  
  ApplyRes_obsvec[paste0("Bin",i),"samplesize"] = sum(1*(!is.na(dataPheno[tempIdxBin, paste0("obsvec_bin", i)])))
  ApplyRes_obsvec[paste0("Bin",i), "age"] = paste0(round(age_obsvec_m$coefficients[2,1]/(10^4),3), " (",
                                                   formatC(age_obsvec_m$coefficients[2,4], format = "e", digits = 2), ")")
  ApplyRes_obsvec[paste0("Bin",i), "sex"] = paste0(round(sex_obsvec_m$coefficients[2,1]/(10^4),3), " (",
                                                   formatC(sex_obsvec_m$coefficients[2,4], format = "e", digits = 2), ")")
  ApplyRes_obsvec[paste0("Bin",i), "infercmv"] = paste0(round(infercmv_obsvec_m$coefficients[2,1]/(10^4),3), " (",
                                                        formatC(infercmv_obsvec_m$coefficients[2,4], format = "e", digits = 2), ")")
  ApplyRes_obsvec[paste0("Bin",i), "cmv"] = paste0(round(cmv_obsvec_m$coefficients[2,1]/(10^4),3), " (",
                                                   formatC(cmv_obsvec_m$coefficients[2,4], format = "e", digits = 2), ")")
  ApplyRes_obsvec[paste0("Bin",i), "totalreads"] = paste0(round(totalreads_obsvec_m$coefficients[2,1]/(10^-2),3), " (",
                                                          formatC(totalreads_obsvec_m$coefficients[2,4], format = "e", digits = 2), ")")
  
  # Shannon index 
  age_shannon_m        = summary(lm(paste0("shannon_bin", i, " ~ yage"), data = dataPheno[tempIdxBin, ]))
  sex_shannon_m        = summary(lm(paste0("shannon_bin", i, " ~ ysex"), data = dataPheno[tempIdxBin, ]))
  infercmv_shannon_m   = summary(lm(paste0("shannon_bin", i, "~ yinfercmv"), data = dataPheno[tempIdxBin, ]))
  cmv_shannon_m        = summary(lm(paste0("shannon_bin", i, "~ ycmv"), data = dataPheno[tempIdxBin, ]))
  totalreads_shannon_m = summary(lm(paste0("shannon_bin", i, "~ ytotalreads"), data = dataPheno[tempIdxBin, ]))
  
  age_shannon[[i]] = age_shannon_m
  sex_shannon[[i]] = sex_shannon_m
  infercmv_shannon[[i]] = infercmv_shannon_m
  cmv_shannon[[i]] = cmv_shannon_m
  totalreads_shannon[[i]] = totalreads_shannon_m
  
  ApplyRes_shannon[paste0("Bin",i),"samplesize"] = sum(1*(!is.na(dataPheno[tempIdxBin, paste0("obsvec_bin", i)])))
  ApplyRes_shannon[paste0("Bin",i), "age"] = paste0(round(age_shannon_m$coefficients[2,1],3), " (",
                                                    formatC(age_shannon_m$coefficients[2,4], format = "e", digits = 2), ")")
  ApplyRes_shannon[paste0("Bin",i), "sex"] = paste0(round(sex_shannon_m$coefficients[2,1],3), " (",
                                                    formatC(sex_shannon_m$coefficients[2,4], format = "e", digits = 2), ")")
  ApplyRes_shannon[paste0("Bin",i), "infercmv"] = paste0(round(infercmv_shannon_m$coefficients[2,1],3), " (",
                                                         formatC(infercmv_shannon_m$coefficients[2,4], format = "e", digits = 2), ")")
  ApplyRes_shannon[paste0("Bin",i), "cmv"] = paste0(round(cmv_shannon_m$coefficients[2,1],3), " (",
                                                    formatC(cmv_shannon_m$coefficients[2,4], format = "e", digits = 2), ")")
  ApplyRes_shannon[paste0("Bin",i), "totalreads"] = paste0(round(totalreads_shannon_m$coefficients[2,1]/(10^-7),3), " (",
                                                           formatC(totalreads_shannon_m$coefficients[2,4], format = "e", digits = 2), ")")
}

##################################################################
################################################################## Multi-bins with Equal and optimal weights
##################################################################

Multibin = function(myfit, vec_size, averagetotalreads){
  
  numfit = length(myfit)
  bClc = vClc = wClc = rep(NA, numfit)
  for(i in 1:numfit){
    bClc[i] = myfit[[i]]$coefficients[2, "Estimate"]
    vClc[i] = myfit[[i]]$coefficients[2, "Std. Error"]^2
  }
  
  # Multi-bin-Equal
  T.param.equal <- mean(bClc)
  T.equal <- sqrt((sum(bClc))^2 / sum(vClc))
  p.equal <- 2 * pnorm(T.equal, lower.tail = FALSE)
  
  # Multi-bin-SSW
  w1Clc = vec_size
  T.param.1 <- sum(w1Clc*bClc)/sum(w1Clc)
  T.optim.1 <- sqrt((sum(w1Clc*bClc)^2) / sum(w1Clc*w1Clc*vClc))
  p.optim.1 <- 2 * pnorm(T.optim.1, lower.tail = FALSE)
  
  # Multi-bin-IVW
  w2Clc = (1/vClc)
  T.param.2 <- sum(w2Clc*bClc)/sum(w2Clc)
  T.optim.2 <- sqrt((sum(w2Clc*bClc))^2/ sum(w2Clc*w2Clc*vClc))
  p.optim.2 <- 2 * pnorm(T.optim.2, lower.tail = FALSE)
  
  return(list(T.param.equal=T.param.equal, T.param.1=T.param.1, T.param.2=T.param.2,
              p.equal=p.equal, p.optim.1=p.optim.1, p.optim.2=p.optim.2))
}

##### Unique sequences counts
vec_size = ApplyRes_obsvec[paste0("Bin", 1:(length(ci)+1)), "samplesize"]

ApplyRes_obsvec["Multi-bin-Equal","samplesize"] = sum(vec_size)
ApplyRes_obsvec["Multi-bin-SSW","samplesize"] = sum(vec_size)
ApplyRes_obsvec["Multi-bin-IVW","samplesize"] = sum(vec_size)

multi_obsvec_age = Multibin(age_obsvec, vec_size, averagetotalreads)
ApplyRes_obsvec["Multi-bin-Equal", "age"] = paste0(round(multi_obsvec_age$T.param.equal/(10^4),3), " (",
                                                 formatC(multi_obsvec_age$p.equal, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-SSW", "age"] = paste0(round(multi_obsvec_age$T.param.1/(10^4),3), " (",
                                                  formatC(multi_obsvec_age$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-IVW", "age"] = paste0(round(multi_obsvec_age$T.param.2/(10^4),3), " (",
                                                      formatC(multi_obsvec_age$p.optim.2, format = "e", digits = 2), ")")

multi_obsvec_sex = Multibin(sex_obsvec, vec_size, averagetotalreads)
ApplyRes_obsvec["Multi-bin-Equal", "sex"] = paste0(round(multi_obsvec_sex$T.param.equal/(10^4),3), " (",
                                                 formatC(multi_obsvec_sex$p.equal, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-SSW", "sex"] = paste0(round(multi_obsvec_sex$T.param.1/(10^4),3), " (",
                                                  formatC(multi_obsvec_sex$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-IVW", "sex"] = paste0(round(multi_obsvec_sex$T.param.2/(10^4),3), " (",
                                                      formatC(multi_obsvec_sex$p.optim.2, format = "e", digits = 2), ")")

multi_obsvec_infercmv = Multibin(infercmv_obsvec, vec_size, averagetotalreads)
ApplyRes_obsvec["Multi-bin-Equal", "infercmv"] = paste0(round(multi_obsvec_infercmv$T.param.equal/(10^4),3), " (",
                                                      formatC(multi_obsvec_infercmv$p.equal, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-SSW", "infercmv"] = paste0(round(multi_obsvec_infercmv$T.param.1/(10^4),3), " (",
                                                       formatC(multi_obsvec_infercmv$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-IVW", "infercmv"] = paste0(round(multi_obsvec_infercmv$T.param.2/(10^4),3), " (",
                                                           formatC(multi_obsvec_infercmv$p.optim.2, format = "e", digits = 2), ")")


multi_obsvec_cmv = Multibin(cmv_obsvec, vec_size, averagetotalreads)
ApplyRes_obsvec["Multi-bin-Equal", "cmv"] = paste0(round(multi_obsvec_cmv$T.param.equal/(10^4),3), " (",
                                                 formatC(multi_obsvec_cmv$p.equal, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-SSW", "cmv"] = paste0(round(multi_obsvec_cmv$T.param.1/(10^4),3), " (",
                                                  formatC(multi_obsvec_cmv$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-IVW", "cmv"] = paste0(round(multi_obsvec_cmv$T.param.2/(10^4),3), " (",
                                                      formatC(multi_obsvec_cmv$p.optim.2, format = "e", digits = 2), ")")

multi_obsvec_totalreads = Multibin(totalreads_obsvec, vec_size, averagetotalreads)
ApplyRes_obsvec["Multi-bin-Equal", "totalreads"] = paste0(round(multi_obsvec_totalreads$T.param.equal/(10^-2),3), " (",
                                                        formatC(multi_obsvec_totalreads$p.equal, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-SSW", "totalreads"] = paste0(round(multi_obsvec_totalreads$T.param.1/(10^-2),3), " (",
                                                         formatC(multi_obsvec_totalreads$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_obsvec["Multi-bin-IVW", "totalreads"] = paste0(round(multi_obsvec_totalreads$T.param.2/(10^-2),3), " (",
                                                             formatC(multi_obsvec_totalreads$p.optim.2, format = "e", digits = 2), ")")

##### Shannon Index 
vec_size = ApplyRes_shannon[paste0("Bin", 1:(length(ci)+1)), "samplesize"]

ApplyRes_shannon["Multi-bin-Equal","samplesize"] = sum(vec_size)
ApplyRes_shannon["Multi-bin-SSW","samplesize"] = sum(vec_size)
ApplyRes_shannon["Multi-bin-IVW","samplesize"] = sum(vec_size)

multi_shannon_age = Multibin(age_shannon, vec_size, averagetotalreads)
ApplyRes_shannon["Multi-bin-Equal", "age"] = paste0(round(multi_shannon_age$T.param.equal,3), " (",
                                                  formatC(multi_shannon_age$p.equal, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-SSW", "age"] = paste0(round(multi_shannon_age$T.param.1,3), " (",
                                                   formatC(multi_shannon_age$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-IVW", "age"] = paste0(round(multi_shannon_age$T.param.2,3), " (",
                                                       formatC(multi_shannon_age$p.optim.2, format = "e", digits = 2), ")")

multi_shannon_sex = Multibin(sex_shannon, vec_size, averagetotalreads)
ApplyRes_shannon["Multi-bin-Equal", "sex"] = paste0(round(multi_shannon_sex$T.param.equal,3), " (",
                                                  formatC(multi_shannon_sex$p.equal, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-SSW", "sex"] = paste0(round(multi_shannon_sex$T.param.1,3), " (",
                                                   formatC(multi_shannon_sex$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-IVW", "sex"] = paste0(round(multi_shannon_sex$T.param.2,3), " (",
                                                       formatC(multi_shannon_sex$p.optim.2, format = "e", digits = 2), ")")

multi_shannon_infercmv = Multibin(infercmv_shannon, vec_size, averagetotalreads)
ApplyRes_shannon["Multi-bin-Equal", "infercmv"] = paste0(round(multi_shannon_infercmv$T.param.equal,3), " (",
                                                       formatC(multi_shannon_infercmv$p.equal, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-SSW", "infercmv"] = paste0(round(multi_shannon_infercmv$T.param.1,3), " (",
                                                        formatC(multi_shannon_infercmv$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-IVW", "infercmv"] = paste0(round(multi_shannon_infercmv$T.param.2,3), " (",
                                                            formatC(multi_shannon_infercmv$p.optim.2, format = "e", digits = 2), ")")

multi_shannon_cmv = Multibin(cmv_shannon, vec_size, averagetotalreads)
ApplyRes_shannon["Multi-bin-Equal", "cmv"] = paste0(round(multi_shannon_cmv$T.param.equal,3), " (",
                                                  formatC(multi_shannon_cmv$p.equal, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-SSW", "cmv"] = paste0(round(multi_shannon_cmv$T.param.1,3), " (",
                                                   formatC(multi_shannon_cmv$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-IVW", "cmv"] = paste0(round(multi_shannon_cmv$T.param.2,3), " (",
                                                       formatC(multi_shannon_cmv$p.optim.2, format = "e", digits = 2), ")")

multi_shannon_totalreads = Multibin(totalreads_shannon, vec_size, averagetotalreads)
ApplyRes_shannon["Multi-bin-Equal", "totalreads"] = paste0(round(multi_shannon_totalreads$T.param.equal/(10^-7),3), " (",
                                                         formatC(multi_shannon_totalreads$p.equal, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-SSW", "totalreads"] = paste0(round(multi_shannon_totalreads$T.param.1/(10^-7),3), " (",
                                                          formatC(multi_shannon_totalreads$p.optim.1, format = "e", digits = 2), ")")
ApplyRes_shannon["Multi-bin-IVW", "totalreads"] = paste0(round(multi_shannon_totalreads$T.param.2/(10^-7),3), " (",
                                                              formatC(multi_shannon_totalreads$p.optim.2, format = "e", digits = 2), ")")


ApplyRes_obsvec
ApplyRes_shannon
