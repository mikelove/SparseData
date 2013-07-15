### R code from vignette source 'SparseData.Rnw'

###################################################
### code chunk number 1: SparseData.Rnw:34-46
###################################################
library(SparseData)
sparse.data.files <- list.files(system.file("extdata",package="SparseData"), 
                                "counts", full=TRUE)
sparse.data.names <- list.files(system.file("extdata",package="SparseData"), 
                                "counts")
feature.file <- list.files(system.file("extdata",package="SparseData"), 
                           "ranges", full=TRUE)
feature.file.name <- list.files(system.file("extdata",package="SparseData"), 
                                "ranges")
# these filenames look like:
sparse.data.names[1]
feature.file.name


###################################################
### code chunk number 2: SparseData.Rnw:64-71
###################################################
sparse.data.conditions <- sub(".+(Fetal.+)\\.counts","\\1",sparse.data.names)
phenoData <- AnnotatedDataFrame(data.frame(filename=sparse.data.names, 
                                           conditions=sparse.data.conditions))
feature.data.frame <- read.delim(feature.file, header=FALSE)
featureData <- AnnotatedDataFrame(data.frame(chr=feature.data.frame[,1], 
                                             start=feature.data.frame[,2], 
                                             end=feature.data.frame[,3]))


###################################################
### code chunk number 3: SparseData.Rnw:76-86
###################################################
# threshold incoming count files at 50%
sparse.data.list <- lapply(sparse.data.files, function(filename) {
  sparseThreshold(Matrix(scan(filename,quiet=TRUE), sparse=TRUE), nzr=.5)
})
# bind the sparse columns together as a sparse matrix
sparse.data <- do.call(cBind, sparse.data.list)
sds <- newSparseDataSet(sparseData = sparse.data, 
                        conditions=phenoData$conditions, 
                        featureData=featureData, 
                        phenoData=phenoData)


###################################################
### code chunk number 4: SparseData.Rnw:91-96
###################################################
logPlusOne <- function(x) log(x + 1)
sparseData(sds) <- applyFunctionSparsely(sparseData(sds), logPlusOne)
norm.mat <- Matrix(diag(1/colMeans(sparseData(sds))),sparse=TRUE)
colnames(norm.mat) <- rownames(pData(sds))
sparseData(sds) <- sparseData(sds) %*% norm.mat


###################################################
### code chunk number 5: SparseData.Rnw:101-104
###################################################
options(mc.cores=1)
sds <- calculateMeans(sds)
sds <- calculateTStats(sds)


###################################################
### code chunk number 6: SparseData.Rnw:109-113
###################################################
# the top five features specific to Fetal_Brain
fData(sds)[head(order(-tStats(sds)[["Fetal_Brain"]]),5),]
# the top five global features
fData(sds)[head(order(-means(sds)[["global"]]),5),]


###################################################
### code chunk number 7: SparseData.Rnw:120-134
###################################################
library(SparseData)
sparse.data.files <- list.files(system.file("extdata",package="SparseData"), 
                                "counts", full=TRUE)
sparse.data.names <- list.files(system.file("extdata",package="SparseData"), 
                                "counts")
sparse.data.conditions <- sub(".+(Fetal.+)\\.counts","\\1",sparse.data.names)
phenoData <- AnnotatedDataFrame(data.frame(filename=sparse.data.names, 
                                           conditions=sparse.data.conditions))
feature.file <- list.files(system.file("extdata",package="SparseData"), 
                           "ranges", full=TRUE)
feature.data.frame <- read.delim(feature.file, header=FALSE)
featureData <- AnnotatedDataFrame(data.frame(chr=feature.data.frame[,1], 
                                             start=feature.data.frame[,2], 
                                             end=feature.data.frame[,3]))


###################################################
### code chunk number 8: SparseData.Rnw:150-155
###################################################
quantile(scan(sparse.data.files[1],quiet=TRUE),0:10/10)
sparse.data.list <- lapply(sparse.data.files, function(filename) {
  sparseThreshold(Matrix(scan(filename,quiet=TRUE), sparse=TRUE), nzr=.5)
})
sparse.data <- do.call(cBind, sparse.data.list)


###################################################
### code chunk number 9: SparseData.Rnw:162-166
###################################################
sds <- newSparseDataSet(sparseData = sparse.data, 
                        conditions=phenoData$conditions, 
                        featureData=featureData, 
                        phenoData=phenoData)


###################################################
### code chunk number 10: SparseData.Rnw:171-179
###################################################
expData <- new("MIAME",
name="2000 ranges of DNase-seq from Roadmap Epigenome Mapping Consortium",
lab="University of Washington",
contact="rharris1@bcm.tmc.edu",
title="Human Reference Epigenome Mapping Project",
url="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18927")
pubMedIds(expData) <- "20944595"
experimentData(sds) <- expData


###################################################
### code chunk number 11: SparseData.Rnw:184-187
###################################################
head(sparseData(sds),10)
head(pData(sds),3)
head(fData(sds),3)


###################################################
### code chunk number 12: SparseData.Rnw:192-193
###################################################
nnzero(sparseData(sds))/prod(dim(sds))


###################################################
### code chunk number 13: SparseData.Rnw:198-203
###################################################
logPlusOne <- function(x) log(x + 1)
sparseData(sds) <- applyFunctionSparsely(sparseData(sds), logPlusOne)
norm.mat <- Matrix(diag(1/colMeans(sparseData(sds))),sparse=TRUE)
colnames(norm.mat) <- rownames(pData(sds))
sparseData(sds) <- sparseData(sds) %*% norm.mat


###################################################
### code chunk number 14: SparseData.Rnw:210-213
###################################################
options(mc.cores=1)
sds <- calculateMeans(sds)
sds <- calculateTStats(sds)


###################################################
### code chunk number 15: SparseData.Rnw:218-222
###################################################
names(means(sds))
means(sds)[["Fetal_Brain"]] <- NULL
sds <- calculateMeans(sds)
names(means(sds))


###################################################
### code chunk number 16: SparseData.Rnw:227-228
###################################################
fData(sds)[head(order(-tStats(sds)[["Fetal_Brain"]]),5),]


###################################################
### code chunk number 17: SparseData.Rnw:234-239
###################################################
par(mfrow=c(3,1),mar=c(1,1,3,1))
for (cond in levels(conditions(sds))) {
  top.feats <- head(order(-tStats(sds)[[cond]]),5)
  image(t(as.matrix(sparseData(sds)[top.feats, order(conditions(sds))])), main=cond, col=grey(10:0/10), xaxt="n", yaxt="n")
}


###################################################
### code chunk number 18: SparseData.Rnw:246-252
###################################################
par(mfrow=c(3,1),mar=c(1,1,3,1))
for (cond in levels(conditions(sds))) {
  bottom.feats <- head(order(tStats(sds)[[cond]]),5)
  image(t(as.matrix(sparseData(sds)[bottom.feats, order(conditions(sds))])), 
        main=cond, col=grey(10:0/10), xaxt="n", yaxt="n")
}


###################################################
### code chunk number 19: SparseData.Rnw:260-261
###################################################
cormat <- sparseCov(sparseData(sds))$cor


###################################################
### code chunk number 20: SparseData.Rnw:265-269
###################################################
par(mar=c(1,1,3,1))
image(t(cormat[order(conditions(sds)),rev(order(conditions(sds)))]), 
      col=colorRampPalette(c("red","white","green"))(21), 
      zlim=c(-1,1), xaxt="n", yaxt="n", main="Correlation matrix of data")


###################################################
### code chunk number 21: SparseData.Rnw:275-278
###################################################
means.matrix <- do.call(cBind, means(sds)[match(
  levels(pData(sds)$condition), names(means(sds)))])
match.cormat <- sparseCov(sparseData(sds), means.matrix)$cor


###################################################
### code chunk number 22: SparseData.Rnw:282-286
###################################################
par(mar=c(1,1,3,1))
image(t(match.cormat[order(conditions(sds)),ncol(match.cormat):1]), 
      col=colorRampPalette(c("red","white","green"))(21), 
      zlim=c(-1,1), xaxt="n", yaxt="n", main="Correlation of data to condition means")


###################################################
### code chunk number 23: SparseData.Rnw:313-323
###################################################
sim.sds <- simulateSparseDataSet(200, c(50, 50, 50), nzg = 0.5, nzs = 0.5)
sim.sds <- calculateMeans(sim.sds, quiet=TRUE)
sim.sds <- calculateTStats(sim.sds, quiet=TRUE)
equalvar.t.stats <- sapply(1:nrow(sim.sds), function(i) t.test(
  x=sparseData(sim.sds)[i,pData(sim.sds)$condition == "c1"], 
  y=sparseData(sim.sds)[i,], var.equal=TRUE)$statistic)
plot(equalvar.t.stats, tStats(sim.sds)[["c1"]], 
     xlab="equal-variance t-statistics",ylab="tStats",
     main="tStats vs. equal-variance t-statistics")
abline(0,1)


###################################################
### code chunk number 24: SparseData.Rnw:333-340
###################################################
x <- simulateSparseDataSet(10,c(2,2,2))
y <- simulateSparseDataSet(10,c(2,2))
sampleNames(y) <- paste("sample",(ncol(x) + 1:ncol(y)),sep="")
pData(y)$sampleID <- sampleNames(y)
z <- combine(x,y)
pData(z)
sparseData(z)


###################################################
### code chunk number 25: SparseData.Rnw:345-356
###################################################
x <- calculateMeans(x)
x <- calculateTStats(x)
y <- calculateMeans(y)
y <- calculateTStats(y)
z <- combine(x,y)
names(means(z))
names(tStats(z))
z <- calculateMeans(z)
z <- calculateTStats(z)
names(means(z))
names(tStats(z))


###################################################
### code chunk number 26: SparseData.Rnw:361-362
###################################################
sessionInfo()


