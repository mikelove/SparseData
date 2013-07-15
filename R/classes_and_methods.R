setClass("SparseDataSet", contains = "eSet",
         representation = representation( 
           means = "list",
           sumSquares = "list",
           tStats = "list"))

newSparseDataSet <- function( sparseData, conditions, phenoData = NULL, featureData = NULL ) {
  phenoData$condition <- factor(conditions)
  sds <- new( "SparseDataSet",
             assayData = assayDataNew( "environment", sparseData=sparseData ),
             phenoData = phenoData, 
             featureData = featureData)
  sds
}

setValidity( "SparseDataSet", function( object ) {
   if( !is(sparseData(object),"dgCMatrix") )
      return( "sparseData slot must be a dgCMatrix" )
   TRUE
} )

setMethod("combine",
          signature=signature(
            x="SparseDataSet", y="SparseDataSet"),
          function(x, y, ...) {
            # this method begins and ends with copied code from combine method in Biobase
            if (class(x) != class(y))
              stop(paste("objects must be the same class, but are ",
                         class(x), ", ", class(y), sep=""))
            if (any(annotation(x) != annotation(y)))
              stop("objects have different annotations: ",
                   annotation(x), ", ", annotation(y))
            
            # here we check to see that the sample names do not overlap
            sharedNames <- intersect(sampleNames(x), sampleNames(y))
            if (length(sharedNames) > 0) {
              stop("shared sample names should be eliminated before calling combine")
            }
            # here we check to see that the feature names are identical
            if (!all(featureNames(x) == featureNames(y))) {
              stop("all feature names should be identical")
            }
            
            # here we ensure that the conditions are defined over the union of levels
            conditions(x) <- factor(conditions(x), levels=union(conditions(x),conditions(y)))
            conditions(y) <- factor(conditions(y), levels=union(conditions(x),conditions(y)))
            
            # here we can concatenate any already calculated statistics
            # which are for unshared conditions.  the statistics for shared
            # conditions need to be recalculated.
            
            sharedConditions <- intersect(conditions(x),conditions(y))
            means(x) <- c(means(x)[names(means(x)) %in% setdiff(levels(conditions(x)),sharedConditions)],
                          means(y)[names(means(y)) %in% setdiff(levels(conditions(y)),sharedConditions)])
            sumSquares(x) <- c(sumSquares(x)[names(sumSquares(x)) %in% setdiff(levels(conditions(x)),sharedConditions)],
                          sumSquares(y)[names(sumSquares(y)) %in% setdiff(levels(conditions(y)),sharedConditions)])

            # we remove t-statistics because these need to be
            # recalculated given a new global mean and global sum of squares
            tStats(x) <- list()
            
            # original code from combine method in Biobase
            assayData(x) <- combine(assayData(x), assayData(y))
            phenoData(x) <- combine(phenoData(x), phenoData(y))
            featureData(x) <- combine(featureData(x), featureData(y))
            experimentData(x) <- combine(experimentData(x),experimentData(y))
            protocolData(x) <- combine(protocolData(x), protocolData(y))
            ## annotation -- constant
            x
          })

setMethod("combine", signature(x="dgCMatrix",y="dgCMatrix"),
          function(x, y, ...) {
            # altered code from combine method in BiocGenerics
            if (length(y) == 0L)
              return(x)
            else if (length(x) == 0L)
              return(y)
            if (mode(x) != mode(y))
              stop("matrix modes ", mode(x), ", ", mode(y), " differ")
            if (typeof(x) != typeof(y))
              warning("matrix typeof ", typeof(x), ", ", typeof(y), " differ")  
            dgcm <- cBind(x,y)
            dgcm
          })

sparseData <- function(object) standardGeneric("sparseData")
setGeneric("sparseData", sparseData)
setMethod("sparseData", signature(object="SparseDataSet"),
          function(object) assayData(object)[["sparseData"]])


setGeneric("sparseData<-", function(object, value) standardGeneric("sparseData<-"))
setReplaceMethod("sparseData", signature(object="SparseDataSet", value="dgCMatrix"),
  function( object, value ) {
   assayData(object)[[ "sparseData" ]] <- value
   validObject(object)
   object
})

setMethod("conditions", signature(object="SparseDataSet"),
          function(object) pData(object)$condition)

setReplaceMethod("conditions", signature(object="SparseDataSet", value="factor"),
  function( object, value ) {
   pData(object)$"condition" <- value
   validObject(object)
   object
})

means <- function(object) standardGeneric("means")
setGeneric("means", means)
setMethod("means", signature(object="SparseDataSet"), function(object) object@means)

setGeneric("means<-", function(object, value) standardGeneric("means<-"))
setReplaceMethod("means", signature(object="SparseDataSet"), function(object, value) {
  object@means <- value
  object
})

calculateMeans <- function(object, nzr=1, quiet=FALSE, recalc=FALSE, ...) standardGeneric("calculateMeans")
setGeneric("calculateMeans", calculateMeans)
setMethod("calculateMeans", signature(object="SparseDataSet"),
          function(object, nzr=1, quiet=FALSE,...) {
            # recalculate statistics for only conditions not in the means and sumSquares list (default)
            # or recalculate all statistics
             if (!recalc) {
              conditions.for.stats <- setdiff(levels(conditions(object)), intersect(names(means(object)),names(sumSquares(object))))
              # if there are some unpaired means and sum of squares for the conditions
              # which are being recalculated, we remove them here
              means(object)[names(means(object)) %in% conditions.for.stats] <- NULL
              sumSquares(object)[names(sumSquares(object)) %in% conditions.for.stats] <- NULL
              if (length(conditions.for.stats) == 0) {
                message("conditions for all statistics have been calculated")
                return(object)
              }
            } else {
              means(object) <- list()
              sumSquares(object) <- list()
              conditions.for.stats <- levels(conditions(object))
            }
             object.stats <- mclapply(conditions.for.stats, function(cond) {
               if (!quiet) message(paste(cond," ",sep=""))
               sparseData.subset <- sparseData(object)[,conditions(object)==cond,drop=FALSE]
               sparseData.mean <- sparseThreshold(rowMeans(sparseData.subset,sparseResult=TRUE),nzr)
               sparseData.sum.squares <- sparseThreshold(rowSums((sparseData.subset - fillMatrixFromColumn(sparseData.mean, sum(conditions(object)==cond)))^2,sparseResult=TRUE),nzr)
               list(mean=sparseData.mean, sum.squares=sparseData.sum.squares)
             },...)
             # extract the means and sum of squares from the object.stats list
             # to construct two lists, append these to the existing lists
             object.means <- lapply(object.stats, function(x) x[["mean"]])
             names(object.means) <- conditions.for.stats
             means(object) <- c(means(object), object.means)
             object.sum.squares <- lapply(object.stats, function(x) x[["sum.squares"]])
             names(object.sum.squares) <- conditions.for.stats
             sumSquares(object) <- c(sumSquares(object), object.sum.squares)
             if (!quiet) message("global")
             means(object)[["global"]] <- sparseThreshold(rowMeans(do.call(cBind, object.means),sparseResult=TRUE),nzr)
             sumSquares(object)[["global"]] <- as(rowSums(do.call(cBind, object.sum.squares),sparseResult=TRUE),"sparseMatrix")
             object
           })

sumSquares <- function(object) standardGeneric("sumSquares")
setGeneric("sumSquares", sumSquares)
setMethod("sumSquares", signature(object="SparseDataSet"), function(object) object@sumSquares)

setGeneric("sumSquares<-", function(object, value) standardGeneric("sumSquares<-"))
setReplaceMethod("sumSquares", signature(object="SparseDataSet"), function(object, value) {
  object@sumSquares <- value
  object
})

tStats <- function(object) standardGeneric("tStats")
setGeneric("tStats", tStats)
setMethod("tStats", signature(object="SparseDataSet"), function(object) object@tStats)

setGeneric("tStats<-", function(object, value) standardGeneric("tStats<-"))
setReplaceMethod("tStats", signature(object="SparseDataSet"), function(object, value) {
  object@tStats <- value
  object
})

calculateTStats <- function(object, offset="mean", quiet=FALSE) standardGeneric("calculateTStats")
setGeneric("calculateTStats", calculateTStats)
setMethod("calculateTStats", signature(object="SparseDataSet"),
          function(object, offset="mean", quiet=FALSE) {
            if (!all(c(levels(conditions),"global") %in% names(means(object))) | !all(c(levels(conditions),"global") %in% names(sumSquares(object)))) {
              stop("not all means or sum of squares have been calculated, use calculateMeans")
            }
            if (!(offset == "mean" | is.numeric(offset))) stop("offset should be the character string 'mean' or a numeric value")
            # using the notation of Tibshirani, Hastie, Narasimhan, & Chu:
            # Diagnosis of multiple cancer types by shrunken centroids of gene expression  2002
            n <- ncol(object)
            K <- nlevels(conditions(object))
            pooled.within.class.sd <- applyFunctionSparsely(sumSquares(object)[["global"]]/(n - K),sqrt)
            # use mean for an offset rather than median because of the zeros
            if (offset == "mean") {
              s_0 <- mean(pooled.within.class.sd)
            } else {
              s_0 <- offset
            }
            object.t.stats <- lapply(levels(conditions(object)), function(cond) {
              if (!quiet) message(paste(cond," ",sep=""))
              n_k <- sum(conditions(object)==cond)
              m_k <- sqrt(1/n + 1/n_k)
              as((means(object)[[cond]] - means(object)[["global"]])/(m_k * (pooled.within.class.sd + s_0)),"sparseMatrix")
            })
            names(object.t.stats) <- levels(conditions(object))
            tStats(object) <- object.t.stats
            object
          })






