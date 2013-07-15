\name{newSparseDataSet}
\alias{newSparseDataSet}
\title{
Create a new SparseDataSet
}
\description{
  Generate a new SparseDataSet from annotation data and a sparse
  matrix.  For more details on the class see
  \code{\link{SparseDataSet-class}}. 
  
  For construction of phenoData and featureData objects see the
  introductory vignette or the ExpressionSet introdcutory vignette of the
  Biobase package.
  
  The sparse data can be constructed using the \code{Matrix} function of
  the Matrix package, see the introductory vignette or check the man files
  from the Matrix package for more information on this class.
  
  For calculation of statistics see \code{\link{calculateMeans}} and
  \code{\link{calculateTStats}}.
}
\usage{
newSparseDataSet(sparseData, conditions, phenoData = NULL, featureData = NULL)
}
\arguments{
  \item{sparseData}{
    an object of class dgCMatrix
  }
  \item{conditions}{
    a factor of the conditions across which statistics will be calculated
  }
  \item{phenoData}{
    the phenotypic data
  }
  \item{featureData}{
    the feature data
  }
}
\value{
a SparseDataSet object
}
\author{
Michael Love
}
\examples{
sds <- simulateSparseDataSet(5, c(5,5))
sds2 <- newSparseDataSet(sparseData(sds), conditions(sds), phenoData(sds), featureData(sds))
}
