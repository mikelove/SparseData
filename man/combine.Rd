\name{combine-methods}
\Rdversion{1.1}
\docType{methods}
\alias{combine}
\alias{combine,SparseDataSet,SparseDataSet-method}
\alias{combine,dgCMatrix,dgCMatrix-method}
\title{ Combine two SparseDataSet objects }
\description{
  This method allows for two \code{\link{SparseDataSet-class}} objects
  to be combined, providing that they are defined over identical
  features and they have no shared sample names.  The function ensures
  that the levels of the conditions are set to the union of the levels
  of conditions of x and y.

  The means and sum of squares are carried over provided they are not in
  the intersection of the two condition levels (as these will need to be
  recalculated).  Also the t-statistics are removed as these will change
  with a new global mean and global sum of squares.

  The method defined for \code{dgCMatrix} objects combines two matrices
  using \code{cBind} 
}
\usage{
  \S4method{combine}{SparseDataSet,SparseDataSet}(x,y,...)
  \S4method{combine}{dgCMatrix,dgCMatrix}(x,y,...)
}
\arguments{
  \item{x}{ a SparseDataSet or dgCMatrix } 
  \item{y}{ a SparseDataSet or dgCMatrix }
  \item{...}{ additional arguments }
}
\value{ a SparseDataSet with combined sparseData, phenoData,
  featureData, experimentData and protocolData.
}
\author{ Michael Love }
\examples{
x <- simulateSparseDataSet(10,c(2,2,2))
y <- simulateSparseDataSet(10,c(2,2))
sampleNames(y) <- paste("sample",(ncol(x) + 1:ncol(y)),sep="")
pData(y)$sampleID <- sampleNames(y)
z <- combine(x,y)
pData(z)
sparseData(z)
}