\name{sparseCov}
\alias{sparseCov}
\alias{sparseCosine}
\alias{sparseEuclid}
\title{
  Sparse covariance/correlation, cosine similarity or Euclidean distance
}
\description{
  Calcultion of matrices of various similarity/distances, which takes use
  of matrix multiplication.  Note: the comparisons are made between
  columns, similar to the \code{cor} function, and unlike the
  \code{dist} function which compares between rows.

  See the timing vignette for comparison of the performance of functions
  against standard dense calculations. 
}
\usage{
sparseCov(X, Y = NULL, XtX = NULL)
sparseCosine(X, Y = NULL, XtX = NULL)
sparseEuclid(X, Y = NULL, XtX = NULL)
}
\arguments{
  \item{X}{
    a matrix X
  }
  \item{Y}{
    optional: a matrix Y, in which case columns of X will be compared to
    columns of Y.  The resulting matrix will be ncol(X) by ncol(Y).
  }
  \item{XtX}{
    optional: if not using Y, supplying XtX, the crossproduct of X with
    itself, will be used rather than recalculating this in the function.
  }
}
\value{
  A list with the covariance and correlation matrix of X, or X and Y;
  The cosine similarity matrix of X, or X and Y;
  The Euclidean distances matrix of X, or X and Y.
}
\references{
  The sparseCov function was developed from a discussion here: \url{http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r}
}
\author{
  Michael Love
}
\examples{
sds1 <- simulateSparseDataSet(100, c(5,5))
sds2 <- simulateSparseDataSet(100, c(2,2))
res <- sparseCov(sparseData(sds1))
res$cor
resxy <- sparseCov(sparseData(sds1),sparseData(sds2))
resxy$cor
}
