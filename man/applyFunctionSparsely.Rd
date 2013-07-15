\name{applyFunctionSparsely}
\alias{applyFunctionSparsely}
\title{ Apply function sparsely }
\description{
  Applies a function to the nonzero elements of a sparse vector or matrix
}
\usage{
  applyFunctionSparsely(z, f)
}
\arguments{
  \item{z}{ a sparse vector or matrix }
  \item{f}{ a function }
}
\value{
  Returns the variable z, with only the value of the nonzero elements changed
}
\author{
  Michael Love
}
\examples{
sds <- simulateSparseDataSet(100, c(5,5))
logPlusOne <- function(x) log(x + 1)
sparseData(sds) <- applyFunctionSparsely(sparseData(sds), logPlusOne)
}
