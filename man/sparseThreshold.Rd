\name{sparseThreshold}
\alias{sparseThreshold}
\title{ sparseThreshold }
\description{
  Pushes the smaller values (in absolute value) to zero to achieve a
  certain nonzero ratio.  In other words, the nonzero values closest to
  zero are assigned zero to achieve at least the specified nonzero ratio.
}
\usage{
  sparseThreshold(z, nzr=.1)
}
\arguments{
  \item{z}{ a sparse vector or matrix }
  \item{nzr}{ the desired nonzero ratio }
}
\value{
  Returns the variable z, with more or equal nonzero elements to achieve
  at least a nonzero ratio of nzr
}
\author{
  Michael Love
}
\examples{
sds <- simulateSparseDataSet(100, c(5,5), nzs=.5, nzg=.3)
nnzero(sparseData(sds))/prod(dim(sds))
sparseData(sds) <- sparseThreshold(sparseData(sds), .1)
nnzero(sparseData(sds))/prod(dim(sds))
}
