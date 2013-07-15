\name{fillMatrixFromColumn}
\alias{fillMatrixFromColumn}
\title{
  fillMatrixFromColumn
}
\description{
  The fillMatrixFromColumn takes a single column and an integer for the
  number of times to repeat this column. 
}
\usage{
fillMatrixFromColumn(col, n)
}
\arguments{
  \item{col}{ a sparse column matrix}
  \item{n}{ integer, how many times to repeat the column }
}
\value{
A sparse matrix
}
\author{
Michael Love
}
\examples{
  sds <- simulateSparseDataSet(5, c(5,5))
  sds <- calculateMeans(sds)
  fillMatrixFromColumn(means(sds)[["c1"]], 3)
}
