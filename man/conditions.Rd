\name{conditions}
\Rdversion{1.1}
\docType{methods}
\alias{conditions}
\alias{conditions,SparseDataSet-method}
\alias{conditions<-,SparseDataSet-method}
\alias{conditions<-,SparseDataSet,factor-method}
\title{ Accessor for conditions }
\description{ Accessor for conditions, as defined for the
  \code{CountDataSet-class}.  A factor which determines the
  conditions across which means, sum of squares and t-statistics will be
  calculated. }
\usage{
  \S4method{conditions}{SparseDataSet}(object)
  \S4method{conditions}{SparseDataSet}(object) <- value
}
\arguments{
  \item{object}{ a SparseDataSet}
  \item{value}{a factor with length equal to the columns of sparseData
    and the rows of phenoData. }
}
\author{ Michael Love }
\examples{
  sds <- simulateSparseDataSet(5, c(5,5))
  conditions(sds)
}