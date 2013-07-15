\name{SparseDataSet-class}
\Rdversion{1.1}
\docType{class}
\alias{SparseDataSet-class}
\title{Class \code{"SparseDataSet"}}
\description{
  An extension of the eSet class for sparse matrix data.  The assay data
  is stored as a dgCMatrix and accessible via the
  \code{\link{sparseData}} function.

  The \code{dgCMatrix-class} is compressed, sparse, column-oriented
  format.  This class is convenient because the typical high-throughput
  assay data is a tall matrix: perhaps hundreds of samples defined over
  thousands or millions of features.  Distances and correlations between
  samples can be calculated quickly using matrix multiplication.
}
\section{Objects from the Class}{
  Objects can be created with calls to \code{\link{newSparseDataSet}}.
}
\section{Slots}{
  \describe{
    \item{\code{means}:}{Object of class \code{"list"} ~~ }
    \item{\code{sumSquares}:}{Object of class \code{"list"} ~~ }
    \item{\code{tStats}:}{Object of class \code{"list"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{sparseData}{\code{signature(object = "SparseDataSet")}: ... }
    \item{conditions}{\code{signature(object = "SparseDataSet")}: ... }
    \item{combine}{\code{signature(object = "SparseDataSet"), signature(object = "SparseDataSet")}: ... }
    \item{calculateMeans}{\code{signature(object = "SparseDataSet")}: ... }
    \item{means}{\code{signature(object = "SparseDataSet")}: ... }
    \item{means<-}{\code{signature(object = "SparseDataSet")}: ... }
    \item{sumSquares}{\code{signature(object = "SparseDataSet")}: ... }
    \item{sumSquares<-}{\code{signature(object = "SparseDataSet")}: ... }
    \item{calculateTStats}{\code{signature(object = "SparseDataSet")}: ... }
    \item{tStats}{\code{signature(object = "SparseDataSet")}: ... }
    \item{tStats<-}{\code{signature(object = "SparseDataSet")}: ... }
  }
}
\author{
  Michael Love
}
