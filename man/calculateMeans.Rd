\name{calculateMeans}
\Rdversion{1.1}
\docType{methods}
\alias{calculateMeans}
\alias{calculateMeans,SparseDataSet-method}
\title{ Calculate the means and sum of squares for a SparseDataSet }
\description{ Calculate the means and sum of squares for a
  SparseDataSet.  The means and sum of squares from samples within a
  condition are calculated and then added as a named list of sparse
  columns (dgCMatrix objects).  An additional sparse column named
  'global' is added to each which is the mean of the condition means and
  the sum of the condition sums of squares.

  This function makes use of the \code{mclapply} function in the
  parallel package, to allow the user to distribute mean and sum of
  squares calculations across cores.
} 
\usage{
  \S4method{calculateMeans}{SparseDataSet}( object, nzr=1, quiet=FALSE,
recalc=FALSE, ... )
}
\arguments{
  \item{object}{ a SparseDataSet }
  \item{nzr}{ the desired nonzero ratio: the means and sum of squares
    can be thresholded by \code{\link{sparseThreshold}}, to achieve a
    certain percent of nonzero elements.  This might be useful as the
    row means of a sparse matrix can be much denser than the individual
    columns. } 
  \item{quiet}{ should the function not print out each condition name as
    it is calculated  }
  \item{recalc}{ should the function recalculate all means and sum of
    squares, or only those which are present in the conditions but
  missing from the means and sumSquares slots  }
  \item{...}{ extra arguments passed to \code{mclapply} of the parallel package }
}
\author{ Michael Love }
\examples{
  sds <- simulateSparseDataSet(5, c(5,5))
  sds <- calculateMeans(sds)
  means(sds)
}