simulateSparseDataSet <- function(n, samples.per.condition,nzg=.1,nzs=.1) {
  m <- sum(samples.per.condition)
  ncond <- length(samples.per.condition)
  cond <- factor(rep(paste("c",1:ncond,sep=""),samples.per.condition))
  mu0 <- ifelse(runif(n) < nzg, rgamma(n,10,1), 0)
  sparseData.list <- lapply(samples.per.condition, function(m_j) {
    mu <- ifelse(runif(n) < nzs, rgamma(n,10,1), 0)
    sample.sparseData <- as(matrix(rnbinom(n * m_j, mu=mu + mu0, size=1),ncol=m_j),"sparseMatrix")
  })
  sm <- do.call(cBind, sparseData.list)
  pd <- new("AnnotatedDataFrame",data.frame(sampleID=paste("sample",1:m,sep=""),row.names=paste("sample",1:m,sep=""),stringsAsFactors=FALSE))
  fd <- new("AnnotatedDataFrame",data.frame(featureID=paste("feat",1:n,sep=""),row.names=paste("feat",1:n,sep=""),stringsAsFactors=FALSE))
  sds <- newSparseDataSet(sm, cond, pd, fd)
  sds
}

# this code was developed from a discussion here:
# http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r
sparseCov <- function(X,Y=NULL,XtX=NULL) {
  if (!is(X,"dgCMatrix")) stop("X should be a dgCMatrix")
  if (is.null(Y)) {
    if (is.null(XtX)) {
      XtX <- crossprod(X)
    } else {
      if (ncol(XtX) != ncol(X) | nrow(XtX) != ncol(X)) stop("XtX should have same number of rows and columns as number of columns of X")
    }
    n <- nrow(X)
    cMeans <- colMeans(X)
    covmat <- (as.matrix(XtX) - n*tcrossprod(cMeans))/(n-1)
    sdvec <- sqrt(diag(covmat))
    cormat <- covmat/crossprod(t(sdvec))
    return(list(cov=covmat,cor=cormat))
  } else {
    if (!is(Y,"dgCMatrix")) stop("Y should be a dgCMatrix")
    if (nrow(X) != nrow(Y)) stop("X and Y should have the same number of rows")
    n <- nrow(X)
    cMeansX <- colMeans(X)
    cMeansY <- colMeans(Y)
    covmat <- (as.matrix(crossprod(X,Y)) - n * tcrossprod(cMeansX,cMeansY))/(n-1)
    sdvecX <- sqrt(diag((as.matrix(crossprod(X)) - n*tcrossprod(cMeansX))/(n-1)))
    sdvecY <- sqrt(diag((as.matrix(crossprod(Y)) - n*tcrossprod(cMeansY))/(n-1)))
    cormat <- covmat/outer(sdvecX,sdvecY)
    return(list(cov=covmat,cor=cormat))
  }
}

sparseCosine <- function(X,Y=NULL,XtX=NULL){
if (!is(X,"dgCMatrix")) stop("X should be a dgCMatrix")
  if (is.null(Y)) {
    if (is.null(XtX)) {
      XtX <- crossprod(X)
    } else {
      if (ncol(XtX) != ncol(X) | nrow(XtX) != ncol(X)) stop("XtX should have same number of rows and columns as number of columns of X")
    }
    lengths <- sqrt(diag(as.matrix(XtX)))
    return(as.matrix(XtX/crossprod(t(lengths))))
  } else {
    if (!is(Y,"dgCMatrix")) stop("Y should be a dgCMatrix")
    if (nrow(X) != nrow(Y)) stop("X and Y should have the same number of rows")
    lengths.X <- sqrt(diag(as.matrix(crossprod(X))))
    lengths.Y <- sqrt(diag(as.matrix(crossprod(Y))))
    XtY <- crossprod(X,Y)
    return(as.matrix(XtY/outer(lengths.X,lengths.Y)))
  }
}

sparseEuclid <- function(X,Y=NULL,XtX=NULL) {
if (!is(X,"dgCMatrix")) stop("X should be a dgCMatrix")
  if (is.null(Y)) {
    if (is.null(XtX)) {
      XtX <- crossprod(X)
    } else {
      if (ncol(XtX) != ncol(X) | nrow(XtX) != ncol(X)) stop("XtX should have same number of rows and columns as number of columns of X")
    }
    X.length <- colSums(X^2)
    length.sum <- outer(X.length,X.length,"+")
    dimnames(length.sum) <- dimnames(XtX)
    return(as.matrix(sqrt(length.sum - 2*XtX)))
  } else {
    if (!is(Y,"dgCMatrix")) stop("Y should be a dgCMatrix")
    if (nrow(X) != nrow(Y)) stop("X and Y should have the same number of rows")
    XtY <- crossprod(X,Y)
    X.length <- colSums(X^2)
    Y.length <- colSums(Y^2)
    length.sum <- outer(X.length,Y.length,"+")
    dimnames(length.sum) <- dimnames(XtY)
    return(as.matrix(sqrt(length.sum - 2*XtY)))
  }
}

applyFunctionSparsely <- function(z, f) {
  if (!(is(z,"dsparseMatrix") | is(z,"dsparseVector"))) stop("z should be a dsparseMatrix or dsparseVector")
  z@x <- f(z@x)
  z
}

sparseThreshold <- function(z, nzr = .1) {
  if (!(is(z,"dsparseMatrix") | is(z,"dsparseVector"))) stop("z should be a dsparseMatrix or dsparseVector")
  if (nzr < 0 | nzr > 1) stop("nzr should be a value between 0 and 1")
  nzratio <- nnzero(z)/length(z)
  if (nzratio > nzr) {
    threshold <- quantile(abs(z@x), 1 - nzr/nzratio)
    z[z != 0 & z <= threshold & z >= -threshold] <- 0
  }
  as(z,"sparseMatrix")
}

fillMatrixFromColumn <- function(col, n) {
  if (!(is(col,"dsparseMatrix") | is(col,"dsparseVector"))) stop("col should be a dsparseMatrix or dsparseVector")
  if (ncol(col) != 1) stop("col should have 1 column")
  x <- rep(col@x, n)
  # note: the i slot is 0-based
  i <- rep(col@i + 1, n)
  j <- rep(seq_len(n),each=length(col@x))
  C <- sparseMatrix(i=i,j=j,x=x,dims=c(length(col),n))
  return(C)
}
