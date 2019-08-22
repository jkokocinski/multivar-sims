toepLeftMult <- function(toepVec,rightVec) {
  # toepVec should be such that the first N entries correspond with the last row
  #   of the Toeplitz matrix (this probably means th reversed ACVF vector).
  stopifnot(class(toepVec)=="numeric", class(rightVec)=="numeric")
  stopifnot(length(toepVec)==2*length(rightVec)-1)
  N <- length(rightVec)
  resVec <- matrix(nrow=N, ncol=1)
  for (j in 1:N) {
    resVec[j,1] <- rbind(toepVec[(N-j+1):(2*N-j)]) %*% cbind(rightVec)
  }
  return(resVec)
}


toepLeftMult2 <- function(toepVec,rightVec) {
  # Performs the matrix multiplication by doing a multiplkication in the
  #   frequency domain, instead of convolving in the time domain, as in
  #   toepLeftMult.
  # toepVec should be such that the first N entries of rev(toepVec) correspond
  #   to the last row of the Toeplitz matrix (just the regular ACVF vector).
  if (class(toepVec)!="numeric") { stop("toepVec must be a numeric vector.") }
  if (class(rightVec)!="numeric") { stop("rightVec must be a numeric vector.") }
  if (length(toepVec)!=2*length(rightVec)-1) {
    stop("lengths of args don't match")
  }
  N <- length(rightVec)
  resVec <- matrix(nrow=N, ncol=1)
  fft.tv <- fft(z=toepVec)
  fft.rv <- fft(z=c(rightVec, rep(0, N-1)))
  fft.res <- as.vector(fft.tv*fft.rv)
  resVec <- Re(fft(z=fft.res, inverse=TRUE)) / (2*N-1)
  resVec <- resVec[seq(N,2*N-1)]
  return(resVec)
}
