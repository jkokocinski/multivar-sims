#
# psumCov : function to compute sum of A^k %*% V %*% t(A)^k, k from 0 to infty,
#   truncated such that successive term difference is less than epsilon.
#
psumCov <- function(A, V=diag(1,dim(A)[1])) {
  if (dim(A)[1]!=dim(A)[2]) stop("Not a square matrix.")
  psum <- ident <- diag(1,dim(A)[1])
  epsilon <- 1
  iter <- 1
  while (epsilon>1e-16) { # threshold for convergence
    psum.old <- psum
    psum <- A %*% psum %*% t(A) + ident
    epsilon <- max(abs(psum-psum.old))
    iter <- iter+1L
  }
  cat(paste0("Series converged after ", iter, " terms.\n"))
  return(psum)
}
# end psumCov
