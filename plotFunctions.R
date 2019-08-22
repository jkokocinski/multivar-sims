#
# vAxisLims : function to compute good plotting limits for sample ACVF plots so
#   that the theoretical line shows up in the plotting range.
#
vAxisLims <- function(vec, val) {
  lims <- vector(mode="numeric",length=2)
  lims[1] <- min( c(min(vec), val) ) - 0.1*abs(diff(range(vec)))
  lims[2] <- max( c(max(vec), val) ) + 0.1*abs(diff(range(vec)))
  return(lims)
}
# end vAxisLims

#
# plotSqGridDims : Function to compute appropriate ylim ordered pair for
#   plotting sample ACVF values and the theretical value.
#
plotSqGridDims <- function(num) {
  if (sqrt(num)==trunc(sqrt(num))) {
    return(sqrt(num)*c(1,1))
  }
  else {
    dimVec <- c(floor(sqrt(num)),ceiling(sqrt(num)))
    while (prod(dimVec) < num) {
      inc1 <- dimVec + c(1,0)
      inc2 <- dimVec + c(0,1)
      if (abs(prod(inc1)-num)<abs(prod(inc2)-num)) {
        dimVec <- inc1
      } else {
        dimVec <- inc2
      }
    } # end while
    return(dimVec)
  } # end else
}
# end plotSqGridDims
