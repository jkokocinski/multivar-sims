#!/usr/bin/env R
library(mAr)
library(multitaper)
library(expm)

numObs <- 2**10

############################ bivariate AR(1) diag A ############################
#
# Bivariate AR(1) with diagonal A (matrix coefficient), so not correlated across
# series.
#   X_{1,t} = 0.7 X_{1,t-1}
#   X_{2,t} = 0.2 X_{2,t-1}
#
phiMat <- diag(c(0.7,0.2))
errCovMat <- diag(1,2,2) # no correlation between components of noise process

X <- mAr.sim(w=rep(0,2), A=phiMat, C=errCovMat, N=numObs)

# calculate covMats, the theoretical covariance matrices, at lags 0, 1, 2
covMats <- list()
baseMat <- matrix(c((1-phiMat[1,1]^2)^(-1), (1-phiMat[1,1]*phiMat[2,2])^(-1),
                    (1-phiMat[1,1]*phiMat[2,2])^(-1), (1-phiMat[2,2]^2)^(-1)),
                  nrow=2, ncol=2, byrow=TRUE) * errCovMat
for (h in 0:3) {
  mtx <- diag(c(phiMat[1,1]^h,phiMat[2,2]^h)) %*% baseMat
  covMats[[h+1]] <- mtx
}
cat("Theoretical covariance matrices:\n")
covMats

# theoretical lag-1 autocovariance function
phiMat %*% solve(diag(1,2,2)-phiMat%*%phiMat)

# sample ACVFs and CCVFs
theACVF1 <- acf(x=X$X1, type="covariance", plot=FALSE)
theACVF2 <- acf(x=X$X2, type="covariance", plot=FALSE)
theCCVF12 <- ccf(x=X$X1, y=X$X2, type="covariance", plot=FALSE)
theCCVF21 <- ccf(x=X$X2, y=X$X1, type="covariance", plot=FALSE)
# end bivariate AR(1) diag A




############################### bivariate AR corr ##############################
#
# Bivariate AR(1), correlated across series.
#
source("psumCov.R") # import psumCov function

phiMat2 <- matrix(c(0.7,-0.5,0.6,0.2), nrow=2, ncol=2, byrow=TRUE)
errCovMat2 <- diag(1,2,2)

bivAR1 <- mAr.sim(w=rep(0,2), A=phiMat2, C=errCovMat2, N=numObs)

bivAR1.ACVF1 <- acf(bivAR1$X1, type="cov", plot=FALSE)
bivAR1.ACVF2 <- acf(bivAR1$X2, type="cov", plot=FALSE)
bivAR1.CCVF12 <- ccf(x=bivAR1$X1, y=bivAR1$X2, type="cov", plot=FALSE)
bivAR1.CCVF21 <- ccf(x=bivAR1$X2, y=bivAR1$X1, type="cov", plot=FALSE)

# compute theoretical cov matrices
matrixSeries <- psumCov(phiMat2)

maxlag <- 50
CCVmats <- matrix(rep(matrixSeries,2*maxlag+1), nrow=2, ncol=2*(2*maxlag+1))
for (h in seq(1,maxlag,1)) {
  startCol <- 2*maxlag+1 + 2*h # left-most column in matrix for lag +h
  CCVmats[,seq(startCol,startCol+1)] <- phiMat2 %*%
    CCVmats[,seq(startCol-2,startCol-1)]
  
  startCol <- dim(CCVmats)[2]-startCol # left-most column in matrix for lag -h
  CCVmats[,seq(startCol,startCol+1)] <- CCVmats[,seq(startCol+2,startCol+3)] %*%
    t(phiMat2)
}

# plots of sample ACCVFs versus theoretical values
par(mar=c(4,4,3,1), mgp=c(2.5, 1, 0))
plot(bivAR1.CCVF12, main="Auto-Cross-Covariance (X[1],X[2])",
     xlim=(maxlag+ceiling(maxlag/5))*c(-1,1))
points(x=(-maxlag:maxlag), y=CCVmats[1,2*seq(0,2*maxlag)+2], col="blue")
dev.off()

par(mar=c(4,4,3,1), mgp=c(2.5, 1, 0))
plot(bivAR1.CCVF21, main="Auto-Cross-Covariance (X[2],X[1])",
     xlim=(maxlag+ceiling(maxlag/5))*c(-1,1))
points(x=(-maxlag:maxlag), y=CCVmats[2,2*seq(0,2*maxlag)+1], col="blue")
dev.off()


# plots of sample ACVFs versus theoretical values
par(mar=c(4,4,3,1), mgp=c(2.5, 1, 0))
plot(bivAR1.ACVF1, main="Autocovariance (X[1])",
     xlim=(maxlag+ceiling(maxlag/3))*c(0,1))
points(x=(0:maxlag), y=CCVmats[1,2*seq(maxlag,2*maxlag)+1], col="blue")
dev.off()

par(mar=c(4,4,3,1), mgp=c(2.5, 1, 0))
plot(bivAR1.ACVF2, main="Autocovariance (X[2])",
     xlim=(maxlag+ceiling(maxlag/3))*c(0,1))
points(x=(0:maxlag), y=CCVmats[2,2*seq(maxlag,2*maxlag)+2], col="blue")
dev.off()
# end bivariate AR corr



######################### ACCVFs sample vs. theoretical ########################
#
# Compare sample ACCVFs to theoretical values, varying the sample size for the
#   simulations. The same process as above is used, defined by phiMat2 and
#   errCovMat2 (correlated across series).
#
source("plotFunctions.R") # import vAxisLims, plotSqGridDims

numObsVec <- 1e3*seq(1,32)
theLags <- seq(0,3)
smplVals <- data.frame(matrix(nrow=length(numObsVec), ncol=1+length(theLags)))
names(smplVals) <- c("N",paste0("lag",theLags))
smplVals$N <- numObsVec
smplVals <- list(smplVals, smplVals, smplVals, smplVals)

cat(paste0("######## start time","\t",Sys.time(),"\n"))
NUM_RLZNS <- 10
for (rlzn in 1:NUM_RLZNS) {
  if (rlzn==1L) { # on 1st iteration, prev.df is the 0 dataframe
    prev.df <- list(
      data.frame(
        matrix(0,nrow=length(numObsVec),ncol=1+length(theLags))
      )
    )
    prev.df[[2]] <- prev.df[[3]] <- prev.df[[4]] <- prev.df[[1]]
  } else {
    for (k in 1:4) {
      prev.df[[k]] <- smplVals[[k]]
    }
  }
  simAR <- mAr.sim(w=rep(0,2), A=phiMat2, C=errCovMat2, N=max(numObsVec))
  for (n in numObsVec) {
    simAR.ACVest <- acf(simAR$X1[1:n], lag.max=max(theLags), type="cov",
                        plot=FALSE)
    smplVals[[1]][which(smplVals[[1]]$N==n),paste0("lag",theLags)] <-
      simAR.ACVest$acf
    
    simAR.ACVest <- acf(simAR$X2[1:n], lag.max=max(theLags), type="cov",
                        plot=FALSE)
    smplVals[[2]][which(smplVals[[2]]$N==n),paste0("lag",theLags)] <-
      simAR.ACVest$acf
    
    simAR.ACVest <- ccf(simAR$X1[1:n], simAR$X2[1:n], lag.max=max(theLags),
                        type="cov", plot=FALSE)
    smplVals[[3]][which(smplVals[[3]]$N==n),paste0("lag",theLags)] <-
      simAR.ACVest$acf[which(simAR.ACVest$lag %in% theLags)]
    smplVals[[4]][which(smplVals[[4]]$N==n),paste0("lag",rev(theLags))] <-
      simAR.ACVest$acf[which(simAR.ACVest$lag %in% ((-1)*theLags))]
  }
  remove(simAR,simAR.ACVest,n)
  for (k in 1:4) { # update dataframes
    smplVals[[k]] <- data.frame(N=numObsVec,
                                smplVals[[k]][,-1]+prev.df[[k]][,-1])
  }
}
cat(paste0("########## end time","\t",Sys.time(),"\n"))

for (k in 1:4) { # divide to get avg.
  smplVals[[k]] <- data.frame(N=numObsVec, smplVals[[k]][,-1]/NUM_RLZNS)
}

# plot results for sample ACVF for X1 series
png(filename="sample-acvf-x1.png", width=640, height=480)
par(mfrow=plotSqGridDims(length(theLags)))
for (l in theLags) {
  par(mar=c(4,4,3,1), mgp=c(2.5, 1, 0))
  plot(x=smplVals[[1]]$N, y=smplVals[[1]][,paste0("lag",l)], type="p",
       xlab="TS realization size",
       ylab=paste0("Lag-", l, " sample ACVF"),
       main=paste0("ACVF of X1 (Avg for ",NUM_RLZNS," realization[s])"),
       ylim=vAxisLims(smplVals[[1]][,paste0("lag",l)], CCVmats[1,2*(maxlag+l)+1])
  )
  abline(h=CCVmats[1,2*(maxlag+l)+1], col="deepskyblue")
}
dev.off()

# plot results for sample ACVF for X2 series
png(filename="sample-acvf-x2.png", width=640, height=480)
par(mfrow=plotSqGridDims(length(theLags)))
for (l in theLags) {
  par(mar=c(4,4,3,1), mgp=c(2.5, 1, 0))
  plot(x=smplVals[[2]]$N, y=smplVals[[2]][,paste0("lag",l)], type="p",
       xlab="TS realization size",
       ylab=paste0("Lag-", l, " sample ACVF"),
       main=paste0("ACVF of X2 (Avg for ",NUM_RLZNS," realization[s])"),
       ylim=vAxisLims(smplVals[[2]][,paste0("lag",l)], CCVmats[2,2*(maxlag+l)+2])
  )
  abline(h=CCVmats[2,2*(maxlag+l)+2], col="deepskyblue")
}
dev.off()


# plot results for sample CCVF of X1 & X2
png(filename="sample-ccvf-x1x2.png", width=640, height=480)
par(mfrow=plotSqGridDims(length(theLags)))
for (l in theLags) {
  par(mar=c(4,4,3,1), mgp=c(2.5, 1, 0))
  plot(x=smplVals[[3]]$N, y=smplVals[[3]][,paste0("lag",l)], type="p",
       xlab="TS realization size",
       ylab=paste0("Lag-", l, " sample CCVF"),
       main=paste0("CCVF of X1 & X2 (Avg for ",NUM_RLZNS," realization[s])"),
       ylim=vAxisLims(smplVals[[3]][,paste0("lag",l)], CCVmats[1,2*(maxlag+l)+2])
  )
  abline(h=CCVmats[1,2*(maxlag+l)+2], col="deepskyblue")
}
dev.off()


# plot results for sample CCVF of X2 & X1
png(filename="sample-ccvf-x2x1.png", width=640, height=480)
par(mfrow=plotSqGridDims(length(theLags)))
for (l in theLags) {
  par(mar=c(4,4,3,1), mgp=c(2.5, 1, 0))
  plot(x=smplVals[[4]]$N, y=smplVals[[4]][,paste0("lag",l)], type="p",
       xlab="TS realization size",
       ylab=paste0("Lag-", l, " sample CCVF"),
       main=paste0("CCVF of X2 & X1 (Avg for ",NUM_RLZNS," realization[s])"),
       ylim=vAxisLims(smplVals[[4]][,paste0("lag",l)], CCVmats[2,2*(maxlag+l)+1])
  )
  abline(h=CCVmats[2,2*(maxlag+l)+1], col="deepskyblue")
}
dev.off()
# end ACCVFs sample vs. theoretical







#################################### Spectra ###################################
#
# Spectra... using bivAR1
#
M <- 3*2^(trunc(log2(numObs))+1) # pad to this length
inds <- seq(0,floor(M/2)-1,1)

#
# bivAR1.th.spec : Computes the theoretical spectra values of the bivAR1 process
#   at specified frequency `l` between -0.5 and -0.5 and auto-cross-covariance
#   matrices `accvf.mats` (specified as a row of m-by-m matrices for lags -h to
#   +h). The returned m-by-m matrix is the result of a matrix multiplication
#   involving complex exponentials.
#
bivAR1.th.spec <- function(l, accvf.mats) {
  stopifnot(abs(l)<=0.5, dim(accvf.mats)[1]==2, dim(accvf.mats)[2] %% 2 == 0)
  
  h <- (dim(accvf.mats)[2]/2 - 1) / 2
  hseq <- seq(-h,h,1)
  expDiag <- diag(as.vector(matrix(exp((0-1i)*2*pi*l*rep(hseq,2)), 2, h*2+1,
                                   byrow=TRUE)
                            )
                  )
  idents <- matrix(rep(diag(1,2),2*h+1), 2*(2*h+1),2, byrow=TRUE)
  specMat <- accvf.mats %*% expDiag %*% idents
  return(specMat)
}
# end bivAR1.th.spec


thSpec <- list(freq=seq(-0.5,0.5,length.out=numf<-101))
thSpec$spec <- matrix(nrow=2, ncol=2*numf)
for (k in 1:numf) {
  startcol <- 2*(k-1) + 1 # left-most column for 2x2 matrix minor
  thSpec$spec[,c(0,1)+startcol] <- bivAR1.th.spec(l=thSpec$freq[k], CCVmats)
}


autospec1 <- spec.mtm(ts(bivAR1$X1), nFFT=M, plot=FALSE)
plot(autospec1)
lines(x=thSpec$freq, y=Re(thSpec$spec[1,(2*seq(1,numf)-1)]),
      lwd=2, col="blue")
dev.off()

autospec2 <- spec.mtm(ts(bivAR1$X2), nFFT=M, plot=FALSE)
plot(autospec2)
lines(x=thSpec$freq, y=Re(thSpec$spec[2,(2*seq(1,numf))]),
      lwd=2, col="blue")
dev.off()


# compute the Slepians
sleps <- multitaper::dpss(numObs, k=7, nw=4)

K <- dim(sleps$v)[2] # number of tapers
X1mtx <- matrix(bivAR1$X1, nrow=length(bivAR1$X1), ncol=K)
X1mtx <- mvfft(z=X1mtx*sleps$v)
X2mtx <- matrix(bivAR1$X2, nrow=length(bivAR1$X2), ncol=K)
X2mtx <- mvfft(z=X2mtx*sleps$v)

crossSpecEst12 <- (X1mtx * Conj(X2mtx)) %*% rep(1,K)/K

# plot (1,2) cross spectrum, real and imaginary parts
par(mar=c(4,4,3,1), mgp=c(2.5, 1, 0))
par(mfcol=c(2,1))
for (prt in c("Re","Im")) {
  yvals <- parse(text=paste0("y=",prt,"(crossSpecEst12)[c(seq(numObs/2+1,numObs), seq(1,numObs/2))]"))
  logy <- ifelse(prt=="Re", "y", "")
  
  plot(x=(seq(0,numObs-1)-numObs/2)/numObs,
       y=eval(yvals), type="l", xlab="Frequency", log=logy,
       ylab=paste0("Cross Spectral Estimate (", prt, ")"))
  
  yvals.l <- parse(text=paste0(prt,"(thSpec$spec[1,(2*seq(1,numf))])"))
  lines(x=thSpec$freq, y=eval(yvals.l), type="l", lwd=2,col="blue")
}
dev.off()

# Note: The following estimate has lags in the order
#     (0,1,2,...,numObs/2 - 1, -numObs/2,...,-2,-1,)
#   so that element 1 is the 0-th lag.
ccv12est <- fft(z=crossSpecEst12, inverse=TRUE) / numObs
plot(x=(-12:12), y=Re(ccv12est)[c(seq(numObs-11,numObs),seq(1,13))], type="h",
     xlab="Lag", ylab="CCVF[12] estimate")
abline(h=0)



# end Spectra

