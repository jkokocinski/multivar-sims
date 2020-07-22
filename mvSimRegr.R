#!/usr/bin/env Rscript

tryCatch(setwd("~/multivar-sims"),
         error=function(e) {
           tryCatch( setwd("~/WORK/Q5/multivar-sims") )
         }
)

library(mAr) # depends on MASS
library(multitaper)
source("helpers_mvSimRegr.R") # import helper functions
source("seasonalFunctions.R")

set.seed(148)

######################## parameters for bivariate AR(1) ########################
phiMat.r <- matrix(c(0.7,-0.5,0.6,0.2), nrow=2, ncol=2, byrow=TRUE)
phiMat.p <- matrix(c(0.7,-0.2,0.5,-0.3), nrow=2, ncol=2, byrow=TRUE)
errCovMat.r <- 1 * (diag(1,2) + 0*matrix(c(0,0.25,0.25,0), 2, 2))
errCovMat.p <- 1 * diag(1,2)

######################## parameters for bivariate AR(2) ########################
# phiMat.r <- matrix(c(0.99,0.1,-0.95,0.3,0.1,0.88,0.1,-0.5), nrow=2, ncol=4, byrow=TRUE)
# phiMat.p <- matrix(c(0.8,0.15,-0.7,0.2,0.3,0.5,0.35,-0.3), nrow=2, ncol=4, byrow=TRUE)
# errCovMat.r <- 1e0 * (diag(c(1,1), 2,2) + matrix(c(0,0.00,0.00,0), 2, 2))
# errCovMat.p <- 1e0 * diag(1,2)

###### parameters for bivariate AR(4) -- based on fit from demand & HOEP #######
# phiMat.p <- rbind(
#   c( 1.5103, 0.0926, -0.9208, -0.1638, 0.5057, 0.1363, -0.1155, -0.0689 ),
#   c( 0.0770, 1.5231, -0.1707, -0.8566, 0.1523, 0.4166, -0.0613, -0.1054 )
# )
# # AR(2) for HOEP
# phiMat.r <- rbind(
#   c(0.6504, 0.024, -0.0708, -0.0085, 0.1016, 0.003, 0.0858, -0.0154),
#   c(0.0095, 0.5936, -0.0067, -0.0267, 0.0048, 0.0826, -0.0058, 0.0847)
# )
# errCovMat.p <- matrix(c(77572,40692,40692,75347), 2, 2)
# errCovMat.r <- matrix(c(364,12,12,371), 2, 2)



################## Run the process simulations and regressions #################
X <- ar.regr.cov(phiMat.p=phiMat.p, phiMat.r=phiMat.r,
                 errCovMat.r=errCovMat.r, errCovMat.p=errCovMat.p,
                 numObsVec=c(4000), NUM_REGR=100,
                 mtmFixed="NW", W=0.01, timeBandProd=4, numTapers=7,
                 adaptWt=TRUE, embedSines=TRUE,
                 linDepY=TRUE, removeLCs=TRUE)

plotCIs(X, stage="",    type="cor", writeImgFile=F)
plotCIs(X, stage="s",   type="cov", writeImgFile=F)
plotCIs(X, stage="s.w", type="cov", writeImgFile=F)

boxplot(cbind(X$betasOverN[[1]]$se.cv.bart,X$betasOverN[[1]]$se.cv.mtap), log="y")
boxplot(cbind(X$betasOverN.s[[1]]$se.cv.bart,X$betasOverN.s[[1]]$se.cv.mtap), log="y")
boxplot(cbind(X$betasOverN.s.w[[1]]$se.cv.bart,X$betasOverN.s.w[[1]]$se.cv.mtap), log="y")

plotCIcompare(X, type="cor", estType="bart")
plotCIcompare(X, type="cor", estType="mtap")
plotCIcompare(X, type="cor", estType="both")

################################## CCVF plots ##################################
# png("img/ccvf-plot.png", width=640, height=320)
# pdf("img/ccvf-plot.pdf", width=6.4, height=3.6)
plotCCVF(resultList=X, plotLags=seq(-40,40), stage="", ave=F)
plotCCVF(resultList=X, plotLags=seq(-40,40), stage="s", ave=F)
plotCCVF(resultList=X, plotLags=seq(-40,40), stage="s.w", ave=F)
# plotLags <- seq(-20,20)
# par(mar=c(4,4,1,1))
# with(X,
#   {
#   N <- max(bhCovCIs$N)
#   plot(x=plotLags, y=theo.ccv.r[plotLags+N], type="h", lwd=2,
#        ylim=range( c(theo.ccv.r[plotLags+N], mtap.ccv.r.ave[plotLags+N],
#                      bart.ccv.r.ave[plotLags+N]) ),
#        ylab="CCVF (Y1,Y2)", xlab="Lag")
#   abline(h=0)
#   points(x=plotLags, y=mtap.ccv.r.ave[plotLags+N], col="blue")
#   points(x=plotLags, y=bart.ccv.r.ave[plotLags+N], col="red")
#   }
# )
dev.off()
# end ACVF plots

################################ covariance plot ###############################
# png("img/cov-plot.png", width=640, height=320)
plotCovB(resultList=X, type="cov", stage="")
plotCovB(resultList=X, type="cov", stage="s")
plotCovB(resultList=X, type="cov", stage="s.w")
# dev.off()
# end cov plot

############################### correlation plot ###############################
# png("img/cor-plot.png", width=640, height=320)
plotCovB(resultList=X, type="cor", stage="")
plotCovB(resultList=X, type="cor", stage="s")
plotCovB(resultList=X, type="cor", stage="s.w")
# dev.off()
# end cor plot




# Plot sample amplitude cross spectrum for Bartlett and MTM for final rlzn and
#   and compare to the theoretical.
N <- X$params.init$numObsVec[1]
nFFT <- length(X$crossSpec)
bart.cspec <- with(X, {
  fft(z=c(respRlzns[,1,params.init$NUM_REGR],rep(0,nFFT-numObs))) *
    Conj(fft(z=c(X$respRlzns[,2,params.init$NUM_REGR],rep(0,nFFT-numObs)))) / nFFT
  }
)
specs <- bivARp.spec(phiMat=phiMat.r, V=errCovMat.r, freqs=c(0,seq(1,nFFT/2))/nFFT)
plot.lims <- range(c(Mod(bart.cspec[-1])**2, Mod(X$crossSpec)**2, Mod(specs[1,2,])**2))
par(mar=c(4,4,1,1))
plot(x=c(0,seq(1,nFFT/2)/nFFT), y=Mod(bart.cspec[1:(nFFT/2+1)])**2,
     type="l", log="y", col="red", ylim=plot.lims,
     xlab="Frequency", ylab="Amplitude Cross Spectrum (Y1,Y2)")
lines(x=c(0,seq(1,nFFT/2)/nFFT), y=Mod(X$crossSpec[1:(nFFT/2+1)])**2, col="blue")
lines(x=c(0,seq(1,nFFT/2)/nFFT), y=Mod(specs[1,2,])**2, col="forestgreen", lwd=2)







