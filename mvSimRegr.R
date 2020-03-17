#!/usr/bin/env Rscript

tryCatch(setwd("~/multivar-sims"),
         error=function(e) {
           tryCatch( setwd("~/WORK/Q5/masters-thesis/code/multivar-sims") )
         }
)

library(mAr) # depends on MASS
library(multitaper)
source("helpers_mvSimRegr.R") # import helper functions
source("seasonalFunctions.R")

######################## parameters for bivariate AR(1) ########################
# phiMat.r <- matrix(c(0.7,-0.5,0.6,0.2), nrow=2, ncol=2, byrow=TRUE)
# phiMat.p <- matrix(c(0.7,-0.2,0.5,-0.3), nrow=2, ncol=2, byrow=TRUE)
# errCovMat.r <- 1 * (diag(1,2) + matrix(c(0,0.55,0.55,0), 2, 2))
# errCovMat.p <- 1 * diag(1,2)

######################## parameters for bivariate AR(2) ########################
phiMat.r <- matrix(c(0.99,0.1,-0.95,0.3,0.1,0.88,0.1,-0.5), nrow=2, ncol=4, byrow=TRUE)
phiMat.p <- matrix(c(0.8,0,-0.7,0,0,0.5,0,-0.3), nrow=2, ncol=4, byrow=TRUE)
errCovMat.r <- 1e1 * (diag(c(1,1), 2,2) + matrix(c(0,0.00,0.00,0), 2, 2))
errCovMat.p <- 1e1 * diag(1,2)

###### parameters for bivariate AR(4) -- based on fit from demand & HOEP #######
# phiMat.p <- rbind(
#   c( 1.5103, 0.0926, -0.9208, -0.1638, 0.5057, 0.1363, -0.1155, -0.0689 ),
#   c( 0.0770, 1.5231, -0.1707, -0.8566, 0.1523, 0.4166, -0.0613, -0.1054 )
# )
# # AR(2) for HOEP
# phiMat.r <- rbind(
#   c(0.6275, 0.0122,  0.1215, -0.0388),
#   c(0.0084, 0.6282, -0.0370,  0.1294)
# )
# errCovMat.p <- matrix(c(77572,40692,40692,75347), 2, 2)
# errCovMat.r <- matrix(c(338,15,15,350), 2, 2)



################## Run the process simulations and regressions #################
X <- ar.regr.cov(phiMat.p=phiMat.p, phiMat.r=phiMat.r,
                 errCovMat.r=errCovMat.r, errCovMat.p=errCovMat.p,
                 numObsVec=c(1024), NUM_REGR=100,
                 mtmFixed="NW", W=0.01, timeBandProd=9, numTapers=17,
                 adaptWt=FALSE, embedSines=TRUE,
                 linDepY=FALSE, computeCorr=TRUE, removeLCs=TRUE)

plotCIs(X, stage="", type="cov", writeImgFile=F)
plotCIs(X, stage="s", type="cov", writeImgFile=F)
plotCIs(X, stage="s.w", type="cov", writeImgFile=F)

boxplot(cbind(X$betasOverN[[1]]$se.cv.bart,X$betasOverN[[1]]$se.cv.mtap), log="y")
boxplot(cbind(X$betasOverN.s[[1]]$se.cv.bart,X$betasOverN.s[[1]]$se.cv.mtap), log="y")
boxplot(cbind(X$betasOverN.s.w[[1]]$se.cv.bart,X$betasOverN.s.w[[1]]$se.cv.mtap), log="y")

plotCIcompare(X, type="cor", estType="bart")
plotCIcompare(X, type="cor", estType="mtap")

################################## CCVF plots ##################################
# png("img/ccvf-plot.png", width=640, height=360)
# pdf("img/ccvf-plot.pdf", width=6.4, height=3.6)
plotCCVF(resultList=X, plotLags=seq(-200,200), stage="", ave=F)
plotCCVF(resultList=X, plotLags=seq(-200,200), stage="s", ave=F)
plotCCVF(resultList=X, plotLags=seq(-200,200), stage="s.w", ave=F)
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
png("img/cov-plot.png", width=640, height=360)
par(mar=c(4,4,1,1))
with(X$betasOverN.s.w[[1]], {
  plot(cv.bart, ylim=range(c(cv.bart, cv.mtap, cv.theo, X$bhCovCIs$smpl.cov[1])),
       col="goldenrod", xlab="Regression #", ylab="cov(b1,b2)")
points(x=1:length(cv.mtap), y=cv.mtap,
       col="blue")
abline(h=cv.theo[1], col="forestgreen", lwd=2)
})
abline(h=X$bhCovCIs.s.w$smpl.cov[1], col="red", lty=2)
dev.off()
# end cov plot

############################### correlation plot ###############################
png("img/cor-plot.png", width=640, height=360)
par(mar=c(4,4,1,1))
with(X$betasOverN[[1]], {
  plot(cor.bart, ylim=range(c(cor.bart, cor.mtap, cor.theo, X$bhCovCIs$smpl.cor[1])),
       col="goldenrod", xlab="Regression #", ylab="cor(b1,b2)")
  points(x=1:length(cor.mtap), y=cor.mtap, col="blue")
  abline(h=cor.theo[1], col="forestgreen", lwd=2)
})
abline(h=X$bhCovCIs$smpl.cor[1], col="red", lty=2)
dev.off()
# end cor plot




# Plot sample amplitude cross spectrum for Bartlett and MTM for final rlzn and
#   and compare to the theoretical.
bart.cspec <- with(X, {
  fft(z=c(respRlzns[,1,params.init$NUM_REGR],rep(0,1024))) *
    Conj(fft(z=c(X$respRlzns[,2,params.init$NUM_REGR],rep(0,1024)))) / 2048
  }
)
plot.lims <- range(c(Mod(bart.cspec[-1])**2, Mod(X$crossSpec)**2, Mod(specs[1,2,])**2))
specs <- bivARp.spec(phiMat=phiMat.r, V=errCovMat.r, freqs=c(0,seq(1,1024))/2048)
par(mar=c(4,4,1,1))
plot(x=c(0,seq(1,1024)/2048), y=Mod(bart.cspec[1:1025])**2,
     type="l", log="y", col="red", ylim=plot.lims,
     xlab="Frequency", ylab="Amplitude Cross Spectrum (Y1,Y2)")
lines(x=c(0,seq(1,1024)/2048), y=Mod(X$crossSpec[1:1025])**2, col="blue")
lines(x=c(0,seq(1,1024)/2048), y=Mod(specs[1,2,])**2, col="forestgreen", lwd=2)
dev.off()







