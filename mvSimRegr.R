#!/usr/bin/env Rscript

tryCatch(setwd("~/multivar-sims"),
         error=function(e) {
           tryCatch( setwd("~/WORK/Q5/masters-thesis/code/multivar-sims") )
         }
)

library(mAr) # depends on MASS
library(multitaper)

source("helpers_mvSimRegr.R") # import helper functions

######################## parameters for bivariate AR(1) ########################
phiMat.r <- matrix(c(0.7,-0.5,0.6,0.2), nrow=2, ncol=2, byrow=TRUE)
phiMat.p <- matrix(c(0.7,-0.2,0.5,-0.3), nrow=2, ncol=2, byrow=TRUE)
errCovMat.r <- 1 * (diag(1,2) + matrix(c(0,0.00,0.00,0), 2, 2))
errCovMat.p <- 1 * diag(1,2)

X <- ar.regr.cov(phiMat.p=phiMat.p, phiMat.r=phiMat.r,
                 errCovMat.r=errCovMat.r, errCovMat.p=errCovMat.p,
                 numObsVec=seq(5,5,5) * 1e3, NUM_REGR=1,
                 mtmFixed="NW", W=0.005, timeBandProd=10, numTapers=19,
                 adaptWt=TRUE, writeImgFile=FALSE, embedSines=FALSE,
                 linDepY=FALSE, computeCorr=TRUE)

mean(X$betasOverN[[1]]$b1)
mean(X$betasOverN[[1]]$b2)
head(X$betasOverN[[1]])
X$result
# end AR(1)


######################## parameters for bivariate AR(2) ########################
phiMat.r <- matrix(c(0.99,0.1,-0.95,0.3,0.1,0.85,0.1,-0.2), nrow=2, ncol=4, byrow=TRUE)
phiMat.p <- matrix(c(0.8,0,-0.7,0,0,0.5,0,-0.3), nrow=2, ncol=4, byrow=TRUE)
errCovMat.r <- 1e1 * (diag(c(1,1), 2,2) + matrix(c(0,0.00,0.00,0), 2, 2))
errCovMat.p <- 1e1 * diag(1,2)

X <- ar.regr.cov(phiMat.p=phiMat.p, phiMat.r=phiMat.r,
                 errCovMat.r=errCovMat.r, errCovMat.p=errCovMat.p,
                 numObsVec=seq(10,10,2) * 1e2, NUM_REGR=100,
                 mtmFixed="NW", W=0.005, timeBandProd=6, numTapers=11,
                 adaptWt=TRUE, writeImgFile=TRUE, embedSines=FALSE,
                 linDepY=FALSE, computeCorr=TRUE)
# end AR(2)


################################## CCVF plots ##################################
png("img/ccvf-plot.png", width=640, height=360)
par(mar=c(4,4,1,1))
N <- max(X$result$N)
plot(x=(-200:200), y=X$theo.ccv.r[(-200:200)+N], type="h",
     ylim=range( c(X$theo.ccv.r, X$mtap.ccv.r, X$bart.ccv.r$acf) ),
     ylab="CCVF (Y1,Y2)", xlab="Lag")
abline(h=0)
points(x=(-200:200), y=X$mtap.ccv.r[(-200:200)+N], col="blue")
points(x=(-200:200), y=X$bart.ccv.r$acf[(-200:200)+N], col="goldenrod")
dev.off()
# end ACVF plots

################################ covariance plot ###############################
png("img/cov-plot.png", width=640, height=360)
par(mar=c(4,4,1,1))
with(X$betasOverN[[1]], {
  plot(cv.bart, ylim=range(c(cv.bart, cv.mtap, cv.theo, X$result$smpl.cov)),
       col="goldenrod", xlab="Regression #", ylab="cov(b1,b2)")
})
points(x=1:100, y=X$betasOverN[[1]]$cv.mtap, col="blue")
abline(h=X$betasOverN[[1]]$cv.theo[1], col="forestgreen", lwd=2)
abline(h=X$result$smpl.cov[1], col="red", lty=2)
dev.off()
# end cov plot

############################### correlation plot ###############################
png("img/cor-plot.png", width=640, height=360)
par(mar=c(4,4,1,1))
with(X$betasOverN[[1]], {
  plot(cor.bart, ylim=range(c(cor.bart, cor.mtap, cor.theo, X$result$smpl.cor)),
         col="goldenrod", xlab="Regression #", ylab="cor(b1,b2)")
})
points(x=1:100, y=X$betasOverN[[1]]$cor.mtap, col="blue")
abline(h=X$betasOverN[[1]]$cor.theo[1], col="forestgreen", lwd=2)
abline(h=X$result$smpl.cor[1], col="red", lty=2)
dev.off()
# end cor plot




# parameters for bivariate AR(6) -- based on fit from demand & HOEP
phiMat.r <- rbind(
  c(2.01,-0.48,-1.49,0.12,0.64,0.6,-0.2,0.25,-0.08,0.11,0.06,0.24),
  c(0.02,0.44,-0.02,0.12,0.02,0.01,-0.01,0.07,0,0.02,0,0.01)
  )
phiMat.p <- cbind(diag(c(0.8,0.5)), diag(c(-0.7,-0.3)), diag(c(0.5,0.5)))
errCovMat.r <- matrix(c(143195,1631,1631,932), 2, 2)
errCovMat.p <- 10000 * diag(1,2)







