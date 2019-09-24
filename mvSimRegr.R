#!/usr/bin/env Rscript

tryCatch(setwd("~/multivar-sims"),
         error=function(e) {
           tryCatch( setwd("~/WORK/Q5/masters-thesis/code/multivar-sims") )
         }
)

library(mAr) # depends on MASS
library(multitaper)

source("helpers_mvSimRegr.R") # import helper functions

################################## set params ##################################
# parameters for bivariate AR(1)
phiMat.r <- matrix(c(0.7,-0.5,0.6,0.2), nrow=2, ncol=2, byrow=TRUE)
phiMat.p <- matrix(c(0.7,-0.2,0.5,-0.3), nrow=2, ncol=2, byrow=TRUE)
errCovMat.r <- 0.5 * (diag(1,2) + matrix(c(0,0.00,0.00,0), 2, 2))
errCovMat.p <- 0.5 * diag(1,2)

X <- ar.regr.cov(phiMat.p=phiMat.p, phiMat.r=phiMat.r,
                 errCovMat.r=errCovMat.r, errCovMat.p=errCovMat.p,
                 numObsVec=seq(5,5,5) * 1e2, NUM_REGR=100,
                 mtmFixed="NW", W=0.005, timeBandProd=6, numTapers=11,
                 adaptWt=TRUE, writeImgFile=FALSE, embedSines=TRUE,
                 linDepY=FALSE, computeCorr=FALSE)

mean(X$betasOverN[[1]]$b1)
mean(X$betasOverN[[1]]$b2)
head(X$betasOverN[[1]])
X$result



# parameters for bivariate AR(p)
phiMat.r <- matrix(c(0.99,0.3,-0.95,0.2,0.1,0.85,0.1,-0.2), nrow=2, ncol=4, byrow=TRUE)
phiMat.p <- matrix(c(0.8,0,-0.7,0,0,0.5,0,-0.3), nrow=2, ncol=2, byrow=TRUE)
errCovMat.r <- 0.5 * (diag(1,2) + matrix(c(0,0.00,0.00,0), 2, 2))
errCovMat.p <- 0.5 * diag(1,2)

X <- ar.regr.cov(phiMat.p=phiMat.p, phiMat.r=phiMat.r,
                 errCovMat.r=errCovMat.r, errCovMat.p=errCovMat.p,
                 numObsVec=seq(4,4,1) * 1e2, NUM_REGR=100,
                 mtmFixed="NW", W=0.005, timeBandProd=6, numTapers=11,
                 adaptWt=FALSE, writeImgFile=FALSE, embedSines=TRUE,
                 linDepY=FALSE, computeCorr=FALSE)







# ACVF plots
par(mar=c(4,4,1,1))
plot(x=(-50:50), y=X$theo.ccv.r[(-50:50)+4e2], type="h",
     ylim=range( c(X$theo.ccv.r, X$mtap.ccv.r, X$bart.ccv.r$acf) ),
     ylab="CCVF (Y1,Y2)", xlab="Lag")
abline(h=0)
points(x=(-50:50), y=X$mtap.ccv.r[(-50:50)+4e2], col="blue")
points(x=(-50:50), y=X$bart.ccv.r$acf[(-50:50)+4e2], col="goldenrod")






