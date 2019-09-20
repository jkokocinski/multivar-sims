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
phiMat.p <- phiMat.r # matrix(c(0.7,-0.2,0.5,-0.3), nrow=2, ncol=2, byrow=TRUE)
errCovMat.r <- 0.5 * (diag(1,2) + matrix(c(0,0.00,0.00,0), 2, 2))
errCovMat.p <- 0.5 * diag(1,2)

X <- ar1.regr.cov(phiMat.p=phiMat.p, phiMat.r=phiMat.r,
                  errCovMat.r=errCovMat.r, errCovMat.p=errCovMat.p,
                  numObsVec=seq(5,5,5) * 1e3, NUM_REGR=100,
                  mtmFixed="NW", W=0.005, timeBandProd=9, numTapers=17,
                  adaptWt=TRUE, writeImgFile=FALSE, embedSines=TRUE,
                  linDepY=FALSE, computeCorr=FALSE)

mean(X$betasOverN[[1]]$b1)
mean(X$betasOverN[[1]]$b2)
head(X$betasOverN[[1]])
X$result



