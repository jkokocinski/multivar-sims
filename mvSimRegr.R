#!/usr/bin/env Rscript

tryCatch(setwd("~/multivar-sims"),
         error=function(e) {
           tryCatch( setwd("~/WORK/Q5/masters-thesis/code/multivar-sims"))
         }
)

source("helpers_mvSimRegr.R") # import helper functions

################################## set params ##################################
# parameters for bivariate AR(1)
phiMat.r <- matrix(c(0.7,-0.5,0.6,0.2), nrow=2, ncol=2, byrow=TRUE)
phiMat.p <- matrix(c(0.4,0,0,0.3), nrow=2, ncol=2, byrow=TRUE)

X <- ar1.regr.cov(phiMat.p=phiMat.p, phiMat.r=phiMat.r, varZ=1,
                  numObsVec=seq(1,5,1) * 1e2,
                  NUM_REGR=100,
                  mtmFixed="NW", timeBandProd=6, numTapers=11,
                  writeImgFile=FALSE, embedSines=FALSE)




# # variance plot
# if (writeImgFile) {
#   png("img/varBetahats.png", width=640, height=360)
# }
# par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
# plot(x=result$N, y=result$var.bart, xlab="Realization Size", ylab="var",
#      log="y", ylim=range(result[,c("var.bart","var.mtap")])*c(1e-1,1e1))
# points(x=result$N, y=result$var.mtap, col="blue", pch=2)
# legend("topright", legend=c("Bartlett","Multitaper"), pch=c(1,2),
#        col=c("black","blue"))
# text(x=min(numObsVec), y=min(result[,c("var.bart","var.mtap")])*0.1, adj=c(0,0),
#      labels=plotComment)
# if (writeImgFile) { dev.off() }



# ggplot(result, aes(x=N, y=mse.mtap)) +
#   geom_point(size = 2) +
#   geom_errorbar(aes(ymin=qSE.02.bart, ymax=qSE.98.bart))





# # do the regressions in freq.domain .......................................
# M <- 3*2^(trunc(log2(numObs))+1) # pad to this length
# inds <- seq(0,floor(M/2)-1,1)
# 
# # compute the Slepians
# sleps <- multitaper::dpss(M, k=11, nw=6)
# 
# K <- dim(sleps$v)[2] # number of tapers
# X1mtx <- matrix(c(bivAR1.p[,1], rep(0,M-numObs)), nrow=M, ncol=K)
# X1mtx <- mvfft(z=X1mtx*sleps$v)
# Y1mtx <- matrix(c(bivAR1.r[,1], rep(0,M-numObs)), nrow=M, ncol=K)
# Y1mtx <- mvfft(z=Y1mtx*sleps$v)
# 
# autoSpecEstXX <- Re(X1mtx * Conj(X1mtx)) %*% rep(1,K)/K
# autoSpecEstYY <- Re(Y1mtx * Conj(Y1mtx)) %*% rep(1,K)/K
# crossSpecEstXY <- (X1mtx * Conj(Y1mtx)) %*% rep(1,K)/K
# 
# mtap.ccv.r <- fft(z=crossSpecEstXY, inverse=TRUE) / numObs
# # put entries in "correct" order, from -numObs to +numObs
# mtap.ccv.r <- Re(mtap.ccv.r[c(seq(M-(numObs-1)+1,M), seq(1,numObs)), 1])
# 
# xferFunEst <- Conj(crossSpecEstXY) / (Re(autoSpecEstXX))
# beta.hat <- Re(fft(z=xferFunEst, inverse=TRUE)) / M
# head(beta.hat)
# beta.hat <- beta.hat[1:(M/2)] # truncate to length M/2
# beta.hat <- beta.hat[1:(bt<-3)] # truncate to length bt
# 
# 
# # compute resids
# fv <- vector(mode="numeric", length=numObs)
# for (j in bt:numObs) {
#   fv[j] <- sum( rev(beta.hat[1:bt]) * bivAR1.p.s[(1:bt)+j-bt,1] )
# }
# rs <- bivAR1.r.s[,1] - fv





