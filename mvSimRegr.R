#!/usr/bin/env Rscript

tryCatch(setwd("~/multivar-sims"),
				 error=function(e) {
								 tryCatch( setwd("~/WORK/Q5/masters-thesis/code/multivar-sims"))
				 }
)


############################ packages and functions ############################
library(mAr)
library(multitaper)
library(expm)

source("psumCov.R") # import psumCov function for computing A/CCVF
source("toepLeftMult.R") # import toepLeftMult function for computing a Toeplitz
                         #   matrix left-multiplying a vector

################################## set params ##################################
writeImgFile <- FALSE # should plots be written to files? or just displayed?
mtmFixed <- "NW" # should 'W' or 'NW' be fixed for mtm? Specify values below.
embedSines <- FALSE # should sinusoids be embedded in process?

if (mtmFixed=="NW") {
  # mtm parameters for fixed NW
  timeBandProd <- 6
  numTapers <- 11
  if ("W" %in% ls()) { rm(W) }
} else if (mtmFixed=="W") {
  # if just W is fixed
  W <- 0.09
  if ("timeBandProd" %in% ls() | "numTapers" %in% ls()) {
    rm(list=intersect(c("timeBandProd","numTapers"), ls()))
  }
} else {
  warning("`mtmFixed` not set properly. Setting to 'NW' and using NW=4, K=7.")
  timeBandProd <- 4
  numTapers <- 7
}

# vector of realization sizes to be used in the simulations
numObsVec <- seq(1,1,1) * 1e2

# number of regressions to be perfomed for each realization size
NUM_REGR <- 100

# parameters for bivariate AR(1)
varZ <- 0.01 # variance of innovations
phiMat.r <- matrix(c(0.4,-0.6,-0.6,0.4), nrow=2, ncol=2, byrow=TRUE)
errCovMat.r <- diag(varZ,2,2)
phiMat.p <- matrix(c(0.4,0,0,0.3), nrow=2, ncol=2, byrow=TRUE)
errCovMat.p <- diag(varZ,2,2)

betasOverN <- list() # list of data.frames for the \hat{beta}s
betas.head <- c("b1","b2","cv.bart","cv.theo","cv.mtap") # headings for betas


############################# do main computations #############################
for (n in seq(1,length(numObsVec))) {
  numObs <- numObsVec[n]

  # compute theoretical A/CCVF matrices for (Y_1, Y_2), the response series
  matrixSeries <- psumCov(phiMat.r, V=errCovMat.r)
  maxlag <- numObs-1
  CCVmats <- matrix(rep(matrixSeries,2*maxlag+1), nrow=2, ncol=2*(2*maxlag+1))
  for (h in seq(1,maxlag,1)) {
    startCol <- 2*maxlag+1 + 2*h # left-most column in matrix for lag +h
    CCVmats[,seq(startCol,startCol+1)] <- phiMat.r %*%
      CCVmats[,seq(startCol-2,startCol-1)]
    
    startCol <- dim(CCVmats)[2]-startCol # left-most column in matrix for lag -h
    CCVmats[,seq(startCol,startCol+1)] <- CCVmats[,seq(startCol+2,startCol+3)] %*%
      t(phiMat.r)
  }
  theo.ccv.r <- CCVmats[1,2*seq(0,2*maxlag)+2] # \cov{Y_1,Y_2} from -maxlag to +maxlag
  # theo.ccv.r[which(theo.ccv.r<1e-18)] <- 0 # truncate at 1e-18
  # theo.ccv.r.mat <- matrix(nrow=numObs, ncol=numObs)
  # theo.ccv.r.mat <- matrix(theo.ccv.r[numObs +
  #                          row(theo.ccv.r.mat)-col(theo.ccv.r.mat)],numObs,numObs)
  
  # values for computing spectra
  M <- 3*2^(trunc(log2(numObs))+1) # pad to this length
  # compute the Slepians
  if (mtmFixed=="NW") {
    sleps <- multitaper::dpss(M, k=numTapers, nw=timeBandProd) # fixed NW
  } else if (mtmFixed=="W") {
    sleps <- multitaper::dpss(M, k=trunc(2*numObs*W), nw=numObs*W) # fixed W
  } else {
    stop("mtmFixed must be either 'NW' or 'W'.")
  }
  
  K <- dim(sleps$v)[2] # number of tapers
  
  # sample df of \hat{beta}s for regular processes and sine-embedded
  betas <- as.data.frame(matrix(nrow=NUM_REGR, ncol=length(betas.head)))
  names(betas) <- betas.head
  # betas.s <- betas
  
  if (embedSines) {
    # generate sinusoids
    sines <- matrix(nrow=numObs, ncol=2)
    sines[,1] <- 10 * varZ * sin(2*pi/30*(1:numObs))
    sines[,2] <- 10 * varZ * sin(2*pi/20*(1:numObs))
  }
  
  bivAR1.p <- mAr.sim(w=rep(0,2), A=phiMat.p, C=errCovMat.p, N=numObs)
    
  cat(paste0("######## start : ",Sys.time()),"\n")
  for (j in 1:NUM_REGR) {
    bivAR1.r <- mAr.sim(w=rep(0,2), A=phiMat.r, C=errCovMat.r, N=numObs)
    
    if (embedSines) {
      # embed the sinusoids
      bivAR1.r <- bivAR1.r + sines
      bivAR1.p <- bivAR1.p + sines
    }
    
    # regular lin regs
    model1 <- lm(bivAR1.r[,1] ~ 0 + bivAR1.p[,1])
    model2 <- lm(bivAR1.r[,2] ~ 0 + bivAR1.p[,2])
    # model1.s <- lm(bivAR1.r.s[,1] ~ 0 + bivAR1.p.s[,1])
    # model2.s <- lm(bivAR1.r.s[,2] ~ 0 + bivAR1.p.s[,2])
    
    
    # Bartlett CCVF of response series
    bart.ccv.r <- ccf(x=bivAR1.r[,1], y=bivAR1.r[,2], type="covariance",
                      lag.max=numObs-1, plot=FALSE)
    # bart.ccv.r$acf[which(bart.ccv.r$acf<1e-18)] <- 0 # truncate at 1e-18
    # put them into a matrix (has a Toeplitz form)
    # bart.ccv.r.mat <- matrix(nrow=numObs, ncol=numObs)
    # bart.ccv.r.mat <- matrix(bart.ccv.r$acf[which(bart.ccv.r$lag==0) +
    #                          row(bart.ccv.r.mat)-col(bart.ccv.r.mat)],
    #                          numObs, numObs)
    
    # \cov(\hat{\beta}_1, \hat{\beta}_2) -- Bartlett-based
    covB.bart <- MASS::ginv(as.matrix(bivAR1.p[,1])) %*%
      toepLeftMult2( as.vector(bart.ccv.r$acf),
                    as.vector(t(MASS::ginv(as.matrix(bivAR1.p[,2]))))
    )
    # covB.bart <- MASS::ginv(as.matrix(bivAR1.p[,1])) %*% bart.ccv.r.mat %*%
    #   t(MASS::ginv(as.matrix(bivAR1.p[,2])))
    
    # \cov(\hat{\beta}_1, \hat{\beta}_2) -- theoretical-based
    covB.theo <- MASS::ginv(as.matrix(bivAR1.p[,1])) %*%
      toepLeftMult2( theo.ccv.r,
                    as.vector(t(MASS::ginv(as.matrix(bivAR1.p[,2]))))
                  )
    # covB.theo <- MASS::ginv(as.matrix(bivAR1.p[,1])) %*% theo.ccv.r.mat %*%
    #   t(MASS::ginv(as.matrix(bivAR1.p[,2])))
    # 
    # compute the cross spectrum of (Y_1,Y_2)
    Y1mtx <- matrix(c(bivAR1.r[,1], rep(0,M-numObs)), nrow=M, ncol=K)
    Y1mtx <- mvfft(z=Y1mtx*sleps$v)
    Y2mtx <- matrix(c(bivAR1.r[,2], rep(0,M-numObs)), nrow=M, ncol=K)
    Y2mtx <- mvfft(z=Y2mtx*sleps$v)
    crossSpecEstY <- (Y1mtx * Conj(Y2mtx)) %*% rep(1,K)/K
    mtap.ccv.r <- fft(z=crossSpecEstY, inverse=TRUE) / numObs
    # put entries in "correct" order, from -numObs to +numObs
    mtap.ccv.r <- Re(mtap.ccv.r[c(seq(M-(numObs-1)+1,M), seq(1,numObs)), 1])
    # mtap.ccv.r[which(mtap.ccv.r<1e-18)] <- 0 # truncate at 1e-18
    # put them into a matrix (has a Toeplitz form)
    # mtap.ccv.r.mat <- matrix(nrow=numObs, ncol=numObs)
    # mtap.ccv.r.mat <- matrix(mtap.ccv.r[numObs +
    #                          row(mtap.ccv.r.mat)-col(mtap.ccv.r.mat)],
    #                          numObs,numObs)
    
    # \cov(\hat{\beta}_1, \hat{\beta}_2) -- multitaper-based
    covB.mtap <- MASS::ginv(as.matrix(bivAR1.p[,1])) %*%
      toepLeftMult2( mtap.ccv.r,
                    as.vector(t(MASS::ginv(as.matrix(bivAR1.p[,2]))))
                  )
    # covB.mtap <- MASS::ginv(as.matrix(bivAR1.p[,1])) %*% mtap.ccv.r.mat %*%
    #   t(MASS::ginv(as.matrix(bivAR1.p[,2])))
    
    # update the betas data.frame
    betas[j,] <- c(as.numeric(model1$coefficients),
                   as.numeric(model2$coefficients),
                   covB.bart, covB.theo, covB.mtap)
    # betas.s[j,] <- c(as.numeric(model1.s$coefficients),
    #                  as.numeric(model2.s$coefficients),
    #                  covB.bart.s, covB.theo.s)
  }
  cat(paste0("########   end : ",Sys.time()),"\n")
  
  betas$se.bart <- (betas$cv.theo-betas$cv.bart)**2
  betas$se.mtap <- (betas$cv.theo-betas$cv.mtap)**2
  
  betasOverN[[n]] <- betas # update the list of dfs
} # end for loop over n in numObsVec



#################################### results ###################################
# produce result df with MSEs, etc.
result <- data.frame(N=numObsVec)
for (k in seq(1,length(betasOverN))) {
  # betasOverN[[k]][,"mse.bart"] <- ( betasOverN[[k]][,"cv.theo"] - betasOverN[[k]][,"cv.bart"] )**2
  # betasOverN[[k]][,"mse.mtap"] <- ( betasOverN[[k]][,"cv.theo"] - betasOverN[[k]][,"cv.mtap"] )**2
  dfObj <- betasOverN[[k]]
  result[k,c("mse.bart","mse.mtap")] <- apply(X=dfObj[,c("se.bart","se.mtap")],MARGIN=2,FUN=mean)
  result[k,c("var.bart","var.mtap")] <- apply(X=dfObj[,c("cv.bart","cv.mtap")],MARGIN=2,FUN=var) * (result$N[k]-1)/result$N[k]
  result[k,c("sqbias.bart","sqbias.mtap")] <- result[k,c("mse.bart","mse.mtap")]-result[k,c("var.bart","var.mtap")]
  
  qSE <- apply(X=dfObj[,c("se.bart","se.mtap")], MARGIN=2, FUN=quantile, c(0.98,0.02)) # quantiles for approximate CIs
  result[k,c("qSE.98.bart","qSE.98.mtap")] <- qSE[1,]
  result[k,c("qSE.02.bart","qSE.02.mtap")] <- qSE[2,]
  
  result[k,"smpl.cov"] <- cov(x=dfObj$b1, y=dfObj$b2)
}

result$CI.bart <- result$qSE.98.bart-result$qSE.02.bart
result$CI.mtap <- result$qSE.98.mtap-result$qSE.02.mtap


################################# plot results #################################
plotComment <- paste0(
  "Fixed ",
  ifelse(mtmFixed=="NW", paste0("NW=",timeBandProd,"; K=", numTapers),
           paste0("W=",W)
         ), "\n",
  ifelse(embedSines, "S", "No s"), "ines embedded\n",
  NUM_REGR, " regressions per N"
  )

# MSE plot
MSE.plotLims <- range( result[,grepl("qSE", names(result))] )
if (writeImgFile) {
  currTime <- gsub(" ", "-", gsub(":","", gsub("-","",Sys.time())))
  pdf(paste0("img/MSEbetahats_", currTime, ".pdf"), width=6, height=6)
}
layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow = TRUE))
par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
plot(x=result$N-5, y=result$mse.bart, xlab="Realization Size", ylab="MSE",
     log="y", ylim=MSE.plotLims*c(1e-1,1e1), col="goldenrod", pch=16)
points(x=result$N+5, y=result$mse.mtap, col="blue", pch=17)
arrows(x0=result$N-5, y0=result$qSE.02.bart, x1=result$N-5, y1=result$qSE.98.bart,
       length=0.05, angle=90, code=3, lwd=2, col="goldenrod")
arrows(x0=result$N+5, y0=result$qSE.02.mtap, x1=result$N+5, y1=result$qSE.98.mtap,
       length=0.05, angle=90, code=3, lwd=2, col="blue")
legend("topright", legend=c("Bartlett","Multitaper"), pch=c(16,17),
       col=c("goldenrod","blue"))
text(x=min(numObsVec), y=MSE.plotLims[1]*0.1, adj=c(0,0),
     labels=plotComment, family="mono")
text(x=numObsVec, y=result$qSE.02.mtap, pos=1,
     labels=sprintf("%.2f", result$CI.bart/result$CI.mtap), family="mono",
     col="blue")
text(x=max(numObsVec), y=MSE.plotLims[1]*0.1, adj=c(1,0), family="mono",
     labels="[ mtm rel. efficiency ]", col="blue")

par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
plot(x=result$N, y=rep(1,dim(result)[1]), xlab="Realization Size",
     ylab="Var. | Sq. Bias", ylim=c(0,1), col="white")
arrows(x0=result$N-5, y0=0, x1=result$N-5, y1=result$var.bart/result$mse.bart,
       length=0.05, angle=90, code=3, lwd=2, col="goldenrod")
arrows(x0=result$N-5, y0=result$var.bart/result$mse.bart, x1=result$N-5, y1=1,
       length=0.05, angle=90, code=3, lwd=2, col="darkgoldenrod")
arrows(x0=result$N+5, y0=0, x1=result$N+5, y1=result$var.mtap/result$mse.mtap,
       length=0.05, angle=90, code=3, lwd=2, col="blue")
arrows(x0=result$N+5, y0=result$var.mtap/result$mse.mtap, x1=result$N+5, y1=1,
       length=0.05, angle=90, code=3, lwd=2, col="darkblue")
if (writeImgFile) { dev.off() }






# variance plot
if (writeImgFile) {
  png("img/varBetahats.png", width=640, height=360)
}
par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
plot(x=result$N, y=result$var.bart, xlab="Realization Size", ylab="var",
     log="y", ylim=range(result[,c("var.bart","var.mtap")])*c(1e-1,1e1))
points(x=result$N, y=result$var.mtap, col="blue", pch=2)
legend("topright", legend=c("Bartlett","Multitaper"), pch=c(1,2),
       col=c("black","blue"))
text(x=min(numObsVec), y=min(result[,c("var.bart","var.mtap")])*0.1, adj=c(0,0),
     labels=plotComment)
if (writeImgFile) { dev.off() }



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





