#################################### psumCov ###################################
#
# psumCov : function to compute sum of A^k %*% V %*% t(A)^k, k from 0 to infty,
#   truncated such that successive term difference is less than epsilon.
#
psumCov <- function(A, V=diag(1,dim(A)[1])) {
  if (dim(A)[1]!=dim(A)[2]) stop("Not a square matrix.")
  
  if (isSymmetric(A)) { # A is symmetric
    cat("A is symmetric.\n")
    if (V[1,1]==V[2,2] & V[1,2]==0 & V[2,1]==0) { # V is scalar mult. of I_2
      cat("V is scalar multiple of identity.\n")
      tryCatch(
        { return( solve(diag(1,2) - A%*%A) * V[1,1] ) },
        error=function(e) {cat("Unstable AR process.\n")}
        )
    }
  }
  
  cat("Using partial sum method.\n")
  psum <- V # initialize the sum; this is the first term
  epsilon <- 1
  iter <- 1
  while (epsilon>1e-16) { # threshold for convergence
    if (iter>1e6) {
      stop("Failed to converge after 1e6 terms.\n")
    }
    psum.old <- psum
    psum <- A %*% psum %*% t(A) + V
    epsilon <- max(abs(psum-psum.old))
    iter <- iter+1L
  }
  cat(paste0("Series converged after ", iter, " terms.\n"))
  return(psum)
}
# end psumCov


################################# toepLeftMult #################################
#
# toepLeftMult : Returns the result of a Toeplitz matrix left-multiplying a
#   vector. Let N be the number of rows and columns on the Toeplitz matrix.
# 
#     * `toepVec` is a vector of the matrix entries such that the first N
#       entries correspond the the first row of the matrix, entries 2 through
#       N+1 correspond to the second row of the matrix, etc. The length of this
#       vector is, thus, 2*N-1.
#     * `rightVec` is the vector to be left-multiplied.
# 
# The function works by performing a series of dot products, rather than having
#   to store a large matrix (stores 2*N-1 values in stead of N^2).
#
toepLeftMult <- function(toepVec,rightVec) {
  # toepVec should be such that the first N entries correspond with the last row
  #   of the Toeplitz matrix (this probably means th reversed ACVF vector).
  stopifnot(class(toepVec)=="numeric", class(rightVec)=="numeric")
  stopifnot(length(toepVec)==2*length(rightVec)-1)
  N <- length(rightVec)
  resVec <- matrix(nrow=N, ncol=1)
  for (j in 1:N) {
    resVec[j,1] <- rbind(toepVec[(N-j+1):(2*N-j)]) %*% cbind(rightVec)
  }
  return(resVec)
}
# end toepLeftMult


################################# toepLeftMult2 ################################
#
# toepLeftMult2 : Has the same function as toepLeftMult, but performs the matrix
#   multiplication by doing a multiplication in the frequency domain, instead of
#   convolving in the time domain, as in toepLeftMult.
#
#     * `toepVec` should be such that the first N entries of rev(toepVec)
#       correspond to the last row of the Toeplitz matrix
#       (If passing an ACVF vector, just use the regular ACVF vector with lags
#       going from negative to positive).
#
toepLeftMult2 <- function(toepVec,rightVec) {
  if (class(toepVec)!="numeric") { stop("toepVec must be a numeric vector.") }
  if (class(rightVec)!="numeric") { stop("rightVec must be a numeric vector.") }
  if (length(toepVec)!=2*length(rightVec)-1) {
    stop("lengths of args don't match")
  }
  N <- length(rightVec)
  resVec <- matrix(nrow=N, ncol=1)
  fft.tv <- fft(z=toepVec)
  fft.rv <- fft(z=c(rightVec, rep(0, N-1)))
  fft.res <- as.vector(fft.tv*fft.rv)
  resVec <- Re(fft(z=fft.res, inverse=TRUE)) / (2*N-1)
  resVec <- resVec[seq(N,2*N-1)]
  return(resVec)
}
# toepLeftMult2


################################## bivAR1.ccvf #################################
#
# bivAR1.ccvf : Compute the CCVF for a bivariate AR(1) process, given the
#   coefficient matrix and the variance of the innovations.
#
bivAR1.ccvf <- function(coefMat, V=diag(1,2), maxlag=1) {
  stopifnot(dim(coefMat)[1]==dim(coefMat)[2], dim(coefMat)[1]==2)
  stopifnot(maxlag>0)
  
  matrixSeries <- psumCov(coefMat, V=V)
  
  CCVmats <- matrix(rep(matrixSeries,2*maxlag+1), nrow=2, ncol=2*(2*maxlag+1))
  for (h in seq(1,maxlag,1)) {
    startCol <- 2*maxlag+1 + 2*h # left-most column in matrix for lag +h
    CCVmats[,seq(startCol,startCol+1)] <- coefMat %*%
      CCVmats[,seq(startCol-2,startCol-1)]
    
    startCol <- dim(CCVmats)[2]-startCol # left-most column in matrix for lag -h
    CCVmats[,seq(startCol,startCol+1)] <- CCVmats[,seq(startCol+2,startCol+3)] %*%
      t(coefMat)
  }
  return(CCVmats)
}


################################# ar1.regr.cov #################################
#
# ar1.regr.cov : Simulates bivariate AR(1) processes for predictor and response
#   time series realizations, and performs ordinary linear regression(s) to get
#   an estimated slope parameter. The covariance of this parameter is estimated
#   using Bartlett ACVF, and again using and inverted cross-spectral estimate of
#   the response series. The resulting MSEs are returned fro each realization
#   size, and are plotted with error bars.
#
#     * `phiMat.p` and `phiMat.r` are, respectively, the predictor and response
#       matrix coefficients for bivariate AR(1) realizations.
#     * `varZ` is the variance of the innovations for both the response and
#       predictor series. Default is 1.
#     * `numObsVec` is a numeric vector containing each realization size to be
#       used.
#     * `NUM_REGR` is the number of regressions to be performed per realization
#       size.
#     * `mtmFixed` should be set to either 'W' or 'NW', according to which
#       should be fixed in the cross-spectral computations.
#     * `M` is the length to zero-pad to in the cross-spectral estimates.
#     * `timeBandProd` and `numTapers` are the fixed time-bandwidth product and
#       number of tapers, respectively, set if mtmFixed=='NW'. Otherwise not
#       used.
#     * `W` is the band limit used in spec.mtm if mtmFixed=='W'. Otherwise
#       ignored.
#     * `writeImgFile` : should image file(s) be written?
#     * `embedSines` : should sinusoids be embedded in AR processes?
#
# The innovations covariance matrices are both the 2-by-2 identity matrix scaled
#   by `varZ`.
#
# Return value is a list containing the list of regression coefficients and
#   their estimated covariances (for each realization size), as well as the
#   estimated MSEs of the regression coefficients in a data frame.
#
ar1.regr.cov <- function(phiMat.p, phiMat.r, varZ=1, numObsVec, NUM_REGR,
                         mtmFixed="NW", timeBandProd=4, numTapers=7, W,
                         writeImgFile=FALSE, embedSines=FALSE) {
  library(mAr)
  library(multitaper)
  
  if (mtmFixed=="W" & missing(W)) {
    # if just W is fixed
    stop("W set to fixed, but value not specified.")
  }
  
  errCovMat.r <- errCovMat.p <- diag(varZ,2,2)
  
  betasOverN <- list() # list of data.frames for the \hat{beta}s
  betas.head <- c("b1","b2","cv.bart","cv.theo","cv.mtap") # headings for betas
  
  
  ############################# do main computations #############################
  for (n in seq(1,length(numObsVec))) {
    numObs <- numObsVec[n]
    
    # \cov{Y_1,Y_2} from -maxlag to +maxlag
    maxlag=numObs-1
    CCVmats <- bivAR1.ccvf(coefMat=phiMat.r, V=errCovMat.r, maxlag=maxlag)
    theo.ccv.r <- CCVmats[1,2*seq(0,2*maxlag)+2]
    # theo.ccv.r[which(theo.ccv.r<1e-18)] <- 0 # truncate at 1e-18
    # theo.ccv.r.mat <- matrix(nrow=numObs, ncol=numObs)
    # theo.ccv.r.mat <- matrix(theo.ccv.r[numObs +
    #                          row(theo.ccv.r.mat)-col(theo.ccv.r.mat)],numObs,numObs)
    
    # compute the Slepians
    M <- 6*2^(trunc(log2(numObs))+1)
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
      Y1mtx <- matrix(c((bivAR1.r[,1]-mean(bivAR1.r[,1])), rep(0,M-numObs)),
                      nrow=M, ncol=K)
      Y1mtx <- mvfft(z=Y1mtx*sleps$v)
      Y2mtx <- matrix(c((bivAR1.r[,2]-mean(bivAR1.r[,2])), rep(0,M-numObs)),
                      nrow=M, ncol=K)
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
  # layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow = TRUE))
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
  
  # par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
  # plot(x=result$N, y=rep(1,dim(result)[1]), xlab="Realization Size",
  #      ylab="Var. | Sq. Bias", ylim=c(0,1), col="white")
  # arrows(x0=result$N-5, y0=0, x1=result$N-5, y1=result$var.bart/result$mse.bart,
  #        length=0.05, angle=90, code=3, lwd=2, col="goldenrod")
  # arrows(x0=result$N-5, y0=result$var.bart/result$mse.bart, x1=result$N-5, y1=1,
  #        length=0.05, angle=90, code=3, lwd=2, col="darkgoldenrod")
  # arrows(x0=result$N+5, y0=0, x1=result$N+5, y1=result$var.mtap/result$mse.mtap,
  #        length=0.05, angle=90, code=3, lwd=2, col="blue")
  # arrows(x0=result$N+5, y0=result$var.mtap/result$mse.mtap, x1=result$N+5, y1=1,
  #        length=0.05, angle=90, code=3, lwd=2, col="darkblue")
  if (writeImgFile) { dev.off() }
  
  
  return(list(betasOverN=betasOverN, result=result))
  
}
# end ar1.regr.cov
