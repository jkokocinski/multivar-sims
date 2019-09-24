#################################### psumCov ###################################
#
# psumCov : function to compute sum of A^k %*% V %*% t(A)^k, k from 0 to infty,
#   truncated such that successive term difference is less than epsilon.
#
psumCov <- function(A, V=diag(1,dim(A)[1])) {
  if (dim(A)[1]!=dim(A)[2]) stop("Not a square matrix.")
  ev <- eigen(x=A, only.values=TRUE)$values
  if (any(Mod(ev) > 1)) stop("Unstable AR model.")
  
  if (isSymmetric(A)) { # A is symmetric
    cat("A is symmetric.\n")
    if (dim(V)[1]==2) {
      if (V[1,1]==V[2,2] & V[1,2]==0 & V[2,1]==0) { # V is scalar mult. of I_2
        cat("V is scalar multiple of identity.\n")
        tryCatch(
          { return( solve(diag(1,2) - A%*%A) * V[1,1] ) },
          error=function(e) {cat("Unstable AR model.\n")}
          )
        }
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
#     * `toepVec` should be such that the order of the vector entries
#       corresponds to traversing down the last column of the Toeplitz matrix,
#       and then traversing down the first column of the Toeplitz matrix, but
#       excluding the duplicated "lag-0" entry.
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
# end toepLeftMult2


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
} # end bivAR1.ccvf


################################## bivARp.ccvf #################################
#
# bivARp.ccvf : Compute the CCVF for a bivariate AR(p) process, given the
#   coefficient matrices (as a 2x2p matrix) and the variance of the innovations.
#
bivARp.ccvf <- function(coefMats, V=diag(1,2), maxlag=1) {
  stopifnot(dim(coefMats)[1]==2, dim(coefMats)[2] %% 2 == 0)
  stopifnot(maxlag>0)
  
  p <- dim(coefMats)[2] / 2
  fullV <- matrix(0, nrow=2*p, ncol=2*p)
  fullV[(1:2),(1:2)] <- V
  
  # block coef matrix with I_2 blocks on the subdiagonal
  fullCoefMat <- rbind(coefMats, diag(1, 2*(p-1), 2*p))
  
  matrixSeries <- psumCov(fullCoefMat, V=fullV)
  
  CCVmats <- matrix(rep(matrixSeries,2*maxlag+1), nrow=2*p, ncol=2*p*(2*maxlag+1))
  for (h in seq(1,maxlag,1)) {
    startCol <- 2*p*maxlag+1 + 2*p*h # left-most column in matrix for lag +h
    CCVmats[,seq(startCol, startCol + 2*p-1)] <- fullCoefMat %*%
      CCVmats[,seq(startCol-2*p,startCol-1)]
    
    startCol <- dim(CCVmats)[2]-(startCol-1+2*p) +1  # left-most column in matrix for lag -h
    CCVmats[,seq(startCol, startCol + 2*p-1)] <- CCVmats[,seq(startCol+2*p,startCol+2*2*p-1)] %*%
      t(fullCoefMat)
  }
  
  CCVmats <- CCVmats[(1:2), sort(outer(seq(0,2*maxlag)*2*p, (1:2), "+"))]
  
  return(CCVmats)
} # end bivARp.ccvf

################################# bivAR1.ccvf.D ################################
#
# bivAR1.ccvf.D : Compute the CCVF for a bivariate AR(1) process, given the
#   coefficient matrix of the dependent AR(1) process, the diagonal matrix
#   coefficient, and the variance of the innovations.
# 
#              cov(y(t+h),y(t)) = D cov(x(t+h),x(t)) D^T + diag(V)
#
#     * `ccvfMats` is a matrix of submatrices (each 2×2), that are
#       autocovariance matrices at lags -h to h. This is the output from the
#       bivAR1.ccvf function, defined above.
#     * `D` is a diagonal matrix which left multiplies the input process to get
#       the output process, not including added noise (see `V`).
#     * `V` is a 2×2 covariance matrix for the bivariate noise which is added to
#       the result of left-multiplying the input by `D`, thus giving the final
#       output process.
#
bivAR1.ccvf.D <- function(ccvfMats, D, V) {
  stopifnot(dim(ccvfMats)[1]==2)
  
  CCVmats.D <- D %*% ccvfMats
  startColInds <- seq(1,dim(CCVmats.D)[2]-1,2)
  for (j in startColInds) {
    CCVmats.D[,((0:1)+j)] <- CCVmats.D[,((0:1)+j)] %*% D +
      as.numeric(2*j==dim(CCVmats.D)[2]) * V    # note that t(D)==D
  }
  
  return(CCVmats.D)
}
# end bivAR1.ccvf.D



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
#     * `errCovMat.r` is the var-cov matrix for the response innovations.
#     * `errCovMat.p` is the var-cov matrix for the predictor innovations.
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
# Return value is a list containing the list of regression coefficients and
#   their estimated covariances (for each realization size), as well as the
#   estimated MSEs of the regression coefficients in a data frame.
#
ar.regr.cov <- function(phiMat.p, phiMat.r, errCovMat.p, errCovMat.r,
                         numObsVec, NUM_REGR,
                         mtmFixed="NW", timeBandProd=4, numTapers=7, W,
                         adaptWt=FALSE,
                         writeImgFile=FALSE, embedSines=TRUE,
                         linDepY=FALSE, computeCorr=FALSE) {
  
  if(length(intersect(c("mAr", "multitaper"), rownames(installed.packages())))<2) {
    stop("Ensure `mAr` and `multitaper` packages are loaded.")
  }
  
  if (mtmFixed=="W") {
    # if just W is fixed
    if (missing(W)) {
      stop("W set to fixed, but value not specified.")
    }
    if (!missing(timeBandProd)) {
      cat("W fixed; timeBandProd not used.\n")
    }
    if (!missing(numTapers)) {
      cat("W fixed; numTapers not used.\n")
    }
  }
  if (linDepY) {
    cat("Responses 'linearly dependent' by design; \`phiMat.r\` not used.\n")
  }
  
  betasOverN <- list() # list of data.frames for the \hat{beta}s
  
  # headings for betas
  betas.head <- c("b1","b2","cv.bart","cv.theo","cv.mtap")
  if (computeCorr) {
    betas.head <- c(betas.head, c("cor.bart","cor.theo","cor.mtap"))
  }
  
  # \cov{Y_1,Y_2} from -(N-1) to +(N-1), where N is max(numObsVec); to subset
  maxmaxlag <- max(numObsVec)-1
  CCVmats.master <- bivARp.ccvf(coefMat=phiMat.r, V=errCovMat.r,
                                maxlag=maxmaxlag)
  
  ############################# do main computations #############################
  for (n in seq(1,length(numObsVec))) {
    numObs <- numObsVec[n]
    cat(paste0("################ N=",numObs,"\n"))
    
    maxlag <- numObs-1
    if (linDepY) {
      phi <- 0.95 # AR coef for Y_1
      a <- 2.5 # scale on Y_2
      varZ.Y1 <- 0.2
      addedNoiseVar.Y2 <- 0.05
      theo.acv1.r <- phi^(abs(seq(-maxlag, maxlag))) * varZ.Y1 / (1-phi^2)
      theo.acv2.r <- a^2 * theo.acv1.r + ((-maxlag:maxlag)==0) * addedNoiseVar.Y2
      theo.ccv.r <- a * theo.acv1.r
    } else {
      # subset CCVmats.master to go from -maxlag to +maxlag
      CCVmats <- CCVmats.master[,2*maxmaxlag + seq(-2*maxlag+1,2*maxlag+2,1)]
      theo.ccv.r <- CCVmats[1,2*seq(0,2*maxlag)+2]
      theo.acv1.r <- CCVmats[1,2*seq(0,2*maxlag)+1]
      theo.acv2.r <- CCVmats[2,2*seq(0,2*maxlag)+2]
    }
    
    # theo.ccv.r[which(theo.ccv.r<1e-18)] <- 0 # truncate at 1e-18
    # theo.ccv.r.mat <- matrix(nrow=numObs, ncol=numObs)
    # theo.ccv.r.mat <- matrix(theo.ccv.r[numObs +
    #                          row(theo.ccv.r.mat)-col(theo.ccv.r.mat)],numObs,numObs)
    
    # compute the Slepians
    # M <- 2^(trunc(log2(numObs))+2)
    if (mtmFixed=="NW") {
      sleps <- multitaper::dpss(numObs, k=numTapers, nw=timeBandProd) # fixed NW
    } else if (mtmFixed=="W") {
      sleps <- multitaper::dpss(numObs, k=trunc(2*numObs*W), nw=numObs*W) # fixed W
    } else {
      stop("mtmFixed must be either 'NW' or 'W'.")
    }
    
    K <- dim(sleps$v)[2] # number of tapers
    slepMat <- sleps$v
    
    # sample df of \hat{beta}s for regular processes and sine-embedded
    betas <- as.data.frame(matrix(nrow=NUM_REGR, ncol=length(betas.head)))
    names(betas) <- betas.head
    # betas.s <- betas
    
    # predictor series (X_1,X_2) realization
    bivAR1.p <- mAr.sim(w=rep(0,2), A=phiMat.p, C=errCovMat.p, N=numObs)
    
    if (embedSines) {
      # generate sinusoids
      sines <- matrix(nrow=numObs, ncol=2)
      sines[,1] <- 0.7*sin(2*pi/360*(1:numObs)) + 0.5*sin(2*pi/150*(1:numObs))
      sines[,2] <- sin(2*pi/240*(1:numObs)) + sin(2*pi/30*(1:numObs))
      sines <- 5 * sines %*% diag(diag(errCovMat.r)) # scale by resp. variances
      
      bivAR1.p <- bivAR1.p + sines # embed sinusoids in predictor series
    }
    
    # Moore-Penrose inverses
    mpi.X1 <- MASS::ginv(as.matrix(bivAR1.p[,1]))
    mpi.X2 <- MASS::ginv(as.matrix(bivAR1.p[,2]))
    
    cat(paste0("######## start : ",Sys.time()),"\n")
    for (j in 1:NUM_REGR) {
      if (linDepY) {
        y1 <- as.numeric(arima.sim(model=list(ar=phi), n=numObs, innov=rnorm(numObs, 0, sd=sqrt(varZ.Y1))))
        bivAR1.r <- data.frame( X1=y1, X2=a*y1 + rnorm(numObs, 0, sd=sqrt(addedNoiseVar.Y2)) )
      } else {
        bivAR1.r <- mAr.sim(w=rep(0,2), A=phiMat.r, C=errCovMat.r, N=numObs)
        # bivAR1.r <- as.matrix(bivAR1.p) %*% t(Dcoef) + MASS::mvrnorm(n=numObs, mu=rep(0,2), Sigma=errCovMat.r)
      }
      
      if (embedSines) {
        bivAR1.r <- bivAR1.r + sines # embed sinusoids in response series
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
      covB.bart <- mpi.X1 %*% toepLeftMult2( as.vector(bart.ccv.r$acf),
                                             as.vector(t(mpi.X2)) )
      # covB.bart <- mpi.X1 %*% bart.ccv.r.mat %*% t(mpi.X2) # brute force
      
      # \cov(\hat{\beta}_1, \hat{\beta}_2) -- theoretical-based
      covB.theo <- mpi.X1 %*% toepLeftMult2( theo.ccv.r, as.vector(t(mpi.X2)) )
      
      # compute the cross spectrum of (Y_1,Y_2) -- averaging of eigencoefs
      # Y1mtx <- matrix(bivAR1.r$X1, nrow=numObs, ncol=K)
      # Y1mtx <- mvfft(z=rbind(Y1mtx*slepMat, matrix(0,nrow=M-numObs, ncol=K)))
      # Y2mtx <- matrix(bivAR1.r$X2, nrow=numObs, ncol=K)
      # Y2mtx <- mvfft(z=rbind(Y2mtx*slepMat, matrix(0,nrow=M-numObs, ncol=K)))
      # crossSpecEstY <- (Y1mtx * Conj(Y2mtx)) %*% rep(1,K)/K
      
      # compute spec.mtm objects for both response series components
      spec.y1 <- multitaper::spec.mtm(timeSeries=ts(bivAR1.r$X1), k=K,
                                      ifelse(mtmFixed=="NW", timeBandProd, 2*numObs*W),
                                      adaptiveWeighting=adaptWt, dpssIN=slepMat,
                                      returnInternals=TRUE, plot=FALSE)
      spec.y2 <- multitaper::spec.mtm(timeSeries=ts(bivAR1.r$X2), k=K,
                                      ifelse(mtmFixed=="NW", timeBandProd, 2*numObs*W),
                                      adaptiveWeighting=adaptWt, dpssIN=slepMat,
                                      returnInternals=TRUE, plot=FALSE)
      
      # wieghts and eigencoefficients
      d1 <- spec.y1$mtm$eigenCoefWt
      d2 <- spec.y2$mtm$eigenCoefWt
      y1 <- spec.y1$mtm$eigenCoefs
      y2 <- spec.y2$mtm$eigenCoefs
      
      if (!adaptWt) {
        d2 <- d1 <- matrix(1, nrow=dim(y1)[1], ncol=dim(y1)[2])
      }
      
      # compute the cross spectrum of (Y_1,Y_2) -- adaptive weighted eigencoefs
      crossSpecEstY <- apply(d1*y1*d2*Conj(y2), MARGIN=1, FUN=sum) /
        apply(d1*d2, MARGIN=1, FUN=sum)
      
      # cross spectrum over full [0,1) interval
      crossSpecEstY <- c(crossSpecEstY, rev(Conj(crossSpecEstY[-1]))[-1])
      
      # CCVF of Y_1 and Y_2
      mtap.ccv.r <- fft(z=crossSpecEstY, inverse=TRUE) / length(crossSpecEstY)
      
      # put entries in "correct" order, from -numObs to +numObs
      mtap.ccv.r <- Re(c(tail(mtap.ccv.r, numObs-1), head(mtap.ccv.r, numObs)))
      
      # \cov(\hat{\beta}_1, \hat{\beta}_2) -- multitaper-based
      covB.mtap <- mpi.X1 %*% toepLeftMult2( mtap.ccv.r, as.vector(t(mpi.X2)) )
      # covB.mtap <- mpi.X1 %*% mtap.ccv.r.mat %*% t(mpi.X2)
      
      if (computeCorr) {
        # Bartlett ACVFs for responses; for use in calculating corr
        bart.acv1.r <- acf(x=bivAR1.r[,1], type="covariance", lag.max=numObs-1,
                           plot=FALSE)
        bart.acv2.r <- acf(x=bivAR1.r[,2], type="covariance", lag.max=numObs-1,
                           plot=FALSE)
        bart.acv1.r <- c(rev(bart.acv1.r$acf), bart.acv1.r$acf[-1])
        bart.acv2.r <- c(rev(bart.acv2.r$acf), bart.acv2.r$acf[-1])
        var.b1.bart <- mpi.X1 %*% toepLeftMult2( bart.acv1.r, as.vector(t(mpi.X1)) )
        var.b2.bart <- mpi.X2 %*% toepLeftMult2( bart.acv2.r, as.vector(t(mpi.X2)) )
        corB.bart <- covB.bart / sqrt(var.b1.bart * var.b2.bart)
        
        # theoretical correlations
        var.b1.theo <- mpi.X1 %*% toepLeftMult2( theo.acv1.r, as.vector(t(mpi.X1)) )
        var.b2.theo <- mpi.X2 %*% toepLeftMult2( theo.acv1.r, as.vector(t(mpi.X2)) )
        corB.theo <- covB.theo / sqrt(var.b1.theo * var.b2.theo)
        
        # autospectra for Y_1 and Y_2
        autoSpecEstY1 <- apply(d1*y1*d1*Conj(y1), MARGIN=1, FUN=sum) /
          apply(d1*d1, MARGIN=1, FUN=sum)
        autoSpecEstY2 <- apply(d2*y2*d2*Conj(y2), MARGIN=1, FUN=sum) /
          apply(d2*d2, MARGIN=1, FUN=sum)
        # autoSpecEstY1 <- (Y1mtx * Conj(Y1mtx)) %*% rep(1,K)/K
        # autoSpecEstY2 <- (Y2mtx * Conj(Y2mtx)) %*% rep(1,K)/K
        
        # autospectra over full [0,1) interval
        autoSpecEstY1 <- c(autoSpecEstY1, rev(Conj(autoSpecEstY1[-1]))[-1])
        autoSpecEstY2 <- c(autoSpecEstY2, rev(Conj(autoSpecEstY2[-1]))[-1])
        
        # ACVFs of Y_1 and Y_2, based on mtm
        mtap.acv1.r <- fft(z=autoSpecEstY1, inverse=TRUE) / length(autoSpecEstY1)
        mtap.acv2.r <- fft(z=autoSpecEstY2, inverse=TRUE) / length(autoSpecEstY2)
        mtap.acv1.r <- Re(c(tail(mtap.acv1.r, numObs-1), head(mtap.acv1.r, numObs)))
        mtap.acv2.r <- Re(c(tail(mtap.acv2.r, numObs-1), head(mtap.acv2.r, numObs)))
        
        var.b1.mtap <- mpi.X1 %*% toepLeftMult2( mtap.acv1.r, as.vector(t(mpi.X1)) )
        var.b2.mtap <- mpi.X2 %*% toepLeftMult2( mtap.acv2.r, as.vector(t(mpi.X2)) )
        corB.mtap <- covB.mtap / sqrt(var.b1.mtap * var.b2.mtap)
      }
      
      # update the betas data.frame
      betas.newRow <- c(as.numeric(model1$coefficients),
                        as.numeric(model2$coefficients),
                        covB.bart, covB.theo, covB.mtap)
      if (computeCorr) {
        betas.newRow <- c(betas.newRow, c(corB.bart, corB.theo, corB.mtap))
      }
      betas[j,] <- betas.newRow
                     
      # betas.s[j,] <- c(as.numeric(model1.s$coefficients),
      #                  as.numeric(model2.s$coefficients),
      #                  covB.bart.s, covB.theo.s)
    }
    cat(paste0("########   end : ",Sys.time()),"\n")
    
    if (!computeCorr) {
      betas$se.bart <- (with(betas, {cv.theo-cv.bart}))**2
      betas$se.mtap <- (with(betas, {cv.theo-cv.mtap}))**2
    } else {
      betas$se.bart <- (with(betas, {cor.theo-cor.bart}))**2
      betas$se.mtap <- (with(betas, {cor.theo-cor.mtap}))**2
    }
    
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
    
    if (!computeCorr) {
      result[k,c("var.bart","var.mtap")] <- apply(X=dfObj[,c("cv.bart","cv.mtap")],MARGIN=2,FUN=var) * (result$N[k]-1)/result$N[k]
    } else {
      result[k,c("var.bart","var.mtap")] <- apply(X=dfObj[,c("cor.bart","cor.mtap")],MARGIN=2,FUN=var) * (result$N[k]-1)/result$N[k]
    }
      
    result[k,c("sqbias.bart","sqbias.mtap")] <- result[k,c("mse.bart","mse.mtap")]-result[k,c("var.bart","var.mtap")]
    
    qSE <- apply(X=dfObj[,c("se.bart","se.mtap")], MARGIN=2, FUN=quantile, c(0.98,0.02)) # quantiles for approximate CIs
    result[k,c("qSE.98.bart","qSE.98.mtap")] <- qSE[1,]
    result[k,c("qSE.02.bart","qSE.02.mtap")] <- qSE[2,]
    
    result[k,"smpl.cov"] <- cov(x=dfObj$b1, y=dfObj$b2)
    result[k,"smpl.cor"] <- cor(x=dfObj$b1, y=dfObj$b2)
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
    NUM_REGR, " regressions per N",
    ifelse(linDepY, "\nY2=Y1+e", "")
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
  
  
  return(list(betasOverN=betasOverN, result=result, theo.ccv.r=theo.ccv.r,
              bart.ccv.r=bart.ccv.r, mtap.ccv.r=mtap.ccv.r, d1=d1, d2=d2,
              theo.acv1.r=theo.acv1.r))
  
}
# end ar1.regr.cov
