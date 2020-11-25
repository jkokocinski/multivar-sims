#################################### psumCov ###################################
#
# psumCov : function to compute sum of A^k %*% V %*% t(A)^k, k from 0 to infty,
#   truncated such that successive term difference is less than epsilon.
#
psumCov <- function(A, V) {
  if (dim(A)[1]!=dim(A)[2]) stop("Not a square matrix.")
  ev <- eigen(x=A, only.values=TRUE)$values
  if (any(Mod(ev) > 1)) stop("Unstable AR model.")
  
  # if (V=="default") V <- diag(c(rep(1,2), rep(0,dim(A)[2]-2)))
  
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
    if (iter>1e12) {
      stop("Failed to converge after 1e12 terms.\n")
    }
    if (iter %in% 10**c(3:11)) { cat(paste0(iter,"\n")) }
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


################################### bilinToep ##################################
#
# bilinToep : Computes the bilinear form of a Toeplitz matrix.
#
#     * `toepVec` should be as in the `toepLeftMult` function, having length
#        2*N-1.
#     * `leftVec` should be a numeric vector the pre-multiplies the Toeplitz
#        matrix.
#     * `rightVec` should be a numeric vector that post-multiplies the Toeplitz
#        matrix.
#
bilinToep <- function(leftVec,toepVec,rightVec) {
  if (class(leftVec)!="numeric") {
    if (class(leftVec)=="matrix") {
      if(any(dim(leftVec)==1)) {
        leftVec <- as.vector(leftVec)
      }
    }
    else {
      stop("leftVec must be a numeric vector.")
    }
  }
  if (class(toepVec)!="numeric") { stop("toepVec must be a numeric vector.") }
  if (class(rightVec)!="numeric") {
    if (class(rightVec)=="matrix") {
      if(any(dim(rightVec)==1)) {
        rightVec <- as.vector(rightVec)
      }
    }
    else {
      stop("rightVec must be a numeric vector.")
    }
  }
  if (length(toepVec)!=2*length(rightVec)-1) {
    stop("non-conformable arguments")
  }
  if (length(leftVec)!=length(rightVec)) {
    stop("non-conformable arguments")
  }
  
  N <- length(leftVec)
  fft.lv <- fft(z=rev(c(rep(0, N-1), leftVec)))
  fft.tv <- fft(z=toepVec)
  fft.rv <- fft(z=c(rightVec, rep(0, N-1)))
  resVec <- Re(fft(z=fft.lv*fft.tv*fft.rv, inverse=TRUE)) / (2*N-1)
  resVal <- tail(resVec, n=1)
  return(resVal)
}
# end bilinToep


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
  
  # matrixSeries <- psumCov(fullCoefMat, V=fullV)
  matrixSeries <- matrix(
    solve( diag(1, (2*p)**2) - kronecker(fullCoefMat, fullCoefMat) ) %*%
      as.vector(fullV),
    nrow=2*p, ncol=2*p)
  
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
                        embedSines=TRUE,
                        linDepY=FALSE, removeLCs=FALSE) {
  
  library(mAr) # depends on MASS
  library(multitaper)
  
  params.init.str <- ls() # vector of initally set parameters
  params.init <- eval(parse(text=paste0("list(",
    paste0(params.init.str,"=",params.init.str, collapse=",\n"),")")))
  
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
  betasOverN.s.w <- betasOverN.s <- betasOverN
  
  # headings for betas
  betas.head <- c("b1","b2","cv.bart","cv.theo","cv.mtap",
                  "cor.bart","cor.theo","cor.mtap",
                  "se.cv.bart","se.cv.mtap","se.cor.bart","se.cor.mtap")
  
  # \cov{Y_1,Y_2} from -(N-1) to +(N-1), where N is max(numObsVec); to subset
  maxmaxlag <- max(numObsVec)-1
  CCVmats.master <- bivARp.ccvf(coefMats=phiMat.r, V=errCovMat.r,
                                maxlag=maxmaxlag)
  
  
  for (n in seq(1,length(numObsVec))) {
    numObs <- numObsVec[n]
    cat(paste0("################ N=", numObs, " ",
               paste0(rep("#",60-nchar(numObs)), collapse=""), "\n"))
    
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
    bivAR1.p <- as.matrix(mAr.sim(w=rep(0,2), A=phiMat.p, C=errCovMat.p, N=numObs))
    
    freqsToEmbed <- (c(24,24))**(-1) # for sinusoids in components 1 & 2
    
    if (embedSines) {
      bivAR1.p.s <- embedSinusoids(input=bivAR1.p, freqs=freqsToEmbed,
                                   amps=sqrt(diag(errCovMat.p)), ampScale=3)
      # bivAR1.p.s <- embedSinusoids(input=bivAR1.p.s, freqs=c(12,12)**(-1),
      #                              amps=sqrt(diag(errCovMat.r)), ampScale=0.7)
      # bivAR1.p.s <- embedSinusoids(input=bivAR1.p.s, freqs=c(21,21)**(-1),
      #                              amps=sqrt(diag(errCovMat.r)), ampScale=0.2)
      # bivAR1.p.s <- embedSinusoids(input=bivAR1.p.s, freqs=c(20,30)**(-1),
      #                              amps=sqrt(diag(errCovMat.r)), ampScale=0.05)
    }
    
    # Moore-Penrose inverses of the predictor realization vectors
    mpi.X1 <- MASS::ginv(bivAR1.p[,1])
    mpi.X2 <- MASS::ginv(bivAR1.p[,2])
    if (embedSines) {
      mpi.X1.s <- MASS::ginv(bivAR1.p.s[,1])
      mpi.X2.s <- MASS::ginv(bivAR1.p.s[,2])
      
      # if (removeLCs) { # deterministic component estimation
      #   seas.x1 <- determineSeasonal(data=bivAR1.p.s[,1], sigCutoff=0.999)
      #   seas.x2 <- determineSeasonal(data=bivAR1.p.s[,2], sigCutoff=0.999)
      # }
    }
      
    # matrices to store CCVF estimates
    bart.ccv.r.mat <- matrix(nrow=NUM_REGR, ncol=2*numObs-1)
    mtap.ccv.r.mat <- matrix(nrow=NUM_REGR, ncol=2*numObs-1)
    
    # initizalize array of response realizations
    respRlzns <- array(data=NA, dim=c(numObs, 2, NUM_REGR))
    respRlzns.s.w <- respRlzns.s <- respRlzns # .s sines embedded; .s.w whitened
    
    ############################# generate realizations ############################
    cat(paste0("# ", Sys.time()) ," | generating realizations...\n")
    pb1 <- txtProgressBar(min=0, max=NUM_REGR, initial=0, char=">", width=80)
    for (j in 1:NUM_REGR) {
      if (linDepY) {
        y1 <- as.numeric(arima.sim(model=list(ar=phi), n=numObs, innov=rnorm(numObs, 0, sd=sqrt(varZ.Y1))))
        bivAR1.r <- cbind( y1, a*y1 + rnorm(numObs, 0, sd=sqrt(addedNoiseVar.Y2)) )
      } else {
        bivAR1.r <- as.matrix(
          mAr.sim(w=rep(0,2), A=phiMat.r, C=errCovMat.r, N=numObs)
        )
        # bivAR1.r <- as.matrix(bivAR1.p) %*% t(Dcoef) + MASS::mvrnorm(n=numObs, mu=rep(0,2), Sigma=errCovMat.r)
      }
      
      respRlzns[,,j] <- bivAR1.r
      
      if (embedSines) {
        bivAR1.r.s <- embedSinusoids(input=bivAR1.r, freqs=freqsToEmbed,
                                     amps=sqrt(diag(errCovMat.r)), ampScale=2)
        # bivAR1.r.s <- embedSinusoids(input=bivAR1.r.s, freqs=c(12,12)**(-1),
        #                              amps=diag(errCovMat.r), ampScale=0.25)
        # bivAR1.r.s <- embedSinusoids(input=bivAR1.r.s, freqs=c(21,21)**(-1),
        #                              amps=diag(errCovMat.r), ampScale=0.1)
        respRlzns.s[,,j] <- bivAR1.r.s
        
        # detect & remove sinusoidal line components
        if (removeLCs) {
          for (p in 1:2) {
            commonSinesObj <- findCommonSines(x=bivAR1.p.s[,p], y=bivAR1.r.s[,p],
              freqThresh=ifelse(mtmFixed=="NW", timeBandProd / numObs, W),
              sigCutoff=1-1/numObs, NW=timeBandProd, K=numTapers)
            detRespSines <- commonSinesObj$fctVals.com.y
            if(is.null(detRespSines)) { detRespSines <- matrix(0,numObs,2) } # to avoid NULLs
            
            respRlzns.s.w[,p,j] <- bivAR1.r.s[,p] - apply(detRespSines, 1, sum)
          }
          
        } # end if(removeLCs)
      } # end if(embedSines)
      setTxtProgressBar(pb=pb1, value=j)
    } # end generate realizations
    setTxtProgressBar(pb=pb1, value=0)
    cat(paste0("# ", Sys.time()) ," | done!\n")
  
  
  
    
    ######################### do regressions and estimation ########################
    cat(paste0("# ", Sys.time()) ," | doing regressions and estimation...\n")
    pb2 <- txtProgressBar(min=0, max=NUM_REGR*(1+embedSines+removeLCs), initial=0, char=">", width=80)
    for (tp in seq(1, 1 + embedSines + removeLCs, 1)) {
      if (tp==1) {
        Yarr <- respRlzns
        predRlzn <- bivAR1.p
        A1 <- mpi.X1
        A2 <- mpi.X2
      } else if (tp==2) {
        Yarr <- respRlzns.s
        predRlzn <- bivAR1.p.s
        A1 <- mpi.X1.s
        A2 <- mpi.X2.s
      } else if (tp==3) {
        Yarr <- respRlzns.s.w # used for the Bartlett est.; YforSpec are respRlzns
        predRlzn <- bivAR1.p
        A1 <- mpi.X1
        A2 <- mpi.X2
      } else { stop("impossible tp") }
      
      for (j in 1:NUM_REGR) {
        Y <- Yarr[,,j]
        # regular lin regs
        if (tp==4) { # do I want to do this for tp==3? Disabled for now.
          model1 <- lm(respRlzns.s.w[,1,j] ~ 0 + predRlzn[,1])
          model2 <- lm(respRlzns.s.w[,2,j] ~ 0 + predRlzn[,2])
        } else {
          model1 <- lm(Y[,1] ~ 0 + predRlzn[,1])
          model2 <- lm(Y[,2] ~ 0 + predRlzn[,2])
        }
        
        # Bartlett CCVF of response series
        bart.ccv.r <- ccf(x=Y[,1], y=Y[,2], type="covariance",
                          lag.max=numObs-1, plot=FALSE)$acf
        
        # \cov(\hat{\beta}_1, \hat{\beta}_2) -- Bartlett-based
        covB.bart <- A1 %*% toepLeftMult2( as.vector(bart.ccv.r),
                                               as.vector(t(A2)) )
        
        # \cov(\hat{\beta}_1, \hat{\beta}_2) -- theoretical-based
        covB.theo <- A1 %*% toepLeftMult2( theo.ccv.r, as.vector(t(A2)) )
        
        if (tp==3) {
          YforSpec <- respRlzns.s.w[,,j]
        } else {
          YforSpec <- Y
        }
        
        # compute spec.mtm objects for both response series components
        spec.y1 <- multitaper::spec.mtm(
          timeSeries=ts(YforSpec[,1]), k=K,
          nw=ifelse(mtmFixed=="NW", timeBandProd, 2*numObs*W),
          adaptiveWeighting=adaptWt, dpssIN=slepMat,
          Ftest=TRUE, returnInternals=TRUE, plot=FALSE
        )
        spec.y2 <- multitaper::spec.mtm(
          timeSeries=ts(YforSpec[,2]), k=K,
          nw=ifelse(mtmFixed=="NW", timeBandProd, 2*numObs*W),
          adaptiveWeighting=adaptWt, dpssIN=slepMat,
          Ftest=TRUE, returnInternals=TRUE, plot=FALSE
        )
        
        # weights and eigencoefficients
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
        covB.mtap <- A1 %*% toepLeftMult2( mtap.ccv.r, as.vector(t(A2)) )
        # covB.mtap <- A1 %*% mtap.ccv.r.mat %*% t(A2)
        
        if (TRUE) {
          # Bartlett ACVFs for responses; for use in calculating corr
          bart.acv1.r <- acf(x=Y[,1], type="covariance", lag.max=numObs-1,
                             plot=FALSE)
          bart.acv2.r <- acf(x=Y[,2], type="covariance", lag.max=numObs-1,
                             plot=FALSE)
          bart.acv1.r <- c(rev(bart.acv1.r$acf), bart.acv1.r$acf[-1])
          bart.acv2.r <- c(rev(bart.acv2.r$acf), bart.acv2.r$acf[-1])
          var.b1.bart <- A1 %*% toepLeftMult2( bart.acv1.r, as.vector(t(A1)) )
          var.b2.bart <- A2 %*% toepLeftMult2( bart.acv2.r, as.vector(t(A2)) )
          corB.bart <- covB.bart / sqrt(var.b1.bart * var.b2.bart)
          
          # theoretical correlations
          var.b1.theo <- A1 %*% toepLeftMult2( theo.acv1.r, as.vector(t(A1)) )
          var.b2.theo <- A2 %*% toepLeftMult2( theo.acv2.r, as.vector(t(A2)) )
          corB.theo <- covB.theo / sqrt(var.b1.theo * var.b2.theo)
          
          # autospectra for Y_1 and Y_2
          autoSpecEstY1 <- apply(d1*y1*d1*Conj(y1), MARGIN=1, FUN=sum) /
            apply(d1*d1, MARGIN=1, FUN=sum)
          autoSpecEstY2 <- apply(d2*y2*d2*Conj(y2), MARGIN=1, FUN=sum) /
            apply(d2*d2, MARGIN=1, FUN=sum)
          
          # autospectra over full [0,1) interval
          autoSpecEstY1 <- c(autoSpecEstY1, rev(Conj(autoSpecEstY1[-1]))[-1])
          autoSpecEstY2 <- c(autoSpecEstY2, rev(Conj(autoSpecEstY2[-1]))[-1])
          
          # ACVFs of Y_1 and Y_2, based on mtm
          mtap.acv1.r <- fft(z=autoSpecEstY1, inverse=TRUE) / length(autoSpecEstY1)
          mtap.acv2.r <- fft(z=autoSpecEstY2, inverse=TRUE) / length(autoSpecEstY2)
          mtap.acv1.r <- Re(c(tail(mtap.acv1.r, numObs-1), head(mtap.acv1.r, numObs)))
          mtap.acv2.r <- Re(c(tail(mtap.acv2.r, numObs-1), head(mtap.acv2.r, numObs)))
          
          var.b1.mtap <- A1 %*% toepLeftMult2( mtap.acv1.r, as.vector(t(A1)) )
          var.b2.mtap <- A2 %*% toepLeftMult2( mtap.acv2.r, as.vector(t(A2)) )
          corB.mtap <- covB.mtap / sqrt(var.b1.mtap * var.b2.mtap)
        }
        
        # update the betas data.frame
        betas.newRow <- c(as.numeric(model1$coefficients),
                          as.numeric(model2$coefficients),
                          covB.bart, covB.theo, covB.mtap)
        betas.newRow <- c(betas.newRow, c(corB.bart, corB.theo, corB.mtap),
                          rep(NA,4))
        betas[j,] <- betas.newRow
        
        # CCVF estimate averages (without dividing by NUM_REGR yet)
        bart.ccv.r.mat[j,] <- bart.ccv.r
        mtap.ccv.r.mat[j,] <- mtap.ccv.r
        
        setTxtProgressBar(pb=pb2, value=(tp-1)*NUM_REGR + j)
      } # end for (j in 1:NUM_REGR)
      setTxtProgressBar(pb=pb2, value=0)
    
      betas$se.cv.bart <- (with(betas, {cv.theo-cv.bart}))**2
      betas$se.cv.mtap <- (with(betas, {cv.theo-cv.mtap}))**2
      betas$se.cor.bart <- (with(betas, {cor.theo-cor.bart}))**2
      betas$se.cor.mtap <- (with(betas, {cor.theo-cor.mtap}))**2
      
      # update the list of dfs
      if (tp==1) {
        betasOverN[[n]] <- betas
        
        bart.ccv.r.last <- bart.ccv.r
        mtap.ccv.r.last <- mtap.ccv.r
        bart.ccv.r.ave <- apply(bart.ccv.r.mat, 2, mean)
        mtap.ccv.r.ave <- apply(mtap.ccv.r.mat, 2, mean)
        B.ccv.r.mat <- bart.ccv.r.mat
        M.ccv.r.mat <- mtap.ccv.r.mat
        crossSpec <- crossSpecEstY
      } else if (tp==2) {
        betasOverN.s[[n]] <- betas
        
        bart.ccv.r.s.last <- bart.ccv.r
        mtap.ccv.r.s.last <- mtap.ccv.r
        bart.ccv.r.s.ave <- apply(bart.ccv.r.mat, 2, mean)
        mtap.ccv.r.s.ave <- apply(mtap.ccv.r.mat, 2, mean)
        B.ccv.r.s.mat <- bart.ccv.r.mat
        M.ccv.r.s.mat <- mtap.ccv.r.mat
        crossSpec.s <- crossSpecEstY
      } else if (tp==3) {
        betasOverN.s.w[[n]] <- betas
        
        bart.ccv.r.s.w.last <- bart.ccv.r
        mtap.ccv.r.s.w.last <- mtap.ccv.r
        bart.ccv.r.s.w.ave <- apply(bart.ccv.r.mat, 2, mean)
        mtap.ccv.r.s.w.ave <- apply(mtap.ccv.r.mat, 2, mean)
        B.ccv.r.s.w.mat <- bart.ccv.r.mat
        M.ccv.r.s.w.mat <- mtap.ccv.r.mat
        crossSpec.s.w <- crossSpecEstY
      } else {
        stop("impossible tp")
      }
      
    } # end for (tp)
    cat(paste0("# ", Sys.time()) ," | done!\n")
  } # end for (n in seq(1,length(numObsVec)))
  
  
  
  #################################### results ###################################
  # produce result df with MSEs, etc.
  for (tp in seq(1, 1 + embedSines + removeLCs, 1)) {
    result <- data.frame(N=numObsVec)
    
    if (tp==1) {
      dfObjList <- betasOverN
    } else if (tp==2) {
      dfObjList <- betasOverN.s
    } else if (tp==3) {
      dfObjList <- betasOverN.s.w
    } else {
      stop("impossible tp")
    }
    
    for (k in seq(1,length(betasOverN))) {
      dfObj <- dfObjList[[k]]
      
      # NOTE: The following three blocks, which do multiple assignments in a single
      #   line do not work on some older versions of R. Works on 3.6.0. Should fix
      #   this.
      
      # estimated MSEs of the cov & cor estimates
      result[k,c("mse.cv.bart","mse.cv.mtap")] <-
        apply(X=dfObj[,c("se.cv.bart","se.cv.mtap")],MARGIN=2,FUN=mean)
      result[k,c("mse.cor.bart","mse.cor.mtap")] <-
        apply(X=dfObj[,c("se.cor.bart","se.cor.mtap")],MARGIN=2,FUN=mean)
      
      # estimated variances of the cov & cor estimates
      result[k, c("var.cv.bart", "var.cv.mtap")] <-
        apply(X=dfObj[,c("cv.bart","cv.mtap")], MARGIN=2, FUN=var) *
        (result$N[k]-1) / result$N[k]
      result[k,c("var.cor.bart","var.cor.mtap")] <-
        apply(X=dfObj[,c("cor.bart","cor.mtap")],MARGIN=2,FUN=var) *
        (result$N[k]-1)/result$N[k]
      
      # estimated squared bias of the cov & cor estimates; MSE subtract var
      result[k,c("sqbias.cv.bart","sqbias.cv.mtap")] <-
        result[k,c("mse.cv.bart","mse.cv.mtap")] -
        result[k,c("var.cv.bart","var.cv.mtap")]
      result[k,c("sqbias.cor.bart","sqbias.cor.mtap")] <-
        result[k,c("mse.cor.bart","mse.cor.mtap")] -
        result[k,c("var.cor.bart","var.cor.mtap")]
      
      # sample quantiles for the squared errors; approximate CIs
      lp <- 0.025; up <- 0.975 # lower and upper probs.
      qSE.cv <- apply(X=dfObj[,c("se.cv.bart","se.cv.mtap")], MARGIN=2,
                      FUN=quantile, c(up,lp))
      result[k,c("qSE.U.cov.bart","qSE.U.cov.mtap")] <- qSE.cv[1,]
      result[k,c("qSE.L.cov.bart","qSE.L.cov.mtap")] <- qSE.cv[2,]
      qSE.cor <- apply(X=dfObj[,c("se.cor.bart","se.cor.mtap")], MARGIN=2,
                       FUN=quantile, c(up,lp))
      result[k,c("qSE.U.cor.bart","qSE.U.cor.mtap")] <- qSE.cor[1,]
      result[k,c("qSE.L.cor.bart","qSE.L.cor.mtap")] <- qSE.cor[2,]
      
      result[k,"smpl.cov"] <- cov(x=dfObj$b1, y=dfObj$b2)
      result[k,"smpl.cor"] <- cor(x=dfObj$b1, y=dfObj$b2)
    }
    
    result$CI.cov.bart <- result$qSE.U.cov.bart-result$qSE.L.cov.bart
    result$CI.cov.mtap <- result$qSE.U.cov.mtap-result$qSE.L.cov.mtap
    result$CI.cor.bart <- result$qSE.U.cor.bart-result$qSE.L.cor.bart
    result$CI.cor.mtap <- result$qSE.U.cor.mtap-result$qSE.L.cor.mtap
    
    if (tp==1) {
      bhCovCIs <- result
    } else if (tp==2) {
      bhCovCIs.s <- result
    } else if (tp==3) {
      bhCovCIs.s.w <- result
    } else {
      stop("impossible tp")
    }
  } # end for tp (results)
  rm(result, betas,betas.head, betas.newRow)
  # rm(list=ls()[which(grepl("bivAR1", ls()))])
  
  localVars <- setdiff(ls(), names(params.init))
  
  returnExpr <- paste0("params.init=params.init,\n",
                       paste0(localVars,"=",localVars, collapse=",\n")
  )
  returnExpr <- paste0("return( list(",returnExpr, ") )")
  # cat(paste0(returnExpr,"\n"))
  
  # return statement
  eval(parse(text=returnExpr))
  
}
# end ar1.regr.cov





################################ embedSinusoids ################################
#
# embedSinusoids : Takes an input bivariate time series and adds a sinusoidal
#   function in each component, with frequencies defined by `freq` and
#   amplitudes defined by `amps`, each numeric vectors of length 2.
#
#     * `input` should be a matrix with 2 columns.
#     * `ampScale` is a scalar to further adjust both amplitudes.
#
embedSinusoids <- function(input, freqs, amps, ampScale) {
  stopifnot(class(input)=="matrix")
  stopifnot(dim(input)[2]==2, length(freqs)==2, length(amps)==2)
  
  numObs <- dim(input)[1]
  
  sines <- matrix(nrow=dim(input)[1], ncol=2)
  sines[,1] <- cos(2*pi*freqs[1]*(1:numObs))
  sines[,2] <- cos(2*pi*freqs[2]*(1:numObs))
  sines <- ampScale * sines %*% diag(amps)
  
  return(input + sines)
} # end embedSinusoids


################################ findCommonSines ###############################
#
# findCommonSines : Run the determineSeasonal function on two time series to
#   find line components in each, then identify those that are 'common' in terms
#   of frequency. Includes jackknifing of amplitude and phase parameters though
#   the determineSeasonal function (by default it is turned off here).
#
findCommonSines <- function(x, y, padFactor=7, freqThresh, sigCutoff,
                            NW=4, K=7, jackknife=FALSE) {
  stopifnot(length(x)==length(y))
  
  source("seasonalFunctions.R")
  
  suppressWarnings(
    {
      seas.x <- determineSeasonal(data=x, padFactor=padFactor, sigCutoff=sigCutoff, NW=NW, K=K, jackknife=jackknife)
      seas.y <- determineSeasonal(data=y, padFactor=padFactor, sigCutoff=sigCutoff, NW=NW, K=K, jackknife=jackknife)
    }
  )
  
  if (is.null(seas.x) | is.null(seas.y)) {
    return( list(fctVals=matrix(0, length(x), 2), paramsX=NULL, paramsY=NULL) )
  }
  
  # commonFreqIndsXY is a matrix where the (i,j)-th entry is TRUE iff the
  #   i-th frequency indentified in seas1 is within the threshold of the
  #   j-th frequency indentified in seas2; threshold defined by `thresh`.
  commonFreqIndsXY <- outer(
    X = seas.x$phaseAmplitudeInfo$freq,
    Y = seas.y$phaseAmplitudeInfo$freq,
    FUN = function(x, y, thresh=freqThresh) {
      abs(x-y) < thresh
    }
  )
  commonFreqInds <- which(commonFreqIndsXY, arr.ind=TRUE)
  incohFreqInds.x <- setdiff(seq(1,dim(commonFreqIndsXY)[1]), commonFreqInds[,1])
  incohFreqInds.y <- setdiff(seq(1,dim(commonFreqIndsXY)[2]), commonFreqInds[,2])
  
  if (dim(commonFreqInds)[1]==0) { # nothing common found
    allSeas.incoh.x <- apply(X=as.matrix(seas.x$sinusoidData), MARGIN=1, FUN=sum)
    allSeas.incoh.y <- apply(X=as.matrix(seas.y$sinusoidData), MARGIN=1, FUN=sum)
    return( list(fctVals.com.x=NULL, fctVals.com.y=NULL,
                 paramsX.com=NULL, paramsY.com=NULL,
                 fctVals.incoh=matrix(c(allSeas.incoh.x, allSeas.incoh.y), ncol=2))
          )
  } else {
    if (length(incohFreqInds.x)<1) {
      allSeas.incoh.x <- rep(0, length(x))
    } else {
      allSeas.incoh.x <- apply(X=as.matrix(seas.x$sinusoidData[,incohFreqInds.x]), MARGIN=1, FUN=sum)
    }
    
    if (length(incohFreqInds.y)<1) {
      allSeas.incoh.y <- rep(0, length(y))
    } else {
      allSeas.incoh.y <- apply(X=as.matrix(seas.y$sinusoidData[,incohFreqInds.y]), MARGIN=1, FUN=sum)
    }
    
    fctVals.com.x <- as.matrix(seas.x$sinusoidData[,commonFreqInds[,1]]) # incl. possible doubles at this point
    fctVals.com.y <- as.matrix(seas.y$sinusoidData[,commonFreqInds[,2]])
    paramsX.com <- seas.x$phaseAmplitudeInfo[commonFreqInds[,1],]
    paramsY.com <- seas.y$phaseAmplitudeInfo[commonFreqInds[,2],]
    fctVals.incoh <- matrix(c(allSeas.incoh.x, allSeas.incoh.y), ncol=2)
    
    # phase-aligned common sinusoids in `x` series
    fctVals.com.x.ph <- matrix(seq(0,length(x)-1,1), nrow=length(x), ncol=dim(fctVals.com.x)[2])
    fctVals.com.x.ph <- 2*pi*fctVals.com.x.ph %*%
      diag(paramsX.com$freq, nrow=dim(paramsX.com)[1], ncol=dim(paramsX.com)[1])
    fctVals.com.x.ph <-
      cos(fctVals.com.x.ph + matrix(
        paramsY.com$phase,
        nrow = length(x),
        ncol = dim(fctVals.com.x)[2],
        byrow=TRUE
      )) %*% diag(paramsX.com$amp, nrow=dim(paramsX.com)[1], ncol=dim(paramsX.com)[1])
    
    keepInds.x <- which(!duplicated(commonFreqInds[,1]))
    keepInds.y <- which(!duplicated(commonFreqInds[,2]))
    fctVals.com.x <- as.matrix(fctVals.com.x[,keepInds.x]) # dupes now removed
    fctVals.com.y <- as.matrix(fctVals.com.y[,keepInds.y])
    paramsX.com <- paramsX.com[keepInds.x,]
    paramsY.com <- paramsY.com[keepInds.y,]
    fctVals.com.x.ph <- as.matrix(fctVals.com.x.ph[,keepInds.x])
    
    # jackknifed amps and phases
    if (jackknife) {
      paramsX.com.jk <- seas.x$phaseAmplitudeInfo.jk[commonFreqInds[,1],,]
      paramsY.com.jk <- seas.y$phaseAmplitudeInfo.jk[commonFreqInds[,2],,]
      
      # below if statements do not execute if params.*.com.jk is a matrix, not an array
      if (length(dim(paramsX.com.jk))==3) {
        paramsX.com.jk <- paramsX.com.jk[keepInds.x,,] # removing dupes
      }
      if (length(dim(paramsY.com.jk))==3) {
        paramsY.com.jk <- paramsY.com.jk[keepInds.y,,]
      }
    } else {
      paramsX.com.jk <- NULL
      paramsY.com.jk <- NULL
    }
    
    return( list(fctVals.com.x=fctVals.com.x,
                 fctVals.com.y=fctVals.com.y,
                 paramsX.com=paramsX.com,
                 paramsY.com=paramsY.com,
                 fctVals.incoh=fctVals.incoh,
                 fctVals.com.x.ph=fctVals.com.x.ph,
                 paramsX.com.jk=paramsX.com.jk,
                 paramsY.com.jk=paramsY.com.jk)
          )
  }
} # end findCommonSines



################################### plotCCVF ###################################
#
# Plot the theoretical CCVF for the specified `plotLags` vector of integers.
#
#     * `stage` : A string specifying the suffix on the relevant data objects.
#        A null string, "", specifies that just the generated AR realizations
#        should be used; ".s" for sinusoids embedded; ".s.w" for sinusoids
#        embedded and then removed using Thompson's F-test. (Note that sinusoids
#        are not removed from for any of the Bartlett estimation.)
#     * `ave` : A logical; should the average estimated CCVFs be plotted? If
#       FALSE, the CCVF estimates from the last generated realizations are used.
#
plotCCVF <- function(resultList, plotLags, stage="", ave=FALSE) {
  if (stage=="") {
    suffix <- ""
  } else if (stage=="s") {
    suffix <- ".s"
  } else if (stage=="s.w") {
    suffix <- ".s.w"
  } else {
    stop("Set a valid `stage`.")
  }
  
  par(mar=c(4,4,1,1))
  plotText <- 
    "with(resultList,
       {
         N <- max(bhCovCIs$N)
         plot(x=plotLags, y=theo.ccv.r[plotLags+N], type=\"h\", lwd=2,
              ylim=range( c(theo.ccv.r[plotLags+N], mtap.ccv.r.last[plotLags+N],
                            bart.ccv.r.last[plotLags+N]) ),
              # ylim=range( c(theo.ccv.r[plotLags+N], mtap.ccv.r.last[plotLags+N]) ),
              ylab=\"CCVF (Y1,Y2)\", xlab=\"Lag\")
         abline(h=0)
         points(x=plotLags, y=mtap.ccv.r.last[plotLags+N], col=\"blue\")
         points(x=plotLags, y=bart.ccv.r.last[plotLags+N], col=\"red\")
       }
  )"
  
  plotText <- gsub("bhCovCIs", paste0("bhCovCIs", suffix), plotText)
  plotText <- gsub("bart.ccv.r", paste0("bart.ccv.r",suffix), plotText)
  plotText <- gsub("mtap.ccv.r", paste0("mtap.ccv.r",suffix), plotText)
  if (ave) {
    plotText <- gsub(".last", ".ave", plotText)
  }
  
  eval(parse(text=plotText))
} # end plotCCVF


#################################### plotCIs ###################################
#
# Plot the estimated confidence intervals for the Bartlett- and MTM-based
#   estimator of the covariance between the \hat{\beta}_0's. The endpoints of
#   these CIs are based on sample quantiles (hard-coded).
#
plotCIs <- function(resultList, type="cov", stage="", writeImgFile=FALSE) {
  with(resultList,
    {
      if (stage=="") {
        result <- bhCovCIs
      } else if (stage=="s") {
        result <- bhCovCIs.s
      } else if (stage=="s.w") {
        result <- bhCovCIs.s.w
      } else {
        stop("Set a valid `stage`.")
      }
        
      # plotComment <- paste0(
      #   "Fixed ",
      #   ifelse(mtmFixed=="NW", paste0("NW=",timeBandProd,"; K=", numTapers),
      #            paste0("W=",W)
      #          ), "\n",
      #   ifelse(embedSines, "S", "No s"), "ines embedded",
      #   ifelse(removeLCs, "; LCs not removed", ""), "\n",
      #   NUM_REGR, " regressions per N",
      #   ifelse(linDepY, "\nY2=Y1+e", "")
      #   )
      
      # MSE plot
      MSE.plotLims <- range(
        result[,(grepl("qSE", names(result)) & grepl(type, names(result)))]
      )
      if (writeImgFile) {
        currTime <- gsub(" ", "-", gsub(":","", gsub("-","",Sys.time())))
        pdf(paste0("img/MSEbetahats_", currTime, ".pdf"), width=7, height=4)
      }
      plotText <- "plot(x=result$N-5, y=result$mse.cv.bart, xlab=\"Realization Size\",
           ylab=\"SE of cov estimator\",
           log=\"y\", ylim=MSE.plotLims*10**c(-0.5,0.5), col=\"red\", pch=16)
      points(x=result$N+5, y=result$mse.cv.mtap, col=\"blue\", pch=17)
      arrows(x0=result$N-5, y0=result$qSE.L.cov.bart, x1=result$N-5, y1=result$qSE.U.cov.bart,
             length=0.05, angle=90, code=3, lwd=2, col=\"red\")
      arrows(x0=result$N+5, y0=result$qSE.L.cov.mtap, x1=result$N+5, y1=result$qSE.U.cov.mtap,
             length=0.05, angle=90, code=3, lwd=2, col=\"blue\")
      legend(\"topright\", legend=c(\"Bartlett\",\"Multitaper\"), pch=c(16,17),
             col=c(\"red\",\"blue\"), horiz=TRUE)
      # text(x=min(params.init$numObsVec), y=par()$yaxp[1], adj=c(0,0),
      #      labels=plotComment, family=\"mono\")
      text(x=params.init$numObsVec, y=result$qSE.L.cov.mtap, pos=1, family=\"mono\", col=\"blue\",
           labels=sprintf(\"%.2f\", result$CI.cov.bart/result$CI.cov.mtap))
      text(x=max(params.init$numObsVec), y=par()$yaxp[1], adj=c(0,0), family=\"mono\",
           labels=\"[ mtm rel. efficiency ]\", col=\"blue\")"
      
      if (type=="cor") {
        plotText <- gsub("cov", "cor", plotText)
        plotText <- gsub("cv", "cor", plotText)
      }
      
      # layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow = TRUE))
      par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
      eval(parse(text=plotText))
      
      if (writeImgFile) { dev.off() }
    }
  )
} # end plotCIs



################################# plotCIcompare ################################
#
# Functions like `plotCIs`, except plots the CIs from each stage (regular, sines
#   embedded, sines embedded & removed) for either the Bartlett- or MTM-based
#   covariance estimator (as specified by `estType` being "bart" or "mtap").
#
#     * `Nind` : the index of the realization size to be used; the row number of
#       the `bhCovCIs*` data.frame.
#
plotCIcompare <- function(resultList, type="cov", estType, Nind=1, writeImgFile=FALSE) {
  stopifnot(length(Nind)==1)
  stopifnot(Nind-floor(Nind)==0)
  
  with(resultList,
    {
      if (estType=="bart") {
        pType <- 16
        plotCol <- "red"
      } else if (estType=="mtap") {
        pType <- 17
        plotCol <- "blue"
      } else if (estType=="both") {
        pType <- 17
        plotCol <- "blue"
      } else {
        stop("Set a valid `estType`: either \"bart\" or \"mtap\".")
      }
      
      stages <- 1:3
        
      # MSE plot
      CIendpointsStr <- paste0(
        paste0(
          rep("bhCovCIs", 6), rep(c("", ".s", ".s.w"), 2), "$qSE.",
          rep(c("L", "U"), each = 3), ".", type, ".", ifelse(estType=="both", "mtap", estType), "[Nind]"
        ),
        collapse = ",\n")
      if (estType=="both") {
        CIendpointsStr <- paste0(CIendpointsStr, ",\n", gsub("\\.mtap", ".bart", CIendpointsStr))
      }
      MSE.plotLims <- eval(parse(text = paste0("range( c(",
                                               CIendpointsStr, ") )")))
      if (writeImgFile) {
        currTime <- gsub(" ", "-", gsub(":","", gsub("-","",Sys.time())))
        pdf(paste0("img/MSEbetahats_", currTime, ".pdf"), width=7, height=4)
      }
      plotText <- "plot(x=stages, y=rep(bhCovCIs$mse.cv.mtap[Nind],3), xlab=\"Realization Type\",
           ylab=\"SE of cov estimator\", xlim=c(0.5,3.5), xaxt=\"n\",
           log=\"y\", ylim=MSE.plotLims*10**c(-0.5,0.5), col=\"white\")"
      pointsText <- "points(x=stages[1], y=bhCovCIs$mse.cv.mtap[Nind], col=plotCol, pch=pType)
      points(x=stages[2], y=bhCovCIs.s$mse.cv.mtap[Nind], col=plotCol, pch=pType)
      points(x=stages[3], y=bhCovCIs.s.w$mse.cv.mtap[Nind], col=plotCol, pch=pType)"
      arrowsText <- "arrows(x0=stages[1], y0=bhCovCIs$qSE.L.cov.mtap[Nind], x1=stages[1], y1=bhCovCIs$qSE.U.cov.mtap[Nind],
             length=0.05, angle=90, code=3, lwd=2, col=plotCol)
      arrows(x0=stages[2], y0=bhCovCIs.s$qSE.L.cov.mtap[Nind], x1=stages[2], y1=bhCovCIs.s$qSE.U.cov.mtap[Nind],
             length=0.05, angle=90, code=3, lwd=2, col=plotCol)
      arrows(x0=stages[3], y0=bhCovCIs.s.w$qSE.L.cov.mtap[Nind], x1=stages[3], y1=bhCovCIs.s.w$qSE.U.cov.mtap[Nind],
             length=0.05, angle=90, code=3, lwd=2, col=plotCol)"
      legendText <- "legend(\"topright\", legend=c(\"Bartlett\",\"Multitaper\"), pch=c(16,17),
            col=c(\"red\",\"blue\"))"
      
      if (type=="cor") {
        plotText <- gsub("cov", "cor", plotText)
        plotText <- gsub("cv", "cor", plotText)
        pointsText <- gsub("cov", "cor", pointsText)
        pointsText <- gsub("cv", "cor", pointsText)
        arrowsText <- gsub("cov", "cor", arrowsText)
        arrowsText <- gsub("cv", "cor", arrowsText)
      }
      if (estType=="bart") {
        plotText <- gsub("\\.mtap", ".bart", plotText)
        pointsText <- gsub("\\.mtap", ".bart", pointsText)
        arrowsText <- gsub("\\.mtap", ".bart", arrowsText)
      } else if (estType=="both") {
        horiOffset <- 0.05
        pointsText <- gsub("x=", paste0("x=+", horiOffset, "+"), pointsText)
        arrowsText <- gsub("x0=", paste0("x0=+", horiOffset, "+"), arrowsText)
        arrowsText <- gsub("x1=", paste0("x1=+", horiOffset, "+"), arrowsText)
        pointsText.b <- gsub("\\.mtap", ".bart", pointsText)
        arrowsText.b <- gsub("\\.mtap", ".bart", arrowsText)
        pointsText.b <- gsub("x=+", "x=-", pointsText.b)
        arrowsText.b <- gsub("x0=\\+", "x0=-", arrowsText.b)
        arrowsText.b <- gsub("x1=\\+", "x1=-", arrowsText.b)
        
        pointsText.b <- gsub("plotCol", "\"red\"", pointsText.b)
        arrowsText.b <- gsub("plotCol", "\"red\"", arrowsText.b)
        pointsText.b <- gsub("pType", "16", pointsText.b)
        
        pointsText <- paste0(pointsText, "\n", pointsText.b)
        arrowsText <- paste0(arrowsText, "\n", arrowsText.b)
      }
      
      # layout(matrix(c(1,1,1,1,2,2), 3, 2, byrow = TRUE))
      par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
      suppressWarnings({
        eval(parse(text=paste(plotText, pointsText, arrowsText, sep="\n"))) # produce the plot
      })
      if (estType=="both") { eval(parse(text=legendText)) }
      axis(side=1, at=1:3, labels=c("orig.","with LCs","LCs removed")) # make x-axis
      
      if (writeImgFile) { dev.off() }
    }
  )
} # end plotCIs


################################### plotCovB ###################################
#
# Plot the covariances of the beta.hat.1 and beta.hat.2 fitted coefficients
#   from the simulations/regressions.
#
plotCovB <- function(resultList, type="cov", stage="", Nind=1) {
  stopifnot(length(Nind)==1)
  stopifnot(Nind-floor(Nind)==0)
  if (stage=="") {
    suffix <- ""
  } else if (stage=="s") {
    suffix <- ".s"
  } else if (stage=="s.w") {
    suffix <- ".s.w"
  } else {
    stop("Set a valid `stage`.")
  }
  
  plotText <- "with(resultList,
  {
    with(betasOverN[[Nind]],
      {
        plot(cv.bart, ylim=range(c(cv.bart, cv.mtap, cv.theo, bhCovCIs$smpl.cov[1])),
             col=\"red\", xlab=\"Regression #\", ylab=\"cov(b1,b2)\")
        points(x=1:length(cv.mtap), y=cv.mtap, col=\"blue\")
        abline(h=cv.theo[1], col=\"forestgreen\", lwd=2)
      })
    abline(h=bhCovCIs$smpl.cov[1], col=\"goldenrod3\", lty=2)
  })"
  
  if (type=="cor") {
    plotText <- gsub("cov", "cor", plotText)
    plotText <- gsub("cv", "cor", plotText)
  }
  plotText <- gsub("bhCovCIs", paste0("bhCovCIs", suffix), plotText)
  plotText <- gsub("betasOverN", paste0("betasOverN", suffix), plotText)
  
  par(mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
  eval(parse(text=plotText)) # produce the plot
} # end plotCovB


################################## bivARp.spec #################################
#
# Compute the spectral matrix for a bivariate VAR(p) process at a given set of
#   frequencies.
#
bivARp.spec <- function(phiMat, V=diag(1,2), freqs=seq(0,0.5,0.01)) {
  p <- dim(phiMat)[2]/2
  V.inv <- solve(V)
  specMtx <- array(0, dim=c(2, 2, length(freqs)))
  for (f.ind in 1:length(freqs)) {
    expSeq <- exp(-1i*2*pi * freqs[f.ind] * seq(1, p)) # to be used in expMat, below
    # expMat is a matrix of complex exponentials; a column of 2x2 identities scaled
    #   by exp(-i*2*pi*f*k), for k=1..p
    expMat <-
      matrix(
        matrix(c(expSeq, rep(0, 2 * p), expSeq), ncol=p, byrow=TRUE),
        nrow = 2 * p,
        ncol = 2,
        byrow = TRUE
      )
    P <- diag(1,2) - phiMat %*% expMat
    P.H <- Conj(t( diag(1,2) - phiMat %*% expMat ))
    specMtx[,,f.ind] <- solve(P.H %*% V.inv %*% P)
  }
  return(specMtx)
} # end bivARp.spec


################################### prewhiten ##################################
#
# Function to pre-whiten a univariate time series.
#
prewhiten <- function(x, sigLevel) {
  stopifnot(length(x)>0)
  N <- length(x)
  if(missing(sigLevel)) { sigLevel <- 1-1/N }
  stopifnot(abs(sigLevel-0.5) < 0.5)
  
  seas <- determineSeasonal(x, sigCutoff=sigLevel)
  
  detectedSines <- apply(seas$sinusoidData, 1, sum)
  
  x.ns <- x - detectedSines - mean(x)
  
  pilotSpec <- spec.mtm(ts(x.ns), nw=10, k=19, plot=FALSE)
  pilotSpec.full <- c(pilotSpec$spec, rev(pilotSpec$spec[-1])[-1])
  theACVF <- Re(fft(pilotSpec.full, inverse=TRUE)) / length(pilotSpec.full)
  
  ARcoefs <- acf2AR(head(theACVF,9))
  x.ns.filt <- stats::filter(x=x.ns, filter=c(0,ARcoefs[1,1]), sides=1)
  x.ns.resid <- x.ns - x.ns.filt
  x.ns.resid <- x.ns.resid[which(!is.na(x.ns.resid))]
  
  varARresid <- var(spec.mtm(ts(x.ns.resid), plot=F)$spec)
  
  p.best <- 1
  for (p in 2:8) {
    new.x.ns.filt <- stats::filter(x=x.ns, filter=c(0,ARcoefs[p,(1:p)]), sides=1)
    new.x.ns.resid <- x.ns - new.x.ns.filt
    new.x.ns.resid <- new.x.ns.resid[which(!is.na(new.x.ns.resid))]
    
    newVar <- var(spec.mtm(ts(new.x.ns.resid), plot=F)$spec)
    
    if (newVar < varARresid) {
      p.best <- p
      x.ns.resid <- new.x.ns.resid
      varARresid <- newVar
    }
  }
  
  x.pw <- c(x.ns.resid, rep(0,p.best)) + detectedSines
  
  return(x.pw)
} # end prewhiten