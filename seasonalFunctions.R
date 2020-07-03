############################# begin Dave R comment #############################
#
# Returns information about the seasonal components in your data
# Depends on sigCutoff (a value between 0 and 1, probably want it around 0.99 ?)
# Set deltat appropriately to get frequency in "cycles per unit"
## The cow example was "cycles per month", deltat was 1 in that case.
#
# Example - if your data is sampled every week, you could use
## deltat = 1, and get "cycles per week"
## OR
## deltat = 7 and get "cycles per day" (I think...) Compare with multitaper to
## make sure
#
############################## end Dave R comment ##############################
#
# Edited to include the padFactor parameter for zero-padding.
#
################################################################################
determineSeasonal.old <- function(data, sigCutoff, padFactor=7, deltat=1, predictNum = 0){
  require('multitaper')
  
  N <- length(data)
  nFFT <- 2^(floor(log2(N)) + padFactor)
  
  spec <- spec.mtm(data, deltat=deltat, Ftest=TRUE, returnInternals=T, plot=F,
                   nFFT=nFFT)
  cmv <- spec$mtm$cmv
  
  freqsIdx <- findLocalFMax(spec, sigCutoff)
  numFreqs <- length(freqsIdx)
  
  if (numFreqs==0) {
    warning("No sinusoids detected.")
    return(NULL)
  }
  
  sines <- matrix(NA, nrow=N, ncol=length(freqsIdx))
  
  phseAmp <- as.data.frame(matrix(NA, numFreqs, 3)) # dataframe for params
  names(phseAmp) <- c("freq", "amp", "phase")
  
  for (i in 1:length(freqsIdx)){
    f.index <- freqsIdx[i]
    # blank <- matrix(data=0,nrow=nFFT,ncol=1)
    # blank[f.index] <- cmv[f.index]
    # blank[nFFT-f.index+2] <- Conj(cmv[f.index])
    # 
    # inv <- fft(blank,inverse=TRUE)
    # sines[, i] <- Re(inv[1:N])
    phseAmp[i, (1:3)] <- c(spec$freq[f.index],
                           2 * Mod(cmv[f.index]),
                           atan2(Im(cmv[f.index]), Re(cmv[f.index])))
    
    sines[,i] <- with(phseAmp,
      { amp[i] * cos(2*pi*freq[i]*seq(0,N-1,deltat) + phase[i]) }
    )
  }
  
  if (predictNum > 0){
    f.index <- freqsIdx
    
    blank <- matrix(data=0,nrow=nFFT,ncol=1)
    blank[f.index] <- cmv[f.index]
    blank[nFFT-f.index+2] <- Conj(cmv[f.index])
    
    inv <- fft(blank,inverse=TRUE)
    
    predicted <- Re(inv[(N+1):(N+predictNum)])
  } else {
    predicted <- NA
  }
  
  list(sinusoidData=sines,
       phaseAmplitudeInfo=phseAmp,
       predicted=predicted)
}

# checks for whether the significant point in the F-test is 
# larger than the previous and next point.  If yes, we have 
# a local max.
# - returns these frequencies
# obj needs to be either a) a spec.mtm object (Ftest class) or 
# b) a driegert.cmv object
findLocalFMax <- function(obj, cutoff){
  # Check whether this is a spec.mtm() object, or from my own CMV code.
  if (any(class(obj) == "Ftest")){
    Fval <- obj$mtm$Ftest
    k <- obj$mtm$k
  } else if (any(class(obj) == "driegert.cmv")){
    Fval <- obj$Ftest$Fval
    k <- obj$k
  } else {
    stop("obj needs to be of class 'driegert.cmv' or 'Ftest'.")
  }
  
  fMaxInd <- which(Fval > qf(cutoff, 2, 2*k-2))
  maxes <- c()
  
  if (length(fMaxInd) == 0){
    return(maxes)
  }
  
  for (i in 1:length(fMaxInd)){
    if (fMaxInd[i] == 1 || fMaxInd[i] == length(Fval)){
      next
    }
    
    if (Fval[fMaxInd[i]] > Fval[fMaxInd[i]-1] && 
        Fval[fMaxInd[i]] > Fval[fMaxInd[i]+1]){
      maxes <- c(maxes, fMaxInd[i])
    }
  }
  
  maxes
}

## Determines the phase and amplitude of a sinusoid 
## fit to data (ie. the data from removePeriod() )
fitSinusoidSingle <- function(series1, dT, f){
  t <- seq(0, by=dT, length.out=length(series1))
  fit <- lm(series1 ~ sin(2*pi*f*t) + cos(2*pi*f*t) -1)
  phse1 <- atan(fit$coef[2] / fit$coef[1])
  amp1 <- fit$coef[1] / cos(phse1)
  
  c(amp1, phse1)
}


################################################################################
# determineSeasonal, but it jackknifes the line components
################################################################################
determineSeasonal <- function(data, sigCutoff, padFactor=7, deltat=1,
                                 predictNum=0, jackknife=FALSE) {
  require('multitaper')
  
  N <- length(data)
  nFFT <- 2^(floor(log2(N)) + padFactor)
  
  # first, find the line components based on using all the eigencoefs
  spec <- spec.mtm(data, deltat=deltat, Ftest=TRUE, returnInternals=T, plot=F,
                   nFFT=nFFT)
  cmv <- spec$mtm$cmv
  
  freqsIdx <- findLocalFMax(spec, sigCutoff)
  numFreqs <- length(freqsIdx)
  
  if (numFreqs==0) {
    warning("No sinusoids detected.")
    return(NULL)
  }
  
  sines <- matrix(NA, nrow=N, ncol=length(freqsIdx))
  
  phseAmp <- as.data.frame(matrix(NA, numFreqs, 3)) # dataframe for params
  names(phseAmp) <- c("freq", "amp", "phase")
  
  for (i in 1:length(freqsIdx)){
    f.index <- freqsIdx[i]
    # blank <- matrix(data=0,nrow=nFFT,ncol=1)
    # blank[f.index] <- cmv[f.index]
    # blank[nFFT-f.index+2] <- Conj(cmv[f.index])
    # 
    # inv <- fft(blank,inverse=TRUE)
    # sines[, i] <- Re(inv[1:N])
    phseAmp[i, (1:3)] <- c(spec$freq[f.index],
                           2 * Mod(cmv[f.index]),
                           atan2(Im(cmv[f.index]), Re(cmv[f.index])))
    
    sines[,i] <- with(phseAmp,
      { amp[i] * cos(2*pi*freq[i]*seq(0,N-1,deltat) + phase[i]) }
    )
  }
  
  if (jackknife) {
    spec.jk <- spec.mtm(data, deltat=deltat, Ftest=TRUE, returnInternals=T, plot=F,
                        nFFT=nFFT)
    K <- spec.jk$mtm$k
    sines.jk <- array(NA, dim=c(dim(sines),K))
    phseAmp.jk <- array(NA, dim=c(numFreqs, 3, K)) # array for params
    
    # eigencoefficients and weights
    # y <- spec.jk$mtm$eigenCoefs
    # d <- ifelse(is.null(spec.jk$mtm$eigenCoefWt),
    #         matrix(1, nrow=dim(y)[1], ncol=dim(y)[2]), spec.jk$mtm$eigenCoefWt)
    dw <- spec.jk$mtm$dpss$v
    swz <- apply(dw, 2, sum)
    swz[seq(2,K,2)] <- 0
    # ssqswz <- as.numeric(t(swz)%*%swz)
    # cft <- mvfft(rbind(dw * data, matrix(0, nFFT-N, K)))
    # cft <- cft[seq(1, nFFT/2+1),]
    # cmv <- (cft %*% swz) / ssqswz
    
    for (k in 1:K) {
      # spec.k <- apply(Mod(d[,-k])**2 * Mod(y[,-k])**2, MARGIN=1, FUN=sum) /
      #   apply(Mod(d[,-k])**2, MARGIN=1, FUN=sum)
      
      cmv.k <- (spec.jk$mtm$eigenCoefs[,-k] %*% swz[-k]) / as.numeric(t(swz[-k])%*%swz[-k])
      
      for (i in 1:length(freqsIdx)) {
        f.index <- freqsIdx[i]
        phseAmp.jk[i, (1:3), k] <- c(spec$freq[f.index],
                                     2 * Mod(cmv.k[f.index]),
                                     atan2(Im(cmv.k[f.index]), Re(cmv.k[f.index])))
        
        sines.jk[,i,k] <- phseAmp.jk[i,2,k] *
          cos(2*pi*phseAmp.jk[i,1,k]*seq(0,N-1,deltat) + phseAmp.jk[i,3,k])
      }
    }
  } else {
    phseAmp.jk <- NA
    sines.jk <- NA
  }
  
  if (predictNum > 0){
    f.index <- freqsIdx
    
    blank <- matrix(data=0,nrow=nFFT,ncol=1)
    blank[f.index] <- cmv[f.index]
    blank[nFFT-f.index+2] <- Conj(cmv[f.index])
    
    inv <- fft(blank,inverse=TRUE)
    
    predicted <- Re(inv[(N+1):(N+predictNum)])
  } else {
    predicted <- NA
  }
  
  list(sinusoidData=sines,
       phaseAmplitudeInfo=phseAmp,
       predicted=predicted,
       phaseAmplitudeInfo.jk=phseAmp.jk, sinusoidData.jk=sines.jk)
}
