# Returns information about the seasonal components in your data
# Depends on sigCutoff (a value between 0 and 1, probably want it around 0.99 ?)
# Set deltat appropriately to get frequency in "cycles per unit"
## The cow example was "cycles per month", deltat was 1 in that case.

# Example - if your data is sampled every week, you could use
## deltat = 1, and get "cycles per week"
## OR
## deltat = 7 and get "cycles per day" (I think... ) Compare with multitaper to make sure
determineSeasonal <- function(data, sigCutoff, deltat=1, predictNum = 0){
  require('multitaper')
  
  N <- length(data)
  nFFT <- 2^(floor(log2(N)) + 7)
  
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
  # sines2 <- sines
  
  phseAmp <- as.data.frame(matrix(NA, numFreqs, 3))
  names(phseAmp) <- c("freq", "amp", "phase")
  # phseAmp$amp2 <- NA; phseAmp$phase2 <- NA
  
  for (i in 1:length(freqsIdx)){
    f.index <- freqsIdx[i]
    blank <- matrix(data=0,nrow=nFFT,ncol=1)
    blank[f.index] <- cmv[f.index]
    blank[nFFT-f.index+2] <- Conj(cmv[f.index])
    
    inv <- fft(blank,inverse=TRUE)
    sines[, i] <- Re(inv[1:N])
    # phseAmp[i,(1:3)] <- c(spec$freq[f.index],
    #                       fitSinusoidSingle(sines[,i], 1, f=spec$freq[f.index]))
    # phseAmp[i, (4:5)] <- c( 2 * Mod(cmv[f.index]),
    #                        atan2(Im(cmv[f.index]), Re(cmv[f.index])) )
    
    # sines2 [,i] <- with(phseAmp,
    #   { amp2[i] * cos(2*pi*freq[i]*seq(1,N,deltat) + phase2[i]) }
    # )
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
       sinusoidData2=sines2,
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
  
  fMaxInd <- which(Fval > qf(cutoff, 2, 2*k))
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
