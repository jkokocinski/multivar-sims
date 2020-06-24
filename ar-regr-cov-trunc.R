################################################################################
ar.regr.cov.trunc <- function(phiMat.p, phiMat.r, errCovMat.p, errCovMat.r,
                        numObsVec, NUM_REGR,
                        mtmFixed="NW", timeBandProd=4, numTapers=7, W,
                        adaptWt=FALSE,
                        embedSines=TRUE,
                        linDepY=FALSE, removeLCs=FALSE,
                        truncs=seq(1,numObsVec-1)) {
  
  library(mAr) # depends on MASS
  library(multitaper)
  
  params.init.str <- ls() # vector of initally set parameters
  params.init <- eval(parse(text=paste0("list(",
    paste0(params.init.str,"=",params.init.str, collapse=",\n"),")")))
  
  if(length(intersect(c("mAr", "multitaper"), rownames(installed.packages())))<2) {
  stop("Ensure `mAr` and `multitaper` packages are loaded.")
  }
  
  stopifnot(length(numObsVec)==1)
  
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
  betas.head <- c("b1","b2")
  
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
      a <- 1 #2.5 # scale on Y_2
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
    covB.mtap.mat <- covB.bart.mat <- array(dim=c(NUM_REGR, length(truncs)))
    var.b1.mtap.mat <- var.b1.bart.mat <- array(dim=c(NUM_REGR, length(truncs)))
    var.b2.mtap.mat <- var.b2.bart.mat <- array(dim=c(NUM_REGR, length(truncs)))
    
    # predictor series (X_1,X_2) realization
    bivAR1.p <- as.matrix(mAr.sim(w=rep(0,2), A=phiMat.p, C=errCovMat.p, N=numObs))
    
    freqsToEmbed <- (c(90,90))**(-1)
    
    if (embedSines) {
      bivAR1.p.s <- embedSinusoids(input=bivAR1.p, freqs=freqsToEmbed,
                                   amps=sqrt(diag(errCovMat.p)), ampScale=8)
      # bivAR1.p.s <- embedSinusoids(input=bivAR1.p.s, freqs=c(90,90)**(-1),
      #                              amps=sqrt(diag(errCovMat.r)), ampScale=0.7)
      # bivAR1.p.s <- embedSinusoids(input=bivAR1.p.s, freqs=c(45,45)**(-1),
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
      }
      
      respRlzns[,,j] <- bivAR1.r
      
      if (embedSines) {
        bivAR1.r.s <- embedSinusoids(input=bivAR1.r, freqs=freqsToEmbed,
                                     amps=sqrt(diag(errCovMat.r)), ampScale=10)
        # bivAR1.r.s <- embedSinusoids(input=bivAR1.r.s, freqs=c(90,90)**(-1),
        #                              amps=diag(errCovMat.r), ampScale=0.25)
        # bivAR1.r.s <- embedSinusoids(input=bivAR1.r.s, freqs=c(45,45)**(-1),
        #                              amps=diag(errCovMat.r), ampScale=0.1)
        respRlzns.s[,,j] <- bivAR1.r.s
        
        # detect & remove sinusoidal line components
        if (removeLCs) {
          for (p in 1:2) {
            commonSinesObj <- findCommonSines(x=bivAR1.p.s[,p], y=bivAR1.r.s[,p],
              freqThresh=ifelse(mtmFixed=="NW", timeBandProd / numObs, W),
              sigCutoff=0.999)
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
      
      # \cov(\hat{\beta}_1, \hat{\beta}_2) -- theoretical-based
      covB.theo <- A1 %*% toepLeftMult2( theo.ccv.r, as.vector(t(A2)) )
      
      # theoretical correlations
      var.b1.theo <- A1 %*% toepLeftMult2( theo.acv1.r, as.vector(t(A1)) )
      var.b2.theo <- A2 %*% toepLeftMult2( theo.acv2.r, as.vector(t(A2)) )
      corB.theo <- covB.theo / sqrt(var.b1.theo * var.b2.theo)
      
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
        
        # Bartlett CCVF of response series (full)
        bart.ccv.r <- ccf(x=Y[,1], y=Y[,2], type="covariance",
                          lag.max=numObs-1, plot=FALSE)$acf
        
        # \cov(\hat{\beta}_1, \hat{\beta}_2) -- Bartlett-based
        covB.barts <- rep(NA, length(truncs)) # initialize the vector
        for (tr.i in 1:length(truncs)) {
          covB.barts[tr.i] <- A1 %*%
            toepLeftMult2( as.vector(bart.ccv.r) *
                             (abs(numObs-(1:(2*numObs-1))) < truncs[tr.i]+1),
                           as.vector(t(A2)) )
        }
        covB.bart.mat[j,] <- covB.barts # update matrix
        
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
        covB.mtaps <- rep(NA, length(truncs)) # initialize the vector
        for (tr.i in 1:length(truncs)) {
          covB.mtaps[tr.i] <- A1 %*%
            toepLeftMult2( mtap.ccv.r *
                             (abs(numObs-(1:(2*numObs-1))) < truncs[tr.i]+1),
                           as.vector(t(A2)) )
        }
        covB.mtap.mat[j,] <- covB.mtaps # update matrix
        
        if (TRUE) {
          # Bartlett ACVFs for responses; for use in calculating corr
          bart.acv1.r <- acf(x=Y[,1], type="covariance", lag.max=numObs-1,
                             plot=FALSE)
          bart.acv2.r <- acf(x=Y[,2], type="covariance", lag.max=numObs-1,
                             plot=FALSE)
          bart.acv1.r <- c(rev(bart.acv1.r$acf), bart.acv1.r$acf[-1])
          bart.acv2.r <- c(rev(bart.acv2.r$acf), bart.acv2.r$acf[-1])
          
          var.b1.barts <- rep(NA, length(truncs)) # initialize the vector
          var.b2.barts <- rep(NA, length(truncs)) # initialize the vector
          for (tr.i in 1:length(truncs)) {
            var.b1.barts[tr.i] <- A1 %*%
              toepLeftMult2( bart.acv1.r *
                               (abs(numObs-(1:(2*numObs-1))) < truncs[tr.i]+1),
                             as.vector(t(A1)) )
            var.b2.barts[tr.i] <- A2 %*%
              toepLeftMult2( bart.acv2.r *
                               (abs(numObs-(1:(2*numObs-1))) < truncs[tr.i]+1),
                             as.vector(t(A2)) )
          }
          var.b1.bart.mat[j,] <- var.b1.barts # update matrix
          var.b2.bart.mat[j,] <- var.b2.barts # update matrix
          
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
          
          var.b1.mtaps <- rep(NA, length(truncs)) # initialize the vector
          var.b2.mtaps <- rep(NA, length(truncs)) # initialize the vector
          for (tr.i in 1:length(truncs)) {
            var.b1.mtaps[tr.i] <- A1 %*%
              toepLeftMult2( mtap.acv1.r *
                               (abs(numObs-(1:(2*numObs-1))) < truncs[tr.i]+1),
                             as.vector(t(A1)) )
            var.b2.mtaps[tr.i] <- A2 %*%
              toepLeftMult2( mtap.acv2.r *
                               (abs(numObs-(1:(2*numObs-1))) < truncs[tr.i]+1),
                             as.vector(t(A2)) )
          }
          var.b1.mtap.mat[j,] <- var.b1.mtaps # update matrix
          var.b2.mtap.mat[j,] <- var.b2.mtaps # update matrix
          # corB.mtap <- covB.mtap / sqrt(var.b1.mtap * var.b2.mtap)
        }
        
        # update the betas data.frame
        betas.newRow <- c(as.numeric(model1$coefficients),
                          as.numeric(model2$coefficients))
        betas[j,] <- betas.newRow
        
        setTxtProgressBar(pb=pb2, value=(tp-1)*NUM_REGR + j)
      } # end for (j in 1:NUM_REGR)
      
      # make the corr matrices over trunc length
      corB.bart.mat <- covB.bart.mat / sqrt(var.b1.bart.mat * var.b2.bart.mat)
      corB.mtap.mat <- covB.mtap.mat / sqrt(var.b1.mtap.mat * var.b2.mtap.mat)
      
      setTxtProgressBar(pb=pb2, value=0)
    
      # update the list of dfs
      if (tp==1) {
        betasOverN[[n]] <- betas
      } else if (tp==2) {
        betasOverN.s[[n]] <- betas
      } else if (tp==3) {
        betasOverN.s.w[[n]] <- betas
      } else {
        stop("impossible tp")
      }
      
    } # end for (tp)
    cat(paste0("# ", Sys.time()) ," | done!\n")
  } # end for (n in seq(1,length(numObsVec)))
  
  
  rm(betas,betas.head, betas.newRow)
  
  localVars <- setdiff(ls(), names(params.init))
  
  returnExpr <- paste0("params.init=params.init,\n",
                       paste0(localVars,"=",localVars, collapse=",\n")
  )
  returnExpr <- paste0("return( list(",returnExpr, ") )")
  # cat(paste0(returnExpr,"\n"))
  
  # return statement
  eval(parse(text=returnExpr))
  
}
