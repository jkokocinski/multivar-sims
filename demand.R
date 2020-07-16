#!/usr/bin/env Rscript

tryCatch(setwd("~/multivar-sims"),
         error=function(e) {
           tryCatch( setwd("~/WORK/Q5/multivar-sims") )
         }
)

library(mAr)
library(multitaper)
library(signal)
library(plot.matrix)
library(DescTools)
# library(parallel); library(foreach); library(doParallel)
source("helpers_mvSimRegr.R")
source("seasonalFunctions.R")

############################### parallel options ###############################
# numCores <- ceiling(parallel::detectCores()/4)
# cl <- parallel::makeCluster(numCores)
# registerDoParallel(cl)

############################ import data and format ############################
importData <- function(prefix, seqInds, suffix="", filext, skip=0) {
  theData <- data.frame()
  
  for (indind in 1:length(seqInds)) {
  	newData <- read.csv(file=paste0(prefix, seqInds[indind], suffix, filext),
  	                    header=TRUE, stringsAsFactors=FALSE, skip=skip)
  	theData <- rbind(theData, newData)
  }
  return(theData)
}

demand <- importData("data/demand/PUB_Demand_", 2003:2020, "", ".csv", skip=3)
hoep <- importData("data/hoep/PUB_PriceHOEPPredispOR_", 2003:2020, "", ".csv", skip=3)

# temperat <- data.frame()
# for (yr in 2003:2012) {
#   temperat.temp <- importData("data/en_climate_hourly_ON_6158733_", sprintf("%02.0f", 1:12), paste0("-", yr, "_P1H"), ".csv", skip=0)
#   temperat <- rbind(temperat, temperat.temp)
#   rm(temperat.temp)
# }

demand$datetime <- as.POSIXct(
  x=paste0(demand$Date, " ", sprintf("%02.0f", demand$Hour), ":00:00"),
  format="%Y-%m-%d %H:%M:%S", tz="UTC"
)

hoep$datetime <- as.POSIXct(
  x=paste0(hoep$Date, " ", sprintf("%02.0f", hoep$Hour), ":00:00"),
  format="%Y-%m-%d %H:%M:%S", tz="UTC"
)

iesoData <- merge.data.frame(x=demand, y=hoep, sort=TRUE)

iesoData$year <- as.integer(substr(iesoData$Date, 1, 4))
iesoData$month <- as.integer(substr(iesoData$Date, 6, 7))
iesoData$day <- as.integer(substr(iesoData$Date, 9, 10))

iesoData <- iesoData[,c("year","month","day","Hour","Ontario.Demand","HOEP",
                        "datetime")]
names(iesoData)[4:6] <- c("hour", "demand", "hoep")
remove(demand, hoep) # clean up

iesoData$demand <- as.numeric(iesoData$demand)
iesoData$hoep <- gsub(",", "", iesoData$hoep) # remove thousands separator
iesoData$hoep <- as.numeric(iesoData$hoep)

iesoData <- iesoData[ with(iesoData, { order(year, month, day, hour) }), ]
row.names(iesoData) <- as.character(seq(1,dim(iesoData)[1])) # correct row.names

# # temperature data
# temperat <- temperat[,(6:10)]
# temperat[,4] <- as.numeric(substr(x=temperat[,4], start=1, stop=2)) # convert time to integer
# names(temperat) <- c("year","month","day","hour","temp")
# 
# # weather data from SSC case study
# weather <- readxl::read_xlsx(path="data/ssc/ssc2020_hourly_weather.xlsx", sheet=2,
#                              col_names=TRUE, trim_ws=TRUE, progress=TRUE)
# weather <- as.data.frame(weather)
# weather$year <- as.integer(substr(weather$time, 1, 4))
# weather$month <- as.integer(substr(weather$time, 6, 7))
# weather$day <- as.integer(substr(weather$time, 9, 10))
# weather$hour <- as.integer(substr(weather$time, 12, 13))+1


################################# plot the data ################################
# png(file="img/demand-hoep_plots.png", width=640, height=720)
# par(mfrow=c(2,1))
# par(mar=c(4,4,1,1))
# plot(x=iesoData$datetime, y=iesoData$demand, type="l",
#      xlab="Date", ylab="Ontario Demand (kW)")
# plot(x=iesoData$datetime, y=iesoData$hoep, type="l",
#      xlab="Date", ylab="HOEP")
# dev.off()

######################## compute auto- and cross-spectra #######################
# png(file="img/3_spec-demand-full.png", width=640, height=320)
# par(mar=c(5,4,1,1))
# multitaper::spec.mtm(ts(iesoData$demand), nw=4, k=7, main="")
# dev.off()
# png(file="img/3_spec-hoep-full.png", width=640, height=320)
# par(mar=c(5,4,1,1))
# multitaper::spec.mtm(ts(iesoData$hoep), nw=4, k=7, main="")
# dev.off()


############################## do the regressions ##############################
hrNum <- 1:(365*3*24) # hour number vector to be passed to lm as the predictor variable; for removing trend and de-meaning

timeBandProd <- 4
numTapers <- 7
adaptWt <- TRUE # adaptive weigting in mtm cross-spec estimate?

timeSegs <- 2003:2016
nts <- length(timeSegs)
covB.mtap.mat <- covB.bart.mat <- matrix(NA, nts, nts)
corB.mtap.mat <- corB.bart.mat <- matrix(NA, nts, nts)
betaHat0Vec <- rep(NA, nts)
betaHatLCs <- matrix(NA, nrow=nts, ncol=3) # 3 columns for persistent LCs (daily, semidaily, octaweekly)
jkLC <- F # jackknife beta-hats for line components?
betaHatLCs.jk <- array(NA, dim=c(dim(betaHatLCs), numTapers))

for (j1 in 2:(nts-1)) {
# foreach(j1 = 2:(5)) %:%
  # foreach(j2 = j1:5) %dopar% {
  for (j2 in j1:(nts-1)) {
    inds1 <- head(which(iesoData$year %in% timeSegs[j1+(-1:1)]), length(hrNum))
    inds2 <- head(which(iesoData$year %in% timeSegs[j2+(-1:1)]), length(hrNum))
    Y1 <- iesoData[inds1, "hoep"]
    Y2 <- iesoData[inds2, "hoep"]
    X1 <- iesoData[inds1, "demand"]
    X2 <- iesoData[inds2, "demand"]
    # X1 <- (weather[which(weather$year %in% timeSegs[j1+(-1:1)]), "temperature"])[hrNum]
    # X2 <- (weather[which(weather$year %in% timeSegs[j2+(-1:1)]), "temperature"])[hrNum]
    
    numObs <- length(Y1)
    
    # demean the time series
    Y1 <- Y1 - mean(Y1)
    Y2 <- Y2 - mean(Y2)
    X1 <- X1 - mean(X1)
    X2 <- X2 - mean(X2)
    
    # prewhitened series
    X1.pw <- prewhiten(X1, sigLevel=0.999)
    Y1.pw <- prewhiten(Y1, sigLevel=0.999)
    X2.pw <- prewhiten(X2, sigLevel=0.999)
    Y2.pw <- prewhiten(Y2, sigLevel=0.999)
    
    # find common sinusoidal components
    cs1Obj <- findCommonSines(x=X1.pw, y=Y1.pw, padFactor=3, freqThresh=timeBandProd/numObs, sigCutoff=0.999,
                              NW=timeBandProd, K=numTapers, jackknife=jkLC)
    if (j1==j2) {
      cs2Obj <- cs1Obj
    } else {
      cs2Obj <- findCommonSines(x=X2.pw, y=Y2.pw, padFactor=3, freqThresh=timeBandProd/numObs, sigCutoff=0.999,
                                NW=timeBandProd, K=numTapers, jackknife=jkLC)
    }
    
    # indices in the params data.frame of the persistent LCs
    LCs.X1.ind <- c(which(abs(cs1Obj$paramsX.com$freq-1/24) < 1/numObs),
                    which(abs(cs1Obj$paramsX.com$freq-1/12) < 1/numObs),
                    which(abs(cs1Obj$paramsX.com$freq-1/21) < 1/numObs))
    LCs.Y1.ind <- c(which(abs(cs1Obj$paramsY.com$freq-1/24) < 1/numObs),
                    which(abs(cs1Obj$paramsY.com$freq-1/12) < 1/numObs),
                    which(abs(cs1Obj$paramsY.com$freq-1/21) < 1/numObs))
    LCs.X2.ind <- c(which(abs(cs2Obj$paramsX.com$freq-1/24) < 1/numObs),
                    which(abs(cs2Obj$paramsX.com$freq-1/12) < 1/numObs),
                    which(abs(cs2Obj$paramsX.com$freq-1/21) < 1/numObs))
    LCs.Y2.ind <- c(which(abs(cs2Obj$paramsY.com$freq-1/24) < 1/numObs),
                    which(abs(cs2Obj$paramsY.com$freq-1/12) < 1/numObs),
                    which(abs(cs2Obj$paramsY.com$freq-1/21) < 1/numObs))
    
    # beta-hat for each persistent LC
    if (j1==j2) {
      for (l in 1:3) {
        betaHatLCs[j1,l] <- cs1Obj$paramsY.com$amp[LCs.Y1.ind[l]] / cs1Obj$paramsX.com$amp[LCs.X1.ind[l]]
        if (jkLC) {
          betaHatLCs.jk[j1,l,] <- cs1Obj$paramsY.com.jk[LCs.Y1.ind[l],2,] / cs1Obj$paramsX.com$amp[LCs.X1.ind[l]]
        }
      }
    }
    
    # remove UNcommon LCs in pred. and resp.; additionally remove common LCs in pred.
    X1.w <- X1 - 0*cs1Obj$fctVals.incoh[,1] - apply(cs1Obj$fctVals.com.x, 1, sum)
    Y1.w <- Y1 - 0*cs1Obj$fctVals.incoh[,2] - apply(cs1Obj$fctVals.com.y, 1, sum)
    X2.w <- X2 - 0*cs2Obj$fctVals.incoh[,1] - apply(cs2Obj$fctVals.com.x, 1, sum)
    Y2.w <- Y2 - 0*cs2Obj$fctVals.incoh[,2] - apply(cs2Obj$fctVals.com.y, 1, sum)
    
    # phase-aligned common LCs in pred. series
    cs1.com.x.ph <- cs1Obj$fctVals.com.x.ph
    cs2.com.x.ph <- cs2Obj$fctVals.com.x.ph
    
    # transform the residuals of the HOEP series
    #   (optimal lambda for entire span of series found to be 0.1)
    Y1.w <- DescTools::BoxCox(Y1.w + abs(min(Y1.w)) + 1, lambda=0.1)
    Y2.w <- DescTools::BoxCox(Y2.w + abs(min(Y2.w)) + 1, lambda=0.1)
    
    mpi.X1.w <- MASS::ginv(X1.w)
    mpi.X2.w <- MASS::ginv(X2.w)
    # mpi.X1.w.ph <- MASS::ginv(cbind(X1.w, cs1.com.x.ph)) # for beta1 vector
    # mpi.X2.w.ph <- MASS::ginv(cbind(X2.w, cs2.com.x.ph)) # for beta2 vector
    
    # store beta-hat (new)
    if (j1==j2) {
      betaHat0Vec[j1] <- mpi.X1.w %*% Y1.w
    }
    
    # compute spec.mtm objects for both response series components
    spec.y1 <- multitaper::spec.mtm(
      timeSeries=ts(Y1.w), nw=timeBandProd, k=numTapers, adaptiveWeighting=TRUE,
      returnInternals=TRUE, plot=FALSE
    )
    spec.y2 <- multitaper::spec.mtm(
      timeSeries=ts(Y2.w), nw=timeBandProd, k=numTapers, adaptiveWeighting=TRUE,
      returnInternals=TRUE, plot=FALSE
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
    
    # F-statistics for demand and HOEP
    # par(mar=c(4,4,1,1))
    # plot(x=spec.y1$freq, y=spec.y1$mtm$Ftest, type="l") # Y1
    # par(mar=c(4,4,1,1))
    # plot(x=spec.y2$freq, y=spec.y2$mtm$Ftest, type="l") # Y2
    
    
    # plot amplitude and phase cross-spectra
    # png(file="img/demand-hoep_cross-spec.png", width=1024, height=576)
    # par(mfrow=c(2,1))
    # plotInds <- seq(1,length(crossSpecEstY)/2) # for plotting on [0,0.5)
    # freqGrid <- (plotInds-1)/length(crossSpecEstY)
    # # amplitude cross spectrum estimate
    # par(mar=c(4,4,1,1))
    # plot(x=freqGrid, y=Mod(crossSpecEstY[plotInds]), type="l", log="y",
    #      xlab="Frequency", ylab="Amplitude Cross Spectrum Estimate")
    # # abline(v=seq(0,0.5, by=1/24), col="blue", lty=2, lwd=1) # daily, bi-daily, etc.
    # # abline(v=seq(0,0.5, by=1/(24*7)), col="orange", lty=2) # weekly, bi-weely, etc.
    # # phase cross spectrum estimate
    # par(mar=c(4,4,1,1))
    # plot(x=freqGrid,
    #      y=signal::unwrap(atan2(Im(crossSpecEstY),Re(crossSpecEstY))[plotInds]),
    #      type="l", xlab="Frequency", ylab="Phase Cross Spectrum Estimate")
    # dev.off()
    
    # CCVF of Y_1 and Y_2
    mtap.ccv.r <- fft(z=crossSpecEstY, inverse=TRUE) / length(crossSpecEstY)
    # CCVF of Y_2 and Y_1
    mtap.ccv.r.21 <- fft(z=Conj(crossSpecEstY), inverse=TRUE) / length(crossSpecEstY)
    
    # put entries in "correct" order, from -numObs to +numObs
    mtap.ccv.r <- Re(c(tail(mtap.ccv.r, numObs-1), head(mtap.ccv.r, numObs)))
    mtap.ccv.r.21 <- Re(c(tail(mtap.ccv.r.21, numObs-1), head(mtap.ccv.r.21, numObs)))
    
    # \cov(\hat{\beta}_1, \hat{\beta}_2) -- multitaper-based
    covB.mtap <- mpi.X1.w %*% toepLeftMult2( mtap.ccv.r, as.vector(t(mpi.X2.w)) )
    # covB.mtap.21 <- mpi.X2.w %*% toepLeftMult2( mtap.ccv.r.21, as.vector(t(mpi.X1.w)) )
    
    # autospectra for Y1.w and Y2.w
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
    
    var.b1.mtap <- bilinToep(mpi.X1.w, mtap.acv1.r, mpi.X1.w)
    var.b2.mtap <- bilinToep(mpi.X2.w, mtap.acv2.r, mpi.X2.w)
    corB.mtap <- covB.mtap / sqrt(var.b1.mtap * var.b2.mtap)
    
    # update covB MTM matrices
    covB.mtap.mat[j1,j2] <- covB.mtap
    covB.mtap.mat[j2,j1] <- covB.mtap # covB.mtap.21
    corB.mtap.mat[j1,j2] <- corB.mtap
    corB.mtap.mat[j2,j1] <- corB.mtap
    
    bart.ccv.r <- ccf(x=Y1.w, y=Y2.w, type="covariance", lag.max=numObs-1, plot=F)
    covB.bart <- bilinToep(mpi.X1.w, as.vector(bart.ccv.r$acf), mpi.X2.w)
    
    # Bartlett ACVFs for responses; for use in calculating corr
    bart.acv1.r <- acf(x=Y1.w, type="covariance", lag.max=numObs-1, plot=FALSE)
    bart.acv2.r <- acf(x=Y2.w, type="covariance", lag.max=numObs-1, plot=FALSE)
    bart.acv1.r <- c(rev(bart.acv1.r$acf), bart.acv1.r$acf[-1])
    bart.acv2.r <- c(rev(bart.acv2.r$acf), bart.acv2.r$acf[-1])
    var.b1.bart <- bilinToep(mpi.X1.w,  bart.acv1.r, mpi.X1.w)
    var.b2.bart <- bilinToep(mpi.X2.w, bart.acv2.r, mpi.X2.w)
    corB.bart <- covB.bart / sqrt(var.b1.bart * var.b2.bart)
    
    # update covB Bartlett matrices
    covB.bart.mat[j1,j2] <- covB.bart
    covB.bart.mat[j2,j1] <- covB.bart # covB.bart.21
    corB.bart.mat[j1,j2] <- corB.bart
    corB.bart.mat[j2,j1] <- corB.bart
    
    # c(j1,j2,covB.bart,corB.bart,covB.mtap,corB.mtap)
    
    cat(paste0("# ", Sys.time()), "\tdone ", j1, ",", j2, "\n")
  }
}

# load(file="demand-result.RData")

numRows <- dim(covB.bart.mat)[1] - 2L
betaHat0Vec <- betaHat0Vec[which(!is.na(betaHat0Vec))]
betaHatLCs <- betaHatLCs[-c(1,dim(betaHatLCs)[1]),]
betaHatLCs.jk <- betaHatLCs.jk[-c(1,dim(betaHatLCs.jk)[1]),,]
BCov <- matrix(covB.bart.mat[which(!is.na(covB.bart.mat))], numRows, numRows)
MCov <- matrix(covB.mtap.mat[which(!is.na(covB.mtap.mat))], numRows, numRows)
BCor <- solve(sqrt(diag(diag(BCov)))) %*% BCov %*% solve(sqrt(diag(diag(BCov))))
MCor <- solve(sqrt(diag(diag(MCov)))) %*% MCov %*% solve(sqrt(diag(diag(MCov))))

diag(BCor) <- rep(1,numRows)
diag(MCor) <- rep(1,numRows)

save(betaHat0Vec, BCov, MCov, BCor, MCor, betaHatLCs, betaHatLCs.jk, file="demand-result.RData")

# stuff for printing in LaTeX
xtable::xtableMatharray(x=BCor, digits=2, display=rep("G", dim(BCov)[2]+1))
xtable::xtableMatharray(x=MCor, digits=2, display=rep("G", dim(BCov)[2]+1))


############################### plots of results ###############################
# plot of the beta-hats over time
png(file="img/3_beta-hats.png", width=640, height=320)
par(mar=c(4,5,1,1))
plot(betaHat0Vec, ylab=expression(hat(beta)^{(j)}),
     xlab="time segment middle year", xaxt="n", xlim=0.5+c(0,numRows))
axis(side=1, at=(1:numRows), labels=rev(rev(timeSegs[-1])[-1]))
dev.off()

# BCor and MCor matrix "heat maps"
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(BCor, breaks=seq(-1, 0.9, by=0.1),  border=T,
  col = c("black", paste(rep("grey", 18), seq(5, 95, 5), sep=""),
          "white")
)
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(MCor, breaks=seq(-1,0.9,by=0.1), border=T,
     col=c("black",paste(rep("grey",18), seq(5,95,5), sep=""),"white"))

# using green and red colours
par(mar=c(4.5, 4, 3.5, 4.1))
plot(BCor, breaks=seq(-1,1,by=0.1), border=T,
     col=c(hsv(2/3,seq(1,0,by=-0.1),1),hsv(0,seq(0.1,1,by=0.1),1)),
     main="Bartlett Correlation Estimates")
par(mar=c(4.5, 4, 3.5, 4.1))
plot(MCor, breaks=seq(-1,1,by=0.1), border=T,
     col=c(hsv(2/3,seq(1,0,by=-0.1),1),hsv(0,seq(0.1,1,by=0.1),1)),
     main="MTM Correlation Estimates")

pdf(file="img/3_BCor-MCor-vals.pdf", width=8, height=6)
par(mar=c(4.5, 4, 2, 4.1))
plot(BCor * lower.tri(BCor, diag=F) + MCor * upper.tri(MCor, diag=T),
     breaks=seq(-1,1,by=0.1), border=T,
     col=c(hsv(2/3,seq(1,0,by=-0.1),1),hsv(0,seq(0.1,1,by=0.1),1)),
     main="Bartlett/MTM Correlation Estimates", digits=2)
dev.off()

# differenced series of beta-hats
b.diff <- diff(betaHat0Vec, lag=1)
b.diff.2 <- diff(betaHat0Vec, lag=2)

MCov.diff <- matrix(NA, length(b.diff), length(b.diff))
MCov.diff.2 <- matrix(NA, length(b.diff.2), length(b.diff.2))

# MCor for the differenced series of beta-hats
MCor.diff <- matrix(NA, length(b.diff), length(b.diff))
for (j in 1:length(b.diff)) {
  for (l in 1:length(b.diff)) {
    MCov.diff[j,l] <- MCov[j+1,l+1] + MCov[j,l] - MCov[j+1,l] - MCov[j,l+1]
    MCor.diff[j,l] <- MCov.diff[j,l] /
      sqrt(MCov[j+1,j+1]+MCov[j,j]-2*MCov[j+1,j]) /
      sqrt(MCov[l+1,l+1]+MCov[l,l]-2*MCov[l+1,l])
  }
}
diag(MCor.diff) <- 1

pdf(file="img/3_MCor-diff-vals.pdf", width=8, height=6)
par(mar=c(4.5, 4, 3.5, 4.1))
plot(MCor.diff, breaks=seq(-1,1,by=0.1), border=T,
     col=c(hsv(2/3,seq(1,0,by=-0.1),1),hsv(0,seq(0.1,1,by=0.1),1)),
     main="MTM Correlation Estimates for Differenced beta-hat series",
     digits=2)
dev.off()

# MCor for the 2nd-differenced series of beta-hats
MCor.diff.2 <- matrix(NA, length(b.diff.2), length(b.diff.2))
for (j in 1:length(b.diff.2)) {
  for (l in 1:length(b.diff.2)) {
    MCov.diff.2[j,l] <- MCov.diff[j+1,l+1] + MCov.diff[j,l] - MCov.diff[j+1,l] - MCov.diff[j,l+1]
    MCor.diff.2[j,l] <- MCov.diff.2[j,l] /
      sqrt(MCov.diff[j+1,j+1]+MCov.diff[j,j]-2*MCov.diff[j+1,j]) /
      sqrt(MCov.diff[l+1,l+1]+MCov.diff[l,l]-2*MCov.diff[l+1,l])
  }
}
diag(MCor.diff.2) <- 1

detach(package:plot.matrix) # sometimes it's annoying when plotting column vectors



############################# confidence intervals #############################
CIs.mtap.names <- c("beta.hat","uncon.var","condl.var","CI.lower","CI.upper",
                    "new.CI.lower","new.CI.upper","CI.len","new.CI.len")
CIs.mtap <- as.data.frame(matrix(NA, numRows, length(CIs.mtap.names)))
names(CIs.mtap) <- CIs.mtap.names

# CIs for a beta-hat, conditional on all the others
CIs.mtap$beta.hat <- betaHat0Vec
CIs.mtap$uncon.var <- diag(MCov)
for (bh in 1:numRows) {
  CIs.mtap$condl.var[bh] <- MCov[bh,bh] - MCov[bh,-bh] %*% solve(MCov[-bh,-bh]) %*% MCov[-bh,bh]
}
CIs.mtap$CI.lower <- CIs.mtap$beta.hat - qnorm(0.975) * sqrt(CIs.mtap$uncon.var)
CIs.mtap$CI.upper <- CIs.mtap$beta.hat + qnorm(0.975) * sqrt(CIs.mtap$uncon.var)
CIs.mtap$new.CI.lower <- CIs.mtap$beta.hat - qnorm(0.975) * sqrt(CIs.mtap$condl.var)
CIs.mtap$new.CI.upper <- CIs.mtap$beta.hat + qnorm(0.975) * sqrt(CIs.mtap$condl.var)
CIs.mtap$CI.len <- CIs.mtap$CI.upper - CIs.mtap$CI.lower
CIs.mtap$new.CI.len <- CIs.mtap$new.CI.upper - CIs.mtap$new.CI.lower

middleyears <- rev(rev(timeSegs[-1])[-1])

# plot just the MTM-based old CIs (no cov info) and new CIs (using MCov)
png(file="img/3_beta-hats-CIs-mtap.png", width=640, height=320)
par(mar=c(4,5,1,1))
plot(x=middleyears, y=betaHat0Vec, type="p", pch=17,
     xlim=range(timeSegs)+c(0.5,-0.5),
     ylim=range(c(CIs.mtap$CI.lower,CIs.mtap$CI.upper)),
     xlab="time segment middle year", ylab=expression(hat(beta)^{(j)}))
arrows(x0=middleyears, y0=CIs.mtap$CI.lower, x1=middleyears,
       y1=CIs.mtap$CI.upper, length=0.05, angle=90, code=3, lwd=2, col="gray60")
arrows(x0=middleyears, y0=CIs.mtap$new.CI.lower, x1=middleyears,
       y1=CIs.mtap$new.CI.upper,length=0.05, angle=90, code=3, lwd=2, col=1)
dev.off()

# print LaTeX table
xtable::xtable(x=CIs.mtap, digits=4, display=rep("G", dim(CIs.mtap)[2]+1))


# create the equivalent of CIs.mtap, but for the Bartlett estimate
CIs.bart.names <- CIs.mtap.names
CIs.bart <- as.data.frame(matrix(NA, numRows, length(CIs.bart.names)))
names(CIs.bart) <- CIs.bart.names

# CIs for a beta-hat, conditional on all the others
CIs.bart$beta.hat <- betaHat0Vec
CIs.bart$uncon.var <- diag(BCov)
for (bh in 1:numRows) {
  CIs.bart$condl.var[bh] <- BCov[bh,bh] - BCov[bh,-bh] %*% solve(BCov[-bh,-bh]) %*% BCov[-bh,bh]
}
CIs.bart$CI.lower <- CIs.bart$beta.hat - qnorm(0.975) * sqrt(CIs.bart$uncon.var)
CIs.bart$CI.upper <- CIs.bart$beta.hat + qnorm(0.975) * sqrt(CIs.bart$uncon.var)
CIs.bart$new.CI.lower <- CIs.bart$beta.hat - qnorm(0.975) * sqrt(CIs.bart$condl.var)
CIs.bart$new.CI.upper <- CIs.bart$beta.hat + qnorm(0.975) * sqrt(CIs.bart$condl.var)
CIs.bart$CI.len <- CIs.bart$CI.upper - CIs.bart$CI.lower
CIs.bart$new.CI.len <- CIs.bart$new.CI.upper - CIs.bart$new.CI.lower

# requires creation of the CIs.bart data.frame
png(file="img/3_beta-hats-CIs-mtap-bart.png", width=640, height=320)
par(mar=c(4,5,1,1))
plot(x=middleyears, y=betaHat0Vec, type="p", pch=17,
     xlim=range(timeSegs)+c(0.5,-0.5),
     ylim=range(c(CIs.bart$CI.lower,CIs.bart$CI.upper)),
     xlab="time segment middle year", ylab=expression(hat(beta)^{(j)}))
arrows(x0=middleyears, y0=CIs.bart$CI.lower, x1=middleyears,
       y1=CIs.bart$CI.upper, length=0.05, angle=90, code=3, lwd=2, col="gray75")
arrows(x0=middleyears, y0=CIs.mtap$CI.lower, x1=middleyears,
       y1=CIs.mtap$CI.upper, length=0.05, angle=90, code=3, lwd=2, col="gray50")
arrows(x0=middleyears, y0=CIs.mtap$new.CI.lower, x1=middleyears,
       y1=CIs.mtap$new.CI.upper,length=0.05, angle=90, code=3, lwd=2, col=1)
legend("topright", lwd=rep(2,3), col=c("gray75","gray50","black"), title="CIs",
       legend=c("Bartlett, indep.", "MTM, indep.", "MTM, cond'l cov."), cex=0.8)
dev.off()

# condl.mean <- (MCov[numRows,numRows])^(-1) * MCov[numRows,-numRows] %*% (cbind(betaHat0Vec[-numRows]))
condl.var.bart <- BCov[numRows,numRows] - BCov[numRows,-numRows] %*% solve(BCov[-numRows,-numRows]) %*% BCov[-numRows,numRows]



