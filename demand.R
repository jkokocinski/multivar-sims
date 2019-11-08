#!/usr/bin/env Rscript

tryCatch(setwd("~/multivar-sims"),
         error=function(e) {
           tryCatch( setwd("~/WORK/Q5/masters-thesis/code/multivar-sims") )
         }
)

require(multitaper)
require(signal)
source("seasonalFunctions.R")

############################ import data and format ############################
importData <- function(prefix, seqInds, suffix="", filext, skip=0) {
  theData <- data.frame()
  
  for (ind in seqInds) {
  	newData <- read.csv(file=paste0(prefix, ind, suffix, filext),
  	                    header=TRUE, stringsAsFactors=FALSE, skip=skip)
  	theData <- rbind(theData, newData)
  }
  return(theData)
}

demand <- importData("data/PUB_Demand_", 2002:2018, "", ".csv", skip=3)
hoep <- importData("data/PUB_PriceHOEPPredispOR_", 2002:2018, "", ".csv", skip=3)

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


################################# plot the data ################################
png(file="img/demand-hoep_plots.png", width=1024, height=576)
par(mfrow=c(2,1))
par(mar=c(4,4,1,1))
plot(x=iesoData$datetime, y=iesoData$demand, type="l",
     xlab="Date", ylab="Ontario Demand (kW)")
plot(x=iesoData$datetime, y=iesoData$hoep, type="l",
     xlab="Date", ylab="HOEP")
dev.off()

######################## compute auto- and cross-spectra #######################
multitaper::spec.mtm(ts(iesoData$demand))
multitaper::spec.mtm(ts(iesoData$hoep))

# compute spec.mtm objects for both response series components
spec.y1 <- multitaper::spec.mtm(
  timeSeries=ts(iesoData$demand), k=6,
  nw=11,
  adaptiveWeighting=TRUE,
  Ftest=TRUE, returnInternals=TRUE, plot=FALSE
)
spec.y2 <- multitaper::spec.mtm(
  timeSeries=ts(iesoData$hoep), k=6,
  nw=11,
  adaptiveWeighting=TRUE,
  Ftest=TRUE, returnInternals=TRUE, plot=FALSE
)

# wieghts and eigencoefficients
d1 <- spec.y1$mtm$eigenCoefWt
d2 <- spec.y2$mtm$eigenCoefWt
y1 <- spec.y1$mtm$eigenCoefs
y2 <- spec.y2$mtm$eigenCoefs

# compute the cross spectrum of (Y_1,Y_2) -- adaptive weighted eigencoefs
crossSpecEstY <- apply(d1*y1*d2*Conj(y2), MARGIN=1, FUN=sum) /
  apply(d1*d2, MARGIN=1, FUN=sum)

# cross spectrum over full [0,1) interval
crossSpecEstY <- c(crossSpecEstY, rev(Conj(crossSpecEstY[-1]))[-1])

# F-statistics for demand and HOEP
par(mar=c(4,4,1,1))
plot(x=spec.y1$freq, y=spec.y1$mtm$Ftest, type="l") # demand
par(mar=c(4,4,1,1))
plot(x=spec.y2$freq, y=spec.y2$mtm$Ftest, type="l") # HOEP


# plot amplitude and phase cross-spectra
png(file="img/demand-hoep_cross-spec.png", width=1024, height=576)
par(mfrow=c(2,1))
plotInds <- seq(1,length(crossSpecEstY)/2) # for plotting on [0,0.5)
freqGrid <- (plotInds-1)/length(crossSpecEstY)
# amplitude cross spectrum estimate
par(mar=c(4,4,1,1))
plot(x=freqGrid, y=Mod(crossSpecEstY[plotInds]), type="l", log="y",
     xlab="Frequency", ylab="Amplitude Cross Spectrum Estimate")
abline(v=seq(0,0.5, by=1/24), col="blue", lty=2, lwd=1) # daily, bi-daily, etc.
abline(v=seq(0,0.5, by=1/(24*7)), col="orange", lty=2) # weekly, bi-weely, etc.
# phase cross spectrum estimate
par(mar=c(4,4,1,1))
plot(x=freqGrid,
     y=signal::unwrap(atan2(Im(crossSpecEstY),Re(crossSpecEstY))[plotInds]),
     type="l", xlab="Frequency", ylab="Phase Cross Spectrum Estimate")
dev.off()



######################## deterministic signal detection ########################
system.time(
  {
    seas.demand <- determineSeasonal(data=iesoData$demand, sigCutoff=0.999,
                                       padFactor=2)
  }
)

system.time(
  {
    seas.hoep <- determineSeasonal(data=iesoData$hoep, sigCutoff=0.999,
                                       padFactor=2)
  }
)



################################# daily dataset ################################
iesoData.d <- subset(x=iesoData, subset=(hour==1))
row.names(iesoData.d) <- as.character(seq(1,dim(iesoData.d)[1])) # correct row.names
spec.y1.d <- multitaper::spec.mtm(
  timeSeries=ts(iesoData.d$demand), k=6, nw=11,
  adaptiveWeighting=TRUE,
  Ftest=TRUE, returnInternals=TRUE, plot=TRUE
)


################################ low-pass filter ###############################
filterSleps <- multitaper::dpss(n=dim(iesoData)[1], nw=trunc(dim(iesoData)[1]*(1/24/7/4)), k=trunc(2*dim(iesoData)[1]*(1/24/7/4)))
filtSlepComb <- filterSleps$v %*% cbind(filterSleps$eigen)
plot(Mod(fft(filtSlepComb))**2, type="l")














