#!/usr/bin/env Rscript

tryCatch(setwd("~/multivar-sims"),
         error=function(e) {
           tryCatch( setwd("~/WORK/Q5/masters-thesis/code/multivar-sims") )
         }
)

require(multitaper)

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

iesoData <- merge.data.frame(x=demand[,-(1:2)], y=hoep[,-(1:2)], sort=TRUE)
iesoData <- iesoData[,c("datetime", "Ontario.Demand", "HOEP")]
names(iesoData)[2:3] <- c("demand", "hoep")
remove(demand, hoep) # clean up
iesoData$hoep <- gsub(",", "", iesoData$hoep) # remove thousands separator
iesoData$hoep <- as.numeric(iesoData$hoep)

plot(x=iesoData$datetime, y=iesoData$demand, type="l")
plot(x=iesoData$datetime, y=iesoData$hoep, type="l")

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


