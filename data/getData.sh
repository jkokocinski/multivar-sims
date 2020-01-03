#!/usr/bin/bash

if [ $1 == "demand" ]; then
  # demand data
  wget http://reports.ieso.ca/public/Demand/PUB_Demand_{2002..2018}.csv
fi

if [ $1 == "hoep" ]; then
  # HOEP data
  wget http://reports.ieso.ca/public/PriceHOEPPredispOR/PUB_PriceHOEPPredispOR_{2002..2018}.csv
fi

# tempearature data
if [ $1 == "temperature" ]; then
  for year in `seq 2002 2018`
  do
  	for month in `seq 1 12`
  	do
  		wget --content-disposition "https://climate.weather.gc.ca/climate_data/bulk_data_e.html?format=csv&stationID=5097&Year=${year}&Month=${month}&Day=14&timeframe=1&submit=Download+Data"
  	done
  done
fi
