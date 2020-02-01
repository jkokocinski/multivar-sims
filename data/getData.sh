#!/usr/bin/bash

if [ $1 == "demand" ]; then
  # demand data
  wget -nc http://reports.ieso.ca/public/Demand/PUB_Demand_{2002..2019}.csv
fi

if [ $1 == "hoep" ]; then
  # HOEP data
  wget -nc http://reports.ieso.ca/public/PriceHOEPPredispOR/PUB_PriceHOEPPredispOR_{2002..2019}.csv
fi

if [ $1 == "temperature" ]; then
  # tempearature data
  for year in `seq 2002 2019`
  do
  	for month in `seq 1 12`
  	do
  		wget -nc --content-disposition "https://climate.weather.gc.ca/climate_data/bulk_data_e.html?format=csv&stationID=5097&Year=${year}&Month=${month}&Day=14&timeframe=1&submit=Download+Data"
  	done
  done
fi

if [ $1 == "weather" ]; then
  # weather data from SSC
  wget -nc "https://ssc.ca/sites/default/files/imce/ssc2020_hourly_weather.xlsx"
fi
