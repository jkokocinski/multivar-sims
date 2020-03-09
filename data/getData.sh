#!/usr/bin/bash

echo Which data?
echo \[ demand \| hoep \| temperature \| weather \]
read dataName

if [ $dataName == "demand" ]; then
  # demand data
  cd demand
  wget -nc http://reports.ieso.ca/public/Demand/PUB_Demand_{2002..2019}.csv
  cd ..
fi

if [ $dataName == "hoep" ]; then
  # HOEP data
  cd hoep
  wget -nc http://reports.ieso.ca/public/PriceHOEPPredispOR/PUB_PriceHOEPPredispOR_{2002..2019}.csv
  cd ..
fi

if [ $dataName == "temperature" ]; then
  # tempearature data
  cd temperature
  for year in `seq 2002 2019`
  do
  	for month in `seq 1 12`
  	do
  		wget -nc --content-disposition "https://climate.weather.gc.ca/climate_data/bulk_data_e.html?format=csv&stationID=5097&Year=${year}&Month=${month}&Day=14&timeframe=1&submit=Download+Data"
  	done
  done
  cd ..
fi

if [ $dataName == "weather" ]; then
  # weather data from SSC
  cd ssc
  wget -nc "https://ssc.ca/sites/default/files/imce/ssc2020_hourly_weather.xlsx"
  cd ..
fi
