#! /bin/bash

RUNNO=$1

# read -p "Please enter the run number: " RUNNO

touch hdtv_openfiles_fitted.sh
chmod +x hdtv_openfiles_fitted.sh

echo "root open /home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/run_${RUNNO}.root" > hdtv_openfiles_fitted.sh
echo "config set fit.quickfit.region 100" >> hdtv_openfiles_fitted.sh


#going to shift all detectors to common values 
for (( n=0; n<=4; n++))
do
    echo "root get CebraE${n}" >> hdtv_openfiles_fitted.sh
    echo "fit read run_${RUNNO}_det_${n}_fits.txt" >> hdtv_openfiles_fitted.sh

done

echo "calibration position list read run_${RUNNO}_cebra_cal.txt" >> hdtv_openfiles_fitted.sh

hdtv -b hdtv_openfiles_fitted.sh

