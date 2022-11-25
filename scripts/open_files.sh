#! /bin/bash

RUNNO=$1

# read -p "Please enter the run number: " RUNNO

touch hdtv_openfiles.sh
chmod +x hdtv_openfiles.sh

touch hdtv_savefiles.sh
chmod +x hdtv_savefiles.sh

echo "root open /home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/run_${RUNNO}.root" > hdtv_openfiles.sh
echo "config set fit.quickfit.region 150" >> hdtv_openfiles.sh

echo "s show all" > hdtv_savefiles.sh

#going to shift all detectors to common values 
for (( n=0; n<=4; n++))
do
    echo "root get CebraE${n}" >> hdtv_openfiles.sh

    echo "s show ${n}" >> hdtv_savefiles.sh
    #12C run 82 det 4 
    # echo "calibration position assign 0.0 435.3 1.0 588.9 2.0 772.9 3.0 966.0" >> hdtv_savefiles.sh 
    #49Ti run 184 det 4 
    echo "calibration position assign 0.0 494.4 0.1 532.2 1.0 647.2 2.0 880.0 3.0 1151.7 4.0 1430.7" >> hdtv_savefiles.sh 
    echo "fit write run_${RUNNO}_det_${n}_fits.txt" >> hdtv_savefiles.sh

done

echo "calibration position list write run_${RUNNO}_cebra_cal.txt" >> hdtv_savefiles.sh
echo "s show all" >> hdtv_savefiles.sh

hdtv -b hdtv_openfiles.sh

