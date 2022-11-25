root open /home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/run_184.root
config set fit.quickfit.region 100
root get CebraE0
fit read run_184_det_0_fits.txt
root get CebraE1
fit read run_184_det_1_fits.txt
root get CebraE2
fit read run_184_det_2_fits.txt
root get CebraE3
fit read run_184_det_3_fits.txt
root get CebraE4
fit read run_184_det_4_fits.txt
calibration position list read run_184_cebra_cal.txt
