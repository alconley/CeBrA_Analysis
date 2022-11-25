root open 49Ti_Histograms.root

root get 49Ti_dpg/x1_only1plane_Protons
fit read 
calibration position list read x1_only1plane_cal.txt
fit read 49Ti_x1_fits.txt
fit show decomposition 0-100

root get 49Ti_dpg/50Ti_Excited_States/Array_0_keV_cebraE
root get 49Ti_dpg/50Ti_Excited_States/Array_0_keV_x1

root get 49Ti_dpg/50Ti_Excited_States/Array_1553_keV_cebraE
fit read Array_1553_keV_cebraE.fit
root get 49Ti_dpg/50Ti_Excited_States/Array_1553_keV_x1
fit read Array_1553_keV_x1.fit

root get 49Ti_dpg/50Ti_Excited_States/Array_2674_keV_cebraE
fit read Array_2674_keV_cebraE.fit
root get 49Ti_dpg/50Ti_Excited_States/Array_2674_keV_x1

root get 49Ti_dpg/50Ti_Excited_States/Array_4100_doublet_keV_cebraE
fit read Array_4100_doublet_keV_cebraE.fit
root get 49Ti_dpg/50Ti_Excited_States/Array_4100_doublet_keV_x1

root get 49Ti_dpg/50Ti_Excited_States/Array_4500_doublet_keV_cebraE
fit read Array_4500_doublet_keV_cebraE.fit
root get 49Ti_dpg/50Ti_Excited_States/Array_4500_doublet_keV_x1

root get 49Ti_dpg/50Ti_Excited_States/Array_4880_keV_cebraE
fit read Array_4880_keV_cebraE.fit
root get 49Ti_dpg/50Ti_Excited_States/Array_4880_keV_x1


# root get 49Ti_dpg/50Ti_Excited_States/Array_5186_keV_cebraE
# root get 49Ti_dpg/50Ti_Excited_States/Array_5379_keV_cebraE
# root get 49Ti_dpg/50Ti_Excited_States/Array_5500_5630_keV_cebraE
# root get 49Ti_dpg/50Ti_Excited_States/Array_5660_5780_keV_cebraE
# root get 49Ti_dpg/50Ti_Excited_States/Array_5837_keV_cebraE
# root get 49Ti_dpg/50Ti_Excited_States/Array_5946_keV_cebraE
# root get 49Ti_dpg/50Ti_Excited_States/Array_6020_6100_keV_cebraE
# root get 49Ti_dpg/50Ti_Excited_States/Array_6123_keV_cebraE
