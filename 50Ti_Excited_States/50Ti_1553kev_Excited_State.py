import ROOT
import numpy as np
from scipy.integrate import quad
from hdtv_xmlfits_to_root import general_xml

inputfile = ROOT.TFile.Open("/home/alconley/Projects/CeBrA_Analysis/49Ti.root","READ")
outputfile = ROOT.TFile.Open("/home/alconley/Projects/CeBrA_Analysis/50Ti_Excited_States.root", "RECREATE")

def hdtv_fit_to_root_guas(kev_per_bin,fit_path,fit_filename):
    
    uncal_list, cal_list = general_xml(f"{fit_path}{fit_filename}",f"{fit_filename}")
    
    fits = []
    fit_paramaters = []
    for i in range(len(cal_list)):
        Energy = cal_list[i][0]
        Energy_Uncertainity = cal_list[i][1]
        FWHM = cal_list[i][2]
        FWHM_Uncertainity = cal_list[i][3]
        Volume = cal_list[i][4]
        Volume_Uncertainity = cal_list[i][5]

        Sigma =  FWHM/2.355
        
        def f(x):
            return np.exp(-0.5*((x-Energy)/Sigma)**2)
        res,err = quad(f,Energy-(3*FWHM),Energy+(3*FWHM))
        
        Constant = (Volume/res)*kev_per_bin
        
        fit = ROOT.TF1('fit','gaus',0,6500)
        fit.SetParameters(Constant,Energy,Sigma)
        fits.append(fit)
        
        hdtv_fit = [Energy,FWHM,Volume]
        fit_paramaters.append(hdtv_fit)
        
    return fits,fit_paramaters

def CeBrA_Efficiency(Gamma_Energy):
    detector0 = 0.21928639275535525 * np.exp(-0.0002981977479889143*Gamma_Energy) + 1.063323409320861 * np.exp(-0.0021900446692233407*Gamma_Energy)
    detector4 = 0.787874293808005 * np.exp(-0.0002592572815672215*Gamma_Energy) + 1.2384884551061637 * np.exp(-0.0021335254959945443*Gamma_Energy)
    Efficency = (4*detector0 )+ (1*detector4)
    return Efficency

def GuasFit(histogram,x1,x2):
    fit = ROOT.TF1('fit','gaus',x1,x2)
    histogram.Fit(fit,'R')
    
    Volume = fit.GetParameter(0)
    Energy = fit.GetParameter(1)
    FWHM = 2.355 * fit.GetParameter(2)

    # print("Fit Energy: ", Energy)
    # print("Fit FWHM: ", FWHM)
    return Volume, Energy, FWHM

#if no fit files, type 0
def excited_state(Excited_State_Name,fit_file_path,cebraE_fit_file,particle_fit_file):
    
    particle_gamma = inputfile.Get("49Ti_Protons/CeBrA_Energy_Calibrated/x1_cebraE_ECal_Array")
    particle = inputfile.Get("49Ti_Protons/CeBrA_Energy_Calibrated/Array_x1_ECal_Protons")

    excited_state_particle = inputfile.Get(f"49Ti_Protons/50Ti_Excited_States/{Excited_State_Name}_keV/{Excited_State_Name}_keV_x1_ECal")
    excited_state_gamma = inputfile.Get(f"49Ti_Protons/50Ti_Excited_States/{Excited_State_Name}_keV/{Excited_State_Name}_keV_CeBrA_ECal")
    
    canvas = ROOT.TCanvas(f"50Ti_{Excited_State_Name}_keV_Excited_State",f"50Ti_{Excited_State_Name}_keV_Excited_State",1400,1000)

    canvas.Divide(2,2,0.001,0.001)
    canvas.Draw()

    #particle-gamma matrix
    canvas.cd(1)
    particle_gamma.Draw("colz")
    particle_gamma.SetTitle("Particle-#gamma Coincidence Matrix; 50Ti Energy [keV]; #gamma-ray Energy [keV];")

    #focal plane spectrum
    canvas.cd(2)
    particle.Draw()
    particle.SetTitle("Partcle Spectrum; 50Ti Energy [keV]; ;")

    #exctied state gamma spectrum
    canvas.cd(3)
    excited_state_gamma.Draw()
    excited_state_gamma.SetTitle(f"{Excited_State_Name} keV Excited State Gamma Spectrum; gamma-ray Energy [keV]; ;")

    if cebraE_fit_file != 0:
        
        cebra_fits,cebra_fit_parameters = hdtv_fit_to_root_guas(13,f"{fit_file_path}",f"{Excited_State_Name}_keV_CeBrA_ECal.fit")
        for i in range(len(cebra_fits)):
            cebra_fits[i].Draw("same")
            print(f"This is index: {i}. [Energy,Width,Volume] = {cebra_fit_parameters[i]}")
    else:
        cebra_fits=0
        cebra_fit_parameters = 0
        

    canvas.cd(4)
    excited_state_particle.Draw()
    excited_state_particle.SetTitle(f"{Excited_State_Name} keV Excited State Particle Spectrum; 50Ti Energy [keV]; ;")

    if particle_fit_file != 0:
        x1_fits,x1_fit_parameters = hdtv_fit_to_root_guas(16,f"{fit_file_path}",f"{Excited_State_Name}_keV_x1_ECal.fit")
        for i in range(len(x1_fits)):
            x1_fits[i].Draw("same")
    else:
        x1_fits = 0
        x1_fit_parameters = 0
    
    canvas.Update()
    
    return canvas,cebra_fits,cebra_fit_parameters,x1_fits,x1_fit_parameters
             
# def I_gamma(relative_gamma_index, other_gamma_index):
#     #Gamma energy of the peak that the intensity are relative to
#     E_j = cebra_fit_parameters[relative_gamma_index][0]
#     V_j = abs(cebra_fit_parameters[relative_gamma_index][2])

#     E_i = cebra_fit_parameters[other_gamma_index][0]
#     V_i = abs(cebra_fit_parameters[other_gamma_index][2])

#     I_gamma = (V_i/CeBrA_Efficiency(E_i)) / (V_j/CeBrA_Efficiency(E_j)) * 100
#     return I_gamma
    
def excited_state_bands(Excited_State_Name, Band_Name, fit_file_path, cebraE_fit_file, particle_fit_file):
    
    canvas = ROOT.TCanvas(f"50Ti_{Excited_State_Name}_keV_Excited_State_{Band_Name}_Band",f"50Ti_{Excited_State_Name}_keV_Excited_State_{Band_Name}_Band",1400,1000)
    canvas.Divide(2,2,0.001,0.001)

    excited_state_particle_gamma_band = inputfile.Get(f"49Ti_Protons/50Ti_Excited_States/{Excited_State_Name}_keV/{Excited_State_Name}_kev_{Band_Name}_Band_cebraE_x1")
    excited_state_particle_band = inputfile.Get(f"49Ti_Protons/50Ti_Excited_States/{Excited_State_Name}_keV/{Excited_State_Name}_kev_{Band_Name}_Band_x1")
    excited_state_gamma_band = inputfile.Get(f"49Ti_Protons/50Ti_Excited_States/{Excited_State_Name}_keV/{Excited_State_Name}_kev_{Band_Name}_Band_cebraE")

    canvas.cd(1)
    excited_state_particle_gamma_band.Draw("colz")
    excited_state_particle_gamma_band.SetTitle("Particle-#gamma Coincidence Matrix; 50Ti Energy [keV]; #gamma-ray Energy [keV];")

    canvas.cd(2)
    excited_state_particle_band.Draw()
    excited_state_particle_band.SetTitle(f"{Excited_State_Name} keV Excited State {Band_Name} Band Particle Spectrum; 50Ti Energy [keV]; ;")
    
    if particle_fit_file != 0:
        x1_band_fits,x1_band_fit_parameters = hdtv_fit_to_root_guas(16,f"{fit_file_path}",f"{Excited_State_Name}_keV_{Band_Name}_Band_x1.fit")
        for i in range(len(x1_band_fits)):
            x1_band_fits[i].Draw("same")
    else:
        x1_band_fits = 0
        x1_band_fit_parameters = 0

    canvas.cd(3)
    excited_state_gamma_band.Draw()
    excited_state_gamma_band.SetTitle(f"{Excited_State_Name} keV Excited State {Band_Name} Band Gamma Spectrum; gamma-ray Energy [keV]; ;")

    if cebraE_fit_file != 0:
        cebra_band_fits, cebra_band_fit_parameters = hdtv_fit_to_root_guas(13,f"{fit_file_path}",f"{Excited_State_Name}_keV_{Band_Name}_Band_cebraE.fit")
        for i in range(len(cebra_band_fits)):
            cebra_band_fits[i].Draw("same")
            print(f"{Band_Name} Band gamma fit index: {i}. [Energy,Width,Volume] = {cebra_band_fit_parameters[i]}")
    else:
        cebra_band_fits = 0
        cebra_band_fit_parameters = 0
            
    canvas.Update()
    
    return canvas, cebra_band_fits, cebra_band_fit_parameters, x1_band_fits, x1_band_fit_parameters

Excited_State_Name = "4100_doublet"

Fit_Path_Excited_State = f"/home/alconley/Projects/CeBrA_Analysis/50Ti_Excited_States/fits/{Excited_State_Name}_keV/"
Excited_State_Canvas, cebra_fits, cebra_fit_para, x1_fits, x1_fit_para = excited_state(Excited_State_Name,Fit_Path_Excited_State,0,0)

Excited_State_GS_Band_Canvas, cebra_GS_Band_fits, cebra_GS_Band_fit_para, x1_GS_Band_fits, x1_GS_Band_fit_para = excited_state_bands(Excited_State_Name,"GS",Fit_Path_Excited_State,0,0)
Excited_State_2_Plus_Band_Canvas, cebra_2_Plus_Band_fits, cebra_2_Plus_Band_fit_para, x1_2_Plus_Band_fits, x1_2_Plus_Band_fit_para = excited_state_bands(Excited_State_Name,"2_Plus",Fit_Path_Excited_State,1,0)
Excited_State_4_Plus_Band_Canvas, cebra_4_Plus_Band_fits, cebra_4_Plus_Band_fit_para, x1_4_Plus_Band_fits, x1_4_Plus_Band_fit_para = excited_state_bands(Excited_State_Name,"4_Plus",Fit_Path_Excited_State,1,0)
# Excited_State_6_Plus_Band_Canvas, cebra_6_Plus_Band_fits, cebra_6_Plus_Band_fit_para, x1_6_Plus_Band_fits, x1_6_Plus_Band_fit_para = excited_state_bands(Excited_State_Name,"6_Plus",Fit_Path_Excited_State,0,0)



print(cebra_4_Plus_Band_fit_para)

input("")