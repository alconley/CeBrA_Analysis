import ROOT
import numpy as np
from scipy.integrate import quad
from hdtv_xmlfits_to_root import general_xml

inputfile = ROOT.TFile.Open("/home/alconley/Projects/CeBrA_Analysis/49Ti_dp_Histograms.root","READ")
# outputfile = ROOT.TFile.Open("/home/alconley/Projects/CeBrA_Analysis/50Ti_Excited_States.root", "RECREATE")

#declare file name energies
x1_E = "4880"

def hdtv_fit_to_root_guas(kev_per_bin,fit_filename):
    
    fit_path = "/home/alconley/Projects/CeBrA_Analysis/50Ti_Excited_States/fits/"
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

particle_gamma = inputfile.Get("49Ti_dpg/Energy_Calibrated/x1_cebraE_ECal_Array")
particle = inputfile.Get("49Ti_dpg/Energy_Calibrated/Array_x1_ECal_Protons")

excited_state_particle = inputfile.Get(f"49Ti_dpg/50Ti_Excited_States/Array_{x1_E}_keV_x1")
excited_state_gamma = inputfile.Get(f"49Ti_dpg/50Ti_Excited_States/Array_{x1_E}_keV_cebraE")

def GuasFit(histogram,x1,x2):
    fit = ROOT.TF1('fit','gaus',x1,x2)
    histogram.Fit(fit,'R')
    
    Volume = fit.GetParameter(0)
    Energy = fit.GetParameter(1)
    FWHM = 2.355 * fit.GetParameter(2)

    # print("Fit Energy: ", Energy)
    # print("Fit FWHM: ", FWHM)
    return Volume, Energy, FWHM

c = ROOT.TCanvas(f"50Ti_{x1_E}_keV_Excited_State",f"50Ti_{x1_E}_keVExcited_State",1400,1000)

c.Divide(2,2,0.001,0.001)

c.cd(1)
particle_gamma.Draw("colz")
particle_gamma.SetTitle("Particle-#gamma Coincidence Matrix; 50Ti Energy [keV]; #gamma-ray Energy [keV];")

c.cd(2)
particle.Draw()
particle.SetTitle("Partcle Spectrum; 50Ti Energy [keV]; ;")

c.cd(3)
excited_state_gamma.Draw()
excited_state_gamma.SetTitle(f"{x1_E} keV Excited State Gamma Spectrum; gamma-ray Energy [keV]; ;")

cebra_fits,cebra_fit_parameters = hdtv_fit_to_root_guas(13,f"{x1_E}_keV_cebra.fit")
for i in range(len(cebra_fits)):
    cebra_fits[i].Draw("same")
    print(f"This is index: {i}. [Energy,Width,Volume] = {cebra_fit_parameters[i]}")

def I_gamma(relative_gamma_index, other_gamma_index):
    #Gamma energy of the peak that the intensity are relative to
    E_j = cebra_fit_parameters[relative_gamma_index][0]
    V_j = abs(cebra_fit_parameters[relative_gamma_index][2])

    E_i = cebra_fit_parameters[other_gamma_index][0]
    V_i = abs(cebra_fit_parameters[other_gamma_index][2])

    I_gamma = (V_i/CeBrA_Efficiency(E_i)) / (V_j/CeBrA_Efficiency(E_j)) * 100
    return I_gamma

I_2205 = I_gamma(0,0)
I_1681 = I_gamma(0,1)
I_733 = I_gamma(0,4)

print(f"\n 50Ti: {x1_E} keV Excited state")
print(f"I(gamma) = I(2205 keV) = {I_2205}")
print(f"I(gamma) = I(1681 keV) = {I_1681}")
print(f"I(gamma) = I(733 keV) = {I_733}")

print(CeBrA_Efficiency(3826))
print(CeBrA_Efficiency(2704))


c.cd(4)
excited_state_particle.Draw()
excited_state_particle.SetTitle(f"{x1_E} keV Excited State Particle Spectrum; 50Ti Energy [keV]; ;")

x1_fits,x1_fit_parameters = hdtv_fit_to_root_guas(16,f"{x1_E}_keV_x1.fit")
for i in range(len(x1_fits)):
    x1_fits[i].Draw("same")
    

c.Update()

# c.Write(f"50Ti_{x1_E[i]}_keV_Excited_State")
c.Draw()

input("hit enter when done")