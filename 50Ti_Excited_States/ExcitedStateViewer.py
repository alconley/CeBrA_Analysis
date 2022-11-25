import ROOT

inputfile = ROOT.TFile.Open("/home/alconley/Projects/CeBrA_Analysis/49Ti_dp_Histograms.root","READ")
outputfile = ROOT.TFile.Open("/home/alconley/Projects/CeBrA_Analysis/50Ti_Excited_States.root", "RECREATE")

#declare file name energies
x1_E = ["0","1553","2674","4100_doublet","4500_doublet","4880","5186",\
        "5379","5500_5630","5660_5780","5837","5946","6020_6100","6123"]

for i in range(0,len(x1_E)):

    particle_gamma = inputfile.Get("49Ti_dpg/Energy_Calibrated/x1_cebraE_ECal_Array")
    particle = inputfile.Get("49Ti_dpg/Energy_Calibrated/Array_x1_ECal_Protons")
    
    excited_state_particle = inputfile.Get(f"49Ti_dpg/50Ti_Excited_States/Array_{x1_E[i]}_keV_x1")
    excited_state_gamma = inputfile.Get(f"49Ti_dpg/50Ti_Excited_States/Array_{x1_E[i]}_keV_cebraE")

    c = ROOT.TCanvas(f"50Ti_{x1_E[i]}_keV_Excited_State",f"50Ti_{x1_E[i]}_keVExcited_State",1400,1000)

    c.Divide(2,2,0.001,0.001)

    c.cd(1)
    particle_gamma.Draw("colz")
    particle_gamma.SetTitle("Particle-#gamma Coincidence Matrix; 50Ti Energy [keV]; #gamma-ray Energy [keV];")

    c.cd(2)
    particle.Draw()
    particle.SetTitle("Partcle Spectrum; 50Ti Energy [keV]; ;")

    c.cd(3)
    excited_state_gamma.Draw()
    excited_state_gamma.SetTitle(f"{x1_E[i]} keV Excited State Gamma Spectrum; gamma-ray Energy [keV]; ;")

    c.cd(4)
    excited_state_particle.Draw()
    excited_state_particle.SetTitle(f"{x1_E[i]} keV Excited State Particle Spectrum; 50Ti Energy [keV]; ;")

    c.Update()
    
    c.Write(f"50Ti_{x1_E[i]}_keV_Excited_State")

    # input("hit enter when done")