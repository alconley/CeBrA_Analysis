import ROOT

inputfile = ROOT.TFile.Open("/home/alconley/Projects/CeBrA_Analysis/49Ti_dp_Histograms.root","READ")
outputfile = ROOT.TFile.Open("/home/alconley/Projects/CeBrA_Analysis/50Ti_Bands.root", "RECREATE")

#declare file name energies
bands = ["GS","2_Plus","4_Plus"]

for i in range(0,len(bands)):

    particle_gamma = inputfile.Get("49Ti_dpg/Energy_Calibrated/x1_cebraE_ECal_Array")
    
    particle_gamma_band = inputfile.Get(f"49Ti_dpg/{bands[i]}_Band/Array_cebraE_x1_ECal_{bands[i]}_Band")
    
    band_particle = inputfile.Get(f"49Ti_dpg/{bands[i]}_Band/Array_x1_ECal_{bands[i]}_Band")
    band_gamma = inputfile.Get(f"49Ti_dpg/{bands[i]}_Band//Array_cebraE_{bands[i]}_Band")

    c = ROOT.TCanvas(f"50Ti_{bands[i]}_Band",f"50Ti_{bands[i]}_Band",1400,1000)

    c.Divide(2,2,0.001,0.001)

    c.cd(1)
    particle_gamma.Draw("colz")
    particle_gamma.SetTitle("Particle-#gamma Coincidence Matrix; 50Ti Energy [keV]; #gamma-ray Energy [keV];")

    c.cd(2)
    particle_gamma_band.Draw("colz")
    particle_gamma_band.SetTitle(f"Particle-#gamma Coincidence Matrix {bands[i]} Band; 50Ti Energy [keV]; #gamma-ray Energy [keV];")

    c.cd(3)
    band_gamma.Draw()
    band_gamma.SetTitle(f"{bands[i]} Band Gamma Spectrum; gamma-ray Energy [keV]; ;")

    c.cd(4)
    band_particle.Draw()
    band_particle.SetTitle(f"{bands[i]} Band Particle Spectrum; 50Ti Energy [keV]; ;")

    c.Update()
    
    c.Write(f"50Ti_{bands[i]} Band")

#     input("hit enter when done")