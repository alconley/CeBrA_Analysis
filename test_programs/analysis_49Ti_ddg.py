import ROOT
import math
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer


Ti49_RunList = [184,185,186,187,188,190,191,192,193,\
    194,195,196,197,198,199,200,201,202,203,204,205,\
        206,207,208,209,210,211,212,213,214,215,216,\
            217,218,219,220,221,222,223,224,225,226,\
                227,228,229,230,231,232]
                

# Ti49_RunList=[184]

for RUNNO in Ti49_RunList:
# for RUNNO in C12_RunList:
        
    print("Current Run: ",RUNNO)
    start = timer()

    #where the analyzed files/tree is located
    # PATH = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/12C_analyzed/"
    PATH = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/"

    OUT = "run_" +str(RUNNO)+"_histo_49Ti_ddg.root"
    output = ROOT.TFile.Open(OUT, "RECREATE")

    #number of detectors
    n_det = 5
    
    treeName = "SPSTree"

    #PLOTS SHOULD BE DRAWN AT THE END

    #optimize the number of thread before you run the loop
    ROOT.ROOT.EnableImplicitMT(5)
    
######################################################################################################################################################
    
    #gain matching calibration numbers, got these with a HDTV Script
    #need to parse the HDTV calibration file
    intercept_GM = []
    slope_GM = []
    # with open("/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/12C_analyzed/run_"+str(RUNNO)+"_cebra_cal.txt","r") as f:
    with open("/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/run_"+str(RUNNO)+"_cebra_cal.txt","r") as f:
        stripped=[s.strip() for s in f]
        for line in stripped:
            name, intercept, slope=line.split()
            intercept_GM.append(intercept)
            slope_GM.append(slope)
    
    #Intitial Frame
    #Need to gain match all the cebra detectors before this program
    Raw_Frame = ROOT.RDataFrame(treeName, PATH + "run_"+str(RUNNO)+".root")\
                    .Define("cebraTime0_toScint","cebraTime0-scintLeftTime")\
                    .Define("cebraTime1_toScint","cebraTime1-scintLeftTime")\
                    .Define("cebraTime2_toScint","cebraTime2-scintLeftTime")\
                    .Define("cebraTime3_toScint","cebraTime3-scintLeftTime")\
                    .Define("cebraTime4_toScint","cebraTime4-scintLeftTime")\
                    .Define("cebraE0_GM","cebraE0*"+slope_GM[0]+"+"+intercept_GM[0])\
                    .Define("cebraE1_GM","cebraE1*"+slope_GM[1]+"+"+intercept_GM[1])\
                    .Define("cebraE2_GM","cebraE2*"+slope_GM[2]+"+"+intercept_GM[2])\
                    .Define("cebraE3_GM","cebraE3*"+slope_GM[3]+"+"+intercept_GM[3])\
                    .Define("cebraE4_GM","cebraE4*"+slope_GM[4]+"+"+intercept_GM[4])


######################################################################################################################################################

    #Gets the Particle Identification
    PID = Raw_Frame.Filter("scintLeftTime != -1 && anodeBackTime != -1").Histo2D(("ScintL_AnodeBack","PID",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
    # c_PID = ROOT.TCanvas()
    # PID.Draw("colz") 

######################################################################################################################################################

    #49Ti, run 184 through 232
    Deuterons = "scintLeft > 200 && scintLeft < 550 && anodeBack > 300 && anodeBack < 1100"
    
    PID_Deuterons = Raw_Frame.Filter(Deuterons).Histo2D(("ScintL_AnodeBack_Deuterons","Deuterons",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
        
######################################################################################################################################################
                                                                                                
    #cebra plots gained matched and raw cebra
    cebraE_noCuts = ["cebraE0_noCuts","cebraE1_noCuts","cebraE2_noCuts","cebraE3_noCuts","cebraE4_noCuts"]
    cebraE_GM_noCuts = ["cebraE0_GM_noCuts","cebraE1_GM_noCuts","cebraE2_GM_noCuts","cebraE3_GM_noCuts","cebraE4_GM_noCuts"]
    cebraTime_toScint_noCuts = ["cebraTime0_toScint_noCuts","cebraTime1_toScint_noCuts","cebraTime2_toScint_noCuts","cebraTime3_toScint_noCuts","cebraTime4_toScint_noCuts"]
    
    for i in range(n_det):
        cebraE_noCuts[i] = Raw_Frame.Filter(f"cebraE{i}>-1").Histo1D((f"cebraE{i}_noCuts",f"cebraE{i}_noCuts",4096,0,4096),f"cebraE{i}")
        cebraE_GM_noCuts[i] = Raw_Frame.Filter(f"cebraE{i}>-1").Histo1D((f"cebraE{i}_GM_noCuts",f"cebraE{i}_GM_noCuts",1024,0,4096),f"cebraE{i}_GM")
        cebraTime_toScint_noCuts[i] = Raw_Frame.Filter(f"scintLeftTime != -1  && anodeBackTime != -1 && cebraE{i} != -1").Histo1D((f"cebraTime{i}_toScint_noCuts",f"cebraTime{i}_toScint_noCuts",6400,-3200,3200),f"cebraTime{i}_toScint")
            
######################################################################################################################################################

    #Creating a new frames for different particle groups
    Deuteron_Frame = Raw_Frame.Filter(Deuterons)
    
    xavg_Deuterons = Deuteron_Frame.Histo1D(("xavg_Deuterons","xavg_Deuterons",600,-300,300),"xavg")
    x1_bothplanes_Deuterons = Deuteron_Frame.Histo1D(("x1_bothplanes_Deuterons","x1_bothplanes_Deuterons",600,-300,300),"x1")
    x2_bothplanes_Deuterons = Deuteron_Frame.Histo1D(("x2_bothplanes_Deuterons","x2_bothplanes_Deuterons",600,-300,300),"x2")
    x1_only1plane_Deuterons = Deuteron_Frame.Filter("x1 != -1e6 && x2 == -1e6").Histo1D(("x1_only1plane_Deuterons","x1_only1plane_Deuterons",600,-300,300),"x1")
    
    
    cebraTime_toScint_Deuterons = ["cebraTime0_toScint_Deuterons","cebraTime1_toScint_Deuterons","cebraTime2_toScint_Deuterons","cebraTime3_toScint_Deuterons","cebraTime4_toScint_Deuterons"]
    for i in range(n_det):
        cebraTime_toScint_Deuterons[i] = Raw_Frame.Filter(Deuterons).Histo1D((f"cebraTime{i}_toScint_Deuterons",f"cebraTime{i}_toScint_Deuterons",6400,-3200,3200),f"cebraTime{i}_toScint")


    #49Ti
    cebraTime0_gate_deuterons = "cebraTime0_toScint > -733 && cebraTime0_toScint < -720"    #Deuterons_TCut0
    cebraTime1_gate_deuterons = "cebraTime1_toScint > -733 && cebraTime1_toScint < -717"    #Deuterons_TCut1
    cebraTime2_gate_deuterons = "cebraTime2_toScint > -731 && cebraTime2_toScint < -717"    #Deuterons_TCut2
    cebraTime3_gate_deuterons = "cebraTime3_toScint > -730 && cebraTime3_toScint < -718"    #Deuterons_TCut3
    cebraTime4_gate_deuterons = "cebraTime4_toScint > -700 && cebraTime4_toScint < -688"    #Deuterons_TCut4
    cebraTimeGate_deuterons = [cebraTime0_gate_deuterons,cebraTime1_gate_deuterons,cebraTime2_gate_deuterons,cebraTime3_gate_deuterons,cebraTime4_gate_deuterons]
    
######################################################################################################################################################
    
    #initilizing name arrays
    det_i_Deuteron_Frame = ["det_0_Deuteron_Frame","det_1_Deuteron_Frame","det_2_Deuteron_Frame","det_3_Deuteron_Frame","det_4_Deuteron_Frame"]
    
    cebraTime_toScint_Deuterons_TCut = ["cebraTime0_toScint_Deuterons_TCut0","cebraTime1_toScint_Deuterons_TCut1","cebraTime2_toScint_Deuterons_TCut2","cebraTime3_toScint_Deuterons_TCut3","cebraTime4_toScint_Deuterons_TCut4"]
    
    cebraE_Deuterons_TCut = ["cebraE0_Deuterons_TCut0","cebraE1_Deuterons_TCut1","cebraE2_Deuterons_TCut2","cebraE3_Deuterons_TCut3","cebraE4_Deuterons_TCut4"]
    
    x1_Deuterons_TCut = ["x1_Deuterons_TCut0","x1_Deuterons_TCut1","x1_Deuterons_TCut2","x1_Deuterons_TCut3","x1_Deuterons_TCut4"]
    
    xavg_Deuterons_TCut = ["xavg_Deuterons_TCut0","xavg_Deuterons_TCut1","xavg_Deuterons_TCut2","xavg_Deuterons_TCut3","xavg_Deuterons_TCut4"]
    
    x1_cebraE_Deuterons_TCut = ["x1_cebraE0_Deuterons_TCut0","x1_cebraE1_Deuterons_TCut1","x1_cebraE2_Deuterons_TCut2","x1_cebraE3_Deuterons_TCut3","x1_cebraE4_Deuterons_TCut4"]
    
    xavg_cebraE_Deuterons_TCut = ["xavg_cebraE0_Deuterons_TCut0","xavg_cebraE1_Deuterons_TCut1","xavg_cebraE2_Deuterons_TCut2","xavg_cebraE3_Deuterons_TCut3","xavg_cebraE4_Deuterons_TCut4"]
    

    for i in range(n_det):
        
        #new frame for each detector with the particles in the cebraTime to Scintillator peak of intrest
        det_i_Deuteron_Frame[i] = Deuteron_Frame.Filter(cebraTimeGate_deuterons[i])
        
        #CeBrA time to scintillator for the desired particle frame
        cebraTime_toScint_Deuterons_TCut[i] = det_i_Deuteron_Frame[i].Histo1D((f"cebraTime{i}_toScint_Deuterons_TCut{i}",f"cebraTime{i}_toScint_Deuterons_TCut{i}",6400,-3200,3200),f"cebraTime{i}_toScint")
        
        #CeBrA Energy plots with the particle gate and the cebraTime to Scintillator cut
        cebraE_Deuterons_TCut[i] = det_i_Deuteron_Frame[i].Histo1D((f"cebraE{i}_Deuterons_TCut",f"cebraE{i}_Deuterons_TCut",4096,0,4096),f"cebraE{i}_GM")
        
        #focal plane x1 with the particle group nd the cebraTime to Scintillator cut
        x1_Deuterons_TCut[i] = det_i_Deuteron_Frame[i].Histo1D((f"x1_Deuterons_TCut{i}",f"x1_Deuterons_TCut{i}",600,-300,300),"x1")
        
        #focal plane x1 with the particle group nd the cebraTime to Scintillator cut
        xavg_Deuterons_TCut[i] = det_i_Deuteron_Frame[i].Histo1D((f"xavg_Deuterons_TCut{i}",f"xavg_Deuterons_TCut{i}",600,-300,300),"xavg")
        
        #Particle-Gamma Concidence Matricies
        x1_cebraE_Deuterons_TCut[i] = det_i_Deuteron_Frame[i].Histo2D((f"x1_cebraE{i}_Deuterons_TCut{i}",f"x1_cebraE{i}_Deuterons_TCut{i}",600,-300,300,4096,0,4096), "x1",f"cebraE{i}_GM")
        
        xavg_cebraE_Deuterons_TCut[i] = det_i_Deuteron_Frame[i].Histo2D((f"xavg_cebraE{i}_Deuterons_TCut{i}",f"xavg_cebraE{i}_Deuterons_TCut{i}",600,-300,300,4096,0,4096), "xavg",f"cebraE{i}_GM")
        
        
    
    #Maybe a cleaner way to do this?
    ######################################################################################################################################################
    #Deuterons
    #combines the gamme detectors to populate focal plane x1 
    
    Array_x1_Deuterons_TCut = ROOT.TH1D('Array_x1_Deuterons_TCut','Array_x1_Deuterons_TCut',600,-300,300)

    x1_Deuterons_TCut_Array_Initial = det_i_Deuteron_Frame[0].Fill(Array_x1_Deuterons_TCut,["x1"])
    x1_Deuterons_TCut1 = det_i_Deuteron_Frame[1].Fill(Array_x1_Deuterons_TCut,["x1"])
    x1_Deuterons_TCut2 = det_i_Deuteron_Frame[2].Fill(Array_x1_Deuterons_TCut,["x1"])
    x1_Deuterons_TCut3 = det_i_Deuteron_Frame[3].Fill(Array_x1_Deuterons_TCut,["x1"])
    x1_Deuterons_TCut4 = det_i_Deuteron_Frame[4].Fill(Array_x1_Deuterons_TCut,["x1"])

    x1_Deuterons_TCut_Array_Initial.Add(x1_Deuterons_TCut1.GetValue())
    x1_Deuterons_TCut_Array_Initial.Add(x1_Deuterons_TCut2.GetValue())
    x1_Deuterons_TCut_Array_Initial.Add(x1_Deuterons_TCut3.GetValue())
    x1_Deuterons_TCut_Array_Initial.Add(x1_Deuterons_TCut4.GetValue())
    
    #combines the gamme detectors to populate the gamma spectrum
    Array_cebraE_Deuterons_TCut = ROOT.TH1D('Array_cebraE_Deuterons_TCut','Array_cebraE_Deuterons_TCut',4096,0,4096)

    cebraE_Deuterons_TCut_Array_Initial = det_i_Deuteron_Frame[0].Fill(Array_cebraE_Deuterons_TCut,["cebraE0_GM"])
    cebraE_Deuterons_TCut1 = det_i_Deuteron_Frame[1].Fill(Array_cebraE_Deuterons_TCut,["cebraE1_GM"])
    cebraE_Deuterons_TCut2 = det_i_Deuteron_Frame[2].Fill(Array_cebraE_Deuterons_TCut,["cebraE2_GM"])
    cebraE_Deuterons_TCut3 = det_i_Deuteron_Frame[3].Fill(Array_cebraE_Deuterons_TCut,["cebraE3_GM"])
    cebraE_Deuterons_TCut4 = det_i_Deuteron_Frame[4].Fill(Array_cebraE_Deuterons_TCut,["cebraE4_GM"])
    
    cebraE_Deuterons_TCut_Array_Initial.Add(cebraE_Deuterons_TCut1.GetValue())
    cebraE_Deuterons_TCut_Array_Initial.Add(cebraE_Deuterons_TCut2.GetValue())
    cebraE_Deuterons_TCut_Array_Initial.Add(cebraE_Deuterons_TCut3.GetValue())
    cebraE_Deuterons_TCut_Array_Initial.Add(cebraE_Deuterons_TCut4.GetValue())

    
    Array_CeBrA_x1_Deuterons_TCut = ROOT.TH2D('Array_x1_CeBrA_Deuterons_TCut','Array_x1_CeBrA_Deuterons_TCut',600,-300,300,4096,0,4096)

    x1_cebraE_Deuterons_TCut_Array_Initial = det_i_Deuteron_Frame[0].Fill(Array_CeBrA_x1_Deuterons_TCut,["x1","cebraE0_GM"])
    x1_cebraE_Deuterons_TCut1 = det_i_Deuteron_Frame[1].Fill(Array_CeBrA_x1_Deuterons_TCut,["x1","cebraE1_GM"])
    x1_cebraE_Deuterons_TCut2 = det_i_Deuteron_Frame[2].Fill(Array_CeBrA_x1_Deuterons_TCut,["x1","cebraE2_GM"])
    x1_cebraE_Deuterons_TCut3 = det_i_Deuteron_Frame[3].Fill(Array_CeBrA_x1_Deuterons_TCut,["x1","cebraE3_GM"])
    x1_cebraE_Deuterons_TCut4 = det_i_Deuteron_Frame[4].Fill(Array_CeBrA_x1_Deuterons_TCut,["x1","cebraE4_GM"])

    x1_cebraE_Deuterons_TCut_Array_Initial.Add(x1_cebraE_Deuterons_TCut1.GetValue())
    x1_cebraE_Deuterons_TCut_Array_Initial.Add(x1_cebraE_Deuterons_TCut2.GetValue())
    x1_cebraE_Deuterons_TCut_Array_Initial.Add(x1_cebraE_Deuterons_TCut3.GetValue())
    x1_cebraE_Deuterons_TCut_Array_Initial.Add(x1_cebraE_Deuterons_TCut4.GetValue())

    
    
    
    
    # c1 = ROOT.TCanvas()

    output.cd()
    
    PID.Write()
    PID_Deuterons.Write()

    xavg_Deuterons.Write()
    x1_only1plane_Deuterons.Write()
    x1_bothplanes_Deuterons.Write()
    x2_bothplanes_Deuterons.Write()

    for i in range(n_det):
        cebraE_noCuts[i].Write()
        cebraE_GM_noCuts[i].Write()
        
        cebraTime_toScint_noCuts[i].Write()
        
        cebraTime_toScint_Deuterons[i].Write()
        
        cebraTime_toScint_Deuterons_TCut[i].Write()
        
        cebraE_Deuterons_TCut[i].Write()
        
        x1_Deuterons_TCut[i].Write()
        
        xavg_Deuterons_TCut[i].Write()
        
        x1_cebraE_Deuterons_TCut[i].Write()
        
        xavg_cebraE_Deuterons_TCut[i].Write()
        
    x1_Deuterons_TCut_Array_Initial.Write()
    cebraE_Deuterons_TCut_Array_Initial.Write()
    x1_cebraE_Deuterons_TCut_Array_Initial.Write()

    output.Close()

    end = timer()
    print("Total Elapsed time for current run: ",end - start)
    # input("Hit Enter when done")


###############
#figure out how to sum the 2d plots in order to get a better energy calibration for the projections
#do the GS Band Code / and band by declaring two points and the width of the band