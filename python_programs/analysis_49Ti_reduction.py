from cmath import pi
import ROOT
import math
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

# Ti49_RunList = [184,185,186,187,188,190,191,192,193,\
#     194,195,196,197,198,199,200,201,202,203,204,205,\
#         206,207,208,209,210,211,212,213,214,215,216,\
#             217,218,219,220,221,222,223,224,225,226,\
#                 227,228,229,230,231,232]
                
Ti49_RunList=[184]

#number of detectors
n_det = 5

for RUNNO in Ti49_RunList:
# for RUNNO in C12_RunList:
        
    print("Current Run: ",RUNNO)
    start = timer()
    
    #where the analyzed files/tree is located
    PATH = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/"
    
    # OUTPATH = f"./run_{RUNNO}_common_histo.root"
    # output = ROOT.TFile.Open(OUTPATH, "RECREATE")
    
    #Do Not have it multi-thread if using snapshot 

    treeName = "SPSTree"
    
    #gain matching calibration numbers, got these with a HDTV Script
    #need to parse the HDTV calibration file
    intercept_GM = []
    slope_GM = []
    # with open("/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/12C_analyzed/run_{RUNNO}_cebra_cal.txt","r") as f:
    with open(f"/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/run_{RUNNO}_cebra_cal.txt","r") as f:
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
                    
    

    cebraE_noCuts = []
    cebraE_GM_noCuts = []
    cebraTime_toScint_noCuts = []
    for i in range(n_det):
        
        cebraE_noCuts.append(f"cebraE{i}_noCuts")
        cebraE_noCuts[i] = Raw_Frame.Filter(f"cebraE{i}>-1").Histo1D((f"cebraE{i}_noCuts",f"cebraE{i}_noCuts",4096,0,4096),f"cebraE{i}")
        
        cebraE_GM_noCuts.append(f"cebraE{i}_GM_noCuts")
        cebraE_GM_noCuts[i] = Raw_Frame.Filter(f"cebraE{i}>-1").Histo1D((f"cebraE{i}_GM_noCuts",f"cebraE{i}_GM_noCuts",4096,0,4096),f"cebraE{i}_GM")
        
        cebraTime_toScint_noCuts.append(f"cebraTime{i}_toScint_noCuts")
        cebraTime_toScint_noCuts[i] = Raw_Frame.Filter(f"scintLeftTime != -1  && anodeBackTime != -1 && cebraE{i} != -1").Histo1D((f"cebraTime{i}_toScint_noCuts",f"cebraTime{i}_toScint_noCuts",6400,-3200,3200),f"cebraTime{i}_toScint")
        

    #Gets the Particle Identification
    PID = Raw_Frame.Filter("scintLeftTime != -1 && anodeBackTime != -1").Histo2D(("ScintL_AnodeBack","PID",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
    c_PID = ROOT.TCanvas()
    PID.Draw("colz") 

######################################################################################################################################################
    #Select the protons for d,p reaction
    
    #12C, run 82 through 123
    # Protons = "scintLeft > 300 && scintLeft < 900 && anodeBack > 100 && anodeBack < 800"
    
    #49Ti, run 184 through 232
    Protons = "scintLeft > 800 && scintLeft < 1600 && anodeBack > 50 && anodeBack < 800"
    Deuterons = "scintLeft > 150 && scintLeft < 550 && anodeBack > 280 && anodeBack < 1200 && (anodeBack - 700 <= 3.558 * (scintLeft - 130)) "
    
    
    # PID_Protons = Raw_Frame.Filter(Protons).Histo2D(("ScintL_AnodeBack_Protons","Protons",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
    PID_Deuterons = Raw_Frame.Filter(Deuterons).Histo2D(("ScintL_AnodeBack_Deuterons","Deuterons",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
        
    # Proton_Frame = Raw_Frame.Filter(Protons)
    Deuteron_Frame = Raw_Frame.Filter(Deuterons)
                                                                      
    # OUT = f"./run_{RUNNO}_protons.root"
    OUT = f"./run_{RUNNO}_deuterons.root"


    # Recreating the tree with only the data that I want
    # Proton_Frame.Snapshot("SPSTree",OUT,{"anodeBack","scintLeft","scintLeftTime","cathode","xavg","x1","x2","theta",\
    #     "cebraE0","cebraE1","cebraE2","cebraE3","cebraE4","cebraTime0","cebraTime1","cebraTime2","cebraTime3","cebraTime4",\
    #         "cebraTime0_toScint","cebraTime1_toScint","cebraTime2_toScint","cebraTime3_toScint","cebraTime4_toScint",\
    #             "cebraE0_GM","cebraE1_GM","cebraE2_GM","cebraE3_GM","cebraE4_GM"})
    
    # Deuteron_Frame.Snapshot("SPSTree",OUT,{"anodeBack","scintLeft","scintLeftTime","cathode","xavg","x1","x2","theta",\
    #     "cebraE0","cebraE1","cebraE2","cebraE3","cebraE4","cebraTime0","cebraTime1","cebraTime2","cebraTime3","cebraTime4",\
    #         "cebraTime0_toScint","cebraTime1_toScint","cebraTime2_toScint","cebraTime3_toScint","cebraTime4_toScint",\
    #             "cebraE0_GM","cebraE1_GM","cebraE2_GM","cebraE3_GM","cebraE4_GM"})
    
    # output.cd()
    
    # output.mkdir("Common")
    # output.cd("Common")
    
    # PID.Write()
    # for i in range(n_det):
        
    #     cebraE_noCuts[i].Write()
    #     cebraE_GM_noCuts[i].Write()
    #     cebraTime_toScint_noCuts[i].Write()
    
    # c1 = ROOT.TCanvas()
    # PID_Deuterons.Draw("Colz")


    # output.Close()

    end = timer()
    print("Total Elapsed time for current run: ",end - start)
    input("Hit Enter when done")
    