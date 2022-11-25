from cmath import pi
import ROOT
import math
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

# C12_RunList = [82, 83, 84 ,85 ,86, 87, 88 ,89, 90,\
#          91, 92, 93, 94, 95, 96, 106, 107, 108,\
#              110, 111, 112, 113, 114, 115, 116,\
#                  117, 119, 120, 121, 122, 123]

# C12_RunList = [82]

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
    # PATH = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingD4ir/12C_analyzed/"
    PATH = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/"
    
    OUT = f"./49Ti/run_{RUNNO}_histo.root"
    output = ROOT.TFile.Open(OUT, "RECREATE")

    #number of detectors
    n_det = 5
    

    #PLOTS SHOULD BE DRAWN AT THE END
    #optimize the number of thread before you run the loop
    ROOT.ROOT.EnableImplicitMT(5)
    treeName = "SPSTree"
    
######################################################################################################################################################
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

######################################################################################################################################################

    #Gets the Particle Identification
    PID = Raw_Frame.Filter("scintLeftTime != -1 && anodeBackTime != -1").Histo2D(("ScintL_AnodeBack","PID",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
    # c_PID = ROOT.TCanvas()
    # PID.Draw("colz") 

######################################################################################################################################################
    #Select the protons for d,p reaction
    
    #12C, run 82 through 123
    # Protons = "scintLeft > 300 && scintLeft < 900 && anodeBack > 100 && anodeBack < 800"
    
    #49Ti, run 184 through 232
    Protons = "scintLeft > 800 && scintLeft < 1600 && anodeBack > 50 && anodeBack < 800"
    Deuterons = "scintLeft > 200 && scintLeft < 550 && anodeBack > 300 && anodeBack < 1100"
    
    PID_Protons = Raw_Frame.Filter(Protons).Histo2D(("ScintL_AnodeBack_Protons","Protons",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
    PID_Deuterons = Raw_Frame.Filter(Deuterons).Histo2D(("ScintL_AnodeBack_Deuterons","Deuterons",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
        
######################################################################################################################################################

    cebraE_noCuts = []
    cebraE_GM_noCuts = []
    cebraTime_toScint_noCuts = []
    for i in range(n_det):
        
        cebraE_noCuts.append(f"cebraE{i}_noCuts")
        cebraE_noCuts[i] = Raw_Frame.Filter(f"cebraE{i}>-1").Histo1D((f"cebraE{i}_noCuts",f"cebraE{i}_noCuts",4096,0,4096),f"cebraE{i}")
        
        cebraE_GM_noCuts.append(f"cebraE{i}_GM_noCuts")
        cebraE_GM_noCuts[i] = Raw_Frame.Filter(f"cebraE{i}>-1").Histo1D((f"cebraE{i}_GM_noCuts",f"cebraE{i}_GM_noCuts",1024,0,4096),f"cebraE{i}_GM")
        
        cebraTime_toScint_noCuts.append(f"cebraTime{i}_toScint_noCuts")
        cebraTime_toScint_noCuts[i] = Raw_Frame.Filter(f"scintLeftTime != -1  && anodeBackTime != -1 && cebraE{i} != -1").Histo1D((f"cebraTime{i}_toScint_noCuts",f"cebraTime{i}_toScint_noCuts",6400,-3200,3200),f"cebraTime{i}_toScint")
        
######################################################################################################################################################

    #ax^2+bx+c
    x1_cal_a = -0.003824652423017848
    x1_cal_b = -21.158024905122485
    x1_cal_c = 2157.116872413482
    
    cebraE_cal_a = -2.221225e-5
    cebraE_cal_b = 1.6342829
    cebraE_cal_c = -63.33630
    
    # 1553.768	1003.36
    # 2205.722	1414.85
    # 2618.6	1680.18
    # 3271	    2100

    #Creating a new frames for different particle groups
    Proton_Frame = Raw_Frame.Filter(Protons).Define("x1_ECal", f"{x1_cal_a}*x1*x1 + {x1_cal_b}*x1 + {x1_cal_c}")\
                            .Define("cebraE0_ECal", f"{cebraE_cal_a}*cebraE0_GM*cebraE0_GM + {cebraE_cal_b}*cebraE0_GM + {cebraE_cal_c}")\
                            .Define("cebraE1_ECal", f"{cebraE_cal_a}*cebraE1_GM*cebraE1_GM + {cebraE_cal_b}*cebraE1_GM + {cebraE_cal_c}")\
                            .Define("cebraE2_ECal", f"{cebraE_cal_a}*cebraE2_GM*cebraE2_GM + {cebraE_cal_b}*cebraE2_GM + {cebraE_cal_c}")\
                            .Define("cebraE3_ECal", f"{cebraE_cal_a}*cebraE3_GM*cebraE3_GM + {cebraE_cal_b}*cebraE3_GM + {cebraE_cal_c}")\
                            .Define("cebraE4_ECal", f"{cebraE_cal_a}*cebraE4_GM*cebraE4_GM + {cebraE_cal_b}*cebraE4_GM + {cebraE_cal_c}")
                                                                      
    # Deuteron_Frame = Raw_Frame.Filter(Deuterons)
    
    xavg_Protons = Proton_Frame.Histo1D(("xavg_Protons","xavg_Protons",600,-300,300),"xavg")
    x1_bothplanes_Protons = Proton_Frame.Filter("x1 != -1e6 && x2 != -1e6").Histo1D(("x1_bothplanes_Protons","x1_bothplanes_Protons",600,-300,300),"x1")
    x2_bothplanes_Protons = Proton_Frame.Filter("x1 != -1e6 && x2 != -1e6").Histo1D(("x2_bothplanes_Protons","x2_bothplanes_Protons",600,-300,300),"x2")
    x1_only1plane_Protons = Proton_Frame.Filter("x1 != -1e6 && x2 == -1e6").Histo1D(("x1_only1plane_Protons","x1_only1plane_Protons",600,-300,300),"x1")
    x1_theta_Protons = Proton_Frame.Histo2D(("x1_theta_Protons","x1_theta_Protons",600,-300,300,100,0,pi/2),"x1","theta")
    xavg_theta_Protons = Proton_Frame.Histo2D(("xavg_theta_Protons","xavg_theta_Protons",600,-300,300,100,0,pi/2),"xavg","theta")
    
    x1_only1plane_Protons_ECal = Proton_Frame.Filter("x1 != -1e6 && x2 == -1e6").Histo1D(("x1_only1plane_Protons_ECal","x1_only1plane_Protons_ECal",512,-100,8092),"x1_ECal")
    x1_bothplane_Protons_ECal = Proton_Frame.Filter("x1 != -1e6 && x2 != -1e6").Histo1D(("x1_bothplanes_Protons_ECal","x1_bothplanes_Protons_ECal",512,-100,8092),"x1_ECal")
    x1_ECal_theta_Protons = Proton_Frame.Histo2D(("x1_ECal_theta_Protons","x1_ECal_theta_Protons",512,-100,8092,600,0,pi/2),"x1_ECal","theta")
    

    
    cebraTime_toScint_Protons = []
    cebraTime_toScint_Deuterons = []
    for i in range(n_det):
        
        cebraTime_toScint_Protons.append(f"cebraTime{i}_toScint_Protons")
        cebraTime_toScint_Protons[i] = Proton_Frame.Histo1D((f"cebraTime{i}_toScint_Protons",f"cebraTime{i}_toScint_Protons",6400,-3200,3200),f"cebraTime{i}_toScint")
        
        # cebraTime_toScint_Deuterons.append(f"cebraTime{i}_toScint_Deuterons")
        # cebraTime_toScint_Deuterons[i] = Deuteron_Frame.Histo1D((f"cebraTime{i}_toScint_Deuterons",f"cebraTime{i}_toScint_Deuterons",6400,-3200,3200),f"cebraTime{i}_toScint")


    #12C
    # cebraTime0_gate_protons = "cebraTime0_toScint > -692 && cebraTime0_toScint < -684"
    # cebraTime1_gate_protons = "cebraTime1_toScint > -689 && cebraTime1_toScint < -681"
    # cebraTime2_gate_protons = "cebraTime2_toScint > -691 && cebraTime2_toScint < -682"
    # cebraTime3_gate_protons = "cebraTime3_toScint > -690 && cebraTime3_toScint < -682"
    # cebraTime4_gate_protons = "cebraTime4_toScint > -660 && cebraTime4_toScint < -652"
    
    #49Ti
    cebraTime0_gate_protons = "cebraTime0_toScint > -678 && cebraTime0_toScint < -670"     #Protons_TCut0
    cebraTime1_gate_protons = "cebraTime1_toScint > -676 && cebraTime1_toScint < -669"     #Protons_TCut1
    cebraTime2_gate_protons = "cebraTime2_toScint > -676 && cebraTime2_toScint < -669"     #Protons_TCut2
    cebraTime3_gate_protons = "cebraTime3_toScint > -676 && cebraTime3_toScint < -668"     #Protons_TCut3
    cebraTime4_gate_protons = "cebraTime4_toScint > -645 && cebraTime4_toScint < -638"     #Protons_TCut4
    cebraTimeGate_protons = [cebraTime0_gate_protons,cebraTime1_gate_protons,cebraTime2_gate_protons,cebraTime3_gate_protons,cebraTime4_gate_protons]
    
    # cebraTime0_gate_deuterons = "cebraTime0_toScint > -733 && cebraTime0_toScint < -720"    #Deuterons_TCut0
    # cebraTime1_gate_deuterons = "cebraTime1_toScint > -733 && cebraTime1_toScint < -717"    #Deuterons_TCut1
    # cebraTime2_gate_deuterons = "cebraTime2_toScint > -731 && cebraTime2_toScint < -717"    #Deuterons_TCut2
    # cebraTime3_gate_deuterons = "cebraTime3_toScint > -730 && cebraTime3_toScint < -718"    #Deuterons_TCut3
    # cebraTime4_gate_deuterons = "cebraTime4_toScint > -700 && cebraTime4_toScint < -688"    #Deuterons_TCut4
    # cebraTimeGate_deuterons = [cebraTime0_gate_deuterons,cebraTime1_gate_deuterons,cebraTime2_gate_deuterons,cebraTime3_gate_deuterons,cebraTime4_gate_deuterons]
    
######################################################################################################################################################

    det_i_Proton_Frame = []
    cebraTime_toScint_Protons_TCut = []
    cebraE_Protons_TCut = []
    x1_Protons_TCut = []
    x1_cebraE_Protons_TCut = []
    xavg_Protons_TCut = []
    xavg_cebraE_Protons_TCut = []
    
    #Energy Calibrated
    x1_ECal_Protons_TCut = []
    cebraE_ECal_Protons_TCut = []
    x1_cebraE_ECal_Protons_TCut = []
    
    #Ground State Band 
    GS_Band_Proton_Frame = []
    x1_GS_Band_Protons_TCut = []
    cebraE_GS_Band_Protons_TCut = []
    x1_cebraE_GS_Band_Protons_TCut = []
    
    #Excited States Energies
    # x1_Energy_Width = 35    
    #Calibration Energy Peaks
    # x1_Energies = [0.0,1553.794,2674.932,4147.21,4789.97,4880.705,5186.103,5379.942,5547.81,5837.2,6136.3]
    
    det_i_Excited_State_Frame = []
    
    x1_GS_Gate = "-100 < x1_ECal && x1_ECal < 100"
    x1_1553_Gate = "1400 < x1_ECal && x1_ECal < 1700"
    x1_2674_Gate = "2520 < x1_ECal && x1_ECal < 2760"
    x1_4100_doublet_Gate = "4050 < x1_ECal && x1_ECal < 4275"
    x1_4500_doublet_Gate = "4450 < x1_ECal && x1_ECal < 4650"
    x1_4880_Gate = "4720 < x1_ECal && x1_ECal < 5000"
    x1_5186_Gate = "5120 < x1_ECal && x1_ECal < 5240"
    x1_5379_Gate = "5322 < x1_ECal && x1_ECal < 5460"
    x1_5500_5630_Gate = "5500 < x1_ECal && x1_ECal < 5630"
    x1_5660_5780_Gate = "5500 < x1_ECal && x1_ECal < 5780"
    x1_5837_Gate = "5780 < x1_ECal && x1_ECal < 5880"
    x1_5946_Gate = "5890 < x1_ECal && x1_ECal < 6020"
    x1_6020_6100_Gate = "6020 < x1_ECal && x1_ECal < 6100"
    x1_6123_Gate = "6100 < x1_ECal && x1_ECal < 6190"
    
    x1_E_Gates = [["0",x1_GS_Gate], ["1553",x1_1553_Gate], ["2674",x1_2674_Gate], ["4100_doublet", x1_4100_doublet_Gate],\
        ["4500_doublet",x1_4500_doublet_Gate], ["4880",x1_4880_Gate], ["5186",x1_5186_Gate], ["5379",x1_5379_Gate],\
            ["5500_5630",x1_5500_5630_Gate], ["5660_5780",x1_5660_5780_Gate], ["5837",x1_5837_Gate], ["5946",x1_5946_Gate],\
                ["6020_6100",x1_6020_6100_Gate], ["6123",x1_6123_Gate]]
        
    x1_E_List_Protons_TCut = []
    cebraE_E_List_Protons_TCut = []
    
    for i in range(n_det):
                        
        det_i_Proton_Frame.append(f"det_{i}_Proton_Frame")
        det_i_Proton_Frame[i] = Proton_Frame.Filter(cebraTimeGate_protons[i])
        
        cebraTime_toScint_Protons_TCut.append(f"cebraTime{i}_toScint_Protons_TCut{i}")
        cebraTime_toScint_Protons_TCut[i] = det_i_Proton_Frame[i].Histo1D((f"cebraTime{i}_toScint_Protons_TCut{i}",f"cebraTime{i}_toScint_Protons_TCut{i}",6400,-3200,3200),f"cebraTime{i}_toScint")
        
        cebraE_Protons_TCut.append(f"cebraE{i}_Protons_TCut{i}")
        cebraE_Protons_TCut[i] = det_i_Proton_Frame[i].Histo1D((f"cebraE{i}_Protons_TCut",f"cebraE{i}_Protons_TCut",4096,0,4096),f"cebraE{i}_GM")
        
        x1_Protons_TCut.append(f"x1_Protons_TCut{i}")
        x1_Protons_TCut[i] = det_i_Proton_Frame[i].Histo1D((f"x1_Protons_TCut{i}",f"x1_Protons_TCut{i}",600,-300,300),"x1")
        
        x1_cebraE_Protons_TCut.append(f"x1_cebraE{i}_Protons_TCut{i}")
        x1_cebraE_Protons_TCut[i] = det_i_Proton_Frame[i].Histo2D((f"x1_cebraE{i}_Protons_TCut{i}",f"x1_cebraE{i}_Protons_TCut{i}",600,-300,300,4096,0,4096), "x1",f"cebraE{i}_GM")
        
        xavg_Protons_TCut.append(f"xavg_Protons_TCut{i}")
        xavg_Protons_TCut[i] = det_i_Proton_Frame[i].Histo1D((f"xavg_Protons_TCut{i}",f"xavg_Protons_TCut{i}",600,-300,300),"xavg")
        
        xavg_cebraE_Protons_TCut.append(f"xavg_cebraE{i}_Protons_TCut{i}")
        xavg_cebraE_Protons_TCut[i] = det_i_Proton_Frame[i].Histo2D((f"xavg_cebraE{i}_Protons_TCut{i}",f"xavg_cebraE{i}_Protons_TCut{i}",600,-300,300,4096,0,4096), "xavg",f"cebraE{i}_GM")
        
        #Energy Calibrated
        x1_ECal_Protons_TCut.append("x1_ECal_Protons_TCut{i}")
        x1_ECal_Protons_TCut[i] = det_i_Proton_Frame[i].Histo1D((f"x1_ECal_Protons_TCut{i}",f"x1_ECal_Protons_TCut{i}",512,-100,8092),"x1_ECal")
        
        cebraE_ECal_Protons_TCut.append(f"cebraE{i}_ECal_Protons_TCut{i}")
        cebraE_ECal_Protons_TCut[i] = det_i_Proton_Frame[i].Histo1D((f"cebraE{i}_ECal_Protons_TCut",f"cebraE{i}_ECal_Protons_TCut",680,0,6800),f"cebraE{i}_ECal")
        
        x1_cebraE_ECal_Protons_TCut.append(f"x1_cebraE{i}_ECal_Protons_TCut{i}")
        x1_cebraE_ECal_Protons_TCut[i] = det_i_Proton_Frame[i].Histo2D((f"x1_cebraE{i}_ECal_Protons_TCut{i}",f"x1_cebraE{i}_ECal_Protons_TCut{i}",512,-100,8092,1024,0,8192), "x1_ECal",f"cebraE{i}_ECal")
        

        Excited_State_Frame = []
        det_i_x1 = []
        det_i_cebraE = []
        for j in range(0,len(x1_E_Gates)):
            
            Excited_State_Frame.append(f"{x1_E_Gates[j][0]}_Frame")
            Excited_State_Frame[j] = det_i_Proton_Frame[i].Filter(x1_E_Gates[j][1])
            
            det_i_x1.append(f"x1_{x1_E_Gates[j][0]}_Protons_TCut{i}")
            # det_i_x1[j] = det_i_Proton_Frame[i].Filter(x1_E_Gates[j][1]).Histo1D((f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",512,-100,8092),"x1_ECal")
            det_i_x1[j] = Excited_State_Frame[j].Histo1D((f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",512,-100,8092),"x1_ECal")
            
            det_i_cebraE.append(f"cebraE_{x1_E_Gates[j][0]}_Protons_TCut{i}")
            # det_i_cebraE[j] = det_i_Proton_Frame[i].Filter(x1_E_Gates[j][1]).Histo1D((f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",512,-100,8092),f"cebraE{i}_ECal")
            det_i_cebraE[j] = Excited_State_Frame[j].Histo1D((f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",1024,0,4096),f"cebraE{i}_GM")
            
            
        x1_E_List_Protons_TCut.append(det_i_x1)
        cebraE_E_List_Protons_TCut.append(det_i_cebraE)
        det_i_Excited_State_Frame.append(Excited_State_Frame)

        
        #Ground State Band (also Energy Calibrated)
        
        #Filter condition
        GS_Band = f"( cebraE{i}_ECal <= (x1_ECal +100) ) && ( cebraE{i}_ECal >= (x1_ECal-100) )"
        
        # Ground State Band Data Frame
        GS_Band_Proton_Frame.append(f"det_{i}_GS_Band_Proton_Frame")
        GS_Band_Proton_Frame[i] = det_i_Proton_Frame[i].Filter(GS_Band)
        
        x1_GS_Band_Protons_TCut.append("x1_GS_Band_Protons_TCut{i}")
        x1_GS_Band_Protons_TCut[i] = GS_Band_Proton_Frame[i].Histo1D((f"x1_GS_Band_Protons_TCut{i}",f"x1_GS_Band_Protons_TCut{i}",512,-100,8092),"x1_ECal")
        
        cebraE_GS_Band_Protons_TCut.append(f"cebraE{i}_GS_Band_Protons_TCut{i}")
        cebraE_GS_Band_Protons_TCut[i] = GS_Band_Proton_Frame[i].Histo1D((f"cebraE{i}_GS_Band_Protons_TCut",f"cebraE{i}_GS_Band_Protons_TCut",1024,0,8192),f"cebraE{i}_ECal")
        
        x1_cebraE_GS_Band_Protons_TCut.append(f"x1_cebraE{i}_GS_Band_Protons_TCut{i}")
        x1_cebraE_GS_Band_Protons_TCut[i] = GS_Band_Proton_Frame[i].Histo2D((f"x1_cebraE{i}_GS_Band_Protons_TCut{i}",f"x1_cebraE{i}_GS_Band_Protons_TCut{i}",512,-100,8092,1024,0,8192), "x1_ECal",f"cebraE{i}_ECal")
        
        
        if i==0:
            #initilizing CeBrA Plots
            
            #Not Calibrated
            Array_x1_Protons_TCut = det_i_Proton_Frame[i].Histo1D(("Array_x1_Protons_TCut","Array_x1_Protons_TCut",600,-300,300),"x1")
            Array_cebraE_Protons_TCut = det_i_Proton_Frame[i].Histo1D(("Array_cebraE_Protons_TCut","Array_cebraE_Protons_TCut",4096,0,4096),f"cebraE{i}_GM")
            Array_x1_cebraE_Protons_TCut = det_i_Proton_Frame[i].Histo2D(("Array_x1_cebraE_Protons_TCut","Array_x1_cebraE_Protons_TCut",600,-300,300,4096,0,4096), "x1",f"cebraE{i}_GM")
            
            #Energy Calibrated
            Array_x1_ECal_Protons = det_i_Proton_Frame[i].Histo1D(("Array_x1_ECal_Protons","Array_x1_ECal_Protons",512,-100,8092),"x1_ECal")
            Array_cebraE_ECal_Protons = det_i_Proton_Frame[i].Histo1D(("Array_cebraE_ECal_Protons","Array_cebraE_ECal_Protons",680,0,6800),f"cebraE{i}_ECal")
            Array_x1_cebraE_ECal_Protons = det_i_Proton_Frame[i].Histo2D(("x1_cebraE_ECal_Array","x1_cebraE_ECal_Array",512,-100,8092,1024,0,8192), "x1_ECal",f"cebraE{i}_ECal")
            
            #GS Band
            Array_x1_GS_Band_Protons = GS_Band_Proton_Frame[i].Histo1D(("Array_x1_GS_Band_Protons","Array_x1_GS_Band_Protons_TCut",512,-100,8092),"x1_ECal")
            Array_cebraE_GS_Band_Protons = GS_Band_Proton_Frame[i].Histo1D(("Array_cebraE_GS_Band_Protons","Array_cebraE_GS_Band_Protons",1024,0,8192),f"cebraE{i}_ECal")
            Array_x1_cebraE_GS_Band_Protons = GS_Band_Proton_Frame[i].Histo2D(("x1_cebraE_GS_Band_Array","x1_cebraE_GS_Band_Array",512,-100,8092,1024,0,8192), "x1_ECal",f"cebraE{i}_ECal")
            
            #Excited States
            Array_x1_E_List = []
            Array_cebraE_E_List = []
            for j in range(0,len(x1_E_Gates)):
                
                                                
                Array_x1_E_List.append(f"Array_{x1_E_Gates[j][0]}_keV_x1")
                # Array_x1_E_List[j] = det_i_Proton_Frame[i].Filter(x1_E_Gates[j][1]).Histo1D((f"Array_{x1_E_Gates[j][0]}_keV_x1",f"Array_{x1_E_Gates[j][0]}_keV_x1",512,-100,8092),"x1_ECal")
                Array_x1_E_List[j] = det_i_Excited_State_Frame[i][j].Histo1D((f"Array_{x1_E_Gates[j][0]}_keV_x1",f"Array_{x1_E_Gates[j][0]}_keV_x1",512,-100,8092),"x1_ECal")


                Array_cebraE_E_List.append(f"Array_{x1_E_Gates[j][0]}_keV_cebraE")
                # Array_cebraE_E_List[j] = det_i_Proton_Frame[i].Filter(x1_E_Gates[j][1]).Histo1D((f"Array_{x1_E_Gates[j][0]}_keV_cebraE",f"Array_{x1_E_Gates[j][0]}_keV_cebraE",512,-100,8092),f"cebraE{i}_ECal")
                Array_cebraE_E_List[j] = det_i_Excited_State_Frame[i][j].Histo1D((f"Array_{x1_E_Gates[j][0]}_keV_cebraE",f"Array_{x1_E_Gates[j][0]}_keV_cebraE",1024,0,4096),f"cebraE{i}_GM")
                
            #     # print(Array_x1_E)
            
            
        if i!=0:
            #Not Calibrated
            Array_x1_Protons_TCut.Add(x1_Protons_TCut[i].GetPtr())
            Array_cebraE_Protons_TCut.Add(cebraE_Protons_TCut[i].GetPtr())
            Array_x1_cebraE_Protons_TCut.Add(x1_cebraE_Protons_TCut[i].GetPtr())
            
            #Energy Calibrated
            Array_x1_ECal_Protons.Add(x1_ECal_Protons_TCut[i].GetPtr())
            Array_cebraE_ECal_Protons.Add(cebraE_ECal_Protons_TCut[i].GetPtr())
            Array_x1_cebraE_ECal_Protons.Add(x1_cebraE_ECal_Protons_TCut[i].GetPtr())
            
            #GS Band
            Array_x1_GS_Band_Protons.Add(x1_GS_Band_Protons_TCut[i].GetPtr())
            Array_cebraE_GS_Band_Protons.Add(cebraE_GS_Band_Protons_TCut[i].GetPtr())
            Array_x1_cebraE_GS_Band_Protons.Add(x1_cebraE_GS_Band_Protons_TCut[i].GetPtr())
            
            for j in range(0,len(x1_E_Gates)):
                Array_x1_E_List[j].Add(x1_E_List_Protons_TCut[i][j].GetPtr())
                Array_cebraE_E_List[j].Add(cebraE_E_List_Protons_TCut[i][j].GetPtr())
                
            
    # print(x1_E_List_Protons_TCut)
    
    # c1 = ROOT.TCanvas()

    output.cd()
    
    output.mkdir("Common")
    
    output.mkdir("49Ti_dpg")
    output.mkdir("49Ti_dpg/Time_Correlation/")
    output.mkdir("49Ti_dpg/Particle_Gamma/")
    output.mkdir("49Ti_dpg/Individual_Spectrums/")
    output.mkdir("49Ti_dpg/Energy_Calibrated/")
    output.mkdir("49Ti_dpg/GS_Band/")
    output.mkdir("49Ti_dpg/50Ti_Excited_States/")
    
    
    output.cd("Common")
    PID.Write()
    
    output.cd("49Ti_dpg")
    PID_Protons.Write()
    xavg_Protons.Write()
    x1_only1plane_Protons.Write()
    x1_bothplanes_Protons.Write()
    x2_bothplanes_Protons.Write()
    x1_theta_Protons.Write()
    xavg_theta_Protons.Write()
    
    output.cd("49Ti_dpg/Energy_Calibrated/")
    x1_only1plane_Protons_ECal.Write()
    x1_bothplane_Protons_ECal.Write()
    x1_ECal_theta_Protons.Write()
    
    for i in range(n_det):
        
        output.cd("Common")
        cebraE_noCuts[i].Write()
        cebraE_GM_noCuts[i].Write()
        cebraTime_toScint_noCuts[i].Write()
        
        output.cd("49Ti_dpg/Time_Correlation/")
        cebraTime_toScint_Protons[i].Write()
        cebraTime_toScint_Protons_TCut[i].Write()
        
        output.cd("49Ti_dpg/Individual_Spectrums/")
        cebraE_Protons_TCut[i].Write()
        x1_Protons_TCut[i].Write()
        xavg_Protons_TCut[i].Write()
        
        output.cd("49Ti_dpg/Particle_Gamma/")
        x1_cebraE_Protons_TCut[i].Write()
        xavg_cebraE_Protons_TCut[i].Write()
        

    # Summed Array Plots
    output.cd("49Ti_dpg/")
    Array_x1_Protons_TCut.Write()
    Array_cebraE_Protons_TCut.Write()
    Array_x1_cebraE_Protons_TCut.Write() 
    
    output.cd("49Ti_dpg/50Ti_Excited_States/")
    for j in range(0,len(x1_E_Gates)):
        Array_x1_E_List[j].Write()
        Array_cebraE_E_List[j].Write()
        
    #Summed Array Energy Calibrated Plots
    output.cd("49Ti_dpg/Energy_Calibrated/")
    Array_x1_ECal_Protons.Write()
    Array_cebraE_ECal_Protons.Write()
    Array_x1_cebraE_ECal_Protons.Write()
    
    #Summed Array GS Band Plots
    output.cd("49Ti_dpg/GS_Band/")
    Array_x1_GS_Band_Protons.Write()
    Array_cebraE_GS_Band_Protons.Write()
    Array_x1_cebraE_GS_Band_Protons.Write()

    output.Close()

    end = timer()
    print("Total Elapsed time for current run: ",end - start)
    # input("Hit Enter when done")