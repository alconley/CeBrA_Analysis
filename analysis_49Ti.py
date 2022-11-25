from cmath import pi
import ROOT
import math
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

start = timer()

ROOT.ROOT.EnableImplicitMT(5)

#This is for data that is reduced to a specific particle group and CeBrA detectors Gain matched
#Current Columns in the SPSTree
#  "anodeBack","scintLeft","scintLeftTime","cathode","xavg","x1","x2","theta",
#  "cebraE0","cebraE1","cebraE2","cebraE3","cebraE4","cebraTime0","cebraTime1","cebraTime2","cebraTime3","cebraTime4",
#  "cebraTime0_toScint","cebraTime1_toScint","cebraTime2_toScint","cebraTime3_toScint","cebraTime4_toScint",
#  "cebraE0_GM","cebraE1_GM","cebraE2_GM","cebraE3_GM","cebraE4_GM"

# Ti49_RunList = [184,185,186,187,188,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,\
#         206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232]
                
# ECal_Gamma_Histo = [1024,0,4096]
ECal_Gamma_Histo = [500,0,6500]
ECal_CeBrA_Para = [0.000034916,1.463915338,49.7375612844]

def analysis_basic(Path_to_Tree, TreeName, Particle, NumDetectors, ECal_FP_List, ECal_Focal_Histo, ECal_FP, output):
    
    Reduced_Frame = ROOT.RDataFrame(TreeName,f"{Path_to_Tree}*.root")\
        .Define(f"{ECal_FP}_ECal", f"{ECal_FP_List[0]}*{ECal_FP}*{ECal_FP} + {ECal_FP_List[1]}*{ECal_FP} + {ECal_FP_List[2]}")\
        .Define("cebraE0_ECal", f"{ECal_CeBrA_Para[0]}*cebraE0_GM*cebraE0_GM + {ECal_CeBrA_Para[1]}*cebraE0_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE1_ECal", f"{ECal_CeBrA_Para[0]}*cebraE1_GM*cebraE1_GM + {ECal_CeBrA_Para[1]}*cebraE1_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE2_ECal", f"{ECal_CeBrA_Para[0]}*cebraE2_GM*cebraE2_GM + {ECal_CeBrA_Para[1]}*cebraE2_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE3_ECal", f"{ECal_CeBrA_Para[0]}*cebraE3_GM*cebraE3_GM + {ECal_CeBrA_Para[1]}*cebraE3_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE4_ECal", f"{ECal_CeBrA_Para[0]}*cebraE4_GM*cebraE4_GM + {ECal_CeBrA_Para[1]}*cebraE4_GM + {ECal_CeBrA_Para[2]}")
    
    PID_particle = Reduced_Frame.Histo2D((f"ScintL_AnodeBack_{Particle}",f"{Particle}",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")                                                                 
    
    xavg_particle = Reduced_Frame.Histo1D((f"xavg_{Particle}",f"xavg_{Particle}",600,-300,300),"xavg")
    x1_bothplanes_particle = Reduced_Frame.Filter("x1 != -1e6 && x2 != -1e6").Histo1D((f"x1_bothplanes_{Particle}",f"x1_bothplanes_{Particle}",600,-300,300),"x1")
    x2_bothplanes_particle = Reduced_Frame.Filter("x1 != -1e6 && x2 != -1e6").Histo1D((f"x2_bothplanes_{Particle}",f"x2_bothplanes_{Particle}",600,-300,300),"x2")
    x1_only1plane_particle = Reduced_Frame.Filter("x1 != -1e6 && x2 == -1e6").Histo1D((f"x1_only1plane_{Particle}",f"x1_only1plane_{Particle}",600,-300,300),"x1")
    
    x1_theta_particle = Reduced_Frame.Histo2D((f"x1_theta_{Particle}",f"x1_theta_{Particle}",600,-300,300,600,0,pi/2),"x1","theta")
    xavg_theta_particle = Reduced_Frame.Histo2D((f"xavg_theta_{Particle}",f"xavg_theta_{Particle}",600,-300,300,600,0,pi/2),"xavg","theta")
    
    FP_particle_ECal = Reduced_Frame.Histo1D((f"{ECal_FP}_{Particle}_ECal",f"{ECal_FP}_{Particle}_ECal",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),f"{ECal_FP}_ECal")
    FP_ECal_theta_particle = Reduced_Frame.Histo2D((f"{ECal_FP}_ECal_theta_{Particle}",f"{ECal_FP}_ECal_theta_{Particle}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],600,0,pi/2),f"{ECal_FP}_ECal","theta")

    cebraTime_toScint_particle = []
    for i in range(NumDetectors):
        cebraTime_toScint_particle.append(Reduced_Frame.Histo1D((f"cebraTime{i}_toScint_{Particle}",f"cebraTime{i}_toScint_{Particle}",6400,-3200,3200),f"cebraTime{i}_toScint"))
        
    output.cd()
    output.mkdir(f"49Ti_{Particle}")
    output.mkdir(f"49Ti_{Particle}/Basic_Particle_Plots")
    output.mkdir(f"49Ti_{Particle}/Basic_Particle_Plots/Time_Correlation")
    output.mkdir(f"49Ti_{Particle}/Basic_Particle_Plots/Energy_Calibrated")
    
    output.cd(f"49Ti_{Particle}/Basic_Particle_Plots")
    PID_particle.Write()
    xavg_particle.Write()
    x1_bothplanes_particle.Write()
    x2_bothplanes_particle.Write()
    x1_only1plane_particle.Write()
    x1_theta_particle.Write()
    xavg_theta_particle.Write()
    
    output.cd(f"49Ti_{Particle}/Basic_Particle_Plots/Time_Correlation")
    for i in range(NumDetectors):
        cebraTime_toScint_particle[i].Write()
        
    output.cd(f"49Ti_{Particle}/Basic_Particle_Plots/Energy_Calibrated")
    FP_particle_ECal.Write()
    FP_ECal_theta_particle.Write()
    
    output.Close()
       
def cebraTime_toScint_Gate(detector,lowvalue,highvalue):
    TimeGate = f"cebraTime{detector}_toScint > {lowvalue} && cebraTime{detector}_toScint < {highvalue}"
    return TimeGate

def analysis(Path_to_Tree, TreeName, Particle, NumDetectors, FP_Variable, ECal_FP_List, ECal_Focal_Histo, cebraTimeGate, output):
    
    Reduced_Frame = ROOT.RDataFrame(TreeName,f"{Path_to_Tree}*.root")\
        .Define(f"{FP_Variable}_ECal", f"{ECal_FP_List[0]}*{FP_Variable}*{FP_Variable} + {ECal_FP_List[1]}*{FP_Variable} + {ECal_FP_List[2]}")\
        .Define("cebraE0_ECal", f"{ECal_CeBrA_Para[0]}*cebraE0_GM*cebraE0_GM + {ECal_CeBrA_Para[1]}*cebraE0_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE1_ECal", f"{ECal_CeBrA_Para[0]}*cebraE1_GM*cebraE1_GM + {ECal_CeBrA_Para[1]}*cebraE1_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE2_ECal", f"{ECal_CeBrA_Para[0]}*cebraE2_GM*cebraE2_GM + {ECal_CeBrA_Para[1]}*cebraE2_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE3_ECal", f"{ECal_CeBrA_Para[0]}*cebraE3_GM*cebraE3_GM + {ECal_CeBrA_Para[1]}*cebraE3_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE4_ECal", f"{ECal_CeBrA_Para[0]}*cebraE4_GM*cebraE4_GM + {ECal_CeBrA_Para[1]}*cebraE4_GM + {ECal_CeBrA_Para[2]}")
    
    #FP_variable would either be x1 or xavg

    det_i_Reduced_Frame = []
    cebraTime_toScint_particle_TCut = []
    cebraE_particle_TCut = []
    FP_Variable_particle_TCut = []
    FP_Variable_cebraE_particle_TCut = []

    #Energy Calibrated
    FP_Variable_ECal_particle_TCut = []
    cebraE_ECal_particle_TCut = []
    FP_Variable_cebraE_ECal_particle_TCut = []
    
    for i in range(NumDetectors):
    
        det_i_Reduced_Frame.append(Reduced_Frame.Filter(cebraTimeGate[i]))
        cebraTime_toScint_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"cebraTime{i}_toScint_{Particle}_TCut{i}",f"cebraTime{i}_toScint_{Particle}_TCut{i}",6400,-3200,3200),f"cebraTime{i}_toScint"))
        cebraE_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"cebraE{i}_{Particle}_TCut",f"cebraE{i}_{Particle}_TCut",4096,0,4096),f"cebraE{i}_GM"))
        FP_Variable_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"{FP_Variable}_{Particle}_TCut{i}",f"{FP_Variable}_{Particle}_TCut{i}",600,-300,300),f"{FP_Variable}"))
        FP_Variable_cebraE_particle_TCut.append(det_i_Reduced_Frame[i].Histo2D((f"{FP_Variable}_cebraE{i}_{Particle}_TCut{i}",f"{FP_Variable}_cebraE{i}_{Particle}_TCut{i}",600,-300,300,4096,0,4096), f"{FP_Variable}",f"cebraE{i}_GM"))
        
        #Energy Calibrated
        FP_Variable_ECal_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"{FP_Variable}_ECal_{Particle}_TCut{i}",f"{FP_Variable}_ECal_{Particle}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),f"{FP_Variable}_ECal"))
        cebraE_ECal_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"cebraE{i}_ECal_{Particle}_TCut",f"cebraE{i}_ECal_{Particle}_TCut",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))
        FP_Variable_cebraE_ECal_particle_TCut.append(det_i_Reduced_Frame[i].Histo2D((f"{FP_Variable}_cebraE{i}_ECal_{Particle}_TCut{i}",f"{FP_Variable}_cebraE{i}_ECal_{Particle}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"{FP_Variable}_ECal",f"cebraE{i}_ECal"))
        
            
        if i==0:
            #initilizing CeBrA Plots
            
            #Not Calibrated
            Array_FP_Variable_particle_TCut = det_i_Reduced_Frame[i].Histo1D((f"Array_{FP_Variable}_{Particle}_TCut",f"Array_{FP_Variable}_{Particle}_TCut",600,-300,300),f"{FP_Variable}")
            Array_cebraE_particle_TCut = det_i_Reduced_Frame[i].Histo1D((f"Array_cebraE_{Particle}_TCut",f"Array_cebraE_{Particle}_TCut",4096,0,4096),f"cebraE{i}_GM")
            Array_FP_Variable_cebraE_particle_TCut = det_i_Reduced_Frame[i].Histo2D((f"Array_{FP_Variable}_cebraE_{Particle}_TCut",f"Array_{FP_Variable}_cebraE_{Particle}_TCut",600,-300,300,4096,0,4096), f"{FP_Variable}",f"cebraE{i}_GM")
            
            #Energy Calibrated
            Array_FP_Variable_ECal_particle = det_i_Reduced_Frame[i].Histo1D((f"Array_{FP_Variable}_ECal_{Particle}",f"Array_{FP_Variable}_ECal_{Particle}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),f"{FP_Variable}_ECal")
            Array_cebraE_ECal_particle = det_i_Reduced_Frame[i].Histo1D((f"Array_cebraE_ECal_{Particle}",f"Array_cebraE_ECal_{Particle}",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal")
            Array_FP_Variable_cebraE_ECal_particle = det_i_Reduced_Frame[i].Histo2D((f"{FP_Variable}_cebraE_ECal_Array",f"{FP_Variable}_cebraE_ECal_Array",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), f"{FP_Variable}_ECal",f"cebraE{i}_ECal")

    
        if i!=0:
            #Not Calibrated
            Array_FP_Variable_particle_TCut.Add(FP_Variable_particle_TCut[i].GetPtr())
            Array_cebraE_particle_TCut.Add(cebraE_particle_TCut[i].GetPtr())
            Array_FP_Variable_cebraE_particle_TCut.Add(FP_Variable_cebraE_particle_TCut[i].GetPtr())
            
            #Energy Calibrated
            Array_FP_Variable_ECal_particle.Add(FP_Variable_ECal_particle_TCut[i].GetPtr())
            Array_cebraE_ECal_particle.Add(cebraE_ECal_particle_TCut[i].GetPtr())
            Array_FP_Variable_cebraE_ECal_particle.Add(FP_Variable_cebraE_ECal_particle_TCut[i].GetPtr())
            
            
    output.cd()
    output.mkdir(f"49Ti_{Particle}")
    output.mkdir(f"49Ti_{Particle}/Particle_Time_Cut_Reduced/")
    output.mkdir(f"49Ti_{Particle}/Particle_Time_Cut_Reduced/Time_Correlation")
    output.mkdir(f"49Ti_{Particle}/Particle_Time_Cut_Reduced/Individual_Spectrums")
    
    for j in range(NumDetectors):
        
        output.cd(f"49Ti_{Particle}/Particle_Time_Cut_Reduced/Time_Correlation")
        cebraTime_toScint_particle_TCut[j].Write()
        
        output.cd(f"49Ti_{Particle}/Particle_Time_Cut_Reduced/Individual_Spectrums")
        cebraE_particle_TCut[j].Write()
        FP_Variable_particle_TCut[j].Write()
        FP_Variable_cebraE_particle_TCut[j].Write()
    
    output.mkdir(f"49Ti_{Particle}/Particle_Time_Cut_Reduced/CeBrA")
    output.cd(f"49Ti_{Particle}/Particle_Time_Cut_Reduced/CeBrA")
    
    Array_FP_Variable_particle_TCut.Write()
    Array_cebraE_particle_TCut.Write()
    Array_FP_Variable_cebraE_particle_TCut.Write()
    
    output.mkdir(f"49Ti_{Particle}/CeBrA_Energy_Calibrated")
    output.cd(f"49Ti_{Particle}/CeBrA_Energy_Calibrated")
    Array_FP_Variable_ECal_particle.Write()
    Array_cebraE_ECal_particle.Write()
    Array_FP_Variable_cebraE_ECal_particle.Write()
    
    output.Close()
                         
def bands(Path_to_Tree,TreeName,Particle,NumDetectors,FP_Variable,ECal_FP_List, ECal_Focal_Histo, cebraTimeGate, Band_Name, x1, y1, x2, y2, Width, output):

    Reduced_Frame = ROOT.RDataFrame(TreeName,f"{Path_to_Tree}*.root")\
        .Define(f"{FP_Variable}_ECal", f"{ECal_FP_List[0]}*{FP_Variable}*{FP_Variable} + {ECal_FP_List[1]}*{FP_Variable} + {ECal_FP_List[2]}")\
        .Define("cebraE0_ECal", f"{ECal_CeBrA_Para[0]}*cebraE0_GM*cebraE0_GM + {ECal_CeBrA_Para[1]}*cebraE0_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE1_ECal", f"{ECal_CeBrA_Para[0]}*cebraE1_GM*cebraE1_GM + {ECal_CeBrA_Para[1]}*cebraE1_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE2_ECal", f"{ECal_CeBrA_Para[0]}*cebraE2_GM*cebraE2_GM + {ECal_CeBrA_Para[1]}*cebraE2_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE3_ECal", f"{ECal_CeBrA_Para[0]}*cebraE3_GM*cebraE3_GM + {ECal_CeBrA_Para[1]}*cebraE3_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE4_ECal", f"{ECal_CeBrA_Para[0]}*cebraE4_GM*cebraE4_GM + {ECal_CeBrA_Para[1]}*cebraE4_GM + {ECal_CeBrA_Para[2]}")
    
    #Band_Name, FP_Variable, and Y_Variable have to be strings
        
    #Basically give the function two points and it will create a filter condition
    #allowed the x variable to change incase we want to use xavg
    #x = "x1_ECal"
    
    #slope
    m = (y2-y1)/(x2-x1)
    
    Band_Frame = []
    Band_X = []
    Band_Y = []
    Band_YX = []
    
    Band_Filter_Condition = []
    
    for i in range(NumDetectors):
    
        Band_Filter_Condition.append(f"( (cebraE{i}_ECal - {y1}) <= ({m}*({FP_Variable}_ECal-{x1}) + {Width} )) && ( (cebraE{i}_ECal - {y1}) >= ({m}*({FP_Variable}_ECal-{x1}) - {Width}) )")
        
        Band_Frame.append(Reduced_Frame.Filter(cebraTimeGate[i]).Filter(Band_Filter_Condition[i]))

        Band_X.append(Band_Frame[i].Histo1D((f"{FP_Variable}_{Band_Name}_TCut{i}",f"{FP_Variable}_{Band_Name}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]), f"{FP_Variable}_ECal"))
        Band_Y.append(Band_Frame[i].Histo1D((f"cebraE{i}_ECal_{Band_Name}_TCut",f"cebraE{i}_ECal_{Band_Name}_TCut",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))
        Band_YX.append(Band_Frame[i].Histo2D((f"{FP_Variable}_cebraE{i}_ECal_{Band_Name}_TCut{i}",f"{FP_Variable}_cebraE{i}_ECal_{Band_Name}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), f"{FP_Variable}_ECal",f"cebraE{i}_ECal"))
        
        if i==0:
            X_Sum_Plot = Band_Frame[i].Histo1D((f"Array_{FP_Variable}_{Band_Name}",f"Array_{FP_Variable}_{Band_Name}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),f"{FP_Variable}_ECal")
            Y_Sum_Plot = Band_Frame[i].Histo1D((f"Array_cebraE_{Band_Name}",f"Array_cebraE_{Band_Name}",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal")
            YX_Sum_Plot = Band_Frame[i].Histo2D((f"Array_cebraE_{FP_Variable}_{Band_Name}",f"Array_cebraE_{FP_Variable}_{Band_Name}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), f"{FP_Variable}_ECal",f"cebraE{i}_ECal")

        if i!=0:
            X_Sum_Plot.Add(Band_X[i].GetPtr())
            Y_Sum_Plot.Add(Band_Y[i].GetPtr())
            YX_Sum_Plot.Add(Band_YX[i].GetPtr())
    
    output.cd()
    output.mkdir(f"49Ti_{Particle}/Bands/")
    output.mkdir(f"49Ti_{Particle}/Bands/{Band_Name}/")
    output.cd(f"49Ti_{Particle}/Bands/{Band_Name}/")
    
    X_Sum_Plot.Write()
    Y_Sum_Plot.Write()
    YX_Sum_Plot.Write()
        
    return Band_Filter_Condition
        
def Excited_State(Path_to_Tree, TreeName, Particle, NumDetectors, FP_Variable, ECal_FP_List, ECal_Focal_Histo, cebraTimeGate, Excited_State_Name, x1, x2, output):
    
    Reduced_Frame = ROOT.RDataFrame(TreeName,f"{Path_to_Tree}*.root")\
        .Define(f"{FP_Variable}_ECal", f"{ECal_FP_List[0]}*{FP_Variable}*{FP_Variable} + {ECal_FP_List[1]}*{FP_Variable} + {ECal_FP_List[2]}")\
        .Define("cebraE0_ECal", f"{ECal_CeBrA_Para[0]}*cebraE0_GM*cebraE0_GM + {ECal_CeBrA_Para[1]}*cebraE0_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE1_ECal", f"{ECal_CeBrA_Para[0]}*cebraE1_GM*cebraE1_GM + {ECal_CeBrA_Para[1]}*cebraE1_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE2_ECal", f"{ECal_CeBrA_Para[0]}*cebraE2_GM*cebraE2_GM + {ECal_CeBrA_Para[1]}*cebraE2_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE3_ECal", f"{ECal_CeBrA_Para[0]}*cebraE3_GM*cebraE3_GM + {ECal_CeBrA_Para[1]}*cebraE3_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE4_ECal", f"{ECal_CeBrA_Para[0]}*cebraE4_GM*cebraE4_GM + {ECal_CeBrA_Para[1]}*cebraE4_GM + {ECal_CeBrA_Para[2]}")
    
    Excited_State_Filter = f"{FP_Variable}_ECal > {x1} && {FP_Variable}_ECal < {x2}"

    Excited_State_Frame = []
    Excited_State_cebra = []
    Excited_State_fp_variable = []
    
    for i in range(NumDetectors):
        Excited_State_Frame.append(Reduced_Frame.Filter(Excited_State_Filter).Filter(cebraTimeGate[i]))
        Excited_State_fp_variable.append(Excited_State_Frame[i].Histo1D((f"{Excited_State_Name}_keV_{FP_Variable}_ECal",f"{Excited_State_Name}_keV_{FP_Variable}_ECal",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]), f"{FP_Variable}_ECal"))
        Excited_State_cebra.append(Excited_State_Frame[i].Histo1D((f"{Excited_State_Name}_keV_cebra",f"{Excited_State_Name}_keV_cebra",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))
        
        if i==0:
            X_Sum_Plot =  Excited_State_Frame[i].Histo1D((f"{Excited_State_Name}_keV_{FP_Variable}_ECal",f"{Excited_State_Name}_keV_{FP_Variable}_ECal",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),f"{FP_Variable}_ECal")
            Y_Sum_Plot =  Excited_State_Frame[i].Histo1D((f"{Excited_State_Name}_keV_CeBrA_ECal",f"{Excited_State_Name}_keV_CeBrA_ECal", ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), f"cebraE{i}_ECal")

        if i!=0:
            X_Sum_Plot.Add(Excited_State_fp_variable[i].GetPtr())
            Y_Sum_Plot.Add(Excited_State_cebra[i].GetPtr())
     
    output.cd()
    output.mkdir(f"49Ti_{Particle}/50Ti_Excited_States/{Excited_State_Name}_keV/")
    output.cd(f"49Ti_{Particle}/50Ti_Excited_States/{Excited_State_Name}_keV/")
    
    X_Sum_Plot.Write()
    Y_Sum_Plot.Write()
    
    return Excited_State_Filter
  
def Excited_States_with_Bands(Path_to_Tree, TreeName, Particle, NumDetectors,FP_Variable, ECal_FP_List, ECal_Focal_Histo, cebraTimeGate,Excited_State_Name, Band_Filter, Excited_State_Filter, Name, output):
    
    Reduced_Frame = ROOT.RDataFrame(TreeName,f"{Path_to_Tree}*.root")\
        .Define(f"{FP_Variable}_ECal", f"{ECal_FP_List[0]}*{FP_Variable}*{FP_Variable} + {ECal_FP_List[1]}*{FP_Variable} + {ECal_FP_List[2]}")\
        .Define("cebraE0_ECal", f"{ECal_CeBrA_Para[0]}*cebraE0_GM*cebraE0_GM + {ECal_CeBrA_Para[1]}*cebraE0_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE1_ECal", f"{ECal_CeBrA_Para[0]}*cebraE1_GM*cebraE1_GM + {ECal_CeBrA_Para[1]}*cebraE1_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE2_ECal", f"{ECal_CeBrA_Para[0]}*cebraE2_GM*cebraE2_GM + {ECal_CeBrA_Para[1]}*cebraE2_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE3_ECal", f"{ECal_CeBrA_Para[0]}*cebraE3_GM*cebraE3_GM + {ECal_CeBrA_Para[1]}*cebraE3_GM + {ECal_CeBrA_Para[2]}")\
        .Define("cebraE4_ECal", f"{ECal_CeBrA_Para[0]}*cebraE4_GM*cebraE4_GM + {ECal_CeBrA_Para[1]}*cebraE4_GM + {ECal_CeBrA_Para[2]}")
    
    Excited_State_Band_Frame = []
    Band_X = []
    Band_Y = []
    Band_YX = []
    
    for i in range(NumDetectors):
            
        Excited_State_Band_Frame.append(Reduced_Frame.Filter(cebraTimeGate[i]).Filter(Band_Filter[i]).Filter(Excited_State_Filter))

        Band_X.append(Excited_State_Band_Frame[i].Histo1D((f"{FP_Variable}_{Name}_TCut{i}",f"{FP_Variable}_{Name}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]), f"{FP_Variable}_ECal"))
        Band_Y.append(Excited_State_Band_Frame[i].Histo1D((f"cebraE{i}_ECal_{Name}_TCut",f"cebraE{i}_ECal_{Name}_TCut",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))
        Band_YX.append(Excited_State_Band_Frame[i].Histo2D((f"{FP_Variable}_cebraE{i}_ECal_{Name}_TCut{i}",f"{FP_Variable}_cebraE{i}_ECal_{Name}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), f"{FP_Variable}_ECal",f"cebraE{i}_ECal"))
        
        if i==0:
            X_Sum_Plot = Excited_State_Band_Frame[i].Histo1D((f"{Name}_{FP_Variable}",f"{Name}_{FP_Variable}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),f"{FP_Variable}_ECal")
            Y_Sum_Plot = Excited_State_Band_Frame[i].Histo1D((f"{Name}_cebraE",f"{Name}_cebraE",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal")
            YX_Sum_Plot = Excited_State_Band_Frame[i].Histo2D((f"{Name}_cebraE_{FP_Variable}",f"{Name}_cebraE_{FP_Variable}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), f"{FP_Variable}_ECal",f"cebraE{i}_ECal")

        if i!=0:
            X_Sum_Plot.Add(Band_X[i].GetPtr())
            Y_Sum_Plot.Add(Band_Y[i].GetPtr())
            YX_Sum_Plot.Add(Band_YX[i].GetPtr())
    
    output.cd(f"49Ti_{Particle}/50Ti_Excited_States/{Excited_State_Name}_keV/")
    
    X_Sum_Plot.Write()
    Y_Sum_Plot.Write()
    YX_Sum_Plot.Write()
        
# Protons
################################################################################################################################################
print("Protons\n")

Proton_Path = f"/home/alconley/Projects/CeBrA_Analysis/Sept_2022_49Ti/49Ti_Analyzed_Protons_Reduced/"
#ax^2+bx+c, [a,b,c]
ECal_Focal_Para_Protons = [-0.003824652423017848,-21.158024905122485,2157.116872413482]
ECal_Focal_Histo_Protons = [406,-100,6396]

cebraTimeGate_Protons = [cebraTime_toScint_Gate(0,-678,-670),\
    cebraTime_toScint_Gate(1,-676,-669),\
    cebraTime_toScint_Gate(2,-676,-669),\
    cebraTime_toScint_Gate(3,-676,-668),\
    cebraTime_toScint_Gate(4,-645,-638)]

print("Basic Plots")
output_basic_path_protons = f"49Ti_Protons_basic_histo.root"
output_basic_protons = ROOT.TFile.Open(output_basic_path_protons, "RECREATE")
analysis_basic(Proton_Path, "SPSTree", "Protons", 5, ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, "x1", output_basic_protons)

print("Plots with the time cuts")
output_timecut_path_protons = f"49Ti_Protons_timecut_histo.root"
output_timecut_protons = ROOT.TFile.Open(output_timecut_path_protons, "RECREATE")
analysis(Proton_Path, "SPSTree", "Protons", 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons, output_timecut_protons)

print("Band Plots")
output_bands_path_protons = f"49Ti_Protons_bands_histo.root"
output_bands_protons = ROOT.TFile.Open(output_bands_path_protons, "RECREATE")
output_bands_protons.mkdir(f"49Ti_Protons")

Ti49_GS_Band_Filter = bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons,ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"GS_Band",0,0,1553,1553,50, output_bands_protons)
Ti49_2_Plus_Band_Filter = bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons,ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"2_Plus_Band",1553,0,4172,2618,50, output_bands_protons)
Ti49_4_Plus_Band_Filter = bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons,ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4_Plus_Band",2674,0,4147,1472,50, output_bands_protons)
Ti49_6_Plus_Band_Filter = bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons,ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"6_Plus_Band",3198,0,6123,2924.9,50, output_bands_protons)
output_bands_protons.Close()

print('Exctied State plots')
output_excited_states_path_protons = f"49Ti_Protons_Excited_States_histo.root"
output_excited_states_protons = ROOT.TFile.Open(output_excited_states_path_protons, "RECREATE")
output_excited_states_protons.mkdir(f"49Ti_Protons")
output_excited_states_protons.mkdir(f"49Ti_Protons/50Ti_Excited_States/")

fp_filter_0_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"0",-100,100, output_excited_states_protons)
fp_filter_1553_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"1553",1400,1700, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"1553",Ti49_GS_Band_Filter,fp_filter_1553_keV,"1553_kev_GS_Band",output_excited_states_protons)

fp_filter_2674_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"2674",2520,2760, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"2674",Ti49_2_Plus_Band_Filter,fp_filter_2674_keV,"2674_kev_2_Plus_Band",output_excited_states_protons)

fp_filter_4100_doublet_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4100_doublet",4050,4275, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4100_doublet",Ti49_GS_Band_Filter,fp_filter_4100_doublet_keV,"4100_doublet_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4100_doublet",Ti49_2_Plus_Band_Filter,fp_filter_4100_doublet_keV,"4100_doublet_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4100_doublet",Ti49_4_Plus_Band_Filter,fp_filter_4100_doublet_keV,"4100_doublet_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_4500_doublet_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4500_doublet",4450,4650, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4500_doublet",Ti49_GS_Band_Filter,fp_filter_4500_doublet_keV,"4500_doublet_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4500_doublet",Ti49_2_Plus_Band_Filter,fp_filter_4500_doublet_keV,"4500_doublet_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4500_doublet",Ti49_4_Plus_Band_Filter,fp_filter_4500_doublet_keV,"4500_doublet_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_4880_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4880",4720,5000, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4880",Ti49_GS_Band_Filter,fp_filter_4880_keV,"4880_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4880",Ti49_2_Plus_Band_Filter,fp_filter_4880_keV,"4880_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"4880",Ti49_4_Plus_Band_Filter,fp_filter_4880_keV,"4880_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_5186_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5186",5120,5240, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5186",Ti49_GS_Band_Filter,fp_filter_5186_keV,"5186_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5186",Ti49_2_Plus_Band_Filter,fp_filter_5186_keV,"5186_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5186",Ti49_4_Plus_Band_Filter,fp_filter_5186_keV,"5186_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_5379_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5379",5322,5460, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5379",Ti49_GS_Band_Filter,fp_filter_5379_keV,"5379_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5379",Ti49_2_Plus_Band_Filter,fp_filter_5379_keV,"5379_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5379",Ti49_4_Plus_Band_Filter,fp_filter_5379_keV,"5379_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_5500_5630_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5500_5630",5500,5630, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5500_5630",Ti49_GS_Band_Filter,fp_filter_5500_5630_keV,"5500_5630_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5500_5630",Ti49_2_Plus_Band_Filter,fp_filter_5500_5630_keV,"5500_5630_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5500_5630",Ti49_4_Plus_Band_Filter,fp_filter_5500_5630_keV,"5500_5630_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_5660_5780_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5660_5780",5660,5780, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5660_5780",Ti49_GS_Band_Filter,fp_filter_5660_5780_keV,"5660_5780_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5660_5780",Ti49_2_Plus_Band_Filter,fp_filter_5660_5780_keV,"5660_5780_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5660_5780",Ti49_4_Plus_Band_Filter,fp_filter_5660_5780_keV,"5660_5780_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_5837_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5837",5780,5880, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5837",Ti49_GS_Band_Filter,fp_filter_5837_keV,"5837_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5837",Ti49_2_Plus_Band_Filter,fp_filter_5837_keV,"5837_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5837",Ti49_4_Plus_Band_Filter,fp_filter_5837_keV,"5837_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_5946_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5946",5890,6020, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5946",Ti49_GS_Band_Filter,fp_filter_5946_keV,"5946_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5946",Ti49_2_Plus_Band_Filter,fp_filter_5946_keV,"5946_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"5946",Ti49_4_Plus_Band_Filter,fp_filter_5946_keV,"5946_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_6020_6100_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"6020_6100",6020,6100, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"6020_6100",Ti49_GS_Band_Filter,fp_filter_6020_6100_keV,"6020_6100_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"6020_6100",Ti49_2_Plus_Band_Filter,fp_filter_6020_6100_keV,"6020_6100_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"6020_6100",Ti49_4_Plus_Band_Filter,fp_filter_6020_6100_keV,"6020_6100_kev_4_Plus_Band",output_excited_states_protons)

fp_filter_6123_keV = Excited_State(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"6123",6100,6190, output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"6123",Ti49_GS_Band_Filter,fp_filter_6123_keV,"6123_kev_GS_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"6123",Ti49_2_Plus_Band_Filter,fp_filter_6123_keV,"6123_kev_2_Plus_Band",output_excited_states_protons)
Excited_States_with_Bands(Proton_Path, "SPSTree", "Protons" , 5, "x1", ECal_Focal_Para_Protons, ECal_Focal_Histo_Protons, cebraTimeGate_Protons,"6123",Ti49_4_Plus_Band_Filter,fp_filter_6123_keV,"6123_kev_4_Plus_Band",output_excited_states_protons)

 
output_excited_states_protons.Close()



# Deuterons
################################################################################################################################################

Deuteron_Path = f"/home/alconley/Projects/CeBrA_Analysis/Sept_2022_49Ti/49Ti_Analyzed_Deuterons_Reduced/"
#ax^2+bx+c, [a,b,c]
ECal_Focal_Para_Deuterons = [0,1,0]
ECal_Focal_Histo_Deuterons = [600,-300,300]

cebraTimeGate_Deuterons = [cebraTime_toScint_Gate(0,-733,-720),\
    cebraTime_toScint_Gate(1,-733,-717),\
    cebraTime_toScint_Gate(2,-731,-717),\
    cebraTime_toScint_Gate(3,-730,-718),\
    cebraTime_toScint_Gate(4,-700,-688)]

# output_basic_path_deuterons = f"49Ti_Deuterons_basic_histo.root"
# output_basic_deuterons = ROOT.TFile.Open(output_basic_path_deuterons, "RECREATE")
# analysis_basic(Deuteron_Path, "SPSTree", "Deuterons", 5, ECal_Focal_Para_Deuterons, ECal_Focal_Histo_Deuterons, "xavg", output_basic_deuterons)

# output_timecut_path_deuterons = f"49Ti_Deuterons_timecut_histo.root"
# output_timecut_path_deuterons = ROOT.TFile.Open(output_timecut_path_deuterons, "RECREATE")
# analysis(Deuteron_Path, "SPSTree", "Deuterons", 5, "xavg", ECal_Focal_Para_Deuterons, ECal_Focal_Histo_Deuterons, cebraTimeGate_Deuterons, output_timecut_path_deuterons)

end = timer()
print("Total Elapsed time: ",end - start)