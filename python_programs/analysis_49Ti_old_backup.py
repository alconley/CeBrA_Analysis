from cmath import pi
import ROOT
import math
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

#This is for data that is reduced to a specific particle group and CeBrA detectors Gain matched
#Current Columns in the SPSTree
#  "anodeBack","scintLeft","scintLeftTime","cathode","xavg","x1","x2","theta",
#  "cebraE0","cebraE1","cebraE2","cebraE3","cebraE4","cebraTime0","cebraTime1","cebraTime2","cebraTime3","cebraTime4",
#  "cebraTime0_toScint","cebraTime1_toScint","cebraTime2_toScint","cebraTime3_toScint","cebraTime4_toScint",
#  "cebraE0_GM","cebraE1_GM","cebraE2_GM","cebraE3_GM","cebraE4_GM"

Ti49_RunList = [184,185,186,187,188,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,\
        206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232]
                
# Ti49_RunList=[184]

#number of detectors
n_det = 5
treeName = "SPSTree"

#[bins,initial,final]
# ECal_Gamma_Histo = [1024,0,4096]
ECal_Gamma_Histo = [500,0,6500]
ECal_Focal_Histo = [406,-100,6396]

#ax^2+bx+c, [a,b,c]
ECal_Focal_Para = [-0.003824652423017848,-21.158024905122485,2157.116872413482]
ECal_CeBrA_Para = [0.000034916,1.463915338,49.7375612844]
# ECal_CeBrA_Para = [0, 1, 0]



def analysis_basic(RunNumber,Path_to_Tree,Particle,ECal_Focalplane):
    
    
    
    PID_particle = Reduced_Frame.Histo2D((f"ScintL_AnodeBack_{Particle}",f"{Particle}",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")                                                                 
    
    xavg_particle = Reduced_Frame.Histo1D((f"xavg_{Particle}",f"xavg_{Particle}",600,-300,300),"xavg")
    x1_bothplanes_particle = Reduced_Frame.Filter("x1 != -1e6 && x2 != -1e6").Histo1D((f"x1_bothplanes_{Particle}",f"x1_bothplanes_{Particle}",600,-300,300),"x1")
    x2_bothplanes_particle = Reduced_Frame.Filter("x1 != -1e6 && x2 != -1e6").Histo1D((f"x2_bothplanes_{Particle}",f"x2_bothplanes_{Particle}",600,-300,300),"x2")
    x1_only1plane_particle = Reduced_Frame.Filter("x1 != -1e6 && x2 == -1e6").Histo1D((f"x1_only1plane_{Particle}",f"x1_only1plane_{Particle}",600,-300,300),"x1")
    
    x1_theta_particle = Reduced_Frame.Histo2D((f"x1_theta_{Particle}",f"x1_theta_{Particle}",600,-300,300,600,0,pi/2),"x1","theta")
    xavg_theta_particle = Reduced_Frame.Histo2D((f"xavg_theta_{Particle}",f"xavg_theta_{Particle}",600,-300,300,600,0,pi/2),"xavg","theta")
    
    focalplane_particle_ECal = Reduced_Frame.Histo1D((f"{ECal_Focalplane}_{Particle}_ECal",f"{ECal_Focalplane}_{Particle}_ECal",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),f"{ECal_Focalplane}_ECal")
    focalplane_ECal_theta_particle = Reduced_Frame.Histo2D((f"{ECal_Focalplane}_ECal_theta_{Particle}",f"{ECal_Focalplane}_ECal_theta_{Particle}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],600,0,pi/2),f"{ECal_Focalplane}_ECal","theta")

    cebraTime_toScint_particle = []
    for i in range(n_det):
        cebraTime_toScint_particle.append(Reduced_Frame.Histo1D((f"cebraTime{i}_toScint_{Particle}",f"cebraTime{i}_toScint_{Particle}",6400,-3200,3200),f"cebraTime{i}_toScint"))
        
    return PID_particle, xavg_particle, x1_bothplanes_particle, x2_bothplanes_particle, x1_only1plane_particle, x1_theta_particle, xavg_theta_particle,\
        focalplane_particle_ECal, focalplane_ECal_theta_particle,cebraTime_toScint_particle
        
def cebraTime_toScint_Gate(detector,lowvalue,highvalue):
    TimeGate = f"cebraTime{detector}_toScint > {lowvalue} && cebraTime{detector}_toScint < {highvalue}"
    return TimeGate

def analysis(NumDetectors,Particle,Focalplane):
    #Focalplane would either be x1 or xavg
    #Particle would either be Protons or Deuterons
    
    det_i_Reduced_Frame = []
    cebraTime_toScint_particle_TCut = []
    cebraE_particle_TCut = []
    focalplane_particle_TCut = []
    focalplane_cebraE_particle_TCut = []

    #Energy Calibrated
    focalplane_ECal_particle_TCut = []
    cebraE_ECal_particle_TCut = []
    focalplane_cebraE_ECal_particle_TCut = []
    
    #############################################
    #Stuff that is protons 49Ti(d,pg),x1_ECal Specific
    def fp_Gate(Energy,FocalPlaneVariable,lowvalue,highvalue):
            PositionGate = [Energy,f"{lowvalue} < {FocalPlaneVariable} && {FocalPlaneVariable} < {highvalue}"]
            
            return PositionGate
        
    x1_E_Gates = [fp_Gate("0","x1_ECal",-100,100),\
            fp_Gate("1553","x1_ECal",1400,1700),\
            fp_Gate("2674","x1_ECal",2520,2760),\
            fp_Gate("4100_doublet","x1_ECal",4050,4275),\
            fp_Gate("4500_doublet","x1_ECal",4450,4650),\
            fp_Gate("4880","x1_ECal",4720,5000),\
            fp_Gate("5186","x1_ECal",5120,5240),\
            fp_Gate("5379","x1_ECal",5322,5460),\
            fp_Gate("5500_5630","x1_ECal",5500,5630),\
            fp_Gate("5660_5780","x1_ECal",5660,5780),\
            fp_Gate("5837","x1_ECal",5780,5880),\
            fp_Gate("5946","x1_ECal",5890,6020),\
            fp_Gate("6020_6100","x1_ECal",6020,6100),\
            fp_Gate("6123","x1_ECal",6100,6190)]
    NumExcited = len(x1_E_Gates)
    
    det_i_Excited_State_Frame = []
    x1_E_List_Protons_TCut = []
    cebraE_E_List_Protons_TCut = []
    
    for i in range(NumDetectors):
                        
        det_i_Reduced_Frame.append(Reduced_Frame.Filter(cebraTimeGate[i]))
        cebraTime_toScint_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"cebraTime{i}_toScint_{Particle}_TCut{i}",f"cebraTime{i}_toScint_{Particle}_TCut{i}",6400,-3200,3200),f"cebraTime{i}_toScint"))
        cebraE_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"cebraE{i}_{Particle}_TCut",f"cebraE{i}_{Particle}_TCut",4096,0,4096),f"cebraE{i}_GM"))
        focalplane_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"{Focalplane}_{Particle}_TCut{i}",f"{Focalplane}_{Particle}_TCut{i}",600,-300,300),f"{Focalplane}"))
        focalplane_cebraE_particle_TCut.append(det_i_Reduced_Frame[i].Histo2D((f"{Focalplane}_cebraE{i}_{Particle}_TCut{i}",f"{Focalplane}_cebraE{i}_{Particle}_TCut{i}",600,-300,300,4096,0,4096), f"{Focalplane}",f"cebraE{i}_GM"))
        
        #Energy Calibrated
        focalplane_ECal_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"{Focalplane}_ECal_{Particle}_TCut{i}",f"{Focalplane}_ECal_{Particle}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),f"{Focalplane}_ECal"))
        cebraE_ECal_particle_TCut.append(det_i_Reduced_Frame[i].Histo1D((f"cebraE{i}_ECal_{Particle}_TCut",f"cebraE{i}_ECal_{Particle}_TCut",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))
        focalplane_cebraE_ECal_particle_TCut.append(det_i_Reduced_Frame[i].Histo2D((f"{Focalplane}_cebraE{i}_ECal_{Particle}_TCut{i}",f"{Focalplane}_cebraE{i}_ECal_{Particle}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"{Focalplane}_ECal",f"cebraE{i}_ECal"))
        # focalplane_cebraE_ECal_particle_TCut.append(det_i_Reduced_Frame[i].Histo2D((f"{Focalplane}_cebraE{i}_ECal_{Particle}_TCut{i}",f"{Focalplane}_cebraE{i}_ECal_{Particle}_TCut{i}",600,-300,300,ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"{Focalplane}",f"cebraE{i}_ECal"))
        
        if Particle=="Protons":
            Excited_State_Frame = []
            det_i_x1 = []
            det_i_cebraE = []
            for j in range(0,len(x1_E_Gates)):
                
                Excited_State_Frame.append(det_i_Reduced_Frame[i].Filter(x1_E_Gates[j][1]))
                det_i_x1.append(Excited_State_Frame[j].Histo1D((f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),"x1_ECal"))
                det_i_cebraE.append(Excited_State_Frame[j].Histo1D((f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))
                
            x1_E_List_Protons_TCut.append(det_i_x1)
            cebraE_E_List_Protons_TCut.append(det_i_cebraE)
            det_i_Excited_State_Frame.append(Excited_State_Frame)
            
        if i==0:
            #initilizing CeBrA Plots
            
            #Not Calibrated
            Array_focalplane_particle_TCut = det_i_Reduced_Frame[i].Histo1D((f"Array_{Focalplane}_{Particle}_TCut",f"Array_{Focalplane}_{Particle}_TCut",600,-300,300),f"{Focalplane}")
            Array_cebraE_particle_TCut = det_i_Reduced_Frame[i].Histo1D((f"Array_cebraE_{Particle}_TCut",f"Array_cebraE_{Particle}_TCut",4096,0,4096),f"cebraE{i}_GM")
            Array_focalplane_cebraE_particle_TCut = det_i_Reduced_Frame[i].Histo2D((f"Array_{Focalplane}_cebraE_{Particle}_TCut",f"Array_{Focalplane}_cebraE_{Particle}_TCut",600,-300,300,4096,0,4096), f"{Focalplane}",f"cebraE{i}_GM")
            
            #Energy Calibrated
            Array_focalplane_ECal_particle = det_i_Reduced_Frame[i].Histo1D((f"Array_{Focalplane}_ECal_{Particle}",f"Array_{Focalplane}_ECal_{Particle}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),f"{Focalplane}_ECal")
            Array_cebraE_ECal_particle = det_i_Reduced_Frame[i].Histo1D((f"Array_cebraE_ECal_{Particle}",f"Array_cebraE_ECal_{Particle}",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal")
            Array_focalplane_cebraE_ECal_particle = det_i_Reduced_Frame[i].Histo2D((f"{Focalplane}_cebraE_ECal_Array",f"{Focalplane}_cebraE_ECal_Array",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), f"{Focalplane}_ECal",f"cebraE{i}_ECal")
            # Array_focalplane_cebraE_ECal_particle = det_i_Reduced_Frame[i].Histo2D((f"{Focalplane}_cebraE_ECal_Array",f"{Focalplane}_cebraE_ECal_Array",600,-300,300,ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), f"{Focalplane}",f"cebraE{i}_ECal")

            if Particle=="Protons":
                #Excited States
                Array_x1_E_List = []
                Array_cebraE_E_List = []
                for j in range(0,NumExcited):
                    Array_x1_E_List.append(det_i_Excited_State_Frame[i][j].Histo1D((f"Array_{x1_E_Gates[j][0]}_keV_x1",f"Array_{x1_E_Gates[j][0]}_keV_x1",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),"x1_ECal"))
                    Array_cebraE_E_List.append(det_i_Excited_State_Frame[i][j].Histo1D((f"Array_{x1_E_Gates[j][0]}_keV_cebraE",f"Array_{x1_E_Gates[j][0]}_keV_cebraE",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))

    
        if i!=0:
            #Not Calibrated
            Array_focalplane_particle_TCut.Add(focalplane_particle_TCut[i].GetPtr())
            Array_cebraE_particle_TCut.Add(cebraE_particle_TCut[i].GetPtr())
            Array_focalplane_cebraE_particle_TCut.Add(focalplane_cebraE_particle_TCut[i].GetPtr())
            
            #Energy Calibrated
            Array_focalplane_ECal_particle.Add(focalplane_ECal_particle_TCut[i].GetPtr())
            Array_cebraE_ECal_particle.Add(cebraE_ECal_particle_TCut[i].GetPtr())
            Array_focalplane_cebraE_ECal_particle.Add(focalplane_cebraE_ECal_particle_TCut[i].GetPtr())
            
            if Particle=="Protons":
                for j in range(0,NumExcited):
                    Array_x1_E_List[j].Add(x1_E_List_Protons_TCut[i][j].GetPtr())
                    Array_cebraE_E_List[j].Add(cebraE_E_List_Protons_TCut[i][j].GetPtr())
    
    if Particle=="Deuterons":
        return cebraTime_toScint_particle_TCut, cebraE_particle_TCut, focalplane_particle_TCut, focalplane_cebraE_particle_TCut,\
            focalplane_ECal_particle_TCut, cebraE_ECal_particle_TCut,focalplane_cebraE_ECal_particle_TCut,\
            Array_focalplane_particle_TCut, Array_cebraE_particle_TCut, Array_focalplane_cebraE_particle_TCut,\
            Array_focalplane_ECal_particle,Array_cebraE_ECal_particle, Array_focalplane_cebraE_ECal_particle
            
    if Particle=="Protons":
        return cebraTime_toScint_particle_TCut, cebraE_particle_TCut, focalplane_particle_TCut, focalplane_cebraE_particle_TCut,\
            focalplane_ECal_particle_TCut, cebraE_ECal_particle_TCut,focalplane_cebraE_ECal_particle_TCut,\
            Array_focalplane_particle_TCut, Array_cebraE_particle_TCut, Array_focalplane_cebraE_particle_TCut,\
            Array_focalplane_ECal_particle,Array_cebraE_ECal_particle, Array_focalplane_cebraE_ECal_particle,\
                Array_x1_E_List,Array_cebraE_E_List,NumExcited
            
def bands(NumDetectors,Band_Name,X_Variable,x1,y1,x2,y2,Width):
    #Band_Name, X_Variable, and Y_Variable have to be strings
        
    #Basically give the function two points and it will create a filter condition
    #allowed the x variable to change incase we want to use xavg
    #x = "x1_ECal"
    
    #slope
    m = (y2-y1)/(x2-x1)
    
    Band_Frame = []
    Band_X = []
    Band_Y = []
    Band_YX = []
    
    for i in range(NumDetectors):
    
        Filter_Condition = f"( (cebraE{i}_ECal - {y1}) <= ({m}*({X_Variable}-{x1}) + {Width} )) && ( (cebraE{i}_ECal - {y1}) >= ({m}*({X_Variable}-{x1}) - {Width}) )"
        
        Band_Frame.append(Reduced_Frame.Filter(cebraTimeGate[i]).Filter(Filter_Condition))

        Band_X.append(Band_Frame[i].Histo1D((f"{X_Variable}_{Band_Name}_TCut{i}",f"{X_Variable}_{Band_Name}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]), X_Variable))
        Band_Y.append(Band_Frame[i].Histo1D((f"cebraE{i}_ECal_{Band_Name}_TCut",f"cebraE{i}_ECal_{Band_Name}_TCut",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))
        Band_YX.append(Band_Frame[i].Histo2D((f"{X_Variable}_cebraE{i}_ECal_{Band_Name}_TCut{i}",f"{X_Variable}_cebraE{i}_ECal_{Band_Name}_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), X_Variable,f"cebraE{i}_ECal"))
        
        if i==0:
            X_Sum_Plot = Band_Frame[i].Histo1D((f"Array_{X_Variable}_{Band_Name}",f"Array_{X_Variable}_{Band_Name}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),X_Variable)
            Y_Sum_Plot = Band_Frame[i].Histo1D((f"Array_cebraE_{Band_Name}",f"Array_cebraE_{Band_Name}",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal")
            YX_Sum_Plot = Band_Frame[i].Histo2D((f"Array_cebraE_{X_Variable}_{Band_Name}",f"Array_cebraE_{X_Variable}_{Band_Name}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2],ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]), X_Variable,f"cebraE{i}_ECal")

        if i!=0:
            X_Sum_Plot.Add(Band_X[i].GetPtr())
            Y_Sum_Plot.Add(Band_Y[i].GetPtr())
            YX_Sum_Plot.Add(Band_YX[i].GetPtr())
        
    return X_Sum_Plot, Y_Sum_Plot, YX_Sum_Plot

# def Excited_States(NumDetectors,Excited_State_Name,FP_Variable,low_fp_gate,high_fp_gate):
    
#     PositionGate = f"{low_fp_gate} < {FP_Variable} && {FP_Variable} < {high_fp_gate}"
    
#     Excited_State_Frame = []
#     det_i_x1 = []
#     det_i_cebraE = []
#     for i in range(NumDetectors):
#         Excited_State_Frame.append(det_i_Reduced_Frame[i].Filter(PositionGate))
#         det_i_x1.append(Excited_State_Frame[j].Histo1D((f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),"x1_ECal"))
#         det_i_cebraE.append(Excited_State_Frame[j].Histo1D((f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))
    
    
    
#     Excited_State_Frame = []
#     det_i_x1 = []
#     det_i_cebraE = []
#     for j in range(0,len(x1_E_Gates)):
        
#         Excited_State_Frame.append(det_i_Reduced_Frame[i].Filter(x1_E_Gates[j][1]))
#         det_i_x1.append(Excited_State_Frame[j].Histo1D((f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"x1_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",ECal_Focal_Histo[0],ECal_Focal_Histo[1],ECal_Focal_Histo[2]),"x1_ECal"))
#         det_i_cebraE.append(Excited_State_Frame[j].Histo1D((f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",f"cebraE_{x1_E_Gates[j][0]}_keV_Protons_TCut{i}",ECal_Gamma_Histo[0],ECal_Gamma_Histo[1],ECal_Gamma_Histo[2]),f"cebraE{i}_ECal"))
        
#     x1_E_List_Protons_TCut.append(det_i_x1)
#     cebraE_E_List_Protons_TCut.append(det_i_cebraE)
#     det_i_Excited_State_Frame.append(Excited_State_Frame)

# For Protons

def saveplots(NumDetectors,Particle,ReactionName):
    
    output.cd()
    
    output.mkdir("Common")
    
    output.mkdir(f"{ReactionName}")
    output.mkdir(f"{ReactionName}/Time_Correlation/")
    output.mkdir(f"{ReactionName}/Particle_Gamma/")
    output.mkdir(f"{ReactionName}/Individual_Spectrums/")
    output.mkdir(f"{ReactionName}/Energy_Calibrated/")
    
    output.cd(f"{ReactionName}")
    
    PID_particle.Write()
    xavg_particle.Write()
    x1_only1plane_particle.Write()
    x1_bothplanes_particle.Write()
    x2_bothplanes_particle.Write()
    x1_theta_particle.Write()
    xavg_theta_particle.Write()
    
    output.cd(f"{ReactionName}/Energy_Calibrated/")
    focalplane_particle_ECal.Write()
    focalplane_ECal_theta_particle.Write()
    
    for i in range(NumDetectors):

        output.cd(f"{ReactionName}/Time_Correlation/")
        cebraTime_toScint_particle[i].Write()
        cebraTime_toScint_particle_TCut[i].Write()
        
        output.cd(f"{ReactionName}/Individual_Spectrums/")
        cebraE_particle_TCut[i].Write()
        focalplane_particle_TCut[i].Write()
        
        output.cd(f"{ReactionName}/Particle_Gamma/")
        focalplane_cebraE_particle_TCut[i].Write()
        
    # Summed Array Plots
    output.cd(f"{ReactionName}/")
    Array_focalplane_particle_TCut.Write()
    Array_cebraE_particle_TCut.Write()
    Array_focalplane_cebraE_particle_TCut.Write() 
    
    #Summed Array Energy Calibrated Plots
    output.cd(f"{ReactionName}/Energy_Calibrated/")
    Array_focalplane_ECal_particle.Write()
    Array_cebraE_ECal_particle.Write()
    Array_focalplane_cebraE_ECal_particle.Write()
    
    #Excited States Plots
    if Particle=="Protons":
        output.mkdir(f"{ReactionName}/50Ti_Excited_States/")
        output.cd("49Ti_dpg/50Ti_Excited_States/")
        for j in range(0,NumExcited):
            Array_x1_E_List[j].Write()
            Array_cebraE_E_List[j].Write()
    
    #Band Plots
    if Particle=="Protons":        
        output.mkdir("49Ti_dpg/GS_Band/")
        output.cd("49Ti_dpg/GS_Band/")
        Array_x1_GS_Band.Write()
        Array_cebraE_GS_Band.Write()
        Array_cebraE_x1_GS_Band.Write()
        
        output.mkdir("49Ti_dpg/2_Plus_Band/")
        output.cd("49Ti_dpg/2_Plus_Band/")
        Array_x1_2_Plus_Band.Write()
        Array_cebraE_2_Plus_Band.Write()
        Array_cebraE_x1_2_Plus_Band.Write()
        
        output.mkdir("49Ti_dpg/4_Plus_Band/")
        output.cd("49Ti_dpg/4_Plus_Band/")
        Array_x1_4_Plus_Band.Write()
        Array_cebraE_4_Plus_Band.Write()
        Array_cebraE_x1_4_Plus_Band.Write()
        
        output.mkdir("49Ti_dpg/6_Plus_Band/")
        output.cd("49Ti_dpg/6_Plus_Band/")
        Array_x1_6_Plus_Band.Write()
        Array_cebraE_6_Plus_Band.Write()
        Array_cebraE_x1_6_Plus_Band.Write()
        
    output.Close()

ROOT.ROOT.EnableImplicitMT(5)
######################################

#Input
#Either "Protons" or "Deuterons" for particle
#Either "x1" (Protons) or "xavg" (Deuterons) for good_focalplane

#ddg
# particle = "Deuterons"
# good_focalplane = "xavg"
# reaction = "49Ti_ddg"

#dpg
particle = "Protons"
good_focalplane = "x1"
reaction = "49Ti_dpg"

######################################
if particle == "Protons":
    cebraTimeGate = [cebraTime_toScint_Gate(0,-678,-670),\
        cebraTime_toScint_Gate(1,-676,-669),\
        cebraTime_toScint_Gate(2,-676,-669),\
        cebraTime_toScint_Gate(3,-676,-668),\
        cebraTime_toScint_Gate(4,-645,-638)]
    
if particle == "Deuterons":
    cebraTimeGate = [cebraTime_toScint_Gate(0,-733,-720),\
        cebraTime_toScint_Gate(1,-733,-717),\
        cebraTime_toScint_Gate(2,-731,-717),\
        cebraTime_toScint_Gate(3,-730,-718),\
        cebraTime_toScint_Gate(4,-700,-688)]

for RUNNO in Ti49_RunList:
        
    print("Current Run: ",RUNNO)
    start = timer()
    
    PATH = f"/home/alconley/Projects/CeBrA_Analysis/49Ti_Analyzed_{particle}_Reduced/"
    
    OUT = f"49Ti_{particle}_Histograms/run_{RUNNO}_{particle}_histo.root"
    output = ROOT.TFile.Open(OUT, "RECREATE")
    
    if particle == "Protons":
        Reduced_Frame = ROOT.RDataFrame(treeName, PATH + f"run_{RUNNO}_protons.root")\
            .Define(f"{good_focalplane}_ECal", f"{ECal_Focal_Para[0]}*{good_focalplane}*{good_focalplane} + {ECal_Focal_Para[1]}*{good_focalplane} + {ECal_Focal_Para[2]}")\
            .Define("cebraE0_ECal", f"{ECal_CeBrA_Para[0]}*cebraE0_GM*cebraE0_GM + {ECal_CeBrA_Para[1]}*cebraE0_GM + {ECal_CeBrA_Para[2]}")\
            .Define("cebraE1_ECal", f"{ECal_CeBrA_Para[0]}*cebraE1_GM*cebraE1_GM + {ECal_CeBrA_Para[1]}*cebraE1_GM + {ECal_CeBrA_Para[2]}")\
            .Define("cebraE2_ECal", f"{ECal_CeBrA_Para[0]}*cebraE2_GM*cebraE2_GM + {ECal_CeBrA_Para[1]}*cebraE2_GM + {ECal_CeBrA_Para[2]}")\
            .Define("cebraE3_ECal", f"{ECal_CeBrA_Para[0]}*cebraE3_GM*cebraE3_GM + {ECal_CeBrA_Para[1]}*cebraE3_GM + {ECal_CeBrA_Para[2]}")\
            .Define("cebraE4_ECal", f"{ECal_CeBrA_Para[0]}*cebraE4_GM*cebraE4_GM + {ECal_CeBrA_Para[1]}*cebraE4_GM + {ECal_CeBrA_Para[2]}")
            
    if particle == "Deuterons":
        Reduced_Frame = ROOT.RDataFrame(treeName, PATH + f"run_{RUNNO}_deuterons.root")\
            .Define(f"{good_focalplane}_ECal", f"{ECal_Focal_Para[0]}*{good_focalplane}*{good_focalplane} + {ECal_Focal_Para[1]}*{good_focalplane} + {ECal_Focal_Para[2]}")\
            .Define("cebraE0_ECal", f"{ECal_CeBrA_Para[0]}*cebraE0_GM*cebraE0_GM + {ECal_CeBrA_Para[1]}*cebraE0_GM + {ECal_CeBrA_Para[2]}")\
            .Define("cebraE1_ECal", f"{ECal_CeBrA_Para[0]}*cebraE1_GM*cebraE1_GM + {ECal_CeBrA_Para[1]}*cebraE1_GM + {ECal_CeBrA_Para[2]}")\
            .Define("cebraE2_ECal", f"{ECal_CeBrA_Para[0]}*cebraE2_GM*cebraE2_GM + {ECal_CeBrA_Para[1]}*cebraE2_GM + {ECal_CeBrA_Para[2]}")\
            .Define("cebraE3_ECal", f"{ECal_CeBrA_Para[0]}*cebraE3_GM*cebraE3_GM + {ECal_CeBrA_Para[1]}*cebraE3_GM + {ECal_CeBrA_Para[2]}")\
            .Define("cebraE4_ECal", f"{ECal_CeBrA_Para[0]}*cebraE4_GM*cebraE4_GM + {ECal_CeBrA_Para[1]}*cebraE4_GM + {ECal_CeBrA_Para[2]}")
    
    #Basic Plots
    PID_particle, xavg_particle, x1_bothplanes_particle, x2_bothplanes_particle, x1_only1plane_particle, x1_theta_particle, xavg_theta_particle,\
        focalplane_particle_ECal, focalplane_ECal_theta_particle,cebraTime_toScint_particle = analysis_basic(particle,good_focalplane)
        
    #Plots with time cuts, coincidence matricies, and summed array
    
    if particle == "Protons":
        cebraTime_toScint_particle_TCut, cebraE_particle_TCut, focalplane_particle_TCut, focalplane_cebraE_particle_TCut,\
            focalplane_ECal_particle_TCut, cebraE_ECal_particle_TCut,focalplane_cebraE_ECal_particle_TCut,\
                Array_focalplane_particle_TCut,Array_cebraE_particle_TCut, Array_focalplane_cebraE_particle_TCut,\
                    Array_focalplane_ECal_particle,Array_cebraE_ECal_particle,Array_focalplane_cebraE_ECal_particle,\
                        Array_x1_E_List,Array_cebraE_E_List,NumExcited = analysis(5,particle,good_focalplane)
    else:
        cebraTime_toScint_particle_TCut, cebraE_particle_TCut, focalplane_particle_TCut, focalplane_cebraE_particle_TCut,\
            focalplane_ECal_particle_TCut, cebraE_ECal_particle_TCut,focalplane_cebraE_ECal_particle_TCut,\
                Array_focalplane_particle_TCut,Array_cebraE_particle_TCut, Array_focalplane_cebraE_particle_TCut,\
                    Array_focalplane_ECal_particle,Array_cebraE_ECal_particle,Array_focalplane_cebraE_ECal_particle = analysis(5,particle,good_focalplane)
    
    if particle == "Protons":
        #Band plots
        Array_x1_GS_Band, Array_cebraE_GS_Band, Array_cebraE_x1_GS_Band = bands(5,"GS_Band",f"{good_focalplane}_ECal",0,0,1553,1553,50)
        Array_x1_2_Plus_Band, Array_cebraE_2_Plus_Band, Array_cebraE_x1_2_Plus_Band = bands(5,"2_Plus_Band",f"{good_focalplane}_ECal",1553,0,4172,2618,50)
        Array_x1_4_Plus_Band, Array_cebraE_4_Plus_Band, Array_cebraE_x1_4_Plus_Band = bands(5,"4_Plus_Band",f"{good_focalplane}_ECal",2674,0,4147,1472,50)
        Array_x1_6_Plus_Band, Array_cebraE_6_Plus_Band, Array_cebraE_x1_6_Plus_Band = bands(5,"6_Plus_Band",f"{good_focalplane}_ECal",3198,0,6123,2924.9,50)
    
    # if particle == "Protons":
    #     #Band plots
    #     Array_x1_GS_Band, Array_cebraE_GS_Band, Array_cebraE_x1_GS_Band = bands(5,"GS_Band",f"{good_focalplane}_ECal",0,0,1553,1000,50)
    #     Array_x1_2_Plus_Band, Array_cebraE_2_Plus_Band, Array_cebraE_x1_2_Plus_Band = bands(5,"2_Plus_Band",f"{good_focalplane}_ECal",1553,0,4172,1700,50)
    #     Array_x1_4_Plus_Band, Array_cebraE_4_Plus_Band, Array_cebraE_x1_4_Plus_Band = bands(5,"4_Plus_Band",f"{good_focalplane}_ECal",2674,0,4880,1427,50)
    #     Array_x1_6_Plus_Band, Array_cebraE_6_Plus_Band, Array_cebraE_x1_6_Plus_Band = bands(5,"6_Plus_Band",f"{good_focalplane}_ECal",3198,0,6123,1770,50)
        
    
    saveplots(5,particle,reaction)
    
    end = timer()
    print("Total Elapsed time for current run: ",end - start)
    # input("Hit Enter when done")