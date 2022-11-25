import ROOT
from timeit import default_timer as timer

#Do Not have it multi-thread if using snapshot 

PATH = "/media/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/49Ti_analyzed/"
treeName = "SPSTree"

#multiple runs (right now it looks at where you have all the analyzed trees at, a good idea in order to get your PID Cut)
Globbed_Raw_Frame = ROOT.RDataFrame(treeName, f"{PATH}*.root")

#Can show what columns are in the tree
#You could just extract what the columns you want and not reduce the data using the PID too
print(Globbed_Raw_Frame.GetColumnNames())

#Shows the Particle Identification
PID = Globbed_Raw_Frame.Filter("scintLeftTime != -1 && anodeBackTime != -1").Histo2D(("ScintL_AnodeBack","PID",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
# c_PID = ROOT.TCanvas()
# PID.Draw("colz") 

#Declare your cut in a string form
Protons = "scintLeft > 800 && scintLeft < 1600 && anodeBack > 50 && anodeBack < 800"

#if you want to make sure the cut looks good
# PID_Protons = Globbed_Raw_Frame.Filter(Protons).Histo2D(("ScintL_AnodeBack_Protons","Protons",2048,0,2048,4096,0,4096), "scintLeft","anodeBack")
# c_PID_Protons = ROOT.TCanvas()
# PID_Protons.Draw("colz") 

input("Hit Enter when done")


Ti49_RunList = [184,185,186,187,188,190,191,192,193,\
    194,195,196,197,198,199,200,201,202,203,204,205,\
        206,207,208,209,210,211,212,213,214,215,216,\
            217,218,219,220,221,222,223,224,225,226,\
                227,228,229,230,231,232]

# for RUNNO in Ti49_RunList:
        
#     print("Current Run: ",RUNNO)
#     start = timer()
    
#     #SPSTree for each run declared here as "Raw_Frame"
#     Raw_Frame = ROOT.RDataFrame(treeName, PATH + f"run_{RUNNO}.root")

#     #Reduce the raw frame with just the particle of choice, in this case protons
#     #You can also declare new columns here too if you want
#     Proton_Frame = Raw_Frame.Filter(Protons)
           
#     #The output file location and name                                                           
#     OUT = f"./run_{RUNNO}_protons.root"
    
#     # Recreating the tree with only the data that I want using Snapshot
#     Proton_Frame.Snapshot("SPSTree",OUT,{"anodeBack","scintLeft","scintLeftTime","cathode","xavg","x1","x2","theta",\
#         "cebraE0","cebraE1","cebraE2","cebraE3","cebraE4","cebraTime0","cebraTime1","cebraTime2","cebraTime3","cebraTime4",\
#             "cebraTime0_toScint","cebraTime1_toScint","cebraTime2_toScint","cebraTime3_toScint","cebraTime4_toScint",\
#                 "cebraE0_GM","cebraE1_GM","cebraE2_GM","cebraE3_GM","cebraE4_GM"})
    
#     end = timer()
#     print("Total Elapsed time for current run: ",end - start)
    