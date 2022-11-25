from cmath import pi
import ROOT
import math
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

start = timer()

#where the analyzed files/tree is located
PATH = "/home/alconley/Projects/CeBrA_Analysis/49Ti_Analyzed_Protons_Reduced/"

OUT = "1st_excited_state_tilt_correction.root"
output = ROOT.TFile.Open(OUT, "RECREATE")

#number of detectors
n_det = 5

treeName = "SPSTree"

#optimize the number of thread before you run the loop
ROOT.ROOT.EnableImplicitMT(5)

Reduced_Frame = ROOT.RDataFrame(treeName, PATH + "*.root")

                                                                               
x1_Protons = Reduced_Frame.Filter("x1 != -1e6 && x2 != -1e6").Histo1D(("x1_Protons","x1_Protons",600,-300,300),"x1")

x1_theta_Protons = Reduced_Frame.Histo2D(("x1_theta_Protons","x1_theta_Protons",600,-300,300,600,0,pi/2),"x1","theta")

                            
#first excited state correction
correcting_frame = Reduced_Frame.Filter("0.6 < theta && theta < 0.9").Filter("22 < x1 && x1 < 35").AsNumpy(columns=["x1","theta"])

#4 MeV state excited state correction uncalibrated
# correcting_frame = Proton_Frame.Filter("0.6 < theta && theta < 0.9").Filter("-102 < x1 && x1 < 92").AsNumpy(columns=["x1","theta"])

y = correcting_frame["x1"]
x = correcting_frame["theta"]

a, stats = np.polynomial.polynomial.polyfit(x,y,3,full=True)
print(a)
print(stats)

poly = np.polynomial.polynomial.Polynomial(a)
print(poly)

f = lambda th: a[1]*th+a[2]*th**2+a[3]*th**3

plotX = np.linspace(0.4,1.4)
plotY = poly(plotX)

plt.plot(x,y,".",plotX,plotY)

# xPrime = y - f(x)
# plt.plot(x,xPrime,".")
# plt.show()


# add 28.3115 for first excited state correction to line up with uncorrected spectrum
Corrected_Frame = Reduced_Frame.Define("x1_Corrected",f"x1 - ({a[0]}+{a[1]}*theta + {a[2]}*theta*theta + {a[3]}*theta*theta*theta)+28.3115")

x1_Protons_Corrected = Corrected_Frame.Filter("theta > 0.6 && theta < 0.9").Histo1D(("x1_Protons_Corrected","x1_Protons_Corrected",600,-300,300),"x1_Corrected")
x1_theta_Protons_Corrected = Corrected_Frame.Histo2D(("x1_theta_Protons_Corrected","x1_theta_Protons_Corrected",600,-300,300,600,0,pi/2),"x1_Corrected","theta")

# c1 = ROOT.TCanvas()
# x1_Protons.Draw()

# c2 = ROOT.TCanvas()
# x1_theta_Protons.Draw("colz")

output.cd()

x1_Protons.Write()
x1_theta_Protons.Write()

x1_Protons_Corrected.Write()
x1_theta_Protons_Corrected.Write()


output.Close()
plt.show()


end = timer()
print("Total Elapsed time for current run: ",end - start)

input("Hit Enter when done")
