
// #include <ROOT/RDataFrame.hxx>
// #include "../../git_workspace/SPS_CEBRA_EventBuilder/src/spsdict/DataStructs.h"


// #define FMT_HEADER_ONLY
// #include </home/alconley/Programs/fmt-9.0.0/include/fmt/format.h>

// void test()
{
    TStopwatch time;

    auto fileName = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/12C_analyzed/run_82.root";
    // std::string outputfile = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/histograms/60Co_spectrum.root";
    auto treeName = "SPSTree";

    ROOT::EnableImplicitMT();
    ROOT::RDataFrame d(treeName, fileName);

    auto PID =d.Filter("scintLeftTime > -1").Histo2D({"ScintL_AnodeBack","Scintilator vs Back Anode",2048,0,2048,4096,0,4096}, "scintLeft","anodeBack");
    auto protons = "scintLeft > 300 && scintLeft < 900 && anodeBack > 100 && anodeBack < 800";

    auto ProtonFrame = d.Filter(protons).Define("cebraTime_toScint_0","cebraTime0-scintLeftTime");

    auto cebraTime0_toScint = ProtonFrame.Histo1D({"cebraTime_toScint_0","cebraTime_toScint_0",6400,-3200,3200},"cebraTime_toScint_0");
    auto time_0_gate = ("cebraTime_toScint_0 > -692 && cebraTime_toScint_0 < -684");
    auto cebraTime0_toScint_TP = ProtonFrame.Filter(time_0_gate).Histo1D({"cebraTime_toScint_0","cebraTime_toScint_0",6400,-3200,3200},"cebraTime_toScint_0");
    auto cebraE0_TP = ProtonFrame.Filter(time_0_gate).Histo1D({"cebraE0","cebraE0_TP",4096,0,4096},"cebraE0");
    auto x1_TP = ProtonFrame.Filter(time_0_gate).Histo1D({"x1","x1",600,-300,300},"x1");
    auto x1_cebraE0_TP = ProtonFrame.Filter(time_0_gate).Histo2D({"x1_cebraE0_TP","x1_cebraE0_TP",600,-300,300,4096,0,4096}, "x1","cebraE0");
    
    auto PID_Protons =d.Filter(protons).Histo2D({"ScintL_AnodeBack","Scintilator vs Back Anode",2048,0,2048,4096,0,4096}, "scintLeft","anodeBack");
    auto xavg_Cut =d.Filter(protons).Histo1D({"xavg","xavg",600,-300,300},"xavg");
    auto x1_Cut =d.Filter(protons).Histo1D({"x1","x1",600,-300,300},"x1");
    auto x2_Cut =d.Filter(protons).Histo1D({"x2","x2",600,-300,300},"x2");





    // auto time4 = "cebraTime_toScint_4>-660 && cebraTime_toScint_4<-650";
    // auto Cebra4_toScint = d2.Filter(time4).Histo1D({"cebraTime_toScint_4","cebraTime_toScint_4",20,-665,-645},"cebraTime_toScint_4");

    // cebraTime0_toScint_Gate->Draw();
    // cebraE0_TP->Draw();
    x1_cebraE0_TP->Draw("colz");


    // ROOT::RDF::RunGraphs({PID,PID_Protons,xavg_Cut,x1_Cut,x2_Cut});
    // PID->Draw("colz");
    // cebraE3_noCuts->Draw();
    // Cebra4_toScint->Draw();
    // x1_Cut->Draw();
    // Cebra0_toScint_Cut->Draw();

    // auto colNames = d.GetColumnNames();
    // // Print columns' names
    // for (auto &&colName : colNames) std::cout << colName << std::endl;
   
    // auto colType = d.GetColumnType("cebraE0");
    // // Print column type
    // std::cout << "Specified Column is " << colType << std::endl;

    // return 0;

    time.Print();
}
