#include <TROOT.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>

#include "../../git_workspace/SPS_CEBRA_EventBuilder/src/spsdict/DataStructs.h"
#include <string>
#include <iostream>

R__LOAD_LIBRARY(../../git_workspace/SPS_CEBRA_EventBuilder/lib/libSPSDict.so);

void analysis()
{

    std::string inputfile = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/analyzed/run_82.root";
    std::string outputfile = "/home/alconley/SanDisk/Research/SPS_CEBRA_Sept_2022/WorkingDir/histograms/run_82_histograms.root";

    TFile* input = TFile::Open(inputfile.c_str(),"read");
    if(!input->IsOpen())
    {
        std::cerr<<"Blergh"<<std::endl;
        return;
    }

    TTree* tree = (TTree*) input->Get("SPSTree");

    TFile* output = TFile::Open(outputfile.c_str(), "recreate");
    if(!output->IsOpen())
    {
        std::cerr<<"Blergh"<<std::endl;
        return;
    }

    ProcessedEvent* eventHandle = new ProcessedEvent();

    tree->SetBranchAddress("event", &eventHandle);

    uint64_t nentries = tree->GetEntries();
    uint64_t counter = 0;
    uint64_t flushVal = 0.05 * nentries;
    uint64_t flushCount = 0;
    

    TH1F* cebraE0_Shifted = new TH1F("cebraE0_noCuts","cebraE0_noCuts", 4096, 0, 4096);
    TH1F* cebraE1_Shifted = new TH1F("cebraE1_noCuts","cebraE1_noCuts", 4096, 0, 4096);
    TH1F* cebraE2_Shifted = new TH1F("cebraE2_noCuts","cebraE2_noCuts", 4096, 0, 4096);
    TH1F* cebraE3_Shifted = new TH1F("cebraE3_noCuts","cebraE3_noCuts", 4096, 0, 4096);
    TH1F* cebraE4_Shifted = new TH1F("cebraE4_noCuts","cebraE4_noCuts", 4096, 0, 4096);

    double m_0 = 1;
    double y_0 = 0;
  
    double cebraE0;
    double cebraE1;
    double cebraE2;
    double cebraE3;
    double cebraE4;

    double cebraE0_shifted;
    double cebraE1_shifted;
    double cebraE2_shifted;
    double cebraE3_shifted;
    double cebraE4_shifted;

    for(uint64_t i=0; i<nentries; i++)
    {
        tree->GetEntry(i);
        counter++;
        if(counter == flushVal)
        {
            counter =0;
            flushCount++;
            std::cout << "Percent of data processed: " << flushCount * 5 << "%" << std::endl;
        }

        // cebraE4_histo->Fill(eventHandle->cebraE[4]);
        // // bigXavg = 2.0 * eventHandle->xavg;


        cebraE4 = m_0*(eventHandle->cebraE[4])+y_0;
        cebraE4_Shifted->Fill(cebraE4);



    }

    input->Close();
    output->cd();

    
    cebraE4_Shifted->Write();




    output->Close();

    delete input;
    delete output;
    delete eventHandle;
}
