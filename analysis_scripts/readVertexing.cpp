#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TMath.h>
#include <TH1.h>
#include <TCanvas.h>


void readVertexing(const char *filename = "../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.8/performance_vertexing.root", 
              const char *outputname = "vertexing",
              const char *pngname = "vertexing")
{
    gROOT->SetBatch(kTRUE); // Turn off graphics

    // written by ChatGPT

    //---------------------------------------
    // open files
    //---------------------------------------
    // Open the ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Get the TTree from the file
    TTree *tree = (TTree *)file->Get("vertexing");

    //---------------------------------------
    // create histos
    //---------------------------------------
    TFile *fout = new TFile(outputname, "recreate");

    TH1D *hresX = new TH1D("hresX", "; x_{PV}^{gen}-x_{PV}^{rec} (mm);entries", 200, -0.06, 0.06);
    TH1D *hresY = new TH1D("hresY", "; y_{PV}^{gen}-y_{PV}^{rec} (mm);entries", 200, -0.06, 0.06);
    TH1D *hresZ = new TH1D("hresZ", "; z_{PV}^{gen}-z_{PV}^{rec} (mm);entries", 200, -0.3, 0.3);
    TH1D *hposX = new TH1D("hposX", "; x_{PV}^{gen} (mm);entries", 200, -0.6, 0.6);
    TH1D *hposY = new TH1D("hposY", "; y_{PV}^{gen} (mm);entries", 200, -0.6, 0.6);
    TH1D *hpos1Z = new TH1D("hpos1Z", "; z_{PV}^{gen} (mm);entries", 200, -2, 2);
    TH1D *hpos2Z = new TH1D("hpos2Z", "; z_{PV}^{gen} (mm);entries", 200, -2-0*(12+1.5), 2-0*(12+1.5));
    TH1D *hpos3Z = new TH1D("hpos3Z", "; z_{PV}^{gen} (mm);entries", 200, -2-0*(12+1.5)*2, 2-0*(12+1.5)*2);
    TH1D *hpos4Z = new TH1D("hpos4Z", "; z_{PV}^{gen} (mm);entries", 200, -2-0*(12+1.5)*3, 2-0*(12+1.5)*3);
    TH1D *hpos5Z = new TH1D("hpos5Z", "; z_{PV}^{gen} (mm);entries", 200, -2-0*(12+1.5)*4, 2-0*(12+1.5)*4);
    TH1D *hpullX = new TH1D("hpullX", "; pull x;entries", 200, -5, 5);
    TH1D *hpullY = new TH1D("hpullY", "; pull y;entries", 200, -5, 5);
    TH1D *hpullZ = new TH1D("hpullZ", "; pull z;entries", 200, -5, 5);
    TH1D *hNTrk = new TH1D("hNTrk", "; N_{trk}^{contr};entries", 200, -0.5, 999.5);
    TH2D *hNTrkVsZ = new TH2D("hNTrkVsZ", "; N_{trk}^{contr};z_{PV}^{rec};entries", 200, -0.5, 999.5, 130, -60,5);
    TH2D *hNTrkvsNTrkTruth = new TH2D("hNTrkvsNTrkTruth", "; N_{trk}^{contr}; N_{trk}^{truth};entries", 200, -0.5, 999.5, 200, -0.5, 999.5);
    TH1D *hNTrkTruth = new TH1D("hNTrkTruth", "; N_{trk}^{truth};entries", 200, -0.5, 999.5);
    TH1D *hDeltaNTrk = new TH1D("hDeltaNTrk", "; N_{trk}^{truth}-N_{trk}^{contr};entries", 400, -1000, 1000);



    TH2D *hresXvsZ = new TH2D("hresXvsZ", "; x_{PV}^{gen}-x_{PV}^{rec} (mm);z;entries", 200, -0.2, 0.2, 130, -60,5);
    TH2D *hresYvsZ = new TH2D("hresYvsZ", "; y_{PV}^{gen}-y_{PV}^{rec} (mm);z;entries", 200, -0.2, 0.2, 130, -60,5);
    TH2D *hresZvsZ = new TH2D("hresZvsZ", "; z_{PV}^{gen}-z_{PV}^{rec} (mm);z;entries", 200, -1, 1, 130, -60,5);
    TH2D *hpullXvsZ = new TH2D("hpullXvsZ", "; pull x;z;entries", 200, -100, 100, 130, -60,5);
    TH2D *hpullYvsZ = new TH2D("hpullYvsZ", "; pull y;z;entries", 200, -100, 100, 130, -60,5);
    TH2D *hpullZvsZ = new TH2D("hpullZvsZ", "; pull z;z;entries", 200, -100, 100, 130, -60,5);

    //---------------------------------------
    // read variables
    //---------------------------------------
    // Define variables to hold kinematic data of true particles
    std::vector<float> *t_resx = new std::vector<float>;
    std::vector<float> *t_resy = nullptr;
    std::vector<float> *t_resz = nullptr;
    std::vector<float> *t_x = nullptr;
    std::vector<float> *t_y = nullptr;
    std::vector<float> *t_z = nullptr;
    std::vector<float> *t_pullx = nullptr;
    std::vector<float> *t_pully = nullptr;
    std::vector<float> *t_pullz = nullptr;
    std::vector<float> *t_nTrk = nullptr;
    std::vector<float> *t_nTrkTruth = nullptr;
    // Set branch addresses for truth particle
    tree->SetBranchAddress("resX", &t_resx);
    tree->SetBranchAddress("resY", &t_resy);
    tree->SetBranchAddress("resZ", &t_resz);

    tree->SetBranchAddress("recoX", &t_y);
    tree->SetBranchAddress("recoY", &t_x);
    tree->SetBranchAddress("recoZ", &t_z);

    tree->SetBranchAddress("pullX", &t_pullx);
    tree->SetBranchAddress("pullY", &t_pully);
    tree->SetBranchAddress("pullZ", &t_pullz);

    tree->SetBranchAddress("nTracksRecoVtx", &t_nTrk);
    tree->SetBranchAddress("nTracksTruthVtx", &t_nTrkTruth);

    //---------------------------------------
    // loop over entries
    //---------------------------------------
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i)
    {
        tree->GetEntry(i);
        if (t_resx->size() < 1)
            continue;

        int sizz = (int)t_resx->size();
        for (int j = 0; j < sizz; j++)
        {
            hposX->Fill(t_x->at(j));
            hposY->Fill(t_y->at(j));
            hresX->Fill(t_resx->at(j));
            hresY->Fill(t_resy->at(j));
            hresZ->Fill(t_resz->at(j));
            hpullX->Fill(t_pullx->at(j));
            hpullY->Fill(t_pully->at(j));
            hpullZ->Fill(t_pullz->at(j));
            hpos1Z->Fill(t_z->at(j)+0.75);
            hpos2Z->Fill(t_z->at(j)+0.75+(12+1.5));
            hpos3Z->Fill(t_z->at(j)+0.75+(12+1.5)*2);
            hpos4Z->Fill(t_z->at(j)+0.75+(12+1.5)*3);
            hpos5Z->Fill(t_z->at(j)+0.75+(12+1.5)*4);
            hNTrk->Fill(t_nTrk->at(j));
            hNTrkTruth->Fill(t_nTrkTruth->at(j));
            hDeltaNTrk->Fill(t_nTrkTruth->at(j)-t_nTrk->at(j));
            hNTrkVsZ->Fill(t_nTrk->at(j), t_z->at(j));
            hNTrkvsNTrkTruth->Fill(t_nTrk->at(j), t_nTrkTruth->at(j));

            hresXvsZ->Fill(t_resx->at(j), t_z->at(j));
            hresYvsZ->Fill(t_resy->at(j), t_z->at(j));
            hresZvsZ->Fill(t_resz->at(j), t_z->at(j));
            hpullXvsZ->Fill(t_pullx->at(j), t_z->at(j));
            hpullYvsZ->Fill(t_pully->at(j), t_z->at(j));
            hpullZvsZ->Fill(t_pullz->at(j), t_z->at(j));
        }
    }

    //--------------------------------------
    // plot histos
    //--------------------------------------
    TCanvas *c0 = new TCanvas("c0", "",20,20, 3000, 2000);

    c0->Divide(3,2);
    c0->cd(1);
    hresX->Draw();
    c0->cd(2);
    hresY->Draw();
    c0->cd(3);
    hresZ->Draw();
    c0->cd(4);
    hpullX->Draw();
    c0->cd(5);
    hpullY->Draw();
    c0->cd(6);
    hpullZ->Draw();
    c0->SaveAs(pngname);

    TCanvas *c1 = new TCanvas("c1", "",20,20, 3000, 1000);
    c1->Divide(3,1);
    c1->cd(1);
    hNTrk->Draw();
    c1->cd(2);
    hNTrkTruth->Draw();
    c1->cd(3);
    hDeltaNTrk->Draw();

    hposX->Write();
    hposY->Write();
    hpos1Z->Write();
    hpos2Z->Write();
    hpos3Z->Write();
    hpos4Z->Write();
    hpos5Z->Write();
    hresX->Write();
    hresY->Write();
    hresZ->Write();
    hpullX->Write();
    hpullY->Write();
    hpullZ->Write();

    hNTrk->Write();
    hNTrkTruth->Write();
    hDeltaNTrk->Write();
    
    hresXvsZ->Write();
    hresYvsZ->Write();
    hresZvsZ->Write();
    hpullXvsZ->Write();
    hpullYvsZ->Write();
    hpullZvsZ->Write();
    
    hNTrkVsZ->Write();
    hNTrkvsNTrkTruth->Write();
    c0->Write();
    c1->Write();
    // Close the file
    file->Close();
}


int runAll(const char *rootfile = "results/summary.root"){
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.01/performance_vertexing_merged.root", "results/vtx_0.01_merged.root", "results/vtx_0.01_merged.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.02/performance_vertexing_merged.root", "results/vtx_0.02_merged.root", "results/vtx_0.02_merged.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.05/performance_vertexing_merged.root", "results/vtx_0.05_merged.root", "results/vtx_0.05_merged.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.1/performance_vertexing_merged.root", "results/vtx_0.1_merged.root", "results/vtx_0.1_merged.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.2/performance_vertexing_merged.root", "results/vtx_0.2_merged.root", "results/vtx_0.2_merged.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.4/performance_vertexing_merged.root", "results/vtx_0.4_merged.root", "results/vtx_0.4_merged.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.6/performance_vertexing_merged.root", "results/vtx_0.6_merged.root", "results/vtx_0.6_merged.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.8/performance_vertexing_merged.root", "results/vtx_0.8_merged.root", "results/vtx_0.8_merged.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1/performance_vertexing_merged.root", "results/vtx_central_merged.root", "results/vtx_central_merged.png");

    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.01/performance_vertexing_truth.root", "results/vtx_0.01_truth.root", "results/vtx_0.01_truth.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.02/performance_vertexing_truth.root", "results/vtx_0.02_truth.root", "results/vtx_0.02_truth.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.05/performance_vertexing_truth.root", "results/vtx_0.05_truth.root", "results/vtx_0.05_truth.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.1/performance_vertexing_truth.root", "results/vtx_0.1_truth.root", "results/vtx_0.1_truth.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.2/performance_vertexing_truth.root", "results/vtx_0.2_truth.root", "results/vtx_0.2_truth.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.4/performance_vertexing_truth.root", "results/vtx_0.4_truth.root", "results/vtx_0.4_truth.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.6/performance_vertexing_truth.root", "results/vtx_0.6_truth.root", "results/vtx_0.6_truth.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.8/performance_vertexing_truth.root", "results/vtx_0.8_truth.root", "results/vtx_0.8_truth.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1/performance_vertexing_truth.root", "results/vtx_central_truth.root", "results/vtx_central_truth.png");

    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.01/performance_vertexing_step1.root", "results/vtx_0.01_step1.root", "results/vtx_0.01_step1.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.02/performance_vertexing_step1.root", "results/vtx_0.02_step1.root", "results/vtx_0.02_step1.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.05/performance_vertexing_step1.root", "results/vtx_0.05_step1.root", "results/vtx_0.05_step1.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.1/performance_vertexing_step1.root", "results/vtx_0.1_step1.root", "results/vtx_0.1_step1.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.2/performance_vertexing_step1.root", "results/vtx_0.2_step1.root", "results/vtx_0.2_step1.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.4/performance_vertexing_step1.root", "results/vtx_0.4_step1.root", "results/vtx_0.4_step1.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.6/performance_vertexing_step1.root", "results/vtx_0.6_step1.root", "results/vtx_0.6_step1.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.8/performance_vertexing_step1.root", "results/vtx_0.8_step1.root", "results/vtx_0.8_step1.png");
    readVertexing("../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1/performance_vertexing_step1.root", "results/vtx_central_step1.root", "results/vtx_central_step1.png");
    
    const char *outputname[] = {
        "results/vtx_0.01_merged.root",
        "results/vtx_0.02_merged.root",
        "results/vtx_0.05_merged.root",
        "results/vtx_0.1_merged.root",
        "results/vtx_0.2_merged.root",
        "results/vtx_0.4_merged.root",
        "results/vtx_0.6_merged.root",
        "results/vtx_0.8_merged.root",
        "results/vtx_central_merged.root"};

    const char *outputnametruth[] = {
        "results/vtx_0.01_truth.root",
        "results/vtx_0.02_truth.root",
        "results/vtx_0.05_truth.root",
        "results/vtx_0.1_truth.root",
        "results/vtx_0.2_truth.root",
        "results/vtx_0.4_truth.root",
        "results/vtx_0.6_truth.root",
        "results/vtx_0.8_truth.root",
        "results/vtx_central_truth.root"};

    const char *outputnamestep1[] = {
        "results/vtx_0.01_step1.root",
        "results/vtx_0.02_step1.root",
        "results/vtx_0.05_step1.root",
        "results/vtx_0.1_step1.root",
        "results/vtx_0.2_step1.root",
        "results/vtx_0.4_step1.root",
        "results/vtx_0.6_step1.root",
        "results/vtx_0.8_step1.root",
        "results/vtx_central_step1.root"};

    float multiplicity[] = {
        0.01,
        0.02,
        0.05,
        0.1,
        0.2,
        0.4,
        0.6,
        0.8,
        1};


    int numElements = sizeof(multiplicity) / sizeof(multiplicity[0]);
    TGraphErrors* sigma_x_graph = new TGraphErrors(numElements);
    sigma_x_graph->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; #sigma_{x} (#mu m)");
    
    TGraphErrors* sigma_y_graph = new TGraphErrors(numElements);
    sigma_y_graph->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; #sigma_{y} (#mu m)");

    TGraphErrors* sigma_z_graph = new TGraphErrors(numElements);
    sigma_z_graph->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; #sigma_{z} (#mu m)");

    TGraphErrors* ntrk_graph = new TGraphErrors(numElements);
    ntrk_graph->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; ntrk (#mu m)");

    TGraphErrors* eff_graph = new TGraphErrors(numElements);
    eff_graph->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; vtx efficiency (\%)");
    ///////////////
    TGraphErrors* sigma_x_graph_truth = new TGraphErrors(numElements);
    sigma_x_graph_truth->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; #sigma_{x} (#mu m)");
    
    TGraphErrors* sigma_y_graph_truth = new TGraphErrors(numElements);
    sigma_y_graph_truth->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; #sigma_{y} (#mu m)");

    TGraphErrors* sigma_z_graph_truth = new TGraphErrors(numElements);
    sigma_z_graph_truth->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; #sigma_{z} (#mu m)");

    TGraphErrors* ntrk_graph_truth = new TGraphErrors(numElements);
    ntrk_graph_truth->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; ntrk (#mu m)");

    TGraphErrors* eff_graph_truth = new TGraphErrors(numElements);
    eff_graph_truth->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; vtx efficiency (\%)");
    ///////////////
    TGraphErrors* sigma_x_graph_step1 = new TGraphErrors(numElements);
    sigma_x_graph_step1->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; #sigma_{x} (#mu m)");
    
    TGraphErrors* sigma_y_graph_step1 = new TGraphErrors(numElements);
    sigma_y_graph_step1->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; #sigma_{y} (#mu m)");

    TGraphErrors* sigma_z_graph_step1 = new TGraphErrors(numElements);
    sigma_z_graph_step1->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; #sigma_{z} (#mu m)");

    TGraphErrors* ntrk_graph_step1 = new TGraphErrors(numElements);
    ntrk_graph_step1->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; ntrk (#mu m)");

    TGraphErrors* eff_graph_step1 = new TGraphErrors(numElements);
    eff_graph_step1->SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5\%>; vtx efficiency (\%)");


    sigma_x_graph->SetName("sigma_x");
    sigma_y_graph->SetName("sigma_y");
    sigma_z_graph->SetName("sigma_z");
    ntrk_graph->SetName("ntrk");
    eff_graph->SetName("eff");

    sigma_x_graph_truth->SetName("sigma_x_truth");
    sigma_y_graph_truth->SetName("sigma_y_truth");
    sigma_z_graph_truth->SetName("sigma_z_truth");
    ntrk_graph_truth->SetName("ntrk_truth");
    eff_graph_truth->SetName("eff_truth");

    sigma_x_graph_step1->SetName("sigma_x_step1");
    sigma_y_graph_step1->SetName("sigma_y_step1");
    sigma_z_graph_step1->SetName("sigma_z_step1");
    ntrk_graph_step1->SetName("ntrk_step1");
    eff_graph_step1->SetName("eff_step1");

    for(int t_i=0; t_i<numElements; t_i++){
        TFile* vtx = new TFile(outputname[t_i], "read");
        TH1D* hresX = (TH1D*) vtx->Get("hresX");
        TH1D* hresY = (TH1D*) vtx->Get("hresY");
        TH1D* hresZ = (TH1D*) vtx->Get("hresZ");
        TH1D* hntrk = (TH1D*) vtx->Get("hNTrk");
        TF1* gaus = new TF1("gaus","gausn(0)");
        //gaus->SetParameter(1,0);
        //gaus->SetParLimits(1,-0.1,0.1);
        //gaus->SetParameter(2,0.009);
        //gaus->SetParLimits(1,0.003,0.009*10);

        float ntrk = hntrk->GetMean();
        hresX->Fit(gaus);
        float sigmaX = gaus->GetParameter(2)*1000.;
        float errsigmaX = gaus->GetParError(2)*1000.;
        hresY->Fit(gaus);
        float sigmaY = gaus->GetParameter(2)*1000.;
        float errsigmaY = gaus->GetParError(2)*1000.;
        hresZ->Fit(gaus);
        float sigmaZ = gaus->GetParameter(2)*1000.;
        float errsigmaZ = gaus->GetParError(2)*1000.;

        sigma_x_graph->SetPoint(t_i, multiplicity[t_i], sigmaX);
        sigma_x_graph->SetPointError(t_i, 0, 0);

        sigma_y_graph->SetPoint(t_i, multiplicity[t_i], sigmaY);
        sigma_y_graph->SetPointError(t_i, 0, errsigmaY);

        sigma_z_graph->SetPoint(t_i, multiplicity[t_i], sigmaZ);
        sigma_z_graph->SetPointError(t_i, 0, errsigmaZ);

        ntrk_graph->SetPoint(t_i, multiplicity[t_i], ntrk);
        ntrk_graph->SetPointError(t_i, 0, 0);

        eff_graph->SetPoint(t_i, multiplicity[t_i], hntrk->GetEntries()/1000.);
        eff_graph->SetPointError(t_i, 0, 0);

        TFile* vtxtruth = new TFile(outputnametruth[t_i], "read");
        TH1D* hresXtruth = (TH1D*) vtxtruth->Get("hresX");
        TH1D* hresYtruth = (TH1D*) vtxtruth->Get("hresY");
        TH1D* hresZtruth = (TH1D*) vtxtruth->Get("hresZ");
        TH1D* hntrktruth = (TH1D*) vtxtruth->Get("hNTrk");
        float ntrktruth = hntrktruth->GetMean();
        hresXtruth->Fit(gaus);
        float sigmaXtruth = gaus->GetParameter(2)*1000.;
        float errsigmaXtruth = gaus->GetParError(2)*1000.;
        hresYtruth->Fit(gaus);
        float sigmaYtruth = gaus->GetParameter(2)*1000.;
        float errsigmaYtruth = gaus->GetParError(2)*1000.;
        hresZtruth->Fit(gaus);
        float sigmaZtruth = gaus->GetParameter(2)*1000.;
        float errsigmaZtruth = gaus->GetParError(2)*1000.;

        sigma_x_graph_truth->SetPoint(t_i, multiplicity[t_i], sigmaXtruth);
        sigma_x_graph_truth->SetPointError(t_i, 0, 0);

        sigma_y_graph_truth->SetPoint(t_i, multiplicity[t_i], sigmaYtruth);
        sigma_y_graph_truth->SetPointError(t_i, 0, errsigmaYtruth);

        sigma_z_graph_truth->SetPoint(t_i, multiplicity[t_i], sigmaZtruth);
        sigma_z_graph_truth->SetPointError(t_i, 0, errsigmaZtruth);

        ntrk_graph_truth->SetPoint(t_i, multiplicity[t_i], ntrktruth);
        ntrk_graph_truth->SetPointError(t_i, 0, 0);

        eff_graph_truth->SetPoint(t_i, multiplicity[t_i], hntrktruth->GetEntries()/1000.);
        eff_graph_truth->SetPointError(t_i, 0, 0);


        TFile* vtxstep1 = new TFile(outputnamestep1[t_i], "read");
        TH1D* hresXstep1 = (TH1D*) vtxstep1->Get("hresX");
        TH1D* hresYstep1 = (TH1D*) vtxstep1->Get("hresY");
        TH1D* hresZstep1 = (TH1D*) vtxstep1->Get("hresZ");
        TH1D* hntrkstep1 = (TH1D*) vtxstep1->Get("hNTrk");
        float ntrkstep1 = hntrkstep1->GetMean();
        hresXstep1->Fit(gaus);
        float sigmaXstep1 = gaus->GetParameter(2)*1000.;
        float errsigmaXstep1 = gaus->GetParError(2)*1000.;
        hresYstep1->Fit(gaus);
        float sigmaYstep1 = gaus->GetParameter(2)*1000.;
        float errsigmaYstep1 = gaus->GetParError(2)*1000.;
        hresZstep1->Fit(gaus);
        float sigmaZstep1 = gaus->GetParameter(2)*1000.;
        float errsigmaZstep1 = gaus->GetParError(2)*1000.;

        sigma_x_graph_step1->SetPoint(t_i, multiplicity[t_i], sigmaXstep1);
        sigma_x_graph_step1->SetPointError(t_i, 0, 0);

        sigma_y_graph_step1->SetPoint(t_i, multiplicity[t_i], sigmaYstep1);
        sigma_y_graph_step1->SetPointError(t_i, 0, errsigmaYstep1);

        sigma_z_graph_step1->SetPoint(t_i, multiplicity[t_i], sigmaZstep1);
        sigma_z_graph_step1->SetPointError(t_i, 0, errsigmaZstep1);

        ntrk_graph_step1->SetPoint(t_i, multiplicity[t_i], ntrkstep1);
        ntrk_graph_step1->SetPointError(t_i, 0, 0);

        eff_graph_step1->SetPoint(t_i, multiplicity[t_i], hntrkstep1->GetEntries()/1000.);
        eff_graph_step1->SetPointError(t_i, 0, 0);

    }

    TFile* summary = new TFile(rootfile,"recreate");
    sigma_x_graph->SetMarkerStyle(20);
    sigma_y_graph->SetMarkerStyle(20);
    sigma_z_graph->SetMarkerStyle(20);
    ntrk_graph->SetMarkerStyle(20);
    eff_graph->SetMarkerStyle(20);

    sigma_x_graph->Write();
    sigma_y_graph->Write();
    sigma_z_graph->Write();
    ntrk_graph->Write();
    eff_graph->Write();

    sigma_x_graph_truth->SetMarkerStyle(20);
    sigma_y_graph_truth->SetMarkerStyle(20);
    sigma_z_graph_truth->SetMarkerStyle(20);
    ntrk_graph_truth->SetMarkerStyle(20);
    eff_graph_truth->SetMarkerStyle(20);

    sigma_x_graph_truth->SetMarkerColor(kRed);
    sigma_y_graph_truth->SetMarkerColor(kRed);
    sigma_z_graph_truth->SetMarkerColor(kRed);
    ntrk_graph_truth->SetMarkerColor(kRed);
    eff_graph_truth->SetMarkerColor(kRed);
    sigma_x_graph_truth->SetLineColor(kRed);
    sigma_y_graph_truth->SetLineColor(kRed);
    sigma_z_graph_truth->SetLineColor(kRed);
    ntrk_graph_truth->SetLineColor(kRed);
    eff_graph_truth->SetLineColor(kRed);

    sigma_x_graph_truth->Write();
    sigma_y_graph_truth->Write();
    sigma_z_graph_truth->Write();
    ntrk_graph_truth->Write();
    eff_graph_truth->Write();

    sigma_x_graph_step1->SetMarkerStyle(20);
    sigma_y_graph_step1->SetMarkerStyle(20);
    sigma_z_graph_step1->SetMarkerStyle(20);
    ntrk_graph_step1->SetMarkerStyle(20);
    eff_graph_step1->SetMarkerStyle(20);

    sigma_x_graph_step1->SetMarkerColor(kBlue);
    sigma_y_graph_step1->SetMarkerColor(kBlue);
    sigma_z_graph_step1->SetMarkerColor(kBlue);
    ntrk_graph_step1->SetMarkerColor(kBlue);
    eff_graph_step1->SetMarkerColor(kBlue);

    sigma_x_graph_step1->SetLineColor(kBlue);
    sigma_y_graph_step1->SetLineColor(kBlue);
    sigma_z_graph_step1->SetLineColor(kBlue);
    ntrk_graph_step1->SetLineColor(kBlue);
    eff_graph_step1->SetLineColor(kBlue);

    sigma_x_graph_step1->Write();
    sigma_y_graph_step1->Write();
    sigma_z_graph_step1->Write();
    ntrk_graph_step1->Write();
    eff_graph_step1->Write();


    auto legend = new TLegend(0.8,0.8,1,1);
    legend->AddEntry(sigma_x_graph,"Tracks from step 1 & 2","lep");
    legend->AddEntry(sigma_x_graph_truth,"Truth info","lep");
    legend->AddEntry(sigma_x_graph_step1,"Tracks from step1","lep");
    
    TCanvas* cv_sigmax = new TCanvas("cv_sigmax","cv_sigmax",1600,1200);
    sigma_x_graph->Draw("ALP");
    //sigma_x_graph_truth->Draw("same ALP");
    //sigma_x_graph_step1->Draw("same ALP");
    //legend->Draw("same");
    cv_sigmax->SaveAs("sigmax.png");

    TCanvas* cv_sigmay = new TCanvas("cv_sigmay","cv_sigmay",1600,1200);    
    sigma_y_graph->Draw();
    // sigma_y_graph_truth->Draw("same");
    //sigma_y_graph_step1->Draw("same");
    //legend->Draw("same");
    cv_sigmay->SaveAs("sigmay.png");

    TCanvas* cv_sigmaz = new TCanvas("cv_sigmaz","cv_sigmaz",1600,1200);

    TF1* na60zres = new TF1("na60zres","[0]/(sqrt(x*[1]+1))",0,1);
    na60zres->SetParameter(0, 1050);
    na60zres->SetParameter(1, ntrk_graph->GetPointY(numElements-1));
    na60zres->Draw();
    sigma_z_graph->Draw("same");
    //sigma_z_graph_truth->Draw("same");
    //sigma_z_graph_step1->Draw("same");
    //legend->Draw("same");
    cv_sigmaz->SaveAs("sigmaz.png");

    TCanvas* cv_ntrk = new TCanvas("cv_ntrk","cv_ntrk",1600,1200);
    ntrk_graph->Draw();
    ntrk_graph_truth->Draw("same");
    ntrk_graph_step1->Draw("same");
    //legend->Draw("same");
    cv_ntrk->SaveAs("ntrk.png");

    TCanvas* cv_eff = new TCanvas("cv_eff","cv_eff",1600,1200);
    eff_graph->Draw();
    eff_graph_truth->Draw("apl same");
    eff_graph_step1->Draw("apl same");
    //legend->Draw("same");
    cv_eff->SaveAs("eff.png");

    summary->Close();

    return 0;
}