#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TMath.h>
#include <TH1.h>
#include <TCanvas.h>

void checkMatchingEfficiency(
    const char *filenameVT = "../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.8/tracksummary_ambi.root",
    const char *filenameMS = "../output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.8/tracksummary_ambi.root")
{
    gROOT->SetBatch(kTRUE); // Turn off graphics

    // written by ChatGPT

    //---------------------------------------
    // open files
    //---------------------------------------
    // Open the ROOT file
    TFile *fileVT = TFile::Open(filenameVT);
    // Open the ROOT file
    TFile *fileMS = TFile::Open(filenameMS);

    // Get the TTree from the file
    TTree *treeVT = (TTree *)fileVT->Get("tracksummary");
    // Get the TTree from the file
    TTree *treeMS = (TTree *)fileMS->Get("tracksummary");

    //---------------------------------------
    // create histos
    //---------------------------------------

    TH1D *hloc0VT = new TH1D("hloc0VT", "; loc0;entries", 200, -5, 5);
    TH1D *hloc1VT = new TH1D("hloc1VT", "; loc0;entries", 200, -80, 15);
    TH1D *hphiVT = new TH1D("hphiVT", "; loc0;entries", 200, -TMath::Pi(), TMath::Pi());
    TH1D *hthetaVT = new TH1D("hthetaVT", "; loc0;entries", 200, 0, 5);
    TH1D *hqopVT = new TH1D("hqopVT", "; loc0;entries", 200, -10, 10);

    TH1D *hloc0MS = new TH1D("hloc0MS", "; loc0;entries", 200, -5, 5);
    TH1D *hloc1MS = new TH1D("hloc1MS", "; loc0;entries", 200, -80, 15);
    TH1D *hphiMS = new TH1D("hphiMS", "; loc0;entries", 200, -TMath::Pi(), TMath::Pi());
    TH1D *hthetaMS = new TH1D("hthetaMS", "; loc0;entries", 200, 0, 5);
    TH1D *hqopMS = new TH1D("hqopMS", "; loc0;entries", 200, -10, 10);

    TH1D *hdloc0 = new TH1D("hdloc0", "; loc0;entries", 200, -1, 1);
    TH1D *hdloc1 = new TH1D("hdloc1", "; loc0;entries", 200, -1, 1);
    TH1D *hdphi = new TH1D("hdphi", "; loc0;entries", 200, -1, 1);
    TH1D *hdtheta = new TH1D("hdtheta", "; loc0;entries", 200, -1, 1);
    TH1D *hdqop = new TH1D("hdqop", "; loc0;entries", 200, -1, 1);
    TH1D *hdsum = new TH1D("hdsum", "; loc0;entries", 200, 0, 100);


    TH1D *hdloc0Match = new TH1D("hdloc0Match", "; loc0;entries", 200, -1, 1);
    TH1D *hdloc1Match = new TH1D("hdloc1Match", "; loc0;entries", 200, -1, 1);
    TH1D *hdphiMatch = new TH1D("hdphiMatch", "; loc0;entries", 200, -1, 1);
    TH1D *hdthetaMatch = new TH1D("hdthetaMatch", "; loc0;entries", 200, -1, 1);
    TH1D *hdqopMatch = new TH1D("hdqopMatch", "; loc0;entries", 200, -1, 1);
    TH1D *hdsumMatch = new TH1D("hdsumMatch", "; loc0;entries", 200, 0, 100);


    TH1D *hdloc0Fake = new TH1D("hdloc0Fake", "; loc0;entries", 200, -1, 1);
    TH1D *hdloc1Fake = new TH1D("hdloc1Fake", "; loc0;entries", 200, -1, 1);
    TH1D *hdphiFake = new TH1D("hdphiFake", "; loc0;entries", 200, -1, 1);
    TH1D *hdthetaFake = new TH1D("hdthetaFake", "; loc0;entries", 200, -1, 1);
    TH1D *hdqopFake = new TH1D("hdqopFake", "; loc0;entries", 200, -1, 1);
    TH1D *hdsumFake = new TH1D("hdsumFake", "; loc0;entries", 200, 0, 100);


    TH1D *hndloc0 = new TH1D("hndloc0", "; loc0;entries", 200, -3, 3);
    TH1D *hndloc1 = new TH1D("hndloc1", "; loc0;entries", 200, -3, 3);
    TH1D *hndphi = new TH1D("hndphi", "; loc0;entries", 200, -3, 3);
    TH1D *hndtheta = new TH1D("hndtheta", "; loc0;entries", 200, -3, 3);
    TH1D *hndqop = new TH1D("hndqop", "; loc0;entries", 200, -3, 3);
    TH1D *hndsum = new TH1D("hndsum", "; loc0;entries", 200, 0, 20);

    TH1D *hndloc0Match = new TH1D("hndloc0Match", "; loc0;entries", 200, -3, 3);
    TH1D *hndloc1Match = new TH1D("hndloc1Match", "; loc0;entries", 200, -3, 3);
    TH1D *hndphiMatch = new TH1D("hndphiMatch", "; loc0;entries", 200, -3, 3);
    TH1D *hndthetaMatch = new TH1D("hndthetaMatch", "; loc0;entries", 200, -3, 3);
    TH1D *hndqopMatch = new TH1D("hndqopMatch", "; loc0;entries", 200, -3, 3);
    TH1D *hndsumMatch = new TH1D("hndsumMatch", "; loc0;entries", 200, 0, 20);

    TH1D *hndloc0Fake = new TH1D("hndloc0Fake", "; loc0;entries", 200, -3, 3);
    TH1D *hndloc1Fake = new TH1D("hndloc1Fake", "; loc0;entries", 200, -3, 3);
    TH1D *hndphiFake = new TH1D("hndphiFake", "; loc0;entries", 200, -3, 3);
    TH1D *hndthetaFake = new TH1D("hndthetaFake", "; loc0;entries", 200, -3, 3);
    TH1D *hndqopFake = new TH1D("hndqopFake", "; loc0;entries", 200, -3, 3);
    TH1D *hndsumFake = new TH1D("hndsumFake", "; loc0;entries", 200, 0, 20);

    TH1D *matchingEff = new TH1D("matchingEff",";;Efficiency",2,-0.5,1.5);
    TH1D *deltasum = new TH1D("deltasum", "; loc0;entries", 1000, -100, 100);
    //---------------------------------------
    // read variables
    //---------------------------------------
    // Define variables to hold kinematic data of true particles
    std::vector<float> *t_loc0VT = new std::vector<float>;
    std::vector<float> *t_loc1VT = nullptr;
    std::vector<float> *t_phiVT = nullptr;
    std::vector<float> *t_thetaVT = nullptr;
    std::vector<float> *t_qopVT = nullptr;
    std::vector<float> *t_eloc0VT = nullptr;
    std::vector<float> *t_eloc1VT = nullptr;
    std::vector<float> *t_ephiVT = nullptr;
    std::vector<float> *t_ethetaVT = nullptr;
    std::vector<float> *t_eqopVT = nullptr;
    std::vector<float> *t_idVT = nullptr;
    // Set branch addresses for truth particle
    treeVT->SetBranchAddress("eLOC0_fit", &t_loc0VT);
    treeVT->SetBranchAddress("eLOC1_fit", &t_loc1VT);
    treeVT->SetBranchAddress("ePHI_fit", &t_phiVT);
    treeVT->SetBranchAddress("eTHETA_fit", &t_thetaVT);
    treeVT->SetBranchAddress("eQOP_fit", &t_qopVT);
    treeVT->SetBranchAddress("err_eLOC0_fit", &t_eloc0VT);
    treeVT->SetBranchAddress("err_eLOC1_fit", &t_eloc1VT);
    treeVT->SetBranchAddress("err_ePHI_fit", &t_ephiVT);
    treeVT->SetBranchAddress("err_eTHETA_fit", &t_ethetaVT);
    treeVT->SetBranchAddress("err_eQOP_fit", &t_eqopVT);
    treeVT->SetBranchAddress("majorityParticleID", &t_idVT);

    // Define variables to hold kinematic data of true particles
    std::vector<float> *t_loc0MS = new std::vector<float>;
    std::vector<float> *t_loc1MS = nullptr;
    std::vector<float> *t_phiMS = nullptr;
    std::vector<float> *t_thetaMS = nullptr;
    std::vector<float> *t_qopMS = nullptr;
    std::vector<float> *t_eloc0MS = nullptr;
    std::vector<float> *t_eloc1MS = nullptr;
    std::vector<float> *t_ephiMS = nullptr;
    std::vector<float> *t_ethetaMS = nullptr;
    std::vector<float> *t_eqopMS = nullptr;
    std::vector<float> *t_idMS = nullptr;
    // Set branch addresses for truth particle
    treeMS->SetBranchAddress("eLOC0_fit", &t_loc0MS);
    treeMS->SetBranchAddress("eLOC1_fit", &t_loc1MS);
    treeMS->SetBranchAddress("ePHI_fit", &t_phiMS);
    treeMS->SetBranchAddress("eTHETA_fit", &t_thetaMS);
    treeMS->SetBranchAddress("eQOP_fit", &t_qopMS);
    treeMS->SetBranchAddress("err_eLOC0_fit", &t_eloc0MS);
    treeMS->SetBranchAddress("err_eLOC1_fit", &t_eloc1MS);
    treeMS->SetBranchAddress("err_ePHI_fit", &t_ephiMS);
    treeMS->SetBranchAddress("err_eTHETA_fit", &t_ethetaMS);
    treeMS->SetBranchAddress("err_eQOP_fit", &t_eqopMS);
    treeMS->SetBranchAddress("majorityParticleID", &t_idMS);

    //---------------------------------------
    // loop over entries
    //---------------------------------------
    Long64_t nentriesVT = treeVT->GetEntries();
    Long64_t nentriesMS = treeMS->GetEntries();
    //loop over the event
    for (Long64_t i = 0; i < nentriesVT; ++i)
    {
        treeVT->GetEntry(i);
        treeMS->GetEntry(i);
        if (t_loc0VT->size() < 1)
            continue;

        int sizzVT = (int)t_loc0VT->size();
        int sizzMS = (int)t_loc0MS->size();

        //loop over the VT tracks
        for (int j = 0; j < sizzVT; j++)
        {
            hloc0VT->Fill(t_loc0VT->at(j));
            hloc1VT->Fill(t_loc1VT->at(j));
            hphiVT->Fill(t_phiVT->at(j));
            hthetaVT->Fill(t_thetaVT->at(j));
            hqopVT->Fill(t_qopVT->at(j));
            Long64_t min_sum = 100000000000;
            Long64_t match_sum = 100000000000;

            //loop over the MS tracks
            for (int k = 0; k < sizzMS; k++)
            {
                hloc0MS->Fill(t_loc0MS->at(j));
                hloc1MS->Fill(t_loc1MS->at(j));
                hphiMS->Fill(t_phiMS->at(j));
                hthetaMS->Fill(t_thetaMS->at(j));
                hqopMS->Fill(t_qopMS->at(j));

                hdloc0->Fill(t_loc0VT->at(j) - t_loc0MS->at(k));
                hdloc1->Fill(t_loc1VT->at(j) - t_loc1MS->at(k));
                hdphi->Fill(t_phiVT->at(j) - t_phiMS->at(k));
                hdtheta->Fill(t_thetaVT->at(j) - t_thetaMS->at(k));
                hdqop->Fill(t_qopVT->at(j) - t_qopMS->at(k));

                float sum = TMath::Sqrt((t_loc0VT->at(j) - t_loc0MS->at(k))*(t_loc0VT->at(j) - t_loc0MS->at(k))/(t_eloc0VT->at(i)*t_eloc0VT->at(i)+t_eloc0VT->at(k)*t_eloc0VT->at(k)) + 
                                        (t_loc1VT->at(j) - t_loc1MS->at(k))*(t_loc1VT->at(j) - t_loc1MS->at(k))/(t_eloc0VT->at(i)*t_eloc0VT->at(i)+t_eloc0VT->at(k)*t_eloc0VT->at(k)) + 
                                        (t_phiVT->at(j) - t_phiMS->at(k))*(t_phiVT->at(j) - t_phiMS->at(k))/(t_eloc0VT->at(i)*t_eloc0VT->at(i)+t_eloc0VT->at(k)*t_eloc0VT->at(k)) + 
                                        (t_thetaVT->at(j) - t_thetaMS->at(k))*(t_thetaVT->at(j) - t_thetaMS->at(k))/(t_eloc0VT->at(i)*t_eloc0VT->at(i)+t_eloc0VT->at(k)*t_eloc0VT->at(k)) + 
                                        (t_qopVT->at(j) - t_qopMS->at(k))*(t_qopVT->at(j) - t_qopMS->at(k))/(t_eloc0VT->at(i)*t_eloc0VT->at(i)+t_eloc0VT->at(k)*t_eloc0VT->at(k)) +);

                hdsum->Fill(sum);

                if (t_idVT->at(j) == t_idMS->at(k))
                {
                    hdloc0Match->Fill(t_loc0VT->at(j) - t_loc0MS->at(k));
                    hdloc1Match->Fill(t_loc1VT->at(j) - t_loc1MS->at(k));
                    hdphiMatch->Fill(t_phiVT->at(j) - t_phiMS->at(k));
                    hdthetaMatch->Fill(t_thetaVT->at(j) - t_thetaMS->at(k));
                    hdqopMatch->Fill(t_qopVT->at(j) - t_qopMS->at(k));
                    hdsumMatch->Fill(sum);
                    match_sum = sum;
                }
                else
                {
                    hdloc0Fake->Fill(t_loc0VT->at(j) - t_loc0MS->at(k));
                    hdloc1Fake->Fill(t_loc1VT->at(j) - t_loc1MS->at(k));
                    hdphiFake->Fill(t_phiVT->at(j) - t_phiMS->at(k));
                    hdthetaFake->Fill(t_thetaVT->at(j) - t_thetaMS->at(k));
                    hdqopFake->Fill(t_qopVT->at(j) - t_qopMS->at(k));
                    hdsumFake->Fill(sum);
                    if(sum<min_sum)
                        min_sum = sum;
                }
            }
            if(min_sum<match_sum)
                matchingEff->Fill(0);
            else
                matchingEff->Fill(1);

            deltasum->Fill(match_sum-match_sum);
        }
    }
    matchingEff->Scale(1./matchingEff->GetEntries());
    /*
    for (Long64_t i = 0; i < nentries; ++i)
    {
        tree->GetEntry(i);
        if (t_loc0->size() < 1)
            continue;

        int sizz = (int)t_loc0->size();
        for (int j = 0; j < sizz; j++)
        {
            for (int k = j + 1; k < sizz; k++)
            {
                hndloc0->Fill((t_loc0->at(j) - t_loc0->at(k) - hdloc0->GetMean()) / hdloc0->GetStdDev());
                hndloc1->Fill((t_loc1->at(j) - t_loc1->at(k) - hdloc1->GetMean()) / hdloc1->GetStdDev());
                hndphi->Fill((t_phi->at(j) - t_phi->at(k) - hdphi->GetMean()) / hdphi->GetStdDev());
                hndtheta->Fill((t_theta->at(j) - t_theta->at(k) - hdtheta->GetMean()) / hdtheta->GetStdDev());
                hndqop->Fill((t_qop->at(j) - t_qop->at(k) - hdqop->GetMean()) / hdqop->GetStdDev());
                hndsum->Fill(
                    TMath::Sqrt((t_loc0->at(j) - t_loc0->at(k) - hdloc0->GetMean()) / hdloc0->GetStdDev()) +
                    TMath::Sqrt((t_loc1->at(j) - t_loc1->at(k) - hdloc1->GetMean()) / hdloc1->GetStdDev()) +
                    TMath::Sqrt((t_phi->at(j) - t_phi->at(k) - hdphi->GetMean()) / hdphi->GetStdDev()) +
                    TMath::Sqrt((t_theta->at(j) - t_theta->at(k) - hdtheta->GetMean()) / hdtheta->GetStdDev()) +
                    TMath::Sqrt((t_qop->at(j) - t_qop->at(k) - hdqop->GetMean()) / hdqop->GetStdDev()));
            }
        }
    }
    */
    TFile *res = new TFile("matching.root", "recreate");
    hloc0VT->Write();
    hloc1VT->Write();
    hphiVT->Write();
    hthetaVT->Write();
    hqopVT->Write();

    hloc0MS->Write();
    hloc1MS->Write();
    hphiMS->Write();
    hthetaMS->Write();
    hqopMS->Write();


    hdloc1->Write();
    hdphi->Write();
    hdtheta->Write();
    hdqop->Write();
    hdsum->Write();

    hdloc1Match->Write();
    hdphiMatch->Write();
    hdthetaMatch->Write();
    hdqopMatch->Write();
    hdsumMatch->Write();

    hdloc1Fake->Write();
    hdphiFake->Write();
    hdthetaFake->Write();
    hdqopFake->Write();
    hdsumFake->Write();
    

    hndloc1->Write();
    hndphi->Write();
    hndtheta->Write();
    hndqop->Write();
    hndsum->Write();

    hndloc1Match->Write();
    hndphiMatch->Write();
    hndthetaMatch->Write();
    hndqopMatch->Write();
    hndsumMatch->Write();

    hndloc1Fake->Write();
    hndphiFake->Write();
    hndthetaFake->Write();
    hndqopFake->Write();
    hndsumFake->Write();

    matchingEff->Write();
    deltasum->Write();
    res->Close();
}