#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cmath>

void deltaE() {
    TFile *file1 = TFile::Open("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output__noAbs/tracksummary_ambims.root");
    TTree *tree1 = (TTree*)file1->Get("tracksummary");

    TFile *file2 = TFile::Open("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output__noAbs/tracksummary_merged.root");
    TTree *tree2 = (TTree*)file2->Get("tracksummary");

    int n_entries = std::min(tree1->GetEntries(), tree2->GetEntries());

    TH1F *hDelta = new TH1F("hDelta", "E_{before Abs} - E_{after Abs} (GeV)", 4*35, 0.5, 3);
    TH1F *hE1 = new TH1F("hE1", "Energia", 100, 0, 100);
    TH1F *hE2 = new TH1F("hE2", "Energia", 100, 0, 100);
    TH1F *hQoP1 = new TH1F("hQoP1", "QoP 1", 100, -0.5, 0.5);
    TH1F *hQoP2 = new TH1F("hQoP2", "QoP 2", 100, -0.5, 0.5);
    TH1F *hTheta1 = new TH1F("hTheta1", "Theta 1", 100, 0, 0.3);
    TH1F *hTheta2 = new TH1F("hTheta2", "Theta 2", 100, 0, 0.3);
    TH1F *hPhi1 = new TH1F("hPhi1", "Phi 1", 100, -3.14, 3.14);
    TH1F *hPhi2 = new TH1F("hPhi2", "Phi 2", 100, -3.14, 3.14);
    TH1F *hLoc01 = new TH1F("hLoc01", "Loc0 1", 100, -600, 600);
    TH1F *hLoc02 = new TH1F("hLoc02", "Loc0 2", 100, -600, 600);
    TH1F *hLoc11 = new TH1F("hLoc11", "Loc1 1", 100, -600, 600);
    TH1F *hLoc12 = new TH1F("hLoc12", "Loc1 2", 100, -600, 600);


    TH1F *herrQoP1 = new TH1F("herrQoP1", "QoP 1", 100, 0, 0.5);
    TH1F *herrQoP2 = new TH1F("herrQoP2", "QoP 2", 100, 0, 0.5);
    TH1F *hDeltaErr = new TH1F("hDeltaErr", "Error in Delta E", 100, 0, 0.4);
    TH1F *hDeltaErrOvE = new TH1F("hDeltaErrOvE", "Error in Delta E", 100, 0, 1);

    float mmuon = 0.1056583745;

    // PUNTATORI A VETTORI
    std::vector<float> *eQOP1 = nullptr;
    std::vector<float> *eQOP2 = nullptr;
    std::vector<float> *eLoc01 = nullptr;
    std::vector<float> *eLoc02 = nullptr;
    std::vector<float> *eLoc11 = nullptr;
    std::vector<float> *eLoc12 = nullptr;
    std::vector<float> *eTheta1 = nullptr;
    std::vector<float> *eTheta2 = nullptr;
    std::vector<float> *ePhi1 = nullptr;
    std::vector<float> *ePhi2 = nullptr;

    std::vector<float> *err_eQOP1 = nullptr;
    std::vector<float> *err_eQOP2 = nullptr;
    std::vector<float> *err_eLoc01 = nullptr;
    std::vector<float> *err_eLoc02 = nullptr;
    std::vector<float> *err_eLoc11 = nullptr;
    std::vector<float> *err_eLoc12 = nullptr;
    std::vector<float> *err_eTheta1 = nullptr;
    std::vector<float> *err_eTheta2 = nullptr;
    std::vector<float> *err_ePhi1 = nullptr;
    std::vector<float> *err_ePhi2 = nullptr;

    std::vector<float> *pid1 = nullptr;
    std::vector<float> *pid2 = nullptr;

    std::vector<float> *t_p = nullptr;
    std::vector<float> *t_eta = nullptr;



    tree1->SetBranchAddress("t_p", &t_p);
    tree1->SetBranchAddress("t_eta", &t_eta);

    tree1->SetBranchAddress("eQOP_fit", &eQOP1);
    tree2->SetBranchAddress("eQOP_fit", &eQOP2);
    tree1->SetBranchAddress("eTHETA_fit", &eTheta1);
    tree2->SetBranchAddress("eTHETA_fit", &eTheta2);
    tree1->SetBranchAddress("ePHI_fit", &ePhi1);
    tree2->SetBranchAddress("ePHI_fit", &ePhi2);
    tree1->SetBranchAddress("eLOC0_fit", &eLoc01);
    tree2->SetBranchAddress("eLOC0_fit", &eLoc02);
    tree1->SetBranchAddress("eLOC1_fit", &eLoc11);
    tree2->SetBranchAddress("eLOC1_fit", &eLoc12);

    tree1->SetBranchAddress("err_eQOP_fit", &err_eQOP1);
    tree2->SetBranchAddress("err_eQOP_fit", &err_eQOP2);
    tree1->SetBranchAddress("err_eTHETA_fit", &err_eTheta1);
    tree2->SetBranchAddress("err_eTHETA_fit", &err_eTheta2);
    tree1->SetBranchAddress("err_ePHI_fit", &err_ePhi1);
    tree2->SetBranchAddress("err_ePHI_fit", &err_ePhi2);
    tree1->SetBranchAddress("err_eLOC0_fit", &err_eLoc01);
    tree2->SetBranchAddress("err_eLOC0_fit", &err_eLoc02);
    tree1->SetBranchAddress("err_eLOC1_fit", &err_eLoc11);
    tree2->SetBranchAddress("err_eLOC1_fit", &err_eLoc12);


    tree1->SetBranchAddress("majorityParticleId", &pid1);
    tree2->SetBranchAddress("majorityParticleId", &pid2);

    for (int i = 0; i < n_entries; ++i) {
        std::cout<< "Processing entry " << i << std::endl;
        tree1->GetEntry(i);
        tree2->GetEntry(i);
        if (pid1->size() != pid2->size()) {
            continue;
        }
        for (size_t j = 0; j < eQOP1->size() && j < eQOP2->size(); ++j) {
            if( std::sqrt((*t_p)[j]*(*t_p)[j]+mmuon*mmuon) > 20 || std::sqrt((*t_p)[j]*(*t_p)[j]+mmuon*mmuon) < 10 || (*t_eta)[j] > 3.2 || (*t_eta)[j] < 2) {
                continue; // Skip tracks with low momentum or outside eta range
            }
            float invQop1 = 1.0 / (*eQOP1)[j];
            float invQop2 = 1.0 / (*eQOP2)[j];
            float E1 = std::sqrt(invQop1 * invQop1 + mmuon * mmuon);
            float E2 = std::sqrt(invQop2 * invQop2 + mmuon * mmuon);

            hQoP1->Fill((*eQOP1)[j]);
            hQoP2->Fill((*eQOP2)[j]);
            hTheta1->Fill((*eTheta1)[j]);
            hTheta2->Fill((*eTheta2)[j]);
            hPhi1->Fill((*ePhi1)[j]);
            hPhi2->Fill((*ePhi2)[j]);
            hLoc01->Fill((*eLoc01)[j]);
            hLoc02->Fill((*eLoc02)[j]);
            hLoc11->Fill((*eLoc11)[j]);
            hLoc12->Fill((*eLoc12)[j]);
            float errP1  = invQop1*invQop1*(*err_eQOP1)[j];
            float errP2  = invQop2*invQop2*(*err_eQOP2)[j];
            herrQoP1->Fill(errP1);
            herrQoP2->Fill(errP2);
            hDeltaErr->Fill(std::sqrt(errP1*errP1 - errP2*errP2));

            float delta = E1 - E2;
            hDeltaErrOvE->Fill(std::sqrt(errP1*errP1 + errP2*errP2) / delta);
            std::cout << "Entry " << i << ", Track " << j 
                      << ": E1 = " << E1 << " GeV, E2 = " << E2 
                      << " GeV, Delta = " << delta << " GeV" << std::endl;
            hE1->Fill(E1);
            hE2->Fill(E2);
            hDelta->Fill(delta);
        }
    }

    gStyle->SetOptStat(11111111);
    TCanvas *cv = new TCanvas("cv", "Differenza di energia", 1600, 1200);
    hDelta->GetXaxis()->SetTitle("#Delta E (GeV)");
    hDelta->GetYaxis()->SetTitle("Entries");
    hDelta->SetLineColor(kBlue);
    hDelta->Draw();
    cv->SaveAs("delta_energy_comparison.png");

    hE1->GetXaxis()->SetTitle("Energia (GeV)");
    hE1->GetYaxis()->SetTitle("Entries");
    hE1->SetLineColor(kRed);
    hE2->SetLineColor(kBlue);

    hE1->Draw();
    hE2->Draw("SAME");
    cv->SaveAs("energy_distribution.png");

    hQoP1->GetXaxis()->SetTitle("QoP");
    hQoP1->GetYaxis()->SetTitle("Entries");
    hQoP1->SetLineColor(kRed);
    hQoP2->SetLineColor(kBlue);

    hQoP1->Draw();
    hQoP2->Draw("SAME");
    cv->SaveAs("QoP_distribution.png");
    hTheta1->GetXaxis()->SetTitle("Theta");
    hTheta1->GetYaxis()->SetTitle("Entries");
    hTheta1->SetLineColor(kRed);
    hTheta2->SetLineColor(kBlue);
    hTheta1->Draw();
    hTheta2->Draw("SAME");
    cv->SaveAs("theta_distribution.png");
    hPhi1->GetXaxis()->SetTitle("Phi");
    hPhi1->GetYaxis()->SetTitle("Entries");
    hPhi1->SetLineColor(kRed);
    hPhi2->SetLineColor(kBlue);
    hPhi1->Draw();
    hPhi2->Draw("SAME");
    cv->SaveAs("phi_distribution.png");
    hLoc01->GetXaxis()->SetTitle("Loc0");
    hLoc01->GetYaxis()->SetTitle("Entries");
    hLoc01->SetLineColor(kRed);
    hLoc02->SetLineColor(kBlue);
    hLoc01->Draw();
    hLoc02->Draw("SAME");
    cv->SaveAs("loc0_distribution.png");
    hLoc11->GetXaxis()->SetTitle("Loc1");
    hLoc11->GetYaxis()->SetTitle("Entries");
    hLoc11->SetLineColor(kRed);
    hLoc12->SetLineColor(kBlue);
    hLoc11->Draw();
    hLoc12->Draw("SAME");
    cv->SaveAs("loc1_distribution.png");

    herrQoP1->GetXaxis()->SetTitle("Error QoP 1");
    herrQoP1->GetYaxis()->SetTitle("Entries");
    herrQoP1->SetLineColor(kRed);
    herrQoP2->SetLineColor(kBlue);
    herrQoP1->Draw();
    herrQoP2->Draw("SAME");
    cv->SaveAs("error_QoP_distribution.png");

    hDeltaErr->GetXaxis()->SetTitle("Error in Delta E");
    hDeltaErr->GetYaxis()->SetTitle("Entries");
    hDeltaErr->Draw();
    cv->SaveAs("error_delta_energy_distribution.png");
    hDeltaErrOvE->GetXaxis()->SetTitle("Error in Delta E / Delta E");
    hDeltaErrOvE->GetYaxis()->SetTitle("Entries");
    hDeltaErrOvE->Draw();
    cv->SaveAs("error_delta_energy_over_delta_energy_distribution.png");
}
