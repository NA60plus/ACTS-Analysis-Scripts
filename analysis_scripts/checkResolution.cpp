#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

gStyle->SetTitleSize(0.05, "XYZ"); // Sets the size for X, Y, and Z axis titles
gStyle->SetTitleFontSize(0.05);    // Sets the size for the histogram title
gStyle->SetLabelSize(0.04, "XYZ");
gStyle->SetPadLeftMargin(0.15);   // Set left margin
gStyle->SetPadRightMargin(0.05);  // Set right margin
gStyle->SetPadTopMargin(0.05);    // Set top margin
gStyle->SetPadBottomMargin(0.15); // Set bottom margin
gStyle->SetHistLineWidth(2);
gStyle->SetOptTitle(0); // Show mean (1), RMS (1), but hide entries (0) and title (0)
// gStyle->SetOptStat(1100); // Show mean (1), RMS (1), but hide entries (0) and title (0)
gStyle->SetCanvasDefW(1800); // Set default canvas width to 800 pixels
gStyle->SetCanvasDefH(1200); // Set default canvas height to 600 pixels

float getPt(float qop, float theta, float phi)
{
    return 1. / std::abs(qop) * std::sin(theta);
}
void checkResolution(const char *filename = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_rubenxprino_40GeV_omega_filtered_ruben/tracksummary_ambi.root",
                     const char *filenameCarbon = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_rubenxprino_40GeV_omega_carbon_filtered_ruben/tracksummary_ambi.root")
{

    TFile *file = TFile::Open(filename);
    TTree *tree = (TTree *)file->Get("tracksummary");

    TH1D *he_resLoc0 = new TH1D("he_resLoc0", ";res d0;Entries", 200, -0.4, 0.4);
    TH1D *he_resLoc1 = new TH1D("he_resLoc1", ";res z0;Entries", 200, -1, 1);
    TH1D *he_resPhi = new TH1D("he_resPhi", ";res phi;Entries", 200, -0.02, 0.02);
    TH1D *he_resTheta = new TH1D("he_resTheta", ";res theta;Entries", 200, -0.1, 0.1);
    TH1D *he_resQOP = new TH1D("he_resQOP", ";res QoP;Entries", 200, -0.1, 0.1);


    TH2D *he_resLoc0_vs_p = new TH2D("he_resLoc0_vs_p", ";#it{p} (GeV/#it{c});res d0;Entries", 200, -0.4, 0.4);
    TH2D *he_resLoc1_vs_p = new TH2D("he_resLoc1_vs_p", ";#it{p} (GeV/#it{c});res z0;Entries", 200, -1, 1);
    TH2D *he_resPhi_vs_p = new TH2D("he_resPhi_vs_p", ";#it{p} (GeV/#it{c});res phi;Entries", 200, -0.02, 0.02);
    TH2D *he_resTheta_vs_p = new TH2D("he_resTheta_vs_p", ";#it{p} (GeV/#it{c});res theta;Entries", 200, -0.1, 0.1);
    TH2D *he_resQOP_vs_p = new TH2D("he_resQOP_vs_p", ";#it{p} (GeV/#it{c});res QoP;Entries", 200, -0.1, 0.1);

    TH2D *he_resLoc0_vs_pt = new TH2D("he_resLoc0_vs_pt", ";#it{p}_{T} (GeV/#it{c});res d0;Entries", 200, -0.4, 0.4);
    TH2D *he_resLoc1_vs_pt = new TH2D("he_resLoc1_vs_pt", ";#it{p}_{T} (GeV/#it{c});res z0;Entries", 200, -1, 1);
    TH2D *he_resPhi_vs_pt = new TH2D("he_resPhi_vs_pt", ";#it{p}_{T} (GeV/#it{c});res phi;Entries", 200, -0.02, 0.02);
    TH2D *he_resTheta_vs_pt = new TH2D("he_resTheta_vs_pt", ";#it{p}_{T} (GeV/#it{c});res theta;Entries", 200, -0.1, 0.1);
    TH2D *he_resQOP_vs_pt = new TH2D("he_resQOP_vs_pt", ";#it{p}_{T} (GeV/#it{c});res QoP;Entries", 200, -0.1, 0.1);

    std::vector<float> *res_eLOC0 = new std::vector<float>;
    std::vector<float> *res_eLOC1 = new std::vector<float>;
    std::vector<float> *res_ePHI = new std::vector<float>;
    std::vector<float> *res_eTHETA = new std::vector<float>;
    std::vector<float> *res_eQOP = new std::vector<float>;

    std::vector<float> *eLOC0 = new std::vector<float>;
    std::vector<float> *eLOC1 = new std::vector<float>;
    std::vector<float> *ePHI = new std::vector<float>;
    std::vector<float> *eTHETA = new std::vector<float>;
    std::vector<float> *eQOP = new std::vector<float>;

    std::vector<float> *t_p = new std::vector<float>;
    std::vector<float> *t_pT = new std::vector<float>;

    tree->SetBranchAddress("res_eLOC0_fit", &res_eLOC0);
    tree->SetBranchAddress("res_eLOC1_fit", &res_eLOC1);
    tree->SetBranchAddress("res_ePHI_fit", &res_ePHI);
    tree->SetBranchAddress("res_eTHETA_fit", &res_eTHETA);
    tree->SetBranchAddress("res_eQOP_fit", &res_eQOP);

    tree->SetBranchAddress("eLOC0_fit", &eLOC0);
    tree->SetBranchAddress("eLOC1_fit", &eLOC1);
    tree->SetBranchAddress("ePHI_fit", &ePHI);
    tree->SetBranchAddress("eTHETA_fit", &eTHETA);
    tree->SetBranchAddress("eQOP_fit", &eQOP);

    tree->SetBranchAddress("t_p", &t_p);
    tree->SetBranchAddress("t_pT", &t_pT);

    Long64_t nentries = tree->GetEntries();

    for (Long64_t i = 0; i < nentries; ++i)
    {
        tree->GetEntry(i);
        for (int im = 0; im < res_eLOC0_->size(); im++)
        {
            auto pt = getPt(t_p->at(im), eTHETA->at(im), ePHI->at(im));
            he_resLoc0_->Fill(res_eLOC0_->at(im));
            he_resLoc1_->Fill(res_eLOC1_->at(im));
            he_resTheta_->Fill(res_ePHI_->at(im));
            he_resPhi_->Fill(res_eTHETA_->at(im));
            he_resQOP_->Fill(res_eQOP_->at(im));

            he_resLoc0_->Fill(res_eLOC0_->at(im));
            he_resLoc1_->Fill(res_eLOC1_->at(im));
            he_resTheta_->Fill(res_ePHI_->at(im));
            he_resPhi_->Fill(res_eTHETA_->at(im));
            he_resQOP_->Fill(res_eQOP_->at(im));

            he_resLoc0_->Fill(res_eLOC0_->at(im));
            he_resLoc1_->Fill(res_eLOC1_->at(im));
            he_resTheta_->Fill(res_ePHI_->at(im));
            he_resPhi_->Fill(res_eTHETA_->at(im));
            he_resQOP_->Fill(res_eQOP_->at(im));
        }
    }

}

void checkResolution(const char *filenameBase = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/tracksummary_ambi.root",
                     const char *filenameCarbon = "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/tracksummary_ambi.root")
{

    TFile *fileBase = TFile::Open(filenameBase);
    TTree *treeBase = (TTree *)fileBase->Get("tracksummary");
    TFile *fileCarbon = TFile::Open(filenameCarbon);
    TTree *treeCarbon = (TTree *)fileCarbon->Get("tracksummary");

    TH1D *he_resLoc0_Base = new TH1D("he_resLoc0_Base", ";res d0;Entries", 200, -0.4, 0.4);
    TH1D *he_resLoc1_Base = new TH1D("he_resLoc1_Base", ";res z0;Entries", 200, -1, 1);
    TH1D *he_resPhi_Base = new TH1D("he_resPhi_Base", ";res phi;Entries", 200, -0.02, 0.02);
    TH1D *he_resTheta_Base = new TH1D("he_resTheta_Base", ";res theta;Entries", 200, -0.1, 0.1);
    TH1D *he_resQOP_Base = new TH1D("he_resQOP_Base", ";res QoP;Entries", 200, -0.1, 0.1);

    std::vector<float> *res_eLOC0_Base = new std::vector<float>;
    std::vector<float> *res_eLOC1_Base = new std::vector<float>;
    std::vector<float> *res_ePHI_Base = new std::vector<float>;
    std::vector<float> *res_eTHETA_Base = new std::vector<float>;
    std::vector<float> *res_eQOP_Base = new std::vector<float>;
    std::vector<float> *eQOP_Base = new std::vector<float>;

    treeBase->SetBranchAddress("res_eLOC0_fit", &res_eLOC0_Base);
    treeBase->SetBranchAddress("res_eLOC1_fit", &res_eLOC1_Base);
    treeBase->SetBranchAddress("res_ePHI_fit", &res_ePHI_Base);
    treeBase->SetBranchAddress("res_eTHETA_fit", &res_eTHETA_Base);
    treeBase->SetBranchAddress("res_eQOP_fit", &res_eQOP_Base);
    treeBase->SetBranchAddress("eQOP_fit", &eQOP_Base);

    TH1D *he_resLoc0_Carbon = new TH1D("he_resLoc0_Carbon", ";res d0;Entries", 200, -0.4, 0.4);
    TH1D *he_resLoc1_Carbon = new TH1D("he_resLoc1_Carbon", ";res z0;Entries", 200, -1, 1);
    TH1D *he_resPhi_Carbon = new TH1D("he_resPhi_Carbon", ";res phi;Entries", 200, -0.02, 0.02);
    TH1D *he_resTheta_Carbon = new TH1D("he_resTheta_Carbon", ";res theta;Entries", 200, -0.1, 0.1);
    TH1D *he_resQOP_Carbon = new TH1D("he_resQOP_Carbon", ";res QoP;Entries", 200, -0.1, 0.1);

    he_resLoc0_Carbon->SetLineColor(kRed);
    he_resLoc1_Carbon->SetLineColor(kRed);
    he_resPhi_Carbon->SetLineColor(kRed);
    he_resTheta_Carbon->SetLineColor(kRed);
    ;
    he_resQOP_Carbon->SetLineColor(kRed);

    std::vector<float> *res_eLOC0_Carbon = new std::vector<float>;
    std::vector<float> *res_eLOC1_Carbon = new std::vector<float>;
    std::vector<float> *res_ePHI_Carbon = new std::vector<float>;
    std::vector<float> *res_eTHETA_Carbon = new std::vector<float>;
    std::vector<float> *res_eQOP_Carbon = new std::vector<float>;
    std::vector<float> *eQOP_Carbon = new std::vector<float>;

    treeCarbon->SetBranchAddress("res_eLOC0_fit", &res_eLOC0_Carbon);
    treeCarbon->SetBranchAddress("res_eLOC1_fit", &res_eLOC1_Carbon);
    treeCarbon->SetBranchAddress("res_ePHI_fit", &res_ePHI_Carbon);
    treeCarbon->SetBranchAddress("res_eTHETA_fit", &res_eTHETA_Carbon);
    treeCarbon->SetBranchAddress("res_eQOP_fit", &res_eQOP_Carbon);
    treeCarbon->SetBranchAddress("eQOP_fit", &eQOP_Carbon);

    Long64_t nentriesBase = treeBase->GetEntries();

    for (Long64_t i = 0; i < nentriesBase; ++i)
    {
        treeBase->GetEntry(i);
        for (int im = 0; im < res_eLOC0_Base->size(); im++)
        {
            // if(TMath::Abs(eQOP_Base->at(im)) < 0.5) continue;
            he_resLoc0_Base->Fill(res_eLOC0_Base->at(im));
            he_resLoc1_Base->Fill(res_eLOC1_Base->at(im));
            he_resTheta_Base->Fill(res_ePHI_Base->at(im));
            he_resPhi_Base->Fill(res_eTHETA_Base->at(im));
            he_resQOP_Base->Fill(res_eQOP_Base->at(im));
        }
    }

    Long64_t nentriesCarbon = treeCarbon->GetEntries();

    for (Long64_t i = 0; i < nentriesCarbon; ++i)
    {
        treeCarbon->GetEntry(i);
        for (int im = 0; im < res_eLOC0_Carbon->size(); im++)
        {
            // if(TMath::Abs(eQOP_Carbon->at(im)) < 0.5) continue;
            he_resLoc0_Carbon->Fill(res_eLOC0_Carbon->at(im));
            he_resLoc1_Carbon->Fill(res_eLOC1_Carbon->at(im));
            he_resTheta_Carbon->Fill(res_ePHI_Carbon->at(im));
            he_resPhi_Carbon->Fill(res_eTHETA_Carbon->at(im));
            he_resQOP_Carbon->Fill(res_eQOP_Carbon->at(im));
        }
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 1600, 1200);
    c1->Divide(3, 2);
    c1->cd(1);
    he_resLoc0_Base->Draw();
    he_resLoc0_Carbon->Draw("same");
    c1->cd(2);
    he_resLoc1_Base->Draw();
    he_resLoc1_Carbon->Draw("same");
    c1->cd(3);
    he_resPhi_Base->Draw();
    he_resPhi_Carbon->Draw("same");
    c1->cd(4);
    he_resTheta_Base->Draw();
    he_resTheta_Carbon->Draw("same");
    c1->cd(5);
    he_resQOP_Base->Draw();
    he_resQOP_Carbon->Draw("same");
    c1->SaveAs("checkResolution_lay12.png");
}