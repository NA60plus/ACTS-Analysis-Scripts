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

TH1D *ProjectMean(const TH2D *h2, TString suffix, TString xtitle = "", TString ytitle = "", TString title = "")
{
    if (!h2)
        return nullptr;

    int nBinsX = h2->GetNbinsX();
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();

    TH1D *hMean = new TH1D(Form("hMean%s", suffix.Data()), Form("%s;%s;%s", title.Data(), xtitle.Data(), ytitle.Data()), nBinsX, xMin, xMax);

    for (int i = 1; i <= nBinsX; ++i)
    {
        TH1D *proj = h2->ProjectionY("_py", i, i);
        if (proj->GetEntries() > 0)
        {
            double mean = proj->GetMean();
            double meanErr = proj->GetMeanError();
            hMean->SetBinContent(i, mean);
            hMean->SetBinError(i, meanErr);
        }
        else
        {
            hMean->SetBinContent(i, 0);
        }
        delete proj; // Clean up the temporary histogram
    }

    return hMean;
}

TH1D *ProjectRMS(const TH2D *h2, TString suffix, TString xtitle = "", TString ytitle = "", TString title = "")
{
    if (!h2)
        return nullptr;

    int nBinsX = h2->GetNbinsX();
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();

    TH1D *hRMS = new TH1D(Form("hRMS%s", suffix.Data()), Form("%s;%s;%s", title.Data(), xtitle.Data(), ytitle.Data()), nBinsX, xMin, xMax);

    for (int i = 1; i <= nBinsX; ++i)
    {
        TH1D *proj = h2->ProjectionY("_py", i, i);
        if (proj->GetEntries() > 0)
        {
            double rms = proj->GetRMS();
            double rmsErr = proj->GetRMSError();
            hRMS->SetBinContent(i, rms);
            hRMS->SetBinError(i, rmsErr);
        }
        else
        {
            hRMS->SetBinContent(i, 0);
        }
        delete proj; // Clean up the temporary histogram
    }

    return hRMS;
}

void CheckRefit(const char *filenameRef = "tracksummary_ambims.root", const char *filenameVT = "tracksummary_ambims.root", const char *filenameMS = "tracksummary_ambims.root", const char *treename = "tracksummary", const char *outputname = "jpsi", int pdgCode = 223)
{

    TFile *fileRef = TFile::Open(filenameRef);
    TTree *treeRef = (TTree *)fileRef->Get(treename);
    TFile *fileVT = TFile::Open(filenameVT);
    TTree *treeVT = (TTree *)fileVT->Get(treename);
    TFile *fileMS = TFile::Open(filenameMS);
    TTree *treeMS = (TTree *)fileMS->Get(treename);

    float ptMin = 0;
    float ptMax = 3;
    float yMin = 0.;
    float yMax = 5.5;
    float pMin = 0;
    float pMax = 30;

    //---------------------------------------
    // create histos
    //---------------------------------------
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    // Retrieve the particle from the database
    TParticlePDG *particle = pdgDB->GetParticle(pdgCode);
    double Mdileton = particle->Mass();
    double deltaM = 0.7;

    std::string rootname = std::string("invmass") + outputname + std::string(".root");
    TFile *fout = new TFile(rootname.c_str(), "recreate");
    TH1D *he_mass_Ref_min = new TH1D("he_mass_Ref_min", ";M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *he_mass_Ref_max = new TH1D("he_mass_Ref_max", ";M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *he_mass_best = new TH1D("he_mass_best", ";M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *he_mass_Ref = new TH1D("he_mass_Ref", ";M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM);
    TH2D *he_mass_vs_pt_Ref = new TH2D("he_mass_vs_pt_Ref", ";#it{p}_{T} (GeV/#it{c};M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 40, 0, 3, 150, Mdileton-deltaM, Mdileton+deltaM);
    TH2D *he_mass_vs_y_Ref = new TH2D("he_mass_vs_y_Ref", ";#it{y};M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 40, 1.5, 5.5, 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *he_mass_VT = new TH1D("he_mass_VT", ";M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 150, Mdileton*0.8, Mdileton*0.8);
    TH2D *he_mass_vs_pt_VT = new TH2D("he_mass_vs_pt_VT", ";#it{p}_{T} (GeV/#it{c};M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 40, 0, 3, 150, Mdileton-deltaM, Mdileton+deltaM);
    TH2D *he_mass_vs_y_VT = new TH2D("he_mass_vs_y_VT", ";#it{y};M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 40, 1.5, 5.5, 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *he_mass_MS = new TH1D("he_mass_MS", ";M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 150, Mdileton*0.8, Mdileton*0.8);
    TH2D *he_mass_vs_pt_MS = new TH2D("he_mass_vs_pt_MS", ";#it{p}_{T} (GeV/#it{c};M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 40, 0, 3, 150, Mdileton-deltaM, Mdileton+deltaM);
    TH2D *he_mass_vs_y_MS = new TH2D("he_mass_vs_y_MS", ";#it{y};M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 40, 1.5, 5.5, 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *hMatchingEff = new TH1D("hMatchingEff",";#it{p}_{T} (GeV/#it{c}); Matching efficiency",10,0,3);
    TH1D *hMatchingRef = new TH1D("hMatchingRef",";#it{p}_{T} (GeV/#it{c}); Matching efficiency",10,0,3);
    TH2D *he_mass_Ref_vs_MS = new TH2D("he_mass_Ref_vs_MS", ";M_{#mu^{+}#mu^{-}}^{Refit} - M_{PDG} (GeV/#it{c}^{2});M_{#mu^{+}#mu^{-}}^{MS} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM, 150, Mdileton-deltaM, Mdileton+deltaM);
    TH2D *he_mass_Ref_vs_VT = new TH2D("he_mass_Ref_vs_VT", ";M_{#mu^{+}#mu^{-}}^{Refit} - M_{PDG} (GeV/#it{c}^{2});M_{#mu^{+}#mu^{-}}^{VT} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM, 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *he_mass_delta = new TH1D("he_mass_delta", ";|M_{#mu^{+}#mu^{-}}^{Ref} - M_{PDG}|-|M_{#mu^{+}#mu^{-}}^{VT} - M_{PDG}| (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *he_mass_delta_min = new TH1D("he_mass_delta_min", ";|M_{#mu^{+}#mu^{-}}^{Ref} - M_{PDG}|-|M_{#mu^{+}#mu^{-}}^{VT} - M_{PDG}| (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *he_mass_delta_maj = new TH1D("he_mass_delta_maj", ";|M_{#mu^{+}#mu^{-}}^{Ref} - M_{PDG}|-|M_{#mu^{+}#mu^{-}}^{VT} - M_{PDG}| (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM);
    TH2D *he_mass_delta_vs_MS = new TH2D("he_mass_delta_vs_MS", ";M_{#mu^{+}#mu^{-}} - M_{PDG} (GeV/#it{c}^{2});dN/dM", 150, Mdileton-deltaM, Mdileton+deltaM, 150, Mdileton-deltaM, Mdileton+deltaM);
    TH1D *hPRes_Ref = new TH1D("hPRes_Ref", ";(#it{p}^{Rec}-#it{p}^{Gen})/#it{p}^{Gen}; a.u.", 100, -0.3, 0.3);
    TH1D *hPRes_Ref_min = new TH1D("hPRes_Ref_min", ";(#it{p}^{Rec}-#it{p}^{Gen})/#it{p}^{Gen}; a.u.", 100, -0.3, 0.3);
    TH1D *hPRes_Ref_max = new TH1D("hPRes_Ref_max", ";(#it{p}^{Rec}-#it{p}^{Gen})/#it{p}^{Gen}; a.u.", 100, -0.3, 0.3);
    TH1D *hPRes_VT = new TH1D("hPRes_VT", ";(#it{p}^{Rec}-#it{p}^{Gen})/#it{p}^{Gen}; a.u.", 100, -0.3, 0.3);
    TH1D *hPRes_MS = new TH1D("hPRes_MS", ";(#it{p}^{Rec}-#it{p}^{Gen})/#it{p}^{Gen}; a.u.", 100, -0.3, 0.3);
    TH2D *hPRes_Ref_vs_MS = new TH2D("hPRes_Ref_vs_MS", ";(#it{p}^{Rec}_{MS}-#it{p}^{Gen})/#it{p}^{Gen}^{MS};(#it{p}^{Rec}_{VT+MS}-#it{p}^{Gen})/#it{p}^{Gen};", 100, -0.3, 0.3, 100, -0.3, 0.3);
    TH2D *hPRes_Ref_vs_VT = new TH2D("hPRes_Ref_vs_VT", ";(#it{p}^{Rec}_{VT}-#it{p}^{Gen})/#it{p}^{Gen}^{VT};(#it{p}^{Rec}_{VT+MS}-#it{p}^{Gen})/#it{p}^{Gen};", 100, -0.3, 0.3, 100, -0.3, 0.3);
    
    TH1D *hD0Res_Ref = new TH1D("hD0Res_Ref", ";(#it{d}_{0}^{Rec}-#it{d}_{0}^{Gen}); a.u.", 100, -0.5, 0.5);
    TH1D *hZ0Res_Ref = new TH1D("hZ0Res_Ref", ";(#it{z}_{0}^{Rec}-#it{z}_{0}^{Gen}); a.u.", 100, -5, 5);
    TH1D *hThetaRes_Ref = new TH1D("hThetaRes_Ref", ";(#theta^{Rec}-#theta^{Gen})/#theta^{Gen}; a.u.", 100, -0.02, 0.02);
    TH1D *hPhiRes_Ref = new TH1D("hPhiRes_Ref", ";(#phi^{Rec}-#phi^{Gen})/#phi^{Gen}; a.u.", 100, -0.02, 0.02);

    TH1D *hD0Res_Ref_min = new TH1D("hD0Res_Ref_min", ";(#it{d}_{0}^{Rec}-#it{d}_{0}^{Gen}); a.u.", 100, -0.5, 0.5);
    TH1D *hZ0Res_Ref_min = new TH1D("hZ0Res_Ref_min", ";(#it{z}_{0}^{Rec}-#it{z}_{0}^{Gen}); a.u.", 100, -5, 5);
    TH1D *hThetaRes_Ref_min = new TH1D("hThetaRes_Ref_min", ";(#theta^{Rec}-#theta^{Gen})/#theta^{Gen}; a.u.", 100, -0.02, 0.02);
    TH1D *hPhiRes_Ref_min = new TH1D("hPhiRes_Ref_min", ";(#phi^{Rec}-#phi^{Gen})/#phi^{Gen}; a.u.", 100, -0.02, 0.02);

    TH1D *hD0Res_Ref_max = new TH1D("hD0Res_Ref_max", ";(#it{d}_{0}^{Rec}-#it{d}_{0}^{Gen}); a.u.", 100, -0.5, 0.5);
    TH1D *hZ0Res_Ref_max = new TH1D("hZ0Res_Ref_max", ";(#it{z}_{0}^{Rec}-#it{z}_{0}^{Gen}); a.u.", 100, -5, 5);
    TH1D *hThetaRes_Ref_max = new TH1D("hThetaRes_Ref_max", ";(#theta^{Rec}-#theta^{Gen})/#theta^{Gen}; a.u.", 100, -0.02, 0.02);
    TH1D *hPhiRes_Ref_max = new TH1D("hPhiRes_Ref_max", ";(#phi^{Rec}-#phi^{Gen})/#phi^{Gen}; a.u.", 100, -0.02, 0.02);

    TH1D *hD0Res_VT = new TH1D("hD0Res_VT", ";(#it{d}_{0}^{Rec}-#it{d}_{0}^{Gen}); a.u.", 100, -0.5, 0.5);
    TH1D *hZ0Res_VT = new TH1D("hZ0Res_VT", ";(#it{z}_{0}^{Rec}-#it{z}_{0}^{Gen}); a.u.", 100, -5, 5);
    TH1D *hThetaRes_VT = new TH1D("hThetaRes_VT", ";(#theta^{Rec}-#theta^{Gen})/#theta^{Gen}; a.u.", 100, -0.02, 0.02);
    TH1D *hPhiRes_VT = new TH1D("hPhiRes_VT", ";(#phi^{Rec}-#phi^{Gen})/#phi^{Gen}; a.u.", 100, -0.02, 0.02);

    //  Define variables to hold kinematic data of true particles
    std::vector<float> *t_px_Ref = new std::vector<float>;
    std::vector<float> *t_py_Ref = nullptr;
    std::vector<float> *t_pz_Ref = nullptr;
    std::vector<float> *t_p_Ref = nullptr;
    std::vector<float> *t_z0_Ref = nullptr;
    std::vector<float> *t_d0_Ref = nullptr;
    std::vector<float> *t_theta_Ref = nullptr;
    std::vector<float> *t_phi_Ref = nullptr;
    std::vector<int> *t_nmu_Ref = new std::vector<int>;
    std::vector<int> *t_majPartId_Ref = new std::vector<int>;

    std::vector<float> *eLOC0_fit_Ref = new std::vector<float>;
    std::vector<float> *eLOC1_fit_Ref = new std::vector<float>;
    std::vector<float> *ePHI_fit_Ref = new std::vector<float>;
    std::vector<float> *eTHETA_fit_Ref = new std::vector<float>;
    std::vector<float> *eQOP_fit_Ref = new std::vector<float>;
    std::vector<float> *err_eQOP_fit_Ref = new std::vector<float>;
    //  Define variables to hold kinematic data of true particles
    std::vector<float> *t_px_VT = new std::vector<float>;
    std::vector<float> *t_py_VT = nullptr;
    std::vector<float> *t_pz_VT = nullptr;
    std::vector<float> *t_p_VT = nullptr;
    std::vector<float> *t_z0_VT = nullptr;
    std::vector<float> *t_d0_VT = nullptr;
    std::vector<float> *t_theta_VT = nullptr;
    std::vector<float> *t_phi_VT = nullptr;
    std::vector<int> *t_nmu_VT = new std::vector<int>;
    std::vector<int> *t_majPartId_VT = new std::vector<int>;

    std::vector<float> *eLOC0_fit_VT = new std::vector<float>;
    std::vector<float> *eLOC1_fit_VT = new std::vector<float>;
    std::vector<float> *ePHI_fit_VT = new std::vector<float>;
    std::vector<float> *eTHETA_fit_VT = new std::vector<float>;
    std::vector<float> *eQOP_fit_VT = new std::vector<float>;
    std::vector<float> *err_eQOP_fit_VT = new std::vector<float>;
    //  Define variables to hold kinematic data of true particles
    std::vector<float> *t_px_MS = new std::vector<float>;
    std::vector<float> *t_py_MS = nullptr;
    std::vector<float> *t_pz_MS = nullptr;
    std::vector<float> *t_p_MS = nullptr;
    std::vector<float> *t_z0_MS = nullptr;
    std::vector<float> *t_d0_MS = nullptr;
    std::vector<float> *t_theta_MS = nullptr;
    std::vector<float> *t_phi_MS = nullptr;
    std::vector<int> *t_nmu_MS = new std::vector<int>;
    std::vector<int> *t_majPartId_MS = new std::vector<int>;

    std::vector<float> *eLOC0_fit_MS = new std::vector<float>;
    std::vector<float> *eLOC1_fit_MS = new std::vector<float>;
    std::vector<float> *ePHI_fit_MS = new std::vector<float>;
    std::vector<float> *eTHETA_fit_MS = new std::vector<float>;
    std::vector<float> *eQOP_fit_MS = new std::vector<float>;
    std::vector<float> *err_eQOP_fit_MS = new std::vector<float>;

    // Set branch addresses for truth particle
    treeRef->SetBranchAddress("t_px", &t_px_Ref);
    treeRef->SetBranchAddress("t_py", &t_py_Ref);
    treeRef->SetBranchAddress("t_pz", &t_pz_Ref);
    treeRef->SetBranchAddress("t_p", &t_p_Ref);
    treeRef->SetBranchAddress("t_phi", &t_phi_Ref);
    treeRef->SetBranchAddress("t_theta", &t_theta_Ref);
    treeRef->SetBranchAddress("t_d0", &t_d0_Ref);
    treeRef->SetBranchAddress("t_z0", &t_z0_Ref);
    treeRef->SetBranchAddress("track_nr", &t_nmu_Ref);
    treeRef->SetBranchAddress("track_nr", &t_nmu_Ref);
    treeRef->SetBranchAddress("majorityParticleId", &t_majPartId_Ref);
    treeRef->SetBranchAddress("eLOC0_fit", &eLOC0_fit_Ref);
    treeRef->SetBranchAddress("eLOC1_fit", &eLOC1_fit_Ref);
    treeRef->SetBranchAddress("ePHI_fit", &ePHI_fit_Ref);
    treeRef->SetBranchAddress("eTHETA_fit", &eTHETA_fit_Ref);
    treeRef->SetBranchAddress("eQOP_fit", &eQOP_fit_Ref);
    treeRef->SetBranchAddress("err_eQOP_fit", &err_eQOP_fit_Ref);
    // Set branch addresses for truth particle
    treeVT->SetBranchAddress("t_px", &t_px_VT);
    treeVT->SetBranchAddress("t_py", &t_py_VT);
    treeVT->SetBranchAddress("t_pz", &t_pz_VT);
    treeVT->SetBranchAddress("t_p", &t_p_VT);
    treeVT->SetBranchAddress("t_phi", &t_phi_VT);
    treeVT->SetBranchAddress("t_theta", &t_theta_VT);
    treeVT->SetBranchAddress("t_d0", &t_d0_VT);
    treeVT->SetBranchAddress("t_z0", &t_z0_VT);
    treeVT->SetBranchAddress("track_nr", &t_nmu_VT);
    treeVT->SetBranchAddress("track_nr", &t_nmu_VT);
    treeVT->SetBranchAddress("majorityParticleId", &t_majPartId_VT);
    treeVT->SetBranchAddress("eLOC0_fit", &eLOC0_fit_VT);
    treeVT->SetBranchAddress("eLOC1_fit", &eLOC1_fit_VT);
    treeVT->SetBranchAddress("ePHI_fit", &ePHI_fit_VT);
    treeVT->SetBranchAddress("eTHETA_fit", &eTHETA_fit_VT);
    treeVT->SetBranchAddress("eQOP_fit", &eQOP_fit_VT);
    treeVT->SetBranchAddress("err_eQOP_fit", &err_eQOP_fit_VT);
    // Set branch addresses for truth particle
    treeMS->SetBranchAddress("t_px", &t_px_MS);
    treeMS->SetBranchAddress("t_py", &t_py_MS);
    treeMS->SetBranchAddress("t_pz", &t_pz_MS);
    treeMS->SetBranchAddress("t_p", &t_p_MS);
    treeMS->SetBranchAddress("t_phi", &t_phi_MS);
    treeMS->SetBranchAddress("t_theta", &t_theta_MS);
    treeMS->SetBranchAddress("t_d0", &t_d0_MS);
    treeMS->SetBranchAddress("t_z0", &t_z0_MS);
    treeMS->SetBranchAddress("track_nr", &t_nmu_MS);
    treeMS->SetBranchAddress("majorityParticleId", &t_majPartId_MS);
    treeMS->SetBranchAddress("eLOC0_fit", &eLOC0_fit_MS);
    treeMS->SetBranchAddress("eLOC1_fit", &eLOC1_fit_MS);
    treeMS->SetBranchAddress("ePHI_fit", &ePHI_fit_MS);
    treeMS->SetBranchAddress("eTHETA_fit", &eTHETA_fit_MS);
    treeMS->SetBranchAddress("eQOP_fit", &eQOP_fit_MS);
    treeMS->SetBranchAddress("err_eQOP_fit", &err_eQOP_fit_MS);

    // Create TLorentzVector objects for the muons
    TLorentzVector t_muon1, t_muon2;
    TLorentzVector e_muon1_fit, e_muon2_fit;

    Double_t mass = 0.1056583755;

    //---------------------------------------
    // loop over entries
    //---------------------------------------

    Long64_t nentries = treeRef->GetEntries();
    Double_t onemuon = 0;
    Double_t zeromuon = 0;
    Double_t twomuons = 0;

    for (Long64_t i = 0; i < nentries; ++i)
    {
        treeRef->GetEntry(i);
        treeVT->GetEntry(i);
        treeMS->GetEntry(i);
        if(t_nmu_MS->size() != t_nmu_Ref->size() || t_nmu_VT->size() != t_nmu_Ref->size()){
            continue;
        }
        for (int im = 0; im < t_nmu_VT->size(); im++){
            double pVT = 1.0 / std::abs(eQOP_fit_VT->at(im));
            double pRef = 1.0 / std::abs(eQOP_fit_Ref->at(im));
            double pMS = 1.0 / std::abs(eQOP_fit_MS->at(im));
            double px = pVT * std::sin(eTHETA_fit_VT->at(im)) * std::cos(ePHI_fit_VT->at(im));
            double py = pVT * std::sin(eTHETA_fit_VT->at(im)) * std::sin(ePHI_fit_VT->at(im));
            double pz = pVT * std::cos(eTHETA_fit_VT->at(im));
            double E = std::sqrt(px * px + py * py + pz * pz + mass * mass);
            hPRes_Ref->Fill((pRef-t_p_Ref->at(im))/t_p_Ref->at(im));
            if((pMS-t_p_Ref->at(im))/t_p_Ref->at(im) > -0.03){
                hPRes_Ref_max->Fill((pRef-t_p_Ref->at(im))/t_p_Ref->at(im));
            
                hD0Res_Ref_max->Fill(eLOC0_fit_Ref->at(im));
                hZ0Res_Ref_max->Fill(eLOC1_fit_Ref->at(im));
                hThetaRes_Ref_max->Fill((eTHETA_fit_Ref->at(im)-t_theta_Ref->at(im))/t_theta_Ref->at(im));
                hPhiRes_Ref_max->Fill((ePHI_fit_Ref->at(im)-t_phi_Ref->at(im))/t_phi_Ref->at(im));
            }
            else{
                hPRes_Ref_min->Fill((pRef-t_p_Ref->at(im))/t_p_Ref->at(im));
            
                hD0Res_Ref_min->Fill(eLOC0_fit_Ref->at(im));
                hZ0Res_Ref_min->Fill(eLOC1_fit_Ref->at(im));
                hThetaRes_Ref_min->Fill((eTHETA_fit_Ref->at(im)-t_theta_Ref->at(im))/t_theta_Ref->at(im));
                hPhiRes_Ref_min->Fill((ePHI_fit_Ref->at(im)-t_phi_Ref->at(im))/t_phi_Ref->at(im));
            }

            hPRes_VT->Fill((pVT-t_p_VT->at(im))/t_p_VT->at(im));
            hPRes_MS->Fill((pMS-t_p_MS->at(im))/t_p_MS->at(im));

            hPRes_Ref_vs_MS->Fill((pMS-t_p_MS->at(im))/t_p_MS->at(im),(pRef-t_p_Ref->at(im))/t_p_Ref->at(im));
            hPRes_Ref_vs_VT->Fill((pVT-t_p_VT->at(im))/t_p_VT->at(im),(pRef-t_p_Ref->at(im))/t_p_Ref->at(im));
        
            e_muon1_fit.SetPxPyPzE(px, py, pz, E);

            hMatchingRef->Fill(TMath::Sqrt(px*px+py*py));
            if (t_majPartId_MS->at(im) == t_majPartId_VT->at(im)){
                hMatchingEff->Fill(TMath::Sqrt(px*px+py*py));
            }

            hD0Res_Ref->Fill(eLOC0_fit_Ref->at(im));
            hZ0Res_Ref->Fill(eLOC1_fit_Ref->at(im));
            hThetaRes_Ref->Fill((eTHETA_fit_Ref->at(im)-t_theta_Ref->at(im))/t_theta_Ref->at(im));
            hPhiRes_Ref->Fill((ePHI_fit_Ref->at(im)-t_phi_Ref->at(im))/t_phi_Ref->at(im));
            hD0Res_VT->Fill((eLOC0_fit_VT->at(im)-t_d0_VT->at(im)));
            hZ0Res_VT->Fill((eLOC1_fit_VT->at(im)-t_z0_VT->at(im)));
            hThetaRes_VT->Fill((eTHETA_fit_VT->at(im)-t_theta_VT->at(im))/t_theta_VT->at(im));
            hPhiRes_VT->Fill((ePHI_fit_VT->at(im)-t_phi_VT->at(im))/t_phi_VT->at(im));
            
        }
        if (t_nmu_Ref->size() < 2)
            continue;

        for (int im = 0; im < t_nmu_Ref->size() - 1; im++)
        {
            double pgen1 = t_p_Ref->at(im);
            double p1_Ref = 1.0 / std::abs(eQOP_fit_Ref->at(im));
            double px1_Ref = p1_Ref * std::sin(eTHETA_fit_Ref->at(im)) * std::cos(ePHI_fit_Ref->at(im));
            double py1_Ref = p1_Ref * std::sin(eTHETA_fit_Ref->at(im)) * std::sin(ePHI_fit_Ref->at(im));
            double pz1_Ref = p1_Ref * std::cos(eTHETA_fit_Ref->at(im));
            double E1_Ref = std::sqrt(px1_Ref * px1_Ref + py1_Ref * py1_Ref + pz1_Ref * pz1_Ref + mass * mass);

            double p1_VT = 1.0 / std::abs(eQOP_fit_VT->at(im));
            double px1_VT = p1_VT * std::sin(eTHETA_fit_VT->at(im)) * std::cos(ePHI_fit_VT->at(im));
            double py1_VT = p1_VT * std::sin(eTHETA_fit_VT->at(im)) * std::sin(ePHI_fit_VT->at(im));
            double pz1_VT = p1_VT * std::cos(eTHETA_fit_VT->at(im));
            double E1_VT = std::sqrt(px1_VT * px1_VT + py1_VT * py1_VT + pz1_VT * pz1_VT + mass * mass);

            double p1_MS = 1.0 / std::abs(eQOP_fit_MS->at(im));
            double px1_MS = p1_MS * std::sin(eTHETA_fit_MS->at(im)) * std::cos(ePHI_fit_MS->at(im));
            double py1_MS = p1_MS * std::sin(eTHETA_fit_MS->at(im)) * std::sin(ePHI_fit_MS->at(im));
            double pz1_MS = p1_MS * std::cos(eTHETA_fit_MS->at(im));
            double E1_MS = std::sqrt(px1_MS * px1_MS + py1_MS * py1_MS + pz1_MS * pz1_MS + mass * mass);

            for (int im2 = im + 1; im2 < t_nmu_Ref->size(); im2++)
            {
                //if(t_majPartId_VT->at(im)!=t_majPartId_MS->at(im)) continue;
                //if(t_majPartId_VT->at(im2)!=t_majPartId_MS->at(im2)) continue;
                double pgen2 = t_p_Ref->at(im2);
                double p2_Ref = 1.0 / std::abs(eQOP_fit_Ref->at(im2));
                double px2_Ref = p2_Ref * std::sin(eTHETA_fit_Ref->at(im2)) * std::cos(ePHI_fit_Ref->at(im2));
                double py2_Ref = p2_Ref * std::sin(eTHETA_fit_Ref->at(im2)) * std::sin(ePHI_fit_Ref->at(im2));
                double pz2_Ref = p2_Ref * std::cos(eTHETA_fit_Ref->at(im2));
                double E2_Ref = std::sqrt(px2_Ref * px2_Ref + py2_Ref * py2_Ref + pz2_Ref * pz2_Ref + mass * mass);

                double p2_VT = 1.0 / std::abs(eQOP_fit_VT->at(im2));
                double px2_VT = p2_VT * std::sin(eTHETA_fit_VT->at(im2)) * std::cos(ePHI_fit_VT->at(im2));
                double py2_VT = p2_VT * std::sin(eTHETA_fit_VT->at(im2)) * std::sin(ePHI_fit_VT->at(im2));
                double pz2_VT = p2_VT * std::cos(eTHETA_fit_VT->at(im2));
                double E2_VT = std::sqrt(px2_VT * px2_VT + py2_VT * py2_VT + pz2_VT * pz2_VT + mass * mass);

                double p2_MS = 1.0 / std::abs(eQOP_fit_MS->at(im2));
                double px2_MS = p2_MS * std::sin(eTHETA_fit_MS->at(im2)) * std::cos(ePHI_fit_MS->at(im2));
                double py2_MS = p2_MS * std::sin(eTHETA_fit_MS->at(im2)) * std::sin(ePHI_fit_MS->at(im2));
                double pz2_MS = p2_MS * std::cos(eTHETA_fit_MS->at(im2));
                double E2_MS = std::sqrt(px2_MS * px2_MS + py2_MS * py2_MS + pz2_MS * pz2_MS + mass * mass);

                e_muon1_fit.SetPxPyPzE(px1_Ref, py1_Ref, pz1_Ref, E1_Ref);
                e_muon2_fit.SetPxPyPzE(px2_Ref, py2_Ref, pz2_Ref, E2_Ref);

                double e_mass_Ref = (e_muon1_fit + e_muon2_fit).M();
                double e_pt_Ref = (e_muon1_fit + e_muon2_fit).Pt();
                double e_y_Ref = (e_muon1_fit + e_muon2_fit).Rapidity();

                e_muon1_fit.SetPxPyPzE(px1_VT, py1_VT, pz1_VT, E1_VT);
                e_muon2_fit.SetPxPyPzE(px2_VT, py2_VT, pz2_VT, E2_VT);

                double e_mass_VT = (e_muon1_fit + e_muon2_fit).M();
                double e_pt_VT = (e_muon1_fit + e_muon2_fit).Pt();
                double e_y_VT = (e_muon1_fit + e_muon2_fit).Rapidity();

                e_muon1_fit.SetPxPyPzE(px1_MS, py1_MS, pz1_MS, E1_MS);
                e_muon2_fit.SetPxPyPzE(px2_MS, py2_MS, pz2_MS, E2_MS);

                double e_mass_MS = (e_muon1_fit + e_muon2_fit).M();
                double e_pt_MS = (e_muon1_fit + e_muon2_fit).Pt();
                double e_y_MS = (e_muon1_fit + e_muon2_fit).Rapidity();

                double pMS = 1.0 / std::abs(eQOP_fit_MS->at(im));
                double pMS2 = 1.0 / std::abs(eQOP_fit_MS->at(im2));
                double pVT = 1.0 / std::abs(eQOP_fit_VT->at(im));
                double pVT2 = 1.0 / std::abs(eQOP_fit_VT->at(im2));
                //if(TMath::Abs(e_mass_MS) < 0.02)

                if((pMS-t_p_Ref->at(im))/t_p_Ref->at(im) > -0.03 && (pMS2-t_p_Ref->at(im2))/t_p_Ref->at(im2) > -0.03){
                    he_mass_Ref_max->Fill(e_mass_Ref);
                    
                }
                else{
                    he_mass_Ref_min->Fill(e_mass_Ref);
                }

                if((pMS-pVT)/pVT > -0.03 && (pMS2-pVT2)/pVT2 > -0.03){
                    he_mass_Ref_max->Fill(e_mass_Ref);
                    he_mass_best->Fill(e_mass_Ref);
                    
                }
                else{
                    he_mass_Ref_min->Fill(e_mass_Ref);
                    he_mass_best->Fill(e_mass_VT);
                }
                he_mass_Ref->Fill(e_mass_Ref);
                he_mass_vs_y_Ref->Fill(e_y_Ref, e_mass_Ref);
                he_mass_vs_pt_Ref->Fill(e_pt_Ref, e_mass_Ref);
                he_mass_VT->Fill(e_mass_VT);
                he_mass_vs_y_VT->Fill(e_y_VT, e_mass_VT);
                he_mass_vs_pt_VT->Fill(e_pt_VT, e_mass_VT);
                he_mass_MS->Fill(e_mass_MS);
                he_mass_vs_y_MS->Fill(e_y_MS, e_mass_MS);
                he_mass_vs_pt_MS->Fill(e_pt_MS, e_mass_MS);

                he_mass_Ref_vs_MS->Fill( e_mass_MS,e_mass_Ref);
                he_mass_Ref_vs_VT->Fill(e_mass_VT,e_mass_Ref);
                he_mass_delta->Fill(TMath::Abs(e_mass_Ref)-TMath::Abs(e_mass_VT));
                he_mass_delta_vs_MS->Fill(TMath::Abs(e_mass_Ref)-TMath::Abs(e_mass_VT),e_mass_MS);

                if((pMS-t_p_Ref->at(im))/t_p_Ref->at(im) > -0.03 && (pMS2-t_p_Ref->at(im2))/t_p_Ref->at(im2) > -0.03)
                    he_mass_delta_maj->Fill(TMath::Abs(e_mass_Ref)-TMath::Abs(e_mass_VT));
                else
                    he_mass_delta_min->Fill(TMath::Abs(e_mass_Ref)-TMath::Abs(e_mass_VT));
            }
        }
    }
    hMatchingEff->Divide(hMatchingRef);

    std::string matchname = outputname + std::string("_matching.root");
    TFile* f = new TFile(matchname.c_str(),"recreate");


    he_mass_Ref_max->Scale(1./he_mass_Ref_max->GetEntries());
    he_mass_Ref_min->Scale(1./he_mass_Ref_min->GetEntries());
    he_mass_best->Scale(1./he_mass_best->GetEntries());
    he_mass_Ref->Scale(1./he_mass_Ref->GetEntries());
    he_mass_VT->Scale(1./he_mass_VT->GetEntries());
    he_mass_MS->Scale(1./he_mass_MS->GetEntries());
    hD0Res_Ref->Write();
    hZ0Res_Ref->Write();
    hThetaRes_Ref->Write();
    hPhiRes_Ref->Write();
    
    hD0Res_Ref_max->Write();
    hZ0Res_Ref_max->Write();
    hThetaRes_Ref_max->Write();
    hPhiRes_Ref_max->Write();
    
    hD0Res_Ref_min->Write();
    hZ0Res_Ref_min->Write();
    hThetaRes_Ref_min->Write();
    hPhiRes_Ref_min->Write();
    
    hD0Res_VT->Write();
    hZ0Res_VT->Write();
    hThetaRes_VT->Write();
    hPhiRes_VT->Write();

    hMatchingEff->Write();
    he_mass_Ref_vs_MS->Write();
    he_mass_Ref_vs_VT->Write();
    he_mass_delta->Write();
    he_mass_delta_vs_MS->Write();
    he_mass_delta_maj->Write();
    he_mass_delta_min->Write();

    he_mass_Ref_max->Write();
    he_mass_Ref_min->Write();
    he_mass_best->Write();
    he_mass_Ref->Write();
    he_mass_vs_y_Ref->Write();
    he_mass_vs_pt_Ref->Write(); 
    he_mass_VT->Write();
    he_mass_vs_y_VT->Write();
    he_mass_vs_pt_VT->Write();
    he_mass_MS->Write();
    he_mass_vs_y_MS->Write();
    he_mass_vs_pt_MS->Write();

    
    hD0Res_Ref_min->Scale(1./hD0Res_Ref_min->GetEntries());
    hD0Res_Ref_max->Scale(1./hD0Res_Ref_max->GetEntries());
    hD0Res_Ref->Scale(1./hD0Res_Ref->GetEntries());
    hD0Res_VT->Scale(1./hD0Res_VT->GetEntries());

    hZ0Res_Ref_min->Scale(1./hZ0Res_Ref_min->GetEntries());
    hZ0Res_Ref_max->Scale(1./hZ0Res_Ref_max->GetEntries());
    hZ0Res_Ref->Scale(1./hZ0Res_Ref->GetEntries());
    hZ0Res_VT->Scale(1./hZ0Res_VT->GetEntries());

    hPhiRes_Ref_min->Scale(1./hPhiRes_Ref_min->GetEntries());
    hPhiRes_Ref_max->Scale(1./hPhiRes_Ref_max->GetEntries());
    hPhiRes_Ref->Scale(1./hPhiRes_Ref->GetEntries());
    hPhiRes_VT->Scale(1./hPhiRes_VT->GetEntries());

    hThetaRes_Ref_min->Scale(1./hThetaRes_Ref_min->GetEntries());
    hThetaRes_Ref_max->Scale(1./hThetaRes_Ref_max->GetEntries());
    hThetaRes_Ref->Scale(1./hThetaRes_Ref->GetEntries());
    hThetaRes_VT->Scale(1./hThetaRes_VT->GetEntries());

    hPRes_Ref_min->Scale(1./hPRes_Ref_min->GetEntries());
    hPRes_Ref_max->Scale(1./hPRes_Ref_max->GetEntries());
    hPRes_Ref->Scale(1./hPRes_Ref->GetEntries());
    hPRes_VT->Scale(1./hPRes_VT->GetEntries());
    hPRes_MS->Scale(1./hPRes_MS->GetEntries());
    hPRes_Ref->Write();
    hPRes_Ref_max->Write();
    hPRes_Ref_min->Write();
    hPRes_VT->Write();
    hPRes_MS->Write();
    hPRes_Ref_vs_MS->Write();
    hPRes_Ref_vs_VT->Write();
    TCanvas* cv = new TCanvas("cv", "", 1800,1600);

    he_mass_Ref->SetLineColor(kRed);
    he_mass_VT->SetLineColor(kBlue);
    he_mass_MS->SetLineColor(kGreen+2);
    he_mass_Ref->Draw();
    he_mass_VT->Draw("same");
    he_mass_MS->Draw("same");
    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
    /*
    he_mass_Ref->Fit("gaus", "", "", -deltaM, deltaM);
    TF1 *fitFunction_Ref = he_mass_Ref->GetFunction("gaus");
    he_mass_VT->Fit("gaus", "", "", -deltaM/2., deltaM/2.);
    TF1 *fitFunction_VT = he_mass_VT->GetFunction("gaus");
    //TF1 *fitFunction_VT = new TF1("fitFunction_VT", "gaus", -deltaM, deltaM);
    //fitFunction_VT->SetParameters(1000,0,7);
    he_mass_MS->Fit("gaus", "", "", -deltaM, deltaM);
    TF1 *fitFunction_MS = he_mass_MS->GetFunction("gaus");
    
    float sigma_Ref = 1000*fitFunction_Ref->GetParameter(2);
    float err_Ref = 1000*fitFunction_Ref->GetParError(2);
    float sigma_VT = 1000*fitFunction_VT->GetParameter(2);
    float err_VT = 1000*fitFunction_VT->GetParError(2);
    float sigma_MS = 1000*fitFunction_MS->GetParameter(2);
    float err_MS = 1000*fitFunction_MS->GetParError(2);

    fitFunction_Ref->SetLineColor(kRed);
    fitFunction_VT->SetLineColor(kBlue);
    fitFunction_MS->SetLineColor(kGreen+2);

    leg->AddEntry(he_mass_Ref,Form("VT+MS RMS = %0.f #pm %0.f MeV/#it{c}^{2}, #sigma = %0.f #pm %0.f MeV/#it{c}^{2}",he_mass_Ref->GetRMS()*1000.,he_mass_Ref->GetRMSError()*1000., sigma_Ref, err_Ref), "l");
    leg->AddEntry(he_mass_VT,Form("VT RMS = %0.f #pm %0.f MeV/#it{c}^{2}, #sigma = %0.f #pm %0.f MeV/#it{c}^{2}",he_mass_VT->GetRMS()*1000.,he_mass_VT->GetRMSError()*1000., sigma_VT, err_VT), "l");
    leg->AddEntry(he_mass_MS,Form("MS RMS = %0.f #pm %0.f MeV/#it{c}^{2}, #sigma = %0.f #pm %0.f MeV/#it{c}^{2}",he_mass_MS->GetRMS()*1000.,he_mass_MS->GetRMSError()*1000., sigma_MS, err_MS), "l");
    leg->Draw();
    */
    std::string pngname = outputname + std::string("_mass_comp.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());
    gStyle->SetOptStat(111111);

    TLegend* leg2 = new TLegend(0.6,0.6,0.9,0.9);
    // p res
    hPRes_Ref_min->SetLineColor(kRed);
    hPRes_Ref_max->SetLineColor(kBlue);
    hPRes_VT->SetLineColor(kGreen+2);
    hPRes_Ref_max->Draw();
    hPRes_Ref_min->Draw("same");
    hPRes_VT->Draw("same");
    leg2->AddEntry(hPRes_Ref_max,"VT+MS with small E loss", "l");
    leg2->AddEntry(hPRes_Ref_min,"VT+MS with large E loss", "l");
    leg2->AddEntry(hPRes_VT,"VT only", "l");
    leg2->Draw();

    pngname = outputname + std::string("_ptres_comp.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());
    // d0 res
    hD0Res_Ref_min->SetLineColor(kRed);
    hD0Res_Ref_max->SetLineColor(kBlue);
    hD0Res_VT->SetLineColor(kGreen+2);
    hD0Res_Ref_max->Draw();
    hD0Res_Ref_min->Draw("same");
    hD0Res_VT->Draw("same");
    leg2->Draw();

    pngname = outputname + std::string("_d0res_comp.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());
    // z0 res
    hZ0Res_Ref_min->SetLineColor(kRed);
    hZ0Res_Ref_max->SetLineColor(kBlue);
    hZ0Res_VT->SetLineColor(kGreen+2);
    hZ0Res_Ref_max->Draw();
    hZ0Res_Ref_min->Draw("same");
    hZ0Res_VT->Draw("same");
    leg2->Draw();

    pngname = outputname + std::string("_z0res_comp.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());
    // phi res
    hPhiRes_Ref_min->SetLineColor(kRed);
    hPhiRes_Ref_max->SetLineColor(kBlue);
    hPhiRes_VT->SetLineColor(kGreen+2);
    hPhiRes_Ref_max->Draw();
    hPhiRes_Ref_min->Draw("same");
    hPhiRes_VT->Draw("same");
    leg2->Draw();

    pngname = outputname + std::string("_phires_comp.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());
    // theta res
    hThetaRes_Ref_min->SetLineColor(kRed);
    hThetaRes_Ref_max->SetLineColor(kBlue);
    hThetaRes_VT->SetLineColor(kGreen+2);
    hThetaRes_Ref_max->Draw();
    hThetaRes_Ref_min->Draw("same");
    hThetaRes_VT->Draw("same");
    leg2->Draw();

    pngname = outputname + std::string("_thetares_comp.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());
    
    he_mass_Ref_min->SetLineColor(kRed);
    he_mass_Ref_max->SetLineColor(kBlue);
    he_mass_VT->SetLineColor(kGreen+2);
    he_mass_Ref_max->Draw("");
    he_mass_Ref_min->Draw("same ");
    he_mass_VT->Draw("same ");
    leg2->Draw();

    pngname = outputname + std::string("_tmassres_comp.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());


    hPRes_Ref_vs_MS->Draw("colz");
    pngname = outputname + std::string("_ref_vs_ms.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());
    /////////////////////////////////

    TLegend* leg3 = new TLegend(0.6,0.6,0.9,0.9);

    hPRes_Ref->SetLineColor(kRed);
    hPRes_VT->SetLineColor(kBlue);
    hPRes_MS->SetLineColor(kGreen+2);
    hPRes_VT->Draw();
    hPRes_Ref->Draw("same");
    hPRes_MS->Draw("same");
    leg3->AddEntry(hPRes_Ref,"VT+MS", "l");
    leg3->AddEntry(hPRes_VT,"VT only", "l");
    leg3->AddEntry(hPRes_MS,"MS only", "l");
    leg3->Draw();

    pngname = outputname + std::string("_pref_comp.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());
    ///////////////////////////////////////////
    TLegend* leg4= new TLegend(0.6,0.6,0.9,0.9);
    // per il proposals
    he_mass_MS->SetLineColor(kRed);
    he_mass_VT->SetLineColor(kBlue);
    he_mass_VT->Draw();
    he_mass_MS->Draw("same");
    leg4->AddEntry(he_mass_VT,"Combined", "l");
    leg4->AddEntry(he_mass_MS,"MS only", "l");
    //leg4->Draw();

    pngname = outputname + std::string("_mbest.png");
    cv->Update();
    cv->SaveAs(pngname.c_str());

    std::cout<<"RMS RESULTS"<<std::endl;
    std::cout<<"Best m resolution: "<<he_mass_best->GetRMS()*1000.<<" +- "<<he_mass_best->GetRMSError()*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"VT m resolution: "<<he_mass_VT->GetRMS()*1000.<<" +- "<<he_mass_VT->GetRMSError()*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"MS m resolution: "<<he_mass_MS->GetRMS()*1000.<<" +- "<<he_mass_MS->GetRMSError()*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"VT+MS m resolution: "<<he_mass_Ref->GetRMS()*1000.<<" +- "<<he_mass_Ref->GetRMSError()*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"VT+MS m resolution with small E loss: "<<he_mass_Ref_max->GetRMS()*1000.<<" +- "<<he_mass_Ref_max->GetRMSError()*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"VT+MS m resolution with large E loss: "<<he_mass_Ref_min->GetRMS()*1000.<<" +- "<<he_mass_Ref_min->GetRMSError()*1000.<<" MeV/c^2"<<std::endl;
    
    he_mass_best->Fit("gaus", "", "", Mdileton-deltaM/2., Mdileton+ deltaM/2.);
    TF1 *fitFunction_best = he_mass_best->GetFunction("gaus");
    he_mass_VT->Fit("gaus", "", "", Mdileton-deltaM/2., Mdileton+ deltaM/2.);
    TF1 *fitFunction_VT = he_mass_VT->GetFunction("gaus");
    he_mass_MS->Fit("gaus", "", "", Mdileton-deltaM/2., Mdileton+ deltaM/2.);
    TF1 *fitFunction_MS = he_mass_MS->GetFunction("gaus");
    he_mass_Ref->Fit("gaus", "", "", Mdileton-deltaM/2., Mdileton+ deltaM/2.);
    TF1 *fitFunction_Ref = he_mass_Ref->GetFunction("gaus");
    he_mass_Ref_max->Fit("gaus", "", "", Mdileton-deltaM/2., Mdileton+ deltaM/2.);
    TF1 *fitFunction_Ref_max = he_mass_Ref_max->GetFunction("gaus");
    he_mass_Ref_min->Fit("gaus", "", "", Mdileton-deltaM/2., Mdileton+ deltaM/2.);
    TF1 *fitFunction_Ref_min = he_mass_Ref_min->GetFunction("gaus");
    std::cout<<"FIT RESULTS"<<std::endl;
    std::cout<<"Best m resolution: "<<fitFunction_best->GetParameter(2)*1000.<<" +- "<<fitFunction_best->GetParError(2)*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"VT m resolution: "<<fitFunction_VT->GetParameter(2)*1000.<<" +- "<<fitFunction_VT->GetParError(2)*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"MS m resolution: "<<fitFunction_MS->GetParameter(2)*1000.<<" +- "<<fitFunction_MS->GetParError(2)*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"VT+MS m resolution: "<<fitFunction_Ref->GetParameter(2)*1000.<<" +- "<<fitFunction_Ref->GetParError(2)*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"VT+MS m resolution with small E loss: "<<fitFunction_Ref_max->GetParameter(2)*1000.<<" +- "<<fitFunction_Ref_max->GetParError(2)*1000.<<" MeV/c^2"<<std::endl;
    std::cout<<"VT+MS m resolution with large E loss: "<<fitFunction_Ref_min->GetParameter(2)*1000.<<" +- "<<fitFunction_Ref_min->GetParError(2)*1000.<<" MeV/c^2"<<std::endl;
    

    f->Close();
    /*
    */
}

void plotRefit()
{
    CheckRefit("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_omega_full_muons_maxSeedsPerSpM_primary_1_maxSeedsPerSpM_secondary_20_numMeasurementsCutOff_1_chi2CutOff_15_fatras/tracksummary_matched.root",
                "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_omega_full_muons_maxSeedsPerSpM_primary_1_maxSeedsPerSpM_secondary_20_numMeasurementsCutOff_1_chi2CutOff_15_fatras/tracksummary_matchedVT.root",
                "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_omega_full_muons_maxSeedsPerSpM_primary_1_maxSeedsPerSpM_secondary_20_numMeasurementsCutOff_1_chi2CutOff_15_fatras/tracksummary_matchedMS.root", "tracksummary", "omega", 223);
    
    CheckRefit("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_jpsi_full_muons_maxSeedsPerSpM_primary_1_maxSeedsPerSpM_secondary_20_numMeasurementsCutOff_1_chi2CutOff_15_fatras/tracksummary_matched.root",
                "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_jpsi_full_muons_maxSeedsPerSpM_primary_1_maxSeedsPerSpM_secondary_20_numMeasurementsCutOff_1_chi2CutOff_15_fatras/tracksummary_matchedVT.root",
                "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_jpsi_full_muons_maxSeedsPerSpM_primary_1_maxSeedsPerSpM_secondary_20_numMeasurementsCutOff_1_chi2CutOff_15_fatras/tracksummary_matchedMS.root", "tracksummary", "jpsi", 443);
    
    
}