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
gStyle->SetHistLineWidth(1);
gStyle->SetOptTitle(0);   // Show mean (1), RMS (1), but hide entries (0) and title (0)
//gStyle->SetOptStat(1100); // Show mean (1), RMS (1), but hide entries (0) and title (0)
//gStyle->SetCanvasDefW(1800);  // Set default canvas width to 800 pixels
//gStyle->SetCanvasDefH(1200);  // Set default canvas height to 600 pixels

TH1D* ProjectMean(const TH2D* h2, TString suffix, TString xtitle="", TString ytitle="", TString title ="")
{
    if (!h2) return nullptr;
    
    int nBinsX = h2->GetNbinsX();
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();
    
    TH1D* hMean = new TH1D(Form("hMean%s",suffix.Data()), Form("%s;%s;%s",title.Data(),xtitle.Data(),ytitle.Data()), nBinsX, xMin, xMax);
    
    for (int i = 1; i <= nBinsX; ++i) {
        TH1D* proj = h2->ProjectionY("_py", i, i);
        if (proj->GetEntries() > 0) {
          double mean = proj->GetMean();
          double meanErr = proj->GetMeanError();
            hMean->SetBinContent(i, mean);
            hMean->SetBinError(i, meanErr);
        } else {
            hMean->SetBinContent(i, 0);
        }
        delete proj; // Clean up the temporary histogram
    }
    
    return hMean;
}

TH1D* ProjectRMS(const TH2D* h2, TString suffix, TString xtitle="", TString ytitle="", TString title ="")
{
    if (!h2) return nullptr;
    
    int nBinsX = h2->GetNbinsX();
    double xMin = h2->GetXaxis()->GetXmin();
    double xMax = h2->GetXaxis()->GetXmax();
    
    TH1D* hRMS = new TH1D(Form("hRMS%s",suffix.Data()), Form("%s;%s;%s",title.Data(),xtitle.Data(),ytitle.Data()), nBinsX, xMin, xMax);
    
    for (int i = 1; i <= nBinsX; ++i) {
        TH1D* proj = h2->ProjectionY("_py", i, i);
        if (proj->GetEntries() > 0) {
          double rms = proj->GetRMS();
          double rmsErr = proj->GetRMSError();
            hRMS->SetBinContent(i, rms);
            hRMS->SetBinError(i, rmsErr);
        } else {
            hRMS->SetBinContent(i, 0);
        }
        delete proj; // Clean up the temporary histogram
    }
    
    return hRMS;
}


void CheckJPsi(const char *filename = "tracksummary_ambims.root", const char *treename = "tracksummary", const char *outputname = "jpsi", int pdgCode = 223)
{

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
  TTree *tree = (TTree *)file->Get(treename);
  if (!tree)
  {
    std::cerr << "Error: Could not find tree " << treename << " in file " << filename << std::endl;
    file->Close();
    return;
  }

  //---------------------------------------
  // create histos
  //---------------------------------------
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  // Retrieve the particle from the database
  TParticlePDG *particle = pdgDB->GetParticle(pdgCode);
  double Mdileton = particle->Mass();
  double deltaM = 0.2;
  if(pdgCode==223)
    deltaM = 0.7;

  std::string rootname = std::string("invmass")+outputname+std::string(".root");
  TFile *fout = new TFile(rootname.c_str(), "recreate");
  TH1D *he_mass_fit = new TH1D("he_mass_fit", ";M_{#mu^{+}#mu^{-}} (GeV/#it{c}^{2});dN/dM", 150, Mdileton*(1-deltaM), Mdileton*(1+deltaM));
  TH2D *he_mass_fit_vs_pt = new TH2D("he_mass_fit_vs_pt", ";#it{p}_{T} (GeV/#it{c};M_{#mu^{+}#mu^{-}} (GeV/#it{c}^{2});dN/dM", 40,0,3, 150, Mdileton*(1-deltaM), Mdileton*(1+deltaM));
  TH2D *he_mass_fit_vs_y = new TH2D("he_mass_fit_vs_y", ";#it{y};M_{#mu^{+}#mu^{-}} (GeV/#it{c}^{2});dN/dM", 40,1.5,5.5, 150, Mdileton*(1-deltaM), Mdileton*(1+deltaM));

  float ptMin = 0;
  float ptMax = 3;
  float yMin = 0.;
  float yMax = 5.5;
  float pMin = 0;
  float pMax = 30;

  TH1D *he_pullLoc0 = new TH1D("he_pullLoc0", ";pull d0;Entries", 200, -10, 10);
  TH1D *he_pullLoc1 = new TH1D("he_pullLoc1", ";pull z0;Entries", 200, -10, 10);
  TH1D *he_pullPhi = new TH1D("he_pullPhi", ";pull phi;Entries", 200, -10, 10);
  TH1D *he_pullTheta = new TH1D("he_pullTheta", ";pull theta;Entries", 200, -10, 10);
  TH1D *he_pullQOP = new TH1D("he_pullQOP", ";pull QoP;Entries", 200, -10, 10);
  TH1D *he_pullP = new TH1D("he_pullP", ";pull P;Entries", 200, -10, 10);
  TH1D *he_resPoverP = new TH1D("he_resPoverP", ";(#it{p}_{rec}-#it{p}_{gen})/#it{p}_{gen};Entries", 200, -1, 1);
  TH1D *he_resPxoverPx = new TH1D("he_resPxoverPx", ";(#it{p}_{x rec}-#it{p}_{x gen})/#it{p}_{x gen};Entries", 200, -1, 1);
  TH1D *he_resPyoverPy = new TH1D("he_resPyoverPy", ";(#it{p}_{y rec}-#it{p}_{y gen})/#it{p}_{y gen};Entries", 200, -1, 1);
  TH1D *he_resPzoverPz = new TH1D("he_resPzoverPz", ";(#it{p}_{z rec}-#it{p}_{z gen})/#it{p}_{z gen};Entries", 200, -1, 1);

  TH2D *he_pullLoc0_vs_p = new TH2D("he_pullLoc0_vs_p", ";#it{p} (GeV/#it{c};pull d0;Entries", 40,pMin,pMax,200, -10, 10);
  TH2D *he_pullLoc1_vs_p = new TH2D("he_pullLoc1_vs_p", ";#it{p} (GeV/#it{c};pull z0;Entries", 40,pMin,pMax,200, -10, 10);
  TH2D *he_pullPhi_vs_p = new TH2D("he_pullPhi_vs_p", ";#it{p} (GeV/#it{c};pull phi;Entries", 40,pMin,pMax,200, -10, 10);
  TH2D *he_pullTheta_vs_p = new TH2D("he_pullTheta_vs_p", ";#it{p} (GeV/#it{c};pull theta;Entries", 40,pMin,pMax,200, -10, 10);
  TH2D *he_pullQOP_vs_p = new TH2D("he_pullQOP_vs_p", ";#it{p} (GeV/#it{c};pull QoP;Entries", 40,pMin,pMax,200, -10, 10);
  TH2D *he_pullP_vs_p = new TH2D("he_pullP_vs_p", ";#it{p} (GeV/#it{c};pull P;Entries", 40,pMin,pMax,200, -10, 10);

  TH2D *he_pullLoc0_vs_pt = new TH2D("he_pullLoc0_vs_pt", ";#it{p}_{T} (GeV/#it{c};pull d0;Entries", 40,ptMin,ptMax,200, -10, 10);
  TH2D *he_pullLoc1_vs_pt = new TH2D("he_pullLoc1_vs_pt", ";#it{p}_{T} (GeV/#it{c};pull z0;Entries", 40,ptMin,ptMax,200, -10, 10);
  TH2D *he_pullPhi_vs_pt = new TH2D("he_pullPhi_vs_pt", ";#it{p}_{T} (GeV/#it{c};pull phi;Entries", 40,ptMin,ptMax,200, -10, 10);
  TH2D *he_pullTheta_vs_pt = new TH2D("he_pullTheta_vs_pt", ";#it{p}_{T} (GeV/#it{c};pull theta;Entries", 40,ptMin,ptMax,200, -10, 10);
  TH2D *he_pullQOP_vs_pt = new TH2D("he_pullQOP_vs_pt", ";#it{p}_{T} (GeV/#it{c};pull QoP;Entries", 40,ptMin,ptMax,200, -10, 10);
  TH2D *he_pullP_vs_pt = new TH2D("he_pullP_vs_pt", ";#it{p}_{T} (GeV/#it{c};pull P;Entries", 40,ptMin,ptMax,200, -10, 10);

  TH2D *he_pullLoc0_vs_y = new TH2D("he_pullLoc0_vs_y", ";#it{y};pull d0;Entries", 40,yMin,yMax, 200, -10, 10);
  TH2D *he_pullLoc1_vs_y = new TH2D("he_pullLoc1_vs_y", ";#it{y};pull z0;Entries", 40,yMin,yMax, 200, -10, 10);
  TH2D *he_pullPhi_vs_y = new TH2D("he_pullPhi_vs_y", ";#it{y};pull phi;Entries", 40,yMin,yMax, 200, -10, 10);
  TH2D *he_pullTheta_vs_y = new TH2D("he_pullTheta_vs_y", ";#it{y};pull theta;Entries", 40,yMin,yMax, 200, -10, 10);
  TH2D *he_pullQOP_vs_y = new TH2D("he_pullQOP_vs_y", ";#it{y};pull QoP;Entries", 40,yMin,yMax, 200, -10, 10);
  TH2D *he_pullP_vs_y = new TH2D("he_pullP_vs_y", ";#it{y};pull P;Entries", 40,yMin,yMax, 200, -10, 10);
  //TH1D *he_mass_fit = new TH1D("he_mass_fit", ";M_{#mu^{+}#mu^{-}} (GeV/#it{c}^{2});dN/dM", 150, 0.2,0.9);
  //---------------------------------------
  // read variables
  //----------------------------//----------
  // Define variables to hold kinematic data of true particles
  std::vector<float> *t_px = new std::vector<float>;
  std::vector<float> *t_py = nullptr;
  std::vector<float> *t_pz = nullptr;
  std::vector<float> *t_p = nullptr;
  std::vector<int> *t_charge = new std::vector<int>;
  std::vector<int> *t_nmu = new std::vector<int>;

  // Define variables to hold kinematic data of fitted particles
  std::vector<float> *eLOC0_pull = new std::vector<float>;
  std::vector<float> *eLOC1_pull = new std::vector<float>;
  std::vector<float> *ePHI_pull = new std::vector<float>;
  std::vector<float> *eTHETA_pull = new std::vector<float>;
  std::vector<float> *eQOP_pull = new std::vector<float>;

  std::vector<float> *eLOC0_fit = new std::vector<float>;
  std::vector<float> *eLOC1_fit = new std::vector<float>;
  std::vector<float> *ePHI_fit = new std::vector<float>;
  std::vector<float> *eTHETA_fit = new std::vector<float>;
  std::vector<float> *eQOP_fit = new std::vector<float>;
  std::vector<float> *err_eQOP_fit = new std::vector<float>;

  // Set branch addresses for truth particle
  tree->SetBranchAddress("t_px", &t_px);
  tree->SetBranchAddress("t_py", &t_py);
  tree->SetBranchAddress("t_pz", &t_pz);
  tree->SetBranchAddress("t_p", &t_p);
  tree->SetBranchAddress("t_charge", &t_charge);
  tree->SetBranchAddress("track_nr", &t_nmu);

  // set branche address for reco particle
  tree->SetBranchAddress("pull_eLOC0_fit", &eLOC0_pull);
  tree->SetBranchAddress("pull_eLOC1_fit", &eLOC1_pull);
  tree->SetBranchAddress("pull_ePHI_fit", &ePHI_pull);
  tree->SetBranchAddress("pull_eTHETA_fit", &eTHETA_pull);
  tree->SetBranchAddress("pull_eQOP_fit", &eQOP_pull);

  tree->SetBranchAddress("eLOC0_fit", &eLOC0_fit);
  tree->SetBranchAddress("eLOC1_fit", &eLOC1_fit);
  tree->SetBranchAddress("ePHI_fit", &ePHI_fit);
  tree->SetBranchAddress("eTHETA_fit", &eTHETA_fit);
  tree->SetBranchAddress("eQOP_fit", &eQOP_fit);
  tree->SetBranchAddress("err_eQOP_fit", &err_eQOP_fit);

  // Create TLorentzVector objects for the muons
  TLorentzVector t_muon1, t_muon2;
  TLorentzVector e_muon1_fit, e_muon2_fit;

  Double_t mass = 0.1056583755;

  //---------------------------------------
  // loop over entries
  //---------------------------------------

  Long64_t nentries = tree->GetEntries();
  Double_t onemuon = 0;
  Double_t zeromuon = 0;
  Double_t twomuons = 0;

  for (Long64_t i = 0; i < nentries; ++i)
  {
    tree->GetEntry(i);
    // std::cout << "Entry " << i << " number of muons= " << t_nmu->size() << std::endl;

    int sizz = (int)t_nmu->size();
    for (int im = 0; im < t_nmu->size(); im++){
      TLorentzVector t_helper;
      float sigmap = err_eQOP_fit->at(im)*(1./TMath::Abs(eQOP_fit->at(im))*1./TMath::Abs(eQOP_fit->at(im)));
      Double_t t_E = TMath::Sqrt(t_px->at(im) * t_px->at(im) + t_py->at(im) * t_py->at(im) + t_pz->at(im) * t_pz->at(im) + mass * mass);
      t_helper.SetPxPyPzE(t_px->at(im), t_py->at(im), t_pz->at(im), t_E);

      //if(t_helper.Pt() < 0.5) continue;
      he_pullLoc0->Fill(eLOC0_pull->at(im));
      he_pullLoc1->Fill(eLOC1_pull->at(im));
      he_pullTheta->Fill(ePHI_pull->at(im));
      he_pullPhi->Fill(eTHETA_pull->at(im));
      he_pullQOP->Fill(eQOP_pull->at(im));
      he_pullP->Fill((1./TMath::Abs(eQOP_fit->at(im))-t_p->at(im))/sigmap);
      he_resPoverP->Fill((1./TMath::Abs(eQOP_fit->at(im))-t_p->at(im))/t_p->at(im));

      float px = 1./std::abs(eQOP_fit->at(im))*std::sin(eTHETA_fit->at(im))*std::cos(ePHI_fit->at(im));
      float py = 1./std::abs(eQOP_fit->at(im))*std::sin(eTHETA_fit->at(im))*std::sin(ePHI_fit->at(im));
      float pz = 1./std::abs(eQOP_fit->at(im))*std::cos(eTHETA_fit->at(im));

      he_resPxoverPx->Fill((px-t_px->at(im))/t_px->at(im));
      he_resPyoverPy->Fill((py-t_py->at(im))/t_py->at(im));
      he_resPzoverPz->Fill((pz-t_pz->at(im))/t_pz->at(im));


      he_pullLoc0_vs_p->Fill(t_helper.P(), eLOC0_pull->at(im));
      he_pullLoc1_vs_p->Fill(t_helper.P(), eLOC1_pull->at(im));
      he_pullTheta_vs_p->Fill(t_helper.P(), ePHI_pull->at(im));
      he_pullPhi_vs_p->Fill(t_helper.P(), eTHETA_pull->at(im));
      he_pullQOP_vs_p->Fill(t_helper.P(), eQOP_pull->at(im));
      he_pullP_vs_p->Fill(t_helper.P(), (1./TMath::Abs(eQOP_fit->at(im))-t_p->at(im))/sigmap);

      he_pullLoc0_vs_pt->Fill(t_helper.Pt(), eLOC0_pull->at(im));
      he_pullLoc1_vs_pt->Fill(t_helper.Pt(), eLOC1_pull->at(im));
      he_pullTheta_vs_pt->Fill(t_helper.Pt(), ePHI_pull->at(im));
      he_pullPhi_vs_pt->Fill(t_helper.Pt(), eTHETA_pull->at(im));
      he_pullQOP_vs_pt->Fill(t_helper.Pt(), eQOP_pull->at(im));
      he_pullP_vs_pt->Fill(t_helper.Pt(), (1./TMath::Abs(eQOP_fit->at(im))-t_p->at(im))/sigmap);

      he_pullLoc0_vs_y->Fill(t_helper.Rapidity(), eLOC0_pull->at(im));
      he_pullLoc1_vs_y->Fill(t_helper.Rapidity(), eLOC1_pull->at(im));
      he_pullTheta_vs_y->Fill(t_helper.Rapidity(), ePHI_pull->at(im));
      he_pullPhi_vs_y->Fill(t_helper.Rapidity(), eTHETA_pull->at(im));
      he_pullQOP_vs_y->Fill(t_helper.Rapidity(), eQOP_pull->at(im));
      he_pullP_vs_y->Fill(t_helper.Rapidity(), (1./TMath::Abs(eQOP_fit->at(im))-t_p->at(im))/sigmap);
    }

    if (t_nmu->size() < 1)
    {
      zeromuon++;
      continue;
    }

    // skip events where there is only one muon
    if (t_nmu->size() < 2)
    {
      if (t_nmu->size() == 1)
        onemuon++;
      continue;
    }
    twomuons++;

    for (int im = 0; im < t_nmu->size() - 1; im++)
    {
      Double_t t_E = TMath::Sqrt(t_px->at(im) * t_px->at(im) + t_py->at(im) * t_py->at(im) + t_pz->at(im) * t_pz->at(im) + mass * mass);

      t_muon1.SetPxPyPzE(t_px->at(im), t_py->at(im), t_pz->at(im), t_E);
      //if(t_muon1.Pt() < 0.5) continue;
      // Calculate momentum for reco particles
      double p1 = 1.0 / std::abs(eQOP_fit->at(im));
      double px1 = p1 * std::sin(eTHETA_fit->at(im)) * std::cos(ePHI_fit->at(im));
      double py1 = p1 * std::sin(eTHETA_fit->at(im)) * std::sin(ePHI_fit->at(im));
      double pz1 = p1 * std::cos(eTHETA_fit->at(im));
      double E1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + mass * mass);
      e_muon1_fit.SetPxPyPzE(px1, py1, pz1, E1);


      for (int im2 = im + 1; im2 < t_nmu->size(); im2++)
      {

        Double_t t_E = TMath::Sqrt(t_px->at(im2) * t_px->at(im2) + t_py->at(im2) * t_py->at(im2) + t_pz->at(im2) * t_pz->at(im2) + mass * mass);
        t_muon2.SetPxPyPzE(t_px->at(im2), t_py->at(im2), t_pz->at(im2), t_E);

        //if(t_muon2.Pt() < 0.5) continue;
        double t_invariant_mass = (t_muon1 + t_muon2).M();

        double p2 = 1.0 / std::abs(eQOP_fit->at(im2));
        double px2 = p2 * std::sin(eTHETA_fit->at(im2)) * std::cos(ePHI_fit->at(im2));
        double py2 = p2 * std::sin(eTHETA_fit->at(im2)) * std::sin(ePHI_fit->at(im2));
        double pz2 = p2 * std::cos(eTHETA_fit->at(im2));
        double E2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + mass * mass);
        e_muon2_fit.SetPxPyPzE(px2, py2, pz2, E2);

        double e_invariant_mass_fit = (e_muon1_fit + e_muon2_fit).M();
        he_mass_fit->Fill(e_invariant_mass_fit);
        std::cout<<e_invariant_mass_fit<<std::endl;
        he_mass_fit_vs_y->Fill((e_muon1_fit + e_muon2_fit).Rapidity(), e_invariant_mass_fit);
        he_mass_fit_vs_pt->Fill((e_muon1_fit + e_muon2_fit).Pt(), e_invariant_mass_fit);
      }
    }
  }

  std::cout << " ****************************************** " << std::endl;
  std::cout << "   NEntries= " << nentries << std::endl;
  std::cout << "   events with zero muons= " << zeromuon << std::endl;
  std::cout << "   events with one muon= " << onemuon << std::endl;
  std::cout << "   events with two muons= " << twomuons << std::endl;
  //--------------------------------------
  // plot histos
  //--------------------------------------

  TCanvas *c10 = new TCanvas("c10", "Inv mass", 1600, 1200);
  c10->Update();
  gStyle->SetOptFit(0);

  gStyle->SetTitleSize(0.05, "XYZ"); // Sets the size for X, Y, and Z axis titles
  gStyle->SetTitleFontSize(0.05);    // Sets the size for the histogram title
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetPadLeftMargin(0.15);   // Set left margin
  gStyle->SetPadRightMargin(0.05);  // Set right margin
  gStyle->SetPadTopMargin(0.05);    // Set top margin
  gStyle->SetPadBottomMargin(0.15); // Set bottom margin
  gStyle->SetHistLineWidth(1);
  gStyle->SetOptTitle(0);   // Show mean (1), RMS (1), but hide entries (0) and title (0)
  //gStyle->SetOptStat(0); // Show mean (1), RMS (1), but hide entries (0) and title (0)

  //he_mass_fit->Rebin(2);
  he_mass_fit->SetMarkerColor(8);
  he_mass_fit->SetMarkerStyle(20);
  he_mass_fit->SetMarkerSize(2);
  he_mass_fit->SetLineWidth(1);
//  he_mass_fit->Fit("gaus", "", "", 3, 3.3);
  //he_mass_fit->Fit("gaus","","MR",Mdileton-0.1,Mdileton+0.1);//, "", "", 2.9, 3.3);
  //TF1 *fitFunction = he_mass_fit->GetFunction("gaus");
  //fitFunction->SetLineWidth(1);
  //fitFunction->SetNpx(1000);
  he_mass_fit->Draw("P");

  fout->cd();
  he_mass_fit->Write();
  he_mass_fit_vs_y->Write();
  he_mass_fit_vs_pt->Write();
  TLatex latex;
  latex.SetTextSize(0.06);         // Sets the text size
  latex.SetTextFont(42);           // Sets the text font
  latex.SetTextColor(kBlack);      // Sets the text color (optional)
  latex.SetTextAlign(13);   // Align at the top-left corner (1=left, 3=top)
  
  // Draw LaTeX text
  //latex.DrawLatex(3.5, he_mass_fit->GetMaximum()*0.85, Form("#sigma_{M} = %.0f MeV/#it{c}^{2}", fitFunction->GetParameter(2)*1000.));
  //latex.DrawLatex(0.52, he_mass_fit->GetMaximum()*0.85, Form("#sigma_{M} = %.0f MeV/#it{c}^{2}", fitFunction->GetParameter(2)*1000.));


  std::string pngname = outputname+std::string(".png");
  c10->Update();
  c10->SaveAs(pngname.c_str());

  TH1D* he_pullLoc0_vs_y_mean = ProjectMean(he_pullLoc0_vs_y, "_loc0_y","#it{y}; Mean pull_{Loc0}"); 
  TH1D* he_pullLoc1_vs_y_mean = ProjectMean(he_pullLoc1_vs_y, "_loc1_y","#it{y}; Mean pull_{Loc1}"); 
  TH1D* he_pullPhi_vs_y_mean = ProjectMean(he_pullPhi_vs_y, "_phi_y","#it{y}; Mean pull_{#phi}"); 
  TH1D* he_pullTheta_vs_y_mean = ProjectMean(he_pullTheta_vs_y, "_theta_y","#it{y}; Mean pull_{#theta}");
  TH1D* he_pullQOP_vs_y_mean = ProjectMean(he_pullQOP_vs_y, "_qop_y","#it{y}; Mean pull_{QoP}"); 

  TH1D* he_pullLoc0_vs_y_rms = ProjectRMS(he_pullLoc0_vs_y, "_loc0_y","#it{y}; RMS pull_{Loc0}"); 
  TH1D* he_pullLoc1_vs_y_rms = ProjectRMS(he_pullLoc1_vs_y, "_loc1_y","#it{y}; RMS pull_{Loc1}"); 
  TH1D* he_pullPhi_vs_y_rms = ProjectRMS(he_pullPhi_vs_y, "_phi_y","#it{y}; RMS pull_{#phi}"); 
  TH1D* he_pullTheta_vs_y_rms = ProjectRMS(he_pullTheta_vs_y, "_theta_y","#it{y}; RMS pull_{#theta}");
  TH1D* he_pullQOP_vs_y_rms = ProjectRMS(he_pullQOP_vs_y, "_qop_y","#it{y}; RMS pull_{QoP}"); 

  TH1D* he_pullLoc0_vs_pt_mean = ProjectMean(he_pullLoc0_vs_pt, "_loc0_pt","#it{p}_{T} (GeV/#it{c}); Mean pull_{Loc0}"); 
  TH1D* he_pullLoc1_vs_pt_mean = ProjectMean(he_pullLoc1_vs_pt, "_loc1_pt","#it{p}_{T} (GeV/#it{c}); Mean pull_{Loc1}"); 
  TH1D* he_pullPhi_vs_pt_mean = ProjectMean(he_pullPhi_vs_pt, "_phi_pt","#it{p}_{T} (GeV/#it{c}); Mean pull_{#phi}"); 
  TH1D* he_pullTheta_vs_pt_mean = ProjectMean(he_pullTheta_vs_pt, "_theta_pt","#it{p}_{T} (GeV/#it{c}); Mean pull_{#theta}");
  TH1D* he_pullQOP_vs_pt_mean = ProjectMean(he_pullQOP_vs_pt, "_qop_pt","#it{p}_{T} (GeV/#it{c}); Mean pull_{QoP}"); 

  TH1D* he_pullLoc0_vs_pt_rms = ProjectRMS(he_pullLoc0_vs_pt, "_loc0_pt","#it{p}_{T} (GeV/#it{c}); RMS pull_{Loc0}"); 
  TH1D* he_pullLoc1_vs_pt_rms = ProjectRMS(he_pullLoc1_vs_pt, "_loc1_pt","#it{p}_{T} (GeV/#it{c}); RMS pull_{Loc1}"); 
  TH1D* he_pullPhi_vs_pt_rms = ProjectRMS(he_pullPhi_vs_pt, "_phi_pt","#it{p}_{T} (GeV/#it{c}); RMS pull_{#phi}"); 
  TH1D* he_pullTheta_vs_pt_rms = ProjectRMS(he_pullTheta_vs_pt, "_theta_pt","#it{p}_{T} (GeV/#it{c}); RMS pull_{#theta}");
  TH1D* he_pullQOP_vs_pt_rms = ProjectRMS(he_pullQOP_vs_pt, "_qop_pt","#it{p}_{T} (GeV/#it{c}); RMS pull_{QoP}"); 

  TH1D* he_pullLoc0_vs_p_mean = ProjectMean(he_pullLoc0_vs_p, "_loc0_p","#it{p} (GeV/#it{c}); Mean pull_{Loc0}"); 
  TH1D* he_pullLoc1_vs_p_mean = ProjectMean(he_pullLoc1_vs_p, "_loc1_p","#it{p} (GeV/#it{c}); Mean pull_{Loc1}"); 
  TH1D* he_pullPhi_vs_p_mean = ProjectMean(he_pullPhi_vs_p, "_phi_p","#it{p} (GeV/#it{c}); Mean pull_{#phi}"); 
  TH1D* he_pullTheta_vs_p_mean = ProjectMean(he_pullTheta_vs_p, "_theta_p","#it{p} (GeV/#it{c}); Mean pull_{#theta}");
  TH1D* he_pullQOP_vs_p_mean = ProjectMean(he_pullQOP_vs_p, "_qop_p","#it{p} (GeV/#it{c}); Mean pull_{QoP}"); 

  TH1D* he_pullLoc0_vs_p_rms = ProjectRMS(he_pullLoc0_vs_p, "_loc0_p","#it{p} (GeV/#it{c}); RMS pull_{Loc0}"); 
  TH1D* he_pullLoc1_vs_p_rms = ProjectRMS(he_pullLoc1_vs_p, "_loc1_p","#it{p} (GeV/#it{c}); RMS pull_{Loc1}"); 
  TH1D* he_pullPhi_vs_p_rms = ProjectRMS(he_pullPhi_vs_p, "_phi_p","#it{p} (GeV/#it{c}); RMS pull_{#phi}"); 
  TH1D* he_pullTheta_vs_p_rms = ProjectRMS(he_pullTheta_vs_p, "_theta_p","#it{p} (GeV/#it{c}); RMS pull_{#theta}");
  TH1D* he_pullQOP_vs_p_rms = ProjectRMS(he_pullQOP_vs_p, "_qop_p","#it{p} (GeV/#it{c}); RMS pull_{QoP}"); 

  he_pullLoc0->Write(); 
  he_pullLoc1->Write(); 
  he_pullPhi->Write(); 
  he_pullTheta->Write();
  he_pullQOP->Write(); 
  he_pullP->Write(); 
  he_resPoverP->Write();
  he_resPzoverPz->Write();
  he_resPxoverPx->Write();
  he_resPyoverPy->Write();

  
  TCanvas *cvRes = new TCanvas("cvRes", "", 1920*2, 1280*2);

  gPad->SetLogy();
  cvRes->Divide(3,1);
  cvRes->Update();
  cvRes->cd(1);
  gPad->SetLogy();
  he_resPxoverPx->Draw();
  cvRes->cd(2);
  gPad->SetLogy();
  he_resPyoverPy->Draw();
  cvRes->cd(3);
  gPad->SetLogy();
  he_resPzoverPz->Draw();

  cvRes->Update();
  pngname = outputname+std::string("_resp.png");
  cvRes->SaveAs(pngname.c_str());

  he_pullLoc0_vs_y->Write(); 
  he_pullLoc1_vs_y->Write(); 
  he_pullPhi_vs_y->Write(); 
  he_pullTheta_vs_y->Write();
  he_pullQOP_vs_y->Write(); 
  he_pullP_vs_y->Write(); 

  he_pullLoc0_vs_y_mean->Write(); 
  he_pullLoc1_vs_y_mean->Write(); 
  he_pullPhi_vs_y_mean->Write(); 
  he_pullTheta_vs_y_mean->Write();
  he_pullQOP_vs_y_mean->Write(); 

  he_pullLoc0_vs_y_rms->Write(); 
  he_pullLoc1_vs_y_rms->Write(); 
  he_pullPhi_vs_y_rms->Write(); 
  he_pullTheta_vs_y_rms->Write();
  he_pullQOP_vs_y_rms->Write(); 

  he_pullLoc0_vs_pt->Write(); 
  he_pullLoc1_vs_pt->Write(); 
  he_pullPhi_vs_pt->Write(); 
  he_pullTheta_vs_pt->Write();
  he_pullQOP_vs_pt->Write(); 
  he_pullP_vs_pt->Write(); 
  

  he_pullLoc0_vs_pt_mean->Write(); 
  he_pullLoc1_vs_pt_mean->Write(); 
  he_pullPhi_vs_pt_mean->Write(); 
  he_pullTheta_vs_pt_mean->Write();
  he_pullQOP_vs_pt_mean->Write(); 

  he_pullLoc0_vs_pt_rms->Write(); 
  he_pullLoc1_vs_pt_rms->Write(); 
  he_pullPhi_vs_pt_rms->Write(); 
  he_pullTheta_vs_pt_rms->Write();
  he_pullQOP_vs_pt_rms->Write(); 

  he_pullLoc0_vs_p->Write(); 
  he_pullLoc1_vs_p->Write(); 
  he_pullPhi_vs_p->Write(); 
  he_pullTheta_vs_p->Write();
  he_pullQOP_vs_p->Write(); 
  he_pullP_vs_p->Write(); 
  

  he_pullLoc0_vs_p_mean->Write(); 
  he_pullLoc1_vs_p_mean->Write(); 
  he_pullPhi_vs_p_mean->Write(); 
  he_pullTheta_vs_p_mean->Write();
  he_pullQOP_vs_p_mean->Write(); 

  he_pullLoc0_vs_p_rms->Write(); 
  he_pullLoc1_vs_p_rms->Write(); 
  he_pullPhi_vs_p_rms->Write(); 
  he_pullTheta_vs_p_rms->Write();
  he_pullQOP_vs_p_rms->Write(); 

  ////////////////////////////////////////////
  TLine *lineMeanPt = new TLine(ptMin, 0, ptMax, 0);
  TLine *lineMeanP = new TLine(pMin, 0, pMax, 0);
  TLine *lineMeanY = new TLine(yMin, 0, yMax, 0);
  TLine *lineRMSPt = new TLine(ptMin, 1, ptMax, 1);
  TLine *lineRMSP = new TLine(pMin, 1, pMax, 1);
  TLine *lineRMSY = new TLine(yMin, 1, yMax, 1);
  lineMeanPt->SetLineColor(kRed);
  lineMeanP->SetLineColor(kRed);
  lineMeanY->SetLineColor(kRed);
  lineRMSPt->SetLineColor(kRed);
  lineRMSP->SetLineColor(kRed);
  lineRMSY->SetLineColor(kRed);
  gStyle->SetOptStat(0);
  
  TCanvas *cvDiv = new TCanvas("cvDiv", "", 1920*2, 1280*2);
  cvDiv->Divide(5,2);
  cvDiv->Update();
  cvDiv->cd(1);
  he_pullLoc0_vs_y_mean->Draw();
  lineMeanY->Draw("same");
  cvDiv->cd(2);
  he_pullLoc1_vs_y_mean->Draw();
  lineMeanY->Draw("same");
  cvDiv->cd(3);
  he_pullPhi_vs_y_mean->Draw();
  lineMeanY->Draw("same");
  cvDiv->cd(4);
  he_pullTheta_vs_y_mean->Draw();
  lineMeanY->Draw("same");
  cvDiv->cd(5);
  he_pullQOP_vs_y_mean->Draw();
  lineMeanY->Draw("same");
  cvDiv->cd(6);
  he_pullLoc0_vs_y_rms->Draw();
  lineRMSY->Draw("same");
  cvDiv->cd(7);
  he_pullLoc1_vs_y_rms->Draw();
  lineRMSY->Draw("same");
  cvDiv->cd(8);
  he_pullPhi_vs_y_rms->Draw();
  lineRMSY->Draw("same");
  cvDiv->cd(9);
  he_pullTheta_vs_y_rms->Draw();
  lineRMSY->Draw("same");
  cvDiv->cd(10);
  he_pullQOP_vs_y_rms->Draw();
  lineRMSY->Draw("same");
  cvDiv->Write();
  pngname = outputname+std::string("_y.png");
  cvDiv->SaveAs(pngname.c_str());
  pngname = outputname+std::string("_y.pdf");
  cvDiv->SaveAs(pngname.c_str());
  cvDiv->cd(1);
  he_pullLoc0_vs_pt_mean->Draw();
  lineMeanPt->Draw("same");
  cvDiv->cd(2);
  he_pullLoc1_vs_pt_mean->Draw();
  lineMeanPt->Draw("same");
  cvDiv->cd(3);
  he_pullPhi_vs_pt_mean->Draw();
  lineMeanPt->Draw("same");
  cvDiv->cd(4);
  he_pullTheta_vs_pt_mean->Draw();
  lineMeanPt->Draw("same");
  cvDiv->cd(5);
  he_pullQOP_vs_pt_mean->Draw();
  lineMeanPt->Draw("same");
  cvDiv->cd(6);
  he_pullLoc0_vs_pt_rms->Draw();
  lineRMSPt->Draw("same");
  cvDiv->cd(7);
  he_pullLoc1_vs_pt_rms->Draw();
  lineRMSPt->Draw("same");
  cvDiv->cd(8);
  he_pullPhi_vs_pt_rms->Draw();
  lineRMSPt->Draw("same");
  cvDiv->cd(9);
  he_pullTheta_vs_pt_rms->Draw();
  lineRMSPt->Draw("same");
  cvDiv->cd(10);
  he_pullQOP_vs_pt_rms->Draw();
  lineRMSPt->Draw("same");
  cvDiv->Write();
  pngname = outputname+std::string("_pt.png");
  cvDiv->SaveAs(pngname.c_str());
  pngname = outputname+std::string("_pt.pdf");
  cvDiv->SaveAs(pngname.c_str());
  cvDiv->cd(1);
  he_pullLoc0_vs_p_mean->Draw();
  lineMeanP->Draw("same");
  cvDiv->cd(2);
  he_pullLoc1_vs_p_mean->Draw();
  lineMeanP->Draw("same");
  cvDiv->cd(3);
  he_pullPhi_vs_p_mean->Draw();
  lineMeanP->Draw("same");
  cvDiv->cd(4);
  he_pullTheta_vs_p_mean->Draw();
  lineMeanP->Draw("same");
  cvDiv->cd(5);
  he_pullQOP_vs_p_mean->Draw();
  lineMeanP->Draw("same");
  cvDiv->cd(6);
  he_pullLoc0_vs_p_rms->Draw();
  lineRMSP->Draw("same");
  cvDiv->cd(7);
  he_pullLoc1_vs_p_rms->Draw();
  lineRMSP->Draw("same");
  cvDiv->cd(8);
  he_pullPhi_vs_p_rms->Draw();
  lineRMSP->Draw("same");
  cvDiv->cd(9);
  he_pullTheta_vs_p_rms->Draw();
  lineRMSP->Draw("same");
  cvDiv->cd(10);
  he_pullQOP_vs_p_rms->Draw();
  lineRMSP->Draw("same");
  cvDiv->Write();
  pngname = outputname+std::string("_p.png");
  cvDiv->SaveAs(pngname.c_str());
  pngname = outputname+std::string("_p.pdf");
  cvDiv->SaveAs(pngname.c_str());
  
  c10->cd();
  he_pullLoc0->Draw();
  pngname = std::string("he_pullLoc0")+outputname+std::string(".png");
  c10->SaveAs(pngname.c_str());

  he_pullLoc1->Draw();
  pngname = std::string("he_pullLoc1")+outputname+std::string(".png");
  c10->SaveAs(pngname.c_str());

  he_pullPhi->Draw();
  pngname = std::string("he_pullPhi")+outputname+std::string(".png");
  c10->SaveAs(pngname.c_str());

  he_pullTheta->Draw();
  pngname = std::string("he_pullTheta")+outputname+std::string(".png");
  c10->SaveAs(pngname.c_str());

  he_pullQOP->Draw();
  pngname = std::string("he_pullQOP")+outputname+std::string(".png");
  c10->SaveAs(pngname.c_str());
  fout->Close();
}

void plotJpsi()
{

  /*
  Omega 223
  Jpsi 443
  Eta 221
  Phi 333
  */
  gStyle->SetTitleSize(0.05, "XYZ"); // Sets the size for X, Y, and Z axis titles
  gStyle->SetTitleFontSize(0.05);    // Sets the size for the histogram title
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetPadLeftMargin(0.15);   // Set left margin
  gStyle->SetPadRightMargin(0.05);  // Set right margin
  gStyle->SetPadTopMargin(0.05);    // Set top margin
  gStyle->SetPadBottomMargin(0.15); // Set bottom margin
  gStyle->SetHistLineWidth(1);
  gStyle->SetOptTitle(0);   // Show mean (1), RMS (1), but hide entries (0) and title (0)
  //gStyle->SetOptStat(1100); // Show mean (1), RMS (1), but hide entries (0) and title (0)
  gROOT->SetBatch(kTRUE);
  //CheckJPsi("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_Sec_jpsi_muons/tracksummary_ambi.root", "tracksummary", "jpsiVT",443);
  //CheckJPsi("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_Sec_jpsi_muons/tracksummary_ambims.root", "tracksummary", "jpsiMS",443);
  //CheckJPsi("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_Sec_omega2Body_muons/tracksummary_ambi.root", "tracksummary", "omegatVT",223);
  //CheckJPsi("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_Sec_omega2Body_muons/tracksummary_ambims.root", "tracksummary", "omegatMS",223);
  
  TFile *foutOmegaMS = new TFile("invmassomegatMS.root", "read");
  TH1D *omegaMS = (TH1D*)foutOmegaMS->Get("he_mass_fit");
  TFile *foutOmegaVT = new TFile("invmassomegatVT.root", "read");
  TH1D *omegaVT = (TH1D*)foutOmegaVT->Get("he_mass_fit");

  TFile *foutjpsiMS = new TFile("invmassjpsiMS.root", "read");
  TH1D *jpsiMS = (TH1D*)foutjpsiMS->Get("he_mass_fit");
  TFile *foutjpsiVT = new TFile("invmassjpsiVT.root", "read");
  TH1D *jpsiVT = (TH1D*)foutjpsiVT->Get("he_mass_fit");

  // Create your own Gaussian functions to ensure they exist
  TF1* gausVT = new TF1("gausVT", "gaus", 3, 3.3);  // Adjust range as needed
  TF1* gausMS = new TF1("gausMS", "gaus", 3, 3.3);
  gStyle->SetOptStat(0);
  // Fit without drawing
  jpsiVT->Fit(gausVT, "RNQ");  // N = no draw, Q = quiet
  jpsiMS->Fit(gausMS, "RNQ");

  TCanvas *c1 = new TCanvas("c1", "Inv mass");
  jpsiVT->SetMarkerColor(kBlue);
  jpsiVT->SetLineColor(kBlue);
  jpsiVT->SetMarkerStyle(20);
  jpsiVT->SetMarkerSize(0.7);
  jpsiMS->SetMarkerColor(kRed);
  jpsiMS->SetLineColor(kRed);
  jpsiMS->SetMarkerStyle(21);
  jpsiMS->SetMarkerSize(0.7);
  jpsiVT->SetLineWidth(1);
  jpsiMS->SetLineWidth(1);

  jpsiMS->Scale(jpsiVT->Integral() / jpsiMS->Integral());
  jpsiVT->Draw("EP");
  jpsiMS->Draw("EP same");
  jpsiVT->Draw("HIST SAME ");
  jpsiMS->Draw("HIST same");

  // Dummy clones for legend with larger markers
  TH1D* jpsiVT_legend = (TH1D*)jpsiVT->Clone("jpsiVT_legend");
  TH1D* jpsiMS_legend = (TH1D*)jpsiMS->Clone("jpsiMS_legend");
  jpsiVT_legend->SetMarkerSize(1);  // larger marker for legend only
  jpsiMS_legend->SetMarkerSize(1);
  TLegend legend(0.55, 0.55, 0.9, 0.9);
  legend.SetTextSize(0.043);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetHeader("J/#psi#rightarrow#mu^{+}#mu^{-}", "C");  // "C" centers the title
  legend.AddEntry(jpsiVT_legend, Form("VS+MS #sigma_{M} = %.0f MeV/#it{c}^{2}", 1000 * gausVT->GetParameter(2)), "p");
  legend.AddEntry(jpsiMS_legend, Form("MS #sigma_{M} = %.0f MeV/#it{c}^{2}", 1000 * gausMS->GetParameter(2)), "p");
  legend.Draw();

  c1->SaveAs("ACTS_InvMass_jpsi.png");
  c1->SaveAs("ACTS_InvMass_jpsi.pdf");

  // Create your own Gaussian functions to ensure they exist
  TF1* gausVTo = new TF1("gausVTo", "gaus", 0.6, 1);  // Adjust range as needed
  TF1* gausMSo = new TF1("gausMSo", "gaus", 0.6, 1);

  // Fit without drawing
  omegaVT->Fit(gausVTo, "RNQ");  // N = no draw, Q = quiet
  omegaMS->Fit(gausMSo, "RNQ");

  omegaVT->SetMarkerColor(kBlue);
  omegaVT->SetLineColor(kBlue);
  omegaVT->SetLineWidth(1);
  omegaVT->SetMarkerStyle(20);
  omegaVT->SetMarkerSize(0.7);
  omegaMS->SetMarkerColor(kRed);
  omegaMS->SetLineColor(kRed);
  omegaMS->SetMarkerStyle(21);
  omegaMS->SetMarkerSize(0.7);
  omegaMS->SetLineWidth(1);

  // Dummy clones for legend with larger markers
  TH1D* omegaVT_legend = (TH1D*)omegaVT->Clone("omegaVT_legend");
  TH1D* omegaMS_legend = (TH1D*)omegaMS->Clone("omegaMS_legend");
  omegaVT_legend->SetMarkerSize(1);  // larger marker for legend only
  omegaMS_legend->SetMarkerSize(1);

  omegaMS->Scale(omegaVT->Integral() / omegaMS->Integral());
  omegaVT->GetXaxis()->SetRangeUser(0.4, 1.2);
  omegaVT->Draw("EP");
  omegaVT->Draw("HIST SAME");
  omegaMS->Draw("EP same");
  omegaMS->Draw("HIST same");
  TLegend legendOmega(0.55, 0.55, 0.9, 0.9);
  legendOmega.SetTextSize(0.045);
  legendOmega.SetBorderSize(0);
  legendOmega.SetFillStyle(0);
  legendOmega.SetHeader("#omega#rightarrow#mu^{+}#mu^{-}", "C");  // "C" centers the title
  legendOmega.AddEntry(omegaVT_legend, Form("VS+MS #sigma_{M} = %.0f MeV/#it{c}^{2}", 1000 * gausVTo->GetParameter(2)), "p");
  legendOmega.AddEntry(omegaMS_legend, Form("MS #sigma_{M} = %.0f MeV/#it{c}^{2}", 1000 * gausMSo->GetParameter(2)), "p");
  legendOmega.Draw();
  c1->SaveAs("ACTS_InvMass_omega.png");
  c1->SaveAs("ACTS_InvMass_omega.pdf");

}