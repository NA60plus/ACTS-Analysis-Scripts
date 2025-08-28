#include <TFile.h>
#include <TTree.h>
#include <iostream>

double GetFWHM(TH1F* hist) {
    int maxBin = hist->GetMaximumBin();  
    double maxVal = hist->GetBinContent(maxBin);  
    double halfMax = maxVal / 2.0;  

    int binLow = maxBin, binHigh = maxBin;

    // Find left half-maximum crossing
    while (binLow > 1 && hist->GetBinContent(binLow) > halfMax) {
        binLow--;
    }

    // Find right half-maximum crossing
    while (binHigh < hist->GetNbinsX() && hist->GetBinContent(binHigh) > halfMax) {
        binHigh++;
    }

    double xLow = hist->GetXaxis()->GetBinCenter(binLow);
    double xHigh = hist->GetXaxis()->GetBinCenter(binHigh);

    return xHigh - xLow;  // FWHM
}

int countOutsideRange(TH1F* hist, double x_min, double x_max) {
    int nBins = hist->GetNbinsX();
    int count = 0;

    for (int bin = 1; bin <= nBins; ++bin) {  // Bins start from 1
        double x_center = hist->GetXaxis()->GetBinCenter(bin);
        double bin_content = hist->GetBinContent(bin);

        if (x_center < x_min || x_center > x_max) {
            count += bin_content;
        }
    }
    
    return count;
}

int readELoss() {
    //gStyle->SetOptStat(0);
    // Open the ROOT file
    TFile *file = TFile::Open("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output__noAbs/particles_simulation_ms.root", "READ");
    // Get the TTree
    TTree *tree = (TTree *)file->Get("particles");
    // Variables to hold branch data
    UInt_t event_id;
    std::vector<unsigned long> *particle_id = nullptr;
    std::vector<float> *particle_type = nullptr;
    std::vector<float> *vx = nullptr;
    std::vector<float> *vy = nullptr;
    std::vector<float> *vz = nullptr;
    std::vector<float> *vt = nullptr;
    std::vector<float> *px = nullptr;
    std::vector<float> *py = nullptr;
    std::vector<float> *pz = nullptr;
    std::vector<float> *eta = nullptr;
    std::vector<float> *e_loss = nullptr;
    std::vector<float> *total_x0 = nullptr;
    std::vector<float> *total_l0 = nullptr;
    std::vector<float> *number_of_hits = nullptr;

    double mmuon = 0.1056583745; // Muon mass in GeV/c^2
    
    // Set branch addresses
    tree->SetBranchAddress("event_id",&event_id);
    tree->SetBranchAddress("particle_type",&particle_type);
    tree->SetBranchAddress("particle_id",&particle_id);
    tree->SetBranchAddress("vx",&vx);
    tree->SetBranchAddress("vy",&vy);
    tree->SetBranchAddress("vz",&vz); 
    tree->SetBranchAddress("vt",&vt);
    tree->SetBranchAddress("px",&px);
    tree->SetBranchAddress("py",&py);
    tree->SetBranchAddress("pz",&pz); 
    tree->SetBranchAddress("eta",&eta); 
    tree->SetBranchAddress("e_loss",&e_loss);
    tree->SetBranchAddress("total_x0",&total_x0);
    tree->SetBranchAddress("total_l0",&total_l0);
    tree->SetBranchAddress("number_of_hits",&number_of_hits);

    TH2F *hHitsVsEloss = new TH2F("hHitsVsEloss",";N_{hits};#DeltaE (GeV/#it{c}^{2});Entries",6,0.5,6.5,400,0,20);
    TH2F *hX0VsEloss = new TH2F("hX0VsEloss",";total X_0;#DeltaE (GeV/#it{c}^{2});Entries",400,0,30,400,0,20);
    TH1F *hL0tot = new TH1F("hL0tot",";total l_0;Entries",300,0,30);
    TH1F *hL0tothits = new TH1F("hL0tothits",";total l_0;Entries",300,0,30);
    TH1F *hElossHits = new TH1F("hElossHits",";#Delta E (GeV/#it{c}^{2});Entries",4*35,0.5,3);
    TH1F *hElossHitsNeg = new TH1F("hElossHitsNeg",";#Delta E (GeV/#it{c}^{2});Entries",4*35,0.5,3);
    TH1F *hElossHitsPos = new TH1F("hElossHitsPos",";#Delta E (GeV/#it{c}^{2});Entries",4*35,0.5,3);
    TH2F *hElossHitsVsP = new TH2F("hElossHitsVsP",";#Delta E (GeV/#it{c}^{2});p (GeV/#it{c});Entries",400,0,20,10,0,20);
    TH1F *hElossHitsGauss = new TH1F("hElossHitsGauss",";#DeltaE (GeV/#it{c}^{2});Entries",400,0,20);
    TH1F *hElossHitsDelta = new TH1F("hElossHitsDelta",";#DeltaE (GeV/#it{c}^{2});Entries",400,0,20);

    // Loop over entries in the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        for (int j = 0; j < particle_id->size(); ++j) {
            hX0VsEloss->Fill((*total_x0)[j], (*e_loss)[j]);
            //if((*total_x0)[j] > 17) continue;
            hL0tot->Fill((*total_l0)[j]);
            hHitsVsEloss->Fill((*number_of_hits)[j], (*e_loss)[j]);

            double energy = std::sqrt((*px)[j]*(*px)[j] + (*py)[j]*(*py)[j] + (*pz)[j]*(*pz)[j] + mmuon*mmuon)+ (*e_loss)[j];
            if ((*number_of_hits)[j] == 10 && (*eta)[j] < 3.2 && (*eta)[j] > 0 && energy > 10 && energy < 20) {
                hElossHits->Fill((*e_loss)[j]);
                hL0tothits->Fill((*total_l0)[j]);
                if ((*particle_type)[j] < 0) {
                    hElossHitsNeg->Fill((*e_loss)[j]);
                } else {
                    hElossHitsPos->Fill((*e_loss)[j]);
                }
                double p = std::sqrt((*px)[j]*(*px)[j] + (*py)[j]*(*py)[j] + (*pz)[j]*(*pz)[j]);
                hElossHitsVsP->Fill((*e_loss)[j], p);
            }
        }
    }
    hElossHits->GetMaximum();
    hElossHits->GetMaximumBin();
    for(int i=1;i<hElossHits->GetMaximumBin()*2;i++){
        if(i<=hElossHits->GetMaximumBin())
            hElossHitsGauss->SetBinContent(i,hElossHits->GetBinContent(i));
        else{
            int j = i-hElossHits->GetMaximumBin();
            hElossHitsGauss->SetBinContent(i,hElossHits->GetBinContent(hElossHits->GetMaximumBin()-j));
        }
        hElossHitsDelta->SetBinContent(i,hElossHits->GetBinContent(i)-hElossHitsGauss->GetBinContent(i));
        
    }
    int integralAll = 0;
    int integralDelta = 0;
    int integralGauss = 0;
    for(int i=1;i<hElossHits->GetNbinsX();i++){
        integralAll += hElossHits->GetBinContent(i);
        integralDelta += hElossHitsDelta->GetBinContent(i);;
        integralGauss += hElossHitsGauss->GetBinContent(i);;
    }
    std::cout<<"integralAll "<<integralAll<<std::endl;
    std::cout<<"integralGauss "<<integralGauss<<std::endl;
    std::cout<<"integralDelta "<<integralDelta<<std::endl;
    std::cout<<"Bias_ratio "<<(double)integralDelta/integralAll<<std::endl;
    TFile *f = new TFile("eloss.root","RECREATE");
    hHitsVsEloss->Write();
    hX0VsEloss->Write();
    hElossHitsGauss->Write();
    hElossHitsDelta->Write();
    hElossHitsVsP->Write();
    hElossHits->Write();
    hElossHitsNeg->Write();
    hElossHitsPos->Write();

    double sigma = GetFWHM(hElossHits)*0.42;



    TCanvas* cv = new TCanvas("cv", "", 1800,1600);
    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);

    //hElossHits->SetLineColor(kRed);
    std::cout<<"sigma+hElossHits->GetMean() "<<sigma+hElossHits->GetMean()<<std::endl;
    std::cout<<"sigma-hElossHits->GetMean() "<<hElossHits->GetMean()-sigma<<std::endl;
    TLine rline(hElossHits->GetMean()+sigma*3,0,hElossHits->GetMean()+sigma*3,800);
    TLine lline(hElossHits->GetMean()-sigma*3,0,hElossHits->GetMean()-sigma*3,800);
    rline.SetLineColor(kRed);
    lline.SetLineColor(kRed);

    hElossHits->Draw();
    rline.Draw();
    lline.Draw();

    std::string pngname =  std::string("eloss.png");
    cv->SaveAs(pngname.c_str());

    int out_of_range_count = countOutsideRange(hElossHits, hElossHits->GetMean()-sigma*3, hElossHits->GetMean()+sigma*3);

    std::cout << "Entries outside [" << hElossHits->GetMean()-sigma*3 << ", " << hElossHits->GetMean()+sigma*3.0 << "]: " 
              << out_of_range_count/hElossHits->GetEntries() << std::endl;

    hL0tothits->SetLineColor(kRed);
    hL0tot->Draw();
    hL0tothits->Draw("same");
    std::cout<<hL0tothits->GetMean()<<std::endl;
    pngname =  std::string("l0tot.png");
    cv->SaveAs(pngname.c_str());

    f->Close();
    return 0;
}
