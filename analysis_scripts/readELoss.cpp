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
    gStyle->SetOptStat(0);
    // Open the ROOT file
    TFile *file = TFile::Open("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_jpsi/particles_simulation.root", "READ");
    // Get the TTree
    TTree *tree = (TTree *)file->Get("particles");
    // Variables to hold branch data
    UInt_t event_id;
    std::vector<unsigned long> *particle_id = nullptr;
    std::vector<float> *vx = nullptr;
    std::vector<float> *vy = nullptr;
    std::vector<float> *vz = nullptr;
    std::vector<float> *vt = nullptr;
    std::vector<float> *px = nullptr;
    std::vector<float> *py = nullptr;
    std::vector<float> *pz = nullptr;
    std::vector<float> *e_loss = nullptr;
    std::vector<float> *total_x0 = nullptr;
    std::vector<float> *total_l0 = nullptr;
    std::vector<float> *number_of_hits = nullptr;
    
    // Set branch addresses
    tree->SetBranchAddress("event_id",&event_id);
    tree->SetBranchAddress("particle_id",&particle_id);
    tree->SetBranchAddress("vx",&vx);
    tree->SetBranchAddress("vy",&vy);
    tree->SetBranchAddress("vz",&vz); 
    tree->SetBranchAddress("vt",&vt);
    tree->SetBranchAddress("px",&px);
    tree->SetBranchAddress("py",&py);
    tree->SetBranchAddress("pz",&pz); 
    tree->SetBranchAddress("e_loss",&e_loss);
    tree->SetBranchAddress("total_x0",&total_x0);
    tree->SetBranchAddress("total_l0",&total_l0);
    tree->SetBranchAddress("number_of_hits",&number_of_hits);

    TH2F *hHitsVsEloss = new TH2F("hHitsVsEloss",";N_{hits};#DeltaE (GeV/#it{c}^{2});Entries",11,0.5,11.5,400,0,20);
    TH2F *hX0VsEloss = new TH2F("hX0VsEloss",";total X_0;#DeltaE (GeV/#it{c}^{2});Entries",400,0,30,400,0,20);
    TH1F *hL0tot = new TH1F("hL0tot",";total l_0;Entries",300,0,30);
    TH1F *hL0tot11hits = new TH1F("hL0tot11hits",";total l_0;Entries",300,0,30);
    TH1F *hEloss11Hits = new TH1F("hEloss11Hits",";#Delta E (GeV/#it{c}^{2});Entries",400,0,20);
    TH1F *hEloss11HitsGauss = new TH1F("hEloss11HitsGauss",";#DeltaE (GeV/#it{c}^{2});Entries",400,0,20);
    TH1F *hEloss11HitsDelta = new TH1F("hEloss11HitsDelta",";#DeltaE (GeV/#it{c}^{2});Entries",400,0,20);

    // Loop over entries in the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        for (int j = 0; j < particle_id->size(); ++j) {
            hX0VsEloss->Fill((*total_x0)[j], (*e_loss)[j]);
            hL0tot->Fill((*total_l0)[j]);
            hHitsVsEloss->Fill((*number_of_hits)[j], (*e_loss)[j]);
            if ((*number_of_hits)[j] == 11) {
                hEloss11Hits->Fill((*e_loss)[j]);
                hL0tot11hits->Fill((*total_l0)[j]);
            }
        }
    }
    hEloss11Hits->GetMaximum();
    hEloss11Hits->GetMaximumBin();
    for(int i=1;i<hEloss11Hits->GetMaximumBin()*2;i++){
        if(i<=hEloss11Hits->GetMaximumBin())
            hEloss11HitsGauss->SetBinContent(i,hEloss11Hits->GetBinContent(i));
        else{
            int j = i-hEloss11Hits->GetMaximumBin();
            hEloss11HitsGauss->SetBinContent(i,hEloss11Hits->GetBinContent(hEloss11Hits->GetMaximumBin()-j));
        }
        hEloss11HitsDelta->SetBinContent(i,hEloss11Hits->GetBinContent(i)-hEloss11HitsGauss->GetBinContent(i));
    }
    int integralAll = 0;
    int integralDelta = 0;
    int integralGauss = 0;
    for(int i=1;i<hEloss11Hits->GetNbinsX();i++){
        integralAll += hEloss11Hits->GetBinContent(i);
        integralDelta += hEloss11HitsDelta->GetBinContent(i);;
        integralGauss += hEloss11HitsGauss->GetBinContent(i);;
    }
    std::cout<<"integralAll "<<integralAll<<std::endl;
    std::cout<<"integralGauss "<<integralGauss<<std::endl;
    std::cout<<"integralDelta "<<integralDelta<<std::endl;
    std::cout<<"Bias_ratio "<<(double)integralDelta/integralAll<<std::endl;
    TFile *f = new TFile("eloss.root","RECREATE");
    hHitsVsEloss->Write();
    hX0VsEloss->Write();
    hEloss11HitsGauss->Write();
    hEloss11HitsDelta->Write();
    
    hEloss11Hits->Write();

    double sigma = GetFWHM(hEloss11Hits)*0.42;



    TCanvas* cv = new TCanvas("cv", "", 1800,1600);
    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);

    //hEloss11Hits->SetLineColor(kRed);
    std::cout<<"sigma+hEloss11Hits->GetMean() "<<sigma+hEloss11Hits->GetMean()<<std::endl;
    std::cout<<"sigma-hEloss11Hits->GetMean() "<<hEloss11Hits->GetMean()-sigma<<std::endl;
    TLine rline(hEloss11Hits->GetMean()+sigma*3,0,hEloss11Hits->GetMean()+sigma*3,800);
    TLine lline(hEloss11Hits->GetMean()-sigma*3,0,hEloss11Hits->GetMean()-sigma*3,800);
    rline.SetLineColor(kRed);
    lline.SetLineColor(kRed);

    //hEloss11HitsGauss->SetLineColor(kBlue);
    //hEloss11HitsDelta->SetLineColor(kGreen+2);
    /*
    leg->AddEntry(he_mass_Ref,Form("VT+MS RMS = %0.f #pm %0.f MeV/#it{c}^{2}, #sigma = %0.f #pm %0.f MeV/#it{c}^{2}",he_mass_Ref->GetRMS()*1000.,he_mass_Ref->GetRMSError()*1000., sigma_Ref, err_Ref), "l");
    leg->AddEntry(he_mass_VT,Form("VT RMS = %0.f #pm %0.f MeV/#it{c}^{2}, #sigma = %0.f #pm %0.f MeV/#it{c}^{2}",he_mass_VT->GetRMS()*1000.,he_mass_VT->GetRMSError()*1000., sigma_VT, err_VT), "l");
    leg->AddEntry(he_mass_MS,Form("MS RMS = %0.f #pm %0.f MeV/#it{c}^{2}, #sigma = %0.f #pm %0.f MeV/#it{c}^{2}",he_mass_MS->GetRMS()*1000.,he_mass_MS->GetRMSError()*1000., sigma_MS, err_MS), "l");
    leg->Draw();
    */
    hEloss11Hits->Draw();
    rline.Draw();
    lline.Draw();

    std::string pngname =  std::string("eloss.png");
    cv->SaveAs(pngname.c_str());

    int out_of_range_count = countOutsideRange(hEloss11Hits, hEloss11Hits->GetMean()-sigma*3, hEloss11Hits->GetMean()+sigma*3);

    std::cout << "Entries outside [" << hEloss11Hits->GetMean()-sigma*3 << ", " << hEloss11Hits->GetMean()+sigma*3.0 << "]: " 
              << out_of_range_count/hEloss11Hits->GetEntries() << std::endl;

    hL0tot11hits->SetLineColor(kRed);
    hL0tot->Draw();
    hL0tot11hits->Draw("same");
    std::cout<<hL0tot11hits->GetMean()<<std::endl;
    pngname =  std::string("l0tot.png");
    cv->SaveAs(pngname.c_str());

    f->Close();
    return 0;
}
