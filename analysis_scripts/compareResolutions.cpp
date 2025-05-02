#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <vector>
#include <string>

int compareResolutions() {
    std::vector<std::string> suffix_list = {"Carbon","Carbon_NoCross"};
    std::vector<int> color_list = {kRed, kBlue};
    std::vector<int> marker_list = {20, 21};
    std::vector<std::string> hist_list = {
        "sigma_d0_vs_p",
        "sigma_d0_vs_pt",
        "sigma_z0_vs_p",
        "sigma_z0_vs_pt",
        "sigma_poverp_vs_p",
        "sigma_ptoverpt_vs_pt"
    };

    TCanvas* cv = new TCanvas("cv", "cv", 1920, 1080);
    for (const auto& hist_name : hist_list) {
        bool first = true;
        for (size_t i = 0; i < suffix_list.size(); ++i) {
            std::string filename = "results/histogram" + suffix_list[i] + ".root";
            TFile* file = TFile::Open(filename.c_str());
            if (!file || file->IsZombie()) {
                std::cerr << "Error opening file: " << filename << std::endl;
                continue;
            }
    
            TH1D* hist = dynamic_cast<TH1D*>(file->Get(hist_name.c_str()));
            if (!hist) {
                std::cerr << "Histogram " << hist_name << " not found in " << filename << std::endl;
                file->Close();
                continue;
            }
    
            std::cout << "Entries in " << hist_name << " (" << filename << "): " << hist->GetEntries() << std::endl;
    
            // Print bin contents to confirm non-empty data
            for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
                std::cout << "Bin " << bin << ": " << hist->GetBinContent(bin) << std::endl;
            }
    
            hist->SetMarkerColor(color_list[i]);
            hist->SetLineColor(color_list[i]);
            hist->SetMarkerStyle(marker_list[i]);
            hist->SetMarkerSize(1.5);
    
            // Ensure axes are set properly
            hist->GetXaxis()->SetRangeUser(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
            hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.2);  // Scale Y-axis for visibility
    
            hist->Draw(first ? "E1" : "E1 SAME");
            first = false;
    
            file->Close();
        }
    
        // Force update of the canvas
        gPad->Modified();
        gPad->Update();
    
        std::string save_name = "results/histogram" + hist_name + ".png";
        cv->SaveAs(save_name.c_str());
    }
    delete cv;
    return 0;
}
