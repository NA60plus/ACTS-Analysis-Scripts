#include <TFile.h>
#include <TTree.h>
#include <iostream>

int readMeas() {
    // Open the ROOT file
    TFile *file = TFile::Open("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_rubenxprino/measurements.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file!" << std::endl;
        return 1;
    }


    // Variables to hold branch data
    Int_t event_nr;
    Int_t volume_id;
    Int_t layer_id;
    Int_t surface_id;

    float rec_loc0;
    float rec_loc1;
    
    // Set branch addresses

    // Get the TTree
    TTree *tree = nullptr;
    TH2D *hHistXYLayer2_s1 = new TH2D("hHistXYLayer2_s1",";x (mm);y (mm)",155,-155,155,155,-155,155);
    std::vector<TH2D*> hHistXYLayers;
    for(int iPlane=0;iPlane<11;iPlane++){
        if(iPlane<5)
            hHistXYLayers.push_back(new TH2D(Form("hHistXYLayer%d",iPlane),";x (mm);y (mm)",155,-155,155,155,-155,155));
        else
            hHistXYLayers.push_back(new TH2D(Form("hHistXYLayer%d",iPlane),";x (mm);y (mm)",200,-3000,3000,200,-3000,3000));


        tree = (TTree *)file->Get(Form("vol1_lay%d",(iPlane+1)*2));
        tree->SetBranchAddress("event_nr",&event_nr);
        tree->SetBranchAddress("volume_id",&volume_id);
        tree->SetBranchAddress("layer_id",&layer_id);
        tree->SetBranchAddress("surface_id",&surface_id);

        tree->SetBranchAddress("rec_loc0",&rec_loc0);
        tree->SetBranchAddress("rec_loc1",&rec_loc1);
        if (!tree) {
            std::cerr << "Error: Cannot find TTree in file!" << std::endl;
            file->Close();
            return 1;
        }
        // Loop over entries in the tree
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            tree->GetEntry(i);
            hHistXYLayers[iPlane]->Fill(rec_loc0, rec_loc1);
            if(layer_id == 2 &&surface_id == 1){
                    hHistXYLayer2_s1->Fill(rec_loc0, rec_loc1);
                
            }
            //std::cout<<tx<<" "<<ty<<" "<<tz<<" "<<tt<<" "<<tpx<<" "<<tpy<<" "<<tpz<<" "<<te<<" "<<deltapx<<" "<<deltapy<<" "<<deltapz<<" "<<deltae<<" "<<index<<" "<<volume_id<<" "<<boundary_id<<" "<<layer_id<<" "<<approach_id<<" "<<surface_id<<std::endl;
        }
    }


    TFile *f = new TFile("outputMeas.root","RECREATE");
    for(int i=0;i<11;i++){
        hHistXYLayers[i]->Write();
    }
    hHistXYLayer2_s1->Write();
    f->Close();
    return 0;
}
