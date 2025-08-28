#include <TFile.h>
#include <TTree.h>
#include <iostream>

int readHits() {
    // Open the ROOT file
    TFile *file = TFile::Open("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/test_geometry_output/hits.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file!" << std::endl;
        return 1;
    }

    // Get the TTree
    TTree *tree = (TTree *)file->Get("hits");
    if (!tree) {
        std::cerr << "Error: Cannot find TTree in file!" << std::endl;
        file->Close();
        return 1;
    }

    // Variables to hold branch data
    UInt_t event_id;
    ULong64_t geometry_id;
    ULong64_t particle_id;
    float tx;
    float ty;
    float tz; 
    float tt;
    float tpx;
    float tpy;
    float tpz; 
    float te;
    float deltapx;
    float deltapy;
    float deltapz; 
    float deltae;
    Int_t index;
    UInt_t volume_id;
    UInt_t boundary_id;
    UInt_t layer_id;
    UInt_t approach_id;
    UInt_t sensitive_id;
    
    // Set branch addresses
    tree->SetBranchAddress("event_id",&event_id);
    tree->SetBranchAddress("geometry_id",&geometry_id);
    tree->SetBranchAddress("particle_id",&particle_id);
    tree->SetBranchAddress("tx",&tx);
    tree->SetBranchAddress("ty",&ty);
    tree->SetBranchAddress("tz",&tz); 
    tree->SetBranchAddress("tt",&tt);
    tree->SetBranchAddress("tpx",&tpx);
    tree->SetBranchAddress("tpy",&tpy);
    tree->SetBranchAddress("tpz",&tpz); 
    tree->SetBranchAddress("te",&te);
    tree->SetBranchAddress("deltapx",&deltapx);
    tree->SetBranchAddress("deltapy",&deltapy);
    tree->SetBranchAddress("deltapz",&deltapz); 
    tree->SetBranchAddress("deltae",&deltae);
    tree->SetBranchAddress("index",&index);
    tree->SetBranchAddress("volume_id",&volume_id);
    tree->SetBranchAddress("boundary_id",&boundary_id);
    tree->SetBranchAddress("layer_id",&layer_id);
    tree->SetBranchAddress("approach_id",&approach_id);
    tree->SetBranchAddress("sensitive_id",&sensitive_id);
    std::vector<TH2D*> hHistXYLayers;
    for(int i=0;i<11;i++){
        if(i<5)
        hHistXYLayers.push_back(new TH2D(Form("hHistXYLayer%d",i),";x (mm);y (mm)",155,-155,155,155,-155,155));
        else
            hHistXYLayers.push_back(new TH2D(Form("hHistXYLayer%d",i),";x (mm);y (mm)",200,-3000,3000,200,-3000,3000));
    }
    TH2D *hHistXYLayer2_s1 = new TH2D("hHistXYLayer2_s1",";x (mm);y (mm)",155,-155,155,155,-155,155);

    // Loop over entries in the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        int station = layer_id/2-1;
        hHistXYLayers[station]->Fill(tx,ty);
        if(layer_id == 2 &&sensitive_id == 1){
                hHistXYLayer2_s1->Fill(tx,ty);
            
        }
        //std::cout<<tx<<" "<<ty<<" "<<tz<<" "<<tt<<" "<<tpx<<" "<<tpy<<" "<<tpz<<" "<<te<<" "<<deltapx<<" "<<deltapy<<" "<<deltapz<<" "<<deltae<<" "<<index<<" "<<volume_id<<" "<<boundary_id<<" "<<layer_id<<" "<<approach_id<<" "<<sensitive_id<<std::endl;
    }


    TFile *f = new TFile("output.root","RECREATE");
    for(int i=0;i<11;i++){
        hHistXYLayers[i]->Write();
    }
    hHistXYLayer2_s1->Write();
    f->Close();
    return 0;
}
