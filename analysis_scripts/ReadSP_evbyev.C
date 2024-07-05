void LoadStyle();

void ReadSP_evbyev(char *fname= "pippo.out"){

//LoadStyle();
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
		 
FILE *ff;
ff= fopen(fname,"r");
if(!ff){
  printf("No file found\n");
  exit(1);
}

//----------------------------------------------------------------------------------
// Create histos
//----------------------------------------------------------------------------------

Double_t y[5]={25.,50.,75.,100.,125};
Double_t r[5];
Double_t xx = 15.;
Double_t zz = 6.5;
for(int i=0;i<5;i++) r[i] = TMath::Sqrt(y[i]*y[i]);
//for(int i=0;i<5;i++) r[i] = TMath::Sqrt(y[i]*y[i]+xx*xx);
TLine *l[5];
for(int i=0;i<5;i++){
  l[i] = new TLine(-zz,r[i],zz,r[i]);
  l[i]->SetLineColor(2);
  l[i]->SetLineWidth(3);
}

const int nzgrid=4;
TLine *lzgrid[nzgrid];
Double_t zgrid[nzgrid]={-6.,-2.,2.,6.};
for(int i=0;i<nzgrid;i++){
  lzgrid[i] = new TLine(zgrid[i],-2,zgrid[i],2);
  lzgrid[i]->SetLineColor(kOrange);
  lzgrid[i]->SetLineWidth(3);
}

const int nn=1;
Double_t phi[nn];
for(int i=0;i<nn;i++){
  phi[i] = -TMath::Pi()+i*(TMath::Pi()*2/nn);
}
TLine *lphi[nn];
for(int i=0;i<nn;i++){
  lphi[i] = new TLine(phi[i],-15.,phi[i],15);
  lphi[i]->SetLineColor(2);
  lphi[i]->SetLineWidth(3);
}

//----------------------------------------------------------------------------------
// Read input file
//----------------------------------------------------------------------------------
char tag[150];

char str[200];

Int_t nsp= 0, nseed= 0;
Int_t countEv=0;
    
for(Int_t i=0;i<20000000;i++){ 
  if(i%1000==0)printf("%d lines read\n",i);  
  if(!fgets(str,200,ff)) {
     printf("END OF FILE\n");
     break;
  }

  Double_t spx, spy, spz;
  Double_t seedx, seedy, seedz;
  Double_t rr, Phi, z;
  Int_t numseed;

  TString tstr(str);
  TCanvas *c0;
  TH2D *hSP;
  TH2D *hSeed;
  TH2D *hSPPhiz;
  TH2D *hSeedPhiz;
  TH2D *hSPxy;
  TH2D *hSeedxy;
  TH2D *hSPyz;
  TH2D *hSeedyz;
  TString newev = kFALSE;
  
  if(tstr.Contains("NA60+_===================NewEvent=============================")) {
     c0 = new TCanvas("c0","c0",20,20,1500,1000);
     c0->Divide(2,2);
     hSP = new TH2D("hSP","SP",1000,-10,10,750,0.,150);
     hSP->GetYaxis()->SetTitle("R (mm)");
     hSP->GetXaxis()->SetTitle("z (mm)");
     hSP->SetMarkerSize(1.6);
     hSP->SetMarkerStyle(20);
     hSP->SetMarkerColor(4);

     hSeed = new TH2D("hSeed","Seed",1000,-10,10,750,0.,150);
     hSeed->GetYaxis()->SetTitle("R (mm)");
     hSeed->GetXaxis()->SetTitle("z (mm)");
     hSeed->SetMarkerStyle(20);
     hSeed->SetMarkerColor(8);

     hSPPhiz= new TH2D("hSPPhiz","SPPhiz",700,-7,7,1000,-2,2);
     hSPPhiz->GetYaxis()->SetTitle("Phi");
     hSPPhiz->GetXaxis()->SetTitle("z (mm)");
     hSPPhiz->SetMarkerColor(4);
     hSPPhiz->SetMarkerSize(1.6);
     hSPPhiz->SetMarkerStyle(20);

     hSeedPhiz= new TH2D("hSeedPhiz","SeedPhiz",700,-7,7,1000,-2,2);
     hSeedPhiz->GetYaxis()->SetTitle("Phi");
     hSeedPhiz->GetXaxis()->SetTitle("z (mm)");
     hSeedPhiz->SetMarkerStyle(20);
     hSeedPhiz->SetMarkerColor(8);

     hSPxy = new TH2D("hSPxy","SPxy",1500,-15,15,750,0.,150);
     hSPxy->GetYaxis()->SetTitle("y (mm)");
     hSPxy->GetXaxis()->SetTitle("x (mm)");
     hSPxy->SetMarkerSize(1.6);
     hSPxy->SetMarkerStyle(20);
     hSPxy->SetMarkerColor(4);

     hSeedxy = new TH2D("hSeedxy","Seedxy",1500,-15,15,750,0.,150);
     hSeedxy->GetYaxis()->SetTitle("y (mm)");
     hSeedxy->GetXaxis()->SetTitle("x (mm)");
     hSeedxy->SetMarkerStyle(20);
     hSeedxy->SetMarkerColor(8);

     hSPyz = new TH2D("hSPxy","SPxy",700,-7,7,750,0.,150);
     hSPyz->GetYaxis()->SetTitle("y (mm)");
     hSPyz->GetXaxis()->SetTitle("z (mm)");
     hSPyz->SetMarkerSize(1.6);
     hSPyz->SetMarkerStyle(20);
     hSPyz->SetMarkerColor(4);

     hSeedyz = new TH2D("hSeedxy","Seedxy",700,-7,7,750,0.,150);
     hSeedyz->GetYaxis()->SetTitle("y (mm)");
     hSeedyz->GetXaxis()->SetTitle("z (mm)");
     hSeedyz->SetMarkerStyle(20);
     hSeedyz->SetMarkerColor(8);
     
     printf("\n==== New event %d ===\n", countEv);
     countEv++; 
     numseed=0;
  }

       
  if(tstr.Contains("NA60+_position= ")) {
     nsp++;
     sscanf(str,"%s %lf %lf %lf",tag,&spx,&spy,&spz);
     printf("sp= %lf %lf %lf\n",spx,spy,spz);    
     hSP->Fill(spz,TMath::Sqrt(spx*spx+spy*spy));
     hSPPhiz->Fill(spz,TMath::ATan(spy/spx));
     hSPxy->Fill(spx,spy);
     hSPyz->Fill(spz,spy);
     
  }  
  if(tstr.Contains("NA60+_EstimateTrackParamsFromSeed_SP(x,y,z)=")) {
     nseed++;
     numseed++;
     sscanf(str,"%s %lf %lf %lf",tag,&seedx,&seedz,&seedy);    // credo di avere un - perche' z' = -y
     printf("seed= %lf %lf %lf\n",seedx,seedy,-seedz);	 
     hSeed->Fill(-seedz,TMath::Sqrt(seedx*seedx+seedy*seedy));
     hSeedPhiz->Fill(-seedz,TMath::ATan(seedy/seedx));
     hSeedxy->Fill(seedx,seedy);
     hSeedyz->Fill(-seedz,seedy);
     
  }  

//----------------------------------------------------------------------------------
// Draw histos
//----------------------------------------------------------------------------------
 
  if(tstr.Contains("finished event")) { 
     c0->cd(1);
     hSP->Draw();
     for(int i=0;i<5;i++) l[i]->Draw();
     hSP->Draw("same");
     hSeed->Draw("same");
     char event[30];
     sprintf(event,"Event n. %d",countEv-1);
     TLatex *lat = new TLatex(-9,140,event);
     lat->SetTextColor(kOrange);
     lat->Draw("same");
     c0->cd(2);
     hSPPhiz->Draw();
     hSeedPhiz->Draw("same");
     for(int i=0;i<nzgrid;i++) lzgrid[i]->Draw();
     c0->cd(3);
     hSPxy->Draw();
     hSeedxy->Draw("same");
     c0->cd(4);
     hSPyz->Draw();
     hSeedyz->Draw("same");
     char lnumseed[30];
     sprintf(lnumseed,"n. seed = %3.0f",numseed/3.);
     TLatex *lat2 = new TLatex(-6,130,lnumseed);
     lat2->SetTextColor(kBlue);
     lat2->SetTextSize(0.05);
     lat2->Draw("same");
     
     gPad->WaitPrimitive();     
   }
}

printf("\n Total number of SP = %d\n",nsp);
printf(" Total number of seeds = %f\n",nseed/3.);

}      
