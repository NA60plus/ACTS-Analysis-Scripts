void LoadStyle();

void ReadActs_Sprotated(TString EvByEv = kFALSE, Int_t nev = -1, char *fname = "testK0S40GeV.out"){

//LoadStyle();
		 
FILE *ff;
ff= fopen(fname,"r");
if(!ff){
  printf("No file found\n");
  exit(1);
}

//----------------------------------------------------------------------------------
// Create histos
//----------------------------------------------------------------------------------
TH2D *hzvsphiM= new TH2D("hzvsphiM","z vs phi Middle SP",100,-TMath::Pi(),TMath::Pi(),5,12.5,137.5);
TH2D *hrvsphiM= new TH2D("hrvsphiM","r vs phi Middle SP",100,-TMath::Pi(),TMath::Pi(),100,0.,7.5);
TH2D *hzvsphiT= new TH2D("hzvsphiT","z vs phi Top SP",100,-TMath::Pi(),TMath::Pi(),5,12.5,137.5);
TH2D *hrvsphiT= new TH2D("hrvsphiT","r vs phi Top SP",100,-TMath::Pi(),TMath::Pi(),100,0.,7.5);
TH2D *hzvsphiB= new TH2D("hzvsphiB","z vs phi Bottom SP",100,-TMath::Pi(),TMath::Pi(),5,12.5,137.5);
TH2D *hrvsphiB= new TH2D("hrvsphiB","r vs phi Bottom SP",100,-TMath::Pi(),TMath::Pi(),100,0.,7.5);

gStyle->SetOptStat(1111);
TH1D *hrM = new TH1D("hrM","hrM",400,0.,400.);
hrM->GetXaxis()->SetTitle("r Middle");
TH1D *hBottomdeltaR = new TH1D("hBottomdeltaR","hBottomdeltaR",300,0.,300.);
hBottomdeltaR->GetXaxis()->SetTitle("doublet cut: R |Bottom-Middle|");
TH1D *hTopdeltaR = new TH1D("hTopdeltaR","hTopdeltaR",300,0.,300.);
hTopdeltaR->GetXaxis()->SetTitle("doublet cut: R |Top-Middle|");

TH1D *hBottomzOrigin = new TH1D("hBottomzOrigin","hBottomzOrigin",3000,-15.,15.);
hBottomzOrigin->GetXaxis()->SetTitle("doublet cut Bottom-Middle: zOrigin= (zM * deltaR - rM * deltaZ)/deltaR");
TH1D *hBottomcotTheta = new TH1D("hBottomcotTheta","hBottomcotTheta",200,-1.,1.);
hBottomcotTheta->GetXaxis()->SetTitle("doublet cut Bottom-Middle: deltaZdeltaR = cottheta");
TH1D *hBottomdeltaZ = new TH1D("hBottomdeltaZ","hBottomdeltaZ",200,-20.,20.);
hBottomdeltaZ->GetXaxis()->SetTitle("doublet cut Bottom-Middle: deltaZ= |zM - ZBottom or Top|");

TH1D *hTopzOrigin = new TH1D("hTopzOrigin","hTopzOrigin",3000,-15.,15.);
hTopzOrigin->GetXaxis()->SetTitle("doublet cut Top-Middle: zOrigin= (zM * deltaR - rM * deltaZ)/deltaR");
TH1D *hTopcotTheta = new TH1D("hTopcotTheta","hTopcotTheta",200,-1.,1.);
hTopcotTheta->GetXaxis()->SetTitle("doublet cut Top-Middle: deltaZdeltaR = cottheta");
TH1D *hTopdeltaZ = new TH1D("hTopdeltaZ","hTopdeltaZ",200,-20.,20.);
hTopdeltaZ->GetXaxis()->SetTitle("doublet cut Top-Middle: deltaZ= |zM - ZTop or Top|");

TH1D *hdeltaCotTheta2PtMin = new TH1D("hdeltaCotTheta2PtMin","hdeltaCotTheta2PtMin-error2-scatteringInRegion",1000,-0.06,0.005);
hdeltaCotTheta2PtMin->GetXaxis()->SetTitle("deltaCotTheta2PtMin-error2-scatteringInRegion");
TH1D *hHelixCut = new TH1D("hHelixCut","hHelixCut",10000,0.,1.e8);
hHelixCut->GetXaxis()->SetTitle("hHelixCut");
TH1D *hdeltaCotTheta2PtEstimate = new TH1D("hdeltaCotTheta2PtEstimate","hdeltaCotTheta2PtEstimate-error2-scatteringInRegion",1000,-0.005,0.005);
hdeltaCotTheta2PtEstimate->GetXaxis()->SetTitle("hdeltaCotTheta2PtEstimate-error2-scatteringInRegion");
TH1D *hIm = new TH1D("hIm","hImpactParameter",1000,0.,10.);
hIm->GetXaxis()->SetTitle("Impact parameter");

//TH1D *hInvHelixDiameter = new TH1D("hInvHelixDiameter","hInvHelixDiameter in filterSeeds_2SpFixed",1000,-1,1);
//hinvHelixDiameter->GetXaxis()->SetTitle("invHelixDiameter");
TH1D *hinvHelixDiameterVec = new TH1D("hinvHelixDiameterVec","invHelixDiameterVec[CompatibleSeeds] in filterSeeds_2SpFixed",400,-0.0002,0.0002);
hinvHelixDiameterVec->GetXaxis()->SetTitle("invHelixDiameterVec[compatible-default]");
TH1D *hdeltaR_in2SpFixed = new TH1D("hdeltaR_in2SpFixed","deltaR in filterSeeds_2SpFixed",600,-300,300);
hdeltaR_in2SpFixed->GetXaxis()->SetTitle("deltaR_in2SpFixed");

TH1D *hCotThetaAvg2 = new TH1D("hCotThetaAvg2","hCotThetaAvg2 in FilterCandidates",350,-0.02,0.05);
hCotThetaAvg2->GetXaxis()->SetTitle("CotThetaAvg2");


TH2D *hSPInSeedsZR= new TH2D("hSPInSeedsZR","hSPInSeedsZR",300,-150,150,300,0,TMath::Sqrt(381.175*381.175+150*150));
hSPInSeedsZR->GetXaxis()->SetTitle("z (mm)");
hSPInSeedsZR->GetYaxis()->SetTitle("r (mm)");

Int_t nlayers= 5;
TH2D *hSPInSeedsXZ[nlayers];
char name[100];
for(int i=0;i<nlayers;i++){
  sprintf(name,"hSPInSeedsXZ%d",i);
  hSPInSeedsXZ[i] = new TH2D(name,name,300,-150,150,300,-150,150);
  hSPInSeedsXZ[i]->GetXaxis()->SetTitle("x (mm)");
  hSPInSeedsXZ[i]->GetYaxis()->SetTitle("z (mm)");
}
TH1D *hSP[nlayers];
for(int i=0;i<nlayers;i++){
  sprintf(name,"hSP%d",i);
  hSP[i] = new TH1D(name,name,250,0.,250);
  hSP[i]->GetXaxis()->SetTitle("R [mm] (in x-z plane)");
  hSP[i]->GetYaxis()->SetTitle("dN/(2Pi*R*DeltaR)/Nev [mm2]");
}

TH2D *hNSeedsvsNSP= new TH2D("hNSeedsvsNSP","hNSeedsvsNSP",500,0.,500,500,0.,500.);
hNSeedsvsNSP->GetXaxis()->SetTitle("N SP");
hNSeedsvsNSP->GetYaxis()->SetTitle("N seeds");

TH1D *hbestSeedQuality = new TH1D("hbestSeedQuality","filterSeeds_1SpFixed_bestSeedQuality",100,0,1000);
hbestSeedQuality->GetXaxis()->SetTitle("bestSeedQuality");

TH1D *hcompatibleSeed = new TH1D("hcompatibleSeed","hcompatibleSeed",10,0.,10);
hcompatibleSeed->GetXaxis()->SetTitle("compatibleSeeds");

TH1D *hParamFromSeed_Loc0 = new TH1D("hParamFromSeed_Loc0","hParamFromSeed_Loc0",200,-150,150);
hParamFromSeed_Loc0->GetXaxis()->SetTitle("ParamFromSeed_Loc0");
TH1D *hParamFromSeed_Loc1 = new TH1D("hParamFromSeed_Loc1","hParamFromSeed_Loc1",200,-150,150);
hParamFromSeed_Loc1->GetXaxis()->SetTitle("ParamFromSeed_Loc1");
TH1D *hParamFromSeed_Phi = new TH1D("hParamFromSeed_Phi","hParamFromSeed_Phi",200,-TMath::Pi(),TMath::Pi());
hParamFromSeed_Phi->GetXaxis()->SetTitle("ParamFromSeed_Phi");
TH1D *hParamFromSeed_Theta = new TH1D("hParamFromSeed_Theta","hParamFromSeed_Theta",200,-1,1);
hParamFromSeed_Theta->GetXaxis()->SetTitle("ParamFromSeed_Theta");
TH1D *hParamFromSeed_QOverP = new TH1D("hParamFromSeed_QOverP","hParamFromSeed_QOverP",200,-10,10);
hParamFromSeed_QOverP->GetXaxis()->SetTitle("ParamFromSeed_QOverP");

TH1D *hChi2CKF = new TH1D("hChi2CKF","hChi2CKF",500,0.,500);
hChi2CKF->SetTitle("Chi2 CKF");

TH1D *hSeedPerformances = new TH1D("hSeedPerformances","Seed Performances per event",5,0,5);
TString namelab[5]= {"nTotSeeds","nTotMatchedSeeds","nTotParticles","nTotMatchedParticles","nTotDuplicatedParticles"};
for(int k=0;k<5;k++) hSeedPerformances->GetXaxis()->SetBinLabel(k+1,namelab[k]);
hSeedPerformances->GetXaxis()->SetLabelSize(0.07);

TH1D *hSeedPerformances2 = new TH1D("hSeedPerformances2","Seed Performances",5,0,5);
TString namelab2[5]= {"Efficiency","Fake rate","Seed purity","DuplicationRate","AvNumberDuplicatedSeeds"};
for(int k=0;k<5;k++) hSeedPerformances2->GetXaxis()->SetBinLabel(k+1,namelab2[k]);
hSeedPerformances2->GetXaxis()->SetLabelSize(0.07);

TH1D *hEtaSeeds = new TH1D("hEtaSeeds","hEtaSeeds per event",80,-1,1);
hEtaSeeds->GetXaxis()->SetTitle("eta seed");
hEtaSeeds->GetYaxis()->SetTitle("eta seed per event");

TH1D *hPtSeeds = new TH1D("hPtSeeds","hPtSeeds per event",250,0,5);
hPtSeeds->GetXaxis()->SetTitle("Pt seed");
hPtSeeds->GetYaxis()->SetTitle("Pt seed per event");

TH1D *hd0Seeds = new TH1D("hd0Seeds","hd0Seeds per event",200,-0,0.4);
hd0Seeds->GetXaxis()->SetTitle("d0 seed");
hd0Seeds->GetYaxis()->SetTitle("d0 seed per event");

TH1D *hCurvatureSeeds = new TH1D("hCurvatureSeeds","hCurvatureSeeds per event",100,-0.0005,0.0005);
hCurvatureSeeds->GetXaxis()->SetTitle("curvature seed");
hCurvatureSeeds->GetYaxis()->SetTitle("curvature seed per event");


//----------------------------------------------------------------------------------
// Read input file
//----------------------------------------------------------------------------------
char tag[150];
Double_t r, phi,z;
Double_t invH;
Double_t spx, spy, spz;
Double_t spx_seed, spy_seed, spz_seed;
Int_t nseed, nSP;
Double_t rM, rMinMiddle, rMaxMiddle;
Double_t BottomdeltaR, BottomdeltaRMinSP;
Double_t TopdeltaR, TopdeltaRMaxSP;
Double_t zOriginBottom, collRegionMinBottom, collRegionMaxBottom;
Double_t deltaZBottom, deltaZMinBottom, deltaZMaxBottom;
Double_t deltaZdeltaRBottom, cotThetaMinBottom, cotThetaMaxBottom;
Double_t zOriginTop, collRegionMinTop, collRegionMaxTop;
Double_t deltaZTop, deltaZMinTop, deltaZMaxTop;
Double_t deltaZdeltaRTop, cotThetaMinTop, cotThetaMaxTop;

Double_t Helix, minHelix;
Double_t deltaTheta, error22, p2scatterSigma;
Double_t Im, impactMax;
Double_t invHelixDiameterVec, deltaHelix, invHelix;
Double_t deltaR_in2SpFixed, Rmin;
Double_t deltaCotTheta2, error2, scatteringInRegion2;
Double_t compatibleSeed;
Double_t ParamLoc0, ParamLoc1, ParamPhi, ParamTheta, ParamQOverP;
Double_t CotThetaAvg2;
Double_t Chi2CKF, Chi2Cut;
Double_t bestSeedQuality;

Double_t Seed_cotThetaAvg;
Double_t Seed_curvature;
Double_t Seed_d0;
Double_t Seed_pT;

Double_t nTotalSeeds, nTotalMatchedSeeds, nTotalParticles, nTotalMatchedParticles, nTotalDuplicatedParticles;
Double_t Efficiency, Fake, Purity, Duplication, AvNumberDuplicatedSeeds;

Int_t countEv=0;

Int_t pippo=0, pippo1=0, pippo2=0, pippo3=0, pippo4=0, pippo5=0, pippo6=0, pippo7=0;
Int_t  pippo4_1=0, pippo4_2=0, pippo4_3=0, pippo4_4=0, pippo4_5=0, pippo4_6=0, pippo4_7=0, pippo4_8=0;
Int_t pippo8=0, pippo9=0, pippo10=0, pippo11=0, pippo12=0, pippo13=0, pippo14=0, pippo15=0;

char str[200];
    
for(Int_t i=0;i<100000000;i++){ 
  if(i%1000000==0)printf("%d lines read\n",i);  
  if(!fgets(str,200,ff)) {
     printf("END OF FILE\n");
     break;
  }
       
  TString tstr(str);
  if(tstr.Contains("NA60+_===================NewEvent=============================")) {
     countEv++; 
  }
  
  if((EvByEv==kTRUE && countEv == (nev+1)) || (EvByEv==kFALSE)) {
  
  if(tstr.Contains("NA60+_position=")) {
     sscanf(str,"%s %lf %lf %lf",tag,&spx,&spy,&spz); 
     Double_t rr = TMath::Sqrt(spx*spx+spz*spz);   
     if(spy==71.175) hSP[0]->Fill(rr);
     if(spy==151.175) hSP[1]->Fill(rr);
     if(spy==201.175) hSP[2]->Fill(rr);
     if(spy==251.175) hSP[3]->Fill(rr);
     if(spy==381.175) hSP[4]->Fill(rr);
  }

  if(tstr.Contains("NA60+_EstimateTrackParamsFromSeed_SP(x,y,z)=")) {
     sscanf(str,"%s %lf %lf %lf",tag,&spx_seed,&spy_seed,&spz_seed);   
     if(spy<250) 
     printf("SP x= %f, y= %f, z= %f\n",spx,spy,spz);
     hSPInSeedsZR->Fill(spz_seed,TMath::Sqrt(spx_seed*spx_seed+spy_seed*spy_seed));
     if(spy==71.175) hSPInSeedsXZ[0]->Fill(spx_seed,spz_seed);
     if(spy==151.175) hSPInSeedsXZ[1]->Fill(spx_seed,spz_seed);
     if(spy==201.175) hSPInSeedsXZ[2]->Fill(spx_seed,spz_seed);
     if(spy==251.175) hSPInSeedsXZ[3]->Fill(spx_seed,spz_seed);
     if(spy==381.175) hSPInSeedsXZ[4]->Fill(spx_seed,spz_seed);
  }  
//   if(tstr.Contains("SeedingAlgorithm_NSeed,NSP: ")) {
//      sscanf(str,"%s %d %d",tag,&nseed,&nSP);    
// //     printf("NSeeds= %d NSP= %d\n",nseed,nSP);
//      hNSeedsvsNSP->Fill(nSP,nseed);
//   }
  if(tstr.Contains("NA60+_SeedFinder_middleSP_rM,rMinMiddle,rMaxMiddle: ")) {
     sscanf(str,"%s %lf %lf %lf",tag,&rM,&rMinMiddle,&rMaxMiddle);    
//     printf("rM= %f ,rMinMiddle= %f ,rMaxMiddle= %f\n",rM,rMinMiddle,rMaxMiddle);
     hrM->Fill(rM);  
   }  
  if(tstr.Contains("NA60+_SeedFinder_getCompatibleDoublets_Bottom_DeltaR_deltRMin= ")) {
     sscanf(str,"%s %lf %lf",tag,&BottomdeltaR,&BottomdeltaRMinSP);    
     hBottomdeltaR->Fill(BottomdeltaR);  
  }
  if(tstr.Contains("NA60+_SeedFinder_getCompatibleDoublets_Top_DeltaR_deltRMax= ")) {
     sscanf(str,"%s %lf %lf",tag,&TopdeltaR,&TopdeltaRMaxSP);    
     hTopdeltaR->Fill(TopdeltaR);  
  }
//   if(tstr.Contains("SeedFinder_getCompatibleDoublets_Top_DeltaR_deltaRMin,deltRMax= ")) {
//      sscanf(str,"%s %lf %lf %lf",tag,&TopdeltaR,&TopdeltaRMinSP,&TopdeltaRMaxSP);    
//      hTopdeltaR->Fill(TopdeltaR);  
//   }
  if(tstr.Contains("NA60+_SeedFinder_getCompatibleDoublets_BOTTOM_zOriginTimesDeltaR/DeltaR_collRegionMin,collRegionMax= ")) {
     sscanf(str,"%s %lf %lf %lf",tag,&zOriginBottom,&collRegionMinBottom,&collRegionMaxBottom);    
//     printf("zOrigin= %f collRegionMin= %f collRegionMax= %f\n",zOrigin,collRegionMin,collRegionMax);
     hBottomzOrigin->Fill(zOriginBottom);  
  }
   if(tstr.Contains("NA60+_SeedFinder_getCompatibleDoublets_BOTTOM_deltaZ/deltaR,cotThetaMin,cotThetaMax: ")) {
     sscanf(str,"%s %lf %lf %lf",tag,&deltaZdeltaRBottom,&cotThetaMinBottom,&cotThetaMaxBottom);    
//     printf("deltaZdeltaR= %f, cotThetaMin= %f ,cotThetaMax= %f\n",deltaZdeltaR,cotThetaMin,cotThetaMax);
     hBottomcotTheta->Fill(deltaZdeltaRBottom);  
   }
   if(tstr.Contains("NA60+_SeedFinder_getCompatibleDoublets_BOTTOM_deltaZ,deltaZMin,deltaZMax: ")) {
     sscanf(str,"%s %lf %lf %lf",tag,&deltaZBottom,&deltaZMinBottom,&deltaZMaxBottom);   
//     printf("deltaZ= %f deltaZMin= %f ,deltaZMax= %f\n",deltaZ,deltaZMin,deltaZMax); 
     hBottomdeltaZ->Fill(deltaZBottom);  
   }  
  if(tstr.Contains("NA60+_SeedFinder_getCompatibleDoublets_TOP_zOriginTimesDeltaR/DeltaR_collRegionMin,collRegionMax= ")) {
     sscanf(str,"%s %lf %lf %lf",tag,&zOriginTop,&collRegionMinTop,&collRegionMaxTop);    
//     printf("zOrigin= %f collRegionMin= %f collRegionMax= %f\n",zOrigin,collRegionMin,collRegionMax);
     hTopzOrigin->Fill(zOriginTop);  
  }
   if(tstr.Contains("NA60+_SeedFinder_getCompatibleDoublets_TOP_deltaZ/deltaR,cotThetaMin,cotThetaMax: ")) {
     sscanf(str,"%s %lf %lf %lf",tag,&deltaZdeltaRTop,&cotThetaMinTop,&cotThetaMaxTop);    
//     printf("deltaZdeltaR= %f, cotThetaMin= %f ,cotThetaMax= %f\n",deltaZdeltaR,cotThetaMin,cotThetaMax);
     hTopcotTheta->Fill(deltaZdeltaRTop);  
   }
   if(tstr.Contains("NA60+_SeedFinder_getCompatibleDoublets_TOP_deltaZ,deltaZMin,deltaZMax: ")) {
     sscanf(str,"%s %lf %lf %lf",tag,&deltaZTop,&deltaZMinTop,&deltaZMaxTop);   
//     printf("deltaZ= %f deltaZMin= %f ,deltaZMax= %f\n",deltaZ,deltaZMin,deltaZMax); 
     hTopdeltaZ->Fill(deltaZTop);  
   }  
   
   
   if(tstr.Contains("NA60+_SeedFinder_filterCandidates_deltaCotTheta2,error2,scatteringInRegion2: ")){
     sscanf(str,"%s %lf %lf %lf",tag,&deltaCotTheta2,&error2,&scatteringInRegion2); 
     hdeltaCotTheta2PtMin->Fill(deltaCotTheta2-(error2+scatteringInRegion2));  
   }   
   if(tstr.Contains("NA60+_SeedFinder_filterCandidates_minHelixDiameter2: ")) {
     sscanf(str,"%s %lf %lf",tag,&Helix,&minHelix);   
     //printf("Helix= %f %f\n",Helix,minHelix);
     hHelixCut->Fill(Helix);  
   }
   if(tstr.Contains("NA60+_SeedFinder_filterCandidates_deltaTheta_scattering,error2,p2scatterSigma: ")) {
     sscanf(str,"%s %lf %lf %lf",tag,&deltaTheta,&error22,&p2scatterSigma);   
     hdeltaCotTheta2PtEstimate->Fill(deltaTheta-error22-p2scatterSigma);  
   }
   if(tstr.Contains("NA60+_SeedFinder_filterCandidates_Im_impactMax: ")) {
     sscanf(str,"%s %lf %lf",tag,&Im,&impactMax);   
     hIm->Fill(Im);  
   }
   
//=========================================================================   
   if(tstr.Contains("NA60+_SeedFilter_2SP_invHelixDiameterVec[compatibleTopSPIndex]_invHelixDiameter[topSPIndex]_m_cfg.deltaInvHelixDiameter: ")) {
     sscanf(str,"%s %lf %lf %lf",tag,&invHelixDiameterVec,&invHelix,&deltaHelix);   
     hinvHelixDiameterVec->Fill(invHelixDiameterVec-invHelix);  
   }
   if(tstr.Contains("NA60+_SeedFilter_2SP_deltaR_in2SpFixed,min: ")) {
     sscanf(str,"%s %lf %lf",tag,&deltaR_in2SpFixed,&Rmin);   
     hdeltaR_in2SpFixed->Fill(deltaR_in2SpFixed);  
   }
   if(tstr.Contains("NA60+_compatibleSeedR.size(): ")) {
     sscanf(str,"%s %lf",tag,&compatibleSeed);   
     hcompatibleSeed->Fill(compatibleSeed);  
   }
   if(tstr.Contains("NA60+_SeedFinder_cotThetaAvg2= ")) {
     sscanf(str,"%s %lf",tag,&CotThetaAvg2);   
     hCotThetaAvg2->Fill(CotThetaAvg2);  
   }
   if(tstr.Contains("NA60+_SeedFilter_filterSeeds_1SpFixed_bestSeedQuality= ")) {
     sscanf(str,"%s %lf",tag,&bestSeedQuality);   
     hbestSeedQuality->Fill(bestSeedQuality);  
   }

//=========================================================================      
   if(tstr.Contains("NA60+_SeedParameter_cotThetaAvg= ")) {
     sscanf(str,"%s %lf",tag,&Seed_cotThetaAvg); 
     Double_t theta = TMath::ATan(1./TMath::Abs(Seed_cotThetaAvg));  
     if(Seed_cotThetaAvg<0) hEtaSeeds->Fill(TMath::Log(TMath::Tan(theta/2.)));  // to allow negative eta
     else hEtaSeeds->Fill(-TMath::Log(TMath::Tan(theta/2.)));  
     //printf("cotThetaAvg2= %f Eta= %f\n",Seed_cotThetaAvg,-TMath::Log(TMath::Tan(theta/2.)));
   }
   if(tstr.Contains("NA60+_SeedParameter_curvature= ")) {
     sscanf(str,"%s %lf",tag,&Seed_curvature);   
     hCurvatureSeeds->Fill(Seed_curvature);  
   }
   if(tstr.Contains("NA60+_SeedParameter_d0= ")) {
     sscanf(str,"%s %lf",tag,&Seed_d0);   
     hd0Seeds->Fill(Seed_d0);  
   }
   if(tstr.Contains("NA60+_SeedParameter_pT= ")) {
     sscanf(str,"%s %lf",tag,&Seed_pT);   
     hPtSeeds->Fill(Seed_pT);  
     //printf("pT= %f\n",Seed_pT);
   }
   
   
   
  if(tstr.Contains("NA60+_EstimateTrackParamsFromSeed_eBoundLoc0= ")) {
    sscanf(str,"%s %lf",tag,&ParamLoc0);   
    printf("ParamLoc0= %f\n",ParamLoc0);
    hParamFromSeed_Loc0->Fill(ParamLoc0);  
  }
  if(tstr.Contains("NA60+_EstimateTrackParamsFromSeed_eBoundLoc1= ")) {
    sscanf(str,"%s %lf",tag,&ParamLoc1);   
    hParamFromSeed_Loc1->Fill(ParamLoc1);  
  }
  if(tstr.Contains("NA60+_EstimateTrackParamsFromSeed_eBoundPhi= ")) {
    sscanf(str,"%s %lf",tag,&ParamPhi);   
    hParamFromSeed_Phi->Fill(ParamPhi);  
  }
  if(tstr.Contains("NA60+_EstimateTrackParamsFromSeed_eBoundTheta= ")) {
    sscanf(str,"%s %lf",tag,&ParamTheta);   
    hParamFromSeed_Theta->Fill(ParamTheta);  
  }
  if(tstr.Contains("NA60+_EstimateTrackParamsFromSeed_eBoundQOverP= ")) {
    sscanf(str,"%s %lf",tag,&ParamQOverP);   
    hParamFromSeed_QOverP->Fill(ParamQOverP);  
  }
    
  if(tstr.Contains("NA60+_MeasurementSelector.hpp_Chi2_Chi2Cut: ")) {
    sscanf(str,"%s %lf %lf",tag,&Chi2CKF,&Chi2Cut);   
    hChi2CKF->Fill(Chi2CKF);  
  }
     
  if(tstr.Contains("NA60+_Summary_nTotalSeeds=")) {
    sscanf(str,"%s %lf",tag,&nTotalSeeds);   
    printf("nTotalSeeds= %f\n",nTotalSeeds);
    hSeedPerformances->SetBinContent(1,nTotalSeeds);  
  }
  if(tstr.Contains("NA60+_Summary_nTotalMatchedSeeds=")) {
    sscanf(str,"%s %lf",tag,&nTotalMatchedSeeds);   
    hSeedPerformances->SetBinContent(2,nTotalMatchedSeeds);  
  }
  if(tstr.Contains("NA60+_Summary_nTotalParticles=")) {
    sscanf(str,"%s %lf",tag,&nTotalParticles);   
    hSeedPerformances->SetBinContent(3,nTotalParticles);  
  }
  if(tstr.Contains("NA60+_Summary_nTotalMatchedParticles=")) {
    sscanf(str,"%s %lf",tag,&nTotalMatchedParticles);   
    hSeedPerformances->SetBinContent(4,nTotalMatchedParticles);  
  }
  if(tstr.Contains("NA60+_Summary_nTotalDuplicatedParticles=")) {
    sscanf(str,"%s %lf",tag,&nTotalDuplicatedParticles);   
    hSeedPerformances->SetBinContent(5,nTotalDuplicatedParticles);  
  }

 if(tstr.Contains("NA60+_Summary_Eff=")) {
    sscanf(str,"%s %lf",tag,&Efficiency);   
    hSeedPerformances2->SetBinContent(1,Efficiency);  
  }
  if(tstr.Contains("NA60+_Summary_Fakerate=")) {
    sscanf(str,"%s %lf",tag,&Fake);   
    hSeedPerformances2->SetBinContent(2,Fake);  
  }
  if(tstr.Contains("NA60+_Summary_Purity=")) {
    sscanf(str,"%s %lf",tag,&Purity);   
    hSeedPerformances2->SetBinContent(3,Purity);  
  }
  if(tstr.Contains("NA60+_Summary_Duplication=")) {
    sscanf(str,"%s %lf",tag,&Duplication);   
    hSeedPerformances2->SetBinContent(4,Duplication);  
  }
  if(tstr.Contains("NA60+_Summary_nDuplicatedSeeds=")) {
    sscanf(str,"%s %lf",tag,&AvNumberDuplicatedSeeds);   
    hSeedPerformances2->SetBinContent(5,AvNumberDuplicatedSeeds);  
  }
//=========================================================  
  
  if(tstr.Contains("NA60+_SeedParameter_curvature=")) {
    pippo++;
  }
  if(tstr.Contains("NA60+_SeedFilter_filterSeeds_2SpFixed")) {
    pippo1++; 
  }
  if(tstr.Contains("NA60+_SeedFilter_topSpVec.size()=")) {
    pippo2++; 
  }
  if(tstr.Contains("NA60+_seed_filter_afterdeltaR")) {
    pippo3++; 
  }
  if(tstr.Contains("NA60+_updated weight= m_cfg.compatSeedWeight")) {
    pippo4++; 
  }
  if(tstr.Contains("NA60+_test1")) {
    pippo4_6++; 
  }
  if(tstr.Contains("NA60+_test2")) {
    pippo4_7++; 
  }
  if(tstr.Contains("NA60+_test3")) {
    pippo4_8++; 
  }
  if(tstr.Contains("NA60+_SeedFilter keep seed with good weight")) {
    pippo5++; 
  }
  if(tstr.Contains("NA60+_SeedFilter_filterSeeds_1SpFixed maxSeeds")) {
    pippo6++; 
  }
  if(tstr.Contains("NA60+_SeedFilter_1SpFixed inside loop ")) {
    pippo4_1++; 
  }
  if(tstr.Contains("NA60+_SeedFilter_filterSeeds_1SpFixed assign quality")) {
    pippo4_2++; 
  }
  if(tstr.Contains("NA60+_SeedFilter_filterSeeds_1SpFixedbestSeedQuality")) {
    pippo7++; 
  }
//   if(tstr.Contains("compatibleSeedR.size():")) {
//     pippo8++; 
//   }
//   if(tstr.Contains("before check")) {
//     pippo9++; 
//   }
//   if(tstr.Contains("updated weight+= m_cfg.seedWeightIncrement")) {
//     pippo10++; 
//   }
//   if(tstr.Contains("before 1")) {
//     pippo11++; 
//   }
//   if(tstr.Contains("before 2")) {
//     pippo12++; 
//   }
//   if(tstr.Contains("before 3")) {
//     pippo13++; 
//   }
//   if(tstr.Contains("SeedFinder_confirmation_weighinconfirmation=")) {  //lose events here
//     pippo14++; 
//   }
//   if(tstr.Contains("SeedFilter check weight ")) {
//     pippo15++; 
//   }
  
  
}
}

printf("\n Number ov events= %d\n",countEv);

printf("1--- N entries in SeedFinder_filtercandidates end of Top loop = %d\n",pippo);
printf("2--- N entries in SeedFinder_filtercandidates end of Bottom loop = %d (ie number of good BottomSP)\n",pippo1);
printf("3--- N entries at the beginning of filterSeeds_2SpFixed = %d  (ie number of good BottomSP)\n",pippo2);
printf("4--- N entries in filterSeeds_2SpFixed fill filterhisto = %d\n",pippo3);
printf("5--- N entries in filterSeeds_2SpFixed end loop othertop = %d\n",pippo4);
// printf("--- test1 = %d\n",pippo8);
// printf("--- test2 = %d\n",pippo9);
// printf("--- test3 = %d\n",pippo10);
// printf("--- test4 = %d\n",pippo11);
// printf("--- test5 = %d\n",pippo12);
// printf("--- test6 = %d\n",pippo13);
// printf("--- test7 = %d\n",pippo14);
// printf("--- test8 = %d\n",pippo15);
printf("6--- N entries in filterSeeds_2SpFixed end loop top = %d \n",pippo5);
printf("7--- N entries in filterSeeds_1SpFixed beginning = %d\n",pippo6);
printf("8--- N entries in filterSeeds_1SpFixed inside loop %d\n",pippo4_1);
printf("9--- N entries in filterSeeds_2SpFixed seedconfirmation %d\n",pippo4_2);
printf("10--- N entries in filterSeeds_1SpFixed end = %d\n",pippo7);

//----------------------------------------------------------------------------------
// Plot histos
//----------------------------------------------------------------------------------
/*
TCanvas *c1 = new TCanvas("c1","z vs phi SP");
c1->Divide(2,2);
c1->cd(1);
hzvsphiM->Draw("colz");
c1->cd(3);
hzvsphiT->Draw("colz");
c1->cd(4);
hzvsphiB->Draw("colz");

TCanvas *c2 = new TCanvas("c2","r vs phi SP");
c2->Divide(2,2);
c2->cd(1);
hrvsphiM->Draw("colz");
c2->cd(3);
hrvsphiT->Draw("colz");
c2->cd(4);
hrvsphiB->Draw("colz");
*/
TCanvas *c4 = new TCanvas("c4","SP in seeds",1400,1000);
c4->Divide(3,2);
c4->cd(1);
hSPInSeedsZR->Draw("colz");
c4->cd(2);
hSPInSeedsXZ[0]->Draw("colz");
c4->cd(3);
hSPInSeedsXZ[1]->Draw("colz");
c4->cd(4);
hSPInSeedsXZ[2]->Draw("colz");
c4->cd(5);
hSPInSeedsXZ[3]->Draw("colz");
c4->cd(6);
hSPInSeedsXZ[4]->Draw("colz");

//TCanvas *c5 = new TCanvas("c5","NSeeds vs N SP");
//hNSeedsvsNSP->Draw("colz");

//=========================================================================
TCanvas *c6 = new TCanvas("c6","SeedFinder",1400,1000);
c6->Divide(3,3);
c6->cd(1);
hrM->Draw();
TLine *lrMinMiddle = new TLine(rMinMiddle,0.,rMinMiddle,hrM->GetMaximum());
TLine *lrMaxMiddle = new TLine(rMaxMiddle,0.,rMaxMiddle,hrM->GetMaximum());
lrMinMiddle->SetLineColor(2);
lrMaxMiddle->SetLineColor(2);
lrMinMiddle->SetLineWidth(2);
lrMaxMiddle->SetLineWidth(2);
lrMinMiddle->Draw("same");
lrMaxMiddle->Draw("same");

c6->cd(2);
hTopdeltaR->Draw();
TLine *lTopdeltaRMaxSP = new TLine(TopdeltaRMaxSP,0.,TopdeltaRMaxSP,1000);
lTopdeltaRMaxSP->SetLineColor(2);
lTopdeltaRMaxSP->SetLineWidth(2);
lTopdeltaRMaxSP->Draw("same");

c6->cd(3);
hTopzOrigin->Draw();
TLine *lTopcollRegionMin = new TLine(collRegionMinTop,0.,collRegionMinTop,hTopzOrigin->GetMaximum());
TLine *lTopcollRegionMax = new TLine(collRegionMaxTop,0.,collRegionMaxTop,hTopzOrigin->GetMaximum());
lTopcollRegionMin->SetLineColor(2);
lTopcollRegionMax->SetLineColor(2);
lTopcollRegionMin->SetLineWidth(2);
lTopcollRegionMax->SetLineWidth(2);
lTopcollRegionMin->Draw("same");
lTopcollRegionMax->Draw("same");

c6->cd(4);
hTopcotTheta->Draw();
TLine *lTopcotThetaMin = new TLine(cotThetaMinTop,0.,cotThetaMinTop,hTopcotTheta->GetMaximum());
TLine *lTopcotThetaMax = new TLine(cotThetaMaxTop,0.,cotThetaMaxTop,hTopcotTheta->GetMaximum());
lTopcotThetaMin->SetLineColor(2);
lTopcotThetaMax->SetLineColor(2);
lTopcotThetaMin->SetLineWidth(2);
lTopcotThetaMax->SetLineWidth(2);
lTopcotThetaMin->Draw("same");
lTopcotThetaMax->Draw("same");

c6->cd(5);
hTopdeltaZ->Draw();
TLine *lTopdeltaZMin = new TLine(deltaZMinTop,0.,deltaZMinTop,hTopdeltaZ->GetMaximum());
TLine *lTopdeltaZMax = new TLine(deltaZMaxTop,0.,deltaZMaxTop,hTopdeltaZ->GetMaximum());
lTopdeltaZMin->SetLineColor(2);
lTopdeltaZMax->SetLineColor(2);
lTopdeltaZMin->SetLineWidth(2);
lTopdeltaZMax->SetLineWidth(2);
lTopdeltaZMin->Draw("same");
lTopdeltaZMax->Draw("same");

c6->cd(6);
hBottomdeltaR->Draw();
TLine *lBottomdeltaRMinSP = new TLine(BottomdeltaRMinSP,0.,BottomdeltaRMinSP,hBottomdeltaR->GetMaximum());
lBottomdeltaRMinSP->SetLineColor(2);
lBottomdeltaRMinSP->SetLineWidth(2);
lBottomdeltaRMinSP->Draw("same");

c6->cd(7);
hBottomzOrigin->Draw();
TLine *lBottomcollRegionMin = new TLine(collRegionMinBottom,0.,collRegionMinBottom,hBottomzOrigin->GetMaximum());
TLine *lBottomcollRegionMax = new TLine(collRegionMaxBottom,0.,collRegionMaxBottom,hBottomzOrigin->GetMaximum());
lBottomcollRegionMin->SetLineColor(2);
lBottomcollRegionMax->SetLineColor(2);
lBottomcollRegionMin->SetLineWidth(2);
lBottomcollRegionMax->SetLineWidth(2);
lBottomcollRegionMin->Draw("same");
lBottomcollRegionMax->Draw("same");

c6->cd(8);
hBottomcotTheta->Draw();
TLine *lBottomcotThetaMin = new TLine(cotThetaMinBottom,0.,cotThetaMinBottom,hBottomcotTheta->GetMaximum());
TLine *lBottomcotThetaMax = new TLine(cotThetaMaxBottom,0.,cotThetaMaxBottom,hBottomcotTheta->GetMaximum());
lBottomcotThetaMin->SetLineColor(2);
lBottomcotThetaMax->SetLineColor(2);
lBottomcotThetaMin->SetLineWidth(2);
lBottomcotThetaMax->SetLineWidth(2);
lBottomcotThetaMin->Draw("same");
lBottomcotThetaMax->Draw("same");

c6->cd(9);
hBottomdeltaZ->Draw();
TLine *lBottomdeltaZMin = new TLine(deltaZMinBottom,0.,deltaZMinBottom,hBottomdeltaZ->GetMaximum());
TLine *lBottomdeltaZMax = new TLine(deltaZMaxBottom,0.,deltaZMaxBottom,hBottomdeltaZ->GetMaximum());
lBottomdeltaZMin->SetLineColor(2);
lBottomdeltaZMax->SetLineColor(2);
lBottomdeltaZMin->SetLineWidth(2);
lBottomdeltaZMax->SetLineWidth(2);
lBottomdeltaZMin->Draw("same");
lBottomdeltaZMax->Draw("same");


//=========================================================================
TCanvas *c7 = new TCanvas("c7","SeedFinder_FilterCandidates",1400,1000);
c7->Divide(3,2);
c7->cd(1);
hCotThetaAvg2->Draw();
TLine *lCotThetaAvg2 = new TLine(0,0.,0,hCotThetaAvg2->GetMaximum());
lCotThetaAvg2->SetLineColor(2);
lCotThetaAvg2->SetLineWidth(2);
lCotThetaAvg2->Draw("same");

c7->cd(2);
hdeltaCotTheta2PtMin->Draw();
TLine *lscatteringInRegion2 = new TLine(0,0.,0,hdeltaCotTheta2PtMin->GetMaximum());
lscatteringInRegion2->SetLineColor(2);
lscatteringInRegion2->SetLineWidth(2);
lscatteringInRegion2->Draw("same");

c7->cd(3);
gPad->SetLogy(1);
hHelixCut->Draw();
TLine *lHelixMin = new TLine(minHelix,0.,minHelix,hHelixCut->GetMaximum());
lHelixMin->SetLineColor(2);
lHelixMin->SetLineWidth(2);
lHelixMin->Draw("same");

c7->cd(4);
hdeltaCotTheta2PtEstimate->Draw();
TLine *lp2scatterSigma = new TLine(0,0.,0,hdeltaCotTheta2PtEstimate->GetMaximum());
lp2scatterSigma->SetLineColor(2);
lp2scatterSigma->SetLineWidth(2);
lp2scatterSigma->Draw("same");

c7->cd(5);
gPad->SetLogy(1);
hIm->Draw();
TLine *limpactMax = new TLine(impactMax,0.,impactMax,hIm->GetMaximum());
limpactMax->SetLineColor(2);
limpactMax->SetLineWidth(2);
limpactMax->Draw("same");

//=========================================================================
TCanvas *c8 = new TCanvas("c8","SeedFilter",1400,1000);
c8->Divide(2,2);


//c8->cd(1);
//hInvHelixDiameter->Draw();

c8->cd(1);
hinvHelixDiameterVec->Draw();
TLine *llowdeltaHelix = new TLine(-deltaHelix,0.,-deltaHelix,hinvHelixDiameterVec->GetMaximum());
TLine *lhighdeltaHelix = new TLine(deltaHelix,0.,deltaHelix,hinvHelixDiameterVec->GetMaximum());
llowdeltaHelix->SetLineColor(2);
lhighdeltaHelix->SetLineColor(2);
llowdeltaHelix->SetLineWidth(2);
lhighdeltaHelix->SetLineWidth(2);
llowdeltaHelix->Draw("same");
lhighdeltaHelix->Draw("same");

c8->cd(2);
hdeltaR_in2SpFixed->Draw();
TLine *lRmin = new TLine(-Rmin,0.,-Rmin,hdeltaR_in2SpFixed->GetMaximum());
lRmin->SetLineColor(2);
lRmin->SetLineWidth(2);
lRmin->Draw("same");
TLine *lRmax = new TLine(Rmin,0.,Rmin,hdeltaR_in2SpFixed->GetMaximum());
lRmax->SetLineColor(2);
lRmax->SetLineWidth(2);
lRmax->Draw("same");

c8->cd(3);
gPad->SetLogy(1);
hbestSeedQuality->Draw();

//=========================================================================
TCanvas *c9 = new TCanvas("c9","CompatibleSeed");
c9->Divide(1,1);
hcompatibleSeed->Draw();

//=========================================================================
TCanvas *c10 = new TCanvas("c10","seed parameters",1400,1000);
c10->Divide(2,2);
c10->cd(1);
hEtaSeeds->Scale(1./countEv);
hEtaSeeds->Draw("hist");
c10->cd(2);
hd0Seeds->Scale(1./countEv);
hd0Seeds->Draw("hist");
c10->cd(3);
hPtSeeds->Scale(1./countEv);
hPtSeeds->Draw("hist");
c10->cd(4);
hCurvatureSeeds->Scale(1./countEv);
hCurvatureSeeds->Draw("hist");

//=========================================================================
TCanvas *c11 = new TCanvas("c11","EstimateTrackParamsFromSeed",1400,1000);
c11->Divide(3,2);
c11->cd(1);
hParamFromSeed_Loc0->Draw();
c11->cd(2);
hParamFromSeed_Loc1->Draw();
c11->cd(3);
hParamFromSeed_Phi->Draw();
c11->cd(4);
hParamFromSeed_Theta->Draw();
c11->cd(5);
hParamFromSeed_QOverP->Draw();

//=========================================================================
/*
TCanvas *c112 = new TCanvas("c12","CKF");
hChi2CKF->Draw();
TLine *lChi2Cut = new TLine(Chi2Cut,0.,Chi2Cut,hChi2CKF->GetMaximum());
lChi2Cut->SetLineColor(2);
lChi2Cut->SetLineWidth(2);
lChi2Cut->Draw("same");
*/
//=========================================================================

for(int i=0;i<nlayers;i++){
  for(int j=1;j<=hSP[i]->GetNbinsX();j++){
    //if (i==2) printf("j=%d i=%d content = %f countEv= %d\n",i, j, hSP[i]->GetBinContent(j),countEv);
    hSP[i]->SetBinContent(j,hSP[i]->GetBinContent(j)/countEv/(2.*TMath::Pi()*hSP[i]->GetBinCenter(j)*hSP[i]->GetBinWidth(j)));
  }
}
TCanvas *c13 = new TCanvas("c13","SP",1400,1000);
c13->Divide(3,2);
for(int i=0;i<nlayers;i++){
  c13->cd(i+1);
  hSP[i]->Draw("");
}
//=========================================================================
TCanvas *c14 = new TCanvas("c14","SeedPerformances",1400,1000);
c14->Divide(1,2);
c14->cd(1);
hSeedPerformances->Draw("hist");
hSeedPerformances->Scale(1./countEv);
hSeedPerformances->SetLineColor(kBlue);
hSeedPerformances->SetFillColor(kOrange);
c14->cd(2);
hSeedPerformances2->Draw();
hSeedPerformances2->SetLineColor(kRed);
hSeedPerformances2->SetFillColor(kMagenta);


c4->SaveAs("SPXZ.png");
c6->SaveAs("SeedFinder.png");
c7->SaveAs("SeedFilterCand.png");
c8->SaveAs("SeedFilter.png");
c9->SaveAs("CompatSeed.png");
c10->SaveAs("SeedParams.png");
c11->SaveAs("EstimateTrackParamsFromSeed.png");
c13->SaveAs("SP.png");
c14->SaveAs("SeedPerformances.png");

}      
