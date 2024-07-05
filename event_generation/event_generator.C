
#include <cmath>
#include <iostream>
#include <string>
#include <bitset>
#include <TF1.h>
#include <TRandom.h>
#include <filesystem>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TParticle.h>
#include <TH3D.h>

#include "TDatabasePDG.h"
#include "AliAODRecoDecay.h"
#include "AliDecayer.h"
#include "AliDecayerEvtGen.h"

#include "strange_particles_data.h"

namespace fs = std::filesystem;

TH3D *h3Dpow;
TH1D *hD0pt;
TH1D *hD0y;
TF1 *fpt;
TF1 *fy;
TF1 *fdNdYPi, *fdNdYK, *fdNdYP, *fdNdPtPi, *fdNdPtK, *fdNdPtP;
double fNChPi, fNChK, fNChP;

void InitBgGenerationPart(double, double, double, double, double,
                          double, double, double, double, double,
                          double, double, double, double,
                          double, double,
                          double, double, double, double, double);
void GenBgEvent(double, double, double, int, bool, bool, bool, AliDecayerEvtGen *, FILE *, int, float);

const double kMassP = 0.938;
const double kMassK = 0.4937;
const double kMassPi = 0.1396;
const double kMassMu = 0.1057;
const double kMassE = 0.0005;

void event_generator(int nev,
                     int Eint = 40,
                     bool addBkg = true,
                     bool addD0 = false,
                     bool addSec = false,
                     double beamSigma = 0,
                     bool fullTargetSystem = false,
                     int onlyT = -1,
                     float periferal = 1)
{

  AliDecayerEvtGen *fDecayer = new AliDecayerEvtGen();
  fDecayer->Init(); // read the default decay table DECAY.DEC and particle table
  const char *privateDecayTable = "decaytables/USERTABD0.DEC";

  TFile *filPow = new TFile("pp0_frag-PtSpectra-Boost.root");
  if (addD0)
  {
    h3Dpow = (TH3D *)filPow->Get("hptyeta421");
    hD0pt = (TH1D *)h3Dpow->ProjectionX("hD0pt");
    hD0y = (TH1D *)h3Dpow->ProjectionY("hD0y");
  }
  else if (addSec)
  {
    privateDecayTable = "decaytables/USERTABK0.DEC";
  }
  fDecayer->SetDecayTablePath((char *)privateDecayTable);
  fDecayer->ReadDecayTable();

  fDecayer->ForceDecay();

  // (1) 40 GeV pi K PRC66 (2002)054902
  // (2) 40 GeV p PRC83 (2011) 014901
  if (Eint == 40)
    printf("--- Background parameters for E=40 GeV/nucleon ---\n");
  double y0BG = 2.22;  // gaussian y mean - 40 GeV
  double sigyBG = 1.2; // .. sigma
  double TBG = 0.17;   // inv.slope of thermal pt distribution

  double y0BGPi = 0.666;       // (1) mean of the two gaussians displaced symmetrically wrt y=0
  double y0BGKplus = 0.694;    // (1)
  double y0BGKminus = 0.569;   // (1)
  double y0BGP = 0.907;        // (2) checked approx
  double sigyBGPi = 0.872;     // (1) sigma of the two gaussians displaced symmetrically wrt y=0
  double sigyBGKplus = 0.725;  // (1)
  double sigyBGKminus = 0.635; // (1)
  double sigyBGP = 0.798;      // (2) checked approx
  double NBGPi = 74.;          //(1) normalization paramters of the sum of two gaussian
  double NBGKplus = 16.2;      // (1)
  double NBGKminus = 6.03;     // (1)
  double NBGP = 37.5;          // (2) checked approx
  double Piratio = 0.91;
  double TBGpi = 0.169; //(1) from a fit in 0.2<M_t-m<0.7 GeV
  double TBGK = 0.229;  // (1) average of K+ (232 MeV) and K- (226 MeV)
  double TBGP = 0.25;   // ?
  double yminBG = -20;  // min y to generate
  double ymaxBG = 20;   // max y to generate
  double ptminBG = 0.;  // min pt to generate
  double ptmaxBG = 100; // max pt to generate

  double dndyBGPi = 615.;
  double dndyBGK = 78.;
  double dndyBGP = 150.;

  if (Eint == 158)
  {
    // pions and Kaons from  NA49 nucl-ex/0205002
    printf("--- Background parameters for E=158 GeV/nucleon ---\n");
    y0BG = 2.9; // gaussian y mean - 160 GeV
    y0BGPi = 0.72;
    y0BGKplus = 0.839;
    y0BGKminus = 0.727;
    y0BGP = 39.8;
    sigyBGPi = 1.18;
    sigyBGKplus = 0.88;
    sigyBGKminus = 0.81;
    sigyBGP = 8.07;
    dndyBGPi = 1258.;
    dndyBGK = 155.;
    dndyBGP = 292.;
    NBGPi = 107.6;
    NBGKplus = 23.4;
    NBGKminus = 12.8;
    NBGP = 2.55e+06;
    Piratio = 0.97;
    TBGpi = 0.18; // inv.slope of thermal pt distribution
    TBGK = 0.23;  // inv.slope of thermal pt distribution
    TBGP = 0.31;  // inv.slope of thermal pt distribution
  }

  InitBgGenerationPart(NBGPi, NBGKplus, NBGKminus, NBGP, Piratio, y0BG, y0BGPi, y0BGKplus, y0BGKminus, y0BGP, sigyBGPi, sigyBGKplus, sigyBGKminus, sigyBGP, yminBG, ymaxBG, TBGpi, TBGK, TBGP, ptminBG, ptmaxBG);
  std::string options = "";
  if (!addBkg)
    options += "_noBkg";
  if (addD0)
    options += "_D0";
  else if (addSec)
    options += "_Sec";
  if (fullTargetSystem)
  {
    options += "_realTarget";
    options += "_beamSigma_" + std::to_string(beamSigma);
  }
  else if (onlyT != -1)
  {
    options += "_onlyTarget_" + std::to_string(onlyT);
    options += "_beamSigma_" + std::to_string(beamSigma);
  }
  if (periferal != 1)
    options += "_periferal_factor_" + std::to_string(periferal);
  std::string directoryName = "events_" + std::to_string(Eint) + "GeV" + options;

  if (!fs::exists(directoryName))
  {                                                     // Check if directory doesn't exist
    bool created = fs::create_directory(directoryName); // Create the directory
    if (created)
    {
      std::cout << "Directory created successfully: " << directoryName << std::endl;
    }
    else
    {
      std::cerr << "Failed to create directory: " << directoryName << std::endl;
      return; // Return error code
    }
  }
  else
  {
    std::cout << "Directory already exists: " << directoryName << std::endl;
  }
  printf("pion   multiplicity in %f<y<%f = %f\n", yminBG, ymaxBG, fNChPi);
  printf("kaon   multiplicity in %f<y<%f = %f\n", yminBG, ymaxBG, fNChK);
  printf("proton multiplicity in %f<y<%f = %f\n", yminBG, ymaxBG, fNChP);

  double vX = 0, vY = 0, vZ = 0;

  // 1.5 mm thick Pb disks spaced by 12
  int nTargets = 5;
  double targetThickness = 1.5;   // mm
  double distBetweenTargets = 12; // mm
  char csvname[80];

  gRandom->SetSeed(1801);
  std::string qa_file_name = directoryName+"/QA_plots.root";
  TFile *check = new TFile(qa_file_name.c_str(), "recreate");
  TH1D *hGenZ = new TH1D("hGenZ", ";z (mm);Entries", 1100, -100, 10);
  TH1D *hGenY = new TH1D("hGenY", ";y (mm);Entries", 100, -3, 3);
  TH1D *hGenX = new TH1D("hGenX", ";x (mm);Entries", 100, -3, 3);
  TH2D *hGenXY = new TH2D("hGenXY", ";x (mm);y (mm);Entries", 100, -3, 3, 100, -3, 3);
  FILE *fpcsv;
  for (int i = 0; i < nev; i++)
  {

    int nzeros = 9;
    int eventNumber = i;
    while (eventNumber != 0)
    {
      eventNumber /= 10;
      nzeros--;
    }
    if (i == 0)
      nzeros = 8;

    std::string eventFile = directoryName + "/event" + std::string(nzeros, '0') + std::to_string(static_cast<int>(i)) + "-particles.csv";
    std::cout << "File: " << eventFile << std::endl;
    sprintf(csvname, eventFile.c_str(), i);

    fpcsv = fopen(csvname, "w");
    printf("\nEvent no. %d\n", i);
    if (beamSigma > 0)
    {
      vX = gRandom->Gaus(0, beamSigma);
      while (TMath::Abs(vX) > 0.5)
        vX = gRandom->Gaus(0, beamSigma);

      vY = gRandom->Gaus(0, beamSigma);
      while (TMath::Abs(vY) > 0.5)
        vY = gRandom->Gaus(0, beamSigma);
    }
    if (fullTargetSystem)
    {
      int target = static_cast<int>(gRandom->Rndm() * nTargets);
      vZ = -(gRandom->Rndm() * targetThickness + target * (distBetweenTargets + targetThickness));
      std::cout << vZ << std::endl;
    }
    else if (onlyT != -1)
    {
      vZ = -(gRandom->Rndm() * targetThickness + (onlyT) * (distBetweenTargets + targetThickness));
    }
    hGenZ->Fill(vZ);
    hGenX->Fill(vX);
    hGenY->Fill(vY);
    hGenXY->Fill(vX, vY);
    GenBgEvent(vX, vY, vZ, i, addBkg, addD0, addSec, fDecayer, fpcsv, Eint, periferal);
    fclose(fpcsv);
  }
  filPow->Close();
  check->cd();
  hGenZ->Write();
  hGenX->Write();
  hGenY->Write();
  hGenXY->Write();
  check->Close();
  // fDecayer->Delete();
}

void InitBgGenerationPart(double NPi, double NKplus, double NKminus, double NP, double Piratio,
                          double y0, double y0Pi, double y0Kplus, double y0Kminus, double y0P,
                          double sigyPi, double sigyKplus, double sigyKminus, double sigyP,
                          double ymin, double ymax,
                          double Tpi, double TK, double TP, double ptmin, double ptmax)
{
  // initialize bg generation routines
  fNChPi = NPi * (1 + Piratio) * sigyPi * TMath::Sqrt(TMath::Pi() / 2.) * (TMath::Erf((ymax - y0 - y0Pi) / sqrt(2.) / sigyPi) + TMath::Erf((ymax - y0 + y0Pi) / sqrt(2.) / sigyPi) - TMath::Erf((ymin - y0 - y0Pi) / sqrt(2.) / sigyPi) - TMath::Erf((ymin - y0 + y0Pi) / sqrt(2.) / sigyPi));
  fNChK = NKplus * sigyKplus * TMath::Sqrt(TMath::Pi() / 2.) * (TMath::Erf((ymax - y0 - y0Kplus) / sqrt(2.) / sigyKplus) + TMath::Erf((ymax - y0 + y0Kplus) / sqrt(2.) / sigyKplus) - TMath::Erf((ymin - y0 - y0Kplus) / sqrt(2.) / sigyKplus) - TMath::Erf((ymin - y0 + y0Kplus) / sqrt(2.) / sigyKplus)) +
          NKminus * sigyKminus * TMath::Sqrt(TMath::Pi() / 2.) * (TMath::Erf((ymax - y0 - y0Kminus) / sqrt(2.) / sigyKminus) + TMath::Erf((ymax - y0 + y0Kminus) / sqrt(2.) / sigyKminus) - TMath::Erf((ymin - y0 - y0Kminus) / sqrt(2.) / sigyKminus) - TMath::Erf((ymin - y0 + y0Kminus) / sqrt(2.) / sigyKminus));
  fNChP = NP * sigyP * TMath::Sqrt(TMath::Pi() / 2.) * (TMath::Erf((ymax - y0 - y0P) / sqrt(2.) / sigyP) + TMath::Erf((ymax - y0 + y0P) / sqrt(2.) / sigyP) - TMath::Erf((ymin - y0 - y0P) / sqrt(2.) / sigyP) - TMath::Erf((ymin - y0 + y0P) / sqrt(2.) / sigyP));
  //
  fdNdYPi = new TF1("dndy", "exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )", ymin, ymax);
  fdNdYPi->SetParameters(y0, y0Pi, sigyPi);
  fdNdYK = new TF1("dndy", "exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )", ymin, ymax);
  fdNdYK->SetParameters(y0, y0Kplus, sigyKplus);
  fdNdYP = new TF1("dndy", "37.45*(exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2)) )", ymin, ymax);
  fdNdYP->SetParameters(y0, y0P, sigyP);

  fdNdPtPi = new TF1("dndptPi", "x*exp(-sqrt(x*x+1.949e-02)/[0])", ptmin, ptmax); // pion
  fdNdPtK = new TF1("dndptK", "x*exp(-sqrt(x*x+0.493*0.493)/[0])", ptmin, ptmax); // kaon
  fdNdPtP = new TF1("dndptP", "x*exp(-sqrt(x*x+0.938*0.938)/[0])", ptmin, ptmax); // proton
  fdNdPtPi->SetParameter(0, Tpi);
  fdNdPtK->SetParameter(0, TK);
  fdNdPtP->SetParameter(0, TP);
}

void GenBgEvent(double x, double y, double z, int event, bool addBkg, bool addD0, bool addSec, AliDecayerEvtGen *fDecayer, FILE *fpcsv, int Eint, float periferal)
{

  if (fNChPi < 0 && fNChK < 0 && fNChP < 0)
    return;
  // generate bg. events from simple thermalPt-gaussian Y parameterization
  if (!fdNdYPi || !fdNdYK || !fdNdYP || !fdNdPtPi || !fdNdPtK || !fdNdPtP)
  {
    printf("Background generation was not initialized\n");
    exit(1);
  }

  //               in pdg pr x  y  z t  px py pz m
  /* from BARCODE
  () const { return level(0); }
  /// Return the secondary vertex identifier.
  () const { return level(1); }
  /// Return the particle identifier.
  () const { return level(2); }
  /// Return the generation identifier.
  () const { return level(3); }
  /// Return the sub-particle identifier.
  subParticle() const { return level(4);

  vertexPrimary(12)-vertexSecondary(12)-particle(16)-generation(8)-subParticle(16)
  classActsFatras::Barcode : public Acts::MultiIndex<uint64_t, 12, 12, 16, 8, 16>
  */
  //
  int ntrTot = 0;
  // pions
  double ntrPi = gRandom->Poisson(fNChPi * periferal);
  printf("fNChPi=%f ntrPi=%f\n", fNChPi, ntrPi);
  printf("vX=%f vY=%f vZ=%f \n", x, y, z);
  // 2|0|14|0|0
  int particleNumber = 1;
  fprintf(fpcsv, "particle_id,particle_type,process,vx,vy,vz,vt,px,py,pz,m,q\n");
  if (addBkg)
  {
    for (int itr = 0; itr < ntrPi; itr++)
    {
      double yrap = fdNdYPi->GetRandom();
      double pt = fdNdPtPi->GetRandom();
      double phi = gRandom->Rndm() * TMath::Pi() * 2;
      float charge = gRandom->Rndm() > 0.52 ? 1 : -1;
      int ptypecsv = 211;
      double masscsv = 0.13957039;
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + kMassPi * kMassPi) * TMath::SinH(yrap)};

      std::string binaryString = std::bitset<12>(event + 1).to_string() + std::string(12, '0') + std::bitset<16>(particleNumber++).to_string() + std::string(24, '0');
      unsigned long long barcode = std::stoull(binaryString, nullptr, 2);
      fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, ptypecsv, 0, x, y, z, 0., pxyz[0], pxyz[1], pxyz[2], masscsv, charge);
    }
    // kaons
    double ntrK = gRandom->Poisson(fNChK * periferal);
    printf("fNChK=%f ntrK=%f\n", fNChK, ntrK);
    for (int itr = 0; itr < ntrK; itr++)
    {
      double yrap = fdNdYK->GetRandom();
      double pt = fdNdPtK->GetRandom();
      double phi = gRandom->Rndm() * TMath::Pi() * 2;
      float charge = gRandom->Rndm() > 0.3 ? 1 : -1;
      int ptypecsv = 321;
      double masscsv = 0.493677;
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + kMassK * kMassK) * TMath::SinH(yrap)};

      std::string binaryString = std::bitset<12>(event + 1).to_string() + std::string(12, '0') + std::bitset<16>(particleNumber++).to_string() + std::string(24, '0');
      unsigned long long barcode = std::stoull(binaryString, nullptr, 2);
      fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, ptypecsv, 0, x, y, z, 0., pxyz[0], pxyz[1], pxyz[2], masscsv, charge);
    }
    // protons
    double ntrP = gRandom->Poisson(fNChP * periferal);
    printf("fNChP=%f ntrP=%f\n", fNChP, ntrP);
    for (int itr = 0; itr < ntrP; itr++)
    {
      double yrap = fdNdYP->GetRandom();
      double pt = fdNdPtP->GetRandom();
      double phi = gRandom->Rndm() * TMath::Pi() * 2;
      double charge = 1;
      int ptypecsv = 2212;
      double masscsv = 0.938272119;
      double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + kMassP * kMassP) * TMath::SinH(yrap)};
      std::string binaryString = std::bitset<12>(event + 1).to_string() + std::string(12, '0') + std::bitset<16>(particleNumber++).to_string() + std::string(24, '0');
      unsigned long long barcode = std::stoull(binaryString, nullptr, 2);
      fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, ptypecsv, 0, x, y, z, 0., pxyz[0], pxyz[1], pxyz[2], masscsv, charge);
    }
  }
  if (addD0 || addSec)
  {
    std::vector<int> vPdg = {310, 3122};

    int nMomGen = 0;
    for (auto pdg : vPdg)
    {
      int avTot = GetMultiplicity(pdg, Eint, true) * GetBRatio(pdg);
      double mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
      int totMoms = gRandom->Poisson(avTot * periferal);

      bool matter = true;

      fpt = new TF1("fpt", "x*exp(-TMath::Sqrt(x**2+[0]**2)/[1])", 0, 100);
      fpt->SetParameter(0, mass);
      fpt->SetParameter(1, GetTslope(TMath::Abs(pdg), Eint, matter) / 1000);

      // rapidity distribution
      fy = new TF1("fy", " exp(-0.5*((x-[0]-[2])/[1])**2)+exp(-0.5*((x+[0]-[2])/[1])**2)", -20, 20);
      fy->SetParameter(0, GetY0Rapidity(TMath::Abs(pdg), Eint, matter));
      fy->SetParameter(1, GetSigmaRapidity(TMath::Abs(pdg), Eint, matter));
      fy->SetParameter(2, GetY0(Eint));

      for (int nMom = 1; nMom <= totMoms; nMom++)
      {
        nMomGen++;
        double ptGenD = addSec ? fpt->GetRandom() : hD0pt->GetRandom();
        double yGenD = addSec ? fy->GetRandom() : hD0y->GetRandom();
        double phi = gRandom->Rndm() * 2 * TMath::Pi();
        double pxGenD = ptGenD * TMath::Cos(phi);
        double pyGenD = ptGenD * TMath::Sin(phi);

        double mt = TMath::Sqrt(ptGenD * ptGenD + mass * mass);
        double pzGenD = mt * TMath::SinH(yGenD);
        double en = mt * TMath::CosH(yGenD);

        std::string binaryString = std::bitset<12>(event + 1).to_string() + std::string(12, '0') + std::bitset<16>(particleNumber++).to_string() + std::string(24, '0');
        unsigned long long barcode = std::stoull(binaryString, nullptr, 2);
        TLorentzVector *mom = new TLorentzVector();
        mom->SetPxPyPzE(pxGenD, pyGenD, pzGenD, en);
        fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, pdg, nMom, x, y, z, 0., pxGenD, pyGenD, pzGenD, mass, 0.0);
        TClonesArray *particles = new TClonesArray("TParticle", 10);
        Int_t np;
        do
        {
          fDecayer->Decay(pdg, mom);
          np = fDecayer->ImportParticles(particles);
        } while (np < 0);

        // loop on decay products
        for (int i = 1; i < np; i++)
        {
          if (i == 3)
            break;
          TParticle *iparticle1 = (TParticle *)particles->At(i);
          double charge = iparticle1->GetPdgCode() / TMath::Abs(iparticle1->GetPdgCode());
          double vX = iparticle1->Vx();
          double vY = iparticle1->Vy();
          double vZ = iparticle1->Vz();
          double px = iparticle1->Px();
          double py = iparticle1->Py();
          double pz = iparticle1->Pz();
          int ptypecsv = TMath::Abs(iparticle1->GetPdgCode());
          double masscsv = TDatabasePDG::Instance()->GetParticle(ptypecsv)->Mass();
          binaryString = std::bitset<12>(event + 1).to_string() + std::bitset<12>(nMomGen).to_string() + std::bitset<16>(particleNumber++).to_string() + std::bitset<8>(0).to_string() + std::bitset<16>(0 * int(i - 1)).to_string();
          barcode = std::stoull(binaryString, nullptr, 2);
          fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, ptypecsv, nMom, x + vX, y + vY, z + vZ, 0., px, py, pz, masscsv, charge);
        }
      }
    }
  }
}
