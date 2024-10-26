#if !defined(__CINT__) || defined(__MAKECINT__)
// #define _NOALIROOT_
#include "KMCDetectorFwd.h"
#include "KMCProbeFwd.h"
#include "KMCFlukaParser.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TROOT.h"
#include "GenMUONLMR.h"
#include "TStopwatch.h"
#include "TTree.h"

// settings for signal generation
double MotherMass;
int generYPtPar, ProcType;
// double y0SG   = 1.9;   // gaussian y mean - 20 GeV
// double y0SG   = 2.08;   // gaussian y mean - 30 GeV
double y0SG = 2.22;   // gaussian y mean - 40 GeV
double sigySG = 1.;   // .. sigma
double yminSG = -10.; // min y to generate
double ymaxSG = 10.;  //
double TSG;           // inv.slope of thermal pt distribution
double ptminSG = 0.01;
double ptmaxSG = 3;

void SetProcessParameters(const Char_t *Process, Double_t E);

void runDiMuGenLMR(const Char_t *Process = "jpsi", Double_t Eint)
{
  SetProcessParameters(Process, TMath::Abs(Eint));
  GenMUONLMR *gener = new GenMUONLMR(7, 1);
  //      0         1        2         3         4        5           6         7
  //  "fPtPion","fPtKaon","fPtEta","fPtRho","fPtOmega","fPtPhi","fPtEtaPrime"  "fPtJPsi"
  gener->SetYParams(generYPtPar, 1., y0SG, sigySG, 0.);
  gener->SetPtParams(generYPtPar, 1., TSG, MotherMass, 0.);
  // 2 = rho
  // 3 = omega 2 body
  // 4 = omega Dalitz
  //
  // Processes eta 2B=0  eta D=1  rho =2   omega 2B=3  omega D=4  phi=5    eta p=6       pi=7      K=8      J/psi=9
  //           kEtaLMR   kEtaLMR  kRhoLMR  kOmegaLMR   kOmegaLMR  kPhiLMR  kEtaPrimeLMR  kPionLMR  kK        kJPsi
  //
  //  BR       5.8e-6    3.1e-4   4.55e-5  7.28e-5     1.3e-4     2.86e-4  1.04e-4       1         0.6344    0.05

  // Added June 2019 for J/psi generation
  if (strcmp(Process, "jpsisch") == 0)
  {
    gener->SetYParams(generYPtPar, 1., p1_xF, p2_xF, p3_xF, p4_xF);
  }
  else
  {
    gener->SetYParams(generYPtPar, 1., y0SG, sigySG, 0., 0.);
  }

  if (strcmp(Process, "jpsi") == 0 || strcmp(Process, "jpsisch") == 0)
  {
    gener->SetPtParams(generYPtPar, 1., par_p0, par_n, par_n2);
  }
  else
  {
    gener->SetPtParams(generYPtPar, 1., TSG, MotherMass, 0.);
  }

  gener->GenerateSingleProcess(ProcType);
  gener->Generate();

  for (int i = 0; i < 2; i++)
  {

    TParticle *iparticle = gener->GetMuon(i);
    double charge = iparticle->GetPdgCode() / TMath::Abs(iparticle->GetPdgCode());
    double vX = iparticle->Vx();
    double vY = iparticle->Vy();
    double vZ = iparticle->Vz();
    double px = iparticle->Px();
    double py = iparticle->Py();
    double pz = iparticle->Pz();
    int ptypecsv = TMath::Abs(iparticle->GetPdgCode());
    double masscsv = TDatabasePDG::Instance()->GetParticle(ptypecsv)->Mass();
    binaryString = std::bitset<12>(event + 1).to_string() + std::bitset<12>(nMomGen).to_string() + std::bitset<16>(particleNumber++).to_string() + std::bitset<8>(0).to_string() + std::bitset<16>(0 * int(i - 1)).to_string();
    barcode = std::stoull(binaryString, nullptr, 2);

    fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, ptypecsv, nMom, x + vX, y + vY, z + vZ, 0., px, py, pz, masscsv, charge);
    fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, ptypecsv, nMom, x + vX, y + vY, z + vZ, 0., px, py, pz, masscsv, charge);
  }

  //    //printf("px0=%f py0=%f pz0=%f e0=%f m0=%f\n",fMu[0]->Px(),fMu[0]->Py(),fMu[0]->Pz(),fMu[0]->Energy(),fMu[0]->GetMass());
  //    //printf("px1=%f py1=%f pz1=%f e1=%f m1=%f\n",fMu[1]->Px(),fMu[1]->Py(),fMu[1]->Pz(),fMu[1]->Energy(),fMu[1]->GetMass());
}

void SetProcessParameters(const Char_t *Process, Double_t E)
{
  // modified June2019 for J/psi parameterization
  if (E == 20.)
  {
    y0SG = 1.9;
    sigySG = 1.;
  }
  else if (E == 30.)
  {
    y0SG = 2.08;
    sigySG = 1.;
  }
  else if (E == 40.)
  {
    y0SG = 2.22;
    if (strcmp(Process, "jpsi") == 0)
      sigySG = 0.19; // PYTHIA 50 GeV
    else
      sigySG = 1.;
  }
  else if (E == 50.)
  {
    y0SG = 2.33;
    if (strcmp(Process, "jpsi") == 0)
      sigySG = 0.19; // PYTHIA 50 GeV
  }
  else if (E == 60.)
  {
    y0SG = 2.42;
    sigySG = 1.;
  }
  else if (E == 70.)
  {
    y0SG = 2.50;
    if (strcmp(Process, "jpsi") == 0)
      sigySG = 0.28; // PYTHIA 70 GeV
  }
  else if (E == 80.)
  {
    y0SG = 2.57;
    sigySG = 1.;
  }
  else if (E == 90.)
  {
    y0SG = 2.63;
    if (strcmp(Process, "jpsi") == 0)
      sigySG = 0.33; // PYTHIA 90 GeV
  }
  else if (E == 110.)
  {
    y0SG = 2.73;
    if (strcmp(Process, "jpsi") == 0)
      sigySG = 0.37; // PYTHIA 110 GeV
  }
  else if (E == 130.)
  {
    y0SG = 2.81;
    if (strcmp(Process, "jpsi") == 0)
      sigySG = 0.40; // PYTHIA 130 GeV
  }
  else if (E == 150.)
  {
    y0SG = 2.88;
    if (strcmp(Process, "jpsi") == 0)
      sigySG = 0.42; // PYTHIA 150 GeV
  }
  else if (E == 160.)
  {
    y0SG = 2.9;
    if (strcmp(Process, "jpsi") == 0)
      sigySG = 0.42; // PYTHIA
    else
      sigySG = 1.2; // NA49
  }
  else if (E == 400.)
  {
    y0SG = 3.37;
    sigySG = 1.;
  }

  if (strcmp(Process, "etaDalitz") == 0)
  { // eta Dalitz
    ProcType = 1;
    generYPtPar = 2;
    MotherMass = 0.549;
    if (E == 40.)
    {
      TSG = 0.225;
    }
    else if (E == 160.)
    {
      TSG = 0.24;
    }
  }
  else if (strcmp(Process, "eta2Body") == 0)
  { // eta Dalitz
    ProcType = 0;
    generYPtPar = 2;
    MotherMass = 0.549;
    if (E == 40.)
    {
      TSG = 0.225;
    }
    else if (E == 160.)
    {
      TSG = 0.24;
    }
  }
  else if (strcmp(Process, "rho") == 0)
  { // rho
    ProcType = 2;
    generYPtPar = 3;
    MotherMass = 0.775;
    if (E == 40.)
    {
      TSG = 0.25;
    }
    else if (E == 160.)
    {
      TSG = 0.29;
    }
  }
  else if (strcmp(Process, "omega2Body") == 0)
  { // omega 2 body
    ProcType = 3;
    generYPtPar = 4;
    MotherMass = 0.781;
    if (E == 40.)
    {
      TSG = 0.25;
    }
    else if (E == 160.)
    {
      TSG = 0.29;
    }
  }
  else if (strcmp(Process, "omegaDalitz") == 0)
  { // omega Dalitz
    ProcType = 4;
    generYPtPar = 4;
    MotherMass = 0.781;
    if (E == 40.)
    {
      TSG = 0.25;
    }
    else if (E == 160.)
    {
      TSG = 0.29;
    }
  }
  else if (strcmp(Process, "phi") == 0)
  { // phi
    ProcType = 5;
    generYPtPar = 5;
    MotherMass = 1.02;
    if (E == 40.)
    {
      TSG = 0.25;
    }
    else if (E == 160.)
    {
      TSG = 0.30; // NA49
    }
  }
  else if (strcmp(Process, "etaPrime") == 0)
  { // eta prime
    ProcType = 6;
    generYPtPar = 6;
    MotherMass = 0.958;
    if (E == 40.)
    {
      TSG = 0.25;
    }
    else if (E == 160.)
    {
      TSG = 0.29;
    }
    // modified June 2019 J/psi parameterizations
  }
  else if (strcmp(Process, "jpsi") == 0 || strcmp(Process, "jpsisch") == 0)
  { // J/psi
    ProcType = 9;
    if (strcmp(Process, "jpsi") == 0)
    {
      ProcType = 9;
      generYPtPar = 7;
    }
    else if (strcmp(Process, "jpsisch") == 0)
    {
      ProcType = 10;
      generYPtPar = 8;
    }
    MotherMass = 3.0969;
    //     if(E == 40.) {
    //       TSG = 0.25;
    //     }
    //     else if(E == 160.){
    //       TSG = 0.284; //NA50
    //     }
    if (E == 150.)
    {
      par_p0 = 5.34;
      par_n = 18.8;
      par_n2 = 2.70;
    }
    else if (E == 130.)
    {
      par_p0 = 13.96;
      par_n = 201.1;
      par_n2 = 2.58;
    }
    else if (E == 110.)
    {
      par_p0 = 3.73;
      par_n = 8.83;
      par_n2 = 2.84;
    }
    else if (E == 90.)
    {
      par_p0 = 4.89;
      par_n = 17.3;
      par_n2 = 2.80;
    }
    else if (E == 70.)
    {
      par_p0 = 3.86;
      par_n = 10.6;
      par_n2 = 2.92;
    }
    else if (E == 50.)
    {
      par_p0 = 11.03;
      par_n = 429.7;
      par_n2 = 3.27;
    }
    else if (E == 40.)
    { // use same parameter as for E=40 GeV
      par_p0 = 11.03;
      par_n = 429.7;
      par_n2 = 3.27;
    }
  }

  if (strcmp(Process, "jpsisch") == 0)
  { // J/psi
    ProcType = 10;
    generYPtPar = 8;
    MotherMass = 3.0969;

    Double_t a = 13.5;
    Double_t b = 44.9;
    p4_xF = y0SG;
    p1_xF = 3.096; // mjpsi
    if (E == 150.)
    {
      p2_xF = 16.8; // sqrts
    }
    else if (E == 130.)
    {
      p2_xF = 15.7;
    }
    else if (E == 110.)
    {
      p2_xF = 14.4;
    }
    else if (E == 90.)
    {
      p2_xF = 13.1;
    }
    else if (E == 70.)
    {
      p2_xF = 11.5;
    }
    else if (E == 50.)
    {
      p2_xF = 9.77;
    }
    else if (E == 40.)
    {
      p2_xF = 8.76;
    }
    p3_xF = a / (1 + b / p2_xF);
    // printf("p1_xF=%f, p2_xF=%f, p3_xF=%f, p4_xF=%f\n", p1_xF, p2_xF, p3_xF, p4_xF);
  }

  // printf("Generating %s\n", Process);
  // printf("Particle mass = %f\n", MotherMass);
  // printf("Inverse slope = %f\n", TSG);
  // printf("rapidity sigma = %f\n", sigySG);
}