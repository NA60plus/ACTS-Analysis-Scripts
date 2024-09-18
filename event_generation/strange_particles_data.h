#ifndef MEASUREMENT_UTILS_H
#define MEASUREMENT_UTILS_H
#endif

#include <TString.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TMath.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TParticle.h"
#include "AliAODRecoDecay.h"
#include "AliDecayer.h"
#include "AliDecayerEvtGen.h"
#include "TDatabasePDG.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

// settings for signal generation
// rapidity range
const double yminSG = -2.; 
const double ymaxSG = 8.;
// pT range
const double ptminSG = 0;
const double ptmaxSG = 10; 

//number of energy of the colliding beam
const int NEnergy = 5;
//number of particles 
const int NParticles = 5;
//array of the energy
double Elab[NEnergy] = {20,30,40,80,158};

//T parameter in the exponential pT distribution [matter/antimatter][particle][beam energy] i
//                                                    phi                        K               Lambda                Omega                    Xi
double Tslope[2][NParticles][NEnergy] = {{{196.8,237.4,244.6,239.8,298.7},{0,0,229, 223.1, 229},{244,249,258,265,301},{0,0,218,0,267},{221,233,222,227,277}},//matter e rapidity distribution [matter/antimatter][particle][beam energy]
                                         {{196.8,237.4,244.6,239.8,298.7},{0,0,226,217,226},{339,284,301,292,303},{0,0,218,0,259},{311,277,255,321,0}}};//antimatter
//sigma parameter of the gaussians of th                    phi                        K                       Lambda                Omega                    Xi
double sigma_rapidity[2][NParticles][NEnergy] = {{{0.425,0.538,0.696,0.658,1.451},{0,0,0.725,0.792,0.88},{0.51,0.66,0.91,0.87,0},{0,0,0.6,0,1.2},{0.45,0.56,0.76,0.71,1.18}},//matter
                                                 {{0.425,0.538,0.696,0.658,1.451},{0,0,0.635,0.705,0.81},{0.62,0.69,0.77,0.83,1.00},{0,0,0.6,0,1.0},{0,0.76,0.65,0.87,0.73}}};//antimatter
//mu paramter of the gaussian of the rapidity distribution [matter/antimatter][particle][beam energy]
//                                                       phi                        K                   Lambda                   Omega                    Xi
double y0_rapidity[2][NParticles][NEnergy] = {{{0.425,0.538,0.487,0.682,0.},{0,0,0.694,0.742,0.839},{0.49,0.59,0.65,0.94,0},{0,0,0,0,0},{0.45,0.47,0.54,0.68,0}},//matter
                                              {{0.425,0.538,0.487,0.682,0.},{0,0,0.569,0.668,0.727},{0,0,0,0,0},{0,0,0.0,0},{0,0,0,0,0}}};//antimatter
//multiplicity for event [matter/antimatter][particle][beam energy]
//                                                       phi                        K                   Lambda                   Omega                    Xi
double multiplicity[2][NParticles][NEnergy] = {{{1.89,1.84,2.55,4.04,8.46},{0,0,59.1,76.9,103.0},{27.1,36.9,43.1,50.1,44.9},{0,0,0.14,0,0.43},{1.50,2.42,2.96,3.80,4.04}},//matter
                                               {{1.89,1.84,2.55,4.04,8.46},{0,0,19.2,32.4,51.9},{0.16,0.39,0.68,1.82,3.07},{0,0,0.14,0,0.19},{0,0.12,0.13,0.58,0.66}}};//antimatter
//branching ratio [particle]
//                            phi     K   Lambda  Omega   Xi
double bratio[NParticles] = {0.489, 0.692, 0.639, 0.433, 0.638};//for Omega (Xi) it's given by 0.678*0.639 (0.999*0.639)
//name of the particles
TString particle_name[NParticles] = {"phi","K0s","Lambda","Omega","Xi"};
TString all_particle_name[NParticles+3] = {"pion","kaon","proton","phi","K0s","Lambda","Omega","Xi"};
//pdg code of the particles
double pdg_mother[NParticles] = {333,310,3122,3334,3312};
//pdg code of the daughters
//                                     K+    K-   pi+  pi-    p    pi-  Lambda  K-  Lambda pi-
double pdg_daughter[NParticles][2] = {{321,-321},{211,-211},{2212,-211},{3122,-321},{3122,-211}};


double GetTslope(int pdgParticle, double Eint, bool matter = true){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E){
      index_E = counter;
      break;
    }
    counter++;
  }
  if(TMath::Abs(pdgParticle)==310)
    return (Tslope[0][index_pdg][index_E]+Tslope[1][index_pdg][index_E])/2.;
  else
    return Tslope[matter ? 0 : 1][index_pdg][index_E];
}

double GetSigmaRapidity(int pdgParticle, double Eint, bool matter = true){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E){
      index_E = counter;
      break;
    }
    counter++;
  }
  
  if(TMath::Abs(pdgParticle)==310)
    return (sigma_rapidity[0][index_pdg][index_E]+sigma_rapidity[1][index_pdg][index_E])/2.;
  else
    return sigma_rapidity[matter ? 0 : 1][index_pdg][index_E];
}

double GetBRatio(int pdgParticle){
  int index_pdg = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  return bratio[index_pdg];
}

double GetY0Rapidity(int pdgParticle, double Eint, bool matter = true){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E){
      index_E = counter;
      break;
    }
    counter++;
  }
  if(TMath::Abs(pdgParticle)==310)
    return (y0_rapidity[0][index_pdg][index_E]+y0_rapidity[1][index_pdg][index_E])/2.;
  else
    return y0_rapidity[matter ? 0 : 1][index_pdg][index_E];
}

double GetMultiplicity(int pdgParticle, double Eint, bool matter = true){
  int index_pdg = 0;
  int index_E = 0;
  int counter = 0;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      break;
    }
    counter++;
  }
  counter = 0;
  for(auto& E: Elab){
    if(Eint == E){
      index_E = counter;
      break;
    }
    counter++;
  }
  
  if(TMath::Abs(pdgParticle)==310)
    return (multiplicity[0][index_pdg][index_E]+multiplicity[1][index_pdg][index_E])/2.;
  else
    return multiplicity[matter ? 0 : 1][index_pdg][index_E];
}

double GetY0(double Eint){
  double y0BG = 2.9;
  if (Eint == 20)
    y0BG = 1.9;   // gaussian y mean - 20 GeV
  else if (Eint == 40)
    y0BG = 2.22; // gaussian y mean - 40 GeV
  else if (Eint == 60)
    y0BG = 2.42;  // gaussian y mean - 60 GeV
  else if (Eint == 80)
    y0BG = 2.57;  // gaussian y mean - 80 GeV
  else if (Eint == 158)
    y0BG = 2.9;   // gaussian y mean - 160 GeV
  else if (Eint == 400)
    y0BG = 3.37;  // gaussian y mean - 400 GeV
  return y0BG;
}


void GetPDGDaughters(int pdgParticle, int pdgDaughters[], bool matter = true){
  int index_pdg = 0;
  int counter = 0;
  int sign = matter ? 1 : -1;
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg){
      index_pdg = counter;
      for(int i = 0; i < 2; i++)
        pdgDaughters[i] = pdg_daughter[index_pdg][i]*sign;
      return;
    }
    counter++;
  }
}

int GetArrayPosition(int pdgParticle){
  int counter = 3;

  if(TMath::Abs(pdgParticle) == 211) 
    return 0;
  if(TMath::Abs(pdgParticle) == 321) 
    return 1;
  if(TMath::Abs(pdgParticle) == 2212) 
    return 2;
  
  for(auto& pdg: pdg_mother){
    if(pdgParticle == pdg)
      return counter;
    counter++;
  }

  return NParticles+3;


}

int GetNBody(int pdgParticle, int& pdg_unstable_dau, bool matter = true){
  pdgParticle = TMath::Abs(pdgParticle);
  int pdgDaughters[2];
  GetPDGDaughters(pdgParticle, pdgDaughters, true);
  int sign = matter ? 1 : -1;
  for(auto& pdg: pdg_mother){
    if(pdgDaughters[0] == pdg){
      pdg_unstable_dau = pdg*sign;
      return 3;
    }
    else if(pdgDaughters[1] == pdg){
      pdg_unstable_dau = pdg*sign;
      return 3;
    }
  }

  return 2;
}
