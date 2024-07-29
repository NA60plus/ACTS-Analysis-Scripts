#include <cmath>
#include <iostream>
#include <string>
#include <bitset>
#include <TF1.h>
#include <TRandom.h>

FILE *fptxt;
FILE *fptxtZRotated;

void BFieldNA60Plus(int nstepx = 100, int nstepy = 100, int nstepz = 400, bool addMNP33 = true, bool addMEP48 = true)
{
  //================================================================================
  // MNP33 field map from https://inspirehep.net/files/4706e9a213757975c85c3fdbedc6f5bc
  // fig 6
  // filed centered in 0,0,0 --> to be moved in the correct position
  Double_t ByMEP48 = 1.5; // T

  const int nn = 56;
  Double_t z[nn] = {
      -12000., -11000., -10000., -9000., -8000., -7000., -6000., -5000.,
      -3983.0, -3548.9, -2987.2, -2485.1, -2025.5, -1795.7, -1600.0, -1404.3,
      -1268.1, -1166.0, -1080.9, -987.2, -868.1, -774.5, -672.3, -544.7,
      -459.6, -374.5, -263.8, -161.7, -8.5, 144.7, 263.8, 323.4,
      391.5, 459.6, 510.6, 570.2, 672.3, 723.4, 791.5, 876.6,
      961.7, 1055.3, 1140.4, 1268.1, 1395.7, 1548.9, 1719.1, 1914.9,
      2093.6, 2383.0, 2383.0, 2748.9, 3131.9, 3506.4, 3966.0, 5000.}; // mm
  Double_t By[nn] = {
      0., 0., 0., 0., 0., 0., 0., 0.,
      0.001, 0.001, 0.006, 0.016, 0.033, 0.047,
      0.069, 0.095, 0.123, 0.148, 0.174, 0.198,
      0.229, 0.262, 0.297, 0.332, 0.351, 0.378,
      0.396, 0.407, 0.414, 0.409, 0.395, 0.381,
      0.368, 0.348, 0.330, 0.313, 0.289, 0.267,
      0.247, 0.221, 0.195, 0.174, 0.152, 0.122,
      0.098, 0.074, 0.055, 0.040, 0.027, 0.017,
      0.017, 0.010, 0.003, 0.001, 0.000, 0.}; // T

  Double_t zorigin = 3850; // center of the magnet[mm]

  for (int i = 0; i < nn; i++)
    z[i] += zorigin;

  TGraph *gg = new TGraph(nn, z, By);
  gg->SetMarkerColor(2);
  gg->SetMarkerStyle(20);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TH2D *hnull = new TH2D("hnull", "hnull", 100, -8000, 8000, 100, 0., 0.5);
  hnull->GetXaxis()->SetTitle("z[mm]");
  hnull->GetYaxis()->SetTitle("By[T]");
  hnull->Draw();
  gg->Draw("P");
  TF1* fitB = new TF1("fitB","gaus(0)+gaus(3)",z[0],5000+zorigin);
  fitB->SetParameter(0,0.298);
  fitB->SetParameter(1,zorigin);
  fitB->SetParameter(2,820);
  fitB->SetParameter(3,0.119);
  fitB->SetParameter(4,zorigin);
  fitB->SetParameter(5,820*2);
  gg->Fit(fitB,"MR+");

  Double_t zz;
  // for(int i=0;i<16000;i+=500){
  //   zz = -8000+i;
  //   printf("i= %f, val = %f\n",zz,gg->Eval(zz));
  // }

  //================================================================================
  // write field map
  // format of field map x, y, z, bx, by, bz

  char txtname[30];
  sprintf(txtname, "BField_MNP33.txt");
  fptxt = fopen(txtname, "w");

  char txtnameZRotated[30];
  sprintf(txtnameZRotated, "BFieldZRotated_MNP33.txt");
  fptxtZRotated = fopen(txtnameZRotated, "w");

  TH1D *hField = new TH1D("hField", "hField", nstepz, -8000., 8000.);
  hField->GetXaxis()->SetTitle("z[mm]");
  hField->GetYaxis()->SetTitle("By[T]");

  // B field map should be defined over the full TGeo volume
  Double_t xmintot = -3080.; // [mm]
  Double_t xmaxtot = 3080.;  // [mm]

  Double_t ymintot = -3080.; // [mm]
  Double_t ymaxtot = 3080.;  // [mm]

  Double_t zmintot = -8000.; // [mm]
  Double_t zmaxtot = 8050.;  // [mm]  // to be sure the map covers up to 8000

  Double_t stepx = (xmaxtot - xmintot) / nstepx;
  printf("number of steps in x = %d, width= %4.3f\n", nstepx, stepx);
  Double_t stepy = (ymaxtot - ymintot) / nstepy;
  printf("number of steps in y = %d, width= %4.3f\n", nstepy, stepy);
  Double_t stepz = (zmaxtot - zmintot) / nstepz;
  printf("number of steps in z = %d, width= %4.3f\n", nstepz, stepz);

  Double_t xval = 0;
  Double_t yval = 0;
  Double_t zval = 0;
  Double_t bxval = 0;
  Double_t byval = 0;
  Double_t bzval = 0;

  for (int ix = 0; ix < nstepx; ix++)
  {
    xval = xmintot + ix * stepx;
    for (int iy = 0; iy < nstepy; iy++)
    {
      yval = ymintot + iy * stepy;
      for (int iz = 0; iz < nstepz; iz++)
      {
        zval = zmintot + iz * stepz;
        bxval = 0;
        bzval = 0;
        byval = (Double_t)fitB->Eval(zval);
        if (zval < 840 && zval > -150)
          byval += ByMEP48;
        fprintf(fptxt, "%f %f %f %f %f %f\n", xval, yval, zval, bxval, byval, bzval);
        hField->SetBinContent(hField->FindBin(zval), byval);
      }
    }
  }

  for (int ix = 0; ix < nstepx; ix++)
  {
    for (int iz = 0; iz < nstepz; iz++)
    {
      for (int iy = nstepy - 1; iy >= 0; iy--)
      {
        xval = xmintot + ix * stepx;
        yval = ymintot + iy * stepy;
        zval = zmintot + iz * stepz;

        bxval = 0;
        bzval = 0;
        byval = (Double_t)gg->Eval(zval);
        fprintf(fptxtZRotated, "%f %f %f %f %f %f\n", xval, zval, -yval, bxval, bzval, -byval);
      }
    }
  }

  TCanvas *cc = new TCanvas("cc", "cc", 20, 20, 400, 400);
  hField->Draw("colz");

  fclose(fptxt);
  fclose(fptxtZRotated);
  TFile *f = new TFile("BFieldNA60plus.root", "recreate");
  hField->Write();
  f->Close();
}
