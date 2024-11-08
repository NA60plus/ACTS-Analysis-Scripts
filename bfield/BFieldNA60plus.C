#include <cmath>
#include <iostream>
#include <string>
#include <bitset>
#include <TF1.h>
#include <TRandom.h>

FILE *fptxt;
FILE *fptxtZRotated;

void BFieldNA60plus(int nstepx = 200, int nstepy = 200, int nstepz = 400, bool addMNP33 = true, bool addMEP48 = true)
{
  //================================================================================
  // B field map should be defined over the full TGeo volume

  Double_t xmintot = -4000.; // [mm]
  Double_t xmaxtot = 4000.;  // [mm]

  Double_t ymintot = -4000.; // [mm]
  Double_t ymaxtot = 4000.;  // [mm]

  Double_t zmintot = -9200.; // [mm]
  Double_t zmaxtot = 9200.;  // [mm]  // to be sure the map covers up to 9000

  // MNP33 field map from https://inspirehep.net/files/4706e9a213757975c85c3fdbedc6f5bc fig 6
  // field centered in 0,0,0 --> to be moved in the correct position

  Double_t ByMEP48 = 1.5;                                     // T
  Double_t ByMNP33 = 0.746 / (5900. - 3000) * (5100. - 3800); // the integrated field in the region of the first 4 chambers (3000-5900) is now applied to the region of the B field (3800-5100)
  Double_t TotByMNP33 = 0.864;                                // over the full range -270, +270

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

  Double_t zorigin = 4450; // center of the magnet[mm] for the long setup
                           //  Double_t zorigin = 3850; // center of the magnet[mm] for the ultrashort setup

  for (int i = 0; i < nn; i++)
    z[i] += zorigin;

  TGraph *gg = new TGraph(nn, z, By);
  gg->SetMarkerColor(2);
  gg->SetMarkerStyle(20);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TH2D *hnull = new TH2D("hnull", "hnull", 100, zmintot, zmaxtot, 100, 0., 0.5);
  hnull->GetXaxis()->SetTitle("z[mm]");
  hnull->GetYaxis()->SetTitle("By[T]");
  hnull->Draw();
  gg->Draw("P");

  TF1 *fitB = new TF1("fitB", "gaus(0)+gaus(3)", zmintot, zmaxtot);
  //  TF1* fitB = new TF1("fitB","gaus(0)+gaus(3)",z[0],5000+zorigin); //modified
  fitB->SetParameter(0, 0.298);
  fitB->SetParameter(1, zorigin);
  fitB->SetParameter(2, 820);
  fitB->SetParameter(3, 0.119);
  fitB->SetParameter(4, zorigin);
  fitB->SetParameter(5, 820 * 2);
  gg->Fit(fitB, "MR+");

  printf("Integral MNP33 (full range)= %f\n", fitB->Integral(0., 8000));
  printf("Integral MNP33 (3000-5300) = %f\n", fitB->Integral(3000., 5300));
  printf("Integral MNP33 (3000-5900) = %f\n", fitB->Integral(3000., 5900));
  printf("Integral MNP33 (3000-8500) = %f T\n", fitB->Integral(3000., 8500));
  printf("Integrated B field in T between 3000-5300 (first 3 chambers) = %f\n", fitB->Integral(3000., 5300) / (5300. - 3000));

  //================================================================================
  // write field map
  // format of field map x, y, z, bx, by, bz

  char txtname[40];
  //  sprintf(txtname, "BFieldNA60plus_MNP33integrated_longsetup.txt");
  sprintf(txtname, "BFieldNA60plus_longsetup.txt");
  fptxt = fopen(txtname, "w");

  char txtnameZRotated[40];
  //  sprintf(txtnameZRotated, "BFieldNA60plus_ZRotated_MNP33integrated_longsetup.txt");
  sprintf(txtnameZRotated, "BFieldNA60plus_ZRotated_longsetup.txt");
  fptxtZRotated = fopen(txtnameZRotated, "w");

  TH1D *hField = new TH1D("hField", "", nstepz, -100., 9000.);
  hField->GetXaxis()->SetTitle("z[mm]");
  hField->GetYaxis()->SetTitle("By[T]");

  TH1D *hFieldRotated = new TH1D("hField", "hField", nstepz, -100., 9000.);
  hFieldRotated->GetXaxis()->SetTitle("y_new = z[mm]");
  hFieldRotated->GetYaxis()->SetTitle("Bz[T]");

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
        // byval = 0;
        byval = (Double_t)fitB->Eval(zval);
        if (zval > -155 && zval < 845)
          byval += ByMEP48;
        // if (zval>3800 && zval <5100) byval= ByMNP33; // if integrated field

        fprintf(fptxt, "%f %f %f %f %f %f\n", xval, yval, zval, bxval, byval, bzval);
        //hField->SetBinContent(hField->FindBin(zval), byval);
      }
    }
  }

  // rotated coordinates
  for (int ix = 0; ix < nstepx; ix++)
  {
    xval = xmintot + ix * stepx;
    for (int iz = 0; iz < nstepz; iz++)
    {
      zval = zmintot + iz * stepz;
      for (int iy = nstepy - 1; iy >= 0; iy--)
      { // originally y were decreasing
        yval = ymintot + iy * stepy;

        bxval = 0;
        bzval = 0;
        // byval = 0;
        byval = (Double_t)fitB->Eval(zval);
        if (zval > -155 && zval < 845)
          byval += ByMEP48;
        // if (zval>3800 && zval <5100) byval= ByMNP33; // if integrated field

        fprintf(fptxtZRotated, "%f %f %f %f %f %f\n", xval, zval, -yval, bxval, bzval, -byval);
        hFieldRotated->SetBinContent(hFieldRotated->FindBin(zval), -byval);
      }
    }
  }


  for (int iz = 1; iz < nstepz+1; iz++)
  {
    zval = hField->GetBinCenter(iz);
    byval = (Double_t)fitB->Eval(zval);
    if (zval > -155 && zval < 845)
      byval += ByMEP48;
    hField->SetBinContent(hField->FindBin(zval), byval);
  }


  //   for (int ix = 0; ix < nstepx; ix++) {
  //     xval = xmintot + ix * stepx;
  //     for (int iz = 0; iz < nstepz; iz++) {
  //       zval = zmintot + iz * stepz;
  //       for (int iy = nstepy - 1; iy >= 0; iy--) {
  //         yval = ymintot + iy * stepy;
  //
  //         bxval = 0;
  //         bzval = 0;
  //
  //         //byval = (Double_t)fitB->Eval(zval); //modified
  //         if (zval > -155 && zval < 845) {
  // 	  byval = ByMEP48;  //ADDED
  // 	} else 	if(zval>3800 && zval <5100) {
  // 	  byval = 0.555/1.3; // integrate B field (B = 0.555 TM)
  // 	} else byval = 0;
  //
  //         //fprintf(fptxtZRotated, "%f %f %f %f %f %f\n", xval, zval, -yval, bxval, bzval, -byval);
  //  	hFieldRotated->SetBinContent(hFieldRotated->FindBin(zval),-byval);
  //       }
  //     }
  //   }

  TCanvas *cc = new TCanvas("cc", "cc", 20, 20, 400, 400);
  cc->Divide(1, 2);
  cc->cd(1);
  hField->Draw("colz");
  cc->cd(2);
  hFieldRotated->Draw("colz");

  TCanvas *cc2 = new TCanvas("cc2", "cc2", 20, 20, 1600, 1200);
  hField->Draw("colz");
  TLine vt1(71.175, 0, 71.175, 1.55);
  TLine vt2(151.175, 0, 151.175, 1.55);
  TLine vt3(201.175, 0, 201.175, 1.55);
  TLine vt4(251.175, 0, 251.175, 1.55);
  TLine vt5(381.175, 0, 381.175, 1.55);

  TLine ms1(2999.695, 0, 2999.695, 1.55);
  TLine ms2(3599.695, 0, 3599.695, 1.55);
  TLine ms3(5299.695, 0, 5299.695, 1.55);
  TLine ms4(5899.695, 0, 5899.695, 1.55);
  TLine ms5(8099.695, 0, 8099.695, 1.55);
  TLine ms6(8499.695, 0, 8499.695, 1.55);

  vt1.SetLineColor(kRed);
  vt2.SetLineColor(kRed);
  vt3.SetLineColor(kRed);
  vt4.SetLineColor(kRed);
  vt5.SetLineColor(kRed);
  ms1.SetLineColor(kBlue);
  ms2.SetLineColor(kBlue);
  ms3.SetLineColor(kBlue);
  ms4.SetLineColor(kBlue);
  ms5.SetLineColor(kBlue);
  ms6.SetLineColor(kBlue);

  vt1.SetLineWidth(3);
  vt2.SetLineWidth(3);
  vt3.SetLineWidth(3);
  vt4.SetLineWidth(3);
  vt5.SetLineWidth(3);
  ms1.SetLineWidth(3);
  ms2.SetLineWidth(3);
  ms3.SetLineWidth(3);
  ms4.SetLineWidth(3);
  ms5.SetLineWidth(3);
  ms6.SetLineWidth(3);

  vt1.Draw("same");
  vt2.Draw("same");
  vt3.Draw("same");
  vt4.Draw("same");
  vt5.Draw("same");

  ms1.Draw("same");
  ms2.Draw("same");
  ms3.Draw("same");
  ms4.Draw("same");
  ms5.Draw("same");
  ms6.Draw("same");
  cc2->SaveAs("BFieldZ.png");
  fclose(fptxt);
  fclose(fptxtZRotated);
  TFile *f = new TFile("BFieldNA60plus_longsetup.root", "recreate");
  hField->Write();
  cc2->Write();
  f->Close();
}
