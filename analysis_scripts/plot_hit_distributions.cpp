TGraphErrors* DivideGraphs(const TGraphErrors* g1, const TGraphErrors* g2) {
    int n1 = g1->GetN();
    int n2 = g2->GetN();

    if (n1 != n2) {
        std::cerr << "Graphs have different number of points!" << std::endl;
        return nullptr;
    }

    auto g_ratio = new TGraphErrors();
    
    for (int i = 0; i < n1; ++i) {
        double x1, y1, x2, y2;
        g1->GetPoint(i, x1, y1);
        g2->GetPoint(i, x2, y2);
        
        if (x1 != x2) {
            std::cerr << "X values don't match at point " << i << "!" << std::endl;
            delete g_ratio;
            return nullptr;
        }

        double ey1 = g1->GetErrorY(i);
        double ey2 = g2->GetErrorY(i);

        if (y2 == 0) {
            std::cerr << "Division by zero at point " << i << "!" << std::endl;
            continue;
        }

        double ratio = y1 / y2;
        double eratio = ratio * std::sqrt( std::pow(ey1 / y1, 2) + std::pow(ey2 / y2, 2) );

        g_ratio->SetPoint(i, x1, ratio);
        g_ratio->SetPointError(i, 0, eratio);  // Assume no x error
    }

    return g_ratio;
}


void plot_hit_distributions()
{
  TFile *file_ex = TFile::Open("/home/giacomo/acts_for_NA60+/teltest/telescope_simulation/fatras/hits.root");    // Replace with your actual ROOT file
  TTree *tree_ex = (TTree *)file_ex->Get("hits");                    // Replace with your TTree name
  TFile *file_c = TFile::Open("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_checks_carbon/hits_vt.root"); // Replace with your actual ROOT file
  TTree *tree_c = (TTree *)file_c->Get("hits");                  // Replace with your TTree name

  const int nLayers = 5;
  double zlayer[nLayers] = {72,
                           152,
                           202,
                           252,
                           382,
                          };
  double layer[nLayers] = {0,1,2,3,4};

  double mean_x_ex[nLayers], rms_x_ex[nLayers], mean_y_ex[nLayers], rms_y_ex[nLayers];
  double mean_x_c[nLayers], rms_x_c[nLayers], mean_y_c[nLayers], rms_y_c[nLayers];
  double err_mean_x_ex[nLayers], err_rms_x_ex[nLayers], err_mean_y_ex[nLayers], err_rms_y_ex[nLayers];
  double err_mean_x_c[nLayers], err_rms_x_c[nLayers], err_mean_y_c[nLayers], err_rms_y_c[nLayers];

  double error_x[nLayers] = {0.5,0.5,0.5,0.5,0.5}; // x-errors (usually 0 for this type of plot)

  TCanvas *c1 = new TCanvas("c1", "Hit Distributions", 1200, 800);
  TFile *outputFile = new TFile("hit_distributions.root", "RECREATE");
  
  for (int i = 0; i < nLayers; ++i)
  {
    TString cut_ex = Form("tz < %f && tz > %f", zlayer[i], zlayer[i]-10);
    TString cut_c = Form("tz < %f && tz > %f", zlayer[i], zlayer[i]-10);

    tree_ex->Draw("0.1*tx>>hx(300,-0.15,0.15)", cut_ex, "goff");
    TH1 *hx = (TH1 *)gDirectory->Get("hx");
    mean_x_ex[i] = hx->GetMean();
    rms_x_ex[i] = hx->GetRMS();
    err_mean_x_ex[i] = hx->GetMeanError();
    err_rms_x_ex[i] = hx->GetRMSError();
    hx->Draw();
    c1->SaveAs(Form("hit_distribution_x_layer_%d_example.png", i));
    tree_c->Draw("0.1*tx>>hx2(300,1.85,2.15)", cut_c, "goff");
    TH1 *hx2 = (TH1 *)gDirectory->Get("hx2");
    mean_x_c[i] = hx2->GetMean();
    rms_x_c[i] = hx2->GetRMS();
    err_mean_x_c[i] = hx2->GetMeanError();
    err_rms_x_c[i] = hx2->GetRMSError();
    hx2->Draw();
    c1->SaveAs(Form("hit_distribution_x_layer_%d_carbon.png", i));

    tree_ex->Draw("0.1*ty>>hy(300,-0.15,0.15)", cut_ex, "goff");
    TH1 *hy = (TH1 *)gDirectory->Get("hy");
    mean_y_ex[i] = hy->GetMean();
    rms_y_ex[i] = hy->GetRMS();
    err_mean_y_ex[i] = hy->GetMeanError();
    err_rms_y_ex[i] = hy->GetRMSError();
    hy->Draw();
    c1->SaveAs(Form("hit_distribution_y_layer_%d_example.png", i));

    tree_c->Draw("0.1*ty>>hy2(300,1.85,2.15)", cut_c, "goff");
    TH1 *hy2 = (TH1 *)gDirectory->Get("hy2");
    mean_y_c[i] = hy2->GetMean();
    rms_y_c[i] = hy2->GetRMS();
    err_mean_y_c[i] = hy2->GetMeanError();
    err_rms_y_c[i] = hy2->GetRMSError();
    hy2->Draw();
    c1->SaveAs(Form("hit_distribution_y_layer_%d_carbon.png", i));

    hx2->SetName(Form("hx2_layer%d", i));
    hy2->SetName(Form("hy2_layer%d", i));
    hx2->Write();
    hy2->Write();
    
  }


  TCanvas *c = new TCanvas("c", "Hit Distributions", 1200, 800);
  c->Divide(2, 2);

  TGraphErrors *gr_mean_x_ex = new TGraphErrors(nLayers, layer, mean_x_ex, error_x, err_mean_x_ex);
  TGraphErrors *gr_mean_x_c = new TGraphErrors(nLayers, layer, mean_x_c, error_x, err_mean_x_c);

  TGraphErrors *gr_rms_x_ex = new TGraphErrors(nLayers, layer, rms_x_ex, error_x, err_rms_x_ex);
  TGraphErrors *gr_rms_x_c = new TGraphErrors(nLayers, layer, rms_x_c, error_x, err_rms_x_c);

  TGraphErrors *gr_mean_y_ex = new TGraphErrors(nLayers, layer, mean_y_ex, error_x, err_mean_y_ex);
  TGraphErrors *gr_mean_y_c = new TGraphErrors(nLayers, layer, mean_y_c, error_x, err_mean_y_c);

  TGraphErrors *gr_rms_y_ex = new TGraphErrors(nLayers, layer, rms_y_ex, error_x, err_rms_y_ex);
  TGraphErrors *gr_rms_y_c = new TGraphErrors(nLayers, layer, rms_y_c, error_x, err_rms_y_c);

  // Style
  gr_mean_x_ex->SetLineColor(kBlack);
  gr_mean_x_c->SetLineColor(kRed);
  gr_rms_x_ex->SetLineColor(kBlack);
  gr_rms_x_c->SetLineColor(kRed);
  gr_mean_y_ex->SetLineColor(kBlack);
  gr_mean_y_c->SetLineColor(kRed);
  gr_rms_y_ex->SetLineColor(kBlack);
  gr_rms_y_c->SetLineColor(kRed);

  // Style
  gr_mean_x_ex->SetLineWidth(2);
  gr_mean_x_c->SetLineWidth(2);
  gr_rms_x_ex->SetLineWidth(2);
  gr_rms_x_c->SetLineWidth(2);
  gr_mean_y_ex->SetLineWidth(2);
  gr_mean_y_c->SetLineWidth(2);
  gr_rms_y_ex->SetLineWidth(2);
  gr_rms_y_c->SetLineWidth(2);

  gr_rms_x_ex->SetName("rms_x_example");
  gr_rms_x_c->SetName("rms_x_converted");
  gr_rms_y_ex->SetName("rms_y_example");
  gr_rms_y_c->SetName("rms_y_converted");

  gr_rms_x_ex->Write();
  gr_rms_x_c->Write();
  gr_rms_y_ex->Write();
  gr_rms_y_c->Write();

  outputFile->Close();


  gr_mean_x_ex->GetXaxis()->SetRangeUser(-0.5, 4.5);
  gr_mean_y_ex->GetXaxis()->SetRangeUser(-0.5, 4.5);
  gr_rms_x_ex->GetXaxis()->SetRangeUser(-0.5, 4.5);
  gr_rms_y_ex->GetXaxis()->SetRangeUser(-0.5, 4.5);
  // Plot 1: Mean x
  c->cd(1);
  gr_mean_x_ex->Draw("AP");
  gr_mean_x_c->Draw("P SAME");
  gr_mean_x_ex->SetTitle("<x> (cm);Layer;<x> (cm)");
  auto leg1 = new TLegend(0.1, 0.75, 0.5, 0.9);
  leg1->AddEntry(gr_mean_x_ex, "simple telescope", "l");
  leg1->AddEntry(gr_mean_x_c, "my mapping", "l");
  leg1->Draw();

  // Plot 2: RMS x
  c->cd(2);
  gr_rms_x_ex->Draw("AP");
  gr_rms_x_c->Draw("P SAME");
  gr_rms_x_ex->SetTitle("rms_{x} (cm);Layer;rms_{x} (cm)");

  // Plot 3: Mean y
  c->cd(3);
  gr_mean_y_ex->Draw("AP");
  gr_mean_y_c->Draw("P SAME");
  gr_mean_y_ex->SetTitle("<y> (cm);Layer;<y> (cm)");

  // Plot 4: RMS y
  c->cd(4);
  gr_rms_y_ex->Draw("AP");
  gr_rms_y_c->Draw("P SAME");
  gr_rms_y_ex->SetTitle("rms_{y} (cm);Layer;rms_{y} (cm)");

  c->SaveAs("HitDistributions.png");

  TGraphErrors* ratio_rms_y = DivideGraphs(gr_rms_x_ex, gr_rms_x_c);
  TGraphErrors* ratio_rms_x = DivideGraphs(gr_rms_y_ex, gr_rms_y_c);

  TCanvas *c_ratio = new TCanvas("c_ratio", "Ratio of RMS Distributions", 800, 600);
  ratio_rms_y->SetTitle("Ratio of RMS with/without Carbon Plates;Layer;Ratio");
  ratio_rms_y->SetLineColor(kBlue);
  ratio_rms_y->SetLineWidth(2);
  ratio_rms_y->Draw("AP");
  ratio_rms_x->SetLineColor(kRed);
  ratio_rms_x->SetLineWidth(2);
  ratio_rms_x->Draw("P SAME");
  ratio_rms_x->GetXaxis()->SetRangeUser(-0.5, 4.5);
  float ratio_x0 = TMath::Sqrt((0.00178936+0.000533842)/0.000533842);

  std::cout << "Ratio x0: " << ratio_x0 << std::endl;
  TLine *line = new TLine(-0.5, ratio_x0, 4.5, ratio_x0);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->Draw("SAME");
  TLatex *latex = new TLatex(2, ratio_x0 - 0.35, Form("#frac{#theta_{rms}^{w}}{#theta_{rms}^{c}} = #sqrt{#frac{(x/X0)_{w}}{(x/X0)_{c}}} = %.2f", ratio_x0));
  latex->SetTextSize(0.04);

  auto legr = new TLegend(0.5, 0.25, 0.8, 0.65);
  legr->AddEntry(ratio_rms_x, "x", "pl");
  legr->AddEntry(ratio_rms_y, "y", "pl");
  legr->AddEntry(line, Form("#frac{#theta_{rms}^{w}}{#theta_{rms}^{c}} = #sqrt{#frac{(x/X0)_{w}}{(x/X0)_{c}}} = %.2f", ratio_x0), "l");
  legr->Draw();


  //latex->Draw("SAME");

  c_ratio->SaveAs("RatioRMSY.png");
  
}
