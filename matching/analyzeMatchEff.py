import ROOT
import math
import array  # Import the array module

import pandas as pd



def ratio(x):
    return 1/x**2-1

def get_eff(x):
    return 1/math.sqrt(x+1)
#print(get_eff(0.26))

# Read the CSV file
df_chi2 = pd.read_csv("chi2NA60.csv")  # Replace with your actual file path
int_chi2_min_3 = df_chi2['y'].sum()

chi2max = 1.5
int_chi2_min_1_5 = 0
for x,y in zip(df_chi2['x'], df_chi2['y']):
    if x < 1.5:
        int_chi2_min_1_5 += y

print("Matching efficiency NA60 (chi2<3):", int_chi2_min_3)
print("Matching efficiency NA60 (chi2<1.5):", int_chi2_min_1_5)


file_list = [
    "omega_full_ruben/matching.root",
    "phi_full_ruben/matching.root",
    "jpsi_full_ruben/matching.root"
]

ROOT.gStyle.SetPaintTextFormat(".2f")

pdg_codes = [
   223,
   333,
   443,
   223,
   333,
   443
]
suffix_list =  [
    "omega_full_ruben",
    "phi_full_ruben",
    "jpsi_full_ruben"
]


def get_equivalent_sel(hist, matching_eff):
    hist.Scale(1./hist.GetEntries())
    eff = 0
    for i in range(1,hist.GetNbinsX()+1):
        eff += hist.GetBinContent(i)
        if eff > matching_eff:
            return hist.GetBinWidth(1)*i

def applyBDTstyle(effVsPBDT):

    effVsPBDT.SetMarkerStyle(20)
    effVsPBDT.SetMarkerSize(1.2)
    effVsPBDT.SetLineWidth(2)

    effVsPBDT.SetMarkerColor(ROOT.kBlue)
    effVsPBDT.SetLineColor(ROOT.kBlue)

def applyChi2style(effVsPChi2):
    effVsPChi2.SetMarkerStyle(20)
    effVsPChi2.SetMarkerSize(1.2)
    effVsPChi2.SetLineWidth(2)
    effVsPChi2.SetMarkerColor(ROOT.kRed)
    effVsPChi2.SetLineColor(ROOT.kRed)

ROOT.gROOT.SetBatch(True)

ROOT.gStyle.SetOptStat(0)
def main():
    #print(get_equivalent_sel(chi2Phi, int_chi2_min_1_5))

    # Example data points
    x = []  # X values
    y = []  # Y values
    ex = [0]*len(file_list)  # X errors
    ey = [0]*len(file_list)  # Y errors

    pdgDB = ROOT.TDatabasePDG.Instance()


    cv3 = ROOT.TCanvas("cv","cv",1920,1920)  
    rootfile_omega = ROOT.TFile("omega_full_ruben/matching.root","read")
    rootfile_phi = ROOT.TFile("phi_full_ruben/matching.root","read")
    rootfile_jpsi = ROOT.TFile("jpsi_full_ruben/matching.root","read")
    
    effVsElossChi2_omega = rootfile_omega.Get("effVsElossChi2")
    effVsElossChi2_phi = rootfile_phi.Get("effVsElossChi2")
    effVsElossChi2_jpsi = rootfile_jpsi.Get("effVsElossChi2")
    applyChi2style(effVsElossChi2_jpsi)
    applyChi2style(effVsElossChi2_phi)
    applyChi2style(effVsElossChi2_omega)

    effVsElossChi2_phi.SetMarkerColor(ROOT.kBlue)
    effVsElossChi2_phi.SetLineColor(ROOT.kBlue)

    effVsElossChi2_omega.SetMarkerColor(ROOT.kGreen+2)
    effVsElossChi2_omega.SetLineColor(ROOT.kGreen+2)

    effVsElossChi2_jpsi.SetTitle(";#Delta E (GeV/#it{c}^{2});Matching purity (%)")
    effVsElossChi2_jpsi.Draw()
    effVsElossChi2_phi.Draw("same")
    effVsElossChi2_omega.Draw("same")

    leg = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    leg.AddEntry(effVsElossChi2_jpsi, "#mu from J/#psi", "lep")
    leg.AddEntry(effVsElossChi2_phi, "#mu from #phi", "lep")
    leg.AddEntry(effVsElossChi2_omega, "#mu from #omega", "lep")
    leg.Draw()
    cv3.SaveAs("effVsElossChi2_comp.png")
    ###########################

    effVsElossRelChi2_omega = rootfile_omega.Get("effVsElossRelChi2")
    effVsElossRelChi2_phi = rootfile_phi.Get("effVsElossRelChi2")
    effVsElossRelChi2_jpsi = rootfile_jpsi.Get("effVsElossRelChi2")
    applyChi2style(effVsElossRelChi2_jpsi)
    applyChi2style(effVsElossRelChi2_phi)
    applyChi2style(effVsElossRelChi2_omega)

    effVsElossRelChi2_phi.SetMarkerColor(ROOT.kBlue)
    effVsElossRelChi2_phi.SetLineColor(ROOT.kBlue)

    effVsElossRelChi2_omega.SetMarkerColor(ROOT.kGreen+2)
    effVsElossRelChi2_omega.SetLineColor(ROOT.kGreen+2)

    effVsElossRelChi2_jpsi.SetTitle(";#Delta E /E;Matching purity (%)")
    effVsElossRelChi2_jpsi.Draw()
    effVsElossRelChi2_phi.Draw("same")
    effVsElossRelChi2_omega.Draw("same")
    leg.Draw()
    cv3.SaveAs("effVsElossRelChi2_comp.png")

    for pdg, file, suffix in zip(pdg_codes, file_list, suffix_list):
        particle = pdgDB.GetParticle(pdg)
        x.append(particle.Mass())

        rootfile = ROOT.TFile(file,"read")
        all = rootfile.Get("hPDistr").GetEntries()
        matched = rootfile.Get("hPMatchedChi2").GetEntries()
        y.append(ratio(matched/all))
        

        hBDTscoreMatched = rootfile.Get("hBDTscoreMatched")
        hBDTscoreFake = rootfile.Get("hBDTscoreFake")
        
        hChi2Matched = rootfile.Get("hChi2Matched")
        hChi2Fake = rootfile.Get("hChi2Fake")
        
        hBDTscoreMatched.Rebin(2)
        hBDTscoreFake.Rebin(2)
        hChi2Matched.Rebin(2)
        hChi2Fake.Rebin(2)

        hBDTscoreMatched.SetFillStyle(3004)
        hBDTscoreFake.SetFillStyle(3004)
        hChi2Matched.SetFillStyle(3004)
        hChi2Fake.SetFillStyle(3004)

        hBDTscoreMatched.GetYaxis().SetTitle("Counts (a.u.)")
        hBDTscoreFake.GetYaxis().SetTitle("Counts (a.u.)")
        hChi2Matched.GetYaxis().SetTitle("Counts (a.u.)")
        hChi2Fake.GetYaxis().SetTitle("Counts (a.u.)")

        hBDTscoreMatched.SetFillColorAlpha(ROOT.kRed, 0.3)
        hBDTscoreFake.SetFillColorAlpha(ROOT.kBlue, 0.3)
        hChi2Matched.SetFillColorAlpha(ROOT.kRed, 0.3)
        hChi2Fake.SetFillColorAlpha(ROOT.kBlue, 0.3)

        cv = ROOT.TCanvas("cv","cv",1920,1920)
        cv.SetLogy(1)
        cv.SetFrameFillColor(0)  # Ensure transparent background

        hBDTscoreMatched.Scale(1./hBDTscoreMatched.GetEntries())
        hBDTscoreFake.Scale(1./hBDTscoreFake.GetEntries())
        hChi2Matched.Scale(1./hChi2Matched.GetEntries())
        hChi2Fake.Scale(1./hChi2Fake.GetEntries())
        hBDTscoreMatched.GetYaxis().SetRangeUser(0.000001,1)
        hChi2Matched.GetYaxis().SetRangeUser(0.000001,1)
        hBDTscoreMatched.Draw("F HIST")
        hBDTscoreFake.Draw("SAME F HIST")
        cv.Update()
        cv.SaveAs(suffix+"/BDTscoreComp.png")
        hChi2Matched.Draw("F HIST")
        hChi2Fake.Draw("SAME F HIST")
        cv.Update()
        cv.SaveAs(suffix+"/Chi2scoreComp.png")

        cv2 = ROOT.TCanvas("cv2","cv2",1920,1920)
        hYVsPVT = rootfile.Get("hYVsPVT")
        hYVsPVT.Draw("colz")
        cv2.SaveAs(suffix+"/allYvsPVT.png")
        hYVsPAllBDT = rootfile.Get("hYVsPAllBDT")
        hYVsPAllBDT.Draw("colz")
        cv2.SaveAs(suffix+"/hYVsPAllBDT.png")


        effVsPBDT = rootfile.Get("effVsPBDT")
        effVsPBDT.SetTitle(";#it{p} (GeV/#it{c});Matching purity (%)")
        effVsPBDT.Draw()
        cv2.SaveAs(suffix+"/effVsPBDT.png")

        effVsPChi2 = rootfile.Get("effVsPChi2")
        effVsPChi2.SetTitle(";#it{p} (GeV/#it{c});Matching purity (%)")
        effVsPChi2.Draw()
        cv2.SaveAs(suffix+"/effVsPChi2.png")

        effVsYBDT = rootfile.Get("effVsYBDT")
        effVsYBDT.SetTitle(";#it{y};Matching purity (%)")
        effVsYBDT.Draw()
        cv2.SaveAs(suffix+"/effVsYBDT.png")

        effVsYChi2 = rootfile.Get("effVsYChi2")
        applyChi2style(effVsYChi2)
        effVsYChi2.SetTitle(";#it{y};Matching purity (%)")
        effVsYChi2.Draw()
        cv2.SaveAs(suffix+"/effVsYChi2.png")

        effVsYVsPBDT = rootfile.Get("effVsYVsPBDT")
        effVsYVsPBDT.Draw("colz text")
        cv2.SaveAs(suffix+"/effVsYVsPBDT.png")
        
        effVsYVsPChi2 = rootfile.Get("effVsYVsPChi2")
        effVsYVsPChi2.Draw("colz text")
        cv2.SaveAs(suffix+"/effVsYVsPChi2.png")
        

        effVsElossBDT = rootfile.Get("effVsElossBDT")
        applyChi2style(effVsElossBDT)
        effVsElossBDT.SetTitle(";#Delta E (GeV/#it{c}^{2});Matching purity (%)")
        effVsElossBDT.Draw()
        cv2.SaveAs(suffix+"/effVsElossBDT.png")

        effVsElossChi2 = rootfile.Get("effVsElossChi2")
        applyChi2style(effVsElossChi2)
        effVsElossChi2.SetTitle(";#Delta E (GeV/#it{c}^{2});Matching purity (%)")
        effVsElossChi2.Draw()
        cv2.SaveAs(suffix+"/effVsElossChi2.png")

        effVsElossRelBDT = rootfile.Get("effVsElossRelBDT")
        applyChi2style(effVsElossRelBDT)
        effVsElossRelBDT.SetTitle(";#Delta E / E;Matching purity (%)")
        effVsElossRelBDT.Draw()
        cv2.SaveAs(suffix+"/effVsElossRelBDT.png")

        effVsElossRelChi2 = rootfile.Get("effVsElossRelChi2")
        applyChi2style(effVsElossRelChi2)
        effVsElossRelChi2.SetTitle(";#Delta E / E;Matching purity (%)")
        effVsElossRelChi2.Draw()
        cv2.SaveAs(suffix+"/effVsElossRelChi2.png")






    # Create the TGraphErrors object
    graph = ROOT.TGraphErrors(len(x), array.array('d', x), array.array('d', y), 
                            array.array('d', ex), array.array('d', ey))

    # Set graph title and axis labels
    graph.SetTitle(";M (GeV#it{c}^{2});fake matches / correct matches")

    # Customize the appearance
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(3)
    graph.SetLineWidth(3)
    graph.SetMarkerColor(ROOT.kBlue)
    graph.SetLineColor(ROOT.kRed)

    # Draw the graph
    canvas = ROOT.TCanvas("c1", "TGraphErrors Example", 1920, 1920)
    graph.Draw("APL")  # "A" to draw axis, "P" for points with errors
    canvas.SaveAs("fakeOvCorrectMatchesVsMass.png")  # Save as image
    canvas.Draw()


if __name__ == "__main__":
    main()
