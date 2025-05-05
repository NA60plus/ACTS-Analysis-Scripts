import ROOT

# Open the ROOT file
root_file = ROOT.TFile("geoRubenCarbon/propagation-material.root", "READ")
# Open the ROOT file
#root_file = ROOT.TFile("/home/giacomo/Downloads/propagation-material_Ruben_VTMS.root", "READ")

# Access the TTree
tree = root_file.Get("material-tracks")

# Check if the TTree exists
if not tree:
    print("TTree not found!")
    root_file.Close()
    exit()

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".0f")
hLay = [ROOT.TH2F("hLay"+str(i), "", 30, -150, 150, 30, -150, 150) for i in range(11)]
hLayX0 = [ROOT.TH2F("hLayX0"+str(i), "", 30, -150, 150, 30, -150, 150) for i in range(11)]
hZ = [ROOT.TH1F("hZ"+str(i), "", 30,0.5,29.5) for i in range(11)]
hLayRef = [ROOT.TH2F("hLayRef"+str(i), "", 30, -150, 150, 30, -150, 150) for i in range(11)]
counter = 0
# Loop over the entries in the TTree
for entry in tree:
    for i in range(0, entry.mat_z.size()):
        k = 5
        if entry.mat_z[i] < 100 and entry.mat_z[i] > 60:
            k = 0
        hZ[k].Fill(entry.mat_Z[i])
        hLay[k].Fill(entry.mat_x[i], entry.mat_y[i], entry.mat_Z[i])
        hLayX0[k].Fill(entry.mat_x[i], entry.mat_y[i], entry.mat_X0[i]/1000)
        hLayRef[k].Fill(entry.mat_x[i], entry.mat_y[i])
    counter += 1
    if counter > 50000:
        break

for i in range(1,hLay[0].GetNbinsX()+1):
    for j in range(1, hLay[0].GetNbinsY()+1):
        for k in range(11):
            hLay[k].SetBinContent(i, j, hLay[k].GetBinContent(i, j) / hLayRef[k].GetBinContent(i, j) if hLayRef[k].GetBinContent(i, j) > 0 else 0)
            hLayX0[k].SetBinContent(i, j, hLayX0[k].GetBinContent(i, j) / hLayRef[k].GetBinContent(i, j) if hLayRef[k].GetBinContent(i, j) > 0 else 0)

cv = ROOT.TCanvas("cv", "cv", 800, 600)
for i in range(11):
    hLay[i].Draw("colz text")
    cv.SaveAs("geoRuben/propagation-material"+str(i)+".png")
    hLayX0[i].Draw("colz text")
    cv.SaveAs("geoRuben/propagation-material"+str(i)+"-x0.png")
    hZ[i].Draw()
    cv.SaveAs("geoRuben/propagation-material"+str(i)+"-z.png")

# Close the ROOT file after processing
root_file.Close()
