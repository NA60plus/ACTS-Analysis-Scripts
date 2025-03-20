import pandas as pd

# Read the CSV file
fakes = pd.read_csv("/home/giacomo/Downloads/data_fake.csv")  # Replace with your actual file path
match = pd.read_csv("/home/giacomo/Downloads/all_match.csv")  # Replace with your actual file path

# Sort data by the 'x' column
fakes = fakes.sort_values(by='x')
match = match.sort_values(by='x')
# Compute the sum of column 'y'
fakes_sum = fakes['y'].sum()
all_match = match['y'].sum()

print("Sum of column y:", fakes_sum)

print("Sum of column y:", all_match)

print("Sum of column y:", (all_match-fakes_sum)/fakes_sum)

import ROOT


hMatch = ROOT.TH1D("hMatch",";M (GeV/#it{c}^{2});countes",len(match['y']),0.4,1.6)
hFake = ROOT.TH1D("hFake",";M (GeV/#it{c}^{2});countes",len(match['y']),0.4,1.6)
index = 0
for y_fake,y_match in zip(fakes['y'],match["y"]):
    index += 1
    hFake.SetBinContent(index,y_fake)
    hMatch.SetBinContent(index,y_match)

cv = ROOT.TCanvas("cv","cv")

hMatch.SetMarkerStyle(24)
hFake.SetMarkerStyle(25)
hMatch.SetMarkerColor(ROOT.kBlack)
hFake.SetMarkerColor(ROOT.kBlue)
cv.SetLogy(1)
hMatch.Draw()
hFake.Draw("same")
cv.SaveAs("omegaOld.png")