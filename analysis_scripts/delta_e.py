import ROOT
import numpy as np
import math
# Apri i file ROOT
file1 = ROOT.TFile.Open("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output__noAbs/tracksummary_ambims.root")
file2 = ROOT.TFile.Open("/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output__noAbs/tracksummary_merged.root") # MS0

# Prendi l'albero tracksummary
tree1 = file1.Get("tracksummary")
tree2 = file2.Get("tracksummary")

# Assumiamo che i due alberi abbiano lo stesso numero di entries
n_entries = min(tree1.GetEntries(), tree2.GetEntries())

hDelta = ROOT.TH1F("hDelta", "Differenza di energia tra due alberi", 100, -0.5, 0.5)
hE1 = ROOT.TH1F("hE1", "Energia", 100, 0, 100)
hE2 = ROOT.TH1F("hE2", "Energia", 100, 0, 100)

mmuon = 0.1056583745  # Massa del muone in GeV/c^2

for i in range(n_entries):
    tree1.GetEntry(i)
    tree2.GetEntry(i)

    # Leggi majorityParticleId da entrambi
    pid1 = tree1.majorityParticleId
    pid2 = tree2.majorityParticleId

    # Controlla se sono uguali
    if pid1 == pid2:
        # Prendi energia, qui uso t_p come esempio
        energy1 = [ ((1/p1)**2 + mmuon**2)**0.5 for p1 in tree1.eQOP_fit ]
        energy2 = [ ((1/p2)**2 + mmuon**2)**0.5 for p2 in tree2.eQOP_fit ]

        energy1_np = np.array(energy1)
        energy2_np = np.array(energy2)
        print(energy1_np, energy2_np)
        diff_energy = energy1_np - energy2_np
        for e in energy1_np:
            hE1.Fill(e)
        for e in energy2_np:
            hE2.Fill(e)
        for delta in diff_energy:
            hDelta.Fill(delta)

cv = ROOT.TCanvas("cv", "Differenza di energia", 1600, 1200)
hDelta.SetXTitle("#Delta E (GeV)")
hDelta.SetYTitle("Entries")
hDelta.SetLineColor(ROOT.kBlue)
hDelta.Draw()
cv.SaveAs("delta_energy_comparison.png")

hE1.SetXTitle("Energia (GeV)")
hE1.SetYTitle("Entries")
hE1.SetLineColor(ROOT.kRed)
hE1.Draw()
hE2.SetLineColor(ROOT.kBlue)
hE2.Draw("SAME")
cv.SaveAs("energy_distribution.png")