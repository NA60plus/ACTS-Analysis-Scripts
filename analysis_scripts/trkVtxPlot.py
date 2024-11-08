import os
import ROOT

ROOT.gStyle.SetTitleSize(0.05, "XYZ") # Sets the size for X, Y, and Z axis titles
ROOT.gStyle.SetTitleFontSize(0.05) # Sets the size for the histogram title
ROOT.gStyle.SetLabelSize(0.04, "XYZ")
ROOT.gStyle.SetPadLeftMargin(0.15)   #Set left margin
ROOT.gStyle.SetPadRightMargin(0.05)  #Set right margin
ROOT.gStyle.SetPadTopMargin(0.05)    #Set top margin
ROOT.gStyle.SetPadBottomMargin(0.15) #Set bottom margin
ROOT.gStyle.SetHistLineWidth(2)

def main():
    #directories containing timing.tsv
    dir_list = [
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.05_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.1_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.15_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.2_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.25_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.3_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.35_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.4_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.45_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.5_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.55_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.6_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.65_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.7_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.75_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.8_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.85_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.9_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_periferal_factor_0.95_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix",
            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix"
        ]
    
    #fraction of multiplicity in central collisions
    mult_list = [
                    0.05,
                    0.1,
                    0.15,
                    0.2,
                    0.25,
                    0.3,
                    0.35,
                    0.4,
                    0.45,
                    0.5,
                    0.55,
                    0.6,
                    0.65,
                    0.7,
                    0.75,
                    0.8,
                    0.85,
                    0.9,
                    0.95,
                    1
                ]


    # Directory name
    directory = "../results/cpu_time_studies"
    # Create the directory
    os.makedirs(directory, exist_ok=True)

    ntrk=685

    hpassed = ROOT.TH1D("hpassed",";Track multiplicity:Target ID efficiency",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    htotal = ROOT.TH1D("htotal",";Track multiplicity;Target ID efficiency",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    
    for i, (dir) in enumerate(dir_list):
        file = ROOT.TFile(dir+"/performance_tracklet_vertexing.root","read")
        eff_vs_zgen = file.Get("eff_vs_zgen")
 
        passed = eff_vs_zgen.GetPassedHistogram().GetEntries()
        total = eff_vs_zgen.GetTotalHistogram().GetEntries()

        hpassed.SetBinContent(i+1, passed)
        htotal.SetBinContent(i+1, total)

    eff = ROOT.TEfficiency(hpassed,htotal)
    eff.SetMarkerStyle(20)
    cv = ROOT.TCanvas("cv","cv",1800,1200)
    eff.Draw()
    cv.SaveAs("trkVtxEffciency.png")
    ROOT.gStyle.SetOptStat(0)
    file = ROOT.TFile(dir_list[-1]+"/performance_tracklet_vertexing.root","read")
    allTrk = file.Get("all/histName4")
    allTrk.SetTitle(";#it{z} (cm);Entries")
    trueTrk = file.Get("peak/histPeakName4")
    trueTrk.SetLineColor(ROOT.kRed)


    allTrk.SetTitleSize(0.05, "XYZ") # Sets the size for X, Y, and Z axis titles
    #allTrk.SetTitleFontSize(0.05) # Sets the size for the histogram title
    allTrk.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetPadLeftMargin(0.15)   #Set left margin
    ROOT.gStyle.SetPadRightMargin(0.05)  #Set right margin
    ROOT.gStyle.SetPadTopMargin(0.05)    #Set top margin
    ROOT.gStyle.SetPadBottomMargin(0.15) #Set bottom margin
    allTrk.SetLineWidth(2)
    trueTrk.SetLineWidth(2)

    legend = ROOT.TLegend(0.2,0.65,0.5,0.85)
    legend.AddEntry(allTrk, "All tracklets","f")
    legend.AddEntry(trueTrk, "True tracklets","f")

    maxRangeY = allTrk.GetMaximum()
    allTrk.GetYaxis().SetRangeUser(0,maxRangeY*1.1)
    allTrk.Draw()
    trueTrk.Draw("same")
    legend.Draw()
    cv.SaveAs("plotTrkZ.png")    

    # Define a rectangle using TBox(x1, y1, x2, y2)
    # (x1, y1) is the bottom-left corner
    # (x2, y2) is the top-right corner

    box1 = ROOT.TBox(-1.5, 0, 0, maxRangeY*1.05)
    box2 = ROOT.TBox(-1.5*2-12, 0, (-1.5-12), maxRangeY*1.05)
    box3 = ROOT.TBox(-1.5*3-12*2, 0, (-1.5-12)*2, maxRangeY*1.05)
    box4 = ROOT.TBox(-1.5*4-12*3, 0, (-1.5-12)*3, maxRangeY*1.05)
    box5 = ROOT.TBox(-1.5*5-12*4, 0, (-1.5-12)*4, maxRangeY*1.05)

    # Set fill color and style
    box1.SetFillColor(ROOT.kRed)
    box1.SetFillStyle(3004)
    box2.SetFillColor(ROOT.kRed)
    box2.SetFillStyle(3004)
    box3.SetFillColor(ROOT.kRed)
    box3.SetFillStyle(3004)
    box4.SetFillColor(ROOT.kRed)
    box4.SetFillStyle(3004)
    box5.SetFillColor(ROOT.kRed)
    box5.SetFillStyle(3004)

    legend.AddEntry(box1, "Targets","f")
    # Draw the box1
    box1.Draw("same")
    box2.Draw("same")
    box3.Draw("same")
    box4.Draw("same")
    box5.Draw("same")
    legend.Draw()

    cv.SaveAs("plotTrkZTargets.png")    
    
if __name__ == "__main__":
    main()