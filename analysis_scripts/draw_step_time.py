import os
import ROOT
import csv

ROOT.gStyle.SetTitleSize(0.05, "XYZ") # Sets the size for X, Y, and Z axis titles
ROOT.gStyle.SetTitleFontSize(0.05) # Sets the size for the histogram title
ROOT.gStyle.SetLabelSize(0.04, "XYZ")
ROOT.gStyle.SetPadLeftMargin(0.15);   #Set left margin
ROOT.gStyle.SetPadRightMargin(0.05);  #Set right margin
ROOT.gStyle.SetPadTopMargin(0.05);    #Set top margin
ROOT.gStyle.SetPadBottomMargin(0.15); #Set bottom margin

def get_time_step(dir):
    keeplist = [
            "Algorithm:SeedingAlgorithmNA60",
            "Algorithm:TrackParamsEstimationAlgorithm",
            "Algorithm:SeedsToPrototracks",
            "Algorithm:TrackFindingAlgorithm",
            "Algorithm:TracksToTrajectories",
            "Algorithm:GreedyAmbiguityResolutionAlgorithm"
            "Algorithm:FilterMeasurementsAlgorithm",
            "Algorithm:SpacePointMaker",
            #"Writer:RootTrackParameterWriter"
        ]
    
    step_counter = -1
    step_times = []
    with open(dir+'/timing.tsv', newline='') as tsvfile:
        # Create a CSV reader object with tab delimiter
        reader = csv.reader(tsvfile, delimiter='\t')
        # Iterate over each row in the file
        for row in reader:
            if row[0] in keeplist:
                if "Algorithm:SpacePointMaker" in row[0]:
                    step_counter += 1
                    step_times.append(float(row[2]))
                step_times[step_counter] += float(row[2])
        return step_times
    

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
    #list of names for the plots
    suffix_list = [
                    "0.1",
                    "0.2",
                    "0.4"
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

    hs = ROOT.THStack("hs",";Track multiplicity; CPU time/event (a.u.) ")
    hsnorm = ROOT.THStack("hsnorm",";Track multiplicity; CPU time /Max CPU time")
    hsself = ROOT.THStack("hsself",";Track multiplicity; fraction of CPU time")
    h1 = ROOT.TH1D("h1",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    h2 = ROOT.TH1D("h2",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    h3 = ROOT.TH1D("h3",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    h4 = ROOT.TH1D("h4",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)

    h1_norm = ROOT.TH1D("h1_norm",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    h2_norm = ROOT.TH1D("h2_norm",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    h3_norm = ROOT.TH1D("h3_norm",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    h4_norm = ROOT.TH1D("h4_norm",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)

    h1_self = ROOT.TH1D("h1_self",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    h2_self = ROOT.TH1D("h2_self",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    h3_self = ROOT.TH1D("h3_self",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    h4_self = ROOT.TH1D("h4_self",";;CPU Time / event",len(mult_list), 0*mult_list[-1]*ntrk, mult_list[-1]*ntrk)
    
    hsself.SetMinimum(0)
    hsself.SetMaximum(1)    
    
    hsnorm.SetMinimum(0)
    hsnorm.SetMaximum(1)    

    legend = ROOT.TLegend(0.15,0.6,0.3,0.8)
    legend.AddEntry(h1, "Step 1","f")
    legend.AddEntry(h2, "Step 2","f")
    legend.AddEntry(h3, "Step 3","f")
    legend.AddEntry(h4, "Step 4","f")

    legend2 = ROOT.TLegend(0.25,0.3,0.4,0.5)
    legend2.AddEntry(h1, "Step 1","f")
    legend2.AddEntry(h2, "Step 2","f")
    legend2.AddEntry(h3, "Step 3","f")
    legend2.AddEntry(h4, "Step 4","f")

    for i, (dir) in enumerate(dir_list):
        tot_time = get_time_step(dir=dir)
        h1.SetBinContent(i+1, tot_time[0])
        h2.SetBinContent(i+1, tot_time[1])
        h3.SetBinContent(i+1, tot_time[2])
        h4.SetBinContent(i+1, tot_time[3])

        sum_tot_time = tot_time[0] + tot_time[1] + tot_time[2] + tot_time[3]
        h1_self.SetBinContent(i+1, tot_time[0]/sum_tot_time)
        h2_self.SetBinContent(i+1, tot_time[1]/sum_tot_time)
        h3_self.SetBinContent(i+1, tot_time[2]/sum_tot_time)
        h4_self.SetBinContent(i+1, tot_time[3]/sum_tot_time)

        h1_norm.SetBinContent(i+1, tot_time[0])
        h2_norm.SetBinContent(i+1, tot_time[1])
        h3_norm.SetBinContent(i+1, tot_time[2])
        h4_norm.SetBinContent(i+1, tot_time[3])

    colors = [
        ROOT.kViolet + 9,
        ROOT.kYellow - 7,
        ROOT.kRed - 10,
        ROOT.kSpring + 9
    ]
    colors = [
        ROOT.kPink - 9,
        ROOT.kAzure - 9,
        ROOT.kTeal - 7,
        ROOT.kOrange - 3
    ]
    
    colors = [
        ROOT.kOrange + 2,
        ROOT.kRed -3 ,
        ROOT.kGreen -2,
        ROOT.kAzure - 2
    ]
    h1.SetFillColor(colors[0])
    h1.SetMarkerColor(colors[0])
    h2.SetFillColor(colors[1])
    h2.SetMarkerColor(colors[1])
    h3.SetFillColor(colors[2])
    h3.SetMarkerColor(colors[2])
    h4.SetFillColor(colors[3])
    h4.SetMarkerColor(colors[3])

    h1.SetMarkerStyle(21)
    h2.SetMarkerStyle(21)
    h3.SetMarkerStyle(21)
    h4.SetMarkerStyle(21)
    hs.Add(h1)
    hs.Add(h2)
    hs.Add(h3)
    hs.Add(h4)


    h1_self.SetFillColor(colors[0])
    h1_self.SetMarkerColor(colors[0])
    h2_self.SetFillColor(colors[1])
    h2_self.SetMarkerColor(colors[1])
    h3_self.SetFillColor(colors[2])
    h3_self.SetMarkerColor(colors[2])
    h4_self.SetFillColor(colors[3])
    h4_self.SetMarkerColor(colors[3])

    h1_self.SetMarkerStyle(21)
    h2_self.SetMarkerStyle(21)
    h3_self.SetMarkerStyle(21)
    h4_self.SetMarkerStyle(21)

    
    h1_norm.SetFillColor(colors[0])
    h1_norm.SetMarkerColor(colors[0])
    h2_norm.SetFillColor(colors[1])
    h2_norm.SetMarkerColor(colors[1])
    h3_norm.SetFillColor(colors[2])
    h3_norm.SetMarkerColor(colors[2])
    h4_norm.SetFillColor(colors[3])
    h4_norm.SetMarkerColor(colors[3])

    h1_norm.SetMarkerStyle(21)
    h2_norm.SetMarkerStyle(21)
    h3_norm.SetMarkerStyle(21)
    h4_norm.SetMarkerStyle(21)

    
    hsself.Add(h1_self)
    hsself.Add(h2_self)
    hsself.Add(h3_self)
    hsself.Add(h4_self)

    cv = ROOT.TCanvas("cv","cv",1800,1200)
    hs.Draw()
    legend.Draw()
    cv.SaveAs("cpu_time.png")

    max_time = h1.GetBinContent(len(mult_list)) + h2.GetBinContent(len(mult_list)) + h3.GetBinContent(len(mult_list)) + h4.GetBinContent(len(mult_list))
    print(max_time)
    h1_norm.Scale(1./max_time)
    h2_norm.Scale(1./max_time)
    h3_norm.Scale(1./max_time)
    h4_norm.Scale(1./max_time)

    
    hsnorm.Add(h1_norm)
    hsnorm.Add(h2_norm)
    hsnorm.Add(h3_norm)
    hsnorm.Add(h4_norm)

    hsnorm.Draw("hist")
    legend.Draw()
    cv.SaveAs("cpu_time_norm.png")

    hsself.Draw("hist")
    legend2.Draw()
    cv.SaveAs("self_cpu_time.png")


if __name__ == "__main__":
    main()