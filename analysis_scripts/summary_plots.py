
import ROOT
import uproot
import os
from cpu_time_studies import get_run_time

def read_outfile(file):
    # Open the .out file for reading
    with open(file, 'r') as file:
        # Initialize a dictionary to store the values
        values = {}

        # Read each line in the file
        for line in file:
            # Check if the line starts with the desired prefix
            if line.startswith('NA60+_Summary_'):
                # Split the line by whitespace
                parts = line.split()

                # Check if the line contains the desired values
                if len(parts) == 2 and parts[0].endswith('='):
                    # Extract the key and value from the line
                    key = parts[0].replace('NA60+_Summary_', "").replace('=', "")
                    value = float(parts[1])

                    # Store the key-value pair in the dictionary
                    if key in values:
                        values[key] += value
                    else:
                        values[key] = value

    return values

    

def get_purity(file_path = "output_standardSeeding_iterativeVertexing/tracksummary_ambi.root"):
    tree_name = "tracksummary"

    # Open the ROOT file and get the TTree
    root_file = uproot.open(file_path)
    tree = root_file[tree_name]

    # Get the branch
    hit_tree = tree["nMajorityHits"]
    meas_tree = tree["nMeasurements"]

    # Get the data from the branch as a numpy array
    hit_data = hit_tree.array(library="np")
    meas_data = meas_tree.array(library="np")

    hMajOverMeas5 = ROOT.TH1F("hMajOverMeas5", ";nMajority/nHits (nHits=5);entries",5, 0.1, 1.1)
    hMajOverMeas4 = ROOT.TH1F("hMajOverMeas4", ";nMajority/nHits (nHits=4);entries",5, 0.1, 1.1)
    hMajOverMeas = ROOT.TH1F("hMajOverMeas", ";nMajority/nHits;entries",5, 0.1, 1.1)

    for hit_ev,meas_ev in zip(hit_data,meas_data):
        for hit,meas in zip(hit_ev,meas_ev):
            hMajOverMeas.Fill(hit/meas)
            if hit ==5:
                hMajOverMeas5.Fill(hit/meas)
            elif hit==4:
                hMajOverMeas4.Fill(hit/meas)
    
    return hMajOverMeas5, hMajOverMeas4, hMajOverMeas



# Set batch mode
ROOT.gROOT.SetBatch(True)


def drawSummary(files, dirlist, cutlist, outputDir):

    ROOT.gStyle.SetPaintTextFormat("0.0f")
    NbinsY = 5
    NbinsY_tot = 5
    th2_summary = ROOT.TH2D("th2_summary","",len(files),0,len(files),NbinsY,0,NbinsY)
    th2_summary_eff = ROOT.TH2D("th2_summary_eff","",len(files),0,len(files),NbinsY_tot,0,NbinsY_tot)
    th2_eff = ROOT.TH2D("th2_eff","",len(files),0,len(files),9,0,9)

    th1_eff_seed = ROOT.TH1D("th1_eff_seed","Seeding;;Efficiency",len(files),0,len(files))
    th1_eff_ckf = ROOT.TH1D("th1_eff_ckf","CKF;;Efficiency",len(files),0,len(files))
    th1_eff_ambi = ROOT.TH1D("th1_eff_ambi","Ambiguity resolution;;Efficiency",len(files),0,len(files))

    th1_fake_seed = ROOT.TH1D("th1_fake_seed","Seeding;;Fake rate",len(files),0,len(files))
    th1_fake_ckf = ROOT.TH1D("th1_fake_ckf","CKF;;Fake rate",len(files),0,len(files))
    th1_fake_ambi = ROOT.TH1D("th1_fake_ambi","Ambiguity resolution;;Fake rate",len(files),0,len(files))

    th1_dup_seed = ROOT.TH1D("th1_dup_seed","Seeding;;Duplicate rate",len(files),0,len(files))
    th1_dup_ckf = ROOT.TH1D("th1_dup_ckf","CKF;;Duplicate rate",len(files),0,len(files))
    th1_dup_ambi = ROOT.TH1D("th1_dup_ambi","Ambiguity resolution;;Duplicate rate",len(files),0,len(files))

    th1_seedPerPart = ROOT.TH1D("th1_seedPerPart",";; Seeds per particle",len(files),0,len(files))

    th1goodHitRate = ROOT.TH1D("th1goodHitRate",";; Fraction of good hits",len(files),0,len(files))

    th2_eff.GetYaxis().SetBinLabel(9, "Seeding efficiency")
    th2_eff.GetYaxis().SetBinLabel(8, "CKF efficiency")
    th2_eff.GetYaxis().SetBinLabel(7, "Ambi res. efficiency")
    th2_eff.GetYaxis().SetBinLabel(6, "Seeding fake rate")
    th2_eff.GetYaxis().SetBinLabel(5, "CKF fake rate")
    th2_eff.GetYaxis().SetBinLabel(4, "Ambi res. fake rate")
    th2_eff.GetYaxis().SetBinLabel(3, "Seeding duplicate rate")
    th2_eff.GetYaxis().SetBinLabel(2, "CKF duplicate rate")
    th2_eff.GetYaxis().SetBinLabel(1, "Ambi res. duplicate rate")

    # Check if the directory already exists
    outputDir = "results/" + outputDir
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        
    ROOT.gStyle.SetOptStat(0)
    output = ROOT.TFile(outputDir+"/summary.root","recreate")
    cv = ROOT.TCanvas("cv","cv",2200,1600)
    i = 0
    for file,dir in zip(files, dirlist):
        a = read_outfile(file)

        i+=1

        # Set the tick labels on the x-axis
        for k in range(1, th2_summary.GetNbinsX() + 1):
            th2_summary.GetXaxis().SetBinLabel(k, cutlist[k-1])
            th2_summary_eff.GetXaxis().SetBinLabel(k, cutlist[k-1])
            th2_eff.GetXaxis().SetBinLabel(k, cutlist[k-1])

            th1_eff_seed.GetXaxis().SetBinLabel(k, cutlist[k-1])
            th1_eff_ckf.GetXaxis().SetBinLabel(k, cutlist[k-1])
            th1_eff_ambi.GetXaxis().SetBinLabel(k, cutlist[k-1])

            th1_fake_seed.GetXaxis().SetBinLabel(k, cutlist[k-1])
            th1_fake_ckf.GetXaxis().SetBinLabel(k, cutlist[k-1])
            th1_fake_ambi.GetXaxis().SetBinLabel(k, cutlist[k-1])

            th1_dup_seed.GetXaxis().SetBinLabel(k, cutlist[k-1])
            th1_dup_ckf.GetXaxis().SetBinLabel(k, cutlist[k-1])
            th1_dup_ambi.GetXaxis().SetBinLabel(k, cutlist[k-1])

            th1_seedPerPart.GetXaxis().SetBinLabel(k, cutlist[k-1])

            th1goodHitRate.GetXaxis().SetBinLabel(k, cutlist[k-1])

        #########################
        #
        #  Set axis labels
        #
        ########################
    
        if i==1:    
            keys = a.keys()
            print(keys)
            # Set the tick labels on the y-axis
            j=1
            jeff=1
            for key in a.keys():
        
                if key == 'Eff':
                    th2_summary_eff.GetYaxis().SetBinLabel(NbinsY_tot+1-jeff, "Seeding efficiency")
                    jeff+=1
                elif key == 'Fakerate':
                    th2_summary_eff.GetYaxis().SetBinLabel(NbinsY_tot+1-jeff, "Seeding fake rate")
                    jeff+=1
                elif key == 'Duplication':
                    th2_summary_eff.GetYaxis().SetBinLabel(NbinsY_tot+1-jeff, "Seeding duplicate rate")
                    jeff+=1
                elif key =='nDuplicatedSeeds':
                    th2_summary_eff.GetYaxis().SetBinLabel(NbinsY_tot+1-jeff, key)
                    jeff+=1
                else:
                    th2_summary.GetYaxis().SetBinLabel(NbinsY+1-j, key)
                    j+=1
            th2_summary_eff.GetYaxis().SetBinLabel(NbinsY_tot+1-jeff, "Run time [s/event]")


        ##################################
        j = 1
        jeff = 1
        for key, value in a.items():
            if key != 'Eff' and key != 'Fakerate' and  key != 'Duplication' and key !='nDuplicatedSeeds':
                th2_summary.SetBinContent(i,NbinsY+1-j,value)
                j+=1
            else:
                th2_summary_eff.SetBinContent(i,NbinsY+1-jeff,value)
                jeff+=1

        eff_ckf = 0.0001
        dup_ckf = 0.0001
        fake_ckf = 0.0001
        eff_ambi = 0.0001
        dup_ambi = 0.0001
        fake_ambi = 0.0001
        step = ""
        while True:
            
            if os.path.exists(dir+"/performance_"+step+"ckf.root"):
                ckf_output = ROOT.TFile(dir+"/performance_"+step+"ckf.root","read")
                ambi_output = ROOT.TFile(dir+"/performance_ambi"+step+".root","read")
                eff_ckf = ckf_output.Get("eff_particles")[0]
                dup_ckf = ckf_output.Get("duplicaterate_particles")[0]
                fake_ckf = ckf_output.Get("fakerate_particles")[0]
                eff_ambi = ambi_output.Get("eff_particles")[0]
                dup_ambi = ambi_output.Get("duplicaterate_particles")[0]
                fake_ambi = ambi_output.Get("fakerate_particles")[0]
                step += "ll"
                ckf_output.Close()
                ambi_output.Close()
            else:
                break

        time,_,_ = get_run_time(dir)
        th2_summary_eff.SetBinContent(i,NbinsY_tot+1-jeff, time)
        th2_eff.SetBinContent(i,9,a["Eff"])   
        th2_eff.SetBinContent(i,8,eff_ckf)
        th2_eff.SetBinContent(i,7,eff_ambi)
        th2_eff.SetBinContent(i,6,a["Fakerate"])   
        th2_eff.SetBinContent(i,5,fake_ckf)
        th2_eff.SetBinContent(i,4,fake_ambi)
        th2_eff.SetBinContent(i,3,a["Duplication"])   
        th2_eff.SetBinContent(i,2,dup_ckf)
        th2_eff.SetBinContent(i,1,dup_ambi)

        th1_seedPerPart.SetBinContent(i, a["nTotalSeeds"]/a["nTotalParticles"])

        th1_eff_seed.SetBinContent(i, a["Eff"])
        th1_eff_ckf.SetBinContent(i, eff_ckf)
        th1_eff_ambi.SetBinContent(i, eff_ambi)

        th1_fake_seed.SetBinContent(i, a["Fakerate"])
        th1_fake_ckf.SetBinContent(i, fake_ckf)
        th1_fake_ambi.SetBinContent(i, fake_ambi)

        th1_dup_seed.SetBinContent(i, a["Duplication"])
        th1_dup_ckf.SetBinContent(i, dup_ckf)
        th1_dup_ambi.SetBinContent(i, dup_ambi)
        
        cv.SetLogy(1)
        hPurity5, hPurity4, hPurity = get_purity(file_path = dir+"/tracksummary_ambi.root")
        output.cd()
        hPurity5.Draw("text hist")
        cv.SaveAs(outputDir+"/hPurity5_"+file+".png")
        hPurity5.Write()

        hPurity4.Draw("text hist")
        cv.SaveAs(outputDir+"/hPurity4_"+file+".png")
        hPurity4.Write()

        hPurity.Draw("text hist")
        cv.SaveAs(outputDir+"/hPurity_"+file+".png")
        hPurity.Write()
        print(hPurity.GetMean())
        th1goodHitRate.SetBinContent(i, hPurity.GetMean())
    output.cd()
    cv.SetLogy(0)
    #ROOT.gStyle.SetPaintTextFormat("0.3f")

    th2_summary.Draw("text")

    ROOT.gPad.SetTopMargin(0.025)
    ROOT.gPad.SetBottomMargin(0.2)
    ROOT.gPad.SetLeftMargin(0.2)
    ROOT.gPad.SetRightMargin(0.1)
    cv.SaveAs(outputDir+"/seedingSummary.png")
    th2_summary.Write()

    ROOT.gStyle.SetPaintTextFormat("0.3f")
    th2_summary_eff.Draw("text")
    #ROOT.gPad.SetTopMargin(0.05)
    #ROOT.gPad.SetLeftMargin(0.2)
    #ROOT.gPad.SetRightMargin(0.04)
    cv.SaveAs(outputDir+"/seedingEffSummary.png")
    th2_summary_eff.Write()

    th2_eff.Draw("text colz")
    th2_eff.GetZaxis().SetRangeUser(0,1)
    cv.SaveAs(outputDir+"/effSummary.png")
    th2_eff.Write()

    # Add a legend
    legend = ROOT.TLegend(0.2, 0.2, 0.4, 0.3)
    ROOT.gStyle.SetOptTitle(0)

    th1_eff_seed.SetLineColor(ROOT.kRed)
    th1_eff_ckf.SetLineColor(ROOT.kBlue)
    th1_eff_ambi.SetLineColor(ROOT.kGreen)

    legend.AddEntry(th1_eff_seed, "Seeding", "lep")
    legend.AddEntry(th1_eff_ckf, "CKF", "lep")
    legend.AddEntry(th1_eff_ambi, "Ambiguity resolution", "lep")
    th1_eff_seed.GetYaxis().SetRangeUser(0,1)

    th1_eff_seed.Draw()
    th1_eff_ckf.Draw("same")
    th1_eff_ambi.Draw("same")
    legend.Draw()
    cv.SaveAs(outputDir+"/efficiency.png")

    th1_fake_seed.SetLineColor(ROOT.kRed)
    th1_fake_ckf.SetLineColor(ROOT.kBlue)
    th1_fake_ambi.SetLineColor(ROOT.kGreen)

    th1_fake_seed.GetYaxis().SetRangeUser(0,1)
    th1_fake_seed.Draw()
    th1_fake_ckf.Draw("same")
    th1_fake_ambi.Draw("same")
    legend.Draw()
    cv.SaveAs(outputDir+"/fakerate.png")

    th1_dup_seed.SetLineColor(ROOT.kRed)
    th1_dup_ckf.SetLineColor(ROOT.kBlue)
    th1_dup_ambi.SetLineColor(ROOT.kGreen)
    th1_dup_seed.GetYaxis().SetRangeUser(0,1)
    th1_dup_seed.Draw()
    th1_dup_ckf.Draw("same")
    th1_dup_ambi.Draw("same")
    legend.Draw()
    cv.SaveAs(outputDir+"/duplicaterate.png")

    th1_eff_seed.Write()
    th1_eff_ckf.Write()
    th1_eff_ambi.Write()

    th1_fake_seed.Write()
    th1_fake_ckf.Write()
    th1_fake_ambi.Write()

    th1_dup_seed.Write()
    th1_dup_ckf.Write()
    th1_dup_ambi.Write()
    th1_seedPerPart.Write()
    th1_seedPerPart.Draw()
    th1_seedPerPart.GetYaxis().SetRangeUser(0,14)
    cv.SaveAs(outputDir+"/SeedPerParticle.png")
    th1goodHitRate.Write()
    th1goodHitRate.Draw("text")
    cv.SaveAs(outputDir+"/goodHitRate.png")
    output.Close()

######
## example
######
files =[
            "testNOprim.out",
            "testNOsec.out"
        ]

dirlist =[
            "output_40GeV_newSeeding_standardSeeding_noPrimary_deadZones_maxSeedSpMPrim1_maxSeedSpMSec5_ImpMax1_dZMax50_branch1",
            "output_40GeV_newSeeding_standardSeeding_noSecondary_deadZones_maxSeedSpMPrim1_maxSeedSpMSec5_ImpMax1_dZMax50_branch1"

]


cutlist =[
            "secondary",
            "primary",
        ]

drawSummary(files, dirlist, cutlist, "New seeding")

