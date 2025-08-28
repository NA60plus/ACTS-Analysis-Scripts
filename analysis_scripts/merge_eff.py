import ROOT
import copy


ROOT.gStyle.SetTitleSize(0.05, "XYZ") # Sets the size for X, Y, and Z axis titles
ROOT.gStyle.SetTitleFontSize(0.05) # Sets the size for the histogram title
ROOT.gStyle.SetLabelSize(0.04, "XYZ")
ROOT.gStyle.SetPadLeftMargin(0.15);   #Set left margin
ROOT.gStyle.SetPadRightMargin(0.05);  #Set right margin
ROOT.gStyle.SetPadTopMargin(0.05);    #Set top margin
ROOT.gStyle.SetPadBottomMargin(0.15); #Set bottom margin

def merge_efficiencies(path_list,suffix,name):
    """
    Function to merge the efficiencies for the different steps of the reconstruction

    Args:
        path_list (list): list of the root files with the efficiencies
        suffix (_type_): string to be attached to output file name
    """
    hs = ROOT.THStack("hs",";#it{p}_{T} (GeV/#it{c}); Tracking efficiency")
    # Create a TCanvas and draw the TEfficiency histograms
    canvas = ROOT.TCanvas("canvas", "", 1800, 1200)
    # Open the ROOT files
    file_list = []
    for path in path_list:
        file_list.append(ROOT.TFile.Open(path))
    
    color = [
        ROOT.kOrange + 2,
        ROOT.kRed -3 ,
        ROOT.kGreen -2,
        ROOT.kAzure - 2
    ]
    index = 0

    legend = ROOT.TLegend(0.45,0.2,0.6,0.4)
    # Get the trees from the files
    obj = "trackeff_vs_pT"
    if "eta" in suffix:
        obj = "trackeff_vs_eta"
    if "phi" in suffix:
        obj = "trackeff_vs_phi"
    eff_list = []
    for file in file_list:
        eff = file.Get(obj)
        eff_list.append(eff)
        eff.Draw()
        canvas.SaveAs(name+"eff_step"+str(len(eff_list))+suffix+".png")
        if len(eff_list)==1:
            passed_sum = eff.GetPassedHistogram()
            
            total = eff.GetTotalHistogram()

            #passed_sum.GetXaxis().SetRangeUser(0,2.5)
            #total.GetXaxis().SetRangeUser(0,2.5)

            passed_sum2 = copy.copy(eff.GetPassedHistogram())
            passed_sum2.Divide(total)

            passed_sum2.SetFillColor(color[index])
            passed_sum2.SetMarkerStyle(21)
            passed_sum2.SetMarkerColor(color[index])
            hs.Add(passed_sum2)
            print(passed_sum)
            print(total)
        else:
            passed_sum2 = copy.copy(eff.GetPassedHistogram())
            #passed_sum2.GetXaxis().SetRangeUser(0,2.5)

            passed_sum.Add(passed_sum2)
            
            #passed_sum2 = copy.copy(eff.GetPassedHistogram())


            passed_sum2.SetFillColor(color[index])
            passed_sum2.SetMarkerStyle(21)
            passed_sum2.SetMarkerColor(color[index])
            passed_sum2.Divide(total)

            hs.Add(passed_sum2)
        index+=1

        legend.AddEntry(passed_sum2, "Step "+str(index),"f")
    hs.Draw()
    legend.Draw()

    canvas.SaveAs(name+"stack.png")

    efficiency_final = ROOT.TEfficiency(passed_sum, total)
    efficiency_final.Draw()
    efficiency_final.SetTitle(";#it{p}_{T} (GeV/#it{c});Tracking efficiency")
    efficiency_final.SetMarkerStyle(20)
    canvas.Draw()
    canvas.SaveAs(name+"eff_merged_"+suffix+".png")
    for file in file_list:
        file.Close()


var_list = ["pT"]



files = [
        "output_40GeV_Sec_jpsi",
        "output_40GeV_Sec_jpsi_primaryOnly",
        "output_40GeV_Sec_jpsi_secondaryOnly"]
names = [
        "all",
        "primaryOnly",
        "secondaryOnly"]

for name,file in zip(names, files):
    for var in var_list:
        merge_efficiencies(["/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/"+file+"/performance_ambi.root",
                            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/"+file+"/performance_ambill.root",
                            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/"+file+"/performance_ambillll.root",
                            "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/"+file+"/performance_ambillllll.root"],
                            var,
                            name+"_ambi_")
    