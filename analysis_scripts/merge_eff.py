import ROOT

def merge_efficiencies(path_list,suffix):
    """
    Function to merge the efficiencies for the different steps of the reconstruction

    Args:
        path_list (list): list of the root files with the efficiencies
        suffix (_type_): string to be attached to output file name
    """
    # Create a TCanvas and draw the TEfficiency histograms
    canvas = ROOT.TCanvas("canvas", "", 1600, 1200)
    # Open the ROOT files
    file_list = []
    for path in path_list:
        file_list.append(ROOT.TFile.Open(path))
        
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
        canvas.Print("eff_step"+str()+suffix+".png")
        if len(eff_list)==1:
            passed_sum = eff.GetPassedHistogram()
            total = eff.GetTotalHistogram()
        else:
            passed_sum.Add(eff)

    efficiency_final = ROOT.TEfficiency(passed_sum, total)
    efficiency_final.Draw()

    canvas.Draw()
    canvas.Print("results/eff_merged_"+suffix+".png")
    for file in file_list:
        file.Close()

var_list = ["pT","eta","phi"]
for var in var_list:
    merge_efficiencies(["output_40GeV_FullEv_doubleStep_standardSeeding_deadZones_maxSeedSpMPrim7_maxSeedSpMSec22_confirmation_primary_ImpMax0.57_dZMax7.37_branch1_eff1.0/performance_ambi.root",
                    "output_40GeV_FullEv_doubleStep_standardSeeding_deadZones_maxSeedSpMPrim7_maxSeedSpMSec22_confirmation_primary_ImpMax0.57_dZMax7.37_branch1_eff1.0/performance_ambiSecondary.root"],
                    var)