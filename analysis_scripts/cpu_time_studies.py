
import os
import ROOT
import csv

def get_run_time(dir):
    """
    Function to get the running times
    Args:
        dir (str): directory with the timing information (timing.tsv)

    Returns:
        float, dictionary, list: total run time per event in s, dictionary with the run times per ev for each algorithm, list of the algorithm and cpu times
    """
    # List of algorithms we want to include in the running time calculation
    keeplist = [
            "Algorithm:TrackletVertexingAlgorithm",
            "Algorithm:SeedingAlgorithmNA60",
            "Algorithm:TrackParamsEstimationAlgorithm",
            "Algorithm:SeedsToPrototracks",
            "Algorithm:TrackFindingAlgorithm",
            "Algorithm:TracksToTrajectories",
            "Algorithm:GreedyAmbiguityResolutionAlgorithm"
            "Algorithm:FilterMeasurementsAlgorithm",
            "Algorithm:SpacePointMaker",
            "Algorithm:TruthSeedingAlgorithm",
            "Writer:RootTrackParameterWriter"
            "Algorithm:IterativeVertexFinder",
            "Writer:VertexPerformanceWriter"
        ]
    dict_time = {}
    list_time = []
    # Open the TSV file for reading
    with open(dir+'/timing.tsv', newline='') as tsvfile:
        # Create a CSV reader object with tab delimiter
        reader = csv.reader(tsvfile, delimiter='\t')
        tot_time = 0
        track_state_counter = 0
        track_summary_counter = 0
        # Iterate over each row in the file
        for row in reader:
            if row[0] in keeplist:
                # These writes are used twice: once during the CKF and once during the ambiguity resolutionA. We want only the ambiguity resolution.
                if "Writer:RootTrackStates" in row[0]:
                    track_state_counter += 1
                    if track_state_counter%2==0 :
                        continue
                if "Writer:RootTrackSummary" in row[0]:
                    track_summary_counter += 1
                    if track_summary_counter%2==0:
                        continue

                tot_time += float(row[2])
                if row[0] in dict_time:
                    dict_time[row[0]] += float(row[2])
                else:
                    dict_time[row[0]] = float(row[2])
                list_time.append(row)
        return tot_time, dict_time, list_time
    



def draw_runtime_pie(dict_time, title, directory = "../results/cpu_time_studies", suffix =""):
    """
    Function to draw a pie chart with the running time for each algorithm

    Args:
        dict_time (dict): dictionary with the running time for each algorithm
        title (str): pie chart title
        suffix (str, optional): ending string of the .png . Defaults to "".
    """
    dict_time = dict(sorted(dict_time.items(), key=lambda item: item[1], reverse=True))

    th1 = ROOT.TH1D("hist",title+";;",len(dict_time),0,len(dict_time))

    i = 1
    for val in dict_time.values():
        th1.SetBinContent(i, val)
        i+=1

    pie = ROOT.TPie(th1)
    pie.SetRadius(0.3)
    i = 0
    for key in dict_time.keys():
        pie.SetEntryFillColor(i, i+1)
        pie.SetEntryLabel(i, key.replace("Algorithm:","").replace("Writer:",""))
        i += 1

    cv = ROOT.TCanvas("cv","cv",1600,1600)
    #cv.SetRightMargin(0.)

    pie.Draw()

    pieleg = pie.MakeLegend()
    pieleg.SetY1(.6)
    pieleg.SetY2(1)
    pieleg.SetX1(.7)
    pieleg.SetX2(1)
 
    cv.SaveAs(directory+"/pie"+suffix+".png")

def draw_step_comparison_pie(list_time, algo = "Algorithm:TrackFindingAlgorithm", directory = "../results/cpu_time_studies", suffix =""):
    """
    Function to draw a pie chart with the running time of an algorithm for each step

    Args:
        list_time (_type_): list of the running times
        algo (str, optional): selected algorithm . Defaults to "Algorithm:TrackFindingAlgorithm".
        suffix (str, optional): ending string of the .png . Defaults to "".
    """
    steps = []
    tot_time = 0
    for al in list_time:
        if algo == al[0]:
            steps.append(al[2])
            tot_time+=float(al[2])
            
    th1 = ROOT.TH1D("hist","Time("+algo.replace("Algorithm:","")+")="+str(round(tot_time,3))+" s;;",len(steps),0,len(steps))

    i = 1
    for step in steps:
        th1.SetBinContent(i, float(step))
        i+=1

    pie = ROOT.TPie(th1)
    pie.SetRadius(0.3)
    i = 0
    for step in steps:
        pie.SetEntryFillColor(i, i+1)
        pie.SetEntryLabel(i, "Step "+str(i+1))
        i += 1

    cv = ROOT.TCanvas("cv","cv",1600,1600)

    pie.Draw()

    pieleg = pie.MakeLegend()
    pieleg.SetY1(.6)
    pieleg.SetY2(.9)
    pieleg.SetX1(.8)
    pieleg.SetX2(.95)
 
    cv.SaveAs(directory+"/pieStep"+suffix+".png")

def print_expected_cpu_time(fit_function, multiplicity_ratio, centrality_max = 50, multiplicity_scale_factor = 1, n_collisions = 10**11):
    """
    This function prints the expected CPU time to analyze one run.
    
    Args:
        fit_function (TF1): fit function of the CPU time vs normalized particle multiplicity
        multiplicity_ratio (list): list <npart>/<npart 0-5%> for each centrality class
        multiplicity_scale_factor (int, optional): _description_. Defaults to 1.
        n_collisions (_type_, optional): number of expected collision in one run. Defaults to 10**11.
    """
    mean_time = 0
    for rat in multiplicity_ratio:
        mean_time += fit_function.Eval(rat*multiplicity_scale_factor)

    mean_time/=len(multiplicity_ratio)
    second_to_year = 60.*60.*24.*365.
    print("Time per event:",mean_time, " s ")
    print("Total run time with "+str(n_collisions)+" MB collisions:",mean_time*n_collisions*centrality_max/100./second_to_year," y")




def main():

    #directories containing timing.tsv
    dir_list = [
                "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_realTarget_beam0.5_twosteps_rej0.114_perigeeZ400_suffix",
                "/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_40GeV_newSeeding_standardSeeding_deadZones_maxSeedSpMPrim1_maxSeedSpMSec20_ImpMax1_dZMax50_branch1_realTarget_beam0.5_fluka_twosteps_rej0.114_perigeeZ400_suffix"
        ]
    #list of names for the plots
    suffix_list = [
                "No fluka",
                "Fluka"
                ]
    #fraction of multiplicity in central collisions
    mult_list = [
                    0.01,
                    0.05
                ]


    # Directory name
    directory = "../results/cpu_time_studies"
    # Create the directory
    os.makedirs(directory, exist_ok=True)

    time_graph = ROOT.TGraphErrors(len(dir_list))
    time_graph.SetTitle("; <Multiplicity>/<Multiplicity 40 GeV, 0-5%>; CPU time / event (s)")

    for i, (dir, suffix, mult) in enumerate(zip(dir_list, suffix_list, mult_list)):
        tot_time, dict_time, list_time = get_run_time(dir=dir)
        print(dict_time)
        draw_runtime_pie(dict_time, "Time/event = " + str(round(tot_time, 3)) + " s", directory, suffix)
        draw_step_comparison_pie(list_time, directory=directory, suffix=suffix)
        time_graph.SetPoint(i, mult, tot_time)
        time_graph.SetPointError(i, 0, 0)  # Assuming no error, if there's an error value, replace 0 with the error

    fit_function = ROOT.TF1("fit_function", "[1]*TMath::Power(x,[0])")  # Linear fit
    time_graph.Fit(fit_function, "M+")  # "R" option for fitting in range
    
    cv = ROOT.TCanvas("cv", "cv", 1600, 1200)
    time_graph.SetMarkerStyle(21)
    time_graph.Draw("AP")
    t = ROOT.TText(.5,.5,"CPU time = C x N_{trk}^"+str(round(fit_function.GetParameter(0),2)))
    t.Draw()
    cv.SaveAs(directory+"/timevsMultiplicity.png")
    cv.SetLogy()
    cv.SaveAs(directory+"/timevsMultiplicityLog.png")

    # multiplicity_energy = pions + kaons + pions + decay dauthers from K0s->pipi and lambda->ppi
    multiplicity158 = 1258.000000 + 155.000000 +  292.000000 + 77.45*2*0.692 + 47.97*2*0.639
    multiplicity40 = 615.000000 + 78.000000 + 150.000000 + 39.15*2*0.692 + 43*2*0.639
    multiplicity20 = 410 + 51 + 148. + (39.15*2*0.692 + 43*2*0.639)/2.
    
    # list of <npart>/<npart 0-5%> for each 5% centrality class from 0 to 80%(50%) from a Glauber simulation with E = 158 AGeV
    ratio80 = [1.0, 0.7301204819277108, 0.6144578313253012, 0.5180722891566265, 0.43614457831325304, 0.363855421686747, 0.29759036144578316, 0.2421686746987952, 0.19518072289156627, 0.15421686746987953, 0.12048192771084337, 0.09156626506024096, 0.06746987951807229, 0.04819277108433735, 0.033734939759036145, 0.02289156626506024]
    ratio50 = [1.0, 0.7301204819277108, 0.6144578313253012, 0.5180722891566265, 0.43614457831325304, 0.363855421686747, 0.29759036144578316, 0.2421686746987952, 0.19518072289156627, 0.15421686746987953]

    print("==========================================================")
    print("E = 20 AGeV, centrality = 0-50%")
    print_expected_cpu_time(fit_function, ratio50, 50, multiplicity20/multiplicity40)
    print("E = 20 AGeV, centrality = 0-80%")
    print_expected_cpu_time(fit_function, ratio80, 80, multiplicity20/multiplicity40)
    print("E = 40 AGeV, centrality = 0-50%")
    print_expected_cpu_time(fit_function, ratio50, 50)
    print("E = 40 AGeV, centrality = 0-80%")
    print_expected_cpu_time(fit_function, ratio80, 80)
    print("E = 158 AGeV, centrality = 0-50%")
    print_expected_cpu_time(fit_function, ratio50, 50, multiplicity158/multiplicity40)
    print("E = 158 AGeV, centrality = 0-80%")
    print_expected_cpu_time(fit_function, ratio80, 80, multiplicity158/multiplicity40)


if __name__ == "__main__":
    main()