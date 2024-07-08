import uproot
import ROOT

def fill_dist(dists, hist, hist2):
    """
    Function to fill the histograms with the distance between the hits

    Args:
        dists (list): List of hits
        hist (TH1D): Histograms with the distance distribution
        hist2 (TH1D): Histograms with the distance distribution
    """
    for i in range(0, len(dists)):
        for j in range(i+1, len(dists)):
            hist.Fill(ROOT.TMath.Sqrt((dists[i][0]-dists[j][0])**2+(dists[i][1]-dists[j][1])**2))
            hist2.Fill(ROOT.TMath.Sqrt((dists[i][0]-dists[j][0])**2+(dists[i][1]-dists[j][1])**2))
            
            
            
def read_meas_tree(file_path, suffix, n_layer=5, station_side=5000):
    """
    Function to read the measurement tree to check the hit distributions in the stations
    Args:
        file_path (str): File path
        suffix (str): Ending string of the saved files
        n_layer (int, optional): Number of detectors. Defaults to 5.
        station_side (int, optional): Range of the histograms. Defaults to 5000.
    """
    tree_name = "vol1"
    # Open the ROOT file and get the TTree
    root_file = uproot.open(file_path)
    tree = root_file[tree_name]
    # Create a ROOT file to store the histograms
    output_file = ROOT.TFile("results/hits_distributions_"+suffix+".root", "RECREATE")

    # Get the branch
    x = tree["true_x"]
    y = tree["true_y"]
    layer = tree["layer_id"]
    event = tree["event_nr"]

    # Get the data from the branch as a numpy array
    x_data = x.array(library="np")
    y_data = y.array(library="np")
    layer_data = layer.array(library="np")
    event_data = event.array(library="np")

    ROOT.gStyle.SetOptStat(0)
    # Create a histogram for the branch
    
    nbins = 1000 #int(300/0.228)
    hist_xy = []
    for i in range(0, n_layer):
        hist_xy.append(ROOT.TH2F("hist_xy_"+str(i), "Plane 0;x [mm]; y [mm];", nbins, -int(station_side),int(station_side), nbins, -int(station_side),int(station_side)))
    
    minimum_dist = ROOT.TH1F("minimum_dist", " ;distance between hits [mm]; entries", 2000,0,200)
    minimum_dist_zoom = ROOT.TH1F("minimum_dist_zoom", " ;distance between hits [mm]; entries", int(1000/20),0,1)
    
    counter =0
    this_event = -1
    hits_list = []
    for x, y, layer, event in zip(x_data, y_data, layer_data, event_data):
        if this_event != event:
            if this_event != -1:
                fill_dist(hits_list, minimum_dist, minimum_dist_zoom)
            this_event = event
            print(event)
            hits_list = []
            
        #if this_event == -1 or this_event!=event_data:
        #    this_event = event
        if layer/2. < n_layer+1:
            counter += 1
            hist_xy[int(layer/2.)].Fill(x,y)
        if layer == 2:
            hits_list.append([x,y])
        
    # Write the histogram to the output file


    ROOT.gStyle.SetOptStat(1)
    minimum_dist.Write()
    minimum_dist_zoom.Write()
    print("mean distance: ",minimum_dist.GetMean())
    ROOT.gStyle.SetOptStat(0)

    cv  = ROOT.TCanvas("cv","cv",1600,1200)
    cv.SetRightMargin(0.2)

    for hist in hist_xy:
        hist.GetXaxis().SetLabelOffset(0.1)
        hist.Write()
        hist.Draw("colz")
        cv.SaveAs(hist.GetName()+"_"+suffix+".png")
        
    # Close the ROOT file and output file
    output_file.Close()
    root_file.close()


def main():
    read_meas_tree(file_path = "/home/giacomo/acts_for_NA60+/acts_na60plus_utils/muon_arm/measurements.root", suffix="muon_40GeV")
    
if __name__ == "__main__":
    main()
