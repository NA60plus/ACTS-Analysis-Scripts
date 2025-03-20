import uproot
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score, confusion_matrix
import matplotlib.pyplot as plt
import awkward as ak
import ROOT
import argparse
import pickle
from scipy.optimize import linear_sum_assignment
import os
import time
import concurrent.futures


def parse_args():
    parser = argparse.ArgumentParser(description="Options")
    parser.add_argument(
        "-t", "--training", action="store_true", required=False, help="Run training"
    )
    parser.add_argument(
        "-a",
        "--application",
        action="store_true",
        required=False,
        help="Run application",
    )
    parser.add_argument(
        "-p",
        "--path",
        type=str,
        default="/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/output/output_omega_full_ruben",
        required=False,
        help="Path",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="omega_ruben",
        required=False,
        help="Name of the output",
    )
    parser.add_argument(
        "-s",
        "--split",
        type=int,
        default=1000,
        required=False,
        help="Events for the training",
    )
    parser.add_argument(
        "-c", "--chi2", type=float, default=6, required=False, help="Chi2 selection"
    )
    parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=4,
        required=False,
        help="Number of worker processes",
    )

    return parser.parse_args()

def compute_rapidity(p, theta, mass = 0.1056583755):
    """
    Compute the rapidity given momentum magnitude (p), polar angle (theta), and mass.
    
    Parameters:
    p (float): Momentum magnitude
    theta (float): Polar angle in radians
    mass (float): Mass of the particle
    
    Returns:
    float: Rapidity y
    """
    pz = p * np.cos(theta)  # Compute p_z
    E = np.sqrt(p**2 + mass**2)  # Compute energy
    
    if E == pz:
        return np.inf  # Avoid division by zero for extreme cases
    
    return 0.5 * np.log((E + pz) / (E - pz))

def load_data(track_params, path, split):
    # Load ROOT files and extract TTrees - use a list to simplify processing
    vt_files = [
        path + "/tracksummary_ambi.root",
        path + "/tracksummary_ambill.root",
        path + "/tracksummary_ambillll.root",
        path + "/tracksummary_ambillllll.root",
    ]
    detectorMS_file = path + "/tracksummary_ambims.root"
    tree_name = "tracksummary"

    # Use a more efficient approach to load and merge data
    data_VT_list = []
    for file_path in vt_files:
        tree = uproot.open(file_path)[tree_name]
        data_VT_list.append(tree.arrays(track_params, library="ak"))

    # Merge the VT data more efficiently
    merged_data = {}
    for key in track_params:
        arrays_to_concat = [data[key] for data in data_VT_list]
        merged_data[key] = ak.concatenate(arrays_to_concat, axis=1)

    data_VT = ak.zip(merged_data)

    # Load MS data
    tree_MS = uproot.open(detectorMS_file)[tree_name]
    data_MS = tree_MS.arrays(track_params, library="ak")

    # Split data for training and application
    data_VT_training = data_VT[:split]
    data_VT_application = data_VT[split:]
    data_MS_training = data_MS[:split]
    data_MS_application = data_MS[split:]

    return data_VT_training, data_VT_application, data_MS_training, data_MS_application


def select_best_pairs(A, B, probabilities, max_to_min=True):
    """
    Seleziona le coppie con la probabilità più alta.
    :param A: Lista o array degli elementi di Detector A
    :param B: Lista o array degli elementi di Detector B
    :param probabilities: Matrice (n x n) con le probabilità di accoppiamento tra ogni elemento di A e ogni elemento di B
    :return: Lista di tuple (indice_A, indice_B, probabilità)
    """
    if max_to_min:
        # Convertiamo il problema di massimizzazione in un problema di minimizzazione
        cost_matrix = -probabilities
    else:
        cost_matrix = probabilities
    # Risolviamo il problema dell'assegnazione ottimale
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Creiamo la lista delle coppie selezionate
    selected_pairs = [
        (A[i], B[j], probabilities[i, j]) for i, j in zip(row_ind, col_ind)
    ]

    return selected_pairs


def process_event_training(event_data):
    event_VT, event_MS = event_data
    pairs = []
    labels = []

    # Vectorize calculations where possible
    vt_params = np.array([event_VT[k] for k in range(2, 7)])
    ms_params = np.array([event_MS[k] for k in range(2, 7)])
    vt_errors = np.array([event_VT[k + 5] for k in range(2, 7)])
    ms_errors = np.array([event_MS[k + 5] for k in range(2, 7)])

    # Extract particle IDs for matching
    vt_ids = event_VT[0]
    ms_ids = event_MS[0]

    for i in range(len(ms_ids)):
        for j in range(len(vt_ids)):
            # Calculate features more efficiently
            diff_squared = np.square(vt_params[:, j] - ms_params[:, i])
            error_squared = np.square(vt_errors[:, j]) + np.square(ms_errors[:, i])
            chi2_components = diff_squared / error_squared

            features = list(chi2_components)
            features.extend(vt_params[:, j])
            features.extend(ms_params[:, i])
            features.append(1 / event_VT[6][j] - 1 / event_MS[6][i])
            features.append(event_MS[6][i] * event_VT[6][j])

            pairs.append(features)

            # Match determination
            if ms_ids[i].size > 0 and vt_ids[j].size > 0:
                is_match = 1 if ms_ids[i] == vt_ids[j] else 0
            else:
                is_match = 0
            labels.append(is_match)

    return pairs, labels


def run_training(
    data_VT, data_MS, track_params, output_dir="", model_name="bdt_model.pkl", workers=4
):
    start_time = time.time()

    # Prepare data for parallel processing
    event_data = list(
        zip(
            zip(*[data_VT[param] for param in track_params]),
            zip(*[data_MS[param] for param in track_params]),
        )
    )

    # Limit to first 1000 events for training
    event_data = event_data[:1000]

    # Process events in parallel
    all_pairs = []
    all_labels = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(process_event_training, event_data))

    # Combine results
    for pairs, labels in results:
        all_pairs.extend(pairs)
        all_labels.extend(labels)

    # Convert to NumPy arrays
    X = np.array(all_pairs)
    y = np.array(all_labels)

    # Split dataset into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )

    # Train XGBoost with optimized parameters
    bdt = xgb.XGBClassifier(
        n_estimators=100,
        max_depth=4,
        learning_rate=0.1,
        subsample=0.8,
        tree_method="hist",  # Faster algorithm
        n_jobs=workers,  # Use parallel processing
    )

    bdt.fit(X_train, y_train)

    # Evaluate the model
    y_pred = bdt.predict(X_test)
    bdt_score = bdt.predict_proba(X_test)[:, 1]

    accuracy = accuracy_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, bdt_score)
    conf_matrix = confusion_matrix(y_test, y_pred)

    print(f"Accuracy: {accuracy:.4f}")
    print(f"ROC AUC Score: {roc_auc:.4f}")
    print("Confusion Matrix:")
    print(conf_matrix)

    # Save model
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, model_name), "wb") as f:
        pickle.dump(bdt, f)

    end_time = time.time()
    print(f"Model training runtime: {end_time - start_time:.4f} seconds")

def run_application(
    data_VT,
    data_MS,
    track_params,
    output_dir="",
    chi2_selection=3,
    model_name="bdt_model.pkl",
):

    start_time = time.time()  # Start the timer
    with open(output_dir + "/" + model_name, "rb") as f:
        model = pickle.load(f)

    hPMatchedBDT = ROOT.TH1D("hPMatchedBDT", ";#it{p} (GeV/#it{c});Counts", 20, 0, 20)
    hPMatchedChi2 = ROOT.TH1D("hPMatchedChi2", ";#it{p} (GeV/#it{c});Counts", 20, 0, 20)
    hPDistr = ROOT.TH1D("hPDistr", ";#it{p} (GeV/#it{c});Counts", 20, 0, 20)
    hPDistrChi2 = ROOT.TH1D("hPDistrChi2", ";#it{p} (GeV/#it{c});Counts", 20, 0, 20)
    hYMatchedBDT = ROOT.TH1D("hYMatchedBDT", ";#it{y};Counts", 20, 0.5, 5.5)
    hYMatchedChi2 = ROOT.TH1D("hYMatchedChi2", ";#it{y};Counts", 20, 0.5, 5.5)
    hYVsPMatchedBDT = ROOT.TH2D("hYVsPMatchedBDT", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20)
    hYVsPMatchedChi2 = ROOT.TH2D("hYVsPMatchedChi2", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20)
    
    hYVsPAllBDT = ROOT.TH2D("hYVsPAllBDT", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20)
    hYVsPAllChi2 = ROOT.TH2D("hYVsPAllChi2", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20)

    hChi2VsBDTP = ROOT.TH2D("hChi2VsBDTP", ";#it{p} (GeV/#it{c});#it{p} (GeV/#it{c});Counts", 20, 0, 20, 20, 0, 20)
    hChi2VsBDTY = ROOT.TH2D("hChi2VsBDTY", ";#it{y};#it{y};Counts", 20, 0.5, 5.5, 20, 0.5, 5.5)
    

    hYDistr = ROOT.TH1D("hYDistr", ";#it{y};Counts", 20, 0.5, 5.5)
    hYDistrChi2 = ROOT.TH1D("hYDistrChi2", ";#it{y};Counts", 20, 0.5, 5.5)

    
    hBDTscoreMatched = ROOT.TH1D("hBDTscoreMatched", ";BDT probability;Counts (a.u.)", 200, 0, 1)
    hChi2Matched = ROOT.TH1D("hChi2Matched", ";Matching #chi^{2};Counts (a.u.)", 200, 0, 40)
    hBDTscoreFake = ROOT.TH1D("hBDTscoreFake", ";BDT probability;Counts (a.u.)", 200, 0, 1)
    hChi2Fake = ROOT.TH1D("hChi2Fake", ";Matching #chi^{2};Counts (a.u.)", 200, 0, 40)

    fake_dileptons = 0
    true_dileptons = 0

    n_muons = 0
    n_muons_chi2 = 0
    n_muons_bdt = 0
    n_matches_chi2 = 0
    n_matches_bdt = 0
    event = 0
    # loop over the events
    for event_num, (event_VT_data, event_MS_data) in enumerate(zip(
        zip(*[data_VT[param] for param in track_params]),
        zip(*[data_MS[param] for param in track_params]),
    ), 1):
        if event_num % 100 == 0:
            print("Processing event ", event_num)
        
        # Extract data once to avoid repeated indexing
        ms_count = len(event_MS_data[0])
        vt_count = len(event_VT_data[0])
        n_muons += ms_count
        
        # Pre-allocate arrays for better performance
        ms_indices = list(range(ms_count))
        vt_indices = list(range(vt_count))
        
        # Initialize tracking arrays
        matches = [0] * ms_count
        fakes = [0] * ms_count
        
        # Pre-compute matrices for better performance
        chi2_matrix = np.zeros((ms_count, vt_count), dtype=float)
        bdt_score_matrix = np.zeros((ms_count, vt_count), dtype=float)
        labels_matrix = np.zeros((ms_count, vt_count), dtype=int)
        
        # Precompute common values
        ms_values = {
            k: [event_MS_data[k][i] for i in range(ms_count)] for k in range(2, 7)
        }
        ms_errors = {
            k: [event_MS_data[k+5][i] for i in range(ms_count)] for k in range(2, 7)
        }
        vt_values = {
            k: [event_VT_data[k][j] for j in range(vt_count)] for k in range(2, 7)
        }
        vt_errors = {
            k: [event_VT_data[k+5][j] for j in range(vt_count)] for k in range(2, 7)
        }
        
        # Vectorized computation of features
        for i in range(ms_count):
            features_list = []
            
            for j in range(vt_count):
                # Compute chi2 value
                chi2 = sum(
                    (vt_values[k][j] - ms_values[k][i]) ** 2 / 
                    (vt_errors[k][j] ** 2 + ms_errors[k][i] ** 2)
                    for k in range(2, 7)
                )
                chi2_matrix[i, j] = chi2
                
                # Build feature vector
                features = [
                    (vt_values[k][j] - ms_values[k][i]) ** 2 / 
                    (vt_errors[k][j] ** 2 + ms_errors[k][i] ** 2)
                    for k in range(2, 7)
                ]
                features += [vt_values[k][j] for k in range(2, 7)]
                features += [ms_values[k][i] for k in range(2, 7)]
                features += [1 / event_VT_data[6][j] - 1 / event_MS_data[6][i]]
                features += [event_MS_data[6][i] * event_VT_data[6][j]]
                
                features_list.append(features)
                
                # Determine if this is a match
                is_match = 1 if event_MS_data[0][i] == event_VT_data[0][j] else 0
                labels_matrix[i, j] = is_match
                
                # Update histograms and counters
                if is_match:
                    if chi2 < chi2_selection:
                        matches[i] += 1
                    hChi2Matched.Fill(chi2)
                else:
                    if chi2 < chi2_selection and event_MS_data[6][i] * event_VT_data[6][j] > 0:
                        fakes[i] += 1
                    hChi2Fake.Fill(chi2)
            
            # Compute BDT scores in a batch
            bdt_scores = model.predict_proba(np.array(features_list))[:, 1]
            bdt_score_matrix[i] = bdt_scores
            
            # Update BDT histograms
            for j, (bdt, label) in enumerate(zip(bdt_scores, labels_matrix[i])):
                if label == 0:
                    hBDTscoreFake.Fill(bdt)
                else:
                    hBDTscoreMatched.Fill(bdt)
        
        # Compute dileptons
        fake_dileptons += sum(fakes[i] * (fakes[j] + matches[j]) for i in range(ms_count) for j in range(i+1, ms_count))
        true_dileptons += sum(matches[i] * matches[j] for i in range(ms_count) for j in range(i+1, ms_count))
        
        # Select best pairs
        chi2_best_pairs = select_best_pairs(ms_indices, vt_indices, chi2_matrix, False)
        bdt_best_pairs = select_best_pairs(ms_indices, vt_indices, bdt_score_matrix)
        
        # Process the selected pairs
        # Make sure both lists have the same length before zipping
        min_pairs = min(len(chi2_best_pairs), len(bdt_best_pairs))
        
        for i in range(min_pairs):
            # Handle each pair separately to avoid unpacking errors
            pair_chi2 = chi2_best_pairs[i]
            pair_bdt = bdt_best_pairs[i]
            
            # Check if the pairs are tuples with 2 elements as expected
            if isinstance(pair_chi2, tuple) and len(pair_chi2) == 2:
                ms_idx, vt_chi2_idx = pair_chi2
            else:
                # If not a proper tuple, skip this iteration
                continue
                
            if isinstance(pair_bdt, tuple) and len(pair_bdt) == 2:
                _, vt_bdt_idx = pair_bdt
            else:
                # If not a proper tuple, skip this iteration
                continue
            
            # Process Chi2 pair
            pChi2 = abs(1.0 / event_VT_data[6][vt_chi2_idx])
            yChi2 = compute_rapidity(pChi2, event_VT_data[5][vt_chi2_idx])
            
            hPDistrChi2.Fill(pChi2)
            hYDistrChi2.Fill(yChi2)
            hYVsPAllChi2.Fill(yChi2, pChi2)
            n_muons_chi2 += 1
            
            if labels_matrix[ms_idx, vt_chi2_idx] == 1:
                n_matches_chi2 += 1
                hPMatchedChi2.Fill(pChi2)
                hYMatchedChi2.Fill(yChi2)
                hYVsPMatchedChi2.Fill(yChi2, pChi2)
            
            # Process BDT pair
            pBDT = abs(1.0 / event_VT_data[6][vt_bdt_idx])
            yBDT = compute_rapidity(pBDT, event_VT_data[5][vt_bdt_idx])
            
            hPDistr.Fill(pBDT)
            hYDistr.Fill(yBDT)
            hYVsPAllBDT.Fill(yBDT, pBDT)
            n_muons_bdt += 1
            
            if labels_matrix[ms_idx, vt_bdt_idx] == 1:
                n_matches_bdt += 1
                hPMatchedBDT.Fill(pBDT)
                hYMatchedBDT.Fill(yBDT)
                hYVsPMatchedBDT.Fill(yBDT, pBDT)
            
            hChi2VsBDTP.Fill(pChi2, pBDT)
    ##########################################

    cv = ROOT.TCanvas("cv", "cv", 1920, 1080)

    effVsP = ROOT.TEfficiency(hPMatchedChi2, hPMatchedBDT)
    effVsY = ROOT.TEfficiency(hYMatchedChi2, hYMatchedBDT)
    effVsPBDT = ROOT.TEfficiency(hPMatchedBDT, hPDistr)
    effVsPChi2 = ROOT.TEfficiency(hPMatchedChi2, hPDistrChi2)
    effVsYBDT = ROOT.TEfficiency(hYMatchedBDT, hYDistr)
    effVsYChi2 = ROOT.TEfficiency(hYMatchedChi2, hYDistrChi2)
    effVsYVsPBDT = ROOT.TEfficiency(hYVsPMatchedBDT, hYVsPAllBDT)
    effVsYVsPChi2 = ROOT.TEfficiency(hYVsPMatchedChi2, hYVsPAllChi2)

    effVsP.SetName("effVsP")
    effVsY.SetName("effVsY")
    effVsPBDT.SetName("effVsPBDT")
    effVsPChi2.SetName("effVsPChi2")
    effVsYBDT.SetName("effVsYBDT")
    effVsYChi2.SetName("effVsYChi2")
    effVsYVsPBDT.SetName("effVsYVsPBDT")
    effVsYVsPChi2.SetName("effVsYVsPChi2")

    effVsPBDT.SetTitle(";#it{p} (GeV/#it{c});Matching efficiency (%)")
    effVsPChi2.SetTitle(";#it{p} (GeV/#it{c});Matching efficiency (%)")
    effVsYBDT.SetTitle(";#it{p} (GeV/#it{c});Matching efficiency (%)")
    effVsYChi2.SetTitle(";#it{p} (GeV/#it{c});Matching efficiency (%)")

    effVsPBDT.SetMarkerStyle(20)
    effVsPBDT.SetMarkerSize(1.2)
    effVsPBDT.SetLineWidth(2)

    effVsPChi2.SetMarkerStyle(20)
    effVsPChi2.SetMarkerSize(1.2)
    effVsPChi2.SetLineWidth(2)

    effVsPBDT.SetMarkerColor(ROOT.kBlue)
    effVsPBDT.SetLineColor(ROOT.kBlue)
    effVsPChi2.SetMarkerColor(ROOT.kRed)
    effVsPChi2.SetLineColor(ROOT.kRed)
    effVsPChi2.Draw()
    effVsPBDT.Draw("same")

    effVsYBDT.SetMarkerStyle(20)
    effVsYBDT.SetMarkerSize(1.2)
    effVsYBDT.SetLineWidth(2)

    effVsYChi2.SetMarkerStyle(20)
    effVsYChi2.SetMarkerSize(1.2)
    effVsYChi2.SetLineWidth(2)

    effVsYBDT.SetMarkerColor(ROOT.kBlue)
    effVsYBDT.SetLineColor(ROOT.kBlue)
    effVsYChi2.SetMarkerColor(ROOT.kRed)
    effVsYChi2.SetLineColor(ROOT.kRed)

    leg = ROOT.TLegend(0.25, 0.2, 0.4, 0.4)
    leg.AddEntry(effVsPChi2, "#chi^{2}", "lep")
    leg.AddEntry(effVsPBDT, "BDT", "lep")
    leg.Draw()
    cv.SaveAs(output_dir + "/matchingEff.png")

    output_file = ROOT.TFile(output_dir + "/matching.root", "recreate")
    effVsP.Write()
    effVsY.Write()
    effVsPBDT.Write()
    effVsPChi2.Write()
    effVsYBDT.Write()
    effVsYChi2.Write()
    effVsYVsPBDT.Write()
    effVsYVsPChi2.Write()
    cv.Write()
    hPMatchedBDT.Write()
    hPMatchedChi2.Write()
    hPDistr.Write()
    hPDistrChi2.Write()
    hYMatchedBDT.Write()
    hYMatchedChi2.Write()
    hYDistr.Write()
    hYDistrChi2.Write()

    hBDTscoreMatched.Write()
    hChi2Matched.Write()
    hBDTscoreFake.Write()
    hChi2Fake.Write()

    hYVsPMatchedBDT.Write()
    hYVsPMatchedChi2.Write()
    
    hYVsPAllBDT.Write()
    hYVsPAllChi2.Write()
    
    hChi2VsBDTP.Write()
    hChi2VsBDTY.Write()

    output_file.Close()
    
    print("all mu bdt: ", n_muons_bdt)
    print("all mu chi2: ", n_muons_chi2)
    print("matches mu bdt: ", n_matches_bdt)
    print("matches mu chi2: ", n_matches_chi2)
    print("efficiency bdt: ", n_matches_bdt / n_muons_bdt)
    print("efficiency chi2: ", n_matches_chi2 / n_muons_chi2)
    print("fakes dileptons = ", fake_dileptons)
    print("true dileptons = ", true_dileptons)
    end_time = time.time()  # End the timer

    print(f"Model application runtime time: {end_time - start_time:.4f} seconds")


def main():
    args = parse_args()
    TRAINING = args.training
    APPLICATION = args.application
    PATH = args.path
    OUTPUT = args.output
    SPLIT = args.split
    CHI2 = args.chi2

    track_params = [
        "majorityParticleId",  # 0
        "nMeasurements",  # 1
        "eLOC0_fit",  # 2
        "eLOC1_fit",  # 3
        "ePHI_fit",  # 4
        "eTHETA_fit",  # 5
        "eQOP_fit",  # 6
        "err_eLOC0_fit",  # 7
        "err_eLOC1_fit",  # 8
        "err_ePHI_fit",  # 9
        "err_eTHETA_fit",  # 10
        "err_eQOP_fit",  # 11
    ]

    data_VT_training, data_VT_application, data_MS_training, data_MS_application = (
        load_data(track_params, PATH, SPLIT)
    )

    # Create the directory
    os.makedirs(
        OUTPUT, exist_ok=True
    )  # 'exist_ok=True' prevents error if the directory already exists

    if TRAINING:
        run_training(data_VT_training, data_MS_training, track_params, OUTPUT)

    if APPLICATION:
        run_application(
            data_VT_application, data_MS_application, track_params, OUTPUT, CHI2
        )


if __name__ == "__main__":
    main()
