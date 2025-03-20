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
import math

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
import numpy as np

def compute_rho(data_VT, data_MS, k1, k2):
    """
    Compute the correlation coefficient rho between two measurement components k1 and k2.
    
    data_VT: List or array containing all VT events
    data_MS: List or array containing all MS events
    k1, k2: Indices of the measurement components (0-based, should map to k=2..6)
    
    Returns: Estimated rho (correlation coefficient)
    """
    
    # Extract measurements across multiple events
    X = np.array([event_VT[k1] - event_MS[k1] for event_VT, event_MS in zip(data_VT, data_MS)])
    Y = np.array([event_VT[k2] - event_MS[k2] for event_VT, event_MS in zip(data_VT, data_MS)])

    # Compute covariance and standard deviations
    covariance = np.cov(X, Y)[0, 1]  # Cov(X, Y)
    sigma_X = np.std(X, ddof=1)  # Standard deviation of X
    sigma_Y = np.std(Y, ddof=1)  # Standard deviation of Y

    # Compute rho, ensuring we don't divide by zero
    if sigma_X > 0 and sigma_Y > 0:
        rho = covariance / (sigma_X * sigma_Y)
    else:
        rho = 0  # No correlation if standard deviation is zero
    
    return rho


def compute_covariance(event_VT, event_MS, i, j, k1, k2, rho=0.1):
    """
    Compute off-diagonal covariance terms.
    
    This should be based on your specific knowledge of correlations between measurements.
    For example, assuming some correlation coefficient rho:
    """

    return rho * np.sqrt(
        (event_VT[k1 + 5][j]**2 + event_MS[k1 + 5][i]**2) *
        (event_VT[k2 + 5][j]**2 + event_MS[k2 + 5][i]**2)
    )

def load_data(track_params, path, split):
    # Load ROOT files and extract TTrees - use a list to simplify processing
    vt_files = [
        path + "/tracksummary_ambi.root",
        path + "/tracksummary_ambill.root",
        #path + "/tracksummary_ambillll.root",
        #path + "/tracksummary_ambillllll.root",
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

    particle_file = path + "/particles_simulation.root"
    tree_name = "particles"
    tree = uproot.open(particle_file)[tree_name]
    params = ["particle_id","e_loss","p"]
    particles = tree.arrays(params, library="ak")

    return data_VT_training, data_VT_application, data_MS_training, data_MS_application, particles[split:]


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
    particles,
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


    hElossMatchedBDT = ROOT.TH1D("hElossMatchedBDT", ";#Delta E (GeV/#it{c}^{2});Counts", 20, 0, 10)
    hElossMatchedChi2 = ROOT.TH1D("hElossMatchedChi2", ";#Delta E (GeV/#it{c}^{2});Counts", 20, 0, 10)
    hElossDistr = ROOT.TH1D("hElossDistr", ";#Delta E (GeV/#it{c}^{2});Counts", 20, 0, 10)
    hElossDistrChi2 = ROOT.TH1D("hElossDistrChi2", ";#Delta E (GeV/#it{c}^{2});Counts", 20, 0, 10)

    hElossRelMatchedBDT = ROOT.TH1D("hElossRelMatchedBDT", ";#Delta E / E ;Counts", 20, 0, 1)
    hElossRelMatchedChi2 = ROOT.TH1D("hElossRelMatchedChi2", ";#Delta E / E ;Counts", 20, 0, 1)
    hElossRelDistr = ROOT.TH1D("hElossRelDistr", ";#Delta E / E ;Counts", 20, 0, 1)
    hElossRelDistrChi2 = ROOT.TH1D("hElossRelDistrChi2", ";#Delta E / E ;Counts", 20, 0, 1)

    hYMatchedBDT = ROOT.TH1D("hYMatchedBDT", ";#it{y};Counts", 20, 0.5, 5.5)
    hYMatchedChi2 = ROOT.TH1D("hYMatchedChi2", ";#it{y};Counts", 20, 0.5, 5.5)
    hYVsPMatchedBDT = ROOT.TH2D("hYVsPMatchedBDT", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20)
    hYVsPMatchedChi2 = ROOT.TH2D("hYVsPMatchedChi2", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20)
    
    hYVsPVT = ROOT.TH2D("hYVsPVT", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20)

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
    vt_missing = 0
    n_muons_chi2 = 0
    n_muons_bdt = 0
    n_matches_chi2 = 0
    n_matches_bdt = 0
    event = 0

    # loop over the events
    for event_VT, event_MS in zip(
        zip(*[data_VT[param] for param in track_params]),
        zip(*[data_MS[param] for param in track_params]),
    ):
        event += 1

        if event % 100 == 0:
            print("Processing event ", event)
            print("efficiency bdt: ", n_matches_bdt / n_muons_bdt)
            print("efficiency chi2: ", n_matches_chi2 / n_muons_chi2)
            print("vt_missing : ", vt_missing / n_muons_chi2)

        y_list = []
        bdt_score_list = []
        chi2_list = []

        matches = [0, 0]
        fakes = [0, 0]
        n_muons += len(event_MS[0])

        for i in range(len(event_MS[0])):
            # Build candidate pairs using all possible combinations
            pairsBDT = []
            labels_list = []
            # loop over the tracks
            for j in range(len(event_VT[0])):

                # Example residual vector
                r = np.array([event_VT[k][j] - event_MS[k][i] for k in range(2, 7)])
                # Example covariance matrix

                # Construct the covariance matrix
                C = np.zeros((5, 5))  # 5x5 matrix since k ranges from 2 to 6

                for k1 in range(5):
                    for k2 in range(5):
                        C[k1, k2] = event_VT[7+k1*5+k2][j] + event_MS[7+k1*5+k2][i]
                # Compute the inverse of the covariance matrix
                C_inv = np.linalg.inv(C)

                # Compute chi-square using matrix operations
                chi2 = r.T @ C_inv @ r
                chi2_list.append(chi2)

                #####

                features = [
                    (event_VT[k+2][j] - event_MS[k+2][i]) ** 2
                    / (event_VT[7+k*5+k][j] ** 2 + event_MS[7+k*5+k][i] ** 2)
                    for k in range(5)
                ]
                features += [event_VT[k][j] for k in range(2, 7)]
                features += [event_MS[k][i] for k in range(2, 7)]
                features += [1 / event_VT[6][j] - 1 / event_MS[6][i]]
                features += [event_MS[6][i] * event_VT[6][j]]

                pairsBDT.append(features)
                if i == 0:
                    pVT= event_VT[32][j]
                    yVT = compute_rapidity(pVT,event_VT[13][j])
                    hYVsPVT.Fill(yVT,pVT)
                    
                if event_MS[0][i] == event_VT[0][j]:
                    if chi2 < chi2_selection:
                        matches[i] += 1
                    is_match = 1
                    hChi2Matched.Fill(chi2)
                else:
                    if chi2 < chi2_selection and event_MS[6][i] * event_VT[6][j] > 0:
                        fakes[i] += 1
                    is_match = 0
                    hChi2Fake.Fill(chi2)

                labels_list.append(is_match)
            if 1 not in labels_list:
                vt_missing += 1
            bdt_score = model.predict_proba(np.array(pairsBDT))[:, 1]

            for bdt, label in zip(bdt_score, labels_list):
                if label == 0:
                    hBDTscoreFake.Fill(bdt)
                else:
                    hBDTscoreMatched.Fill(bdt)

            y_list.append(np.array(labels_list))
            bdt_score_list.append(bdt_score)

        fake_dileptons += (
            fakes[0] * fakes[1] + fakes[0] * matches[1] + fakes[1] * matches[0]
        )
        true_dileptons += matches[0] * matches[1]

        ms = [i for i in range(len(event_MS[0]))]
        vt = [i for i in range(len(event_VT[0]))]

        # chi2_list = [(0 if i in check_vt_tracks) for i in range(len(event_MS[0]))]

        chi2_matrix = np.array(chi2_list, dtype=float)

        if chi2_matrix.ndim == 1:
            chi2_matrix = chi2_matrix.reshape(len(ms), len(vt))

        chi2_best_pairs = select_best_pairs(ms, vt, chi2_matrix, False)

        bdt_matrix = np.array(bdt_score_list, dtype=float)
        if bdt_matrix.ndim == 1:
            bdt_matrix = bdt_matrix.reshape(len(ms), len(vt))
        bdt_best_pairs = select_best_pairs(ms, vt, bdt_matrix)


        for pairChi2, pairBDT in zip(chi2_best_pairs,bdt_best_pairs):
            # Extract the particles array for the given event i
            index = list(particles[event-1]["particle_id"]).index(event_MS[0][pairBDT[0]])
            
            # Access the value of p for this particle
            eloss = particles[event-1]["e_loss"][index]
            pgen = particles[event-1]["p"][index]
            eloss_rel = eloss/math.sqrt(pgen**2+0.1056583755**2)

            hElossDistr.Fill(eloss)
            hElossDistrChi2.Fill(eloss)
            
            hElossRelDistr.Fill(eloss_rel)
            hElossRelDistrChi2.Fill(eloss_rel)

            pChi2 = event_MS[32][pairChi2[0]]
            yChi2 = compute_rapidity(pChi2,event_MS[33][pairChi2[0]])
            hPDistrChi2.Fill(pChi2)
            hYDistrChi2.Fill(yChi2)
            hYVsPAllChi2.Fill(yChi2,pChi2)
            n_muons_chi2 += 1
            if y_list[pairChi2[0]][pairChi2[1]] == 1:
                n_matches_chi2 += 1
                hPMatchedChi2.Fill(pChi2)
                hYMatchedChi2.Fill(yChi2)
                hYVsPMatchedChi2.Fill(yChi2,pChi2)
                hElossMatchedChi2.Fill(eloss)
                hElossRelMatchedChi2.Fill(eloss_rel)


            pBDT= event_MS[32][pairBDT[0]]
            yBDT = compute_rapidity(pBDT,event_MS[33][pairBDT[0]])
            hPDistr.Fill(pBDT)
            hYDistr.Fill(yBDT)
            hYVsPAllBDT.Fill(yBDT,pBDT)
            n_muons_bdt += 1
            if y_list[pairBDT[0]][pairBDT[1]] == 1:
                n_matches_bdt += 1
                hPMatchedBDT.Fill(pBDT)
                hYMatchedBDT.Fill(yBDT)
                hYVsPMatchedBDT.Fill(yBDT,pBDT)
                hElossMatchedBDT.Fill(eloss)
                hElossRelMatchedBDT.Fill(eloss_rel)

            hChi2VsBDTP.Fill(pChi2, pBDT)
            hChi2VsBDTY.Fill(yChi2, yBDT)
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

    effVsElossBDT = ROOT.TEfficiency(hElossMatchedBDT, hElossDistr)
    effVsElossChi2 = ROOT.TEfficiency(hElossMatchedChi2, hElossDistrChi2)
    effVsElossRelBDT = ROOT.TEfficiency(hElossRelMatchedBDT, hElossRelDistr)
    effVsElossRelChi2 = ROOT.TEfficiency(hElossRelMatchedChi2, hElossRelDistrChi2)

    effVsP.SetName("effVsP")
    effVsY.SetName("effVsY")
    effVsPBDT.SetName("effVsPBDT")
    effVsPChi2.SetName("effVsPChi2")
    effVsYBDT.SetName("effVsYBDT")
    effVsYChi2.SetName("effVsYChi2")
    effVsYVsPBDT.SetName("effVsYVsPBDT")
    effVsYVsPChi2.SetName("effVsYVsPChi2")

    effVsElossBDT.SetName("effVsElossBDT")
    effVsElossChi2.SetName("effVsElossChi2")
    effVsElossRelBDT.SetName("effVsElossRelBDT")
    effVsElossRelChi2.SetName("effVsElossRelChi2")

    effVsPBDT.SetTitle(";#it{p} (GeV/#it{c});Matching efficiency (%)")
    effVsPChi2.SetTitle(";#it{p} (GeV/#it{c});Matching efficiency (%)")
    effVsYBDT.SetTitle(";#it{y};Matching efficiency (%)")
    effVsYChi2.SetTitle(";#it{y};Matching efficiency (%)")

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

    effVsElossBDT.Write()
    effVsElossChi2.Write()
    effVsElossRelBDT.Write()
    effVsElossRelChi2.Write()

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
    hYVsPVT.Write()
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

        "cov_eLOC0_eLOC0",  # 7 1
        "cov_eLOC0_eLOC1",  # 8
        "cov_eLOC0_ePHI",  # 9
        "cov_eLOC0_eTHETA",  # 10
        "cov_eLOC0_eQOP",  # 11

        "cov_eLOC1_eLOC0",  # 12
        "cov_eLOC1_eLOC1",  # 13
        "cov_eLOC1_ePHI",  # 14
        "cov_eLOC1_eTHETA",  # 15
        "cov_eLOC1_eQOP",  # 16

        "cov_ePHI_eLOC0",  # 7
        "cov_ePHI_eLOC1",  # 8
        "cov_ePHI_ePHI",  # 9
        "cov_ePHI_eTHETA",  # 10
        "cov_ePHI_eQOP",  # 21

        "cov_eTHETA_eLOC0",  # 7
        "cov_eTHETA_eLOC1",  # 8
        "cov_eTHETA_ePHI",  # 9
        "cov_eTHETA_eTHETA",  # 10
        "cov_eTHETA_eQOP",  # 26

        "cov_eQOP_eLOC0",  # 7
        "cov_eQOP_eLOC1",  # 8
        "cov_eQOP_ePHI",  # 9
        "cov_eQOP_eTHETA",  # 10
        "cov_eQOP_eQOP",  #  31
        "t_p",  # 32
        "t_theta",  # 33
    ]

    data_VT_training, data_VT_application, data_MS_training, data_MS_application, particles = (
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
            data_VT_application, data_MS_application, particles, track_params, OUTPUT, CHI2
        )


if __name__ == "__main__":
    main()
