#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <tuple>

// Compute rapidity given momentum magnitude (p), polar angle (theta), and mass
double compute_rapidity(double p, double theta, double mass = 0.1056583755) {
    double pz = p * std::cos(theta);
    double E = std::sqrt(p * p + mass * mass);
    
    if (E == pz) {
        return std::numeric_limits<double>::infinity();
    }
    return 0.5 * std::log((E + pz) / (E - pz));
}

def select_best_pairs(A, B, probabilities, max_to_min=True):
    """
    Seleziona le coppie con la probabilità più alta.
    :param A: Lista o array degli elementi di Detector A
    :param B: Lista o array degli elementi di Detector B
    :param probabilities: Matrice (n x n) con le probabilità di accoppiamento tra ogni elemento di A e ogni elemento di B
    :return: Lista di tuple (indice_A, indice_B, probabilità)
    """
    if max_to_min:
        // Convertiamo il problema di massimizzazione in un problema di minimizzazione
        cost_matrix = -probabilities
    else:
        cost_matrix = probabilities
    // Risolviamo il problema dell'assegnazione ottimale
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    // Creiamo la lista delle coppie selezionate
    selected_pairs = [
        (A[i], B[j], probabilities[i, j]) for i, j in zip(row_ind, col_ind)
    ]

    return selected_pairs


void run_application(
    const std::vector<double>& data_VT,
    const std::vector<double>& data_MS,
    const std::vector<double>& particles,
    const std::vector<std::string>& track_params,
    const std::string& input_dir = "",
    const std::string& output_dir = "",
    double chi2_selection = 3.0) {
    /*
    */
    std::string file_name = input_dir + std::string("/tracksummary_ambi.root");
    TFile *fileVT_step1 = TFile::Open(file_name.c_str());
    TTree *treeVT_step1 = (TTree *)fileVT->Get("tracksummary");

    file_name = input_dir + std::string("/tracksummary_ambill.root");
    TFile *fileVT_step2 = TFile::Open(file_name.c_str());
    TTree *treeVT_step2 = (TTree *)fileVT->Get("tracksummary");

    file_name = input_dir + std::string("/tracksummary_ambillll.root");
    TFile *fileVT_step3 = TFile::Open(file_name.c_str());
    TTree *treeVT_step3 = (TTree *)fileVT->Get("tracksummary");

    file_name = input_dir + std::string("/tracksummary_ambillllll.root");
    TFile *fileVT_step4 = TFile::Open(file_name.c_str());
    TTree *treeVT_step4 = (TTree *)fileVT->Get("tracksummary");

    file_name = input_dir + std::string("/tracksummary_ambims.root");
    TFile *fileMS = TFile::Open(file_name.c_str());
    TTree *treeMS = (TTree *)fileMS->Get("tracksummary");
    
    // Allocazione dei vettori per ciascuna variabile
    std::vector<float>* majorityParticleId = new std::vector<float>;
    std::vector<float>* nMeasurements     = new std::vector<float>;
    std::vector<float>* eLOC0_fit         = new std::vector<float>;
    std::vector<float>* eLOC1_fit         = new std::vector<float>;
    std::vector<float>* ePHI_fit          = new std::vector<float>;
    std::vector<float>* eTHETA_fit        = new std::vector<float>;
    std::vector<float>* eQOP_fit          = new std::vector<float>;

    std::vector<float>* cov_eLOC0_eLOC0   = new std::vector<float>;
    std::vector<float>* cov_eLOC0_eLOC1   = new std::vector<float>;
    std::vector<float>* cov_eLOC0_ePHI    = new std::vector<float>;
    std::vector<float>* cov_eLOC0_eTHETA  = new std::vector<float>;
    std::vector<float>* cov_eLOC0_eQOP    = new std::vector<float>;

    std::vector<float>* cov_eLOC1_eLOC0   = new std::vector<float>;
    std::vector<float>* cov_eLOC1_eLOC1   = new std::vector<float>;
    std::vector<float>* cov_eLOC1_ePHI    = new std::vector<float>;
    std::vector<float>* cov_eLOC1_eTHETA  = new std::vector<float>;
    std::vector<float>* cov_eLOC1_eQOP    = new std::vector<float>;

    std::vector<float>* cov_ePHI_eLOC0    = new std::vector<float>;
    std::vector<float>* cov_ePHI_eLOC1    = new std::vector<float>;
    std::vector<float>* cov_ePHI_ePHI     = new std::vector<float>;
    std::vector<float>* cov_ePHI_eTHETA   = new std::vector<float>;
    std::vector<float>* cov_ePHI_eQOP     = new std::vector<float>;

    std::vector<float>* cov_eTHETA_eLOC0  = new std::vector<float>;
    std::vector<float>* cov_eTHETA_eLOC1  = new std::vector<float>;
    std::vector<float>* cov_eTHETA_ePHI   = new std::vector<float>;
    std::vector<float>* cov_eTHETA_eTHETA = new std::vector<float>;
    std::vector<float>* cov_eTHETA_eQOP   = new std::vector<float>;

    std::vector<float>* cov_eQOP_eLOC0    = new std::vector<float>;
    std::vector<float>* cov_eQOP_eLOC1    = new std::vector<float>;
    std::vector<float>* cov_eQOP_ePHI     = new std::vector<float>;
    std::vector<float>* cov_eQOP_eTHETA   = new std::vector<float>;
    std::vector<float>* cov_eQOP_eQOP     = new std::vector<float>;

    std::vector<float>* t_p               = new std::vector<float>;
    std::vector<float>* t_theta           = new std::vector<float>;

    // Imposta gli indirizzi dei branch
    treeRef->SetBranchAddress("majorityParticleId", &majorityParticleId);
    treeRef->SetBranchAddress("nMeasurements",     &nMeasurements);
    treeRef->SetBranchAddress("eLOC0_fit",         &eLOC0_fit);
    treeRef->SetBranchAddress("eLOC1_fit",         &eLOC1_fit);
    treeRef->SetBranchAddress("ePHI_fit",          &ePHI_fit);
    treeRef->SetBranchAddress("eTHETA_fit",        &eTHETA_fit);
    treeRef->SetBranchAddress("eQOP_fit",          &eQOP_fit);

    treeRef->SetBranchAddress("cov_eLOC0_eLOC0",   &cov_eLOC0_eLOC0);
    treeRef->SetBranchAddress("cov_eLOC0_eLOC1",   &cov_eLOC0_eLOC1);
    treeRef->SetBranchAddress("cov_eLOC0_ePHI",    &cov_eLOC0_ePHI);
    treeRef->SetBranchAddress("cov_eLOC0_eTHETA",  &cov_eLOC0_eTHETA);
    treeRef->SetBranchAddress("cov_eLOC0_eQOP",    &cov_eLOC0_eQOP);

    treeRef->SetBranchAddress("cov_eLOC1_eLOC0",   &cov_eLOC1_eLOC0);
    treeRef->SetBranchAddress("cov_eLOC1_eLOC1",   &cov_eLOC1_eLOC1);
    treeRef->SetBranchAddress("cov_eLOC1_ePHI",    &cov_eLOC1_ePHI);
    treeRef->SetBranchAddress("cov_eLOC1_eTHETA",  &cov_eLOC1_eTHETA);
    treeRef->SetBranchAddress("cov_eLOC1_eQOP",    &cov_eLOC1_eQOP);

    treeRef->SetBranchAddress("cov_ePHI_eLOC0",    &cov_ePHI_eLOC0);
    treeRef->SetBranchAddress("cov_ePHI_eLOC1",    &cov_ePHI_eLOC1);
    treeRef->SetBranchAddress("cov_ePHI_ePHI",     &cov_ePHI_ePHI);
    treeRef->SetBranchAddress("cov_ePHI_eTHETA",   &cov_ePHI_eTHETA);
    treeRef->SetBranchAddress("cov_ePHI_eQOP",     &cov_ePHI_eQOP);

    treeRef->SetBranchAddress("cov_eTHETA_eLOC0",  &cov_eTHETA_eLOC0);
    treeRef->SetBranchAddress("cov_eTHETA_eLOC1",  &cov_eTHETA_eLOC1);
    treeRef->SetBranchAddress("cov_eTHETA_ePHI",   &cov_eTHETA_ePHI);
    treeRef->SetBranchAddress("cov_eTHETA_eTHETA", &cov_eTHETA_eTHETA);
    treeRef->SetBranchAddress("cov_eTHETA_eQOP",   &cov_eTHETA_eQOP);

    treeRef->SetBranchAddress("cov_eQOP_eLOC0",    &cov_eQOP_eLOC0);
    treeRef->SetBranchAddress("cov_eQOP_eLOC1",    &cov_eQOP_eLOC1);
    treeRef->SetBranchAddress("cov_eQOP_ePHI",     &cov_eQOP_ePHI);
    treeRef->SetBranchAddress("cov_eQOP_eTHETA",   &cov_eQOP_eTHETA);
    treeRef->SetBranchAddress("cov_eQOP_eQOP",     &cov_eQOP_eQOP);

    treeRef->SetBranchAddress("t_p",       &t_p);
    treeRef->SetBranchAddress("t_theta",   &t_theta);




    TH1D hPMatchedChi2("hPMatchedChi2", ";#it{p} (GeV/#it{c});Counts", 20, 0, 20);
    TH1D hPDistr("hPDistr", ";#it{p} (GeV/#it{c});Counts", 20, 0, 20);
    TH1D hYDistr("hYDistr", ";#it{y} (GeV/#it{c});Counts", 20, 0, 20);

    TH1D hElossMatchedChi2("hElossMatchedChi2", ";#Delta E (GeV/#it{c}^{2});Counts", 20, 0, 10);
    TH1D hElossDistr("hElossDistr", ";#Delta E (GeV/#it{c}^{2});Counts", 20, 0, 10);

    TH1D hElossRelMatchedChi2("hElossRelMatchedChi2", ";#Delta E / E ;Counts", 20, 0, 1);
    TH1D hElossRelDistr("hElossRelDistr", ";#Delta E / E ;Counts", 20, 0, 1);

    TH1D hYMatchedChi2("hYMatchedChi2", ";#it{y};Counts", 20, 0.5, 5.5);
    TH2D hYVsPMatchedChi2("hYVsPMatchedChi2", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20);
    
    TH2D hYVsPVT("hYVsPVT", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20);
    TH2D hYVsPAllChi2("hYVsPAllChi2", ";#it{y};#it{p} (GeV/#it{c});Counts", 20, 0.5, 5.5, 20, 0, 20);

    TH2D hChi2VsBDTY("hChi2VsBDTY", ";#it{y};#it{y};Counts", 20, 0.5, 5.5, 20, 0.5, 5.5);

    TH1D hYDistr("hYDistr", ";#it{y};Counts", 20, 0.5, 5.5);
    TH1D hYDistrChi2("hYDistrChi2", ";#it{y};Counts", 20, 0.5, 5.5);

    TH1D hChi2Matched("hChi2Matched", ";Matching #chi^{2};Counts (a.u.)", 200, 0, 40);
    TH1D hChi2Fake("hChi2Fake", ";Matching #chi^{2};Counts (a.u.)", 200, 0, 40);

    int fake_dileptons = 0;
    int true_dileptons = 0;
    int n_muons = 0;
    int vt_missing = 0;
    int n_muons_chi2 = 0;
    int n_muons_bdt = 0;
    int n_matches_chi2 = 0;
    int n_matches_bdt = 0;
    int event = 0;

    // loop over the events
    for event_VT, event_MS in zip(
        zip(*[data_VT[param] for param in track_params]),
        zip(*[data_MS[param] for param in track_params]),
    ):
        event += 1

        if (event % 100 == 0){
            std::cout<<"Processing event "<< event<<std::endl;
            std::cout<<"efficiency chi2: "<< n_matches_chi2 / n_muons_chi2<<std::endl;
            std::cout<<"vt_missing : "<< vt_missing / n_muons_chi2<<std::endl;
        }
        std::vector<int> y_list;
        std::vector<float> chi2_list;

        int matches = {0, 0};
        int fakes = {0, 0};
        n_muons += len(event_MS[0]);

        for i in range(len(event_MS[0])):
            // Build candidate pairs using all possible combinations
            pairsBDT = []
            labels_list = []
            // loop over the tracks
            for j in range(len(event_VT[0])):

                // Example residual vector
                r = np.array([event_VT[k][j] - event_MS[k][i] for k in range(2, 7)])
                // Example covariance matrix

                // Construct the covariance matrix
                C = np.zeros((5, 5))  // 5x5 matrix since k ranges from 2 to 6

                for k1 in range(5):
                    for k2 in range(5):
                        C[k1, k2] = event_VT[7+k1*5+k2][j] + event_MS[7+k1*5+k2][i]
                // Compute the inverse of the covariance matrix
                C_inv = np.linalg.inv(C)

                // Compute chi-square using matrix operations
                chi2 = r.T @ C_inv @ r
                chi2_list.append(chi2)

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

            y_list.append(np.array(labels_list))

        fake_dileptons += (
            fakes[0] * fakes[1] + fakes[0] * matches[1] + fakes[1] * matches[0]
        )
        true_dileptons += matches[0] * matches[1]

        ms = [i for i in range(len(event_MS[0]))]
        vt = [i for i in range(len(event_VT[0]))]

        chi2_matrix = np.array(chi2_list, dtype=float)

        if chi2_matrix.ndim == 1:
            chi2_matrix = chi2_matrix.reshape(len(ms), len(vt))

        chi2_best_pairs = select_best_pairs(ms, vt, chi2_matrix, False)


        for pairChi2 in zip(chi2_best_pairs):
            // Extract the particles array for the given event i
            index = list(particles[event-1]["particle_id"]).index(event_MS[0][pairChi2[0]])
            
            // Access the value of p for this particle
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

    
    TCanvas cv("cv", "cv", 1920, 1080);

        TEfficiency effVsPChi2(hPMatchedChi2, hPDistrChi2);
        TEfficiency effVsYChi2(hYMatchedChi2, hYDistrChi2);
        TEfficiency effVsYVsPChi2(hYVsPMatchedChi2, hYVsPAllChi2);
    
        TEfficiency effVsElossChi2(hElossMatchedChi2, hElossDistrChi2);
        TEfficiency effVsElossRelChi2(hElossRelMatchedChi2, hElossRelDistrChi2);
    
        effVsPChi2.SetName("effVsPChi2");
        effVsYChi2.SetName("effVsYChi2");
        effVsYVsPChi2.SetName("effVsYVsPChi2");
        effVsElossChi2.SetName("effVsElossChi2");
        effVsElossRelChi2.SetName("effVsElossRelChi2");
    
        effVsPChi2.SetTitle(";#it{p} (GeV/#it{c});Matching efficiency (%)");
        effVsYChi2.SetTitle(";#it{y};Matching efficiency (%)");
    
        TFile output_file((output_dir + "/matching.root").c_str(), "RECREATE");
    
        effVsPChi2.Write();
        effVsYChi2.Write();
        effVsYVsPChi2.Write();
        effVsElossChi2.Write();
        effVsElossRelChi2.Write();
    
        cv.Write();
        hPMatchedChi2.Write();
        hPDistr.Write();
        hYMatchedChi2.Write();
        hYDistr.Write();
        hBDTscoreMatched.Write();
        hChi2Matched.Write();
        hBDTscoreFake.Write();
        hChi2Fake.Write();
        hYVsPMatchedChi2.Write();
        hYVsPAllChi2.Write();
        hYVsPVT.Write();
    
        output_file.Close();
    std::cout<<"all mu bdt: "<< n_muons_bdt << std::endl;
    std::cout<<"all mu chi2: "<< n_muons_chi2 << std::endl;
    std::cout<<"matches mu bdt: "<< n_matches_bdt << std::endl;
    std::cout<<"matches mu chi2: "<< n_matches_chi2 << std::endl;
    std::cout<<"efficiency bdt: "<< n_matches_bdt / n_muons_bdt << std::endl;
    std::cout<<"efficiency chi2: "<< n_matches_chi2 / n_muons_chi2 << std::endl;
    std::cout<<"fakes dileptons = "<< fake_dileptons << std::endl;
    std::cout<<"true dileptons = "<< true_dileptons << std::endl;


    // Pulizia (delete) dei vettori allocati, se non gestiti automaticamente
    delete majorityParticleId;
    delete nMeasurements;
    delete eLOC0_fit;
    delete eLOC1_fit;
    delete ePHI_fit;
    delete eTHETA_fit;
    delete eQOP_fit;

    delete cov_eLOC0_eLOC0;
    delete cov_eLOC0_eLOC1;
    delete cov_eLOC0_ePHI;
    delete cov_eLOC0_eTHETA;
    delete cov_eLOC0_eQOP;

    delete cov_eLOC1_eLOC0;
    delete cov_eLOC1_eLOC1;
    delete cov_eLOC1_ePHI;
    delete cov_eLOC1_eTHETA;
    delete cov_eLOC1_eQOP;

    delete cov_ePHI_eLOC0;
    delete cov_ePHI_eLOC1;
    delete cov_ePHI_ePHI;
    delete cov_ePHI_eTHETA;
    delete cov_ePHI_eQOP;

    delete cov_eTHETA_eLOC0;
    delete cov_eTHETA_eLOC1;
    delete cov_eTHETA_ePHI;
    delete cov_eTHETA_eTHETA;
    delete cov_eTHETA_eQOP;

    delete cov_eQOP_eLOC0;
    delete cov_eQOP_eLOC1;
    delete cov_eQOP_ePHI;
    delete cov_eQOP_eTHETA;
    delete cov_eQOP_eQOP;

    delete t_p;
    delete t_theta;

    }

int main() {
    std::vector<std::string> track_params = {
        "majorityParticleId", "nMeasurements", "eLOC0_fit", "eLOC1_fit", "ePHI_fit", "eTHETA_fit", "eQOP_fit",
        "cov_eLOC0_eLOC0", "cov_eLOC0_eLOC1", "cov_eLOC0_ePHI", "cov_eLOC0_eTHETA", "cov_eLOC0_eQOP",
        "cov_eLOC1_eLOC0", "cov_eLOC1_eLOC1", "cov_eLOC1_ePHI", "cov_eLOC1_eTHETA", "cov_eLOC1_eQOP",
        "cov_ePHI_eLOC0", "cov_ePHI_eLOC1", "cov_ePHI_ePHI", "cov_ePHI_eTHETA", "cov_ePHI_eQOP",
        "cov_eTHETA_eLOC0", "cov_eTHETA_eLOC1", "cov_eTHETA_ePHI", "cov_eTHETA_eTHETA", "cov_eTHETA_eQOP",
        "cov_eQOP_eLOC0", "cov_eQOP_eLOC1", "cov_eQOP_ePHI", "cov_eQOP_eTHETA", "cov_eQOP_eQOP",
        "t_p", "t_theta"
    };
    return 0;
}
