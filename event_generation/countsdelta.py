def count_lines_in_file(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for line in file)

import os

def list_files_in_directory(directory_path):
    try:
        return [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]
    except FileNotFoundError:
        return f"Directory '{directory_path}' not found."

# Example usage
directory_path = '/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/event_generation/bkghits_40GeV_vt'
files_list = list_files_in_directory(directory_path)

import ROOT

th1 = ROOT.TH1D("th1",";N_{#delta};entries",100,0,2000)
# Example usage
for file_path in files_list:
    lines_count = count_lines_in_file(directory_path+"/"+file_path)-1
    th1.Fill(lines_count)

cv = ROOT.TCanvas("cv","cv",1600,1200)
th1.Draw()
cv.SaveAs("delta multiplicity.png")