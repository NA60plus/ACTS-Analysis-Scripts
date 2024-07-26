import ROOT
import os
import csv

dict_part = {
    "PROTON": [1, 2212],
    "APROTON": [2, -2212],
    "ELECTRON": [3, 11],
    "POSITRON": [4, -11],
    "NEUTRIE": [5, 12],
    "ANEUTRIE": [6, -12],
    "PHOTON": [7, 22],
    "NEUTRON": [8, 2112],
    "ANEUTRON": [9, -2112],
    "MUON+": [10, -13],
    "MUON-": [11, 13],
    "KAONLONG": [12, 130],
    "PION+": [13, 211],
    "PION-": [14, -211],
    "KAON+": [15, 321],
    "KAON-": [16, -321],
    "LAMBDA": [17, 3122],
    "ALAMBDA": [18, -3122],
    "KAONSHRT": [19, 310]
}

# Inverting the key and the first item
inverted_dict = {v[0]: [k,v[1]] for k, v in dict_part.items()}

# Specify the directory
directory = 'FLUKA/150GeV/'
outdir = 'FLUKA/csv-150GeV/'
suffix = '150GeV'
def create_csv_from_fluka(directory, outdir, suffix):
    # List all files and directories in the specified directory
    files_and_dirs = os.listdir(directory)
    # Filter out directories, keeping only files
    files = [f for f in files_and_dirs if os.path.isfile(os.path.join(directory, f))]
    histPDG = ROOT.TH1D("histPDG",";;Entries",19,0.5,19.5)
    for key,item in dict_part.items():
        histPDG.GetXaxis().SetBinLabel(item[0], inverted_dict[item[0]][0])

    histxy = ROOT.TH2D("histxy",";x;y;Entries",200,-400,400,200,-400,400)
    histp = ROOT.TH1D("histp",";#it{p}_{T} ;Entries",100,0,10)
    histe = ROOT.TH1D("histe","; E (GeV/#it{c}^{2});Entries",100,0,10)
    histpt = ROOT.TH1D("histpt","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)

    hist_dict = {
        "PROTON": [
                    ROOT.TH1D("histp_PROTON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_PROTON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_PROTON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)
                ],
        "APROTON": [
                    ROOT.TH1D("histp_APROTON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_APROTON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_APROTON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "ELECTRON": [
                    ROOT.TH1D("histp_ELECTRON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_ELECTRON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_ELECTRON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "POSITRON": [
                    ROOT.TH1D("histp_POSITRON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_POSITRON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_POSITRON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "NEUTRIE": [
                    ROOT.TH1D("histp_NEUTRIE","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_NEUTRIE","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_NEUTRIE","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "ANEUTRIE": [
                    ROOT.TH1D("histp_ANEUTRIE","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_ANEUTRIE","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_ANEUTRIE","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "PHOTON": [
                    ROOT.TH1D("histp_PHOTON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_PHOTON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_PHOTON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "NEUTRON": [
                    ROOT.TH1D("histp_NEUTRON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_NEUTRON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_NEUTRON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "ANEUTRON": [
                    ROOT.TH1D("histp_ANEUTRON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_ANEUTRON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_ANEUTRON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "MUON+": [
                    ROOT.TH1D("histp_MUON+","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_MUON+","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_MUON+","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "MUON-": [
                    ROOT.TH1D("histp_MUON-","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_MUON-","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_MUON-","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "KAONLONG": [
                    ROOT.TH1D("histp_KAONLONG","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_KAONLONG","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_KAONLONG","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "PION+": [
                    ROOT.TH1D("histp_PION+","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_PION+","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_PION+","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "PION-": [
                    ROOT.TH1D("histp_PION-","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_PION-","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_PION-","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "KAON+": [
                    ROOT.TH1D("histp_KAON+","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_KAON+","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_KAON+","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "KAON-": [
                    ROOT.TH1D("histp_KAON-","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_KAON-","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_KAON-","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "LAMBDA": [
                    ROOT.TH1D("histp_LAMBDA","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_LAMBDA","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_LAMBDA","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "ALAMBDA": [
                    ROOT.TH1D("histp_ALAMBDA","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_ALAMBDA","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_ALAMBDA","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)],
        "KAONSHRT": [
                    ROOT.TH1D("histp_KAONSHRT","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_KAONSHRT","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_KAONSHRT","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)]
    }

    event = 0
    for input_file in files:
        with open(directory+"/"+input_file, 'r') as infile:
                reader = csv.reader(infile, delimiter=' ', skipinitialspace=True)
                # Data to be written
                data = [
                    ["particle_id","particle_type","process","vx","vy","vz","vt","px","py","pz","m","q"]
                ]
                particleNumber = 0
                for line in reader:
                    for i in line:
                        if "event" in i:
                            
                            outfile = open(outdir+"event"+str(event)+"-csv.csv", 'w', newline='')
                            writer = csv.writer(outfile)
                            writer.writerows(data)
                            outfile.close()

                            particleNumber = 0
                            event += 1
                            data = [
                                ["particle_id","particle_type","process","vx","vy","vz","vt","px","py","pz","m","q"]
                            ]

                        
                    """
                    "elemento del set-up","particle id","luogo dove l'interaz. adronica ha avuto luogo","energia totale","x,y,z","coseni direttori"

                    PixStn2            3 AbsoPlu1   9.07811627E-04   14.2597504  -10.3486423       20.1179695       3.07518914E-02 6.16198704E-02             0.997625828
                    """

                    #line = line.split()
                    if line[0] == "MS0":
                        if line[-1] == "":  
                            popped = line.pop(-1)
                        pdg_code = inverted_dict[int(line[1])][1]
                        mass = ROOT.TDatabasePDG.Instance().GetParticle(pdg_code).Mass()
                        shift = 0
                        if len(line) == 9:
                            shift = 1
                        p = ROOT.TMath.Sqrt(float(line[3-shift])**2-mass**2)
                        histe.Fill(float(line[3-shift]))
                        vx = float(line[4-shift])
                        vy = float(line[5-shift])
                        vz = float(line[6-shift])
                        px = p*float(line[7-shift])
                        py = p*float(line[8-shift])
                        pz = p*float(line[9-shift])
                        histPDG.Fill(int(line[1]))
                        histpt.Fill(ROOT.TMath.Sqrt(px**2+py**2))
                        histp.Fill(ROOT.TMath.Sqrt(px**2+py**2+pz**2))
                        hist_dict[inverted_dict[int(line[1])][0]][0].Fill(ROOT.TMath.Sqrt(px**2+py**2))
                        hist_dict[inverted_dict[int(line[1])][0]][1].Fill(ROOT.TMath.Sqrt(px**2+py**2+pz**2))
                        hist_dict[inverted_dict[int(line[1])][0]][2].Fill(float(line[3-shift]))
                        histxy.Fill(vx,vy)
                        particleNumber += 1  # Exampl          e particleNumber value

                        # Convert event + 1 to a 12-bit binary string
                        event_binary = format(event + 1, '012b')

                        # Create a string of 12 zeros
                        zeros_12 = '0' * 12

                        # Convert particleNumber to a 16-bit binary string and increment particleNumber
                        particleNumber_binary = format(particleNumber, '016b')
                        particleNumber += 1

                        # Create a string of 24 zeros
                        zeros_24 = '0' * 24

                        # Concatenate all binary strings
                        binaryString = event_binary + zeros_12 + particleNumber_binary + zeros_24

                        # Convert the binary string to an unsigned long long integer
                        barcode = int(binaryString, 2)

                        # Write the data to the CSV file
                        data.append([barcode,pdg_code,0,vx,vy,vz,px,py,pz,mass,pdg_code/ROOT.TMath.Abs(pdg_code)])
            
    outfile.close()

    print("nevent: ",event) 
    print("average molteplicity: ",histPDG.GetEntries()/event) 
    cv = ROOT.TCanvas("cv","cv")
    histPDG.Draw()
    cv.SaveAs("pdgcode"+suffix+".png")

    output = ROOT.TFile("flukaout"+suffix+".root","recreate")
    histPDG.Write()
    histp.Write()
    histpt.Write()
    histe.Write()
    histxy.Write()
    for key,item in hist_dict.items():

        item[0].Write()
        item[1].Write()
        item[2].Write()
    output.Close()
    

#create_csv_from_fluka('FLUKA/40GeV/', 'FLUKA/csv-40GeV/', '40GeV')
#create_csv_from_fluka('FLUKA/80GeV/', 'FLUKA/csv-80GeV/', '80GeV')
#create_csv_from_fluka('FLUKA/120GeV/', 'FLUKA/csv-120GeV/', '120GeV')
#create_csv_from_fluka('FLUKA/150GeV/', 'FLUKA/csv-150GeV/', '150GeV')


def merge_csv(file1, file2, output_file):
    infile1 = open(file1, 'r')
    infile2 = open(file2, 'r')
    
    outfile = open(output_file, 'w', newline='')
    
    for line1 in infile1:
        outfile.write(line1)
    
    first = True
    for line2 in infile2:
        if first:
            first = False
            continue
        outfile.write(line2)
    
    outfile.close()
    infile1.close()
    infile2.close()
    
    print(f"Merged file saved as {output_file}")

# Example usage
directory_jpsi = ""
directory_bkg = ""

files_and_dirs_jpsi = os.listdir(directory_jpsi)
files_jpsi = [f for f in files_and_dirs_jpsi if os.path.isfile(os.path.join(directory_jpsi, f))]
files_and_dirs_bkg = os.listdir(directory_bkg)
files_bkg = [f for f in files_and_dirs_bkg if os.path.isfile(os.path.join(directory_bkg, f))]



for i, (jpi, bkg) in enumerate(zip(files_jpsi, files_bkg)):
    event = str(i)
    event = "0"*(9-len(event))+event
    merge_csv('FLUKA/csv-40GeV/event1-csv.csv', 'FLUKA/csv-40GeV/event2-csv.csv', 'FLUKA/event'+event+'-particles.csv.csv')