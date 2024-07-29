import ROOT
import os
import csv
import re

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
directory = '/home/giacomo/FLUKA/150GeV/'
outdir = '/home/giacomo/FLUKA/csv-150GeV/'
suffix = '150GeV'
def create_csv_from_fluka(directory, outdir, suffix):
    # List all files and directories in the specified directory
    files_and_dirs = os.listdir(directory)
    # Filter out directories, keeping only files
    files = [f for f in files_and_dirs if os.path.isfile(os.path.join(directory, f))]
    histPDG = ROOT.TH1D("histPDG",";;Entries",19,0.5,19.5)
    for key,item in dict_part.items():
        histPDG.GetXaxis().SetBinLabel(item[0], inverted_dict[item[0]][0])


    histStation = ROOT.TH1D("histStation",";;Entries",6,-0.5,5.5)
    histSource = ROOT.TH1D("histSource",";;Entries",3,0.5,3.5)
    histSource.GetXaxis().SetBinLabel(1, "Target")
    histSource.GetXaxis().SetBinLabel(2, "Absorber")
    histSource.GetXaxis().SetBinLabel(3, "None")

    histxy = ROOT.TH2D("histxy",";x;y;Entries",200,-90,90,200,-90,90)
    histp = ROOT.TH1D("histp",";#it{p}_{T} ;Entries",100,0,10)
    histe = ROOT.TH1D("histe","; E (GeV/#it{c}^{2});Entries",100,0,10)
    histpt = ROOT.TH1D("histpt","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3)

    hist_dict = {
        "PROTON": [
                    ROOT.TH1D("histp_PROTON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_PROTON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_PROTON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_PROTON",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_PROTON",";;Entries",3,0.5,3.5)],
        "APROTON": [
                    ROOT.TH1D("histp_APROTON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_APROTON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_APROTON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_APROTON",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_APROTON",";;Entries",3,0.5,3.5)],
        "ELECTRON": [
                    ROOT.TH1D("histp_ELECTRON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_ELECTRON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_ELECTRON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_ELECTRON",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_ELECTRON",";;Entries",3,0.5,3.5)],
        "POSITRON": [
                    ROOT.TH1D("histp_POSITRON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_POSITRON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_POSITRON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_POSITRON",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_POSITRON",";;Entries",3,0.5,3.5)],
        "NEUTRIE": [
                    ROOT.TH1D("histp_NEUTRIE","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_NEUTRIE","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_NEUTRIE","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_NEUTRIE",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_NEUTRIE",";;Entries",3,0.5,3.5)],
        "ANEUTRIE": [
                    ROOT.TH1D("histp_ANEUTRIE","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_ANEUTRIE","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_ANEUTRIE","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_ANEUTRIE",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_ANEUTRIE",";;Entries",3,0.5,3.5)],
        "PHOTON": [
                    ROOT.TH1D("histp_PHOTON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_PHOTON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_PHOTON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_PHOTON",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_PHOTON",";;Entries",3,0.5,3.5)],
        "NEUTRON": [
                    ROOT.TH1D("histp_NEUTRON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_NEUTRON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_NEUTRON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_NEUTRON",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_NEUTRON",";;Entries",3,0.5,3.5)],
        "ANEUTRON": [
                    ROOT.TH1D("histp_ANEUTRON","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_ANEUTRON","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_ANEUTRON","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_ANEUTRON",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_ANEUTRON",";;Entries",3,0.5,3.5)],
        "MUON+": [
                    ROOT.TH1D("histp_MUON+","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_MUON+","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_MUON+","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_MUON+",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_MUON+",";;Entries",3,0.5,3.5)],
        "MUON-": [
                    ROOT.TH1D("histp_MUON-","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_MUON-","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_MUON-","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_MUON-",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_MUON-",";;Entries",3,0.5,3.5)],
        "KAONLONG": [
                    ROOT.TH1D("histp_KAONLONG","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_KAONLONG","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_KAONLONG","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_KAONLONG",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_KAONLONG",";;Entries",3,0.5,3.5)],
        "PION+": [
                    ROOT.TH1D("histp_PION+","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_PION+","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_PION+","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_PION+",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_PION+",";;Entries",3,0.5,3.5)],
        "PION-": [
                    ROOT.TH1D("histp_PION-","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_PION-","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_PION-","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_PION-",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_PION-",";;Entries",3,0.5,3.5)],
        "KAON+": [
                    ROOT.TH1D("histp_KAON+","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_KAON+","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_KAON+","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_KAON+",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_KAON+",";;Entries",3,0.5,3.5)],
        "KAON-": [
                    ROOT.TH1D("histp_KAON-","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_KAON-","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_KAON-","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_KAON-",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_KAON-",";;Entries",3,0.5,3.5)],
        "LAMBDA": [
                    ROOT.TH1D("histp_LAMBDA","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_LAMBDA","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_LAMBDA","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_LAMBDA",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_LAMBDA",";;Entries",3,0.5,3.5)],
        "ALAMBDA": [
                    ROOT.TH1D("histp_ALAMBDA","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_ALAMBDA","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_ALAMBDA","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_ALAMBDA",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSource_ALAMBDA",";;Entries",3,0.5,3.5)],
        "KAONSHRT": [
                    ROOT.TH1D("histp_KAONSHRT","; #it{p} (GeV/#it{c});Entries",100,0,10),
                    ROOT.TH1D("histe_KAONSHRT","; E (GeV/#it{c}^{2});Entries",100,0,10),
                    ROOT.TH1D("histpt_KAONSHRT","; #it{p}_{T} (GeV/#it{c});Entries",100,0,3),
                    ROOT.TH2D("histxy_KAONSHRT",";x;y;Entries",200,-90,90,200,-90,90),
                    ROOT.TH1D("histSourcey_KAONSHRT",";;Entries",3,0.5,3.5)]
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
                            if source == 2:
                                event += 1
                            data = [
                                ["particle_id","particle_type","process","vx","vy","vz","vt","px","py","pz","m","q"]
                            ]
                            break

                        
                    """
                    "elemento del set-up","particle id","luogo dove l'interaz. adronica ha avuto luogo","energia totale","x,y,z","coseni direttori"

                    PixStn2            3 AbsoPlu1   9.07811627E-04   14.2597504  -10.3486423       20.1179695       3.07518914E-02 6.16198704E-02             0.997625828
                    """

                    source = 3
                    for i in line:
                        if "Target" in i:
                            source = 2
                        if "Abso" in i:
                            source = 1
                    
                    if source != 2:
                        continue
                    #line = line.split()
                    if line[0] == "MS0":
                        
                        match = re.search(r'\d+', line[0])
                        if match:
                            number = match.group()
                            histStation.Fill(int(number))

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
                        hist_dict[inverted_dict[int(line[1])][0]][3].Fill(vx,vy)
                        
                        hist_dict[inverted_dict[int(line[1])][0]][4].Fill(source)
                        histSource.Fill(source)
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

    try:
        os.mkdir(suffix)
    except FileExistsError:
        print(f"Directory '{suffix}' already exists.")
    except PermissionError:
        print(f"Permission denied: Cannot create directory '{suffix}'.")

    print("nevent: ",event) 
    print("average molteplicity: ",histPDG.GetEntries()/event) 
    cv = ROOT.TCanvas("cv","cv",2000,2000)
    histPDG.Draw()
    cv.SaveAs(suffix+"/pdgcode"+suffix+".png")
    histp.Draw()
    cv.SaveAs(suffix+"/histp"+suffix+".png")
    histpt.Draw()
    cv.SaveAs(suffix+"/histt"+suffix+".png")
    histe.Draw()
    cv.SaveAs(suffix+"/histe"+suffix+".png")
    histxy.Draw()
    cv.SaveAs(suffix+"/histxy"+suffix+".png")
    histSource.Draw()
    cv.SaveAs(suffix+"/histSource"+suffix+".png")
    histStation.Draw()
    cv.SaveAs(suffix+"/histStation"+suffix+".png")

    output = ROOT.TFile(suffix+"/flukaout"+suffix+".root","recreate")
    histPDG.Write()
    histp.Write()
    histpt.Write()
    histe.Write()
    histxy.Write()
    histSource.Write()
    histStation.Write()
    for key, item in hist_dict.items():

        item[0].Write()
        item[1].Write()
        item[2].Write()
        item[3].Write()
        
        item[4].GetXaxis().SetBinLabel(1, "Target")
        item[4].GetXaxis().SetBinLabel(2, "Absorber")
        item[4].GetXaxis().SetBinLabel(3, "None")
        item[4].Write()
        if item[0].GetEntries() > 0:
            item[0].Draw()
            cv.SaveAs(suffix+"/"+item[0].GetName()+suffix+".png")
            item[1].Draw()
            cv.SaveAs(suffix+"/"+item[1].GetName()+suffix+".png")
            item[2].Draw()
            cv.SaveAs(suffix+"/"+item[2].GetName()+suffix+".png")
            item[3].Draw()
            cv.SaveAs(suffix+"/"+item[3].GetName()+suffix+".png")
        
    output.Close()
    

create_csv_from_fluka('/home/giacomo/FLUKA/40GeV/', '/home/giacomo/FLUKA/csv-40GeV/', '40GeV')
create_csv_from_fluka('/home/giacomo/FLUKA/80GeV/', '/home/giacomo/FLUKA/csv-80GeV/', '80GeV')
create_csv_from_fluka('/home/giacomo/FLUKA/120GeV/', '/home/giacomo/FLUKA/csv-120GeV/', '120GeV')
create_csv_from_fluka('/home/giacomo/FLUKA/150GeV/', '/home/giacomo/FLUKA/csv-150GeV/', '150GeV')


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
    merge_csv('/home/giacomo/FLUKA/csv-40GeV/event1-csv.csv', '/home/giacomo/FLUKA/csv-40GeV/event2-csv.csv', '/home/giacomo/FLUKA/event'+event+'-particles.csv.csv')