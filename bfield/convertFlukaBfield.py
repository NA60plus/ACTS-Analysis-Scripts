import re
import ROOT

ROOT.gStyle.SetTitleSize(0.05, "XYZ") # Sets the size for X, Y, and Z axis titles
ROOT.gStyle.SetTitleFontSize(0.05) # Sets the size for the histogram title
ROOT.gStyle.SetLabelSize(0.04, "XYZ")
ROOT.gStyle.SetPadLeftMargin(0.15)   #Set left margin
ROOT.gStyle.SetPadRightMargin(0.05)  #Set right margin
ROOT.gStyle.SetPadTopMargin(0.05)    #Set top margin
ROOT.gStyle.SetPadBottomMargin(0.15) #Set bottom margin
ROOT.gStyle.SetHistLineWidth(2)

def get_b_values(path="data/MEP48_field_map.inp"):
    # Read the file
    with open(path, "r") as file:
        lines = file.readlines()

    # Collect numbers from lines starting with MGNDATA
    bins = []
    range_x = []
    range_y = []
    range_z = []
    counter = 0
    numbers = []
    for line in lines:
        if line.startswith("MGNDATA"):
            # Extract all numbers, including scientific notation
            extracted_numbers = re.findall(
                r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?", line
            )
            if len(extracted_numbers) == 7:
                _ = extracted_numbers.pop(6)

                for i in range(0, len(extracted_numbers)):
                    extracted_numbers[i] = float(extracted_numbers[i])
            numbers.extend(extracted_numbers)

        if line.startswith("MGNCREAT"):
            counter += 1
            if counter == 2:
                bins = re.findall(
                    r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?", line
                )
                for i in range(0, len(bins)):
                    bins[i] = int(bins[i])
            elif counter == 3:
                extracted_numbers = re.findall(
                    r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?", line
                )
                range_x = [float(extracted_numbers[0]), float(extracted_numbers[3])]
                range_y = [float(extracted_numbers[1]), float(extracted_numbers[4])]
                range_z = [float(extracted_numbers[2]), float(extracted_numbers[5])]

    # Convert the extracted numbers to floats
    numbers = [float(num) for num in numbers]

    return bins, range_x, range_y, range_z, numbers


def from_line_to_matrix(line, nbin):
    matrix = []
    counter_z = 0
    z_bin = []
    for element in line:
        if counter_z == nbin:
            matrix.append(z_bin)
            z_bin = []
            counter_z = 0
        z_bin.append(element)
        counter_z += 1

    matrix.append(z_bin)

    return matrix


def rebin_matrix(matrix, bins):
    new_matrix = []
    for m in range(0,bins[2]):
        if m%2==0:#keep z ending with 5
            new_matrix.append(matrix[m])

    newest_matrix = []
    for m in range(0,len(new_matrix)):
        new_line = []
        for x in range(0,bins[0]):
            for y in range(0,bins[1]):
                if y%2==0 and x%2==0:
                    new_line.append(new_matrix[m][x+y*bins[0]])

        newest_matrix.append(new_line)
    return newest_matrix


"""

3
2
1 
0 1 2 3  


 4
 3
 2
 1 
 0  
-1 0 1 2 3 4 




"""


def extend_formatted_field(
    field, xmin, xmax, ymin, ymax, zmin, zmax, range_x, range_y, range_z, step, bins, zero = False
):
    
    print("field len", len(field))
    print("field[0] len", len(field[0]))
    # expand in z
    if range_z[0] > zmin:
        nstep = ROOT.TMath.Abs(range_z[0] - zmin) / step
        bins[2] += nstep
        range_z[0] = zmin
        for _ in range(0, int(nstep)):
            if zero:
                field.insert(0, [[0,0,0]] * int(bins[0]*bins[1]))
            else:
                field.insert(0, [[0,0,0]] * int(bins[0]*bins[1]))
    if range_z[1] < zmax:
        nstep = ROOT.TMath.Abs(range_z[1] - zmax) / step
        bins[2] += nstep
        range_z[1] = zmax
        for _ in range(0, int(nstep)):
            if zero:
                field.append([[0,0,0]] * int(bins[0]*bins[1]))
            else:
                field.append([[0,0,0]] * int(bins[0]*bins[1]))
    # expand in y
    new_field = []
    if range_y[0] > ymin:
        nstep = ROOT.TMath.Abs(range_y[0] - ymin) / step
        bins[1] += nstep
        range_y[0] = ymin
        for line in field:
            for _ in range(0, int(nstep)):
                if zero:
                    for _ in range(0, bins[0]):
                            line.insert(0, [0,0,0])
            if not zero:
                first_line = line[0:bins[0]]
                for _ in range(0, int(nstep)):
                    line = first_line + line
                new_field.append(line)
    if not zero:
        field = new_field
    new_field = []
    if range_y[1] < ymax:
        nstep = ROOT.TMath.Abs(range_y[1] - ymax) / step
        bins[1] += nstep
        range_y[1] = ymax
        for line in field:
            for _ in range(0, int(nstep)):
                for _ in range(0, bins[0]):
                    if zero:
                        line.append([0,0,0])

            if not zero:
                last_line = line[len(line)-bins[0]-1:-1]
                for _ in range(0, int(nstep)):
                    line = line+last_line

                new_field.append(line)
    if not zero:
        field = new_field

    print("field len", len(field))
    print("field[0] len", len(field[0]))

    # expand in x
    nstepmin = 0
    if range_x[0] > xmin:
        nstepmin = ROOT.TMath.Abs(range_x[0] - xmin) / step
        range_x[0] = xmin

    nstepmax = 0
    if range_x[1] < xmax:
        nstepmax = ROOT.TMath.Abs(range_x[1] - xmax) / step
        range_x[1] = xmax

    """
    for line in field:
        index = 0
        index_break = 0
        while index_break < bins[1]:
            index_break += 1
            for _ in range(0, int(nstepmin)):
                if zero:
                    line.insert(int(index), [0,0,0])
                else:
                    line.insert(int(index), line[int(index)])
            index += nstepmin + bins[0]
            for _ in range(0, int(nstepmax)):
                if zero:
                    line.insert(int(index), [0,0,0])
                else:
                    line.insert(int(index), line[int(index-1)])
            index += nstepmax
    """

    for line in field:
        index = 0
        for i in range(0, int(bins[1])):
            for n in range(0, int(nstepmin)):
                if zero:
                    line.insert(int(index), [0,0,0])
                else:
                    line.insert(int(index), line[int(index)-1])
            index += nstepmin + bins[0]

        index_end = nstepmin + bins[0]
        for _ in range(0, int(bins[1])):
            for _ in range(0, int(nstepmax)):
                if zero:
                    line.insert(int(index_end), [0,0,0])
                else:
                    line.insert(int(index_end), line[int(index_end)-1])
            index_end += nstepmax + nstepmin + bins[0]

    bins[0] += nstepmax + nstepmin
    bins[0] = int(bins[0])
    bins[1] = int(bins[1])
    bins[2] = int(bins[2])
    return field, range_x, range_y, range_z, bins


def convert_map(magnet="MEP48", zshift=0):
    path = "data/" + magnet + "_field_map.inp"

    bins, range_x, range_y, range_z, field = get_b_values(path)
    binsize = int((range_y[1]-range_y[0])/bins[0]+1)
    matrix = from_line_to_matrix(field, 3)
    matrix = from_line_to_matrix(matrix, len(matrix) / bins[2])
    range_z[0] = range_z[0] + zshift
    range_z[1] = range_z[1] + zshift


    matrix, range_x, range_y, range_z, bins = extend_formatted_field(matrix, -300, 300, -300, 300, -155, 1005, range_x, range_y, range_z, binsize, bins)

        
    hFieldxMax = ROOT.TH1D("hFieldxMax", ";#it{z} (cm);|#it{B}_{x}|^{max} (T)", bins[2], range_z[0], range_z[1])
    hFieldyMax = ROOT.TH1D("hFieldyMax", ";#it{z} (cm);|#it{B}_{y}|^{max} (T)", bins[2], range_z[0], range_z[1])
    hFieldzMax = ROOT.TH1D("hFieldzMax", ";#it{z} (cm);|#it{B}_{z}|^{max} (T)", bins[2], range_z[0], range_z[1])
    if magnet == "MEP48":

            
        hFieldxMax = ROOT.TH1D("hFieldxMax", ";#it{z} (cm);|#it{B}_{x}|^{max} (T)", 10, 0, 50)
        hFieldyMax = ROOT.TH1D("hFieldyMax", ";#it{z} (cm);|#it{B}_{y}|^{max} (T)", 10, 0, 50)
        hFieldzMax = ROOT.TH1D("hFieldzMax", ";#it{z} (cm);|#it{B}_{z}|^{max} (T)", 10, 0, 50)

    stepx = (range_x[1] - range_x[0] + binsize) / bins[0]
    stepy = (range_y[1] - range_y[0] + binsize) / bins[1]
    stepz = (range_z[1] - range_z[0] + binsize) / bins[2]


    for iz in range(0, bins[2]):

        zval = range_z[0] + iz * stepz
        bxvalMax = -1
        byvalMax = -1
        bzvalMax = -1

        for ix in range(0, bins[0]):
            xval = range_x[0] + ix * stepx
            if (xval>15 or xval<-15) and magnet == "MEP48":
                continue
            elif xval>100 or xval<-100:
                continue

            for iy in range(0, bins[1]):
                yval = range_y[0] + iy * stepy
                if (yval>15 or yval<-15) and magnet == "MEP48":
                    continue
                elif yval>100 or yval<-100:
                    continue
                if ROOT.TMath.Abs(matrix[iz][ix+iy*bins[0]][0]) > bxvalMax:
                    bxvalMax = ROOT.TMath.Abs(matrix[iz][ix+iy*bins[0]][0])
                if ROOT.TMath.Abs(matrix[iz][ix+iy*bins[0]][1]) > byvalMax:
                    byvalMax = ROOT.TMath.Abs(matrix[iz][ix+iy*bins[0]][1])
                if ROOT.TMath.Abs(matrix[iz][ix+iy*bins[0]][2]) > bzvalMax:
                    bzvalMax = ROOT.TMath.Abs(matrix[iz][ix+iy*bins[0]][2])
        hFieldxMax.Fill(zval, bxvalMax)
        hFieldyMax.Fill(zval, byvalMax)
        hFieldzMax.Fill(zval, bzvalMax)

    maxRangeX = hFieldxMax.GetMaximum()
    maxRangeY = hFieldyMax.GetMaximum()
    maxRangeZ = hFieldzMax.GetMaximum()
    
    if magnet == "MEP48":

        hFieldxMax.GetYaxis().SetRangeUser(0, 0.15)
        hFieldyMax.GetYaxis().SetRangeUser(0, 1.55)
        hFieldzMax.GetYaxis().SetRangeUser(0, 0.45)

        matrix = rebin_matrix(matrix,bins)
        bins[0] = int((bins[0]-1)/2+1)
        bins[1] = int((bins[1]-1)/2+1)
        bins[2] = int((bins[2]-1)/2+1)
        binsize =  2*binsize

    hFieldx = ROOT.TH1D("hFieldx", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldy = ROOT.TH1D("hFieldy", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldz = ROOT.TH1D("hFieldz", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldxmin = ROOT.TH1D("hFieldxmin", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldymin = ROOT.TH1D("hFieldymin", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldzmin = ROOT.TH1D("hFieldzmin", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldxplus = ROOT.TH1D("hFieldxplus", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldyplus = ROOT.TH1D("hFieldyplus", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldzplus = ROOT.TH1D("hFieldzplus", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldxminplus = ROOT.TH1D("hFieldxminplus", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldyminplus = ROOT.TH1D("hFieldyminplus", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldzminplus = ROOT.TH1D("hFieldzminplus", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldxplusmin = ROOT.TH1D("hFieldxplusmin", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldyplusmin = ROOT.TH1D("hFieldyplusmin", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldzplusmin = ROOT.TH1D("hFieldzplusmin", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldx.SetLineColor(ROOT.kBlack)
    hFieldy.SetLineColor(ROOT.kBlack)
    hFieldz.SetLineColor(ROOT.kBlack)
    
    hFieldxmin.SetLineColor(ROOT.kRed)
    hFieldxplus.SetLineColor(ROOT.kBlue)
    hFieldxminplus.SetLineColor(ROOT.kGreen + 1)
    hFieldxplusmin.SetLineColor(ROOT.kOrange + 1)

    hFieldymin.SetLineColor(ROOT.kRed)
    hFieldyplus.SetLineColor(ROOT.kBlue)
    hFieldyminplus.SetLineColor(ROOT.kGreen + 1)
    hFieldyplusmin.SetLineColor(ROOT.kOrange + 1)

    hFieldzmin.SetLineColor(ROOT.kRed)
    hFieldzplus.SetLineColor(ROOT.kBlue)
    hFieldzminplus.SetLineColor(ROOT.kGreen + 1)
    hFieldzplusmin.SetLineColor(ROOT.kOrange + 1)

    stepx = (range_x[1] - range_x[0] + binsize) / bins[0]
    stepy = (range_y[1] - range_y[0] + binsize) / bins[1]
    stepz = (range_z[1] - range_z[0] + binsize) / bins[2]

    file = open(magnet + ".txt", "w")
    plotrange = 10 if magnet == "MEP48" else 70
        
    for ix in range(0, bins[0]):
        xval = range_x[0] + ix * stepx
        #if xval % 10 == 5:
            #continue
        for iy in range(0, bins[1]):
            yval = range_y[0] + iy * stepy
            #if yval % 10 == 5:
                #continue
            for iz in range(0, bins[2]):
                zval = range_z[0] + iz * stepz
                #if zval % 10 == 0:
                    #continue

                bxval = matrix[iz][ix + iy * bins[0]][0]
                byval = matrix[iz][ix + iy * bins[0]][1]
                bzval = matrix[iz][ix + iy * bins[0]][2]
                file.write(f"{xval*10} {yval*10} {zval*10} {bxval} {byval} {bzval}\n")

                if yval == 0 and xval == 0 and magnet == "MEP48":
                    hFieldx.Fill(zval, bxval)
                    hFieldy.Fill(zval, byval)
                    hFieldz.Fill(zval, bzval)
                if yval == 40 and xval == 10 and magnet == "MNP33":
                    hFieldx.Fill(zval, bxval)
                    hFieldy.Fill(zval, byval)
                    hFieldz.Fill(zval, bzval)
                if yval == -plotrange and xval == -plotrange:
                    hFieldxmin.Fill(zval, bxval)
                    hFieldymin.Fill(zval, byval)
                    hFieldzmin.Fill(zval, bzval)
                if yval == +plotrange and xval == +plotrange:
                    hFieldxplus.Fill(zval, bxval)
                    hFieldyplus.Fill(zval, byval)
                    hFieldzplus.Fill(zval, bzval)
                if yval == -plotrange and xval == +plotrange:
                    hFieldxminplus.Fill(zval, bxval)
                    hFieldyminplus.Fill(zval, byval)
                    hFieldzminplus.Fill(zval, bzval)
                if yval == +plotrange and xval == -plotrange:
                    hFieldxplusmin.Fill(zval, bxval)
                    hFieldyplusmin.Fill(zval, byval)
                    hFieldzplusmin.Fill(zval, bzval)

    fileRot = open(magnet + "Rotated.txt", "w")
    for ix in range(0, bins[0]):
        xval = range_x[0] + ix * stepx
        for iz in range(0, bins[2]):
            zval = range_z[0] + iz * stepz
            for iy in range(bins[1] - 1, -1, -1):
                yval = range_y[0] + iy * stepy
                bxval = matrix[iz][ix + iy * bins[0]][0]
                byval = matrix[iz][ix + iy * bins[0]][1]
                bzval = matrix[iz][ix + iy * bins[0]][2]
                fileRot.write(
                    f"{xval*10} {zval*10} {-yval*10} {bxval} {bzval} {-byval}\n"
                )

    file.close()
    fileRot.close()

    cv = ROOT.TCanvas("cv", "cv", 3000, 1100)
    ROOT.gStyle.SetOptStat(0)
    cv.Divide(3)
    for i in range(1, 4):
        pad = cv.cd(i)
        pad.SetMargin(
            0.2, 0.1, 0.1, 0.1
        )  # Set left, right, bottom, top margins (in fraction of pad)
    
#    if magnet == "MEP48":

    cv.cd(1)
    hFieldxMax.Draw("hist") 

    vt1_x = ROOT.TLine(7.1175, 0, 7.1175, maxRangeX)
    vt2_x = ROOT.TLine(15.1175, 0, 15.1175, maxRangeX)
    vt3_x = ROOT.TLine(20.1175, 0, 20.1175, maxRangeX)
    vt4_x = ROOT.TLine(25.1175, 0, 25.1175, maxRangeX)
    vt5_x = ROOT.TLine(38.1175, 0, 38.1175, maxRangeX)

    ms1_x = ROOT.TLine(299.9695, 0, 299.9695, maxRangeX)
    ms2_x = ROOT.TLine(359.9695, 0, 359.9695, maxRangeX)
    ms3_x = ROOT.TLine(529.9695, 0, 529.9695, maxRangeX)
    ms4_x = ROOT.TLine(589.9695, 0, 589.9695, maxRangeX)
    ms5_x = ROOT.TLine(809.9695, 0, 809.9695, maxRangeX)
    ms6_x = ROOT.TLine(849.9695, 0, 849.9695, maxRangeX)

    vt1_x.SetLineColor(ROOT.kRed)
    vt2_x.SetLineColor(ROOT.kRed)
    vt3_x.SetLineColor(ROOT.kRed)
    vt4_x.SetLineColor(ROOT.kRed)
    vt5_x.SetLineColor(ROOT.kRed)
    ms1_x.SetLineColor(ROOT.kRed)
    ms2_x.SetLineColor(ROOT.kRed)
    ms3_x.SetLineColor(ROOT.kRed)
    ms4_x.SetLineColor(ROOT.kRed)
    ms5_x.SetLineColor(ROOT.kRed)
    ms6_x.SetLineColor(ROOT.kRed)

    vt1_x.SetLineWidth(2)
    vt2_x.SetLineWidth(2)
    vt3_x.SetLineWidth(2)
    vt4_x.SetLineWidth(2)
    vt5_x.SetLineWidth(2)
    ms1_x.SetLineWidth(2)
    ms2_x.SetLineWidth(2)
    ms3_x.SetLineWidth(2)
    ms4_x.SetLineWidth(2)
    ms5_x.SetLineWidth(2)
    ms6_x.SetLineWidth(2)

    if magnet == "MEP48":
        vt1_x.Draw("same")
        vt2_x.Draw("same")
        vt3_x.Draw("same")
        vt4_x.Draw("same")
        vt5_x.Draw("same")
    else:
        ms1_x.Draw("same")
        ms2_x.Draw("same")
        ms3_x.Draw("same")
        ms4_x.Draw("same")
        ms5_x.Draw("same")
        ms6_x.Draw("same")

    cv.cd(2)
    hFieldyMax.Draw("hist") 
    vt1_y = ROOT.TLine(7.1175, 0, 7.1175, maxRangeY)
    vt2_y = ROOT.TLine(15.1175, 0, 15.1175, maxRangeY)
    vt3_y = ROOT.TLine(20.1175, 0, 20.1175, maxRangeY)
    vt4_y = ROOT.TLine(25.1175, 0, 25.1175, maxRangeY)
    vt5_y = ROOT.TLine(38.1175, 0, 38.1175, maxRangeY)

    ms1_y = ROOT.TLine(299.9695, 0, 299.9695, maxRangeY)
    ms2_y = ROOT.TLine(359.9695, 0, 359.9695, maxRangeY)
    ms3_y = ROOT.TLine(529.9695, 0, 529.9695, maxRangeY)
    ms4_y = ROOT.TLine(589.9695, 0, 589.9695, maxRangeY)
    ms5_y = ROOT.TLine(809.9695, 0, 809.9695, maxRangeY)
    ms6_y = ROOT.TLine(849.9695, 0, 849.9695, maxRangeY)

    vt1_y.SetLineColor(ROOT.kRed)
    vt2_y.SetLineColor(ROOT.kRed)
    vt3_y.SetLineColor(ROOT.kRed)
    vt4_y.SetLineColor(ROOT.kRed)
    vt5_y.SetLineColor(ROOT.kRed)
    ms1_y.SetLineColor(ROOT.kRed)
    ms2_y.SetLineColor(ROOT.kRed)
    ms3_y.SetLineColor(ROOT.kRed)
    ms4_y.SetLineColor(ROOT.kRed)
    ms5_y.SetLineColor(ROOT.kRed)
    ms6_y.SetLineColor(ROOT.kRed)

    vt1_y.SetLineWidth(2)
    vt2_y.SetLineWidth(2)
    vt3_y.SetLineWidth(2)
    vt4_y.SetLineWidth(2)
    vt5_y.SetLineWidth(2)
    ms1_y.SetLineWidth(2)
    ms2_y.SetLineWidth(2)
    ms3_y.SetLineWidth(2)
    ms4_y.SetLineWidth(2)
    ms5_y.SetLineWidth(2)
    ms6_y.SetLineWidth(2)
    if magnet == "MEP48":
        vt1_y.Draw("same")
        vt2_y.Draw("same")
        vt3_y.Draw("same")
        vt4_y.Draw("same")
        vt5_y.Draw("same")
    else:
        ms1_y.Draw("same")
        ms2_y.Draw("same")
        ms3_y.Draw("same")
        ms4_y.Draw("same")
        ms5_y.Draw("same")
        ms6_y.Draw("same")

    cv.cd(3)
    hFieldzMax.Draw("hist") 

    vt1 = ROOT.TLine(7.1175, 0, 7.1175, maxRangeZ)
    vt2 = ROOT.TLine(15.1175, 0, 15.1175, maxRangeZ)
    vt3 = ROOT.TLine(20.1175, 0, 20.1175, maxRangeZ)
    vt4 = ROOT.TLine(25.1175, 0, 25.1175, maxRangeZ)
    vt5 = ROOT.TLine(38.1175, 0, 38.1175, maxRangeZ)

    ms1 = ROOT.TLine(299.9695, 0, 299.9695, maxRangeZ)
    ms2 = ROOT.TLine(359.9695, 0, 359.9695, maxRangeZ)
    ms3 = ROOT.TLine(529.9695, 0, 529.9695, maxRangeZ)
    ms4 = ROOT.TLine(589.9695, 0, 589.9695, maxRangeZ)
    ms5 = ROOT.TLine(809.9695, 0, 809.9695, maxRangeZ)
    ms6 = ROOT.TLine(849.9695, 0, 849.9695, maxRangeZ)

    vt1.SetLineColor(ROOT.kRed)
    vt2.SetLineColor(ROOT.kRed)
    vt3.SetLineColor(ROOT.kRed)
    vt4.SetLineColor(ROOT.kRed)
    vt5.SetLineColor(ROOT.kRed)
    ms1.SetLineColor(ROOT.kRed)
    ms2.SetLineColor(ROOT.kRed)
    ms3.SetLineColor(ROOT.kRed)
    ms4.SetLineColor(ROOT.kRed)
    ms5.SetLineColor(ROOT.kRed)
    ms6.SetLineColor(ROOT.kRed)

    vt1.SetLineWidth(2)
    vt2.SetLineWidth(2)
    vt3.SetLineWidth(2)
    vt4.SetLineWidth(2)
    vt5.SetLineWidth(2)
    ms1.SetLineWidth(2)
    ms2.SetLineWidth(2)
    ms3.SetLineWidth(2)
    ms4.SetLineWidth(2)
    ms5.SetLineWidth(2)
    ms6.SetLineWidth(2)
    if magnet == "MEP48":
        vt1.Draw("same")
        vt2.Draw("same")
        vt3.Draw("same")
        vt4.Draw("same")
        vt5.Draw("same")
    else:
        ms1.Draw("same")
        ms2.Draw("same")
        ms3.Draw("same")
        ms4.Draw("same")
        ms5.Draw("same")
        ms6.Draw("same")
    cv.SaveAs("figures/" + magnet + "_max_fieldvsz.png")

    cv.cd(1)
    hFieldx.Draw("hist") 
    cv.cd(2)
    hFieldy.Draw("hist") 
    cv.cd(3)
    hFieldz.Draw("hist") 

    cv.SaveAs("figures/" + magnet + "_fieldvsz.png")

    legend = ROOT.TLegend(0.4, 0.7, 0.7, 0.9)
    if magnet == "MEP48":
        legend.AddEntry(hFieldx, "x = 0 cm, y = 0 cm", "f")
        legend.AddEntry(hFieldxmin, "x = -15 cm, y = -15 cm", "f")
        legend.AddEntry(hFieldxplus, "x = 15 cm, y = 15 cm", "f")
        legend.AddEntry(hFieldxminplus, "x = -15 cm, y = 15 cm", "f")
        legend.AddEntry(hFieldxplusmin, "x = 15 cm, y = -15 cm", "f")
        hFieldx.GetYaxis().SetRangeUser(-0.2, 0.2)
        hFieldz.GetYaxis().SetRangeUser(-0.7, 0.7)
    else:
        legend.AddEntry(hFieldx, "x = 10 cm, y = 40 cm", "f")
        legend.AddEntry(hFieldxmin, "x = -70 cm, y = -70 cm", "f")
        legend.AddEntry(hFieldxplus, "x = 70 cm, y = 70 cm", "f")
        legend.AddEntry(hFieldxminplus, "x = -70 cm, y = 70 cm", "f")
        legend.AddEntry(hFieldxplusmin, "x = 70 cm, y = -70 cm", "f")
        hFieldx.GetYaxis().SetRangeUser(-0.025, 0.045)
        hFieldy.GetYaxis().SetRangeUser(0, 0.6)
        hFieldz.GetYaxis().SetRangeUser(-0.25, 0.25)

    cv.cd(1)
    hFieldxmin.Draw("hist same")
    hFieldxplus.Draw("hist same")
    hFieldxminplus.Draw("hist same")
    hFieldxplusmin.Draw("hist same")

    if magnet == "MNP33":
        legend.Draw("hist") 
    cv.cd(2)
    hFieldymin.Draw("hist same")
    hFieldyplus.Draw("hist same")
    hFieldyminplus.Draw("hist same")
    hFieldyplusmin.Draw("hist same")

    if magnet == "MEP48":
        legend.Draw("hist") 

    cv.cd(3)
    hFieldzmin.Draw("hist same")
    hFieldzplus.Draw("hist same")
    hFieldzminplus.Draw("hist same")
    hFieldzplusmin.Draw("hist same")
    cv.Update()
    cv.SaveAs("figures/" + magnet + "_fieldvsz_all.png")
    return matrix, range_x, range_y, range_z, bins


def merge_map(maps,binsize,range_x,range_y,range_z,bins):
    mep48 = maps[0]
    mnp33 = maps[1]

    # Initialize the result matrix with zeros
    matrix = [[[0,0,0] for _ in range(len(mep48[0]))] for _ in range(len(mep48))]

    # Sum element by element
    for i in range(len(mep48)):
        for j in range(len(mep48[0])):
            matrix[i][j][0] = mep48[i][j][0] + mnp33[i][j][0]
            matrix[i][j][1] = mep48[i][j][1] + mnp33[i][j][1]
            matrix[i][j][2] = mep48[i][j][2] + mnp33[i][j][2]

    hFieldxMax = ROOT.TH1D("hFieldxMax", ";#it{z} (cm);|#it{B}_{x}|^{max} (T)", bins[2], range_z[0], range_z[1])
    hFieldyMax = ROOT.TH1D("hFieldyMax", ";#it{z} (cm);|#it{B}_{y}|^{max} (T)", bins[2], range_z[0], range_z[1])
    hFieldzMax = ROOT.TH1D("hFieldzMax", ";#it{z} (cm);|#it{B}_{z}|^{max} (T)", bins[2], range_z[0], range_z[1])

    hFieldx = ROOT.TH1D("hFieldx", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldy = ROOT.TH1D("hFieldy", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldz = ROOT.TH1D("hFieldz", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldxmin = ROOT.TH1D("hFieldxmin", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldymin = ROOT.TH1D("hFieldymin", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldzmin = ROOT.TH1D("hFieldzmin", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldxplus = ROOT.TH1D("hFieldxplus", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldyplus = ROOT.TH1D("hFieldyplus", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldzplus = ROOT.TH1D("hFieldzplus", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldxminplus = ROOT.TH1D("hFieldxminplus", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldyminplus = ROOT.TH1D("hFieldyminplus", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldzminplus = ROOT.TH1D("hFieldzminplus", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldxplusmin = ROOT.TH1D("hFieldxplusmin", ";#it{z} (cm);Bx (T)", bins[2], range_z[0], range_z[1])
    hFieldyplusmin = ROOT.TH1D("hFieldyplusmin", ";#it{z} (cm);By (T)", bins[2], range_z[0], range_z[1])
    hFieldzplusmin = ROOT.TH1D("hFieldzplusmin", ";#it{z} (cm);Bz (T)", bins[2], range_z[0], range_z[1])

    hFieldx.SetLineColor(ROOT.kBlack)
    hFieldy.SetLineColor(ROOT.kBlack)
    hFieldz.SetLineColor(ROOT.kBlack)
    
    hFieldxmin.SetLineColor(ROOT.kRed)
    hFieldxplus.SetLineColor(ROOT.kBlue)
    hFieldxminplus.SetLineColor(ROOT.kGreen + 1)
    hFieldxplusmin.SetLineColor(ROOT.kOrange + 1)

    hFieldymin.SetLineColor(ROOT.kRed)
    hFieldyplus.SetLineColor(ROOT.kBlue)
    hFieldyminplus.SetLineColor(ROOT.kGreen + 1)
    hFieldyplusmin.SetLineColor(ROOT.kOrange + 1)

    hFieldzmin.SetLineColor(ROOT.kRed)
    hFieldzplus.SetLineColor(ROOT.kBlue)
    hFieldzminplus.SetLineColor(ROOT.kGreen + 1)
    hFieldzplusmin.SetLineColor(ROOT.kOrange + 1)

    stepx = (range_x[1] - range_x[0] + binsize) / bins[0]
    stepy = (range_y[1] - range_y[0] + binsize) / bins[1]
    stepz = (range_z[1] - range_z[0] + binsize) / bins[2]

    file = open("NewBFieldNA60plus_longsetup.txt", "w")
    plotrange = 10


    for iz in range(0, bins[2]):
        zval = range_z[0] + iz * stepz
        bxvalMax = ROOT.TMath.Abs(matrix[iz][0][0])
        byvalMax = ROOT.TMath.Abs(matrix[iz][0][1])
        bzvalMax = ROOT.TMath.Abs(matrix[iz][0][2])
        for m in matrix[iz]:
            if ROOT.TMath.Abs(m[0]) > bxvalMax:
                bxvalMax = ROOT.TMath.Abs(m[0])
            if ROOT.TMath.Abs(m[1]) > byvalMax:
                byvalMax = ROOT.TMath.Abs(m[1])
            if ROOT.TMath.Abs(m[2]) > bzvalMax:
                bzvalMax = ROOT.TMath.Abs(m[2])
        hFieldxMax.Fill(zval, bxvalMax)
        hFieldyMax.Fill(zval, byvalMax)
        hFieldzMax.Fill(zval, bzvalMax)

    for ix in range(0, bins[0]):
        xval = range_x[0] + ix * stepx
        #if xval % 10 == 5:
            #continue
        for iy in range(0, bins[1]):
            yval = range_y[0] + iy * stepy
            #if yval % 10 == 5:
                #continue
            for iz in range(0, bins[2]):
                zval = range_z[0] + iz * stepz
                #if zval % 10 == 0:
                    #continue

                bxval = matrix[iz][ix + iy * bins[0]][0]
                byval = matrix[iz][ix + iy * bins[0]][1]
                bzval = matrix[iz][ix + iy * bins[0]][2]

                file.write(f"{xval*10} {yval*10} {zval*10} {bxval} {byval} {bzval}\n")


                if yval == 0 and xval == 0:
                    hFieldx.Fill(zval, bxval)
                    hFieldy.Fill(zval, byval)
                    hFieldz.Fill(zval, bzval)
                if yval == -plotrange and xval == -plotrange:
                    hFieldxmin.Fill(zval, bxval)
                    hFieldymin.Fill(zval, byval)
                    hFieldzmin.Fill(zval, bzval)
                if yval == +plotrange and xval == +plotrange:
                    hFieldxplus.Fill(zval, bxval)
                    hFieldyplus.Fill(zval, byval)
                    hFieldzplus.Fill(zval, bzval)
                if yval == -plotrange and xval == +plotrange:
                    hFieldxminplus.Fill(zval, bxval)
                    hFieldyminplus.Fill(zval, byval)
                    hFieldzminplus.Fill(zval, bzval)
                if yval == +plotrange and xval == -plotrange:
                    hFieldxplusmin.Fill(zval, bxval)
                    hFieldyplusmin.Fill(zval, byval)
                    hFieldzplusmin.Fill(zval, bzval)

    fileRot = open("NewBFieldNA60plus_longsetupRotated.txt", "w")
    for ix in range(0, bins[0]):
        xval = range_x[0] + ix * stepx
        for iz in range(0, bins[2]):
            zval = range_z[0] + iz * stepz
            for iy in range(bins[1] - 1, -1, -1):
                yval = range_y[0] + iy * stepy
                bxval = matrix[iz][ix + iy * bins[0]][0]
                byval = matrix[iz][ix + iy * bins[0]][1]
                bzval = matrix[iz][ix + iy * bins[0]][2]
                fileRot.write(
                    f"{xval*10} {zval*10} {-yval*10} {bxval} {bzval} {-byval}\n"
                )

    file.close()
    fileRot.close()

    cv = ROOT.TCanvas("cv", "cv", 3000, 1100)
    ROOT.gStyle.SetOptStat(0)
    cv.Divide(3)
    for i in range(1, 4):
        pad = cv.cd(i)
        pad.SetMargin(
            0.2, 0.1, 0.1, 0.1
        )  # Set left, right, bottom, top margins (in fraction of pad)

    cv.cd(1)
    hFieldxMax.Draw("hist") 
    cv.cd(2)
    hFieldyMax.Draw("hist") 
    cv.cd(3)
    hFieldzMax.Draw("hist") 

    cv.SaveAs("figures/max_fieldvsz.png")

    cv.cd(1)
    hFieldx.Draw("hist") 
    cv.cd(2)
    hFieldy.Draw("hist") 
    cv.cd(3)
    hFieldz.Draw("hist") 

    cv.SaveAs("figures/all_fieldvsz.png")

    legend = ROOT.TLegend(0.4, 0.7, 0.7, 0.9)
    legend.AddEntry(hFieldx, "x = 0 cm, y = 0 cm", "f")
    legend.AddEntry(hFieldxmin, "x = -15 cm, y = -15 cm", "f")
    legend.AddEntry(hFieldxplus, "x = 15 cm, y = 15 cm", "f")
    legend.AddEntry(hFieldxminplus, "x = -15 cm, y = 15 cm", "f")
    legend.AddEntry(hFieldxplusmin, "x = 15 cm, y = -15 cm", "f")
    hFieldx.GetYaxis().SetRangeUser(-0.2, 0.2)
    hFieldz.GetYaxis().SetRangeUser(-0.7, 0.7)

    cv.cd(1)
    hFieldxmin.Draw("hist same")
    hFieldxplus.Draw("hist same")
    hFieldxminplus.Draw("hist same")
    hFieldxplusmin.Draw("hist same")

    legend.Draw("hist") 
    cv.cd(2)
    hFieldymin.Draw("hist same")
    hFieldyplus.Draw("hist same")
    hFieldyminplus.Draw("hist same")
    hFieldyplusmin.Draw("hist same")
    cv.cd(3)
    hFieldzmin.Draw("hist same")
    hFieldzplus.Draw("hist same")
    hFieldzminplus.Draw("hist same")
    hFieldzplusmin.Draw("hist same")
    cv.Update()
    cv.SaveAs("figures/all_fieldvsz_all.png")
    return matrix, range_x, range_y, range_z, bins

def main():
    MNP33_matrix, _, _, _, _ = convert_map("MNP33", 445)
    MEP48_matrix, MEP48_range_x, MEP48_range_y, MEP48_range_z, MEP48_bins = convert_map("MEP48", 0)
    merge_map([MEP48_matrix,MNP33_matrix],10,MEP48_range_x, MEP48_range_y, MEP48_range_z, MEP48_bins)


if __name__ == "__main__":
    main()