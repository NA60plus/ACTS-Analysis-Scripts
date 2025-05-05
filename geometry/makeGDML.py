import ROOT

dirs = ["geoRubenNewCarbon"]
for dir in dirs:
    ROOT.TGeoManager.Import(dir+"/geometry_Ruben.root")
    ROOT.gGeoManager.Export(dir+"/geometry_Ruben.gdml")