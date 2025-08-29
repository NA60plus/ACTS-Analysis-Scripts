import ROOT
import argparse
import os

# Set up argument parser
parser = argparse.ArgumentParser(description="Convert ROOT geometry files to GDML.")
parser.add_argument(
    "dirs", metavar="DIR", type=str, nargs="+",
    help="List of directories containing geometry.root files"
)
args = parser.parse_args()

# Loop through provided directories
for dir in args.dirs:
    root_file = os.path.join(dir, "geometry.root")
    gdml_file = os.path.join(dir, "geometry.gdml")

    print(f"Converting {root_file} to {gdml_file}...")

    ROOT.TGeoManager.Import(root_file)
    ROOT.gGeoManager.Export(gdml_file)
