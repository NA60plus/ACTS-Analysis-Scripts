#!/bin/bash

# Define an array of directory names
directories=("geoRubenNewCarbon")
python_dir=/home/giacomo/acts_for_NA60+/acts_old/Examples/Scripts/Python
geo_dir=/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry
for dir in "${directories[@]}"; do
    echo "Processing $dir..."
    
    # Create the directory if it doesn't exist
    mkdir -p "$dir"
    
    # Convert response to JSON
    python3 /home/giacomo/acts_for_NA60+/acts_old/Examples/Scripts/Detectors/TGeoDetector/tgeo-response2json.py tgeo-json-config_muons_Ruben_VTMS.response > "$dir/tgeoRubenVol.json"

    # Clean up previous runs
    rm -rf obj/* csv/*

    # Run geometry conversion
    python3 $python_dir/geometry.py --jsonFile="$geo_dir/$dir/tgeoRubenVol.json" --tgeo_fileName="$geo_dir/$dir/geometry_Ruben.root"

    # Run material recording
    python3 $python_dir/material_recording.py --input="$geo_dir/$dir/geometry_Ruben.gdml" --tracks=100 -n 10000

    # Run material mapping
    python3 $python_dir/material_mapping.py --jsonFile="$geo_dir/$dir/tgeoRubenVol.json" --tgeo_fileName="$geo_dir/$dir/geometry_Ruben.root"

    # Run material validation
    python3 $python_dir/material_validation.py -o propagation-material -n 10000 --jsonFile="$geo_dir/$dir/tgeoRubenVol.json" --tgeo_fileName="$geo_dir/$dir/geometry_Ruben.root" > "$dir/test.out"

    # Move output files
    mv propagation-material.root "$dir/"
    mv geometry-map.json "$dir/"
    mv geant4_material_tracks.root "$dir/"
    mv material-map.json "$dir/"
    mv material-map_tracks.root "$dir/"

    echo "Finished processing $dir."
done
