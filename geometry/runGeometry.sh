#!/bin/bash

# Define an array of directory names
acts_dir=/home/giacomo/acts_for_NA60+/acts/Examples/Scripts
run_mapping=true
plot_mapping=true
removeVS=false
removeMS=false
nevts=300
ntracks=100

opts=""
if $removeVS; then
  opts="$opts --remove-vs"
fi
if $removeMS; then
  opts="$opts --remove-ms"
fi

geo_dirs=("fullgeo")

mkdir -p obj
mkdir -p csv

for geo_dir in "${geo_dirs[@]}"; do
    echo "Processing $geo_dir..."
    
    # Set python_dir to the current directory being processed
    python_dir=$(pwd)/$geo_dir
    
    # Create the directory if it doesn't exist
    mkdir -p "$geo_dir"
    
    python3 makeGDML.py $geo_dir

    # Clean up previous runs
    rm -rf obj/* csv/*

    if $run_mapping; then
        # Run geometry conversion
        opts_tmp="${opts} --geometry-file $geo_dirs/geometry.root"
        python3 $acts_dir/Python/dice_geometry.py $opts_tmp
        
        python3 geomap_modifier.py -o "geometry-map.json" -n "geometry-map-new.json"

        # Run material recording
        python3 $acts_dir/Python/material_recording_dice.py --input="$geo_dir/geometry.gdml" --tracks=$ntracks -n $nevts

        # Run material mapping
        python3 $acts_dir/Python/material_mapping_dice.py -i "geometry-map-new.json" -o "material-map.json" $opts_tmp  -n $nevts

        # Run material validation
        python3 $acts_dir/Python/material_validation_dice.py -o propagation-material -n $nevts  -m "material-map.json" -t $ntracks  $opts_tmp
        # Move output files
        mv propagation-material.root "$geo_dir/"
        mv geometry-map.json "$geo_dir/"
        mv geometry-map-new.json "$geo_dir/"
        mv geant4_material_tracks.root "$geo_dir/"
        mv material-map.json "$geo_dir/"
        mv material-map_tracks.root "$geo_dir/"
    fi

    if $plot_mapping; then
        
        if $use_old_acts; then
            acts_dir=/home/giacomo/acts_for_NA60+/acts_old/Examples/Scripts
        fi
        mkdir -p "$python_dir/Validation"

        echo ${acts_dir}/MaterialMapping/Mat_map.C
        
         root -l -b $acts_dir/MaterialMapping/Mat_map.C'("'"$python_dir/propagation-material.root"'","'"$python_dir/material-map_tracks.root"'","'"$python_dir/Validation"'")' <<< '.q'
        
        mkdir -p "$python_dir/Surfaces"
        mkdir -p "$python_dir/Surfaces/prop_plot"
        mkdir -p "$python_dir/Surfaces/map_plot"
        mkdir -p "$python_dir/Surfaces/ratio_plot"
        mkdir -p "$python_dir/Surfaces/dist_plot"
        mkdir -p "$python_dir/Surfaces/1D_plot"
        
        root -l -b $acts_dir/MaterialMapping/Mat_map_surface_plot_ratio.C'("'"$python_dir/propagation-material.root"'","'"$python_dir/material-map_tracks.root"'",10000,"'"$python_dir/Surfaces/ratio_plot"'","'"$python_dir/Surfaces/prop_plot"'","'"$python_dir/Surfaces/map_plot"'")' <<< '.q'
        
        root -l -b $acts_dir/MaterialMapping/Mat_map_surface_plot_dist.C'("'"$python_dir/propagation-material.root"'",-1,"'"$python_dir/Surfaces/dist_plot"'")' <<< '.q'
        
        root -l -b $acts_dir/MaterialMapping/Mat_map_surface_plot_1D.C'("'"$python_dir/propagation-material.root"'",100000,"'"$python_dir/Surfaces/1D_plot"'")' <<< '.q'
        

    fi

    #echo "Finished processing $geo_dir."
done
