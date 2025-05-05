mkdir Validation

root -l -b /home/giacomo/acts_for_NA60+/acts_old/Examples/Scripts/MaterialMapping/Mat_map.C'("geoRuben/propagation-material.root","geoRuben/material-map_tracks.root","Validation")'
.q

mkdir Surfaces
cd Surfaces
mkdir prop_plot
mkdir map_plot
mkdir ratio_plot
mkdir dist_plot
mkdir 1D_plot
cd ..

root -l -b /home/giacomo/acts_for_NA60+/acts_old/Examples/Scripts/MaterialMapping/Mat_map_surface_plot_ratio.C'("geoRuben/propagation-material.root","geoRuben/material-map_tracks.root",10000,"Surfaces/ratio_plot","Surfaces/prop_plot","Surfaces/map_plot")'
.q
root -l -b /home/giacomo/acts_for_NA60+/acts_old/Examples/Scripts/MaterialMapping/Mat_map_surface_plot_dist.C'("geoRuben/propagation-material.root",-1,"Surfaces/dist_plot")'
.q
root -l -b /home/giacomo/acts_for_NA60+/acts_old/Examples/Scripts/MaterialMapping/Mat_map_surface_plot_1D.C'("geoRuben/propagation-material.root",100000,"Surfaces/1D_plot")'
.q

root -l 
.x /home/giacomo/acts_for_NA60+/acts_old/Examples/Scripts/MaterialMapping/Mat_map.C("/home/arnaldi/cernbox/ACTS/ACTSsource_hits/Examples/Scripts/Python/propagation-material_Ruben_VTMS.root","/home/arnaldi/cernbox/ACTS/ACTSsource_hits/Examples/Scripts/Python/Mapping/material-map_tracks_Ruben_VTMS.root","Validation")
