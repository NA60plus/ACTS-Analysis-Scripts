# ACTS-Analysis-Scripts
Repository for the script to run the NA60+ reconstruction chain and analyze the results
# Geometry
The vertex telescope geometry can be downloaded from https://cernbox.cern.ch/s/z2HWIp9xKg1Tmxg.
# Analysis of the data
The scripts the event generation are contained in event_generation.
# Run the reconstruction
The reconstruction of the generated events with the vertex telescope can be run with:

```
python3 full_chain_vt.py > printout.out
```
# Analysis of the data
The script for the analysis are contained in analysis_scripts.