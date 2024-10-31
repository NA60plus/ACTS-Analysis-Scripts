# Event generator
This directory contains the script to generate events in csv files that can be read by ACTS. The events contains prompt proton, pions, and kaons generated according to the NA49 measurements in Pb-Pb collisions at 40 AGeV and 158 AGeV[^1][^2]. It can also generate $\text{K}^0_S \rightarrow \pi^+ + \pi^-$, $\Lambda^0 \rightarrow \text{K}^- + \pi^+$ from NA49 measurements[^3], and $\text{D}^0 \rightarrow \pi^- + \text{p}$ from  POWHEG-BOX event generator[^4]. The decays are simulated using AliDecayerEvtGen. 

To run the script:

```
.L event_generator.C+


event_generator(1000)
```

You can download a set of events with prompt particles, $\text{K}^0_S \rightarrow \pi^+ + \pi^-$, and $\Lambda^0 \rightarrow \text{K}^- + \pi^+$ at 40 AGeV from https://cernbox.cern.ch/s/Ch80tIdItovM4R6.

# Background for the muon spectrometer
The background hits are taken from FLUKA simulations. The FLUKA simulations are MB events and contains both the particles created from the collision of a Pb ion with one of the targets, their decays daughters, the delta rays and the particles produced in the absorbers.
The hits are converted into csv files for ACTS. In the csv file for the VT only electrons and muons are saved.

The simulation and csv files can be found on cernbox at /eos/experiment/na60plus/MonteCarlo/Fluka .

To run the script:

```
.L bkg_generator.C+


bkg_generator()
```


[^1]: Afanasiev, Sergey V., et al. "Energy dependence of pion and kaon production in central Pb+ Pb collisions." Physical Review C 66.5 (2002): 054902.
[^2]: Antičić, Tome, et al. "Energy and centrality dependence of anti-p and p production and the anti-Lambda/anti-p ratio in Pb+ Pb collisions between 20A GeV and 158A GeV." (2006).
[^3]: Mischke, André, et al. "Lambda production in central Pb+ Pb collisions at CERN-SPS energies." Journal of Physics G: Nuclear and Particle Physics 28.7 (2002): 1761.
[^4]: Alioli, Simone, et al. "A general framework for implementing NLO calculations in shower Monte Carlo programs: the POWHEG BOX." Journal of High Energy Physics 2010.6 (2010): 1-58.