# Event generator
This directory contains the script to generate events in csv files that can be read by ACTS. The events contains prompt proton, pions, and kaons generated according to the NA49 measurements in Pb-Pb collisions at 40 AGeV and 158 AGeV. It can also generate K$^0_S\rightarrow \pi^+ + \pi^-$, $\Lambda^0 \rightarrow \text{K}^- + \pi^+$, and D$^0\rightarrow \pi^- + \text{p}$. The decays are simulated using AliDecayerEvtGen. 

To run the script:

```
.L event_generator.C+


event_generator(1000)
```
