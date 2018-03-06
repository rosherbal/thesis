# swG2M

This directory contains all the necessary files to reproduce the analysis of the result's chapter "Twists of kinases and phosphatases for a robust control of the G2/M transition".

- swCC-main.R -> Main code to call and run the functions to produce the results presented in the result's chapter.
- swCC-general.R -> General parameters required to run the functions called by swCC-main.R.
- swCC-models.R -> Contains the initial conditions for state variables and parameters, as well as the ODEs.
- MIsystem.lbs -> Contains tle LBS code of the Mutual Inhibitions (MI) system presented in the network emulation figure.
- G2Msystem.lbs -> Contains tle LBS code of the Gap 2 phase to mitosis (G2M) system presented in the network emulation figure.
