# README

## Excitable Tissue
Vertex model and effective model for an excitable tissue, as detailed in "Excitable dynamics driven by mechanical feedback in epithelial tissues".

## Authors / Contributors
- Fernanda Pérez-Verdugo (Carnegie Mellon University)
- Samuel Banks (Yale University)
- Shiladitya Banerjee (Carnegie Mellon University)

## System Requirements
- Simulation code is run in Python 3
- Numerical analisis is run in Mathematica 12

## QUICKSTART GUIDE
- /Simulations/ directory contains three directories with all the examples needed to re-do all the simulations ran in the study. Specifically, the examples includes the ones shown in Fig1D, Fig4A, Fig5I.
Details:

--- MainFile.py: Main file to run the simulation. Results are saved in the directories "length_tension" and "results".

--- snapshots.py: This script generate tissue snapshots using the data from "/length_tension/" and "/results/". Snapshots are saved in the directory: "/saved_snapshots/".

--- make_plot.py: Make the plot shown in Fig1D.

- EffectiveModel.nb: Contains the numerical analysis of the effective model leading to the boundaries shown in Fig3.

## Contribution guidelines
- Email: (shiladtb@andrew.cmu.edu)
  
## Who do I talk to?
- Fernanda Pérez-Verdugo (fverdugo@andrew.cmu.edu)
- Samuel Banks (sam.banks@yale.edu)
- Shiladitya Banerjee (shiladtb@andrew.cmu.edu)
