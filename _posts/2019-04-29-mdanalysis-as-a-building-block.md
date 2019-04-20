---
layout: post
title: MDAnalysis as a building block
---

While our documentation is mostly focused on using MDAnalysis for exploratory
analysis it is equally well suited to build your own analysis library on top of
it. Below is a list of all projects we know about that use MDAnalysis.

# Visualization tools

1. [nglview](https://github.com/arose/nglview): nglview is a tool to visualize
   trajectories in jupyter notebooks.
2. [mda-pymol](https://nms.kcl.ac.uk/lorenz.lab/wp/index.php/2019/04/17/mda-pymol/):
   MDAnalysis has been embedded into PyMOL to read many different MD formats
   directly

# Analysis tools

1. [pydiffusion](https://github.com/bio-phys/pydiffusion): Analyze the
   rotational diffusion of your molecules.
2. [pytim](https://marcello-sega.github.io/pytim/): Pytim is a package based on
   MDAnalysis for the identification and analysis of surface molecules in
   configuration files or in trajectories from molecular dynamics simulations.
3. [pycontact](https://github.com/maxscheurer/pycontact): Analysis of
   non-covalent interactions in MD trajectories.
4. [pyPcazip](https://www.sciencedirect.com/science/article/pii/S2352711016300036):
   A PCA-based toolkit for compression and analysis of molecular simulation data
5. [RotamerConvolveMD](https://github.com/MDAnalysis/RotamerConvolveMD):
   Analysis of molecular dynamics trajectories or conformational ensembles in
   terms of spin-label distances as probed in double electron-electron resonance
   (DEER) experiments.
6. [PBxplore](https://github.com/pierrepo/PBxplore): PBxplore is a suite of
   tools dedicated to Protein Block (PB) analysis.
7. [cgheliparm](https://github.com/ifaust83/cgheliparm): Scripts used to analyze
   dsDNA structures from Martini MD simulations.
8. [accelerated_sampling_with_autoencoder](https://github.com/weiHelloWorld/accelerated_sampling_with_autoencoder):
   This is the framework for running accelerated sampling with data-augmented
   autoencoders.

# Distributions

MDAnalysis is also included in [NMRBox](https://nmrbox.org/) a distribution of
common software to analyze NMR measurements.


If you know of other tools that are build on MDAnalysis please share them with
us on twitter.
