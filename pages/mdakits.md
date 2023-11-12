---
layout: page
title: MDAKits and MDA-based tools
---

MDAnalysis is developed with extensibility in mind, allowing
scientists to create various tools with the components
that the library provides. To help developers in this process,
and to advertise these wonderful tools to the MDAnalysis user community,
we have created an ecosystem of [MDAnalysis Toolkits, or "MDAKits"](https://mdakits.mdanalysis.org/about.html).

## MDAKits 

Please look at the [registry of MDAKits](https://mdakits.mdanalysis.org/mdakits.html)
for a list of toolkits that [meet the minimum MDAKit requirements](https://mdakits.mdanalysis.org/about.html#requirements).
These are regularly tested against the latest and development versions of
MDAnalysis. Whilst we cannot guarantee the completeness of the MDAKit's test coverage, if you see a green badge, the code should be minimally behaving as the MDAKit authors intended it!


*If you have your own MDAKit and you would like to add it to the
registry, please look at the [contribution instructions](https://mdakits.mdanalysis.org/add.html)
or share it with on our [mailing list]({{ site.mailinglists.discussion.url }})!*


## Other (non-MDAKit registered) tools using MDAnalysis

Below we list projects that use MDAnalysis and are not registered as MDAKits:


### Visualization tools

-  [nglview](https://github.com/arose/nglview): nglview is a tool to visualize
   trajectories in jupyter notebooks.
-  [MDSrv](https://github.com/nglviewer/mdsrv) streams and visualizes MD trajectories interactively within web browsers
-  [MD Trajectories in
   PyMOL](https://nms.kcl.ac.uk/lorenz.lab/wp/?p=1768): MDAnalysis has
   been embedded into Open Source PyMOL to read many different MD formats
   directly. Source code is available from 
   [bieniekmateusz/pymol-mdanalysis](https://github.com/bieniekmateusz/pymol-mdanalysis) on GitHub.
-  [MolecularNodes](https://github.com/BradyAJohnston/MolecularNodes) plugin for the [Blender](https://www.blender.org/) rendering package; the plugin provides a convenient method for importing structural biology files, including MD trajectories, into Blender, and several nodes for working with atomic data inside of Blender's Geometry Nodes.

### Analysis tools

-  [pydiffusion](https://github.com/bio-phys/pydiffusion): Analyze the
   rotational diffusion of your molecules.
-  [pytim](https://marcello-sega.github.io/pytim/): Pytim is a package based on
   MDAnalysis for the identification and analysis of surface molecules in
   configuration files or in trajectories from molecular dynamics simulations.
-  [pycontact](https://github.com/maxscheurer/pycontact): Analysis of
   non-covalent interactions in MD trajectories.
-  [RotamerConvolveMD](https://github.com/MDAnalysis/RotamerConvolveMD):
   Analysis of molecular dynamics trajectories or conformational ensembles in
   terms of spin-label distances as probed in double electron-electron resonance
   (DEER) experiments.
-  [PBxplore](https://github.com/pierrepo/PBxplore): PBxplore is a suite of
   tools dedicated to Protein Block (PB) analysis.
-  [cgheliparm](https://github.com/ifaust83/cgheliparm): Scripts used to analyze
   dsDNA structures from Martini MD simulations.
-  [accelerated_sampling_with_autoencoder](https://github.com/weiHelloWorld/accelerated_sampling_with_autoencoder):
   This is the framework for running accelerated sampling with data-augmented
   autoencoders.
-  [BioEn](https://github.com/bio-phys/BioEn): BioEn integrates a broad range of experimental data to refine ensembles of structures.
-  [kugupu](https://github.com/kugupu/kugupu): A molecular network generator to study charge transport pathways in amorphous materials
-  [PyInteraph](https://github.com/ELELAB/pyinteraph): A software tool
   for the analysis of structural communication in protein ensembles,
   including a PyMOL plugin and an InteractionPlotter.
-  [MAICoS](https://gitlab.com/netzlab/maicos): Analyze molecular dynamics simulations of 
   interfacial and confined systems.
-  [taurenmd](https://taurenmd.readthedocs.io/en/latest/): A command-line interface for analysis of Molecular Dynamics simulations.
-  [PENSA](https://github.com/drorlab/pensa): A toolkit for exploratory analysis and comparison of protein structural ensembles
-  [LiPyphilic](https://lipyphilic.readthedocs.io/en/latest/): A Python package for the analysis of lipid membrane simulations.
-  [MDVoxelSegmentation](https://github.com/marrink-lab/MDVoxelSegmentation): A voxel-based approach for dynamic cluster analysis of molecular dynamics trajectories. 

### Molecular modeling tools

-  [Swarm-CG](https://github.com/GMPavanLab/Swarm-CG): Automatically optimizes the bonded terms of a MARTINI-like coarse-grained (CG) molecular model with respect to its reference all-atom (AA) trajectory, via [FST-PSO](https://github.com/aresio/fst-pso).

### Simulation packages
-  [ESPResSo](http://espressomd.org/) is a software package for
   performing and analyzing Molecular Dynamics many-particle
   simulations of coarse-grained atomistic or bead-spring models as
   they are used in soft matter research in physics, chemistry and
   molecular biology.

### Distributions

MDAnalysis is included in [NMRBox](https://nmrbox.org/) a distribution of
common software to analyze NMR measurements. 

Some generic distributions such as [MacPorts:
py-MDAnalysis](https://ports.macports.org/port/py-MDAnalysis/summary)
(macOS), [Fedora:
python-MDAnalysis](https://src.fedoraproject.org/rpms/python-MDAnalysis/)
and [archlinux:
python-mdanalysis](https://aur.archlinux.org/packages/python-mdanalysis/)
(Linux) maintain packages of MDAnalysis.


