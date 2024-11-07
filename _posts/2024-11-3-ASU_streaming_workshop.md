---
layout: post
title: "MDAnalysis Streaming Online Developer Workshop (Dec 4, 2024)"
---


Have you ever wanted to analyze sub-picosecond dynamics in your trajectories? Trajectory file sizes too large? Want to sync up your analysis and trajectory production? Lucky for you MDAnalysis, in conjunction with [Arizona State University][ASU] (ASU) and with the support of a [CSSI Elements][CSSI] grant from the [National Science Foundation][NSF], is holding a **free, online developer workshop** focused on streaming and inline analysis of molecular simulations on **December 4th** 2024.

The general idea of streaming, just like with Netflix, is to transfer data piece-by-piece as needed instead of transferring entire files. In our case, the data generated during a running simulation is transmitted to MDAnalysis for processing without ever being stored on disk.

Our streaming interface is built on top of the TCP/IP socket protocol and can transmit data between distinct processes on: A) the same computer; B) different computers in a local network; C) via the internet. This allows analyzing MD simulation trajectories live while they are being generated. As a result, the streaming interface allows analyzing data at femtosecond-scale time intervals which would create massive trajectories and slow down the simulation engine if written to disk.

This online workshop is intended to introduce participants to streaming of trajectories directly from simulation engines, inline analysis 
of simulations, and all the awesome science you can do with streaming. This workshop is suitable for students, developers, and researchers in the broad area of computational (bio)chemistry, materials science, and chemical engineering. It is designed for those who have some familiarity with MDAnalysis and are comfortable working with [Python](https://www.python.org/), [Jupyter
Notebooks](https://jupyter-notebook.readthedocs.io/en/stable/) and a molecular simulation engine such as [LAMMPS][LAMMPS], [GROMACS][GROMACS] or [NAMD][NAMD].



## Workshop Overview

The program will run from 8:00 am to 12:00 pm Pacific time on Wednesday, December 4th.
In the workshop, we will focus on contextualizing MD streaming, showing you some of its use cases from working as basic connective tissue to advanced, high-time-resolution analyses, and getting your hands dirty with streaming in a live-coding activity in an easy-to-use workshop environment.

| Topic | Duration |
| --- | --- |
| 👋 Welcome  | 5 min |
| 📦 MDAnalysis mission & ecosystem | 15 min |
| 🖼️ Streaming: big picture  | 15 min |
| 👀 Streaming: first look | 10 min |
| ❓ Q&A: Streaming overview  | 5 min |
| 📦Streaming: MD packages, IMDClient | 15 min |
| 👀 Demo: Multiple analyses on NAMD simulation stream | 10 min |
| 💤 Break | 10 min |
| 🎯Activity: Write your own stream analysis  | 40 min |
| 📦 Streaming: MDAnalysis functionality | 10 min |
| ❓Q&A: Streaming with MDAnalysis | 5 min |
| 👀 Application: Velocity correlation functions and 2PT | 10 min |
| 👀 Application: Ion channel permeation | 10 min |
| ❓ Q&A: Applications | 5 min |
| 🔮 Future direction | 5 min |
| 📖 Open Forum | 20 min |
| 🚪 Closing | 5 min |

## Registration

Attendance at this workshop will be *free*, and we encourage anyone with an interest in attending to register below. 

<a href="https://docs.google.com/forms/d/e/1FAIpQLSfSOmPEcV3uLBLFEo1EvQGPh1CwpWyKxChPZp_VSW9rNJLTgw/viewform" target="_blank" style="background:#FF9200;padding:10px;margin:10px 0px;text-align:center;text-decoration:none;font-size:12pt;color:#000000;display:inline-block;border-radius:3px">Register</a>



## Workshop materials
All materials are made available in the https://github.com/Becksteinlab/imd-workshop-2024 repository.

Prepare for the interactive workshop activities by following the [set-up instructions](https://github.com/Becksteinlab/imd-workshop-2024).

## Who to Contact

If you have any questions or special requests related to this workshop, you may [contact the organizing committee](mailto:workshops@mdanalysis.org).

- @hmacdope @yuxuanzhuang @IAlibay @jaclark5 (MDAnalysis) @orbeckst @ljwoods2 @HeydenLabASU @amruthesht @hcho38 (ASU)


[ASU]: https://www.asu.edu
[CSSI]: https://new.nsf.gov/funding/opportunities/cssi-cyberinfrastructure-sustained-scientific-innovation
[NSF]: https://www.nsf.gov/
[LAMMPS]: https://www.lammps.org/#gsc.tab=0
[GROMACS]: https://www.gromacs.org/
[NAMD]: https://www.ks.uiuc.edu/Research/namd/
