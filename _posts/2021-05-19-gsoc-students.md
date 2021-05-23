---
layout: post
title: Google Summer of Code Students 2021
---

We are happy to announce that MDAnalysis is hosting two [GSoC][gsoc]
students this year -- @ojeda-e, and @orioncohen. MDAnalysis has been accepted as its own
[organization with GSoC][mda-gsoc] for a second year running and we are grateful to Google for granting us
two student slots for two exciting projects. Both the students and mentors have a very exciting few months ahead!

## Estefania Barreto-Ojeda: Curvature analysis of biological membranes 

<!-- <img
src="https://ojeda-e.github.io/assets/images/profile-photo.jpg"
title="Estefania Barreto-Ojeda" alt="Estefania Barreto-Ojeda"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" /> -->

Interested in contributing to an open-source initiative, Estefania will expand the capabilities of MDAnalysis by integrating a new [MDAnalysis module to calculate membrane curvature][ojeda-project] to derive and visualize membrane curvature profiles of protein-membrane/membrane-only systems obtained from Molecular Dynamics (MD) simulations. With the introduction of this analysis module, users will rapidly extract mean and gaussian curvature of biological membranes and their respective visualization in 2D-profile maps.

Estefania is a Ph.D. candidate in Biophysical Chemistry at The University of Calgary, Canada. Her research work is focused on membrane curvature induced by ABC transporters, a superfamily of transmembrane proteins involved in cancer and antibiotic resistance. A typical day for Estefania includes running Coarse-Grained (CG) MD simulations using the [Martini force field][martini-url], reading literature on ABCs, and working on cool data visualization workflows. In her free time, she enjoys camping and road tripping in the Canadian Rockie Mountains and going for long bike rides. 

Estefania can be found on github as [@ojeda-e][ojeda-git] and on twitter as [@ebojeda][ojeda-twitter].

Her journey will be documented on the blog [Le Mirroir][ojeda-blog].


## Orion Cohen: A Solvation Module for MDAnalysis

<!-- <img
src="https://perssongroup.lbl.gov/img/ocohen.jpg"
title="Orion Cohen" alt="Orion Cohen"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" /> -->

The macroscopic behavior of a liquid is determined by its microscopic
structure. For ionic systems, like batteries and many enzymes, the solvation
environment surrounding ions is especially important. By studying the solvation
of interesting materials, scientists can better understand, engineer, and
design new technologies. The aim of this project is to implement a robust
and cohesive set of methods for solvation analysis that would be widely
useful in both biomolecular and battery electrolyte simulations. The core of
the solvation module will be a cohesive set of functions for easily working
with ionic solvation shells. Building from that core functionality, the
module will implement several analysis methods for analyzing ion pairing,
ion speciation, residence times, and shell association and dissociation.

Orion is a Ph.D. student at the University of California Berkeley, working
with Dr. Kristin Persson at the Lawrence Berkeley National Laboratory. His
work leverages high-throughput chemical simulations and machine learning to
discover new materials for Lithium-ion batteries. Orion is passionate about
making science more accessible, reproducible, and efficient with powerful
open-source software. If he isn't toiling at his computer, Orion is probably
playing board games, camping, or relaxing in the temperate Berkeley sun.

— @richardjgowers @IAlibay @fiona-naughton @orbeckst @lilyminium @hmacdope (mentors)

[gsoc]: https://summerofcode.withgoogle.com
[mda-gsoc]: https://summerofcode.withgoogle.com/organizations/6414449348444160/
[martini-url]: http://cgmartini.nl/
[ojeda-git]: https://github.com/ojeda-e
[ojeda-blog]: https://ojeda-e.github.io/
[ojeda-twitter]: https://twitter.com/ebojeda
[ojeda-project]: https://summerofcode.withgoogle.com/projects/#5098282306502

