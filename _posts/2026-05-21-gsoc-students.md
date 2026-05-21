---
layout: post
title: Google Summer of Code Students 2026
---

We’re excited to welcome _three_ contributors to MDAnalysis through [Google Summer of Code][gsoc] this year:
@jauy123, @kunjsinha, and @PardhavMaradani.

This marks our seventh consecutive year participating as an independent [GSoC organization][mda-gsoc].
A big thank you goes out to Google for supporting these *three* exciting projects.
We're looking forward to getting started!

## Kunj Sinha: [Interface for post-simulation analysis (“crawling”) of WESTPA simulations ](https://summerofcode.withgoogle.com/programs/2026/projects/VJMrp2UW)

<img
src="https://avatars.githubusercontent.com/kunjsinha"
title="Kunj Sinha" alt="Kunj Sinha"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

This project will implement WESTPAParser and WESTPAReader inside [westpa/westpa](https://github.com/westpa/westpa), exposing [WESTPA's HDF5 Framework](https://westpa.readthedocs.io/en/latest/users_guide/hdf5.html) simulation data as a standard [MDAnalysis Universe](https://userguide.mdanalysis.org/stable/universe.html). Post-simulation analysis currently requires custom boilerplate code via [w_crawl](https://westpa.readthedocs.io/en/latest/users_guide/command_line_tools/w_crawl.html) which this project will replace with a single command, making the entire MDAnalysis toolkit accessible on WESTPA data.

Kunj is an undergraduate student at PES University, pursuing a Bachelor of Technology in Computer Science and Engineering. He has always had an interest in various fields of science and technology since his early school days. In his free time, he listens to music, plays different musical instruments and enjoys cooking as well.

You can find Kunj on [Github](https://github.com/kunjsinha) and [LinkedIn](https://www.linkedin.com/in/kunjsinha).

To keep up with his work, you can check out his [blog](https://kunjsinha.github.io/blog).

## Pardhav Maradani: [Dashboard for tracking MD simulation progress with the new streaming interface](https://summerofcode.withgoogle.com/programs/2026/projects/DzTMshtu)

<img
src="https://avatars.githubusercontent.com/pardhavmaradani"
title="Pardhav Maradani" alt="Pardhav Maradani"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

This project will develop a browser-based, real-time dashboard for molecular dynamics simulations using the new IMDv3 streaming protocol, enabling live monitoring, analysis, and visualization of running simulations. The dashboard will leverage [imdclient](https://imdclient.readthedocs.io/) and MDAnalysis to receive and analyze data streams from running simulations.

Pardhav is an undergraduate student from India pursuing a Bachelors in Computer Science and Engineering from Vellore Institute of Technology (Vellore) and a BS in Data Science and Applications from Indian Institute of Technology (IIT) Madras.

You can find Pardhav on GitHub [@PardhavMaradani](https://github.com/PardhavMaradani).

To see updates on this project, you can check out his [blog](https://pardhavmaradani.github.io/categories/gsoc-2026/).

## Josh Uy: [Adding Additional Functionality and Enhancements to the Fetcher Module](https://summerofcode.withgoogle.com/programs/2026/projects/Jz9S7jxP)

<img
src=""{{site.images}}/joshuy.jpg"
title="Josh Uy" alt="Josh Uy"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

This project intends to add web retrieval functionality to MDAnalysis by augmenting the existing fetcher module. By introducing new fetchers,
it is intended to allow MDAnalysis to download structures and trajectories from external databases such as [AlphaFold](https://alphafold.ebi.ac.uk/) and the [Molecular Dynamics Data Bank](https://alphafold.ebi.ac.uk/) by passing downloaded data to the in-built [Universe](https://userguide.mdanalysis.org/stable/universe.html) class. Additionally, this project intends to introduce a new suite of functions that is capable of retrieving any molecular data from any database supporting REST API.

Joshua Raphael Uy is a Physics Ph.D student at Arizona State University. He earned a Bachelor of Science in Physics, a Bachelor of Arts in
Mathematics at Miami University in Oxford, Ohio, and a Master in Science in Physics at Arizona State University. He has a [bio page](https://becksteinlab.physics.asu.edu/people/155/joshua-raphael-uy) and a [Linkedin](https://www.linkedin.com/in/joshua-raphael-uy-9540a0202).

To follow this project, you can check out Josh's [blog](https://jauy123.github.io/) 

— @amruthesht @BradyAJohnston @HeydenLabASU @jeremyleung521 @ltchong @nilay-v3rma @orbeckst @talagayev @yuxuanzhuang @IAlibay @jennaswa (@MDAnalysis/gsoc-mentors and org admins)

[gsoc]: https://summerofcode.withgoogle.com
[mda-gsoc]: https://summerofcode.withgoogle.com/programs/2026/organizations/mdanalysis
