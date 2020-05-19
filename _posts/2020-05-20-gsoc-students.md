---
layout: post
title: Google Summer of Code Students 2020
---

We are happy to anounce that MDAnalysis is hosting three [GSoC][gsoc]
students this year -- @hmacdope, @cbouy, and @yuxuanzhuang. This is
the first year that MDAnalysis has been accepted as its own
[organization with GSoC][mda-gsoc] and we are grateful to Google for granting us
three student slots so that we can have _three exciting GSoC
projects_.

## Hugo MacDermott-Opeskin: Trajectory New Generation: the trajectory format for the future of simulation

<img
src="https://pbs.twimg.com/profile_images/1140166815919665152/D-x2Tvae_400x400.jpg"
title="Hugo MacDermott-Opeskin" alt="Hugo MacDermott-Opeskin"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Trajectory storage has always proved problematic for the molecular
simulation community, as large volumes of data can be generated
quickly. Traditional trajectory formats suffer from poor portability,
large file sizes and limited ability to include metadata relevant to
simulation. The Trajectory New Generation (TNG) format developed by
the [GROMACS][] team represents the first trajectory format with small
file sizes, metadata storage, archive integrity verification and
user/software signatures. The primary goal of [this
project](https://summerofcode.withgoogle.com/projects/#5116604104310784)
is for @hmacdope to refactor the existing TNG code into C++ to provide
clarity and usability for GROMACS, other simulation packages and
analysis tools. Thin FORTRAN and Python layers are also desirable to
encourage widespread adoption and are a secondary goal of the
project. An efficient and transferable implementation of the TNG
format will represent a major step forward for the computational
molecular sciences community, enabling easy storage and replication of
simulations.

This project is a collaboration with the [GROMACS][] developer team
with @acmnpv from GROMACS serving as a co-mentor.

Hugo MacDermott-Opeskin is a PhD student in computational chemistry at
the Australian National University. His work focuses on studying
membrane biophysics through molecular dynamics simulations coupled
with enhanced sampling techniques.  Hugo can be found on github as
@hmacdope and on twitter as [@hugomacdermott][hmacdope-twitter]. When
not hard at work Hugo can be found running or mountain biking in the
Canberra hills.

Through GSoC Hugo aims to bring the TNG next generation trajectory
format to the simulation community and he will document his experience at
[his "Biophysics Bonanza" blog][hmacdope-blog].


## Cédric Bouysset: From RDKit to the Universe and back

<img
src="https://cbouy.github.io/assets/img/photo-CV.jpeg"
title="Cédric Bouysset" alt="Cédric Bouysset"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

The aim of the [RDKit interoperability
project](https://summerofcode.withgoogle.com/projects/#6750913248624640)
is to give MDAnalysis the ability to use [RDKit][]’s `Chem.Mol`
structure as an input to an MDAnalysis `Universe`, but also to
convert a `Universe` or `AtomGroup` to an RDKit `molecule`. RDKit is
one of the most complete and one of the most commonly used
chemoinformatics package, yet it lacks file readers for formats
typically encountered in MD simulations. @cbouy will implement in
MDAnalysis the ability to switch back and forth between a `Universe`
and an RDKit `molecule` to perform typical chemoinformatics
calculations and so add a lot of value to both packages.

Cédric is a PhD student in molecular modelling at Université Côte
D'Azur, France. His research aims to decipher the molecular basis
of chemosensory perception (smell and taste) using computational
tools. His day-to-day work includes; modelling bitter taste receptors,
building machine-learning models to search for molecules with
interesting olfactive or sapid properties, maintaining the website
of the Global Consortium of Chemosensory Researchers, and a bit of
teaching. In his free time he enjoys cooking and playing video games. 

Cédric will describe his progress in [his blog][cbouy-blog].


## Yuxuan Zhuang: Serialize Universes for parallel

<img
src="https://i0.wp.com/www.biophysics.se/wp-content/uploads/2018/06/IMG_4767.jpg"
title="Yuxuan Zhuang" alt="Yuxuan Zhuang"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

MDAnalysis is a flexible and relatively fast framework for complex
analysis tasks in molecular dynamics (MD) simulations. To achieve a
flawless implementation of parallelism, @yuxuanzhuang will implement
[serialization support for
`Universe`](https://summerofcode.withgoogle.com/projects/#5812065073102848),
the core of MDAnalysis. Furthermore, he will adapt this new
serialization functionality into Dask, multiprocessing, or
MPI. Additionally, he will run tests, write documentation, and run
benchmarks.

Yuxuan is a PhD student at Stockholm University. He mainly works on
understanding pentameric ligand-gated ion channels from MD simulations.
His daily workflow involves setting up and running simulations,
on lab clusters or HPC centers, and performing various analyses on the
MD trajectories in his jupyter notebook.

Yuxuan will chronicle his work on [his blog
"Simulacrum"][yuxuanzhuang-blog].


— @richardjgowers @IAlibay @acmnpv @fiona-naughton @orbeckst (mentors)

[gsoc]: https://summerofcode.withgoogle.com
[mda-gsoc]: https://summerofcode.withgoogle.com/organizations/4891814374408192/
[GROMACS]: http://www.gromacs.org
[RDKit]: http://rdkit.org/
[hmacdope-blog]: https://hmacdope.github.io
[hmacdope-twitter]: https://twitter.com/hugomacdermott
[cbouy-blog]: https://cbouy.github.io/blog/
[yuxuanzhuang-blog]: https://yuxuanzhuang.github.io/
