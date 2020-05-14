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
simulation. The Trajectory New Generation format developed by the
[Gromacs][] team represents the first trajectory format with small
file sizes, metadata storage, archive integrity verification and
user/software signatures. In [this
project](https://summerofcode.withgoogle.com/projects/#5116604104310784),
@hmacdope will refactor this code into C++ to provide clarity and
usability for both the Gromacs project and for other simulation
packages and analysis tools and is the primary goal of this
project. Thin FORTRAN and Python layers are also desirable to
encourage widespread adoption and are a secondary goal of the
project. An efficient and transferable implementation of the TNG
format will represent a major step forward for the computational
molecular sciences community, enabling easy storage and replication of
simulations.

This project is a collaboration with the [Gromacs][] developer team
with @acmnpv from Gromacs serving as a co-mentor.

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
structure as an input to an MDAnalysis’ `Universe`, but also to
convert a `Universe` or `AtomGroup` to an RDKit `molecule`. RDKit is
one of the most complete and one of the most commonly used
chemoinformatics package, yet it lacks file readers for formats
typically encountered in MD simulations. @cbouy will implement in
MDAnalysis the ability to switch back and forth between a `Universe`
and an RDKit `molecule` to perform typical chemoinformatics
calculations and so add a lot of value to both packages.

Cédric is a PhD student at the Institute of Chemistry of Nice, in the
south of France. His research project tackles taste perception (and
sometimes olfaction) using computational approaches. His day-to-day
work includes homology modeling of bitter taste receptors, building
machine-learning models to predict taste modalities from chemical
structures, writing Python or Bash scripts, and a bit of teaching.

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

> BIOGRAPHY HERE

Yuxuan will chronicle his work on [his blog
"Simulacrum"][yuxuanzhuang-blog].


— @richardjgowers @IAlibay @acmnpv @fiona-naughton @orbeckst (mentors)

[gsoc]: https://summerofcode.withgoogle.com
[mda-gsoc]: https://summerofcode.withgoogle.com/organizations/4891814374408192/
[Gromacs]: https://www.gromacs.org
[RDKit]: http://rdkit.org/
[hmacdope-gsoc]: https://summerofcode.withgoogle.com/projects/
[hmacdope-blog]: https://hmacdope.github.io
[hmacdope-twitter]: https://twitter.com/hugomacdermott
[cbouy-gsoc]: https://summerofcode.withgoogle.com/projects/
[cbouy-blog]: https://cbouy.github.io/blog/
[yuxuanzhuang-gsoc]: https://summerofcode.withgoogle.com/projects/
[yuxuanzhuang-blog]: http://wsygzyx.com/
