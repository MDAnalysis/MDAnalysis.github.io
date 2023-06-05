---
layout: post
title: Station1 Frontiers Fellowship Students 2023
---
<p>
<img
src="https://images.squarespace-cdn.com/content/v1/5992e994bebafb41ac4baa11/1512655046354-ZVQXATKKVTLFLLOX2JJB/icon++newwwwww2+%283%29.png"
title="Station1 Condensed Logo" alt="Station1 Condensed Logo"
style="float: right; height: 8em; " />
</p>

MDAnalysis is excited to be participating this year as a host organization for the [Station1 Frontiers Fellowship][sff] program for the first time. We thank Station1 for supporting two students, John Vinh Ong and Karen Bekhazi, who will be working with us to improve and extend MDAnalysis's interoperability. Stay tuned for updates on the project, and keep reading to learn more about John and Karen!

## John Vinh Ong

<img
src="John's Picture Here"
title="John Vinh Ong" alt="John Vinh Ong"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Write something about yourself here!

Add where we can find you here! (github, twitter, others)

Add a link to your blog here!

## Karen Bekhazi

<img
src="Karen's Picture Here"
title="Karen Bekhazi" alt="Karen Bekhazi"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Write something about yourself here!

Add where we can find you here! (github, twitter, others)

Add a link to your blog here!

## Project Summary

The development of MDAnalysis has always been driven by the growing need for standardized, accessible analysis tools for open, reproducible, and collaborative research. This project therefore consists of three main objectives: improving
the interoperability between MDAnalysis and the Protein Data Bank (PDB) format; extending MDAnalysis’ interoperability; and adding benchmarking cases to the MDAnalysis library to optimize performance. MDAnalysis allows users to read particle-based trajectories like the ones produced by MD simulations or individual coordinate frames - such as biomolecules in the PDB format - and access the atomic coordinates through NumPy arrays. This project will include fixing [open issues related to PDB formatting][PDBissues]. MDAnalysis is always aiming to not only expand its number of supported package-specific data formats, but also to [expand its direct interoperability][interoperability] with other popular packages for molecular analysis by becoming API compatible. We have already added converters for the ParmEd and RDKit libraries. We aim to expand the range of our converters framework to include packages in three categories: widely-used analysis libraries, such as [MDTraj][mdtraj] and [pytraj][pytraj]; libraries that can expand the range of supported formats, such as [OpenBabel][OpenBabel]; and direct interfaces with computational chemistry engines, such as [OpenMM][openmm] and [Psi4][psi4]. This project would thus give students the opportunity to create converter classes to and from MDAnalysis to their chosen package(s). The performance of the MDAnalysis library is assessed by [automated benchmarks][benchmarks] with [ASV][asv]. The goal of this project is to increase the performance assessment coverage and identify code that should be improved. Specifically, this project will involve writing benchmark cases; analyzing the performance history to identify code that needs to be improved; and optimizing the code for at least one of the discovered performance bottlenecks.

— @hmacdope @richardjgowers @jennaswa (mentors and org admin)

[sff]: https://www.station1.org/sff
[PDBissues]: https://github.com/MDAnalysis/mdanalysis/labels/Format-PDB
[interoperability]: https://www.mdanalysis.org/2020/08/03/roadmap/
[mdtraj]: https://www.mdtraj.org/1.9.8.dev0/index.html
[pytraj]: https://amber-md.github.io/pytraj/latest/index.html
[OpenBabel]: https://openbabel.org/wiki/Main_Page
[openmm]: https://openmm.org/
[psi4]: https://psicode.org/
[benchmarks]: https://www.mdanalysis.org/benchmarks/
[asv]: https://github.com/airspeed-velocity/asv/