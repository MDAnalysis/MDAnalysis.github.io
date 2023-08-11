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

## John Ong

<img
src="https://jong9559.github.io/assets/img/prof_pic.jpg"
title="John Ong" alt="John Ong"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Hello everyone, my name is John Ong, and I am one of the Station1 interns at MDAnalysis this summer. I am a soon-to-be third-year undergraduate student at the University of Utah studying mathematics and computer science. Please forgive me for not having introduced myself earlier.

MDAnalysis is core to many professional applications, including treatments for various diseases. The Protein Data Bank (PDB) is an international repository for organic structures and contributes greatly to understanding and creating treatments for diseases. However, data representation techniques for structures on the PDB often differ, despite the same master format (PDBx/mmCIF). Hence, there exists "dialects" within the PDB's files. 

This summer, Karen and I explored the reasons for the existence of these dialects. Through literature and interviews we've carried out with both key members within the PDB and PDB user(s), we've found experimental methods, knowledge of proteins, data processing methods, and user error to be the most common contributors to these dialects. We also studied the motivations for format changes and regulations and the lack thereof within the PDB and its user base. 

All in all, this was a rewarding summer. While we did not achieve as many technical results as we would have liked, studying change has been interesting. I'm honored to have been a mentee of Richard Gowers, Ph.D., and Hugo MacDermott-Opeskin, Ph.D., and to have had the support of MDAnalysis and its community this summer. 

You can find me on [LinkedIn][john-linkedin] and you can find me on [my website][john-website]! 

I hope to be able to continue contributing to this incredible community!


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
[john-linkedin]: https://www.linkedin.com/in/john-ong/
[john-website]: https://jong9559.github.io/