---
layout: post
title: "GSoC Report: From RDKit to the Universe and back 🚀"
---

With the end of the summer comes the end of my awesome Google Summer of Code adventure with the MDAnalysis team. It's also the occasion for me to report on the code I've implemented, the things I've learned and the challenges I've faced.

## Summary of the project
---

The goal of my project, "[From RDKit to the Universe and back][project-url]", was to provide interoperability between MDAnalysis and RDKit. *i.e.* to be able to:
* Read an RDKit molecule and create an MDAnalysis `Universe` with it,
* Convert an MDAnalysis `AtomGroup` to an RDKit `Molecule`,
* Leverage RDKit's functionalities directly from MDAnalysis (descriptors, fingerprints, aromaticity perception...*etc.*)

With this in mind, the project was easily cut down in 2 main deliverables (the RDKit Reader/Parser and the Converter) and a few smaller ones for each "wrapper" functionality.


## Contributions
---

### Merged PRs

* **Converting an RDKit molecule to MDAnalysis**: [#2707](https://github.com/MDAnalysis/mdanalysis/pull/2707)  
  This was the first part of the project and it helped me get acquainted with MDAnalysis as I had never used it before GSoC. The goal here was to be able to "parse" an RDKit molecule and build an MDAnalysis Universe from it.  
  It was also the occasion for me to start interfacing MDAnalysis with RDKit as I got to implement a new classmethod for the Universe, the `Universe.from_smiles` method which allows us to build a Universe from a SMILES string, but also to add atom selection based on aromaticity.  
  At the end of this first PR, I was familiar with the different core objects that compose the Universe as well as writing tests with pytest and documentation with sphinx.

* **Converting an MDAnalysis universe to RDKit**: [#2775](https://github.com/MDAnalysis/mdanalysis/pull/2775)  
  The second goal of the project was to be able to convert an MDAnalysis AtomGroup to an RDKit molecule, allowing users to analyse MD trajectories in RDKit. While this should be trivial once we know how to do the opposite operation, it was actually a real challenge to get a molecule with the correct bond order and charges out of it.  
  As you may know, most MD topology file formats don't keep track of bond orders and formal charges so we have to find a way to infer this information from what we have. In our case, we require all hydrogen atoms to be explicit in the topology file, as well as elements and bonds (although these two can be guessed). Then the bond orders and formal charges are inferred based on atomic valencies and the number of unpaired electrons, followed by a standardization step of functional groups and conjugated systems. This last step is needed because the algorithm implemented to guess bond orders and charges is dependent on the order in which atoms are read.  
  This took a more time than originally planned in the project timeline but was well worth it.

### In progress

* **SMARTS selection**: [#2883](https://github.com/MDAnalysis/mdanalysis/pull/2883/files/55691558811a1091bd8ff7a1868bbde140f63d72)  
  [SMARTS][smarts-url] is an extension of the SMILES language that is used for substructure searching. Being able to select atoms based on SMARTS queries, and combine these selections with those already available in MDAnalysis might be one of the key features that will come out of this project.  
  This is currently missing a few tests and clarification in the docs.

* **Wrap RDKit drawing code for AtomGroups**: [#2900](https://github.com/MDAnalysis/mdanalysis/pull/2900/files/1f9d3d04c14643a54aaf0cc89511f2811b51f35c)  
  This PR allows us to draw images (SVG, PNG, and GIF) of AtomGroups using RDKit. It also adds rich displays to AtomGroups in notebooks (a.k.a. `__repr__` methods).
  Before working on this, I thought the only way to use alternative representations for Python objects was to define a `_repr_*_` method for your class, where * is a MIME type such as png, html...*etc.* There is actually a second way, where you tell IPython directly how it's supposed to represent an object. This allows funky representations of any object, even python built-in types, *i.e.* `int` as roman numerals and so on. I also wrote my first metaclass here, to register different "viewer" classes when more become available in the future. It's also not straightforward to write tests for the images as different versions of RDKit or other packages will lead to slightly different outputs.
  More discussion, code review, and tests are needed before this is ready.

* **Wrap RDKit descriptors and fingerprints**: [#2912](https://github.com/MDAnalysis/mdanalysis/pull/2912/files/76a51ece7b676a5f5ba3d683fdf942706f87c141)  
  This PR adds new kinds of analysis that are typically performed on small molecules in the chemoinformatics field. Fingerprints are mostly used for calculating similarity metrics between molecules and a reference. Molecular descriptors could be used to describe all the sampled conformations of a ligand in a binding pocket during a simulation, and given as input to a machine-learning model for clustering, scoring binding poses...*etc.*  
  This is currently missing more discussion, documentation, and a few tests.

* **Change the Converters API**: [#2882](https://github.com/MDAnalysis/mdanalysis/pull/2882/files/43f6b34a3f7cfe57ed5fa296dbc800f304e95198)   
  While developing the RDKitConverter, some interesting points were made about the current API used to convert AtomGroups, *i.e.* `u.atoms.convert_to("RDKIT")`. This method is case sensitive, it doesn't allow to pass arguments to the converter class, and it requires users to read the documentation to know which converters are available. This PR corrects the two first points and adds the possibility to either use the previous syntax or tab-completion to find the available converters *i.e.* `u.atoms.convert_to.rdkit()`.  
  This was inspired by pandas `df.plot(kind="scatter", ...)` which is also accessible as `df.plot.scatter(...)`.  
  Since this is an API change, more discussion with core developers is needed for now.

### Left to do

* Guessers for aromaticity and Gasteiger charges through RDKit
* Tutorial on the reader/converter and wrapped RDKit functionalities in the UserGuide. This will make it easier for users to know about the features I implemented and how to use them properly, rather than searching for every single feature in the documentation.
* Documentation on the RDKit format in the UserGuide: [#69](https://github.com/MDAnalysis/UserGuide/pull/69)


## Conclusion
---

This project taught me a lot of things on software development and the Google Summer of Code experience has been an incredible and valuable to me. All of this wouldn't have been possible without the help of many people, including my amazing mentors (@IAlibay, @fiona-naughton and @richardjgowers), but also the rest of the MDAnalysis team as a lot of them got involved and gave me great feedback, so thank you to all of them!


— @cbouy

[project-url]: https://summerofcode.withgoogle.com/projects/#6750913248624640
[smarts-url]: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
