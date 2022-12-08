---
layout: post
title: GSoC Report: "Adding a new guessing methodology for molecular systems properties in MDAnalysis"
---



## Project motivation

In [MDAnalysis], molecular topologies come from various file formats and forcefields, each of which has its own rules regarding formatting and molecular properties representation. Not all molecular properties are read from files, either because they are simply not provided in the file or because they can be inferred from the nature of the molecular system or forcefield definitions (e.g., for atomic systems masses can be inferred from atom type, while in the Martini forcefield there are only three definitions of masses for its known beads). MDAnalysis had different guessing methods that can infer missing molecular properties for the [Universe]'s topology (eg. [guessing Atom type, mass, and bond]). But those methods suffered from the main issue of being a general-prpose methods, that doesn't cover any kind of special definitions for the various forcefields and file formats, which causes the guessing output to be inaccurate and not reliable for some topologies; as not all topologies speak the same language. That’s why we need to develop tailored context-aware guesser classes to serve different molecular dynamics worlds, in addition, we need to establish a way to pass this context to the [Universe], so the [Universe] uses the appropriate guesser class in guessing topology attributes. 

So, my [GSoC project] was about developing a new guessing API, that will make the guessing process more convenient and context-specific, plus developing the first two context-specific guessers, which are [PDB] Guesser and             [Martini forcefield] Guesser.


### 1- Developing the new guesser API

The first thing to start with was the development of the new guesser API and the building of the [DefaultGuesser] class (which will represent a general-purpose guesser that has the current general-purpose guesser methods). But before that, we must first remove all the usage of the old guessing methods inside the MDAnalysis library to prepare for the new stage of having context-aware guessers and the guessing API. But this can lead to unexpected wrong behavior if we did without being careful. Guessing mainly takes place inside 16 topology parsers (`types` and `mass` guessing specifically), and not all guessing inside those parsers happens the same way. So, to not break the current behavior, we must be careful with removing all the guessing from those parsers and yet keep the same default output while developing the new API.

N.B.: it was important to have all the work on developing both the guesser API and the [DefaultGuesser] plus changing parsers behavior to be done inside one pull request here: [#3753]; because all thoses changes have a direct effect on each others, for example, the optmization of the new guesser API must be measured by how far the topology output changed after removing all the guessings processes from parsers. This was the most challenging thing about the project; as all these required changes are interconnected and must be done simultaneously, to see how it reflects and interacts with each other, and accordingly we can achieve optimization. So, working on this PR was not sequential as mentioned below, it was rather a continuous back-and-forth update on all touched parts of the library.

#### - Removing types and mass guessing from Topology Parsers ([commit 0cbc497]) 
First, I navigated through all parsers to see where and how types and mass guessing take place. The below table shows the behavior of types and masses guessing inside parsers.
![Guessing inside Parsers](public/images/final_report_aya/parsers_guessing.png)
As we see, not all parsers use the same attribute in guessing masses. The [DefaultGuesser] guesses masses by first looking for the elements attribute and if not available, it looks for the types attribute, this behavior preserves the default behavior of all parsers except the ones that guess masses from names. Parsers that guess masses from types are [FHIAIMSParser], [TXYZParser], and [XYZParser]. For both FHIAIMSParser and XYZParser `element` attributes are provided, and luckily it is just a copy of the `names` attribute, so the default behavior of these parsers is not affected. For TXYZParser, I added a feature of publishing Elements attribute to its topology output in case all read names are valid elements ([merged PR #3836]).

#### - Developing the `guess_topologyAttributes` API
As described in issue [# 3704], the new guesser API will have the tasks of guessing and adding new topology attributes to the [Universe]. The user only passes the context (eg. "default", "pdb", "martini") and the attributes of interest to be guessed, either through the `to_guess` or `force_guess` parameters and the API does the rest. Since I was writing my proposal, I was keen to develop the guesser API with this high level of abstraction and make guessing attributes easy and straightforward to the user without bothering about which method should be called and which attribute(s) is used as a reference in guessing other attributes, so all this work is handled inside the API and the guessers classes. The `guess_topologyAttributes()` make the following processes:

a- gets the appropriate Guesser class that matches the passed context.

b- checks if the Guesser class support guessing the attributes passed to the`to_guess` and/or `force_guess` parameters.

c- check if any attribute passed to the `to_guess` parameter already exists in the topology attributes of the [Universe]. If so, warn the user that only empty values will be filled for this attribute, if any exists and in case the user wishes to override all the attribute values, he must pass it to the `force_guess` parameter instead of the `to_guess` one.

d- guessing attributes is handled by the `guess_attrs method()`, which is declared inside the parent guesser class GuesserBase. It manages partial and complete guessing of attributes and calls the appropriate guesser method for each attribute.

e- each attribute guesser method searches for the reference attribute to begin guessing from it, and if not found in the [Universe], it calls the `guess_topologyAttributes`  to try guessing this reference attribute.

f- after guessing the attribute, the API adds it to the [Universe] with the help of `add_topologyAttr` or `_add_topology_objects` [Universe]’s methods

Example of using the `guess_topologyAttributes` at [Universe] initiation:

```python
# to guess bonds for a [Universe]:
 
 import MDAnalysis as mda
 from MDAnalysisTests.datafiles import two_water_gro
 
 u = mda.Universe(two_water_gro, context='default', to_guess=['bonds'])
```

Example of using the `guess_topologyAttributes` directly:

```python
# guess masses and types attribute::
 
u.guess_TopologyAttributes(context='default', to_guess=['masses', 'types'])
```
Example of passing empty `to_guess` list to `guess_topologyAttributes` so no guessing takes place at universe creation:

```python
# silencing masses and types guessing at universe creation::
 
 import MDAnalysis as mda
 from MDAnalysisTests.datafiles import two_water_gro
 
 u = mda.Universe(two_water_gro, to_guess=())
```
More explanation is found in the user guide here: [guess_topologyAttributes], [Guessing]

#### - [DefaultGuesser] class
The `DefaultGuesser` class holds the same old guessing methods but with modifying them to be compatible with the new guess_topologyAttrsibutes API. Moreover, I added a new feature to type guessing method so that it now can guess types from masses if names are not available ([commit 348f62d])

#### - Testing the old parser behavior is preserved
The `to_guess` parameters of the [Universe] have a default value of `("types", "masses")` to maintain the default behavior of the parsers. To make sure of not breaking old behavior, I added three tests in the parser's base.py module to check three things: 

a- types and masses are guessed as expected in all [Universe]'s created with the parsers after removing guessing from them.

b- The values of the guessed types with the `guess_topologyAttributes` API are the same as those from the old behavior. 

c- The values of the guessed masses with the `guess_topologyAttributes` API are the same as those from the old behavior. 
Once those three tests are passed, then it's safe to say that we are not breaking the default behavior of the code ([commit af84927]).

I also discovered a bug in Topology's methods `gussed_attributes()` and `read_attributes()`  while working on developing the guesser API and fixed it ([merged PR # 3779])


### 2- Working on PDBGuesser

Currently, I’m working on developing the PDBGuesser [issue #3856]. The generic guessing methods can’t deal optimally with PDB files, which makes guessing processes slow and not reliable for pdb-generated topologies. So, if we had a PDB-aware guesser, this process could improve significantly, especially that PDB has a huge archive called the [chemical component dictionary (CCD)], which describes every single residue that exists in the PDB database (its atom names, atom elements, bonds, bond orders, charges, aromaticity, etc.), plus that PDB has a well-formatted structure, that makes it easy to infer topology properties from it. I’m working on developing guesser methods for elements, masses, bonds, and aromaticity.

#### a. Elements guessing
PDB has a well-defined format for names, from which we can get the atomic symbol easily.
Atom names are found in columns 13-16. The first two characters represent the atom symbol, and if the symbol consists of one character, then the first character is blank. At the third character comes the remoteness indicator code for amino acid residues `['A', 'B', 'G', 'D', 'E', 'Z', 'H']`. Then the last character is a branching factor if needed.

The above rules are the standard rules but there are some exceptions to them:
- If the first character is blank and the second character is not recognized as an atomic symbol, we check if the third character contains "H", "C", "N", "O", "P" or "S", then it is considered the atomic symbol.

- If the first character is a digit, ”, ’, or  *, then the second character is the atomic symbol.

- If the first character in 'H' and the residues are a standard amino acid, nucleic acid, or known hetero groups (found in pdb_tables.py), then the atom element is 'H'.

- If the first two characters are not recognized as an atomic symbol and the first character is 'H', then the element is H.

Based on these rules, I developed the `guess_types` method for PDBGuesser [pr #3866].

#### b. Masses 
Mass guessing is the same as generic mass guessing methods, I just added a more detailed message about how many successful guessings happened and how many failed, in addition to which atom type/element the guesser failed to guess mass to.

PDBGuesser is still under discussion, so some of the current implementations may be updated in the future.


## Future work
I’m currently working on the PDBGuesser and plan to implement the MartiniGuesser after the GSoC period. Completing those context-aware guessers is not just important for having more tailored guessers, but also crucial for testing how the new guesser methodology is flexible and can work fine with different types of Guessers classes that will be implemented in the future.


## Lessons learned
GSoC was my first software internship, and I feel lucky that I got the chance to participate in such a wonderful program under the MDAnalysis organization, where everyone is supportive and friendly. I have learned about contributing to the open-source community, and the value of being a part of it. I also learned more about software engineering principles and how to better estimate the time and effort needed for your work. Additionally, I learned how to set priorities in developing new features, and the importance of test-driven software development. Finally, I learned how to effectively communicate and represent my work, and the importance of clean code and well-documented steps.

I’m happy with my experience with MDAnalysis and looking forward to increasing my contribution to the library, especially in the guessing and topology parts which I gained lots of experience at.

## Acknowledgements
I'd like to thank all my mentors for their effort and valubale lesson they gave to me through the program period, and I'm 
specially gratful for @jbarnoud (Jonathan Barnoud) for his endless guidance and patience through every step in the project.

-- @aya9aladdin


[MDAnalysis]: https://www.mdanalysis.org/
[guessing Atom type, mass, and bond]: https://userguide.mdanalysis.org/stable/formats/guessing.html
[GSoC project]: https://summerofcode.withgoogle.com/programs/2022/projects/B1Y0nTh2
[PDB]: https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction#:~:text=These%20files%20are%20available%20in,the%20atoms%20and%20their%20coordinates.
[Martini forcefield]: http://cgmartini.nl/index.php/tutorials
[FHIAIMSParser]: https://userguide.mdanalysis.org/stable/formats/reference/in.html#in-format
[TXYZParser]: https://userguide.mdanalysis.org/stable/formats/reference/txyz.html#txyz-format
[XYZParser]: https://userguide.mdanalysis.org/stable/formats/reference/xyz.html#xyz-format
[guess_topologyAttributes]: https://aya9aladdin.github.io/UserGuide/universe.html#guessing-topology-attributes
[Guessing]: https://aya9aladdin.github.io/UserGuide/formats/guessing.html
[chemical component dictionary (CCD)]: https://www.wwpdb.org/data/ccd
[Universe]: https://userguide.mdanalysis.org/stable/universe.html
[DefaultGuesser]: https://aya9aladdin.github.io/UserGuide/formats/guessers/default_guesser.html#default-guesser

[#3753]: https://github.com/MDAnalysis/mdanalysis/pull/3753
[commit 0cbc497]: https://github.com/MDAnalysis/mdanalysis/pull/3753/commits/0cbc4971bae924bd4d030022e5ff06faeb8c0aa6
[merged PR #3836]: https://github.com/MDAnalysis/mdanalysis/pull/3826
[# 3704]: https://github.com/MDAnalysis/mdanalysis/issues/3704
[commit 348f62d]: https://github.com/MDAnalysis/mdanalysis/pull/3753/commits/348f62dfb2707f10387a5ddab37426ccb7abba00
[commit af84927]: https://github.com/MDAnalysis/mdanalysis/pull/3753/commits/af849272dd120c04219fbd558d02c25afa24caf8
[merged PR # 3779]: https://github.com/MDAnalysis/mdanalysis/pull/3779
[issue #3856]: https://github.com/MDAnalysis/mdanalysis/issues/3856
[pr #3866]: https://github.com/MDAnalysis/mdanalysis/pull/3866
