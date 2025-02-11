---
layout: post
title: Luna Morrow GSoC 2024 Final Submission Blog
---

# Luna Morrow GSoC 2024 Final Submission Blog

I can’t believe I am at the end of [GSoC](https://summerofcode.withgoogle.com/)! This past 6 months has absolutely flown by, but I have learnt so much. I am really grateful for the amazing support I have had along the way from my mentors Hugo ( @hmacdope ), Cédric ( @cbouy ) and Xu ( @xhgchen ). My body of work can be located on the MDAnalysis GitHub at [mda-openbabel-converter](https://github.com/MDAnalysis/mda-openbabel-converter).


# Why a converter

Direct interoperability between molecular dynamics software is critical for enabling collaboration, data transfer and straightforward use by scientists. [OpenBabel](https://openbabel.org/) is a popular toolbox for chemical molecular modelling research as it enables searching, conversions, analysis and data storage. The ability to interconvert chemical file formats with OpenBabel, in particular, opens up the ability to utilise and work with other packages, as OpenBabel enables input and writing of over 100 chemical data file formats. Therefore, enabling an MDAnalysis Universe to be inter-converted with an OpenBabel OBMol would greatly increase the data formats available to MDAnalysis. Furthermore, it would encourage greater adoption of MDAnalysis as an “all-in-one” package for molecular dynamics analysis.


# From OpenBabel to MDAnalysis

The first aim of my GSoC project was to convert OpenBabel OBMols to MDAnalysis Universes. This required the extraction of both molecule and atom data, alongside their 3D positions. I hit many roadblocks in this section, with the two most difficult being the installation of OpenBabel 3.1.1.1 Python bindings and understanding some of the more complex OpenBabel API. Around half of the coding period was spent writing the two classes for converting the atom attributes and positions respectively. Time was also committed to developing tests for both classes.

First, I created the OpenBabelParser Class (converts the atom attributes) and some basic unit tests for it in [#12](https://github.com/MDAnalysis/mda-openbabel-converter/pull/12). Once this PR was opened, it revealed that the CI was not correctly setup. With additional help from Irfan ( @IAlibay ), this was corrected in [#13](https://github.com/MDAnalysis/mda-openbabel-converter/pull/13). Next, I made the OpenBabelReader Class (converts the atom positions) and a whole suite of tests for it in [#16](https://github.com/MDAnalysis/mda-openbabel-converter/pull/16), and then applied some of these ‘extended’ test conditions to the OpenBabelParser in [#18](https://github.com/MDAnalysis/mda-openbabel-converter/pull/18). During this time some issues were found with compatibility between, and accessing of, attributes from OpenBabel. I opened an issue on the OpenBabel GitHub to obtain assistance (see [here](https://github.com/openbabel/openbabel/issues/2708)) and also made a issue to flag an attribute that was proving difficult but will need to be supported later in [#17](https://github.com/MDAnalysis/mda-openbabel-converter/issues/17). 


# Documentation

Once the OpenBabelParser and OpenBabelReader were developed, I setup ReadTheDocs [#19](https://github.com/MDAnalysis/mda-openbabel-converter/pull/19) and wrote documentation for these two classes [#20](https://github.com/MDAnalysis/mda-openbabel-converter/pull/20). The documentation is readily available on ReadTheDocs [here](https://mda-openbabel-converter.readthedocs.io/en/latest/). The documentation is functional but currently quite sparse, I have created an issue to improve the landing page and getting started page, once functionality is complete, at [#21](https://github.com/MDAnalysis/mda-openbabel-converter/issues/21).


# From MDAnalysis to OpenBabel

The next stage of this this project was to implement the OpenBabelConverter Class to convert an MDAnalysis Universe to OpenBabel OBMol. Unfortunately, I have run out of time to implement this class during the allotted time for GSoC. There is an open issue detailing this class at [#22](https://github.com/MDAnalysis/mda-openbabel-converter/issues/22), which I will be expanding on shortly. While my time with GSoC has ended, I am keen to complete this project in my own time and continue on with MDAnalysis as a developer.


# Next Steps

I will be finishing the OpenBabelConverter, including tests and documentation. The next steps will then be to add attributes that were left for later, further development of the documentation with worked examples and deploying a release of this package. I would also like to integrate this converter into the MDAnalysis package so that it can be used alongside the other converters, as this would increase its visibility and usage.


# What can we do now?

With what I currently implemented at the end of the GSoC coding period, users can very easily convert an OBMol to an MDAnalysis universe.

An example of how to use the converter to convert an OBMol to an MDAnalysis universe is shown below:



```python

from openbabel import openbabel as ob
import MDAnalysis as mda
obconversion = ob.OBConversion()
obconversion.SetInFormat("pdb")
mol = ob.OBMol()
obconversion.ReadFile(mol, "1crn.pdb")
u = mda.Universe(mol)

```

# What I’ve Learned

Things I have learned during GSoC include:
* the importance of good documentation
* how to develop the ‘backbone’ of a python package so it can be installed, tested and used
* developing tests with pytest
* managing myself and being able to make decisions about when something is out of scope or unviable to implement
* how to use git and GitHub
* developing classes that inherit from classes that were not developed by me


# My experience

I have had an amazing time during GSoC. I had a lot of support and felt very welcomed and encouraged. I can’t wait for this converter to be up and running, so that the community can benefit from it. I am also really grateful for the experience and technical growth GSoC granted me, as I know it will be incredibly beneficial for my future in this field. 