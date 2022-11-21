---
layout: post
title: GSoC 2022 - Adding Energy Readers to MDAnalysis
---

## Motivation
In molecular dynamics simulations, users frequently have to inspect energy-like terms such as potential or kinetic energy, temperature, or pressure. This is so common a task that even small inefficiencies add up. Currently, users have to create intermediate files from their MD simulation’s output files to obtain plot-able data, and this quickly becomes cumbersome when multiple terms are to be inspected. Being able to read in the energy output files directly would make this more convenient.

Therefore, I wanted to add readers for energy-type files (output files containing information on potential and kinetic energy, temperature, pressure, and other such terms) from a number of MD engines to the auxiliary module of MDAnalysis in this project. This would make quality control of MD simulations much more convenient, and allow users to analyse the energy data without the need for switching windows or writing intermediate files directly from within their scripts or jupyter notebooks.

In a first instance, I focussed on a reader for [EDR files], which are energy files written by
GROMACS during simulations. EDR files are binary files which follow the [XDR protocol]. To read these
files, @jbarnoud had previously written the [panedr] Python package, which was the
foundation of my work this summer.

## Adapting Panedr for use in MDAnalysis
The panedr package makes use of the [xdrlib] Python module to parse EDR files
and return the data in the form of a [pandas] DataFrame. My GSoC project started out
adapting this package for use in MDAnalysis. In particular, we wanted to avoid making
pandas a dependency in MDAnalyis. This necessitated some refactoring of panedr ([PR #33]),
which ultimately led to a restructuring of the code into two distinct packages: [panedr
and pyedr] (PRs [#42] and [#50]). Both packages read EDR files, but one returns the
data as a pandas DataFrame, the other as a dictionary of NumPy arrays. Both also
expose a function to return a dictionary of units of the energy terms found in the file ([PR #56]).

Example:

```python
import pyedr
file = "path/to/edr/file.edr"
energy_dictionary = pyedr.edr_to_dict(file)
unit_dictionary = pyedr.get_unit_dictionary(file)
```


## EDRReader
With Pyedr available, I started work on the implementation of an EDRReader in MDAnalyis ([PR #3749]).
Here, I benefited hugely from the existing [AuxReader framework].
However, while working on the reader, it quickly became apparent that the auxiliary API would need to be changed to accommodate
the large number of terms found in EDR files. While the [XVGReader] still works as previously, the new base case for adding auxiliary data assumes a dictionary to be passed. The dictionary maps the name to be used in MDAnalysis to the names read from the EDR file. This is shown in the following minimal working example:

```python
import MDAnalysis as mda
from MDAnalysisTests.datafiles import AUX_EDR, AUX_EDR_TPR, AUX_EDR_XTC
term_dict = {"temp": "Temperature", "epot": "Potential"}
aux = mda.auxiliary.EDR.EDRReader(AUX_EDR)
u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
u.trajectory.add_auxiliary(term_dict, aux)
```

Aside from this API change, the EDRReader can do everything the XVGReader can. In addition to that, it
has some new functionality.
* Because EDR files can become reasonably large, a memory warning will be issued when more than a gigabyte of storage is used by the auxiliary data. This default value of 1 GB can be changed by passing a value as `memory_limit` when creating the EDRReader object.
* EDR files store data of a large number of different quantities, so it is important to know their units as well. The EDRReader therefore has a `unit_dict` attribute that contains this information. By default, units found in the EDR file will be converted to [MDAnalysis base units] on reading. This can be disabled by setting `convert_units` to False on creation of the reader.
* In addition to associating data with trajectories, the EDRReader can also return the NumPy arrays of selected data, which is useful for plotting, for example. This is done via the EDRReader's `get_data` method.

Additionally, the new auxiliary readers allow the selection of frames based on the
values of the auxiliary data. For example, it is possible to select only frames
with a potential energy below a certain threshold as follows:
```python
u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
term_dict = {"epot": "Potential"}
u.trajectory.add_auxiliary(term_dict, aux)
selected_frames = np.array([ts.frame for ts in u.trajectory if ts.aux.epot < -524600])
```
More details on the EDRReader's functionality can be found in the [MDAnalysis User Guide].

## Outlook
Through this project, the AuxReader framework was expanded, and handling of EDR
files was made more convenient with pyedr and EDRReaders. I am continually making
improvements to these contributions, and will include an auxiliary reader for
NumPy arrays in the future. This NumPyReader will be very useful, because many
analysis methods in MDAnalysis return their results in the form of NumPy arrays.
Having the option of associating these results with trajectories will facilitate
further analyses, for example allowing the slicing of trajectories by RMSD to a reference
structure.


## Lessons learned
Participating in the Summer of Code was a great opportunity for me. I learned a lot, from small things like individual code patterns to larger points concerning overall best practices, the value of test-driven development, and package management. This is thanks in large part
to the mentorship and advice I have received from @hmacdope, @ialibay, @orbeckst, and @fiona-naughton.
Thanks very much to you all, and to @jbarnoud.


-- @bfedder



---
[EDR files]: https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#edr
[XDR protocol]: https://docs.oracle.com/cd/E19683-01/816-1435/xdrproto-1/index.html
[panedr]: https://github.com/mdanalysis/panedr
[xdrlib]: https://docs.python.org/3/library/xdrlib.html
[pandas]: https://pandas.pydata.org
[panedr and pyedr]: https://pypi.org/project/panedr/
[PR #33]: https://github.com/MDAnalysis/panedr/pull/33
[#42]: https://github.com/MDAnalysis/panedr/pull/42
[#50]: https://github.com/MDAnalysis/panedr/pull/50
[PR #56]: https://github.com/MDAnalysis/panedr/pull/56
[PR #3749]: https://github.com/MDAnalysis/mdanalysis/pull/3749
[AuxReader framework]: https://userguide.mdanalysis.org/stable/formats/auxiliary.html
[XVGReader]: https://docs.mdanalysis.org/stable/documentation_pages/auxiliary/XVG.html
[MDAnalysis base units]: https://docs.mdanalysis.org/stable/documentation_pages/units.html#id68
[MDAnalysis User Guide]: https://userguide.mdanalysis.org/2.4.0-dev0/formats/auxiliary.html#edr-files
