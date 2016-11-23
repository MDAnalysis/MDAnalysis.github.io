---
layout: post
title: Support for MMTF has arrived in MDAnalysis!
---

The upcoming 0.16.0 release of MDAnalysis will have
[support for MMTF](https://twitter.com/mmtf_spec/status/799704395046760448)!
[MMTF](http://mmtf.rcsb.org/) is a new format designed to provide compact, efficient and fast browsing of the [Protein
Data Bank](www.rcsb.org/).
Support for MMTF within MDAnalysis is offered through reading locally stored MMTF files,
or through fetching the file directly from the PDB archive through the new
`MDAnalysis.fetch_mmtf` function.

```python
import MDAnalysis as mda

# Load a local file
u = mda.Universe('myfile.mmtf')

# Or download directly from PDB by providing PBD id
u = mda.fetch_mmtf('3J3Q')

```

The performance of loading MMTF files is a large improvement over traditional ascii PDB files,
with the above system of approximately 2.4M atoms taking under 10 seconds to load.
The compressed format and [efficient algorithms](https://github.com/rcsb/mmtf/blob/v1.0/spec.md)
for storing the data mean that downloading structures will also require much less bandwidth,
making this possible even on slow connections.

MMTF files can support many different models for a given structure and this is made available
through the `.models` attribute of a MDAnalysis Universe.  This provides a list of AtomGroup
objects each representing a different model.  These models are able to each have a
different topology.

```python
from __future__ import print_function
import MDAnalysis as mda

u = mda.fetch_mmtf('4P3R')

print("This file has {} models".format(len(u.models)))

# Iterate over all models
for model in u.models:
    # analyse each model!

# Select atoms in a given molecule
ag = u.select_atoms('model 4 and name Ca')
```

Finally, full interoperability between different formats is provided in MDAnalysis, allowing
MMTF files to be written to any of our
[supported formats](http://pythonhosted.org/MDAnalysis/documentation_pages/coordinates/init.html#id1).
For example to download a MMTF file and write out to Gromacs GRO file:

```python
import MDAnalysis as mda

u = mda.fetch_mmtf('4AKE')

u.atoms.write('4ake.gro')

```

— @richardjgowers
