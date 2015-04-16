---
layout: otherpage
title: Paper Errata
---

The paper on MDAnalysis was published as

  * N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
    O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular
    Dynamics Simulations. _J. Comput. Chem._ **32** (2011), 2319–2327,
    DOI: [10.1002/jcc.21787](http://dx.doi.org/10.1002/jcc.21787)

The [HTML version of the
paper](http://onlinelibrary.wiley.com/doi/10.1002/jcc.21787/full) has
so many indentation errors in the code example to make it almost
completely useless for learning MDAnalysis. The
[PDF](http://onlinelibrary.wiley.com/doi/10.1002/jcc.21787/pdf)
version is mostly correct although the final formatted version ended
up with a few line breaks that are illegal in the Python syntax. The
free PubmedCentral version  [PMID
21500218](http://www.ncbi.nlm.nih.gov/pubmed/21500218), however,
contains perfectly formatted Python code.

Corrected code is listed here and can also be found in the
[examples]({{site.github.repo}}/tree/master/package/examples) that are
part of the MDAnalysis package. We are sorry for the inconvenience but
assure you that we tried our best to impart the importance of spaces
in Python code on the typesetter.

We are also collecting any further errata (which correct errors or
omissions in both the PDF and the HTML version). If you find further
problems then please let us know by emailing one of the authors or
filing an issue in the MDAnalysis [Issue Tracker]({{site.github.issues}}).

## Errata ##
### page 4 ###
#### Writing the coordinates of a selection ####

```python
universe.selectAtoms("byres (resname SOL and around 3.5 protein)").write("solvation-shell.pdb")
```

### page 5 ###
#### Block average example ####
The **PDF** version also contains a crucial indentation error.

```python
def blocked(universe, nblocks, analyze): 
    size = universe.trajectory.numframes/nblocks
    blocks = [] 
    for block in xrange(nblocks): 
        a = [] 
        for ts in u.trajectory[block*size:(block+1)*size]: 
            a.append(analyze(universe)) 
        blocks.append(numpy.average(a)) 
    blockaverage = numpy.average(blocks) 
    blockstd = numpy.std(blocks) 
    return nblocks, size, blockaverage, blockstd
```

See [examples/blocks.py]({{site.github.repo}}/blob/master/package/examples/blocks.py)
for a working example. 

#### Radius of gyration analysis function ####
```python
def rgyr(universe): 
    return universe.selectAtoms("protein").radiusOfGyration()
```

#### Time series analysis ####
```python
MDAnalysis.collection.addTimeseries(\ 
    Timeseries.Dihedral(atomselection)) 
data = universe.trajectory.correl(\
    MDAnalysis.collection, skip=10)
```

### page 6 ###
#### radial distribution function ####

The normalisation of the radial distribution function of water is not
correct in the simplified example in the paper; use
[examples/radial\_distribution\_function.py]({{site.github.repo}}/blob/master/package/examples/radial_distribution_function.py)
from the distribution:

```python
# Normalize RDF
radii = 0.5 * (edges[1:] + edges[:-1])
vol = (4./3.)*numpy.pi*(numpy.power(edges[1:], 3)-numpy.power(edges[:-1], 3))
# get water density
import MDAnalysis.core.units
density = MDAnalysis.core.units.convert(1.0, 'water', 'Angstrom^{-3}')
norm = density * (n - 1) / 2 * numframes
# normalize the RDF
rdf /= norm * vol
```

### Table 1 ###

[CHARMM](http://www.charmm.org) has a way to calculate 3D densities
using the [COORdinates ANALysis
... IHIS](http://www.charmm.org/html/documentation/c35b1/corman.html)
facility. The table states incorrectly that CHARMM can not presently
carry out this kind of calculations.
