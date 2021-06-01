
## Writing Custom Readers and Parsers for MDAnalysis


molecular simulation full of custom formats, or variants on existing formats

these can be tricky to read as columns get rearranged

similarly, the number of fields used can vary greatly between different formats

in this blog, we show how MDAnalysis can be used to read any format,
including adding custom attributes


```python
myfile = """\
7
My custom file thingy
C red      1.0 1.0 1.0
H blue     1.0 1.5 1.5
H blue     1.0 0.5 0.5
C red      2.0 2.0 2.0
H yellow   2.0 2.5 2.5
H yellow   2.0 1.5 1.5
O octarine 3.0 3.0 3.0
"""
with open('myfile.xyzc', 'w') as out:
    out.write(myfile)
```

In the above example, we have something that looks mostly like an XYZ file, with an extra column which indicates the colour of particles.  This will trip up the existing XYZ parser as the position of the coordinate information has moved.


```python
import MDAnalysis as mda

u = mda.Universe('myfile.xyzc', format='XYZ')
```


    ---------------------------------------------------------------------------

    EOFError                                  Traceback (most recent call last)

    <ipython-input-2-19807941e521> in <module>()
          1 import MDAnalysis as mda
          2 
    ----> 3 u = mda.Universe('myfile.xyzc', format='XYZ')
    

    /home/richard/miniconda2/envs/mydev/lib/python2.7/site-packages/MDAnalysis/core/universe.pyc in __init__(self, *args, **kwargs)
        268             else:
        269                 coordinatefile = args[1:]
    --> 270             self.load_new(coordinatefile, **kwargs)
        271 
        272         # Check for guess_bonds


    /home/richard/miniconda2/envs/mydev/lib/python2.7/site-packages/MDAnalysis/core/universe.pyc in load_new(self, filename, format, in_memory, **kwargs)
        416         kwargs['n_atoms'] = self.atoms.n_atoms
        417 
    --> 418         self.trajectory = reader(filename, **kwargs)
        419         if self.trajectory.n_atoms != len(self.atoms):
        420             raise ValueError("The topology and {form} trajectory files don't"


    /home/richard/miniconda2/envs/mydev/lib/python2.7/site-packages/MDAnalysis/coordinates/XYZ.pyc in __init__(self, filename, **kwargs)
        309         # (Also cannot just use seek() or reset() because that would break
        310         # with urllib2.urlopen() streams)
    --> 311         self._read_next_timestep()
        312 
        313     @property


    /home/richard/miniconda2/envs/mydev/lib/python2.7/site-packages/MDAnalysis/coordinates/XYZ.pyc in _read_next_timestep(self, ts)
        369             return ts
        370         except (ValueError, IndexError) as err:
    --> 371             raise EOFError(err)
        372 
        373     def _reopen(self):


    EOFError: could not convert string to float: red


Here, as expected, the XYZReader chokes on the word 'red' where it expects a floating point number.

We can rewrite the XYZReader to anticipate the extra column:


```python
import numpy as np
from MDAnalysis.coordinates.XYZ import XYZReader
from MDAnalysis.topology.XYZParser import XYZParser

class XYZCReader(XYZReader):
    format = 'XYZC'
    
    def _read_next_timestep(self, ts=None):
        # check that the timestep object exists
        if ts is None:
            ts = self.ts

        f = self.xyzfile

        try:
            # we assume that there are only two header lines per frame
            f.readline()
            f.readline()
            for i in range(self.n_atoms):
                self.ts._pos[i] = np.float32(f.readline().split()[2:5])
            ts.frame += 1
            return ts
        except (ValueError, IndexError) as err:
            raise EOFError(err)
            
class XYZCParser(XYZParser):
    format = 'XYZC'
```


```python
u = mda.Universe('myfile.xyzc')
```


```python
print u.atoms.names
print u.atoms.positions
```

    ['C' 'H' 'H' 'C' 'H' 'H' 'O']
    [[ 1.   1.   1. ]
     [ 1.   1.5  1.5]
     [ 1.   0.5  0.5]
     [ 2.   2.   2. ]
     [ 2.   2.5  2.5]
     [ 2.   1.5  1.5]
     [ 3.   3.   3. ]]


Now it seems to work.

But what about our colour field?


```python
# to import all the stuff that the original XYZParser used
from MDAnalysis.topology.XYZParser import *
from MDAnalysis.core.topologyattrs import AtomAttr
from MDAnalysis.core.groups import AtomGroup
from collections import defaultdict


class Colours(AtomAttr):
    attrname = 'colours'
    singular = 'colour'
    transplants = defaultdict(list)

    @property
    def complimentary(self):
        COMP = {
            'red': 'green', 'green': 'red',
            'yellow': 'purple', 'purple': 'yellow',
            'blue': 'orange', 'orange': 'blue',
            'octarine': '???', '???': 'octarine'
        }
        return np.array([COMP[c] for c in self.colours], dtype=object)

    transplants[AtomGroup].append(('complimentary', complimentary))

    
class XYZCParser(XYZParser):
    format='XYZC'
    
    def parse(self):
        with openany(self.filename, 'r') as inf:
            natoms = int(inf.readline().strip())
            inf.readline()

            names = np.zeros(natoms, dtype=object)
            colours = np.zeros(natoms, dtype=object)

            for i, line in enumerate(inf):
                if i == natoms:
                    break
                names[i] = line.split()[0]
                colours[i] = line.split()[1]

        # Guessing time
        atomtypes = guessers.guess_types(names)
        masses = guessers.guess_masses(atomtypes)

        attrs = [Atomnames(names),
                 Atomids(np.arange(natoms) + 1),
                 Atomtypes(atomtypes, guessed=True),
                 Masses(masses, guessed=True),
                 Resids(np.array([1])),
                 Resnums(np.array([1])),
                 Segids(np.array(['SYSTEM'], dtype=object)),
                 Colours(colours),
                 ]

        top = Topology(natoms, 1, 1,
                       attrs=attrs)

        return top
```


```python
u = mda.Universe('myfile.xyzc')
```


```python
print u.atoms.colours
print u.atoms.complimentary
```

    ['red' 'blue' 'blue' 'red' 'yellow' 'yellow' 'octarine']
    ['green' 'orange' 'orange' 'green' 'purple' 'purple' '???']

