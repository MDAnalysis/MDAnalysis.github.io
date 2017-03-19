---
layout: page
title: Basic example
---

A typical usage pattern is to iterate through a trajectory and analyze
coordinates for every frame. In the following example the end-to-end
distance of a protein and the radius of gyration of the backbone atoms
are calculated:

{% highlight python %}
import MDAnalysis
from MDAnalysis.tests.datafiles import PSF, DCD   # test trajectory
import numpy.linalg

u = MDAnalysis.Universe(PSF,DCD)  # always start with a Universe
nterm = u.s4AKE.atoms.N[0]   # can access via segid (s4AKE) and atom name
cterm = u.s4AKE.atoms.C[-1]  # ... takes the last atom named 'C'
bb = u.select_atoms('protein and backbone')  # a selection (AtomGroup)

for ts in u.trajectory:     # iterate through all frames
    r = cterm.position - nterm.position # end-to-end vector from atom positions
    d = numpy.linalg.norm(r)  # end-to-end distance
    rgyr = bb.radius_of_gyration()  # method of AtomGroup
    print("frame = {0}: d = {1} A, Rgyr = {2} A".format(
          ts.frame, d, rgyr))
{% endhighlight %}

To find out what else you can do, head over to [learning
MDAnalysis]({{site.baseurl}}pages/learning_MDAnalysis) to have a look
at the tutorials and the docs.
