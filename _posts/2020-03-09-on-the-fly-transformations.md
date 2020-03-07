---
layout: post
title: On-the-fly transformations
---

**On-the-fly transformations** have been introduced in version 0.19.0 of MDAnalysis.
This feature is part of @davidercruz 's [Google Summer of Code
2018 project]({{ site.baseurl }}{% post_url 2018-04-26-gsoc-students %}) and brings to 
MDAnalysis a whole new level of functionality, allowing for new and
more efficient workflows when analyzing and visualizing simulation trajectories. 
The documentation for these new functions can be found in the docs for 
[`MDAnalysis.transformations`][otf-docs]

## Why do we need transformations?
When visualizing and analyzing trajectories from molecular dynamics simulations, some
prior modifications are often required.
Examples of the most usual modifications or transformations are removing artifacts
from periodic boundary conditions, which cause some issues with some molecular
viewers (PyMol for example), removing the rotation and translation of a particular
molecule and/or centering it in the unit cell, which helps focus on the its actual
conformational changes by removing their natural movement in solution.
These transformations help us better identify patterns in the behavior of our
biological systems, and, more importantly, show them to the world.

## The advantage of using MDAnalysis for trajectory transformations
Many simulation packages often contain tools to transform and analyze trajectories, such
as [Gromacs][] `gmx trjconv` command. However, most of the times, the user is required to
apply all the intended transformations to the whole trajectory (or the portion of
interest) prior to visualization and analysis. This often requires processing huge files,
sometimes more than once. Moreover, some tools such as `trjconv` do not support frame
indexing for the most popular trajectory formats, requiring iterating over frames that are
not needed for that particular analysis.  Trajectory transformations in MDAnalysis, on the
other end, have one great advantage - they are performed on-the-fly for each frame that is
read. Transformations are added to a universe as a transformation workflow containing one
or more transformations. The API also makes it easy to add new transformations for your
own projects.  Another things that really makes the "on-the-fly" aspect of the MDAnalysis
transformations shine is coupling it to a visualization widget such as [NGL Viewer][].

## Using MDAnalysis transformations
Now it's time to learn how to use the trajectory transformations in MDAnalysis. During the
following steps, we will apply some transformations on a 1 ns trajectory of a simple
19-residue peptide embeded in a 128-DMPC membrane, showing the
[Gromacs][] `gmx trjconv` command and the equivalent MDAnalysis code
and output. To keep things lightweight, frames are were taken every 100 ps, and water
molecules were removed. This can be easily done with MDAnalysis.

### Preparation: Example trajectory
We get our example trajectory from
[`MDAnalysisData.membrane_peptide`](https://www.mdanalysis.org/MDAnalysisData/membrane_peptide.html):

```python
import MDAnalysis as mda
import MDAnalysisData

peptide = MDAnalysisData.datasets.fetch_membrane_peptide()
u = mda.Universe(peptide.topology, peptide.trajectory)
u.transfer_to_memory(step=10)
```

The above commands will download the `peptide.topology` (a Gromacs TPR file name
"memb_pept.tpr") and the `peptide.trajectory` "memb_pept.xtc" in XTC format.

The original trajectory has 1000 frames but for making the visualizations in this post
shorter, we will only keep every 10th frame by using an in-memory representation (see
[Universe.transfer_to_memory()]({{ site.docs.mdanalysis.url
}}/documentation_pages/core/universe.html#MDAnalysis.core.universe.Universe.transfer_to_memory));
when trying these examples yourself you can omit the line
`u.transfer_to_memory(step=10)`. In the following we just write
```python
u = mda.Universe(peptide.topology, peptide.trajectory)
```

### Visualization
We use nglview for visualizing our trajectory in the jupyter notebook. In all cases we add
a unit cell representation and rotate the view with commands such as 
```python
import nglview as nv
import numpy as np

view = nv.show_mdanalysis(u)
view.add_unitcell()
view.control.rotate(
    mda.lib.transformations.quaternion_from_euler(
        -np.pi/2, np.pi/3, np.pi/6, 'rzyz').tolist())
view.control.zoom(-0.3)
view
```
but for simplicity, in the following we only write
```
nv.show_mdanalysis(u)
```

The movies were rendered as animated GIFs with 
```python
from nglview.contrib.movie import MovieMaker
movie = MovieMaker(view, fps=24, output=movie.gif')
movie.make()
```

### Example 1: making everything whole again
When performing MD simulations using periodic boundary conditions, molecules will often
cross the limits of the unit cell. When this happens, some atoms of the molecule will
show up on the the opposing side of the unit cell and some molecular viewers will show
stretched bonds and other visual artifacts depending on the visual representation of the
system. 
This is the case of our system. Without any modifications, when we look at the trajectory
of our system, things become more cluttered and confusing:


```python
import warnings
warnings.filterwarnings('ignore') # nglview is missing some PDB-only attributes and complains 

import MDAnalysis as mda
import nglview as nv

u = mda.Universe(peptide.topology, peptide.trajectory)

nv.show_mdanalysis(u)
```

<img src="{{ site.baseurl }}{{ site.images }}/otf/peptide_raw.gif" title="raw trajectory"
alt="raw trajectory" />

Using `trjconv`, one way to make every molecules whole again would be:

    gmx trjconv -f pept_in_memb.xtc -s pept_in_memb.tpr -pbc mol -o output.xtc
    
In MDAnalysis this can be done with the `unwrap` transformation, which takes an AtomGroup as argument.
This can be done as follows:


```python
from MDAnalysis import transformations

# a custom atom group can be passed as an argument. In this case we will use all the atoms
# in the Universe u
u = mda.Universe(peptide.topology, peptide.trajectory)

# we define the transformation
workflow = [transformations.unwrap(u.atoms)]
```

Now that we have a workflow - in this case it is only a single transformation - we add
it to the `trajectory` object so it can be applied in each frame that we want to read.

```python
u.trajectory.add_transformations(*workflow)
```

If we want to, we can do other things with the trajectory without having to generate a new file
with the transformed trajectory.

This is how our trajectory looks like:

```python
nv.show_mdanalysis(u)
```

<img src="{{ site.baseurl }}{{ site.images }}/otf/peptide_wrapped.gif" title="unwrapped trajectory"
alt="unwrapped trajectory" />


As you can see, the artifacts caused by the atoms crossing the boundaries of the unit cell are now gone.

### Example 2: what if we also want to center the peptide in the unit cell?
In that case, using `trjconv` we would do something like this:
   
    gmx trjconv -f pept_in_memb.xtc -s pept_in_memb.tpr -pbc mol -center -o output.xtc

And we choose `Protein` as the group to be centered.

In MDAnalysis we use the `center_in_box` transformation. As the name says, this transformation will move all the
atoms of the frame, so that a given AtomGroup is centered in the unit cell.
`center_in_box` takes an AtomGroup as a mandatory argument. Optional arguments include `weights`, which is used
to calculate the weighted center of the given AtomGroup (if weights='mass' then the center of mass is calculated),
`center_to` which is used when the user needs to center the AtomGroup in a custom point instead of the center of
the unit cell, and `wrap` which, if `True`, causes all the atoms of the AtomGroup to be moved to the unit cell
before calculating the weighted center.

You can see that the `transformations` workflow below has three steps:
 - make everything molecule whole again with `unwrap` ;
 - center the protein in the unit cell with `center_in_box` - this causes some of the phospholipids to
 fall outside the unit cell ;
 - shift the molecules (`compound='fragments'`) back to the unit cell using `wrap`

This is how it looks:

```python
u = mda.Universe(peptide.topology, peptide.trajectory)
prot = u.select_atoms("protein")
ag = u.atoms
# we will use mass as weights for the center calculation
workflow = (transformations.unwrap(ag),
                   transformations.center_in_box(prot, center='mass'),
                   transformations.wrap(ag, compound='fragments'))
u.trajectory.add_transformations(*workflow)
nv.show_mdanalysis(u)
```
 
<img src="{{ site.baseurl }}{{ site.images }}/otf/peptide_centered.gif" title="centered
and unwrapped trajectory" alt="centered and unwrapped trajectory" />
 
### Example 3: what if we want to do a fitting of the protein?
Fitting is useful when processing trajectories for visualization and analyses - it removes the translations
and rotations of the molecule, allowing us to have a better look at the structural changes that happen in
our simulations.
If we want to do this using `trjconv` we would do have to do this in two steps:

    gmx trjconv -f pept_in_memb.xtc -s pept_in_memb.tpr -pbc mol -center -o midstep.xtc
    gmx trjconv -f midstep.xtc -s pept_in_memb.tpr -fit rot+trans -o output.xtc

And we choose `Protein` as the group to be centered and for the least squares fitting.

In MDAnalysis we just add another transformation to our workflow - `fit_rot_trans`. This transformation takes
the AtomGroup to be fitted as argument, an AtomGroup to be used as reference and, by default, it behaves just
as the option `-fit rot+trans`. If given a `plane` argument, the fitting is performed on a given plane. If
`plane=xy` then the transformation will behave as `-fit rotxy+transxy`, but the `xz` and `yz` planes are also
supported. Just as in `center_in_box`, a `weights` argument can be passed to the function, and it will dictate
how much each atom of the molecule contributes to the least squares fitting. 
Here's what the workflow looks like:


```python
u = mda.Universe(peptide.topology, peptide.trajectory)
prot = u.select_atoms("protein")
# we load another universe to define the reference
# it uses the same input files, but this doesn't have to be always the case
ref_u = u.copy()
reference = ref_u.select_atoms("protein")
ag = u.atoms
workflow = (transformations.unwrap(ag),
                   transformations.center_in_box(prot, center='mass'),
                   transformations.wrap(ag, compound='fragments'),
                   transformations.fit_rot_trans(prot, reference))
u.trajectory.add_transformations(*workflow)
nv.show_mdanalysis(u)
```

<img src="{{ site.baseurl }}{{ site.images }}/otf/peptide_fitted.gif" title="fitted on protein
and unwrapped trajectory" alt="fitted on protein and unwrapped trajectory" />

It looks a bit confusing with the membrane so we can also look at only the protein


```python
view = nv.show_mdanalysis(prot)
view.w.add_line()
view
```

<img src="{{ site.baseurl }}{{ site.images }}/otf/peptideonly_fitted.gif" title="fitted on protein
and unwrapped trajectory (protein only)" alt="fitted on protein and unwrapped trajectory
(protein only)" />


This transformation is good when we want to see how the conformation of the protein evolves with time. 

But, in this case, we also have a membrane. How does the protein behave in the membrane? Doing a least
squares fitting in the `xy` plane can help us have a better look. Here's how it goes:


```python
u = mda.Universe(peptide.topology, peptide.trajectory)
prot = u.select_atoms("protein")
ref_u = u.copy()
reference = ref_u.select_atoms("protein")
ag = u.atoms
workflow = (transformations.unwrap(ag),
                   transformations.center_in_box(prot),
                   transformations.wrap(ag, compound='fragments'),
                   transformations.fit_rot_trans(prot, reference, plane='xy', weights="mass"))
u.trajectory.add_transformations(*workflow)
```

For the visualization we will hide the lipid tails and only indicate the phosphorous
atoms:

```python
protein_P = u.select_atoms("protein or name P")
view = nv.show_mdanalysis(protein_P)
view.add_line()
view
```

<img src="{{ site.baseurl }}{{ site.images }}/otf/peptide_P_fitted_xy.gif" title="fitted on
protein in x-y plane and unwrapped trajectory, protein and lipid phosphorous atoms are
shown" alt="fitted on protein in x-y plane and unwrapped trajectory" />


This transformation keeps the membrane horizontal, while the protein rotation in the z-axis is removed, and
it becomes particularly useful when observing protein insertion.

### Example 4: I want to do my own transformations...
The beauty of MDAnalysis transformations is the ability to easily create custom transformations.
All transformations must have the following structure:
    
```python
def custom_transform(args): # arguments at this point are not mandatory
    #do some things
        
    def wrapped(ts): 
        # This wrapped function must only take a Timestep as argument
        # and perform the actual changes to the timestep
        
        return ts
        
    return wrapped
```
       
Let's create one here:


```python
def up_by_2():
    def wrapped(ts):
        # here's where the magic happens 
        # we create a numpy float32 array to avoid reduce floating
        # point errors
        ts.positions += np.asarray([0,0,20])
        return ts
    return wrapped
```

Now lets add our transformation to a workflow.



```python
import numpy as np
u = mda.Universe(peptide.topology, peptide.trajectory)

# loading another universe to better see the changes made by our transformation
previous = u.copy()
# making the unmodified universe whole accross the trajectory
previous.trajectory.add_transformations(mda.transformations.unwrap(previous.atoms))

ag = u.atoms

workflow = (transformations.unwrap(ag),
                   up_by_2())
u.trajectory.add_transformations(*workflow)
```


All atoms in the `u` Universe are shifted up by 20 angstroms but this does not look
much different from what we have seen before. So let's do something more interesting and
just move a selection such as the peptide.

The transformations can accept arguments. Let's modify `up_by_2` so that only the peptide is translated
in the z coordinate:


```python
def protein_up_by_2(ag):
    def wrapped(ts):
        # here's where the magic happens 
        # we create a numpy float32 array to avoid reduce floating
        # point errors
        ag.positions += np.asarray([0,0,20])
        return ts
    return wrapped
```

We'll add the new transformation to the workflow and see what happens.


```python
u = mda.Universe(peptide.topology, peptide.trajectory)
ag = u.atoms
prot = u.select_atoms("protein")
workflow = (transformations.unwrap(ag),
                   protein_up_by_2(prot),
                   transformations.wrap(ag, compound='fragments'))
u.trajectory.add_transformations(*workflow)
nv.show_mdanalysis(u)
```


<img src="{{ site.baseurl }}{{ site.images }}/otf/peptide_up2.gif" title="peptide
translated by 20 Å upwards" alt="peptide translated by 20 Å upwards" />


The two examples of custom transformations shown here are very simple. But 
more complex things can be done, and we encourage you to try them!

## Final remarks

These transformations are a new feature in MDAnalysis and some transformations such as
wrap/unwrap are still comparatively slow but this post should have given you some ideas
what you will now be able to do. In particular, one can transform a trajectory before any
analysis code sees it so one could implement trajectory smoothing or projections and then
directly analyze the pre-processed trajectory without having to write any intermediate
files or change any of the existing analysis functions.

Transformations behave differently when used with "out of core" trajectories (the normal
approach in MDAnalysis, where each trajectory frame is read from disk into memory when
needed) and "in core" trajectories (generated with `Universe.transfer_to_memory()`, also
known as the "MemoryReader"). For on-disk trajectories, the transformations are performed
whenever a frame is read from disk. For in-memory trajectories, the transformations are
applied once to and the modified trajectory is stored in memory. Therefore, in-memory
trajectories with transformations can appear to take a long time to load because all
calculations are done immediately.

This has been a quick demonstration of the power of the new on-the-fly transformations of
MDAnalysis. There are more transformations available for you to explore and a whole lot
more for you to create for your own molecular system. More [information on trajectory
transformations][otf-docs] can be found in the online docs of MDAnalysis. 




— @davidercruz

[Gromacs]: http://www.gromacs.org
[NGL Viewer]: http://nglviewer.org/ngl/api/index.html
[otf-docs]: {{ site.docs.mdanalysis.url }}/documentation_pages/trajectory_transformations.html
