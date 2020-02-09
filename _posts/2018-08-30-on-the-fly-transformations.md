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
[`MDAnalysis.transformations`](https://www.mdanalysis.org/docs/documentation_pages/trajectory_transformations.html).

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
Many simulation packages often contain tools to transform and analyze trajectories,
such as GROMACS' `trjconv`. However, most of the times, the user is required to apply
all the intended transformations to the whole trajectory (or the portion of interest)
prior to visualization and analysis. This often requires processing huge files, sometimes
more than once. Moreover, some tools such as `trjconv` do not support frame indexing
for the most popular trajectory formats, requiring iterating over frames that are not
needed for that particular analysis.
Trajectory transformations in MDAnalysis, on the other end, have one great advantage - they
are performed on-the-fly for each frame that is read. Transformations are added to a universe
as a transformation workflow containing one or more transformations. The API also makes it easy
to add new transformations for your own projects.
Another things that really makes the "on-the-fly" aspect of the MDAnalysis transformations
shine is coupling it to a visualization widget such as
[NGL Viewer](http://nglviewer.org/ngl/api/index.html).

## Using MDAnalysis transformations
Now it's time to learn how to use the trajectory transformations in MDAnalysis. During the
following steps, we will apply some transformations on a 1 ns trajectory of a simple
19-residue peptide embeded in a 128-DMPC membrane, showing the GROMACS `trjconv` command
and the equivalent MDAnalysis code and output. To keep things lightweight, frames are
were taken every 100 ps, and water molecules were removed. This can be easily done with
MDAnalysis.

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
warnings.filterwarnings('ignore') # some attributes are missing 
import MDAnalysis as mda
from MDAnalysis import transformations
import nglview as nv
u = mda.Universe('pept_in_memb.tpr', 'pept_in_memb.xtc')
nv.show_mdanalysis(u)
```


Using `trjconv`, one way to may every molecules whole again would be:

    gmx trjconv -f pept_in_memb.xtc -s pept_in_memb.tpr -pbc mol -o output.xtc
    
In MDAnalysis this can be done with the `unwrap` transformation, which takes an AtomGroup as argument.
This can be done as follows:


```python
# a custom atom group can be passed as an argument. In this case we will use all the atoms
# in the Universe u
u = mda.Universe('pept_in_memb.tpr', 'pept_in_memb.xtc')
ag = u.atoms
# we define the transformation
transformations = transformations.unwrap(ag)
```

Now that we have a workflow - in this case it is only a single transformation - we add
it to the `trajectory` object so it can be applied in each frame that we want to read.

```python
u.trajectory.add_transformations(transformations)
```

If we want to, we can do other things with the trajectory without having to generate a new file
with the transformed trajectory.

This is how our trajectory looks like:

```python
nv.show_mdanalysis(u)
```

as you can see, the artifacts caused by the atoms crossing the boundaries of the unit cell are now gone.

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
u = mda.Universe('pept_in_memb.tpr', 'pept_in_memb.xtc')
prot = u.select_atoms("protein")
ag = u.atoms
# we will use mass as weights for the center calculation
transformations = (transformations.unwrap(ag),
                   transformations.center_in_box(prot, center='mass'),
                   transformations.wrap(ag, compound='fragments'))
u.trajectory.add_transformations(*transformations)
nv.show_mdanalysis(u)
```
 
 
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
u = mda.Universe('pept_in_memb.tpr', 'pept_in_memb.xtc')
prot = u.select_atoms("protein")
# we load another universe to define the reference
# it uses the same input files, but this doesn't have to be always the case
ref_u = u.copy()
reference = ref_u.select_atoms("protein")
ag = u.atoms
transformations = (transformations.unwrap(ag),
                   transformations.center_in_box(prot, center='mass'),
                   transformations.wrap(ag, compound='fragments'),
                   transformations.fit_rot_trans(prot, reference))
u.trajectory.add_transformations(*transformations)
nv.show_mdanalysis(u)

```


It looks a bit confusing with the membrane...


```python
t = nv.MDAnalysisTrajectory(prot)
w = nv.NGLWidget(t)
w.add_line()
w
```


This transformation is good when we want to see how the conformation of the protein evolves with time. 

But, in this case, we also have a membrane. How does the protein behave in the membrane? Doing a least
squares fitting in the `xy` plane can help us have a better look. Here's how it goes:


```python
u = mda.Universe('pept_in_memb.tpr', 'pept_in_memb.xtc')
prot = u.select_atoms("protein")
ref_u = u.copy()
reference = ref_u.select_atoms("protein")
ag = u.atoms
transformations = (transformations.unwrap(ag),
                   transformations.center_in_box(prot),
                   transformations.wrap(ag, compound='fragments'),
                   transformations.fit_rot_trans(prot, reference, plane='xy', weights="mass"))
u.trajectory.add_transformations(*transformations)
# let's hide the lipid tails to have a better view
view_selection = u.select_atoms("protein or name P")
t = nv.MDAnalysisTrajectory(view_selection)
w = nv.NGLWidget(t)
w.add_line()
w
```


This transformation keeps the membrane horizontal, while the protein rotation in the z-axis is removed, and
it becomes particularly useful when observing protein insertion.

### Example 4: I want to do my own transformations...
The beauty of MDAnalysis transformations is the ability to easily create custom transformations.
All transformations must have the following structure:
    
```Python
def custom_transform(args): # arguments at this point are not mandatory
    #do some things
        
    def wrapped(ts): 
        # this wrapped function must only take a Timestep as argument
        # this is where the actual changes to the timestep must be done
            
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
u = mda.Universe('pept_in_memb.tpr', 'pept_in_memb.xtc')

# loading another universe to better see the changes made by our transformation
previous = mda.Universe('pept_in_memb.tpr', 'pept_in_memb.xtc')
# making the unmodified universe whole accross the trajectory
previous.trajectory.add_transformations(mda.transformations.unwrap(previous.atoms))

ag = u.atoms

transformations = (transformations.unwrap(ag),
                   up_by_2())
u.trajectory.add_transformations(*transformations)
view_selection = previous.select_atoms("protein or name P")
t = nv.MDAnalysisTrajectory(view_selection)
w = nv.NGLWidget(t)
w.add_line()
w.add_trajectory(u)
w
```


As you can see, the atoms in the `u` Universe have been shifted up by 20 angstroms.

The transformations can accept arguments. Let's modify `up_by_2` so that only the peptide is translated
in the z coordinate:


```python
def protein_up_by_2(agroup):
    
    def wrapped(ts):
        # here's where the magic happens 
        # we create a numpy float32 array to avoid reduce floating
        # point errors
        agroup.positions += np.asarray([0,0,20])
        
        return ts
    
    return wrapped
```

We'll add the new transformation to the workflow and see what happens.


```python
u = mda.Universe('pept_in_memb.tpr', 'pept_in_memb.xtc')
ag = u.atoms
prot = u.select_atoms("protein")
transformations = (transformations.unwrap(ag),
                   protein_up_by_2(prot),
                   transformations.wrap(ag, compound='fragments'))
u.trajectory.add_transformations(*transformations)
nv.show_mdanalysis(u)
```


The two examples of custom transformations shown here are very simple. But 
more complex things can be done, and we encourage you to try them!

This has been a quick demonstration on the power of the new on-the-fly transformations
of MDAnalysis. There are more transformations available for you to explore and a whole
lot more for you to create for your own molecular system. More information on trajectory
transformations can be found in the [online docs of MDAnalysis](https://www.mdanalysis.org/mdanalysis).
