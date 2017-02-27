---
layout: post
title: A shiny, new topology system 
---

With MDAnalysis 0.16.0 on the horizon, we wanted to showcase a major development that most users will probably not notice if we've done our job well.
In fall 2015, @richardjgowers and I set to work on redesigning the topology system from scratch.
This system determines how atom, residue, and segment information is internally represented and exposed to everything in the API (``Universe``, ``AtomGroup``, etc.), and the old scheme had issues with data duplication, maintaining consistency between atom and residue attributes, and performance for large systems.
We hoped to resolve all of these issues with our new design.

The starting point of this work was (the now infamous) [issue 363](https://github.com/MDAnalysis/mdanalysis/issues/363), which floated the idea of holding all atom, residue, and segment attributes in arrays instead of lists of ``Atom``, ``Residue``, and ``Segment`` objects.
This approach turned the way topology data such as atom names, resids, masses, etc. are stored in a ``Universe`` on its head, going from an array of structs (list of ``Atom`` objects with individual attributes) to a struct of arrays (an array for each attribute, one entry per ``Atom``).

Now, over a year later, the finishing touches on this work are being prepared for release.
This post is meant to serve as a brief view to what has changed internally, what has changed externally, and what benefits this gives us looking forward to the future.

## Internal changes that shouldn't affect external behavior

In the new system, each atom is a member of exactly one residue, and each residue is a member of exactly one segment.
The new `Topology` object keeps an array giving the residue membership of each atom, and likewise an array giving segment membership of each residue.
Getting the resname of the residue of a group of atoms, then, is achieved by taking the indices of these atoms to fancy-index the `Atoms->Residues` array, and then using the result of this to fancy-index the `Resnames` array.
For example, if the  `Topology` has 5 atoms and 3 residues, with membership (`Atoms->Residues`) and `Resnames` arrays as below:

```
       Atoms->Residues           Resnames
 index ---------------     index --------
     0 0                       0 GLU
     1 2                       1 LYS
     2 1                       2 ALA
     3 1
     4 2
```

calling `AtomGroup.resnames` for an `AtomGroup` with atoms [2, 0, 1, 2] will yield (pseudocode):

```
"Atoms->Residues"[[2, 0, 1, 2]] --> [1, 0, 2, 1]
"Resnames"[[1, 0, 2, 1]]        --> ['LYS', 'GLU', 'ALA', 'LYS']
```

This scheme only works if each atom is a member of one and only one residue, and likewise if residues are members of one and only one segment.
Furthermore, `AtomGroup`s, `ResidueGroup`s, and `SegmentGroup`s are very thin, storing only the indices of their members as a `numpy` array.
This gives a number of advantages:

1. **Performance**. We get up to an 8x speedup over the old scheme when accessing attributes. Setting attributes can give up to a 40x speedup.
2. **Memory**. We don't store, for example, a resname for each atom, but instead store attributes at the level they make sense for.
3. **Consistency**. Since attributes are stored in one place, we avoid cases where the topology is in an inconsistent state, e.g. two atoms in the same residue give a different resname.
4. **No staleness**. Because e.g. `ResidueGroup`s are only an array of indices, not a list of `Residue` objects generated upon creation of the group, changes of resiude-level properties by another `ResidueGroup` are always reflected consistently by every other one. Data is not duplicated anywhere in this scheme, and is all contained in the `Topology` object.

For further performance comparisons, check out this [notebook](http://nbviewer.jupyter.org/gist/dotsdl/0e0fbd409e3e102d0458).

## External changes that may affect how you use MDAnalysis

Previously, every object except ``Atom`` subclassed from ``AtomGroup``.
This meant that calling `.positions` of would give you the positions of the ``Atom``s contained within that group.

Previous class structure:
```
Atom

AtomGroup  -> Residue
           -> ResidueGroup -> Segment
                           -> SegmentGroup
```
New class structure:
```
Group    -> AtomGroup
         -> ResidueGroup
         -> SegmentGroup

Atom
Residue
Segment
```

Now each object only contains information pertaining to that particular object.
A ``Residue`` object only yields information about the residue; to get to the atoms, use ``Residue.atoms``.

### Why this was changed

Previously everything inheriting from ``AtomGroup`` made it unclear at what level of topology a given method or attribute was working on.
For example, does ``ResidueGroup.charges`` give the charge of the residues or the atoms?
Also, it was unclear what size a given output would be (see [issue 411](https://github.com/MDAnalysis/mdanalysis/issues/411)).

### How to work around this

To access atom-level information from anything that isn't an ``AtomGroup``, use the `.atoms` level accessor.
For example, changing all `.positions` calls on anything that isn't an `AtomGroup` to `.atoms.positions`.


## Going forward: what does this mean for MDAnalysis as a project?

A major benefit of the new topology system is that information about the topology of a ``Universe`` is now completely encapsulated in the ``Topology`` object.
This not only makes development and maintenance easier, but also opens the door to some exciting new possibilities as simulation systems grow larger.
A single ``Topology`` object can now be cleanly shared by multiple ``Universe`` instances, each with their own trajectory reader(s), making common operations such as fitting a trajectory to a reference structure or doing parallel analysis of many trajectories more feasible for large systems.
The ``Topology`` object can also be serialized more easily, making parallelization on workers without shared memory (using libraries such as [``distributed``](http://distributed.readthedocs.io/en/latest/)) within the realm of possibility out-of-the-box.

Making these things work is an ongoing effort, but the MDAnalysis [coredevs](https://github.com/orgs/MDAnalysis/teams/coredevs) are working to take advantage of all these possibilities.
We look forward to the benefits this brings not only to the project, but also to all our users going forward.
We hope you like what we've done here.

-- @dotsdl
