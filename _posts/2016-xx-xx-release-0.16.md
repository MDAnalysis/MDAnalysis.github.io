---
layout: post
title: Release 0.16.0
---

We have just released MDAnalysis version 0.16.0. This release contains new
features as well as bug fixes. Highlights are listed below but for more details
see the [release notes](https://github.com/MDanalysis/mdanalysis/wiki/...).

This release includes the work of our GSoC students Fiona Naughton
(@fiona-naughton) and John Detlefs (@jdetle). We have a big release this time
with tones of new features, performance improvements and bug fixes. In total
**XX** People contributed to this release.

# Upgrade

You can upgrade with `pip install --upgrade MDAnalysis`

# Noticable Changes

## Attach arbitrary time series to your trajectories

Our GSoC student @fiona-naughton has implemented an auxillary reader to add
arbitrary time series to a universe. The time series are kept in sync with the
trajectory so it is possible to iterate through the trajectory and access the
auxiliary data corresponding to the current time step.

```python
import MDAnalysis as mda

# Create your universe as usual
universe = mda.Universe(SIMULATION)
# Attach an auxiliary time serie with the name `pull_force`
universe.trajectory.add_auxiliaty('pull_force', 'md_f.xvg')
# Itarete through your trajectory, the time serie is kept in sync
for time_step in universe.trajectory:
    print(time_step.aux.pull_force)
```

@fiona-naugthon worked at offering several convenient way to iterate through your
data. Read the [documentation](aux-doc) or [Fiona's blog posts](fiona-blog) to learn more about the feature.

This feature is still in it's beginning and will be expanded in future releases. You can
follow the conversation on the [initial issue](issue785) or on the [pull request](pr868).
So far, only the XVG format used by [gromacs] and [grace] are supported. Open an issue
if you need support for other time series formats.

## Do a dimension reduction with PCA and Dmaps

@jdetle has implemented two new dimension reduction algorithms,
[Principal Component Analysis](pca) and [Diffusion Maps](dmaps-paper). Both can
be found in the analysis submodule. As an example lets look at the first two PCA
dimensions of ADK from our test files.

```python
import matplotlib.pyplot as plt
import MDAnalyis as mda
from MDAnalysis.analysis.pca import PCA
from MDAnalyisTests.datafiles import PSF, DCD

plt.style.use('ggplot')

u = mda.Universe(PSF, DCD)
ca = u.select_atoms('protein and name CA')

pca = PCA(u, select='protein and name CA', quiet=False).run()
reduced_data = pca.transform(ca, n_components=2)

f, ax = plt.subplots()
ax.plot(d[:, 0], d[:, 1], 'o')
ax.set(xlabel=r'PC$_1$ [$\AA$]', ylabel=r'PC$_2$ [$\AA$]', title='PCA of ADK')
```

![PCA projection]({{site.images}}pca-release-0.16.png)

## Convenience functions to create a new analysis

A while back we introduced a new frame work for analysis to unify the API for
the different analysis methods we offer. With this release we also add a new
class `AnalysisFromFunction` to make it easier to calculate observables from a
simulation. Now code like this with a handwritten loop.

```python
result = []
for ts in u.trajectory:
    result.append(u.atoms.center_of_geometry())
results = np.asarray(results)
```

Can now be converted into this.

```python
from MDAnalyis.analysis.base import AnalysisFromFunction
cog = AnalysisFromFunction(lambda ag : ag.center_of_geometry(), u.atoms).run()
cog.results
```

This class also takes arguments to adjust the iteration (`start`,`stop`,`step`)
and you can add verbosity with `quiet=False`. You will also profit from any
performance improvements in the analysis class in the future without changing
your code. If you have a specific observable that you want to calculate several
times you can also create a new analysis class with `analysis_class` like this.

```python
from MDAnalyis.analysis.base import analysis_class

def cog(ag):
    return ag.center_of_geometry()

COG = analysis_class(cog)

cog_results = COG(u.atoms, step=2, quiet=False).run()
```

## Speed improvements in RMSD

Thanks for work from @rbrtdlgd our RMSD calculations are about 40% faster now.
If you are using the low-level qcprot algorithm your self intead of our provided
wrappers you have to change your code since the API has changed. For more see
the [CHANGELOG].

## MemoryReader: Reading trajectories from memory

MDAnalysis typically reads trajectories from files on-demand, so that it can efficiently deal with large trajectories - even those that do not fit in memory. However, in some cases, both for convenience and for efficiency, it can be an advantage to work with trajectories directly in memory. In this release, we have introduced a MemoryReader, which makes this possible.

The MemoryReader works with numpy arrays, using the same format as that used by for instance `DCDReader.timeseries()`. You can create a Universe directly from such an array:

```python
import numpy as np
from MDAnalysis import Universe
from MDAnalysisTests.datafiles import DCD, PSF
from MDAnalysis.coordinates.memory import MemoryReader

# Create a Universe using a DCD reader
universe = Universe(PSF, DCD)

# Create a numpy array with random coordinates (100 frames) for the same topology
coordinates = np.random.uniform(size=(universe.atoms.n_atoms, 100, 3)).cumsum(0)

# Create a new Universe directly from these coordinates
universe2 = Universe(PSF, coordinates, format=MemoryReader)
```

The MemoryReader will work just as any other reader. In particular, you can iterate over it as usual, or use the `.timeseries()` method to retrieve a reference to the raw array:

```python
coordinates_fac = universe2.trajectory.timeseries(format='fac')
```

Certain operations can be speeded up by moving a trajectory to memory, and we have therefore
added functionality to directly transfer any existing trajectory to a MemoryReader using `Universe.transfer_to_memory`:

```python
universe = Universe(PSF, DCD)
universe.transfer_to_memory()     # Switches to a MemoryReader representation
```

You can also do this directly upon construction of a Universe, by using the `in_memory` flag:

```python
universe = Universe(PSF, DCD, in_memory=True)
```

Likewise, the `AlignTraj` class in the analysis/align.py module also has an `in_memory` flag, allowing it to do in-place alignments in memory.


## Incorporation of the ENCORE ensemble similarity library

The **ENCORE** ensemble similarity library has been integrated with MDAnalysis as [MDAnalysis.analysis.encore](http://docs.mdanalysis.org/documentation_pages/analysis/encore.html). It implements a variety of techniques for calculating similarities between structural ensembles (trajectories), as described in this publication:

Tiberti M, Papaleo E, Bengtsen T, Boomsma W, Lindorff-Larsen K (2015), ENCORE: Software for Quantitative Ensemble Comparison. PLoS Comput Biol 11(10): e1004415. doi:[10.1371/journal.pcbi.1004415](http://doi.org/10.1371/journal.pcbi.1004415).

Using the similarity measures is simply a matter of loading the trajectories or experimental ensembles that one would like to compare as MDAnalysis.Universe objects:

```python
from MDAnalysis import Universe
import MDAnalysis.analysis.encore as encore
from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
u1 = Universe(PSF, DCD)
u2 = Universe(PSF, DCD2)
```

and running the similarity measures on them, choosing among 1) the Harmonic Ensemble Similarity measure:

```python
hes_similarities, details = encore.hes([u1, u2])
print hes_similarities
```
```
[[        0.         38279683.9587939]
 [ 38279683.9587939         0.       ]]
```

2) the Clustering Ensemble Similarity measure:

```python
ces_similarities, details = encore.ces([u1, u2])
print ces_similarities
```
```
[[ 0.          0.68070702]
 [ 0.68070702  0.        ]]
```

or 3) the Dimensionality Reduction Ensemble Similarity measure:

```python
dres_similarities, details = encore.dres([u1, u2])
print dres_similarities
```
```
[[ 0.          0.65434461]
 [ 0.65434461  0.        ]]
```
Similarities are written in a square symmetric matrix having the same dimensions and ordering as the input list, with each element being the similarity value for a pair of the input ensembles. 

The encore library includes a general interface to various clustering and dimensionality reduction algorithms (through the [scikit-learn](http://scikit-learn.org/) package), which makes it easy to switch between clustering and dimensionality reduction algorithms when using the `ces` and `dres` functions. The clustering and dimensionality reduction functionality is also directly available through the `cluster` and `reduce_dimensionality` functions. For instance, to cluster the conformations from the two universes defined above, we can write:
```python
cluster_collection = encore.cluster([u1,u2])
print cluster_collection
```
```
0 (size:5,centroid:1): array([ 0,  1,  2,  3, 98])
1 (size:5,centroid:6): array([4, 5, 6, 7, 8])
2 (size:7,centroid:12): array([ 9, 10, 11, 12, 13, 14, 15])
…
```
In addition to standard cluster membership information, the `cluster_collection` output keep track of the origin of each conformation, so you check how the different trajectories are represented in each cluster. For further details, see the documentation of the individual functions within Encore.


# Minor Enhancements

- No more deprecation warning spam when MDAnalyis is imported
- analysis.align has a new AlignTraj class following the analysis class style
- all new analysis classes now print additional information with the `quiet=False`.
- RMSD has been ported to the new analysis class style

# Other Changes

A list of all changes can be found in the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/develop/package/CHANGELOG).

[dmaps-paper]: dx.doi.org/10.1073/pnas.0500334102
[pca]: https://en.wikipedia.org/wiki/Principal_component_analysis
[aux-doc]: http://www.mdanalysis.org/MDAnalysis/documentation_pages/auxiliary/init.html
[fiona-blog]: http://fiona-naughton.github.io/blog/
[isue785]: https://github.com/MDAnalysis/mdanalysis/issues/785
[pr868]: https://github.com/MDAnalysis/mdanalysis/pull/868
[gromacs]: http://www.gromacs.org
[grace]: http://plasma-gate.weizmann.ac.il/Grace/
