---
layout: post
title: Release 0.16.0
---

We have just released MDAnalysis version 0.16.0. This release contains new
features as well as bug fixes. Highlights are listed below but for more details
see the [release notes](https://github.com/MDanalysis/mdanalysis/wiki/...).

13 People contributed to this release.

# Upgrade

You can upgrade with `pip install --upgrade MDAnalysis`

# Noticable Changes

## Attach arbitraty time series to your trajectories

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

This feature is still in it's beginning and will be expanded in future releases.
So far, only the XVG format used by [gromacs] and [grace] are supported. Open an issue
if you need support for other time series formats.

## Do a dimension reduction with PCA and Dmaps

@jdetle has implemented two new dimension reduction algorithms,
[Principal Component Analysis](pca) and [Diffusion Maps](dmaps-paper). Both can
be found in the analysis submodule. As an example lets look at the first two PCA
dimensions of [ADK](adk-sim).

```python
import MDAnalyis as mda
from MDAnalysis.analysis.pca import PCA

u = mda.Universe(ADK_SIMULATION)
backbone = u.select_atoms('backbone)

pca = PCA(backbone).run()
reduced_data = pca.transform(backbone, n_components=2)

f, ax = plt.subplots()
ax.scatter(reduced_data[:, 0], reduced_data[:, 1], edgecolor='none')
ax.set(xlabel='PC_1', ylabel='PC_2', title='PCA of ADK')
```

**TODO**: add picture

## Speed improvements in RMSD

Thanks for work from @rbrtdlgd our RMSD calculations are about 40% faster now.
If you are using the low-level qcprot algorithm your self intead of our provided
wrappers you have to change your code since the API has changed. For more see
the [CHANGELOG].

# Minor Enhancements

- No more deprecation warning spam when MDAnalyis is imported
- analysis.align has a new AlignTraj class following the analysis class style
- all new analysis classes now print additional information with the `quiet=False`.
- new analysis ported to new style: AlignTraj, RMSD,


# Other Changes

A list of all changes can be found in the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/develop/package/CHANGELOG).

[dmaps-paper]: http://look.me.up.a.clementi.md.paper
[pca]: http://wikipedia?
[adk-sim]: link to download instructions
[aux-doc]: link to auxiliary doc
[fiona-blog]: http://fiona-naughton.github.io/blog/
[gromacs]: http://www.gromacs.org
[grace]: http://plasma-gate.weizmann.ac.il/Grace/
