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

## Access Force/pressure/... data in a timestep

Our GSoC student @fiona-naughton has implemented a auxillary reader to add
arbitrary data to a universe.

**TODO** Code examples of usage

This feature is still in it's beginning and will be expanded in future releases

## Do a dimension reduction with PCA and Dmaps

@jdetle has implemented two new dimension reduction algorithms,
[Principal Component Analysis](pca) and [Diffusion Maps](dmaps-paper). Both can
be found in the analysis submodule. As an example lets look at the first two PCA
dimensions of [ADK](adk-sim).

    import MDAnalyis as mda
    from MDAnalysis.analysis.pca import PCA

    u = mda.Universe(ADK_SIMULATION)
    backbone = u.select_atoms('backbone)

    pca = PCA(backbone).run()
    reduced_data = pca.transform(backbone, n_components=2)

    f, ax = plt.subplots()
    ax.scatter(reduced_data[:, 0], reduced_data[:, 1], edgecolor='none')
    ax.set(xlabel='PC_1', ylabel='PC_2', title='PCA of ADK')

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
