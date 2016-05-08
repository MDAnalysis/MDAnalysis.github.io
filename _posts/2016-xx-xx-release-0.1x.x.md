---
layout: post
title: Release 0.1x.x
---

We have just released MDAnalysis version 0.1x.x

# Upgrade

You can upgrade with `pip install --upgrade MDAnalysis`

# Noticable Changes

## Rewrite of TRR and XTC file handling

Our implementation of the Gromacs
[xdrlib](http://www.gromacs.org/Developer_Zone/Programming_Guide/XTC_Library)
has been completely rewritten in cython this release. This changes brings us one
step closer towards supporting Python 3. The only user facing API change is that
we don't save persistent frame offsets with the pickle module anymore but
with numpy's 'npz' format. This improves reopening of xtc/trr files.

## Going strong towards Python 3

We are planning to support Python 3 as well in the future. This release has
started our process of adapting MDAnalysis to be compatible with Python 2.7 and
Python 3.

## Deprecations in anticipation of new topology system

The next release (0.16.0) will bring a very big change to the internal workings of MDAnalysis.
The topology system (how atom, residue, segment attributes are stored and manipulated) has been
completely redesigned, giving many advantages and resolving a number of longstanding issues.
You can read more about this at the [original issue](https://github.com/MDAnalysis/mdanalysis/issues/363),
or see [a short summary of the new system](https://github.com/MDAnalysis/mdanalysis/wiki/Issue363-Changes)
on the wiki.

In preparation for this change, we have introduced [deprecation
warnings](https://github.com/MDAnalysis/mdanalysis/issues/599) for all
components of the existing topology system suggesting the corresponding usage
under the new system.  Please adjust existing code as you encounter these
warnings. 

## Others

This release contains other performance enhancements and fixes. For a detailed
list see the [release notes]()
