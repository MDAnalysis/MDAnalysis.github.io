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

## Others

This release contains other performance enhancements and fixes. For a detailed
list see the [release notes]()
