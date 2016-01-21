---
layout: post
title: Release 0.14.0
---

We have just released MDAnalysis version 0.14.0

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

## Impicit OR in selections

Many long selection strings used in select_atoms have been simplified through
allowing implicit OR in the arguments.  For example to select all atoms
with one of a few names previous required lots of ORs

```python
# OLD
u.select_atoms('name Ca or name N or name Ch')

# NEW
u.select_atoms('name Ca N Ch')
```

The new syntax allows multiple arguments after the keyword `name`.
The selection will keep eating arguments till it hits a keyword.
The use of wildcards is still possible too, making the selection
of all atoms with a type beginning with 'C' or 'N' as simple as:

```python
u.select_atoms('type C* N*')
```

Similarly, for selecting ranges of resids

```python
# OLD
u.select_atoms('resid 1:10 or resid 40:50 or resid 56 or resid 67')

# NEW
u.select_atoms('resid 1:10 40:50 56 67')
```

This new behaviour works for name, type, resname, segid, altLoc, resid,
resnum and bynum selections!  The old behaviour will still work,
but we feel this should save a lot of typing!

## Going strong towards Python 3

We are planning to support Python 3 as well in the future. This release has
started our process of adapting MDAnalysis to be compatible with Python 2.7 and
Python 3.

## Others

This release contains other performance enhancements and fixes. For a detailed
list see the [release notes]()
