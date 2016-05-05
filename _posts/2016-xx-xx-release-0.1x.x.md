---
layout: post
title: Release 0.15.0
---

We have just release MDAnalysis version 0.15.0. This release contains new
features as well as bug fixes. Highlights are listed below but for more details
see the [release notes](https://github.com/MDanalysis/mdanalysis/wiki/...).

We also had a lot of contributions from GSoC applicants. Thanks to our GSoC
students @fiona-naughton and @jdetle as well as the other applicants @Saxenauts,
@Endle, @abhinavgupta94 and @pedrishi.

# Upgrade

You can upgrade with `pip install --upgrade MDAnalysis`

# Noticable Changes

## Revamped Contact Analysis

Contact Analysis has been completely rewritten in the new *Contacts* class. This
class offers a standard native contact analysis as well as a contact analysis
developed by [Best & Hummer][best-hummer-paper]. A Q1-Q2 analysis is not
available directly as *q1q2*. We have also made the *Contacts* extendable so
that you can pass it your own cut-off functions, the *q1q2* analysis is actually
only a wrapper of *Contacts* that makes use of this flexibility. More
information can be found in the [documentation][contacts-docs].

The old `ContactAnalysis1` and `ContactAnalysis` classes will be removed in the
next release.

## RMSD Calculation

The [*rmsd*][rmsd-docs] function now doesn't super position the given
coordinates by default. The coordinates aren't changed now by default, instead
you can control it with the new *center* and *superposition* keywords.

[contacts-docs]: http://www.mdanalysis.org/mdanalysis/documentation_pages/analysis/contacts.html
[best-hummer-paper]: http://www.pnas.org/content/110/44/17874
[rmsd-docs]: http://www.mdanalysis.org/mdanalysis/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.rmsd
