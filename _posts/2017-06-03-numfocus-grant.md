---
layout: post
title: NumFOCUS small grant for Python 3 support
---

We have generously been awarded a [small development grant][grant-post] by
NumFOCUS to finally add
Python 3 support to MDAnalysis. To do this [Richard Gowers][rg] and [Tyler Reddy][tr] will
be hosted at [Oliver Beckstein's][ob] [lab at Arizona State University][beckstlab]
in the summer for a week of hacking.

MDAnalysis started almost 10 years ago
when Python was around version 2.4 and interfacing with existing C code was
mostly done with writing C-wrappers that directly used CPython. This legacy code
has hampered a speedy full transition to Python 3 and consequently MDAnalysis
lags behind the rest of the scientific Python community in fully supporting
Python 3.
Although about 80% of code passes unit tests in Python 3, we urgently need to
close the remaining 20% gap in order to support our user base and to safeguard
the long term viability of the project.

In the meantime we are busy porting our last Python 2.7 only C-extension, the
DCD Reader and Writer, to Cython. We now have a working cython version that can
be used without MDAnalysis, similar to our XTC and TRR readers, only cleaning up
the code and writing documentation is left to do. You can check our
progress [here][dcd-pr].

[rg]: https://github.com/richardjgowers
[tr]: https://github.com/tylerjereddy
[ob]: https://github.com/orbeckst
[grant-post]: https://www.numfocus.org/blog/numfocus-awards-small-development-grants-to-projects/
[beckstlab]: https://becksteinlab.physics.asu.edu/
[dcd-pr]: https://github.com/MDAnalysis/mdanalysis/pull/1155