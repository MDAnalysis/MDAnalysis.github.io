---
layout: post
title: NumFOCUS small grant for Python 3 support
---

We received a [small development grant][grant-post] by NumFOCUS to finally add
python 3 support to MDAnalysis. To do this Richard Gowers and Tyler Reddy will
visit Oliver Becksteins lab in the summer for a week. 

Although about 80% of code passes unit tests in Python 3, we urgently need to
close the remaining 20% gap in order to support our user base and to safeguard
the long term viability of the project. MDAnalysis started almost 10 years ago
when Python was around version 2.4 and interfacing with existing C code was
mostly done with writing C-wrappers that directly used CPython. This legacy code
has hampered a speedy full transition to Python 3 and consequently MDAnalysis
lags behind the rest of the scientific Python community in fully supporting
Python 3.

In the meantime we are busy porting our last python 2.7 only c-extension, the
DCD reader and writer, to Cython. We now have a working cython version that can
be used without MDAnalysis, similar to our XTC and TRR readers, only cleaning up
the code and writing documentation is left to do. You can check our
progress [here][dcd-pr].


[grant-post]: https://www.numfocus.org/blog/numfocus-awards-small-development-grants-to-projects/
[dcd-pr]: https://github.com/MDAnalysis/mdanalysis/pull/1155
