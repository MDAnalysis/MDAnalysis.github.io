---
layout: default
title: MDAnalysis
---


<img src="{{ site.baseurl }}/public/mdanalysis-logo_square.png"
style="float: right" alt="MDAnalysis" width="30%"/>

**MDAnalysis** is an object-oriented Python library to analyze
trajectories from molecular dynamics (MD) simulations in [many popular
formats]({{site.pypi.docs}}/documentation_pages/coordinates/init.html#id1). It
can write most of these formats, too, together with [atom
selections]({{site.pypi.docs}}/documentation_pages/selections_modules.html#selection-exporters)
suitable for visualization or native analysis tools.

MDAnalysis allows one to read particle-based trajectories (including
individual coordinate frames such as biomolecules in the PDB format)
and access the atomic coordinates through
[NumPy](http://numpy.scipy.org/) arrays. This provides a flexible and
relatively fast framework for complex analysis tasks. In addition,
powerful atom
[selection commands]({{site.pypi.docs}}/documentation_pages/selections.html)
are implemented. Trajectories can also be manipulated (for instance,
fit to a reference structure) and written out. The
[basic example]({{ site.baseurl }}/pages/basic_example) demonstrates some
of these features.

Read more:

* [installation quick start]({{ site.baseurl }}/pages/installation_quick_start)
* [learning MDAnalysis]({{ site.baseurl }}/pages/learning_MDAnalysis)

Also, check out the [blog]({{ site.baseurl }}/blog) or subscribe to our 
[news feed]({{ site.baseurl }}/{{site.feed.path}}) to follow development
updates and events.

### Availability

All source code is available under the
[GNU General Public License, version 2](https://www.gnu.org/licenses/gpl-2.0.html)
(or any later version at your choice) from
[github.com/MDAnalysis/mdanalysis](https://github.com/MDAnalysis/mdanalysis)
and the Python Package index
[pypi.python.org/pypi/MDAnalysis](http://pypi.python.org/pypi/MDAnalysis).


### Participating

Ask **questions** on the [{{site.mailinglists.discussion.name}}
mailing list]({{site.mailinglists.discussion.url}}) and join the
discussion.

Please report **bugs** or **enhancement requests** through the [Issue
Tracker]({{site.github.issues}}). 

MDAnalysis is **open source** and welcomes *your* contributions. [Fork
the repository on
GitHub](https://github.com/MDAnalysis/mdanalysis#fork-destination-box)
and submit a pull request. Participate on the [{{site.mailinglists.developer.name}}
mailing list]({{site.mailinglists.developer.url}}).

<a href="https://github.com/MDAnalysis/mdanalysis"><img style="position: absolute; top:
0; right: 0; border: 0;"
src="https://camo.githubusercontent.com/a6677b08c955af8400f44c6298f40e7d19cc5b2d/68747470733a2f2f73332e616d617a6f6e6177732e636f6d2f6769746875622f726962626f6e732f666f726b6d655f72696768745f677261795f3664366436642e706e67"
alt="Fork me on GitHub"
data-canonical-src="https://s3.amazonaws.com/github/ribbons/forkme_right_gray_6d6d6d.png"></a>
