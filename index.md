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
[NumPy](https://numpy.scipy.org/) arrays. This provides a flexible and
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

MDAnalysis can be [**easily
installed**]({{site.baseurl}}/pages/installation_quick_start/) with
its dependencies using the ``pip`` or ``conda`` package managers;
MDAnalysis is also available in some recent Linux distributions.

All **source code** is available under the
[GNU General Public License, version 2](https://www.gnu.org/licenses/gpl-2.0.html)
(or any later version at your choice) from
[github.com/MDAnalysis/mdanalysis](https://github.com/MDAnalysis/mdanalysis)
and the Python Package index
[pypi.org/project/MDAnalysis](https://pypi.org/project/MDAnalysis).



### Participating

Ask **questions** on the [{{ site.mailinglists.discussion.name }}
mailing list]({{ site.mailinglists.discussion.url }}) and **be part of the
conversation**. You can also join the [{{ site.discord.name }}
Discord Server]({{ site.discord.url }})  to talk with other users and
developers. (In order to join our Discord server, use the invitation link
[{{ site.discord.invite }}]({{ site.discord.invite }}).)


Please report **bugs** or **enhancement requests** through the [Issue
Tracker]({{ site.github.issues }}). 

MDAnalysis is **open source** and welcomes *your* contributions. [Fork
the repository on
GitHub](https://github.com/MDAnalysis/mdanalysis#fork-destination-box)
and submit a pull request. Participate on the [{{
site.mailinglists.developer.name }}
mailing list]({{ site.mailinglists.developer.url }}).



### Supporting

If you like MDAnalysis and want to support our mission, please
consider making a donation to support our efforts.

{{ site.numfocus.donate_button }}

(Donations are made through [our fiscal sponsor][], [NumFOCUS][], which is
a 501(c)(3) non-profit charity in the United States; as such,
donations to NumFOCUS are tax-deductible as allowed by law.  As with
any donation, you should consult with your personal tax adviser or the
IRS about your particular tax situation.)

[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)
Additionally, consider displaying our badge in your projects that use
MDAnalysis. We provide several [embedding markup examples]({{ site.baseurl }}/pages/citations/#powered-by-mdanalysis).

<a href="https://github.com/MDAnalysis/mdanalysis"><img style="position: absolute; top:
0; right: 0; border: 0;"
src="https://camo.githubusercontent.com/a6677b08c955af8400f44c6298f40e7d19cc5b2d/68747470733a2f2f73332e616d617a6f6e6177732e636f6d2f6769746875622f726962626f6e732f666f726b6d655f72696768745f677261795f3664366436642e706e67"
alt="Fork me on GitHub"
data-canonical-src="https://s3.amazonaws.com/github/ribbons/forkme_right_gray_6d6d6d.png"></a>


[NumFOCUS]: https://www.numfocus.org
[our fiscal sponsor]: {{site.baseurl}}/about#partners
