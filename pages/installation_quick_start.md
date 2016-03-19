---
layout: page
title: Installation Quick Start
---

# Python Package Index

To install the latest stable release with
[pip][pip] do:

{% highlight bash %}
pip install --upgrade MDAnalysis
{% endhighlight %}

To run the test cases and examples, also install the unit tests (about 20 MiB
in size):

{% highlight bash %}
pip install --upgrade MDAnalysisTests
{% endhighlight %}

# Conda

To install the lastest stable release with [conda][conda] do:

{% highlight bash %}
conda config --add channels MDAnalysis
conda install mdanalysis
{% endhighlight %}

To upgrade to the latest stable release.

{% highlight bash %}
conda update mdanalysis
{% endhighlight %}

The conda packages currently only support serial calculations. If you
plan to use the parallel [OpenMP][OpenMP] algorithms you need to install
MDAnalysis with [pip][pip] and have a working OpenMP installation.

# More

For more details on installation and alternative ways to install MDAnalysis
(e.g. through your distribution's package manager, see [Installing
MDAnalysis]({{site.github.wiki}}/Install)).

If you have questions with the installation, please ask on the
[{{site.mailinglists.discussion.name}}]({{site.mailinglists.discussion.url}})
mailing list.

[pip]: http://www.pip-installer.org/en/latest/index.html
[conda]: http://conda.pydata.org/docs/
[OpenMP]: http://openmp.org
