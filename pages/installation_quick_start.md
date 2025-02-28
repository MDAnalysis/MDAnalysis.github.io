---
layout: page
title: Installation Quick Start
order: 2
---

MDAnalysis can be installed using [mamba][] (recommended) or [pip][].

## Mamba ##

If you have [conda][] or [mamba][] installed, run:

{% highlight bash %}
conda install -c conda-forge mamba
mamba install -c conda-forge mdanalysis
{% endhighlight %}

To install [test cases][] cases (about 90 MB
in size) run: 

{% highlight bash %}
mamba install -c conda-forge MDAnalysisTests
{% endhighlight %}

MDAnalysis via [mamba][] supports only serial calculations. 
For OpenMP-based parallelism, use [pip][] and ensure you have 
a working [OpenMP][] installation.

## Python Package Index ##

To install with [pip][] run:

{% highlight bash %}
pip install --upgrade MDAnalysis
{% endhighlight %}

To install [test cases][] run:

{% highlight bash %}
pip install --upgrade MDAnalysisTests
{% endhighlight %}

## More ##

For detailed installation methods, including setting up [mamba][] from scratch, see the [Installation instuctions in the {{ site.docs.userguide.name }}]({{ site.docs.userguide.url }}/stable/installation.html).

If you have questions regarding the installation, please ask on
[{{site.mailinglists.discussion.name}}]({{site.mailinglists.discussion.url}}).

[pip]: https://pip.pypa.io/en/latest/
[mamba]:https://anaconda.org/conda-forge/mamba
[conda]: https://conda.io/
[installation instuctions]: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
[OpenMP]: https://www.openmp.org/
[test cases]: {{ site.docs.userguide.url }}/stable/installation.html#testing
