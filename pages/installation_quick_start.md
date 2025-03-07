---
layout: page
title: Installation Quick Start
order: 2
---

MDAnalysis can be installed using [mamba][], a faster drop-in replacement for [conda][], or [pip][].

## mamba/conda ##

If you don't have [mamba][] installed you can follow the [mamba installation instructions][]. 
To install MDAnalysis run:

{% highlight bash %}
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

To install MDAnalysis with [pip][] run:

{% highlight bash %}
pip install --upgrade MDAnalysis
{% endhighlight %}

To install [test cases][] run:

{% highlight bash %}
pip install --upgrade MDAnalysisTests
{% endhighlight %}

## More ##

For more details about the installation see the [installation instructions in the {{ site.docs.userguide.name }}]({{ site.docs.userguide.url }}/stable/installation.html).

If you have questions regarding the installation, please ask on
[{{site.mailinglists.discussion.name}}]({{site.mailinglists.discussion.url}}).

[pip]: https://pip.pypa.io/en/latest/
[mamba]:https://anaconda.org/conda-forge/mamba
[conda]: https://conda.io/
[mamba installation instructions]: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
[OpenMP]: https://www.openmp.org/
[test cases]: {{ site.docs.userguide.url }}/stable/installation.html#testing
