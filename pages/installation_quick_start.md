---
layout: page
title: Installation Quick Start
order: 2
---

MDAnalysis offers multiple installation methods. For most
users we recommend [mamba][] installation method, a faster
drop-in replacement for [conda][] .

## Mamba ##

If you don't have [conda][] installed yet, follow the [installation
instructions for conda][]. Following that, install [mamba][] with:

{% highlight bash %}
conda install -c conda-forge mamba
{% endhighlight %}

Then, install the latest stable release of MDAnalysis with:

{% highlight bash %}
mamba install -c conda-forge mdanalysis
{% endhighlight %}

To upgrade use:

{% highlight bash %}
mamba update mdanalysis
{% endhighlight %}

To [run the test cases][run_tests] and examples, also install the unit tests (about 90 MB
in size):

{% highlight bash %}
mamba install -c conda-forge MDAnalysisTests
{% endhighlight %}

MDAnalysis installed via mamba supports only serial calculations. 
If you need OpenMP-based parallelism, install MDAnalysis via [pip](#python-package-index) 
and ensure you have a working [OpenMP][] installation.

## Python Package Index ##

To install the latest stable release with
[pip][pip] (which should be available in all Python installations) do:

{% highlight bash %}
pip install --upgrade MDAnalysis
{% endhighlight %}

To [run the test cases][run_tests] and examples, also install the unit tests (about 53 MiB
in size):

{% highlight bash %}
pip install --upgrade MDAnalysisTests
{% endhighlight %}

## More ##

For more installation options and additional details see [Installing
MDAnalysis in the {{ site.docs.userguide.name }}]({{ site.docs.userguide.url }}/stable/installation.html).

If you have questions with the installation, please ask on
[{{site.mailinglists.discussion.name}}]({{site.mailinglists.discussion.url}}).

[pip]: https://pip.pypa.io/en/latest/
[mamba]:https://anaconda.org/conda-forge/mamba
[conda]: https://conda.io/
[installation instructions for conda]: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
[OpenMP]: https://www.openmp.org/
[run_tests]: {{ site.docs.userguide.url }}/stable/installation.html#testing
