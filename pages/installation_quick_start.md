---
layout: page
title: Installation Quick Start
order: 2
---

MDAnalysis offers two methods to install the released version. For most
users we recommend the [*Conda*](#conda) installation method.

## Conda ##

If you don't have [conda][] installed yet, follow the [installation
instructions for conda][].

To install the lastest stable release with [conda][conda] do:

{% highlight bash %}
conda config --add channels conda-forge
conda install mdanalysis
{% endhighlight %}

To upgrade to the latest stable release.

{% highlight bash %}
conda update mdanalysis
{% endhighlight %}

To [run the test cases][run_tests] and examples, also install the unit tests (about 53 MiB
in size):

{% highlight bash %}
conda install MDAnalysisTests
{% endhighlight %}

The conda packages currently only support serial calculations. If you
plan to use the parallel [OpenMP][OpenMP] algorithms you need to
install MDAnalysis from the [Python Package
Index](#python-package-index) and have a working OpenMP installation.

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

For more details on installation and alternative ways to install MDAnalysis see [Installing
MDAnalysis in the {{ site.docs.userguide.name }}]({{ site.docs.userguide.url }}/stable/installation.html).

If you have questions with the installation, please ask on
[{{site.mailinglists.discussion.name}}]({{site.mailinglists.discussion.url}}).

[pip]: https://pip.pypa.io/en/latest/
[conda]: https://conda.io/
[installation instructions for conda]: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
[OpenMP]: https://www.openmp.org/
[run_tests]: {{ site.docs.userguide.url }}/stable/installation.html#testing
