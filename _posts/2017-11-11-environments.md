---
layout: post
title: Working with different package versions
---

In science projects often take years until they are complete. During that time
your the libraries you are using will be updated and most have new features you
like to use. But sometimes the upgrade of a library means that old programs and
scripts will stop working. It can also be that you need a new library version
for another project you are working on. The solution to this problem is to
install the same programs multiple times.

In this post I will explain how conda and python virtual envs can be used to
manage different package versions.

# Conda Environments

[conda]() is a general package manager for scientific applications. It is mostly
used for python packages but the system can be used with any programs. The
[conda-forge] community also provides a large collection of scientific software
for python, R and perl. Conda should be your first choice to manage different
software versions.

In this guide we will only concentrate on creating and managing environments
with conda. For more information on general installation of package please refer
to the [official documentation]().

As a preparation add our conda-channel to your configuration

{% highlight bash %}
conda config  --add channels MDAnalysis
{% endhighlight %}

To create a new environment for your next project that uses MDAnalysis in version
0.15.0 run:

{% highlight bash %}
conda create -n mda-15 MDAnalysis=0.15.0 -y
{% endhighlight %}

This only create a new environment. You still have to activate it to be able to
use it.

{% highlight bash %}
source activate mda-15
{% endhighlight %}

To list your environments

{% highlight bash %}
conda env list
{% endhighlight %}

A nice feature of using conda-environments is that they are easy to share with
colleagues or transferred to other computers. To store the state of the
environment we created in a file called `mda-15-environment`

{% highlight bash %}
conda list --explicit --md5 -n mda-15 > mda-15-environment
{% endhighlight %}

You can now copy this file to a colleague or onto another computer. The first 3
lines also contain instructions how this file can be used with conda.

{% highlight bash %}
# This file may be used to create an environment using:
# $ conda create --name <env> --file <this file>
# platform: linux-64
{% endhighlight %}

Instead of just creating environments for specific versions of software packages
you can of course also create a new environment for each project you are working
on.


More information about conda environments can be found in
the [official documentation]().

# Python Virtual Environments

Virtual environments will only work for python packages.

# Automatically Change Environment Based On Folder

After you created an isolated environment for a project it would be nice if it
would be activated automatically when ever you work on the project. If you use
bash/zsh this can be done on a per folder basis.




