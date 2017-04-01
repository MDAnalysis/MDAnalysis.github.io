---
layout: post
title: Managing software versioning using Conda environments
---

Research projects can often take months to years to complete, however the
precise version of software they use will often have a shorter lifetime than
this.  These new versions of software will often include new features which
might be of great use, but they might also introduce changes which break
your existing work and introduce compatibility issues with other pieces
of software.  So whilst for existing projects we might wish to freeze
all pieces of software installed, for newer projects we instead want to
use the most up-to-date versions, leaving us needing to install multiple
versions of multiple different pieces of software.

In this post we will try to explain how conda and Python virtual environments
can be used to precisely manage the software used across many independent
research projects.


# Conda Environments

[Conda](https://conda.io/docs/index.html) is a general package manager for scientific
applications. It is mostly used for Python packages but the system can be used with
any programs. The [conda-forge] community also provides a large collection of scientific software
for Python, R and perl. Conda should be your first choice to manage different
software versions.

In this guide we will concentrate only on creating and managing environments
with conda. For more information on general installation of package please refer
to the [official documentation]().

Software is made available through different conda channels, which each act as a
source for different software.  When attempting to install packages into a conda
environment, these channels are searched.  In this post we will be using the
MDAnalysis channel which you can add to your configuration like so:

{% highlight bash %}
conda config  --add channels MDAnalysis
{% endhighlight %}

For each research project, it is advised that you create a new environment so that
the software used in each project does interfere across different projects.
To create a new environment for your next project that uses MDAnalysis in version
0.15.0 run:

{% highlight bash %}
conda create -n myproject MDAnalysis=0.15.0 -y
{% endhighlight %}

This has created a new software environment called `myproject` but has not affected
anything currently!  To have access to all the software installed within it we
must first activate it

{% highlight bash %}
source activate myproject
{% endhighlight %}

To list your available environments

{% highlight bash %}
conda env list
{% endhighlight %}

A nice feature of using conda-environments is that they are easy to share with
colleagues or transferred to other computers.  This allows all team members on
a project to use an identical set of software and makes your research projects
to be reproducible.  To store the state of the environment we created in a file
called `myproject-environment`

{% highlight bash %}
conda list --explicit --md5 -n myproject > myproject-environment
{% endhighlight %}

You can now copy this file to a colleague or onto another computer. The first 3
lines also contain instructions how this file can be used with conda.

{% highlight bash %}
# This file may be used to create an environment using:
# $ conda create --name <env> --file <this file>
# platform: linux-64
{% endhighlight %}

More information about conda environments can be found in
the [official documentation]().

# Python Virtual Environments

Virtual environments will only work for python packages.

# Automatically Change Environment Based On Folder

After you created an isolated environment for a project it would be nice if it
would be activated automatically when ever you work on the project. If you use
bash/zsh this can be done on a per folder basis.




