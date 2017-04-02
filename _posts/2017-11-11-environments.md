---
layout: post
title: Managing software versioning using Conda environments
---

Research projects can often take months to years to complete, but the precise
version of software they use will often have a shorter lifetime than this. These
new versions of software will often include new features which might be of great
use, but they might also introduce changes which break your existing work and
introduce compatibility issues with other pieces of software. So whilst for
newer projects we might want to use the most up-to-date versions of software,
for existing projects we want to be able to freeze the versions of software that
are in use. This leads to us needing to install and manage multiple versions of
software across our various research projects.

With the upcoming release of 0.16.0 of MDAnalysis, alongside various
improvements, we are also introducing some changes which could break existing
code. In this post we will explain how [conda environments](#conda-environments)
or [Python virtual environments](#python-virtual-environments) can be used to
manage this transition, allowing you to finish existing projects with version
0.15.0, while also enjoying the benefits provided in version 0.16.0.

# Conda Environments

[Conda](https://conda.io/docs/index.html) is a general package manager for
scientific applications. It is mostly used for Python packages but the system
can be used with any program. The [conda-forge](https://conda-forge.github.io/)
community also provides a large collection of scientific software for Python, R
, Perl and others.

In this guide we will concentrate only on creating and managing environments
with conda. For more general information on installing conda please refer
to the [official documentation](https://conda.io/docs/using/pkgs.html).

Software is made available through different conda channels, which each act as a
source for different software.  When attempting to install packages into a conda
environment, these channels are searched.  In this post we will be using the
MDAnalysis channel which you can add to your configuration like so:

{% highlight bash %}
conda config  --add channels MDAnalysis
{% endhighlight %}

For each research project, it is advised that you create a new environment so that
the software used in each project does not interfere across different projects.
To create a new environment for your next project that uses MDAnalysis in version
0.15.0 run:

{% highlight bash %}
conda create -n myproject mdanalysis=0.15.0 -y
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
colleagues or transferred to other computers. This allows all collaborators on a
project to use an identical set of software and makes your research projects
reproducible. To store the state of the environment we created in a file called
`myproject-environment`

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
the [official documentation](https://conda.io/docs/using/envs.html).

# Python Virtual Environments

Like conda, virtual environments managed with virtualenv allow you to use
different versions of python and python packages for your different project. On
the contrary to conda, however, virtualenv is not a general-purpose package
manager; it leverages what is available on your system, and let you install
python packages using pip.

To use virtual environments you have to install the virtualenv package first.
This can be done with either pip or the package manager of your system:

{% highlight bash %}
pip install virtualenv
# or on ubuntu
sudo apt install virtualenv
# or on fedora
sudo dnf install python-virtualenv
{% endhighlight %}

Virtual environments can be created per project directory.

{% highlight bash %}
cd myproject
# if you want to use the default python interpretor
virtualenv myproject-env
# or if you want to explicitly use python 2.7
virtualenv -p /usr/bin/python2.7 myproject-env
# or if you want to use python 3.5
virtualenv -p /usr/bin/python3.5 myproject-env
# or maybe you want to use pypy
virtualenv -p /usr/bin/pypy myproject-env
{% endhighlight %}

This will create a new folder `myproject-env`. This folder contains the virtual
environment and all packages you have installed in it. To activate it in the
current terminal run:

{% highlight bash %}
source myproject-env/bin/activate
{% endhighlight %}

Now you can install packages via `pip` without affecting your global
environment. The packages you install when the environment is activated will be
available in terminal sessions that have the environment activated. You can
deactivate the virtual environment by running

{% highlight bash %}
deactivate
{% endhighlight %}

The [`virtualenvwrapper` package](https://virtualenvwrapper.readthedocs.io)
makes virtual environment easier to use. It provides some very useful features:

* it organizes the virtual environment in a single directory, avoiding to have
  them scattered throughout the file system;
* it defines command to easy the creation, deletion, and copy of virtual
  environments;
* it defines a command to activate a virtual environment using its name;
* all commands defined by `virtualenvwrapper` have tab-completion for virtual
  environment names.

To use `virtualenvwrapper` you first need to install it *outside* of a virtual
environment:

{% highlight bash %}
pip install virtualenvwrapper
# or on ubuntu
sudo apt install virtualenvwrapper
# or on fedora
sudo dnf install python-virtualenvwrapper
{% endhighlight %}

Then, you need to have it loaded in your terminal session. Add the following
lines in `~/.bashrc`, they will be executed every time you open a new terminal
session:

{% highlight bash %}
# Decide where to store the virtual environments
export WORKON_HOME=~/Envs
# Make sure the directory exists
mkdir -p ${WORKON_HOME}
# Load virtualenvwrapper
source /usr/local/bin/virtualenvwrapper.sh
{% endhighlight %}

Open a new terminal or run `source ~/.bashrc` to update your session. You can
now create a virtual environment with

{% highlight bash %}
mkvirtualenv my-project
{% endhighlight %}

Like the `virtualenv` command we saw earlier, `mkvirtualenv` lets you choose
your python interpretor with the `-p` option. Regardless of your current
working directory, we created the virtual environment in `~/Envs/` and it is
now loaded in our terminal session.

You can load your virtual environments by running `workon my-project`, and
exit them by running `deactivate`.

Virtual environments, especially with `virtualenvwrapper`, can do much more.
The Hitchhikers Guide to Python has a good
[tutorial](http://docs.python-guide.org/en/latest/dev/virtualenvs/) that
gives a more in depth explanation of virtual environments. The
[`virtualenvwrapper` documentation](https://virtualenvwrapper.readthedocs.io)
is also a good resource to read.
