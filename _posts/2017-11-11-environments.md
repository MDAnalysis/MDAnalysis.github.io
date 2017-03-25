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

On this post I will explain several ways to keep maintain different versions of
the same program on a computer. For this I'm introducing several ways to create
isolated environments. The guide will mostly focus on python packages but large
parts can independently be applied to other languages and programs as well.

# Conda Environments

[conda]() is a general package manager for scientific applications. It is mostly
used for python packages but the system can be used with any programs. The
[conda-forge] community also provides a large collection of scientific software
for python, R and perl. Conda should be your first choice to manage different
software versions.

In this guide we will only concentrate on creating and managing environments
with conda. For more information on general installation of package please refer
to the [official documentation]().

To create a new environment for your next project that uses numpy in version
0.10.0 run:

{% highlight bash %}
conda create -n awesome_project numpy=0.10.0 -y
{% endhighlight %}

This only create a new environment 

{% highlight bash %}
source activate awesome_project
{% endhighlight %}

To list your environments

{% highlight bash %}
conda env list
{% endhighlight %}

To remove an environment

{% highlight bash %}
conda remove -all -n awesome_project
{% endhighlight %}


More information about conda environments can be found in
the [official documentation]().

# Linux Environment Modules

You likely already use environment modules if you ever used a super computer.
They have been a part of most unix systems for decades now. To see if you
currently have any modules installed run.

{% highlight bash %}
module avail
{% endhighlight %}

Since the command is called `module` we will refer to environments as modules in
this section. You can show information about a module with

{% highlight bash %}
module info <module-name>
{% endhighlight %}

Modules work by changing environment variables like `PATH` and `LDD_PATH` to
point to software packages you want to use. This means that you should be
familiar with common unix environment variables and installing software in
custom locations. The advantage of this approach though is that you don't have
to install anything new and it can be easier to use software that hasn't been
packaged with conda.


## Create a Environment Module

To start we will create some example environments and later we will show how to
install different versions of vmd. Modules are defined using module-files. So we
first need to create a folder and change into it.

{% highlight bash %}
mkdir -p ~/modules/modulefiles
cd ~/modules/modulefiles
{% endhighlight %}

We also need to tell the module system that it should look for module-files in
that folder. Assuming you are using bash run

{% highlight bash %}
echo "export MODULEPATH=$MODULEPATH:~/modules/modulefiles" >> .bashrc
source .bashrc
{% endhighlight %}

Our first environment will just set a new environment variable called `TEST`

{% highlight bash %}
{% endhighlight %}

More information on linux modules can be found with `man module` and `man modulefile`. 


## Installing software into an environment

## Handling conflicts

The module command allows you to have several modules active at the same time.
To avoid potential version conflicts you have to tell specify in the module file
which other modules are incompatible.


# Python Virtual Environments

Virtual environments will only work for python packages.

# Automatically Change Environment Based On Folder

After you created an isolated environment for a project it would be nice if it
would be activated automatically when ever you work on the project. If you use
bash/zsh this can be done on a per folder basis.




