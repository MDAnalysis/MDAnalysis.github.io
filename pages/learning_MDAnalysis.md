---
layout: page
title: Learning MDAnalysis
---

Once you had a look at the 
[basic example]({{ site.baseurl }}/pages/basic_example) 
you might want to learn more about how to use
MDAnalysis. MDAnalysis is primarily a library that helps you to build
your own tools but it also works very well for **interactive data
exploration** of MD data in [IPython](http://ipython.org/), in
particular within [Jupyter notebooks](https://jupyter.org/) and in
conjunction with [pandas](http://pandas.pydata.org/). MDAnalysis is
well suited for a *rapid development* approach.

The resources below should help you to quickly find out to best use
MDAnalysis for your own specific uses.


## Tutorials ##

The [MDAnalysis
Tutorial](http://www.mdanalysis.org/MDAnalysisTutorial/) serves as an
introduction to the library and there are [other
tutorials]({{site.github.wiki}}/Tutorials)
available, too.

[Interactive Jupyter
notebooks](http://nbviewer.jupyter.org/github/MDAnalysis/binder-notebook/tree/master/notebooks/)
show how to accomplish specific tasks (including visualizing
trajectories with [nglview](http://nglviewer.org/nglview/latest/));
these notebooks can be run in the cloud on Binder (click the "launch
binder" button to start a notebook server).

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/MDAnalysis/binder-notebook/master?filepath=notebooks)


## Documentation ##

See the [Online Documentation]({{site.pypi.docs}})
for more information on how to use MDAnalysis and the available
documentation on the Wiki. The [paper on
MDAnalysis]({{ site.baseurl }}/pages/citations#Gowers2016) contains a
high-level description of the structure and philosophy of the library
together with examples of its use.

The [FAQ]({{ site.github.wiki }}/FAQ) contains a 
growing list of specific (frequently asked) questions and answers.

## Mailing list ##

Finally, you can also ask for advice or help on the
[mdnalysis-discussion mailing
list](http://groups.google.com/group/mdnalysis-discussion). If you
find *bugs* or want to *request enhancements* please [file a
report]({{site.github.wiki}}/ReportingProblems)
in the [Issue Tracker]({{sitemap.github.issues}}).

## Videos ##

### Introductory ###

#### The universe as balls and springs: molecular dynamics in Python
@lilyminium's talk at [PyCon AU 2019](https://2019.pycon-au.org/) *The universe as balls and
springs: molecular dynamics in Python* gives a general introduction to
molecular dynamics and shows how to use MDAnalysis (and other tools
such as [OpenMM](http://openmm.org/), [nglviewer](nglviewer.org/nglview/latest/),
[pandas](https://pandas.pydata.org/),
[plotly](https://pandas.pydata.org/)). If you want to better
understand what MD simulations are and how scientists can make use of
the vast Python eco-system to analyze (and run) MD simulations, start here:

<div class="js-video">
	<iframe src="https://www.youtube.com/embed/X5umNQDqfqQ" frameborder="0"
	allowfullscreen class="video"></iframe>
</div>

#### MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations

@orbeckst's talk at [SciPy 2016](http://scipy2016.scipy.org/) provides
an introduction to the MDAnalysis library, its uses, and underlying philosophy:

<div class="js-video">
	<iframe src="https://www.youtube.com/embed/zVQGFysYDew" frameborder="0"
	allowfullscreen class="video"></iframe>
</div>

Also read the paper [MDAnalysis: A Python package for the rapid
analysis of molecular dynamics
simulations](http://conference.scipy.org/proceedings/scipy2016/oliver_beckstein.html)
which adds detail to the concepts outlined in this talk.



### Intermediate ###

#### Looking at molecules using Python
@jbarnoud presented at the PyGrunn 2017 conference _Looking at
molecules using Python_ where he shows how to use a whole range of
MDAnalysis from the simple to the advanced in Jupyter notebooks (he
also shows off [nglview](http://nglviewer.org/nglview/latest/) for
visualization and [datreant](http://datreant.org) for organizing his
data):

<div class="js-video">
	<iframe src="https://www.youtube.com/embed/RWgt1WMwMUs" frameborder="0"
	allowfullscreen class="video"></iframe>
</div>


