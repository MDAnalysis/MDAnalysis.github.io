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

Start with the [{{ site.docs.quickstart.name }}]({{
site.docs.quickstart.url }}) when you are new MDAnalysis.

Then browse the [{{ site.docs.userguide.name }}]({{ site.docs.userguide.url }}), which contains detailed documentation for all the important parts of MDAnalysis and many self-contained tutorials.


There are a number of [older tutorials]({{site.github.wiki}}/Tutorials) available, too, although we recommend new users start  with [{{ site.docs.quickstart.name }}]({{
site.docs.quickstart.url }}) and then start reading the [{{ site.docs.userguide.name }}]({{ site.docs.userguide.url }}).


## Documentation ##

The [{{ site.docs.userguide.name }}]({{ site.docs.userguide.url }})
contains installation instructions, the [{{ site.docs.quickstart.name
}}]({{ site.docs.quickstart.url }}), and comprehensive description of
the functionality of MDAnalysis from a user's perspective. **New users
should start here!**

The [{{ site.docs.mdanalysis.name }}]({{ site.docs.mdanalysis.url }})
contains technical information on how to use MDAnalysis. 

The [paper on MDAnalysis]({{ site.baseurl
}}/pages/citations#Gowers2016) contains a high-level description of
the structure and philosophy of the library together with examples of
its use.

The [FAQ]({{ site.github.wiki }}/FAQ) contains a 
growing list of specific (frequently asked) questions and answers.

## Mailing list ##

You can ask for advice or help on the [mdnalysis-discussion mailing
list]({{ site.mailinglists.discussion.url }}). If you find *bugs* or
want to *request enhancements* please [file a report]({{
site.github.wiki }}/ReportingProblems) in the [Issue Tracker]({{
sitemap.github.issues }}).

## Videos ##

The videos listed below were given by core developers at
conferences. They highlight various aspects of MDAnalysis and show how
to use it in a research context.

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


#### BioExcel Webinar: MDAnalysis: Interoperable analysis of biomolecular simulations in Python

In this [BioExcel](https://bioexcel.eu/) webinar, three of the MDAnalysis Core
Developers (@orbeckst, @lilyminium, @IAlibay) summarize the **basics of
MDAnalysis**, show more advanced ways to **hack MDAnalysis** and outline
**future developments**.

<div class="js-video">
	<iframe src="https://www.youtube.com/embed/1Wot83DSt4E" frameborder="0"
	allowfullscreen class="video"></iframe>
</div>

