---
layout: page
title: Learning MDAnalysis
order: 3
---

MDAnalysis is a powerful Python library for analyzing MD simulations. While primarily designed to help you build custom analysis tools, it also supports interactive data exploration in environments like [IPython](http://ipython.org/) and [Jupyter notebooks](https://jupyter.org/), especially when combined with [pandas](http://pandas.pydata.org/). This makes MDAnalysis an excellent choice for rapid prototyping and exploratory analysis.

Whether you're new to MDAnalysis or looking to deepen your expertise, this page will guide you through our learning resources.

## Tutorials ##

If you are new to MDAnalysis, start with:

- [{{ site.docs.quickstart.name }}]({{ site.docs.quickstart.url }}) – A quick introduction to MDAnalysis.
- [{{ site.docs.userguide.name }}]({{ site.docs.userguide.url }}) – Detailed documentation with in-depth tutorials.

For full documentation, visit the [Documentation page]({{ site.docs.url }}).

You can ask for advice or help on [{{ site.mailinglists.discussion.name }}]({{
site.mailinglists.discussion.url }}). If you find *bugs* or want to *request enhancements* please [file a report]({{site.github.wiki }}/ReportingProblems) in the [Issue Tracker]({{sitemap.github.issues }}).

## Videos ##

The following videos, presented by core developers at conferences, highlight various aspects of MDAnalysis and demonstrate its use in research.

In addition to conference talks, we have recorded and published a number of MDAnalysis workshop sessions on our [YouTube channel](https://www.youtube.com/channel/UC3TCuK-z_bJNdwWCvsH9D3Q), which cover hands-on tutorials and more in-depth discussions of the library's features. 

### Introductory ###

#### The universe as balls and springs: molecular dynamics in Python
@lilyminium's talk at [PyCon AU 2019](https://2019.pycon-au.org/) *The universe as balls and
springs: molecular dynamics in Python* gives a general introduction to
molecular dynamics and shows how to use MDAnalysis (and other tools
such as [OpenMM](http://openmm.org/), [nglviewer](https://nglviewer.org/nglview/latest/),
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


