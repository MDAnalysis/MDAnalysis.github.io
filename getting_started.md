---
layout: page
title: Getting Started 
order: 3
---

MDAnalysis is a Python library for analyzing molecular dynamics (MD) simulations. This page will guide you through quickly installing MDAnalysis, exploring a basic example, and accessing learning resources to get started.

## Installation Instructions

Our [Installation Quick Start]({{ site.baseurl }}/pages/installation_quick_start) guide contains instructions to get MDAnalysis up and running in a few minutes. Try this guide first. 

For a more thorough installation with advanced options, including dependencies and environment configurations, see the [installation instructions in the User Guide]({{ site.docs.userguide.url }}/stable/installation.html). 

## Basic Example

Once MDAnalysis is installed, you can load a trajectory and perform a simple analysis. A typical usage pattern is to iterate through a trajectory and analyze
coordinates for every frame.

In the following example, the end-to-end distance of a protein and the radius of gyration of the backbone atoms are calculated:

<div class="wide-code">
    {% highlight python %}
    import MDAnalysis
    from MDAnalysis.tests.datafiles import PSF, DCD   # test trajectory
    import numpy.linalg

    # load trajectory and topology into a Universe
    u = MDAnalysis.Universe(PSF,DCD)  
    # from the 4AKE segid, select
    # the first atom named N and the last atom named C
    nterm = u.select_atoms('segid 4AKE and name N')[0]
    cterm = u.select_atoms('segid 4AKE and name C')[-1]
    # select the backbone atoms (AtomGroup)
    bb = u.select_atoms('protein and backbone') 

    for ts in u.trajectory: # iterate through all frames
        r = cterm.position - nterm.position # end-to-end vector from atom positions
        d = numpy.linalg.norm(r) # end-to-end distance
        rgyr = bb.radius_of_gyration() # method of AtomGroup
        print("frame = {0}: d = {1} A, Rgyr = {2} A".format(
              ts.frame, d, rgyr))
    {% endhighlight %}
</div>

## Learning Resources

To find out what else you can do, head over to [Learning
MDAnalysis]({{ site.baseurl }}/pages/learning_MDAnalysis) to explore
tutorials and documentation.

If you have questions, visit our [Community]({{ site.baseurl }}/pages/community) page to learn about available discussion channels. Happy coding!

