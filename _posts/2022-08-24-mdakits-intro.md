---
layout: post
title: MDAKits
---

As part of our [CZI EOSS 4 grant]({% post_url 2021-08-31-CZI-EOSS4 %})
we announced our plans to create an **MDAKit ecosystem**. With this post we aim to make
our plans more concrete and solicit feedback from the community.

Beyond the outline provided here, the complete details of our plans can be found in our *white paper* named **MDAKits: Supporting and promoting
the development of community packages leveraging the MDAnalysis
library [v0.1.0]**, available as a PDF at DOI
[10.6084/m9.figshare.20520726.v1][].


## What is an MDAKit?

MDAKits are standalone packages containing code using MDAnalysis components that solves a specific
scientific problems or in some form enhances the functionality of MDAnalysis core
library. An MDAKit can be written by anyone and hosted
anywhere.

A MDAKit can be **registered** in the MDAKits registry. In this case,
it has to fulfill a number of additional requirements such as
open-source licensed, hosted in a version control system, clear
designation of authors/maintainers, documentation, and tests and
continuous integration. Registered MDAKits will be *listed publicly*
and thus be advertised to the whole MDAnalysis community. They will
also be *continuosly tested* against the latest released version and
the current development version of the core MDAnalysis library so that
users and developers have an up-to-date view of the code health of an
MDAKit.


## Why?

The open sharing of code that abides by the basic principles of
[FAIR](https://doi.org/10.15497/RDA00068) (findability, accessibility,
interoperability, and reusability) is essential to robust,
reproducible, and transparent science. However, scientists typically
are not supported in making the substantial effort required to make
software FAIR-compliant, or incentivized with academic recognition or
reward.  

Our goal with MDAKits is to **lower the barrier for researchers to
produce FAIR software.** 

We support developers in creating new packages, guiding them through
the process of achieving best practices and FAIR compliance. At the
same time, we hope to make MDAnalysis useful to a broader community.


## How to develop an MDAKit?

We are producing tools for creating MDAKits to help
developers and we are working on infrastructure to publicize
MDAKits. Our work on MDAKits is an ongoing process but you can now get
started creating your own MDAKit:


### MDAKit project template

Our first tool is the
[**cookiecutter-mdakit**](https://github.com/MDAnalysis/cookiecutter-mdakit),
a [cookiecutter](https://github.com/audreyr/cookiecutter) template
that generates a **skeleton project** that implements our recommended
best practices. With [`cookiecutter`
installed](https://cookiecutter.readthedocs.io/en/latest/installation.html#install-cookiecutter),
execute the following command inside the folder you want to create the
skeletal repository

```bash
cookiecutter gh:MDAnalysis/cookiecutter-mdakit
```

Follow the prompts or hit enter for the default options. 

Then add your own code to the project. Add tests --- you can extend
the example tests in the template that are suitable for
MDAnalysis-based. Commit and push your changes.


<small>(The MDAKit cookiecutter is based off the [Cookiecutter for
Computational Molecular Sciences (CMS) Python
Packages](https://github.com/MolSSI/cookiecutter-cms) by Levi N. Naden
and Jessica A. Nash from the [Molecular Sciences Software Institute
(MolSSI)](http://molssi.org/) and Daniel G. A. Smith of
[ENTOS](https://www.entos.ai/). Thank you!)</small>


### Registering an MDAKit

If you want to **register your MDAKit** then create a pull request to
add a meta data entry `metadata.yaml` to [MDAnalysis:
MDAKits/mdakits/*{YOUR_MDAKIT_NAME}*](https://github.com/MDAnalysis/MDAKits/tree/main/mdakits)
(where you will also find a template to get you started). Your PR will
be reviewed for compliance with the requirements (for right now, see the [white
paper][10.6084/m9.figshare.20520726.v1] for specifics). Once
registered, your MDAKit will be continuously tested.


### Towards publication

The best practices that we encourage MDAKits to fulfill essentially
amount to the majority of the contribution criteria for submissions to
software-focused journals such as the [Journal Open Source
Software](https://joss.theoj.org) (JOSS). We encourage MDAKits to
consider submission to a journal such as JOSS once they meet the
required levels of best practices. We are working towards streamlining
the submission process.


## Give us feedback!

We are looking for **feedback from the community**: please let us know
via our [malinglist or discord]({{ site.baseurl }}/#participating)
channels or via the [MDAKits issue
tracker](https://github.com/MDAnalysis/MDAKits/issues) what your
thoughts are:

* As a **user**: What do you like or dislike about the MDAKits
  approach? Would you want to use an MDAKit?
* As a **developer**: Would you be interested in creating an MDAKit?
  What should we do to make it easy for you?
  
Get in touch! MDAKits are new and we look forward to adapting the
initial (v0.1.0!) approach based on what we hear from the community.

--- @IAlibay @jbarnoud @orbeckst @richardjgowers  @fiona-naughton @lilyminium 


[10.6084/m9.figshare.20520726.v1]: https://doi.org/10.6084/m9.figshare.20520726.v1
