---
layout: post
title: Google Summer of Code and Outreachy Students 2024
---

We are happy to announce that MDAnalysis is hosting three [Google Summer of Code][gsoc] (GSoC) contributors - [@ljwoods2](https://github.com/ljwoods2), [@talagayev](https://github.com/talagayev) and [@lunamorrow](https://github.com/lunamorrow) - and an [Outreachy][outreachy] intern - [@adetutudeborah](https://github.com/adetutudeborah). MDAnalysis has been accepted as its own [organization with GSoC][mda-gsoc] for the fifth year running, and we are excited to be hosting our second-ever [MDAnalysis Outreachy internship](https://www.outreachy.org/alums/2024-05/). We are grateful to Google and Outreachy for granting us the opportunity to get started on _four_ very exciting projects!

# GSoC

<a href="https://summerofcode.withgoogle.com/"> <img
    src="https://developers.google.com/open-source/gsoc/images/gsoc2016-sun-373x373.png"
    title="Google Summer of Code" alt="Google Summer of Code"
    style="display: inline; float: right; height: 4em; margin: 0
    0.5em" /></a>

## Lawson Woods: [Enhancing the Interoperability and Efficiency of the Zarrtraj MDAKit](https://summerofcode.withgoogle.com/programs/2024/projects/BYYAE9MR)

<img
src="https://avatars.githubusercontent.com/ljwoods2"
title="Lawson Woods" alt="Lawson Woods"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Streaming trajectory data directly from cloud storage will allow researchers to share their trajectories more easily than ever before, however, rewriting existing trajectory file readers for cloud-streaming would be a slow and complex process. [zarrtraj](https://github.com/Becksteinlab/zarrtraj), a trajectory file format based on [Zarr](https://zarr.readthedocs.io/en/stable/) storage, presents an alternative which natively supports AWS S3, Google Cloud Buckets, and Azure Blob Storage in addition to reading from a local directory without much added complexity. Lawson's GSoC 2024 project will be to continue enhancing the ease-of-use, robustness, and efficiency of this format.

Lawson is a computer science student at Arizona State University (ASU) graduating in December 2024. You can find Lawson on GitHub as [@ljwoods2](https://github.com/ljwoods2) and on [LinkedIn](https://www.linkedin.com/in/lawson-woods/). To see updates on the `zarrtraj` project as well as his other projects you can check out [his blog](https://ljwoods2.github.io/).

## Valerij Talagayev: [2D Visualization for Small Molecules](https://summerofcode.withgoogle.com/programs/2024/projects/sfy3kuqc)

<img
src="https://avatars.githubusercontent.com/talagayev"
title="Valerij Talagayev" alt="Valerij Talagayev"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

2D Visualization of small molecules is quite a helpful way to recognize the different parts of the small molecule. The main objective of Valerij's GSoC project will revolve around the 2D visualization of small molecules through the application of the MDAnalysis RDKit converter, which would allow to convert ligands into RDKit objects and thus allow the user to look at the ligand in his MDAnalysis Universe in an easy to use interactive display in [Jupyter Notebook](https://jupyter-notebook.readthedocs.io/en/stable/).

Valerij is a pharmacy student who graduated from the University of Marburg, where he worked at the [Kolb Lab](https://www.uni-marburg.de/en/fb16/ipc/kolb-group). Valerij is currently a PhD student at the Free University of Berlin in the [Wolber Lab](https://www.bcp.fu-berlin.de/en/pharmazie/faecher/pharmazeutische_chemie/wolber/index.html) with his research revolving around Toll-like receptor antagonists with addition to method development such as [OpenMMDL](https://github.com/wolberlab/OpenMMDL).

You can find Valerij on Github as [talagayev](https://github.com/talagayev) and on [LinkedIn](https://www.linkedin.com/in/valerij-talagayev-260bb820b). To see the updates on the 2D visualization of small molecules project follow his [blog](https://talagayev.github.io/). 

## Luna Morrow: [Extend MDAnalysis Interoperability with OpenBabel](https://summerofcode.withgoogle.com/programs/2024/projects/yLzX6MjS)

<img
src="Picture Here"
title="Luna Morrow" alt="Luna Morrow"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Write a short description of your project here!

Write something about yourself here! (affiliation, background, interests - work and/or personal, etc.)

Add where we can find you here! (github, linkedin, twitter, others, etc.)

Add a link to your blog here!

# Outreachy

<a href="https://www.outreachy.org/"><img
    src="{{ site.baseurl }}{{ site.images }}/Outreachy-logo.svg"
    title="Outreachy" alt="Outreachy"
    style="display: inline; float: right; height: 4em; margin: 0 0.5em" /></a>

## Adetutu Oluwasanmi: [Develop a Communications Strategy for a Growing MDAnalysis User and Contributor Base](https://www.outreachy.org/alums/2024-05/)

<img
src="https://avatars.githubusercontent.com/u/69110554?s=400&u=e20bfb5b20f86b27359b5443084c96016a9817ac&v=4"
title="Adetutu Oluwasanmi" alt="Adetutu Oluwasanmi"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />


MDAnalysis aims to establish clear communication channels to effectively disseminate information, engage the community, and encourage participation. The main objective of Adetutu's [Outreachy](http://outreachy.org/) project is to develop a social media and communications strategy for MDAnalysis to grow its user and contributor base.

Adetutu will conduct surveys to understand the communication preferences of MDAnalysis developers and users. Based on the survey data, she will recommend the best platforms for discussions and increasing engagement.
Additionally, she will set up an email newsletter to keep the community updated on MDAnalysis news, announcements, and events.

Adetutu Oluwasanmi is a First-Class Microbiology graduate from Lagos State University. She recently completed a 1-year diploma in software engineering at AltSchool Africa, where she gained coding and software development skills.
Her previous roles include social media manager, communications intern, and freelance science writer, giving her a diverse skill set that she brings to her current project.

Adetutu is on GitHub as [@adetutudeborah](https://github.com/adetutudeborah), and on Twitter as [@adetutuoluwa2](https://twitter.com/adetutuoluwa2).

She will be documenting her internship journey on her [Medium blog](https://medium.com/@adetutuoluwasanmi). Feel free to check it out for updates on her experience!


— @BradyAJohnston @cbouy @hmacdope @IAlibay @jennaswa @micaela-matta @orbeckst @richardgowers @xhgchen @yuxuanzhuang (mentors and org admins)

[gsoc]: https://summerofcode.withgoogle.com
[outreachy]: https://www.outreachy.org/
[mda-gsoc]: https://summerofcode.withgoogle.com/programs/2024/organizations/mdanalysis
