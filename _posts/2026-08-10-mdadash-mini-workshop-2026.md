---
layout: post
title: "MDAnalysis Dashboard - Mini Workshop (Aug 20, 2026)"
---

[Release 2.10.0]({% post_url 2025-10-26-release-2.10.0 %}) of MDAnalysis introduced support for interactive MD (IMDv3) stream reading using the IMDReader and the [imdclient][] package. This allows analyzing MD simulation trajectories live while they are being generated and thus enabling analysis of sub-picosecond dynamics in trajectories.

As part of a Google Summer of Code (GSoC) 2026 [project](https://summerofcode.withgoogle.com/programs/2026/projects/DzTMshtu), we built a new tool called [mdadash](https://github.com/MDAnalysis/mdadash) (short for 'MDAnalysis Dashboard') that provides a browser-based real-time dashboard for monitoring, tracking and analyzing running MD simulations using this new [IMDv3 streaming][imdclient] interface.

The dashboard provides an interactive GUI for defining and running real-time analyses. There is support for per-frame observables and buffered, time-dependent analyses along with live visualization of analysis results as real-time plots. Custom user-defined analyses are also supported through built-in Notebook functionality.

<img src="{{ site.images }}/mdadash-mini-workshop-2026.png" alt="MDAnalysis Dashboard"/>

## Workshop Overview

The program will run from 9:30 am to 11:30 am Pacific time (UTC 16:30 – 18:30) on Thursday, August 20th. The workshop will introduce [mdadash](https://github.com/MDAnalysis/mdadash), demonstrate its various features, use cases and how it can be customized for your own analysis needs. The interactive activities will allow participants to try out the dashboard in an easy-to-use workshop environment.

## Registration

Attendance at this workshop will be *free*, and we encourage anyone with an interest in attending to register below. 

[Registraion link to be added]

## Workshop materials

All materials are made available in the [github.com/PardhavMaradani/mdadash-mini-workshop-2026](https://github.com/PardhavMaradani/mdadash-mini-workshop-2026) repository.

Prepare for the interactive workshop activities by following the [setup instructions](https://github.com/PardhavMaradani/mdadash-mini-workshop-2026?tab=readme-ov-file#setup-instructions).

## Who to Contact

If you have any questions or special requests related to this workshop, you may [contact the organizing committee](mailto:workshops@mdanalysis.org).

- @PardhavMaradani
- @orbeckst (MDAnalysis)
- @amruthesht @HeydenLabASU  (ASU)
- @jeremyleung521 (University of Pittsburgh)


[imdclient]: https://imdclient.readthedocs.io/
