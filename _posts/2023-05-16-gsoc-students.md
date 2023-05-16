---
layout: post
title: Google Summer of Code Students 2023
---

We are happy to announce that MDAnalysis is hosting two [GSoC][gsoc] students this year - [@xhgchen](https://github.com/xhgchen) and [@marinegor](https://github.com/marinegor). MDAnalysis has been accepted as its own [organization with GSoC][mda-gsoc] for the fourth year running, and we are grateful to Google for granting us the opportunity to get started on _two_ very exciting GSoC projects!

## Xu Hong Chen: Add calculations for self-diffusivity and viscosity

<img
src="https://avatars.githubusercontent.com/u/110699064?v=4"
title="Xu Hong Chen" alt="Xu Hong Chen"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Transport properties are values defining mass, momentum, heat, or charge transfer and are vital to biomolecular research and chemical engineering. MDAnalysis currently provides little support for these calculations and adding them to MDAnalysis would greatly benefit researchers and engineers. Xu Hong will implement a class to calculate self-diffusivity via the Green-Kubo method and utility functions to compute shear viscosity via the Einstein and Green-Kubo methods. He will test and revise his algorithms rigorously to ensure they are accurate and efficient. To maximize user accessibility, Xu Hong will write clear and understandable documentation linked to the best practices in the literature. Xu Hong's project page can be found [here](https://summerofcode.withgoogle.com/programs/2023/projects/4vt9npUg).

Xu Hong is an incoming computer science student at the University of British Columbia and a biochemistry graduate from the University of Alberta. He has 2 years of experience in molecular dynamics from his research on the molecular mechanism of Hsp90 in the [Spyracopoulos Lab](https://lspy.biochem.ualberta.ca/index.php). He is interested in scientific software, high-performance computing, computational research, and open-source software. To rest and relax, Xu Hong enjoys playing piano and listening to Chopin, Schubert, and Bach.

Xu Hong is on GitHub as [@xhgchen](https://github.com/xhgchen) and on LinkedIn as [Xu Hong Chen](https://www.linkedin.com/in/xu-hong-chen/). He will be sharing his experience on [his blog](https://xhgchen.github.io/).

## Egor Marin: Implementation of parallel analysis in MDAnalysis

<img
src="https://pbs.twimg.com/profile_images/1261341752281182208/K0-9mHNm_400x400.jpg"
title="Egor Marin" alt="Egor Marin"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Even though MDAnalysis (since 2.0) allows for serialization of almost all fundamental components, it still lacks a seamless way to run analysis of trajectories in a parallel fashion. For his [GSoC project](https://summerofcode.withgoogle.com/programs/2023/projects/cOTjpLid), Egor will be implementing a parallel-ready backend for MDAnalysis, hoping that it will increase the speed and accessibility of large-scale molecular dynamics analysis.

Egor will hopefully graduate this year fom the PhD program of the University of Groningen (Netherlands). His main expertise lies in structural biology, as can be seen by his [work](https://scholar.google.com/citations?user=FJbv9XcAAAAJ) during his bachelor's and master's. Apart from that, he does a lot of computer administration work, and enjoys occasional bouldering and snowboarding.

You can find Egor on [github](https://github.com/marinegor) and [twitter](https://twitter.com/egor__marin). You can follow Egor's GSoC work throughout the summer on [his blog](https://marinegor.github.io/).

— @hmacdope @orionarcher @yuxuanzhuang @RMeli @orbeckst @richardjgowers @ianmkenney @IAlibay @jennaswa (mentors and org admins)

[gsoc]: https://summerofcode.withgoogle.com
[mda-gsoc]: https://summerofcode.withgoogle.com/programs/2023/organizations/mdanalysis
