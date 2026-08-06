---
layout: post
title: MDAnalysis Awarded Research Software Maintenance Fund Grant
---

![Software Sustainability Institute Research Software Maintenance Fund Round 2 logo](https://www.software.ac.uk/sites/default/files/2026-04/SSI%20-%202025-12-12T100432.105.png)

We are excited to announce that MDAnalysis has been awarded a grant from the [Software Sustainability Institute's][SSI] [Research Software Maintenance Fund][RSMF] (RSMF), a £5.8M investment from [UKRI][] to support the long-term sustainability of research software that matters to the UK research community. MDAnalysis is one of [19 projects funded in Round 2][RSMF round 2 announcement], out of 229 initial expressions of interest.

This funding will allow us to hire a dedicated Research Software Engineer (RSE) to carry out essential maintenance work which is key for long-term sustainability but difficult to achieve through volunteer effort alone. It will also contribute to retaining our Community Manager to continue the community engagement activities and mentorship programs that are central to the project.

## Why this grant matters

Despite MDAnalysis's scale of impact the project has, until now, relied entirely on volunteer effort to keep the library maintained. This is the **first grant MDAnalysis has ever received that is dedicated specifically to maintenance** rather than new feature development, and it comes at an important time: previous funding streams that supported MDAnalysis, such as the Chan Zuckerberg Initiative's Essential Open Source Software for Science program (which funded our [EOSS4][] and [EOSS5][] awards), have been discontinued, narrowing the available opportunities.

## What the funding will support

Over the 12-month grant period, our RSE and Community Manager will work across three areas:

**1. Sustainable maintenance.** The RSE will audit and document our continuous integration and release workflows, lead regular releases, and keep MDAnalysis compatible with fast-moving upstream dependencies such as NumPy, Python, and Cython. Other planned work includes modernising our test suite, maintaining the [MDAKits][] ecosystem, and adding support for new file formats, alongside new "Guesser" components for PDB structures and MARTINI coarse-grained simulations.

**2. Documentation and user experience.** We will overhaul the [User Guide][], cleanly separating conceptual/user-facing documentation from the API reference, consolidate and add continuous testing to our library of tutorials, and mine our [community support channels][] (GitHub Discussions, Discord, mailing lists) to build out a FAQ section.

**3. Community building and mentorship.** Our Community Manager will continue coordinating our participation in mentorship programs such as Google Summer of Code, maintain onboarding pathways for new Core Developers, and keep up our long-standing training partnerships with UK-based communities including the [Thomas Young Centre][], [CCPBioSim][], and [CCP5][]. We will also continue to grow relationships with industry stakeholders to help diversify our funding base beyond this grant.

We will be advertising for the RSE position soon — if you're passionate about open-source scientific software and want to help sustain a tool used by thousands of researchers worldwide, keep an eye on our channels for the job posting.

## Acknowledgements

This grant recognises that maintaining widely used research software requires sustained, professional investment, and we're grateful to the SSI and UKRI for supporting that vision. Thank you also to King's College London and our partners at the Thomas Young Centre, CCPBioSim, and CCP5, and to our fiscal sponsor [NumFOCUS][] for its continued support. Finally, thank you to everyone who uses, contributes to, and advocates for MDAnalysis — this grant is as much a recognition of your work as it is an investment in our future.

This project will be led by Micaela Matta @micaela-matta (Lead) and Jenna Swarthout Goddard @jennaswa (Co-Lead and Community Manager).

## About MDAnalysis

The MDAnalysis package is one of the most widely used Python-based libraries for the analysis and manipulation of molecular simulations. The objective of MDAnalysis is to provide simple, flexible, and efficient means of handling molecular structure data from simulations and experiments and support research in biophysics, biochemistry, materials science and beyond.

[MDAnalysis][] is a fiscally-sponsored project of [NumFOCUS][], a nonprofit dedicated to supporting the open source scientific computing community.


[SSI]: https://www.software.ac.uk/
[RSMF]: https://www.software.ac.uk/programmes/research-software-maintenance-fund
[RSMF round 2 announcement]: https://www.software.ac.uk/news/research-software-maintenance-fund-awards-funding-19-projects-round-2
[UKRI]: https://www.ukri.org/
[EOSS4]: {{ site.baseurl }}{% post_url 2021-08-31-CZI-EOSS4 %}
[EOSS5]: {{ site.baseurl }}{% post_url 2022-11-10-CZI_EOSS5 %}
[community support channels]: {{ site.baseurl }}/community/#ask-questions--get-help
[MDAKits]: https://mdakits.mdanalysis.org/
[User Guide]: https://userguide.mdanalysis.org/
[Thomas Young Centre]: https://www.thomasyoungcentre.org/
[CCPBioSim]: https://ccpbiosim.ac.uk/
[CCP5]: https://www.ccp5.ac.uk/
[MDAnalysis]: {{ site.baseurl }}/about
[NumFOCUS]: https://www.numfocus.org/
