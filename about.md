---
layout: page
title: About MDAnalysis
---

<img src="{{ site.baseurl }}/public/mdanalysis-logo_square.png"
style="float: right" alt="MDAnalysis" width="30%"/>


## Community

**MDAnalysis** is developed and maintained as a freely available, open-source
project by a global community of scientists. The MDAnalysis community adheres
to our [Code of Conduct]({{site.baseurl}}/pages/conduct/) and invites everyone
to [participate]({{site.baseurl}}/#participating) --- be it on GitHub Discussions,
through issue reports, or code contributions.

All *contributors* to the MDAnalysis library and its subprojects are acknowledged
in a file called `AUTHORS` in each source code repository and in the list of
contributions; as examples see the [`AUTHORS` file for
mdanalysis](https://github.com/MDAnalysis/mdanalysis/blob/develop/package/AUTHORS)
and the [contributors for
MDAnalysis/mdanalysis](https://github.com/MDAnalysis/mdanalysis/graphs/contributors).

MDAnalysis and the algorithms implemented in the library and the subprojects are
scientific software that are described in [academic
publications]({{site.baseurl}}/pages/citations/).

MDAnalysis is used in a variety of [other molecular dynamics tools]({{site.baseurl}}/pages/mdakits/).


## Open source

All our [code]({{site.baseurl}}/#availability) and our [teaching
materials]({{site.baseurl}}/pages/learning_MDAnalysis/) are available
under open source licenses from repositories at
[https://github.com/MDAnalysis/](https://github.com/MDAnalysis/). The
MDAnalysis library itself is published under the [GNU General Public
License, version 2](https://www.gnu.org/licenses/gpl-2.0.html); other
supporting libraries are published under the MIT or the BSD-3 clause
licence. 

Installable packages are available through the popular ``pip`` and
``conda`` package managers as well as some Linux distributions.


## Governance

Project leadership is provided by a subset of contributors, the
*MDAnalysis Core Developers*
([@MDAnalysis/coredevs](https://github.com/orgs/MDAnalysis/teams/coredevs))
who have produced substantial contributions over extended lengths of
time and who remain active in reviewing issues and discussions on
GitHub Discussions and our Discord server. 

### MDAnalysis Core Developers

The **Core Developers** lead the MDAnalysis project and are responsible to
the community and to NumFOCUS, our fiscal sponsor. They **represent
the project publicly** and **vote to make decisions for the project**.

PIs on a grant submitted by MDAnalysis via NumFOCUS must be Core Developers
while co-PIs do not have to be Core Developers.

Core Developers are granted commit rights (write access) to the [GitHub source
code repositories][orgrepo] and thus can approve pull requests for merges.

The current
[@MDAnalysis/coredevs](https://github.com/orgs/MDAnalysis/teams/coredevs/members)
team ("MDAnalysis Core Developers") consists of:


- @fiona-naughton
- @hmacdope
- @lilyminium
- @micaela-matta
- @orbeckst
- @richardjgowers
- @RMeli
- @tylerjereddy

### MDAnalysis Emeriti Core Developers

**Emerita/Emeritus Core Developers** are former Core Developers who remain
connected to the project but have stepped back from the day-to-day
decision making. Emeriti Core Developers can reinstate themselves to
Core Developer status.

Emeriti Core Developers maintain commit rights (write access) to the
[GitHub source code repositories][orgrepo] and can approve pull requests for
merges.

The current *Emeriti Core Developers* are:

- @dotsdl
- Elizabeth Denning
- @IAlibay
- @jandom
- @jbarnoud
- @kain88-de
- @mnmelo
- @mtiberti
- @nmichaud
- @PicoCentauri
- @seb-buch
- @zemanj

### Decision Making Process and Membership

1. All decisions are made by *simple majority*[^1] of the [MDAnalysis Core
   Developers](#mdanalysis-core-developers).
2. New *Core Developers* are elected with a simple majority of current
   MDAnalysis Core Developers.
3. Current Core Developers are polled annually to *opt-in* to remain
   a Core Developer; otherwise they transition to [Emerita/Emeritus
   Core Developer](#mdanalysis-emeriti-core-developers) status.

[^1]: A [simple majority][] is defined as *more than half the votes
     cast*. Abstentions or blanks are excluded in calculating a
     majority vote. Totals do not include votes cast by someone not
     entitled to vote[^2] or improper multiple votes by a single
     member.  Illegal votes[^3] are counted as votes cast;  
     if only two choices (such as a binary "yes"/"no" vote) are
     possible, a majority vote is more "yes" than "no" votes.
	 	
[^2]: See [MDAnalysis Core Developers](#mdanalysis-core-developers)
    for the list of *individuals entitled to vote*.
	
[^3]: *Illegal votes* are votes that were cast for ineligible choices.

### Code of Conduct

Everyone in the MDAnalysis community adheres to our [Code of
Conduct]({{site.baseurl}}/pages/conduct/).  A rotating subset of three
MDAnalysis Core Developers is tasked to respond to and to investigate
[Code of Conduct]({{site.baseurl}}/pages/conduct/) violations.


## Partners

MDAnalysis is a [fiscally sponsored
project]({{site.numfocus.sponsored_project}}) of [NumFOCUS][], a nonprofit
dedicated to supporting the open source scientific computing
community. 

If you like MDAnalysis and want to support our mission, please
consider making a [donation]({{site.numfocus.donate}}) to support our
efforts. NumFOCUS is a 501(c)(3) non-profit charity in the United
States; as such, donations to NumFOCUS are tax-deductible as allowed
by law.  As with any donation, you should consult with your personal
tax adviser or the IRS about your particular tax situation.

{{ site.numfocus.donate_button }}


## Funding

We are grateful for financial support from the following organizations, which have supported MDAnalysis either through direct funding or indirectly by funding MDAnalysis contributors.

### [Chan Zuckerberg Initiative][] (CZI)
MDAnalysis has been supported by the [Essential Open Source for Science](https://chanzuckerberg.com/rfa/essential-open-source-software-for-science/) (EOSS) program from the CZI Donor-Advised Fund (DAF), an advised fund of Silicon Valley Community Foundation (funder DOI 10.13039/100014989)

- EOSS5, 2022-253062 (**2022**): [MDAnalysis: Outreach and Project Manager](https://chanzuckerberg.com/eoss/proposals/mdanalysis-outreach-and-project-manager/) (**Personnel**: @IAlibay, @jennaswa, @micaela-matta, @orbeckst, @richardjgowers (*PI*))
- EOSS4, DAF2021-237663, DOI [https://doi.org/10.37921/426590wiobus](https://doi.org/10.37921/426590wiobus) (**2021**): [MDAnalysis: Faster, Extensible Molecular Analysis for Reproducible Science](https://chanzuckerberg.com/eoss/proposals/mdanalysis-faster-extensible-molecular-analysis-for-reproducible-science/) (**Personnel**: @fiona-naughton, @hmacdope, @IAlibay, @ianmkenney, @lilyminium, @orbeckst (*PI*), @richardjgowers)

<a href="https://chanzuckerberg.com/"><img
	src="{{site.images}}/CZI_Logo.jpg" title="Chan Zuckerberg
	Initiative" alt="Chan Zuckerberg Initiative" style="display:
	inline; float: left; height: 4em; margin: 0 0.5em" /></a>

### [Google](https://opensource.google/)
The following contributors were sponsored to work on MDAnalysis through the [Google Summer of Code](https://summerofcode.withgoogle.com/) program.

- **2023**: @marinegor, @xhgchen 
- **2022**: @aya9aladdin, @BFedder
- **2021**: @ojeda-e, @orioncohen
- **2020**: @cbouy, @hmacdope, @yuxuanzhuang
- **2019**: @NinadBhat
- **2018**: @ayushsuhane, @davidercruz
- **2017**: @utkbansal
- **2016**: @fiona-naughton, @jdetle

The following technical writers were sponsored to work on MDAnalysis through the [Google Season of Docs](https://developers.google.com/season-of-docs) program.

- **2019**: @lilyminium

<a href="https://summerofcode.withgoogle.com/"> <img
    src="https://developers.google.com/open-source/gsoc/images/gsoc2016-sun-373x373.png"
    title="Google Summer of Code" alt="Google Summer of Code"
    style="display: inline; float: left; height: 4em; margin: 0
    0.5em" /></a>

<a href="https://developers.google.com/season-of-docs"> <img
    src="https://developers.google.com/season-of-docs/images/SeasonofDocs_Icon_Grey_300ppi_trimmed.png"
    title="Google Season of Docs" alt="Google Season of Docs"
    style="display: inline; float: left; height: 4em; margin: 0 0.5em" /></a>

### [National Science Foundation](https://www.nsf.gov/)
Earlier work was partially supported by the NSF (as part of award ACI-1443054).

- NSF DIBBS award, ACI-1443054 (**2014**): [CIF21 DIBBs: Middleware and High Performance Analytics Libraries for Scalable Data Science](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1443054) (**MDAnalysis Personnel**: @orbeckst)

The following students were sponsored to work on MDAnalysis through the [NSF Research Experience for Undergraduates](https://www.nsf.gov/crssprgm/reu/) (REU) program.

- **2021**: @ALescoulie, @edisj
- **2020**: @edisj
- **2019**: @nawtrey
- **2018**: @hfmull
- **2017**: @kaceyaurum
- **2016**: @rbrtdlgd
- **2015**: @ianmkenney

<a href="https://nsf.gov/">
<img src="{{site.images}}/nsf.jpg" title="National Science
	Foundation" alt="National Science Foundation" style="display:
	inline; float: left; height: 4em; margin: 0 0.5em" /></a>

### [NumFOCUS][]
MDAnalysis thanks NumFOCUS for its continued support as our fiscal sponsor, as well as through its [Small Development Grants](https://numfocus.org/programs/small-development-grants) (SDG) program.

- SDG Round 2 (**2023**): Unified and comprehensive documentation and learning resources for MDAnalysis (**Personnel**: @IAlibay, @jennaswa, @lilyminium, @micaela-matta, @orbeckst)
- SDG Round 2 (**2022**): Improving the organization and content of MDAnalysis teaching materials (**Personnel**: @micaela-matta, @pgbarletta)
- SDG Round 1 (**2020**): Periodic boundary handling and on the fly transformations
- SDG Round 2 (**2018**): MDAnalysis tutorial and hackathon
- SDG Round 1 (**2017**): Widening platform availability for MDAnalysis: Full Python 3 Support

<a href="{{site.numfocus.sponsored_project}}"><img
    src="{{site.images}}/numfocus-sponsored.png" title="NumFOCUS
    sponsored project" alt="NumFOCUS Sponsored" style="display:
    inline; float: left; height: 4em; margin: 0 0.5em" /></a>

### [Outreachy](https://www.outreachy.org/)
The following contributors were sponsored to work on MDAnalysis through the [Outreachy](https://www.outreachy.org/) program.

- **2022**: @umak1106

<a href="https://www.outreachy.org/"><img
    src="{{ site.baseurl }}{{ site.images }}/Outreachy-logo.svg"
    title="Outreachy" alt="Outreachy"
    style="display: inline; float: left; height: 4em; margin: 0 0.5em" /></a>

### [Station1](https://www.station1.org/)
The following contributors were sponsored to work on MDAnalysis through the [Station1 Frontiers Fellowship](https://www.station1.org/sff) program.

- **2023**: @jong9559, @KarenBekhazi

<a href="https://www.station1.org/"><img
    src="/public/images/station1_condensed_logo.png"
    title="Station1 Logo" alt="Station1 Logo"
    style="display: inline; float: left; height: 4em; margin: 0 0.5em " /></a>

<div style="clear: both"></div>

## Feedback

MDAnalysis welcomes feedback from its users and community to help improve. If you have any general feedback or comments to make about the MDAnalysis, the community, events, or other aspects, please [let us know at this form here](https://forms.gle/n8GLe2QsL2hW2QiDA)!

------


[NumFOCUS]: https://www.numfocus.org
[simple majority]: https://en.wikipedia.org/wiki/Majority#Majority_vote
[orgrepo]: https://github.com/MDAnalysis
[Chan Zuckerberg Initiative]: https://chanzuckerberg.com/
