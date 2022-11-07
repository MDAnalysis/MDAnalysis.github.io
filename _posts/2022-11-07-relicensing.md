---
layout: post
title: Relicensing MDAnalysis
---

<a href="https://www.gnu.org/licenses/lgpl-3.0.en.html">
<img src="https://www.gnu.org/graphics/lgplv3-with-text-154x68.png"
title="LGPLv3" alt="LGPLv3 logo" style="float: right"/>
</a>

This blog post outlines MDAnalysis' proposal to change its license
to the [Lesser GNU Public License][LGPL] (LGPL v3+).

A summary of our [reasons for proposing this license
change](#Rationale-for-changing-licenses), alongside upcoming actions
for [community members and library
contributors](#How-will-the-relicensing-process-work?) are provided.


> ⚠️ **Disclaimer**
> The MDAnalysis core team members are not
> lawyers. As such the information provided here __does not, and is not
> intended to, constitute legal advice__. This blog post also does not
> represent MDAnalysis' full legal position on software licensing; it
> simply aims to inform MDAnalysis developers and users on why
> we believe the library should be relicensed.
>
> Further information on open-source software licensing can be found
> from sources such as the [Open Source Initiative][OSI],
> [tl;drLegal][tldr legal] and the [Software Sustainability Insitute][SSI licensing].
>
> Should you have any concerns about licensing, we always strongly
> recommend getting legal advice before making any decisions on how
> licensing changes may affect you.


## Overview

We want to **change the license of MDAnalysis** from the [GNU General
Public License v2 (or any later versions)][GPLv2] (GPL v2+) to the less
restrictive [Lesser GNU Lesser Public License v3 (or any later versions)][LGPL]
(LGPL v3+) license. Both are [_open source licenses_][OSI] but it is
our view that the [LGPL v3+][LGPL] will give developers more freedom
in how they license any of their own codes that make use of MDAnalysis.

As detailed by the [Open Source Definition][OSD], licenses are core to
the definition of open source. "Open source doesn't just mean access
to the source code". The license defines how code can be used, copied,
changed, and incorporated into other code.

**License changes will affect how people interact with the MDAnalysis code
base going forward. We need the agreement of our contributors and
community members to change from GPL v2+ to LGPL v3+.**

In this post we want to share our motivation, outline the relicensing
process, and invite comments / questions from the community.


## Rationale for license change


### Why is GPL v2+ no longer the best choice?

Since its initial release in 2008, MDAnalysis has grown from a small
Python package used by a handful of enthusiastic graduate students and
postdocs to a mature library that is used by thousands of researchers
in the molecular sciences. The MDAnalysis library was published under
an open source license from the start so that anyone could freely use
it, contribute to it, and build on it. We chose the GNU Public License
version 2+ for this purpose. Thanks to the "copy-left" portion of the
GPL v2+, which requires anyone using MDAnalysis in their own code to
also adopt a compatible version of the GPL for their code, contributors
could feel that any time and work that they donated and invested into
MDAnalysis would not end up contributing to software under non
open-source licensing.

However, the GPL v2+ has also created barriers to adoption of MDAnalysis.
It has, under our (and many other's) interpretation of the GPL v2+
license, prevented developers from making their own code available under
non-GPL licenses. Ultimately, it is the MDAnalysis core team's view that
we do not want to dictate how our developers and users should license
their code, but we also wish to ensure that work on the MDAnalysis
library remains open and free.

Changing to a less restrictive license would benefit the MDAnalysis community, increasing
the number of codes which can use MDAnalysis, and enabling users in
corporate environments to use the library with more certainty. The
reduced licensing complexity also paves the way for our proposed
[MDAKit ecosystem]({{ site.baseurl }}{% post_url 2022-08-24-mdakits-intro %}).


### Why now the LGPL v3+?

We therefore propose to undergo the process of **relicensing MDAnalysis
under the [Lesser GNU Public License v3 (or any later versions)][LGPL]**.
This open source license fulfills a number of important requirements for us:

1. Downstream codes are able to freely import or link to MDAnalysis
   library components without impacting the license choice of the
   downstream code.
   
2. Downstream codes are able to use and subclass any MDAnalysis components
   under its application programming interface (classes, methods, and
   data objects), without impacting the license choice of the
   downstream code.
   
3. Codes that either copy or extend the MDAnalysis library should
   fall under the copyleft license requirements of the MDAnalysis
   library license.
   
Thus, it is our view that the LGPL v3+ license gives people the freedom
to choose any license for their own code that *makes use* of the MDAnalysis
library as a whole (namely ``import MDAnalysis`` or subclassing). This
includes closed / commercial licenses (although we encourage the use of
open source licenses). However, one would not be able to just take parts of
the MDAnalysis code and add it into another codebase unless this
other code is then *also* licensed under a compatible copyleft license
(e.g. GPLv3+/LGPLv3+).

We considered other popular licenses but none fulfilled the requirements listed above.

   
## How will the relicensing process work?

As of writing, MDAnalysis has [over 160 contributors][contributors],
all of whom have contributed code under the terms of the GPL v2+
license. We also have a large user community that uses the library
for many wonderful scientific applications, including several
downstream libraries.

Ultimately, the final decision on relicensing [rests with code
authors](#Contacting-contributors).  However, we fully recognise that
this is a big change for the MDAnalysis user base and the wider
molecular sciences community. As always, we are fully invested in
ensuring that our actions reflect the needs of our community. We
therefore want to give everyone an opportunity to [ask questions about
or comment on the relicensing effort](#Consultation-period) as part of
this process.


### Consultation period (7th November until 5th December 2022)

We will start the process with an open consultation period lasting 28 days from 7th November to 5th December 2022 (anywhere on earth).

During this period we encourage members of the community,
both developers and users, to comment on and ask questions about the
proposed relicensing efforts. The aim is to ensure that relicensing is
indeed in the interest of the community. We will do our best to account
for any concerns raised before attempting to continue with the long and
time-consuming process of relicensing.

We wish to open this conversation on our [public forums]({% link
index.md %}#participating) (mailing lists, discord, twitter). As legal
matters such as licensing can sometimes be sensitive in nature we have
also set up an email address (*licensing@mdanalysis.org*) monitored
solely by the MDAnalysis Core Developers for any private queries that
you may have.

A summary of open discussions and frequently asked questions will be
made available on the [MDAnalysis wiki][faq wiki].

_Note: Whilst the consultation will only last 28 days, we will continue
to engage with conversations on this topic for the entire length of the
relicensing process._


### Contacting contributors (6th December onwards)

After the consultation period, we will contact every code contributor to
the core MDAnalysis library with a request to agree to changing their
contributitions' license from the current "GPL v2 or any later version"
to "LGPL v3 or any later version".

It is important that we hear back from as many contributors as possible.
If you have contributed to MDAnalysis in the past but have since changed
your git-linked contact details, we would kindly ask if you could email
*licensing@mdanalysis.org* to let us know how best to contact you.


### License change

We do not know how long relicensing will take, especially as contacting
historical contributors will likely be a very slow process. Nevertheless,
our aim is to change the license as quickly as possible. We will keep the
community regularly updated on our progress.


## Acknowledgments

We are very grateful for the administrative and legal support from our
fiscal sponsor, [NumFOCUS][]. 


-- The MDAnalysis [Core Developers][]

[OSI]: https://opensource.org/osd
[tldr legal]: https://tldrlegal.com/
[SSI licensing]: https://www.software.ac.uk/resources/guides/choosing-open-source-licence
[OSD]: https://opensource.org/osd
[GPLv2]: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
[LGPL]: https://www.gnu.org/licenses/lgpl-3.0.en.html
[contributors]: https://github.com/MDAnalysis/mdanalysis/blob/develop/package/AUTHORS
[faq wiki]: https://github.com/MDAnalysis/mdanalysis/wiki/GPLv2--to-LGPLv3--relicensing-summary-and-FAQ
[NumFOCUS]: https://www.numfocus.org
[Core Developers]: {% link about.md %}#mdanalysis-core-developers
