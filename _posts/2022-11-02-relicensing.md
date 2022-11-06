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
> lawyers, as such the information provided here __does not, and is not
> intended to, constitute legal advice__. This blog post also does not
> represent MDAnalysis' full legal position on software licensing, it
> purely aims to inform MDAnalysis developers and users on why
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
(LGPL v3+) license. Both are [_open source licenses_][OSI] but it is our view
that the [LGPL][] will give developers more freedom in how they license any of
their own codes that makes use of MDAnalysis.

The license is the legal basis that determines how code can be used, copied,
changed, and incorporated into other code. The license is the only
reason that one can use someone else's code --- it does not matter if
the code is publicly available somewhere. If it does not have a
license you cannot use it legally. Therefore, getting the license
right is very important.

**In order to change the MDAnalysis license, we need the agreement of
all our contributors to change from GPL to LGPL.**

In this post we want to share the motivation, outline the relicensing
process, and invite questions from the community.


## Rationale for license change


### Why not GPL any more?

Since its initial release in 2008, MDAnalysis has grown from a small
Python package used by a handful of enthusiastic graduate students and
postdocs to a mature library that is used by thousands of researchers
in the molecular sciences. The MDAnalysis library was published under
an open source license from the start so that anyone could freely use
it, contribute to it, and build on it. We chose the GNU Public License
(GPL) for this purpose. Thanks to the "copy-left" portion of the GPL,
which requires anyone using MDAnalysis in their own code to also adopt
the GPL for their code, contributors could feel that any time and work
that they donated and invested into MDAnalysis would not end up
contributing to software that they could not use themselves in the
future.

However, the GPL has also created barriers to adoption of MDAnalysis
in that it prevented developers to make their own code available under
other permissive open source licenses. Ultimately, we do not want to
tell our developers and users how to license their code but we also
want to clearly signal that the work on MDAnalysis itself is
guaranteed to remain open and free.


### Why now the LGPL?

We therefore decided to start the process to **relicense MDAnalysis
under the Lesser GNU Public License v3 ([LGPL][])**. This open source
license fulfills a number of important requirements for us:

1. Downstream codes are able to freely import or link to MDAnalysis
   library components without impacting the license choice of the
   downstream code.
   
2. Downstream codes are able to use and subclass any MDAnalysis components
   under its application programming interface (classes, methods, and
   data objects), without impacting the license choice of the
   downstream code.
   
3. Codes which either copy or extend the MDAnalysis library should
   fall under the copyleft license requirements of the MDAnalysis
   library license.
   
Thus, the LGPL gives people the freedom to choose any license for
their own code that *makes use* of the MDAnalysis library as a whole
(namely ``import MDAnalysis`` or subclassing). They can also make a
commercial product out of it with a restrictive license (although we
encourage open source licenses). However, one cannot just take parts
of the MDAnalysis code and put it into another piece of code unless
the other code is then *also* licensed under the GPL/LGPL. The LGPL
ensures that the MDAnalysis library itself remains open and free.

We considered other popular licenses but none fulfilled our
requirements.

   
## How will the relicensing process work?

As of writing, MDAnalysis has [over 160 contributors][contributors]. All of which have
contributed code under the terms of the GPL v2+ license. We also have
a large user community which uses the library for many wonderful scientific
applications, including several downstream libraries.

Ultimately, the final decision on relicensing [rests with code
authors](#Contacting-contributors).  However, we fully recognise that
this is a big change for the MDAnalysis user base and the wider
molecular sciences community. We are, as always, fully invested in
ensuring that our actions reflect the needs of our community. We
therefore want to give everyone an opportunity to [ask questions about
or comment on the relicensing effort](#Consultation-period) as part of
this process.


### Consultation period (Today until XXth December 2022)

We will start the process with a 28 days open consultation period.

During this period we would like to encourage members of the community, both
developers and users, to comment on and ask questions about the proposed
relicensing efforts. The aim is to ensure that relicensing is indeed in the
interest of the community. We will do our best to account for any concerns
raised before attempting to continue with the long and time-consuming
process of relicensing.

We wish to open this converstation on our [public forums]({% link
index.md %}#participating) (mailing lists, discord, twitter). As legal matters,
such as licensing, can sometimes be sensitive in nature we have also set up
an email address (*licensing@mdanalysis.org*) monitored solely by the MDAnalysis
Core Developers for any private queries which you may have.

_Note: Whilst the consultation will only last 28 days, we will continue
to engage with conversations on this topic for the entire length of the
relicensing process._


### Contacting contributors (XXth December onwards)

After the consultation period, we will contact every code contributor to
the core MDAnalysis library with a request to agree to changing their
contributitions' license from the current "GPL v2 or any later version" to 
"LGPL v3 or any later version".

It is important that we hear back from as many contributors as possible. If you have
contributed to MDAnalysis in the past but have since changed your git-linked contact
details, we would kindly ask if you could email *licensing@mdanalysis.org* to let us
know how best to contact you.


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
[GPLv2]: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
[LGPL]: https://www.gnu.org/licenses/lgpl-3.0.en.html
[contributors]: https://github.com/MDAnalysis/mdanalysis/blob/develop/package/AUTHORS
[NumFOCUS]: https://www.numfocus.org
[Core Developers]: {% link about.md %}#mdanalysis-core-developers
