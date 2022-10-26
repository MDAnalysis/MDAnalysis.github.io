---
layout: post
title: Relicensing MDAnalysis
---

<img src="https://www.gnu.org/graphics/lgplv3-with-text-154x68.png"
title="LGPLv3" alt="LGPLv3 logo" style="float: right"/>

We want to **change the license of MDAnalysis** from the [GNU General
Public License][GPLv2] (GPL v2+) to the more permissive [Lesser GNU
Lesser Public License][LGPL] (LGPL v3+) license. Both are *open source
licenses* but the [LGPL][] gives developers more freedom in how they
license any of their own code that makes use of MDAnalysis.

The license is the legal basis that determines how code can be copied,
changed, and incorporated into other code. The license is the only
reason that one can use someone else's code --- it does not matter if
the code is publicly available somewhere. If it does not have a
license you cannot use it legally. Therefore, getting the license
right is very important.

**In order to change the MDAnalysis license, we need the agreement of
all our contributors to change from GPL to LGPL.**

In this post we want to share the motivation, outline the relicensing
process, and invite questions from the community.



## Why not GPL any more?

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


## Why now the LGPL?

We therefore decided to start the process to **relicense MDAnalysis
under the Lesser GNU Public License v3 (LGPL)**. This open source
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
(namely ``import MDAnalysis`` or subclassing). They can even make a
commercial product out of it with a restrictive license (although we
encourage open source licenses). However, one cannot just take parts
of the MDAnalysis code and put it into another piece of code unless
the other code is then *also* licensed under the GPL/LGPL. The LGPL
ensures that the MDAnalysis library itself is guaranteed to remain
open and free.

We considered other licenses but none fulfilled our requirements.

   
## How will the relicensing process work?

MDAnalysis has over 150 contributors. They all contributed code under
the GPL v2+. In order to relicense, we need *all our
contributors to agree to the license change*. 

We know that this is a big change for the MDAnalysis community (and
it's a big administrative and legal undertaking for us). We therefore
want to give everyone an opportunity to ask questions so that you can
make an informed decision.

### Consultation period

For the next month, the MDAnalysis Core Developers will answer your
questions regarding the relicensing in our [public forums]({% link
index.md %}#participating) (mailing lists, discord, twitter).

### Contacting contributors

After the consultation period, we will contact every contributor to
MDAnalysis with a request to agree to changing the license from the
current "GPL v2 or any later version" to "LGPL v3 or any later
version".

### License change

Once we have the necessary agreement, we will switch over the license
of MDAnalysis.

### Timeline

* November 2022: consultation period
* December 2022: start contacting contributors
* 2023?: We do not know how long the process will take. Our aim is to
change the license as quickly as possible. We will keep the community
updated.


## Acknowledgments

We are very grateful for the administrative and legal support from our
fiscal sponsor, [NumFOCUS][]. 


-- The MDAnalysis [Core Developers][]

[GPLv2]: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
[LGPL]: https://www.gnu.org/licenses/lgpl-3.0.en.html
[NumFOCUS]: https://www.numfocus.org
[Core Developers]: {% link about.md %}#mdanalysis-core-developers
