---
layout: post
title: An update on relicensing the MDAnalysis library (towards LGPLv2.1+)
---

> [!WARNING]
> **Disclaimer**
> The MDAnalysis core team members are not
> lawyers. As such the information provided here __does not, and is not
> intended to constitute legal advice__. This blog post also does not
> represent MDAnalysis' full legal position on software licensing; it
> simply aims to inform MDAnalysis developers and users on certain
> aspects of the MDAnalysis relicensing efforts and how we think it
> may impact them going forward.
>
> Further information on open-source software licensing can be found
> from sources such as the [Open Source Initiative][OSI],
> [tl;dr Legal][tldr legal] and the [Software Sustainability Institute][SSI licensing].
>
> Should you have any concerns about licensing, we always strongly
> recommend getting legal advice before making any decisions on how
> licensing changes may affect you.

## Short description

This blog post provides an update on MDAnalysis’ proposed license change. Specifically, we outline our
decision to opt for a [GNU Lesser General Public License v2.1+ (LGPL v2.1+)][LGPL v2.1] instead of
the [previously communicated]() [GNU Lesser General Public License v3+ (LGPL v3+)][LGPL v3].

We also detail some [interim changes](#why-the-interim-gpl-v3-and-lgpl-v3-package-licenses-the-issue-with-apache-20)
to the MDAnalysis [package license](#package-license), alongside some
[upcoming actions for library contributors](#package-license).

## Overview

In November of 2022 the MDAnalysis team [announced a proposal to change the license of its core library](https://www.mdanalysis.org/2022/11/07/relicensing/#relicensing-mdanalysis)
from the [GNU General Public License v2+][GPL v2] to the [GNU Lesser General Public License v3+][LGPL v3].
The full rationale for this change can be found [here](https://www.mdanalysis.org/2022/11/07/relicensing/#rationale-for-license-change)
and mainly centered around three main arguments:

1. Allowing downstream codes to freely import or link to the MDAnalysis library with minimal worry about the license they choose to use.
2. Allowing downstream codes to use and subclass any MDAnalysis components under its application programming interface
   (classes, methods, and data objects), without incurring copyleft and restricting the license choice of the downstream code.
3. Ensuring that codes that either directly copy or extend the MDAnalysis library’s source code would fall under the copyleft
   license requirements of the MDAnalysis library license. Doing so would ensure that the eventual license would more closely match the
   original intent of the historical [GPL v2+][GPL v2] contributions which were kindly donated to MDAnalysis over the years.

It is our view that a weak copyleft license, such as the [GNU Lesser General Public License][LGPL v3] would best match these requirements.

Upon announcing the proposed license change, we engaged in a [consultation period](https://www.mdanalysis.org/2022/11/07/relicensing/#consultation-period-7th-november-until-5th-december-2022)
where we encouraged members of the community to comment on and ask questions about the proposed relicensing efforts.
Overall, the community response was overwhelmingly positive, and led to some very fruitful discussions. We would like
to thank all members of the community who participated. Some interesting points were brought to our attention, particularly
when it came to the range of licenses compatible with [LGPL v3+][LGPL v3]. This raised a lot of questions about the
suitability of [LGPL v3+][LGPL v3], and the current status of the MDAnalysis [GPL v2+][GPL v2] package license, so we
decided to pause the relicensing process whilst we worked with our legal counsel to find an appropriate path forward.

After careful consideration we have decided to aim for a license change to the [LGPL v2.1+][LGPL v2.1] via interim
[package license](#package-license) changes to [GPL v3+][GPL v3] and [LGPL v3+][LGPL v3].

In this post we will provide our rationale for this change, alongside key information for how we will be
enacting these changes.

**License changes will affect how users and developers will interact with the MDAnalysis code base going
forward. We urge all members of the community to read this post and [get in touch](#getting-in-touch) with
the MDAnalysis core developer team should you have any questions or concerns**.

## Rationale for our change in target license

*Table 1: MDAnalysis’ view of the dynamic linking license compatibility matrix. License pairs
detailed with a compatibility of “No” have claimed incompatibilities between their license terms.
License pairs detailed as “Yes ([copyleft](#copyleft))” are compatible but may incur copyleft through the
interacting GNU GPL license. License pairs detailed with the compatibility “Yes” are, to our knowledge,
fully compatible without incurring any [copyleft](#copyleft) penalty.*

| License         | [GPL v2][]          | [GPL v3][]          | [LGPL v2.1][]       | [LGPL v3][]         | [Apache 2.0]()     | [MIT][]            | [BSD 3.0][]        |
|-----------------|---------------------|---------------------|---------------------|---------------------|--------------------|--------------------|--------------------|
| [GPL v2][]      | Yes                 | No \[1\]            | Yes ([copyleft]())  | No \[1\]            | No                 | Yes ([copyleft]()) | Yes ([copyleft]()) |
| [GPL v3][]      | No \[1\]            | Yes                 | No \[2\]            | Yes ([copyleft]())  | Yes ([copyleft]()) | Yes ([copyleft]()) | Yes ([copyleft]()) |
| [LGPL v2.1][]   | Yes ([copyleft]())  | No                  | Yes                 | No \[2\]            | No                 | Yes                | Yes                |
| [LGPL v3][]     | No \[1\]            | Yes ([copyleft]())  | No \[2\]            | Yes                 | Yes                | Yes                | Yes                |
| [Apache 2.0][]  | No                  | Yes ([copyleft]())  | No                  | Yes                 | Yes                | Yes                | Yes                |
| [MIT][]         | Yes ([copyleft]())  | Yes ([copyleft]())  | Yes                 | Yes                 | Yes                | Yes                | Yes                |
| [BSD 3.0][]     | Yes ([copyleft]())  | Yes ([copyleft]())  | Yes                 | Yes                 | Yes                | Yes                | Yes                |

*\[1\] Compatible if the [GPL v2][] license is provided under the “or any later version” terms.
\[2\] Compatible if the [LGPL v2.1][] license is provided under the “or any later version” terms.*

### Why LGPL v2.1+ instead of LGPLv3+?

During the community consultation, it was brought to our attention that by choosing an [LGPL v3+][LGPL v3]
license, it would no longer be compatible with older versions of the GNU General Public License, namely
[GPL v2][] and [LGPL v2.1][] (see Table 1). The reasons for this are complex and we will not attempt to
explain them here, we will instead refer the reader to the [FSF GPL license compatibility FAQ entry](https://www.gnu.org/licenses/gpl-faq.html#AllCompatibility)
for further details.

Overall the impact of this incompatibility is likely very low. Few codes are solely licensed under
[GPL v2][] or [LGPL v2.1][], with most under “or any later version” terms, allowing the codes to relicense
under a later version and avoiding this issue entirely. However, this would still require several downstream
packages to update their [package license](#package-license), which is not an insignificant amount of work.
It is our view that it is best for the community that we try to avoid this issue, especially given that the
MDAnalysis library itself has historically been released under [GPL v2+][GPL v2].

Our solution is instead to compromise and try to ultimately relicense MDAnalysis under the terms of the [LGPL v2.1+][LGPL v2.1].
This would allow other libraries to interact with MDAnalysis under either the terms of [LGPL v2.1][], [LGPL v3][],
or any as of yet unreleased future versions of the GNU Lesser General Public License. This does mean that the
MDAnalysis library would not fully benefit from some of the newer protections and clarifications introduced in
the [LGPL v3][] license, but we do not believe this to have a significant impact on how we expect the library
to be developed and used.

### Why the interim GPL v3+ and LGPL v3+ package licenses? (The issue with Apache 2.0).

Whilst reviewing the compatibility of GNU GPL licenses, we were made aware of a likely incompatibility between
the [Apache 2.0][] license and either [GPL v2][] or [LGPL v2.1][]. The exact details are again quite complex and
murky, and we would instead refer you to [expert opinions][] on the matter. Nevertheless, it is our understanding
that the prevalent opinion is that (L)GPLv2 codes cannot dynamically link to [Apache 2.0][] codes, something which
the MDAnalysis library has historically been doing by depending on the [fasteners][] and [mmtf-python][] packages.

Both of these packages bring important core functionality to MDAnalysis and cannot be quickly removed without impacting
the user experience. We believe that this is something we can do as we move towards a v3.0 release of MDAnalysis,
replacing the use of [fasteners][] with alternatives such as [py-filelock][] and removing the ability to read
[MMTF][] files in favor [newer PDB file formats](https://mmcif.wwpdb.org/). As this will likely take some time,
we will temporarily change our [package license](#package-license) to [GPL v3+][GPL v3], which we have already done
as part of the [2.6.0 release of MDAnalysis](https://www.mdanalysis.org/blog/#release-2-6-0-and-2-6-1-of-mdanalysis).
We may also possibly go through an interim [LGPL v3+][LGPL v3] license, before reaching our final [LGPL v2.1+][LGPL v2.1]
[package license](#package-license). Further details how we see the license process to occur can be found
[later in this blog](#detailing-our-proposed-relicensing-process).

We note that we can make these changes to the MDAnalysis [package license](#package-license) without requiring
sign offs from all developers as the library’s [source code license](#source-code-license) will still remain
under [GPL v2+][GPL v2] (and eventually [LGPL v2.1+][LGPL v2.1]). This allows us to release the package under
any appropriate license version of the GNU GPL. It does however mean that we will be temporarily at risk of
downstream packages encountering the [previously mentioned incompatibilities](#why-lgpl-v21-instead-of-lgplv3)
between [GPL v2][] and [GPL v3][]. **Should this impact you, please [get in touch with the core developer team](#getting-in-touch)
and we will attempt to explore alternative packaging solutions**.

### Why not opt for a more permissive (and likely less complicated) license such as MIT?

As detailed in our [original justification](https://www.mdanalysis.org/2022/11/07/relicensing/#rationale-for-license-change)
we believe that a weak copyleft license, i.e. the Lesser GNU GPL licenses, best respects the original intent of historical
contributions to MDAnalysis. Additionally, we are also aware of various portions of the MDAnalysis library which likely have
derivative code from other [LGPL v2.1+][LGPL v2.1] licensed codes, such as [GROMACS](https://www.gromacs.org/about.html#license).
We do not believe it to be possible at this point to remove these portions of the library without significantly impacting the
user experience.

## Detailing our proposed relicensing process

**_Table 2. A stepwise description of the proposed MDAnalysis relicensing process._**

| Step | [Package License](#package-license) | [Source Code License](#source-code-license) from Historical Contributions | [Source Code License](#source-code-license) from New Contributions | Details & Conditions of Change                                                        | MDAnalysis Versions |
|------|-------------------------------------|---------------------------------------------------------------------------|-------------------------------------------------------------------|---------------------------------------------------------------------------------------|---------------------|
| 1.   | [GPL v2+][GPL v2]                   | [GPL v2+][GPL v2]                                                         | [GPL v2+][GPL v2]                                                 | This is the initial state of the MDAnalysis library license.                          | < v2.6.0            |
| 2.   | [GPL v3+][GPL v3]                   | [GPL v2+][GPL v2]                                                         | [LGPL v2.1+][LGPL v2.1]                                           | As of commit [44733fc][] all new contributions are made under [LGPL v2.1+][LGPL v2.1] | v2.6.0+             |
| 3.   | [LGPL v3+][LGPL v3]                 | [LGPL v2.1+][LGPL v2.1]                                                   | [LGPL v2.1+][LGPL v2.1]                                           | Once historical contributors have agreed to relicense their contributions.            | N/A                 |
| 4.   | [LGPL v2.1+][LGPL v2.1]             | [LGPL v2.1+][LGPL v2.1]                                                   | [LGPL v2.1+][LGPL v2.1]                                           | Once Apache 2.0 dependencies have been removed from the library.                      | N/A                 |

In this section, we will briefly go over the various steps we propose to take (Table 2) as we go about relicensing the
[MDAnalysis core library][].

### 1. A GPL v2+ package and source code (MDAnalysis v2.5 and lower)

This is the original state of the license up until the [2.5.x releases of MDAnalysis](https://www.mdanalysis.org/blog/page2/#releases-2-5-0-and-2-4-x-of-mdanalysis).
The package and source code were released under a [GPL v2+][GPL v2] license.

### 2. A GPL v3+ package, and a mix of historical GPL v2+ and new LGPL v2.1+ source code (MDAnalysis v2.6 onwards)

As of commit [44733fc][], all new contributions to MDAnalysis have been made under the terms of the
[LGPL v2.1+][LGPL v2.1] and the package is released under [GPL v3+][GPL v3] to better reflect our
understanding of the [Apache 2.0 incompatibility issues](#why-the-interim-gpl-v3-and-lgpl-v3-package-licenses-the-issue-with-apache-20).
Historical source code contributions remain under the terms of [GPL v2+][GPL v2].

This will remain the state of MDAnalysis until we obtain approval from our developers to relicense
their code to [LGPL v2.1+][LGPL v2.1]. As we obtain contributor agreements, these historical
source code contributions will slowly turn to [LGPL v2.1+][LGPL v2.1].

### 3. An LGPL v3+ package, and an LGPL v2.1+ source code (Optional)

Should we obtain approval to relicense from our developers before we remove our [Apache 2.0][]
dependencies, we will then start releasing the MDAnalysis package under a [LGPL v3+][LGPL v3] license.

### 4. An LGPL v2.1+ package and source code

This is the state which we aim for the library to eventually be licensed under. Once we both obtain
the agreement of our developers to relicense and remove [Apache 2.0][] dependencies, we will be
able to release the MDAnalysis package under [LGPL v2.1+][LGPL v2.1].

## How does this affect me?

### As a contributor

If you have any historical contributions you will be asked to relicense them to [LGPL v2.1+][LGPL v2.1]
instead of the initially proposed [LGPL v3+][LGPL v3].

### As a developer

Ultimately we believe that the proposed change to [LGPL v2.1+][LGPL v2.1] should offer a better compatibility
(Table 1) with whatever license you may choose to use for your own code. The interim package license changes
may however impact you, if you believe this to be the case [please do get in touch with us](#getting-in-touch).

### As a user

We do not anticipate this to have a significant impact on users. There may be some slight changes to how you
[might access MMTF files which are discussed below](#removing-apache-20-dependencies).

## Next steps

### Contacting contributors

We will shortly be contacting every historical contributor to MDAnalysis and ask them to agree to changing
their contribution’s license from the current “GPL v2 or any later version” to “LGPL v2.1 or any later version”.

Specifically we will be asking all contributors to respond with the following statement:

```
I, <your name>, GitHub handle <your GitHub handle>, am the copyright owner of all the code I submitted
to the MDAnalysis library using my GitHub handle. I agree to relicense all my code contributions to the
MDAnalysis library under the LGPLv2.1+ license from this day going forward.

Typed Signature
```

It is important that we hear back from as many contributors as possible. If you have contributed to
MDAnalysis in the past but have since changed your git-linked contact details, we would kindly ask if
you could email *licensing@mdanalysis.org* to let us know how best to contact you.

### Removing Apache 2.0 dependencies

As we move towards a v3.0 release of MDAnalysis, we will be looking to replace existing [Apache 2.0][]
core dependencies. Of particular impact, this may mean that access to reading [MMTF][] files may no longer
be available by default through the MDAnalysis library, but will instead be offered through an optionally
installable [MDAKit](https://mdakits.mdanalysis.org/about.html).

## Getting in touch

We understand that license changes are a complex subject that can have a large impact on how users and
developers may interact with the MDAnalysis library. Should you have any concerns, questions, or comments,
we urge you to get in touch with the core developer team. You can do so over our
[public forums]({% link index.md %}#participating) or privately by emailing *licensing@mdanalysis.org*.

## Definitions

### package license

The license under which the MDAnalysis library is packaged and released. This may be different from
but still compatible with the [source code license](#source-code-license).

### source code license

The license under which the MDAnalysis source code is released. This only accounts for the code itself
and not any impact in licensing incurred from linking (i.e. `python importing`) to other libraries.

### license compatibility

Whether or not the terms of two licenses allows them to be compatible with each other.

### copyleft

Whether the terms of a given license enforces that other codes (even under other licenses) must be bound
by the same conditions. See the [FSF definition](https://www.gnu.org/licenses/copyleft.en.html) for more information.

## Acknowledgements

We are very grateful for the administrative and legal support from our fiscal sponsor, [NumFOCUS][].

This relicensing effort has been made possible in part by the CZI grant DAF2021-237663,
grant DOI https://doi.org/10.37921/426590wiobus from the Chan Zuckerberg Initiative DAF,
an advised fund of Silicon Valley Community Foundation (funder DOI 10.13039/100014989).

– The MDAnalysis [Core Developers][]

[OSI]: https://opensource.org/osd
[tldr legal]: https://tldrlegal.com/
[SSI licensing]: https://www.software.ac.uk/resources/guides/choosing-open-source-licence
[OSD]: https://opensource.org/osd
[GPL v2]: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
[GPL v3]: https://www.gnu.org/licenses/gpl-3.0.en.html
[LGPL v2.1]: https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html
[LGPL v3]: https://www.gnu.org/licenses/lgpl-3.0.en.html
[MIT]: https://opensource.org/license/mit/
[BSD 3.0]: https://opensource.org/license/bsd-3-clause/
[Apache 2.0]: https://www.apache.org/licenses/LICENSE-2.0
[contributors]: https://github.com/MDAnalysis/mdanalysis/blob/develop/package/AUTHORS
[NumFOCUS]: https://www.numfocus.org
[expert opinions]: https://scholar.google.com/scholar?hl=en&as_sdt=0%2C5&q=heather+meeker+open+source+for+business&btnG=
[fasteners]: https://fasteners.readthedocs.io/en/latest/
[mmtf-python]: https://github.com/rcsb/mmtf-python
[py-filelock]: https://github.com/tox-dev/filelock
[MMTF]: https://mmtf.rcsb.org/
[44733fc]: https://github.com/MDAnalysis/mdanalysis/commit/44733fc214dcfdcc2b7cb3e3705258781bb491bd
[MDAnalysis core library]: https://github.com/MDAnalysis/mdanalysis
[Core Developers]: {% link about.md %}#mdanalysis-core-developers

