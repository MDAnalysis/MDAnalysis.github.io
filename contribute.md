---
layout: page
title: Contribute
order: 8
---

MDAnalysis is an open-source project that welcomes and encourages contributions from the community. Whether you're a developer, documentation writer, or interested in improving the website, there are many ways to get involved.

You can contribute by **developing new features, fixing bugs, improving documentation, enhancing tutorials, contributing to MDAKits, improving the website, or engaging in discussions**. 

Contributions also extend to our **mentoring and outreach programs**, where participants gain hands-on experience in open-source development and scientific software. You can also contribute by **mentoring new contributors**, helping them navigate MDAnalysis, review their work, and share best practices.

If you're new, you can start by participating in community conversations and mentoring programs. For more information on how you can participate, check out the [Community]({{ site.baseurl }}/community/) and [Events]({{ site.baseurl }}/events/) pages.

The [MDAnalysis Project][] hosts several open-source repositories beyond the main library and MDAKits. You are welcome to contribute to any of them! Each repository typically has a `CONTRIBUTING.md` file to help you get started. If you’re unsure where to begin, you can:

- Check the repository’s issue tracker for open tasks.
- Submit a pull request for improvements.
- Contact maintainers via [GitHub discussions][] or MDAnalysis [Community]({{ site.baseurl }}/community/) forums.

This page provides an overview of how you can contribute to the [Codebase](#contributing-to-the-codebase), [MDAKits](#contributing-to-mdakits), [Documentation](#contributing-to-the-documentation), or [Mentoring & Outreach Programs](#mentoring--outreach-programs).

## Contributing to the Codebase
If you are interested in developing new features, fixing bugs, or improving performance, contributing to the MDAnalysis main codebase is a great way to get involved. 

### Main Codebase Repository
MDAnalysis is developed and maintained in the [MDAnalysis Main Code Repository]({{ site.github.repo }}), which contains:

- The core library for molecular dynamics analysis.
- The [API documentation]({{ site.docs.mdanalysis.url }}/stable/index.html), generated from in-code docstrings.
- [Open issues]({{ site.github.issues }}), where contributions such as bug fixes and feature implementations are tracked.

#### How to Get Started
- Read the [Contributing to the Main Codebase]({{ site.docs.userguide.url }}/stable/contributing_code.html#working-with-mdanalysis-code) guide, which outlines all the necessary steps.
- Browse [GitHub Issues]({{ site.github.issues }}) to find open tasks. Look for "good first issue" if you're a beginner. You can also report bugs or suggest improvements there.
- Follow best practices, including
    - [Writing and Running Tests]({{ site.docs.userguide.url }}/stable/testing.html#testing) to ensure your changes work correctly.
    - [Updating the Documentation]({{ site.docs.userguide.url }}/stable/contributing_code.html#building-code-documentation) when making modifications.

## Contributing to MDAKits
[MDAKits]({{ site.mdakits.registry }}) are **community-driven projects** that extend MDAnalysis functionality. They allow researchers and developers to build tools tailored for specific workflows while staying within the MDAnalysis ecosystem.

### MDAKits Registry
The [MDAKits Registry]({{ site.mdakits.registry }}) lists active projects contributed by the community. These projects cover various applications, such as enhanced trajectory analysis, molecular simulations, and integration with machine learning frameworks.

#### How to Get Started
- Find open issues in [MDAKits GitHub Issues]({{ site.mdakits.issues }}).
- Propose a new MDAKit by following the [MDAKits Guide]({{ site.mdakits.guide }}).
- Contribute to an existing MDAKit by collaborating with maintainers to improve functionality, fix bugs, or add new features.
- Explore existing MDAKits in the [MDAKits Registry]({{ site.mdakits.registry }}) to find projects that match your interests. Each entry includes functionality, installation instructions, and contribution guidelines.

MDAKits follow MDAnalysis coding standards but operate independently. Contributions to these projects are highly encouraged and help expand the MDAnalysis ecosystem.

## Contributing to the Documentation
Good documentation is essential to making MDAnalysis accessible to users and developers. MDAnalysis maintains two types of documentation: the [User Guide]({{ site.docs.userguide.url }}) and the [API Reference]({{ site.docs.mdanalysis.url }}/stable/index.html).

### User Guide
The [User Guide]({{ site.docs.userguide.url }}/stable/index.html) provides tutorials and explanatory content to help users install and learn MDAnalysis. You can contribute by updating unclear explanations, adding new tutorials or examples, or fixing outdated information. 

#### How to Get Started
- Work on the [User Guide Repository]({{ site.docs.userguide.repo }}).
- Follow the [Contributing to the User Guide]({{ site.docs.userguide.url }}/stable/contributing_docs.html).
- Check for open documentation-related issues in the [User Guide GitHub Issues]({{ site.docs.userguide.issues }}).

### API Reference
The [API reference]({{ site.docs.mdanalysis.url }}/stable/index.html) is **automatically generated** from *docstrings* in the codebase. It provides a structured reference for developers working with MDAnalysis. You can contribute by improving *docstrings* for functions, classes, and modules or updating outdated documentation.

#### How to Get Started
- Work on the [MDAnalysis Main Code Repository]({{ site.github.repo }}).
- Follow the [API Documentation Contribution Guide]({{ site.docs.userguide.url }}/stable/contributing_code.html#working-with-mdanalysis-docs).
- Browse open API documentation issues in [GitHub Issues]({{ site.github.issues }}).

## Mentoring & Outreach Programs
MDAnalysis has actively participated in various mentoring and outreach programs **to help students, researchers, and early-career contributors gain experience in open-source development and scientific software**. Through these programs, participants receive mentorship, contribute to real-world projects, and develop valuable skills in computational science and software engineering.

### Stay informed
For the latest announcements, deadlines, and project ideas, check out our [News]({{ site.baseurl }}/blog/) where we share updates on participation and application guidelines.  

We have been involved in several mentoring initiatives:

- **[Google Summer of Code (GSoC)](https://summerofcode.withgoogle.com/)** \
    A global, remote program where students contribute to open-source projects under mentorship, enhancing their programming skills while engaging with the open-source community.

- **[Google Season of Docs](https://developers.google.com/season-of-docs)** \
    A program that connects technical writers with open-source projects to improve documentation quality and accessibility.

- **[Outreachy](https://www.outreachy.org/)** \
    A diversity-focused program offering paid remote internships to individuals underrepresented in tech, providing opportunities to contribute to open-source software.

- **[Station1 Frontiers Fellowship](https://www.station1.org/sff)**  \
    A unique research and internship program designed for students in science and technology, promoting inclusive and socially impactful research.

- **[CompChemURG](https://www.bindingsites.co.uk/home)**  \
    A mentoring initiative focused on computational chemistry, supporting undergraduates and early-career researchers in gaining expertise in the field.

Many of our [core developers]({{ site.baseurl }}/pages/governance/#mdanalysis-core-developers) started through these programs &mdash; consider joining and becoming part of the MDAnalysis community!

### GSoC
MDAnalysis has participated in Google Summer of Code for several years, mentoring students on a variety of open-source projects. Explore our past projects:

- [GSoC 2025]({% post_url 2025-05-23-gsoc-students %})
- GSoC 2024 — [Small Projects]({% post_url 2024-08-30-gsoc2024-final-sm %}); [Large Projects]({% post_url 2025-02-11-gsoc2024-final-lg %})
- [GSoC 2023]({% post_url 2024-01-18-gsoc2023-final %})
- [GSoC 2022]({% post_url 2022-12-15-gsoc_outreachy2022-final %})
- [GSoC 2021]({% post_url 2021-9-14-gsoc2021-final %})
- [GSoC 2020]({% post_url 2020-09-10-gsoc2020-final %})
- [GSoC 2019]({% post_url 2019-05-24-gsoc-students %})
- [GSoC 2018]({% post_url 2018-04-26-gsoc-students %})
- [GSoC 2017]({% post_url 2017-06-04-gsoc-students %})
- [GSoC 2016]({% post_url 2016-04-26-gsoc %})

## Other Contributions
If you are interested in contributing in other ways, such as writing a **blog post**, we welcome contributions that document your experience with MDAnalysis &mdash; whether it's learning the software, participating in a user group meeting, or working on a project. Feel free to reach out via the [Community]({{ site.baseurl }}/community/) page for guidance on submitting a post.

If you are interested in improving our website, you can refer to the [README file][] in the [Website Repository][]. While the website is already maintained by the MDAnalysis team, any fixes and updates are always welcome.

[MDAnalysis Project]: https://github.com/MDAnalysis/
[GitHub discussions]: https://github.com/MDAnalysis/mdanalysis/discussions
[Website Repository]: https://github.com/MDAnalysis/MDAnalysis.github.io
[README file]: https://github.com/namiroues/MDAnalysis.github.io/blob/master/README.md

