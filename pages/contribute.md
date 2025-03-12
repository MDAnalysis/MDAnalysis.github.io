---
layout: page
title: Contribute
order: 6
---

MDAnalysis is an open-source project that welcomes and encourages contributions from the community. Whether you're a developer, documentation writer, or interested in improving the website, there are many ways to get involved.

You can contribute by **developing new features, fixing bugs, improving documentation, enhancing tutorials, contributing to MDAKits, improving the website, or engaging in discussions**. 

If you're new, you can start by participating in community conversations and mentoring programs. For more information on how you can participate, check out the [Community]({{ site.baseurl }}/pages/community/) and [Events]({{ site.baseurl }}/pages/events/) pages.

This page provides an overview of how you can contribute to the [Codebase](#contributing-to-the-codebase), [MDAKits](#contributing-to-mdakits) or [Documentation](#contributing-to-the-documentation).

## Contributing to the Codebase

If you are interested in developing new features, fixing bugs, or improving performance, contributing to the MDAnalysis main codebase is a great way to get involved. 

### Main Codebase Repository
MDAnalysis is developed and maintained in the [MDAnalysis Main Code Repository]({{ site.github.repo }}), which contains:

- The core library for molecular dynamics analysis.
- The [API documentation]({{ site.docs.mdanalysis.url }}/stable/index.html), generated from in-code docstrings.
- [Open issues]({{ site.github.issues }}), where contributions such as bug fixes and feature implementations are tracked.

### How to Get Started
- Read the [Contributing to the Main Codebase]({{ site.docs.userguide.url }}/stable/contributing_code.html#working-with-mdanalysis-code) guide, which outlines all the necessary steps.
- Browse [GitHub Issues]({{ site.github.issues }}) to find open tasks. Look for "good first issue" if you're a beginner. You can also report bugs or suggest improvements there.
- Follow best practices, including
    - [Writing and Running Tests]({{ site.docs.userguide.url }}/stable/testing.html#testing) to ensure your changes work correctly.
    - [Updating the Documentation]({{ site.docs.userguide.url }}/stable/contributing_code.html#building-code-documentation) when making modifications.

## Contributing to MDAKits

[MDAKits]({{ site.mdakits.registry }}) are **community-driven projects** that extend MDAnalysis functionality. They allow researchers and developers to build tools tailored for specific workflows while staying within the MDAnalysis ecosystem.

### MDAKits Registry
The [MDAKits Registry]({{ site.mdakits.registry }}) lists active projects contributed by the community. These projects cover various applications, such as enhanced trajectory analysis, molecular simulations, and integration with machine learning frameworks.

### How to Get Started
- **Explore existing MDAKits** in the [MDAKits Registry]({{ site.mdakits.registry }}) to find projects that match your interests. Each entry includes functionality, installation instructions, and contribution guidelines.
- **Find open issues** in [MDAKits GitHub Issues]({{ site.mdakits.issues }}).
- **Propose a new MDAKit** by following the [MDAKits Guide]({{ site.mdakits.guide }}).
- **Contribute to an existing MDAKit** by collaborating with maintainers to improve functionality, fix bugs, or add new features.

MDAKits follow MDAnalysis coding standards but operate independently. Contributions to these projects are highly encouraged and help expand the MDAnalysis ecosystem.

## Contributing to the Documentation

Good documentation is essential to making MDAnalysis accessible to users and developers. MDAnalysis maintains two types of documentation: the [User Guide]({{ site.docs.userguide.url }}) and the [API Reference]({{ site.docs.mdanalysis.url }}/stable/index.html).

### User Guide

The [User Guide]({{ site.docs.userguide.url }}/stable/index.html) provides tutorials and explanatory content to help users install and learn MDAnalysis. You can contribute by updating unclear explanations, adding new tutorials or examples, or fixing outdated information. 

### How to Get Started
- Work on the [User Guide Repository]({{ site.docs.userguide.repo }}).
- Follow the [Contributing to the User Guide]({{ site.docs.userguide.url }}/stable/contributing_docs.html).
- Check for open documentation-related issues in the [User Guide GitHub Issues]({{ site.docs.userguide.issues }}).

### API Reference

The [API reference]({{ site.docs.mdanalysis.url }}/stable/index.html) is **automatically generated** from *docstrings* in the codebase. It provides a structured reference for developers working with MDAnalysis. You can contribute by improving *docstrings* for functions, classes, and modules or updating outdated documentation.

To get started:
- Work on the [MDAnalysis Main Code Repository]({{ site.github.repo }}).
- Follow the [API Documentation Contribution Guide]({{ site.docs.userguide.url }}/stable/contributing_code.html#working-with-mdanalysis-docs).
- Browse open API documentation issues in [GitHub Issues]({{ site.github.issues }}).

## Other Contributions

If you are interested in contributing in other ways, such as improving the website or blog, you can refer to the [README file][] in the [Website Repository][]. While the website is already maintained by the MDAnalysis team, minor fixes and updates are always welcome. 

If you'd like to write a **blog post**, feel free to reach out via the [Community]({{ site.baseurl }}/pages/community/) page.

[Website Repository]: https://github.com/MDAnalysis/MDAnalysis.github.io
[README file]: https://github.com/namiroues/MDAnalysis.github.io/blob/master/README.md


