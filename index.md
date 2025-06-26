---
layout: default
title: MDAnalysis
---


<img src="{{ site.baseurl }}/public/mdanalysis-logo_square.png"
style="float: right" alt="MDAnalysis" width="30%"/>

 The **MDAnalysis Project** develops open-source tools for the analysis of molecular simulation data, providing researchers with efficient and accessible solutions for studying molecular structures and dynamics.

 MDAnalysis is a fiscally sponsored project of [NumFOCUS][], a nonprofit that promotes open practices in research, data, and scientific computing.
 
At the core of the project is *MDAnalysis*, **an open-source Python library** for analyzing molecular dynamics (MD) trajectories across [multiple formats]({{ site.docs.userguide.url }}/formats/index.html#formats). *MDAnalysis* enables seamless reading and writing of simulation data, allowing users to efficiently analyze molecular structures and dynamics, including particle-based trajectories and individual coordinate frames (e.g., biomolecules in the [PDB format][]).

With *MDAnalysis*, you can access atomic coordinates as [NumPy][] arrays, providing a flexible and efficient framework for complex analysis tasks. The library includes powerful [atom selection commands]({{ site.docs.userguide.url }}/selections.html) for extracting subsets of structures and supports trajectory transformations (e.g., fitting to a reference structure).

Learn more about our **mission, development, team, and governance** on the [About MDAnalysis]({{ site.baseurl }}/about/) page.

## Get Started  

If you're new to *MDAnalysis*, explore the following resources to install the `MDAnalysis` package, access tutorials, and dive into the documentation:

- [Getting Started]({{ site.baseurl }}/getting_started/)  
- [Learning MDAnalysis]({{ site.baseurl }}/learning_MDAnalysis/)  
- [Documentation]({{ site.baseurl }}/documentation/)  

## Community 

MDAnalysis is driven by an active **community of users and contributors**. Stay updated and get involved through the following pages:  

<!-- TODO: Add link [Get Involved]({{ site.baseurl }}/pages/about/#get-involved)
 -->
- [Community]({{ site.baseurl }}/community/) &mdash; Get involved, ask questions, and collaborate.
- [News]({{ site.baseurl }}/blog) &mdash; Stay updated with development news and community highlights.
- [Events]({{ site.baseurl }}/events/) &mdash; Join workshops, conferences, and mentoring programs.

## Contributing 

Want to contribute to MDAnalysis? Here’s how:  

- [Contribute]({{ site.baseurl }}/contribute/) &mdash; Learn how to help improve MDAnalysis through coding, documentation, or discussions.  
- [MDAKits and MDA-based tools]({{ site.baseurl }}/mdakits/) &mdash; Build and extend MDAnalysis with MDAKits.  
- [Support MDAnalysis financially](#funding--support) &mdash; Donate via NumFOCUS to sustain the project.


## Citing MDAnalysis

If you use MDAnalysis in your research, please cite it appropriately and consider displaying our MDAnalysis badge in your projects. All citation formats and badge instructions are available in our [Citation guidelines]({{ site.baseurl }}/citations/).

Additionally, consider displaying our badge in your projects that use MDAnalysis. We provide several [embedding markup examples]({{ site.baseurl }}/citations/#powered-by-mdanalysis).

[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)

## Funding & Support

If you find MDAnalysis useful and would like to support its continued development, please consider making a donation. You can learn more about out Sponsors in our [Funding]({{ site.baseurl }}/pages/funding/) page.

{{ site.numfocus.donate_button }}

<small>
    Donations are made through [our fiscal sponsor][], [NumFOCUS][], which is a 501(c)(3) non-profit charity in the United States; as such, donations to NumFOCUS are tax-deductible as allowed by law.  As with any donation, you should consult with your personal tax adviser or the IRS about your particular tax situation.
</small>

[NumFOCUS]: https://www.numfocus.org
[our fiscal sponsor]: {{site.baseurl}}/about#partners
[NumPy]: https://numpy.org/
[PDB format]: https://www.wwpdb.org/documentation/file-format
