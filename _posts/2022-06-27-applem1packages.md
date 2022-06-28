---
layout: post
title: Apple M1 conda packages for MDAnalysis 2.2.0
---

We now also have [conda-forge
packages](https://anaconda.org/conda-forge/mdanalysis/files) for our
[MDAnalysis 2.2.0 release]({{ site.baseurl }}{% post_url
2022-06-02-release-2.2.0 %}) that directly support the [Apple
M1](https://en.wikipedia.org/wiki/Apple_M1) ARM architecture (labelled
*osx-arm64*).

On all

* supported Python versions (**3.8, 3.9, 3.10**)
* supported Operating Systems (**Linux**,  **Windows**, **MacOS**)

you are now able to install and upgrade with `conda`

```bash
conda update -c conda-forge mdanalysis
```

For everything else about the new release, read our blog post about
[MDAnalysis 2.2]({{ site.baseurl }}{% post_url
2022-06-02-release-2.2.0 %}).


— The MDAnalysis Team


