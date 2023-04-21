---
layout: post
title: SolvationAnalysis Published in JOSS
---

We're happy to announce that the [SolvationAnalysis][solvation] Python package is now [published][paper] in the [Journal of Open Source Software][joss] (JOSS):

Orion Archer Cohen, Hugo Macdermott-Opeskin, Lauren Lee, Tingzheng Hou, Kara D. Fong, Ryan Kingsbury, Jingyang Wang, and Kristin A. Persson, (2023). *SolvationAnalysis: A Python toolkit for understanding liquid solvation structure in classical molecular dynamics simulations.* **Journal of Open Source Software**, 8(84), 5183, https://doi.org/10.21105/joss.05183



Originally developed by @orioncohen as a [GSoC 2021 project][gsocblog], SolvationAnalysis builds on MDAnalysis and pandas to make analyzing solvation structure much, much easier. With a few lines of code, you can extract key solvation properties from any molecular dynamics simulation. With a few more, you can visualize solvation trends within and between different solutions. Further, SolvationAnalysis is designed with extensibility in mind, exposing a core representation of solvation that users can use for brand new analyses. To get started with using SolvationAnalysis in your workflow, check out the documentation and tutorials on [readthedocs][docs].

![SolvationAnalysis summary figure](https://github.com/MDAnalysis/solvation-analysis/raw/main/joss_paper/summary_figure.jpg)

[solvation]: https://github.com/MDAnalysis/solvation-analysis
[paper]: https://joss.theoj.org/papers/10.21105/joss.05183
[joss]: https://joss.theoj.org/
[gsocblog]: https://www.mdanalysis.org/2021/09/02/gsoc-final-report-orion/
[docs]: https://solvation-analysis.readthedocs.io/en/latest/
