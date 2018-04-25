---
layout: post
title: Google Summer of Code Students 2018
---

We are happy to anounce that MDAnalysis is hosting two [GSoC][gsoc]
students for [NumFOCUS][numfocus] this year, [Ayush
Suhane][ayush-gsoc] (@ayushsuhane) on GitHub) and [Davide
Cruz][davide-gsoc] (@davidercruz).

# Ayush Suhane: Improve Distance Search Methods in MDAnalysis

<img
src="https://avatars0.githubusercontent.com/u/34154224?s=400&v=4"
title="Ayush Suhane" alt="Ayush Suhane"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

With the capability of multiple MD codes to easily handle milions of
atoms, a major roadblock to analysis of this vast amount of data
corresponding to positions of each atoms at every timestep is the time
to evaluate pairwise distance between multiple atoms. Almost every
operation requires the distance between the pair of atoms, fast
calculation of pairwise distance is of utmost importance. Multiple
basic analysis functions like Radial Distribution Function, Contact
Matrices, depepend very heavily on fast distance evaluations. Apart
from naive approach for pairwise calculations which scale as
$$\mathcal{O}(N^2)$$, other forms of data structures like KDTree,
Octree are sugested for faster calulations based on the
requirements. Based on the MDAnalysis, two use cases are identified as
highly used in majority of the analysis algorithms. The goal of the
project is to identify the data structure based on the requirements of
the use case and implement in the MDAnalysis library along with clear
documentations and test cases.

Ayush is ... (ADD SHORT BIO)

Ayush will describe his progress on [his blog][ayush-blog].

# Davide Cruz: Implementing on-the-fly coordinate transformations

<img
src="https://avatars0.githubusercontent.com/u/37750297?s=400&v=4"
title="Davide Cruz" alt="Davide Cruz"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Implement trajectory transformations on the MDAnalysis API, to be
called on-the-fly by the user, eliminating the requirement for
multiple intermediate steps of modifying and saving the trajectory,
and giving users a more efficient and simple workflow for simulation
data analysis.

Davide is currently on the last year of his PhD on Molecular Biosciences
at ITQB-NOVA in Lisbon, Portugal. For his thesis he is using MDAnalysis to
analyse the results of molecular dynamics simulationations and this GSoC
project is an opportunity to contribute to the community. He expects to
learn a lot about python and software development during this summer.

Davide will describe his progress on [his blog](https://davidercruz.github.io).


# Other NumFOCUS students

NumFOCUS is hosting 45 students this year for several of their supported and
affiliated projects. You can find out about the other
students
[here](https://github.com/numfocus/gsoc/blob/master/2018/accepted_student_blogs.md).

[gsoc]: https://summerofcode.withgoogle.com
[numfocus]: https://www.numfocus.org/
[ayush-gsoc]: https://summerofcode.withgoogle.com/projects/#5050592943144960
[ayush-blog]: XXX
[davide-gsoc]: https://summerofcode.withgoogle.com/projects/#5194121086500864
[davide-blog]: https://davidercruz.github.io
