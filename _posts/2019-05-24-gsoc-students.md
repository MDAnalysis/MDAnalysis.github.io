---
layout: post
title: Google Summer of Code Student 2019
---

We are happy to anounce that MDAnalysis is hosting one [GSoC][gsoc]
students for [NumFOCUS][numfocus] this year, [Ninad Bhat][ninad-gsoc] (@NinadBhat) on GitHub).

# Ninad Bhat: Better Periodic Boundary Handling

<img
src="https://ninadbhat.github.io/images/profile.png"
title="NinadBhat Suhane" alt="Ninad Bhat"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Molecular simulations are predominantly ran under periodic boundary
conditions, i.e., upon leaving one face of the simulation volume, you
re-enter in the opposite face. This can lead to molecules being split
over the periodic boundary, which requires rectification before
performing calculations. In this project, Ninad will implement
wrapping and unwrapping functionality in the various AtomGroup methods
that use the position of particles, e.g., the calculation of the
center of mass. In order to improve performance, the wrapping and
unwrapping methods will be implemented in Cython.


Ninad is a researcher in materials science at IIT Bombay (Mumbai).


Ninad will describe his progress on [his blog][ninad-blog].

— @jbarnoud @richardjgowers @micaela-matta @orbeckst (mentors)

[gsoc]: https://summerofcode.withgoogle.com
[numfocus]: https://www.numfocus.org/
[ninad-gsoc]: https://summerofcode.withgoogle.com/projects/#5319625758212096
[ninad-blog]: https://ninadbhat.github.io/
