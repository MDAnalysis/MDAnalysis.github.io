---
layout: post
title: Google Summer of Code Students 2022
---

We are happy to announce that MDAnalysis is hosting two [GSoC][gsoc]
students this year -- @aya9aladdin, and @BFedder. MDAnalysis has been accepted as its own
[organization with GSoC][mda-gsoc] for the third year running and we are grateful to Google for granting us
the opportunity to undertake _two very exciting GSoC projects_!

## Aya Mohamed Alaa: Context-aware guessers

<img
src="Aya's Picture Here"
title="Aya Mohamed Alaa" alt="Aya Mohamed Alaa"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Write something about your project and yourself here (don't forget to link to the project page on the GSoC website)!

Add where we can find you here! (github, twitter, others)

Add a link to your blog here!


## Bjarne Feddersen: Adding Energy Readers to MDAnalysis

<img
src="https://bfedder.github.io/assets/images/profile_photo.jpeg"
title="Bjarne Feddersen" alt="Bjarne Feddersen"
style="float: left; width: 110px; height: 110px; border-radius: 20px; border: 15px solid white" />

Bjarne is joining us to work on the [auxiliary data framework][aux-guide]. This framework allows the association of non-trajectory data with the frames of a trajectory in MDAnalysis, which could for example be used for further analyses or to filter frames of a trajectory based on non-trajectory information. The framework is currently limited to working with XVG files due to the fact that no other file readers are implemented. This will change with Bjarne's project. In particular, he will implement [new AuxReaders to read energy files][bjarne-project]. These files are produced by MD engines during a simulation and contain information on, for example, the system's kinetic and potential energy, temperature, or pressure. These terms describe important quantitites of the system and as such help with evaluation of simulations. Implementing new AuxReaders for these files in MDAnalysis will make this part of the analysis more convenient and at the same time broaden the scope of the auxiliary data framework. To further aid this second objective, Bjarne will compile the lessons he will learn while writing the new AuxReaders into a comprehensive tutorial to make future additions to the framework easier.

Bjarne is a DPhil student in [Phil Biggin's group][bigginlab] at the University of Oxford, where he uses MD simulations and other computational tools to investigate the mechanism of action of voltage-gated sodium channels. To balance his desk job he likes to spend time outdoors, and especially enjoys cycling through the green English countryside. 

Bjarne is on GitHub as @BFedder and on [LinkedIn][bjarne-linkedin]. He will be reporting on his project [on his blog][bjarne-blog].


— @jbarnoud @hmacdope @ojeda-e @IAlibay @fiona-naughton @orbeckst @lilyminium (mentors)

[gsoc]: https://summerofcode.withgoogle.com
[mda-gsoc]: https://summerofcode.withgoogle.com/programs/2022/organizations/mdanalysis
[bjarne-project]: https://summerofcode.withgoogle.com/programs/2022/projects/wbLbZmGk
[bjarne-blog]: https://bfedder.github.io
[aux-guide]: https://userguide.mdanalysis.org/stable/formats/auxiliary.html
[bigginlab]: https://bigginlab.web.ox.ac.uk
[bjarne-linkedin]: https://www.linkedin.com/in/bjarne-feddersen-407184187/
