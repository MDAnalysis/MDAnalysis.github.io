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
src="https://user-images.githubusercontent.com/27581535/170899701-bdb01612-8764-4d61-97f8-73777f2f1d44.jpg"
title="Aya Mohamed Alaa" alt="Aya Mohamed Alaa"
style="float: left; width: 90px; height: 110px; border-radius: 20px; border: 15px solid white"/>

MDAnalysis supports many molecular dynamics simulation files and force fields, which brings the challenge of tailoring the guessing methods to each type of format. Aya is going to work on building a context-aware guesser API [Context-aware guesser (implementing implementing PDB and Martini guessers)](https://summerofcode.withgoogle.com/programs/2022/projects/B1Y0nTh2) that can provide a way for the user to pass the context of his work conveniently through an API so that the guessing of missing attributes would be more accurate and reliable. 

Aya is pharmaceutical sciences graduate from Ain-Shams University (Egypt) and a premaster student at Arab Academy for science and technology (Egypt). She is interested in working on the link between computer science and biological sciences as she believes in the importance of finding a way to look at biological questions from the computational and engineering point of view. On the weekends, Aya enjoys outdoor sports, especially football, and tennis, additionally, Aya enjoys swimming as it helps reduce the feeling of stress (Egypt's sun will force you to spend the whole day at the pool though).

You can find Aya on [GitHub](https://github.com/aya9aladdin), [Twitter](https://twitter.com/AyaSalim0909), and [Facebook](https://www.facebook.com/)

This is Aya's personal blog for documenting her summer with MDAnalysis [Blog](https://sites.google.com/pharma.asu.edu.eg/aya-gsoc/home) 






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
