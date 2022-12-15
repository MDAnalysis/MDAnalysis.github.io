---
layout: post
title: Outreachy Report - 2022 Improve MDAnalysis by implementing type hinting
---
## About Me
I am Uma Kadam , a Computer Science and Engineering undergraduate at Indian Institute Of Information Technology Guwahati. My interests primarily lie in exploring ML & AI which is evident from my research internship experience on NLP . Participating in numerous hackathons, some of which I won, allowed me to explore new technologies and domains in computer science. My involvement with Outreachy provided me with a first-hand introduction to Open Source and set me on the path to a successful career in technology.

[![Linkedin](https://i.stack.imgur.com/gVE0j.png) LinkedIn](https://in.linkedin.com/in/uma-kadam-7885341b0)
&nbsp;
[![GitHub](https://i.stack.imgur.com/tskMh.png) GitHub](https://github.com/umak1106)


## Outreachy

[Outreachy](https://www.outreachy.org/) is a 12+ week internship program where contributors work with an open-source organization under the guidance of experienced mentors. By providing opportunities to work with participating organizations, Outreachy supports people from underrepresented groups in technological sector. Providing mentorship to build technical skills and establishing an inclusive community that has no room for systematic bias or discrimination, it aims to help minority members pave the way into the tech industry.

## Motivation 

The object-oriented Python library MDAnalysis analyzes trajectory data derived from molecular dynamics (MD) simulations in many popular formats.
Python's dynamic nature enables us to develop with speed, flexibility, and ease of use but if you are not careful, you may trade short-term expedience for long-term lack of maintainability because of the dynamic nature of Python. Type hints are implemented with the help of [typing](https://docs.python.org/3/library/typing.html) module.

While type hints and type annotations do hint towards or indicate the appropriate types they do not enforce them. By utilizing typecheckers such as [mypy](http://mypy-lang.org/), I made sure that the code is performing what it should regarding the types passed around between functions, and the annotated function signatures only boosted the readability of the code as well as improved communication within it.

### Introducing type annotations and type hints provided a multitude of benefits like:

* The need to document the type in the docstring got eliminated.

* Datatypes got clearly defined in the code, removing any potential datatype ambiguity.

* It was beneficial in catching errors, improving linting and providing a more cleaner architecture for MDanalysis .

* It was helpful in making the code more organized and speeding up the debugging process .

* It provided us with optional static typing to leverage the best of both static and dynamic typing.

## Contributions made during Outreachy Project


* **Addition of Mypy requirements to CI pipeline**: [#3705](https://github.com/MDAnalysis/mdanalysis/pull/3705/files)

    Addition of mypy to github actions workflow provided us the ability to run a mypy check on every new pull request made and to raise appropriate warnings whenever the type hints provided were incorrect or erroneous. Added customizations like running mypy checks only on more prioritized modules and so on . The type checking in CI was made to be blocking wherever required.
    
* **Providing Type Hints for lib module and annotations for init file**: [#3823](https://github.com/MDAnalysis/mdanalysis/pull/3823/files) [#3729](https://github.com/MDAnalysis/mdanalysis/pull/3729/files) 
    
    The first PR deals with task of providing type annotations for the init file.
    
    The second PR deals with the task of providing type annotations and type hints for the lib module and usage of mypy type checker to ensure all the functions and variables are correctly type hinted. Usage of typing module and numpy.typing module for providing type hints.
    In the lib module I provided type hints for [NeighbourSearch](https://docs.mdanalysis.org/2.0.0/documentation_pages/lib/NeighborSearch.html), [PeriodicKDTree](https://docs.mdanalysis.org/2.0.0/documentation_pages/lib/pkdtree.html), [Mathematical_helper_functions](https://docs.mdanalysis.org/2.0.0/documentation_pages/lib/mdamath.html) and so on.
    
* **Providing type hints for core module**: [#3719](https://github.com/MDAnalysis/mdanalysis/pull/3719/files) [Docs](https://docs.mdanalysis.org/stable/documentation_pages/core_modules.html)
    
    Dealt with the issue caused by circular imports and provided type hints for the core module.
    
* **Providing type hints for visualization module, auxiliary module, topology module, analysis module and converters module**

    [#3781](https://github.com/MDAnalysis/mdanalysis/pull/3781)  [#3746](https://github.com/MDAnalysis/mdanalysis/pull/3746) [#3782](https://github.com/MDAnalysis/mdanalysis/pull/3782) [#3744](https://github.com/MDAnalysis/mdanalysis/pull/3774/files) [#3752](https://github.com/MDAnalysis/mdanalysis/pull/3752) [#3784](https://github.com/MDAnalysis/mdanalysis/pull/3784)
    
    Usage of np.ndarray and other type hints from the typing module for type hinting and annotating the streamlines and streamlines_3D files in the visualization module.Provide type hints, as well as use mypy to efficiently perform type checks for the different parsers present in the Topology module, and see if any unit tests have been broken as a result.Providing type hints for converters, auxiliary and other modules.


The annotations are not all merged yet, and they do not cover the full module. However, spending efforts on type annotations gave us an idea of the challenges it represents for MDAnalysis and will lead to better code in the future. Already, we identified corner cases that could have led to bugs.

## Example :

```python
    def search(self, atoms, radius, level='A'):
```

With the help of type hints this is tranformed to:

```python
    def search(self, atoms: AtomGroup, radius: float, level: str = 'A') -> Optional[Union[AtomGroup, ResidueGroup, SegmentGroup]]:
```
Here the type hints used for the input parameters of the function search inform us that  the type of Atoms is [Atomgroup](https://userguide.mdanalysis.org/1.1.1/atomgroup.html), radius is a float value and level is a string which takes default value 'A'.

The return type of the function can be None or  [AtomGroup](https://userguide.mdanalysis.org/1.1.1/atomgroup.html) or [ResidueGroup](https://docs.mdanalysis.org/1.1.1/documentation_pages/core/groups.html#MDAnalysis.core.groups.ResidueGroup) or [SegmentGroup](https://docs.mdanalysis.org/2.3.0/documentation_pages/core/groups.html#MDAnalysis.core.groups.SegmentGroup).
Here the Optional keyword suggests that the type can be None or the type enclosed within its brackets and the Union keyword suggests that the type could be anything from the types conatined within its brackets.


    
## Conclusion :
In my opinion, Outreachy provided an incredible learning opportunity for a newcomer to open source development like me. There is no doubt that the internship was an exceptional experience and I would like to extend my sincere gratitude to Outreachy, my mentors, and the MDAnalysis community for that experience. During my time working with the MDAnalysis library, I have learned a great deal not only from the technical side of working with the library, but also from experiencing the rich set of policies and guidelines that are created by the MDAnalysis community, and how they shape MDAnalysis' inner workings.
    
