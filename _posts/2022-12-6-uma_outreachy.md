---
layout: post
title: Outreachy Report - 2022 Improve MDAnalysis by implementing type hinting
---

## Motivation 

The object-oriented Python library MDAnalysis analyzes trajectory data derived from molecular dynamics (MD) simulations in many popular formats.
Python's dynamic nature enables us to develop with speed, flexibility, and ease of use, but if you are not careful, you can result in short-term expedience at the expense of long-term maintainability as a result of the dynamic nature of Python.

While type hints and type annotations do hint towards or indicate the appropriate types they do not enforce them. By utilizing typecheckers such as mypy, I made sure that the code is performing what it should regarding the types passed around between functions, and the annotated function signatures only boosted the readability of the code as well as improved communication within it.

### Introducing type annotations and type hints provided a multitude of benefits like:

The need to document the type in the docstring got eliminated.

Datatypes got clearly defined in the code, removing any potential datatype ambiguity.

It was beneficial in catching errors, improving linting and providing a more cleaner architecture for MDanalysis .

It was helpful in making the code more organized and speeding up the debugging process .

It provided us with optional static typing to leverage the best of both static and dynamic typing.

## Contributions made during Outreachy Project


* **Addition of Mypy requirements to CI pipeline**: [#3705](https://github.com/MDAnalysis/mdanalysis/pull/3705/files)

    Addition of mypy to github actions workflow provided us the ability to run a mypy check on every new pull request made and to raise appropriate warnings whenever the type hints provided were incorrect or errorneous. Added customizations like running mypy checks only on more prioritized modules and so on . The type checking in CI was made to be blocking wherever required.
    
* **Providing Type Hints for lib module and annotations for init file**: [#3823](https://github.com/MDAnalysis/mdanalysis/pull/3823/files) [#3729](https://github.com/MDAnalysis/mdanalysis/pull/3729/files)
    
    The first PR deals with task of providing type annotations for the init file.
    
    The second PR deals with the task of providing type annotations and type hints for the lib module and usage of mypy type checker to ensure all the functions and variables are correctly type hinted. Usage of typing module and numpy.typing module for providing type hints.
    
* **Providing type hints for core module**: [#3719](https://github.com/MDAnalysis/mdanalysis/pull/3719/files)
    
    Dealt with the issue caused by circular imports and provided type hints for the core module.
    
* **Providing type hints for visualization module, auxiliary module, topology module, analysis module and converters module**

    [#3781](https://github.com/MDAnalysis/mdanalysis/pull/3781)  [#3746](https://github.com/MDAnalysis/mdanalysis/pull/3746) [#3782](https://github.com/MDAnalysis/mdanalysis/pull/3782) [#3744](https://github.com/MDAnalysis/mdanalysis/pull/3774/files) [#3752](https://github.com/MDAnalysis/mdanalysis/pull/3752) [#3784](https://github.com/MDAnalysis/mdanalysis/pull/3784)
    
    Usage of np.ndarray and other type hints from the typing module for type hinting and annotating the streamlines and streamlines_3D files in the visualization module.Provide type hints, as well as use mypy to efficiently perform type checks for the different parsers present in the Topology module, and see if any unit tests have been broken as a result.Providing type hints for converters, auxiliary and other modules.
    
    
## Conclusion :
In my opinion, Outreachy provided an incredible learning opportunity for a newcomer to open source development like me. There is no doubt that the internship was an exceptional experience and I would like to extend my sincere gratitude to Outreachy, my mentors, and the MDAnalysis community for that experience. During my time working with the MDAnalysis library, I have learned a great deal not only from the technical side of working with the library, but also from experiencing the rich set of policies and guidelines that are created by the MDAnalysis community, and how they shape MDAnalysis' inner workings.
    

  








