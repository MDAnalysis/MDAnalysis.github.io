---
layout: post
title: Improvements in distance search methods
---

We are pleased to announce another successful year of Google Summer of Code with [NumFOCUS][] organization. Thanks to [Richard Gowers][] and [Jonathan Barnoud][] for mentoring the GSoC students. This year, one of the projects, was to improve the performance of pairwise distance computations, which is used quite frequently in MDAnalysis in different forms. MDAnalysis v 0.19.0 and higher will include an internal function to speed up the neighbor search computations using automatic method selection. As a result, a user need not bother with the technical details of implementations of different search algorithms. A flexible functionality with easily extendible interface ``capped_distance`` is introduced in ``MDAnalysis.lib.distances`` for this purpose. Aside from the major highlight of ``capped_distance`` and its improvements, additional functionality of ``augment_coordinates`` is also implemented for any curious user to implement different neighbor search algorithm with periodic boundary conditions.

One of the major bottleneck in various analysis routines in MDAnalysis (and typically in Molecular Dynamics studies) is the evaluation of pairwise distances among the particles. The primary problem revolves around fixed radius neighbor search algorithms. MDAnalysis offers a suite of algorithms including brute force method, tree-based binary search algorithms to solve such problems. While these methods are suitable for a variety of analysis functions using pairwise distances in MDAnalysis, one of the question was whether one can improve the performance of distance calculations using other established neighbor search methods.

This question led to the inception of Google Summer of Code [project][] with [NumFOCUS][]. [Ayush Suhane][] completed the project and was able to demonstrate performance improvements for specific cases of distance selections, identification of bonds and Radial distribution function in the analysis module of MDAnalysis. More details on the commit history, PR's and blog posts can be found in the final [report][] submitted to GSoC. Real-time benchmarks for specific modules in MDAnalysis can be found [here](https://www.mdanalysis.org/benchmarks/). 

The major highlight of the project is the introduction of ``capped_distance`` which allows automatic selection of methods based on predefined set of rules to evaluate pairs of atoms in the neighborhood of any particle. It allows a user-friendly interface for the developers to quickly implement any new algorithm throughout MDAnalysis modules. To test any new algorithm, one must comply with the following protocol:

```python
def newmethod_capped(reference, configuration, max_cutoff, min_cutoff=None, box=None, return_distance=True):
    """
        An Algorithm to evaluate pairs between reference and configuration atoms
        and corresponding distances
    """
    return pairs, distances
```

Once the method is defined, register the function name in ``_determine_method`` in ``MDAnalysis.lib.distances`` as:

```python
methods = {'bruteforce': _bruteforce_capped,
           'pkdtree': _pkdtree_capped,
           'nsgrid': _nsgrid_capped,
           'newmethod': newmethod_capped}
```
That's it. The new method is ready to be tested across functions which use ``capped_distance``. For any specific application, it can be called as ``capped_distance(ref, conf, max_dist, method=newmethod)`` from the function.

As mentioned above, MDAnalysis offers support of three different algorithms namely [bruteforce][] which is a naive pairwise distance calculation algorithm and implemented in MDAnalysis even for parallel execution, [pkdtree][] is a wrapper method around binary tree search algorithm, [nsgrid][] is an implementation of cell-list algorithm. During the tenure of GSoC'18, an additional method ``nsgrid`` is implemented in MDAnalysis with the help of [Sebastien Buchoux][]. For more information, the reader is encouraged to read the [blog], which include detailed information about different algorithms and their implementation.

While implementing any new algorithm for Molecular dynamics trajectories, one additional requirement is to handle the periodic boundary conditions. A combination of versatile function ``augment_coordinates`` and ``undo_augment`` can be used with any algorithm to handle PBC. The main idea is to extend the box by generating duplicate particles in the vicinity of the box by ``augment_coordinates``. These duplicates, as well as the original particles, can now be used with any algorithm to evaluate the nearest neighbors. After the operation, the duplicate particles can be reverted back to their original particle indices using ``undo_augment``. These functions are available in ``MDAnalysis.lib._augment``. We encourage the interested readers to try different algorithms using these functions. Hopefully, you can help us improve the performance further with your feedback. As a starting point, the skeleton to enable PBC would take the following form:

```python
def newmethod_search(coords, centers, radius, box=None):
    aug, mapping = augment_coordinates(coords, box, radius)
    all_coords = no.concatenate([coords, aug])
    """
        Perform operations for distance evaluations
        with **all_coords** using the new algorithm 
        and obtain the result in indices
    """
    indices = undo_augment(indices, mapping, len(coords))
    return indices
```

Finally, this function can be tested with ``capped_distance`` to check the performance against already implemented algorithms in MDAnalysis.

This was a flavor of what work was done during GSoC'18. Apart from performance improvements, it is envisioned that this internal functionality will reduce the burden from the user to understand all the technical details of distance search algorithms and instead allow a user to focus on their analysis, as well as allow future developers to easily implement any new algorithm which can exceed the present performance benchmarks.

As a final note, we managed to get an improvement of ~ 2-3 times in Radial Distribution Function computation, ~ 10 times in identification of bonds, and ~ 10 times in distance based selections for the already existing benchmarks in MDAnalysis. The performance is also found to improve with larger datasets but is not reported in benchmarks. Any motivated reader is welcome to submit their feedbacks about the performance of the above-mentioned functions on their data, and/or a benchmark which we would be happy to showcase to the world.


[project]: https://summerofcode.withgoogle.com/projects/#5050592943144960 
[NumFOCUS]: https://numfocus.org/
[Ayush Suhane]: https://github.com/ayushsuhane
[report]: https://gist.github.com/ayushsuhane/fd114cda20e93b0f61a8acb6d25d3276
[bruteforce]: http://www.csl.mtu.edu/cs4321/www/Lectures/Lecture%206%20-%20Brute%20Force%20Closest%20Pair%20and%20Convex%20and%20Exhausive%20Search.htm
[pkdtree]: https://en.wikipedia.org/wiki/K-d_tree
[nsgrid]: https://en.wikipedia.org/wiki/Cell_lists
[blog]: https://ayushsuhane.github.io/
[Sebastien Buchoux]: https://github.com/seb-buch
[Richard Gowers]: https://github.com/richardjgowers
[Jonathan Barnoud]: https://github.com/jbarnoud