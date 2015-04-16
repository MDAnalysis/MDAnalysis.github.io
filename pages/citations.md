---
layout: page
title: Citations
---

When using MDAnalysis in published work, please cite

 * <a name="MichaudAgrawal2011"></a>N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
   O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics
   Simulations. *J. Comput. Chem.* **32** (2011), 2319-2327,
   doi:[10.1002/jcc.21787](http://dx.doi.org/10.1002/jcc.21787). 
   PMCID:[PMC3144279](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3144279/)

   (If you are reading the [HTML version of the
   paper](http://onlinelibrary.wiley.com/doi/10.1002/jcc.21787/full),
   have a look at the [paper
   errata]({{site.baseurl}}pages/errata).
   The free PubmedCentral manuscript
   [PMC3144279](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3144279/)
   has correct code, which can be copied and pasted.)

## <a name="IncludedAlgorithms"></a>Included algorithms ##

If you use the [RMSD alignment
code](https://pythonhosted.org/MDAnalysis/documentation_pages/analysis/align.html)
that uses the [QCProt
module](https://pythonhosted.org/MDAnalysis/documentation_pages/core/qcprot.html)
please also cite

 * Douglas L. Theobald. Rapid calculation of RMSD using a quaternion-based
   characteristic polynomial. *Acta Crystallographica A* **61** (2005),
   478–480.

 * Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald. Fast determination
   of the optimal rotational matrix for macromolecular
   superpositions. *J. Comput. Chem.* **31** (2010), 1561-1563.

If you use the helix analysis algorithm HELANAL in
[MDAnalysis.analysis.helanal](https://pythonhosted.org/MDAnalysis/documentation_pages/analysis/helanal.html)
please cite

 * Bansal M, Kumar S, Velavan R. HELANAL - A program to characterise helix
   geometry in proteins. *J. Biomol. Struct. Dyn.* **17** (2000), 811-819

If you use the GNM trajectory analysis code in
[MDAnalysis.analysis.gnm](https://pythonhosted.org/MDAnalysis/documentation_pages/analysis/gnm.html)
please cite

 * Benjamin A. Hall, Samantha L. Kaye, Andy Pang, Rafael Perera, and Philip
   C. Biggin. Characterization of Protein Conformational States by Normal-Mode
   Frequencies. *J. Am. Chem. Soc.* **129** (2007), 11394-11401.

Thanks!


## Acknowledgements

MDAnalysis was originally inspired by the Schulten Group's
[MDTools](http://www.ks.uiuc.edu/Development/MDTools/) for Python, and the DCD
reading code is derived from VMD's
[catdcd](http://www.ks.uiuc.edu/Development/MDTools/catdcd/). MDAnalysis is GPL
licensed, except for some 3rd party code that is included under GPL-compatible
licenses; for instance the dcd reading code is under the [UIUC Open Source
Licence](http://www.ks.uiuc.edu/Development/MDTools/catdcd/license.html). See
the files AUTHORS and LICENSE in the distribution for details.

Some time-critical routines are written in C or [cython](http://cython.org) and
require a working C compiler. The minimum required version of Python is
2.6. MDAnalysis has been successfully used on Linux and Mac OS X.

The **MDAnalysis 'Atom' Logo** was designed by **Christian Beckstein** and is
licensed under a [Creative Commons Attribution-NoDerivs 3.0 Unported
License](http://creativecommons.org/licenses/by-nd/3.0/).

