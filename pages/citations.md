---
layout: page
title: Citations
---

MDAnalysis and the included algorithms are scientific software that
are described in academic publications. Please cite them when you use
them in published work. 

It is possible to [automatically generate a list of
references]({{site.pypi.docs}}/documentation_pages/references.html#citations-using-duecredit)
for any program that uses MDAnalysis. This list (in common reference
manager formats) contains the citations associated with the specific
algorithms and libraries that were used in the program.



## MDAnalysis (library) ##

When using MDAnalysis in published work, please cite the following two papers:

 * <a name="Gowers2016"></a>R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy,
   M. N. Melo, S. L. Seyler, D. L. Dotson, J. Domanski, S. Buchoux,
   I. M. Kenney, and
   O. Beckstein. [MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations](http://conference.scipy.org/proceedings/scipy2016/oliver_beckstein.html). In
   S. Benthall and S. Rostrup, editors, *Proceedings of the 15th Python in
   Science Conference*, pages 102-109, Austin, TX, 2016. SciPy.

 * <a name="MichaudAgrawal2011"></a>N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
   O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics
   Simulations. *J. Comput. Chem.* **32** (2011), 2319-2327,
   doi:[10.1002/jcc.21787](http://dx.doi.org/10.1002/jcc.21787).
   PMCID:[PMC3144279](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3144279/)

   (If you are reading the [HTML version of the
   paper](http://onlinelibrary.wiley.com/doi/10.1002/jcc.21787/full),
   have a look at the [paper
   errata]({{ site.baseurl }}/pages/errata).
   The free PubmedCentral manuscript
   [PMC3144279](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3144279/)
   has correct code, which can be copied and pasted.)

## <a name="IncludedAlgorithms"></a>Included algorithms ##

If you use the [RMSD alignment
code](http://docs.mdanalysis.org/documentation_pages/analysis/align.html)
that uses the [QCProt
module](http://docs.mdanalysis.org/documentation_pages/core/qcprot.html)
please also cite

 * Douglas L. Theobald. Rapid calculation of RMSD using a quaternion-based
   characteristic polynomial. *Acta Crystallographica A* **61** (2005),
   478-480. doi: [10.1107/S0108767305015266](http://doi.org/10.1107/S0108767305015266)

 * Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald. Fast determination
   of the optimal rotational matrix for macromolecular
   superpositions. *J. Comput. Chem.* **31** (2010), 1561-1563. doi:
   [10.1002/jcc.21439](http://doi.org/10.1002/jcc.21439)

If you use the helix analysis algorithm HELANAL in
[MDAnalysis.analysis.helanal](http://docs.mdanalysis.org/documentation_pages/analysis/helanal.html)
please cite

 * Bansal M, Kumar S, Velavan R. HELANAL - A program to characterise helix
   geometry in proteins. *J. Biomol. Struct. Dyn.* **17** (2000),
   811-819. doi:
   [10.1080/07391102.2000.10506570](http://doi.org/10.1080/07391102.2000.10506570)

If you use the GNM trajectory analysis code in
[MDAnalysis.analysis.gnm](http://docs.mdanalysis.org/documentation_pages/analysis/gnm.html)
please cite

 * Benjamin A. Hall, Samantha L. Kaye, Andy Pang, Rafael Perera, and Philip
   C. Biggin. Characterization of Protein Conformational States by Normal-Mode
   Frequencies. *J. Am. Chem. Soc.* **129** (2007), 11394-11401. doi:
   [10.1021/ja071797y](http://doi.org/10.1021/ja071797y)

If you use the water analysis code in
[MDAnalysis.analysis.waterdynamics](http://docs.mdanalysis.org/documentation_pages/analysis/waterdynamics.html)
please cite

 * Araya-Secchi, R., Tomas Perez-Acle, Seung-gu Kang, Tien Huynh,
   Alejandro Bernardin, Yerko Escalona, Jose-Antonio Garate, Agustin
   D. Martinez, Isaac E. Garcia, Juan C. Saez, Ruhong
   Zhou. Characterization of a novel water pocket inside the human
   Cx26 hemichannel structure. *Biophysical Journal* ***107*** (2014),
   599-612. doi: [10.1016/j.bpj.2014.05.037](http://doi.org/10.1016/j.bpj.2014.05.037)

If you use the Path Similarity Analysis (PSA) code in
[MDAnalysis.analysis.psa](http://docs.mdanalysis.org/documentation_pages/analysis/psa.html)
please cite

 * Seyler SL, Kumar A, Thorpe MF, Beckstein O (2015) Path Similarity
   Analysis: A Method for Quantifying Macromolecular Pathways. *PLoS
   Comput Biol* **11**(10): e1004568. doi: [10.1371/journal.pcbi.1004568](http://dx.doi.org/10.1371/journal.pcbi.1004568)

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
2.7. MDAnalysis has been successfully used on Linux and Mac OS X.

The **MDAnalysis 'Atom' Logo** was designed by **Christian Beckstein** and is
licensed under a [Creative Commons Attribution-NoDerivs 3.0 Unported
License](http://creativecommons.org/licenses/by-nd/3.0/).
