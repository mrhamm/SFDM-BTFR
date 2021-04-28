***Contact Info*** 
Please contact Bih@duke.edu if you have any questions in regards to 
    this repository.
***Contact Info***

***Publication Info***
This repository contains a set of specialized mathematical routines that 
    correspond to the publication at https://arxiv.org/abs/2004.07792
***Publication Info***

***Libraries and Routines***
SSSStateGenerator.py  contains the necessary functionality to produce solutions
    to the spherically symmetric and static Einstein Klein Gordon Equations.  
    EDIT THIS FILE AT YOUR OWN RISK, the majority of the other routines depend 
    on the format of this file.  We STRONGLY SUGGEST read chapter 3 and the appendix
    of the arxiv publication listed above before making any edits.  

SPARCRotationCurves.py and SPARCTullyFisherRelation.py provide a PANDAS dataframe 
    of the relevant observational data from the SPARC survey - further details of 
    SPARC can be found in these files.  

DiskFractionSolver.py is a routine that corresponds to Chapter 3 of the arxiv publication.  
    This will take the dark matter only states from the SSSStateSolver package, and 
    use a continuation method to include the effects of ordinary visible matter.  

BoundaryProblemSolver.py will apply a scaling relationship to a family of solutions
    and output a new set of solutions with a tully-fisher like trend.  This corersponds
    to the Amplitude-Wavelength boundary condtion from the publication.  

TullyFisherFittingRoutine.py contains a simplified routine that combines the entirety 
   of this repository to produce a basic fit to the observed BTFR.  
***Libraries and Routines***


*************RUNTIME NOTICE***************
 The routines contained in this repository can take an extremely long time to 
    resolve and converge to their final results.  The runtime will increase
    quadratically with the total number of static states calculated.  Computations
    of the first 100 states can take upwards of 48 hours
*************RUNTIME NOTICE*****************