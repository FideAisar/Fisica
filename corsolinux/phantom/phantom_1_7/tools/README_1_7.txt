In this version 1_7 the possibility has been introduced of choosing  to define 
  boson polarizations (when requested) in the center of mass system of four 
  outgoing particles, or using the previous lab definition. The four particles 
  have to be indicated by the user (in r.in) and are normally the decay 
  particles of the bosons. 
  We refer to our paper arXiv:2007.07133 for an analysis of the differences 
   between the two definitions in vector boson scattering processes.
  In 1_7 it is necessary to define in r.in a constant which multiply the 
  computed pdf scale (which of course can be 1.d0). 
  For details, see the r.in in the distribution 
    
Version 1_6 is very similar to the beta version 1_5_1 b.
   It only differs for the correction of a bug in fxn.f which could manifest 
    itself only for processes with all quarks and 4 identical fermions in the 
    final state (as dd -> ddddd~d~). 
    In the tools directory the scripts gendir.scr and setupdir*.pl 
    have been modified so that they can be used forthe CONDOR queue system 
    at Cern  

    The version 1_5_1 b was a beta version.
    Compared to the  previous 1_5_b
    - the problem,  which gave some rare events in which 
    the color flow was inconsistent with the mothers assignment, has been solved.
    - the integration of processes with b bbar inistial states has been made
     more efficient.

Version 1_6 (and of course 1_5_1_b) differs from the previous one  1_3_2 in 
particular for the introduction of the possibility  of 
   - computing only resonant W and/or Z diagrams 
   - computing polarized cross sections for W and/or Z diagrams 
   - projecting on shell the above contributions

We refer to our papers  arXiv:1710.09339 and arXiv:1907.04722  
for the  correct meaning of the above extensions and their use in physical 
studies for W's.

THe above computations can be used only for alpha_ew^6.

In the last versions there are also some additional possible cuts implemented 
and another possible scale of renormalization and factorization.


It is mandatory that one uses the r.in file present in this distribution 
(see the tools directory), modifying it, but keeping its structure. In this 
r.in there are explanations and warnings that should be carefully taken into 
account. 
We have taken away in this version some options that have not been used in 
recent years (as for instance unitarization models, user defined cuts or cuts 
stronger than those used in the cross section calculations), in order to 
simplify the use of the r.in

In the tools directory one can find several versions of the setupdir.pl file, 
perl script used to create directories, input files and batch-run files for 
computing the various processes in the first step (ionesh=0). 
setupdir2.pl has to be used when i_ccfam=1 as it produces only one representative
for processes which differ by charge and/or family conjugation. 
setupdirall.pl when one sets i_ccfam=0 as it produces all processes, also the 
charge and family conjiugated.  
The corresponding scripts setupdir2_nob.pl and setupdirall_nob.pl do not produce 
processes with external b quarks. 

iccfam=1 cannot be set if one has resonant W's. It can be set if one has only 
resonant Z's decaying to leptons (the decay particles are in this case ccfam 
invariiant), but also in this case not for the polarized cross sections.   

setupdirall.pl or setupdirall_nob.pl  has normally to be chosen (together with 
i_ccfam=0) when producing polarized polarized cross section in order not to mix 
W+ and W- polarizations. If one has final resonant W's and b's in the final 
state the ew top contributions may be very large and they can be excluded not 
using the b's or using the appropriate deltacuttop cut.   

As for all versions from 1_3, this version can be used with intel compiler ifort 
or with gfortran. We recommend the use of intel compiler when available 
as in our experience the executable compiled with gfortran is slower than the 
corresponding compiled with intel (by at least a factor two). 

The user must compile the code using the makefile.
He has to substitute the address of  PDFLIBDIR with the one he uses.
Moreover he has to choose (uncomment) the flags F77, FLAG, FLA4z... 
corresponding to the compiler he uses and comment the corresponding ones for 
the others.
One will find two possibilities for gfortran: cern lxplus and cern lxplus7. 
The default gfortran on lxplus is old (2010) and does not support for our code 
the best optimizations. In lxplus 7 there is a more recent version (2015) which 
allows them and consequently runs faster. If you use gfortran we suggest, 
if possible, to use this second compilation. If you have any doubt, just try to 
compile it this way. If it arrives at the end and produces the file phantom.exe,
 then it works. 

For other information on PHANTOM and how to use it please refer to the directory
/afs/cern.ch/work/b/ballest/public/phantom
and the readme.txt you find there.



 
