#  $*  is name of the target with the suffix deleted
#  $@  is name of the target



#PDFLIBDIR = /afs/cern.ch/work/b/ballest/public/phantom/LHAPDF-6.1.5_work/lib
#PDFLIBDIR = /home/phantom/phantom/LHAPDF-6.1.5_work/lib
PDFLIBDIR = /Users/maina/local/LHAPDF-6.2.0_work/


#The following  definitions are to be used for compilation on
# a recent version Intel compiler. (cern, torino 03)
#  They must be commented for other compilers.

#03
#F77 = ifort
#FLAG = -c -save -traceback
#FLAG4z = -c -save   -traceback
#FLAGggzww = -c -save  -traceback
#FLAGgg3z = -c -save  -traceback
#FLAG2z2wqcd = -c -save -O1   -traceback
#FLAG4zqcd = -c -save  -traceback
#LDFLAG =   -L$(PDFLIBDIR) -lLHAPDF

#The following  definitions are to be used for compilation on
# an older version Intel compiler. (torino 02)
#  They must be commented for other compilers.

#02
# F77 = ifort
# FLAG = -c -save -traceback
# FLAG4z = -c -save -O1  -traceback
# FLAGggzww = -c -save -O1 -traceback
# FLAGgg3z = -c -save -O1 -traceback
# FLAG2z2wqcd = -c -save -O1   -traceback
# FLAG4zqcd = -c -save -O1 -traceback
# LDFLAG =   -L$(PDFLIBDIR) -lLHAPDF



# The following definitions are to be used for a recent version
#  of gfortran  (e.g. lxplus at cern)

F77 = gfortran
FLAG = -c -fno-automatic -O3  
FLAG4z = -c -fno-automatic -O3  
FLAG2z2wqcd=-c -fno-automatic -O3 
FLAG4zqcd = -c -fno-automatic -O3 
FLAGggzww = -c -fno-automatic -O3 
FLAGgg3z = -c -fno-automatic -O3 
LDFLAG =   -L$(PDFLIBDIR) -lLHAPDF

# For the same gfortran, with also debugging options

#F77 = gfortran
#FLAG = -c -fno-automatic -O3  -g   -ffpe-trap=invalid,zero,overflow,underflow
#FLAG4z = -c -fno-automatic -O3  -g  -ffpe-trap=invalid,zero,overflow,underflow
#FLAG2z2wqcd=-c -fno-automatic -O3 -g  -ffpe-trap=invalid,zero,overflow,underflow
#FLAG4zqcd = -c -fno-automatic -O3 -g  -ffpe-trap=invalid,zero,overflow,underflow
#FLAGggzww = -c -fno-automatic -O3 -g  -ffpe-trap=invalid,zero,overflow,underflow
#FLAGgg3z = -c -fno-automatic -O3 -g   -ffpe-trap=invalid,zero,overflow,underflow
#LDFLAG =   -L$(PDFLIBDIR) -lLHAPDF



# The following definitions are to be used for an old version
#  of gfortran 

#F77 = gfortran
#FLAG = -c -fno-automatic   -O3
#FLAG4z = -c -fno-automatic  -O1
#FLAGggzww = -c -fno-automatic -O3  
#FLAGgg3z = -c -fno-automatic -O0 
#FLAG2z2wqcd = -c -fno-automatic -O3
#FLAG4zqcd = -c -fno-automatic -O0  
#LDFLAG =   -L$(PDFLIBDIR) -lLHAPDF



#LDFLAG =
#.f.o: ; $(F77) $(FLAG)  $*.f

OBJS =  phantom.o phavegas.o oneshot.o phread.o util.o \
	phsp_7.o  phsp_ini7.o \
	phsp1_1_4_multi7_c.o  phsp1_1_4_multi7_cjac.o \
	phsp2_4_multi7.o phsp2_4_multi7jac.o \
	phsp1_1_31_multi_c.o phsp1_1_31_multi_cjac.o \
	phsp1_2_3_multi5_c.o phsp1_2_3_multi5_cjac.o\
	phsp3_3_multi5.o phsp3_3_multi5jac.o \
	phsp1_5to1_4to31_multi_c.o phsp1_5to1_4to31_multi_cjac.o \
	phsp2_4to31_multi5.o phsp2_4to31_multi5jac.o \
	coupling.o bernoul.o readinput.o\
	procini.o procextraini.o proc.o \
	integ.o fxn.o extrema.o cuts.o \
	ccfcsym.o perm.o coleval.o colevalew.o isign.o LHAFileInit.o\
	fourw.o twoztwow.o fourz.o \
	fourwqcd.o twoztwowqcd.o amp8fqcd.o\
	perm2g.o amp2g.o ampem.o ggzww.o isr.o \
	pythia_interface.o pythia6225.o \
        circe.o storeLH.o \
	fourz_massless.o twoztwow_massless.o \
	fourzqcd_massless.o twoztwowqcd_massless.o \
	gg3z_massless.o ggzww_massless.o \
        unit_amp.o unit.o pol.o dpa_proj.o singleosp.o


OBJS1 =  fourzqcd.o gg3z.o

# This pseudo target prevents from stopping for some error code different
#   from 0 . In particular here I am not linking circe and it would stop
#   before executing the last target wphact:

.IGNORE:

#
### this target is there just to make wphact.exe executable
#

phantom: phantom.exe
	chmod +x phantom.exe

#
### the main command for the makefile
#

phantom.exe : $(OBJS) $(OBJS1)
	$(F77) -o $@ $(LDFLAG) $(OBJS) $(OBJS1) $(LIB) $(LDFLAG)

extrema.o  :  extrema.f common.h common_subproc.h common_cut.h
	$(F77)  $(FLAG)  $*.f
fxn.o :  fxn.f common.h common_subproc.h common_cut.h
	$(F77)  $(FLAG)  $*.f
readinput.o :  readinput.f common.h common_subproc.h common_cut.h
	$(F77)  $(FLAG)  $*.f

cuts.o :  cuts.f common.h  common_subproc.h
	$(F77)  $(FLAG)  $*.f
integ.o  :  integ.f  common.h  common_subproc.h
	$(F77)  $(FLAG)  $*.f
oneshot.o  :  oneshot.f  common.h  common_subproc.h
	$(F77)  $(FLAG)  $*.f
procextraini.o  :  procextraini.f  common.h  common_subproc.h
	$(F77)  $(FLAG)  $*.f
proc.o   :  proc.f  common.h  common_subproc.h
	$(F77)  $(FLAG)  $*.f
procini.o   :  procini.f  common.h  common_subproc.h
	$(F77)  $(FLAG)  $*.f
pythia_interface.o  :  pythia_interface.f common.h
	$(F77)  $(FLAG)  $*.f
pythia6225.o  :  pythia6225.f
	$(F77)  $(FLAG)  $*.f
phantom.o  :  phantom.f common.h
	$(F77)  $(FLAG)  $*.f
coupling.o  :  coupling.f common.h
	$(F77)  $(FLAG)  $*.f
LHAFileInit.o  :  LHAFileInit.f common.h
	$(F77)  $(FLAG)  $*.f
fourw.o  :  fourw.f common.h
	$(F77)  $(FLAG)  $*.f
twoztwow.o  :  twoztwow.f common.h
	$(F77)  $(FLAG)  $*.f
fourwqcd.o  :  fourwqcd.f  common.h
	$(F77)  $(FLAG)  $*.f
amp8fqcd.o :  amp8fqcd.f common.h
	$(F77)  $(FLAG)  $*.f
amp2g.o  :  amp2g.f common.h
	$(F77)  $(FLAG)  $*.f
phavegas.o :  phavegas.f
	$(F77)  $(FLAG)  $*.f
phread.o   :  phread.f
	$(F77)  $(FLAG)  $*.f
util.o     :  util.f
	$(F77)  $(FLAG)  $*.f
phsp1_1_4_multi7_c.o  :  phsp1_1_4_multi7_c.f
	$(F77)  $(FLAG)  $*.f
phsp1_1_4_multi7_cjac.o  :  phsp1_1_4_multi7_cjac.f
	$(F77)  $(FLAG)  $*.f
phsp2_4_multi7.o  :  phsp2_4_multi7.f
	$(F77)  $(FLAG)  $*.f
phsp2_4_multi7jac.o  :  phsp2_4_multi7jac.f
	$(F77)  $(FLAG)  $*.f
phsp1_1_31_multi_c.o  :  phsp1_1_31_multi_c.f
	$(F77)  $(FLAG)  $*.f
phsp1_1_31_multi_cjac.o  :  phsp1_1_31_multi_cjac.f
	$(F77)  $(FLAG)  $*.f
phsp1_2_3_multi5_c.o   :  phsp1_2_3_multi5_c.f
	$(F77)  $(FLAG)  $*.f
phsp1_2_3_multi5_cjac.o  :  phsp1_2_3_multi5_cjac.f
	$(F77)  $(FLAG)  $*.f
phsp3_3_multi5.o  :  phsp3_3_multi5.f
	$(F77)  $(FLAG)  $*.f
phsp3_3_multi5jac.o  :  phsp3_3_multi5jac.f
	$(F77)  $(FLAG)  $*.f
phsp1_5to1_4to31_multi_c.o  :  phsp1_5to1_4to31_multi_c.f
	$(F77)  $(FLAG)  $*.f
phsp1_5to1_4to31_multi_cjac.o  :  phsp1_5to1_4to31_multi_cjac.f
	$(F77)  $(FLAG)  $*.f
phsp2_4to31_multi5.o  :  phsp2_4to31_multi5.f
	$(F77)  $(FLAG)  $*.f
phsp2_4to31_multi5jac.o  :  phsp2_4to31_multi5jac.f
	$(F77)  $(FLAG)  $*.f
phsp_7.o  :  phsp_7.f
	$(F77)  $(FLAG)  $*.f
phsp_ini7.o  :  phsp_ini7.f
	$(F77)  $(FLAG)  $*.f
bernoul.o   :  bernoul.f
	$(F77)  $(FLAG)  $*.f
ccfcsym.o  :  ccfcsym.f
	$(F77)  $(FLAG)  $*.f
perm.o  :  perm.f
	$(F77)  $(FLAG)  $*.f
coleval.o  :  coleval.f
	$(F77)  $(FLAG)  $*.f
colevalew.o  :  colevalew.f
	$(F77)  $(FLAG)  $*.f
isign.o  :  isign.f
	$(F77)  $(FLAG)  $*.f
perm2g.o  :  perm2g.f
	$(F77)  $(FLAG)  $*.f
ampem.o  :  ampem.f common.h common_unitarization.h
	$(F77)  $(FLAG)  $*.f
isr.o  :  isr.f
	$(F77)  $(FLAG)  $*.f
circe.o :  circe.f
	$(F77)  $(FLAG)  $*.f
storeLH.o :  storeLH.f common.h common_subproc.h
	$(F77)  $(FLAG)  $*.f

fourz_massless.o  :  fourz_massless.f common.h
	$(F77)  $(FLAG)  $*.f
twoztwow_massless.o  :  twoztwow_massless.f common.h
	$(F77)  $(FLAG)  $*.f
fourzqcd_massless.o  :  fourzqcd_massless.f common.h
	$(F77)  $(FLAG)  $*.f
twoztwowqcd_massless.o  :  twoztwowqcd_massless.f common.h
	$(F77)  $(FLAG)  $*.f


gg3z_massless.o  :  gg3z_massless.f common.h
	$(F77)  $(FLAGgg3z)  $*.f
ggzww_massless.o  :  ggzww_massless.f common.h
	$(F77)  $(FLAGggzww)  $*.f


fourz.o  :  fourz.f common.h
	$(F77)  $(FLAG4z)  $*.f
ggzww.o   :  ggzww.f  common.h
	$(F77)  $(FLAGggzww)  $*.f
gg3z.o  : gg3z.f  common.h phact_data_types.inc
	$(F77)  $(FLAGgg3z)  $*.f
twoztwowqcd.o  :  twoztwowqcd.f common.h
	$(F77)  $(FLAG2z2wqcd)  $*.f
fourzqcd.o : fourzqcd.f  common.h
	$(F77)  $(FLAG4zqcd)  $*.f

unit_amp.o : unit_amp.f  common.h common_unitarization.h
	$(F77)  $(FLAG)  $*.f
unit.o : unit.f  common.h common_unitarization.h
	$(F77)  $(FLAG)  $*.f
pol.o :   pol.f
	$(F77)  $(FLAG)  $*.f

dpa_proj.o :   dpa_proj.f common.h
	$(F77)  $(FLAG)  $*.f
singleosp.o :   singleosp.f common.h
	$(F77)  $(FLAG)  $*.f

clean:
	$(RM) -f $(OBJS) phantom.exe


