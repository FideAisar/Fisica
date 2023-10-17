C*********************************************************************
 
C...UPINIT
C...Routine to be called by user to set up user-defined processes.
C...Code below only intended as example, without any claim of realism.
C...Especially it shows what info needs to be put in HEPRUP.
 
      SUBROUTINE UPINIT
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)

C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/

C...Extra commonblock to transfer run info.
      COMMON/PRIV/MODE,NLIM   
      SAVE/PRIV/

C PHANTOM commonblocks
*sandro 27/7
c      COMMON/fxn_main/s,flux,factor
      INCLUDE 'common.h'     
*sandro 27/7 end

C....Pythia commonblock - needed for setting PDF's; see below.
C      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C      SAVE /PYPARS/      

*sandro 27/7
c      ecoll=sqrt(s)
*sandro 27/7 end

C...Set incoming beams: LHC.
      IDBMUP(1)=2212
      IDBMUP(2)=2212
      EBMUP(1)=ecoll/2.d0
      EBMUP(2)=ecoll/2.d0

C...Set PDF's of incoming beams: CTEQ 5L.
C...Note that Pythia will not look at PDFGUP and PDFSUP.  
      PDFGUP(1)=4
      PDFSUP(1)=46
      PDFGUP(2)=PDFGUP(1)
      PDFSUP(2)=PDFSUP(1)
      
C...If you want Pythia to use PDFLIB, you have to set it by hand.
C...(You also have to ensure that the dummy routines
C...PDFSET, STRUCTM and STRUCTP in Pythia are not linked.)      
C      MSTP(52)=2
C      MSTP(51)=1000*PDFGUP(1)+PDFSUP(1)

C...Decide on weighting strategy: we are generating unweighted events.
      IDWTUP=3

C...Number of external processes. 
      NPRUP=1

C...Set up process.
C Cross section
      XSECUP(1)=1.D0
C Error on cross section
      XERRUP(1)=0.d0
C Maximum weight
      XMAXUP(1)=1.D0
C Communication code
      LPRUP(1)=661

      RETURN
      END
 
C*********************************************************************
 
C...UPEVNT
C...Sample routine to generate events of various kinds.
C...Not intended to be realistic, but rather to show in closed
C...and understandable form what such a routine might look like.
C...Especially it shows what info needs to be put in HEPEUP.
 
      SUBROUTINE UPEVNT
 
C...Double precision and integer declarations.
      IMPLICIT REAL*8(A-H, O-Z)
      IMPLICIT INTEGER(I-N)

C...Extra commonblock to transfer run info.
      COMMON/PRIV/MODE,NLIM   
      SAVE/PRIV/

C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500) 
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      REAL*8 XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP  
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)  
      SAVE /HEPEUP/ 

      RETURN
      END 
 
