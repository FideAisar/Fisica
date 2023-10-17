************************************************************************
* LHAFileInit.f
*
* Last update: Jun 25, 2006
************************************************************************

      SUBROUTINE LHAFileInit
 
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

C PHANTOM commonblocks
      INCLUDE 'common.h'     

C...Set incoming beams:
      IF(i_coll.eq.1)THEN       ! p-p collider
        IDBMUP(1)=2212
        IDBMUP(2)=2212
      ELSEIF(i_coll.eq.2)THEN   ! p-pbar collider
        IDBMUP(1)=2212
        IDBMUP(2)=-2212
c giuseppe ILC 25/06/2007
      ELSEIF(i_coll.eq.3)THEN   ! e+e- collider
        IDBMUP(1)=-11
        IDBMUP(2)=11
c end giuseppe ILC 25/06/2007
      ENDIF
      EBMUP(1)=ecoll/2.d0
      EBMUP(2)=ecoll/2.d0

c giuseppe ILC 25/06/2007
      IF(i_coll.eq.1 .OR. i_coll.eq.2)THEN ! hadron colliders
C...Set PDF's of incoming beams: CTEQ 5L.
        PDFGUP(1)=4
        PDFSUP(1)=46
        PDFGUP(2)=PDFGUP(1)
        PDFSUP(2)=PDFSUP(1)
      ELSEIF(i_coll.eq.3)THEN   ! e+e- collider
        PDFGUP(1)=-1
        PDFSUP(1)=-1
        PDFGUP(2)=PDFGUP(1)
        PDFSUP(2)=PDFSUP(1)
      ENDIF
c end giuseppe ILC 25/06/2007
      
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

      IF(iwrite_event.EQ.1)THEN
        WRITE(23,'(A)') '<init>'
        WRITE(23,'(2(I6,1x),2(E13.5,1x),6(I4,1x))')
     &              IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),
     &              PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),
     &              IDWTUP,NPRUP
        WRITE(23,'(3(E11.1,1x),I4)') 
     &              XSECUP(1),XERRUP(1),XMAXUP(1),LPRUP(1)
        WRITE(23,'(A)') '</init>'
      ENDIF

      RETURN
      END
 
