************************************************************************
*                                                                      *
*                              PHANTOM                                 *
*                                                                      *
*                     PHAct New TOrino Montecarlo                      *
*                                                                      *
*     a Monte Carlo program for six-fermion final state processes      *
*                             at the LHC                               *
*                                                                      *
*                   based on PHACT matrix elements                     *
*                                and                                   *
*                  adaptive+multichannel integration                   *
*                                                                      *
*                              Authors:                                *
*               Alessandro Ballestrero, Aissa Belhouari,               *
*                   Giuseppe Bevilacqua, Ezio Maina                    *
*                                                                      *
*                            version 1.0                               *
*                                                                      *
************************************************************************
*                                                                      *
* This program, its steering file, explanations and examples,          *
* can be obtained upon request from A. Ballestrero and E. Maina        *
*                                                                      *
************************************************************************

      Program PHANTOM


      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT COMPLEX*16 (c)

      INCLUDE 'common.h'
*sandro 27/7
c      COMMON/collenergy/ecoll
c      COMMON/fxn_main/s,flux,factor
      COMMON/fxn_main/flux,factor
*sandro 27/7 end

C for the reading system
      CHARACTER*80 pluto
      COMMON/abread/pluto
C How many Pythia event get printed
      COMMON/HADOUTPUT/IPRDEB
      DATA IPRDEB/5/
*ONESHOT 
      COMMON/phaones/ionesh
*PDF's
      CHARACTER*200 PDFname
      COMMON/phpdfli/ PDFname

      DATA pb_conv/0.38937966d9/


**Ezio Per SGE
**      PRINT*,  '  name of the input file ?'
**      READ (*,'(a)') pluto       ! name of the file
      pluto='r.in'
**Ezio-end Per SGE
      PRINT*, '                          INPUT FROM FILE '//pluto

      CALL readinput
	
      IF(ihadronize.EQ.1)THEN
C...Initialize with external process.
        CALL PYINIT('USER',' ',' ',0D0)
      ENDIF

c giuseppe ILC 25/06/2007
      IF(i_coll.EQ.1 .OR. i_coll.EQ.2)THEN ! hadron colliders

C...Initialize PDF's
c       PDFname='/home/usr04/ballestr/PDFsets/cteq5l.LHgrid'
*        CALL InitPDFset(PDFname)
        CALL InitPDFsetByName(PDFname)
        CALL InitPDF(0)

      ELSEIF(i_coll.EQ.3)THEN   ! e+e- collider
c giuseppe 30/04/2008: LHAPDF must be initialized also for e+e- 
c collisions if alphasPDF has to be used
        CALL InitPDFset(PDFname)
        CALL InitPDF(0)
c end giuseppe 30/04/2008
        WRITE(*,*),'*************************************'
        WRITE(*,*),'*          e+e- collider            *'
        WRITE(*,*),'*************************************'
        WRITE(*,*),''
        IF(i_beamstrahlung.EQ.1)THEN
          WRITE(*,*),'BEAMSTRAHLUNG: ON'
          WRITE(*,*),'>>>>>> Beamstrahlung description: <<<<<<'
          WRITE(*,*),'CIRCE 1.0'
          WRITE(*,*),'Reference:'
          WRITE(*,*),'T. Ohl'
          WRITE(*,*),'hep-ph/9607454'
          WRITE(*,*),'>>>>>>                            <<<<<<'
          WRITE(*,*),''
        ELSE
          WRITE(*,*),'BEAMSTRAHLUNG: OFF'
          WRITE(*,*),''
        ENDIF
        IF(i_isr.EQ.1)THEN
          WRITE(*,*),'ISR: ON'
          WRITE(*,*),'>>>>>> ISR description: <<<<<<'
          WRITE(*,*),'Reference:'
          WRITE(*,*),'W. Beenakker et al.'
          WRITE(*,*),'hep-ph/9602351'
          WRITE(*,*),'>>>>>>                  <<<<<<'
          WRITE(*,*),''
        ELSE
          WRITE(*,*),'ISR: OFF'
          WRITE(*,*),''
        ENDIF
        WRITE(*,*),''

      ENDIF                     ! IF(i_coll.EQ.1 .OR. i_coll.EQ.2)THEN
c end giuseppe ILC 25/06/2007

      CALL coupling

      CALL procini
      
      s=ecoll*ecoll
      flux=0.5d0/s              ! flux factor
      factor=pb_conv      ! overall conversion in pb



* Routine Vegas parameter:

c giuseppe ILC 25/06/2007
      IF(i_coll.EQ.1)THEN       ! p-p collider
        ndim=15
      ELSEIF(i_coll.EQ.2)THEN   ! p-pbar collider
        ndim=15
      ELSEIF(i_coll.EQ.3)THEN   ! e+e- collider
        IF(i_isr.EQ.1)THEN
          ndim=15
        ELSEIF(i_isr.EQ.0)THEN
          ndim=13
        ELSE
          print*,'***ERROR: i_isr must be either 0 or 1'
          STOP
        ENDIF
      ELSE
        print*,'***ERROR: i_coll must be either 1, 2 or 3'
        STOP
      ENDIF                     ! IF(i_coll.EQ.1)THEN

* Initialization of CIRCE parameterization for e+e- beamstrahlung
c Parameterization available only for center of mass energy = 350, 500, 
c 800, 1000 and 1600 GeV
      IF(i_coll.EQ.3)THEN
        IF(i_beamstrahlung.NE.0 .AND. i_beamstrahlung.NE.1)THEN
          print*,'***ERROR: i_beamstrahlung must be either 0 or 1'
          STOP
        ENDIF
        IF(i_beamstrahlung.EQ.1)THEN
          call CIRCES(0.d0,0.d0,ecoll,2,1,1997 04 17,0)
        ENDIF                   !i_beamstrahlung
      ENDIF                     !i_coll
c end giuseppe ILC 25/06/2007


* Integration limits:
      DO i=1,ndim
        region(i)=0.d0
        region(ndim+i)=1.d0
      ENDDO !i

*sandro 25/07
c      CALL proc  ! proc  determines the variables for every phase space
                 ! that can be used and the variables for the process
                 ! at hand
                 ! If ioneshot=1 the latter are read from phavegas.dat
                 !  and only the former are determined in proc

c      if (ionesh.eq.0) then

      if (ionesh.eq.0) then

      CALL proc  ! proc  determines the variables for every phase space
                 ! that can be used and the variables for the process
                 ! at hand
                 ! If ioneshot=1 the latter are read from phavegas.dat
                 !  and only the former are determined in proc
*sandro 25/07 end 

        CALL procextraini
                
        CALL extrema
            
        CALL integration

      else
        
        CALL oneshot    ! the call to extrema 
*sandro 25/07 
                        ! and to proc   
*sandro 25/07 end       
                        ! is in oneshot,
                        !  once the process has been chosen

      endif

      STOP
      END
      

       









