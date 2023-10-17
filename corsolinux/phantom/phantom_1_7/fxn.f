c
c Last update: Jul 18, 2017
c
***********************************************************************
*     STANDARD VEGAS INPUT:
*     PARAMETERS
*     mxnphs          : total number of possible phase space mappings
*     THROGH COMMON BLOCKS:
*     iphs_ind    : phase space to be used for generating momenta
*     ialfa(mxnphs)  : 0/1 inactive/active phase space
*     alfa(mxnphs)   : phase space weight in multichannel
*     idp(8)       : particle identities: ASSUMED OUTGOING
*     idp_inout(8) : idp_inout(i)=1(-1) for outgoing(incoming) particles
*     iorder(8,mxnphs)
*
***********************************************************************

      REAL*8 FUNCTION fxn(x,wgt)

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT DOUBLE COMPLEX (c)

*************************OVERALL COMMONS *****************************

      INCLUDE 'common.h'
      INCLUDE 'common_cut.h'
      INCLUDE 'common_subproc.h'
c giuseppe 07/07/2007
      DIMENSION pt5min(3),pt6min(3)
c end giuseppe 07/07/2007
      COMMON/fxn_extrema/pt5min,pt6min,shatmin
*sandro 27/7
c      COMMON/fxn_main/s,flux,factor
      COMMON/fxn_main/flux,factor
*sandro 27/7 end
*     INTEGRATION COMMONS
c giuseppe 21/05/2007
      DIMENSION ncall_therm(2)
c end giuseppe 21/05/2007
      COMMON/iinteg_readinput/ncall_therm,itmx_therm,ncall,itmx
      COMMON/rinteg_readinput/acc_therm,acc
c giuseppe 07/02/2007
      COMMON/ifxn_readinput/i_PDFscale
*1_7
      COMMON/rfxn_readinput/pdfconst
*1_7end
c end giuseppe 07/02/2007
c scale choice
      COMMON/rfxn_readinput/fixed_PDFscale
c end scale choice

      COMMON/phspinput_proc/
     &  rm0a(mxnphs),gam0a(mxnphs),spm0alo(mxnphs),spm0ahi(mxnphs),
     &  rm0b(mxnphs),gam0b(mxnphs),spm0blo(mxnphs),spm0bhi(mxnphs),
     &  rm1(mxnphs),gam1(mxnphs),spm1lo(mxnphs),spm1hi(mxnphs),
     &  rm1a(mxnphs),gam1a(mxnphs),spm1alo(mxnphs),spm1ahi(mxnphs),
     &  rm1b(mxnphs),gam1b(mxnphs),spm1blo(mxnphs),spm1bhi(mxnphs),
     &  rm11(mxnphs),gam11(mxnphs),spm11lo(mxnphs),spm11hi(mxnphs),
     &  rm11a(mxnphs),gam11a(mxnphs),spm11alo(mxnphs),spm11ahi(mxnphs),
     &  rm11b(mxnphs),gam11b(mxnphs),spm11blo(mxnphs),spm11bhi(mxnphs),
     &  rm11c(mxnphs),gam11c(mxnphs),spm11clo(mxnphs),spm11chi(mxnphs),
     &  rm12a(mxnphs),gam12a(mxnphs),spm12alo(mxnphs),spm12ahi(mxnphs),
     &  rm12b(mxnphs),gam12b(mxnphs),spm12blo(mxnphs),spm12bhi(mxnphs),
     &  rm12c(mxnphs),gam12c(mxnphs),spm12clo(mxnphs),spm12chi(mxnphs),
     &  rm2(mxnphs),gam2(mxnphs),spm2lo(mxnphs),spm2hi(mxnphs),
     &  rm2a(mxnphs),gam2a(mxnphs),spm2alo(mxnphs),spm2ahi(mxnphs),
     &  rm2b(mxnphs),gam2b(mxnphs),spm2blo(mxnphs),spm2bhi(mxnphs),
     &  rm2c(mxnphs),gam2c(mxnphs),spm2clo(mxnphs),spm2chi(mxnphs),
     &  rm21a(mxnphs),gam21a(mxnphs),spm21alo(mxnphs),spm21ahi(mxnphs),
     &  rm21b(mxnphs),gam21b(mxnphs),spm21blo(mxnphs),spm21bhi(mxnphs),
     &  rm111(mxnphs),gam111(mxnphs),spm111lo(mxnphs),spm111hi(mxnphs),
     &  rm111a(mxnphs),gam111a(mxnphs),spm111alo(mxnphs),
     &     spm111ahi(mxnphs),
     &  rm111b(mxnphs),gam111b(mxnphs),spm111blo(mxnphs),
     &     spm111bhi(mxnphs),
     &  rm1111(mxnphs),gam1111(mxnphs),spm1111lo(mxnphs),
     &     spm1111hi(mxnphs),
     &  rmta(mxnphs),expa(mxnphs),spta,rmtb(mxnphs),expb(mxnphs),sptb

      COMMON /pharand/ idum

***********************VEGAS *****************************************
      DIMENSION x(15)
*6
c      COMMON/abresl/resl(10),standdevl(10)
*6end
      COMMON/abstat/ncall_eff
******************CONTROL FLAGS***************************************
      COMMON/i_norm/i_normalize
      COMMON/phaones/ionesh
*************************EXTERNAL PARTICLES **************************
      DIMENSION  xm(8),xmsq(8)
*************************DISTRIBUTION FUNCTIONS **********************
      DIMENSION distfun(8)
********LESHOUCHES USER PROCESS EVENT COMMON BLOCK *******************
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      REAL*8 XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &   ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &   VTIMUP(MAXNUP),SPINUP(MAXNUP)
********PDF's in LESHOUCHES EVENT FILE FORMAT *******************
      INTEGER id1_LH(8),id2_LH(8)
      REAL*8 xpdf1_LH(8),xpdf2_LH(8)
********LHAPDF *******************
      REAL*8 ff1(-6:6),ff2(-6:6)
c giuseppe ILC 26/05/2007
************BEAMSTRAHLUNG AND ISR (for e+e- initial state) ************
      REAL*8 ffISR, random
      EXTERNAL ffISR
      EXTERNAL random
c end giuseppe ILC 26/05/2007
c scale choice
      DIMENSION nlept(4),nquark(4),qucos(4),nqucentr(2)
c

*6
c******************************PYTHIA *********************************
c      PARAMETER(NMXHEP=4000)
c      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
c     &   JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
*6end
***********************Hadronization output **************************
      INTEGER IPRDEB
      COMMON/HADOUTPUT/IPRDEB
*******************PHASE SPACE and MULTICHANNEL **********************
*lab
      DIMENSION p(0:3,8),pcoll(0:3,8),plab(0:3,8)

*4cmpol
      DIMENSION p4cm(0:3),pcall(0:3,8)
*4cmpolend

*     res
c     giovanni 24/4/18 begin
      REAL*8 offactor
c     include the common with idosp(4)
c     include the common with idp(8)
c     giovanni 24/4/18 end
*zw
*      common/osp/i_osp
*      common/polww/i_ww,i_pol,i_polwp,i_polwm,i_polwel,i_polwmu     
*zwend
*resend

 
*labend
c giuseppe 29/09/2007
      INTEGER iphsp
c end giuseppe 29/09/2007
*************************PERMUTATIONS ********************************
      INTEGER idin0(2),ipos0(3,-16:21),index(4,24),ifact(0:4),
     &     nidentical(-16:21),idin0_pdf(2),
     &     iordaux(8),iorderaux(8),iordtemp(3),nmultiple,
     &     i_multzero,i_identical(3),index_pythia(8),iordamp(8),
     &     idp0(8),idp0_inout(8)
      REAL*8 paux(0:3,8),symfact
      DATA index/
     &           1,2,3,4,         !even
     &           2,1,3,4,
     &           3,1,2,4,         !even
     &           1,3,2,4,
     &           2,3,1,4,         !even
     &           3,2,1,4,         
     &  !permutations of first 3 indices end here
     &           4,2,1,3,         !even
     &           2,4,1,3,     
     &           1,4,2,3,         !even
     &           4,1,2,3,
     &           2,1,4,3,         !even
     &           1,2,4,3,
     &
     &           3,4,1,2,         !even
     &           4,3,1,2,
     &           4,1,3,2,         !even
     &           1,4,3,2,
     &           1,3,4,2,         !even
     &           3,1,4,2,
     &
     &           4,3,2,1,         !even
     &           3,4,2,1,
     &           2,4,3,1,         !even
     &           4,2,3,1,
     &           4,2,4,1,         !even
     &           2,3,4,1
     &          /,
     &     ifact/1,1,2,6,24/
***   SYMMETRY BETWEEN FIRST AND SECOND FAMILY, C CONJUGATION **********
      INTEGER ifamexchange12(-16:21)
      DATA    ifamexchange12/
     &     -16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-2,-1,-4,-3,0,
     &     3,4,1,2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21/
****************************CUTS *************************************
      INTEGER ipasscuts
      EXTERNAL ipasscuts
**********************RANDOM NUMBERS *********************************
      REAL*8 ran2
      INTEGER idum
*************************UTILITY *************************************
      INTEGER i,j,l,k,i0,i1,i2,i3

*************************UTILITY *************************************
      DATA ifirst / 1/

***********************************************************************
*************************EXECUTABLE CODE *****************************
***********************************************************************

c giuseppe ILC 26/06/2007
* INITIAL-STATE KINEMATICS
      IF(i_coll.EQ.1 .OR. i_coll.EQ.2)THEN ! hadron colliders

        i_startout=3

        e_cm_min=sqrt(shatmin)
        e_cm_max=sqrt(s)
        spm1=3.d3
        call PHSP_INI2(e_cm_max,e_cm_min,spm1,
     &                           alfa1,
     &               xf1_min,xf1_max,xf2_min,xf2_max)
        IF(x(1).LE.alfa1.AND.alfa1.GT.0.d0)THEN
          vegas_x=x(1)/alfa1
          e_cm=(xf1_max-xf1_min)*vegas_x+xf1_min
          rjacz=(xf1_max-xf1_min)*2.d0*e_cm/alfa1
        ELSE
          vegas_x=(x(1)-alfa1)/(1.d0-alfa1)
          e_cm=(xf2_max-xf2_min)*vegas_x+xf2_min
          rjacz=(xf2_max-xf2_min)*2.d0*e_cm/(1.d0-alfa1)
        ENDIF
        z=e_cm*e_cm/s
        rjacz=rjacz/s           ! Jacobian of dz
        IF(z.LE.0.d0.OR.z.GT.1.d0)THEN
          fxn=0.d0
          RETURN
        ENDIF
        ys=-2.d0*log(sqrt(z))*x(2)+log(sqrt(z))
        rjacys=-2.d0*log(sqrt(z)) ! Jacobian of dy

      ELSEIF(i_coll.EQ.3)THEN   ! e+e- collider

        i_startout=1

        IF(i_beamstrahlung.EQ.1)THEN
* BEAMSTRAHLUNG on: x1beam, x2beam generated unweighted by GIRCEE
          call GIRCEE(x1beam,x2beam,random)
          z1=x1beam*x2beam
c Note: with current definition of b1 and b2, ys is MINUS the rapidity
          ys1=-0.5d0*log(x1beam/x2beam)
        ELSE
* BEAMSTRAHLUNG off
          z1=1.d0
          ys1=0.d0
        ENDIF                   ! IF(i_beamstrahlung.EQ.1)THEN

        IF(i_isr.EQ.1)THEN
* ISR on: the number of integration variables is 15
          i_startout=i_startout+2

          saux=z1*s

          call map_ISR(x(1),x(2),saux,bb1,bb2)
          z2=bb1*bb2
          ys2=-0.5d0*log(bb1/bb2)
          rjacz=1.d0
          rjacys=1.d0

        ELSE
* ISR off
          z2=1.d0
          ys2=0.d0
          rjacz=1.d0
          rjacys=1.d0
        ENDIF                   ! IF(i_isr.EQ.1)THEN

        z=z1*z2
        ys=ys1+ys2

      ENDIF                     ! IF(i_coll.EQ.1 .OR. i_coll.EQ.2)THEN
c end giuseppe ILC 26/06/2007


c giuseppe 19/12/2006: iflip is a new flag for the exchange x1<->x2.
c Each exchange of initial state momenta makes the flag to change sign
      iflip=1
c end giuseppe 19/12/2006

      fxn=1.d0

      shat=z*s
      IF(shat.LE.shatmin)THEN
        fxn=0.d0
        RETURN
      ENDIF
      ecm=dsqrt(shat)
      fluxhat=flux/z

* SOME USEFUL QUANTITIES FOR BOOSTS FROM HARD SCATTERING REST FRAME
*     TO COLLIDER FRAME
      b1=sqrt(z)*exp(-ys)
      b2=sqrt(z)*exp(ys)
c giuseppe ILC 26/06/2007
c      IF(b1.LE.1.d-5 .OR .b2.LE.1.d-5 .OR. b1.GE.1.d0
c     &     .OR. b2.GE.1.d0)THEN
      IF(b1.LE.1.d-5 .OR. b2.LE.1.d-5 .OR. b1.GT.1.d0
     &     .OR. b2.GT.1.d0)THEN
c end giuseppe ILC 26/06/2007
        fxn=0.d0
        RETURN
      ENDIF
*     boost variables
      IF((b1+b2).LE.0.d0 .OR. (b1*b2).LE.0.d0)THEN
        fxn=0.d0
        RETURN
      ENDIF
      bcm=(b1-b2)/(b1+b2)
      gcm=(b1+b2)/(2.d0*sqrt(b1*b2))

*     MASSES
      DO i=1,8
        xm(i)=rmass(idp(i))
        xmsq(i)=rmass2(idp(i))
      ENDDO

*     EXTRACT INITIAL STATE PARTICLES
      idin0(1)=idp(iorder(1,iphs_ind))
      idin0(2)=idp(iorder(2,iphs_ind))
      IF(idin0(1).EQ.idin0(2))THEN
        nidenticalinitial=2
      ELSE
        nidenticalinitial=1
      ENDIF

      xmasq=xmsq(iorder(1,iphs_ind))
      xmbsq=xmsq(iorder(2,iphs_ind))

c giuseppe ILC 26/06/2007
      IF(iphs_ind.le.13)THEN
c giuseppe 29/09/2007
        iphsp=1
c end giuseppe 29/09/2007
        CALL phsp1_1_4_multi7_c(x(i_startout),shat,
     &       rm1a(iphs_ind),gam1a(iphs_ind),spm1alo(iphs_ind),
     &       spm1ahi(iphs_ind),
     &       rm1b(iphs_ind),gam1b(iphs_ind),spm1blo(iphs_ind),
     &       spm1bhi(iphs_ind),
     &       rm11a(iphs_ind),gam11a(iphs_ind),
     &       spm11alo(iphs_ind),spm11ahi(iphs_ind),
     &       rm11b(iphs_ind),gam11b(iphs_ind),spm11blo(iphs_ind),
     &       spm11bhi(iphs_ind),
     &       rm11c(iphs_ind),gam11c(iphs_ind),spm11clo(iphs_ind),
     &       spm11chi(iphs_ind),
     &       rm12a(iphs_ind),gam12a(iphs_ind),
     &       spm12alo(iphs_ind),spm12ahi(iphs_ind),
     &       rm12b(iphs_ind),gam12b(iphs_ind),spm12blo(iphs_ind),
     &       spm12bhi(iphs_ind),
     &       rm12c(iphs_ind),gam12c(iphs_ind),spm12clo(iphs_ind),
     &       spm12chi(iphs_ind),
     &       rmta(iphs_ind),expa(iphs_ind),
     &       spta,rmtb(iphs_ind),expb(iphs_ind),sptb,
     &       pt5min,pt6min,
c giuseppe 07/07/2007
     &       idp(iorder(7,iphs_ind)),idp(iorder(8,iphs_ind)),
c end giuseppe 07/07/2007
     &       xmasq,xmbsq,
     &       xm,xmsq,iorder(3,iphs_ind),
     &       p,xjac)
      ELSE IF(iphs_ind.le.24)THEN
c giuseppe 29/09/2007
        iphsp=2
c end giuseppe 29/09/2007
        CALL phsp2_4_multi7(x(i_startout),shat,
     &       rm1a(iphs_ind),gam1a(iphs_ind),spm1alo(iphs_ind),
     &       spm1ahi(iphs_ind),
     &       rm1b(iphs_ind),gam1b(iphs_ind),spm1blo(iphs_ind),
     &       spm1bhi(iphs_ind),
     &       rm11a(iphs_ind),gam11a(iphs_ind),spm11alo(iphs_ind),
     &       spm11ahi(iphs_ind),
     &       rm11b(iphs_ind),gam11b(iphs_ind),spm11blo(iphs_ind),
     &       spm11bhi(iphs_ind),
     &       rm11c(iphs_ind),gam11c(iphs_ind),spm11clo(iphs_ind),
     &       spm11chi(iphs_ind),
     &       rm12a(iphs_ind),gam12a(iphs_ind),
     &       spm12alo(iphs_ind),spm12ahi(iphs_ind),
     &       rm12b(iphs_ind),gam12b(iphs_ind),spm12blo(iphs_ind),
     &       spm12bhi(iphs_ind),
     &       rm12c(iphs_ind),gam12c(iphs_ind),spm12clo(iphs_ind),
     &       spm12chi(iphs_ind),
     &       rm2a(iphs_ind),gam2a(iphs_ind),
     &       spm2alo(iphs_ind),spm2ahi(iphs_ind),
     &       rm2b(iphs_ind),gam2b(iphs_ind),
     &       spm2blo(iphs_ind),spm2bhi(iphs_ind),
     &       rm2c(iphs_ind),gam2c(iphs_ind),
     &       spm2clo(iphs_ind),spm2chi(iphs_ind),
     &       xm,xmsq,iorder(3,iphs_ind),
     &       p,xjac)
      ELSE IF(iphs_ind.le.32)THEN
c giuseppe 29/09/2007
        iphsp=3
c end giuseppe 29/09/2007
        CALL phsp1_1_31_multi_c(x(i_startout),shat,
     &       rm1a(iphs_ind),gam1a(iphs_ind),spm1alo(iphs_ind),
     &        spm1ahi(iphs_ind),
     &       rm1b(iphs_ind),gam1b(iphs_ind),spm1blo(iphs_ind),
     &        spm1bhi(iphs_ind),
     &       rm11a(iphs_ind),gam11a(iphs_ind),spm11alo(iphs_ind),
     &        spm11ahi(iphs_ind),
     &       rm11b(iphs_ind),gam11b(iphs_ind),spm11blo(iphs_ind),
     &        spm11bhi(iphs_ind),
     &       rm111a(iphs_ind),gam111a(iphs_ind),spm111alo(iphs_ind),
     &        spm111ahi(iphs_ind),
     &       rm111b(iphs_ind),gam111b(iphs_ind),spm111blo(iphs_ind),
     &        spm111bhi(iphs_ind),
     &       rmta(iphs_ind),expa(iphs_ind),spta,
     &       rmtb(iphs_ind),expb(iphs_ind),sptb,
     &       pt5min,pt6min,
c giuseppe 07/07/2007
     &       idp(iorder(7,iphs_ind)),idp(iorder(8,iphs_ind)),
c end giuseppe 07/07/2007
     &       xmasq,xmbsq,
     &       xm,xmsq,iorder(3,iphs_ind),
     &       p,xjac)
      ELSE IF(iphs_ind.le.54)THEN
c giuseppe 29/09/2007
        iphsp=4
c end giuseppe 29/09/2007
        CALL phsp1_2_3_multi5_c(x(i_startout),shat,
     &       rm0a(iphs_ind),gam0a(iphs_ind),spm0alo(iphs_ind),
     &        spm0ahi(iphs_ind),
     &       rm0b(iphs_ind),gam0b(iphs_ind),spm0blo(iphs_ind),
     &        spm0bhi(iphs_ind),
     &       rm1a(iphs_ind),gam1a(iphs_ind),spm1alo(iphs_ind),
     &        spm1ahi(iphs_ind),
     &       rm1b(iphs_ind),gam1b(iphs_ind),spm1blo(iphs_ind),
     &        spm1bhi(iphs_ind),
     &       rm11a(iphs_ind),gam11a(iphs_ind),spm11alo(iphs_ind),
     &        spm11ahi(iphs_ind),
     &       rm11b(iphs_ind),gam11b(iphs_ind),spm11blo(iphs_ind),
     &        spm11bhi(iphs_ind),
     &       rm2a(iphs_ind),gam2a(iphs_ind),spm2alo(iphs_ind),
     &        spm2ahi(iphs_ind),
     &       rm2b(iphs_ind),gam2b(iphs_ind),spm2blo(iphs_ind),
     &        spm2bhi(iphs_ind),
     &       rmta(iphs_ind),expa(iphs_ind),spta,
     &       pt6min,
     &       idin0(1),idin0(2),idp(iorder(8,iphs_ind)),
     &       xmasq,xmbsq,
     &       xm,xmsq,iorder(3,iphs_ind),
     &       p,xjac)
      ELSE IF(iphs_ind.le.60)THEN
c giuseppe 29/09/2007
        iphsp=5
c end giuseppe 29/09/2007
        CALL phsp3_3_multi5(x(i_startout),shat,
     &       rm1a(iphs_ind),gam1a(iphs_ind),spm1alo(iphs_ind),
     &        spm1ahi(iphs_ind),
     &       rm1b(iphs_ind),gam1b(iphs_ind),spm1blo(iphs_ind),
     &        spm1bhi(iphs_ind),
     &       rm11a(iphs_ind),gam11a(iphs_ind),spm11alo(iphs_ind),
     &        spm11ahi(iphs_ind),
     &       rm11b(iphs_ind),gam11b(iphs_ind),spm11blo(iphs_ind),
     &        spm11bhi(iphs_ind),
     &       rm2a(iphs_ind),gam2a(iphs_ind),spm2alo(iphs_ind),
     &        spm2ahi(iphs_ind),
     &       rm2b(iphs_ind),gam2b(iphs_ind),spm2blo(iphs_ind),
     &        spm2bhi(iphs_ind),
     &       rm21a(iphs_ind),gam21a(iphs_ind),spm21alo(iphs_ind),
     &        spm21ahi(iphs_ind),
     &       rm21b(iphs_ind),gam21b(iphs_ind),spm21blo(iphs_ind),
     &        spm21bhi(iphs_ind),
     &       xm,xmsq,iorder(3,iphs_ind),
     &       p,xjac)
      ELSE IF(iphs_ind.le.68)THEN
c giuseppe 29/09/2007
        iphsp=6
c end giuseppe 29/09/2007
        CALL phsp2_4to31_multi5(x(i_startout),shat,
     &       rm1a(iphs_ind),gam1a(iphs_ind),spm1alo(iphs_ind),
     &        spm1ahi(iphs_ind),
     &       rm1b(iphs_ind),gam1b(iphs_ind),spm1blo(iphs_ind),
     &        spm1bhi(iphs_ind),
     &       rm11a(iphs_ind),gam11a(iphs_ind),spm11alo(iphs_ind),
     &        spm11ahi(iphs_ind),
     &       rm11b(iphs_ind),gam11b(iphs_ind),spm11blo(iphs_ind),
     &        spm11bhi(iphs_ind),
     &       rm111a(iphs_ind),gam111a(iphs_ind),spm111alo(iphs_ind),
     &        spm111ahi(iphs_ind),
     &       rm111b(iphs_ind),gam111b(iphs_ind),spm111blo(iphs_ind),
     &        spm111bhi(iphs_ind),
     &       rm2a(iphs_ind),gam2a(iphs_ind),spm2alo(iphs_ind),
     &        spm2ahi(iphs_ind),
     &       rm2b(iphs_ind),gam2b(iphs_ind),spm2blo(iphs_ind),
     &        spm2bhi(iphs_ind),
     &       xm,xmsq,iorder(3,iphs_ind),
     &       p,xjac)
      ELSE IF(iphs_ind.le.76)THEN
c giuseppe 29/09/2007
        iphsp=7
c end giuseppe 29/09/2007
        CALL phsp1_5to1_4to31_multi_c(x(i_startout),shat,
     &       rm1(iphs_ind),gam1(iphs_ind),spm1lo(iphs_ind),
     &        spm1hi(iphs_ind),
     &       rm11(iphs_ind),gam11(iphs_ind),spm11lo(iphs_ind),
     &        spm11hi(iphs_ind),
     &       rm111(iphs_ind),gam111(iphs_ind),spm111lo(iphs_ind),
     &        spm111hi(iphs_ind),
     &       rm1111(iphs_ind),gam1111(iphs_ind),spm1111lo(iphs_ind),
     &        spm1111hi(iphs_ind),
     &       rmta(iphs_ind),expa(iphs_ind),spta,
     &       pt6min,
     &       idin0(1),idin0(2),idp(iorder(8,iphs_ind)),
     &       xmasq,xmbsq,
     &       xm,xmsq,iorder(3,iphs_ind),
     &       p,xjac)
      ELSE
        print*,'***fxn.f ERROR'
        print*,'   iphs_ind:',iphs_ind,
     &       'exceeding the maximum number of channels'
        STOP
      ENDIF
c end giuseppe ILC 26/06/2007

      IF(xjac.LE.0.d0)THEN
        fxn=0.d0
        RETURN
      ENDIF

* INCOMING MOMENTA IN HARD SCATTERING CENTER OF MASS FRAME
*  THEY ARE IN A FIXED POSITION IN THE ARRAY p EVEN WHEN IDENTICAL
*  PARTICLES ARE PRESENT
      app=((shat-xmasq-xmbsq)**2
     &     -4.d0*xmasq*xmbsq)/(4.d0*shat)
      IF(app.LE.0.d0)THEN
        fxn=0.d0
        RETURN
      ENDIF
      pmod=sqrt(app)
      p(0,iorder(1,iphs_ind))=sqrt(pmod*pmod+xmasq)
      p(1,iorder(1,iphs_ind))=0.d0
      p(2,iorder(1,iphs_ind))=0.d0
      p(3,iorder(1,iphs_ind))=pmod
      p(0,iorder(2,iphs_ind))=sqrt(pmod*pmod+xmbsq)
      p(1,iorder(2,iphs_ind))=0.d0
      p(2,iorder(2,iphs_ind))=0.d0
      p(3,iorder(2,iphs_ind))=-pmod

* BOOST TO COLLIDER FRAME
      DO i=1,8
        CALL ZBOOST(gcm,bcm,p(0,i),pcoll(0,i))
      ENDDO

* CUTS: EXPECTS pcoll((0:3,8) in one-to-one correspondence with idp and
*     idp_inout


      ipass1=ipasscuts(pcoll,
     &     e_min_lep,pt_min_lep,
c giuseppe 01/03/2007
     &     eta_max_onelep,
c end giuseppe 01/03/2007
     &     eta_max_lep,ptmiss_min,
     &     e_min_j,pt_min_j,eta_max_j,
     &     eta_def_jf_min,eta_def_jb_max,eta_def_jc_max,
     &     pt_min_jcjc,rm_min_jj,rm_min_ll,rm_min_jlep,
     &     rm_min_jcjc,rm_max_jcjc,rm_min_jfjb,
* six
     &     rm_min_4l,rm_min_2l2cq,
* YR
     &     rm_max_4l,rm_max_2l2cq,
* YRend
* sixend
* cuttop
     &     deltacuttop, 
* cuttopend
*ww
     &     rcutwlept,rcutptwelectron,
*wwend
*zw
     &     rcutzlept,
*zwend

     &     eta_min_jfjb,
     &     d_ar_jj,d_ar_jlep,d_ar_leplep,
     &     thetamin_jj,thetamin_jlep,
c giuseppe 01/03/2007
     &     thetamin_leplep,
c end giuseppe 01/03/2007
     &     i_e_min_lep,i_pt_min_lep,
     &     i_eta_max_onelep,
     &     i_eta_max_lep,i_ptmiss_min,
     &     i_e_min_j,i_pt_min_j,i_eta_max_j,i_eta_jf_jb_jc,
     &     i_pt_min_jcjc,i_rm_min_jj,i_rm_min_ll,i_rm_min_jlep,
     &     i_rm_min_jcjc,i_rm_max_jcjc,i_rm_min_jfjb,
* six
     &     i_rm_min_4l,i_rm_min_2l2cq,
* YR
     &     i_rm_max_4l,i_rm_max_2l2cq,
* YRend
* sixend
* cuttop
     &     i_deltacuttop, 
* cuttopend
*ww
     &     i_cutwlept, i_cutptwelectron,
*wwend
*zw
     &     i_cutzlept, 
*zwend
     &     i_eta_min_jfjb,
     &     i_d_ar_jj,i_d_ar_jlep,i_d_ar_leplep,
     &     i_thetamin_jj,i_thetamin_jlep,
c giuseppe 01/03/2007
     &     i_thetamin_leplep,
c end giuseppe 01/03/2007
     &     i_usercuts, 0)

      IF (ipass1.eq.0)then
        fxn=0.d0
        RETURN
      ENDIF


      if (ionesh.eq.1.and.iextracuts.eq.1) then
        ipass2=ipasscuts(pcoll,
     &      e_min_lepos,pt_min_lepos,
c giuseppe 01/03/2007
     &      eta_max_onelepos,
c end giuseppe 01/03/2007
     &      eta_max_lepos,ptmiss_minos,
     &      e_min_jos,pt_min_jos,eta_max_jos,
     &      eta_def_jf_minos,eta_def_jb_maxos,eta_def_jc_maxos,
     &      pt_min_jcjcos,rm_min_jjos,rm_min_llos,rm_min_jlepos,
     &      rm_min_jcjcos,rm_max_jcjcos,rm_min_jfjbos,
* six
     &      rm_min_4los,rm_min_2l2cqos,
* sixend
* cuttop
     &     deltacuttopos, 
* cuttopend
     &      eta_min_jfjbos,
     &      d_ar_jjos,d_ar_jlepos,d_ar_leplepos,
     &      thetamin_jjos,thetamin_jlepos,
c giuseppe 01/03/2007
     &      thetamin_leplepos,
c end giuseppe 01/03/2007
     &      i_e_min_lepos,i_pt_min_lepos,
     &      i_eta_max_onelepos,
     &      i_eta_max_lepos,i_ptmiss_minos,
     &      i_e_min_jos,i_pt_min_jos,i_eta_max_jos,i_eta_jf_jb_jcos,
     &      i_pt_min_jcjcos,i_rm_min_jjos,i_rm_min_llos,i_rm_min_jlepos,
     &      i_rm_min_jcjcos,i_rm_max_jcjcos,i_rm_min_jfjbos,
* six
     &      i_rm_min_4los,i_rm_min_2l2cqos,
* sixend
* cuttop
     &     i_deltacuttopos, 
* cuttopend
     &      i_eta_min_jfjbos,
     &      i_d_ar_jjos,i_d_ar_jlepos,i_d_ar_leplepos,
     &      i_thetamin_jjos,i_thetamin_jlepos,
c giuseppe 01/03/2007
     &      i_thetamin_leplepos,
c end giuseppe 01/03/2007
     &      0,i_usercutsos)

        IF (ipass2.eq.0)then
          fxn=0.d0
          RETURN
        ENDIF
      endif

* CALL DISTRIBUTION FUNCTIONS

*  FACTORIZATION SCALE FOR PDF'S:

c giuseppe 01/03/2007
* The input flag i_PDFscale allows the user to choose the way of
* calculating the PDF scale.
*  i_PDFscale = 1: for all processes, the pT's of ALL outgoing particles
*                  is considered, so that Q**2= rmw**2+average(pT**2);
*  i_PDFscale = 2: process by process, the flavour content is analysed
*                  and the pT of the RECONSTRUCTED TOP(ANTITOP) quark is
*                  considered, so that Q**2= rmt**2+pT_top**2
*                  If this is impossible due to the flavour content, the
*                  PDF scale is calculated as in the case i_PDFscale=1.
*  i_PDFscale = 3: fixed scale given in input
*
*  i_PDFscale = 4  Q= m4l/sqrt(2) (invariant mass of the 4 leptons)/sqrt(2)
*                  (valid only for prosesses with four outgoing leptons)
*
c giovanni 2017/10/24 begin
*  i_PDFscale = 5  Q= sqrt(ptj1*ptj2) square-root of the product of
*                  transverse momenta of the 2 jets with largest pt
*                  (for all processes with at least 2 final state jets)
c giovanni 2017/10/24 end


      i_iftop=0                 ! top production flag is OFF by default

      IF(i_PDFscale.EQ.2)THEN

        n_igr=0
        IF(igr(2).EQ.1)THEN
          n_igr=2
        ELSEIF(igr(5).EQ.1)THEN
          n_igr=5
        ENDIF

        IF(n_igr.NE.0)THEN      ! 2Z2W or 2g1Z2W

          n_bout=0
          n_bbarout=0
          DO i=1,8
            IF(idp_inout(i).GT.0)THEN
              IF(idp(i).EQ.5)THEN
                n_bout=n_bout+1
              ELSEIF(idp(i).EQ.-5)THEN
                n_bbarout=n_bbarout+1
              ENDIF
            ENDIF
          ENDDO

          pT2_top=0.d0
          rm2ref=0.d0

          IF(n_bout.GT.0
     &         .AND. idp_inout(igo(5,n_igr)).GT.0
     &         .AND. idp_inout(igo(6,n_igr)).GT.0)THEN

            i_iftop=1           ! set top production flag ON
            px_wplus=p(1,igo(5,n_igr))+p(1,igo(6,n_igr))
            py_wplus=p(2,igo(5,n_igr))+p(2,igo(6,n_igr))
            pz_wplus=p(3,igo(5,n_igr))+p(3,igo(6,n_igr))
            E_wplus=p(0,igo(5,n_igr))+p(0,igo(6,n_igr))

            ncounter=0
            DO i=1,3,2
              IF(idp(igo(i,n_igr)).EQ.5
     &             .AND. idp_inout(igo(i,n_igr)).GT.0)THEN
                ncounter=ncounter+1
                px_bout=p(1,igo(i,n_igr))
                py_bout=p(2,igo(i,n_igr))
                pz_bout=p(3,igo(i,n_igr))
                E_bout=p(0,igo(i,n_igr))

                pt2aux=(px_bout+px_wplus)**2
     &               +(py_bout+py_wplus)**2
                rm2aux=(E_bout+E_wplus)**2-pt2aux
     &               -(pz_bout+pz_wplus)**2

                IF(ncounter.EQ.1)THEN
                  rm2ref=rm2aux
                  pT2_top=pt2aux
                ELSE ! solve the ambiguity on top pT reconstruction
                  IF(abs(rm2aux-rmt2).LT.abs(rm2ref-rmt2))THEN
                    pT2_top=pt2aux
                  ENDIF
                ENDIF
              ENDIF
            ENDDO

          ELSEIF(n_bbarout.GT.0
     &           .AND. idp_inout(igo(7,n_igr)).GT.0
     &           .AND. idp_inout(igo(8,n_igr)).GT.0)THEN

            i_iftop=1           ! set top production flag ON
            px_wminus=p(1,igo(7,n_igr))+p(1,igo(8,n_igr))
            py_wminus=p(2,igo(7,n_igr))+p(2,igo(8,n_igr))
            pz_wminus=p(3,igo(7,n_igr))+p(3,igo(8,n_igr))
            E_wminus=p(0,igo(7,n_igr))+p(0,igo(8,n_igr))

            ncounter=0
            DO i=2,4,2
              IF(idp(igo(i,n_igr)).EQ.-5
     &             .AND. idp_inout(igo(i,n_igr)).GT.0)THEN
                ncounter=ncounter+1
                px_bbarout=p(1,igo(i,n_igr))
                py_bbarout=p(2,igo(i,n_igr))
                pz_bbarout=p(3,igo(i,n_igr))
                E_bbarout=p(0,igo(i,n_igr))

                pt2aux=(px_bbarout+px_wminus)**2
     &               +(py_bbarout+py_wminus)**2
                rm2aux=(E_bbarout+E_wminus)**2-pt2aux
     &               -(pz_bbarout+pz_wminus)**2

                IF(ncounter.EQ.1)THEN
                  rm2ref=rm2aux
                  pT2_top=pt2aux
                ELSE ! solve the ambiguity on top pT reconstruction
                  IF(abs(rm2aux-rmt2).LT.abs(rm2ref-rmt2))THEN
                    pT2_top=pt2aux
                  ENDIF
                ENDIF
              ENDIF
            ENDDO

          ENDIF                 ! IF(n_bout.GT.0 ...
        ENDIF                   ! IF(n_igr.NE.0)THEN

        IF(i_iftop.EQ.1)THEN    !top production => use top pT scale
          Q=rmt**2+pT2_top
ctest giuseppe
          if(pT2_top.EQ.0.d0)then
            print*,'***fxn.f: ERROR in PDF scale calculation !!!'
            print*,'n_bout,n_bbarout',n_bout,n_bbarout
            print*,'pt2aux,pT2_top',pt2aux,pT2_top
            print*,'rm2aux,rm2ref,rmt2',rm2aux,rm2ref,rmt2
            print*,'n_igr',n_igr
            print*,'idp(igo(i,n_igr))',(idp(igo(i,n_igr)),i=1,8)
            STOP
          endif
c          if (ionesh.ne.1 .and. ifirst.eq.1) then
c            print*,' '
c            print*, 'PDF scale: Q**2= m_top**2+pT_top**2'
c            print*,' '
c          endif
ctestend giuseppe
          PDFscale=sqrt(Q)
        ENDIF
      ENDIF                     ! IF(i_PDFscale.EQ.2)THEN


      IF(i_PDFscale.EQ.1 .OR.(i_PDFscale.EQ.2.and.i_iftop.EQ.0))THEN
        Q=rmw**2
        DO i=1,8
          IF(idp_inout(i).GT.0)THEN
            Q=Q+(p(1,i)**2+p(2,i)**2)/6.d0
          ENDIF
        ENDDO
ctest giuseppe
c        if (ionesh.ne.1 .and. ifirst.eq.1) then
c          print*,' '
c          print*, 'PDF scale: Q**2= m_W**2+average(pT**2)'
c          print*,' '
c        endif
ctestend giuseppe

        PDFscale=sqrt(Q)

      ENDIF

c scale choice

      IF(i_PDFscale.EQ.3)THEN
        PDFscale=fixed_PDFscale
      ENDIF

      IF(i_PDFscale.EQ.4)THEN
        ntotlept=0
        do i=1,8
          if (idp_inout(i).gt.0.and.abs(idp(i)).gt.10.and.
     &               abs(idp(i)).lt.17)then
            ntotlept=ntotlept+1
            nlept(ntotlept)=i
          endif
        enddo

        if (ntotlept.eq.4) then

          ptot0=pcoll(0,nlept(1))+pcoll(0,nlept(2))
     &         +pcoll(0,nlept(3))+pcoll(0,nlept(4))
          ptot1=pcoll(1,nlept(1))+pcoll(1,nlept(2))
     &         +pcoll(1,nlept(3))+pcoll(1,nlept(4))
          ptot2=pcoll(2,nlept(1))+pcoll(2,nlept(2))
     &         +pcoll(2,nlept(3))+pcoll(2,nlept(4))
          ptot3=pcoll(3,nlept(1))+pcoll(3,nlept(2))
     &         +pcoll(3,nlept(3))+pcoll(3,nlept(4))

        elseif (ntotlept.eq.2) then

* find outgoing quarks or gluons
          ntotquark=0
          do i=1,8
            if (idp_inout(i).gt.0.and.(abs(idp(i)).lt.6.or.abs(idp(i)).
     &           eq.21))then
              ntotquark=ntotquark+1
              nquark(ntotquark)=i
            endif
          enddo
* compute  cosines
          DO i=1,4
            pmodq=pcoll(1,nquark(i))**2+pcoll(2,nquark(i))**2
     &           +pcoll(3,nquark(i))**2
            pmodq=sqrt(pmodq)
            qucos(i)=pcoll(3,nquark(i))/pmodq
            IF(qucos(i).GT.1.d0.OR.qucos(i).LT.-1.d0) RETURN
          ENDDO
* FIND MOST FORWARD AND MOST BACKWARD JETS
          qucos_max=-1.1d0
          qucos_min=1.1d0
          do i=1,4
            IF(qucos(i).LT.qucos_min)THEN
              qucos_min=qucos(i)
              i_qucos_min=i
            ENDIF
            IF(qucos(i).GT.qucos_max)THEN
              qucos_max=qucos(i)
              i_qucos_max=i
            ENDIF
          enddo
          icentr=1
          do i=1,4
            if (i.ne.i_qucos_min.and.i.ne.i_qucos_max) then
              nqucentr(icentr)=i
              icentr=icentr+1
            endif
          enddo
          ptot0=pcoll(0,nlept(1))+pcoll(0,nlept(2))
     &     +pcoll(0,nquark(nqucentr(1)))+pcoll(0,nquark(nqucentr(2)))
          ptot1=pcoll(1,nlept(1))+pcoll(1,nlept(2))
     &     +pcoll(1,nquark(nqucentr(1)))+pcoll(1,nquark(nqucentr(2)))
          ptot2=pcoll(2,nlept(1))+pcoll(2,nlept(2))
     &     +pcoll(2,nquark(nqucentr(1)))+pcoll(2,nquark(nqucentr(2)))
          ptot3=pcoll(3,nlept(1))+pcoll(3,nlept(2))
     &     +pcoll(3,nquark(nqucentr(1)))+pcoll(3,nquark(nqucentr(2)))
        else

          print*, '                     *****ERROR*****'
          print*, 'this PDF scale is for 4 or 2 leptons in final state'
          print*, 'this process seems to have', ntotlept,'leptons'
          stop
        endif



        Q=(ptot0**2-ptot1**2-ptot2**2-ptot3**2)/2.d0
        if (Q.gt.50d0) then
          PDFscale=sqrt(Q)
        else
          PDFscale=sqrt(50.d0)
        endif
      ENDIF
      
c giovanni 2017/10/24 begin
      IF(i_PDFscale.EQ.5)THEN           
        ntotquark=0
        do i=1,8
            if (idp_inout(i).gt.0.and.(abs(idp(i)).lt.6.or.abs(idp(i)).
     &           eq.21))then
              ntotquark=ntotquark+1
              nquark(ntotquark)=i
            endif
        enddo
        ptjL = 0.d0
        iptjL = 0
        ptjNL = 0.d0
        iptjNL = 0
        if (ntotquark.lt.2) then
          print*, '***********************ERROR*********************'
          print*, 'PDF scale requires at least 2 jets in final state'
          print*, 'this process seems to have', ntotquarks,'jets'
          stop
        else          
          do i=1,ntotquark
           ptjtemp = sqrt(pcoll(1,nquark(i))**2+pcoll(2,nquark(i))**2) 
           if (ptjtemp .gt. ptjL) then
            iptjL = nquark(i)
            ptjL = ptjtemp
           endif               
          enddo
c          print*,'leading jet pt= ',ptjL          
          do i=1,ntotquark
           ptjtemp = sqrt(pcoll(1,nquark(i))**2+pcoll(2,nquark(i))**2)
           if (ptjtemp .gt. ptjNL .and. nquark(i) .NE. iptjL) then
            iptjNL =  nquark(i)
            ptjNL = ptjtemp
           endif           
          enddo
c          print*,'next-to-leading jet pt= ',ptjNL
          PDFscale = sqrt(ptjL*ptjNL)
c          print*,'PDFscale = ',PDFscale   
        endif  
      ENDIF
c giovanni 2017/10/24 end
*1_7
      PDFscale=PDFscale*pdfconst
*1_7end           
c end scale choice
c end giuseppe 01/03/2007


* Use LHAPDF default for points outside the allowed range: frozen at
* boundaries

c      IF(PDFScale.GE.10000.d0.OR.PDFScale.LE.1.d0)THEN
c        fxn=0.d0
c        RETURN
c      ENDIF

*  RENORMALIZATION SCALE SQUARED FOR ALFAS:

* Alfa_s from Pythia
*      Q2_QCD=Q
*      alfa_s=pyalps(Q2_QCD)

* Alfa_s from LHAPDF: factorization scale = renormalization scale
      alfa_s=alphasPDF(PDFScale)
      gs2=alfa_s*fourpi


c giuseppe ILC 26/06/2007
      IF(i_coll.EQ.1 .OR. i_coll.EQ.2) THEN ! hadron colliders
        CALL evolvePDF(b1,PDFScale,ff1)
        CALL evolvePDF(b2,PDFScale,ff2)
      ENDIF
c end giuseppe ILC 26/06/2007


*** accounting for the difference between PDG and LHAPDF convention
      DO i=1,2
        IF(ABS(idin0(i)).EQ.21)THEN
          idin0_pdf(i)=0
        ELSE
          idin0_pdf(i)=idin0(i)
        ENDIF
      ENDDO                     !i


      IF(i_coll.eq.1)THEN       ! p-p collider

        IF(i_ccfam.EQ.0)THEN
          maxdistfun=1
          id1_LH(1)=-idin0_pdf(1)
          id2_LH(1)=-idin0_pdf(2)
          xpdf1_LH(1)=ff1(id1_LH(1))
          xpdf2_LH(1)=ff2(id2_LH(1))
          distfun(1)=xpdf1_LH(1)*xpdf2_LH(1)/b1/b2
        ELSE                    ! sum
          maxdistfun=4
          i1_pdf=ifamexchange12(-idin0_pdf(1))
          i2_pdf=ifamexchange12(-idin0_pdf(2))

          id1_LH(1)=-idin0_pdf(1)
          id2_LH(1)=-idin0_pdf(2)
          xpdf1_LH(1)=ff1(id1_LH(1))
          xpdf2_LH(1)=ff2(id2_LH(1))
          distfun(1)=xpdf1_LH(1)*xpdf2_LH(1)/b1/b2

c* correction for autoconjugate CC and/or FAM processes Sandro

          if (ionesh.eq.1.or.ifirst.eq.1) then
            call ccfcsym(iproc,iccsum,ifcsum,iccfcsum)
ctest
            if (ionesh.ne.1) then
              print*, 'iccsum,ifcsum,iccfcsum'
              print*, iccsum,ifcsum,iccfcsum
            endif
ctestend
            ifirst=0
          endif

* FAMILY_1 <--> FAMILY_2
          if (ifcsum.eq.0) then
            distfun(2)=0.d0
          else
            id1_LH(2)=i1_pdf
            id2_LH(2)=i2_pdf
            xpdf1_LH(2)=ff1(id1_LH(2))
            xpdf2_LH(2)=ff2(id2_LH(2))
            distfun(2)=xpdf1_LH(2)*xpdf2_LH(2)/b1/b2
          endif

* C CONJUGATION
          if (iccsum.eq.0) then
            distfun(3)=0.d0
          else
            id1_LH(3)=idin0_pdf(1)
            id2_LH(3)=idin0_pdf(2)
            xpdf1_LH(3)=ff1(id1_LH(3))
            xpdf2_LH(3)=ff2(id2_LH(3))
            distfun(3)=xpdf1_LH(3)*xpdf2_LH(3)/b1/b2
          endif

*     FAMILY_1 <--> FAMILY_2 + C CONJUGATION
          if (iccfcsum.eq.0) then
            distfun(4)=0.d0
          else
            id1_LH(4)=-i1_pdf
            id2_LH(4)=-i2_pdf
            xpdf1_LH(4)=ff1(id1_LH(4))
            xpdf2_LH(4)=ff2(id2_LH(4))
            distfun(4)=xpdf1_LH(4)*xpdf2_LH(4)/b1/b2
          endif
c giuseppe 16/11/2006
c*     SUM
c        distfuntot=distfun(1)
c        DO i=2,4
c          distfuntot=distfuntot+distfun(i)
c        ENDDO
c      ENDIF      !i_ccfam.EQ.0

        ENDIF                   !i_ccfam.EQ.0
* SUM
        distfuntot=distfun(1)
        DO i=2,maxdistfun
          distfuntot=distfuntot+distfun(i)
        ENDDO
c end giuseppe 16/11/2006

      ELSEIF(i_coll.eq.2)THEN   ! p-pbar collider

c giuseppe 20/11/2006: no bug here, only change in the order for more
c clarity/ease of maintenance. First the whole case 'i_ccfam.EQ.0' is
c accomplished, then the whole case 'i_ccfam.EQ.1'
        IF(i_ccfam.EQ.0)THEN
          maxdistfun=1
          id1_LH(1)=-idin0_pdf(1)
          id2_LH(1)=-idin0_pdf(2)
          xpdf1_LH(1)=ff1(id1_LH(1))
          xpdf2_LH(1)=ff2(-id2_LH(1))
          distfun(1)=xpdf1_LH(1)*xpdf2_LH(1)/b1/b2

          IF(nidenticalinitial.EQ.1)THEN !initial partons non-identical
            maxdistfun=2
            id1_LH(2)=-idin0_pdf(2)
            id2_LH(2)=-idin0_pdf(1)
            xpdf1_LH(2)=ff1(id1_LH(2))
            xpdf2_LH(2)=ff2(-id2_LH(2))
            distfun(2)=xpdf1_LH(2)*xpdf2_LH(2)/b1/b2
          ENDIF

        ELSEIF(i_ccfam.EQ.1)THEN
          maxdistfun=4
          i1_pdf=ifamexchange12(-idin0_pdf(1))
          i2_pdf=ifamexchange12(-idin0_pdf(2))
          id1_LH(1)=-idin0_pdf(1)
          id2_LH(1)=-idin0_pdf(2)
          xpdf1_LH(1)=ff1(id1_LH(1))
          xpdf2_LH(1)=ff2(-id2_LH(1))
          distfun(1)=xpdf1_LH(1)*xpdf2_LH(1)/b1/b2
c* correction for autoconjugate CC and/or FAM processes Sandro

          if (ionesh.eq.1.or.ifirst.eq.1) then
            call ccfcsym(iproc,iccsum,ifcsum,iccfcsum)
ctest
            if (ionesh.ne.1) then
              print*, 'iccsum,ifcsum,iccfcsum'
              print*, iccsum,ifcsum,iccfcsum
            endif
ctestend
            ifirst=0
          endif

* FAMILY_1 <--> FAMILY_2
          if (ifcsum.eq.0) then
            distfun(2)=0.d0
          else
            id1_LH(2)=i1_pdf
            id2_LH(2)=i2_pdf
            xpdf1_LH(2)=ff1(id1_LH(2))
            xpdf2_LH(2)=ff2(-id2_LH(2))
            distfun(2)=xpdf1_LH(2)*xpdf2_LH(2)/b1/b2
          endif

* C CONJUGATION
          if (iccsum.eq.0) then
            distfun(3)=0.d0
          else
            id1_LH(3)=idin0_pdf(1)
            id2_LH(3)=idin0_pdf(2)
            xpdf1_LH(3)=ff1(id1_LH(3))
            xpdf2_LH(3)=ff2(-id2_LH(3))
            distfun(3)=xpdf1_LH(3)*xpdf2_LH(3)/b1/b2
          endif

*     FAMILY_1 <--> FAMILY_2 + C CONJUGATION
          if (iccfcsum.eq.0) then
            distfun(4)=0.d0
          else
            id1_LH(4)=-i1_pdf
            id2_LH(4)=-i2_pdf
            xpdf1_LH(4)=ff1(id1_LH(4))
            xpdf2_LH(4)=ff2(-id2_LH(4))
            distfun(4)=xpdf1_LH(4)*xpdf2_LH(4)/b1/b2
          endif

          IF(nidenticalinitial.EQ.1)THEN
            maxdistfun=8
            i1_pdf=ifamexchange12(-idin0_pdf(2))
            i2_pdf=ifamexchange12(-idin0_pdf(1))
            id1_LH(5)=-idin0_pdf(2)
            id2_LH(5)=-idin0_pdf(1)
            xpdf1_LH(5)=ff1(id1_LH(5))
            xpdf2_LH(5)=ff2(-id2_LH(5))
            distfun(5)=xpdf1_LH(5)*xpdf2_LH(5)/b1/b2
c* correction for autoconjugate CC and/or FAM processes Sandro Aug05

* FAMILY_1 <--> FAMILY_2
            if (ifcsum.eq.0) then
              distfun(6)=0.d0
            else
              id1_LH(6)=i1_pdf
              id2_LH(6)=i2_pdf
              xpdf1_LH(6)=ff1(id1_LH(6))
              xpdf2_LH(6)=ff2(-id2_LH(6))
              distfun(6)=xpdf1_LH(6)*xpdf2_LH(6)/b1/b2
            endif

* C CONJUGATION
            if (iccsum.eq.0) then
              distfun(7)=0.d0
            else
              id1_LH(7)=idin0_pdf(2)
              id2_LH(7)=idin0_pdf(1)
              xpdf1_LH(7)=ff1(id1_LH(7))
              xpdf2_LH(7)=ff2(-id2_LH(7))
              distfun(7)=xpdf1_LH(7)*xpdf2_LH(7)/b1/b2
            endif

*     FAMILY_1 <--> FAMILY_2 + C CONJUGATION
            if (iccfcsum.eq.0) then
              distfun(8)=0.d0
            else
              id1_LH(8)=-i1_pdf
              id2_LH(8)=-i2_pdf
              xpdf1_LH(8)=ff1(id1_LH(8))
              xpdf2_LH(8)=ff2(-id2_LH(8))
              distfun(8)=xpdf1_LH(8)*xpdf2_LH(8)/b1/b2
            endif

          ENDIF                 ! IF(nidenticalinitial.EQ.1)THEN

        ENDIF                   ! IF(i_ccfam.EQ.0)THEN
c end giuseppe 20/11/2006

* correctionend
*     SUM
        distfuntot=distfun(1)
        DO i=2,maxdistfun
          distfuntot=distfuntot+distfun(i)
        ENDDO

c giuseppe ILC 26/06/2007
      ELSEIF(i_coll.EQ.3)THEN   ! e+e- collider

        id1_LH(1)=-idin0_pdf(1)
        id2_LH(1)=-idin0_pdf(2)

        IF(i_isr.EQ.1)THEN
* ISR on
          QISRscale=sqrt(saux)
          xpdf1_LH(1)=ffISR(bb1,QISRscale)
          xpdf2_LH(1)=ffISR(bb2,QISRscale)
          distfun(1)=xpdf1_LH(1)*xpdf2_LH(1)
        ELSE
* ISR off
          distfun(1)=1.d0
        ENDIF                   ! IF(i_isr.EQ.1)THEN

        IF(i_ccfam.EQ.0)THEN
          maxdistfun=1
        ELSEIF(i_ccfam.EQ.1)THEN
          maxdistfun=4
          if (ionesh.eq.1.or.ifirst.eq.1) then
            call ccfcsym(iproc,iccsum,ifcsum,iccfcsum)
ctest
            if (ionesh.ne.1) then
              print*, 'iccsum,ifcsum,iccfcsum'
              print*, iccsum,ifcsum,iccfcsum
            endif
ctestend
            ifirst=0
          endif

* FAMILY_1 <--> FAMILY_2
          if (ifcsum.eq.0) then
            distfun(2)=0.d0
          else
            distfun(2)=distfun(1)
          endif

* C CONJUGATION
          if (iccsum.eq.0) then
            distfun(3)=0.d0
          else
            distfun(3)=distfun(1)
          endif

* FAMILY_1 <--> FAMILY_2  +  C CONJUGATION
          if (iccfcsum.eq.0) then
            distfun(4)=0.d0
          else
            distfun(4)=distfun(1)
          endif
        ENDIF                   ! IF(i_ccfam.EQ.0)THEN

* SUM
        distfuntot=distfun(1)
        DO i=2,maxdistfun
          distfuntot=distfuntot+distfun(i)
        ENDDO
c end giuseppe ILC 26/06/2007

      ENDIF                     !i_coll

      IF(distfuntot.LE.0.d0)THEN
        fxn=0.d0
        RETURN
      ENDIF

* SYMMETRY BETWEEN INCOMING PROTONS
c      IF(idin0(1).NE.idin0(2)) distfuntot=2.d0*distfuntot
*** no more symmetry between incoming protons for single run .
***  for one shot the factor 2 is accounted for in oneshot sampling

*     IF DOING NORMALIZATION RETURN
      IF(i_normalize.EQ.1)THEN
        fxn=distfuntot*rjacys*rjacz
**** number of effective calls
        ncall_eff=ncall_eff+1
        RETURN
      ENDIF
* CHECK FOR IDENTICAL FINAL STATE PARTICLES
      DO i=-16,21
        nidentical(i)=0
      ENDDO
      DO i=1,8
        IF(idp_inout(i).GT.0)THEN
          nidentical(idp(i))=nidentical(idp(i))+1
          ipos0(nidentical(idp(i)),idp(i))=i
        ENDIF
      ENDDO
* FIND MULTIPLICITIES > 1:THERE CAN BE AT MOST 3
* FIND ALSO A PARTICLE WITH MULTIPLICITY ZERO TO BE USED
*     IN REDUNDANT DO-LOOPS
      nmultiple=0
      DO i=-16,21
        IF(nidentical(i).GT.1)THEN
          nmultiple=nmultiple+1
          i_identical(nmultiple)=i
        ENDIF
        IF(nidentical(i).EQ.0) i_multzero=i
      ENDDO
      IF(nmultiple.GT.3)THEN
        PRINT*,'nmultiple=',nmultiple,'    SOMETHING IS WRONG!'
        STOP
      ENDIF
      DO i= nmultiple+1,3
        i_identical(i)=i_multzero
      ENDDO
      symfact= ifact(nidentical(i_identical(3)))
     &     *ifact(nidentical(i_identical(2)))
     &     *ifact(nidentical(i_identical(1)))

      denfac=0.d0
** Ezio 18/03/04
      izerojac=0
** Ezio-end 18/03/04

* COMPUTE JACOBIANS, SYMMETRIZE ON IDENTICAL PARTICLES
      DO i3=1,ifact(nidentical(i_identical(3)))
        DO i2=1,ifact(nidentical(i_identical(2)))
          DO i1=1,ifact(nidentical(i_identical(1)))
* LOAD BASIC ORDERING
            DO l=1,8
              iordaux(l)=l
            ENDDO
* PERMUTE INDICES OF IDENTICAL PARTICLES
            DO l=1,nidentical(i_identical(1))
              iordtemp(index(l,i1))=ipos0(l,i_identical(1))
            ENDDO
            DO l=1,nidentical(i_identical(1))
              iordaux(ipos0(l,i_identical(1)))=iordtemp(l)
            ENDDO
*test
*            print*, ' '
*            print*, i_identical(2),nidentical(i_identical(2))
*testend
            DO l=1,nidentical(i_identical(2))
*test
*              print*,index(l,i2),ipos0(l,i_identical(2))
*testend
              iordtemp(index(l,i2))=ipos0(l,i_identical(2))
            ENDDO
            DO l=1,nidentical(i_identical(2))
              iordaux(ipos0(l,i_identical(2)))=iordtemp(l)
            ENDDO
            DO l=1,nidentical(i_identical(3))
              iordtemp(index(l,i3))=ipos0(l,i_identical(3))
            ENDDO
            DO l=1,nidentical(i_identical(3))
              iordaux(ipos0(l,i_identical(3)))=iordtemp(l)
            ENDDO
* LOOP ON PHASE SPACES
            DO k=1,mxnphs
              IF(ialfa(k).NE.0)THEN
* PERMUTE IORDER
                DO l=1,8
                  iorderaux(l)=iordaux(iorder(l,k))
                ENDDO

                IF(k.le.13)THEN
                  CALL phsp1_1_4_multi7_cjac(shat,
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm11c(k),gam11c(k),spm11clo(k),spm11chi(k),
     &                 rm12a(k),gam12a(k),spm12alo(k),spm12ahi(k),
     &                 rm12b(k),gam12b(k),spm12blo(k),spm12bhi(k),
     &                 rm12c(k),gam12c(k),spm12clo(k),spm12chi(k),
     &                 rmta(k),expa(k),spta,rmtb(k),expb(k),sptb,
     &                 pt5min,pt6min,
c giuseppe 07/07/2007
     &                 idp(iorder(7,k)),idp(iorder(8,k)),
c end giuseppe 07/07/2007
     &                 xmasq,xmbsq,
     &                 xm,xmsq,iorderaux(3),
     &                 p,xjac)
                ELSE IF(k.le.24)THEN
                  CALL phsp2_4_multi7jac(shat,
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm11c(k),gam11c(k),spm11clo(k),spm11chi(k),
     &                 rm12a(k),gam12a(k),spm12alo(k),spm12ahi(k),
     &                 rm12b(k),gam12b(k),spm12blo(k),spm12bhi(k),
     &                 rm12c(k),gam12c(k),spm12clo(k),spm12chi(k),
     &                 rm2a(k),gam2a(k),spm2alo(k),spm2ahi(k),
     &                 rm2b(k),gam2b(k),spm2blo(k),spm2bhi(k),
     &                 rm2c(k),gam2c(k),spm2clo(k),spm2chi(k),
     &                 xm,xmsq,iorderaux(3),
     &                 p,xjac)
                ELSE IF(k.le.32)THEN
                  CALL phsp1_1_31_multi_cjac(shat,
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm111a(k),gam111a(k),spm111alo(k),spm111ahi(k),
     &                 rm111b(k),gam111b(k),spm111blo(k),spm111bhi(k),
     &                 rmta(k),expa(k),spta,rmtb(k),expb(k),sptb,
     &                 pt5min,pt6min,
c giuseppe 07/07/2007
     &                 idp(iorder(7,k)),idp(iorder(8,k)),
c end giuseppe 07/07/2007
     &                 xmasq,xmbsq,
     &                 xm,xmsq,iorderaux(3),
     &                 p,xjac)
                ELSE IF(k.le.54)THEN
                  CALL phsp1_2_3_multi5_cjac(shat,
     &                 rm0a(k),gam0a(k),spm0alo(k),spm0ahi(k),
     &                 rm0b(k),gam0b(k),spm0blo(k),spm0bhi(k),
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm2a(k),gam2a(k),spm2alo(k),spm2ahi(k),
     &                 rm2b(k),gam2b(k),spm2blo(k),spm2bhi(k),
     &                 rmta(k),expa(k),spta,
     &                 pt6min,
     &                 idin0(1),idin0(2),idp(iorder(8,k)),
     &                 xmasq,xmbsq,
     &                 xm,xmsq,iorderaux(3),
     &                 p,xjac)
                ELSE IF(k.le.60)THEN
                  CALL phsp3_3_multi5jac(shat,
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm2a(k),gam2a(k),spm2alo(k),spm2ahi(k),
     &                 rm2b(k),gam2b(k),spm2blo(k),spm2bhi(k),
     &                 rm21a(k),gam21a(k),spm21alo(k),spm21ahi(k),
     &                 rm21b(k),gam21b(k),spm21blo(k),spm21bhi(k),
     &                 xm,xmsq,iorderaux(3),
     &                 p,xjac)
                ELSE IF(k.le.68)THEN
                  CALL phsp2_4to31_multi5jac(shat,
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm111a(k),gam111a(k),spm111alo(k),spm111ahi(k),
     &                 rm111b(k),gam111b(k),spm111blo(k),spm111bhi(k),
     &                 rm2a(k),gam2a(k),spm2alo(k),spm2ahi(k),
     &                 rm2b(k),gam2b(k),spm2blo(k),spm2bhi(k),
     &                 xm,xmsq,iorderaux(3),
     &                 p,xjac)
                ELSE IF(k.le.76)THEN
                  CALL phsp1_5to1_4to31_multi_cjac(shat,
     &                 rm1(k),gam1(k),spm1lo(k),spm1hi(k),
     &                 rm11(k),gam11(k),spm11lo(k),spm11hi(k),
     &                 rm111(k),gam111(k),spm111lo(k),spm111hi(k),
     &                 rm1111(k),gam1111(k),spm1111lo(k),spm1111hi(k),
     &                 rmta(k),expa(k),spta,
     &                 pt6min,
     &                 idin0(1),idin0(2),idp(iorder(8,k)),
     &                 xmasq,xmbsq,
     &                 xm,xmsq,iorderaux(3),
     &                 p,xjac)
                ELSE
                  print*,'***fxn.f ERROR'
                  print*,'   iphs_ind:',iphs_ind,
     &                 'exceeding the maximum number of channels'
                  STOP
                ENDIF

************************************************************************
                IF(xjac.GT.0.d0)THEN

                  IF(nidenticalinitial.EQ.2  .AND.
     &                 (k.EQ.2 .OR. k.EQ.3
     &                 .OR. k.EQ.25 .OR. k.EQ.26
     &                 .OR. k.EQ.27 .OR. k.EQ.28
     &                 .OR. k.EQ.34 .OR. k.EQ.35
     &                 .OR. k.EQ.43 .OR. k.EQ.44
     &                 .OR. k.EQ.45 .OR. k.EQ.46
     &                 .OR. k.EQ.51 .OR. k.EQ.52
     &                 .OR. k.EQ.53 .OR. k.EQ.54)
     &                 )THEN
                    denfac=denfac+alfa(k)/xjac/symfact/rnormal(k)
     &                   /2.d0
                  ELSE
                    denfac=denfac+alfa(k)/xjac/symfact/rnormal(k)
                  ENDIF

                ELSE

                  izerojac=1

                ENDIF           ! ENDIF xjac.GT.0.d0
              ENDIF             ! ENDIF(ialfa(k).NE.0)
            ENDDO               ! ENDDO k PHASE SPACES
          ENDDO                 ! ENDDO i1 PERMUTATIONS2
        ENDDO                   ! ENDDO i2 PERMUTATIONS2
      ENDDO                     ! ENDDO i3 PERMUTATIONS1

* FLIP FINAL STATE MOMENTA FOR PHASE SPACES ,,,,,,,,,,,,
* IN ORDER TO ACCOUNT FOR SYMMETRY BETWEEN INCOMING PARTICLES
* ...
***********************************************************************
c giuseppe 01/06/2007
c      IF(nidenticalinitial.EQ.2)THEN
      IF(nidenticalinitial.EQ.2 .AND. izerojac.EQ.0)THEN
c end giuseppe 01/06/2007
        DO i=1,8
          DO j=0,3
            paux(j,i)=p(j,i)
          ENDDO
        ENDDO
* FLIP
        DO i=3,8
          paux(3,iorder(i,iphs_ind))=-paux(3,iorder(i,iphs_ind))
          paux(2,iorder(i,iphs_ind))=-paux(2,iorder(i,iphs_ind))
	ENDDO
*     COMPUTE EXTRA JACOBIANS
        DO i3=1,ifact(nidentical(i_identical(3)))
          DO i2=1,ifact(nidentical(i_identical(2)))
            DO i1=1,ifact(nidentical(i_identical(1)))
* LOAD BASIC ORDERING
              DO l=1,8
                iordaux(l)=l
              ENDDO
* PERMUTE INDICES OF IDENTICAL PARTICLES
              DO l=1,nidentical(i_identical(1))
                iordtemp(index(l,i1))=ipos0(l,i_identical(1))
              ENDDO
              DO l=1,nidentical(i_identical(1))
                iordaux(ipos0(l,i_identical(1)))=iordtemp(l)
              ENDDO
              DO l=1,nidentical(i_identical(2))
                iordtemp(index(l,i2))=ipos0(l,i_identical(2))
              ENDDO
              DO l=1,nidentical(i_identical(2))
                iordaux(ipos0(l,i_identical(2)))=iordtemp(l)
              ENDDO
              DO l=1,nidentical(i_identical(3))
                iordtemp(index(l,i3))=ipos0(l,i_identical(3))
              ENDDO
              DO l=1,nidentical(i_identical(3))
                iordaux(ipos0(l,i_identical(3)))=iordtemp(l)
              ENDDO
* LOOP ON PHASE SPACES
              DO k=2,3
                IF(ialfa(k).NE.0)THEN
* PERMUTE IORDER
                  DO l=1,8
                    iorderaux(l)=iordaux(iorder(l,k))
                  ENDDO
                  CALL phsp1_1_4_multi7_cjac(shat,
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm11c(k),gam11c(k),spm11clo(k),spm11chi(k),
     &                 rm12a(k),gam12a(k),spm12alo(k),spm12ahi(k),
     &                 rm12b(k),gam12b(k),spm12blo(k),spm12bhi(k),
     &                 rm12c(k),gam12c(k),spm12clo(k),spm12chi(k),
     &                 rmta(k),expa(k),spta,rmtb(k),expb(k),sptb,
     &                 pt5min,pt6min,
c giuseppe 07/07/2007
     &                 idp(iorder(7,k)),idp(iorder(8,k)),
c end giuseppe 07/07/2007
     &                 xmasq,xmbsq,
     &                 xm,xmsq,iorderaux(3),
     &                 paux,xjac)
                  IF(xjac.GT.0.d0)THEN
                    denfac=denfac+alfa(k)/xjac/symfact/rnormal(k)
     &                   /2.d0
                  ELSE
                    izerojac=1
                  ENDIF         !
                ENDIF           ! ENDIF(ialfa(k).NE.0)
              ENDDO             ! ENDDO k PHASE SPACES
              DO k=25,28
                IF(ialfa(k).NE.0)THEN
* PERMUTE IORDER
                  DO l=1,8
                    iorderaux(l)=iordaux(iorder(l,k))
                  ENDDO
                  CALL phsp1_1_31_multi_cjac(shat,
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm111a(k),gam111a(k),spm111alo(k),spm111ahi(k),
     &                 rm111b(k),gam111b(k),spm111blo(k),spm111bhi(k),
     &                 rmta(k),expa(k),spta,rmtb(k),expb(k),sptb,
     &                 pt5min,pt6min,
c giuseppe 07/07/2007
     &                 idp(iorder(7,k)),idp(iorder(8,k)),
c end giuseppe 07/07/2007
     &                 xmasq,xmbsq,
     &                 xm,xmsq,iorderaux(3),
c Giuseppe 2/8/06     &                 p,xjac)
     &                 paux,xjac)

                  IF(xjac.GT.0.d0)THEN
                    denfac=denfac+alfa(k)/xjac/symfact/rnormal(k)
     &                   /2.d0
                  ELSE
                    izerojac=1
                  ENDIF         !
                ENDIF           ! ENDIF(ialfa(k).NE.0)
              ENDDO             ! ENDDO k PHASE SPACES
              DO k=34,35
                IF(ialfa(k).NE.0)THEN
* PERMUTE IORDER
                  DO l=1,8
                    iorderaux(l)=iordaux(iorder(l,k))
                  ENDDO
                  CALL phsp1_2_3_multi5_cjac(shat,
     &                 rm0a(k),gam0a(k),spm0alo(k),spm0ahi(k),
     &                 rm0b(k),gam0b(k),spm0blo(k),spm0bhi(k),
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm2a(k),gam2a(k),spm2alo(k),spm2ahi(k),
     &                 rm2b(k),gam2b(k),spm2blo(k),spm2bhi(k),
     &                 rmta(k),expa(k),spta,
     &                 pt6min,
     &                 idin0(1),idin0(2),idp(iorder(8,k)),
     &                 xmasq,xmbsq,
     &                 xm,xmsq,iorderaux(3),
c Giuseppe 2/8/06     &                 p,xjac)
     &                 paux,xjac)

                  IF(xjac.GT.0.d0)THEN
                    denfac=denfac+alfa(k)/xjac/symfact/rnormal(k)
     &                   /2.d0
                  ELSE
                    izerojac=1
                  ENDIF         !
                ENDIF           ! ENDIF(ialfa(k).NE.0)
              ENDDO             ! ENDDO k PHASE SPACES
              DO k=43,46
                IF(ialfa(k).NE.0)THEN
* PERMUTE IORDER
                  DO l=1,8
                    iorderaux(l)=iordaux(iorder(l,k))
                  ENDDO
                  CALL phsp1_2_3_multi5_cjac(shat,
     &                 rm0a(k),gam0a(k),spm0alo(k),spm0ahi(k),
     &                 rm0b(k),gam0b(k),spm0blo(k),spm0bhi(k),
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm2a(k),gam2a(k),spm2alo(k),spm2ahi(k),
     &                 rm2b(k),gam2b(k),spm2blo(k),spm2bhi(k),
     &                 rmta(k),expa(k),spta,
     &                 pt6min,
     &                 idin0(1),idin0(2),idp(iorder(8,k)),
     &                 xmasq,xmbsq,
     &                 xm,xmsq,iorderaux(3),
     &                 paux,xjac)

                  IF(xjac.GT.0.d0)THEN
                    denfac=denfac+alfa(k)/xjac/symfact/rnormal(k)
     &                   /2.d0
                  ELSE
                    izerojac=1
                  ENDIF         !
                ENDIF           ! ENDIF(ialfa(k).NE.0)
              ENDDO             ! ENDDO k PHASE SPACES
              DO k=51,54
                IF(ialfa(k).NE.0)THEN
* PERMUTE IORDER
                  DO l=1,8
                    iorderaux(l)=iordaux(iorder(l,k))
                  ENDDO
                  CALL phsp1_2_3_multi5_cjac(shat,
     &                 rm0a(k),gam0a(k),spm0alo(k),spm0ahi(k),
     &                 rm0b(k),gam0b(k),spm0blo(k),spm0bhi(k),
     &                 rm1a(k),gam1a(k),spm1alo(k),spm1ahi(k),
     &                 rm1b(k),gam1b(k),spm1blo(k),spm1bhi(k),
     &                 rm11a(k),gam11a(k),spm11alo(k),spm11ahi(k),
     &                 rm11b(k),gam11b(k),spm11blo(k),spm11bhi(k),
     &                 rm2a(k),gam2a(k),spm2alo(k),spm2ahi(k),
     &                 rm2b(k),gam2b(k),spm2blo(k),spm2bhi(k),
     &                 rmta(k),expa(k),spta,
     &                 pt6min,
     &                 idin0(1),idin0(2),idp(iorder(8,k)),
     &                 xmasq,xmbsq,
     &                 xm,xmsq,iorderaux(3),
     &                 paux,xjac)
                  IF(xjac.GT.0.d0)THEN
                    denfac=denfac+alfa(k)/xjac/symfact/rnormal(k)
     &                   /2.d0
                  ELSE
                    izerojac=1
                  ENDIF         !
                ENDIF           ! ENDIF(ialfa(k).NE.0)
              ENDDO             ! ENDDO k PHASE SPACES
            ENDDO               ! ENDDO i1 PERMUTATIONS2
          ENDDO                 ! ENDDO i2 PERMUTATIONS2
        ENDDO                   ! ENDDO i3 PERMUTATIONS1
      ENDIF                     ! IF(nidenticalinitial.EQ.2 ...)

      IF(izerojac.EQ.1)THEN
        fxn=0.d0
        RETURN
      ENDIF
** Ezio-end 18/03/04

***elena (24/02/04)
      IF(denfac.LE.0.d0)THEN
        fxn=0.d0
        RETURN
      ENDIF
***elenaend

* ROTATE WITH PROBABILITY 0.5 ALL MOMENTA BY 180 DEGREES AROUND x-AXIS
      IF(i_coll.eq.1)THEN       ! p-p collider

        IF (nidenticalinitial.EQ.2 .OR.
     &       (ionesh.EQ.1.AND.i_exchincoming.eq.1)) THEN
          IF(ran2(idum).GT.0.5d0)THEN
            DO i=1,8
              pcoll(3,i)=-pcoll(3,i)
              pcoll(2,i)=-pcoll(2,i)
            ENDDO
c giuseppe 19/12/2006
            iflip=-iflip
c end giuseppe 19/12/2006
          ENDIF
        ENDIF

      ELSEIF(i_coll.eq.2)THEN   ! p-pbar collider

        IF (nidenticalinitial.EQ.2 ) THEN
          IF(ran2(idum).GT.0.5d0)THEN
            DO i=1,8
              pcoll(3,i)=-pcoll(3,i)
              pcoll(2,i)=-pcoll(2,i)
            ENDDO
c giuseppe 19/12/2006
            iflip=-iflip
c end giuseppe 19/12/2006
          ENDIF
        ENDIF

c In case of e+e- collider, initial-state particles are never
c identical. Hence no rotation has to be performed.
c     ELSEIF(i_coll.eq.3)THEN ! e+e- collider

      ENDIF                     ! i_coll.eq.1

*lab
*     CALL AMPLITUDE: USE CM MOMENTA
c
c      DO i=1,8
c        DO mu=0,3
c          p(mu,i)=p(mu,i)*idp_inout(i)
c        ENDDO                   !mu
c      ENDDO                     !i
cctot
c*sandro2/3/07
ccc giuseppe 26/02/2007
ccc      if (i_nofullew.eq.-1.)then
cc      if (i_nofullew.eq.-1)then
ccc end giuseppe 26/02/2007
c      if (i_pertorder.eq.1)then
c*sandro2/3/07end
c        CALL amp(p,idp,igr,igo,ionesh,symfact,amod,iordamp)
c        amod=amod*elcharge2**6
c      else
cctotend
c        IF(ngluon.EQ.0)THEN
c          CALL ampqcd(p,idp,igr,igo,ionesh,symfact,amod,iordamp)
c        ELSEIF(ngluon.EQ.2)THEN
c          CALL amp2g(p,idp,igr(5),igo(1,5),ionesh,symfact,ngin,amod,
c     &       iordamp)
c        ENDIF
cctot
c      endif
cctotend

*     CALL AMPLITUDE: USE lab (=coll) MOMENTA

      DO i=1,8
        CALL ZBOOST(gcm,bcm,p(0,i),plab(0,i))
      ENDDO

      DO i=1,8
        DO mu=0,3
          plab(mu,i)=plab(mu,i)*idp_inout(i)
        ENDDO                   !mu
      ENDDO                     !i


c     giovanni 24/4/18 begin
************************** ON SHELL PROJECTIONS: BEGIN **************************             
        if (i_osp .eq. 0) then
           offactor = 1.d0
        elseif (i_osp .eq. 1) then
           call singleosp(plab,offactor)
c           print*,'fxn.f: offactor (osp 1) = ', offactor
        elseif (i_osp .eq. 2) then
           call doubleosp(plab,offactor)
        endif
************************** ON SHELL PROJECTIONS: END **************************
c     giovanni 24/4/18 end

*4cmpol
      if (i_4cmpol.eq.0) then
        do i=1,8
          do mu=0,3
            pcall(mu,i)=plab(mu,i)
          enddo
        enddo
      elseif (i_4cmpol.eq.1) then
        do mu=0,3
          p4cm(mu)=0.d0
        enddo
        irep=0
        do i=1,8
          if (idp_inout(i).gt.0) then
            do j=1,4
              if (idp(i).eq.id4cmpol(j)) then
                irep=irep+1
                do mu=0,3
                  p4cm(mu)=p4cm(mu)+plab(mu,i)
                enddo
              endif
            enddo
          endif
        enddo
        if (irep.ne.4) then
          print*, 'ERROR   irep different from 4 !'
          stop
        endif
        do i=1,8
          call boostinv(plab(0,i),p4cm,pcall(0,i))
        enddo
      else
        print *, 'ERROR i4cmpol value not admitted !'
        STOP
      endif
        
c      if (i_pertorder.eq.1)then
c        CALL amp(plab,idp,igr,igo,ionesh,symfact,amod,iordamp)
c        amod=amod*elcharge2**6
c      else
c        IF(ngluon.EQ.0)THEN
c          CALL ampqcd(plab,idp,igr,igo,ionesh,symfact,amod,iordamp)
c        ELSEIF(ngluon.EQ.2)THEN
c          CALL amp2g(plab,idp,igr(5),igo(1,5),ionesh,symfact,ngin,amod,
c     &       iordamp)
c        ENDIF
c      endif

      if (i_pertorder.eq.1)then
        CALL amp(pcall,idp,igr,igo,ionesh,symfact,amod,iordamp)
        amod=amod*elcharge2**6
      else
        IF(ngluon.EQ.0)THEN
          CALL ampqcd(pcall,idp,igr,igo,ionesh,symfact,amod,iordamp)
        ELSEIF(ngluon.EQ.2)THEN
          CALL amp2g(pcall,idp,igr(5),igo(1,5),ionesh,symfact,ngin,amod,
     &       iordamp)
        ENDIF
      endif

*4cmpolend

c     giovanni 24/4/18 begin
c      print*,'fxn.f: amod (osp 1) = ', amod
      amod=amod*offactor
c      print*,'fxn.f: amod rescaled (osp 1) = ', amod
c     giovanni 24/4/18 end
      
      
      if (ISNAN(amod)) then
         print*,'fxn.f: ISNAN(amod) (osp 1) .... STOP'
         STOP
      endif
       
     

*labend

      fxn=alfa(iphs_ind)*distfuntot*amod*fluxhat*factor*rjacys*rjacz/
     &     denfac/rnormal(iphs_ind)


* REORDER MOMENTA AND IDENTITIES AS DETERMINED BY THE AMPLITUDE
      DO i=1,8
        DO mu=0,3
          paux(mu,iordamp(i))=pcoll(mu,i)
        ENDDO                   !mu
      ENDDO                     !i
      DO i=1,8
        DO mu=0,3
          pcoll(mu,i)=paux(mu,i)
        ENDDO                   !mu
      ENDDO                     !i
* IDP AND IDP_INOUT MUST BE PRESERVED
      DO i=1,8
        idp0(i)=idp(i)
        idp0_inout(i)=idp_inout(i)
      ENDDO                     !i
      DO i=1,8
        iorderaux(iordamp(i))=idp0(i)
        iordaux(iordamp(i))=idp0_inout(i)
      ENDDO                     !i
      DO i=1,8
        idp0(i)=iorderaux(i)
        idp0_inout(i)=iordaux(i)
      ENDDO                     !i

* DETERMINE THE FAMILY AND CC PARTICLE SET IF NECESSARY
      IF(i_coll.eq.1)THEN       ! p-p collider

        IF(i_ccfam.EQ.1)THEN
c giuseppe 20/11/2006
c          DO i=1,4
          DO i=1,maxdistfun
c end giuseppe 20/11/2006
            distfun(i)=distfun(i)/distfuntot
          ENDDO                 !i
          i=1
          distaux=distfun(1)
          zz=ran2(idum)
          DO WHILE (zz.gt.distaux .and. i.LT.maxdistfun)
            i=i+1
            distaux=distaux+distfun(i)
          ENDDO
*
          IF(i.EQ.2)THEN
            DO j=1,8
              IF(abs(idp0(j)).lt.6) THEN
                idp0(j)=ifamexchange12(idp0(j))
              ENDIF
            ENDDO
          ELSEIF(i.EQ.3)THEN
            DO j=1,8
              IF(idp0(j).LE.16) THEN
                idp0(j)=-idp0(j)
              ENDIF
            ENDDO
          ELSEIF(i.EQ.4)THEN
            DO j=1,8
              IF(abs(idp0(j)).lt.6) THEN
                idp0(j)=-ifamexchange12(idp0(j))
              ELSE
                IF(idp0(j).LE.16) THEN
                  idp0(j)=-idp0(j)
                ENDIF
              ENDIF
            ENDDO
          ENDIF
* CHANGE FOR PARITY IF CC EXCHANGE (amplitude must remain the same)
          IF(i.EQ.1.OR.i.EQ.2)THEN
            ichargeconj=0
          ELSEIF(i.EQ.3.OR.i.EQ.4)THEN
            ichargeconj=1
            DO j=1,8
              DO k=1,3
                pcoll(k,j)=-pcoll(k,j)
              ENDDO
            ENDDO
c giuseppe 19/12/2006
            iflip=-iflip
c end giuseppe 19/12/2006
          ENDIF
* Save i for LHFile output
          i_selected=i
        ELSE
          i_selected=1
        ENDIF                   ! i_ccfam.EQ.1

      ELSEIF(i_coll.eq.2)THEN   ! p-pbar collider

        DO i=1,maxdistfun
          distfun(i)=distfun(i)/distfuntot
        ENDDO                   !i
        i=1
        distaux=distfun(1)
        zz=ran2(idum)
        DO WHILE (zz.GT.distaux .and. i.LT.maxdistfun)
          i=i+1
          distaux=distaux+distfun(i)
        ENDDO
*
c giuseppe 20/11/2006: added IF(i_ccfam.EQ.1) ... ELSEIF(i_ccfam.EQ.0)
        IF(i_ccfam.EQ.1)THEN
          IF(i.EQ.2)THEN
            DO j=1,8
              IF(abs(idp0(j)).lt.6) THEN
                idp0(j)=ifamexchange12(idp0(j))
              ENDIF
            ENDDO
          ELSEIF(i.EQ.3)THEN
            DO j=1,8
              IF(idp0(j).LE.16) THEN
                idp0(j)=-idp0(j)
              ENDIF
            ENDDO
          ELSEIF(i.EQ.4)THEN
            DO j=1,8
              IF(abs(idp0(j)).lt.6) THEN
                idp0(j)=-ifamexchange12(idp0(j))
              ELSE
                IF(idp0(j).LE.16) THEN
                  idp0(j)=-idp0(j)
                ENDIF
              ENDIF
            ENDDO
          ELSEIF(i.EQ.5)THEN    ! no flavour change
          ELSEIF(i.EQ.6)THEN
            DO j=1,8
              IF(abs(idp0(j)).lt.6) THEN
                idp0(j)=ifamexchange12(idp0(j))
              ENDIF
            ENDDO
          ELSEIF(i.EQ.7)THEN
            DO j=1,8
              IF(idp0(j).LE.16) THEN
                idp0(j)=-idp0(j)
              ENDIF
            ENDDO
          ELSEIF(i.EQ.8)THEN
            DO j=1,8
              IF(abs(idp0(j)).lt.6) THEN
                idp0(j)=-ifamexchange12(idp0(j))
              ELSE
                IF(idp0(j).LE.16) THEN
                  idp0(j)=-idp0(j)
                ENDIF
              ENDIF
            ENDDO
          ENDIF                 ! IF(i.EQ.2)THEN

* CHANGE FOR PARITY IF CC EXCHANGE (amplitude must remain the same)
* MUST ALSO ROTATE by 180 deg AROUND x TO BRING PROTON ALONG +z
* THE TWO OPERATION COMBINE TO x --> -x

          IF(i.EQ.1.OR.i.EQ.2)THEN
            ichargeconj=0
          ELSEIF(i.EQ.3.OR.i.EQ.4)THEN
* k --> -k + Rot_x
            ichargeconj=1
            DO j=1,8
              pcoll(1,j)=-pcoll(1,j)
            ENDDO
          ELSEIF(i.EQ.5.OR.i.EQ.6)THEN
*  Rot_x
            ichargeconj=0
            DO j=1,8
              pcoll(2,j)=-pcoll(2,j)
              pcoll(3,j)=-pcoll(3,j)
            ENDDO
c giuseppe 19/12/2006
            iflip=-iflip
c end giuseppe 19/12/2006
          ELSEIF(i.EQ.7.OR.i.EQ.8)THEN
* k --> -k (+ Rot_x + Rot_x)
            ichargeconj=1
            DO j=1,8
              DO k=1,3
                pcoll(k,j)=-pcoll(k,j)
              ENDDO
            ENDDO
c giuseppe 19/12/2006
            iflip=-iflip
c end giuseppe 19/12/2006
          ENDIF

        ELSEIF(i_ccfam.EQ.0)THEN
          ichargeconj=0
          IF(i.EQ.2)THEN
*  Rot_x
            DO j=1,8
              pcoll(2,j)=-pcoll(2,j)
              pcoll(3,j)=-pcoll(3,j)
            ENDDO
c giuseppe 19/12/2006
            iflip=-iflip
c end giuseppe 19/12/2006
          ENDIF

        ENDIF                   ! IF(i_ccfam.EQ.1)THEN
c end giuseppe 20/11/2006

* Save i for LHFile output
        i_selected=i

c giuseppe ILC 26/06/2007
      ELSEIF(i_coll.eq.3)THEN   ! e+e- collider

        IF(i_ccfam.EQ.1)THEN

          DO i=1,maxdistfun
            distfun(i)=distfun(i)/distfuntot
          ENDDO                 !i
          i=1
          distaux=distfun(1)
          zz=ran2(idum)
          DO WHILE (zz.gt.distaux .and. i.LT.maxdistfun)
            i=i+1
            distaux=distaux+distfun(i)
          ENDDO
*
          IF(i.EQ.2)THEN        ! FAM
            DO j=1,8
              IF(abs(idp0(j)).lt.6) THEN
                idp0(j)=ifamexchange12(idp0(j))
              ENDIF
            ENDDO
          ELSEIF(i.EQ.3)THEN    ! CC
            DO j=1,8
              IF(idp0(j).LE.16) THEN
                idp0(j)=-idp0(j)
              ENDIF
            ENDDO
          ELSEIF(i.EQ.4)THEN    ! FAM+CC
            DO j=1,8
              IF(abs(idp0(j)).lt.6) THEN
                idp0(j)=-ifamexchange12(idp0(j))
              ELSE
                IF(idp0(j).LE.16) THEN
                  idp0(j)=-idp0(j)
                ENDIF
              ENDIF
            ENDDO
          ENDIF


c giuseppe 30/10/2007

c* CHANGE FOR PARITY IF CC EXCHANGE (amplitude must remain the same)
c* MUST ALSO ROTATE by 180 deg AROUND x TO BRING ELECTRON ALONG +z
c* THE TWO OPERATION COMBINE TO x --> -x
c          IF(i.EQ.1.OR.i.EQ.2)THEN
c            ichargeconj=0
c          ELSEIF(i.EQ.3.OR.i.EQ.4)THEN
c* k --> -k + Rot_x
c            ichargeconj=1
c            DO j=1,8
c              pcoll(1,j)=-pcoll(1,j)
c            ENDDO
c          ENDIF

* CHANGE FOR PARITY IF CC EXCHANGE (amplitude must remain the same)
          IF(i.EQ.1.OR.i.EQ.2)THEN
            ichargeconj=0
          ELSEIF(i.EQ.3.OR.i.EQ.4)THEN
            ichargeconj=1
            DO j=1,8
              DO k=1,3
                pcoll(k,j)=-pcoll(k,j)
              ENDDO
            ENDDO
          ENDIF

c end giuseppe 30/10/2007

* Save i for LHFile output
          i_selected=i
        ELSE
          ichargeconj=0
          i_selected=1
        ENDIF                   ! i_ccfam.EQ.1
c end giuseppe ILC 26/06/2007

      ENDIF                     ! i_coll

      IF(iflat.EQ.1)THEN
* FIND MAXIMUM WEIGHT
*6
c        IF(it.EQ.(itmx-1) .AND. init.GE.1)THEN
c          rmaxfxn=max(rmaxfxn,fxn*wgt)
c        ENDIF
c        IF(it.EQ.itmx .AND. init.GE.1)THEN
c          rmaxfxn_2=max(rmaxfxn_2,fxn*wgt)
c        ENDIF

        rmaxfxnnew=max(rmaxfxnnew,fxn*wgt)

*6end
      ENDIF

      IF((iflat.EQ.1 .OR. ionesh.EQ.1)
*6
c     &     .AND. it.EQ.itmx.AND.init.GE.1)THEN
     &     .AND.init.GE.1)THEN
*6end
        IF(fxn*wgt.GT.rmaxfxn*scalemax)THEN
          novermax=novermax+1
        ENDIF
        ry=scalemax*rmaxfxn*ran2(idum)
        IF(ry.LE.fxn*wgt)THEN   ! EVENT ACCEPTED
*6
c          NEVHEP=NEVHEP+1
c          nevent=NEVHEP
          nevent=nevent+1
*6end

* STORE INFORMATION IN LESHOUCHES-PYTHIA FORMAT
c giuseppe 29/09/2007
          CALL store_LH_information(iphsp,idp0,idp0_inout,xm,pcoll,
     &         iordamp,ichargeconj,PDFscale)
c end giuseppe 29/09/2007

* CALL PS PROGRAM
          IF(ihadronize.EQ.1)THEN
            CALL PYEVNT
* WRITE OUTPUT OF PS PROGRAM
*  HERE EVENTUALLY CALL A USER PROVIDED ROUTINE FOR WRITING RESULTS,
*  OR SIMPLY USE PYLIST BELOW FOR IPRDEB (GIVEN IN INPUT) EVENTS.
            IF (NEVENT.LE.IPRDEB) CALL PYLIST(2)
          ENDIF  !ihadronize

*  WRITE LESHOUCHES EVENT: UNIT=23 is phamom.dat
*old*          IF(iwrite_event.EQ.1)THEN
*old*            WRITE(23,99999)(IDUP(i),i=1,8)
*old*            WRITE(23,99998)(ISTUP(i),i=1,8)
*old*            DO i=1,8
*old*              WRITE(23,99997)(PUP(j,i),j=1,4)
*old*            ENDDO
*old*            WRITE(23,99996)((ICOLUP(j,i),j=1,2),i=1,8)
*old*	    WRITE(23,99995)SCALUP,AQEDUP,AQCDUP
*old*          ENDIF
          IF(iwrite_event.EQ.1)THEN
            WRITE(23,'(A)')'<event>'
c old            WRITE(23,'(2(I4,1x),4(E13.6,1x))')
            WRITE(23,'(2(I4,1x),4(E21.14,1x))')
     &              NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
            DO i=1,NUP
c old              WRITE(23,'(6(I4,1x),5(E13.6,1x),2(E10.1))')
              WRITE(23,'(6(I4,1x),5(E21.14,1x),2(E10.1))')
     &                    IDUP(i),ISTUP(i),MOTHUP(1,i),MOTHUP(2,i),
     &                    ICOLUP(1,i),ICOLUP(2,i),(PUP(j,i),j=1,5),
     &                    VTIMUP(i),SPINUP(i)
            ENDDO

c giuseppe ILC 26/06/2007
            IF(i_coll.EQ.1 .OR. i_coll.EQ.2)THEN ! hadron colliders
* PDF information

c giuseppe 19/12/2006: exchange momentum fractions of initial state
c partons and related pdf's if initial state momenta have been exchanged
              IF(iflip.LT.0)THEN
                DO i=1,maxdistfun
                  idaux=id1_LH(i)
                  id1_LH(i)=id2_LH(i)
                  id2_LH(i)=idaux
                  xpdfaux=xpdf1_LH(i)
                  xpdf1_LH(i)=xpdf2_LH(i)
                  xpdf2_LH(i)=xpdfaux
                ENDDO
                baux=b1
                b1=b2
                b2=baux
              ENDIF
c end giuseppe 19/12/2006

* Revert to PDG numbering for gluon
              IF(id1_LH(i_selected).EQ.0) id1_LH(i_selected)=21
              IF(id2_LH(i_selected).EQ.0) id2_LH(i_selected)=21
              WRITE(23,'(A5,2(I4,1x),5(E13.6,1x))')
     &             '#pdf ',id1_LH(i_selected),id2_LH(i_selected),
     &             b1,b2,PDFScale,
     &             xpdf1_LH(i_selected),xpdf2_LH(i_selected)

c tetamp
c              WRITE(23,*) '#ampmodsquared', amod
c testampend

              WRITE(23,'(A)')'</event>'

            ELSEIF(i_coll.EQ.3)THEN ! e+e- collider
c* Beamstrahlung information (only if active)
c              IF(i_beamstrahlung.EQ.1)THEN
c                WRITE(23,'(A5,2(I4,1x),5(E13.6,1x))')
c     &               '#beamstrahlung ',id1_LH(1),id2_LH(1),
c     &               x1beam,x2beam,ecoll
c              ENDIF
c              IF(i_isr.EQ.1)THEN
c*     ISR information (only if active)
c                IF(i_beamstrahlung.EQ.0)THEN
c                  x1isr=b1
c                  x2isr=b2
c                ELSEIF(i_beamstrahlung.EQ.1)THEN
c                  x1isr=bb1
c                  x2isr=bb2
c                ENDIF
c                WRITE(23,'(A5,2(I4,1x),5(E13.6,1x))')
c     &               '#isr ',id1_LH(1),id2_LH(1),
c     &               x1isr,x2isr,QISRscale,
c     &               xpdf1_LH(1),xpdf2_LH(1)
c              ENDIF
              WRITE(23,'(A)')'</event>'

            ENDIF               ! IF(i_coll.EQ.1 .OR. i_coll.EQ.2)THEN
c end giuseppe ILC 26/06/2007
          ENDIF                 ! IF(iwrite_event.EQ.1)THEN
        ENDIF                   ! ENDIF EVENT ACCEPTED
      ENDIF                     ! ENDIF iflat.EQ.1 .OR. ioneshot.EQ.1

****  numero di chiamate effettive
      ncall_eff=ncall_eff+1

99999 FORMAT('IDUP  ',8I4)
99998 FORMAT('ISTUP ',8I4)
99997 FORMAT(4E17.9)
99996 FORMAT('ICOLUP ',16I4)
99995 FORMAT('SCALUP ',E17.9,' AQEDUP ',E17.9,' AQCDUP ',E17.9)

      RETURN
      END



      SUBROUTINE ZBOOST(gcm,bcm,q1,q2)
      IMPLICIT NONE
      REAL*8 gcm,bcm,q1(0:3),q2(0:3)
      q2(0)=gcm*(q1(0)+bcm*q1(3))
      q2(1)=q1(1)
      q2(2)=q1(2)
      q2(3)=q1(3)/gcm+bcm*q2(0)
      RETURN
      END


