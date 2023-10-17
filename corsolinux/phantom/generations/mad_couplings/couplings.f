C
C-----------------------------------------------------------------------------
C
      double precision function alfa(alfa0,qsq )
C
C-----------------------------------------------------------------------------
C
C	This function returns the 1-loop value of alpha.
C
C	INPUT: 
C		qsq   = Q^2
C
C-----------------------------------------------------------------------------
C
      implicit none
      double precision  qsq,alfa0
c
c constants
c
      double precision  One, Three, Pi,zmass
      parameter( One = 1.0d0, Three = 3.0d0 )
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( zmass = 91.188d0 )
cc
      alfa = alfa0 / ( 1.0d0 - alfa0*dlog( qsq/zmass**2 ) /Three /Pi )
ccc
      return
      end


      subroutine ewcoup
c****************************************************************     
c     EW PARAMETERS CALCULATIONS 
c****************************************************************
      implicit none
c
c constants
c
      double complex  ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
      double precision  Zero, One, Two, Three, Four, Half, Rt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
      parameter( Three = 3.0d0, Four = 4.0d0, Half = 0.5d0 )
      parameter( Rt2   = 1.414213562d0 )
      double precision  Pi, Fourpi
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Fourpi = Four * Pi )
c
c     include file
c     
      include 'coupl.inc'
      include 'input.inc'


      if ( idef.eq.0 ) then   
c------------------------------------------------------------
c     option=0 : MadEvent default (= AlpGen with iewopt=2) 
c------------------------------------------------------------

c-- equal to the input values
         sin2w  = hss   
         alpha  = hale  
         zmass  = zm    
c-- derived
         wmass  = zmass * sqrt( One - sin2w )
         gfermi = alpha * Fourpi/sin2w/(8d0*wmass**2/Rt2)     

      else if ( idef.eq.1 ) then
c-----------------------------------------------------
c     option=1 : LUSIFER and AlpGen (iewopt=3) default 
c-----------------------------------------------------

c-- equal to the input values
         zmass  = zm
         wmass  = wm
         gfermi = hgf
c-- derived
         sin2w  = One-(wmass/zmass)**2
         alpha  = Rt2*gfermi*wmass**2*sin2w/pi

      else if ( idef.eq.2 ) then
c-------------------------------------------------------------------
c     option=2 : W and Z mass are derived from couplings
c-------------------------------------------------------------------

c-- equal to the input values
         gfermi = hgf
         alpha  = hale
         sin2w  = hss
c-- derived
         wmass  = dsqrt(alpha*pi/sin2w/gfermi/Rt2)
         zmass  = wmass/dsqrt(One-sin2w)  
         

      else if ( idef.eq.3 ) then
c-----------------------------------------------------------------
c     option=3 : USER choice : you should know what you're doing!!
c-----------------------------------------------------------------

c-- equal to the input values (remember: only 3 of them are independent!)
         gfermi = hgf
         alpha  = hale   
         sin2w  = hss    
         wmass  = wm     
         zmass  = zm     
c-- derived

c--example:MCFM parameters,
c  rho=1+3d0*alpha/16d0/pi/sin2w*(mtop/wm)**2, 
c  wmass=mz*dsqrt(rho)*cosw
c  mtop= 159 GeV
c--uncomment the following lines 
c         hgf    = 1.16638d-5
c         hale   = 1/128.9d0 
c         wm     = 80.41d0     
c         zm     = 91.187d0     
c
c         gfermi = hgf
c         alpha  = hale   
c         sin2w  = fourpi*alpha/(8d0*wm**2*hgf/rt2)
c         wmass  = wm     
c         zmass  = zm     

         

      else
c---------------------------
c else error output and halt
c---------------------------
         write(6,*) '*** EWcoup called w/ invalid option - halting ***'
         stop
      end if

      return
      end

      
C
C-----------------------------------------------------------------------------
C
      double precision function alfaw( qsq,nh )
C
C-----------------------------------------------------------------------------
C
C	This function returns the 1-loop value of alpha_w.
C
C	INPUT: 
C		qsq = Q^2
C               nh  = # of Higgs doublets
C
C-----------------------------------------------------------------------------
C
      implicit none
      double precision  qsq, alphaw, dum
      integer  nh, nq
      include 'coupl.inc'
c
c constants
c
      double precision  Two, Four, Pi, Twpi
      parameter( Two = 2.0d0, Four = 4.0d0 )
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Twpi = 3.0d0*Four*Pi )
cc
      if ( qsq.ge.tmass**2 ) then
         nq = 6
      else
         nq = 5
      end if
      alphaw = gw**2 / Four / Pi
      dum = ( 22.0d0 - Four*nq - nh/Two ) / Twpi
      alfaw = alphaw / ( 1.0d0 + dum*alphaw*dlog( qsq/zmass**2 ) )
ccc
      return
      end

C-----------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION ALPHAS(Q,AMZ,NLOOP)
c
c     Evaluation of strong coupling constant alpha_S
c     Author: R.K. Ellis
c
c     q -- scale at which alpha_s is to be evaluated
c     amz -- value of alpha_s at the mass of the Z-boson
c     nloop -- the number of loops (1,2, or 3) at which beta 
c     function is evaluated to determine running.
c     the values of the cmass and the bmass should be set
c     in common block qmass.
C-----------------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION Q,T,AMZ,AMZ0,AMB,AMC
      DOUBLE PRECISION AS_OUT
      INTEGER NLOOP,NLOOP0,NF3,NF4,NF5
      PARAMETER(NF5=5,NF4=4,NF3=3)
C
      REAL*8       CMASS,BMASS
      COMMON/QMASS/CMASS,BMASS
      DATA CMASS,BMASS/1.42D0,4.7D0/  ! HEAVY QUARK MASSES FOR THRESHOLDS
C
      REAL*8 ZMASS
      DATA ZMASS/91.188D0/
C
      SAVE AMZ0,NLOOP0,AMB,AMC
      DATA AMZ0,NLOOP0/0D0,0/
      IF (Q .LE. 0D0) THEN 
         WRITE(6,*) 'q .le. 0 in alphas'
         WRITE(6,*) 'q= ',Q
         PAUSE
      ENDIF
      IF (AMZ .LE. 0D0) THEN 
         WRITE(6,*) 'amz .le. 0 in alphas',AMZ
         WRITE(6,*) 'continue with amz=0.1185?'
         PAUSE
         AMZ=0.1185D0
      ENDIF
      IF (CMASS .LE. 0.3D0) THEN 
         WRITE(6,*) 'cmass .le. 0.3GeV in alphas',CMASS
         WRITE(6,*) 'continue with cmass=1.5GeV?'
         PAUSE
         CMASS=1.42D0
      ENDIF
      IF (BMASS .LE. 0D0) THEN 
         WRITE(6,*) 'bmass .le. 0 in alphas',BMASS
         WRITE(6,*) 'COMMON/QMASS/CMASS,BMASS'
         PAUSE
         BMASS=4.7D0
      ENDIF
c--- establish value of coupling at b- and c-mass and save
      IF ((AMZ .NE. AMZ0) .OR. (NLOOP .NE. NLOOP0)) THEN
         AMZ0=AMZ
         NLOOP0=NLOOP
         T=2D0*DLOG(BMASS/ZMASS)
         CALL NEWTON1(T,AMZ,AMB,NLOOP,NF5)
         T=2D0*DLOG(CMASS/BMASS)
         CALL NEWTON1(T,AMB,AMC,NLOOP,NF4)
      ENDIF

c--- evaluate strong coupling at scale q
      IF (Q  .LT. BMASS) THEN
           IF (Q  .LT. CMASS) THEN
             T=2D0*DLOG(Q/CMASS)
             CALL NEWTON1(T,AMC,AS_OUT,NLOOP,NF3)
           ELSE
             T=2D0*DLOG(Q/BMASS)
             CALL NEWTON1(T,AMB,AS_OUT,NLOOP,NF4)
           ENDIF
      ELSE
      T=2D0*DLOG(Q/ZMASS)
      CALL NEWTON1(T,AMZ,AS_OUT,NLOOP,NF5)
      ENDIF
      ALPHAS=AS_OUT
      RETURN
      END


      SUBROUTINE DIFF(Q,AMZ,NLOOP)
      IMPLICIT NONE
      DOUBLE PRECISION BETA(3:5),B0(3:5),C1(3:5),C2(3:5)
      INTEGER NLOOP,J
      DOUBLE PRECISION Q,QP,QM,AMZ
      DOUBLE PRECISION X1,X2,X3,EP,DIFF1,ALPHAS
      REAL*8       CMASS,BMASS
      COMMON/QMASS/CMASS,BMASS	
C---     B0=(11.-2.*F/3.)/4./PI
      DATA B0/0.716197243913527D0,0.66314559621623D0,0.61009394851893D0/
C---     C1=(102.D0-38.D0/3.D0*F)/4.D0/PI/(11.D0-2.D0/3.D0*F)
      DATA C1/.565884242104515D0,0.49019722472304D0,0.40134724779695D0/
C---     C2=(2857.D0/2.D0-5033*F/18.D0+325*F**2/54)
C---     /16.D0/PI**2/(11.D0-2.D0/3.D0*F)
      DATA C2/0.453013579178645D0,0.30879037953664D0,0.14942733137107D0/
C---     DEL=SQRT(4*C2-C1**2)

      X1=ALPHAS(Q,AMZ,1)
      X2=ALPHAS(Q,AMZ,2)
      X3=ALPHAS(Q,AMZ,3)
      J=3
      IF (Q .GT. CMASS) J=4
      IF (Q .GT. BMASS) J=5
      EP=.001D0
      QP=Q*(1.D0+EP)
      QM=Q*(1.D0-EP)
      IF (NLOOP .EQ.1) THEN 
      BETA(J)=-B0(J)*X1**2
      DIFF1=(ALPHAS(QP,AMZ,1)-ALPHAS(QM,AMZ,1))/4d0/EP/BETA(J)
      ENDIF
      IF (NLOOP .EQ.2) THEN 
      BETA(J)=-B0(J)*X2**2*(1.D0+C1(J)*X2)
      DIFF1=(ALPHAS(QP,AMZ,2)-ALPHAS(QM,AMZ,2))/4d0/EP/BETA(J)
      ENDIF
      IF (NLOOP .EQ.3) THEN 
      BETA(J)=-B0(J)*X3**2*(1.D0+C1(J)*X3+C2(J)*X3**2)
      DIFF1=(ALPHAS(QP,AMZ,3)-ALPHAS(QM,AMZ,3))/4d0/EP/BETA(J)
      ENDIF
      WRITE(6,*) Q,DIFF1,NLOOP
      RETURN
      END

      SUBROUTINE NEWTON1(T,A_IN,A_OUT,NLOOP,NF)
C     Author: R.K. Ellis

c---  calculate a_out using nloop beta-function evolution 
c---  with nf flavours, given starting value as-in
c---  given as_in and logarithmic separation between 
c---  input scale and output scale t.
c---  Evolution is performed using Newton's method,
c---  with a precision given by tol.

      IMPLICIT NONE
      INTEGER NLOOP,NF
      REAL*8 T,A_IN,A_OUT,AS,TOL,F2,F3,F,FP,DELTA
      REAL*8 B0(3:5),C1(3:5),C2(3:5),DEL(3:5)
      PARAMETER(TOL=5.D-4)
C---     B0=(11.-2.*NF/3.)/4./PI
      DATA B0/0.716197243913527D0,0.66314559621623D0,0.61009394851893D0/
C---     C1=(102.D0-38.D0/3.D0*NF)/4.D0/PI/(11.D0-2.D0/3.D0*NF)
      DATA C1/.565884242104515D0,0.49019722472304D0,0.40134724779695D0/
C---     C2=(2857.D0/2.D0-5033*NF/18.D0+325*NF**2/54)
C---     /16.D0/PI**2/(11.D0-2.D0/3.D0*NF)
      DATA C2/0.453013579178645D0,0.30879037953664D0,0.14942733137107D0/
C---     DEL=SQRT(4*C2-C1**2)
      DATA DEL/1.22140465909230D0,0.99743079911360D0,0.66077962451190D0/
      F2(AS)=1D0/AS+C1(NF)*LOG((C1(NF)*AS)/(1D0+C1(NF)*AS))
      F3(AS)=1D0/AS+0.5D0*C1(NF)
     & *LOG((C2(NF)*AS**2)/(1D0+C1(NF)*AS+C2(NF)*AS**2))
     & -(C1(NF)**2-2D0*C2(NF))/DEL(NF)
     & *ATAN((2D0*C2(NF)*AS+C1(NF))/DEL(NF))

           
      A_OUT=A_IN/(1D0+A_IN*B0(NF)*T)
      IF (NLOOP .EQ. 1) RETURN
      A_OUT=A_IN/(1D0+B0(NF)*A_IN*T+C1(NF)*A_IN*LOG(1D0+A_IN*B0(NF)*T))
      IF (A_OUT .LT. 0D0) AS=0.3D0
 30   AS=A_OUT

      IF (NLOOP .EQ. 2) THEN
      F=B0(NF)*T+F2(A_IN)-F2(AS)
      FP=1D0/(AS**2*(1D0+C1(NF)*AS))
      ENDIF
      IF (NLOOP .EQ. 3) THEN
      F=B0(NF)*T+F3(A_IN)-F3(AS)
      FP=1D0/(AS**2*(1D0+C1(NF)*AS+C2(NF)*AS**2))
      ENDIF
      A_OUT=AS-F/FP
      DELTA=ABS(F/FP/AS)
      IF (DELTA .GT. TOL) GO TO 30
      RETURN
      END


C-----------------------------------------------------------------------------
C
      double precision function mfrun(mf,scale,asmz,nloop)
C
C-----------------------------------------------------------------------------
C
C	This function returns the 2-loop value of a MSbar fermion mass
C       at a given scale.
C
C	INPUT: mf    = MSbar mass of fermion at MSbar fermion mass scale 
C	       scale = scale at which the running mass is evaluated
C	       asmz  = AS(MZ) : this is passed to alphas(scale,asmz,nloop)
C              nloop = # of loops in the evolution
C       
C
C
C	EXTERNAL:      double precision alphas(scale,asmz,nloop)
C                      
C-----------------------------------------------------------------------------
C
      implicit none
C
C     ARGUMENTS
C
      double precision  mf,scale,asmz
      integer           nloop
C
C     LOCAL
C
      double precision  beta0, beta1,gamma0,gamma1
      double precision  A1,as,asmf,l2
      integer  nf
C
C     EXTERNAL
C
      double precision  alphas
      external          alphas
C
C     COMMON
C
      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
c
c     CONSTANTS
c
      double precision  One, Two, Three, Pi
      parameter( One = 1.0d0, Two = 2.0d0, Three = 3.0d0 )
      parameter( Pi = 3.14159265358979323846d0) 
cc
C
C
      if ( mf.gt.tmass ) then
         nf = 6
      else
         nf = 5
      end if

      beta0 = ( 11.0d0 - Two/Three *nf )/4d0
      beta1 = ( 102d0  - 38d0/Three*nf )/16d0
      gamma0= 1d0
      gamma1= ( 202d0/3d0  - 20d0/9d0*nf )/16d0
      A1    = -beta1*gamma0/beta0**2+gamma1/beta0
      as    = alphas(scale,asmz,nloop)
      asmf  = alphas(mf   ,asmz,nloop)
      l2    = (1+ A1*as/Pi)/(1+ A1*asmf/Pi)
      
      
      mfrun = mf * (as/asmf)**(gamma0/beta0)

      if(nloop.eq.2) mfrun =mfrun*l2
ccc
      return
      end

      SUBROUTINE TOPWID(RMT,RMW,RMB,RGW,GW,RGT)
c*************************************************************************
c     THE TOTAL WEAK DECAY WIDTH OF THE TOP QUARK, INCLUDING
c     THE EFFECTS OF BOTTOM MASS AND, IF IGW=1,  A FINITE W WIDTH.
c     From James Stirling 6-10-94
c
c     RMT=TOP MASS
c     RMW=W   MASS
c     RMB=B   MASS
c     RGW=W   WIDTH
c     GW =WEAK COUPLING
c
c     RGT=TOP WIDTH
c
c*************************************************************************
      IMPLICIT COMPLEX*16(A-H,O-Z)
      REAL*8 RMT,RMB,RMW,XW,XB,RGW,RGT,GW
*
      PI=4.*DATAN(1.D0)
      XGW=dcmplx(GW/2d0/dsqrt(2d0))
*
      XB=RMB/RMT
      XW=RMW/RMT
*
      RM=XB**2
      OM=1.+RM-DCMPLX(RMW**2,RMW*RGW)/RMT**2
      Y1=OM+CDSQRT(OM*OM-4.*RM)
      Y0=OM-CDSQRT(OM*OM-4.*RM)
      Z1=2.
      Z0=2.*CDSQRT(RM)
*
      D0=(-Y0**8+3.*Y0**7*RM+3.*Y0**7-8.*Y0**6*RM-12.*Y0**5*RM**
     . 2-12.*Y0**5*RM+96.*Y0**4*RM**2-48.*Y0**3*RM**3-48.*Y0**3*
     . RM**2-128.*Y0**2*RM**3+192.*Y0*RM**4+192.*Y0*RM**3-256.*
     . RM**4)/(24.*Y0**4*(Y1-Y0))
      D1=(-Y1**8+3.*Y1**7*RM+3.*Y1**7-8.*Y1**6*RM-12.*Y1**5*RM**
     . 2-12.*Y1**5*RM+96.*Y1**4*RM**2-48.*Y1**3*RM**3-48.*Y1**3*
     . RM**2-128.*Y1**2*RM**3+192.*Y1*RM**4+192.*Y1*RM**3-256.*
     . RM**4)/(24.*Y1**4*(Y1-Y0))
      A4=(32.*RM**4*(Y1-Y0))/(3.*Y1*Y0*(Y1-Y0))
      A3=(8.*RM**3*(-3.*Y1**2*Y0*RM-3.*Y1**2*Y0+4.*Y1**2*RM+3.*
     . Y1*Y0**2*RM+3.*Y1*Y0**2-4.*Y0**2*RM))/(3.*Y1**2*Y0**2*(Y1
     . -Y0))
      A2=(8.*RM**3*(2.*Y1**3*Y0**2-3.*Y1**3*Y0*RM-3.*Y1**3*Y0+4.
     . *Y1**3*RM-2.*Y1**2*Y0**3+3.*Y1*Y0**3*RM+3.*Y1*Y0**3-4.*Y0
     . **3*RM))/(3.*Y1**3*Y0**3*(Y1-Y0))
      A1=(2.*RM**2*(3.*Y1**4*Y0**3*RM+3.*Y1**4*Y0**3+8.*Y1**4*Y0
     . **2*RM-12.*Y1**4*Y0*RM**2-12.*Y1**4*Y0*RM+16.*Y1**4*RM**2
     . -3.*Y1**3*Y0**4*RM-3.*Y1**3*Y0**4-8.*Y1**2*Y0**4*RM+12.*
     . Y1*Y0**4*RM**2+12.*Y1*Y0**4*RM-16.*Y0**4*RM**2))/(3.*Y1**
     . 4*Y0**4*(Y1-Y0))
      B0=(Y1**3-3.*Y1**2*RM-3.*Y1**2+8.*Y1*RM-Y0**3+3.*Y0**2*RM+
     . 3.*Y0**2-8.*Y0*RM)/(24.*(Y1-Y0))
      B1=(Y1+Y0-3.*RM-3.)/24.
      B2=1./24.
*
      RINT=D0*CDLOG((Z1-Y0)/(Z0-Y0))
     .    -D1*CDLOG((Y1-Z1)/(Y1-Z0))
     .    -A4/3.*(1./Z1**3-1./Z0**3)
     .    -A3/2.*(1./Z1**2-1./Z0**2)
     .    -A2   *(1./Z1   -1./Z0   )
     .    +A1*CDLOG(Z1/Z0)
     .    +B0   *(Z1   -Z0   )
     .    +B1/2.*(Z1**2-Z0**2)
     .    +B2/3.*(Z1**3-Z0**3)
*
      XGW4=XGW**4
*
* TOTAL WIDTH INCLUDES FLAVOUR & COLOUR FACTORS
      RGT=RMT**3/(RMW*RGW)*XGW4/(8.*PI**3)*DIMAG(RINT)
      RGT=9.*RGT
      RETURN
      END
