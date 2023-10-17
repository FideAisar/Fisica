      PROGRAM A_S
c     Interactive interface to ALPHAS(Q,AMZ,NLOOP) by R.K. Ellis
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
      DOUBLE PRECISION Q,AMZ,AS,ALPHAS,Q2,ALFA
      INTEGER NLOOP

      WRITE(*,*) 'SCALE ? (GeV)'
      READ(*,*) Q
      WRITE(*,*) 'ALPHA_S at MZ ?'
      READ(*,*) AMZ
      WRITE(*,*) 'NLOOP ?'
      READ(*,*) NLOOP
      AS=ALPHAS(Q,AMZ,NLOOP)
      WRITE(*,*)'ALPHA_S at',Q,': ',AS,' with Ellis routine'
      Q2=Q*Q
      AS=alfa(AMZ,Q2)
      WRITE(*,*)'ALPHA_S at',Q,': ',AS,' with basic routine'
      END
      
