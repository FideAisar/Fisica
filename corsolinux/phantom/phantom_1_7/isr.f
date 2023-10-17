************************************************************************
* isr.f
*
* Last update: Jul 06, 2007
*
* This file contains the implementation of the structure function
* ffISR for Initial State Radiation (ISR).
************************************************************************

      real*8 function ffISR(xx,Qscale)

      implicit real*8 (a-h,o-z)

* implementation of electron/positron structure function:
c QED Leading Log ISR structure function O(alpha_EM**3)
c Based on results by the "WW cross sections and distributions" working 
c group, CERN Workshop "Physics at LEP2" (1994/1995)  [hep-ph/9602351]

      parameter(nterm=10)       ! Bernoulli
      external rli2             ! dilogarithm function
      parameter(rme=0.51099906d-3,alpha_me=1.d0/137.0359895d0,
     &     ge=0.5772156649d0,pi=3.141592653589793238462643d0)


      s=QScale**2
      rl=dlog(s/rme/rme)
      betah=alpha_me*(rl-1.d0)/pi
      gamma=exp(gammln(1.d0+betah))
      rlim=exp(betah*(0.75d0-ge))*betah/gamma

* Implementation of the structure function.
c Note: the structure function is multiplied by a factor
c (1.d0/betah)*(1.d0-xx)**(1.d0-betah) which represents the jacobian
c associated to the mapping: x_vegas = (1.d0 - x_map)**betah

      CALL Bernoulli(nterm)
      if(xx.eq.1.d0)then
        ffISR = rlim
      else
        ffISR = rlim + (1.d0-xx)**(1.d0-betah)*(
c O(betah) contribution
     &       -betah*(1.d0+xx)/2.d0
c O(betah**2) contribution
     &       +(betah**2/8.d0)*(-4.d0*(1.d0+xx)*log(1.d0-xx)+3.d0
     &       *(1.d0+xx)*log(xx)-4.d0*log(xx)/(1.
     &       d0-xx)-5.d0-xx)
c O(betah**3) contribution
     &       -(betah**3/48.d0)*((1+xx)*(6.d0*rli2(xx,nterm)+
     &       12.d0*log(1-xx)*log(1-xx)-3.d0*pi**2)+(1.d0/(1.d0-xx))
     &       *(1.5d0*(1.d0+8.d0*xx+3.d0*xx**2)*log(xx)+
     &       6.d0*(xx+5)*(1-xx)*log(1-xx)+
     &       12.d0*(1+xx**2)*log(xx)*log(1-xx)-
     &       0.5d0*(1.d0+7.d0*xx**2)*log(xx)*log(xx)+
     &       0.25d0*(39.d0-24.d0*xx-15.d0*xx**2)))
     &       )
      endif

      ffISR=ffISR/betah

      return
      end



      subroutine map_ISR(x1,x2,Q2scale,xx1,xx2)

      implicit real*8 (a-h,o-z)
      
      parameter(rme=0.51099906d-3,alpha_me=1.d0/137.0359895d0,
     &     pi=3.141592653589793238462643d0)

* implementation of mapping
      rl=dlog(Q2scale/rme/rme)
      betah=alpha_me*(rl-1.d0)/pi

      xx1=1.d0-x1**(1.d0/betah)
      xx2=1.d0-x2**(1.d0/betah)

      return
      end




      double precision FUNCTION gammln(xx)

      IMPLICIT NONE
      REAL*8 xx
      INTEGER j
      double precision ser,stp,tmp,xa,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      xa=xx
      y=xa
      tmp=xa+5.5d0
      tmp=(xa+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      DO 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
   11 CONTINUE
      gammln=tmp+log(stp*ser/xa)
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #>,1')5c).
 
