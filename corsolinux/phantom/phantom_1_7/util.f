      subroutine boost(q,pboost,qprime)
      implicit real*8 (a-h,o-z)

c this subroutine performs the boost according to the vector pboost of 
c  the fourvector q to the fourvector qprime. The vector q and the 
c   resulting qprime might also have the same name in the call       
c  Let's call rmboost the mass of pboost. With this transformation
c  if q= (rmboost,0,0,0), qprime=pboost 
    

      dimension q(0:3), pboost(0:3), qprime (0:3)

      rmboost=pboost(0)**2-pboost(1)**2-pboost(2)**2-pboost(3)**2
      rmboost=sqrt(max(rmboost,0.d0))
      aux=(q(0)*pboost(0)+q(1)*pboost(1)+q(2)*pboost(2)+q(3)*pboost(3))
     &      /rmboost
      aaux=(aux+q(0))/(pboost(0)+rmboost)
      qprime(0)=aux
      qprime(1)=q(1)+aaux*pboost(1)
      qprime(2)=q(2)+aaux*pboost(2)
      qprime(3)=q(3)+aaux*pboost(3)

      return
      end

      subroutine boostinv(q,pboost,qprime)
      implicit real*8 (a-h,o-z)

c this subroutine performs the inverse boost according to the vector 
c  pboost  ( which corresponds to the boost according to pboost(0), 
c  -pboost(i) ) of the fourvector q to the fourvector qprime. 
c  The vector q and the resulting qprime might also have the same 
c  name in the call       
c  Let's call rmboost the mass of pboost. With this transformation
c  if q=pboost, qprime= (rmboost,0,0,0),   

      dimension q(0:3), pboost(0:3), qprime (0:3)

      rmboost=pboost(0)**2-pboost(1)**2-pboost(2)**2-pboost(3)**2
      rmboost=sqrt(max(rmboost,0.d0))
      aux=(q(0)*pboost(0)-q(1)*pboost(1)-q(2)*pboost(2)-q(3)*pboost(3))
     &      /rmboost
      aaux=(aux+q(0))/(pboost(0)+rmboost)
      qprime(0)=aux
      qprime(1)=q(1)-aaux*pboost(1)
      qprime(2)=q(2)-aaux*pboost(2)
      qprime(3)=q(3)-aaux*pboost(3)

      return
      end




      subroutine rot(p,q,pp)

* the subroutine performs a rotation such that the 3-vector p goes into
* the 3-vector pp. The rotation is such that the vector q', of the same
*  lenght as q but along the z axis, becomes q after the rotation
      
      implicit real*8 (a-h,o-z)

      dimension pp(3), p(3), q(3) 


      qmod=sqrt(q(1)**2+q(2)**2+q(3)**2)
      if (qmod.eq.0.d0) then
        print*, ' ERROR in subroutine rot '
        print*, ' spatial q components are 0.d0 ! '
        stop
      endif

      cth=q(3)/qmod
      sth=1.d0-cth**2
      if (sth.eq.0.d0) then
        pp(1)=p(1)
        pp(2)=p(2)
        pp(3)=p(3)
        return
      endif
      sth=sqrt(sth)

      qmodt=sqrt(q(1)**2+q(2)**2)
      if (qmodt.eq.0.d0) then
        pp(1)=p(1)
        pp(2)=p(2)
        pp(3)=p(3)
        return
      endif
        
      cfi=q(1)/qmodt
      sfi=q(2)/qmodt

c  avoid possible problems if p and pp are the same vector:
      
      p1=p(1)
      p2=p(2)
      p3=p(3)

c  perform the rotation

      pp(1)=cth*cfi*p1-sfi*p2+sth*cfi*p3
      pp(2)=cth*sfi*p1+cfi*p2+sth*sfi*p3
      pp(3)=-sth   *p1       +cth*    p3


      return
      end

      subroutine rotinv(p,q,pp)
* the subroutine performs a rotation such that the 3-vector p goes into
* the 3-vector pp. The rotation is such that the vector q goes to the z axis.
* It performs first a rotation around z which brings the vectot in the xz-plane and
* then a rotaion around y which brings the vector along the z axis
      
      implicit real*8 (a-h,o-z)
      dimension pp(3), p(3), q(3) 

      qmodt=q(1)**2+q(2)**2
      qmod=qmodt+q(3)**2
      qmodt=sqrt(qmodt)
      qmod=sqrt(qmod)

      if (qmod.eq.0.d0) then
        print*, ' ERROR in subroutine rot '
        print*, ' spatial q components are 0.d0 ! '
        stop
      endif

      cth=q(3)/qmod
      sth=1.d0-cth**2
      if (sth.eq.0.d0) then
        pp(1)=p(1)
        pp(2)=p(2)
        pp(3)=p(3)
        return
      endif
      sth=sqrt(sth)

      if (qmodt.eq.0.d0) then
        pp(1)=p(1)
        pp(2)=p(2)
        pp(3)=p(3)
        return
      endif
        
      cfi=q(1)/qmodt
      sfi=q(2)/qmodt

c  avoid possible problems if p and pp are the same vector:     
      p1=p(1)
      p2=p(2)
      p3=p(3)

c  perform the rotation
      pp(1)=cth*cfi*p1+cth*sfi*p2-sth*p3
      pp(2)=-sfi*p1   +cfi*p2
      pp(3)=+sth*cfi*p1+sth*sfi*p2+cth*    p3

      return
      end


      
*****************************************************************************
* Auxiliary subroutine
* INPUT:
*       pbst        : intermediate particle momentum
*       shat        : the corresponding mass (assumed positive)
*       pmod        : the absolute value of the 3-momentum of the two
*                     daughter particles in the decay CM
*       xm1sq,xm2sq : the squared masses of the daughter particles
*       x1,x2       : two random numbers uniformly distributed in [0,1]
* OUTPUT:
*       p1,p2       : the momenta of the two daughters in the frame in
*                     which the intermediate particle has momentum pbst 
*                             p1 + p2 = pbst
*
* If pbst=(sqrt(shat),0,0,0) the boost reduces to the identity
* Ezio May 2 2007
* A rotation of the generated momentum has been added in such a way to assume as
* z-axis the boost momentum
* 
*****************************************************************************
      SUBROUTINE TWOBODY(pbst,shat,pmod,xm1sq,xm2sq,x1,x2,p1,p2)
      IMPLICIT NONE
      REAL*8 pbst(0:3),shat,pmod,xm1sq,xm2sq,x1,x2,p1(0:3),p2(0:3)
* Local variables
      REAL*8 app,costh,sinth,ph,rmboost,aapp,pi,twopi
      INTEGER i
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi)
       
      costh=min(x1*2.d0-1.d0,1.d0)
      costh=max(costh,-1.d0)
      sinth=sqrt(1.d0-costh*costh)
      ph=x2*twopi
* p1
      p1(0)=sqrt(pmod*pmod+xm1sq)
      p1(1)=pmod*sinth*cos(ph)
      p1(2)=pmod*sinth*sin(ph)
      p1(3)=pmod*costh
* rotate in such a way that angles are relative to the direction of pbst
      call rot(p1(1),pbst(1),p1(1))
* boost back
      rmboost=sqrt(shat)
      app=(p1(0)*pbst(0)+p1(1)*pbst(1)+p1(2)*pbst(2)+p1(3)*pbst(3))
     &      /rmboost
      aapp=(app+p1(0))/(pbst(0)+rmboost)
      p1(0)=app
      p1(1)=p1(1)+aapp*pbst(1)
      p1(2)=p1(2)+aapp*pbst(2)
      p1(3)=p1(3)+aapp*pbst(3)
      DO i=0,3
        p2(i)=pbst(i)-p1(i)
      ENDDO
      
      RETURN
      END
      


*****************************************************************************
* Auxiliary subroutine. Defines extrema for flat+Breit-Wigner+flat mapping
* of resonances. Checks for position of estrsup,estrinf with respect to
* spm1lo,spm1hi
* INPUT:
*   estrsup,estrinf: upper and lower limit for the generated mass
*   spm1lo,spm1hi  : lower and upper limit of the Breit-Wigner region
*   rmx1,gamx1       : mass and width of the resonance
* OUTPUT:
*  alfa1,alfa2       : separation points in [0,1] between the three regions
*                      Can also be 0 or 1. Some care necessary in this cases
*  auxalfa2,xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max
*                    : auxiliary quantities used for generating the mass
*                      and computing the jacobian
*****************************************************************************
      SUBROUTINE PHSP_INI(estrsup,estrinf,spm1lo,spm1hi,
     &                       rmx1,gamx1,
     &                           alfa1,alfa2,auxalfa2,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max)
     
      IMPLICIT REAL*8(a-h,o-z)
      
      IF(estrsup.LE.spm1lo)THEN ! only flat1
        alfa2=1.d0
        alfa1=1.d0
        auxalfa2=alfa2
        xf1_max=estrsup
        xf1_min=estrinf
      ELSEIF(estrsup.LE.spm1hi)THEN ! flat1 + W/Z ?
        alfa2=1.d0
        auxalfa2=alfa2
        IF((estrinf.GE.spm1lo))THEN  ! W/Z only
          alfa1=0.d0
          xw1_max=atan((estrsup**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((estrinf**2-rmx1**2)/(gamx1*rmx1))
        ELSE !  flat1 + W/Z !
          xw1_min=atan((spm1lo**2-rmx1**2)/(gamx1*rmx1))
          xw1_max=atan((estrsup**2-rmx1**2)/(gamx1*rmx1))

          rkappa2low=(xw1_max-xw1_min)*
     &          ((spm1lo**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          xf1_max=spm1lo
          xf1_min=estrinf
          rkappa1=(xf1_max-xf1_min)*2.d0*spm1lo
* If rkappa2low=0 then estrsup is numerically indistinguishable from spm1lo
          IF(rkappa2low.LE.0.d0)THEN    
            alfa1=1.d0
            auxalfa2=alfa2
            xf1_max=estrsup
            xf1_min=estrinf
	    RETURN
	  ENDIF
          rkappa=rkappa1/rkappa2low
          alfa1=rkappa/(rkappa+1.d0)
        ENDIF
      ELSE ! estrsup.GT.spm1hi: flat1 + W/Z + flat2 ?
        IF(estrinf.GE.spm1hi)THEN  ! only flat2
          alfa1=0.d0
          alfa2=0.d0
          auxalfa2=alfa2
          xf2_max=estrsup
          xf2_min=estrinf
        ELSEIF((estrinf.GE.spm1lo))THEN !  W/Z + flat2 
          alfa1=0.d0
          xf2_max=estrsup
          xf2_min=spm1hi
          rkappa3=(xf2_max-xf2_min)*2.d0*spm1hi
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((estrinf**2-rmx1**2)/(gamx1*rmx1))
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3
c          alfa2=rk2/(rk2+1.d0)
          alfa2=0.5d0
          auxalfa2=alfa2
        ELSE !  flat + W/Z + flat
          xf2_max=estrsup
          xf2_min=spm1hi
          rkappa3=(xf2_max-xf2_min)*2.d0*spm1hi
          xw1_min=atan((spm1lo**2-rmx1**2)/(gamx1*rmx1))
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rkappa2low=(xw1_max-xw1_min)*
     &          ((spm1lo**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3
c          auxalfa2=rk2/(rk2+1.d0)
          auxalfa2=0.5d0
c          auxalfa2=8.d0/9.d0
          xf1_min=estrinf
          xf1_max=spm1lo
          rkappa1=(xf1_max-xf1_min)*2.d0*spm1lo
          rk1=rkappa1/rkappa2low*auxalfa2
c          alfa1=rk1/(rk1+1.d0)
          alfa1=1.d0/3.d0
c          alfa1=1.d0/10.d0
          alfa2=auxalfa2*(1.d0-alfa1)+alfa1
        ENDIF
      ENDIF
      
      RETURN
      END

*****************************************************************************
* Auxiliary subroutine. Defines extrema for flat+Breit-Wigner+flat
*                                              +Breit-Wigner+flat mapping
* of resonances. Checks for position of estrsup,estrinf with respect to
* spm1lo,spm1hi
* INPUT:
*   estrsup,estrinf: upper and lower limit for the generated mass
*   spm1lo,spm1hi  : lower and upper limit of first Breit-Wigner region
*   rmx1,gamx1       : mass and width of first resonance
*   spm2lo,spm2hi  : lower and upper limit of second Breit-Wigner region
*   rmx2,gamx2       : mass and width of second resonance
* OUTPUT:
*  alfa1,alfa2,alfa3,alfa4: separation points in [0,1] between the five regions
*                      Can also be 0 or 1. Some care necessary in this cases
*  auxalfa2,xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,xw2_min,xw2_max,
*  xf3_min,xf3_max,auxalfa3,auxalfa4
*                    : auxiliary quantities used for generating the mass
*                      and computing the jacobian
*****************************************************************************
      SUBROUTINE PHSP_INI5(estrsup,estrinf,
     &                     spm1lo,spm1hi,rmx1,gamx1,
     &                     spm2lo,spm2hi,rmx2,gamx2,
     &                     alfa1,alfa2,auxalfa2,
     &                     alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max)
     
      IMPLICIT REAL*8(a-h,o-z)
      
      IF(estrsup.LE.spm1lo)THEN ! only flat1
        alfa1=1.d0
        alfa2=1.d0
        alfa3=1.d0
        alfa4=1.d0
        auxalfa4=alfa4
        auxalfa3=alfa3
        auxalfa2=alfa2
        xf1_max=estrsup
        xf1_min=estrinf
      ELSEIF(estrsup.LE.spm1hi)THEN ! flat1 + res1 ?
        alfa2=1.d0
        alfa3=1.d0
        alfa4=1.d0
        auxalfa4=alfa4
        auxalfa3=alfa3
        auxalfa2=alfa2
        IF((estrinf.GE.spm1lo))THEN  ! res1 only
          alfa1=0.d0
          xw1_max=atan((estrsup**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((estrinf**2-rmx1**2)/(gamx1*rmx1))
        ELSE !  flat1 + res1 !
          xw1_min=atan((spm1lo**2-rmx1**2)/(gamx1*rmx1))
          xw1_max=atan((estrsup**2-rmx1**2)/(gamx1*rmx1))
          rkappa2low=(xw1_max-xw1_min)*
     &          ((spm1lo**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          xf1_max=spm1lo
          xf1_min=estrinf
          rkappa1=(xf1_max-xf1_min)*2.d0*spm1lo
* If rkappa2low=0 then estrsup is numerically indistinguishable from spm1lo
          IF(rkappa2low.LE.0.d0)THEN    
            alfa1=1.d0
            xf1_max=estrsup
            xf1_min=estrinf
            RETURN
          ENDIF
          rk1=rkappa1/rkappa2low
          alfa1=rk1/(rk1+1.d0)
        ENDIF
      ELSEIF(estrsup.LE.spm2lo)THEN ! estrsup.GT.spm1hi: flat1 + res1 + flat2 ?
        alfa3=1.d0
        alfa4=1.d0
        auxalfa4=alfa4
        auxalfa3=alfa3
        IF(estrinf.GE.spm1hi)THEN  ! only flat2
          alfa1=0.d0
          alfa2=0.d0
          auxalfa2=alfa2
          xf2_max=estrsup
          xf2_min=estrinf
        ELSEIF((estrinf.GE.spm1lo))THEN !  res1 + flat2 
          alfa1=0.d0
          xf2_max=estrsup
          xf2_min=spm1hi
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((estrinf**2-rmx1**2)/(gamx1*rmx1))
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3low
          alfa2=rk2/(rk2+1.d0)
          auxalfa2=alfa2
        ELSE !  flat1 + res1 + flat2
          xf2_max=estrsup
          xf2_min=spm1hi
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          xw1_min=atan((spm1lo**2-rmx1**2)/(gamx1*rmx1))
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          rkappa2low=(xw1_max-xw1_min)*
     &          ((spm1lo**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3low
ccc          auxalfa2=rk2/(rk2+1.d0)
          auxalfa2=0.5d0
          xf1_min=estrinf
          xf1_max=spm1lo
          rkappa1=(xf1_max-xf1_min)*2.d0*spm1lo
          rk1=rkappa1/rkappa2low*auxalfa2
CCC          alfa1=rk1/(rk1+1.d0)
          alfa1=1.d0/3.d0
          alfa2=auxalfa2*(1.d0-alfa1)+alfa1
        ENDIF
      ELSEIF(estrsup.LE.spm2hi)THEN 
             ! estrsup.GT.spm2lo: flat1 + res1 + flat2 + res2?
        alfa4=1.d0
        auxalfa4=alfa4
        IF(estrinf.GE.spm2lo)THEN  ! only res2
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          auxalfa3=alfa3
          auxalfa2=alfa2
          xw2_max=atan((estrsup**2-rmx2**2)/(gamx2*rmx2))
          xw2_min=atan((estrinf**2-rmx2**2)/(gamx2*rmx2))
        ELSEIF(estrinf.GE.spm1hi)THEN  ! flat2 + res2
          alfa1=0.d0
          alfa2=0.d0
          auxalfa2=alfa2
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          xw2_max=atan((estrsup**2-rmx2**2)/(gamx2*rmx2))
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
* If rkappa4low=0 then estrsup is numerically indistinguishable from spm2lo
          IF(rkappa4low.LE.0.d0)THEN    
            alfa3=1.d0
            auxalfa3=alfa3
            xf2_max=estrsup
            xf2_min=estrinf
            RETURN
          ENDIF
          xf2_max=spm2lo
          xf2_min=estrinf
          rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
          rk3=rkappa3high/rkappa4low
          alfa3=rk3/(rk3+1.d0)
          auxalfa3=alfa3
        ELSEIF(estrinf.GE.spm1lo)THEN  ! res1 + flat2 + res2
          alfa1=0.d0
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          xw2_max=atan((estrsup**2-rmx2**2)/(gamx2*rmx2))
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
* If rkappa4low=0 then estrsup is numerically indistinguishable from spm2lo
          IF(rkappa4low.LE.0.d0)THEN    
            auxalfa3=1.d0
            xf2_max=estrsup
          ELSE
            xf2_max=spm2lo
            rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
            rk3=rkappa3high/rkappa4low
ccc            auxalfa3=rk3/(rk3+1.d0)
            auxalfa3=0.5d0
          ENDIF
          xf2_min=spm1hi
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((estrinf**2-rmx1**2)/(gamx1*rmx1))
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3low*auxalfa3
ccc          alfa2=rk2/(rk2+1.d0)
          alfa2=1.d0/3.d0
          auxalfa2=alfa2
          alfa3=auxalfa3*(1.d0-alfa2)+alfa2
          IF(rkappa4low.LE.0.d0) alfa3=1.d0    
        ELSE  ! flat1 + res1 + flat2 + res2
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          xw2_max=atan((estrsup**2-rmx2**2)/(gamx2*rmx2))
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
* If rkappa4low=0 then estrsup is numerically indistinguishable from spm2lo
          IF(rkappa4low.LE.0.d0)THEN    
            auxalfa3=1.d0
            xf2_max=estrsup
          ELSE
            xf2_max=spm2lo
            rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
            rk3=rkappa3high/rkappa4low
ccc            auxalfa3=rk3/(rk3+1.d0)
            auxalfa3=0.5d0
          ENDIF
          xf2_min=spm1hi
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((spm1lo**2-rmx1**2)/(gamx1*rmx1))
          rkappa2low=(xw1_max-xw1_min)*
     &          ((spm1lo**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3low*auxalfa3
ccc          auxalfa2=rk2/(rk2+1.d0)
          auxalfa2=1.d0/3.d0
          xf1_max=spm1lo
          xf1_min=estrinf
          rkappa1=(xf1_max-xf1_min)*2.d0*spm1lo
          rk1=rkappa1/rkappa2low*auxalfa2
ccc          alfa1=rk1/(rk1+1.d0)
          alfa1=1.d0/4.d0
          alfa2=auxalfa2*(1.d0-alfa1)+alfa1
          alfa3=auxalfa3*(1.d0-auxalfa2)*(1.d0-alfa1)+alfa2
          IF(rkappa4low.LE.0.d0) alfa3=1.d0    
        ENDIF
      ELSE ! estrsup.GT.spm2hi: flat1 + res1 + flat2 + res2 + flat3?
        IF(estrinf.GE.spm2hi)THEN  ! only flat3
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          alfa4=0.d0
          auxalfa4=alfa4
          auxalfa3=alfa3
          auxalfa2=alfa2
          xf3_max=estrsup
          xf3_min=estrinf
        ELSEIF(estrinf.GE.spm2lo)THEN  ! res2 + flat3
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          auxalfa3=alfa3
          auxalfa2=alfa2
          xf3_max=estrsup
          xf3_min=spm2hi
          rkappa5=(xf3_max-xf3_min)*2.d0*spm2hi
          xw2_min=atan((estrinf**2-rmx2**2)/(gamx2*rmx2))
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5
ccc          alfa4=rk4/(rk4+1.d0)
          alfa4=0.5d0
          auxalfa4=alfa4
        ELSEIF(estrinf.GE.spm1hi)THEN  ! flat2 + res2 + flat3
          alfa1=0.d0
          alfa2=0.d0
          auxalfa2=alfa2
          xf3_max=estrsup
          xf3_min=spm2hi
          rkappa5=(xf3_max-xf3_min)*2.d0*spm2hi
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5
ccc          auxalfa4=rk4/(rk4+1.d0)
          auxalfa4=0.5d0
          xf2_max=spm2lo
          xf2_min=estrinf
          rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
          rk3=rkappa3high/rkappa4low*auxalfa4
ccc          alfa3=rk3/(rk3+1.d0)
          alfa3=1.d0/3.d0
          auxalfa3=alfa3
          alfa4=auxalfa4*(1.d0-alfa3)+alfa3
        ELSEIF(estrinf.GE.spm1lo)THEN  ! res1 + flat2 + res2 + flat3
          alfa1=0.d0
          xf3_max=estrsup
          xf3_min=spm2hi
          rkappa5=(xf3_max-xf3_min)*2.d0*spm2hi
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5
ccc          auxalfa4=rk4/(rk4+1.d0)
          auxalfa4=0.5d0
          xf2_max=spm2lo
          xf2_min=spm1hi
          rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          rk3=rkappa3high/rkappa4low*auxalfa4
ccc          auxalfa3=rk3/(rk3+1.d0)
          auxalfa3=1.d0/3.d0
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((estrinf**2-rmx1**2)/(gamx1*rmx1))
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3low*auxalfa3
ccc          alfa2=rk2/(rk2+1.d0)
          alfa2=1.d0/4.d0
          auxalfa2=alfa2
          alfa3=auxalfa3*(1.d0-alfa2)+alfa2
          alfa4=auxalfa4*(1.d0-auxalfa3)*(1.d0-alfa2)+alfa3
        ELSE ! estrinf.GT.spm1lo: flat1 + res1 + flat2 + res2+ flat3
          xf3_max=estrsup
          xf3_min=spm2hi
          rkappa5=(xf3_max-xf3_min)*2.d0*spm2hi
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5
ccc          auxalfa4=rk4/(rk4+1.d0)
          auxalfa4=0.5d0
          xf2_max=spm2lo
          xf2_min=spm1hi
          rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          rk3=rkappa3high/rkappa4low*auxalfa4
ccc          auxalfa3=rk3/(rk3+1.d0)
          auxalfa3=1.d0/3.d0
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((spm1lo**2-rmx1**2)/(gamx1*rmx1))
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rkappa2low=(xw1_max-xw1_min)*
     &          ((spm1lo**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3low*auxalfa3
ccc          auxalfa2=rk2/(rk2+1.d0)
          auxalfa2=1.d0/4.d0
          xf1_max=spm1lo
          xf1_min=estrinf
          rkappa1=(xf1_max-xf1_min)*2.d0*spm1lo
          rk1=rkappa1/rkappa2low*auxalfa2
ccc          alfa1=rk1/(rk1+1.d0)
          alfa1=1.d0/5.d0
          alfa2=auxalfa2*(1.d0-alfa1)+alfa1
          alfa3=auxalfa3*(1.d0-auxalfa2)*(1.d0-alfa1)+alfa2
          alfa4=auxalfa4*(1.d0-auxalfa3)*(1.d0-auxalfa2)*(1.d0-alfa1)
     &                               +alfa3
        ENDIF
      ENDIF
      
      RETURN
      END

      
      
*****************************************************************************
* Auxiliary subroutine. Defines extrema for flat+t**alfa mapping of t-channel
* Propagator. Checks for position of t_sup,t_inf with respect to spt1
* INPUT:
*   t_sup,t_inf      : upper and lower limit for t
*   spt1             : separation point between flat and mapped regions
*   exp1             : exponent of the mapping
* OUTPUT:
*  alfa              : separation points in [0,1] between flat and mapped regions
*                      Can also be 0 or 1. Some care necessary in this cases
*  xt1_min,xt1_max,xt2_min,xt2_max
*                    : auxiliary quantities used for generating t
*                      and computing the jacobian
*****************************************************************************
* a --> c, b --> d

      SUBROUTINE T_INI(t_inf,t_sup,spt1,exp1,
     &                 alfa,xt1_min,xt1_max,xt2_min,xt2_max)
      IMPLICIT REAL*8(a-h,o-z)

      IF(t_sup.LE.spt1)THEN      ! only flat
        alfa=1.d0
        xt1_min=t_inf
        xt1_max=t_sup
      ELSE                       ! flat + mapped?
        IF(t_inf.GE.spt1)THEN    ! only mapped
          alfa=0.d0
          IF(exp1.eq.0.d0)THEN
            xt2_min=t_inf
            xt2_max=t_sup   
          ELSEIF(exp1.ne.1.d0)THEN 
            xt2_min=(-t_sup)**(1.d0-exp1)
            xt2_max=(-t_inf)**(1.d0-exp1)
          ELSE
            xt2_min=log((-t_sup))
	    xt2_max=log((-t_inf))
          ENDIF
	ELSE                     ! flat + mapped!
          xt1_min=t_inf
          xt1_max=spt1
          rkappa1=(xt1_max-xt1_min)
          IF(exp1.eq.0.d0)THEN
            xt2_min=spt1
            xt2_max=t_sup   
            rkappa2=(xt2_max-xt2_min)
          ELSEIF(exp1.ne.1.d0)THEN 
            xt2_min=(-t_sup)**(1.d0-exp1)
            xt2_max=(-spt1)**(1.d0-exp1)
	    rkappa2=(xt2_max-xt2_min)/(1.d0-exp1)*(-spt1)**exp1
          ELSE
            xt2_min=log((-t_sup))
	    xt2_max=log((-spt1))
	    rkappa2=(xt2_max-xt2_min)*(-spt1)
          ENDIF
          rkappa=rkappa1/rkappa2
          alfa=rkappa/(rkappa+1.d0)
	ENDIF
      ENDIF 
	
      RETURN
      END
        

*****************************************************************************
* Auxiliary subroutine. Defines extrema for flat+flat mapping
* of shat. Checks for position of estrsup,estrinf with respect to
* spm1lo,spm1hi
* INPUT:
*   estrsup,estrinf: upper and lower limit for the generated mass
*   spm1  : separation point of the two regions region
* OUTPUT:
*  alfa1             : separation point in [0,1] between the two regions
*                      Can also be 0 or 1. Some care necessary in this cases.
*                      At present set to 2/3.
*  xf1_min,xf1_max,xf2_min,xf2_max
*                    : auxiliary quantities used for generating the mass
*                      and computing the jacobian
*****************************************************************************
      SUBROUTINE PHSP_INI2(estrsup,estrinf,spm1,
     &                           alfa1,
     &               xf1_min,xf1_max,xf2_min,xf2_max)
     
      IMPLICIT REAL*8(a-h,o-z)
      
      IF(estrsup.LE.spm1)THEN ! only flat1
        alfa1=1.d0
        xf1_max=estrsup
        xf1_min=estrinf
      ELSE ! estrsup.GT.spm1: flat1 + flat2 ?
        IF(estrinf.GE.spm1)THEN  ! only flat2
          alfa1=0.d0
          xf2_max=estrsup
          xf2_min=estrinf
        ELSE !  flat1 + flat2
          xf2_max=estrsup
          xf2_min=spm1
          rkappa2=(xf2_max-xf2_min)*2.d0*spm1
          xf1_min=estrinf
          xf1_max=spm1
          rkappa1=(xf1_max-xf1_min)*2.d0*spm1
          rk1=rkappa1/rkappa2
c          alfa1=rk1/(rk1+1.d0)
          alfa1=0.85d0
        ENDIF
      ENDIF
      
      RETURN
      END

      
*****************************************************************************
* Auxiliary subroutine
* INPUT:
*       pbst        : intermediate particle momentum
*       shat        : the corresponding mass (assumed positive)
*       pmod        : the absolute value of the 3-momentum of the two
*                     daughter particles in the decay CM
*       xm1sq,xm2sq : the squared masses of the daughter particles
*       x1,x2       : two random numbers uniformly distributed in [0,1]
* OUTPUT:
*       p1,p2       : the momenta of the two daughters in the frame in
*                     which the intermediate particle has momentum pbst 
*                             p1 + p2 = pbst
*
* If pbst=(sqrt(shat),0,0,0) the boost reduces to the identity
* Ezio May 7 2007
* A rotation of the generated momentum has been added in such a way to assume as
* z-axis the direction of a massless particle along the z-axis in the overall
* center of mass
* 
*****************************************************************************
      SUBROUTINE TWOBODY_Z(pbst,shat,pmod,xm1sq,xm2sq,x1,x2,p1,p2)
      IMPLICIT NONE
      REAL*8 pbst(0:3),shat,pmod,xm1sq,xm2sq,x1,x2,p1(0:3),p2(0:3)
* Local variables
      REAL*8 app,costh,sinth,ph,rmboost,aapp,pi,twopi
      REAL*8 pz(0:3)
      INTEGER i
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi)
       
      costh=min(x1*2.d0-1.d0,1.d0)
      costh=max(costh,-1.d0)
      sinth=sqrt(1.d0-costh*costh)
      ph=x2*twopi
* p1
      p1(0)=sqrt(pmod*pmod+xm1sq)
      p1(1)=pmod*sinth*cos(ph)
      p1(2)=pmod*sinth*sin(ph)
      p1(3)=pmod*costh
* pz
      pz(0)=1.d0
      pz(1)=0.d0      
      pz(2)=0.d0      
      pz(3)=1.d0
* boost to pbst center of mass
      call boostinv(pz,pbst,pz)      
* rotate in such a way that angles are relative to the direction of pz (bbosted)
      call rot(p1(1),pz(1),p1(1))
* boost back
      rmboost=sqrt(shat)
      app=(p1(0)*pbst(0)+p1(1)*pbst(1)+p1(2)*pbst(2)+p1(3)*pbst(3))
     &      /rmboost
      aapp=(app+p1(0))/(pbst(0)+rmboost)
      p1(0)=app
      p1(1)=p1(1)+aapp*pbst(1)
      p1(2)=p1(2)+aapp*pbst(2)
      p1(3)=p1(3)+aapp*pbst(3)
      DO i=0,3
        p2(i)=pbst(i)-p1(i)
      ENDDO
      
      RETURN
      END

*****************************************************************************
* Auxiliary subroutine. Determines the invariant mass squared of a branching.
* The mass is supposed to resonate at two different masses rmx1 and rmx2.
* INPUT:
*   xvar: a random number in [0,1]
*   rmx1,gamx1     : mass and width of first resonance
*   rmx2,gamx2     : mass and width of second resonance
*   iregion        : index for filling the iphsp_fxn COMMON
* 
*  alfa1,alfa2,alfa3,alfa4: separation points in [0,1] between the five regions
*                      Can also be 0 or 1. Some care necessary in this cases
*  auxalfa2,xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,xw2_min,xw2_max,
*  xf3_min,xf3_max,auxalfa3,auxalfa4
*                    : auxiliary quantities used for generating the mass
*                      and computing the jacobian
* OUTPUT:
*   svar: the generated invariant mass squared
*   xjacvar: the jacobian
*****************************************************************************
      SUBROUTINE PHSP_5(xvar,
     &                     rmx1,gamx1,rmx2,gamx2,
     &                     alfa1,alfa2,auxalfa2,
     &                     alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &                     iregion,
     &                     svar,xjacvar)
     
      IMPLICIT REAL*8(a-h,o-z)
c giuseppe 29/09/2007
      INTEGER maxnresonant
      PARAMETER (maxnresonant=4)
      INTEGER isresonant(maxnresonant),imapregion(maxnresonant)
      COMMON /iphsp_fxn/ isresonant,imapregion      
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi)

* Compute the invariant mass squared and the corresponding jacobian

      xjacvar=1.d0

      IF(xvar.LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=xvar/alfa1
        evar=(xf1_max-xf1_min)*vegas_x+xf1_min
        svar=evar*evar
        xjacvar=xjacvar*(xf1_max-xf1_min)*2.d0*evar/alfa1
      ELSEIF(xvar.LE.alfa2.AND.alfa2.GT.0.d0)THEN
        vegas_x=(xvar-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        svar=(gamx1*rmx1*tan((x_max-x_min)*vegas_x+x_min)+rmx1**2)
        IF(svar.LE.0.d0)THEN
          xjacvar=0.d0
          RETURN
        ENDIF
        xjacvar=xjacvar*(x_max-x_min)*((svar-rmx1**2)**2+
     &      (gamx1*rmx1)**2)/(gamx1*rmx1)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(iregion)=1
        imapregion(iregion)=2
c end giuseppe 29/09/2007
      ELSEIF(xvar.LE.alfa3.AND.alfa3.GT.0.d0)THEN
        vegas_x=(xvar-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        evar=(xf2_max-xf2_min)*vegas_x+xf2_min
        svar=evar*evar
        xjacvar=xjacvar*(xf2_max-xf2_min)*2.d0*evar
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(xvar.LE.alfa4.AND.alfa4.GT.0.d0)THEN
        vegas_x=(xvar-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/auxalfa4
        x_max=xw2_max
        x_min=xw2_min
        svar=(gamx2*rmx2*tan((x_max-x_min)*vegas_x+x_min)+rmx2**2)
        IF(svar.LE.0.d0)THEN
          xjacvar=0.d0
          RETURN
        ENDIF
        xjacvar=xjacvar*(x_max-x_min)*((svar-rmx2**2)**2+
     &      (gamx2*rmx2)**2)/(gamx2*rmx2)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
c giuseppe 29/09/2007
        isresonant(iregion)=1
        imapregion(iregion)=4
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(xvar-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
        evar=(xf3_max-xf3_min)*vegas_x+xf3_min
        svar=evar*evar
        xjacvar=xjacvar*(xf3_max-xf3_min)*2.d0*evar
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF      
      xjacvar=xjacvar/twopi
      
      RETURN
      END



      
