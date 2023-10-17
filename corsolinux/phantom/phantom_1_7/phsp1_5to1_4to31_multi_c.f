c
c Last update: Sep 29, 2007
c
************************************************************************
* Subroutine to compute the momenta p1,....,p6 of six particles of mass 
* xm1,...,xm6 given a point in the 14 dimensional hypercube [0,1]^14
* Returns also xjac, the weight due to the jacobian.
*
* The phase space mapping corresponds to 
*      [1t,a-->p6]+[5(rm1) --> 4(rm11)+p5 --> 3(rm111)+p4+p5 
*                                         --> 2(rm1111,p1,p2)+p3+p4+p5] 
*      ("t-channel process")
*   and
*      [5(rm1) --> 4(rm11)+p5 --> 3(rm111)+p4+p5
*                             --> 2(rm1111,p1,p2)+p3+p4+p5]
*      ("s-channel process")
*
* where rm1,rm11,rm111, rm1111 are the masses of the intermediate 
*  resonances
*
* rm1,rm11,rm111 can resonate at top mass, rm1111 can resonate at W mass
*
* One exponent (expa) is available for mapping the t-channel propagator.
* In case of an "s-channel process", expa is set equal to 0 in order the
* mapping on the t-channel to be flat.
* rmta is the corresponding t-channel mass.
* This mapping is applied for 1.d0>cos(theta)>cos_sep. cos_sep is a 
*   PARAMETER
* Outside this range the cosine is generated flat. 
*
* If in the main x had dimension 15 with x(1) and x(2) fixing the 
*  initial state particle momentum fraction the subrotine can be called 
*  as call phsp1_5to1_4to31(x(3),shat,.......)
*
* INPUT:
*    x,shat        : total energy squared
*    xm(8),xmsq(8) : masses and masses squared of the particles
*                        in the order used in the main program      
*    rmN,gamN      : N=1,11,111,1111 mass and width for intermediate 
*                        resonance N.
*    spmNlo,spmNhi : N=1,11,111,1111 low and high limit of resonance 
*                       mapping. outside this range the mapping is flat
*    pt6min_in(3)  : array of possible minimum pt's of particle 6,
*                    depending on its flavour. The first element refers 
*                    to jets, the second one to charged leptons and the 
*                    third one to neutrinos.
*    iorder(6)     : specifies the permutation that relates internal and
*                    external momenta: p = p(i,iorder(1))
*    id1,id2       : idp(iorder(1)),idp(iorder(2)) used to determine the
*                    direction of the incoming particles involved in 
*                    the t-channel 
*    id6           : idp(iorder(8)) used to determine which external 
*                    particle 
*
* OUTPUT:
*    p(0:3,8)      : momenta in the hard scattering CM.
*                        Includes initial momenta.
*    xjac          : full jacobian
*
* INTERNAL VARIABLES
* xmasq,xmbsq are the masses squared of the two incoming particles 
* (a --> 6  in the t)
* The routine determines from iorder whether
* pa=sqrt(shat)/2*(1,0,0,1), pb=sqrt(shat)/2*(1,0,0,-1) or viceversa
* in the massless limit.
************************************************************************

      subroutine phsp1_5to1_4to31_multi_c(x,shat,
     &                         rm1,gam1,spm1lo,spm1hi,
     &                         rm11,gam11,spm11lo,spm11hi,
     &                         rm111,gam111,spm111lo,spm111hi,
     &                         rm1111,gam1111,spm1111lo,spm1111hi,
     &                         rmta,expa_in,spta,
     &                         pt6min_in,
     &                         id1,id2,id6,
     &                         xmasq_in,xmbsq_in,
     &                         xm,xmsq,iorder,
     &                         p,xjac)

      IMPLICIT NONE
      REAL*8 x(13),shat,
     &     rm1,gam1,spm1lo,spm1hi,
     &     rm11,gam11,spm11lo,spm11hi,
     &     rm111,gam111,spm111lo,spm111hi,
     &     rm1111,gam1111,spm1111lo,spm1111hi,
     &     rmta,expa,expa_in,
     &     pt6min_in(3),
     &     pt6min,
     &     xmasq_in,xmbsq_in,
     &     xm(8),xmsq(8),
     &     p(0:3,8),xjac
      INTEGER iorder(6),id1,id2,id6
* Local variables
      REAL*8 e_sup,e_inf,x_max,x_min,app,
     &       s1111,s111,s11,s1,
     &       e1111,e111,e11,e1,
     &       q1111(0:3),q111(0:3),q11(0:3),q1(0:3),
     &       p1111mod,p111mod,p11mod,p1mod,pmod,p0mod,
     &       t_inf,t_sup,t1,
     &       pi,twopi,fourpi,
     &       costh,sinth,ph,ran2,
     &       xmasq,xmbsq,rmtasq,
     &       signfactor,smallmass

      INTEGER idum,io1,io2,io3,io4,io5,io6,iaux,ischannel,i

      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi,smallmass=1.d-3)

* Multimapping
      REAL*8 alfa1,alfa2,auxalfa2,
     &       xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &       vegas_x
* Multimapping t
      REAL*8 alfa,xt1_min,xt1_max,xt2_min,xt2_max,spta
      REAL*8 cos_sep
      PARAMETER(cos_sep=0.4d0)

      COMMON /pharand/ idum
c giuseppe 29/09/2007
      INTEGER maxnresonant
      PARAMETER (maxnresonant=4)
      INTEGER isresonant(maxnresonant),imapregion(maxnresonant)
      COMMON /iphsp_fxn/ isresonant,imapregion

      DO i=1,maxnresonant
         isresonant(i)=0
         imapregion(i)=0
      ENDDO
c end giuseppe 29/09/2007

      io1=iorder(1)
      io2=iorder(2)
      io3=iorder(3)
      io4=iorder(4)
      io5=iorder(5)
      io6=iorder(6)


c Determine which minimum pt must be assigned to t-channel, depending
c on final-state fermion's flavour
      IF(abs(id6).LE.5 .OR. id6.EQ.21)THEN
        pt6min=pt6min_in(1)
      ELSEIF(abs(id6).EQ.11)THEN
        pt6min=pt6min_in(2)
      ELSEIF(abs(id6).EQ.12)THEN
        pt6min=pt6min_in(3)
      ELSE
        print*,'***phsp1_5_to1_4to31 ERROR'
        STOP
      ENDIF



      ischannel=0

c In case of an s-channel process, the exponential mapping has to be
c suppressed. The following statement allows the input parameter
c (expa_in) NOT to be modified (thus avoiding possible undesired 
c effects).
      expa=expa_in

* Determine which particle is paired with 6 in t-channel
      IF(id6.EQ.21)THEN   ! 6 is a  gluon. It couples to another gluon
        iaux=id6
        IF(id1.EQ.iaux)THEN
          xmasq=xmasq_in
          xmbsq=xmbsq_in
          signfactor=1.d0
        ELSEIF(id2.EQ.iaux)THEN
          xmasq=xmbsq_in
          xmbsq=xmasq_in
          signfactor=-1.d0
        ELSE
          PRINT*,'WRONG ASSIGNEMENT IN PHASE SPACE 1_5to1_4to31'
          STOP
        ENDIF
      ELSE                      ! 6 is a quark. Use rmta to find partner
        IF(rmta.LT.smallmass)THEN ! t channel boson is neutral
          iaux=-id6
          IF(id1.EQ.iaux)THEN
            xmasq=xmasq_in
            xmbsq=xmbsq_in
            signfactor=1.d0
          ELSEIF(id2.EQ.iaux)THEN
            xmasq=xmbsq_in
            xmbsq=xmasq_in
            signfactor=-1.d0
c Notice: when the phase space 1_5to1_4to31 is used to cover an 
c s-channel process, the value of the parameter rmta in procini.f is 
c assumed to be zero. This is the reason for the following check holds
c for rmta.LT.smallmass only.
          ELSEIF( (id1+id2).EQ.0 .OR. (id1+id2).EQ.1 .OR. 
     &           (id1+id2).EQ.-1)THEN
            ischannel=1         ! s-channel process
          ELSE
            PRINT*,'WRONG ASSIGNEMENT IN PHASE SPACE 1_5to1_4to31'
            STOP
          ENDIF
        ELSE                    ! t channel boson is charged
	  IF(mod(abs(id6),2).EQ.0)THEN ! id6 even: 6 is up-type
	    IF(id6.GT.0)THEN    ! quark
              iaux=-(id6-1)
	    ELSE                ! antiquark
	      iaux=-(id6+1)
	    ENDIF 
          ELSE                  ! id6 odd: 6 is down-type
	    IF(id6.GT.0)THEN    ! quark
              iaux=-(id6+1)
	    ELSE                ! antiquark
	      iaux=-(id6-1)
	    ENDIF 
          ENDIF
          IF(id1.EQ.iaux)THEN
            xmasq=xmasq_in
            xmbsq=xmbsq_in
            signfactor=1.d0
          ELSEIF(id2.EQ.iaux)THEN
            xmasq=xmbsq_in
            xmbsq=xmasq_in
            signfactor=-1.d0
          ELSE
            PRINT*,'WRONG ASSIGNEMENT IN PHASE SPACE 1_5to1_4to31'
            STOP
	  ENDIF  
        ENDIF                   !IF(rmta.LT.smallmass)
      ENDIF                     ! IF(id6.EQ.21)

      IF(ischannel.EQ.1)THEN
        xmasq=xmasq_in
        xmbsq=xmbsq_in
        signfactor=1.d0
* In case of s-channel process, set expa=0 in order to avoid mapping on
* the t-channel variable
        expa=0.d0
      ENDIF                     ! IF(ischannel.EQ.1)

      xjac=1.d0

* s1 = (q11+p5)**2 = q1**2
      e_sup = sqrt(shat)-xm(io6)
      e_inf = xm(io1)+xm(io2)+xm(io3)+xm(io4)+xm(io5)
      IF(e_sup.LE.e_inf)THEN    
        xjac=0.d0
        RETURN
      ENDIF
      CALL PHSP_INI(e_sup,e_inf,
     &              spm1lo,spm1hi,rm1,gam1,
     &              alfa1,alfa2,auxalfa2,
     &              xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max)
      IF(x(1).LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=x(1)/alfa1
        e1=(xf1_max-xf1_min)*vegas_x+xf1_min
        s1=e1*e1
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e1/alfa1
      ELSEIF(x(1).LE.alfa2.AND.alfa2.GT.0.d0)THEN
        vegas_x=(x(1)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        s1=(gam1*rm1*tan((x_max-x_min)*vegas_x+x_min)+rm1**2)
        IF(s1.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s1-rm1**2)**2+
     &       (gam1*rm1)**2)/(gam1*rm1)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(1)=1
        imapregion(1)=2
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(1)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)
        e1=(xf2_max-xf2_min)*vegas_x+xf2_min
        s1=e1*e1
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e1
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
      ENDIF
      xjac=xjac/twopi


* pa+pb --> p6 + q1; (pa-p6)**2 = (pb-q1)**2 = t1
      rmtasq=rmta*rmta
      app=((shat-s1-xmsq(io6))**2-4.d0*s1*xmsq(io6))/(4.d0*shat)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      pmod=sqrt(app)
      app=((shat-xmasq-xmbsq)**2-4.d0*xmasq*xmbsq)/(4.d0*shat)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p0mod=sqrt(app)

      IF(pt6min.LE.0.d0)THEN
        t_inf = xmbsq+s1-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                (shat+s1-xmsq(io6))
     &                -2.d0*p0mod*pmod
        app=(xmbsq-s1)*(xmasq-xmsq(io6))
     &                +(s1-xmbsq+xmasq-xmsq(io6))*
     &                (s1*xmasq-xmbsq*xmsq(io6))/shat
        t_sup = app/t_inf
      ELSE
        IF(pt6min.GE.pmod)THEN
          xjac=0.d0
          RETURN
        ENDIF
        app=sqrt((pmod-pt6min)*(pmod+pt6min))
        t_inf = xmbsq+s1-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                (shat+s1-xmsq(io6))
     &               -2.d0*p0mod*app
        t_sup = xmbsq+s1-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                (shat+s1-xmsq(io6))
     &               +2.d0*p0mod*app
      ENDIF

      spta=xmbsq+s1-0.5d0/shat*(shat+xmbsq-xmasq)*
     &     (shat+s1-xmsq(io6))
     &     +2.d0*p0mod*pmod*cos_sep

* From now on we use t1 = (pa-p6)**2 - rmtasq = (pb-q1)**2 - rmtasq
      t_inf=t_inf-rmtasq
      t_sup=t_sup-rmtasq
      spta=spta-rmtasq
      CALL T_INI(t_inf,t_sup,spta,expa,
     &           alfa,xt1_min,xt1_max,xt2_min,xt2_max)
      IF(x(2).LE.alfa .AND. alfa.GT.0.d0)THEN
        vegas_x=x(2)/alfa
        t1=(xt1_max-xt1_min)*vegas_x+xt1_min
        xjac=xjac*0.5d0*(xt1_max-xt1_min)/p0mod/pmod
        xjac=xjac/alfa
      ELSE
        vegas_x=(x(2)-alfa)/(1.d0-alfa)
        IF(expa.eq.0.d0)THEN 
          t1=(xt2_max-xt2_min)*vegas_x+xt2_min
          xjac=xjac*0.5d0*(xt2_max-xt2_min)/p0mod/pmod
        ELSEIF(expa.ne.1.d0)THEN
          t1=-((xt2_max-xt2_min)*vegas_x+xt2_min)**(1.d0/(1.d0-expa))
          xjac=xjac*0.5d0*(xt2_max-xt2_min)
     &         /(1.d0-expa)/p0mod/pmod*(-t1)**expa 
        ELSE
          t1=-exp((xt2_max-xt2_min)*vegas_x+xt2_min)
          xjac=xjac*0.5d0*(xt2_max-xt2_min)/p0mod/pmod*(-t1)
        ENDIF
        xjac=xjac/(1.d0-alfa)
      ENDIF

      costh=0.5*(t1+rmtasq-xmbsq-s1
     &     +0.5d0/shat*(shat+xmbsq-xmasq)*(shat+s1-xmsq(io6)))
     &     /p0mod/pmod
      costh=costh*signfactor

* q1 and p6
      costh=min(costh,1.d0)
      costh=max(costh,-1.d0)
      sinth=sqrt(1.d0-costh*costh)
      ph=ran2(idum)*twopi
      p(0,io6)=sqrt(pmod*pmod+xmsq(io6))
      p(1,io6)=pmod*sinth*cos(ph)
      p(2,io6)=pmod*sinth*sin(ph)
      p(3,io6)=pmod*costh
      q1(0)=sqrt(pmod*pmod+s1)
      q1(1)=-p(1,io6)
      q1(2)=-p(2,io6)
      q1(3)=-p(3,io6)
* beta/2/(4*pi)**2
      xjac=xjac*pmod/sqrt(shat)/(fourpi)**2
      xjac=xjac*twopi

* s11 = (q111+p4)**2 = q11**2
      e_sup = sqrt(s1)-xm(io5)
      e_inf = xm(io1)+xm(io2)+xm(io3)+xm(io4)
      IF(e_sup.LE.e_inf)THEN    
        xjac=0.d0
        RETURN
      ENDIF
      CALL PHSP_INI(e_sup,e_inf,
     &              spm11lo,spm11hi,rm11,gam11,
     &              alfa1,alfa2,auxalfa2,
     &              xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max)
      IF(x(3).LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=x(3)/alfa1
        e11=(xf1_max-xf1_min)*vegas_x+xf1_min
        s11=e11*e11
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e11/alfa1
      ELSEIF(x(3).LE.alfa2.AND.alfa2.GT.0.d0)THEN
        vegas_x=(x(3)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        s11=(gam11*rm11*tan((x_max-x_min)*vegas_x+x_min)+rm11**2)
        IF(s11.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s11-rm11**2)**2+
     &       (gam11*rm11)**2)/(gam11*rm11)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(2)=1
        imapregion(2)=2
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(3)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)
        e11=(xf2_max-xf2_min)*vegas_x+xf2_min
        s11=e11*e11
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e11
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
      ENDIF
      xjac=xjac/twopi

* s111 = (q1111+p5)**2 = q111**2
      e_sup = sqrt(s11)-xm(io4)
      e_inf = xm(io1)+xm(io2)+xm(io3)
      IF(e_sup.LE.e_inf)THEN    
        xjac=0.d0
        RETURN
      ENDIF
      CALL PHSP_INI(e_sup,e_inf,
     &              spm111lo,spm111hi,rm111,gam111,
     &              alfa1,alfa2,auxalfa2,
     &              xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max)
      IF(x(4).LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=x(4)/alfa1
        e111=(xf1_max-xf1_min)*vegas_x+xf1_min
        s111=e111*e111
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e111/alfa1
      ELSEIF(x(4).LE.alfa2.AND.alfa2.GT.0.d0)THEN
        vegas_x=(x(4)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        s111=(gam111*rm111*tan((x_max-x_min)*vegas_x+x_min)+rm111**2)
        IF(s111.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s111-rm111**2)**2+
     &       (gam111*rm111)**2)/(gam111*rm111)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(3)=1
        imapregion(3)=2
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(4)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)
        e111=(xf2_max-xf2_min)*vegas_x+xf2_min
        s111=e111*e111
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e111
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
      ENDIF
      xjac=xjac/twopi

* s1111 = (p1+p2)**2 = q1111**2
      e_sup = sqrt(s111)-xm(io3)
      e_inf = xm(io1)+xm(io2)
      IF(e_sup.LE.e_inf)THEN
        xjac=0.d0
        RETURN
      ENDIF
      CALL PHSP_INI(e_sup,e_inf,
     &              spm1111lo,spm1111hi,rm1111,gam1111,
     &              alfa1,alfa2,auxalfa2,
     &              xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max)
      IF(x(5).LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=x(5)/alfa1
        e1111=(xf1_max-xf1_min)*vegas_x+xf1_min
        s1111=e1111*e1111
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e1111/alfa1
      ELSEIF(x(5).LE.alfa2.AND.alfa2.GT.0.d0)THEN
        vegas_x=(x(5)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        s1111=(gam1111*rm1111*tan((x_max-x_min)*vegas_x+x_min)
     &       +rm1111**2)
        IF(s1111.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s1111-rm1111**2)**2+
     &      (gam1111*rm1111)**2)/(gam1111*rm1111)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(4)=1
        imapregion(4)=2
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(5)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)
        e1111=(xf2_max-xf2_min)*vegas_x+xf2_min
        s1111=e1111*e1111
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e1111
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
      ENDIF
      xjac=xjac/twopi


* q11 and p5
      app=((s1-s11-xmsq(io5))**2-4.d0*s11*xmsq(io5))/(4.d0*s1)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p1mod=sqrt(app)
      CALL TWOBODY(q1,s1,p1mod,s11,xmsq(io5),x(6),x(7),q11,p(0,io5))
* beta/2/(4*pi)**2
      xjac=xjac*p1mod/sqrt(s1)/(fourpi)**2
      xjac=xjac*fourpi

* q111 and p4
      app=((s11-s111-xmsq(io4))**2-4.d0*s111*xmsq(io4))/(4.d0*s11)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p11mod=sqrt(app)
      CALL TWOBODY(q11,s11,p11mod,s111,xmsq(io4),x(8),x(9),q111,
     &     p(0,io4))
* beta/2/(4*pi)**2
      xjac=xjac*p11mod/sqrt(s11)/(fourpi)**2
      xjac=xjac*fourpi

* q1111 and p3
      app=((s111-s1111-xmsq(io3))**2-4.d0*s1111*xmsq(io3))/(4.d0*s111)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p111mod=sqrt(app)
      CALL TWOBODY(q111,s111,p111mod,s1111,xmsq(io3),x(10),x(11),q1111,
     &     p(0,io3))
* beta/2/(4*pi)**2
      xjac=xjac*p111mod/sqrt(s111)/(fourpi)**2
      xjac=xjac*fourpi

* q1111 --> p1 and p2
      app=((s1111-xmsq(io1)-xmsq(io2))**2-4.d0*xmsq(io1)*xmsq(io2))
     &     /(4.d0*s1111)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p1111mod=sqrt(app)
      CALL TWOBODY(q1111,s1111,p1111mod,xmsq(io1),xmsq(io2),x(12),
     &     x(13),p(0,io1),p(0,io2))
* beta/2/(4*pi)**2
      xjac=xjac*p1111mod/sqrt(s1111)/(fourpi)**2
      xjac=xjac*fourpi

      RETURN 
      END

