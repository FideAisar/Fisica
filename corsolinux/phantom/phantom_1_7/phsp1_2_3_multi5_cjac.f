c
c Last update: Jul 07, 2007
c
************************************************************************
* Subroutine to compute the jacobian corresponding to the phase space 
*   mapping:
*      [1t,a-->p6]+ 5(rm0[a,b])
*                    -->  [ [3(rm1[a,b]) --> 2(rm11[a,b],p1,p2)+p3] +
*                          + 2(rm2[a,b],p4,p5) ]
* where rm0[a,b],rm1[a,b],rm11[a,b],rm2[a,b] are the masses of the 
*  intermediate resonances
* rm0,rm1,rm11 and rm2 can resonate at two different masses 
*  (e.g. Z and H) with rm0a < rm0b, rm1a < rm1b, rm11a < rm11b, 
*  rm2a < rm2b.
* If only one resonance is needed set the corresponding 
*  spmXahi=spmXalo=0.d0
*
* One exponent (expa) is available for mapping the t-channel propagator.
* rmta is the corresponding t-channel mass.
* This mapping is applied for 1.d0>cos(theta)>cos_sep. cos_sep is a 
*   PARAMETER
* Outside this range the cosine is generated flat. 
*
* INPUT:
*    x,shat        : vegas random points, total energy squared
*    xm(8),xmsq(8) : masses and masses squared of the particles
*                        in the order used in the main program      
*    rmN,gamN      : N=0[a,b],1[a,b],11[a,b],2[a,b] mass and width 
*                     for intermediate resonance N.
*    spmNlo,spmNhi : N=0[a,b],1[a,b],11[a,b],2[a,b] low and high 
*                     limit of resonance mapping. Outside this range 
*                     the mapping is flat
*    pt6min_in(3)  : array of possible minimum pt's of particle 6,
*                    depending on its flavour. The first element refers 
*                    to jets, the second one to charged leptons and the 
*                    third one to neutrinos.
*    iorder(6)     : specifies the permutation that relates internal and
*                        external momenta: p = p(i,iorder(1))
*    id1,id2       : idp(iorder(1)),idp(iorder(2)) used to determine the
*                    direction of the incoming particles involved in the
*                    t-channel 
*    id6           : idp(iorder(8)) used to determine which external 
*                    particle
*                    particles in the final state.
* xmasq_in,xmbsq_in : masses squared of the two incoming particles
*    p(0:3,8)      : momenta in the hard scattering CM.
*                        Includes initial momenta.
* OUTPUT:
*    xjac          : full jacobian
*
* INTERNAL VARIABLES
* xmasq,xmbsq are the masses squared of the two incoming particles 
* (a --> 6  in the t)
* The routine determines from iorder whether
* pa=sqrt(shat)/2*(1,0,0,1), pb=sqrt(shat)/2*(1,0,0,-1) or viceversa
* in the massless limit.
************************************************************************

        subroutine phsp1_2_3_multi5_cjac(shat,
     &                    rm0a,gam0a,spm0alo,spm0ahi,
     &                    rm0b,gam0b,spm0blo,spm0bhi,
     &                    rm1a,gam1a,spm1alo,spm1ahi,
     &                    rm1b,gam1b,spm1blo,spm1bhi,
     &                    rm11a,gam11a,spm11alo,spm11ahi,
     &                    rm11b,gam11b,spm11blo,spm11bhi,
     &                    rm2a,gam2a,spm2alo,spm2ahi,
     &                    rm2b,gam2b,spm2blo,spm2bhi,
     &                    rmta,expa,spta,
c giuseppe 07/07/2007
c     &                    pt6min,
     &                    pt6min_in,
c end giuseppe 07/07/2007
     &                    id1,id2,id6,
     &                    xmasq_in,xmbsq_in,
     &                    xm,xmsq,iorder,
     &                    p,xjac)
     
      IMPLICIT NONE
      REAL*8 shat,
     &       rm0a,gam0a,spm0alo,spm0ahi,
     &       rm0b,gam0b,spm0blo,spm0bhi,
     &       rm1a,gam1a,spm1alo,spm1ahi,
     &       rm1b,gam1b,spm1blo,spm1bhi,
     &       rm11a,gam11a,spm11alo,spm11ahi,
     &       rm11b,gam11b,spm11blo,spm11bhi,
     &       rm2a,gam2a,spm2alo,spm2ahi,
     &       rm2b,gam2b,spm2blo,spm2bhi,
     &       rmta,expa,
c giuseppe 07/07/2007
     &       pt6min_in(3),
c end giuseppe 07/07/2007
     &       pt6min,
     &       xmasq_in,xmbsq_in,
     &       xm(8),xmsq(8),
     &       p(0:3,8),xjac
      INTEGER iorder(6),id1,id2,id6
* Local variables
      REAL*8 e_sup,e_inf,x_max,x_min,app,
     &       s11,s2,s1,
     &       e11,e2,e1,p11mod,p2mod,p1mod,pmod,
     &       saux,eaux,t_inf,t_sup,t1,pauxmod,p0mod,
     &       pi,twopi,fourpi,
     &       xmasq,xmbsq,rmtasq,
     &       signfactor,smallmass
     
      INTEGER io1,io2,io3,io4,io5,io6,i,iaux
     
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi,smallmass=1.d-3)
      
* Mutimapping
      REAL*8 alfa1,alfa2,auxalfa2,
     &       alfa3,auxalfa3,alfa4,auxalfa4,
     &       xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &       xw2_min,xw2_max,xf3_min,xf3_max
* Mutimapping t
      REAL*8 alfa,xt1_min,xt1_max,xt2_min,xt2_max,spta
      REAL*8 cos_sep
      PARAMETER(cos_sep=0.4d0)
      

      io1=iorder(1)
      io2=iorder(2)
      io3=iorder(3)
      io4=iorder(4)
      io5=iorder(5)
      io6=iorder(6)


c giuseppe 07/07/2007
c Determine which minimum pt must be assigned to t-channel, depending
c on final-state fermion's flavour
      IF(abs(id6).LE.5 .OR. id6.EQ.21)THEN
        pt6min=pt6min_in(1)
      ELSEIF(abs(id6).EQ.11)THEN
        pt6min=pt6min_in(2)
      ELSEIF(abs(id6).EQ.12)THEN
        pt6min=pt6min_in(3)
      ELSE
        print*,'***phsp1_2_3_jac ERROR'
        STOP
      ENDIF
c end giuseppe 07/07/2007


* Determine which particle is paired with 6 in t-channel
      IF(id6.EQ.21)THEN    ! 6 is a  gluon. It couples to another gluon
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
          PRINT*,'WRONG ASSIGNEMENT IN PHASE SPACE 1_2_3'
          STOP
        ENDIF
      ELSE                     ! 6 is a quark. Use rmta to find partner                             
        IF(rmta.LT.smallmass)THEN  ! t channel boson is neutral
          iaux=-id6
          IF(id1.EQ.iaux)THEN
            xmasq=xmasq_in
            xmbsq=xmbsq_in
            signfactor=1.d0
          ELSEIF(id2.EQ.iaux)THEN
            xmasq=xmbsq_in
            xmbsq=xmasq_in
            signfactor=-1.d0
          ELSE
            PRINT*,'WRONG ASSIGNEMENT IN PHASE SPACE 1_2_3'
            STOP
          ENDIF
        ELSE                   ! t channel boson is charged
	  IF(mod(abs(id6),2).EQ.0)THEN  ! id6 even: 6 is up-type
	    IF(id6.GT.0)THEN                ! quark
              iaux=-(id6-1)
	    ELSE                        ! antiquark
	      iaux=-(id6+1)
	    ENDIF 
          ELSE                          ! id6 odd: 6 is down-type
	    IF(id6.GT.0)THEN                ! quark
              iaux=-(id6+1)
	    ELSE                        ! antiquark
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
            PRINT*,'WRONG ASSIGNEMENT IN PHASE SPACE 1_2_3'
            STOP
	  ENDIF  
        ENDIF       !IF(rmta.LT.smallmass)
      ENDIF      ! IF(abs(id6).EQ.5)
        
      xjac=1.d0

C All these invariant masses should be computed once and for all and
C stored in some array.

      s11=(p(0,io1)+p(0,io2))**2
      s1=(p(0,io1)+p(0,io2)+p(0,io3))**2
      s2=(p(0,io4)+p(0,io5))**2
      saux=(p(0,io1)+p(0,io2)+p(0,io3)+p(0,io4)+p(0,io5))**2
      DO i=1,3
        s11=s11-(p(i,io1)+p(i,io2))**2
        s1=s1-(p(i,io1)+p(i,io2)+p(i,io3))**2
        s2=s2-(p(i,io4)+p(i,io5))**2
        saux=saux-(p(i,io1)+p(i,io2)+p(i,io3)+p(i,io4)+p(i,io5))**2
      ENDDO 
      IF( s11.LE.1.d-8 .OR. s1.LE.1.d-8 .OR. s2.LE.1.d-8
     &          .OR. saux.LE.1.d-8 )THEN
        xjac=0.d0
        RETURN
      ENDIF

* saux = (q1+q2)**2 = qaux**2
      e_sup = sqrt(shat)-xm(io6)
      e_inf = xm(io1)+xm(io2)+xm(io3)+xm(io4)+xm(io5)
      IF(e_sup.LE.e_inf)THEN    
        xjac=0.d0
      RETURN
      ENDIF
c 5-region multimapping for saux
c      eaux=sqrt(saux)
c      xjac=xjac*(e_sup-e_inf)*2.d0*eaux
      CALL PHSP_INI5(e_sup,e_inf,
     &               spm0alo,spm0ahi,rm0a,gam0a,
     &               spm0blo,spm0bhi,rm0b,gam0b,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max)
      eaux=sqrt(saux)
      IF(eaux.LE.spm0alo .AND. alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*eaux/alfa1
      ELSEIF(eaux.LE.spm0ahi .AND. alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((saux-rm0a**2)**2+
     &       (gam0a*rm0a)**2)/(gam0a*rm0a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(eaux.LE.spm0blo .AND. alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*eaux
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(eaux.LE.spm0bhi .AND. alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((saux-rm0b**2)**2+
     &       (gam0b*rm0b)**2)/(gam0b*rm0b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
        xjac=xjac*(xf3_max-xf3_min)*2.d0*eaux
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF
      xjac=xjac/twopi
      
* pa+pb --> p6 + qaux; (pa-p6)**2 = (pb-qaux)**2 = t1
      rmtasq=rmta*rmta
      app=((shat-saux-xmsq(io6))**2-4.d0*saux*xmsq(io6))/(4.d0*shat)
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
        t_inf =xmbsq+saux-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                                    (shat+saux-xmsq(io6))
     &               -2.d0*p0mod*pmod
        app=(xmbsq-saux)*(xmasq-xmsq(io6))
     &           +(saux-xmbsq+xmasq-xmsq(io6))*
     &                        (saux*xmasq-xmbsq*xmsq(io6))/shat
        t_sup=app/t_inf
      ELSE
        IF(pt6min.GE.pmod)THEN
          xjac=0.d0
          RETURN
        ENDIF
        app=sqrt((pmod-pt6min)*(pmod+pt6min))
        t_inf =xmbsq+saux-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                                    (shat+saux-xmsq(io6))
     &               -2.d0*p0mod*app
        t_sup =xmbsq+saux-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                                    (shat+saux-xmsq(io6))
     &               +2.d0*p0mod*app
      ENDIF
      
      spta=xmbsq+saux-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                                    (shat+saux-xmsq(io6))
     &               +2.d0*p0mod*pmod*cos_sep

* From now on we use t1 = (pa-p6)**2 - rmtasq = (pb-qaux)**2 - rmtasq
      t_inf=t_inf-rmtasq
      t_sup=t_sup-rmtasq
      spta=spta-rmtasq
      CALL T_INI(t_inf,t_sup,spta,expa,
     &                 alfa,xt1_min,xt1_max,xt2_min,xt2_max)
      t1=xmasq+xmsq(io6)-0.5d0/shat*(shat+xmasq-xmbsq)*
     &                                    (shat+xmsq(io6)-saux)
     &               +2.d0*p0mod*p(3,io6)*signfactor-rmtasq
      IF(t1.LE.spta)THEN
        xjac=xjac*0.5d0*(xt1_max-xt1_min)/p0mod/pmod
        xjac=xjac/alfa
      ELSE
        IF(expa.eq.0.d0)THEN 
          xjac=xjac*0.5d0*(xt2_max-xt2_min)/p0mod/pmod
        ELSEIF(expa.ne.1.d0)THEN 
          xjac=xjac*0.5d0*(xt2_max-xt2_min)
     &       /(1.d0-expa)/p0mod/pmod*(-t1)**expa
        ELSE
          xjac=xjac*0.5d0*(xt2_max-xt2_min)/p0mod/pmod*(-t1)
        ENDIF
        xjac=xjac/(1.d0-alfa)
      ENDIF

* qaux and p6
* beta/2/(4*pi)**2
      xjac=xjac*pmod/sqrt(shat)/(fourpi)**2
      xjac=xjac*twopi

* s1 = (q11+p3)**2 = q1**2
      e_sup = sqrt(saux)-xm(io4)-xm(io5)
      e_inf = xm(io1)+xm(io2)+xm(io3)    
      IF(e_sup.LE.e_inf)THEN    
        xjac=0.d0
        RETURN
      ENDIF
      CALL PHSP_INI5(e_sup,e_inf,
     &               spm1alo,spm1ahi,rm1a,gam1a,
     &               spm1blo,spm1bhi,rm1b,gam1b,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max)
      e1=sqrt(s1)
      IF(e1.LE.spm1alo .AND. alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e1/alfa1
      ELSEIF(e1.LE.spm1ahi.AND.alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s1-rm1a**2)**2+
     &       (gam1a*rm1a)**2)/(gam1a*rm1a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e1.LE.spm1blo.AND.alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e1
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e1.LE.spm1bhi.AND.alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s1-rm1b**2)**2+
     &       (gam1b*rm1b)**2)/(gam1b*rm1b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e1
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF
      xjac=xjac/twopi
      
* s11 = (p1+p2)**2 = q11**2
      e_sup = sqrt(s1)-xm(io3)
      e_inf = xm(io1)+xm(io2)      
      IF(e_sup.LE.e_inf)THEN    
        xjac=0.d0
        RETURN
      ENDIF
      CALL PHSP_INI5(e_sup,e_inf,
     &               spm11alo,spm11ahi,rm11a,gam11a,
     &               spm11blo,spm11bhi,rm11b,gam11b,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max)
      e11=sqrt(s11)
      IF(e11.LE.spm11alo .AND. alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e11/alfa1
      ELSEIF(e11.LE.spm11ahi .AND. alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s11-rm11a**2)**2+
     &       (gam11a*rm11a)**2)/(gam11a*rm11a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e11.LE.spm11blo .AND. alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e11
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e11.LE.spm11bhi .AND. alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s11-rm11b**2)**2+
     &       (gam11b*rm11b)**2)/(gam11b*rm11b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e11
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF
      xjac=xjac/twopi

* s2 = (p4+p5)**2 = q2**2
      e_sup = sqrt(saux)-sqrt(s1)
      e_inf = xm(io4)+xm(io5)      
      IF(e_sup.LE.e_inf)THEN    
        xjac=0.d0
        RETURN
      ENDIF
      CALL PHSP_INI5(e_sup,e_inf,
     &               spm2alo,spm2ahi,rm2a,gam2a,
     &               spm2blo,spm2bhi,rm2b,gam2b,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max)
      e2=sqrt(s2)
      IF(e2.LE.spm2alo.AND.alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e2/alfa1
      ELSEIF(e2.LE.spm2ahi.AND.alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s2-rm2a**2)**2+
     &      (gam2a*rm2a)**2)/(gam2a*rm2a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e2.LE.spm2blo.AND.alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e2
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e2.LE.spm2bhi.AND.alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s2-rm2b**2)**2+
     &      (gam2b*rm2b)**2)/(gam2b*rm2b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e2
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF      
      xjac=xjac/twopi


* q1 and q2
      app=((saux-s1-s2)**2-4.d0*s1*s2)/(4.d0*saux)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      pauxmod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*pauxmod/sqrt(saux)/(fourpi)**2
      xjac=xjac*fourpi

* q11 and p3
      app=((s1-s11-xmsq(io3))**2-4.d0*s11*xmsq(io3))/(4.d0*s1)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p1mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p1mod/sqrt(s1)/(fourpi)**2
      xjac=xjac*fourpi

* q11 --> p1 and p2      
      app=((s11-xmsq(io1)-xmsq(io2))**2
     &              -4.d0*xmsq(io1)*xmsq(io2))/(4.d0*s11)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p11mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p11mod/sqrt(s11)/(fourpi)**2
      xjac=xjac*fourpi

* q2 --> p4 and p5
      app=((s2-xmsq(io4)-xmsq(io5))**2
     &              -4.d0*xmsq(io4)*xmsq(io5))/(4.d0*s2)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p2mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p2mod/sqrt(s2)/(fourpi)**2
      xjac=xjac*fourpi

      RETURN
      END
      

