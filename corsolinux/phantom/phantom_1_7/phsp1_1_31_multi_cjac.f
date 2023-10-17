c
c Last update: Jul 07, 2007
c
************************************************************************
* Subroutine to compute the jacobian corresponding to the phase space 
*  mapping:
*          [1t]+[1t]+[4(rm1[a,b]) --> 3(rm11[a,b])+p4 
*                                      --> 2(rm111[a,b])+p3+p4] 
* where rm1[a,b],rm11[a,b],rm111[a,b] are the masses of the intermediate
*  resonances.
* rm1, rm11 and rm111 can resonate at two different masses (e.g. Z and 
*  H) with rm1a < rm1b, rm11a < rm11b and rm111a < rm111b.
* Two exponents (expa,2) are available for mapping the two t-channel 
*  propagators.
* rmta, rmtb are the corresponding t-channel masses.
* These mapping are applied for 1.d0>cos(theta)>cos_sep. cos_sep is a 
*  PARAMETER
* Outside this range the cosine is generated flat. 
* xmasq,xmbsq are the masses squared of the two incoming particles 
*     (a --> 5, b --> 6 in the t's)
* pa=sqrt(shat)/2*(1,0,0,1) pb=sqrt(shat)/2*(1,0,0,-1) in the massless 
*   limit.
*
* INPUT:
*  shat          : total energy squared
*  xm(8),xmsq(8) : masses and masses squared of the particles
*                        in the order used in the main program      
*  rmN,gamN      : N=1,11,111 mass and width for intermediate 
*                   resonance N
*  spmNlo,spmNhi : N=111,11,1 low and high limit of resonance mapping.
*                        outside this range the mapping is flat
*  pt5min_in(3)  : array of possible minimum pt's of particle 5, 
*                  depending on its flavour. The first element refers 
*                  to jets, the second one to charged leptons and the 
*                  third one to neutrinos.
*  pt6min_in(3)  : array of possible minimum pt's of particle 6,
*                  depending on its flavour. The first element refers 
*                  to jets, the second one to charged leptons and the 
*                  third one to neutrinos.
*  id5           : idp(iorder(7)) used to determine which minimum pt
*                  must be considered in the t-channel.
*  id6           : idp(iorder(8)) used to determine which minimum pt
*                  must be considered in the t-channel.
*  iorder(6)     : specifies the permutation that relates internal and
*                        external momenta: p = p(i,iorder(1))
*  xmasq,xmbsq   : masses squared of the two incoming particles 
*                        (a --> 5, b --> 6 in the t's)
*  p(0:3,8)      : momenta in the hard scattering CM.
*                        Includes initial momenta.
* OUTPUT:
*  xjac          : full jacobian
************************************************************************

        subroutine phsp1_1_31_multi_cjac(shat,
     &                    rm1a,gam1a,spm1alo,spm1ahi,
     &                    rm1b,gam1b,spm1blo,spm1bhi,
     &                    rm11a,gam11a,spm11alo,spm11ahi,
     &                    rm11b,gam11b,spm11blo,spm11bhi,
     &                    rm111a,gam111a,spm111alo,spm111ahi,
     &                    rm111b,gam111b,spm111blo,spm111bhi,
     &                    rmta,expa,spta,rmtb,expb,sptb,
c giuseppe 07/07/2007
c     &                    pt5min,pt6min,
     &                    pt5min_in,pt6min_in,
     &                    id5,id6,
c end giuseppe 07/07/2007
     &                    xmasq,xmbsq,
     &                    xm,xmsq,iorder,
     &                    p,xjac)
     
      IMPLICIT NONE
      REAL*8 shat,
     &       rm1a,gam1a,spm1alo,spm1ahi,
     &       rm1b,gam1b,spm1blo,spm1bhi,
     &       rm11a,gam11a,spm11alo,spm11ahi,
     &       rm11b,gam11b,spm11blo,spm11bhi,
     &       rm111a,gam111a,spm111alo,spm111ahi,
     &       rm111b,gam111b,spm111blo,spm111bhi,
     &       rmta,expa,rmtb,expb,
c giuseppe 07/07/2007
     &       pt5min_in(3),pt6min_in(3),
c end giuseppe 07/07/2007
     &       pt5min,pt6min,
     &       xmasq,xmbsq,
     &       xm(8),xmsq(8),
     &       p(0:3,8),xjac
c giuseppe 07/07/2007
      INTEGER iorder(6),id5,id6
c end giuseppe 07/07/2007

* Local variables
      REAL*8 e_sup,e_inf,x_max,x_min,app,
     &       s111,s11,s1,pmod_15,
     &       e111,e11,e1,p111mod,p11mod,p1mod,pmod,
     &       saux,eaux,t_inf,t_sup,t1,t2,pauxmod,pbauxmod,p0mod,
     &       t_inf0,t_sup0,
     &       pi,twopi,fourpi,
     &       rmtasq,rmtbsq,t1aux
     
      INTEGER io1,io2,io3,io4,io5,io6,i
      

      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi)
      
* Mutimapping BW
      REAL*8 alfa1,alfa2,auxalfa2,
     &       alfa3,auxalfa3,alfa4,auxalfa4,
     &       xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &       xw2_min,xw2_max,xf3_min,xf3_max

* Mutimapping t
      REAL*8 alfa,xt1_min,xt1_max,xt2_min,xt2_max,spta,sptb
      REAL*8 cos_sep
      PARAMETER(cos_sep=0.4d0)
      
      io1=iorder(1)
      io2=iorder(2)
      io3=iorder(3)
      io4=iorder(4)
      io5=iorder(5)
      io6=iorder(6)

      xjac=1.d0


c giuseppe 07/07/2007
c Determine which minimum pt must be assigned to t-channel, depending
c on final-state fermion's flavour
      IF(abs(id5).LE.5 .OR. id5.EQ.21)THEN
        pt5min=pt5min_in(1)
      ELSEIF(abs(id5).EQ.11)THEN
        pt5min=pt5min_in(2)
      ELSEIF(abs(id5).EQ.12)THEN
        pt5min=pt5min_in(3)
      ELSE
        print*,'***phsp1_1_31_jac ERROR'
        STOP
      ENDIF
      IF(abs(id6).LE.5 .OR. id6.EQ.21)THEN
        pt6min=pt6min_in(1)
      ELSEIF(abs(id6).EQ.11)THEN
        pt6min=pt6min_in(2)
      ELSEIF(abs(id6).EQ.12)THEN
        pt6min=pt6min_in(3)
      ELSE
        print*,'***phsp1_1_31_jac ERROR'
        STOP
      ENDIF
c end giuseppe 07/07/2007


C All these invariant masses should be computed once and for all and
C stored in some array.

      s111=(p(0,io1)+p(0,io2))**2
      s11=(p(0,io1)+p(0,io2)+p(0,io3))**2
      s1=(p(0,io1)+p(0,io2)+p(0,io3)+p(0,io4))**2
      saux=(p(0,io1)+p(0,io2)+p(0,io3)+p(0,io4)+p(0,io6))**2
      DO i=1,3
        s111=s111-(p(i,io1)+p(i,io2))**2
        s11=s11-(p(i,io1)+p(i,io2)+p(i,io3))**2
        s1=s1-(p(i,io1)+p(i,io2)+p(i,io3)+p(i,io4))**2
        saux=saux-(p(i,io1)+p(i,io2)+p(i,io3)+p(i,io4)+p(i,io6))**2
      ENDDO 
      IF( s111.LE.1.d-8 .OR. s11.LE.1.d-8 .OR. s1.LE.1.d-8
     &             .OR. saux.LE.1.d-8 )THEN
        xjac=0.d0
        RETURN
      ENDIF

* saux = (q1+p6)**2 = qaux**2: flat
      e_sup = sqrt(shat)-xm(io5)
      e_inf = xm(io1)+xm(io2)+xm(io3)+xm(io4)+xm(io6)  
c      e_inf = sqrt(s1)+xm(io6)    
      eaux=sqrt(saux)
      xjac=xjac*(e_sup-e_inf)*2.d0*eaux
      xjac=xjac/twopi
      
* pa+pb --> p5 + qaux; (pa-p5)**2 = (pb-qaux)**2 = t1
      rmtasq=rmta*rmta
      rmtbsq=rmtb*rmtb    
      app=((shat-saux-xmsq(io5))**2-4.d0*saux*xmsq(io5))/(4.d0*shat)
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
     
      IF(pt5min.LE.0.d0)THEN
        t_inf =xmbsq+saux-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                                    (shat+saux-xmsq(io5))
     &               -2.d0*p0mod*pmod
        app=(xmbsq-saux)*(xmasq-xmsq(io5))
     &           +(saux-xmbsq+xmasq-xmsq(io5))*
     &                        (saux*xmasq-xmbsq*xmsq(io5))/shat
        t_sup=app/t_inf
      ELSE
        IF(pt5min.GE.pmod)THEN
          xjac=0.d0
          RETURN
        ENDIF
        app=sqrt((pmod-pt5min)*(pmod+pt5min))
        t_inf =xmbsq+saux-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                                    (shat+saux-xmsq(io5))
     &               -2.d0*p0mod*app
        t_sup =xmbsq+saux-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                                    (shat+saux-xmsq(io5))
     &               +2.d0*p0mod*app
      ENDIF
      
      spta=xmbsq+saux-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                                    (shat+saux-xmsq(io5))
     &               +2.d0*p0mod*pmod*cos_sep

* From now on we use t1 = (pa-p5)**2 - rmtasq = (pb-qaux)**2 - rmtasq
      t_inf=t_inf-rmtasq
      t_sup=t_sup-rmtasq
      spta=spta-rmtasq
      CALL T_INI(t_inf,t_sup,spta,expa,
     &                 alfa,xt1_min,xt1_max,xt2_min,xt2_max)
      t1=xmasq+xmsq(io5)-0.5d0/shat*(shat+xmasq-xmbsq)*
     &                                    (shat+xmsq(io5)-saux)
     &               +2.d0*p0mod*p(3,io5)-rmtasq
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

* qaux and p5
* beta/2/(4*pi)**2
      xjac=xjac*pmod/sqrt(shat)/(fourpi)**2
      xjac=xjac*twopi

* s1 = (q11+q12)**2 = q1**2
c      e_sup = sqrt(shat)-xm(io5)-xm(io6)
c      e_inf = sqrt(s11)+ sqrt(s12)    
      e_sup = sqrt(saux)-xm(io6)
      e_inf = xm(io1)+xm(io2)+xm(io3)+xm(io4)   
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
      ELSEIF(e1.LE.spm1ahi .AND. alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s1-rm1a**2)**2+
     &       (gam1a*rm1a)**2)/(gam1a*rm1a)/(1.d0-alfa1)/auxalfa2        
      ELSEIF(e1.LE.spm1blo .AND. alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e1
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e1.LE.spm1bhi .AND. alfa4.GT.0.d0)THEN
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
      
* q1 and q6
      app=((saux-s1-xmsq(io6))**2-4.d0*s1*xmsq(io6))/(4.d0*saux)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      pauxmod=sqrt(app)
      t1aux=t1+rmtasq
      app=((saux-t1aux-xmbsq)**2-4.d0*t1aux*xmbsq)/(4.d0*saux)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      pbauxmod=sqrt(app)
      
*Define theta, phi as a function of t2 = (pb-p6)**2 = ((pa-p5)-q1)**2
      t_inf=t1aux+s1-0.5d0/saux*(saux+t1aux-xmbsq)*(saux+s1-xmsq(io6))
     &        -2.d0*pauxmod*pbauxmod
      app=(t1aux-s1)*(xmbsq-xmsq(io6))
     &     +(s1-t1aux+xmbsq-xmsq(io6))*(s1*xmbsq-t1aux*xmsq(io6))/saux
      t_sup=app/t_inf
* In order to get the exact limit of integration for t2 from pt6min we 
* should know the invariant mass (q1+p5)**2. 
* However (q1+p5)**2 > s1+xmsq(io5)**2. Hence we can establish bounds 
* on t2 which are looser than the exact ones.
      IF(pt6min.gt.0.d0)THEN
        app=((shat-(s1+xmsq(io5))-xmsq(io6))**2
     &                   -4.d0*(s1+xmsq(io5))*xmsq(io6))/(4.d0*shat)
        IF(app.LE.pt6min*pt6min)THEN
          xjac=0.d0
          RETURN
        ENDIF
        pmod_15=sqrt(app)
        app=sqrt((pmod_15-pt6min)*(pmod_15+pt6min))
        
        t_inf0 =xmbsq+xmsq(io6)-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                             (shat+xmsq(io6)-(s1+xmsq(io5)))
     &               -2.d0*p0mod*app
        t_sup0 =xmbsq+xmsq(io6)-0.5d0/shat*(shat+xmbsq-xmasq)*
     &                             (shat+xmsq(io6)-(s1+xmsq(io5)))
     &               +2.d0*p0mod*app
        t_inf=max(t_inf,t_inf0)
        t_sup=min(t_sup,t_sup0)
      ENDIF

      sptb=t1aux+s1-0.5d0/saux*(saux+t1aux-xmbsq)*(saux+s1-xmsq(io6))
     &        +2.d0*pauxmod*pbauxmod*cos_sep
     
* From now on we use t2 = (pb-p6)**2 - rmtbsq = ((pa-p5)-q1)**2 - rmtbsq
      t_inf=t_inf-rmtbsq
      t_sup=t_sup-rmtbsq
      sptb=sptb-rmtbsq
      CALL T_INI(t_inf,t_sup,sptb,expb,
     &                 alfa,xt1_min,xt1_max,xt2_min,xt2_max)
      t2=xmbsq+xmsq(io6)
     &      -2.d0*sqrt(p0mod*p0mod+xmbsq)*p(0,io6)
     &                    -2.d0*p0mod*p(3,io6)-rmtbsq
      IF(t2.LE.sptb)THEN
        xjac=xjac*0.5d0*(xt1_max-xt1_min)/pauxmod/pbauxmod
        xjac=xjac/alfa
      ELSE 
        IF(expb.eq.0.d0)THEN 
          xjac=xjac*0.5d0*(xt2_max-xt2_min)/pauxmod/pbauxmod
        ELSEIF(expb.ne.1.d0)THEN
          xjac=xjac*0.5d0*(xt2_max-xt2_min)
     &       /(1.d0-expb)/pauxmod/pbauxmod*(-t2)**expb
        ELSE
          xjac=xjac*0.5d0*(xt2_max-xt2_min)/pauxmod/pbauxmod*(-t2)
        ENDIF
        xjac=xjac/(1.d0-alfa)
      ENDIF
* beta/2/(4*pi)**2
      xjac=xjac*pauxmod/sqrt(saux)/(fourpi)**2
      xjac=xjac*twopi

* s11 = (p1+p2+p3)**2 = (q111+p3)**2 = q11**2
      e_sup = sqrt(s1)-xm(io4)
      e_inf = xm(io1)+xm(io2)+xm(io3)     
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

* s111 = (p1+p2)**2 = q111**2
      e_sup = sqrt(s11)-xm(io3)
      e_inf = xm(io1)+xm(io2)      
      IF(e_sup.LE.e_inf)THEN    
        xjac=0.d0
        RETURN
      ENDIF
      CALL PHSP_INI5(e_sup,e_inf,
     &               spm111alo,spm111ahi,rm111a,gam111a,
     &               spm111blo,spm111bhi,rm111b,gam111b,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max)
      e111=sqrt(s111)
      IF(e111.LE.spm111alo .AND. alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e111/alfa1
      ELSEIF(e111.LE.spm111ahi .AND. alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s111-rm111a**2)**2+
     &       (gam111a*rm111a)**2)/(gam111a*rm111a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e111.LE.spm111blo .AND. alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e111
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e111.LE.spm111bhi .AND. alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s111-rm111b**2)**2+
     &       (gam111b*rm111b)**2)/(gam111b*rm111b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e111
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF
      xjac=xjac/twopi
      
* q11 and p4
      app=((s1-s11-xmsq(io4))**2-4.d0*s11*xmsq(io4))/(4.d0*s1)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p1mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p1mod/sqrt(s1)/(fourpi)**2
      xjac=xjac*fourpi

* q11 --> p111 and p3      
      app=((s11-s111-xmsq(io3))**2
     &              -4.d0*s111*xmsq(io3))/(4.d0*s11)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p11mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p11mod/sqrt(s11)/(fourpi)**2
      xjac=xjac*fourpi

* q111 --> p1 and p2
      app=((s111-xmsq(io1)-xmsq(io2))**2
     &              -4.d0*xmsq(io1)*xmsq(io2))/(4.d0*s111)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p111mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p111mod/sqrt(s111)/(fourpi)**2
      xjac=xjac*fourpi

      RETURN
      END
      

