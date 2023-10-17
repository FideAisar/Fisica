c
c Last update:  May 29, 2014
c
************************************************************************
* Subroutine to compute the momenta p1,....,p6 of six particles of mass
* xm1,...,xm6 given a point in the 14 dimensional hypercube [0,1]^14
* Returns also xjac, the weight due to the jacobian.
*
* The phase space mapping corresponds to
*          [1t]+[1t]+[4(rm1[a,b]) --> 2(rm11[a,b,c],p1,p2)+2(rm12[a,b,c],p3,p4)]
* where rm1,rm11,rm12 are the masses of the intermediate resonances
* rm1 can resonate at two different masses
*    (e.g. H, H') with rm1a < rm1b.
* rm11 and rm12 can resonate at three different masses
*    (e.g. Z, H, H') with
* rmXa < rmXb < rmXc.
* If only one resonance is needed set the corresponding
*   spmXahi=spmXalo=0.d0
*
* Two exponents (expa,b) are available for mapping the t-channel
*   propagators.
* rmta, rmtb are the corresponding t-channel masses.
* These mapping are applied for 1.d0>cos(theta)>cos_sep. cos_sep is
*   a PARAMETER
* Outside this range the cosine is generated flat.
*y xmasq,xmbsq are the masses squared of the two incoming particles
*   (a --> 5, b --> 6 in  the t's)
* pa=sqrt(shat)/2*(1,0,0,1) pb=sqrt(shat)/2*(1,0,0,-1) in the massless
*  limit.
*
* If in the main x had dimension 15 with x(1) and x(2) fixing the
*  initial state particle momentum fraction the subroutine can be
*  called as  call phsp1_1_4(x(3),shat,.......)
*

* INPUT:

*  x,shat        : vegas random points, total energy squared
*  xm(8),xmsq(8) : masses and masses squared of the particles
*                        in the order used in the main program
*  rmN,gamN      : N=1[a,b],11[a,b,c],12[a,b,c] mass and width for intermediate
*                        resonance N.
*  spmNlo,spmNhi : N=11[a,b,c],12[a,b,c],1[a,b] low and high limit of resonance
*                        mapping. Outside this range the mapping is flat
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
*  xmasq,xmbsq   : masses of the two incoming particles
*                        (a --> 5, b --> 6 in the t's)

* OUTPUT:

*  p(0:3,8)      : momenta in the hard scattering CM.
*                        Includes initial momenta.
*  xjac          : full jacobian
************************************************************************

        subroutine phsp1_1_4_multi7_c(x,shat,
     &                    rm1a,gam1a,spm1alo,spm1ahi,
     &                    rm1b,gam1b,spm1blo,spm1bhi,
     &                    rm11a,gam11a,spm11alo,spm11ahi,
     &                    rm11b,gam11b,spm11blo,spm11bhi,
     &                    rm11c,gam11c,spm11clo,spm11chi,
     &                    rm12a,gam12a,spm12alo,spm12ahi,
     &                    rm12b,gam12b,spm12blo,spm12bhi,
     &                    rm12c,gam12c,spm12clo,spm12chi,
     &                    rmta,expa,spta,rmtb,expb,sptb,
     &                    pt5min_in,pt6min_in,
     &                    id5,id6,
     &                    xmasq,xmbsq,
     &                    xm,xmsq,iorder,
     &                    p,xjac)

      IMPLICIT NONE
      REAL*8 x(13),shat,
     &       rm1a,gam1a,spm1alo,spm1ahi,
     &       rm1b,gam1b,spm1blo,spm1bhi,
     &       rm11a,gam11a,spm11alo,spm11ahi,
     &       rm11b,gam11b,spm11blo,spm11bhi,
     &       rm11c,gam11c,spm11clo,spm11chi,
     &       rm12a,gam12a,spm12alo,spm12ahi,
     &       rm12b,gam12b,spm12blo,spm12bhi,
     &       rm12c,gam12c,spm12clo,spm12chi,
     &       rmta,expa,rmtb,expb,
     &       pt5min_in(3),pt6min_in(3),
     &       pt5min,pt6min,
     &       xmasq,xmbsq,
     &       xm(8),xmsq(8),
     &       p(0:3,8),xjac,xjacvar,svar
      INTEGER iorder(6),id5,id6

* Local variables
      REAL*8 e_sup,e_inf,x_max,x_min,app,
     &       s11,s12,s1,q11(0:3),q12(0:3),q1(0:3),qaux(0:3),
     &       pbaux(0:3),p6aux(0:3),
     &       e11,e12,e1,p11mod,p12mod,p1mod,pmod,pmod_15,
     &       saux,eaux,t_inf,t_sup,t1,t2,pauxmod,pbauxmod,p0mod,
     &       t_inf0,t_sup0,
     &       pi,twopi,fourpi,
     &       rmtasq,rmtbsq,costh,sinth,ph,ran2,t1aux

      INTEGER idum,io1,io2,io3,io4,io5,io6,i,iregion
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi)

* Mutimapping
      REAL*8 alfa1,alfa2,auxalfa2,
     &       alfa3,auxalfa3,alfa4,auxalfa4,
     &       alfa5,auxalfa5,alfa6,auxalfa6,
     &       xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &       xw2_min,xw2_max,xf3_min,xf3_max,
     &       xw3_min,xw3_max,xf4_min,xf4_max,
     &       vegas_x
* Mutimapping t
      REAL*8 alfa,xt1_min,xt1_max,xt2_min,xt2_max,spta,sptb
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

      xjac=1.d0


c Determine which minimum pt must be assigned to t-channel, depending
c on final-state fermion's flavour
      IF(abs(id5).LE.5 .OR. id5.EQ.21)THEN
        pt5min=pt5min_in(1)
      ELSEIF(abs(id5).EQ.11)THEN
        pt5min=pt5min_in(2)
      ELSEIF(abs(id5).EQ.12)THEN
        pt5min=pt5min_in(3)
      ELSE
        print*,'***phsp1_1_4 ERROR'
        STOP
      ENDIF
      IF(abs(id6).LE.5 .OR. id6.EQ.21)THEN
        pt6min=pt6min_in(1)
      ELSEIF(abs(id6).EQ.11)THEN
        pt6min=pt6min_in(2)
      ELSEIF(abs(id6).EQ.12)THEN
        pt6min=pt6min_in(3)
      ELSE
        print*,'***phsp1_1_4 ERROR'
        STOP
      ENDIF


* saux = (q1+p6)**2 = qaux**2: flat
      e_sup = sqrt(shat)-xm(io5)
c      e_inf = sqrt(s1)+xm(io6)
      e_inf = xm(io1)+xm(io2)+xm(io3)+xm(io4)+xm(io6)
c      eaux=(e_sup-e_inf)*x(4)+e_inf
      eaux=(e_sup-e_inf)*x(1)+e_inf
      saux=eaux*eaux
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
        app=(pmod-pt5min)*(pmod+pt5min)
        IF(app.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        app=sqrt(app)
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
     &       /(1.d0-expa)/p0mod/pmod*(-t1)**expa
        ELSE
          t1=-exp((xt2_max-xt2_min)*vegas_x+xt2_min)
	  xjac=xjac*0.5d0*(xt2_max-xt2_min)/p0mod/pmod*(-t1)
        ENDIF
        xjac=xjac/(1.d0-alfa)
      ENDIF

      costh=0.5*(t1+rmtasq-xmbsq-saux
     &        +0.5d0/shat*(shat+xmbsq-xmasq)*(shat+saux-xmsq(io5)))
     &       /p0mod/pmod

* qaux and p5
      costh=min(costh,1.d0)
      costh=max(costh,-1.d0)
      sinth=sqrt(1.d0-costh*costh)
      ph=ran2(idum)*twopi
      p(0,io5)=sqrt(pmod*pmod+xmsq(io5))
      p(1,io5)=pmod*sinth*cos(ph)
      p(2,io5)=pmod*sinth*sin(ph)
      p(3,io5)=pmod*costh
      qaux(0)=sqrt(pmod*pmod+saux)
      qaux(1)=-p(1,io5)
      qaux(2)=-p(2,io5)
      qaux(3)=-p(3,io5)
* beta/2/(4*pi)**2
      xjac=xjac*pmod/sqrt(shat)/(fourpi)**2
      xjac=xjac*twopi

* s1 = (q11+q12)**2 = q1**2
      e_sup = sqrt(saux)-xm(io6)
      e_inf = xm(io1)+xm(io2)+xm(io3)+xm(io4)
      IF( e_inf.GE.e_sup)THEN
        xjac=0.d0
        RETURN
      ENDIF
      iregion=1
      CALL PHSP_INI5(e_sup,e_inf,
     &               spm1alo,spm1ahi,rm1a,gam1a,
     &               spm1blo,spm1bhi,rm1b,gam1b,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max)
      CALL PHSP_5(x(3),
     &               rm1a,gam1a,rm1b,gam1b,
     &                     alfa1,alfa2,auxalfa2,
     &                     alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &                     iregion,
     &                     svar,xjacvar)
      s1=svar
      xjac=xjac*xjacvar

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

**** boost pb to paux  cm: pb=(eb,0,0,-p0mod)
      pbaux(0)=0.5d0*(shat-xmasq+xmbsq)/sqrt(shat)
      pbaux(1)=0.d0
      pbaux(2)=0.d0
      pbaux(3)=-p0mod
      CALL boostinv(pbaux,qaux,pbaux)
*Define theta, phi as a function of t2 = (pb-p6)**2 = ((pa-p5)-q1)**2
      t_inf=t1aux+s1-0.5d0/saux*(saux+t1aux-xmbsq)*(saux+s1-xmsq(io6))
     &        -2.d0*pauxmod*pbauxmod
      app=(t1aux-s1)*(xmbsq-xmsq(io6))
     &     +(s1-t1aux+xmbsq-xmsq(io6))*(s1*xmbsq-t1aux*xmsq(io6))/saux
      t_sup=app/t_inf
* In order to get the exact limit of integration for t2 from pt6min we should
* know the invariant mass (q1+p5)**2. However (q1+p5)**2 > s1+xmsq(io5)**2.
* Hence we can establish bounds on t2 which are looser than the exact ones.
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

* From now on we use t2 =  (pb-p6)**2 - rmtbsq = ((pa-p5)-q1)**2 - rmtbsq
      t_inf=t_inf-rmtbsq
      t_sup=t_sup-rmtbsq
      sptb=sptb-rmtbsq
      CALL T_INI(t_inf,t_sup,sptb,expb,
     &                 alfa,xt1_min,xt1_max,xt2_min,xt2_max)
      IF(x(4).LE.alfa .AND. alfa.GT.0.d0)THEN
        vegas_x=x(4)/alfa
        t2=(xt1_max-xt1_min)*vegas_x+xt1_min
        xjac=xjac*0.5d0*(xt1_max-xt1_min)/pauxmod/pbauxmod
        xjac=xjac/alfa
      ELSE
        vegas_x=(x(4)-alfa)/(1.d0-alfa)
        IF(expb.eq.0.d0)THEN
          t2=(xt2_max-xt2_min)*vegas_x+xt2_min
          xjac=xjac*0.5d0*(xt2_max-xt2_min)/pauxmod/pbauxmod
        ELSEIF(expb.ne.1.d0)THEN
          t2=-((xt2_max-xt2_min)*vegas_x+xt2_min)**(1.d0/(1.d0-expb))
          xjac=xjac*0.5d0*(xt2_max-xt2_min)
     &       /(1.d0-expb)/pauxmod/pbauxmod*(-t2)**expb
        ELSE
          t2=-exp((xt2_max-xt2_min)*vegas_x+xt2_min)
	  xjac=xjac*0.5d0*(xt2_max-xt2_min)/pauxmod/pbauxmod*(-t2)
        ENDIF
        xjac=xjac/(1.d0-alfa)
      ENDIF
      costh=0.5*(t2+rmtbsq-t1aux-s1
     &        +0.5d0/saux*(saux+t1aux-xmbsq)*(saux+s1-xmsq(io6)))
     &       /pauxmod/pbauxmod
      costh=min(costh,1.d0)
      costh=max(costh,-1.d0)
      sinth=sqrt(1.d0-costh*costh)
      ph=x(5)*twopi

* p6aux = q6 in qaux cm
      p6aux(0)=sqrt(pauxmod*pauxmod+xmsq(io6))
      p6aux(1)=pauxmod*sinth*cos(ph)
      p6aux(2)=pauxmod*sinth*sin(ph)
      p6aux(3)=pauxmod*costh

* Rotate p6aux to use pbaux as new z axis
      CALL rot(p6aux(1),pbaux(1),p6aux(1))
* Boost p6aux to overall cm
      CALL boost(p6aux,qaux,p(0,io6))
      DO i=0,3
        q1(i)=qaux(i)-p(i,io6)
      ENDDO
* beta/2/(4*pi)**2
      xjac=xjac*pauxmod/sqrt(saux)/(fourpi)**2
      xjac=xjac*twopi

* s11 = (p1+p2)**2 = q11**2
      e_sup = sqrt(s1)-xm(io3)-xm(io4)
      e_inf = xm(io1)+xm(io2)
      IF( e_inf.GE.e_sup)THEN
        xjac=0.d0
        RETURN
      ENDIF
      iregion=3
      CALL PHSP_INI7(e_sup,e_inf,
     &               spm11alo,spm11ahi,rm11a,gam11a,
     &               spm11blo,spm11bhi,rm11b,gam11b,
     &               spm11clo,spm11chi,rm11c,gam11c,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               alfa5,auxalfa5,alfa6,auxalfa6,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &               xw3_min,xw3_max,xf4_min,xf4_max)
      CALL PHSP_7(x(6),
     &               rm11a,gam11a,rm11b,gam11b,rm11c,gam11c,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               alfa5,auxalfa5,alfa6,auxalfa6,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &               xw3_min,xw3_max,xf4_min,xf4_max,
     &                     iregion,
     &                     svar,xjacvar)
      s11=svar
      xjac=xjac*xjacvar

* s12 = (p3+p4)**2 = q12**2
      e_sup = sqrt(s1)-sqrt(s11)
      e_inf = xm(io3)+xm(io4)
      IF( e_inf.GE.e_sup)THEN
        xjac=0.d0
        RETURN
      ENDIF
      iregion=4
      CALL PHSP_INI7(e_sup,e_inf,
     &               spm12alo,spm12ahi,rm12a,gam12a,
     &               spm12blo,spm12bhi,rm12b,gam12b,
     &               spm12clo,spm12chi,rm12c,gam12c,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               alfa5,auxalfa5,alfa6,auxalfa6,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &               xw3_min,xw3_max,xf4_min,xf4_max)
      CALL PHSP_7(x(7),
     &               rm12a,gam12a,rm12b,gam12b,rm12c,gam12c,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               alfa5,auxalfa5,alfa6,auxalfa6,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &               xw3_min,xw3_max,xf4_min,xf4_max,
     &                     iregion,
     &                     svar,xjacvar)
      s12=svar
      xjac=xjac*xjacvar

* q11 and q12
      app=((s1-s11-s12)**2-4.d0*s11*s12)/(4.d0*s1)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p1mod=sqrt(app)
      CALL TWOBODY(q1,s1,p1mod,s11,s12,x(8),x(9),q11,q12)
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
      CALL TWOBODY(q11,s11,p11mod,xmsq(io1),xmsq(io2),x(10),x(11),
     &             p(0,io1),p(0,io2))
* beta/2/(4*pi)**2
      xjac=xjac*p11mod/sqrt(s11)/(fourpi)**2
      xjac=xjac*fourpi

* q12 --> p3 and p4
      app=((s12-xmsq(io3)-xmsq(io4))**2
     &              -4.d0*xmsq(io3)*xmsq(io4))/(4.d0*s12)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p12mod=sqrt(app)
      CALL TWOBODY(q12,s12,p12mod,xmsq(io3),xmsq(io4),x(12),x(13),
     &             p(0,io3),p(0,io4))
* beta/2/(4*pi)**2
      xjac=xjac*p12mod/sqrt(s12)/(fourpi)**2
      xjac=xjac*fourpi

      RETURN
      END


