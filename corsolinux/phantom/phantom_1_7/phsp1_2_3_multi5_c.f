c
c Last update: Sep 29, 2007
c
************************************************************************
* Subroutine to compute the momenta p1,....,p6 of six particles of mass 
* xm1,...,xm6 given a point in the 14 dimensional hypercube [0,1]^14
* Returns also xjac, the weight due to the jacobian.
*
* The phase space mapping corresponds to 
*      [1t,a-->p6]+ 5(rm0[a,b])
*                    -->  [ [3(rm1[a,b]) --> 2(rm11[a,b],p1,p2)+p3] +
*                          + 2(rm2[a,b],p4,p5) ]
* where rm0[a,b],rm1[a,b],rm11[a,b],rm2[a,b] are the masses of the 
*  intermediate resonances
*
* rm0,rm1, rm11 and rm2 can resonate at two different masses (e.g. Z 
*  and H) with rm0a < rm0b, rm1a < rm1b, rm11a < rm11b, rm2a < rm2b.
*
* One exponent (expa) is available for mapping the t-channel propagator.
* rmta is the corresponding t-channel mass.
* This mapping is applied for 1.d0>cos(theta)>cos_sep. cos_sep is 
*  a PARAMETER
* Outside this range the cosine is generated flat. 
*
* If in the main x had dimension 15 with x(1) and x(2) fixing the 
*  initial state particle momentum fraction the subrotine can be called 
*  as call phsp1_1_4(x(3),shat,.......)
* 

* INPUT:

*  x,shat        : vegas random points, total energy squared
*  xm(8),xmsq(8) : masses and masses squared of the particles
*                   in the order used in the main program      
*  rmN,gamN      : N=0[a,b],1[a,b],11[a,b],2[a,b] mass and width for 
*                   intermediate resonance N.
*  spmNlo,spmNhi : N=0[a,b],1[a,b],11[a,b],2[a,b] low and high limit of 
*                   resonance mapping.
*                   Outside this range the mapping is flat
*  pt6min_in(3)  : array of possible minimum pt's of particle 6,
*                  depending on its flavour. The first element refers 
*                  to jets, the second one to charged leptons and the 
*                  third one to neutrinos.
*  iorder(6)     : specifies the permutation that relates internal and
*                   external momenta: p = p(i,iorder(1))
*  id1,id2       : idp(iorder(1)),idp(iorder(2)) used to determine the
*                   direction of the incoming particles involved in the
*                   t-channel 
*  id6           : idp(iorder(8)) used to determine which external
*                   particle enters in the t-channel exchange.
*                   Useful in case of identical particles in the final 
*                   state.
* xmasq_in,xmbsq_in : masses squared of the two incoming particles
 
* OUTPUT:

*  p(0:3,8)      : momenta in the hard scattering CM.
*                        Includes initial momenta.
*  xjac          : full jacobian
*
* INTERNAL VARIABLES
* xmasq,xmbsq are the masses squared of the two incoming particles 
* (a --> 6  in the t)
* The routine determines from iorder whether
* pa=sqrt(shat)/2*(1,0,0,1), pb=sqrt(shat)/2*(1,0,0,-1) or viceversa
* in the massless limit.
************************************************************************

        subroutine phsp1_2_3_multi5_c(x,shat,
     &                    rm0a,gam0a,spm0alo,spm0ahi,
     &                    rm0b,gam0b,spm0blo,spm0bhi,
     &                    rm1a,gam1a,spm1alo,spm1ahi,
     &                    rm1b,gam1b,spm1blo,spm1bhi,
     &                    rm11a,gam11a,spm11alo,spm11ahi,
     &                    rm11b,gam11b,spm11blo,spm11bhi,
     &                    rm2a,gam2a,spm2alo,spm2ahi,
     &                    rm2b,gam2b,spm2blo,spm2bhi,
     &                    rmta,expa,spta,
     &                    pt6min_in,
     &                    id1,id2,id6,
     &                    xmasq_in,xmbsq_in,
     &                    xm,xmsq,iorder,
     &                    p,xjac)
     
      IMPLICIT NONE
      REAL*8 x(13),shat,
     &       rm0a,gam0a,spm0alo,spm0ahi,
     &       rm0b,gam0b,spm0blo,spm0bhi,
     &       rm1a,gam1a,spm1alo,spm1ahi,
     &       rm1b,gam1b,spm1blo,spm1bhi,
     &       rm11a,gam11a,spm11alo,spm11ahi,
     &       rm11b,gam11b,spm11blo,spm11bhi,
     &       rm2a,gam2a,spm2alo,spm2ahi,
     &       rm2b,gam2b,spm2blo,spm2bhi,
     &       rmta,expa,
     &       pt6min_in(3),
     &       pt6min,
     &       xmasq_in,xmbsq_in,
     &       xm(8),xmsq(8),
     &       p(0:3,8),xjac
      INTEGER iorder(6),id1,id2,id6
* Local variables
      REAL*8 e_sup,e_inf,x_max,x_min,app,
     &       s11,s2,s1,q11(0:3),q2(0:3),q1(0:3),qaux(0:3),
     &       e11,e2,e1,p11mod,p2mod,p1mod,pmod,
     &       saux,eaux,t_inf,t_sup,t1,pauxmod,p0mod,
     &       pi,twopi,fourpi,
     &       xmasq,xmbsq,rmtasq,costh,sinth,ph,ran2,
     &       signfactor,smallmass
     
      INTEGER idum,io1,io2,io3,io4,io5,io6,iaux,i
     
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi,smallmass=1.d-3)
      
* Mutimapping
      REAL*8 alfa1,alfa2,auxalfa2,
     &       alfa3,auxalfa3,alfa4,auxalfa4,
     &       xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &       xw2_min,xw2_max,xf3_min,xf3_max,
     &       vegas_x
* Mutimapping t
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
        print*,'***phsp1_2_3 ERROR'
        STOP
      ENDIF


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

* saux = (q1+q2)**2 = qaux**2
      e_sup = sqrt(shat)-xm(io6)
      e_inf = xm(io1)+xm(io2)+xm(io3)+xm(io4)+xm(io5)
      IF(e_sup.LE.e_inf)THEN    
        xjac=0.d0
        RETURN
      ENDIF
c 5-region multimapping for saux
c      eaux=(e_sup-e_inf)*x(1)+e_inf
c      saux=eaux*eaux
c      xjac=xjac*(e_sup-e_inf)*2.d0*eaux
      CALL PHSP_INI5(e_sup,e_inf,
     &               spm0alo,spm0ahi,rm0a,gam0a,
     &               spm0blo,spm0bhi,rm0b,gam0b,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max)
      IF(x(1).LE.alfa1 .AND. alfa1.GT.0.d0)THEN
        vegas_x=x(1)/alfa1
        eaux=(xf1_max-xf1_min)*vegas_x+xf1_min
        saux=eaux*eaux
        xjac=xjac*(xf1_max-xf1_min)*2.d0*eaux/alfa1
      ELSEIF(x(1).LE.alfa2 .AND. alfa2.GT.0.d0)THEN
        vegas_x=(x(1)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        saux=(gam0a*rm0a*tan((x_max-x_min)*vegas_x+x_min)+rm0a**2)
        IF(saux.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((saux-rm0a**2)**2+
     &       (gam0a*rm0a)**2)/(gam0a*rm0a)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(1)=1
        imapregion(1)=2
c end giuseppe 29/09/2007
      ELSEIF(x(1).LE.alfa3 .AND. alfa3.GT.0.d0)THEN
        vegas_x=(x(1)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        eaux=(xf2_max-xf2_min)*vegas_x+xf2_min
        saux=eaux*eaux
        xjac=xjac*(xf2_max-xf2_min)*2.d0*eaux
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(x(1).LE.alfa4 .AND. alfa4.GT.0.d0)THEN
        vegas_x=(x(1)-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/auxalfa4
        x_max=xw2_max
        x_min=xw2_min
        saux=(gam0b*rm0b*tan((x_max-x_min)*vegas_x+x_min)+rm0b**2)
        IF(saux.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((saux-rm0b**2)**2+
     &       (gam0b*rm0b)**2)/(gam0b*rm0b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
c giuseppe 29/09/2007
        isresonant(1)=1
        imapregion(1)=4
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(1)-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/(1.d0-auxalfa4)
        eaux=(xf3_max-xf3_min)*vegas_x+xf3_min
        saux=eaux*eaux
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
     &        +0.5d0/shat*(shat+xmbsq-xmasq)*(shat+saux-xmsq(io6)))
     &       /p0mod/pmod  
      costh=costh*signfactor 

* qaux and p6
      costh=min(costh,1.d0)
      costh=max(costh,-1.d0)
      sinth=sqrt(1.d0-costh*costh)
      ph=ran2(idum)*twopi
      p(0,io6)=sqrt(pmod*pmod+xmsq(io6))
      p(1,io6)=pmod*sinth*cos(ph)
      p(2,io6)=pmod*sinth*sin(ph)
      p(3,io6)=pmod*costh
      qaux(0)=sqrt(pmod*pmod+saux)
      qaux(1)=-p(1,io6)
      qaux(2)=-p(2,io6)
      qaux(3)=-p(3,io6)
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
      IF(x(3).LE.alfa1 .AND. alfa1.GT.0.d0)THEN
        vegas_x=x(3)/alfa1
        e1=(xf1_max-xf1_min)*vegas_x+xf1_min
        s1=e1*e1
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e1/alfa1
      ELSEIF(x(3).LE.alfa2 .AND. alfa2.GT.0.d0)THEN
        vegas_x=(x(3)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        s1=(gam1a*rm1a*tan((x_max-x_min)*vegas_x+x_min)+rm1a**2)
        IF(s1.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s1-rm1a**2)**2+
     &       (gam1a*rm1a)**2)/(gam1a*rm1a)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(2)=1
        imapregion(2)=2
c end giuseppe 29/09/2007
      ELSEIF(x(3).LE.alfa3 .AND. alfa3.GT.0.d0)THEN
        vegas_x=(x(3)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        e1=(xf2_max-xf2_min)*vegas_x+xf2_min
        s1=e1*e1
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e1
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(x(3).LE.alfa4 .AND. alfa4.GT.0.d0)THEN
        vegas_x=(x(3)-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/auxalfa4
        x_max=xw2_max
        x_min=xw2_min
        s1=(gam1b*rm1b*tan((x_max-x_min)*vegas_x+x_min)+rm1b**2)
        IF(s1.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s1-rm1b**2)**2+
     &       (gam1b*rm1b)**2)/(gam1b*rm1b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
c giuseppe 29/09/2007
        isresonant(2)=1
        imapregion(2)=4
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(3)-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/(1.d0-auxalfa4)
        e1=(xf3_max-xf3_min)*vegas_x+xf3_min
        s1=e1*e1
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
      IF(x(4).LE.alfa1 .AND. alfa1.GT.0.d0)THEN
        vegas_x=x(4)/alfa1
        e11=(xf1_max-xf1_min)*vegas_x+xf1_min
        s11=e11*e11
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e11/alfa1
      ELSEIF(x(4).LE.alfa2 .AND. alfa2.GT.0.d0)THEN
        vegas_x=(x(4)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        s11=(gam11a*rm11a*tan((x_max-x_min)*vegas_x+x_min)+rm11a**2)
        IF(s11.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s11-rm11a**2)**2+
     &       (gam11a*rm11a)**2)/(gam11a*rm11a)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(4)=1
        imapregion(4)=2
c end giuseppe 29/09/2007
      ELSEIF(x(4).LE.alfa3 .AND. alfa3.GT.0.d0)THEN
        vegas_x=(x(4)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        e11=(xf2_max-xf2_min)*vegas_x+xf2_min
        s11=e11*e11
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e11
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(x(4).LE.alfa4 .AND. alfa4.GT.0.d0)THEN
        vegas_x=(x(4)-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/auxalfa4
        x_max=xw2_max
        x_min=xw2_min
        s11=(gam11b*rm11b*tan((x_max-x_min)*vegas_x+x_min)+rm11b**2)
        IF(s11.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s11-rm11b**2)**2+
     &       (gam11b*rm11b)**2)/(gam11b*rm11b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
c giuseppe 29/09/2007
        isresonant(4)=1
        imapregion(4)=4
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(4)-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/(1.d0-auxalfa4)
        e11=(xf3_max-xf3_min)*vegas_x+xf3_min
        s11=e11*e11
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
      IF(x(5).LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=x(5)/alfa1
        e2=(xf1_max-xf1_min)*vegas_x+xf1_min
        s2=e2*e2
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e2/alfa1
      ELSEIF(x(5).LE.alfa2.AND.alfa2.GT.0.d0)THEN
        vegas_x=(x(5)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        s2=(gam2a*rm2a*tan((x_max-x_min)*vegas_x+x_min)+rm2a**2)
        IF(s2.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s2-rm2a**2)**2+
     &      (gam2a*rm2a)**2)/(gam2a*rm2a)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(3)=1
        imapregion(3)=2
c end giuseppe 29/09/2007
      ELSEIF(x(5).LE.alfa3.AND.alfa3.GT.0.d0)THEN
        vegas_x=(x(5)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        e2=(xf2_max-xf2_min)*vegas_x+xf2_min
        s2=e2*e2
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e2
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(x(5).LE.alfa4.AND.alfa4.GT.0.d0)THEN
        vegas_x=(x(5)-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/auxalfa4
        x_max=xw2_max
        x_min=xw2_min
        s2=(gam2b*rm2b*tan((x_max-x_min)*vegas_x+x_min)+rm2b**2)
        IF(s2.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s2-rm2b**2)**2+
     &      (gam2b*rm2b)**2)/(gam2b*rm2b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
c giuseppe 29/09/2007
        isresonant(3)=1
        imapregion(3)=4
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(5)-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
        e2=(xf3_max-xf3_min)*vegas_x+xf3_min
        s2=e2*e2
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
      CALL TWOBODY(qaux,saux,pauxmod,s1,s2,x(6),x(7),q1,q2)
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
      CALL TWOBODY(q1,s1,p1mod,s11,xmsq(io3),x(8),x(9),q11,p(0,io3))
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

* q2 --> p4 and p5
      app=((s2-xmsq(io4)-xmsq(io5))**2
     &              -4.d0*xmsq(io4)*xmsq(io5))/(4.d0*s2)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p2mod=sqrt(app)
      CALL TWOBODY(q2,s2,p2mod,xmsq(io4),xmsq(io5),x(12),x(13),
     &             p(0,io4),p(0,io5))
* beta/2/(4*pi)**2
      xjac=xjac*p2mod/sqrt(s2)/(fourpi)**2
      xjac=xjac*fourpi

      RETURN
      END
      

