c
c Last update: Sep 29, 2007
c
************************************************************************
* Subroutine to compute the momenta p1,....,p6 of six particles of mass 
* xm1,...,xm6 given a point in the 14 dimensional hypercube [0,1]^14
* Returns also xjac, the weight due to the jacobian.
*
* The phase space mapping corresponds to 
*  [4(rm1[a,b]) --> 3(rm11[a,b])+p4 --> 2(rm111[a,b],p1,p2)+p3+p4] +
*                   [2(rm2[a,b],p5,p6)]
* where rm1[a,b], rm11[a,b],rm111[a,b],rm2[a,b] are the masses of the 
*  intermediate resonances
* rm1, rm11 and rm111 can resonate at two different masses (e.g. Z and 
*  H) with rm1a < rm1b, rm11a < rm11b, rm111a < rm111b and rm2a < rm2b.
* If only one resonance is needed set the corresponding 
*    spmXahi=spmXalo=0.d0
* 
* If in the main x had dimension 15 with x(1) and x(2) fixing the 
*  initial state particle momentum fraction the subroutine can be called
*   as     call phsp2_4(x(3),shat,.......)
*
* INPUT:
*  x,shat        : vegas random points, total energy squared
*  xm(8),xmsq(8) : masses and masses squared of the particles
*                        in the order used in the main program      
*  rmN,gamN      : N=1[a,b],11[a,b],111[a,b],2[a,b] mass and width for 
*                        intermediate resonance N. 
*  spmNlo,spmNhi : N=1[a,b],11[a,b],111[a,b],2[a,b] low and high limit 
*                        of resonance mapping. Outside this range the 
*                        mapping is flat
*  iorder(6)     : specifies the permutation that relates internal and
*                      external momenta: p = p(i,iorder(1))
* OUTPUT:
*  p(0:3,8)      : momenta in the hard scattering CM.
*                        Includes initial momenta.
*  xjac          : full jacobian
************************************************************************

      subroutine phsp2_4to31_multi5(x,shat,
     &                    rm1a,gam1a,spm1alo,spm1ahi,
     &                    rm1b,gam1b,spm1blo,spm1bhi, 
     &                    rm11a,gam11a,spm11alo,spm11ahi,
     &                    rm11b,gam11b,spm11blo,spm11bhi,
     &                    rm111a,gam111a,spm111alo,spm111ahi,
     &                    rm111b,gam111b,spm111blo,spm111bhi,
     &                    rm2a,gam2a,spm2alo,spm2ahi,
     &                    rm2b,gam2b,spm2blo,spm2bhi,
     &                    xm,xmsq,iorder,
     &                    p,xjac)                

      IMPLICIT NONE
      REAL*8 x(13),shat,
     &       rm1a,gam1a,spm1alo,spm1ahi,
     &       rm1b,gam1b,spm1blo,spm1bhi, 
     &       rm11a,gam11a,spm11alo,spm11ahi,
     &       rm11b,gam11b,spm11blo,spm11bhi,
     &       rm111a,gam111a,spm111alo,spm111ahi,
     &       rm111b,gam111b,spm111blo,spm111bhi,
     &       rm2a,gam2a,spm2alo,spm2ahi,
     &       rm2b,gam2b,spm2blo,spm2bhi,
     &       xm(8),xmsq(8),
     &       p(0:3,8),xjac
      INTEGER iorder(6)
* Local variables
      REAL*8 e_sup,e_inf,x_max,x_min,app,
     &       s111,s11,s1,s2,q111(0:3),q11(0:3),q1(0:3),q2(0:3),
     &       e111,e11,e1,e2,p111mod,p11mod,p1mod,p2mod,pmod,
     &       pi,twopi,fourpi,
     &       costh,sinth,ph,ran2
      INTEGER idum,io1,io2,io3,io4,io5,io6,i
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi)

* Multimapping
      REAL*8 alfa1,alfa2,auxalfa2,
     &       alfa3,auxalfa3,alfa4,auxalfa4,
     &       xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &       xw2_min,xw2_max,xf3_min,xf3_max,
     &       vegas_x  

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


* s111 = (p1+p2)**2 = q111**2
      e_sup = sqrt(shat)-xm(io3)-xm(io4)-xm(io5)-xm(io6)
      e_inf = xm(io1)+xm(io2)
      IF(e_inf.GE.e_sup)THEN
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
      IF(x(1).LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=x(1)/alfa1
        e111=(xf1_max-xf1_min)*vegas_x+xf1_min
        s111=e111*e111
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e111/alfa1
      ELSEIF(x(1).LE.alfa2.AND.alfa2.GT.0.d0)THEN
        vegas_x=(x(1)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        s111=(gam111a*rm111a*tan((x_max-x_min)*vegas_x+x_min)+rm111a**2)
        IF(s111.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s111-rm111a**2)**2+
     &       (gam111a*rm111a)**2)/(gam111a*rm111a)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(4)=1
        imapregion(4)=2
c end giuseppe 29/09/2007
      ELSEIF(x(1).LE.alfa3.AND.alfa3.GT.0.d0)THEN
        vegas_x=(x(1)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        e111=(xf2_max-xf2_min)*vegas_x+xf2_min
        s111=e111*e111
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e111
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(x(1).LE.alfa4.AND.alfa4.GT.0.d0)THEN
        vegas_x=(x(1)-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/auxalfa4
        x_max=xw2_max
        x_min=xw2_min
        s111=(gam111b*rm111b*tan((x_max-x_min)*vegas_x+x_min)+rm111b**2)
        IF(s111.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s111-rm111b**2)**2+
     &       (gam111b*rm111b)**2)/(gam111b*rm111b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
c giuseppe 29/09/2007
        isresonant(4)=1
        imapregion(4)=4
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(1)-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
        e111=(xf3_max-xf3_min)*vegas_x+xf3_min
        s111=e111*e111
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e111
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF
      xjac=xjac/twopi

* s11 = (q111+p3)**2 = q11**2
      e_sup = sqrt(shat)-xm(io4)-xm(io5)-xm(io6)
      e_inf = sqrt(s111)+xm(io3)
      IF(e_inf.GE.e_sup)THEN
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
      IF(x(2).LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=x(2)/alfa1
        e11=(xf1_max-xf1_min)*vegas_x+xf1_min
        s11=e11*e11
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e11/alfa1
      ELSEIF(x(2).LE.alfa2.AND.alfa2.GT.0.d0)THEN
        vegas_x=(x(2)-alfa1)/(1.d0-alfa1)/auxalfa2
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
        isresonant(3)=1
        imapregion(3)=2
c end giuseppe 29/09/2007
      ELSEIF(x(2).LE.alfa3.AND.alfa3.GT.0.d0)THEN
        vegas_x=(x(2)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        e11=(xf2_max-xf2_min)*vegas_x+xf2_min
        s11=e11*e11
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e11
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(x(2).LE.alfa4.AND.alfa4.GT.0.d0)THEN
        vegas_x=(x(2)-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/auxalfa4
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
        isresonant(3)=1
        imapregion(3)=4
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(2)-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
        e11=(xf3_max-xf3_min)*vegas_x+xf3_min
        s11=e11*e11
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e11
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF
      xjac=xjac/twopi

* s1 = (q11+p4)**2 = q1**2
      e_sup = sqrt(shat)-xm(io5)-xm(io6)
      e_inf = sqrt(s11)+xm(io4)
      IF(e_inf.GE.e_sup)THEN
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
      IF(x(3).LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=x(3)/alfa1
        e1=(xf1_max-xf1_min)*vegas_x+xf1_min
        s1=e1*e1
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e1/alfa1
      ELSEIF(x(3).LE.alfa2.AND.alfa2.GT.0.d0)THEN
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
        isresonant(1)=1
        imapregion(1)=2
c end giuseppe 29/09/2007
      ELSEIF(x(3).LE.alfa3.AND.alfa3.GT.0.d0)THEN
        vegas_x=(x(3)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        e1=(xf2_max-xf2_min)*vegas_x+xf2_min
        s1=e1*e1
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e1
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(x(3).LE.alfa4.AND.alfa4.GT.0.d0)THEN
        vegas_x=(x(3)-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/auxalfa4
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
        isresonant(1)=1
        imapregion(1)=4
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(3)-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
        e1=(xf3_max-xf3_min)*vegas_x+xf3_min
        s1=e1*e1
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e1
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF
      xjac=xjac/twopi

* s2 = (p5+p6)**2 = q2**2
      e_sup = sqrt(shat)-sqrt(s1)
      e_inf = xm(io5)+xm(io6)
      IF(e_inf.GE.e_sup)THEN
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
      IF(x(4).LE.alfa1 .AND. alfa1.GT.0.d0)THEN
        vegas_x=x(4)/alfa1
        e2=(xf1_max-xf1_min)*vegas_x+xf1_min
        s2=e2*e2
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e2/alfa1
      ELSEIF(x(4).LE.alfa2 .AND. alfa2.GT.0.d0)THEN 
        vegas_x=(x(4)-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        s2=(gam2a*rm2a*tan((x_max-x_min)*vegas_x+x_min)+rm2a**2)
        IF(s2.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s2-rm2a**2)**2+
     &       (gam2a*rm2a)**2)/(gam2a*rm2a)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(2)=1
        imapregion(2)=2
c end giuseppe 29/09/2007
      ELSEIF(x(4).LE.alfa3 .AND. alfa3.GT.0.d0)THEN
        vegas_x=(x(4)-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        e2=(xf2_max-xf2_min)*vegas_x+xf2_min
        s2=e2*e2
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e2
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(x(4).LE.alfa4 .AND. alfa4.GT.0.d0)THEN
        vegas_x=(x(4)-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/auxalfa4
        x_max=xw2_max
        x_min=xw2_min
        s2=(gam2b*rm2b*tan((x_max-x_min)*vegas_x+x_min)+rm2b**2)
        IF(s2.LE.0.d0)THEN
          xjac=0.d0
          RETURN
        ENDIF
        xjac=xjac*(x_max-x_min)*((s2-rm2b**2)**2+
     &       (gam2b*rm2b)**2)/(gam2b*rm2b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
c giuseppe 29/09/2007
        isresonant(2)=1
        imapregion(2)=4
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(x(4)-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/(1.d0-auxalfa4)
        e2=(xf3_max-xf3_min)*vegas_x+xf3_min
        s2=e2*e2
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e2
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)
     &       /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF
      xjac=xjac/twopi


* q1 and q2
      app=((shat-s1-s2)**2-4.d0*s1*s2)/(4.d0*shat)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      pmod=sqrt(app)
      costh=min(x(5)*2.d0-1.d0,1.d0)
      costh=max(costh,-1.d0)
      sinth=sqrt(1.d0-costh*costh)
      ph=ran2(idum)*twopi
      q1(0)=sqrt(pmod*pmod+s1)
      q1(1)=pmod*sinth*cos(ph)
      q1(2)=pmod*sinth*sin(ph)
      q1(3)=pmod*costh
      q2(0)=sqrt(pmod*pmod+s2)
      q2(1)=-q1(1)
      q2(2)=-q1(2)
      q2(3)=-q1(3)
* beta/2/(4*pi)**2
      xjac=xjac*pmod/sqrt(shat)/(fourpi)**2
      xjac=xjac*fourpi

* q11 and p4
      app=((s1-s11-xmsq(io4))**2-4.d0*s11*xmsq(io4))/(4.d0*s1)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p1mod=sqrt(app)
      CALL TWOBODY(q1,s1,p1mod,s11,xmsq(io4),x(6),x(7),q11,p(0,io4))
* beta/2/(4*pi)**2
      xjac=xjac*p1mod/sqrt(s1)/(fourpi)**2
      xjac=xjac*fourpi

* q111 and p3
      app=((s11-s111-xmsq(io3))**2-4.d0*s111*xmsq(io3))/(4.d0*s11)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p11mod=sqrt(app)
      CALL TWOBODY(q11,s11,p11mod,s111,xmsq(io3),x(8),x(9),q111,
     &     p(0,io3))
* beta/2/(4*pi)**2
      xjac=xjac*p11mod/sqrt(s11)/(fourpi)**2
      xjac=xjac*fourpi

* q111 --> p1 and p2
      app=((s111-xmsq(io1)-xmsq(io2))**2-4.d0*xmsq(io1)*xmsq(io2))
     &     /(4.d0*s111)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p111mod=sqrt(app)
      CALL TWOBODY(q111,s111,p111mod,xmsq(io1),xmsq(io2),x(10),x(11),
     &     p(0,io1),p(0,io2))
* beta/2/(4*pi)**2
      xjac=xjac*p111mod/sqrt(s111)/(fourpi)**2
      xjac=xjac*fourpi

*q2 --> p5 and p6
      app=((s2-xmsq(io5)-xmsq(io6))**2-4.d0*xmsq(io5)*xmsq(io6))
     &     /(4.d0*s2)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p2mod=sqrt(app)
      CALL TWOBODY(q2,s2,p2mod,xmsq(io5),xmsq(io6),x(12),x(13),
     &     p(0,io5),p(0,io6))
* beta/2/(4*pi)**2
      xjac=xjac*p2mod/sqrt(s2)/(fourpi)**2
      xjac=xjac*fourpi


      RETURN
      END
