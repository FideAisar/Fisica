************************************************************************
* Subroutine to compute the jacobian corresponding to the phase space 
*  mapping:
* [4(rm1[a,b]) --> 3(rm11[a,b])+p4 --> 2(rm111[a,b],p1,p2)+p3+p4] +
*                  [2(rm2[a,b],p5,p6)]
* where rm1[a,b], rm11[a,b],rm111[a,b],rm2[a,b] are the masses of the 
*  intermediate resonances
* rm1, rm11 and rm111 can resonate at two different masses (e.g. Z and 
*  H) with rm1a < rm1b, rm11a < rm11b, rm111a < rm111b and rm2a < rm2b.
* If only one resonance is needed set the corresponding 
*    spmXahi=spmXalo=0.d0
* 
* INPUT:
*    shat          : total energy squared
*    xm(8),xmsq(8) : masses and masses squared of the particles
*                        in the order used in the main program      
*    rmN,gamN      : N=1[a,b],11[a,b],111[a,b],2[a,b] mass and width for
*                        intermediate resonance N. 
*    spmNlo,spmNhi : N=1[a,b],11[a,b],111[a,b],2[a,b] low and high limit
*                        of resonance mapping. Outside this range the 
*                        mapping is flat
*    iorder(6)     : specifies the permutation that relates internal and
*                        external momenta: p = p(i,iorder(1))
*    p(0:3,8)      : momenta in the hard scattering CM.
*                        Includes initial momenta.
* OUTPUT:
*    xjac          : full jacobian
************************************************************************

      subroutine phsp2_4to31_multi5jac(shat,
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
      REAL*8 shat,
     &     rm1a,gam1a,spm1alo,spm1ahi,
     &     rm1b,gam1b,spm1blo,spm1bhi, 
     &     rm11a,gam11a,spm11alo,spm11ahi,
     &     rm11b,gam11b,spm11blo,spm11bhi,
     &     rm111a,gam111a,spm111alo,spm111ahi,
     &     rm111b,gam111b,spm111blo,spm111bhi,
     &     rm2a,gam2a,spm2alo,spm2ahi,
     &     rm2b,gam2b,spm2blo,spm2bhi,
     &     xm(8),xmsq(8),
     &     p(0:3,8),xjac
      INTEGER iorder(6)
* Local variables
      REAL*8 e_sup,e_inf,x_max,x_min,app,
     &     s111,s11,s1,s2,
     &     e111,e11,e1,e2,p111mod,p11mod,p1mod,p2mod,pmod,
     &     pi,twopi,fourpi
      INTEGER io1,io2,io3,io4,io5,io6,i
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi)

* Multimapping
      REAL*8 alfa1,alfa2,auxalfa2,
     &       alfa3,auxalfa3,alfa4,auxalfa4,
     &       xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &       xw2_min,xw2_max,xf3_min,xf3_max

      io1=iorder(1)
      io2=iorder(2)
      io3=iorder(3)
      io4=iorder(4)
      io5=iorder(5)
      io6=iorder(6)

      xjac=1.d0

C All these invariant masses should be computed once and for all and
C stored in some array.

      s111=(p(0,io1)+p(0,io2))**2
      s11=(p(0,io1)+p(0,io2)+p(0,io3))**2
      s1=(p(0,io1)+p(0,io2)+p(0,io3)+p(0,io4))**2
      s2=(p(0,io5)+p(0,io6))**2
      DO i=1,3
        s111=s111-(p(i,io1)+p(i,io2))**2
        s11=s11-(p(i,io1)+p(i,io2)+p(i,io3))**2
        s1=s1-(p(i,io1)+p(i,io2)+p(i,io3)+p(i,io4))**2
        s2=s2-(p(i,io5)+p(i,io6))**2
      ENDDO
      IF( s111.LE.1.d-8 .OR. s11.LE.1.d-8 .OR. s1.LE.1.d-8
     &          .OR. s2.LE.1.d-8 )THEN
	xjac=0.d0
        RETURN
      ENDIF

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
      e111=sqrt(s111)
      IF(e111.LE.spm111alo.AND.alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e111/alfa1
      ELSEIF(e111.LE.spm111ahi.AND.alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s111-rm111a**2)**2+
     &       (gam111a*rm111a)**2)/(gam111a*rm111a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e111.LE.spm111blo.AND.alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e111
     &                   /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e111.LE.spm111bhi.AND.alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s111-rm111b**2)**2+
     &      (gam111b*rm111b)**2)/(gam111b*rm111b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e111
     &                   /(1.d0-alfa1)/(1.d0-auxalfa2)
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
      e11=sqrt(s11)
      IF(e11.LE.spm11alo.AND.alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e11/alfa1
      ELSEIF(e11.LE.spm11ahi.AND.alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s11-rm11a**2)**2+
     &       (gam11a*rm11a)**2)/(gam11a*rm11a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e11.LE.spm11blo.AND.alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e11
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e11.LE.spm11bhi.AND.alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s11-rm11b**2)**2+
     &       (gam11b*rm11b)**2)/(gam11b*rm11b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
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
      e1=sqrt(s1)
      IF(e1.LE.spm1alo.AND.alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e1/alfa1
      ELSEIF(e1.LE.spm1ahi.AND.alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s1-rm1a**2)**2+
     &       (gam1a*rm1a)**2)/(gam1a*rm1a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e1.LE.spm1blo.AND.alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e1
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e1.LE.spm1bhi.AND.alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s1-rm1b**2)**2+
     &      (gam1b*rm1b)**2)/(gam1b*rm1b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
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
      e2=sqrt(s2)
      IF(e2.LE.spm2alo .AND. alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e2/alfa1
      ELSEIF(e2.LE.spm2ahi .AND. alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s2-rm2a**2)**2+
     &       (gam2a*rm2a)**2)/(gam2a*rm2a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e2.LE.spm2blo .AND. alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e2
     &       /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e2.LE.spm2bhi .AND. alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s2-rm2b**2)**2+
     &       (gam2b*rm2b)**2)/(gam2b*rm2b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
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
* beta/2/(4*pi)**2
      xjac=xjac*p11mod/sqrt(s11)/(fourpi)**2
      xjac=xjac*fourpi

* q111 -> p1 and p2
      app=((s111-xmsq(io1)-xmsq(io2))**2-4.d0*xmsq(io1)*xmsq(io2))
     &     /(4.d0*s111)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p111mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p111mod/sqrt(s111)/(fourpi)**2
      xjac=xjac*fourpi

* q2 --> p5 and p6
      app=((s2-xmsq(io5)-xmsq(io6))**2-4.d0*xmsq(io5)*xmsq(io6))
     &     /(4.d0*s2)
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

