************************************************************************
* Subroutine to compute the jacobian corresponding to the phase space 
* mapping:
*         [3(rm1[a,b]) --> 2(rm11[a,b],p1,p2)+p5] 
*               + [3(rm2[a,b]) --> 2(rm21[a,b],p3,p4)+p6]
* where rm1[a,b],rm11[a,b],rm2[a,b],rm21[a,b] are the masses of the 
* intermediate resonances
* rm1,rm11,rm2 and rm21 can resonate at two different masses (e.g.
* Z and H) with
* rm1a < rm1b, rm11a < rm11b, rm2a < rm2b  and rm21a < rm21b.
* If only one resonance is needed set the corresponding 
*  spmXahi=spmXalo=0.d0
*
* INPUT:
*     shat          : total energy squared
*     xm(8),xmsq(8) : masses and masses squared of the particles
*                        in the order used in the main program      
*     rmN,gamN      : N=1[a,b],2[a,b],11[a,b],21[a,b] mass, width for 
*                        intermediate resonance N.
*     spmNlo,spmNhi : N=1[a,b],2[a,b],11[a,b],21[a,b] low and high limit
*                         of resonance mapping. Outside this range the
*                         mapping is flat
*     iorder(6)     : specifies the permutation that relates internal 
*                       and  external momenta: p = p(i,iorder(1))
*     p(0:3,8)      : momenta in the hard scattering CM.
*                        Includes initial momenta.
* OUTPUT:
*     xjac          : full jacobian
************************************************************************

        subroutine phsp3_3_multi5jac(shat,
     &                    rm1a,gam1a,spm1alo,spm1ahi,
     &                    rm1b,gam1b,spm1blo,spm1bhi,
     &                    rm11a,gam11a,spm11alo,spm11ahi,
     &                    rm11b,gam11b,spm11blo,spm11bhi,
     &                    rm2a,gam2a,spm2alo,spm2ahi,
     &                    rm2b,gam2b,spm2blo,spm2bhi,
     &                    rm21a,gam21a,spm21alo,spm21ahi,
     &                    rm21b,gam21b,spm21blo,spm21bhi,
     &                    xm,xmsq,iorder,
     &                    p,xjac)
     
      IMPLICIT NONE
      REAL*8 shat,
     &       rm1a,gam1a,spm1alo,spm1ahi,
     &       rm1b,gam1b,spm1blo,spm1bhi,
     &       rm11a,gam11a,spm11alo,spm11ahi,
     &       rm11b,gam11b,spm11blo,spm11bhi,
     &       rm2a,gam2a,spm2alo,spm2ahi,
     &       rm2b,gam2b,spm2blo,spm2bhi,
     &       rm21a,gam21a,spm21alo,spm21ahi,
     &       rm21b,gam21b,spm21blo,spm21bhi,
     &       xm(8),xmsq(8),
     &       p(0:3,8),xjac
      INTEGER iorder(6)
* Local variables
      REAL*8 e_sup,e_inf,x_max,x_min,app,
     &       s11,s21,s1,s2,
     &       e11,e21,e1,e2,p11mod,p21mod,p1mod,p2mod,pmod,
     &       pi,twopi,fourpi
      INTEGER io1,io2,io3,io4,io5,io6,i
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi)

* Mutimapping
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

      s11=(p(0,io1)+p(0,io2))**2
      s21=(p(0,io3)+p(0,io4))**2
      s1=(p(0,io1)+p(0,io2)+p(0,io5))**2
      s2=(p(0,io3)+p(0,io4)+p(0,io6))**2
      DO i=1,3
        s11=s11-(p(i,io1)+p(i,io2))**2
        s21=s21-(p(i,io3)+p(i,io4))**2
        s1=s1-(p(i,io1)+p(i,io2)+p(i,io5))**2
        s2=s2-(p(i,io3)+p(i,io4)+p(i,io6))**2
      ENDDO
      IF( s11.LE.1.d-8 .OR. s21.LE.1.d-8 .OR. s1.LE.1.d-8 
     &         .OR. s2.LE.1.d-8)THEN
        xjac=0.d0
        RETURN
      ENDIF
       
* s11 = (p1+p2)**2 = q11**2
      e_sup = sqrt(shat)-xm(io3)-xm(io4)-xm(io5)-xm(io6)
      e_inf = xm(io1)+xm(io2)      
      IF( e_inf.GE.e_sup)THEN
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
     &      (gam11a*rm11a)**2)/(gam11a*rm11a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e11.LE.spm11blo.AND.alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e11
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e11.LE.spm11bhi.AND.alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s11-rm11b**2)**2+
     &      (gam11b*rm11b)**2)/(gam11b*rm11b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e11
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF      
      xjac=xjac/twopi
      
* s21 = (p3+p4)**2 = q21**2
      e_sup = sqrt(shat)-xm(io5)-xm(io6)-sqrt(s11)
      e_inf = xm(io3)+xm(io4)      
      IF( e_inf.GE.e_sup)THEN
        xjac=0.d0
        RETURN
      ENDIF     
      CALL PHSP_INI5(e_sup,e_inf,
     &               spm21alo,spm21ahi,rm21a,gam21a,
     &               spm21blo,spm21bhi,rm21b,gam21b,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max)
      e21=sqrt(s21)
      IF(e21.LE.spm21alo.AND.alfa1.GT.0.d0)THEN
        xjac=xjac*(xf1_max-xf1_min)*2.d0*e21/alfa1
      ELSEIF(e21.LE.spm21ahi.AND.alfa2.GT.0.d0)THEN
        x_max=xw1_max
        x_min=xw1_min
        xjac=xjac*(x_max-x_min)*((s21-rm21a**2)**2+
     &      (gam21a*rm21a)**2)/(gam21a*rm21a)/(1.d0-alfa1)/auxalfa2
      ELSEIF(e21.LE.spm21blo.AND.alfa3.GT.0.d0)THEN
        xjac=xjac*(xf2_max-xf2_min)*2.d0*e21
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(e21.LE.spm21bhi.AND.alfa4.GT.0.d0)THEN
        x_max=xw2_max
        x_min=xw2_min
        xjac=xjac*(x_max-x_min)*((s21-rm21b**2)**2+
     &      (gam21b*rm21b)**2)/(gam21b*rm21b)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
      ELSE
        xjac=xjac*(xf3_max-xf3_min)*2.d0*e21
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
      ENDIF      
      xjac=xjac/twopi

* s1 = (q11+p5)**2 = q1**2
      e_sup = sqrt(shat)-xm(io6)-sqrt(s21)
      e_inf = sqrt(s11)+xm(io5)      
      IF( e_inf.GE.e_sup)THEN
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
      
* s2 = (q21+p6)**2 = q2**2
      e_sup = sqrt(shat)-sqrt(s1)
      e_inf = sqrt(s21)+xm(io6)      
      IF( e_inf.GE.e_sup)THEN
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
      
* q11 and p5
      app=((s1-s11-xmsq(io5))**2-4.d0*s11*xmsq(io5))/(4.d0*s1)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p1mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p1mod/sqrt(s1)/(fourpi)**2
      xjac=xjac*fourpi

* q21 and q6
      app=((s2-s21-xmsq(io6))**2-4.d0*s21*xmsq(io6))/(4.d0*s2)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p2mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p2mod/sqrt(s2)/(fourpi)**2
      xjac=xjac*fourpi

* p1 and p2      
      app=((s11-xmsq(io1)-xmsq(io2))**2
     &       -4.d0*xmsq(io1)*xmsq(io2))/(4.d0*s11)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p11mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p11mod/sqrt(s11)/(fourpi)**2
      xjac=xjac*fourpi
      
* p3 and p4      
      app=((s21-xmsq(io3)-xmsq(io4))**2
     &        -4.d0*xmsq(io3)*xmsq(io4))/(4.d0*s21)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p21mod=sqrt(app)
* beta/2/(4*pi)**2
      xjac=xjac*p21mod/sqrt(s21)/(fourpi)**2
      xjac=xjac*fourpi
      
      RETURN
      END
      
