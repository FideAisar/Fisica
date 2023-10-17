c
c Last update: Mar 29, 2014
c
************************************************************************
* Subroutine to compute the momenta p1,....,p6 of six particles of mass
* xm1,...,xm6 given a point in the 14 dimensional hypercube [0,1]^14
* Returns also xjac, the weight due to the jacobian.
*
* The phase space mapping corresponds to
*   [4(rm1[a,b]) --> 2(rm11[a,b,c],p1,p2)+2(rm12[a,b,c],p3,p4)] +
*                     [2(rm2[a,b,c],p5,p6)]
* where rm1[a,b], rm11[a,b,c],rm2[a,b,c],rm12[a,b,c] are the masses of the
*   intermediate resonances
* rm1 can resonate at two different masses
*    (e.g. H, H') with rm1a < rm1b.
* rm11, rm12 and rm2 can resonate at three different masses
*    (e.g. Z, H, H') with
* rmXa < rmXb < rmXc.
* CHECK
* If only one resonance is needed set the corresponding
*    spmXahi=spmXalo=0.d0
* ENDCHECK
*
* If in the main x had dimension 15 with x(1) and x(2) fixing the
*  initial state particle momentum fraction the subroutine can be called
*   as     call phsp2_4(x(3),shat,.......)
*

* INPUT:

*  x,shat        : vegas random points, total energy squared
*  xm(8),xmsq(8) : masses and masses squared of the particles
*                        in the order used in the main program
*  rmN,gamN      : N=1[a,b],2[a,b,c],11[a,b,c],12[a,b,c] mass and width for
*                        intermediate resonance N.
*  spmNlo,spmNhi : N=1[a,b],11[a,b,c],12[a,b,c],1,2[a,b,c] low and high limit of
*                        resonance mapping. Outside this range the
*                         mapping is flat
*  iorder(6)     : specifies the permutation that relates internal and
*
*                      external momenta: p = p(i,iorder(1))
* OUTPUT:

*  p(0:3,8)      : momenta in the hard scattering CM.
*                        Includes initial momenta.
*  xjac          : full jacobian
************************************************************************

        subroutine phsp2_4_multi7(x,shat,
     &                    rm1a,gam1a,spm1alo,spm1ahi,
     &                    rm1b,gam1b,spm1blo,spm1bhi,
     &                    rm11a,gam11a,spm11alo,spm11ahi,
     &                    rm11b,gam11b,spm11blo,spm11bhi,
     &                    rm11c,gam11c,spm11clo,spm11chi,
     &                    rm12a,gam12a,spm12alo,spm12ahi,
     &                    rm12b,gam12b,spm12blo,spm12bhi,
     &                    rm12c,gam12c,spm12clo,spm12chi,
     &                    rm2a,gam2a,spm2alo,spm2ahi,
     &                    rm2b,gam2b,spm2blo,spm2bhi,
     &                    rm2c,gam2c,spm2clo,spm2chi,
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
     &       rm2a,gam2a,spm2alo,spm2ahi,
     &       rm2b,gam2b,spm2blo,spm2bhi,
     &       rm2c,gam2c,spm2clo,spm2chi,
     &       xm(8),xmsq(8),
     &       p(0:3,8),xjac,xjacvar,svar
      INTEGER iorder(6)
* Local variables
      REAL*8 e_sup,e_inf,x_max,x_min,app,
     &       s11,s12,s1,s2,q11(0:3),q12(0:3),q1(0:3),q2(0:3),
     &       e11,e12,e1,e2,p11mod,p12mod,p1mod,p2mod,pmod,
     &       pi,twopi,fourpi,
     &       costh,sinth,ph,ran2
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

* s11 = (p1+p2)**2 = q11**2
      e_sup = sqrt(shat)-xm(io3)-xm(io4)-xm(io5)-xm(io6)
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
      CALL PHSP_7(x(1),
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
      e_sup = sqrt(shat)-xm(io5)-xm(io6)-sqrt(s11)
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
      CALL PHSP_7(x(2),
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

* s1 = (q11+q12)**2 = q1**2
      e_sup = sqrt(shat)-xm(io5)-xm(io6)
      e_inf = sqrt(s11)+ sqrt(s12)
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

* s2 = (p5+p6)**2 = q2**2
      e_sup = sqrt(shat)-sqrt(s1)
      e_inf = xm(io5)+xm(io6)
      IF( e_inf.GE.e_sup)THEN
        xjac=0.d0
        RETURN
      ENDIF
      iregion=2
      CALL PHSP_INI7(e_sup,e_inf,
     &               spm2alo,spm2ahi,rm2a,gam2a,
     &               spm2blo,spm2bhi,rm2b,gam2b,
     &               spm2clo,spm2chi,rm2c,gam2c,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               alfa5,auxalfa5,alfa6,auxalfa6,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &               xw3_min,xw3_max,xf4_min,xf4_max)
      CALL PHSP_7(x(4),
     &               rm2a,gam2a,rm2b,gam2b,rm2c,gam2c,
     &               alfa1,alfa2,auxalfa2,
     &               alfa3,auxalfa3,alfa4,auxalfa4,
     &               alfa5,auxalfa5,alfa6,auxalfa6,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &               xw3_min,xw3_max,xf4_min,xf4_max,
     &                     iregion,
     &                     svar,xjacvar)
      s2=svar
      xjac=xjac*xjacvar

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

* q11 and q12
      app=((s1-s11-s12)**2-4.d0*s11*s12)/(4.d0*s1)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p1mod=sqrt(app)
      CALL TWOBODY_Z(q1,s1,p1mod,s11,s12,x(6),x(7),q11,q12)
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
      CALL TWOBODY(q11,s11,p11mod,xmsq(io1),xmsq(io2),x(8),x(9),
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
      CALL TWOBODY(q12,s12,p12mod,xmsq(io3),xmsq(io4),x(10),x(11),
     &             p(0,io3),p(0,io4))
* beta/2/(4*pi)**2
      xjac=xjac*p12mod/sqrt(s12)/(fourpi)**2
      xjac=xjac*fourpi

* q2 --> p5 and p6
      app=((s2-xmsq(io5)-xmsq(io6))**2
     &             -4.d0*xmsq(io5)*xmsq(io6))/(4.d0*s2)
      IF(app.LE.0.d0)THEN
        xjac=0.d0
        RETURN
      ENDIF
      p2mod=sqrt(app)
      CALL TWOBODY(q2,s2,p2mod,xmsq(io5),xmsq(io6),x(12),x(13),
     &             p(0,io5),p(0,io6))
* beta/2/(4*pi)**2
      xjac=xjac*p2mod/sqrt(s2)/(fourpi)**2
      xjac=xjac*fourpi

      RETURN
      END


