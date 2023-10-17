*****************************************************************************
* Auxiliary subroutine. Defines extrema for flat+Breit-Wigner+flat
*                                              +Breit-Wigner+flat
*                                              +Breit-Wigner+flat mapping
* of resonances. Checks for position of estrsup,estrinf with respect to
* spm1lo,spm1hi
* INPUT:
*   estrsup,estrinf: upper and lower limit for the generated mass
*   spm1lo,spm1hi  : lower and upper limit of first Breit-Wigner region
*   rmx1,gamx1       : mass and width of first resonance
*   spm2lo,spm2hi  : lower and upper limit of second Breit-Wigner region
*   rmx2,gamx2       : mass and width of second resonance
*   spm3lo,spm3hi  : lower and upper limit of second Breit-Wigner region
*   rmx3,gamx3       : mass and width of second resonance
* OUTPUT:
*  alfa1,alfa2,alfa3,alfa4,alfa5,alfa6,: separation points in [0,1]
*                                        between the seven regions
*                      Can also be 0 or 1. Some care necessary in this cases
*  auxalfa2,xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,xw2_min,xw2_max,
*  xf3_min,xf3_max,auxalfa3,auxalfa4,xw3_min,xw3_max,
*  xf4_min,xf4_max,auxalfa5,auxalfa6
*                    : auxiliary quantities used for generating the mass
*                      and computing the jacobian
*****************************************************************************
      SUBROUTINE PHSP_INI7(estrsup,estrinf,
     &                     spm1lo,spm1hi,rmx1,gamx1,
     &                     spm2lo,spm2hi,rmx2,gamx2,
     &                     spm3lo,spm3hi,rmx3,gamx3,
     &                     alfa1,alfa2,auxalfa2,
     &                     alfa3,auxalfa3,alfa4,auxalfa4,
     &                     alfa5,auxalfa5,alfa6,auxalfa6,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &               xw3_min,xw3_max,xf4_min,xf4_max
     &                    )

      IMPLICIT REAL*8(a-h,o-z)

      IF(estrsup.LE.spm1lo)THEN ! only flat1
        alfa1=1.d0
        alfa2=1.d0
        alfa3=1.d0
        alfa4=1.d0
        alfa5=1.d0
        alfa6=1.d0
        auxalfa6=alfa6
        auxalfa5=alfa5
        auxalfa4=alfa4
        auxalfa3=alfa3
        auxalfa2=alfa2
        xf1_max=estrsup
        xf1_min=estrinf
        RETURN
      ELSEIF(estrsup.LE.spm1hi)THEN ! flat1 + res1 ?
        alfa2=1.d0
        alfa3=1.d0
        alfa4=1.d0
        alfa5=1.d0
        alfa6=1.d0
        auxalfa6=alfa6
        auxalfa5=alfa5
        auxalfa4=alfa4
        auxalfa3=alfa3
        auxalfa2=alfa2
        IF((estrinf.GE.spm1lo))THEN  ! res1 only
          alfa1=0.d0
          xw1_max=atan((estrsup**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((estrinf**2-rmx1**2)/(gamx1*rmx1))
        ELSE !  flat1 + res1 !
          xw1_max=atan((estrsup**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((spm1lo**2-rmx1**2)/(gamx1*rmx1))
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
        alfa5=1.d0
        alfa6=1.d0
        auxalfa6=alfa6
        auxalfa5=alfa5
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
             ! estrsup.GT.spm2hi: flat1 + res1 + flat2 + res2?
        alfa4=1.d0
        alfa5=1.d0
        alfa6=1.d0
        auxalfa6=alfa6
        auxalfa5=alfa5
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
          xf2_min=spm1hi
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
      ELSEIF(estrsup.LE.spm3lo)THEN
                 ! estrsup.LE.spm3lo: flat1 + res1 + flat2 + res2 + flat3?
        alfa5=1.d0
        alfa6=1.d0
        auxalfa6=alfa6
        auxalfa5=alfa5
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
        RETURN
      ELSEIF(estrsup.LE.spm3hi)THEN
                 ! estrsup.LE.spm3hi: flat1 + res1 + flat2 + res2 + flat3 + res3?
        alfa6=1.d0
        auxalfa6=alfa6
        IF(estrinf.GE.spm3lo)THEN  ! only res3
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          alfa4=0.d0
          alfa5=0.d0
          auxalfa5=alfa5
          auxalfa4=alfa4
          auxalfa3=alfa3
          auxalfa2=alfa2
          xw3_max=atan((estrsup**2-rmx3**2)/(gamx3*rmx3))
          xw3_min=atan((estrinf**2-rmx3**2)/(gamx3*rmx3))
        ELSEIF(estrinf.GE.spm2hi)THEN  ! flat3 + res3
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          alfa4=0.d0
          auxalfa4=alfa4
          auxalfa3=alfa3
          auxalfa2=alfa2
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((estrsup**2-rmx3**2)/(gamx3*rmx3))

          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
* If rkappa6low=0 then estrsup is numerically indistinguishable from spm3lo
          IF(rkappa6low.LE.0.d0)THEN
            alfa5=1.d0
            auxalfa5=alfa5
            xf3_max=estrsup
            xf3_min=estrinf
            RETURN
          ENDIF
          xf3_max=spm3lo
          xf3_min=estrinf
          rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
          rk5=rkappa5high/rkappa6low
          alfa5=rk5/(rk5+1.d0)
          auxalfa5=alfa5
! end flat3 + res3
        ELSEIF(estrinf.GE.spm2lo)THEN  ! res2 + flat3 + res3
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          auxalfa3=alfa3
          auxalfa2=alfa2
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((estrsup**2-rmx3**2)/(gamx3*rmx3))
          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          xf3_min=spm2hi    !! ERRORE IN FILE ORIGINALE
* If rkappa6low=0 then estrsup is numerically indistinguishable from spm3lo
C          IF(rkappa6low.LE.0.d0)THEN
C            alfa5=1.d0
C            auxalfa5=alfa5
C            xf3_max=estrsup
C          ELSE
            xf3_max=spm3lo
            rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
            rk5=rkappa5high/rkappa6low
ccc            auxalfa5=rk5/(rk5+1.d0)
            auxalfa5=0.5d0
C          ENDIF
          rkappa5low=(xf3_max-xf3_min)*2.d0*spm2hi
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          xw2_min=atan((estrinf**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5low*auxalfa5
ccc          alfa4=rk4/(rk4+1.d0)
          alfa4=1.d0/3.d0
          auxalfa4=alfa4
          alfa5=auxalfa5*(1.d0-alfa4)+alfa4
        ELSEIF(estrinf.GE.spm1hi)THEN    ! flat2 + res2 + flat3 + res3
          alfa1=0.d0
          alfa2=0.d0
          auxalfa2=alfa2
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((estrsup**2-rmx3**2)/(gamx3*rmx3))
          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          xf3_min=spm2hi
* If rkappa6low=0 then estrsup is numerically indistinguishable from spm2lo
          IF(rkappa6low.LE.0.d0)THEN
            auxalfa5=1.d0
            xf3_max=estrsup
          ELSE
            xf3_max=spm3lo
            rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
            rk5=rkappa5high/rkappa6low
ccc            auxalfa5=rk5/(rk5+1.d0)
            auxalfa5=0.5d0
          ENDIF
          rkappa5low=(xf3_max-xf3_min)*2.d0*spm2hi
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5low*auxalfa5
ccc          auxalfa4=rk4/(rk4+1.d0)
          auxalfa4=1.d0/3.d0
          xf2_max=spm2lo
          xf2_min=estrinf
          rkappa3=(xf2_max-xf2_min)*2.d0*spm2lo
          rk3=rkappa3/rkappa4low*auxalfa4
ccc          alfa3=rk3/(rk3+1.d0)
          alfa3=1.d0/4.d0
          auxalfa3=alfa3
          alfa4=auxalfa4*(1.d0-alfa3)+alfa3
          alfa5=auxalfa5*(1.d0-auxalfa4)*(1.d0-alfa3)+alfa4
          IF(rkappa6low.LE.0.d0) alfa5=1.d0
        ELSEIF(estrinf.GE.spm1lo)THEN    ! res1 + flat2 + res2 + flat3 + res3
          alfa1=0.d0
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((estrsup**2-rmx3**2)/(gamx3*rmx3))
          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          xf3_max=spm3lo
          xf3_min=spm2hi
          rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
          rkappa5low=(xf3_max-xf3_min)*2.d0*spm2hi
          rk5=rkappa5high/rkappa6low
ccc          auxalfa5=rk5/(rk5+1.d0)
          auxalfa5=0.5d0
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5low*auxalfa5
ccc          auxalfa4=rk4/(rk4+1.d0)
          auxalfa4=1.d0/3.d0
          xf2_max=spm2lo
          xf2_min=spm1hi
          rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          rk3=rkappa3high/rkappa4low*auxalfa4
ccc          auxalfa3=rk3/(rk3+1.d0)
          auxalfa3=1.d0/4.d0
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((estrinf**2-rmx1**2)/(gamx1*rmx1))
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3low*auxalfa3
ccc          auxalfa2=rk2/(rk2+1.d0)
          auxalfa2=1.d0/5.d0
          alfa2=auxalfa2
          alfa3=auxalfa3*(1.d0-auxalfa2)+alfa2
          alfa4=auxalfa4*(1.d0-auxalfa3)*(1.d0-auxalfa2)+alfa3
          alfa5=auxalfa5*(1.d0-auxalfa4)*(1.d0-auxalfa3)
     &                          *(1.d0-auxalfa2)+alfa4
        ELSE        ! estrinf.LT.spm1lo: flat1 + res1 + flat2 + res2 + flat3 + res3
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((estrsup**2-rmx3**2)/(gamx3*rmx3))
          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          xf3_max=spm3lo
          xf3_min=spm2hi
          rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
          rkappa5low=(xf3_max-xf3_min)*2.d0*spm2hi
          rk5=rkappa5high/rkappa6low
ccc          auxalfa5=rk5/(rk5+1.d0)
          auxalfa5=0.5d0
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5low*auxalfa5
ccc          auxalfa4=rk4/(rk4+1.d0)
          auxalfa4=1.d0/3.d0
          xf2_max=spm2lo
          xf2_min=spm1hi
          rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          rk3=rkappa3high/rkappa4low*auxalfa4
ccc          auxalfa3=rk3/(rk3+1.d0)
          auxalfa3=1.d0/4.d0
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
          auxalfa2=1.d0/5.d0
          xf1_max=spm1lo
          xf1_min=estrinf
          rkappa1=(xf1_max-xf1_min)*2.d0*spm1lo
          rk1=rkappa1/rkappa2low*auxalfa2
ccc          alfa1=rk1/(rk1+1.d0)
          alfa1=1.d0/6.d0
          alfa2=auxalfa2*(1.d0-alfa1)+alfa1
          alfa3=auxalfa3*(1.d0-auxalfa2)*(1.d0-alfa1)+alfa2
          alfa4=auxalfa4*(1.d0-auxalfa3)*(1.d0-auxalfa2)*(1.d0-alfa1)
     &                      +alfa3
          alfa5=auxalfa5*(1.d0-auxalfa4)*(1.d0-auxalfa3)*
     &               (1.d0-auxalfa2)*(1.d0-alfa1)+alfa4
        ENDIF
      ELSE
                 ! estrsup.GE.spm3hi: flat1 + res1 + flat2 + res2 + flat3 +
                 ! res3 + flat4?
        IF(estrinf.GE.spm3hi)THEN  ! only flat4
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          alfa4=0.d0
          alfa5=0.d0
          alfa6=0.d0
          auxalfa6=alfa6
          auxalfa5=alfa5
          auxalfa4=alfa4
          auxalfa3=alfa3
          auxalfa2=alfa2
          xf4_max=estrsup
          xf4_min=estrinf
          RETURN
         ELSEIF(estrinf.GE.spm3lo)THEN  ! res3 + flat4
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          alfa4=0.d0
          alfa5=0.d0
          auxalfa5=alfa5
          auxalfa4=alfa4
          auxalfa3=alfa3
          auxalfa2=alfa2
          xf4_max=estrsup
          xf4_min=spm3hi
          rkappa7=(xf4_max-xf4_min)*2.d0*spm3hi
          xw3_min=atan((estrinf**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((spm3hi**2-rmx3**2)/(gamx3*rmx3))
          rkappa6high=(xw3_max-xw3_min)*
     &          ((spm3hi**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rk6=rkappa6high/rkappa7
ccc          alfa6=rk6/(rk6+1.d0)
          alfa6=0.5d0
          auxalfa6=alfa6
        ELSEIF(estrinf.GE.spm2hi)THEN  ! flat3 + res3 + flat4
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          alfa4=0.d0
          auxalfa4=alfa4
          auxalfa3=alfa3
          auxalfa2=alfa2
          xf4_max=estrsup
          xf4_min=spm3hi
          rkappa7=(xf4_max-xf4_min)*2.d0*spm3hi
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((spm3hi**2-rmx3**2)/(gamx3*rmx3))
          rkappa6high=(xw3_max-xw3_min)*
     &          ((spm3hi**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rk6=rkappa6high/rkappa7
ccc          auxalfa6=rk6/(rk6+1.d0)
          auxalfa6=0.5d0
          xf3_max=spm3lo
          xf3_min=estrinf
          rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
          rk5=rkappa5high/rkappa6low*auxalfa6
ccc          alfa5=rk5/(rk5+1.d0)
          alfa5=1.d0/3.d0
          auxalfa5=alfa5
          alfa6=auxalfa6*(1.d0-alfa5)+alfa5
        ELSEIF(estrinf.GE.spm2lo)THEN  ! res2 + flat3 + res3 + flat4
          alfa1=0.d0
          alfa2=0.d0
          alfa3=0.d0
          auxalfa3=alfa3
          auxalfa2=alfa2
          xf4_max=estrsup
          xf4_min=spm3hi
          rkappa7=(xf4_max-xf4_min)*2.d0*spm3hi
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((spm3hi**2-rmx3**2)/(gamx3*rmx3))
          rkappa6high=(xw3_max-xw3_min)*
     &          ((spm3hi**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rk6=rkappa6high/rkappa7
ccc          auxalfa6=rk6/(rk6+1.d0)
          auxalfa6=0.5d0
          xf3_max=spm3lo
          xf3_min=spm2hi
          rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
          rkappa5low=(xf3_max-xf3_min)*2.d0*spm2hi
          rk5=rkappa5high/rkappa6low*auxalfa6
ccc          auxalfa5=rk5/(rk5+1.d0)
          auxalfa5=1.d0/3.d0
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          xw2_min=atan((estrinf**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5low*auxalfa5
ccc          alfa4=rk4/(rk4+1.d0)
          alfa4=1.d0/4.d0
          auxalfa4=alfa4
          alfa5=auxalfa5*(1.d0-alfa4)+alfa4
          alfa6=auxalfa6*(1.d0-auxalfa5)*(1.d0-alfa4)+alfa5
        ELSEIF(estrinf.GE.spm1hi)THEN ! flat2 + res2 + flat3 + res3 + flat4
          alfa1=0.d0
          alfa2=0.d0
          auxalfa2=alfa2
          xf4_max=estrsup
          xf4_min=spm3hi
          rkappa7=(xf4_max-xf4_min)*2.d0*spm3hi
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((spm3hi**2-rmx3**2)/(gamx3*rmx3))
          rkappa6high=(xw3_max-xw3_min)*
     &          ((spm3hi**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rk6=rkappa6high/rkappa7
ccc          auxalfa6=rk6/(rk6+1.d0)
          auxalfa6=0.5d0
          xf3_max=spm3lo
          xf3_min=spm2hi
          rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
          rkappa5low=(xf3_max-xf3_min)*2.d0*spm2hi
          rk5=rkappa5high/rkappa6low*auxalfa6
ccc          auxalfa5=rk5/(rk5+1.d0)
          auxalfa5=1.d0/3.d0
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5low*auxalfa5
ccc          auxalfa4=rk4/(rk4+1.d0)
          auxalfa4=1.d0/4.d0
          xf2_max=spm2lo
          xf2_min=estrinf
          rkappa3=(xf2_max-xf2_min)*2.d0*spm2lo
          rk3=rkappa3/rkappa4low*auxalfa4
ccc          alfa3=rk3/(rk3+1.d0)
          auxalfa3=1.d0/5.d0
          alfa3=auxalfa3
          alfa4=auxalfa4*(1.d0-alfa3)+alfa3
          alfa5=auxalfa5*(1.d0-auxalfa4)*(1.d0-alfa3)+alfa4
          alfa6=auxalfa6*(1.d0-auxalfa5)*(1.d0-auxalfa4)*(1.d0-alfa3)
     &                               +alfa5
        ELSEIF(estrinf.GE.spm1lo)THEN ! res1 + flat2 + res2 + flat3 + res3+ flat4
          alfa1=0.d0
          xf4_max=estrsup
          xf4_min=spm3hi
          rkappa7=(xf4_max-xf4_min)*2.d0*spm3hi
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((spm3hi**2-rmx3**2)/(gamx3*rmx3))
          rkappa6high=(xw3_max-xw3_min)*
     &          ((spm3hi**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rk6=rkappa6high/rkappa7
ccc          auxalfa6=rk6/(rk6+1.d0)
          auxalfa6=0.5d0
          xf3_max=spm3lo
          xf3_min=spm2hi
          rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
          rkappa5low=(xf3_max-xf3_min)*2.d0*spm2hi
          rk5=rkappa5high/rkappa6low*auxalfa6
ccc          auxalfa5=rk5/(rk5+1.d0)
          auxalfa5=1.d0/3.d0
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5low*auxalfa5
ccc          auxalfa4=rk4/(rk4+1.d0)
          auxalfa4=1.d0/4.d0
          xf2_max=spm2lo
          xf2_min=spm1hi
          rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          rk3=rkappa3high/rkappa4low*auxalfa4
ccc          auxalfa3=rk3/(rk3+1.d0)
          auxalfa3=1.d0/5.d0
          xw1_max=atan((spm1hi**2-rmx1**2)/(gamx1*rmx1))
          xw1_min=atan((estrinf**2-rmx1**2)/(gamx1*rmx1))
          rkappa2high=(xw1_max-xw1_min)*
     &          ((spm1hi**2-rmx1**2)**2+(gamx1*rmx1)**2)
     &                /(gamx1*rmx1)
          rk2=rkappa2high/rkappa3low*auxalfa3
ccc          auxalfa2=rk2/(rk2+1.d0)
          auxalfa2=1.d0/6.d0
          alfa2=auxalfa2
          alfa3=auxalfa3*(1.d0-auxalfa2)+alfa2
          alfa4=auxalfa4*(1.d0-auxalfa3)*(1.d0-auxalfa2)+alfa3
          alfa5=auxalfa5*(1.d0-auxalfa4)*(1.d0-auxalfa3)*
     &                          (1.d0-auxalfa2)+alfa4
          alfa6=auxalfa6*(1.d0-auxalfa5)*(1.d0-auxalfa4)*
     &                       (1.d0-auxalfa3)*(1.d0-auxalfa2)+alfa5
        ELSE  ! estrinf.LT.spm1lo  flat1 + res1 + flat2 + res2 + flat3 + res3 + flat4
          xf4_max=estrsup
          xf4_min=spm3hi
          rkappa7=(xf4_max-xf4_min)*2.d0*spm3hi
          xw3_min=atan((spm3lo**2-rmx3**2)/(gamx3*rmx3))
          xw3_max=atan((spm3hi**2-rmx3**2)/(gamx3*rmx3))
          rkappa6high=(xw3_max-xw3_min)*
     &          ((spm3hi**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rkappa6low=(xw3_max-xw3_min)*
     &          ((spm3lo**2-rmx3**2)**2+(gamx3*rmx3)**2)
     &                /(gamx3*rmx3)
          rk6=rkappa6high/rkappa7
ccc          auxalfa6=rk6/(rk6+1.d0)
          auxalfa6=0.5d0
          xf3_max=spm3lo
          xf3_min=spm2hi
          rkappa5high=(xf3_max-xf3_min)*2.d0*spm3lo
          rkappa5low=(xf3_max-xf3_min)*2.d0*spm2hi
          rk5=rkappa5high/rkappa6low*auxalfa6
ccc          auxalfa5=rk5/(rk5+1.d0)
          auxalfa5=1.d0/3.d0
          xw2_max=atan((spm2hi**2-rmx2**2)/(gamx2*rmx2))
          xw2_min=atan((spm2lo**2-rmx2**2)/(gamx2*rmx2))
          rkappa4high=(xw2_max-xw2_min)*
     &          ((spm2hi**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rkappa4low=(xw2_max-xw2_min)*
     &          ((spm2lo**2-rmx2**2)**2+(gamx2*rmx2)**2)
     &                /(gamx2*rmx2)
          rk4=rkappa4high/rkappa5low*auxalfa5
ccc          auxalfa4=rk4/(rk4+1.d0)
          auxalfa4=1.d0/4.d0
          xf2_max=spm2lo
          xf2_min=spm1hi
          rkappa3high=(xf2_max-xf2_min)*2.d0*spm2lo
          rkappa3low=(xf2_max-xf2_min)*2.d0*spm1hi
          rk3=rkappa3high/rkappa4low*auxalfa4
ccc          auxalfa3=rk3/(rk3+1.d0)
          auxalfa3=1.d0/5.d0
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
          auxalfa2=1.d0/6.d0
          xf1_max=spm1lo
          xf1_min=estrinf
          rkappa1=(xf1_max-xf1_min)*2.d0*spm1lo
          rk1=rkappa1/rkappa2low*auxalfa2
ccc          alfa1=rk1/(rk1+1.d0)
          alfa1=1.d0/7.d0
          alfa2=auxalfa2*(1.d0-alfa1)+alfa1
          alfa3=auxalfa3*(1.d0-auxalfa2)*(1.d0-alfa1)+alfa2
          alfa4=auxalfa4*(1.d0-auxalfa3)*(1.d0-auxalfa2)*(1.d0-alfa1)
     &                      +alfa3
          alfa5=auxalfa5*(1.d0-auxalfa4)*(1.d0-auxalfa3)*
     &               (1.d0-auxalfa2)*(1.d0-alfa1)+alfa4
          alfa6=auxalfa6*(1.d0-auxalfa5)*(1.d0-auxalfa4)*
     &         (1.d0-auxalfa3)*(1.d0-auxalfa2)*(1.d0-alfa1)+alfa5
        ENDIF



      ENDIF

      RETURN
      END



