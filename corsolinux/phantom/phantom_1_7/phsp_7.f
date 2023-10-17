*****************************************************************************
* Auxiliary subroutine. Determines the invariant mass squared of a branching.
* The mass is supposed to resonate at three different masses rmx1, rmx2 and rmx3.
* INPUT:
*   xvar: a random number in [0,1]
*   rmx1,gamx1     : mass and width of first resonance
*   rmx2,gamx2     : mass and width of second resonance
*   rmx3,gamx3     : mass and width of third resonance
*   iregion        : index for filling the iphsp_fxn COMMON
*
*  alfa1,alfa2,alfa3,alfa4,alfa5,alfa6:
*                      separation points in [0,1] between the five regions
*                      Can also be 0 or 1. Some care necessary in this cases
*  auxalfa2,xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,xw2_min,xw2_max,
*  xf3_min,xf3_max,auxalfa3,auxalfa4,xw3_min,xw3_max,
*  xf4_min,xf4_max,auxalfa5,auxalfa6
*                    : auxiliary quantities used for generating the mass
*                      and computing the jacobian
* OUTPUT:
*   svar: the generated invariant mass squared
*   xjacvar: the jacobian
*****************************************************************************
      SUBROUTINE PHSP_7(xvar,
     &                     rmx1,gamx1,rmx2,gamx2,rmx3,gamx3,
     &                     alfa1,alfa2,auxalfa2,
     &                     alfa3,auxalfa3,alfa4,auxalfa4,
     &                     alfa5,auxalfa5,alfa6,auxalfa6,
     &               xf1_min,xf1_max,xw1_min,xw1_max,xf2_min,xf2_max,
     &               xw2_min,xw2_max,xf3_min,xf3_max,
     &               xw3_min,xw3_max,xf4_min,xf4_max,
     &                     iregion,
     &                     svar,xjacvar)

      IMPLICIT REAL*8(a-h,o-z)
c giuseppe 29/09/2007
      INTEGER maxnresonant
      PARAMETER (maxnresonant=4)
      INTEGER isresonant(maxnresonant),imapregion(maxnresonant)
      COMMON /iphsp_fxn/ isresonant,imapregion
      PARAMETER (pi=3.141592653589793238462643d0,twopi=2.d0*pi,
     &           fourpi=4.d0*pi)

* Compute the invariant mass squared and the corresponding jacobian

      xjacvar=1.d0

      IF(xvar.LE.alfa1.AND.alfa1.GT.0.d0)THEN
        vegas_x=xvar/alfa1
        evar=(xf1_max-xf1_min)*vegas_x+xf1_min
        svar=evar*evar
        IF(svar.LE.0.d0)THEN
          xjacvar=0.d0
          RETURN
        ENDIF
        xjacvar=xjacvar*(xf1_max-xf1_min)*2.d0*evar/alfa1
      ELSEIF(xvar.LE.alfa2.AND.alfa2.GT.0.d0)THEN
        vegas_x=(xvar-alfa1)/(1.d0-alfa1)/auxalfa2
        x_max=xw1_max
        x_min=xw1_min
        svar=(gamx1*rmx1*tan((x_max-x_min)*vegas_x+x_min)+rmx1**2)
        IF(svar.LE.0.d0)THEN
          xjacvar=0.d0
          RETURN
        ENDIF
        xjacvar=xjacvar*(x_max-x_min)*((svar-rmx1**2)**2+
     &      (gamx1*rmx1)**2)/(gamx1*rmx1)/(1.d0-alfa1)/auxalfa2
c giuseppe 29/09/2007
        isresonant(iregion)=1
        imapregion(iregion)=2
c end giuseppe 29/09/2007
      ELSEIF(xvar.LE.alfa3.AND.alfa3.GT.0.d0)THEN
        vegas_x=(xvar-alfa2)/(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
        evar=(xf2_max-xf2_min)*vegas_x+xf2_min
        svar=evar*evar
        IF(svar.LE.0.d0)THEN
          xjacvar=0.d0
          RETURN
        ENDIF
        xjacvar=xjacvar*(xf2_max-xf2_min)*2.d0*evar
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)/auxalfa3
      ELSEIF(xvar.LE.alfa4.AND.alfa4.GT.0.d0)THEN
        vegas_x=(xvar-alfa3)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/auxalfa4
        x_max=xw2_max
        x_min=xw2_min
        svar=(gamx2*rmx2*tan((x_max-x_min)*vegas_x+x_min)+rmx2**2)
        IF(svar.LE.0.d0)THEN
          xjacvar=0.d0
          RETURN
        ENDIF
        xjacvar=xjacvar*(x_max-x_min)*((svar-rmx2**2)**2+
     &      (gamx2*rmx2)**2)/(gamx2*rmx2)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/auxalfa4
c giuseppe 29/09/2007
        isresonant(iregion)=1
        imapregion(iregion)=4
c end giuseppe 29/09/2007
      ELSEIF(xvar.LE.alfa5.AND.alfa5.GT.0.d0)THEN
        vegas_x=(xvar-alfa4)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)/auxalfa5
        evar=(xf3_max-xf3_min)*vegas_x+xf3_min
        svar=evar*evar
        IF(svar.LE.0.d0)THEN
          xjacvar=0.d0
          RETURN
        ENDIF
        xjacvar=xjacvar*(xf3_max-xf3_min)*2.d0*evar
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)/auxalfa5
      ELSEIF(xvar.LE.alfa6.AND.alfa6.GT.0.d0)THEN
        vegas_x=(xvar-alfa5)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
     &                      /(1.d0-auxalfa5)/auxalfa6
        x_max=xw3_max
        x_min=xw3_min
        svar=(gamx3*rmx3*tan((x_max-x_min)*vegas_x+x_min)+rmx3**2)
        IF(svar.LE.0.d0)THEN
          xjacvar=0.d0
          RETURN
        ENDIF
        xjacvar=xjacvar*(x_max-x_min)*((svar-rmx3**2)**2+
     &      (gamx3*rmx3)**2)/(gamx3*rmx3)/(1.d0-alfa1)
     &       /(1.d0-auxalfa2)/(1.d0-auxalfa3)/(1.d0-auxalfa4)
     &                      /(1.d0-auxalfa5)/auxalfa6
c giuseppe 29/09/2007
        isresonant(iregion)=1
        imapregion(iregion)=6
c end giuseppe 29/09/2007
      ELSE
        vegas_x=(xvar-alfa6)/(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
     &                      /(1.d0-auxalfa5)/(1.d0-auxalfa6)
        evar=(xf4_max-xf4_min)*vegas_x+xf4_min
        svar=evar*evar
        IF(svar.LE.0.d0)THEN
          xjacvar=0.d0
          RETURN
        ENDIF
        xjacvar=xjacvar*(xf4_max-xf4_min)*2.d0*evar
     &                  /(1.d0-alfa1)/(1.d0-auxalfa2)
     &                      /(1.d0-auxalfa3)/(1.d0-auxalfa4)
     &                      /(1.d0-auxalfa5)/(1.d0-auxalfa6)
      ENDIF
      xjacvar=xjacvar/twopi

      RETURN
      END


