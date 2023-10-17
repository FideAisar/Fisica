      IMPLICIT REAL*8 (a-b,d-h,o-z)
      READ(*,*)gamh
      xx=xhiggswidth(gamh)
      WRITE(*,*)'gamh=',gamh,'   xhiggswidth=',xx
      END

      
      REAL*8 FUNCTION xhiggswidth(gamh)
      IMPLICIT REAL*8 (a-b,d-h,o-z)
c Gamh(120 GeV) = 3.5d-3, Gamh(150 GeV) = 1.7d-2, Gamh(250 GeV) = 4.d0, Gamh(1000 GeV) = 650
      IF(gamh.LE.0.0d0)THEN
        xhiggswidth=-1.d0
      ELSEIF(gamh.LT.1.7d-2)THEN
        xhiggswidth=30.d0+(gamh-3.5d-3)/(1.7d-2-3.5d-3)*(20.d0-30.d0)
      ELSEIF(gamh.LT.4.0d0)THEN
        xhiggswidth=20.d0+(gamh-1.7d-2)/(4.0d0-1.7d-2)*(5.d0-20.d0)
      ELSEIF(gamh.LT.650.0d0)THEN
        xhiggswidth=5.d0+(gamh-4.0d0)/(650.0d0-4.0d0)*(1.d0-5.d0)
      ELSE
        xhiggswidth=-1.d0
      END IF
      RETURN
      END
