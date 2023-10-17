***********************************************************************
*                        PROGRAM BR_H2hh                              *
*                                                                     *
*                                                                     *
*  Purpose:  It contains DATA (for masses, W and Z widths and         *
*                                                     alfa_couplings) *
*            It computes the Branching Ratio of a Heavy Higgs into    *
*            a pair of light Higgses in the Singlet Model.            *
*            Uses Passarino's parametrization of the SM light Higgs   *
*            multiplied by cos(alpha)**2 and adds H-->hh at tree      *
*            level.                                                   *
*                                                                     *
***********************************************************************


      PROGRAM BR_H2hh

      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)

      DATA rmw/80.40d0/, rmz/91.187d0/, rmt/175.d0/, rmb/4.8d0/,
     &     rmc/0.75d0/,rmtau/1.78d0/,
     &     gamw/0.2042774d+01/, gamz/0.25007d+01/
      DATA gf/1.16639d-05/, alfainv_me/137.0359895d0/,
     &  alfas_w/0.1225d0/, alfas_z/0.123d0/,
     &  alfas_t/0.100d0/, alfas_h/0.100d0/

      DATA pi/3.141592653589793238462643d0/

      DATA czero/(0.d0,0.d0)/, cuno/(1.d0,0.d0)/, cim/(0.d0,1.d0)/

      rmh=125.d0
      tgbeta=1.d0


* sin(alpha)=0.1
      ghhfactor=0.1d0
      ghhfactor2=ghhfactor*ghhfactor
      ghfactor2=1.d0-ghhfactor2
      ghfactor = sqrt(ghfactor2)

      OPEN(UNIT=51,FILE='BRH2hh_01.dat')
      WRITE(51,*)'sin(alpha):',ghhfactor,'  tg(beta):',tgbeta
      WRITE(51,*)'rmh:',rmh

      DO i=250,1000
        rmhh=i
        call gridHt(rmhh,gamhh)
        gamhh=gamhh*ghhfactor2

* Add H->hh channel
        alfa=1.d0/alfainv_me
        esquared=alfa*4.d0*pi
        gamH2hh=esquared*rmhh**3/
     &      (128.d0*pi*rmw**2*(1.d0-rmw**2/rmz**2))*
     &      sqrt(1.d0-4.d0*rmh**2/rmhh**2)*
     &      (1.d0+2.d0*rmh**2/rmhh**2)**2*
     &      ghfactor**2*ghhfactor**2*(ghfactor+ghhfactor*tgbeta)**2
        gamhh=gamhh+gamH2hh
        BR=gamH2hh/gamhh
        WRITE(51,*)rmhh,BR,gamH2hh,gamhh
      ENDDO
      CLOSE(51)


* sin(alpha)=0.2
      ghhfactor=0.2d0
      ghhfactor2=ghhfactor*ghhfactor
      ghfactor2=1.d0-ghhfactor2
      ghfactor = sqrt(ghfactor2)

      OPEN(UNIT=51,FILE='BRH2hh_02.dat')
      WRITE(51,*)'sin(alpha):',ghhfactor,'  tg(beta):',tgbeta
      WRITE(51,*)'rmh:',rmh

      DO i=250,1000
        rmhh=i
        call gridHt(rmhh,gamhh)
        gamhh=gamhh*ghhfactor2

* Add H->hh channel
        alfa=1.d0/alfainv_me
        esquared=alfa*4.d0*pi
        gamH2hh=esquared*rmhh**3/
     &      (128.d0*pi*rmw**2*(1.d0-rmw**2/rmz**2))*
     &      sqrt(1.d0-4.d0*rmh**2/rmhh**2)*
     &      (1.d0+2.d0*rmh**2/rmhh**2)**2*
     &      ghfactor**2*ghhfactor**2*(ghfactor+ghhfactor*tgbeta)**2
        gamhh=gamhh+gamH2hh
        BR=gamH2hh/gamhh
        WRITE(51,*)rmhh,BR,gamH2hh,gamhh
      ENDDO
      CLOSE(51)



* sin(alpha)=0.3
      ghhfactor=0.3d0
      ghhfactor2=ghhfactor*ghhfactor
      ghfactor2=1.d0-ghhfactor2
      ghfactor = sqrt(ghfactor2)

      OPEN(UNIT=51,FILE='BRH2hh_03.dat')
      WRITE(51,*)'sin(alpha):',ghhfactor,'  tg(beta):',tgbeta
      WRITE(51,*)'rmh:',rmh

      DO i=250,1000
        rmhh=i
        call gridHt(rmhh,gamhh)
        gamhh=gamhh*ghhfactor2

* Add H->hh channel
        alfa=1.d0/alfainv_me
        esquared=alfa*4.d0*pi
        gamH2hh=esquared*rmhh**3/
     &      (128.d0*pi*rmw**2*(1.d0-rmw**2/rmz**2))*
     &      sqrt(1.d0-4.d0*rmh**2/rmhh**2)*
     &      (1.d0+2.d0*rmh**2/rmhh**2)**2*
     &      ghfactor**2*ghhfactor**2*(ghfactor+ghhfactor*tgbeta)**2
        gamhh=gamhh+gamH2hh
        BR=gamH2hh/gamhh
        WRITE(51,*)rmhh,BR,gamH2hh,gamhh
      ENDDO
      CLOSE(51)



      END



*    The following routines have been written by Giampiero Passarino
*
*
*
*----------------------------------------------------------------
* H total grids
*----------------------------------------------------------------
*
      SUBROUTINE gridHt(mass,evalue)
* 
      IMPLICIT NONE
*
      INTEGER i,top,gdim
      REAL*8 u,value,evalue,mass
      REAL*8, dimension(105) :: bc,cc,dc
* 
* u value of M_H at which the spline is to be evaluated
* top= -1,0,1 lower, central, upper value for m_top
*
      gdim= 105
*
      CALL FMMsplineSingleHt(bc,cc,dc,top,gdim)
*
      u= mass
      CALL Seval3SingleHt(u,bc,cc,dc,top,gdim,value)
*
      evalue= value
*
      RETURN
*
*-----------------------------------------------------------------------
*
      CONTAINS
*
      SUBROUTINE FMMsplineSingleHt(b,c,d,top,gdim)
*
*---------------------------------------------------------------------------
*
      INTEGER k,n,i,top,gdim,l
*
      REAL*8, dimension(105) :: xc,yc
      REAL*8, dimension(105) :: x,y
*
      REAL*8, DIMENSION(gdim) :: b 
* linear coeff
*
      REAL*8, DIMENSION(gdim) :: c 
* quadratic coeff.
*
      REAL*8, DIMENSION(gdim) :: d 
* cubic coeff.
*
      REAL*8 :: t
      REAL*8,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
*
* The grid
*
*
      DATA (xc(i),i=1,105)/
     #   90.d0,95.d0,100.d0,105.d0,110.d0,115.d0,120.d0,
     #   125.d0,130.d0,135.d0,140.d0,145.d0,150.d0,155.d0,160.d0,165.d0,
     #   170.d0,175.d0,180.d0,185.d0,190.d0,195.d0,200.d0,210.d0,220.d0,
     #   230.d0,240.d0,250.d0,260.d0,270.d0,280.d0,290.d0,300.d0,310.d0,
     #   320.d0,330.d0,340.d0,350.d0,360.d0,370.d0,380.d0,390.d0,400.d0,
     #   410.d0,420.d0,430.d0,440.d0,450.d0,460.d0,470.d0,480.d0,490.d0,
     #   500.d0,510.d0,520.d0,530.d0,540.d0,550.d0,560.d0,570.d0,580.d0,
     #   590.d0,600.d0,610.d0,620.d0,630.d0,640.d0,650.d0,660.d0,670.d0,
     #   680.d0,690.d0,700.d0,710.d0,720.d0,730.d0,740.d0,750.d0,760.d0,
     #   770.d0,780.d0,790.d0,800.d0,810.d0,820.d0,830.d0,840.d0,850.d0,
     #   860.d0,870.d0,880.d0,890.d0,900.d0,910.d0,920.d0,930.d0,940.d0,
     #   950.d0,960.d0,970.d0,980.d0,990.d0,1000.d0,1500.d0,2000.d0/
*
      DATA (yc(i),i=1,105)/
     # 2.20d-3,2.32d-3,2.46d-3,2.62d-3,2.82d-3,3.09d-3,3.47d-3,4.03d-3,
     # 4.87d-3,6.14d-3,8.12d-3,1.14d-2,1.73d-2,3.02d-2,8.29d-2,2.46d-1,
     # 3.80d-1,5.00d-1,6.31d-1,8.32d-1,1.04d0,1.24d0,1.43d0,1.85d0,
     # 2.31d0,2.82d0,3.40d0,4.04d0,4.76d0,5.55d0,6.43d0,7.39d0,8.43d0,
     # 9.57d0,10.8d0,12.1d0,13.5d0,15.2d0,17.6d0,20.2d0,23.1d0,26.1d0,
     # 29.2d0,32.5d0,35.9d0,39.4d0,43.1d0,46.9d0,50.8d0,54.9d0,59.1d0,
     # 63.5d0,68.0d0,72.7d0,77.6d0,82.6d0,87.7d0,93.1d0,98.7d0,104.d0,
     # 110.d0,116.d0,123.d0,129.d0,136.d0,143.d0,150.d0,158.d0,166.d0,
     # 174.d0,182.d0,190.d0,199.d0,208.d0,218.d0,227.d0,237.d0,248.d0,
     # 258.d0,269.d0,281.d0,292.d0,304.d0,317.d0,330.d0,343.d0,357.d0,
     # 371.d0,386.d0,401.d0,416.d0,432.d0,449.d0,466.d0,484.d0,502.d0,
     # 521.d0,540.d0,560.d0,581.d0,602.d0,624.d0,647.d0,3.38d3,1.58d4/
*
      n= 105
      FORALL(l=1:105)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL

*.....Set up tridiagonal system.........................................
*     b=diagonal, d=offdiagonal, c=right-hand side
*
      d(1)= x(2)-x(1)
      c(2)= (y(2)-y(1))/d(1)
      DO k= 2,n-1
       d(k)= x(k+1)-x(k)
       b(k)= TWO*(d(k-1)+d(k))
       c(k+1)= (y(k+1)-y(k))/d(k)
       c(k)= c(k+1)-c(k)
      END DO
*
*.....End conditions.  third derivatives at x(1) and x(n) obtained
*     from divided differences.......................................
*
      b(1)= -d(1)
      b(n)= -d(n-1)
      c(1)= ZERO
      c(n)= ZERO
      IF (n > 3) THEN
       c(1)= c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
       c(n)= c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
       c(1)= c(1)*d(1)*d(1)/(x(4)-x(1))
       c(n)= -c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
      END IF
*
      DO k=2,n    ! forward elimination
       t= d(k-1)/b(k-1)
       b(k)= b(k)-t*d(k-1)
       c(k)= c(k)-t*c(k-1)
      END DO
*
      c(n)= c(n)/b(n)   
*
* back substitution ( makes c the sigma of text)
*
      DO k=n-1,1,-1
       c(k)= (c(k)-d(k)*c(k+1))/b(k)
      END DO
*
*.....Compute polynomial coefficients...................................
*
      b(n)= (y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
      DO k=1,n-1
       b(k)= (y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
       d(k)= (c(k+1)-c(k))/d(k)
       c(k)= THREE*c(k)
      END DO
      c(n)= THREE*c(n)
      d(n)= d(n-1)
*
      RETURN
*
      END Subroutine FMMsplineSingleHt   
*
*------------------------------------------------------------------------
*
      SUBROUTINE Seval3SingleHt(u,b,c,d,top,gdim,f,fp,fpp,fppp)
*
* ---------------------------------------------------------------------------
*
      REAL*8,INTENT(IN) :: u 
* abscissa at which the spline is to be evaluated
*
      INTEGER j,k,n,l,top,gdim
*
      REAL*8, dimension(105) :: xc,yc
      REAL*8, dimension(105) :: x,y
      REAL*8, DIMENSION(gdim) :: b,c,d 
* linear,quadratic,cubic coeff
*
      REAL*8,INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp 
* function, 1st,2nd,3rd deriv
*
      INTEGER, SAVE :: i=1
      REAL*8    :: dx
      REAL*8,PARAMETER:: TWO=2.0, THREE=3.0, SIX=6.0
*
* The grid
*
      DATA (xc(l),l=1,105)/
     #   90.d0,95.d0,100.d0,105.d0,110.d0,115.d0,120.d0,
     #   125.d0,130.d0,135.d0,140.d0,145.d0,150.d0,155.d0,160.d0,165.d0,
     #   170.d0,175.d0,180.d0,185.d0,190.d0,195.d0,200.d0,210.d0,220.d0,
     #   230.d0,240.d0,250.d0,260.d0,270.d0,280.d0,290.d0,300.d0,310.d0,
     #   320.d0,330.d0,340.d0,350.d0,360.d0,370.d0,380.d0,390.d0,400.d0,
     #   410.d0,420.d0,430.d0,440.d0,450.d0,460.d0,470.d0,480.d0,490.d0,
     #   500.d0,510.d0,520.d0,530.d0,540.d0,550.d0,560.d0,570.d0,580.d0,
     #   590.d0,600.d0,610.d0,620.d0,630.d0,640.d0,650.d0,660.d0,670.d0,
     #   680.d0,690.d0,700.d0,710.d0,720.d0,730.d0,740.d0,750.d0,760.d0,
     #   770.d0,780.d0,790.d0,800.d0,810.d0,820.d0,830.d0,840.d0,850.d0,
     #   860.d0,870.d0,880.d0,890.d0,900.d0,910.d0,920.d0,930.d0,940.d0,
     #   950.d0,960.d0,970.d0,980.d0,990.d0,1000.d0,1500.0d0,2000.d0/
*
      DATA (yc(l),l=1,105)/
     # 2.20d-3,2.32d-3,2.46d-3,2.62d-3,2.82d-3,3.09d-3,3.47d-3,4.03d-3,
     # 4.87d-3,6.14d-3,8.12d-3,1.14d-2,1.73d-2,3.02d-2,8.29d-2,2.46d-1,
     # 3.80d-1,5.00d-1,6.31d-1,8.32d-1,1.04d0,1.24d0,1.43d0,1.85d0,
     # 2.31d0,2.82d0,3.40d0,4.04d0,4.76d0,5.55d0,6.43d0,7.39d0,8.43d0,
     # 9.57d0,10.8d0,12.1d0,13.5d0,15.2d0,17.6d0,20.2d0,23.1d0,26.1d0,
     # 29.2d0,32.5d0,35.9d0,39.4d0,43.1d0,46.9d0,50.8d0,54.9d0,59.1d0,
     # 63.5d0,68.0d0,72.7d0,77.6d0,82.6d0,87.7d0,93.1d0,98.7d0,104.d0,
     # 110.d0,116.d0,123.d0,129.d0,136.d0,143.d0,150.d0,158.d0,166.d0,
     # 174.d0,182.d0,190.d0,199.d0,208.d0,218.d0,227.d0,237.d0,248.d0,
     # 258.d0,269.d0,281.d0,292.d0,304.d0,317.d0,330.d0,343.d0,357.d0,
     # 371.d0,386.d0,401.d0,416.d0,432.d0,449.d0,466.d0,484.d0,502.d0,
     # 521.d0,540.d0,560.d0,581.d0,602.d0,624.d0,647.d0,3.38d3,1.58d4/
*
      n= 105
      FORALL(l=1:105)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL
*
*.....First check if u is in the same interval found on the
*     last call to Seval.............................................
*
      IF (  (i<1) .OR. (i >= n) ) i=1
      IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
       i=1   
*
* binary search
*
       j= n+1
       DO
        k= (i+j)/2
        IF (u < x(k)) THEN
         j= k
        ELSE
         i= k
        ENDIF
        IF (j <= i+1) EXIT
       ENDDO
      ENDIF
*
      dx= u-x(i)   
*
* evaluate the spline
*
      IF (Present(f))    f= y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      IF (Present(fp))   fp= b(i)+dx*(TWO*c(i) + dx*THREE*d(i))
      IF (Present(fpp))  fpp= TWO*c(i) + dx*SIX*d(i)
      IF (Present(fppp)) fppp= SIX*d(i)
*
      RETURN
*
      END Subroutine Seval3SingleHt  
*
      END SUBROUTINE gridHt
*
*----------------------------------------------------------------
* ZZ YR1 grids
*----------------------------------------------------------------
*
      SUBROUTINE gridYZZ(mass,evalue)
* 
      IMPLICIT NONE
*
      INTEGER i,top,gdim
      REAL*8 u,value,evalue,mass
      REAL*8, dimension(109) :: bc,cc,dc
* 
* u value of M_H at which the spline is to be evaluated
* top= -1,0,1 lower, central, upper value for m_top
*
      gdim= 109
*
      CALL FMMsplineSingleYZZ(bc,cc,dc,top,gdim)
*
      u= mass
      CALL Seval3SingleYZZ(u,bc,cc,dc,top,gdim,value)
*
      evalue= value
*
      RETURN
*
*-----------------------------------------------------------------------
*
      CONTAINS
*
      SUBROUTINE FMMsplineSingleYZZ(b,c,d,top,gdim)
*
*---------------------------------------------------------------------------
*
      INTEGER k,n,i,top,gdim,l
*
      REAL*8, dimension(109) :: xc,yc
      REAL*8, dimension(109) :: x,y
*
      REAL*8, DIMENSION(gdim) :: b 
* linear coeff
*
      REAL*8, DIMENSION(gdim) :: c 
* quadratic coeff.
*
      REAL*8, DIMENSION(gdim) :: d 
* cubic coeff.
*
      REAL*8 :: t
      REAL*8,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
*
* The grid
*
*
      DATA (xc(i),i=1,109)/
     #   90.d0,95.d0,100.d0,105.d0,110.d0,115.d0,120.d0,
     #   125.d0,130.d0,135.d0,140.d0,145.d0,150.d0,155.d0,160.d0,165.d0,
     #   170.d0,175.d0,180.d0,185.d0,190.d0,195.d0,200.d0,210.d0,220.d0,
     #   230.d0,240.d0,250.d0,260.d0,270.d0,280.d0,290.d0,300.d0,310.d0,
     #   320.d0,330.d0,340.d0,350.d0,360.d0,370.d0,380.d0,390.d0,400.d0,
     #   410.d0,420.d0,430.d0,440.d0,450.d0,460.d0,470.d0,480.d0,490.d0,
     #   500.d0,510.d0,520.d0,530.d0,540.d0,550.d0,560.d0,570.d0,580.d0,
     #   590.d0,600.d0,610.d0,620.d0,630.d0,640.d0,650.d0,660.d0,670.d0,
     #   680.d0,690.d0,700.d0,710.d0,720.d0,730.d0,740.d0,750.d0,760.d0,
     #   770.d0,780.d0,790.d0,800.d0,810.d0,820.d0,830.d0,840.d0,850.d0,
     #   860.d0,870.d0,880.d0,890.d0,900.d0,910.d0,920.d0,930.d0,940.d0,
     #   950.d0,960.d0,970.d0,980.d0,990.d0,1000.d0,
     #   1.2d3,1.4d3,1.5d3,1.6d3,1.8d3,2.0d3/
*
      DATA (yc(i),i=1,109)/
     # 9.27d-7,1.56d-6,2.79d-6,5.63d-6,1.25d-5,2.70d-5,5.57d-5,1.07d-4,
     # 1.95d-4,3.38d-4,5.62d-4,9.06d-4,1.43d-3,2.22d-3,3.44d-3,5.47d-3,
     # 8.98d-3,1.62d-2,3.80d-2,1.25d-1,2.18d-1,2.95d-1,3.66d-1,5.06d-1,
     # 6.54d-1,8.16d-1,9.97d-1,1.20d0,1.42d0,1.67d0,1.95d0,2.25d0,
     # 2.59d0,2.95d0,3.34d0,3.76d0,4.19d0,4.66d0,5.22d0,5.81d0,6.45d0,
     # 7.13d0,7.85d0,8.62d0,9.43d0,10.3d0,11.2d0,12.2d0,13.2d0,14.2d0,
     # 15.3d0,16.5d0,17.7d0,19.0d0,20.4d0,21.8d0,23.2d0,24.7d0,26.3d0,
     # 28.0d0,29.7d0,31.5d0,33.4d0,35.4d0,37.4d0,39.5d0,41.7d0,44.0d0,
     # 46.3d0,48.8d0,51.3d0,53.9d0,56.7d0,59.5d0,62.4d0,65.5d0,68.6d0,
     # 71.9d0,75.2d0,78.7d0,82.3d0,86.1d0,89.9d0,93.9d0,98.0d0,102.d0,
     # 107.d0,111.d0,116.d0,121.d0,126.d0,131.d0,137.d0,142.d0,148.d0,
     # 154.d0,160.d0,166.d0,173.d0,180.d0,187.d0,194.d0,201.2d0,
     # 413.5d0,813.5d0,1128.6d0,1555.6d0,2898.5d0,5253.4d0/
*
      n= 109
      FORALL(l=1:109)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL

*.....Set up tridiagonal system.........................................
*     b=diagonal, d=offdiagonal, c=right-hand side
*
      d(1)= x(2)-x(1)
      c(2)= (y(2)-y(1))/d(1)
      DO k= 2,n-1
       d(k)= x(k+1)-x(k)
       b(k)= TWO*(d(k-1)+d(k))
       c(k+1)= (y(k+1)-y(k))/d(k)
       c(k)= c(k+1)-c(k)
      END DO
*
*.....End conditions.  third derivatives at x(1) and x(n) obtained
*     from divided differences.......................................
*
      b(1)= -d(1)
      b(n)= -d(n-1)
      c(1)= ZERO
      c(n)= ZERO
      IF (n > 3) THEN
       c(1)= c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
       c(n)= c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
       c(1)= c(1)*d(1)*d(1)/(x(4)-x(1))
       c(n)= -c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
      END IF
*
      DO k=2,n    ! forward elimination
       t= d(k-1)/b(k-1)
       b(k)= b(k)-t*d(k-1)
       c(k)= c(k)-t*c(k-1)
      END DO
*
      c(n)= c(n)/b(n)   
*
* back substitution ( makes c the sigma of text)
*
      DO k=n-1,1,-1
       c(k)= (c(k)-d(k)*c(k+1))/b(k)
      END DO
*
*.....Compute polynomial coefficients...................................
*
      b(n)= (y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
      DO k=1,n-1
       b(k)= (y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
       d(k)= (c(k+1)-c(k))/d(k)
       c(k)= THREE*c(k)
      END DO
      c(n)= THREE*c(n)
      d(n)= d(n-1)
*
      RETURN
*
      END Subroutine FMMsplineSingleYZZ   
*
*------------------------------------------------------------------------
*
      SUBROUTINE Seval3SingleYZZ(u,b,c,d,top,gdim,f,fp,fpp,fppp)
*
* ---------------------------------------------------------------------------
*
      REAL*8,INTENT(IN) :: u 
* abscissa at which the spline is to be evaluated
*
      INTEGER j,k,n,l,top,gdim
*
      REAL*8, dimension(109) :: xc,yc
      REAL*8, dimension(109) :: x,y
      REAL*8, DIMENSION(gdim) :: b,c,d 
* linear,quadratic,cubic coeff
*
      REAL*8,INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp 
* function, 1st,2nd,3rd deriv
*
      INTEGER, SAVE :: i=1
      REAL*8    :: dx
      REAL*8,PARAMETER:: TWO=2.0, THREE=3.0, SIX=6.0
*
* The grid
*
      DATA (xc(l),l=1,109)/
     #   90.d0,95.d0,100.d0,105.d0,110.d0,115.d0,120.d0,
     #   125.d0,130.d0,135.d0,140.d0,145.d0,150.d0,155.d0,160.d0,165.d0,
     #   170.d0,175.d0,180.d0,185.d0,190.d0,195.d0,200.d0,210.d0,220.d0,
     #   230.d0,240.d0,250.d0,260.d0,270.d0,280.d0,290.d0,300.d0,310.d0,
     #   320.d0,330.d0,340.d0,350.d0,360.d0,370.d0,380.d0,390.d0,400.d0,
     #   410.d0,420.d0,430.d0,440.d0,450.d0,460.d0,470.d0,480.d0,490.d0,
     #   500.d0,510.d0,520.d0,530.d0,540.d0,550.d0,560.d0,570.d0,580.d0,
     #   590.d0,600.d0,610.d0,620.d0,630.d0,640.d0,650.d0,660.d0,670.d0,
     #   680.d0,690.d0,700.d0,710.d0,720.d0,730.d0,740.d0,750.d0,760.d0,
     #   770.d0,780.d0,790.d0,800.d0,810.d0,820.d0,830.d0,840.d0,850.d0,
     #   860.d0,870.d0,880.d0,890.d0,900.d0,910.d0,920.d0,930.d0,940.d0,
     #   950.d0,960.d0,970.d0,980.d0,990.d0,1000.d0,
     #   1.2d3,1.4d3,1.5d3,1.6d3,1.8d3,2.0d3/
*
      DATA (yc(l),l=1,109)/
     # 9.27d-7,1.56d-6,2.79d-6,5.63d-6,1.25d-5,2.70d-5,5.57d-5,1.07d-4,
     # 1.95d-4,3.38d-4,5.62d-4,9.06d-4,1.43d-3,2.22d-3,3.44d-3,5.47d-3,
     # 8.98d-3,1.62d-2,3.80d-2,1.25d-1,2.18d-1,2.95d-1,3.66d-1,5.06d-1,
     # 6.54d-1,8.16d-1,9.97d-1,1.20d0,1.42d0,1.67d0,1.95d0,2.25d0,
     # 2.59d0,2.95d0,3.34d0,3.76d0,4.19d0,4.66d0,5.22d0,5.81d0,6.45d0,
     # 7.13d0,7.85d0,8.62d0,9.43d0,10.3d0,11.2d0,12.2d0,13.2d0,14.2d0,
     # 15.3d0,16.5d0,17.7d0,19.0d0,20.4d0,21.8d0,23.2d0,24.7d0,26.3d0,
     # 28.0d0,29.7d0,31.5d0,33.4d0,35.4d0,37.4d0,39.5d0,41.7d0,44.0d0,
     # 46.3d0,48.8d0,51.3d0,53.9d0,56.7d0,59.5d0,62.4d0,65.5d0,68.6d0,
     # 71.9d0,75.2d0,78.7d0,82.3d0,86.1d0,89.9d0,93.9d0,98.0d0,102.d0,
     # 107.d0,111.d0,116.d0,121.d0,126.d0,131.d0,137.d0,142.d0,148.d0,
     # 154.d0,160.d0,166.d0,173.d0,180.d0,187.d0,194.d0,201.2d0,
     # 413.5d0,813.5d0,1128.6d0,1555.6d0,2898.5d0,5253.4d0/
*
      n= 109
      FORALL(l=1:109)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL
*
*.....First check if u is in the same interval found on the
*     last call to Seval.............................................
*
      IF (  (i<1) .OR. (i >= n) ) i=1
      IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
       i=1   
*
* binary search
*
       j= n+1
       DO
        k= (i+j)/2
        IF (u < x(k)) THEN
         j= k
        ELSE
         i= k
        ENDIF
        IF (j <= i+1) EXIT
       ENDDO
      ENDIF
*
      dx= u-x(i)   
*
* evaluate the spline
*
      IF (Present(f))    f= y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      IF (Present(fp))   fp= b(i)+dx*(TWO*c(i) + dx*THREE*d(i))
      IF (Present(fpp))  fpp= TWO*c(i) + dx*SIX*d(i)
      IF (Present(fppp)) fppp= SIX*d(i)
*
      RETURN
*
      END Subroutine Seval3SingleYZZ  
*
      END SUBROUTINE gridYZZ
*
*----------------------------------------------------------------
* WW YR1 grids
*----------------------------------------------------------------
*
      SUBROUTINE gridYWW(mass,evalue)
* 
      IMPLICIT NONE
*
      INTEGER i,top,gdim
      REAL*8 u,value,evalue,mass
      REAL*8, dimension(109) :: bc,cc,dc
* 
* u value of M_H at which the spline is to be evaluated
* top= -1,0,1 lower, central, upper value for m_top
*
      gdim= 109
*
      CALL FMMsplineSingleYWW(bc,cc,dc,top,gdim)
*
      u= mass
      CALL Seval3SingleYWW(u,bc,cc,dc,top,gdim,value)
*
      evalue= value
*
      RETURN
*
*-----------------------------------------------------------------------
*
      CONTAINS
*
      SUBROUTINE FMMsplineSingleYWW(b,c,d,top,gdim)
*
*---------------------------------------------------------------------------
*
      INTEGER k,n,i,top,gdim,l
*
      REAL*8, dimension(109) :: xc,yc
      REAL*8, dimension(109) :: x,y
*
      REAL*8, DIMENSION(gdim) :: b 
* linear coeff
*
      REAL*8, DIMENSION(gdim) :: c 
* quadratic coeff.
*
      REAL*8, DIMENSION(gdim) :: d 
* cubic coeff.
*
      REAL*8 :: t
      REAL*8,PARAMETER:: ZERO=0.0, TWO=2.0, THREE=3.0
*
* The grid
*
*
      DATA (xc(i),i=1,109)/
     #   90.d0,95.d0,100.d0,105.d0,110.d0,115.d0,120.d0,
     #   125.d0,130.d0,135.d0,140.d0,145.d0,150.d0,155.d0,160.d0,165.d0,
     #   170.d0,175.d0,180.d0,185.d0,190.d0,195.d0,200.d0,210.d0,220.d0,
     #   230.d0,240.d0,250.d0,260.d0,270.d0,280.d0,290.d0,300.d0,310.d0,
     #   320.d0,330.d0,340.d0,350.d0,360.d0,370.d0,380.d0,390.d0,400.d0,
     #   410.d0,420.d0,430.d0,440.d0,450.d0,460.d0,470.d0,480.d0,490.d0,
     #   500.d0,510.d0,520.d0,530.d0,540.d0,550.d0,560.d0,570.d0,580.d0,
     #   590.d0,600.d0,610.d0,620.d0,630.d0,640.d0,650.d0,660.d0,670.d0,
     #   680.d0,690.d0,700.d0,710.d0,720.d0,730.d0,740.d0,750.d0,760.d0,
     #   770.d0,780.d0,790.d0,800.d0,810.d0,820.d0,830.d0,840.d0,850.d0,
     #   860.d0,870.d0,880.d0,890.d0,900.d0,910.d0,920.d0,930.d0,940.d0,
     #   950.d0,960.d0,970.d0,980.d0,990.d0,1000.d0,
     #   1.2d3,1.4d3,1.5d3,1.6d3,1.8d3,2.0d3/
*
      DATA (yc(i),i=1,109)/
     # 4.60d-6,1.10d-5,2.72d-5,6.36d-5,1.36d-4,2.68d-4,4.95d-4,8.73d-4,
     # 1.49d-3,2.47d-3,4.10d-3,6.86d-3,1.21d-2,2.41d-2,7.53d-2,2.36d-1,
     # 3.67d-1,4.80d-1,5.88d-1,7.02d-1,8.17d-1,9.36d-1,1.06d0,1.33d0,
     # 1.65d0,2.00d0,2.39d0,2.83d0,3.32d0,3.87d0,4.47d0,5.12d0,5.83d0,
     # 6.60d0,7.44d0,8.32d0,9.26d0,10.3d0,11.4d0,12.7d0,14.1d0,15.5d0,
     # 17.0d0,18.6d0,20.2d0,22.0d0,23.9d0,25.8d0,27.9d0,30.0d0,32.3d0,
     # 34.6d0,37.1d0,39.7d0,42.4d0,45.2d0,48.1d0,51.2d0,54.4d0,57.7d0,
     # 61.2d0,64.8d0,68.5d0,72.4d0,76.5d0,80.7d0,85.0d0,89.6d0,94.3d0,
     # 99.1d0,104.d0,109.d0,115.d0,120.d0,126.d0,132.d0,139.d0,145.d0,
     # 152.d0,159.d0,166.d0,173.d0,181.d0,189.d0,197.d0,206.d0,214.d0,
     # 223.d0,233.d0,242.d0,252.d0,263.d0,273.d0,284.d0,296.d0,308.d0,
     # 320.d0,332.d0,345.d0,359.d0,373.d0,387.d0,401.9d0,
     # 824.7d0,1623.7d0,2253.9d0,3108.6d0,5797.9d0,10517.2d0/
*
      n= 109
      FORALL(l=1:109)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL

*.....Set up tridiagonal system.........................................
*     b=diagonal, d=offdiagonal, c=right-hand side
*
      d(1)= x(2)-x(1)
      c(2)= (y(2)-y(1))/d(1)
      DO k= 2,n-1
       d(k)= x(k+1)-x(k)
       b(k)= TWO*(d(k-1)+d(k))
       c(k+1)= (y(k+1)-y(k))/d(k)
       c(k)= c(k+1)-c(k)
      END DO
*
*.....End conditions.  third derivatives at x(1) and x(n) obtained
*     from divided differences.......................................
*
      b(1)= -d(1)
      b(n)= -d(n-1)
      c(1)= ZERO
      c(n)= ZERO
      IF (n > 3) THEN
       c(1)= c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
       c(n)= c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
       c(1)= c(1)*d(1)*d(1)/(x(4)-x(1))
       c(n)= -c(n)*d(n-1)*d(n-1)/(x(n)-x(n-3))
      END IF
*
      DO k=2,n    ! forward elimination
       t= d(k-1)/b(k-1)
       b(k)= b(k)-t*d(k-1)
       c(k)= c(k)-t*c(k-1)
      END DO
*
      c(n)= c(n)/b(n)   
*
* back substitution ( makes c the sigma of text)
*
      DO k=n-1,1,-1
       c(k)= (c(k)-d(k)*c(k+1))/b(k)
      END DO
*
*.....Compute polynomial coefficients...................................
*
      b(n)= (y(n)-y(n-1))/d(n-1)+d(n-1)*(c(n-1)+c(n)+c(n))
      DO k=1,n-1
       b(k)= (y(k+1)-y(k))/d(k)-d(k)*(c(k+1)+c(k)+c(k))
       d(k)= (c(k+1)-c(k))/d(k)
       c(k)= THREE*c(k)
      END DO
      c(n)= THREE*c(n)
      d(n)= d(n-1)
*
      RETURN
*
      END Subroutine FMMsplineSingleYWW   
*
*------------------------------------------------------------------------
*
      SUBROUTINE Seval3SingleYWW(u,b,c,d,top,gdim,f,fp,fpp,fppp)
*
* ---------------------------------------------------------------------------
*
      REAL*8,INTENT(IN) :: u 
* abscissa at which the spline is to be evaluated
*
      INTEGER j,k,n,l,top,gdim
*
      REAL*8, dimension(109) :: xc,yc
      REAL*8, dimension(109) :: x,y
      REAL*8, DIMENSION(gdim) :: b,c,d 
* linear,quadratic,cubic coeff
*
      REAL*8,INTENT(OUT),OPTIONAL:: f,fp,fpp,fppp 
* function, 1st,2nd,3rd deriv
*
      INTEGER, SAVE :: i=1
      REAL*8    :: dx
      REAL*8,PARAMETER:: TWO=2.0, THREE=3.0, SIX=6.0
*
* The grid
*
      DATA (xc(l),l=1,109)/
     #   90.d0,95.d0,100.d0,105.d0,110.d0,115.d0,120.d0,
     #   125.d0,130.d0,135.d0,140.d0,145.d0,150.d0,155.d0,160.d0,165.d0,
     #   170.d0,175.d0,180.d0,185.d0,190.d0,195.d0,200.d0,210.d0,220.d0,
     #   230.d0,240.d0,250.d0,260.d0,270.d0,280.d0,290.d0,300.d0,310.d0,
     #   320.d0,330.d0,340.d0,350.d0,360.d0,370.d0,380.d0,390.d0,400.d0,
     #   410.d0,420.d0,430.d0,440.d0,450.d0,460.d0,470.d0,480.d0,490.d0,
     #   500.d0,510.d0,520.d0,530.d0,540.d0,550.d0,560.d0,570.d0,580.d0,
     #   590.d0,600.d0,610.d0,620.d0,630.d0,640.d0,650.d0,660.d0,670.d0,
     #   680.d0,690.d0,700.d0,710.d0,720.d0,730.d0,740.d0,750.d0,760.d0,
     #   770.d0,780.d0,790.d0,800.d0,810.d0,820.d0,830.d0,840.d0,850.d0,
     #   860.d0,870.d0,880.d0,890.d0,900.d0,910.d0,920.d0,930.d0,940.d0,
     #   950.d0,960.d0,970.d0,980.d0,990.d0,1000.d0,
     #   1.2d3,1.4d3,1.5d3,1.6d3,1.8d3,2.0d3/
*
      DATA (yc(l),l=1,109)/
     # 4.60d-6,1.10d-5,2.72d-5,6.36d-5,1.36d-4,2.68d-4,4.95d-4,8.73d-4,
     # 1.49d-3,2.47d-3,4.10d-3,6.86d-3,1.21d-2,2.41d-2,7.53d-2,2.36d-1,
     # 3.67d-1,4.80d-1,5.88d-1,7.02d-1,8.17d-1,9.36d-1,1.06d0,1.33d0,
     # 1.65d0,2.00d0,2.39d0,2.83d0,3.32d0,3.87d0,4.47d0,5.12d0,5.83d0,
     # 6.60d0,7.44d0,8.32d0,9.26d0,10.3d0,11.4d0,12.7d0,14.1d0,15.5d0,
     # 17.0d0,18.6d0,20.2d0,22.0d0,23.9d0,25.8d0,27.9d0,30.0d0,32.3d0,
     # 34.6d0,37.1d0,39.7d0,42.4d0,45.2d0,48.1d0,51.2d0,54.4d0,57.7d0,
     # 61.2d0,64.8d0,68.5d0,72.4d0,76.5d0,80.7d0,85.0d0,89.6d0,94.3d0,
     # 99.1d0,104.d0,109.d0,115.d0,120.d0,126.d0,132.d0,139.d0,145.d0,
     # 152.d0,159.d0,166.d0,173.d0,181.d0,189.d0,197.d0,206.d0,214.d0,
     # 223.d0,233.d0,242.d0,252.d0,263.d0,273.d0,284.d0,296.d0,308.d0,
     # 320.d0,332.d0,345.d0,359.d0,373.d0,387.d0,401.9d0,
     # 824.7d0,1623.7d0,2253.9d0,3108.6d0,5797.9d0,10517.2d0/
*
      n= 109
      FORALL(l=1:109)
       x(l)= xc(l)
       y(l)= yc(l)
      ENDFORALL
*
*.....First check if u is in the same interval found on the
*     last call to Seval.............................................
*
      IF (  (i<1) .OR. (i >= n) ) i=1
      IF ( (u < x(i))  .OR.  (u >= x(i+1)) ) THEN
       i=1   
*
* binary search
*
       j= n+1
       DO
        k= (i+j)/2
        IF (u < x(k)) THEN
         j= k
        ELSE
         i= k
        ENDIF
        IF (j <= i+1) EXIT
       ENDDO
      ENDIF
*
      dx= u-x(i)   
*
* evaluate the spline
*
      IF (Present(f))    f= y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
      IF (Present(fp))   fp= b(i)+dx*(TWO*c(i) + dx*THREE*d(i))
      IF (Present(fpp))  fpp= TWO*c(i) + dx*SIX*d(i)
      IF (Present(fppp)) fppp= SIX*d(i)
*
      RETURN
*
      END Subroutine Seval3SingleYWW  
*
      END SUBROUTINE gridYWW
*
*-------------------------------------------------------------------------------
*
