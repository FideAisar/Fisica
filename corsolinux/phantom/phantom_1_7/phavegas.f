C
C OpenMP  MXDIM=15 
C
C Must be compiled with the additional flag -omp on True64
C in order to use OpenMP features.
C The number of threads/processors is specified in the calling c-shell
C script (The input script) with:
C setenv OMP_NUM_THREADS n where n is the number of cpu's to be used
C
 
      SUBROUTINE vegas(region,ndim,fxn,init,ncall,itmx,nprn,tgral,sd,
     *  chi2a,acc,xi,it,ndo,si,swgt,schi)
      IMPLICIT NONE
Cc ho aggiunto acc
      INTEGER init,itmx,ncall,ndim,nprn,ndmx,mxdim,ncall_eff
      REAL*8 tgral,chi2a,sd,acc,region(2*ndim),fxn,alph,tiny
      PARAMETER (alph=1.5,ndmx=50,mxdim=15,tiny=1.d-30)
C      PARAMETER (ALPH=1.5,NDMX=100,MXDIM=10,TINY=1.e-30)
      EXTERNAL fxn
CU    USES fxn,ran2,rebin
      INTEGER i,idum,it,j,k,mds,nd,ndo,ng,npg,ia(mxdim),kg(mxdim)
      REAL*8 rcalls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo,
     *d(ndmx,mxdim),di(ndmx,mxdim),dt(mxdim),dx(mxdim),r(ndmx),x(mxdim),
     *xi(ndmx,mxdim),xin(ndmx),ran2
      REAL*8 schi,si,swgt,resl,standdevl
*6
c      COMMON/abresl/resl(10),standdevl(10)
*6end
      COMMON/abchia/rcalls
      COMMON/abstat/ncall_eff
      INTEGER nevent,nflevts
      COMMON/phaiflav/nevent,nflevts

*for PHASE
      INTEGER ionesh
      COMMON/phaones/ionesh
      COMMON /pharand/ idum
*6
      INTEGER ivegasstop
      COMMON /phavestop/ivegasstop
*6end
*
      INTEGER ntotcells,naux,nlabel,nlocal
      REAL*8 sigma_sq_inv
      REAL*8 ran_numb(mxdim)

      DATA mds/1/           ! mds =  1 importance sampling only
                            ! mds = -1 importance + stratified sampling
      SAVE
  
      IF(init.LE.0)THEN     ! Normal entry. Start from scratch
        mds=1
        ndo=1
C**
        it=1
C**
        DO 11 j=1,ndim
          xi(1,j)=1.d0     ! Initialization of variables passed to rebin
   11   CONTINUE            ! The actual values are totally irrelevant
      ENDIF   ! End IF (init.LE.0)
      
      
      IF (init.LE.1)THEN    ! Enter here to inherit grid but not results
        si=0.d0
        swgt=0.d0
        schi=0.d0
C**
        it=1
C**
      ENDIF    ! End IF (init.LE.1)

      IF (init.LE.2)THEN    ! Enter here to inherit grid and results
        nd=ndmx             ! nd = number of elementary cell in the grid
        ! they are grouped ndmx/ng at a time for generating points
        ! but they are kept separated for grid refinements (di(i,j) and 
                     ! possibly d(i,j) run on nd values per dimension) 
        ng=1         ! ng = number of supercells in each dimension
        IF(mds.NE.0)THEN
          ng=(ncall/2.d0+0.25d0)**(1.d0/ndim)
          mds=1
          IF((2*ng-ndmx).ge.0)then
            mds=-1
            naux=ng/ndmx+1
            nd=ng/naux
            ng=naux*nd
          ENDIF
        ENDIF
        ntotcells=ng**ndim            ! total number of supercells
        npg=max(ncall/ntotcells,2)    ! number of points per supercell
        rcalls=npg*ntotcells          
                       ! total number of effective calls per iteration
        dxg=1.d0/ng
        dv2g=(rcalls*dxg**ndim)**2/npg/npg/(npg-1.d0)
        xnd=nd
        dxg=dxg*xnd                   ! number of cells per supercell
        xjac=1.d0/rcalls

        DO 12 j=1,ndim
          dx(j)=region(j+ndim)-region(j)
          xjac=xjac*dx(j)
   12   CONTINUE

        IF(nd.NE.ndo)THEN   !  Do binning if necessary
c protect ionesh from redefining the grid

          if (ionesh.eq.1) then
            print*, '   '
            print*,' *************ERROR***********'
            print*, 'nd.ne ndo in wvegas!'
            stop
          endif
          

          DO 13 i=1,nd
            r(i)=1.d0
   13     CONTINUE
          DO 14 j=1,ndim
            CALL rebin(ndo/xnd,nd,r,xin,xi(1,j))     
                              ! here all cells are created equal
   14     CONTINUE
          ndo=nd
        ENDIF
        IF(nprn.GE.0) WRITE(*,200) ndim,rcalls,it,itmx
      ENDIF    ! ENF IF (init.LE.2)
      
      DO 28 it=it,itmx      ! Main iteration loop
*6
c        IF(it.GE.2.AND.acc*abs(tgral).ge.sd) RETURN
*6end
        ti=0.d0
        tsi=0.d0
        DO 16 j=1,ndim
          kg(j)=1
          DO 15 i=1,nd
            d(i,j)=0.d0
            di(i,j)=0.d0
   15     CONTINUE
   16   CONTINUE

c for oneshot a cell is chosen at random and 
c   the function is evaluated only once

       IF (ionesh.EQ.1) THEN
          wgt=xjac
          DO  j=1,ndim
            kg(j)=ran2(idum)*ng+1

            xn=(kg(j)-ran2(idum))*dxg+1.d0
            ia(j)=max(min(int(xn),ndmx),1)
            IF(ia(j).GT.1)THEN     ! xi(i,j)= upper limit of i-th cell 
                                        !  in j-th direction
              xo=xi(ia(j),j)-xi(ia(j)-1,j)
              rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
            ELSE                        ! there is no xi(0,j)
              xo=xi(ia(j),j)
              rc=(xn-ia(j))*xo
            ENDIF
            x(j)=region(j)+rc*dx(j)
            wgt=wgt*xo*xnd
          END DO
          f=wgt*fxn(x,wgt)
          RETURN
        ENDIF

c normal vegas continue: loop on supercells

* Return point for loop on cells 
C   10   CONTINUE
c The inner loop has been made explicit and exposed to parallelization
c
c Given ndim dimension and ng cells per dimension we have ng**ndim cells
c To identify a cell we can take a number ncell between 0 and ng**ndim-1
c and interpret it as a number in base ng with ndim digits.
c Taking ndim times the mod[ng,ncell] and then integer-dividing by ng
c (and adding one to all mod's) one gets an ndim unique vector with 
c entries between 1 and ng.
c 
!    $OMP PARALLEL DEFAULT(SHARED)
!    $OMP+         PRIVATE(NLOCAL,KG,FB,F2B,WGT,XN,IA,XO,RC,X,F,F2)
!$OMP PARALLEL DEFAULT(PRIVATE)
!$OMP+ SHARED(NG,NDIM,NPG,XJAC,IDUM,DXG,XI,REGION,XND,DX,DI,MDS,D,
!$OMP+        TI,TSI)

!   $OMP DO REDUCTION(+:TI,TSI)
!$OMP DO
        DO nlabel=0,ng**ndim-1
          nlocal=nlabel
          DO k=1,ndim
            kg(k)=mod(nlocal,ng)+1
            nlocal=nlocal/ng
          ENDDO   
          fb=0.d0
          f2b=0.d0
          DO 19 k=1,npg      ! Loop on points in a supercell
            wgt=xjac
!$OMP CRITICAL(idum_updating)
            DO j=1,ndim
              ran_numb(j)=ran2(idum)        ! Flat random variables
            ENDDO
!$OMP END CRITICAL(idum_updating)
            DO 17 j=1,ndim
             xn=(kg(j)-ran_numb(j))*dxg+1.d0  
                            ! Uniform coordinate in a supercell
              ia(j)=max(min(int(xn),ndmx),1)  
                           ! Uniform cell index: integer part of xn
              IF(ia(j).GT.1)THEN
                xo=xi(ia(j),j)-xi(ia(j)-1,j)  
                         ! width of the ia(j)-th cell in the j direction
                rc=xi(ia(j)-1,j)+(xn-ia(j))*xo 
                               ! upper of (ia(j)-1)-th cell + 
                               ! uniformdistribution  inside the cell
              ELSE
                xo=xi(ia(j),j)
                rc=(xn-ia(j))*xo
              ENDIF
* EM 10/3/03 additional check to avoid x(j)>1.d0
              rc=min(rc,1.d0)
* EMend
              x(j)=region(j)+rc*dx(j)
              wgt=wgt*xo*xnd    ! wgt includes width of cell dxg=xnd/ng
                                ! ng**(-ndmx) is in xjac 
   17       CONTINUE
            f=wgt*fxn(x,wgt)
            f2=f*f
            fb=fb+f
            f2b=f2b+f2
!$OMP CRITICAL(point_downloading)
            DO 18 j=1,ndim
              di(ia(j),j)=di(ia(j),j)+f
              IF(mds.GE.0) d(ia(j),j)=d(ia(j),j)+f2
   18       CONTINUE
!$OMP END CRITICAL(point_downloading)
   19     CONTINUE             ! end do on k=1,npg
          f2b=sqrt(f2b*npg)
          f2b=(f2b-fb)*(f2b+fb)
          IF (f2b.LE.0.d0) f2b=tiny
!$OMP CRITICAL(cell_result_downloading)
          ti=ti+fb            ! ti accumulates fxn*wgt
          tsi=tsi+f2b    
          ! tsi accumulates   sum(npg*(fxn*wgt)**2) - (sum(f*wgt))**2
       ! the extra factor of npg is later eliminated when multiplying by
       ! dv2g leading to the correct averages
       ! NOTICE THAT THE COMPUTATION OF ti AND tsi IS THE
       ! ONLY LINK BETWEEN THE DIFFERENT CYCLES IN THE LOOP
       ! IF ONE USES THE NEW METHOD OF DETERMINING kg(k)
          IF(mds.LT.0)THEN    ! Use stratified sampling
            DO 21 j=1,ndim
              d(ia(j),j)=d(ia(j),j)+f2b
   21       CONTINUE
          ENDIF
!$OMP END CRITICAL(cell_result_downloading)
       ENDDO   ! end do on nlabel
!$OMP END DO       
!$OMP END PARALLEL

! Compute final result for this iteration: ti already contains the 
!    average value of f*wgt
        tsi=tsi*dv2g        
        sigma_sq_inv=1.d0/tsi  ! 1/sigma(itmx)**2
        si=si+sigma_sq_inv*ti  ! Accumulates  ti(itmx)/sigma(itmx)**2
        schi=schi+sigma_sq_inv*ti**2    
                               ! Accumulates  ti(itmx)**2/sigma(itmx)**2
        swgt=swgt+sigma_sq_inv ! Accumulates  1/sigma(itmx)**2
        tgral=si/swgt
        chi2a=max((schi-si*tgral)/(it-.99d0),0.d0)
        sd=sqrt(1.d0/swgt)     ! sd= 1/sqrt(sum of 1/sigma(itmx)**2)
        tsi=sqrt(tsi)          ! sd for this iteration
        IF(nprn.GE.0)THEN
C**  aggiunta di ncall_eff e sua inizializzazione dopo ogni iterazione
C**  ho modificato anche il FORMAT 201
          WRITE(*,201) it,ncall_eff,it,ti,tsi,tgral,sd,chi2a
**for pgf to write immeditely output
*          call flush(6)

          ncall_eff=0
*6
c          resl(it)=ti
c          standdevl(it)=tsi
*6end
          IF(nprn.NE.0)THEN
            DO 23 j=1,ndim
              WRITE(*,202) j,(xi(i,j),di(i,j),i=1+nprn/2,nd,nprn)
   23       CONTINUE
          ENDIF
        ENDIF

        DO 25 j=1,ndim      ! refine the grid
          xo=d(1,j)
          xn=d(2,j)
          d(1,j)=(xo+xn)/2.d0  
             ! average of (f*wgt)**2 aver the first two cells --> d(1,j)
          dt(j)=d(1,j)
          DO 24 i=2,nd-1
            rc=xo+xn
            xo=xn
            xn=d(i+1,j)
            d(i,j)=(rc+xn)/3.d0   
           ! average of (f*wgt)**2 aver the cells (i-1,i,i+1) --> d(i,j)
            dt(j)=dt(j)+d(i,j)    
       ! dt sums over the averaged cells. Notice that the sums does not 
       ! coincide with the sum over the cells before averaging
   24     CONTINUE
          d(nd,j)=(xo+xn)/2.d0    
       ! average of (f*wgt)**2 aver the last two cells --> d(nd,j)
          dt(j)=dt(j)+d(nd,j)     
   25   CONTINUE
        DO 27 j=1,ndim
          rc=0.d0
          DO 26 i=1,nd
            IF(d(i,j).lt.tiny) d(i,j)=tiny
            r(i)=((1.-d(i,j)/dt(j))/(log(dt(j))-log(d(i,j))))**alph    
                                           ! black magic
            rc=rc+r(i)
   26     CONTINUE
          CALL rebin(rc/xnd,nd,r,xin,xi(1,j))
   27   CONTINUE
*6
        IF(it.GE.2.AND.acc*abs(tgral).ge.sd) THEN
          ivegasstop=1
          RETURN
        ENDIF
*6end
   28 CONTINUE
      RETURN
  200 FORMAT(/' input parameters for vegas:  ndim=',i3,'  ncall=',
     *f12.0/28x,'  it=',i5,'  itmx=',i5)
  201 FORMAT(/' iteration no.',i3,':',12x,'effective ncall=',i11/
     *' iteration no.',i3,': ','integral =',g14.7,'+/- ',g9.2/
     *' all iterations:   integral =',g14.7,'+/- ',g9.3,' chi**2/it'
     *'n =',g9.2)
  202 FORMAT(/' data for axis ',i2/'    X       delta i       ',
     *'   x       delta i       ','    x       delta i       ',/(1x,
     *f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #>,1')5c).

      FUNCTION ran2(idum)
      INTEGER idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      REAL*8 ran2,am,eps,rnmx
      PARAMETER (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     *ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,
     *ntab=32,ndiv=1+imm1/ntab,eps=1.2d-7,rnmx=1.-eps)
      INTEGER idum2,j,k,iv(ntab),iy
      DATA idum2/123456789/, iv/ntab*0/, iy/0/
      IF (idum.LE.0) THEN
        idum=max(-idum,1)
        idum2=idum
        DO 11 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          IF (idum.LT.0) idum=idum+im1
          IF (j.LE.ntab) iv(j)=idum
   11   CONTINUE
        iy=iv(1)
      ENDIF
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      IF (idum.LT.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      IF (idum2.LT.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      IF(iy.LT.1)iy=iy+imm1
      ran2=min(am*iy,rnmx)
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #>,1')5c).

      SUBROUTINE rebin(rc,nd,r,xin,xi)
      INTEGER nd
      REAL*8 rc,r(*),xi(*),xin(*)
C xin is purely internal to rebin. It enters the call to rebin in order to
C be an assumed dimension array
C rc is the sum of all bin height divided by the numer of bins = content of one new bin
C r is the distribution which has to be flattened
C xi is the array of the upper bounds of each cell which is recreated anew each time 
      INTEGER i,k
      REAL*8 dr,xn,xo
      k=0
      xn=0.d0
      dr=0.d0
      DO 11 i=1,nd-1
    1   IF(rc.GT.dr)THEN
          k=k+1
          dr=dr+r(k)
          xo=xn
          xn=xi(k)
          GOTO 1
        ENDIF
        dr=dr-rc
        xin(i)=xn-(xn-xo)*dr/r(k)
   11 CONTINUE
      DO 12 i=1,nd-1
        xi(i)=xin(i)
   12 CONTINUE
      xi(nd)=1.d0
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #>,1')5c).



* EM additional random number routines are needed if random numbers
* have to be generated in the full program OUTSIDE of VEGAS ( e.g for the
* overall azimuthal angle of an event) in order to avoid conflict between
* threads when using more than one processor

      SUBROUTINE random (r)
      IMPLICIT NONE
      REAL*8 r
      INTEGER m, a, c
      PARAMETER (M = 259200, A = 7141, C = 54773)
      INTEGER n
      SAVE n
      DATA n /0/
      n = mod(n*a+c,m)
      r = dble (n) / dble (m)
      END


      FUNCTION rnd_1()
      IMPLICIT NONE
      REAL*8 rnd_1
      INTEGER m, a, c
      PARAMETER (M = 259200, A = 7141, C = 54773)
      INTEGER n
      SAVE n
      DATA n /0/
      n = mod(n*a+c,m)
      rnd_1 = dble (n) / dble (m)
      END
