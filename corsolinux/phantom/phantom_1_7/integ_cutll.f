c
c Last update: Jun 25, 2007
c
c In case of e+e- collisions without ISR, the number of integration 
c variables (ndim) is 13, otherwise it is 15 as in the case of hadronic 
c colliders. Therefore in DO loops mxdim has been replaced by ndim 
c (I-search 'giuseppe ILC' to see changes).

************************************************************************
*                        SUBROUTINE INTEGRATION                        *
*                                                                      *
*                                                                      *
* Purpose:  integration run based on a modified multichannel approach. *
*         It determines the alfa(i) during thermalization.             *
*                                                                      *
*  Call to subroutine VEGAS                                            *
*  Fill phavegas0i.dat if iflat.eq.1 ( see WRITE(24,*) )               *
*                                                                      *
* 10/07/07  Small alfas (< 10**-3) are eliminated in third warm-up     *
*           iteration.                                                 *
*                                                                      *
* 21/05/07  All channels which survive after the third warm-up         *
*           iteration are processed with a larger number of calls for  *
*           the remaining iterations.                                  *
*                                                                      *
* variable acc for integration introduced depending on alfas           *
*  (see comment *varacc and variable accur                             *
*                                                                      *
************************************************************************


      SUBROUTINE integration


      IMPLICIT REAL*8 (a-b,d-h,o-z)


      INCLUDE 'common.h'
      INCLUDE 'common_subproc.h'
      PARAMETER (ndmx=50,mxdim=15)
c giuseppe 21/05/2007
      DIMENSION ncall_therm(2)
c end giuseppe 21/05/2007
      COMMON/iinteg_readinput/ncall_therm,itmx_therm,ncall,itmx
      COMMON/rinteg_readinput/acc_therm,acc
      COMMON/i_norm/i_normalize

      DIMENSION avgi(mxnphs),sd(mxnphs),y(ndmx,mxdim,mxnphs),
     &          ndo(mxnphs),si(mxnphs),swgt(mxnphs),schi(mxnphs)

*** variables for subroutine Vegas

      COMMON /pharand/ idum
      COMMON/abchia/rcalls
*6
      COMMON /phavestop/ivegasstop
      DIMENSION ytmp(ndmx,mxdim)
*6end

*** variables for phavegas/i/.dat, with i =sequential phase_space index 
      CHARACTER*80 phavegas
      character*1 charint1
      character*2 charint2

*** flag to use (1) or not (0) normalization in alfas determination 
      DATA inormal/1/

      EXTERNAL fxn

*** normalization : determines rnormal which is used in fxn
      IF (inormal.EQ.1) then

      PRINT*,' '
      PRINT*,'------------------------------------------------------'
      PRINT*,' '
      PRINT*,'NORMALIZATION'


* i_normalize is a flag that tells fxn not to compute amplitudes
        i_normalize=1

        itmx_norm=3
*        ncall_norm=200000
        ncall_norm=2000000
        acc_norm=0.003d0
        DO j=1,mxnphs
          rnormal(j)=1.d0
        ENDDO

        i=0
        DO j=1,mxnphs
          IF(ialfa(j).NE.0)THEN
            i=i+1
            init_norm=0
            iphs_ind=j
            print*,' '
            print*,'iphs_ind=',iphs_ind
            CALL vegas(region,ndim,fxn,init_norm,ncall_norm,itmx_norm,
     &             nprn,rnormal(iphs_ind),sd(i),rchi2a,acc_norm,
     &             y(1,1,i),it,ndo(i),si(i),swgt(i),schi(i))
          ENDIF
        ENDDO !j

        i_normalize=0

      ELSE

        DO j=1,mxnphs
          rnormal(j)=1.d0
        ENDDO

      ENDIF


*** initialization alfa(i)
      DO j=1,mxnphs
        IF(ialfa(j).NE.0)THEN
          alfa(j)=1.d0/n_kphs
c          print*,'iphs_ind and alfa=',j,alfa(j)
        ELSE
          alfa(j)=0.d0
        ENDIF
      ENDDO !j


      PRINT*,' '
      PRINT*,'------------------------------------------------------'
      PRINT*,' '
      PRINT*,'ALFA(i) DETERMINATION'
      PRINT*,' '

*** THERMALIZATION
*6
c      DO nit=1,itmx_therm

* ivegasstoptherm is the product of the ivegasstop of the various 
*   channels, so that all channel have eventually to stop at the same iteration

      ivegasstoptherm=0
      nit=1
      iweed=0

c giuseppe 21/05/2007
      ncall_ther=ncall_therm(1)
c end giuseppe 21/05/2007

      ichanged=0

      DO WHILE (nit.LE.itmx_therm.AND.ivegasstoptherm.EQ.0)
*6end

        IF(nit.eq.4.and.ichanged.eq.0)THEN
          ncall_ther=ncall_therm(2)
          itmx_therm=itmx_therm-3
          PRINT*,' '
          PRINT*,'----------------------------------------------------'
     &         //'------'
          PRINT*,' '
          PRINT*,'THERMALIZATION'
          PRINT*,' '
        ENDIF


*7
      avgi_tot=0.d0
      sd_tot=0.d0
*7end

        i=0
*6
        ivegasstoptherm=1        
*6end

        DO j=1,mxnphs
          IF(ialfa(j).NE.0)THEN
            i=i+1
            iphs_ind=j
            print*, '   '
            print*,'iphs_ind=',iphs_ind

            IF(nit.EQ.1.and.ichanged.eq.0)THEN
              init_therm=0
c giuseppe 21/05/2007
            ELSEIF(nit.EQ.4.and.ichanged.eq.0.or.nit.eq.1.
     &             and.ichanged.eq.1)THEN
              init_therm=1
              nit=1
              IF(ichanged.eq.0)THEN
              ENDIF
              ichanged=1
c end giuseppe 21/05/2007
            ELSE
              init_therm=2
            ENDIF
            it=nit
*6

            ivegasstop=0        ! flag which becomes =1 when vegas reach
                                !      the requested accuracy acc.    
*6end

c giuseppe 21/05/2007
c            CALL vegas(region,ndim,fxn,init_therm,ncall_therm,nit,
            CALL vegas(region,ndim,fxn,init_therm,ncall_ther,nit,
c end giuseppe 21/05/2007
     &                 nprn,avgi(i),sd(i),rchi2a,acc_therm,y(1,1,i),
     &                 it,ndo(i),si(i),swgt(i),schi(i))
*6
            IF(iweed.EQ.1) ivegasstop=0    ! if some phase space has been
                                           ! eliminated force all itmx_therm
                                           ! iterations
            ivegasstoptherm=ivegasstoptherm*ivegasstop
*6end

c results of various channels are put together only for the last 
c   iteration
*7
c            IF (nit.EQ.itmx_therm.OR.ivegasstoptherm.NE.0) THEN
c              avgi_tot=0.d0
c              sd_tot=0.d0
c            ENDIF
*7end
            IF (nit.EQ.itmx_therm.OR.ivegasstoptherm.NE.0) THEN
              avgi_tot=avgi_tot+avgi(i)
              sd_tot=sd_tot+sd(i)**2
            ENDIF
c
          ENDIF
        ENDDO !j
        
        app=0.d0
        i=0
        DO j=1,mxnphs
          IF(ialfa(j).ne.0)then
            i=i+1
            app=app+avgi(i)
          ENDIF
        ENDDO !j

        i=0
        DO j=1,mxnphs
          IF(ialfa(j).ne.0) then
            i=i+1
            iphs_ind=j
            alfa(iphs_ind)=avgi(i)/app
            print*,'iphs_ind and alfa=',iphs_ind,alfa(j)
          ENDIF
        ENDDO !j
*6

c Eliminate alfa's which are smaller than alfamin=1.d-3
c Performed only after the third iteration
c If some phase space is eliminated two more iteration are forced
        alfamin=1.d-3
        IF (nit.EQ.3.and.ichanged.eq.0) THEN
          app=0.d0
          i=0
          k=0
          DO j=1,mxnphs
            IF(ialfa(j).NE.0)THEN
              i=i+1
              IF(alfa(j).LT.alfamin) THEN
                ialfa(j)=0
              ELSE         
                k=k+1
                app=app+alfa(j)
                IF(k.NE.i)THEN  ! shift integration data
                  avgi(k)=avgi(i)
                  sd(k)=sd(i)
                  ndo(k)=ndo(i)
                  si(k)=si(i)
                  swgt(k)=swgt(i)
                  schi(k)=schi(i)
c giuseppe ILC 25/06/2007
c                  DO n2=1,mxdim
                  DO n2=1,ndim
c end giuseppe ILC 25/06/2007
                    DO n1=1,ndmx
                      y(n1,n2,k)=y(n1,n2,i)
                    ENDDO
                  ENDDO
                ENDIF ! IF(k.NE.i)
              ENDIF
            ENDIF
          ENDDO !j
c Renormalize
          IF(k.NE.i)THEN
            iweed=1
            DO j=1,mxnphs
              IF(ialfa(j).NE.0) THEN
                alfa(j)=alfa(j)/app
                print*,'NEW iphs_ind and alfa=',j,alfa(j)
              ENDIF
            ENDDO
          ENDIF
        ENDIF ! End check on small alfa's
          
        nit=nit+1
*6end

      ENDDO !nit
      sd_tot=sqrt(sd_tot)
      PRINT*,' '

      PRINT 200
  200 FORMAT('------------------------------------------------------')
      PRINT*,' '
      PRINT 201,avgi_tot,sd_tot
  201 FORMAT(' Sigma = ',d13.7,' +/-',d10.3,' (pb)')

      do j=1,mxnphs
        if(ialfa(j).ne.0)then
          print*,'iphs_ind=',j,'  alfa=',alfa(j)
        endif
      enddo

*** integration

      PRINT*,' '
      PRINT*,'------------------------------------------------------'
      PRINT*,' '
      PRINT*,'INTEGRATION'
      PRINT*,' '

      avgi_tot=0.d0
      sd_tot=0.d0

      i=0
      DO j=1,mxnphs
        IF(ialfa(j).NE.0)THEN
          i=i+1
          init=1
          iphs_ind=j
          print*,'iphs_ind=',iphs_ind

*varacc
          if (alfa(j).gt.0.1d0) then
            accur=acc
          else
            factacc=(.1d0/alfa(j))**(2.d0/3.d0)
            accur=acc*factacc
          endif
          print*,'alfa,accur',alfa(j),accur
*varaccend

          IF(iflat.EQ.0)THEN
            CALL vegas(region,ndim,fxn,init,ncall,itmx,nprn,avgi(i),
*varacc     &               sd(i),rchi2a,acc,y(1,1,i),it,ndo(i),si(i),
     &               sd(i),rchi2a,accur,y(1,1,i),it,ndo(i),si(i),
*varaccend
     &               swgt(i),schi(i))
          ELSE IF(iflat.EQ.1)THEN

* initialization of the variables used if IFLAT=1

            ivegasstop=0       ! flag which becomes =1 when vegas reach
                                !      the requested accuracy acc.    
            nevent=0            ! number of generated events
            novermax=0          ! number of points exceeding the maximum
            rmaxfxn=0.d0        ! maximum value of fxn*wgt after last 
                                !  but one Vegas iteration 
            rmaxfxnnew=0        ! maximum value of fxn*wgt after last
                                 ! Vegas iteration

            nit=1

            DO WHILE (nit.LE.itmx.AND.ivegasstop.EQ.0)

              IF(nit.GE.2)init=2

              CALL vegas(region,ndim,fxn,init,ncall,nit,nprn,avgi(i),
*varacc     &             sd(i),rchi2a,acc,y(1,1,i),it,ndo(i),si(i),
     &             sd(i),rchi2a,accur,y(1,1,i),it,ndo(i),si(i),
*varaccend
     &             swgt(i),schi(i))

* ivegasstop in phavegas can exit with value 1 only after the first iteration
              IF (ivegasstop.EQ.0.AND.nit.NE.itmx) THEN
                rmaxfxn=rmaxfxnnew
                rmaxfxnnew=0
                novermax=0
                nevent=0
              ENDIF

              IF(istorvegas.EQ.1)THEN
                avgitmp=avgi(i)
                sdtmp=sd(i)
                rchi2atmp=rchi2a
                DO n=1,ndmx
c giuseppe ILC 25/06/2007
c                  DO m=1,mxdim
                  DO m=1,ndim
c end giuseppe ILC 25/06/2007
                    ytmp(n,m)=y(n,m,i)
                  ENDDO         !m
                ENDDO           !n
                ndotmp=ndo(i)
                sitmp=si(i)
                swgttmp=swgt(i)
                schitmp=schi(i)
              ENDIF

              nit=nit+1

            ENDDO  !nit and ivegasstop

            IF(istorvegas.EQ.1)THEN
              if (i.le.9) then
                write (charint1,'(i1)') i
                charint2='0'//charint1
              elseif (i.le.99) then
                write (charint2,'(i2)') i
              endif
              phavegas='phavegas'//charint2//'.dat'

              OPEN(unit=24,file=phavegas,status='new')
              
              do k=1,8
                WRITE(24,*)iproc(k)
              enddo  
              WRITE(24,*) igroup,n_kphs,iphs_ind,nidenticalinitial
              do jj=1,mxnphs
                WRITE(24,*)ialfa(jj),rnormal(jj),alfa(jj)
              enddo
              do k=1,8
                WRITE(24,*)idp(k),idp_inout(k)
              enddo
              do jj=1,mxnphs
                do k=1,8
                  WRITE(24,*)iorder(k,jj)
                enddo
              enddo
              WRITE(24,*)avgitmp,sdtmp,rchi2atmp
              DO n=1,ndmx
c giuseppe ILC 25/06/2007
c                DO m=1,mxdim
                DO m=1,ndim
c end giuseppe ILC 25/06/2007
                  WRITE(24,*)ytmp(n,m)
                ENDDO           !m
              ENDDO             !n
              WRITE(24,*)ndotmp,sitmp,swgttmp,schitmp

              WRITE(24,*)rmaxfxnnew,rcalls
              WRITE(24,*)avgi(i)
              WRITE(24,*)ncall,itmx
              CLOSE(24)
            ENDIF

*6end
          ENDIF       !iflat

          avgi_tot=avgi_tot+avgi(i)
          sd_tot=sd_tot+sd(i)**2

          PRINT 202
 202    FORMAT('------------------------------------------------------')
          PRINT*,' '
          PRINT 203,avgi(i),sd(i)
 203      FORMAT(' Sigma = ',d13.7,' +/-',d10.3,' (pb)')

          IF(iflat.EQ.1)THEN
            PRINT*,'Informations about flat events generation:'
            PRINT*,'         ----------------------           '

           PRINT 204,rmaxfxn
 204       FORMAT(' Maximum after next-to-last VEGAS iteration = ',d9.3)
*6
c            PRINT 205,rmaxfxn_2
            PRINT 205,rmaxfxnnew
*6end
 205        FORMAT(' Maximum after last VEGAS iteration = ',d9.3)
            PRINT 206,nevent
 206        FORMAT(' Flat events number = ',i9)
            PRINT 207,novermax
 207        FORMAT(' number of function values over maximum = ',i9)
          ENDIF

        ENDIF ! (ialfa(j).NE.0)

      ENDDO !j

      sd_tot=sqrt(sd_tot)

      PRINT 208,avgi_tot,sd_tot
 208  FORMAT(' Sigma = ',d13.7,' +/-',d10.3,' (pb)')

      RETURN
      END
