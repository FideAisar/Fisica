************************************************************************
* oneshot.f
*
* Last update: Jun 25, 2007
*
c In case of e+e- collisions without ISR, the number of integration 
c variables is 13, otherwise it is 15 as for the hadronic colliders.
c Therefore in DO loops mxdim has been replaced by ndim.
c (I-search 'giuseppe ILC' to see the changes).
************************************************************************

      subroutine oneshot

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT COMPLEX*16 (c)

c giuseppe 08/09/2006
c      CHARACTER*80 charaux
      CHARACTER*200 charaux
c end giuseppe 08/09/2006
      INCLUDE 'common.h'
      INCLUDE 'common_subproc.h'
      PARAMETER (ndmx=50,mxdim=15)
      PARAMETER (nmaxfiles=500)
      COMMON/phaones/ionesh

      COMMON /pharand/ idum

      dimension igroupos(nmaxfiles),n_kphsos(nmaxfiles)
      integer iphs_indos(nmaxfiles),nidenticalinitialos(nmaxfiles),
     &     ialfaos(mxnphs,nmaxfiles),idpos(8,nmaxfiles),
     &     idp_inoutos(8,nmaxfiles),iorderos(8,mxnphs,nmaxfiles)
*sandro 25/07
     &     ,iprocos(8,nmaxfiles)
*sandro 25/07 end
      dimension rnormalos(mxnphs,nmaxfiles),alfaos(mxnphs,nmaxfiles)

      dimension avgi1os(nmaxfiles),sd1os(nmaxfiles),rchi2aos(nmaxfiles),
     &     y1os(ndmx,mxdim,nmaxfiles)
      dimension ndo1os(nmaxfiles)
      dimension si1os(nmaxfiles),swgt1os(nmaxfiles),schi1os(nmaxfiles),
     &     rmaxfxnos(nmaxfiles), rcallos(nmaxfiles),
     &     avgi_totos(nmaxfiles)

      dimension ncallos(nmaxfiles),itmxos(nmaxfiles)
      dimension prob(nmaxfiles)

c      dimension y1(ndmx,mxdim)

c rnacc counts the number of generated events for every input file,   
c rntotproc the total number of events for every input file
  
      dimension rnacc(nmaxfiles), rntotproc(nmaxfiles)
      dimension noverm (nmaxfiles)

      EXTERNAL fxn


      IF(iwrite_event.EQ.1)THEN
        OPEN(unit=23,file='phamom.dat',status='new')
        WRITE(23,'(A)')'<LesHouchesEvents version="1.0">'
        CALL LHAFileInit      
      ENDIF


      call ireadnoclose('nfiles',nfiles,1)

      do i=1,nfiles
	  
        read (30,'(a)') charaux

        print *, charaux

        open (unit=24,file=charaux, status='old')
*sandro 25/07
        do k=1,8
          READ(24,*)iprocos(k,i)
        enddo
*sandro 25/07 end
        READ(24,*) igroupos(i),n_kphsos(i),iphs_indos(i),
     &       nidenticalinitialos(i)
        do j=1,mxnphs
          READ(24,*)ialfaos(j,i),rnormalos(j,i),alfaos(j,i)
        enddo
        do k=1,8
          READ(24,*)idpos(k,i),idp_inoutos(k,i)
        enddo
        do j=1,mxnphs
          do k=1,8
            READ(24,*)iorderos(k,j,i)
          enddo
        enddo
        READ(24,*)avgi1os(i),sd1os(i),rchi2aos(i)
        DO k=1,ndmx
c giuseppe ILC 25/06/2007
c          DO j=1,mxdim
          DO j=1,ndim
c end giuseppe ILC 25/06/2007
            READ(24,*)y1os(k,j,i)
          ENDDO  !j
        ENDDO    !k
        READ(24,*)ndo1os(i),si1os(i),swgt1os(i),schi1os(i)
        READ(24,*)rmaxfxnos(i), rcallos(i)
        READ(24,*)avgi_totos(i)
        READ(24,*)ncallos(i),itmxos(i)
        close (unit=24)
      enddo

      close (unit=30)

c  determine the probabilities of the various processes
c   prob(i)= Sum_1^i prob(i)  
      rmaxtot=0.d0

      IF(i_coll.eq.1)THEN    ! p-p collider

        if (i_exchincoming.eq.0) then
          do i=1,nfiles
            rmaxtot=rmaxtot+rmaxfxnos(i)*rcallos(i)  
          enddo
        else if (i_exchincoming.eq.1) then
          do i=1,nfiles
            rmaxtot=rmaxtot+
     &      (1+mod(nidenticalinitialos(i),2))*rmaxfxnos(i)*rcallos(i)
          enddo
        else
          print *, ' ERROR: value of i_exchincoming not valid!!! '
          stop
        endif
        if (i_exchincoming.eq.0) then
          prob(1)=rmaxfxnos(1)*rcallos(1)/rmaxtot
          do i=2,nfiles
            prob(i)=rmaxfxnos(i)*rcallos(i)/rmaxtot+prob(i-1)
          enddo
        else
          prob(1)=(1+mod(nidenticalinitialos(1),2))*rmaxfxnos(1)*
     &       rcallos(1)/rmaxtot
          do i=2,nfiles
            prob(i)=(1+mod(nidenticalinitialos(i),2))*rmaxfxnos(i)*
     &         rcallos(i)/rmaxtot+prob(i-1)
          enddo
        endif
      
      ELSEIF(i_coll.eq.2)THEN    ! p-pbar collider

        do i=1,nfiles
          rmaxtot=rmaxtot+rmaxfxnos(i)*rcallos(i)  
        enddo
        prob(1)=rmaxfxnos(1)*rcallos(1)/rmaxtot
        do i=2,nfiles
          prob(i)=rmaxfxnos(i)*rcallos(i)/rmaxtot+prob(i-1)
        enddo

      ELSEIF(i_coll.eq.3)THEN   ! e+e- collider

        do i=1,nfiles
          rmaxtot=rmaxtot+rmaxfxnos(i)*rcallos(i)  
        enddo
        prob(1)=rmaxfxnos(1)*rcallos(1)/rmaxtot
        do i=2,nfiles
          prob(i)=rmaxfxnos(i)*rcallos(i)/rmaxtot+prob(i-1)
        enddo
      
      ENDIF                     ! i_coll

c  now start unweighted sampling 

      nevent=0

      do while (nevent.lt.nunwevts)

        neventstart=nevent
        novermaxstart=novermax
             
c     choose the process

        x=ran2(idum)
        i=1
        do while (x.gt.prob(i))
          i=i+1
        enddo
        

*once the process has been chosen

*sandro 25/07

c define the proc and call proc (need iproc also for the part with
c    ionesh=1

        do k=1,8
          iproc(k)=iprocos(k,i)
        enddo
*test proc
c        CALL proc ! proc  determines the variables for every phase space
c                 ! that can be used and the variables for the process
c                 ! at hand
c                 ! If ioneshot=1 the latter are read from phavegas.dat
c                 !  and only the former are determined in proc
        CALL proc ! proc  determines the variables for every phase space
* test proc end
*sandro 25/07 end
        CALL procextraini

c     define the input variables
          
        igroup=igroupos(i)
        n_kphs=n_kphsos(i)
        iphs_ind=iphs_indos(i)
        nidenticalinitial=nidenticalinitialos(i)
        do j=1,mxnphs
          ialfa(j)=ialfaos(j,i)
          rnormal(j)=rnormalos(j,i)
          alfa(j)=alfaos(j,i)
        enddo
        do k=1,8
          idp(k)=idpos(k,i)
          idp_inout(k)=idp_inoutos(k,i)
          do j=1,mxnphs
            iorder(k,j)=iorderos(k,j,i)
          enddo
        enddo
        

*eventually choose at random the flavour of leptons 
*  - up to now only for  one couple
        if (i_emutau.ge.1) then
          kk=0
          do k=1,8
            indabs=abs(idp(k))
            if (indabs.gt.10.and.indabs.lt.17) then
              kk=kk+1
              if (kk.eq.1) k1=k
              if (kk.eq.2) k2=k
            endif
          enddo
          if (kk.eq.2) then
            kk1=abs(idp(k1))
            kk2=abs(idp(k2))
            if (abs(kk1-kk2).gt.1) then 
              print*, 'difference between kk1 and kk2 not allowed'
              stop
            endif

            if (i_emutau.eq.1) then 
              irand=int(ran2(idum)*2)+1
            elseif (i_emutau.eq.2) then
              irand=int(ran2(idum)*3)+1
            else
              print*, 'value of i_emutau not allowed'
              stop
            endif
              
            if (kk1.eq.11.or.kk1.eq.13.or.kk1.eq.15) then
              if (irand.eq.1) then
                idp(k1)=sign(11,idp(k1))
                if(kk1.eq.kk2) then
                  idp(k2)=sign(11,idp(k2))
                else
                  idp(k2)=sign(12,idp(k2))
                endif
              elseif (irand.eq.2) then
                idp(k1)=sign(13,idp(k1))
                if(kk1.eq.kk2) then
                  idp(k2)=sign(13,idp(k2))
                else
                  idp(k2)=sign(14,idp(k2))
                endif
              elseif (irand.eq.3.and.i_emutau.eq.2) then
                idp(k1)=sign(15,idp(k1))
                if(kk1.eq.kk2) then
                  idp(k2)=sign(15,idp(k2))
                else
                  idp(k2)=sign(16,idp(k2))
                endif
              else
                print*, 'irand not allowed'
                stop
              endif
            else if (kk1.eq.12.or.kk1.eq.14.or.kk1.eq.16) then
              if (irand.eq.1) then
                idp(k1)=sign(12,idp(k1))
                if(kk1.eq.kk2) then
                  idp(k2)=sign(12,idp(k2))
                else
                  idp(k2)=sign(11,idp(k2))
                endif
              elseif (irand.eq.2) then
                idp(k1)=sign(14,idp(k1))
                if(kk1.eq.kk2) then
                  idp(k2)=sign(14,idp(k2))
                else
                  idp(k2)=sign(13,idp(k2))
                endif
              elseif (irand.eq.3.and.i_emutau.eq.2) then
                idp(k1)=sign(16,idp(k1))
                if(kk1.eq.kk2) then
                  idp(k2)=sign(16,idp(k2))
                else
                  idp(k2)=sign(15,idp(k2))
                endif
              else
                print*, 'irand not allowed'
                stop
              endif
            endif
          else
            print*, 'i_emutau.gt.1 implemented only for 2 leptons'
            stop
          endif
        endif


c   define vegas variables
        ncall=ncallos(i)
        itmx=itmxos(i)
        avgi1=avgi1os(i)
        sd1=sd1os(i)
        rchi2a=rchi2aos(i)
c        DO k=1,ndmx
c          DO j=1,mxdim
c            y1(k,j)=y1os(k,j,i)
c          ENDDO                 !j
c        ENDDO                   !k
        ndo1=ndo1os(i)
        si1=si1os(i)
        swgt1=swgt1os(i)
        schi1=schi1os(i)
        
        rmaxfxn=rmaxfxnos(i)
        
        rcalls=rcallos(i)
        
          

        rntotproc(i)=rntotproc(i)+1.d0
          
        CALL extrema
          
c  integration
* ricordare che per quanto scritto in fxn nel generare eventi 
* it deve essere uguale itmx and init g.e.1  
        nit=itmx
        init=2
        nprn=-1
        it1=itmx

        CALL vegas(region,ndim,fxn,init,ncall,nit,nprn,avgi1,
     &       sd1,rchi2a,acc,y1os(1,1,i),it1,ndo1,si1,swgt1,schi1)

        if (nevent.gt.neventstart) then
          rnacc(i)=rnacc(i)+1.d0
          rnacctot=rnacctot+1.d0
      if (mod(nevent,100).eq.0) print *,'Event =',nevent
        endif         
        
        if (novermax.gt.novermaxstart) then
          noverm(i)=noverm(i)+1
        endif
        
        rntot=rntot+1.d0
        
      enddo
      


c compute the integral
	
      res=rmaxtot*rnacctot/rntot*scalemax

c compute the error on the integral
	
      errres=rmaxtot*scalemax*sqrt(((rnacctot/rntot)-(rnacctot/rntot)
     &     **2)/rntot)


c print results


      print*, 'number of events over max.', novermax
        
      print*, 'total integral=',res*(i_emutau+1)

      print*, 'total error=',errres*(i_emutau+1)
        
      print*, 'total unweighted events',int(rnacctot) 
        
*23/11/04
c      print*, 'total weighted events',int(rntot)
      if (rntot.le.2.d9) then
        print*, 'total weighted events',int(rntot)
      else
        print*, 'total weighted events',rntot
      endif
*23/11/04end
        
      print*, '  '
      print*, '  '
      print*, ' number of events over max,' 
      print*, ' unweighted evts, weighted evts and'
      print*, ' integrals for various proc and chan  '
        

      do i=1,nfiles
        print *, '     '
        print*, '         '
        print *,'proc=',i  
        print *,'noverm', noverm(i)
*23/11/04
c        print *,'unweighted and weighted events',
c     &       int(rnacc(i)), int(rntotproc(i)) 
      if (rntotproc(i).le.2.d9) then
        print *,'unweighted and weighted events',
     &       int(rnacc(i)), int(rntotproc(i))
      else
        print *,'unweighted and weighted events',
     &       int(rnacc(i)), rntotproc(i)
      endif
*23/11/04end
        exint=rmaxtot/rntot*rnacc(i)*scalemax
        if (int(rnacc(i)).gt.0) then
          exerr=exint/sqrt(rnacc(i))
        else
          exerr=0.d0
        endif
        
        exerr = rmaxtot*scalemax*sqrt(((rnacc(i)/rntot)-(rnacc(i)
     &       /rntot)**2)/rntot)
        
        print*, 'extimated integral', exint*(i_emutau+1)
        print*, 'extimated error', exerr*(i_emutau+1)
      enddo 

      IF(iwrite_event.EQ.1)THEN
        WRITE(23,'(A)')'</LesHouchesEvents>'      
        CLOSE(23)
      ENDIF

      stop
      end
