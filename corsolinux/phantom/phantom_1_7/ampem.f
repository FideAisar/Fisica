c amp -- Optimized version: reorders indeces in DO loops and avoids 
c redundant double-counting of interference terms in |Amp|**2
*
* Last update: Jun 26, 2007
*
* This routine computes the amplitude for all processes with
* 2 incoming fermions and  6 outgoing, making use of the  
* three routines fourw.f , twoztwow.f and fourz.f
* It makes use of the four momenta p(0:3,8) (all considered
* outgoing) and of the particle identities idp(8) (all outgoing)

*  idp_inout = -1/+1 for incoming/outgoing particles NOT PASSED 
*    as  symfact is given

*  igr(i) = 1/0 if amplitude to be computed 
*   ( i=1 4W,  2 2Z2W, 3 another 2Z2W, 4 4Z )
*  igo (8,4) order of the eight particles for the call of the
*     four igr amplitudes above. for instance the second 2Z2W,
*     corresponding to igr(3),if active, must be called as
*          CALL twoztwow(p(0,igo(1,3)),p(0,igo(2,3)),
*     &               p(0,igo(3,3)),p(0,igo(4,3)),
*     &               p(0,igo(5,3)),p(0,igo(6,3)),
*     &               p(0,igo(7,3)),p(0,igo(8,3)),
*     &               idp(igo(1,3)),idp(igo(2,3)),
*     &               idp(igo(3,3)),idp(igo(4,3)),
*     &               idp(igo(5,3)),idp(igo(6,3)),
*     &               idp(igo(7,3)),idp(igo(8,3)),
*     &               cresw(..)...)
*
*  ionesh 1/0 if oneshot mode is on/off
*  res  result
*  iord order of the chosen LH configuration

* It must be used together with the function isignperm(igo)
*  and  perm(idp, igr, igo, isign,
*     &     n4w,n2z2w,n4z, ig4w, ig2z2w,ig4z,isign4w,isign2z2w,isign4z)

      subroutine amp(p,idp,igr,igo,ionesh,symfact,res,iord)

      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)

      dimension p(0:3,8)
      dimension idp(8),igr(4),igo(8,4),iord(8)
      
      dimension isign(4)
      dimension ii(8)

      parameter (maxperm4w=4,maxperm2z2w=18*2,max4z=24)
      dimension ig4w(8,maxperm4w), ig2z2w(8,maxperm2z2w),ig4z(8,max4z)
      dimension isign4w(maxperm4w),isign2z2w(maxperm2z2w),isign4z(max4z)
      
      parameter (namp2max=64)

      dimension isig(namp2max),ipos(8,namp2max)
       
      dimension amp2(namp2max),factcol(namp2max,namp2max)

c giuseppe ILC 26/06/2007: already defined in coupling.f
c      parameter (czero=(0.d0,0.d0))
c end giuseppe ILC 26/06/2007
     
* declaration for record of double structure in fow 
*    at the end we have cres(namp2max).z(2,2,2,2).z(2,2,2,2)
      type zz
         double complex z(2,2,2,2)
      end type

      type z4
         double complex za(2,2)
      end type

c diogo -- include za2(2,2) and lbquark
      type fow
         double complex w,zaux(2,2,2,2),za2(2,2)
         type(zz) z(2,2,2,2)
         type(z4) za(2,2,2,2,2,2)
      end type
      type(fow) cres(namp2max)

      logical lbquark

c diogoend
      
* the choiche cres(i).za(2,2,2,2,2,2).za(2,2)
*  is to conform to fourz .


      COMMON /pharand/ idum

ctot sandro 4/2/07 
c      DATA ifirst/1/
      DATA ifirst/1/,igluons/0/
ctotend sandroend

c giuseppe ILC 26/06/2007
      INCLUDE 'common.h'
c end giuseppe ILC 26/06/2007

c diogo 27/04/2010
      include 'common_unitarization.h'
c end diogo

      if (ionesh.eq.1.or.ifirst.eq.1.and.ionesh.eq.0) then
 
* determine the sign of the available igo's

        do i=1,4
          if (igr(i).eq.1) then
ctot sandro 4/2/07
c            isign(i)=isignperm(igo(1,i))
            isign(i)=isignperm(igo(1,i),igluons)
ctotend sandroend
          endif
        enddo

* determine the number of calls for fourw, twoztwow, fourz
*  and their orders and signs

        call perm(idp, igr, igo, isign,
     &       n4w,n2z2w,n4z, ig4w, ig2z2w,ig4z,isign4w,isign2z2w,isign4z)

        namp2=n4z+n2z2w+n4w
       
*        print*, 'namp2=',namp2
*        print*, 'n4z=',n4z
*        print*, 'n2z2w=',n2z2w
*        print*, 'n4w=',n4w
			
        do i=1,n4z
          isig(i)=isign4z(i)
          do j=1,8
            ipos(j,i)=ig4z(j,i)
          enddo
        enddo

        do i=1,n2z2w
          k=n4z+i
          isig(k)=isign2z2w(i)
          do j=1,8
            ipos(j,k)=ig2z2w(j,i)
          enddo
        enddo

        n4ze2z2w=n4z+n2z2w
        do i=1,n4w
          k=n4ze2z2w+i
          isig(k)=isign4w(i)
          do j=1,8
            ipos(j,k)=ig4w(j,i)
          enddo
        enddo

* determine color factors
*  

        do j=1,namp2
          do i=j,namp2
            factcol(i,j)=rcolor(ipos(1,i),ipos(1,j),idp)
            factcol(j,i)=factcol(i,j)
          enddo
        enddo

        ifirst=0
      endif

* set to zero to avoid previous unwanted components


      do i=1,namp2
        do i4=1,2
        do i3=1,2
        do i2=1,2
        do i1=1,2
        do i8=1,2
        do i7=1,2
        do i6=1,2
        do i5=1,2
	    cres(i)%z(i1,i2,i3,i4)%z(i5,i6,i7,i8)=czero
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
      enddo

c diogo
*  see if there are b-quarks around (since there are no top quarks in the 
*  final state we need a pair 5 -5, so we can test just particles
      lbquark =.false.
      do i=1,7,2
	if (idp(i).eq.5) then
	 lbquark=.true.
       endif   
      enddo
      
* b-quark present OR flag mass=1     
* sig
*      if (lbquark.or.i_massive.eq.1) then
      if (lbquark.and.i_signal.eq.0.or.i_massive.eq.1) then
* sigend                 
      res=0.d0

      do i=1,n4z
	
	amp2(i)=0.d0      
        call fourz(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres(i)%za(1,1,1,1,1,1)%za(1,1))

        do i6=1,2
          do i5=1,2
            do i4=1,2
              do i3=1,2
                do i2=1,2
                do i1=1,2
                do i8=1,2
                do i7=1,2
                 ii(ipos(1,i))=i1
                 ii(ipos(2,i))=i2
                 ii(ipos(3,i))=i3
                 ii(ipos(4,i))=i4
                 ii(ipos(5,i))=i5
                 ii(ipos(6,i))=i6
                 ii(ipos(7,i))=i7
                 ii(ipos(8,i))=i8

                 cres(i)%z(ii(1),ii(2),ii(3),
     &               ii(4))%z(ii(5),ii(6),
     &               ii(7),ii(8))=
     &               cres(i)%za(i1,i2,i3,i4,i5,i6)%za(i7,i8)*isig(i)
     
     
* diagonal terms      
            amp2(i)=amp2(i)+
     &              cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8))*
     &              conjg(cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8)))*
     &              factcol(i,i) 

* interference terms

        do j=1,i-1
	    
              res = res +
     &             2.d0*cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8))*
     &             conjg(cres(j)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8)))*
     &             factcol(i,j) 

        enddo
	     
                enddo
                enddo
                enddo
                enddo
              enddo
            enddo
          enddo
        enddo
	
	res=res+amp2(i)

      enddo


      do i=n4z+1,n4z+n2z2w
      
        amp2(i)=0.d0
        call twoztwow(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres(i)%zaux)

        do i4=1,2
          do i3=1,2
            do i2=1,2
              do i1=1,2
                ii(ipos(1,i))=i1
                ii(ipos(2,i))=i2
                ii(ipos(3,i))=i3
                ii(ipos(4,i))=i4
                ii(ipos(5,i))=2
                ii(ipos(6,i))=2
                ii(ipos(7,i))=2
                ii(ipos(8,i))=2
                cres(i)%z(ii(1),ii(2),ii(3),
     &               ii(4))%z(ii(5),ii(6),
     &               ii(7),ii(8))=
     &               cres(i)%zaux(i1,i2,i3,i4)*isig(i)
     
* diagonal terms      
            amp2(i)=amp2(i)+
     &              cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8))*
     &              conjg(cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8)))*
     &              factcol(i,i) 

* interference terms

        do j=1,i-1

              res = res +
     &             2.d0*cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8))*
     &             conjg(cres(j)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8)))*
     &             factcol(i,j) 

        enddo
		     
              enddo
            enddo
          enddo
        enddo
	
	res=res+amp2(i)

      enddo

* should be zero, because that isn't 4W with b quark
      do i=n4ze2z2w+1,namp2
        
        amp2(i)=0.d0
        call fourw(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres(i)%w)

        cres(i)%z(2,2,2,2)%z(2,2,2,2)=cres(i)%w*isig(i)
      
* diagonal terms      
            amp2(i)=amp2(i)+
     &              cres(i)%z(2,2,2,2)%z(2,2,2,2)*
     &              conjg(cres(i)%z(2,2,2,2)%z(2,2,2,2))*
     &              factcol(i,i) 

* interference terms

        do j=1,i-1
	    
              res = res +
     &             2.d0*cres(i)%z(2,2,2,2)%z(2,2,2,2)*
     &             conjg(cres(j)%z(2,2,2,2)%z(2,2,2,2))*
     &             factcol(i,j) 

        enddo
        
	res=res+amp2(i)
	
      enddo

* b-quark absent AND imass=0      
      else
      
*       print*, 'ids in order of process i=3:', 
*     &   idp(ipos(1,3)),idp(ipos(2,3)),idp(ipos(3,3)),idp(ipos(4,3)),
*     &   idp(ipos(5,3)),idp(ipos(6,3)),idp(ipos(7,3)),idp(ipos(8,3))
      
      
       res=0.d0
            
      do i=1,n4z
      
        amp2(i)=0.d0
        call fourz_massless(
     &         p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres(i)%zaux)

          do i7=1,2
          do i5=1,2
          do i3=1,2
          do i1=1,2

            ii(ipos(1,i))=i1
            ii(ipos(2,i))=i1
            ii(ipos(3,i))=i3
            ii(ipos(4,i))=i3
            ii(ipos(5,i))=i5
            ii(ipos(6,i))=i5
            ii(ipos(7,i))=i7
            ii(ipos(8,i))=i7

            cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8))=
     &               cres(i)%zaux(i1,i3,i5,i7)*isig(i)
     
* diagonal terms      
            amp2(i)=amp2(i)+
     &              cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8))*
     &              conjg(cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8)))*
     &              factcol(i,i) 

* interference terms

        do j=1,i-1
	    
              res = res +
     &             2.d0*cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8))*
     &             conjg(cres(j)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8)))*
     &             factcol(i,j) 

        enddo
	    
	  enddo
	  enddo
	  enddo
	  enddo
        
*	print*, 'amp2(',i,')=',amp2(i)
	res=res+amp2(i)
	 
      enddo

      do i=n4z+1,n4z+n2z2w
             
        amp2(i)=0.d0
        call twoztwow_massless(
     &         p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres(i)%za2)

          do i3=1,2
          do i1=1,2
                ii(ipos(1,i))=i1
                ii(ipos(2,i))=i1
                ii(ipos(3,i))=i3
                ii(ipos(4,i))=i3
                ii(ipos(5,i))=2
                ii(ipos(6,i))=2
                ii(ipos(7,i))=2
                ii(ipos(8,i))=2
                cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                      ii(5),ii(6),ii(7),ii(8))=
     &               cres(i)%za2(i1,i3)*isig(i)

* diagonal terms 
            amp2(i)=amp2(i)+
     &              cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8))*
     &              conjg(cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8)))*
     &              factcol(i,i)     

* interference terms

        do j=1,i-1
	    
              res = res +
     &             2.d0*cres(i)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8))*
     &             conjg(cres(j)%z(ii(1),ii(2),ii(3),ii(4))%z(
     &                  ii(5),ii(6),ii(7),ii(8)))*
     &             factcol(i,j) 

        enddo

          enddo
          enddo
*        print*, 'amp2(',i,')=',amp2(i)
	res=res+amp2(i)
	
      enddo

      do i=n4ze2z2w+1,namp2
	     
        amp2(i)=0.d0	
        call fourw(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres(i)%w)

        cres(i)%z(2,2,2,2)%z(2,2,2,2)=cres(i)%w*isig(i)

* diagonal terms  
            amp2(i)=amp2(i)+
     &              cres(i)%z(2,2,2,2)%z(2,2,2,2)*
     &              conjg(cres(i)%z(2,2,2,2)%z(2,2,2,2))*
     &              factcol(i,i)     
 
* interference terms

        do j=1,i-1
	    
              res = res +
     &             2.d0*cres(i)%z(2,2,2,2)%z(2,2,2,2)*
     &             conjg(cres(j)%z(2,2,2,2)%z(2,2,2,2))*
     &             factcol(i,j) 
     
        enddo
	
*	print*, 'amp2(',i,')=',amp2(i)
	res=res+amp2(i)
        
      enddo

      endif
c diogoend


c giuseppe ILC 26/06/2007
* spin average 1/4
* color average 1/9 (only in case of initial state quarks)
      if(i_coll.eq.1 .or. i_coll.eq.2)then ! hadron colliders
        rinvclav=9.d0
      elseif(i_coll.eq.3)then   ! e+e- collider
        rinvclav=1.d0
      endif

      res=res/4.d0/rinvclav/symfact

c end giuseppe ILC 26/06/2007

*extract the configuration nchosen for color flow

      amp2tot=0.d0
      do i=1,namp2
	amp2tot=amp2tot+amp2(i)
      enddo

      do i=1,namp2
	amp2(i)=amp2(i)/amp2tot
      enddo

      i=1
      am2=amp2(1)
      zzran=ran2(idum)
      do while (zzran.gt.am2.and.i.lt.namp2)
	i=i+1
	am2=am2+amp2(i)
      enddo	

      nchosen=i

*compute iord as the inverse of ipos(i,nchosen) 

      do i=1,8
	iord(ipos(i,nchosen))=i
      enddo


      
      return
      end

