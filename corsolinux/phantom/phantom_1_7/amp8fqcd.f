c amp8fqcd -- Optimized version: avoids loops on identically zero color
c configurations; avoids redundant double counting in calculating 
c interference terms. 
c In structure cres containing the final amplitude, the indeces 
c representing the colour configuration are the first to the right.
*
* Last update: Jun 25, 2007
*
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
*     &               cres(..)...)
*
*  ionesh 1/0 if oneshot mode is on/off
*  res  result
*  iord order of the chosen LH configuration

* It must be used together with the function isignperm(igo,igluons)
*  and the routine  perm(idp, igr, igo, isign,
*     &     n4w,n2z2w,n4z, ig4w, ig2z2w,ig4z,isign4w,isign2z2w,isign4z)

      subroutine ampqcd(p,idp,igr,igo,ionesh,symfact,res,iord)

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
       
      dimension amp2(namp2max),amp2k(7,namp2max)
      dimension factcol(7,namp2max,7,namp2max)
      dimension iactive(7,namp2max)
      dimension rcoup(7)

* declaration for record of double structure in fow 
*    at the end we have cresh(2,2,2,2).h(2,2,2,2).col(7,namp2max)

      type tot
         double complex col(7,namp2max)
      end type

      type ftot
         type(tot) h(2,2,2,2)
      end type
      type(ftot) cresh(2,2,2,2)


*declaration for fourw

      type fww
         double complex w
      end type
      type(fww) cresww(7)

*declaration for twoztwow

      type fzw
         double complex zw(2,2,2,2)
      end type
      type(fzw) creszw(6)

*diogo massless 30/12/08
      dimension cres4zaux(2,2,2,2,7),cres2z2waux(2,2,6)
      
      logical lbquark 
*diogoend      

* the choiche creszz.zz(2,2,2,2,2,2).zz(2,2,7)
*  is to conform to fourz .
      type zz
         double complex zz(2,2,7)
      end type

      type fzz
         type(zz) zz(2,2,2,2,2,2)
      end type
      type(fzz) creszz

* color structure

      parameter (maxcouple=4,maxtriple=2)

      type tre
         integer ithree(3)
      end type

      type due
         integer itwo(2)
      end type

      type uno 
         integer ncouple,ntriple
         type(due)  couple(maxcouple)
         type(tre)  triple(maxtriple)
      end type
      type(uno)  col(7,namp2max)



      COMMON /pharand/ idum

      DATA ifirst/1/,igluons/0/

      INCLUDE 'common.h'

      e6=elcharge2**3
      g2e4=elcharge2**2*gs2
*sandro2/3/07
c      if (i_nofullew.eq.1) then
      if (i_pertorder.eq.2) then
*sandro2/3/07end
        rcoup(1)=0.d0
        kmin=2
      else
        rcoup(1)=e6
        kmin=1
      endif
ctot
c      do i=2,7
c        rcoup(i)=g2e4
c      enddo

*sandro2/3/07
c      if (i_nofullew.eq.-2) then
      if (i_pertorder.eq.0) then
*sandro2/3/07end
        do i=2,7
          rcoup(i)=0.d0
        enddo
        kmax=1
      else
        do i=2,7
          rcoup(i)=g2e4
        enddo
        kmax=7
      endif
ctotend

      if (ionesh.eq.1.or.ifirst.eq.1.and.ionesh.eq.0) then
 
* determine the sign of the available igo's

        do i=1,4
          if (igr(i).eq.1) then
            isign(i)=isignperm(igo(1,i),igluons)
          endif
        enddo


* determine the number of calls for fourw, twoztwow, fourz
*  and their orders and signs

        call perm(idp, igr, igo, isign,
     &       n4w,n2z2w,n4z, ig4w, ig2z2w,ig4z,isign4w,isign2z2w,isign4z)

        namp2=n4z+n2z2w+n4w


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


* determine the number of hadronic couples
        nhadcoup=0
        do i=1,8
          if (abs(idp(i)).le.6)  nhadcoup=nhadcoup+1
        enddo
        nhadcoup=nhadcoup/2

* fill for every configuration 1...namp2 the seven color structures
        do i=1,namp2
ctot
c          do k=kmin,7
          do k=kmin,kmax
ctotend

* l'ordine della chiamata e' ipos(1,i)....(ipos(8,i)
            if (k.eq.1) then
              iactive(k,i)=1
              col(k,i)%ncouple=nhadcoup
              col(k,i)%ntriple=0
              icount=1
              do j=1,8,2
                if (idp(ipos(j,i)).le.6) then
                  col(k,i)%couple(icount)%itwo(1)=ipos(j,i)
                  col(k,i)%couple(icount)%itwo(2)=ipos(j+1,i)
                  icount=icount+1
                endif
              enddo  
            else
              col(k,i)%ncouple=nhadcoup-2
              col(k,i)%ntriple=2
              if (k.eq.2) then  
                if (idp(ipos(1,i)).le.6.and.idp(ipos(3,i)).le.6) then
                  iactive(k,i)=1  
                  col(k,i)%triple(1)%ithree(1)=ipos(1,i)
                  col(k,i)%triple(1)%ithree(2)=20
                  col(k,i)%triple(1)%ithree(3)=ipos(2,i)
                  col(k,i)%triple(2)%ithree(1)=ipos(3,i)
                  col(k,i)%triple(2)%ithree(2)=20
                  col(k,i)%triple(2)%ithree(3)=ipos(4,i)
                  icount=1
                  if (idp(ipos(5,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(5,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(6,i)
                    icount=icount+1
                  endif
                  if (idp(ipos(7,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(7,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(8,i)
                  endif
                else
                  iactive(k,i)=0
                endif
              elseif (k.eq.3) then  
                if (idp(ipos(1,i)).le.6.and.idp(ipos(5,i)).le.6) then
                  iactive(k,i)=1  
                  col(k,i)%triple(1)%ithree(1)=ipos(1,i)
                  col(k,i)%triple(1)%ithree(2)=20
                  col(k,i)%triple(1)%ithree(3)=ipos(2,i)
                  col(k,i)%triple(2)%ithree(1)=ipos(5,i)
                  col(k,i)%triple(2)%ithree(2)=20
                  col(k,i)%triple(2)%ithree(3)=ipos(6,i)
                  icount=1
                  if (idp(ipos(3,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(3,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(4,i)
                    icount=icount+1
                  endif
                  if (idp(ipos(7,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(7,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(8,i)
                  endif
                else
                  iactive(k,i)=0
                endif
              elseif (k.eq.4) then  
                if (idp(ipos(3,i)).le.6.and.idp(ipos(5,i)).le.6) then
                  iactive(k,i)=1  
                  col(k,i)%triple(1)%ithree(1)=ipos(3,i)
                  col(k,i)%triple(1)%ithree(2)=20
                  col(k,i)%triple(1)%ithree(3)=ipos(4,i)
                  col(k,i)%triple(2)%ithree(1)=ipos(5,i)
                  col(k,i)%triple(2)%ithree(2)=20
                  col(k,i)%triple(2)%ithree(3)=ipos(6,i)
                  icount=1
                  if (idp(ipos(1,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(1,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(2,i)
                    icount=icount+1
                  endif
                  if (idp(ipos(7,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(7,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(8,i)
                  endif
                else
                  iactive(k,i)=0
                endif
              elseif (k.eq.5) then  
                if (idp(ipos(1,i)).le.6.and.idp(ipos(7,i)).le.6) then
                  iactive(k,i)=1  
                  col(k,i)%triple(1)%ithree(1)=ipos(1,i)
                  col(k,i)%triple(1)%ithree(2)=20
                  col(k,i)%triple(1)%ithree(3)=ipos(2,i)
                  col(k,i)%triple(2)%ithree(1)=ipos(7,i)
                  col(k,i)%triple(2)%ithree(2)=20
                  col(k,i)%triple(2)%ithree(3)=ipos(8,i)
                  icount=1
                  if (idp(ipos(3,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(3,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(4,i)
                    icount=icount+1
                  endif
                  if (idp(ipos(5,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(5,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(6,i)
                  endif
                else
                  iactive(k,i)=0
                endif
              elseif (k.eq.6) then  
                if (idp(ipos(3,i)).le.6.and.idp(ipos(7,i)).le.6) then
                  iactive(k,i)=1  
                  col(k,i)%triple(1)%ithree(1)=ipos(3,i)
                  col(k,i)%triple(1)%ithree(2)=20
                  col(k,i)%triple(1)%ithree(3)=ipos(4,i)
                  col(k,i)%triple(2)%ithree(1)=ipos(7,i)
                  col(k,i)%triple(2)%ithree(2)=20
                  col(k,i)%triple(2)%ithree(3)=ipos(8,i)
                  icount=1
                  if (idp(ipos(1,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(1,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(2,i)
                    icount=icount+1
                  endif
                  if (idp(ipos(5,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(5,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(6,i)
                  endif
                else
                  iactive(k,i)=0
                endif
              elseif (k.eq.7) then  
                if (idp(ipos(5,i)).le.6.and.idp(ipos(7,i)).le.6) then
                  iactive(k,i)=1  
                  col(k,i)%triple(1)%ithree(1)=ipos(5,i)
                  col(k,i)%triple(1)%ithree(2)=20
                  col(k,i)%triple(1)%ithree(3)=ipos(6,i)
                  col(k,i)%triple(2)%ithree(1)=ipos(7,i)
                  col(k,i)%triple(2)%ithree(2)=20
                  col(k,i)%triple(2)%ithree(3)=ipos(8,i)
                  icount=1
                  if (idp(ipos(1,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(1,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(2,i)
                    icount=icount+1
                  endif
                  if (idp(ipos(3,i)).le.6) then
                    col(k,i)%couple(icount)%itwo(1)=ipos(3,i)
                    col(k,i)%couple(icount)%itwo(2)=ipos(4,i)
                  endif
                else
                  iactive(k,i)=0
                endif
              endif
            endif

          enddo
        enddo

* determine color factors
*  
ctot
c        kmax=7
ctotend

csandro31/1/07
c        lim=namp2*kmax
c        do m1=0,lim-1
c          do m2=m1,lim-1
c            ki=mod(m1,kmax)+1
c            i=(m1+1-ki)/kmax +1
c            kj=mod(m2,kmax)+1
c            j=(m2+1-kj)/kmax +1


c diogo : avoid computing color factor for interference between 
c same fermion lines with different color configuration, which 
c are always zero

       do i=1,namp2
        
	do ki=kmin,kmax
         if (iactive(ki,i).eq.1) then       
* diagonal terms
          call coleval(col(ki,i),col(ki,i),rescolfact)
	  factcol(ki,i,ki,i)=rescolfact
* interference terms
         do j=1,i-1
	  do kj=kmin,kmax
	   if (iactive(kj,j).eq.1) then
	    call coleval(col(kj,j),col(ki,i),rescolfact)
	    factcol(ki,i,kj,j)=rescolfact
	   else
	    factcol(ki,i,kj,j)=0.d0
	   endif
	   
	   factcol(kj,j,ki,i)=factcol(ki,i,kj,j)
	   
	  enddo
	 enddo
	  
        endif	
       enddo
      enddo 
c diogoend

*      krange=kmax-kmin+1
*      lim=namp2*krange
*      do m1=0,lim-1
*        ki=mod(m1,krange)+kmin
*        i=m1/krange+1
*          do m2=m1,lim-1
*            kj=mod(m2,krange)+kmin
*            j=m2/krange +1
csandro31/1/07end


* fill color structure
*            if (iactive(kj,j).eq.1.and.iactive(ki,i).eq.1) then
*              call coleval(col(kj,j),col(ki,i),rescolfact)
*              factcol(ki,i,kj,j)=rescolfact
*            else
*              factcol(ki,i,kj,j)=0.d0
*            endif
*            factcol(kj,j,ki,i)=factcol(ki,i,kj,j)
*          enddo
*        enddo


        ifirst=0
      endif   !(ionesh.eq.1.or.ifirst.eq.1.and.ionesh.eq.0)


* set to zero to avoid previous unwanted components

      do i4=1,2
        do i3=1,2
          do i2=1,2
            do i1=1,2
              do i8=1,2
                do i7=1,2
                  do i6=1,2
                    do i5=1,2
                      do i=1,namp2
c here the maximum limit must remain 7 and not kmax because 
c   a process 2z2w fills up only to k=6 and the seventh k item
c can remain different from zero in oneshot if not set to zero here
                        do k=kmin,7
                          cresh(i1,i2,i3,i4)%h(i5,i6,i7,i8)%col(k,i)
     &                         = czero
                        enddo
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
*diogoend

* b-quark present OR flag mass=1     
      if (lbquark.or.i_massive.eq.1) then
      
      res=0.d0
      
      do i=1,n4z
        amp2(i)=0.d0
	
          call fourzqcd(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),creszz%zz(1,1,1,1,1,1)%zz(1,1,1))

ctot
c         do k=kmin,7
         do k=kmin,kmax
ctotend
         amp2k(k,i)=0.d0
	 
           if(iactive(k,i).eq.1)then
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

                              cresh(ii(1),ii(2),ii(3),ii(4))
     &                             %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &                             = creszz%zz(i1,i2,i3,i4,i5,i6)
     &                             %zz(i7,i8,k)*isig(i)*rcoup(k)

* diagonal terms      
            amp2k(k,i)=amp2k(k,i)+
     &              cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &                 *conjg(cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i))
     &                 *factcol(k,i,k,i)
     
     
* interference terms 

c diogo avoid computing interference between different color 
c configuration of same fermion line structure
        do j=1,i-1
	do kj=kmin,kmax
	  if (iactive(kj,j).eq.1) then
	     
              res = res +
     &             2.d0*cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &                *conjg(cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(kj,j))
     &                *factcol(k,i,kj,j)
          endif   

        enddo
	enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          endif
	  amp2(i)=amp2(i)+amp2k(k,i)
	  
        enddo
        res=res+amp2(i)
	
      enddo

      do i=n4z+1,n4z+n2z2w
      
      amp2(i)=0.d0
      
        call twoztwowqcd(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),creszw(1)%zw(1,1,1,1))

ctot
c        do k=kmin,6 
c          if(iactive(k,i).eq.1)then
        do k=kmin,kmax 
	
	amp2k(k,i)=0.d0 
	
          if(iactive(k,i).eq.1.and.k.ne.7)then
ctotend
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

                    cresh(ii(1),ii(2),ii(3),ii(4))
     &                   %h(ii(5),ii(6),ii(7),ii(8))%col(k,i) 
     &                   = creszw(k)%zw(i1,i2,i3,i4)*isig(i)*rcoup(k)

* diagonal terms      
            amp2k(k,i)=amp2k(k,i)+
     &              cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &                 *conjg(cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i))
     &                 *factcol(k,i,k,i)

* interference terms !!!!!!!!!!!!!!!!!!!!!!!!!!!????

        do j=1,i-1
	do kj=kmin,kmax
	  if (iactive(kj,j).eq.1) then
	     
              res = res +
     &             2.d0*cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &                *conjg(cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(kj,j))
     &                *factcol(k,i,kj,j)
          endif   

        enddo
	enddo
	
                  enddo
                enddo
              enddo
            enddo
          endif
	  amp2(i)=amp2(i)+amp2k(k,i)
	  
        enddo
        res=res+amp2(i)
	
      enddo

      do i=n4ze2z2w+1,namp2
      amp2(i)=0.d0
      
        call fourwqcd(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cresww(1)%w)

ctot
c        do k=kmin,7
        do k=kmin,kmax
ctotend
         amp2k(k,i)=0.d0
          if(iactive(k,i).eq.1)then
            cresh(2,2,2,2)%h(2,2,2,2)%col(k,i) 
     &           = cresww(k)%w*isig(i)*rcoup(k)
     
*diagonal terms      
            amp2k(k,i)=amp2k(k,i)+
     &       cresh(2,2,2,2)%h(2,2,2,2)%col(k,i)
     &       *conjg(cresh(2,2,2,2)%h(2,2,2,2)%col(k,i))
     &                 *factcol(k,i,k,i)

* interference terms !!!!!!!!!!!!!!!!!!!!!!!!!!!????

        do j=1,i-1
	do kj=kmin,kmax
	  if (iactive(kj,j).eq.1) then
	     
            res = res +
     &        2.d0*cresh(2,2,2,2)%h(2,2,2,2)%col(k,i)
     &        *conjg(cresh(2,2,2,2)%h(2,2,2,2)%col(kj,j))
     &                *factcol(k,i,kj,j)
          endif   

        enddo
	enddo
	
          endif
	  amp2(i)=amp2(i)+amp2k(k,i)
	  
        enddo
        res=res+amp2(i)
	
      enddo

*b-quark ausent AND imass=0
      else
      
      res=0.d0

      do i=1,n4z
       amp2(i)=0.d0

          call fourzqcd_massless(
     &         p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres4zaux)

ctot
c         do k=kmin,7
         do k=kmin,kmax
ctotend
           amp2k(k,i)=0.d0
           if(iactive(k,i).eq.1)then
            
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

                cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &   	    = cres4zaux(i1,i3,i5,i7,k)*isig(i)*rcoup(k)
     
     
* diagonal terms      
            amp2k(k,i)=amp2k(k,i)+
     &              cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &                 *conjg(cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i))
     &                 *factcol(k,i,k,i)
     
* interference terms !!!!!!!!!!!!!!!!!!!!!!!!!!!????

        do j=1,i-1
	do kj=kmin,kmax
	  if (iactive(kj,j).eq.1) then
	     
              res = res +
     &             2.d0*cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &                *conjg(cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(kj,j))
     &                *factcol(k,i,kj,j)
          endif   

        enddo
	enddo     

	    enddo
	    enddo
	    enddo
	    enddo	    
	    
          endif
	  amp2(i)=amp2(i)+amp2k(k,i)  
	  
        enddo
	
	res=res+amp2(i)

      enddo




      do i=n4z+1,n4z+n2z2w
      
       amp2(i)=0.d0
	      
        call twoztwowqcd_massless(
     &         p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres2z2waux)

ctot
c        do k=kmin,6 
c          if(iactive(k,i).eq.1)then
        do k=kmin,kmax
	  
	  amp2k(k,i)=0.d0 
          if(iactive(k,i).eq.1.and.k.ne.7)then
ctotend
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

                    cresh(ii(1),ii(2),ii(3),ii(4))
     &                   %h(ii(5),ii(6),ii(7),ii(8))%col(k,i) 
     &                   = cres2z2waux(i1,i3,k)*isig(i)*rcoup(k)

* diagonal terms      
            amp2k(k,i)=amp2k(k,i)+
     &              cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &                 *conjg(cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i))
     &                 *factcol(k,i,k,i)

* interference terms !!!!!!!!!!!!!!!!!!!!!!!!!!!????

        do j=1,i-1
	do kj=kmin,kmax
	  if (iactive(kj,j).eq.1) then
	     
              res = res +
     &             2.d0*cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(k,i)
     &                *conjg(cresh(ii(1),ii(2),ii(3),ii(4))
     &	            %h(ii(5),ii(6),ii(7),ii(8))%col(kj,j))
     &                *factcol(k,i,kj,j)
          endif   

        enddo
	enddo
	
            enddo
            enddo
	    
          endif	  
	  amp2(i)=amp2(i)+amp2k(k,i)
	  
        enddo
	res=res+amp2(i)

      enddo

      do i=n4ze2z2w+1,namp2
      
        amp2(i)=0.d0
      
        call fourwqcd(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cresww(1)%w)

ctot
c        do k=kmin,7
        do k=kmin,kmax
ctotend
          amp2k(k,i)=0.d0
          if(iactive(k,i).eq.1)then
            cresh(2,2,2,2)%h(2,2,2,2)%col(k,i) 
     &           = cresww(k)%w*isig(i)*rcoup(k)
     
* diagonal terms      
            amp2k(k,i)=amp2k(k,i)+
     &       cresh(2,2,2,2)%h(2,2,2,2)%col(k,i)
     &       *conjg(cresh(2,2,2,2)%h(2,2,2,2)%col(k,i))
     &                 *factcol(k,i,k,i)

* interference terms !!!!!!!!!!!!!!!!!!!!!!!!!!!????

        do j=1,i-1
	do kj=kmin,kmax
	  if (iactive(kj,j).eq.1) then
	     
            res = res +
     &        2.d0*cresh(2,2,2,2)%h(2,2,2,2)%col(k,i)
     &        *conjg(cresh(2,2,2,2)%h(2,2,2,2)%col(kj,j))
     &                *factcol(k,i,kj,j)
          endif   

        enddo
	enddo
     
          endif
	  amp2(i)=amp2(i)+amp2k(k,i)
	  
        enddo
        res=res+amp2(i)
	
      enddo
      
      endif  !if b-quark      

c giuseppe ILC 25/06/2007
* spin average 1/4
* color average 1/9 (only in case of initial state quarks)
      if(i_coll.eq.1 .or. i_coll.eq.2)then ! hadron colliders
        rinvclav=9.d0
      elseif(i_coll.eq.3)then   ! e+e- collider
        rinvclav=1.d0
      endif
      res=res/4.d0/rinvclav/symfact
c end giuseppe ILC 25/06/2007


** extract the configuration ichosen for color flow

      amp2tot=0.d0
      do i=1,namp2
	amp2tot=amp2tot+amp2(i)
      enddo
*ctest
      if (amp2tot.eq.0.d0) then
	print*, '   '
	print*, 'amp2tot = 0!!!!!!!!!' 
	print*, 'PLEASE CHECK THAT IT IS CORRECT FOR THIS PROCESS'
	stop
      endif
*ctestend

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

      ichosen=i


* extract the configuration kchosen for color flow
*sandro 22/1/07
c      amp2tot=0.d0
c      do k=kmin,7
c	 amp2tot=amp2tot+amp2k(k,i)
c      enddo
c
c      do k=kmin,7
c	 amp2(k)=amp2k(k,i)/amp2tot
c      enddo
c
c      k=kmin
c      am2=amp2(1)
c      zzran=ran2(idum)
c      do while (zzran.gt.am2.and.k.lt.7)
c	 k=k+1
c	 am2=am2+amp2(k)
c      enddo	  

      amp2tot=0.d0
ctot
c      do k=kmin,7
      do k=kmin,kmax
ctotend
	amp2tot=amp2tot+amp2k(k,ichosen)
      enddo

ctot
c      do k=kmin,7
      do k=kmin,kmax
ctotend
	amp2k(k,ichosen)=amp2k(k,ichosen)/amp2tot
      enddo

      k=kmin
      am2=amp2k(kmin,ichosen)
      zzran=ran2(idum)
ctot
c      do while (zzran.gt.am2.and.k.lt.7)
      do while (zzran.gt.am2.and.k.lt.kmax)
ctotend
	k=k+1
	am2=am2+amp2k(k,ichosen)
      enddo	 
*sandro 22/1/07 end
      kchosen=k


* compute iord as the inverse of ipos(i,nchosen) 
c
c old ampem
c      do i=1,8
c	 iord(ipos(i,nchosen))=i
c      enddo

      if (kchosen.eq.1) then
	do i=1,nhadcoup
	  iord(col(kchosen,ichosen)%couple(i)%itwo(1))=2*i-1	  
	  iord(col(kchosen,ichosen)%couple(i)%itwo(2))=2*i 
	enddo

      else
	iord(col(kchosen,ichosen)%triple(1)%ithree(1))=1
	iord(col(kchosen,ichosen)%triple(2)%ithree(3))=2
	iord(col(kchosen,ichosen)%triple(2)%ithree(1))=3
	iord(col(kchosen,ichosen)%triple(1)%ithree(3))=4
	do i=1,col(kchosen,ichosen)%ncouple
	  iord(col(kchosen,ichosen)%couple(i)%itwo(1))=4+2*i-1      
	  iord(col(kchosen,ichosen)%couple(i)%itwo(2))=4+2*i 
	enddo
      endif

***************

* fill the non hadronic positions if needed
      icount=2*nhadcoup+1
      do j=1,8
	if (abs(idp(ipos(j,ichosen))).gt.6) then
	  iord(ipos(j,ichosen))=icount
	  icount=icount+1
	endif
      enddo  

     
      return
      end

