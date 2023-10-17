c amp2g -- Optimized version: avoids loops on identically zero color 
c configurations; avoids redundant double counting in calculating 
c interference terms. 
c In structure cres containing the final amplitude, the indeces 
c representing the colour configuration are the first to the right.
*
* Last update: Jun 25, 2007
*
*

* This routine computes the amplitude for all processes with            
* 2 gluons and  6 fermions as  outgoing particles, making use of the    
* two  routines ggzww.f and gg3z.f.                                     
* It makes use of the four momenta p(0:3,8) (all considered             
* outgoing) and of the particle identities idp(8) (all outgoing)        
  
*  idp_inout = -1/+1 for incoming/outgoing particles NOT PASSED         
*    as  symfact is given                                               
  
*  igr(i) = 1/0 if amplitude to be computed                             
*   ( i=1 4W,  2 2Z2W, 3 another 2Z2W, 4 4Z ,) 5 2g1Z2W , 6 2g3Z       
  
*   but we use only igr(2) and igo(8,2) with a call igr(5) igo(1,5)     
*       from fxn.f                                                      
  
*  igo (8,2) order of the eight particles for the call of the           
*     four igr amplitudes above. For instance the process  1Z2W,        
*     corresponding to igr(1) locally,if active, must be called as      
*        call  ggzww(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),      
*     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),           
*     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),           
*     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),           
*     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),           
*     &       idp(ipos(8,i)),creszww)                                   
*    where ipos(1,i),ipos(2,i)... are the equivalent of igo(1,1),       
*    igo(2,1).... for the ith permutation                               
*                                                                       
*  ionesh 1/0 if oneshot mode is on/off                                 
*  res  result                                                          
*  iord order of the chosen LH configuration                            
  
* It must be used together with the function isignperm(igo)             
*  and the routine  perm(idp, igr, igo, isign,                          
*     &     n4w,n2z2w,n4z, ig4w, ig2z2w,ig4z,isign4w,isign2z2w,isign4z) 
  
********from ggzww                                                      
* cres  is on output the complex helicity amplitude                     
* cres=cres(2,2,2,2,12). The first two indeces refer to the             
* helicities of the gluons, the third and fourth one to the chiralities 
* of particles 3 and 4, while all other fermions have only chirality    
* index 2 (corresponding to -)                                          
* cres(i1,i2,i3,i4,1)  --> T[3,1,4]T[5,2,6][7,8]                        
* cres(i1,i2,i3,i4,2)  --> T[3,2,4]T[5,1,6][7,8]                        
* cres(i1,i2,i3,i4,3)  --> T[3,1,4]T[5,6][7,2,8]                        
* cres(i1,i2,i3,i4,4)  --> T[3,2,4]T[5,6][7,1,8]                        
* cres(i1,i2,i3,i4,5)  --> T[3,4]T[5,1,6][7,2,8]                        
* cres(i1,i2,i3,i4,6)  --> T[3,4]T[5,2,6][7,1,8]                        
* cres(i1,i2,i3,i4,7)  --> T[3,1,2,4]T[5,6][7,8]                        
* cres(i1,i2,i3,i4,8)  --> T[3,2,1,4]T[5,6][7,8]                        
* cres(i1,i2,i3,i4,9)  --> T[3,4]T[5,1,2,6][7,8]                        
* cres(i1,i2,i3,i4,10) --> T[3,4]T[5,2,1,6][7,8]                        
* cres(i1,i2,i3,i4,11) --> T[3,4]T[5,6][7,1,2,8]                        
* cres(i1,i2,i3,i4,12) --> T[3,4]T[5,6][7,2,1,8]                        
  
*      subroutine ggzww(p1,p2,p3,p4,p5,p6,p7,p8,                        
*     &              id1,id2,id3,id4,id5,id6,id7,id8,cres)              
  
********from gg3z                                                       
*      type gg                                                    
*         double complex gg(2,2,12)                                     
*      end type                                                    
*                                                                       
*      type zzz                                                   
*        type(gg) z(2,2,2,2,2,2)                                       
*      end type                                                    
*      type(zzz) cres                                                  
*                                                                       
*      subroutine ggzzz(p1,p2,p3,p4,p5,p6,p7,p8,                        
*     &              id1,id2,id3,id4,id5,id6,id7,id8,cres)              
  
  
  
* ingl= number of  incoming gluons for mean values.                     
  
      subroutine amp2g(p,idp,igr,igo,ionesh,
     &     symfact,ingl,res,iord)
  
      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)
  
      dimension p(0:3,8)
      dimension idp(8),igr(2),igo(8,2),iord(8)
  
      dimension isign(2)
  
      dimension ii(8)
  
      parameter (maxperm1z2w=4,max3z=6)
      dimension ig1z2w(8,maxperm1z2w),ig3z(8,max3z)
      dimension isign1z2w(maxperm1z2w),isign3z(max3z)
  
      parameter (namp2max=10)
  
      dimension isig(namp2max),ipos(8,namp2max)
  
      dimension amp2(namp2max),amp2k(12,namp2max)
      dimension factcol(12,namp2max,12,namp2max)
      dimension iactive(12,namp2max)
  
* declaration for record of double structure for final cres             
*    at the end we have cres(2,2).q(2,2,2,2,2,2).col(12,namp2max)
      type tot
         double complex col(12,namp2max)
      end type
  
      type ftot
         type(tot) q(2,2,2,2,2,2)
      end type
      type(ftot) cres(2,2)

*declaration for ggzww                                                  
  
      dimension creszww(2,2,2,2,12),creszww_massless(2,2,2,12)
  
*declaration for gg3z                                                   
      type res3z
        double complex caux(2,2,12)
      end type
      type(res3z) cres3z(2,2,2,2,2,2)
  
      dimension cres3z_massless(2,2,2,2,2,12)
      
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
      type(uno) col(12,namp2max)
  
      logical lbquark 
  
      COMMON /pharand/ idum
  
      DATA ifirst/1/,igluons/1/
  
      INCLUDE 'common.h'
  
      g4e8=elcharge2**4*gs2**2

  
      if (ionesh.eq.1.or.ifirst.eq.1.and.ionesh.eq.0) then  

* determine the sign of the available igo's                             
  
        do i=1,2
          if (igr(i).eq.1) then
            isign(i)=isignperm(igo(1,i),igluons)
          endif
        enddo
  
  
* determine the number of calls for 1z2w and 3z                         
*  and their orders and signs                                           
  
        call perm2g(idp, igr, igo, isign,
     &       n1z2w,n3z, ig1z2w, ig3z,isign1z2w,isign3z)

  
        namp2=n1z2w+n3z
  
        do i=1,n1z2w
          isig(i)=isign1z2w(i)
          do j=1,8
            ipos(j,i)=ig1z2w(j,i)
          enddo
        enddo
  
        do i=1,n3z
          k=n1z2w+i
          isig(k)=isign3z(i)
          do j=1,8
            ipos(j,k)=ig3z(j,i)
          enddo
        enddo

* determine the number of hadronic couples                              
        nhadcoup=0
        do i=1,8
          if (abs(idp(i)).le.6)  nhadcoup=nhadcoup+1
        enddo
        nhadcoup=nhadcoup/2
  
* fill for every configuration 1...namp2 the twelve color structures  
        do i=1,namp2
          do k=1,12
* l'ordine della chiamata e' ipos(1,i)....(ipos(8,i)                    
  
            if (k.eq.1) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-2
              if (idp(ipos(3,i)).le.6.and.idp(ipos(5,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(3,i)
                col(k,i)%triple(1)%ithree(2)=1
                col(k,i)%triple(1)%ithree(3)=ipos(4,i)
                col(k,i)%triple(2)%ithree(1)=ipos(5,i)
                col(k,i)%triple(2)%ithree(2)=2
                col(k,i)%triple(2)%ithree(3)=ipos(6,i)
                if (idp(ipos(7,i)).le.6) then
                  col(k,i)%couple(1)%itwo(1)=ipos(7,i)
                  col(k,i)%couple(1)%itwo(2)=ipos(8,i)
                endif
              else
                iactive(k,i)=0
              endif
            elseif (k.eq.2) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-2
              if (idp(ipos(3,i)).le.6.and.idp(ipos(5,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(3,i)
                col(k,i)%triple(1)%ithree(2)=2
                col(k,i)%triple(1)%ithree(3)=ipos(4,i)
                col(k,i)%triple(2)%ithree(1)=ipos(5,i)
                col(k,i)%triple(2)%ithree(2)=1
                col(k,i)%triple(2)%ithree(3)=ipos(6,i)
                if (idp(ipos(7,i)).le.6) then
                  col(k,i)%couple(1)%itwo(1)=ipos(7,i)
                  col(k,i)%couple(1)%itwo(2)=ipos(8,i)
                endif
              else
                iactive(k,i)=0
              endif
            elseif (k.eq.3) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-2
              if (idp(ipos(3,i)).le.6.and.idp(ipos(7,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(3,i)
                col(k,i)%triple(1)%ithree(2)=1
                col(k,i)%triple(1)%ithree(3)=ipos(4,i)
                col(k,i)%triple(2)%ithree(1)=ipos(7,i)
                col(k,i)%triple(2)%ithree(2)=2
                col(k,i)%triple(2)%ithree(3)=ipos(8,i)
                if (idp(ipos(7,i)).le.6) then
                  col(k,i)%couple(1)%itwo(1)=ipos(5,i)
                  col(k,i)%couple(1)%itwo(2)=ipos(6,i)
                endif
              else
                iactive(k,i)=0
              endif
            elseif (k.eq.4) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-2
              if (idp(ipos(3,i)).le.6.and.idp(ipos(7,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(3,i)
                col(k,i)%triple(1)%ithree(2)=2
                col(k,i)%triple(1)%ithree(3)=ipos(4,i)
                col(k,i)%triple(2)%ithree(1)=ipos(7,i)
                col(k,i)%triple(2)%ithree(2)=1
                col(k,i)%triple(2)%ithree(3)=ipos(8,i)
                if (idp(ipos(7,i)).le.6) then
                  col(k,i)%couple(1)%itwo(1)=ipos(5,i)
                  col(k,i)%couple(1)%itwo(2)=ipos(6,i)
                endif
              else
                iactive(k,i)=0
              endif
            elseif (k.eq.5) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-2
              if (idp(ipos(5,i)).le.6.and.idp(ipos(7,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(5,i)
                col(k,i)%triple(1)%ithree(2)=1
                col(k,i)%triple(1)%ithree(3)=ipos(6,i)
                col(k,i)%triple(2)%ithree(1)=ipos(7,i)
                col(k,i)%triple(2)%ithree(2)=2
                col(k,i)%triple(2)%ithree(3)=ipos(8,i)
                if (idp(ipos(7,i)).le.6) then
                  col(k,i)%couple(1)%itwo(1)=ipos(3,i)
                  col(k,i)%couple(1)%itwo(2)=ipos(4,i)
                endif
              else
                iactive(k,i)=0
              endif
            elseif (k.eq.6) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-2
              if (idp(ipos(5,i)).le.6.and.idp(ipos(7,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(5,i)
                col(k,i)%triple(1)%ithree(2)=2
                col(k,i)%triple(1)%ithree(3)=ipos(6,i)
                col(k,i)%triple(2)%ithree(1)=ipos(7,i)
                col(k,i)%triple(2)%ithree(2)=1
                col(k,i)%triple(2)%ithree(3)=ipos(8,i)
                if (idp(ipos(7,i)).le.6) then
                  col(k,i)%couple(1)%itwo(1)=ipos(3,i)
                  col(k,i)%couple(1)%itwo(2)=ipos(4,i)
                endif
              else
                iactive(k,i)=0
              endif
  
            elseif (k.eq.7) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-1
             if (idp(ipos(3,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(3,i)
                col(k,i)%triple(1)%ithree(2)=1
                col(k,i)%triple(1)%ithree(3)=21
                col(k,i)%triple(2)%ithree(1)=21
                col(k,i)%triple(2)%ithree(2)=2
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
            elseif (k.eq.8) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-1
              if (idp(ipos(3,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(3,i)
                col(k,i)%triple(1)%ithree(2)=2
                col(k,i)%triple(1)%ithree(3)=21
                col(k,i)%triple(2)%ithree(1)=21
                col(k,i)%triple(2)%ithree(2)=1
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
            elseif (k.eq.9) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-1
              if (idp(ipos(5,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(5,i)
                col(k,i)%triple(1)%ithree(2)=1
                col(k,i)%triple(1)%ithree(3)=21
                col(k,i)%triple(2)%ithree(1)=21
                col(k,i)%triple(2)%ithree(2)=2
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
            elseif (k.eq.10) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-1
              if (idp(ipos(5,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(5,i)
                col(k,i)%triple(1)%ithree(2)=2
                col(k,i)%triple(1)%ithree(3)=21
                col(k,i)%triple(2)%ithree(1)=21
                col(k,i)%triple(2)%ithree(2)=1
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
            elseif (k.eq.11) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-1
              if (idp(ipos(7,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(7,i)
                col(k,i)%triple(1)%ithree(2)=1
                col(k,i)%triple(1)%ithree(3)=21
                col(k,i)%triple(2)%ithree(1)=21
                col(k,i)%triple(2)%ithree(2)=2
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
            elseif (k.eq.12) then
              col(k,i)%ntriple=2
              col(k,i)%ncouple=nhadcoup-1
              if (idp(ipos(7,i)).le.6) then
                iactive(k,i)=1
                col(k,i)%triple(1)%ithree(1)=ipos(7,i)
                col(k,i)%triple(1)%ithree(2)=2
                col(k,i)%triple(1)%ithree(3)=21
                col(k,i)%triple(2)%ithree(1)=21
                col(k,i)%triple(2)%ithree(2)=1
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
            endif
  
          enddo
        enddo

 
* determine color factors                                               
*                                                                       
  
        kmax=12
        lim=namp2*kmax
        do m1=0,lim-1
          do m2=m1,lim-1
            ki=mod(m1,kmax)+1
            i=(m1+1-ki)/kmax +1
            kj=mod(m2,kmax)+1
            j=(m2+1-kj)/kmax +1
  
* fill color structure                                                  
            if (iactive(kj,j).eq.1.and.iactive(ki,i).eq.1) then
               call coleval(col(kj,j),col(ki,i),rescolfact)
              factcol(ki,i,kj,j)=rescolfact
            else
              factcol(ki,i,kj,j)=0.d0
            endif
            factcol(kj,j,ki,i)=factcol(ki,i,kj,j)
          enddo
        enddo

  
        ifirst=0
      endif   !(ionesh.eq.1.or.ifirst.eq.1.and.ionesh.eq.0)
  
* set to zero to avoid previous unwanted components                     
  
      do i2=1,2
        do i1=1,2
          do i8=1,2
            do i7=1,2
              do i6=1,2
                do i5=1,2
                  do i4=1,2
                    do i3=1,2
                      do i=1,namp2
                        do k=1,12
                          cres(i1,i2)%q(i3,i4,i5,i6,i7,i8)%col(k,i)
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
*  see if there are b-quarks around (since there are no top quarks in 
* the final state we need a pair 5 -5, so we can test just particles
      lbquark =.false.
      do i=3,7,2
	if (idp(i).eq.5) then
	 lbquark=.true.
       endif   
      enddo
*diogoend

* b-quark present OR flag mass=1     
      if (lbquark.or.i_massive.eq.1) then
      
      do i=1,n1z2w
        call  ggzww(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &       idp(ipos(8,i)),creszww)

        do k=1,12
          if(iactive(k,i).eq.1)then
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
                    cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &                   ii(6),ii(7),ii(8))%col(k,i)=
     &                   creszww(i1,i2,i3,i4,k)*isig(i)
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo

      enddo
  
  
      do i=n1z2w+1,namp2
  
          call ggzzz(p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres3z)

        do k=1,12
          if(iactive(k,i).eq.1)then
            do i8=1,2
              do i7=1,2
                do i6=1,2
                  do i5=1,2
                    do i4=1,2
                      do i3=1,2
                        do i2=1,2
                          do i1=1,2
                            ii(ipos(1,i))=i1
                            ii(ipos(2,i))=i2
                            ii(ipos(3,i))=i3
                            ii(ipos(4,i))=i4
                            ii(ipos(5,i))=i5
                            ii(ipos(6,i))=i6
                            ii(ipos(7,i))=i7
                            ii(ipos(8,i))=i8
                            cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &                           ii(6),ii(7),ii(8))%col(k,i)=
     &                           cres3z(i3,i4,i5,i6,i7,i8)%caux(i1,i2,k)
     &                           *isig(i)
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo

      enddo
  
  
  
* compute res                                                           
      res=0.d0  

* first consider diagonal terms
* compute amp2(k,i) (k=kmin,12 i=1,namp2) summed over all helicity 
* configurations

      do i=1,namp2
        amp2(i)=0.d0
        do k=1,12
          amp2k(k,i)=0.d0
          if(iactive(k,i).eq.1)then
            do i2=1,2
              do i1=1,2
                do i8=1,2
                do i7=1,2
                do i6=1,2
                do i5=1,2
                do i4=1,2
                do i3=1,2
                  amp2k(k,i)=amp2k(k,i)+
     &             cres(i1,i2)%q(i3,i4,i5,i6,i7,i8)%col(k,i)* 
     &             conjg(cres(i1,i2)%q(i3,i4,i5,i6,i7,i8)%col(k,i))*
     &             factcol(k,i,k,i)
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
      enddo

      do i=1,namp2
        res=res+amp2(i)
      enddo

* then consider interference terms

      kmax=12
      lim=namp2*kmax
      do m1=0,lim-1
        ki=mod(m1,kmax)+1
        i=(m1+1-ki)/kmax +1
        if(iactive(ki,i).eq.1)then
          do m2=m1+1,lim-1
            kj=mod(m2,kmax)+1
            j=(m2+1-kj)/kmax +1
            if(iactive(kj,j).eq.1)then
              do i2=1,2
              do i1=1,2
              do i8=1,2
              do i7=1,2
              do i6=1,2
              do i5=1,2
              do i4=1,2
              do i3=1,2
                res = res +
     &               2.d0*cres(i1,i2)%q(i3,i4,i5,i6,i7,i8)%col(ki,i)*
     &               conjg(cres(i1,i2)%q(i3,i4,i5,i6,i7,i8)%col(kj,j))
     &               *factcol(ki,i,kj,j)  
              enddo
              enddo
              enddo
              enddo
              enddo
              enddo
              enddo
              enddo
            endif
          enddo
        endif
      enddo


*b-quark absent AND imass=0
      else
      
      res=0.d0
      
      do i=1,n1z2w      
      amp2(i)=0.d0
      
        call  ggzww_massless(
     &         p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),creszww_massless)

        do k=1,12
	amp2k(k,i)=0.d0
	
          if(iactive(k,i).eq.1)then
	   do i3=1,2
	   do i2=1,2
	   do i1=1,2
	     ii(ipos(1,i))=i1
	     ii(ipos(2,i))=i2
	     ii(ipos(3,i))=i3
	     ii(ipos(4,i))=i3
	     ii(ipos(5,i))=2
	     ii(ipos(6,i))=2
	     ii(ipos(7,i))=2
	     ii(ipos(8,i))=2
	     cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &		  ii(6),ii(7),ii(8))%col(k,i)=
     &		  creszww_massless(i1,i2,i3,k)*isig(i)
     
* diagonal terms
            amp2k(k,i)=amp2k(k,i)+
     &          cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(k,i)* 
     &          conjg(cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(k,i))*
     &          factcol(k,i,k,i)

* interference terms 

* Same master amp but different color            
	     do kj=1,k-1
	       if (iactive(kj,i).eq.1) then
	          res = res +
     &               2.d0*cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(k,i)*
     &               conjg(cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(kj,i))
     &               *factcol(k,i,kj,i)                  
	       endif	    
	     enddo

* Other master amp	     
             do j=1,i-1             
	     do kj=1,12
	       if (iactive(kj,j).eq.1) then	       
	          res = res +
     &               2.d0*cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(k,i)*
     &               conjg(cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(kj,j))
     &               *factcol(k,i,kj,j)                  
	       endif	    
	     enddo
	     enddo
     
           enddo
           enddo
           enddo
          
	  endif
	  amp2(i)=amp2(i)+amp2k(k,i)  
	  
        enddo
	
	
*        print*,'diagonal(',i,')=',amp2(i)
	
        res=res+amp2(i)
	
*	print*, 'diag+interf(',i,')=',res
	
      enddo
  
      
  
      do i=n1z2w+1,namp2
      amp2(i)=0.d0
      
          call ggzzz_massless(
     &         p(0,ipos(1,i)),p(0,ipos(2,i)),p(0,ipos(3,i)),
     &         p(0,ipos(4,i)),p(0,ipos(5,i)),p(0,ipos(6,i)),
     &         p(0,ipos(7,i)),p(0,ipos(8,i)),idp(ipos(1,i)),
     &         idp(ipos(2,i)),idp(ipos(3,i)),idp(ipos(4,i)),
     &         idp(ipos(5,i)),idp(ipos(6,i)),idp(ipos(7,i)),
     &         idp(ipos(8,i)),cres3z_massless)

        do k=1,12
	amp2k(k,i)=0.d0
	
          if(iactive(k,i).eq.1)then
            do i7=1,2
            do i5=1,2
            do i3=1,2
            do i2=1,2
            do i1=1,2
              ii(ipos(1,i))=i1
              ii(ipos(2,i))=i2
              ii(ipos(3,i))=i3
              ii(ipos(4,i))=i3
              ii(ipos(5,i))=i5
              ii(ipos(6,i))=i5
              ii(ipos(7,i))=i7
              ii(ipos(8,i))=i7
              cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &          ii(6),ii(7),ii(8))%col(k,i)=
     &          cres3z_massless(i1,i2,i3,i5,i7,k)
     &          *isig(i)

* diagonal terms
            amp2k(k,i)=amp2k(k,i)+
     &          cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(k,i)* 
     &          conjg(cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(k,i))*
     &          factcol(k,i,k,i)

* interference terms 

* Same master amp but different color            
	     do kj=1,k-1
	       if (iactive(kj,i).eq.1) then
	          res = res +
     &               2.d0*cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(k,i)*
     &               conjg(cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(kj,i))
     &               *factcol(k,i,kj,i)                  
	       endif	    
	     enddo

* Other master amp
             do j=1,i-1             
	     do kj=1,12
	       if (iactive(kj,j).eq.1) then
	          res = res +
     &               2.d0*cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(k,i)*
     &               conjg(cres(ii(1),ii(2))%q(ii(3),ii(4),ii(5),
     &               ii(6),ii(7),ii(8))%col(kj,j))
     &               *factcol(k,i,kj,j)
	       endif	    
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

      endif  !if b-quark        

*      print*,'res=',res
       
c giuseppe ILC 25/06/2007
* spin average 1/4
* color average in case of initial state partons only
      if(i_coll.eq.1 .or. i_coll.eq.2)then ! hadron colliders
        if(ingl.eq.0)then
          rinvclav=9.d0
        elseif(ingl.eq.1)then
          rinvclav=24.d0
        elseif(ingl.eq.2)then
          rinvclav=64.d0
        endif
      elseif(i_coll.eq.3)then   ! e+e- collider
        rinvclav=1.d0
      endif
c end giuseppe ILC 25/06/2007


*coupling g4e8
      res=res/4.d0/rinvclav/symfact*g4e8
      
* extract the configuration ichosen for color flow                      
  
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
  
      ichosen=i
  
  
* extract the configuration kchosen for color flow                      
  
      amp2tot=0.d0
      do k=1,12
	amp2tot=amp2tot+amp2k(k,i)
      enddo
  
      do k=1,12
	amp2k(k,i)=amp2k(k,i)/amp2tot
      enddo
      k=1
      am2=amp2k(1,i)
      zzran=ran2(idum)
      do while (zzran.gt.am2.and.k.lt.12)
	k=k+1
	am2=am2+amp2k(k,i)
      enddo
  
      kchosen=k
      
* compute iord as the inverse of ipos(i,nchosen)                        
*                                                                       
* old ampem                                                             
*      do i=1,8                                                         
*        iord(ipos(i,nchosen))=i                                        
*      enddo                                                            
*  
      if (kchosen.le.6)then
	iord(col(kchosen,ichosen)%triple(1)%ithree(1))=1
	iord(col(kchosen,ichosen)%triple(1)%ithree(2))=2
	iord(col(kchosen,ichosen)%triple(1)%ithree(3))=3
	iord(col(kchosen,ichosen)%triple(2)%ithree(1))=4
	iord(col(kchosen,ichosen)%triple(2)%ithree(2))=5
	iord(col(kchosen,ichosen)%triple(2)%ithree(3))=6
  
	do i=1,col(kchosen,ichosen)%ncouple
	  iord(col(kchosen,ichosen)%couple(i)%itwo(1))=6+2*i-1
	  iord(col(kchosen,ichosen)%couple(i)%itwo(2))=6+2*i
	enddo
  
      else
	iord(col(kchosen,ichosen)%triple(1)%ithree(1))=1
	iord(col(kchosen,ichosen)%triple(1)%ithree(2))=2
	iord(col(kchosen,ichosen)%triple(2)%ithree(2))=3
	iord(col(kchosen,ichosen)%triple(2)%ithree(3))=4
  
	do i=1,col(kchosen,ichosen)%ncouple
	  iord(col(kchosen,ichosen)%couple(i)%itwo(1))=4+2*i-1
	  iord(col(kchosen,ichosen)%couple(i)%itwo(2))=4+2*i
	enddo
      endif
  
  
***************                                                         
*  
* fill the non hadronic positions if needed                             
      icount=2*nhadcoup+3
      do j=1,8
	if (abs(idp(ipos(j,ichosen))).gt.6.and.
     &       idp(ipos(j,ichosen)).ne.21) then
	  iord(ipos(j,ichosen))=icount
	  icount=icount+1
	endif
      enddo
   
 
      return
      end
  
