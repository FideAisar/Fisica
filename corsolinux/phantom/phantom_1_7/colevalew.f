      double precision function rcolor (iiord1,iiord2,idp)

* This function evaluates the color factor of two amplitudes
* with eight external fermions, given respectively in the order 
* iiord1 and iiord2
* This means that, given the order 1 2 3 4 5 6 7 8 to which
* corresponds  particle identity idp(1) ...idp(8) and the momenta
* p1,p2,...p8,
* the first amplitude will have the order
*     iiord1(1), iiord1(2)....iiord1(8)
* and identities
*     idp(iiord1(1)), idp(iiord1(2))....idp(iiord1(8))    
* 
* This means that the particle iiord1(1) will form a fermion line
* with iiord1(2), iiord1(3) with iiord1(4) ...
* Analogous considerations for the second amplitude

* It is assumed that particles are given always in the order
* particle antiparticle, particle antiparticle, ....
* and the various couples form a W+, a W- or a Z
 

      implicit real *8 (a-h,o-z)
      
      dimension iiord1(8), iiord2(8), idp(8), iord1(8), iord2(8),icol(8)

* copy the non quarks in iord1 and iord2 and count them

      ncolor=0
      do i=1,8
        if (abs(idp(iiord1(i))).le.5) then
          ncolor=ncolor+1
          iord1(ncolor)=iiord1(i) 
        endif
      enddo

      ncolor2=0
      do i=1,8
        if (abs(idp(iiord2(i))).le.5) then
          ncolor2=ncolor2+1
          iord2(ncolor2)=iiord2(i) 
        endif
      enddo

      if (ncolor.ne.ncolor2) then
        print*, 'ERROR in color assignment'
        stop
      endif  

ctest      print*, 'ncolor', ncolor

ctest      print*, 'iord1', iord1
ctest      print*, 'iord2', iord2
      

      do j=2,ncolor,2

        icol(j)=iord1(j)

        i=1

        do while (iord2(i).ne.iord1(j-1))
          i=i+2
        enddo
      
        icol(j-1)=iord2(i+1)

      enddo

ctest      print*, icol

      n=ncolor
      
      rcolor=1.d0
      
      do while (n.gt.0)
      
        if (icol(n-1).eq.icol(n)) then
          rcolor=rcolor*3.d0
        else
          i=2
          do while (icol(i).ne.icol(n-1))
            i=i+2
          enddo
          icol(i)=icol(n)
        endif
        n=n-2
ctest        print*, 'rcolor',rcolor
ctest        print*, 'n',n
ctest        print*, 'icol',icol
      enddo

      return
      end
