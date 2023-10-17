
      integer function isignperm(igo,igluons)
      implicit real*8 (a-h,o-z)
      integer igo(8),ignew(8)

* this function determines whether the permutation igo is odd
*  or even with respect to the fundamental 1 2 3 4 5 6 7 8
* If odd isignper is -1, if even it is 1
      

      do k=1,8
        ignew(k)=igo(k)
      enddo

ctest      print*,ignew

      if (igluons.eq.1) then
        jstart=3
      else
        jstart=1
      endif

      ntot=0
      do j=jstart,7

        i=j
        do while (ignew(i).ne.j)
          i=i+1
        enddo

ctest        print*, 'j',j,'i',i

        if (i.gt.j) then
        
          n=i-j
          ntot=ntot+n

          igsalt=ignew(i)
        
          do k=i-1,j,-1
            ignew(k+1)=ignew(k)
          enddo
          ignew(j)=igsalt

ctest          print*, ignew

        elseif (i.lt.j) then
          print*, 'ERROR in Reorder'
          stop
        endif

            
ctest        print*,n,ntot
ctest        print*, '  '

      enddo

      isignperm =(-1)** mod(ntot,2)
        
      return
      end
