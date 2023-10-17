      subroutine perm(idp, igrp, igo, isign,
     &     n4w,n2z2w,n4z, ig4w, ig2z2w,ig4z,isign4w,isign2z2w,isign4z)
      implicit real*8 (a-h,o-z)
      dimension idp(8), igrp(4),igo(8,4),isign(4)
      parameter (maxperm4w=4,maxperm2z2w=18*2,max4z=24)
      dimension ig4w(8,maxperm4w), ig2z2w(8,maxperm2z2w),ig4z(8,max4z)
      dimension isign4w(maxperm4w),isign2z2w(maxperm2z2w),isign4z(max4z)
      dimension numexc4w(2,2),numexc2z2w(2,9),numexc4z(2,6)
      DATA ifirst/1/


      if (ifirst.eq.1) then
*possible exchanges for 4 W's
        numexc4w(1,1)=1
        numexc4w(2,1)=5
        numexc4w(1,2)=3
        numexc4w(2,2)=7

*possible exchanges for 2z2w

* possible exchange between 2 z's
        numexc2z2w(1,1)=1
        numexc2z2w(2,1)=3

* possible exch on the previous
        numexc2z2w(1,2)=1
        numexc2z2w(2,2)=5
        numexc2z2w(1,3)=3
        numexc2z2w(2,3)=5

* possible exch on the previous
        numexc2z2w(1,4)=1
        numexc2z2w(2,4)=7
        numexc2z2w(1,5)=3
        numexc2z2w(2,5)=7
        
* possible exch on the previous
        numexc2z2w(1,6)=2
        numexc2z2w(2,6)=6
        numexc2z2w(1,7)=4
        numexc2z2w(2,7)=6
        
* possible exch on the previous
        numexc2z2w(1,8)=2
        numexc2z2w(2,8)=8
        numexc2z2w(1,9)=4
        numexc2z2w(2,9)=8

*possible exchanges for 4 Z's

* exchange among 1 and others 
        numexc4z(1,1)=1
        numexc4z(2,1)=3
        numexc4z(1,2)=1
        numexc4z(2,2)=5
        numexc4z(1,3)=1
        numexc4z(2,3)=7

* exchange among 3 and higher to be performed on the previous 
        numexc4z(1,4)=3
        numexc4z(2,4)=5
        numexc4z(1,5)=3
        numexc4z(2,5)=7

* exchange among 5 and 7 to be performed on the previous 
        numexc4z(1,6)=5
        numexc4z(2,6)=7

      ifirst = 0
      endif

      n4w=0
      n2z2w=0
      n4z=0



      if (igrp(1).eq.1) then
        n4w=n4w+1
        do i=1,8 
          ig4w(i,1)=igo(i,1)
        enddo
        isign4w(1)=isign(1)



        do k=1,2
          if (idp(igo(numexc4w(1,k),1)).eq.idp(igo(numexc4w(2,k),1))) 
     &         then
            n4wold=n4w
            do j=1,n4wold
              n4w=n4w+1
              do i=1,8
                ig4w(i,n4wold+j)=ig4w(i,j)
              end do
              ig4w(numexc4w(1,k),n4wold+j)=ig4w(numexc4w(2,k),j)
              ig4w(numexc4w(2,k),n4wold+j)=ig4w(numexc4w(1,k),j)
              isign4w(n4wold+j)=-isign4w(j)            
            enddo
          endif
        enddo

c        do i=1,n4w
c          print*, (idp(ig4w(j,i)),j=1,8)
c          print*, (ig4w(j,i),j=1,8)
c          ilsegno=isignperm(ig4w(1,i))
c          print*, 'isign4w,ilsegno',isign4w(i),ilsegno
c          print*, '   '
c        enddo
      
      endif

      do nn=1,2   ! number of the 2z2w channel
        if (nn.eq.1.and.igrp(2).eq.1.or.nn.eq.2.and.igrp(3).eq.1) then
          n2z2w=n2z2w+1
          n2z2w_loc=1
          do i=1,8 
            ig2z2w(i,n2z2w)=igo(i,nn+1)
          enddo
          isign2z2w(n2z2w)=isign(nn+1)

          if (idp(igo(numexc2z2w(1,1),nn+1)).eq.
     &         idp(igo(numexc2z2w(2,1),nn+1))) then
            n2z2w=n2z2w+1
            n2z2w_loc=2
            do i=1,8
              ig2z2w(i,n2z2w)=ig2z2w(i,n2z2w-1)
            end do
            ig2z2w(numexc2z2w(1,1),n2z2w)=
     &           ig2z2w(numexc2z2w(2,1),n2z2w-1)
            ig2z2w(numexc2z2w(2,1),n2z2w)=
     &           ig2z2w(numexc2z2w(1,1),n2z2w-1)
            isign2z2w(n2z2w)=-isign2z2w(n2z2w-1)
          endif


          do k=2,8,2

            n2z2wold=n2z2w_loc
            do m=0,1
              if (idp(igo(numexc2z2w(1,k+m),nn+1)).eq.
     &             idp(igo(numexc2z2w(2,k+m),nn+1))) then
                do j=1,n2z2wold
                  n2z2w_loc=n2z2w_loc+1
                  n2z2w=n2z2w+1
c                  do i=1,8
c                    ig2z2w(i,n2z2w)=ig2z2w(i,n2z2w-n2z2wold)
c                  end do
c                  ig2z2w(numexc2z2w(1,k+m),n2z2w)=
c     &                 ig2z2w(numexc2z2w(2,k+m),n2z2w-n2z2wold)
c                  ig2z2w(numexc2z2w(2,k+m),n2z2w)=
c     &                 ig2z2w(numexc2z2w(1,k+m),n2z2w-n2z2wold)
c                  isign2z2w(n2z2w)=-isign2z2w(n2z2w-n2z2wold)
                  do i=1,8
                    ig2z2w(i,n2z2w)=ig2z2w(i,j)
                  end do
                  ig2z2w(numexc2z2w(1,k+m),n2z2w)=
     &                 ig2z2w(numexc2z2w(2,k+m),j)
                  ig2z2w(numexc2z2w(2,k+m),n2z2w)=
     &                 ig2z2w(numexc2z2w(1,k+m),j)
                  isign2z2w(n2z2w)=-isign2z2w(j)
                enddo
              endif
            enddo

          enddo

        endif

      enddo



      if (igrp(4).eq.1) then
        n4z=n4z+1
        do i=1,8 
          ig4z(i,1)=igo(i,4)
        enddo
        isign4z(1)=isign(4)


        do m=3,1,-1
          if (m.eq.3) iprec=0
          if (m.eq.2) iprec=3
          if (m.eq.1) iprec=5

          n4zold=n4z
          do k=1,m
            if (idp(igo(numexc4z(1,k+iprec),4)).eq.
     &           idp(igo(numexc4z(2,k+iprec),4))) then
              do j=1,n4zold
                n4z=n4z+1
                do i=1,8
                  ig4z(i,n4z)=ig4z(i,j)
                enddo
                ig4z(numexc4z(1,k+iprec),n4z)=
     &               ig4z(numexc4z(2,k+iprec),j)
                ig4z(numexc4z(2,k+iprec),n4z)=
     &               ig4z(numexc4z(1,k+iprec),j)
                isign4z(n4z)=-isign4z(j)

c                print*, 'm,k,k+iprec,n4z'
c                print*,m,k,k+iprec,n4z
c                print*, (ig4z(jj,n4z),jj=1,8)
              enddo  
              
            endif
          enddo 
        enddo

      endif



      return
      end


      
