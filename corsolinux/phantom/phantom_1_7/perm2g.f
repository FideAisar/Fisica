      subroutine perm2g(idp, igrp, igo, isign,
     &       n1z2w,n3z, ig1z2w, ig3z,isign1z2w,isign3z)

      implicit real*8 (a-h,o-z)
      dimension idp(8), igrp(2),igo(8,2),isign(2)
      parameter (maxperm1z2w=4,max3z=6)
      dimension  ig1z2w(8,maxperm1z2w),ig3z(8,max3z)
      dimension  isign1z2w(maxperm1z2w),isign3z(max3z)
      dimension  numexc1z2w(2,4),numexc3z(2,3)
      DATA ifirst/1/

      if (ifirst.eq.1) then

*possible exchanges for 1z2w
        numexc1z2w(1,1)=3
        numexc1z2w(2,1)=5
        numexc1z2w(1,2)=3
        numexc1z2w(2,2)=7
        numexc1z2w(1,3)=4
        numexc1z2w(2,3)=6
        numexc1z2w(1,4)=4
        numexc1z2w(2,4)=8

*possible exchanges for 3 Z's

* exchange among 3 and higher 
        numexc3z(1,1)=3
        numexc3z(2,1)=5
        numexc3z(1,2)=3
        numexc3z(2,2)=7

* exchange among 5 and 7 to be performed on the previous 
        numexc3z(1,3)=5
        numexc3z(2,3)=7

      ifirst = 0
      endif

      n1z2w=0
      n3z=0

 
      if (igrp(1).eq.1) then
        n1z2w=n1z2w+1
        do i=1,8 
          ig1z2w(i,1)=igo(i,1)
        enddo
        isign1z2w(1)=isign(1)

        do k=1,4
          if (idp(igo(numexc1z2w(1,k),1)).eq.
     &         idp(igo(numexc1z2w(2,k),1))) then
            n1z2wold=n1z2w
            do j=1,n1z2wold
              n1z2w=n1z2w+1
              do i=1,8
                ig1z2w(i,n1z2wold+j)=ig1z2w(i,j)
              end do
              ig1z2w(numexc1z2w(1,k),n1z2wold+j)=
     &             ig1z2w(numexc1z2w(2,k),j)
              ig1z2w(numexc1z2w(2,k),n1z2wold+j)=
     &             ig1z2w(numexc1z2w(1,k),j)
              isign1z2w(n1z2wold+j)=-isign1z2w(j)            
            enddo
          endif
        enddo
      endif

      if (igrp(2).eq.1) then
        n3z=n3z+1
        do i=1,8 
          ig3z(i,1)=igo(i,2)
        enddo
        isign3z(1)=isign(2)

        do m=2,1,-1
          if (m.eq.2) iprec=0
          if (m.eq.1) iprec=2

          n3zold=n3z
          do k=1,m
            if (idp(igo(numexc3z(1,k+iprec),2)).eq.
     &           idp(igo(numexc3z(2,k+iprec),2))) then
              do j=1,n3zold
                n3z=n3z+1
                do i=1,8
                  ig3z(i,n3z)=ig3z(i,j)
                enddo
                ig3z(numexc3z(1,k+iprec),n3z)=
     &               ig3z(numexc3z(2,k+iprec),j)
                ig3z(numexc3z(2,k+iprec),n3z)=
     &               ig3z(numexc3z(1,k+iprec),j)
                isign3z(n3z)=-isign3z(j)
              enddo  
            endif
          enddo 
        enddo
      endif

      return
      end


      
