
      subroutine coleval(col1,col2,rescolfact)

      implicit real*8 (a-h,o-z)

      parameter (maxcouple=4,maxtriple=2)
      parameter (maxcoupletot=2*maxcouple+maxtriple)
      parameter (maxconfigtot=maxtriple**2)
      parameter (twoinv=1.d0/2.d0,sixinvm=-1.d0/6.d0)

*structure for the two incoming color configurations

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
      type(uno)  col1,col2

*structure for the sum of the triples

      type tre3
         integer ithree(3)
      end type

      type uno3
         integer ntriple
         type(tre)  triple(2*maxtriple)
      end type
      type(uno3)  colthree



c      type tretot
c         integer ithree(3)
c      end type

      type duetot
         integer itwo(2)
      end type

      type unotot
         real*8 fact
c         integer ncouple,ntriple
         integer ncouple
         type(duetot)  couple(maxcoupletot)
      end type
      type(unotot)  col(maxconfigtot)

*      print*, 'col1'
*      print*, col1


*      print*, 'col2'
*      print*, col2

* put together the triples changing signes to the "bar" dummy gluons
*  or quarks.  dummy gluons or quarks must have  a number ge  20

      colthree%ntriple=col1%ntriple+col2%ntriple

* the number of final configurations will be :
      if (colthree%ntriple.eq.0) then
        ntotalconfig=1
      else
        ntotalconfig=2**(colthree%ntriple/2)
      endif


* initalize col()%fact
      do n=1,ntotalconfig
        col(n)%fact=1.d0
      enddo

      do k=1,col1%ntriple
        do j=1,3
          colthree%triple(k)%ithree(j)=col1%triple(k)%ithree(j)
        enddo
      enddo
      do k=1,col2%ntriple 
* invert the bar order
        if (col2%triple(k)%ithree(3).ge.20) then
          colthree%triple(k+col1%ntriple)%ithree(1)=
     &         -col2%triple(k)%ithree(3)
        else
          colthree%triple(k+col1%ntriple)%ithree(1)=
     &         col2%triple(k)%ithree(3)
        endif
        if (col2%triple(k)%ithree(1).ge.20) then
          colthree%triple(k+col1%ntriple)%ithree(3)=
     &         -col2%triple(k)%ithree(1)
          else
          colthree%triple(k+col1%ntriple)%ithree(3)=
     &         col2%triple(k)%ithree(1)
        endif
        if (col2%triple(k)%ithree(2).ge.20) then
          colthree%triple(k+col1%ntriple)%ithree(2)=
     &         -col2%triple(k)%ithree(2)
        else
          colthree%triple(k+col1%ntriple)%ithree(2)=
     &         col2%triple(k)%ithree(2)
        endif
      enddo

*      print*, 'colthree'
*      print*, colthree


* put together the doubles in all col()% The number of col will
*  be equal to the total numer of ntriples

      col(1)%ncouple=col1%ncouple+col2%ncouple
      do k=1,col1%ncouple
        do i=1,2
          col(1)%couple(k)%itwo(i)=col1%couple(k)%itwo(i)
        enddo
      enddo

*      print*, 'first col(1)',col(1)

      do k=1,col2%ncouple
*invert the bar order
c        do i=1,2
c          col(1)%couple(k+col1%ncouple)%itwo(i)=col2%couple(k)%itwo(i)
c        enddo
        col(1)%couple(k+col1%ncouple)%itwo(1)=col2%couple(k)%itwo(2)
        col(1)%couple(k+col1%ncouple)%itwo(2)=col2%couple(k)%itwo(1)
      enddo

      nfilled=1

*      print*, 'col(1)',col(1)


* fill all other couples from the triples with the usual rule
*       1/2 X - 1/6 ||

      do nn=1,colthree%ntriple/2
        k=1
        do while (colthree%triple(k)%ithree(2).eq.0.and. 
     &       k.lt.colthree%ntriple-1) 
          k=k+1
        enddo
        i1=colthree%triple(k)%ithree(2)
        colthree%triple(k)%ithree(2)=0
        n=k+1
        do while (colthree%triple(n)%ithree(2).eq.0.or.
     &       i1.ne.colthree%triple(n)%ithree(2))
          n=n+1
        enddo
        if (n.gt.colthree%ntriple) print*, 'ERROR: n too big in coleval'
        colthree%triple(n)%ithree(2)=0

* at this point k and n should indicate the couple connected by
*  the gluon
c        print*, 'k, n'
c        print*, k, n

        do l=1,nfilled

          col(l)%couple(col(l)%ncouple+1)%itwo(1)=
     &         colthree%triple(k)%ithree(1)
          col(l)%couple(col(l)%ncouple+1)%itwo(2)=
     &         colthree%triple(n)%ithree(3)

          col(l)%couple(col(l)%ncouple+2)%itwo(1)=
     &         colthree%triple(n)%ithree(1)
          col(l)%couple(col(l)%ncouple+2)%itwo(2)=
     &         colthree%triple(k)%ithree(3)
          
          

          col(l+nfilled)%fact=col(l)%fact
          do m=1,col(l)%ncouple
            do ii=1,2
              col(l+nfilled)%couple(m)%itwo(ii)=
     &             col(l)%couple(m)%itwo(ii)
            enddo
          enddo
          col(l+nfilled)%couple(col(l)%ncouple+1)%itwo(1)=
     &         colthree%triple(k)%ithree(1)
          col(l+nfilled)%couple(col(l)%ncouple+1)%itwo(2)=
     &         colthree%triple(k)%ithree(3)

          col(l+nfilled)%couple(col(l)%ncouple+2)%itwo(1)=
     &         colthree%triple(n)%ithree(1)
          col(l+nfilled)%couple(col(l)%ncouple+2)%itwo(2)=
     &         colthree%triple(n)%ithree(3)
 
          col(l)%ncouple=col(l)%ncouple+2
          col(l+nfilled)%ncouple=col(l)%ncouple

* adjust the the factors
          col(l)%fact=col(l)%fact*twoinv
          col(l+nfilled)%fact=col(l+nfilled)%fact*sixinvm
            
        enddo

        nfilled=2*nfilled
      enddo


*      print*, 'col'

*      print*, col

*      print*,'nfilled',nfilled

      rescolfact=0.d0
      do l=1,nfilled
        n=col(l)%ncouple
*        print*,'n',n
        do while (n.gt.0)
          if (col(l)%couple(n)%itwo(2).eq.col(l)%couple(n)%itwo(1)) then
            col(l)%fact=col(l)%fact*3.d0
          else
            i=1
            do while (col(l)%couple(i)%itwo(2).ne.
     &           col(l)%couple(n)%itwo(1))
              i=i+1
*              print *, i
            enddo
            col(l)%couple(i)%itwo(2)=col(l)%couple(n)%itwo(2)
          endif
          n=n-1

        enddo

c        print*,'col(l)%fact',col(l)%fact
        rescolfact=rescolfact+col(l)%fact

      enddo


      return
      end



