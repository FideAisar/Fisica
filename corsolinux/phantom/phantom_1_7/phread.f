      subroutine convert(ii,charint1,charint2,charint3)

      character*1 charint1
      character*2 charint2
      character*3 charint3

      if (ii.le.9) then
        write (charint1,'(i1)') ii
      elseif (ii.le.99) then
        write (charint2,'(i2)') ii
      elseif (ii.le.999) then
        write (charint3,'(i3)') ii
      endif
      
      return
      end

      subroutine rread(name,var,num)
      implicit real*8 (a-h, o-z)
      character *(*) name
      character *81 line
      character *80 pluto
      dimension var(num)
      common/abread/pluto
      common/phwrite/iwritevar
      
      open (unit=20,file=pluto, status='old')
      le=len(name)
    1 read (20,'(a)',end=100) line
      
      if (line(1:40).ne.'                                        '.or.
     &   line(41:80).ne.'                                        ') 
     &    then
        if (line(1:1).ne.'C'.and.line(1:1).ne.'*') then
          i=0
          do while (line(i+1:i+1).eq.' '.and.(i+le).lt.80) 
             i=i+1
          end do
          if (line(i+1:i+le).eq.name.and.line(i+le+1:i+le+1).eq.' ') 
     &                     then
            n=0
            l=i+le
    2       do while (l.lt.80.and.n.lt.num)
              do while (line(l+1:l+1).eq.' '.and.l.lt.81) 
                l=l+1
              enddo 
              if (l.le.80) then
                 n=n+1
                 read (line(l:), *,end=200,err=100) var(n)
              endif
              do while (line(l+1:l+1).ne.' '.and.l.le.79)
                l=l+1
              enddo
  200         continue
            enddo
            
            if (n.lt.num) then
              read (20,'(a)',end=100) line 
              l=1
              goto 2
            endif
            close (unit=20)
            if (iwritevar.ne.1) print*, name//'   ', var
            return 
          endif
        endif
      endif

      goto 1
  100 print*, '                     '
      print*, '                     '
      print*, '!!!!! INPUT ERROR AT THE VARIABLE '//name//' !!!!!!'
      print*, '                     '
      print*, '                     '
      stop
      end


      subroutine iread(name,ivar,num)
      implicit real*8 (a-h, o-z)
      character *(*) name
      character *81 line
      character *80 pluto
      dimension ivar(num)
      common/abread/pluto
      common/phwrite/iwritevar

      open (unit=20,file=pluto, status='old')
      le=len(name)
    1 read (20,'(a)',end=100) line
      
      if (line(1:40).ne.'                                        '.or.
     &   line(41:80).ne.'                                        ') 
     &    then
        if (line(1:1).ne.'C'.and.line(1:1).ne.'*') then
          i=0
          do while (line(i+1:i+1).eq.' '.and.(i+le).lt.80) 
             i=i+1
          end do
          if (line(i+1:i+le).eq.name.and.line(i+le+1:i+le+1).eq.' ')
     &                     then
            n=0
            l=i+le
    2       do while (l.lt.80.and.n.lt.num)
              do while (line(l+1:l+1).eq.' '.and.l.lt.81) 
                l=l+1
              enddo 
              if (l.le.80) then
                 n=n+1
                 read (line(l:), *,end=200,err=100) ivar(n)
              endif
              do while (line(l+1:l+1).ne.' '.and.l.le.79)
                l=l+1
              enddo
  200         continue
            enddo
            
            if (n.lt.num) then
              read (20,'(a)',end=100) line 
              l=1
              goto 2
            endif
            close (unit=20)
            if (iwritevar.ne.1) print*, name//'   ', ivar
            return 
          endif
        endif
      endif

      goto 1
  100 print*, '                     '
      print*, '                     '
      print*, '!!!!! INPUT ERROR AT THE VARIABLE '//name//' !!!!!!'
      print*, '                     '
      print*, '                     '
      stop
      end


conesh 

      subroutine ireadnoclose(name,ivar,num)
c  this routine is completely similar to iread, but it does not close
c  the input file  to which it is assigned unit 30 and it does not 
c  print the name and the values of the variables.
c   It can be useful if one has to read sequentially something after
c   having read some variable.

      implicit real*8 (a-h, o-z)
      character *(*) name
      character *81 line
      character *80 pluto
      dimension ivar(num)
      common/abread/pluto
      common/phwrite/iwritevar
      

      open (unit=30,file=pluto, status='old')
      le=len(name)
    1 read (30,'(a)',end=100) line
      
      if (line(1:40).ne.'                                        '.or.
     &   line(41:80).ne.'                                        ') 
     &    then
        if (line(1:1).ne.'C'.and.line(1:1).ne.'*') then
          i=0
          do while (line(i+1:i+1).eq.' '.and.(i+le).lt.80) 
             i=i+1
          end do
          if (line(i+1:i+le).eq.name.and.line(i+le+1:i+le+1).eq.' ')
     &                     then
            n=0
            l=i+le
    2       do while (l.lt.80.and.n.lt.num)
              do while (line(l+1:l+1).eq.' '.and.l.lt.81) 
                l=l+1
              enddo 
              if (l.le.80) then
                 n=n+1
                 read (line(l:), *,end=200,err=100) ivar(n)
              endif
              do while (line(l+1:l+1).ne.' '.and.l.le.79)
                l=l+1
              enddo
  200         continue
            enddo
            
            if (n.lt.num) then
              read (30,'(a)',end=100) line 
              l=1
              goto 2
            endif
c            close (unit=30)
            if (iwritevar.ne.1) print*, name//'   ', ivar
            return 
          endif
        endif
      endif

      goto 1
  100 print*, '                     '
      print*, '                     '
      print*, '!!!!! INPUT ERROR AT THE VARIABLE '//name//' !!!!!!'
      print*, '                     '
      print*, '                     '
      stop
      end

coneshend



      subroutine cread4(name,cvar,num)
      implicit real*8 (a-h, o-z)
      character *(*) name
      character *81 line
      character *80 pluto
      character *4 cvar
      dimension cvar(num)
      common/abread/pluto
      common/phwrite/iwritevar

      open (unit=20,file=pluto, status='old')
      le=len(name)
    1 read (20,'(a)',end=100) line
      
      if (line(1:40).ne.'                                        '.or.
     &   line(41:80).ne.'                                        ') 
     &    then
        if (line(1:1).ne.'C'.and.line(1:1).ne.'*') then
          i=0
          do while (line(i+1:i+1).eq.' '.and.(i+le).lt.80) 
             i=i+1
          end do
          if (line(i+1:i+le).eq.name.and.line(i+le+1:i+le+1).eq.' ')
     &                     then
            n=0
            l=i+le
    2       do while (l.lt.80.and.n.lt.num)
              do while (line(l+1:l+1).eq.' '.and.l.lt.81) 
                l=l+1
              enddo 
              if (l.le.80) then
                 n=n+1
                 read (line(l+1:),'(a)',end=200,err=100) cvar(n)
              endif
              do while (line(l+1:l+1).ne.' '.and.l.le.79)
                l=l+1
              enddo
  200         continue
            enddo
            
            if (n.lt.num) then
              read (20,'(a)',end=100) line 
              l=1
              goto 2
            endif
            close (unit=20)
c            do i=1,num
            if (iwritevar.ne.1) print*, name//'   ', cvar 
c            enddo
            return 
          endif
        endif
      endif

      goto 1
  100 print*, '                     '
      print*, '                     '
      print*, '!!!!! INPUT ERROR AT THE VARIABLE '//name//' !!!!!!'
      print*, '                     '
      print*, '                     '
      stop
      end

      subroutine cread(name,cvar)
      implicit real*8 (a-h, o-z)
      character *(*) name
      character *300 line
      character *80 pluto
      character *(*) cvar
      common/abread/pluto
      common/phwrite/iwritevar

      open (unit=20,file=pluto, status='old')
      le=len(name)
    1 read (20,'(a)',end=100) line
      if (line(1:40).ne.'                                        '.or.
     &   line(41:80).ne.'                                        ') 
     &    then
        if (line(1:1).ne.'C'.and.line(1:1).ne.'*') then
          i=0
          do while (line(i+1:i+1).eq.' '.and.(i+le).lt.80) 
             i=i+1
          end do
          if (line(i+1:i+le).eq.name.and.line(i+le+1:i+le+1).eq.' ')
     &                     then
            l=i+le
              do while (line(l+1:l+1).eq.' '.and.l.lt.201) 
                l=l+1
              enddo 
             if (l.le.200) then
                 read (line(l+1:),'(a)',end=200,err=100) cvar
              endif
  200         continue
            
            close (unit=20)
            if (iwritevar.ne.1) print*, name//'   ', cvar 
            return 
          endif
        endif
      endif

      goto 1
  100 print*, '                     '
      print*, '                     '
      print*, '!!!!! INPUT ERROR AT THE VARIABLE '//name//' !!!!!!'
      print*, '                     '
      print*, '                     '
      stop
      end

