* Adds two gnuplot datafiles

      implicit real*8 (a-h,o-z)
      character * 200 filein1,filein2,fileout,pippo

      print*, 'absolute address of input file 1?'
      read (5,'(a)') filein1
      print*, 'absolute address of input file 2?'
      read (5,'(a)') filein2
      print*, 'absolute address of output file?'
      read (5,'(a)') fileout

      open (unit=21,status='old',file=filein1)     
      open (unit=22,status='old',file=filein2)    
      open (unit=23,status='new',file=fileout)

* Assume 5 lines to be skipped
      do i=1,5
        read(21,*)pippo
        read(22,*)pippo
        write(23,*)pippo(1:(index(pippo,' ')-1))
      enddo

* Assume 1000 datapoints
      do i=1,1000
        read(21,*)rx1,ry1
        read(22,*)rx2,ry2
        if(rx1.eq.rx2)then
          write(23,*)rx1,ry1+ry2
        else
          write(*,*)'x coordinates do not match:',rx1,rx2
          stop
        endif
      enddo

      close(21)
      close(22)
      close(23)

      end

