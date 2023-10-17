* This program computes the total integral of a generation
* The input file can be prepared from a directory, in which all partial
* generations gen1 gen2 .... are stored, with the command:
*     grep -A 1 total\ integral gen*/run.o*  > res 
* Then run the program with the command .../gentotint.exe > result and you will
*  get in the file result the cross sections evaluated by all partial generations
*  and the average which is the total cross section with its extimated error.


      implicit real*8 (a-h,o-z)
      character * 80 pippo
      character * 80 sigma
      character * 80 error

c      print*, 'nome del file?'
c      read (5,'(a)') pippo
      pippo='res'
      npi=index(pippo,' ')-1

      open (unit=20,status='old',file=pippo(1:npi))     

      totint=0.d0
      toterr=0.d0
      nentries=0

    1 read (20,'(a)',err=101,end=101) sigma
      read (20,'(a)',err=101,end=101) error
      npi=index(sigma,'=')     
      read (sigma(npi+1:),*) x
      npi=index(error,'=')     
      read (error(npi+1:),*) y
      print*,x,' +/-',y
      totint=totint+x
      toterr=toterr+y**2
      nentries=nentries+1

      read (20,'(a)',err=101,end=101) pippo

      goto 1

  101 continue

      toterr=sqrt(toterr)/nentries
      totint=totint/nentries

      print*, '    '
      print*,'nentries=',nentries
      print*, '    '
      print*, 'total cross section=', totint, ' +/-', toterr
      
      stop
      end
