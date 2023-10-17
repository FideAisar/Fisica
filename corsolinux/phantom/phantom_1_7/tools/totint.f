* This program computes the total integral of a grid computation 
*  taking into account the factor two which has to be given
* to the processes in which the incoming particles are not identical.
* The input file can be prepared from a directory, in which are computed
*  all processes for the grids, with the command:
*  grep SIGMA */run.out > res  if run.out is the name of the files containing
*  the results of the integration of all * processes
*       (grep SIGMA */run.o* > res for Torino queues)
* In this way in the file res there will be ne name of all processes followed by
*  the corresponding cross section SIGMA

* Then you run the program .../totint.exe > result
*  In the results file in which you sent the output you will get
*  the name of all processes, their multiplication factor, the cross section
*   and at the end the total cross section with its extimated error  
 

      implicit real*8 (a-h,o-z)
      character * 80 pippo
      character * 80 sigma

c      print*, 'nome del file?'
c      read (5,'(a)') pippo
      pippo='res'
      npi=index(pippo,' ')-1

      open (unit=20,status='old',file=pippo(1:npi))     

      totint=0.d0
      toterr=0.d0
      nentries=0

    1 read (20,'(a)',err=101,end=101) pippo
      np1=index(pippo,'/')
      np2=index(pippo,'SIGMA')
      sigma=pippo(np2:)
      read (sigma(9:),*) x
      read (sigma(27:),*) y
      print *, pippo(1:np1)
c      print*,x
c      print*,y
      read (sigma(9:),*) x
      read (sigma(27:),*) y
      npi=index(pippo,'-')-1
c      print*,npi
c      print*, pippo(1:npi)
c      print*,x,y
      if (mod(npi,2).eq.0.and.
     &     (pippo(1:(npi)/2).eq.pippo(npi/2+1:npi))) then
        totint=totint+x
        toterr=toterr+y**2
        print*, '*1'
      else
        totint=totint+2.d0*x
        toterr=toterr+2.d0*y**2
        print*, '*2'
      endif

      print *,sigma

      nentries=nentries+1
c      print*, 'totint, toterr' ,totint, toterr

****
      goto 1
  101 continue

      toterr=sqrt(toterr)

      print*,'nentries=',nentries
      print*, '    '
      print*, 'total cross section=', totint, ' +/-', toterr

      stop
      end

