***********************************************************************
*                       SUBROUTINE CCFCSYM                            *
*                                                                     *
*                                                                     *
*  Purpose:  It determines whether the Charge Conjugate (CC),         *
*   Family Conjugate (FC) and CC+FC of a given processes have to      *
*   be summed upon by fxn.f (with the appropriate PDF) or some of     * 
*   them are just the same process and must not be summed             *
*                                                                     *
*  The process is determined by the particle identities iproc(8)      *
*     iproc(1) and iproc(2) are incoming, the other six are outgoing  *
*                                                                     *
*      The CC, FC, CCFC processes are generated.                      *
*      Every one is compared  with the previous ones (including       *
*      iproc itself) and if different the corresponding flag is set   *
*      to 1                                                           *
*                                                                     *
*      iccsum =1 means  CC has to be summed                           *
*      ifcsum =1 means  FC has to be summed                           *
*      iccfcsum =1 means  CCFC has to be summed                       * 
*                                                                     *
***********************************************************************


      subroutine CCFCSYM(iproc,iccsum,ifcsum,iccfcsum)
      
      implicit real*8 (a-h,o-z)
      
      dimension iproc(8),icc(8), ifc(8), iccfc(8) 

      COMMON/phaones/ionesh

c   PDG  convention:
c   1=d, 2=u, 3=s, 4=c, 5=b, 6=t
c   11=e-, 12=v_e, 13=mu-, 14=v_mu, 15=tau-, 16=v_tau
c   all antiparticles have the same number but opposite sign
c   moreover: 21=gluon, 22=gamma, 23=Z0, 24=W+, 25=h


* charge conjugation array
*       Careful CC of non existing particles are set =0 !!!!
      INTEGER iccconj(-25:25)
      DATA    iccconj/
     &     0,24,0,0,0,0,0,0,0,16,15,14,13,12,11,0,0,0,0,6,5,4,3,2,1,0,
     &     -1,-2,-3,-4,-5,-6,0,0,0,0,-11,-12,-13,-14,-15,-16,
     &     0,0,0,0,21,22,23,-24,25 /

* family conjugation array
      INTEGER ifamconj(-25:25)
      DATA    ifamconj/
     &  -25,-24,-23,-22,-21,-20,-19,-18,-17,-16,
     &  -15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,
     &         -2,-1,-4,-3,0,3,4,1,2,
     &  5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 /


c initialization

      iccsum=0
      ifcsum=0
      iccfcsum=0


      do i=1,8
        icc(i)=iccconj(iproc(i))
        ifc(i)=ifamconj(iproc(i))
        iccfc(i)=iccconj(ifc(i))
      enddo

ctest
      if (ionesh.ne.1) then
        print*, 'iproc',iproc
        print*, 'icc  ',icc
        print*, 'ifc  ',ifc
        print*, 'iccfc',iccfc
      endif
ctestend      

      call comptot(icc,iproc,idiv)
      if (idiv.eq.1) iccsum=1

      call comptot(ifc,iproc,idiv)
      if (idiv.eq.1.and.iccsum.eq.1) then
        call comptot(ifc,icc,idiv)
      endif
      if (idiv.eq.1) ifcsum=1

      call comptot(iccfc,iproc,idiv)
      if (idiv.eq.1.and.iccsum.eq.1) then
        call comptot(iccfc,icc,idiv)
      endif
      if (idiv.eq.1.and.ifcsum.eq.1) then
        call comptot(iccfc,ifc,idiv)
      endif
      if (idiv.eq.1) iccfcsum=1

      return
      end


***********************************************************************
*                     SUBROUTINE COMPTOT                              *
*                                                                     *
*  Purpose:  It determines whether two group of particles 8 particles *
*   iproc1 and iproc2  represent the same process apart from ordering.*
*  The first two particles of each group are assumed to be incoming,  *
*   the other six are outgoing. Therefore the two processes iproc1    *
*   and iproc2 are the same if the content of the first two is the    *
*   same AND the content of the second six is the same                *
*                                                                     *
*     idiv=1 if they two processes are not identical                  *
*     idiv=0 if they are identical                                    *
*                                                                     *
***********************************************************************

      subroutine comptot(iproc1,iproc2,idiv)
      dimension iproc1(8),iproc2(8)
      
      idiv=0
      call comp(iproc1,iproc2,2,idiv1)
      if (idiv1.eq.1) then
        idiv=1
      else
        call comp(iproc1(3),iproc2(3),6,idiv)
      endif

      return
      end



***********************************************************************
*                         SUBROUTINE COMP                             *
*                                                                     *
*  Purpose:  It determines whether two group of particles idp and ip  *
*     of dimension n (actually max n = 8, but it may be changed       *
*     in the parameter) are the same apart from their order.          *
*                                                                     *
*     idiv=1 if they are not identical                               *
*     idiv=0 if they are identical                                   *
*                                                                     *
***********************************************************************


      subroutine comp(idp,ip,n,idiv)
      parameter (maxn=8)
      dimension idp(maxn),ip(maxn),idploc(maxn)

      idiv=0

      do j=1,n
        idploc(j)=idp(j)
      enddo

      do j=1,n
        i=j
c   one compares ip(j) with  idploc(i) (i=j,n) . If none is 
c    equal, i become greater than n  ad idiv is set to 1
c    if the i-th is equal, this is exchanged whith the one in position
c    j and one starts again 
        do while(idploc(i).ne.ip(j).and.i.le.n)
          i=i+1
        enddo
        if (i.gt.n) then
          idiv=1
          return
        else
          iaux=idploc(j)
          idploc(j)=idploc(i)
          idploc(i)=iaux
        endif
      enddo

      return
      end
