************************************************************************
*                        SUBROUTINE READINPUT                          *
*                                                                      *
*                                                                      *
*     Purpose:  reading input_file                                     *
*     Call to subroutines in wread.f                                   *
*                                                                      *
************************************************************************

      SUBROUTINE readinput

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT DOUBLE COMPLEX (c)

      INCLUDE 'common.h'
      INCLUDE 'common_cut.h'
      INCLUDE 'common_subproc.h'
      INCLUDE 'common_unitarization.h'

*     INTEGRATION COMMONS
c giuseppe 21/05/2007
      DIMENSION ncall_therm(2)
c end giuseppe 21/05/2007


      COMMON/iinteg_readinput/ncall_therm,itmx_therm,ncall,itmx
      COMMON/rinteg_readinput/acc_therm,acc
c giuseppe 07/02/2007
      COMMON/ifxn_readinput/i_PDFscale
*1_7
      COMMON/rfxn_readinput/pdfconst
*1_7end
c end giuseppe 07/02/2007
c s14
* heavh
*      COMMON/coupling_readinput/ghfactor
      COMMON/coupling_readinput/ghfactor,ghhfactor,rcosa,tgbeta
      COMMON/int_coupling_readinput/i_singlet,i_hh
*hevhend 
c s14end
c scale choice
      COMMON/rfxn_readinput/fixed_PDFscale
c end scale choice 
      COMMON/phaones/ionesh

      COMMON /pharand/ idum
*PDF's
      CHARACTER*200 PDFname
      COMMON/phpdfli/ PDFname
*zw
* *wwpol
*       common/polww/i_ww,i_pol,i_polwp,i_polwm,i_polwel,i_polwmu
* c DPA (giovanni, 5/12/2016) 
* *res
* c      common/dpa/i_dpa
*       common/osp/i_osp
* *resend
* c DPA end (giovanni, 5/12/2016) 
* *wwpolend      
*zwend
      CALL iread('idum',idum,1) ! idum=initialization random number seed
                                !  must be a large negative number

      CALL cread('PDFname',PDFname)

c giuseppe 06/02/2007
      CALL iread('i_PDFscale',i_PDFscale,1) ! selects the way of 
                                            ! calculating the PDF scale:
                                            ! =1 for all processes, 
                                            !     based on pT's of ALL
                                            !     OUTGOING PARTICLES
                                            ! =2 process by process,
                                            !     based on pT of the
                                            !     (RECONSTRUCTED) TOP
                                            !     if possible, otherwise
                                            !     as done in option 1
c end giuseppe 06/02/2007
c scale choice
                                            ! =3 Q a fixed numerical scale 
                                            !      given in r.in
                                            ! =4 Q= m4l/sqrt(2)(invariant mass 
                                            !         of the 4 leptons)/sqrt(2) 
                                            !  (valid only for prosesses with 
                                            !       four outgoing leptons)
c giovanni 2017/10/24 begin
                                            ! =5  Q= sqrt(ptj1*ptj2)
                                            !     square-root of the product of
                                            !     pt's of the 2 jets with largest pt
                                            !     (for all processes with
                                            !     at least 2 final statejets)
c giovanni 2017/10/24 end


      if (i_PDFscale.eq.3) then
        CALL rread('fixed_PDFscale',fixed_PDFscale,1)
      endif
c end scale choice 
*1_7
      CALL rread('pdfconst',pdfconst,1) ! Fix the numerical constant by which 
                                        ! the PDFscale is multiplied
*1_7end

c giuseppe ILC 25/06/2007
      CALL iread('i_coll',i_coll,1) ! determines the type of accelerator
                                    !  i_coll=1 => p-p  
                                    !  i_coll=2 => p-pbar
                                    !  i_coll=3 => e+e- 
      IF(i_coll.EQ.3)THEN
        CALL iread('i_isr',i_isr,1) ! yes/no initial state radiation
                                    ! (ISR) for e+e- collider only
        CALL iread('i_beamstrahlung',i_beamstrahlung,1) 
                                    ! yes/no beamstrahlung 
                                    ! for e+e- collider only
      ENDIF                     ! IF(i_coll.EQ.3)THEN
c end giuseppe ILC 25/06/2007

      CALL iread('perturbativeorder',i_pertorder,1 ) 
              !i_pertorder = 1 alpha_em^6 with dedicated amp
              !           = 2 alpha_s^2alpha_em^4
              !           = 3 alpha_em^6 + alpha_s^2alpha_em^4
              !           = 0 alpha_em^6 with amp8fqcd (for test only)
 

*sandro2/3/07end

*diogo29/3/09
      CALL iread('i_massive',i_massive,1 ) 
              !i_massive = 0 use massless amp unless there is at least a b quark
              !          = 1 always use massive amplitudes (massive Z-lines)
*diogoend

* oldread

      i_unitarize=0

c$$$*diogo24/4/10 Unitarization
c$$$
c$$$C     iunittype= 0 no unitarization
c$$$C            1 kmatrix
c$$$C            2 pade
c$$$C            3 nd
c$$$C            4 largen     
c$$$      CALL iread('i_unitarize',i_unitarize,1)
c$$$            !i_unitarize = 0   no changes in amplitudes
c$$$	    !            = 1   use unitarization routines
c$$$      if (i_unitarize.eq.1) then
c$$$        CALL iread('iunittype',iunittype,1)
c$$$C         iunittype= 0 no unitarization
c$$$C                    1 kmatrix
c$$$C                    2 pade
c$$$C                    3 nd
c$$$C                    4 largen
c$$$      
c$$$        CALL iread('inlo',inlo,1) !yes/no NLO 
c$$$        CALL iread('ireson',ireson,1) ! yes/no resonances
c$$$        CALL iread('isigma',isigma,1) 
c$$$C         yes/no sigma isosinglet scalar resonance IJ=00.
c$$$        CALL iread('irho',irho,1) 
c$$$C         yes/no rho isotriplet vector resonance IJ=11.      
c$$$        CALL iread('iphi',iphi,1) 
c$$$C         yes/no phi isosinglet scalar resonance IJ=20.
c$$$        CALL iread('iff',iff,1) 
c$$$C         yes/no f isosinglet tensor resonance IJ=02.      
c$$$        CALL iread('ia',ia,1) 
c$$$C         yes/no a isoquintet tensor resonance IJ=22.            
c$$$
c$$$C rg_R  resonances couplings
c$$$C rm_R  resonances masses      
c$$$C gam_R resonances widths, must be specified in case you don't 
c$$$C       use a unitarization scheme.
c$$$            
c$$$        CALL rread('rg_sigma',rg_sigma,1)
c$$$        CALL rread('rm_sigma',rm_sigma,1)
c$$$        CALL rread('gam_sigma',gam_sigma,1)
c$$$      
c$$$        CALL rread('rg_rho',rg_rho,1)
c$$$        CALL rread('rm_rho',rm_rho,1)
c$$$        CALL rread('gam_rho',gam_rho,1)
c$$$      
c$$$        CALL rread('rg_phi',rg_phi,1)
c$$$        CALL rread('rm_phi',rm_phi,1)
c$$$        CALL rread('gam_phi',gam_phi,1)
c$$$      
c$$$        CALL rread('rg_ff',rg_ff,1)
c$$$        CALL rread('rm_ff',rm_ff,1)
c$$$        CALL rread('gam_ff',gam_ff,1)      
c$$$      
c$$$        CALL rread('rg_a',rg_a,1)
c$$$        CALL rread('rm_a',rm_a,1)
c$$$        CALL rread('gam_a',gam_a,1)      
c$$$
c$$$C NLO parameter, alphas are higher order operators parameters
c$$$C rmu is the renormalization scale.
c$$$        CALL rread('rmu',rmu,1)
c$$$        CALL rread('alpha4',alpha4,1)
c$$$        CALL rread('alpha5',alpha5,1)
c$$$
c$$$C N/D scheme: mass parameter       
c$$$        CALL rread('rmnd',rmnd,1)      
c$$$
c$$$      endif
c$$$
c$$$*end diogo

*oldreadend


      CALL iread('ionesh',ionesh,1) 
                     ! 0= normal run of one process   
                     ! 1= one shot generation of all indicated processes
*zw
      if (ionesh.eq.0) then
        CALL iread('iproc',iproc,8) ! process
      endif
*zwend 

     
      CALL rread('ecoll',ecoll,1) ! collider energy

      CALL rread('rmh',rmh,1)   ! Higgs mass (GeV) 
	                         ! <0 means no Higgs  
 
c s14
      CALL rread('ghfactor',ghfactor,1)   ! factor which multiplies SM higgs
                                          ! couplings
                                          ! for SM use    ghfactor    1.d0
                  ! when  i_singlet=1 (see below) ghfactor not considered !!!

      if (ghfactor.ne.1.d0.and.i_pertorder.ne.1) then
        print*, 'ERROR: '
        print*, ' higgs couplings different from SM '
        print*, ' (ghfactor different from 1.d0) '
        print*, ' implemented only for  alpha_em^6  (i_pertorder = 1)'
        stop
      endif
      
c s14end

      CALL rread('gamh',gamh,1)   ! Higgs width (GeV)
*                                  ! <0 means computed by phantom
*	                           ! in this last case SM gamh is multiplied
*                                  ! by ghfactor**2 
*                                  ! or by rcosa**2 if i_singlet=1 (see below)

*zw

      CALL iread('i_ww',i_ww,1) ! i_ww= 0 full computation
*                               ! i_ww= 1 only 1 resonant  w diagrams 
*                               ! i_ww= 2 only 2 resonant  w diagrams 


       if (i_ww.ge.1) then
         CALL iread('idw',idw,4)!(four numbers must be given, 
*                                ! but only the first two are considered
*                                !  if i_ww=1, all of them if i_ww=2)
*                                ! the first two correspond to the decay of the
*                                ! first w , the second two eventually 
*                                ! to the decay of the second w
*                                ! The first number of any couple must 
*                                ! correspond to the particle, the second to 
*                                ! the antiparticle (negative) of the decay.  

         jend=0
         ireap=0
         if (ionesh.eq.0) then
           if (i_ww.eq.1) then
             jend=2
             if (idw(1).le.0.or.idw(2).ge.0) then
               print*, 'ERROR:' 
               print*, 'the order particle antiparticle not respected'
               stop
             endif
           endif
           if (i_ww.eq.2) then
             jend=4
             if (idw(1).le.0.or.idw(2).ge.0.or.
     &            idw(3).le.0.or.idw(4).ge.0) then
               print*, 'ERROR:' 
               print*, 'the order particle antiparticle not respected'
               stop
             endif
           endif
           do j=1,jend
             do k=j+1,jend
               if (idw(j).eq.idw(k)) ireap=1
             enddo
               
             neq=0
             do i=3,8
               if (iproc(i).eq.idw(j)) neq=neq+1
             enddo
             if (neq.ne.1.or.ireap.eq.1) then
               print*, 'ERROR'
               print*, 'some resonant particle is not in the outgoing'
               print*, 'ones or it appears more than once !'
               stop
             endif
           enddo
         endif
     


         CALL iread('ipolw',ipolw,2) ! the first index refers to the 
                                      ! polarization of the first w, 
                                      ! the second of the second if i_ww=2
                                      ! the indexes can be: 
                                      ! 0 no polarization, 1 longitudinal,
                                      ! 2 left, 3 right ,4 transverse 
       endif
       

      CALL iread('i_zz',i_zz,1) ! i_zz= 0 full computation
*                               ! i_zz= 1 only 1 resonant  z diagrams 
*                               ! i_zz= 2 only 2 resonant  z diagrams 

       if (i_zz.ge.1) then
         CALL iread('idz',idz,4)!(four numbers must be given, 
*                                ! but only the first two are considered
*                                !  if i_zz=1, all of them if i_zz=2)
*                                ! the first two correspond to the decay of the
*                                ! first z , the second two eventually 
*                                ! to the decay of the second z
*                                ! The first number of any couple must 
*                                ! correspond to the particle, the second to 
*                                ! the antiparticle (negative) of the decay  


         jend=0
         ireap=0
         if (ionesh.eq.0) then
           if (i_zz.eq.1) then
             jend=2
             if (idz(1).le.0.or.idz(2).ge.0) then
               print*, 'ERROR:' 
               print*, 'the order particle antiparticle not respected'
               stop
             endif
           endif
           if (i_zz.eq.2) then
             jend=4
             if (idz(1).le.0.or.idz(2).ge.0.or.
     &            idz(3).le.0.or.idz(4).ge.0) then
               print*, 'ERROR:' 
               print*, 'the order particle antiparticle not respected'
               stop
             endif
           endif
           do j=1,jend
             do k=j+1,jend
               if (idz(j).eq.idz(k)) ireap=1
             enddo
               
             neq=0
             do i=3,8
               if (iproc(i).eq.idz(j)) neq=neq+1
             enddo
             if (ionesh.eq.0) then
               if (neq.ne.1.or.ireap.eq.1) then
                 print*, 'ERROR'
                 print*, 'some resonant particle is not in the outgoing'
                 print*, 'ones or it appears more than once !'
                 stop
               endif
             endif
           enddo
         endif
     


         CALL iread('ipolz',ipolz,2) ! the first index refers to the 
                                      ! polarization of the first z, 
                                      ! the second of the second if i_zz=2
                                      ! the indexes can be: 
                                      ! 0 no polarization, 1 longitudinal,
                                      ! 2 left, 3 right ,4 transverse 
       endif
       

      if ((i_ww.gt.0.or.i_zz.gt.0).and.i_pertorder.gt.1) then
        print *, '                                     ERROR : ' 
        print*, ' RESONANT COMPUTATIONS, ON SHELL PROJECTIONS AND'
        print*, ' POLARIZATIONS CAN BE USED ONLY FOR  i_pertorder = 1'
        stop
      endif

* comment 1 res w  and 1 res z diagrams are obtained fixing i_ww=1 and i_zz=1

*4cmpol
      if (i_ww.ge.1.or.i_zz.ge.1) then
        if (ipolw(1).gt.0.or.ipolw(2).gt.0.
     &       or.ipolz(1).gt.0.or.ipolz(2).gt.0) then
          CALL iread('i_4cmpol',i_4cmpol,1) 
                               ! i_4cmpol = 0  polarizations defined in the lab
                               ! i_4cmpol = 1  polarizations defined in cm of
                               !    four particles to be indicted below 

          if (i_4cmpol.gt.0) then
            CALL iread('id4cmpol',id4cmpol,4)  
                                       ! identity of the particles which form
                                       ! the cm in which the polarizations 
*check                                       ! are defined
            ireap=0
            do j=1,4
              do k=j+1,4
                if (id4cmpol(j).eq.id4cmpol(k)) ireap=1
              enddo
               
              neq=0
              do i=3,8
                if (iproc(i).eq.id4cmpol(j)) neq=neq+1
              enddo
              if (neq.ne.1.or.ireap.eq.1) then
                if (ionesh.eq.0) then
                  print*, 'ERROR'
                print*,'some  particle for 4 cm is not in the outgoing'
                  print*, 'ones or it appears more than once !'
                  stop
                endif
              endif
*check end
            enddo
          endif
        endif
      endif

*4cmpolend      

      if (i_ww.ge.1.or.i_zz.ge.1 ) then
        CALL iread('i_osp',i_osp,1) ! i_osp = 0  no kinematics change
                                    ! i_osp = 1  on shell projection scheme for
                                    !             1 boson decaying
                                    ! i_osp = 2  on shell projection scheme for
                                    !             2 bosons decaying to letptons

         if (i_osp.gt.0) then
           CALL iread('idosp',idosp,4) ! identity of the particles which must
                                       ! be projected. Only the first couple
                                       ! counts if i_osp.eq.1.
                                       ! For every couple the first is the 
                                       ! particle, the second the antiparticle
         endif
       endif
*zwend

* sig
      CALL iread('i_signal',i_signal,1) ! i_signal= 0 full computation 
                                        ! i_signal> 0 Higgs signal only 
                                        !        (only for i_pertorder = 1 
                                        !                       alpha_em^6
                                        !              and i_unitarize = 0)
                                        ! i_signal= 1 s channel contributions 
                                        !        to boson boson scattering 
                                        !   (boson boson-> Higgs -> boson boson)
                                        ! i_signal= 2 all contributions 
                                        !        (s+t+u channels)
                                        !        to boson boson scattering 
                                        ! i_signal= 3 all contributions 
                                        !        (s+t+u channels)
                                        !        to boson boson scattering 
                                        !        and Higgstrahlung with 
                                        !           H -> boson boson

      if (i_signal.gt.0.and.i_pertorder.ne.1) then
         print*, 'ERROR: for i_signal.gt.0 i_pertorder must be 1'
         stop
      endif
      if (i_signal.gt.0.and.i_unitarize.ne.0) then
        print*,  'ERROR: for i_signal.gt.0 i_unitarize must be 0)'
        stop
      endif
 
* sigend

* heavh  

* SINGLET MODEL OPTION

* singlet model implementation (see e.g. Pruna Roberts arXiv:1303.1150) 

      CALL iread('i_singlet',i_singlet,1) ! yes/no singlet implementation 

      if (i_singlet.eq.1.and.i_pertorder.ne.1) then
        print*, 'WARNING: '
        print*, ' higgs couplings for Singlet Model '
        print*, ' implemented only for  alpha_em^6 '
      endif

* SINGLET MODEL PARAMETERS
      if (i_singlet.eq.1) then

        CALL rread('rmhh',rmhh,1) ! mass of the heavier higgs

        CALL rread('rcosa',rcosa,1) ! parameter cos alfa of arXiv:1303.1150

        CALL rread('tgbeta',tgbeta,1) ! parameter  tg beta of arXiv:1303.1150

        CALL rread('gamhh',gamhh,1) ! heavier Higgs width (GeV)
                                    ! <0 means computed by phantom
                                    ! in this last case SM gamh is multiplied
                                    ! by (1-rcosa**2) + decay of heavy higgs 
                                    !     to 2 light ones: hh-> h+h 
      endif


* HEAVY HIGGS NOT IN THE SINGLET CONTEST (for test only)

      CALL iread('i_hh',i_hh,1) ! yes/no heavy higgs ( not singlet) 


      if (i_hh.eq. 1) then

        if (i_singlet.eq.1) then
          print*,'ERROR'
          print*,'i_hh and i_singlet cannot be both = 1'
          stop
        endif

*from now on parameters regarding a second heavier higgs scalar
****    hh stays for heavy higgs. 


        CALL rread('rmhh_ns',rmhh,1) ! heavier higgs mass (GeV) (not singlet)

        if (i_pertorder.ne.1) then
          print*, 'ERROR: '
          print*, ' higgs couplings different from SM '
          print*, ' i_hh=1 '
          print*, ' implemented only for  alpha_em^6  (i_pertorder = 1)'
          stop
        endif

        CALL rread('ghhfactor',ghhfactor,1) ! factor for second higgs, 
                                           ! which multiplies SM higgs
                                          ! couplings 
                 ! if ghhfactor is negative, ghhfactor=sqrt(1-ghfactor**2)

        IF (ghhfactor.lt.0.d0) ghhfactor=sqrt(1.d0-ghfactor**2)
     
        if (ghhfactor.le.0.d0) then
          print*, ' WARNING '
          print*, ' a hevy higgs has been requested  '
          print*, ' but the couplings are zero'
*          stop
        endif

        CALL rread('gamhh_ns',gamhh,1) ! heavier Higgs width (GeV)
                                    ! <0 means computed by phantom
                                    ! in this last case SM gamh is multiplied
                                    ! by ghhfactor**2 
      endif
      
* set rmhh <0.d0 ad a flag for no heavy higgs (no i_singlet nor i_hh =1)
      if (i_singlet.ne.1.and.i_hh.ne.1) rmhh=-1.d0

* heavend


      CALL iread('i_ccfam',i_ccfam,1)           ! family+CC conjugate

*zw
      if (i_ccfam.eq.1) then
        if (i_ww.gt.0) then
          print*, 'ERROR: i_ccfam cannot be 1 if i_ww>0'
          stop
        endif
        if(i_zz.eq.1) then
          if (idz(1).lt.11.or.idz(1).gt.16) then   
            print*, 'ERROR: i_ccfam cannot be 1' 
            print*, '       if resonant Z particles are not leptons'
            stop
          endif
        elseif(i_zz.eq.2) then
          if (idz(1).lt.11.or.idz(1).gt.16.or.idz(3).lt.11.
     &         or.idz(3).gt.16) then
            print*, 'ERROR: i_ccfam cannot be 1' 
            print*, '       if resonant Z particles are not leptons'
            stop
          endif
        endif
      endif
*zwend
       

      if (ionesh.eq.0) then

*     READ INPUT FOR THE INTEGRATION
        CALL rread('acc_therm',acc_therm,1) ! thermalization accuracy
c giuseppe 21/05/2007
c        CALL iread('ncall_therm',ncall_therm,1)
c                                ! thermalization calls per iteration
        CALL iread('ncall_therm',ncall_therm,2)
                                ! thermalization calls per iteration.
                                ! The first component refers to the 
                                ! number of calls for the first 3 
                                ! iterations, the second one to the 
                                ! calls for the remaining iterations.
c end giuseppe 21/05/2007
        CALL iread('itmx_therm',itmx_therm,1) !thermalization iterations
        CALL rread('acc',acc,1) ! integration accuracy
        CALL iread('ncall',ncall,1) ! integration calls per iteration
        CALL iread('itmx',itmx,1) ! integration iterations 

*     READ INPUT FOR FLAT EVENT GENERATION
        CALL iread('iflat',iflat,1) ! yes/no flat event generation
        IF(iflat.EQ.1)THEN
c          CALL rread('scalemax0',scalemax0,1) 
c                                !scale factor for the maximum
c          scalemax=scalemax0
          scalemax=1.1

c          CALL iread('istorvegas',istorvegas,1)
                                ! yes/no VEGAS data stored
          istorvegas=1

c          CALL iread('iwrite_event',iwrite_event,1)
                  ! yes/no momenta of flat events written in .dat files

          iwrite_event=0

        ENDIF   !iflat

      elseif (ionesh.eq.1) then

*     scale factor for generation

        CALL rread('scalemax',scalemax,1) 
                            !scale factor for the maximum

        CALL iread ('nunwevts',nunwevts,1)	 
               !  number of unweighted events to be produced

        CALL iread('iwrite_event',iwrite_event,1)
                  ! yes/no momenta of flat events written in .dat files

c giuseppe 29/09/2007
        CALL iread('iwrite_mothers',iwrite_mothers,1)
                                ! yes/no information about intermediate 
                                ! particles (mothers) in .dat files
c end giuseppe 29/09/2007

        CALL iread('ihadronize',ihadronize,1) 
                                ! yes/no call to hadronization

        CALL iread('i_exchincoming',i_exchincoming,1)  

c        CALL iread('i_emutau',i_emutau,1)  
        i_emutau=0
      endif

*     READ INPUT FOR CUTS
      CALL iread('i_e_min_lep',i_e_min_lep,1) ! lepton energy lower cuts 
                                !     (GeV)
      IF(i_e_min_lep.EQ.1) CALL rread('e_min_lep',e_min_lep,1) 
      CALL iread('i_pt_min_lep',i_pt_min_lep,1) 
                                ! lepton pt lower cuts (GeV)
      IF(i_pt_min_lep.EQ.1) CALL rread('pt_min_lep',pt_min_lep,1)

c giuseppe 06/02/2007
      CALL iread('i_eta_max_onelep',i_eta_max_onelep,1)
                                ! maximum rapidity (absolute value) for
                                ! AT LEAST one lepton, i.e. at least one
                                ! final state lepton is required to be
                                ! central
c end giuseppe 06/02/2007

      CALL iread('i_eta_max_lep',i_eta_max_lep,1) 
                                ! lepton rapidity upper cuts (absolute
                                ! value), i.e. ALL final state leptons
                                ! are required to be central
c giuseppe 01/03/2007
      IF(i_eta_max_onelep.EQ.1)THEN
        CALL rread('eta_max_onelep',eta_max_onelep,1)
      ENDIF
      IF(i_eta_max_lep.EQ.1)THEN
        CALL rread('eta_max_lep',eta_max_lep,1)      
      ENDIF
c end giuseppe 01/03/2007

      CALL iread('i_ptmiss_min',i_ptmiss_min,1) 
                                ! missing pt lower cuts (GeV)
      IF(i_ptmiss_min.EQ.1) CALL rread('ptmiss_min',ptmiss_min,1)
      CALL iread('i_e_min_j',i_e_min_j,1) ! jet energy lower cuts (GeV)
      IF(i_e_min_j.EQ.1) CALL rread('e_min_j',e_min_j,1)
      CALL iread('i_pt_min_j',i_pt_min_j,1) ! jet pt lower cuts (GeV)
      IF(i_pt_min_j.EQ.1) CALL rread('pt_min_j',pt_min_j,1)
      CALL iread('i_eta_max_j',i_eta_max_j,1) 
                                ! jet rapidity upper cuts  (absolute value)
      IF(i_eta_max_j.EQ.1) CALL rread('eta_max_j',eta_max_j,1)

      CALL iread('i_eta_jf_jb_jc',i_eta_jf_jb_jc,1)
                                ! rapidity of forward, backward and central jets
      IF(i_eta_jf_jb_jc.EQ.1) THEN
        CALL rread('eta_def_jf_min',eta_def_jf_min,1) 
                                ! min rapidity for a jet to be called forward
        CALL rread('eta_def_jb_max',eta_def_jb_max,1) 
                                ! max rapidity for a jet to be called backward
        CALL rread('eta_def_jc_max',eta_def_jc_max,1) 
                                ! max rapidity for a jet to be called central (absolute value)
      ENDIF ! (i_eta_jf_jb_jc.EQ.1)

      CALL iread('i_pt_min_jcjc',i_pt_min_jcjc,1) ! pt lower cuts 
                                ! on two centraljets (GeV)
      IF(i_pt_min_jcjc.EQ.1) 
     &     CALL rread('pt_min_jcjc',pt_min_jcjc,1)

      CALL iread('i_rm_min_jj',i_rm_min_jj,1) 
                                ! minimum invariant mass  between jets (GeV)
      IF(i_rm_min_jj.EQ.1) CALL rread('rm_min_jj',rm_min_jj,1)

      CALL iread('i_rm_min_ll',i_rm_min_ll,1) 
                         ! minimum invariant mass  between charged leptons (GeV)
      IF(i_rm_min_ll.EQ.1) CALL rread('rm_min_ll',rm_min_ll,1)

      CALL iread('i_rm_min_jlep',i_rm_min_jlep,1) 
                                ! minimum invariant mass between jets and lepton (GeV)
      IF(i_rm_min_jlep.EQ.1) CALL rread('rm_min_jlep',rm_min_jlep,1)
      CALL iread('i_rm_min_jcjc',i_rm_min_jcjc,1)
                                ! minimum invariant mass between central jets (GeV)
      IF(i_rm_min_jcjc.EQ.1) 
     &     CALL rread('rm_min_jcjc',rm_min_jcjc,1) 
      CALL iread('i_rm_max_jcjc',i_rm_max_jcjc,1) 
                                ! maximum invariant mass between central jets (GeV)
      IF(i_rm_max_jcjc.EQ.1) 
     &     CALL rread('rm_max_jcjc',rm_max_jcjc,1) 
      CALL iread('i_rm_min_jfjb',i_rm_min_jfjb,1) 
                                ! minimum invariant mass between forward and backward jet
      IF(i_rm_min_jfjb.EQ.1) 
     &     CALL rread('rm_min_jfjb',rm_min_jfjb,1)

* six
      CALL iread('i_rm_min_4l',i_rm_min_4l,1) 
                     ! minimum invariant mass of 4 leptons for processes with 4l
      IF(i_rm_min_4l.EQ.1) 
     &     CALL rread('rm_min_4l',rm_min_4l,1)

      if (i_osp.eq.2.and.rm_min_4l.lt.2*rmw) then
        print*, '      i_osp=2 requires rm_min_4l > 2 mW'
        stop
      endif

      CALL iread('i_rm_min_2l2cq',i_rm_min_2l2cq,1) 
                    ! minimum invariant mass of 2 leptons and 
                    ! 2 central quark for processes with 2l
      IF(i_rm_min_2l2cq.EQ.1) 
     &     CALL rread('rm_min_2l2cq',rm_min_2l2cq,1)

* YR
      CALL iread('i_rm_max_4l',i_rm_max_4l,1) 
                     ! maximum invariant mass of 4 leptons for processes with 4l
      IF(i_rm_max_4l.EQ.1) 
     &     CALL rread('rm_max_4l',rm_max_4l,1)
      CALL iread('i_rm_max_2l2cq',i_rm_max_2l2cq,1) 
                    ! maximum invariant mass of 2 leptons and 2 central quark for processes with 2l
      IF(i_rm_max_2l2cq.EQ.1) 
     &     CALL rread('rm_max_2l2cq',rm_max_2l2cq,1)
* YRend

* sixend

* cuttop

      CALL iread('i_deltacuttop',i_deltacuttop,1) 
                     ! yes/no cut on the invariant mass of all triplets of
                     ! particles who could form a top 
                     ! to avoid top contributions
                     ! The corresponding deltacuttop fixes the 
                     ! interval excluded : topmas +- deltatacuttop
 
      IF(i_deltacuttop.EQ.1) 
     &     CALL rread('deltacuttop',deltacuttop,1)

* cuttopend

*ww
* rcutwlept

      CALL iread('i_cutwlept',i_cutwlept,1) 
*                     ! yes/no cut on the invariant mass of all couples
*                     ! of letpons who could form a W+ or W- 
*                     ! The corresponding rcutwlept fixes the 
*                     ! interval remaining : Wmass +- rcutwlept
 
      IF(i_cutwlept.EQ.1) 
     &     CALL rread('rcutwlept',rcutwlept,1)

* rcutptwelectron

      CALL iread('i_cutptwelectron',i_cutptwelectron,1) 
*                     ! yes/no cut on the pT of the electron-antineutrino
*                     ! electron pair
 
      IF(i_cutptwelectron.EQ.1) 
     &     CALL rread('rcutptwelectron',rcutptwelectron,1)

*wwend

*zw
      CALL iread('i_cutzlept',i_cutzlept,1) 
*                     ! yes/no cut on the invariant mass of all couples
*                     ! of letpons who could form a Z 
*                     ! The corresponding rcutzlept fixes the 
*                     ! interval remaining : Zmass +- rcutwlept
 
      IF(i_cutzlept.EQ.1) 
     &     CALL rread('rcutzlept',rcutzlept,1)
*zwend

      CALL iread('i_eta_min_jfjb',i_eta_min_jfjb,1) 
                 ! minimum rapidity difference between forward and backward jet
      IF(i_eta_min_jfjb.EQ.1) 
     &     CALL rread('eta_min_jfjb',eta_min_jfjb,1)
      CALL iread('i_d_ar_jj',i_d_ar_jj,1) 
                                ! minimum delta_R separation between jets
      IF(i_d_ar_jj.EQ.1) CALL rread('d_ar_jj',d_ar_jj,1)
      CALL iread('i_d_ar_jlep',i_d_ar_jlep,1)
                                ! minimum delta_R separation between jets and
				! lepton
      IF(i_d_ar_jlep.EQ.1) CALL rread('d_ar_jlep',d_ar_jlep,1)
***aggiunta 1.2.1 su richiesta Pietro
      CALL iread('i_d_ar_leplep',i_d_ar_leplep,1)
                      ! minimum delta_R separation between two charged leptons
				! 
      IF(i_d_ar_leplep.EQ.1) CALL rread('d_ar_leplep',d_ar_leplep,1)
***
      CALL iread('i_thetamin_jj',i_thetamin_jj,1) 
                                ! minimum angle separation between jets (cosine)
      IF(i_thetamin_jj.EQ.1) 
     &     CALL rread('thetamin_jj',thetamin_jj,1)
      CALL iread('i_thetamin_jlep',i_thetamin_jlep,1) 
                                ! minimum angle separation between jets and lepton (cosine)
      IF(i_thetamin_jlep.EQ.1) 
     &     CALL rread('thetamin_jlep',thetamin_jlep,1)
c giuseppe 01/03/2007
      CALL iread('i_thetamin_leplep',i_thetamin_leplep,1) 
                                ! minimum angle separation between charged leptons (cosine)
      IF(i_thetamin_leplep.EQ.1)THEN
        CALL rread('thetamin_leplep',thetamin_leplep,1)
      ENDIF
c end giuseppe 01/03/2007
      CALL iread('i_usercuts',i_usercuts,1) ! yes/no (1/0) additional
                                ! user-defined cuts

      if (ionesh.eq.1) then

*oldread
       iextracuts=0
c$$$
c$$$*     eventual more restrictive cuts
c$$$
c$$$
c$$$        CALL iread('iextracuts',iextracuts,1)
c$$$
c$$$        if (iextracuts. eq.1) then
c$$$
c$$$          CALL iread('i_e_min_lepos',i_e_min_lepos,1) 
c$$$                                ! lepton energy lower cuts (GeV)
c$$$          IF(i_e_min_lepos.EQ.1) CALL rread('e_min_lepos',e_min_lepos,1) 
c$$$          CALL iread('i_pt_min_lepos',i_pt_min_lepos,1) 
c$$$                                ! lepton pt lower cuts (GeV)
c$$$          IF(i_pt_min_lepos.EQ.1) 
c$$$     &         CALL rread('pt_min_lepos',pt_min_lepos,1)
c$$$
c$$$c giuseppe 06/02/2007
c$$$          CALL iread('i_eta_max_onelepos',i_eta_max_onelepos,1) 
c$$$                                ! maximum rapidity (absolute value) for
c$$$                                ! AT LEAST one lepton, i.e. at least one
c$$$                                ! final state lepton is required to be 
c$$$                                ! central
c$$$c end giuseppe 06/02/2007
c$$$
c$$$          CALL iread('i_eta_max_lepos',i_eta_max_lepos,1) 
c$$$                                ! lepton rapidity upper cuts (absolute 
c$$$                                ! value), i.e. ALL final state leptons 
c$$$                                ! are required to be central
c$$$c giuseppe 01/03/2007
c$$$          IF(i_eta_max_onelepos.EQ.1)THEN
c$$$            CALL rread('eta_max_onelepos',eta_max_onelepos,1)
c$$$          ENDIF
c$$$          IF(i_eta_max_lepos.EQ.1)THEN
c$$$            CALL rread('eta_max_lepos',eta_max_lepos,1)
c$$$          ENDIF
c$$$c end giuseppe 01/03/2007
c$$$
c$$$          CALL iread('i_ptmiss_minos',i_ptmiss_minos,1) 
c$$$                                ! missing pt lower cuts (GeV)
c$$$          IF(i_ptmiss_minos.EQ.1) 
c$$$     &         CALL rread('ptmiss_minos',ptmiss_minos,1)
c$$$          CALL iread('i_e_min_jos',i_e_min_jos,1) 
c$$$                                ! jet energy lower cuts (GeV)
c$$$          IF(i_e_min_jos.EQ.1) CALL rread('e_min_jos',e_min_jos,1)
c$$$          CALL iread('i_pt_min_jos',i_pt_min_jos,1) 
c$$$                                !jet pt lower cuts (GeV)
c$$$          IF(i_pt_min_jos.EQ.1) CALL rread('pt_min_jos',pt_min_jos,1)
c$$$          CALL iread('i_eta_max_jos',i_eta_max_jos,1) 
c$$$                             ! jet rapidity upper cuts  (absolute value)
c$$$          IF(i_eta_max_jos.EQ.1) CALL rread('eta_max_jos',eta_max_jos,1)
c$$$
c$$$
c$$$          CALL iread('i_eta_jf_jb_jcos',i_eta_jf_jb_jcos,1)
c$$$                             ! rapidity of forward, backward and central jets
c$$$          IF(i_eta_jf_jb_jcos.EQ.1) THEN
c$$$            CALL rread('eta_def_jf_minos',eta_def_jf_minos,1) 
c$$$                          ! min rapidity for a jet to be called forward
c$$$            CALL rread('eta_def_jb_maxos',eta_def_jb_maxos,1) 
c$$$                          ! max rapidity for a jet to be called backward
c$$$            CALL rread('eta_def_jc_maxos',eta_def_jc_maxos,1)
c$$$                ! max rapidity for a jet to be called central (absolute value)
c$$$          ENDIF   !  (i_eta_jf_jb_jcos.EQ.1)
c$$$
c$$$          CALL iread('i_pt_min_jcjcos',i_pt_min_jcjcos,1) 
c$$$                         ! pt lower cuts  on two centraljets (GeV)
c$$$          IF(i_pt_min_jcjcos.EQ.1) 
c$$$     &         CALL rread('pt_min_jcjcos',pt_min_jcjcos,1)
c$$$
c$$$          CALL iread('i_rm_min_jjos',i_rm_min_jjos,1) 
c$$$                           ! minimum invariant mass  between jets (GeV)
c$$$          IF(i_rm_min_jjos.EQ.1) CALL rread('rm_min_jjos',rm_min_jjos,1)
c$$$
c$$$          CALL iread('i_rm_min_llos',i_rm_min_llos,1) 
c$$$                   ! minimum invariant mass  between charged leptons (GeV)
c$$$          IF(i_rm_min_llos.EQ.1) CALL rread('rm_min_llos',rm_min_llos,1)
c$$$
c$$$          CALL iread('i_rm_min_jlepos',i_rm_min_jlepos,1) 
c$$$                      ! minimum invariant mass between jets and lepton (GeV)
c$$$          IF(i_rm_min_jlepos.EQ.1) 
c$$$     &         CALL rread('rm_min_jlepos',rm_min_jlepos,1)
c$$$          CALL iread('i_rm_min_jcjcos',i_rm_min_jcjcos,1)
c$$$                     ! minimum invariant mass between central jets (GeV)
c$$$          IF(i_rm_min_jcjcos.EQ.1) 
c$$$     &         CALL rread('rm_min_jcjcos',rm_min_jcjcos,1) 
c$$$          CALL iread('i_rm_max_jcjcos',i_rm_max_jcjcos,1) 
c$$$                     ! maximum invariant mass between central jets (GeV)
c$$$          IF(i_rm_max_jcjcos.EQ.1) 
c$$$     &         CALL rread('rm_max_jcjcos',rm_max_jcjcos,1) 
c$$$          CALL iread('i_rm_min_jfjbos',i_rm_min_jfjbos,1) 
c$$$               ! minimum invariant mass between forward and backward jet
c$$$          IF(i_rm_min_jfjbos.EQ.1) 
c$$$     &         CALL rread('rm_min_jfjbos',rm_min_jfjbos,1) 
c$$$
c$$$* six
c$$$      CALL iread('i_rm_min_4los',i_rm_min_4los,1) 
c$$$                     ! minimum invariant mass of 4 leptons for processes with 4l
c$$$      IF(i_rm_min_4los.EQ.1) 
c$$$     &     CALL rread('rm_min_4los',rm_min_4los,1)
c$$$      CALL iread('i_rm_min_2l2cqos',i_rm_min_2l2cqos,1) 
c$$$                    ! minimum invariant mass of 2 leptons and 2 central quark for processes with 2l
c$$$      IF(i_rm_min_2l2cqos.EQ.1) 
c$$$     &     CALL rread('rm_min_2l2cqos',rm_min_2l2cqos,1)
c$$$* sixend
c$$$
c$$$* cuttop
c$$$
c$$$          CALL iread('i_deltacuttopos',i_deltacuttopos,1) 
c$$$                     ! yes/no cut on the invariant mass of all triplets of
c$$$                     ! particles who could form a top 
c$$$                     ! to avoid top contributions
c$$$                     ! The corresponding deltacuttop fixes the 
c$$$                     ! interval excluded : topmas +- deltatacuttop
c$$$          
c$$$          IF(i_deltacuttopos.EQ.1) 
c$$$     &     CALL rread('deltacuttopos',deltacuttopos,1)
c$$$
c$$$* cuttopend
c$$$
c$$$          CALL iread('i_eta_min_jfjbos',i_eta_min_jfjbos,1) 
c$$$          ! minimum rapidity difference between forward and backward jet
c$$$          IF(i_eta_min_jfjbos.EQ.1) 
c$$$     &         CALL rread('eta_min_jfjbos',eta_min_jfjbos,1)
c$$$          CALL iread('i_d_ar_jjos',i_d_ar_jjos,1) 
c$$$                    ! minimum delta_R separation between jets
c$$$          IF(i_d_ar_jjos.EQ.1) CALL rread('d_ar_jjos',d_ar_jjos,1)
c$$$          CALL iread('i_d_ar_jlepos',i_d_ar_jlepos,1)
c$$$                   ! minimum delta_R separation between jets and lepton
c$$$          IF(i_d_ar_jlepos.EQ.1) CALL rread('d_ar_jlepos',d_ar_jlepos,1)
c$$$          CALL iread('i_thetamin_jjos',i_thetamin_jjos,1) 
c$$$                   ! minimum angle separation between jets (cosine)
c$$$          IF(i_thetamin_jjos.EQ.1) 
c$$$     &         CALL rread('thetamin_jjos',thetamin_jjos,1)
c$$$          CALL iread('i_thetamin_jlepos',i_thetamin_jlepos,1) 
c$$$             ! minimum angle separation between jets and lepton (cosine)
c$$$          IF(i_thetamin_jlepos.EQ.1) 
c$$$     &         CALL rread('thetamin_jlepos',thetamin_jlepos,1)
c$$$c giuseppe 01/03/2007
c$$$          CALL iread('i_thetamin_leplepos',i_thetamin_leplepos,1) 
c$$$          IF(i_thetamin_leplepos.EQ.1) 
c$$$     &         CALL rread('thetamin_leplepos',thetamin_leplepos,1)
c$$$c end giuseppe 01/03/2007
c$$$          CALL iread('i_usercutsos',i_usercutsos,1) 
c$$$                   ! yes/no (1/0) additional  user-defined cuts
c$$$        endif

*oldreadend

c the number of files to consider and their names are read in oneshot

      endif

      RETURN
      END







