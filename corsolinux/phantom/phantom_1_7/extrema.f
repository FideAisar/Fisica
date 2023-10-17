c
c Last update: Jul 20, 2007
c

* EXTREMA.F computes pt5min,pt6min,shatmin using the cuts selected 
*   by the user and stores them in the common/min_extrema/
*  
      SUBROUTINE extrema

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT COMPLEX*16 (c)

      include 'common.h'
      include 'common_cut.h'
      include 'common_subproc.h'
c giuseppe 07/07/2007
      DIMENSION pt5min(3),pt6min(3)
c end giuseppe 07/07/2007
      COMMON/fxn_extrema/pt5min,pt6min,shatmin

*      LOCAL VARIABLES
      DIMENSION idpout(6),idpin(2)
         
* EXTRACT FINAL AND INITIAL STATE PARTICLES
* COUNT FINAL STATE JETS, MUONS, ELECTRONS, NEUTRINOS
* COMPUTE SUM OF FINAL AND INITIAL STATE MASSES
      njets=0
      nleptons=0
      nneutrinos=0
      totalmassout=0.d0
      totalmassin=0.d0
      iout=1
      iin=1
      DO i=1,8
        IF(idp_inout(i).GT.0)THEN
          idpout(iout)=idp(i)
          IF(abs(idpout(iout)).LE.5
     &         .or.abs(idpout(iout)).EQ.21) njets=njets+1
          IF((abs(idpout(iout)).EQ.11).OR.(abs(idpout(iout)).EQ.13)
     &        .OR.(abs(idpout(iout)).EQ.15)) nleptons=nleptons+1
          IF(abs(idpout(iout)).EQ.12 .OR.
     &       abs(idpout(iout)).EQ.14) nneutrinos=nneutrinos+1
          totalmassout=totalmassout+rmass(idpout(iout))
          iout=iout+1
        ELSE
          idpin(iin)=idp(i)
          totalmassin=totalmassin+rmass(idpin(iin))
          iin=iin+1
        ENDIF
      ENDDO
      IF(iout.NE.7  .OR. iin.NE.3)THEN
        WRITE(*,*) "WRONG INITIAL/FINAL STATE"
        STOP
      ENDIF 

C   PDG        parton label (6, 5, 4, 3, 2, 1, 0, -1, ......,-5,-6)
C                       for (t, b, c, s, u, d, g, d~, ....., b~,t~)
C   PDG        parton label (11, 12,  13,  14,  15,   16)
C                       for (e-, ve, mu-, vmu, tau, vtau)
      shatmin=0.d0

* MINIMUM MASS SEPARATION BETWEEN JETS
      IF(i_rm_min_jj.EQ.1) THEN
        shatmin=shatmin+njets*(njets-1.d0)/2.d0*rm_min_jj*rm_min_jj
* SUBTRACT EXTRA MASSES: s= SUM(i>j)(2*p_i*p_j)+SUM(i)(m^2_i)
*                        m_ij=2*p_i*p_j+m^2_i+m^2_j
*                        each momentum appears in (njet-1) m_ij
        DO i=1,6
c giuseppe 26/01/2007
c          IF(abs(idpout(i)).LE.5.or.abs(idpout(iout)).EQ.21)THEN
          IF(abs(idpout(i)).LE.5.or.abs(idpout(i)).EQ.21)THEN
c end giuseppe 26/01/2007
            shatmin=shatmin-(njets-2)*rmass2(idpout(i))
          ENDIF
        ENDDO
*      MINIMUM MASS SEPARATION BETWEEN CENTRAL JETS IF rm_min_jcjc > rm_min_jj
        IF(i_rm_min_jcjc.EQ.1) THEN
          IF(rm_min_jcjc.GT.rm_min_jj)THEN
             shatmin=shatmin+
     &        (njets-2.d0)*(njets-3.d0)/2.d0*(rm_min_jcjc*rm_min_jcjc
     &                                   -rm_min_jj*rm_min_jj)
          ENDIF
        ENDIF
*      MINIMUM INVARIANT MASS OF FORWARD AND BACKWARD JETS
        IF(i_rm_min_jfjb.EQ.1) THEN
          IF(rm_min_jfjb.GT.rm_min_jj)THEN
             shatmin=shatmin+rm_min_jfjb*rm_min_jfjb
     &                                  -rm_min_jj*rm_min_jj
          ENDIF
        ENDIF
      ENDIF

* ASK FOR MINIMUM MASS SEPARATION BETWEEN CHARGED LEPTONS AND JETS
* ASSUME ONLY ONE IS PRESENT
      IF(i_rm_min_jlep.EQ.1) THEN
        IF(nleptons.GT.0)THEN
          shatmin=shatmin+njets*rm_min_jlep*rm_min_jlep
        ENDIF
* SUBTRACT EXTRA MASSES IF NECESSARY: Muons considered massless
        IF(i_rm_min_jj.EQ.1)THEN
          DO i=1,6
c giuseppe 26/01/2007
c            IF(abs(idpout(i)).LE.5.or.abs(idpout(iout)).EQ.21)THEN
            IF(abs(idpout(i)).LE.5.or.abs(idpout(i)).EQ.21)THEN
c end giuseppe 26/01/2007
              shatmin=shatmin-rmass2(idpout(i))
            ENDIF
          ENDDO
        ENDIF
      ENDIF

* ADDITIONAL CONSTRAINTS ON shatmin
      emin=nleptons*pt_min_lep+min(nneutrinos,1)*ptmiss_min
      DO i=1,6
c giuseppe 26/01/2007
c        IF(abs(idpout(i)).LE.5.or.abs(idpout(iout)).EQ.21)THEN
        IF(abs(idpout(i)).LE.5.or.abs(idpout(i)).EQ.21)THEN
c end giuseppe 26/01/2007
          emin=emin+sqrt(rmass2(idpout(i))+pt_min_j*pt_min_j)
        ENDIF
      ENDDO
      shatmin=max(shatmin,totalmassout*totalmassout,
     &                 totalmassin*totalmassin,emin*emin)
      

c giuseppe 07/07/2007
*old      pt5min=max(pt_min_lep,ptmiss_min,pt_min_j)
*old      IF(pt_min_lep.GT.0.d0) pt5min=min(pt5min,pt_min_lep)
*old      IF(ptmiss_min.GT.0.d0) pt5min=min(pt5min,ptmiss_min)
*old      IF(pt_min_j.GT.0.d0) pt5min=min(pt5min,pt_min_j)
*old      pt6min=pt5min

c Look for possible t-channels
c giuseppe 20/07/2007
      IF(idpin(1).EQ.21)THEN
        iauxz1=idpin(1)
        iauxw1=0
      ELSE
        iauxz1=-idpin(1)
        IF(mod(abs(idpin(1)),2).EQ.0)THEN ! idpin(1) even: up-type
          IF(idpin(1).GT.0)THEN
            iauxw1=-(idpin(1)-1)
          ELSE                  ! idpin(1).LT.0
            iauxw1=-(idpin(1)+1)
          ENDIF
        ELSE                    ! idpin(1) odd: down-type
          IF(idpin(1).GT.0)THEN
            iauxw1=-(idpin(1)+1)
          ELSE                  ! idpin(1).LT.0
            iauxw1=-(idpin(1)-1)
          ENDIF
        ENDIF
      ENDIF

      IF(idpin(2).EQ.21)THEN
        iauxz2=idpin(2)
        iauxw2=0
      ELSE
        iauxz2=-idpin(2)
        IF(mod(abs(idpin(2)),2).EQ.0)THEN ! idpin(2) even: up-type
          IF(idpin(2).GT.0)THEN
            iauxw2=-(idpin(2)-1)
          ELSE                  ! idpin(2).LT.0
            iauxw2=-(idpin(2)+1)
          ENDIF
        ELSE                    ! idpin(2) odd: down-type
          IF(idpin(2).GT.0)THEN
            iauxw2=-(idpin(2)+1)
          ELSE                  ! idpin(2).LT.0
            iauxw2=-(idpin(2)-1)
          ENDIF
        ENDIF
      ENDIF
c end giuseppe 20/07/2007

c Compute minimum value of transverse momentum allowed for t-channel by 
c cuts on energy and pseudorapidity for jets, charged leptons and 
c neutrinos
      DO i=1,3
        pt5min(i)=0.d0
        pt6min(i)=0.d0
      ENDDO

      ncountert_j=0
      ncountert_lep=0
      ncountert_nu=0
      pt_min_auxj=0.d0
      pt_min_auxlep=0.d0
      ptmiss_min_aux=0.d0

      DO i=1,6
        IF(idpout(i).EQ.iauxz1 .OR. idpout(i).EQ.iauxw1
     &       .OR.  idpout(i).EQ.iauxz2 .OR. idpout(i).EQ.iauxw2)THEN
          IF((abs(idpout(i)).LE.5 .OR. idpout(i).EQ.21) 
     &         .AND. (i_coll.EQ.1 .OR. i_coll.EQ.2) )THEN
            ncountert_j=ncountert_j+1
            IF(i_e_min_j.EQ.1 .AND. i_eta_max_j.EQ.1)THEN
              pt_min_auxj=sqrt(e_min_j*e_min_j-rmass2(idpout(i)))
     &             *sin(2.d0*atan(exp(-eta_max_j)))
              IF(ncountert_j.EQ.1)THEN
                pt5min(1)=pt_min_auxj
              ELSE
                pt5min(1)=min(pt5min(1),pt_min_auxj)
              ENDIF
            ENDIF
          ELSEIF(abs(idpout(i)).EQ.11 .AND. i_coll.EQ.3)THEN
            ncountert_lep=ncountert_lep+1
            IF(i_e_min_lep.EQ.1 .AND. i_eta_max_lep.EQ.1)THEN
              pt_min_auxlep=sqrt(e_min_lep*e_min_lep
     &             -rmass2(idpout(i)))*sin(2.d0
     &             *atan(exp(-eta_max_lep)))
              IF(ncountert_lep.EQ.1)THEN
                pt5min(2)=pt_min_auxlep
              ELSE
                pt5min(2)=min(pt5min(2),pt_min_auxlep)
              ENDIF
            ENDIF
          ELSEIF(abs(idpout(i)).EQ.12 .AND. i_coll.EQ.3)THEN
c no cuts on minimum energy / maximum pseudorapidity can be imposed to 
c neutrinos
            ncountert_nu=ncountert_nu+1
          ELSE
            print*,'***extrema.f ERROR'
            STOP
          ENDIF
        ENDIF
      ENDDO


      IF(i_pt_min_j.EQ.1)THEN
        pt5min(1)=max(pt5min(1),pt_min_j)
      ENDIF

      IF(i_pt_min_lep.EQ.1)THEN
        pt5min(2)=max(pt5min(2),pt_min_lep)
      ENDIF

      IF(i_ptmiss_min.EQ.1)THEN
c Note: ptmiss_min represents the TOTAL missing pT, hence it can be 
c used safely as minimum pT for t-channels only in presence of ONE 
c final-state neutrino.
        IF(ncountert_nu.LT.2)THEN
          pt5min(3)=max(pt5min(3),ptmiss_min)
        ENDIF
      ENDIF

      DO i=1,3
        pt6min(i)=pt5min(i)
      ENDDO
c end giuseppe 07/07/2007

      RETURN
      END


