c
c Last update: Mar 01, 2007
c
* CUT FUNCTION: RETURNS 1 IF ALL CUTS ARE PASSED, 0 IF NOT 

      INTEGER FUNCTION IPASSCUTS(p,
     &  e_min_lep,pt_min_lep,
c giuseppe 01/03/2007
     &  eta_max_onelep,
c end giuseppe 01/03/2007
     &  eta_max_lep,ptmiss_min,
     &  e_min_j,pt_min_j,eta_max_j,
     &  eta_def_jf_min,eta_def_jb_max,eta_def_jc_max,
     &  pt_min_jcjc,rm_min_jj,rm_min_ll,rm_min_jlep,
     &  rm_min_jcjc,rm_max_jcjc,rm_min_jfjb,
* six
     &     rm_min_4l,rm_min_2l2cq,
* YR
     &     rm_max_4l,rm_max_2l2cq,
* YRend
* sixend
* cuttop
     &   deltacuttop, 
* cuttopend
*ww
     &     rcutwlept,rcutptwelectron,
*wwend
*zw
     &     rcutzlept,
*zwend
     &  eta_min_jfjb,
     &  d_ar_jj,d_ar_jlep,d_ar_leplep,
     &  thetamin_jj,thetamin_jlep,
c giuseppe 01/03/2007
     &  thetamin_leplep,
c end giuseppe 01/03/2007
     &  i_e_min_lep,i_pt_min_lep,
     &  i_eta_max_onelep,
     &  i_eta_max_lep,i_ptmiss_min,
     &  i_e_min_j,i_pt_min_j,i_eta_max_j,i_eta_jf_jb_jc,
     &  i_pt_min_jcjc,i_rm_min_jj,i_rm_min_ll,i_rm_min_jlep,
     &  i_rm_min_jcjc,i_rm_max_jcjc,i_rm_min_jfjb,
* six
     &     i_rm_min_4l,i_rm_min_2l2cq,
* YR
     &     i_rm_max_4l,i_rm_max_2l2cq,
* YR end
* sixend
* cuttop
     &     i_deltacuttop, 
* cuttopend
*ww
     &     i_cutwlept,i_cutptwelectron,
*wwend
*zw
     &     i_cutzlept,
*zwend
     &  i_eta_min_jfjb,
     &  i_d_ar_jj,i_d_ar_jlep,i_d_ar_leplep,
     &  i_thetamin_jj,i_thetamin_jlep,
c giuseppe 01/03/2007
     &  i_thetamin_leplep,
c end giuseppe 01/03/2007
     &  i_usercuts, i_usercutsos)

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      include 'common.h'
      INCLUDE 'common_subproc.h'

      DIMENSION p(0:3,8)
* cuttop
      EXTERNAL itoppass
* cuttopend
     

*      LOCAL VARIABLES
      DIMENSION p0(0:3,6),eta0(6),pt0(6),rmsq(6,6),pmod(6),
     &          delphi(6,6),deltheta(6,6),paux(2),pmiss(0:3),
     &          pTWelectron(2)

      DIMENSION idp0(6),invord0(8),
     &          ifjet(6),ifclepton(6),ifneutrino(6)

* six

      DIMENSION nlept(4),nquark(4),qucos(4),nqucentr(2)

* sixend

      ipasscuts=0

      DO i=1,6
        ifjet(i)=0
        ifclepton(i)=0
        ifneutrino(i)=0
      ENDDO
      DO i=0,3
        pmiss(i)=0.d0
      ENDDO

*sandro6/3/07
      nfinallept=0
      nfinalneut=0
*sandro6/3/07end

* six
* YR
c      IF(i_rm_min_4l.eq.1.or.i_rm_min_2l2cq.eq.1)THEN
      IF(i_rm_min_4l.eq.1.or.i_rm_min_2l2cq.eq.1.or.
     &    i_rm_max_4l.eq.1.or.i_rm_max_2l2cq.eq.1)THEN
* YRend
        ntotlept=0
        do i=1,8
          if (idp_inout(i).gt.0.and.abs(idp(i)).gt.10.and.
     &               abs(idp(i)).lt.17)then
            ntotlept=ntotlept+1
            nlept(ntotlept)=i
          endif
        enddo
* YR
c        if (ntotlept.eq.4) then 
        if (ntotlept.eq.4.and.
     &         (i_rm_min_4l.eq.1.or.i_rm_max_4l.eq.1)) then 
* YRend         
          ptot0=p(0,nlept(1))+p(0,nlept(2))
     &         +p(0,nlept(3))+p(0,nlept(4))
          ptot1=p(1,nlept(1))+p(1,nlept(2))
     &         +p(1,nlept(3))+p(1,nlept(4))
          ptot2=p(2,nlept(1))+p(2,nlept(2))
     &         +p(2,nlept(3))+p(2,nlept(4))
          ptot3=p(3,nlept(1))+p(3,nlept(2))
     &         +p(3,nlept(3))+p(3,nlept(4))


          rrmass2=(ptot0**2-ptot1**2-ptot2**2-ptot3**2)
          if (rrmass2.lt.0.d0) RETURN
          
c          Q=sqrt(ptot0**2-ptot1**2-ptot2**2-ptot3**2)

          Q=sqrt(rrmass2)

* YR
c          if (Q.lt.rm_min_4l) RETURN
          if (i_rm_min_4l.eq.1.and.Q.lt.rm_min_4l) RETURN
          if (i_rm_max_4l.eq.1.and.Q.gt.rm_max_4l) RETURN

c        elseif (ntotlept.eq.2) then
        elseif (ntotlept.eq.2.and.
     &         (i_rm_min_2l2cq.eq.1.or.i_rm_max_2l2cq.eq.1)) then
* YRend

* find outgoing quarks or gluons
          ntotquark=0
          do i=1,8
            if (idp_inout(i).gt.0.and.(abs(idp(i)).lt.6.or.abs(idp(i)).
     &           eq.21))then
              ntotquark=ntotquark+1
              nquark(ntotquark)=i
            endif
          enddo
* compute  cosines          
          DO i=1,4
            pmodq=p(1,nquark(i))**2+p(2,nquark(i))**2
     &           +p(3,nquark(i))**2
            pmodq=sqrt(pmodq)
            qucos(i)=p(3,nquark(i))/pmodq
            IF(qucos(i).GT.1.d0.OR.qucos(i).LT.-1.d0) RETURN
          ENDDO
* FIND MOST FORWARD AND MOST BACKWARD JETS
          qucos_max=-1.1d0
          qucos_min=1.1d0
          do i=1,4
            IF(qucos(i).LT.qucos_min)THEN
              qucos_min=qucos(i)
              i_qucos_min=i
            ENDIF
            IF(qucos(i).GT.qucos_max)THEN
              qucos_max=qucos(i)
              i_qucos_max=i
            ENDIF
          enddo
          icentr=1
          do i=1,4
            if (i.ne.i_qucos_min.and.i.ne.i_qucos_max) then
              nqucentr(icentr)=i
              icentr=icentr+1
            endif
          enddo
          ptot0=p(0,nlept(1))+p(0,nlept(2))
     &     +p(0,nquark(nqucentr(1)))+p(0,nquark(nqucentr(2)))
          ptot1=p(1,nlept(1))+p(1,nlept(2))
     &     +p(1,nquark(nqucentr(1)))+p(1,nquark(nqucentr(2)))
          ptot2=p(2,nlept(1))+p(2,nlept(2))
     &     +p(2,nquark(nqucentr(1)))+p(2,nquark(nqucentr(2)))
          ptot3=p(3,nlept(1))+p(3,nlept(2))
     &     +p(3,nquark(nqucentr(1)))+p(3,nquark(nqucentr(2)))

          rrmass2=ptot0**2-ptot1**2-ptot2**2-ptot3**2
          if (rrmass2.lt.0.d0) RETURN

c          Q=sqrt(ptot0**2-ptot1**2-ptot2**2-ptot3**2)

          Q=sqrt(rrmass2)

* YR
c         if (Q.lt.rm_min_2l2cq) RETURN
         if (i_rm_min_2l2cq.eq.1.and.Q.lt.rm_min_2l2cq) RETURN
         if (i_rm_max_2l2cq.eq.1.and.Q.gt.rm_max_2l2cq) RETURN
* YRend 
        endif
      ENDIF

* sixend

      
          
* EXTRACT FINAL STATE PARTICLES AND RECORD THEIR LOCATION IN p0
      i0=1
      DO i=1,8
        IF(idp_inout(i).GT.0)THEN
          idp0(i0)=idp(i)
          DO j=0,3
            p0(j,i0)=p(j,i)
          ENDDO
          invord0(i)=i0
          IF(abs(idp0(i0)).LE.5 .OR. idp0(i0).EQ.21)THEN
            ifjet(i0)=1
          ELSEIF(abs(idp0(i0)).EQ.11 .OR. abs(idp0(i0)).EQ.13 .OR.
     &      abs(idp0(i0)).EQ.15) THEN
            ifclepton(i0)=1
*sandro6/3/07
            nfinallept=nfinallept+1
*sandro6/3/07end
          ELSEIF(abs(idp0(i0)).EQ.12 .OR. abs(idp0(i0)).EQ.14 .OR.
     &      abs(idp0(i0)).EQ.16) THEN
            ifneutrino(i0)=1
            nfinalneut=nfinalneut+1
          ELSE
            WRITE(*,*) 'ERROR IN cuts.f: WRONG IDP:',idp0(i0)
          ENDIF   

          i0=i0+1
        ENDIF
      ENDDO

C   PDG        parton label (6, 5, 4, 3, 2, 1, 0, -1, ......,-5,-6)
C                       for (t, b, c, s, u, d, g, d~, ....., b~,t~)
C   PDG        parton label (11, 12,  13,  14,  15,   16)
C                       for (e-, ve, mu-, vmu, tau, vtau)

* cuttop

      if (i_deltacuttop.eq.1) then
        do i=1,6
          if(idp0(i).eq.5) then
* search for a W+ which can form a top with the b
* e+ ve
            itoppasse=itoppass(idp0,i,-11,12,rmt,deltacuttop,p0)
            if (itoppasse.eq.0) return
* mu+ vmu
            itoppassmu=itoppass(idp0,i,-13,14,rmt,deltacuttop,p0)
            if (itoppassmu.eq.0) return
* u d~
            itoppassu=itoppass(idp0,i,2,-1,rmt,deltacuttop,p0)
            if (itoppassu.eq.0) return
* c s~
            itoppassc=itoppass(idp0,i,4,-3,rmt,deltacuttop,p0)
            if (itoppassc.eq.0) return
          endif
          if(idp0(i).eq.-5) then
* search for a W- which can form a top with the b~
* e- ve~
            itoppasse=itoppass(idp0,i,11,-12,rmt,deltacuttop,p0)
            if (itoppasse.eq.0) return
* mu- vmu~
            itoppassmu=itoppass(idp0,i,13,-14,rmt,deltacuttop,p0)
            if (itoppassmu.eq.0) return
* u~ d
            itoppassu=itoppass(idp0,i,-2,1,rmt,deltacuttop,p0)
            if (itoppassu.eq.0) return
* c~ s
            itoppassc=itoppass(idp0,i,-4,3,rmt,deltacuttop,p0)
            if (itoppassc.eq.0) return
          endif          
        enddo
      endif

* cuttopend


* COMPUTE PT'S AND RAPIDITIES
      DO i=1,6
        part_ptsq=p0(1,i)*p0(1,i)+p0(2,i)*p0(2,i)
        pt0(i)=sqrt(part_ptsq)
        part_3psq=part_ptsq+p0(3,i)*p0(3,i)
        pmod(i)=sqrt(part_3psq)
* elmod
        app=p0(3,i)/pmod(i)
        IF(app.GT.1.d0.OR.app.LT.-1.d0) RETURN
        eta0(i)=-log(tan(0.5d0*acos(app)))
* elmodend
        IF(ifneutrino(i).EQ.1)THEN
          DO j=0,3
            pmiss(j)=pmiss(j)+p0(j,i)
          ENDDO
        ENDIF
      ENDDO

* CHARGED LEPTONS
c giuseppe 06/02/2007
      n_good=0
c end giuseppe 06/02/2007
      DO i=1,6
        IF(ifclepton(i).EQ.1)THEN
* ENERGY
          IF(i_e_min_lep.EQ.1) THEN
            IF(p0(0,i).LT.e_min_lep) RETURN
          ENDIF
* PT
          IF(i_pt_min_lep.EQ.1) THEN
            IF(pt0(i).LT.pt_min_lep) RETURN
          ENDIF
* RAPIDITY
c giuseppe 01/03/2007
          IF(i_eta_max_onelep.EQ.1) THEN
            IF(abs(eta0(i)).LT.eta_max_onelep) n_good=n_good+1
          ENDIF
c end giuseppe 01/03/2007
          IF(i_eta_max_lep.EQ.1) THEN
            IF(abs(eta0(i)).GT.eta_max_lep) RETURN
          ENDIF
        ENDIF
      ENDDO
c giuseppe 06/02/2007
      IF(i_eta_max_onelep.EQ.1) THEN
*sandro6/3/07
        IF(n_good.EQ.0.and.nfinallept.gt.0) RETURN
*sandro6/3/07end
      ENDIF
c end giuseppe 06/02/2007


* NEUTRINOS
* PT
      IF(i_ptmiss_min.EQ.1.and.nfinalneut.ne.0) THEN
        part_ptsq=pmiss(1)*pmiss(1)+pmiss(2)*pmiss(2)
        IF(part_ptsq.LT.ptmiss_min*ptmiss_min) RETURN
      ENDIF

* JETS
      eta_max_j0=-1.d3
      eta_min_j0=1.d3
      DO i=1,6
        IF(ifjet(i).EQ.1)THEN
* ENERGY
          IF(i_e_min_j.EQ.1) THEN
            IF(p0(0,i).LT.e_min_j) RETURN
          ENDIF
* PT
          IF(i_pt_min_j.EQ.1) THEN
            IF(pt0(i).LT.pt_min_j) RETURN
          ENDIF
* RAPIDITY
          IF(i_eta_max_j.EQ.1) THEN
            IF(abs(eta0(i)).GT.eta_max_j) RETURN
          ENDIF
* FIND MOST FORWARD AND MOST BACKWARD JETS
          IF(eta0(i).LT.eta_min_j0)THEN
            eta_min_j0=eta0(i)
            i_eta_min_j0=i
          ENDIF
          IF(eta0(i).GT.eta_max_j0)THEN
            eta_max_j0=eta0(i)
            i_eta_max_j0=i
          ENDIF
        ENDIF
      ENDDO

      IF(i_eta_jf_jb_jc.EQ.1) THEN
* ASK FOR THE MOST BACKWARD JET TO BE BACKWARD
        IF(eta_min_j0.GT.eta_def_jb_max) RETURN
* ASK FOR THE MOST FORWARD JET TO BE FORWARD
        IF(eta_max_j0.LT.eta_def_jf_min) RETURN
* ASK FOR OTHER TWO JETS TO BE CENTRAL
        DO i=1,6 
          IF(ifjet(i).EQ.1 .AND. 
     &       i.NE.i_eta_min_j0 .AND. i.NE.i_eta_max_j0)THEN
            IF(abs(eta0(i)).GT.eta_def_jc_max) RETURN
          ENDIF
        ENDDO
      ENDIF   !  (i_eta_jf_jb_jc.EQ.1)

* ASK FOR MINIMUM RAPIDITY SEPARATION BETWEEN MOST FORWARD AND MOST BACKWARD 
* JET
      IF(i_eta_min_jfjb.EQ.1) THEN
        IF(eta_max_j0-eta_min_j0.LT.eta_min_jfjb) RETURN
      ENDIF
          

* ASK FOR MINIMUM PT OF TWO CENTRAL JETS
      IF(i_pt_min_jcjc.EQ.1)THEN
        paux(1)=0.d0
        paux(2)=0.d0
        DO i=1,6 
          IF(ifjet(i).EQ.1 .AND. 
     &       i.NE.i_eta_min_j0 .AND. i.NE.i_eta_max_j0)THEN
            paux(1)=paux(1)+p0(1,i)
            paux(2)=paux(2)+p0(2,i)
          ENDIF
        ENDDO
        part_ptsq=paux(1)*paux(1)+paux(2)*paux(2)
        IF(part_ptsq.LT.pt_min_jcjc*pt_min_jcjc) RETURN
      ENDIF
      
* COMPUTE INVARIANT MASSES
      DO i=1,6
        DO j=i+1,6
          rmsq(i,j)=(p0(0,i)+p0(0,j))**2
          DO i0=1,3
            rmsq(i,j)=rmsq(i,j)-(p0(i0,i)+p0(i0,j))**2
            if (rmsq(i,j).lt.0.d0) RETURN
          ENDDO
          rmsq(j,i)=rmsq(i,j)
        ENDDO
      ENDDO

*ww
* ASK FOR INVARIANT MASSES OF LEPTONIC W'S BETWEEN rmw+-rcutwlept 
      if(i_cutwlept.eq.1) then
        uprmwlept=rmw+rcutwlept
        rlowrmwlept=rmw-rcutwlept
        do i=1,6
          do k=11,13,2
            if(idp0(i).eq.k)then
              do j=1,6
                if(idp0(j).eq.-(k+1))then
                  rmtmp=sqrt(rmsq(i,j))
                  if (rmtmp.gt.uprmwlept.or.
     &                 rmtmp.lt.rlowrmwlept) then
c                    print*, 1,idp0(i),idp0(j), rmtmp
                    RETURN
c                  else
c                    print*, 1, rmtmp 
                  endif
                endif
              enddo
            endif
            if(idp0(i).eq.-k)then
              do j=1,6
                if(idp0(j).eq.k+1)then
                  rmtmp=sqrt(rmsq(i,j))
                  if (rmtmp.gt.uprmwlept.or.
     &                 rmtmp.lt.rlowrmwlept) then
c                    print*, 2,idp0(i),idp0(j), rmtmp
                    RETURN
c                  else
c                    print*, 2, rmtmp
                  endif
                endif
               enddo
            endif
          enddo
        enddo
      endif
*zw
* ASK FOR INVARIANT MASSES OF LEPTONIC Z'S BETWEEN rmz+-rcutzlept 
      if(i_cutzlept.eq.1) then
        uprmzlept=rmz+rcutzlept
        rlowrmzlept=rmz-rcutzlept
        do i=1,6
          do k=11,13,2
            if(idp0(i).eq.k)then
              do j=1,6
                if(idp0(j).eq.-(k))then
                  rmtmp=sqrt(rmsq(i,j))
                  if (rmtmp.gt.uprmzlept.or.
     &                 rmtmp.lt.rlowrmzlept) then
                    RETURN
c                  else
c                    print*, 1, rmtmp 
                  endif
                endif
              enddo
            endif
          enddo
        enddo
      endif
*zwend
* ASK FOR MINIMUM pT OF ELECTRON-ANTINEUTRINOELECTRON W 
      if(i_cutptwelectron.eq.1) then
        do i=1,6
          do k=11,11
            if(idp0(i).eq.k)then
	      do jj=1,2
	        pTWelectron(jj)=p0(jj,i)
	      enddo
              do j=1,6
                if(idp0(j).eq.-(k+1))then
	          do jj=1,2
	            pTWelectron(jj)=pTWelectron(jj)+p0(jj,j)
	          enddo
		  pTWel=0.d0
		  do jj=1,2
		    pTWel=pTWel+pTWelectron(jj)*pTWelectron(jj)
		  enddo
		  pTWel=sqrt(pTWel)
                  if (pTWel.lt.rcutptwelectron) then
                    RETURN
                  endif
                endif
              enddo
            endif
          enddo
        enddo
      endif
*wwend

* ASK FOR MINIMUM INVARIANT MASS OF FORWARD AND BACKWARD JETS
       IF(i_rm_min_jfjb.EQ.1) THEN
         IF(rmsq(i_eta_min_j0,i_eta_max_j0).LT.
     &                                rm_min_jfjb*rm_min_jfjb)
     &                          RETURN
       ENDIF
       
* ASK FOR MINIMUM MASS SEPARATION BETWEEN JETS
      IF(i_rm_min_jj.EQ.1) THEN
        DO i=1,6
          IF(ifjet(i).EQ.1)THEN
            DO j=i+1,6
              IF(ifjet(j).EQ.1)THEN
                IF(rmsq(i,j).LT.rm_min_jj*rm_min_jj) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      
* ASK FOR MINIMUM MASS SEPARATION BETWEEN CENTRAL JETS
      IF(i_rm_min_jcjc.EQ.1) THEN
        DO i=1,6
          IF(ifjet(i).EQ.1 .AND. 
     &       i.NE.i_eta_min_j0 .AND. i.NE.i_eta_max_j0)THEN
            DO j=i+1,6
              IF(ifjet(j).EQ.1 .AND. 
     &         j.NE.i_eta_min_j0 .AND. j.NE.i_eta_max_j0)THEN
                IF(rmsq(i,j).LT.rm_min_jcjc*rm_min_jcjc) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      
* ASK FOR MAXIMUM MASS SEPARATION BETWEEN CENTRAL JETS
      IF(i_rm_max_jcjc.EQ.1) THEN
        DO i=1,6
          IF(ifjet(i).EQ.1 .AND. 
     &       i.NE.i_eta_min_j0 .AND. i.NE.i_eta_max_j0)THEN
            DO j=i+1,6
              IF(ifjet(j).EQ.1 .AND. 
     &         j.NE.i_eta_min_j0 .AND. j.NE.i_eta_max_j0)THEN
                IF(rmsq(i,j).GT.rm_max_jcjc*rm_max_jcjc) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

* giovanni 2017/10/31 : reintrodotta la modifica
* sandro 25/7 toglie questa modifica
* diogo 18/6
* modification of M(ll) cut for unitarization models to avoid that
* s goes to zero in any kind of process. For this reason, it is 
* asked that any lepton pair, same sign, different family too,
* must be heavier than rm_min_ll.


* ASK FOR MINIMUM MASS SEPARATION FOR OPPOSITE SIGN SAME FLAVOUR 
*        CHARGED LEPTONS
      IF(i_rm_min_ll.EQ.1) THEN
        DO i=1,6
          IF(ifclepton(i).EQ.1)THEN
            DO j=i+1,6
*              IF(ifclepton(j).EQ.1)THEN
              IF(ifclepton(j).EQ.1.AND.idp0(i).eq.-idp0(j))THEN
                IF(rmsq(i,j).LT.rm_min_ll*rm_min_ll) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      
* end diogo      
*end sandro
      
* ASK FOR MINIMUM MASS SEPARATION BETWEEN CHARGED LEPTON AND JETS
      IF(i_rm_min_jlep.EQ.1) THEN
        DO i=1,6
          IF(ifclepton(i).EQ.1)THEN
            DO j=1,6
              IF(ifjet(j).EQ.1)THEN
                IF(rmsq(i,j).LT.rm_min_jlep*rm_min_jlep) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
         

*COMPUTE SEPARATION IN ANGLE AND IN AZIMUTHAL ANGLE PHI
      DO i=1,6
        DO j=i+1,6
          delphi(i,j)=(p0(1,i)*p0(1,j)+p0(2,i)*p0(2,j))/pt0(i)/pt0(j)
          deltheta(i,j)=(p0(1,i)*p0(1,j)+p0(2,i)*p0(2,j)
     &                          +p0(3,i)*p0(3,j))/pmod(i)/pmod(j)
* elmod
          IF(delphi(i,j).GT.1.d0.OR.delphi(i,j).LT.-1.d0) RETURN
          IF(deltheta(i,j).GT.1.d0.OR.deltheta(i,j).LT.-1.d0) RETURN
* elmodend
          delphi(i,j)=acos(delphi(i,j))
          delphi(j,i)=delphi(i,j)
c giuseppe 06/02/2007
c          deltheta(i,j)=acos(deltheta(i,j))
c  sandro per gfortran
c          deltheta(i,j)=acosd(deltheta(i,j))
          deltheta(i,j)=acos(deltheta(i,j))/pi*180.d0
c end giuseppe 06/02/2007
          deltheta(j,i)=deltheta(i,j)
        ENDDO
      ENDDO

* ASK FOR MINIMUM DELTA_R SEPARATION BETWEEN JETS
      IF(i_d_ar_jj.EQ.1) THEN
        DO i=1,6
          IF(ifjet(i).EQ.1)THEN
            DO j=i+1,6
              IF(ifjet(j).EQ.1)THEN
                delr0=delphi(i,j)*delphi(i,j)
     &                          +(eta0(i)-eta0(j))**2
                IF(delr0.LT.d_ar_jj*d_ar_jj) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
         
* ASK FOR MINIMUM DELTA_R SEPARATION BETWEEN JETS AND CHARGED LEPTON
      IF(i_d_ar_jlep.EQ.1) THEN
        DO i=1,6
          IF(ifclepton(i).EQ.1)THEN
            DO j=1,6
              IF(ifjet(j).EQ.1)THEN
                delr0=delphi(i,j)*delphi(i,j)
     &                          +(eta0(i)-eta0(j))**2
                IF(delr0.LT.d_ar_jlep*d_ar_jlep) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

* ASK FOR MINIMUM DELTA_R SEPARATION BETWEEN TWO CHARGED LEPTONS
      IF(i_d_ar_leplep.EQ.1) THEN
c        print*, 'check'
        DO i=1,6
          IF(ifclepton(i).EQ.1)THEN
            DO j=i+1,6
              IF(ifclepton(j).EQ.1)THEN
                delr0=delphi(i,j)*delphi(i,j)
     &                          +(eta0(i)-eta0(j))**2
                IF(delr0.LT.d_ar_leplep*d_ar_leplep) RETURN
ctest                IF(delr0.LT.d_ar_leplep*d_ar_leplep) then
ctest     t             print *, sqrt(delr0)
ctest                  do kkk=1,8
ctest                    print*,idp(kkk)
ctest                    print*,p(0,kkk),p(1,kkk),p(2,kkk),p(3,kkk)
ctest                  enddo
ctest                    print*, '  '
ctest                  RETURN
ctest                endif
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

* ASK FOR MINIMUM ANGULAR SEPARATION BETWEEN JETS
      IF(i_thetamin_jj.EQ.1) THEN
        DO i=1,6
          IF(ifjet(i).EQ.1)THEN
            DO j=i+1,6
              IF(ifjet(j).EQ.1)THEN
                IF(deltheta(i,j).LT.thetamin_jj) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
         
* ASK FOR MINIMUM ANGULAR SEPARATION BETWEEN JETS AND CHARGED LEPTON
      IF(i_thetamin_jlep.EQ.1) THEN
        DO i=1,6
          IF(ifclepton(i).EQ.1)THEN
            DO j=1,6
              IF(ifjet(j).EQ.1)THEN
                IF(deltheta(i,j).LT.thetamin_jlep) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

c giuseppe 01/03/2007
* ASK FOR MINIMUM ANGULAR SEPARATION BETWEEN CHARGED LEPTONS
c giuseppe 26/06/2007
c      IF(i_thetamin_leplep.EQ.1)THEN
c        DO i=1,6
c          IF(ifclepton(i).EQ.1)THEN
c            DO j=1,6
c              IF(ifclepton(j).EQ.1)THEN
c                IF(deltheta(i,j).LT.thetamin_leplep) RETURN
c              ENDIF
c            ENDDO
c          ENDIF
c        ENDDO
c      ENDIF
      IF(i_thetamin_leplep.EQ.1)THEN
        DO i=1,6
          IF(ifclepton(i).EQ.1)THEN
            DO j=i+1,6
              IF(ifclepton(j).EQ.1)THEN
                IF(deltheta(i,j).LT.thetamin_leplep) RETURN
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
c end giuseppe 26/06/2007
c end giuseppe 01/03/2007

* USER SUPPLIED CUTS
      IF(i_usercuts.EQ.1) THEN
        IF(iuserfunc(p0,idp0,eta0,pt0,pmod,rmsq,delphi,deltheta,
     &        i_eta_max_j0,eta_max_j0,i_eta_min_j0,eta_min_j0,
     &        invord0)
     &             .EQ.0) RETURN
      ENDIF

* USER SUPPLIED CUTS in oneshot
      IF(i_usercutsos.EQ.1) THEN
        IF(iuserfuncos(p0,idp0,eta0,pt0,pmod,rmsq,delphi,deltheta,
     &        i_eta_max_j0,eta_max_j0,i_eta_min_j0,eta_min_j0,
     &        invord0)
     &             .EQ.0) RETURN
      ENDIF
      
* IF WE GET HERE ALL CUTS HAVE BEEN PASSED
      ipasscuts=1

      RETURN

      END


* SAMPLE USER CUTS ROUTINE

      INTEGER FUNCTION IUSERFUNC(p0,idp0,eta0,pt0,pmod,rmsq,delphi,
     &        deltheta,i_eta_max_j0,eta_max_j0,i_eta_min_j0,eta_min_j0,
     &        invord0)

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      include 'common.h'
      DIMENSION
     &   p0(0:3,6),eta0(6),pt0(6),pmod(6),rmsq(6,6),
     &  delphi(6,6),deltheta(6,6)
      DIMENSION
     &   idp0(6),invord0(8)
* LOCAL VARIABLES

      iuserfunc=0 
       
c* ASK FOR MINIMUM TRANSVERSE MOMENTUM OF FOUR CENTRAL PARTICLES
c      pperpmin=100.d0
c      px=0.d0
c      py=0.d0
c      DO i=1,6
c        IF(i.NE.i_eta_max_j0 .AND. i.NE.i_eta_min_j0) THEN
c          px=px+p0(1,i)
c          py=py+p0(2,i)
c        ENDIF
c      ENDDO
c      IF(sqrt(px*px+py*py).LT.pperpmin) RETURN

cccccccccccccccccccccccccc diogo 28/04/2010
* ASK FOR VBS cuts. For reactions of type 
* qq-qq munu enu.

*      rM_vv_min=500.d0

*      rM_vv=0.d0
*      rE_vv=0.d0
*      rpx_vv=0.d0
*      rpy_vv=0.d0
*      rpz_vv=0.d0

*      rM_v1=0.d0
*      rE_v1=0.d0
*      rpx_v1=0.d0
*      rpy_v1=0.d0
*      rpz_v1=0.d0

*      rM_v2=0.d0
*      rE_v2=0.d0
*      rpx_v2=0.d0
*      rpy_v2=0.d0
*      rpz_v2=0.d0

*      DO i=1,6
*	IF(abs(idp0(i)).eq.11 .OR. abs(idp0(i)).eq.12)THEN ! elec
*	  rE_vv=rE_vv+p0(0,i)
*	  rpx_vv=rpx_vv+p0(1,i)
*	  rpy_vv=rpy_vv+p0(2,i)
*	  rpz_vv=rpz_vv+p0(3,i)
*	  
*	  rE_v1=rE_v1+p0(0,i)
*	  rpx_v1=rpx_v1+p0(1,i)
*	  rpy_v1=rpy_v1+p0(2,i)
*	  rpz_v1=rpz_v1+p0(3,i)
*	ELSE IF(abs(idp0(i)).eq.13 .OR. abs(idp0(i)).eq.14)THEN !mu 
*	  rE_vv=rE_vv+p0(0,i)
*	  rpx_vv=rpx_vv+p0(1,i)
*	  rpy_vv=rpy_vv+p0(2,i)
*	  rpz_vv=rpz_vv+p0(3,i)
*	  
*	  rE_v2=rE_v2+p0(0,i)
*	  rpx_v2=rpx_v2+p0(1,i)
*	  rpy_v2=rpy_v2+p0(2,i)
*	  rpz_v2=rpz_v2+p0(3,i)	
*	  
*	ENDIF
*      ENDDO

*      rM_vv=sqrt(rE_vv*rE_vv-rpx_vv*rpx_vv-rpy_vv*rpy_vv-rpz_vv*rpz_vv)
*      rM_v1=sqrt(rE_v1*rE_v1-rpx_v1*rpx_v1-rpy_v1*rpy_v1-rpz_v1*rpz_v1)
*      rM_v2=sqrt(rE_v2*rE_v2-rpx_v2*rpx_v2-rpy_v2*rpy_v2-rpz_v2*rpz_v2)
*      
*      IF(rM_vv.LT.rM_vv_min) RETURN
*      IF(rM_v1.LT.70.d0.OR.rM_v1.GT.100.d0) RETURN
*      IF(rM_v2.LT.70.d0.OR.rM_v2.GT.100.d0) RETURN

      
cccccccccccccccccccccccccc end giuseppe 25/07/2007

* IF WE GET HERE ALL USER CUTS HAVE BEEN PASSED
      iuserfunc=1

      RETURN

      END
      
      
* Oneshot second SAMPLE USER CUTS ROUTINE

      INTEGER FUNCTION IUSERFUNCOS(p0,idp0,eta0,pt0,pmod,rmsq,delphi,
     &        deltheta,i_eta_max_j0,eta_max_j0,i_eta_min_j0,eta_min_j0,
     &        invord0)

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      include 'common.h'
      DIMENSION
     &   p0(0:3,6),eta0(6),pt0(6),pmod(6),rmsq(6,6),
     &  delphi(6,6),deltheta(6,6)
      DIMENSION
     &   idp0(6),invord0(8)
* LOCAL VARIABLES

      iuserfuncos=0 
       
* ASK FOR MINIMUM TRANSVERSE MOMENTUM OF FOUR CENTRAL PARTICLES
      pperpmin=100.d0
      px=0.d0
      py=0.d0
      DO i=1,6
        IF(i.NE.i_eta_max_j0 .AND. i.NE.i_eta_min_j0) THEN
          px=px+p0(1,i)
          py=py+p0(2,i)
        ENDIF
      ENDDO
      IF(sqrt(px*px+py*py).LT.pperpmin) RETURN

* IF WE GET HERE ALL USER CUTS HAVE BEEN PASSED
      iuserfuncos=1
      RETURN

      END
      
* cuttop
      INTEGER FUNCTION itoppass(idp,indb,ipw1,ipw2,rmt,deltacuttop,p)
* this function controls if among the 6 external particles there are the two
* ipw1 and ipw2 and if in this case they form with the particle in position
* indb an invariant mass within the interval rmt+-deltacuttop. 
* If this happens itoppass is 0 , otherwise it is 1
* idp(6) are the identities of the particles and p(0:3,6) their momenta

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      DIMENSION p(0:3,6)
      DIMENSION idp(6)
   
      itoppass=0
* search for particle ipw1
      do i=1,6
        if (idp(i).eq.ipw1) then
* search for particle ipw2
          do j=1,6
            if (idp(j).eq.ipw2) then
              rinvmass=(p(0,i)+p(0,j)+p(0,indb))**2-(p(1,i)+p(1,j)+
     &             p(1,indb))**2-(p(2,i)+p(2,j)+p(2,indb))**2
     &             -(p(3,i)+p(3,j)+p(3,indb))**2
              if (rinvmass.lt.0.d0) RETURN
              rinvmass=sqrt(rinvmass)
              rinf=rmt-deltacuttop
              rsup=rmt+deltacuttop
              if(rinvmass.gt.rinf.and.rinvmass.lt.rsup) return
            endif
          enddo
        endif
      enddo        

      itoppass=1
      return

      end
* cuttopend
