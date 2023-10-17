c
c Last update: Oct 30, 2007
c
c Added check that EW resonances are color singlet when decaying hadronically
c
************************************************************************
* SUBROUTINE store_LH_information                                      *
*                                                                      *
* Purpose: store information in Les Houches Event File (LHEF) format   *
************************************************************************

      subroutine store_LH_information(iphsp,idp0,idp0_inout,xm,pcoll,
     &     iordamp,ichargeconj,PDFscale)

      implicit real*8 (a-b,d-h,o-z)

      include 'common.h'
      include 'common_subproc.h'


      integer maxNmothers,maxNparticles
      parameter (maxNmothers=4) ! up to 4 intermediate particles
      parameter (maxNparticles=12)

      integer idp0(8),idp0_inout(8),iordamp(8),index_pythia(8),
     &     ichargeconj,i_mother,isconsistent,ichain,iskip,
     &     i,k,mu,i_outgoing,icolorline,iaux,i_odd,
     &     ismother,index,NUP_aux,IDUP_aux(maxNparticles),
     &     ISTUP_aux(maxNparticles),MOTHUP_aux(2,maxNparticles),
     &     ICOLUP_aux(2,maxNparticles),iord_daughters_LH(6),
     &     itemp(8)

      real *8 xm(8),pcoll(0:3,8),PDFscale,
     &     PUP_aux(5,maxNparticles),VTIMUP_aux(maxNparticles),
     &     SPINUP_aux(maxNparticles)


* Parameters of phase-space multimapping
      COMMON/phspinput_proc/
     &  rm0a(mxnphs),gam0a(mxnphs),spm0alo(mxnphs),spm0ahi(mxnphs),
     &  rm0b(mxnphs),gam0b(mxnphs),spm0blo(mxnphs),spm0bhi(mxnphs),
     &  rm1(mxnphs),gam1(mxnphs),spm1lo(mxnphs),spm1hi(mxnphs),
     &  rm1a(mxnphs),gam1a(mxnphs),spm1alo(mxnphs),spm1ahi(mxnphs),
     &  rm1b(mxnphs),gam1b(mxnphs),spm1blo(mxnphs),spm1bhi(mxnphs),
     &  rm11(mxnphs),gam11(mxnphs),spm11lo(mxnphs),spm11hi(mxnphs),
     &  rm11a(mxnphs),gam11a(mxnphs),spm11alo(mxnphs),spm11ahi(mxnphs),
     &  rm11b(mxnphs),gam11b(mxnphs),spm11blo(mxnphs),spm11bhi(mxnphs),
     &  rm11c(mxnphs),gam11c(mxnphs),spm11clo(mxnphs),spm11chi(mxnphs),
     &  rm12a(mxnphs),gam12a(mxnphs),spm12alo(mxnphs),spm12ahi(mxnphs),
     &  rm12b(mxnphs),gam12b(mxnphs),spm12blo(mxnphs),spm12bhi(mxnphs),
     &  rm12c(mxnphs),gam12c(mxnphs),spm12clo(mxnphs),spm12chi(mxnphs),
     &  rm2(mxnphs),gam2(mxnphs),spm2lo(mxnphs),spm2hi(mxnphs),
     &  rm2a(mxnphs),gam2a(mxnphs),spm2alo(mxnphs),spm2ahi(mxnphs),
     &  rm2b(mxnphs),gam2b(mxnphs),spm2blo(mxnphs),spm2bhi(mxnphs),
     &  rm2c(mxnphs),gam2c(mxnphs),spm2clo(mxnphs),spm2chi(mxnphs),
     &  rm21a(mxnphs),gam21a(mxnphs),spm21alo(mxnphs),spm21ahi(mxnphs),
     &  rm21b(mxnphs),gam21b(mxnphs),spm21blo(mxnphs),spm21bhi(mxnphs),
     &  rm111(mxnphs),gam111(mxnphs),spm111lo(mxnphs),spm111hi(mxnphs),
     &  rm111a(mxnphs),gam111a(mxnphs),spm111alo(mxnphs),
     &     spm111ahi(mxnphs),
     &  rm111b(mxnphs),gam111b(mxnphs),spm111blo(mxnphs),
     &     spm111bhi(mxnphs),
     &  rm1111(mxnphs),gam1111(mxnphs),spm1111lo(mxnphs),
     &     spm1111hi(mxnphs),
     &  rmta(mxnphs),expa(mxnphs),spta,rmtb(mxnphs),expb(mxnphs),sptb

      parameter (maxnresonant=4)
      COMMON/iphsp_fxn/isresonant(maxnresonant),imapregion(maxnresonant)


* Les Houches User Process Event common block
      integer MAXNUP
      parameter (MAXNUP=500)
      integer NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      real*8 XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      common/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &   ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &   VTIMUP(MAXNUP),SPINUP(MAXNUP)





* Store information in Les Houches format (LHEF).
* Incoming particles are expected in the first two positions.

c At most 4 intermediate resonances are possible
      NUP_aux=maxNparticles

c Initialize to zero LHEF commons
      DO i=1,NUP_aux
         IDUP_aux(i)=0
         ISTUP_aux(i)=0
         DO mu=1,5
            PUP_aux(mu,i)=0.d0
         ENDDO
         DO k=1,2
            MOTHUP_aux(k,i)=0
            ICOLUP_aux(k,i)=0
         ENDDO
         VTIMUP_aux(i)=0.d0
         SPINUP_aux(i)=0.d0
      ENDDO


c The LHEF list of final-state particles is intended to begin with
c (possible) mothers.
c First, assign Les Houches commons to EXTERNAL final-state particles.
c LH information about mothers will be reconstructed by their daughters.

      i_outgoing = 3 + maxNmothers
      DO i=1,8
         IF(idp0_inout(i).EQ.-1)THEN
            IF(pcoll(3,i).GT.0.d0)THEN
               index_pythia(i)=1
            ELSE
               index_pythia(i)=2
            ENDIF
         ELSE
            index_pythia(i)=i_outgoing
            i_outgoing=i_outgoing+1
         ENDIF
      ENDDO

c Assign colour flow to all external particles
      icolorline=500

c     all momenta outgoing, fermionic flow from right to left
c     q ----<-(g,gg)----q~
c     start from the quark, moving against fermionic flow
      i_odd=1


      DO i=1,8
         IF(idp0(i).LE.16)THEN
            IDUP_aux(index_pythia(i))=idp0(i)*idp0_inout(i)
         ELSE
            IDUP_aux(index_pythia(i))=idp0(i)
         ENDIF
         ISTUP_aux(index_pythia(i))=idp0_inout(i)
         DO j=1,3
            PUP_aux(j,index_pythia(i))=pcoll(j,i)
         ENDDO
         PUP_aux(4,index_pythia(i))=pcoll(0,i)
c Ezio 16/06/2018
c         PUP_aux(5,index_pythia(i))=xm(i)
         PUP_aux(5,index_pythia(i))=rmass(IDUP_aux(index_pythia(i)))

         IF(abs(idp0(i)).le.5)THEN ! quark or antiquark

            IF(i_odd.EQ.1)THEN  ! quark, if outgoing
               icolorline=icolorline+1
               IF(idp0_inout(i).EQ.1)THEN ! = (q)-
                  ICOLUP_aux(1,index_pythia(i))=icolorline
                  ICOLUP_aux(2,index_pythia(i))=0
               ELSE             ! = (q~)-
                  ICOLUP_aux(1,index_pythia(i))=0
                  ICOLUP_aux(2,index_pythia(i))=icolorline
               ENDIF
               i_odd=0          !next fermion will be an antiquark

            ELSE                ! antiquark, if outgoing
               IF(idp0_inout(i).EQ.1)THEN ! = -(q~)
                  ICOLUP_aux(1,index_pythia(i))=0
                  ICOLUP_aux(2,index_pythia(i))=icolorline
               ELSE             ! = -(q)
                  ICOLUP_aux(1,index_pythia(i))=icolorline
                  ICOLUP_aux(2,index_pythia(i))=0
               ENDIF
               i_odd=1          !next fermion will be a quark
            ENDIF

         ELSEIF(idp0(i).EQ.21)THEN ! gluon = -(q~)|(q)- if outgoing

            IF(idp0_inout(i).EQ.1)THEN ! = -(q~)|(q)-
               ICOLUP_aux(2,index_pythia(i))=icolorline
               icolorline=icolorline+1
               ICOLUP_aux(1,index_pythia(i))=icolorline
            ELSE                ! = -(q)|(q~)-
               ICOLUP_aux(1,index_pythia(i))=icolorline
               icolorline=icolorline+1
               ICOLUP_aux(2,index_pythia(i))=icolorline
            ENDIF

         ELSE                   ! lepton
            ICOLUP_aux(1,index_pythia(i))=0
            ICOLUP_aux(2,index_pythia(i))=0
         ENDIF
      ENDDO

      IF (ichargeconj.eq.1) THEN
         DO i=1,8
            iaux=ICOLUP_aux(1,index_pythia(i))
            ICOLUP_aux(1,index_pythia(i))=ICOLUP_aux(2,index_pythia(i))
            ICOLUP_aux(2,index_pythia(i))=iaux
         ENDDO
      ENDIF

c Start by assigning the relationship of all external particles to the
c initial state (i.e. 1,2 as mothers). This relationship will be
c eventually updated later on if intermediate resonances (mothers) are
c found.
      DO k=1,2
c incoming particles have no mother
         MOTHUP_aux(1,k)=0
         MOTHUP_aux(2,k)=0
      ENDDO
      DO k=3+maxNmothers,NUP_aux
         MOTHUP_aux(1,k)=1
         MOTHUP_aux(2,k)=2
      ENDDO

      DO k=3+maxNmothers,NUP_aux
         SPINUP_aux(k)=9.d0
         VTIMUP_aux(k)=0.d0
      ENDDO

* Complete Les Houches commons
      SCALUP=PDFscale
      AQEDUP=1.d0/alfainv
      AQCDUP=alfa_s
      IDPRUP=661
      XWGTUP=1.d0



      if(iwrite_mothers.eq.1)then

* Now look for possible intermediate particles (interpreted as
* resonances in multiparticle invariant mass). The identification
* procedure is stricly dependent on the phase-space routine called.

      i_mother=2
      ichain=1
      iskip=0
      isconsistent=1

      DO i=1,8
        itemp(i)=index_pythia(iordamp(iorder(i,iphs_ind)))
      ENDDO

      if(iphsp.eq.1)then        ! phsp1_1_4

        igranmoth1=1
        igranmoth2=2
        if(isresonant(1).eq.1)then ! rm1[a,b]
          i_mother=i_mother+1
          ndaughters=4
          iord_daughters_LH(1)=itemp(3)
          iord_daughters_LH(2)=itemp(4)
          iord_daughters_LH(3)=itemp(5)
          iord_daughters_LH(4)=itemp(6)
          if(imapregion(1).eq.2)then
              rmass_mapping=rm1a(iphs_ind)
          elseif(imapregion(1).eq.4)then
              rmass_mapping=rm1b(iphs_ind)
          endif

          call add_LH_mother2(ndaughters,iord_daughters_LH,
     &          rmass_mapping,i_mother,igranmoth1,igranmoth2,
     &          IDUP_aux,ISTUP_aux,
     &          MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &          isconsistent)

          igranmoth1=i_mother
          igranmoth2=i_mother
        endif


c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(3))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(3)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(3))
      ELSE
       ico1=ICOLUP_aux(2,itemp(3))
      ENDIF

      IF(IDUP_aux(itemp(4)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(4))
      ELSE
       ico2=ICOLUP_aux(2,itemp(4))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(3).eq.1 .AND. issinglet.EQ.1)then ! rm11[a,b,c]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1)=itemp(3)
            iord_daughters_LH(2)=itemp(4)
            if(imapregion(3).eq.2)then
               rmass_mapping=rm11a(iphs_ind)
            elseif(imapregion(3).eq.4)then
               rmass_mapping=rm11b(iphs_ind)
            elseif(imapregion(3).eq.6)then
               rmass_mapping=rm11c(iphs_ind)
            endif

            call add_LH_mother2(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,igranmoth1,igranmoth2,
     &           IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(5))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(5)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(5))
      ELSE
       ico1=ICOLUP_aux(2,itemp(5))
      ENDIF

      IF(IDUP_aux(itemp(6)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(6))
      ELSE
       ico2=ICOLUP_aux(2,itemp(6))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(4).eq.1 .AND. issinglet.EQ.1)then ! rm12[a,b,c]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1)=itemp(5)
            iord_daughters_LH(2)=itemp(6)
            if(imapregion(4).eq.2)then
               rmass_mapping=rm12a(iphs_ind)
            elseif(imapregion(4).eq.4)then
               rmass_mapping=rm12b(iphs_ind)
            elseif(imapregion(4).eq.6)then
               rmass_mapping=rm12c(iphs_ind)
            endif

            call add_LH_mother2(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,igranmoth1,igranmoth2,
     &           IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif


      elseif(iphsp.eq.2)then    ! phsp2_4

c      write(23,'(A12,4(4I,1X),A13,4I)') 'isresonant: ',
c     &              (isresonant(k),k=1,4),'  iphs_ind : ',iphs_ind
cccc      Write(23,'(A9,8(I4,1X))')'iorder  :',(iorder(k,iphs_ind),k=1,8)
        igranmoth1=1
        igranmoth2=2
        if(isresonant(1).eq.1)then ! rm1
          i_mother=i_mother+1
          ndaughters=4
          iord_daughters_LH(1)=itemp(3)
          iord_daughters_LH(2)=itemp(4)
          iord_daughters_LH(3)=itemp(5)
          iord_daughters_LH(4)=itemp(6)
          if(imapregion(1).eq.2)then
             rmass_mapping=rm1a(iphs_ind)
          elseif(imapregion(1).eq.4)then
             rmass_mapping=rm1b(iphs_ind)
          endif

          call add_LH_mother2(ndaughters,iord_daughters_LH,
     &          rmass_mapping,i_mother,igranmoth1,igranmoth2,
     &          IDUP_aux,ISTUP_aux,
     &          MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &          isconsistent)

          igranmoth1=i_mother
          igranmoth2=i_mother
         endif

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(3))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(3)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(3))
      ELSE
       ico1=ICOLUP_aux(2,itemp(3))
      ENDIF

      IF(IDUP_aux(itemp(4)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(4))
      ELSE
       ico2=ICOLUP_aux(2,itemp(4))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(3).eq.1 .AND. issinglet.EQ.1)then ! rm11[a,b,c]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1)=itemp(3)
            iord_daughters_LH(2)=itemp(4)
            if(imapregion(3).eq.2)then
               rmass_mapping=rm11a(iphs_ind)
            elseif(imapregion(3).eq.4)then
               rmass_mapping=rm11b(iphs_ind)
            elseif(imapregion(3).eq.6)then
               rmass_mapping=rm11c(iphs_ind)
            endif

            call add_LH_mother2(ndaughters,iord_daughters_LH,
     &          rmass_mapping,i_mother,igranmoth1,igranmoth2,
     &          IDUP_aux,ISTUP_aux,
     &          MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &          isconsistent)
          endif

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(5))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(5)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(5))
      ELSE
       ico1=ICOLUP_aux(2,itemp(5))
      ENDIF

      IF(IDUP_aux(itemp(6)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(6))
      ELSE
       ico2=ICOLUP_aux(2,itemp(6))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(4).eq.1 .AND. issinglet.EQ.1)then ! rm12[a,b]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1)=itemp(5)
            iord_daughters_LH(2)=itemp(6)
            if(imapregion(4).eq.2)then
               rmass_mapping=rm12a(iphs_ind)
            elseif(imapregion(4).eq.4)then
               rmass_mapping=rm12b(iphs_ind)
            elseif(imapregion(4).eq.6)then
               rmass_mapping=rm12c(iphs_ind)
            endif

            call add_LH_mother2(ndaughters,iord_daughters_LH,
     &          rmass_mapping,i_mother,igranmoth1,igranmoth2,
     &          IDUP_aux,ISTUP_aux,
     &          MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &          isconsistent)
          endif

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(7))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(7)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(7))
      ELSE
       ico1=ICOLUP_aux(2,itemp(7))
      ENDIF

      IF(IDUP_aux(itemp(8)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(8))
      ELSE
       ico2=ICOLUP_aux(2,itemp(8))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(2).eq.1 .AND. issinglet.EQ.1)then ! rm2[a,b,c]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1)=itemp(7)
            iord_daughters_LH(2)=itemp(8)
            if(imapregion(2).eq.2)then
               rmass_mapping=rm2a(iphs_ind)
            elseif(imapregion(2).eq.4)then
               rmass_mapping=rm2b(iphs_ind)
            elseif(imapregion(2).eq.6)then
               rmass_mapping=rm2c(iphs_ind)
            endif

            call add_LH_mother(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif

      elseif(iphsp.eq.3)then    ! phsp1_1_31

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(3))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(3)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(3))
      ELSE
       ico1=ICOLUP_aux(2,itemp(3))
      ENDIF

      IF(IDUP_aux(itemp(4)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(4))
      ELSE
       ico2=ICOLUP_aux(2,itemp(4))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(3).eq.1 .AND. issinglet.EQ.1)then ! rm111[a,b]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1) = itemp(3)
            iord_daughters_LH(2) = itemp(4)
            if(imapregion(3).eq.2)then
               rmass_mapping=rm111a(iphs_ind)
            elseif(imapregion(3).eq.4)then
               rmass_mapping=rm111b(iphs_ind)
            endif

            call add_LH_mother(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons (e.g. both {q-qbar} and
c {q-qbar-g} resonant in the W/Z range)
         if(rm111a(iphs_ind).eq.rm11a(iphs_ind) .and.
     &        rm111b(iphs_ind).eq.rm11b(iphs_ind))then
           ichain=ichain+1
         else
           ichain=1
           iskip=0
         endif
         if(ichain.gt.1 .and.
     &        (isresonant(3).eq.1 .and. isconsistent.eq.1))then
           iskip=1
         endif
c giuseppe 08/09/2007
         if(isresonant(3).ne.1 .or. isconsistent.eq.0)then
           iskip=1
         endif
c end giuseppe 08/09/2007

         if(isresonant(2).eq.1)then ! rm11[a,b]
           i_mother=i_mother+1
           if(iskip.eq.0)then
             ndaughters=3
             iord_daughters_LH(1) =  itemp(3)
             iord_daughters_LH(2) =  itemp(4)
             iord_daughters_LH(3) =  itemp(5)
             if(imapregion(2).eq.2)then
               rmass_mapping=rm11a(iphs_ind)
             elseif(imapregion(2).eq.4)then
               rmass_mapping=rm11b(iphs_ind)
             endif

             call add_LH_mother(ndaughters,iord_daughters_LH,
     &            rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &            MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &            isconsistent)
           endif
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons
         if(rm11a(iphs_ind).eq.rm1a(iphs_ind) .and.
     &        rm11b(iphs_ind).eq.rm1b(iphs_ind))then
           ichain=ichain+1
         else
           ichain=1
c giuseppe 08/09/2007
c            iskip=0
c end giuseppe 08/09/2007
         endif
         if(ichain.gt.1 .and. (iskip.eq.1 .or.
     &        (isresonant(2).eq.1 .and. isconsistent.eq.1)))then
           iskip=1
         endif

         if(isresonant(1).eq.1)then ! rm1[a,b]
           i_mother=i_mother+1
           if(iskip.eq.0)then
             ndaughters=4
             iord_daughters_LH(1) =  itemp(3)
             iord_daughters_LH(2) =  itemp(4)
             iord_daughters_LH(3) =  itemp(5)
             iord_daughters_LH(4) =  itemp(6)
             if(imapregion(1).eq.2)then
               rmass_mapping=rm1a(iphs_ind)
             elseif(imapregion(1).eq.4)then
               rmass_mapping=rm1b(iphs_ind)
             endif

             call add_LH_mother(ndaughters,iord_daughters_LH,
     &            rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &            MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &            isconsistent)
           endif
         endif



      elseif(iphsp.eq.4)then    ! phsp1_2_3

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(3))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(3)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(3))
      ELSE
       ico1=ICOLUP_aux(2,itemp(3))
      ENDIF

      IF(IDUP_aux(itemp(4)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(4))
      ELSE
       ico2=ICOLUP_aux(2,itemp(4))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(4).eq.1 .AND. issinglet.EQ.1)then ! rm11[a,b]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1) = itemp(3)
            iord_daughters_LH(2) = itemp(4)
           if(imapregion(4).eq.2)then
               rmass_mapping=rm11a(iphs_ind)
            elseif(imapregion(4).eq.4)then
               rmass_mapping=rm11b(iphs_ind)
            endif

            call add_LH_mother(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons
         if(rm11a(iphs_ind).eq.rm1a(iphs_ind) .and.
     &        rm11b(iphs_ind).eq.rm1b(iphs_ind))then
            ichain=ichain+1
         else
            ichain=1
            iskip=0
         endif
         if(ichain.gt.1 .and.
     &        (isresonant(4).eq.1 .and. isconsistent.eq.1))then
            iskip=1
         endif
c giuseppe 08/09/2007
         if(isresonant(4).ne.1 .or. isconsistent.eq.0)then
           iskip=1
         endif
c end giuseppe 08/09/2007

         if(isresonant(2).eq.1)then ! rm1[a,b]
            i_mother=i_mother+1
            if(iskip.eq.0)then
               ndaughters=3
               iord_daughters_LH(1) = itemp(3)
               iord_daughters_LH(2) = itemp(4)
               iord_daughters_LH(3) = itemp(5)
               if(imapregion(2).eq.2)then
                  rmass_mapping=rm1a(iphs_ind)
               elseif(imapregion(2).eq.4)then
                  rmass_mapping=rm1b(iphs_ind)
               endif

               call add_LH_mother(ndaughters,iord_daughters_LH,
     &              rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &              MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &              isconsistent)
            endif
         endif

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(7))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(7)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(7))
      ELSE
       ico1=ICOLUP_aux(2,itemp(7))
      ENDIF

      IF(IDUP_aux(itemp(6)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(6))
      ELSE
       ico2=ICOLUP_aux(2,itemp(6))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(3).eq.1 .AND. issinglet.EQ.1)then ! rm2[a,b]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1) = itemp(6)
            iord_daughters_LH(2) = itemp(7)
           if(imapregion(3).eq.2)then
               rmass_mapping=rm2a(iphs_ind)
            elseif(imapregion(3).eq.4)then
               rmass_mapping=rm2b(iphs_ind)
            endif

            call add_LH_mother(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif

c giuseppe 08/09/2007
         if(isresonant(3).ne.1 .or. isconsistent.eq.0)then
           iskip=1
         endif
c end giuseppe 08/09/2007

         if(isresonant(1).eq.1)then ! rm0[a,b]
           i_mother=i_mother+1
c giuseppe 08/09/2007
           if(iskip.eq.0)then
c end giuseppe 08/09/2007
             ndaughters=5
             iord_daughters_LH(1) = itemp(3)
             iord_daughters_LH(2) = itemp(4)
             iord_daughters_LH(3) = itemp(5)
             iord_daughters_LH(4) = itemp(6)
             iord_daughters_LH(5) = itemp(7)
             if(imapregion(1).eq.2)then
               rmass_mapping=rm0a(iphs_ind)
             elseif(imapregion(1).eq.4)then
               rmass_mapping=rm0b(iphs_ind)
             endif

             call add_LH_mother(ndaughters,iord_daughters_LH,
     &            rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &            MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &            isconsistent)
c giuseppe 08/09/2007
           endif
c end giuseppe 08/09/2007
         endif


      elseif(iphsp.eq.5)then    ! phsp3_3

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(3))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(3)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(3))
      ELSE
       ico1=ICOLUP_aux(2,itemp(3))
      ENDIF

      IF(IDUP_aux(itemp(4)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(4))
      ELSE
       ico2=ICOLUP_aux(2,itemp(4))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(3).eq.1 .AND. issinglet.EQ.1)then ! rm11[a,b]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1) = itemp(3)
            iord_daughters_LH(2) = itemp(4)
            if(imapregion(3).eq.2)then
               rmass_mapping=rm11a(iphs_ind)
            elseif(imapregion(3).eq.4)then
               rmass_mapping=rm11b(iphs_ind)
            endif

            call add_LH_mother(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons
         if(rm11a(iphs_ind).eq.rm1a(iphs_ind) .and.
     &        rm11b(iphs_ind).eq.rm1b(iphs_ind))then
           ichain=ichain+1
         else
           ichain=1
           iskip=0
         endif
         if(ichain.gt.1 .and.
     &        (isresonant(3).eq.1 .and. isconsistent.eq.1))then
           iskip=1
         endif
c giuseppe 08/09/2007
         if(isresonant(3).ne.1 .or. isconsistent.eq.0)then
           iskip=1
         endif
c end giuseppe 08/09/2007

         if(isresonant(1).eq.1)then ! rm1[a,b]
           i_mother=i_mother+1
           if(iskip.eq.0)then
             ndaughters=3
             iord_daughters_LH(1) = itemp(3)
             iord_daughters_LH(2) = itemp(4)
             iord_daughters_LH(3) = itemp(7)
             if(imapregion(1).eq.2)then
               rmass_mapping=rm1a(iphs_ind)
             elseif(imapregion(1).eq.4)then
               rmass_mapping=rm1b(iphs_ind)
             endif

             call add_LH_mother(ndaughters,iord_daughters_LH,
     &            rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &            MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &            isconsistent)
           endif
         endif

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(5))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(5)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(5))
      ELSE
       ico1=ICOLUP_aux(2,itemp(5))
      ENDIF

      IF(IDUP_aux(itemp(6)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(6))
      ELSE
       ico2=ICOLUP_aux(2,itemp(6))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(4).eq.1 .AND. issinglet.EQ.1)then ! rm21[a,b]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1) = itemp(5)
            iord_daughters_LH(2) = itemp(6)
            if(imapregion(4).eq.2)then
               rmass_mapping=rm21a(iphs_ind)
            elseif(imapregion(4).eq.4)then
               rmass_mapping=rm21b(iphs_ind)
            endif

            call add_LH_mother(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons
         ichain=1
         iskip=0
         if(rm21a(iphs_ind).eq.rm2a(iphs_ind) .and.
     &        rm21b(iphs_ind).eq.rm2b(iphs_ind))then
           ichain=ichain+1
         else
           ichain=1
           iskip=0
         endif
         if(ichain.gt.1 .and.
     &        (isresonant(4).eq.1  .and. isconsistent.eq.1))then
           iskip=1
         endif
c giuseppe 08/09/2007
         if(isresonant(4).ne.1 .or. isconsistent.eq.0)then
           iskip=1
         endif
c end giuseppe 08/09/2007

         if(isresonant(2).eq.1)then ! rm2[a,b]
            i_mother=i_mother+1
            if(iskip.eq.0)then
               ndaughters=3
               iord_daughters_LH(1) = itemp(5)
               iord_daughters_LH(2) = itemp(6)
               iord_daughters_LH(3) = itemp(8)
               if(imapregion(2).eq.2)then
                  rmass_mapping=rm2a(iphs_ind)
               elseif(imapregion(2).eq.4)then
                  rmass_mapping=rm2b(iphs_ind)
               endif

               call add_LH_mother(ndaughters,iord_daughters_LH,
     &              rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &              MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &              isconsistent)
            endif
         endif



      elseif(iphsp.eq.6)then    ! phsp2_4to31

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(3))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(3)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(3))
      ELSE
       ico1=ICOLUP_aux(2,itemp(3))
      ENDIF

      IF(IDUP_aux(itemp(4)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(4))
      ELSE
       ico2=ICOLUP_aux(2,itemp(4))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(4).eq.1)then ! rm111[a,b]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1) = itemp(3)
            iord_daughters_LH(2) = itemp(4)
            if(imapregion(4).eq.2)then
               rmass_mapping=rm111a(iphs_ind)
            elseif(imapregion(4).eq.4)then
               rmass_mapping=rm111b(iphs_ind)
            endif

            call add_LH_mother(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons
         if(rm111a(iphs_ind).eq.rm11a(iphs_ind) .and.
     &        rm111b(iphs_ind).eq.rm11b(iphs_ind))then
           ichain=ichain+1
         else
           ichain=1
           iskip=0
         endif
         if(ichain.gt.1 .and.
     &        (isresonant(4).eq.1 .and. isconsistent.eq.1))then
           iskip=1
         endif
c giuseppe 08/09/2007
         if(isresonant(4).ne.1 .or. isconsistent.eq.0)then
           iskip=1
         endif
c end giuseppe 08/09/2007

         if(isresonant(3).eq.1)then ! rm11[a,b]
            i_mother=i_mother+1
            if(iskip.eq.0)then
               ndaughters=3
               iord_daughters_LH(1) = itemp(3)
               iord_daughters_LH(2) = itemp(4)
               iord_daughters_LH(3) = itemp(5)
               if(imapregion(3).eq.2)then
                  rmass_mapping=rm11a(iphs_ind)
               elseif(imapregion(3).eq.4)then
                  rmass_mapping=rm11b(iphs_ind)
               endif

               call add_LH_mother(ndaughters,iord_daughters_LH,
     &              rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &              MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &              isconsistent)
            endif
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons
         if(rm11a(iphs_ind).eq.rm1a(iphs_ind) .and.
     &        rm11b(iphs_ind).eq.rm1b(iphs_ind))then
            ichain=ichain+1
         else
            ichain=1
c giuseppe 08/09/2007
c            iskip=0
c end giuseppe 08/09/2007
         endif
         if(ichain.gt.1 .and. (iskip.eq.1 .or.
     &        (isresonant(3).eq.1 .and. isconsistent.eq.1)))then
           iskip=1
         endif

         if(isresonant(1).eq.1)then ! rm1[a,b]
            i_mother=i_mother+1
            if(iskip.eq.0)then
               ndaughters=4
               iord_daughters_LH(1) = itemp(3)
               iord_daughters_LH(2) = itemp(4)
               iord_daughters_LH(3) = itemp(5)
               iord_daughters_LH(4) = itemp(6)
               if(imapregion(1).eq.2)then
                  rmass_mapping=rm1a(iphs_ind)
               elseif(imapregion(1).eq.4)then
                  rmass_mapping=rm1b(iphs_ind)
               endif

               call add_LH_mother(ndaughters,iord_daughters_LH,
     &              rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &              MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &              isconsistent)
            endif
         endif

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(7))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(7)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(7))
      ELSE
       ico1=ICOLUP_aux(2,itemp(7))
      ENDIF

      IF(IDUP_aux(itemp(8)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(8))
      ELSE
       ico2=ICOLUP_aux(2,itemp(8))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(2).eq.1 .AND. issinglet.EQ.1)then ! rm2[a,b]
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1) = itemp(7)
            iord_daughters_LH(2) = itemp(8)
            if(imapregion(2).eq.2)then
               rmass_mapping=rm2a(iphs_ind)
            elseif(imapregion(2).eq.4)then
               rmass_mapping=rm2b(iphs_ind)
            endif

            call add_LH_mother(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif



      elseif(iphsp.eq.7)then    ! phsp1_5to1_4to31

c Check if colorflow is consistent with presence of a resonance
      issinglet=1
      IF(abs(IDUP_aux(itemp(3))).le.5)THEN ! quark or antiquark

      IF(IDUP_aux(itemp(3)).GT.0)THEN 
       ico1=ICOLUP_aux(1,itemp(3))
      ELSE
       ico1=ICOLUP_aux(2,itemp(3))
      ENDIF

      IF(IDUP_aux(itemp(4)).GT.0)THEN 
       ico2=ICOLUP_aux(1,itemp(4))
      ELSE
       ico2=ICOLUP_aux(2,itemp(4))
      ENDIF

      IF(ico1.NE.ico2) issinglet=0

      ENDIF

         if(isresonant(4).eq.1 .AND. issinglet.EQ.1)then ! rm1111
            i_mother=i_mother+1
            ndaughters=2
            iord_daughters_LH(1) = itemp(3)
            iord_daughters_LH(2) = itemp(4)
            rmass_mapping=rm1111(iphs_ind)

            call add_LH_mother(ndaughters,iord_daughters_LH,
     &           rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &           MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &           isconsistent)
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons
         if(rm1111(iphs_ind).eq.rm111(iphs_ind))then
            ichain=ichain+1
         else
            ichain=1
            iskip=0
         endif
         if(ichain.gt.1 .and.
     &        (isresonant(4).eq.1 .and. isconsistent.eq.1))then
            iskip=1
         endif
c giuseppe 08/09/2007
         if(isresonant(4).ne.1 .or. isconsistent.eq.0)then
           iskip=1
         endif
c end giuseppe 08/09/2007

         if(isresonant(3).eq.1)then ! rm111
            i_mother=i_mother+1
            if(iskip.eq.0)then
               ndaughters=3
               iord_daughters_LH(1) = itemp(3)
               iord_daughters_LH(2) = itemp(4)
               iord_daughters_LH(3) = itemp(5)
               rmass_mapping=rm111(iphs_ind)

               call add_LH_mother(ndaughters,iord_daughters_LH,
     &              rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &              MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &              isconsistent)
            endif
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons
         if(rm111(iphs_ind).eq.rm11(iphs_ind))then
           ichain=ichain+1
         else
           ichain=1
c giuseppe 08/09/2007
c            iskip=0
c end giuseppe 08/09/2007
         endif
         if(ichain.gt.1 .and. (iskip.eq.1 .or.
     &        (isresonant(3).eq.1 .and. isconsistent.eq.1)))then
           iskip=1
         endif

         if(isresonant(2).eq.1)then ! rm11
            i_mother=i_mother+1
            if(iskip.eq.0)then
               ndaughters=4
               iord_daughters_LH(1) = itemp(3)
               iord_daughters_LH(2) = itemp(4)
               iord_daughters_LH(3) = itemp(5)
               iord_daughters_LH(4) = itemp(6)
              rmass_mapping=rm11(iphs_ind)

               call add_LH_mother(ndaughters,iord_daughters_LH,
     &              rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &              MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &              isconsistent)
            endif
         endif

c Avoid double counting of intermediate particles due to possible
c multiple resonances in presence of gluons
         if(rm11(iphs_ind).eq.rm1(iphs_ind))then
            ichain=ichain+1
         else
            ichain=1
c giuseppe 08/09/2007
c            iskip=0
c end giuseppe 08/09/2007
         endif
         if(ichain.gt.1 .and. (iskip.eq.1 .or.
     &        (isresonant(2).eq.1 .and. isconsistent.eq.1)))then
            iskip=1
         endif

         if(isresonant(1).eq.1)then ! rm1
            i_mother=i_mother+1
            if(iskip.eq.0)then
               ndaughters=5
               iord_daughters_LH(1) = itemp(3)
               iord_daughters_LH(2) = itemp(4)
               iord_daughters_LH(3) = itemp(5)
               iord_daughters_LH(4) = itemp(6)
               iord_daughters_LH(5) = itemp(7)
               rmass_mapping=rm1(iphs_ind)

               call add_LH_mother(ndaughters,iord_daughters_LH,
     &              rmass_mapping,i_mother,IDUP_aux,ISTUP_aux,
     &              MOTHUP_aux,PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux,
     &              isconsistent)
            endif
         endif

      endif                     ! if(iphsp.eq.1)then

      endif                     ! if(iwrite_mothers.eq.1)then



* Finally initialize Les Houches commons.

      if(iwrite_mothers.eq.1)then

         NUP=maxNparticles

         do i=1,NUP_aux
            do while(IDUP_aux(i).eq.0)
               NUP=NUP-1
               call shift_LH_stack(i,IDUP_aux,ISTUP_aux,MOTHUP_aux,
     &              PUP_aux,ICOLUP_aux,VTIMUP_aux,SPINUP_aux)
            enddo
         enddo

         do i=1,NUP
            IDUP(i)=IDUP_aux(i)
            ISTUP(i)=ISTUP_aux(i)
            do mu=1,5
               PUP(mu,i)=PUP_aux(mu,i)
            enddo
            do k=1,2
               ICOLUP(k,i)=ICOLUP_aux(k,i)
            enddo
            VTIMUP(i)=VTIMUP_aux(i)
            SPINUP(i)=SPINUP_aux(i)
            do k=1,2
               MOTHUP(k,i)=MOTHUP_aux(k,i)
            enddo
         enddo


      else

         NUP=maxNparticles-maxNmothers

         do i=1,NUP
            if(i.le.2)then
               index=i
            else
               index=i+maxNmothers
            endif

            IDUP(i)=IDUP_aux(index)
            ISTUP(i)=ISTUP_aux(index)
            do mu=1,5
               PUP(mu,i)=PUP_aux(mu,index)
            enddo
            do k=1,2
               ICOLUP(k,i)=ICOLUP_aux(k,index)
            enddo
            VTIMUP(i)=VTIMUP_aux(index)
            SPINUP(i)=SPINUP_aux(index)
            do k=1,2
               MOTHUP(k,i)=MOTHUP_aux(k,index)
            enddo
         enddo

      endif                     ! if(iwrite_mothers.eq.1)then

      return
      end










***********************************************************************
* SUBROUTINE add_LH_mother                                            *
*                                                                     *
* Purpuse: fill LH commons for the candidate mother (if consistent    *
*          with colour flow configuration chosen for the event at     *
*          hand)                                                      *
***********************************************************************

      subroutine add_LH_mother(ndaughters,iord_daughters_LH,
     &     rmass_mapping,i_mother,IDUP,ISTUP,MOTHUP,PUP,ICOLUP,
     &     VTIMUP,SPINUP,isconsistent)

      implicit real*8 (a-b,d-h,o-z)
c common.h is required to know the value of boson and top masses
      include 'common.h'

      integer maxNUP
      parameter (maxNUP=12)
      integer ndaughters,iord_daughters_LH(6),
     &     i,k,idaux,nf,nb,mu,imoth1,imoth2,
     &     i_mother,IDUP(maxNUP),ISTUP(maxNUP),MOTHUP(2,maxNUP),
     &     ICOLUP(2,maxNUP)
      real*8 rmass_mapping,PUP(5,maxNUP),VTIMUP(maxNUP),SPINUP(maxNUP)

c Determine the ID of the candidate mother (W+/-, Z, H or top) combining
c information from flavour of fermions and multimapping.
c
c   PDG  convention:
c   1=d, 2=u, 3=s, 4=c, 5=b, 6=t
c   11=e-, 12=v_e, 13=mu-, 14=v_mu, 15=tau-, 16=v_tau
c   all antiparticles have the same number but opposite sign
c   moreover: 21=gluon, 22=gamma, 23=Z0, 24=W+, 25=h

      idaux=0
      nf=0
      nb=0
      do i=1,ndaughters
         if(IDUP(iord_daughters_LH(i)).eq.0)then
            print*,'***add_LH_mother: ERROR 0'
            STOP
         endif
         if(abs(IDUP(iord_daughters_LH(i))).le.16)then ! fermion
            nf=nf+1
            idaux = idaux + IDUP(iord_daughters_LH(i))
            if(abs(IDUP(iord_daughters_LH(i))).eq.5)then ! b
               nb=nb+1
            endif
         endif
      enddo

      if(nf.gt.0)then
         if(idaux.eq.0)then     ! neutral boson
            if(rmass_mapping.eq.rmz)then ! Z
               IDUP(i_mother)=23
            elseif(rmass_mapping.eq.rmh .and. rmh.ge.0.d0)then ! H
               IDUP(i_mother)=25
            elseif(rmass_mapping.eq.rmhh .and. rmhh.ge.0.d0)then ! Heavy H
               IDUP(i_mother)=625
            else
               print*,'***add_LH_mother: ERROR 11'
               STOP
            endif
         elseif(idaux.eq.1 .and. rmass_mapping.eq.rmw)then ! W+
            IDUP(i_mother)=24
         elseif(idaux.eq.-1 .and. rmass_mapping.eq.rmw)then ! W-
            IDUP(i_mother)=-24
         elseif(nb.gt.0 .and. rmass_mapping.eq.rmt)then
            if(idaux.eq.6)then  ! t
               IDUP(i_mother)=6
            elseif(idaux.eq.-6)then ! tbar
               IDUP(i_mother)=-6
            else
               print*,'***add_LH_mother: ERROR 12'
               STOP
            endif
         else
            print*,'***add_LH_mother: ERROR 13'
            STOP
         endif
      else
         print*,'***add_LH_mother: ERROR 1'
      endif


c giuseppe 30/10/2007
c In case of a resonance candidate for Higgs decaying in two fermions
c (+ n gluons possibly), one must check for consistency that the
c daughter fermions are b quarks.
      isconsistent=1
      if( (IDUP(i_mother).eq.25 .or. IDUP(i_mother).eq.625)
     &                     .and. nf.eq.2)then
        if(nb.ne.2)then
          isconsistent=0
        endif
      endif
c end giuseppe 30/10/2007


c giuseppe 30/10/2007
      if(isconsistent.eq.1)then
c end giuseppe 30/10/2007

c Check consistency with the colour flow configuration chosen in the
c amplitude
        call set_LH_color(ndaughters,iord_daughters_LH,i_mother,
     &       IDUP,ICOLUP,isconsistent)

        if(isconsistent.eq.1)then

c Reconstruct momentum of the intermediate particle
          do mu=1,5
            PUP(mu,i_mother)=0.d0
          enddo
          do i=1,ndaughters
            PUP(1,i_mother)=PUP(1,i_mother)+PUP(1,iord_daughters_LH(i))
            PUP(2,i_mother)=PUP(2,i_mother)+PUP(2,iord_daughters_LH(i))
            PUP(3,i_mother)=PUP(3,i_mother)+PUP(3,iord_daughters_LH(i))
            PUP(4,i_mother)=PUP(4,i_mother)+PUP(4,iord_daughters_LH(i))
          enddo
          auxpx2 = PUP(1,i_mother)*PUP(1,i_mother)
          auxpy2 = PUP(2,i_mother)*PUP(2,i_mother)
          auxpz2 = PUP(3,i_mother)*PUP(3,i_mother)
          auxE2  = PUP(4,i_mother)*PUP(4,i_mother)
          PUP(5,i_mother) = sqrt(auxE2-auxpx2-auxpy2-auxpz2)

          MOTHUP(1,i_mother)=1
          MOTHUP(2,i_mother)=2

c Assign correct relationship with this mother
          do i=1,ndaughters

            imoth1=MOTHUP(1,iord_daughters_LH(i))
            imoth2=MOTHUP(2,iord_daughters_LH(i))

            if(imoth1.eq.1 .and. imoth2.eq.2)then
c External particle has not yet a mother. Assign to this.
              MOTHUP(1,iord_daughters_LH(i))=i_mother
              MOTHUP(2,iord_daughters_LH(i))=i_mother
            else
c External particle has already a mother. Then follow backwards the
c chain of decays to eventually update the relationship of intermediate
c particles involved.
              k=0
              do while(imoth1.ne.1 .and. imoth2.ne.2)
                k=imoth1
                if(imoth1.eq.0 .or. imoth2.eq.0 .or.
     &               imoth1.ne.imoth2)then
                  print*,'***add_LH_mother: ERROR 2'
                  STOP
                endif
                imoth1=MOTHUP(1,imoth1)
                imoth2=MOTHUP(2,imoth2)
              enddo
              if(k.lt.i_mother)then
                MOTHUP(1,k)=i_mother
                MOTHUP(2,k)=i_mother
              endif
            endif
          enddo


c Complete Les Houches information for LHEF
c Notice: ICOLUP has been already determined inside subroutine
c set_LH_color.
          ISTUP(i_mother)=2
          VTIMUP(i_mother)=0.d0
          SPINUP(i_mother)=9.d0

        else
          IDUP(i_mother)=0
        endif                   ! if(isconsistent.eq.1)then
c giuseppe 30/10/2007
      else
        IDUP(i_mother)=0
      endif                     ! if(isconsistent.eq.1)then
c end giuseppe 30/10/2007

      return
      end




***********************************************************************
* SUBROUTINE add_LH_mother2                                           *
*                                                                     *
* Purpuse: fill LH commons for the candidate mother (if consistent    *
*          with colour flow configuration chosen for the event at     *
*          hand)                                                      *
*          Grandmothers are given in input                            *
***********************************************************************

      subroutine add_LH_mother2(ndaughters,iord_daughters_LH,
     &     rmass_mapping,i_mother,igranmoth1,igranmoth2,
     &     IDUP,ISTUP,MOTHUP,PUP,ICOLUP,
     &     VTIMUP,SPINUP,isconsistent)

      implicit real*8 (a-b,d-h,o-z)
c common.h is required to know the value of boson and top masses
      include 'common.h'

      integer maxNUP
      parameter (maxNUP=12)
      integer ndaughters,iord_daughters_LH(6),
     &     i,k,idaux,nf,nb,mu,imoth1,imoth2,igranmoth1,igranmoth2,
     &     i_mother,IDUP(maxNUP),ISTUP(maxNUP),MOTHUP(2,maxNUP),
     &     ICOLUP(2,maxNUP)
      real*8 rmass_mapping,PUP(5,maxNUP),VTIMUP(maxNUP),SPINUP(maxNUP)

c Determine the ID of the candidate mother (W+/-, Z, H or top) combining
c information from flavour of fermions and multimapping.
c
c   PDG  convention:
c   1=d, 2=u, 3=s, 4=c, 5=b, 6=t
c   11=e-, 12=v_e, 13=mu-, 14=v_mu, 15=tau-, 16=v_tau
c   all antiparticles have the same number but opposite sign
c   moreover: 21=gluon, 22=gamma, 23=Z0, 24=W+, 25=h

      idaux=0
      nf=0
      nb=0
      do i=1,ndaughters
         if(IDUP(iord_daughters_LH(i)).eq.0)then
            print*,'***add_LH_mother: ERROR 0'
            STOP
         endif
         if(abs(IDUP(iord_daughters_LH(i))).le.16)then ! fermion
            nf=nf+1
            idaux = idaux + IDUP(iord_daughters_LH(i))
            if(abs(IDUP(iord_daughters_LH(i))).eq.5)then ! b
               nb=nb+1
            endif
         endif
      enddo

      if(nf.gt.0)then
         if(idaux.eq.0)then     ! neutral boson
            if(rmass_mapping.eq.rmz)then ! Z
               IDUP(i_mother)=23
            elseif(rmass_mapping.eq.rmh .and. rmh.ge.0.d0)then ! H
               IDUP(i_mother)=25
            elseif(rmass_mapping.eq.rmhh .and. rmhh.ge.0.d0)then ! Heavy H
               IDUP(i_mother)=625
            else
               print*,'***add_LH_mother: ERROR 11'
               STOP
            endif
         elseif(idaux.eq.1 .and. rmass_mapping.eq.rmw)then ! W+
            IDUP(i_mother)=24
         elseif(idaux.eq.-1 .and. rmass_mapping.eq.rmw)then ! W-
            IDUP(i_mother)=-24
         elseif(nb.gt.0 .and. rmass_mapping.eq.rmt)then
            if(idaux.eq.6)then  ! t
               IDUP(i_mother)=6
            elseif(idaux.eq.-6)then ! tbar
               IDUP(i_mother)=-6
            else
               print*,'***add_LH_mother: ERROR 12'
               STOP
            endif
         else
            print*,'***add_LH_mother: ERROR 13'
            STOP
         endif
      else
         print*,'***add_LH_mother: ERROR 1'
      endif


c giuseppe 30/10/2007
c In case of a resonance candidate for Higgs decaying in two fermions
c (+ n gluons possibly), one must check for consistency that the
c daughter fermions are b quarks.
      isconsistent=1
      if( (IDUP(i_mother).eq.25 .or. IDUP(i_mother).eq.625)
     &                     .and. nf.eq.2)then
        if(nb.ne.2)then
          isconsistent=0
        endif
      endif
c end giuseppe 30/10/2007


c giuseppe 30/10/2007
      if(isconsistent.eq.1)then
c end giuseppe 30/10/2007

c Check consistency with the colour flow configuration chosen in the
c amplitude
        call set_LH_color(ndaughters,iord_daughters_LH,i_mother,
     &       IDUP,ICOLUP,isconsistent)

c Reconstruct momentum of the intermediate particle
        do mu=1,5
          PUP(mu,i_mother)=0.d0
        enddo
        do i=1,ndaughters
          PUP(1,i_mother)=PUP(1,i_mother)+PUP(1,iord_daughters_LH(i))
          PUP(2,i_mother)=PUP(2,i_mother)+PUP(2,iord_daughters_LH(i))
          PUP(3,i_mother)=PUP(3,i_mother)+PUP(3,iord_daughters_LH(i))
          PUP(4,i_mother)=PUP(4,i_mother)+PUP(4,iord_daughters_LH(i))
        enddo
        auxpx2 = PUP(1,i_mother)*PUP(1,i_mother)
        auxpy2 = PUP(2,i_mother)*PUP(2,i_mother)
        auxpz2 = PUP(3,i_mother)*PUP(3,i_mother)
        auxE2  = PUP(4,i_mother)*PUP(4,i_mother)
        PUP(5,i_mother) = sqrt(auxE2-auxpx2-auxpy2-auxpz2)

        MOTHUP(1,i_mother)=igranmoth1
        MOTHUP(2,i_mother)=igranmoth2

c Assign correct relationship with this mother
        do i=1,ndaughters
          MOTHUP(1,iord_daughters_LH(i))=i_mother
          MOTHUP(2,iord_daughters_LH(i))=i_mother
        enddo

c Complete Les Houches information for LHEF
c Notice: ICOLUP has been already determined inside subroutine
c set_LH_color.
        ISTUP(i_mother)=2
        VTIMUP(i_mother)=0.d0
        SPINUP(i_mother)=9.d0

      else
        IDUP(i_mother)=0
      endif                     ! if(isconsistent.eq.1)then
c end giuseppe 30/10/2007

      return
      end







************************************************************************
* SUBROUTINE set_LH_color                                              *
*                                                                      *
* Purpose: Check colour consistency and assign colour flow information *
*          (ICOLUP) to a candidate mother                              *
************************************************************************

      subroutine set_LH_color(ndaughters,iord_daughters_LH,i_mother,
     &     IDUP,ICOLUP,isconsistent)

      implicit none

      integer maxNUP
      parameter (maxNUP=12)

      integer ndaughters,iord_daughters_LH(6),i_mother,
     &     IDUP(maxNUP),ICOLUP(2,maxNUP),isconsistent,i,j,ifound,
     &     nfound,istop,imatch


      isconsistent=1

      if(IDUP(i_mother).eq.23 .or. abs(IDUP(i_mother)).eq.24 .or.
     &     IDUP(i_mother).eq.25)then ! EW boson

c Check if the mother has a colour configuration consistent with the
c one chosen in the amplitude: the daughters must form a colour singlet.
c All nonzero LEFT and RIGHT colour flow lines (i.e. ICOLUP(1,.) and
c ICOLUP(2,.)) must be shared between daughter particles.
c This means that for each nonzero ICOLUP(1,i) there must be an
c ICOLUP(2,j) with the same value, where i,j indicate two particles in
c the set of daughters. The same yelds for 2<->1.

         nfound=0
         istop=0
         i=1
         do while(i.le.ndaughters .and. istop.eq.0)
            j=1
            if(ICOLUP(1,iord_daughters_LH(i)).ne.0)then
               imatch=0
               do while(imatch.eq.0 .and. j.le.ndaughters)
                  if(i.ne.j .and. ICOLUP(2,iord_daughters_LH(j)).eq.
     &                 ICOLUP(1,iord_daughters_LH(i)))then
                     imatch=1
                  endif
                  j=j+1
               enddo
               if(j.eq.(ndaughters+1) .and. imatch.eq.0)then
                  nfound=nfound+1
                  if(nfound.gt.0)then
                     istop=1
                  endif
               endif
            endif

            j=1
            if(ICOLUP(2,iord_daughters_LH(i)).ne.0)then
               imatch=0
               do while(imatch.eq.0 .and. j.le.ndaughters)
                  if(i.ne.j .and. ICOLUP(1,iord_daughters_LH(j)).eq.
     &                 ICOLUP(2,iord_daughters_LH(i)))then
                     imatch=1
                  endif
                  j=j+1
               enddo
               if(j.eq.(ndaughters+1) .and. imatch.eq.0)then
                  nfound=nfound+1
                  if(nfound.gt.0)then
                     istop=1
                  endif
               endif
            endif
            i=i+1
         enddo

         if(nfound.ne.0)then
c If there is one or more colour flow lines not shared between
c daughters, then consistency with the colour configuration chosen in
c the amplitude is lost.
            isconsistent=0
         endif

         if(isconsistent.eq.1)then
            ICOLUP(1,i_mother)=0
            ICOLUP(2,i_mother)=0
         endif


      elseif(IDUP(i_mother).eq.6)then ! top

c In case of 3 daughters only [i.e. t -> W b], the colour configuration
c (ICOLUP) of the top quark is identical to that of the b quark:
c
c                              \ /
c                               |
c                               |W
c                           t___|___b
c
c In the general case of n daughters [i.e. t-> W b + (n-3)gluons],
c top ICOLUP is reconstructed by finding the (only) daughter which
c does not share a nonzero LEFT colour flow line (i.e. ICOLUP(1,.))
c with other daughter particles. It is indicated as g* in the picture
c below:
c
c                g*  g  g     g    \ /   g     g
c                 |  |  |     |     |    |     |
c                 |  |  | ... |     |W   | ... |
c             t___|__|__|_____|_____|____|_____|___b
c

         ifound=0
         nfound=0
         istop=0
         i=1
         do while(i.le.ndaughters .and. istop.eq.0)
            j=1
            if(ICOLUP(1,iord_daughters_LH(i)).ne.0)then
               imatch=0
               do while(imatch.eq.0 .and. j.le.ndaughters)
                  if(i.ne.j .and. ICOLUP(2,iord_daughters_LH(j)).eq.
     &                 ICOLUP(1,iord_daughters_LH(i)))then
                     imatch=1
                  endif
                  j=j+1
               enddo
               if(j.eq.(ndaughters+1) .and. imatch.eq.0)then
c a candidate for top ICOLUP has been found.
                  nfound=nfound+1
                  ifound=i
                  if(nfound.gt.1)then
                     istop=1
                  endif
               endif
            endif
            i=i+1
         enddo

         if(nfound.ne.1)then
c If there is more than one candidate for top colour flow, i.e. more
c than one colour flow line not shared between daughters, then
c consistency with the colour configuration chosen in the amplitude
c is lost.
            isconsistent=0
         endif

         if(isconsistent.eq.1)then
            ICOLUP(1,i_mother)=ICOLUP(1,iord_daughters_LH(ifound))
            ICOLUP(2,i_mother)=ICOLUP(2,iord_daughters_LH(ifound))
         endif


      elseif(IDUP(i_mother).eq.-6)then ! antitop

c In the general case of n daughters [i.e. tbar-> W bbar + (n-3)gluons],
c top ICOLUP is reconstructed by finding the (only) daughter which
c does not share a nonzero RIGHT colour flow line (i.e. ICOLUP(2,.))
c with other daughter particles. It is indicated as g* in the picture
c below:
c
c                g*  g  g     g    \ /   g     g
c                 |  |  |     |     |    |     |
c                 |  |  | ... |     |W   | ... |
c          tbar___|__|__|_____|_____|____|_____|___bbar
c

         ifound=0
         nfound=0
         istop=0
         i=1
         do while(i.le.ndaughters .and. istop.eq.0)
            j=1
            if(ICOLUP(2,iord_daughters_LH(i)).ne.0)then
               imatch=0
               do while(imatch.eq.0 .and. j.le.ndaughters)
                  if(i.ne.j .and. ICOLUP(1,iord_daughters_LH(j)).eq.
     &                 ICOLUP(2,iord_daughters_LH(i)))then
                     imatch=1
                  endif
                  j=j+1
               enddo
               if(j.eq.(ndaughters+1) .and. imatch.eq.0)then
c a candidate for top ICOLUP has been found.
                  nfound=nfound+1
                  ifound=i
                  if(nfound.gt.1)then
                     istop=1
                  endif
               endif
            endif
            i=i+1
         enddo

         if(nfound.ne.1)then
c If there is more than one candidate for antitop colour flow, i.e. more
c than one colour flow line not shared between daughters, then
c consistency with the colour configuration chosen in the amplitude
c is lost.
            isconsistent=0
         endif

         if(isconsistent.eq.1)then
            ICOLUP(1,i_mother)=ICOLUP(1,iord_daughters_LH(ifound))
            ICOLUP(2,i_mother)=ICOLUP(2,iord_daughters_LH(ifound))
         endif

      endif                     ! if(IDUP(i_mother).eq. ...)

      return
      end








************************************************************************
* SUBROUTINE shift_LH_stack                                            *
*                                                                      *
* Purpose: shift of one unit, from index_shift on, a stack of LHEF-like*
* arrays with fixed dimension.                                         *
************************************************************************

      subroutine shift_LH_stack(index_shift,IDUP,ISTUP,MOTHUP,PUP,
     &     ICOLUP,VTIMUP,SPINUP)

      implicit none

      integer maxNUP
      parameter (maxNUP=12)

      integer index_shift,NUP,IDUP(maxNUP),ISTUP(maxNUP),
     &     MOTHUP(2,maxNUP),ICOLUP(2,maxNUP),i,mu,k

      real*8 PUP(5,maxNUP),VTIMUP(maxNUP),SPINUP(maxNUP)


      do i=index_shift+1,maxNUP
         IDUP(i-1)=IDUP(i)
         ISTUP(i-1)=ISTUP(i)
         do mu=1,5
            PUP(mu,i-1)=PUP(mu,i)
         enddo
         do k=1,2
            ICOLUP(k,i-1)=ICOLUP(k,i)
         enddo
         VTIMUP(i-1)=VTIMUP(i)
         SPINUP(i-1)=SPINUP(i)
         do k=1,2
            MOTHUP(k,i-1)=MOTHUP(k,i)
         enddo
      enddo

c At this point the shift has been accomplished.
c As MOTHUP information is connected with the ordering in LHEF particle
c list, after a shift one must update the values of MOTHUP's when
c required.

      do i=index_shift,maxNUP
         do k=1,2
            if(MOTHUP(k,i).GE.index_shift)then
               MOTHUP(k,i)=MOTHUP(k,i)-1
            endif
         enddo
      enddo

      return
      end









