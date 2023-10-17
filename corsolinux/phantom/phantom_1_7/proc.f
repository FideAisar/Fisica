c
c Last update: Jun 25, 2007
c

C For 8,9/6,7/19,20,21 add check on identical W+/W-/Z

C INPUT (integer vector)
* iproc (e.g.  1 1 -> 2 2 3 -4 13 -14  for   d d -> u u s c~ mu- vm~)


C OUTPUT(igroup,n_kphs,ialfa,idp,idp_inout,iorder,{input_phasespace})
* igroup=group number
* n_kphs=total number of phase spaces
* ialfa(i)=0/1 (off/on phase space)
* idp(i)= particle identity according to the amplitude order
* idp_inout(j,i)=-1/+1  (incoming/outgoing particles)
* iorder(j,i) gives the momenta order as passed to the amplitude 
* i.e. iorder(j_phsp,iphs_ind)=amplitude_index 
* {input_phasespace} = set of phase_space input 
* e.g. rm1(i),gam1(i) etc.

C   PDG  convention:
C   1=d, 2=u, 3=s, 4=c, 5=b, 6=t
C   11=e-, 12=v_e, 13=mu-, 14=v_mu, 15=tau-, 16=v_tau
C   all antiparticles have the same number but opposite sign
C   moreover: 21=gluon, 22=gamma, 23=Z0, 24=W+, 25=h


* Phase space list and numbering:

* t(00) W+ W-                                     ialfa(1)  1_1_4
* t(0W) Z W-  and t(0W) Z W+                      ialfa(2)  1_1_4
* t(W0) Z W-  and t(W0) Z W+                      ialfa(3)  1_1_4
* t(WW) W- W- and t(WW) W+ W+                     ialfa(4)  1_1_4
* t(WW) W+ W-                                     ialfa(5)  1_1_4
* t(00) Z Z                                       ialfa(6)  1_1_4
* t(WW) Z Z                                       ialfa(7)  1_1_4
* t(g0) W+ W-                                     ialfa(8)  1_1_4
* t(0g) W+ W-                                     ialfa(9)  1_1_4
* t(gW) Z W-   and  t(gW) Z W+                    ialfa(10) 1_1_4
* t(Wg) Z W-   and  t(Wg) Z W+                    ialfa(11) 1_1_4
* t(g0) Z Z                                       ialfa(12) 1_1_4
* t(0g) Z Z                                       ialfa(13) 1_1_4
* W- -> [W- W+] W-                                ialfa(14) 2_4 
* W- -> W- [W+ W-]                                ialfa(15) 2_4 
* W+ -> [W+ W-] W+                                ialfa(16) 2_4 
* W+ -> W+ [W- W+]                                ialfa(17) 2_4 
* Z  -> Z  W+ W-   and  gg -> Z  W+ W-            ialfa(18) 2_4 
* W- -> Z  Z  W- and W+ -> Z  Z  W+               ialfa(19) 2_4 
* Z  -> [Z Z] Z   and gg -> [Z Z] Z               ialfa(20) 2_4 
* Z  -> [Z Z Z]                                   ialfa(21) 2_4 
* Z  -> Z [Z Z]                                   ialfa(22) 2_4 
* qq -> {gq} {gqW-}   and  qq -> {gq} {gqW+}      ialfa(23) 2_4
* qq -> {gq} {gqZ}                                ialfa(24) 2_4
* t(W0) t                                         ialfa(25) 1_1_31
* t(0W) t                                         ialfa(26) 1_1_31
* t(0W) {ggW-}   and   t(0W) {ggW+}               ialfa(27) 1_1_31 
* t(W0) {ggW-}   and   t(W0) {ggW+}               ialfa(28) 1_1_31
* t(WW) {ggZ}                                     ialfa(29) 1_1_31
* t(gW) t                                         ialfa(30) 1_1_31
* t(Wg) t                                         ialfa(31) 1_1_31
* t(00) {ggZ}                                     ialfa(32) 1_1_31
* t(W ) Z  t                                      ialfa(33) 1_2_3
* t(0 ) W- t                                      ialfa(34) 1_2_3
* t(0 ) W+ t                                      ialfa(35) 1_2_3
* t(g ) W- t    and  t(g ) W+ t                   ialfa(36) 1_2_3
* g t(W ) Z {gW-}   and   g t(W ) Z {gW+}         ialfa(37) 1_2_3
* g t(W ) {gZ} W-   and   g t(W ) {gZ} W+         ialfa(38) 1_2_3
* g t(0 ) W+ {gW-}                                ialfa(39) 1_2_3
* g t(0 ) {gW+} W-                                ialfa(40) 1_2_3
* g t(0 ) {gZ} Z                                  ialfa(41) 1_2_3
* g t(0 ) Z {gZ}                                  ialfa(42) 1_2_3
* t(0 ) W {ggq} -> t(0 ) W g{gq}                  ialfa(43) 1_2_3
* t(0 ) W {ggq} -> t(0 ) W q{gg}                  ialfa(44) 1_2_3
* t(W ) W  {ggq} -> t(W ) W  g{gq}                ialfa(45) 1_2_3
* t(W ) W  {ggq} -> t(W ) W  q{gg}                ialfa(46) 1_2_3
* t(W+ ) Z {ggq} -> t(W+ ) Z g{gq}                ialfa(47) 1_2_3
* t(W+ ) Z {ggq} -> t(W+ ) Z q{gg}                ialfa(48) 1_2_3
* t(W- ) Z {ggq} -> t(W- ) Z g{gq}                ialfa(49) 1_2_3
* t(W- ) Z {ggq} -> t(W- ) Z q{gg}                ialfa(50) 1_2_3
* t(0[1] ) Z {ggq} -> t(0[1] ) Z g{gq}            ialfa(51) 1_2_3
* t(0[1] ) Z {ggq} -> t(0[1] ) Z q{gg}            ialfa(52) 1_2_3
* t(0[2] ) Z {ggq} -> t(0[2] ) Z g{gq}            ialfa(53) 1_2_3
* t(0[2] ) Z {ggq} -> t(0[2] ) Z q{gg}            ialfa(54) 1_2_3
* W- -> t  Z  and  W+ -> t  Z                     ialfa(55) 3_3 
* Z  -> t  t   and  gg  -> t  t                   ialfa(56) 3_3 
* W- -> {gZ} {gW-}  and  W+ -> {gZ} {gW+}         ialfa(57) 3_3
* Z  -> {gW+} {gW-}                               ialfa(58) 3_3
* Z  -> {gZ} {gZ}                                 ialfa(59) 3_3
* gq -> {gW+} t    and  gq -> {gW-} t             ialfa(60) 3_3
* W- -> {ggZ} W-     and  W+ -> {ggZ} W+          ialfa(61) 2_4to31
* W- -> Z {ggW-}    and  W+ -> Z {ggW+}           ialfa(62) 2_4to31
* Z  -> {ggW+} W-                                 ialfa(63) 2_4to31
* Z  -> W+ {ggW-}                                 ialfa(64) 2_4to31
* gq -> W+ {gt} -> W+ g W-
*            and   gq -> W- {gt} -> W- g W+       ialfa(65) 2_4to31
* gq -> W+ g t -> W+ {gW-}
*            and   gq -> W- g t -> W- {gW+}       ialfa(66) 2_4to31
* Z  -> {ggZ} Z                                   ialfa(67) 2_4to31
* Z  -> Z {ggZ}                                   ialfa(68) 2_4to31
* W- -> {ggtbar} -> gg W-   
*            and  W+ -> {ggt} -> gg W+            ialfa(69) 1_5to1_4to31
* W- -> g {gtbar} -> g {gW-}                           
*            and  W+ -> g {gt} -> g {gW+}         ialfa(70) 1_5to1_4to31
* W- -> gg tbar -> {ggW-}   
*            and  W+ -> gg t -> {ggW+}            ialfa(71) 1_5to1_4to31
* t(W ) {ggtbar} -> t(W ) gg W-                              
*            and  t(W ) {ggt} -> t(W ) gg W+      ialfa(72) 1_5to1_4to31
* t(W ) g {gtbar} -> t(W ) g {gW-}  
*            and  t(W ) g {gt} -> t(W ) g {gW+}   ialfa(73) 1_5to1_4to31
* t(W ) gg tbar -> t(W ) {ggW-} 
*            and  t(W ) gg t -> t(W ) {ggW+}      ialfa(74) 1_5to1_4to31
* g t(W ) {gtbar} -> g t(W ) g W-  
*            and   g t(W ) {gt} -> g t(W ) g W+   ialfa(75) 1_5to1_4to31
* g t(W ) g tbar -> g t(W ) {gW-}  
*            and   g t(W ) g t -> g t(W ) {gW+}   ialfa(76) 1_5to1_4to31




* idout(i) i=1,8 particle/antiparticle id all assumed outgoing
* idpart(i) i=1,4  particle id
* idpartinout(i) i=1,4  particle initial(-1) or final(1) status
* idantipart(i) i=1,4  antiparticle id
* idantipartinout(i) i=1,4  antiparticle initial(-1) or final(1) status
* ipartinitial(i),i=1,2 tells whether the initial particles are 
*   particles (1) or antiparticles (-1)
* nwm(j),nwp(j),nz0(j) number of W-/W+/Z in the j_th match 
* idwm(i,j),idwp(i,j),idz0(i,j): identity of the particle in i-th pair 
*   in the j_th match 
* iwminout(i,j),iwpinout(i,j),iz0inout(i,j): outgoing(2), t-type(0) or 
*   incoming(-2) status of i-th pair in the j_th match 
* itype(i,j)=-1/+1/0 tells whether the i-th pair in the j_th match is a 
*   W-/W+/Z
* itypeinitial(i)=-1/+1/0 i=1,2 gives type of boson which contain 
*   initial particles
*
* igr(i)=0/1 indicates whether the i-th type of phase space is present, 
*   where:
*    1 --> 4W   2 --> 2Z2W-1  3 --> 2Z2W-2   4 --> 4Z   5 --> 2g1Z2W  
*    6 --> 2g3Z
*  NB: for 2Z2W there can be two different possibilities if the state
*      can be reconstructed as 4W AND as 2Z2W AND 4Z
* igo(8,6) gives the ordering in which the momenta are passed to the 
*   amplitudes
*
* Momenta are ordered as follows:
* 4W            W+W-W+W-
* 2Z2W          Z Z W+W-  (beware ud~du~cs~sc~: two possibilities exist)
* 4Z            Z Z Z Z
* 2gZ2W         ggZ W+W-
* 2g3Z          ggZ Z Z


      SUBROUTINE  proc

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT DOUBLE COMPLEX (c)
      
      INCLUDE 'common.h'
      INCLUDE 'common_subproc.h'
      DIMENSION ifact(0:4),index(4,24),invindex(4,24)
      DATA ifact/1,1,2,6,24/
      DATA index/
     &     1,2,3,4,               !even
     &     2,1,3,4,
     &     2,3,1,4,               !even
     &     1,3,2,4,
     &     3,1,2,4,               !even
     &     3,2,1,4,
     
     &     3,2,4,1,               !even
     &     3,1,4,2,               
     &     1,3,4,2,               !even
     &     2,3,4,1,               
     &     2,1,4,3,               !even
     &     1,2,4,3,               

     &     1,4,2,3,               !even
     &     2,4,1,3,
     &     2,4,3,1,               !even
     &     1,4,3,2,
     &     3,4,1,2,               !even
     &     3,4,2,1,

     &     4,3,2,1,               !even
     &     4,3,1,2,               
     &     4,1,3,2,               !even
     &     4,2,3,1,               
     &     4,2,1,3,               !even
     &     4,1,2,3               
     &         /
      DATA invindex/
     &     1,2,3,4,               !even
     &     2,1,3,4,
     &     3,1,2,4,               !even
     &     1,3,2,4,
     &     2,3,1,4,               !even
     &     3,2,1,4,

     &     4,2,1,3,
     &     2,4,1,3,
     &     1,4,2,3,
     &     4,1,2,3,
     &     2,1,4,3,
     &     1,2,4,3,
     &     1,3,4,2,
     &     3,1,4,2,
     &     4,1,3,2,
     &     1,4,3,2,
     &     3,4,1,2,
     &     4,3,1,2,
     &     4,3,2,1,
     &     3,4,2,1,
     &     2,4,3,1,
     &     4,2,3,1,
     &     3,2,4,1,
     &     2,3,4,1
     &         /
      DIMENSION inout(8,mxnphs)
      DIMENSION idhep(8),idout(8)
      DIMENSION ipartinout(4),iantipartinout(4),
     &          idpart(4), idantipart(4), !igluon(4),
     &          idwplus(2),iwplusinout(2),iposwplus(2),
     &          idwminus(2),iwminusinout(2),iposwminus(2),
     &          idzzero(4),izzeroinout(4),iposzzero(4),
     &          igluoninout(4),
     &          nwp(8),nwm(8),nz0(8),iperm(8),
     &          idwp(2,8),iwpinout(2,8),iposwp(2,8),
     &          idwm(2,8),iwminout(2,8),iposwm(2,8),
     &          idz0(4,8),iz0inout(4,8),iposz0(4,8),
     &          idz0_aux(2),
     &          ipartinitial(2),itype(4,8),
     &          nwpf(8),nwmf(8),nz0f(8),nfinal(8),
     &          iz0_b(4,8),
     &          idpartlocal(4),idantipartlocal(4),
     &          ipartiolocal(4),iantipartiolocal(4),
     &          idp_iolocal(8),idplocal(8),
     &          iswap(2),itypeinitial(2),
c final state bosons which decay into leptons
     &          nwpflep(8),nwmflep(8),nz0flep(8)

* Information for procextraini.f
      COMMON/info_proc/
c number of leptonic final state bosons for each channel
     &     num_wpflep(mxnphs),num_wmflep(mxnphs),num_z0flep(mxnphs)


* initializing output variables:
      do j=1,mxnphs
        ialfa(j)=0
      enddo
      DO i=1,mxnphs
        DO j=1,8 
          iorder(j,i)=0
        ENDDO
      ENDDO
      DO j=1,8 
        idp_inout(j)=1
      ENDDO

* particle identity
      DO i=1,8
        idhep(i)=iproc(i)
      ENDDO

      DO i=1,6
        igr(i)=0
      ENDDO

* 'outgoing process' identity number: the sign of the incoming fermions
* is reversed in order to consider all outgoing particles. 
      DO i=1,2
        IF(idhep(i).NE.21)THEN
          idout(i)=-idhep(i)
        ELSE
          idout(i)=idhep(i)
        ENDIF
      ENDDO
      DO i=3,8
        idout(i)=idhep(i)
      ENDDO

      n_kphs=0
      

* Gluons particles and antiparticles

      ipartinitial(1)=0
      ipartinitial(2)=0

      ngluon=0
      ngf=0
      npart=0
      nantipart=0
      DO i=1,8
        IF(idout(i).eq.21)THEN
           ngluon=ngluon+1
c           igluon(ngluon)=idout(i)
            IF(i.LT.3)THEN
              igluoninout(ngluon)=-1
            ELSE
              igluoninout(ngluon)=1
              ngf=ngf+1
            ENDIF
        ELSE
          IF(idout(i).LT.0)THEN
            nantipart=nantipart+1
            idantipart(nantipart)=idout(i)
            IF(i.LT.3)THEN
              iantipartinout(nantipart)=-1
              ipartinitial(i)=-1
            ELSE
              iantipartinout(nantipart)=1
            ENDIF
          ELSE
            npart=npart+1
            idpart(npart)=idout(i)
            IF(i.LT.3)THEN
              ipartinout(npart)=-1
              ipartinitial(i)=1
            ELSE
              ipartinout(npart)=1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      ngin=ngluon-ngf
      
      IF(npart.NE.nantipart)THEN
        WRITE(*,*) 'INPUT ERROR: Fermion number is not conserved'
        STOP
c        RETURN
      ENDIF

* Add gluons
      DO i=1,ngluon/2
        idpart(npart+i)=21       
        ipartinout(npart+i)=igluoninout(2*i-1)
        idantipart(npart+i)=21
        iantipartinout(npart+i)=igluoninout(2*i)
      ENDDO
      
* Check f-fbar pairs 
      nmatch=0
      DO ipermute=1,ifact(npart)
        imatch=1
        j=1
        nwplus=0
        nwminus=0
        nzzero=0
        DO WHILE(imatch.EQ.1 .AND. j.LE.npart)
          idtemp=idpart(index(j,ipermute))
          idiff=idtemp+idantipart(j)
          IF(abs(idiff).GT.1)THEN
            imatch=0
          ELSE
            iinouttemp=ipartinout(index(j,ipermute))
     &                    +iantipartinout(j)
            ifdown=mod(idtemp,2)  ! 0=up 1=down
            IF(idiff.EQ.1 .AND. ifdown.EQ.0)THEN
              jtemp=1
              ifound=0

              DO WHILE(jtemp.LE.nwplus .AND. ifound.EQ.0)
                IF(iinouttemp.LT.iwplusinout(jtemp))THEN
                  ifound=1
                ELSEIF(iinouttemp.GT.iwplusinout(jtemp))THEN
                  jtemp=jtemp+1
                ELSE
                  IF(idtemp.LE.idwplus(jtemp))THEN
                    ifound=1
                  ELSE
                    jtemp=jtemp+1
                  ENDIF
                ENDIF
              ENDDO

*  Shift rest of the stack
              DO k=nwplus,jtemp,-1
                idwplus(k+1)=idwplus(k)
                iwplusinout(k+1)=iwplusinout(k)
                iposwplus(k+1)=iposwplus(k)
              ENDDO
              idwplus(jtemp)=idtemp
              iwplusinout(jtemp)=iinouttemp
              iposwplus(jtemp)=j
              nwplus=nwplus+1
            ELSEIF(idiff.EQ.-1 .AND. ifdown.EQ.1)THEN
              jtemp=1
              ifound=0

              DO WHILE(jtemp.LE.nwminus .AND. ifound.EQ.0)
                IF(iinouttemp.LT.iwminusinout(jtemp))THEN
                  ifound=1
                ELSEIF(iinouttemp.GT.iwminusinout(jtemp))THEN
                  jtemp=jtemp+1
                ELSE
                  IF(idtemp.LE.idwminus(jtemp))THEN
                    ifound=1
                  ELSE
                    jtemp=jtemp+1
                  ENDIF
                ENDIF
              ENDDO

*  Shift rest of the stack
              DO k=nwminus,jtemp,-1
                idwminus(k+1)=idwminus(k)
                iwminusinout(k+1)=iwminusinout(k)
                iposwminus(k+1)=iposwminus(k)
              ENDDO
              idwminus(jtemp)=idtemp
              iwminusinout(jtemp)=iinouttemp
              iposwminus(jtemp)=j
              nwminus=nwminus+1
            ELSEIF(idiff.EQ.0)THEN
              jtemp=1
              ifound=0

              DO WHILE(jtemp.LE.nzzero .AND. ifound.EQ.0)
                IF(iinouttemp.LT.izzeroinout(jtemp))THEN
                  ifound=1
                ELSEIF(iinouttemp.GT.izzeroinout(jtemp))THEN
                  jtemp=jtemp+1
                ELSE
                  IF(idtemp.LE.idzzero(jtemp))THEN
                    ifound=1
                  ELSE
                    jtemp=jtemp+1
                  ENDIF
                ENDIF
              ENDDO

*  Shift rest of the stack
              DO k=nzzero,jtemp,-1
                idzzero(k+1)=idzzero(k)
                izzeroinout(k+1)=izzeroinout(k)
                iposzzero(k+1)=iposzzero(k)
              ENDDO
              idzzero(jtemp)=idtemp
              izzeroinout(jtemp)=iinouttemp
              iposzzero(jtemp)=j
              nzzero=nzzero+1
            ELSE  ! nomatch
              imatch=0
            ENDIF
          ENDIF
          j=j+1
        ENDDO          
            
        IF(imatch.EQ.1)THEN
          IF(nwplus.NE.nwminus)THEN
            WRITE(*,*) 'INPUT ERROR: Charge is not conserved'
            STOP
c            RETURN
          ENDIF
* Check if already found
          inew=1
          DO i=1,nmatch
* Check if same number of all types of EW vector bosons
              IF(nwplus.NE.nwp(i) .OR. nwminus.NE.nwm(i)
     &           .OR. nzzero.NE.nz0(i))THEN
*  Let inew=1
              ELSE
* Check all EW vector bosons
                ifound=1
                j=1
                DO WHILE( ifound.EQ.1 .AND. j.LE.nwplus)
                  IF(idwplus(j).NE.idwp(j,i) .OR.
     &                iwplusinout(j).NE.iwpinout(j,i)) ifound=0
                  j=j+1
                ENDDO 
                j=1
                DO WHILE( ifound.EQ.1 .AND. j.LE.nwminus)
                  IF(idwminus(j).NE.idwm(j,i) .OR.
     &                iwminusinout(j).NE.iwminout(j,i)) ifound=0
                  j=j+1
                ENDDO 
                j=1
                DO WHILE( ifound.EQ.1 .AND. j.LE.nzzero)
                  IF(idzzero(j).NE.idz0(j,i) .OR.
     &                izzeroinout(j).NE.iz0inout(j,i)) ifound=0
                  j=j+1
                ENDDO
                IF(ifound.EQ.1) inew=0
            ENDIF
          ENDDO  ! on i=1,nmatch
* Record new match. Matches are ordered by decreasing number of Z's.
* Therefore 4Z comes before 2Z2W which comes before 4W
* ggZZZ comes before ggZWW
          IF(inew.EQ.1)THEN
* nmatch now is the number matches already in the stack   
             ifound=0
             jtemp=1
             DO WHILE(jtemp.LE.nmatch .AND. ifound.EQ.0)
               IF(nzzero.GE.nz0(jtemp))THEN
                 ifound=1
               ELSE
                 jtemp=jtemp+1
               ENDIF
             ENDDO
*  Shift rest of the stack
             DO k=nmatch,jtemp,-1
               nwp(k+1)=nwp(k)
               nwm(k+1)=nwm(k)
               nz0(k+1)=nz0(k)
               iperm(k+1)=iperm(k)
               DO ii=1,nwp(k)
                 idwp(ii,k+1)=idwp(ii,k)
                 iwpinout(ii,k+1)=iwpinout(ii,k)
                 iposwp(ii,k+1)=iposwp(ii,k)
c giuseppe 03/01/2007
c                 itype(iposwplus(ii),k+1)=itype(iposwplus(ii),k)
                 itype(iposwp(ii,k+1),k+1)=itype(iposwp(ii,k),k)
c end giuseppe 03/01/2006
               ENDDO
               nwpf(k+1)=nwpf(k)
               DO ii=1,nwm(k)
                 idwm(ii,k+1)=idwm(ii,k)
                 iwminout(ii,k+1)=iwminout(ii,k)
                 iposwm(ii,k+1)=iposwm(ii,k)
c giuseppe 03/01/2007
c                 itype(iposwminus(ii),k+1)=itype(iposwminus(ii),k)
                 itype(iposwm(ii,k+1),k+1)=itype(iposwm(ii,k),k)
c end giuseppe 03/01/2006
               ENDDO
               nwmf(k+1)=nwmf(k)
               DO ii=1,nz0(k)
                 idz0(ii,k+1)=idz0(ii,k)
                 iz0inout(ii,k+1)=iz0inout(ii,k)
                 iposz0(ii,k+1)=iposz0(ii,k)
c giuseppe 03/01/2007
c                 itype(iposzzero(ii),k+1)=itype(iposzzero(ii),k)
                 itype(iposz0(ii,k+1),k+1)=itype(iposz0(ii,k),k)
c end giuseppe 03/01/2006
                 iz0_b(ii,k+1)=iz0_b(ii,k)
               ENDDO
               nz0f(k+1)=nz0f(k)
               nfinal(k+1)=nfinal(k)
             ENDDO
* Insert new entry             
             nwp(jtemp)=nwplus
             nwm(jtemp)=nwminus
             nz0(jtemp)=nzzero
             iperm(jtemp)=ipermute
             nf=0
             nfaux=0
             DO ii=1,nwplus
               idwp(ii,jtemp)=idwplus(ii)
               iwpinout(ii,jtemp)=iwplusinout(ii)
               iposwp(ii,jtemp)=iposwplus(ii)
               itype(iposwplus(ii),jtemp)=1
               IF(iwplusinout(ii).EQ.2)THEN
                 nf=nf+1
                 nfaux=nfaux+1
               ENDIF
             ENDDO
             nwpf(jtemp)=nfaux
             nfaux=0
             DO ii=1,nwminus
               idwm(ii,jtemp)=idwminus(ii)
               iwminout(ii,jtemp)=iwminusinout(ii)
               iposwm(ii,jtemp)=iposwminus(ii)
               itype(iposwminus(ii),jtemp)=-1
               IF(iwminusinout(ii).EQ.2)THEN
                 nf=nf+1
                 nfaux=nfaux+1
               ENDIF
             ENDDO
             nwmf(jtemp)=nfaux
             nfaux=0
             DO ii=1,nzzero
               idz0(ii,jtemp)=idzzero(ii)
               iz0inout(ii,jtemp)=izzeroinout(ii)
               iposz0(ii,jtemp)=iposzzero(ii)
               itype(iposzzero(ii),jtemp)=0
               IF(idzzero(ii).EQ.5)THEN
                 iz0_b(ii,jtemp)=1
               ELSE
                 iz0_b(ii,jtemp)=0
               ENDIF
               IF(izzeroinout(ii).EQ.2)THEN
                 nf=nf+1
                 nfaux=nfaux+1
               ENDIF
             ENDDO
             nz0f(jtemp)=nfaux
             nfinal(jtemp)=nf

             nmatch=nmatch+1
          ENDIF  !  on  (inew.EQ.1)   

        ENDIF    ! on (imatch.EQ.1)
      ENDDO      ! i=1,ifact(npart)
 
      IF(nmatch.EQ.0)THEN
        WRITE(*,*) 
     &       'INPUT ERROR: The particle content does not correspond'
        WRITE(*,*) 'to a valid SM process'
        STOP
c        RETURN
      ENDIF
     
* Count final state b's and antib's
      nb=0
      nantib=0
      DO i=3,8
        IF(idout(i).EQ.5) nb=nb+1
        IF(idout(i).EQ.-5) nantib=nantib+1
      ENDDO

* Count final state bosons which decay into LEPTONS
      DO i=1,nmatch
        nwpflep(i)=0
        nwmflep(i)=0
        nz0flep(i)=0
      ENDDO

      DO i=1,nmatch

        DO k=1,nwp(i)
          IF(iwpinout(k,i).EQ.2 .AND. 
     &         (idwp(k,i).GE.11 .AND. idwp(k,i).LE.16) )THEN
            nwpflep(i)=nwpflep(i)+1 ! final state leptonic W+'s
          ENDIF
        ENDDO

        DO k=1,nwm(i)
          IF(iwminout(k,i).EQ.2 .AND. 
     &         (idwm(k,i).GE.11 .AND. idwm(k,i).LE.16) )THEN
            nwmflep(i)=nwmflep(i)+1 ! final state leptonic W-'s
          ENDIF
        ENDDO

        DO k=1,nz0(i)
          IF(iz0inout(k,i).EQ.2 .AND.
     &         (idz0(k,i).GE.11 .AND. idz0(k,i).LE.16) )THEN
            nz0flep(i)=nz0flep(i)+1 ! final state leptonic Z's
          ENDIF
        ENDDO

      ENDDO

* Information for procextraini.f: number of final state leptonic bosons
* for each phase space channel (to be assigned if required)
      DO k=1,mxnphs
        num_wpflep(k)=0
        num_wmflep(k)=0
        num_z0flep(k)=0
      ENDDO

c The following flag refers to phase space channel qq -> {gq} {gqZ} 
c only (see below).
      ifound_phsp=0


c      PRINT*,'ng:',ngluon,' ngf:',ngf
c      DO i=1,nmatch
c        PRINT*,i,'nwpf:',nwpf(i),'nwmf:',nwmf(i),'nz0f:',nz0f(i),
c     &            'nfinal:',nfinal(i)
c       PRINT*,'idz0',(idz0(j,i),j=1,nz0(i))
c      ENDDO



       
C  No gluons
      IF(ngluon.EQ.0)THEN
        DO i=1,nmatch
c TEST
c        print*,(idpart(index(j,iperm(i))),idantipart(j),j=1,4)
c        print*,(ipartinout(index(j,iperm(i))),iantipartinout(j),j=1,4)
c        PRINT*,i,'nwpf:',nwpf(i),'nwmf:',nwmf(i),'nz0f:',nz0f(i),
c     &            'nfinal:',nfinal(i)
c        PRINT*,'iz0_b,iz0inout:',(iz0_b(j,i),iz0inout(j,i),j=1,nz0(i))
c ENDTEST
          IF(nfinal(i).EQ.3)THEN
c 2_4 or 3_3
            IF(nwpf(i).EQ.2)THEN
*            WRITE(*,*) i,' W+ -> (W+ W-) W+ and   W+ (W- W+)  : 16,17'
              local_index=16
              ialfa(local_index)=1
              n_kphs=n_kphs+1
              
              IF(ipartinitial(1).EQ.1)THEN
                iorder(1,local_index)=7    ! W+ in --> W- 
                iorder(2,local_index)=8
              ELSE
                iorder(1,local_index)=8
                iorder(2,local_index)=7
              ENDIF
              iorder(3,local_index)=1    ! W+
              iorder(4,local_index)=2
              iorder(5,local_index)=3    ! W-
              iorder(6,local_index)=4  
              iorder(7,local_index)=5    ! W+
              iorder(8,local_index)=6
                
              IF(ipartinitial(1).EQ.1)THEN
                idplocal(iorder(1,local_index))=idwm(1,i)
                idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
              ELSE
                idplocal(iorder(1,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(2,local_index))=idwm(1,i)
              ENDIF
              idplocal(iorder(5,local_index))=idwm(2,i)
              idplocal(iorder(6,local_index))=-(idwm(2,i)+1)
              idplocal(iorder(3,local_index))=idwp(1,i)
              idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
              idplocal(iorder(7,local_index))=idwp(2,i)
              idplocal(iorder(8,local_index))=-(idwp(2,i)-1)

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(1)=1
                DO j=1,8
                  igo(j,1)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(1).EQ.0)THEN
                  igr(1)=1
                  DO j=1,4
                    igo(2*j-1,1)=2*index(j,ipp)-1
                    igo(2*j,1)=2*index(j,ipa)
                  ENDDO
                ENDIF 
              ENDIF
                
              IF(idwp(1,i).NE.idwp(2,i))THEN
                local_index=17
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=iorder(1,local_index-1)
                iorder(2,local_index)=iorder(2,local_index-1)
                iorder(3,local_index)=iorder(7,local_index-1)
                iorder(4,local_index)=iorder(8,local_index-1)
                iorder(5,local_index)=iorder(5,local_index-1)
                iorder(6,local_index)=iorder(6,local_index-1)
                iorder(7,local_index)=iorder(3,local_index-1)
                iorder(8,local_index)=iorder(4,local_index-1)
              ENDIF
            ELSEIF(nwmf(i).EQ.2)THEN
*            WRITE(*,*) i,' W- -> (W- W+) W- and   W- (W+ W-)  : 14,15'
              local_index=14
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              IF(ipartinitial(1).EQ.1)THEN
                iorder(1,local_index)=1    ! W- in --> W+ 
                iorder(2,local_index)=2
              ELSE
                iorder(1,local_index)=2 
                iorder(2,local_index)=1
              ENDIF
              iorder(3,local_index)=3    ! W-
              iorder(4,local_index)=4  
              iorder(5,local_index)=5    ! W+
              iorder(6,local_index)=6
              iorder(7,local_index)=7    ! W-
              iorder(8,local_index)=8
                
              IF(ipartinitial(1).EQ.1)THEN
                idplocal(iorder(1,local_index))=idwp(1,i)
                idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
              ELSE
                idplocal(iorder(1,local_index))=-(idwp(1,i)-1)
                idplocal(iorder(2,local_index))=idwp(1,i)
              ENDIF
              idplocal(iorder(5,local_index))=idwp(2,i)
              idplocal(iorder(6,local_index))=-(idwp(2,i)-1)
              idplocal(iorder(3,local_index))=idwm(1,i)
              idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
              idplocal(iorder(7,local_index))=idwm(2,i)
              idplocal(iorder(8,local_index))=-(idwm(2,i)+1)

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(1)=1
                DO j=1,8
                  igo(j,1)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(1).EQ.0)THEN
                  igr(1)=1
                  DO j=1,4
                    igo(2*j-1,1)=2*index(j,ipp)-1
                    igo(2*j,1)=2*index(j,ipa)
                  ENDDO
                ENDIF 
              ENDIF
              
              IF(idwm(1,i).NE.idwm(2,i))THEN
                local_index=15
                ialfa(local_index)=1
                n_kphs=n_kphs+1
                
                iorder(1,local_index)=iorder(1,local_index-1) !W-in ->W+
                iorder(2,local_index)=iorder(2,local_index-1)
                iorder(3,local_index)=iorder(7,local_index-1) ! W-
                iorder(4,local_index)=iorder(8,local_index-1)  
                iorder(5,local_index)=iorder(5,local_index-1) ! W+
                iorder(6,local_index)=iorder(6,local_index-1)
                iorder(7,local_index)=iorder(3,local_index-1) ! W-
                iorder(8,local_index)=iorder(4,local_index-1)
              ENDIF
            ELSEIF(nz0f(i).EQ.3)THEN
c Z2 and Z3 make a Higgs
*              WRITE(*,*) i,' Z -> Z Z Z  : 20,21,22'
              local_index=20
              ialfa(local_index)=1  
              n_kphs=n_kphs+1
              
              IF(ipartinitial(1).EQ.1)THEN
                iorder(1,local_index)=1    ! Z 
                iorder(2,local_index)=2
              ELSE
                iorder(1,local_index)=2 
                iorder(2,local_index)=1
              ENDIF
              iorder(3,local_index)=3    ! Z2
              iorder(4,local_index)=4  
              iorder(5,local_index)=5    ! Z3
              iorder(6,local_index)=6
              iorder(7,local_index)=7    ! Z4
              iorder(8,local_index)=8
              
              IF(ipartinitial(1).EQ.1)THEN
                idplocal(iorder(1,local_index))=idz0(1,i)
                idplocal(iorder(2,local_index))=-idz0(1,i)
              ELSE
                idplocal(iorder(1,local_index))=-idz0(1,i)
                idplocal(iorder(2,local_index))=idz0(1,i)
              ENDIF
              idplocal(iorder(3,local_index))=idz0(2,i)    ! Z
              idplocal(iorder(4,local_index))=-idz0(2,i)  
              idplocal(iorder(5,local_index))=idz0(3,i)    ! Z
              idplocal(iorder(6,local_index))=-idz0(3,i)
              idplocal(iorder(7,local_index))=idz0(4,i)    ! Z
              idplocal(iorder(8,local_index))=-idz0(4,i)

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(4)=1
                DO j=1,8
                  igo(j,4)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                     idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(4).EQ.0)THEN
                  igr(4)=1
                  DO j=1,4
                    igo(2*j-1,4)=2*index(j,ipp)-1
                    igo(2*j,4)=2*index(j,ipa)
                  ENDDO
                ENDIF 
              ENDIF
              
              IF(idz0(3,i).NE.idz0(4,i))THEN
c Z2 and Z4 make a Higgs
                local_index=21
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=iorder(1,local_index-1) 
                iorder(2,local_index)=iorder(2,local_index-1)
                iorder(3,local_index)=iorder(3,local_index-1)    ! Z2
                iorder(4,local_index)=iorder(4,local_index-1)  
                iorder(5,local_index)=iorder(7,local_index-1)    ! Z4
                iorder(6,local_index)=iorder(8,local_index-1)
                iorder(7,local_index)=iorder(5,local_index-1)    ! Z3
                iorder(8,local_index)=iorder(6,local_index-1)
              ENDIF

              IF(idz0(2,i).NE.idz0(3,i))THEN
c Z3 and Z4 make a Higgs
                local_index=22
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=iorder(1,local_index-2) 
                iorder(2,local_index)=iorder(2,local_index-2)
                iorder(3,local_index)=iorder(5,local_index-2)    ! Z3
                iorder(4,local_index)=iorder(6,local_index-2)  
                iorder(5,local_index)=iorder(7,local_index-2)    ! Z4
                iorder(6,local_index)=iorder(8,local_index-2)
                iorder(7,local_index)=iorder(3,local_index-2)    ! Z2
                iorder(8,local_index)=iorder(4,local_index-2)
              ENDIF
            ELSEIF(nz0f(i).EQ.2)THEN 
*              WRITE(*,*) i,' W -> Z Z W  : 19'
              local_index=19
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(3,local_index)=1 ! Z
              iorder(4,local_index)=2  
              iorder(5,local_index)=3 ! Z
              iorder(6,local_index)=4
              IF(iwpinout(1,i).EQ.-2)THEN
                IF(ipartinitial(1).EQ.1)THEN
                  iorder(1,local_index)=5 ! W+
                  iorder(2,local_index)=6
                ELSE  
                  iorder(1,local_index)=6    
                  iorder(2,local_index)=5
                ENDIF  
                iorder(7,local_index)=7 ! W-
                iorder(8,local_index)=8  
              ELSE
                iorder(7,local_index)=5 ! W+
                iorder(8,local_index)=6  
                IF(ipartinitial(1).EQ.1)THEN
                  iorder(1,local_index)=7 ! W-
                  iorder(2,local_index)=8
                ELSE  
                  iorder(1,local_index)=8
                  iorder(2,local_index)=7
                ENDIF  
              ENDIF 

              idplocal(iorder(3,local_index))=idz0(1,i)
              idplocal(iorder(4,local_index))=-idz0(1,i)
              idplocal(iorder(5,local_index))=idz0(2,i)
              idplocal(iorder(6,local_index))=-idz0(2,i)
              IF(iwpinout(1,i).EQ.-2)THEN
                IF(ipartinitial(1).EQ.1)THEN
                  idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(2,local_index))=-(idwp(1,i)-1)  
                ELSE  
                  idplocal(iorder(1,local_index))=-(idwp(1,i)-1)
                  idplocal(iorder(2,local_index))=idwp(1,i)  
                ENDIF  
                idplocal(iorder(7,local_index))=idwm(1,i) ! W-
                idplocal(iorder(8,local_index))=-(idwm(1,i)+1)  
              ELSE
                idplocal(iorder(7,local_index))=idwp(1,i) ! W+
                idplocal(iorder(8,local_index))=-(idwp(1,i)-1)  
                IF(ipartinitial(1).EQ.1)THEN
                  idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
                ELSE  
                  idplocal(iorder(1,local_index))=-(idwm(1,i)+1)
                  idplocal(iorder(2,local_index))=idwm(1,i)
                ENDIF  
              ENDIF 

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(2)=1
                DO j=1,8
                  igo(j,2)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
                idwm_ref=idwm(1,i)
                idwp_ref=idwp(1,i)
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(2).EQ.0)THEN
                  igr(2)=1
                  DO j=1,4
                    igo(2*j-1,2)=2*index(j,ipp)-1
                    igo(2*j,2)=2*index(j,ipa)
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                  IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                    DO j=1,4
                      igo(2*j-1,3)=2*index(j,ipp)-1
                      igo(2*j,3)=2*index(j,ipa)
                    ENDDO
                    igr(3)=1       
                  ENDIF
                ENDIF 
              ENDIF
            ELSE  ! Z+Wm+Wp
*             WRITE(*,*) i,' Z  -> Z  W+ W-  : 18'
              local_index=18
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              IF(ipartinitial(1).EQ.1)THEN
                iorder(1,local_index)=1    ! Z 
                iorder(2,local_index)=2
              ELSE
                iorder(1,local_index)=2 
                iorder(2,local_index)=1
              ENDIF
              iorder(3,local_index)=5    ! W+
              iorder(4,local_index)=6  
              iorder(5,local_index)=7    ! W-
              iorder(6,local_index)=8
              iorder(7,local_index)=3    ! Z
              iorder(8,local_index)=4
              
              IF(ipartinitial(1).EQ.1)THEN
                idplocal(iorder(1,local_index))=idz0(1,i)
                idplocal(iorder(2,local_index))=-idz0(1,i)
              ELSE
                idplocal(iorder(1,local_index))=-idz0(1,i)
                idplocal(iorder(2,local_index))=idz0(1,i)
              ENDIF
              idplocal(iorder(7,local_index))=idz0(2,i)
              idplocal(iorder(8,local_index))=-idz0(2,i)
              idplocal(iorder(3,local_index))=idwp(1,i)    ! W+
              idplocal(iorder(4,local_index))=-(idwp(1,i)-1)  
              idplocal(iorder(5,local_index))=idwm(1,i)    ! W-
              idplocal(iorder(6,local_index))=-(idwm(1,i)+1)

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(2)=1
                DO j=1,8
                  igo(j,2)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
                idwm_ref=idwm(1,i)
                idwp_ref=idwp(1,i)
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(2).EQ.0)THEN
                  igr(2)=1
                  DO j=1,4
                    igo(2*j-1,2)=2*index(j,ipp)-1
                    igo(2*j,2)=2*index(j,ipa)
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                  IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                    DO j=1,4
                      igo(2*j-1,3)=2*index(j,ipp)-1
                      igo(2*j,3)=2*index(j,ipa)
                    ENDDO
                    igr(3)=1       
                  ENDIF
                ENDIF 
              ENDIF
            ENDIF
c Check for top
            IF(nwmf(i).EQ.1 .AND. nantib.GE.1)THEN 
              IF(nwpf(i).EQ.1 .AND. nb.GE.1)THEN
*               WRITE(*,*) i,' Z  -> t  t  : 56'
* The initial Z appears first, The final (bb) one second.
                local_index=56
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(ipartinitial(1).EQ.1)THEN
                  iorder(1,local_index)=1    ! Z 
                  iorder(2,local_index)=2
                ELSE
                  iorder(1,local_index)=2 
                  iorder(2,local_index)=1
                ENDIF
                iorder(7,local_index)=3    ! Z  b
                iorder(8,local_index)=4    !    bbar
                iorder(3,local_index)=5    ! W+
                iorder(4,local_index)=6  
                iorder(5,local_index)=7    ! W-
                iorder(6,local_index)=8

                IF(ipartinitial(1).EQ.1)THEN
                  idplocal(iorder(1,local_index))=idz0(1,i)
                  idplocal(iorder(2,local_index))=-idz0(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=-idz0(1,i)
                  idplocal(iorder(2,local_index))=idz0(1,i)
                ENDIF
                idplocal(iorder(7,local_index))=idz0(2,i)    ! Z  b
                idplocal(iorder(8,local_index))= -idz0(2,i)   !    bbar
                idplocal(iorder(3,local_index))=idwp(1,i)    ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1) 
                idplocal(iorder(5,local_index))=idwm(1,i)    ! W-
                idplocal(iorder(6,local_index))=-(idwm(1,i)+1)

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1      
                    ENDIF
                  ENDIF 
                ENDIF
              ELSE ! Two final Z's, one in b-bbar
*                WRITE(*,*) i,' W- -> tbar  Z b  : 55'
* There are two final state Z's. The bb one has to be located
                local_index=55
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(5,local_index)=1 ! Z 
                iorder(6,local_index)=2
                iorder(8,local_index)=3 ! Z  b
                iorder(7,local_index)=4 !    bbar
                iorder(3,local_index)=7 ! W-
                iorder(4,local_index)=8  
                IF(ipartinitial(1).EQ.1)THEN
                  iorder(1,local_index)=5 ! W- -> W+
                  iorder(2,local_index)=6
                ELSE
                  iorder(1,local_index)=6    
                  iorder(2,local_index)=5
                ENDIF

                IF(idz0(2,i).EQ.5)THEN
                  idaux1=1
                  idaux2=2
                ELSE
                  idaux1=2       
                  idaux2=1
                ENDIF
                idplocal(iorder(5,local_index))=idz0(idaux1,i)
                idplocal(iorder(6,local_index))=-idz0(idaux1,i)
                idplocal(iorder(8,local_index))=idz0(idaux2,i) ! Z  b
                idplocal(iorder(7,local_index))=-idz0(idaux2,i) ! bbar
                idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                IF(ipartinitial(1).EQ.1)THEN
                  idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
                ELSE
                  idplocal(iorder(1,local_index))=-(idwp(1,i)-1)
                  idplocal(iorder(2,local_index))=idwp(1,i)
                ENDIF

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN 
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1       
                    ENDIF 
                  ENDIF
                ENDIF
              ENDIF
            ELSEIF(nwpf(i).EQ.1 .AND. nb.GE.1)THEN
*             WRITE(*,*) i,' W+ -> t  Z bbar  : 55'
              local_index=55
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(5,local_index)=1 ! Z 
              iorder(6,local_index)=2
              iorder(7,local_index)=3 ! Z  b
              iorder(8,local_index)=4 !    bbar
              iorder(3,local_index)=5 ! W+
              iorder(4,local_index)=6  
              IF(ipartinitial(1).EQ.1)THEN
                iorder(1,local_index)=7 ! W+ -> W-
                iorder(2,local_index)=8
              ELSE
                iorder(1,local_index)=8
                iorder(2,local_index)=7
              ENDIF

              IF(idz0(2,i).EQ.5)THEN
                idaux1=1
                idaux2=2
              ELSE
                idaux1=2       
                idaux2=1
              ENDIF
              idplocal(iorder(5,local_index))=idz0(idaux1,i) ! Z 
              idplocal(iorder(6,local_index))=-idz0(idaux1,i)
              idplocal(iorder(7,local_index))=idz0(idaux2,i) ! Z  b
              idplocal(iorder(8,local_index))=-idz0(idaux2,i) !    bbar
              idplocal(iorder(3,local_index))=idwp(1,i) ! W+
              idplocal(iorder(4,local_index))=-(idwp(1,i)-1)  
              IF(ipartinitial(1).EQ.1)THEN
                idplocal(iorder(1,local_index))=idwm(1,i) ! W+ -> W-
                idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
              ELSE
                idplocal(iorder(1,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(2,local_index))=idwm(1,i)
              ENDIF

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(2)=1
                DO j=1,8
                  igo(j,2)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
                idwm_ref=idwm(1,i)
                idwp_ref=idwp(1,i)
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(2).EQ.0)THEN
                  igr(2)=1
                  DO j=1,4
                    igo(2*j-1,2)=2*index(j,ipp)-1
                    igo(2*j,2)=2*index(j,ipa)
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                  IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                    DO j=1,4
                      igo(2*j-1,3)=2*index(j,ipp)-1
                      igo(2*j,3)=2*index(j,ipa)
                    ENDDO
                    igr(3)=1       
                  ENDIF
                ENDIF 
              ENDIF
            ENDIF               ! end Check for top
          ELSE                  ! nfinal=2
c 1_1_4 or 1_1_31 or 1_2_3

c Check for type of the two t-channel bosons.
c I am sure that the two initial particles are in different bosons.
c Determine if W+ is u-->W+d or  d~-->W+u~ etc
c iswap(i)=1,2 records the order in which the two bosons which include 
c the initial particles appear in the list. Necessary only when the two 
c bosons are identical. We first check the order in which they appear 
c in the basic order. Then whether they have been interchanged by the 
c ordering by flavour.

c Notice that when pairing particles and antiparticles only particles 
c are rotated.
            IF(ipartinitial(1).EQ.-1)THEN  ! 1=antiparticle #1
              itypeinitial(1)=itype(1,i)
              IF(ipartinitial(2).EQ.-1)THEN  ! 2=antiparticle #2
                itypeinitial(2)=itype(2,i)
              ELSE                          ! 2=particle #1  
                itypeinitial(2)=itype(invindex(1,iperm(i)),i)
              ENDIF
              iflagswap=0
            ELSE                          ! 1=particle #1
              itypeinitial(1)=itype(invindex(1,iperm(i)),i)
              IF(ipartinitial(2).EQ.-1)THEN  ! 2=antiparticle #1
                itypeinitial(2)=itype(1,i)
                iflagswap=1
              ELSE                          ! 2=particle #2  
                itypeinitial(2)=itype(invindex(2,iperm(i)),i)
                IF(invindex(1,iperm(i)).LT.invindex(2,iperm(i)))THEN
                  iflagswap=0
                ELSE
                  iflagswap=1
                ENDIF
              ENDIF
            ENDIF   ! end  Check for type of t-channel boson        

c If iflagswap=0 the boson containing initial-particle-1 comes before
c the one containing initial-particle-2, after if iflagswap=1
c Now check if flavour ordering has interchanged the two bosons.
           
            IF(itypeinitial(1).EQ.itypeinitial(2))THEN
              IF(itypeinitial(1).EQ.1)THEN
                IF(iposwp(2,i).LT.iposwp(1,i)) iflagswap=1-iflagswap
              ELSEIF(itypeinitial(1).EQ.-1)THEN
                IF(iposwm(2,i).LT.iposwm(1,i)) iflagswap=1-iflagswap
              ELSE
                IF(iposz0(2,i).LT.iposz0(1,i)) iflagswap=1-iflagswap
              ENDIF
            ENDIF
            
            IF(iflagswap.EQ.0)THEN
              iswap(1)=1
              iswap(2)=2
            ELSE
              iswap(1)=2
              iswap(2)=1
            ENDIF
            
            IF(nwpf(i).EQ.2)THEN
*              WRITE(*,*) i,' t(WW) W+ W+  : 4'
              local_index=4
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(3,local_index)=1    ! W+
              iorder(4,local_index)=2 
              iorder(5,local_index)=5    ! W+
              iorder(6,local_index)=6 
              IF(ipartinitial(1).EQ.-1)THEN
                iorder(1,local_index)=4
                iorder(7,local_index)=3
              ELSE
                iorder(1,local_index)=3
                iorder(7,local_index)=4
              ENDIF
              IF(ipartinitial(2).EQ.-1)THEN
                iorder(2,local_index)=8
                iorder(8,local_index)=7
              ELSE
                iorder(2,local_index)=7
                iorder(8,local_index)=8
              ENDIF

              idplocal(iorder(3,local_index))=idwp(1,i)    ! W+
              idplocal(iorder(4,local_index))=-(idwp(1,i)-1) 
              idplocal(iorder(5,local_index))=idwp(2,i)    ! W+
              idplocal(iorder(6,local_index))=-(idwp(2,i)-1)
              IF(ipartinitial(1).EQ.-1)THEN
                idplocal(iorder(1,local_index))=-(idwm(iswap(1),i)+1)
                idplocal(iorder(7,local_index))=idwm(iswap(1),i)
              ELSE
                idplocal(iorder(1,local_index))=idwm(iswap(1),i)
                idplocal(iorder(7,local_index))=-(idwm(iswap(1),i)+1)
              ENDIF
              IF(ipartinitial(2).EQ.-1)THEN
                idplocal(iorder(2,local_index))=-(idwm(iswap(2),i)+1)
                idplocal(iorder(8,local_index))=idwm(iswap(2),i)
              ELSE
                idplocal(iorder(2,local_index))=idwm(iswap(2),i)
                idplocal(iorder(8,local_index))=-(idwm(iswap(2),i)+1)
              ENDIF

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(1)=1
                DO j=1,8
                  igo(j,1)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(1).EQ.0)THEN
                  igr(1)=1
                  DO j=1,4
                    igo(2*j-1,1)=2*index(j,ipp)-1
                    igo(2*j,1)=2*index(j,ipa)
                  ENDDO
                ENDIF 
              ENDIF
            ELSEIF(nwmf(i).EQ.2)THEN
*              WRITE(*,*) i,' t(WW) W- W-  : 4'
              local_index=4
              ialfa(local_index)=1
              n_kphs=n_kphs+1
              
              iorder(3,local_index)=3    ! W-
              iorder(4,local_index)=4 
              iorder(5,local_index)=7    ! W-
              iorder(6,local_index)=8 
              IF(ipartinitial(1).EQ.-1)THEN
                iorder(1,local_index)=2
                iorder(7,local_index)=1
              ELSE
                iorder(1,local_index)=1
                iorder(7,local_index)=2
              ENDIF
              IF(ipartinitial(2).EQ.-1)THEN
                iorder(2,local_index)=6
                iorder(8,local_index)=5
              ELSE
                iorder(2,local_index)=5
                iorder(8,local_index)=6
              ENDIF
                
              idplocal(iorder(3,local_index))=idwm(1,i)    ! W-
              idplocal(iorder(4,local_index))=-(idwm(1,i)+1) 
              idplocal(iorder(5,local_index))=idwm(2,i)    ! W-
              idplocal(iorder(6,local_index))=-(idwm(2,i)+1) 
              IF(ipartinitial(1).EQ.-1)THEN
                idplocal(iorder(1,local_index))=-(idwp(iswap(1),i)-1)
                idplocal(iorder(7,local_index))=idwp(iswap(1),i)
              ELSE
                idplocal(iorder(1,local_index))=idwp(iswap(1),i)
                idplocal(iorder(7,local_index))=-(idwp(iswap(1),i)-1)
              ENDIF
              IF(ipartinitial(2).EQ.-1)THEN
                idplocal(iorder(2,local_index))=-(idwp(iswap(2),i)-1)
                idplocal(iorder(8,local_index))=idwp(iswap(2),i)
              ELSE
                idplocal(iorder(2,local_index))=idwp(iswap(2),i)
                idplocal(iorder(8,local_index))=-(idwp(iswap(2),i)-1)
              ENDIF
                
              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(1)=1
                DO j=1,8
                  igo(j,1)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(1).EQ.0)THEN
                  igr(1)=1
                  DO j=1,4
                    igo(2*j-1,1)=2*index(j,ipp)-1
                    igo(2*j,1)=2*index(j,ipa)
                  ENDDO
                ENDIF 
              ENDIF
            ELSEIF(nwmf(i).EQ.1 .AND. nwpf(i).EQ.1)THEN
              IF(nz0(i).EQ.2)THEN
*                WRITE(*,*) i,' t(00) W+ W-  : 1'
                local_index=1
                ialfa(local_index)=1
                n_kphs=n_kphs+1
                
                iorder(3,local_index)=5    ! W+
                iorder(4,local_index)=6 
                iorder(5,local_index)=7    ! W-
                iorder(6,local_index)=8 
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=2
                  iorder(7,local_index)=1
                ELSE
                  iorder(1,local_index)=1
                  iorder(7,local_index)=2
                ENDIF
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=4
                  iorder(8,local_index)=3
                ELSE
                  iorder(2,local_index)=3
                  iorder(8,local_index)=4
                ENDIF
                  
                idplocal(iorder(3,local_index))=idwp(1,i)   ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1) 
                idplocal(iorder(5,local_index))=idwm(1,i)    ! W-
                idplocal(iorder(6,local_index))=-(idwm(1,i)+1) 
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-idz0(iswap(1),i)
                  idplocal(iorder(7,local_index))=idz0(iswap(1),i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(iswap(1),i)
                  idplocal(iorder(7,local_index))=-idz0(iswap(1),i)
                ENDIF
                IF(ipartinitial(2).EQ.-1)THEN
                  idplocal(iorder(2,local_index))=-idz0(iswap(2),i)
                  idplocal(iorder(8,local_index))=idz0(iswap(2),i)
                ELSE
                  idplocal(iorder(2,local_index))=idz0(iswap(2),i)
                  idplocal(iorder(8,local_index))=-idz0(iswap(2),i)
                ENDIF
                  
                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1       
                    ENDIF
                  ENDIF 
                ENDIF
              ELSE
*                WRITE(*,*) i,' t(WW) W+ W-  : 5'
                local_index=5
                ialfa(local_index)=1
                n_kphs=n_kphs+1
                
                iorder(3,local_index)=5    ! W+
                iorder(4,local_index)=6 
                iorder(5,local_index)=7    ! W-
                iorder(6,local_index)=8 
                IF(itypeinitial(1).EQ.1)THEN  ! First W is +
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=2
                    iorder(7,local_index)=1
                  ELSE
                    iorder(1,local_index)=1
                    iorder(7,local_index)=2
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN
                    iorder(2,local_index)=4
                    iorder(8,local_index)=3
                  ELSE
                    iorder(2,local_index)=3
                    iorder(8,local_index)=4
                  ENDIF
                ELSE  ! First W is -
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=4
                    iorder(7,local_index)=3
                  ELSE
                    iorder(1,local_index)=3
                    iorder(7,local_index)=4
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN
                    iorder(2,local_index)=2
                    iorder(8,local_index)=1
                  ELSE
                    iorder(2,local_index)=1
                    iorder(8,local_index)=2
                  ENDIF
                ENDIF
                  
                idplocal(iorder(3,local_index))=idwp(2,i)    ! W+
                idplocal(iorder(4,local_index))=-(idwp(2,i)-1) 
                idplocal(iorder(5,local_index))=idwm(2,i)    ! W-
                idplocal(iorder(6,local_index))=-(idwm(2,i)+1) 
                IF(itypeinitial(1).EQ.1)THEN  ! First W is +
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwp(1,i)-1)
                    idplocal(iorder(7,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwp(1,i)
                    idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
                    idplocal(iorder(8,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idwm(1,i)
                    idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                  ENDIF
                ELSE  ! First W is -
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwm(1,i)+1)
                    idplocal(iorder(7,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwm(1,i)
                    idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
                    idplocal(iorder(8,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idwp(1,i)
                    idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ENDIF

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(1)=1
                  DO j=1,8
                    igo(j,1)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                   IF(igr(1).EQ.0)THEN
                     igr(1)=1
                     DO j=1,4
                       igo(2*j-1,1)=2*index(j,ipp)-1
                       igo(2*j,1)=2*index(j,ipa)
                     ENDDO
                   ENDIF 
                ENDIF
              ENDIF
            ELSEIF(nwmf(i).EQ.1)THEN
              IF(itypeinitial(1).EQ.0 .OR. idout(1).EQ.idout(2))THEN
*                WRITE(*,*) i,' t(0W) Z W-  : 2'
                local_index=2
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(3,local_index)=3    ! Z
                iorder(4,local_index)=4 
                iorder(5,local_index)=7    ! W-
                iorder(6,local_index)=8 
                IF(ipartinitial(1).EQ.-1)THEN  ! Z
                  iorder(1,local_index)=2
                  iorder(7,local_index)=1
                ELSE
                  iorder(1,local_index)=1
                  iorder(7,local_index)=2
                ENDIF
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=6
                  iorder(8,local_index)=5
                ELSE
                  iorder(2,local_index)=5
                  iorder(8,local_index)=6
                ENDIF


                idplocal(iorder(3,local_index))=idz0(2,i)    ! Z
                idplocal(iorder(4,local_index))=-idz0(2,i) 
                idplocal(iorder(5,local_index))=idwm(1,i)    ! W-
                idplocal(iorder(6,local_index))=-(idwm(1,i)+1) 
                IF(ipartinitial(1).EQ.-1)THEN      ! Z
                  idplocal(iorder(1,local_index))=-idz0(1,i)
                  idplocal(iorder(7,local_index))=idz0(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(1,i)
                  idplocal(iorder(7,local_index))=-idz0(1,i)
                ENDIF
                IF(ipartinitial(2).EQ.-1)THEN      ! W+
                  idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
                  idplocal(iorder(8,local_index))=idwp(1,i)
                ELSE
                  idplocal(iorder(2,local_index))=idwp(1,i)
                  idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                ENDIF

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1       
                    ENDIF
                  ENDIF 
                ENDIF
              ENDIF
              IF(itypeinitial(2).EQ.0 .OR. idout(1).EQ.idout(2))THEN
*                WRITE(*,*) i,' t(W0) Z W-  : 3'
                local_index=3
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(idout(1).NE.idout(2))THEN
                  iorder(3,local_index)=3    ! Z
                  iorder(4,local_index)=4 
                  iorder(5,local_index)=7    ! W-
                  iorder(6,local_index)=8 
                  IF(ipartinitial(1).EQ.-1)THEN  ! W+
                    iorder(1,local_index)=6
                    iorder(7,local_index)=5
                  ELSE
                    iorder(1,local_index)=5
                    iorder(7,local_index)=6
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN  ! Z
                    iorder(2,local_index)=2
                    iorder(8,local_index)=1
                  ELSE
                    iorder(2,local_index)=1
                    iorder(8,local_index)=2
                  ENDIF


                  idplocal(iorder(3,local_index))=idz0(2,i)    ! Z
                  idplocal(iorder(4,local_index))=-idz0(2,i) 
                  idplocal(iorder(5,local_index))=idwm(1,i)    ! W-
                  idplocal(iorder(6,local_index))=-(idwm(1,i)+1) 
                  IF(ipartinitial(1).EQ.-1)THEN      ! W+
                    idplocal(iorder(1,local_index))=-(idwp(1,i)-1)
                    idplocal(iorder(7,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwp(1,i)
                    idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN      ! Z
                    idplocal(iorder(2,local_index))=-idz0(1,i)
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idz0(1,i)
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF

                  IF(n_kphs.eq.1)THEN
                   iabase=local_index
                    igr(2)=1
                    DO j=1,8
                      igo(j,2)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(2).EQ.0)THEN
                      igr(2)=1
                      DO j=1,4
                        igo(2*j-1,2)=2*index(j,ipp)-1
                        igo(2*j,2)=2*index(j,ipa)
                      ENDDO
                      idwm_ref=idwm(1,i)
                      idwp_ref=idwp(1,i)
                    ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                      IF(idwm(1,i).NE.idwm_ref .OR. 
     &                                     idwp(1,i).NE.idwp_ref)THEN
                        DO j=1,4
                          igo(2*j-1,3)=2*index(j,ipp)-1
                          igo(2*j,3)=2*index(j,ipa)
                        ENDDO
                        igr(3)=1      
                      ENDIF 
                    ENDIF 
                  ENDIF
                ELSE  ! iabase should be 2
                  local_index=3
                  iorder(1,local_index)=iorder(2,local_index-1)
                  iorder(2,local_index)=iorder(1,local_index-1)
                  iorder(3,local_index)=iorder(3,local_index-1)
                  iorder(4,local_index)=iorder(4,local_index-1)
                  iorder(5,local_index)=iorder(5,local_index-1)
                  iorder(6,local_index)=iorder(6,local_index-1)
                  iorder(7,local_index)=iorder(8,local_index-1)
                  iorder(8,local_index)=iorder(7,local_index-1)
                ENDIF
              ENDIF
            ELSEIF(nwpf(i).EQ.1)THEN
              IF(itypeinitial(1).EQ.0.OR. idout(1).EQ.idout(2))THEN
*                WRITE(*,*) i,' t(0W) Z W+  : 2'
                local_index=2
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(3,local_index)=3    ! Z
                iorder(4,local_index)=4 
                iorder(5,local_index)=5    ! W+
                iorder(6,local_index)=6 
                IF(ipartinitial(1).EQ.-1)THEN  ! Z
                  iorder(1,local_index)=2
                  iorder(7,local_index)=1
                ELSE
                  iorder(1,local_index)=1
                  iorder(7,local_index)=2
                ENDIF
                IF(ipartinitial(2).EQ.-1)THEN  ! W-
                  iorder(2,local_index)=8
                  iorder(8,local_index)=7
                ELSE
                  iorder(2,local_index)=7
                  iorder(8,local_index)=8
                ENDIF


                idplocal(iorder(3,local_index))=idz0(2,i)    ! Z
                idplocal(iorder(4,local_index))=-idz0(2,i) 
                idplocal(iorder(5,local_index))=idwp(1,i)    ! W+
                idplocal(iorder(6,local_index))=-(idwp(1,i)-1) 
                IF(ipartinitial(1).EQ.-1)THEN      ! Z
                  idplocal(iorder(1,local_index))=-idz0(1,i)
                  idplocal(iorder(7,local_index))=idz0(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(1,i)
                  idplocal(iorder(7,local_index))=-idz0(1,i)
                ENDIF
                IF(ipartinitial(2).EQ.-1)THEN      ! W-
                  idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
                  idplocal(iorder(8,local_index))=idwm(1,i)
                ELSE
                  idplocal(iorder(2,local_index))=idwm(1,i)
                  idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                ENDIF

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1       
                    ENDIF
                  ENDIF 
                ENDIF
              ENDIF
              IF(itypeinitial(2).EQ.0 .OR. idout(1).EQ.idout(2))THEN
*                WRITE(*,*) i,' t(W0) Z W+  : 3'
                local_index=3
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(idout(1).NE.idout(2))THEN
                  iorder(3,local_index)=3    ! Z
                  iorder(4,local_index)=4 
                  iorder(5,local_index)=5    ! W+
                  iorder(6,local_index)=6 
                  IF(ipartinitial(1).EQ.-1)THEN  ! W-
                    iorder(1,local_index)=8
                    iorder(7,local_index)=7
                  ELSE
                    iorder(1,local_index)=7
                    iorder(7,local_index)=8
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN  ! Z
                    iorder(2,local_index)=2
                    iorder(8,local_index)=1
                  ELSE
                    iorder(2,local_index)=1
                    iorder(8,local_index)=2
                  ENDIF


                  idplocal(iorder(3,local_index))=idz0(2,i)    ! Z
                  idplocal(iorder(4,local_index))=-idz0(2,i) 
                  idplocal(iorder(5,local_index))=idwp(1,i)    ! W+
                  idplocal(iorder(6,local_index))=-(idwp(1,i)-1) 
                  IF(ipartinitial(1).EQ.-1)THEN      ! W-
                    idplocal(iorder(1,local_index))=-(idwm(1,i)+1)
                    idplocal(iorder(7,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwm(1,i)
                    idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN      ! Z
                    idplocal(iorder(2,local_index))=-idz0(1,i)
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idz0(1,i)
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF

                  IF(n_kphs.eq.1)THEN
                   iabase=local_index
                    igr(2)=1
                    DO j=1,8
                      igo(j,2)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(2).EQ.0)THEN
                      igr(2)=1
                      DO j=1,4
                        igo(2*j-1,2)=2*index(j,ipp)-1
                        igo(2*j,2)=2*index(j,ipa)
                      ENDDO
                      idwm_ref=idwm(1,i)
                      idwp_ref=idwp(1,i)
                    ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                      IF(idwm(1,i).NE.idwm_ref .OR. 
     &                                     idwp(1,i).NE.idwp_ref)THEN
                        DO j=1,4
                          igo(2*j-1,3)=2*index(j,ipp)-1
                          igo(2*j,3)=2*index(j,ipa)
                        ENDDO
                        igr(3)=1       
                      ENDIF
                    ENDIF 
                  ENDIF
                ELSE  ! iabase should be 2
                  local_index=3
                  iorder(1,local_index)=iorder(2,local_index-1)
                  iorder(2,local_index)=iorder(1,local_index-1)
                  iorder(3,local_index)=iorder(3,local_index-1)
                  iorder(4,local_index)=iorder(4,local_index-1)
                  iorder(5,local_index)=iorder(5,local_index-1)
                  iorder(6,local_index)=iorder(6,local_index-1)
                  iorder(7,local_index)=iorder(8,local_index-1)
                  iorder(8,local_index)=iorder(7,local_index-1)
                ENDIF
              ENDIF
            ELSE    !(nz0f(i).EQ.2)
              IF(nz0(i).eq.4)THEN
*                WRITE(*,*) i,' t(00) Z Z  : 6'
                local_index=6
                ialfa(local_index)=1
                n_kphs=n_kphs+1
                
                iorder(3,local_index)=5    ! Z
                iorder(4,local_index)=6 
                iorder(5,local_index)=7    ! Z
                iorder(6,local_index)=8 
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=2
                  iorder(7,local_index)=1
                ELSE
                  iorder(1,local_index)=1
                  iorder(7,local_index)=2
                ENDIF
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=4
                  iorder(8,local_index)=3
                ELSE
                  iorder(2,local_index)=3
                  iorder(8,local_index)=4
                ENDIF
                  
                idplocal(iorder(3,local_index))=idz0(3,i)   ! Z
                idplocal(iorder(4,local_index))=-idz0(3,i) 
                idplocal(iorder(5,local_index))=idz0(4,i)    ! Z
                idplocal(iorder(6,local_index))=-idz0(4,i) 
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-idz0(iswap(1),i)
                  idplocal(iorder(7,local_index))=idz0(iswap(1),i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(iswap(1),i)
                  idplocal(iorder(7,local_index))=-idz0(iswap(1),i)
                ENDIF
                IF(ipartinitial(2).EQ.-1)THEN
                  idplocal(iorder(2,local_index))=-idz0(iswap(2),i)
                  idplocal(iorder(8,local_index))=idz0(iswap(2),i)
                ELSE
                  idplocal(iorder(2,local_index))=idz0(iswap(2),i)
                  idplocal(iorder(8,local_index))=-idz0(iswap(2),i)
                ENDIF
                  
                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(4)=1
                  DO j=1,8
                    igo(j,4)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(4).EQ.0)THEN
                    igr(4)=1
                    DO j=1,4
                      igo(2*j-1,4)=2*index(j,ipp)-1
                      igo(2*j,4)=2*index(j,ipa)
                    ENDDO
                  ENDIF 
                ENDIF
              ELSE 
*                WRITE(*,*) i,' t(WW) Z Z  : 7'
                local_index=7
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(3,local_index)=1    ! Z
                iorder(4,local_index)=2 
                iorder(5,local_index)=3    ! Z
                iorder(6,local_index)=4 
                IF(itypeinitial(1).EQ.1)THEN  ! First W is +
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=6
                    iorder(7,local_index)=5
                  ELSE
                    iorder(1,local_index)=5
                    iorder(7,local_index)=6
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN
                    iorder(2,local_index)=8
                    iorder(8,local_index)=7
                  ELSE
                    iorder(2,local_index)=7
                    iorder(8,local_index)=8
                  ENDIF
                ELSE  ! First W is -
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=8
                    iorder(7,local_index)=7
                  ELSE
                    iorder(1,local_index)=7
                    iorder(7,local_index)=8
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN
                    iorder(2,local_index)=6
                    iorder(8,local_index)=5
                  ELSE
                    iorder(2,local_index)=5
                    iorder(8,local_index)=6
                  ENDIF
                ENDIF
                  
                idplocal(iorder(3,local_index))=idz0(1,i)    ! Z
                idplocal(iorder(4,local_index))=-idz0(1,i) 
                idplocal(iorder(5,local_index))=idz0(2,i)    ! Z
                idplocal(iorder(6,local_index))=-idz0(2,i) 
                IF(itypeinitial(1).EQ.1)THEN  ! First W is +
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwp(1,i)-1)
                    idplocal(iorder(7,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwp(1,i)
                    idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
                    idplocal(iorder(8,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idwm(1,i)
                    idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                  ENDIF
                ELSE  ! First W is -
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwm(1,i)+1)
                    idplocal(iorder(7,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwm(1,i)
                    idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
                    idplocal(iorder(8,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idwp(1,i)
                    idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ENDIF

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
* giuseppe bug correction 21/2/07
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
* giuseppe bug correction 21/2/07 end
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1       
                    ENDIF
                  ENDIF 
                ENDIF
              ENDIF
            ENDIF
c Check for top. These possibilities add to those previously found
c ZW+ in final state 
            IF(nwpf(i).EQ.1 .AND. nz0f(i).EQ.1)THEN
c Case 1 tbbar or tbarb --> W+Z(bbbar) in the final state         
              ifound=0
              DO j=1,nz0(i)
                IF(iz0_b(j,i).EQ.1 .AND. iz0inout (j,i).EQ.2) ifound=1
              ENDDO
              IF(ifound.EQ.1)THEN
                IF(itypeinitial(1).EQ.0 .OR. idout(1).EQ.idout(2))THEN
*                WRITE(*,*) i,' t(0W) bbar t -> Z W+  : 26'
                  local_index=26
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  iorder(5,local_index)=3    ! Z
                  iorder(6,local_index)=4 
                  iorder(3,local_index)=5    ! W+
                  iorder(4,local_index)=6 
                  IF(ipartinitial(1).EQ.-1)THEN  ! Z
                    iorder(1,local_index)=2
                    iorder(7,local_index)=1
                  ELSE
                    iorder(1,local_index)=1
                    iorder(7,local_index)=2
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN  ! W-
                    iorder(2,local_index)=8
                    iorder(8,local_index)=7
                  ELSE
                    iorder(2,local_index)=7
                    iorder(8,local_index)=8
                  ENDIF
                  
                  idplocal(iorder(5,local_index))=idz0(2,i)    ! Z
                  idplocal(iorder(6,local_index))=-idz0(2,i) 
                  idplocal(iorder(3,local_index))=idwp(1,i)    ! W+
                  idplocal(iorder(4,local_index))=-(idwp(1,i)-1) 
                  IF(ipartinitial(1).EQ.-1)THEN      ! Z
                    idplocal(iorder(1,local_index))=-idz0(1,i)
                    idplocal(iorder(7,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idz0(1,i)
                    idplocal(iorder(7,local_index))=-idz0(1,i)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN      ! W-
                    idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
                    idplocal(iorder(8,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idwm(1,i)
                    idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                  ENDIF

                  IF(n_kphs.eq.1)THEN
                    iabase=local_index
                    igr(2)=1
                    DO j=1,8
                      igo(j,2)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(2).EQ.0)THEN
                      igr(2)=1
                      DO j=1,4
                        igo(2*j-1,2)=2*index(j,ipp)-1
                        igo(2*j,2)=2*index(j,ipa)
                      ENDDO
                      idwm_ref=idwm(1,i)
                      idwp_ref=idwp(1,i)
                    ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                      IF(idwm(1,i).NE.idwm_ref .OR. 
     &                                     idwp(1,i).NE.idwp_ref)THEN
                        DO j=1,4
                          igo(2*j-1,3)=2*index(j,ipp)-1
                          igo(2*j,3)=2*index(j,ipa)
                        ENDDO
                        igr(3)=1       
                      ENDIF
                    ENDIF 
                  ENDIF
                ENDIF
                IF(itypeinitial(2).EQ.0 .OR. idout(1).EQ.idout(2))
     &                                             THEN   
*                 WRITE(*,*) i,' t(W0) bbar t -> Z W+  : 25'
                  local_index=25
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  iorder(5,local_index)=3    ! Z
                  iorder(6,local_index)=4 
                  iorder(3,local_index)=5    ! W+
                  iorder(4,local_index)=6 
                  IF(ipartinitial(1).EQ.-1)THEN  ! W-
                    iorder(1,local_index)=8
                    iorder(7,local_index)=7
                  ELSE
                    iorder(1,local_index)=7
                    iorder(7,local_index)=8
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN  ! Z
                    iorder(2,local_index)=2
                    iorder(8,local_index)=1
                  ELSE
                    iorder(2,local_index)=1
                    iorder(8,local_index)=2
                  ENDIF
                  
                  idplocal(iorder(5,local_index))=idz0(2,i)    ! Z
                  idplocal(iorder(6,local_index))=-idz0(2,i) 
                  idplocal(iorder(3,local_index))=idwp(1,i)    ! W+
                  idplocal(iorder(4,local_index))=-(idwp(1,i)-1) 
                  IF(ipartinitial(1).EQ.-1)THEN      ! W-
                    idplocal(iorder(1,local_index))=-(idwm(1,i)+1)
                    idplocal(iorder(7,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwm(1,i)
                    idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN      ! Z
                    idplocal(iorder(2,local_index))=-idz0(1,i)
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idz0(1,i)
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF

                  IF(n_kphs.eq.1)THEN
                    iabase=local_index
                    igr(2)=1
                    DO j=1,8
                      igo(j,2)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(2).EQ.0)THEN
                      igr(2)=1
                      DO j=1,4
                        igo(2*j-1,2)=2*index(j,ipp)-1
                        igo(2*j,2)=2*index(j,ipa)
                      ENDDO
                      idwm_ref=idwm(1,i)
                      idwp_ref=idwp(1,i)
                    ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                      IF(idwm(1,i).NE.idwm_ref .OR. 
     &                                     idwp(1,i).NE.idwp_ref)THEN
                        DO j=1,4
                          igo(2*j-1,3)=2*index(j,ipp)-1
                          igo(2*j,3)=2*index(j,ipa)
                        ENDDO
                        igr(3)=1       
                      ENDIF
                    ENDIF 
                  ENDIF
                ENDIF  
              ENDIF  ! end case 1
c Case 2  b T(W) --> t Z --> b W Z
              IF((nantib-nb).EQ.-1)THEN
                                        ! need b in final state  
*                WRITE(*,*) i,' t(W+) b -> Z t -> Z W+ b  : 33'
                local_index=33
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(6,local_index)=3    ! Z final
                iorder(7,local_index)=4
                iorder(3,local_index)=5    ! W+
                iorder(4,local_index)=6
                iorder(5,local_index)=1    ! b

                DO j=1,2
                  IF(itypeinitial(j).EQ.-1) iwin=j
                ENDDO

                IF(abs(idout(1)).EQ.5)THEN
                  iorder(1,local_index)=2
                  IF(ipartinitial(iwin).EQ.-1)THEN
                    iorder(2,local_index)=8
                    iorder(8,local_index)=7
                  ELSE
                    iorder(2,local_index)=7
                    iorder(8,local_index)=8
                  ENDIF
                ELSE
                  iorder(2,local_index)=2
                  IF(ipartinitial(iwin).EQ.-1)THEN
                    iorder(1,local_index)=8
                    iorder(8,local_index)=7
                  ELSE
                    iorder(1,local_index)=7
                    iorder(8,local_index)=8
                  ENDIF
                ENDIF
                
                idplocal(iorder(6,local_index))=idz0(2,i)    ! Z final
                idplocal(iorder(7,local_index))=-idz0(2,i)
                idplocal(iorder(3,local_index))=idwp(1,i)    ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                idplocal(iorder(5,local_index))=5    ! b

                IF(abs(idout(1)).EQ.5)THEN
                  idplocal(iorder(1,local_index))=-5
                  IF(ipartinitial(iwin).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
                    idplocal(iorder(8,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idwm(1,i)
                    idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                  ENDIF
                ELSE
                  idplocal(iorder(2,local_index))=-5
                  IF(ipartinitial(iwin).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwm(1,i)+1)
                    idplocal(iorder(8,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwm(1,i)
                    idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                  ENDIF
                ENDIF

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1       
                    ENDIF
                  ENDIF 
                ENDIF
              ENDIF  ! end case 2
            ENDIF  ! end ZW+ in final state
c ZW- in final state 
            IF(nwmf(i).EQ.1 .AND. nz0f(i).EQ.1)THEN
c Case 1 tbbar or tbarb --> W+Z(bbbar) in the final state         
              ifound=0
              DO j=1,nz0(i)
                IF(iz0_b(j,i).EQ.1 .AND. iz0inout (j,i).EQ.2) ifound=1
              ENDDO
              IF(ifound.EQ.1)THEN
                IF(itypeinitial(1).EQ.0 .OR. idout(1).EQ.idout(2))THEN
*                WRITE(*,*) i,' t(0W) b tbar -> Z W-  : 26'
                  local_index=26
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  iorder(5,local_index)=4    ! Z  bbar
                  iorder(6,local_index)=3 
                  iorder(3,local_index)=7    ! W-
                  iorder(4,local_index)=8 
                  IF(ipartinitial(1).EQ.-1)THEN  ! Z
                    iorder(1,local_index)=2
                    iorder(7,local_index)=1
                  ELSE
                    iorder(1,local_index)=1
                    iorder(7,local_index)=2
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN  ! W+
                    iorder(2,local_index)=6
                    iorder(8,local_index)=5
                  ELSE
                    iorder(2,local_index)=5
                    iorder(8,local_index)=6
                  ENDIF

                  idplocal(iorder(5,local_index))=-idz0(2,i) ! Z  bbar
                  idplocal(iorder(6,local_index))=idz0(2,i) 
                  idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(4,local_index))=-(idwm(1,i)+1) 
                  IF(ipartinitial(1).EQ.-1)THEN ! Z
                    idplocal(iorder(1,local_index))=-idz0(1,i)
                    idplocal(iorder(7,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idz0(1,i)
                    idplocal(iorder(7,local_index))=-idz0(1,i)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN ! W+
                    idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
                    idplocal(iorder(8,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idwp(1,i)
                    idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                  ENDIF

                  IF(n_kphs.eq.1)THEN
                    iabase=local_index
                    igr(2)=1
                    DO j=1,8
                      igo(j,2)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(2).EQ.0)THEN
                      igr(2)=1
                      DO j=1,4
                        igo(2*j-1,2)=2*index(j,ipp)-1
                        igo(2*j,2)=2*index(j,ipa)
                      ENDDO
                      idwm_ref=idwm(1,i)
                      idwp_ref=idwp(1,i)
                    ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                      IF(idwm(1,i).NE.idwm_ref .OR. 
     &                                     idwp(1,i).NE.idwp_ref)THEN
                        DO j=1,4
                          igo(2*j-1,3)=2*index(j,ipp)-1
                          igo(2*j,3)=2*index(j,ipa)
                        ENDDO
                        igr(3)=1       
                      ENDIF
                    ENDIF 
                  ENDIF
                ENDIF
                IF(itypeinitial(2).EQ.0 .OR. idout(1).EQ.idout(2))
     &                                             THEN   
*                  WRITE(*,*) i,' t(W0) b tbar -> Z W-  : 25'
                  local_index=25
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  iorder(5,local_index)=4    ! Z  bbar
                  iorder(6,local_index)=3 
                  iorder(3,local_index)=7    ! W-
                  iorder(4,local_index)=8 
                  IF(ipartinitial(1).EQ.-1)THEN  ! W+
                    iorder(1,local_index)=6
                    iorder(7,local_index)=5
                  ELSE
                    iorder(1,local_index)=5
                    iorder(7,local_index)=6
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN  ! Z
                    iorder(2,local_index)=2
                    iorder(8,local_index)=1
                  ELSE
                    iorder(2,local_index)=1
                    iorder(8,local_index)=2
                  ENDIF

                  idplocal(iorder(5,local_index))=-idz0(2,i) ! Z  bbar
                  idplocal(iorder(6,local_index))=idz0(2,i) 
                  idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(4,local_index))=-(idwm(1,i)+1) 
                  IF(ipartinitial(1).EQ.-1)THEN ! W+
                    idplocal(iorder(1,local_index))=-(idwp(1,i)-1)
                    idplocal(iorder(7,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwp(1,i)
                    idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                  ENDIF
                  IF(ipartinitial(2).EQ.-1)THEN ! Z
                    idplocal(iorder(2,local_index))=-idz0(1,i)
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idz0(1,i)
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF

                  IF(n_kphs.eq.1)THEN
                    iabase=local_index
                    igr(2)=1
                    DO j=1,8
                      igo(j,2)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(2).EQ.0)THEN
                      igr(2)=1
                      DO j=1,4
                        igo(2*j-1,2)=2*index(j,ipp)-1
                        igo(2*j,2)=2*index(j,ipa)
                      ENDDO
                      idwm_ref=idwm(1,i)
                      idwp_ref=idwp(1,i)
                    ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                      IF(idwm(1,i).NE.idwm_ref .OR. 
     &                                     idwp(1,i).NE.idwp_ref)THEN
                        DO j=1,4
                          igo(2*j-1,3)=2*index(j,ipp)-1
                          igo(2*j,3)=2*index(j,ipa)
                        ENDDO
                        igr(3)=1       
                      ENDIF
                    ENDIF 
                  ENDIF
                ENDIF  
              ENDIF  ! end case 1
c Case 2  b T(W) --> t Z --> b W Z
              IF((nantib-nb).EQ.1)THEN  ! need bbar in final
                                                         ! state
*               WRITE(*,*) i,' t(W-) bbar -> Z tbar -> Z W- bbar  : 33'
                local_index=33
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(6,local_index)=3    ! Z final
                iorder(7,local_index)=4
                iorder(3,local_index)=7    ! W-
                iorder(4,local_index)=8
                iorder(5,local_index)=2    ! bbar

                DO j=1,2
                  IF(itypeinitial(j).EQ.1) iwin=j
                ENDDO
              
                IF(abs(idout(1)).EQ.5)THEN
                  iorder(1,local_index)=1
                  IF(ipartinitial(iwin).EQ.-1)THEN
                    iorder(2,local_index)=6
                    iorder(8,local_index)=5
                  ELSE
                    iorder(2,local_index)=5
                    iorder(8,local_index)=6
                  ENDIF
                ELSE
                  iorder(2,local_index)=1
                  IF(ipartinitial(iwin).EQ.-1)THEN
                    iorder(1,local_index)=6
                    iorder(8,local_index)=5
                  ELSE
                    iorder(1,local_index)=5
                    iorder(8,local_index)=6
                  ENDIF
                ENDIF
                

                idplocal(iorder(6,local_index))=idz0(2,i)    ! Z final
                idplocal(iorder(7,local_index))=-idz0(2,i)
                idplocal(iorder(3,local_index))=idwm(1,i)    ! W-
                idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(5,local_index))=-5    ! bbar

                IF(abs(idout(1)).EQ.5)THEN
                  idplocal(iorder(1,local_index))=5
                  IF(ipartinitial(iwin).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
                    idplocal(iorder(8,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idwp(1,i)
                    idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ELSE
                  idplocal(iorder(2,local_index))=5
                  IF(ipartinitial(iwin).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwp(1,i)-1)
                    idplocal(iorder(8,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwp(1,i)
                    idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ENDIF
                
                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1       
                    ENDIF
                  ENDIF 
                ENDIF
              ENDIF  ! end case 2
            ENDIF  ! end ZW- in final state
c WW in final state 
            IF(nwmf(i).EQ.1 .AND. nwpf(i).EQ.1 .AND. nz0(i).EQ.2)THEN
              IF(nantib.GE.1)THEN
*               WRITE(*,*) i,'  t(0 ) bbar -> W+ tbar(W- bbar) : 35'
                local_index=35
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(6,local_index)=5    ! W+
                iorder(7,local_index)=6 
                iorder(3,local_index)=7    ! W-
                iorder(4,local_index)=8
c locate Z with initial state bbar. If there are two identical Z's 
c symmetrization will take care of it.
                IF(idz0(1,i).NE.idz0(2,i))THEN  ! only one Z is b-bbar
                  IF(idz0(1,i).EQ.5)THEN
                    iz2top=1
                  ELSE
                    iz2top=2
                  ENDIF
                ELSE            ! Both Z's are b-bbar
                  IF(ipartinitial(1).EQ.ipartinitial(2))THEN ! bbar-bbar
                               ! symmetrization will take care of things
                    iz2top=1
                  ELSE
                    IF(ipartinitial(1).EQ.1)THEN ! the initial bbar 
                                                 ! is part-1
                      iz2top=iflagswap+1
                    ELSE        ! the initial bbar is part-2
                      iz2top=2-iflagswap
                    ENDIF
                  ENDIF
                ENDIF

c If iswap(i)=i first Z includes first incoming particle
                IF(iz2top.eq.1)THEN
                  iorder(iswap(1),local_index)=1      
                  iorder(5,local_index)=2
                  IF(ipartinitial(iswap(2)).EQ.-1)THEN
                    iorder(iswap(2),local_index)=4
                    iorder(8,local_index)=3
                  ELSE
                    iorder(iswap(2),local_index)=3
                    iorder(8,local_index)=4
                  ENDIF
                ELSE
                  iorder(iswap(2),local_index)=3      
                  iorder(5,local_index)=4
                  IF(ipartinitial(iswap(1)).EQ.-1)THEN
                    iorder(iswap(1),local_index)=2
                    iorder(8,local_index)=1
                  ELSE
                    iorder(iswap(1),local_index)=1
                    iorder(8,local_index)=2
                  ENDIF
                ENDIF
                
                idplocal(iorder(6,local_index))=idwp(1,i) ! W+
                idplocal(iorder(7,local_index))=-(idwp(1,i)-1) 
                idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                IF(iz2top.eq.1)THEN
                  idplocal(iorder(iswap(1),local_index))=5      
                  idplocal(iorder(5,local_index))=-5
                  IF(ipartinitial(iswap(2)).EQ.-1)THEN
                    idplocal(iorder(iswap(2),local_index))=-idz0(2,i)
                    idplocal(iorder(8,local_index))=idz0(2,i)
                  ELSE
                    idplocal(iorder(iswap(2),local_index))=idz0(2,i)
                    idplocal(iorder(8,local_index))=-idz0(2,i)
                  ENDIF
                ELSE
                  idplocal(iorder(iswap(2),local_index))=5      
                  idplocal(iorder(5,local_index))=-5
                  IF(ipartinitial(iswap(1)).EQ.-1)THEN
                    idplocal(iorder(iswap(1),local_index))=-idz0(1,i)
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(iswap(1),local_index))=idz0(1,i)
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF
                ENDIF


                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1       
                    ENDIF
                  ENDIF 
                ENDIF
              ENDIF
              IF(nb.GE.1)THEN
*               WRITE(*,*) i,'  t(0 ) b -> W- t(W+ b)  : 34'
                local_index=34
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(6,local_index)=7    ! W-
                iorder(7,local_index)=8 
                iorder(3,local_index)=5    ! W+
                iorder(4,local_index)=6
c locate Z with initial state b. If there are two identical Z's 
c symmetrization will take care of it.
                IF(idz0(1,i).NE.idz0(2,i))THEN  ! only one Z is b-bbar
                  IF(idz0(1,i).EQ.5)THEN
                    iz2top=1
                  ELSE
                    iz2top=2
                  ENDIF
                ELSE                           ! Both Z's are b-bbar
                  IF(ipartinitial(1).EQ.ipartinitial(2))THEN  ! b-b
                              ! symmetrization will take care of things
                    iz2top=1
                  ELSE
                    IF(ipartinitial(1).EQ.-1)THEN ! the initial b 
                                                  ! is part-1
                      iz2top=iflagswap+1
                    ELSE                     ! the initial b is part-2
                      iz2top=2-iflagswap
                    ENDIF
                  ENDIF
                ENDIF

                IF(iz2top.eq.1)THEN
                  iorder(iswap(1),local_index)=2      
                  iorder(5,local_index)=1
                  IF(ipartinitial(iswap(2)).EQ.-1)THEN
                    iorder(iswap(2),local_index)=4
                    iorder(8,local_index)=3
                  ELSE
                    iorder(iswap(2),local_index)=3
                    iorder(8,local_index)=4
                  ENDIF
                ELSE
                  iorder(iswap(2),local_index)=4      
                  iorder(5,local_index)=3
                  IF(ipartinitial(iswap(1)).EQ.-1)THEN
                    iorder(iswap(1),local_index)=2
                    iorder(8,local_index)=1
                  ELSE
                    iorder(iswap(1),local_index)=1
                    iorder(8,local_index)=2
                  ENDIF
                ENDIF


                idplocal(iorder(6,local_index))=idwm(1,i)    ! W-
                idplocal(iorder(7,local_index))=-(idwm(1,i)+1) 
                idplocal(iorder(3,local_index))=idwp(1,i)    ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                IF(iz2top.eq.1)THEN
                  idplocal(iorder(iswap(1),local_index))=-5      
                  idplocal(iorder(5,local_index))=5
                  IF(ipartinitial(iswap(2)).EQ.-1)THEN
                    idplocal(iorder(iswap(2),local_index))=-idz0(2,i)
                    idplocal(iorder(8,local_index))=idz0(2,i)
                  ELSE
                    idplocal(iorder(iswap(2),local_index))=idz0(2,i)
                    idplocal(iorder(8,local_index))=-idz0(2,i)
                  ENDIF
                ELSE
                  idplocal(iorder(iswap(2),local_index))=-5      
                  idplocal(iorder(5,local_index))=5
                  IF(ipartinitial(iswap(1)).EQ.-1)THEN
                    idplocal(iorder(iswap(1),local_index))=-idz0(1,i)
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(iswap(1),local_index))=idz0(1,i)
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF
                ENDIF

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(2)=1
                  DO j=1,8
                    igo(j,2)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                  idwm_ref=idwm(1,i)
                  idwp_ref=idwp(1,i)
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart,
     &                       idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(2).EQ.0)THEN
                    igr(2)=1
                    DO j=1,4
                      igo(2*j-1,2)=2*index(j,ipp)-1
                      igo(2*j,2)=2*index(j,ipa)
                    ENDDO
                    idwm_ref=idwm(1,i)
                    idwp_ref=idwp(1,i)
                  ELSEIF(igr(2).EQ.1 .AND. igr(3).EQ.0)THEN
                    IF(idwm(1,i).NE.idwm_ref .OR. idwp(1,i).NE.idwp_ref)
     &                                     THEN
                      DO j=1,4
                        igo(2*j-1,3)=2*index(j,ipp)-1
                        igo(2*j,3)=2*index(j,ipa)
                      ENDDO
                      igr(3)=1       
                    ENDIF
                  ENDIF 
                ENDIF
              ENDIF
            ENDIF  ! end WW in final state
          ENDIF
        ENDDO
      ELSEIF(ngluon.EQ.2 .AND. ngf.EQ.0)THEN
        DO i=1,nmatch
          IF(nz0f(i).EQ.1)THEN
*            WRITE(*,*) i,' gg -> Z  W+ W-  : 18'
            local_index=18
            ialfa(local_index)=1
            n_kphs=n_kphs+1

            iorder(1,local_index)=1 
            iorder(2,local_index)=2
            iorder(3,local_index)=5 ! W+
            iorder(4,local_index)=6  
            iorder(5,local_index)=7 ! W-
            iorder(6,local_index)=8
            iorder(7,local_index)=3 ! Z
            iorder(8,local_index)=4
          
            idplocal(iorder(1,local_index))=21
            idplocal(iorder(2,local_index))=21
            idplocal(iorder(7,local_index))=idz0(1,i)
            idplocal(iorder(8,local_index))=-idz0(1,i)
            idplocal(iorder(3,local_index))=idwp(1,i) ! W+
            idplocal(iorder(4,local_index))=-(idwp(1,i)-1)  
            idplocal(iorder(5,local_index))=idwm(1,i) ! W-
            idplocal(iorder(6,local_index))=-(idwm(1,i)+1)

            IF(n_kphs.eq.1)THEN
              iabase=local_index
              igr(5)=1
              DO j=1,8
                igo(j,5)=j
              ENDDO
              DO j=1,8
                idp(j)=idplocal(j)
              ENDDO
              DO j=1,2
                idp_inout(iorder(j,iabase))=-1
              ENDDO
            ELSE
              ial=local_index
              CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                   idplocal,idp,idp_inout,ipp,ipa)
              IF(igr(5).EQ.0)THEN
                igr(5)=1
                DO j=1,4
                  igo(2*j-1,5)=2*index(j,ipp)-1
                  igo(2*j,5)=2*index(j,ipa)
                ENDDO
              ENDIF 
            ENDIF
          ELSE   ! nz0f(i).EQ.3
c Z1 and Z2 make a Higgs if Z3 is a bbbar pair. The latter must be 
c put in position 7-8 of the phase space 2_4.
c No other phase space is necessary because of the symmetrization 
c performed in fxn.
*           WRITE(*,*) i,' gg -> Z Z Z  : 20'
            local_index=20
            ialfa(local_index)=1
            n_kphs=n_kphs+1
c Find Z(b-bbar) [if there are more than one, take the first 
c encountered] and put it at position 7-8 of the phase space 2_4; the 
c remaining two Z's are those which make a Higgs and their ID 
c information is kept in idz0_aux(2).
            j_aux=1
            k_aux=1
            DO WHILE(idz0(k_aux,i).NE.5 .AND. k_aux.LT.3)
              k_aux=k_aux+1
            ENDDO
c now k_aux represents the index of Z(b-bbar) [the first one, in 
c  presence of more than one b-bbar pairs]. If there are no b-bbar 
c pairs, then k_aux=3 .
            DO i_aux=1,k_aux-1
              idz0_aux(j_aux)=idz0(i_aux,i)
              j_aux=j_aux+1
            ENDDO
            DO i_aux=k_aux+1,3
              idz0_aux(j_aux)=idz0(i_aux,i)
              j_aux=j_aux+1
            ENDDO
c
            iorder(1,local_index)=1 
            iorder(2,local_index)=2
            iorder(3,local_index)=3 ! Z1
            iorder(4,local_index)=4  
            iorder(5,local_index)=5 ! Z2
            iorder(6,local_index)=6
            iorder(7,local_index)=7 ! Z3
            iorder(8,local_index)=8

            idplocal(iorder(1,local_index))=21
            idplocal(iorder(2,local_index))=21
            idplocal(iorder(3,local_index))=idz0_aux(1) ! Z1
            idplocal(iorder(4,local_index))=-idz0_aux(1)
            idplocal(iorder(5,local_index))=idz0_aux(2) ! Z2
            idplocal(iorder(6,local_index))=-idz0_aux(2)
            idplocal(iorder(7,local_index))=idz0(k_aux,i) ! Z3
            idplocal(iorder(8,local_index))=-idz0(k_aux,i)

            IF(n_kphs.eq.1)THEN
              iabase=local_index
              igr(6)=1
              DO j=1,8
                igo(j,6)=j
              ENDDO
              DO j=1,8
                idp(j)=idplocal(j)
              ENDDO
              DO j=1,2
                idp_inout(iorder(j,iabase))=-1
              ENDDO
            ELSE
              ial=local_index
              CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &             idplocal,idp,idp_inout,ipp,ipa)
              IF(igr(6).EQ.0)THEN
                igr(6)=1
                DO j=1,4
                  igo(2*j-1,6)=2*index(j,ipp)-1
                  igo(2*j,6)=2*index(j,ipa)
                ENDDO
              ENDIF 
            ENDIF
          ENDIF                 ! IF(nz0(i).EQ.1) ELSE
c LOOK for top 
          IF(nwpf(i).EQ.1 .AND. nb.GE.1)THEN
*           WRITE(*,*) i,' gg  -> t  t  : 56'
            local_index=56
            ialfa(local_index)=1
            n_kphs=n_kphs+1

            iorder(1,local_index)=1 
            iorder(2,local_index)=2
            iorder(7,local_index)=3    ! Z  b
            iorder(8,local_index)=4    !    bbar
            iorder(3,local_index)=5    ! W+
            iorder(4,local_index)=6  
            iorder(5,local_index)=7    ! W-
            iorder(6,local_index)=8

            idplocal(iorder(1,local_index))=21
            idplocal(iorder(2,local_index))=21
            idplocal(iorder(7,local_index))=idz0(1,i)    ! Z  b
            idplocal(iorder(8,local_index))= -idz0(1,i)   !    bbar
            idplocal(iorder(3,local_index))=idwp(1,i)    ! W+
            idplocal(iorder(4,local_index))=-(idwp(1,i)-1) 
            idplocal(iorder(5,local_index))=idwm(1,i)    ! W-
            idplocal(iorder(6,local_index))=-(idwm(1,i)+1)

            IF(n_kphs.eq.1)THEN
              iabase=local_index
              igr(5)=1
              DO j=1,8
                igo(j,5)=j
              ENDDO
              DO j=1,8
                idp(j)=idplocal(j)
              ENDDO
              DO j=1,2
                idp_inout(iorder(j,iabase))=-1
              ENDDO
            ELSE
              ial=local_index
              CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                   idplocal,idp,idp_inout,ipp,ipa)
              IF(igr(5).EQ.0)THEN
                igr(5)=1
                DO j=1,4
                  igo(2*j-1,5)=2*index(j,ipp)-1
                  igo(2*j,5)=2*index(j,ipa)
                ENDDO
              ENDIF 
            ENDIF
          ENDIF
        ENDDO    ! i=1,nmatch




      ELSEIF(ngluon.EQ.2 .AND. ngf.EQ.1)THEN

c This check prevents from DIS-like processes, i.e. e p -> X, which are
c not covered by Phantom. In this case the program stops with a warning.
        IF( (abs(iproc(1)).GE.11 .AND. abs(iproc(1)).LE. 16) .OR.
     &       (abs(iproc(2)).GE.11 .AND. abs(iproc(2)).LE. 16) )THEN
          WRITE(*,*) '***proc.f ERROR: this version of the code does '
          WRITE(*,*) '                 not cover processes of type '
          WRITE(*,*) '                 e g -> X'
          PRINT*,'iproc:',(iproc(j),j=1,8)
          STOP
        ENDIF

***********************************************************************
*AVAILABLE PHASE SPACE CHANNELS:
* t(g0) W+ W-                                     ialfa(8)  1_1_4
* t(0g) W+ W-                                     ialfa(9)  1_1_4
* t(gW) Z W-   and  t(gW) Z W+                    ialfa(10) 1_1_4
* t(Wg) Z W-   and  t(Wg) Z W+                    ialfa(11) 1_1_4
* t(g0) Z Z                                       ialfa(12) 1_1_4
* t(0g) Z Z                                       ialfa(13) 1_1_4
* t(gW) t                                         ialfa(30) 1_1_31
* t(Wg) t                                         ialfa(31) 1_1_31
* t(g ) W- t    and  t(g ) W+ t                   ialfa(36) 1_2_3
* g t(W ) Z {gW-}   and   g t(W ) Z {gW+}         ialfa(37) 1_2_3
* g t(W ) {gZ} W-   and   g t(W ) {gZ} W+         ialfa(38) 1_2_3
* g t(0 ) W+ {gW-}                                ialfa(39) 1_2_3
* g t(0 ) {gW+} W-                                ialfa(40) 1_2_3
* g t(0 ) {gZ} Z                                  ialfa(41) 1_2_3
* g t(0 ) Z {gZ}                                  ialfa(42) 1_2_3
* gq -> {gW+} t    and  gq -> {gW-} t             ialfa(60) 3_3
* g t(W ) {gtbar} -> g t(W ) g W-  
*            and   g t(W ) {gt} -> g t(W ) g W+   ialfa(75) 1_5to1_4to31
* g t(W ) g tbar -> g t(W ) {gW-}  
*            and   g t(W ) g t -> g t(W ) {gW+}   ialfa(76) 1_5to1_4to31
***********************************************************************
        DO i=1,nmatch

c Check for type of t-channel boson.
c ipartinitial(i)=0 means that incoming particle #i is a GLUON. In this 
c case itypeinitial(i) is conventionally set to 21
          IF(ipartinitial(1).EQ.-1)THEN ! 1=antiparticle #1, particle #2
                                        ! must be a gluon
            itypeinitial(1)=itype(1,i)
            itypeinitial(2)=21
          ELSEIF(ipartinitial(1).EQ.1)THEN ! 1=particle #1
            itypeinitial(1)=itype(invindex(1,iperm(i)),i) 
            itypeinitial(2)=21
          ELSE                  ! 1=gluon #1
            itypeinitial(1)=21
            IF(ipartinitial(2).EQ.-1)THEN ! 2=antiparticle #2
              itypeinitial(2)=itype(1,i)  
            ELSE                ! 2=particle #2
              itypeinitial(2)=itype(invindex(2,iperm(i)),i)
            ENDIF
          ENDIF
c end Check for type of t-channel boson


c nfinal(i).EQ.2 for all processes of type qg -> gX
          IF(nwpf(i).EQ.1 .AND. nwmf(i).EQ.1)THEN
***********************************************************************
* t(g0) W+ W-                                     ialfa(8)  1_1_4
* t(0g) W+ W-                                     ialfa(9)  1_1_4
* g t(0 ) W+ {gW-}                                ialfa(39) 1_2_3
* g t(0 ) {gW+} W-                                ialfa(40) 1_2_3
***********************************************************************

            IF(iproc(1).EQ.21)THEN
*              WRITE(*,*) i,' t(g0) W+ W-  : 8'
              local_index=8
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=1 ! g
              iorder(7,local_index)=2 ! g
              iorder(3,local_index)=5 ! W+
              iorder(4,local_index)=6
              iorder(5,local_index)=7 ! W-
              iorder(6,local_index)=8
              IF(ipartinitial(2).EQ.-1)THEN
                iorder(2,local_index)=4 ! Z
                iorder(8,local_index)=3
              ELSE
                iorder(2,local_index)=3 ! Z
                iorder(8,local_index)=4
              ENDIF

              idplocal(iorder(1,local_index))=21 ! g
              idplocal(iorder(7,local_index))=21 ! g
              idplocal(iorder(3,local_index))=idwp(1,i) ! W+
              idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
              idplocal(iorder(5,local_index))=idwm(1,i) ! W-
              idplocal(iorder(6,local_index))=-(idwm(1,i)+1)
              IF(ipartinitial(2).EQ.-1)THEN
                idplocal(iorder(2,local_index))=-idz0(1,i) ! Z
                idplocal(iorder(8,local_index))=idz0(1,i)
              ELSE
                idplocal(iorder(2,local_index))=idz0(1,i) ! Z
                idplocal(iorder(8,local_index))=-idz0(1,i)
              ENDIF

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(5)=1
                DO j=1,8
                  igo(j,5)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &               idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(5).EQ.0)THEN
                  igr(5)=1
                  DO j=1,4
                    igo(2*j-1,5)=2*index(j,ipp)-1
                    igo(2*j,5)=2*index(j,ipa)
                  ENDDO
                ENDIF 
              ENDIF


            ELSE                ! iproc(2).EQ.21
*              WRITE(*,*) i,' t(0g) W+ W-  : 9'
              local_index=9
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(2,local_index)=1 ! g
              iorder(8,local_index)=2 ! g
              iorder(3,local_index)=5 ! W+
              iorder(4,local_index)=6
              iorder(5,local_index)=7 ! W-
              iorder(6,local_index)=8
              IF(ipartinitial(1).EQ.-1)THEN
                iorder(1,local_index)=4 ! Z
                iorder(7,local_index)=3
              ELSE
                iorder(1,local_index)=3 ! Z
                iorder(7,local_index)=4
              ENDIF

              idplocal(iorder(2,local_index))=21 ! g
              idplocal(iorder(8,local_index))=21 ! g
              idplocal(iorder(3,local_index))=idwp(1,i) ! W+
              idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
              idplocal(iorder(5,local_index))=idwm(1,i) ! W-
              idplocal(iorder(6,local_index))=-(idwm(1,i)+1)
              IF(ipartinitial(1).EQ.-1)THEN
                idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                idplocal(iorder(7,local_index))=idz0(1,i)
              ELSE
                idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                idplocal(iorder(7,local_index))=-idz0(1,i)
              ENDIF

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(5)=1
                DO j=1,8
                  igo(j,5)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &               idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(5).EQ.0)THEN
                  igr(5)=1
                  DO j=1,4
                    igo(2*j-1,5)=2*index(j,ipp)-1
                    igo(2*j,5)=2*index(j,ipa)
                  ENDDO
                ENDIF 
              ENDIF
            ENDIF               ! IF(iproc(1).EQ.21)THEN

c Channels 'g t(0 ) {gW+} W-' and 'g t(0 ) W+ {gW-}' are ignored in case
c of two leptonic bosons in the final state (no set of q-qbar-g can 
c resonate at W mass, hence channels 't(g0) W+ W-' and 't(0g) W+ W-' are
c enough)
            IF((nwpflep(i)+nwmflep(i)).LT.2)THEN

              IF(nwpflep(i).EQ.0)THEN ! W- may be leptonic or not
*                WRITE(*,*) i,'  g t(0 ) {gW+} W- : 40'
                local_index=40
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(iproc(1).EQ.21)THEN
                  iorder(1,local_index)=1 ! g
                  IF(ipartinitial(2).EQ.-1)THEN
                    iorder(2,local_index)=4 ! Z
                    iorder(8,local_index)=3
                  ELSE
                    iorder(2,local_index)=3 ! Z
                    iorder(8,local_index)=4
                  ENDIF
                ELSE            ! iproc(2).EQ.21
                  iorder(2,local_index)=1 ! g
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=4 ! Z
                    iorder(8,local_index)=3
                  ELSE
                    iorder(1,local_index)=3 ! Z
                    iorder(8,local_index)=4
                  ENDIF
                ENDIF
                iorder(3,local_index)=5 ! W+
                iorder(4,local_index)=6
                iorder(5,local_index)=2 ! g
                iorder(6,local_index)=7 ! W-
                iorder(7,local_index)=8

                IF(iproc(1).EQ.21)THEN
                  idplocal(iorder(1,local_index))=21 ! g
                  IF(ipartinitial(2).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF
                ELSE            ! iproc(2).EQ.21
                  idplocal(iorder(2,local_index))=21 ! g
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF
                ENDIF
                idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                idplocal(iorder(5,local_index))=21 ! g
                idplocal(iorder(6,local_index))=idwm(1,i) ! W-
                idplocal(iorder(7,local_index))=-(idwm(1,i)+1)

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF
              ENDIF             ! IF(nwpflep(i).EQ.0)THEN


              IF(nwmflep(i).EQ.0)THEN ! W+ may be leptonic or not
*                WRITE(*,*) i,'  g t(0 ) W+ {gW-} : 39'
                local_index=39
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(iproc(1).EQ.21)THEN
                  iorder(1,local_index)=1 ! g
                  IF(ipartinitial(2).EQ.-1)THEN
                    iorder(2,local_index)=4 ! Z
                    iorder(8,local_index)=3
                  ELSE
                    iorder(2,local_index)=3 ! Z
                    iorder(8,local_index)=4
                  ENDIF
                ELSE            ! iproc(2).EQ.21
                  iorder(2,local_index)=1 ! g
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=4 ! Z
                    iorder(8,local_index)=3
                  ELSE
                    iorder(1,local_index)=3 ! Z
                    iorder(8,local_index)=4
                  ENDIF
                ENDIF
                iorder(3,local_index)=7 ! W-
                iorder(4,local_index)=8
                iorder(5,local_index)=2 ! g
                iorder(6,local_index)=5 ! W+
                iorder(7,local_index)=6

                IF(iproc(1).EQ.21)THEN
                  idplocal(iorder(1,local_index))=21 ! g
                  IF(ipartinitial(2).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF
                ELSE            ! iproc(2).EQ.21
                  idplocal(iorder(2,local_index))=21 ! g
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF
                ENDIF
                idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(5,local_index))=21 ! g
                idplocal(iorder(6,local_index))=idwp(1,i) ! W+
                idplocal(iorder(7,local_index))=-(idwp(1,i)-1)

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF
              ENDIF             ! IF(nwmflep(i).EQ.0)THEN

            ENDIF               ! IF((nwpflep(i)+nwmflep(i)).LT.2)THEN


          ELSEIF((nwmf(i)+nwpf(i)).EQ.1 .AND. nz0f(i).EQ.1)THEN
***********************************************************************
* t(gW) Z W-   and  t(gW) Z W+                    ialfa(10) 1_1_4
* t(Wg) Z W-   and  t(Wg) Z W+                    ialfa(11) 1_1_4
* g t(W ) Z {gW-}   and   g t(W ) Z {gW+}         ialfa(37) 1_2_3
* g t(W ) {gZ} W-   and   g t(W ) {gZ} W+         ialfa(38) 1_2_3
***********************************************************************
c Merging cases with one final W+ or W-. Requires determination of 
c final-state W's sign. 
c Notice that initial state W line is assumed to be partonic.

            IF(iproc(1).EQ.21)THEN
*             WRITE(*,*) i,' t(gW) Z W : 10'
              local_index=10
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=1 ! g
              iorder(7,local_index)=2 ! g
              iorder(3,local_index)=3 ! Z
              iorder(4,local_index)=4
              IF(nwmf(i).EQ.1)THEN
                iorder(5,local_index)=7 ! W-
                iorder(6,local_index)=8
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=6 ! W+
                  iorder(8,local_index)=5
                ELSE
                  iorder(2,local_index)=5 ! W+
                  iorder(8,local_index)=6
                ENDIF
              ELSE              ! nwpf(i).EQ.1
                iorder(5,local_index)=5 ! W+
                iorder(6,local_index)=6
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=8 ! W-
                  iorder(8,local_index)=7
                ELSE
                  iorder(2,local_index)=7 ! W-
                  iorder(8,local_index)=8
                ENDIF
              ENDIF

              idplocal(iorder(1,local_index))=21 ! g
              idplocal(iorder(7,local_index))=21 ! g
              idplocal(iorder(3,local_index))=idz0(1,i) ! Z
              idplocal(iorder(4,local_index))=-idz0(1,i)
              IF(nwmf(i).EQ.1)THEN
                idplocal(iorder(5,local_index))=idwm(1,i) ! W-
                idplocal(iorder(6,local_index))=-(idwm(1,i)+1)
                IF(ipartinitial(2).EQ.-1)THEN
                  idplocal(iorder(2,local_index))=-(idwp(1,i)-1) ! W+
                  idplocal(iorder(8,local_index))=idwp(1,i)
                ELSE
                  idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                ENDIF
              ELSE
                idplocal(iorder(5,local_index))=idwp(1,i)! W+
                idplocal(iorder(6,local_index))=-(idwp(1,i)-1)
                IF(ipartinitial(2).EQ.-1)THEN
                  idplocal(iorder(2,local_index))=-(idwm(1,i)+1) ! W-
                  idplocal(iorder(8,local_index))=idwm(1,i)
                ELSE
                  idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                ENDIF
              ENDIF

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(5)=1
                DO j=1,8
                  igo(j,5)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &               idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(5).EQ.0)THEN
                  igr(5)=1
                  DO j=1,4
                    igo(2*j-1,5)=2*index(j,ipp)-1
                    igo(2*j,5)=2*index(j,ipa)
                  ENDDO
                ENDIF
              ENDIF


            ELSE                ! iproc(2).EQ.21
*             WRITE(*,*) i,' t(Wg) Z W : 11'
              local_index=11
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(2,local_index)=1 ! g
              iorder(8,local_index)=2 ! g
              iorder(3,local_index)=3 ! Z
              iorder(4,local_index)=4
              IF(nwmf(i).EQ.1)THEN
                iorder(5,local_index)=7 ! W-
                iorder(6,local_index)=8
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=6 ! W+
                  iorder(7,local_index)=5
                ELSE
                  iorder(1,local_index)=5 ! W+
                  iorder(7,local_index)=6
                ENDIF
              ELSE              ! nwpf(i).EQ.1
                iorder(5,local_index)=5 ! W+
                iorder(6,local_index)=6
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=8 ! W-
                  iorder(7,local_index)=7
                ELSE
                  iorder(1,local_index)=7 ! W-
                  iorder(7,local_index)=8
                ENDIF
              ENDIF

              idplocal(iorder(2,local_index))=21 ! g
              idplocal(iorder(8,local_index))=21 ! g
              idplocal(iorder(3,local_index))=idz0(1,i) ! Z
              idplocal(iorder(4,local_index))=-idz0(1,i)
              IF(nwmf(i).EQ.1)THEN
                idplocal(iorder(5,local_index))=idwm(1,i) ! W-
                idplocal(iorder(6,local_index))=-(idwm(1,i)+1)
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
                  idplocal(iorder(7,local_index))=idwp(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                ENDIF
              ELSE              ! nwpf(i).EQ.1
                idplocal(iorder(5,local_index))=idwp(1,i) ! W+
                idplocal(iorder(6,local_index))=-(idwp(1,i)-1)
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
                  idplocal(iorder(7,local_index))=idwm(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                ENDIF
              ENDIF

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(5)=1
                DO j=1,8
                  igo(j,5)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &               idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(5).EQ.0)THEN
                  igr(5)=1
                  DO j=1,4
                    igo(2*j-1,5)=2*index(j,ipp)-1
                    igo(2*j,5)=2*index(j,ipa)
                  ENDDO
                ENDIF
              ENDIF

            ENDIF               ! IF(iproc(1).EQ.21)THEN


c Channels 'g t(W ) Z {gW}' and 'g t(W ) {gZ} W' are ignored in case of 
c two leptonic bosons in the final state (no set of q-qbar-g can 
c resonate at V mass, hence channels 't(gW) Z W' and 't(Wg) Z W' are 
c enough)
            IF((nwpflep(i)+nwmflep(i)+nz0flep(i)).LT.2)THEN

              IF((nwpflep(i)+nwmflep(i)).EQ.0)THEN ! Z may be leptonic
*               WRITE(*,*) i,' g t(W ) Z {gW} : 37'
                local_index=37
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(nwmf(i).EQ.1)THEN

                  IF(iproc(1).EQ.21)THEN
                    iorder(1,local_index)=1 ! g
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=6 ! W+
                      iorder(8,local_index)=5
                    ELSE
                      iorder(2,local_index)=5 ! W+
                      iorder(8,local_index)=6
                    ENDIF
                  ELSE          ! iproc(2).EQ.21
                    iorder(2,local_index)=1 ! g
                    IF(ipartinitial(1).EQ.-1)THEN
                      iorder(1,local_index)=6 ! W+
                      iorder(8,local_index)=5
                    ELSE
                      iorder(1,local_index)=5 ! W+
                      iorder(8,local_index)=6
                    ENDIF
                  ENDIF
                  iorder(3,local_index)=7 ! W-
                  iorder(4,local_index)=8
                  iorder(5,local_index)=2 ! g
                  iorder(6,local_index)=3 ! Z
                  iorder(7,local_index)=4

                ELSE            ! nwpf(i).EQ.1

                  IF(iproc(1).EQ.21)THEN
                    iorder(1,local_index)=1 ! g
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=8 ! W-
                      iorder(8,local_index)=7
                    ELSE
                      iorder(2,local_index)=7 ! W-
                      iorder(8,local_index)=8
                    ENDIF
                  ELSE          ! iproc(2).EQ.21
                    iorder(2,local_index)=1 ! g
                    IF(ipartinitial(1).EQ.-1)THEN
                      iorder(1,local_index)=8 ! W-
                      iorder(8,local_index)=7
                    ELSE
                      iorder(1,local_index)=7 ! W-
                      iorder(8,local_index)=8
                    ENDIF
                  ENDIF         ! IF(iproc(1).EQ.21)THEN
                  iorder(3,local_index)=5 ! W+
                  iorder(4,local_index)=6
                  iorder(5,local_index)=2 ! g
                  iorder(6,local_index)=3 ! Z
                  iorder(7,local_index)=4

                ENDIF           ! IF(nwmf(i).EQ.1)THEN

                IF(nwmf(i).EQ.1)THEN
                  IF(iproc(1).EQ.21)THEN
                    idplocal(iorder(1,local_index))=21 ! g
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwp(1,i)-1) !W+
                      idplocal(iorder(8,local_index))=idwp(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                      idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                    ENDIF
                  ELSE          ! iproc(2).EQ.21
                    idplocal(iorder(2,local_index))=21 ! g
                    IF(ipartinitial(1).EQ.-1)THEN
                      idplocal(iorder(1,local_index))=-(idwp(1,i)-1) !W+
                      idplocal(iorder(8,local_index))=idwp(1,i)
                    ELSE
                      idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                      idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                    ENDIF
                  ENDIF
                  idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                  idplocal(iorder(5,local_index))=21 ! g
                  idplocal(iorder(6,local_index))=idz0(1,i) ! Z
                  idplocal(iorder(7,local_index))=-idz0(1,i)
                ELSE            ! nwpf(i).EQ.1

                  IF(iproc(1).EQ.21)THEN
                    idplocal(iorder(1,local_index))=21 ! g
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwm(1,i)+1) !W-
                      idplocal(iorder(8,local_index))=idwm(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                      idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                    ENDIF
                  ELSE          ! iproc(2).EQ.21
                    idplocal(iorder(2,local_index))=21 ! g
                    IF(ipartinitial(1).EQ.-1)THEN
                      idplocal(iorder(1,local_index))=-(idwm(1,i)+1) !W-
                      idplocal(iorder(8,local_index))=idwm(1,i)
                    ELSE
                      idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                      idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                    ENDIF
                  ENDIF         ! IF(iproc(1).EQ.21)THEN
                  idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                  idplocal(iorder(5,local_index))=21 ! g
                  idplocal(iorder(6,local_index))=idz0(1,i) ! Z
                  idplocal(iorder(7,local_index))=-idz0(1,i)

                ENDIF           ! IF(nwmf(i).EQ.1)THEN

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF((nwpflep(i)+nwmflep(i)).EQ.0)THEN


              IF(nz0flep(i).EQ.0)THEN ! final state W may be leptonic
*               WRITE(*,*) i,' g t(W ) {gZ} W : 38'
                local_index=38
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(nwmf(i).EQ.1)THEN

                  IF(iproc(1).EQ.21)THEN
                    iorder(1,local_index)=1 ! g
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=6 ! W+
                      iorder(8,local_index)=5
                    ELSE
                      iorder(2,local_index)=5 ! W+
                      iorder(8,local_index)=6
                    ENDIF
                  ELSE          ! iproc(2).EQ.21
                    iorder(2,local_index)=1 ! g
                    IF(ipartinitial(1).EQ.-1)THEN
                      iorder(1,local_index)=6 ! W+
                      iorder(8,local_index)=5
                    ELSE
                      iorder(1,local_index)=5 ! W+
                      iorder(8,local_index)=6
                    ENDIF
                  ENDIF 
                  iorder(3,local_index)=3 ! Z
                  iorder(4,local_index)=4
                  iorder(5,local_index)=2 ! g
                  iorder(6,local_index)=7 ! W-
                  iorder(7,local_index)=8

                ELSE            ! nwpf(i).EQ.1

                  IF(iproc(1).EQ.21)THEN
                    iorder(1,local_index)=1 ! g
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=8 ! W-
                      iorder(8,local_index)=7
                    ELSE
                      iorder(2,local_index)=7 ! W-
                      iorder(8,local_index)=8
                    ENDIF
                  ELSE          ! iproc(2).EQ.21
                    iorder(2,local_index)=1 ! g
                    IF(ipartinitial(1).EQ.-1)THEN
                      iorder(1,local_index)=8 ! W-
                      iorder(8,local_index)=7
                    ELSE
                      iorder(1,local_index)=7 ! W-
                      iorder(8,local_index)=8
                    ENDIF
                  ENDIF
                  iorder(3,local_index)=3 ! Z
                  iorder(4,local_index)=4
                  iorder(5,local_index)=2 ! g
                  iorder(6,local_index)=5 ! W+
                  iorder(7,local_index)=6

                ENDIF           ! IF(nwmf(i).EQ.1)THEN

                IF(nwmf(i).EQ.1)THEN

                  IF(iproc(1).EQ.21)THEN
                    idplocal(iorder(1,local_index))=21 ! g
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwp(1,i)-1) !W+
                      idplocal(iorder(8,local_index))=idwp(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                      idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                    ENDIF
                  ELSE          ! iproc(2).EQ.21
                    idplocal(iorder(2,local_index))=21 ! g
                    IF(ipartinitial(1).EQ.-1)THEN
                      idplocal(iorder(1,local_index))=-(idwp(1,i)-1) !W+
                      idplocal(iorder(8,local_index))=idwp(1,i)
                    ELSE
                      idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                      idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                    ENDIF
                  ENDIF
                  idplocal(iorder(3,local_index))=idz0(1,i) ! Z
                  idplocal(iorder(4,local_index))=-idz0(1,i)
                  idplocal(iorder(5,local_index))=21 ! g
                  idplocal(iorder(6,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(7,local_index))=-(idwm(1,i)+1)

                ELSE            ! nwpf(i).EQ.1

                  IF(iproc(1).EQ.21)THEN
                    idplocal(iorder(1,local_index))=21 ! g
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwm(1,i)+1) !W-
                      idplocal(iorder(8,local_index))=idwm(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                      idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                    ENDIF
                  ELSE          ! iproc(2).EQ.21
                    idplocal(iorder(2,local_index))=21 ! g
                    IF(ipartinitial(1).EQ.-1)THEN
                      idplocal(iorder(1,local_index))=-(idwm(1,i)+1) !W-
                      idplocal(iorder(8,local_index))=idwm(1,i)
                    ELSE
                      idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                      idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                    ENDIF
                  ENDIF
                  idplocal(iorder(3,local_index))=idz0(1,i) ! Z
                  idplocal(iorder(4,local_index))=-idz0(1,i)
                  idplocal(iorder(5,local_index))=21 ! g
                  idplocal(iorder(6,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(7,local_index))=-(idwp(1,i)-1)

                ENDIF           ! IF(nwmf(i).EQ.1)THEN

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF(nz0flep(i).EQ.0)THEN

            ENDIF      ! IF((nwmflep(i)+nwpflep(i)+nz0flep(i).LT.2)THEN


          ELSEIF(nz0f(i).EQ.2)THEN
***********************************************************************
* t(g0) Z Z                                       ialfa(12) 1_1_4
* t(0g) Z Z                                       ialfa(13) 1_1_4
* g t(0 ) {gZ} Z                                  ialfa(41) 1_2_3
* g t(0 ) Z {gZ}                                  ialfa(42) 1_2_3
***********************************************************************

            IF(iproc(1).EQ.21)THEN
*             WRITE(*,*) i,' t(g0) Z Z  : 12'
              local_index=12
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=1 ! g
              iorder(7,local_index)=2 ! g
              IF(ipartinitial(2).EQ.-1)THEN
                iorder(2,local_index)=4 ! Z
                iorder(8,local_index)=3
              ELSE
                iorder(2,local_index)=3 ! Z
                iorder(8,local_index)=4
              ENDIF
              iorder(3,local_index)=5 ! Z
              iorder(4,local_index)=6
              iorder(5,local_index)=7 ! Z
              iorder(6,local_index)=8

              idplocal(iorder(1,local_index))=21 ! g
              idplocal(iorder(7,local_index))=21 ! g
              IF(ipartinitial(2).EQ.-1)THEN
                idplocal(iorder(2,local_index))=-idz0(1,i) ! Z
                idplocal(iorder(8,local_index))=idz0(1,i)
              ELSE
                idplocal(iorder(2,local_index))=idz0(1,i) ! Z
                idplocal(iorder(8,local_index))=-idz0(1,i)
              ENDIF
              idplocal(iorder(3,local_index))=idz0(2,i) !Z 
              idplocal(iorder(4,local_index))=-idz0(2,i)
              idplocal(iorder(5,local_index))=idz0(3,i) ! Z
              idplocal(iorder(6,local_index))=-idz0(3,i)


              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(6)=1
                DO j=1,8
                  igo(j,6)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &               idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(6).EQ.0)THEN
                  igr(6)=1
                  DO j=1,4
                    igo(2*j-1,6)=2*index(j,ipp)-1
                    igo(2*j,6)=2*index(j,ipa)
                  ENDDO
                ENDIF
              ENDIF


            ELSE                !iproc(2).EQ.21
*              WRITE(*,*) i,' t(0g) Z Z  : 13'
              local_index=13
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(2,local_index)=1 ! g
              iorder(8,local_index)=2 ! g
              IF(ipartinitial(1).EQ.-1)THEN
                iorder(1,local_index)=4 ! Z
                iorder(7,local_index)=3
              ELSE
                iorder(1,local_index)=3 ! Z
                iorder(7,local_index)=4
              ENDIF
              iorder(3,local_index)=5 ! Z
              iorder(4,local_index)=6
              iorder(5,local_index)=7 ! Z
              iorder(6,local_index)=8

              idplocal(iorder(2,local_index))=21 ! g
              idplocal(iorder(8,local_index))=21 ! g
              IF(ipartinitial(1).EQ.-1)THEN
                idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                idplocal(iorder(7,local_index))=idz0(1,i)
              ELSE
                idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                idplocal(iorder(7,local_index))=-idz0(1,i)
              ENDIF
              idplocal(iorder(3,local_index))=idz0(2,i) ! Z
              idplocal(iorder(4,local_index))=-idz0(2,i)
              idplocal(iorder(5,local_index))=idz0(3,i) ! Z
              idplocal(iorder(6,local_index))=-idz0(3,i)

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(6)=1
                DO j=1,8
                  igo(j,6)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &               idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(6).EQ.0)THEN
                  igr(6)=1
                  DO j=1,4
                    igo(2*j-1,6)=2*index(j,ipp)-1
                    igo(2*j,6)=2*index(j,ipa)
                  ENDDO
                ENDIF
              ENDIF

            ENDIF               ! IF(iproc(1).EQ.21)THEN



c Channels 'g t(0 ) {gZ} Z' and 'g t(0 ) Z {gZ}' are ignored in case of 
c two leptonic bosons in the final state (no set of q-qbar-g can 
c resonate at Z mass, hence channels 't(g0) Z Z' and 't(0g) Z Z' are 
c enough).

            IF(nz0flep(i) .LT. 2)THEN
c Careful in assignment of iorders. Must determine which Z is leptonic.
c In case of only one leptonic Z (if any), I assume it is represented by
c idz0(3,i), whereas idz0(1,i) represents the t-channel Z boson.

*             WRITE(*,*) i,' g t(0 ) {gZ} Z  : 41'
              local_index=41
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              IF(iproc(1).EQ.21)THEN
                iorder(1,local_index)=1 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=4 ! Z
                  iorder(8,local_index)=3
                ELSE
                  iorder(2,local_index)=3 ! Z
                  iorder(8,local_index)=4
                ENDIF
              ELSE              ! iproc(2).EQ.21
                iorder(2,local_index)=1 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=4 ! Z
                  iorder(8,local_index)=3
                ELSE
                  iorder(1,local_index)=3 ! Z
                  iorder(8,local_index)=4
                ENDIF
              ENDIF
              iorder(3,local_index)=5 ! Z
              iorder(4,local_index)=6
              iorder(5,local_index)=2 ! g
              iorder(6,local_index)=7 ! Z leptonic [if any]
              iorder(7,local_index)=8

              IF(iproc(1).EQ.21)THEN
                idplocal(iorder(1,local_index))=21 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  idplocal(iorder(2,local_index))=-idz0(1,i) ! Z
                  idplocal(iorder(8,local_index))=idz0(1,i)
                ELSE
                  idplocal(iorder(2,local_index))=idz0(1,i) ! Z
                  idplocal(iorder(8,local_index))=-idz0(1,i)
                ENDIF
              ELSE              ! iproc(2).EQ.21
                idplocal(iorder(2,local_index))=21 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                  idplocal(iorder(8,local_index))=idz0(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                  idplocal(iorder(8,local_index))=-idz0(1,i)
                ENDIF
              ENDIF
              idplocal(iorder(3,local_index))=idz0(2,i) ! Z
              idplocal(iorder(4,local_index))=-idz0(2,i)
              idplocal(iorder(5,local_index))=21 ! g
              idplocal(iorder(6,local_index))=idz0(3,i) !Z lept [if any]
              idplocal(iorder(7,local_index))=-idz0(3,i)

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(6)=1
                DO j=1,8
                  igo(j,6)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &               idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(6).EQ.0)THEN
                  igr(6)=1
                  DO j=1,4
                    igo(2*j-1,6)=2*index(j,ipp)-1
                    igo(2*j,6)=2*index(j,ipa)
                  ENDDO
                ENDIF
              ENDIF


              IF(nz0flep(i).EQ.0)THEN
                IF(idz0(2,i).NE.idz0(3,i))THEN
c Careful in assignment of iorders. idz0(3,i) has to be associated with 
c {gZ}.

*                 WRITE(*,*) i,' g t(0 ) Z {gZ}  : 42'
                  local_index=42
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  iorder(1,local_index)=iorder(1,local_index-1) 
                  iorder(2,local_index)=iorder(2,local_index-1)
                  iorder(3,local_index)=iorder(6,local_index-1)
                  iorder(4,local_index)=iorder(7,local_index-1)
                  iorder(5,local_index)=iorder(5,local_index-1)  
                  iorder(6,local_index)=iorder(3,local_index-1)  
                  iorder(7,local_index)=iorder(4,local_index-1)  
                  iorder(8,local_index)=iorder(8,local_index-1)
                ENDIF           ! IF(idz0(2,i).NE.idz0(3,i))THEN

              ENDIF             ! IF(nz0flep(i).EQ.0)THEN

            ENDIF               ! IF(nz0flep(i) .LT. 2)THEN


          ELSE
            print*,'***proc.f : ERROR in gq -> g X !!!'
            print*,'   Check boson matching'
            STOP
          ENDIF                 ! IF(nwpf(i).EQ.1 .AND. nwmf(i).EQ.1)

c Check for top
***********************************************************************
* t(gW) t                                         ialfa(30) 1_1_31
* t(Wg) t                                         ialfa(31) 1_1_31
* t(g ) W- t    and  t(g ) W+ t                   ialfa(36) 1_2_3
* gq -> {gW+} t    and  gq -> {gW-} t             ialfa(60) 3_3
* gq -> W+ {gt} -> W+ g W-
*            and   gq -> W- {gt} -> W- g W+       ialfa(65) 2_4to31
* gq -> W+ g t -> W+ {gW-}
*            and   gq -> W- g t -> W- {gW+}       ialfa(66) 2_4to31
* g t(W ) {gtbar} -> g t(W ) g W-  
*            and   g t(W ) {gt} -> g t(W ) g W+   ialfa(75) 1_5to1_4to31
* g t(W ) g tbar -> g t(W ) {gW-}  
*            and   g t(W ) g t -> g t(W ) {gW+}   ialfa(76) 1_5to1_4to31
***********************************************************************
          ifound=0
          DO j=1,nz0(i)
            IF(iz0_b(j,i).EQ.1) ifound=1
          ENDDO
          IF(ifound.EQ.1)THEN
            IF(nwmf(i).EQ.1 .AND. nwpf(i).EQ.1)THEN
***********************************************************************
* t(g ) W- t    and  t(g ) W+ tbar                ialfa(36) 1_2_3
* g bbar -> {gW+} tbar   and   g b -> {gW-} t     ialfa(60) 3_3
* g bbar -> W+ {gtbar} -> W+ g W-
*            and   g b -> W- {gt} -> W- g W+      ialfa(65) 2_4to31
* g bbar -> W+ g tbar -> W+ {gW-}
*            and   g b -> W- g t -> W- {gW+}      ialfa(66) 2_4to31
***********************************************************************
              IF((nantib-nb).EQ.-1)THEN ! g b -> b X
*               WRITE(*,*) i,' t(g ) W- t  : 36'
                local_index=36
                ialfa(local_index)=1
                n_kphs=n_kphs+1

c giuseppe 05/01/2007
c                IF(iproc(1).EQ.1)THEN
                IF(iproc(1).EQ.21)THEN
c end giuseppe 05/01/2007
                  iorder(1,local_index)=1 ! g
                  iorder(2,local_index)=4 ! bbar
                ELSE            ! iproc(2).EQ.21
                  iorder(1,local_index)=4 ! bbar
                  iorder(2,local_index)=1 ! g
                ENDIF
                iorder(3,local_index)=5 ! W+
                iorder(4,local_index)=6 
                iorder(5,local_index)=3 ! b
                iorder(6,local_index)=7 ! W-
                iorder(7,local_index)=8 
                iorder(8,local_index)=2 ! g

c giuseppe 05/01/2007
c                IF(iproc(1).EQ.1)THEN
                IF(iproc(1).EQ.21)THEN
c end giuseppe 05/01/2007
                  idplocal(iorder(1,local_index))=21 ! g
                  idplocal(iorder(2,local_index))=-idz0(1,i) ! bbar
                ELSE            ! iproc(2).EQ.21
                  idplocal(iorder(1,local_index))=-idz0(1,i) ! bbar
                  idplocal(iorder(2,local_index))=21 ! g
                ENDIF
                idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1) 
                idplocal(iorder(5,local_index))=idz0(1,i) ! b
                idplocal(iorder(6,local_index))=idwm(1,i) ! W-
                idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(8,local_index))=21 ! g

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF


*               WRITE(*,*) i,' g b -> W- {gt} -> W- g W+  : 65'
                local_index=65
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(iproc(1).EQ.21)THEN
                  iorder(1,local_index)=1 ! g
                  iorder(2,local_index)=4 ! bbar
                ELSE            ! iproc(2).EQ.21
                  iorder(1,local_index)=4 ! bbar
                  iorder(2,local_index)=1 ! g
                ENDIF
                iorder(3,local_index)=5 ! W+
                iorder(4,local_index)=6
                iorder(5,local_index)=3 ! b
                iorder(6,local_index)=2 ! g
                iorder(7,local_index)=7 ! W-
                iorder(8,local_index)=8

                IF(iproc(1).EQ.21)THEN
                  idplocal(iorder(1,local_index))=21 ! g
                  idplocal(iorder(2,local_index))=-idz0(1,i) ! bbar
                ELSE            ! iproc(2).EQ.21
                  idplocal(iorder(1,local_index))=-idz0(1,i) ! bbar
                  idplocal(iorder(2,local_index))=21 ! g
                ENDIF
                idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                idplocal(iorder(5,local_index))=idz0(1,i) ! b
                idplocal(iorder(6,local_index))=21 ! g
                idplocal(iorder(7,local_index))=idwm(1,i) ! W-
                idplocal(iorder(8,local_index))=-(idwm(1,i)+1)

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF


                IF(nwmflep(i).EQ.0)THEN
*                 WRITE(*,*) i,' g b -> {gW-} t  : 60'
                  local_index=60
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  IF(iproc(1).EQ.21)THEN
                    iorder(1,local_index)=1 ! g
                    iorder(2,local_index)=4 ! bbar
                  ELSE          ! iproc(2).EQ.21
                    iorder(1,local_index)=4 ! bbar
                    iorder(2,local_index)=1 ! g
                  ENDIF
                  iorder(3,local_index)=7 ! W-
                  iorder(4,local_index)=8
                  iorder(7,local_index)=2 ! g
                  iorder(5,local_index)=5 ! W+
                  iorder(6,local_index)=6
                  iorder(8,local_index)=3 ! b

                  IF(iproc(1).EQ.21)THEN
                    idplocal(iorder(1,local_index))=21 ! g
                    idplocal(iorder(2,local_index))=-idz0(1,i) ! bbar
                  ELSE          ! iproc(2).EQ.21
                    idplocal(iorder(1,local_index))=-idz0(1,i) ! bbar
                    idplocal(iorder(2,local_index))=21 ! g
                  ENDIF
                  idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                  idplocal(iorder(7,local_index))=21 ! g
                  idplocal(iorder(5,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(6,local_index))=-(idwp(1,i)-1)
                  idplocal(iorder(8,local_index))=idz0(1,i) ! b

                  IF(n_kphs.eq.1)THEN
                    iabase=local_index
                    igr(5)=1
                    DO j=1,8
                      igo(j,5)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                   idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(5).EQ.0)THEN
                      igr(5)=1
                      DO j=1,4
                        igo(2*j-1,5)=2*index(j,ipp)-1
                        igo(2*j,5)=2*index(j,ipa)
                      ENDDO
                    ENDIF
                  ENDIF

                ENDIF           ! IF(nwmflep(i).EQ.0)THEN


                IF(nwpflep(i).EQ.0)THEN
*                 WRITE(*,*) i,' g b -> W- g t -> W- {gW+}  : 66'
                  local_index=66
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  IF(iproc(1).EQ.21)THEN
                    iorder(1,local_index)=1 ! g
                    iorder(2,local_index)=4 ! bbar
                  ELSE          ! iproc(2).EQ.21
                    iorder(1,local_index)=4 ! bbar
                    iorder(2,local_index)=1 ! g
                  ENDIF
                  iorder(3,local_index)=5 ! W+
                  iorder(4,local_index)=6
                  iorder(5,local_index)=2 ! g
                  iorder(6,local_index)=3 ! b
                  iorder(7,local_index)=7 ! W-
                  iorder(8,local_index)=8

                  IF(iproc(1).EQ.21)THEN
                    idplocal(iorder(1,local_index))=21 ! g
                    idplocal(iorder(2,local_index))=-idz0(1,i) ! bbar
                  ELSE          ! iproc(2).EQ.21
                    idplocal(iorder(1,local_index))=-idz0(1,i) ! bbar
                    idplocal(iorder(2,local_index))=21 ! g
                  ENDIF
                  idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                  idplocal(iorder(5,local_index))=21 ! g
                  idplocal(iorder(6,local_index))=idz0(1,i) ! b
                  idplocal(iorder(7,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(8,local_index))=-(idwm(1,i)+1)

                  IF(n_kphs.eq.1)THEN
                    iabase=local_index
                    igr(5)=1
                    DO j=1,8
                      igo(j,5)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                   idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(5).EQ.0)THEN
                      igr(5)=1
                      DO j=1,4
                        igo(2*j-1,5)=2*index(j,ipp)-1
                        igo(2*j,5)=2*index(j,ipa)
                      ENDDO
                    ENDIF
                  ENDIF

                ENDIF           ! IF(nwpflep(i).EQ.0)THEN


              ELSEIF((nantib-nb).EQ.1)THEN ! g bbar -> bbar X
*               WRITE(*,*) i,' t(g ) W+ tbar  : 36'
                local_index=36
                ialfa(local_index)=1
                n_kphs=n_kphs+1

c giuseppe 05/01/2007
c                IF(iproc(1).EQ.1)THEN
                IF(iproc(1).EQ.21)THEN
c end giuseppe 05/01/2007
                  iorder(1,local_index)=1 ! g
                  iorder(2,local_index)=3 ! b
                ELSE            ! iproc(2).EQ.21
                  iorder(1,local_index)=3 ! b
                  iorder(2,local_index)=1 ! g
                ENDIF
                iorder(3,local_index)=7 ! W-
                iorder(4,local_index)=8
                iorder(5,local_index)=4 ! bbar
                iorder(6,local_index)=5 ! W+
                iorder(7,local_index)=6
                iorder(8,local_index)=2 ! g

c giuseppe 05/01/2007
c                IF(iproc(1).EQ.1)THEN
                IF(iproc(1).EQ.21)THEN
c end giuseppe 05/01/2007
                  idplocal(iorder(1,local_index))=21 ! g
                  idplocal(iorder(2,local_index))=idz0(1,i) ! b
                ELSE            ! iproc(2).EQ.21
                  idplocal(iorder(1,local_index))=idz0(1,i) ! b
                  idplocal(iorder(2,local_index))=21 ! g
                ENDIF
                idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(5,local_index))=-idz0(1,i) ! bbar
                idplocal(iorder(6,local_index))=idwp(1,i) ! W+
                idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                idplocal(iorder(8,local_index))=21 ! g

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF


*               WRITE(*,*) i,' g bbar -> W+ {gtbar} -> W+ g W-  : 65'
                local_index=65
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(iproc(1).EQ.21)THEN
                  iorder(1,local_index)=1 ! g
                  iorder(2,local_index)=3 ! b
                ELSE            ! iproc(2).EQ.21
                  iorder(1,local_index)=3 ! b
                  iorder(2,local_index)=1 ! g
                ENDIF
                iorder(3,local_index)=7 ! W-
                iorder(4,local_index)=8
                iorder(5,local_index)=4 ! bbar
                iorder(6,local_index)=2 ! g
                iorder(7,local_index)=5 ! W+
                iorder(8,local_index)=6

                IF(iproc(1).EQ.21)THEN
                  idplocal(iorder(1,local_index))=21 ! g
                  idplocal(iorder(2,local_index))=idz0(1,i) ! b
                ELSE            ! iproc(2).EQ.21
                  idplocal(iorder(1,local_index))=idz0(1,i) ! b
                  idplocal(iorder(2,local_index))=21 ! g
                ENDIF
                idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(5,local_index))=-idz0(1,i) ! bbar
                idplocal(iorder(6,local_index))=21 ! g
                idplocal(iorder(7,local_index))=idwp(1,i) ! W+
                idplocal(iorder(8,local_index))=-(idwp(1,i)-1)

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF


                IF(nwpflep(i).EQ.0)THEN
*                 WRITE(*,*) i,' g bbar -> {gW+} tbar  : 60'
                  local_index=60
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  IF(iproc(1).EQ.21)THEN
                    iorder(1,local_index)=1 ! g
                    iorder(2,local_index)=3 ! b
                  ELSE          ! iproc(2).EQ.21
                    iorder(1,local_index)=3 ! b
                    iorder(2,local_index)=1 ! g
                  ENDIF
                  iorder(3,local_index)=5 ! W+
                  iorder(4,local_index)=6
                  iorder(7,local_index)=2 ! g
                  iorder(5,local_index)=7 ! W-
                  iorder(6,local_index)=8
                  iorder(8,local_index)=4 ! bbar

                  IF(iproc(1).EQ.21)THEN
                    idplocal(iorder(1,local_index))=21 ! g
                    idplocal(iorder(2,local_index))=idz0(1,i) ! b
                  ELSE          ! iproc(2).EQ.21
                    idplocal(iorder(1,local_index))=idz0(1,i) ! b
                    idplocal(iorder(2,local_index))=21 ! g
                  ENDIF
                  idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                  idplocal(iorder(7,local_index))=21 ! g
                  idplocal(iorder(5,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(6,local_index))=-(idwm(1,i)+1)
                  idplocal(iorder(8,local_index))=-idz0(1,i) ! bbar

                  IF(n_kphs.eq.1)THEN
                    iabase=local_index
                    igr(5)=1
                    DO j=1,8
                      igo(j,5)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                   idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(5).EQ.0)THEN
                      igr(5)=1
                      DO j=1,4
                        igo(2*j-1,5)=2*index(j,ipp)-1
                        igo(2*j,5)=2*index(j,ipa)
                      ENDDO
                    ENDIF
                  ENDIF

                ENDIF           ! IF(nwpflep(i).EQ.0)THEN


                IF(nwmflep(i).EQ.0)THEN
*                 WRITE(*,*) i,' g bbar -> W+ g tbar -> W+ {gW-}  : 66'
                  local_index=66
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  IF(iproc(1).EQ.21)THEN
                    iorder(1,local_index)=1 ! g
                    iorder(2,local_index)=3 ! b
                  ELSE          ! iproc(2).EQ.21
                    iorder(1,local_index)=3 ! b
                    iorder(2,local_index)=1 ! g
                  ENDIF
                  iorder(3,local_index)=7 ! W-
                  iorder(4,local_index)=8
                  iorder(5,local_index)=2 ! g
                  iorder(6,local_index)=4 ! bbar
                  iorder(7,local_index)=5 ! W+
                  iorder(8,local_index)=6

                  IF(iproc(1).EQ.21)THEN
                    idplocal(iorder(1,local_index))=21 ! g
                    idplocal(iorder(2,local_index))=idz0(1,i) ! b
                  ELSE          ! iproc(2).EQ.21
                    idplocal(iorder(1,local_index))=idz0(1,i) ! b
                    idplocal(iorder(2,local_index))=21 ! g
                  ENDIF
                  idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                  idplocal(iorder(5,local_index))=21 ! g
                  idplocal(iorder(6,local_index))=-idz0(1,i) ! bbar
                  idplocal(iorder(7,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(8,local_index))=-(idwp(1,i)-1)

                  IF(n_kphs.eq.1)THEN
                    iabase=local_index
                    igr(5)=1
                    DO j=1,8
                      igo(j,5)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                   idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(5).EQ.0)THEN
                      igr(5)=1
                      DO j=1,4
                        igo(2*j-1,5)=2*index(j,ipp)-1
                        igo(2*j,5)=2*index(j,ipa)
                      ENDDO
                    ENDIF
                  ENDIF

                ENDIF           ! IF(nwmflep(i).EQ.0)THEN

              ENDIF             ! IF((nantib-nb).EQ.-1)THEN


            ELSEIF(nwmf(i).EQ.1 .AND. nz0f(i).EQ.1)THEN
***********************************************************************
* t(gW) t                                         ialfa(30) 1_1_31
* t(Wg) t                                         ialfa(31) 1_1_31
* g t(W ) {gtbar} -> g t(W ) g W-  
*            and   g t(W ) {gt} -> g t(W ) g W+   ialfa(75) 1_5to1_4to31
* g t(W ) g tbar -> g t(W ) {gW-}  
*            and   g t(W ) g t -> g t(W ) {gW+}   ialfa(76) 1_5to1_4to31
***********************************************************************
c Notice that here Z line is assumed to be b-bar. Only final state W- 
c can be leptonic.

              IF(iproc(1).EQ.21)THEN
*               WRITE(*,*) i,' t(gW) tbar  : 30'
                local_index=30
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=1 ! g
                iorder(7,local_index)=2 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=6 ! W+
                  iorder(8,local_index)=5
                ELSE
                  iorder(2,local_index)=5 ! W+
                  iorder(8,local_index)=6
                ENDIF
                iorder(3,local_index)=7 ! W-
                iorder(4,local_index)=8
                iorder(5,local_index)=4 ! bbar
                iorder(6,local_index)=3 ! b

                idplocal(iorder(1,local_index))=21 ! g
                idplocal(iorder(7,local_index))=21 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  idplocal(iorder(2,local_index))=-(idwp(1,i)-1) ! W+
                  idplocal(iorder(8,local_index))=idwp(1,i)
                ELSE
                  idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                ENDIF
                idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(5,local_index))=-idz0(1,i) ! bbar
                idplocal(iorder(6,local_index))=idz0(1,i) ! b

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF


              ELSE              ! iproc(2).EQ.21
*               WRITE(*,*) i,' t(Wg) tbar  : 31'
                local_index=31
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(2,local_index)=1 ! g
                iorder(8,local_index)=2 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=6 ! W+
                  iorder(7,local_index)=5
                ELSE
                  iorder(1,local_index)=5 ! W+
                  iorder(7,local_index)=6
                ENDIF
                iorder(3,local_index)=7 ! W-
                iorder(4,local_index)=8
                iorder(5,local_index)=4 ! bbar
                iorder(6,local_index)=3 ! b

                idplocal(iorder(2,local_index))=21 ! g
                idplocal(iorder(8,local_index))=21 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
                  idplocal(iorder(7,local_index))=idwp(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                ENDIF
                idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(5,local_index))=-idz0(1,i) ! bbar
                idplocal(iorder(6,local_index))=idz0(1,i) ! b

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF(iproc(1).EQ.21)THEN


*         WRITE(*,*) i,' g t(W ) {gtbar} -> g t(W ) g W-  : 75'
              local_index=75
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              IF(iproc(1).EQ.21)THEN
                iorder(1,local_index)=1 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=6 ! W+
                  iorder(8,local_index)=5
                ELSE
                  iorder(2,local_index)=5 ! W+
                  iorder(8,local_index)=6
                ENDIF
              ELSE              ! iproc(2).EQ.21
                iorder(2,local_index)=1 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=6 ! W+
                  iorder(8,local_index)=5
                ELSE
                  iorder(1,local_index)=5 ! W+
                  iorder(8,local_index)=6
                ENDIF
              ENDIF
              iorder(3,local_index)=7 ! W-
              iorder(4,local_index)=8
              iorder(5,local_index)=4 ! bbar
              iorder(6,local_index)=2 ! g
              iorder(7,local_index)=3 ! b

              IF(iproc(1).EQ.21)THEN
                idplocal(iorder(1,local_index))=21 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  idplocal(iorder(2,local_index))=-(idwp(1,i)-1) ! W+
                  idplocal(iorder(8,local_index))=idwp(1,i)
                ELSE
                  idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                ENDIF
              ELSE              ! iproc(2).EQ.21
                idplocal(iorder(2,local_index))=21 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
                  idplocal(iorder(8,local_index))=idwp(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                ENDIF
              ENDIF
              idplocal(iorder(3,local_index))=idwm(1,i) ! W-
              idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
              idplocal(iorder(5,local_index))=-idz0(1,i) ! bbar
              idplocal(iorder(6,local_index))=21 ! g
              idplocal(iorder(7,local_index))=idz0(1,i) ! b

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(5)=1
                DO j=1,8
                  igo(j,5)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &               idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(5).EQ.0)THEN
                  igr(5)=1
                  DO j=1,4
                    igo(2*j-1,5)=2*index(j,ipp)-1
                    igo(2*j,5)=2*index(j,ipa)
                  ENDDO
                ENDIF
              ENDIF


              IF(nwmflep(i).EQ.0)THEN
c This channel differs from the previous one by the simple exchange
c  g <-> bbar in the phase space 1_5to1_4to31.

*               WRITE(*,*) i,' g t(W ) g tbar -> g t(W ) {gW-}  : 76'
                local_index=76
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=iorder(1,local_index-1)
                iorder(2,local_index)=iorder(2,local_index-1)
                iorder(3,local_index)=iorder(3,local_index-1) ! W-
                iorder(4,local_index)=iorder(4,local_index-1)  
                iorder(5,local_index)=iorder(6,local_index-1) ! g
                iorder(6,local_index)=iorder(5,local_index-1) ! bbar
                iorder(7,local_index)=iorder(7,local_index-1) ! b
                iorder(8,local_index)=iorder(8,local_index-1)

              ENDIF             ! IF(nwmflep(i).EQ.0)THEN


            ELSEIF(nwpf(i).EQ.1 .AND. nz0f(i).EQ.1)THEN
***********************************************************************
* t(gW) t                                         ialfa(30) 1_1_31
* t(Wg) t                                         ialfa(31) 1_1_31
* g t(W ) {gtbar} -> g t(W ) g W-  
*            and   g t(W ) {gt} -> g t(W ) g W+   ialfa(75) 1_5to1_4to31
* g t(W ) g tbar -> g t(W ) {gW-}  
*            and   g t(W ) g t -> g t(W ) {gW+}   ialfa(76) 1_5to1_4to31
***********************************************************************
c Notice that here Z line is assumed to be b-bar. Only final state W+
c can be leptonic.

              IF(iproc(1).EQ.21)THEN
*               WRITE(*,*) i,' t(gW) t  : 30'
                local_index=30
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=1 ! g
                iorder(7,local_index)=2 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=8 ! W-
                  iorder(8,local_index)=7
                ELSE
                  iorder(2,local_index)=7 ! W-
                  iorder(8,local_index)=8
                ENDIF
                iorder(3,local_index)=5 ! W+
                iorder(4,local_index)=6
                iorder(5,local_index)=3 ! b
                iorder(6,local_index)=4 ! bbar

                idplocal(iorder(1,local_index))=21 ! g
                idplocal(iorder(7,local_index))=21 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  idplocal(iorder(2,local_index))=-(idwm(1,i)+1) ! W-
                  idplocal(iorder(8,local_index))=idwm(1,i)
                ELSE
                  idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                ENDIF
                idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                idplocal(iorder(5,local_index))=idz0(1,i) ! b
                idplocal(iorder(6,local_index))=-idz0(1,i) ! bbar

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF


              ELSE              ! iproc(2).EQ.21
*               WRITE(*,*) i,' t(Wg) t  : 31'
                local_index=31
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(2,local_index)=1 ! g
                iorder(8,local_index)=2 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=8 ! W-
                  iorder(7,local_index)=7
                ELSE
                  iorder(1,local_index)=7 ! W-
                  iorder(7,local_index)=8
                ENDIF
                iorder(3,local_index)=5 ! W+
                iorder(4,local_index)=6
                iorder(5,local_index)=3 ! b
                iorder(6,local_index)=4 ! bbar

                idplocal(iorder(2,local_index))=21 ! g
                idplocal(iorder(8,local_index))=21 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
                  idplocal(iorder(7,local_index))=idwm(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                ENDIF
                idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                idplocal(iorder(5,local_index))=idz0(1,i) ! b
                idplocal(iorder(6,local_index))=-idz0(1,i) ! bbar

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF(iproc(1).EQ.21)THEN


*             WRITE(*,*) i,' g t(W ) {gt} -> g t(W ) g W+  : 75'
              local_index=75
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              IF(iproc(1).EQ.21)THEN
                iorder(1,local_index)=1 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  iorder(2,local_index)=8 ! W-
                  iorder(8,local_index)=7
                ELSE
                  iorder(2,local_index)=7 ! W-
                  iorder(8,local_index)=8
                ENDIF
              ELSE              ! iproc(2).EQ.21
                iorder(2,local_index)=1 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=8 ! W-
                  iorder(8,local_index)=7
                ELSE
                  iorder(1,local_index)=7 ! W-
                  iorder(8,local_index)=8
                ENDIF
              ENDIF
              iorder(3,local_index)=5 ! W+
              iorder(4,local_index)=6
              iorder(5,local_index)=3 ! b
              iorder(6,local_index)=2 ! g
              iorder(7,local_index)=4 ! bbar

              IF(iproc(1).EQ.21)THEN
                idplocal(iorder(1,local_index))=21 ! g
                IF(ipartinitial(2).EQ.-1)THEN
                  idplocal(iorder(2,local_index))=-(idwm(1,i)+1) ! W-
                  idplocal(iorder(8,local_index))=idwm(1,i)
                ELSE
                  idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                ENDIF
              ELSE              ! iproc(2).EQ.21
                idplocal(iorder(2,local_index))=21 ! g
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
                  idplocal(iorder(8,local_index))=idwm(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                ENDIF
              ENDIF
              idplocal(iorder(3,local_index))=idwp(1,i) ! W+
              idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
              idplocal(iorder(5,local_index))=idz0(1,i) ! b
              idplocal(iorder(6,local_index))=21 ! g
              idplocal(iorder(7,local_index))=-idz0(1,i) ! bbar

              IF(n_kphs.eq.1)THEN
                iabase=local_index
                igr(5)=1
                DO j=1,8
                  igo(j,5)=j
                ENDDO
                DO j=1,8
                  idp(j)=idplocal(j)
                ENDDO
                DO j=1,2
                  idp_inout(iorder(j,iabase))=-1
                ENDDO
              ELSE
                ial=local_index
                CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &               idplocal,idp,idp_inout,ipp,ipa)
                IF(igr(5).EQ.0)THEN
                  igr(5)=1
                  DO j=1,4
                    igo(2*j-1,5)=2*index(j,ipp)-1
                    igo(2*j,5)=2*index(j,ipa)
                  ENDDO
                ENDIF
              ENDIF


              IF(nwpflep(i).EQ.0)THEN
c This channel differs from the previous one by the simple exchange
c  g <-> b in the phase space 1_5to1_4to31.

*               WRITE(*,*) i,' g t(W ) g t -> g t(W ) {gW+}  : 76'
                local_index=76
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=iorder(1,local_index-1)
                iorder(2,local_index)=iorder(2,local_index-1)
                iorder(3,local_index)=iorder(3,local_index-1) ! W+
                iorder(4,local_index)=iorder(4,local_index-1)  
                iorder(5,local_index)=iorder(6,local_index-1) ! g
                iorder(6,local_index)=iorder(5,local_index-1) ! b
                iorder(7,local_index)=iorder(7,local_index-1) ! bbar
                iorder(8,local_index)=iorder(8,local_index-1)

              ENDIF             ! IF(nwpflep(i).EQ.0)THEN

            ENDIF               !IF(nwmf(i).EQ.1 .AND. nwpf(i).EQ.1)THEN

          ENDIF                 ! IF(ifound.EQ.1)THEN
c end Check for top

        ENDDO                   ! DO i=1,nmatch




      ELSEIF(ngluon.EQ.2 .AND. ngf.EQ.2)THEN

c Check existence of at least one pair of final-state quarks
        DO i=1,nmatch
          IF( (abs(iproc(1)).GE.11 .AND. abs(iproc(1)).LE. 16) .AND.
     &         (abs(iproc(2)).GE.11 .AND. abs(iproc(2)).LE. 16) )THEN
            IF( (nwpflep(i)+nwmflep(i)+nz0flep(i)).EQ.nfinal(i) )THEN
              WRITE(*,*) '***proc.f ERROR: check iproc'
              PRINT*,'iproc:',(iproc(j),j=1,8)
              STOP
            ENDIF
          ENDIF
        ENDDO

***********************************************************************
*AVAILABLE PHASE SPACE CHANNELS:
* t(0W) {ggW-}   and   t(0W) {ggW+}               ialfa(27) 1_1_31 
* t(W0) {ggW-}   and   t(W0) {ggW+}               ialfa(28) 1_1_31
* t(WW) {ggZ}                                     ialfa(29) 1_1_31
* t(00) {ggZ}                                     ialfa(32) 1_1_31
* W- -> {gZ} {gW-}  and  W+ -> {gZ} {gW+}         ialfa(57) 3_3
* Z  -> {gW+} {gW-}                               ialfa(58) 3_3
* Z  -> {gZ} {gZ}                                 ialfa(59) 3_3
* W- -> {ggZ} W-     and  W+ -> {ggZ} W+          ialfa(61) 2_4to31
* W- -> Z {ggW-}    and  W+ -> Z {ggW+}           ialfa(62) 2_4to31
* Z  -> {ggW+} W-                                 ialfa(63) 2_4to31
* Z  -> W+ {ggW-}                                 ialfa(64) 2_4to31
* Z  -> {ggZ} Z                                   ialfa(67) 2_4to31
* Z  -> Z {ggZ}                                   ialfa(68) 2_4to31
* W- -> {ggtbar} -> gg W-   
*            and  W+ -> {ggt} -> gg W+            ialfa(69) 1_5to1_4to31
* W- -> g {gtbar} -> g {gW-}                           
*            and  W+ -> g {gt} -> g {gW+}         ialfa(70) 1_5to1_4to31
* W- -> gg tbar -> {ggW-}   
*            and  W+ -> gg t -> {ggW+}            ialfa(71) 1_5to1_4to31
* t(W ) {ggtbar} -> t(W ) gg W- 
*            and  t(W ) {ggt} -> t(W ) gg W+      ialfa(72) 1_5to1_4to31
* t(W ) g {gtbar} -> t(W ) g {gW-} 
*            and  t(W ) g {gt} -> t(W ) g {gW+}   ialfa(73) 1_5to1_4to31
* t(W ) gg tbar -> t(W ) {ggW-} 
*            and  t(W ) gg t -> t(W ) {ggW+}      ialfa(74) 1_5to1_4to31
* t(0 ) W- {ggq} -> t(0 ) W- g{gq}
*          and  t(0 ) W+ {ggq} -> t(0 ) W+ g{gq}  ialfa(43) 1_2_3
* t(0 ) W- {ggq} -> t(0 ) W- q{gg}
*          and  t(0 ) W+ {ggq} -> t(0 ) W+ q{gg}  ialfa(44) 1_2_3
* t(W ) Z  {ggq} -> t(W ) Z  g{gq}                ialfa(45) 1_2_3
* t(W ) Z  {ggq} -> t(W ) Z  q{gg}                ialfa(46) 1_2_3
* qq -> {gq} {gqW-}   and  qq -> {gq} {gqW+}      ialfa(23) 2_4
* t(W+ ) Z {ggq} -> t(W+ ) Z g{gq}                ialfa(47) 1_2_3
* t(W+ ) Z {ggq} -> t(W+ ) Z q{gg}                ialfa(48) 1_2_3
* t(W- ) Z {ggq} -> t(W- ) Z g{gq}                ialfa(49) 1_2_3
* t(W- ) Z {ggq} -> t(W- ) Z q{gg}                ialfa(50) 1_2_3
* t(0[1] ) Z {ggq} -> t(0[1] ) Z g{gq}            ialfa(51) 1_2_3
* t(0[1] ) Z {ggq} -> t(0[1] ) Z q{gg}            ialfa(52) 1_2_3
* t(0[2] ) Z {ggq} -> t(0[2] ) Z g{gq}            ialfa(53) 1_2_3
* t(0[2] ) Z {ggq} -> t(0[2] ) Z q{gg}            ialfa(54) 1_2_3
* qq -> {gq} {gqZ}                                ialfa(24) 2_4
***********************************************************************
        DO i=1,nmatch
c nfinal(i) can be 1 or 2 for processes of type qq -> gg X
          IF(nfinal(i).EQ.2)THEN
***********************************************************************
* W- -> {gZ} {gW-}  and  W+ -> {gZ} {gW+}         ialfa(57) 3_3
* Z  -> {gW+} {gW-}                               ialfa(58) 3_3
* Z  -> {gZ} {gZ}                                 ialfa(59) 3_3
* W- -> {ggZ} W-     and  W+ -> {ggZ} W+          ialfa(61) 2_4to31
* W- -> Z {ggW-}    and  W+ -> Z {ggW+}           ialfa(62) 2_4to31
* Z  -> {ggW+} W-                                 ialfa(63) 2_4to31
* Z  -> W+ {ggW-}                                 ialfa(64) 2_4to31
* Z  -> {ggZ} Z                                   ialfa(67) 2_4to31
* Z  -> Z {ggZ}                                   ialfa(68) 2_4to31
***********************************************************************
            IF(nwpf(i).EQ.1 .AND. nwmf(i).EQ.1)THEN
***********************************************************************
* Z  -> {gW+} {gW-}                               ialfa(58) 3_3
* Z  -> {ggW+} W-                                 ialfa(63) 2_4to31
* Z  -> W+ {ggW-}                                 ialfa(64) 2_4to31
***********************************************************************
c In case of two final state bosons decaying into leptons, just one of 
c the following three channels is necessary: we choose 
c 'Z  -> {gW+} {gW-}   [3_3]'. As no set of q-qbar-g can resonate at W
c mass, the value of invariant mass parameters rm1[a,b],rm2[a,b] in phsp
c 3_3 must be set to zero [see procextraini.f].

              IF((nwpflep(i).EQ.0 .AND. nwmflep(i).EQ.0) .OR. 
     &             (nwpflep(i).EQ.1 .AND. nwmflep(i).EQ.1))THEN
*               WRITE(*,*) i,' Z  -> {gW+} {gW-}  : 58'
                local_index=58
                ialfa(local_index)=1
                n_kphs=n_kphs+1

c information (common) for procextraini.f
                num_wpflep(local_index)=nwpflep(i)
                num_wmflep(local_index)=nwmflep(i)
                num_z0flep(local_index)=nz0flep(i)
c
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=4 ! Z
                  iorder(2,local_index)=3
                ELSE
                  iorder(1,local_index)=3 ! Z
                  iorder(2,local_index)=4
                ENDIF
                iorder(3,local_index)=5 ! W+
                iorder(4,local_index)=6
                iorder(7,local_index)=1 ! g
                iorder(5,local_index)=7 ! W-
                iorder(6,local_index)=8
                iorder(8,local_index)=2 ! g

                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                  idplocal(iorder(2,local_index))=idz0(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                  idplocal(iorder(2,local_index))=-idz0(1,i)
                ENDIF
                idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                idplocal(iorder(7,local_index))=21 ! g
                idplocal(iorder(5,local_index))=idwm(1,i) ! W-
                idplocal(iorder(6,local_index))=-(idwm(1,i)+1)
                idplocal(iorder(8,local_index))=21 ! g

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF((nwpflep(i).EQ.0 .AND. ...))


              IF(nwpflep(i).EQ.0)THEN
*               WRITE(*,*) i,' Z  -> {ggW+} W-  : 63'
                local_index=63
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=4 ! Z
                  iorder(2,local_index)=3
                ELSE
                  iorder(1,local_index)=3 ! Z
                  iorder(2,local_index)=4
                ENDIF
                iorder(3,local_index)=5 ! W+
                iorder(4,local_index)=6
                iorder(5,local_index)=1 ! g
                iorder(6,local_index)=2 ! g
                iorder(7,local_index)=7 ! W-
                iorder(8,local_index)=8

                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                  idplocal(iorder(2,local_index))=idz0(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                  idplocal(iorder(2,local_index))=-idz0(1,i)
                ENDIF
                idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                idplocal(iorder(5,local_index))=21 ! g
                idplocal(iorder(6,local_index))=21 ! g
                idplocal(iorder(7,local_index))=idwm(1,i) ! W-
                idplocal(iorder(8,local_index))=-(idwm(1,i)+1)

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF(nwpflep(i).EQ.0)THEN


              IF(nwmflep(i).EQ.0)THEN
c This channel differs from the previous one by the simple exchange 
c  W+ <-> W- in the phase space 2_4to31.

*               WRITE(*,*) i,' Z  -> W+ {ggW-}  : 64'
                local_index=64
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=iorder(1,local_index-1)
                iorder(2,local_index)=iorder(2,local_index-1)
                iorder(3,local_index)=iorder(7,local_index-1) ! W-
                iorder(4,local_index)=iorder(8,local_index-1)  
                iorder(5,local_index)=iorder(5,local_index-1) ! g
                iorder(6,local_index)=iorder(6,local_index-1) ! g
                iorder(7,local_index)=iorder(3,local_index-1) ! W+
                iorder(8,local_index)=iorder(4,local_index-1)

              ENDIF             ! IF(nwmflep(i).EQ.0)THEN


            ELSEIF((nwmf(i)+nwpf(i)).EQ.1 .AND. nz0f(i).EQ.1)THEN
***********************************************************************
* W- -> {gZ} {gW-}  and  W+ -> {gZ} {gW+}         ialfa(57) 3_3
* W- -> {ggZ} W-     and  W+ -> {ggZ} W+          ialfa(61) 2_4to31
* W- -> Z {ggW-}    and  W+ -> Z {ggW+}           ialfa(62) 2_4to31
***********************************************************************
c Merging cases with one final W+ or W-. Requires determination of 
c final-state W's sign

c In case of two final state bosons decaying into leptons, just one of 
c the following three channels is necessary: we choose 
c 'W  -> {gZ} {gW}   [3_3]'. As no set of q-qbar-g can resonate at W/Z
c mass, the value of invariant mass parameters rm1[a,b],rm2[a,b] in phsp
c 3_3 must be set to zero [see procextraini.f].

              IF((nz0flep(i).EQ.0 .AND. (nwpflep(i)+nwmflep(i)).EQ.0) 
     &             .OR. 
     &             (nz0flep(i).EQ.1 .AND. (nwpflep(i)+nwmflep(i)).EQ.1))
     &             THEN
*               WRITE(*,*) i,' W -> {gZ} {gW}  : 57'
                local_index=57
                ialfa(local_index)=1
                n_kphs=n_kphs+1

c information (common) for procextraini.f
                num_wpflep(local_index)=nwpflep(i)
                num_wmflep(local_index)=nwmflep(i)
                num_z0flep(local_index)=nz0flep(i)
c
                IF(nwmf(i).EQ.1)THEN
                  iorder(5,local_index)=7 ! W-
                  iorder(6,local_index)=8
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=6 ! W+
                    iorder(2,local_index)=5
                  ELSE
                    iorder(1,local_index)=5 ! W+
                    iorder(2,local_index)=6
                  ENDIF
                ELSE            ! nwpf(i).EQ.1
                  iorder(5,local_index)=5 ! W+
                  iorder(6,local_index)=6
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=8 ! W-
                    iorder(2,local_index)=7
                  ELSE
                    iorder(1,local_index)=7 ! W-
                    iorder(2,local_index)=8
                  ENDIF
                ENDIF
                iorder(3,local_index)=3 ! Z
                iorder(4,local_index)=4
                iorder(7,local_index)=1 ! g
                iorder(8,local_index)=2 ! g

                IF(nwmf(i).EQ.1)THEN
                  idplocal(iorder(5,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(6,local_index))=-(idwm(1,i)+1)
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
                    idplocal(iorder(2,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                    idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ELSE            ! nwpf(i).EQ.1
                  idplocal(iorder(5,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(6,local_index))=-(idwp(1,i)-1)
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
                    idplocal(iorder(2,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                    idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
                  ENDIF
                ENDIF
                idplocal(iorder(3,local_index))=idz0(1,i) ! Z
                idplocal(iorder(4,local_index))=-idz0(1,i)
                idplocal(iorder(7,local_index))=21 ! g
                idplocal(iorder(8,local_index))=21 ! g

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF((nz0flep(i).EQ.0 .AND. ... ))


              IF(nz0flep(i).EQ.0)THEN
*               WRITE(*,*) i,' W -> {ggZ} W  : 61'
                local_index=61
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(nwmf(i).EQ.1)THEN
                  iorder(7,local_index)=7 ! W-
                  iorder(8,local_index)=8
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=6 ! W+
                    iorder(2,local_index)=5
                  ELSE
                    iorder(1,local_index)=5 ! W+
                    iorder(2,local_index)=6
                  ENDIF
                ELSE            ! nwpf(i).EQ.1
                  iorder(7,local_index)=5 ! W+
                  iorder(8,local_index)=6
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=8 ! W-
                    iorder(2,local_index)=7
                  ELSE
                    iorder(1,local_index)=7 ! W-
                    iorder(2,local_index)=8
                  ENDIF
                ENDIF
                iorder(3,local_index)=3 ! Z
                iorder(4,local_index)=4
                iorder(5,local_index)=1 ! g
                iorder(6,local_index)=2 ! g

                IF(nwmf(i).EQ.1)THEN
                  idplocal(iorder(7,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
                    idplocal(iorder(2,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                    idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ELSE            ! nwpf(i).EQ.1
                  idplocal(iorder(7,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
                    idplocal(iorder(2,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                    idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
                  ENDIF
                ENDIF
                idplocal(iorder(3,local_index))=idz0(1,i) ! Z
                idplocal(iorder(4,local_index))=-idz0(1,i)
                idplocal(iorder(5,local_index))=21 ! g
                idplocal(iorder(6,local_index))=21 ! g

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF(nz0flep(i).EQ.0)THEN


              IF((nwpflep(i)+nwmflep(i)).EQ.0)THEN
*               WRITE(*,*) i,' W -> Z {ggW}  : 62'
                local_index=62
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                IF(nwmf(i).EQ.1)THEN
                  iorder(3,local_index)=7 ! W-
                  iorder(4,local_index)=8
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=6 ! W+
                    iorder(2,local_index)=5
                  ELSE
                    iorder(1,local_index)=5 ! W+
                    iorder(2,local_index)=6
                  ENDIF
                ELSE            ! nwpf(i).EQ.1
                  iorder(3,local_index)=5 ! W+
                  iorder(4,local_index)=6
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=8 ! W-
                    iorder(2,local_index)=7
                  ELSE
                    iorder(1,local_index)=7 ! W-
                    iorder(2,local_index)=8
                  ENDIF
                ENDIF
                iorder(7,local_index)=3 ! Z
                iorder(8,local_index)=4
                iorder(5,local_index)=1 ! g
                iorder(6,local_index)=2 ! g

                IF(nwmf(i).EQ.1)THEN
                  idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                  idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
                    idplocal(iorder(2,local_index))=idwp(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                    idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ELSE            ! nwpf(i).EQ.1
                  idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                  idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
                    idplocal(iorder(2,local_index))=idwm(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                    idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
                  ENDIF
                ENDIF
                idplocal(iorder(7,local_index))=idz0(1,i) ! Z
                idplocal(iorder(8,local_index))=-idz0(1,i)
                idplocal(iorder(5,local_index))=21 ! g
                idplocal(iorder(6,local_index))=21 ! g

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF((nwpflep(i)+nwmflep(i)).EQ.0)THEN


            ELSEIF(nz0f(i).EQ.2)THEN
***********************************************************************
* Z  -> {gZ} {gZ}                                 ialfa(59) 3_3
* Z  -> {ggZ} Z                                   ialfa(67) 2_4to31
* Z  -> Z {ggZ}                                   ialfa(68) 2_4to31
***********************************************************************
c Ezio 2018_08_02
c Added channel 6: 
* t(00) Z Z                                       ialfa(6)  1_1_4
c to account for the possibility of the two gluons radiating from the initial
c q-qbar state and the possibility of having H --> ZZ if the quarks are b's
c
c In case of two final state bosons decaying into leptons, just one of 
c the following three channels is necessary: we choose 
c 'Z  -> {gZ} {gZ}   [3_3]'. As no set of q-qbar-g can resonate at Z
c mass, the value of invariant mass parameters rm1[a,b],rm2[a,b] in phsp
c 3_3 must be set to zero [see procextraini.f].

*               WRITE(*,*) i,' t(00) Z Z  :6  1_1_4
c       if(npippo.EQ.0) GOTO 666
                local_index=6
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(3,local_index)=5    ! Z
                iorder(4,local_index)=6 
                iorder(5,local_index)=7    ! Z
                iorder(6,local_index)=8 
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=4
                  iorder(2,local_index)=3
                ELSE
                  iorder(1,local_index)=3
                  iorder(2,local_index)=4
                ENDIF
                iorder(7,local_index)=1   ! g
                iorder(8,local_index)=2
                  
                idplocal(iorder(3,local_index))=idz0(2,i)   ! Z
                idplocal(iorder(4,local_index))=-idz0(2,i) 
                idplocal(iorder(5,local_index))=idz0(3,i)    ! Z
                idplocal(iorder(6,local_index))=-idz0(3,i) 
                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-idz0(1,i)
                  idplocal(iorder(2,local_index))=idz0(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(1,i)
                  idplocal(iorder(2,local_index))=-idz0(1,i)
                ENDIF
                idplocal(iorder(7,local_index))=21 ! g
                idplocal(iorder(8,local_index))=21 ! g
                  
                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(6)=1
                  DO j=1,8
                    igo(j,6)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(6).EQ.0)THEN
                    igr(6)=1
                    DO j=1,4
                      igo(2*j-1,6)=2*index(j,ipp)-1
                      igo(2*j,6)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF
c666    CONTINUE
c end Ezio 2018_08_02


              IF(nz0flep(i).EQ.0 .OR. nz0flep(i).EQ.2)THEN
*               WRITE(*,*) i,' Z  -> {gZ} {gZ}  : 59'
                local_index=59
                ialfa(local_index)=1
                n_kphs=n_kphs+1

c information (common) for procextraini.f
                num_wpflep(local_index)=nwpflep(i)
                num_wmflep(local_index)=nwmflep(i)
                num_z0flep(local_index)=nz0flep(i)
c
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=4 ! Z
                  iorder(2,local_index)=3
                ELSE
                  iorder(1,local_index)=3 ! Z
                  iorder(2,local_index)=4
                ENDIF
                iorder(3,local_index)=5 ! Z
                iorder(4,local_index)=6
                iorder(7,local_index)=1 ! g
                iorder(5,local_index)=7 ! Z
                iorder(6,local_index)=8
                iorder(8,local_index)=2 ! g

                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                  idplocal(iorder(2,local_index))=idz0(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                  idplocal(iorder(2,local_index))=-idz0(1,i)
                ENDIF
                idplocal(iorder(3,local_index))=idz0(2,i) ! Z
                idplocal(iorder(4,local_index))=-idz0(2,i)
                idplocal(iorder(7,local_index))=21 ! g
                idplocal(iorder(5,local_index))=idz0(3,i) ! Z
                idplocal(iorder(6,local_index))=-idz0(3,i)
                idplocal(iorder(8,local_index))=21 ! g

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(6)=1
                  DO j=1,8
                    igo(j,6)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(6).EQ.0)THEN
                    igr(6)=1
                    DO j=1,4
                      igo(2*j-1,6)=2*index(j,ipp)-1
                      igo(2*j,6)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF(nz0flep(i).EQ.0 .OR. ...)


              IF(nz0flep(i).EQ.0 .OR. nz0flep(i).EQ.1)THEN
*               WRITE(*,*) i,' Z  -> {ggZ} Z  : 67'
                local_index=67
                ialfa(local_index)=1
                n_kphs=n_kphs+1

c Careful in assigning iorders of Z. Should know which Z is leptonic:
c if any, it is represented by idz0(3,i) [last in the list].
c idz0(1,i) represents the t-channel Z boson.
                IF(ipartinitial(1).EQ.-1)THEN
                  iorder(1,local_index)=4 ! Z in
                  iorder(2,local_index)=3
                ELSE
                  iorder(1,local_index)=3 ! Z in
                  iorder(2,local_index)=4
                ENDIF
                iorder(3,local_index)=5 ! Z1
                iorder(4,local_index)=6
                iorder(5,local_index)=1 ! g
                iorder(6,local_index)=2 ! g
                iorder(7,local_index)=7 ! Z2 (possibly leptonic)
                iorder(8,local_index)=8

                IF(ipartinitial(1).EQ.-1)THEN
                  idplocal(iorder(1,local_index))=-idz0(1,i) ! Z in
                  idplocal(iorder(2,local_index))=idz0(1,i)
                ELSE
                  idplocal(iorder(1,local_index))=idz0(1,i) ! Z in
                  idplocal(iorder(2,local_index))=-idz0(1,i)
                ENDIF
                idplocal(iorder(3,local_index))=idz0(2,i) ! Z1
                idplocal(iorder(4,local_index))=-idz0(2,i)
                idplocal(iorder(5,local_index))=21 ! g
                idplocal(iorder(6,local_index))=21 ! g
                idplocal(iorder(7,local_index))=idz0(3,i) ! Z2(possibly 
                                                          !   leptonic)
                idplocal(iorder(8,local_index))=-idz0(3,i)

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(6)=1
                  DO j=1,8
                    igo(j,6)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(6).EQ.0)THEN
                    igr(6)=1
                    DO j=1,4
                      igo(2*j-1,6)=2*index(j,ipp)-1
                      igo(2*j,6)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

              ENDIF             ! IF(nz0flep(i).EQ.0 .OR. ...)


              IF(nz0flep(i).EQ.0)THEN
                IF(idz0(2,i).NE.idz0(3,i))THEN
c This channel differs from the previous one by the simple exchange 
c of the two final state Z's in the phase space 2_4to31.

*                 WRITE(*,*) i,' Z  -> Z {ggZ}  : 68'
                  local_index=68
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  iorder(1,local_index)=iorder(1,local_index-1) ! Z in
                  iorder(2,local_index)=iorder(2,local_index-1)
                  iorder(3,local_index)=iorder(7,local_index-1) ! Z2
                  iorder(4,local_index)=iorder(8,local_index-1)
                  iorder(5,local_index)=iorder(5,local_index-1) ! g
                  iorder(6,local_index)=iorder(6,local_index-1) ! g
                  iorder(7,local_index)=iorder(3,local_index-1) ! Z1
                  iorder(8,local_index)=iorder(4,local_index-1)
                ENDIF           ! IF(idz0(2,i).NE.idz0(3,i))THEN

              ENDIF             ! IF(nz0flep(i).EQ.0)THEN


            ELSE
              print*,'***proc.f : ERROR in qq -> gg X !!!'
              print*,'   Uncorrect matching in nfinal(i).eq.2'
              STOP

            ENDIF               !IF(nwpf(i).EQ.1 .AND. nwmf(i).EQ.1)THEN


          ELSE                  ! nfinal(i).EQ.1
***********************************************************************
* t(0W) {ggW-}   and   t(0W) {ggW+}               ialfa(27) 1_1_31 
* t(W0) {ggW-}   and   t(W0) {ggW+}               ialfa(28) 1_1_31
* t(WW) {ggZ}                                     ialfa(29) 1_1_31
* t(00) {ggZ}                                     ialfa(32) 1_1_31
* t(0 ) W- {ggq} -> t(0 ) W- g{gq}
*          and  t(0 ) W+ {ggq} -> t(0 ) W+ g{gq}  ialfa(43) 1_2_3
* t(0 ) W- {ggq} -> t(0 ) W- q{gg}
*          and  t(0 ) W+ {ggq} -> t(0 ) W+ q{gg}  ialfa(44) 1_2_3
* t(W ) Z  {ggq} -> t(W ) Z  g{gq}                ialfa(45) 1_2_3
* t(W ) Z  {ggq} -> t(W ) Z  q{gg}                ialfa(46) 1_2_3
* qq -> {gq} {gqW-}   and  qq -> {gq} {gqW+}      ialfa(23) 2_4
* t(W+ ) Z {ggq} -> t(W+ ) Z g{gq}                ialfa(47) 1_2_3
* t(W+ ) Z {ggq} -> t(W+ ) Z q{gg}                ialfa(48) 1_2_3
* t(W- ) Z {ggq} -> t(W- ) Z g{gq}                ialfa(49) 1_2_3
* t(W- ) Z {ggq} -> t(W- ) Z q{gg}                ialfa(50) 1_2_3
* t(0[1] ) Z {ggq} -> t(0[1] ) Z g{gq}            ialfa(51) 1_2_3
* t(0[1] ) Z {ggq} -> t(0[1] ) Z q{gg}            ialfa(52) 1_2_3
* t(0[2] ) Z {ggq} -> t(0[2] ) Z g{gq}            ialfa(53) 1_2_3
* t(0[2] ) Z {ggq} -> t(0[2] ) Z q{gg}            ialfa(54) 1_2_3
* qq -> {gq} {gqZ}                                ialfa(24) 2_4
***********************************************************************

c Check for type of t-channel boson.
c Determine if W+ is u-->W+d or  d~-->W+u~ etc
c iswap(i)=1,2 records the order in which the two bosons which include 
c the initial particles appear in the list. Necessary when the two 
c bosons are identical. We first check the order in which they appear in
c the basic order. Then whether they have been interchanged by the 
c ordering by flavour
            IF(ipartinitial(1).EQ.-1)THEN ! 1=antiparticle #1
              itypeinitial(1)=itype(1,i)
              IF(ipartinitial(2).EQ.-1)THEN ! 2=antiparticle #2
                itypeinitial(2)=itype(2,i)
              ELSE              ! 2=particle #1  
                itypeinitial(2)=itype(invindex(1,iperm(i)),i)
              ENDIF
              iflagswap=0
            ELSE                ! 1=particle #1
              itypeinitial(1)=itype(invindex(1,iperm(i)),i)
              IF(ipartinitial(2).EQ.-1)THEN ! 2=antiparticle #1
                itypeinitial(2)=itype(1,i)
                iflagswap=1
              ELSE              ! 2=particle #2  
                itypeinitial(2)=itype(invindex(2,iperm(i)),i)
                IF(invindex(1,iperm(i)).LT.invindex(2,iperm(i)))THEN
                  iflagswap=0
                ELSE
                  iflagswap=1
                ENDIF
              ENDIF
            ENDIF               ! end  Check for type of t-channel boson
            
            IF(itypeinitial(1).EQ.itypeinitial(2))THEN
              IF(itypeinitial(1).EQ.1)THEN
                IF(iposwp(2,i).LT.iposwp(1,i)) iflagswap=1-iflagswap
              ELSEIF(itypeinitial(1).EQ.-1)THEN
                IF(iposwm(2,i).LT.iposwm(1,i)) iflagswap=1-iflagswap
              ELSE
                IF(iposz0(2,i).LT.iposz0(1,i)) iflagswap=1-iflagswap
              ENDIF
            ENDIF
            
            IF(iflagswap.EQ.0)THEN
              iswap(1)=1
              iswap(2)=2
            ELSE
              iswap(1)=2
              iswap(2)=1
            ENDIF


            IF((nwmf(i)+nwpf(i)).EQ.1)THEN 
***********************************************************************
* t(0W) {ggW-}   and   t(0W) {ggW+}               ialfa(27) 1_1_31 
* t(W0) {ggW-}   and   t(W0) {ggW+}               ialfa(28) 1_1_31
* t(0 ) W- {ggq} -> t(0 ) W- g{gq}
*          and  t(0 ) W+ {ggq} -> t(0 ) W+ g{gq}  ialfa(43) 1_2_3
* t(0 ) W- {ggq} -> t(0 ) W- q{gg}
*          and  t(0 ) W+ {ggq} -> t(0 ) W+ q{gg}  ialfa(44) 1_2_3
* t(W ) Z  {ggq} -> t(W ) Z  g{gq}                ialfa(45) 1_2_3
* t(W ) Z  {ggq} -> t(W ) Z  q{gg}                ialfa(46) 1_2_3
* qq -> {gq} {gqW-}   and  qq -> {gq} {gqW+}      ialfa(23) 2_4
***********************************************************************
c Merging cases with one final W+ or W-. Requires determination of 
c final-state W's sign.

c giuseppe ILC 25/06/2007
              IF(abs(iproc(1)).LT.11 .AND. abs(iproc(2)).LT.11)THEN

*               WRITE(*,*) i,'  t(0 ) W {ggq} -> t(0 ) W g{gq} : 43'
                local_index=43
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(4,local_index)=1 ! g
                iorder(5,local_index)=2 ! g
                IF(itypeinitial(1).EQ.0)THEN
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=4 ! Z
                    iorder(8,local_index)=3
                  ELSE
                    iorder(1,local_index)=3
                    iorder(8,local_index)=4
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=6 ! W+
                      iorder(3,local_index)=5
                    ELSE
                      iorder(2,local_index)=5
                      iorder(3,local_index)=6
                    ENDIF
                    iorder(6,local_index)=7 ! W-
                    iorder(7,local_index)=8
                  ELSE          ! nwpf(i).EQ.1
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=8 ! W-
                      iorder(3,local_index)=7
                    ELSE
                      iorder(2,local_index)=7
                      iorder(3,local_index)=8
                    ENDIF
                    iorder(6,local_index)=5 ! W+
                    iorder(7,local_index)=6
                  ENDIF
                ELSE            ! itypeinitial(2).EQ.0
                  IF(ipartinitial(2).EQ.-1)THEN
                    iorder(2,local_index)=4 ! Z
                    iorder(8,local_index)=3
                  ELSE
                    iorder(2,local_index)=3
                    iorder(8,local_index)=4
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(1).EQ.-1)THEN
                      iorder(1,local_index)=6 ! W+
                      iorder(3,local_index)=5
                    ELSE
                      iorder(1,local_index)=5
                      iorder(3,local_index)=6
                    ENDIF
                    iorder(6,local_index)=7 ! W-
                    iorder(7,local_index)=8
                  ELSE          ! nwpf(i).EQ.1
                    IF(ipartinitial(1).EQ.-1)THEN
                      iorder(1,local_index)=8 ! W-
                      iorder(3,local_index)=7
                    ELSE
                      iorder(1,local_index)=7
                      iorder(3,local_index)=8
                    ENDIF
                    iorder(6,local_index)=5 ! W+
                    iorder(7,local_index)=6
                  ENDIF
                ENDIF           ! IF(itypeinitial(1).EQ.0)THEN

                idplocal(iorder(4,local_index))=21 ! g
                idplocal(iorder(5,local_index))=21 ! g
                IF(itypeinitial(1).EQ.0)THEN
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idz0(1,i)
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwp(1,i)-1) !W+
                      idplocal(iorder(3,local_index))=idwp(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwp(1,i)
                      idplocal(iorder(3,local_index))=-(idwp(1,i)-1)
                    ENDIF
                    idplocal(iorder(6,local_index))=idwm(1,i) ! W-
                    idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                  ELSE          ! nwpf(i).EQ.1
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwm(1,i)+1) !W-
                      idplocal(iorder(3,local_index))=idwm(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwm(1,i)
                      idplocal(iorder(3,local_index))=-(idwm(1,i)+1)
                    ENDIF
                    idplocal(iorder(6,local_index))=idwp(1,i) ! W+
                    idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ELSE            ! itypeinitial(2).EQ.0
                  IF(ipartinitial(2).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-idz0(1,i) ! Z
                    idplocal(iorder(8,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idz0(1,i)
                    idplocal(iorder(8,local_index))=-idz0(1,i)
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(1).EQ.-1)THEN
                      idplocal(iorder(1,local_index))=-(idwp(1,i)-1) !W+
                      idplocal(iorder(3,local_index))=idwp(1,i)
                    ELSE
                      idplocal(iorder(1,local_index))=idwp(1,i)
                      idplocal(iorder(3,local_index))=-(idwp(1,i)-1)
                    ENDIF
                    idplocal(iorder(6,local_index))=idwm(1,i) ! W-
                    idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                  ELSE          ! nwpf(i).EQ.1
                    IF(ipartinitial(1).EQ.-1)THEN
                      idplocal(iorder(1,local_index))=-(idwm(1,i)+1) !W-
                      idplocal(iorder(3,local_index))=idwm(1,i)
                    ELSE
                      idplocal(iorder(1,local_index))=idwm(1,i)
                      idplocal(iorder(3,local_index))=-(idwm(1,i)+1)
                    ENDIF
                    idplocal(iorder(6,local_index))=idwp(1,i) ! W+
                    idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ENDIF           ! IF(itypeinitial(1).EQ.0)THEN

                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF

c Channel ialfa(44) differs from ialfa(43)  by the simple exchange 
c g <-> q in the phase space 1_2_3.
*               WRITE(*,*) i,'  t(0 ) W {ggq} -> t(0 ) W q{gg} : 44'
                local_index=44
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=iorder(1,local_index-1)
                iorder(2,local_index)=iorder(2,local_index-1)
                iorder(3,local_index)=iorder(5,local_index-1) ! g
                iorder(4,local_index)=iorder(4,local_index-1) ! g
                iorder(5,local_index)=iorder(3,local_index-1)
                iorder(6,local_index)=iorder(6,local_index-1) ! W
                iorder(7,local_index)=iorder(7,local_index-1)
                iorder(8,local_index)=iorder(8,local_index-1)

c Channel ialfa(45) differs from ialfa(43)  by the simple exchange
c 1<->2 and 3<->8 in the phase space 1_2_3.
*               WRITE(*,*) i,'  t(W ) W  {ggq} -> t(W ) W  g{gq} : 45'
                local_index=45
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=iorder(1,local_index-2)
                iorder(2,local_index)=iorder(2,local_index-2)
                iorder(3,local_index)=iorder(8,local_index-2) !g
                iorder(4,local_index)=iorder(4,local_index-2) !g
                iorder(5,local_index)=iorder(5,local_index-2)
                iorder(6,local_index)=iorder(6,local_index-2) !W
                iorder(7,local_index)=iorder(7,local_index-2)
                iorder(8,local_index)=iorder(3,local_index-2)

c Channel ialfa(46) differs from ialfa(45)  by the simple exchange
c 3<->5 in the phase space 1_2_3.
*               WRITE(*,*) i,'  t(W ) W  {ggq} -> t(W ) W  q{gg} : 46'
                local_index=46
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(1,local_index)=iorder(1,local_index-1)
                iorder(2,local_index)=iorder(2,local_index-1)
                iorder(3,local_index)=iorder(5,local_index-1) ! g
                iorder(4,local_index)=iorder(4,local_index-1) ! g
                iorder(5,local_index)=iorder(3,local_index-1)
                iorder(6,local_index)=iorder(6,local_index-1) ! W
                iorder(7,local_index)=iorder(7,local_index-1)
                iorder(8,local_index)=iorder(8,local_index-1)


*               WRITE(*,*) i,'  qq -> {gq} {gqW}  : 23'
                local_index=23
                ialfa(local_index)=1
                n_kphs=n_kphs+1

                iorder(3,local_index)=1 ! g
                iorder(7,local_index)=2 ! g
                IF(itypeinitial(1).EQ.0)THEN
                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=4 ! Z
                    iorder(4,local_index)=3
                  ELSE
                    iorder(1,local_index)=3 ! Z
                    iorder(4,local_index)=4
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=6 ! W+
                      iorder(8,local_index)=5
                    ELSE
                      iorder(2,local_index)=5 ! W+
                      iorder(8,local_index)=6
                    ENDIF
                    iorder(5,local_index)=7 ! W-
                    iorder(6,local_index)=8
                  ELSE
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=8 ! W-
                      iorder(8,local_index)=7
                    ELSE
                      iorder(2,local_index)=7 ! W-
                      iorder(8,local_index)=8
                    ENDIF
                    iorder(5,local_index)=5 ! W+
                    iorder(6,local_index)=6
                  ENDIF
                ELSE            ! itypeinitial(2).EQ.0
                  IF(ipartinitial(2).EQ.-1)THEN
                    iorder(2,local_index)=4 ! Z
                    iorder(4,local_index)=3
                  ELSE
                    iorder(2,local_index)=3 ! Z
                    iorder(4,local_index)=4
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(1).EQ.-1)THEN
                      iorder(1,local_index)=6 ! W+
                      iorder(8,local_index)=5
                    ELSE
                      iorder(1,local_index)=5 ! W+
                      iorder(8,local_index)=6
                    ENDIF
                    iorder(5,local_index)=7 ! W-
                    iorder(6,local_index)=8
                  ELSE
                    IF(ipartinitial(1).EQ.-1)THEN
                      iorder(1,local_index)=8 ! W-
                      iorder(8,local_index)=7
                    ELSE
                      iorder(1,local_index)=7 ! W-
                      iorder(8,local_index)=8
                    ENDIF
                    iorder(5,local_index)=5 ! W+
                    iorder(6,local_index)=6
                  ENDIF
                ENDIF

                idplocal(iorder(3,local_index))=21 ! g
                idplocal(iorder(7,local_index))=21 ! g
                IF(itypeinitial(1).EQ.0)THEN
                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                    idplocal(iorder(4,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                    idplocal(iorder(4,local_index))=-idz0(1,i)
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwp(1,i)-1) !W+
                      idplocal(iorder(8,local_index))=idwp(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                      idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                    ENDIF
                    idplocal(iorder(5,local_index))=idwm(1,i) ! W-
                    idplocal(iorder(6,local_index))=-(idwm(1,i)+1)
                  ELSE
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwm(1,i)+1) !W-
                      idplocal(iorder(8,local_index))=idwm(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                      idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                    ENDIF
                    idplocal(iorder(5,local_index))=idwp(1,i) ! W+
                    idplocal(iorder(6,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ELSE            ! itypeinitial(2).EQ.0
                  IF(ipartinitial(2).EQ.-1)THEN
                    idplocal(iorder(2,local_index))=-idz0(1,i) ! Z
                    idplocal(iorder(4,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(2,local_index))=idz0(1,i) ! Z
                    idplocal(iorder(4,local_index))=-idz0(1,i)
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(1).EQ.-1)THEN
                      idplocal(iorder(1,local_index))=-(idwp(1,i)-1) !W+
                      idplocal(iorder(8,local_index))=idwp(1,i)
                    ELSE
                      idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                      idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                    ENDIF
                    idplocal(iorder(5,local_index))=idwm(1,i) ! W-
                    idplocal(iorder(6,local_index))=-(idwm(1,i)+1)
                  ELSE
                    IF(ipartinitial(1).EQ.-1)THEN
                      idplocal(iorder(1,local_index))=-(idwm(1,i)+1) !W-
                      idplocal(iorder(8,local_index))=idwm(1,i)
                    ELSE
                      idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                      idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                    ENDIF
                    idplocal(iorder(5,local_index))=idwp(1,i) ! W+
                    idplocal(iorder(6,local_index))=-(idwp(1,i)-1)
                  ENDIF
                ENDIF
              
                IF(n_kphs.eq.1)THEN
                  iabase=local_index
                  igr(5)=1
                  DO j=1,8
                    igo(j,5)=j
                  ENDDO
                  DO j=1,8
                    idp(j)=idplocal(j)
                  ENDDO
                  DO j=1,2
                    idp_inout(iorder(j,iabase))=-1
                  ENDDO
                ELSE
                  ial=local_index
                  CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                 idplocal,idp,idp_inout,ipp,ipa)
                  IF(igr(5).EQ.0)THEN
                    igr(5)=1
                    DO j=1,4
                      igo(2*j-1,5)=2*index(j,ipp)-1
                      igo(2*j,5)=2*index(j,ipa)
                    ENDDO
                  ENDIF
                ENDIF
              ENDIF             ! IF( abs(iproc(1)).LT.11 .AND. ...
c end giuseppe ILC 25/06/2007


              IF(nwpflep(i).EQ.0 .AND. nwmflep(i).EQ.0)THEN
                IF(itypeinitial(1).EQ.0 .OR. idout(1).EQ.idout(2))THEN
*                WRITE(*,*) i,'  t(0W) {ggW} : 27'
                  local_index=27
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  IF(ipartinitial(1).EQ.-1)THEN
                    iorder(1,local_index)=4 ! Z
                    iorder(7,local_index)=3
                  ELSE
                    iorder(1,local_index)=3 ! Z
                    iorder(7,local_index)=4
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=6 ! W+
                      iorder(8,local_index)=5
                    ELSE
                      iorder(2,local_index)=5 ! W+
                      iorder(8,local_index)=6
                    ENDIF
                    iorder(3,local_index)=7 ! W-
                    iorder(4,local_index)=8
                  ELSE          ! nwpf(i).EQ.1
                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=8 ! W-
                      iorder(8,local_index)=7
                    ELSE
                      iorder(2,local_index)=7 ! W-
                      iorder(8,local_index)=8
                    ENDIF
                    iorder(3,local_index)=5 ! W+
                    iorder(4,local_index)=6
                  ENDIF
                  iorder(5,local_index)=1 ! g
                  iorder(6,local_index)=2 ! g

                  IF(ipartinitial(1).EQ.-1)THEN
                    idplocal(iorder(1,local_index))=-idz0(1,i) ! Z
                    idplocal(iorder(7,local_index))=idz0(1,i)
                  ELSE
                    idplocal(iorder(1,local_index))=idz0(1,i) ! Z
                    idplocal(iorder(7,local_index))=-idz0(1,i)
                  ENDIF
                  IF(nwmf(i).EQ.1)THEN
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwp(1,i)-1) !W+
                      idplocal(iorder(8,local_index))=idwp(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                      idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                    ENDIF
                    idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                    idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                  ELSE          ! nwpf(i).EQ.1
                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-(idwm(1,i)+1) !W-
                      idplocal(iorder(8,local_index))=idwm(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                      idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                    ENDIF
                    idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                    idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                  ENDIF
                  idplocal(iorder(5,local_index))=21 ! g
                  idplocal(iorder(6,local_index))=21 ! g

                  IF(n_kphs.eq.1)THEN
                    iabase=local_index
                    igr(5)=1
                    DO j=1,8
                      igo(j,5)=j
                    ENDDO
                    DO j=1,8
                      idp(j)=idplocal(j)
                    ENDDO
                    DO j=1,2
                      idp_inout(iorder(j,iabase))=-1
                    ENDDO
                  ELSE
                    ial=local_index
                    CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &                   idplocal,idp,idp_inout,ipp,ipa)
                    IF(igr(5).EQ.0)THEN
                      igr(5)=1
                      DO j=1,4
                        igo(2*j-1,5)=2*index(j,ipp)-1
                        igo(2*j,5)=2*index(j,ipa)
                      ENDDO
                    ENDIF
                  ENDIF

                ENDIF           ! IF(itypeinitial(1).EQ.0 .OR. ...)
                
                IF(itypeinitial(2).EQ.0 .OR. idout(1).EQ.idout(2))THEN
*                WRITE(*,*) i,' t(W0) {ggW}  : 28'
                  local_index=28
                  ialfa(local_index)=1
                  n_kphs=n_kphs+1

                  IF(idout(1).NE.idout(2))THEN

                    IF(ipartinitial(2).EQ.-1)THEN
                      iorder(2,local_index)=4 ! Z
                      iorder(8,local_index)=3
                    ELSE
                      iorder(2,local_index)=3 ! Z
                      iorder(8,local_index)=4
                    ENDIF
                    IF(nwmf(i).EQ.1)THEN
                      IF(ipartinitial(1).EQ.-1)THEN
                        iorder(1,local_index)=6 ! W+
                        iorder(7,local_index)=5
                      ELSE
                        iorder(1,local_index)=5 ! W+
                        iorder(7,local_index)=6
                      ENDIF
                      iorder(3,local_index)=7 ! W-
                      iorder(4,local_index)=8
                    ELSE        ! nwpf(i).EQ.1
                      IF(ipartinitial(1).EQ.-1)THEN
                        iorder(1,local_index)=8 ! W-
                        iorder(7,local_index)=7
                      ELSE
                        iorder(1,local_index)=7 ! W-
                        iorder(7,local_index)=8
                      ENDIF
                      iorder(3,local_index)=5 ! W+
                      iorder(4,local_index)=6
                    ENDIF
                    iorder(5,local_index)=1 ! g
                    iorder(6,local_index)=2 ! g

                    IF(ipartinitial(2).EQ.-1)THEN
                      idplocal(iorder(2,local_index))=-idz0(1,i) ! Z
                      idplocal(iorder(8,local_index))=idz0(1,i)
                    ELSE
                      idplocal(iorder(2,local_index))=idz0(1,i) ! Z
                      idplocal(iorder(8,local_index))=-idz0(1,i)
                    ENDIF
                    IF(nwmf(i).EQ.1)THEN
                      IF(ipartinitial(1).EQ.-1)THEN
                       idplocal(iorder(1,local_index))=-(idwp(1,i)-1)!W+
                       idplocal(iorder(7,local_index))=idwp(1,i)
                     ELSE
                       idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                       idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                     ENDIF
                     idplocal(iorder(3,local_index))=idwm(1,i) ! W-
                     idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                   ELSE         ! nwpf(i).EQ.1
                     IF(ipartinitial(1).EQ.-1)THEN
                       idplocal(iorder(1,local_index))=-(idwm(1,i)+1)!W-
                       idplocal(iorder(7,local_index))=idwm(1,i)
                     ELSE
                       idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                       idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                     ENDIF
                     idplocal(iorder(3,local_index))=idwp(1,i) ! W+
                     idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                   ENDIF
                   idplocal(iorder(5,local_index))=21 ! g
                   idplocal(iorder(6,local_index))=21 ! g

                   IF(n_kphs.eq.1)THEN
                     iabase=local_index
                     igr(5)=1
                     DO j=1,8
                       igo(j,5)=j
                     ENDDO
                     DO j=1,8
                       idp(j)=idplocal(j)
                     ENDDO
                     DO j=1,2
                       idp_inout(iorder(j,iabase))=-1
                     ENDDO
                   ELSE
                     ial=local_index
                     CALL REORDER(iorder(1,ial),iorder(1,iabase),
     &                    npart+1,idplocal,idp,idp_inout,ipp,ipa)
                     IF(igr(5).EQ.0)THEN
                       igr(5)=1
                       DO j=1,4
                         igo(2*j-1,5)=2*index(j,ipp)-1
                         igo(2*j,5)=2*index(j,ipa)
                       ENDDO
                     ENDIF
                   ENDIF

                 ELSE           ! iabase should be 27
                   local_index=28
                   iorder(1,local_index)=iorder(2,local_index-1)
                   iorder(2,local_index)=iorder(1,local_index-1)
                   iorder(3,local_index)=iorder(3,local_index-1)
                   iorder(4,local_index)=iorder(4,local_index-1)
                   iorder(5,local_index)=iorder(5,local_index-1)
                   iorder(6,local_index)=iorder(6,local_index-1)
                   iorder(7,local_index)=iorder(8,local_index-1)
                   iorder(8,local_index)=iorder(7,local_index-1)
 
                 ENDIF          ! IF(idout(1).NE.idout(2))THEN

               ENDIF            ! IF(itypeinitial(2).EQ.0 .OR. ...)
             ENDIF              ! IF(nwpflep(i).EQ.0 .AND. ...)


           ELSEIF(nz0f(i).EQ.1)THEN
***********************************************************************
* t(WW) {ggZ}                                     ialfa(29) 1_1_31
* t(00) {ggZ}                                     ialfa(32) 1_1_31
* t(W+ ) Z {ggq} -> t(W+ ) Z g{gq}                ialfa(47) 1_2_3
* t(W+ ) Z {ggq} -> t(W+ ) Z q{gg}                ialfa(48) 1_2_3
* t(W- ) Z {ggq} -> t(W- ) Z g{gq}                ialfa(49) 1_2_3
* t(W- ) Z {ggq} -> t(W- ) Z q{gg}                ialfa(50) 1_2_3
* t(0[1] ) Z {ggq} -> t(0[1] ) Z g{gq}            ialfa(51) 1_2_3
* t(0[1] ) Z {ggq} -> t(0[1] ) Z q{gg}            ialfa(52) 1_2_3
* t(0[2] ) Z {ggq} -> t(0[2] ) Z g{gq}            ialfa(53) 1_2_3
* t(0[2] ) Z {ggq} -> t(0[2] ) Z q{gg}            ialfa(54) 1_2_3
* qq -> {gq} {gqZ}                                ialfa(24) 2_4
***********************************************************************
             IF(nwm(i).EQ.1 .AND. nwp(i).EQ.1)THEN

c giuseppe ILC 25/06/2007
               IF(abs(iproc(1)).LT.11 .AND. abs(iproc(2)).LT.11)THEN

*                WRITE(*,*) i,'t(W+ ) Z {ggq} -> t(W+ ) Z g{gq}  : 47'
                 local_index=47
                 ialfa(local_index)=1
                 n_kphs=n_kphs+1

                 iorder(4,local_index)=1 ! g
                 iorder(5,local_index)=2 ! g
                 iorder(6,local_index)=3 ! Z
                 iorder(7,local_index)=4
                 IF(itypeinitial(1).EQ.1)THEN
                   IF(ipartinitial(1).EQ.-1)THEN
                     iorder(1,local_index)=6 ! W+
                     iorder(8,local_index)=5
                   ELSE
                     iorder(1,local_index)=5 ! W+
                     iorder(8,local_index)=6
                   ENDIF
                   IF(ipartinitial(2).EQ.-1)THEN
                     iorder(2,local_index)=8 ! W-
                     iorder(3,local_index)=7
                   ELSE
                     iorder(2,local_index)=7 ! W-
                     iorder(3,local_index)=8
                   ENDIF
                 ELSE           ! itypeinitial(2).EQ.1
                   IF(ipartinitial(2).EQ.-1)THEN
                     iorder(2,local_index)=6 ! W+
                     iorder(8,local_index)=5
                   ELSE
                     iorder(2,local_index)=5 ! W+
                     iorder(8,local_index)=6
                   ENDIF
                   IF(ipartinitial(1).EQ.-1)THEN
                     iorder(1,local_index)=8 ! W-
                     iorder(3,local_index)=7
                   ELSE
                     iorder(1,local_index)=7 ! W-
                     iorder(3,local_index)=8
                   ENDIF
                 ENDIF

                 idplocal(iorder(4,local_index))=21 ! g
                 idplocal(iorder(5,local_index))=21 ! g
                 idplocal(iorder(6,local_index))=idz0(1,i) ! Z
                 idplocal(iorder(7,local_index))=-idz0(1,i)
                 IF(itypeinitial(1).EQ.1)THEN
                   IF(ipartinitial(1).EQ.-1)THEN
                     idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
                     idplocal(iorder(8,local_index))=idwp(1,i)
                   ELSE
                     idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                     idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                   ENDIF
                   IF(ipartinitial(2).EQ.-1)THEN
                     idplocal(iorder(2,local_index))=-(idwm(1,i)+1) ! W-
                     idplocal(iorder(3,local_index))=idwm(1,i)
                   ELSE
                     idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                     idplocal(iorder(3,local_index))=-(idwm(1,i)+1)
                   ENDIF
                 ELSE           ! itypeinitial(2).EQ.1
                   IF(ipartinitial(2).EQ.-1)THEN
                     idplocal(iorder(2,local_index))=-(idwp(1,i)-1) ! W+
                     idplocal(iorder(8,local_index))=idwp(1,i)
                   ELSE
                     idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                     idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                   ENDIF
                   IF(ipartinitial(1).EQ.-1)THEN
                     idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
                     idplocal(iorder(3,local_index))=idwm(1,i)
                   ELSE
                     idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                     idplocal(iorder(3,local_index))=-(idwm(1,i)+1)
                   ENDIF
                 ENDIF

                 IF(n_kphs.eq.1)THEN
                   iabase=local_index
                   igr(5)=1
                   DO j=1,8
                     igo(j,5)=j
                   ENDDO
                   DO j=1,8
                     idp(j)=idplocal(j)
                   ENDDO
                   DO j=1,2
                     idp_inout(iorder(j,iabase))=-1
                   ENDDO
                 ELSE
                   ial=local_index
                   CALL REORDER(iorder(1,ial),iorder(1,iabase),
     &                  npart+1,idplocal,idp,idp_inout,ipp,ipa)
                   IF(igr(5).EQ.0)THEN
                     igr(5)=1
                     DO j=1,4
                       igo(2*j-1,5)=2*index(j,ipp)-1
                       igo(2*j,5)=2*index(j,ipa)
                     ENDDO
                   ENDIF
                 ENDIF


*     WRITE(*,*) i,'t(W+ ) Z {ggq} -> t(W+ ) Z q{gg}  : 48'
c Channel ialfa(48) differs from ialfa(47)  by the simple exchange 
c 3<->5 (g<->q) in the phase space 1_2_3.
                 local_index=48
                 ialfa(local_index)=1
                 n_kphs=n_kphs+1

                 iorder(1,local_index)=iorder(1,local_index-1)
                 iorder(2,local_index)=iorder(2,local_index-1)
                 iorder(3,local_index)=iorder(5,local_index-1) ! g
                 iorder(4,local_index)=iorder(4,local_index-1) ! g
                 iorder(5,local_index)=iorder(3,local_index-1) ! q
                 iorder(6,local_index)=iorder(6,local_index-1) ! Z
                 iorder(7,local_index)=iorder(7,local_index-1)
                 iorder(8,local_index)=iorder(8,local_index-1)

*     WRITE(*,*) i,'t(W- ) Z {ggq} -> t(W- ) Z g{gq}  : 49'
c Channel ialfa(49) differs from ialfa(47)  by the simple exchange 
c 3<->8 in the phase space 1_2_3.
                 local_index=49
                 ialfa(local_index)=1
                 n_kphs=n_kphs+1

                 iorder(1,local_index)=iorder(1,local_index-2)
                 iorder(2,local_index)=iorder(2,local_index-2)
                 iorder(3,local_index)=iorder(8,local_index-2) ! q
                 iorder(4,local_index)=iorder(4,local_index-2) ! g
                 iorder(5,local_index)=iorder(5,local_index-2) ! g
                 iorder(6,local_index)=iorder(6,local_index-2) ! Z
                 iorder(7,local_index)=iorder(7,local_index-2)
                 iorder(8,local_index)=iorder(3,local_index-2)

*     WRITE(*,*) i,'t(W- ) Z {ggq} -> t(W- ) Z q{gg}  : 50'
c Channel ialfa(50) differs from ialfa(49)  by the simple exchange 
c 3<->5 (g<->q) in the phase space 1_2_3.
                 local_index=50
                 ialfa(local_index)=1
                 n_kphs=n_kphs+1

                 iorder(1,local_index)=iorder(1,local_index-1)
                 iorder(2,local_index)=iorder(2,local_index-1)
                 iorder(3,local_index)=iorder(5,local_index-1) ! g
                 iorder(4,local_index)=iorder(4,local_index-1) ! g
                 iorder(5,local_index)=iorder(3,local_index-1) ! q
                 iorder(6,local_index)=iorder(6,local_index-1) ! Z
                 iorder(7,local_index)=iorder(7,local_index-1)
                 iorder(8,local_index)=iorder(8,local_index-1)


c The following channel merges two possibilities: the four quarks 
c denoted by 'q' may be associated to two W or two Z bosons. 
c If the process can be both of type 2gZWW and 2g3Z, calling two 
c different channels of type qq -> {gq} {gqZ} is redundant because the 
c difference in the ordering is irrelevant. As in this case the channel 
c qq -> {gq} {gqZ} would be evaluated twice (it fulfills two different 
c matchings), the flag ifound_phsp is necessary to avoid double counting
c on the variable n_kphs.

                 IF(ifound_phsp.EQ.0)THEN
                   ifound_phsp=1

*     WRITE(*,*) i,'qq -> {gq} {gqZ}  : 24'
                   local_index=24
                   ialfa(local_index)=1
                   n_kphs=n_kphs+1

                   iorder(3,local_index)=1 ! g
                   iorder(7,local_index)=2 ! g
                   iorder(5,local_index)=3 ! Z
                   iorder(6,local_index)=4
                   IF(itypeinitial(1).EQ.1)THEN ! First W is +
                     IF(ipartinitial(1).EQ.-1)THEN
                       iorder(1,local_index)=6 ! W+
                       iorder(4,local_index)=5
                     ELSE
                       iorder(1,local_index)=5 ! W+
                       iorder(4,local_index)=6
                     ENDIF
                     IF(ipartinitial(2).EQ.-1)THEN
                       iorder(2,local_index)=8 ! W-
                       iorder(8,local_index)=7
                     ELSE
                       iorder(2,local_index)=7 ! W-
                       iorder(8,local_index)=8
                     ENDIF
                   ELSE         ! First W is -
                     IF(ipartinitial(1).EQ.-1)THEN
                       iorder(1,local_index)=8 ! W-
                       iorder(4,local_index)=7
                     ELSE
                       iorder(1,local_index)=7 ! W-
                       iorder(4,local_index)=8
                     ENDIF
                     IF(ipartinitial(2).EQ.-1)THEN
                       iorder(2,local_index)=6 ! W+
                       iorder(8,local_index)=5
                     ELSE
                       iorder(2,local_index)=5 ! W+
                       iorder(8,local_index)=6
                     ENDIF
                   ENDIF

                   idplocal(iorder(3,local_index))=21 ! g
                   idplocal(iorder(7,local_index))=21 ! g
                   idplocal(iorder(5,local_index))=idz0(1,i) ! Z
                   idplocal(iorder(6,local_index))=-idz0(1,i)
                   IF(itypeinitial(1).EQ.1)THEN ! First W is +
                     IF(ipartinitial(1).EQ.-1)THEN
                       idplocal(iorder(1,local_index))=-(idwp(1,i)-1)!W+
                       idplocal(iorder(4,local_index))=idwp(1,i)
                     ELSE
                       idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                       idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
                     ENDIF
                     IF(ipartinitial(2).EQ.-1)THEN
                       idplocal(iorder(2,local_index))=-(idwm(1,i)+1)!W-
                       idplocal(iorder(8,local_index))=idwm(1,i)
                     ELSE
                       idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                       idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                     ENDIF
                   ELSE         ! First W is -
                     IF(ipartinitial(1).EQ.-1)THEN
                       idplocal(iorder(1,local_index))=-(idwm(1,i)+1)!W-
                       idplocal(iorder(4,local_index))=idwm(1,i)
                     ELSE
                       idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                       idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
                     ENDIF
                     IF(ipartinitial(2).EQ.-1)THEN
                       idplocal(iorder(2,local_index))=-(idwp(1,i)-1)!W+
                       idplocal(iorder(8,local_index))=idwp(1,i)
                     ELSE
                       idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                       idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                     ENDIF
                   ENDIF

                   IF(n_kphs.eq.1)THEN
                     iabase=local_index
                     igr(5)=1
                     DO j=1,8
                       igo(j,5)=j
                     ENDDO
                     DO j=1,8
                       idp(j)=idplocal(j)
                     ENDDO
                     DO j=1,2
                       idp_inout(iorder(j,iabase))=-1
                     ENDDO
                   ELSE
                     ial=local_index
                     CALL REORDER(iorder(1,ial),iorder(1,iabase),
     &                    npart+1,idplocal,idp,idp_inout,ipp,ipa)
                     IF(igr(5).EQ.0)THEN
                       igr(5)=1
                       DO j=1,4
                         igo(2*j-1,5)=2*index(j,ipp)-1
                         igo(2*j,5)=2*index(j,ipa)
                       ENDDO
                     ENDIF
                   ENDIF
                 ENDIF          ! IF(ifound_phsp.EQ.0)THEN
               ENDIF            ! IF(abs(iproc(1)).LT.11 .AND. ...
c end giuseppe ILC 25/06/2007

               IF(nz0flep(i).EQ.0)THEN
*                WRITE(*,*) i,' t(WW) {ggZ}  : 29'
                 local_index=29
                 ialfa(local_index)=1
                 n_kphs=n_kphs+1

                 iorder(3,local_index)=3 ! Z
                 iorder(4,local_index)=4
                 iorder(5,local_index)=1 ! g
                 iorder(6,local_index)=2 ! g
                 IF(itypeinitial(1).EQ.1)THEN ! First W is +
                   IF(ipartinitial(1).EQ.-1)THEN
                     iorder(1,local_index)=6 ! W+
                     iorder(7,local_index)=5
                   ELSE
                     iorder(1,local_index)=5 ! W+
                     iorder(7,local_index)=6
                   ENDIF
                   IF(ipartinitial(2).EQ.-1)THEN
                     iorder(2,local_index)=8 ! W-
                     iorder(8,local_index)=7
                   ELSE
                     iorder(2,local_index)=7 ! W-
                     iorder(8,local_index)=8
                   ENDIF
                 ELSEIF(itypeinitial(1).EQ.-1)THEN ! First W is -
                   IF(ipartinitial(1).EQ.-1)THEN
                     iorder(1,local_index)=8 ! W-
                     iorder(7,local_index)=7
                   ELSE
                     iorder(1,local_index)=7 ! W-
                     iorder(7,local_index)=8
                   ENDIF
                   IF(ipartinitial(2).EQ.-1)THEN
                     iorder(2,local_index)=6 ! W+
                     iorder(8,local_index)=5
                   ELSE
                     iorder(2,local_index)=5 ! W+
                     iorder(8,local_index)=6
                   ENDIF
                 ELSE
                   print*,'***proc.f:ERROR in itypeinitial'
     &                  //' determination'
                   print*,'itypeinitial(1)',itypeinitial(1)
                   print*,'itypeinitial(2)',itypeinitial(2)
                   STOP
                 ENDIF

                 idplocal(iorder(3,local_index))=idz0(1,i) ! Z
                 idplocal(iorder(4,local_index))=-idz0(1,i)
                 idplocal(iorder(5,local_index))=21 ! g
                 idplocal(iorder(6,local_index))=21 ! g
                 IF(itypeinitial(1).EQ.1)THEN ! First W is +
                   IF(ipartinitial(1).EQ.-1)THEN
                     idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
                     idplocal(iorder(7,local_index))=idwp(1,i)
                   ELSE
                     idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                     idplocal(iorder(7,local_index))=-(idwp(1,i)-1)
                   ENDIF
                   IF(ipartinitial(2).EQ.-1)THEN
                     idplocal(iorder(2,local_index))=-(idwm(1,i)+1) ! W-
                     idplocal(iorder(8,local_index))=idwm(1,i)
                   ELSE
                     idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                     idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
                   ENDIF
                 ELSE           ! First W is -
                   IF(ipartinitial(1).EQ.-1)THEN
                     idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
                     idplocal(iorder(7,local_index))=idwm(1,i)
                   ELSE
                     idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                     idplocal(iorder(7,local_index))=-(idwm(1,i)+1)
                   ENDIF
                   IF(ipartinitial(2).EQ.-1)THEN
                     idplocal(iorder(2,local_index))=-(idwp(1,i)-1) ! W+
                     idplocal(iorder(8,local_index))=idwp(1,i)
                   ELSE
                     idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                     idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
                   ENDIF
                 ENDIF
                   
                 IF(n_kphs.eq.1)THEN
                   iabase=local_index
                   igr(5)=1
                   DO j=1,8
                     igo(j,5)=j
                   ENDDO
                   DO j=1,8
                     idp(j)=idplocal(j)
                   ENDDO
                   DO j=1,2
                     idp_inout(iorder(j,iabase))=-1
                   ENDDO
                 ELSE
                   ial=local_index
                   CALL REORDER(iorder(1,ial),iorder(1,iabase),
     &                  npart+1,idplocal,idp,idp_inout,ipp,ipa)
                   IF(igr(5).EQ.0)THEN
                     igr(5)=1
                     DO j=1,4
                       igo(2*j-1,5)=2*index(j,ipp)-1
                       igo(2*j,5)=2*index(j,ipa)
                     ENDDO
                   ENDIF
                 ENDIF
               ENDIF            ! IF(nz0flep(i).EQ.0)THEN


             ELSE               ! two neutral bosons in t-channel

c giuseppe ILC 25/06/2007
               IF(abs(iproc(1)).LT.11 .AND. abs(iproc(2)).LT.11)THEN

*     WRITE(*,*) i,' t(0[1] ) Z {ggq} -> t(0[1] ) Z g{gq}  : 51'
                 local_index=51
                 ialfa(local_index)=1
                 n_kphs=n_kphs+1

                 IF(ipartinitial(iswap(1)).EQ.-1)THEN
                   iorder(iswap(1),local_index)=4 ! Z
                   iorder(8,local_index)=3
                 ELSE
                   iorder(iswap(1),local_index)=3
                   iorder(8,local_index)=4
                 ENDIF
                 IF(ipartinitial(iswap(2)).EQ.-1)THEN
                   iorder(iswap(2),local_index)=6 ! Z
                   iorder(3,local_index)=5
                 ELSE
                   iorder(iswap(2),local_index)=5
                   iorder(3,local_index)=6
                 ENDIF
                 iorder(4,local_index)=1 ! g
                 iorder(5,local_index)=2 ! g
                 iorder(6,local_index)=7 ! Z
                 iorder(7,local_index)=8

                 IF(ipartinitial(iswap(1)).EQ.-1)THEN
                   idplocal(iorder(iswap(1),local_index))=-idz0(1,i) ! Z
                   idplocal(iorder(8,local_index))=idz0(1,i)
                 ELSE
                   idplocal(iorder(iswap(1),local_index))=idz0(1,i)
                   idplocal(iorder(8,local_index))=-idz0(1,i)
                 ENDIF
                 IF(ipartinitial(iswap(2)).EQ.-1)THEN
                   idplocal(iorder(iswap(2),local_index))=-idz0(2,i) ! Z
                   idplocal(iorder(3,local_index))=idz0(2,i)
                 ELSE
                   idplocal(iorder(iswap(2),local_index))=idz0(2,i)
                   idplocal(iorder(3,local_index))=-idz0(2,i)
                 ENDIF
                 idplocal(iorder(4,local_index))=21 ! g
                 idplocal(iorder(5,local_index))=21 ! g
                 idplocal(iorder(6,local_index))=idz0(3,i) ! Z
                 idplocal(iorder(7,local_index))=-idz0(3,i)

                 IF(n_kphs.eq.1)THEN
                   iabase=local_index
                   igr(6)=1
                   DO j=1,8
                     igo(j,6)=j
                   ENDDO
                   DO j=1,8
                     idp(j)=idplocal(j)
                   ENDDO
                   DO j=1,2
                     idp_inout(iorder(j,iabase))=-1
                   ENDDO
                 ELSE
                   ial=local_index
                   CALL REORDER(iorder(1,ial),iorder(1,iabase),
     &                  npart+1,idplocal,idp,idp_inout,ipp,ipa)
                   IF(igr(6).EQ.0)THEN
                     igr(6)=1
                     DO j=1,4
                       igo(2*j-1,6)=2*index(j,ipp)-1
                       igo(2*j,6)=2*index(j,ipa)
                     ENDDO
                   ENDIF
                 ENDIF

*     WRITE(*,*) i,'t(0[1] ) Z {ggq} -> t(0[1] ) Z q{gg}  : 52'
c Channel ialfa(52) differs from ialfa(51)  by the simple exchange 
c 3<->5 (g<->q) in the phase space 1_2_3.
                 local_index=52
                 ialfa(local_index)=1
                 n_kphs=n_kphs+1

                 iorder(1,local_index)=iorder(1,local_index-1)
                 iorder(2,local_index)=iorder(2,local_index-1)
                 iorder(3,local_index)=iorder(5,local_index-1) ! g
                 iorder(4,local_index)=iorder(4,local_index-1) ! g
                 iorder(5,local_index)=iorder(3,local_index-1) ! q
                 iorder(6,local_index)=iorder(6,local_index-1) ! Z
                 iorder(7,local_index)=iorder(7,local_index-1)
                 iorder(8,local_index)=iorder(8,local_index-1)


                 IF(idout(1).NE.idout(2))THEN
*     WRITE(*,*) i,'t(0[2] ) Z {ggq} -> t(0[2] ) Z g{gq}  : 53'
c Channel ialfa(53) differs from ialfa(51)  by the simple exchange 
c 3<->8 in the phase space 1_2_3. t(0[2] ) means that the second neutral
c boson is mapped in t-channel.
                   local_index=53
                   ialfa(local_index)=1
                   n_kphs=n_kphs+1

                   iorder(1,local_index)=iorder(1,local_index-2)
                   iorder(2,local_index)=iorder(2,local_index-2)
                   iorder(3,local_index)=iorder(8,local_index-2) ! q
                   iorder(4,local_index)=iorder(4,local_index-2) ! g
                   iorder(5,local_index)=iorder(5,local_index-2) ! g
                   iorder(6,local_index)=iorder(6,local_index-2) ! Z
                   iorder(7,local_index)=iorder(7,local_index-2)
                   iorder(8,local_index)=iorder(3,local_index-2)

*     WRITE(*,*) i,'t(0[2] ) Z {ggq} -> t(0[2] ) Z q{gg}  : 54'
c Channel ialfa(54) differs from ialfa(53)  by the simple exchange 
c 3<->5 (g<->q) in the phase space 1_2_3.
                   local_index=54
                   ialfa(local_index)=1
                   n_kphs=n_kphs+1

                   iorder(1,local_index)=iorder(1,local_index-1)
                   iorder(2,local_index)=iorder(2,local_index-1)
                   iorder(3,local_index)=iorder(5,local_index-1) ! g
                   iorder(4,local_index)=iorder(4,local_index-1) ! g
                   iorder(5,local_index)=iorder(3,local_index-1) ! q
                   iorder(6,local_index)=iorder(6,local_index-1) ! Z
                   iorder(7,local_index)=iorder(7,local_index-1)
                   iorder(8,local_index)=iorder(8,local_index-1)
                 ENDIF          ! IF(idout(1).NE.idout(2))THEN


c The following channel merges two possibilities: the four quarks 
c denoted by 'q' may be associated to two W or two Z bosons. 
c If the process can be both of type 2gZWW and 2g3Z, calling two 
c different channels of type qq -> {gq} {gqZ} is redundant because the 
c difference in the ordering is irrelevant. As in this case the channel
c qq -> {gq} {gqZ} would be evaluated twice (it fulfills two different 
c matchings), the flag ifound_phsp is necessary to avoid double counting
c on the variable n_kphs.

                 IF(ifound_phsp.EQ.0)THEN
                   ifound_phsp=1

*     WRITE(*,*) i,'qq -> {gq} {gqZ}  : 24'
                   local_index=24
                   ialfa(local_index)=1
                   n_kphs=n_kphs+1

                   iorder(3,local_index)=1 ! g
                   iorder(7,local_index)=2 ! g
                   iorder(5,local_index)=7 ! Z
                   iorder(6,local_index)=8
                   IF(ipartinitial(iswap(1)).EQ.-1)THEN
                     iorder(iswap(1),local_index)=4 ! Z
                     iorder(4,local_index)=3
                   ELSE
                     iorder(iswap(1),local_index)=3 ! Z
                     iorder(4,local_index)=4
                   ENDIF
                   IF(ipartinitial(iswap(2)).EQ.-1)THEN
                     iorder(iswap(2),local_index)=6 ! Z
                     iorder(8,local_index)=5
                   ELSE
                     iorder(iswap(2),local_index)=5 ! Z
                     iorder(8,local_index)=6
                   ENDIF

                   idplocal(iorder(3,local_index))=21 ! g
                   idplocal(iorder(7,local_index))=21 ! g
                   idplocal(iorder(5,local_index))=idz0(3,i) ! Z
                   idplocal(iorder(6,local_index))=-idz0(3,i)
                   IF(ipartinitial(iswap(1)).EQ.-1)THEN
                     idplocal(iorder(iswap(1),local_index))=-idz0(1,i)!Z
                     idplocal(iorder(4,local_index))=idz0(1,i)
                   ELSE
                     idplocal(iorder(iswap(1),local_index))=idz0(1,i) !Z
                     idplocal(iorder(4,local_index))=-idz0(1,i)
                   ENDIF
                   IF(ipartinitial(iswap(2)).EQ.-1)THEN
                     idplocal(iorder(iswap(2),local_index))=-idz0(2,i)!Z
                     idplocal(iorder(8,local_index))=idz0(2,i)
                   ELSE
                     idplocal(iorder(iswap(2),local_index))=idz0(2,i) !Z
                     idplocal(iorder(8,local_index))=-idz0(2,i)
                   ENDIF

                   IF(n_kphs.eq.1)THEN
                     iabase=local_index
                     igr(6)=1
                     DO j=1,8
                       igo(j,6)=j
                     ENDDO
                     DO j=1,8
                       idp(j)=idplocal(j)
                     ENDDO
                     DO j=1,2
                       idp_inout(iorder(j,iabase))=-1
                     ENDDO
                   ELSE
                     ial=local_index
                     CALL REORDER(iorder(1,ial),iorder(1,iabase),
     &                    npart+1,idplocal,idp,idp_inout,ipp,ipa)
                     IF(igr(6).EQ.0)THEN
                       igr(6)=1
                       DO j=1,4
                         igo(2*j-1,6)=2*index(j,ipp)-1
                         igo(2*j,6)=2*index(j,ipa)
                       ENDDO
                     ENDIF
                   ENDIF
                 ENDIF          ! IF(ifound_phsp.EQ.0)THEN
               ENDIF            ! IF(abs(iproc(1)).LT.11 .AND. ...
c end giuseppe ILC 25/06/2007

               IF(nz0flep(i).EQ.0)THEN
*     WRITE(*,*) i,' t(00) {ggZ}  : 32'
                 local_index=32
                 ialfa(local_index)=1
                 n_kphs=n_kphs+1

                 IF(ipartinitial(1).EQ.-1)THEN
                   iorder(1,local_index)=4 ! Z
                   iorder(7,local_index)=3
                 ELSE
                   iorder(1,local_index)=3 ! Z
                   iorder(7,local_index)=4
                 ENDIF
                 IF(ipartinitial(2).EQ.-1)THEN
                   iorder(2,local_index)=6 ! Z
                   iorder(8,local_index)=5
                 ELSE
                   iorder(2,local_index)=5 ! Z
                   iorder(8,local_index)=6
                 ENDIF
                 iorder(3,local_index)=7 ! Z
                 iorder(4,local_index)=8
                 iorder(5,local_index)=1 ! g
                 iorder(6,local_index)=2 ! g

                 IF(ipartinitial(1).EQ.-1)THEN
                   idplocal(iorder(1,local_index))=-idz0(iswap(1),i) ! Z
                   idplocal(iorder(7,local_index))=idz0(iswap(1),i)
                 ELSE
                   idplocal(iorder(1,local_index))=idz0(iswap(1),i) ! Z
                   idplocal(iorder(7,local_index))=-idz0(iswap(1),i)
                 ENDIF
                 IF(ipartinitial(2).EQ.-1)THEN
                   idplocal(iorder(2,local_index))=-idz0(iswap(2),i) ! Z
                   idplocal(iorder(8,local_index))=idz0(iswap(2),i)
                 ELSE
                   idplocal(iorder(2,local_index))=idz0(iswap(2),i) ! Z
                   idplocal(iorder(8,local_index))=-idz0(iswap(2),i)
                 ENDIF
                 idplocal(iorder(3,local_index))=idz0(3,i) ! Z
                 idplocal(iorder(4,local_index))=-idz0(3,i)
                 idplocal(iorder(5,local_index))=21 ! g
                 idplocal(iorder(6,local_index))=21 ! g

                 IF(n_kphs.eq.1)THEN
                   iabase=local_index
                   igr(6)=1
                   DO j=1,8
                     igo(j,6)=j
                   ENDDO
                   DO j=1,8
                     idp(j)=idplocal(j)
                   ENDDO
                   DO j=1,2
                     idp_inout(iorder(j,iabase))=-1
                   ENDDO
                 ELSE
                   ial=local_index
                   CALL REORDER(iorder(1,ial),iorder(1,iabase),
     &                  npart+1,idplocal,idp,idp_inout,ipp,ipa)
                   IF(igr(6).EQ.0)THEN
                     igr(6)=1
                     DO j=1,4
                       igo(2*j-1,6)=2*index(j,ipp)-1
                       igo(2*j,6)=2*index(j,ipa)
                     ENDDO
                   ENDIF
                 ENDIF
               ENDIF            ! IF(nz0flep(i).EQ.0)THEN

             ENDIF              ! IF(nwm(i).EQ.1 .AND. nwp(i).EQ.1)THEN

           ELSE
             print*,'***proc.f : ERROR in qq -> gg X !!!'
             print*,'   Uncorrect matching in nfinal(i).eq.1'
             STOP

           ENDIF                ! IF((nwmf(i)+nwpf(i)).EQ.1)THEN

         ENDIF                  ! IF(nfinal(i).EQ.2)THEN

c Check for top
***********************************************************************
* W- -> {ggtbar} -> gg W-   
*            and  W+ -> {ggt} -> gg W+            ialfa(69) 1_5to1_4to31
* W- -> g {gtbar} -> g {gW-} 
*            and  W+ -> g {gt} -> g {gW+}         ialfa(70) 1_5to1_4to31
* W- -> gg tbar -> {ggW-}   
*            and  W+ -> gg t -> {ggW+}            ialfa(71) 1_5to1_4to31
* t(W ) {ggtbar} -> t(W ) gg W-                              
*            and  t(W ) {ggt} -> t(W ) gg W+      ialfa(72) 1_5to1_4to31
* t(W ) g {gtbar} -> t(W ) g {gW-}  
*            and  t(W ) g {gt} -> t(W ) g {gW+}   ialfa(73) 1_5to1_4to31
* t(W ) gg tbar -> t(W ) {ggW-} 
*            and  t(W ) gg t -> t(W ) {ggW+}      ialfa(74) 1_5to1_4to31
***********************************************************************
c gg t in final state

c Case 1: "W+ -> gg t" (s-channel process)
          IF(nwpf(i).EQ.1 .AND. nb.EQ.1 .AND. nantib.EQ.1)THEN
*           WRITE(*,*) i,' W+ -> {ggt} -> gg W+  : 69'
            local_index=69
            ialfa(local_index)=1
            n_kphs=n_kphs+1

            IF(ipartinitial(1).EQ.-1)THEN
              iorder(1,local_index)=8 ! W-
              iorder(2,local_index)=7
            ELSE
              iorder(1,local_index)=7 ! W-
              iorder(2,local_index)=8
            ENDIF
            iorder(3,local_index)=5 ! W+
            iorder(4,local_index)=6
            iorder(5,local_index)=3 ! b
            iorder(6,local_index)=1 ! g
            iorder(7,local_index)=2 ! g
            iorder(8,local_index)=4 ! bbar

            IF(ipartinitial(1).EQ.-1)THEN
              idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
              idplocal(iorder(2,local_index))=idwm(1,i)
            ELSE
              idplocal(iorder(1,local_index))=idwm(1,i) ! W-
              idplocal(iorder(2,local_index))=-(idwm(1,i)+1)
            ENDIF
            idplocal(iorder(3,local_index))=idwp(1,i) ! W+
            idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
            idplocal(iorder(5,local_index))=idz0(1,i) ! b
            idplocal(iorder(6,local_index))=21 ! g
            idplocal(iorder(7,local_index))=21 ! g
            idplocal(iorder(8,local_index))=-idz0(1,i) ! bbar

            IF(n_kphs.eq.1)THEN
              iabase=local_index
              igr(5)=1
              DO j=1,8
                igo(j,5)=j
              ENDDO
              DO j=1,8
                idp(j)=idplocal(j)
              ENDDO
              DO j=1,2
                idp_inout(iorder(j,iabase))=-1
              ENDDO
            ELSE
              ial=local_index
              CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &             idplocal,idp,idp_inout,ipp,ipa)
              IF(igr(5).EQ.0)THEN
                igr(5)=1
                DO j=1,4
                  igo(2*j-1,5)=2*index(j,ipp)-1
                  igo(2*j,5)=2*index(j,ipa)
                ENDDO
              ENDIF
            ENDIF


            IF(nwpflep(i).EQ.0)THEN

c Channel ialfa(70) differs from ialfa(69) by the simple exchange 
c g <-> b in the phase space 1_5to1_4to31.

*             WRITE(*,*) i,' W+ -> g {gt} -> g {gW+}  : 70'
              local_index=70
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=iorder(1,local_index-1) ! W-
              iorder(2,local_index)=iorder(2,local_index-1)
              iorder(3,local_index)=iorder(3,local_index-1) ! W+
              iorder(4,local_index)=iorder(4,local_index-1)
              iorder(5,local_index)=iorder(6,local_index-1) ! g
              iorder(6,local_index)=iorder(5,local_index-1) ! b
              iorder(7,local_index)=iorder(7,local_index-1) ! g
              iorder(8,local_index)=iorder(8,local_index-1) ! bbar


c Channel ialfa(71) differs from ialfa(69) by the simple exchange 
c g1 <-> b  et b <-> g2 [cyclic permutation] in the phase space 
c 1_5to1_4to31.

*             WRITE(*,*) i,' W+ -> gg t -> {ggW+}  : 71'
              local_index=71
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=iorder(1,local_index-2) ! W-
              iorder(2,local_index)=iorder(2,local_index-2)
              iorder(3,local_index)=iorder(3,local_index-2) ! W+
              iorder(4,local_index)=iorder(4,local_index-2)
              iorder(5,local_index)=iorder(6,local_index-2) ! g
              iorder(6,local_index)=iorder(7,local_index-2) ! g
              iorder(7,local_index)=iorder(5,local_index-2) ! b
              iorder(8,local_index)=iorder(8,local_index-2) ! bbar

            ENDIF               ! IF(nwpflep(i).EQ.0)THEN

          ENDIF ! IF(nwpf(i).EQ.1 .AND. nb.EQ.1 .AND. nantib.EQ.1)THEN


c Case 2:  "b t(W ) -> gg t" (t-channel process)
          IF(nwpf(i).EQ.1 .AND. (nantib-nb).EQ.-1)THEN 
*           WRITE(*,*) i,' t(W ) {ggt} -> t(W ) gg W+  : 72'
            local_index=72
            ialfa(local_index)=1
            n_kphs=n_kphs+1

            IF(abs(idout(1)).EQ.5)THEN
              iorder(1,local_index)=4 ! bbar
              IF(ipartinitial(2).EQ.-1)THEN
                iorder(2,local_index)=8 ! W-
                iorder(8,local_index)=7
              ELSE
                iorder(2,local_index)=7 ! W-
                iorder(8,local_index)=8
              ENDIF
            ELSE                ! abs(idout(2)).EQ.5
              iorder(2,local_index)=4 ! bbar
              IF(ipartinitial(1).EQ.-1)THEN
                iorder(1,local_index)=8 ! W-
                iorder(8,local_index)=7
              ELSE
                iorder(1,local_index)=7 ! W-
                iorder(8,local_index)=8
              ENDIF
            ENDIF
            iorder(3,local_index)=5 ! W+
            iorder(4,local_index)=6
            iorder(5,local_index)=3 ! b
            iorder(6,local_index)=1 ! g
            iorder(7,local_index)=2 ! g

            IF(abs(idout(1)).EQ.5)THEN
              idplocal(iorder(1,local_index))=-idz0(1,i) ! bbar
              IF(ipartinitial(2).EQ.-1)THEN
                idplocal(iorder(2,local_index))=-(idwm(1,i)+1) ! W-
                idplocal(iorder(8,local_index))=idwm(1,i)
              ELSE
                idplocal(iorder(2,local_index))=idwm(1,i) ! W-
                idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
              ENDIF
            ELSE                ! abs(idout(2)).EQ.5
              idplocal(iorder(2,local_index))=-idz0(1,i) ! bbar
              IF(ipartinitial(1).EQ.-1)THEN
                idplocal(iorder(1,local_index))=-(idwm(1,i)+1) ! W-
                idplocal(iorder(8,local_index))=idwm(1,i)
              ELSE
                idplocal(iorder(1,local_index))=idwm(1,i) ! W-
                idplocal(iorder(8,local_index))=-(idwm(1,i)+1)
              ENDIF
            ENDIF
            idplocal(iorder(3,local_index))=idwp(1,i) ! W+
            idplocal(iorder(4,local_index))=-(idwp(1,i)-1)
            idplocal(iorder(5,local_index))=idz0(1,i) ! b
            idplocal(iorder(6,local_index))=21 ! g
            idplocal(iorder(7,local_index))=21 ! g

            IF(n_kphs.eq.1)THEN
              iabase=local_index
              igr(5)=1
              DO j=1,8
                igo(j,5)=j
              ENDDO
              DO j=1,8
                idp(j)=idplocal(j)
              ENDDO
              DO j=1,2
                idp_inout(iorder(j,iabase))=-1
              ENDDO
            ELSE
              ial=local_index
              CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &             idplocal,idp,idp_inout,ipp,ipa)
              IF(igr(5).EQ.0)THEN
                igr(5)=1
                DO j=1,4
                  igo(2*j-1,5)=2*index(j,ipp)-1
                  igo(2*j,5)=2*index(j,ipa)
                ENDDO
              ENDIF
            ENDIF


            IF(nwpflep(i).EQ.0)THEN

c Channel ialfa(73) differs from ialfa(72) by the simple exchange
c g <-> b in the phase space 1_5to1_4to31.

*             WRITE(*,*) i,' t(W ) g {gt} -> t(W ) g {gW+}  : 73'
              local_index=73
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=iorder(1,local_index-1)
              iorder(2,local_index)=iorder(2,local_index-1)
              iorder(3,local_index)=iorder(3,local_index-1) ! W+
              iorder(4,local_index)=iorder(4,local_index-1)
              iorder(5,local_index)=iorder(6,local_index-1) ! g
              iorder(6,local_index)=iorder(5,local_index-1) ! b
              iorder(7,local_index)=iorder(7,local_index-1) ! g
              iorder(8,local_index)=iorder(8,local_index-1)


c Channel ialfa(74) differs from ialfa(72) by the simple exchange
c g1 <-> b  et b <-> g2 [cyclic permutation] in the phase space 
c 1_5to1_4to31.

*             WRITE(*,*) i,' t(W ) gg t -> t(W ) {ggW+}  : 74'
              local_index=74
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=iorder(1,local_index-2)
              iorder(2,local_index)=iorder(2,local_index-2)
              iorder(3,local_index)=iorder(3,local_index-2) ! W+
              iorder(4,local_index)=iorder(4,local_index-2)
              iorder(5,local_index)=iorder(6,local_index-2) ! g
              iorder(6,local_index)=iorder(7,local_index-2) ! g
              iorder(7,local_index)=iorder(5,local_index-2) ! b
              iorder(8,local_index)=iorder(8,local_index-2)

            ENDIF               ! IF(nwpflep(i).EQ.0)THEN

          ENDIF ! IF(nwpf(i).EQ.1 .AND. (nantib-nb).EQ.-1)THEN 


c gg tbar in final state

c Case 1: "W- -> gg tbar"
          IF(nwmf(i).EQ.1 .AND. nb.EQ.1 .AND. nantib.EQ.1)THEN
*           WRITE(*,*) i,' W- -> {ggtbar} -> gg W-  : 69'
            local_index=69
            ialfa(local_index)=1
            n_kphs=n_kphs+1

            IF(ipartinitial(1).EQ.-1)THEN
              iorder(1,local_index)=6 ! W+
              iorder(2,local_index)=5
            ELSE
              iorder(1,local_index)=5 ! W+
              iorder(2,local_index)=6
            ENDIF
            iorder(3,local_index)=7 ! W-
            iorder(4,local_index)=8
            iorder(5,local_index)=4 ! bbar
            iorder(6,local_index)=1 ! g
            iorder(7,local_index)=2 ! g
            iorder(8,local_index)=3 ! b

            IF(ipartinitial(1).EQ.-1)THEN
              idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
              idplocal(iorder(2,local_index))=idwp(1,i)
            ELSE
              idplocal(iorder(1,local_index))=idwp(1,i) ! W+
              idplocal(iorder(2,local_index))=-(idwp(1,i)-1)
            ENDIF
            idplocal(iorder(3,local_index))=idwm(1,i) ! W-
            idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
            idplocal(iorder(5,local_index))=-idz0(1,i) ! bbar
            idplocal(iorder(6,local_index))=21 ! g
            idplocal(iorder(7,local_index))=21 ! g
            idplocal(iorder(8,local_index))=idz0(1,i) ! b

            IF(n_kphs.eq.1)THEN
              iabase=local_index
              igr(5)=1
              DO j=1,8
                igo(j,5)=j
              ENDDO
              DO j=1,8
                idp(j)=idplocal(j)
              ENDDO
              DO j=1,2
                idp_inout(iorder(j,iabase))=-1
              ENDDO
            ELSE
              ial=local_index
              CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &             idplocal,idp,idp_inout,ipp,ipa)
              IF(igr(5).EQ.0)THEN
                igr(5)=1
                DO j=1,4
                  igo(2*j-1,5)=2*index(j,ipp)-1
                  igo(2*j,5)=2*index(j,ipa)
                ENDDO
              ENDIF
            ENDIF


            IF(nwmflep(i).EQ.0)THEN

c Channel ialfa(70) differs from ialfa(69) by the simple exchange 
c g <-> bbar in the phase space 1_5to1_4to31.

*             WRITE(*,*) i,' W- -> g {gtbar} -> g {gW-}  : 70'
              local_index=70
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=iorder(1,local_index-1) ! W+
              iorder(2,local_index)=iorder(2,local_index-1)
              iorder(3,local_index)=iorder(3,local_index-1) ! W-
              iorder(4,local_index)=iorder(4,local_index-1)
              iorder(5,local_index)=iorder(6,local_index-1) ! g
              iorder(6,local_index)=iorder(5,local_index-1) ! bbar
              iorder(7,local_index)=iorder(7,local_index-1) ! g
              iorder(8,local_index)=iorder(8,local_index-1) ! b


c Channel ialfa(71) differs from ialfa(69) by the simple exchange 
c g1 <-> bbar  et bbar <-> g2 [cyclic permutation] in the phase space 
c 1_5to1_4to31.

*             WRITE(*,*) i,' W- -> gg tbar -> {ggW-}  : 71'
              local_index=71
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=iorder(1,local_index-2) ! W+
              iorder(2,local_index)=iorder(2,local_index-2)
              iorder(3,local_index)=iorder(3,local_index-2) ! W-
              iorder(4,local_index)=iorder(4,local_index-2)
              iorder(5,local_index)=iorder(6,local_index-2) ! g
              iorder(6,local_index)=iorder(7,local_index-2) ! g
              iorder(7,local_index)=iorder(5,local_index-2) ! bbar
              iorder(8,local_index)=iorder(8,local_index-2) ! b

            ENDIF               ! IF(nwmflep(i).EQ.0)THEN

          ENDIF ! IF(nwmf(i).EQ.1 .AND. nb.EQ.1 .AND. nantib.EQ.1)THEN


c Case 2: "bbar t(W ) -> gg tbar"
          IF(nwmf(i).EQ.1 .AND. (nantib-nb).EQ.1)THEN
*           WRITE(*,*) i,' t(W ) {ggtbar} -> t(W ) gg W-  : 72'
            local_index=72
            ialfa(local_index)=1
            n_kphs=n_kphs+1

            IF(abs(idout(1)).EQ.5)THEN
              iorder(1,local_index)=3 ! b
              IF(ipartinitial(2).EQ.-1)THEN
                iorder(2,local_index)=6 ! W+
                iorder(8,local_index)=5
              ELSE
                iorder(2,local_index)=5 ! W+
                iorder(8,local_index)=6
              ENDIF
            ELSE                ! abs(idout(2)).EQ.5
              iorder(2,local_index)=3 ! b
              IF(ipartinitial(1).EQ.-1)THEN
                iorder(1,local_index)=6 ! W+
                iorder(8,local_index)=5
              ELSE
                iorder(1,local_index)=5 ! W+
                iorder(8,local_index)=6
              ENDIF
            ENDIF
            iorder(3,local_index)=7 ! W-
            iorder(4,local_index)=8
            iorder(5,local_index)=4 ! bbar
            iorder(6,local_index)=1 ! g
            iorder(7,local_index)=2 ! g

            IF(abs(idout(1)).EQ.5)THEN
              idplocal(iorder(1,local_index))=idz0(1,i) ! b
              IF(ipartinitial(2).EQ.-1)THEN
                idplocal(iorder(2,local_index))=-(idwp(1,i)-1) ! W+
                idplocal(iorder(8,local_index))=idwp(1,i)
              ELSE
                idplocal(iorder(2,local_index))=idwp(1,i) ! W+
                idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
              ENDIF
            ELSE                ! abs(idout(2)).EQ.5
              idplocal(iorder(2,local_index))=idz0(1,i) ! b
              IF(ipartinitial(1).EQ.-1)THEN
                idplocal(iorder(1,local_index))=-(idwp(1,i)-1) ! W+
                idplocal(iorder(8,local_index))=idwp(1,i)
              ELSE
                idplocal(iorder(1,local_index))=idwp(1,i) ! W+
                idplocal(iorder(8,local_index))=-(idwp(1,i)-1)
              ENDIF
            ENDIF
            idplocal(iorder(3,local_index))=idwm(1,i) ! W-
            idplocal(iorder(4,local_index))=-(idwm(1,i)+1)
            idplocal(iorder(5,local_index))=-idz0(1,i) ! bbar
            idplocal(iorder(6,local_index))=21 ! g
            idplocal(iorder(7,local_index))=21 ! g

            IF(n_kphs.eq.1)THEN
              iabase=local_index
              igr(5)=1
              DO j=1,8
                igo(j,5)=j
              ENDDO
              DO j=1,8
                idp(j)=idplocal(j)
              ENDDO
              DO j=1,2
                idp_inout(iorder(j,iabase))=-1
              ENDDO
            ELSE
              ial=local_index
              CALL REORDER(iorder(1,ial),iorder(1,iabase),npart+1,
     &             idplocal,idp,idp_inout,ipp,ipa)
              IF(igr(5).EQ.0)THEN
                igr(5)=1
                DO j=1,4
                  igo(2*j-1,5)=2*index(j,ipp)-1
                  igo(2*j,5)=2*index(j,ipa)
                ENDDO
              ENDIF
            ENDIF


            IF(nwmflep(i).EQ.0)THEN

c Channel ialfa(73) differs from ialfa(72) by the simple exchange
c g <-> bbar in the phase space 1_5to1_4to31.

*             WRITE(*,*) i,' t(W ) g {gtbar} -> t(W ) g {gW-}  : 73'
              local_index=73
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=iorder(1,local_index-1)
              iorder(2,local_index)=iorder(2,local_index-1)
              iorder(3,local_index)=iorder(3,local_index-1) ! W-
              iorder(4,local_index)=iorder(4,local_index-1)
              iorder(5,local_index)=iorder(6,local_index-1) ! g
              iorder(6,local_index)=iorder(5,local_index-1) ! bbar
              iorder(7,local_index)=iorder(7,local_index-1) ! g
              iorder(8,local_index)=iorder(8,local_index-1)


c Channel ialfa(74) differs from ialfa(72) by the simple exchange 
c g1 <-> bbar  et bbar <-> g2 [cyclic permutation] in the phase space 
c 1_5to1_4to31.

*             WRITE(*,*) i,' t(W ) gg tbar -> t(W ) {ggW-}  : 74'
              local_index=74
              ialfa(local_index)=1
              n_kphs=n_kphs+1

              iorder(1,local_index)=iorder(1,local_index-2)
              iorder(2,local_index)=iorder(2,local_index-2)
              iorder(3,local_index)=iorder(3,local_index-2) ! W-
              iorder(4,local_index)=iorder(4,local_index-2)
              iorder(5,local_index)=iorder(6,local_index-2) ! g
              iorder(6,local_index)=iorder(7,local_index-2) ! g
              iorder(7,local_index)=iorder(5,local_index-2) ! bbar
              iorder(8,local_index)=iorder(8,local_index-2)

            ENDIF               ! IF(nwmflep(i).EQ.0)THEN

          ENDIF ! IF(nwmf(i).EQ.1 .AND. (nantib-nb).EQ.1)THEN
c end Check for top

        ENDDO                   ! DO i=1,nmatch

      ENDIF                     ! IF(ngluon.EQ.0)THEN



      
c       PRINT*,'      idp:',idp
c       PRINT*,'idp_inout:',idp_inout
c       PRINT*,'iabase',iabase

c      DO i=1,mxnphs
c        IF(ialfa(i).eq.1)THEN
c          PRINT*,i,': iorder:',(iorder(j,i),j=1,8)
c        ENDIF
c      ENDDO
c
c      DO i=1,6
c        IF(igr(i).EQ.1)THEN
c          PRINT*,'igr',i,'igo',(igo(j,i),j=1,8)
c        ENDIF
c      ENDDO



      RETURN
      END





      SUBROUTINE REORDER(iorder,iorderbase,npart,idplocal,
     &                   idp,idp_inout,ipp,ipa)
      DIMENSION ifact(0:4),index(4,24)
      DATA ifact/1,1,2,6,24/
      DATA index/
     &     1,2,3,4,               !even
     &     2,1,3,4,
     &     2,3,1,4,               !even
     &     1,3,2,4,
     &     3,1,2,4,               !even
     &     3,2,1,4,
     
     &     3,2,4,1,               !even
     &     3,1,4,2,               
     &     1,3,4,2,               !even
     &     2,3,4,1,               
     &     2,1,4,3,               !even
     &     1,2,4,3,               

     &     1,4,2,3,               !even
     &     2,4,1,3,
     &     2,4,3,1,               !even
     &     1,4,3,2,
     &     3,4,1,2,               !even
     &     3,4,2,1,

     &     4,3,2,1,               !even
     &     4,3,1,2,               
     &     4,1,3,2,               !even
     &     4,2,3,1,               
     &     4,2,1,3,               !even
     &     4,1,2,3               
     &         /

       DIMENSION iorder(8),iorderbase(8),idp(8),idp_inout(8),
     &           idpartamp(4),idantipartamp(4),
     &           ipartampinout(4),iantipartampinout(4),
     &           idplocal(8),idp_iolocal(8),
     &           idpartlocal(4),ipartiolocal(4),
     &           idantipartlocal(4),iantipartiolocal(4),
     &           invorder(8)

      DO j=1,8
        IF(mod(iorderbase(j),2).EQ.1)THEN
          idpartamp(iorderbase(j)/2+1)=
     &                idp(iorderbase(j))
          ipartampinout(iorderbase(j)/2+1)=
     &                idp_inout(iorderbase(j))
        ELSE
          idantipartamp(iorderbase(j)/2)=
     &                idp(iorderbase(j))
          iantipartampinout(iorderbase(j)/2)=
     &                idp_inout(iorderbase(j))
        ENDIF
      ENDDO
      
      DO j=1,8
        invorder(iorder(j))=j
      ENDDO
      
      DO j=1,2
        idp_iolocal(iorder(j))=-1
      ENDDO
      DO j=3,8
        idp_iolocal(iorder(j))=1
      ENDDO
      DO j=1,8
        IF(mod(iorder(j),2).EQ.1)THEN
          idpartlocal(iorder(j)/2+1)=idplocal(iorder(j))
          ipartiolocal(iorder(j)/2+1)=idp_iolocal(iorder(j))
        ELSE
          idantipartlocal(iorder(j)/2)=idplocal(iorder(j))
          iantipartiolocal(iorder(j)/2)=idp_iolocal(iorder(j))
        ENDIF
      ENDDO
      
c      print*,'idpartamp',idpartamp
c      print*,'ipartampinout',ipartampinout
c      print*,'idantipartamp',idantipartamp
c      print*,'iantipartampinout',iantipartampinout
c      print*,'idplocal',idplocal
c      print*,'idp_iolocal',idp_iolocal
c      print*,'idpartlocal',idpartlocal
c      print*,'ipartiolocal',ipartiolocal
c      print*,'idantipartlocal',idantipartlocal
c      print*,'iantipartiolocal',iantipartiolocal

      ipp=1
      ifound=0
      DO WHILE(ipp.LE.ifact(npart) .AND.
     &           ifound.EQ.0)
        imatch=0
        DO j=1,npart
          imatch=imatch
     &      +abs(idpartamp(index(j,ipp))-idpartlocal(j))
     &      +abs(ipartampinout(index(j,ipp))-ipartiolocal(j))
        ENDDO
        IF(imatch.EQ.0)THEN
          ifound=1
        ELSE
          ipp=ipp+1
        ENDIF
      ENDDO
      
c      print*,'ip part:',ipp

      DO j=1,4
cc        iorder(2*j-1)=iorderbase(2*invindex(j,ip)-1)
cc        iorder(2*invindex(j,ip)-1)=iorderbase(2*j-1)
cc        iorder(2*invindex(j,ip)-1)=2*j-1
cc        iorder(2*j-1)=2*invindex(j,ip)-1
        iorder(invorder(2*j-1))=2*index(j,ipp)-1
      ENDDO
      ipa=1
      ifound=0
      DO WHILE(ipa.LE.ifact(npart) .AND.
     &           ifound.EQ.0)
        imatch=0
        DO j=1,npart
          imatch=imatch
     &    +abs(idantipartamp(index(j,ipa))-idantipartlocal(j))
     &    +abs(iantipartampinout(index(j,ipa))
     &         -iantipartiolocal(j))
        ENDDO
        IF(imatch.EQ.0)THEN
          ifound=1
        ELSE
          ipa=ipa+1
        ENDIF
      ENDDO

c      print*,'ip antipart:',ipa

      DO j=1,4
cc        iorder(2*j)=iorderbase(2*invindex(j,ip))
cc        iorder(2*invindex(j,ip))=iorderbase(2*j)
cc        iorder(2*invindex(j,ip))=2*j
cc        iorder(2*j)=2*invindex(j,ip)
        iorder(invorder(2*j))=2*index(j,ipa)
      ENDDO
      
c      print*,'ip antipart:',ipa

      RETURN
      END
      
