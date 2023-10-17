*
* procextraini.f
*
* Last update: Oct 30, 2007
*
*
************************************************************************
*NOTICE: using NEW versions of phase space routines 1_1_3, 1_2_3, 3_3:
*
*    phsp                PHASE       --->         PHANTOM
*_______________________________________________________________________
*
*   1_1_3            rm1,rm11,rm111     rm1[a,b],rm11[a,b],rm111[a,b]
*
*   1_2_3               rm1,rm11         rm0[a,b],rm1[a,b],rm11[a,b]
*
*    3_3                rm1,rm2               rm1[a,b],rm2[a,b]
*
************************************************************************


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
************************************************************************



      SUBROUTINE procextraini

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      implicit double complex (c)


      INCLUDE 'common.h'

      INCLUDE 'common_subproc.h'
      COMMON/phaones/ionesh
* the following variables represent the input of the various phase
* spaces
* rmN=invariant mass of the N-fork
* gamN=width of the N-fork
* spmNlo=lower bound of the resonance in the N-fork
* spmNhi=upper bound of the resonance in the N-fork
* in the multi-mapping strategy,
* [spmNlo,spmNhi] defines the range where the mapping on the resonance
* is on
* rmtN=boson mass in the "N" t-channel mapping
* expN=value of the exponent used in the "N" t-channel mapping
* sptN=separaction between mapped and flat regions in the "N" t-channel
* see the specific phase space for details.

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


      COMMON/phspinput_data/spmhlo,spmhhi,spmhhlo,spmhhhi,
     &                                   spmzlo,spmzhi,spmwlo,spmwhi

* Information from proc.f
      COMMON/info_proc/
c number of leptonic final state bosons for each channel. This
c information is loaded in proc.f only for those channels which do
c require it in procextraini.
     &     num_wpflep(mxnphs),num_wmflep(mxnphs),num_z0flep(mxnphs)


*** define phase space inputs which are process dependent

      n5f=0
      DO j=3,8
        IF(ABS(iproc(j)).EQ.5)n5f=n5f+1
      ENDDO !j
      n5=n5f
      IF(ABS(iproc(1)).EQ.5)n5=n5-1
      IF(ABS(iproc(2)).EQ.5)n5=n5-1
c giuseppe 30/10/2007
      n5tot=0
      DO j=1,8
        IF(ABS(iproc(j)).EQ.5)n5tot=n5tot+1
      ENDDO
c n5_t is the number of possible t-channel-like b quark lines
      n5_t=n5f-n5
      n5i=n5f-n5
c end giuseppe 30/10/2007

c Ezio 2018_08_02
c Nel caso dello spazio fasi 6 dvo sapere se ci sono due gluoni nello stato
c finale oppure no
      ngf=0
      DO j=3,8
        IF(iproc(j).EQ.21)ngf=ngf+1
      ENDDO !j
c end Ezio 2018_08_02

c giuseppe 30/10/2007
* t(00) W+ W-                                     ialfa(1)  1_1_4
      local_index=1
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007


* t(0W) Z W-  and t(0W) Z W+                 ialfa(2)  1_1_4
      local_index=2
      IF(n5.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
      ENDIF
      IF(n5.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
      ELSE
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
      ENDIF

* t(W0) Z W-  and t(W0) Z W+                 ialfa(3)  1_1_4
      local_index=3
      IF(n5.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
      ENDIF
      IF(n5.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
      ELSE
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
      ENDIF

c giuseppe 30/10/2007
* t(WW) W+ W-                                     ialfa(5)  1_1_4
      local_index=5
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007

* t(00) Z Z                                       ialfa(6)  1_1_4
      local_index=6
c Ezio 2018_08_02
       IF(ngf.EQ.0) THEN
c giuseppe 30/10/2007
        IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
          rm1a(local_index)=0.d0
          gam1a(local_index)=0.d0
          spm1alo(local_index)=0.d0
          spm1ahi(local_index)=0.d0
        ELSE
          rm1a(local_index)=rmh
          gam1a(local_index)=gamh
          spm1alo(local_index)=spmhlo
          spm1ahi(local_index)=spmhhi
        ENDIF
        IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
          rm1b(local_index)=0.d0
          gam1b(local_index)=0.d0
          spm1blo(local_index)=ecoll
          spm1bhi(local_index)=ecoll
        ELSE
          rm1b(local_index)=rmhh
          gam1b(local_index)=gamhh
          spm1blo(local_index)=spmhhlo
          spm1bhi(local_index)=spmhhhi
        ENDIF
c end giuseppe 30/10/2007
        IF(n5.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
          rm11a(local_index)=0.d0
          gam11a(local_index)=0.d0
          spm11alo(local_index)=0.d0
          spm11ahi(local_index)=0.d0
          rm11b(local_index)=rmz
          gam11b(local_index)=gamz
          spm11blo(local_index)=spmzlo
          spm11bhi(local_index)=spmzhi
          rm12a(local_index)=0.d0
          gam12a(local_index)=0.d0
          spm12alo(local_index)=0.d0
          spm12ahi(local_index)=0.d0
          rm12b(local_index)=rmz
          gam12b(local_index)=gamz
          spm12blo(local_index)=spmzlo
          spm12bhi(local_index)=spmzhi
        ELSE
          rm11a(local_index)=rmz
          gam11a(local_index)=gamz
          spm11alo(local_index)=spmzlo
          spm11ahi(local_index)=spmzhi
          rm11b(local_index)=rmh
          gam11b(local_index)=gamh
          spm11blo(local_index)=spmhlo
          spm11bhi(local_index)=spmhhi
          rm12a(local_index)=rmz
          gam12a(local_index)=gamz
          spm12alo(local_index)=spmzlo
          spm12ahi(local_index)=spmzhi
          rm12b(local_index)=rmh
          gam12b(local_index)=gamh
          spm12blo(local_index)=spmhlo
          spm12bhi(local_index)=spmhhi
        ENDIF
        IF(n5.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
          rm11c(local_index)=0.d0
          gam11c(local_index)=0.d0
          spm11clo(local_index)=ecoll
          spm11chi(local_index)=ecoll
          rm12c(local_index)=0.d0
          gam12c(local_index)=0.d0
          spm12clo(local_index)=ecoll
          spm12chi(local_index)=ecoll
        ELSE
          rm11c(local_index)=rmhh
          gam11c(local_index)=gamhh
          spm11clo(local_index)=spmhhlo
          spm11chi(local_index)=spmhhhi
          rm12c(local_index)=rmhh
          gam12c(local_index)=gamhh
          spm12clo(local_index)=spmhhlo
          spm12chi(local_index)=spmhhhi
        ENDIF
      ELSEIF(ngf.EQ.2) THEN
        IF(rmh.GE.1.d3.OR.rmh.LT.0.d0.OR.n5i.LT.2)THEN
          rm1a(local_index)=0.d0
          gam1a(local_index)=0.d0
          spm1alo(local_index)=0.d0
          spm1ahi(local_index)=0.d0
        ELSE
          rm1a(local_index)=rmh
          gam1a(local_index)=gamh
          spm1alo(local_index)=spmhlo
          spm1ahi(local_index)=spmhhi
        ENDIF
        IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0.OR.n5i.LT.2)THEN
          rm1b(local_index)=0.d0
          gam1b(local_index)=0.d0
          spm1blo(local_index)=ecoll
          spm1bhi(local_index)=ecoll
        ELSE
          rm1b(local_index)=rmhh
          gam1b(local_index)=gamhh
          spm1blo(local_index)=spmhhlo
          spm1bhi(local_index)=spmhhhi
        ENDIF
c end giuseppe 30/10/2007
        IF(n5.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
          rm11a(local_index)=0.d0
          gam11a(local_index)=0.d0
          spm11alo(local_index)=0.d0
          spm11ahi(local_index)=0.d0
          rm11b(local_index)=rmz
          gam11b(local_index)=gamz
          spm11blo(local_index)=spmzlo
          spm11bhi(local_index)=spmzhi
          rm12a(local_index)=0.d0
          gam12a(local_index)=0.d0
          spm12alo(local_index)=0.d0
          spm12ahi(local_index)=0.d0
          rm12b(local_index)=rmz
          gam12b(local_index)=gamz
          spm12blo(local_index)=spmzlo
          spm12bhi(local_index)=spmzhi
        ELSE
          rm11a(local_index)=rmz
          gam11a(local_index)=gamz
          spm11alo(local_index)=spmzlo
          spm11ahi(local_index)=spmzhi
          rm11b(local_index)=rmh
          gam11b(local_index)=gamh
          spm11blo(local_index)=spmhlo
          spm11bhi(local_index)=spmhhi
          rm12a(local_index)=rmz
          gam12a(local_index)=gamz
          spm12alo(local_index)=spmzlo
          spm12ahi(local_index)=spmzhi
          rm12b(local_index)=rmh
          gam12b(local_index)=gamh
          spm12blo(local_index)=spmhlo
          spm12bhi(local_index)=spmhhi
        ENDIF
        IF(n5.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
          rm11c(local_index)=0.d0
          gam11c(local_index)=0.d0
          spm11clo(local_index)=ecoll
          spm11chi(local_index)=ecoll
          rm12c(local_index)=0.d0
          gam12c(local_index)=0.d0
          spm12clo(local_index)=ecoll
          spm12chi(local_index)=ecoll
        ELSE
          rm11c(local_index)=rmhh
          gam11c(local_index)=gamhh
          spm11clo(local_index)=spmhhlo
          spm11chi(local_index)=spmhhhi
          rm12c(local_index)=rmhh
          gam12c(local_index)=gamhh
          spm12clo(local_index)=spmhhlo
          spm12chi(local_index)=spmhhhi
        ENDIF
      ENDIF  ! end IF(ngf.EQ.0)
ccccc here 
* t(WW) Z Z                                       ialfa(7)  1_1_4
      local_index=7
c giuseppe 30/10/2007
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007
      IF(n5.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm12a(local_index)=rmz
        gam12a(local_index)=gamz
        spm12alo(local_index)=spmzlo
        spm12ahi(local_index)=spmzhi
        rm12b(local_index)=rmh
        gam12b(local_index)=gamh
        spm12blo(local_index)=spmhlo
        spm12bhi(local_index)=spmhhi
      ENDIF
      IF(n5.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
      ELSE
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
        rm12c(local_index)=rmhh
        gam12c(local_index)=gamhh
        spm12clo(local_index)=spmhhlo
        spm12chi(local_index)=spmhhhi
      ENDIF

* t(g0) W+ W-                                     ialfa(8)  1_1_4
      local_index=8
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF

* t(0g) W+ W-                                     ialfa(9)  1_1_4
      local_index=9
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF

* t(gW) Z W-   and  t(gW) Z W+                    ialfa(10) 1_1_4
      local_index=10
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
      ELSE
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
      ENDIF

* t(Wg) Z W-   and  t(Wg) Z W+                    ialfa(11) 1_1_4
      local_index=11
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
      ELSE
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
      ENDIF

* t(g0) Z Z                                       ialfa(12) 1_1_4
      local_index=12
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
      ELSEIF(n5f.LT.2)THEN
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm12a(local_index)=rmz
        gam12a(local_index)=gamz
        spm12alo(local_index)=spmzlo
        spm12ahi(local_index)=spmzhi
        rm12b(local_index)=rmh
        gam12b(local_index)=gamh
        spm12blo(local_index)=spmhlo
        spm12bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
      ELSEIF(n5f.LT.2)THEN
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
        rm12c(local_index)=rmhh
        gam12c(local_index)=gamhh
        spm12clo(local_index)=spmhhlo
        spm12chi(local_index)=spmhhhi
      ENDIF

* t(0g) Z Z                                       ialfa(13) 1_1_4
      local_index=13
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
      ELSEIF(n5f.LT.2)THEN
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm12a(local_index)=rmz
        gam12a(local_index)=gamz
        spm12alo(local_index)=spmzlo
        spm12ahi(local_index)=spmzhi
        rm12b(local_index)=rmh
        gam12b(local_index)=gamh
        spm12blo(local_index)=spmhlo
        spm12bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
      ELSEIF(n5f.LT.2)THEN
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
        rm12c(local_index)=rmhh
        gam12c(local_index)=gamhh
        spm12clo(local_index)=spmhhlo
        spm12chi(local_index)=spmhhhi
      ENDIF


c giuseppe 30/10/2007
* W- -> [W- W+] W-                                ialfa(14) 2_4
      local_index=14
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007

c giuseppe 30/10/2007
* W- -> W- [W+ W-]                                ialfa(15) 2_4
      local_index=15
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007

c giuseppe 30/10/2007
* W+ -> [W+ W-] W+                                ialfa(16) 2_4
      local_index=16
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007

c giuseppe 30/10/2007
* W+ -> W+ [W- W+]                                ialfa(17) 2_4
      local_index=17
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007

* Z  -> Z  W+ W-   and  gg -> Z  W+ W-       ialfa(18) 2_4
      local_index=18
c giuseppe 30/10/2007
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5tot.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5_t.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm2c(local_index)=0.d0
        gam2c(local_index)=0.d0
        spm2clo(local_index)=ecoll
        spm2chi(local_index)=ecoll
      ELSE
        rm2c(local_index)=rmhh
        gam2c(local_index)=gamhh
        spm2clo(local_index)=spmhhlo
        spm2chi(local_index)=spmhhhi
      ENDIF

* W- -> Z  Z  W- and W+ -> Z  Z  W+          ialfa(19) 2_4
      local_index=19
c giuseppe 30/10/2007
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5tot.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5tot.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm12a(local_index)=rmz
        gam12a(local_index)=gamz
        spm12alo(local_index)=spmzlo
        spm12ahi(local_index)=spmzhi
        rm12b(local_index)=rmh
        gam12b(local_index)=gamh
        spm12blo(local_index)=spmhlo
        spm12bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
      ELSE
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
        rm12c(local_index)=rmhh
        gam12c(local_index)=gamhh
        spm12clo(local_index)=spmhhlo
        spm12chi(local_index)=spmhhhi
      ENDIF

* Z  -> [Z Z] Z   and gg -> [Z Z] Z          ialfa(20) 2_4
      local_index=20
c giuseppe 30/10/2007
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5tot.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5tot.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm12a(local_index)=rmz
        gam12a(local_index)=gamz
        spm12alo(local_index)=spmzlo
        spm12ahi(local_index)=spmzhi
        rm12b(local_index)=rmh
        gam12b(local_index)=gamh
        spm12blo(local_index)=spmhlo
        spm12bhi(local_index)=spmhhi
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
        rm2c(local_index)=0.d0
        gam2c(local_index)=0.d0
        spm2clo(local_index)=ecoll
        spm2chi(local_index)=ecoll
      ELSE
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
        rm12c(local_index)=rmhh
        gam12c(local_index)=gamhh
        spm12clo(local_index)=spmhhlo
        spm12chi(local_index)=spmhhhi
        rm2c(local_index)=rmhh
        gam2c(local_index)=gamhh
        spm2clo(local_index)=spmhhlo
        spm2chi(local_index)=spmhhhi
      ENDIF

* Z  -> [Z Z Z]                              ialfa(21) 2_4
      local_index=21
c giuseppe 30/10/2007
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5tot.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5tot.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm12a(local_index)=rmz
        gam12a(local_index)=gamz
        spm12alo(local_index)=spmzlo
        spm12ahi(local_index)=spmzhi
        rm12b(local_index)=rmh
        gam12b(local_index)=gamh
        spm12blo(local_index)=spmhlo
        spm12bhi(local_index)=spmhhi
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
        rm2c(local_index)=0.d0
        gam2c(local_index)=0.d0
        spm2clo(local_index)=ecoll
        spm2chi(local_index)=ecoll
      ELSE
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
        rm12c(local_index)=rmhh
        gam12c(local_index)=gamhh
        spm12clo(local_index)=spmhhlo
        spm12chi(local_index)=spmhhhi
        rm2c(local_index)=rmhh
        gam2c(local_index)=gamhh
        spm2clo(local_index)=spmhhlo
        spm2chi(local_index)=spmhhhi
      ENDIF

* Z  -> Z [Z Z]                              ialfa(22) 2_4
      local_index=22
c giuseppe 30/10/2007
      IF(rmh.GE.1.d3.OR.rmh.LT.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5tot.EQ.0))THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
      ELSE
        rm1a(local_index)=rmh
        gam1a(local_index)=gamh
        spm1alo(local_index)=spmhlo
        spm1ahi(local_index)=spmhhi
      ENDIF
      IF(rmhh.GE.1.d3.OR.rmhh.LE.0.d0 .OR.
     &     ((i_pertorder.EQ.2 .OR. ngluon.GT.0) .AND. n5tot.EQ.0))THEN
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=ecoll
        spm1bhi(local_index)=ecoll
      ELSE
        rm1b(local_index)=rmhh
        gam1b(local_index)=gamhh
        spm1blo(local_index)=spmhhlo
        spm1bhi(local_index)=spmhhhi
      ENDIF
c end giuseppe 30/10/2007
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm12a(local_index)=rmz
        gam12a(local_index)=gamz
        spm12alo(local_index)=spmzlo
        spm12ahi(local_index)=spmzhi
        rm12b(local_index)=rmh
        gam12b(local_index)=gamh
        spm12blo(local_index)=spmhlo
        spm12bhi(local_index)=spmhhi
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm11c(local_index)=0.d0
        gam11c(local_index)=0.d0
        spm11clo(local_index)=ecoll
        spm11chi(local_index)=ecoll
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
        rm2c(local_index)=0.d0
        gam2c(local_index)=0.d0
        spm2clo(local_index)=ecoll
        spm2chi(local_index)=ecoll
      ELSE
        rm11c(local_index)=rmhh
        gam11c(local_index)=gamhh
        spm11clo(local_index)=spmhhlo
        spm11chi(local_index)=spmhhhi
        rm12c(local_index)=rmhh
        gam12c(local_index)=gamhh
        spm12clo(local_index)=spmhhlo
        spm12chi(local_index)=spmhhhi
        rm2c(local_index)=rmhh
        gam2c(local_index)=gamhh
        spm2clo(local_index)=spmhhlo
        spm2chi(local_index)=spmhhhi
      ENDIF

* qq -> {gq} {gqZ}                                ialfa(24) 2_4
      local_index=24
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm12a(local_index)=0.d0
        gam12a(local_index)=0.d0
        spm12alo(local_index)=0.d0
        spm12ahi(local_index)=0.d0
        rm12b(local_index)=rmz
        gam12b(local_index)=gamz
        spm12blo(local_index)=spmzlo
        spm12bhi(local_index)=spmzhi
      ELSE
        rm12a(local_index)=rmz
        gam12a(local_index)=gamz
        spm12alo(local_index)=spmzlo
        spm12ahi(local_index)=spmzhi
        rm12b(local_index)=rmh
        gam12b(local_index)=gamh
        spm12blo(local_index)=spmhlo
        spm12bhi(local_index)=spmhhi
      ENDIF
      IF(n5f.EQ.0 .OR. rmhh.GE.1.d3.OR.rmhh.LE.0.d0)THEN
        rm12c(local_index)=0.d0
        gam12c(local_index)=0.d0
        spm12clo(local_index)=ecoll
        spm12chi(local_index)=ecoll
      ELSE
        rm12c(local_index)=rmhh
        gam12c(local_index)=gamhh
        spm12clo(local_index)=spmhhlo
        spm12chi(local_index)=spmhhhi
      ENDIF

* t(WW) {ggZ}                                     ialfa(29) 1_1_31
      local_index=29
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm111a(local_index)=0.d0
        gam111a(local_index)=0.d0
        spm111alo(local_index)=0.d0
        spm111ahi(local_index)=0.d0
        rm111b(local_index)=rmz
        gam111b(local_index)=gamz
        spm111blo(local_index)=spmzlo
        spm111bhi(local_index)=spmzhi
      ELSE
        rm1a(local_index)=rmz
        gam1a(local_index)=gamz
        spm1alo(local_index)=spmzlo
        spm1ahi(local_index)=spmzhi
        rm1b(local_index)=rmh
        gam1b(local_index)=gamh
        spm1blo(local_index)=spmhlo
        spm1bhi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm111a(local_index)=rmz
        gam111a(local_index)=gamz
        spm111alo(local_index)=spmzlo
        spm111ahi(local_index)=spmzhi
        rm111b(local_index)=rmh
        gam111b(local_index)=gamh
        spm111blo(local_index)=spmhlo
        spm111bhi(local_index)=spmhhi
      ENDIF

* t(00) {ggZ}                                     ialfa(32) 1_1_31
      local_index=32
      IF(n5.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm111a(local_index)=0.d0
        gam111a(local_index)=0.d0
        spm111alo(local_index)=0.d0
        spm111ahi(local_index)=0.d0
        rm111b(local_index)=rmz
        gam111b(local_index)=gamz
        spm111blo(local_index)=spmzlo
        spm111bhi(local_index)=spmzhi
      ELSE
        rm1a(local_index)=rmz
        gam1a(local_index)=gamz
        spm1alo(local_index)=spmzlo
        spm1ahi(local_index)=spmzhi
        rm1b(local_index)=rmh
        gam1b(local_index)=gamh
        spm1blo(local_index)=spmhlo
        spm1bhi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm111a(local_index)=rmz
        gam111a(local_index)=gamz
        spm111alo(local_index)=spmzlo
        spm111ahi(local_index)=spmzhi
        rm111b(local_index)=rmh
        gam111b(local_index)=gamh
        spm111blo(local_index)=spmhlo
        spm111bhi(local_index)=spmhhi
      ENDIF

* t(W ) Z  t                                 ialfa(33) 1_2_3
      local_index=33
      IF(n5.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* g t(W ) Z {gW-}   and   g t(W ) Z {gW+}         ialfa(37) 1_2_3
      local_index=37
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* g t(W ) {gZ} W-   and   g t(W ) {gZ} W+         ialfa(38) 1_2_3
      local_index=38
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
      ELSE
        rm1a(local_index)=rmz
        gam1a(local_index)=gamz
        spm1alo(local_index)=spmzlo
        spm1ahi(local_index)=spmzhi
        rm1b(local_index)=rmh
        gam1b(local_index)=gamh
        spm1blo(local_index)=spmhlo
        spm1bhi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
      ENDIF

* g t(0 ) W+ {gW-}                                ialfa(39) 1_2_3
      local_index=39
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=0.d0
        gam0b(local_index)=0.d0
        spm0blo(local_index)=0.d0
        spm0bhi(local_index)=0.d0
      ELSE
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=rmh
        gam0b(local_index)=gamh
        spm0blo(local_index)=spmhlo
        spm0bhi(local_index)=spmhhi
      ENDIF

* g t(0 ) {gW+} W-                                ialfa(40) 1_2_3
      local_index=40
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=0.d0
        gam0b(local_index)=0.d0
        spm0blo(local_index)=0.d0
        spm0bhi(local_index)=0.d0
      ELSE
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=rmh
        gam0b(local_index)=gamh
        spm0blo(local_index)=spmhlo
        spm0bhi(local_index)=spmhhi
      ENDIF

* g t(0 ) {gZ} Z                                  ialfa(41) 1_2_3
      local_index=41
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=0.d0
        gam0b(local_index)=0.d0
        spm0blo(local_index)=0.d0
        spm0bhi(local_index)=0.d0
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSEIF(n5f.LT.2)THEN
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=rmh
        gam0b(local_index)=gamh
        spm0blo(local_index)=spmhlo
        spm0bhi(local_index)=spmhhi
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=rmh
        gam0b(local_index)=gamh
        spm0blo(local_index)=spmhlo
        spm0bhi(local_index)=spmhhi
        rm1a(local_index)=rmz
        gam1a(local_index)=gamz
        spm1alo(local_index)=spmzlo
        spm1ahi(local_index)=spmzhi
        rm1b(local_index)=rmh
        gam1b(local_index)=gamh
        spm1blo(local_index)=spmhlo
        spm1bhi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* g t(0 ) Z {gZ}                                  ialfa(42) 1_2_3
      local_index=42
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=0.d0
        gam0b(local_index)=0.d0
        spm0blo(local_index)=0.d0
        spm0bhi(local_index)=0.d0
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSEIF(n5f.LT.2)THEN
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=rmh
        gam0b(local_index)=gamh
        spm0blo(local_index)=spmhlo
        spm0bhi(local_index)=spmhhi
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm0a(local_index)=0.d0
        gam0a(local_index)=0.d0
        spm0alo(local_index)=0.d0
        spm0ahi(local_index)=0.d0
        rm0b(local_index)=rmh
        gam0b(local_index)=gamh
        spm0blo(local_index)=spmhlo
        spm0bhi(local_index)=spmhhi
        rm1a(local_index)=rmz
        gam1a(local_index)=gamz
        spm1alo(local_index)=spmzlo
        spm1ahi(local_index)=spmzhi
        rm1b(local_index)=rmh
        gam1b(local_index)=gamh
        spm1blo(local_index)=spmhlo
        spm1bhi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* t(W+ ) Z {ggq} -> t(W+ ) Z g{gq}                ialfa(47) 1_2_3
      local_index=47
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* t(W+ ) Z {ggq} -> t(W+ ) Z q{gg}                ialfa(48) 1_2_3
      local_index=48
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* t(W- ) Z {ggq} -> t(W- ) Z g{gq}                ialfa(49) 1_2_3
      local_index=49
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* t(W- ) Z {ggq} -> t(W- ) Z q{gg}                ialfa(50) 1_2_3
      local_index=50
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* t(0[1] ) Z {ggq} -> t(0[1] ) Z g{gq}            ialfa(51) 1_2_3
      local_index=51
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* t(0[1] ) Z {ggq} -> t(0[1] ) Z q{gg}            ialfa(52) 1_2_3
      local_index=52
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* t(0[2] ) Z {ggq} -> t(0[2] ) Z g{gq}            ialfa(53) 1_2_3
      local_index=53
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* t(0[2] ) Z {ggq} -> t(0[2] ) Z q{gg}            ialfa(54) 1_2_3
      local_index=54
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* W- -> t  Z  and  W+ -> t  Z                ialfa(55) 3_3
      local_index=55
      IF(n5f.NE.4 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm21a(local_index)=0.d0
        gam21a(local_index)=0.d0
        spm21alo(local_index)=0.d0
        spm21ahi(local_index)=0.d0
        rm21b(local_index)=rmz
        gam21b(local_index)=gamz
        spm21blo(local_index)=spmzlo
        spm21bhi(local_index)=spmzhi
      ELSE
        rm21a(local_index)=rmz
        gam21a(local_index)=gamz
        spm21alo(local_index)=spmzlo
        spm21ahi(local_index)=spmzhi
        rm21b(local_index)=rmh
        gam21b(local_index)=gamh
        spm21blo(local_index)=spmhlo
        spm21bhi(local_index)=spmhhi
      ENDIF

* W- -> {gZ} {gW-}  and  W+ -> {gZ} {gW+}         ialfa(57) 3_3
      local_index=57
      IF(num_z0flep(local_index).EQ.0)THEN
        IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
          rm1a(local_index)=0.d0
          gam1a(local_index)=0.d0
          spm1alo(local_index)=0.d0
          spm1ahi(local_index)=0.d0
          rm1b(local_index)=rmz
          gam1b(local_index)=gamz
          spm1blo(local_index)=spmzlo
          spm1bhi(local_index)=spmzhi
        ELSE
          rm1a(local_index)=rmz
          gam1a(local_index)=gamz
          spm1alo(local_index)=spmzlo
          spm1ahi(local_index)=spmzhi
          rm1b(local_index)=rmh
          gam1b(local_index)=gamh
          spm1blo(local_index)=spmhlo
          spm1bhi(local_index)=spmhhi
        ENDIF
      ELSE        ! Z leptonic, {gZ} cannot be resonant at Z/H mass
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=0.d0
        spm1bhi(local_index)=0.d0
      ENDIF
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
      ELSE
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
      ENDIF
      IF((num_wpflep(local_index)+num_wmflep(local_index)).EQ.0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmw
        gam2b(local_index)=gamw
        spm2blo(local_index)=spmwlo
        spm2bhi(local_index)=spmwhi
      ELSE        ! W leptonic, {gW} cannot be resonant at W mass
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=0.d0
        gam2b(local_index)=0.d0
        spm2blo(local_index)=0.d0
        spm2bhi(local_index)=0.d0
      ENDIF

* Z  -> {gW+} {gW-}                               ialfa(58) 3_3
      local_index=58
      IF(num_wpflep(local_index).EQ.0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmw
        gam1b(local_index)=gamw
        spm1blo(local_index)=spmwlo
        spm1bhi(local_index)=spmwhi
      ELSE            ! W+ leptonic, {gW+} cannot be resonant at W mass
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=0.d0
        gam1b(local_index)=0.d0
        spm1blo(local_index)=0.d0
        spm1bhi(local_index)=0.d0
      ENDIF
      IF(num_wmflep(local_index).EQ.0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmw
        gam2b(local_index)=gamw
        spm2blo(local_index)=spmwlo
        spm2bhi(local_index)=spmwhi
      ELSE            ! W- leptonic, {gW-} cannot be resonant at W mass
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=0.d0
        gam2b(local_index)=0.d0
        spm2blo(local_index)=0.d0
        spm2bhi(local_index)=0.d0
      ENDIF

* Z  -> {gZ} {gZ}                                 ialfa(59) 3_3
      local_index=59
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        IF(num_z0flep(local_index).EQ.0)THEN
          rm1b(local_index)=rmz
          gam1b(local_index)=gamz
          spm1blo(local_index)=spmzlo
          spm1bhi(local_index)=spmzhi
        ELSE          ! Z leptonic, {gZ} cannot be resonant at Z mass
          rm1b(local_index)=0.d0
          gam1b(local_index)=0.d0
          spm1blo(local_index)=0.d0
          spm1bhi(local_index)=0.d0
        ENDIF
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        IF(num_z0flep(local_index).EQ.0)THEN
          rm2b(local_index)=rmz
          gam2b(local_index)=gamz
          spm2blo(local_index)=spmzlo
          spm2bhi(local_index)=spmzhi
        ELSE          ! Z leptonic, {gZ} cannot be resonant at Z mass
          rm2b(local_index)=0.d0
          gam2b(local_index)=0.d0
          spm2blo(local_index)=0.d0
          spm2bhi(local_index)=0.d0
        ENDIF
        rm21a(local_index)=0.d0
        gam21a(local_index)=0.d0
        spm21alo(local_index)=0.d0
        spm21ahi(local_index)=0.d0
        rm21b(local_index)=rmz
        gam21b(local_index)=gamz
        spm21blo(local_index)=spmzlo
        spm21bhi(local_index)=spmzhi
      ELSE
c Notice: 'Z  -> {gZ} {gZ}' is called only if the two final Z's
c are both leptonic or non-leptonic (see proc.f). Here the two Z's in
c the final state are surely non-leptonic, hence the IF condition
c on num_z0flep(local_index) can be avoided.
        rm1a(local_index)=rmz
        gam1a(local_index)=gamz
        spm1alo(local_index)=spmzlo
        spm1ahi(local_index)=spmzhi
        rm1b(local_index)=rmh
        gam1b(local_index)=gamh
        spm1blo(local_index)=spmhlo
        spm1bhi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
        rm21a(local_index)=rmz
        gam21a(local_index)=gamz
        spm21alo(local_index)=spmzlo
        spm21ahi(local_index)=spmzhi
        rm21b(local_index)=rmh
        gam21b(local_index)=gamh
        spm21blo(local_index)=spmhlo
        spm21bhi(local_index)=spmhhi
      ENDIF

* W- -> {ggZ} W-     and  W+ -> {ggZ} W+          ialfa(61) 2_4to31
      local_index=61
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm111a(local_index)=0.d0
        gam111a(local_index)=0.d0
        spm111alo(local_index)=0.d0
        spm111ahi(local_index)=0.d0
        rm111b(local_index)=rmz
        gam111b(local_index)=gamz
        spm111blo(local_index)=spmzlo
        spm111bhi(local_index)=spmzhi
      ELSE
        rm1a(local_index)=rmz
        gam1a(local_index)=gamz
        spm1alo(local_index)=spmzlo
        spm1ahi(local_index)=spmzhi
        rm1b(local_index)=rmh
        gam1b(local_index)=gamh
        spm1blo(local_index)=spmhlo
        spm1bhi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm111a(local_index)=rmz
        gam111a(local_index)=gamz
        spm111alo(local_index)=spmzlo
        spm111ahi(local_index)=spmzhi
        rm111b(local_index)=rmh
        gam111b(local_index)=gamh
        spm111blo(local_index)=spmhlo
        spm111bhi(local_index)=spmhhi
      ENDIF

* W- -> Z {ggW-}    and  W+ -> Z {ggW+}           ialfa(62) 2_4to31
      local_index=62
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* Z  -> {ggZ} Z                                   ialfa(67) 2_4to31
      local_index=67
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm111a(local_index)=0.d0
        gam111a(local_index)=0.d0
        spm111alo(local_index)=0.d0
        spm111ahi(local_index)=0.d0
        rm111b(local_index)=rmz
        gam111b(local_index)=gamz
        spm111blo(local_index)=spmzlo
        spm111bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm1a(local_index)=rmz
        gam1a(local_index)=gamz
        spm1alo(local_index)=spmzlo
        spm1ahi(local_index)=spmzhi
        rm1b(local_index)=rmh
        gam1b(local_index)=gamh
        spm1blo(local_index)=spmhlo
        spm1bhi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm111a(local_index)=rmz
        gam111a(local_index)=gamz
        spm111alo(local_index)=spmzlo
        spm111ahi(local_index)=spmzhi
        rm111b(local_index)=rmh
        gam111b(local_index)=gamh
        spm111blo(local_index)=spmhlo
        spm111bhi(local_index)=spmhhi
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF

* Z  -> Z {ggZ}                                   ialfa(68) 2_4to31
      local_index=68
      IF(n5f.EQ.0 .OR. rmh.GE.1.d3.OR.rmh.LT.0.d0)THEN
        rm1a(local_index)=0.d0
        gam1a(local_index)=0.d0
        spm1alo(local_index)=0.d0
        spm1ahi(local_index)=0.d0
        rm1b(local_index)=rmz
        gam1b(local_index)=gamz
        spm1blo(local_index)=spmzlo
        spm1bhi(local_index)=spmzhi
        rm11a(local_index)=0.d0
        gam11a(local_index)=0.d0
        spm11alo(local_index)=0.d0
        spm11ahi(local_index)=0.d0
        rm11b(local_index)=rmz
        gam11b(local_index)=gamz
        spm11blo(local_index)=spmzlo
        spm11bhi(local_index)=spmzhi
        rm111a(local_index)=0.d0
        gam111a(local_index)=0.d0
        spm111alo(local_index)=0.d0
        spm111ahi(local_index)=0.d0
        rm111b(local_index)=rmz
        gam111b(local_index)=gamz
        spm111blo(local_index)=spmzlo
        spm111bhi(local_index)=spmzhi
        rm2a(local_index)=0.d0
        gam2a(local_index)=0.d0
        spm2alo(local_index)=0.d0
        spm2ahi(local_index)=0.d0
        rm2b(local_index)=rmz
        gam2b(local_index)=gamz
        spm2blo(local_index)=spmzlo
        spm2bhi(local_index)=spmzhi
      ELSE
        rm1a(local_index)=rmz
        gam1a(local_index)=gamz
        spm1alo(local_index)=spmzlo
        spm1ahi(local_index)=spmzhi
        rm1b(local_index)=rmh
        gam1b(local_index)=gamh
        spm1blo(local_index)=spmhlo
        spm1bhi(local_index)=spmhhi
        rm11a(local_index)=rmz
        gam11a(local_index)=gamz
        spm11alo(local_index)=spmzlo
        spm11ahi(local_index)=spmzhi
        rm11b(local_index)=rmh
        gam11b(local_index)=gamh
        spm11blo(local_index)=spmhlo
        spm11bhi(local_index)=spmhhi
        rm111a(local_index)=rmz
        gam111a(local_index)=gamz
        spm111alo(local_index)=spmzlo
        spm111ahi(local_index)=spmzhi
        rm111b(local_index)=rmh
        gam111b(local_index)=gamh
        spm111blo(local_index)=spmhlo
        spm111bhi(local_index)=spmhhi
        rm2a(local_index)=rmz
        gam2a(local_index)=gamz
        spm2alo(local_index)=spmzlo
        spm2ahi(local_index)=spmzhi
        rm2b(local_index)=rmh
        gam2b(local_index)=gamh
        spm2blo(local_index)=spmhlo
        spm2bhi(local_index)=spmhhi
      ENDIF


* end of phase space input


      RETURN
      END
