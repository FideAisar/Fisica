*
* procini.f
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


      SUBROUTINE procini

      IMPLICIT REAL*8 (a-b,d-h,o-z)
      implicit double complex (c)

C   Cteq convention:
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C
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
************************************************************************


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


      INCLUDE 'common.h'

      INCLUDE 'common_subproc.h'
      COMMON/phaones/ionesh
* the following variables represent the input of the various phase
* spaces:
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


* these variables are used in the phase space. Their value is specified
* once for all in the DATA, and it is not computed anymore
      DATA expa/mxnphs*0.8d0/, expb/mxnphs*0.8d0/

*** define phase space input
* Coefficient of the Higgs width in spmXlo and spmXhi: linear
* interpolation.
* Here we define the number of Higgs widths to be taken around the Higgs
* mass, for defining the phase space region where the mapping on the
* Higgs resonance must be switched on.
* variables {xhiggswidth, spmNlo, spmNhi where N=h,w,z,t} are defined
* just once, independently on the process.

* IF gamh > 650 GeV or rmh < 0 the Higgs resonance is not mapped
* When both a Z and a Higgs resonance are possible (n5 > 0) only the
* first is activated
      xxh=xhiggswidth(gamh)
      IF(xxh.LT.0.0)THEN
        spmhlo=0.d0
        spmhhi=0.d0
      ELSE
        spmhlo=rmh-xxh*gamh
        spmhhi=rmh+xxh*gamh
      ENDIF

      spmhlo=max(spmhlo,0.d0)
      spmhhi=min(spmhhi,ecoll)
* Heavy Higgs
      xxhh=xhiggswidth(gamhh)
      IF(xxh.LT.0.0)THEN
        spmhhlo=ecoll
        spmhhhi=ecoll
      ELSE
        spmhhlo=rmhh-xxhh*gamhh
        spmhhhi=rmhh+xxhh*gamhh
      ENDIF

      spmhhlo=max(spmhhlo,0.d0)
      spmhhhi=min(spmhhhi,ecoll)

      spmwlo=rmw-5.d0*gamw
      spmwhi=rmw+5.d0*gamw
      spmzlo=rmz-5.d0*gamz
      spmzhi=rmz+5.d0*gamz
      spmtlo=rmt-5.d0*gamt
      spmthi=rmt+5.d0*gamt

      IF(
     &    (spmzhi.GE.spmhlo .AND. rmh.GE.0.d0 .AND. spmhlo.GT.0.d0)
     &        .OR. 
     &    (spmhhi.GE.spmhhlo .AND. rmh.GE.0.d0 .AND. rmhh.GE.0.d0)
     &        .OR. 
     &    (rmh.GE.rmhh .AND. rmh.GE.0.d0 .AND. rmhh.GE.0.d0)
     &   )THEN
       WRITE(*,*)
     & '***************************************************************'
      WRITE(*,*)
     & '                          WARNING                              '
      WRITE(*,*)
     & ' PHANTOM requires rmz<rmh<rmhh and maps the Breit-Wigner peaks '
      WRITE(*,*)
     & ' in intervals whose width is computed in procini.f.            '
      WRITE(*,*)
     & ' The values given in input lead to overlapping intervals.      '
      WRITE(*,*)
     & ' This would produce wrong results therefore PHANTOM stops.     '
      WRITE(*,*)
     & '***************************************************************'
      stop
      ENDIF

* Continue defining phase space inputs.
* Some phase spaces need iproc to know how many b's are
* in the process. The number of b's determines whether or not the
* mapping on the Higgs resonance must be switched on in the phase space.
* Parameters which depend upon the Higgs resonance are declared in
* PROCEXTRAINI subroutine (see procextraini.f).



* t(00) W+ W-                                     ialfa(1)  1_1_4
      local_index=1
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b] declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rmta(local_index)=0.d0
      rmtb(local_index)=0.d0


* t(0W) Z W-  and t(0W) Z W+                      ialfa(2)  1_1_4
      local_index=2
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
c rm11[a,b,c] declared inside PROCEXTRAINI
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rmta(local_index)=0.d0
      rmtb(local_index)=rmw


* t(W0) Z W-  and t(W0) Z W+                      ialfa(3)  1_1_4
      local_index=3
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
c rm11[a,b,c] declared inside PROCEXTRAINI
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rmta(local_index)=rmw
      rmtb(local_index)=0.d0


* t(WW) W- W- and t(WW) W+ W+                     ialfa(4)  1_1_4
      local_index=4
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rmta(local_index)=rmw
      rmtb(local_index)=rmw


* t(WW) W+ W-                                     ialfa(5)  1_1_4
      local_index=5
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b] declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rmta(local_index)=rmw
      rmtb(local_index)=rmw


* t(00) Z Z                                       ialfa(6)  1_1_4
      local_index=6
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b], rm11[a,b,c], rm12[a,b,c] declared inside PROCEXTRAINI
      rmta(local_index)=0.d0
      rmtb(local_index)=0.d0


* t(WW) Z Z                                       ialfa(7)  1_1_4
      local_index=7
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b], rm11[a,b,c], rm12[a,b,c] declared inside PROCEXTRAINI
      rmta(local_index)=rmw
      rmtb(local_index)=rmw


* t(g0) W+ W-                                     ialfa(8)  1_1_4
      local_index=8
c rm1[a,b] declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rmta(local_index)=0.d0
      rmtb(local_index)=0.d0


* t(0g) W+ W-                                     ialfa(9)  1_1_4
      local_index=9
c rm1[a,b] declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rmta(local_index)=0.d0
      rmtb(local_index)=0.d0


* t(gW) Z W-   and  t(gW) Z W+                    ialfa(10) 1_1_4
      local_index=10
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
c rm11[a,b,c] declared inside PROCEXTRAINI
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rmta(local_index)=0.d0
      rmtb(local_index)=rmw


* t(Wg) Z W-   and  t(Wg) Z W+                    ialfa(11) 1_1_4
      local_index=11
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
c rm11[a,b,c] declared inside PROCEXTRAINI
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rmta(local_index)=rmw
      rmtb(local_index)=0.d0


* t(g0) Z Z                                       ialfa(12) 1_1_4
      local_index=12
c rm1[a,b], rm11[a,b,c], rm12[a,b,c] declared inside PROCEXTRAINI
      rmta(local_index)=0.d0
      rmtb(local_index)=0.d0


* t(0g) Z Z                                       ialfa(13) 1_1_4
      local_index=13
c rm1[a,b], rm11[a,b,c], rm12[a,b,c] declared inside PROCEXTRAINI
      rmta(local_index)=0.d0
      rmtb(local_index)=0.d0


* W- -> [W- W+] W-                                ialfa(14) 2_4
      local_index=14
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1 declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rm2c(local_index)=0.d0
      gam2c(local_index)=0.d0
      spm2clo(local_index)=ecoll
      spm2chi(local_index)=ecoll


* W- -> W- [W+ W-]                                ialfa(15) 2_4
      local_index=15
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b] declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rm2c(local_index)=0.d0
      gam2c(local_index)=0.d0
      spm2clo(local_index)=ecoll
      spm2chi(local_index)=ecoll


* W+ -> [W+ W-] W+                                ialfa(16) 2_4
      local_index=16
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b] declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rm2c(local_index)=0.d0
      gam2c(local_index)=0.d0
      spm2clo(local_index)=ecoll
      spm2chi(local_index)=ecoll


* W+ -> W+ [W- W+]                                ialfa(17) 2_4
      local_index=17
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1 declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rm2c(local_index)=0.d0
      gam2c(local_index)=0.d0
      spm2clo(local_index)=ecoll
      spm2chi(local_index)=ecoll


* Z  -> Z  W+ W-   and  gg -> Z  W+ W-            ialfa(18) 2_4
      local_index=18
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b] declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
c rm2[a,b,c] declared inside PROCEXTRAINI


* W- -> Z  Z  W- and W+ -> Z  Z  W+               ialfa(19) 2_4
      local_index=19
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b], rm11[a,b,c], rm12[a,b,c] declared inside PROCEXTRAINI
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rm2c(local_index)=0.d0
      gam2c(local_index)=0.d0
      spm2clo(local_index)=ecoll
      spm2chi(local_index)=ecoll


* Z  -> [Z Z] Z   and gg -> [Z Z] Z               ialfa(20) 2_4
      local_index=20
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b], rm11[a,b,c], rm12[a,b,c], rm2[a,b,c] declared inside PROCEXTRAINI


* Z  -> [Z Z Z]                                   ialfa(21) 2_4
      local_index=21
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b], rm11[a,b,c], rm12[a,b,c], rm2[a,b,c] declared inside PROCEXTRAINI


* Z  -> Z [Z Z]                                   ialfa(22) 2_4
      local_index=22
c giuseppe 30/10/2007
cx      rm1(local_index)=rmh
cx      gam1(local_index)=gamh
cx      spm1lo(local_index)=spmhlo
cx      spm1hi(local_index)=spmhhi
c end giuseppe 30/10/2007
c rm1[a,b], rm11[a,b,c], rm12[a,b,c], rm2[a,b,c] declared inside PROCEXTRAINI


* qq -> {gq} {gqW-}   and  qq -> {gq} {gqW+}      ialfa(23) 2_4
      local_index=23
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
      rm12a(local_index)=0.d0
      gam12a(local_index)=0.d0
      spm12alo(local_index)=0.d0
      spm12ahi(local_index)=0.d0
      rm12b(local_index)=rmw
      gam12b(local_index)=gamw
      spm12blo(local_index)=spmwlo
      spm12bhi(local_index)=spmwhi
      rm12c(local_index)=0.d0
      gam12c(local_index)=0.d0
      spm12clo(local_index)=ecoll
      spm12chi(local_index)=ecoll
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=0.d0
      gam2b(local_index)=0.d0
      spm2blo(local_index)=0.d0
      spm2bhi(local_index)=0.d0
      rm2c(local_index)=0.d0
      gam2c(local_index)=0.d0
      spm2clo(local_index)=ecoll
      spm2chi(local_index)=ecoll


* qq -> {gq} {gqZ}                                ialfa(24) 2_4
      local_index=24
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rm11c(local_index)=0.d0
      gam11c(local_index)=0.d0
      spm11clo(local_index)=ecoll
      spm11chi(local_index)=ecoll
c rm12[a,b,c] declared inside PROCEXTRAINI
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=0.d0
      gam2b(local_index)=0.d0
      spm2blo(local_index)=0.d0
      spm2bhi(local_index)=0.d0
      rm2c(local_index)=0.d0
      gam2c(local_index)=0.d0
      spm2clo(local_index)=ecoll
      spm2chi(local_index)=ecoll


* t(W0) t                                         ialfa(25) 1_1_31
      local_index=25
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmt
      gam11b(local_index)=gamt
      spm11blo(local_index)=spmtlo
      spm11bhi(local_index)=spmthi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rmta(local_index)=rmw
      rmtb(local_index)=0.d0


* t(0W) t                                         ialfa(26) 1_1_31
      local_index=26
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmt
      gam11b(local_index)=gamt
      spm11blo(local_index)=spmtlo
      spm11bhi(local_index)=spmthi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rmta(local_index)=0.d0
      rmtb(local_index)=rmw


* t(0W) {ggW-}   and   t(0W) {ggW+}               ialfa(27) 1_1_31
      local_index=27
c This phase space channel is activated in case of HADRONIC W only,
c hence {gW} and {ggW} can be resonant at W mass
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmw
      gam1b(local_index)=gamw
      spm1blo(local_index)=spmwlo
      spm1bhi(local_index)=spmwhi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rmta(local_index)=0.d0
      rmtb(local_index)=rmw


* t(W0) {ggW-}   and   t(W0) {ggW+}               ialfa(28) 1_1_31
      local_index=28
c This phase space channel is activated in case of HADRONIC W only,
c hence {gW} and {ggW} can be resonant at W mass
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmw
      gam1b(local_index)=gamw
      spm1blo(local_index)=spmwlo
      spm1bhi(local_index)=spmwhi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rmta(local_index)=rmw
      rmtb(local_index)=0.d0


* t(WW) {ggZ}                                     ialfa(29) 1_1_31
      local_index=29
c rm1[a,b], rm11[a,b], rm111[a,b] declared inside PROCEXTRAINI
      rmta(local_index)=rmw
      rmtb(local_index)=rmw


* t(gW) t                                         ialfa(30) 1_1_31
      local_index=30
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmt
      gam11b(local_index)=gamt
      spm11blo(local_index)=spmtlo
      spm11bhi(local_index)=spmthi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rmta(local_index)=0.d0
      rmtb(local_index)=rmw


* t(Wg) t                                         ialfa(31) 1_1_31
      local_index=31
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmt
      gam11b(local_index)=gamt
      spm11blo(local_index)=spmtlo
      spm11bhi(local_index)=spmthi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rmta(local_index)=rmw
      rmtb(local_index)=0.d0


* t(00) {ggZ}                                     ialfa(32) 1_1_31
      local_index=32
c rm1[a,b], rm11[a,b], rm111[a,b] declared inside PROCEXTRAINI
      rmta(local_index)=0.d0
      rmtb(local_index)=0.d0


* t(W ) Z  t                                      ialfa(33) 1_2_3
      local_index=33
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
      rm1b(local_index)=rmt
      gam1b(local_index)=gamt
      spm1blo(local_index)=spmtlo
      spm1bhi(local_index)=spmthi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
c rm2[a,b] declared inside PROCEXTRAINI
      rmta(local_index)=rmw


* t(0 ) W- t                                      ialfa(34) 1_2_3
      local_index=34
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
      rm1b(local_index)=rmt
      gam1b(local_index)=gamt
      spm1blo(local_index)=spmtlo
      spm1bhi(local_index)=spmthi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=0.d0


* t(0 ) W+ t                                      ialfa(35) 1_2_3
      local_index=35
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
      rm1b(local_index)=rmt
      gam1b(local_index)=gamt
      spm1blo(local_index)=spmtlo
      spm1bhi(local_index)=spmthi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=0.d0


* t(g ) W- t    and  t(g ) W+ t                   ialfa(36) 1_2_3
      local_index=36
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
      rm1b(local_index)=rmt
      gam1b(local_index)=gamt
      spm1blo(local_index)=spmtlo
      spm1bhi(local_index)=spmthi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=0.d0


* g t(W ) Z {gW-}   and   g t(W ) Z {gW+}         ialfa(37) 1_2_3
      local_index=37
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
      rm1b(local_index)=rmw
      gam1b(local_index)=gamw
      spm1blo(local_index)=spmwlo
      spm1bhi(local_index)=spmwhi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
c rm2[a,b] declared inside PROCEXTRAINI
      rmta(local_index)=rmw


* g t(W ) {gZ} W-   and   g t(W ) {gZ} W+         ialfa(38) 1_2_3
      local_index=38
      rm0a(local_index)=0.d0
      gam0a(local_index)=0.d0
      spm0alo(local_index)=0.d0
      spm0ahi(local_index)=0.d0
      rm0b(local_index)=0.d0
      gam0b(local_index)=0.d0
      spm0blo(local_index)=0.d0
      spm0bhi(local_index)=0.d0
c rm1[a,b], rm11[a,b] declared inside PROCEXTRAINI
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=rmw


* g t(0 ) W+ {gW-}                                ialfa(39) 1_2_3
      local_index=39
c rm0[a,b] declared inside PROCEXTRAINI
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmw
      gam1b(local_index)=gamw
      spm1blo(local_index)=spmwlo
      spm1bhi(local_index)=spmwhi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=0.d0


* g t(0 ) {gW+} W-                                ialfa(40) 1_2_3
      local_index=40
c rm0[a,b] declared inside PROCEXTRAINI
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmw
      gam1b(local_index)=gamw
      spm1blo(local_index)=spmwlo
      spm1bhi(local_index)=spmwhi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=0.d0


* g t(0 ) {gZ} Z                                  ialfa(41) 1_2_3
      local_index=41
c rm0[a,b], rm1[a,b], rm11[a,b], rm2[a,b] declared inside PROCEXTRAINI
      rmta(local_index)=0.d0


* g t(0 ) Z {gZ}                                  ialfa(42) 1_2_3
      local_index=42
c rm0[a,b], rm1[a,b], rm11[a,b], rm2[a,b] declared inside PROCEXTRAINI
      rmta(local_index)=0.d0


* t(0 ) W- {ggq} -> t(0 ) W- g{gq}
*          and  t(0 ) W+ {ggq} -> t(0 ) W+ g{gq}  ialfa(43) 1_2_3
      local_index=43
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=0.d0


* t(0 ) W- {ggq} -> t(0 ) W- q{gg}
*          and  t(0 ) W+ {ggq} -> t(0 ) W+ q{gg}  ialfa(44) 1_2_3
      local_index=44
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=0.d0


* t(W ) W  {ggq} -> t(W ) W  g{gq}                ialfa(45) 1_2_3
      local_index=45
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=rmw


* t(W ) W  {ggq} -> t(W ) W  q{gg}                ialfa(46) 1_2_3
      local_index=46
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi
      rmta(local_index)=rmw


* t(W+ ) Z {ggq} -> t(W+ ) Z g{gq}                ialfa(47) 1_2_3
      local_index=47
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rmta(local_index)=rmw
c rm2[a,b] declared inside PROCEXTRAINI


* t(W+ ) Z {ggq} -> t(W+ ) Z q{gg}                ialfa(48) 1_2_3
      local_index=48
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rmta(local_index)=rmw
c rm2[a,b] declared inside PROCEXTRAINI


* t(W- ) Z {ggq} -> t(W- ) Z g{gq}                ialfa(49) 1_2_3
      local_index=49
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rmta(local_index)=rmw
c rm2[a,b] declared inside PROCEXTRAINI


* t(W- ) Z {ggq} -> t(W- ) Z q{gg}                ialfa(50) 1_2_3
      local_index=50
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rmta(local_index)=rmw
c rm2[a,b] declared inside PROCEXTRAINI


* t(0[1] ) Z {ggq} -> t(0[1] ) Z g{gq}            ialfa(51) 1_2_3
      local_index=51
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rmta(local_index)=0.d0
c rm2[a,b] declared inside PROCEXTRAINI


* t(0[1] ) Z {ggq} -> t(0[1] ) Z q{gg}            ialfa(52) 1_2_3
      local_index=52
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rmta(local_index)=0.d0
c rm2[a,b] declared inside PROCEXTRAINI


* t(0[2] ) Z {ggq} -> t(0[2] ) Z g{gq}            ialfa(53) 1_2_3
      local_index=53
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rmta(local_index)=0.d0
c rm2[a,b] declared inside PROCEXTRAINI


* t(0[2] ) Z {ggq} -> t(0[2] ) Z q{gg}            ialfa(54) 1_2_3
      local_index=54
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
      rm1b(local_index)=0.d0
      gam1b(local_index)=0.d0
      spm1blo(local_index)=0.d0
      spm1bhi(local_index)=0.d0
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=0.d0
      gam11b(local_index)=0.d0
      spm11blo(local_index)=0.d0
      spm11bhi(local_index)=0.d0
      rmta(local_index)=0.d0
c rm2[a,b] declared inside PROCEXTRAINI


* W- -> t  Z  and  W+ -> t  Z                     ialfa(55) 3_3
      local_index=55
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmt
      gam1b(local_index)=gamt
      spm1blo(local_index)=spmtlo
      spm1bhi(local_index)=spmthi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=0.d0
      gam2b(local_index)=0.d0
      spm2blo(local_index)=0.d0
      spm2bhi(local_index)=0.d0
c rm21[a,b] declared inside PROCEXTRAINI


* Z  -> t  t   and  gg  -> t  t                   ialfa(56) 3_3
      local_index=56
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmt
      gam1b(local_index)=gamt
      spm1blo(local_index)=spmtlo
      spm1bhi(local_index)=spmthi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmt
      gam2b(local_index)=gamt
      spm2blo(local_index)=spmtlo
      spm2bhi(local_index)=spmthi
      rm21a(local_index)=0.d0
      gam21a(local_index)=0.d0
      spm21alo(local_index)=0.d0
      spm21ahi(local_index)=0.d0
      rm21b(local_index)=rmw
      gam21b(local_index)=gamw
      spm21blo(local_index)=spmwlo
      spm21bhi(local_index)=spmwhi


* W- -> {gZ} {gW-}  and  W+ -> {gZ} {gW+}         ialfa(57) 3_3
      local_index=57
c rm1[a,b], rm11[a,b], rm2[a,b] declared inside PROCEXTRAINI
      rm21a(local_index)=0.d0
      gam21a(local_index)=0.d0
      spm21alo(local_index)=0.d0
      spm21ahi(local_index)=0.d0
      rm21b(local_index)=rmw
      gam21b(local_index)=gamw
      spm21blo(local_index)=spmwlo
      spm21bhi(local_index)=spmwhi


* Z  -> {gW+} {gW-}                               ialfa(58) 3_3
      local_index=58
c rm1[a,b], rm2[a,b] declared inside PROCEXTRAINI
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm21a(local_index)=0.d0
      gam21a(local_index)=0.d0
      spm21alo(local_index)=0.d0
      spm21ahi(local_index)=0.d0
      rm21b(local_index)=rmw
      gam21b(local_index)=gamw
      spm21blo(local_index)=spmwlo
      spm21bhi(local_index)=spmwhi


* Z  -> {gZ} {gZ}                                 ialfa(59) 3_3
      local_index=59
c rm1[a,b], rm11[a,b], rm2[a,b], rm21[a,b] declared inside PROCEXTRAINI


* gq -> {gW+} t    and  gq -> {gW-} t             ialfa(60) 3_3
      local_index=60
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmw
      gam1b(local_index)=gamw
      spm1blo(local_index)=spmwlo
      spm1bhi(local_index)=spmwhi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmt
      gam2b(local_index)=gamt
      spm2blo(local_index)=spmtlo
      spm2bhi(local_index)=spmthi
      rm21a(local_index)=0.d0
      gam21a(local_index)=0.d0
      spm21alo(local_index)=0.d0
      spm21ahi(local_index)=0.d0
      rm21b(local_index)=rmw
      gam21b(local_index)=gamw
      spm21blo(local_index)=spmwlo
      spm21bhi(local_index)=spmwhi


* W- -> {ggZ} W-     and  W+ -> {ggZ} W+          ialfa(61) 2_4to31
      local_index=61
c rm1[a,b], rm11[a,b], rm111[a,b] declared inside PROCEXTRAINI
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi


* W- -> Z {ggW-}    and  W+ -> Z {ggW+}           ialfa(62) 2_4to31
      local_index=62
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmw
      gam1b(local_index)=gamw
      spm1blo(local_index)=spmwlo
      spm1bhi(local_index)=spmwhi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
c rm2[a,b] declared inside PROCEXTRAINI


* Z  -> {ggW+} W-                                 ialfa(63) 2_4to31
      local_index=63
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmw
      gam1b(local_index)=gamw
      spm1blo(local_index)=spmwlo
      spm1bhi(local_index)=spmwhi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi


* Z  -> W+ {ggW-}                                 ialfa(64) 2_4to31
      local_index=64
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmw
      gam1b(local_index)=gamw
      spm1blo(local_index)=spmwlo
      spm1bhi(local_index)=spmwhi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi


* gq -> W+ {gt} -> W+ g W-
*            and   gq -> W- {gt} -> W- g W+       ialfa(65) 2_4to31
      local_index=65
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmt
      gam1b(local_index)=gamt
      spm1blo(local_index)=spmtlo
      spm1bhi(local_index)=spmthi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmt
      gam11b(local_index)=gamt
      spm11blo(local_index)=spmtlo
      spm11bhi(local_index)=spmthi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi


* gq -> W+ g t -> W+ {gW-}
*            and   gq -> W- g t -> W- {gW+}       ialfa(66) 2_4to31
      local_index=66
      rm1a(local_index)=0.d0
      gam1a(local_index)=0.d0
      spm1alo(local_index)=0.d0
      spm1ahi(local_index)=0.d0
      rm1b(local_index)=rmt
      gam1b(local_index)=gamt
      spm1blo(local_index)=spmtlo
      spm1bhi(local_index)=spmthi
      rm11a(local_index)=0.d0
      gam11a(local_index)=0.d0
      spm11alo(local_index)=0.d0
      spm11ahi(local_index)=0.d0
      rm11b(local_index)=rmw
      gam11b(local_index)=gamw
      spm11blo(local_index)=spmwlo
      spm11bhi(local_index)=spmwhi
      rm111a(local_index)=0.d0
      gam111a(local_index)=0.d0
      spm111alo(local_index)=0.d0
      spm111ahi(local_index)=0.d0
      rm111b(local_index)=rmw
      gam111b(local_index)=gamw
      spm111blo(local_index)=spmwlo
      spm111bhi(local_index)=spmwhi
      rm2a(local_index)=0.d0
      gam2a(local_index)=0.d0
      spm2alo(local_index)=0.d0
      spm2ahi(local_index)=0.d0
      rm2b(local_index)=rmw
      gam2b(local_index)=gamw
      spm2blo(local_index)=spmwlo
      spm2bhi(local_index)=spmwhi


* Z  -> {ggZ} Z                                   ialfa(67) 2_4to31
      local_index=67
c rm1[a,b], rm11[a,b], rm111[a,b], rm2[a,b] declared in PROCEXTRAINI


* Z  -> Z {ggZ}                                   ialfa(68) 2_4to31
      local_index=68
c rm1[a,b], rm11[a,b], rm111[a,b], rm2[a,b] declared in PROCEXTRAINI


* W- -> {ggtbar} -> gg W-
*            and  W+ -> {ggt} -> gg W+            ialfa(69) 1_5to1_4to31
      local_index=69
      rm1(local_index)=rmt
      gam1(local_index)=gamt
      spm1lo(local_index)=spmtlo
      spm1hi(local_index)=spmthi
      rm11(local_index)=rmt
      gam11(local_index)=gamt
      spm11lo(local_index)=spmtlo
      spm11hi(local_index)=spmthi
      rm111(local_index)=rmt
      gam111(local_index)=gamt
      spm111lo(local_index)=spmtlo
      spm111hi(local_index)=spmthi
      rm1111(local_index)=rmw
      gam1111(local_index)=gamw
      spm1111lo(local_index)=spmwlo
      spm1111hi(local_index)=spmwhi
      rmta(local_index)=0.d0


* W- -> g {gtbar} -> g {gW-}
*            and  W+ -> g {gt} -> g {gW+}         ialfa(70) 1_5to1_4to31
      local_index=70
      rm1(local_index)=rmt
      gam1(local_index)=gamt
      spm1lo(local_index)=spmtlo
      spm1hi(local_index)=spmthi
      rm11(local_index)=rmt
      gam11(local_index)=gamt
      spm11lo(local_index)=spmtlo
      spm11hi(local_index)=spmthi
      rm111(local_index)=rmw
      gam111(local_index)=gamw
      spm111lo(local_index)=spmwlo
      spm111hi(local_index)=spmwhi
      rm1111(local_index)=rmw
      gam1111(local_index)=gamw
      spm1111lo(local_index)=spmwlo
      spm1111hi(local_index)=spmwhi
      rmta(local_index)=0.d0


* W- -> gg tbar -> {ggW-}
*            and  W+ -> gg t -> {ggW+}            ialfa(71) 1_5to1_4to31
      local_index=71
      rm1(local_index)=rmt
      gam1(local_index)=gamt
      spm1lo(local_index)=spmtlo
      spm1hi(local_index)=spmthi
      rm11(local_index)=rmw
      gam11(local_index)=gamw
      spm11lo(local_index)=spmwlo
      spm11hi(local_index)=spmwhi
      rm111(local_index)=rmw
      gam111(local_index)=gamw
      spm111lo(local_index)=spmwlo
      spm111hi(local_index)=spmwhi
      rm1111(local_index)=rmw
      gam1111(local_index)=gamw
      spm1111lo(local_index)=spmwlo
      spm1111hi(local_index)=spmwhi
      rmta(local_index)=0.d0


* t(W ) {ggtbar} -> t(W ) gg W-
*            and  t(W ) {ggt} -> t(W ) gg W+      ialfa(72) 1_5to1_4to31
      local_index=72
      rm1(local_index)=rmt
      gam1(local_index)=gamt
      spm1lo(local_index)=spmtlo
      spm1hi(local_index)=spmthi
      rm11(local_index)=rmt
      gam11(local_index)=gamt
      spm11lo(local_index)=spmtlo
      spm11hi(local_index)=spmthi
      rm111(local_index)=rmt
      gam111(local_index)=gamt
      spm111lo(local_index)=spmtlo
      spm111hi(local_index)=spmthi
      rm1111(local_index)=rmw
      gam1111(local_index)=gamw
      spm1111lo(local_index)=spmwlo
      spm1111hi(local_index)=spmwhi
      rmta(local_index)=rmw


* t(W ) g {gtbar} -> t(W ) g {gW-}
*            and  t(W ) g {gt} -> t(W ) g {gW+}   ialfa(73) 1_5to1_4to31
      local_index=73
      rm1(local_index)=rmt
      gam1(local_index)=gamt
      spm1lo(local_index)=spmtlo
      spm1hi(local_index)=spmthi
      rm11(local_index)=rmt
      gam11(local_index)=gamt
      spm11lo(local_index)=spmtlo
      spm11hi(local_index)=spmthi
      rm111(local_index)=rmw
      gam111(local_index)=gamw
      spm111lo(local_index)=spmwlo
      spm111hi(local_index)=spmwhi
      rm1111(local_index)=rmw
      gam1111(local_index)=gamw
      spm1111lo(local_index)=spmwlo
      spm1111hi(local_index)=spmwhi
      rmta(local_index)=rmw


* t(W ) gg tbar -> t(W ) {ggW-}
*            and  t(W ) gg t -> t(W ) {ggW+}      ialfa(74) 1_5to1_4to31
      local_index=74
      rm1(local_index)=rmt
      gam1(local_index)=gamt
      spm1lo(local_index)=spmtlo
      spm1hi(local_index)=spmthi
      rm11(local_index)=rmw
      gam11(local_index)=gamw
      spm11lo(local_index)=spmwlo
      spm11hi(local_index)=spmwhi
      rm111(local_index)=rmw
      gam111(local_index)=gamw
      spm111lo(local_index)=spmwlo
      spm111hi(local_index)=spmwhi
      rm1111(local_index)=rmw
      gam1111(local_index)=gamw
      spm1111lo(local_index)=spmwlo
      spm1111hi(local_index)=spmwhi
      rmta(local_index)=rmw


* g t(W ) {gtbar} -> g t(W ) g W-
*            and   g t(W ) {gt} -> g t(W ) g W+   ialfa(75) 1_5to1_4to31
      local_index=75
      rm1(local_index)=0.d0
      gam1(local_index)=0.d0
      spm1lo(local_index)=0.d0
      spm1hi(local_index)=0.d0
      rm11(local_index)=rmt
      gam11(local_index)=gamt
      spm11lo(local_index)=spmtlo
      spm11hi(local_index)=spmthi
      rm111(local_index)=rmt
      gam111(local_index)=gamt
      spm111lo(local_index)=spmtlo
      spm111hi(local_index)=spmthi
      rm1111(local_index)=rmw
      gam1111(local_index)=gamw
      spm1111lo(local_index)=spmwlo
      spm1111hi(local_index)=spmwhi
      rmta(local_index)=rmw


* g t(W ) g tbar -> g t(W ) {gW-}
*            and   g t(W ) g t -> g t(W ) {gW+}   ialfa(76) 1_5to1_4to31
      local_index=76
      rm1(local_index)=0.d0
      gam1(local_index)=0.d0
      spm1lo(local_index)=0.d0
      spm1hi(local_index)=0.d0
      rm11(local_index)=rmt
      gam11(local_index)=gamt
      spm11lo(local_index)=spmtlo
      spm11hi(local_index)=spmthi
      rm111(local_index)=rmw
      gam111(local_index)=gamw
      spm111lo(local_index)=spmwlo
      spm111hi(local_index)=spmwhi
      rm1111(local_index)=rmw
      gam1111(local_index)=gamw
      spm1111lo(local_index)=spmwlo
      spm1111hi(local_index)=spmwhi
      rmta(local_index)=rmw


* end of phase space input


      RETURN
      END
      
      
      
      
      REAL*8 FUNCTION xhiggswidth(gamh)
      IMPLICIT REAL*8 (a-b,d-h,o-z)
c Gamh(120 GeV) = 3.5d-3, Gamh(150 GeV) = 1.7d-2, Gamh(250 GeV) = 4.d0, Gamh(1000 GeV) = 650
      IF(gamh.LE.0.0d0)THEN
        xhiggswidth=-1.d0
      ELSEIF(gamh.LT.1.7d-2)THEN
        xhiggswidth=30.d0+(gamh-3.5d-3)/(1.7d-2-3.5d-3)*(20.d0-30.d0)
      ELSEIF(gamh.LT.4.0d0)THEN
        xhiggswidth=20.d0+(gamh-1.7d-2)/(4.0d0-1.7d-2)*(5.d0-20.d0)
      ELSEIF(gamh.LT.650.0d0)THEN
        xhiggswidth=5.d0+(gamh-4.0d0)/(650.0d0-4.0d0)*(1.d0-5.d0)
      ELSE
        xhiggswidth=-1.d0
      END IF
      RETURN
      END

