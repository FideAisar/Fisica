************************************************************************
* common.h
*
* Last update: Jul 30, 2014
*
c Version extended to e+e- initial state.
c Added two common flags for BEAMSTRAHLUNG and ISR
c Added flag for printing intermediate particles (mothers) in event file
c Added flag for massless amplitudes
************************************************************************

c these commons are evaluated once in coupling.f and are never changed
*cmass

       COMPLEX*16    sw,s2w,rcw,rc2w,rcotw,rcot2w
     
       COMPLEX*16  fcr,fcl,zcr,zcl,wcl
     &     ,rhbb,rhtt,rhww,rhzz,r2h2w,r2h2z,r3h,r4h
     &     ,rhhbb,rhhtt,rhhww,rhhzz,r2hh2w,r2hh2z,r3hh,r4hh
     &     ,r1h1hh2w,r1h1hh2z,r1h2hh,r2h1hh,r1h3hh,r2h2hh,r3h1hh
 
*cmassend
*cmass
c      common/phacou1/gf,sw,s2w,rcw,rc2w,rcotw,rcot2w,
c     &     alfainv,elcharge2,alfainv_me,alfa_s,gs2
c      common/phacou2/fcr(-21:21),fcl(-21:21),zcr(-21:21),zcl(-21:21),wcl
cc s14
cc     &     ,rhbb,rhtt    
c     &     ,rhbb,rhtt,rhww,rhzz,r2h2w,r2h2z,r3h,r4h
c* heavh
c     &     ,rhhbb,rhhtt,rhhww,rhhzz,r2hh2w,r2hh2z,r3hh,r4hh
c     &     ,r1h1hh2w,r1h1hh2z,r1h2hh,r2h1hh,r1h3hh,r2h2hh,r3h1hh
c* heavhend
cc s14end    

      common/phacoulreal/alfainv,elcharge2,alfa_s,gs2
      common/phacou1/sw,s2w,rcw,rc2w,rcotw,rcot2w
      common/phacou2/fcr(-21:21),fcl(-21:21),zcr(-21:21),zcl(-21:21),wcl
     &     ,rhbb,rhtt,rhww,rhzz,r2h2w,r2h2z,r3h,r4h
     &     ,rhhbb,rhhtt,rhhww,rhhzz,r2hh2w,r2hh2z,r3hh,r4hh
     &     ,r1h1hh2w,r1h1hh2z,r1h2hh,r2h1hh,r1h3hh,r2h2hh,r3h1hh
*cmassend
      common/phamawi/rmw,rmw2,gamw,rmz,rmz2,gamz,rmt,rmt2,gamt,rmb,rmb2,
* heavh
*     &     rmh,rmh2,gamh,rmass(-21:21),rmass2(-21:21)
     &     rmh,rmh2,gamh,rmhh,rmhh2,gamhh,rmass(-21:21),rmass2(-21:21)
* heavhend
      common/phapar/pi,twopi,fourpi
* heavh
*      COMPLEX*16 czero,cuno,cim,cmw2,cmz2,cmh2,cmt2,cmass2
*      common/phacomp/czero,cuno,cim,cmw2,cmz2,cmh2,cmt2,cmass2(-21:21)
      COMPLEX*16 czero,cuno,cim,cmw2,cmz2,cmh2,cmhh2,cmt2,cmass2
      common/phacomp/czero,cuno,cim,cmw2,cmz2,cmh2,
     &                               cmhh2,cmt2,cmass2(-21:21)
* heavhend
      common/phaint/ilept(-21:21),ineutri(-21:21),iup(-21:21)
     &     ,imass(-21:21)


c collider energy common
      COMMON/collenergy/ecoll,s

c these commons are necessary for flat event generation. 
c  /phaiflav/ is the one present in phavegas
 
      COMMON/phaifla/iflat,istorvegas,iwrite_event,iwrite_mothers,
     &     ihadronize,novermax
      COMMON/phaiflav/nevent,nflevts
*6
c      COMMON/pharfla/scalemax,rmaxfxn,rmaxfxn_2
      COMMON/pharfla/scalemax,rmaxfxn,rmaxfxnnew
*6end
c this common decides how many events are generated in oneshot
      COMMON/phaonev/nunwevts

c this common gives integration limits defined in main to integration.f
c   and oneshot.f

      COMMON/phareg/region(30),ndim  

c this common determines whether one wants to sum over cc and fam 
c  symmetries, sum over the exchange of the 
c  incoming particles if not identical, change at random the couple
c  lepton neutrino.
c  i_ccfam must be the same in single run and oneshot generation
c  i_exchincoming and i_emutau are valid only for oneshot. 
c  In single run one does not sum over exhange of incoming and change
c  the leptonic couple.
  
*sandro2/3/07
cc  the variable i_nofullew  eliminates from 8 fermion qcd amplitudes the 
cc    order alpha_em^6 when set =1 , it computes all when set =0 
c
c      COMMON/icontrol/i_ccfam,i_exchincoming,i_emutau,i_nofullew  

c i_pertorder = 1 alpha_em^6 with dedicated amp
c              = 2 alpha_s^2alpha_em^4
c              = 3 alpha_em^6 + alpha_s^2alpha_em^4
c              = 0 alpha_em^6 with amp8fqcd (for test only)

c diogo 29/03/2009
c   i_massive = 0 use massless amp unless there is at least a b quark
c             = 1 always use massive amplitudes (massive Z-lines)

c      COMMON/icontrol/i_ccfam,i_exchincoming,i_emutau,i_pertorder
      COMMON/icontrol/i_ccfam,i_exchincoming,i_emutau,i_pertorder,
* sig
c     &     i_massive
     &     i_massive,i_signal
*zw
**res
*     &     ,nleptfromw1,nleptfromw2
**resend
     &     ,i_ww,idw(4),i_zz,idz(4)
     &     ,ipolw(2),ipolz(2),i_osp,idosp(4)
*4cmpol
     &     ,i_4cmpol,id4cmpol(4)  
*4cmpolend

*zwend
* sigend

c diogoend


*sandro2/3/07end


c giuseppe ILC 25/06/2007
c this common determines the type of accelerator
c icoll=1 => p-p  i_coll=2 => p-pbar  i_coll=3 => e+e-
 
      COMMON/collider/i_coll

c this common contains the input flags for the inclusion of 
c BEAMSTRAHLUNG and ISR in case of e+e- initial state

      COMMON/extrarad/i_isr,i_beamstrahlung
c end giuseppe ILC 25/06/2007
