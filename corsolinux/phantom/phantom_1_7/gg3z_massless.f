      subroutine ggzzz_massless(p1,p2,p3,p4,p5,p6,p7,p8,
     &              id1,id2,id3,id4,id5,id6,id7,id8,cres)
  
      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)
  
      include 'phact_data_types_massless.inc'
      include 'common.h'
  
      dimension cres(2,2,2,2,2,12)
  
* VARIABLE DEFINITIONS                                                  
*four momenta                                                           
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
     &         ,p5(0:3),p6(0:3),p7(0:3),p8(0:3)
* auxiliaries                                                           
      type(tl0) laux_imu(2,0:3),laux_iii(2,2,2)
*      type(t0) uaux_mu(0:3)                                            
  
* -> for forks                                                          
*    -simple fork                                                       
*      dimension p34(0;3)                                               
      type(polcom) cz34(2),cf34(2)
*      dimension p56(0;3)                                               
      type(polcom) cz56(2),cf56(2)
*      dimension p78(0;3)                                               
      type(polcom) cz78(2),cf78(2)
  
      dimension p12(0:3)
      type(polcom) ce1(2),ce2(2),ctrip12(2,2)
  
*    -one-gluon fork                                                    
      dimension p31(0:3),p41(0:3)
      type(polcom) cz314(2,2),cf314(2,2)
      type(tl0) l3_1(2),lz341(0:3),lf341(0:3)
      type(tr0) r4_1(2),rz431(0:3),rf431(0:3)
      dimension p32(0:3),p42(0:3)
      type(polcom) cz324(2,2),cf324(2,2)
      type(tl0) l3_2(2),lz342(0:3),lf342(0:3)
      type(tr0) r4_2(2),rz432(0:3),rf432(0:3)
      dimension p51(0:3),p61(0:3)
      type(polcom) cz516(2,2),cf516(2,2)
      type(tl0) l5_1(2),lz561(0:3),lf561(0:3)
      type(tr0) r6_1(2),rz651(0:3),rf651(0:3)
      dimension p52(0:3),p62(0:3)
      type(polcom) cz526(2,2),cf526(2,2)
      type(tl0) l5_2(2),lz562(0:3),lf562(0:3)
      type(tr0) r6_2(2),rz652(0:3),rf652(0:3)
      dimension p71(0:3),p81(0:3)
      type(polcom) cz718(2,2),cf718(2,2)
      type(tl0) l7_1(2),lz781(0:3),lf781(0:3)
      type(tr0) r8_1(2),rz871(0:3),rf871(0:3)
      dimension p72(0:3),p82(0:3)
      type(polcom) cz728(2,2),cf728(2,2)
      type(tl0) l7_2(2),lz782(0:3),lf782(0:3)
      type(tr0) r8_2(2),rz872(0:3),rf872(0:3)
  
*     -two-gluons fork                                                  
      dimension p312(0:3),p412(0:3)
      type(tl0) l3_gg(2,2),lz3412(0:3),lf3412(0:3)
      type(tr0) r4_gg(2,2),rz4312(0:3),rf4312(0:3)
      type(polcom) cz3124(2,2,2),cf3124(2,2,2)
      type(t0) u31_2(2),u41_2(2),uz31(0:3),uf31(0:3)
      type(tl0) l3_12(2,2)
      type(tr0) r4_12(2,2)
      type(polcom) cz3214(2,2,2),cf3214(2,2,2)
      type(t0) u32_1(2),u42_1(2),uz32(0:3),uf32(0:3)
      type(tl0) l3_21(2,2)
      type(tr0) r4_21(2,2)
      dimension p512(0:3),p612(0:3)
      type(tl0) l5_gg(2,2),lz5612(0:3),lf5612(0:3)
      type(tr0) r6_gg(2,2),rz6512(0:3),rf6512(0:3)
      type(polcom) cz5126(2,2,2),cf5126(2,2,2)
      type(t0) u51_2(2),u61_2(2),uz51(0:3),uf51(0:3)
      type(tl0) l5_12(2,2)
      type(tr0) r6_12(2,2)
      type(polcom) cz5216(2,2,2),cf5216(2,2,2)
      type(t0) u52_1(2),u62_1(2),uz52(0:3),uf52(0:3)
      type(tl0) l5_21(2,2)
      type(tr0) r6_21(2,2)
      dimension p712(0:3),p812(0:3)
      type(tl0) l7_gg(2,2),lz7812(0:3),lf7812(0:3)
      type(tr0) r8_gg(2,2),rz8712(0:3),rf8712(0:3)
      type(polcom) cz7128(2,2,2),cf7128(2,2,2)
      type(t0) u71_2(2),u81_2(2),uz71(0:3),uf71(0:3)
      type(tl0) l7_12(2,2)
      type(tr0) r8_12(2,2)
      type(polcom) cz7218(2,2,2),cf7218(2,2,2)
      type(t0) u72_1(2),u82_1(2),uz72(0:3),uf72(0:3)
      type(tl0) l7_21(2,2)
      type(tr0) r8_21(2,2)
  
* -> lines without gluon                                                
      dimension p356(0:3),p456(0:3)
      type(tl0) l3_56(2)
      type(tr0) r4_56(2)
  
      dimension p3156(0:3)
      type(tl0) l3_516(2,2),l3_5126(2,2,2)
      type(tr0) r4_516(2,2),r4_5126(2,2,2)
  
      dimension p3256(0:3)
      type(tl0) l3_526(2,2),l3_5216(2,2,2)
      type(tr0) r4_526(2,2),r4_5216(2,2,2)
      dimension p378(0:3),p478(0:3)
      type(tl0) l3_78(2)
      type(tr0) r4_78(2)
  
      dimension p3178(0:3)
      type(tl0) l3_718(2,2),l3_7128(2,2,2)
      type(tr0) r4_718(2,2),r4_7128(2,2,2)
  
      dimension p3278(0:3)
      type(tl0) l3_728(2,2),l3_7218(2,2,2)
      type(tr0) r4_728(2,2),r4_7218(2,2,2)
      dimension p534(0:3),p634(0:3)
      type(tl0) l5_34(2)
      type(tr0) r6_34(2)
  
      dimension p5134(0:3)
      type(tl0) l5_314(2,2),l5_3124(2,2,2)
      type(tr0) r6_314(2,2),r6_3124(2,2,2)
  
      dimension p5234(0:3)
      type(tl0) l5_324(2,2),l5_3214(2,2,2)
      type(tr0) r6_324(2,2),r6_3214(2,2,2)
      dimension p578(0:3),p678(0:3)
      type(tl0) l5_78(2)
      type(tr0) r6_78(2)
  
      dimension p5178(0:3)
      type(tl0) l5_718(2,2),l5_7128(2,2,2)
      type(tr0) r6_718(2,2),r6_7128(2,2,2)
  
      dimension p5278(0:3)
      type(tl0) l5_728(2,2),l5_7218(2,2,2)
      type(tr0) r6_728(2,2),r6_7218(2,2,2)
      dimension p734(0:3),p834(0:3)
      type(tl0) l7_34(2)
      type(tr0) r8_34(2)
  
      dimension p7134(0:3)
      type(tl0) l7_314(2,2),l7_3124(2,2,2)
      type(tr0) r8_314(2,2),r8_3124(2,2,2)
  
      dimension p7234(0:3)
      type(tl0) l7_324(2,2),l7_3214(2,2,2)
      type(tr0) r8_324(2,2),r8_3214(2,2,2)
      dimension p756(0:3),p856(0:3)
      type(tl0) l7_56(2)
      type(tr0) r8_56(2)
  
      dimension p7156(0:3)
      type(tl0) l7_516(2,2),l7_5126(2,2,2)
      type(tr0) r8_516(2,2),r8_5126(2,2,2)
  
      dimension p7256(0:3)
      type(tl0) l7_526(2,2),l7_5216(2,2,2)
      type(tr0) r8_526(2,2),r8_5216(2,2,2)
  
* -> lines with one gluon                                               
      type(tl0) l3_156(2,2),l3_1526(2,2,2),
     &  l3_56718(2,2,2)
      type(t0) u31_56(2),u31_526(2,2),u356_1(2),
     &  u3516_2(2),u356_718(2,2),u3516_78(2)
  
      type(tl0) l3_256(2,2),l3_2516(2,2,2),
     &  l3_56728(2,2,2)
      type(t0) u32_56(2),u32_516(2,2),u356_2(2),
     &  u3526_1(2),u356_728(2,2),u3526_78(2)
  
      type(tl0) l3_178(2,2),l3_1728(2,2,2),
     &  l3_78516(2,2,2)
      type(t0) u31_78(2),u31_728(2,2),u378_1(2),
     &  u3718_2(2),u378_516(2,2),u3718_56(2)
  
      type(tl0) l3_278(2,2),l3_2718(2,2,2),
     &  l3_78526(2,2,2)
      type(t0) u32_78(2),u32_718(2,2),u378_2(2),
     &  u3728_1(2),u378_526(2,2),u3728_56(2)
  
      type(tl0) l5_134(2,2),l5_1324(2,2,2),
     &  l5_34718(2,2,2)
      type(t0) u51_34(2),u51_324(2,2),u534_1(2),
     &  u5314_2(2),u534_718(2,2),u5314_78(2)
  
      type(tl0) l5_234(2,2),l5_2314(2,2,2),
     &  l5_34728(2,2,2)
      type(t0) u52_34(2),u52_314(2,2),u534_2(2),
     &  u5324_1(2),u534_728(2,2),u5324_78(2)
  
      type(tl0) l5_178(2,2),l5_1728(2,2,2),
     &  l5_78314(2,2,2)
      type(t0) u51_78(2),u51_728(2,2),u578_1(2),
     &  u5718_2(2),u578_314(2,2),u5718_34(2)
  
      type(tl0) l5_278(2,2),l5_2718(2,2,2),
     &  l5_78324(2,2,2)
      type(t0) u52_78(2),u52_718(2,2),u578_2(2),
     &  u5728_1(2),u578_324(2,2),u5728_34(2)
  
      type(tl0) l7_134(2,2),l7_1324(2,2,2),
     &  l7_34516(2,2,2)
      type(t0) u71_34(2),u71_324(2,2),u734_1(2),
     &  u7314_2(2),u734_516(2,2),u7314_56(2)
  
      type(tl0) l7_234(2,2),l7_2314(2,2,2),
     &  l7_34526(2,2,2)
      type(t0) u72_34(2),u72_314(2,2),u734_2(2),
     &  u7324_1(2),u734_526(2,2),u7324_56(2)
  
      type(tl0) l7_156(2,2),l7_1526(2,2,2),
     &  l7_56314(2,2,2)
      type(t0) u71_56(2),u71_526(2,2),u756_1(2),
     &  u7516_2(2),u756_314(2,2),u7516_34(2)
  
      type(tl0) l7_256(2,2),l7_2516(2,2,2),
     &  l7_56324(2,2,2)
      type(t0) u72_56(2),u72_516(2,2),u756_2(2),
     &  u7526_1(2),u756_324(2,2),u7526_34(2)
  
  
* -> lines with two gluons                                              
      type(t0) u312_56(2),u312_78(2),u356_78(2),u378_56(2),
     &  u356_gg(2,2),u378_gg(2,2)
      type(tr0) r4_5678(2,2),r4_156(2,2),r4_256(2,2),
     &  r4_178(2,2),r4_278(2,2)
      type(tl0) l3_5678(2,2)
      type(t0) u512_34(2),u512_78(2),u534_78(2),u578_34(2),
     &  u534_gg(2,2),u578_gg(2,2)
      type(tr0) r6_3478(2,2),r6_134(2,2),r6_234(2,2),
     &  r6_178(2,2),r6_278(2,2)
      type(tl0) l5_3478(2,2)
      type(t0) u712_34(2),u712_56(2),u734_56(2),u756_34(2),
     &  u734_gg(2,2),u756_gg(2,2)
      type(tr0) r8_3456(2,2),r8_134(2,2),r8_234(2,2),
     &  r8_156(2,2),r8_256(2,2)
      type(tl0) l7_3456(2,2)
  
* COMPUTATION --------------------------------------------------------  
*four momenta k0                                                        
* pk0 -- p=p3
      p3k0=p3(0)-p3(1)
* pk0 -- p=p4
      p4k0=p4(0)-p4(1)
* pk0 -- p=p5
      p5k0=p5(0)-p5(1)
* pk0 -- p=p6
      p6k0=p6(0)-p6(1)
* pk0 -- p=p7
      p7k0=p7(0)-p7(1)
* pk0 -- p=p8
      p8k0=p8(0)-p8(1)
  
  
* -> for forks                                                          
*    -simple fork                                                       
  
* quqd -- p=p3,q=p4
      quqd=p3(0)*p4(0)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3)
      s34=2.d0*quqd
  
      ccr=zcr(id3)/(-s34+cmz2)
      ccl=zcl(id3)/(-s34+cmz2)
* T10 -- qu=p3,qd=p4,v=0,a=cz34(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p4(3)+p4(2)*p3(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p3k0*p4(0)+p4k0*p3(0)
      cz34(1)%e(0)=ccr*(auxa+ceps_0)
      cz34(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=1,a=cz34(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p4(1)+p4k0*p3(1)
      cz34(1)%e(1)=ccr*(auxa+ceps_0)
      cz34(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=2,a=cz34(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p4(3)+p4k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(2)+p4k0*p3(2)
      cz34(1)%e(2)=ccr*(auxa+ceps_0)
      cz34(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=3,a=cz34(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p4(2)-p4k0*p3(2)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(3)+p4k0*p3(3)
      cz34(1)%e(3)=ccr*(auxa+ceps_0)
      cz34(2)%e(3)=ccl*(auxa-ceps_0)
  
      do i3=1,2
* pk0 -- p=cz34(i3)%e
      cz34(i3)%ek0=cz34(i3)%e(0)-cz34(i3)%e(1)
      end do
  
      if (ineutri(id3).ne.1) then
      ccr=fcr(id3)/(-s34)
      ccl=fcl(id3)/(-s34)
* T10 -- qu=p3,qd=p4,v=0,a=cf34(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p4(3)+p4(2)*p3(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p3k0*p4(0)+p4k0*p3(0)
      cf34(1)%e(0)=ccr*(auxa+ceps_0)
      cf34(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=1,a=cf34(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p4(1)+p4k0*p3(1)
      cf34(1)%e(1)=ccr*(auxa+ceps_0)
      cf34(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=2,a=cf34(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p4(3)+p4k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(2)+p4k0*p3(2)
      cf34(1)%e(2)=ccr*(auxa+ceps_0)
      cf34(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=3,a=cf34(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p4(2)-p4k0*p3(2)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(3)+p4k0*p3(3)
      cf34(1)%e(3)=ccr*(auxa+ceps_0)
      cf34(2)%e(3)=ccl*(auxa-ceps_0)
  
      do i3=1,2
* pk0 -- p=cf34(i3)%e
      cf34(i3)%ek0=cf34(i3)%e(0)-cf34(i3)%e(1)
      end do
      endif
  
  
* quqd -- p=p5,q=p6
      quqd=p5(0)*p6(0)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3)
      s56=2.d0*quqd
  
      ccr=zcr(id5)/(-s56+cmz2)
      ccl=zcl(id5)/(-s56+cmz2)
* T10 -- qu=p5,qd=p6,v=0,a=cz56(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p5(2)*p6(3)+p6(2)*p5(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p5k0*p6(0)+p6k0*p5(0)
      cz56(1)%e(0)=ccr*(auxa+ceps_0)
      cz56(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=1,a=cz56(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p6(1)+p6k0*p5(1)
      cz56(1)%e(1)=ccr*(auxa+ceps_0)
      cz56(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=2,a=cz56(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p5k0*p6(3)+p6k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(2)+p6k0*p5(2)
      cz56(1)%e(2)=ccr*(auxa+ceps_0)
      cz56(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=3,a=cz56(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p5k0*p6(2)-p6k0*p5(2)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(3)+p6k0*p5(3)
      cz56(1)%e(3)=ccr*(auxa+ceps_0)
      cz56(2)%e(3)=ccl*(auxa-ceps_0)
  
      do i5=1,2
* pk0 -- p=cz56(i5)%e
      cz56(i5)%ek0=cz56(i5)%e(0)-cz56(i5)%e(1)
      end do
  
      if (ineutri(id5).ne.1) then
      ccr=fcr(id5)/(-s56)
      ccl=fcl(id5)/(-s56)
* T10 -- qu=p5,qd=p6,v=0,a=cf56(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p5(2)*p6(3)+p6(2)*p5(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p5k0*p6(0)+p6k0*p5(0)
      cf56(1)%e(0)=ccr*(auxa+ceps_0)
      cf56(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=1,a=cf56(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p6(1)+p6k0*p5(1)
      cf56(1)%e(1)=ccr*(auxa+ceps_0)
      cf56(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=2,a=cf56(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p5k0*p6(3)+p6k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(2)+p6k0*p5(2)
      cf56(1)%e(2)=ccr*(auxa+ceps_0)
      cf56(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=3,a=cf56(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p5k0*p6(2)-p6k0*p5(2)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(3)+p6k0*p5(3)
      cf56(1)%e(3)=ccr*(auxa+ceps_0)
      cf56(2)%e(3)=ccl*(auxa-ceps_0)
  
      do i5=1,2
* pk0 -- p=cf56(i5)%e
      cf56(i5)%ek0=cf56(i5)%e(0)-cf56(i5)%e(1)
      end do
      endif
  
  
* quqd -- p=p7,q=p8
      quqd=p7(0)*p8(0)-p7(1)*p8(1)-p7(2)*p8(2)-p7(3)*p8(3)
      s78=2.d0*quqd
  
      ccr=zcr(id7)/(-s78+cmz2)
      ccl=zcl(id7)/(-s78+cmz2)
* T10 -- qu=p7,qd=p8,v=0,a=cz78(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p7(2)*p8(3)+p8(2)*p7(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p7k0*p8(0)+p8k0*p7(0)
      cz78(1)%e(0)=ccr*(auxa+ceps_0)
      cz78(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=1,a=cz78(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p8(1)+p8k0*p7(1)
      cz78(1)%e(1)=ccr*(auxa+ceps_0)
      cz78(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=2,a=cz78(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p7k0*p8(3)+p8k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(2)+p8k0*p7(2)
      cz78(1)%e(2)=ccr*(auxa+ceps_0)
      cz78(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=3,a=cz78(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p7k0*p8(2)-p8k0*p7(2)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(3)+p8k0*p7(3)
      cz78(1)%e(3)=ccr*(auxa+ceps_0)
      cz78(2)%e(3)=ccl*(auxa-ceps_0)
  
      do i7=1,2
* pk0 -- p=cz78(i7)%e
      cz78(i7)%ek0=cz78(i7)%e(0)-cz78(i7)%e(1)
      end do
  
      if (ineutri(id7).ne.1) then
      ccr=fcr(id7)/(-s78)
      ccl=fcl(id7)/(-s78)
* T10 -- qu=p7,qd=p8,v=0,a=cf78(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p7(2)*p8(3)+p8(2)*p7(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p7k0*p8(0)+p8k0*p7(0)
      cf78(1)%e(0)=ccr*(auxa+ceps_0)
      cf78(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=1,a=cf78(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p8(1)+p8k0*p7(1)
      cf78(1)%e(1)=ccr*(auxa+ceps_0)
      cf78(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=2,a=cf78(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p7k0*p8(3)+p8k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(2)+p8k0*p7(2)
      cf78(1)%e(2)=ccr*(auxa+ceps_0)
      cf78(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=3,a=cf78(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p7k0*p8(2)-p8k0*p7(2)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(3)+p8k0*p7(3)
      cf78(1)%e(3)=ccr*(auxa+ceps_0)
      cf78(2)%e(3)=ccl*(auxa-ceps_0)
  
      do i7=1,2
* pk0 -- p=cf78(i7)%e
      cf78(i7)%ek0=cf78(i7)%e(0)-cf78(i7)%e(1)
      end do
      endif
  
  
* gluon polarizations                                                   
  
* pol -- p=p1,m=0,e1=ce1(1).e,e2=ce1(2).e,e3=0,nhz=0                    
      rkt=sqrt(p1(1)**2+p1(2)**2)
      IF(abs(rkt/p1(0)).GT.1.d-10)THEN
      ce1(1)%e(1)=p1(1)*p1(3)/p1(0)/rkt
      ce1(1)%e(2)=p1(2)*p1(3)/p1(0)/rkt
      ce1(1)%e(3)=-rkt/p1(0)
      ce1(2)%e(1)=-p1(2)/rkt
      ce1(2)%e(2)=p1(1)/rkt
      ELSE
      ce1(1)%e(1)=1.d0
      ce1(1)%e(2)=0.d0
      ce1(1)%e(3)=0.d0
      ce1(2)%e(1)=0.d0
      ce1(2)%e(2)=1.d0
      ENDIF
* pol -- p=p2,m=0,e1=ce2(1)%e,e2=ce2(2)%e,e3=0,nhz=0                    
      rkt=sqrt(p2(1)**2+p2(2)**2)
      IF(abs(rkt/p2(0)).GT.1.d-10)THEN
      ce2(1)%e(1)=p2(1)*p2(3)/p2(0)/rkt
      ce2(1)%e(2)=p2(2)*p2(3)/p2(0)/rkt
      ce2(1)%e(3)=-rkt/p2(0)
      ce2(2)%e(1)=-p2(2)/rkt
      ce2(2)%e(2)=p2(1)/rkt
      ELSE
      ce2(1)%e(1)=1.d0
      ce2(1)%e(2)=0.d0
      ce2(1)%e(3)=0.d0
      ce2(2)%e(1)=0.d0
      ce2(2)%e(2)=-1.d0
      ENDIF
  
      do i=1,2
* pk0 -- p=ce1(i)%e
      ce1(i)%ek0=ce1(i)%e(0)-ce1(i)%e(1)
      end do
      do i=1,2
* pk0 -- p=ce2(i)%e
      ce2(i)%ek0=ce2(i)%e(0)-ce2(i)%e(1)
      end do
  
      do mu=0,3
         p12(mu)=p1(mu)+p2(mu)
      enddo
* p.q -- p.q=s12,p=p12,q=p12,bef=,aft=
      s12=(p12(0)*p12(0)-p12(1)*p12(1)-p12(2)*p12(2)-p12(3)*p12(
     & 3))
* triple vertex                                                         
*** triple vertex -- pfz(mu)=p2(mu),pwm(mu)=p1(mu),pwp(mu)=(-p12(mu)),ef
* z=ce2(i2),ewm=ce1(i1),ewp=#,res=ctrip12(i1?,i2?)%e(mu#)
      do mu=0,3
      vfz(mu)=p1(mu)-(-p12(mu))
      vwm(mu)=(-p12(mu))-p2(mu)
      vwp(mu)=p2(mu)-p1(mu)
      end do !mu
* vfz%efz
      do i2=1,2
* p.q -- p.q=ce2(i2)%v,p=ce2(i2)%e,q=vfz,bef=,aft=
      ce2(i2)%v=(ce2(i2)%e(0)*vfz(0)-ce2(i2)%e(1)*vfz(1)-ce2(i2)
     & %e(2)*vfz(2)-ce2(i2)%e(3)*vfz(3))
      end do
* vwm%ewm
      do i1=1,2
* p.q -- p.q=ce1(i1)%v,p=ce1(i1)%e,q=vwm,bef=,aft=
      ce1(i1)%v=(ce1(i1)%e(0)*vwm(0)-ce1(i1)%e(1)*vwm(1)-ce1(i1)
     & %e(2)*vwm(2)-ce1(i1)%e(3)*vwm(3))
      end do
* efz%ewm
      do i2=1,2
      do i1=1,2
* p.q -- p.q=caux,p=ce2(i2)%e,q=ce1(i1)%e,bef=,aft=
      caux=(ce2(i2)%e(0)*ce1(i1)%e(0)-ce2(i2)%e(1)*ce1(i1)%e(1)-
     & ce2(i2)%e(2)*ce1(i1)%e(2)-ce2(i2)%e(3)*ce1(i1)%e(3))
      do mu=0,3
      ctrip12(i1,i2)%e(mu)=ce2(i2)%v*ce1(i1)%e(mu)+ce1(i1)%v*ce2
     & (i2)%e(mu)+vwp(mu)*caux
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
* pk0 -- p=ctrip12(i1,i2)%e
      ctrip12(i1,i2)%ek0=ctrip12(i1,i2)%e(0)-ctrip12(i1,i2)%e(1)
      end do
      end do
  
  
*    -one-gluon fork                                                    
  
  
      do mu=0,3
         p31(mu)=p3(mu)+p1(mu)
*	 p341(mu)=p31(mu)+p4(mu)                                              
      enddo
* pk0 -- p=p31
      p31k0=p31(0)-p31(1)
* p.q -- p.q=s31,p=p31,q=p31,bef=,aft=
      s31=(p31(0)*p31(0)-p31(1)*p31(1)-p31(2)*p31(2)-p31(3)*p31(
     & 3))
      f31=s31*p31k0
  
* p.q -- p.q=s341,p=p31,q=p4,bef=2.d0*,aft=+s31
      s341=2.d0*(p31(0)*p4(0)-p31(1)*p4(1)-p31(2)*p4(2)-p31(3)*p
     & 4(3))+s31
      do mu=0,3
         p41(mu)=-p1(mu)-p4(mu)
      enddo
* pk0 -- p=p41
      p41k0=p41(0)-p41(1)
* p.q -- p.q=s41,p=p41,q=p41,bef=,aft=
      s41=(p41(0)*p41(0)-p41(1)*p41(1)-p41(2)*p41(2)-p41(3)*p41(
     & 3))
      f41=s41*p41k0
  
      if (ilept(id3).ne.1) then
      quqd=s31/2d0
      ccr=1.d0/(f31)
      ccl=1.d0/(f31)
      do i1=1,2
* TL0 -- qu=p3,qd=p31,v=ce1(i1)%e,a=l3_1(i1)%a,c=l3_1(i1)%c,cr=ccr,cl=cc
* l,nsum=0
      ceps_0=-ce1(i1)%ek0*(p3(2)*p31(3)-p31(2)*p3(3))+p3k0*(ce1(
     & i1)%e(2)*p31(3)-p31(2)*ce1(i1)%e(3))-p31k0*(ce1(i1)%e(2)*
     & p3(3)-p3(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p3k0+p3(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=ce1(i1)%e(0)*p3(0)-ce1(i1)%e(1)*p3(1)-ce1(i1)%e(2)*p3
     & (2)-ce1(i1)%e(3)*p3(3)
      cvqd=ce1(i1)%e(0)*p31(0)-ce1(i1)%e(1)*p31(1)-ce1(i1)%e(2)*
     & p31(2)-ce1(i1)%e(3)*p31(3)
      cauxa=-ce1(i1)%ek0*quqd+p3k0*cvqd+p31k0*cvqu
      cauxc=+ce1(i1)%ek0*p3(2)-p3k0*ce1(i1)%e(2)
      l3_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      l3_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      l3_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s41/2d0
      do i1=1,2
* TR0 -- qu=p41,qd=p4,v=ce1(i1)%e,a=r4_1(i1)%a,b=r4_1(i1)%b,cr=1.d0,cl=1
* .d0,nsum=0
      ceps_0=-ce1(i1)%ek0*(p41(2)*p4(3)-p4(2)*p41(3))+p41k0*(ce1
     & (i1)%e(2)*p4(3)-p4(2)*ce1(i1)%e(3))-p4k0*(ce1(i1)%e(2)*p4
     & 1(3)-p41(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-ce1(i1)%e(3)*p4k0+p4(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p41(0)-ce1(i1)%e(1)*p41(1)-ce1(i1)%e(2)*
     & p41(2)-ce1(i1)%e(3)*p41(3)
      cvqd=ce1(i1)%e(0)*p4(0)-ce1(i1)%e(1)*p4(1)-ce1(i1)%e(2)*p4
     & (2)-ce1(i1)%e(3)*p4(3)
      cauxa=-ce1(i1)%ek0*quqd+p41k0*cvqd+p4k0*cvqu
      cauxb=-ce1(i1)%ek0*p4(2)+p4k0*ce1(i1)%e(2)
      r4_1(i1)%a(1)=(cauxa+ceps_0)
      r4_1(i1)%a(2)=(cauxa-ceps_0)
      r4_1(i1)%b(1)=(cauxb-ceps_2)
      r4_1(i1)%b(2)=(-cauxb-ceps_2)
      end do
  
  
      cden= (-s341+cmz2)
* quqd -- p=p31,q=p4
      quqd=p31(0)*p4(0)-p31(1)*p4(1)-p31(2)*p4(2)-p31(3)*p4(3)
      ccr=zcr(id4)/cden
      ccl=zcl(id4)/cden
* TR0 -- qu=p31,qd=p4,v=0,a=rz431(0)%a,b=rz431(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p31(2)*p4(3)+p4(2)*p31(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p31k0*p4(0)+p4k0*p31(0)
      rz431(0)%a(1)=ccr*(auxa+ceps_0)
      rz431(0)%a(2)=ccl*(auxa-ceps_0)
      rz431(0)%b(1)=-ccl*(p4(2)+ceps_2)
      rz431(0)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p31,qd=p4,v=1,a=rz431(1)%a,b=rz431(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p31k0*p4(1)+p4k0*p31(1)
      rz431(1)%a(1)=ccr*(auxa+ceps_0)
      rz431(1)%a(2)=ccl*(auxa-ceps_0)
      rz431(1)%b(1)=-ccl*(p4(2)+ceps_2)
      rz431(1)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p31,qd=p4,v=2,a=rz431(2)%a,b=rz431(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p31k0*p4(3)+p4k0*p31(3)
      ceps_0=eps_0*cim
      auxa=p31k0*p4(2)+p4k0*p31(2)
      rz431(2)%a(1)=ccr*(auxa+ceps_0)
      rz431(2)%a(2)=ccl*(auxa-ceps_0)
      rz431(2)%b(1)=-ccl*p4k0
      rz431(2)%b(2)=ccr*p4k0
* TR0 -- qu=p31,qd=p4,v=3,a=rz431(3)%a,b=rz431(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p31k0*p4(2)-p4k0*p31(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p31k0*p4(3)+p4k0*p31(3)
      rz431(3)%a(1)=ccr*(auxa+ceps_0)
      rz431(3)%a(2)=ccl*(auxa-ceps_0)
      rz431(3)%b(1)=-ccl*ceps_2
      rz431(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz314(&,i1)%e(m),a1=l3_1(i1)%a,c1=l3_1(i1)%c,a2=rz431(m)%a
* ,b2=rz431(m)%b,prq=s31,bef=,aft=
      cz314(1,i1)%e(m)=(l3_1(i1)%a(1)*rz431(m)%a(1)+l3_1(i1)%c(1
     & )*s31*rz431(m)%b(2))
      cz314(2,i1)%e(m)=(l3_1(i1)%c(2)*s31*rz431(m)%b(1)+l3_1(i1)
     & %a(2)*rz431(m)%a(2))
      end do
      end do
  
      cden=(-s341+cmz2)*f41
* quqd -- p=p3,q=p41
      quqd=p3(0)*p41(0)-p3(1)*p41(1)-p3(2)*p41(2)-p3(3)*p41(3)
      ccr=zcr(id3)/cden
      ccl=zcl(id3)/cden
* TL0 -- qu=p3,qd=p41,v=0,a=lz341(0)%a,c=lz341(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p41(3)+p41(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p41(0)+p41k0*p3(0)
      lz341(0)%a(1)=ccr*(auxa+ceps_0)
      lz341(0)%a(2)=ccl*(auxa-ceps_0)
      lz341(0)%c(1)=ccr*(p3(2)+ceps_1)
      lz341(0)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p41,v=1,a=lz341(1)%a,c=lz341(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p41(1)+p41k0*p3(1)
      lz341(1)%a(1)=ccr*(auxa+ceps_0)
      lz341(1)%a(2)=ccl*(auxa-ceps_0)
      lz341(1)%c(1)=ccr*(p3(2)+ceps_1)
      lz341(1)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p41,v=2,a=lz341(2)%a,c=lz341(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p41(3)+p41k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p41(2)+p41k0*p3(2)
      lz341(2)%a(1)=ccr*(auxa+ceps_0)
      lz341(2)%a(2)=ccl*(auxa-ceps_0)
      lz341(2)%c(1)=ccr*p3k0
      lz341(2)%c(2)=-ccl*p3k0
* TL0 -- qu=p3,qd=p41,v=3,a=lz341(3)%a,c=lz341(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p41(2)-p41k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p41(3)+p41k0*p3(3)
      lz341(3)%a(1)=ccr*(auxa+ceps_0)
      lz341(3)%a(2)=ccl*(auxa-ceps_0)
      lz341(3)%c(1)=ccr*ceps_1
      lz341(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz314(&,i1)%e(m),a1=lz341(m)%a,c1=lz341(m)%c,a2=r4_1(i1)%a
* ,b2=r4_1(i1)%b,prq=s41,bef=cz314(&,i1)%e(m)+,aft=
      cz314(1,i1)%e(m)=cz314(1,i1)%e(m)+(lz341(m)%a(1)*r4_1(i1)%
     & a(1)+lz341(m)%c(1)*s41*r4_1(i1)%b(2))
      cz314(2,i1)%e(m)=cz314(2,i1)%e(m)+(lz341(m)%c(2)*s41*r4_1(
     & i1)%b(1)+lz341(m)%a(2)*r4_1(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cz314(i3,i1)%e
      cz314(i3,i1)%ek0=cz314(i3,i1)%e(0)-cz314(i3,i1)%e(1)
      end do
      end do
  
  
      cden= (-s341)
* quqd -- p=p31,q=p4
      quqd=p31(0)*p4(0)-p31(1)*p4(1)-p31(2)*p4(2)-p31(3)*p4(3)
      ccr=fcr(id4)/cden
      ccl=fcl(id4)/cden
* TR0 -- qu=p31,qd=p4,v=0,a=rf431(0)%a,b=rf431(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p31(2)*p4(3)+p4(2)*p31(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p31k0*p4(0)+p4k0*p31(0)
      rf431(0)%a(1)=ccr*(auxa+ceps_0)
      rf431(0)%a(2)=ccl*(auxa-ceps_0)
      rf431(0)%b(1)=-ccl*(p4(2)+ceps_2)
      rf431(0)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p31,qd=p4,v=1,a=rf431(1)%a,b=rf431(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p31k0*p4(1)+p4k0*p31(1)
      rf431(1)%a(1)=ccr*(auxa+ceps_0)
      rf431(1)%a(2)=ccl*(auxa-ceps_0)
      rf431(1)%b(1)=-ccl*(p4(2)+ceps_2)
      rf431(1)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p31,qd=p4,v=2,a=rf431(2)%a,b=rf431(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p31k0*p4(3)+p4k0*p31(3)
      ceps_0=eps_0*cim
      auxa=p31k0*p4(2)+p4k0*p31(2)
      rf431(2)%a(1)=ccr*(auxa+ceps_0)
      rf431(2)%a(2)=ccl*(auxa-ceps_0)
      rf431(2)%b(1)=-ccl*p4k0
      rf431(2)%b(2)=ccr*p4k0
* TR0 -- qu=p31,qd=p4,v=3,a=rf431(3)%a,b=rf431(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p31k0*p4(2)-p4k0*p31(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p31k0*p4(3)+p4k0*p31(3)
      rf431(3)%a(1)=ccr*(auxa+ceps_0)
      rf431(3)%a(2)=ccl*(auxa-ceps_0)
      rf431(3)%b(1)=-ccl*ceps_2
      rf431(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf314(&,i1)%e(m),a1=l3_1(i1)%a,c1=l3_1(i1)%c,a2=rf431(m)%a
* ,b2=rf431(m)%b,prq=s31,bef=,aft=
      cf314(1,i1)%e(m)=(l3_1(i1)%a(1)*rf431(m)%a(1)+l3_1(i1)%c(1
     & )*s31*rf431(m)%b(2))
      cf314(2,i1)%e(m)=(l3_1(i1)%c(2)*s31*rf431(m)%b(1)+l3_1(i1)
     & %a(2)*rf431(m)%a(2))
      end do
      end do
  
      cden=(-s341)*f41
* quqd -- p=p3,q=p41
      quqd=p3(0)*p41(0)-p3(1)*p41(1)-p3(2)*p41(2)-p3(3)*p41(3)
      ccr=fcr(id3)/cden
      ccl=fcl(id3)/cden
* TL0 -- qu=p3,qd=p41,v=0,a=lf341(0)%a,c=lf341(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p41(3)+p41(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p41(0)+p41k0*p3(0)
      lf341(0)%a(1)=ccr*(auxa+ceps_0)
      lf341(0)%a(2)=ccl*(auxa-ceps_0)
      lf341(0)%c(1)=ccr*(p3(2)+ceps_1)
      lf341(0)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p41,v=1,a=lf341(1)%a,c=lf341(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p41(1)+p41k0*p3(1)
      lf341(1)%a(1)=ccr*(auxa+ceps_0)
      lf341(1)%a(2)=ccl*(auxa-ceps_0)
      lf341(1)%c(1)=ccr*(p3(2)+ceps_1)
      lf341(1)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p41,v=2,a=lf341(2)%a,c=lf341(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p41(3)+p41k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p41(2)+p41k0*p3(2)
      lf341(2)%a(1)=ccr*(auxa+ceps_0)
      lf341(2)%a(2)=ccl*(auxa-ceps_0)
      lf341(2)%c(1)=ccr*p3k0
      lf341(2)%c(2)=-ccl*p3k0
* TL0 -- qu=p3,qd=p41,v=3,a=lf341(3)%a,c=lf341(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p41(2)-p41k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p41(3)+p41k0*p3(3)
      lf341(3)%a(1)=ccr*(auxa+ceps_0)
      lf341(3)%a(2)=ccl*(auxa-ceps_0)
      lf341(3)%c(1)=ccr*ceps_1
      lf341(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf314(&,i1)%e(m),a1=lf341(m)%a,c1=lf341(m)%c,a2=r4_1(i1)%a
* ,b2=r4_1(i1)%b,prq=s41,bef=cf314(&,i1)%e(m)+,aft=
      cf314(1,i1)%e(m)=cf314(1,i1)%e(m)+(lf341(m)%a(1)*r4_1(i1)%
     & a(1)+lf341(m)%c(1)*s41*r4_1(i1)%b(2))
      cf314(2,i1)%e(m)=cf314(2,i1)%e(m)+(lf341(m)%c(2)*s41*r4_1(
     & i1)%b(1)+lf341(m)%a(2)*r4_1(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cf314(i3,i1)%e
      cf314(i3,i1)%ek0=cf314(i3,i1)%e(0)-cf314(i3,i1)%e(1)
      end do
      end do
  
      endif  ! id3 = quark
  
      do mu=0,3
         p32(mu)=p3(mu)+p2(mu)
*	 p342(mu)=p32(mu)+p4(mu)                                              
      enddo
* pk0 -- p=p32
      p32k0=p32(0)-p32(1)
* p.q -- p.q=s32,p=p32,q=p32,bef=,aft=
      s32=(p32(0)*p32(0)-p32(1)*p32(1)-p32(2)*p32(2)-p32(3)*p32(
     & 3))
      f32=s32*p32k0
  
* p.q -- p.q=s342,p=p32,q=p4,bef=2.d0*,aft=+s32
      s342=2.d0*(p32(0)*p4(0)-p32(1)*p4(1)-p32(2)*p4(2)-p32(3)*p
     & 4(3))+s32
      do mu=0,3
         p42(mu)=-p2(mu)-p4(mu)
      enddo
* pk0 -- p=p42
      p42k0=p42(0)-p42(1)
* p.q -- p.q=s42,p=p42,q=p42,bef=,aft=
      s42=(p42(0)*p42(0)-p42(1)*p42(1)-p42(2)*p42(2)-p42(3)*p42(
     & 3))
      f42=s42*p42k0
  
      if (ilept(id3).ne.1) then
      quqd=s32/2d0
      ccr=1.d0/(f32)
      ccl=1.d0/(f32)
      do i1=1,2
* TL0 -- qu=p3,qd=p32,v=ce2(i1)%e,a=l3_2(i1)%a,c=l3_2(i1)%c,cr=ccr,cl=cc
* l,nsum=0
      ceps_0=-ce2(i1)%ek0*(p3(2)*p32(3)-p32(2)*p3(3))+p3k0*(ce2(
     & i1)%e(2)*p32(3)-p32(2)*ce2(i1)%e(3))-p32k0*(ce2(i1)%e(2)*
     & p3(3)-p3(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p3k0+p3(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=ce2(i1)%e(0)*p3(0)-ce2(i1)%e(1)*p3(1)-ce2(i1)%e(2)*p3
     & (2)-ce2(i1)%e(3)*p3(3)
      cvqd=ce2(i1)%e(0)*p32(0)-ce2(i1)%e(1)*p32(1)-ce2(i1)%e(2)*
     & p32(2)-ce2(i1)%e(3)*p32(3)
      cauxa=-ce2(i1)%ek0*quqd+p3k0*cvqd+p32k0*cvqu
      cauxc=+ce2(i1)%ek0*p3(2)-p3k0*ce2(i1)%e(2)
      l3_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      l3_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      l3_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s42/2d0
      do i1=1,2
* TR0 -- qu=p42,qd=p4,v=ce2(i1)%e,a=r4_2(i1)%a,b=r4_2(i1)%b,cr=1.d0,cl=1
* .d0,nsum=0
      ceps_0=-ce2(i1)%ek0*(p42(2)*p4(3)-p4(2)*p42(3))+p42k0*(ce2
     & (i1)%e(2)*p4(3)-p4(2)*ce2(i1)%e(3))-p4k0*(ce2(i1)%e(2)*p4
     & 2(3)-p42(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-ce2(i1)%e(3)*p4k0+p4(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p42(0)-ce2(i1)%e(1)*p42(1)-ce2(i1)%e(2)*
     & p42(2)-ce2(i1)%e(3)*p42(3)
      cvqd=ce2(i1)%e(0)*p4(0)-ce2(i1)%e(1)*p4(1)-ce2(i1)%e(2)*p4
     & (2)-ce2(i1)%e(3)*p4(3)
      cauxa=-ce2(i1)%ek0*quqd+p42k0*cvqd+p4k0*cvqu
      cauxb=-ce2(i1)%ek0*p4(2)+p4k0*ce2(i1)%e(2)
      r4_2(i1)%a(1)=(cauxa+ceps_0)
      r4_2(i1)%a(2)=(cauxa-ceps_0)
      r4_2(i1)%b(1)=(cauxb-ceps_2)
      r4_2(i1)%b(2)=(-cauxb-ceps_2)
      end do
  
  
      cden= (-s342+cmz2)
* quqd -- p=p32,q=p4
      quqd=p32(0)*p4(0)-p32(1)*p4(1)-p32(2)*p4(2)-p32(3)*p4(3)
      ccr=zcr(id4)/cden
      ccl=zcl(id4)/cden
* TR0 -- qu=p32,qd=p4,v=0,a=rz432(0)%a,b=rz432(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p32(2)*p4(3)+p4(2)*p32(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p32k0*p4(0)+p4k0*p32(0)
      rz432(0)%a(1)=ccr*(auxa+ceps_0)
      rz432(0)%a(2)=ccl*(auxa-ceps_0)
      rz432(0)%b(1)=-ccl*(p4(2)+ceps_2)
      rz432(0)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p32,qd=p4,v=1,a=rz432(1)%a,b=rz432(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p32k0*p4(1)+p4k0*p32(1)
      rz432(1)%a(1)=ccr*(auxa+ceps_0)
      rz432(1)%a(2)=ccl*(auxa-ceps_0)
      rz432(1)%b(1)=-ccl*(p4(2)+ceps_2)
      rz432(1)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p32,qd=p4,v=2,a=rz432(2)%a,b=rz432(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p32k0*p4(3)+p4k0*p32(3)
      ceps_0=eps_0*cim
      auxa=p32k0*p4(2)+p4k0*p32(2)
      rz432(2)%a(1)=ccr*(auxa+ceps_0)
      rz432(2)%a(2)=ccl*(auxa-ceps_0)
      rz432(2)%b(1)=-ccl*p4k0
      rz432(2)%b(2)=ccr*p4k0
* TR0 -- qu=p32,qd=p4,v=3,a=rz432(3)%a,b=rz432(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p32k0*p4(2)-p4k0*p32(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p32k0*p4(3)+p4k0*p32(3)
      rz432(3)%a(1)=ccr*(auxa+ceps_0)
      rz432(3)%a(2)=ccl*(auxa-ceps_0)
      rz432(3)%b(1)=-ccl*ceps_2
      rz432(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz324(&,i1)%e(m),a1=l3_2(i1)%a,c1=l3_2(i1)%c,a2=rz432(m)%a
* ,b2=rz432(m)%b,prq=s32,bef=,aft=
      cz324(1,i1)%e(m)=(l3_2(i1)%a(1)*rz432(m)%a(1)+l3_2(i1)%c(1
     & )*s32*rz432(m)%b(2))
      cz324(2,i1)%e(m)=(l3_2(i1)%c(2)*s32*rz432(m)%b(1)+l3_2(i1)
     & %a(2)*rz432(m)%a(2))
      end do
      end do
  
      cden=(-s342+cmz2)*f42
* quqd -- p=p3,q=p42
      quqd=p3(0)*p42(0)-p3(1)*p42(1)-p3(2)*p42(2)-p3(3)*p42(3)
      ccr=zcr(id3)/cden
      ccl=zcl(id3)/cden
* TL0 -- qu=p3,qd=p42,v=0,a=lz342(0)%a,c=lz342(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p42(3)+p42(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p42(0)+p42k0*p3(0)
      lz342(0)%a(1)=ccr*(auxa+ceps_0)
      lz342(0)%a(2)=ccl*(auxa-ceps_0)
      lz342(0)%c(1)=ccr*(p3(2)+ceps_1)
      lz342(0)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p42,v=1,a=lz342(1)%a,c=lz342(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p42(1)+p42k0*p3(1)
      lz342(1)%a(1)=ccr*(auxa+ceps_0)
      lz342(1)%a(2)=ccl*(auxa-ceps_0)
      lz342(1)%c(1)=ccr*(p3(2)+ceps_1)
      lz342(1)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p42,v=2,a=lz342(2)%a,c=lz342(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p42(3)+p42k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p42(2)+p42k0*p3(2)
      lz342(2)%a(1)=ccr*(auxa+ceps_0)
      lz342(2)%a(2)=ccl*(auxa-ceps_0)
      lz342(2)%c(1)=ccr*p3k0
      lz342(2)%c(2)=-ccl*p3k0
* TL0 -- qu=p3,qd=p42,v=3,a=lz342(3)%a,c=lz342(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p42(2)-p42k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p42(3)+p42k0*p3(3)
      lz342(3)%a(1)=ccr*(auxa+ceps_0)
      lz342(3)%a(2)=ccl*(auxa-ceps_0)
      lz342(3)%c(1)=ccr*ceps_1
      lz342(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz324(&,i1)%e(m),a1=lz342(m)%a,c1=lz342(m)%c,a2=r4_2(i1)%a
* ,b2=r4_2(i1)%b,prq=s42,bef=cz324(&,i1)%e(m)+,aft=
      cz324(1,i1)%e(m)=cz324(1,i1)%e(m)+(lz342(m)%a(1)*r4_2(i1)%
     & a(1)+lz342(m)%c(1)*s42*r4_2(i1)%b(2))
      cz324(2,i1)%e(m)=cz324(2,i1)%e(m)+(lz342(m)%c(2)*s42*r4_2(
     & i1)%b(1)+lz342(m)%a(2)*r4_2(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cz324(i3,i1)%e
      cz324(i3,i1)%ek0=cz324(i3,i1)%e(0)-cz324(i3,i1)%e(1)
      end do
      end do
  
  
      cden= (-s342)
* quqd -- p=p32,q=p4
      quqd=p32(0)*p4(0)-p32(1)*p4(1)-p32(2)*p4(2)-p32(3)*p4(3)
      ccr=fcr(id4)/cden
      ccl=fcl(id4)/cden
* TR0 -- qu=p32,qd=p4,v=0,a=rf432(0)%a,b=rf432(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p32(2)*p4(3)+p4(2)*p32(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p32k0*p4(0)+p4k0*p32(0)
      rf432(0)%a(1)=ccr*(auxa+ceps_0)
      rf432(0)%a(2)=ccl*(auxa-ceps_0)
      rf432(0)%b(1)=-ccl*(p4(2)+ceps_2)
      rf432(0)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p32,qd=p4,v=1,a=rf432(1)%a,b=rf432(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p32k0*p4(1)+p4k0*p32(1)
      rf432(1)%a(1)=ccr*(auxa+ceps_0)
      rf432(1)%a(2)=ccl*(auxa-ceps_0)
      rf432(1)%b(1)=-ccl*(p4(2)+ceps_2)
      rf432(1)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p32,qd=p4,v=2,a=rf432(2)%a,b=rf432(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p32k0*p4(3)+p4k0*p32(3)
      ceps_0=eps_0*cim
      auxa=p32k0*p4(2)+p4k0*p32(2)
      rf432(2)%a(1)=ccr*(auxa+ceps_0)
      rf432(2)%a(2)=ccl*(auxa-ceps_0)
      rf432(2)%b(1)=-ccl*p4k0
      rf432(2)%b(2)=ccr*p4k0
* TR0 -- qu=p32,qd=p4,v=3,a=rf432(3)%a,b=rf432(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p32k0*p4(2)-p4k0*p32(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p32k0*p4(3)+p4k0*p32(3)
      rf432(3)%a(1)=ccr*(auxa+ceps_0)
      rf432(3)%a(2)=ccl*(auxa-ceps_0)
      rf432(3)%b(1)=-ccl*ceps_2
      rf432(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf324(&,i1)%e(m),a1=l3_2(i1)%a,c1=l3_2(i1)%c,a2=rf432(m)%a
* ,b2=rf432(m)%b,prq=s32,bef=,aft=
      cf324(1,i1)%e(m)=(l3_2(i1)%a(1)*rf432(m)%a(1)+l3_2(i1)%c(1
     & )*s32*rf432(m)%b(2))
      cf324(2,i1)%e(m)=(l3_2(i1)%c(2)*s32*rf432(m)%b(1)+l3_2(i1)
     & %a(2)*rf432(m)%a(2))
      end do
      end do
  
      cden=(-s342)*f42
* quqd -- p=p3,q=p42
      quqd=p3(0)*p42(0)-p3(1)*p42(1)-p3(2)*p42(2)-p3(3)*p42(3)
      ccr=fcr(id3)/cden
      ccl=fcl(id3)/cden
* TL0 -- qu=p3,qd=p42,v=0,a=lf342(0)%a,c=lf342(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p42(3)+p42(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p42(0)+p42k0*p3(0)
      lf342(0)%a(1)=ccr*(auxa+ceps_0)
      lf342(0)%a(2)=ccl*(auxa-ceps_0)
      lf342(0)%c(1)=ccr*(p3(2)+ceps_1)
      lf342(0)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p42,v=1,a=lf342(1)%a,c=lf342(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p42(1)+p42k0*p3(1)
      lf342(1)%a(1)=ccr*(auxa+ceps_0)
      lf342(1)%a(2)=ccl*(auxa-ceps_0)
      lf342(1)%c(1)=ccr*(p3(2)+ceps_1)
      lf342(1)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p42,v=2,a=lf342(2)%a,c=lf342(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p42(3)+p42k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p42(2)+p42k0*p3(2)
      lf342(2)%a(1)=ccr*(auxa+ceps_0)
      lf342(2)%a(2)=ccl*(auxa-ceps_0)
      lf342(2)%c(1)=ccr*p3k0
      lf342(2)%c(2)=-ccl*p3k0
* TL0 -- qu=p3,qd=p42,v=3,a=lf342(3)%a,c=lf342(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p42(2)-p42k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p42(3)+p42k0*p3(3)
      lf342(3)%a(1)=ccr*(auxa+ceps_0)
      lf342(3)%a(2)=ccl*(auxa-ceps_0)
      lf342(3)%c(1)=ccr*ceps_1
      lf342(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf324(&,i1)%e(m),a1=lf342(m)%a,c1=lf342(m)%c,a2=r4_2(i1)%a
* ,b2=r4_2(i1)%b,prq=s42,bef=cf324(&,i1)%e(m)+,aft=
      cf324(1,i1)%e(m)=cf324(1,i1)%e(m)+(lf342(m)%a(1)*r4_2(i1)%
     & a(1)+lf342(m)%c(1)*s42*r4_2(i1)%b(2))
      cf324(2,i1)%e(m)=cf324(2,i1)%e(m)+(lf342(m)%c(2)*s42*r4_2(
     & i1)%b(1)+lf342(m)%a(2)*r4_2(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cf324(i3,i1)%e
      cf324(i3,i1)%ek0=cf324(i3,i1)%e(0)-cf324(i3,i1)%e(1)
      end do
      end do
  
      endif  ! id3 = quark
  
  
      do mu=0,3
         p51(mu)=p5(mu)+p1(mu)
*	 p561(mu)=p51(mu)+p6(mu)                                              
      enddo
* pk0 -- p=p51
      p51k0=p51(0)-p51(1)
* p.q -- p.q=s51,p=p51,q=p51,bef=,aft=
      s51=(p51(0)*p51(0)-p51(1)*p51(1)-p51(2)*p51(2)-p51(3)*p51(
     & 3))
      f51=s51*p51k0
  
* p.q -- p.q=s561,p=p51,q=p6,bef=2.d0*,aft=+s51
      s561=2.d0*(p51(0)*p6(0)-p51(1)*p6(1)-p51(2)*p6(2)-p51(3)*p
     & 6(3))+s51
      do mu=0,3
         p61(mu)=-p1(mu)-p6(mu)
      enddo
* pk0 -- p=p61
      p61k0=p61(0)-p61(1)
* p.q -- p.q=s61,p=p61,q=p61,bef=,aft=
      s61=(p61(0)*p61(0)-p61(1)*p61(1)-p61(2)*p61(2)-p61(3)*p61(
     & 3))
      f61=s61*p61k0
  
      if (ilept(id5).ne.1) then
      quqd=s51/2d0
      ccr=1.d0/(f51)
      ccl=1.d0/(f51)
      do i1=1,2
* TL0 -- qu=p5,qd=p51,v=ce1(i1)%e,a=l5_1(i1)%a,c=l5_1(i1)%c,cr=ccr,cl=cc
* l,nsum=0
      ceps_0=-ce1(i1)%ek0*(p5(2)*p51(3)-p51(2)*p5(3))+p5k0*(ce1(
     & i1)%e(2)*p51(3)-p51(2)*ce1(i1)%e(3))-p51k0*(ce1(i1)%e(2)*
     & p5(3)-p5(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p5k0+p5(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=ce1(i1)%e(0)*p5(0)-ce1(i1)%e(1)*p5(1)-ce1(i1)%e(2)*p5
     & (2)-ce1(i1)%e(3)*p5(3)
      cvqd=ce1(i1)%e(0)*p51(0)-ce1(i1)%e(1)*p51(1)-ce1(i1)%e(2)*
     & p51(2)-ce1(i1)%e(3)*p51(3)
      cauxa=-ce1(i1)%ek0*quqd+p5k0*cvqd+p51k0*cvqu
      cauxc=+ce1(i1)%ek0*p5(2)-p5k0*ce1(i1)%e(2)
      l5_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      l5_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      l5_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s61/2d0
      do i1=1,2
* TR0 -- qu=p61,qd=p6,v=ce1(i1)%e,a=r6_1(i1)%a,b=r6_1(i1)%b,cr=1.d0,cl=1
* .d0,nsum=0
      ceps_0=-ce1(i1)%ek0*(p61(2)*p6(3)-p6(2)*p61(3))+p61k0*(ce1
     & (i1)%e(2)*p6(3)-p6(2)*ce1(i1)%e(3))-p6k0*(ce1(i1)%e(2)*p6
     & 1(3)-p61(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-ce1(i1)%e(3)*p6k0+p6(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p61(0)-ce1(i1)%e(1)*p61(1)-ce1(i1)%e(2)*
     & p61(2)-ce1(i1)%e(3)*p61(3)
      cvqd=ce1(i1)%e(0)*p6(0)-ce1(i1)%e(1)*p6(1)-ce1(i1)%e(2)*p6
     & (2)-ce1(i1)%e(3)*p6(3)
      cauxa=-ce1(i1)%ek0*quqd+p61k0*cvqd+p6k0*cvqu
      cauxb=-ce1(i1)%ek0*p6(2)+p6k0*ce1(i1)%e(2)
      r6_1(i1)%a(1)=(cauxa+ceps_0)
      r6_1(i1)%a(2)=(cauxa-ceps_0)
      r6_1(i1)%b(1)=(cauxb-ceps_2)
      r6_1(i1)%b(2)=(-cauxb-ceps_2)
      end do
  
  
      cden= (-s561+cmz2)
* quqd -- p=p51,q=p6
      quqd=p51(0)*p6(0)-p51(1)*p6(1)-p51(2)*p6(2)-p51(3)*p6(3)
      ccr=zcr(id6)/cden
      ccl=zcl(id6)/cden
* TR0 -- qu=p51,qd=p6,v=0,a=rz651(0)%a,b=rz651(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p51(2)*p6(3)+p6(2)*p51(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p51k0*p6(0)+p6k0*p51(0)
      rz651(0)%a(1)=ccr*(auxa+ceps_0)
      rz651(0)%a(2)=ccl*(auxa-ceps_0)
      rz651(0)%b(1)=-ccl*(p6(2)+ceps_2)
      rz651(0)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p51,qd=p6,v=1,a=rz651(1)%a,b=rz651(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p51k0*p6(1)+p6k0*p51(1)
      rz651(1)%a(1)=ccr*(auxa+ceps_0)
      rz651(1)%a(2)=ccl*(auxa-ceps_0)
      rz651(1)%b(1)=-ccl*(p6(2)+ceps_2)
      rz651(1)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p51,qd=p6,v=2,a=rz651(2)%a,b=rz651(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p51k0*p6(3)+p6k0*p51(3)
      ceps_0=eps_0*cim
      auxa=p51k0*p6(2)+p6k0*p51(2)
      rz651(2)%a(1)=ccr*(auxa+ceps_0)
      rz651(2)%a(2)=ccl*(auxa-ceps_0)
      rz651(2)%b(1)=-ccl*p6k0
      rz651(2)%b(2)=ccr*p6k0
* TR0 -- qu=p51,qd=p6,v=3,a=rz651(3)%a,b=rz651(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p51k0*p6(2)-p6k0*p51(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p51k0*p6(3)+p6k0*p51(3)
      rz651(3)%a(1)=ccr*(auxa+ceps_0)
      rz651(3)%a(2)=ccl*(auxa-ceps_0)
      rz651(3)%b(1)=-ccl*ceps_2
      rz651(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz516(&,i1)%e(m),a1=l5_1(i1)%a,c1=l5_1(i1)%c,a2=rz651(m)%a
* ,b2=rz651(m)%b,prq=s51,bef=,aft=
      cz516(1,i1)%e(m)=(l5_1(i1)%a(1)*rz651(m)%a(1)+l5_1(i1)%c(1
     & )*s51*rz651(m)%b(2))
      cz516(2,i1)%e(m)=(l5_1(i1)%c(2)*s51*rz651(m)%b(1)+l5_1(i1)
     & %a(2)*rz651(m)%a(2))
      end do
      end do
  
      cden=(-s561+cmz2)*f61
* quqd -- p=p5,q=p61
      quqd=p5(0)*p61(0)-p5(1)*p61(1)-p5(2)*p61(2)-p5(3)*p61(3)
      ccr=zcr(id5)/cden
      ccl=zcl(id5)/cden
* TL0 -- qu=p5,qd=p61,v=0,a=lz561(0)%a,c=lz561(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p5(2)*p61(3)+p61(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p61(0)+p61k0*p5(0)
      lz561(0)%a(1)=ccr*(auxa+ceps_0)
      lz561(0)%a(2)=ccl*(auxa-ceps_0)
      lz561(0)%c(1)=ccr*(p5(2)+ceps_1)
      lz561(0)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p61,v=1,a=lz561(1)%a,c=lz561(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p61(1)+p61k0*p5(1)
      lz561(1)%a(1)=ccr*(auxa+ceps_0)
      lz561(1)%a(2)=ccl*(auxa-ceps_0)
      lz561(1)%c(1)=ccr*(p5(2)+ceps_1)
      lz561(1)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p61,v=2,a=lz561(2)%a,c=lz561(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p5k0*p61(3)+p61k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p61(2)+p61k0*p5(2)
      lz561(2)%a(1)=ccr*(auxa+ceps_0)
      lz561(2)%a(2)=ccl*(auxa-ceps_0)
      lz561(2)%c(1)=ccr*p5k0
      lz561(2)%c(2)=-ccl*p5k0
* TL0 -- qu=p5,qd=p61,v=3,a=lz561(3)%a,c=lz561(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p5k0*p61(2)-p61k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p61(3)+p61k0*p5(3)
      lz561(3)%a(1)=ccr*(auxa+ceps_0)
      lz561(3)%a(2)=ccl*(auxa-ceps_0)
      lz561(3)%c(1)=ccr*ceps_1
      lz561(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz516(&,i1)%e(m),a1=lz561(m)%a,c1=lz561(m)%c,a2=r6_1(i1)%a
* ,b2=r6_1(i1)%b,prq=s61,bef=cz516(&,i1)%e(m)+,aft=
      cz516(1,i1)%e(m)=cz516(1,i1)%e(m)+(lz561(m)%a(1)*r6_1(i1)%
     & a(1)+lz561(m)%c(1)*s61*r6_1(i1)%b(2))
      cz516(2,i1)%e(m)=cz516(2,i1)%e(m)+(lz561(m)%c(2)*s61*r6_1(
     & i1)%b(1)+lz561(m)%a(2)*r6_1(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cz516(i3,i1)%e
      cz516(i3,i1)%ek0=cz516(i3,i1)%e(0)-cz516(i3,i1)%e(1)
      end do
      end do
  
  
      cden= (-s561)
* quqd -- p=p51,q=p6
      quqd=p51(0)*p6(0)-p51(1)*p6(1)-p51(2)*p6(2)-p51(3)*p6(3)
      ccr=fcr(id6)/cden
      ccl=fcl(id6)/cden
* TR0 -- qu=p51,qd=p6,v=0,a=rf651(0)%a,b=rf651(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p51(2)*p6(3)+p6(2)*p51(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p51k0*p6(0)+p6k0*p51(0)
      rf651(0)%a(1)=ccr*(auxa+ceps_0)
      rf651(0)%a(2)=ccl*(auxa-ceps_0)
      rf651(0)%b(1)=-ccl*(p6(2)+ceps_2)
      rf651(0)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p51,qd=p6,v=1,a=rf651(1)%a,b=rf651(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p51k0*p6(1)+p6k0*p51(1)
      rf651(1)%a(1)=ccr*(auxa+ceps_0)
      rf651(1)%a(2)=ccl*(auxa-ceps_0)
      rf651(1)%b(1)=-ccl*(p6(2)+ceps_2)
      rf651(1)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p51,qd=p6,v=2,a=rf651(2)%a,b=rf651(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p51k0*p6(3)+p6k0*p51(3)
      ceps_0=eps_0*cim
      auxa=p51k0*p6(2)+p6k0*p51(2)
      rf651(2)%a(1)=ccr*(auxa+ceps_0)
      rf651(2)%a(2)=ccl*(auxa-ceps_0)
      rf651(2)%b(1)=-ccl*p6k0
      rf651(2)%b(2)=ccr*p6k0
* TR0 -- qu=p51,qd=p6,v=3,a=rf651(3)%a,b=rf651(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p51k0*p6(2)-p6k0*p51(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p51k0*p6(3)+p6k0*p51(3)
      rf651(3)%a(1)=ccr*(auxa+ceps_0)
      rf651(3)%a(2)=ccl*(auxa-ceps_0)
      rf651(3)%b(1)=-ccl*ceps_2
      rf651(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf516(&,i1)%e(m),a1=l5_1(i1)%a,c1=l5_1(i1)%c,a2=rf651(m)%a
* ,b2=rf651(m)%b,prq=s51,bef=,aft=
      cf516(1,i1)%e(m)=(l5_1(i1)%a(1)*rf651(m)%a(1)+l5_1(i1)%c(1
     & )*s51*rf651(m)%b(2))
      cf516(2,i1)%e(m)=(l5_1(i1)%c(2)*s51*rf651(m)%b(1)+l5_1(i1)
     & %a(2)*rf651(m)%a(2))
      end do
      end do
  
      cden=(-s561)*f61
* quqd -- p=p5,q=p61
      quqd=p5(0)*p61(0)-p5(1)*p61(1)-p5(2)*p61(2)-p5(3)*p61(3)
      ccr=fcr(id5)/cden
      ccl=fcl(id5)/cden
* TL0 -- qu=p5,qd=p61,v=0,a=lf561(0)%a,c=lf561(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p5(2)*p61(3)+p61(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p61(0)+p61k0*p5(0)
      lf561(0)%a(1)=ccr*(auxa+ceps_0)
      lf561(0)%a(2)=ccl*(auxa-ceps_0)
      lf561(0)%c(1)=ccr*(p5(2)+ceps_1)
      lf561(0)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p61,v=1,a=lf561(1)%a,c=lf561(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p61(1)+p61k0*p5(1)
      lf561(1)%a(1)=ccr*(auxa+ceps_0)
      lf561(1)%a(2)=ccl*(auxa-ceps_0)
      lf561(1)%c(1)=ccr*(p5(2)+ceps_1)
      lf561(1)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p61,v=2,a=lf561(2)%a,c=lf561(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p5k0*p61(3)+p61k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p61(2)+p61k0*p5(2)
      lf561(2)%a(1)=ccr*(auxa+ceps_0)
      lf561(2)%a(2)=ccl*(auxa-ceps_0)
      lf561(2)%c(1)=ccr*p5k0
      lf561(2)%c(2)=-ccl*p5k0
* TL0 -- qu=p5,qd=p61,v=3,a=lf561(3)%a,c=lf561(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p5k0*p61(2)-p61k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p61(3)+p61k0*p5(3)
      lf561(3)%a(1)=ccr*(auxa+ceps_0)
      lf561(3)%a(2)=ccl*(auxa-ceps_0)
      lf561(3)%c(1)=ccr*ceps_1
      lf561(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf516(&,i1)%e(m),a1=lf561(m)%a,c1=lf561(m)%c,a2=r6_1(i1)%a
* ,b2=r6_1(i1)%b,prq=s61,bef=cf516(&,i1)%e(m)+,aft=
      cf516(1,i1)%e(m)=cf516(1,i1)%e(m)+(lf561(m)%a(1)*r6_1(i1)%
     & a(1)+lf561(m)%c(1)*s61*r6_1(i1)%b(2))
      cf516(2,i1)%e(m)=cf516(2,i1)%e(m)+(lf561(m)%c(2)*s61*r6_1(
     & i1)%b(1)+lf561(m)%a(2)*r6_1(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cf516(i3,i1)%e
      cf516(i3,i1)%ek0=cf516(i3,i1)%e(0)-cf516(i3,i1)%e(1)
      end do
      end do
  
      endif  ! id5 = quark
  
      do mu=0,3
         p52(mu)=p5(mu)+p2(mu)
*	 p562(mu)=p52(mu)+p6(mu)                                              
      enddo
* pk0 -- p=p52
      p52k0=p52(0)-p52(1)
* p.q -- p.q=s52,p=p52,q=p52,bef=,aft=
      s52=(p52(0)*p52(0)-p52(1)*p52(1)-p52(2)*p52(2)-p52(3)*p52(
     & 3))
      f52=s52*p52k0
  
* p.q -- p.q=s562,p=p52,q=p6,bef=2.d0*,aft=+s52
      s562=2.d0*(p52(0)*p6(0)-p52(1)*p6(1)-p52(2)*p6(2)-p52(3)*p
     & 6(3))+s52
      do mu=0,3
         p62(mu)=-p2(mu)-p6(mu)
      enddo
* pk0 -- p=p62
      p62k0=p62(0)-p62(1)
* p.q -- p.q=s62,p=p62,q=p62,bef=,aft=
      s62=(p62(0)*p62(0)-p62(1)*p62(1)-p62(2)*p62(2)-p62(3)*p62(
     & 3))
      f62=s62*p62k0
  
      if (ilept(id5).ne.1) then
      quqd=s52/2d0
      ccr=1.d0/(f52)
      ccl=1.d0/(f52)
      do i1=1,2
* TL0 -- qu=p5,qd=p52,v=ce2(i1)%e,a=l5_2(i1)%a,c=l5_2(i1)%c,cr=ccr,cl=cc
* l,nsum=0
      ceps_0=-ce2(i1)%ek0*(p5(2)*p52(3)-p52(2)*p5(3))+p5k0*(ce2(
     & i1)%e(2)*p52(3)-p52(2)*ce2(i1)%e(3))-p52k0*(ce2(i1)%e(2)*
     & p5(3)-p5(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p5k0+p5(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=ce2(i1)%e(0)*p5(0)-ce2(i1)%e(1)*p5(1)-ce2(i1)%e(2)*p5
     & (2)-ce2(i1)%e(3)*p5(3)
      cvqd=ce2(i1)%e(0)*p52(0)-ce2(i1)%e(1)*p52(1)-ce2(i1)%e(2)*
     & p52(2)-ce2(i1)%e(3)*p52(3)
      cauxa=-ce2(i1)%ek0*quqd+p5k0*cvqd+p52k0*cvqu
      cauxc=+ce2(i1)%ek0*p5(2)-p5k0*ce2(i1)%e(2)
      l5_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      l5_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      l5_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s62/2d0
      do i1=1,2
* TR0 -- qu=p62,qd=p6,v=ce2(i1)%e,a=r6_2(i1)%a,b=r6_2(i1)%b,cr=1.d0,cl=1
* .d0,nsum=0
      ceps_0=-ce2(i1)%ek0*(p62(2)*p6(3)-p6(2)*p62(3))+p62k0*(ce2
     & (i1)%e(2)*p6(3)-p6(2)*ce2(i1)%e(3))-p6k0*(ce2(i1)%e(2)*p6
     & 2(3)-p62(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-ce2(i1)%e(3)*p6k0+p6(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p62(0)-ce2(i1)%e(1)*p62(1)-ce2(i1)%e(2)*
     & p62(2)-ce2(i1)%e(3)*p62(3)
      cvqd=ce2(i1)%e(0)*p6(0)-ce2(i1)%e(1)*p6(1)-ce2(i1)%e(2)*p6
     & (2)-ce2(i1)%e(3)*p6(3)
      cauxa=-ce2(i1)%ek0*quqd+p62k0*cvqd+p6k0*cvqu
      cauxb=-ce2(i1)%ek0*p6(2)+p6k0*ce2(i1)%e(2)
      r6_2(i1)%a(1)=(cauxa+ceps_0)
      r6_2(i1)%a(2)=(cauxa-ceps_0)
      r6_2(i1)%b(1)=(cauxb-ceps_2)
      r6_2(i1)%b(2)=(-cauxb-ceps_2)
      end do
  
  
      cden= (-s562+cmz2)
* quqd -- p=p52,q=p6
      quqd=p52(0)*p6(0)-p52(1)*p6(1)-p52(2)*p6(2)-p52(3)*p6(3)
      ccr=zcr(id6)/cden
      ccl=zcl(id6)/cden
* TR0 -- qu=p52,qd=p6,v=0,a=rz652(0)%a,b=rz652(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p52(2)*p6(3)+p6(2)*p52(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p52k0*p6(0)+p6k0*p52(0)
      rz652(0)%a(1)=ccr*(auxa+ceps_0)
      rz652(0)%a(2)=ccl*(auxa-ceps_0)
      rz652(0)%b(1)=-ccl*(p6(2)+ceps_2)
      rz652(0)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p52,qd=p6,v=1,a=rz652(1)%a,b=rz652(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p52k0*p6(1)+p6k0*p52(1)
      rz652(1)%a(1)=ccr*(auxa+ceps_0)
      rz652(1)%a(2)=ccl*(auxa-ceps_0)
      rz652(1)%b(1)=-ccl*(p6(2)+ceps_2)
      rz652(1)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p52,qd=p6,v=2,a=rz652(2)%a,b=rz652(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p52k0*p6(3)+p6k0*p52(3)
      ceps_0=eps_0*cim
      auxa=p52k0*p6(2)+p6k0*p52(2)
      rz652(2)%a(1)=ccr*(auxa+ceps_0)
      rz652(2)%a(2)=ccl*(auxa-ceps_0)
      rz652(2)%b(1)=-ccl*p6k0
      rz652(2)%b(2)=ccr*p6k0
* TR0 -- qu=p52,qd=p6,v=3,a=rz652(3)%a,b=rz652(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p52k0*p6(2)-p6k0*p52(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p52k0*p6(3)+p6k0*p52(3)
      rz652(3)%a(1)=ccr*(auxa+ceps_0)
      rz652(3)%a(2)=ccl*(auxa-ceps_0)
      rz652(3)%b(1)=-ccl*ceps_2
      rz652(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz526(&,i1)%e(m),a1=l5_2(i1)%a,c1=l5_2(i1)%c,a2=rz652(m)%a
* ,b2=rz652(m)%b,prq=s52,bef=,aft=
      cz526(1,i1)%e(m)=(l5_2(i1)%a(1)*rz652(m)%a(1)+l5_2(i1)%c(1
     & )*s52*rz652(m)%b(2))
      cz526(2,i1)%e(m)=(l5_2(i1)%c(2)*s52*rz652(m)%b(1)+l5_2(i1)
     & %a(2)*rz652(m)%a(2))
      end do
      end do
  
      cden=(-s562+cmz2)*f62
* quqd -- p=p5,q=p62
      quqd=p5(0)*p62(0)-p5(1)*p62(1)-p5(2)*p62(2)-p5(3)*p62(3)
      ccr=zcr(id5)/cden
      ccl=zcl(id5)/cden
* TL0 -- qu=p5,qd=p62,v=0,a=lz562(0)%a,c=lz562(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p5(2)*p62(3)+p62(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p62(0)+p62k0*p5(0)
      lz562(0)%a(1)=ccr*(auxa+ceps_0)
      lz562(0)%a(2)=ccl*(auxa-ceps_0)
      lz562(0)%c(1)=ccr*(p5(2)+ceps_1)
      lz562(0)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p62,v=1,a=lz562(1)%a,c=lz562(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p62(1)+p62k0*p5(1)
      lz562(1)%a(1)=ccr*(auxa+ceps_0)
      lz562(1)%a(2)=ccl*(auxa-ceps_0)
      lz562(1)%c(1)=ccr*(p5(2)+ceps_1)
      lz562(1)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p62,v=2,a=lz562(2)%a,c=lz562(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p5k0*p62(3)+p62k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p62(2)+p62k0*p5(2)
      lz562(2)%a(1)=ccr*(auxa+ceps_0)
      lz562(2)%a(2)=ccl*(auxa-ceps_0)
      lz562(2)%c(1)=ccr*p5k0
      lz562(2)%c(2)=-ccl*p5k0
* TL0 -- qu=p5,qd=p62,v=3,a=lz562(3)%a,c=lz562(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p5k0*p62(2)-p62k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p62(3)+p62k0*p5(3)
      lz562(3)%a(1)=ccr*(auxa+ceps_0)
      lz562(3)%a(2)=ccl*(auxa-ceps_0)
      lz562(3)%c(1)=ccr*ceps_1
      lz562(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz526(&,i1)%e(m),a1=lz562(m)%a,c1=lz562(m)%c,a2=r6_2(i1)%a
* ,b2=r6_2(i1)%b,prq=s62,bef=cz526(&,i1)%e(m)+,aft=
      cz526(1,i1)%e(m)=cz526(1,i1)%e(m)+(lz562(m)%a(1)*r6_2(i1)%
     & a(1)+lz562(m)%c(1)*s62*r6_2(i1)%b(2))
      cz526(2,i1)%e(m)=cz526(2,i1)%e(m)+(lz562(m)%c(2)*s62*r6_2(
     & i1)%b(1)+lz562(m)%a(2)*r6_2(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cz526(i3,i1)%e
      cz526(i3,i1)%ek0=cz526(i3,i1)%e(0)-cz526(i3,i1)%e(1)
      end do
      end do
  
  
      cden= (-s562)
* quqd -- p=p52,q=p6
      quqd=p52(0)*p6(0)-p52(1)*p6(1)-p52(2)*p6(2)-p52(3)*p6(3)
      ccr=fcr(id6)/cden
      ccl=fcl(id6)/cden
* TR0 -- qu=p52,qd=p6,v=0,a=rf652(0)%a,b=rf652(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p52(2)*p6(3)+p6(2)*p52(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p52k0*p6(0)+p6k0*p52(0)
      rf652(0)%a(1)=ccr*(auxa+ceps_0)
      rf652(0)%a(2)=ccl*(auxa-ceps_0)
      rf652(0)%b(1)=-ccl*(p6(2)+ceps_2)
      rf652(0)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p52,qd=p6,v=1,a=rf652(1)%a,b=rf652(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p52k0*p6(1)+p6k0*p52(1)
      rf652(1)%a(1)=ccr*(auxa+ceps_0)
      rf652(1)%a(2)=ccl*(auxa-ceps_0)
      rf652(1)%b(1)=-ccl*(p6(2)+ceps_2)
      rf652(1)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p52,qd=p6,v=2,a=rf652(2)%a,b=rf652(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p52k0*p6(3)+p6k0*p52(3)
      ceps_0=eps_0*cim
      auxa=p52k0*p6(2)+p6k0*p52(2)
      rf652(2)%a(1)=ccr*(auxa+ceps_0)
      rf652(2)%a(2)=ccl*(auxa-ceps_0)
      rf652(2)%b(1)=-ccl*p6k0
      rf652(2)%b(2)=ccr*p6k0
* TR0 -- qu=p52,qd=p6,v=3,a=rf652(3)%a,b=rf652(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p52k0*p6(2)-p6k0*p52(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p52k0*p6(3)+p6k0*p52(3)
      rf652(3)%a(1)=ccr*(auxa+ceps_0)
      rf652(3)%a(2)=ccl*(auxa-ceps_0)
      rf652(3)%b(1)=-ccl*ceps_2
      rf652(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf526(&,i1)%e(m),a1=l5_2(i1)%a,c1=l5_2(i1)%c,a2=rf652(m)%a
* ,b2=rf652(m)%b,prq=s52,bef=,aft=
      cf526(1,i1)%e(m)=(l5_2(i1)%a(1)*rf652(m)%a(1)+l5_2(i1)%c(1
     & )*s52*rf652(m)%b(2))
      cf526(2,i1)%e(m)=(l5_2(i1)%c(2)*s52*rf652(m)%b(1)+l5_2(i1)
     & %a(2)*rf652(m)%a(2))
      end do
      end do
  
      cden=(-s562)*f62
* quqd -- p=p5,q=p62
      quqd=p5(0)*p62(0)-p5(1)*p62(1)-p5(2)*p62(2)-p5(3)*p62(3)
      ccr=fcr(id5)/cden
      ccl=fcl(id5)/cden
* TL0 -- qu=p5,qd=p62,v=0,a=lf562(0)%a,c=lf562(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p5(2)*p62(3)+p62(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p62(0)+p62k0*p5(0)
      lf562(0)%a(1)=ccr*(auxa+ceps_0)
      lf562(0)%a(2)=ccl*(auxa-ceps_0)
      lf562(0)%c(1)=ccr*(p5(2)+ceps_1)
      lf562(0)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p62,v=1,a=lf562(1)%a,c=lf562(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p62(1)+p62k0*p5(1)
      lf562(1)%a(1)=ccr*(auxa+ceps_0)
      lf562(1)%a(2)=ccl*(auxa-ceps_0)
      lf562(1)%c(1)=ccr*(p5(2)+ceps_1)
      lf562(1)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p62,v=2,a=lf562(2)%a,c=lf562(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p5k0*p62(3)+p62k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p62(2)+p62k0*p5(2)
      lf562(2)%a(1)=ccr*(auxa+ceps_0)
      lf562(2)%a(2)=ccl*(auxa-ceps_0)
      lf562(2)%c(1)=ccr*p5k0
      lf562(2)%c(2)=-ccl*p5k0
* TL0 -- qu=p5,qd=p62,v=3,a=lf562(3)%a,c=lf562(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p5k0*p62(2)-p62k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p62(3)+p62k0*p5(3)
      lf562(3)%a(1)=ccr*(auxa+ceps_0)
      lf562(3)%a(2)=ccl*(auxa-ceps_0)
      lf562(3)%c(1)=ccr*ceps_1
      lf562(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf526(&,i1)%e(m),a1=lf562(m)%a,c1=lf562(m)%c,a2=r6_2(i1)%a
* ,b2=r6_2(i1)%b,prq=s62,bef=cf526(&,i1)%e(m)+,aft=
      cf526(1,i1)%e(m)=cf526(1,i1)%e(m)+(lf562(m)%a(1)*r6_2(i1)%
     & a(1)+lf562(m)%c(1)*s62*r6_2(i1)%b(2))
      cf526(2,i1)%e(m)=cf526(2,i1)%e(m)+(lf562(m)%c(2)*s62*r6_2(
     & i1)%b(1)+lf562(m)%a(2)*r6_2(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cf526(i3,i1)%e
      cf526(i3,i1)%ek0=cf526(i3,i1)%e(0)-cf526(i3,i1)%e(1)
      end do
      end do
  
      endif  ! id5 = quark
  
  
      do mu=0,3
         p71(mu)=p7(mu)+p1(mu)
*	 p781(mu)=p71(mu)+p8(mu)                                              
      enddo
* pk0 -- p=p71
      p71k0=p71(0)-p71(1)
* p.q -- p.q=s71,p=p71,q=p71,bef=,aft=
      s71=(p71(0)*p71(0)-p71(1)*p71(1)-p71(2)*p71(2)-p71(3)*p71(
     & 3))
      f71=s71*p71k0
  
* p.q -- p.q=s781,p=p71,q=p8,bef=2.d0*,aft=+s71
      s781=2.d0*(p71(0)*p8(0)-p71(1)*p8(1)-p71(2)*p8(2)-p71(3)*p
     & 8(3))+s71
      do mu=0,3
         p81(mu)=-p1(mu)-p8(mu)
      enddo
* pk0 -- p=p81
      p81k0=p81(0)-p81(1)
* p.q -- p.q=s81,p=p81,q=p81,bef=,aft=
      s81=(p81(0)*p81(0)-p81(1)*p81(1)-p81(2)*p81(2)-p81(3)*p81(
     & 3))
      f81=s81*p81k0
  
      if (ilept(id7).ne.1) then
      quqd=s71/2d0
      ccr=1.d0/(f71)
      ccl=1.d0/(f71)
      do i1=1,2
* TL0 -- qu=p7,qd=p71,v=ce1(i1)%e,a=l7_1(i1)%a,c=l7_1(i1)%c,cr=ccr,cl=cc
* l,nsum=0
      ceps_0=-ce1(i1)%ek0*(p7(2)*p71(3)-p71(2)*p7(3))+p7k0*(ce1(
     & i1)%e(2)*p71(3)-p71(2)*ce1(i1)%e(3))-p71k0*(ce1(i1)%e(2)*
     & p7(3)-p7(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p7k0+p7(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=ce1(i1)%e(0)*p7(0)-ce1(i1)%e(1)*p7(1)-ce1(i1)%e(2)*p7
     & (2)-ce1(i1)%e(3)*p7(3)
      cvqd=ce1(i1)%e(0)*p71(0)-ce1(i1)%e(1)*p71(1)-ce1(i1)%e(2)*
     & p71(2)-ce1(i1)%e(3)*p71(3)
      cauxa=-ce1(i1)%ek0*quqd+p7k0*cvqd+p71k0*cvqu
      cauxc=+ce1(i1)%ek0*p7(2)-p7k0*ce1(i1)%e(2)
      l7_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      l7_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      l7_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s81/2d0
      do i1=1,2
* TR0 -- qu=p81,qd=p8,v=ce1(i1)%e,a=r8_1(i1)%a,b=r8_1(i1)%b,cr=1.d0,cl=1
* .d0,nsum=0
      ceps_0=-ce1(i1)%ek0*(p81(2)*p8(3)-p8(2)*p81(3))+p81k0*(ce1
     & (i1)%e(2)*p8(3)-p8(2)*ce1(i1)%e(3))-p8k0*(ce1(i1)%e(2)*p8
     & 1(3)-p81(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-ce1(i1)%e(3)*p8k0+p8(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p81(0)-ce1(i1)%e(1)*p81(1)-ce1(i1)%e(2)*
     & p81(2)-ce1(i1)%e(3)*p81(3)
      cvqd=ce1(i1)%e(0)*p8(0)-ce1(i1)%e(1)*p8(1)-ce1(i1)%e(2)*p8
     & (2)-ce1(i1)%e(3)*p8(3)
      cauxa=-ce1(i1)%ek0*quqd+p81k0*cvqd+p8k0*cvqu
      cauxb=-ce1(i1)%ek0*p8(2)+p8k0*ce1(i1)%e(2)
      r8_1(i1)%a(1)=(cauxa+ceps_0)
      r8_1(i1)%a(2)=(cauxa-ceps_0)
      r8_1(i1)%b(1)=(cauxb-ceps_2)
      r8_1(i1)%b(2)=(-cauxb-ceps_2)
      end do
  
  
      cden= (-s781+cmz2)
* quqd -- p=p71,q=p8
      quqd=p71(0)*p8(0)-p71(1)*p8(1)-p71(2)*p8(2)-p71(3)*p8(3)
      ccr=zcr(id8)/cden
      ccl=zcl(id8)/cden
* TR0 -- qu=p71,qd=p8,v=0,a=rz871(0)%a,b=rz871(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p71(2)*p8(3)+p8(2)*p71(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p71k0*p8(0)+p8k0*p71(0)
      rz871(0)%a(1)=ccr*(auxa+ceps_0)
      rz871(0)%a(2)=ccl*(auxa-ceps_0)
      rz871(0)%b(1)=-ccl*(p8(2)+ceps_2)
      rz871(0)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p71,qd=p8,v=1,a=rz871(1)%a,b=rz871(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p71k0*p8(1)+p8k0*p71(1)
      rz871(1)%a(1)=ccr*(auxa+ceps_0)
      rz871(1)%a(2)=ccl*(auxa-ceps_0)
      rz871(1)%b(1)=-ccl*(p8(2)+ceps_2)
      rz871(1)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p71,qd=p8,v=2,a=rz871(2)%a,b=rz871(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p71k0*p8(3)+p8k0*p71(3)
      ceps_0=eps_0*cim
      auxa=p71k0*p8(2)+p8k0*p71(2)
      rz871(2)%a(1)=ccr*(auxa+ceps_0)
      rz871(2)%a(2)=ccl*(auxa-ceps_0)
      rz871(2)%b(1)=-ccl*p8k0
      rz871(2)%b(2)=ccr*p8k0
* TR0 -- qu=p71,qd=p8,v=3,a=rz871(3)%a,b=rz871(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p71k0*p8(2)-p8k0*p71(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p71k0*p8(3)+p8k0*p71(3)
      rz871(3)%a(1)=ccr*(auxa+ceps_0)
      rz871(3)%a(2)=ccl*(auxa-ceps_0)
      rz871(3)%b(1)=-ccl*ceps_2
      rz871(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz718(&,i1)%e(m),a1=l7_1(i1)%a,c1=l7_1(i1)%c,a2=rz871(m)%a
* ,b2=rz871(m)%b,prq=s71,bef=,aft=
      cz718(1,i1)%e(m)=(l7_1(i1)%a(1)*rz871(m)%a(1)+l7_1(i1)%c(1
     & )*s71*rz871(m)%b(2))
      cz718(2,i1)%e(m)=(l7_1(i1)%c(2)*s71*rz871(m)%b(1)+l7_1(i1)
     & %a(2)*rz871(m)%a(2))
      end do
      end do
  
      cden=(-s781+cmz2)*f81
* quqd -- p=p7,q=p81
      quqd=p7(0)*p81(0)-p7(1)*p81(1)-p7(2)*p81(2)-p7(3)*p81(3)
      ccr=zcr(id7)/cden
      ccl=zcl(id7)/cden
* TL0 -- qu=p7,qd=p81,v=0,a=lz781(0)%a,c=lz781(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p7(2)*p81(3)+p81(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p81(0)+p81k0*p7(0)
      lz781(0)%a(1)=ccr*(auxa+ceps_0)
      lz781(0)%a(2)=ccl*(auxa-ceps_0)
      lz781(0)%c(1)=ccr*(p7(2)+ceps_1)
      lz781(0)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p81,v=1,a=lz781(1)%a,c=lz781(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p81(1)+p81k0*p7(1)
      lz781(1)%a(1)=ccr*(auxa+ceps_0)
      lz781(1)%a(2)=ccl*(auxa-ceps_0)
      lz781(1)%c(1)=ccr*(p7(2)+ceps_1)
      lz781(1)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p81,v=2,a=lz781(2)%a,c=lz781(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p7k0*p81(3)+p81k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p81(2)+p81k0*p7(2)
      lz781(2)%a(1)=ccr*(auxa+ceps_0)
      lz781(2)%a(2)=ccl*(auxa-ceps_0)
      lz781(2)%c(1)=ccr*p7k0
      lz781(2)%c(2)=-ccl*p7k0
* TL0 -- qu=p7,qd=p81,v=3,a=lz781(3)%a,c=lz781(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p7k0*p81(2)-p81k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p81(3)+p81k0*p7(3)
      lz781(3)%a(1)=ccr*(auxa+ceps_0)
      lz781(3)%a(2)=ccl*(auxa-ceps_0)
      lz781(3)%c(1)=ccr*ceps_1
      lz781(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz718(&,i1)%e(m),a1=lz781(m)%a,c1=lz781(m)%c,a2=r8_1(i1)%a
* ,b2=r8_1(i1)%b,prq=s81,bef=cz718(&,i1)%e(m)+,aft=
      cz718(1,i1)%e(m)=cz718(1,i1)%e(m)+(lz781(m)%a(1)*r8_1(i1)%
     & a(1)+lz781(m)%c(1)*s81*r8_1(i1)%b(2))
      cz718(2,i1)%e(m)=cz718(2,i1)%e(m)+(lz781(m)%c(2)*s81*r8_1(
     & i1)%b(1)+lz781(m)%a(2)*r8_1(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cz718(i3,i1)%e
      cz718(i3,i1)%ek0=cz718(i3,i1)%e(0)-cz718(i3,i1)%e(1)
      end do
      end do
  
  
      cden= (-s781)
* quqd -- p=p71,q=p8
      quqd=p71(0)*p8(0)-p71(1)*p8(1)-p71(2)*p8(2)-p71(3)*p8(3)
      ccr=fcr(id8)/cden
      ccl=fcl(id8)/cden
* TR0 -- qu=p71,qd=p8,v=0,a=rf871(0)%a,b=rf871(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p71(2)*p8(3)+p8(2)*p71(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p71k0*p8(0)+p8k0*p71(0)
      rf871(0)%a(1)=ccr*(auxa+ceps_0)
      rf871(0)%a(2)=ccl*(auxa-ceps_0)
      rf871(0)%b(1)=-ccl*(p8(2)+ceps_2)
      rf871(0)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p71,qd=p8,v=1,a=rf871(1)%a,b=rf871(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p71k0*p8(1)+p8k0*p71(1)
      rf871(1)%a(1)=ccr*(auxa+ceps_0)
      rf871(1)%a(2)=ccl*(auxa-ceps_0)
      rf871(1)%b(1)=-ccl*(p8(2)+ceps_2)
      rf871(1)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p71,qd=p8,v=2,a=rf871(2)%a,b=rf871(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p71k0*p8(3)+p8k0*p71(3)
      ceps_0=eps_0*cim
      auxa=p71k0*p8(2)+p8k0*p71(2)
      rf871(2)%a(1)=ccr*(auxa+ceps_0)
      rf871(2)%a(2)=ccl*(auxa-ceps_0)
      rf871(2)%b(1)=-ccl*p8k0
      rf871(2)%b(2)=ccr*p8k0
* TR0 -- qu=p71,qd=p8,v=3,a=rf871(3)%a,b=rf871(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p71k0*p8(2)-p8k0*p71(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p71k0*p8(3)+p8k0*p71(3)
      rf871(3)%a(1)=ccr*(auxa+ceps_0)
      rf871(3)%a(2)=ccl*(auxa-ceps_0)
      rf871(3)%b(1)=-ccl*ceps_2
      rf871(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf718(&,i1)%e(m),a1=l7_1(i1)%a,c1=l7_1(i1)%c,a2=rf871(m)%a
* ,b2=rf871(m)%b,prq=s71,bef=,aft=
      cf718(1,i1)%e(m)=(l7_1(i1)%a(1)*rf871(m)%a(1)+l7_1(i1)%c(1
     & )*s71*rf871(m)%b(2))
      cf718(2,i1)%e(m)=(l7_1(i1)%c(2)*s71*rf871(m)%b(1)+l7_1(i1)
     & %a(2)*rf871(m)%a(2))
      end do
      end do
  
      cden=(-s781)*f81
* quqd -- p=p7,q=p81
      quqd=p7(0)*p81(0)-p7(1)*p81(1)-p7(2)*p81(2)-p7(3)*p81(3)
      ccr=fcr(id7)/cden
      ccl=fcl(id7)/cden
* TL0 -- qu=p7,qd=p81,v=0,a=lf781(0)%a,c=lf781(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p7(2)*p81(3)+p81(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p81(0)+p81k0*p7(0)
      lf781(0)%a(1)=ccr*(auxa+ceps_0)
      lf781(0)%a(2)=ccl*(auxa-ceps_0)
      lf781(0)%c(1)=ccr*(p7(2)+ceps_1)
      lf781(0)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p81,v=1,a=lf781(1)%a,c=lf781(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p81(1)+p81k0*p7(1)
      lf781(1)%a(1)=ccr*(auxa+ceps_0)
      lf781(1)%a(2)=ccl*(auxa-ceps_0)
      lf781(1)%c(1)=ccr*(p7(2)+ceps_1)
      lf781(1)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p81,v=2,a=lf781(2)%a,c=lf781(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p7k0*p81(3)+p81k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p81(2)+p81k0*p7(2)
      lf781(2)%a(1)=ccr*(auxa+ceps_0)
      lf781(2)%a(2)=ccl*(auxa-ceps_0)
      lf781(2)%c(1)=ccr*p7k0
      lf781(2)%c(2)=-ccl*p7k0
* TL0 -- qu=p7,qd=p81,v=3,a=lf781(3)%a,c=lf781(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p7k0*p81(2)-p81k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p81(3)+p81k0*p7(3)
      lf781(3)%a(1)=ccr*(auxa+ceps_0)
      lf781(3)%a(2)=ccl*(auxa-ceps_0)
      lf781(3)%c(1)=ccr*ceps_1
      lf781(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf718(&,i1)%e(m),a1=lf781(m)%a,c1=lf781(m)%c,a2=r8_1(i1)%a
* ,b2=r8_1(i1)%b,prq=s81,bef=cf718(&,i1)%e(m)+,aft=
      cf718(1,i1)%e(m)=cf718(1,i1)%e(m)+(lf781(m)%a(1)*r8_1(i1)%
     & a(1)+lf781(m)%c(1)*s81*r8_1(i1)%b(2))
      cf718(2,i1)%e(m)=cf718(2,i1)%e(m)+(lf781(m)%c(2)*s81*r8_1(
     & i1)%b(1)+lf781(m)%a(2)*r8_1(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cf718(i3,i1)%e
      cf718(i3,i1)%ek0=cf718(i3,i1)%e(0)-cf718(i3,i1)%e(1)
      end do
      end do
  
      endif  ! id7 = quark
  
      do mu=0,3
         p72(mu)=p7(mu)+p2(mu)
*	 p782(mu)=p72(mu)+p8(mu)                                              
      enddo
* pk0 -- p=p72
      p72k0=p72(0)-p72(1)
* p.q -- p.q=s72,p=p72,q=p72,bef=,aft=
      s72=(p72(0)*p72(0)-p72(1)*p72(1)-p72(2)*p72(2)-p72(3)*p72(
     & 3))
      f72=s72*p72k0
  
* p.q -- p.q=s782,p=p72,q=p8,bef=2.d0*,aft=+s72
      s782=2.d0*(p72(0)*p8(0)-p72(1)*p8(1)-p72(2)*p8(2)-p72(3)*p
     & 8(3))+s72
      do mu=0,3
         p82(mu)=-p2(mu)-p8(mu)
      enddo
* pk0 -- p=p82
      p82k0=p82(0)-p82(1)
* p.q -- p.q=s82,p=p82,q=p82,bef=,aft=
      s82=(p82(0)*p82(0)-p82(1)*p82(1)-p82(2)*p82(2)-p82(3)*p82(
     & 3))
      f82=s82*p82k0
  
      if (ilept(id7).ne.1) then
      quqd=s72/2d0
      ccr=1.d0/(f72)
      ccl=1.d0/(f72)
      do i1=1,2
* TL0 -- qu=p7,qd=p72,v=ce2(i1)%e,a=l7_2(i1)%a,c=l7_2(i1)%c,cr=ccr,cl=cc
* l,nsum=0
      ceps_0=-ce2(i1)%ek0*(p7(2)*p72(3)-p72(2)*p7(3))+p7k0*(ce2(
     & i1)%e(2)*p72(3)-p72(2)*ce2(i1)%e(3))-p72k0*(ce2(i1)%e(2)*
     & p7(3)-p7(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p7k0+p7(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=ce2(i1)%e(0)*p7(0)-ce2(i1)%e(1)*p7(1)-ce2(i1)%e(2)*p7
     & (2)-ce2(i1)%e(3)*p7(3)
      cvqd=ce2(i1)%e(0)*p72(0)-ce2(i1)%e(1)*p72(1)-ce2(i1)%e(2)*
     & p72(2)-ce2(i1)%e(3)*p72(3)
      cauxa=-ce2(i1)%ek0*quqd+p7k0*cvqd+p72k0*cvqu
      cauxc=+ce2(i1)%ek0*p7(2)-p7k0*ce2(i1)%e(2)
      l7_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      l7_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      l7_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s82/2d0
      do i1=1,2
* TR0 -- qu=p82,qd=p8,v=ce2(i1)%e,a=r8_2(i1)%a,b=r8_2(i1)%b,cr=1.d0,cl=1
* .d0,nsum=0
      ceps_0=-ce2(i1)%ek0*(p82(2)*p8(3)-p8(2)*p82(3))+p82k0*(ce2
     & (i1)%e(2)*p8(3)-p8(2)*ce2(i1)%e(3))-p8k0*(ce2(i1)%e(2)*p8
     & 2(3)-p82(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-ce2(i1)%e(3)*p8k0+p8(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p82(0)-ce2(i1)%e(1)*p82(1)-ce2(i1)%e(2)*
     & p82(2)-ce2(i1)%e(3)*p82(3)
      cvqd=ce2(i1)%e(0)*p8(0)-ce2(i1)%e(1)*p8(1)-ce2(i1)%e(2)*p8
     & (2)-ce2(i1)%e(3)*p8(3)
      cauxa=-ce2(i1)%ek0*quqd+p82k0*cvqd+p8k0*cvqu
      cauxb=-ce2(i1)%ek0*p8(2)+p8k0*ce2(i1)%e(2)
      r8_2(i1)%a(1)=(cauxa+ceps_0)
      r8_2(i1)%a(2)=(cauxa-ceps_0)
      r8_2(i1)%b(1)=(cauxb-ceps_2)
      r8_2(i1)%b(2)=(-cauxb-ceps_2)
      end do
  
  
      cden= (-s782+cmz2)
* quqd -- p=p72,q=p8
      quqd=p72(0)*p8(0)-p72(1)*p8(1)-p72(2)*p8(2)-p72(3)*p8(3)
      ccr=zcr(id8)/cden
      ccl=zcl(id8)/cden
* TR0 -- qu=p72,qd=p8,v=0,a=rz872(0)%a,b=rz872(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p72(2)*p8(3)+p8(2)*p72(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p72k0*p8(0)+p8k0*p72(0)
      rz872(0)%a(1)=ccr*(auxa+ceps_0)
      rz872(0)%a(2)=ccl*(auxa-ceps_0)
      rz872(0)%b(1)=-ccl*(p8(2)+ceps_2)
      rz872(0)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p72,qd=p8,v=1,a=rz872(1)%a,b=rz872(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p72k0*p8(1)+p8k0*p72(1)
      rz872(1)%a(1)=ccr*(auxa+ceps_0)
      rz872(1)%a(2)=ccl*(auxa-ceps_0)
      rz872(1)%b(1)=-ccl*(p8(2)+ceps_2)
      rz872(1)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p72,qd=p8,v=2,a=rz872(2)%a,b=rz872(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p72k0*p8(3)+p8k0*p72(3)
      ceps_0=eps_0*cim
      auxa=p72k0*p8(2)+p8k0*p72(2)
      rz872(2)%a(1)=ccr*(auxa+ceps_0)
      rz872(2)%a(2)=ccl*(auxa-ceps_0)
      rz872(2)%b(1)=-ccl*p8k0
      rz872(2)%b(2)=ccr*p8k0
* TR0 -- qu=p72,qd=p8,v=3,a=rz872(3)%a,b=rz872(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p72k0*p8(2)-p8k0*p72(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p72k0*p8(3)+p8k0*p72(3)
      rz872(3)%a(1)=ccr*(auxa+ceps_0)
      rz872(3)%a(2)=ccl*(auxa-ceps_0)
      rz872(3)%b(1)=-ccl*ceps_2
      rz872(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz728(&,i1)%e(m),a1=l7_2(i1)%a,c1=l7_2(i1)%c,a2=rz872(m)%a
* ,b2=rz872(m)%b,prq=s72,bef=,aft=
      cz728(1,i1)%e(m)=(l7_2(i1)%a(1)*rz872(m)%a(1)+l7_2(i1)%c(1
     & )*s72*rz872(m)%b(2))
      cz728(2,i1)%e(m)=(l7_2(i1)%c(2)*s72*rz872(m)%b(1)+l7_2(i1)
     & %a(2)*rz872(m)%a(2))
      end do
      end do
  
      cden=(-s782+cmz2)*f82
* quqd -- p=p7,q=p82
      quqd=p7(0)*p82(0)-p7(1)*p82(1)-p7(2)*p82(2)-p7(3)*p82(3)
      ccr=zcr(id7)/cden
      ccl=zcl(id7)/cden
* TL0 -- qu=p7,qd=p82,v=0,a=lz782(0)%a,c=lz782(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p7(2)*p82(3)+p82(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p82(0)+p82k0*p7(0)
      lz782(0)%a(1)=ccr*(auxa+ceps_0)
      lz782(0)%a(2)=ccl*(auxa-ceps_0)
      lz782(0)%c(1)=ccr*(p7(2)+ceps_1)
      lz782(0)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p82,v=1,a=lz782(1)%a,c=lz782(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p82(1)+p82k0*p7(1)
      lz782(1)%a(1)=ccr*(auxa+ceps_0)
      lz782(1)%a(2)=ccl*(auxa-ceps_0)
      lz782(1)%c(1)=ccr*(p7(2)+ceps_1)
      lz782(1)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p82,v=2,a=lz782(2)%a,c=lz782(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p7k0*p82(3)+p82k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p82(2)+p82k0*p7(2)
      lz782(2)%a(1)=ccr*(auxa+ceps_0)
      lz782(2)%a(2)=ccl*(auxa-ceps_0)
      lz782(2)%c(1)=ccr*p7k0
      lz782(2)%c(2)=-ccl*p7k0
* TL0 -- qu=p7,qd=p82,v=3,a=lz782(3)%a,c=lz782(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p7k0*p82(2)-p82k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p82(3)+p82k0*p7(3)
      lz782(3)%a(1)=ccr*(auxa+ceps_0)
      lz782(3)%a(2)=ccl*(auxa-ceps_0)
      lz782(3)%c(1)=ccr*ceps_1
      lz782(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz728(&,i1)%e(m),a1=lz782(m)%a,c1=lz782(m)%c,a2=r8_2(i1)%a
* ,b2=r8_2(i1)%b,prq=s82,bef=cz728(&,i1)%e(m)+,aft=
      cz728(1,i1)%e(m)=cz728(1,i1)%e(m)+(lz782(m)%a(1)*r8_2(i1)%
     & a(1)+lz782(m)%c(1)*s82*r8_2(i1)%b(2))
      cz728(2,i1)%e(m)=cz728(2,i1)%e(m)+(lz782(m)%c(2)*s82*r8_2(
     & i1)%b(1)+lz782(m)%a(2)*r8_2(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cz728(i3,i1)%e
      cz728(i3,i1)%ek0=cz728(i3,i1)%e(0)-cz728(i3,i1)%e(1)
      end do
      end do
  
  
      cden= (-s782)
* quqd -- p=p72,q=p8
      quqd=p72(0)*p8(0)-p72(1)*p8(1)-p72(2)*p8(2)-p72(3)*p8(3)
      ccr=fcr(id8)/cden
      ccl=fcl(id8)/cden
* TR0 -- qu=p72,qd=p8,v=0,a=rf872(0)%a,b=rf872(0)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p72(2)*p8(3)+p8(2)*p72(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p72k0*p8(0)+p8k0*p72(0)
      rf872(0)%a(1)=ccr*(auxa+ceps_0)
      rf872(0)%a(2)=ccl*(auxa-ceps_0)
      rf872(0)%b(1)=-ccl*(p8(2)+ceps_2)
      rf872(0)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p72,qd=p8,v=1,a=rf872(1)%a,b=rf872(1)%b,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p72k0*p8(1)+p8k0*p72(1)
      rf872(1)%a(1)=ccr*(auxa+ceps_0)
      rf872(1)%a(2)=ccl*(auxa-ceps_0)
      rf872(1)%b(1)=-ccl*(p8(2)+ceps_2)
      rf872(1)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p72,qd=p8,v=2,a=rf872(2)%a,b=rf872(2)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=-p72k0*p8(3)+p8k0*p72(3)
      ceps_0=eps_0*cim
      auxa=p72k0*p8(2)+p8k0*p72(2)
      rf872(2)%a(1)=ccr*(auxa+ceps_0)
      rf872(2)%a(2)=ccl*(auxa-ceps_0)
      rf872(2)%b(1)=-ccl*p8k0
      rf872(2)%b(2)=ccr*p8k0
* TR0 -- qu=p72,qd=p8,v=3,a=rf872(3)%a,b=rf872(3)%b,cr=ccr,cl=ccl,nsum=0
      eps_0=p72k0*p8(2)-p8k0*p72(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p72k0*p8(3)+p8k0*p72(3)
      rf872(3)%a(1)=ccr*(auxa+ceps_0)
      rf872(3)%a(2)=ccl*(auxa-ceps_0)
      rf872(3)%b(1)=-ccl*ceps_2
      rf872(3)%b(2)=-ccr*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf728(&,i1)%e(m),a1=l7_2(i1)%a,c1=l7_2(i1)%c,a2=rf872(m)%a
* ,b2=rf872(m)%b,prq=s72,bef=,aft=
      cf728(1,i1)%e(m)=(l7_2(i1)%a(1)*rf872(m)%a(1)+l7_2(i1)%c(1
     & )*s72*rf872(m)%b(2))
      cf728(2,i1)%e(m)=(l7_2(i1)%c(2)*s72*rf872(m)%b(1)+l7_2(i1)
     & %a(2)*rf872(m)%a(2))
      end do
      end do
  
      cden=(-s782)*f82
* quqd -- p=p7,q=p82
      quqd=p7(0)*p82(0)-p7(1)*p82(1)-p7(2)*p82(2)-p7(3)*p82(3)
      ccr=fcr(id7)/cden
      ccl=fcl(id7)/cden
* TL0 -- qu=p7,qd=p82,v=0,a=lf782(0)%a,c=lf782(0)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p7(2)*p82(3)+p82(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p82(0)+p82k0*p7(0)
      lf782(0)%a(1)=ccr*(auxa+ceps_0)
      lf782(0)%a(2)=ccl*(auxa-ceps_0)
      lf782(0)%c(1)=ccr*(p7(2)+ceps_1)
      lf782(0)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p82,v=1,a=lf782(1)%a,c=lf782(1)%c,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p82(1)+p82k0*p7(1)
      lf782(1)%a(1)=ccr*(auxa+ceps_0)
      lf782(1)%a(2)=ccl*(auxa-ceps_0)
      lf782(1)%c(1)=ccr*(p7(2)+ceps_1)
      lf782(1)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p82,v=2,a=lf782(2)%a,c=lf782(2)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=-p7k0*p82(3)+p82k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p82(2)+p82k0*p7(2)
      lf782(2)%a(1)=ccr*(auxa+ceps_0)
      lf782(2)%a(2)=ccl*(auxa-ceps_0)
      lf782(2)%c(1)=ccr*p7k0
      lf782(2)%c(2)=-ccl*p7k0
* TL0 -- qu=p7,qd=p82,v=3,a=lf782(3)%a,c=lf782(3)%c,cr=ccr,cl=ccl,nsum=0
      eps_0=p7k0*p82(2)-p82k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p82(3)+p82k0*p7(3)
      lf782(3)%a(1)=ccr*(auxa+ceps_0)
      lf782(3)%a(2)=ccl*(auxa-ceps_0)
      lf782(3)%c(1)=ccr*ceps_1
      lf782(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cf728(&,i1)%e(m),a1=lf782(m)%a,c1=lf782(m)%c,a2=r8_2(i1)%a
* ,b2=r8_2(i1)%b,prq=s82,bef=cf728(&,i1)%e(m)+,aft=
      cf728(1,i1)%e(m)=cf728(1,i1)%e(m)+(lf782(m)%a(1)*r8_2(i1)%
     & a(1)+lf782(m)%c(1)*s82*r8_2(i1)%b(2))
      cf728(2,i1)%e(m)=cf728(2,i1)%e(m)+(lf782(m)%c(2)*s82*r8_2(
     & i1)%b(1)+lf782(m)%a(2)*r8_2(i1)%a(2))
      end do
      end do
  
      do i3=1,2
      do i1=1,2
* pk0 -- p=cf728(i3,i1)%e
      cf728(i3,i1)%ek0=cf728(i3,i1)%e(0)-cf728(i3,i1)%e(1)
      end do
      end do
  
      endif  ! id7 = quark
  
*     -two-gluons fork                                                  
  
      do m=0,3
       p312(m)=p31(m) + p2(m)
      enddo
* pk0 -- p=p312
      p312k0=p312(0)-p312(1)
* p.q -- p.q=s312,p=p312,q=p312,bef=,aft=
      s312=(p312(0)*p312(0)-p312(1)*p312(1)-p312(2)*p312(2)-p312
     & (3)*p312(3))
      f312=s312*p312k0
  
      do m=0,3
       p412(m)=p41(m) - p2(m)
      enddo
* pk0 -- p=p412
      p412k0=p412(0)-p412(1)
* p.q -- p.q=s412,p=p412,q=p412,bef=,aft=
      s412=(p412(0)*p412(0)-p412(1)*p412(1)-p412(2)*p412(2)-p412
     & (3)*p412(3))
      f412=s412*p412k0
  
* p.q -- p.q=s3412,p=p312,q=p4,bef=2.d0*,aft=+s312
      s3412=2.d0*(p312(0)*p4(0)-p312(1)*p4(1)-p312(2)*p4(2)-p312
     & (3)*p4(3))+s312
  
      if (ilept(id3).ne.1) then
  
* triple vertex in left/right line                                      
* quqd -- p=p412,q=p4
      quqd=p412(0)*p4(0)-p412(1)*p4(1)-p412(2)*p4(2)-p412(3)*p4(
     & 3)
      ccr=1.d0/(s12)
      ccl=1.d0/(s12)
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p412,qd=p4,v=ctrip12(i1,i2)%e,a=r4_gg(i1,i2)%a,b=r4_gg(i1,i2
* )%b,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p412(2)*p4(3)-p4(2)*p412(3))+p
     & 412k0*(ctrip12(i1,i2)%e(2)*p4(3)-p4(2)*ctrip12(i1,i2)%e(3
     & ))-p4k0*(ctrip12(i1,i2)%e(2)*p412(3)-p412(2)*ctrip12(i1,i
     & 2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-ctrip12(i1,i2)%e(3)*p4k0+p4(3)*ctrip12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=ctrip12(i1,i2)%e(0)*p412(0)-ctrip12(i1,i2)%e(1)*p412(
     & 1)-ctrip12(i1,i2)%e(2)*p412(2)-ctrip12(i1,i2)%e(3)*p412(3
     & )
      cvqd=ctrip12(i1,i2)%e(0)*p4(0)-ctrip12(i1,i2)%e(1)*p4(1)-c
     & trip12(i1,i2)%e(2)*p4(2)-ctrip12(i1,i2)%e(3)*p4(3)
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p412k0*cvqd+p4k0*cvqu
      cauxb=-ctrip12(i1,i2)%ek0*p4(2)+p4k0*ctrip12(i1,i2)%e(2)
      r4_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      r4_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      r4_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      r4_gg(i1,i2)%b(2)=ccr*(-cauxb-ceps_2)
      end do
      end do
  
* quqd -- p=p3,q=p312
      quqd=p3(0)*p312(0)-p3(1)*p312(1)-p3(2)*p312(2)-p3(3)*p312(
     & 3)
      ccr=1.d0/(f312*s12)
      ccl=1.d0/(f312*s12)
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p3,qd=p312,v=ctrip12(i1,i2)%e,a=l3_gg(i1,i2)%a,c=l3_gg(i1,i2
* )%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p3(2)*p312(3)-p312(2)*p3(3))+p
     & 3k0*(ctrip12(i1,i2)%e(2)*p312(3)-p312(2)*ctrip12(i1,i2)%e
     & (3))-p312k0*(ctrip12(i1,i2)%e(2)*p3(3)-p3(2)*ctrip12(i1,i
     & 2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ctrip12(i1,i2)%e(3)*p3k0+p3(3)*ctrip12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=ctrip12(i1,i2)%e(0)*p3(0)-ctrip12(i1,i2)%e(1)*p3(1)-c
     & trip12(i1,i2)%e(2)*p3(2)-ctrip12(i1,i2)%e(3)*p3(3)
      cvqd=ctrip12(i1,i2)%e(0)*p312(0)-ctrip12(i1,i2)%e(1)*p312(
     & 1)-ctrip12(i1,i2)%e(2)*p312(2)-ctrip12(i1,i2)%e(3)*p312(3
     & )
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p3k0*cvqd+p312k0*cvqu
      cauxc=+ctrip12(i1,i2)%ek0*p3(2)-p3k0*ctrip12(i1,i2)%e(2)
      l3_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l3_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l3_gg(i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l3_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      cden=(-s3412+cmz2)
* quqd -- p=p312,q=p4
      quqd=p312(0)*p4(0)-p312(1)*p4(1)-p312(2)*p4(2)-p312(3)*p4(
     & 3)
      ccr=zcr(id4)/cden
      ccl=zcl(id4)/cden
* TR0 -- qu=p312,qd=p4,v=0,a=rz4312(0)%a,b=rz4312(0)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rz4312(0)%a(1)=ccr*(auxa+ceps_0)
      rz4312(0)%a(2)=ccl*(auxa-ceps_0)
      rz4312(0)%b(1)=-ccl*(p4(2)+ceps_2)
      rz4312(0)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=1,a=rz4312(1)%a,b=rz4312(1)%b,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rz4312(1)%a(1)=ccr*(auxa+ceps_0)
      rz4312(1)%a(2)=ccl*(auxa-ceps_0)
      rz4312(1)%b(1)=-ccl*(p4(2)+ceps_2)
      rz4312(1)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=2,a=rz4312(2)%a,b=rz4312(2)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rz4312(2)%a(1)=ccr*(auxa+ceps_0)
      rz4312(2)%a(2)=ccl*(auxa-ceps_0)
      rz4312(2)%b(1)=-ccl*p4k0
      rz4312(2)%b(2)=ccr*p4k0
* TR0 -- qu=p312,qd=p4,v=3,a=rz4312(3)%a,b=rz4312(3)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p312k0*p4(3)+p4k0*p312(3)
      rz4312(3)%a(1)=ccr*(auxa+ceps_0)
      rz4312(3)%a(2)=ccl*(auxa-ceps_0)
      rz4312(3)%b(1)=-ccl*ceps_2
      rz4312(3)%b(2)=-ccr*ceps_2
  
      cden=(-s3412+cmz2)*f412
* quqd -- p=p3,q=p412
      quqd=p3(0)*p412(0)-p3(1)*p412(1)-p3(2)*p412(2)-p3(3)*p412(
     & 3)
      ccr=zcr(id3)/cden
      ccl=zcl(id3)/cden
* TL0 -- qu=p3,qd=p412,v=0,a=lz3412(0)%a,c=lz3412(0)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lz3412(0)%a(1)=ccr*(auxa+ceps_0)
      lz3412(0)%a(2)=ccl*(auxa-ceps_0)
      lz3412(0)%c(1)=ccr*(p3(2)+ceps_1)
      lz3412(0)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=1,a=lz3412(1)%a,c=lz3412(1)%c,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lz3412(1)%a(1)=ccr*(auxa+ceps_0)
      lz3412(1)%a(2)=ccl*(auxa-ceps_0)
      lz3412(1)%c(1)=ccr*(p3(2)+ceps_1)
      lz3412(1)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=2,a=lz3412(2)%a,c=lz3412(2)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lz3412(2)%a(1)=ccr*(auxa+ceps_0)
      lz3412(2)%a(2)=ccl*(auxa-ceps_0)
      lz3412(2)%c(1)=ccr*p3k0
      lz3412(2)%c(2)=-ccl*p3k0
* TL0 -- qu=p3,qd=p412,v=3,a=lz3412(3)%a,c=lz3412(3)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p412(3)+p412k0*p3(3)
      lz3412(3)%a(1)=ccr*(auxa+ceps_0)
      lz3412(3)%a(2)=ccl*(auxa-ceps_0)
      lz3412(3)%c(1)=ccr*ceps_1
      lz3412(3)%c(2)=ccl*ceps_1
      cden=(-s3412)
* quqd -- p=p312,q=p4
      quqd=p312(0)*p4(0)-p312(1)*p4(1)-p312(2)*p4(2)-p312(3)*p4(
     & 3)
      ccr=fcr(id4)/cden
      ccl=fcl(id4)/cden
* TR0 -- qu=p312,qd=p4,v=0,a=rf4312(0)%a,b=rf4312(0)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rf4312(0)%a(1)=ccr*(auxa+ceps_0)
      rf4312(0)%a(2)=ccl*(auxa-ceps_0)
      rf4312(0)%b(1)=-ccl*(p4(2)+ceps_2)
      rf4312(0)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=1,a=rf4312(1)%a,b=rf4312(1)%b,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rf4312(1)%a(1)=ccr*(auxa+ceps_0)
      rf4312(1)%a(2)=ccl*(auxa-ceps_0)
      rf4312(1)%b(1)=-ccl*(p4(2)+ceps_2)
      rf4312(1)%b(2)=ccr*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=2,a=rf4312(2)%a,b=rf4312(2)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rf4312(2)%a(1)=ccr*(auxa+ceps_0)
      rf4312(2)%a(2)=ccl*(auxa-ceps_0)
      rf4312(2)%b(1)=-ccl*p4k0
      rf4312(2)%b(2)=ccr*p4k0
* TR0 -- qu=p312,qd=p4,v=3,a=rf4312(3)%a,b=rf4312(3)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p312k0*p4(3)+p4k0*p312(3)
      rf4312(3)%a(1)=ccr*(auxa+ceps_0)
      rf4312(3)%a(2)=ccl*(auxa-ceps_0)
      rf4312(3)%b(1)=-ccl*ceps_2
      rf4312(3)%b(2)=-ccr*ceps_2
  
      cden=(-s3412)*f412
* quqd -- p=p3,q=p412
      quqd=p3(0)*p412(0)-p3(1)*p412(1)-p3(2)*p412(2)-p3(3)*p412(
     & 3)
      ccr=fcr(id3)/cden
      ccl=fcl(id3)/cden
* TL0 -- qu=p3,qd=p412,v=0,a=lf3412(0)%a,c=lf3412(0)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lf3412(0)%a(1)=ccr*(auxa+ceps_0)
      lf3412(0)%a(2)=ccl*(auxa-ceps_0)
      lf3412(0)%c(1)=ccr*(p3(2)+ceps_1)
      lf3412(0)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=1,a=lf3412(1)%a,c=lf3412(1)%c,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lf3412(1)%a(1)=ccr*(auxa+ceps_0)
      lf3412(1)%a(2)=ccl*(auxa-ceps_0)
      lf3412(1)%c(1)=ccr*(p3(2)+ceps_1)
      lf3412(1)%c(2)=ccl*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=2,a=lf3412(2)%a,c=lf3412(2)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lf3412(2)%a(1)=ccr*(auxa+ceps_0)
      lf3412(2)%a(2)=ccl*(auxa-ceps_0)
      lf3412(2)%c(1)=ccr*p3k0
      lf3412(2)%c(2)=-ccl*p3k0
* TL0 -- qu=p3,qd=p412,v=3,a=lf3412(3)%a,c=lf3412(3)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p412(3)+p412k0*p3(3)
      lf3412(3)%a(1)=ccr*(auxa+ceps_0)
      lf3412(3)%a(2)=ccl*(auxa-ceps_0)
      lf3412(3)%c(1)=ccr*ceps_1
      lf3412(3)%c(2)=ccl*ceps_1
  
* quqd -- p=p31,q=p312
      quqd=p31(0)*p312(0)-p31(1)*p312(1)-p31(2)*p312(2)-p31(3)*p
     & 312(3)
      ccr=1.d0/(f312)
      ccl=1.d0/(f312)
      do i2=1,2
* T0 -- qu=p31,qd=p312,v=ce2(i2)%e,a=u31_2(i2)%a,b=u31_2(i2)%b,c=u31_2(i
* 2)%c,d=u31_2(i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i2)%ek0*(p31(2)*p312(3)-p312(2)*p31(3))+p31k0*
     & (ce2(i2)%e(2)*p312(3)-p312(2)*ce2(i2)%e(3))-p312k0*(ce2(i
     & 2)%e(2)*p31(3)-p31(2)*ce2(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i2)%e(3)*p31k0+p31(3)*ce2(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i2)%e(3)*p312k0+p312(3)*ce2(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i2)%e(0)*p31(0)-ce2(i2)%e(1)*p31(1)-ce2(i2)%e(2)*
     & p31(2)-ce2(i2)%e(3)*p31(3)
      cvqd=ce2(i2)%e(0)*p312(0)-ce2(i2)%e(1)*p312(1)-ce2(i2)%e(2
     & )*p312(2)-ce2(i2)%e(3)*p312(3)
      cauxa=-ce2(i2)%ek0*quqd+p31k0*cvqd+p312k0*cvqu
      cauxb=-ce2(i2)%ek0*p312(2)+p312k0*ce2(i2)%e(2)
      cauxc=+ce2(i2)%ek0*p31(2)-p31k0*ce2(i2)%e(2)
      u31_2(i2)%a(1)=ccr*(cauxa+ceps_0)
      u31_2(i2)%a(2)=ccl*(cauxa-ceps_0)
      u31_2(i2)%b(1)=ccl*(cauxb-ceps_2)
      u31_2(i2)%b(2)=ccr*(-cauxb-ceps_2)
      u31_2(i2)%c(1)=ccr*(cauxc+ceps_1)
      u31_2(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u31_2(i2)%d(1)=ccl*ce2(i2)%ek0
      u31_2(i2)%d(2)=ccr*ce2(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=l3_12(i1,i2)%a,cc=l3_12(i1,i2)%c,a1=l3_1(i1)%a,c1=l3_1(i1)%
* c,a2=u31_2(i2)%a,b2=u31_2(i2)%b,c2=u31_2(i2)%c,d2=u31_2(i2)%d,prq=s31,
* nsum=0
      l3_12(i1,i2)%a(1)=l3_1(i1)%a(1)*u31_2(i2)%a(1)+l3_1(i1)%c(
     & 1)*s31*u31_2(i2)%b(2)
      l3_12(i1,i2)%c(1)=l3_1(i1)%a(1)*u31_2(i2)%c(1)+l3_1(i1)%c(
     & 1)*s31*u31_2(i2)%d(2)
      l3_12(i1,i2)%c(2)=l3_1(i1)%c(2)*s31*u31_2(i2)%d(1)+l3_1(i1
     & )%a(2)*u31_2(i2)%c(2)
      l3_12(i1,i2)%a(2)=l3_1(i1)%c(2)*s31*u31_2(i2)%b(1)+l3_1(i1
     & )%a(2)*u31_2(i2)%a(2)
      end do
      end do
  
* quqd -- p=p412,q=p42
      quqd=p412(0)*p42(0)-p412(1)*p42(1)-p412(2)*p42(2)-p412(3)*
     & p42(3)
      ccr=1.d0/(f42)
      ccl=1.d0/(f42)
      do i1=1,2
* T0 -- qu=p412,qd=p42,v=ce1(i1)%e,a=u42_1(i1)%a,b=u42_1(i1)%b,c=u42_1(i
* 1)%c,d=u42_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p412(2)*p42(3)-p42(2)*p412(3))+p412k0
     & *(ce1(i1)%e(2)*p42(3)-p42(2)*ce1(i1)%e(3))-p42k0*(ce1(i1)
     & %e(2)*p412(3)-p412(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p412k0+p412(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p42k0+p42(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p412(0)-ce1(i1)%e(1)*p412(1)-ce1(i1)%e(2
     & )*p412(2)-ce1(i1)%e(3)*p412(3)
      cvqd=ce1(i1)%e(0)*p42(0)-ce1(i1)%e(1)*p42(1)-ce1(i1)%e(2)*
     & p42(2)-ce1(i1)%e(3)*p42(3)
      cauxa=-ce1(i1)%ek0*quqd+p412k0*cvqd+p42k0*cvqu
      cauxb=-ce1(i1)%ek0*p42(2)+p42k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p412(2)-p412k0*ce1(i1)%e(2)
      u42_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u42_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u42_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u42_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u42_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u42_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u42_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u42_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0 -- aa=r4_12(i1,i2)%a,bb=r4_12(i1,i2)%b,a1=u42_1(i1)%a,b1=u42_1(i1
* )%b,c1=u42_1(i1)%c,d1=u42_1(i1)%d,a2=r4_2(i2)%a,b2=r4_2(i2)%b,prq=s42,
* nsum=0
      r4_12(i1,i2)%a(1)=u42_1(i1)%a(1)*r4_2(i2)%a(1)+u42_1(i1)%c
     & (1)*s42*r4_2(i2)%b(2)
      r4_12(i1,i2)%b(1)=u42_1(i1)%d(1)*s42*r4_2(i2)%b(1)+u42_1(i
     & 1)%b(1)*r4_2(i2)%a(2)
      r4_12(i1,i2)%b(2)=u42_1(i1)%b(2)*r4_2(i2)%a(1)+u42_1(i1)%d
     & (2)*s42*r4_2(i2)%b(2)
      r4_12(i1,i2)%a(2)=u42_1(i1)%c(2)*s42*r4_2(i2)%b(1)+u42_1(i
     & 1)%a(2)*r4_2(i2)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l3_12(i1,i2)%a(1)=l3_12(i1,i2)%a(1) +
     &   l3_gg(i1,i2)%a(1)
       l3_12(i1,i2)%a(2)=l3_12(i1,i2)%a(2) +
     &   l3_gg(i1,i2)%a(2)
       l3_12(i1,i2)%c(1)=l3_12(i1,i2)%c(1) +
     &   l3_gg(i1,i2)%c(1)
       l3_12(i1,i2)%c(2)=l3_12(i1,i2)%c(2) +
     &   l3_gg(i1,i2)%c(2)
  
       r4_12(i1,i2)%a(1)=r4_12(i1,i2)%a(1) +
     &   r4_gg(i1,i2)%a(1)
       r4_12(i1,i2)%a(2)=r4_12(i1,i2)%a(2) +
     &   r4_gg(i1,i2)%a(2)
       r4_12(i1,i2)%b(1)=r4_12(i1,i2)%b(1) +
     &   r4_gg(i1,i2)%b(1)
       r4_12(i1,i2)%b(2)=r4_12(i1,i2)%b(2) +
     &   r4_gg(i1,i2)%b(2)
  
      enddo
      enddo
  
  
      cden=(-s3412+cmz2)*f42
* quqd -- p=p31,q=p42
      quqd=p31(0)*p42(0)-p31(1)*p42(1)-p31(2)*p42(2)-p31(3)*p42(
     & 3)
      ccr=zcr(id3)/cden
      ccl=zcl(id3)/cden
* T0 -- qu=p31,qd=p42,v=0,a=uz31(0)%a,b=uz31(0)%b,c=uz31(0)%c,d=uz31(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p31(2)*p42(3)+p42(2)*p31(3)
      ceps_0=eps_0*cim
      ceps_1=p31(3)*cim
      ceps_2=p42(3)*cim
      auxa=-quqd+p31k0*p42(0)+p42k0*p31(0)
      uz31(0)%a(1)=ccr*(auxa+ceps_0)
      uz31(0)%a(2)=ccl*(auxa-ceps_0)
      uz31(0)%b(1)=-ccl*(p42(2)+ceps_2)
      uz31(0)%b(2)=ccr*(p42(2)-ceps_2)
      uz31(0)%c(1)=ccr*(p31(2)+ceps_1)
      uz31(0)%c(2)=ccl*(-p31(2)+ceps_1)
      uz31(0)%d(1)=ccl
      uz31(0)%d(2)=ccr
* T0 -- qu=p31,qd=p42,v=1,a=uz31(1)%a,b=uz31(1)%b,c=uz31(1)%c,d=uz31(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p31k0*p42(1)+p42k0*p31(1)
      uz31(1)%a(1)=ccr*(auxa+ceps_0)
      uz31(1)%a(2)=ccl*(auxa-ceps_0)
      uz31(1)%b(1)=-ccl*(p42(2)+ceps_2)
      uz31(1)%b(2)=ccr*(p42(2)-ceps_2)
      uz31(1)%c(1)=ccr*(p31(2)+ceps_1)
      uz31(1)%c(2)=ccl*(-p31(2)+ceps_1)
      uz31(1)%d(1)=ccl
      uz31(1)%d(2)=ccr
* T0 -- qu=p31,qd=p42,v=2,a=uz31(2)%a,b=uz31(2)%b,c=uz31(2)%c,d=uz31(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p31k0*p42(3)+p42k0*p31(3)
      ceps_0=eps_0*cim
      auxa=p31k0*p42(2)+p42k0*p31(2)
      uz31(2)%a(1)=ccr*(auxa+ceps_0)
      uz31(2)%a(2)=ccl*(auxa-ceps_0)
      uz31(2)%b(1)=-ccl*p42k0
      uz31(2)%b(2)=ccr*p42k0
      uz31(2)%c(1)=ccr*p31k0
      uz31(2)%c(2)=-ccl*p31k0
* T0 -- qu=p31,qd=p42,v=3,a=uz31(3)%a,b=uz31(3)%b,c=uz31(3)%c,d=uz31(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p31k0*p42(2)-p42k0*p31(2)
      ceps_0=eps_0*cim
      ceps_1=p31k0*cim
      ceps_2=p42k0*cim
      auxa=+p31k0*p42(3)+p42k0*p31(3)
      uz31(3)%a(1)=ccr*(auxa+ceps_0)
      uz31(3)%a(2)=ccl*(auxa-ceps_0)
      uz31(3)%b(1)=-ccl*ceps_2
      uz31(3)%b(2)=-ccr*ceps_2
      uz31(3)%c(1)=ccr*ceps_1
      uz31(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l3_1(i1)%a,c1=l3_
* 1(i1)%c,a2=uz31(mu)%a,b2=uz31(mu)%b,c2=uz31(mu)%c,d2=uz31(mu)%d,prq=s3
* 1,nsum=0
      laux_imu(i1,mu)%a(1)=l3_1(i1)%a(1)*uz31(mu)%a(1)+l3_1(i1)%
     & c(1)*s31*uz31(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l3_1(i1)%a(1)*uz31(mu)%c(1)+l3_1(i1)%
     & c(1)*s31*uz31(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l3_1(i1)%c(2)*s31*uz31(mu)%d(1)+l3_1(
     & i1)%a(2)*uz31(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l3_1(i1)%c(2)*s31*uz31(mu)%b(1)+l3_1(
     & i1)%a(2)*uz31(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz3124(&,i1,i2)%e(mu),a1=laux_imu(i1,mu)%a,c1=laux_imu(i1,
* mu)%c,a2=r4_2(i2)%a,b2=r4_2(i2)%b,prq=s42,bef=,aft=
      cz3124(1,i1,i2)%e(mu)=(laux_imu(i1,mu)%a(1)*r4_2(i2)%a(1)+
     & laux_imu(i1,mu)%c(1)*s42*r4_2(i2)%b(2))
      cz3124(2,i1,i2)%e(mu)=(laux_imu(i1,mu)%c(2)*s42*r4_2(i2)%b
     & (1)+laux_imu(i1,mu)%a(2)*r4_2(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz3124(&,i1,i2)%e(mu),a1=l3_12(i1,i2)%a,c1=l3_12(i1,i2)%c,
* a2=rz4312(mu)%a,b2=rz4312(mu)%b,prq=s312,bef=cz3124(&,i1,i2)%e(mu)+,af
* t=
      cz3124(1,i1,i2)%e(mu)=cz3124(1,i1,i2)%e(mu)+(l3_12(i1,i2)%
     & a(1)*rz4312(mu)%a(1)+l3_12(i1,i2)%c(1)*s312*rz4312(mu)%b(
     & 2))
      cz3124(2,i1,i2)%e(mu)=cz3124(2,i1,i2)%e(mu)+(l3_12(i1,i2)%
     & c(2)*s312*rz4312(mu)%b(1)+l3_12(i1,i2)%a(2)*rz4312(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz3124(&,i1,i2)%e(mu),a1=lz3412(mu)%a,c1=lz3412(mu)%c,a2=r
* 4_12(i1,i2)%a,b2=r4_12(i1,i2)%b,prq=s412,bef=cz3124(&,i1,i2)%e(mu)+,af
* t=
      cz3124(1,i1,i2)%e(mu)=cz3124(1,i1,i2)%e(mu)+(lz3412(mu)%a(
     & 1)*r4_12(i1,i2)%a(1)+lz3412(mu)%c(1)*s412*r4_12(i1,i2)%b(
     & 2))
      cz3124(2,i1,i2)%e(mu)=cz3124(2,i1,i2)%e(mu)+(lz3412(mu)%c(
     & 2)*s412*r4_12(i1,i2)%b(1)+lz3412(mu)%a(2)*r4_12(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cz3124(i3,i1,i2)%e
      cz3124(i3,i1,i2)%ek0=cz3124(i3,i1,i2)%e(0)-cz3124(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
      cden=(-s3412)*f42
* quqd -- p=p31,q=p42
      quqd=p31(0)*p42(0)-p31(1)*p42(1)-p31(2)*p42(2)-p31(3)*p42(
     & 3)
      ccr=fcr(id3)/cden
      ccl=fcl(id3)/cden
* T0 -- qu=p31,qd=p42,v=0,a=uf31(0)%a,b=uf31(0)%b,c=uf31(0)%c,d=uf31(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p31(2)*p42(3)+p42(2)*p31(3)
      ceps_0=eps_0*cim
      ceps_1=p31(3)*cim
      ceps_2=p42(3)*cim
      auxa=-quqd+p31k0*p42(0)+p42k0*p31(0)
      uf31(0)%a(1)=ccr*(auxa+ceps_0)
      uf31(0)%a(2)=ccl*(auxa-ceps_0)
      uf31(0)%b(1)=-ccl*(p42(2)+ceps_2)
      uf31(0)%b(2)=ccr*(p42(2)-ceps_2)
      uf31(0)%c(1)=ccr*(p31(2)+ceps_1)
      uf31(0)%c(2)=ccl*(-p31(2)+ceps_1)
      uf31(0)%d(1)=ccl
      uf31(0)%d(2)=ccr
* T0 -- qu=p31,qd=p42,v=1,a=uf31(1)%a,b=uf31(1)%b,c=uf31(1)%c,d=uf31(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p31k0*p42(1)+p42k0*p31(1)
      uf31(1)%a(1)=ccr*(auxa+ceps_0)
      uf31(1)%a(2)=ccl*(auxa-ceps_0)
      uf31(1)%b(1)=-ccl*(p42(2)+ceps_2)
      uf31(1)%b(2)=ccr*(p42(2)-ceps_2)
      uf31(1)%c(1)=ccr*(p31(2)+ceps_1)
      uf31(1)%c(2)=ccl*(-p31(2)+ceps_1)
      uf31(1)%d(1)=ccl
      uf31(1)%d(2)=ccr
* T0 -- qu=p31,qd=p42,v=2,a=uf31(2)%a,b=uf31(2)%b,c=uf31(2)%c,d=uf31(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p31k0*p42(3)+p42k0*p31(3)
      ceps_0=eps_0*cim
      auxa=p31k0*p42(2)+p42k0*p31(2)
      uf31(2)%a(1)=ccr*(auxa+ceps_0)
      uf31(2)%a(2)=ccl*(auxa-ceps_0)
      uf31(2)%b(1)=-ccl*p42k0
      uf31(2)%b(2)=ccr*p42k0
      uf31(2)%c(1)=ccr*p31k0
      uf31(2)%c(2)=-ccl*p31k0
* T0 -- qu=p31,qd=p42,v=3,a=uf31(3)%a,b=uf31(3)%b,c=uf31(3)%c,d=uf31(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p31k0*p42(2)-p42k0*p31(2)
      ceps_0=eps_0*cim
      ceps_1=p31k0*cim
      ceps_2=p42k0*cim
      auxa=+p31k0*p42(3)+p42k0*p31(3)
      uf31(3)%a(1)=ccr*(auxa+ceps_0)
      uf31(3)%a(2)=ccl*(auxa-ceps_0)
      uf31(3)%b(1)=-ccl*ceps_2
      uf31(3)%b(2)=-ccr*ceps_2
      uf31(3)%c(1)=ccr*ceps_1
      uf31(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l3_1(i1)%a,c1=l3_
* 1(i1)%c,a2=uf31(mu)%a,b2=uf31(mu)%b,c2=uf31(mu)%c,d2=uf31(mu)%d,prq=s3
* 1,nsum=0
      laux_imu(i1,mu)%a(1)=l3_1(i1)%a(1)*uf31(mu)%a(1)+l3_1(i1)%
     & c(1)*s31*uf31(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l3_1(i1)%a(1)*uf31(mu)%c(1)+l3_1(i1)%
     & c(1)*s31*uf31(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l3_1(i1)%c(2)*s31*uf31(mu)%d(1)+l3_1(
     & i1)%a(2)*uf31(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l3_1(i1)%c(2)*s31*uf31(mu)%b(1)+l3_1(
     & i1)%a(2)*uf31(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf3124(&,i1,i2)%e(mu),a1=laux_imu(i1,mu)%a,c1=laux_imu(i1,
* mu)%c,a2=r4_2(i2)%a,b2=r4_2(i2)%b,prq=s42,bef=,aft=
      cf3124(1,i1,i2)%e(mu)=(laux_imu(i1,mu)%a(1)*r4_2(i2)%a(1)+
     & laux_imu(i1,mu)%c(1)*s42*r4_2(i2)%b(2))
      cf3124(2,i1,i2)%e(mu)=(laux_imu(i1,mu)%c(2)*s42*r4_2(i2)%b
     & (1)+laux_imu(i1,mu)%a(2)*r4_2(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf3124(&,i1,i2)%e(mu),a1=l3_12(i1,i2)%a,c1=l3_12(i1,i2)%c,
* a2=rf4312(mu)%a,b2=rf4312(mu)%b,prq=s312,bef=cf3124(&,i1,i2)%e(mu)+,af
* t=
      cf3124(1,i1,i2)%e(mu)=cf3124(1,i1,i2)%e(mu)+(l3_12(i1,i2)%
     & a(1)*rf4312(mu)%a(1)+l3_12(i1,i2)%c(1)*s312*rf4312(mu)%b(
     & 2))
      cf3124(2,i1,i2)%e(mu)=cf3124(2,i1,i2)%e(mu)+(l3_12(i1,i2)%
     & c(2)*s312*rf4312(mu)%b(1)+l3_12(i1,i2)%a(2)*rf4312(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf3124(&,i1,i2)%e(mu),a1=lf3412(mu)%a,c1=lf3412(mu)%c,a2=r
* 4_12(i1,i2)%a,b2=r4_12(i1,i2)%b,prq=s412,bef=cf3124(&,i1,i2)%e(mu)+,af
* t=
      cf3124(1,i1,i2)%e(mu)=cf3124(1,i1,i2)%e(mu)+(lf3412(mu)%a(
     & 1)*r4_12(i1,i2)%a(1)+lf3412(mu)%c(1)*s412*r4_12(i1,i2)%b(
     & 2))
      cf3124(2,i1,i2)%e(mu)=cf3124(2,i1,i2)%e(mu)+(lf3412(mu)%c(
     & 2)*s412*r4_12(i1,i2)%b(1)+lf3412(mu)%a(2)*r4_12(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cf3124(i3,i1,i2)%e
      cf3124(i3,i1,i2)%ek0=cf3124(i3,i1,i2)%e(0)-cf3124(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
* quqd -- p=p32,q=p312
      quqd=p32(0)*p312(0)-p32(1)*p312(1)-p32(2)*p312(2)-p32(3)*p
     & 312(3)
      ccr=1.d0/(f312)
      ccl=1.d0/(f312)
      do i2=1,2
* T0 -- qu=p32,qd=p312,v=ce1(i2)%e,a=u32_1(i2)%a,b=u32_1(i2)%b,c=u32_1(i
* 2)%c,d=u32_1(i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i2)%ek0*(p32(2)*p312(3)-p312(2)*p32(3))+p32k0*
     & (ce1(i2)%e(2)*p312(3)-p312(2)*ce1(i2)%e(3))-p312k0*(ce1(i
     & 2)%e(2)*p32(3)-p32(2)*ce1(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i2)%e(3)*p32k0+p32(3)*ce1(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i2)%e(3)*p312k0+p312(3)*ce1(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i2)%e(0)*p32(0)-ce1(i2)%e(1)*p32(1)-ce1(i2)%e(2)*
     & p32(2)-ce1(i2)%e(3)*p32(3)
      cvqd=ce1(i2)%e(0)*p312(0)-ce1(i2)%e(1)*p312(1)-ce1(i2)%e(2
     & )*p312(2)-ce1(i2)%e(3)*p312(3)
      cauxa=-ce1(i2)%ek0*quqd+p32k0*cvqd+p312k0*cvqu
      cauxb=-ce1(i2)%ek0*p312(2)+p312k0*ce1(i2)%e(2)
      cauxc=+ce1(i2)%ek0*p32(2)-p32k0*ce1(i2)%e(2)
      u32_1(i2)%a(1)=ccr*(cauxa+ceps_0)
      u32_1(i2)%a(2)=ccl*(cauxa-ceps_0)
      u32_1(i2)%b(1)=ccl*(cauxb-ceps_2)
      u32_1(i2)%b(2)=ccr*(-cauxb-ceps_2)
      u32_1(i2)%c(1)=ccr*(cauxc+ceps_1)
      u32_1(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u32_1(i2)%d(1)=ccl*ce1(i2)%ek0
      u32_1(i2)%d(2)=ccr*ce1(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=l3_21(i1,i2)%a,cc=l3_21(i1,i2)%c,a1=l3_2(i2)%a,c1=l3_2(i2)%
* c,a2=u32_1(i1)%a,b2=u32_1(i1)%b,c2=u32_1(i1)%c,d2=u32_1(i1)%d,prq=s32,
* nsum=0
      l3_21(i1,i2)%a(1)=l3_2(i2)%a(1)*u32_1(i1)%a(1)+l3_2(i2)%c(
     & 1)*s32*u32_1(i1)%b(2)
      l3_21(i1,i2)%c(1)=l3_2(i2)%a(1)*u32_1(i1)%c(1)+l3_2(i2)%c(
     & 1)*s32*u32_1(i1)%d(2)
      l3_21(i1,i2)%c(2)=l3_2(i2)%c(2)*s32*u32_1(i1)%d(1)+l3_2(i2
     & )%a(2)*u32_1(i1)%c(2)
      l3_21(i1,i2)%a(2)=l3_2(i2)%c(2)*s32*u32_1(i1)%b(1)+l3_2(i2
     & )%a(2)*u32_1(i1)%a(2)
      end do
      end do
  
* quqd -- p=p412,q=p41
      quqd=p412(0)*p41(0)-p412(1)*p41(1)-p412(2)*p41(2)-p412(3)*
     & p41(3)
      ccr=1.d0/(f41)
      ccl=1.d0/(f41)
      do i1=1,2
* T0 -- qu=p412,qd=p41,v=ce2(i1)%e,a=u41_2(i1)%a,b=u41_2(i1)%b,c=u41_2(i
* 1)%c,d=u41_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p412(2)*p41(3)-p41(2)*p412(3))+p412k0
     & *(ce2(i1)%e(2)*p41(3)-p41(2)*ce2(i1)%e(3))-p41k0*(ce2(i1)
     & %e(2)*p412(3)-p412(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p412k0+p412(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p41k0+p41(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p412(0)-ce2(i1)%e(1)*p412(1)-ce2(i1)%e(2
     & )*p412(2)-ce2(i1)%e(3)*p412(3)
      cvqd=ce2(i1)%e(0)*p41(0)-ce2(i1)%e(1)*p41(1)-ce2(i1)%e(2)*
     & p41(2)-ce2(i1)%e(3)*p41(3)
      cauxa=-ce2(i1)%ek0*quqd+p412k0*cvqd+p41k0*cvqu
      cauxb=-ce2(i1)%ek0*p41(2)+p41k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p412(2)-p412k0*ce2(i1)%e(2)
      u41_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u41_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u41_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u41_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u41_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u41_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u41_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u41_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0 -- aa=r4_21(i1,i2)%a,bb=r4_21(i1,i2)%b,a1=u41_2(i2)%a,b1=u41_2(i2
* )%b,c1=u41_2(i2)%c,d1=u41_2(i2)%d,a2=r4_1(i1)%a,b2=r4_1(i1)%b,prq=s41,
* nsum=0
      r4_21(i1,i2)%a(1)=u41_2(i2)%a(1)*r4_1(i1)%a(1)+u41_2(i2)%c
     & (1)*s41*r4_1(i1)%b(2)
      r4_21(i1,i2)%b(1)=u41_2(i2)%d(1)*s41*r4_1(i1)%b(1)+u41_2(i
     & 2)%b(1)*r4_1(i1)%a(2)
      r4_21(i1,i2)%b(2)=u41_2(i2)%b(2)*r4_1(i1)%a(1)+u41_2(i2)%d
     & (2)*s41*r4_1(i1)%b(2)
      r4_21(i1,i2)%a(2)=u41_2(i2)%c(2)*s41*r4_1(i1)%b(1)+u41_2(i
     & 2)%a(2)*r4_1(i1)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l3_21(i1,i2)%a(1)=l3_21(i1,i2)%a(1) -
     &   l3_gg(i1,i2)%a(1)
       l3_21(i1,i2)%a(2)=l3_21(i1,i2)%a(2) -
     &   l3_gg(i1,i2)%a(2)
       l3_21(i1,i2)%c(1)=l3_21(i1,i2)%c(1) -
     &   l3_gg(i1,i2)%c(1)
       l3_21(i1,i2)%c(2)=l3_21(i1,i2)%c(2) -
     &   l3_gg(i1,i2)%c(2)
  
       r4_21(i1,i2)%a(1)=r4_21(i1,i2)%a(1) -
     &   r4_gg(i1,i2)%a(1)
       r4_21(i1,i2)%a(2)=r4_21(i1,i2)%a(2) -
     &   r4_gg(i1,i2)%a(2)
       r4_21(i1,i2)%b(1)=r4_21(i1,i2)%b(1) -
     &   r4_gg(i1,i2)%b(1)
       r4_21(i1,i2)%b(2)=r4_21(i1,i2)%b(2) -
     &   r4_gg(i1,i2)%b(2)
  
      enddo
      enddo
  
  
      cden=(-s3412+cmz2)*f41
* quqd -- p=p32,q=p41
      quqd=p32(0)*p41(0)-p32(1)*p41(1)-p32(2)*p41(2)-p32(3)*p41(
     & 3)
      ccr=zcr(id3)/cden
      ccl=zcl(id3)/cden
* T0 -- qu=p32,qd=p41,v=0,a=uz32(0)%a,b=uz32(0)%b,c=uz32(0)%c,d=uz32(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p32(2)*p41(3)+p41(2)*p32(3)
      ceps_0=eps_0*cim
      ceps_1=p32(3)*cim
      ceps_2=p41(3)*cim
      auxa=-quqd+p32k0*p41(0)+p41k0*p32(0)
      uz32(0)%a(1)=ccr*(auxa+ceps_0)
      uz32(0)%a(2)=ccl*(auxa-ceps_0)
      uz32(0)%b(1)=-ccl*(p41(2)+ceps_2)
      uz32(0)%b(2)=ccr*(p41(2)-ceps_2)
      uz32(0)%c(1)=ccr*(p32(2)+ceps_1)
      uz32(0)%c(2)=ccl*(-p32(2)+ceps_1)
      uz32(0)%d(1)=ccl
      uz32(0)%d(2)=ccr
* T0 -- qu=p32,qd=p41,v=1,a=uz32(1)%a,b=uz32(1)%b,c=uz32(1)%c,d=uz32(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p32k0*p41(1)+p41k0*p32(1)
      uz32(1)%a(1)=ccr*(auxa+ceps_0)
      uz32(1)%a(2)=ccl*(auxa-ceps_0)
      uz32(1)%b(1)=-ccl*(p41(2)+ceps_2)
      uz32(1)%b(2)=ccr*(p41(2)-ceps_2)
      uz32(1)%c(1)=ccr*(p32(2)+ceps_1)
      uz32(1)%c(2)=ccl*(-p32(2)+ceps_1)
      uz32(1)%d(1)=ccl
      uz32(1)%d(2)=ccr
* T0 -- qu=p32,qd=p41,v=2,a=uz32(2)%a,b=uz32(2)%b,c=uz32(2)%c,d=uz32(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p32k0*p41(3)+p41k0*p32(3)
      ceps_0=eps_0*cim
      auxa=p32k0*p41(2)+p41k0*p32(2)
      uz32(2)%a(1)=ccr*(auxa+ceps_0)
      uz32(2)%a(2)=ccl*(auxa-ceps_0)
      uz32(2)%b(1)=-ccl*p41k0
      uz32(2)%b(2)=ccr*p41k0
      uz32(2)%c(1)=ccr*p32k0
      uz32(2)%c(2)=-ccl*p32k0
* T0 -- qu=p32,qd=p41,v=3,a=uz32(3)%a,b=uz32(3)%b,c=uz32(3)%c,d=uz32(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p32k0*p41(2)-p41k0*p32(2)
      ceps_0=eps_0*cim
      ceps_1=p32k0*cim
      ceps_2=p41k0*cim
      auxa=+p32k0*p41(3)+p41k0*p32(3)
      uz32(3)%a(1)=ccr*(auxa+ceps_0)
      uz32(3)%a(2)=ccl*(auxa-ceps_0)
      uz32(3)%b(1)=-ccl*ceps_2
      uz32(3)%b(2)=-ccr*ceps_2
      uz32(3)%c(1)=ccr*ceps_1
      uz32(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l3_2(i1)%a,c1=l3_
* 2(i1)%c,a2=uz32(mu)%a,b2=uz32(mu)%b,c2=uz32(mu)%c,d2=uz32(mu)%d,prq=s3
* 2,nsum=0
      laux_imu(i1,mu)%a(1)=l3_2(i1)%a(1)*uz32(mu)%a(1)+l3_2(i1)%
     & c(1)*s32*uz32(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l3_2(i1)%a(1)*uz32(mu)%c(1)+l3_2(i1)%
     & c(1)*s32*uz32(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l3_2(i1)%c(2)*s32*uz32(mu)%d(1)+l3_2(
     & i1)%a(2)*uz32(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l3_2(i1)%c(2)*s32*uz32(mu)%b(1)+l3_2(
     & i1)%a(2)*uz32(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz3214(&,i1,i2)%e(mu),a1=laux_imu(i2,mu)%a,c1=laux_imu(i2,
* mu)%c,a2=r4_1(i1)%a,b2=r4_1(i1)%b,prq=s41,bef=,aft=
      cz3214(1,i1,i2)%e(mu)=(laux_imu(i2,mu)%a(1)*r4_1(i1)%a(1)+
     & laux_imu(i2,mu)%c(1)*s41*r4_1(i1)%b(2))
      cz3214(2,i1,i2)%e(mu)=(laux_imu(i2,mu)%c(2)*s41*r4_1(i1)%b
     & (1)+laux_imu(i2,mu)%a(2)*r4_1(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz3214(&,i1,i2)%e(mu),a1=l3_21(i1,i2)%a,c1=l3_21(i1,i2)%c,
* a2=rz4312(mu)%a,b2=rz4312(mu)%b,prq=s312,bef=cz3214(&,i1,i2)%e(mu)+,af
* t=
      cz3214(1,i1,i2)%e(mu)=cz3214(1,i1,i2)%e(mu)+(l3_21(i1,i2)%
     & a(1)*rz4312(mu)%a(1)+l3_21(i1,i2)%c(1)*s312*rz4312(mu)%b(
     & 2))
      cz3214(2,i1,i2)%e(mu)=cz3214(2,i1,i2)%e(mu)+(l3_21(i1,i2)%
     & c(2)*s312*rz4312(mu)%b(1)+l3_21(i1,i2)%a(2)*rz4312(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz3214(&,i1,i2)%e(mu),a1=lz3412(mu)%a,c1=lz3412(mu)%c,a2=r
* 4_21(i1,i2)%a,b2=r4_21(i1,i2)%b,prq=s412,bef=cz3214(&,i1,i2)%e(mu)+,af
* t=
      cz3214(1,i1,i2)%e(mu)=cz3214(1,i1,i2)%e(mu)+(lz3412(mu)%a(
     & 1)*r4_21(i1,i2)%a(1)+lz3412(mu)%c(1)*s412*r4_21(i1,i2)%b(
     & 2))
      cz3214(2,i1,i2)%e(mu)=cz3214(2,i1,i2)%e(mu)+(lz3412(mu)%c(
     & 2)*s412*r4_21(i1,i2)%b(1)+lz3412(mu)%a(2)*r4_21(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cz3214(i3,i1,i2)%e
      cz3214(i3,i1,i2)%ek0=cz3214(i3,i1,i2)%e(0)-cz3214(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
      cden=(-s3412)*f41
* quqd -- p=p32,q=p41
      quqd=p32(0)*p41(0)-p32(1)*p41(1)-p32(2)*p41(2)-p32(3)*p41(
     & 3)
      ccr=fcr(id3)/cden
      ccl=fcl(id3)/cden
* T0 -- qu=p32,qd=p41,v=0,a=uf32(0)%a,b=uf32(0)%b,c=uf32(0)%c,d=uf32(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p32(2)*p41(3)+p41(2)*p32(3)
      ceps_0=eps_0*cim
      ceps_1=p32(3)*cim
      ceps_2=p41(3)*cim
      auxa=-quqd+p32k0*p41(0)+p41k0*p32(0)
      uf32(0)%a(1)=ccr*(auxa+ceps_0)
      uf32(0)%a(2)=ccl*(auxa-ceps_0)
      uf32(0)%b(1)=-ccl*(p41(2)+ceps_2)
      uf32(0)%b(2)=ccr*(p41(2)-ceps_2)
      uf32(0)%c(1)=ccr*(p32(2)+ceps_1)
      uf32(0)%c(2)=ccl*(-p32(2)+ceps_1)
      uf32(0)%d(1)=ccl
      uf32(0)%d(2)=ccr
* T0 -- qu=p32,qd=p41,v=1,a=uf32(1)%a,b=uf32(1)%b,c=uf32(1)%c,d=uf32(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p32k0*p41(1)+p41k0*p32(1)
      uf32(1)%a(1)=ccr*(auxa+ceps_0)
      uf32(1)%a(2)=ccl*(auxa-ceps_0)
      uf32(1)%b(1)=-ccl*(p41(2)+ceps_2)
      uf32(1)%b(2)=ccr*(p41(2)-ceps_2)
      uf32(1)%c(1)=ccr*(p32(2)+ceps_1)
      uf32(1)%c(2)=ccl*(-p32(2)+ceps_1)
      uf32(1)%d(1)=ccl
      uf32(1)%d(2)=ccr
* T0 -- qu=p32,qd=p41,v=2,a=uf32(2)%a,b=uf32(2)%b,c=uf32(2)%c,d=uf32(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p32k0*p41(3)+p41k0*p32(3)
      ceps_0=eps_0*cim
      auxa=p32k0*p41(2)+p41k0*p32(2)
      uf32(2)%a(1)=ccr*(auxa+ceps_0)
      uf32(2)%a(2)=ccl*(auxa-ceps_0)
      uf32(2)%b(1)=-ccl*p41k0
      uf32(2)%b(2)=ccr*p41k0
      uf32(2)%c(1)=ccr*p32k0
      uf32(2)%c(2)=-ccl*p32k0
* T0 -- qu=p32,qd=p41,v=3,a=uf32(3)%a,b=uf32(3)%b,c=uf32(3)%c,d=uf32(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p32k0*p41(2)-p41k0*p32(2)
      ceps_0=eps_0*cim
      ceps_1=p32k0*cim
      ceps_2=p41k0*cim
      auxa=+p32k0*p41(3)+p41k0*p32(3)
      uf32(3)%a(1)=ccr*(auxa+ceps_0)
      uf32(3)%a(2)=ccl*(auxa-ceps_0)
      uf32(3)%b(1)=-ccl*ceps_2
      uf32(3)%b(2)=-ccr*ceps_2
      uf32(3)%c(1)=ccr*ceps_1
      uf32(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l3_2(i1)%a,c1=l3_
* 2(i1)%c,a2=uf32(mu)%a,b2=uf32(mu)%b,c2=uf32(mu)%c,d2=uf32(mu)%d,prq=s3
* 2,nsum=0
      laux_imu(i1,mu)%a(1)=l3_2(i1)%a(1)*uf32(mu)%a(1)+l3_2(i1)%
     & c(1)*s32*uf32(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l3_2(i1)%a(1)*uf32(mu)%c(1)+l3_2(i1)%
     & c(1)*s32*uf32(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l3_2(i1)%c(2)*s32*uf32(mu)%d(1)+l3_2(
     & i1)%a(2)*uf32(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l3_2(i1)%c(2)*s32*uf32(mu)%b(1)+l3_2(
     & i1)%a(2)*uf32(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf3214(&,i1,i2)%e(mu),a1=laux_imu(i2,mu)%a,c1=laux_imu(i2,
* mu)%c,a2=r4_1(i1)%a,b2=r4_1(i1)%b,prq=s41,bef=,aft=
      cf3214(1,i1,i2)%e(mu)=(laux_imu(i2,mu)%a(1)*r4_1(i1)%a(1)+
     & laux_imu(i2,mu)%c(1)*s41*r4_1(i1)%b(2))
      cf3214(2,i1,i2)%e(mu)=(laux_imu(i2,mu)%c(2)*s41*r4_1(i1)%b
     & (1)+laux_imu(i2,mu)%a(2)*r4_1(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf3214(&,i1,i2)%e(mu),a1=l3_21(i1,i2)%a,c1=l3_21(i1,i2)%c,
* a2=rf4312(mu)%a,b2=rf4312(mu)%b,prq=s312,bef=cf3214(&,i1,i2)%e(mu)+,af
* t=
      cf3214(1,i1,i2)%e(mu)=cf3214(1,i1,i2)%e(mu)+(l3_21(i1,i2)%
     & a(1)*rf4312(mu)%a(1)+l3_21(i1,i2)%c(1)*s312*rf4312(mu)%b(
     & 2))
      cf3214(2,i1,i2)%e(mu)=cf3214(2,i1,i2)%e(mu)+(l3_21(i1,i2)%
     & c(2)*s312*rf4312(mu)%b(1)+l3_21(i1,i2)%a(2)*rf4312(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf3214(&,i1,i2)%e(mu),a1=lf3412(mu)%a,c1=lf3412(mu)%c,a2=r
* 4_21(i1,i2)%a,b2=r4_21(i1,i2)%b,prq=s412,bef=cf3214(&,i1,i2)%e(mu)+,af
* t=
      cf3214(1,i1,i2)%e(mu)=cf3214(1,i1,i2)%e(mu)+(lf3412(mu)%a(
     & 1)*r4_21(i1,i2)%a(1)+lf3412(mu)%c(1)*s412*r4_21(i1,i2)%b(
     & 2))
      cf3214(2,i1,i2)%e(mu)=cf3214(2,i1,i2)%e(mu)+(lf3412(mu)%c(
     & 2)*s412*r4_21(i1,i2)%b(1)+lf3412(mu)%a(2)*r4_21(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cf3214(i3,i1,i2)%e
      cf3214(i3,i1,i2)%ek0=cf3214(i3,i1,i2)%e(0)-cf3214(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
      endif !id3 = quark
  
      do m=0,3
       p512(m)=p51(m) + p2(m)
      enddo
* pk0 -- p=p512
      p512k0=p512(0)-p512(1)
* p.q -- p.q=s512,p=p512,q=p512,bef=,aft=
      s512=(p512(0)*p512(0)-p512(1)*p512(1)-p512(2)*p512(2)-p512
     & (3)*p512(3))
      f512=s512*p512k0
  
      do m=0,3
       p612(m)=p61(m) - p2(m)
      enddo
* pk0 -- p=p612
      p612k0=p612(0)-p612(1)
* p.q -- p.q=s612,p=p612,q=p612,bef=,aft=
      s612=(p612(0)*p612(0)-p612(1)*p612(1)-p612(2)*p612(2)-p612
     & (3)*p612(3))
      f612=s612*p612k0
  
* p.q -- p.q=s5612,p=p512,q=p6,bef=2.d0*,aft=+s512
      s5612=2.d0*(p512(0)*p6(0)-p512(1)*p6(1)-p512(2)*p6(2)-p512
     & (3)*p6(3))+s512
  
      if (ilept(id5).ne.1) then
  
* triple vertex in left/right line                                      
* quqd -- p=p612,q=p6
      quqd=p612(0)*p6(0)-p612(1)*p6(1)-p612(2)*p6(2)-p612(3)*p6(
     & 3)
      ccr=1.d0/(s12)
      ccl=1.d0/(s12)
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p612,qd=p6,v=ctrip12(i1,i2)%e,a=r6_gg(i1,i2)%a,b=r6_gg(i1,i2
* )%b,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p612(2)*p6(3)-p6(2)*p612(3))+p
     & 612k0*(ctrip12(i1,i2)%e(2)*p6(3)-p6(2)*ctrip12(i1,i2)%e(3
     & ))-p6k0*(ctrip12(i1,i2)%e(2)*p612(3)-p612(2)*ctrip12(i1,i
     & 2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-ctrip12(i1,i2)%e(3)*p6k0+p6(3)*ctrip12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=ctrip12(i1,i2)%e(0)*p612(0)-ctrip12(i1,i2)%e(1)*p612(
     & 1)-ctrip12(i1,i2)%e(2)*p612(2)-ctrip12(i1,i2)%e(3)*p612(3
     & )
      cvqd=ctrip12(i1,i2)%e(0)*p6(0)-ctrip12(i1,i2)%e(1)*p6(1)-c
     & trip12(i1,i2)%e(2)*p6(2)-ctrip12(i1,i2)%e(3)*p6(3)
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p612k0*cvqd+p6k0*cvqu
      cauxb=-ctrip12(i1,i2)%ek0*p6(2)+p6k0*ctrip12(i1,i2)%e(2)
      r6_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      r6_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      r6_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      r6_gg(i1,i2)%b(2)=ccr*(-cauxb-ceps_2)
      end do
      end do
  
* quqd -- p=p5,q=p512
      quqd=p5(0)*p512(0)-p5(1)*p512(1)-p5(2)*p512(2)-p5(3)*p512(
     & 3)
      ccr=1.d0/(f512*s12)
      ccl=1.d0/(f512*s12)
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p5,qd=p512,v=ctrip12(i1,i2)%e,a=l5_gg(i1,i2)%a,c=l5_gg(i1,i2
* )%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p5(2)*p512(3)-p512(2)*p5(3))+p
     & 5k0*(ctrip12(i1,i2)%e(2)*p512(3)-p512(2)*ctrip12(i1,i2)%e
     & (3))-p512k0*(ctrip12(i1,i2)%e(2)*p5(3)-p5(2)*ctrip12(i1,i
     & 2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ctrip12(i1,i2)%e(3)*p5k0+p5(3)*ctrip12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=ctrip12(i1,i2)%e(0)*p5(0)-ctrip12(i1,i2)%e(1)*p5(1)-c
     & trip12(i1,i2)%e(2)*p5(2)-ctrip12(i1,i2)%e(3)*p5(3)
      cvqd=ctrip12(i1,i2)%e(0)*p512(0)-ctrip12(i1,i2)%e(1)*p512(
     & 1)-ctrip12(i1,i2)%e(2)*p512(2)-ctrip12(i1,i2)%e(3)*p512(3
     & )
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p5k0*cvqd+p512k0*cvqu
      cauxc=+ctrip12(i1,i2)%ek0*p5(2)-p5k0*ctrip12(i1,i2)%e(2)
      l5_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l5_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_gg(i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l5_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      cden=(-s5612+cmz2)
* quqd -- p=p512,q=p6
      quqd=p512(0)*p6(0)-p512(1)*p6(1)-p512(2)*p6(2)-p512(3)*p6(
     & 3)
      ccr=zcr(id6)/cden
      ccl=zcl(id6)/cden
* TR0 -- qu=p512,qd=p6,v=0,a=rz6512(0)%a,b=rz6512(0)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p512(2)*p6(3)+p6(2)*p512(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p512k0*p6(0)+p6k0*p512(0)
      rz6512(0)%a(1)=ccr*(auxa+ceps_0)
      rz6512(0)%a(2)=ccl*(auxa-ceps_0)
      rz6512(0)%b(1)=-ccl*(p6(2)+ceps_2)
      rz6512(0)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p512,qd=p6,v=1,a=rz6512(1)%a,b=rz6512(1)%b,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p512k0*p6(1)+p6k0*p512(1)
      rz6512(1)%a(1)=ccr*(auxa+ceps_0)
      rz6512(1)%a(2)=ccl*(auxa-ceps_0)
      rz6512(1)%b(1)=-ccl*(p6(2)+ceps_2)
      rz6512(1)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p512,qd=p6,v=2,a=rz6512(2)%a,b=rz6512(2)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p512k0*p6(3)+p6k0*p512(3)
      ceps_0=eps_0*cim
      auxa=p512k0*p6(2)+p6k0*p512(2)
      rz6512(2)%a(1)=ccr*(auxa+ceps_0)
      rz6512(2)%a(2)=ccl*(auxa-ceps_0)
      rz6512(2)%b(1)=-ccl*p6k0
      rz6512(2)%b(2)=ccr*p6k0
* TR0 -- qu=p512,qd=p6,v=3,a=rz6512(3)%a,b=rz6512(3)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p512k0*p6(2)-p6k0*p512(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p512k0*p6(3)+p6k0*p512(3)
      rz6512(3)%a(1)=ccr*(auxa+ceps_0)
      rz6512(3)%a(2)=ccl*(auxa-ceps_0)
      rz6512(3)%b(1)=-ccl*ceps_2
      rz6512(3)%b(2)=-ccr*ceps_2
  
      cden=(-s5612+cmz2)*f612
* quqd -- p=p5,q=p612
      quqd=p5(0)*p612(0)-p5(1)*p612(1)-p5(2)*p612(2)-p5(3)*p612(
     & 3)
      ccr=zcr(id5)/cden
      ccl=zcl(id5)/cden
* TL0 -- qu=p5,qd=p612,v=0,a=lz5612(0)%a,c=lz5612(0)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p5(2)*p612(3)+p612(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p612(0)+p612k0*p5(0)
      lz5612(0)%a(1)=ccr*(auxa+ceps_0)
      lz5612(0)%a(2)=ccl*(auxa-ceps_0)
      lz5612(0)%c(1)=ccr*(p5(2)+ceps_1)
      lz5612(0)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p612,v=1,a=lz5612(1)%a,c=lz5612(1)%c,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p5k0*p612(1)+p612k0*p5(1)
      lz5612(1)%a(1)=ccr*(auxa+ceps_0)
      lz5612(1)%a(2)=ccl*(auxa-ceps_0)
      lz5612(1)%c(1)=ccr*(p5(2)+ceps_1)
      lz5612(1)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p612,v=2,a=lz5612(2)%a,c=lz5612(2)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p5k0*p612(3)+p612k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p612(2)+p612k0*p5(2)
      lz5612(2)%a(1)=ccr*(auxa+ceps_0)
      lz5612(2)%a(2)=ccl*(auxa-ceps_0)
      lz5612(2)%c(1)=ccr*p5k0
      lz5612(2)%c(2)=-ccl*p5k0
* TL0 -- qu=p5,qd=p612,v=3,a=lz5612(3)%a,c=lz5612(3)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p5k0*p612(2)-p612k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p612(3)+p612k0*p5(3)
      lz5612(3)%a(1)=ccr*(auxa+ceps_0)
      lz5612(3)%a(2)=ccl*(auxa-ceps_0)
      lz5612(3)%c(1)=ccr*ceps_1
      lz5612(3)%c(2)=ccl*ceps_1
      cden=(-s5612)
* quqd -- p=p512,q=p6
      quqd=p512(0)*p6(0)-p512(1)*p6(1)-p512(2)*p6(2)-p512(3)*p6(
     & 3)
      ccr=fcr(id6)/cden
      ccl=fcl(id6)/cden
* TR0 -- qu=p512,qd=p6,v=0,a=rf6512(0)%a,b=rf6512(0)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p512(2)*p6(3)+p6(2)*p512(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p512k0*p6(0)+p6k0*p512(0)
      rf6512(0)%a(1)=ccr*(auxa+ceps_0)
      rf6512(0)%a(2)=ccl*(auxa-ceps_0)
      rf6512(0)%b(1)=-ccl*(p6(2)+ceps_2)
      rf6512(0)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p512,qd=p6,v=1,a=rf6512(1)%a,b=rf6512(1)%b,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p512k0*p6(1)+p6k0*p512(1)
      rf6512(1)%a(1)=ccr*(auxa+ceps_0)
      rf6512(1)%a(2)=ccl*(auxa-ceps_0)
      rf6512(1)%b(1)=-ccl*(p6(2)+ceps_2)
      rf6512(1)%b(2)=ccr*(p6(2)-ceps_2)
* TR0 -- qu=p512,qd=p6,v=2,a=rf6512(2)%a,b=rf6512(2)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p512k0*p6(3)+p6k0*p512(3)
      ceps_0=eps_0*cim
      auxa=p512k0*p6(2)+p6k0*p512(2)
      rf6512(2)%a(1)=ccr*(auxa+ceps_0)
      rf6512(2)%a(2)=ccl*(auxa-ceps_0)
      rf6512(2)%b(1)=-ccl*p6k0
      rf6512(2)%b(2)=ccr*p6k0
* TR0 -- qu=p512,qd=p6,v=3,a=rf6512(3)%a,b=rf6512(3)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p512k0*p6(2)-p6k0*p512(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p512k0*p6(3)+p6k0*p512(3)
      rf6512(3)%a(1)=ccr*(auxa+ceps_0)
      rf6512(3)%a(2)=ccl*(auxa-ceps_0)
      rf6512(3)%b(1)=-ccl*ceps_2
      rf6512(3)%b(2)=-ccr*ceps_2
  
      cden=(-s5612)*f612
* quqd -- p=p5,q=p612
      quqd=p5(0)*p612(0)-p5(1)*p612(1)-p5(2)*p612(2)-p5(3)*p612(
     & 3)
      ccr=fcr(id5)/cden
      ccl=fcl(id5)/cden
* TL0 -- qu=p5,qd=p612,v=0,a=lf5612(0)%a,c=lf5612(0)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p5(2)*p612(3)+p612(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p612(0)+p612k0*p5(0)
      lf5612(0)%a(1)=ccr*(auxa+ceps_0)
      lf5612(0)%a(2)=ccl*(auxa-ceps_0)
      lf5612(0)%c(1)=ccr*(p5(2)+ceps_1)
      lf5612(0)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p612,v=1,a=lf5612(1)%a,c=lf5612(1)%c,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p5k0*p612(1)+p612k0*p5(1)
      lf5612(1)%a(1)=ccr*(auxa+ceps_0)
      lf5612(1)%a(2)=ccl*(auxa-ceps_0)
      lf5612(1)%c(1)=ccr*(p5(2)+ceps_1)
      lf5612(1)%c(2)=ccl*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p612,v=2,a=lf5612(2)%a,c=lf5612(2)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p5k0*p612(3)+p612k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p612(2)+p612k0*p5(2)
      lf5612(2)%a(1)=ccr*(auxa+ceps_0)
      lf5612(2)%a(2)=ccl*(auxa-ceps_0)
      lf5612(2)%c(1)=ccr*p5k0
      lf5612(2)%c(2)=-ccl*p5k0
* TL0 -- qu=p5,qd=p612,v=3,a=lf5612(3)%a,c=lf5612(3)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p5k0*p612(2)-p612k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p612(3)+p612k0*p5(3)
      lf5612(3)%a(1)=ccr*(auxa+ceps_0)
      lf5612(3)%a(2)=ccl*(auxa-ceps_0)
      lf5612(3)%c(1)=ccr*ceps_1
      lf5612(3)%c(2)=ccl*ceps_1
  
* quqd -- p=p51,q=p512
      quqd=p51(0)*p512(0)-p51(1)*p512(1)-p51(2)*p512(2)-p51(3)*p
     & 512(3)
      ccr=1.d0/(f512)
      ccl=1.d0/(f512)
      do i2=1,2
* T0 -- qu=p51,qd=p512,v=ce2(i2)%e,a=u51_2(i2)%a,b=u51_2(i2)%b,c=u51_2(i
* 2)%c,d=u51_2(i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i2)%ek0*(p51(2)*p512(3)-p512(2)*p51(3))+p51k0*
     & (ce2(i2)%e(2)*p512(3)-p512(2)*ce2(i2)%e(3))-p512k0*(ce2(i
     & 2)%e(2)*p51(3)-p51(2)*ce2(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i2)%e(3)*p51k0+p51(3)*ce2(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i2)%e(3)*p512k0+p512(3)*ce2(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i2)%e(0)*p51(0)-ce2(i2)%e(1)*p51(1)-ce2(i2)%e(2)*
     & p51(2)-ce2(i2)%e(3)*p51(3)
      cvqd=ce2(i2)%e(0)*p512(0)-ce2(i2)%e(1)*p512(1)-ce2(i2)%e(2
     & )*p512(2)-ce2(i2)%e(3)*p512(3)
      cauxa=-ce2(i2)%ek0*quqd+p51k0*cvqd+p512k0*cvqu
      cauxb=-ce2(i2)%ek0*p512(2)+p512k0*ce2(i2)%e(2)
      cauxc=+ce2(i2)%ek0*p51(2)-p51k0*ce2(i2)%e(2)
      u51_2(i2)%a(1)=ccr*(cauxa+ceps_0)
      u51_2(i2)%a(2)=ccl*(cauxa-ceps_0)
      u51_2(i2)%b(1)=ccl*(cauxb-ceps_2)
      u51_2(i2)%b(2)=ccr*(-cauxb-ceps_2)
      u51_2(i2)%c(1)=ccr*(cauxc+ceps_1)
      u51_2(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u51_2(i2)%d(1)=ccl*ce2(i2)%ek0
      u51_2(i2)%d(2)=ccr*ce2(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=l5_12(i1,i2)%a,cc=l5_12(i1,i2)%c,a1=l5_1(i1)%a,c1=l5_1(i1)%
* c,a2=u51_2(i2)%a,b2=u51_2(i2)%b,c2=u51_2(i2)%c,d2=u51_2(i2)%d,prq=s51,
* nsum=0
      l5_12(i1,i2)%a(1)=l5_1(i1)%a(1)*u51_2(i2)%a(1)+l5_1(i1)%c(
     & 1)*s51*u51_2(i2)%b(2)
      l5_12(i1,i2)%c(1)=l5_1(i1)%a(1)*u51_2(i2)%c(1)+l5_1(i1)%c(
     & 1)*s51*u51_2(i2)%d(2)
      l5_12(i1,i2)%c(2)=l5_1(i1)%c(2)*s51*u51_2(i2)%d(1)+l5_1(i1
     & )%a(2)*u51_2(i2)%c(2)
      l5_12(i1,i2)%a(2)=l5_1(i1)%c(2)*s51*u51_2(i2)%b(1)+l5_1(i1
     & )%a(2)*u51_2(i2)%a(2)
      end do
      end do
  
* quqd -- p=p612,q=p62
      quqd=p612(0)*p62(0)-p612(1)*p62(1)-p612(2)*p62(2)-p612(3)*
     & p62(3)
      ccr=1.d0/(f62)
      ccl=1.d0/(f62)
      do i1=1,2
* T0 -- qu=p612,qd=p62,v=ce1(i1)%e,a=u62_1(i1)%a,b=u62_1(i1)%b,c=u62_1(i
* 1)%c,d=u62_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p612(2)*p62(3)-p62(2)*p612(3))+p612k0
     & *(ce1(i1)%e(2)*p62(3)-p62(2)*ce1(i1)%e(3))-p62k0*(ce1(i1)
     & %e(2)*p612(3)-p612(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p612k0+p612(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p62k0+p62(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p612(0)-ce1(i1)%e(1)*p612(1)-ce1(i1)%e(2
     & )*p612(2)-ce1(i1)%e(3)*p612(3)
      cvqd=ce1(i1)%e(0)*p62(0)-ce1(i1)%e(1)*p62(1)-ce1(i1)%e(2)*
     & p62(2)-ce1(i1)%e(3)*p62(3)
      cauxa=-ce1(i1)%ek0*quqd+p612k0*cvqd+p62k0*cvqu
      cauxb=-ce1(i1)%ek0*p62(2)+p62k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p612(2)-p612k0*ce1(i1)%e(2)
      u62_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u62_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u62_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u62_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u62_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u62_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u62_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u62_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0 -- aa=r6_12(i1,i2)%a,bb=r6_12(i1,i2)%b,a1=u62_1(i1)%a,b1=u62_1(i1
* )%b,c1=u62_1(i1)%c,d1=u62_1(i1)%d,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s62,
* nsum=0
      r6_12(i1,i2)%a(1)=u62_1(i1)%a(1)*r6_2(i2)%a(1)+u62_1(i1)%c
     & (1)*s62*r6_2(i2)%b(2)
      r6_12(i1,i2)%b(1)=u62_1(i1)%d(1)*s62*r6_2(i2)%b(1)+u62_1(i
     & 1)%b(1)*r6_2(i2)%a(2)
      r6_12(i1,i2)%b(2)=u62_1(i1)%b(2)*r6_2(i2)%a(1)+u62_1(i1)%d
     & (2)*s62*r6_2(i2)%b(2)
      r6_12(i1,i2)%a(2)=u62_1(i1)%c(2)*s62*r6_2(i2)%b(1)+u62_1(i
     & 1)%a(2)*r6_2(i2)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l5_12(i1,i2)%a(1)=l5_12(i1,i2)%a(1) +
     &   l5_gg(i1,i2)%a(1)
       l5_12(i1,i2)%a(2)=l5_12(i1,i2)%a(2) +
     &   l5_gg(i1,i2)%a(2)
       l5_12(i1,i2)%c(1)=l5_12(i1,i2)%c(1) +
     &   l5_gg(i1,i2)%c(1)
       l5_12(i1,i2)%c(2)=l5_12(i1,i2)%c(2) +
     &   l5_gg(i1,i2)%c(2)
  
       r6_12(i1,i2)%a(1)=r6_12(i1,i2)%a(1) +
     &   r6_gg(i1,i2)%a(1)
       r6_12(i1,i2)%a(2)=r6_12(i1,i2)%a(2) +
     &   r6_gg(i1,i2)%a(2)
       r6_12(i1,i2)%b(1)=r6_12(i1,i2)%b(1) +
     &   r6_gg(i1,i2)%b(1)
       r6_12(i1,i2)%b(2)=r6_12(i1,i2)%b(2) +
     &   r6_gg(i1,i2)%b(2)
  
      enddo
      enddo
  
  
      cden=(-s5612+cmz2)*f62
* quqd -- p=p51,q=p62
      quqd=p51(0)*p62(0)-p51(1)*p62(1)-p51(2)*p62(2)-p51(3)*p62(
     & 3)
      ccr=zcr(id5)/cden
      ccl=zcl(id5)/cden
* T0 -- qu=p51,qd=p62,v=0,a=uz51(0)%a,b=uz51(0)%b,c=uz51(0)%c,d=uz51(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p51(2)*p62(3)+p62(2)*p51(3)
      ceps_0=eps_0*cim
      ceps_1=p51(3)*cim
      ceps_2=p62(3)*cim
      auxa=-quqd+p51k0*p62(0)+p62k0*p51(0)
      uz51(0)%a(1)=ccr*(auxa+ceps_0)
      uz51(0)%a(2)=ccl*(auxa-ceps_0)
      uz51(0)%b(1)=-ccl*(p62(2)+ceps_2)
      uz51(0)%b(2)=ccr*(p62(2)-ceps_2)
      uz51(0)%c(1)=ccr*(p51(2)+ceps_1)
      uz51(0)%c(2)=ccl*(-p51(2)+ceps_1)
      uz51(0)%d(1)=ccl
      uz51(0)%d(2)=ccr
* T0 -- qu=p51,qd=p62,v=1,a=uz51(1)%a,b=uz51(1)%b,c=uz51(1)%c,d=uz51(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p51k0*p62(1)+p62k0*p51(1)
      uz51(1)%a(1)=ccr*(auxa+ceps_0)
      uz51(1)%a(2)=ccl*(auxa-ceps_0)
      uz51(1)%b(1)=-ccl*(p62(2)+ceps_2)
      uz51(1)%b(2)=ccr*(p62(2)-ceps_2)
      uz51(1)%c(1)=ccr*(p51(2)+ceps_1)
      uz51(1)%c(2)=ccl*(-p51(2)+ceps_1)
      uz51(1)%d(1)=ccl
      uz51(1)%d(2)=ccr
* T0 -- qu=p51,qd=p62,v=2,a=uz51(2)%a,b=uz51(2)%b,c=uz51(2)%c,d=uz51(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p51k0*p62(3)+p62k0*p51(3)
      ceps_0=eps_0*cim
      auxa=p51k0*p62(2)+p62k0*p51(2)
      uz51(2)%a(1)=ccr*(auxa+ceps_0)
      uz51(2)%a(2)=ccl*(auxa-ceps_0)
      uz51(2)%b(1)=-ccl*p62k0
      uz51(2)%b(2)=ccr*p62k0
      uz51(2)%c(1)=ccr*p51k0
      uz51(2)%c(2)=-ccl*p51k0
* T0 -- qu=p51,qd=p62,v=3,a=uz51(3)%a,b=uz51(3)%b,c=uz51(3)%c,d=uz51(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p51k0*p62(2)-p62k0*p51(2)
      ceps_0=eps_0*cim
      ceps_1=p51k0*cim
      ceps_2=p62k0*cim
      auxa=+p51k0*p62(3)+p62k0*p51(3)
      uz51(3)%a(1)=ccr*(auxa+ceps_0)
      uz51(3)%a(2)=ccl*(auxa-ceps_0)
      uz51(3)%b(1)=-ccl*ceps_2
      uz51(3)%b(2)=-ccr*ceps_2
      uz51(3)%c(1)=ccr*ceps_1
      uz51(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l5_1(i1)%a,c1=l5_
* 1(i1)%c,a2=uz51(mu)%a,b2=uz51(mu)%b,c2=uz51(mu)%c,d2=uz51(mu)%d,prq=s5
* 1,nsum=0
      laux_imu(i1,mu)%a(1)=l5_1(i1)%a(1)*uz51(mu)%a(1)+l5_1(i1)%
     & c(1)*s51*uz51(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l5_1(i1)%a(1)*uz51(mu)%c(1)+l5_1(i1)%
     & c(1)*s51*uz51(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l5_1(i1)%c(2)*s51*uz51(mu)%d(1)+l5_1(
     & i1)%a(2)*uz51(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l5_1(i1)%c(2)*s51*uz51(mu)%b(1)+l5_1(
     & i1)%a(2)*uz51(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz5126(&,i1,i2)%e(mu),a1=laux_imu(i1,mu)%a,c1=laux_imu(i1,
* mu)%c,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s62,bef=,aft=
      cz5126(1,i1,i2)%e(mu)=(laux_imu(i1,mu)%a(1)*r6_2(i2)%a(1)+
     & laux_imu(i1,mu)%c(1)*s62*r6_2(i2)%b(2))
      cz5126(2,i1,i2)%e(mu)=(laux_imu(i1,mu)%c(2)*s62*r6_2(i2)%b
     & (1)+laux_imu(i1,mu)%a(2)*r6_2(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz5126(&,i1,i2)%e(mu),a1=l5_12(i1,i2)%a,c1=l5_12(i1,i2)%c,
* a2=rz6512(mu)%a,b2=rz6512(mu)%b,prq=s512,bef=cz5126(&,i1,i2)%e(mu)+,af
* t=
      cz5126(1,i1,i2)%e(mu)=cz5126(1,i1,i2)%e(mu)+(l5_12(i1,i2)%
     & a(1)*rz6512(mu)%a(1)+l5_12(i1,i2)%c(1)*s512*rz6512(mu)%b(
     & 2))
      cz5126(2,i1,i2)%e(mu)=cz5126(2,i1,i2)%e(mu)+(l5_12(i1,i2)%
     & c(2)*s512*rz6512(mu)%b(1)+l5_12(i1,i2)%a(2)*rz6512(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz5126(&,i1,i2)%e(mu),a1=lz5612(mu)%a,c1=lz5612(mu)%c,a2=r
* 6_12(i1,i2)%a,b2=r6_12(i1,i2)%b,prq=s612,bef=cz5126(&,i1,i2)%e(mu)+,af
* t=
      cz5126(1,i1,i2)%e(mu)=cz5126(1,i1,i2)%e(mu)+(lz5612(mu)%a(
     & 1)*r6_12(i1,i2)%a(1)+lz5612(mu)%c(1)*s612*r6_12(i1,i2)%b(
     & 2))
      cz5126(2,i1,i2)%e(mu)=cz5126(2,i1,i2)%e(mu)+(lz5612(mu)%c(
     & 2)*s612*r6_12(i1,i2)%b(1)+lz5612(mu)%a(2)*r6_12(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cz5126(i3,i1,i2)%e
      cz5126(i3,i1,i2)%ek0=cz5126(i3,i1,i2)%e(0)-cz5126(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
      cden=(-s5612)*f62
* quqd -- p=p51,q=p62
      quqd=p51(0)*p62(0)-p51(1)*p62(1)-p51(2)*p62(2)-p51(3)*p62(
     & 3)
      ccr=fcr(id5)/cden
      ccl=fcl(id5)/cden
* T0 -- qu=p51,qd=p62,v=0,a=uf51(0)%a,b=uf51(0)%b,c=uf51(0)%c,d=uf51(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p51(2)*p62(3)+p62(2)*p51(3)
      ceps_0=eps_0*cim
      ceps_1=p51(3)*cim
      ceps_2=p62(3)*cim
      auxa=-quqd+p51k0*p62(0)+p62k0*p51(0)
      uf51(0)%a(1)=ccr*(auxa+ceps_0)
      uf51(0)%a(2)=ccl*(auxa-ceps_0)
      uf51(0)%b(1)=-ccl*(p62(2)+ceps_2)
      uf51(0)%b(2)=ccr*(p62(2)-ceps_2)
      uf51(0)%c(1)=ccr*(p51(2)+ceps_1)
      uf51(0)%c(2)=ccl*(-p51(2)+ceps_1)
      uf51(0)%d(1)=ccl
      uf51(0)%d(2)=ccr
* T0 -- qu=p51,qd=p62,v=1,a=uf51(1)%a,b=uf51(1)%b,c=uf51(1)%c,d=uf51(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p51k0*p62(1)+p62k0*p51(1)
      uf51(1)%a(1)=ccr*(auxa+ceps_0)
      uf51(1)%a(2)=ccl*(auxa-ceps_0)
      uf51(1)%b(1)=-ccl*(p62(2)+ceps_2)
      uf51(1)%b(2)=ccr*(p62(2)-ceps_2)
      uf51(1)%c(1)=ccr*(p51(2)+ceps_1)
      uf51(1)%c(2)=ccl*(-p51(2)+ceps_1)
      uf51(1)%d(1)=ccl
      uf51(1)%d(2)=ccr
* T0 -- qu=p51,qd=p62,v=2,a=uf51(2)%a,b=uf51(2)%b,c=uf51(2)%c,d=uf51(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p51k0*p62(3)+p62k0*p51(3)
      ceps_0=eps_0*cim
      auxa=p51k0*p62(2)+p62k0*p51(2)
      uf51(2)%a(1)=ccr*(auxa+ceps_0)
      uf51(2)%a(2)=ccl*(auxa-ceps_0)
      uf51(2)%b(1)=-ccl*p62k0
      uf51(2)%b(2)=ccr*p62k0
      uf51(2)%c(1)=ccr*p51k0
      uf51(2)%c(2)=-ccl*p51k0
* T0 -- qu=p51,qd=p62,v=3,a=uf51(3)%a,b=uf51(3)%b,c=uf51(3)%c,d=uf51(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p51k0*p62(2)-p62k0*p51(2)
      ceps_0=eps_0*cim
      ceps_1=p51k0*cim
      ceps_2=p62k0*cim
      auxa=+p51k0*p62(3)+p62k0*p51(3)
      uf51(3)%a(1)=ccr*(auxa+ceps_0)
      uf51(3)%a(2)=ccl*(auxa-ceps_0)
      uf51(3)%b(1)=-ccl*ceps_2
      uf51(3)%b(2)=-ccr*ceps_2
      uf51(3)%c(1)=ccr*ceps_1
      uf51(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l5_1(i1)%a,c1=l5_
* 1(i1)%c,a2=uf51(mu)%a,b2=uf51(mu)%b,c2=uf51(mu)%c,d2=uf51(mu)%d,prq=s5
* 1,nsum=0
      laux_imu(i1,mu)%a(1)=l5_1(i1)%a(1)*uf51(mu)%a(1)+l5_1(i1)%
     & c(1)*s51*uf51(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l5_1(i1)%a(1)*uf51(mu)%c(1)+l5_1(i1)%
     & c(1)*s51*uf51(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l5_1(i1)%c(2)*s51*uf51(mu)%d(1)+l5_1(
     & i1)%a(2)*uf51(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l5_1(i1)%c(2)*s51*uf51(mu)%b(1)+l5_1(
     & i1)%a(2)*uf51(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf5126(&,i1,i2)%e(mu),a1=laux_imu(i1,mu)%a,c1=laux_imu(i1,
* mu)%c,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s62,bef=,aft=
      cf5126(1,i1,i2)%e(mu)=(laux_imu(i1,mu)%a(1)*r6_2(i2)%a(1)+
     & laux_imu(i1,mu)%c(1)*s62*r6_2(i2)%b(2))
      cf5126(2,i1,i2)%e(mu)=(laux_imu(i1,mu)%c(2)*s62*r6_2(i2)%b
     & (1)+laux_imu(i1,mu)%a(2)*r6_2(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf5126(&,i1,i2)%e(mu),a1=l5_12(i1,i2)%a,c1=l5_12(i1,i2)%c,
* a2=rf6512(mu)%a,b2=rf6512(mu)%b,prq=s512,bef=cf5126(&,i1,i2)%e(mu)+,af
* t=
      cf5126(1,i1,i2)%e(mu)=cf5126(1,i1,i2)%e(mu)+(l5_12(i1,i2)%
     & a(1)*rf6512(mu)%a(1)+l5_12(i1,i2)%c(1)*s512*rf6512(mu)%b(
     & 2))
      cf5126(2,i1,i2)%e(mu)=cf5126(2,i1,i2)%e(mu)+(l5_12(i1,i2)%
     & c(2)*s512*rf6512(mu)%b(1)+l5_12(i1,i2)%a(2)*rf6512(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf5126(&,i1,i2)%e(mu),a1=lf5612(mu)%a,c1=lf5612(mu)%c,a2=r
* 6_12(i1,i2)%a,b2=r6_12(i1,i2)%b,prq=s612,bef=cf5126(&,i1,i2)%e(mu)+,af
* t=
      cf5126(1,i1,i2)%e(mu)=cf5126(1,i1,i2)%e(mu)+(lf5612(mu)%a(
     & 1)*r6_12(i1,i2)%a(1)+lf5612(mu)%c(1)*s612*r6_12(i1,i2)%b(
     & 2))
      cf5126(2,i1,i2)%e(mu)=cf5126(2,i1,i2)%e(mu)+(lf5612(mu)%c(
     & 2)*s612*r6_12(i1,i2)%b(1)+lf5612(mu)%a(2)*r6_12(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cf5126(i3,i1,i2)%e
      cf5126(i3,i1,i2)%ek0=cf5126(i3,i1,i2)%e(0)-cf5126(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
* quqd -- p=p52,q=p512
      quqd=p52(0)*p512(0)-p52(1)*p512(1)-p52(2)*p512(2)-p52(3)*p
     & 512(3)
      ccr=1.d0/(f512)
      ccl=1.d0/(f512)
      do i2=1,2
* T0 -- qu=p52,qd=p512,v=ce1(i2)%e,a=u52_1(i2)%a,b=u52_1(i2)%b,c=u52_1(i
* 2)%c,d=u52_1(i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i2)%ek0*(p52(2)*p512(3)-p512(2)*p52(3))+p52k0*
     & (ce1(i2)%e(2)*p512(3)-p512(2)*ce1(i2)%e(3))-p512k0*(ce1(i
     & 2)%e(2)*p52(3)-p52(2)*ce1(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i2)%e(3)*p52k0+p52(3)*ce1(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i2)%e(3)*p512k0+p512(3)*ce1(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i2)%e(0)*p52(0)-ce1(i2)%e(1)*p52(1)-ce1(i2)%e(2)*
     & p52(2)-ce1(i2)%e(3)*p52(3)
      cvqd=ce1(i2)%e(0)*p512(0)-ce1(i2)%e(1)*p512(1)-ce1(i2)%e(2
     & )*p512(2)-ce1(i2)%e(3)*p512(3)
      cauxa=-ce1(i2)%ek0*quqd+p52k0*cvqd+p512k0*cvqu
      cauxb=-ce1(i2)%ek0*p512(2)+p512k0*ce1(i2)%e(2)
      cauxc=+ce1(i2)%ek0*p52(2)-p52k0*ce1(i2)%e(2)
      u52_1(i2)%a(1)=ccr*(cauxa+ceps_0)
      u52_1(i2)%a(2)=ccl*(cauxa-ceps_0)
      u52_1(i2)%b(1)=ccl*(cauxb-ceps_2)
      u52_1(i2)%b(2)=ccr*(-cauxb-ceps_2)
      u52_1(i2)%c(1)=ccr*(cauxc+ceps_1)
      u52_1(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u52_1(i2)%d(1)=ccl*ce1(i2)%ek0
      u52_1(i2)%d(2)=ccr*ce1(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=l5_21(i1,i2)%a,cc=l5_21(i1,i2)%c,a1=l5_2(i2)%a,c1=l5_2(i2)%
* c,a2=u52_1(i1)%a,b2=u52_1(i1)%b,c2=u52_1(i1)%c,d2=u52_1(i1)%d,prq=s52,
* nsum=0
      l5_21(i1,i2)%a(1)=l5_2(i2)%a(1)*u52_1(i1)%a(1)+l5_2(i2)%c(
     & 1)*s52*u52_1(i1)%b(2)
      l5_21(i1,i2)%c(1)=l5_2(i2)%a(1)*u52_1(i1)%c(1)+l5_2(i2)%c(
     & 1)*s52*u52_1(i1)%d(2)
      l5_21(i1,i2)%c(2)=l5_2(i2)%c(2)*s52*u52_1(i1)%d(1)+l5_2(i2
     & )%a(2)*u52_1(i1)%c(2)
      l5_21(i1,i2)%a(2)=l5_2(i2)%c(2)*s52*u52_1(i1)%b(1)+l5_2(i2
     & )%a(2)*u52_1(i1)%a(2)
      end do
      end do
  
* quqd -- p=p612,q=p61
      quqd=p612(0)*p61(0)-p612(1)*p61(1)-p612(2)*p61(2)-p612(3)*
     & p61(3)
      ccr=1.d0/(f61)
      ccl=1.d0/(f61)
      do i1=1,2
* T0 -- qu=p612,qd=p61,v=ce2(i1)%e,a=u61_2(i1)%a,b=u61_2(i1)%b,c=u61_2(i
* 1)%c,d=u61_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p612(2)*p61(3)-p61(2)*p612(3))+p612k0
     & *(ce2(i1)%e(2)*p61(3)-p61(2)*ce2(i1)%e(3))-p61k0*(ce2(i1)
     & %e(2)*p612(3)-p612(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p612k0+p612(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p61k0+p61(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p612(0)-ce2(i1)%e(1)*p612(1)-ce2(i1)%e(2
     & )*p612(2)-ce2(i1)%e(3)*p612(3)
      cvqd=ce2(i1)%e(0)*p61(0)-ce2(i1)%e(1)*p61(1)-ce2(i1)%e(2)*
     & p61(2)-ce2(i1)%e(3)*p61(3)
      cauxa=-ce2(i1)%ek0*quqd+p612k0*cvqd+p61k0*cvqu
      cauxb=-ce2(i1)%ek0*p61(2)+p61k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p612(2)-p612k0*ce2(i1)%e(2)
      u61_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u61_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u61_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u61_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u61_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u61_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u61_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u61_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0 -- aa=r6_21(i1,i2)%a,bb=r6_21(i1,i2)%b,a1=u61_2(i2)%a,b1=u61_2(i2
* )%b,c1=u61_2(i2)%c,d1=u61_2(i2)%d,a2=r6_1(i1)%a,b2=r6_1(i1)%b,prq=s61,
* nsum=0
      r6_21(i1,i2)%a(1)=u61_2(i2)%a(1)*r6_1(i1)%a(1)+u61_2(i2)%c
     & (1)*s61*r6_1(i1)%b(2)
      r6_21(i1,i2)%b(1)=u61_2(i2)%d(1)*s61*r6_1(i1)%b(1)+u61_2(i
     & 2)%b(1)*r6_1(i1)%a(2)
      r6_21(i1,i2)%b(2)=u61_2(i2)%b(2)*r6_1(i1)%a(1)+u61_2(i2)%d
     & (2)*s61*r6_1(i1)%b(2)
      r6_21(i1,i2)%a(2)=u61_2(i2)%c(2)*s61*r6_1(i1)%b(1)+u61_2(i
     & 2)%a(2)*r6_1(i1)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l5_21(i1,i2)%a(1)=l5_21(i1,i2)%a(1) -
     &   l5_gg(i1,i2)%a(1)
       l5_21(i1,i2)%a(2)=l5_21(i1,i2)%a(2) -
     &   l5_gg(i1,i2)%a(2)
       l5_21(i1,i2)%c(1)=l5_21(i1,i2)%c(1) -
     &   l5_gg(i1,i2)%c(1)
       l5_21(i1,i2)%c(2)=l5_21(i1,i2)%c(2) -
     &   l5_gg(i1,i2)%c(2)
  
       r6_21(i1,i2)%a(1)=r6_21(i1,i2)%a(1) -
     &   r6_gg(i1,i2)%a(1)
       r6_21(i1,i2)%a(2)=r6_21(i1,i2)%a(2) -
     &   r6_gg(i1,i2)%a(2)
       r6_21(i1,i2)%b(1)=r6_21(i1,i2)%b(1) -
     &   r6_gg(i1,i2)%b(1)
       r6_21(i1,i2)%b(2)=r6_21(i1,i2)%b(2) -
     &   r6_gg(i1,i2)%b(2)
  
      enddo
      enddo
  
  
      cden=(-s5612+cmz2)*f61
* quqd -- p=p52,q=p61
      quqd=p52(0)*p61(0)-p52(1)*p61(1)-p52(2)*p61(2)-p52(3)*p61(
     & 3)
      ccr=zcr(id5)/cden
      ccl=zcl(id5)/cden
* T0 -- qu=p52,qd=p61,v=0,a=uz52(0)%a,b=uz52(0)%b,c=uz52(0)%c,d=uz52(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p52(2)*p61(3)+p61(2)*p52(3)
      ceps_0=eps_0*cim
      ceps_1=p52(3)*cim
      ceps_2=p61(3)*cim
      auxa=-quqd+p52k0*p61(0)+p61k0*p52(0)
      uz52(0)%a(1)=ccr*(auxa+ceps_0)
      uz52(0)%a(2)=ccl*(auxa-ceps_0)
      uz52(0)%b(1)=-ccl*(p61(2)+ceps_2)
      uz52(0)%b(2)=ccr*(p61(2)-ceps_2)
      uz52(0)%c(1)=ccr*(p52(2)+ceps_1)
      uz52(0)%c(2)=ccl*(-p52(2)+ceps_1)
      uz52(0)%d(1)=ccl
      uz52(0)%d(2)=ccr
* T0 -- qu=p52,qd=p61,v=1,a=uz52(1)%a,b=uz52(1)%b,c=uz52(1)%c,d=uz52(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p52k0*p61(1)+p61k0*p52(1)
      uz52(1)%a(1)=ccr*(auxa+ceps_0)
      uz52(1)%a(2)=ccl*(auxa-ceps_0)
      uz52(1)%b(1)=-ccl*(p61(2)+ceps_2)
      uz52(1)%b(2)=ccr*(p61(2)-ceps_2)
      uz52(1)%c(1)=ccr*(p52(2)+ceps_1)
      uz52(1)%c(2)=ccl*(-p52(2)+ceps_1)
      uz52(1)%d(1)=ccl
      uz52(1)%d(2)=ccr
* T0 -- qu=p52,qd=p61,v=2,a=uz52(2)%a,b=uz52(2)%b,c=uz52(2)%c,d=uz52(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p52k0*p61(3)+p61k0*p52(3)
      ceps_0=eps_0*cim
      auxa=p52k0*p61(2)+p61k0*p52(2)
      uz52(2)%a(1)=ccr*(auxa+ceps_0)
      uz52(2)%a(2)=ccl*(auxa-ceps_0)
      uz52(2)%b(1)=-ccl*p61k0
      uz52(2)%b(2)=ccr*p61k0
      uz52(2)%c(1)=ccr*p52k0
      uz52(2)%c(2)=-ccl*p52k0
* T0 -- qu=p52,qd=p61,v=3,a=uz52(3)%a,b=uz52(3)%b,c=uz52(3)%c,d=uz52(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p52k0*p61(2)-p61k0*p52(2)
      ceps_0=eps_0*cim
      ceps_1=p52k0*cim
      ceps_2=p61k0*cim
      auxa=+p52k0*p61(3)+p61k0*p52(3)
      uz52(3)%a(1)=ccr*(auxa+ceps_0)
      uz52(3)%a(2)=ccl*(auxa-ceps_0)
      uz52(3)%b(1)=-ccl*ceps_2
      uz52(3)%b(2)=-ccr*ceps_2
      uz52(3)%c(1)=ccr*ceps_1
      uz52(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l5_2(i1)%a,c1=l5_
* 2(i1)%c,a2=uz52(mu)%a,b2=uz52(mu)%b,c2=uz52(mu)%c,d2=uz52(mu)%d,prq=s5
* 2,nsum=0
      laux_imu(i1,mu)%a(1)=l5_2(i1)%a(1)*uz52(mu)%a(1)+l5_2(i1)%
     & c(1)*s52*uz52(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l5_2(i1)%a(1)*uz52(mu)%c(1)+l5_2(i1)%
     & c(1)*s52*uz52(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l5_2(i1)%c(2)*s52*uz52(mu)%d(1)+l5_2(
     & i1)%a(2)*uz52(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l5_2(i1)%c(2)*s52*uz52(mu)%b(1)+l5_2(
     & i1)%a(2)*uz52(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz5216(&,i1,i2)%e(mu),a1=laux_imu(i2,mu)%a,c1=laux_imu(i2,
* mu)%c,a2=r6_1(i1)%a,b2=r6_1(i1)%b,prq=s61,bef=,aft=
      cz5216(1,i1,i2)%e(mu)=(laux_imu(i2,mu)%a(1)*r6_1(i1)%a(1)+
     & laux_imu(i2,mu)%c(1)*s61*r6_1(i1)%b(2))
      cz5216(2,i1,i2)%e(mu)=(laux_imu(i2,mu)%c(2)*s61*r6_1(i1)%b
     & (1)+laux_imu(i2,mu)%a(2)*r6_1(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz5216(&,i1,i2)%e(mu),a1=l5_21(i1,i2)%a,c1=l5_21(i1,i2)%c,
* a2=rz6512(mu)%a,b2=rz6512(mu)%b,prq=s512,bef=cz5216(&,i1,i2)%e(mu)+,af
* t=
      cz5216(1,i1,i2)%e(mu)=cz5216(1,i1,i2)%e(mu)+(l5_21(i1,i2)%
     & a(1)*rz6512(mu)%a(1)+l5_21(i1,i2)%c(1)*s512*rz6512(mu)%b(
     & 2))
      cz5216(2,i1,i2)%e(mu)=cz5216(2,i1,i2)%e(mu)+(l5_21(i1,i2)%
     & c(2)*s512*rz6512(mu)%b(1)+l5_21(i1,i2)%a(2)*rz6512(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz5216(&,i1,i2)%e(mu),a1=lz5612(mu)%a,c1=lz5612(mu)%c,a2=r
* 6_21(i1,i2)%a,b2=r6_21(i1,i2)%b,prq=s612,bef=cz5216(&,i1,i2)%e(mu)+,af
* t=
      cz5216(1,i1,i2)%e(mu)=cz5216(1,i1,i2)%e(mu)+(lz5612(mu)%a(
     & 1)*r6_21(i1,i2)%a(1)+lz5612(mu)%c(1)*s612*r6_21(i1,i2)%b(
     & 2))
      cz5216(2,i1,i2)%e(mu)=cz5216(2,i1,i2)%e(mu)+(lz5612(mu)%c(
     & 2)*s612*r6_21(i1,i2)%b(1)+lz5612(mu)%a(2)*r6_21(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cz5216(i3,i1,i2)%e
      cz5216(i3,i1,i2)%ek0=cz5216(i3,i1,i2)%e(0)-cz5216(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
      cden=(-s5612)*f61
* quqd -- p=p52,q=p61
      quqd=p52(0)*p61(0)-p52(1)*p61(1)-p52(2)*p61(2)-p52(3)*p61(
     & 3)
      ccr=fcr(id5)/cden
      ccl=fcl(id5)/cden
* T0 -- qu=p52,qd=p61,v=0,a=uf52(0)%a,b=uf52(0)%b,c=uf52(0)%c,d=uf52(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p52(2)*p61(3)+p61(2)*p52(3)
      ceps_0=eps_0*cim
      ceps_1=p52(3)*cim
      ceps_2=p61(3)*cim
      auxa=-quqd+p52k0*p61(0)+p61k0*p52(0)
      uf52(0)%a(1)=ccr*(auxa+ceps_0)
      uf52(0)%a(2)=ccl*(auxa-ceps_0)
      uf52(0)%b(1)=-ccl*(p61(2)+ceps_2)
      uf52(0)%b(2)=ccr*(p61(2)-ceps_2)
      uf52(0)%c(1)=ccr*(p52(2)+ceps_1)
      uf52(0)%c(2)=ccl*(-p52(2)+ceps_1)
      uf52(0)%d(1)=ccl
      uf52(0)%d(2)=ccr
* T0 -- qu=p52,qd=p61,v=1,a=uf52(1)%a,b=uf52(1)%b,c=uf52(1)%c,d=uf52(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p52k0*p61(1)+p61k0*p52(1)
      uf52(1)%a(1)=ccr*(auxa+ceps_0)
      uf52(1)%a(2)=ccl*(auxa-ceps_0)
      uf52(1)%b(1)=-ccl*(p61(2)+ceps_2)
      uf52(1)%b(2)=ccr*(p61(2)-ceps_2)
      uf52(1)%c(1)=ccr*(p52(2)+ceps_1)
      uf52(1)%c(2)=ccl*(-p52(2)+ceps_1)
      uf52(1)%d(1)=ccl
      uf52(1)%d(2)=ccr
* T0 -- qu=p52,qd=p61,v=2,a=uf52(2)%a,b=uf52(2)%b,c=uf52(2)%c,d=uf52(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p52k0*p61(3)+p61k0*p52(3)
      ceps_0=eps_0*cim
      auxa=p52k0*p61(2)+p61k0*p52(2)
      uf52(2)%a(1)=ccr*(auxa+ceps_0)
      uf52(2)%a(2)=ccl*(auxa-ceps_0)
      uf52(2)%b(1)=-ccl*p61k0
      uf52(2)%b(2)=ccr*p61k0
      uf52(2)%c(1)=ccr*p52k0
      uf52(2)%c(2)=-ccl*p52k0
* T0 -- qu=p52,qd=p61,v=3,a=uf52(3)%a,b=uf52(3)%b,c=uf52(3)%c,d=uf52(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p52k0*p61(2)-p61k0*p52(2)
      ceps_0=eps_0*cim
      ceps_1=p52k0*cim
      ceps_2=p61k0*cim
      auxa=+p52k0*p61(3)+p61k0*p52(3)
      uf52(3)%a(1)=ccr*(auxa+ceps_0)
      uf52(3)%a(2)=ccl*(auxa-ceps_0)
      uf52(3)%b(1)=-ccl*ceps_2
      uf52(3)%b(2)=-ccr*ceps_2
      uf52(3)%c(1)=ccr*ceps_1
      uf52(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l5_2(i1)%a,c1=l5_
* 2(i1)%c,a2=uf52(mu)%a,b2=uf52(mu)%b,c2=uf52(mu)%c,d2=uf52(mu)%d,prq=s5
* 2,nsum=0
      laux_imu(i1,mu)%a(1)=l5_2(i1)%a(1)*uf52(mu)%a(1)+l5_2(i1)%
     & c(1)*s52*uf52(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l5_2(i1)%a(1)*uf52(mu)%c(1)+l5_2(i1)%
     & c(1)*s52*uf52(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l5_2(i1)%c(2)*s52*uf52(mu)%d(1)+l5_2(
     & i1)%a(2)*uf52(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l5_2(i1)%c(2)*s52*uf52(mu)%b(1)+l5_2(
     & i1)%a(2)*uf52(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf5216(&,i1,i2)%e(mu),a1=laux_imu(i2,mu)%a,c1=laux_imu(i2,
* mu)%c,a2=r6_1(i1)%a,b2=r6_1(i1)%b,prq=s61,bef=,aft=
      cf5216(1,i1,i2)%e(mu)=(laux_imu(i2,mu)%a(1)*r6_1(i1)%a(1)+
     & laux_imu(i2,mu)%c(1)*s61*r6_1(i1)%b(2))
      cf5216(2,i1,i2)%e(mu)=(laux_imu(i2,mu)%c(2)*s61*r6_1(i1)%b
     & (1)+laux_imu(i2,mu)%a(2)*r6_1(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf5216(&,i1,i2)%e(mu),a1=l5_21(i1,i2)%a,c1=l5_21(i1,i2)%c,
* a2=rf6512(mu)%a,b2=rf6512(mu)%b,prq=s512,bef=cf5216(&,i1,i2)%e(mu)+,af
* t=
      cf5216(1,i1,i2)%e(mu)=cf5216(1,i1,i2)%e(mu)+(l5_21(i1,i2)%
     & a(1)*rf6512(mu)%a(1)+l5_21(i1,i2)%c(1)*s512*rf6512(mu)%b(
     & 2))
      cf5216(2,i1,i2)%e(mu)=cf5216(2,i1,i2)%e(mu)+(l5_21(i1,i2)%
     & c(2)*s512*rf6512(mu)%b(1)+l5_21(i1,i2)%a(2)*rf6512(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf5216(&,i1,i2)%e(mu),a1=lf5612(mu)%a,c1=lf5612(mu)%c,a2=r
* 6_21(i1,i2)%a,b2=r6_21(i1,i2)%b,prq=s612,bef=cf5216(&,i1,i2)%e(mu)+,af
* t=
      cf5216(1,i1,i2)%e(mu)=cf5216(1,i1,i2)%e(mu)+(lf5612(mu)%a(
     & 1)*r6_21(i1,i2)%a(1)+lf5612(mu)%c(1)*s612*r6_21(i1,i2)%b(
     & 2))
      cf5216(2,i1,i2)%e(mu)=cf5216(2,i1,i2)%e(mu)+(lf5612(mu)%c(
     & 2)*s612*r6_21(i1,i2)%b(1)+lf5612(mu)%a(2)*r6_21(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cf5216(i3,i1,i2)%e
      cf5216(i3,i1,i2)%ek0=cf5216(i3,i1,i2)%e(0)-cf5216(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
      endif !id3 = quark
  
      do m=0,3
       p712(m)=p71(m) + p2(m)
      enddo
* pk0 -- p=p712
      p712k0=p712(0)-p712(1)
* p.q -- p.q=s712,p=p712,q=p712,bef=,aft=
      s712=(p712(0)*p712(0)-p712(1)*p712(1)-p712(2)*p712(2)-p712
     & (3)*p712(3))
      f712=s712*p712k0
  
      do m=0,3
       p812(m)=p81(m) - p2(m)
      enddo
* pk0 -- p=p812
      p812k0=p812(0)-p812(1)
* p.q -- p.q=s812,p=p812,q=p812,bef=,aft=
      s812=(p812(0)*p812(0)-p812(1)*p812(1)-p812(2)*p812(2)-p812
     & (3)*p812(3))
      f812=s812*p812k0
  
* p.q -- p.q=s7812,p=p712,q=p8,bef=2.d0*,aft=+s712
      s7812=2.d0*(p712(0)*p8(0)-p712(1)*p8(1)-p712(2)*p8(2)-p712
     & (3)*p8(3))+s712
  
      if (ilept(id7).ne.1) then
  
* triple vertex in left/right line                                      
* quqd -- p=p812,q=p8
      quqd=p812(0)*p8(0)-p812(1)*p8(1)-p812(2)*p8(2)-p812(3)*p8(
     & 3)
      ccr=1.d0/(s12)
      ccl=1.d0/(s12)
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p812,qd=p8,v=ctrip12(i1,i2)%e,a=r8_gg(i1,i2)%a,b=r8_gg(i1,i2
* )%b,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p812(2)*p8(3)-p8(2)*p812(3))+p
     & 812k0*(ctrip12(i1,i2)%e(2)*p8(3)-p8(2)*ctrip12(i1,i2)%e(3
     & ))-p8k0*(ctrip12(i1,i2)%e(2)*p812(3)-p812(2)*ctrip12(i1,i
     & 2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-ctrip12(i1,i2)%e(3)*p8k0+p8(3)*ctrip12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=ctrip12(i1,i2)%e(0)*p812(0)-ctrip12(i1,i2)%e(1)*p812(
     & 1)-ctrip12(i1,i2)%e(2)*p812(2)-ctrip12(i1,i2)%e(3)*p812(3
     & )
      cvqd=ctrip12(i1,i2)%e(0)*p8(0)-ctrip12(i1,i2)%e(1)*p8(1)-c
     & trip12(i1,i2)%e(2)*p8(2)-ctrip12(i1,i2)%e(3)*p8(3)
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p812k0*cvqd+p8k0*cvqu
      cauxb=-ctrip12(i1,i2)%ek0*p8(2)+p8k0*ctrip12(i1,i2)%e(2)
      r8_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      r8_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      r8_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      r8_gg(i1,i2)%b(2)=ccr*(-cauxb-ceps_2)
      end do
      end do
  
* quqd -- p=p7,q=p712
      quqd=p7(0)*p712(0)-p7(1)*p712(1)-p7(2)*p712(2)-p7(3)*p712(
     & 3)
      ccr=1.d0/(f712*s12)
      ccl=1.d0/(f712*s12)
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p7,qd=p712,v=ctrip12(i1,i2)%e,a=l7_gg(i1,i2)%a,c=l7_gg(i1,i2
* )%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p7(2)*p712(3)-p712(2)*p7(3))+p
     & 7k0*(ctrip12(i1,i2)%e(2)*p712(3)-p712(2)*ctrip12(i1,i2)%e
     & (3))-p712k0*(ctrip12(i1,i2)%e(2)*p7(3)-p7(2)*ctrip12(i1,i
     & 2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ctrip12(i1,i2)%e(3)*p7k0+p7(3)*ctrip12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=ctrip12(i1,i2)%e(0)*p7(0)-ctrip12(i1,i2)%e(1)*p7(1)-c
     & trip12(i1,i2)%e(2)*p7(2)-ctrip12(i1,i2)%e(3)*p7(3)
      cvqd=ctrip12(i1,i2)%e(0)*p712(0)-ctrip12(i1,i2)%e(1)*p712(
     & 1)-ctrip12(i1,i2)%e(2)*p712(2)-ctrip12(i1,i2)%e(3)*p712(3
     & )
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p7k0*cvqd+p712k0*cvqu
      cauxc=+ctrip12(i1,i2)%ek0*p7(2)-p7k0*ctrip12(i1,i2)%e(2)
      l7_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l7_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_gg(i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l7_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      cden=(-s7812+cmz2)
* quqd -- p=p712,q=p8
      quqd=p712(0)*p8(0)-p712(1)*p8(1)-p712(2)*p8(2)-p712(3)*p8(
     & 3)
      ccr=zcr(id8)/cden
      ccl=zcl(id8)/cden
* TR0 -- qu=p712,qd=p8,v=0,a=rz8712(0)%a,b=rz8712(0)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p712(2)*p8(3)+p8(2)*p712(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p712k0*p8(0)+p8k0*p712(0)
      rz8712(0)%a(1)=ccr*(auxa+ceps_0)
      rz8712(0)%a(2)=ccl*(auxa-ceps_0)
      rz8712(0)%b(1)=-ccl*(p8(2)+ceps_2)
      rz8712(0)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p712,qd=p8,v=1,a=rz8712(1)%a,b=rz8712(1)%b,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p712k0*p8(1)+p8k0*p712(1)
      rz8712(1)%a(1)=ccr*(auxa+ceps_0)
      rz8712(1)%a(2)=ccl*(auxa-ceps_0)
      rz8712(1)%b(1)=-ccl*(p8(2)+ceps_2)
      rz8712(1)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p712,qd=p8,v=2,a=rz8712(2)%a,b=rz8712(2)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p712k0*p8(3)+p8k0*p712(3)
      ceps_0=eps_0*cim
      auxa=p712k0*p8(2)+p8k0*p712(2)
      rz8712(2)%a(1)=ccr*(auxa+ceps_0)
      rz8712(2)%a(2)=ccl*(auxa-ceps_0)
      rz8712(2)%b(1)=-ccl*p8k0
      rz8712(2)%b(2)=ccr*p8k0
* TR0 -- qu=p712,qd=p8,v=3,a=rz8712(3)%a,b=rz8712(3)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p712k0*p8(2)-p8k0*p712(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p712k0*p8(3)+p8k0*p712(3)
      rz8712(3)%a(1)=ccr*(auxa+ceps_0)
      rz8712(3)%a(2)=ccl*(auxa-ceps_0)
      rz8712(3)%b(1)=-ccl*ceps_2
      rz8712(3)%b(2)=-ccr*ceps_2
  
      cden=(-s7812+cmz2)*f812
* quqd -- p=p7,q=p812
      quqd=p7(0)*p812(0)-p7(1)*p812(1)-p7(2)*p812(2)-p7(3)*p812(
     & 3)
      ccr=zcr(id7)/cden
      ccl=zcl(id7)/cden
* TL0 -- qu=p7,qd=p812,v=0,a=lz7812(0)%a,c=lz7812(0)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p7(2)*p812(3)+p812(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p812(0)+p812k0*p7(0)
      lz7812(0)%a(1)=ccr*(auxa+ceps_0)
      lz7812(0)%a(2)=ccl*(auxa-ceps_0)
      lz7812(0)%c(1)=ccr*(p7(2)+ceps_1)
      lz7812(0)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p812,v=1,a=lz7812(1)%a,c=lz7812(1)%c,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p7k0*p812(1)+p812k0*p7(1)
      lz7812(1)%a(1)=ccr*(auxa+ceps_0)
      lz7812(1)%a(2)=ccl*(auxa-ceps_0)
      lz7812(1)%c(1)=ccr*(p7(2)+ceps_1)
      lz7812(1)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p812,v=2,a=lz7812(2)%a,c=lz7812(2)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p7k0*p812(3)+p812k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p812(2)+p812k0*p7(2)
      lz7812(2)%a(1)=ccr*(auxa+ceps_0)
      lz7812(2)%a(2)=ccl*(auxa-ceps_0)
      lz7812(2)%c(1)=ccr*p7k0
      lz7812(2)%c(2)=-ccl*p7k0
* TL0 -- qu=p7,qd=p812,v=3,a=lz7812(3)%a,c=lz7812(3)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p7k0*p812(2)-p812k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p812(3)+p812k0*p7(3)
      lz7812(3)%a(1)=ccr*(auxa+ceps_0)
      lz7812(3)%a(2)=ccl*(auxa-ceps_0)
      lz7812(3)%c(1)=ccr*ceps_1
      lz7812(3)%c(2)=ccl*ceps_1
      cden=(-s7812)
* quqd -- p=p712,q=p8
      quqd=p712(0)*p8(0)-p712(1)*p8(1)-p712(2)*p8(2)-p712(3)*p8(
     & 3)
      ccr=fcr(id8)/cden
      ccl=fcl(id8)/cden
* TR0 -- qu=p712,qd=p8,v=0,a=rf8712(0)%a,b=rf8712(0)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p712(2)*p8(3)+p8(2)*p712(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p712k0*p8(0)+p8k0*p712(0)
      rf8712(0)%a(1)=ccr*(auxa+ceps_0)
      rf8712(0)%a(2)=ccl*(auxa-ceps_0)
      rf8712(0)%b(1)=-ccl*(p8(2)+ceps_2)
      rf8712(0)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p712,qd=p8,v=1,a=rf8712(1)%a,b=rf8712(1)%b,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p712k0*p8(1)+p8k0*p712(1)
      rf8712(1)%a(1)=ccr*(auxa+ceps_0)
      rf8712(1)%a(2)=ccl*(auxa-ceps_0)
      rf8712(1)%b(1)=-ccl*(p8(2)+ceps_2)
      rf8712(1)%b(2)=ccr*(p8(2)-ceps_2)
* TR0 -- qu=p712,qd=p8,v=2,a=rf8712(2)%a,b=rf8712(2)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p712k0*p8(3)+p8k0*p712(3)
      ceps_0=eps_0*cim
      auxa=p712k0*p8(2)+p8k0*p712(2)
      rf8712(2)%a(1)=ccr*(auxa+ceps_0)
      rf8712(2)%a(2)=ccl*(auxa-ceps_0)
      rf8712(2)%b(1)=-ccl*p8k0
      rf8712(2)%b(2)=ccr*p8k0
* TR0 -- qu=p712,qd=p8,v=3,a=rf8712(3)%a,b=rf8712(3)%b,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p712k0*p8(2)-p8k0*p712(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p712k0*p8(3)+p8k0*p712(3)
      rf8712(3)%a(1)=ccr*(auxa+ceps_0)
      rf8712(3)%a(2)=ccl*(auxa-ceps_0)
      rf8712(3)%b(1)=-ccl*ceps_2
      rf8712(3)%b(2)=-ccr*ceps_2
  
      cden=(-s7812)*f812
* quqd -- p=p7,q=p812
      quqd=p7(0)*p812(0)-p7(1)*p812(1)-p7(2)*p812(2)-p7(3)*p812(
     & 3)
      ccr=fcr(id7)/cden
      ccl=fcl(id7)/cden
* TL0 -- qu=p7,qd=p812,v=0,a=lf7812(0)%a,c=lf7812(0)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p7(2)*p812(3)+p812(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p812(0)+p812k0*p7(0)
      lf7812(0)%a(1)=ccr*(auxa+ceps_0)
      lf7812(0)%a(2)=ccl*(auxa-ceps_0)
      lf7812(0)%c(1)=ccr*(p7(2)+ceps_1)
      lf7812(0)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p812,v=1,a=lf7812(1)%a,c=lf7812(1)%c,cr=ccr,cl=ccl,nsu
* m=0
      auxa=-quqd+p7k0*p812(1)+p812k0*p7(1)
      lf7812(1)%a(1)=ccr*(auxa+ceps_0)
      lf7812(1)%a(2)=ccl*(auxa-ceps_0)
      lf7812(1)%c(1)=ccr*(p7(2)+ceps_1)
      lf7812(1)%c(2)=ccl*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p812,v=2,a=lf7812(2)%a,c=lf7812(2)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=-p7k0*p812(3)+p812k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p812(2)+p812k0*p7(2)
      lf7812(2)%a(1)=ccr*(auxa+ceps_0)
      lf7812(2)%a(2)=ccl*(auxa-ceps_0)
      lf7812(2)%c(1)=ccr*p7k0
      lf7812(2)%c(2)=-ccl*p7k0
* TL0 -- qu=p7,qd=p812,v=3,a=lf7812(3)%a,c=lf7812(3)%c,cr=ccr,cl=ccl,nsu
* m=0
      eps_0=p7k0*p812(2)-p812k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p812(3)+p812k0*p7(3)
      lf7812(3)%a(1)=ccr*(auxa+ceps_0)
      lf7812(3)%a(2)=ccl*(auxa-ceps_0)
      lf7812(3)%c(1)=ccr*ceps_1
      lf7812(3)%c(2)=ccl*ceps_1
  
* quqd -- p=p71,q=p712
      quqd=p71(0)*p712(0)-p71(1)*p712(1)-p71(2)*p712(2)-p71(3)*p
     & 712(3)
      ccr=1.d0/(f712)
      ccl=1.d0/(f712)
      do i2=1,2
* T0 -- qu=p71,qd=p712,v=ce2(i2)%e,a=u71_2(i2)%a,b=u71_2(i2)%b,c=u71_2(i
* 2)%c,d=u71_2(i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i2)%ek0*(p71(2)*p712(3)-p712(2)*p71(3))+p71k0*
     & (ce2(i2)%e(2)*p712(3)-p712(2)*ce2(i2)%e(3))-p712k0*(ce2(i
     & 2)%e(2)*p71(3)-p71(2)*ce2(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i2)%e(3)*p71k0+p71(3)*ce2(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i2)%e(3)*p712k0+p712(3)*ce2(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i2)%e(0)*p71(0)-ce2(i2)%e(1)*p71(1)-ce2(i2)%e(2)*
     & p71(2)-ce2(i2)%e(3)*p71(3)
      cvqd=ce2(i2)%e(0)*p712(0)-ce2(i2)%e(1)*p712(1)-ce2(i2)%e(2
     & )*p712(2)-ce2(i2)%e(3)*p712(3)
      cauxa=-ce2(i2)%ek0*quqd+p71k0*cvqd+p712k0*cvqu
      cauxb=-ce2(i2)%ek0*p712(2)+p712k0*ce2(i2)%e(2)
      cauxc=+ce2(i2)%ek0*p71(2)-p71k0*ce2(i2)%e(2)
      u71_2(i2)%a(1)=ccr*(cauxa+ceps_0)
      u71_2(i2)%a(2)=ccl*(cauxa-ceps_0)
      u71_2(i2)%b(1)=ccl*(cauxb-ceps_2)
      u71_2(i2)%b(2)=ccr*(-cauxb-ceps_2)
      u71_2(i2)%c(1)=ccr*(cauxc+ceps_1)
      u71_2(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u71_2(i2)%d(1)=ccl*ce2(i2)%ek0
      u71_2(i2)%d(2)=ccr*ce2(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=l7_12(i1,i2)%a,cc=l7_12(i1,i2)%c,a1=l7_1(i1)%a,c1=l7_1(i1)%
* c,a2=u71_2(i2)%a,b2=u71_2(i2)%b,c2=u71_2(i2)%c,d2=u71_2(i2)%d,prq=s71,
* nsum=0
      l7_12(i1,i2)%a(1)=l7_1(i1)%a(1)*u71_2(i2)%a(1)+l7_1(i1)%c(
     & 1)*s71*u71_2(i2)%b(2)
      l7_12(i1,i2)%c(1)=l7_1(i1)%a(1)*u71_2(i2)%c(1)+l7_1(i1)%c(
     & 1)*s71*u71_2(i2)%d(2)
      l7_12(i1,i2)%c(2)=l7_1(i1)%c(2)*s71*u71_2(i2)%d(1)+l7_1(i1
     & )%a(2)*u71_2(i2)%c(2)
      l7_12(i1,i2)%a(2)=l7_1(i1)%c(2)*s71*u71_2(i2)%b(1)+l7_1(i1
     & )%a(2)*u71_2(i2)%a(2)
      end do
      end do
  
* quqd -- p=p812,q=p82
      quqd=p812(0)*p82(0)-p812(1)*p82(1)-p812(2)*p82(2)-p812(3)*
     & p82(3)
      ccr=1.d0/(f82)
      ccl=1.d0/(f82)
      do i1=1,2
* T0 -- qu=p812,qd=p82,v=ce1(i1)%e,a=u82_1(i1)%a,b=u82_1(i1)%b,c=u82_1(i
* 1)%c,d=u82_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p812(2)*p82(3)-p82(2)*p812(3))+p812k0
     & *(ce1(i1)%e(2)*p82(3)-p82(2)*ce1(i1)%e(3))-p82k0*(ce1(i1)
     & %e(2)*p812(3)-p812(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p812k0+p812(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p82k0+p82(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p812(0)-ce1(i1)%e(1)*p812(1)-ce1(i1)%e(2
     & )*p812(2)-ce1(i1)%e(3)*p812(3)
      cvqd=ce1(i1)%e(0)*p82(0)-ce1(i1)%e(1)*p82(1)-ce1(i1)%e(2)*
     & p82(2)-ce1(i1)%e(3)*p82(3)
      cauxa=-ce1(i1)%ek0*quqd+p812k0*cvqd+p82k0*cvqu
      cauxb=-ce1(i1)%ek0*p82(2)+p82k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p812(2)-p812k0*ce1(i1)%e(2)
      u82_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u82_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u82_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u82_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u82_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u82_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u82_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u82_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0 -- aa=r8_12(i1,i2)%a,bb=r8_12(i1,i2)%b,a1=u82_1(i1)%a,b1=u82_1(i1
* )%b,c1=u82_1(i1)%c,d1=u82_1(i1)%d,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s82,
* nsum=0
      r8_12(i1,i2)%a(1)=u82_1(i1)%a(1)*r8_2(i2)%a(1)+u82_1(i1)%c
     & (1)*s82*r8_2(i2)%b(2)
      r8_12(i1,i2)%b(1)=u82_1(i1)%d(1)*s82*r8_2(i2)%b(1)+u82_1(i
     & 1)%b(1)*r8_2(i2)%a(2)
      r8_12(i1,i2)%b(2)=u82_1(i1)%b(2)*r8_2(i2)%a(1)+u82_1(i1)%d
     & (2)*s82*r8_2(i2)%b(2)
      r8_12(i1,i2)%a(2)=u82_1(i1)%c(2)*s82*r8_2(i2)%b(1)+u82_1(i
     & 1)%a(2)*r8_2(i2)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l7_12(i1,i2)%a(1)=l7_12(i1,i2)%a(1) +
     &   l7_gg(i1,i2)%a(1)
       l7_12(i1,i2)%a(2)=l7_12(i1,i2)%a(2) +
     &   l7_gg(i1,i2)%a(2)
       l7_12(i1,i2)%c(1)=l7_12(i1,i2)%c(1) +
     &   l7_gg(i1,i2)%c(1)
       l7_12(i1,i2)%c(2)=l7_12(i1,i2)%c(2) +
     &   l7_gg(i1,i2)%c(2)
  
       r8_12(i1,i2)%a(1)=r8_12(i1,i2)%a(1) +
     &   r8_gg(i1,i2)%a(1)
       r8_12(i1,i2)%a(2)=r8_12(i1,i2)%a(2) +
     &   r8_gg(i1,i2)%a(2)
       r8_12(i1,i2)%b(1)=r8_12(i1,i2)%b(1) +
     &   r8_gg(i1,i2)%b(1)
       r8_12(i1,i2)%b(2)=r8_12(i1,i2)%b(2) +
     &   r8_gg(i1,i2)%b(2)
  
      enddo
      enddo
  
  
      cden=(-s7812+cmz2)*f82
* quqd -- p=p71,q=p82
      quqd=p71(0)*p82(0)-p71(1)*p82(1)-p71(2)*p82(2)-p71(3)*p82(
     & 3)
      ccr=zcr(id7)/cden
      ccl=zcl(id7)/cden
* T0 -- qu=p71,qd=p82,v=0,a=uz71(0)%a,b=uz71(0)%b,c=uz71(0)%c,d=uz71(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p71(2)*p82(3)+p82(2)*p71(3)
      ceps_0=eps_0*cim
      ceps_1=p71(3)*cim
      ceps_2=p82(3)*cim
      auxa=-quqd+p71k0*p82(0)+p82k0*p71(0)
      uz71(0)%a(1)=ccr*(auxa+ceps_0)
      uz71(0)%a(2)=ccl*(auxa-ceps_0)
      uz71(0)%b(1)=-ccl*(p82(2)+ceps_2)
      uz71(0)%b(2)=ccr*(p82(2)-ceps_2)
      uz71(0)%c(1)=ccr*(p71(2)+ceps_1)
      uz71(0)%c(2)=ccl*(-p71(2)+ceps_1)
      uz71(0)%d(1)=ccl
      uz71(0)%d(2)=ccr
* T0 -- qu=p71,qd=p82,v=1,a=uz71(1)%a,b=uz71(1)%b,c=uz71(1)%c,d=uz71(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p71k0*p82(1)+p82k0*p71(1)
      uz71(1)%a(1)=ccr*(auxa+ceps_0)
      uz71(1)%a(2)=ccl*(auxa-ceps_0)
      uz71(1)%b(1)=-ccl*(p82(2)+ceps_2)
      uz71(1)%b(2)=ccr*(p82(2)-ceps_2)
      uz71(1)%c(1)=ccr*(p71(2)+ceps_1)
      uz71(1)%c(2)=ccl*(-p71(2)+ceps_1)
      uz71(1)%d(1)=ccl
      uz71(1)%d(2)=ccr
* T0 -- qu=p71,qd=p82,v=2,a=uz71(2)%a,b=uz71(2)%b,c=uz71(2)%c,d=uz71(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p71k0*p82(3)+p82k0*p71(3)
      ceps_0=eps_0*cim
      auxa=p71k0*p82(2)+p82k0*p71(2)
      uz71(2)%a(1)=ccr*(auxa+ceps_0)
      uz71(2)%a(2)=ccl*(auxa-ceps_0)
      uz71(2)%b(1)=-ccl*p82k0
      uz71(2)%b(2)=ccr*p82k0
      uz71(2)%c(1)=ccr*p71k0
      uz71(2)%c(2)=-ccl*p71k0
* T0 -- qu=p71,qd=p82,v=3,a=uz71(3)%a,b=uz71(3)%b,c=uz71(3)%c,d=uz71(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p71k0*p82(2)-p82k0*p71(2)
      ceps_0=eps_0*cim
      ceps_1=p71k0*cim
      ceps_2=p82k0*cim
      auxa=+p71k0*p82(3)+p82k0*p71(3)
      uz71(3)%a(1)=ccr*(auxa+ceps_0)
      uz71(3)%a(2)=ccl*(auxa-ceps_0)
      uz71(3)%b(1)=-ccl*ceps_2
      uz71(3)%b(2)=-ccr*ceps_2
      uz71(3)%c(1)=ccr*ceps_1
      uz71(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l7_1(i1)%a,c1=l7_
* 1(i1)%c,a2=uz71(mu)%a,b2=uz71(mu)%b,c2=uz71(mu)%c,d2=uz71(mu)%d,prq=s7
* 1,nsum=0
      laux_imu(i1,mu)%a(1)=l7_1(i1)%a(1)*uz71(mu)%a(1)+l7_1(i1)%
     & c(1)*s71*uz71(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l7_1(i1)%a(1)*uz71(mu)%c(1)+l7_1(i1)%
     & c(1)*s71*uz71(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l7_1(i1)%c(2)*s71*uz71(mu)%d(1)+l7_1(
     & i1)%a(2)*uz71(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l7_1(i1)%c(2)*s71*uz71(mu)%b(1)+l7_1(
     & i1)%a(2)*uz71(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz7128(&,i1,i2)%e(mu),a1=laux_imu(i1,mu)%a,c1=laux_imu(i1,
* mu)%c,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s82,bef=,aft=
      cz7128(1,i1,i2)%e(mu)=(laux_imu(i1,mu)%a(1)*r8_2(i2)%a(1)+
     & laux_imu(i1,mu)%c(1)*s82*r8_2(i2)%b(2))
      cz7128(2,i1,i2)%e(mu)=(laux_imu(i1,mu)%c(2)*s82*r8_2(i2)%b
     & (1)+laux_imu(i1,mu)%a(2)*r8_2(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz7128(&,i1,i2)%e(mu),a1=l7_12(i1,i2)%a,c1=l7_12(i1,i2)%c,
* a2=rz8712(mu)%a,b2=rz8712(mu)%b,prq=s712,bef=cz7128(&,i1,i2)%e(mu)+,af
* t=
      cz7128(1,i1,i2)%e(mu)=cz7128(1,i1,i2)%e(mu)+(l7_12(i1,i2)%
     & a(1)*rz8712(mu)%a(1)+l7_12(i1,i2)%c(1)*s712*rz8712(mu)%b(
     & 2))
      cz7128(2,i1,i2)%e(mu)=cz7128(2,i1,i2)%e(mu)+(l7_12(i1,i2)%
     & c(2)*s712*rz8712(mu)%b(1)+l7_12(i1,i2)%a(2)*rz8712(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz7128(&,i1,i2)%e(mu),a1=lz7812(mu)%a,c1=lz7812(mu)%c,a2=r
* 8_12(i1,i2)%a,b2=r8_12(i1,i2)%b,prq=s812,bef=cz7128(&,i1,i2)%e(mu)+,af
* t=
      cz7128(1,i1,i2)%e(mu)=cz7128(1,i1,i2)%e(mu)+(lz7812(mu)%a(
     & 1)*r8_12(i1,i2)%a(1)+lz7812(mu)%c(1)*s812*r8_12(i1,i2)%b(
     & 2))
      cz7128(2,i1,i2)%e(mu)=cz7128(2,i1,i2)%e(mu)+(lz7812(mu)%c(
     & 2)*s812*r8_12(i1,i2)%b(1)+lz7812(mu)%a(2)*r8_12(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cz7128(i3,i1,i2)%e
      cz7128(i3,i1,i2)%ek0=cz7128(i3,i1,i2)%e(0)-cz7128(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
      cden=(-s7812)*f82
* quqd -- p=p71,q=p82
      quqd=p71(0)*p82(0)-p71(1)*p82(1)-p71(2)*p82(2)-p71(3)*p82(
     & 3)
      ccr=fcr(id7)/cden
      ccl=fcl(id7)/cden
* T0 -- qu=p71,qd=p82,v=0,a=uf71(0)%a,b=uf71(0)%b,c=uf71(0)%c,d=uf71(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p71(2)*p82(3)+p82(2)*p71(3)
      ceps_0=eps_0*cim
      ceps_1=p71(3)*cim
      ceps_2=p82(3)*cim
      auxa=-quqd+p71k0*p82(0)+p82k0*p71(0)
      uf71(0)%a(1)=ccr*(auxa+ceps_0)
      uf71(0)%a(2)=ccl*(auxa-ceps_0)
      uf71(0)%b(1)=-ccl*(p82(2)+ceps_2)
      uf71(0)%b(2)=ccr*(p82(2)-ceps_2)
      uf71(0)%c(1)=ccr*(p71(2)+ceps_1)
      uf71(0)%c(2)=ccl*(-p71(2)+ceps_1)
      uf71(0)%d(1)=ccl
      uf71(0)%d(2)=ccr
* T0 -- qu=p71,qd=p82,v=1,a=uf71(1)%a,b=uf71(1)%b,c=uf71(1)%c,d=uf71(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p71k0*p82(1)+p82k0*p71(1)
      uf71(1)%a(1)=ccr*(auxa+ceps_0)
      uf71(1)%a(2)=ccl*(auxa-ceps_0)
      uf71(1)%b(1)=-ccl*(p82(2)+ceps_2)
      uf71(1)%b(2)=ccr*(p82(2)-ceps_2)
      uf71(1)%c(1)=ccr*(p71(2)+ceps_1)
      uf71(1)%c(2)=ccl*(-p71(2)+ceps_1)
      uf71(1)%d(1)=ccl
      uf71(1)%d(2)=ccr
* T0 -- qu=p71,qd=p82,v=2,a=uf71(2)%a,b=uf71(2)%b,c=uf71(2)%c,d=uf71(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p71k0*p82(3)+p82k0*p71(3)
      ceps_0=eps_0*cim
      auxa=p71k0*p82(2)+p82k0*p71(2)
      uf71(2)%a(1)=ccr*(auxa+ceps_0)
      uf71(2)%a(2)=ccl*(auxa-ceps_0)
      uf71(2)%b(1)=-ccl*p82k0
      uf71(2)%b(2)=ccr*p82k0
      uf71(2)%c(1)=ccr*p71k0
      uf71(2)%c(2)=-ccl*p71k0
* T0 -- qu=p71,qd=p82,v=3,a=uf71(3)%a,b=uf71(3)%b,c=uf71(3)%c,d=uf71(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p71k0*p82(2)-p82k0*p71(2)
      ceps_0=eps_0*cim
      ceps_1=p71k0*cim
      ceps_2=p82k0*cim
      auxa=+p71k0*p82(3)+p82k0*p71(3)
      uf71(3)%a(1)=ccr*(auxa+ceps_0)
      uf71(3)%a(2)=ccl*(auxa-ceps_0)
      uf71(3)%b(1)=-ccl*ceps_2
      uf71(3)%b(2)=-ccr*ceps_2
      uf71(3)%c(1)=ccr*ceps_1
      uf71(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l7_1(i1)%a,c1=l7_
* 1(i1)%c,a2=uf71(mu)%a,b2=uf71(mu)%b,c2=uf71(mu)%c,d2=uf71(mu)%d,prq=s7
* 1,nsum=0
      laux_imu(i1,mu)%a(1)=l7_1(i1)%a(1)*uf71(mu)%a(1)+l7_1(i1)%
     & c(1)*s71*uf71(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l7_1(i1)%a(1)*uf71(mu)%c(1)+l7_1(i1)%
     & c(1)*s71*uf71(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l7_1(i1)%c(2)*s71*uf71(mu)%d(1)+l7_1(
     & i1)%a(2)*uf71(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l7_1(i1)%c(2)*s71*uf71(mu)%b(1)+l7_1(
     & i1)%a(2)*uf71(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf7128(&,i1,i2)%e(mu),a1=laux_imu(i1,mu)%a,c1=laux_imu(i1,
* mu)%c,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s82,bef=,aft=
      cf7128(1,i1,i2)%e(mu)=(laux_imu(i1,mu)%a(1)*r8_2(i2)%a(1)+
     & laux_imu(i1,mu)%c(1)*s82*r8_2(i2)%b(2))
      cf7128(2,i1,i2)%e(mu)=(laux_imu(i1,mu)%c(2)*s82*r8_2(i2)%b
     & (1)+laux_imu(i1,mu)%a(2)*r8_2(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf7128(&,i1,i2)%e(mu),a1=l7_12(i1,i2)%a,c1=l7_12(i1,i2)%c,
* a2=rf8712(mu)%a,b2=rf8712(mu)%b,prq=s712,bef=cf7128(&,i1,i2)%e(mu)+,af
* t=
      cf7128(1,i1,i2)%e(mu)=cf7128(1,i1,i2)%e(mu)+(l7_12(i1,i2)%
     & a(1)*rf8712(mu)%a(1)+l7_12(i1,i2)%c(1)*s712*rf8712(mu)%b(
     & 2))
      cf7128(2,i1,i2)%e(mu)=cf7128(2,i1,i2)%e(mu)+(l7_12(i1,i2)%
     & c(2)*s712*rf8712(mu)%b(1)+l7_12(i1,i2)%a(2)*rf8712(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf7128(&,i1,i2)%e(mu),a1=lf7812(mu)%a,c1=lf7812(mu)%c,a2=r
* 8_12(i1,i2)%a,b2=r8_12(i1,i2)%b,prq=s812,bef=cf7128(&,i1,i2)%e(mu)+,af
* t=
      cf7128(1,i1,i2)%e(mu)=cf7128(1,i1,i2)%e(mu)+(lf7812(mu)%a(
     & 1)*r8_12(i1,i2)%a(1)+lf7812(mu)%c(1)*s812*r8_12(i1,i2)%b(
     & 2))
      cf7128(2,i1,i2)%e(mu)=cf7128(2,i1,i2)%e(mu)+(lf7812(mu)%c(
     & 2)*s812*r8_12(i1,i2)%b(1)+lf7812(mu)%a(2)*r8_12(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cf7128(i3,i1,i2)%e
      cf7128(i3,i1,i2)%ek0=cf7128(i3,i1,i2)%e(0)-cf7128(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
* quqd -- p=p72,q=p712
      quqd=p72(0)*p712(0)-p72(1)*p712(1)-p72(2)*p712(2)-p72(3)*p
     & 712(3)
      ccr=1.d0/(f712)
      ccl=1.d0/(f712)
      do i2=1,2
* T0 -- qu=p72,qd=p712,v=ce1(i2)%e,a=u72_1(i2)%a,b=u72_1(i2)%b,c=u72_1(i
* 2)%c,d=u72_1(i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i2)%ek0*(p72(2)*p712(3)-p712(2)*p72(3))+p72k0*
     & (ce1(i2)%e(2)*p712(3)-p712(2)*ce1(i2)%e(3))-p712k0*(ce1(i
     & 2)%e(2)*p72(3)-p72(2)*ce1(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i2)%e(3)*p72k0+p72(3)*ce1(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i2)%e(3)*p712k0+p712(3)*ce1(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i2)%e(0)*p72(0)-ce1(i2)%e(1)*p72(1)-ce1(i2)%e(2)*
     & p72(2)-ce1(i2)%e(3)*p72(3)
      cvqd=ce1(i2)%e(0)*p712(0)-ce1(i2)%e(1)*p712(1)-ce1(i2)%e(2
     & )*p712(2)-ce1(i2)%e(3)*p712(3)
      cauxa=-ce1(i2)%ek0*quqd+p72k0*cvqd+p712k0*cvqu
      cauxb=-ce1(i2)%ek0*p712(2)+p712k0*ce1(i2)%e(2)
      cauxc=+ce1(i2)%ek0*p72(2)-p72k0*ce1(i2)%e(2)
      u72_1(i2)%a(1)=ccr*(cauxa+ceps_0)
      u72_1(i2)%a(2)=ccl*(cauxa-ceps_0)
      u72_1(i2)%b(1)=ccl*(cauxb-ceps_2)
      u72_1(i2)%b(2)=ccr*(-cauxb-ceps_2)
      u72_1(i2)%c(1)=ccr*(cauxc+ceps_1)
      u72_1(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u72_1(i2)%d(1)=ccl*ce1(i2)%ek0
      u72_1(i2)%d(2)=ccr*ce1(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=l7_21(i1,i2)%a,cc=l7_21(i1,i2)%c,a1=l7_2(i2)%a,c1=l7_2(i2)%
* c,a2=u72_1(i1)%a,b2=u72_1(i1)%b,c2=u72_1(i1)%c,d2=u72_1(i1)%d,prq=s72,
* nsum=0
      l7_21(i1,i2)%a(1)=l7_2(i2)%a(1)*u72_1(i1)%a(1)+l7_2(i2)%c(
     & 1)*s72*u72_1(i1)%b(2)
      l7_21(i1,i2)%c(1)=l7_2(i2)%a(1)*u72_1(i1)%c(1)+l7_2(i2)%c(
     & 1)*s72*u72_1(i1)%d(2)
      l7_21(i1,i2)%c(2)=l7_2(i2)%c(2)*s72*u72_1(i1)%d(1)+l7_2(i2
     & )%a(2)*u72_1(i1)%c(2)
      l7_21(i1,i2)%a(2)=l7_2(i2)%c(2)*s72*u72_1(i1)%b(1)+l7_2(i2
     & )%a(2)*u72_1(i1)%a(2)
      end do
      end do
  
* quqd -- p=p812,q=p81
      quqd=p812(0)*p81(0)-p812(1)*p81(1)-p812(2)*p81(2)-p812(3)*
     & p81(3)
      ccr=1.d0/(f81)
      ccl=1.d0/(f81)
      do i1=1,2
* T0 -- qu=p812,qd=p81,v=ce2(i1)%e,a=u81_2(i1)%a,b=u81_2(i1)%b,c=u81_2(i
* 1)%c,d=u81_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p812(2)*p81(3)-p81(2)*p812(3))+p812k0
     & *(ce2(i1)%e(2)*p81(3)-p81(2)*ce2(i1)%e(3))-p81k0*(ce2(i1)
     & %e(2)*p812(3)-p812(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p812k0+p812(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p81k0+p81(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p812(0)-ce2(i1)%e(1)*p812(1)-ce2(i1)%e(2
     & )*p812(2)-ce2(i1)%e(3)*p812(3)
      cvqd=ce2(i1)%e(0)*p81(0)-ce2(i1)%e(1)*p81(1)-ce2(i1)%e(2)*
     & p81(2)-ce2(i1)%e(3)*p81(3)
      cauxa=-ce2(i1)%ek0*quqd+p812k0*cvqd+p81k0*cvqu
      cauxb=-ce2(i1)%ek0*p81(2)+p81k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p812(2)-p812k0*ce2(i1)%e(2)
      u81_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u81_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u81_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u81_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u81_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u81_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u81_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u81_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0 -- aa=r8_21(i1,i2)%a,bb=r8_21(i1,i2)%b,a1=u81_2(i2)%a,b1=u81_2(i2
* )%b,c1=u81_2(i2)%c,d1=u81_2(i2)%d,a2=r8_1(i1)%a,b2=r8_1(i1)%b,prq=s81,
* nsum=0
      r8_21(i1,i2)%a(1)=u81_2(i2)%a(1)*r8_1(i1)%a(1)+u81_2(i2)%c
     & (1)*s81*r8_1(i1)%b(2)
      r8_21(i1,i2)%b(1)=u81_2(i2)%d(1)*s81*r8_1(i1)%b(1)+u81_2(i
     & 2)%b(1)*r8_1(i1)%a(2)
      r8_21(i1,i2)%b(2)=u81_2(i2)%b(2)*r8_1(i1)%a(1)+u81_2(i2)%d
     & (2)*s81*r8_1(i1)%b(2)
      r8_21(i1,i2)%a(2)=u81_2(i2)%c(2)*s81*r8_1(i1)%b(1)+u81_2(i
     & 2)%a(2)*r8_1(i1)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l7_21(i1,i2)%a(1)=l7_21(i1,i2)%a(1) -
     &   l7_gg(i1,i2)%a(1)
       l7_21(i1,i2)%a(2)=l7_21(i1,i2)%a(2) -
     &   l7_gg(i1,i2)%a(2)
       l7_21(i1,i2)%c(1)=l7_21(i1,i2)%c(1) -
     &   l7_gg(i1,i2)%c(1)
       l7_21(i1,i2)%c(2)=l7_21(i1,i2)%c(2) -
     &   l7_gg(i1,i2)%c(2)
  
       r8_21(i1,i2)%a(1)=r8_21(i1,i2)%a(1) -
     &   r8_gg(i1,i2)%a(1)
       r8_21(i1,i2)%a(2)=r8_21(i1,i2)%a(2) -
     &   r8_gg(i1,i2)%a(2)
       r8_21(i1,i2)%b(1)=r8_21(i1,i2)%b(1) -
     &   r8_gg(i1,i2)%b(1)
       r8_21(i1,i2)%b(2)=r8_21(i1,i2)%b(2) -
     &   r8_gg(i1,i2)%b(2)
  
      enddo
      enddo
  
  
      cden=(-s7812+cmz2)*f81
* quqd -- p=p72,q=p81
      quqd=p72(0)*p81(0)-p72(1)*p81(1)-p72(2)*p81(2)-p72(3)*p81(
     & 3)
      ccr=zcr(id7)/cden
      ccl=zcl(id7)/cden
* T0 -- qu=p72,qd=p81,v=0,a=uz72(0)%a,b=uz72(0)%b,c=uz72(0)%c,d=uz72(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p72(2)*p81(3)+p81(2)*p72(3)
      ceps_0=eps_0*cim
      ceps_1=p72(3)*cim
      ceps_2=p81(3)*cim
      auxa=-quqd+p72k0*p81(0)+p81k0*p72(0)
      uz72(0)%a(1)=ccr*(auxa+ceps_0)
      uz72(0)%a(2)=ccl*(auxa-ceps_0)
      uz72(0)%b(1)=-ccl*(p81(2)+ceps_2)
      uz72(0)%b(2)=ccr*(p81(2)-ceps_2)
      uz72(0)%c(1)=ccr*(p72(2)+ceps_1)
      uz72(0)%c(2)=ccl*(-p72(2)+ceps_1)
      uz72(0)%d(1)=ccl
      uz72(0)%d(2)=ccr
* T0 -- qu=p72,qd=p81,v=1,a=uz72(1)%a,b=uz72(1)%b,c=uz72(1)%c,d=uz72(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p72k0*p81(1)+p81k0*p72(1)
      uz72(1)%a(1)=ccr*(auxa+ceps_0)
      uz72(1)%a(2)=ccl*(auxa-ceps_0)
      uz72(1)%b(1)=-ccl*(p81(2)+ceps_2)
      uz72(1)%b(2)=ccr*(p81(2)-ceps_2)
      uz72(1)%c(1)=ccr*(p72(2)+ceps_1)
      uz72(1)%c(2)=ccl*(-p72(2)+ceps_1)
      uz72(1)%d(1)=ccl
      uz72(1)%d(2)=ccr
* T0 -- qu=p72,qd=p81,v=2,a=uz72(2)%a,b=uz72(2)%b,c=uz72(2)%c,d=uz72(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p72k0*p81(3)+p81k0*p72(3)
      ceps_0=eps_0*cim
      auxa=p72k0*p81(2)+p81k0*p72(2)
      uz72(2)%a(1)=ccr*(auxa+ceps_0)
      uz72(2)%a(2)=ccl*(auxa-ceps_0)
      uz72(2)%b(1)=-ccl*p81k0
      uz72(2)%b(2)=ccr*p81k0
      uz72(2)%c(1)=ccr*p72k0
      uz72(2)%c(2)=-ccl*p72k0
* T0 -- qu=p72,qd=p81,v=3,a=uz72(3)%a,b=uz72(3)%b,c=uz72(3)%c,d=uz72(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p72k0*p81(2)-p81k0*p72(2)
      ceps_0=eps_0*cim
      ceps_1=p72k0*cim
      ceps_2=p81k0*cim
      auxa=+p72k0*p81(3)+p81k0*p72(3)
      uz72(3)%a(1)=ccr*(auxa+ceps_0)
      uz72(3)%a(2)=ccl*(auxa-ceps_0)
      uz72(3)%b(1)=-ccl*ceps_2
      uz72(3)%b(2)=-ccr*ceps_2
      uz72(3)%c(1)=ccr*ceps_1
      uz72(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l7_2(i1)%a,c1=l7_
* 2(i1)%c,a2=uz72(mu)%a,b2=uz72(mu)%b,c2=uz72(mu)%c,d2=uz72(mu)%d,prq=s7
* 2,nsum=0
      laux_imu(i1,mu)%a(1)=l7_2(i1)%a(1)*uz72(mu)%a(1)+l7_2(i1)%
     & c(1)*s72*uz72(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l7_2(i1)%a(1)*uz72(mu)%c(1)+l7_2(i1)%
     & c(1)*s72*uz72(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l7_2(i1)%c(2)*s72*uz72(mu)%d(1)+l7_2(
     & i1)%a(2)*uz72(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l7_2(i1)%c(2)*s72*uz72(mu)%b(1)+l7_2(
     & i1)%a(2)*uz72(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz7218(&,i1,i2)%e(mu),a1=laux_imu(i2,mu)%a,c1=laux_imu(i2,
* mu)%c,a2=r8_1(i1)%a,b2=r8_1(i1)%b,prq=s81,bef=,aft=
      cz7218(1,i1,i2)%e(mu)=(laux_imu(i2,mu)%a(1)*r8_1(i1)%a(1)+
     & laux_imu(i2,mu)%c(1)*s81*r8_1(i1)%b(2))
      cz7218(2,i1,i2)%e(mu)=(laux_imu(i2,mu)%c(2)*s81*r8_1(i1)%b
     & (1)+laux_imu(i2,mu)%a(2)*r8_1(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz7218(&,i1,i2)%e(mu),a1=l7_21(i1,i2)%a,c1=l7_21(i1,i2)%c,
* a2=rz8712(mu)%a,b2=rz8712(mu)%b,prq=s712,bef=cz7218(&,i1,i2)%e(mu)+,af
* t=
      cz7218(1,i1,i2)%e(mu)=cz7218(1,i1,i2)%e(mu)+(l7_21(i1,i2)%
     & a(1)*rz8712(mu)%a(1)+l7_21(i1,i2)%c(1)*s712*rz8712(mu)%b(
     & 2))
      cz7218(2,i1,i2)%e(mu)=cz7218(2,i1,i2)%e(mu)+(l7_21(i1,i2)%
     & c(2)*s712*rz8712(mu)%b(1)+l7_21(i1,i2)%a(2)*rz8712(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cz7218(&,i1,i2)%e(mu),a1=lz7812(mu)%a,c1=lz7812(mu)%c,a2=r
* 8_21(i1,i2)%a,b2=r8_21(i1,i2)%b,prq=s812,bef=cz7218(&,i1,i2)%e(mu)+,af
* t=
      cz7218(1,i1,i2)%e(mu)=cz7218(1,i1,i2)%e(mu)+(lz7812(mu)%a(
     & 1)*r8_21(i1,i2)%a(1)+lz7812(mu)%c(1)*s812*r8_21(i1,i2)%b(
     & 2))
      cz7218(2,i1,i2)%e(mu)=cz7218(2,i1,i2)%e(mu)+(lz7812(mu)%c(
     & 2)*s812*r8_21(i1,i2)%b(1)+lz7812(mu)%a(2)*r8_21(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cz7218(i3,i1,i2)%e
      cz7218(i3,i1,i2)%ek0=cz7218(i3,i1,i2)%e(0)-cz7218(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
      cden=(-s7812)*f81
* quqd -- p=p72,q=p81
      quqd=p72(0)*p81(0)-p72(1)*p81(1)-p72(2)*p81(2)-p72(3)*p81(
     & 3)
      ccr=fcr(id7)/cden
      ccl=fcl(id7)/cden
* T0 -- qu=p72,qd=p81,v=0,a=uf72(0)%a,b=uf72(0)%b,c=uf72(0)%c,d=uf72(0)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p72(2)*p81(3)+p81(2)*p72(3)
      ceps_0=eps_0*cim
      ceps_1=p72(3)*cim
      ceps_2=p81(3)*cim
      auxa=-quqd+p72k0*p81(0)+p81k0*p72(0)
      uf72(0)%a(1)=ccr*(auxa+ceps_0)
      uf72(0)%a(2)=ccl*(auxa-ceps_0)
      uf72(0)%b(1)=-ccl*(p81(2)+ceps_2)
      uf72(0)%b(2)=ccr*(p81(2)-ceps_2)
      uf72(0)%c(1)=ccr*(p72(2)+ceps_1)
      uf72(0)%c(2)=ccl*(-p72(2)+ceps_1)
      uf72(0)%d(1)=ccl
      uf72(0)%d(2)=ccr
* T0 -- qu=p72,qd=p81,v=1,a=uf72(1)%a,b=uf72(1)%b,c=uf72(1)%c,d=uf72(1)%
* d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p72k0*p81(1)+p81k0*p72(1)
      uf72(1)%a(1)=ccr*(auxa+ceps_0)
      uf72(1)%a(2)=ccl*(auxa-ceps_0)
      uf72(1)%b(1)=-ccl*(p81(2)+ceps_2)
      uf72(1)%b(2)=ccr*(p81(2)-ceps_2)
      uf72(1)%c(1)=ccr*(p72(2)+ceps_1)
      uf72(1)%c(2)=ccl*(-p72(2)+ceps_1)
      uf72(1)%d(1)=ccl
      uf72(1)%d(2)=ccr
* T0 -- qu=p72,qd=p81,v=2,a=uf72(2)%a,b=uf72(2)%b,c=uf72(2)%c,d=uf72(2)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p72k0*p81(3)+p81k0*p72(3)
      ceps_0=eps_0*cim
      auxa=p72k0*p81(2)+p81k0*p72(2)
      uf72(2)%a(1)=ccr*(auxa+ceps_0)
      uf72(2)%a(2)=ccl*(auxa-ceps_0)
      uf72(2)%b(1)=-ccl*p81k0
      uf72(2)%b(2)=ccr*p81k0
      uf72(2)%c(1)=ccr*p72k0
      uf72(2)%c(2)=-ccl*p72k0
* T0 -- qu=p72,qd=p81,v=3,a=uf72(3)%a,b=uf72(3)%b,c=uf72(3)%c,d=uf72(3)%
* d,cr=ccr,cl=ccl,nsum=0
      eps_0=p72k0*p81(2)-p81k0*p72(2)
      ceps_0=eps_0*cim
      ceps_1=p72k0*cim
      ceps_2=p81k0*cim
      auxa=+p72k0*p81(3)+p81k0*p72(3)
      uf72(3)%a(1)=ccr*(auxa+ceps_0)
      uf72(3)%a(2)=ccl*(auxa-ceps_0)
      uf72(3)%b(1)=-ccl*ceps_2
      uf72(3)%b(2)=-ccr*ceps_2
      uf72(3)%c(1)=ccr*ceps_1
      uf72(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0 -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l7_2(i1)%a,c1=l7_
* 2(i1)%c,a2=uf72(mu)%a,b2=uf72(mu)%b,c2=uf72(mu)%c,d2=uf72(mu)%d,prq=s7
* 2,nsum=0
      laux_imu(i1,mu)%a(1)=l7_2(i1)%a(1)*uf72(mu)%a(1)+l7_2(i1)%
     & c(1)*s72*uf72(mu)%b(2)
      laux_imu(i1,mu)%c(1)=l7_2(i1)%a(1)*uf72(mu)%c(1)+l7_2(i1)%
     & c(1)*s72*uf72(mu)%d(2)
      laux_imu(i1,mu)%c(2)=l7_2(i1)%c(2)*s72*uf72(mu)%d(1)+l7_2(
     & i1)%a(2)*uf72(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l7_2(i1)%c(2)*s72*uf72(mu)%b(1)+l7_2(
     & i1)%a(2)*uf72(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf7218(&,i1,i2)%e(mu),a1=laux_imu(i2,mu)%a,c1=laux_imu(i2,
* mu)%c,a2=r8_1(i1)%a,b2=r8_1(i1)%b,prq=s81,bef=,aft=
      cf7218(1,i1,i2)%e(mu)=(laux_imu(i2,mu)%a(1)*r8_1(i1)%a(1)+
     & laux_imu(i2,mu)%c(1)*s81*r8_1(i1)%b(2))
      cf7218(2,i1,i2)%e(mu)=(laux_imu(i2,mu)%c(2)*s81*r8_1(i1)%b
     & (1)+laux_imu(i2,mu)%a(2)*r8_1(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf7218(&,i1,i2)%e(mu),a1=l7_21(i1,i2)%a,c1=l7_21(i1,i2)%c,
* a2=rf8712(mu)%a,b2=rf8712(mu)%b,prq=s712,bef=cf7218(&,i1,i2)%e(mu)+,af
* t=
      cf7218(1,i1,i2)%e(mu)=cf7218(1,i1,i2)%e(mu)+(l7_21(i1,i2)%
     & a(1)*rf8712(mu)%a(1)+l7_21(i1,i2)%c(1)*s712*rf8712(mu)%b(
     & 2))
      cf7218(2,i1,i2)%e(mu)=cf7218(2,i1,i2)%e(mu)+(l7_21(i1,i2)%
     & c(2)*s712*rf8712(mu)%b(1)+l7_21(i1,i2)%a(2)*rf8712(mu)%a(
     & 2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0 -- aa=cf7218(&,i1,i2)%e(mu),a1=lf7812(mu)%a,c1=lf7812(mu)%c,a2=r
* 8_21(i1,i2)%a,b2=r8_21(i1,i2)%b,prq=s812,bef=cf7218(&,i1,i2)%e(mu)+,af
* t=
      cf7218(1,i1,i2)%e(mu)=cf7218(1,i1,i2)%e(mu)+(lf7812(mu)%a(
     & 1)*r8_21(i1,i2)%a(1)+lf7812(mu)%c(1)*s812*r8_21(i1,i2)%b(
     & 2))
      cf7218(2,i1,i2)%e(mu)=cf7218(2,i1,i2)%e(mu)+(lf7812(mu)%c(
     & 2)*s812*r8_21(i1,i2)%b(1)+lf7812(mu)%a(2)*r8_21(i1,i2)%a(
     & 2))
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* pk0 -- p=cf7218(i3,i1,i2)%e
      cf7218(i3,i1,i2)%ek0=cf7218(i3,i1,i2)%e(0)-cf7218(i3,i1,i2
     & )%e(1)
      end do
      end do
      end do
  
      endif !id3 = quark
  
* -> lines without gluon                                                
  
       do m=0,3
        p356(m)=p3(m)+p5(m)+p6(m)
       enddo
* pk0 -- p=p356
      p356k0=p356(0)-p356(1)
* p.q -- p.q=s356,p=p356,q=p356,bef=,aft=
      s356=(p356(0)*p356(0)-p356(1)*p356(1)-p356(2)*p356(2)-p356
     & (3)*p356(3))
      f356=s356*p356k0
  
       do m=0,3
        p478(m)= -p4(m)-p7(m)-p8(m)
       enddo
* pk0 -- p=p478
      p478k0=p478(0)-p478(1)
* p.q -- p.q=s478,p=p478,q=p478,bef=,aft=
      s478=(p478(0)*p478(0)-p478(1)*p478(1)-p478(2)*p478(2)-p478
     & (3)*p478(3))
      f478=s478*p478k0
  
      if (ilept(id3).ne.1.or.ilept(id7).ne.1) then
      quqd=(s356-s56)/2d0
      ccr=zcr(id3)/(f356)
      ccl=zcl(id3)/(f356)
      do i5=1,2
* TL0 -- qu=p3,qd=p356,v=cz56(i5)%e,a=l3_56(i5)%a,c=l3_56(i5)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p3(2)*p356(3)-p356(2)*p3(3))+p3k0*(c
     & z56(i5)%e(2)*p356(3)-p356(2)*cz56(i5)%e(3))-p356k0*(cz56(
     & i5)%e(2)*p3(3)-p3(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p3k0+p3(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz56(i5)%e(0)*p3(0)-cz56(i5)%e(1)*p3(1)-cz56(i5)%e(2)
     & *p3(2)-cz56(i5)%e(3)*p3(3)
      cvqd=cz56(i5)%e(0)*p356(0)-cz56(i5)%e(1)*p356(1)-cz56(i5)%
     & e(2)*p356(2)-cz56(i5)%e(3)*p356(3)
      cauxa=-cz56(i5)%ek0*quqd+p3k0*cvqd+p356k0*cvqu
      cauxc=+cz56(i5)%ek0*p3(2)-p3k0*cz56(i5)%e(2)
      l3_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      l3_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      l3_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      l3_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      ccr=fcr(id3)/(f356)
      ccl=fcl(id3)/(f356)
      do i5=1,2
* TL0 -- qu=p3,qd=p356,v=cf56(i5)%e,a=l3_56(i5)%a,c=l3_56(i5)%c,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p3(2)*p356(3)-p356(2)*p3(3))+p3k0*(c
     & f56(i5)%e(2)*p356(3)-p356(2)*cf56(i5)%e(3))-p356k0*(cf56(
     & i5)%e(2)*p3(3)-p3(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p3k0+p3(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf56(i5)%e(0)*p3(0)-cf56(i5)%e(1)*p3(1)-cf56(i5)%e(2)
     & *p3(2)-cf56(i5)%e(3)*p3(3)
      cvqd=cf56(i5)%e(0)*p356(0)-cf56(i5)%e(1)*p356(1)-cf56(i5)%
     & e(2)*p356(2)-cf56(i5)%e(3)*p356(3)
      cauxa=-cf56(i5)%ek0*quqd+p3k0*cvqd+p356k0*cvqu
      cauxc=+cf56(i5)%ek0*p3(2)-p3k0*cf56(i5)%e(2)
      l3_56(i5)%a(1)=l3_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      l3_56(i5)%a(2)=l3_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      l3_56(i5)%c(1)=l3_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      l3_56(i5)%c(2)=l3_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      endif
  
      if (ilept(id4).ne.1.or.ilept(id5).ne.1) then
      quqd=(-s478+s78)/2d0
      do i7=1,2
* TR0 -- qu=p478,qd=p4,v=cz78(i7)%e,a=r4_78(i7)%a,b=r4_78(i7)%b,cr=zcr(i
* d4),cl=zcl(id4),nsum=0
      ceps_0=-cz78(i7)%ek0*(p478(2)*p4(3)-p4(2)*p478(3))+p478k0*
     & (cz78(i7)%e(2)*p4(3)-p4(2)*cz78(i7)%e(3))-p4k0*(cz78(i7)%
     & e(2)*p478(3)-p478(2)*cz78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz78(i7)%e(3)*p4k0+p4(3)*cz78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i7)%e(0)*p478(0)-cz78(i7)%e(1)*p478(1)-cz78(i7)%
     & e(2)*p478(2)-cz78(i7)%e(3)*p478(3)
      cvqd=cz78(i7)%e(0)*p4(0)-cz78(i7)%e(1)*p4(1)-cz78(i7)%e(2)
     & *p4(2)-cz78(i7)%e(3)*p4(3)
      cauxa=-cz78(i7)%ek0*quqd+p478k0*cvqd+p4k0*cvqu
      cauxb=-cz78(i7)%ek0*p4(2)+p4k0*cz78(i7)%e(2)
      r4_78(i7)%a(1)=zcr(id4)*(cauxa+ceps_0)
      r4_78(i7)%a(2)=zcl(id4)*(cauxa-ceps_0)
      r4_78(i7)%b(1)=zcl(id4)*(cauxb-ceps_2)
      r4_78(i7)%b(2)=zcr(id4)*(-cauxb-ceps_2)
      end do
  
      do i7=1,2
* TR0 -- qu=p478,qd=p4,v=cf78(i7)%e,a=r4_78(i7)%a,b=r4_78(i7)%b,cr=fcr(i
* d4),cl=fcl(id4),nsum=1
      ceps_0=-cf78(i7)%ek0*(p478(2)*p4(3)-p4(2)*p478(3))+p478k0*
     & (cf78(i7)%e(2)*p4(3)-p4(2)*cf78(i7)%e(3))-p4k0*(cf78(i7)%
     & e(2)*p478(3)-p478(2)*cf78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf78(i7)%e(3)*p4k0+p4(3)*cf78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i7)%e(0)*p478(0)-cf78(i7)%e(1)*p478(1)-cf78(i7)%
     & e(2)*p478(2)-cf78(i7)%e(3)*p478(3)
      cvqd=cf78(i7)%e(0)*p4(0)-cf78(i7)%e(1)*p4(1)-cf78(i7)%e(2)
     & *p4(2)-cf78(i7)%e(3)*p4(3)
      cauxa=-cf78(i7)%ek0*quqd+p478k0*cvqd+p4k0*cvqu
      cauxb=-cf78(i7)%ek0*p4(2)+p4k0*cf78(i7)%e(2)
      r4_78(i7)%a(1)=r4_78(i7)%a(1)+fcr(id4)*(cauxa+ceps_0)
      r4_78(i7)%a(2)=r4_78(i7)%a(2)+fcl(id4)*(cauxa-ceps_0)
      r4_78(i7)%b(1)=r4_78(i7)%b(1)+fcl(id4)*(cauxb-ceps_2)
      r4_78(i7)%b(2)=r4_78(i7)%b(2)+fcr(id4)*(-cauxb-ceps_2)
      end do
      endif
  
  
       do m=0,3
        p3156(m)=p356(m) +p1(m)
       enddo
* pk0 -- p=p3156
      p3156k0=p3156(0)-p3156(1)
* p.q -- p.q=s3156,p=p3156,q=p3156,bef=,aft=
      s3156=(p3156(0)*p3156(0)-p3156(1)*p3156(1)-p3156(2)*p3156(
     & 2)-p3156(3)*p3156(3))
      f3156=s3156*p3156k0
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p3,q=p478
      quqd=p3(0)*p478(0)-p3(1)*p478(1)-p3(2)*p478(2)-p3(3)*p478(
     & 3)
      ccr=zcr(id3)/(f478)
      ccl=zcl(id3)/(f478)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p3,qd=p478,v=cz5126(i5,i1,i2)%e,a=l3_5126(i5,i1,i2)%a,c=l3_5
* 126(i5,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz5126(i5,i1,i2)%ek0*(p3(2)*p478(3)-p478(2)*p3(3))
     & +p3k0*(cz5126(i5,i1,i2)%e(2)*p478(3)-p478(2)*cz5126(i5,i1
     & ,i2)%e(3))-p478k0*(cz5126(i5,i1,i2)%e(2)*p3(3)-p3(2)*cz51
     & 26(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz5126(i5,i1,i2)%e(3)*p3k0+p3(3)*cz5126(i5,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz5126(i5,i1,i2)%e(0)*p3(0)-cz5126(i5,i1,i2)%e(1)*p3(
     & 1)-cz5126(i5,i1,i2)%e(2)*p3(2)-cz5126(i5,i1,i2)%e(3)*p3(3
     & )
      cvqd=cz5126(i5,i1,i2)%e(0)*p478(0)-cz5126(i5,i1,i2)%e(1)*p
     & 478(1)-cz5126(i5,i1,i2)%e(2)*p478(2)-cz5126(i5,i1,i2)%e(3
     & )*p478(3)
      cauxa=-cz5126(i5,i1,i2)%ek0*quqd+p3k0*cvqd+p478k0*cvqu
      cauxc=+cz5126(i5,i1,i2)%ek0*p3(2)-p3k0*cz5126(i5,i1,i2)%e(
     & 2)
      l3_5126(i5,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l3_5126(i5,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l3_5126(i5,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l3_5126(i5,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id3)/(f478)
      ccl=fcl(id3)/(f478)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p3,qd=p478,v=cf5126(i5,i1,i2)%e,a=l3_5126(i5,i1,i2)%a,c=l3_5
* 126(i5,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf5126(i5,i1,i2)%ek0*(p3(2)*p478(3)-p478(2)*p3(3))
     & +p3k0*(cf5126(i5,i1,i2)%e(2)*p478(3)-p478(2)*cf5126(i5,i1
     & ,i2)%e(3))-p478k0*(cf5126(i5,i1,i2)%e(2)*p3(3)-p3(2)*cf51
     & 26(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf5126(i5,i1,i2)%e(3)*p3k0+p3(3)*cf5126(i5,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf5126(i5,i1,i2)%e(0)*p3(0)-cf5126(i5,i1,i2)%e(1)*p3(
     & 1)-cf5126(i5,i1,i2)%e(2)*p3(2)-cf5126(i5,i1,i2)%e(3)*p3(3
     & )
      cvqd=cf5126(i5,i1,i2)%e(0)*p478(0)-cf5126(i5,i1,i2)%e(1)*p
     & 478(1)-cf5126(i5,i1,i2)%e(2)*p478(2)-cf5126(i5,i1,i2)%e(3
     & )*p478(3)
      cauxa=-cf5126(i5,i1,i2)%ek0*quqd+p3k0*cvqd+p478k0*cvqu
      cauxc=+cf5126(i5,i1,i2)%ek0*p3(2)-p3k0*cf5126(i5,i1,i2)%e(
     & 2)
      l3_5126(i5,i1,i2)%a(1)=l3_5126(i5,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l3_5126(i5,i1,i2)%a(2)=l3_5126(i5,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l3_5126(i5,i1,i2)%c(1)=l3_5126(i5,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l3_5126(i5,i1,i2)%c(2)=l3_5126(i5,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id3).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p3,q=p3156
      quqd=p3(0)*p3156(0)-p3(1)*p3156(1)-p3(2)*p3156(2)-p3(3)*p3
     & 156(3)
      ccr=zcr(id3)/(f3156)
      ccl=zcl(id3)/(f3156)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p3,qd=p3156,v=cz516(i5,i1)%e,a=l3_516(i5,i1)%a,c=l3_516(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz516(i5,i1)%ek0*(p3(2)*p3156(3)-p3156(2)*p3(3))+p
     & 3k0*(cz516(i5,i1)%e(2)*p3156(3)-p3156(2)*cz516(i5,i1)%e(3
     & ))-p3156k0*(cz516(i5,i1)%e(2)*p3(3)-p3(2)*cz516(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz516(i5,i1)%e(3)*p3k0+p3(3)*cz516(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz516(i5,i1)%e(0)*p3(0)-cz516(i5,i1)%e(1)*p3(1)-cz516
     & (i5,i1)%e(2)*p3(2)-cz516(i5,i1)%e(3)*p3(3)
      cvqd=cz516(i5,i1)%e(0)*p3156(0)-cz516(i5,i1)%e(1)*p3156(1)
     & -cz516(i5,i1)%e(2)*p3156(2)-cz516(i5,i1)%e(3)*p3156(3)
      cauxa=-cz516(i5,i1)%ek0*quqd+p3k0*cvqd+p3156k0*cvqu
      cauxc=+cz516(i5,i1)%ek0*p3(2)-p3k0*cz516(i5,i1)%e(2)
      l3_516(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l3_516(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_516(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l3_516(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id3)/(f3156)
      ccl=fcl(id3)/(f3156)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p3,qd=p3156,v=cf516(i5,i1)%e,a=l3_516(i5,i1)%a,c=l3_516(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf516(i5,i1)%ek0*(p3(2)*p3156(3)-p3156(2)*p3(3))+p
     & 3k0*(cf516(i5,i1)%e(2)*p3156(3)-p3156(2)*cf516(i5,i1)%e(3
     & ))-p3156k0*(cf516(i5,i1)%e(2)*p3(3)-p3(2)*cf516(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf516(i5,i1)%e(3)*p3k0+p3(3)*cf516(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf516(i5,i1)%e(0)*p3(0)-cf516(i5,i1)%e(1)*p3(1)-cf516
     & (i5,i1)%e(2)*p3(2)-cf516(i5,i1)%e(3)*p3(3)
      cvqd=cf516(i5,i1)%e(0)*p3156(0)-cf516(i5,i1)%e(1)*p3156(1)
     & -cf516(i5,i1)%e(2)*p3156(2)-cf516(i5,i1)%e(3)*p3156(3)
      cauxa=-cf516(i5,i1)%ek0*quqd+p3k0*cvqd+p3156k0*cvqu
      cauxc=+cf516(i5,i1)%ek0*p3(2)-p3k0*cf516(i5,i1)%e(2)
      l3_516(i5,i1)%a(1)=l3_516(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l3_516(i5,i1)%a(2)=l3_516(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l3_516(i5,i1)%c(1)=l3_516(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l3_516(i5,i1)%c(2)=l3_516(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p356,q=p4
      quqd=p356(0)*p4(0)-p356(1)*p4(1)-p356(2)*p4(2)-p356(3)*p4(
     & 3)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p356,qd=p4,v=cz7128(i7,i1,i2)%e,a=r4_7128(i7,i1,i2)%a,b=r4_7
* 128(i7,i1,i2)%b,cr=zcr(id4),cl=zcl(id4),nsum=0
      ceps_0=-cz7128(i7,i1,i2)%ek0*(p356(2)*p4(3)-p4(2)*p356(3))
     & +p356k0*(cz7128(i7,i1,i2)%e(2)*p4(3)-p4(2)*cz7128(i7,i1,i
     & 2)%e(3))-p4k0*(cz7128(i7,i1,i2)%e(2)*p356(3)-p356(2)*cz71
     & 28(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz7128(i7,i1,i2)%e(3)*p4k0+p4(3)*cz7128(i7,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz7128(i7,i1,i2)%e(0)*p356(0)-cz7128(i7,i1,i2)%e(1)*p
     & 356(1)-cz7128(i7,i1,i2)%e(2)*p356(2)-cz7128(i7,i1,i2)%e(3
     & )*p356(3)
      cvqd=cz7128(i7,i1,i2)%e(0)*p4(0)-cz7128(i7,i1,i2)%e(1)*p4(
     & 1)-cz7128(i7,i1,i2)%e(2)*p4(2)-cz7128(i7,i1,i2)%e(3)*p4(3
     & )
      cauxa=-cz7128(i7,i1,i2)%ek0*quqd+p356k0*cvqd+p4k0*cvqu
      cauxb=-cz7128(i7,i1,i2)%ek0*p4(2)+p4k0*cz7128(i7,i1,i2)%e(
     & 2)
      r4_7128(i7,i1,i2)%a(1)=zcr(id4)*(cauxa+ceps_0)
      r4_7128(i7,i1,i2)%a(2)=zcl(id4)*(cauxa-ceps_0)
      r4_7128(i7,i1,i2)%b(1)=zcl(id4)*(cauxb-ceps_2)
      r4_7128(i7,i1,i2)%b(2)=zcr(id4)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p356,qd=p4,v=cf7128(i7,i1,i2)%e,a=r4_7128(i7,i1,i2)%a,b=r4_7
* 128(i7,i1,i2)%b,cr=fcr(id4),cl=fcl(id4),nsum=1
      ceps_0=-cf7128(i7,i1,i2)%ek0*(p356(2)*p4(3)-p4(2)*p356(3))
     & +p356k0*(cf7128(i7,i1,i2)%e(2)*p4(3)-p4(2)*cf7128(i7,i1,i
     & 2)%e(3))-p4k0*(cf7128(i7,i1,i2)%e(2)*p356(3)-p356(2)*cf71
     & 28(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf7128(i7,i1,i2)%e(3)*p4k0+p4(3)*cf7128(i7,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf7128(i7,i1,i2)%e(0)*p356(0)-cf7128(i7,i1,i2)%e(1)*p
     & 356(1)-cf7128(i7,i1,i2)%e(2)*p356(2)-cf7128(i7,i1,i2)%e(3
     & )*p356(3)
      cvqd=cf7128(i7,i1,i2)%e(0)*p4(0)-cf7128(i7,i1,i2)%e(1)*p4(
     & 1)-cf7128(i7,i1,i2)%e(2)*p4(2)-cf7128(i7,i1,i2)%e(3)*p4(3
     & )
      cauxa=-cf7128(i7,i1,i2)%ek0*quqd+p356k0*cvqd+p4k0*cvqu
      cauxb=-cf7128(i7,i1,i2)%ek0*p4(2)+p4k0*cf7128(i7,i1,i2)%e(
     & 2)
      r4_7128(i7,i1,i2)%a(1)=r4_7128(i7,i1,i2)%a(1)+fcr(id4)*(ca
     & uxa+ceps_0)
      r4_7128(i7,i1,i2)%a(2)=r4_7128(i7,i1,i2)%a(2)+fcl(id4)*(ca
     & uxa-ceps_0)
      r4_7128(i7,i1,i2)%b(1)=r4_7128(i7,i1,i2)%b(1)+fcl(id4)*(ca
     & uxb-ceps_2)
      r4_7128(i7,i1,i2)%b(2)=r4_7128(i7,i1,i2)%b(2)+fcr(id4)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id4).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p3156,q=p4
      quqd=p3156(0)*p4(0)-p3156(1)*p4(1)-p3156(2)*p4(2)-p3156(3)
     & *p4(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p3156,qd=p4,v=cz728(i7,i2)%e,a=r4_728(i7,i2)%a,b=r4_728(i7,i
* 2)%b,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cz728(i7,i2)%ek0*(p3156(2)*p4(3)-p4(2)*p3156(3))+p
     & 3156k0*(cz728(i7,i2)%e(2)*p4(3)-p4(2)*cz728(i7,i2)%e(3))-
     & p4k0*(cz728(i7,i2)%e(2)*p3156(3)-p3156(2)*cz728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz728(i7,i2)%e(3)*p4k0+p4(3)*cz728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz728(i7,i2)%e(0)*p3156(0)-cz728(i7,i2)%e(1)*p3156(1)
     & -cz728(i7,i2)%e(2)*p3156(2)-cz728(i7,i2)%e(3)*p3156(3)
      cvqd=cz728(i7,i2)%e(0)*p4(0)-cz728(i7,i2)%e(1)*p4(1)-cz728
     & (i7,i2)%e(2)*p4(2)-cz728(i7,i2)%e(3)*p4(3)
      cauxa=-cz728(i7,i2)%ek0*quqd+p3156k0*cvqd+p4k0*cvqu
      cauxb=-cz728(i7,i2)%ek0*p4(2)+p4k0*cz728(i7,i2)%e(2)
      r4_728(i7,i2)%a(1)=zcr(id3)*(cauxa+ceps_0)
      r4_728(i7,i2)%a(2)=zcl(id3)*(cauxa-ceps_0)
      r4_728(i7,i2)%b(1)=zcl(id3)*(cauxb-ceps_2)
      r4_728(i7,i2)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p3156,qd=p4,v=cf728(i7,i2)%e,a=r4_728(i7,i2)%a,b=r4_728(i7,i
* 2)%b,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cf728(i7,i2)%ek0*(p3156(2)*p4(3)-p4(2)*p3156(3))+p
     & 3156k0*(cf728(i7,i2)%e(2)*p4(3)-p4(2)*cf728(i7,i2)%e(3))-
     & p4k0*(cf728(i7,i2)%e(2)*p3156(3)-p3156(2)*cf728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf728(i7,i2)%e(3)*p4k0+p4(3)*cf728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf728(i7,i2)%e(0)*p3156(0)-cf728(i7,i2)%e(1)*p3156(1)
     & -cf728(i7,i2)%e(2)*p3156(2)-cf728(i7,i2)%e(3)*p3156(3)
      cvqd=cf728(i7,i2)%e(0)*p4(0)-cf728(i7,i2)%e(1)*p4(1)-cf728
     & (i7,i2)%e(2)*p4(2)-cf728(i7,i2)%e(3)*p4(3)
      cauxa=-cf728(i7,i2)%ek0*quqd+p3156k0*cvqd+p4k0*cvqu
      cauxb=-cf728(i7,i2)%ek0*p4(2)+p4k0*cf728(i7,i2)%e(2)
      r4_728(i7,i2)%a(1)=r4_728(i7,i2)%a(1)+fcr(id3)*(cauxa+ceps
     & _0)
      r4_728(i7,i2)%a(2)=r4_728(i7,i2)%a(2)+fcl(id3)*(cauxa-ceps
     & _0)
      r4_728(i7,i2)%b(1)=r4_728(i7,i2)%b(1)+fcl(id3)*(cauxb-ceps
     & _2)
      r4_728(i7,i2)%b(2)=r4_728(i7,i2)%b(2)+fcr(id3)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p3256(m)=p356(m) +p2(m)
       enddo
* pk0 -- p=p3256
      p3256k0=p3256(0)-p3256(1)
* p.q -- p.q=s3256,p=p3256,q=p3256,bef=,aft=
      s3256=(p3256(0)*p3256(0)-p3256(1)*p3256(1)-p3256(2)*p3256(
     & 2)-p3256(3)*p3256(3))
      f3256=s3256*p3256k0
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p3,q=p478
      quqd=p3(0)*p478(0)-p3(1)*p478(1)-p3(2)*p478(2)-p3(3)*p478(
     & 3)
      ccr=zcr(id3)/(f478)
      ccl=zcl(id3)/(f478)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p3,qd=p478,v=cz5216(i5,i1,i2)%e,a=l3_5216(i5,i1,i2)%a,c=l3_5
* 216(i5,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz5216(i5,i1,i2)%ek0*(p3(2)*p478(3)-p478(2)*p3(3))
     & +p3k0*(cz5216(i5,i1,i2)%e(2)*p478(3)-p478(2)*cz5216(i5,i1
     & ,i2)%e(3))-p478k0*(cz5216(i5,i1,i2)%e(2)*p3(3)-p3(2)*cz52
     & 16(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz5216(i5,i1,i2)%e(3)*p3k0+p3(3)*cz5216(i5,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz5216(i5,i1,i2)%e(0)*p3(0)-cz5216(i5,i1,i2)%e(1)*p3(
     & 1)-cz5216(i5,i1,i2)%e(2)*p3(2)-cz5216(i5,i1,i2)%e(3)*p3(3
     & )
      cvqd=cz5216(i5,i1,i2)%e(0)*p478(0)-cz5216(i5,i1,i2)%e(1)*p
     & 478(1)-cz5216(i5,i1,i2)%e(2)*p478(2)-cz5216(i5,i1,i2)%e(3
     & )*p478(3)
      cauxa=-cz5216(i5,i1,i2)%ek0*quqd+p3k0*cvqd+p478k0*cvqu
      cauxc=+cz5216(i5,i1,i2)%ek0*p3(2)-p3k0*cz5216(i5,i1,i2)%e(
     & 2)
      l3_5216(i5,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l3_5216(i5,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l3_5216(i5,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l3_5216(i5,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id3)/(f478)
      ccl=fcl(id3)/(f478)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p3,qd=p478,v=cf5216(i5,i1,i2)%e,a=l3_5216(i5,i1,i2)%a,c=l3_5
* 216(i5,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf5216(i5,i1,i2)%ek0*(p3(2)*p478(3)-p478(2)*p3(3))
     & +p3k0*(cf5216(i5,i1,i2)%e(2)*p478(3)-p478(2)*cf5216(i5,i1
     & ,i2)%e(3))-p478k0*(cf5216(i5,i1,i2)%e(2)*p3(3)-p3(2)*cf52
     & 16(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf5216(i5,i1,i2)%e(3)*p3k0+p3(3)*cf5216(i5,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf5216(i5,i1,i2)%e(0)*p3(0)-cf5216(i5,i1,i2)%e(1)*p3(
     & 1)-cf5216(i5,i1,i2)%e(2)*p3(2)-cf5216(i5,i1,i2)%e(3)*p3(3
     & )
      cvqd=cf5216(i5,i1,i2)%e(0)*p478(0)-cf5216(i5,i1,i2)%e(1)*p
     & 478(1)-cf5216(i5,i1,i2)%e(2)*p478(2)-cf5216(i5,i1,i2)%e(3
     & )*p478(3)
      cauxa=-cf5216(i5,i1,i2)%ek0*quqd+p3k0*cvqd+p478k0*cvqu
      cauxc=+cf5216(i5,i1,i2)%ek0*p3(2)-p3k0*cf5216(i5,i1,i2)%e(
     & 2)
      l3_5216(i5,i1,i2)%a(1)=l3_5216(i5,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l3_5216(i5,i1,i2)%a(2)=l3_5216(i5,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l3_5216(i5,i1,i2)%c(1)=l3_5216(i5,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l3_5216(i5,i1,i2)%c(2)=l3_5216(i5,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id3).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p3,q=p3256
      quqd=p3(0)*p3256(0)-p3(1)*p3256(1)-p3(2)*p3256(2)-p3(3)*p3
     & 256(3)
      ccr=zcr(id3)/(f3256)
      ccl=zcl(id3)/(f3256)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p3,qd=p3256,v=cz526(i5,i1)%e,a=l3_526(i5,i1)%a,c=l3_526(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz526(i5,i1)%ek0*(p3(2)*p3256(3)-p3256(2)*p3(3))+p
     & 3k0*(cz526(i5,i1)%e(2)*p3256(3)-p3256(2)*cz526(i5,i1)%e(3
     & ))-p3256k0*(cz526(i5,i1)%e(2)*p3(3)-p3(2)*cz526(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz526(i5,i1)%e(3)*p3k0+p3(3)*cz526(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz526(i5,i1)%e(0)*p3(0)-cz526(i5,i1)%e(1)*p3(1)-cz526
     & (i5,i1)%e(2)*p3(2)-cz526(i5,i1)%e(3)*p3(3)
      cvqd=cz526(i5,i1)%e(0)*p3256(0)-cz526(i5,i1)%e(1)*p3256(1)
     & -cz526(i5,i1)%e(2)*p3256(2)-cz526(i5,i1)%e(3)*p3256(3)
      cauxa=-cz526(i5,i1)%ek0*quqd+p3k0*cvqd+p3256k0*cvqu
      cauxc=+cz526(i5,i1)%ek0*p3(2)-p3k0*cz526(i5,i1)%e(2)
      l3_526(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l3_526(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_526(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l3_526(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id3)/(f3256)
      ccl=fcl(id3)/(f3256)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p3,qd=p3256,v=cf526(i5,i1)%e,a=l3_526(i5,i1)%a,c=l3_526(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf526(i5,i1)%ek0*(p3(2)*p3256(3)-p3256(2)*p3(3))+p
     & 3k0*(cf526(i5,i1)%e(2)*p3256(3)-p3256(2)*cf526(i5,i1)%e(3
     & ))-p3256k0*(cf526(i5,i1)%e(2)*p3(3)-p3(2)*cf526(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf526(i5,i1)%e(3)*p3k0+p3(3)*cf526(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf526(i5,i1)%e(0)*p3(0)-cf526(i5,i1)%e(1)*p3(1)-cf526
     & (i5,i1)%e(2)*p3(2)-cf526(i5,i1)%e(3)*p3(3)
      cvqd=cf526(i5,i1)%e(0)*p3256(0)-cf526(i5,i1)%e(1)*p3256(1)
     & -cf526(i5,i1)%e(2)*p3256(2)-cf526(i5,i1)%e(3)*p3256(3)
      cauxa=-cf526(i5,i1)%ek0*quqd+p3k0*cvqd+p3256k0*cvqu
      cauxc=+cf526(i5,i1)%ek0*p3(2)-p3k0*cf526(i5,i1)%e(2)
      l3_526(i5,i1)%a(1)=l3_526(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l3_526(i5,i1)%a(2)=l3_526(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l3_526(i5,i1)%c(1)=l3_526(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l3_526(i5,i1)%c(2)=l3_526(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p356,q=p4
      quqd=p356(0)*p4(0)-p356(1)*p4(1)-p356(2)*p4(2)-p356(3)*p4(
     & 3)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p356,qd=p4,v=cz7218(i7,i1,i2)%e,a=r4_7218(i7,i1,i2)%a,b=r4_7
* 218(i7,i1,i2)%b,cr=zcr(id4),cl=zcl(id4),nsum=0
      ceps_0=-cz7218(i7,i1,i2)%ek0*(p356(2)*p4(3)-p4(2)*p356(3))
     & +p356k0*(cz7218(i7,i1,i2)%e(2)*p4(3)-p4(2)*cz7218(i7,i1,i
     & 2)%e(3))-p4k0*(cz7218(i7,i1,i2)%e(2)*p356(3)-p356(2)*cz72
     & 18(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz7218(i7,i1,i2)%e(3)*p4k0+p4(3)*cz7218(i7,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz7218(i7,i1,i2)%e(0)*p356(0)-cz7218(i7,i1,i2)%e(1)*p
     & 356(1)-cz7218(i7,i1,i2)%e(2)*p356(2)-cz7218(i7,i1,i2)%e(3
     & )*p356(3)
      cvqd=cz7218(i7,i1,i2)%e(0)*p4(0)-cz7218(i7,i1,i2)%e(1)*p4(
     & 1)-cz7218(i7,i1,i2)%e(2)*p4(2)-cz7218(i7,i1,i2)%e(3)*p4(3
     & )
      cauxa=-cz7218(i7,i1,i2)%ek0*quqd+p356k0*cvqd+p4k0*cvqu
      cauxb=-cz7218(i7,i1,i2)%ek0*p4(2)+p4k0*cz7218(i7,i1,i2)%e(
     & 2)
      r4_7218(i7,i1,i2)%a(1)=zcr(id4)*(cauxa+ceps_0)
      r4_7218(i7,i1,i2)%a(2)=zcl(id4)*(cauxa-ceps_0)
      r4_7218(i7,i1,i2)%b(1)=zcl(id4)*(cauxb-ceps_2)
      r4_7218(i7,i1,i2)%b(2)=zcr(id4)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p356,qd=p4,v=cf7218(i7,i1,i2)%e,a=r4_7218(i7,i1,i2)%a,b=r4_7
* 218(i7,i1,i2)%b,cr=fcr(id4),cl=fcl(id4),nsum=1
      ceps_0=-cf7218(i7,i1,i2)%ek0*(p356(2)*p4(3)-p4(2)*p356(3))
     & +p356k0*(cf7218(i7,i1,i2)%e(2)*p4(3)-p4(2)*cf7218(i7,i1,i
     & 2)%e(3))-p4k0*(cf7218(i7,i1,i2)%e(2)*p356(3)-p356(2)*cf72
     & 18(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf7218(i7,i1,i2)%e(3)*p4k0+p4(3)*cf7218(i7,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf7218(i7,i1,i2)%e(0)*p356(0)-cf7218(i7,i1,i2)%e(1)*p
     & 356(1)-cf7218(i7,i1,i2)%e(2)*p356(2)-cf7218(i7,i1,i2)%e(3
     & )*p356(3)
      cvqd=cf7218(i7,i1,i2)%e(0)*p4(0)-cf7218(i7,i1,i2)%e(1)*p4(
     & 1)-cf7218(i7,i1,i2)%e(2)*p4(2)-cf7218(i7,i1,i2)%e(3)*p4(3
     & )
      cauxa=-cf7218(i7,i1,i2)%ek0*quqd+p356k0*cvqd+p4k0*cvqu
      cauxb=-cf7218(i7,i1,i2)%ek0*p4(2)+p4k0*cf7218(i7,i1,i2)%e(
     & 2)
      r4_7218(i7,i1,i2)%a(1)=r4_7218(i7,i1,i2)%a(1)+fcr(id4)*(ca
     & uxa+ceps_0)
      r4_7218(i7,i1,i2)%a(2)=r4_7218(i7,i1,i2)%a(2)+fcl(id4)*(ca
     & uxa-ceps_0)
      r4_7218(i7,i1,i2)%b(1)=r4_7218(i7,i1,i2)%b(1)+fcl(id4)*(ca
     & uxb-ceps_2)
      r4_7218(i7,i1,i2)%b(2)=r4_7218(i7,i1,i2)%b(2)+fcr(id4)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id4).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p3256,q=p4
      quqd=p3256(0)*p4(0)-p3256(1)*p4(1)-p3256(2)*p4(2)-p3256(3)
     & *p4(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p3256,qd=p4,v=cz718(i7,i2)%e,a=r4_718(i7,i2)%a,b=r4_718(i7,i
* 2)%b,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cz718(i7,i2)%ek0*(p3256(2)*p4(3)-p4(2)*p3256(3))+p
     & 3256k0*(cz718(i7,i2)%e(2)*p4(3)-p4(2)*cz718(i7,i2)%e(3))-
     & p4k0*(cz718(i7,i2)%e(2)*p3256(3)-p3256(2)*cz718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz718(i7,i2)%e(3)*p4k0+p4(3)*cz718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz718(i7,i2)%e(0)*p3256(0)-cz718(i7,i2)%e(1)*p3256(1)
     & -cz718(i7,i2)%e(2)*p3256(2)-cz718(i7,i2)%e(3)*p3256(3)
      cvqd=cz718(i7,i2)%e(0)*p4(0)-cz718(i7,i2)%e(1)*p4(1)-cz718
     & (i7,i2)%e(2)*p4(2)-cz718(i7,i2)%e(3)*p4(3)
      cauxa=-cz718(i7,i2)%ek0*quqd+p3256k0*cvqd+p4k0*cvqu
      cauxb=-cz718(i7,i2)%ek0*p4(2)+p4k0*cz718(i7,i2)%e(2)
      r4_718(i7,i2)%a(1)=zcr(id3)*(cauxa+ceps_0)
      r4_718(i7,i2)%a(2)=zcl(id3)*(cauxa-ceps_0)
      r4_718(i7,i2)%b(1)=zcl(id3)*(cauxb-ceps_2)
      r4_718(i7,i2)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p3256,qd=p4,v=cf718(i7,i2)%e,a=r4_718(i7,i2)%a,b=r4_718(i7,i
* 2)%b,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cf718(i7,i2)%ek0*(p3256(2)*p4(3)-p4(2)*p3256(3))+p
     & 3256k0*(cf718(i7,i2)%e(2)*p4(3)-p4(2)*cf718(i7,i2)%e(3))-
     & p4k0*(cf718(i7,i2)%e(2)*p3256(3)-p3256(2)*cf718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf718(i7,i2)%e(3)*p4k0+p4(3)*cf718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf718(i7,i2)%e(0)*p3256(0)-cf718(i7,i2)%e(1)*p3256(1)
     & -cf718(i7,i2)%e(2)*p3256(2)-cf718(i7,i2)%e(3)*p3256(3)
      cvqd=cf718(i7,i2)%e(0)*p4(0)-cf718(i7,i2)%e(1)*p4(1)-cf718
     & (i7,i2)%e(2)*p4(2)-cf718(i7,i2)%e(3)*p4(3)
      cauxa=-cf718(i7,i2)%ek0*quqd+p3256k0*cvqd+p4k0*cvqu
      cauxb=-cf718(i7,i2)%ek0*p4(2)+p4k0*cf718(i7,i2)%e(2)
      r4_718(i7,i2)%a(1)=r4_718(i7,i2)%a(1)+fcr(id3)*(cauxa+ceps
     & _0)
      r4_718(i7,i2)%a(2)=r4_718(i7,i2)%a(2)+fcl(id3)*(cauxa-ceps
     & _0)
      r4_718(i7,i2)%b(1)=r4_718(i7,i2)%b(1)+fcl(id3)*(cauxb-ceps
     & _2)
      r4_718(i7,i2)%b(2)=r4_718(i7,i2)%b(2)+fcr(id3)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p378(m)=p3(m)+p7(m)+p8(m)
       enddo
* pk0 -- p=p378
      p378k0=p378(0)-p378(1)
* p.q -- p.q=s378,p=p378,q=p378,bef=,aft=
      s378=(p378(0)*p378(0)-p378(1)*p378(1)-p378(2)*p378(2)-p378
     & (3)*p378(3))
      f378=s378*p378k0
  
       do m=0,3
        p456(m)= -p4(m)-p5(m)-p6(m)
       enddo
* pk0 -- p=p456
      p456k0=p456(0)-p456(1)
* p.q -- p.q=s456,p=p456,q=p456,bef=,aft=
      s456=(p456(0)*p456(0)-p456(1)*p456(1)-p456(2)*p456(2)-p456
     & (3)*p456(3))
      f456=s456*p456k0
  
      if (ilept(id3).ne.1.or.ilept(id5).ne.1) then
      quqd=(s378-s78)/2d0
      ccr=zcr(id3)/(f378)
      ccl=zcl(id3)/(f378)
      do i5=1,2
* TL0 -- qu=p3,qd=p378,v=cz78(i5)%e,a=l3_78(i5)%a,c=l3_78(i5)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p3(2)*p378(3)-p378(2)*p3(3))+p3k0*(c
     & z78(i5)%e(2)*p378(3)-p378(2)*cz78(i5)%e(3))-p378k0*(cz78(
     & i5)%e(2)*p3(3)-p3(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p3k0+p3(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz78(i5)%e(0)*p3(0)-cz78(i5)%e(1)*p3(1)-cz78(i5)%e(2)
     & *p3(2)-cz78(i5)%e(3)*p3(3)
      cvqd=cz78(i5)%e(0)*p378(0)-cz78(i5)%e(1)*p378(1)-cz78(i5)%
     & e(2)*p378(2)-cz78(i5)%e(3)*p378(3)
      cauxa=-cz78(i5)%ek0*quqd+p3k0*cvqd+p378k0*cvqu
      cauxc=+cz78(i5)%ek0*p3(2)-p3k0*cz78(i5)%e(2)
      l3_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      l3_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      l3_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      l3_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      ccr=fcr(id3)/(f378)
      ccl=fcl(id3)/(f378)
      do i5=1,2
* TL0 -- qu=p3,qd=p378,v=cf78(i5)%e,a=l3_78(i5)%a,c=l3_78(i5)%c,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p3(2)*p378(3)-p378(2)*p3(3))+p3k0*(c
     & f78(i5)%e(2)*p378(3)-p378(2)*cf78(i5)%e(3))-p378k0*(cf78(
     & i5)%e(2)*p3(3)-p3(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p3k0+p3(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf78(i5)%e(0)*p3(0)-cf78(i5)%e(1)*p3(1)-cf78(i5)%e(2)
     & *p3(2)-cf78(i5)%e(3)*p3(3)
      cvqd=cf78(i5)%e(0)*p378(0)-cf78(i5)%e(1)*p378(1)-cf78(i5)%
     & e(2)*p378(2)-cf78(i5)%e(3)*p378(3)
      cauxa=-cf78(i5)%ek0*quqd+p3k0*cvqd+p378k0*cvqu
      cauxc=+cf78(i5)%ek0*p3(2)-p3k0*cf78(i5)%e(2)
      l3_78(i5)%a(1)=l3_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      l3_78(i5)%a(2)=l3_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      l3_78(i5)%c(1)=l3_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      l3_78(i5)%c(2)=l3_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      endif
  
      if (ilept(id4).ne.1.or.ilept(id7).ne.1) then
      quqd=(-s456+s56)/2d0
      do i7=1,2
* TR0 -- qu=p456,qd=p4,v=cz56(i7)%e,a=r4_56(i7)%a,b=r4_56(i7)%b,cr=zcr(i
* d4),cl=zcl(id4),nsum=0
      ceps_0=-cz56(i7)%ek0*(p456(2)*p4(3)-p4(2)*p456(3))+p456k0*
     & (cz56(i7)%e(2)*p4(3)-p4(2)*cz56(i7)%e(3))-p4k0*(cz56(i7)%
     & e(2)*p456(3)-p456(2)*cz56(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz56(i7)%e(3)*p4k0+p4(3)*cz56(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i7)%e(0)*p456(0)-cz56(i7)%e(1)*p456(1)-cz56(i7)%
     & e(2)*p456(2)-cz56(i7)%e(3)*p456(3)
      cvqd=cz56(i7)%e(0)*p4(0)-cz56(i7)%e(1)*p4(1)-cz56(i7)%e(2)
     & *p4(2)-cz56(i7)%e(3)*p4(3)
      cauxa=-cz56(i7)%ek0*quqd+p456k0*cvqd+p4k0*cvqu
      cauxb=-cz56(i7)%ek0*p4(2)+p4k0*cz56(i7)%e(2)
      r4_56(i7)%a(1)=zcr(id4)*(cauxa+ceps_0)
      r4_56(i7)%a(2)=zcl(id4)*(cauxa-ceps_0)
      r4_56(i7)%b(1)=zcl(id4)*(cauxb-ceps_2)
      r4_56(i7)%b(2)=zcr(id4)*(-cauxb-ceps_2)
      end do
  
      do i7=1,2
* TR0 -- qu=p456,qd=p4,v=cf56(i7)%e,a=r4_56(i7)%a,b=r4_56(i7)%b,cr=fcr(i
* d4),cl=fcl(id4),nsum=1
      ceps_0=-cf56(i7)%ek0*(p456(2)*p4(3)-p4(2)*p456(3))+p456k0*
     & (cf56(i7)%e(2)*p4(3)-p4(2)*cf56(i7)%e(3))-p4k0*(cf56(i7)%
     & e(2)*p456(3)-p456(2)*cf56(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf56(i7)%e(3)*p4k0+p4(3)*cf56(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i7)%e(0)*p456(0)-cf56(i7)%e(1)*p456(1)-cf56(i7)%
     & e(2)*p456(2)-cf56(i7)%e(3)*p456(3)
      cvqd=cf56(i7)%e(0)*p4(0)-cf56(i7)%e(1)*p4(1)-cf56(i7)%e(2)
     & *p4(2)-cf56(i7)%e(3)*p4(3)
      cauxa=-cf56(i7)%ek0*quqd+p456k0*cvqd+p4k0*cvqu
      cauxb=-cf56(i7)%ek0*p4(2)+p4k0*cf56(i7)%e(2)
      r4_56(i7)%a(1)=r4_56(i7)%a(1)+fcr(id4)*(cauxa+ceps_0)
      r4_56(i7)%a(2)=r4_56(i7)%a(2)+fcl(id4)*(cauxa-ceps_0)
      r4_56(i7)%b(1)=r4_56(i7)%b(1)+fcl(id4)*(cauxb-ceps_2)
      r4_56(i7)%b(2)=r4_56(i7)%b(2)+fcr(id4)*(-cauxb-ceps_2)
      end do
      endif
  
  
       do m=0,3
        p3178(m)=p378(m) +p1(m)
       enddo
* pk0 -- p=p3178
      p3178k0=p3178(0)-p3178(1)
* p.q -- p.q=s3178,p=p3178,q=p3178,bef=,aft=
      s3178=(p3178(0)*p3178(0)-p3178(1)*p3178(1)-p3178(2)*p3178(
     & 2)-p3178(3)*p3178(3))
      f3178=s3178*p3178k0
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p3,q=p456
      quqd=p3(0)*p456(0)-p3(1)*p456(1)-p3(2)*p456(2)-p3(3)*p456(
     & 3)
      ccr=zcr(id3)/(f456)
      ccl=zcl(id3)/(f456)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p3,qd=p456,v=cz7128(i7,i1,i2)%e,a=l3_7128(i7,i1,i2)%a,c=l3_7
* 128(i7,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz7128(i7,i1,i2)%ek0*(p3(2)*p456(3)-p456(2)*p3(3))
     & +p3k0*(cz7128(i7,i1,i2)%e(2)*p456(3)-p456(2)*cz7128(i7,i1
     & ,i2)%e(3))-p456k0*(cz7128(i7,i1,i2)%e(2)*p3(3)-p3(2)*cz71
     & 28(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz7128(i7,i1,i2)%e(3)*p3k0+p3(3)*cz7128(i7,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz7128(i7,i1,i2)%e(0)*p3(0)-cz7128(i7,i1,i2)%e(1)*p3(
     & 1)-cz7128(i7,i1,i2)%e(2)*p3(2)-cz7128(i7,i1,i2)%e(3)*p3(3
     & )
      cvqd=cz7128(i7,i1,i2)%e(0)*p456(0)-cz7128(i7,i1,i2)%e(1)*p
     & 456(1)-cz7128(i7,i1,i2)%e(2)*p456(2)-cz7128(i7,i1,i2)%e(3
     & )*p456(3)
      cauxa=-cz7128(i7,i1,i2)%ek0*quqd+p3k0*cvqd+p456k0*cvqu
      cauxc=+cz7128(i7,i1,i2)%ek0*p3(2)-p3k0*cz7128(i7,i1,i2)%e(
     & 2)
      l3_7128(i7,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l3_7128(i7,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l3_7128(i7,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l3_7128(i7,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id3)/(f456)
      ccl=fcl(id3)/(f456)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p3,qd=p456,v=cf7128(i7,i1,i2)%e,a=l3_7128(i7,i1,i2)%a,c=l3_7
* 128(i7,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf7128(i7,i1,i2)%ek0*(p3(2)*p456(3)-p456(2)*p3(3))
     & +p3k0*(cf7128(i7,i1,i2)%e(2)*p456(3)-p456(2)*cf7128(i7,i1
     & ,i2)%e(3))-p456k0*(cf7128(i7,i1,i2)%e(2)*p3(3)-p3(2)*cf71
     & 28(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf7128(i7,i1,i2)%e(3)*p3k0+p3(3)*cf7128(i7,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf7128(i7,i1,i2)%e(0)*p3(0)-cf7128(i7,i1,i2)%e(1)*p3(
     & 1)-cf7128(i7,i1,i2)%e(2)*p3(2)-cf7128(i7,i1,i2)%e(3)*p3(3
     & )
      cvqd=cf7128(i7,i1,i2)%e(0)*p456(0)-cf7128(i7,i1,i2)%e(1)*p
     & 456(1)-cf7128(i7,i1,i2)%e(2)*p456(2)-cf7128(i7,i1,i2)%e(3
     & )*p456(3)
      cauxa=-cf7128(i7,i1,i2)%ek0*quqd+p3k0*cvqd+p456k0*cvqu
      cauxc=+cf7128(i7,i1,i2)%ek0*p3(2)-p3k0*cf7128(i7,i1,i2)%e(
     & 2)
      l3_7128(i7,i1,i2)%a(1)=l3_7128(i7,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l3_7128(i7,i1,i2)%a(2)=l3_7128(i7,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l3_7128(i7,i1,i2)%c(1)=l3_7128(i7,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l3_7128(i7,i1,i2)%c(2)=l3_7128(i7,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id3).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p3,q=p3178
      quqd=p3(0)*p3178(0)-p3(1)*p3178(1)-p3(2)*p3178(2)-p3(3)*p3
     & 178(3)
      ccr=zcr(id3)/(f3178)
      ccl=zcl(id3)/(f3178)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p3,qd=p3178,v=cz718(i5,i1)%e,a=l3_718(i5,i1)%a,c=l3_718(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz718(i5,i1)%ek0*(p3(2)*p3178(3)-p3178(2)*p3(3))+p
     & 3k0*(cz718(i5,i1)%e(2)*p3178(3)-p3178(2)*cz718(i5,i1)%e(3
     & ))-p3178k0*(cz718(i5,i1)%e(2)*p3(3)-p3(2)*cz718(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz718(i5,i1)%e(3)*p3k0+p3(3)*cz718(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz718(i5,i1)%e(0)*p3(0)-cz718(i5,i1)%e(1)*p3(1)-cz718
     & (i5,i1)%e(2)*p3(2)-cz718(i5,i1)%e(3)*p3(3)
      cvqd=cz718(i5,i1)%e(0)*p3178(0)-cz718(i5,i1)%e(1)*p3178(1)
     & -cz718(i5,i1)%e(2)*p3178(2)-cz718(i5,i1)%e(3)*p3178(3)
      cauxa=-cz718(i5,i1)%ek0*quqd+p3k0*cvqd+p3178k0*cvqu
      cauxc=+cz718(i5,i1)%ek0*p3(2)-p3k0*cz718(i5,i1)%e(2)
      l3_718(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l3_718(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_718(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l3_718(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id3)/(f3178)
      ccl=fcl(id3)/(f3178)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p3,qd=p3178,v=cf718(i5,i1)%e,a=l3_718(i5,i1)%a,c=l3_718(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf718(i5,i1)%ek0*(p3(2)*p3178(3)-p3178(2)*p3(3))+p
     & 3k0*(cf718(i5,i1)%e(2)*p3178(3)-p3178(2)*cf718(i5,i1)%e(3
     & ))-p3178k0*(cf718(i5,i1)%e(2)*p3(3)-p3(2)*cf718(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf718(i5,i1)%e(3)*p3k0+p3(3)*cf718(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf718(i5,i1)%e(0)*p3(0)-cf718(i5,i1)%e(1)*p3(1)-cf718
     & (i5,i1)%e(2)*p3(2)-cf718(i5,i1)%e(3)*p3(3)
      cvqd=cf718(i5,i1)%e(0)*p3178(0)-cf718(i5,i1)%e(1)*p3178(1)
     & -cf718(i5,i1)%e(2)*p3178(2)-cf718(i5,i1)%e(3)*p3178(3)
      cauxa=-cf718(i5,i1)%ek0*quqd+p3k0*cvqd+p3178k0*cvqu
      cauxc=+cf718(i5,i1)%ek0*p3(2)-p3k0*cf718(i5,i1)%e(2)
      l3_718(i5,i1)%a(1)=l3_718(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l3_718(i5,i1)%a(2)=l3_718(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l3_718(i5,i1)%c(1)=l3_718(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l3_718(i5,i1)%c(2)=l3_718(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p378,q=p4
      quqd=p378(0)*p4(0)-p378(1)*p4(1)-p378(2)*p4(2)-p378(3)*p4(
     & 3)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p378,qd=p4,v=cz5126(i5,i1,i2)%e,a=r4_5126(i5,i1,i2)%a,b=r4_5
* 126(i5,i1,i2)%b,cr=zcr(id4),cl=zcl(id4),nsum=0
      ceps_0=-cz5126(i5,i1,i2)%ek0*(p378(2)*p4(3)-p4(2)*p378(3))
     & +p378k0*(cz5126(i5,i1,i2)%e(2)*p4(3)-p4(2)*cz5126(i5,i1,i
     & 2)%e(3))-p4k0*(cz5126(i5,i1,i2)%e(2)*p378(3)-p378(2)*cz51
     & 26(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz5126(i5,i1,i2)%e(3)*p4k0+p4(3)*cz5126(i5,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz5126(i5,i1,i2)%e(0)*p378(0)-cz5126(i5,i1,i2)%e(1)*p
     & 378(1)-cz5126(i5,i1,i2)%e(2)*p378(2)-cz5126(i5,i1,i2)%e(3
     & )*p378(3)
      cvqd=cz5126(i5,i1,i2)%e(0)*p4(0)-cz5126(i5,i1,i2)%e(1)*p4(
     & 1)-cz5126(i5,i1,i2)%e(2)*p4(2)-cz5126(i5,i1,i2)%e(3)*p4(3
     & )
      cauxa=-cz5126(i5,i1,i2)%ek0*quqd+p378k0*cvqd+p4k0*cvqu
      cauxb=-cz5126(i5,i1,i2)%ek0*p4(2)+p4k0*cz5126(i5,i1,i2)%e(
     & 2)
      r4_5126(i5,i1,i2)%a(1)=zcr(id4)*(cauxa+ceps_0)
      r4_5126(i5,i1,i2)%a(2)=zcl(id4)*(cauxa-ceps_0)
      r4_5126(i5,i1,i2)%b(1)=zcl(id4)*(cauxb-ceps_2)
      r4_5126(i5,i1,i2)%b(2)=zcr(id4)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p378,qd=p4,v=cf5126(i5,i1,i2)%e,a=r4_5126(i5,i1,i2)%a,b=r4_5
* 126(i5,i1,i2)%b,cr=fcr(id4),cl=fcl(id4),nsum=1
      ceps_0=-cf5126(i5,i1,i2)%ek0*(p378(2)*p4(3)-p4(2)*p378(3))
     & +p378k0*(cf5126(i5,i1,i2)%e(2)*p4(3)-p4(2)*cf5126(i5,i1,i
     & 2)%e(3))-p4k0*(cf5126(i5,i1,i2)%e(2)*p378(3)-p378(2)*cf51
     & 26(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf5126(i5,i1,i2)%e(3)*p4k0+p4(3)*cf5126(i5,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf5126(i5,i1,i2)%e(0)*p378(0)-cf5126(i5,i1,i2)%e(1)*p
     & 378(1)-cf5126(i5,i1,i2)%e(2)*p378(2)-cf5126(i5,i1,i2)%e(3
     & )*p378(3)
      cvqd=cf5126(i5,i1,i2)%e(0)*p4(0)-cf5126(i5,i1,i2)%e(1)*p4(
     & 1)-cf5126(i5,i1,i2)%e(2)*p4(2)-cf5126(i5,i1,i2)%e(3)*p4(3
     & )
      cauxa=-cf5126(i5,i1,i2)%ek0*quqd+p378k0*cvqd+p4k0*cvqu
      cauxb=-cf5126(i5,i1,i2)%ek0*p4(2)+p4k0*cf5126(i5,i1,i2)%e(
     & 2)
      r4_5126(i5,i1,i2)%a(1)=r4_5126(i5,i1,i2)%a(1)+fcr(id4)*(ca
     & uxa+ceps_0)
      r4_5126(i5,i1,i2)%a(2)=r4_5126(i5,i1,i2)%a(2)+fcl(id4)*(ca
     & uxa-ceps_0)
      r4_5126(i5,i1,i2)%b(1)=r4_5126(i5,i1,i2)%b(1)+fcl(id4)*(ca
     & uxb-ceps_2)
      r4_5126(i5,i1,i2)%b(2)=r4_5126(i5,i1,i2)%b(2)+fcr(id4)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id4).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p3178,q=p4
      quqd=p3178(0)*p4(0)-p3178(1)*p4(1)-p3178(2)*p4(2)-p3178(3)
     & *p4(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p3178,qd=p4,v=cz526(i7,i2)%e,a=r4_526(i7,i2)%a,b=r4_526(i7,i
* 2)%b,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cz526(i7,i2)%ek0*(p3178(2)*p4(3)-p4(2)*p3178(3))+p
     & 3178k0*(cz526(i7,i2)%e(2)*p4(3)-p4(2)*cz526(i7,i2)%e(3))-
     & p4k0*(cz526(i7,i2)%e(2)*p3178(3)-p3178(2)*cz526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz526(i7,i2)%e(3)*p4k0+p4(3)*cz526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz526(i7,i2)%e(0)*p3178(0)-cz526(i7,i2)%e(1)*p3178(1)
     & -cz526(i7,i2)%e(2)*p3178(2)-cz526(i7,i2)%e(3)*p3178(3)
      cvqd=cz526(i7,i2)%e(0)*p4(0)-cz526(i7,i2)%e(1)*p4(1)-cz526
     & (i7,i2)%e(2)*p4(2)-cz526(i7,i2)%e(3)*p4(3)
      cauxa=-cz526(i7,i2)%ek0*quqd+p3178k0*cvqd+p4k0*cvqu
      cauxb=-cz526(i7,i2)%ek0*p4(2)+p4k0*cz526(i7,i2)%e(2)
      r4_526(i7,i2)%a(1)=zcr(id3)*(cauxa+ceps_0)
      r4_526(i7,i2)%a(2)=zcl(id3)*(cauxa-ceps_0)
      r4_526(i7,i2)%b(1)=zcl(id3)*(cauxb-ceps_2)
      r4_526(i7,i2)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p3178,qd=p4,v=cf526(i7,i2)%e,a=r4_526(i7,i2)%a,b=r4_526(i7,i
* 2)%b,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cf526(i7,i2)%ek0*(p3178(2)*p4(3)-p4(2)*p3178(3))+p
     & 3178k0*(cf526(i7,i2)%e(2)*p4(3)-p4(2)*cf526(i7,i2)%e(3))-
     & p4k0*(cf526(i7,i2)%e(2)*p3178(3)-p3178(2)*cf526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf526(i7,i2)%e(3)*p4k0+p4(3)*cf526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf526(i7,i2)%e(0)*p3178(0)-cf526(i7,i2)%e(1)*p3178(1)
     & -cf526(i7,i2)%e(2)*p3178(2)-cf526(i7,i2)%e(3)*p3178(3)
      cvqd=cf526(i7,i2)%e(0)*p4(0)-cf526(i7,i2)%e(1)*p4(1)-cf526
     & (i7,i2)%e(2)*p4(2)-cf526(i7,i2)%e(3)*p4(3)
      cauxa=-cf526(i7,i2)%ek0*quqd+p3178k0*cvqd+p4k0*cvqu
      cauxb=-cf526(i7,i2)%ek0*p4(2)+p4k0*cf526(i7,i2)%e(2)
      r4_526(i7,i2)%a(1)=r4_526(i7,i2)%a(1)+fcr(id3)*(cauxa+ceps
     & _0)
      r4_526(i7,i2)%a(2)=r4_526(i7,i2)%a(2)+fcl(id3)*(cauxa-ceps
     & _0)
      r4_526(i7,i2)%b(1)=r4_526(i7,i2)%b(1)+fcl(id3)*(cauxb-ceps
     & _2)
      r4_526(i7,i2)%b(2)=r4_526(i7,i2)%b(2)+fcr(id3)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p3278(m)=p378(m) +p2(m)
       enddo
* pk0 -- p=p3278
      p3278k0=p3278(0)-p3278(1)
* p.q -- p.q=s3278,p=p3278,q=p3278,bef=,aft=
      s3278=(p3278(0)*p3278(0)-p3278(1)*p3278(1)-p3278(2)*p3278(
     & 2)-p3278(3)*p3278(3))
      f3278=s3278*p3278k0
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p3,q=p456
      quqd=p3(0)*p456(0)-p3(1)*p456(1)-p3(2)*p456(2)-p3(3)*p456(
     & 3)
      ccr=zcr(id3)/(f456)
      ccl=zcl(id3)/(f456)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p3,qd=p456,v=cz7218(i7,i1,i2)%e,a=l3_7218(i7,i1,i2)%a,c=l3_7
* 218(i7,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz7218(i7,i1,i2)%ek0*(p3(2)*p456(3)-p456(2)*p3(3))
     & +p3k0*(cz7218(i7,i1,i2)%e(2)*p456(3)-p456(2)*cz7218(i7,i1
     & ,i2)%e(3))-p456k0*(cz7218(i7,i1,i2)%e(2)*p3(3)-p3(2)*cz72
     & 18(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz7218(i7,i1,i2)%e(3)*p3k0+p3(3)*cz7218(i7,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz7218(i7,i1,i2)%e(0)*p3(0)-cz7218(i7,i1,i2)%e(1)*p3(
     & 1)-cz7218(i7,i1,i2)%e(2)*p3(2)-cz7218(i7,i1,i2)%e(3)*p3(3
     & )
      cvqd=cz7218(i7,i1,i2)%e(0)*p456(0)-cz7218(i7,i1,i2)%e(1)*p
     & 456(1)-cz7218(i7,i1,i2)%e(2)*p456(2)-cz7218(i7,i1,i2)%e(3
     & )*p456(3)
      cauxa=-cz7218(i7,i1,i2)%ek0*quqd+p3k0*cvqd+p456k0*cvqu
      cauxc=+cz7218(i7,i1,i2)%ek0*p3(2)-p3k0*cz7218(i7,i1,i2)%e(
     & 2)
      l3_7218(i7,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l3_7218(i7,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l3_7218(i7,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l3_7218(i7,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id3)/(f456)
      ccl=fcl(id3)/(f456)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p3,qd=p456,v=cf7218(i7,i1,i2)%e,a=l3_7218(i7,i1,i2)%a,c=l3_7
* 218(i7,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf7218(i7,i1,i2)%ek0*(p3(2)*p456(3)-p456(2)*p3(3))
     & +p3k0*(cf7218(i7,i1,i2)%e(2)*p456(3)-p456(2)*cf7218(i7,i1
     & ,i2)%e(3))-p456k0*(cf7218(i7,i1,i2)%e(2)*p3(3)-p3(2)*cf72
     & 18(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf7218(i7,i1,i2)%e(3)*p3k0+p3(3)*cf7218(i7,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf7218(i7,i1,i2)%e(0)*p3(0)-cf7218(i7,i1,i2)%e(1)*p3(
     & 1)-cf7218(i7,i1,i2)%e(2)*p3(2)-cf7218(i7,i1,i2)%e(3)*p3(3
     & )
      cvqd=cf7218(i7,i1,i2)%e(0)*p456(0)-cf7218(i7,i1,i2)%e(1)*p
     & 456(1)-cf7218(i7,i1,i2)%e(2)*p456(2)-cf7218(i7,i1,i2)%e(3
     & )*p456(3)
      cauxa=-cf7218(i7,i1,i2)%ek0*quqd+p3k0*cvqd+p456k0*cvqu
      cauxc=+cf7218(i7,i1,i2)%ek0*p3(2)-p3k0*cf7218(i7,i1,i2)%e(
     & 2)
      l3_7218(i7,i1,i2)%a(1)=l3_7218(i7,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l3_7218(i7,i1,i2)%a(2)=l3_7218(i7,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l3_7218(i7,i1,i2)%c(1)=l3_7218(i7,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l3_7218(i7,i1,i2)%c(2)=l3_7218(i7,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id3).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p3,q=p3278
      quqd=p3(0)*p3278(0)-p3(1)*p3278(1)-p3(2)*p3278(2)-p3(3)*p3
     & 278(3)
      ccr=zcr(id3)/(f3278)
      ccl=zcl(id3)/(f3278)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p3,qd=p3278,v=cz728(i5,i1)%e,a=l3_728(i5,i1)%a,c=l3_728(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz728(i5,i1)%ek0*(p3(2)*p3278(3)-p3278(2)*p3(3))+p
     & 3k0*(cz728(i5,i1)%e(2)*p3278(3)-p3278(2)*cz728(i5,i1)%e(3
     & ))-p3278k0*(cz728(i5,i1)%e(2)*p3(3)-p3(2)*cz728(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz728(i5,i1)%e(3)*p3k0+p3(3)*cz728(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz728(i5,i1)%e(0)*p3(0)-cz728(i5,i1)%e(1)*p3(1)-cz728
     & (i5,i1)%e(2)*p3(2)-cz728(i5,i1)%e(3)*p3(3)
      cvqd=cz728(i5,i1)%e(0)*p3278(0)-cz728(i5,i1)%e(1)*p3278(1)
     & -cz728(i5,i1)%e(2)*p3278(2)-cz728(i5,i1)%e(3)*p3278(3)
      cauxa=-cz728(i5,i1)%ek0*quqd+p3k0*cvqd+p3278k0*cvqu
      cauxc=+cz728(i5,i1)%ek0*p3(2)-p3k0*cz728(i5,i1)%e(2)
      l3_728(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l3_728(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_728(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l3_728(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id3)/(f3278)
      ccl=fcl(id3)/(f3278)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p3,qd=p3278,v=cf728(i5,i1)%e,a=l3_728(i5,i1)%a,c=l3_728(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf728(i5,i1)%ek0*(p3(2)*p3278(3)-p3278(2)*p3(3))+p
     & 3k0*(cf728(i5,i1)%e(2)*p3278(3)-p3278(2)*cf728(i5,i1)%e(3
     & ))-p3278k0*(cf728(i5,i1)%e(2)*p3(3)-p3(2)*cf728(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf728(i5,i1)%e(3)*p3k0+p3(3)*cf728(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf728(i5,i1)%e(0)*p3(0)-cf728(i5,i1)%e(1)*p3(1)-cf728
     & (i5,i1)%e(2)*p3(2)-cf728(i5,i1)%e(3)*p3(3)
      cvqd=cf728(i5,i1)%e(0)*p3278(0)-cf728(i5,i1)%e(1)*p3278(1)
     & -cf728(i5,i1)%e(2)*p3278(2)-cf728(i5,i1)%e(3)*p3278(3)
      cauxa=-cf728(i5,i1)%ek0*quqd+p3k0*cvqd+p3278k0*cvqu
      cauxc=+cf728(i5,i1)%ek0*p3(2)-p3k0*cf728(i5,i1)%e(2)
      l3_728(i5,i1)%a(1)=l3_728(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l3_728(i5,i1)%a(2)=l3_728(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l3_728(i5,i1)%c(1)=l3_728(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l3_728(i5,i1)%c(2)=l3_728(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p378,q=p4
      quqd=p378(0)*p4(0)-p378(1)*p4(1)-p378(2)*p4(2)-p378(3)*p4(
     & 3)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p378,qd=p4,v=cz5216(i5,i1,i2)%e,a=r4_5216(i5,i1,i2)%a,b=r4_5
* 216(i5,i1,i2)%b,cr=zcr(id4),cl=zcl(id4),nsum=0
      ceps_0=-cz5216(i5,i1,i2)%ek0*(p378(2)*p4(3)-p4(2)*p378(3))
     & +p378k0*(cz5216(i5,i1,i2)%e(2)*p4(3)-p4(2)*cz5216(i5,i1,i
     & 2)%e(3))-p4k0*(cz5216(i5,i1,i2)%e(2)*p378(3)-p378(2)*cz52
     & 16(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz5216(i5,i1,i2)%e(3)*p4k0+p4(3)*cz5216(i5,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz5216(i5,i1,i2)%e(0)*p378(0)-cz5216(i5,i1,i2)%e(1)*p
     & 378(1)-cz5216(i5,i1,i2)%e(2)*p378(2)-cz5216(i5,i1,i2)%e(3
     & )*p378(3)
      cvqd=cz5216(i5,i1,i2)%e(0)*p4(0)-cz5216(i5,i1,i2)%e(1)*p4(
     & 1)-cz5216(i5,i1,i2)%e(2)*p4(2)-cz5216(i5,i1,i2)%e(3)*p4(3
     & )
      cauxa=-cz5216(i5,i1,i2)%ek0*quqd+p378k0*cvqd+p4k0*cvqu
      cauxb=-cz5216(i5,i1,i2)%ek0*p4(2)+p4k0*cz5216(i5,i1,i2)%e(
     & 2)
      r4_5216(i5,i1,i2)%a(1)=zcr(id4)*(cauxa+ceps_0)
      r4_5216(i5,i1,i2)%a(2)=zcl(id4)*(cauxa-ceps_0)
      r4_5216(i5,i1,i2)%b(1)=zcl(id4)*(cauxb-ceps_2)
      r4_5216(i5,i1,i2)%b(2)=zcr(id4)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p378,qd=p4,v=cf5216(i5,i1,i2)%e,a=r4_5216(i5,i1,i2)%a,b=r4_5
* 216(i5,i1,i2)%b,cr=fcr(id4),cl=fcl(id4),nsum=1
      ceps_0=-cf5216(i5,i1,i2)%ek0*(p378(2)*p4(3)-p4(2)*p378(3))
     & +p378k0*(cf5216(i5,i1,i2)%e(2)*p4(3)-p4(2)*cf5216(i5,i1,i
     & 2)%e(3))-p4k0*(cf5216(i5,i1,i2)%e(2)*p378(3)-p378(2)*cf52
     & 16(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf5216(i5,i1,i2)%e(3)*p4k0+p4(3)*cf5216(i5,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf5216(i5,i1,i2)%e(0)*p378(0)-cf5216(i5,i1,i2)%e(1)*p
     & 378(1)-cf5216(i5,i1,i2)%e(2)*p378(2)-cf5216(i5,i1,i2)%e(3
     & )*p378(3)
      cvqd=cf5216(i5,i1,i2)%e(0)*p4(0)-cf5216(i5,i1,i2)%e(1)*p4(
     & 1)-cf5216(i5,i1,i2)%e(2)*p4(2)-cf5216(i5,i1,i2)%e(3)*p4(3
     & )
      cauxa=-cf5216(i5,i1,i2)%ek0*quqd+p378k0*cvqd+p4k0*cvqu
      cauxb=-cf5216(i5,i1,i2)%ek0*p4(2)+p4k0*cf5216(i5,i1,i2)%e(
     & 2)
      r4_5216(i5,i1,i2)%a(1)=r4_5216(i5,i1,i2)%a(1)+fcr(id4)*(ca
     & uxa+ceps_0)
      r4_5216(i5,i1,i2)%a(2)=r4_5216(i5,i1,i2)%a(2)+fcl(id4)*(ca
     & uxa-ceps_0)
      r4_5216(i5,i1,i2)%b(1)=r4_5216(i5,i1,i2)%b(1)+fcl(id4)*(ca
     & uxb-ceps_2)
      r4_5216(i5,i1,i2)%b(2)=r4_5216(i5,i1,i2)%b(2)+fcr(id4)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id4).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p3278,q=p4
      quqd=p3278(0)*p4(0)-p3278(1)*p4(1)-p3278(2)*p4(2)-p3278(3)
     & *p4(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p3278,qd=p4,v=cz516(i7,i2)%e,a=r4_516(i7,i2)%a,b=r4_516(i7,i
* 2)%b,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cz516(i7,i2)%ek0*(p3278(2)*p4(3)-p4(2)*p3278(3))+p
     & 3278k0*(cz516(i7,i2)%e(2)*p4(3)-p4(2)*cz516(i7,i2)%e(3))-
     & p4k0*(cz516(i7,i2)%e(2)*p3278(3)-p3278(2)*cz516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz516(i7,i2)%e(3)*p4k0+p4(3)*cz516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz516(i7,i2)%e(0)*p3278(0)-cz516(i7,i2)%e(1)*p3278(1)
     & -cz516(i7,i2)%e(2)*p3278(2)-cz516(i7,i2)%e(3)*p3278(3)
      cvqd=cz516(i7,i2)%e(0)*p4(0)-cz516(i7,i2)%e(1)*p4(1)-cz516
     & (i7,i2)%e(2)*p4(2)-cz516(i7,i2)%e(3)*p4(3)
      cauxa=-cz516(i7,i2)%ek0*quqd+p3278k0*cvqd+p4k0*cvqu
      cauxb=-cz516(i7,i2)%ek0*p4(2)+p4k0*cz516(i7,i2)%e(2)
      r4_516(i7,i2)%a(1)=zcr(id3)*(cauxa+ceps_0)
      r4_516(i7,i2)%a(2)=zcl(id3)*(cauxa-ceps_0)
      r4_516(i7,i2)%b(1)=zcl(id3)*(cauxb-ceps_2)
      r4_516(i7,i2)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p3278,qd=p4,v=cf516(i7,i2)%e,a=r4_516(i7,i2)%a,b=r4_516(i7,i
* 2)%b,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cf516(i7,i2)%ek0*(p3278(2)*p4(3)-p4(2)*p3278(3))+p
     & 3278k0*(cf516(i7,i2)%e(2)*p4(3)-p4(2)*cf516(i7,i2)%e(3))-
     & p4k0*(cf516(i7,i2)%e(2)*p3278(3)-p3278(2)*cf516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf516(i7,i2)%e(3)*p4k0+p4(3)*cf516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf516(i7,i2)%e(0)*p3278(0)-cf516(i7,i2)%e(1)*p3278(1)
     & -cf516(i7,i2)%e(2)*p3278(2)-cf516(i7,i2)%e(3)*p3278(3)
      cvqd=cf516(i7,i2)%e(0)*p4(0)-cf516(i7,i2)%e(1)*p4(1)-cf516
     & (i7,i2)%e(2)*p4(2)-cf516(i7,i2)%e(3)*p4(3)
      cauxa=-cf516(i7,i2)%ek0*quqd+p3278k0*cvqd+p4k0*cvqu
      cauxb=-cf516(i7,i2)%ek0*p4(2)+p4k0*cf516(i7,i2)%e(2)
      r4_516(i7,i2)%a(1)=r4_516(i7,i2)%a(1)+fcr(id3)*(cauxa+ceps
     & _0)
      r4_516(i7,i2)%a(2)=r4_516(i7,i2)%a(2)+fcl(id3)*(cauxa-ceps
     & _0)
      r4_516(i7,i2)%b(1)=r4_516(i7,i2)%b(1)+fcl(id3)*(cauxb-ceps
     & _2)
      r4_516(i7,i2)%b(2)=r4_516(i7,i2)%b(2)+fcr(id3)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p534(m)=p5(m)+p3(m)+p4(m)
       enddo
* pk0 -- p=p534
      p534k0=p534(0)-p534(1)
* p.q -- p.q=s534,p=p534,q=p534,bef=,aft=
      s534=(p534(0)*p534(0)-p534(1)*p534(1)-p534(2)*p534(2)-p534
     & (3)*p534(3))
      f534=s534*p534k0
  
       do m=0,3
        p678(m)= -p6(m)-p7(m)-p8(m)
       enddo
* pk0 -- p=p678
      p678k0=p678(0)-p678(1)
* p.q -- p.q=s678,p=p678,q=p678,bef=,aft=
      s678=(p678(0)*p678(0)-p678(1)*p678(1)-p678(2)*p678(2)-p678
     & (3)*p678(3))
      f678=s678*p678k0
  
      if (ilept(id5).ne.1.or.ilept(id7).ne.1) then
      quqd=(s534-s34)/2d0
      ccr=zcr(id5)/(f534)
      ccl=zcl(id5)/(f534)
      do i5=1,2
* TL0 -- qu=p5,qd=p534,v=cz34(i5)%e,a=l5_34(i5)%a,c=l5_34(i5)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*(c
     & z34(i5)%e(2)*p534(3)-p534(2)*cz34(i5)%e(3))-p534k0*(cz34(
     & i5)%e(2)*p5(3)-p5(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p5k0+p5(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i5)%e(0)*p5(0)-cz34(i5)%e(1)*p5(1)-cz34(i5)%e(2)
     & *p5(2)-cz34(i5)%e(3)*p5(3)
      cvqd=cz34(i5)%e(0)*p534(0)-cz34(i5)%e(1)*p534(1)-cz34(i5)%
     & e(2)*p534(2)-cz34(i5)%e(3)*p534(3)
      cauxa=-cz34(i5)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cz34(i5)%ek0*p5(2)-p5k0*cz34(i5)%e(2)
      l5_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      l5_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      l5_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      l5_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      ccr=fcr(id5)/(f534)
      ccl=fcl(id5)/(f534)
      do i5=1,2
* TL0 -- qu=p5,qd=p534,v=cf34(i5)%e,a=l5_34(i5)%a,c=l5_34(i5)%c,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*(c
     & f34(i5)%e(2)*p534(3)-p534(2)*cf34(i5)%e(3))-p534k0*(cf34(
     & i5)%e(2)*p5(3)-p5(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p5k0+p5(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf34(i5)%e(0)*p5(0)-cf34(i5)%e(1)*p5(1)-cf34(i5)%e(2)
     & *p5(2)-cf34(i5)%e(3)*p5(3)
      cvqd=cf34(i5)%e(0)*p534(0)-cf34(i5)%e(1)*p534(1)-cf34(i5)%
     & e(2)*p534(2)-cf34(i5)%e(3)*p534(3)
      cauxa=-cf34(i5)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cf34(i5)%ek0*p5(2)-p5k0*cf34(i5)%e(2)
      l5_34(i5)%a(1)=l5_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      l5_34(i5)%a(2)=l5_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      l5_34(i5)%c(1)=l5_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      l5_34(i5)%c(2)=l5_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      endif
  
      if (ilept(id6).ne.1.or.ilept(id3).ne.1) then
      quqd=(-s678+s78)/2d0
      do i7=1,2
* TR0 -- qu=p678,qd=p6,v=cz78(i7)%e,a=r6_78(i7)%a,b=r6_78(i7)%b,cr=zcr(i
* d6),cl=zcl(id6),nsum=0
      ceps_0=-cz78(i7)%ek0*(p678(2)*p6(3)-p6(2)*p678(3))+p678k0*
     & (cz78(i7)%e(2)*p6(3)-p6(2)*cz78(i7)%e(3))-p6k0*(cz78(i7)%
     & e(2)*p678(3)-p678(2)*cz78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz78(i7)%e(3)*p6k0+p6(3)*cz78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i7)%e(0)*p678(0)-cz78(i7)%e(1)*p678(1)-cz78(i7)%
     & e(2)*p678(2)-cz78(i7)%e(3)*p678(3)
      cvqd=cz78(i7)%e(0)*p6(0)-cz78(i7)%e(1)*p6(1)-cz78(i7)%e(2)
     & *p6(2)-cz78(i7)%e(3)*p6(3)
      cauxa=-cz78(i7)%ek0*quqd+p678k0*cvqd+p6k0*cvqu
      cauxb=-cz78(i7)%ek0*p6(2)+p6k0*cz78(i7)%e(2)
      r6_78(i7)%a(1)=zcr(id6)*(cauxa+ceps_0)
      r6_78(i7)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_78(i7)%b(1)=zcl(id6)*(cauxb-ceps_2)
      r6_78(i7)%b(2)=zcr(id6)*(-cauxb-ceps_2)
      end do
  
      do i7=1,2
* TR0 -- qu=p678,qd=p6,v=cf78(i7)%e,a=r6_78(i7)%a,b=r6_78(i7)%b,cr=fcr(i
* d6),cl=fcl(id6),nsum=1
      ceps_0=-cf78(i7)%ek0*(p678(2)*p6(3)-p6(2)*p678(3))+p678k0*
     & (cf78(i7)%e(2)*p6(3)-p6(2)*cf78(i7)%e(3))-p6k0*(cf78(i7)%
     & e(2)*p678(3)-p678(2)*cf78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf78(i7)%e(3)*p6k0+p6(3)*cf78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i7)%e(0)*p678(0)-cf78(i7)%e(1)*p678(1)-cf78(i7)%
     & e(2)*p678(2)-cf78(i7)%e(3)*p678(3)
      cvqd=cf78(i7)%e(0)*p6(0)-cf78(i7)%e(1)*p6(1)-cf78(i7)%e(2)
     & *p6(2)-cf78(i7)%e(3)*p6(3)
      cauxa=-cf78(i7)%ek0*quqd+p678k0*cvqd+p6k0*cvqu
      cauxb=-cf78(i7)%ek0*p6(2)+p6k0*cf78(i7)%e(2)
      r6_78(i7)%a(1)=r6_78(i7)%a(1)+fcr(id6)*(cauxa+ceps_0)
      r6_78(i7)%a(2)=r6_78(i7)%a(2)+fcl(id6)*(cauxa-ceps_0)
      r6_78(i7)%b(1)=r6_78(i7)%b(1)+fcl(id6)*(cauxb-ceps_2)
      r6_78(i7)%b(2)=r6_78(i7)%b(2)+fcr(id6)*(-cauxb-ceps_2)
      end do
      endif
  
  
       do m=0,3
        p5134(m)=p534(m) +p1(m)
       enddo
* pk0 -- p=p5134
      p5134k0=p5134(0)-p5134(1)
* p.q -- p.q=s5134,p=p5134,q=p5134,bef=,aft=
      s5134=(p5134(0)*p5134(0)-p5134(1)*p5134(1)-p5134(2)*p5134(
     & 2)-p5134(3)*p5134(3))
      f5134=s5134*p5134k0
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p5,q=p678
      quqd=p5(0)*p678(0)-p5(1)*p678(1)-p5(2)*p678(2)-p5(3)*p678(
     & 3)
      ccr=zcr(id5)/(f678)
      ccl=zcl(id5)/(f678)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p5,qd=p678,v=cz3124(i3,i1,i2)%e,a=l5_3124(i3,i1,i2)%a,c=l5_3
* 124(i3,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz3124(i3,i1,i2)%ek0*(p5(2)*p678(3)-p678(2)*p5(3))
     & +p5k0*(cz3124(i3,i1,i2)%e(2)*p678(3)-p678(2)*cz3124(i3,i1
     & ,i2)%e(3))-p678k0*(cz3124(i3,i1,i2)%e(2)*p5(3)-p5(2)*cz31
     & 24(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz3124(i3,i1,i2)%e(3)*p5k0+p5(3)*cz3124(i3,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz3124(i3,i1,i2)%e(0)*p5(0)-cz3124(i3,i1,i2)%e(1)*p5(
     & 1)-cz3124(i3,i1,i2)%e(2)*p5(2)-cz3124(i3,i1,i2)%e(3)*p5(3
     & )
      cvqd=cz3124(i3,i1,i2)%e(0)*p678(0)-cz3124(i3,i1,i2)%e(1)*p
     & 678(1)-cz3124(i3,i1,i2)%e(2)*p678(2)-cz3124(i3,i1,i2)%e(3
     & )*p678(3)
      cauxa=-cz3124(i3,i1,i2)%ek0*quqd+p5k0*cvqd+p678k0*cvqu
      cauxc=+cz3124(i3,i1,i2)%ek0*p5(2)-p5k0*cz3124(i3,i1,i2)%e(
     & 2)
      l5_3124(i3,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l5_3124(i3,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_3124(i3,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l5_3124(i3,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id5)/(f678)
      ccl=fcl(id5)/(f678)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p5,qd=p678,v=cf3124(i3,i1,i2)%e,a=l5_3124(i3,i1,i2)%a,c=l5_3
* 124(i3,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf3124(i3,i1,i2)%ek0*(p5(2)*p678(3)-p678(2)*p5(3))
     & +p5k0*(cf3124(i3,i1,i2)%e(2)*p678(3)-p678(2)*cf3124(i3,i1
     & ,i2)%e(3))-p678k0*(cf3124(i3,i1,i2)%e(2)*p5(3)-p5(2)*cf31
     & 24(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf3124(i3,i1,i2)%e(3)*p5k0+p5(3)*cf3124(i3,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf3124(i3,i1,i2)%e(0)*p5(0)-cf3124(i3,i1,i2)%e(1)*p5(
     & 1)-cf3124(i3,i1,i2)%e(2)*p5(2)-cf3124(i3,i1,i2)%e(3)*p5(3
     & )
      cvqd=cf3124(i3,i1,i2)%e(0)*p678(0)-cf3124(i3,i1,i2)%e(1)*p
     & 678(1)-cf3124(i3,i1,i2)%e(2)*p678(2)-cf3124(i3,i1,i2)%e(3
     & )*p678(3)
      cauxa=-cf3124(i3,i1,i2)%ek0*quqd+p5k0*cvqd+p678k0*cvqu
      cauxc=+cf3124(i3,i1,i2)%ek0*p5(2)-p5k0*cf3124(i3,i1,i2)%e(
     & 2)
      l5_3124(i3,i1,i2)%a(1)=l5_3124(i3,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l5_3124(i3,i1,i2)%a(2)=l5_3124(i3,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l5_3124(i3,i1,i2)%c(1)=l5_3124(i3,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l5_3124(i3,i1,i2)%c(2)=l5_3124(i3,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id5).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p5,q=p5134
      quqd=p5(0)*p5134(0)-p5(1)*p5134(1)-p5(2)*p5134(2)-p5(3)*p5
     & 134(3)
      ccr=zcr(id5)/(f5134)
      ccl=zcl(id5)/(f5134)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p5,qd=p5134,v=cz314(i5,i1)%e,a=l5_314(i5,i1)%a,c=l5_314(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz314(i5,i1)%ek0*(p5(2)*p5134(3)-p5134(2)*p5(3))+p
     & 5k0*(cz314(i5,i1)%e(2)*p5134(3)-p5134(2)*cz314(i5,i1)%e(3
     & ))-p5134k0*(cz314(i5,i1)%e(2)*p5(3)-p5(2)*cz314(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i5,i1)%e(3)*p5k0+p5(3)*cz314(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz314(i5,i1)%e(0)*p5(0)-cz314(i5,i1)%e(1)*p5(1)-cz314
     & (i5,i1)%e(2)*p5(2)-cz314(i5,i1)%e(3)*p5(3)
      cvqd=cz314(i5,i1)%e(0)*p5134(0)-cz314(i5,i1)%e(1)*p5134(1)
     & -cz314(i5,i1)%e(2)*p5134(2)-cz314(i5,i1)%e(3)*p5134(3)
      cauxa=-cz314(i5,i1)%ek0*quqd+p5k0*cvqd+p5134k0*cvqu
      cauxc=+cz314(i5,i1)%ek0*p5(2)-p5k0*cz314(i5,i1)%e(2)
      l5_314(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l5_314(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_314(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l5_314(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id5)/(f5134)
      ccl=fcl(id5)/(f5134)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p5,qd=p5134,v=cf314(i5,i1)%e,a=l5_314(i5,i1)%a,c=l5_314(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf314(i5,i1)%ek0*(p5(2)*p5134(3)-p5134(2)*p5(3))+p
     & 5k0*(cf314(i5,i1)%e(2)*p5134(3)-p5134(2)*cf314(i5,i1)%e(3
     & ))-p5134k0*(cf314(i5,i1)%e(2)*p5(3)-p5(2)*cf314(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i5,i1)%e(3)*p5k0+p5(3)*cf314(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf314(i5,i1)%e(0)*p5(0)-cf314(i5,i1)%e(1)*p5(1)-cf314
     & (i5,i1)%e(2)*p5(2)-cf314(i5,i1)%e(3)*p5(3)
      cvqd=cf314(i5,i1)%e(0)*p5134(0)-cf314(i5,i1)%e(1)*p5134(1)
     & -cf314(i5,i1)%e(2)*p5134(2)-cf314(i5,i1)%e(3)*p5134(3)
      cauxa=-cf314(i5,i1)%ek0*quqd+p5k0*cvqd+p5134k0*cvqu
      cauxc=+cf314(i5,i1)%ek0*p5(2)-p5k0*cf314(i5,i1)%e(2)
      l5_314(i5,i1)%a(1)=l5_314(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l5_314(i5,i1)%a(2)=l5_314(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l5_314(i5,i1)%c(1)=l5_314(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l5_314(i5,i1)%c(2)=l5_314(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p534,q=p6
      quqd=p534(0)*p6(0)-p534(1)*p6(1)-p534(2)*p6(2)-p534(3)*p6(
     & 3)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p534,qd=p6,v=cz7128(i7,i1,i2)%e,a=r6_7128(i7,i1,i2)%a,b=r6_7
* 128(i7,i1,i2)%b,cr=zcr(id6),cl=zcl(id6),nsum=0
      ceps_0=-cz7128(i7,i1,i2)%ek0*(p534(2)*p6(3)-p6(2)*p534(3))
     & +p534k0*(cz7128(i7,i1,i2)%e(2)*p6(3)-p6(2)*cz7128(i7,i1,i
     & 2)%e(3))-p6k0*(cz7128(i7,i1,i2)%e(2)*p534(3)-p534(2)*cz71
     & 28(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz7128(i7,i1,i2)%e(3)*p6k0+p6(3)*cz7128(i7,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz7128(i7,i1,i2)%e(0)*p534(0)-cz7128(i7,i1,i2)%e(1)*p
     & 534(1)-cz7128(i7,i1,i2)%e(2)*p534(2)-cz7128(i7,i1,i2)%e(3
     & )*p534(3)
      cvqd=cz7128(i7,i1,i2)%e(0)*p6(0)-cz7128(i7,i1,i2)%e(1)*p6(
     & 1)-cz7128(i7,i1,i2)%e(2)*p6(2)-cz7128(i7,i1,i2)%e(3)*p6(3
     & )
      cauxa=-cz7128(i7,i1,i2)%ek0*quqd+p534k0*cvqd+p6k0*cvqu
      cauxb=-cz7128(i7,i1,i2)%ek0*p6(2)+p6k0*cz7128(i7,i1,i2)%e(
     & 2)
      r6_7128(i7,i1,i2)%a(1)=zcr(id6)*(cauxa+ceps_0)
      r6_7128(i7,i1,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_7128(i7,i1,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      r6_7128(i7,i1,i2)%b(2)=zcr(id6)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p534,qd=p6,v=cf7128(i7,i1,i2)%e,a=r6_7128(i7,i1,i2)%a,b=r6_7
* 128(i7,i1,i2)%b,cr=fcr(id6),cl=fcl(id6),nsum=1
      ceps_0=-cf7128(i7,i1,i2)%ek0*(p534(2)*p6(3)-p6(2)*p534(3))
     & +p534k0*(cf7128(i7,i1,i2)%e(2)*p6(3)-p6(2)*cf7128(i7,i1,i
     & 2)%e(3))-p6k0*(cf7128(i7,i1,i2)%e(2)*p534(3)-p534(2)*cf71
     & 28(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf7128(i7,i1,i2)%e(3)*p6k0+p6(3)*cf7128(i7,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf7128(i7,i1,i2)%e(0)*p534(0)-cf7128(i7,i1,i2)%e(1)*p
     & 534(1)-cf7128(i7,i1,i2)%e(2)*p534(2)-cf7128(i7,i1,i2)%e(3
     & )*p534(3)
      cvqd=cf7128(i7,i1,i2)%e(0)*p6(0)-cf7128(i7,i1,i2)%e(1)*p6(
     & 1)-cf7128(i7,i1,i2)%e(2)*p6(2)-cf7128(i7,i1,i2)%e(3)*p6(3
     & )
      cauxa=-cf7128(i7,i1,i2)%ek0*quqd+p534k0*cvqd+p6k0*cvqu
      cauxb=-cf7128(i7,i1,i2)%ek0*p6(2)+p6k0*cf7128(i7,i1,i2)%e(
     & 2)
      r6_7128(i7,i1,i2)%a(1)=r6_7128(i7,i1,i2)%a(1)+fcr(id6)*(ca
     & uxa+ceps_0)
      r6_7128(i7,i1,i2)%a(2)=r6_7128(i7,i1,i2)%a(2)+fcl(id6)*(ca
     & uxa-ceps_0)
      r6_7128(i7,i1,i2)%b(1)=r6_7128(i7,i1,i2)%b(1)+fcl(id6)*(ca
     & uxb-ceps_2)
      r6_7128(i7,i1,i2)%b(2)=r6_7128(i7,i1,i2)%b(2)+fcr(id6)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id6).ne.1.or.ilept(id3).ne.1) then
  
* quqd -- p=p5134,q=p6
      quqd=p5134(0)*p6(0)-p5134(1)*p6(1)-p5134(2)*p6(2)-p5134(3)
     & *p6(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p5134,qd=p6,v=cz728(i7,i2)%e,a=r6_728(i7,i2)%a,b=r6_728(i7,i
* 2)%b,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cz728(i7,i2)%ek0*(p5134(2)*p6(3)-p6(2)*p5134(3))+p
     & 5134k0*(cz728(i7,i2)%e(2)*p6(3)-p6(2)*cz728(i7,i2)%e(3))-
     & p6k0*(cz728(i7,i2)%e(2)*p5134(3)-p5134(2)*cz728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz728(i7,i2)%e(3)*p6k0+p6(3)*cz728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz728(i7,i2)%e(0)*p5134(0)-cz728(i7,i2)%e(1)*p5134(1)
     & -cz728(i7,i2)%e(2)*p5134(2)-cz728(i7,i2)%e(3)*p5134(3)
      cvqd=cz728(i7,i2)%e(0)*p6(0)-cz728(i7,i2)%e(1)*p6(1)-cz728
     & (i7,i2)%e(2)*p6(2)-cz728(i7,i2)%e(3)*p6(3)
      cauxa=-cz728(i7,i2)%ek0*quqd+p5134k0*cvqd+p6k0*cvqu
      cauxb=-cz728(i7,i2)%ek0*p6(2)+p6k0*cz728(i7,i2)%e(2)
      r6_728(i7,i2)%a(1)=zcr(id5)*(cauxa+ceps_0)
      r6_728(i7,i2)%a(2)=zcl(id5)*(cauxa-ceps_0)
      r6_728(i7,i2)%b(1)=zcl(id5)*(cauxb-ceps_2)
      r6_728(i7,i2)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p5134,qd=p6,v=cf728(i7,i2)%e,a=r6_728(i7,i2)%a,b=r6_728(i7,i
* 2)%b,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cf728(i7,i2)%ek0*(p5134(2)*p6(3)-p6(2)*p5134(3))+p
     & 5134k0*(cf728(i7,i2)%e(2)*p6(3)-p6(2)*cf728(i7,i2)%e(3))-
     & p6k0*(cf728(i7,i2)%e(2)*p5134(3)-p5134(2)*cf728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf728(i7,i2)%e(3)*p6k0+p6(3)*cf728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf728(i7,i2)%e(0)*p5134(0)-cf728(i7,i2)%e(1)*p5134(1)
     & -cf728(i7,i2)%e(2)*p5134(2)-cf728(i7,i2)%e(3)*p5134(3)
      cvqd=cf728(i7,i2)%e(0)*p6(0)-cf728(i7,i2)%e(1)*p6(1)-cf728
     & (i7,i2)%e(2)*p6(2)-cf728(i7,i2)%e(3)*p6(3)
      cauxa=-cf728(i7,i2)%ek0*quqd+p5134k0*cvqd+p6k0*cvqu
      cauxb=-cf728(i7,i2)%ek0*p6(2)+p6k0*cf728(i7,i2)%e(2)
      r6_728(i7,i2)%a(1)=r6_728(i7,i2)%a(1)+fcr(id5)*(cauxa+ceps
     & _0)
      r6_728(i7,i2)%a(2)=r6_728(i7,i2)%a(2)+fcl(id5)*(cauxa-ceps
     & _0)
      r6_728(i7,i2)%b(1)=r6_728(i7,i2)%b(1)+fcl(id5)*(cauxb-ceps
     & _2)
      r6_728(i7,i2)%b(2)=r6_728(i7,i2)%b(2)+fcr(id5)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p5234(m)=p534(m) +p2(m)
       enddo
* pk0 -- p=p5234
      p5234k0=p5234(0)-p5234(1)
* p.q -- p.q=s5234,p=p5234,q=p5234,bef=,aft=
      s5234=(p5234(0)*p5234(0)-p5234(1)*p5234(1)-p5234(2)*p5234(
     & 2)-p5234(3)*p5234(3))
      f5234=s5234*p5234k0
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p5,q=p678
      quqd=p5(0)*p678(0)-p5(1)*p678(1)-p5(2)*p678(2)-p5(3)*p678(
     & 3)
      ccr=zcr(id5)/(f678)
      ccl=zcl(id5)/(f678)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p5,qd=p678,v=cz3214(i3,i1,i2)%e,a=l5_3214(i3,i1,i2)%a,c=l5_3
* 214(i3,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz3214(i3,i1,i2)%ek0*(p5(2)*p678(3)-p678(2)*p5(3))
     & +p5k0*(cz3214(i3,i1,i2)%e(2)*p678(3)-p678(2)*cz3214(i3,i1
     & ,i2)%e(3))-p678k0*(cz3214(i3,i1,i2)%e(2)*p5(3)-p5(2)*cz32
     & 14(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz3214(i3,i1,i2)%e(3)*p5k0+p5(3)*cz3214(i3,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz3214(i3,i1,i2)%e(0)*p5(0)-cz3214(i3,i1,i2)%e(1)*p5(
     & 1)-cz3214(i3,i1,i2)%e(2)*p5(2)-cz3214(i3,i1,i2)%e(3)*p5(3
     & )
      cvqd=cz3214(i3,i1,i2)%e(0)*p678(0)-cz3214(i3,i1,i2)%e(1)*p
     & 678(1)-cz3214(i3,i1,i2)%e(2)*p678(2)-cz3214(i3,i1,i2)%e(3
     & )*p678(3)
      cauxa=-cz3214(i3,i1,i2)%ek0*quqd+p5k0*cvqd+p678k0*cvqu
      cauxc=+cz3214(i3,i1,i2)%ek0*p5(2)-p5k0*cz3214(i3,i1,i2)%e(
     & 2)
      l5_3214(i3,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l5_3214(i3,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_3214(i3,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l5_3214(i3,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id5)/(f678)
      ccl=fcl(id5)/(f678)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p5,qd=p678,v=cf3214(i3,i1,i2)%e,a=l5_3214(i3,i1,i2)%a,c=l5_3
* 214(i3,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf3214(i3,i1,i2)%ek0*(p5(2)*p678(3)-p678(2)*p5(3))
     & +p5k0*(cf3214(i3,i1,i2)%e(2)*p678(3)-p678(2)*cf3214(i3,i1
     & ,i2)%e(3))-p678k0*(cf3214(i3,i1,i2)%e(2)*p5(3)-p5(2)*cf32
     & 14(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf3214(i3,i1,i2)%e(3)*p5k0+p5(3)*cf3214(i3,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf3214(i3,i1,i2)%e(0)*p5(0)-cf3214(i3,i1,i2)%e(1)*p5(
     & 1)-cf3214(i3,i1,i2)%e(2)*p5(2)-cf3214(i3,i1,i2)%e(3)*p5(3
     & )
      cvqd=cf3214(i3,i1,i2)%e(0)*p678(0)-cf3214(i3,i1,i2)%e(1)*p
     & 678(1)-cf3214(i3,i1,i2)%e(2)*p678(2)-cf3214(i3,i1,i2)%e(3
     & )*p678(3)
      cauxa=-cf3214(i3,i1,i2)%ek0*quqd+p5k0*cvqd+p678k0*cvqu
      cauxc=+cf3214(i3,i1,i2)%ek0*p5(2)-p5k0*cf3214(i3,i1,i2)%e(
     & 2)
      l5_3214(i3,i1,i2)%a(1)=l5_3214(i3,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l5_3214(i3,i1,i2)%a(2)=l5_3214(i3,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l5_3214(i3,i1,i2)%c(1)=l5_3214(i3,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l5_3214(i3,i1,i2)%c(2)=l5_3214(i3,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id5).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p5,q=p5234
      quqd=p5(0)*p5234(0)-p5(1)*p5234(1)-p5(2)*p5234(2)-p5(3)*p5
     & 234(3)
      ccr=zcr(id5)/(f5234)
      ccl=zcl(id5)/(f5234)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p5,qd=p5234,v=cz324(i5,i1)%e,a=l5_324(i5,i1)%a,c=l5_324(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz324(i5,i1)%ek0*(p5(2)*p5234(3)-p5234(2)*p5(3))+p
     & 5k0*(cz324(i5,i1)%e(2)*p5234(3)-p5234(2)*cz324(i5,i1)%e(3
     & ))-p5234k0*(cz324(i5,i1)%e(2)*p5(3)-p5(2)*cz324(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i5,i1)%e(3)*p5k0+p5(3)*cz324(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz324(i5,i1)%e(0)*p5(0)-cz324(i5,i1)%e(1)*p5(1)-cz324
     & (i5,i1)%e(2)*p5(2)-cz324(i5,i1)%e(3)*p5(3)
      cvqd=cz324(i5,i1)%e(0)*p5234(0)-cz324(i5,i1)%e(1)*p5234(1)
     & -cz324(i5,i1)%e(2)*p5234(2)-cz324(i5,i1)%e(3)*p5234(3)
      cauxa=-cz324(i5,i1)%ek0*quqd+p5k0*cvqd+p5234k0*cvqu
      cauxc=+cz324(i5,i1)%ek0*p5(2)-p5k0*cz324(i5,i1)%e(2)
      l5_324(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l5_324(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_324(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l5_324(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id5)/(f5234)
      ccl=fcl(id5)/(f5234)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p5,qd=p5234,v=cf324(i5,i1)%e,a=l5_324(i5,i1)%a,c=l5_324(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf324(i5,i1)%ek0*(p5(2)*p5234(3)-p5234(2)*p5(3))+p
     & 5k0*(cf324(i5,i1)%e(2)*p5234(3)-p5234(2)*cf324(i5,i1)%e(3
     & ))-p5234k0*(cf324(i5,i1)%e(2)*p5(3)-p5(2)*cf324(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i5,i1)%e(3)*p5k0+p5(3)*cf324(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf324(i5,i1)%e(0)*p5(0)-cf324(i5,i1)%e(1)*p5(1)-cf324
     & (i5,i1)%e(2)*p5(2)-cf324(i5,i1)%e(3)*p5(3)
      cvqd=cf324(i5,i1)%e(0)*p5234(0)-cf324(i5,i1)%e(1)*p5234(1)
     & -cf324(i5,i1)%e(2)*p5234(2)-cf324(i5,i1)%e(3)*p5234(3)
      cauxa=-cf324(i5,i1)%ek0*quqd+p5k0*cvqd+p5234k0*cvqu
      cauxc=+cf324(i5,i1)%ek0*p5(2)-p5k0*cf324(i5,i1)%e(2)
      l5_324(i5,i1)%a(1)=l5_324(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l5_324(i5,i1)%a(2)=l5_324(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l5_324(i5,i1)%c(1)=l5_324(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l5_324(i5,i1)%c(2)=l5_324(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p534,q=p6
      quqd=p534(0)*p6(0)-p534(1)*p6(1)-p534(2)*p6(2)-p534(3)*p6(
     & 3)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p534,qd=p6,v=cz7218(i7,i1,i2)%e,a=r6_7218(i7,i1,i2)%a,b=r6_7
* 218(i7,i1,i2)%b,cr=zcr(id6),cl=zcl(id6),nsum=0
      ceps_0=-cz7218(i7,i1,i2)%ek0*(p534(2)*p6(3)-p6(2)*p534(3))
     & +p534k0*(cz7218(i7,i1,i2)%e(2)*p6(3)-p6(2)*cz7218(i7,i1,i
     & 2)%e(3))-p6k0*(cz7218(i7,i1,i2)%e(2)*p534(3)-p534(2)*cz72
     & 18(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz7218(i7,i1,i2)%e(3)*p6k0+p6(3)*cz7218(i7,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz7218(i7,i1,i2)%e(0)*p534(0)-cz7218(i7,i1,i2)%e(1)*p
     & 534(1)-cz7218(i7,i1,i2)%e(2)*p534(2)-cz7218(i7,i1,i2)%e(3
     & )*p534(3)
      cvqd=cz7218(i7,i1,i2)%e(0)*p6(0)-cz7218(i7,i1,i2)%e(1)*p6(
     & 1)-cz7218(i7,i1,i2)%e(2)*p6(2)-cz7218(i7,i1,i2)%e(3)*p6(3
     & )
      cauxa=-cz7218(i7,i1,i2)%ek0*quqd+p534k0*cvqd+p6k0*cvqu
      cauxb=-cz7218(i7,i1,i2)%ek0*p6(2)+p6k0*cz7218(i7,i1,i2)%e(
     & 2)
      r6_7218(i7,i1,i2)%a(1)=zcr(id6)*(cauxa+ceps_0)
      r6_7218(i7,i1,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_7218(i7,i1,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      r6_7218(i7,i1,i2)%b(2)=zcr(id6)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p534,qd=p6,v=cf7218(i7,i1,i2)%e,a=r6_7218(i7,i1,i2)%a,b=r6_7
* 218(i7,i1,i2)%b,cr=fcr(id6),cl=fcl(id6),nsum=1
      ceps_0=-cf7218(i7,i1,i2)%ek0*(p534(2)*p6(3)-p6(2)*p534(3))
     & +p534k0*(cf7218(i7,i1,i2)%e(2)*p6(3)-p6(2)*cf7218(i7,i1,i
     & 2)%e(3))-p6k0*(cf7218(i7,i1,i2)%e(2)*p534(3)-p534(2)*cf72
     & 18(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf7218(i7,i1,i2)%e(3)*p6k0+p6(3)*cf7218(i7,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf7218(i7,i1,i2)%e(0)*p534(0)-cf7218(i7,i1,i2)%e(1)*p
     & 534(1)-cf7218(i7,i1,i2)%e(2)*p534(2)-cf7218(i7,i1,i2)%e(3
     & )*p534(3)
      cvqd=cf7218(i7,i1,i2)%e(0)*p6(0)-cf7218(i7,i1,i2)%e(1)*p6(
     & 1)-cf7218(i7,i1,i2)%e(2)*p6(2)-cf7218(i7,i1,i2)%e(3)*p6(3
     & )
      cauxa=-cf7218(i7,i1,i2)%ek0*quqd+p534k0*cvqd+p6k0*cvqu
      cauxb=-cf7218(i7,i1,i2)%ek0*p6(2)+p6k0*cf7218(i7,i1,i2)%e(
     & 2)
      r6_7218(i7,i1,i2)%a(1)=r6_7218(i7,i1,i2)%a(1)+fcr(id6)*(ca
     & uxa+ceps_0)
      r6_7218(i7,i1,i2)%a(2)=r6_7218(i7,i1,i2)%a(2)+fcl(id6)*(ca
     & uxa-ceps_0)
      r6_7218(i7,i1,i2)%b(1)=r6_7218(i7,i1,i2)%b(1)+fcl(id6)*(ca
     & uxb-ceps_2)
      r6_7218(i7,i1,i2)%b(2)=r6_7218(i7,i1,i2)%b(2)+fcr(id6)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id6).ne.1.or.ilept(id3).ne.1) then
  
* quqd -- p=p5234,q=p6
      quqd=p5234(0)*p6(0)-p5234(1)*p6(1)-p5234(2)*p6(2)-p5234(3)
     & *p6(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p5234,qd=p6,v=cz718(i7,i2)%e,a=r6_718(i7,i2)%a,b=r6_718(i7,i
* 2)%b,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cz718(i7,i2)%ek0*(p5234(2)*p6(3)-p6(2)*p5234(3))+p
     & 5234k0*(cz718(i7,i2)%e(2)*p6(3)-p6(2)*cz718(i7,i2)%e(3))-
     & p6k0*(cz718(i7,i2)%e(2)*p5234(3)-p5234(2)*cz718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz718(i7,i2)%e(3)*p6k0+p6(3)*cz718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz718(i7,i2)%e(0)*p5234(0)-cz718(i7,i2)%e(1)*p5234(1)
     & -cz718(i7,i2)%e(2)*p5234(2)-cz718(i7,i2)%e(3)*p5234(3)
      cvqd=cz718(i7,i2)%e(0)*p6(0)-cz718(i7,i2)%e(1)*p6(1)-cz718
     & (i7,i2)%e(2)*p6(2)-cz718(i7,i2)%e(3)*p6(3)
      cauxa=-cz718(i7,i2)%ek0*quqd+p5234k0*cvqd+p6k0*cvqu
      cauxb=-cz718(i7,i2)%ek0*p6(2)+p6k0*cz718(i7,i2)%e(2)
      r6_718(i7,i2)%a(1)=zcr(id5)*(cauxa+ceps_0)
      r6_718(i7,i2)%a(2)=zcl(id5)*(cauxa-ceps_0)
      r6_718(i7,i2)%b(1)=zcl(id5)*(cauxb-ceps_2)
      r6_718(i7,i2)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p5234,qd=p6,v=cf718(i7,i2)%e,a=r6_718(i7,i2)%a,b=r6_718(i7,i
* 2)%b,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cf718(i7,i2)%ek0*(p5234(2)*p6(3)-p6(2)*p5234(3))+p
     & 5234k0*(cf718(i7,i2)%e(2)*p6(3)-p6(2)*cf718(i7,i2)%e(3))-
     & p6k0*(cf718(i7,i2)%e(2)*p5234(3)-p5234(2)*cf718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf718(i7,i2)%e(3)*p6k0+p6(3)*cf718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf718(i7,i2)%e(0)*p5234(0)-cf718(i7,i2)%e(1)*p5234(1)
     & -cf718(i7,i2)%e(2)*p5234(2)-cf718(i7,i2)%e(3)*p5234(3)
      cvqd=cf718(i7,i2)%e(0)*p6(0)-cf718(i7,i2)%e(1)*p6(1)-cf718
     & (i7,i2)%e(2)*p6(2)-cf718(i7,i2)%e(3)*p6(3)
      cauxa=-cf718(i7,i2)%ek0*quqd+p5234k0*cvqd+p6k0*cvqu
      cauxb=-cf718(i7,i2)%ek0*p6(2)+p6k0*cf718(i7,i2)%e(2)
      r6_718(i7,i2)%a(1)=r6_718(i7,i2)%a(1)+fcr(id5)*(cauxa+ceps
     & _0)
      r6_718(i7,i2)%a(2)=r6_718(i7,i2)%a(2)+fcl(id5)*(cauxa-ceps
     & _0)
      r6_718(i7,i2)%b(1)=r6_718(i7,i2)%b(1)+fcl(id5)*(cauxb-ceps
     & _2)
      r6_718(i7,i2)%b(2)=r6_718(i7,i2)%b(2)+fcr(id5)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p578(m)=p5(m)+p7(m)+p8(m)
       enddo
* pk0 -- p=p578
      p578k0=p578(0)-p578(1)
* p.q -- p.q=s578,p=p578,q=p578,bef=,aft=
      s578=(p578(0)*p578(0)-p578(1)*p578(1)-p578(2)*p578(2)-p578
     & (3)*p578(3))
      f578=s578*p578k0
  
       do m=0,3
        p634(m)= -p6(m)-p3(m)-p4(m)
       enddo
* pk0 -- p=p634
      p634k0=p634(0)-p634(1)
* p.q -- p.q=s634,p=p634,q=p634,bef=,aft=
      s634=(p634(0)*p634(0)-p634(1)*p634(1)-p634(2)*p634(2)-p634
     & (3)*p634(3))
      f634=s634*p634k0
  
      if (ilept(id5).ne.1.or.ilept(id3).ne.1) then
      quqd=(s578-s78)/2d0
      ccr=zcr(id5)/(f578)
      ccl=zcl(id5)/(f578)
      do i5=1,2
* TL0 -- qu=p5,qd=p578,v=cz78(i5)%e,a=l5_78(i5)%a,c=l5_78(i5)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p5(2)*p578(3)-p578(2)*p5(3))+p5k0*(c
     & z78(i5)%e(2)*p578(3)-p578(2)*cz78(i5)%e(3))-p578k0*(cz78(
     & i5)%e(2)*p5(3)-p5(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p5k0+p5(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz78(i5)%e(0)*p5(0)-cz78(i5)%e(1)*p5(1)-cz78(i5)%e(2)
     & *p5(2)-cz78(i5)%e(3)*p5(3)
      cvqd=cz78(i5)%e(0)*p578(0)-cz78(i5)%e(1)*p578(1)-cz78(i5)%
     & e(2)*p578(2)-cz78(i5)%e(3)*p578(3)
      cauxa=-cz78(i5)%ek0*quqd+p5k0*cvqd+p578k0*cvqu
      cauxc=+cz78(i5)%ek0*p5(2)-p5k0*cz78(i5)%e(2)
      l5_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      l5_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      l5_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      l5_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      ccr=fcr(id5)/(f578)
      ccl=fcl(id5)/(f578)
      do i5=1,2
* TL0 -- qu=p5,qd=p578,v=cf78(i5)%e,a=l5_78(i5)%a,c=l5_78(i5)%c,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p5(2)*p578(3)-p578(2)*p5(3))+p5k0*(c
     & f78(i5)%e(2)*p578(3)-p578(2)*cf78(i5)%e(3))-p578k0*(cf78(
     & i5)%e(2)*p5(3)-p5(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p5k0+p5(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf78(i5)%e(0)*p5(0)-cf78(i5)%e(1)*p5(1)-cf78(i5)%e(2)
     & *p5(2)-cf78(i5)%e(3)*p5(3)
      cvqd=cf78(i5)%e(0)*p578(0)-cf78(i5)%e(1)*p578(1)-cf78(i5)%
     & e(2)*p578(2)-cf78(i5)%e(3)*p578(3)
      cauxa=-cf78(i5)%ek0*quqd+p5k0*cvqd+p578k0*cvqu
      cauxc=+cf78(i5)%ek0*p5(2)-p5k0*cf78(i5)%e(2)
      l5_78(i5)%a(1)=l5_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      l5_78(i5)%a(2)=l5_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      l5_78(i5)%c(1)=l5_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      l5_78(i5)%c(2)=l5_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      endif
  
      if (ilept(id6).ne.1.or.ilept(id7).ne.1) then
      quqd=(-s634+s34)/2d0
      do i7=1,2
* TR0 -- qu=p634,qd=p6,v=cz34(i7)%e,a=r6_34(i7)%a,b=r6_34(i7)%b,cr=zcr(i
* d6),cl=zcl(id6),nsum=0
      ceps_0=-cz34(i7)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k0*
     & (cz34(i7)%e(2)*p6(3)-p6(2)*cz34(i7)%e(3))-p6k0*(cz34(i7)%
     & e(2)*p634(3)-p634(2)*cz34(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i7)%e(3)*p6k0+p6(3)*cz34(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i7)%e(0)*p634(0)-cz34(i7)%e(1)*p634(1)-cz34(i7)%
     & e(2)*p634(2)-cz34(i7)%e(3)*p634(3)
      cvqd=cz34(i7)%e(0)*p6(0)-cz34(i7)%e(1)*p6(1)-cz34(i7)%e(2)
     & *p6(2)-cz34(i7)%e(3)*p6(3)
      cauxa=-cz34(i7)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cz34(i7)%ek0*p6(2)+p6k0*cz34(i7)%e(2)
      r6_34(i7)%a(1)=zcr(id6)*(cauxa+ceps_0)
      r6_34(i7)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_34(i7)%b(1)=zcl(id6)*(cauxb-ceps_2)
      r6_34(i7)%b(2)=zcr(id6)*(-cauxb-ceps_2)
      end do
  
      do i7=1,2
* TR0 -- qu=p634,qd=p6,v=cf34(i7)%e,a=r6_34(i7)%a,b=r6_34(i7)%b,cr=fcr(i
* d6),cl=fcl(id6),nsum=1
      ceps_0=-cf34(i7)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k0*
     & (cf34(i7)%e(2)*p6(3)-p6(2)*cf34(i7)%e(3))-p6k0*(cf34(i7)%
     & e(2)*p634(3)-p634(2)*cf34(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf34(i7)%e(3)*p6k0+p6(3)*cf34(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i7)%e(0)*p634(0)-cf34(i7)%e(1)*p634(1)-cf34(i7)%
     & e(2)*p634(2)-cf34(i7)%e(3)*p634(3)
      cvqd=cf34(i7)%e(0)*p6(0)-cf34(i7)%e(1)*p6(1)-cf34(i7)%e(2)
     & *p6(2)-cf34(i7)%e(3)*p6(3)
      cauxa=-cf34(i7)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cf34(i7)%ek0*p6(2)+p6k0*cf34(i7)%e(2)
      r6_34(i7)%a(1)=r6_34(i7)%a(1)+fcr(id6)*(cauxa+ceps_0)
      r6_34(i7)%a(2)=r6_34(i7)%a(2)+fcl(id6)*(cauxa-ceps_0)
      r6_34(i7)%b(1)=r6_34(i7)%b(1)+fcl(id6)*(cauxb-ceps_2)
      r6_34(i7)%b(2)=r6_34(i7)%b(2)+fcr(id6)*(-cauxb-ceps_2)
      end do
      endif
  
  
       do m=0,3
        p5178(m)=p578(m) +p1(m)
       enddo
* pk0 -- p=p5178
      p5178k0=p5178(0)-p5178(1)
* p.q -- p.q=s5178,p=p5178,q=p5178,bef=,aft=
      s5178=(p5178(0)*p5178(0)-p5178(1)*p5178(1)-p5178(2)*p5178(
     & 2)-p5178(3)*p5178(3))
      f5178=s5178*p5178k0
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p5,q=p634
      quqd=p5(0)*p634(0)-p5(1)*p634(1)-p5(2)*p634(2)-p5(3)*p634(
     & 3)
      ccr=zcr(id5)/(f634)
      ccl=zcl(id5)/(f634)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p5,qd=p634,v=cz7128(i7,i1,i2)%e,a=l5_7128(i7,i1,i2)%a,c=l5_7
* 128(i7,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz7128(i7,i1,i2)%ek0*(p5(2)*p634(3)-p634(2)*p5(3))
     & +p5k0*(cz7128(i7,i1,i2)%e(2)*p634(3)-p634(2)*cz7128(i7,i1
     & ,i2)%e(3))-p634k0*(cz7128(i7,i1,i2)%e(2)*p5(3)-p5(2)*cz71
     & 28(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz7128(i7,i1,i2)%e(3)*p5k0+p5(3)*cz7128(i7,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz7128(i7,i1,i2)%e(0)*p5(0)-cz7128(i7,i1,i2)%e(1)*p5(
     & 1)-cz7128(i7,i1,i2)%e(2)*p5(2)-cz7128(i7,i1,i2)%e(3)*p5(3
     & )
      cvqd=cz7128(i7,i1,i2)%e(0)*p634(0)-cz7128(i7,i1,i2)%e(1)*p
     & 634(1)-cz7128(i7,i1,i2)%e(2)*p634(2)-cz7128(i7,i1,i2)%e(3
     & )*p634(3)
      cauxa=-cz7128(i7,i1,i2)%ek0*quqd+p5k0*cvqd+p634k0*cvqu
      cauxc=+cz7128(i7,i1,i2)%ek0*p5(2)-p5k0*cz7128(i7,i1,i2)%e(
     & 2)
      l5_7128(i7,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l5_7128(i7,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_7128(i7,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l5_7128(i7,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id5)/(f634)
      ccl=fcl(id5)/(f634)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p5,qd=p634,v=cf7128(i7,i1,i2)%e,a=l5_7128(i7,i1,i2)%a,c=l5_7
* 128(i7,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf7128(i7,i1,i2)%ek0*(p5(2)*p634(3)-p634(2)*p5(3))
     & +p5k0*(cf7128(i7,i1,i2)%e(2)*p634(3)-p634(2)*cf7128(i7,i1
     & ,i2)%e(3))-p634k0*(cf7128(i7,i1,i2)%e(2)*p5(3)-p5(2)*cf71
     & 28(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf7128(i7,i1,i2)%e(3)*p5k0+p5(3)*cf7128(i7,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf7128(i7,i1,i2)%e(0)*p5(0)-cf7128(i7,i1,i2)%e(1)*p5(
     & 1)-cf7128(i7,i1,i2)%e(2)*p5(2)-cf7128(i7,i1,i2)%e(3)*p5(3
     & )
      cvqd=cf7128(i7,i1,i2)%e(0)*p634(0)-cf7128(i7,i1,i2)%e(1)*p
     & 634(1)-cf7128(i7,i1,i2)%e(2)*p634(2)-cf7128(i7,i1,i2)%e(3
     & )*p634(3)
      cauxa=-cf7128(i7,i1,i2)%ek0*quqd+p5k0*cvqd+p634k0*cvqu
      cauxc=+cf7128(i7,i1,i2)%ek0*p5(2)-p5k0*cf7128(i7,i1,i2)%e(
     & 2)
      l5_7128(i7,i1,i2)%a(1)=l5_7128(i7,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l5_7128(i7,i1,i2)%a(2)=l5_7128(i7,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l5_7128(i7,i1,i2)%c(1)=l5_7128(i7,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l5_7128(i7,i1,i2)%c(2)=l5_7128(i7,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id5).ne.1.or.ilept(id3).ne.1) then
  
* quqd -- p=p5,q=p5178
      quqd=p5(0)*p5178(0)-p5(1)*p5178(1)-p5(2)*p5178(2)-p5(3)*p5
     & 178(3)
      ccr=zcr(id5)/(f5178)
      ccl=zcl(id5)/(f5178)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p5,qd=p5178,v=cz718(i5,i1)%e,a=l5_718(i5,i1)%a,c=l5_718(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz718(i5,i1)%ek0*(p5(2)*p5178(3)-p5178(2)*p5(3))+p
     & 5k0*(cz718(i5,i1)%e(2)*p5178(3)-p5178(2)*cz718(i5,i1)%e(3
     & ))-p5178k0*(cz718(i5,i1)%e(2)*p5(3)-p5(2)*cz718(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz718(i5,i1)%e(3)*p5k0+p5(3)*cz718(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz718(i5,i1)%e(0)*p5(0)-cz718(i5,i1)%e(1)*p5(1)-cz718
     & (i5,i1)%e(2)*p5(2)-cz718(i5,i1)%e(3)*p5(3)
      cvqd=cz718(i5,i1)%e(0)*p5178(0)-cz718(i5,i1)%e(1)*p5178(1)
     & -cz718(i5,i1)%e(2)*p5178(2)-cz718(i5,i1)%e(3)*p5178(3)
      cauxa=-cz718(i5,i1)%ek0*quqd+p5k0*cvqd+p5178k0*cvqu
      cauxc=+cz718(i5,i1)%ek0*p5(2)-p5k0*cz718(i5,i1)%e(2)
      l5_718(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l5_718(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_718(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l5_718(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id5)/(f5178)
      ccl=fcl(id5)/(f5178)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p5,qd=p5178,v=cf718(i5,i1)%e,a=l5_718(i5,i1)%a,c=l5_718(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf718(i5,i1)%ek0*(p5(2)*p5178(3)-p5178(2)*p5(3))+p
     & 5k0*(cf718(i5,i1)%e(2)*p5178(3)-p5178(2)*cf718(i5,i1)%e(3
     & ))-p5178k0*(cf718(i5,i1)%e(2)*p5(3)-p5(2)*cf718(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf718(i5,i1)%e(3)*p5k0+p5(3)*cf718(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf718(i5,i1)%e(0)*p5(0)-cf718(i5,i1)%e(1)*p5(1)-cf718
     & (i5,i1)%e(2)*p5(2)-cf718(i5,i1)%e(3)*p5(3)
      cvqd=cf718(i5,i1)%e(0)*p5178(0)-cf718(i5,i1)%e(1)*p5178(1)
     & -cf718(i5,i1)%e(2)*p5178(2)-cf718(i5,i1)%e(3)*p5178(3)
      cauxa=-cf718(i5,i1)%ek0*quqd+p5k0*cvqd+p5178k0*cvqu
      cauxc=+cf718(i5,i1)%ek0*p5(2)-p5k0*cf718(i5,i1)%e(2)
      l5_718(i5,i1)%a(1)=l5_718(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l5_718(i5,i1)%a(2)=l5_718(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l5_718(i5,i1)%c(1)=l5_718(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l5_718(i5,i1)%c(2)=l5_718(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p578,q=p6
      quqd=p578(0)*p6(0)-p578(1)*p6(1)-p578(2)*p6(2)-p578(3)*p6(
     & 3)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p578,qd=p6,v=cz3124(i3,i1,i2)%e,a=r6_3124(i3,i1,i2)%a,b=r6_3
* 124(i3,i1,i2)%b,cr=zcr(id6),cl=zcl(id6),nsum=0
      ceps_0=-cz3124(i3,i1,i2)%ek0*(p578(2)*p6(3)-p6(2)*p578(3))
     & +p578k0*(cz3124(i3,i1,i2)%e(2)*p6(3)-p6(2)*cz3124(i3,i1,i
     & 2)%e(3))-p6k0*(cz3124(i3,i1,i2)%e(2)*p578(3)-p578(2)*cz31
     & 24(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz3124(i3,i1,i2)%e(3)*p6k0+p6(3)*cz3124(i3,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz3124(i3,i1,i2)%e(0)*p578(0)-cz3124(i3,i1,i2)%e(1)*p
     & 578(1)-cz3124(i3,i1,i2)%e(2)*p578(2)-cz3124(i3,i1,i2)%e(3
     & )*p578(3)
      cvqd=cz3124(i3,i1,i2)%e(0)*p6(0)-cz3124(i3,i1,i2)%e(1)*p6(
     & 1)-cz3124(i3,i1,i2)%e(2)*p6(2)-cz3124(i3,i1,i2)%e(3)*p6(3
     & )
      cauxa=-cz3124(i3,i1,i2)%ek0*quqd+p578k0*cvqd+p6k0*cvqu
      cauxb=-cz3124(i3,i1,i2)%ek0*p6(2)+p6k0*cz3124(i3,i1,i2)%e(
     & 2)
      r6_3124(i3,i1,i2)%a(1)=zcr(id6)*(cauxa+ceps_0)
      r6_3124(i3,i1,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_3124(i3,i1,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      r6_3124(i3,i1,i2)%b(2)=zcr(id6)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p578,qd=p6,v=cf3124(i3,i1,i2)%e,a=r6_3124(i3,i1,i2)%a,b=r6_3
* 124(i3,i1,i2)%b,cr=fcr(id6),cl=fcl(id6),nsum=1
      ceps_0=-cf3124(i3,i1,i2)%ek0*(p578(2)*p6(3)-p6(2)*p578(3))
     & +p578k0*(cf3124(i3,i1,i2)%e(2)*p6(3)-p6(2)*cf3124(i3,i1,i
     & 2)%e(3))-p6k0*(cf3124(i3,i1,i2)%e(2)*p578(3)-p578(2)*cf31
     & 24(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf3124(i3,i1,i2)%e(3)*p6k0+p6(3)*cf3124(i3,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf3124(i3,i1,i2)%e(0)*p578(0)-cf3124(i3,i1,i2)%e(1)*p
     & 578(1)-cf3124(i3,i1,i2)%e(2)*p578(2)-cf3124(i3,i1,i2)%e(3
     & )*p578(3)
      cvqd=cf3124(i3,i1,i2)%e(0)*p6(0)-cf3124(i3,i1,i2)%e(1)*p6(
     & 1)-cf3124(i3,i1,i2)%e(2)*p6(2)-cf3124(i3,i1,i2)%e(3)*p6(3
     & )
      cauxa=-cf3124(i3,i1,i2)%ek0*quqd+p578k0*cvqd+p6k0*cvqu
      cauxb=-cf3124(i3,i1,i2)%ek0*p6(2)+p6k0*cf3124(i3,i1,i2)%e(
     & 2)
      r6_3124(i3,i1,i2)%a(1)=r6_3124(i3,i1,i2)%a(1)+fcr(id6)*(ca
     & uxa+ceps_0)
      r6_3124(i3,i1,i2)%a(2)=r6_3124(i3,i1,i2)%a(2)+fcl(id6)*(ca
     & uxa-ceps_0)
      r6_3124(i3,i1,i2)%b(1)=r6_3124(i3,i1,i2)%b(1)+fcl(id6)*(ca
     & uxb-ceps_2)
      r6_3124(i3,i1,i2)%b(2)=r6_3124(i3,i1,i2)%b(2)+fcr(id6)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id6).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p5178,q=p6
      quqd=p5178(0)*p6(0)-p5178(1)*p6(1)-p5178(2)*p6(2)-p5178(3)
     & *p6(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p5178,qd=p6,v=cz324(i7,i2)%e,a=r6_324(i7,i2)%a,b=r6_324(i7,i
* 2)%b,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cz324(i7,i2)%ek0*(p5178(2)*p6(3)-p6(2)*p5178(3))+p
     & 5178k0*(cz324(i7,i2)%e(2)*p6(3)-p6(2)*cz324(i7,i2)%e(3))-
     & p6k0*(cz324(i7,i2)%e(2)*p5178(3)-p5178(2)*cz324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz324(i7,i2)%e(3)*p6k0+p6(3)*cz324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i7,i2)%e(0)*p5178(0)-cz324(i7,i2)%e(1)*p5178(1)
     & -cz324(i7,i2)%e(2)*p5178(2)-cz324(i7,i2)%e(3)*p5178(3)
      cvqd=cz324(i7,i2)%e(0)*p6(0)-cz324(i7,i2)%e(1)*p6(1)-cz324
     & (i7,i2)%e(2)*p6(2)-cz324(i7,i2)%e(3)*p6(3)
      cauxa=-cz324(i7,i2)%ek0*quqd+p5178k0*cvqd+p6k0*cvqu
      cauxb=-cz324(i7,i2)%ek0*p6(2)+p6k0*cz324(i7,i2)%e(2)
      r6_324(i7,i2)%a(1)=zcr(id5)*(cauxa+ceps_0)
      r6_324(i7,i2)%a(2)=zcl(id5)*(cauxa-ceps_0)
      r6_324(i7,i2)%b(1)=zcl(id5)*(cauxb-ceps_2)
      r6_324(i7,i2)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p5178,qd=p6,v=cf324(i7,i2)%e,a=r6_324(i7,i2)%a,b=r6_324(i7,i
* 2)%b,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cf324(i7,i2)%ek0*(p5178(2)*p6(3)-p6(2)*p5178(3))+p
     & 5178k0*(cf324(i7,i2)%e(2)*p6(3)-p6(2)*cf324(i7,i2)%e(3))-
     & p6k0*(cf324(i7,i2)%e(2)*p5178(3)-p5178(2)*cf324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf324(i7,i2)%e(3)*p6k0+p6(3)*cf324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i7,i2)%e(0)*p5178(0)-cf324(i7,i2)%e(1)*p5178(1)
     & -cf324(i7,i2)%e(2)*p5178(2)-cf324(i7,i2)%e(3)*p5178(3)
      cvqd=cf324(i7,i2)%e(0)*p6(0)-cf324(i7,i2)%e(1)*p6(1)-cf324
     & (i7,i2)%e(2)*p6(2)-cf324(i7,i2)%e(3)*p6(3)
      cauxa=-cf324(i7,i2)%ek0*quqd+p5178k0*cvqd+p6k0*cvqu
      cauxb=-cf324(i7,i2)%ek0*p6(2)+p6k0*cf324(i7,i2)%e(2)
      r6_324(i7,i2)%a(1)=r6_324(i7,i2)%a(1)+fcr(id5)*(cauxa+ceps
     & _0)
      r6_324(i7,i2)%a(2)=r6_324(i7,i2)%a(2)+fcl(id5)*(cauxa-ceps
     & _0)
      r6_324(i7,i2)%b(1)=r6_324(i7,i2)%b(1)+fcl(id5)*(cauxb-ceps
     & _2)
      r6_324(i7,i2)%b(2)=r6_324(i7,i2)%b(2)+fcr(id5)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p5278(m)=p578(m) +p2(m)
       enddo
* pk0 -- p=p5278
      p5278k0=p5278(0)-p5278(1)
* p.q -- p.q=s5278,p=p5278,q=p5278,bef=,aft=
      s5278=(p5278(0)*p5278(0)-p5278(1)*p5278(1)-p5278(2)*p5278(
     & 2)-p5278(3)*p5278(3))
      f5278=s5278*p5278k0
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p5,q=p634
      quqd=p5(0)*p634(0)-p5(1)*p634(1)-p5(2)*p634(2)-p5(3)*p634(
     & 3)
      ccr=zcr(id5)/(f634)
      ccl=zcl(id5)/(f634)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p5,qd=p634,v=cz7218(i7,i1,i2)%e,a=l5_7218(i7,i1,i2)%a,c=l5_7
* 218(i7,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz7218(i7,i1,i2)%ek0*(p5(2)*p634(3)-p634(2)*p5(3))
     & +p5k0*(cz7218(i7,i1,i2)%e(2)*p634(3)-p634(2)*cz7218(i7,i1
     & ,i2)%e(3))-p634k0*(cz7218(i7,i1,i2)%e(2)*p5(3)-p5(2)*cz72
     & 18(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz7218(i7,i1,i2)%e(3)*p5k0+p5(3)*cz7218(i7,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz7218(i7,i1,i2)%e(0)*p5(0)-cz7218(i7,i1,i2)%e(1)*p5(
     & 1)-cz7218(i7,i1,i2)%e(2)*p5(2)-cz7218(i7,i1,i2)%e(3)*p5(3
     & )
      cvqd=cz7218(i7,i1,i2)%e(0)*p634(0)-cz7218(i7,i1,i2)%e(1)*p
     & 634(1)-cz7218(i7,i1,i2)%e(2)*p634(2)-cz7218(i7,i1,i2)%e(3
     & )*p634(3)
      cauxa=-cz7218(i7,i1,i2)%ek0*quqd+p5k0*cvqd+p634k0*cvqu
      cauxc=+cz7218(i7,i1,i2)%ek0*p5(2)-p5k0*cz7218(i7,i1,i2)%e(
     & 2)
      l5_7218(i7,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l5_7218(i7,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_7218(i7,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l5_7218(i7,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id5)/(f634)
      ccl=fcl(id5)/(f634)
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p5,qd=p634,v=cf7218(i7,i1,i2)%e,a=l5_7218(i7,i1,i2)%a,c=l5_7
* 218(i7,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf7218(i7,i1,i2)%ek0*(p5(2)*p634(3)-p634(2)*p5(3))
     & +p5k0*(cf7218(i7,i1,i2)%e(2)*p634(3)-p634(2)*cf7218(i7,i1
     & ,i2)%e(3))-p634k0*(cf7218(i7,i1,i2)%e(2)*p5(3)-p5(2)*cf72
     & 18(i7,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf7218(i7,i1,i2)%e(3)*p5k0+p5(3)*cf7218(i7,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf7218(i7,i1,i2)%e(0)*p5(0)-cf7218(i7,i1,i2)%e(1)*p5(
     & 1)-cf7218(i7,i1,i2)%e(2)*p5(2)-cf7218(i7,i1,i2)%e(3)*p5(3
     & )
      cvqd=cf7218(i7,i1,i2)%e(0)*p634(0)-cf7218(i7,i1,i2)%e(1)*p
     & 634(1)-cf7218(i7,i1,i2)%e(2)*p634(2)-cf7218(i7,i1,i2)%e(3
     & )*p634(3)
      cauxa=-cf7218(i7,i1,i2)%ek0*quqd+p5k0*cvqd+p634k0*cvqu
      cauxc=+cf7218(i7,i1,i2)%ek0*p5(2)-p5k0*cf7218(i7,i1,i2)%e(
     & 2)
      l5_7218(i7,i1,i2)%a(1)=l5_7218(i7,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l5_7218(i7,i1,i2)%a(2)=l5_7218(i7,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l5_7218(i7,i1,i2)%c(1)=l5_7218(i7,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l5_7218(i7,i1,i2)%c(2)=l5_7218(i7,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id5).ne.1.or.ilept(id3).ne.1) then
  
* quqd -- p=p5,q=p5278
      quqd=p5(0)*p5278(0)-p5(1)*p5278(1)-p5(2)*p5278(2)-p5(3)*p5
     & 278(3)
      ccr=zcr(id5)/(f5278)
      ccl=zcl(id5)/(f5278)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p5,qd=p5278,v=cz728(i5,i1)%e,a=l5_728(i5,i1)%a,c=l5_728(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz728(i5,i1)%ek0*(p5(2)*p5278(3)-p5278(2)*p5(3))+p
     & 5k0*(cz728(i5,i1)%e(2)*p5278(3)-p5278(2)*cz728(i5,i1)%e(3
     & ))-p5278k0*(cz728(i5,i1)%e(2)*p5(3)-p5(2)*cz728(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz728(i5,i1)%e(3)*p5k0+p5(3)*cz728(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz728(i5,i1)%e(0)*p5(0)-cz728(i5,i1)%e(1)*p5(1)-cz728
     & (i5,i1)%e(2)*p5(2)-cz728(i5,i1)%e(3)*p5(3)
      cvqd=cz728(i5,i1)%e(0)*p5278(0)-cz728(i5,i1)%e(1)*p5278(1)
     & -cz728(i5,i1)%e(2)*p5278(2)-cz728(i5,i1)%e(3)*p5278(3)
      cauxa=-cz728(i5,i1)%ek0*quqd+p5k0*cvqd+p5278k0*cvqu
      cauxc=+cz728(i5,i1)%ek0*p5(2)-p5k0*cz728(i5,i1)%e(2)
      l5_728(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l5_728(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_728(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l5_728(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id5)/(f5278)
      ccl=fcl(id5)/(f5278)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p5,qd=p5278,v=cf728(i5,i1)%e,a=l5_728(i5,i1)%a,c=l5_728(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf728(i5,i1)%ek0*(p5(2)*p5278(3)-p5278(2)*p5(3))+p
     & 5k0*(cf728(i5,i1)%e(2)*p5278(3)-p5278(2)*cf728(i5,i1)%e(3
     & ))-p5278k0*(cf728(i5,i1)%e(2)*p5(3)-p5(2)*cf728(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf728(i5,i1)%e(3)*p5k0+p5(3)*cf728(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf728(i5,i1)%e(0)*p5(0)-cf728(i5,i1)%e(1)*p5(1)-cf728
     & (i5,i1)%e(2)*p5(2)-cf728(i5,i1)%e(3)*p5(3)
      cvqd=cf728(i5,i1)%e(0)*p5278(0)-cf728(i5,i1)%e(1)*p5278(1)
     & -cf728(i5,i1)%e(2)*p5278(2)-cf728(i5,i1)%e(3)*p5278(3)
      cauxa=-cf728(i5,i1)%ek0*quqd+p5k0*cvqd+p5278k0*cvqu
      cauxc=+cf728(i5,i1)%ek0*p5(2)-p5k0*cf728(i5,i1)%e(2)
      l5_728(i5,i1)%a(1)=l5_728(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l5_728(i5,i1)%a(2)=l5_728(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l5_728(i5,i1)%c(1)=l5_728(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l5_728(i5,i1)%c(2)=l5_728(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p578,q=p6
      quqd=p578(0)*p6(0)-p578(1)*p6(1)-p578(2)*p6(2)-p578(3)*p6(
     & 3)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p578,qd=p6,v=cz3214(i3,i1,i2)%e,a=r6_3214(i3,i1,i2)%a,b=r6_3
* 214(i3,i1,i2)%b,cr=zcr(id6),cl=zcl(id6),nsum=0
      ceps_0=-cz3214(i3,i1,i2)%ek0*(p578(2)*p6(3)-p6(2)*p578(3))
     & +p578k0*(cz3214(i3,i1,i2)%e(2)*p6(3)-p6(2)*cz3214(i3,i1,i
     & 2)%e(3))-p6k0*(cz3214(i3,i1,i2)%e(2)*p578(3)-p578(2)*cz32
     & 14(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz3214(i3,i1,i2)%e(3)*p6k0+p6(3)*cz3214(i3,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz3214(i3,i1,i2)%e(0)*p578(0)-cz3214(i3,i1,i2)%e(1)*p
     & 578(1)-cz3214(i3,i1,i2)%e(2)*p578(2)-cz3214(i3,i1,i2)%e(3
     & )*p578(3)
      cvqd=cz3214(i3,i1,i2)%e(0)*p6(0)-cz3214(i3,i1,i2)%e(1)*p6(
     & 1)-cz3214(i3,i1,i2)%e(2)*p6(2)-cz3214(i3,i1,i2)%e(3)*p6(3
     & )
      cauxa=-cz3214(i3,i1,i2)%ek0*quqd+p578k0*cvqd+p6k0*cvqu
      cauxb=-cz3214(i3,i1,i2)%ek0*p6(2)+p6k0*cz3214(i3,i1,i2)%e(
     & 2)
      r6_3214(i3,i1,i2)%a(1)=zcr(id6)*(cauxa+ceps_0)
      r6_3214(i3,i1,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_3214(i3,i1,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      r6_3214(i3,i1,i2)%b(2)=zcr(id6)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p578,qd=p6,v=cf3214(i3,i1,i2)%e,a=r6_3214(i3,i1,i2)%a,b=r6_3
* 214(i3,i1,i2)%b,cr=fcr(id6),cl=fcl(id6),nsum=1
      ceps_0=-cf3214(i3,i1,i2)%ek0*(p578(2)*p6(3)-p6(2)*p578(3))
     & +p578k0*(cf3214(i3,i1,i2)%e(2)*p6(3)-p6(2)*cf3214(i3,i1,i
     & 2)%e(3))-p6k0*(cf3214(i3,i1,i2)%e(2)*p578(3)-p578(2)*cf32
     & 14(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf3214(i3,i1,i2)%e(3)*p6k0+p6(3)*cf3214(i3,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf3214(i3,i1,i2)%e(0)*p578(0)-cf3214(i3,i1,i2)%e(1)*p
     & 578(1)-cf3214(i3,i1,i2)%e(2)*p578(2)-cf3214(i3,i1,i2)%e(3
     & )*p578(3)
      cvqd=cf3214(i3,i1,i2)%e(0)*p6(0)-cf3214(i3,i1,i2)%e(1)*p6(
     & 1)-cf3214(i3,i1,i2)%e(2)*p6(2)-cf3214(i3,i1,i2)%e(3)*p6(3
     & )
      cauxa=-cf3214(i3,i1,i2)%ek0*quqd+p578k0*cvqd+p6k0*cvqu
      cauxb=-cf3214(i3,i1,i2)%ek0*p6(2)+p6k0*cf3214(i3,i1,i2)%e(
     & 2)
      r6_3214(i3,i1,i2)%a(1)=r6_3214(i3,i1,i2)%a(1)+fcr(id6)*(ca
     & uxa+ceps_0)
      r6_3214(i3,i1,i2)%a(2)=r6_3214(i3,i1,i2)%a(2)+fcl(id6)*(ca
     & uxa-ceps_0)
      r6_3214(i3,i1,i2)%b(1)=r6_3214(i3,i1,i2)%b(1)+fcl(id6)*(ca
     & uxb-ceps_2)
      r6_3214(i3,i1,i2)%b(2)=r6_3214(i3,i1,i2)%b(2)+fcr(id6)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id6).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p5278,q=p6
      quqd=p5278(0)*p6(0)-p5278(1)*p6(1)-p5278(2)*p6(2)-p5278(3)
     & *p6(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p5278,qd=p6,v=cz314(i7,i2)%e,a=r6_314(i7,i2)%a,b=r6_314(i7,i
* 2)%b,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cz314(i7,i2)%ek0*(p5278(2)*p6(3)-p6(2)*p5278(3))+p
     & 5278k0*(cz314(i7,i2)%e(2)*p6(3)-p6(2)*cz314(i7,i2)%e(3))-
     & p6k0*(cz314(i7,i2)%e(2)*p5278(3)-p5278(2)*cz314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz314(i7,i2)%e(3)*p6k0+p6(3)*cz314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i7,i2)%e(0)*p5278(0)-cz314(i7,i2)%e(1)*p5278(1)
     & -cz314(i7,i2)%e(2)*p5278(2)-cz314(i7,i2)%e(3)*p5278(3)
      cvqd=cz314(i7,i2)%e(0)*p6(0)-cz314(i7,i2)%e(1)*p6(1)-cz314
     & (i7,i2)%e(2)*p6(2)-cz314(i7,i2)%e(3)*p6(3)
      cauxa=-cz314(i7,i2)%ek0*quqd+p5278k0*cvqd+p6k0*cvqu
      cauxb=-cz314(i7,i2)%ek0*p6(2)+p6k0*cz314(i7,i2)%e(2)
      r6_314(i7,i2)%a(1)=zcr(id5)*(cauxa+ceps_0)
      r6_314(i7,i2)%a(2)=zcl(id5)*(cauxa-ceps_0)
      r6_314(i7,i2)%b(1)=zcl(id5)*(cauxb-ceps_2)
      r6_314(i7,i2)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p5278,qd=p6,v=cf314(i7,i2)%e,a=r6_314(i7,i2)%a,b=r6_314(i7,i
* 2)%b,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cf314(i7,i2)%ek0*(p5278(2)*p6(3)-p6(2)*p5278(3))+p
     & 5278k0*(cf314(i7,i2)%e(2)*p6(3)-p6(2)*cf314(i7,i2)%e(3))-
     & p6k0*(cf314(i7,i2)%e(2)*p5278(3)-p5278(2)*cf314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf314(i7,i2)%e(3)*p6k0+p6(3)*cf314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i7,i2)%e(0)*p5278(0)-cf314(i7,i2)%e(1)*p5278(1)
     & -cf314(i7,i2)%e(2)*p5278(2)-cf314(i7,i2)%e(3)*p5278(3)
      cvqd=cf314(i7,i2)%e(0)*p6(0)-cf314(i7,i2)%e(1)*p6(1)-cf314
     & (i7,i2)%e(2)*p6(2)-cf314(i7,i2)%e(3)*p6(3)
      cauxa=-cf314(i7,i2)%ek0*quqd+p5278k0*cvqd+p6k0*cvqu
      cauxb=-cf314(i7,i2)%ek0*p6(2)+p6k0*cf314(i7,i2)%e(2)
      r6_314(i7,i2)%a(1)=r6_314(i7,i2)%a(1)+fcr(id5)*(cauxa+ceps
     & _0)
      r6_314(i7,i2)%a(2)=r6_314(i7,i2)%a(2)+fcl(id5)*(cauxa-ceps
     & _0)
      r6_314(i7,i2)%b(1)=r6_314(i7,i2)%b(1)+fcl(id5)*(cauxb-ceps
     & _2)
      r6_314(i7,i2)%b(2)=r6_314(i7,i2)%b(2)+fcr(id5)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p734(m)=p7(m)+p3(m)+p4(m)
       enddo
* pk0 -- p=p734
      p734k0=p734(0)-p734(1)
* p.q -- p.q=s734,p=p734,q=p734,bef=,aft=
      s734=(p734(0)*p734(0)-p734(1)*p734(1)-p734(2)*p734(2)-p734
     & (3)*p734(3))
      f734=s734*p734k0
  
       do m=0,3
        p856(m)= -p8(m)-p5(m)-p6(m)
       enddo
* pk0 -- p=p856
      p856k0=p856(0)-p856(1)
* p.q -- p.q=s856,p=p856,q=p856,bef=,aft=
      s856=(p856(0)*p856(0)-p856(1)*p856(1)-p856(2)*p856(2)-p856
     & (3)*p856(3))
      f856=s856*p856k0
  
      if (ilept(id7).ne.1.or.ilept(id5).ne.1) then
      quqd=(s734-s34)/2d0
      ccr=zcr(id7)/(f734)
      ccl=zcl(id7)/(f734)
      do i5=1,2
* TL0 -- qu=p7,qd=p734,v=cz34(i5)%e,a=l7_34(i5)%a,c=l7_34(i5)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*(c
     & z34(i5)%e(2)*p734(3)-p734(2)*cz34(i5)%e(3))-p734k0*(cz34(
     & i5)%e(2)*p7(3)-p7(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p7k0+p7(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i5)%e(0)*p7(0)-cz34(i5)%e(1)*p7(1)-cz34(i5)%e(2)
     & *p7(2)-cz34(i5)%e(3)*p7(3)
      cvqd=cz34(i5)%e(0)*p734(0)-cz34(i5)%e(1)*p734(1)-cz34(i5)%
     & e(2)*p734(2)-cz34(i5)%e(3)*p734(3)
      cauxa=-cz34(i5)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cz34(i5)%ek0*p7(2)-p7k0*cz34(i5)%e(2)
      l7_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      l7_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      l7_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      l7_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      ccr=fcr(id7)/(f734)
      ccl=fcl(id7)/(f734)
      do i5=1,2
* TL0 -- qu=p7,qd=p734,v=cf34(i5)%e,a=l7_34(i5)%a,c=l7_34(i5)%c,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*(c
     & f34(i5)%e(2)*p734(3)-p734(2)*cf34(i5)%e(3))-p734k0*(cf34(
     & i5)%e(2)*p7(3)-p7(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p7k0+p7(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf34(i5)%e(0)*p7(0)-cf34(i5)%e(1)*p7(1)-cf34(i5)%e(2)
     & *p7(2)-cf34(i5)%e(3)*p7(3)
      cvqd=cf34(i5)%e(0)*p734(0)-cf34(i5)%e(1)*p734(1)-cf34(i5)%
     & e(2)*p734(2)-cf34(i5)%e(3)*p734(3)
      cauxa=-cf34(i5)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cf34(i5)%ek0*p7(2)-p7k0*cf34(i5)%e(2)
      l7_34(i5)%a(1)=l7_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      l7_34(i5)%a(2)=l7_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      l7_34(i5)%c(1)=l7_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      l7_34(i5)%c(2)=l7_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      endif
  
      if (ilept(id8).ne.1.or.ilept(id3).ne.1) then
      quqd=(-s856+s56)/2d0
      do i7=1,2
* TR0 -- qu=p856,qd=p8,v=cz56(i7)%e,a=r8_56(i7)%a,b=r8_56(i7)%b,cr=zcr(i
* d8),cl=zcl(id8),nsum=0
      ceps_0=-cz56(i7)%ek0*(p856(2)*p8(3)-p8(2)*p856(3))+p856k0*
     & (cz56(i7)%e(2)*p8(3)-p8(2)*cz56(i7)%e(3))-p8k0*(cz56(i7)%
     & e(2)*p856(3)-p856(2)*cz56(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz56(i7)%e(3)*p8k0+p8(3)*cz56(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i7)%e(0)*p856(0)-cz56(i7)%e(1)*p856(1)-cz56(i7)%
     & e(2)*p856(2)-cz56(i7)%e(3)*p856(3)
      cvqd=cz56(i7)%e(0)*p8(0)-cz56(i7)%e(1)*p8(1)-cz56(i7)%e(2)
     & *p8(2)-cz56(i7)%e(3)*p8(3)
      cauxa=-cz56(i7)%ek0*quqd+p856k0*cvqd+p8k0*cvqu
      cauxb=-cz56(i7)%ek0*p8(2)+p8k0*cz56(i7)%e(2)
      r8_56(i7)%a(1)=zcr(id8)*(cauxa+ceps_0)
      r8_56(i7)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_56(i7)%b(1)=zcl(id8)*(cauxb-ceps_2)
      r8_56(i7)%b(2)=zcr(id8)*(-cauxb-ceps_2)
      end do
  
      do i7=1,2
* TR0 -- qu=p856,qd=p8,v=cf56(i7)%e,a=r8_56(i7)%a,b=r8_56(i7)%b,cr=fcr(i
* d8),cl=fcl(id8),nsum=1
      ceps_0=-cf56(i7)%ek0*(p856(2)*p8(3)-p8(2)*p856(3))+p856k0*
     & (cf56(i7)%e(2)*p8(3)-p8(2)*cf56(i7)%e(3))-p8k0*(cf56(i7)%
     & e(2)*p856(3)-p856(2)*cf56(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf56(i7)%e(3)*p8k0+p8(3)*cf56(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i7)%e(0)*p856(0)-cf56(i7)%e(1)*p856(1)-cf56(i7)%
     & e(2)*p856(2)-cf56(i7)%e(3)*p856(3)
      cvqd=cf56(i7)%e(0)*p8(0)-cf56(i7)%e(1)*p8(1)-cf56(i7)%e(2)
     & *p8(2)-cf56(i7)%e(3)*p8(3)
      cauxa=-cf56(i7)%ek0*quqd+p856k0*cvqd+p8k0*cvqu
      cauxb=-cf56(i7)%ek0*p8(2)+p8k0*cf56(i7)%e(2)
      r8_56(i7)%a(1)=r8_56(i7)%a(1)+fcr(id8)*(cauxa+ceps_0)
      r8_56(i7)%a(2)=r8_56(i7)%a(2)+fcl(id8)*(cauxa-ceps_0)
      r8_56(i7)%b(1)=r8_56(i7)%b(1)+fcl(id8)*(cauxb-ceps_2)
      r8_56(i7)%b(2)=r8_56(i7)%b(2)+fcr(id8)*(-cauxb-ceps_2)
      end do
      endif
  
  
       do m=0,3
        p7134(m)=p734(m) +p1(m)
       enddo
* pk0 -- p=p7134
      p7134k0=p7134(0)-p7134(1)
* p.q -- p.q=s7134,p=p7134,q=p7134,bef=,aft=
      s7134=(p7134(0)*p7134(0)-p7134(1)*p7134(1)-p7134(2)*p7134(
     & 2)-p7134(3)*p7134(3))
      f7134=s7134*p7134k0
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p7,q=p856
      quqd=p7(0)*p856(0)-p7(1)*p856(1)-p7(2)*p856(2)-p7(3)*p856(
     & 3)
      ccr=zcr(id7)/(f856)
      ccl=zcl(id7)/(f856)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p7,qd=p856,v=cz3124(i3,i1,i2)%e,a=l7_3124(i3,i1,i2)%a,c=l7_3
* 124(i3,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz3124(i3,i1,i2)%ek0*(p7(2)*p856(3)-p856(2)*p7(3))
     & +p7k0*(cz3124(i3,i1,i2)%e(2)*p856(3)-p856(2)*cz3124(i3,i1
     & ,i2)%e(3))-p856k0*(cz3124(i3,i1,i2)%e(2)*p7(3)-p7(2)*cz31
     & 24(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz3124(i3,i1,i2)%e(3)*p7k0+p7(3)*cz3124(i3,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz3124(i3,i1,i2)%e(0)*p7(0)-cz3124(i3,i1,i2)%e(1)*p7(
     & 1)-cz3124(i3,i1,i2)%e(2)*p7(2)-cz3124(i3,i1,i2)%e(3)*p7(3
     & )
      cvqd=cz3124(i3,i1,i2)%e(0)*p856(0)-cz3124(i3,i1,i2)%e(1)*p
     & 856(1)-cz3124(i3,i1,i2)%e(2)*p856(2)-cz3124(i3,i1,i2)%e(3
     & )*p856(3)
      cauxa=-cz3124(i3,i1,i2)%ek0*quqd+p7k0*cvqd+p856k0*cvqu
      cauxc=+cz3124(i3,i1,i2)%ek0*p7(2)-p7k0*cz3124(i3,i1,i2)%e(
     & 2)
      l7_3124(i3,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l7_3124(i3,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_3124(i3,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l7_3124(i3,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id7)/(f856)
      ccl=fcl(id7)/(f856)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p7,qd=p856,v=cf3124(i3,i1,i2)%e,a=l7_3124(i3,i1,i2)%a,c=l7_3
* 124(i3,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf3124(i3,i1,i2)%ek0*(p7(2)*p856(3)-p856(2)*p7(3))
     & +p7k0*(cf3124(i3,i1,i2)%e(2)*p856(3)-p856(2)*cf3124(i3,i1
     & ,i2)%e(3))-p856k0*(cf3124(i3,i1,i2)%e(2)*p7(3)-p7(2)*cf31
     & 24(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf3124(i3,i1,i2)%e(3)*p7k0+p7(3)*cf3124(i3,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf3124(i3,i1,i2)%e(0)*p7(0)-cf3124(i3,i1,i2)%e(1)*p7(
     & 1)-cf3124(i3,i1,i2)%e(2)*p7(2)-cf3124(i3,i1,i2)%e(3)*p7(3
     & )
      cvqd=cf3124(i3,i1,i2)%e(0)*p856(0)-cf3124(i3,i1,i2)%e(1)*p
     & 856(1)-cf3124(i3,i1,i2)%e(2)*p856(2)-cf3124(i3,i1,i2)%e(3
     & )*p856(3)
      cauxa=-cf3124(i3,i1,i2)%ek0*quqd+p7k0*cvqd+p856k0*cvqu
      cauxc=+cf3124(i3,i1,i2)%ek0*p7(2)-p7k0*cf3124(i3,i1,i2)%e(
     & 2)
      l7_3124(i3,i1,i2)%a(1)=l7_3124(i3,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l7_3124(i3,i1,i2)%a(2)=l7_3124(i3,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l7_3124(i3,i1,i2)%c(1)=l7_3124(i3,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l7_3124(i3,i1,i2)%c(2)=l7_3124(i3,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id7).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p7,q=p7134
      quqd=p7(0)*p7134(0)-p7(1)*p7134(1)-p7(2)*p7134(2)-p7(3)*p7
     & 134(3)
      ccr=zcr(id7)/(f7134)
      ccl=zcl(id7)/(f7134)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p7,qd=p7134,v=cz314(i5,i1)%e,a=l7_314(i5,i1)%a,c=l7_314(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz314(i5,i1)%ek0*(p7(2)*p7134(3)-p7134(2)*p7(3))+p
     & 7k0*(cz314(i5,i1)%e(2)*p7134(3)-p7134(2)*cz314(i5,i1)%e(3
     & ))-p7134k0*(cz314(i5,i1)%e(2)*p7(3)-p7(2)*cz314(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i5,i1)%e(3)*p7k0+p7(3)*cz314(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz314(i5,i1)%e(0)*p7(0)-cz314(i5,i1)%e(1)*p7(1)-cz314
     & (i5,i1)%e(2)*p7(2)-cz314(i5,i1)%e(3)*p7(3)
      cvqd=cz314(i5,i1)%e(0)*p7134(0)-cz314(i5,i1)%e(1)*p7134(1)
     & -cz314(i5,i1)%e(2)*p7134(2)-cz314(i5,i1)%e(3)*p7134(3)
      cauxa=-cz314(i5,i1)%ek0*quqd+p7k0*cvqd+p7134k0*cvqu
      cauxc=+cz314(i5,i1)%ek0*p7(2)-p7k0*cz314(i5,i1)%e(2)
      l7_314(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l7_314(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_314(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l7_314(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id7)/(f7134)
      ccl=fcl(id7)/(f7134)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p7,qd=p7134,v=cf314(i5,i1)%e,a=l7_314(i5,i1)%a,c=l7_314(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf314(i5,i1)%ek0*(p7(2)*p7134(3)-p7134(2)*p7(3))+p
     & 7k0*(cf314(i5,i1)%e(2)*p7134(3)-p7134(2)*cf314(i5,i1)%e(3
     & ))-p7134k0*(cf314(i5,i1)%e(2)*p7(3)-p7(2)*cf314(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i5,i1)%e(3)*p7k0+p7(3)*cf314(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf314(i5,i1)%e(0)*p7(0)-cf314(i5,i1)%e(1)*p7(1)-cf314
     & (i5,i1)%e(2)*p7(2)-cf314(i5,i1)%e(3)*p7(3)
      cvqd=cf314(i5,i1)%e(0)*p7134(0)-cf314(i5,i1)%e(1)*p7134(1)
     & -cf314(i5,i1)%e(2)*p7134(2)-cf314(i5,i1)%e(3)*p7134(3)
      cauxa=-cf314(i5,i1)%ek0*quqd+p7k0*cvqd+p7134k0*cvqu
      cauxc=+cf314(i5,i1)%ek0*p7(2)-p7k0*cf314(i5,i1)%e(2)
      l7_314(i5,i1)%a(1)=l7_314(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l7_314(i5,i1)%a(2)=l7_314(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l7_314(i5,i1)%c(1)=l7_314(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l7_314(i5,i1)%c(2)=l7_314(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p734,q=p8
      quqd=p734(0)*p8(0)-p734(1)*p8(1)-p734(2)*p8(2)-p734(3)*p8(
     & 3)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p734,qd=p8,v=cz5126(i5,i1,i2)%e,a=r8_5126(i5,i1,i2)%a,b=r8_5
* 126(i5,i1,i2)%b,cr=zcr(id8),cl=zcl(id8),nsum=0
      ceps_0=-cz5126(i5,i1,i2)%ek0*(p734(2)*p8(3)-p8(2)*p734(3))
     & +p734k0*(cz5126(i5,i1,i2)%e(2)*p8(3)-p8(2)*cz5126(i5,i1,i
     & 2)%e(3))-p8k0*(cz5126(i5,i1,i2)%e(2)*p734(3)-p734(2)*cz51
     & 26(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz5126(i5,i1,i2)%e(3)*p8k0+p8(3)*cz5126(i5,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz5126(i5,i1,i2)%e(0)*p734(0)-cz5126(i5,i1,i2)%e(1)*p
     & 734(1)-cz5126(i5,i1,i2)%e(2)*p734(2)-cz5126(i5,i1,i2)%e(3
     & )*p734(3)
      cvqd=cz5126(i5,i1,i2)%e(0)*p8(0)-cz5126(i5,i1,i2)%e(1)*p8(
     & 1)-cz5126(i5,i1,i2)%e(2)*p8(2)-cz5126(i5,i1,i2)%e(3)*p8(3
     & )
      cauxa=-cz5126(i5,i1,i2)%ek0*quqd+p734k0*cvqd+p8k0*cvqu
      cauxb=-cz5126(i5,i1,i2)%ek0*p8(2)+p8k0*cz5126(i5,i1,i2)%e(
     & 2)
      r8_5126(i5,i1,i2)%a(1)=zcr(id8)*(cauxa+ceps_0)
      r8_5126(i5,i1,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_5126(i5,i1,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      r8_5126(i5,i1,i2)%b(2)=zcr(id8)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p734,qd=p8,v=cf5126(i5,i1,i2)%e,a=r8_5126(i5,i1,i2)%a,b=r8_5
* 126(i5,i1,i2)%b,cr=fcr(id8),cl=fcl(id8),nsum=1
      ceps_0=-cf5126(i5,i1,i2)%ek0*(p734(2)*p8(3)-p8(2)*p734(3))
     & +p734k0*(cf5126(i5,i1,i2)%e(2)*p8(3)-p8(2)*cf5126(i5,i1,i
     & 2)%e(3))-p8k0*(cf5126(i5,i1,i2)%e(2)*p734(3)-p734(2)*cf51
     & 26(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf5126(i5,i1,i2)%e(3)*p8k0+p8(3)*cf5126(i5,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf5126(i5,i1,i2)%e(0)*p734(0)-cf5126(i5,i1,i2)%e(1)*p
     & 734(1)-cf5126(i5,i1,i2)%e(2)*p734(2)-cf5126(i5,i1,i2)%e(3
     & )*p734(3)
      cvqd=cf5126(i5,i1,i2)%e(0)*p8(0)-cf5126(i5,i1,i2)%e(1)*p8(
     & 1)-cf5126(i5,i1,i2)%e(2)*p8(2)-cf5126(i5,i1,i2)%e(3)*p8(3
     & )
      cauxa=-cf5126(i5,i1,i2)%ek0*quqd+p734k0*cvqd+p8k0*cvqu
      cauxb=-cf5126(i5,i1,i2)%ek0*p8(2)+p8k0*cf5126(i5,i1,i2)%e(
     & 2)
      r8_5126(i5,i1,i2)%a(1)=r8_5126(i5,i1,i2)%a(1)+fcr(id8)*(ca
     & uxa+ceps_0)
      r8_5126(i5,i1,i2)%a(2)=r8_5126(i5,i1,i2)%a(2)+fcl(id8)*(ca
     & uxa-ceps_0)
      r8_5126(i5,i1,i2)%b(1)=r8_5126(i5,i1,i2)%b(1)+fcl(id8)*(ca
     & uxb-ceps_2)
      r8_5126(i5,i1,i2)%b(2)=r8_5126(i5,i1,i2)%b(2)+fcr(id8)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id8).ne.1.or.ilept(id3).ne.1) then
  
* quqd -- p=p7134,q=p8
      quqd=p7134(0)*p8(0)-p7134(1)*p8(1)-p7134(2)*p8(2)-p7134(3)
     & *p8(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p7134,qd=p8,v=cz526(i7,i2)%e,a=r8_526(i7,i2)%a,b=r8_526(i7,i
* 2)%b,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cz526(i7,i2)%ek0*(p7134(2)*p8(3)-p8(2)*p7134(3))+p
     & 7134k0*(cz526(i7,i2)%e(2)*p8(3)-p8(2)*cz526(i7,i2)%e(3))-
     & p8k0*(cz526(i7,i2)%e(2)*p7134(3)-p7134(2)*cz526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz526(i7,i2)%e(3)*p8k0+p8(3)*cz526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz526(i7,i2)%e(0)*p7134(0)-cz526(i7,i2)%e(1)*p7134(1)
     & -cz526(i7,i2)%e(2)*p7134(2)-cz526(i7,i2)%e(3)*p7134(3)
      cvqd=cz526(i7,i2)%e(0)*p8(0)-cz526(i7,i2)%e(1)*p8(1)-cz526
     & (i7,i2)%e(2)*p8(2)-cz526(i7,i2)%e(3)*p8(3)
      cauxa=-cz526(i7,i2)%ek0*quqd+p7134k0*cvqd+p8k0*cvqu
      cauxb=-cz526(i7,i2)%ek0*p8(2)+p8k0*cz526(i7,i2)%e(2)
      r8_526(i7,i2)%a(1)=zcr(id7)*(cauxa+ceps_0)
      r8_526(i7,i2)%a(2)=zcl(id7)*(cauxa-ceps_0)
      r8_526(i7,i2)%b(1)=zcl(id7)*(cauxb-ceps_2)
      r8_526(i7,i2)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p7134,qd=p8,v=cf526(i7,i2)%e,a=r8_526(i7,i2)%a,b=r8_526(i7,i
* 2)%b,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cf526(i7,i2)%ek0*(p7134(2)*p8(3)-p8(2)*p7134(3))+p
     & 7134k0*(cf526(i7,i2)%e(2)*p8(3)-p8(2)*cf526(i7,i2)%e(3))-
     & p8k0*(cf526(i7,i2)%e(2)*p7134(3)-p7134(2)*cf526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf526(i7,i2)%e(3)*p8k0+p8(3)*cf526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf526(i7,i2)%e(0)*p7134(0)-cf526(i7,i2)%e(1)*p7134(1)
     & -cf526(i7,i2)%e(2)*p7134(2)-cf526(i7,i2)%e(3)*p7134(3)
      cvqd=cf526(i7,i2)%e(0)*p8(0)-cf526(i7,i2)%e(1)*p8(1)-cf526
     & (i7,i2)%e(2)*p8(2)-cf526(i7,i2)%e(3)*p8(3)
      cauxa=-cf526(i7,i2)%ek0*quqd+p7134k0*cvqd+p8k0*cvqu
      cauxb=-cf526(i7,i2)%ek0*p8(2)+p8k0*cf526(i7,i2)%e(2)
      r8_526(i7,i2)%a(1)=r8_526(i7,i2)%a(1)+fcr(id7)*(cauxa+ceps
     & _0)
      r8_526(i7,i2)%a(2)=r8_526(i7,i2)%a(2)+fcl(id7)*(cauxa-ceps
     & _0)
      r8_526(i7,i2)%b(1)=r8_526(i7,i2)%b(1)+fcl(id7)*(cauxb-ceps
     & _2)
      r8_526(i7,i2)%b(2)=r8_526(i7,i2)%b(2)+fcr(id7)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p7234(m)=p734(m) +p2(m)
       enddo
* pk0 -- p=p7234
      p7234k0=p7234(0)-p7234(1)
* p.q -- p.q=s7234,p=p7234,q=p7234,bef=,aft=
      s7234=(p7234(0)*p7234(0)-p7234(1)*p7234(1)-p7234(2)*p7234(
     & 2)-p7234(3)*p7234(3))
      f7234=s7234*p7234k0
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p7,q=p856
      quqd=p7(0)*p856(0)-p7(1)*p856(1)-p7(2)*p856(2)-p7(3)*p856(
     & 3)
      ccr=zcr(id7)/(f856)
      ccl=zcl(id7)/(f856)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p7,qd=p856,v=cz3214(i3,i1,i2)%e,a=l7_3214(i3,i1,i2)%a,c=l7_3
* 214(i3,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz3214(i3,i1,i2)%ek0*(p7(2)*p856(3)-p856(2)*p7(3))
     & +p7k0*(cz3214(i3,i1,i2)%e(2)*p856(3)-p856(2)*cz3214(i3,i1
     & ,i2)%e(3))-p856k0*(cz3214(i3,i1,i2)%e(2)*p7(3)-p7(2)*cz32
     & 14(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz3214(i3,i1,i2)%e(3)*p7k0+p7(3)*cz3214(i3,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz3214(i3,i1,i2)%e(0)*p7(0)-cz3214(i3,i1,i2)%e(1)*p7(
     & 1)-cz3214(i3,i1,i2)%e(2)*p7(2)-cz3214(i3,i1,i2)%e(3)*p7(3
     & )
      cvqd=cz3214(i3,i1,i2)%e(0)*p856(0)-cz3214(i3,i1,i2)%e(1)*p
     & 856(1)-cz3214(i3,i1,i2)%e(2)*p856(2)-cz3214(i3,i1,i2)%e(3
     & )*p856(3)
      cauxa=-cz3214(i3,i1,i2)%ek0*quqd+p7k0*cvqd+p856k0*cvqu
      cauxc=+cz3214(i3,i1,i2)%ek0*p7(2)-p7k0*cz3214(i3,i1,i2)%e(
     & 2)
      l7_3214(i3,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l7_3214(i3,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_3214(i3,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l7_3214(i3,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id7)/(f856)
      ccl=fcl(id7)/(f856)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p7,qd=p856,v=cf3214(i3,i1,i2)%e,a=l7_3214(i3,i1,i2)%a,c=l7_3
* 214(i3,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf3214(i3,i1,i2)%ek0*(p7(2)*p856(3)-p856(2)*p7(3))
     & +p7k0*(cf3214(i3,i1,i2)%e(2)*p856(3)-p856(2)*cf3214(i3,i1
     & ,i2)%e(3))-p856k0*(cf3214(i3,i1,i2)%e(2)*p7(3)-p7(2)*cf32
     & 14(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf3214(i3,i1,i2)%e(3)*p7k0+p7(3)*cf3214(i3,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf3214(i3,i1,i2)%e(0)*p7(0)-cf3214(i3,i1,i2)%e(1)*p7(
     & 1)-cf3214(i3,i1,i2)%e(2)*p7(2)-cf3214(i3,i1,i2)%e(3)*p7(3
     & )
      cvqd=cf3214(i3,i1,i2)%e(0)*p856(0)-cf3214(i3,i1,i2)%e(1)*p
     & 856(1)-cf3214(i3,i1,i2)%e(2)*p856(2)-cf3214(i3,i1,i2)%e(3
     & )*p856(3)
      cauxa=-cf3214(i3,i1,i2)%ek0*quqd+p7k0*cvqd+p856k0*cvqu
      cauxc=+cf3214(i3,i1,i2)%ek0*p7(2)-p7k0*cf3214(i3,i1,i2)%e(
     & 2)
      l7_3214(i3,i1,i2)%a(1)=l7_3214(i3,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l7_3214(i3,i1,i2)%a(2)=l7_3214(i3,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l7_3214(i3,i1,i2)%c(1)=l7_3214(i3,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l7_3214(i3,i1,i2)%c(2)=l7_3214(i3,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id7).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p7,q=p7234
      quqd=p7(0)*p7234(0)-p7(1)*p7234(1)-p7(2)*p7234(2)-p7(3)*p7
     & 234(3)
      ccr=zcr(id7)/(f7234)
      ccl=zcl(id7)/(f7234)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p7,qd=p7234,v=cz324(i5,i1)%e,a=l7_324(i5,i1)%a,c=l7_324(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz324(i5,i1)%ek0*(p7(2)*p7234(3)-p7234(2)*p7(3))+p
     & 7k0*(cz324(i5,i1)%e(2)*p7234(3)-p7234(2)*cz324(i5,i1)%e(3
     & ))-p7234k0*(cz324(i5,i1)%e(2)*p7(3)-p7(2)*cz324(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i5,i1)%e(3)*p7k0+p7(3)*cz324(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz324(i5,i1)%e(0)*p7(0)-cz324(i5,i1)%e(1)*p7(1)-cz324
     & (i5,i1)%e(2)*p7(2)-cz324(i5,i1)%e(3)*p7(3)
      cvqd=cz324(i5,i1)%e(0)*p7234(0)-cz324(i5,i1)%e(1)*p7234(1)
     & -cz324(i5,i1)%e(2)*p7234(2)-cz324(i5,i1)%e(3)*p7234(3)
      cauxa=-cz324(i5,i1)%ek0*quqd+p7k0*cvqd+p7234k0*cvqu
      cauxc=+cz324(i5,i1)%ek0*p7(2)-p7k0*cz324(i5,i1)%e(2)
      l7_324(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l7_324(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_324(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l7_324(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id7)/(f7234)
      ccl=fcl(id7)/(f7234)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p7,qd=p7234,v=cf324(i5,i1)%e,a=l7_324(i5,i1)%a,c=l7_324(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf324(i5,i1)%ek0*(p7(2)*p7234(3)-p7234(2)*p7(3))+p
     & 7k0*(cf324(i5,i1)%e(2)*p7234(3)-p7234(2)*cf324(i5,i1)%e(3
     & ))-p7234k0*(cf324(i5,i1)%e(2)*p7(3)-p7(2)*cf324(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i5,i1)%e(3)*p7k0+p7(3)*cf324(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf324(i5,i1)%e(0)*p7(0)-cf324(i5,i1)%e(1)*p7(1)-cf324
     & (i5,i1)%e(2)*p7(2)-cf324(i5,i1)%e(3)*p7(3)
      cvqd=cf324(i5,i1)%e(0)*p7234(0)-cf324(i5,i1)%e(1)*p7234(1)
     & -cf324(i5,i1)%e(2)*p7234(2)-cf324(i5,i1)%e(3)*p7234(3)
      cauxa=-cf324(i5,i1)%ek0*quqd+p7k0*cvqd+p7234k0*cvqu
      cauxc=+cf324(i5,i1)%ek0*p7(2)-p7k0*cf324(i5,i1)%e(2)
      l7_324(i5,i1)%a(1)=l7_324(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l7_324(i5,i1)%a(2)=l7_324(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l7_324(i5,i1)%c(1)=l7_324(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l7_324(i5,i1)%c(2)=l7_324(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p734,q=p8
      quqd=p734(0)*p8(0)-p734(1)*p8(1)-p734(2)*p8(2)-p734(3)*p8(
     & 3)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p734,qd=p8,v=cz5216(i5,i1,i2)%e,a=r8_5216(i5,i1,i2)%a,b=r8_5
* 216(i5,i1,i2)%b,cr=zcr(id8),cl=zcl(id8),nsum=0
      ceps_0=-cz5216(i5,i1,i2)%ek0*(p734(2)*p8(3)-p8(2)*p734(3))
     & +p734k0*(cz5216(i5,i1,i2)%e(2)*p8(3)-p8(2)*cz5216(i5,i1,i
     & 2)%e(3))-p8k0*(cz5216(i5,i1,i2)%e(2)*p734(3)-p734(2)*cz52
     & 16(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz5216(i5,i1,i2)%e(3)*p8k0+p8(3)*cz5216(i5,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz5216(i5,i1,i2)%e(0)*p734(0)-cz5216(i5,i1,i2)%e(1)*p
     & 734(1)-cz5216(i5,i1,i2)%e(2)*p734(2)-cz5216(i5,i1,i2)%e(3
     & )*p734(3)
      cvqd=cz5216(i5,i1,i2)%e(0)*p8(0)-cz5216(i5,i1,i2)%e(1)*p8(
     & 1)-cz5216(i5,i1,i2)%e(2)*p8(2)-cz5216(i5,i1,i2)%e(3)*p8(3
     & )
      cauxa=-cz5216(i5,i1,i2)%ek0*quqd+p734k0*cvqd+p8k0*cvqu
      cauxb=-cz5216(i5,i1,i2)%ek0*p8(2)+p8k0*cz5216(i5,i1,i2)%e(
     & 2)
      r8_5216(i5,i1,i2)%a(1)=zcr(id8)*(cauxa+ceps_0)
      r8_5216(i5,i1,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_5216(i5,i1,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      r8_5216(i5,i1,i2)%b(2)=zcr(id8)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p734,qd=p8,v=cf5216(i5,i1,i2)%e,a=r8_5216(i5,i1,i2)%a,b=r8_5
* 216(i5,i1,i2)%b,cr=fcr(id8),cl=fcl(id8),nsum=1
      ceps_0=-cf5216(i5,i1,i2)%ek0*(p734(2)*p8(3)-p8(2)*p734(3))
     & +p734k0*(cf5216(i5,i1,i2)%e(2)*p8(3)-p8(2)*cf5216(i5,i1,i
     & 2)%e(3))-p8k0*(cf5216(i5,i1,i2)%e(2)*p734(3)-p734(2)*cf52
     & 16(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf5216(i5,i1,i2)%e(3)*p8k0+p8(3)*cf5216(i5,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf5216(i5,i1,i2)%e(0)*p734(0)-cf5216(i5,i1,i2)%e(1)*p
     & 734(1)-cf5216(i5,i1,i2)%e(2)*p734(2)-cf5216(i5,i1,i2)%e(3
     & )*p734(3)
      cvqd=cf5216(i5,i1,i2)%e(0)*p8(0)-cf5216(i5,i1,i2)%e(1)*p8(
     & 1)-cf5216(i5,i1,i2)%e(2)*p8(2)-cf5216(i5,i1,i2)%e(3)*p8(3
     & )
      cauxa=-cf5216(i5,i1,i2)%ek0*quqd+p734k0*cvqd+p8k0*cvqu
      cauxb=-cf5216(i5,i1,i2)%ek0*p8(2)+p8k0*cf5216(i5,i1,i2)%e(
     & 2)
      r8_5216(i5,i1,i2)%a(1)=r8_5216(i5,i1,i2)%a(1)+fcr(id8)*(ca
     & uxa+ceps_0)
      r8_5216(i5,i1,i2)%a(2)=r8_5216(i5,i1,i2)%a(2)+fcl(id8)*(ca
     & uxa-ceps_0)
      r8_5216(i5,i1,i2)%b(1)=r8_5216(i5,i1,i2)%b(1)+fcl(id8)*(ca
     & uxb-ceps_2)
      r8_5216(i5,i1,i2)%b(2)=r8_5216(i5,i1,i2)%b(2)+fcr(id8)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id8).ne.1.or.ilept(id3).ne.1) then
  
* quqd -- p=p7234,q=p8
      quqd=p7234(0)*p8(0)-p7234(1)*p8(1)-p7234(2)*p8(2)-p7234(3)
     & *p8(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p7234,qd=p8,v=cz516(i7,i2)%e,a=r8_516(i7,i2)%a,b=r8_516(i7,i
* 2)%b,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cz516(i7,i2)%ek0*(p7234(2)*p8(3)-p8(2)*p7234(3))+p
     & 7234k0*(cz516(i7,i2)%e(2)*p8(3)-p8(2)*cz516(i7,i2)%e(3))-
     & p8k0*(cz516(i7,i2)%e(2)*p7234(3)-p7234(2)*cz516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz516(i7,i2)%e(3)*p8k0+p8(3)*cz516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz516(i7,i2)%e(0)*p7234(0)-cz516(i7,i2)%e(1)*p7234(1)
     & -cz516(i7,i2)%e(2)*p7234(2)-cz516(i7,i2)%e(3)*p7234(3)
      cvqd=cz516(i7,i2)%e(0)*p8(0)-cz516(i7,i2)%e(1)*p8(1)-cz516
     & (i7,i2)%e(2)*p8(2)-cz516(i7,i2)%e(3)*p8(3)
      cauxa=-cz516(i7,i2)%ek0*quqd+p7234k0*cvqd+p8k0*cvqu
      cauxb=-cz516(i7,i2)%ek0*p8(2)+p8k0*cz516(i7,i2)%e(2)
      r8_516(i7,i2)%a(1)=zcr(id7)*(cauxa+ceps_0)
      r8_516(i7,i2)%a(2)=zcl(id7)*(cauxa-ceps_0)
      r8_516(i7,i2)%b(1)=zcl(id7)*(cauxb-ceps_2)
      r8_516(i7,i2)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p7234,qd=p8,v=cf516(i7,i2)%e,a=r8_516(i7,i2)%a,b=r8_516(i7,i
* 2)%b,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cf516(i7,i2)%ek0*(p7234(2)*p8(3)-p8(2)*p7234(3))+p
     & 7234k0*(cf516(i7,i2)%e(2)*p8(3)-p8(2)*cf516(i7,i2)%e(3))-
     & p8k0*(cf516(i7,i2)%e(2)*p7234(3)-p7234(2)*cf516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf516(i7,i2)%e(3)*p8k0+p8(3)*cf516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf516(i7,i2)%e(0)*p7234(0)-cf516(i7,i2)%e(1)*p7234(1)
     & -cf516(i7,i2)%e(2)*p7234(2)-cf516(i7,i2)%e(3)*p7234(3)
      cvqd=cf516(i7,i2)%e(0)*p8(0)-cf516(i7,i2)%e(1)*p8(1)-cf516
     & (i7,i2)%e(2)*p8(2)-cf516(i7,i2)%e(3)*p8(3)
      cauxa=-cf516(i7,i2)%ek0*quqd+p7234k0*cvqd+p8k0*cvqu
      cauxb=-cf516(i7,i2)%ek0*p8(2)+p8k0*cf516(i7,i2)%e(2)
      r8_516(i7,i2)%a(1)=r8_516(i7,i2)%a(1)+fcr(id7)*(cauxa+ceps
     & _0)
      r8_516(i7,i2)%a(2)=r8_516(i7,i2)%a(2)+fcl(id7)*(cauxa-ceps
     & _0)
      r8_516(i7,i2)%b(1)=r8_516(i7,i2)%b(1)+fcl(id7)*(cauxb-ceps
     & _2)
      r8_516(i7,i2)%b(2)=r8_516(i7,i2)%b(2)+fcr(id7)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p756(m)=p7(m)+p5(m)+p6(m)
       enddo
* pk0 -- p=p756
      p756k0=p756(0)-p756(1)
* p.q -- p.q=s756,p=p756,q=p756,bef=,aft=
      s756=(p756(0)*p756(0)-p756(1)*p756(1)-p756(2)*p756(2)-p756
     & (3)*p756(3))
      f756=s756*p756k0
  
       do m=0,3
        p834(m)= -p8(m)-p3(m)-p4(m)
       enddo
* pk0 -- p=p834
      p834k0=p834(0)-p834(1)
* p.q -- p.q=s834,p=p834,q=p834,bef=,aft=
      s834=(p834(0)*p834(0)-p834(1)*p834(1)-p834(2)*p834(2)-p834
     & (3)*p834(3))
      f834=s834*p834k0
  
      if (ilept(id7).ne.1.or.ilept(id3).ne.1) then
      quqd=(s756-s56)/2d0
      ccr=zcr(id7)/(f756)
      ccl=zcl(id7)/(f756)
      do i5=1,2
* TL0 -- qu=p7,qd=p756,v=cz56(i5)%e,a=l7_56(i5)%a,c=l7_56(i5)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p7(2)*p756(3)-p756(2)*p7(3))+p7k0*(c
     & z56(i5)%e(2)*p756(3)-p756(2)*cz56(i5)%e(3))-p756k0*(cz56(
     & i5)%e(2)*p7(3)-p7(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p7k0+p7(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz56(i5)%e(0)*p7(0)-cz56(i5)%e(1)*p7(1)-cz56(i5)%e(2)
     & *p7(2)-cz56(i5)%e(3)*p7(3)
      cvqd=cz56(i5)%e(0)*p756(0)-cz56(i5)%e(1)*p756(1)-cz56(i5)%
     & e(2)*p756(2)-cz56(i5)%e(3)*p756(3)
      cauxa=-cz56(i5)%ek0*quqd+p7k0*cvqd+p756k0*cvqu
      cauxc=+cz56(i5)%ek0*p7(2)-p7k0*cz56(i5)%e(2)
      l7_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      l7_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      l7_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      l7_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      ccr=fcr(id7)/(f756)
      ccl=fcl(id7)/(f756)
      do i5=1,2
* TL0 -- qu=p7,qd=p756,v=cf56(i5)%e,a=l7_56(i5)%a,c=l7_56(i5)%c,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p7(2)*p756(3)-p756(2)*p7(3))+p7k0*(c
     & f56(i5)%e(2)*p756(3)-p756(2)*cf56(i5)%e(3))-p756k0*(cf56(
     & i5)%e(2)*p7(3)-p7(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p7k0+p7(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf56(i5)%e(0)*p7(0)-cf56(i5)%e(1)*p7(1)-cf56(i5)%e(2)
     & *p7(2)-cf56(i5)%e(3)*p7(3)
      cvqd=cf56(i5)%e(0)*p756(0)-cf56(i5)%e(1)*p756(1)-cf56(i5)%
     & e(2)*p756(2)-cf56(i5)%e(3)*p756(3)
      cauxa=-cf56(i5)%ek0*quqd+p7k0*cvqd+p756k0*cvqu
      cauxc=+cf56(i5)%ek0*p7(2)-p7k0*cf56(i5)%e(2)
      l7_56(i5)%a(1)=l7_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      l7_56(i5)%a(2)=l7_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      l7_56(i5)%c(1)=l7_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      l7_56(i5)%c(2)=l7_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      endif
  
      if (ilept(id8).ne.1.or.ilept(id5).ne.1) then
      quqd=(-s834+s34)/2d0
      do i7=1,2
* TR0 -- qu=p834,qd=p8,v=cz34(i7)%e,a=r8_34(i7)%a,b=r8_34(i7)%b,cr=zcr(i
* d8),cl=zcl(id8),nsum=0
      ceps_0=-cz34(i7)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834k0*
     & (cz34(i7)%e(2)*p8(3)-p8(2)*cz34(i7)%e(3))-p8k0*(cz34(i7)%
     & e(2)*p834(3)-p834(2)*cz34(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i7)%e(3)*p8k0+p8(3)*cz34(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i7)%e(0)*p834(0)-cz34(i7)%e(1)*p834(1)-cz34(i7)%
     & e(2)*p834(2)-cz34(i7)%e(3)*p834(3)
      cvqd=cz34(i7)%e(0)*p8(0)-cz34(i7)%e(1)*p8(1)-cz34(i7)%e(2)
     & *p8(2)-cz34(i7)%e(3)*p8(3)
      cauxa=-cz34(i7)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cz34(i7)%ek0*p8(2)+p8k0*cz34(i7)%e(2)
      r8_34(i7)%a(1)=zcr(id8)*(cauxa+ceps_0)
      r8_34(i7)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_34(i7)%b(1)=zcl(id8)*(cauxb-ceps_2)
      r8_34(i7)%b(2)=zcr(id8)*(-cauxb-ceps_2)
      end do
  
      do i7=1,2
* TR0 -- qu=p834,qd=p8,v=cf34(i7)%e,a=r8_34(i7)%a,b=r8_34(i7)%b,cr=fcr(i
* d8),cl=fcl(id8),nsum=1
      ceps_0=-cf34(i7)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834k0*
     & (cf34(i7)%e(2)*p8(3)-p8(2)*cf34(i7)%e(3))-p8k0*(cf34(i7)%
     & e(2)*p834(3)-p834(2)*cf34(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf34(i7)%e(3)*p8k0+p8(3)*cf34(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i7)%e(0)*p834(0)-cf34(i7)%e(1)*p834(1)-cf34(i7)%
     & e(2)*p834(2)-cf34(i7)%e(3)*p834(3)
      cvqd=cf34(i7)%e(0)*p8(0)-cf34(i7)%e(1)*p8(1)-cf34(i7)%e(2)
     & *p8(2)-cf34(i7)%e(3)*p8(3)
      cauxa=-cf34(i7)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cf34(i7)%ek0*p8(2)+p8k0*cf34(i7)%e(2)
      r8_34(i7)%a(1)=r8_34(i7)%a(1)+fcr(id8)*(cauxa+ceps_0)
      r8_34(i7)%a(2)=r8_34(i7)%a(2)+fcl(id8)*(cauxa-ceps_0)
      r8_34(i7)%b(1)=r8_34(i7)%b(1)+fcl(id8)*(cauxb-ceps_2)
      r8_34(i7)%b(2)=r8_34(i7)%b(2)+fcr(id8)*(-cauxb-ceps_2)
      end do
      endif
  
  
       do m=0,3
        p7156(m)=p756(m) +p1(m)
       enddo
* pk0 -- p=p7156
      p7156k0=p7156(0)-p7156(1)
* p.q -- p.q=s7156,p=p7156,q=p7156,bef=,aft=
      s7156=(p7156(0)*p7156(0)-p7156(1)*p7156(1)-p7156(2)*p7156(
     & 2)-p7156(3)*p7156(3))
      f7156=s7156*p7156k0
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p7,q=p834
      quqd=p7(0)*p834(0)-p7(1)*p834(1)-p7(2)*p834(2)-p7(3)*p834(
     & 3)
      ccr=zcr(id7)/(f834)
      ccl=zcl(id7)/(f834)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p7,qd=p834,v=cz5126(i5,i1,i2)%e,a=l7_5126(i5,i1,i2)%a,c=l7_5
* 126(i5,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz5126(i5,i1,i2)%ek0*(p7(2)*p834(3)-p834(2)*p7(3))
     & +p7k0*(cz5126(i5,i1,i2)%e(2)*p834(3)-p834(2)*cz5126(i5,i1
     & ,i2)%e(3))-p834k0*(cz5126(i5,i1,i2)%e(2)*p7(3)-p7(2)*cz51
     & 26(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz5126(i5,i1,i2)%e(3)*p7k0+p7(3)*cz5126(i5,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz5126(i5,i1,i2)%e(0)*p7(0)-cz5126(i5,i1,i2)%e(1)*p7(
     & 1)-cz5126(i5,i1,i2)%e(2)*p7(2)-cz5126(i5,i1,i2)%e(3)*p7(3
     & )
      cvqd=cz5126(i5,i1,i2)%e(0)*p834(0)-cz5126(i5,i1,i2)%e(1)*p
     & 834(1)-cz5126(i5,i1,i2)%e(2)*p834(2)-cz5126(i5,i1,i2)%e(3
     & )*p834(3)
      cauxa=-cz5126(i5,i1,i2)%ek0*quqd+p7k0*cvqd+p834k0*cvqu
      cauxc=+cz5126(i5,i1,i2)%ek0*p7(2)-p7k0*cz5126(i5,i1,i2)%e(
     & 2)
      l7_5126(i5,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l7_5126(i5,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_5126(i5,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l7_5126(i5,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id7)/(f834)
      ccl=fcl(id7)/(f834)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p7,qd=p834,v=cf5126(i5,i1,i2)%e,a=l7_5126(i5,i1,i2)%a,c=l7_5
* 126(i5,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf5126(i5,i1,i2)%ek0*(p7(2)*p834(3)-p834(2)*p7(3))
     & +p7k0*(cf5126(i5,i1,i2)%e(2)*p834(3)-p834(2)*cf5126(i5,i1
     & ,i2)%e(3))-p834k0*(cf5126(i5,i1,i2)%e(2)*p7(3)-p7(2)*cf51
     & 26(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf5126(i5,i1,i2)%e(3)*p7k0+p7(3)*cf5126(i5,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf5126(i5,i1,i2)%e(0)*p7(0)-cf5126(i5,i1,i2)%e(1)*p7(
     & 1)-cf5126(i5,i1,i2)%e(2)*p7(2)-cf5126(i5,i1,i2)%e(3)*p7(3
     & )
      cvqd=cf5126(i5,i1,i2)%e(0)*p834(0)-cf5126(i5,i1,i2)%e(1)*p
     & 834(1)-cf5126(i5,i1,i2)%e(2)*p834(2)-cf5126(i5,i1,i2)%e(3
     & )*p834(3)
      cauxa=-cf5126(i5,i1,i2)%ek0*quqd+p7k0*cvqd+p834k0*cvqu
      cauxc=+cf5126(i5,i1,i2)%ek0*p7(2)-p7k0*cf5126(i5,i1,i2)%e(
     & 2)
      l7_5126(i5,i1,i2)%a(1)=l7_5126(i5,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l7_5126(i5,i1,i2)%a(2)=l7_5126(i5,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l7_5126(i5,i1,i2)%c(1)=l7_5126(i5,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l7_5126(i5,i1,i2)%c(2)=l7_5126(i5,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id7).ne.1.or.ilept(id3).ne.1) then
  
* quqd -- p=p7,q=p7156
      quqd=p7(0)*p7156(0)-p7(1)*p7156(1)-p7(2)*p7156(2)-p7(3)*p7
     & 156(3)
      ccr=zcr(id7)/(f7156)
      ccl=zcl(id7)/(f7156)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p7,qd=p7156,v=cz516(i5,i1)%e,a=l7_516(i5,i1)%a,c=l7_516(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz516(i5,i1)%ek0*(p7(2)*p7156(3)-p7156(2)*p7(3))+p
     & 7k0*(cz516(i5,i1)%e(2)*p7156(3)-p7156(2)*cz516(i5,i1)%e(3
     & ))-p7156k0*(cz516(i5,i1)%e(2)*p7(3)-p7(2)*cz516(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz516(i5,i1)%e(3)*p7k0+p7(3)*cz516(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz516(i5,i1)%e(0)*p7(0)-cz516(i5,i1)%e(1)*p7(1)-cz516
     & (i5,i1)%e(2)*p7(2)-cz516(i5,i1)%e(3)*p7(3)
      cvqd=cz516(i5,i1)%e(0)*p7156(0)-cz516(i5,i1)%e(1)*p7156(1)
     & -cz516(i5,i1)%e(2)*p7156(2)-cz516(i5,i1)%e(3)*p7156(3)
      cauxa=-cz516(i5,i1)%ek0*quqd+p7k0*cvqd+p7156k0*cvqu
      cauxc=+cz516(i5,i1)%ek0*p7(2)-p7k0*cz516(i5,i1)%e(2)
      l7_516(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l7_516(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_516(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l7_516(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id7)/(f7156)
      ccl=fcl(id7)/(f7156)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p7,qd=p7156,v=cf516(i5,i1)%e,a=l7_516(i5,i1)%a,c=l7_516(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf516(i5,i1)%ek0*(p7(2)*p7156(3)-p7156(2)*p7(3))+p
     & 7k0*(cf516(i5,i1)%e(2)*p7156(3)-p7156(2)*cf516(i5,i1)%e(3
     & ))-p7156k0*(cf516(i5,i1)%e(2)*p7(3)-p7(2)*cf516(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf516(i5,i1)%e(3)*p7k0+p7(3)*cf516(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf516(i5,i1)%e(0)*p7(0)-cf516(i5,i1)%e(1)*p7(1)-cf516
     & (i5,i1)%e(2)*p7(2)-cf516(i5,i1)%e(3)*p7(3)
      cvqd=cf516(i5,i1)%e(0)*p7156(0)-cf516(i5,i1)%e(1)*p7156(1)
     & -cf516(i5,i1)%e(2)*p7156(2)-cf516(i5,i1)%e(3)*p7156(3)
      cauxa=-cf516(i5,i1)%ek0*quqd+p7k0*cvqd+p7156k0*cvqu
      cauxc=+cf516(i5,i1)%ek0*p7(2)-p7k0*cf516(i5,i1)%e(2)
      l7_516(i5,i1)%a(1)=l7_516(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l7_516(i5,i1)%a(2)=l7_516(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l7_516(i5,i1)%c(1)=l7_516(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l7_516(i5,i1)%c(2)=l7_516(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p756,q=p8
      quqd=p756(0)*p8(0)-p756(1)*p8(1)-p756(2)*p8(2)-p756(3)*p8(
     & 3)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p756,qd=p8,v=cz3124(i3,i1,i2)%e,a=r8_3124(i3,i1,i2)%a,b=r8_3
* 124(i3,i1,i2)%b,cr=zcr(id8),cl=zcl(id8),nsum=0
      ceps_0=-cz3124(i3,i1,i2)%ek0*(p756(2)*p8(3)-p8(2)*p756(3))
     & +p756k0*(cz3124(i3,i1,i2)%e(2)*p8(3)-p8(2)*cz3124(i3,i1,i
     & 2)%e(3))-p8k0*(cz3124(i3,i1,i2)%e(2)*p756(3)-p756(2)*cz31
     & 24(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz3124(i3,i1,i2)%e(3)*p8k0+p8(3)*cz3124(i3,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz3124(i3,i1,i2)%e(0)*p756(0)-cz3124(i3,i1,i2)%e(1)*p
     & 756(1)-cz3124(i3,i1,i2)%e(2)*p756(2)-cz3124(i3,i1,i2)%e(3
     & )*p756(3)
      cvqd=cz3124(i3,i1,i2)%e(0)*p8(0)-cz3124(i3,i1,i2)%e(1)*p8(
     & 1)-cz3124(i3,i1,i2)%e(2)*p8(2)-cz3124(i3,i1,i2)%e(3)*p8(3
     & )
      cauxa=-cz3124(i3,i1,i2)%ek0*quqd+p756k0*cvqd+p8k0*cvqu
      cauxb=-cz3124(i3,i1,i2)%ek0*p8(2)+p8k0*cz3124(i3,i1,i2)%e(
     & 2)
      r8_3124(i3,i1,i2)%a(1)=zcr(id8)*(cauxa+ceps_0)
      r8_3124(i3,i1,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_3124(i3,i1,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      r8_3124(i3,i1,i2)%b(2)=zcr(id8)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p756,qd=p8,v=cf3124(i3,i1,i2)%e,a=r8_3124(i3,i1,i2)%a,b=r8_3
* 124(i3,i1,i2)%b,cr=fcr(id8),cl=fcl(id8),nsum=1
      ceps_0=-cf3124(i3,i1,i2)%ek0*(p756(2)*p8(3)-p8(2)*p756(3))
     & +p756k0*(cf3124(i3,i1,i2)%e(2)*p8(3)-p8(2)*cf3124(i3,i1,i
     & 2)%e(3))-p8k0*(cf3124(i3,i1,i2)%e(2)*p756(3)-p756(2)*cf31
     & 24(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf3124(i3,i1,i2)%e(3)*p8k0+p8(3)*cf3124(i3,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf3124(i3,i1,i2)%e(0)*p756(0)-cf3124(i3,i1,i2)%e(1)*p
     & 756(1)-cf3124(i3,i1,i2)%e(2)*p756(2)-cf3124(i3,i1,i2)%e(3
     & )*p756(3)
      cvqd=cf3124(i3,i1,i2)%e(0)*p8(0)-cf3124(i3,i1,i2)%e(1)*p8(
     & 1)-cf3124(i3,i1,i2)%e(2)*p8(2)-cf3124(i3,i1,i2)%e(3)*p8(3
     & )
      cauxa=-cf3124(i3,i1,i2)%ek0*quqd+p756k0*cvqd+p8k0*cvqu
      cauxb=-cf3124(i3,i1,i2)%ek0*p8(2)+p8k0*cf3124(i3,i1,i2)%e(
     & 2)
      r8_3124(i3,i1,i2)%a(1)=r8_3124(i3,i1,i2)%a(1)+fcr(id8)*(ca
     & uxa+ceps_0)
      r8_3124(i3,i1,i2)%a(2)=r8_3124(i3,i1,i2)%a(2)+fcl(id8)*(ca
     & uxa-ceps_0)
      r8_3124(i3,i1,i2)%b(1)=r8_3124(i3,i1,i2)%b(1)+fcl(id8)*(ca
     & uxb-ceps_2)
      r8_3124(i3,i1,i2)%b(2)=r8_3124(i3,i1,i2)%b(2)+fcr(id8)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id8).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p7156,q=p8
      quqd=p7156(0)*p8(0)-p7156(1)*p8(1)-p7156(2)*p8(2)-p7156(3)
     & *p8(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p7156,qd=p8,v=cz324(i7,i2)%e,a=r8_324(i7,i2)%a,b=r8_324(i7,i
* 2)%b,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cz324(i7,i2)%ek0*(p7156(2)*p8(3)-p8(2)*p7156(3))+p
     & 7156k0*(cz324(i7,i2)%e(2)*p8(3)-p8(2)*cz324(i7,i2)%e(3))-
     & p8k0*(cz324(i7,i2)%e(2)*p7156(3)-p7156(2)*cz324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz324(i7,i2)%e(3)*p8k0+p8(3)*cz324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i7,i2)%e(0)*p7156(0)-cz324(i7,i2)%e(1)*p7156(1)
     & -cz324(i7,i2)%e(2)*p7156(2)-cz324(i7,i2)%e(3)*p7156(3)
      cvqd=cz324(i7,i2)%e(0)*p8(0)-cz324(i7,i2)%e(1)*p8(1)-cz324
     & (i7,i2)%e(2)*p8(2)-cz324(i7,i2)%e(3)*p8(3)
      cauxa=-cz324(i7,i2)%ek0*quqd+p7156k0*cvqd+p8k0*cvqu
      cauxb=-cz324(i7,i2)%ek0*p8(2)+p8k0*cz324(i7,i2)%e(2)
      r8_324(i7,i2)%a(1)=zcr(id7)*(cauxa+ceps_0)
      r8_324(i7,i2)%a(2)=zcl(id7)*(cauxa-ceps_0)
      r8_324(i7,i2)%b(1)=zcl(id7)*(cauxb-ceps_2)
      r8_324(i7,i2)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p7156,qd=p8,v=cf324(i7,i2)%e,a=r8_324(i7,i2)%a,b=r8_324(i7,i
* 2)%b,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cf324(i7,i2)%ek0*(p7156(2)*p8(3)-p8(2)*p7156(3))+p
     & 7156k0*(cf324(i7,i2)%e(2)*p8(3)-p8(2)*cf324(i7,i2)%e(3))-
     & p8k0*(cf324(i7,i2)%e(2)*p7156(3)-p7156(2)*cf324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf324(i7,i2)%e(3)*p8k0+p8(3)*cf324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i7,i2)%e(0)*p7156(0)-cf324(i7,i2)%e(1)*p7156(1)
     & -cf324(i7,i2)%e(2)*p7156(2)-cf324(i7,i2)%e(3)*p7156(3)
      cvqd=cf324(i7,i2)%e(0)*p8(0)-cf324(i7,i2)%e(1)*p8(1)-cf324
     & (i7,i2)%e(2)*p8(2)-cf324(i7,i2)%e(3)*p8(3)
      cauxa=-cf324(i7,i2)%ek0*quqd+p7156k0*cvqd+p8k0*cvqu
      cauxb=-cf324(i7,i2)%ek0*p8(2)+p8k0*cf324(i7,i2)%e(2)
      r8_324(i7,i2)%a(1)=r8_324(i7,i2)%a(1)+fcr(id7)*(cauxa+ceps
     & _0)
      r8_324(i7,i2)%a(2)=r8_324(i7,i2)%a(2)+fcl(id7)*(cauxa-ceps
     & _0)
      r8_324(i7,i2)%b(1)=r8_324(i7,i2)%b(1)+fcl(id7)*(cauxb-ceps
     & _2)
      r8_324(i7,i2)%b(2)=r8_324(i7,i2)%b(2)+fcr(id7)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
       do m=0,3
        p7256(m)=p756(m) +p2(m)
       enddo
* pk0 -- p=p7256
      p7256k0=p7256(0)-p7256(1)
* p.q -- p.q=s7256,p=p7256,q=p7256,bef=,aft=
      s7256=(p7256(0)*p7256(0)-p7256(1)*p7256(1)-p7256(2)*p7256(
     & 2)-p7256(3)*p7256(3))
      f7256=s7256*p7256k0
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p7,q=p834
      quqd=p7(0)*p834(0)-p7(1)*p834(1)-p7(2)*p834(2)-p7(3)*p834(
     & 3)
      ccr=zcr(id7)/(f834)
      ccl=zcl(id7)/(f834)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p7,qd=p834,v=cz5216(i5,i1,i2)%e,a=l7_5216(i5,i1,i2)%a,c=l7_5
* 216(i5,i1,i2)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz5216(i5,i1,i2)%ek0*(p7(2)*p834(3)-p834(2)*p7(3))
     & +p7k0*(cz5216(i5,i1,i2)%e(2)*p834(3)-p834(2)*cz5216(i5,i1
     & ,i2)%e(3))-p834k0*(cz5216(i5,i1,i2)%e(2)*p7(3)-p7(2)*cz52
     & 16(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz5216(i5,i1,i2)%e(3)*p7k0+p7(3)*cz5216(i5,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cz5216(i5,i1,i2)%e(0)*p7(0)-cz5216(i5,i1,i2)%e(1)*p7(
     & 1)-cz5216(i5,i1,i2)%e(2)*p7(2)-cz5216(i5,i1,i2)%e(3)*p7(3
     & )
      cvqd=cz5216(i5,i1,i2)%e(0)*p834(0)-cz5216(i5,i1,i2)%e(1)*p
     & 834(1)-cz5216(i5,i1,i2)%e(2)*p834(2)-cz5216(i5,i1,i2)%e(3
     & )*p834(3)
      cauxa=-cz5216(i5,i1,i2)%ek0*quqd+p7k0*cvqd+p834k0*cvqu
      cauxc=+cz5216(i5,i1,i2)%ek0*p7(2)-p7k0*cz5216(i5,i1,i2)%e(
     & 2)
      l7_5216(i5,i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      l7_5216(i5,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_5216(i5,i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      l7_5216(i5,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
      ccr=fcr(id7)/(f834)
      ccl=fcl(id7)/(f834)
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TL0 -- qu=p7,qd=p834,v=cf5216(i5,i1,i2)%e,a=l7_5216(i5,i1,i2)%a,c=l7_5
* 216(i5,i1,i2)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf5216(i5,i1,i2)%ek0*(p7(2)*p834(3)-p834(2)*p7(3))
     & +p7k0*(cf5216(i5,i1,i2)%e(2)*p834(3)-p834(2)*cf5216(i5,i1
     & ,i2)%e(3))-p834k0*(cf5216(i5,i1,i2)%e(2)*p7(3)-p7(2)*cf52
     & 16(i5,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf5216(i5,i1,i2)%e(3)*p7k0+p7(3)*cf5216(i5,i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      cvqu=cf5216(i5,i1,i2)%e(0)*p7(0)-cf5216(i5,i1,i2)%e(1)*p7(
     & 1)-cf5216(i5,i1,i2)%e(2)*p7(2)-cf5216(i5,i1,i2)%e(3)*p7(3
     & )
      cvqd=cf5216(i5,i1,i2)%e(0)*p834(0)-cf5216(i5,i1,i2)%e(1)*p
     & 834(1)-cf5216(i5,i1,i2)%e(2)*p834(2)-cf5216(i5,i1,i2)%e(3
     & )*p834(3)
      cauxa=-cf5216(i5,i1,i2)%ek0*quqd+p7k0*cvqd+p834k0*cvqu
      cauxc=+cf5216(i5,i1,i2)%ek0*p7(2)-p7k0*cf5216(i5,i1,i2)%e(
     & 2)
      l7_5216(i5,i1,i2)%a(1)=l7_5216(i5,i1,i2)%a(1)+ccr*(cauxa+c
     & eps_0)
      l7_5216(i5,i1,i2)%a(2)=l7_5216(i5,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l7_5216(i5,i1,i2)%c(1)=l7_5216(i5,i1,i2)%c(1)+ccr*(cauxc+c
     & eps_1)
      l7_5216(i5,i1,i2)%c(2)=l7_5216(i5,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
  
      if (ilept(id7).ne.1.or.ilept(id3).ne.1) then
  
* quqd -- p=p7,q=p7256
      quqd=p7(0)*p7256(0)-p7(1)*p7256(1)-p7(2)*p7256(2)-p7(3)*p7
     & 256(3)
      ccr=zcr(id7)/(f7256)
      ccl=zcl(id7)/(f7256)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p7,qd=p7256,v=cz526(i5,i1)%e,a=l7_526(i5,i1)%a,c=l7_526(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz526(i5,i1)%ek0*(p7(2)*p7256(3)-p7256(2)*p7(3))+p
     & 7k0*(cz526(i5,i1)%e(2)*p7256(3)-p7256(2)*cz526(i5,i1)%e(3
     & ))-p7256k0*(cz526(i5,i1)%e(2)*p7(3)-p7(2)*cz526(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz526(i5,i1)%e(3)*p7k0+p7(3)*cz526(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz526(i5,i1)%e(0)*p7(0)-cz526(i5,i1)%e(1)*p7(1)-cz526
     & (i5,i1)%e(2)*p7(2)-cz526(i5,i1)%e(3)*p7(3)
      cvqd=cz526(i5,i1)%e(0)*p7256(0)-cz526(i5,i1)%e(1)*p7256(1)
     & -cz526(i5,i1)%e(2)*p7256(2)-cz526(i5,i1)%e(3)*p7256(3)
      cauxa=-cz526(i5,i1)%ek0*quqd+p7k0*cvqd+p7256k0*cvqu
      cauxc=+cz526(i5,i1)%ek0*p7(2)-p7k0*cz526(i5,i1)%e(2)
      l7_526(i5,i1)%a(1)=ccr*(cauxa+ceps_0)
      l7_526(i5,i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_526(i5,i1)%c(1)=ccr*(cauxc+ceps_1)
      l7_526(i5,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      ccr=fcr(id7)/(f7256)
      ccl=fcl(id7)/(f7256)
      do i5=1,2
      do i1=1,2
* TL0 -- qu=p7,qd=p7256,v=cf526(i5,i1)%e,a=l7_526(i5,i1)%a,c=l7_526(i5,i
* 1)%c,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf526(i5,i1)%ek0*(p7(2)*p7256(3)-p7256(2)*p7(3))+p
     & 7k0*(cf526(i5,i1)%e(2)*p7256(3)-p7256(2)*cf526(i5,i1)%e(3
     & ))-p7256k0*(cf526(i5,i1)%e(2)*p7(3)-p7(2)*cf526(i5,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf526(i5,i1)%e(3)*p7k0+p7(3)*cf526(i5,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf526(i5,i1)%e(0)*p7(0)-cf526(i5,i1)%e(1)*p7(1)-cf526
     & (i5,i1)%e(2)*p7(2)-cf526(i5,i1)%e(3)*p7(3)
      cvqd=cf526(i5,i1)%e(0)*p7256(0)-cf526(i5,i1)%e(1)*p7256(1)
     & -cf526(i5,i1)%e(2)*p7256(2)-cf526(i5,i1)%e(3)*p7256(3)
      cauxa=-cf526(i5,i1)%ek0*quqd+p7k0*cvqd+p7256k0*cvqu
      cauxc=+cf526(i5,i1)%ek0*p7(2)-p7k0*cf526(i5,i1)%e(2)
      l7_526(i5,i1)%a(1)=l7_526(i5,i1)%a(1)+ccr*(cauxa+ceps_0)
      l7_526(i5,i1)%a(2)=l7_526(i5,i1)%a(2)+ccl*(cauxa-ceps_0)
      l7_526(i5,i1)%c(1)=l7_526(i5,i1)%c(1)+ccr*(cauxc+ceps_1)
      l7_526(i5,i1)%c(2)=l7_526(i5,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
      endif
      endif
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p756,q=p8
      quqd=p756(0)*p8(0)-p756(1)*p8(1)-p756(2)*p8(2)-p756(3)*p8(
     & 3)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p756,qd=p8,v=cz3214(i3,i1,i2)%e,a=r8_3214(i3,i1,i2)%a,b=r8_3
* 214(i3,i1,i2)%b,cr=zcr(id8),cl=zcl(id8),nsum=0
      ceps_0=-cz3214(i3,i1,i2)%ek0*(p756(2)*p8(3)-p8(2)*p756(3))
     & +p756k0*(cz3214(i3,i1,i2)%e(2)*p8(3)-p8(2)*cz3214(i3,i1,i
     & 2)%e(3))-p8k0*(cz3214(i3,i1,i2)%e(2)*p756(3)-p756(2)*cz32
     & 14(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz3214(i3,i1,i2)%e(3)*p8k0+p8(3)*cz3214(i3,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cz3214(i3,i1,i2)%e(0)*p756(0)-cz3214(i3,i1,i2)%e(1)*p
     & 756(1)-cz3214(i3,i1,i2)%e(2)*p756(2)-cz3214(i3,i1,i2)%e(3
     & )*p756(3)
      cvqd=cz3214(i3,i1,i2)%e(0)*p8(0)-cz3214(i3,i1,i2)%e(1)*p8(
     & 1)-cz3214(i3,i1,i2)%e(2)*p8(2)-cz3214(i3,i1,i2)%e(3)*p8(3
     & )
      cauxa=-cz3214(i3,i1,i2)%ek0*quqd+p756k0*cvqd+p8k0*cvqu
      cauxb=-cz3214(i3,i1,i2)%ek0*p8(2)+p8k0*cz3214(i3,i1,i2)%e(
     & 2)
      r8_3214(i3,i1,i2)%a(1)=zcr(id8)*(cauxa+ceps_0)
      r8_3214(i3,i1,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_3214(i3,i1,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      r8_3214(i3,i1,i2)%b(2)=zcr(id8)*(-cauxb-ceps_2)
      end do
      end do
      end do
  
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p756,qd=p8,v=cf3214(i3,i1,i2)%e,a=r8_3214(i3,i1,i2)%a,b=r8_3
* 214(i3,i1,i2)%b,cr=fcr(id8),cl=fcl(id8),nsum=1
      ceps_0=-cf3214(i3,i1,i2)%ek0*(p756(2)*p8(3)-p8(2)*p756(3))
     & +p756k0*(cf3214(i3,i1,i2)%e(2)*p8(3)-p8(2)*cf3214(i3,i1,i
     & 2)%e(3))-p8k0*(cf3214(i3,i1,i2)%e(2)*p756(3)-p756(2)*cf32
     & 14(i3,i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf3214(i3,i1,i2)%e(3)*p8k0+p8(3)*cf3214(i3,i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=cf3214(i3,i1,i2)%e(0)*p756(0)-cf3214(i3,i1,i2)%e(1)*p
     & 756(1)-cf3214(i3,i1,i2)%e(2)*p756(2)-cf3214(i3,i1,i2)%e(3
     & )*p756(3)
      cvqd=cf3214(i3,i1,i2)%e(0)*p8(0)-cf3214(i3,i1,i2)%e(1)*p8(
     & 1)-cf3214(i3,i1,i2)%e(2)*p8(2)-cf3214(i3,i1,i2)%e(3)*p8(3
     & )
      cauxa=-cf3214(i3,i1,i2)%ek0*quqd+p756k0*cvqd+p8k0*cvqu
      cauxb=-cf3214(i3,i1,i2)%ek0*p8(2)+p8k0*cf3214(i3,i1,i2)%e(
     & 2)
      r8_3214(i3,i1,i2)%a(1)=r8_3214(i3,i1,i2)%a(1)+fcr(id8)*(ca
     & uxa+ceps_0)
      r8_3214(i3,i1,i2)%a(2)=r8_3214(i3,i1,i2)%a(2)+fcl(id8)*(ca
     & uxa-ceps_0)
      r8_3214(i3,i1,i2)%b(1)=r8_3214(i3,i1,i2)%b(1)+fcl(id8)*(ca
     & uxb-ceps_2)
      r8_3214(i3,i1,i2)%b(2)=r8_3214(i3,i1,i2)%b(2)+fcr(id8)*(-c
     & auxb-ceps_2)
      end do
      end do
      end do
  
      if (ilept(id8).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p7256,q=p8
      quqd=p7256(0)*p8(0)-p7256(1)*p8(1)-p7256(2)*p8(2)-p7256(3)
     & *p8(3)
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p7256,qd=p8,v=cz314(i7,i2)%e,a=r8_314(i7,i2)%a,b=r8_314(i7,i
* 2)%b,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cz314(i7,i2)%ek0*(p7256(2)*p8(3)-p8(2)*p7256(3))+p
     & 7256k0*(cz314(i7,i2)%e(2)*p8(3)-p8(2)*cz314(i7,i2)%e(3))-
     & p8k0*(cz314(i7,i2)%e(2)*p7256(3)-p7256(2)*cz314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz314(i7,i2)%e(3)*p8k0+p8(3)*cz314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i7,i2)%e(0)*p7256(0)-cz314(i7,i2)%e(1)*p7256(1)
     & -cz314(i7,i2)%e(2)*p7256(2)-cz314(i7,i2)%e(3)*p7256(3)
      cvqd=cz314(i7,i2)%e(0)*p8(0)-cz314(i7,i2)%e(1)*p8(1)-cz314
     & (i7,i2)%e(2)*p8(2)-cz314(i7,i2)%e(3)*p8(3)
      cauxa=-cz314(i7,i2)%ek0*quqd+p7256k0*cvqd+p8k0*cvqu
      cauxb=-cz314(i7,i2)%ek0*p8(2)+p8k0*cz314(i7,i2)%e(2)
      r8_314(i7,i2)%a(1)=zcr(id7)*(cauxa+ceps_0)
      r8_314(i7,i2)%a(2)=zcl(id7)*(cauxa-ceps_0)
      r8_314(i7,i2)%b(1)=zcl(id7)*(cauxb-ceps_2)
      r8_314(i7,i2)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      end do
      end do
  
      do i7=1,2
      do i2=1,2
* TR0 -- qu=p7256,qd=p8,v=cf314(i7,i2)%e,a=r8_314(i7,i2)%a,b=r8_314(i7,i
* 2)%b,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cf314(i7,i2)%ek0*(p7256(2)*p8(3)-p8(2)*p7256(3))+p
     & 7256k0*(cf314(i7,i2)%e(2)*p8(3)-p8(2)*cf314(i7,i2)%e(3))-
     & p8k0*(cf314(i7,i2)%e(2)*p7256(3)-p7256(2)*cf314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf314(i7,i2)%e(3)*p8k0+p8(3)*cf314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i7,i2)%e(0)*p7256(0)-cf314(i7,i2)%e(1)*p7256(1)
     & -cf314(i7,i2)%e(2)*p7256(2)-cf314(i7,i2)%e(3)*p7256(3)
      cvqd=cf314(i7,i2)%e(0)*p8(0)-cf314(i7,i2)%e(1)*p8(1)-cf314
     & (i7,i2)%e(2)*p8(2)-cf314(i7,i2)%e(3)*p8(3)
      cauxa=-cf314(i7,i2)%ek0*quqd+p7256k0*cvqd+p8k0*cvqu
      cauxb=-cf314(i7,i2)%ek0*p8(2)+p8k0*cf314(i7,i2)%e(2)
      r8_314(i7,i2)%a(1)=r8_314(i7,i2)%a(1)+fcr(id7)*(cauxa+ceps
     & _0)
      r8_314(i7,i2)%a(2)=r8_314(i7,i2)%a(2)+fcl(id7)*(cauxa-ceps
     & _0)
      r8_314(i7,i2)%b(1)=r8_314(i7,i2)%b(1)+fcl(id7)*(cauxb-ceps
     & _2)
      r8_314(i7,i2)%b(2)=r8_314(i7,i2)%b(2)+fcr(id7)*(-cauxb-cep
     & s_2)
      end do
      end do
  
      endif
      endif
  
* COLOR STRUCTURE                                                       
*                                                                       
* col            col                                                    
* 1   3-1 5-2    7    3-12                                              
* 2   3-2 5-1    8    3-21                                              
* 3   3-1 7-2    9    5-12                                              
* 4   3-2 7-1   10    5-21                                              
* 5   5-1 7-2   11    7-12                                              
* 6   5-2 7-1   12    7-21                                              
  
* line without gluon; col 1-6                                           
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,1),a1=l7_314(i3,i1)%a,c1=l7_314(i3,i1)%
* c,a2=r8_526(i5,i2)%a,b2=r8_526(i5,i2)%b,prq=s7134,bef=,aft=
      cres(i1,i2,i3,i5,1,1)=(l7_314(i3,i1)%a(1)*r8_526(i5,i2)%a(
     & 1)+l7_314(i3,i1)%c(1)*s7134*r8_526(i5,i2)%b(2))
      cres(i1,i2,i3,i5,2,1)=(l7_314(i3,i1)%c(2)*s7134*r8_526(i5,
     & i2)%b(1)+l7_314(i3,i1)%a(2)*r8_526(i5,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,1),a1=l7_526(i5,i2)%a,c1=l7_526(i5,i2)%
* c,a2=r8_314(i3,i1)%a,b2=r8_314(i3,i1)%b,prq=s7256,bef=cres(i1,i2,i3,i5
* ,&,1)+,aft=
      cres(i1,i2,i3,i5,1,1)=cres(i1,i2,i3,i5,1,1)+(l7_526(i5,i2)
     & %a(1)*r8_314(i3,i1)%a(1)+l7_526(i5,i2)%c(1)*s7256*r8_314(
     & i3,i1)%b(2))
      cres(i1,i2,i3,i5,2,1)=cres(i1,i2,i3,i5,2,1)+(l7_526(i5,i2)
     & %c(2)*s7256*r8_314(i3,i1)%b(1)+l7_526(i5,i2)%a(2)*r8_314(
     & i3,i1)%a(2))
      end do
      end do
      end do
      end do
      endif
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,2),a1=l7_324(i3,i2)%a,c1=l7_324(i3,i2)%
* c,a2=r8_516(i5,i1)%a,b2=r8_516(i5,i1)%b,prq=s7234,bef=,aft=
      cres(i1,i2,i3,i5,1,2)=(l7_324(i3,i2)%a(1)*r8_516(i5,i1)%a(
     & 1)+l7_324(i3,i2)%c(1)*s7234*r8_516(i5,i1)%b(2))
      cres(i1,i2,i3,i5,2,2)=(l7_324(i3,i2)%c(2)*s7234*r8_516(i5,
     & i1)%b(1)+l7_324(i3,i2)%a(2)*r8_516(i5,i1)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,2),a1=l7_516(i5,i1)%a,c1=l7_516(i5,i1)%
* c,a2=r8_324(i3,i2)%a,b2=r8_324(i3,i2)%b,prq=s7156,bef=cres(i1,i2,i3,i5
* ,&,2)+,aft=
      cres(i1,i2,i3,i5,1,2)=cres(i1,i2,i3,i5,1,2)+(l7_516(i5,i1)
     & %a(1)*r8_324(i3,i2)%a(1)+l7_516(i5,i1)%c(1)*s7156*r8_324(
     & i3,i2)%b(2))
      cres(i1,i2,i3,i5,2,2)=cres(i1,i2,i3,i5,2,2)+(l7_516(i5,i1)
     & %c(2)*s7156*r8_324(i3,i2)%b(1)+l7_516(i5,i1)%a(2)*r8_324(
     & i3,i2)%a(2))
      end do
      end do
      end do
      end do
      endif
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,3),a1=l5_314(i3,i1)%a,c1=l5_314(i3,i1)%
* c,a2=r6_728(i7,i2)%a,b2=r6_728(i7,i2)%b,prq=s5134,bef=,aft=
      cres(i1,i2,i3,1,i7,3)=(l5_314(i3,i1)%a(1)*r6_728(i7,i2)%a(
     & 1)+l5_314(i3,i1)%c(1)*s5134*r6_728(i7,i2)%b(2))
      cres(i1,i2,i3,2,i7,3)=(l5_314(i3,i1)%c(2)*s5134*r6_728(i7,
     & i2)%b(1)+l5_314(i3,i1)%a(2)*r6_728(i7,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,3),a1=l5_728(i7,i2)%a,c1=l5_728(i7,i2)%
* c,a2=r6_314(i3,i1)%a,b2=r6_314(i3,i1)%b,prq=s5278,bef=cres(i1,i2,i3,&,
* i7,3)+,aft=
      cres(i1,i2,i3,1,i7,3)=cres(i1,i2,i3,1,i7,3)+(l5_728(i7,i2)
     & %a(1)*r6_314(i3,i1)%a(1)+l5_728(i7,i2)%c(1)*s5278*r6_314(
     & i3,i1)%b(2))
      cres(i1,i2,i3,2,i7,3)=cres(i1,i2,i3,2,i7,3)+(l5_728(i7,i2)
     & %c(2)*s5278*r6_314(i3,i1)%b(1)+l5_728(i7,i2)%a(2)*r6_314(
     & i3,i1)%a(2))
      end do
      end do
      end do
      end do
      endif
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,4),a1=l5_324(i3,i2)%a,c1=l5_324(i3,i2)%
* c,a2=r6_718(i7,i1)%a,b2=r6_718(i7,i1)%b,prq=s5234,bef=,aft=
      cres(i1,i2,i3,1,i7,4)=(l5_324(i3,i2)%a(1)*r6_718(i7,i1)%a(
     & 1)+l5_324(i3,i2)%c(1)*s5234*r6_718(i7,i1)%b(2))
      cres(i1,i2,i3,2,i7,4)=(l5_324(i3,i2)%c(2)*s5234*r6_718(i7,
     & i1)%b(1)+l5_324(i3,i2)%a(2)*r6_718(i7,i1)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,4),a1=l5_718(i7,i1)%a,c1=l5_718(i7,i1)%
* c,a2=r6_324(i3,i2)%a,b2=r6_324(i3,i2)%b,prq=s5178,bef=cres(i1,i2,i3,&,
* i7,4)+,aft=
      cres(i1,i2,i3,1,i7,4)=cres(i1,i2,i3,1,i7,4)+(l5_718(i7,i1)
     & %a(1)*r6_324(i3,i2)%a(1)+l5_718(i7,i1)%c(1)*s5178*r6_324(
     & i3,i2)%b(2))
      cres(i1,i2,i3,2,i7,4)=cres(i1,i2,i3,2,i7,4)+(l5_718(i7,i1)
     & %c(2)*s5178*r6_324(i3,i2)%b(1)+l5_718(i7,i1)%a(2)*r6_324(
     & i3,i2)%a(2))
      end do
      end do
      end do
      end do
      endif
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,5),a1=l3_516(i5,i1)%a,c1=l3_516(i5,i1)%
* c,a2=r4_728(i7,i2)%a,b2=r4_728(i7,i2)%b,prq=s3156,bef=,aft=
      cres(i1,i2,1,i5,i7,5)=(l3_516(i5,i1)%a(1)*r4_728(i7,i2)%a(
     & 1)+l3_516(i5,i1)%c(1)*s3156*r4_728(i7,i2)%b(2))
      cres(i1,i2,2,i5,i7,5)=(l3_516(i5,i1)%c(2)*s3156*r4_728(i7,
     & i2)%b(1)+l3_516(i5,i1)%a(2)*r4_728(i7,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,5),a1=l3_728(i7,i2)%a,c1=l3_728(i7,i2)%
* c,a2=r4_516(i5,i1)%a,b2=r4_516(i5,i1)%b,prq=s3278,bef=cres(i1,i2,&,i5,
* i7,5)+,aft=
      cres(i1,i2,1,i5,i7,5)=cres(i1,i2,1,i5,i7,5)+(l3_728(i7,i2)
     & %a(1)*r4_516(i5,i1)%a(1)+l3_728(i7,i2)%c(1)*s3278*r4_516(
     & i5,i1)%b(2))
      cres(i1,i2,2,i5,i7,5)=cres(i1,i2,2,i5,i7,5)+(l3_728(i7,i2)
     & %c(2)*s3278*r4_516(i5,i1)%b(1)+l3_728(i7,i2)%a(2)*r4_516(
     & i5,i1)%a(2))
      end do
      end do
      end do
      end do
      endif
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,6),a1=l3_526(i5,i2)%a,c1=l3_526(i5,i2)%
* c,a2=r4_718(i7,i1)%a,b2=r4_718(i7,i1)%b,prq=s3256,bef=,aft=
      cres(i1,i2,1,i5,i7,6)=(l3_526(i5,i2)%a(1)*r4_718(i7,i1)%a(
     & 1)+l3_526(i5,i2)%c(1)*s3256*r4_718(i7,i1)%b(2))
      cres(i1,i2,2,i5,i7,6)=(l3_526(i5,i2)%c(2)*s3256*r4_718(i7,
     & i1)%b(1)+l3_526(i5,i2)%a(2)*r4_718(i7,i1)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,6),a1=l3_718(i7,i1)%a,c1=l3_718(i7,i1)%
* c,a2=r4_526(i5,i2)%a,b2=r4_526(i5,i2)%b,prq=s3178,bef=cres(i1,i2,&,i5,
* i7,6)+,aft=
      cres(i1,i2,1,i5,i7,6)=cres(i1,i2,1,i5,i7,6)+(l3_718(i7,i1)
     & %a(1)*r4_526(i5,i2)%a(1)+l3_718(i7,i1)%c(1)*s3178*r4_526(
     & i5,i2)%b(2))
      cres(i1,i2,2,i5,i7,6)=cres(i1,i2,2,i5,i7,6)+(l3_718(i7,i1)
     & %c(2)*s3178*r4_526(i5,i2)%b(1)+l3_718(i7,i1)%a(2)*r4_526(
     & i5,i2)%a(2))
      end do
      end do
      end do
      end do
      endif
  
* line without gluon; col 7-12                                          
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,7),a1=l7_3124(i3,i1,i2)%a,c1=l7_3124(i3
* ,i1,i2)%c,a2=r8_56(i5)%a,b2=r8_56(i5)%b,prq=s856,bef=,aft=
      cres(i1,i2,i3,i5,1,7)=(l7_3124(i3,i1,i2)%a(1)*r8_56(i5)%a(
     & 1)+l7_3124(i3,i1,i2)%c(1)*s856*r8_56(i5)%b(2))
      cres(i1,i2,i3,i5,2,7)=(l7_3124(i3,i1,i2)%c(2)*s856*r8_56(i
     & 5)%b(1)+l7_3124(i3,i1,i2)%a(2)*r8_56(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,7),a1=l7_56(i5)%a,c1=l7_56(i5)%c,a2=r8_
* 3124(i3,i1,i2)%a,b2=r8_3124(i3,i1,i2)%b,prq=s756,bef=cres(i1,i2,i3,i5,
* &,7)+,aft=
      cres(i1,i2,i3,i5,1,7)=cres(i1,i2,i3,i5,1,7)+(l7_56(i5)%a(1
     & )*r8_3124(i3,i1,i2)%a(1)+l7_56(i5)%c(1)*s756*r8_3124(i3,i
     & 1,i2)%b(2))
      cres(i1,i2,i3,i5,2,7)=cres(i1,i2,i3,i5,2,7)+(l7_56(i5)%c(2
     & )*s756*r8_3124(i3,i1,i2)%b(1)+l7_56(i5)%a(2)*r8_3124(i3,i
     & 1,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,7),a1=l5_3124(i3,i1,i2)%a,c1=l5_3124(i3
* ,i1,i2)%c,a2=r6_78(i7)%a,b2=r6_78(i7)%b,prq=s678,bef=cres(i1,i2,i3,&,i
* 7,7)+,aft=
      cres(i1,i2,i3,1,i7,7)=cres(i1,i2,i3,1,i7,7)+(l5_3124(i3,i1
     & ,i2)%a(1)*r6_78(i7)%a(1)+l5_3124(i3,i1,i2)%c(1)*s678*r6_7
     & 8(i7)%b(2))
      cres(i1,i2,i3,2,i7,7)=cres(i1,i2,i3,2,i7,7)+(l5_3124(i3,i1
     & ,i2)%c(2)*s678*r6_78(i7)%b(1)+l5_3124(i3,i1,i2)%a(2)*r6_7
     & 8(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,7),a1=l5_78(i7)%a,c1=l5_78(i7)%c,a2=r6_
* 3124(i3,i1,i2)%a,b2=r6_3124(i3,i1,i2)%b,prq=s578,bef=cres(i1,i2,i3,&,i
* 7,7)+,aft=
      cres(i1,i2,i3,1,i7,7)=cres(i1,i2,i3,1,i7,7)+(l5_78(i7)%a(1
     & )*r6_3124(i3,i1,i2)%a(1)+l5_78(i7)%c(1)*s578*r6_3124(i3,i
     & 1,i2)%b(2))
      cres(i1,i2,i3,2,i7,7)=cres(i1,i2,i3,2,i7,7)+(l5_78(i7)%c(2
     & )*s578*r6_3124(i3,i1,i2)%b(1)+l5_78(i7)%a(2)*r6_3124(i3,i
     & 1,i2)%a(2))
      end do
      end do
      end do
      end do
      endif
  
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,8),a1=l7_3214(i3,i1,i2)%a,c1=l7_3214(i3
* ,i1,i2)%c,a2=r8_56(i5)%a,b2=r8_56(i5)%b,prq=s856,bef=,aft=
      cres(i1,i2,i3,i5,1,8)=(l7_3214(i3,i1,i2)%a(1)*r8_56(i5)%a(
     & 1)+l7_3214(i3,i1,i2)%c(1)*s856*r8_56(i5)%b(2))
      cres(i1,i2,i3,i5,2,8)=(l7_3214(i3,i1,i2)%c(2)*s856*r8_56(i
     & 5)%b(1)+l7_3214(i3,i1,i2)%a(2)*r8_56(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,8),a1=l7_56(i5)%a,c1=l7_56(i5)%c,a2=r8_
* 3214(i3,i1,i2)%a,b2=r8_3214(i3,i1,i2)%b,prq=s756,bef=cres(i1,i2,i3,i5,
* &,8)+,aft=
      cres(i1,i2,i3,i5,1,8)=cres(i1,i2,i3,i5,1,8)+(l7_56(i5)%a(1
     & )*r8_3214(i3,i1,i2)%a(1)+l7_56(i5)%c(1)*s756*r8_3214(i3,i
     & 1,i2)%b(2))
      cres(i1,i2,i3,i5,2,8)=cres(i1,i2,i3,i5,2,8)+(l7_56(i5)%c(2
     & )*s756*r8_3214(i3,i1,i2)%b(1)+l7_56(i5)%a(2)*r8_3214(i3,i
     & 1,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,8),a1=l5_3214(i3,i1,i2)%a,c1=l5_3214(i3
* ,i1,i2)%c,a2=r6_78(i7)%a,b2=r6_78(i7)%b,prq=s678,bef=cres(i1,i2,i3,&,i
* 7,8)+,aft=
      cres(i1,i2,i3,1,i7,8)=cres(i1,i2,i3,1,i7,8)+(l5_3214(i3,i1
     & ,i2)%a(1)*r6_78(i7)%a(1)+l5_3214(i3,i1,i2)%c(1)*s678*r6_7
     & 8(i7)%b(2))
      cres(i1,i2,i3,2,i7,8)=cres(i1,i2,i3,2,i7,8)+(l5_3214(i3,i1
     & ,i2)%c(2)*s678*r6_78(i7)%b(1)+l5_3214(i3,i1,i2)%a(2)*r6_7
     & 8(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,8),a1=l5_78(i7)%a,c1=l5_78(i7)%c,a2=r6_
* 3214(i3,i1,i2)%a,b2=r6_3214(i3,i1,i2)%b,prq=s578,bef=cres(i1,i2,i3,&,i
* 7,8)+,aft=
      cres(i1,i2,i3,1,i7,8)=cres(i1,i2,i3,1,i7,8)+(l5_78(i7)%a(1
     & )*r6_3214(i3,i1,i2)%a(1)+l5_78(i7)%c(1)*s578*r6_3214(i3,i
     & 1,i2)%b(2))
      cres(i1,i2,i3,2,i7,8)=cres(i1,i2,i3,2,i7,8)+(l5_78(i7)%c(2
     & )*s578*r6_3214(i3,i1,i2)%b(1)+l5_78(i7)%a(2)*r6_3214(i3,i
     & 1,i2)%a(2))
      end do
      end do
      end do
      end do
      endif
  
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,9),a1=l7_5126(i5,i1,i2)%a,c1=l7_5126(i5
* ,i1,i2)%c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=,aft=
      cres(i1,i2,i3,i5,1,9)=(l7_5126(i5,i1,i2)%a(1)*r8_34(i3)%a(
     & 1)+l7_5126(i5,i1,i2)%c(1)*s834*r8_34(i3)%b(2))
      cres(i1,i2,i3,i5,2,9)=(l7_5126(i5,i1,i2)%c(2)*s834*r8_34(i
     & 3)%b(1)+l7_5126(i5,i1,i2)%a(2)*r8_34(i3)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,9),a1=l7_34(i3)%a,c1=l7_34(i3)%c,a2=r8_
* 5126(i5,i1,i2)%a,b2=r8_5126(i5,i1,i2)%b,prq=s734,bef=cres(i1,i2,i3,i5,
* &,9)+,aft=
      cres(i1,i2,i3,i5,1,9)=cres(i1,i2,i3,i5,1,9)+(l7_34(i3)%a(1
     & )*r8_5126(i5,i1,i2)%a(1)+l7_34(i3)%c(1)*s734*r8_5126(i5,i
     & 1,i2)%b(2))
      cres(i1,i2,i3,i5,2,9)=cres(i1,i2,i3,i5,2,9)+(l7_34(i3)%c(2
     & )*s734*r8_5126(i5,i1,i2)%b(1)+l7_34(i3)%a(2)*r8_5126(i5,i
     & 1,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,9),a1=l3_5126(i5,i1,i2)%a,c1=l3_5126(i5
* ,i1,i2)%c,a2=r4_78(i7)%a,b2=r4_78(i7)%b,prq=s478,bef=cres(i1,i2,&,i5,i
* 7,9)+,aft=
      cres(i1,i2,1,i5,i7,9)=cres(i1,i2,1,i5,i7,9)+(l3_5126(i5,i1
     & ,i2)%a(1)*r4_78(i7)%a(1)+l3_5126(i5,i1,i2)%c(1)*s478*r4_7
     & 8(i7)%b(2))
      cres(i1,i2,2,i5,i7,9)=cres(i1,i2,2,i5,i7,9)+(l3_5126(i5,i1
     & ,i2)%c(2)*s478*r4_78(i7)%b(1)+l3_5126(i5,i1,i2)%a(2)*r4_7
     & 8(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,9),a1=l3_78(i7)%a,c1=l3_78(i7)%c,a2=r4_
* 5126(i5,i1,i2)%a,b2=r4_5126(i5,i1,i2)%b,prq=s378,bef=cres(i1,i2,&,i5,i
* 7,9)+,aft=
      cres(i1,i2,1,i5,i7,9)=cres(i1,i2,1,i5,i7,9)+(l3_78(i7)%a(1
     & )*r4_5126(i5,i1,i2)%a(1)+l3_78(i7)%c(1)*s378*r4_5126(i5,i
     & 1,i2)%b(2))
      cres(i1,i2,2,i5,i7,9)=cres(i1,i2,2,i5,i7,9)+(l3_78(i7)%c(2
     & )*s378*r4_5126(i5,i1,i2)%b(1)+l3_78(i7)%a(2)*r4_5126(i5,i
     & 1,i2)%a(2))
      end do
      end do
      end do
      end do
      endif
  
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,10),a1=l7_5216(i5,i1,i2)%a,c1=l7_5216(i
* 5,i1,i2)%c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=,aft=
      cres(i1,i2,i3,i5,1,10)=(l7_5216(i5,i1,i2)%a(1)*r8_34(i3)%a
     & (1)+l7_5216(i5,i1,i2)%c(1)*s834*r8_34(i3)%b(2))
      cres(i1,i2,i3,i5,2,10)=(l7_5216(i5,i1,i2)%c(2)*s834*r8_34(
     & i3)%b(1)+l7_5216(i5,i1,i2)%a(2)*r8_34(i3)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,10),a1=l7_34(i3)%a,c1=l7_34(i3)%c,a2=r8
* _5216(i5,i1,i2)%a,b2=r8_5216(i5,i1,i2)%b,prq=s734,bef=cres(i1,i2,i3,i5
* ,&,10)+,aft=
      cres(i1,i2,i3,i5,1,10)=cres(i1,i2,i3,i5,1,10)+(l7_34(i3)%a
     & (1)*r8_5216(i5,i1,i2)%a(1)+l7_34(i3)%c(1)*s734*r8_5216(i5
     & ,i1,i2)%b(2))
      cres(i1,i2,i3,i5,2,10)=cres(i1,i2,i3,i5,2,10)+(l7_34(i3)%c
     & (2)*s734*r8_5216(i5,i1,i2)%b(1)+l7_34(i3)%a(2)*r8_5216(i5
     & ,i1,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,10),a1=l3_5216(i5,i1,i2)%a,c1=l3_5216(i
* 5,i1,i2)%c,a2=r4_78(i7)%a,b2=r4_78(i7)%b,prq=s478,bef=cres(i1,i2,&,i5,
* i7,10)+,aft=
      cres(i1,i2,1,i5,i7,10)=cres(i1,i2,1,i5,i7,10)+(l3_5216(i5,
     & i1,i2)%a(1)*r4_78(i7)%a(1)+l3_5216(i5,i1,i2)%c(1)*s478*r4
     & _78(i7)%b(2))
      cres(i1,i2,2,i5,i7,10)=cres(i1,i2,2,i5,i7,10)+(l3_5216(i5,
     & i1,i2)%c(2)*s478*r4_78(i7)%b(1)+l3_5216(i5,i1,i2)%a(2)*r4
     & _78(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,10),a1=l3_78(i7)%a,c1=l3_78(i7)%c,a2=r4
* _5216(i5,i1,i2)%a,b2=r4_5216(i5,i1,i2)%b,prq=s378,bef=cres(i1,i2,&,i5,
* i7,10)+,aft=
      cres(i1,i2,1,i5,i7,10)=cres(i1,i2,1,i5,i7,10)+(l3_78(i7)%a
     & (1)*r4_5216(i5,i1,i2)%a(1)+l3_78(i7)%c(1)*s378*r4_5216(i5
     & ,i1,i2)%b(2))
      cres(i1,i2,2,i5,i7,10)=cres(i1,i2,2,i5,i7,10)+(l3_78(i7)%c
     & (2)*s378*r4_5216(i5,i1,i2)%b(1)+l3_78(i7)%a(2)*r4_5216(i5
     & ,i1,i2)%a(2))
      end do
      end do
      end do
      end do
      endif
  
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,11),a1=l5_7128(i7,i1,i2)%a,c1=l5_7128(i
* 7,i1,i2)%c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=,aft=
      cres(i1,i2,i3,1,i7,11)=(l5_7128(i7,i1,i2)%a(1)*r6_34(i3)%a
     & (1)+l5_7128(i7,i1,i2)%c(1)*s634*r6_34(i3)%b(2))
      cres(i1,i2,i3,2,i7,11)=(l5_7128(i7,i1,i2)%c(2)*s634*r6_34(
     & i3)%b(1)+l5_7128(i7,i1,i2)%a(2)*r6_34(i3)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,11),a1=l5_34(i3)%a,c1=l5_34(i3)%c,a2=r6
* _7128(i7,i1,i2)%a,b2=r6_7128(i7,i1,i2)%b,prq=s534,bef=cres(i1,i2,i3,&,
* i7,11)+,aft=
      cres(i1,i2,i3,1,i7,11)=cres(i1,i2,i3,1,i7,11)+(l5_34(i3)%a
     & (1)*r6_7128(i7,i1,i2)%a(1)+l5_34(i3)%c(1)*s534*r6_7128(i7
     & ,i1,i2)%b(2))
      cres(i1,i2,i3,2,i7,11)=cres(i1,i2,i3,2,i7,11)+(l5_34(i3)%c
     & (2)*s534*r6_7128(i7,i1,i2)%b(1)+l5_34(i3)%a(2)*r6_7128(i7
     & ,i1,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,11),a1=l3_7128(i7,i1,i2)%a,c1=l3_7128(i
* 7,i1,i2)%c,a2=r4_56(i5)%a,b2=r4_56(i5)%b,prq=s456,bef=cres(i1,i2,&,i5,
* i7,11)+,aft=
      cres(i1,i2,1,i5,i7,11)=cres(i1,i2,1,i5,i7,11)+(l3_7128(i7,
     & i1,i2)%a(1)*r4_56(i5)%a(1)+l3_7128(i7,i1,i2)%c(1)*s456*r4
     & _56(i5)%b(2))
      cres(i1,i2,2,i5,i7,11)=cres(i1,i2,2,i5,i7,11)+(l3_7128(i7,
     & i1,i2)%c(2)*s456*r4_56(i5)%b(1)+l3_7128(i7,i1,i2)%a(2)*r4
     & _56(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,11),a1=l3_56(i5)%a,c1=l3_56(i5)%c,a2=r4
* _7128(i7,i1,i2)%a,b2=r4_7128(i7,i1,i2)%b,prq=s356,bef=cres(i1,i2,&,i5,
* i7,11)+,aft=
      cres(i1,i2,1,i5,i7,11)=cres(i1,i2,1,i5,i7,11)+(l3_56(i5)%a
     & (1)*r4_7128(i7,i1,i2)%a(1)+l3_56(i5)%c(1)*s356*r4_7128(i7
     & ,i1,i2)%b(2))
      cres(i1,i2,2,i5,i7,11)=cres(i1,i2,2,i5,i7,11)+(l3_56(i5)%c
     & (2)*s356*r4_7128(i7,i1,i2)%b(1)+l3_56(i5)%a(2)*r4_7128(i7
     & ,i1,i2)%a(2))
      end do
      end do
      end do
      end do
      endif
  
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,12),a1=l5_7218(i7,i1,i2)%a,c1=l5_7218(i
* 7,i1,i2)%c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=,aft=
      cres(i1,i2,i3,1,i7,12)=(l5_7218(i7,i1,i2)%a(1)*r6_34(i3)%a
     & (1)+l5_7218(i7,i1,i2)%c(1)*s634*r6_34(i3)%b(2))
      cres(i1,i2,i3,2,i7,12)=(l5_7218(i7,i1,i2)%c(2)*s634*r6_34(
     & i3)%b(1)+l5_7218(i7,i1,i2)%a(2)*r6_34(i3)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,12),a1=l5_34(i3)%a,c1=l5_34(i3)%c,a2=r6
* _7218(i7,i1,i2)%a,b2=r6_7218(i7,i1,i2)%b,prq=s534,bef=cres(i1,i2,i3,&,
* i7,12)+,aft=
      cres(i1,i2,i3,1,i7,12)=cres(i1,i2,i3,1,i7,12)+(l5_34(i3)%a
     & (1)*r6_7218(i7,i1,i2)%a(1)+l5_34(i3)%c(1)*s534*r6_7218(i7
     & ,i1,i2)%b(2))
      cres(i1,i2,i3,2,i7,12)=cres(i1,i2,i3,2,i7,12)+(l5_34(i3)%c
     & (2)*s534*r6_7218(i7,i1,i2)%b(1)+l5_34(i3)%a(2)*r6_7218(i7
     & ,i1,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,12),a1=l3_7218(i7,i1,i2)%a,c1=l3_7218(i
* 7,i1,i2)%c,a2=r4_56(i5)%a,b2=r4_56(i5)%b,prq=s456,bef=cres(i1,i2,&,i5,
* i7,12)+,aft=
      cres(i1,i2,1,i5,i7,12)=cres(i1,i2,1,i5,i7,12)+(l3_7218(i7,
     & i1,i2)%a(1)*r4_56(i5)%a(1)+l3_7218(i7,i1,i2)%c(1)*s456*r4
     & _56(i5)%b(2))
      cres(i1,i2,2,i5,i7,12)=cres(i1,i2,2,i5,i7,12)+(l3_7218(i7,
     & i1,i2)%c(2)*s456*r4_56(i5)%b(1)+l3_7218(i7,i1,i2)%a(2)*r4
     & _56(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,12),a1=l3_56(i5)%a,c1=l3_56(i5)%c,a2=r4
* _7218(i7,i1,i2)%a,b2=r4_7218(i7,i1,i2)%b,prq=s356,bef=cres(i1,i2,&,i5,
* i7,12)+,aft=
      cres(i1,i2,1,i5,i7,12)=cres(i1,i2,1,i5,i7,12)+(l3_56(i5)%a
     & (1)*r4_7218(i7,i1,i2)%a(1)+l3_56(i5)%c(1)*s356*r4_7218(i7
     & ,i1,i2)%b(2))
      cres(i1,i2,2,i5,i7,12)=cres(i1,i2,2,i5,i7,12)+(l3_56(i5)%c
     & (2)*s356*r4_7218(i7,i1,i2)%b(1)+l3_56(i5)%a(2)*r4_7218(i7
     & ,i1,i2)%a(2))
      end do
      end do
      end do
      end do
      endif
  
  
  
* -> lines with one gluon                                               
  
      if (ilept(id3).ne.1) then
* I                                                                     
* quqd -- p=p31,q=p3156
      quqd=p31(0)*p3156(0)-p31(1)*p3156(1)-p31(2)*p3156(2)-p31(3
     & )*p3156(3)
      ccr=zcr(id3)/(f3156)
      ccl=zcl(id3)/(f3156)
      do i5=1,2
* T0 -- qu=p31,qd=p3156,v=cz56(i5)%e,a=u31_56(i5)%a,b=u31_56(i5)%b,c=u31
* _56(i5)%c,d=u31_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p31(2)*p3156(3)-p3156(2)*p31(3))+p31
     & k0*(cz56(i5)%e(2)*p3156(3)-p3156(2)*cz56(i5)%e(3))-p3156k
     & 0*(cz56(i5)%e(2)*p31(3)-p31(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p31k0+p31(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p3156k0+p3156(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p31(0)-cz56(i5)%e(1)*p31(1)-cz56(i5)%e(
     & 2)*p31(2)-cz56(i5)%e(3)*p31(3)
      cvqd=cz56(i5)%e(0)*p3156(0)-cz56(i5)%e(1)*p3156(1)-cz56(i5
     & )%e(2)*p3156(2)-cz56(i5)%e(3)*p3156(3)
      cauxa=-cz56(i5)%ek0*quqd+p31k0*cvqd+p3156k0*cvqu
      cauxb=-cz56(i5)%ek0*p3156(2)+p3156k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p31(2)-p31k0*cz56(i5)%e(2)
      u31_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u31_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u31_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u31_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u31_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u31_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u31_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u31_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f3156)
      ccl=fcl(id3)/(f3156)
      do i5=1,2
* T0 -- qu=p31,qd=p3156,v=cf56(i5)%e,a=u31_56(i5)%a,b=u31_56(i5)%b,c=u31
* _56(i5)%c,d=u31_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p31(2)*p3156(3)-p3156(2)*p31(3))+p31
     & k0*(cf56(i5)%e(2)*p3156(3)-p3156(2)*cf56(i5)%e(3))-p3156k
     & 0*(cf56(i5)%e(2)*p31(3)-p31(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p31k0+p31(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p3156k0+p3156(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p31(0)-cf56(i5)%e(1)*p31(1)-cf56(i5)%e(
     & 2)*p31(2)-cf56(i5)%e(3)*p31(3)
      cvqd=cf56(i5)%e(0)*p3156(0)-cf56(i5)%e(1)*p3156(1)-cf56(i5
     & )%e(2)*p3156(2)-cf56(i5)%e(3)*p3156(3)
      cauxa=-cf56(i5)%ek0*quqd+p31k0*cvqd+p3156k0*cvqu
      cauxb=-cf56(i5)%ek0*p3156(2)+p3156k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p31(2)-p31k0*cf56(i5)%e(2)
      u31_56(i5)%a(1)=u31_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u31_56(i5)%a(2)=u31_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u31_56(i5)%b(1)=u31_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u31_56(i5)%b(2)=u31_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u31_56(i5)%c(1)=u31_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u31_56(i5)%c(2)=u31_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u31_56(i5)%d(1)=u31_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u31_56(i5)%d(2)=u31_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l3_156(i1,i5)%a,cc=l3_156(i1,i5)%c,a1=l3_1(i1)%a,c1=l3_1(i1
* )%c,a2=u31_56(i5)%a,b2=u31_56(i5)%b,c2=u31_56(i5)%c,d2=u31_56(i5)%d,pr
* q=s31,nsum=0
      l3_156(i1,i5)%a(1)=l3_1(i1)%a(1)*u31_56(i5)%a(1)+l3_1(i1)%
     & c(1)*s31*u31_56(i5)%b(2)
      l3_156(i1,i5)%c(1)=l3_1(i1)%a(1)*u31_56(i5)%c(1)+l3_1(i1)%
     & c(1)*s31*u31_56(i5)%d(2)
      l3_156(i1,i5)%c(2)=l3_1(i1)%c(2)*s31*u31_56(i5)%d(1)+l3_1(
     & i1)%a(2)*u31_56(i5)%c(2)
      l3_156(i1,i5)%a(2)=l3_1(i1)%c(2)*s31*u31_56(i5)%b(1)+l3_1(
     & i1)%a(2)*u31_56(i5)%a(2)
      end do
      end do
  
* quqd -- p=p356,q=p3156
      quqd=p356(0)*p3156(0)-p356(1)*p3156(1)-p356(2)*p3156(2)-p3
     & 56(3)*p3156(3)
      ccr=1.d0/(f3156)
      ccl=1.d0/(f3156)
      do i1=1,2
* T0 -- qu=p356,qd=p3156,v=ce1(i1)%e,a=u356_1(i1)%a,b=u356_1(i1)%b,c=u35
* 6_1(i1)%c,d=u356_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p356(2)*p3156(3)-p3156(2)*p356(3))+p3
     & 56k0*(ce1(i1)%e(2)*p3156(3)-p3156(2)*ce1(i1)%e(3))-p3156k
     & 0*(ce1(i1)%e(2)*p356(3)-p356(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p356k0+p356(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p3156k0+p3156(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p356(0)-ce1(i1)%e(1)*p356(1)-ce1(i1)%e(2
     & )*p356(2)-ce1(i1)%e(3)*p356(3)
      cvqd=ce1(i1)%e(0)*p3156(0)-ce1(i1)%e(1)*p3156(1)-ce1(i1)%e
     & (2)*p3156(2)-ce1(i1)%e(3)*p3156(3)
      cauxa=-ce1(i1)%ek0*quqd+p356k0*cvqd+p3156k0*cvqu
      cauxb=-ce1(i1)%ek0*p3156(2)+p3156k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p356(2)-p356k0*ce1(i1)%e(2)
      u356_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u356_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u356_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u356_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u356_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u356_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u356_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u356_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l3_156(i1,i5)%a,cc=l3_156(i1,i5)%c,a1=l3_56(i5)%a,c1=l3_56(
* i5)%c,a2=u356_1(i1)%a,b2=u356_1(i1)%b,c2=u356_1(i1)%c,d2=u356_1(i1)%d,
* prq=s356,nsum=1
      l3_156(i1,i5)%a(1)=l3_156(i1,i5)%a(1)+l3_56(i5)%a(1)*u356_
     & 1(i1)%a(1)+l3_56(i5)%c(1)*s356*u356_1(i1)%b(2)
      l3_156(i1,i5)%c(1)=l3_156(i1,i5)%c(1)+l3_56(i5)%a(1)*u356_
     & 1(i1)%c(1)+l3_56(i5)%c(1)*s356*u356_1(i1)%d(2)
      l3_156(i1,i5)%c(2)=l3_156(i1,i5)%c(2)+l3_56(i5)%c(2)*s356*
     & u356_1(i1)%d(1)+l3_56(i5)%a(2)*u356_1(i1)%c(2)
      l3_156(i1,i5)%a(2)=l3_156(i1,i5)%a(2)+l3_56(i5)%c(2)*s356*
     & u356_1(i1)%b(1)+l3_56(i5)%a(2)*u356_1(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id7).ne.1) then
* quqd -- p=p31,q=p456
      quqd=p31(0)*p456(0)-p31(1)*p456(1)-p31(2)*p456(2)-p31(3)*p
     & 456(3)
      ccr=zcr(id3)/(f456)
      ccl=zcl(id3)/(f456)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p31,qd=p456,v=cz728(i7,i2)%e,a=u31_728(i7,i2)%a,b=u31_728(i7,
* i2)%b,c=u31_728(i7,i2)%c,d=u31_728(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz728(i7,i2)%ek0*(p31(2)*p456(3)-p456(2)*p31(3))+p
     & 31k0*(cz728(i7,i2)%e(2)*p456(3)-p456(2)*cz728(i7,i2)%e(3)
     & )-p456k0*(cz728(i7,i2)%e(2)*p31(3)-p31(2)*cz728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz728(i7,i2)%e(3)*p31k0+p31(3)*cz728(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz728(i7,i2)%e(3)*p456k0+p456(3)*cz728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz728(i7,i2)%e(0)*p31(0)-cz728(i7,i2)%e(1)*p31(1)-cz7
     & 28(i7,i2)%e(2)*p31(2)-cz728(i7,i2)%e(3)*p31(3)
      cvqd=cz728(i7,i2)%e(0)*p456(0)-cz728(i7,i2)%e(1)*p456(1)-c
     & z728(i7,i2)%e(2)*p456(2)-cz728(i7,i2)%e(3)*p456(3)
      cauxa=-cz728(i7,i2)%ek0*quqd+p31k0*cvqd+p456k0*cvqu
      cauxb=-cz728(i7,i2)%ek0*p456(2)+p456k0*cz728(i7,i2)%e(2)
      cauxc=+cz728(i7,i2)%ek0*p31(2)-p31k0*cz728(i7,i2)%e(2)
      u31_728(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u31_728(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u31_728(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u31_728(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u31_728(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u31_728(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u31_728(i7,i2)%d(1)=ccl*cz728(i7,i2)%ek0
      u31_728(i7,i2)%d(2)=ccr*cz728(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id3)/(f456)
      ccl=fcl(id3)/(f456)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p31,qd=p456,v=cf728(i7,i2)%e,a=u31_728(i7,i2)%a,b=u31_728(i7,
* i2)%b,c=u31_728(i7,i2)%c,d=u31_728(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf728(i7,i2)%ek0*(p31(2)*p456(3)-p456(2)*p31(3))+p
     & 31k0*(cf728(i7,i2)%e(2)*p456(3)-p456(2)*cf728(i7,i2)%e(3)
     & )-p456k0*(cf728(i7,i2)%e(2)*p31(3)-p31(2)*cf728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf728(i7,i2)%e(3)*p31k0+p31(3)*cf728(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf728(i7,i2)%e(3)*p456k0+p456(3)*cf728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf728(i7,i2)%e(0)*p31(0)-cf728(i7,i2)%e(1)*p31(1)-cf7
     & 28(i7,i2)%e(2)*p31(2)-cf728(i7,i2)%e(3)*p31(3)
      cvqd=cf728(i7,i2)%e(0)*p456(0)-cf728(i7,i2)%e(1)*p456(1)-c
     & f728(i7,i2)%e(2)*p456(2)-cf728(i7,i2)%e(3)*p456(3)
      cauxa=-cf728(i7,i2)%ek0*quqd+p31k0*cvqd+p456k0*cvqu
      cauxb=-cf728(i7,i2)%ek0*p456(2)+p456k0*cf728(i7,i2)%e(2)
      cauxc=+cf728(i7,i2)%ek0*p31(2)-p31k0*cf728(i7,i2)%e(2)
      u31_728(i7,i2)%a(1)=u31_728(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u31_728(i7,i2)%a(2)=u31_728(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u31_728(i7,i2)%b(1)=u31_728(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u31_728(i7,i2)%b(2)=u31_728(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u31_728(i7,i2)%c(1)=u31_728(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u31_728(i7,i2)%c(2)=u31_728(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u31_728(i7,i2)%d(1)=u31_728(i7,i2)%d(1)+ccl*cf728(i7,i2)%e
     & k0
      u31_728(i7,i2)%d(2)=u31_728(i7,i2)%d(2)+ccr*cf728(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_1728(i1,i7,i2)%a,cc=l3_1728(i1,i7,i2)%c,a1=l3_1(i1)%a,c1
* =l3_1(i1)%c,a2=u31_728(i7,i2)%a,b2=u31_728(i7,i2)%b,c2=u31_728(i7,i2)%
* c,d2=u31_728(i7,i2)%d,prq=s31,nsum=0
      l3_1728(i1,i7,i2)%a(1)=l3_1(i1)%a(1)*u31_728(i7,i2)%a(1)+l
     & 3_1(i1)%c(1)*s31*u31_728(i7,i2)%b(2)
      l3_1728(i1,i7,i2)%c(1)=l3_1(i1)%a(1)*u31_728(i7,i2)%c(1)+l
     & 3_1(i1)%c(1)*s31*u31_728(i7,i2)%d(2)
      l3_1728(i1,i7,i2)%c(2)=l3_1(i1)%c(2)*s31*u31_728(i7,i2)%d(
     & 1)+l3_1(i1)%a(2)*u31_728(i7,i2)%c(2)
      l3_1728(i1,i7,i2)%a(2)=l3_1(i1)%c(2)*s31*u31_728(i7,i2)%b(
     & 1)+l3_1(i1)%a(2)*u31_728(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u3278_1                                                
* quqd -- p=p3278,q=p456
      quqd=p3278(0)*p456(0)-p3278(1)*p456(1)-p3278(2)*p456(2)-p3
     & 278(3)*p456(3)
      ccr=1.d0/(f456)
      ccl=1.d0/(f456)
      do i1=1,2
* T0 -- qu=p3278,qd=p456,v=ce1(i1)%e,a=u3728_1(i1)%a,b=u3728_1(i1)%b,c=u
* 3728_1(i1)%c,d=u3728_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p3278(2)*p456(3)-p456(2)*p3278(3))+p3
     & 278k0*(ce1(i1)%e(2)*p456(3)-p456(2)*ce1(i1)%e(3))-p456k0*
     & (ce1(i1)%e(2)*p3278(3)-p3278(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p3278k0+p3278(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p456k0+p456(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p3278(0)-ce1(i1)%e(1)*p3278(1)-ce1(i1)%e
     & (2)*p3278(2)-ce1(i1)%e(3)*p3278(3)
      cvqd=ce1(i1)%e(0)*p456(0)-ce1(i1)%e(1)*p456(1)-ce1(i1)%e(2
     & )*p456(2)-ce1(i1)%e(3)*p456(3)
      cauxa=-ce1(i1)%ek0*quqd+p3278k0*cvqd+p456k0*cvqu
      cauxb=-ce1(i1)%ek0*p456(2)+p456k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p3278(2)-p3278k0*ce1(i1)%e(2)
      u3728_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u3728_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u3728_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u3728_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u3728_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u3728_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u3728_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u3728_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_1728(i1,i7,i2)%a,cc=l3_1728(i1,i7,i2)%c,a1=l3_728(i7,i2)
* %a,c1=l3_728(i7,i2)%c,a2=u3728_1(i1)%a,b2=u3728_1(i1)%b,c2=u3728_1(i1)
* %c,d2=u3728_1(i1)%d,prq=s3278,nsum=1
      l3_1728(i1,i7,i2)%a(1)=l3_1728(i1,i7,i2)%a(1)+l3_728(i7,i2
     & )%a(1)*u3728_1(i1)%a(1)+l3_728(i7,i2)%c(1)*s3278*u3728_1(
     & i1)%b(2)
      l3_1728(i1,i7,i2)%c(1)=l3_1728(i1,i7,i2)%c(1)+l3_728(i7,i2
     & )%a(1)*u3728_1(i1)%c(1)+l3_728(i7,i2)%c(1)*s3278*u3728_1(
     & i1)%d(2)
      l3_1728(i1,i7,i2)%c(2)=l3_1728(i1,i7,i2)%c(2)+l3_728(i7,i2
     & )%c(2)*s3278*u3728_1(i1)%d(1)+l3_728(i7,i2)%a(2)*u3728_1(
     & i1)%c(2)
      l3_1728(i1,i7,i2)%a(2)=l3_1728(i1,i7,i2)%a(2)+l3_728(i7,i2
     & )%c(2)*s3278*u3728_1(i1)%b(1)+l3_728(i7,i2)%a(2)*u3728_1(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p356,q=p41
      quqd=p356(0)*p41(0)-p356(1)*p41(1)-p356(2)*p41(2)-p356(3)*
     & p41(3)
      ccr=zcr(id3)/(f41)
      ccl=zcl(id3)/(f41)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p356,qd=p41,v=cz728(i7,i2)%e,a=u356_728(i7,i2)%a,b=u356_728(i
* 7,i2)%b,c=u356_728(i7,i2)%c,d=u356_728(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz728(i7,i2)%ek0*(p356(2)*p41(3)-p41(2)*p356(3))+p
     & 356k0*(cz728(i7,i2)%e(2)*p41(3)-p41(2)*cz728(i7,i2)%e(3))
     & -p41k0*(cz728(i7,i2)%e(2)*p356(3)-p356(2)*cz728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz728(i7,i2)%e(3)*p356k0+p356(3)*cz728(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz728(i7,i2)%e(3)*p41k0+p41(3)*cz728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz728(i7,i2)%e(0)*p356(0)-cz728(i7,i2)%e(1)*p356(1)-c
     & z728(i7,i2)%e(2)*p356(2)-cz728(i7,i2)%e(3)*p356(3)
      cvqd=cz728(i7,i2)%e(0)*p41(0)-cz728(i7,i2)%e(1)*p41(1)-cz7
     & 28(i7,i2)%e(2)*p41(2)-cz728(i7,i2)%e(3)*p41(3)
      cauxa=-cz728(i7,i2)%ek0*quqd+p356k0*cvqd+p41k0*cvqu
      cauxb=-cz728(i7,i2)%ek0*p41(2)+p41k0*cz728(i7,i2)%e(2)
      cauxc=+cz728(i7,i2)%ek0*p356(2)-p356k0*cz728(i7,i2)%e(2)
      u356_728(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u356_728(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u356_728(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u356_728(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u356_728(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u356_728(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u356_728(i7,i2)%d(1)=ccl*cz728(i7,i2)%ek0
      u356_728(i7,i2)%d(2)=ccr*cz728(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id3)/(f41)
      ccl=fcl(id3)/(f41)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p356,qd=p41,v=cf728(i7,i2)%e,a=u356_728(i7,i2)%a,b=u356_728(i
* 7,i2)%b,c=u356_728(i7,i2)%c,d=u356_728(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf728(i7,i2)%ek0*(p356(2)*p41(3)-p41(2)*p356(3))+p
     & 356k0*(cf728(i7,i2)%e(2)*p41(3)-p41(2)*cf728(i7,i2)%e(3))
     & -p41k0*(cf728(i7,i2)%e(2)*p356(3)-p356(2)*cf728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf728(i7,i2)%e(3)*p356k0+p356(3)*cf728(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf728(i7,i2)%e(3)*p41k0+p41(3)*cf728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf728(i7,i2)%e(0)*p356(0)-cf728(i7,i2)%e(1)*p356(1)-c
     & f728(i7,i2)%e(2)*p356(2)-cf728(i7,i2)%e(3)*p356(3)
      cvqd=cf728(i7,i2)%e(0)*p41(0)-cf728(i7,i2)%e(1)*p41(1)-cf7
     & 28(i7,i2)%e(2)*p41(2)-cf728(i7,i2)%e(3)*p41(3)
      cauxa=-cf728(i7,i2)%ek0*quqd+p356k0*cvqd+p41k0*cvqu
      cauxb=-cf728(i7,i2)%ek0*p41(2)+p41k0*cf728(i7,i2)%e(2)
      cauxc=+cf728(i7,i2)%ek0*p356(2)-p356k0*cf728(i7,i2)%e(2)
      u356_728(i7,i2)%a(1)=u356_728(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u356_728(i7,i2)%a(2)=u356_728(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u356_728(i7,i2)%b(1)=u356_728(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u356_728(i7,i2)%b(2)=u356_728(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u356_728(i7,i2)%c(1)=u356_728(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u356_728(i7,i2)%c(2)=u356_728(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u356_728(i7,i2)%d(1)=u356_728(i7,i2)%d(1)+ccl*cf728(i7,i2)
     & %ek0
      u356_728(i7,i2)%d(2)=u356_728(i7,i2)%d(2)+ccr*cf728(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_56728(i5,i7,i2)%a,cc=l3_56728(i5,i7,i2)%c,a1=l3_56(i5)%a
* ,c1=l3_56(i5)%c,a2=u356_728(i7,i2)%a,b2=u356_728(i7,i2)%b,c2=u356_728(
* i7,i2)%c,d2=u356_728(i7,i2)%d,prq=s356,nsum=0
      l3_56728(i5,i7,i2)%a(1)=l3_56(i5)%a(1)*u356_728(i7,i2)%a(1
     & )+l3_56(i5)%c(1)*s356*u356_728(i7,i2)%b(2)
      l3_56728(i5,i7,i2)%c(1)=l3_56(i5)%a(1)*u356_728(i7,i2)%c(1
     & )+l3_56(i5)%c(1)*s356*u356_728(i7,i2)%d(2)
      l3_56728(i5,i7,i2)%c(2)=l3_56(i5)%c(2)*s356*u356_728(i7,i2
     & )%d(1)+l3_56(i5)%a(2)*u356_728(i7,i2)%c(2)
      l3_56728(i5,i7,i2)%a(2)=l3_56(i5)%c(2)*s356*u356_728(i7,i2
     & )%b(1)+l3_56(i5)%a(2)*u356_728(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u3278_56                                               
* quqd -- p=p3278,q=p41
      quqd=p3278(0)*p41(0)-p3278(1)*p41(1)-p3278(2)*p41(2)-p3278
     & (3)*p41(3)
      ccr=zcr(id3)/(f41)
      ccl=zcl(id3)/(f41)
      do i5=1,2
* T0 -- qu=p3278,qd=p41,v=cz56(i5)%e,a=u3728_56(i5)%a,b=u3728_56(i5)%b,c
* =u3728_56(i5)%c,d=u3728_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p3278(2)*p41(3)-p41(2)*p3278(3))+p32
     & 78k0*(cz56(i5)%e(2)*p41(3)-p41(2)*cz56(i5)%e(3))-p41k0*(c
     & z56(i5)%e(2)*p3278(3)-p3278(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p3278k0+p3278(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p41k0+p41(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p3278(0)-cz56(i5)%e(1)*p3278(1)-cz56(i5
     & )%e(2)*p3278(2)-cz56(i5)%e(3)*p3278(3)
      cvqd=cz56(i5)%e(0)*p41(0)-cz56(i5)%e(1)*p41(1)-cz56(i5)%e(
     & 2)*p41(2)-cz56(i5)%e(3)*p41(3)
      cauxa=-cz56(i5)%ek0*quqd+p3278k0*cvqd+p41k0*cvqu
      cauxb=-cz56(i5)%ek0*p41(2)+p41k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p3278(2)-p3278k0*cz56(i5)%e(2)
      u3728_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u3728_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u3728_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u3728_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u3728_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u3728_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u3728_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u3728_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f41)
      ccl=fcl(id3)/(f41)
      do i5=1,2
* T0 -- qu=p3278,qd=p41,v=cf56(i5)%e,a=u3728_56(i5)%a,b=u3728_56(i5)%b,c
* =u3728_56(i5)%c,d=u3728_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p3278(2)*p41(3)-p41(2)*p3278(3))+p32
     & 78k0*(cf56(i5)%e(2)*p41(3)-p41(2)*cf56(i5)%e(3))-p41k0*(c
     & f56(i5)%e(2)*p3278(3)-p3278(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p3278k0+p3278(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p41k0+p41(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p3278(0)-cf56(i5)%e(1)*p3278(1)-cf56(i5
     & )%e(2)*p3278(2)-cf56(i5)%e(3)*p3278(3)
      cvqd=cf56(i5)%e(0)*p41(0)-cf56(i5)%e(1)*p41(1)-cf56(i5)%e(
     & 2)*p41(2)-cf56(i5)%e(3)*p41(3)
      cauxa=-cf56(i5)%ek0*quqd+p3278k0*cvqd+p41k0*cvqu
      cauxb=-cf56(i5)%ek0*p41(2)+p41k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p3278(2)-p3278k0*cf56(i5)%e(2)
      u3728_56(i5)%a(1)=u3728_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u3728_56(i5)%a(2)=u3728_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u3728_56(i5)%b(1)=u3728_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u3728_56(i5)%b(2)=u3728_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u3728_56(i5)%c(1)=u3728_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u3728_56(i5)%c(2)=u3728_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u3728_56(i5)%d(1)=u3728_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u3728_56(i5)%d(2)=u3728_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_56728(i5,i7,i2)%a,cc=l3_56728(i5,i7,i2)%c,a1=l3_728(i7,i
* 2)%a,c1=l3_728(i7,i2)%c,a2=u3728_56(i5)%a,b2=u3728_56(i5)%b,c2=u3728_5
* 6(i5)%c,d2=u3728_56(i5)%d,prq=s3278,nsum=1
      l3_56728(i5,i7,i2)%a(1)=l3_56728(i5,i7,i2)%a(1)+l3_728(i7,
     & i2)%a(1)*u3728_56(i5)%a(1)+l3_728(i7,i2)%c(1)*s3278*u3728
     & _56(i5)%b(2)
      l3_56728(i5,i7,i2)%c(1)=l3_56728(i5,i7,i2)%c(1)+l3_728(i7,
     & i2)%a(1)*u3728_56(i5)%c(1)+l3_728(i7,i2)%c(1)*s3278*u3728
     & _56(i5)%d(2)
      l3_56728(i5,i7,i2)%c(2)=l3_56728(i5,i7,i2)%c(2)+l3_728(i7,
     & i2)%c(2)*s3278*u3728_56(i5)%d(1)+l3_728(i7,i2)%a(2)*u3728
     & _56(i5)%c(2)
      l3_56728(i5,i7,i2)%a(2)=l3_56728(i5,i7,i2)%a(2)+l3_728(i7,
     & i2)%c(2)*s3278*u3728_56(i5)%b(1)+l3_728(i7,i2)%a(2)*u3728
     & _56(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id3).ne.1) then
* I                                                                     
* quqd -- p=p32,q=p3256
      quqd=p32(0)*p3256(0)-p32(1)*p3256(1)-p32(2)*p3256(2)-p32(3
     & )*p3256(3)
      ccr=zcr(id3)/(f3256)
      ccl=zcl(id3)/(f3256)
      do i5=1,2
* T0 -- qu=p32,qd=p3256,v=cz56(i5)%e,a=u32_56(i5)%a,b=u32_56(i5)%b,c=u32
* _56(i5)%c,d=u32_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p32(2)*p3256(3)-p3256(2)*p32(3))+p32
     & k0*(cz56(i5)%e(2)*p3256(3)-p3256(2)*cz56(i5)%e(3))-p3256k
     & 0*(cz56(i5)%e(2)*p32(3)-p32(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p32k0+p32(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p3256k0+p3256(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p32(0)-cz56(i5)%e(1)*p32(1)-cz56(i5)%e(
     & 2)*p32(2)-cz56(i5)%e(3)*p32(3)
      cvqd=cz56(i5)%e(0)*p3256(0)-cz56(i5)%e(1)*p3256(1)-cz56(i5
     & )%e(2)*p3256(2)-cz56(i5)%e(3)*p3256(3)
      cauxa=-cz56(i5)%ek0*quqd+p32k0*cvqd+p3256k0*cvqu
      cauxb=-cz56(i5)%ek0*p3256(2)+p3256k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p32(2)-p32k0*cz56(i5)%e(2)
      u32_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u32_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u32_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u32_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u32_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u32_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u32_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u32_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f3256)
      ccl=fcl(id3)/(f3256)
      do i5=1,2
* T0 -- qu=p32,qd=p3256,v=cf56(i5)%e,a=u32_56(i5)%a,b=u32_56(i5)%b,c=u32
* _56(i5)%c,d=u32_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p32(2)*p3256(3)-p3256(2)*p32(3))+p32
     & k0*(cf56(i5)%e(2)*p3256(3)-p3256(2)*cf56(i5)%e(3))-p3256k
     & 0*(cf56(i5)%e(2)*p32(3)-p32(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p32k0+p32(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p3256k0+p3256(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p32(0)-cf56(i5)%e(1)*p32(1)-cf56(i5)%e(
     & 2)*p32(2)-cf56(i5)%e(3)*p32(3)
      cvqd=cf56(i5)%e(0)*p3256(0)-cf56(i5)%e(1)*p3256(1)-cf56(i5
     & )%e(2)*p3256(2)-cf56(i5)%e(3)*p3256(3)
      cauxa=-cf56(i5)%ek0*quqd+p32k0*cvqd+p3256k0*cvqu
      cauxb=-cf56(i5)%ek0*p3256(2)+p3256k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p32(2)-p32k0*cf56(i5)%e(2)
      u32_56(i5)%a(1)=u32_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u32_56(i5)%a(2)=u32_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u32_56(i5)%b(1)=u32_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u32_56(i5)%b(2)=u32_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u32_56(i5)%c(1)=u32_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u32_56(i5)%c(2)=u32_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u32_56(i5)%d(1)=u32_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u32_56(i5)%d(2)=u32_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l3_256(i1,i5)%a,cc=l3_256(i1,i5)%c,a1=l3_2(i1)%a,c1=l3_2(i1
* )%c,a2=u32_56(i5)%a,b2=u32_56(i5)%b,c2=u32_56(i5)%c,d2=u32_56(i5)%d,pr
* q=s32,nsum=0
      l3_256(i1,i5)%a(1)=l3_2(i1)%a(1)*u32_56(i5)%a(1)+l3_2(i1)%
     & c(1)*s32*u32_56(i5)%b(2)
      l3_256(i1,i5)%c(1)=l3_2(i1)%a(1)*u32_56(i5)%c(1)+l3_2(i1)%
     & c(1)*s32*u32_56(i5)%d(2)
      l3_256(i1,i5)%c(2)=l3_2(i1)%c(2)*s32*u32_56(i5)%d(1)+l3_2(
     & i1)%a(2)*u32_56(i5)%c(2)
      l3_256(i1,i5)%a(2)=l3_2(i1)%c(2)*s32*u32_56(i5)%b(1)+l3_2(
     & i1)%a(2)*u32_56(i5)%a(2)
      end do
      end do
  
* quqd -- p=p356,q=p3256
      quqd=p356(0)*p3256(0)-p356(1)*p3256(1)-p356(2)*p3256(2)-p3
     & 56(3)*p3256(3)
      ccr=1.d0/(f3256)
      ccl=1.d0/(f3256)
      do i1=1,2
* T0 -- qu=p356,qd=p3256,v=ce2(i1)%e,a=u356_2(i1)%a,b=u356_2(i1)%b,c=u35
* 6_2(i1)%c,d=u356_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p356(2)*p3256(3)-p3256(2)*p356(3))+p3
     & 56k0*(ce2(i1)%e(2)*p3256(3)-p3256(2)*ce2(i1)%e(3))-p3256k
     & 0*(ce2(i1)%e(2)*p356(3)-p356(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p356k0+p356(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p3256k0+p3256(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p356(0)-ce2(i1)%e(1)*p356(1)-ce2(i1)%e(2
     & )*p356(2)-ce2(i1)%e(3)*p356(3)
      cvqd=ce2(i1)%e(0)*p3256(0)-ce2(i1)%e(1)*p3256(1)-ce2(i1)%e
     & (2)*p3256(2)-ce2(i1)%e(3)*p3256(3)
      cauxa=-ce2(i1)%ek0*quqd+p356k0*cvqd+p3256k0*cvqu
      cauxb=-ce2(i1)%ek0*p3256(2)+p3256k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p356(2)-p356k0*ce2(i1)%e(2)
      u356_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u356_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u356_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u356_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u356_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u356_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u356_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u356_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l3_256(i1,i5)%a,cc=l3_256(i1,i5)%c,a1=l3_56(i5)%a,c1=l3_56(
* i5)%c,a2=u356_2(i1)%a,b2=u356_2(i1)%b,c2=u356_2(i1)%c,d2=u356_2(i1)%d,
* prq=s356,nsum=1
      l3_256(i1,i5)%a(1)=l3_256(i1,i5)%a(1)+l3_56(i5)%a(1)*u356_
     & 2(i1)%a(1)+l3_56(i5)%c(1)*s356*u356_2(i1)%b(2)
      l3_256(i1,i5)%c(1)=l3_256(i1,i5)%c(1)+l3_56(i5)%a(1)*u356_
     & 2(i1)%c(1)+l3_56(i5)%c(1)*s356*u356_2(i1)%d(2)
      l3_256(i1,i5)%c(2)=l3_256(i1,i5)%c(2)+l3_56(i5)%c(2)*s356*
     & u356_2(i1)%d(1)+l3_56(i5)%a(2)*u356_2(i1)%c(2)
      l3_256(i1,i5)%a(2)=l3_256(i1,i5)%a(2)+l3_56(i5)%c(2)*s356*
     & u356_2(i1)%b(1)+l3_56(i5)%a(2)*u356_2(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id7).ne.1) then
* quqd -- p=p32,q=p456
      quqd=p32(0)*p456(0)-p32(1)*p456(1)-p32(2)*p456(2)-p32(3)*p
     & 456(3)
      ccr=zcr(id3)/(f456)
      ccl=zcl(id3)/(f456)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p32,qd=p456,v=cz718(i7,i2)%e,a=u32_718(i7,i2)%a,b=u32_718(i7,
* i2)%b,c=u32_718(i7,i2)%c,d=u32_718(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz718(i7,i2)%ek0*(p32(2)*p456(3)-p456(2)*p32(3))+p
     & 32k0*(cz718(i7,i2)%e(2)*p456(3)-p456(2)*cz718(i7,i2)%e(3)
     & )-p456k0*(cz718(i7,i2)%e(2)*p32(3)-p32(2)*cz718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz718(i7,i2)%e(3)*p32k0+p32(3)*cz718(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz718(i7,i2)%e(3)*p456k0+p456(3)*cz718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz718(i7,i2)%e(0)*p32(0)-cz718(i7,i2)%e(1)*p32(1)-cz7
     & 18(i7,i2)%e(2)*p32(2)-cz718(i7,i2)%e(3)*p32(3)
      cvqd=cz718(i7,i2)%e(0)*p456(0)-cz718(i7,i2)%e(1)*p456(1)-c
     & z718(i7,i2)%e(2)*p456(2)-cz718(i7,i2)%e(3)*p456(3)
      cauxa=-cz718(i7,i2)%ek0*quqd+p32k0*cvqd+p456k0*cvqu
      cauxb=-cz718(i7,i2)%ek0*p456(2)+p456k0*cz718(i7,i2)%e(2)
      cauxc=+cz718(i7,i2)%ek0*p32(2)-p32k0*cz718(i7,i2)%e(2)
      u32_718(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u32_718(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u32_718(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u32_718(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u32_718(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u32_718(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u32_718(i7,i2)%d(1)=ccl*cz718(i7,i2)%ek0
      u32_718(i7,i2)%d(2)=ccr*cz718(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id3)/(f456)
      ccl=fcl(id3)/(f456)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p32,qd=p456,v=cf718(i7,i2)%e,a=u32_718(i7,i2)%a,b=u32_718(i7,
* i2)%b,c=u32_718(i7,i2)%c,d=u32_718(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf718(i7,i2)%ek0*(p32(2)*p456(3)-p456(2)*p32(3))+p
     & 32k0*(cf718(i7,i2)%e(2)*p456(3)-p456(2)*cf718(i7,i2)%e(3)
     & )-p456k0*(cf718(i7,i2)%e(2)*p32(3)-p32(2)*cf718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf718(i7,i2)%e(3)*p32k0+p32(3)*cf718(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf718(i7,i2)%e(3)*p456k0+p456(3)*cf718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf718(i7,i2)%e(0)*p32(0)-cf718(i7,i2)%e(1)*p32(1)-cf7
     & 18(i7,i2)%e(2)*p32(2)-cf718(i7,i2)%e(3)*p32(3)
      cvqd=cf718(i7,i2)%e(0)*p456(0)-cf718(i7,i2)%e(1)*p456(1)-c
     & f718(i7,i2)%e(2)*p456(2)-cf718(i7,i2)%e(3)*p456(3)
      cauxa=-cf718(i7,i2)%ek0*quqd+p32k0*cvqd+p456k0*cvqu
      cauxb=-cf718(i7,i2)%ek0*p456(2)+p456k0*cf718(i7,i2)%e(2)
      cauxc=+cf718(i7,i2)%ek0*p32(2)-p32k0*cf718(i7,i2)%e(2)
      u32_718(i7,i2)%a(1)=u32_718(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u32_718(i7,i2)%a(2)=u32_718(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u32_718(i7,i2)%b(1)=u32_718(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u32_718(i7,i2)%b(2)=u32_718(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u32_718(i7,i2)%c(1)=u32_718(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u32_718(i7,i2)%c(2)=u32_718(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u32_718(i7,i2)%d(1)=u32_718(i7,i2)%d(1)+ccl*cf718(i7,i2)%e
     & k0
      u32_718(i7,i2)%d(2)=u32_718(i7,i2)%d(2)+ccr*cf718(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_2718(i1,i7,i2)%a,cc=l3_2718(i1,i7,i2)%c,a1=l3_2(i1)%a,c1
* =l3_2(i1)%c,a2=u32_718(i7,i2)%a,b2=u32_718(i7,i2)%b,c2=u32_718(i7,i2)%
* c,d2=u32_718(i7,i2)%d,prq=s32,nsum=0
      l3_2718(i1,i7,i2)%a(1)=l3_2(i1)%a(1)*u32_718(i7,i2)%a(1)+l
     & 3_2(i1)%c(1)*s32*u32_718(i7,i2)%b(2)
      l3_2718(i1,i7,i2)%c(1)=l3_2(i1)%a(1)*u32_718(i7,i2)%c(1)+l
     & 3_2(i1)%c(1)*s32*u32_718(i7,i2)%d(2)
      l3_2718(i1,i7,i2)%c(2)=l3_2(i1)%c(2)*s32*u32_718(i7,i2)%d(
     & 1)+l3_2(i1)%a(2)*u32_718(i7,i2)%c(2)
      l3_2718(i1,i7,i2)%a(2)=l3_2(i1)%c(2)*s32*u32_718(i7,i2)%b(
     & 1)+l3_2(i1)%a(2)*u32_718(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u3178_2                                                
* quqd -- p=p3178,q=p456
      quqd=p3178(0)*p456(0)-p3178(1)*p456(1)-p3178(2)*p456(2)-p3
     & 178(3)*p456(3)
      ccr=1.d0/(f456)
      ccl=1.d0/(f456)
      do i1=1,2
* T0 -- qu=p3178,qd=p456,v=ce2(i1)%e,a=u3718_2(i1)%a,b=u3718_2(i1)%b,c=u
* 3718_2(i1)%c,d=u3718_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p3178(2)*p456(3)-p456(2)*p3178(3))+p3
     & 178k0*(ce2(i1)%e(2)*p456(3)-p456(2)*ce2(i1)%e(3))-p456k0*
     & (ce2(i1)%e(2)*p3178(3)-p3178(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p3178k0+p3178(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p456k0+p456(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p3178(0)-ce2(i1)%e(1)*p3178(1)-ce2(i1)%e
     & (2)*p3178(2)-ce2(i1)%e(3)*p3178(3)
      cvqd=ce2(i1)%e(0)*p456(0)-ce2(i1)%e(1)*p456(1)-ce2(i1)%e(2
     & )*p456(2)-ce2(i1)%e(3)*p456(3)
      cauxa=-ce2(i1)%ek0*quqd+p3178k0*cvqd+p456k0*cvqu
      cauxb=-ce2(i1)%ek0*p456(2)+p456k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p3178(2)-p3178k0*ce2(i1)%e(2)
      u3718_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u3718_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u3718_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u3718_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u3718_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u3718_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u3718_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u3718_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_2718(i1,i7,i2)%a,cc=l3_2718(i1,i7,i2)%c,a1=l3_718(i7,i2)
* %a,c1=l3_718(i7,i2)%c,a2=u3718_2(i1)%a,b2=u3718_2(i1)%b,c2=u3718_2(i1)
* %c,d2=u3718_2(i1)%d,prq=s3178,nsum=1
      l3_2718(i1,i7,i2)%a(1)=l3_2718(i1,i7,i2)%a(1)+l3_718(i7,i2
     & )%a(1)*u3718_2(i1)%a(1)+l3_718(i7,i2)%c(1)*s3178*u3718_2(
     & i1)%b(2)
      l3_2718(i1,i7,i2)%c(1)=l3_2718(i1,i7,i2)%c(1)+l3_718(i7,i2
     & )%a(1)*u3718_2(i1)%c(1)+l3_718(i7,i2)%c(1)*s3178*u3718_2(
     & i1)%d(2)
      l3_2718(i1,i7,i2)%c(2)=l3_2718(i1,i7,i2)%c(2)+l3_718(i7,i2
     & )%c(2)*s3178*u3718_2(i1)%d(1)+l3_718(i7,i2)%a(2)*u3718_2(
     & i1)%c(2)
      l3_2718(i1,i7,i2)%a(2)=l3_2718(i1,i7,i2)%a(2)+l3_718(i7,i2
     & )%c(2)*s3178*u3718_2(i1)%b(1)+l3_718(i7,i2)%a(2)*u3718_2(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p356,q=p42
      quqd=p356(0)*p42(0)-p356(1)*p42(1)-p356(2)*p42(2)-p356(3)*
     & p42(3)
      ccr=zcr(id3)/(f42)
      ccl=zcl(id3)/(f42)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p356,qd=p42,v=cz718(i7,i2)%e,a=u356_718(i7,i2)%a,b=u356_718(i
* 7,i2)%b,c=u356_718(i7,i2)%c,d=u356_718(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz718(i7,i2)%ek0*(p356(2)*p42(3)-p42(2)*p356(3))+p
     & 356k0*(cz718(i7,i2)%e(2)*p42(3)-p42(2)*cz718(i7,i2)%e(3))
     & -p42k0*(cz718(i7,i2)%e(2)*p356(3)-p356(2)*cz718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz718(i7,i2)%e(3)*p356k0+p356(3)*cz718(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz718(i7,i2)%e(3)*p42k0+p42(3)*cz718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz718(i7,i2)%e(0)*p356(0)-cz718(i7,i2)%e(1)*p356(1)-c
     & z718(i7,i2)%e(2)*p356(2)-cz718(i7,i2)%e(3)*p356(3)
      cvqd=cz718(i7,i2)%e(0)*p42(0)-cz718(i7,i2)%e(1)*p42(1)-cz7
     & 18(i7,i2)%e(2)*p42(2)-cz718(i7,i2)%e(3)*p42(3)
      cauxa=-cz718(i7,i2)%ek0*quqd+p356k0*cvqd+p42k0*cvqu
      cauxb=-cz718(i7,i2)%ek0*p42(2)+p42k0*cz718(i7,i2)%e(2)
      cauxc=+cz718(i7,i2)%ek0*p356(2)-p356k0*cz718(i7,i2)%e(2)
      u356_718(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u356_718(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u356_718(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u356_718(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u356_718(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u356_718(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u356_718(i7,i2)%d(1)=ccl*cz718(i7,i2)%ek0
      u356_718(i7,i2)%d(2)=ccr*cz718(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id3)/(f42)
      ccl=fcl(id3)/(f42)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p356,qd=p42,v=cf718(i7,i2)%e,a=u356_718(i7,i2)%a,b=u356_718(i
* 7,i2)%b,c=u356_718(i7,i2)%c,d=u356_718(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf718(i7,i2)%ek0*(p356(2)*p42(3)-p42(2)*p356(3))+p
     & 356k0*(cf718(i7,i2)%e(2)*p42(3)-p42(2)*cf718(i7,i2)%e(3))
     & -p42k0*(cf718(i7,i2)%e(2)*p356(3)-p356(2)*cf718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf718(i7,i2)%e(3)*p356k0+p356(3)*cf718(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf718(i7,i2)%e(3)*p42k0+p42(3)*cf718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf718(i7,i2)%e(0)*p356(0)-cf718(i7,i2)%e(1)*p356(1)-c
     & f718(i7,i2)%e(2)*p356(2)-cf718(i7,i2)%e(3)*p356(3)
      cvqd=cf718(i7,i2)%e(0)*p42(0)-cf718(i7,i2)%e(1)*p42(1)-cf7
     & 18(i7,i2)%e(2)*p42(2)-cf718(i7,i2)%e(3)*p42(3)
      cauxa=-cf718(i7,i2)%ek0*quqd+p356k0*cvqd+p42k0*cvqu
      cauxb=-cf718(i7,i2)%ek0*p42(2)+p42k0*cf718(i7,i2)%e(2)
      cauxc=+cf718(i7,i2)%ek0*p356(2)-p356k0*cf718(i7,i2)%e(2)
      u356_718(i7,i2)%a(1)=u356_718(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u356_718(i7,i2)%a(2)=u356_718(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u356_718(i7,i2)%b(1)=u356_718(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u356_718(i7,i2)%b(2)=u356_718(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u356_718(i7,i2)%c(1)=u356_718(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u356_718(i7,i2)%c(2)=u356_718(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u356_718(i7,i2)%d(1)=u356_718(i7,i2)%d(1)+ccl*cf718(i7,i2)
     & %ek0
      u356_718(i7,i2)%d(2)=u356_718(i7,i2)%d(2)+ccr*cf718(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_56718(i5,i7,i2)%a,cc=l3_56718(i5,i7,i2)%c,a1=l3_56(i5)%a
* ,c1=l3_56(i5)%c,a2=u356_718(i7,i2)%a,b2=u356_718(i7,i2)%b,c2=u356_718(
* i7,i2)%c,d2=u356_718(i7,i2)%d,prq=s356,nsum=0
      l3_56718(i5,i7,i2)%a(1)=l3_56(i5)%a(1)*u356_718(i7,i2)%a(1
     & )+l3_56(i5)%c(1)*s356*u356_718(i7,i2)%b(2)
      l3_56718(i5,i7,i2)%c(1)=l3_56(i5)%a(1)*u356_718(i7,i2)%c(1
     & )+l3_56(i5)%c(1)*s356*u356_718(i7,i2)%d(2)
      l3_56718(i5,i7,i2)%c(2)=l3_56(i5)%c(2)*s356*u356_718(i7,i2
     & )%d(1)+l3_56(i5)%a(2)*u356_718(i7,i2)%c(2)
      l3_56718(i5,i7,i2)%a(2)=l3_56(i5)%c(2)*s356*u356_718(i7,i2
     & )%b(1)+l3_56(i5)%a(2)*u356_718(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u3178_56                                               
* quqd -- p=p3178,q=p42
      quqd=p3178(0)*p42(0)-p3178(1)*p42(1)-p3178(2)*p42(2)-p3178
     & (3)*p42(3)
      ccr=zcr(id3)/(f42)
      ccl=zcl(id3)/(f42)
      do i5=1,2
* T0 -- qu=p3178,qd=p42,v=cz56(i5)%e,a=u3718_56(i5)%a,b=u3718_56(i5)%b,c
* =u3718_56(i5)%c,d=u3718_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p3178(2)*p42(3)-p42(2)*p3178(3))+p31
     & 78k0*(cz56(i5)%e(2)*p42(3)-p42(2)*cz56(i5)%e(3))-p42k0*(c
     & z56(i5)%e(2)*p3178(3)-p3178(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p3178k0+p3178(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p42k0+p42(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p3178(0)-cz56(i5)%e(1)*p3178(1)-cz56(i5
     & )%e(2)*p3178(2)-cz56(i5)%e(3)*p3178(3)
      cvqd=cz56(i5)%e(0)*p42(0)-cz56(i5)%e(1)*p42(1)-cz56(i5)%e(
     & 2)*p42(2)-cz56(i5)%e(3)*p42(3)
      cauxa=-cz56(i5)%ek0*quqd+p3178k0*cvqd+p42k0*cvqu
      cauxb=-cz56(i5)%ek0*p42(2)+p42k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p3178(2)-p3178k0*cz56(i5)%e(2)
      u3718_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u3718_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u3718_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u3718_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u3718_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u3718_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u3718_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u3718_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f42)
      ccl=fcl(id3)/(f42)
      do i5=1,2
* T0 -- qu=p3178,qd=p42,v=cf56(i5)%e,a=u3718_56(i5)%a,b=u3718_56(i5)%b,c
* =u3718_56(i5)%c,d=u3718_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p3178(2)*p42(3)-p42(2)*p3178(3))+p31
     & 78k0*(cf56(i5)%e(2)*p42(3)-p42(2)*cf56(i5)%e(3))-p42k0*(c
     & f56(i5)%e(2)*p3178(3)-p3178(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p3178k0+p3178(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p42k0+p42(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p3178(0)-cf56(i5)%e(1)*p3178(1)-cf56(i5
     & )%e(2)*p3178(2)-cf56(i5)%e(3)*p3178(3)
      cvqd=cf56(i5)%e(0)*p42(0)-cf56(i5)%e(1)*p42(1)-cf56(i5)%e(
     & 2)*p42(2)-cf56(i5)%e(3)*p42(3)
      cauxa=-cf56(i5)%ek0*quqd+p3178k0*cvqd+p42k0*cvqu
      cauxb=-cf56(i5)%ek0*p42(2)+p42k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p3178(2)-p3178k0*cf56(i5)%e(2)
      u3718_56(i5)%a(1)=u3718_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u3718_56(i5)%a(2)=u3718_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u3718_56(i5)%b(1)=u3718_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u3718_56(i5)%b(2)=u3718_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u3718_56(i5)%c(1)=u3718_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u3718_56(i5)%c(2)=u3718_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u3718_56(i5)%d(1)=u3718_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u3718_56(i5)%d(2)=u3718_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_56718(i5,i7,i2)%a,cc=l3_56718(i5,i7,i2)%c,a1=l3_718(i7,i
* 2)%a,c1=l3_718(i7,i2)%c,a2=u3718_56(i5)%a,b2=u3718_56(i5)%b,c2=u3718_5
* 6(i5)%c,d2=u3718_56(i5)%d,prq=s3178,nsum=1
      l3_56718(i5,i7,i2)%a(1)=l3_56718(i5,i7,i2)%a(1)+l3_718(i7,
     & i2)%a(1)*u3718_56(i5)%a(1)+l3_718(i7,i2)%c(1)*s3178*u3718
     & _56(i5)%b(2)
      l3_56718(i5,i7,i2)%c(1)=l3_56718(i5,i7,i2)%c(1)+l3_718(i7,
     & i2)%a(1)*u3718_56(i5)%c(1)+l3_718(i7,i2)%c(1)*s3178*u3718
     & _56(i5)%d(2)
      l3_56718(i5,i7,i2)%c(2)=l3_56718(i5,i7,i2)%c(2)+l3_718(i7,
     & i2)%c(2)*s3178*u3718_56(i5)%d(1)+l3_718(i7,i2)%a(2)*u3718
     & _56(i5)%c(2)
      l3_56718(i5,i7,i2)%a(2)=l3_56718(i5,i7,i2)%a(2)+l3_718(i7,
     & i2)%c(2)*s3178*u3718_56(i5)%b(1)+l3_718(i7,i2)%a(2)*u3718
     & _56(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id3).ne.1) then
* I                                                                     
* quqd -- p=p31,q=p3178
      quqd=p31(0)*p3178(0)-p31(1)*p3178(1)-p31(2)*p3178(2)-p31(3
     & )*p3178(3)
      ccr=zcr(id3)/(f3178)
      ccl=zcl(id3)/(f3178)
      do i5=1,2
* T0 -- qu=p31,qd=p3178,v=cz78(i5)%e,a=u31_78(i5)%a,b=u31_78(i5)%b,c=u31
* _78(i5)%c,d=u31_78(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p31(2)*p3178(3)-p3178(2)*p31(3))+p31
     & k0*(cz78(i5)%e(2)*p3178(3)-p3178(2)*cz78(i5)%e(3))-p3178k
     & 0*(cz78(i5)%e(2)*p31(3)-p31(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p31k0+p31(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i5)%e(3)*p3178k0+p3178(3)*cz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i5)%e(0)*p31(0)-cz78(i5)%e(1)*p31(1)-cz78(i5)%e(
     & 2)*p31(2)-cz78(i5)%e(3)*p31(3)
      cvqd=cz78(i5)%e(0)*p3178(0)-cz78(i5)%e(1)*p3178(1)-cz78(i5
     & )%e(2)*p3178(2)-cz78(i5)%e(3)*p3178(3)
      cauxa=-cz78(i5)%ek0*quqd+p31k0*cvqd+p3178k0*cvqu
      cauxb=-cz78(i5)%ek0*p3178(2)+p3178k0*cz78(i5)%e(2)
      cauxc=+cz78(i5)%ek0*p31(2)-p31k0*cz78(i5)%e(2)
      u31_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      u31_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      u31_78(i5)%b(1)=ccl*(cauxb-ceps_2)
      u31_78(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u31_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      u31_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u31_78(i5)%d(1)=ccl*cz78(i5)%ek0
      u31_78(i5)%d(2)=ccr*cz78(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f3178)
      ccl=fcl(id3)/(f3178)
      do i5=1,2
* T0 -- qu=p31,qd=p3178,v=cf78(i5)%e,a=u31_78(i5)%a,b=u31_78(i5)%b,c=u31
* _78(i5)%c,d=u31_78(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p31(2)*p3178(3)-p3178(2)*p31(3))+p31
     & k0*(cf78(i5)%e(2)*p3178(3)-p3178(2)*cf78(i5)%e(3))-p3178k
     & 0*(cf78(i5)%e(2)*p31(3)-p31(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p31k0+p31(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i5)%e(3)*p3178k0+p3178(3)*cf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i5)%e(0)*p31(0)-cf78(i5)%e(1)*p31(1)-cf78(i5)%e(
     & 2)*p31(2)-cf78(i5)%e(3)*p31(3)
      cvqd=cf78(i5)%e(0)*p3178(0)-cf78(i5)%e(1)*p3178(1)-cf78(i5
     & )%e(2)*p3178(2)-cf78(i5)%e(3)*p3178(3)
      cauxa=-cf78(i5)%ek0*quqd+p31k0*cvqd+p3178k0*cvqu
      cauxb=-cf78(i5)%ek0*p3178(2)+p3178k0*cf78(i5)%e(2)
      cauxc=+cf78(i5)%ek0*p31(2)-p31k0*cf78(i5)%e(2)
      u31_78(i5)%a(1)=u31_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      u31_78(i5)%a(2)=u31_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      u31_78(i5)%b(1)=u31_78(i5)%b(1)+ccl*(cauxb-ceps_2)
      u31_78(i5)%b(2)=u31_78(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u31_78(i5)%c(1)=u31_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      u31_78(i5)%c(2)=u31_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u31_78(i5)%d(1)=u31_78(i5)%d(1)+ccl*cf78(i5)%ek0
      u31_78(i5)%d(2)=u31_78(i5)%d(2)+ccr*cf78(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l3_178(i1,i5)%a,cc=l3_178(i1,i5)%c,a1=l3_1(i1)%a,c1=l3_1(i1
* )%c,a2=u31_78(i5)%a,b2=u31_78(i5)%b,c2=u31_78(i5)%c,d2=u31_78(i5)%d,pr
* q=s31,nsum=0
      l3_178(i1,i5)%a(1)=l3_1(i1)%a(1)*u31_78(i5)%a(1)+l3_1(i1)%
     & c(1)*s31*u31_78(i5)%b(2)
      l3_178(i1,i5)%c(1)=l3_1(i1)%a(1)*u31_78(i5)%c(1)+l3_1(i1)%
     & c(1)*s31*u31_78(i5)%d(2)
      l3_178(i1,i5)%c(2)=l3_1(i1)%c(2)*s31*u31_78(i5)%d(1)+l3_1(
     & i1)%a(2)*u31_78(i5)%c(2)
      l3_178(i1,i5)%a(2)=l3_1(i1)%c(2)*s31*u31_78(i5)%b(1)+l3_1(
     & i1)%a(2)*u31_78(i5)%a(2)
      end do
      end do
  
* quqd -- p=p378,q=p3178
      quqd=p378(0)*p3178(0)-p378(1)*p3178(1)-p378(2)*p3178(2)-p3
     & 78(3)*p3178(3)
      ccr=1.d0/(f3178)
      ccl=1.d0/(f3178)
      do i1=1,2
* T0 -- qu=p378,qd=p3178,v=ce1(i1)%e,a=u378_1(i1)%a,b=u378_1(i1)%b,c=u37
* 8_1(i1)%c,d=u378_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p378(2)*p3178(3)-p3178(2)*p378(3))+p3
     & 78k0*(ce1(i1)%e(2)*p3178(3)-p3178(2)*ce1(i1)%e(3))-p3178k
     & 0*(ce1(i1)%e(2)*p378(3)-p378(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p378k0+p378(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p3178k0+p3178(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p378(0)-ce1(i1)%e(1)*p378(1)-ce1(i1)%e(2
     & )*p378(2)-ce1(i1)%e(3)*p378(3)
      cvqd=ce1(i1)%e(0)*p3178(0)-ce1(i1)%e(1)*p3178(1)-ce1(i1)%e
     & (2)*p3178(2)-ce1(i1)%e(3)*p3178(3)
      cauxa=-ce1(i1)%ek0*quqd+p378k0*cvqd+p3178k0*cvqu
      cauxb=-ce1(i1)%ek0*p3178(2)+p3178k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p378(2)-p378k0*ce1(i1)%e(2)
      u378_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u378_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u378_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u378_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u378_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u378_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u378_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u378_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l3_178(i1,i5)%a,cc=l3_178(i1,i5)%c,a1=l3_78(i5)%a,c1=l3_78(
* i5)%c,a2=u378_1(i1)%a,b2=u378_1(i1)%b,c2=u378_1(i1)%c,d2=u378_1(i1)%d,
* prq=s378,nsum=1
      l3_178(i1,i5)%a(1)=l3_178(i1,i5)%a(1)+l3_78(i5)%a(1)*u378_
     & 1(i1)%a(1)+l3_78(i5)%c(1)*s378*u378_1(i1)%b(2)
      l3_178(i1,i5)%c(1)=l3_178(i1,i5)%c(1)+l3_78(i5)%a(1)*u378_
     & 1(i1)%c(1)+l3_78(i5)%c(1)*s378*u378_1(i1)%d(2)
      l3_178(i1,i5)%c(2)=l3_178(i1,i5)%c(2)+l3_78(i5)%c(2)*s378*
     & u378_1(i1)%d(1)+l3_78(i5)%a(2)*u378_1(i1)%c(2)
      l3_178(i1,i5)%a(2)=l3_178(i1,i5)%a(2)+l3_78(i5)%c(2)*s378*
     & u378_1(i1)%b(1)+l3_78(i5)%a(2)*u378_1(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id5).ne.1) then
* quqd -- p=p31,q=p478
      quqd=p31(0)*p478(0)-p31(1)*p478(1)-p31(2)*p478(2)-p31(3)*p
     & 478(3)
      ccr=zcr(id3)/(f478)
      ccl=zcl(id3)/(f478)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p31,qd=p478,v=cz526(i7,i2)%e,a=u31_526(i7,i2)%a,b=u31_526(i7,
* i2)%b,c=u31_526(i7,i2)%c,d=u31_526(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz526(i7,i2)%ek0*(p31(2)*p478(3)-p478(2)*p31(3))+p
     & 31k0*(cz526(i7,i2)%e(2)*p478(3)-p478(2)*cz526(i7,i2)%e(3)
     & )-p478k0*(cz526(i7,i2)%e(2)*p31(3)-p31(2)*cz526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz526(i7,i2)%e(3)*p31k0+p31(3)*cz526(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz526(i7,i2)%e(3)*p478k0+p478(3)*cz526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz526(i7,i2)%e(0)*p31(0)-cz526(i7,i2)%e(1)*p31(1)-cz5
     & 26(i7,i2)%e(2)*p31(2)-cz526(i7,i2)%e(3)*p31(3)
      cvqd=cz526(i7,i2)%e(0)*p478(0)-cz526(i7,i2)%e(1)*p478(1)-c
     & z526(i7,i2)%e(2)*p478(2)-cz526(i7,i2)%e(3)*p478(3)
      cauxa=-cz526(i7,i2)%ek0*quqd+p31k0*cvqd+p478k0*cvqu
      cauxb=-cz526(i7,i2)%ek0*p478(2)+p478k0*cz526(i7,i2)%e(2)
      cauxc=+cz526(i7,i2)%ek0*p31(2)-p31k0*cz526(i7,i2)%e(2)
      u31_526(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u31_526(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u31_526(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u31_526(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u31_526(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u31_526(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u31_526(i7,i2)%d(1)=ccl*cz526(i7,i2)%ek0
      u31_526(i7,i2)%d(2)=ccr*cz526(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id3)/(f478)
      ccl=fcl(id3)/(f478)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p31,qd=p478,v=cf526(i7,i2)%e,a=u31_526(i7,i2)%a,b=u31_526(i7,
* i2)%b,c=u31_526(i7,i2)%c,d=u31_526(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf526(i7,i2)%ek0*(p31(2)*p478(3)-p478(2)*p31(3))+p
     & 31k0*(cf526(i7,i2)%e(2)*p478(3)-p478(2)*cf526(i7,i2)%e(3)
     & )-p478k0*(cf526(i7,i2)%e(2)*p31(3)-p31(2)*cf526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf526(i7,i2)%e(3)*p31k0+p31(3)*cf526(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf526(i7,i2)%e(3)*p478k0+p478(3)*cf526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf526(i7,i2)%e(0)*p31(0)-cf526(i7,i2)%e(1)*p31(1)-cf5
     & 26(i7,i2)%e(2)*p31(2)-cf526(i7,i2)%e(3)*p31(3)
      cvqd=cf526(i7,i2)%e(0)*p478(0)-cf526(i7,i2)%e(1)*p478(1)-c
     & f526(i7,i2)%e(2)*p478(2)-cf526(i7,i2)%e(3)*p478(3)
      cauxa=-cf526(i7,i2)%ek0*quqd+p31k0*cvqd+p478k0*cvqu
      cauxb=-cf526(i7,i2)%ek0*p478(2)+p478k0*cf526(i7,i2)%e(2)
      cauxc=+cf526(i7,i2)%ek0*p31(2)-p31k0*cf526(i7,i2)%e(2)
      u31_526(i7,i2)%a(1)=u31_526(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u31_526(i7,i2)%a(2)=u31_526(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u31_526(i7,i2)%b(1)=u31_526(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u31_526(i7,i2)%b(2)=u31_526(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u31_526(i7,i2)%c(1)=u31_526(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u31_526(i7,i2)%c(2)=u31_526(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u31_526(i7,i2)%d(1)=u31_526(i7,i2)%d(1)+ccl*cf526(i7,i2)%e
     & k0
      u31_526(i7,i2)%d(2)=u31_526(i7,i2)%d(2)+ccr*cf526(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_1526(i1,i7,i2)%a,cc=l3_1526(i1,i7,i2)%c,a1=l3_1(i1)%a,c1
* =l3_1(i1)%c,a2=u31_526(i7,i2)%a,b2=u31_526(i7,i2)%b,c2=u31_526(i7,i2)%
* c,d2=u31_526(i7,i2)%d,prq=s31,nsum=0
      l3_1526(i1,i7,i2)%a(1)=l3_1(i1)%a(1)*u31_526(i7,i2)%a(1)+l
     & 3_1(i1)%c(1)*s31*u31_526(i7,i2)%b(2)
      l3_1526(i1,i7,i2)%c(1)=l3_1(i1)%a(1)*u31_526(i7,i2)%c(1)+l
     & 3_1(i1)%c(1)*s31*u31_526(i7,i2)%d(2)
      l3_1526(i1,i7,i2)%c(2)=l3_1(i1)%c(2)*s31*u31_526(i7,i2)%d(
     & 1)+l3_1(i1)%a(2)*u31_526(i7,i2)%c(2)
      l3_1526(i1,i7,i2)%a(2)=l3_1(i1)%c(2)*s31*u31_526(i7,i2)%b(
     & 1)+l3_1(i1)%a(2)*u31_526(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u3256_1                                                
* quqd -- p=p3256,q=p478
      quqd=p3256(0)*p478(0)-p3256(1)*p478(1)-p3256(2)*p478(2)-p3
     & 256(3)*p478(3)
      ccr=1.d0/(f478)
      ccl=1.d0/(f478)
      do i1=1,2
* T0 -- qu=p3256,qd=p478,v=ce1(i1)%e,a=u3526_1(i1)%a,b=u3526_1(i1)%b,c=u
* 3526_1(i1)%c,d=u3526_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p3256(2)*p478(3)-p478(2)*p3256(3))+p3
     & 256k0*(ce1(i1)%e(2)*p478(3)-p478(2)*ce1(i1)%e(3))-p478k0*
     & (ce1(i1)%e(2)*p3256(3)-p3256(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p3256k0+p3256(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p478k0+p478(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p3256(0)-ce1(i1)%e(1)*p3256(1)-ce1(i1)%e
     & (2)*p3256(2)-ce1(i1)%e(3)*p3256(3)
      cvqd=ce1(i1)%e(0)*p478(0)-ce1(i1)%e(1)*p478(1)-ce1(i1)%e(2
     & )*p478(2)-ce1(i1)%e(3)*p478(3)
      cauxa=-ce1(i1)%ek0*quqd+p3256k0*cvqd+p478k0*cvqu
      cauxb=-ce1(i1)%ek0*p478(2)+p478k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p3256(2)-p3256k0*ce1(i1)%e(2)
      u3526_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u3526_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u3526_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u3526_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u3526_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u3526_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u3526_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u3526_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_1526(i1,i7,i2)%a,cc=l3_1526(i1,i7,i2)%c,a1=l3_526(i7,i2)
* %a,c1=l3_526(i7,i2)%c,a2=u3526_1(i1)%a,b2=u3526_1(i1)%b,c2=u3526_1(i1)
* %c,d2=u3526_1(i1)%d,prq=s3256,nsum=1
      l3_1526(i1,i7,i2)%a(1)=l3_1526(i1,i7,i2)%a(1)+l3_526(i7,i2
     & )%a(1)*u3526_1(i1)%a(1)+l3_526(i7,i2)%c(1)*s3256*u3526_1(
     & i1)%b(2)
      l3_1526(i1,i7,i2)%c(1)=l3_1526(i1,i7,i2)%c(1)+l3_526(i7,i2
     & )%a(1)*u3526_1(i1)%c(1)+l3_526(i7,i2)%c(1)*s3256*u3526_1(
     & i1)%d(2)
      l3_1526(i1,i7,i2)%c(2)=l3_1526(i1,i7,i2)%c(2)+l3_526(i7,i2
     & )%c(2)*s3256*u3526_1(i1)%d(1)+l3_526(i7,i2)%a(2)*u3526_1(
     & i1)%c(2)
      l3_1526(i1,i7,i2)%a(2)=l3_1526(i1,i7,i2)%a(2)+l3_526(i7,i2
     & )%c(2)*s3256*u3526_1(i1)%b(1)+l3_526(i7,i2)%a(2)*u3526_1(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p378,q=p41
      quqd=p378(0)*p41(0)-p378(1)*p41(1)-p378(2)*p41(2)-p378(3)*
     & p41(3)
      ccr=zcr(id3)/(f41)
      ccl=zcl(id3)/(f41)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p378,qd=p41,v=cz526(i7,i2)%e,a=u378_526(i7,i2)%a,b=u378_526(i
* 7,i2)%b,c=u378_526(i7,i2)%c,d=u378_526(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz526(i7,i2)%ek0*(p378(2)*p41(3)-p41(2)*p378(3))+p
     & 378k0*(cz526(i7,i2)%e(2)*p41(3)-p41(2)*cz526(i7,i2)%e(3))
     & -p41k0*(cz526(i7,i2)%e(2)*p378(3)-p378(2)*cz526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz526(i7,i2)%e(3)*p378k0+p378(3)*cz526(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz526(i7,i2)%e(3)*p41k0+p41(3)*cz526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz526(i7,i2)%e(0)*p378(0)-cz526(i7,i2)%e(1)*p378(1)-c
     & z526(i7,i2)%e(2)*p378(2)-cz526(i7,i2)%e(3)*p378(3)
      cvqd=cz526(i7,i2)%e(0)*p41(0)-cz526(i7,i2)%e(1)*p41(1)-cz5
     & 26(i7,i2)%e(2)*p41(2)-cz526(i7,i2)%e(3)*p41(3)
      cauxa=-cz526(i7,i2)%ek0*quqd+p378k0*cvqd+p41k0*cvqu
      cauxb=-cz526(i7,i2)%ek0*p41(2)+p41k0*cz526(i7,i2)%e(2)
      cauxc=+cz526(i7,i2)%ek0*p378(2)-p378k0*cz526(i7,i2)%e(2)
      u378_526(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u378_526(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u378_526(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u378_526(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u378_526(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u378_526(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u378_526(i7,i2)%d(1)=ccl*cz526(i7,i2)%ek0
      u378_526(i7,i2)%d(2)=ccr*cz526(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id3)/(f41)
      ccl=fcl(id3)/(f41)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p378,qd=p41,v=cf526(i7,i2)%e,a=u378_526(i7,i2)%a,b=u378_526(i
* 7,i2)%b,c=u378_526(i7,i2)%c,d=u378_526(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf526(i7,i2)%ek0*(p378(2)*p41(3)-p41(2)*p378(3))+p
     & 378k0*(cf526(i7,i2)%e(2)*p41(3)-p41(2)*cf526(i7,i2)%e(3))
     & -p41k0*(cf526(i7,i2)%e(2)*p378(3)-p378(2)*cf526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf526(i7,i2)%e(3)*p378k0+p378(3)*cf526(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf526(i7,i2)%e(3)*p41k0+p41(3)*cf526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf526(i7,i2)%e(0)*p378(0)-cf526(i7,i2)%e(1)*p378(1)-c
     & f526(i7,i2)%e(2)*p378(2)-cf526(i7,i2)%e(3)*p378(3)
      cvqd=cf526(i7,i2)%e(0)*p41(0)-cf526(i7,i2)%e(1)*p41(1)-cf5
     & 26(i7,i2)%e(2)*p41(2)-cf526(i7,i2)%e(3)*p41(3)
      cauxa=-cf526(i7,i2)%ek0*quqd+p378k0*cvqd+p41k0*cvqu
      cauxb=-cf526(i7,i2)%ek0*p41(2)+p41k0*cf526(i7,i2)%e(2)
      cauxc=+cf526(i7,i2)%ek0*p378(2)-p378k0*cf526(i7,i2)%e(2)
      u378_526(i7,i2)%a(1)=u378_526(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u378_526(i7,i2)%a(2)=u378_526(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u378_526(i7,i2)%b(1)=u378_526(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u378_526(i7,i2)%b(2)=u378_526(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u378_526(i7,i2)%c(1)=u378_526(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u378_526(i7,i2)%c(2)=u378_526(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u378_526(i7,i2)%d(1)=u378_526(i7,i2)%d(1)+ccl*cf526(i7,i2)
     & %ek0
      u378_526(i7,i2)%d(2)=u378_526(i7,i2)%d(2)+ccr*cf526(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_78526(i5,i7,i2)%a,cc=l3_78526(i5,i7,i2)%c,a1=l3_78(i5)%a
* ,c1=l3_78(i5)%c,a2=u378_526(i7,i2)%a,b2=u378_526(i7,i2)%b,c2=u378_526(
* i7,i2)%c,d2=u378_526(i7,i2)%d,prq=s378,nsum=0
      l3_78526(i5,i7,i2)%a(1)=l3_78(i5)%a(1)*u378_526(i7,i2)%a(1
     & )+l3_78(i5)%c(1)*s378*u378_526(i7,i2)%b(2)
      l3_78526(i5,i7,i2)%c(1)=l3_78(i5)%a(1)*u378_526(i7,i2)%c(1
     & )+l3_78(i5)%c(1)*s378*u378_526(i7,i2)%d(2)
      l3_78526(i5,i7,i2)%c(2)=l3_78(i5)%c(2)*s378*u378_526(i7,i2
     & )%d(1)+l3_78(i5)%a(2)*u378_526(i7,i2)%c(2)
      l3_78526(i5,i7,i2)%a(2)=l3_78(i5)%c(2)*s378*u378_526(i7,i2
     & )%b(1)+l3_78(i5)%a(2)*u378_526(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u3256_78                                               
* quqd -- p=p3256,q=p41
      quqd=p3256(0)*p41(0)-p3256(1)*p41(1)-p3256(2)*p41(2)-p3256
     & (3)*p41(3)
      ccr=zcr(id3)/(f41)
      ccl=zcl(id3)/(f41)
      do i5=1,2
* T0 -- qu=p3256,qd=p41,v=cz78(i5)%e,a=u3526_78(i5)%a,b=u3526_78(i5)%b,c
* =u3526_78(i5)%c,d=u3526_78(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p3256(2)*p41(3)-p41(2)*p3256(3))+p32
     & 56k0*(cz78(i5)%e(2)*p41(3)-p41(2)*cz78(i5)%e(3))-p41k0*(c
     & z78(i5)%e(2)*p3256(3)-p3256(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p3256k0+p3256(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i5)%e(3)*p41k0+p41(3)*cz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i5)%e(0)*p3256(0)-cz78(i5)%e(1)*p3256(1)-cz78(i5
     & )%e(2)*p3256(2)-cz78(i5)%e(3)*p3256(3)
      cvqd=cz78(i5)%e(0)*p41(0)-cz78(i5)%e(1)*p41(1)-cz78(i5)%e(
     & 2)*p41(2)-cz78(i5)%e(3)*p41(3)
      cauxa=-cz78(i5)%ek0*quqd+p3256k0*cvqd+p41k0*cvqu
      cauxb=-cz78(i5)%ek0*p41(2)+p41k0*cz78(i5)%e(2)
      cauxc=+cz78(i5)%ek0*p3256(2)-p3256k0*cz78(i5)%e(2)
      u3526_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      u3526_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      u3526_78(i5)%b(1)=ccl*(cauxb-ceps_2)
      u3526_78(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u3526_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      u3526_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u3526_78(i5)%d(1)=ccl*cz78(i5)%ek0
      u3526_78(i5)%d(2)=ccr*cz78(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f41)
      ccl=fcl(id3)/(f41)
      do i5=1,2
* T0 -- qu=p3256,qd=p41,v=cf78(i5)%e,a=u3526_78(i5)%a,b=u3526_78(i5)%b,c
* =u3526_78(i5)%c,d=u3526_78(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p3256(2)*p41(3)-p41(2)*p3256(3))+p32
     & 56k0*(cf78(i5)%e(2)*p41(3)-p41(2)*cf78(i5)%e(3))-p41k0*(c
     & f78(i5)%e(2)*p3256(3)-p3256(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p3256k0+p3256(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i5)%e(3)*p41k0+p41(3)*cf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i5)%e(0)*p3256(0)-cf78(i5)%e(1)*p3256(1)-cf78(i5
     & )%e(2)*p3256(2)-cf78(i5)%e(3)*p3256(3)
      cvqd=cf78(i5)%e(0)*p41(0)-cf78(i5)%e(1)*p41(1)-cf78(i5)%e(
     & 2)*p41(2)-cf78(i5)%e(3)*p41(3)
      cauxa=-cf78(i5)%ek0*quqd+p3256k0*cvqd+p41k0*cvqu
      cauxb=-cf78(i5)%ek0*p41(2)+p41k0*cf78(i5)%e(2)
      cauxc=+cf78(i5)%ek0*p3256(2)-p3256k0*cf78(i5)%e(2)
      u3526_78(i5)%a(1)=u3526_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      u3526_78(i5)%a(2)=u3526_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      u3526_78(i5)%b(1)=u3526_78(i5)%b(1)+ccl*(cauxb-ceps_2)
      u3526_78(i5)%b(2)=u3526_78(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u3526_78(i5)%c(1)=u3526_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      u3526_78(i5)%c(2)=u3526_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u3526_78(i5)%d(1)=u3526_78(i5)%d(1)+ccl*cf78(i5)%ek0
      u3526_78(i5)%d(2)=u3526_78(i5)%d(2)+ccr*cf78(i5)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_78526(i5,i7,i2)%a,cc=l3_78526(i5,i7,i2)%c,a1=l3_526(i7,i
* 2)%a,c1=l3_526(i7,i2)%c,a2=u3526_78(i5)%a,b2=u3526_78(i5)%b,c2=u3526_7
* 8(i5)%c,d2=u3526_78(i5)%d,prq=s3256,nsum=1
      l3_78526(i5,i7,i2)%a(1)=l3_78526(i5,i7,i2)%a(1)+l3_526(i7,
     & i2)%a(1)*u3526_78(i5)%a(1)+l3_526(i7,i2)%c(1)*s3256*u3526
     & _78(i5)%b(2)
      l3_78526(i5,i7,i2)%c(1)=l3_78526(i5,i7,i2)%c(1)+l3_526(i7,
     & i2)%a(1)*u3526_78(i5)%c(1)+l3_526(i7,i2)%c(1)*s3256*u3526
     & _78(i5)%d(2)
      l3_78526(i5,i7,i2)%c(2)=l3_78526(i5,i7,i2)%c(2)+l3_526(i7,
     & i2)%c(2)*s3256*u3526_78(i5)%d(1)+l3_526(i7,i2)%a(2)*u3526
     & _78(i5)%c(2)
      l3_78526(i5,i7,i2)%a(2)=l3_78526(i5,i7,i2)%a(2)+l3_526(i7,
     & i2)%c(2)*s3256*u3526_78(i5)%b(1)+l3_526(i7,i2)%a(2)*u3526
     & _78(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id3).ne.1) then
* I                                                                     
* quqd -- p=p32,q=p3278
      quqd=p32(0)*p3278(0)-p32(1)*p3278(1)-p32(2)*p3278(2)-p32(3
     & )*p3278(3)
      ccr=zcr(id3)/(f3278)
      ccl=zcl(id3)/(f3278)
      do i5=1,2
* T0 -- qu=p32,qd=p3278,v=cz78(i5)%e,a=u32_78(i5)%a,b=u32_78(i5)%b,c=u32
* _78(i5)%c,d=u32_78(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p32(2)*p3278(3)-p3278(2)*p32(3))+p32
     & k0*(cz78(i5)%e(2)*p3278(3)-p3278(2)*cz78(i5)%e(3))-p3278k
     & 0*(cz78(i5)%e(2)*p32(3)-p32(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p32k0+p32(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i5)%e(3)*p3278k0+p3278(3)*cz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i5)%e(0)*p32(0)-cz78(i5)%e(1)*p32(1)-cz78(i5)%e(
     & 2)*p32(2)-cz78(i5)%e(3)*p32(3)
      cvqd=cz78(i5)%e(0)*p3278(0)-cz78(i5)%e(1)*p3278(1)-cz78(i5
     & )%e(2)*p3278(2)-cz78(i5)%e(3)*p3278(3)
      cauxa=-cz78(i5)%ek0*quqd+p32k0*cvqd+p3278k0*cvqu
      cauxb=-cz78(i5)%ek0*p3278(2)+p3278k0*cz78(i5)%e(2)
      cauxc=+cz78(i5)%ek0*p32(2)-p32k0*cz78(i5)%e(2)
      u32_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      u32_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      u32_78(i5)%b(1)=ccl*(cauxb-ceps_2)
      u32_78(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u32_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      u32_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u32_78(i5)%d(1)=ccl*cz78(i5)%ek0
      u32_78(i5)%d(2)=ccr*cz78(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f3278)
      ccl=fcl(id3)/(f3278)
      do i5=1,2
* T0 -- qu=p32,qd=p3278,v=cf78(i5)%e,a=u32_78(i5)%a,b=u32_78(i5)%b,c=u32
* _78(i5)%c,d=u32_78(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p32(2)*p3278(3)-p3278(2)*p32(3))+p32
     & k0*(cf78(i5)%e(2)*p3278(3)-p3278(2)*cf78(i5)%e(3))-p3278k
     & 0*(cf78(i5)%e(2)*p32(3)-p32(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p32k0+p32(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i5)%e(3)*p3278k0+p3278(3)*cf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i5)%e(0)*p32(0)-cf78(i5)%e(1)*p32(1)-cf78(i5)%e(
     & 2)*p32(2)-cf78(i5)%e(3)*p32(3)
      cvqd=cf78(i5)%e(0)*p3278(0)-cf78(i5)%e(1)*p3278(1)-cf78(i5
     & )%e(2)*p3278(2)-cf78(i5)%e(3)*p3278(3)
      cauxa=-cf78(i5)%ek0*quqd+p32k0*cvqd+p3278k0*cvqu
      cauxb=-cf78(i5)%ek0*p3278(2)+p3278k0*cf78(i5)%e(2)
      cauxc=+cf78(i5)%ek0*p32(2)-p32k0*cf78(i5)%e(2)
      u32_78(i5)%a(1)=u32_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      u32_78(i5)%a(2)=u32_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      u32_78(i5)%b(1)=u32_78(i5)%b(1)+ccl*(cauxb-ceps_2)
      u32_78(i5)%b(2)=u32_78(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u32_78(i5)%c(1)=u32_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      u32_78(i5)%c(2)=u32_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u32_78(i5)%d(1)=u32_78(i5)%d(1)+ccl*cf78(i5)%ek0
      u32_78(i5)%d(2)=u32_78(i5)%d(2)+ccr*cf78(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l3_278(i1,i5)%a,cc=l3_278(i1,i5)%c,a1=l3_2(i1)%a,c1=l3_2(i1
* )%c,a2=u32_78(i5)%a,b2=u32_78(i5)%b,c2=u32_78(i5)%c,d2=u32_78(i5)%d,pr
* q=s32,nsum=0
      l3_278(i1,i5)%a(1)=l3_2(i1)%a(1)*u32_78(i5)%a(1)+l3_2(i1)%
     & c(1)*s32*u32_78(i5)%b(2)
      l3_278(i1,i5)%c(1)=l3_2(i1)%a(1)*u32_78(i5)%c(1)+l3_2(i1)%
     & c(1)*s32*u32_78(i5)%d(2)
      l3_278(i1,i5)%c(2)=l3_2(i1)%c(2)*s32*u32_78(i5)%d(1)+l3_2(
     & i1)%a(2)*u32_78(i5)%c(2)
      l3_278(i1,i5)%a(2)=l3_2(i1)%c(2)*s32*u32_78(i5)%b(1)+l3_2(
     & i1)%a(2)*u32_78(i5)%a(2)
      end do
      end do
  
* quqd -- p=p378,q=p3278
      quqd=p378(0)*p3278(0)-p378(1)*p3278(1)-p378(2)*p3278(2)-p3
     & 78(3)*p3278(3)
      ccr=1.d0/(f3278)
      ccl=1.d0/(f3278)
      do i1=1,2
* T0 -- qu=p378,qd=p3278,v=ce2(i1)%e,a=u378_2(i1)%a,b=u378_2(i1)%b,c=u37
* 8_2(i1)%c,d=u378_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p378(2)*p3278(3)-p3278(2)*p378(3))+p3
     & 78k0*(ce2(i1)%e(2)*p3278(3)-p3278(2)*ce2(i1)%e(3))-p3278k
     & 0*(ce2(i1)%e(2)*p378(3)-p378(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p378k0+p378(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p3278k0+p3278(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p378(0)-ce2(i1)%e(1)*p378(1)-ce2(i1)%e(2
     & )*p378(2)-ce2(i1)%e(3)*p378(3)
      cvqd=ce2(i1)%e(0)*p3278(0)-ce2(i1)%e(1)*p3278(1)-ce2(i1)%e
     & (2)*p3278(2)-ce2(i1)%e(3)*p3278(3)
      cauxa=-ce2(i1)%ek0*quqd+p378k0*cvqd+p3278k0*cvqu
      cauxb=-ce2(i1)%ek0*p3278(2)+p3278k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p378(2)-p378k0*ce2(i1)%e(2)
      u378_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u378_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u378_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u378_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u378_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u378_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u378_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u378_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l3_278(i1,i5)%a,cc=l3_278(i1,i5)%c,a1=l3_78(i5)%a,c1=l3_78(
* i5)%c,a2=u378_2(i1)%a,b2=u378_2(i1)%b,c2=u378_2(i1)%c,d2=u378_2(i1)%d,
* prq=s378,nsum=1
      l3_278(i1,i5)%a(1)=l3_278(i1,i5)%a(1)+l3_78(i5)%a(1)*u378_
     & 2(i1)%a(1)+l3_78(i5)%c(1)*s378*u378_2(i1)%b(2)
      l3_278(i1,i5)%c(1)=l3_278(i1,i5)%c(1)+l3_78(i5)%a(1)*u378_
     & 2(i1)%c(1)+l3_78(i5)%c(1)*s378*u378_2(i1)%d(2)
      l3_278(i1,i5)%c(2)=l3_278(i1,i5)%c(2)+l3_78(i5)%c(2)*s378*
     & u378_2(i1)%d(1)+l3_78(i5)%a(2)*u378_2(i1)%c(2)
      l3_278(i1,i5)%a(2)=l3_278(i1,i5)%a(2)+l3_78(i5)%c(2)*s378*
     & u378_2(i1)%b(1)+l3_78(i5)%a(2)*u378_2(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id5).ne.1) then
* quqd -- p=p32,q=p478
      quqd=p32(0)*p478(0)-p32(1)*p478(1)-p32(2)*p478(2)-p32(3)*p
     & 478(3)
      ccr=zcr(id3)/(f478)
      ccl=zcl(id3)/(f478)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p32,qd=p478,v=cz516(i7,i2)%e,a=u32_516(i7,i2)%a,b=u32_516(i7,
* i2)%b,c=u32_516(i7,i2)%c,d=u32_516(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz516(i7,i2)%ek0*(p32(2)*p478(3)-p478(2)*p32(3))+p
     & 32k0*(cz516(i7,i2)%e(2)*p478(3)-p478(2)*cz516(i7,i2)%e(3)
     & )-p478k0*(cz516(i7,i2)%e(2)*p32(3)-p32(2)*cz516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz516(i7,i2)%e(3)*p32k0+p32(3)*cz516(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz516(i7,i2)%e(3)*p478k0+p478(3)*cz516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz516(i7,i2)%e(0)*p32(0)-cz516(i7,i2)%e(1)*p32(1)-cz5
     & 16(i7,i2)%e(2)*p32(2)-cz516(i7,i2)%e(3)*p32(3)
      cvqd=cz516(i7,i2)%e(0)*p478(0)-cz516(i7,i2)%e(1)*p478(1)-c
     & z516(i7,i2)%e(2)*p478(2)-cz516(i7,i2)%e(3)*p478(3)
      cauxa=-cz516(i7,i2)%ek0*quqd+p32k0*cvqd+p478k0*cvqu
      cauxb=-cz516(i7,i2)%ek0*p478(2)+p478k0*cz516(i7,i2)%e(2)
      cauxc=+cz516(i7,i2)%ek0*p32(2)-p32k0*cz516(i7,i2)%e(2)
      u32_516(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u32_516(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u32_516(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u32_516(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u32_516(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u32_516(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u32_516(i7,i2)%d(1)=ccl*cz516(i7,i2)%ek0
      u32_516(i7,i2)%d(2)=ccr*cz516(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id3)/(f478)
      ccl=fcl(id3)/(f478)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p32,qd=p478,v=cf516(i7,i2)%e,a=u32_516(i7,i2)%a,b=u32_516(i7,
* i2)%b,c=u32_516(i7,i2)%c,d=u32_516(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf516(i7,i2)%ek0*(p32(2)*p478(3)-p478(2)*p32(3))+p
     & 32k0*(cf516(i7,i2)%e(2)*p478(3)-p478(2)*cf516(i7,i2)%e(3)
     & )-p478k0*(cf516(i7,i2)%e(2)*p32(3)-p32(2)*cf516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf516(i7,i2)%e(3)*p32k0+p32(3)*cf516(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf516(i7,i2)%e(3)*p478k0+p478(3)*cf516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf516(i7,i2)%e(0)*p32(0)-cf516(i7,i2)%e(1)*p32(1)-cf5
     & 16(i7,i2)%e(2)*p32(2)-cf516(i7,i2)%e(3)*p32(3)
      cvqd=cf516(i7,i2)%e(0)*p478(0)-cf516(i7,i2)%e(1)*p478(1)-c
     & f516(i7,i2)%e(2)*p478(2)-cf516(i7,i2)%e(3)*p478(3)
      cauxa=-cf516(i7,i2)%ek0*quqd+p32k0*cvqd+p478k0*cvqu
      cauxb=-cf516(i7,i2)%ek0*p478(2)+p478k0*cf516(i7,i2)%e(2)
      cauxc=+cf516(i7,i2)%ek0*p32(2)-p32k0*cf516(i7,i2)%e(2)
      u32_516(i7,i2)%a(1)=u32_516(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u32_516(i7,i2)%a(2)=u32_516(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u32_516(i7,i2)%b(1)=u32_516(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u32_516(i7,i2)%b(2)=u32_516(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u32_516(i7,i2)%c(1)=u32_516(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u32_516(i7,i2)%c(2)=u32_516(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u32_516(i7,i2)%d(1)=u32_516(i7,i2)%d(1)+ccl*cf516(i7,i2)%e
     & k0
      u32_516(i7,i2)%d(2)=u32_516(i7,i2)%d(2)+ccr*cf516(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_2516(i1,i7,i2)%a,cc=l3_2516(i1,i7,i2)%c,a1=l3_2(i1)%a,c1
* =l3_2(i1)%c,a2=u32_516(i7,i2)%a,b2=u32_516(i7,i2)%b,c2=u32_516(i7,i2)%
* c,d2=u32_516(i7,i2)%d,prq=s32,nsum=0
      l3_2516(i1,i7,i2)%a(1)=l3_2(i1)%a(1)*u32_516(i7,i2)%a(1)+l
     & 3_2(i1)%c(1)*s32*u32_516(i7,i2)%b(2)
      l3_2516(i1,i7,i2)%c(1)=l3_2(i1)%a(1)*u32_516(i7,i2)%c(1)+l
     & 3_2(i1)%c(1)*s32*u32_516(i7,i2)%d(2)
      l3_2516(i1,i7,i2)%c(2)=l3_2(i1)%c(2)*s32*u32_516(i7,i2)%d(
     & 1)+l3_2(i1)%a(2)*u32_516(i7,i2)%c(2)
      l3_2516(i1,i7,i2)%a(2)=l3_2(i1)%c(2)*s32*u32_516(i7,i2)%b(
     & 1)+l3_2(i1)%a(2)*u32_516(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u3156_2                                                
* quqd -- p=p3156,q=p478
      quqd=p3156(0)*p478(0)-p3156(1)*p478(1)-p3156(2)*p478(2)-p3
     & 156(3)*p478(3)
      ccr=1.d0/(f478)
      ccl=1.d0/(f478)
      do i1=1,2
* T0 -- qu=p3156,qd=p478,v=ce2(i1)%e,a=u3516_2(i1)%a,b=u3516_2(i1)%b,c=u
* 3516_2(i1)%c,d=u3516_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p3156(2)*p478(3)-p478(2)*p3156(3))+p3
     & 156k0*(ce2(i1)%e(2)*p478(3)-p478(2)*ce2(i1)%e(3))-p478k0*
     & (ce2(i1)%e(2)*p3156(3)-p3156(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p3156k0+p3156(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p478k0+p478(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p3156(0)-ce2(i1)%e(1)*p3156(1)-ce2(i1)%e
     & (2)*p3156(2)-ce2(i1)%e(3)*p3156(3)
      cvqd=ce2(i1)%e(0)*p478(0)-ce2(i1)%e(1)*p478(1)-ce2(i1)%e(2
     & )*p478(2)-ce2(i1)%e(3)*p478(3)
      cauxa=-ce2(i1)%ek0*quqd+p3156k0*cvqd+p478k0*cvqu
      cauxb=-ce2(i1)%ek0*p478(2)+p478k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p3156(2)-p3156k0*ce2(i1)%e(2)
      u3516_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u3516_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u3516_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u3516_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u3516_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u3516_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u3516_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u3516_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_2516(i1,i7,i2)%a,cc=l3_2516(i1,i7,i2)%c,a1=l3_516(i7,i2)
* %a,c1=l3_516(i7,i2)%c,a2=u3516_2(i1)%a,b2=u3516_2(i1)%b,c2=u3516_2(i1)
* %c,d2=u3516_2(i1)%d,prq=s3156,nsum=1
      l3_2516(i1,i7,i2)%a(1)=l3_2516(i1,i7,i2)%a(1)+l3_516(i7,i2
     & )%a(1)*u3516_2(i1)%a(1)+l3_516(i7,i2)%c(1)*s3156*u3516_2(
     & i1)%b(2)
      l3_2516(i1,i7,i2)%c(1)=l3_2516(i1,i7,i2)%c(1)+l3_516(i7,i2
     & )%a(1)*u3516_2(i1)%c(1)+l3_516(i7,i2)%c(1)*s3156*u3516_2(
     & i1)%d(2)
      l3_2516(i1,i7,i2)%c(2)=l3_2516(i1,i7,i2)%c(2)+l3_516(i7,i2
     & )%c(2)*s3156*u3516_2(i1)%d(1)+l3_516(i7,i2)%a(2)*u3516_2(
     & i1)%c(2)
      l3_2516(i1,i7,i2)%a(2)=l3_2516(i1,i7,i2)%a(2)+l3_516(i7,i2
     & )%c(2)*s3156*u3516_2(i1)%b(1)+l3_516(i7,i2)%a(2)*u3516_2(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p378,q=p42
      quqd=p378(0)*p42(0)-p378(1)*p42(1)-p378(2)*p42(2)-p378(3)*
     & p42(3)
      ccr=zcr(id3)/(f42)
      ccl=zcl(id3)/(f42)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p378,qd=p42,v=cz516(i7,i2)%e,a=u378_516(i7,i2)%a,b=u378_516(i
* 7,i2)%b,c=u378_516(i7,i2)%c,d=u378_516(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz516(i7,i2)%ek0*(p378(2)*p42(3)-p42(2)*p378(3))+p
     & 378k0*(cz516(i7,i2)%e(2)*p42(3)-p42(2)*cz516(i7,i2)%e(3))
     & -p42k0*(cz516(i7,i2)%e(2)*p378(3)-p378(2)*cz516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz516(i7,i2)%e(3)*p378k0+p378(3)*cz516(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz516(i7,i2)%e(3)*p42k0+p42(3)*cz516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz516(i7,i2)%e(0)*p378(0)-cz516(i7,i2)%e(1)*p378(1)-c
     & z516(i7,i2)%e(2)*p378(2)-cz516(i7,i2)%e(3)*p378(3)
      cvqd=cz516(i7,i2)%e(0)*p42(0)-cz516(i7,i2)%e(1)*p42(1)-cz5
     & 16(i7,i2)%e(2)*p42(2)-cz516(i7,i2)%e(3)*p42(3)
      cauxa=-cz516(i7,i2)%ek0*quqd+p378k0*cvqd+p42k0*cvqu
      cauxb=-cz516(i7,i2)%ek0*p42(2)+p42k0*cz516(i7,i2)%e(2)
      cauxc=+cz516(i7,i2)%ek0*p378(2)-p378k0*cz516(i7,i2)%e(2)
      u378_516(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u378_516(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u378_516(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u378_516(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u378_516(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u378_516(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u378_516(i7,i2)%d(1)=ccl*cz516(i7,i2)%ek0
      u378_516(i7,i2)%d(2)=ccr*cz516(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id3)/(f42)
      ccl=fcl(id3)/(f42)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p378,qd=p42,v=cf516(i7,i2)%e,a=u378_516(i7,i2)%a,b=u378_516(i
* 7,i2)%b,c=u378_516(i7,i2)%c,d=u378_516(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf516(i7,i2)%ek0*(p378(2)*p42(3)-p42(2)*p378(3))+p
     & 378k0*(cf516(i7,i2)%e(2)*p42(3)-p42(2)*cf516(i7,i2)%e(3))
     & -p42k0*(cf516(i7,i2)%e(2)*p378(3)-p378(2)*cf516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf516(i7,i2)%e(3)*p378k0+p378(3)*cf516(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf516(i7,i2)%e(3)*p42k0+p42(3)*cf516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf516(i7,i2)%e(0)*p378(0)-cf516(i7,i2)%e(1)*p378(1)-c
     & f516(i7,i2)%e(2)*p378(2)-cf516(i7,i2)%e(3)*p378(3)
      cvqd=cf516(i7,i2)%e(0)*p42(0)-cf516(i7,i2)%e(1)*p42(1)-cf5
     & 16(i7,i2)%e(2)*p42(2)-cf516(i7,i2)%e(3)*p42(3)
      cauxa=-cf516(i7,i2)%ek0*quqd+p378k0*cvqd+p42k0*cvqu
      cauxb=-cf516(i7,i2)%ek0*p42(2)+p42k0*cf516(i7,i2)%e(2)
      cauxc=+cf516(i7,i2)%ek0*p378(2)-p378k0*cf516(i7,i2)%e(2)
      u378_516(i7,i2)%a(1)=u378_516(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u378_516(i7,i2)%a(2)=u378_516(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u378_516(i7,i2)%b(1)=u378_516(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u378_516(i7,i2)%b(2)=u378_516(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u378_516(i7,i2)%c(1)=u378_516(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u378_516(i7,i2)%c(2)=u378_516(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u378_516(i7,i2)%d(1)=u378_516(i7,i2)%d(1)+ccl*cf516(i7,i2)
     & %ek0
      u378_516(i7,i2)%d(2)=u378_516(i7,i2)%d(2)+ccr*cf516(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_78516(i5,i7,i2)%a,cc=l3_78516(i5,i7,i2)%c,a1=l3_78(i5)%a
* ,c1=l3_78(i5)%c,a2=u378_516(i7,i2)%a,b2=u378_516(i7,i2)%b,c2=u378_516(
* i7,i2)%c,d2=u378_516(i7,i2)%d,prq=s378,nsum=0
      l3_78516(i5,i7,i2)%a(1)=l3_78(i5)%a(1)*u378_516(i7,i2)%a(1
     & )+l3_78(i5)%c(1)*s378*u378_516(i7,i2)%b(2)
      l3_78516(i5,i7,i2)%c(1)=l3_78(i5)%a(1)*u378_516(i7,i2)%c(1
     & )+l3_78(i5)%c(1)*s378*u378_516(i7,i2)%d(2)
      l3_78516(i5,i7,i2)%c(2)=l3_78(i5)%c(2)*s378*u378_516(i7,i2
     & )%d(1)+l3_78(i5)%a(2)*u378_516(i7,i2)%c(2)
      l3_78516(i5,i7,i2)%a(2)=l3_78(i5)%c(2)*s378*u378_516(i7,i2
     & )%b(1)+l3_78(i5)%a(2)*u378_516(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u3156_78                                               
* quqd -- p=p3156,q=p42
      quqd=p3156(0)*p42(0)-p3156(1)*p42(1)-p3156(2)*p42(2)-p3156
     & (3)*p42(3)
      ccr=zcr(id3)/(f42)
      ccl=zcl(id3)/(f42)
      do i5=1,2
* T0 -- qu=p3156,qd=p42,v=cz78(i5)%e,a=u3516_78(i5)%a,b=u3516_78(i5)%b,c
* =u3516_78(i5)%c,d=u3516_78(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p3156(2)*p42(3)-p42(2)*p3156(3))+p31
     & 56k0*(cz78(i5)%e(2)*p42(3)-p42(2)*cz78(i5)%e(3))-p42k0*(c
     & z78(i5)%e(2)*p3156(3)-p3156(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p3156k0+p3156(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i5)%e(3)*p42k0+p42(3)*cz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i5)%e(0)*p3156(0)-cz78(i5)%e(1)*p3156(1)-cz78(i5
     & )%e(2)*p3156(2)-cz78(i5)%e(3)*p3156(3)
      cvqd=cz78(i5)%e(0)*p42(0)-cz78(i5)%e(1)*p42(1)-cz78(i5)%e(
     & 2)*p42(2)-cz78(i5)%e(3)*p42(3)
      cauxa=-cz78(i5)%ek0*quqd+p3156k0*cvqd+p42k0*cvqu
      cauxb=-cz78(i5)%ek0*p42(2)+p42k0*cz78(i5)%e(2)
      cauxc=+cz78(i5)%ek0*p3156(2)-p3156k0*cz78(i5)%e(2)
      u3516_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      u3516_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      u3516_78(i5)%b(1)=ccl*(cauxb-ceps_2)
      u3516_78(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u3516_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      u3516_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u3516_78(i5)%d(1)=ccl*cz78(i5)%ek0
      u3516_78(i5)%d(2)=ccr*cz78(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f42)
      ccl=fcl(id3)/(f42)
      do i5=1,2
* T0 -- qu=p3156,qd=p42,v=cf78(i5)%e,a=u3516_78(i5)%a,b=u3516_78(i5)%b,c
* =u3516_78(i5)%c,d=u3516_78(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p3156(2)*p42(3)-p42(2)*p3156(3))+p31
     & 56k0*(cf78(i5)%e(2)*p42(3)-p42(2)*cf78(i5)%e(3))-p42k0*(c
     & f78(i5)%e(2)*p3156(3)-p3156(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p3156k0+p3156(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i5)%e(3)*p42k0+p42(3)*cf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i5)%e(0)*p3156(0)-cf78(i5)%e(1)*p3156(1)-cf78(i5
     & )%e(2)*p3156(2)-cf78(i5)%e(3)*p3156(3)
      cvqd=cf78(i5)%e(0)*p42(0)-cf78(i5)%e(1)*p42(1)-cf78(i5)%e(
     & 2)*p42(2)-cf78(i5)%e(3)*p42(3)
      cauxa=-cf78(i5)%ek0*quqd+p3156k0*cvqd+p42k0*cvqu
      cauxb=-cf78(i5)%ek0*p42(2)+p42k0*cf78(i5)%e(2)
      cauxc=+cf78(i5)%ek0*p3156(2)-p3156k0*cf78(i5)%e(2)
      u3516_78(i5)%a(1)=u3516_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      u3516_78(i5)%a(2)=u3516_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      u3516_78(i5)%b(1)=u3516_78(i5)%b(1)+ccl*(cauxb-ceps_2)
      u3516_78(i5)%b(2)=u3516_78(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u3516_78(i5)%c(1)=u3516_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      u3516_78(i5)%c(2)=u3516_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u3516_78(i5)%d(1)=u3516_78(i5)%d(1)+ccl*cf78(i5)%ek0
      u3516_78(i5)%d(2)=u3516_78(i5)%d(2)+ccr*cf78(i5)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l3_78516(i5,i7,i2)%a,cc=l3_78516(i5,i7,i2)%c,a1=l3_516(i7,i
* 2)%a,c1=l3_516(i7,i2)%c,a2=u3516_78(i5)%a,b2=u3516_78(i5)%b,c2=u3516_7
* 8(i5)%c,d2=u3516_78(i5)%d,prq=s3156,nsum=1
      l3_78516(i5,i7,i2)%a(1)=l3_78516(i5,i7,i2)%a(1)+l3_516(i7,
     & i2)%a(1)*u3516_78(i5)%a(1)+l3_516(i7,i2)%c(1)*s3156*u3516
     & _78(i5)%b(2)
      l3_78516(i5,i7,i2)%c(1)=l3_78516(i5,i7,i2)%c(1)+l3_516(i7,
     & i2)%a(1)*u3516_78(i5)%c(1)+l3_516(i7,i2)%c(1)*s3156*u3516
     & _78(i5)%d(2)
      l3_78516(i5,i7,i2)%c(2)=l3_78516(i5,i7,i2)%c(2)+l3_516(i7,
     & i2)%c(2)*s3156*u3516_78(i5)%d(1)+l3_516(i7,i2)%a(2)*u3516
     & _78(i5)%c(2)
      l3_78516(i5,i7,i2)%a(2)=l3_78516(i5,i7,i2)%a(2)+l3_516(i7,
     & i2)%c(2)*s3156*u3516_78(i5)%b(1)+l3_516(i7,i2)%a(2)*u3516
     & _78(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p51,q=p5134
      quqd=p51(0)*p5134(0)-p51(1)*p5134(1)-p51(2)*p5134(2)-p51(3
     & )*p5134(3)
      ccr=zcr(id5)/(f5134)
      ccl=zcl(id5)/(f5134)
      do i5=1,2
* T0 -- qu=p51,qd=p5134,v=cz34(i5)%e,a=u51_34(i5)%a,b=u51_34(i5)%b,c=u51
* _34(i5)%c,d=u51_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p51(2)*p5134(3)-p5134(2)*p51(3))+p51
     & k0*(cz34(i5)%e(2)*p5134(3)-p5134(2)*cz34(i5)%e(3))-p5134k
     & 0*(cz34(i5)%e(2)*p51(3)-p51(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p51k0+p51(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p5134k0+p5134(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p51(0)-cz34(i5)%e(1)*p51(1)-cz34(i5)%e(
     & 2)*p51(2)-cz34(i5)%e(3)*p51(3)
      cvqd=cz34(i5)%e(0)*p5134(0)-cz34(i5)%e(1)*p5134(1)-cz34(i5
     & )%e(2)*p5134(2)-cz34(i5)%e(3)*p5134(3)
      cauxa=-cz34(i5)%ek0*quqd+p51k0*cvqd+p5134k0*cvqu
      cauxb=-cz34(i5)%ek0*p5134(2)+p5134k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p51(2)-p51k0*cz34(i5)%e(2)
      u51_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u51_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u51_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u51_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u51_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u51_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u51_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u51_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f5134)
      ccl=fcl(id5)/(f5134)
      do i5=1,2
* T0 -- qu=p51,qd=p5134,v=cf34(i5)%e,a=u51_34(i5)%a,b=u51_34(i5)%b,c=u51
* _34(i5)%c,d=u51_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p51(2)*p5134(3)-p5134(2)*p51(3))+p51
     & k0*(cf34(i5)%e(2)*p5134(3)-p5134(2)*cf34(i5)%e(3))-p5134k
     & 0*(cf34(i5)%e(2)*p51(3)-p51(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p51k0+p51(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p5134k0+p5134(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p51(0)-cf34(i5)%e(1)*p51(1)-cf34(i5)%e(
     & 2)*p51(2)-cf34(i5)%e(3)*p51(3)
      cvqd=cf34(i5)%e(0)*p5134(0)-cf34(i5)%e(1)*p5134(1)-cf34(i5
     & )%e(2)*p5134(2)-cf34(i5)%e(3)*p5134(3)
      cauxa=-cf34(i5)%ek0*quqd+p51k0*cvqd+p5134k0*cvqu
      cauxb=-cf34(i5)%ek0*p5134(2)+p5134k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p51(2)-p51k0*cf34(i5)%e(2)
      u51_34(i5)%a(1)=u51_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u51_34(i5)%a(2)=u51_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u51_34(i5)%b(1)=u51_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u51_34(i5)%b(2)=u51_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u51_34(i5)%c(1)=u51_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u51_34(i5)%c(2)=u51_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u51_34(i5)%d(1)=u51_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u51_34(i5)%d(2)=u51_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l5_134(i1,i5)%a,cc=l5_134(i1,i5)%c,a1=l5_1(i1)%a,c1=l5_1(i1
* )%c,a2=u51_34(i5)%a,b2=u51_34(i5)%b,c2=u51_34(i5)%c,d2=u51_34(i5)%d,pr
* q=s51,nsum=0
      l5_134(i1,i5)%a(1)=l5_1(i1)%a(1)*u51_34(i5)%a(1)+l5_1(i1)%
     & c(1)*s51*u51_34(i5)%b(2)
      l5_134(i1,i5)%c(1)=l5_1(i1)%a(1)*u51_34(i5)%c(1)+l5_1(i1)%
     & c(1)*s51*u51_34(i5)%d(2)
      l5_134(i1,i5)%c(2)=l5_1(i1)%c(2)*s51*u51_34(i5)%d(1)+l5_1(
     & i1)%a(2)*u51_34(i5)%c(2)
      l5_134(i1,i5)%a(2)=l5_1(i1)%c(2)*s51*u51_34(i5)%b(1)+l5_1(
     & i1)%a(2)*u51_34(i5)%a(2)
      end do
      end do
  
* quqd -- p=p534,q=p5134
      quqd=p534(0)*p5134(0)-p534(1)*p5134(1)-p534(2)*p5134(2)-p5
     & 34(3)*p5134(3)
      ccr=1.d0/(f5134)
      ccl=1.d0/(f5134)
      do i1=1,2
* T0 -- qu=p534,qd=p5134,v=ce1(i1)%e,a=u534_1(i1)%a,b=u534_1(i1)%b,c=u53
* 4_1(i1)%c,d=u534_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p534(2)*p5134(3)-p5134(2)*p534(3))+p5
     & 34k0*(ce1(i1)%e(2)*p5134(3)-p5134(2)*ce1(i1)%e(3))-p5134k
     & 0*(ce1(i1)%e(2)*p534(3)-p534(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p534k0+p534(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p5134k0+p5134(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p534(0)-ce1(i1)%e(1)*p534(1)-ce1(i1)%e(2
     & )*p534(2)-ce1(i1)%e(3)*p534(3)
      cvqd=ce1(i1)%e(0)*p5134(0)-ce1(i1)%e(1)*p5134(1)-ce1(i1)%e
     & (2)*p5134(2)-ce1(i1)%e(3)*p5134(3)
      cauxa=-ce1(i1)%ek0*quqd+p534k0*cvqd+p5134k0*cvqu
      cauxb=-ce1(i1)%ek0*p5134(2)+p5134k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p534(2)-p534k0*ce1(i1)%e(2)
      u534_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u534_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u534_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u534_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u534_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u534_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u534_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u534_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l5_134(i1,i5)%a,cc=l5_134(i1,i5)%c,a1=l5_34(i5)%a,c1=l5_34(
* i5)%c,a2=u534_1(i1)%a,b2=u534_1(i1)%b,c2=u534_1(i1)%c,d2=u534_1(i1)%d,
* prq=s534,nsum=1
      l5_134(i1,i5)%a(1)=l5_134(i1,i5)%a(1)+l5_34(i5)%a(1)*u534_
     & 1(i1)%a(1)+l5_34(i5)%c(1)*s534*u534_1(i1)%b(2)
      l5_134(i1,i5)%c(1)=l5_134(i1,i5)%c(1)+l5_34(i5)%a(1)*u534_
     & 1(i1)%c(1)+l5_34(i5)%c(1)*s534*u534_1(i1)%d(2)
      l5_134(i1,i5)%c(2)=l5_134(i1,i5)%c(2)+l5_34(i5)%c(2)*s534*
     & u534_1(i1)%d(1)+l5_34(i5)%a(2)*u534_1(i1)%c(2)
      l5_134(i1,i5)%a(2)=l5_134(i1,i5)%a(2)+l5_34(i5)%c(2)*s534*
     & u534_1(i1)%b(1)+l5_34(i5)%a(2)*u534_1(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id7).ne.1) then
* quqd -- p=p51,q=p634
      quqd=p51(0)*p634(0)-p51(1)*p634(1)-p51(2)*p634(2)-p51(3)*p
     & 634(3)
      ccr=zcr(id5)/(f634)
      ccl=zcl(id5)/(f634)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p51,qd=p634,v=cz728(i7,i2)%e,a=u51_728(i7,i2)%a,b=u51_728(i7,
* i2)%b,c=u51_728(i7,i2)%c,d=u51_728(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz728(i7,i2)%ek0*(p51(2)*p634(3)-p634(2)*p51(3))+p
     & 51k0*(cz728(i7,i2)%e(2)*p634(3)-p634(2)*cz728(i7,i2)%e(3)
     & )-p634k0*(cz728(i7,i2)%e(2)*p51(3)-p51(2)*cz728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz728(i7,i2)%e(3)*p51k0+p51(3)*cz728(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz728(i7,i2)%e(3)*p634k0+p634(3)*cz728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz728(i7,i2)%e(0)*p51(0)-cz728(i7,i2)%e(1)*p51(1)-cz7
     & 28(i7,i2)%e(2)*p51(2)-cz728(i7,i2)%e(3)*p51(3)
      cvqd=cz728(i7,i2)%e(0)*p634(0)-cz728(i7,i2)%e(1)*p634(1)-c
     & z728(i7,i2)%e(2)*p634(2)-cz728(i7,i2)%e(3)*p634(3)
      cauxa=-cz728(i7,i2)%ek0*quqd+p51k0*cvqd+p634k0*cvqu
      cauxb=-cz728(i7,i2)%ek0*p634(2)+p634k0*cz728(i7,i2)%e(2)
      cauxc=+cz728(i7,i2)%ek0*p51(2)-p51k0*cz728(i7,i2)%e(2)
      u51_728(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u51_728(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u51_728(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u51_728(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u51_728(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u51_728(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u51_728(i7,i2)%d(1)=ccl*cz728(i7,i2)%ek0
      u51_728(i7,i2)%d(2)=ccr*cz728(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id5)/(f634)
      ccl=fcl(id5)/(f634)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p51,qd=p634,v=cf728(i7,i2)%e,a=u51_728(i7,i2)%a,b=u51_728(i7,
* i2)%b,c=u51_728(i7,i2)%c,d=u51_728(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf728(i7,i2)%ek0*(p51(2)*p634(3)-p634(2)*p51(3))+p
     & 51k0*(cf728(i7,i2)%e(2)*p634(3)-p634(2)*cf728(i7,i2)%e(3)
     & )-p634k0*(cf728(i7,i2)%e(2)*p51(3)-p51(2)*cf728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf728(i7,i2)%e(3)*p51k0+p51(3)*cf728(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf728(i7,i2)%e(3)*p634k0+p634(3)*cf728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf728(i7,i2)%e(0)*p51(0)-cf728(i7,i2)%e(1)*p51(1)-cf7
     & 28(i7,i2)%e(2)*p51(2)-cf728(i7,i2)%e(3)*p51(3)
      cvqd=cf728(i7,i2)%e(0)*p634(0)-cf728(i7,i2)%e(1)*p634(1)-c
     & f728(i7,i2)%e(2)*p634(2)-cf728(i7,i2)%e(3)*p634(3)
      cauxa=-cf728(i7,i2)%ek0*quqd+p51k0*cvqd+p634k0*cvqu
      cauxb=-cf728(i7,i2)%ek0*p634(2)+p634k0*cf728(i7,i2)%e(2)
      cauxc=+cf728(i7,i2)%ek0*p51(2)-p51k0*cf728(i7,i2)%e(2)
      u51_728(i7,i2)%a(1)=u51_728(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u51_728(i7,i2)%a(2)=u51_728(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u51_728(i7,i2)%b(1)=u51_728(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u51_728(i7,i2)%b(2)=u51_728(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u51_728(i7,i2)%c(1)=u51_728(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u51_728(i7,i2)%c(2)=u51_728(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u51_728(i7,i2)%d(1)=u51_728(i7,i2)%d(1)+ccl*cf728(i7,i2)%e
     & k0
      u51_728(i7,i2)%d(2)=u51_728(i7,i2)%d(2)+ccr*cf728(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_1728(i1,i7,i2)%a,cc=l5_1728(i1,i7,i2)%c,a1=l5_1(i1)%a,c1
* =l5_1(i1)%c,a2=u51_728(i7,i2)%a,b2=u51_728(i7,i2)%b,c2=u51_728(i7,i2)%
* c,d2=u51_728(i7,i2)%d,prq=s51,nsum=0
      l5_1728(i1,i7,i2)%a(1)=l5_1(i1)%a(1)*u51_728(i7,i2)%a(1)+l
     & 5_1(i1)%c(1)*s51*u51_728(i7,i2)%b(2)
      l5_1728(i1,i7,i2)%c(1)=l5_1(i1)%a(1)*u51_728(i7,i2)%c(1)+l
     & 5_1(i1)%c(1)*s51*u51_728(i7,i2)%d(2)
      l5_1728(i1,i7,i2)%c(2)=l5_1(i1)%c(2)*s51*u51_728(i7,i2)%d(
     & 1)+l5_1(i1)%a(2)*u51_728(i7,i2)%c(2)
      l5_1728(i1,i7,i2)%a(2)=l5_1(i1)%c(2)*s51*u51_728(i7,i2)%b(
     & 1)+l5_1(i1)%a(2)*u51_728(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u5278_1                                                
* quqd -- p=p5278,q=p634
      quqd=p5278(0)*p634(0)-p5278(1)*p634(1)-p5278(2)*p634(2)-p5
     & 278(3)*p634(3)
      ccr=1.d0/(f634)
      ccl=1.d0/(f634)
      do i1=1,2
* T0 -- qu=p5278,qd=p634,v=ce1(i1)%e,a=u5728_1(i1)%a,b=u5728_1(i1)%b,c=u
* 5728_1(i1)%c,d=u5728_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p5278(2)*p634(3)-p634(2)*p5278(3))+p5
     & 278k0*(ce1(i1)%e(2)*p634(3)-p634(2)*ce1(i1)%e(3))-p634k0*
     & (ce1(i1)%e(2)*p5278(3)-p5278(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p5278k0+p5278(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p634k0+p634(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p5278(0)-ce1(i1)%e(1)*p5278(1)-ce1(i1)%e
     & (2)*p5278(2)-ce1(i1)%e(3)*p5278(3)
      cvqd=ce1(i1)%e(0)*p634(0)-ce1(i1)%e(1)*p634(1)-ce1(i1)%e(2
     & )*p634(2)-ce1(i1)%e(3)*p634(3)
      cauxa=-ce1(i1)%ek0*quqd+p5278k0*cvqd+p634k0*cvqu
      cauxb=-ce1(i1)%ek0*p634(2)+p634k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p5278(2)-p5278k0*ce1(i1)%e(2)
      u5728_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u5728_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u5728_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u5728_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u5728_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u5728_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u5728_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u5728_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_1728(i1,i7,i2)%a,cc=l5_1728(i1,i7,i2)%c,a1=l5_728(i7,i2)
* %a,c1=l5_728(i7,i2)%c,a2=u5728_1(i1)%a,b2=u5728_1(i1)%b,c2=u5728_1(i1)
* %c,d2=u5728_1(i1)%d,prq=s5278,nsum=1
      l5_1728(i1,i7,i2)%a(1)=l5_1728(i1,i7,i2)%a(1)+l5_728(i7,i2
     & )%a(1)*u5728_1(i1)%a(1)+l5_728(i7,i2)%c(1)*s5278*u5728_1(
     & i1)%b(2)
      l5_1728(i1,i7,i2)%c(1)=l5_1728(i1,i7,i2)%c(1)+l5_728(i7,i2
     & )%a(1)*u5728_1(i1)%c(1)+l5_728(i7,i2)%c(1)*s5278*u5728_1(
     & i1)%d(2)
      l5_1728(i1,i7,i2)%c(2)=l5_1728(i1,i7,i2)%c(2)+l5_728(i7,i2
     & )%c(2)*s5278*u5728_1(i1)%d(1)+l5_728(i7,i2)%a(2)*u5728_1(
     & i1)%c(2)
      l5_1728(i1,i7,i2)%a(2)=l5_1728(i1,i7,i2)%a(2)+l5_728(i7,i2
     & )%c(2)*s5278*u5728_1(i1)%b(1)+l5_728(i7,i2)%a(2)*u5728_1(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p534,q=p61
      quqd=p534(0)*p61(0)-p534(1)*p61(1)-p534(2)*p61(2)-p534(3)*
     & p61(3)
      ccr=zcr(id5)/(f61)
      ccl=zcl(id5)/(f61)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p534,qd=p61,v=cz728(i7,i2)%e,a=u534_728(i7,i2)%a,b=u534_728(i
* 7,i2)%b,c=u534_728(i7,i2)%c,d=u534_728(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz728(i7,i2)%ek0*(p534(2)*p61(3)-p61(2)*p534(3))+p
     & 534k0*(cz728(i7,i2)%e(2)*p61(3)-p61(2)*cz728(i7,i2)%e(3))
     & -p61k0*(cz728(i7,i2)%e(2)*p534(3)-p534(2)*cz728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz728(i7,i2)%e(3)*p534k0+p534(3)*cz728(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz728(i7,i2)%e(3)*p61k0+p61(3)*cz728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz728(i7,i2)%e(0)*p534(0)-cz728(i7,i2)%e(1)*p534(1)-c
     & z728(i7,i2)%e(2)*p534(2)-cz728(i7,i2)%e(3)*p534(3)
      cvqd=cz728(i7,i2)%e(0)*p61(0)-cz728(i7,i2)%e(1)*p61(1)-cz7
     & 28(i7,i2)%e(2)*p61(2)-cz728(i7,i2)%e(3)*p61(3)
      cauxa=-cz728(i7,i2)%ek0*quqd+p534k0*cvqd+p61k0*cvqu
      cauxb=-cz728(i7,i2)%ek0*p61(2)+p61k0*cz728(i7,i2)%e(2)
      cauxc=+cz728(i7,i2)%ek0*p534(2)-p534k0*cz728(i7,i2)%e(2)
      u534_728(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u534_728(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u534_728(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u534_728(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u534_728(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u534_728(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u534_728(i7,i2)%d(1)=ccl*cz728(i7,i2)%ek0
      u534_728(i7,i2)%d(2)=ccr*cz728(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id5)/(f61)
      ccl=fcl(id5)/(f61)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p534,qd=p61,v=cf728(i7,i2)%e,a=u534_728(i7,i2)%a,b=u534_728(i
* 7,i2)%b,c=u534_728(i7,i2)%c,d=u534_728(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf728(i7,i2)%ek0*(p534(2)*p61(3)-p61(2)*p534(3))+p
     & 534k0*(cf728(i7,i2)%e(2)*p61(3)-p61(2)*cf728(i7,i2)%e(3))
     & -p61k0*(cf728(i7,i2)%e(2)*p534(3)-p534(2)*cf728(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf728(i7,i2)%e(3)*p534k0+p534(3)*cf728(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf728(i7,i2)%e(3)*p61k0+p61(3)*cf728(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf728(i7,i2)%e(0)*p534(0)-cf728(i7,i2)%e(1)*p534(1)-c
     & f728(i7,i2)%e(2)*p534(2)-cf728(i7,i2)%e(3)*p534(3)
      cvqd=cf728(i7,i2)%e(0)*p61(0)-cf728(i7,i2)%e(1)*p61(1)-cf7
     & 28(i7,i2)%e(2)*p61(2)-cf728(i7,i2)%e(3)*p61(3)
      cauxa=-cf728(i7,i2)%ek0*quqd+p534k0*cvqd+p61k0*cvqu
      cauxb=-cf728(i7,i2)%ek0*p61(2)+p61k0*cf728(i7,i2)%e(2)
      cauxc=+cf728(i7,i2)%ek0*p534(2)-p534k0*cf728(i7,i2)%e(2)
      u534_728(i7,i2)%a(1)=u534_728(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u534_728(i7,i2)%a(2)=u534_728(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u534_728(i7,i2)%b(1)=u534_728(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u534_728(i7,i2)%b(2)=u534_728(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u534_728(i7,i2)%c(1)=u534_728(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u534_728(i7,i2)%c(2)=u534_728(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u534_728(i7,i2)%d(1)=u534_728(i7,i2)%d(1)+ccl*cf728(i7,i2)
     & %ek0
      u534_728(i7,i2)%d(2)=u534_728(i7,i2)%d(2)+ccr*cf728(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_34728(i5,i7,i2)%a,cc=l5_34728(i5,i7,i2)%c,a1=l5_34(i5)%a
* ,c1=l5_34(i5)%c,a2=u534_728(i7,i2)%a,b2=u534_728(i7,i2)%b,c2=u534_728(
* i7,i2)%c,d2=u534_728(i7,i2)%d,prq=s534,nsum=0
      l5_34728(i5,i7,i2)%a(1)=l5_34(i5)%a(1)*u534_728(i7,i2)%a(1
     & )+l5_34(i5)%c(1)*s534*u534_728(i7,i2)%b(2)
      l5_34728(i5,i7,i2)%c(1)=l5_34(i5)%a(1)*u534_728(i7,i2)%c(1
     & )+l5_34(i5)%c(1)*s534*u534_728(i7,i2)%d(2)
      l5_34728(i5,i7,i2)%c(2)=l5_34(i5)%c(2)*s534*u534_728(i7,i2
     & )%d(1)+l5_34(i5)%a(2)*u534_728(i7,i2)%c(2)
      l5_34728(i5,i7,i2)%a(2)=l5_34(i5)%c(2)*s534*u534_728(i7,i2
     & )%b(1)+l5_34(i5)%a(2)*u534_728(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u5278_34                                               
* quqd -- p=p5278,q=p61
      quqd=p5278(0)*p61(0)-p5278(1)*p61(1)-p5278(2)*p61(2)-p5278
     & (3)*p61(3)
      ccr=zcr(id5)/(f61)
      ccl=zcl(id5)/(f61)
      do i5=1,2
* T0 -- qu=p5278,qd=p61,v=cz34(i5)%e,a=u5728_34(i5)%a,b=u5728_34(i5)%b,c
* =u5728_34(i5)%c,d=u5728_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p5278(2)*p61(3)-p61(2)*p5278(3))+p52
     & 78k0*(cz34(i5)%e(2)*p61(3)-p61(2)*cz34(i5)%e(3))-p61k0*(c
     & z34(i5)%e(2)*p5278(3)-p5278(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p5278k0+p5278(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p61k0+p61(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p5278(0)-cz34(i5)%e(1)*p5278(1)-cz34(i5
     & )%e(2)*p5278(2)-cz34(i5)%e(3)*p5278(3)
      cvqd=cz34(i5)%e(0)*p61(0)-cz34(i5)%e(1)*p61(1)-cz34(i5)%e(
     & 2)*p61(2)-cz34(i5)%e(3)*p61(3)
      cauxa=-cz34(i5)%ek0*quqd+p5278k0*cvqd+p61k0*cvqu
      cauxb=-cz34(i5)%ek0*p61(2)+p61k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p5278(2)-p5278k0*cz34(i5)%e(2)
      u5728_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u5728_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u5728_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u5728_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u5728_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u5728_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u5728_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u5728_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f61)
      ccl=fcl(id5)/(f61)
      do i5=1,2
* T0 -- qu=p5278,qd=p61,v=cf34(i5)%e,a=u5728_34(i5)%a,b=u5728_34(i5)%b,c
* =u5728_34(i5)%c,d=u5728_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p5278(2)*p61(3)-p61(2)*p5278(3))+p52
     & 78k0*(cf34(i5)%e(2)*p61(3)-p61(2)*cf34(i5)%e(3))-p61k0*(c
     & f34(i5)%e(2)*p5278(3)-p5278(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p5278k0+p5278(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p61k0+p61(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p5278(0)-cf34(i5)%e(1)*p5278(1)-cf34(i5
     & )%e(2)*p5278(2)-cf34(i5)%e(3)*p5278(3)
      cvqd=cf34(i5)%e(0)*p61(0)-cf34(i5)%e(1)*p61(1)-cf34(i5)%e(
     & 2)*p61(2)-cf34(i5)%e(3)*p61(3)
      cauxa=-cf34(i5)%ek0*quqd+p5278k0*cvqd+p61k0*cvqu
      cauxb=-cf34(i5)%ek0*p61(2)+p61k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p5278(2)-p5278k0*cf34(i5)%e(2)
      u5728_34(i5)%a(1)=u5728_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u5728_34(i5)%a(2)=u5728_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u5728_34(i5)%b(1)=u5728_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u5728_34(i5)%b(2)=u5728_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u5728_34(i5)%c(1)=u5728_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u5728_34(i5)%c(2)=u5728_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u5728_34(i5)%d(1)=u5728_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u5728_34(i5)%d(2)=u5728_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_34728(i5,i7,i2)%a,cc=l5_34728(i5,i7,i2)%c,a1=l5_728(i7,i
* 2)%a,c1=l5_728(i7,i2)%c,a2=u5728_34(i5)%a,b2=u5728_34(i5)%b,c2=u5728_3
* 4(i5)%c,d2=u5728_34(i5)%d,prq=s5278,nsum=1
      l5_34728(i5,i7,i2)%a(1)=l5_34728(i5,i7,i2)%a(1)+l5_728(i7,
     & i2)%a(1)*u5728_34(i5)%a(1)+l5_728(i7,i2)%c(1)*s5278*u5728
     & _34(i5)%b(2)
      l5_34728(i5,i7,i2)%c(1)=l5_34728(i5,i7,i2)%c(1)+l5_728(i7,
     & i2)%a(1)*u5728_34(i5)%c(1)+l5_728(i7,i2)%c(1)*s5278*u5728
     & _34(i5)%d(2)
      l5_34728(i5,i7,i2)%c(2)=l5_34728(i5,i7,i2)%c(2)+l5_728(i7,
     & i2)%c(2)*s5278*u5728_34(i5)%d(1)+l5_728(i7,i2)%a(2)*u5728
     & _34(i5)%c(2)
      l5_34728(i5,i7,i2)%a(2)=l5_34728(i5,i7,i2)%a(2)+l5_728(i7,
     & i2)%c(2)*s5278*u5728_34(i5)%b(1)+l5_728(i7,i2)%a(2)*u5728
     & _34(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p52,q=p5234
      quqd=p52(0)*p5234(0)-p52(1)*p5234(1)-p52(2)*p5234(2)-p52(3
     & )*p5234(3)
      ccr=zcr(id5)/(f5234)
      ccl=zcl(id5)/(f5234)
      do i5=1,2
* T0 -- qu=p52,qd=p5234,v=cz34(i5)%e,a=u52_34(i5)%a,b=u52_34(i5)%b,c=u52
* _34(i5)%c,d=u52_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p52(2)*p5234(3)-p5234(2)*p52(3))+p52
     & k0*(cz34(i5)%e(2)*p5234(3)-p5234(2)*cz34(i5)%e(3))-p5234k
     & 0*(cz34(i5)%e(2)*p52(3)-p52(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p52k0+p52(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p5234k0+p5234(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p52(0)-cz34(i5)%e(1)*p52(1)-cz34(i5)%e(
     & 2)*p52(2)-cz34(i5)%e(3)*p52(3)
      cvqd=cz34(i5)%e(0)*p5234(0)-cz34(i5)%e(1)*p5234(1)-cz34(i5
     & )%e(2)*p5234(2)-cz34(i5)%e(3)*p5234(3)
      cauxa=-cz34(i5)%ek0*quqd+p52k0*cvqd+p5234k0*cvqu
      cauxb=-cz34(i5)%ek0*p5234(2)+p5234k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p52(2)-p52k0*cz34(i5)%e(2)
      u52_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u52_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u52_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u52_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u52_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u52_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u52_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u52_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f5234)
      ccl=fcl(id5)/(f5234)
      do i5=1,2
* T0 -- qu=p52,qd=p5234,v=cf34(i5)%e,a=u52_34(i5)%a,b=u52_34(i5)%b,c=u52
* _34(i5)%c,d=u52_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p52(2)*p5234(3)-p5234(2)*p52(3))+p52
     & k0*(cf34(i5)%e(2)*p5234(3)-p5234(2)*cf34(i5)%e(3))-p5234k
     & 0*(cf34(i5)%e(2)*p52(3)-p52(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p52k0+p52(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p5234k0+p5234(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p52(0)-cf34(i5)%e(1)*p52(1)-cf34(i5)%e(
     & 2)*p52(2)-cf34(i5)%e(3)*p52(3)
      cvqd=cf34(i5)%e(0)*p5234(0)-cf34(i5)%e(1)*p5234(1)-cf34(i5
     & )%e(2)*p5234(2)-cf34(i5)%e(3)*p5234(3)
      cauxa=-cf34(i5)%ek0*quqd+p52k0*cvqd+p5234k0*cvqu
      cauxb=-cf34(i5)%ek0*p5234(2)+p5234k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p52(2)-p52k0*cf34(i5)%e(2)
      u52_34(i5)%a(1)=u52_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u52_34(i5)%a(2)=u52_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u52_34(i5)%b(1)=u52_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u52_34(i5)%b(2)=u52_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u52_34(i5)%c(1)=u52_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u52_34(i5)%c(2)=u52_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u52_34(i5)%d(1)=u52_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u52_34(i5)%d(2)=u52_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l5_234(i1,i5)%a,cc=l5_234(i1,i5)%c,a1=l5_2(i1)%a,c1=l5_2(i1
* )%c,a2=u52_34(i5)%a,b2=u52_34(i5)%b,c2=u52_34(i5)%c,d2=u52_34(i5)%d,pr
* q=s52,nsum=0
      l5_234(i1,i5)%a(1)=l5_2(i1)%a(1)*u52_34(i5)%a(1)+l5_2(i1)%
     & c(1)*s52*u52_34(i5)%b(2)
      l5_234(i1,i5)%c(1)=l5_2(i1)%a(1)*u52_34(i5)%c(1)+l5_2(i1)%
     & c(1)*s52*u52_34(i5)%d(2)
      l5_234(i1,i5)%c(2)=l5_2(i1)%c(2)*s52*u52_34(i5)%d(1)+l5_2(
     & i1)%a(2)*u52_34(i5)%c(2)
      l5_234(i1,i5)%a(2)=l5_2(i1)%c(2)*s52*u52_34(i5)%b(1)+l5_2(
     & i1)%a(2)*u52_34(i5)%a(2)
      end do
      end do
  
* quqd -- p=p534,q=p5234
      quqd=p534(0)*p5234(0)-p534(1)*p5234(1)-p534(2)*p5234(2)-p5
     & 34(3)*p5234(3)
      ccr=1.d0/(f5234)
      ccl=1.d0/(f5234)
      do i1=1,2
* T0 -- qu=p534,qd=p5234,v=ce2(i1)%e,a=u534_2(i1)%a,b=u534_2(i1)%b,c=u53
* 4_2(i1)%c,d=u534_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p534(2)*p5234(3)-p5234(2)*p534(3))+p5
     & 34k0*(ce2(i1)%e(2)*p5234(3)-p5234(2)*ce2(i1)%e(3))-p5234k
     & 0*(ce2(i1)%e(2)*p534(3)-p534(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p534k0+p534(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p5234k0+p5234(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p534(0)-ce2(i1)%e(1)*p534(1)-ce2(i1)%e(2
     & )*p534(2)-ce2(i1)%e(3)*p534(3)
      cvqd=ce2(i1)%e(0)*p5234(0)-ce2(i1)%e(1)*p5234(1)-ce2(i1)%e
     & (2)*p5234(2)-ce2(i1)%e(3)*p5234(3)
      cauxa=-ce2(i1)%ek0*quqd+p534k0*cvqd+p5234k0*cvqu
      cauxb=-ce2(i1)%ek0*p5234(2)+p5234k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p534(2)-p534k0*ce2(i1)%e(2)
      u534_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u534_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u534_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u534_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u534_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u534_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u534_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u534_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l5_234(i1,i5)%a,cc=l5_234(i1,i5)%c,a1=l5_34(i5)%a,c1=l5_34(
* i5)%c,a2=u534_2(i1)%a,b2=u534_2(i1)%b,c2=u534_2(i1)%c,d2=u534_2(i1)%d,
* prq=s534,nsum=1
      l5_234(i1,i5)%a(1)=l5_234(i1,i5)%a(1)+l5_34(i5)%a(1)*u534_
     & 2(i1)%a(1)+l5_34(i5)%c(1)*s534*u534_2(i1)%b(2)
      l5_234(i1,i5)%c(1)=l5_234(i1,i5)%c(1)+l5_34(i5)%a(1)*u534_
     & 2(i1)%c(1)+l5_34(i5)%c(1)*s534*u534_2(i1)%d(2)
      l5_234(i1,i5)%c(2)=l5_234(i1,i5)%c(2)+l5_34(i5)%c(2)*s534*
     & u534_2(i1)%d(1)+l5_34(i5)%a(2)*u534_2(i1)%c(2)
      l5_234(i1,i5)%a(2)=l5_234(i1,i5)%a(2)+l5_34(i5)%c(2)*s534*
     & u534_2(i1)%b(1)+l5_34(i5)%a(2)*u534_2(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id7).ne.1) then
* quqd -- p=p52,q=p634
      quqd=p52(0)*p634(0)-p52(1)*p634(1)-p52(2)*p634(2)-p52(3)*p
     & 634(3)
      ccr=zcr(id5)/(f634)
      ccl=zcl(id5)/(f634)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p52,qd=p634,v=cz718(i7,i2)%e,a=u52_718(i7,i2)%a,b=u52_718(i7,
* i2)%b,c=u52_718(i7,i2)%c,d=u52_718(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz718(i7,i2)%ek0*(p52(2)*p634(3)-p634(2)*p52(3))+p
     & 52k0*(cz718(i7,i2)%e(2)*p634(3)-p634(2)*cz718(i7,i2)%e(3)
     & )-p634k0*(cz718(i7,i2)%e(2)*p52(3)-p52(2)*cz718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz718(i7,i2)%e(3)*p52k0+p52(3)*cz718(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz718(i7,i2)%e(3)*p634k0+p634(3)*cz718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz718(i7,i2)%e(0)*p52(0)-cz718(i7,i2)%e(1)*p52(1)-cz7
     & 18(i7,i2)%e(2)*p52(2)-cz718(i7,i2)%e(3)*p52(3)
      cvqd=cz718(i7,i2)%e(0)*p634(0)-cz718(i7,i2)%e(1)*p634(1)-c
     & z718(i7,i2)%e(2)*p634(2)-cz718(i7,i2)%e(3)*p634(3)
      cauxa=-cz718(i7,i2)%ek0*quqd+p52k0*cvqd+p634k0*cvqu
      cauxb=-cz718(i7,i2)%ek0*p634(2)+p634k0*cz718(i7,i2)%e(2)
      cauxc=+cz718(i7,i2)%ek0*p52(2)-p52k0*cz718(i7,i2)%e(2)
      u52_718(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u52_718(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u52_718(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u52_718(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u52_718(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u52_718(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u52_718(i7,i2)%d(1)=ccl*cz718(i7,i2)%ek0
      u52_718(i7,i2)%d(2)=ccr*cz718(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id5)/(f634)
      ccl=fcl(id5)/(f634)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p52,qd=p634,v=cf718(i7,i2)%e,a=u52_718(i7,i2)%a,b=u52_718(i7,
* i2)%b,c=u52_718(i7,i2)%c,d=u52_718(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf718(i7,i2)%ek0*(p52(2)*p634(3)-p634(2)*p52(3))+p
     & 52k0*(cf718(i7,i2)%e(2)*p634(3)-p634(2)*cf718(i7,i2)%e(3)
     & )-p634k0*(cf718(i7,i2)%e(2)*p52(3)-p52(2)*cf718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf718(i7,i2)%e(3)*p52k0+p52(3)*cf718(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf718(i7,i2)%e(3)*p634k0+p634(3)*cf718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf718(i7,i2)%e(0)*p52(0)-cf718(i7,i2)%e(1)*p52(1)-cf7
     & 18(i7,i2)%e(2)*p52(2)-cf718(i7,i2)%e(3)*p52(3)
      cvqd=cf718(i7,i2)%e(0)*p634(0)-cf718(i7,i2)%e(1)*p634(1)-c
     & f718(i7,i2)%e(2)*p634(2)-cf718(i7,i2)%e(3)*p634(3)
      cauxa=-cf718(i7,i2)%ek0*quqd+p52k0*cvqd+p634k0*cvqu
      cauxb=-cf718(i7,i2)%ek0*p634(2)+p634k0*cf718(i7,i2)%e(2)
      cauxc=+cf718(i7,i2)%ek0*p52(2)-p52k0*cf718(i7,i2)%e(2)
      u52_718(i7,i2)%a(1)=u52_718(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u52_718(i7,i2)%a(2)=u52_718(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u52_718(i7,i2)%b(1)=u52_718(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u52_718(i7,i2)%b(2)=u52_718(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u52_718(i7,i2)%c(1)=u52_718(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u52_718(i7,i2)%c(2)=u52_718(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u52_718(i7,i2)%d(1)=u52_718(i7,i2)%d(1)+ccl*cf718(i7,i2)%e
     & k0
      u52_718(i7,i2)%d(2)=u52_718(i7,i2)%d(2)+ccr*cf718(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_2718(i1,i7,i2)%a,cc=l5_2718(i1,i7,i2)%c,a1=l5_2(i1)%a,c1
* =l5_2(i1)%c,a2=u52_718(i7,i2)%a,b2=u52_718(i7,i2)%b,c2=u52_718(i7,i2)%
* c,d2=u52_718(i7,i2)%d,prq=s52,nsum=0
      l5_2718(i1,i7,i2)%a(1)=l5_2(i1)%a(1)*u52_718(i7,i2)%a(1)+l
     & 5_2(i1)%c(1)*s52*u52_718(i7,i2)%b(2)
      l5_2718(i1,i7,i2)%c(1)=l5_2(i1)%a(1)*u52_718(i7,i2)%c(1)+l
     & 5_2(i1)%c(1)*s52*u52_718(i7,i2)%d(2)
      l5_2718(i1,i7,i2)%c(2)=l5_2(i1)%c(2)*s52*u52_718(i7,i2)%d(
     & 1)+l5_2(i1)%a(2)*u52_718(i7,i2)%c(2)
      l5_2718(i1,i7,i2)%a(2)=l5_2(i1)%c(2)*s52*u52_718(i7,i2)%b(
     & 1)+l5_2(i1)%a(2)*u52_718(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u5178_2                                                
* quqd -- p=p5178,q=p634
      quqd=p5178(0)*p634(0)-p5178(1)*p634(1)-p5178(2)*p634(2)-p5
     & 178(3)*p634(3)
      ccr=1.d0/(f634)
      ccl=1.d0/(f634)
      do i1=1,2
* T0 -- qu=p5178,qd=p634,v=ce2(i1)%e,a=u5718_2(i1)%a,b=u5718_2(i1)%b,c=u
* 5718_2(i1)%c,d=u5718_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p5178(2)*p634(3)-p634(2)*p5178(3))+p5
     & 178k0*(ce2(i1)%e(2)*p634(3)-p634(2)*ce2(i1)%e(3))-p634k0*
     & (ce2(i1)%e(2)*p5178(3)-p5178(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p5178k0+p5178(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p634k0+p634(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p5178(0)-ce2(i1)%e(1)*p5178(1)-ce2(i1)%e
     & (2)*p5178(2)-ce2(i1)%e(3)*p5178(3)
      cvqd=ce2(i1)%e(0)*p634(0)-ce2(i1)%e(1)*p634(1)-ce2(i1)%e(2
     & )*p634(2)-ce2(i1)%e(3)*p634(3)
      cauxa=-ce2(i1)%ek0*quqd+p5178k0*cvqd+p634k0*cvqu
      cauxb=-ce2(i1)%ek0*p634(2)+p634k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p5178(2)-p5178k0*ce2(i1)%e(2)
      u5718_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u5718_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u5718_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u5718_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u5718_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u5718_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u5718_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u5718_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_2718(i1,i7,i2)%a,cc=l5_2718(i1,i7,i2)%c,a1=l5_718(i7,i2)
* %a,c1=l5_718(i7,i2)%c,a2=u5718_2(i1)%a,b2=u5718_2(i1)%b,c2=u5718_2(i1)
* %c,d2=u5718_2(i1)%d,prq=s5178,nsum=1
      l5_2718(i1,i7,i2)%a(1)=l5_2718(i1,i7,i2)%a(1)+l5_718(i7,i2
     & )%a(1)*u5718_2(i1)%a(1)+l5_718(i7,i2)%c(1)*s5178*u5718_2(
     & i1)%b(2)
      l5_2718(i1,i7,i2)%c(1)=l5_2718(i1,i7,i2)%c(1)+l5_718(i7,i2
     & )%a(1)*u5718_2(i1)%c(1)+l5_718(i7,i2)%c(1)*s5178*u5718_2(
     & i1)%d(2)
      l5_2718(i1,i7,i2)%c(2)=l5_2718(i1,i7,i2)%c(2)+l5_718(i7,i2
     & )%c(2)*s5178*u5718_2(i1)%d(1)+l5_718(i7,i2)%a(2)*u5718_2(
     & i1)%c(2)
      l5_2718(i1,i7,i2)%a(2)=l5_2718(i1,i7,i2)%a(2)+l5_718(i7,i2
     & )%c(2)*s5178*u5718_2(i1)%b(1)+l5_718(i7,i2)%a(2)*u5718_2(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p534,q=p62
      quqd=p534(0)*p62(0)-p534(1)*p62(1)-p534(2)*p62(2)-p534(3)*
     & p62(3)
      ccr=zcr(id5)/(f62)
      ccl=zcl(id5)/(f62)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p534,qd=p62,v=cz718(i7,i2)%e,a=u534_718(i7,i2)%a,b=u534_718(i
* 7,i2)%b,c=u534_718(i7,i2)%c,d=u534_718(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz718(i7,i2)%ek0*(p534(2)*p62(3)-p62(2)*p534(3))+p
     & 534k0*(cz718(i7,i2)%e(2)*p62(3)-p62(2)*cz718(i7,i2)%e(3))
     & -p62k0*(cz718(i7,i2)%e(2)*p534(3)-p534(2)*cz718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz718(i7,i2)%e(3)*p534k0+p534(3)*cz718(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz718(i7,i2)%e(3)*p62k0+p62(3)*cz718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz718(i7,i2)%e(0)*p534(0)-cz718(i7,i2)%e(1)*p534(1)-c
     & z718(i7,i2)%e(2)*p534(2)-cz718(i7,i2)%e(3)*p534(3)
      cvqd=cz718(i7,i2)%e(0)*p62(0)-cz718(i7,i2)%e(1)*p62(1)-cz7
     & 18(i7,i2)%e(2)*p62(2)-cz718(i7,i2)%e(3)*p62(3)
      cauxa=-cz718(i7,i2)%ek0*quqd+p534k0*cvqd+p62k0*cvqu
      cauxb=-cz718(i7,i2)%ek0*p62(2)+p62k0*cz718(i7,i2)%e(2)
      cauxc=+cz718(i7,i2)%ek0*p534(2)-p534k0*cz718(i7,i2)%e(2)
      u534_718(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u534_718(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u534_718(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u534_718(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u534_718(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u534_718(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u534_718(i7,i2)%d(1)=ccl*cz718(i7,i2)%ek0
      u534_718(i7,i2)%d(2)=ccr*cz718(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id5)/(f62)
      ccl=fcl(id5)/(f62)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p534,qd=p62,v=cf718(i7,i2)%e,a=u534_718(i7,i2)%a,b=u534_718(i
* 7,i2)%b,c=u534_718(i7,i2)%c,d=u534_718(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf718(i7,i2)%ek0*(p534(2)*p62(3)-p62(2)*p534(3))+p
     & 534k0*(cf718(i7,i2)%e(2)*p62(3)-p62(2)*cf718(i7,i2)%e(3))
     & -p62k0*(cf718(i7,i2)%e(2)*p534(3)-p534(2)*cf718(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf718(i7,i2)%e(3)*p534k0+p534(3)*cf718(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf718(i7,i2)%e(3)*p62k0+p62(3)*cf718(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf718(i7,i2)%e(0)*p534(0)-cf718(i7,i2)%e(1)*p534(1)-c
     & f718(i7,i2)%e(2)*p534(2)-cf718(i7,i2)%e(3)*p534(3)
      cvqd=cf718(i7,i2)%e(0)*p62(0)-cf718(i7,i2)%e(1)*p62(1)-cf7
     & 18(i7,i2)%e(2)*p62(2)-cf718(i7,i2)%e(3)*p62(3)
      cauxa=-cf718(i7,i2)%ek0*quqd+p534k0*cvqd+p62k0*cvqu
      cauxb=-cf718(i7,i2)%ek0*p62(2)+p62k0*cf718(i7,i2)%e(2)
      cauxc=+cf718(i7,i2)%ek0*p534(2)-p534k0*cf718(i7,i2)%e(2)
      u534_718(i7,i2)%a(1)=u534_718(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u534_718(i7,i2)%a(2)=u534_718(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u534_718(i7,i2)%b(1)=u534_718(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u534_718(i7,i2)%b(2)=u534_718(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u534_718(i7,i2)%c(1)=u534_718(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u534_718(i7,i2)%c(2)=u534_718(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u534_718(i7,i2)%d(1)=u534_718(i7,i2)%d(1)+ccl*cf718(i7,i2)
     & %ek0
      u534_718(i7,i2)%d(2)=u534_718(i7,i2)%d(2)+ccr*cf718(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_34718(i5,i7,i2)%a,cc=l5_34718(i5,i7,i2)%c,a1=l5_34(i5)%a
* ,c1=l5_34(i5)%c,a2=u534_718(i7,i2)%a,b2=u534_718(i7,i2)%b,c2=u534_718(
* i7,i2)%c,d2=u534_718(i7,i2)%d,prq=s534,nsum=0
      l5_34718(i5,i7,i2)%a(1)=l5_34(i5)%a(1)*u534_718(i7,i2)%a(1
     & )+l5_34(i5)%c(1)*s534*u534_718(i7,i2)%b(2)
      l5_34718(i5,i7,i2)%c(1)=l5_34(i5)%a(1)*u534_718(i7,i2)%c(1
     & )+l5_34(i5)%c(1)*s534*u534_718(i7,i2)%d(2)
      l5_34718(i5,i7,i2)%c(2)=l5_34(i5)%c(2)*s534*u534_718(i7,i2
     & )%d(1)+l5_34(i5)%a(2)*u534_718(i7,i2)%c(2)
      l5_34718(i5,i7,i2)%a(2)=l5_34(i5)%c(2)*s534*u534_718(i7,i2
     & )%b(1)+l5_34(i5)%a(2)*u534_718(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u5178_34                                               
* quqd -- p=p5178,q=p62
      quqd=p5178(0)*p62(0)-p5178(1)*p62(1)-p5178(2)*p62(2)-p5178
     & (3)*p62(3)
      ccr=zcr(id5)/(f62)
      ccl=zcl(id5)/(f62)
      do i5=1,2
* T0 -- qu=p5178,qd=p62,v=cz34(i5)%e,a=u5718_34(i5)%a,b=u5718_34(i5)%b,c
* =u5718_34(i5)%c,d=u5718_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p5178(2)*p62(3)-p62(2)*p5178(3))+p51
     & 78k0*(cz34(i5)%e(2)*p62(3)-p62(2)*cz34(i5)%e(3))-p62k0*(c
     & z34(i5)%e(2)*p5178(3)-p5178(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p5178k0+p5178(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p62k0+p62(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p5178(0)-cz34(i5)%e(1)*p5178(1)-cz34(i5
     & )%e(2)*p5178(2)-cz34(i5)%e(3)*p5178(3)
      cvqd=cz34(i5)%e(0)*p62(0)-cz34(i5)%e(1)*p62(1)-cz34(i5)%e(
     & 2)*p62(2)-cz34(i5)%e(3)*p62(3)
      cauxa=-cz34(i5)%ek0*quqd+p5178k0*cvqd+p62k0*cvqu
      cauxb=-cz34(i5)%ek0*p62(2)+p62k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p5178(2)-p5178k0*cz34(i5)%e(2)
      u5718_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u5718_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u5718_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u5718_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u5718_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u5718_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u5718_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u5718_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f62)
      ccl=fcl(id5)/(f62)
      do i5=1,2
* T0 -- qu=p5178,qd=p62,v=cf34(i5)%e,a=u5718_34(i5)%a,b=u5718_34(i5)%b,c
* =u5718_34(i5)%c,d=u5718_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p5178(2)*p62(3)-p62(2)*p5178(3))+p51
     & 78k0*(cf34(i5)%e(2)*p62(3)-p62(2)*cf34(i5)%e(3))-p62k0*(c
     & f34(i5)%e(2)*p5178(3)-p5178(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p5178k0+p5178(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p62k0+p62(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p5178(0)-cf34(i5)%e(1)*p5178(1)-cf34(i5
     & )%e(2)*p5178(2)-cf34(i5)%e(3)*p5178(3)
      cvqd=cf34(i5)%e(0)*p62(0)-cf34(i5)%e(1)*p62(1)-cf34(i5)%e(
     & 2)*p62(2)-cf34(i5)%e(3)*p62(3)
      cauxa=-cf34(i5)%ek0*quqd+p5178k0*cvqd+p62k0*cvqu
      cauxb=-cf34(i5)%ek0*p62(2)+p62k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p5178(2)-p5178k0*cf34(i5)%e(2)
      u5718_34(i5)%a(1)=u5718_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u5718_34(i5)%a(2)=u5718_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u5718_34(i5)%b(1)=u5718_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u5718_34(i5)%b(2)=u5718_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u5718_34(i5)%c(1)=u5718_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u5718_34(i5)%c(2)=u5718_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u5718_34(i5)%d(1)=u5718_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u5718_34(i5)%d(2)=u5718_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_34718(i5,i7,i2)%a,cc=l5_34718(i5,i7,i2)%c,a1=l5_718(i7,i
* 2)%a,c1=l5_718(i7,i2)%c,a2=u5718_34(i5)%a,b2=u5718_34(i5)%b,c2=u5718_3
* 4(i5)%c,d2=u5718_34(i5)%d,prq=s5178,nsum=1
      l5_34718(i5,i7,i2)%a(1)=l5_34718(i5,i7,i2)%a(1)+l5_718(i7,
     & i2)%a(1)*u5718_34(i5)%a(1)+l5_718(i7,i2)%c(1)*s5178*u5718
     & _34(i5)%b(2)
      l5_34718(i5,i7,i2)%c(1)=l5_34718(i5,i7,i2)%c(1)+l5_718(i7,
     & i2)%a(1)*u5718_34(i5)%c(1)+l5_718(i7,i2)%c(1)*s5178*u5718
     & _34(i5)%d(2)
      l5_34718(i5,i7,i2)%c(2)=l5_34718(i5,i7,i2)%c(2)+l5_718(i7,
     & i2)%c(2)*s5178*u5718_34(i5)%d(1)+l5_718(i7,i2)%a(2)*u5718
     & _34(i5)%c(2)
      l5_34718(i5,i7,i2)%a(2)=l5_34718(i5,i7,i2)%a(2)+l5_718(i7,
     & i2)%c(2)*s5178*u5718_34(i5)%b(1)+l5_718(i7,i2)%a(2)*u5718
     & _34(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p51,q=p5178
      quqd=p51(0)*p5178(0)-p51(1)*p5178(1)-p51(2)*p5178(2)-p51(3
     & )*p5178(3)
      ccr=zcr(id5)/(f5178)
      ccl=zcl(id5)/(f5178)
      do i5=1,2
* T0 -- qu=p51,qd=p5178,v=cz78(i5)%e,a=u51_78(i5)%a,b=u51_78(i5)%b,c=u51
* _78(i5)%c,d=u51_78(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p51(2)*p5178(3)-p5178(2)*p51(3))+p51
     & k0*(cz78(i5)%e(2)*p5178(3)-p5178(2)*cz78(i5)%e(3))-p5178k
     & 0*(cz78(i5)%e(2)*p51(3)-p51(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p51k0+p51(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i5)%e(3)*p5178k0+p5178(3)*cz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i5)%e(0)*p51(0)-cz78(i5)%e(1)*p51(1)-cz78(i5)%e(
     & 2)*p51(2)-cz78(i5)%e(3)*p51(3)
      cvqd=cz78(i5)%e(0)*p5178(0)-cz78(i5)%e(1)*p5178(1)-cz78(i5
     & )%e(2)*p5178(2)-cz78(i5)%e(3)*p5178(3)
      cauxa=-cz78(i5)%ek0*quqd+p51k0*cvqd+p5178k0*cvqu
      cauxb=-cz78(i5)%ek0*p5178(2)+p5178k0*cz78(i5)%e(2)
      cauxc=+cz78(i5)%ek0*p51(2)-p51k0*cz78(i5)%e(2)
      u51_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      u51_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      u51_78(i5)%b(1)=ccl*(cauxb-ceps_2)
      u51_78(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u51_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      u51_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u51_78(i5)%d(1)=ccl*cz78(i5)%ek0
      u51_78(i5)%d(2)=ccr*cz78(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f5178)
      ccl=fcl(id5)/(f5178)
      do i5=1,2
* T0 -- qu=p51,qd=p5178,v=cf78(i5)%e,a=u51_78(i5)%a,b=u51_78(i5)%b,c=u51
* _78(i5)%c,d=u51_78(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p51(2)*p5178(3)-p5178(2)*p51(3))+p51
     & k0*(cf78(i5)%e(2)*p5178(3)-p5178(2)*cf78(i5)%e(3))-p5178k
     & 0*(cf78(i5)%e(2)*p51(3)-p51(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p51k0+p51(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i5)%e(3)*p5178k0+p5178(3)*cf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i5)%e(0)*p51(0)-cf78(i5)%e(1)*p51(1)-cf78(i5)%e(
     & 2)*p51(2)-cf78(i5)%e(3)*p51(3)
      cvqd=cf78(i5)%e(0)*p5178(0)-cf78(i5)%e(1)*p5178(1)-cf78(i5
     & )%e(2)*p5178(2)-cf78(i5)%e(3)*p5178(3)
      cauxa=-cf78(i5)%ek0*quqd+p51k0*cvqd+p5178k0*cvqu
      cauxb=-cf78(i5)%ek0*p5178(2)+p5178k0*cf78(i5)%e(2)
      cauxc=+cf78(i5)%ek0*p51(2)-p51k0*cf78(i5)%e(2)
      u51_78(i5)%a(1)=u51_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      u51_78(i5)%a(2)=u51_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      u51_78(i5)%b(1)=u51_78(i5)%b(1)+ccl*(cauxb-ceps_2)
      u51_78(i5)%b(2)=u51_78(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u51_78(i5)%c(1)=u51_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      u51_78(i5)%c(2)=u51_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u51_78(i5)%d(1)=u51_78(i5)%d(1)+ccl*cf78(i5)%ek0
      u51_78(i5)%d(2)=u51_78(i5)%d(2)+ccr*cf78(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l5_178(i1,i5)%a,cc=l5_178(i1,i5)%c,a1=l5_1(i1)%a,c1=l5_1(i1
* )%c,a2=u51_78(i5)%a,b2=u51_78(i5)%b,c2=u51_78(i5)%c,d2=u51_78(i5)%d,pr
* q=s51,nsum=0
      l5_178(i1,i5)%a(1)=l5_1(i1)%a(1)*u51_78(i5)%a(1)+l5_1(i1)%
     & c(1)*s51*u51_78(i5)%b(2)
      l5_178(i1,i5)%c(1)=l5_1(i1)%a(1)*u51_78(i5)%c(1)+l5_1(i1)%
     & c(1)*s51*u51_78(i5)%d(2)
      l5_178(i1,i5)%c(2)=l5_1(i1)%c(2)*s51*u51_78(i5)%d(1)+l5_1(
     & i1)%a(2)*u51_78(i5)%c(2)
      l5_178(i1,i5)%a(2)=l5_1(i1)%c(2)*s51*u51_78(i5)%b(1)+l5_1(
     & i1)%a(2)*u51_78(i5)%a(2)
      end do
      end do
  
* quqd -- p=p578,q=p5178
      quqd=p578(0)*p5178(0)-p578(1)*p5178(1)-p578(2)*p5178(2)-p5
     & 78(3)*p5178(3)
      ccr=1.d0/(f5178)
      ccl=1.d0/(f5178)
      do i1=1,2
* T0 -- qu=p578,qd=p5178,v=ce1(i1)%e,a=u578_1(i1)%a,b=u578_1(i1)%b,c=u57
* 8_1(i1)%c,d=u578_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p578(2)*p5178(3)-p5178(2)*p578(3))+p5
     & 78k0*(ce1(i1)%e(2)*p5178(3)-p5178(2)*ce1(i1)%e(3))-p5178k
     & 0*(ce1(i1)%e(2)*p578(3)-p578(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p578k0+p578(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p5178k0+p5178(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p578(0)-ce1(i1)%e(1)*p578(1)-ce1(i1)%e(2
     & )*p578(2)-ce1(i1)%e(3)*p578(3)
      cvqd=ce1(i1)%e(0)*p5178(0)-ce1(i1)%e(1)*p5178(1)-ce1(i1)%e
     & (2)*p5178(2)-ce1(i1)%e(3)*p5178(3)
      cauxa=-ce1(i1)%ek0*quqd+p578k0*cvqd+p5178k0*cvqu
      cauxb=-ce1(i1)%ek0*p5178(2)+p5178k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p578(2)-p578k0*ce1(i1)%e(2)
      u578_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u578_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u578_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u578_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u578_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u578_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u578_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u578_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l5_178(i1,i5)%a,cc=l5_178(i1,i5)%c,a1=l5_78(i5)%a,c1=l5_78(
* i5)%c,a2=u578_1(i1)%a,b2=u578_1(i1)%b,c2=u578_1(i1)%c,d2=u578_1(i1)%d,
* prq=s578,nsum=1
      l5_178(i1,i5)%a(1)=l5_178(i1,i5)%a(1)+l5_78(i5)%a(1)*u578_
     & 1(i1)%a(1)+l5_78(i5)%c(1)*s578*u578_1(i1)%b(2)
      l5_178(i1,i5)%c(1)=l5_178(i1,i5)%c(1)+l5_78(i5)%a(1)*u578_
     & 1(i1)%c(1)+l5_78(i5)%c(1)*s578*u578_1(i1)%d(2)
      l5_178(i1,i5)%c(2)=l5_178(i1,i5)%c(2)+l5_78(i5)%c(2)*s578*
     & u578_1(i1)%d(1)+l5_78(i5)%a(2)*u578_1(i1)%c(2)
      l5_178(i1,i5)%a(2)=l5_178(i1,i5)%a(2)+l5_78(i5)%c(2)*s578*
     & u578_1(i1)%b(1)+l5_78(i5)%a(2)*u578_1(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id3).ne.1) then
* quqd -- p=p51,q=p678
      quqd=p51(0)*p678(0)-p51(1)*p678(1)-p51(2)*p678(2)-p51(3)*p
     & 678(3)
      ccr=zcr(id5)/(f678)
      ccl=zcl(id5)/(f678)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p51,qd=p678,v=cz324(i7,i2)%e,a=u51_324(i7,i2)%a,b=u51_324(i7,
* i2)%b,c=u51_324(i7,i2)%c,d=u51_324(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz324(i7,i2)%ek0*(p51(2)*p678(3)-p678(2)*p51(3))+p
     & 51k0*(cz324(i7,i2)%e(2)*p678(3)-p678(2)*cz324(i7,i2)%e(3)
     & )-p678k0*(cz324(i7,i2)%e(2)*p51(3)-p51(2)*cz324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i7,i2)%e(3)*p51k0+p51(3)*cz324(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz324(i7,i2)%e(3)*p678k0+p678(3)*cz324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i7,i2)%e(0)*p51(0)-cz324(i7,i2)%e(1)*p51(1)-cz3
     & 24(i7,i2)%e(2)*p51(2)-cz324(i7,i2)%e(3)*p51(3)
      cvqd=cz324(i7,i2)%e(0)*p678(0)-cz324(i7,i2)%e(1)*p678(1)-c
     & z324(i7,i2)%e(2)*p678(2)-cz324(i7,i2)%e(3)*p678(3)
      cauxa=-cz324(i7,i2)%ek0*quqd+p51k0*cvqd+p678k0*cvqu
      cauxb=-cz324(i7,i2)%ek0*p678(2)+p678k0*cz324(i7,i2)%e(2)
      cauxc=+cz324(i7,i2)%ek0*p51(2)-p51k0*cz324(i7,i2)%e(2)
      u51_324(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u51_324(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u51_324(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u51_324(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u51_324(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u51_324(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u51_324(i7,i2)%d(1)=ccl*cz324(i7,i2)%ek0
      u51_324(i7,i2)%d(2)=ccr*cz324(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id5)/(f678)
      ccl=fcl(id5)/(f678)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p51,qd=p678,v=cf324(i7,i2)%e,a=u51_324(i7,i2)%a,b=u51_324(i7,
* i2)%b,c=u51_324(i7,i2)%c,d=u51_324(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf324(i7,i2)%ek0*(p51(2)*p678(3)-p678(2)*p51(3))+p
     & 51k0*(cf324(i7,i2)%e(2)*p678(3)-p678(2)*cf324(i7,i2)%e(3)
     & )-p678k0*(cf324(i7,i2)%e(2)*p51(3)-p51(2)*cf324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i7,i2)%e(3)*p51k0+p51(3)*cf324(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf324(i7,i2)%e(3)*p678k0+p678(3)*cf324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i7,i2)%e(0)*p51(0)-cf324(i7,i2)%e(1)*p51(1)-cf3
     & 24(i7,i2)%e(2)*p51(2)-cf324(i7,i2)%e(3)*p51(3)
      cvqd=cf324(i7,i2)%e(0)*p678(0)-cf324(i7,i2)%e(1)*p678(1)-c
     & f324(i7,i2)%e(2)*p678(2)-cf324(i7,i2)%e(3)*p678(3)
      cauxa=-cf324(i7,i2)%ek0*quqd+p51k0*cvqd+p678k0*cvqu
      cauxb=-cf324(i7,i2)%ek0*p678(2)+p678k0*cf324(i7,i2)%e(2)
      cauxc=+cf324(i7,i2)%ek0*p51(2)-p51k0*cf324(i7,i2)%e(2)
      u51_324(i7,i2)%a(1)=u51_324(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u51_324(i7,i2)%a(2)=u51_324(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u51_324(i7,i2)%b(1)=u51_324(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u51_324(i7,i2)%b(2)=u51_324(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u51_324(i7,i2)%c(1)=u51_324(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u51_324(i7,i2)%c(2)=u51_324(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u51_324(i7,i2)%d(1)=u51_324(i7,i2)%d(1)+ccl*cf324(i7,i2)%e
     & k0
      u51_324(i7,i2)%d(2)=u51_324(i7,i2)%d(2)+ccr*cf324(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_1324(i1,i7,i2)%a,cc=l5_1324(i1,i7,i2)%c,a1=l5_1(i1)%a,c1
* =l5_1(i1)%c,a2=u51_324(i7,i2)%a,b2=u51_324(i7,i2)%b,c2=u51_324(i7,i2)%
* c,d2=u51_324(i7,i2)%d,prq=s51,nsum=0
      l5_1324(i1,i7,i2)%a(1)=l5_1(i1)%a(1)*u51_324(i7,i2)%a(1)+l
     & 5_1(i1)%c(1)*s51*u51_324(i7,i2)%b(2)
      l5_1324(i1,i7,i2)%c(1)=l5_1(i1)%a(1)*u51_324(i7,i2)%c(1)+l
     & 5_1(i1)%c(1)*s51*u51_324(i7,i2)%d(2)
      l5_1324(i1,i7,i2)%c(2)=l5_1(i1)%c(2)*s51*u51_324(i7,i2)%d(
     & 1)+l5_1(i1)%a(2)*u51_324(i7,i2)%c(2)
      l5_1324(i1,i7,i2)%a(2)=l5_1(i1)%c(2)*s51*u51_324(i7,i2)%b(
     & 1)+l5_1(i1)%a(2)*u51_324(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u5234_1                                                
* quqd -- p=p5234,q=p678
      quqd=p5234(0)*p678(0)-p5234(1)*p678(1)-p5234(2)*p678(2)-p5
     & 234(3)*p678(3)
      ccr=1.d0/(f678)
      ccl=1.d0/(f678)
      do i1=1,2
* T0 -- qu=p5234,qd=p678,v=ce1(i1)%e,a=u5324_1(i1)%a,b=u5324_1(i1)%b,c=u
* 5324_1(i1)%c,d=u5324_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p5234(2)*p678(3)-p678(2)*p5234(3))+p5
     & 234k0*(ce1(i1)%e(2)*p678(3)-p678(2)*ce1(i1)%e(3))-p678k0*
     & (ce1(i1)%e(2)*p5234(3)-p5234(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p5234k0+p5234(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p678k0+p678(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p5234(0)-ce1(i1)%e(1)*p5234(1)-ce1(i1)%e
     & (2)*p5234(2)-ce1(i1)%e(3)*p5234(3)
      cvqd=ce1(i1)%e(0)*p678(0)-ce1(i1)%e(1)*p678(1)-ce1(i1)%e(2
     & )*p678(2)-ce1(i1)%e(3)*p678(3)
      cauxa=-ce1(i1)%ek0*quqd+p5234k0*cvqd+p678k0*cvqu
      cauxb=-ce1(i1)%ek0*p678(2)+p678k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p5234(2)-p5234k0*ce1(i1)%e(2)
      u5324_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u5324_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u5324_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u5324_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u5324_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u5324_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u5324_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u5324_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_1324(i1,i7,i2)%a,cc=l5_1324(i1,i7,i2)%c,a1=l5_324(i7,i2)
* %a,c1=l5_324(i7,i2)%c,a2=u5324_1(i1)%a,b2=u5324_1(i1)%b,c2=u5324_1(i1)
* %c,d2=u5324_1(i1)%d,prq=s5234,nsum=1
      l5_1324(i1,i7,i2)%a(1)=l5_1324(i1,i7,i2)%a(1)+l5_324(i7,i2
     & )%a(1)*u5324_1(i1)%a(1)+l5_324(i7,i2)%c(1)*s5234*u5324_1(
     & i1)%b(2)
      l5_1324(i1,i7,i2)%c(1)=l5_1324(i1,i7,i2)%c(1)+l5_324(i7,i2
     & )%a(1)*u5324_1(i1)%c(1)+l5_324(i7,i2)%c(1)*s5234*u5324_1(
     & i1)%d(2)
      l5_1324(i1,i7,i2)%c(2)=l5_1324(i1,i7,i2)%c(2)+l5_324(i7,i2
     & )%c(2)*s5234*u5324_1(i1)%d(1)+l5_324(i7,i2)%a(2)*u5324_1(
     & i1)%c(2)
      l5_1324(i1,i7,i2)%a(2)=l5_1324(i1,i7,i2)%a(2)+l5_324(i7,i2
     & )%c(2)*s5234*u5324_1(i1)%b(1)+l5_324(i7,i2)%a(2)*u5324_1(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p578,q=p61
      quqd=p578(0)*p61(0)-p578(1)*p61(1)-p578(2)*p61(2)-p578(3)*
     & p61(3)
      ccr=zcr(id5)/(f61)
      ccl=zcl(id5)/(f61)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p578,qd=p61,v=cz324(i7,i2)%e,a=u578_324(i7,i2)%a,b=u578_324(i
* 7,i2)%b,c=u578_324(i7,i2)%c,d=u578_324(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz324(i7,i2)%ek0*(p578(2)*p61(3)-p61(2)*p578(3))+p
     & 578k0*(cz324(i7,i2)%e(2)*p61(3)-p61(2)*cz324(i7,i2)%e(3))
     & -p61k0*(cz324(i7,i2)%e(2)*p578(3)-p578(2)*cz324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i7,i2)%e(3)*p578k0+p578(3)*cz324(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz324(i7,i2)%e(3)*p61k0+p61(3)*cz324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i7,i2)%e(0)*p578(0)-cz324(i7,i2)%e(1)*p578(1)-c
     & z324(i7,i2)%e(2)*p578(2)-cz324(i7,i2)%e(3)*p578(3)
      cvqd=cz324(i7,i2)%e(0)*p61(0)-cz324(i7,i2)%e(1)*p61(1)-cz3
     & 24(i7,i2)%e(2)*p61(2)-cz324(i7,i2)%e(3)*p61(3)
      cauxa=-cz324(i7,i2)%ek0*quqd+p578k0*cvqd+p61k0*cvqu
      cauxb=-cz324(i7,i2)%ek0*p61(2)+p61k0*cz324(i7,i2)%e(2)
      cauxc=+cz324(i7,i2)%ek0*p578(2)-p578k0*cz324(i7,i2)%e(2)
      u578_324(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u578_324(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u578_324(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u578_324(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u578_324(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u578_324(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u578_324(i7,i2)%d(1)=ccl*cz324(i7,i2)%ek0
      u578_324(i7,i2)%d(2)=ccr*cz324(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id5)/(f61)
      ccl=fcl(id5)/(f61)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p578,qd=p61,v=cf324(i7,i2)%e,a=u578_324(i7,i2)%a,b=u578_324(i
* 7,i2)%b,c=u578_324(i7,i2)%c,d=u578_324(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf324(i7,i2)%ek0*(p578(2)*p61(3)-p61(2)*p578(3))+p
     & 578k0*(cf324(i7,i2)%e(2)*p61(3)-p61(2)*cf324(i7,i2)%e(3))
     & -p61k0*(cf324(i7,i2)%e(2)*p578(3)-p578(2)*cf324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i7,i2)%e(3)*p578k0+p578(3)*cf324(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf324(i7,i2)%e(3)*p61k0+p61(3)*cf324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i7,i2)%e(0)*p578(0)-cf324(i7,i2)%e(1)*p578(1)-c
     & f324(i7,i2)%e(2)*p578(2)-cf324(i7,i2)%e(3)*p578(3)
      cvqd=cf324(i7,i2)%e(0)*p61(0)-cf324(i7,i2)%e(1)*p61(1)-cf3
     & 24(i7,i2)%e(2)*p61(2)-cf324(i7,i2)%e(3)*p61(3)
      cauxa=-cf324(i7,i2)%ek0*quqd+p578k0*cvqd+p61k0*cvqu
      cauxb=-cf324(i7,i2)%ek0*p61(2)+p61k0*cf324(i7,i2)%e(2)
      cauxc=+cf324(i7,i2)%ek0*p578(2)-p578k0*cf324(i7,i2)%e(2)
      u578_324(i7,i2)%a(1)=u578_324(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u578_324(i7,i2)%a(2)=u578_324(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u578_324(i7,i2)%b(1)=u578_324(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u578_324(i7,i2)%b(2)=u578_324(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u578_324(i7,i2)%c(1)=u578_324(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u578_324(i7,i2)%c(2)=u578_324(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u578_324(i7,i2)%d(1)=u578_324(i7,i2)%d(1)+ccl*cf324(i7,i2)
     & %ek0
      u578_324(i7,i2)%d(2)=u578_324(i7,i2)%d(2)+ccr*cf324(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_78324(i5,i7,i2)%a,cc=l5_78324(i5,i7,i2)%c,a1=l5_78(i5)%a
* ,c1=l5_78(i5)%c,a2=u578_324(i7,i2)%a,b2=u578_324(i7,i2)%b,c2=u578_324(
* i7,i2)%c,d2=u578_324(i7,i2)%d,prq=s578,nsum=0
      l5_78324(i5,i7,i2)%a(1)=l5_78(i5)%a(1)*u578_324(i7,i2)%a(1
     & )+l5_78(i5)%c(1)*s578*u578_324(i7,i2)%b(2)
      l5_78324(i5,i7,i2)%c(1)=l5_78(i5)%a(1)*u578_324(i7,i2)%c(1
     & )+l5_78(i5)%c(1)*s578*u578_324(i7,i2)%d(2)
      l5_78324(i5,i7,i2)%c(2)=l5_78(i5)%c(2)*s578*u578_324(i7,i2
     & )%d(1)+l5_78(i5)%a(2)*u578_324(i7,i2)%c(2)
      l5_78324(i5,i7,i2)%a(2)=l5_78(i5)%c(2)*s578*u578_324(i7,i2
     & )%b(1)+l5_78(i5)%a(2)*u578_324(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u5234_78                                               
* quqd -- p=p5234,q=p61
      quqd=p5234(0)*p61(0)-p5234(1)*p61(1)-p5234(2)*p61(2)-p5234
     & (3)*p61(3)
      ccr=zcr(id5)/(f61)
      ccl=zcl(id5)/(f61)
      do i5=1,2
* T0 -- qu=p5234,qd=p61,v=cz78(i5)%e,a=u5324_78(i5)%a,b=u5324_78(i5)%b,c
* =u5324_78(i5)%c,d=u5324_78(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p5234(2)*p61(3)-p61(2)*p5234(3))+p52
     & 34k0*(cz78(i5)%e(2)*p61(3)-p61(2)*cz78(i5)%e(3))-p61k0*(c
     & z78(i5)%e(2)*p5234(3)-p5234(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p5234k0+p5234(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i5)%e(3)*p61k0+p61(3)*cz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i5)%e(0)*p5234(0)-cz78(i5)%e(1)*p5234(1)-cz78(i5
     & )%e(2)*p5234(2)-cz78(i5)%e(3)*p5234(3)
      cvqd=cz78(i5)%e(0)*p61(0)-cz78(i5)%e(1)*p61(1)-cz78(i5)%e(
     & 2)*p61(2)-cz78(i5)%e(3)*p61(3)
      cauxa=-cz78(i5)%ek0*quqd+p5234k0*cvqd+p61k0*cvqu
      cauxb=-cz78(i5)%ek0*p61(2)+p61k0*cz78(i5)%e(2)
      cauxc=+cz78(i5)%ek0*p5234(2)-p5234k0*cz78(i5)%e(2)
      u5324_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      u5324_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      u5324_78(i5)%b(1)=ccl*(cauxb-ceps_2)
      u5324_78(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u5324_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      u5324_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u5324_78(i5)%d(1)=ccl*cz78(i5)%ek0
      u5324_78(i5)%d(2)=ccr*cz78(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f61)
      ccl=fcl(id5)/(f61)
      do i5=1,2
* T0 -- qu=p5234,qd=p61,v=cf78(i5)%e,a=u5324_78(i5)%a,b=u5324_78(i5)%b,c
* =u5324_78(i5)%c,d=u5324_78(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p5234(2)*p61(3)-p61(2)*p5234(3))+p52
     & 34k0*(cf78(i5)%e(2)*p61(3)-p61(2)*cf78(i5)%e(3))-p61k0*(c
     & f78(i5)%e(2)*p5234(3)-p5234(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p5234k0+p5234(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i5)%e(3)*p61k0+p61(3)*cf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i5)%e(0)*p5234(0)-cf78(i5)%e(1)*p5234(1)-cf78(i5
     & )%e(2)*p5234(2)-cf78(i5)%e(3)*p5234(3)
      cvqd=cf78(i5)%e(0)*p61(0)-cf78(i5)%e(1)*p61(1)-cf78(i5)%e(
     & 2)*p61(2)-cf78(i5)%e(3)*p61(3)
      cauxa=-cf78(i5)%ek0*quqd+p5234k0*cvqd+p61k0*cvqu
      cauxb=-cf78(i5)%ek0*p61(2)+p61k0*cf78(i5)%e(2)
      cauxc=+cf78(i5)%ek0*p5234(2)-p5234k0*cf78(i5)%e(2)
      u5324_78(i5)%a(1)=u5324_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      u5324_78(i5)%a(2)=u5324_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      u5324_78(i5)%b(1)=u5324_78(i5)%b(1)+ccl*(cauxb-ceps_2)
      u5324_78(i5)%b(2)=u5324_78(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u5324_78(i5)%c(1)=u5324_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      u5324_78(i5)%c(2)=u5324_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u5324_78(i5)%d(1)=u5324_78(i5)%d(1)+ccl*cf78(i5)%ek0
      u5324_78(i5)%d(2)=u5324_78(i5)%d(2)+ccr*cf78(i5)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_78324(i5,i7,i2)%a,cc=l5_78324(i5,i7,i2)%c,a1=l5_324(i7,i
* 2)%a,c1=l5_324(i7,i2)%c,a2=u5324_78(i5)%a,b2=u5324_78(i5)%b,c2=u5324_7
* 8(i5)%c,d2=u5324_78(i5)%d,prq=s5234,nsum=1
      l5_78324(i5,i7,i2)%a(1)=l5_78324(i5,i7,i2)%a(1)+l5_324(i7,
     & i2)%a(1)*u5324_78(i5)%a(1)+l5_324(i7,i2)%c(1)*s5234*u5324
     & _78(i5)%b(2)
      l5_78324(i5,i7,i2)%c(1)=l5_78324(i5,i7,i2)%c(1)+l5_324(i7,
     & i2)%a(1)*u5324_78(i5)%c(1)+l5_324(i7,i2)%c(1)*s5234*u5324
     & _78(i5)%d(2)
      l5_78324(i5,i7,i2)%c(2)=l5_78324(i5,i7,i2)%c(2)+l5_324(i7,
     & i2)%c(2)*s5234*u5324_78(i5)%d(1)+l5_324(i7,i2)%a(2)*u5324
     & _78(i5)%c(2)
      l5_78324(i5,i7,i2)%a(2)=l5_78324(i5,i7,i2)%a(2)+l5_324(i7,
     & i2)%c(2)*s5234*u5324_78(i5)%b(1)+l5_324(i7,i2)%a(2)*u5324
     & _78(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p52,q=p5278
      quqd=p52(0)*p5278(0)-p52(1)*p5278(1)-p52(2)*p5278(2)-p52(3
     & )*p5278(3)
      ccr=zcr(id5)/(f5278)
      ccl=zcl(id5)/(f5278)
      do i5=1,2
* T0 -- qu=p52,qd=p5278,v=cz78(i5)%e,a=u52_78(i5)%a,b=u52_78(i5)%b,c=u52
* _78(i5)%c,d=u52_78(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p52(2)*p5278(3)-p5278(2)*p52(3))+p52
     & k0*(cz78(i5)%e(2)*p5278(3)-p5278(2)*cz78(i5)%e(3))-p5278k
     & 0*(cz78(i5)%e(2)*p52(3)-p52(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p52k0+p52(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i5)%e(3)*p5278k0+p5278(3)*cz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i5)%e(0)*p52(0)-cz78(i5)%e(1)*p52(1)-cz78(i5)%e(
     & 2)*p52(2)-cz78(i5)%e(3)*p52(3)
      cvqd=cz78(i5)%e(0)*p5278(0)-cz78(i5)%e(1)*p5278(1)-cz78(i5
     & )%e(2)*p5278(2)-cz78(i5)%e(3)*p5278(3)
      cauxa=-cz78(i5)%ek0*quqd+p52k0*cvqd+p5278k0*cvqu
      cauxb=-cz78(i5)%ek0*p5278(2)+p5278k0*cz78(i5)%e(2)
      cauxc=+cz78(i5)%ek0*p52(2)-p52k0*cz78(i5)%e(2)
      u52_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      u52_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      u52_78(i5)%b(1)=ccl*(cauxb-ceps_2)
      u52_78(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u52_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      u52_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u52_78(i5)%d(1)=ccl*cz78(i5)%ek0
      u52_78(i5)%d(2)=ccr*cz78(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f5278)
      ccl=fcl(id5)/(f5278)
      do i5=1,2
* T0 -- qu=p52,qd=p5278,v=cf78(i5)%e,a=u52_78(i5)%a,b=u52_78(i5)%b,c=u52
* _78(i5)%c,d=u52_78(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p52(2)*p5278(3)-p5278(2)*p52(3))+p52
     & k0*(cf78(i5)%e(2)*p5278(3)-p5278(2)*cf78(i5)%e(3))-p5278k
     & 0*(cf78(i5)%e(2)*p52(3)-p52(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p52k0+p52(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i5)%e(3)*p5278k0+p5278(3)*cf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i5)%e(0)*p52(0)-cf78(i5)%e(1)*p52(1)-cf78(i5)%e(
     & 2)*p52(2)-cf78(i5)%e(3)*p52(3)
      cvqd=cf78(i5)%e(0)*p5278(0)-cf78(i5)%e(1)*p5278(1)-cf78(i5
     & )%e(2)*p5278(2)-cf78(i5)%e(3)*p5278(3)
      cauxa=-cf78(i5)%ek0*quqd+p52k0*cvqd+p5278k0*cvqu
      cauxb=-cf78(i5)%ek0*p5278(2)+p5278k0*cf78(i5)%e(2)
      cauxc=+cf78(i5)%ek0*p52(2)-p52k0*cf78(i5)%e(2)
      u52_78(i5)%a(1)=u52_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      u52_78(i5)%a(2)=u52_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      u52_78(i5)%b(1)=u52_78(i5)%b(1)+ccl*(cauxb-ceps_2)
      u52_78(i5)%b(2)=u52_78(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u52_78(i5)%c(1)=u52_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      u52_78(i5)%c(2)=u52_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u52_78(i5)%d(1)=u52_78(i5)%d(1)+ccl*cf78(i5)%ek0
      u52_78(i5)%d(2)=u52_78(i5)%d(2)+ccr*cf78(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l5_278(i1,i5)%a,cc=l5_278(i1,i5)%c,a1=l5_2(i1)%a,c1=l5_2(i1
* )%c,a2=u52_78(i5)%a,b2=u52_78(i5)%b,c2=u52_78(i5)%c,d2=u52_78(i5)%d,pr
* q=s52,nsum=0
      l5_278(i1,i5)%a(1)=l5_2(i1)%a(1)*u52_78(i5)%a(1)+l5_2(i1)%
     & c(1)*s52*u52_78(i5)%b(2)
      l5_278(i1,i5)%c(1)=l5_2(i1)%a(1)*u52_78(i5)%c(1)+l5_2(i1)%
     & c(1)*s52*u52_78(i5)%d(2)
      l5_278(i1,i5)%c(2)=l5_2(i1)%c(2)*s52*u52_78(i5)%d(1)+l5_2(
     & i1)%a(2)*u52_78(i5)%c(2)
      l5_278(i1,i5)%a(2)=l5_2(i1)%c(2)*s52*u52_78(i5)%b(1)+l5_2(
     & i1)%a(2)*u52_78(i5)%a(2)
      end do
      end do
  
* quqd -- p=p578,q=p5278
      quqd=p578(0)*p5278(0)-p578(1)*p5278(1)-p578(2)*p5278(2)-p5
     & 78(3)*p5278(3)
      ccr=1.d0/(f5278)
      ccl=1.d0/(f5278)
      do i1=1,2
* T0 -- qu=p578,qd=p5278,v=ce2(i1)%e,a=u578_2(i1)%a,b=u578_2(i1)%b,c=u57
* 8_2(i1)%c,d=u578_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p578(2)*p5278(3)-p5278(2)*p578(3))+p5
     & 78k0*(ce2(i1)%e(2)*p5278(3)-p5278(2)*ce2(i1)%e(3))-p5278k
     & 0*(ce2(i1)%e(2)*p578(3)-p578(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p578k0+p578(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p5278k0+p5278(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p578(0)-ce2(i1)%e(1)*p578(1)-ce2(i1)%e(2
     & )*p578(2)-ce2(i1)%e(3)*p578(3)
      cvqd=ce2(i1)%e(0)*p5278(0)-ce2(i1)%e(1)*p5278(1)-ce2(i1)%e
     & (2)*p5278(2)-ce2(i1)%e(3)*p5278(3)
      cauxa=-ce2(i1)%ek0*quqd+p578k0*cvqd+p5278k0*cvqu
      cauxb=-ce2(i1)%ek0*p5278(2)+p5278k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p578(2)-p578k0*ce2(i1)%e(2)
      u578_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u578_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u578_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u578_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u578_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u578_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u578_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u578_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l5_278(i1,i5)%a,cc=l5_278(i1,i5)%c,a1=l5_78(i5)%a,c1=l5_78(
* i5)%c,a2=u578_2(i1)%a,b2=u578_2(i1)%b,c2=u578_2(i1)%c,d2=u578_2(i1)%d,
* prq=s578,nsum=1
      l5_278(i1,i5)%a(1)=l5_278(i1,i5)%a(1)+l5_78(i5)%a(1)*u578_
     & 2(i1)%a(1)+l5_78(i5)%c(1)*s578*u578_2(i1)%b(2)
      l5_278(i1,i5)%c(1)=l5_278(i1,i5)%c(1)+l5_78(i5)%a(1)*u578_
     & 2(i1)%c(1)+l5_78(i5)%c(1)*s578*u578_2(i1)%d(2)
      l5_278(i1,i5)%c(2)=l5_278(i1,i5)%c(2)+l5_78(i5)%c(2)*s578*
     & u578_2(i1)%d(1)+l5_78(i5)%a(2)*u578_2(i1)%c(2)
      l5_278(i1,i5)%a(2)=l5_278(i1,i5)%a(2)+l5_78(i5)%c(2)*s578*
     & u578_2(i1)%b(1)+l5_78(i5)%a(2)*u578_2(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id3).ne.1) then
* quqd -- p=p52,q=p678
      quqd=p52(0)*p678(0)-p52(1)*p678(1)-p52(2)*p678(2)-p52(3)*p
     & 678(3)
      ccr=zcr(id5)/(f678)
      ccl=zcl(id5)/(f678)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p52,qd=p678,v=cz314(i7,i2)%e,a=u52_314(i7,i2)%a,b=u52_314(i7,
* i2)%b,c=u52_314(i7,i2)%c,d=u52_314(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz314(i7,i2)%ek0*(p52(2)*p678(3)-p678(2)*p52(3))+p
     & 52k0*(cz314(i7,i2)%e(2)*p678(3)-p678(2)*cz314(i7,i2)%e(3)
     & )-p678k0*(cz314(i7,i2)%e(2)*p52(3)-p52(2)*cz314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i7,i2)%e(3)*p52k0+p52(3)*cz314(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz314(i7,i2)%e(3)*p678k0+p678(3)*cz314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i7,i2)%e(0)*p52(0)-cz314(i7,i2)%e(1)*p52(1)-cz3
     & 14(i7,i2)%e(2)*p52(2)-cz314(i7,i2)%e(3)*p52(3)
      cvqd=cz314(i7,i2)%e(0)*p678(0)-cz314(i7,i2)%e(1)*p678(1)-c
     & z314(i7,i2)%e(2)*p678(2)-cz314(i7,i2)%e(3)*p678(3)
      cauxa=-cz314(i7,i2)%ek0*quqd+p52k0*cvqd+p678k0*cvqu
      cauxb=-cz314(i7,i2)%ek0*p678(2)+p678k0*cz314(i7,i2)%e(2)
      cauxc=+cz314(i7,i2)%ek0*p52(2)-p52k0*cz314(i7,i2)%e(2)
      u52_314(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u52_314(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u52_314(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u52_314(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u52_314(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u52_314(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u52_314(i7,i2)%d(1)=ccl*cz314(i7,i2)%ek0
      u52_314(i7,i2)%d(2)=ccr*cz314(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id5)/(f678)
      ccl=fcl(id5)/(f678)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p52,qd=p678,v=cf314(i7,i2)%e,a=u52_314(i7,i2)%a,b=u52_314(i7,
* i2)%b,c=u52_314(i7,i2)%c,d=u52_314(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf314(i7,i2)%ek0*(p52(2)*p678(3)-p678(2)*p52(3))+p
     & 52k0*(cf314(i7,i2)%e(2)*p678(3)-p678(2)*cf314(i7,i2)%e(3)
     & )-p678k0*(cf314(i7,i2)%e(2)*p52(3)-p52(2)*cf314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i7,i2)%e(3)*p52k0+p52(3)*cf314(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf314(i7,i2)%e(3)*p678k0+p678(3)*cf314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i7,i2)%e(0)*p52(0)-cf314(i7,i2)%e(1)*p52(1)-cf3
     & 14(i7,i2)%e(2)*p52(2)-cf314(i7,i2)%e(3)*p52(3)
      cvqd=cf314(i7,i2)%e(0)*p678(0)-cf314(i7,i2)%e(1)*p678(1)-c
     & f314(i7,i2)%e(2)*p678(2)-cf314(i7,i2)%e(3)*p678(3)
      cauxa=-cf314(i7,i2)%ek0*quqd+p52k0*cvqd+p678k0*cvqu
      cauxb=-cf314(i7,i2)%ek0*p678(2)+p678k0*cf314(i7,i2)%e(2)
      cauxc=+cf314(i7,i2)%ek0*p52(2)-p52k0*cf314(i7,i2)%e(2)
      u52_314(i7,i2)%a(1)=u52_314(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u52_314(i7,i2)%a(2)=u52_314(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u52_314(i7,i2)%b(1)=u52_314(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u52_314(i7,i2)%b(2)=u52_314(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u52_314(i7,i2)%c(1)=u52_314(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u52_314(i7,i2)%c(2)=u52_314(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u52_314(i7,i2)%d(1)=u52_314(i7,i2)%d(1)+ccl*cf314(i7,i2)%e
     & k0
      u52_314(i7,i2)%d(2)=u52_314(i7,i2)%d(2)+ccr*cf314(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_2314(i1,i7,i2)%a,cc=l5_2314(i1,i7,i2)%c,a1=l5_2(i1)%a,c1
* =l5_2(i1)%c,a2=u52_314(i7,i2)%a,b2=u52_314(i7,i2)%b,c2=u52_314(i7,i2)%
* c,d2=u52_314(i7,i2)%d,prq=s52,nsum=0
      l5_2314(i1,i7,i2)%a(1)=l5_2(i1)%a(1)*u52_314(i7,i2)%a(1)+l
     & 5_2(i1)%c(1)*s52*u52_314(i7,i2)%b(2)
      l5_2314(i1,i7,i2)%c(1)=l5_2(i1)%a(1)*u52_314(i7,i2)%c(1)+l
     & 5_2(i1)%c(1)*s52*u52_314(i7,i2)%d(2)
      l5_2314(i1,i7,i2)%c(2)=l5_2(i1)%c(2)*s52*u52_314(i7,i2)%d(
     & 1)+l5_2(i1)%a(2)*u52_314(i7,i2)%c(2)
      l5_2314(i1,i7,i2)%a(2)=l5_2(i1)%c(2)*s52*u52_314(i7,i2)%b(
     & 1)+l5_2(i1)%a(2)*u52_314(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u5134_2                                                
* quqd -- p=p5134,q=p678
      quqd=p5134(0)*p678(0)-p5134(1)*p678(1)-p5134(2)*p678(2)-p5
     & 134(3)*p678(3)
      ccr=1.d0/(f678)
      ccl=1.d0/(f678)
      do i1=1,2
* T0 -- qu=p5134,qd=p678,v=ce2(i1)%e,a=u5314_2(i1)%a,b=u5314_2(i1)%b,c=u
* 5314_2(i1)%c,d=u5314_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p5134(2)*p678(3)-p678(2)*p5134(3))+p5
     & 134k0*(ce2(i1)%e(2)*p678(3)-p678(2)*ce2(i1)%e(3))-p678k0*
     & (ce2(i1)%e(2)*p5134(3)-p5134(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p5134k0+p5134(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p678k0+p678(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p5134(0)-ce2(i1)%e(1)*p5134(1)-ce2(i1)%e
     & (2)*p5134(2)-ce2(i1)%e(3)*p5134(3)
      cvqd=ce2(i1)%e(0)*p678(0)-ce2(i1)%e(1)*p678(1)-ce2(i1)%e(2
     & )*p678(2)-ce2(i1)%e(3)*p678(3)
      cauxa=-ce2(i1)%ek0*quqd+p5134k0*cvqd+p678k0*cvqu
      cauxb=-ce2(i1)%ek0*p678(2)+p678k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p5134(2)-p5134k0*ce2(i1)%e(2)
      u5314_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u5314_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u5314_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u5314_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u5314_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u5314_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u5314_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u5314_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_2314(i1,i7,i2)%a,cc=l5_2314(i1,i7,i2)%c,a1=l5_314(i7,i2)
* %a,c1=l5_314(i7,i2)%c,a2=u5314_2(i1)%a,b2=u5314_2(i1)%b,c2=u5314_2(i1)
* %c,d2=u5314_2(i1)%d,prq=s5134,nsum=1
      l5_2314(i1,i7,i2)%a(1)=l5_2314(i1,i7,i2)%a(1)+l5_314(i7,i2
     & )%a(1)*u5314_2(i1)%a(1)+l5_314(i7,i2)%c(1)*s5134*u5314_2(
     & i1)%b(2)
      l5_2314(i1,i7,i2)%c(1)=l5_2314(i1,i7,i2)%c(1)+l5_314(i7,i2
     & )%a(1)*u5314_2(i1)%c(1)+l5_314(i7,i2)%c(1)*s5134*u5314_2(
     & i1)%d(2)
      l5_2314(i1,i7,i2)%c(2)=l5_2314(i1,i7,i2)%c(2)+l5_314(i7,i2
     & )%c(2)*s5134*u5314_2(i1)%d(1)+l5_314(i7,i2)%a(2)*u5314_2(
     & i1)%c(2)
      l5_2314(i1,i7,i2)%a(2)=l5_2314(i1,i7,i2)%a(2)+l5_314(i7,i2
     & )%c(2)*s5134*u5314_2(i1)%b(1)+l5_314(i7,i2)%a(2)*u5314_2(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p578,q=p62
      quqd=p578(0)*p62(0)-p578(1)*p62(1)-p578(2)*p62(2)-p578(3)*
     & p62(3)
      ccr=zcr(id5)/(f62)
      ccl=zcl(id5)/(f62)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p578,qd=p62,v=cz314(i7,i2)%e,a=u578_314(i7,i2)%a,b=u578_314(i
* 7,i2)%b,c=u578_314(i7,i2)%c,d=u578_314(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz314(i7,i2)%ek0*(p578(2)*p62(3)-p62(2)*p578(3))+p
     & 578k0*(cz314(i7,i2)%e(2)*p62(3)-p62(2)*cz314(i7,i2)%e(3))
     & -p62k0*(cz314(i7,i2)%e(2)*p578(3)-p578(2)*cz314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i7,i2)%e(3)*p578k0+p578(3)*cz314(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz314(i7,i2)%e(3)*p62k0+p62(3)*cz314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i7,i2)%e(0)*p578(0)-cz314(i7,i2)%e(1)*p578(1)-c
     & z314(i7,i2)%e(2)*p578(2)-cz314(i7,i2)%e(3)*p578(3)
      cvqd=cz314(i7,i2)%e(0)*p62(0)-cz314(i7,i2)%e(1)*p62(1)-cz3
     & 14(i7,i2)%e(2)*p62(2)-cz314(i7,i2)%e(3)*p62(3)
      cauxa=-cz314(i7,i2)%ek0*quqd+p578k0*cvqd+p62k0*cvqu
      cauxb=-cz314(i7,i2)%ek0*p62(2)+p62k0*cz314(i7,i2)%e(2)
      cauxc=+cz314(i7,i2)%ek0*p578(2)-p578k0*cz314(i7,i2)%e(2)
      u578_314(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u578_314(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u578_314(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u578_314(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u578_314(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u578_314(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u578_314(i7,i2)%d(1)=ccl*cz314(i7,i2)%ek0
      u578_314(i7,i2)%d(2)=ccr*cz314(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id5)/(f62)
      ccl=fcl(id5)/(f62)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p578,qd=p62,v=cf314(i7,i2)%e,a=u578_314(i7,i2)%a,b=u578_314(i
* 7,i2)%b,c=u578_314(i7,i2)%c,d=u578_314(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf314(i7,i2)%ek0*(p578(2)*p62(3)-p62(2)*p578(3))+p
     & 578k0*(cf314(i7,i2)%e(2)*p62(3)-p62(2)*cf314(i7,i2)%e(3))
     & -p62k0*(cf314(i7,i2)%e(2)*p578(3)-p578(2)*cf314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i7,i2)%e(3)*p578k0+p578(3)*cf314(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf314(i7,i2)%e(3)*p62k0+p62(3)*cf314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i7,i2)%e(0)*p578(0)-cf314(i7,i2)%e(1)*p578(1)-c
     & f314(i7,i2)%e(2)*p578(2)-cf314(i7,i2)%e(3)*p578(3)
      cvqd=cf314(i7,i2)%e(0)*p62(0)-cf314(i7,i2)%e(1)*p62(1)-cf3
     & 14(i7,i2)%e(2)*p62(2)-cf314(i7,i2)%e(3)*p62(3)
      cauxa=-cf314(i7,i2)%ek0*quqd+p578k0*cvqd+p62k0*cvqu
      cauxb=-cf314(i7,i2)%ek0*p62(2)+p62k0*cf314(i7,i2)%e(2)
      cauxc=+cf314(i7,i2)%ek0*p578(2)-p578k0*cf314(i7,i2)%e(2)
      u578_314(i7,i2)%a(1)=u578_314(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u578_314(i7,i2)%a(2)=u578_314(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u578_314(i7,i2)%b(1)=u578_314(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u578_314(i7,i2)%b(2)=u578_314(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u578_314(i7,i2)%c(1)=u578_314(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u578_314(i7,i2)%c(2)=u578_314(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u578_314(i7,i2)%d(1)=u578_314(i7,i2)%d(1)+ccl*cf314(i7,i2)
     & %ek0
      u578_314(i7,i2)%d(2)=u578_314(i7,i2)%d(2)+ccr*cf314(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_78314(i5,i7,i2)%a,cc=l5_78314(i5,i7,i2)%c,a1=l5_78(i5)%a
* ,c1=l5_78(i5)%c,a2=u578_314(i7,i2)%a,b2=u578_314(i7,i2)%b,c2=u578_314(
* i7,i2)%c,d2=u578_314(i7,i2)%d,prq=s578,nsum=0
      l5_78314(i5,i7,i2)%a(1)=l5_78(i5)%a(1)*u578_314(i7,i2)%a(1
     & )+l5_78(i5)%c(1)*s578*u578_314(i7,i2)%b(2)
      l5_78314(i5,i7,i2)%c(1)=l5_78(i5)%a(1)*u578_314(i7,i2)%c(1
     & )+l5_78(i5)%c(1)*s578*u578_314(i7,i2)%d(2)
      l5_78314(i5,i7,i2)%c(2)=l5_78(i5)%c(2)*s578*u578_314(i7,i2
     & )%d(1)+l5_78(i5)%a(2)*u578_314(i7,i2)%c(2)
      l5_78314(i5,i7,i2)%a(2)=l5_78(i5)%c(2)*s578*u578_314(i7,i2
     & )%b(1)+l5_78(i5)%a(2)*u578_314(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u5134_78                                               
* quqd -- p=p5134,q=p62
      quqd=p5134(0)*p62(0)-p5134(1)*p62(1)-p5134(2)*p62(2)-p5134
     & (3)*p62(3)
      ccr=zcr(id5)/(f62)
      ccl=zcl(id5)/(f62)
      do i5=1,2
* T0 -- qu=p5134,qd=p62,v=cz78(i5)%e,a=u5314_78(i5)%a,b=u5314_78(i5)%b,c
* =u5314_78(i5)%c,d=u5314_78(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i5)%ek0*(p5134(2)*p62(3)-p62(2)*p5134(3))+p51
     & 34k0*(cz78(i5)%e(2)*p62(3)-p62(2)*cz78(i5)%e(3))-p62k0*(c
     & z78(i5)%e(2)*p5134(3)-p5134(2)*cz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i5)%e(3)*p5134k0+p5134(3)*cz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i5)%e(3)*p62k0+p62(3)*cz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i5)%e(0)*p5134(0)-cz78(i5)%e(1)*p5134(1)-cz78(i5
     & )%e(2)*p5134(2)-cz78(i5)%e(3)*p5134(3)
      cvqd=cz78(i5)%e(0)*p62(0)-cz78(i5)%e(1)*p62(1)-cz78(i5)%e(
     & 2)*p62(2)-cz78(i5)%e(3)*p62(3)
      cauxa=-cz78(i5)%ek0*quqd+p5134k0*cvqd+p62k0*cvqu
      cauxb=-cz78(i5)%ek0*p62(2)+p62k0*cz78(i5)%e(2)
      cauxc=+cz78(i5)%ek0*p5134(2)-p5134k0*cz78(i5)%e(2)
      u5314_78(i5)%a(1)=ccr*(cauxa+ceps_0)
      u5314_78(i5)%a(2)=ccl*(cauxa-ceps_0)
      u5314_78(i5)%b(1)=ccl*(cauxb-ceps_2)
      u5314_78(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u5314_78(i5)%c(1)=ccr*(cauxc+ceps_1)
      u5314_78(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u5314_78(i5)%d(1)=ccl*cz78(i5)%ek0
      u5314_78(i5)%d(2)=ccr*cz78(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f62)
      ccl=fcl(id5)/(f62)
      do i5=1,2
* T0 -- qu=p5134,qd=p62,v=cf78(i5)%e,a=u5314_78(i5)%a,b=u5314_78(i5)%b,c
* =u5314_78(i5)%c,d=u5314_78(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i5)%ek0*(p5134(2)*p62(3)-p62(2)*p5134(3))+p51
     & 34k0*(cf78(i5)%e(2)*p62(3)-p62(2)*cf78(i5)%e(3))-p62k0*(c
     & f78(i5)%e(2)*p5134(3)-p5134(2)*cf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i5)%e(3)*p5134k0+p5134(3)*cf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i5)%e(3)*p62k0+p62(3)*cf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i5)%e(0)*p5134(0)-cf78(i5)%e(1)*p5134(1)-cf78(i5
     & )%e(2)*p5134(2)-cf78(i5)%e(3)*p5134(3)
      cvqd=cf78(i5)%e(0)*p62(0)-cf78(i5)%e(1)*p62(1)-cf78(i5)%e(
     & 2)*p62(2)-cf78(i5)%e(3)*p62(3)
      cauxa=-cf78(i5)%ek0*quqd+p5134k0*cvqd+p62k0*cvqu
      cauxb=-cf78(i5)%ek0*p62(2)+p62k0*cf78(i5)%e(2)
      cauxc=+cf78(i5)%ek0*p5134(2)-p5134k0*cf78(i5)%e(2)
      u5314_78(i5)%a(1)=u5314_78(i5)%a(1)+ccr*(cauxa+ceps_0)
      u5314_78(i5)%a(2)=u5314_78(i5)%a(2)+ccl*(cauxa-ceps_0)
      u5314_78(i5)%b(1)=u5314_78(i5)%b(1)+ccl*(cauxb-ceps_2)
      u5314_78(i5)%b(2)=u5314_78(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u5314_78(i5)%c(1)=u5314_78(i5)%c(1)+ccr*(cauxc+ceps_1)
      u5314_78(i5)%c(2)=u5314_78(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u5314_78(i5)%d(1)=u5314_78(i5)%d(1)+ccl*cf78(i5)%ek0
      u5314_78(i5)%d(2)=u5314_78(i5)%d(2)+ccr*cf78(i5)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l5_78314(i5,i7,i2)%a,cc=l5_78314(i5,i7,i2)%c,a1=l5_314(i7,i
* 2)%a,c1=l5_314(i7,i2)%c,a2=u5314_78(i5)%a,b2=u5314_78(i5)%b,c2=u5314_7
* 8(i5)%c,d2=u5314_78(i5)%d,prq=s5134,nsum=1
      l5_78314(i5,i7,i2)%a(1)=l5_78314(i5,i7,i2)%a(1)+l5_314(i7,
     & i2)%a(1)*u5314_78(i5)%a(1)+l5_314(i7,i2)%c(1)*s5134*u5314
     & _78(i5)%b(2)
      l5_78314(i5,i7,i2)%c(1)=l5_78314(i5,i7,i2)%c(1)+l5_314(i7,
     & i2)%a(1)*u5314_78(i5)%c(1)+l5_314(i7,i2)%c(1)*s5134*u5314
     & _78(i5)%d(2)
      l5_78314(i5,i7,i2)%c(2)=l5_78314(i5,i7,i2)%c(2)+l5_314(i7,
     & i2)%c(2)*s5134*u5314_78(i5)%d(1)+l5_314(i7,i2)%a(2)*u5314
     & _78(i5)%c(2)
      l5_78314(i5,i7,i2)%a(2)=l5_78314(i5,i7,i2)%a(2)+l5_314(i7,
     & i2)%c(2)*s5134*u5314_78(i5)%b(1)+l5_314(i7,i2)%a(2)*u5314
     & _78(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p71,q=p7134
      quqd=p71(0)*p7134(0)-p71(1)*p7134(1)-p71(2)*p7134(2)-p71(3
     & )*p7134(3)
      ccr=zcr(id7)/(f7134)
      ccl=zcl(id7)/(f7134)
      do i5=1,2
* T0 -- qu=p71,qd=p7134,v=cz34(i5)%e,a=u71_34(i5)%a,b=u71_34(i5)%b,c=u71
* _34(i5)%c,d=u71_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p71(2)*p7134(3)-p7134(2)*p71(3))+p71
     & k0*(cz34(i5)%e(2)*p7134(3)-p7134(2)*cz34(i5)%e(3))-p7134k
     & 0*(cz34(i5)%e(2)*p71(3)-p71(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p71k0+p71(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p7134k0+p7134(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p71(0)-cz34(i5)%e(1)*p71(1)-cz34(i5)%e(
     & 2)*p71(2)-cz34(i5)%e(3)*p71(3)
      cvqd=cz34(i5)%e(0)*p7134(0)-cz34(i5)%e(1)*p7134(1)-cz34(i5
     & )%e(2)*p7134(2)-cz34(i5)%e(3)*p7134(3)
      cauxa=-cz34(i5)%ek0*quqd+p71k0*cvqd+p7134k0*cvqu
      cauxb=-cz34(i5)%ek0*p7134(2)+p7134k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p71(2)-p71k0*cz34(i5)%e(2)
      u71_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u71_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u71_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u71_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u71_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u71_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u71_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u71_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f7134)
      ccl=fcl(id7)/(f7134)
      do i5=1,2
* T0 -- qu=p71,qd=p7134,v=cf34(i5)%e,a=u71_34(i5)%a,b=u71_34(i5)%b,c=u71
* _34(i5)%c,d=u71_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p71(2)*p7134(3)-p7134(2)*p71(3))+p71
     & k0*(cf34(i5)%e(2)*p7134(3)-p7134(2)*cf34(i5)%e(3))-p7134k
     & 0*(cf34(i5)%e(2)*p71(3)-p71(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p71k0+p71(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p7134k0+p7134(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p71(0)-cf34(i5)%e(1)*p71(1)-cf34(i5)%e(
     & 2)*p71(2)-cf34(i5)%e(3)*p71(3)
      cvqd=cf34(i5)%e(0)*p7134(0)-cf34(i5)%e(1)*p7134(1)-cf34(i5
     & )%e(2)*p7134(2)-cf34(i5)%e(3)*p7134(3)
      cauxa=-cf34(i5)%ek0*quqd+p71k0*cvqd+p7134k0*cvqu
      cauxb=-cf34(i5)%ek0*p7134(2)+p7134k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p71(2)-p71k0*cf34(i5)%e(2)
      u71_34(i5)%a(1)=u71_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u71_34(i5)%a(2)=u71_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u71_34(i5)%b(1)=u71_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u71_34(i5)%b(2)=u71_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u71_34(i5)%c(1)=u71_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u71_34(i5)%c(2)=u71_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u71_34(i5)%d(1)=u71_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u71_34(i5)%d(2)=u71_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l7_134(i1,i5)%a,cc=l7_134(i1,i5)%c,a1=l7_1(i1)%a,c1=l7_1(i1
* )%c,a2=u71_34(i5)%a,b2=u71_34(i5)%b,c2=u71_34(i5)%c,d2=u71_34(i5)%d,pr
* q=s71,nsum=0
      l7_134(i1,i5)%a(1)=l7_1(i1)%a(1)*u71_34(i5)%a(1)+l7_1(i1)%
     & c(1)*s71*u71_34(i5)%b(2)
      l7_134(i1,i5)%c(1)=l7_1(i1)%a(1)*u71_34(i5)%c(1)+l7_1(i1)%
     & c(1)*s71*u71_34(i5)%d(2)
      l7_134(i1,i5)%c(2)=l7_1(i1)%c(2)*s71*u71_34(i5)%d(1)+l7_1(
     & i1)%a(2)*u71_34(i5)%c(2)
      l7_134(i1,i5)%a(2)=l7_1(i1)%c(2)*s71*u71_34(i5)%b(1)+l7_1(
     & i1)%a(2)*u71_34(i5)%a(2)
      end do
      end do
  
* quqd -- p=p734,q=p7134
      quqd=p734(0)*p7134(0)-p734(1)*p7134(1)-p734(2)*p7134(2)-p7
     & 34(3)*p7134(3)
      ccr=1.d0/(f7134)
      ccl=1.d0/(f7134)
      do i1=1,2
* T0 -- qu=p734,qd=p7134,v=ce1(i1)%e,a=u734_1(i1)%a,b=u734_1(i1)%b,c=u73
* 4_1(i1)%c,d=u734_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p734(2)*p7134(3)-p7134(2)*p734(3))+p7
     & 34k0*(ce1(i1)%e(2)*p7134(3)-p7134(2)*ce1(i1)%e(3))-p7134k
     & 0*(ce1(i1)%e(2)*p734(3)-p734(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p734k0+p734(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p7134k0+p7134(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p734(0)-ce1(i1)%e(1)*p734(1)-ce1(i1)%e(2
     & )*p734(2)-ce1(i1)%e(3)*p734(3)
      cvqd=ce1(i1)%e(0)*p7134(0)-ce1(i1)%e(1)*p7134(1)-ce1(i1)%e
     & (2)*p7134(2)-ce1(i1)%e(3)*p7134(3)
      cauxa=-ce1(i1)%ek0*quqd+p734k0*cvqd+p7134k0*cvqu
      cauxb=-ce1(i1)%ek0*p7134(2)+p7134k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p734(2)-p734k0*ce1(i1)%e(2)
      u734_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u734_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u734_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u734_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u734_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u734_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u734_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u734_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l7_134(i1,i5)%a,cc=l7_134(i1,i5)%c,a1=l7_34(i5)%a,c1=l7_34(
* i5)%c,a2=u734_1(i1)%a,b2=u734_1(i1)%b,c2=u734_1(i1)%c,d2=u734_1(i1)%d,
* prq=s734,nsum=1
      l7_134(i1,i5)%a(1)=l7_134(i1,i5)%a(1)+l7_34(i5)%a(1)*u734_
     & 1(i1)%a(1)+l7_34(i5)%c(1)*s734*u734_1(i1)%b(2)
      l7_134(i1,i5)%c(1)=l7_134(i1,i5)%c(1)+l7_34(i5)%a(1)*u734_
     & 1(i1)%c(1)+l7_34(i5)%c(1)*s734*u734_1(i1)%d(2)
      l7_134(i1,i5)%c(2)=l7_134(i1,i5)%c(2)+l7_34(i5)%c(2)*s734*
     & u734_1(i1)%d(1)+l7_34(i5)%a(2)*u734_1(i1)%c(2)
      l7_134(i1,i5)%a(2)=l7_134(i1,i5)%a(2)+l7_34(i5)%c(2)*s734*
     & u734_1(i1)%b(1)+l7_34(i5)%a(2)*u734_1(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id5).ne.1) then
* quqd -- p=p71,q=p834
      quqd=p71(0)*p834(0)-p71(1)*p834(1)-p71(2)*p834(2)-p71(3)*p
     & 834(3)
      ccr=zcr(id7)/(f834)
      ccl=zcl(id7)/(f834)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p71,qd=p834,v=cz526(i7,i2)%e,a=u71_526(i7,i2)%a,b=u71_526(i7,
* i2)%b,c=u71_526(i7,i2)%c,d=u71_526(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz526(i7,i2)%ek0*(p71(2)*p834(3)-p834(2)*p71(3))+p
     & 71k0*(cz526(i7,i2)%e(2)*p834(3)-p834(2)*cz526(i7,i2)%e(3)
     & )-p834k0*(cz526(i7,i2)%e(2)*p71(3)-p71(2)*cz526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz526(i7,i2)%e(3)*p71k0+p71(3)*cz526(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz526(i7,i2)%e(3)*p834k0+p834(3)*cz526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz526(i7,i2)%e(0)*p71(0)-cz526(i7,i2)%e(1)*p71(1)-cz5
     & 26(i7,i2)%e(2)*p71(2)-cz526(i7,i2)%e(3)*p71(3)
      cvqd=cz526(i7,i2)%e(0)*p834(0)-cz526(i7,i2)%e(1)*p834(1)-c
     & z526(i7,i2)%e(2)*p834(2)-cz526(i7,i2)%e(3)*p834(3)
      cauxa=-cz526(i7,i2)%ek0*quqd+p71k0*cvqd+p834k0*cvqu
      cauxb=-cz526(i7,i2)%ek0*p834(2)+p834k0*cz526(i7,i2)%e(2)
      cauxc=+cz526(i7,i2)%ek0*p71(2)-p71k0*cz526(i7,i2)%e(2)
      u71_526(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u71_526(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u71_526(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u71_526(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u71_526(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u71_526(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u71_526(i7,i2)%d(1)=ccl*cz526(i7,i2)%ek0
      u71_526(i7,i2)%d(2)=ccr*cz526(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id7)/(f834)
      ccl=fcl(id7)/(f834)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p71,qd=p834,v=cf526(i7,i2)%e,a=u71_526(i7,i2)%a,b=u71_526(i7,
* i2)%b,c=u71_526(i7,i2)%c,d=u71_526(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf526(i7,i2)%ek0*(p71(2)*p834(3)-p834(2)*p71(3))+p
     & 71k0*(cf526(i7,i2)%e(2)*p834(3)-p834(2)*cf526(i7,i2)%e(3)
     & )-p834k0*(cf526(i7,i2)%e(2)*p71(3)-p71(2)*cf526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf526(i7,i2)%e(3)*p71k0+p71(3)*cf526(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf526(i7,i2)%e(3)*p834k0+p834(3)*cf526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf526(i7,i2)%e(0)*p71(0)-cf526(i7,i2)%e(1)*p71(1)-cf5
     & 26(i7,i2)%e(2)*p71(2)-cf526(i7,i2)%e(3)*p71(3)
      cvqd=cf526(i7,i2)%e(0)*p834(0)-cf526(i7,i2)%e(1)*p834(1)-c
     & f526(i7,i2)%e(2)*p834(2)-cf526(i7,i2)%e(3)*p834(3)
      cauxa=-cf526(i7,i2)%ek0*quqd+p71k0*cvqd+p834k0*cvqu
      cauxb=-cf526(i7,i2)%ek0*p834(2)+p834k0*cf526(i7,i2)%e(2)
      cauxc=+cf526(i7,i2)%ek0*p71(2)-p71k0*cf526(i7,i2)%e(2)
      u71_526(i7,i2)%a(1)=u71_526(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u71_526(i7,i2)%a(2)=u71_526(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u71_526(i7,i2)%b(1)=u71_526(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u71_526(i7,i2)%b(2)=u71_526(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u71_526(i7,i2)%c(1)=u71_526(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u71_526(i7,i2)%c(2)=u71_526(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u71_526(i7,i2)%d(1)=u71_526(i7,i2)%d(1)+ccl*cf526(i7,i2)%e
     & k0
      u71_526(i7,i2)%d(2)=u71_526(i7,i2)%d(2)+ccr*cf526(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_1526(i1,i7,i2)%a,cc=l7_1526(i1,i7,i2)%c,a1=l7_1(i1)%a,c1
* =l7_1(i1)%c,a2=u71_526(i7,i2)%a,b2=u71_526(i7,i2)%b,c2=u71_526(i7,i2)%
* c,d2=u71_526(i7,i2)%d,prq=s71,nsum=0
      l7_1526(i1,i7,i2)%a(1)=l7_1(i1)%a(1)*u71_526(i7,i2)%a(1)+l
     & 7_1(i1)%c(1)*s71*u71_526(i7,i2)%b(2)
      l7_1526(i1,i7,i2)%c(1)=l7_1(i1)%a(1)*u71_526(i7,i2)%c(1)+l
     & 7_1(i1)%c(1)*s71*u71_526(i7,i2)%d(2)
      l7_1526(i1,i7,i2)%c(2)=l7_1(i1)%c(2)*s71*u71_526(i7,i2)%d(
     & 1)+l7_1(i1)%a(2)*u71_526(i7,i2)%c(2)
      l7_1526(i1,i7,i2)%a(2)=l7_1(i1)%c(2)*s71*u71_526(i7,i2)%b(
     & 1)+l7_1(i1)%a(2)*u71_526(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u7256_1                                                
* quqd -- p=p7256,q=p834
      quqd=p7256(0)*p834(0)-p7256(1)*p834(1)-p7256(2)*p834(2)-p7
     & 256(3)*p834(3)
      ccr=1.d0/(f834)
      ccl=1.d0/(f834)
      do i1=1,2
* T0 -- qu=p7256,qd=p834,v=ce1(i1)%e,a=u7526_1(i1)%a,b=u7526_1(i1)%b,c=u
* 7526_1(i1)%c,d=u7526_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p7256(2)*p834(3)-p834(2)*p7256(3))+p7
     & 256k0*(ce1(i1)%e(2)*p834(3)-p834(2)*ce1(i1)%e(3))-p834k0*
     & (ce1(i1)%e(2)*p7256(3)-p7256(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p7256k0+p7256(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p834k0+p834(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p7256(0)-ce1(i1)%e(1)*p7256(1)-ce1(i1)%e
     & (2)*p7256(2)-ce1(i1)%e(3)*p7256(3)
      cvqd=ce1(i1)%e(0)*p834(0)-ce1(i1)%e(1)*p834(1)-ce1(i1)%e(2
     & )*p834(2)-ce1(i1)%e(3)*p834(3)
      cauxa=-ce1(i1)%ek0*quqd+p7256k0*cvqd+p834k0*cvqu
      cauxb=-ce1(i1)%ek0*p834(2)+p834k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p7256(2)-p7256k0*ce1(i1)%e(2)
      u7526_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u7526_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u7526_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u7526_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u7526_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u7526_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u7526_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u7526_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_1526(i1,i7,i2)%a,cc=l7_1526(i1,i7,i2)%c,a1=l7_526(i7,i2)
* %a,c1=l7_526(i7,i2)%c,a2=u7526_1(i1)%a,b2=u7526_1(i1)%b,c2=u7526_1(i1)
* %c,d2=u7526_1(i1)%d,prq=s7256,nsum=1
      l7_1526(i1,i7,i2)%a(1)=l7_1526(i1,i7,i2)%a(1)+l7_526(i7,i2
     & )%a(1)*u7526_1(i1)%a(1)+l7_526(i7,i2)%c(1)*s7256*u7526_1(
     & i1)%b(2)
      l7_1526(i1,i7,i2)%c(1)=l7_1526(i1,i7,i2)%c(1)+l7_526(i7,i2
     & )%a(1)*u7526_1(i1)%c(1)+l7_526(i7,i2)%c(1)*s7256*u7526_1(
     & i1)%d(2)
      l7_1526(i1,i7,i2)%c(2)=l7_1526(i1,i7,i2)%c(2)+l7_526(i7,i2
     & )%c(2)*s7256*u7526_1(i1)%d(1)+l7_526(i7,i2)%a(2)*u7526_1(
     & i1)%c(2)
      l7_1526(i1,i7,i2)%a(2)=l7_1526(i1,i7,i2)%a(2)+l7_526(i7,i2
     & )%c(2)*s7256*u7526_1(i1)%b(1)+l7_526(i7,i2)%a(2)*u7526_1(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p734,q=p81
      quqd=p734(0)*p81(0)-p734(1)*p81(1)-p734(2)*p81(2)-p734(3)*
     & p81(3)
      ccr=zcr(id7)/(f81)
      ccl=zcl(id7)/(f81)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p734,qd=p81,v=cz526(i7,i2)%e,a=u734_526(i7,i2)%a,b=u734_526(i
* 7,i2)%b,c=u734_526(i7,i2)%c,d=u734_526(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz526(i7,i2)%ek0*(p734(2)*p81(3)-p81(2)*p734(3))+p
     & 734k0*(cz526(i7,i2)%e(2)*p81(3)-p81(2)*cz526(i7,i2)%e(3))
     & -p81k0*(cz526(i7,i2)%e(2)*p734(3)-p734(2)*cz526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz526(i7,i2)%e(3)*p734k0+p734(3)*cz526(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz526(i7,i2)%e(3)*p81k0+p81(3)*cz526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz526(i7,i2)%e(0)*p734(0)-cz526(i7,i2)%e(1)*p734(1)-c
     & z526(i7,i2)%e(2)*p734(2)-cz526(i7,i2)%e(3)*p734(3)
      cvqd=cz526(i7,i2)%e(0)*p81(0)-cz526(i7,i2)%e(1)*p81(1)-cz5
     & 26(i7,i2)%e(2)*p81(2)-cz526(i7,i2)%e(3)*p81(3)
      cauxa=-cz526(i7,i2)%ek0*quqd+p734k0*cvqd+p81k0*cvqu
      cauxb=-cz526(i7,i2)%ek0*p81(2)+p81k0*cz526(i7,i2)%e(2)
      cauxc=+cz526(i7,i2)%ek0*p734(2)-p734k0*cz526(i7,i2)%e(2)
      u734_526(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u734_526(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u734_526(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u734_526(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u734_526(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u734_526(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u734_526(i7,i2)%d(1)=ccl*cz526(i7,i2)%ek0
      u734_526(i7,i2)%d(2)=ccr*cz526(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id7)/(f81)
      ccl=fcl(id7)/(f81)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p734,qd=p81,v=cf526(i7,i2)%e,a=u734_526(i7,i2)%a,b=u734_526(i
* 7,i2)%b,c=u734_526(i7,i2)%c,d=u734_526(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf526(i7,i2)%ek0*(p734(2)*p81(3)-p81(2)*p734(3))+p
     & 734k0*(cf526(i7,i2)%e(2)*p81(3)-p81(2)*cf526(i7,i2)%e(3))
     & -p81k0*(cf526(i7,i2)%e(2)*p734(3)-p734(2)*cf526(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf526(i7,i2)%e(3)*p734k0+p734(3)*cf526(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf526(i7,i2)%e(3)*p81k0+p81(3)*cf526(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf526(i7,i2)%e(0)*p734(0)-cf526(i7,i2)%e(1)*p734(1)-c
     & f526(i7,i2)%e(2)*p734(2)-cf526(i7,i2)%e(3)*p734(3)
      cvqd=cf526(i7,i2)%e(0)*p81(0)-cf526(i7,i2)%e(1)*p81(1)-cf5
     & 26(i7,i2)%e(2)*p81(2)-cf526(i7,i2)%e(3)*p81(3)
      cauxa=-cf526(i7,i2)%ek0*quqd+p734k0*cvqd+p81k0*cvqu
      cauxb=-cf526(i7,i2)%ek0*p81(2)+p81k0*cf526(i7,i2)%e(2)
      cauxc=+cf526(i7,i2)%ek0*p734(2)-p734k0*cf526(i7,i2)%e(2)
      u734_526(i7,i2)%a(1)=u734_526(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u734_526(i7,i2)%a(2)=u734_526(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u734_526(i7,i2)%b(1)=u734_526(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u734_526(i7,i2)%b(2)=u734_526(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u734_526(i7,i2)%c(1)=u734_526(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u734_526(i7,i2)%c(2)=u734_526(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u734_526(i7,i2)%d(1)=u734_526(i7,i2)%d(1)+ccl*cf526(i7,i2)
     & %ek0
      u734_526(i7,i2)%d(2)=u734_526(i7,i2)%d(2)+ccr*cf526(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_34526(i5,i7,i2)%a,cc=l7_34526(i5,i7,i2)%c,a1=l7_34(i5)%a
* ,c1=l7_34(i5)%c,a2=u734_526(i7,i2)%a,b2=u734_526(i7,i2)%b,c2=u734_526(
* i7,i2)%c,d2=u734_526(i7,i2)%d,prq=s734,nsum=0
      l7_34526(i5,i7,i2)%a(1)=l7_34(i5)%a(1)*u734_526(i7,i2)%a(1
     & )+l7_34(i5)%c(1)*s734*u734_526(i7,i2)%b(2)
      l7_34526(i5,i7,i2)%c(1)=l7_34(i5)%a(1)*u734_526(i7,i2)%c(1
     & )+l7_34(i5)%c(1)*s734*u734_526(i7,i2)%d(2)
      l7_34526(i5,i7,i2)%c(2)=l7_34(i5)%c(2)*s734*u734_526(i7,i2
     & )%d(1)+l7_34(i5)%a(2)*u734_526(i7,i2)%c(2)
      l7_34526(i5,i7,i2)%a(2)=l7_34(i5)%c(2)*s734*u734_526(i7,i2
     & )%b(1)+l7_34(i5)%a(2)*u734_526(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u7256_34                                               
* quqd -- p=p7256,q=p81
      quqd=p7256(0)*p81(0)-p7256(1)*p81(1)-p7256(2)*p81(2)-p7256
     & (3)*p81(3)
      ccr=zcr(id7)/(f81)
      ccl=zcl(id7)/(f81)
      do i5=1,2
* T0 -- qu=p7256,qd=p81,v=cz34(i5)%e,a=u7526_34(i5)%a,b=u7526_34(i5)%b,c
* =u7526_34(i5)%c,d=u7526_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p7256(2)*p81(3)-p81(2)*p7256(3))+p72
     & 56k0*(cz34(i5)%e(2)*p81(3)-p81(2)*cz34(i5)%e(3))-p81k0*(c
     & z34(i5)%e(2)*p7256(3)-p7256(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p7256k0+p7256(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p81k0+p81(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p7256(0)-cz34(i5)%e(1)*p7256(1)-cz34(i5
     & )%e(2)*p7256(2)-cz34(i5)%e(3)*p7256(3)
      cvqd=cz34(i5)%e(0)*p81(0)-cz34(i5)%e(1)*p81(1)-cz34(i5)%e(
     & 2)*p81(2)-cz34(i5)%e(3)*p81(3)
      cauxa=-cz34(i5)%ek0*quqd+p7256k0*cvqd+p81k0*cvqu
      cauxb=-cz34(i5)%ek0*p81(2)+p81k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p7256(2)-p7256k0*cz34(i5)%e(2)
      u7526_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u7526_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u7526_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u7526_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u7526_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u7526_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u7526_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u7526_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f81)
      ccl=fcl(id7)/(f81)
      do i5=1,2
* T0 -- qu=p7256,qd=p81,v=cf34(i5)%e,a=u7526_34(i5)%a,b=u7526_34(i5)%b,c
* =u7526_34(i5)%c,d=u7526_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p7256(2)*p81(3)-p81(2)*p7256(3))+p72
     & 56k0*(cf34(i5)%e(2)*p81(3)-p81(2)*cf34(i5)%e(3))-p81k0*(c
     & f34(i5)%e(2)*p7256(3)-p7256(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p7256k0+p7256(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p81k0+p81(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p7256(0)-cf34(i5)%e(1)*p7256(1)-cf34(i5
     & )%e(2)*p7256(2)-cf34(i5)%e(3)*p7256(3)
      cvqd=cf34(i5)%e(0)*p81(0)-cf34(i5)%e(1)*p81(1)-cf34(i5)%e(
     & 2)*p81(2)-cf34(i5)%e(3)*p81(3)
      cauxa=-cf34(i5)%ek0*quqd+p7256k0*cvqd+p81k0*cvqu
      cauxb=-cf34(i5)%ek0*p81(2)+p81k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p7256(2)-p7256k0*cf34(i5)%e(2)
      u7526_34(i5)%a(1)=u7526_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u7526_34(i5)%a(2)=u7526_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u7526_34(i5)%b(1)=u7526_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u7526_34(i5)%b(2)=u7526_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u7526_34(i5)%c(1)=u7526_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u7526_34(i5)%c(2)=u7526_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u7526_34(i5)%d(1)=u7526_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u7526_34(i5)%d(2)=u7526_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_34526(i5,i7,i2)%a,cc=l7_34526(i5,i7,i2)%c,a1=l7_526(i7,i
* 2)%a,c1=l7_526(i7,i2)%c,a2=u7526_34(i5)%a,b2=u7526_34(i5)%b,c2=u7526_3
* 4(i5)%c,d2=u7526_34(i5)%d,prq=s7256,nsum=1
      l7_34526(i5,i7,i2)%a(1)=l7_34526(i5,i7,i2)%a(1)+l7_526(i7,
     & i2)%a(1)*u7526_34(i5)%a(1)+l7_526(i7,i2)%c(1)*s7256*u7526
     & _34(i5)%b(2)
      l7_34526(i5,i7,i2)%c(1)=l7_34526(i5,i7,i2)%c(1)+l7_526(i7,
     & i2)%a(1)*u7526_34(i5)%c(1)+l7_526(i7,i2)%c(1)*s7256*u7526
     & _34(i5)%d(2)
      l7_34526(i5,i7,i2)%c(2)=l7_34526(i5,i7,i2)%c(2)+l7_526(i7,
     & i2)%c(2)*s7256*u7526_34(i5)%d(1)+l7_526(i7,i2)%a(2)*u7526
     & _34(i5)%c(2)
      l7_34526(i5,i7,i2)%a(2)=l7_34526(i5,i7,i2)%a(2)+l7_526(i7,
     & i2)%c(2)*s7256*u7526_34(i5)%b(1)+l7_526(i7,i2)%a(2)*u7526
     & _34(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p72,q=p7234
      quqd=p72(0)*p7234(0)-p72(1)*p7234(1)-p72(2)*p7234(2)-p72(3
     & )*p7234(3)
      ccr=zcr(id7)/(f7234)
      ccl=zcl(id7)/(f7234)
      do i5=1,2
* T0 -- qu=p72,qd=p7234,v=cz34(i5)%e,a=u72_34(i5)%a,b=u72_34(i5)%b,c=u72
* _34(i5)%c,d=u72_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p72(2)*p7234(3)-p7234(2)*p72(3))+p72
     & k0*(cz34(i5)%e(2)*p7234(3)-p7234(2)*cz34(i5)%e(3))-p7234k
     & 0*(cz34(i5)%e(2)*p72(3)-p72(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p72k0+p72(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p7234k0+p7234(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p72(0)-cz34(i5)%e(1)*p72(1)-cz34(i5)%e(
     & 2)*p72(2)-cz34(i5)%e(3)*p72(3)
      cvqd=cz34(i5)%e(0)*p7234(0)-cz34(i5)%e(1)*p7234(1)-cz34(i5
     & )%e(2)*p7234(2)-cz34(i5)%e(3)*p7234(3)
      cauxa=-cz34(i5)%ek0*quqd+p72k0*cvqd+p7234k0*cvqu
      cauxb=-cz34(i5)%ek0*p7234(2)+p7234k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p72(2)-p72k0*cz34(i5)%e(2)
      u72_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u72_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u72_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u72_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u72_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u72_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u72_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u72_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f7234)
      ccl=fcl(id7)/(f7234)
      do i5=1,2
* T0 -- qu=p72,qd=p7234,v=cf34(i5)%e,a=u72_34(i5)%a,b=u72_34(i5)%b,c=u72
* _34(i5)%c,d=u72_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p72(2)*p7234(3)-p7234(2)*p72(3))+p72
     & k0*(cf34(i5)%e(2)*p7234(3)-p7234(2)*cf34(i5)%e(3))-p7234k
     & 0*(cf34(i5)%e(2)*p72(3)-p72(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p72k0+p72(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p7234k0+p7234(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p72(0)-cf34(i5)%e(1)*p72(1)-cf34(i5)%e(
     & 2)*p72(2)-cf34(i5)%e(3)*p72(3)
      cvqd=cf34(i5)%e(0)*p7234(0)-cf34(i5)%e(1)*p7234(1)-cf34(i5
     & )%e(2)*p7234(2)-cf34(i5)%e(3)*p7234(3)
      cauxa=-cf34(i5)%ek0*quqd+p72k0*cvqd+p7234k0*cvqu
      cauxb=-cf34(i5)%ek0*p7234(2)+p7234k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p72(2)-p72k0*cf34(i5)%e(2)
      u72_34(i5)%a(1)=u72_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u72_34(i5)%a(2)=u72_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u72_34(i5)%b(1)=u72_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u72_34(i5)%b(2)=u72_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u72_34(i5)%c(1)=u72_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u72_34(i5)%c(2)=u72_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u72_34(i5)%d(1)=u72_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u72_34(i5)%d(2)=u72_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l7_234(i1,i5)%a,cc=l7_234(i1,i5)%c,a1=l7_2(i1)%a,c1=l7_2(i1
* )%c,a2=u72_34(i5)%a,b2=u72_34(i5)%b,c2=u72_34(i5)%c,d2=u72_34(i5)%d,pr
* q=s72,nsum=0
      l7_234(i1,i5)%a(1)=l7_2(i1)%a(1)*u72_34(i5)%a(1)+l7_2(i1)%
     & c(1)*s72*u72_34(i5)%b(2)
      l7_234(i1,i5)%c(1)=l7_2(i1)%a(1)*u72_34(i5)%c(1)+l7_2(i1)%
     & c(1)*s72*u72_34(i5)%d(2)
      l7_234(i1,i5)%c(2)=l7_2(i1)%c(2)*s72*u72_34(i5)%d(1)+l7_2(
     & i1)%a(2)*u72_34(i5)%c(2)
      l7_234(i1,i5)%a(2)=l7_2(i1)%c(2)*s72*u72_34(i5)%b(1)+l7_2(
     & i1)%a(2)*u72_34(i5)%a(2)
      end do
      end do
  
* quqd -- p=p734,q=p7234
      quqd=p734(0)*p7234(0)-p734(1)*p7234(1)-p734(2)*p7234(2)-p7
     & 34(3)*p7234(3)
      ccr=1.d0/(f7234)
      ccl=1.d0/(f7234)
      do i1=1,2
* T0 -- qu=p734,qd=p7234,v=ce2(i1)%e,a=u734_2(i1)%a,b=u734_2(i1)%b,c=u73
* 4_2(i1)%c,d=u734_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p734(2)*p7234(3)-p7234(2)*p734(3))+p7
     & 34k0*(ce2(i1)%e(2)*p7234(3)-p7234(2)*ce2(i1)%e(3))-p7234k
     & 0*(ce2(i1)%e(2)*p734(3)-p734(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p734k0+p734(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p7234k0+p7234(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p734(0)-ce2(i1)%e(1)*p734(1)-ce2(i1)%e(2
     & )*p734(2)-ce2(i1)%e(3)*p734(3)
      cvqd=ce2(i1)%e(0)*p7234(0)-ce2(i1)%e(1)*p7234(1)-ce2(i1)%e
     & (2)*p7234(2)-ce2(i1)%e(3)*p7234(3)
      cauxa=-ce2(i1)%ek0*quqd+p734k0*cvqd+p7234k0*cvqu
      cauxb=-ce2(i1)%ek0*p7234(2)+p7234k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p734(2)-p734k0*ce2(i1)%e(2)
      u734_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u734_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u734_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u734_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u734_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u734_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u734_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u734_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l7_234(i1,i5)%a,cc=l7_234(i1,i5)%c,a1=l7_34(i5)%a,c1=l7_34(
* i5)%c,a2=u734_2(i1)%a,b2=u734_2(i1)%b,c2=u734_2(i1)%c,d2=u734_2(i1)%d,
* prq=s734,nsum=1
      l7_234(i1,i5)%a(1)=l7_234(i1,i5)%a(1)+l7_34(i5)%a(1)*u734_
     & 2(i1)%a(1)+l7_34(i5)%c(1)*s734*u734_2(i1)%b(2)
      l7_234(i1,i5)%c(1)=l7_234(i1,i5)%c(1)+l7_34(i5)%a(1)*u734_
     & 2(i1)%c(1)+l7_34(i5)%c(1)*s734*u734_2(i1)%d(2)
      l7_234(i1,i5)%c(2)=l7_234(i1,i5)%c(2)+l7_34(i5)%c(2)*s734*
     & u734_2(i1)%d(1)+l7_34(i5)%a(2)*u734_2(i1)%c(2)
      l7_234(i1,i5)%a(2)=l7_234(i1,i5)%a(2)+l7_34(i5)%c(2)*s734*
     & u734_2(i1)%b(1)+l7_34(i5)%a(2)*u734_2(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id5).ne.1) then
* quqd -- p=p72,q=p834
      quqd=p72(0)*p834(0)-p72(1)*p834(1)-p72(2)*p834(2)-p72(3)*p
     & 834(3)
      ccr=zcr(id7)/(f834)
      ccl=zcl(id7)/(f834)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p72,qd=p834,v=cz516(i7,i2)%e,a=u72_516(i7,i2)%a,b=u72_516(i7,
* i2)%b,c=u72_516(i7,i2)%c,d=u72_516(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz516(i7,i2)%ek0*(p72(2)*p834(3)-p834(2)*p72(3))+p
     & 72k0*(cz516(i7,i2)%e(2)*p834(3)-p834(2)*cz516(i7,i2)%e(3)
     & )-p834k0*(cz516(i7,i2)%e(2)*p72(3)-p72(2)*cz516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz516(i7,i2)%e(3)*p72k0+p72(3)*cz516(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz516(i7,i2)%e(3)*p834k0+p834(3)*cz516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz516(i7,i2)%e(0)*p72(0)-cz516(i7,i2)%e(1)*p72(1)-cz5
     & 16(i7,i2)%e(2)*p72(2)-cz516(i7,i2)%e(3)*p72(3)
      cvqd=cz516(i7,i2)%e(0)*p834(0)-cz516(i7,i2)%e(1)*p834(1)-c
     & z516(i7,i2)%e(2)*p834(2)-cz516(i7,i2)%e(3)*p834(3)
      cauxa=-cz516(i7,i2)%ek0*quqd+p72k0*cvqd+p834k0*cvqu
      cauxb=-cz516(i7,i2)%ek0*p834(2)+p834k0*cz516(i7,i2)%e(2)
      cauxc=+cz516(i7,i2)%ek0*p72(2)-p72k0*cz516(i7,i2)%e(2)
      u72_516(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u72_516(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u72_516(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u72_516(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u72_516(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u72_516(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u72_516(i7,i2)%d(1)=ccl*cz516(i7,i2)%ek0
      u72_516(i7,i2)%d(2)=ccr*cz516(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id7)/(f834)
      ccl=fcl(id7)/(f834)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p72,qd=p834,v=cf516(i7,i2)%e,a=u72_516(i7,i2)%a,b=u72_516(i7,
* i2)%b,c=u72_516(i7,i2)%c,d=u72_516(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf516(i7,i2)%ek0*(p72(2)*p834(3)-p834(2)*p72(3))+p
     & 72k0*(cf516(i7,i2)%e(2)*p834(3)-p834(2)*cf516(i7,i2)%e(3)
     & )-p834k0*(cf516(i7,i2)%e(2)*p72(3)-p72(2)*cf516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf516(i7,i2)%e(3)*p72k0+p72(3)*cf516(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf516(i7,i2)%e(3)*p834k0+p834(3)*cf516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf516(i7,i2)%e(0)*p72(0)-cf516(i7,i2)%e(1)*p72(1)-cf5
     & 16(i7,i2)%e(2)*p72(2)-cf516(i7,i2)%e(3)*p72(3)
      cvqd=cf516(i7,i2)%e(0)*p834(0)-cf516(i7,i2)%e(1)*p834(1)-c
     & f516(i7,i2)%e(2)*p834(2)-cf516(i7,i2)%e(3)*p834(3)
      cauxa=-cf516(i7,i2)%ek0*quqd+p72k0*cvqd+p834k0*cvqu
      cauxb=-cf516(i7,i2)%ek0*p834(2)+p834k0*cf516(i7,i2)%e(2)
      cauxc=+cf516(i7,i2)%ek0*p72(2)-p72k0*cf516(i7,i2)%e(2)
      u72_516(i7,i2)%a(1)=u72_516(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u72_516(i7,i2)%a(2)=u72_516(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u72_516(i7,i2)%b(1)=u72_516(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u72_516(i7,i2)%b(2)=u72_516(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u72_516(i7,i2)%c(1)=u72_516(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u72_516(i7,i2)%c(2)=u72_516(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u72_516(i7,i2)%d(1)=u72_516(i7,i2)%d(1)+ccl*cf516(i7,i2)%e
     & k0
      u72_516(i7,i2)%d(2)=u72_516(i7,i2)%d(2)+ccr*cf516(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_2516(i1,i7,i2)%a,cc=l7_2516(i1,i7,i2)%c,a1=l7_2(i1)%a,c1
* =l7_2(i1)%c,a2=u72_516(i7,i2)%a,b2=u72_516(i7,i2)%b,c2=u72_516(i7,i2)%
* c,d2=u72_516(i7,i2)%d,prq=s72,nsum=0
      l7_2516(i1,i7,i2)%a(1)=l7_2(i1)%a(1)*u72_516(i7,i2)%a(1)+l
     & 7_2(i1)%c(1)*s72*u72_516(i7,i2)%b(2)
      l7_2516(i1,i7,i2)%c(1)=l7_2(i1)%a(1)*u72_516(i7,i2)%c(1)+l
     & 7_2(i1)%c(1)*s72*u72_516(i7,i2)%d(2)
      l7_2516(i1,i7,i2)%c(2)=l7_2(i1)%c(2)*s72*u72_516(i7,i2)%d(
     & 1)+l7_2(i1)%a(2)*u72_516(i7,i2)%c(2)
      l7_2516(i1,i7,i2)%a(2)=l7_2(i1)%c(2)*s72*u72_516(i7,i2)%b(
     & 1)+l7_2(i1)%a(2)*u72_516(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u7156_2                                                
* quqd -- p=p7156,q=p834
      quqd=p7156(0)*p834(0)-p7156(1)*p834(1)-p7156(2)*p834(2)-p7
     & 156(3)*p834(3)
      ccr=1.d0/(f834)
      ccl=1.d0/(f834)
      do i1=1,2
* T0 -- qu=p7156,qd=p834,v=ce2(i1)%e,a=u7516_2(i1)%a,b=u7516_2(i1)%b,c=u
* 7516_2(i1)%c,d=u7516_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p7156(2)*p834(3)-p834(2)*p7156(3))+p7
     & 156k0*(ce2(i1)%e(2)*p834(3)-p834(2)*ce2(i1)%e(3))-p834k0*
     & (ce2(i1)%e(2)*p7156(3)-p7156(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p7156k0+p7156(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p834k0+p834(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p7156(0)-ce2(i1)%e(1)*p7156(1)-ce2(i1)%e
     & (2)*p7156(2)-ce2(i1)%e(3)*p7156(3)
      cvqd=ce2(i1)%e(0)*p834(0)-ce2(i1)%e(1)*p834(1)-ce2(i1)%e(2
     & )*p834(2)-ce2(i1)%e(3)*p834(3)
      cauxa=-ce2(i1)%ek0*quqd+p7156k0*cvqd+p834k0*cvqu
      cauxb=-ce2(i1)%ek0*p834(2)+p834k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p7156(2)-p7156k0*ce2(i1)%e(2)
      u7516_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u7516_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u7516_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u7516_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u7516_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u7516_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u7516_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u7516_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_2516(i1,i7,i2)%a,cc=l7_2516(i1,i7,i2)%c,a1=l7_516(i7,i2)
* %a,c1=l7_516(i7,i2)%c,a2=u7516_2(i1)%a,b2=u7516_2(i1)%b,c2=u7516_2(i1)
* %c,d2=u7516_2(i1)%d,prq=s7156,nsum=1
      l7_2516(i1,i7,i2)%a(1)=l7_2516(i1,i7,i2)%a(1)+l7_516(i7,i2
     & )%a(1)*u7516_2(i1)%a(1)+l7_516(i7,i2)%c(1)*s7156*u7516_2(
     & i1)%b(2)
      l7_2516(i1,i7,i2)%c(1)=l7_2516(i1,i7,i2)%c(1)+l7_516(i7,i2
     & )%a(1)*u7516_2(i1)%c(1)+l7_516(i7,i2)%c(1)*s7156*u7516_2(
     & i1)%d(2)
      l7_2516(i1,i7,i2)%c(2)=l7_2516(i1,i7,i2)%c(2)+l7_516(i7,i2
     & )%c(2)*s7156*u7516_2(i1)%d(1)+l7_516(i7,i2)%a(2)*u7516_2(
     & i1)%c(2)
      l7_2516(i1,i7,i2)%a(2)=l7_2516(i1,i7,i2)%a(2)+l7_516(i7,i2
     & )%c(2)*s7156*u7516_2(i1)%b(1)+l7_516(i7,i2)%a(2)*u7516_2(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p734,q=p82
      quqd=p734(0)*p82(0)-p734(1)*p82(1)-p734(2)*p82(2)-p734(3)*
     & p82(3)
      ccr=zcr(id7)/(f82)
      ccl=zcl(id7)/(f82)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p734,qd=p82,v=cz516(i7,i2)%e,a=u734_516(i7,i2)%a,b=u734_516(i
* 7,i2)%b,c=u734_516(i7,i2)%c,d=u734_516(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz516(i7,i2)%ek0*(p734(2)*p82(3)-p82(2)*p734(3))+p
     & 734k0*(cz516(i7,i2)%e(2)*p82(3)-p82(2)*cz516(i7,i2)%e(3))
     & -p82k0*(cz516(i7,i2)%e(2)*p734(3)-p734(2)*cz516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz516(i7,i2)%e(3)*p734k0+p734(3)*cz516(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz516(i7,i2)%e(3)*p82k0+p82(3)*cz516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz516(i7,i2)%e(0)*p734(0)-cz516(i7,i2)%e(1)*p734(1)-c
     & z516(i7,i2)%e(2)*p734(2)-cz516(i7,i2)%e(3)*p734(3)
      cvqd=cz516(i7,i2)%e(0)*p82(0)-cz516(i7,i2)%e(1)*p82(1)-cz5
     & 16(i7,i2)%e(2)*p82(2)-cz516(i7,i2)%e(3)*p82(3)
      cauxa=-cz516(i7,i2)%ek0*quqd+p734k0*cvqd+p82k0*cvqu
      cauxb=-cz516(i7,i2)%ek0*p82(2)+p82k0*cz516(i7,i2)%e(2)
      cauxc=+cz516(i7,i2)%ek0*p734(2)-p734k0*cz516(i7,i2)%e(2)
      u734_516(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u734_516(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u734_516(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u734_516(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u734_516(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u734_516(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u734_516(i7,i2)%d(1)=ccl*cz516(i7,i2)%ek0
      u734_516(i7,i2)%d(2)=ccr*cz516(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id7)/(f82)
      ccl=fcl(id7)/(f82)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p734,qd=p82,v=cf516(i7,i2)%e,a=u734_516(i7,i2)%a,b=u734_516(i
* 7,i2)%b,c=u734_516(i7,i2)%c,d=u734_516(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf516(i7,i2)%ek0*(p734(2)*p82(3)-p82(2)*p734(3))+p
     & 734k0*(cf516(i7,i2)%e(2)*p82(3)-p82(2)*cf516(i7,i2)%e(3))
     & -p82k0*(cf516(i7,i2)%e(2)*p734(3)-p734(2)*cf516(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf516(i7,i2)%e(3)*p734k0+p734(3)*cf516(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf516(i7,i2)%e(3)*p82k0+p82(3)*cf516(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf516(i7,i2)%e(0)*p734(0)-cf516(i7,i2)%e(1)*p734(1)-c
     & f516(i7,i2)%e(2)*p734(2)-cf516(i7,i2)%e(3)*p734(3)
      cvqd=cf516(i7,i2)%e(0)*p82(0)-cf516(i7,i2)%e(1)*p82(1)-cf5
     & 16(i7,i2)%e(2)*p82(2)-cf516(i7,i2)%e(3)*p82(3)
      cauxa=-cf516(i7,i2)%ek0*quqd+p734k0*cvqd+p82k0*cvqu
      cauxb=-cf516(i7,i2)%ek0*p82(2)+p82k0*cf516(i7,i2)%e(2)
      cauxc=+cf516(i7,i2)%ek0*p734(2)-p734k0*cf516(i7,i2)%e(2)
      u734_516(i7,i2)%a(1)=u734_516(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u734_516(i7,i2)%a(2)=u734_516(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u734_516(i7,i2)%b(1)=u734_516(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u734_516(i7,i2)%b(2)=u734_516(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u734_516(i7,i2)%c(1)=u734_516(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u734_516(i7,i2)%c(2)=u734_516(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u734_516(i7,i2)%d(1)=u734_516(i7,i2)%d(1)+ccl*cf516(i7,i2)
     & %ek0
      u734_516(i7,i2)%d(2)=u734_516(i7,i2)%d(2)+ccr*cf516(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_34516(i5,i7,i2)%a,cc=l7_34516(i5,i7,i2)%c,a1=l7_34(i5)%a
* ,c1=l7_34(i5)%c,a2=u734_516(i7,i2)%a,b2=u734_516(i7,i2)%b,c2=u734_516(
* i7,i2)%c,d2=u734_516(i7,i2)%d,prq=s734,nsum=0
      l7_34516(i5,i7,i2)%a(1)=l7_34(i5)%a(1)*u734_516(i7,i2)%a(1
     & )+l7_34(i5)%c(1)*s734*u734_516(i7,i2)%b(2)
      l7_34516(i5,i7,i2)%c(1)=l7_34(i5)%a(1)*u734_516(i7,i2)%c(1
     & )+l7_34(i5)%c(1)*s734*u734_516(i7,i2)%d(2)
      l7_34516(i5,i7,i2)%c(2)=l7_34(i5)%c(2)*s734*u734_516(i7,i2
     & )%d(1)+l7_34(i5)%a(2)*u734_516(i7,i2)%c(2)
      l7_34516(i5,i7,i2)%a(2)=l7_34(i5)%c(2)*s734*u734_516(i7,i2
     & )%b(1)+l7_34(i5)%a(2)*u734_516(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u7156_34                                               
* quqd -- p=p7156,q=p82
      quqd=p7156(0)*p82(0)-p7156(1)*p82(1)-p7156(2)*p82(2)-p7156
     & (3)*p82(3)
      ccr=zcr(id7)/(f82)
      ccl=zcl(id7)/(f82)
      do i5=1,2
* T0 -- qu=p7156,qd=p82,v=cz34(i5)%e,a=u7516_34(i5)%a,b=u7516_34(i5)%b,c
* =u7516_34(i5)%c,d=u7516_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p7156(2)*p82(3)-p82(2)*p7156(3))+p71
     & 56k0*(cz34(i5)%e(2)*p82(3)-p82(2)*cz34(i5)%e(3))-p82k0*(c
     & z34(i5)%e(2)*p7156(3)-p7156(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p7156k0+p7156(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p82k0+p82(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p7156(0)-cz34(i5)%e(1)*p7156(1)-cz34(i5
     & )%e(2)*p7156(2)-cz34(i5)%e(3)*p7156(3)
      cvqd=cz34(i5)%e(0)*p82(0)-cz34(i5)%e(1)*p82(1)-cz34(i5)%e(
     & 2)*p82(2)-cz34(i5)%e(3)*p82(3)
      cauxa=-cz34(i5)%ek0*quqd+p7156k0*cvqd+p82k0*cvqu
      cauxb=-cz34(i5)%ek0*p82(2)+p82k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p7156(2)-p7156k0*cz34(i5)%e(2)
      u7516_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u7516_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u7516_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u7516_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u7516_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u7516_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u7516_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u7516_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f82)
      ccl=fcl(id7)/(f82)
      do i5=1,2
* T0 -- qu=p7156,qd=p82,v=cf34(i5)%e,a=u7516_34(i5)%a,b=u7516_34(i5)%b,c
* =u7516_34(i5)%c,d=u7516_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p7156(2)*p82(3)-p82(2)*p7156(3))+p71
     & 56k0*(cf34(i5)%e(2)*p82(3)-p82(2)*cf34(i5)%e(3))-p82k0*(c
     & f34(i5)%e(2)*p7156(3)-p7156(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p7156k0+p7156(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p82k0+p82(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p7156(0)-cf34(i5)%e(1)*p7156(1)-cf34(i5
     & )%e(2)*p7156(2)-cf34(i5)%e(3)*p7156(3)
      cvqd=cf34(i5)%e(0)*p82(0)-cf34(i5)%e(1)*p82(1)-cf34(i5)%e(
     & 2)*p82(2)-cf34(i5)%e(3)*p82(3)
      cauxa=-cf34(i5)%ek0*quqd+p7156k0*cvqd+p82k0*cvqu
      cauxb=-cf34(i5)%ek0*p82(2)+p82k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p7156(2)-p7156k0*cf34(i5)%e(2)
      u7516_34(i5)%a(1)=u7516_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u7516_34(i5)%a(2)=u7516_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u7516_34(i5)%b(1)=u7516_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u7516_34(i5)%b(2)=u7516_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u7516_34(i5)%c(1)=u7516_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u7516_34(i5)%c(2)=u7516_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u7516_34(i5)%d(1)=u7516_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u7516_34(i5)%d(2)=u7516_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_34516(i5,i7,i2)%a,cc=l7_34516(i5,i7,i2)%c,a1=l7_516(i7,i
* 2)%a,c1=l7_516(i7,i2)%c,a2=u7516_34(i5)%a,b2=u7516_34(i5)%b,c2=u7516_3
* 4(i5)%c,d2=u7516_34(i5)%d,prq=s7156,nsum=1
      l7_34516(i5,i7,i2)%a(1)=l7_34516(i5,i7,i2)%a(1)+l7_516(i7,
     & i2)%a(1)*u7516_34(i5)%a(1)+l7_516(i7,i2)%c(1)*s7156*u7516
     & _34(i5)%b(2)
      l7_34516(i5,i7,i2)%c(1)=l7_34516(i5,i7,i2)%c(1)+l7_516(i7,
     & i2)%a(1)*u7516_34(i5)%c(1)+l7_516(i7,i2)%c(1)*s7156*u7516
     & _34(i5)%d(2)
      l7_34516(i5,i7,i2)%c(2)=l7_34516(i5,i7,i2)%c(2)+l7_516(i7,
     & i2)%c(2)*s7156*u7516_34(i5)%d(1)+l7_516(i7,i2)%a(2)*u7516
     & _34(i5)%c(2)
      l7_34516(i5,i7,i2)%a(2)=l7_34516(i5,i7,i2)%a(2)+l7_516(i7,
     & i2)%c(2)*s7156*u7516_34(i5)%b(1)+l7_516(i7,i2)%a(2)*u7516
     & _34(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p71,q=p7156
      quqd=p71(0)*p7156(0)-p71(1)*p7156(1)-p71(2)*p7156(2)-p71(3
     & )*p7156(3)
      ccr=zcr(id7)/(f7156)
      ccl=zcl(id7)/(f7156)
      do i5=1,2
* T0 -- qu=p71,qd=p7156,v=cz56(i5)%e,a=u71_56(i5)%a,b=u71_56(i5)%b,c=u71
* _56(i5)%c,d=u71_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p71(2)*p7156(3)-p7156(2)*p71(3))+p71
     & k0*(cz56(i5)%e(2)*p7156(3)-p7156(2)*cz56(i5)%e(3))-p7156k
     & 0*(cz56(i5)%e(2)*p71(3)-p71(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p71k0+p71(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p7156k0+p7156(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p71(0)-cz56(i5)%e(1)*p71(1)-cz56(i5)%e(
     & 2)*p71(2)-cz56(i5)%e(3)*p71(3)
      cvqd=cz56(i5)%e(0)*p7156(0)-cz56(i5)%e(1)*p7156(1)-cz56(i5
     & )%e(2)*p7156(2)-cz56(i5)%e(3)*p7156(3)
      cauxa=-cz56(i5)%ek0*quqd+p71k0*cvqd+p7156k0*cvqu
      cauxb=-cz56(i5)%ek0*p7156(2)+p7156k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p71(2)-p71k0*cz56(i5)%e(2)
      u71_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u71_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u71_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u71_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u71_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u71_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u71_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u71_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f7156)
      ccl=fcl(id7)/(f7156)
      do i5=1,2
* T0 -- qu=p71,qd=p7156,v=cf56(i5)%e,a=u71_56(i5)%a,b=u71_56(i5)%b,c=u71
* _56(i5)%c,d=u71_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p71(2)*p7156(3)-p7156(2)*p71(3))+p71
     & k0*(cf56(i5)%e(2)*p7156(3)-p7156(2)*cf56(i5)%e(3))-p7156k
     & 0*(cf56(i5)%e(2)*p71(3)-p71(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p71k0+p71(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p7156k0+p7156(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p71(0)-cf56(i5)%e(1)*p71(1)-cf56(i5)%e(
     & 2)*p71(2)-cf56(i5)%e(3)*p71(3)
      cvqd=cf56(i5)%e(0)*p7156(0)-cf56(i5)%e(1)*p7156(1)-cf56(i5
     & )%e(2)*p7156(2)-cf56(i5)%e(3)*p7156(3)
      cauxa=-cf56(i5)%ek0*quqd+p71k0*cvqd+p7156k0*cvqu
      cauxb=-cf56(i5)%ek0*p7156(2)+p7156k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p71(2)-p71k0*cf56(i5)%e(2)
      u71_56(i5)%a(1)=u71_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u71_56(i5)%a(2)=u71_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u71_56(i5)%b(1)=u71_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u71_56(i5)%b(2)=u71_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u71_56(i5)%c(1)=u71_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u71_56(i5)%c(2)=u71_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u71_56(i5)%d(1)=u71_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u71_56(i5)%d(2)=u71_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l7_156(i1,i5)%a,cc=l7_156(i1,i5)%c,a1=l7_1(i1)%a,c1=l7_1(i1
* )%c,a2=u71_56(i5)%a,b2=u71_56(i5)%b,c2=u71_56(i5)%c,d2=u71_56(i5)%d,pr
* q=s71,nsum=0
      l7_156(i1,i5)%a(1)=l7_1(i1)%a(1)*u71_56(i5)%a(1)+l7_1(i1)%
     & c(1)*s71*u71_56(i5)%b(2)
      l7_156(i1,i5)%c(1)=l7_1(i1)%a(1)*u71_56(i5)%c(1)+l7_1(i1)%
     & c(1)*s71*u71_56(i5)%d(2)
      l7_156(i1,i5)%c(2)=l7_1(i1)%c(2)*s71*u71_56(i5)%d(1)+l7_1(
     & i1)%a(2)*u71_56(i5)%c(2)
      l7_156(i1,i5)%a(2)=l7_1(i1)%c(2)*s71*u71_56(i5)%b(1)+l7_1(
     & i1)%a(2)*u71_56(i5)%a(2)
      end do
      end do
  
* quqd -- p=p756,q=p7156
      quqd=p756(0)*p7156(0)-p756(1)*p7156(1)-p756(2)*p7156(2)-p7
     & 56(3)*p7156(3)
      ccr=1.d0/(f7156)
      ccl=1.d0/(f7156)
      do i1=1,2
* T0 -- qu=p756,qd=p7156,v=ce1(i1)%e,a=u756_1(i1)%a,b=u756_1(i1)%b,c=u75
* 6_1(i1)%c,d=u756_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p756(2)*p7156(3)-p7156(2)*p756(3))+p7
     & 56k0*(ce1(i1)%e(2)*p7156(3)-p7156(2)*ce1(i1)%e(3))-p7156k
     & 0*(ce1(i1)%e(2)*p756(3)-p756(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p756k0+p756(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p7156k0+p7156(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p756(0)-ce1(i1)%e(1)*p756(1)-ce1(i1)%e(2
     & )*p756(2)-ce1(i1)%e(3)*p756(3)
      cvqd=ce1(i1)%e(0)*p7156(0)-ce1(i1)%e(1)*p7156(1)-ce1(i1)%e
     & (2)*p7156(2)-ce1(i1)%e(3)*p7156(3)
      cauxa=-ce1(i1)%ek0*quqd+p756k0*cvqd+p7156k0*cvqu
      cauxb=-ce1(i1)%ek0*p7156(2)+p7156k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p756(2)-p756k0*ce1(i1)%e(2)
      u756_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u756_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u756_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u756_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u756_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u756_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u756_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u756_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l7_156(i1,i5)%a,cc=l7_156(i1,i5)%c,a1=l7_56(i5)%a,c1=l7_56(
* i5)%c,a2=u756_1(i1)%a,b2=u756_1(i1)%b,c2=u756_1(i1)%c,d2=u756_1(i1)%d,
* prq=s756,nsum=1
      l7_156(i1,i5)%a(1)=l7_156(i1,i5)%a(1)+l7_56(i5)%a(1)*u756_
     & 1(i1)%a(1)+l7_56(i5)%c(1)*s756*u756_1(i1)%b(2)
      l7_156(i1,i5)%c(1)=l7_156(i1,i5)%c(1)+l7_56(i5)%a(1)*u756_
     & 1(i1)%c(1)+l7_56(i5)%c(1)*s756*u756_1(i1)%d(2)
      l7_156(i1,i5)%c(2)=l7_156(i1,i5)%c(2)+l7_56(i5)%c(2)*s756*
     & u756_1(i1)%d(1)+l7_56(i5)%a(2)*u756_1(i1)%c(2)
      l7_156(i1,i5)%a(2)=l7_156(i1,i5)%a(2)+l7_56(i5)%c(2)*s756*
     & u756_1(i1)%b(1)+l7_56(i5)%a(2)*u756_1(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id3).ne.1) then
* quqd -- p=p71,q=p856
      quqd=p71(0)*p856(0)-p71(1)*p856(1)-p71(2)*p856(2)-p71(3)*p
     & 856(3)
      ccr=zcr(id7)/(f856)
      ccl=zcl(id7)/(f856)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p71,qd=p856,v=cz324(i7,i2)%e,a=u71_324(i7,i2)%a,b=u71_324(i7,
* i2)%b,c=u71_324(i7,i2)%c,d=u71_324(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz324(i7,i2)%ek0*(p71(2)*p856(3)-p856(2)*p71(3))+p
     & 71k0*(cz324(i7,i2)%e(2)*p856(3)-p856(2)*cz324(i7,i2)%e(3)
     & )-p856k0*(cz324(i7,i2)%e(2)*p71(3)-p71(2)*cz324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i7,i2)%e(3)*p71k0+p71(3)*cz324(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz324(i7,i2)%e(3)*p856k0+p856(3)*cz324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i7,i2)%e(0)*p71(0)-cz324(i7,i2)%e(1)*p71(1)-cz3
     & 24(i7,i2)%e(2)*p71(2)-cz324(i7,i2)%e(3)*p71(3)
      cvqd=cz324(i7,i2)%e(0)*p856(0)-cz324(i7,i2)%e(1)*p856(1)-c
     & z324(i7,i2)%e(2)*p856(2)-cz324(i7,i2)%e(3)*p856(3)
      cauxa=-cz324(i7,i2)%ek0*quqd+p71k0*cvqd+p856k0*cvqu
      cauxb=-cz324(i7,i2)%ek0*p856(2)+p856k0*cz324(i7,i2)%e(2)
      cauxc=+cz324(i7,i2)%ek0*p71(2)-p71k0*cz324(i7,i2)%e(2)
      u71_324(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u71_324(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u71_324(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u71_324(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u71_324(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u71_324(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u71_324(i7,i2)%d(1)=ccl*cz324(i7,i2)%ek0
      u71_324(i7,i2)%d(2)=ccr*cz324(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id7)/(f856)
      ccl=fcl(id7)/(f856)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p71,qd=p856,v=cf324(i7,i2)%e,a=u71_324(i7,i2)%a,b=u71_324(i7,
* i2)%b,c=u71_324(i7,i2)%c,d=u71_324(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf324(i7,i2)%ek0*(p71(2)*p856(3)-p856(2)*p71(3))+p
     & 71k0*(cf324(i7,i2)%e(2)*p856(3)-p856(2)*cf324(i7,i2)%e(3)
     & )-p856k0*(cf324(i7,i2)%e(2)*p71(3)-p71(2)*cf324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i7,i2)%e(3)*p71k0+p71(3)*cf324(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf324(i7,i2)%e(3)*p856k0+p856(3)*cf324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i7,i2)%e(0)*p71(0)-cf324(i7,i2)%e(1)*p71(1)-cf3
     & 24(i7,i2)%e(2)*p71(2)-cf324(i7,i2)%e(3)*p71(3)
      cvqd=cf324(i7,i2)%e(0)*p856(0)-cf324(i7,i2)%e(1)*p856(1)-c
     & f324(i7,i2)%e(2)*p856(2)-cf324(i7,i2)%e(3)*p856(3)
      cauxa=-cf324(i7,i2)%ek0*quqd+p71k0*cvqd+p856k0*cvqu
      cauxb=-cf324(i7,i2)%ek0*p856(2)+p856k0*cf324(i7,i2)%e(2)
      cauxc=+cf324(i7,i2)%ek0*p71(2)-p71k0*cf324(i7,i2)%e(2)
      u71_324(i7,i2)%a(1)=u71_324(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u71_324(i7,i2)%a(2)=u71_324(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u71_324(i7,i2)%b(1)=u71_324(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u71_324(i7,i2)%b(2)=u71_324(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u71_324(i7,i2)%c(1)=u71_324(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u71_324(i7,i2)%c(2)=u71_324(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u71_324(i7,i2)%d(1)=u71_324(i7,i2)%d(1)+ccl*cf324(i7,i2)%e
     & k0
      u71_324(i7,i2)%d(2)=u71_324(i7,i2)%d(2)+ccr*cf324(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_1324(i1,i7,i2)%a,cc=l7_1324(i1,i7,i2)%c,a1=l7_1(i1)%a,c1
* =l7_1(i1)%c,a2=u71_324(i7,i2)%a,b2=u71_324(i7,i2)%b,c2=u71_324(i7,i2)%
* c,d2=u71_324(i7,i2)%d,prq=s71,nsum=0
      l7_1324(i1,i7,i2)%a(1)=l7_1(i1)%a(1)*u71_324(i7,i2)%a(1)+l
     & 7_1(i1)%c(1)*s71*u71_324(i7,i2)%b(2)
      l7_1324(i1,i7,i2)%c(1)=l7_1(i1)%a(1)*u71_324(i7,i2)%c(1)+l
     & 7_1(i1)%c(1)*s71*u71_324(i7,i2)%d(2)
      l7_1324(i1,i7,i2)%c(2)=l7_1(i1)%c(2)*s71*u71_324(i7,i2)%d(
     & 1)+l7_1(i1)%a(2)*u71_324(i7,i2)%c(2)
      l7_1324(i1,i7,i2)%a(2)=l7_1(i1)%c(2)*s71*u71_324(i7,i2)%b(
     & 1)+l7_1(i1)%a(2)*u71_324(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u7234_1                                                
* quqd -- p=p7234,q=p856
      quqd=p7234(0)*p856(0)-p7234(1)*p856(1)-p7234(2)*p856(2)-p7
     & 234(3)*p856(3)
      ccr=1.d0/(f856)
      ccl=1.d0/(f856)
      do i1=1,2
* T0 -- qu=p7234,qd=p856,v=ce1(i1)%e,a=u7324_1(i1)%a,b=u7324_1(i1)%b,c=u
* 7324_1(i1)%c,d=u7324_1(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce1(i1)%ek0*(p7234(2)*p856(3)-p856(2)*p7234(3))+p7
     & 234k0*(ce1(i1)%e(2)*p856(3)-p856(2)*ce1(i1)%e(3))-p856k0*
     & (ce1(i1)%e(2)*p7234(3)-p7234(2)*ce1(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce1(i1)%e(3)*p7234k0+p7234(3)*ce1(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce1(i1)%e(3)*p856k0+p856(3)*ce1(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce1(i1)%e(0)*p7234(0)-ce1(i1)%e(1)*p7234(1)-ce1(i1)%e
     & (2)*p7234(2)-ce1(i1)%e(3)*p7234(3)
      cvqd=ce1(i1)%e(0)*p856(0)-ce1(i1)%e(1)*p856(1)-ce1(i1)%e(2
     & )*p856(2)-ce1(i1)%e(3)*p856(3)
      cauxa=-ce1(i1)%ek0*quqd+p7234k0*cvqd+p856k0*cvqu
      cauxb=-ce1(i1)%ek0*p856(2)+p856k0*ce1(i1)%e(2)
      cauxc=+ce1(i1)%ek0*p7234(2)-p7234k0*ce1(i1)%e(2)
      u7324_1(i1)%a(1)=ccr*(cauxa+ceps_0)
      u7324_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u7324_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u7324_1(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u7324_1(i1)%c(1)=ccr*(cauxc+ceps_1)
      u7324_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u7324_1(i1)%d(1)=ccl*ce1(i1)%ek0
      u7324_1(i1)%d(2)=ccr*ce1(i1)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_1324(i1,i7,i2)%a,cc=l7_1324(i1,i7,i2)%c,a1=l7_324(i7,i2)
* %a,c1=l7_324(i7,i2)%c,a2=u7324_1(i1)%a,b2=u7324_1(i1)%b,c2=u7324_1(i1)
* %c,d2=u7324_1(i1)%d,prq=s7234,nsum=1
      l7_1324(i1,i7,i2)%a(1)=l7_1324(i1,i7,i2)%a(1)+l7_324(i7,i2
     & )%a(1)*u7324_1(i1)%a(1)+l7_324(i7,i2)%c(1)*s7234*u7324_1(
     & i1)%b(2)
      l7_1324(i1,i7,i2)%c(1)=l7_1324(i1,i7,i2)%c(1)+l7_324(i7,i2
     & )%a(1)*u7324_1(i1)%c(1)+l7_324(i7,i2)%c(1)*s7234*u7324_1(
     & i1)%d(2)
      l7_1324(i1,i7,i2)%c(2)=l7_1324(i1,i7,i2)%c(2)+l7_324(i7,i2
     & )%c(2)*s7234*u7324_1(i1)%d(1)+l7_324(i7,i2)%a(2)*u7324_1(
     & i1)%c(2)
      l7_1324(i1,i7,i2)%a(2)=l7_1324(i1,i7,i2)%a(2)+l7_324(i7,i2
     & )%c(2)*s7234*u7324_1(i1)%b(1)+l7_324(i7,i2)%a(2)*u7324_1(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p756,q=p81
      quqd=p756(0)*p81(0)-p756(1)*p81(1)-p756(2)*p81(2)-p756(3)*
     & p81(3)
      ccr=zcr(id7)/(f81)
      ccl=zcl(id7)/(f81)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p756,qd=p81,v=cz324(i7,i2)%e,a=u756_324(i7,i2)%a,b=u756_324(i
* 7,i2)%b,c=u756_324(i7,i2)%c,d=u756_324(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz324(i7,i2)%ek0*(p756(2)*p81(3)-p81(2)*p756(3))+p
     & 756k0*(cz324(i7,i2)%e(2)*p81(3)-p81(2)*cz324(i7,i2)%e(3))
     & -p81k0*(cz324(i7,i2)%e(2)*p756(3)-p756(2)*cz324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i7,i2)%e(3)*p756k0+p756(3)*cz324(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz324(i7,i2)%e(3)*p81k0+p81(3)*cz324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i7,i2)%e(0)*p756(0)-cz324(i7,i2)%e(1)*p756(1)-c
     & z324(i7,i2)%e(2)*p756(2)-cz324(i7,i2)%e(3)*p756(3)
      cvqd=cz324(i7,i2)%e(0)*p81(0)-cz324(i7,i2)%e(1)*p81(1)-cz3
     & 24(i7,i2)%e(2)*p81(2)-cz324(i7,i2)%e(3)*p81(3)
      cauxa=-cz324(i7,i2)%ek0*quqd+p756k0*cvqd+p81k0*cvqu
      cauxb=-cz324(i7,i2)%ek0*p81(2)+p81k0*cz324(i7,i2)%e(2)
      cauxc=+cz324(i7,i2)%ek0*p756(2)-p756k0*cz324(i7,i2)%e(2)
      u756_324(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u756_324(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u756_324(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u756_324(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u756_324(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u756_324(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u756_324(i7,i2)%d(1)=ccl*cz324(i7,i2)%ek0
      u756_324(i7,i2)%d(2)=ccr*cz324(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id7)/(f81)
      ccl=fcl(id7)/(f81)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p756,qd=p81,v=cf324(i7,i2)%e,a=u756_324(i7,i2)%a,b=u756_324(i
* 7,i2)%b,c=u756_324(i7,i2)%c,d=u756_324(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf324(i7,i2)%ek0*(p756(2)*p81(3)-p81(2)*p756(3))+p
     & 756k0*(cf324(i7,i2)%e(2)*p81(3)-p81(2)*cf324(i7,i2)%e(3))
     & -p81k0*(cf324(i7,i2)%e(2)*p756(3)-p756(2)*cf324(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i7,i2)%e(3)*p756k0+p756(3)*cf324(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf324(i7,i2)%e(3)*p81k0+p81(3)*cf324(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i7,i2)%e(0)*p756(0)-cf324(i7,i2)%e(1)*p756(1)-c
     & f324(i7,i2)%e(2)*p756(2)-cf324(i7,i2)%e(3)*p756(3)
      cvqd=cf324(i7,i2)%e(0)*p81(0)-cf324(i7,i2)%e(1)*p81(1)-cf3
     & 24(i7,i2)%e(2)*p81(2)-cf324(i7,i2)%e(3)*p81(3)
      cauxa=-cf324(i7,i2)%ek0*quqd+p756k0*cvqd+p81k0*cvqu
      cauxb=-cf324(i7,i2)%ek0*p81(2)+p81k0*cf324(i7,i2)%e(2)
      cauxc=+cf324(i7,i2)%ek0*p756(2)-p756k0*cf324(i7,i2)%e(2)
      u756_324(i7,i2)%a(1)=u756_324(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u756_324(i7,i2)%a(2)=u756_324(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u756_324(i7,i2)%b(1)=u756_324(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u756_324(i7,i2)%b(2)=u756_324(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u756_324(i7,i2)%c(1)=u756_324(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u756_324(i7,i2)%c(2)=u756_324(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u756_324(i7,i2)%d(1)=u756_324(i7,i2)%d(1)+ccl*cf324(i7,i2)
     & %ek0
      u756_324(i7,i2)%d(2)=u756_324(i7,i2)%d(2)+ccr*cf324(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_56324(i5,i7,i2)%a,cc=l7_56324(i5,i7,i2)%c,a1=l7_56(i5)%a
* ,c1=l7_56(i5)%c,a2=u756_324(i7,i2)%a,b2=u756_324(i7,i2)%b,c2=u756_324(
* i7,i2)%c,d2=u756_324(i7,i2)%d,prq=s756,nsum=0
      l7_56324(i5,i7,i2)%a(1)=l7_56(i5)%a(1)*u756_324(i7,i2)%a(1
     & )+l7_56(i5)%c(1)*s756*u756_324(i7,i2)%b(2)
      l7_56324(i5,i7,i2)%c(1)=l7_56(i5)%a(1)*u756_324(i7,i2)%c(1
     & )+l7_56(i5)%c(1)*s756*u756_324(i7,i2)%d(2)
      l7_56324(i5,i7,i2)%c(2)=l7_56(i5)%c(2)*s756*u756_324(i7,i2
     & )%d(1)+l7_56(i5)%a(2)*u756_324(i7,i2)%c(2)
      l7_56324(i5,i7,i2)%a(2)=l7_56(i5)%c(2)*s756*u756_324(i7,i2
     & )%b(1)+l7_56(i5)%a(2)*u756_324(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u7234_56                                               
* quqd -- p=p7234,q=p81
      quqd=p7234(0)*p81(0)-p7234(1)*p81(1)-p7234(2)*p81(2)-p7234
     & (3)*p81(3)
      ccr=zcr(id7)/(f81)
      ccl=zcl(id7)/(f81)
      do i5=1,2
* T0 -- qu=p7234,qd=p81,v=cz56(i5)%e,a=u7324_56(i5)%a,b=u7324_56(i5)%b,c
* =u7324_56(i5)%c,d=u7324_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p7234(2)*p81(3)-p81(2)*p7234(3))+p72
     & 34k0*(cz56(i5)%e(2)*p81(3)-p81(2)*cz56(i5)%e(3))-p81k0*(c
     & z56(i5)%e(2)*p7234(3)-p7234(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p7234k0+p7234(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p81k0+p81(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p7234(0)-cz56(i5)%e(1)*p7234(1)-cz56(i5
     & )%e(2)*p7234(2)-cz56(i5)%e(3)*p7234(3)
      cvqd=cz56(i5)%e(0)*p81(0)-cz56(i5)%e(1)*p81(1)-cz56(i5)%e(
     & 2)*p81(2)-cz56(i5)%e(3)*p81(3)
      cauxa=-cz56(i5)%ek0*quqd+p7234k0*cvqd+p81k0*cvqu
      cauxb=-cz56(i5)%ek0*p81(2)+p81k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p7234(2)-p7234k0*cz56(i5)%e(2)
      u7324_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u7324_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u7324_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u7324_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u7324_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u7324_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u7324_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u7324_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f81)
      ccl=fcl(id7)/(f81)
      do i5=1,2
* T0 -- qu=p7234,qd=p81,v=cf56(i5)%e,a=u7324_56(i5)%a,b=u7324_56(i5)%b,c
* =u7324_56(i5)%c,d=u7324_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p7234(2)*p81(3)-p81(2)*p7234(3))+p72
     & 34k0*(cf56(i5)%e(2)*p81(3)-p81(2)*cf56(i5)%e(3))-p81k0*(c
     & f56(i5)%e(2)*p7234(3)-p7234(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p7234k0+p7234(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p81k0+p81(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p7234(0)-cf56(i5)%e(1)*p7234(1)-cf56(i5
     & )%e(2)*p7234(2)-cf56(i5)%e(3)*p7234(3)
      cvqd=cf56(i5)%e(0)*p81(0)-cf56(i5)%e(1)*p81(1)-cf56(i5)%e(
     & 2)*p81(2)-cf56(i5)%e(3)*p81(3)
      cauxa=-cf56(i5)%ek0*quqd+p7234k0*cvqd+p81k0*cvqu
      cauxb=-cf56(i5)%ek0*p81(2)+p81k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p7234(2)-p7234k0*cf56(i5)%e(2)
      u7324_56(i5)%a(1)=u7324_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u7324_56(i5)%a(2)=u7324_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u7324_56(i5)%b(1)=u7324_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u7324_56(i5)%b(2)=u7324_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u7324_56(i5)%c(1)=u7324_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u7324_56(i5)%c(2)=u7324_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u7324_56(i5)%d(1)=u7324_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u7324_56(i5)%d(2)=u7324_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_56324(i5,i7,i2)%a,cc=l7_56324(i5,i7,i2)%c,a1=l7_324(i7,i
* 2)%a,c1=l7_324(i7,i2)%c,a2=u7324_56(i5)%a,b2=u7324_56(i5)%b,c2=u7324_5
* 6(i5)%c,d2=u7324_56(i5)%d,prq=s7234,nsum=1
      l7_56324(i5,i7,i2)%a(1)=l7_56324(i5,i7,i2)%a(1)+l7_324(i7,
     & i2)%a(1)*u7324_56(i5)%a(1)+l7_324(i7,i2)%c(1)*s7234*u7324
     & _56(i5)%b(2)
      l7_56324(i5,i7,i2)%c(1)=l7_56324(i5,i7,i2)%c(1)+l7_324(i7,
     & i2)%a(1)*u7324_56(i5)%c(1)+l7_324(i7,i2)%c(1)*s7234*u7324
     & _56(i5)%d(2)
      l7_56324(i5,i7,i2)%c(2)=l7_56324(i5,i7,i2)%c(2)+l7_324(i7,
     & i2)%c(2)*s7234*u7324_56(i5)%d(1)+l7_324(i7,i2)%a(2)*u7324
     & _56(i5)%c(2)
      l7_56324(i5,i7,i2)%a(2)=l7_56324(i5,i7,i2)%a(2)+l7_324(i7,
     & i2)%c(2)*s7234*u7324_56(i5)%b(1)+l7_324(i7,i2)%a(2)*u7324
     & _56(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p72,q=p7256
      quqd=p72(0)*p7256(0)-p72(1)*p7256(1)-p72(2)*p7256(2)-p72(3
     & )*p7256(3)
      ccr=zcr(id7)/(f7256)
      ccl=zcl(id7)/(f7256)
      do i5=1,2
* T0 -- qu=p72,qd=p7256,v=cz56(i5)%e,a=u72_56(i5)%a,b=u72_56(i5)%b,c=u72
* _56(i5)%c,d=u72_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p72(2)*p7256(3)-p7256(2)*p72(3))+p72
     & k0*(cz56(i5)%e(2)*p7256(3)-p7256(2)*cz56(i5)%e(3))-p7256k
     & 0*(cz56(i5)%e(2)*p72(3)-p72(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p72k0+p72(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p7256k0+p7256(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p72(0)-cz56(i5)%e(1)*p72(1)-cz56(i5)%e(
     & 2)*p72(2)-cz56(i5)%e(3)*p72(3)
      cvqd=cz56(i5)%e(0)*p7256(0)-cz56(i5)%e(1)*p7256(1)-cz56(i5
     & )%e(2)*p7256(2)-cz56(i5)%e(3)*p7256(3)
      cauxa=-cz56(i5)%ek0*quqd+p72k0*cvqd+p7256k0*cvqu
      cauxb=-cz56(i5)%ek0*p7256(2)+p7256k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p72(2)-p72k0*cz56(i5)%e(2)
      u72_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u72_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u72_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u72_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u72_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u72_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u72_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u72_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f7256)
      ccl=fcl(id7)/(f7256)
      do i5=1,2
* T0 -- qu=p72,qd=p7256,v=cf56(i5)%e,a=u72_56(i5)%a,b=u72_56(i5)%b,c=u72
* _56(i5)%c,d=u72_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p72(2)*p7256(3)-p7256(2)*p72(3))+p72
     & k0*(cf56(i5)%e(2)*p7256(3)-p7256(2)*cf56(i5)%e(3))-p7256k
     & 0*(cf56(i5)%e(2)*p72(3)-p72(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p72k0+p72(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p7256k0+p7256(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p72(0)-cf56(i5)%e(1)*p72(1)-cf56(i5)%e(
     & 2)*p72(2)-cf56(i5)%e(3)*p72(3)
      cvqd=cf56(i5)%e(0)*p7256(0)-cf56(i5)%e(1)*p7256(1)-cf56(i5
     & )%e(2)*p7256(2)-cf56(i5)%e(3)*p7256(3)
      cauxa=-cf56(i5)%ek0*quqd+p72k0*cvqd+p7256k0*cvqu
      cauxb=-cf56(i5)%ek0*p7256(2)+p7256k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p72(2)-p72k0*cf56(i5)%e(2)
      u72_56(i5)%a(1)=u72_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u72_56(i5)%a(2)=u72_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u72_56(i5)%b(1)=u72_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u72_56(i5)%b(2)=u72_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u72_56(i5)%c(1)=u72_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u72_56(i5)%c(2)=u72_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u72_56(i5)%d(1)=u72_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u72_56(i5)%d(2)=u72_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l7_256(i1,i5)%a,cc=l7_256(i1,i5)%c,a1=l7_2(i1)%a,c1=l7_2(i1
* )%c,a2=u72_56(i5)%a,b2=u72_56(i5)%b,c2=u72_56(i5)%c,d2=u72_56(i5)%d,pr
* q=s72,nsum=0
      l7_256(i1,i5)%a(1)=l7_2(i1)%a(1)*u72_56(i5)%a(1)+l7_2(i1)%
     & c(1)*s72*u72_56(i5)%b(2)
      l7_256(i1,i5)%c(1)=l7_2(i1)%a(1)*u72_56(i5)%c(1)+l7_2(i1)%
     & c(1)*s72*u72_56(i5)%d(2)
      l7_256(i1,i5)%c(2)=l7_2(i1)%c(2)*s72*u72_56(i5)%d(1)+l7_2(
     & i1)%a(2)*u72_56(i5)%c(2)
      l7_256(i1,i5)%a(2)=l7_2(i1)%c(2)*s72*u72_56(i5)%b(1)+l7_2(
     & i1)%a(2)*u72_56(i5)%a(2)
      end do
      end do
  
* quqd -- p=p756,q=p7256
      quqd=p756(0)*p7256(0)-p756(1)*p7256(1)-p756(2)*p7256(2)-p7
     & 56(3)*p7256(3)
      ccr=1.d0/(f7256)
      ccl=1.d0/(f7256)
      do i1=1,2
* T0 -- qu=p756,qd=p7256,v=ce2(i1)%e,a=u756_2(i1)%a,b=u756_2(i1)%b,c=u75
* 6_2(i1)%c,d=u756_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p756(2)*p7256(3)-p7256(2)*p756(3))+p7
     & 56k0*(ce2(i1)%e(2)*p7256(3)-p7256(2)*ce2(i1)%e(3))-p7256k
     & 0*(ce2(i1)%e(2)*p756(3)-p756(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p756k0+p756(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p7256k0+p7256(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p756(0)-ce2(i1)%e(1)*p756(1)-ce2(i1)%e(2
     & )*p756(2)-ce2(i1)%e(3)*p756(3)
      cvqd=ce2(i1)%e(0)*p7256(0)-ce2(i1)%e(1)*p7256(1)-ce2(i1)%e
     & (2)*p7256(2)-ce2(i1)%e(3)*p7256(3)
      cauxa=-ce2(i1)%ek0*quqd+p756k0*cvqd+p7256k0*cvqu
      cauxb=-ce2(i1)%ek0*p7256(2)+p7256k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p756(2)-p756k0*ce2(i1)%e(2)
      u756_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u756_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u756_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u756_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u756_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u756_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u756_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u756_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i5=1,2
* TLT0 -- aa=l7_256(i1,i5)%a,cc=l7_256(i1,i5)%c,a1=l7_56(i5)%a,c1=l7_56(
* i5)%c,a2=u756_2(i1)%a,b2=u756_2(i1)%b,c2=u756_2(i1)%c,d2=u756_2(i1)%d,
* prq=s756,nsum=1
      l7_256(i1,i5)%a(1)=l7_256(i1,i5)%a(1)+l7_56(i5)%a(1)*u756_
     & 2(i1)%a(1)+l7_56(i5)%c(1)*s756*u756_2(i1)%b(2)
      l7_256(i1,i5)%c(1)=l7_256(i1,i5)%c(1)+l7_56(i5)%a(1)*u756_
     & 2(i1)%c(1)+l7_56(i5)%c(1)*s756*u756_2(i1)%d(2)
      l7_256(i1,i5)%c(2)=l7_256(i1,i5)%c(2)+l7_56(i5)%c(2)*s756*
     & u756_2(i1)%d(1)+l7_56(i5)%a(2)*u756_2(i1)%c(2)
      l7_256(i1,i5)%a(2)=l7_256(i1,i5)%a(2)+l7_56(i5)%c(2)*s756*
     & u756_2(i1)%b(1)+l7_56(i5)%a(2)*u756_2(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id3).ne.1) then
* quqd -- p=p72,q=p856
      quqd=p72(0)*p856(0)-p72(1)*p856(1)-p72(2)*p856(2)-p72(3)*p
     & 856(3)
      ccr=zcr(id7)/(f856)
      ccl=zcl(id7)/(f856)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p72,qd=p856,v=cz314(i7,i2)%e,a=u72_314(i7,i2)%a,b=u72_314(i7,
* i2)%b,c=u72_314(i7,i2)%c,d=u72_314(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz314(i7,i2)%ek0*(p72(2)*p856(3)-p856(2)*p72(3))+p
     & 72k0*(cz314(i7,i2)%e(2)*p856(3)-p856(2)*cz314(i7,i2)%e(3)
     & )-p856k0*(cz314(i7,i2)%e(2)*p72(3)-p72(2)*cz314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i7,i2)%e(3)*p72k0+p72(3)*cz314(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz314(i7,i2)%e(3)*p856k0+p856(3)*cz314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i7,i2)%e(0)*p72(0)-cz314(i7,i2)%e(1)*p72(1)-cz3
     & 14(i7,i2)%e(2)*p72(2)-cz314(i7,i2)%e(3)*p72(3)
      cvqd=cz314(i7,i2)%e(0)*p856(0)-cz314(i7,i2)%e(1)*p856(1)-c
     & z314(i7,i2)%e(2)*p856(2)-cz314(i7,i2)%e(3)*p856(3)
      cauxa=-cz314(i7,i2)%ek0*quqd+p72k0*cvqd+p856k0*cvqu
      cauxb=-cz314(i7,i2)%ek0*p856(2)+p856k0*cz314(i7,i2)%e(2)
      cauxc=+cz314(i7,i2)%ek0*p72(2)-p72k0*cz314(i7,i2)%e(2)
      u72_314(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u72_314(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u72_314(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u72_314(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u72_314(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u72_314(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u72_314(i7,i2)%d(1)=ccl*cz314(i7,i2)%ek0
      u72_314(i7,i2)%d(2)=ccr*cz314(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id7)/(f856)
      ccl=fcl(id7)/(f856)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p72,qd=p856,v=cf314(i7,i2)%e,a=u72_314(i7,i2)%a,b=u72_314(i7,
* i2)%b,c=u72_314(i7,i2)%c,d=u72_314(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf314(i7,i2)%ek0*(p72(2)*p856(3)-p856(2)*p72(3))+p
     & 72k0*(cf314(i7,i2)%e(2)*p856(3)-p856(2)*cf314(i7,i2)%e(3)
     & )-p856k0*(cf314(i7,i2)%e(2)*p72(3)-p72(2)*cf314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i7,i2)%e(3)*p72k0+p72(3)*cf314(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf314(i7,i2)%e(3)*p856k0+p856(3)*cf314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i7,i2)%e(0)*p72(0)-cf314(i7,i2)%e(1)*p72(1)-cf3
     & 14(i7,i2)%e(2)*p72(2)-cf314(i7,i2)%e(3)*p72(3)
      cvqd=cf314(i7,i2)%e(0)*p856(0)-cf314(i7,i2)%e(1)*p856(1)-c
     & f314(i7,i2)%e(2)*p856(2)-cf314(i7,i2)%e(3)*p856(3)
      cauxa=-cf314(i7,i2)%ek0*quqd+p72k0*cvqd+p856k0*cvqu
      cauxb=-cf314(i7,i2)%ek0*p856(2)+p856k0*cf314(i7,i2)%e(2)
      cauxc=+cf314(i7,i2)%ek0*p72(2)-p72k0*cf314(i7,i2)%e(2)
      u72_314(i7,i2)%a(1)=u72_314(i7,i2)%a(1)+ccr*(cauxa+ceps_0)
      u72_314(i7,i2)%a(2)=u72_314(i7,i2)%a(2)+ccl*(cauxa-ceps_0)
      u72_314(i7,i2)%b(1)=u72_314(i7,i2)%b(1)+ccl*(cauxb-ceps_2)
      u72_314(i7,i2)%b(2)=u72_314(i7,i2)%b(2)+ccr*(-cauxb-ceps_2
     & )
      u72_314(i7,i2)%c(1)=u72_314(i7,i2)%c(1)+ccr*(cauxc+ceps_1)
      u72_314(i7,i2)%c(2)=u72_314(i7,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u72_314(i7,i2)%d(1)=u72_314(i7,i2)%d(1)+ccl*cf314(i7,i2)%e
     & k0
      u72_314(i7,i2)%d(2)=u72_314(i7,i2)%d(2)+ccr*cf314(i7,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_2314(i1,i7,i2)%a,cc=l7_2314(i1,i7,i2)%c,a1=l7_2(i1)%a,c1
* =l7_2(i1)%c,a2=u72_314(i7,i2)%a,b2=u72_314(i7,i2)%b,c2=u72_314(i7,i2)%
* c,d2=u72_314(i7,i2)%d,prq=s72,nsum=0
      l7_2314(i1,i7,i2)%a(1)=l7_2(i1)%a(1)*u72_314(i7,i2)%a(1)+l
     & 7_2(i1)%c(1)*s72*u72_314(i7,i2)%b(2)
      l7_2314(i1,i7,i2)%c(1)=l7_2(i1)%a(1)*u72_314(i7,i2)%c(1)+l
     & 7_2(i1)%c(1)*s72*u72_314(i7,i2)%d(2)
      l7_2314(i1,i7,i2)%c(2)=l7_2(i1)%c(2)*s72*u72_314(i7,i2)%d(
     & 1)+l7_2(i1)%a(2)*u72_314(i7,i2)%c(2)
      l7_2314(i1,i7,i2)%a(2)=l7_2(i1)%c(2)*s72*u72_314(i7,i2)%b(
     & 1)+l7_2(i1)%a(2)*u72_314(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u7134_2                                                
* quqd -- p=p7134,q=p856
      quqd=p7134(0)*p856(0)-p7134(1)*p856(1)-p7134(2)*p856(2)-p7
     & 134(3)*p856(3)
      ccr=1.d0/(f856)
      ccl=1.d0/(f856)
      do i1=1,2
* T0 -- qu=p7134,qd=p856,v=ce2(i1)%e,a=u7314_2(i1)%a,b=u7314_2(i1)%b,c=u
* 7314_2(i1)%c,d=u7314_2(i1)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ce2(i1)%ek0*(p7134(2)*p856(3)-p856(2)*p7134(3))+p7
     & 134k0*(ce2(i1)%e(2)*p856(3)-p856(2)*ce2(i1)%e(3))-p856k0*
     & (ce2(i1)%e(2)*p7134(3)-p7134(2)*ce2(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ce2(i1)%e(3)*p7134k0+p7134(3)*ce2(i1)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-ce2(i1)%e(3)*p856k0+p856(3)*ce2(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=ce2(i1)%e(0)*p7134(0)-ce2(i1)%e(1)*p7134(1)-ce2(i1)%e
     & (2)*p7134(2)-ce2(i1)%e(3)*p7134(3)
      cvqd=ce2(i1)%e(0)*p856(0)-ce2(i1)%e(1)*p856(1)-ce2(i1)%e(2
     & )*p856(2)-ce2(i1)%e(3)*p856(3)
      cauxa=-ce2(i1)%ek0*quqd+p7134k0*cvqd+p856k0*cvqu
      cauxb=-ce2(i1)%ek0*p856(2)+p856k0*ce2(i1)%e(2)
      cauxc=+ce2(i1)%ek0*p7134(2)-p7134k0*ce2(i1)%e(2)
      u7314_2(i1)%a(1)=ccr*(cauxa+ceps_0)
      u7314_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u7314_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u7314_2(i1)%b(2)=ccr*(-cauxb-ceps_2)
      u7314_2(i1)%c(1)=ccr*(cauxc+ceps_1)
      u7314_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u7314_2(i1)%d(1)=ccl*ce2(i1)%ek0
      u7314_2(i1)%d(2)=ccr*ce2(i1)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_2314(i1,i7,i2)%a,cc=l7_2314(i1,i7,i2)%c,a1=l7_314(i7,i2)
* %a,c1=l7_314(i7,i2)%c,a2=u7314_2(i1)%a,b2=u7314_2(i1)%b,c2=u7314_2(i1)
* %c,d2=u7314_2(i1)%d,prq=s7134,nsum=1
      l7_2314(i1,i7,i2)%a(1)=l7_2314(i1,i7,i2)%a(1)+l7_314(i7,i2
     & )%a(1)*u7314_2(i1)%a(1)+l7_314(i7,i2)%c(1)*s7134*u7314_2(
     & i1)%b(2)
      l7_2314(i1,i7,i2)%c(1)=l7_2314(i1,i7,i2)%c(1)+l7_314(i7,i2
     & )%a(1)*u7314_2(i1)%c(1)+l7_314(i7,i2)%c(1)*s7134*u7314_2(
     & i1)%d(2)
      l7_2314(i1,i7,i2)%c(2)=l7_2314(i1,i7,i2)%c(2)+l7_314(i7,i2
     & )%c(2)*s7134*u7314_2(i1)%d(1)+l7_314(i7,i2)%a(2)*u7314_2(
     & i1)%c(2)
      l7_2314(i1,i7,i2)%a(2)=l7_2314(i1,i7,i2)%a(2)+l7_314(i7,i2
     & )%c(2)*s7134*u7314_2(i1)%b(1)+l7_314(i7,i2)%a(2)*u7314_2(
     & i1)%a(2)
      end do
      end do
      end do
  
* III *cont.                                                            
* quqd -- p=p756,q=p82
      quqd=p756(0)*p82(0)-p756(1)*p82(1)-p756(2)*p82(2)-p756(3)*
     & p82(3)
      ccr=zcr(id7)/(f82)
      ccl=zcl(id7)/(f82)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p756,qd=p82,v=cz314(i7,i2)%e,a=u756_314(i7,i2)%a,b=u756_314(i
* 7,i2)%b,c=u756_314(i7,i2)%c,d=u756_314(i7,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz314(i7,i2)%ek0*(p756(2)*p82(3)-p82(2)*p756(3))+p
     & 756k0*(cz314(i7,i2)%e(2)*p82(3)-p82(2)*cz314(i7,i2)%e(3))
     & -p82k0*(cz314(i7,i2)%e(2)*p756(3)-p756(2)*cz314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i7,i2)%e(3)*p756k0+p756(3)*cz314(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz314(i7,i2)%e(3)*p82k0+p82(3)*cz314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i7,i2)%e(0)*p756(0)-cz314(i7,i2)%e(1)*p756(1)-c
     & z314(i7,i2)%e(2)*p756(2)-cz314(i7,i2)%e(3)*p756(3)
      cvqd=cz314(i7,i2)%e(0)*p82(0)-cz314(i7,i2)%e(1)*p82(1)-cz3
     & 14(i7,i2)%e(2)*p82(2)-cz314(i7,i2)%e(3)*p82(3)
      cauxa=-cz314(i7,i2)%ek0*quqd+p756k0*cvqd+p82k0*cvqu
      cauxb=-cz314(i7,i2)%ek0*p82(2)+p82k0*cz314(i7,i2)%e(2)
      cauxc=+cz314(i7,i2)%ek0*p756(2)-p756k0*cz314(i7,i2)%e(2)
      u756_314(i7,i2)%a(1)=ccr*(cauxa+ceps_0)
      u756_314(i7,i2)%a(2)=ccl*(cauxa-ceps_0)
      u756_314(i7,i2)%b(1)=ccl*(cauxb-ceps_2)
      u756_314(i7,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u756_314(i7,i2)%c(1)=ccr*(cauxc+ceps_1)
      u756_314(i7,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u756_314(i7,i2)%d(1)=ccl*cz314(i7,i2)%ek0
      u756_314(i7,i2)%d(2)=ccr*cz314(i7,i2)%ek0
      end do
      end do
  
      ccr=fcr(id7)/(f82)
      ccl=fcl(id7)/(f82)
      do i7=1,2
      do i2=1,2
* T0 -- qu=p756,qd=p82,v=cf314(i7,i2)%e,a=u756_314(i7,i2)%a,b=u756_314(i
* 7,i2)%b,c=u756_314(i7,i2)%c,d=u756_314(i7,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf314(i7,i2)%ek0*(p756(2)*p82(3)-p82(2)*p756(3))+p
     & 756k0*(cf314(i7,i2)%e(2)*p82(3)-p82(2)*cf314(i7,i2)%e(3))
     & -p82k0*(cf314(i7,i2)%e(2)*p756(3)-p756(2)*cf314(i7,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i7,i2)%e(3)*p756k0+p756(3)*cf314(i7,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf314(i7,i2)%e(3)*p82k0+p82(3)*cf314(i7,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i7,i2)%e(0)*p756(0)-cf314(i7,i2)%e(1)*p756(1)-c
     & f314(i7,i2)%e(2)*p756(2)-cf314(i7,i2)%e(3)*p756(3)
      cvqd=cf314(i7,i2)%e(0)*p82(0)-cf314(i7,i2)%e(1)*p82(1)-cf3
     & 14(i7,i2)%e(2)*p82(2)-cf314(i7,i2)%e(3)*p82(3)
      cauxa=-cf314(i7,i2)%ek0*quqd+p756k0*cvqd+p82k0*cvqu
      cauxb=-cf314(i7,i2)%ek0*p82(2)+p82k0*cf314(i7,i2)%e(2)
      cauxc=+cf314(i7,i2)%ek0*p756(2)-p756k0*cf314(i7,i2)%e(2)
      u756_314(i7,i2)%a(1)=u756_314(i7,i2)%a(1)+ccr*(cauxa+ceps_
     & 0)
      u756_314(i7,i2)%a(2)=u756_314(i7,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u756_314(i7,i2)%b(1)=u756_314(i7,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u756_314(i7,i2)%b(2)=u756_314(i7,i2)%b(2)+ccr*(-cauxb-ceps
     & _2)
      u756_314(i7,i2)%c(1)=u756_314(i7,i2)%c(1)+ccr*(cauxc+ceps_
     & 1)
      u756_314(i7,i2)%c(2)=u756_314(i7,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u756_314(i7,i2)%d(1)=u756_314(i7,i2)%d(1)+ccl*cf314(i7,i2)
     & %ek0
      u756_314(i7,i2)%d(2)=u756_314(i7,i2)%d(2)+ccr*cf314(i7,i2)
     & %ek0
      end do
      end do
  
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_56314(i5,i7,i2)%a,cc=l7_56314(i5,i7,i2)%c,a1=l7_56(i5)%a
* ,c1=l7_56(i5)%c,a2=u756_314(i7,i2)%a,b2=u756_314(i7,i2)%b,c2=u756_314(
* i7,i2)%c,d2=u756_314(i7,i2)%d,prq=s756,nsum=0
      l7_56314(i5,i7,i2)%a(1)=l7_56(i5)%a(1)*u756_314(i7,i2)%a(1
     & )+l7_56(i5)%c(1)*s756*u756_314(i7,i2)%b(2)
      l7_56314(i5,i7,i2)%c(1)=l7_56(i5)%a(1)*u756_314(i7,i2)%c(1
     & )+l7_56(i5)%c(1)*s756*u756_314(i7,i2)%d(2)
      l7_56314(i5,i7,i2)%c(2)=l7_56(i5)%c(2)*s756*u756_314(i7,i2
     & )%d(1)+l7_56(i5)%a(2)*u756_314(i7,i2)%c(2)
      l7_56314(i5,i7,i2)%a(2)=l7_56(i5)%c(2)*s756*u756_314(i7,i2
     & )%b(1)+l7_56(i5)%a(2)*u756_314(i7,i2)%a(2)
      end do
      end do
      end do
       endif
  
* To use also as u7134_56                                               
* quqd -- p=p7134,q=p82
      quqd=p7134(0)*p82(0)-p7134(1)*p82(1)-p7134(2)*p82(2)-p7134
     & (3)*p82(3)
      ccr=zcr(id7)/(f82)
      ccl=zcl(id7)/(f82)
      do i5=1,2
* T0 -- qu=p7134,qd=p82,v=cz56(i5)%e,a=u7314_56(i5)%a,b=u7314_56(i5)%b,c
* =u7314_56(i5)%c,d=u7314_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p7134(2)*p82(3)-p82(2)*p7134(3))+p71
     & 34k0*(cz56(i5)%e(2)*p82(3)-p82(2)*cz56(i5)%e(3))-p82k0*(c
     & z56(i5)%e(2)*p7134(3)-p7134(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p7134k0+p7134(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p82k0+p82(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p7134(0)-cz56(i5)%e(1)*p7134(1)-cz56(i5
     & )%e(2)*p7134(2)-cz56(i5)%e(3)*p7134(3)
      cvqd=cz56(i5)%e(0)*p82(0)-cz56(i5)%e(1)*p82(1)-cz56(i5)%e(
     & 2)*p82(2)-cz56(i5)%e(3)*p82(3)
      cauxa=-cz56(i5)%ek0*quqd+p7134k0*cvqd+p82k0*cvqu
      cauxb=-cz56(i5)%ek0*p82(2)+p82k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p7134(2)-p7134k0*cz56(i5)%e(2)
      u7314_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u7314_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u7314_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u7314_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u7314_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u7314_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u7314_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u7314_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f82)
      ccl=fcl(id7)/(f82)
      do i5=1,2
* T0 -- qu=p7134,qd=p82,v=cf56(i5)%e,a=u7314_56(i5)%a,b=u7314_56(i5)%b,c
* =u7314_56(i5)%c,d=u7314_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p7134(2)*p82(3)-p82(2)*p7134(3))+p71
     & 34k0*(cf56(i5)%e(2)*p82(3)-p82(2)*cf56(i5)%e(3))-p82k0*(c
     & f56(i5)%e(2)*p7134(3)-p7134(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p7134k0+p7134(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p82k0+p82(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p7134(0)-cf56(i5)%e(1)*p7134(1)-cf56(i5
     & )%e(2)*p7134(2)-cf56(i5)%e(3)*p7134(3)
      cvqd=cf56(i5)%e(0)*p82(0)-cf56(i5)%e(1)*p82(1)-cf56(i5)%e(
     & 2)*p82(2)-cf56(i5)%e(3)*p82(3)
      cauxa=-cf56(i5)%ek0*quqd+p7134k0*cvqd+p82k0*cvqu
      cauxb=-cf56(i5)%ek0*p82(2)+p82k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p7134(2)-p7134k0*cf56(i5)%e(2)
      u7314_56(i5)%a(1)=u7314_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u7314_56(i5)%a(2)=u7314_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u7314_56(i5)%b(1)=u7314_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u7314_56(i5)%b(2)=u7314_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u7314_56(i5)%c(1)=u7314_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u7314_56(i5)%c(2)=u7314_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u7314_56(i5)%d(1)=u7314_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u7314_56(i5)%d(2)=u7314_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i5=1,2
      do i7=1,2
      do i2=1,2
* TLT0 -- aa=l7_56314(i5,i7,i2)%a,cc=l7_56314(i5,i7,i2)%c,a1=l7_314(i7,i
* 2)%a,c1=l7_314(i7,i2)%c,a2=u7314_56(i5)%a,b2=u7314_56(i5)%b,c2=u7314_5
* 6(i5)%c,d2=u7314_56(i5)%d,prq=s7134,nsum=1
      l7_56314(i5,i7,i2)%a(1)=l7_56314(i5,i7,i2)%a(1)+l7_314(i7,
     & i2)%a(1)*u7314_56(i5)%a(1)+l7_314(i7,i2)%c(1)*s7134*u7314
     & _56(i5)%b(2)
      l7_56314(i5,i7,i2)%c(1)=l7_56314(i5,i7,i2)%c(1)+l7_314(i7,
     & i2)%a(1)*u7314_56(i5)%c(1)+l7_314(i7,i2)%c(1)*s7134*u7314
     & _56(i5)%d(2)
      l7_56314(i5,i7,i2)%c(2)=l7_56314(i5,i7,i2)%c(2)+l7_314(i7,
     & i2)%c(2)*s7134*u7314_56(i5)%d(1)+l7_314(i7,i2)%a(2)*u7314
     & _56(i5)%c(2)
      l7_56314(i5,i7,i2)%a(2)=l7_56314(i5,i7,i2)%a(2)+l7_314(i7,
     & i2)%c(2)*s7134*u7314_56(i5)%b(1)+l7_314(i7,i2)%a(2)*u7314
     & _56(i5)%a(2)
      end do
      end do
      end do
      endif
      endif
  
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,1),a1=l3_178(i1,i7)%a,c1=l3_178(i1,i7)%
* c,a2=r4_526(i5,i2)%a,b2=r4_526(i5,i2)%b,prq=s3178,bef=cres(i1,i2,&,i5,
* i7,1)+,aft=
      cres(i1,i2,1,i5,i7,1)=cres(i1,i2,1,i5,i7,1)+(l3_178(i1,i7)
     & %a(1)*r4_526(i5,i2)%a(1)+l3_178(i1,i7)%c(1)*s3178*r4_526(
     & i5,i2)%b(2))
      cres(i1,i2,2,i5,i7,1)=cres(i1,i2,2,i5,i7,1)+(l3_178(i1,i7)
     & %c(2)*s3178*r4_526(i5,i2)%b(1)+l3_178(i1,i7)%a(2)*r4_526(
     & i5,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,1),a1=l3_1526(i1,i5,i2)%a,c1=l3_1526(i1
* ,i5,i2)%c,a2=r4_78(i7)%a,b2=r4_78(i7)%b,prq=s478,bef=cres(i1,i2,&,i5,i
* 7,1)+,aft=
      cres(i1,i2,1,i5,i7,1)=cres(i1,i2,1,i5,i7,1)+(l3_1526(i1,i5
     & ,i2)%a(1)*r4_78(i7)%a(1)+l3_1526(i1,i5,i2)%c(1)*s478*r4_7
     & 8(i7)%b(2))
      cres(i1,i2,2,i5,i7,1)=cres(i1,i2,2,i5,i7,1)+(l3_1526(i1,i5
     & ,i2)%c(2)*s478*r4_78(i7)%b(1)+l3_1526(i1,i5,i2)%a(2)*r4_7
     & 8(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,1),a1=l3_78526(i7,i5,i2)%a,c1=l3_78526(
* i7,i5,i2)%c,a2=r4_1(i1)%a,b2=r4_1(i1)%b,prq=s41,bef=cres(i1,i2,&,i5,i7
* ,1)+,aft=
      cres(i1,i2,1,i5,i7,1)=cres(i1,i2,1,i5,i7,1)+(l3_78526(i7,i
     & 5,i2)%a(1)*r4_1(i1)%a(1)+l3_78526(i7,i5,i2)%c(1)*s41*r4_1
     & (i1)%b(2))
      cres(i1,i2,2,i5,i7,1)=cres(i1,i2,2,i5,i7,1)+(l3_78526(i7,i
     & 5,i2)%c(2)*s41*r4_1(i1)%b(1)+l3_78526(i7,i5,i2)%a(2)*r4_1
     & (i1)%a(2))
      end do
      end do
      end do
      end do
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,1),a1=l5_278(i2,i7)%a,c1=l5_278(i2,i7)%
* c,a2=r6_314(i3,i1)%a,b2=r6_314(i3,i1)%b,prq=s5278,bef=cres(i1,i2,i3,&,
* i7,1)+,aft=
      cres(i1,i2,i3,1,i7,1)=cres(i1,i2,i3,1,i7,1)+(l5_278(i2,i7)
     & %a(1)*r6_314(i3,i1)%a(1)+l5_278(i2,i7)%c(1)*s5278*r6_314(
     & i3,i1)%b(2))
      cres(i1,i2,i3,2,i7,1)=cres(i1,i2,i3,2,i7,1)+(l5_278(i2,i7)
     & %c(2)*s5278*r6_314(i3,i1)%b(1)+l5_278(i2,i7)%a(2)*r6_314(
     & i3,i1)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,1),a1=l5_2314(i2,i3,i1)%a,c1=l5_2314(i2
* ,i3,i1)%c,a2=r6_78(i7)%a,b2=r6_78(i7)%b,prq=s678,bef=cres(i1,i2,i3,&,i
* 7,1)+,aft=
      cres(i1,i2,i3,1,i7,1)=cres(i1,i2,i3,1,i7,1)+(l5_2314(i2,i3
     & ,i1)%a(1)*r6_78(i7)%a(1)+l5_2314(i2,i3,i1)%c(1)*s678*r6_7
     & 8(i7)%b(2))
      cres(i1,i2,i3,2,i7,1)=cres(i1,i2,i3,2,i7,1)+(l5_2314(i2,i3
     & ,i1)%c(2)*s678*r6_78(i7)%b(1)+l5_2314(i2,i3,i1)%a(2)*r6_7
     & 8(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,1),a1=l5_78314(i7,i3,i1)%a,c1=l5_78314(
* i7,i3,i1)%c,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s62,bef=cres(i1,i2,i3,&,i7
* ,1)+,aft=
      cres(i1,i2,i3,1,i7,1)=cres(i1,i2,i3,1,i7,1)+(l5_78314(i7,i
     & 3,i1)%a(1)*r6_2(i2)%a(1)+l5_78314(i7,i3,i1)%c(1)*s62*r6_2
     & (i2)%b(2))
      cres(i1,i2,i3,2,i7,1)=cres(i1,i2,i3,2,i7,1)+(l5_78314(i7,i
     & 3,i1)%c(2)*s62*r6_2(i2)%b(1)+l5_78314(i7,i3,i1)%a(2)*r6_2
     & (i2)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,2),a1=l3_278(i2,i7)%a,c1=l3_278(i2,i7)%
* c,a2=r4_516(i5,i1)%a,b2=r4_516(i5,i1)%b,prq=s3278,bef=cres(i1,i2,&,i5,
* i7,2)+,aft=
      cres(i1,i2,1,i5,i7,2)=cres(i1,i2,1,i5,i7,2)+(l3_278(i2,i7)
     & %a(1)*r4_516(i5,i1)%a(1)+l3_278(i2,i7)%c(1)*s3278*r4_516(
     & i5,i1)%b(2))
      cres(i1,i2,2,i5,i7,2)=cres(i1,i2,2,i5,i7,2)+(l3_278(i2,i7)
     & %c(2)*s3278*r4_516(i5,i1)%b(1)+l3_278(i2,i7)%a(2)*r4_516(
     & i5,i1)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,2),a1=l3_2516(i2,i5,i1)%a,c1=l3_2516(i2
* ,i5,i1)%c,a2=r4_78(i7)%a,b2=r4_78(i7)%b,prq=s478,bef=cres(i1,i2,&,i5,i
* 7,2)+,aft=
      cres(i1,i2,1,i5,i7,2)=cres(i1,i2,1,i5,i7,2)+(l3_2516(i2,i5
     & ,i1)%a(1)*r4_78(i7)%a(1)+l3_2516(i2,i5,i1)%c(1)*s478*r4_7
     & 8(i7)%b(2))
      cres(i1,i2,2,i5,i7,2)=cres(i1,i2,2,i5,i7,2)+(l3_2516(i2,i5
     & ,i1)%c(2)*s478*r4_78(i7)%b(1)+l3_2516(i2,i5,i1)%a(2)*r4_7
     & 8(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,2),a1=l3_78516(i7,i5,i1)%a,c1=l3_78516(
* i7,i5,i1)%c,a2=r4_2(i2)%a,b2=r4_2(i2)%b,prq=s42,bef=cres(i1,i2,&,i5,i7
* ,2)+,aft=
      cres(i1,i2,1,i5,i7,2)=cres(i1,i2,1,i5,i7,2)+(l3_78516(i7,i
     & 5,i1)%a(1)*r4_2(i2)%a(1)+l3_78516(i7,i5,i1)%c(1)*s42*r4_2
     & (i2)%b(2))
      cres(i1,i2,2,i5,i7,2)=cres(i1,i2,2,i5,i7,2)+(l3_78516(i7,i
     & 5,i1)%c(2)*s42*r4_2(i2)%b(1)+l3_78516(i7,i5,i1)%a(2)*r4_2
     & (i2)%a(2))
      end do
      end do
      end do
      end do
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,2),a1=l5_178(i1,i7)%a,c1=l5_178(i1,i7)%
* c,a2=r6_324(i3,i2)%a,b2=r6_324(i3,i2)%b,prq=s5178,bef=cres(i1,i2,i3,&,
* i7,2)+,aft=
      cres(i1,i2,i3,1,i7,2)=cres(i1,i2,i3,1,i7,2)+(l5_178(i1,i7)
     & %a(1)*r6_324(i3,i2)%a(1)+l5_178(i1,i7)%c(1)*s5178*r6_324(
     & i3,i2)%b(2))
      cres(i1,i2,i3,2,i7,2)=cres(i1,i2,i3,2,i7,2)+(l5_178(i1,i7)
     & %c(2)*s5178*r6_324(i3,i2)%b(1)+l5_178(i1,i7)%a(2)*r6_324(
     & i3,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,2),a1=l5_1324(i1,i3,i2)%a,c1=l5_1324(i1
* ,i3,i2)%c,a2=r6_78(i7)%a,b2=r6_78(i7)%b,prq=s678,bef=cres(i1,i2,i3,&,i
* 7,2)+,aft=
      cres(i1,i2,i3,1,i7,2)=cres(i1,i2,i3,1,i7,2)+(l5_1324(i1,i3
     & ,i2)%a(1)*r6_78(i7)%a(1)+l5_1324(i1,i3,i2)%c(1)*s678*r6_7
     & 8(i7)%b(2))
      cres(i1,i2,i3,2,i7,2)=cres(i1,i2,i3,2,i7,2)+(l5_1324(i1,i3
     & ,i2)%c(2)*s678*r6_78(i7)%b(1)+l5_1324(i1,i3,i2)%a(2)*r6_7
     & 8(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,2),a1=l5_78324(i7,i3,i2)%a,c1=l5_78324(
* i7,i3,i2)%c,a2=r6_1(i1)%a,b2=r6_1(i1)%b,prq=s61,bef=cres(i1,i2,i3,&,i7
* ,2)+,aft=
      cres(i1,i2,i3,1,i7,2)=cres(i1,i2,i3,1,i7,2)+(l5_78324(i7,i
     & 3,i2)%a(1)*r6_1(i1)%a(1)+l5_78324(i7,i3,i2)%c(1)*s61*r6_1
     & (i1)%b(2))
      cres(i1,i2,i3,2,i7,2)=cres(i1,i2,i3,2,i7,2)+(l5_78324(i7,i
     & 3,i2)%c(2)*s61*r6_1(i1)%b(1)+l5_78324(i7,i3,i2)%a(2)*r6_1
     & (i1)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,3),a1=l3_156(i1,i5)%a,c1=l3_156(i1,i5)%
* c,a2=r4_728(i7,i2)%a,b2=r4_728(i7,i2)%b,prq=s3156,bef=cres(i1,i2,&,i5,
* i7,3)+,aft=
      cres(i1,i2,1,i5,i7,3)=cres(i1,i2,1,i5,i7,3)+(l3_156(i1,i5)
     & %a(1)*r4_728(i7,i2)%a(1)+l3_156(i1,i5)%c(1)*s3156*r4_728(
     & i7,i2)%b(2))
      cres(i1,i2,2,i5,i7,3)=cres(i1,i2,2,i5,i7,3)+(l3_156(i1,i5)
     & %c(2)*s3156*r4_728(i7,i2)%b(1)+l3_156(i1,i5)%a(2)*r4_728(
     & i7,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,3),a1=l3_1728(i1,i7,i2)%a,c1=l3_1728(i1
* ,i7,i2)%c,a2=r4_56(i5)%a,b2=r4_56(i5)%b,prq=s456,bef=cres(i1,i2,&,i5,i
* 7,3)+,aft=
      cres(i1,i2,1,i5,i7,3)=cres(i1,i2,1,i5,i7,3)+(l3_1728(i1,i7
     & ,i2)%a(1)*r4_56(i5)%a(1)+l3_1728(i1,i7,i2)%c(1)*s456*r4_5
     & 6(i5)%b(2))
      cres(i1,i2,2,i5,i7,3)=cres(i1,i2,2,i5,i7,3)+(l3_1728(i1,i7
     & ,i2)%c(2)*s456*r4_56(i5)%b(1)+l3_1728(i1,i7,i2)%a(2)*r4_5
     & 6(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,3),a1=l3_56728(i5,i7,i2)%a,c1=l3_56728(
* i5,i7,i2)%c,a2=r4_1(i1)%a,b2=r4_1(i1)%b,prq=s41,bef=cres(i1,i2,&,i5,i7
* ,3)+,aft=
      cres(i1,i2,1,i5,i7,3)=cres(i1,i2,1,i5,i7,3)+(l3_56728(i5,i
     & 7,i2)%a(1)*r4_1(i1)%a(1)+l3_56728(i5,i7,i2)%c(1)*s41*r4_1
     & (i1)%b(2))
      cres(i1,i2,2,i5,i7,3)=cres(i1,i2,2,i5,i7,3)+(l3_56728(i5,i
     & 7,i2)%c(2)*s41*r4_1(i1)%b(1)+l3_56728(i5,i7,i2)%a(2)*r4_1
     & (i1)%a(2))
      end do
      end do
      end do
      end do
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,3),a1=l7_256(i2,i5)%a,c1=l7_256(i2,i5)%
* c,a2=r8_314(i3,i1)%a,b2=r8_314(i3,i1)%b,prq=s7256,bef=cres(i1,i2,i3,i5
* ,&,3)+,aft=
      cres(i1,i2,i3,i5,1,3)=cres(i1,i2,i3,i5,1,3)+(l7_256(i2,i5)
     & %a(1)*r8_314(i3,i1)%a(1)+l7_256(i2,i5)%c(1)*s7256*r8_314(
     & i3,i1)%b(2))
      cres(i1,i2,i3,i5,2,3)=cres(i1,i2,i3,i5,2,3)+(l7_256(i2,i5)
     & %c(2)*s7256*r8_314(i3,i1)%b(1)+l7_256(i2,i5)%a(2)*r8_314(
     & i3,i1)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,3),a1=l7_2314(i2,i3,i1)%a,c1=l7_2314(i2
* ,i3,i1)%c,a2=r8_56(i5)%a,b2=r8_56(i5)%b,prq=s856,bef=cres(i1,i2,i3,i5,
* &,3)+,aft=
      cres(i1,i2,i3,i5,1,3)=cres(i1,i2,i3,i5,1,3)+(l7_2314(i2,i3
     & ,i1)%a(1)*r8_56(i5)%a(1)+l7_2314(i2,i3,i1)%c(1)*s856*r8_5
     & 6(i5)%b(2))
      cres(i1,i2,i3,i5,2,3)=cres(i1,i2,i3,i5,2,3)+(l7_2314(i2,i3
     & ,i1)%c(2)*s856*r8_56(i5)%b(1)+l7_2314(i2,i3,i1)%a(2)*r8_5
     & 6(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,3),a1=l7_56314(i5,i3,i1)%a,c1=l7_56314(
* i5,i3,i1)%c,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s82,bef=cres(i1,i2,i3,i5,&
* ,3)+,aft=
      cres(i1,i2,i3,i5,1,3)=cres(i1,i2,i3,i5,1,3)+(l7_56314(i5,i
     & 3,i1)%a(1)*r8_2(i2)%a(1)+l7_56314(i5,i3,i1)%c(1)*s82*r8_2
     & (i2)%b(2))
      cres(i1,i2,i3,i5,2,3)=cres(i1,i2,i3,i5,2,3)+(l7_56314(i5,i
     & 3,i1)%c(2)*s82*r8_2(i2)%b(1)+l7_56314(i5,i3,i1)%a(2)*r8_2
     & (i2)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,4),a1=l3_256(i2,i5)%a,c1=l3_256(i2,i5)%
* c,a2=r4_718(i7,i1)%a,b2=r4_718(i7,i1)%b,prq=s3256,bef=cres(i1,i2,&,i5,
* i7,4)+,aft=
      cres(i1,i2,1,i5,i7,4)=cres(i1,i2,1,i5,i7,4)+(l3_256(i2,i5)
     & %a(1)*r4_718(i7,i1)%a(1)+l3_256(i2,i5)%c(1)*s3256*r4_718(
     & i7,i1)%b(2))
      cres(i1,i2,2,i5,i7,4)=cres(i1,i2,2,i5,i7,4)+(l3_256(i2,i5)
     & %c(2)*s3256*r4_718(i7,i1)%b(1)+l3_256(i2,i5)%a(2)*r4_718(
     & i7,i1)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,4),a1=l3_2718(i2,i7,i1)%a,c1=l3_2718(i2
* ,i7,i1)%c,a2=r4_56(i5)%a,b2=r4_56(i5)%b,prq=s456,bef=cres(i1,i2,&,i5,i
* 7,4)+,aft=
      cres(i1,i2,1,i5,i7,4)=cres(i1,i2,1,i5,i7,4)+(l3_2718(i2,i7
     & ,i1)%a(1)*r4_56(i5)%a(1)+l3_2718(i2,i7,i1)%c(1)*s456*r4_5
     & 6(i5)%b(2))
      cres(i1,i2,2,i5,i7,4)=cres(i1,i2,2,i5,i7,4)+(l3_2718(i2,i7
     & ,i1)%c(2)*s456*r4_56(i5)%b(1)+l3_2718(i2,i7,i1)%a(2)*r4_5
     & 6(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,4),a1=l3_56718(i5,i7,i1)%a,c1=l3_56718(
* i5,i7,i1)%c,a2=r4_2(i2)%a,b2=r4_2(i2)%b,prq=s42,bef=cres(i1,i2,&,i5,i7
* ,4)+,aft=
      cres(i1,i2,1,i5,i7,4)=cres(i1,i2,1,i5,i7,4)+(l3_56718(i5,i
     & 7,i1)%a(1)*r4_2(i2)%a(1)+l3_56718(i5,i7,i1)%c(1)*s42*r4_2
     & (i2)%b(2))
      cres(i1,i2,2,i5,i7,4)=cres(i1,i2,2,i5,i7,4)+(l3_56718(i5,i
     & 7,i1)%c(2)*s42*r4_2(i2)%b(1)+l3_56718(i5,i7,i1)%a(2)*r4_2
     & (i2)%a(2))
      end do
      end do
      end do
      end do
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,4),a1=l7_156(i1,i5)%a,c1=l7_156(i1,i5)%
* c,a2=r8_324(i3,i2)%a,b2=r8_324(i3,i2)%b,prq=s7156,bef=cres(i1,i2,i3,i5
* ,&,4)+,aft=
      cres(i1,i2,i3,i5,1,4)=cres(i1,i2,i3,i5,1,4)+(l7_156(i1,i5)
     & %a(1)*r8_324(i3,i2)%a(1)+l7_156(i1,i5)%c(1)*s7156*r8_324(
     & i3,i2)%b(2))
      cres(i1,i2,i3,i5,2,4)=cres(i1,i2,i3,i5,2,4)+(l7_156(i1,i5)
     & %c(2)*s7156*r8_324(i3,i2)%b(1)+l7_156(i1,i5)%a(2)*r8_324(
     & i3,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,4),a1=l7_1324(i1,i3,i2)%a,c1=l7_1324(i1
* ,i3,i2)%c,a2=r8_56(i5)%a,b2=r8_56(i5)%b,prq=s856,bef=cres(i1,i2,i3,i5,
* &,4)+,aft=
      cres(i1,i2,i3,i5,1,4)=cres(i1,i2,i3,i5,1,4)+(l7_1324(i1,i3
     & ,i2)%a(1)*r8_56(i5)%a(1)+l7_1324(i1,i3,i2)%c(1)*s856*r8_5
     & 6(i5)%b(2))
      cres(i1,i2,i3,i5,2,4)=cres(i1,i2,i3,i5,2,4)+(l7_1324(i1,i3
     & ,i2)%c(2)*s856*r8_56(i5)%b(1)+l7_1324(i1,i3,i2)%a(2)*r8_5
     & 6(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,4),a1=l7_56324(i5,i3,i2)%a,c1=l7_56324(
* i5,i3,i2)%c,a2=r8_1(i1)%a,b2=r8_1(i1)%b,prq=s81,bef=cres(i1,i2,i3,i5,&
* ,4)+,aft=
      cres(i1,i2,i3,i5,1,4)=cres(i1,i2,i3,i5,1,4)+(l7_56324(i5,i
     & 3,i2)%a(1)*r8_1(i1)%a(1)+l7_56324(i5,i3,i2)%c(1)*s81*r8_1
     & (i1)%b(2))
      cres(i1,i2,i3,i5,2,4)=cres(i1,i2,i3,i5,2,4)+(l7_56324(i5,i
     & 3,i2)%c(2)*s81*r8_1(i1)%b(1)+l7_56324(i5,i3,i2)%a(2)*r8_1
     & (i1)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,5),a1=l5_134(i1,i3)%a,c1=l5_134(i1,i3)%
* c,a2=r6_728(i7,i2)%a,b2=r6_728(i7,i2)%b,prq=s5134,bef=cres(i1,i2,i3,&,
* i7,5)+,aft=
      cres(i1,i2,i3,1,i7,5)=cres(i1,i2,i3,1,i7,5)+(l5_134(i1,i3)
     & %a(1)*r6_728(i7,i2)%a(1)+l5_134(i1,i3)%c(1)*s5134*r6_728(
     & i7,i2)%b(2))
      cres(i1,i2,i3,2,i7,5)=cres(i1,i2,i3,2,i7,5)+(l5_134(i1,i3)
     & %c(2)*s5134*r6_728(i7,i2)%b(1)+l5_134(i1,i3)%a(2)*r6_728(
     & i7,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,5),a1=l5_1728(i1,i7,i2)%a,c1=l5_1728(i1
* ,i7,i2)%c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=cres(i1,i2,i3,&,i
* 7,5)+,aft=
      cres(i1,i2,i3,1,i7,5)=cres(i1,i2,i3,1,i7,5)+(l5_1728(i1,i7
     & ,i2)%a(1)*r6_34(i3)%a(1)+l5_1728(i1,i7,i2)%c(1)*s634*r6_3
     & 4(i3)%b(2))
      cres(i1,i2,i3,2,i7,5)=cres(i1,i2,i3,2,i7,5)+(l5_1728(i1,i7
     & ,i2)%c(2)*s634*r6_34(i3)%b(1)+l5_1728(i1,i7,i2)%a(2)*r6_3
     & 4(i3)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,5),a1=l5_34728(i3,i7,i2)%a,c1=l5_34728(
* i3,i7,i2)%c,a2=r6_1(i1)%a,b2=r6_1(i1)%b,prq=s61,bef=cres(i1,i2,i3,&,i7
* ,5)+,aft=
      cres(i1,i2,i3,1,i7,5)=cres(i1,i2,i3,1,i7,5)+(l5_34728(i3,i
     & 7,i2)%a(1)*r6_1(i1)%a(1)+l5_34728(i3,i7,i2)%c(1)*s61*r6_1
     & (i1)%b(2))
      cres(i1,i2,i3,2,i7,5)=cres(i1,i2,i3,2,i7,5)+(l5_34728(i3,i
     & 7,i2)%c(2)*s61*r6_1(i1)%b(1)+l5_34728(i3,i7,i2)%a(2)*r6_1
     & (i1)%a(2))
      end do
      end do
      end do
      end do
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,5),a1=l7_234(i2,i3)%a,c1=l7_234(i2,i3)%
* c,a2=r8_516(i5,i1)%a,b2=r8_516(i5,i1)%b,prq=s7234,bef=cres(i1,i2,i3,i5
* ,&,5)+,aft=
      cres(i1,i2,i3,i5,1,5)=cres(i1,i2,i3,i5,1,5)+(l7_234(i2,i3)
     & %a(1)*r8_516(i5,i1)%a(1)+l7_234(i2,i3)%c(1)*s7234*r8_516(
     & i5,i1)%b(2))
      cres(i1,i2,i3,i5,2,5)=cres(i1,i2,i3,i5,2,5)+(l7_234(i2,i3)
     & %c(2)*s7234*r8_516(i5,i1)%b(1)+l7_234(i2,i3)%a(2)*r8_516(
     & i5,i1)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,5),a1=l7_2516(i2,i5,i1)%a,c1=l7_2516(i2
* ,i5,i1)%c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=cres(i1,i2,i3,i5,
* &,5)+,aft=
      cres(i1,i2,i3,i5,1,5)=cres(i1,i2,i3,i5,1,5)+(l7_2516(i2,i5
     & ,i1)%a(1)*r8_34(i3)%a(1)+l7_2516(i2,i5,i1)%c(1)*s834*r8_3
     & 4(i3)%b(2))
      cres(i1,i2,i3,i5,2,5)=cres(i1,i2,i3,i5,2,5)+(l7_2516(i2,i5
     & ,i1)%c(2)*s834*r8_34(i3)%b(1)+l7_2516(i2,i5,i1)%a(2)*r8_3
     & 4(i3)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,5),a1=l7_34516(i3,i5,i1)%a,c1=l7_34516(
* i3,i5,i1)%c,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s82,bef=cres(i1,i2,i3,i5,&
* ,5)+,aft=
      cres(i1,i2,i3,i5,1,5)=cres(i1,i2,i3,i5,1,5)+(l7_34516(i3,i
     & 5,i1)%a(1)*r8_2(i2)%a(1)+l7_34516(i3,i5,i1)%c(1)*s82*r8_2
     & (i2)%b(2))
      cres(i1,i2,i3,i5,2,5)=cres(i1,i2,i3,i5,2,5)+(l7_34516(i3,i
     & 5,i1)%c(2)*s82*r8_2(i2)%b(1)+l7_34516(i3,i5,i1)%a(2)*r8_2
     & (i2)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,6),a1=l5_234(i2,i3)%a,c1=l5_234(i2,i3)%
* c,a2=r6_718(i7,i1)%a,b2=r6_718(i7,i1)%b,prq=s5234,bef=cres(i1,i2,i3,&,
* i7,6)+,aft=
      cres(i1,i2,i3,1,i7,6)=cres(i1,i2,i3,1,i7,6)+(l5_234(i2,i3)
     & %a(1)*r6_718(i7,i1)%a(1)+l5_234(i2,i3)%c(1)*s5234*r6_718(
     & i7,i1)%b(2))
      cres(i1,i2,i3,2,i7,6)=cres(i1,i2,i3,2,i7,6)+(l5_234(i2,i3)
     & %c(2)*s5234*r6_718(i7,i1)%b(1)+l5_234(i2,i3)%a(2)*r6_718(
     & i7,i1)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,6),a1=l5_2718(i2,i7,i1)%a,c1=l5_2718(i2
* ,i7,i1)%c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=cres(i1,i2,i3,&,i
* 7,6)+,aft=
      cres(i1,i2,i3,1,i7,6)=cres(i1,i2,i3,1,i7,6)+(l5_2718(i2,i7
     & ,i1)%a(1)*r6_34(i3)%a(1)+l5_2718(i2,i7,i1)%c(1)*s634*r6_3
     & 4(i3)%b(2))
      cres(i1,i2,i3,2,i7,6)=cres(i1,i2,i3,2,i7,6)+(l5_2718(i2,i7
     & ,i1)%c(2)*s634*r6_34(i3)%b(1)+l5_2718(i2,i7,i1)%a(2)*r6_3
     & 4(i3)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,6),a1=l5_34718(i3,i7,i1)%a,c1=l5_34718(
* i3,i7,i1)%c,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s62,bef=cres(i1,i2,i3,&,i7
* ,6)+,aft=
      cres(i1,i2,i3,1,i7,6)=cres(i1,i2,i3,1,i7,6)+(l5_34718(i3,i
     & 7,i1)%a(1)*r6_2(i2)%a(1)+l5_34718(i3,i7,i1)%c(1)*s62*r6_2
     & (i2)%b(2))
      cres(i1,i2,i3,2,i7,6)=cres(i1,i2,i3,2,i7,6)+(l5_34718(i3,i
     & 7,i1)%c(2)*s62*r6_2(i2)%b(1)+l5_34718(i3,i7,i1)%a(2)*r6_2
     & (i2)%a(2))
      end do
      end do
      end do
      end do
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,6),a1=l7_134(i1,i3)%a,c1=l7_134(i1,i3)%
* c,a2=r8_526(i5,i2)%a,b2=r8_526(i5,i2)%b,prq=s7134,bef=cres(i1,i2,i3,i5
* ,&,6)+,aft=
      cres(i1,i2,i3,i5,1,6)=cres(i1,i2,i3,i5,1,6)+(l7_134(i1,i3)
     & %a(1)*r8_526(i5,i2)%a(1)+l7_134(i1,i3)%c(1)*s7134*r8_526(
     & i5,i2)%b(2))
      cres(i1,i2,i3,i5,2,6)=cres(i1,i2,i3,i5,2,6)+(l7_134(i1,i3)
     & %c(2)*s7134*r8_526(i5,i2)%b(1)+l7_134(i1,i3)%a(2)*r8_526(
     & i5,i2)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,6),a1=l7_1526(i1,i5,i2)%a,c1=l7_1526(i1
* ,i5,i2)%c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=cres(i1,i2,i3,i5,
* &,6)+,aft=
      cres(i1,i2,i3,i5,1,6)=cres(i1,i2,i3,i5,1,6)+(l7_1526(i1,i5
     & ,i2)%a(1)*r8_34(i3)%a(1)+l7_1526(i1,i5,i2)%c(1)*s834*r8_3
     & 4(i3)%b(2))
      cres(i1,i2,i3,i5,2,6)=cres(i1,i2,i3,i5,2,6)+(l7_1526(i1,i5
     & ,i2)%c(2)*s834*r8_34(i3)%b(1)+l7_1526(i1,i5,i2)%a(2)*r8_3
     & 4(i3)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,6),a1=l7_34526(i3,i5,i2)%a,c1=l7_34526(
* i3,i5,i2)%c,a2=r8_1(i1)%a,b2=r8_1(i1)%b,prq=s81,bef=cres(i1,i2,i3,i5,&
* ,6)+,aft=
      cres(i1,i2,i3,i5,1,6)=cres(i1,i2,i3,i5,1,6)+(l7_34526(i3,i
     & 5,i2)%a(1)*r8_1(i1)%a(1)+l7_34526(i3,i5,i2)%c(1)*s81*r8_1
     & (i1)%b(2))
      cres(i1,i2,i3,i5,2,6)=cres(i1,i2,i3,i5,2,6)+(l7_34526(i3,i
     & 5,i2)%c(2)*s81*r8_1(i1)%b(1)+l7_34526(i3,i5,i2)%a(2)*r8_1
     & (i1)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
  
* -> lines with two gluons                                              
  
      if (ilept(id3).ne.1) then
* I                                                                     
* quqd -- p=p312,q=p478
      quqd=p312(0)*p478(0)-p312(1)*p478(1)-p312(2)*p478(2)-p312(
     & 3)*p478(3)
      ccr=zcr(id3)/(f478)
      ccl=zcl(id3)/(f478)
      do i5=1,2
* T0 -- qu=p312,qd=p478,v=cz56(i5)%e,a=u312_56(i5)%a,b=u312_56(i5)%b,c=u
* 312_56(i5)%c,d=u312_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p312(2)*p478(3)-p478(2)*p312(3))+p31
     & 2k0*(cz56(i5)%e(2)*p478(3)-p478(2)*cz56(i5)%e(3))-p478k0*
     & (cz56(i5)%e(2)*p312(3)-p312(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p312k0+p312(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p478k0+p478(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p312(0)-cz56(i5)%e(1)*p312(1)-cz56(i5)%
     & e(2)*p312(2)-cz56(i5)%e(3)*p312(3)
      cvqd=cz56(i5)%e(0)*p478(0)-cz56(i5)%e(1)*p478(1)-cz56(i5)%
     & e(2)*p478(2)-cz56(i5)%e(3)*p478(3)
      cauxa=-cz56(i5)%ek0*quqd+p312k0*cvqd+p478k0*cvqu
      cauxb=-cz56(i5)%ek0*p478(2)+p478k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p312(2)-p312k0*cz56(i5)%e(2)
      u312_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u312_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u312_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u312_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u312_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u312_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u312_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u312_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f478)
      ccl=fcl(id3)/(f478)
      do i5=1,2
* T0 -- qu=p312,qd=p478,v=cf56(i5)%e,a=u312_56(i5)%a,b=u312_56(i5)%b,c=u
* 312_56(i5)%c,d=u312_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p312(2)*p478(3)-p478(2)*p312(3))+p31
     & 2k0*(cf56(i5)%e(2)*p478(3)-p478(2)*cf56(i5)%e(3))-p478k0*
     & (cf56(i5)%e(2)*p312(3)-p312(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p312k0+p312(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p478k0+p478(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p312(0)-cf56(i5)%e(1)*p312(1)-cf56(i5)%
     & e(2)*p312(2)-cf56(i5)%e(3)*p312(3)
      cvqd=cf56(i5)%e(0)*p478(0)-cf56(i5)%e(1)*p478(1)-cf56(i5)%
     & e(2)*p478(2)-cf56(i5)%e(3)*p478(3)
      cauxa=-cf56(i5)%ek0*quqd+p312k0*cvqd+p478k0*cvqu
      cauxb=-cf56(i5)%ek0*p478(2)+p478k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p312(2)-p312k0*cf56(i5)%e(2)
      u312_56(i5)%a(1)=u312_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u312_56(i5)%a(2)=u312_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u312_56(i5)%b(1)=u312_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u312_56(i5)%b(2)=u312_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u312_56(i5)%c(1)=u312_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u312_56(i5)%c(2)=u312_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u312_56(i5)%d(1)=u312_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u312_56(i5)%d(2)=u312_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TTR0 -- aa=r4_5678(i5,i7)%a,bb=r4_5678(i5,i7)%b,a1=u312_56(i5)%a,b1=u3
* 12_56(i5)%b,c1=u312_56(i5)%c,d1=u312_56(i5)%d,a2=r4_78(i7)%a,b2=r4_78(
* i7)%b,prq=s478,nsum=0
      r4_5678(i5,i7)%a(1)=u312_56(i5)%a(1)*r4_78(i7)%a(1)+u312_5
     & 6(i5)%c(1)*s478*r4_78(i7)%b(2)
      r4_5678(i5,i7)%b(1)=u312_56(i5)%d(1)*s478*r4_78(i7)%b(1)+u
     & 312_56(i5)%b(1)*r4_78(i7)%a(2)
      r4_5678(i5,i7)%b(2)=u312_56(i5)%b(2)*r4_78(i7)%a(1)+u312_5
     & 6(i5)%d(2)*s478*r4_78(i7)%b(2)
      r4_5678(i5,i7)%a(2)=u312_56(i5)%c(2)*s478*r4_78(i7)%b(1)+u
     & 312_56(i5)%a(2)*r4_78(i7)%a(2)
      end do
      end do
  
  
* quqd -- p=p312,q=p456
      quqd=p312(0)*p456(0)-p312(1)*p456(1)-p312(2)*p456(2)-p312(
     & 3)*p456(3)
      ccr=zcr(id3)/(f456)
      ccl=zcl(id3)/(f456)
      do i7=1,2
* T0 -- qu=p312,qd=p456,v=cz78(i7)%e,a=u312_78(i7)%a,b=u312_78(i7)%b,c=u
* 312_78(i7)%c,d=u312_78(i7)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i7)%ek0*(p312(2)*p456(3)-p456(2)*p312(3))+p31
     & 2k0*(cz78(i7)%e(2)*p456(3)-p456(2)*cz78(i7)%e(3))-p456k0*
     & (cz78(i7)%e(2)*p312(3)-p312(2)*cz78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i7)%e(3)*p312k0+p312(3)*cz78(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i7)%e(3)*p456k0+p456(3)*cz78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i7)%e(0)*p312(0)-cz78(i7)%e(1)*p312(1)-cz78(i7)%
     & e(2)*p312(2)-cz78(i7)%e(3)*p312(3)
      cvqd=cz78(i7)%e(0)*p456(0)-cz78(i7)%e(1)*p456(1)-cz78(i7)%
     & e(2)*p456(2)-cz78(i7)%e(3)*p456(3)
      cauxa=-cz78(i7)%ek0*quqd+p312k0*cvqd+p456k0*cvqu
      cauxb=-cz78(i7)%ek0*p456(2)+p456k0*cz78(i7)%e(2)
      cauxc=+cz78(i7)%ek0*p312(2)-p312k0*cz78(i7)%e(2)
      u312_78(i7)%a(1)=ccr*(cauxa+ceps_0)
      u312_78(i7)%a(2)=ccl*(cauxa-ceps_0)
      u312_78(i7)%b(1)=ccl*(cauxb-ceps_2)
      u312_78(i7)%b(2)=ccr*(-cauxb-ceps_2)
      u312_78(i7)%c(1)=ccr*(cauxc+ceps_1)
      u312_78(i7)%c(2)=ccl*(-cauxc+ceps_1)
      u312_78(i7)%d(1)=ccl*cz78(i7)%ek0
      u312_78(i7)%d(2)=ccr*cz78(i7)%ek0
      end do
  
      ccr=fcr(id3)/(f456)
      ccl=fcl(id3)/(f456)
      do i7=1,2
* T0 -- qu=p312,qd=p456,v=cf78(i7)%e,a=u312_78(i7)%a,b=u312_78(i7)%b,c=u
* 312_78(i7)%c,d=u312_78(i7)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i7)%ek0*(p312(2)*p456(3)-p456(2)*p312(3))+p31
     & 2k0*(cf78(i7)%e(2)*p456(3)-p456(2)*cf78(i7)%e(3))-p456k0*
     & (cf78(i7)%e(2)*p312(3)-p312(2)*cf78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i7)%e(3)*p312k0+p312(3)*cf78(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i7)%e(3)*p456k0+p456(3)*cf78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i7)%e(0)*p312(0)-cf78(i7)%e(1)*p312(1)-cf78(i7)%
     & e(2)*p312(2)-cf78(i7)%e(3)*p312(3)
      cvqd=cf78(i7)%e(0)*p456(0)-cf78(i7)%e(1)*p456(1)-cf78(i7)%
     & e(2)*p456(2)-cf78(i7)%e(3)*p456(3)
      cauxa=-cf78(i7)%ek0*quqd+p312k0*cvqd+p456k0*cvqu
      cauxb=-cf78(i7)%ek0*p456(2)+p456k0*cf78(i7)%e(2)
      cauxc=+cf78(i7)%ek0*p312(2)-p312k0*cf78(i7)%e(2)
      u312_78(i7)%a(1)=u312_78(i7)%a(1)+ccr*(cauxa+ceps_0)
      u312_78(i7)%a(2)=u312_78(i7)%a(2)+ccl*(cauxa-ceps_0)
      u312_78(i7)%b(1)=u312_78(i7)%b(1)+ccl*(cauxb-ceps_2)
      u312_78(i7)%b(2)=u312_78(i7)%b(2)+ccr*(-cauxb-ceps_2)
      u312_78(i7)%c(1)=u312_78(i7)%c(1)+ccr*(cauxc+ceps_1)
      u312_78(i7)%c(2)=u312_78(i7)%c(2)+ccl*(-cauxc+ceps_1)
      u312_78(i7)%d(1)=u312_78(i7)%d(1)+ccl*cf78(i7)%ek0
      u312_78(i7)%d(2)=u312_78(i7)%d(2)+ccr*cf78(i7)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TTR0 -- aa=r4_5678(i5,i7)%a,bb=r4_5678(i5,i7)%b,a1=u312_78(i7)%a,b1=u3
* 12_78(i7)%b,c1=u312_78(i7)%c,d1=u312_78(i7)%d,a2=r4_56(i5)%a,b2=r4_56(
* i5)%b,prq=s456,nsum=1
      r4_5678(i5,i7)%a(1)=r4_5678(i5,i7)%a(1)+u312_78(i7)%a(1)*r
     & 4_56(i5)%a(1)+u312_78(i7)%c(1)*s456*r4_56(i5)%b(2)
      r4_5678(i5,i7)%b(1)=r4_5678(i5,i7)%b(1)+u312_78(i7)%d(1)*s
     & 456*r4_56(i5)%b(1)+u312_78(i7)%b(1)*r4_56(i5)%a(2)
      r4_5678(i5,i7)%b(2)=r4_5678(i5,i7)%b(2)+u312_78(i7)%b(2)*r
     & 4_56(i5)%a(1)+u312_78(i7)%d(2)*s456*r4_56(i5)%b(2)
      r4_5678(i5,i7)%a(2)=r4_5678(i5,i7)%a(2)+u312_78(i7)%c(2)*s
     & 456*r4_56(i5)%b(1)+u312_78(i7)%a(2)*r4_56(i5)%a(2)
      end do
      end do
  
* II                                                                    
* quqd -- p=p356,q=p412
      quqd=p356(0)*p412(0)-p356(1)*p412(1)-p356(2)*p412(2)-p356(
     & 3)*p412(3)
      ccr=zcr(id3)/(f412)
      ccl=zcl(id3)/(f412)
      do i7=1,2
* T0 -- qu=p356,qd=p412,v=cz78(i7)%e,a=u356_78(i7)%a,b=u356_78(i7)%b,c=u
* 356_78(i7)%c,d=u356_78(i7)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i7)%ek0*(p356(2)*p412(3)-p412(2)*p356(3))+p35
     & 6k0*(cz78(i7)%e(2)*p412(3)-p412(2)*cz78(i7)%e(3))-p412k0*
     & (cz78(i7)%e(2)*p356(3)-p356(2)*cz78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i7)%e(3)*p356k0+p356(3)*cz78(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i7)%e(3)*p412k0+p412(3)*cz78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i7)%e(0)*p356(0)-cz78(i7)%e(1)*p356(1)-cz78(i7)%
     & e(2)*p356(2)-cz78(i7)%e(3)*p356(3)
      cvqd=cz78(i7)%e(0)*p412(0)-cz78(i7)%e(1)*p412(1)-cz78(i7)%
     & e(2)*p412(2)-cz78(i7)%e(3)*p412(3)
      cauxa=-cz78(i7)%ek0*quqd+p356k0*cvqd+p412k0*cvqu
      cauxb=-cz78(i7)%ek0*p412(2)+p412k0*cz78(i7)%e(2)
      cauxc=+cz78(i7)%ek0*p356(2)-p356k0*cz78(i7)%e(2)
      u356_78(i7)%a(1)=ccr*(cauxa+ceps_0)
      u356_78(i7)%a(2)=ccl*(cauxa-ceps_0)
      u356_78(i7)%b(1)=ccl*(cauxb-ceps_2)
      u356_78(i7)%b(2)=ccr*(-cauxb-ceps_2)
      u356_78(i7)%c(1)=ccr*(cauxc+ceps_1)
      u356_78(i7)%c(2)=ccl*(-cauxc+ceps_1)
      u356_78(i7)%d(1)=ccl*cz78(i7)%ek0
      u356_78(i7)%d(2)=ccr*cz78(i7)%ek0
      end do
  
      ccr=fcr(id3)/(f412)
      ccl=fcl(id3)/(f412)
      do i7=1,2
* T0 -- qu=p356,qd=p412,v=cf78(i7)%e,a=u356_78(i7)%a,b=u356_78(i7)%b,c=u
* 356_78(i7)%c,d=u356_78(i7)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i7)%ek0*(p356(2)*p412(3)-p412(2)*p356(3))+p35
     & 6k0*(cf78(i7)%e(2)*p412(3)-p412(2)*cf78(i7)%e(3))-p412k0*
     & (cf78(i7)%e(2)*p356(3)-p356(2)*cf78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i7)%e(3)*p356k0+p356(3)*cf78(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i7)%e(3)*p412k0+p412(3)*cf78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i7)%e(0)*p356(0)-cf78(i7)%e(1)*p356(1)-cf78(i7)%
     & e(2)*p356(2)-cf78(i7)%e(3)*p356(3)
      cvqd=cf78(i7)%e(0)*p412(0)-cf78(i7)%e(1)*p412(1)-cf78(i7)%
     & e(2)*p412(2)-cf78(i7)%e(3)*p412(3)
      cauxa=-cf78(i7)%ek0*quqd+p356k0*cvqd+p412k0*cvqu
      cauxb=-cf78(i7)%ek0*p412(2)+p412k0*cf78(i7)%e(2)
      cauxc=+cf78(i7)%ek0*p356(2)-p356k0*cf78(i7)%e(2)
      u356_78(i7)%a(1)=u356_78(i7)%a(1)+ccr*(cauxa+ceps_0)
      u356_78(i7)%a(2)=u356_78(i7)%a(2)+ccl*(cauxa-ceps_0)
      u356_78(i7)%b(1)=u356_78(i7)%b(1)+ccl*(cauxb-ceps_2)
      u356_78(i7)%b(2)=u356_78(i7)%b(2)+ccr*(-cauxb-ceps_2)
      u356_78(i7)%c(1)=u356_78(i7)%c(1)+ccr*(cauxc+ceps_1)
      u356_78(i7)%c(2)=u356_78(i7)%c(2)+ccl*(-cauxc+ceps_1)
      u356_78(i7)%d(1)=u356_78(i7)%d(1)+ccl*cf78(i7)%ek0
      u356_78(i7)%d(2)=u356_78(i7)%d(2)+ccr*cf78(i7)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TLT0 -- aa=l3_5678(i5,i7)%a,cc=l3_5678(i5,i7)%c,a1=l3_56(i5)%a,c1=l3_5
* 6(i5)%c,a2=u356_78(i7)%a,b2=u356_78(i7)%b,c2=u356_78(i7)%c,d2=u356_78(
* i7)%d,prq=s356,nsum=0
      l3_5678(i5,i7)%a(1)=l3_56(i5)%a(1)*u356_78(i7)%a(1)+l3_56(
     & i5)%c(1)*s356*u356_78(i7)%b(2)
      l3_5678(i5,i7)%c(1)=l3_56(i5)%a(1)*u356_78(i7)%c(1)+l3_56(
     & i5)%c(1)*s356*u356_78(i7)%d(2)
      l3_5678(i5,i7)%c(2)=l3_56(i5)%c(2)*s356*u356_78(i7)%d(1)+l
     & 3_56(i5)%a(2)*u356_78(i7)%c(2)
      l3_5678(i5,i7)%a(2)=l3_56(i5)%c(2)*s356*u356_78(i7)%b(1)+l
     & 3_56(i5)%a(2)*u356_78(i7)%a(2)
      end do
      end do
  
  
* quqd -- p=p378,q=p412
      quqd=p378(0)*p412(0)-p378(1)*p412(1)-p378(2)*p412(2)-p378(
     & 3)*p412(3)
      ccr=zcr(id3)/(f412)
      ccl=zcl(id3)/(f412)
      do i5=1,2
* T0 -- qu=p378,qd=p412,v=cz56(i5)%e,a=u378_56(i5)%a,b=u378_56(i5)%b,c=u
* 378_56(i5)%c,d=u378_56(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i5)%ek0*(p378(2)*p412(3)-p412(2)*p378(3))+p37
     & 8k0*(cz56(i5)%e(2)*p412(3)-p412(2)*cz56(i5)%e(3))-p412k0*
     & (cz56(i5)%e(2)*p378(3)-p378(2)*cz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i5)%e(3)*p378k0+p378(3)*cz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i5)%e(3)*p412k0+p412(3)*cz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i5)%e(0)*p378(0)-cz56(i5)%e(1)*p378(1)-cz56(i5)%
     & e(2)*p378(2)-cz56(i5)%e(3)*p378(3)
      cvqd=cz56(i5)%e(0)*p412(0)-cz56(i5)%e(1)*p412(1)-cz56(i5)%
     & e(2)*p412(2)-cz56(i5)%e(3)*p412(3)
      cauxa=-cz56(i5)%ek0*quqd+p378k0*cvqd+p412k0*cvqu
      cauxb=-cz56(i5)%ek0*p412(2)+p412k0*cz56(i5)%e(2)
      cauxc=+cz56(i5)%ek0*p378(2)-p378k0*cz56(i5)%e(2)
      u378_56(i5)%a(1)=ccr*(cauxa+ceps_0)
      u378_56(i5)%a(2)=ccl*(cauxa-ceps_0)
      u378_56(i5)%b(1)=ccl*(cauxb-ceps_2)
      u378_56(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u378_56(i5)%c(1)=ccr*(cauxc+ceps_1)
      u378_56(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u378_56(i5)%d(1)=ccl*cz56(i5)%ek0
      u378_56(i5)%d(2)=ccr*cz56(i5)%ek0
      end do
  
      ccr=fcr(id3)/(f412)
      ccl=fcl(id3)/(f412)
      do i5=1,2
* T0 -- qu=p378,qd=p412,v=cf56(i5)%e,a=u378_56(i5)%a,b=u378_56(i5)%b,c=u
* 378_56(i5)%c,d=u378_56(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i5)%ek0*(p378(2)*p412(3)-p412(2)*p378(3))+p37
     & 8k0*(cf56(i5)%e(2)*p412(3)-p412(2)*cf56(i5)%e(3))-p412k0*
     & (cf56(i5)%e(2)*p378(3)-p378(2)*cf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i5)%e(3)*p378k0+p378(3)*cf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i5)%e(3)*p412k0+p412(3)*cf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i5)%e(0)*p378(0)-cf56(i5)%e(1)*p378(1)-cf56(i5)%
     & e(2)*p378(2)-cf56(i5)%e(3)*p378(3)
      cvqd=cf56(i5)%e(0)*p412(0)-cf56(i5)%e(1)*p412(1)-cf56(i5)%
     & e(2)*p412(2)-cf56(i5)%e(3)*p412(3)
      cauxa=-cf56(i5)%ek0*quqd+p378k0*cvqd+p412k0*cvqu
      cauxb=-cf56(i5)%ek0*p412(2)+p412k0*cf56(i5)%e(2)
      cauxc=+cf56(i5)%ek0*p378(2)-p378k0*cf56(i5)%e(2)
      u378_56(i5)%a(1)=u378_56(i5)%a(1)+ccr*(cauxa+ceps_0)
      u378_56(i5)%a(2)=u378_56(i5)%a(2)+ccl*(cauxa-ceps_0)
      u378_56(i5)%b(1)=u378_56(i5)%b(1)+ccl*(cauxb-ceps_2)
      u378_56(i5)%b(2)=u378_56(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u378_56(i5)%c(1)=u378_56(i5)%c(1)+ccr*(cauxc+ceps_1)
      u378_56(i5)%c(2)=u378_56(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u378_56(i5)%d(1)=u378_56(i5)%d(1)+ccl*cf56(i5)%ek0
      u378_56(i5)%d(2)=u378_56(i5)%d(2)+ccr*cf56(i5)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TLT0 -- aa=l3_5678(i5,i7)%a,cc=l3_5678(i5,i7)%c,a1=l3_78(i7)%a,c1=l3_7
* 8(i7)%c,a2=u378_56(i5)%a,b2=u378_56(i5)%b,c2=u378_56(i5)%c,d2=u378_56(
* i5)%d,prq=s378,nsum=1
      l3_5678(i5,i7)%a(1)=l3_5678(i5,i7)%a(1)+l3_78(i7)%a(1)*u37
     & 8_56(i5)%a(1)+l3_78(i7)%c(1)*s378*u378_56(i5)%b(2)
      l3_5678(i5,i7)%c(1)=l3_5678(i5,i7)%c(1)+l3_78(i7)%a(1)*u37
     & 8_56(i5)%c(1)+l3_78(i7)%c(1)*s378*u378_56(i5)%d(2)
      l3_5678(i5,i7)%c(2)=l3_5678(i5,i7)%c(2)+l3_78(i7)%c(2)*s37
     & 8*u378_56(i5)%d(1)+l3_78(i7)%a(2)*u378_56(i5)%c(2)
      l3_5678(i5,i7)%a(2)=l3_5678(i5,i7)%a(2)+l3_78(i7)%c(2)*s37
     & 8*u378_56(i5)%b(1)+l3_78(i7)%a(2)*u378_56(i5)%a(2)
      end do
      end do
  
* III                                                                   
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r4_278(i2,i7)%a,bb=r4_278(i2,i7)%b,a1=u3516_2(i2)%a,b1=u351
* 6_2(i2)%b,c1=u3516_2(i2)%c,d1=u3516_2(i2)%d,a2=r4_78(i7)%a,b2=r4_78(i7
* )%b,prq=s478,nsum=0
      r4_278(i2,i7)%a(1)=u3516_2(i2)%a(1)*r4_78(i7)%a(1)+u3516_2
     & (i2)%c(1)*s478*r4_78(i7)%b(2)
      r4_278(i2,i7)%b(1)=u3516_2(i2)%d(1)*s478*r4_78(i7)%b(1)+u3
     & 516_2(i2)%b(1)*r4_78(i7)%a(2)
      r4_278(i2,i7)%b(2)=u3516_2(i2)%b(2)*r4_78(i7)%a(1)+u3516_2
     & (i2)%d(2)*s478*r4_78(i7)%b(2)
      r4_278(i2,i7)%a(2)=u3516_2(i2)%c(2)*s478*r4_78(i7)%b(1)+u3
     & 516_2(i2)%a(2)*r4_78(i7)%a(2)
      end do
      end do
  
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r4_278(i2,i7)%a,bb=r4_278(i2,i7)%b,a1=u3516_78(i7)%a,b1=u35
* 16_78(i7)%b,c1=u3516_78(i7)%c,d1=u3516_78(i7)%d,a2=r4_2(i2)%a,b2=r4_2(
* i2)%b,prq=s42,nsum=1
      r4_278(i2,i7)%a(1)=r4_278(i2,i7)%a(1)+u3516_78(i7)%a(1)*r4
     & _2(i2)%a(1)+u3516_78(i7)%c(1)*s42*r4_2(i2)%b(2)
      r4_278(i2,i7)%b(1)=r4_278(i2,i7)%b(1)+u3516_78(i7)%d(1)*s4
     & 2*r4_2(i2)%b(1)+u3516_78(i7)%b(1)*r4_2(i2)%a(2)
      r4_278(i2,i7)%b(2)=r4_278(i2,i7)%b(2)+u3516_78(i7)%b(2)*r4
     & _2(i2)%a(1)+u3516_78(i7)%d(2)*s42*r4_2(i2)%b(2)
      r4_278(i2,i7)%a(2)=r4_278(i2,i7)%a(2)+u3516_78(i7)%c(2)*s4
     & 2*r4_2(i2)%b(1)+u3516_78(i7)%a(2)*r4_2(i2)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r4_256(i2,i5)%a,bb=r4_256(i2,i5)%b,a1=u3718_2(i2)%a,b1=u371
* 8_2(i2)%b,c1=u3718_2(i2)%c,d1=u3718_2(i2)%d,a2=r4_56(i5)%a,b2=r4_56(i5
* )%b,prq=s456,nsum=0
      r4_256(i2,i5)%a(1)=u3718_2(i2)%a(1)*r4_56(i5)%a(1)+u3718_2
     & (i2)%c(1)*s456*r4_56(i5)%b(2)
      r4_256(i2,i5)%b(1)=u3718_2(i2)%d(1)*s456*r4_56(i5)%b(1)+u3
     & 718_2(i2)%b(1)*r4_56(i5)%a(2)
      r4_256(i2,i5)%b(2)=u3718_2(i2)%b(2)*r4_56(i5)%a(1)+u3718_2
     & (i2)%d(2)*s456*r4_56(i5)%b(2)
      r4_256(i2,i5)%a(2)=u3718_2(i2)%c(2)*s456*r4_56(i5)%b(1)+u3
     & 718_2(i2)%a(2)*r4_56(i5)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r4_256(i2,i5)%a,bb=r4_256(i2,i5)%b,a1=u3718_56(i5)%a,b1=u37
* 18_56(i5)%b,c1=u3718_56(i5)%c,d1=u3718_56(i5)%d,a2=r4_2(i2)%a,b2=r4_2(
* i2)%b,prq=s42,nsum=1
      r4_256(i2,i5)%a(1)=r4_256(i2,i5)%a(1)+u3718_56(i5)%a(1)*r4
     & _2(i2)%a(1)+u3718_56(i5)%c(1)*s42*r4_2(i2)%b(2)
      r4_256(i2,i5)%b(1)=r4_256(i2,i5)%b(1)+u3718_56(i5)%d(1)*s4
     & 2*r4_2(i2)%b(1)+u3718_56(i5)%b(1)*r4_2(i2)%a(2)
      r4_256(i2,i5)%b(2)=r4_256(i2,i5)%b(2)+u3718_56(i5)%b(2)*r4
     & _2(i2)%a(1)+u3718_56(i5)%d(2)*s42*r4_2(i2)%b(2)
      r4_256(i2,i5)%a(2)=r4_256(i2,i5)%a(2)+u3718_56(i5)%c(2)*s4
     & 2*r4_2(i2)%b(1)+u3718_56(i5)%a(2)*r4_2(i2)%a(2)
      end do
      end do
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r4_178(i2,i7)%a,bb=r4_178(i2,i7)%b,a1=u3526_1(i2)%a,b1=u352
* 6_1(i2)%b,c1=u3526_1(i2)%c,d1=u3526_1(i2)%d,a2=r4_78(i7)%a,b2=r4_78(i7
* )%b,prq=s478,nsum=0
      r4_178(i2,i7)%a(1)=u3526_1(i2)%a(1)*r4_78(i7)%a(1)+u3526_1
     & (i2)%c(1)*s478*r4_78(i7)%b(2)
      r4_178(i2,i7)%b(1)=u3526_1(i2)%d(1)*s478*r4_78(i7)%b(1)+u3
     & 526_1(i2)%b(1)*r4_78(i7)%a(2)
      r4_178(i2,i7)%b(2)=u3526_1(i2)%b(2)*r4_78(i7)%a(1)+u3526_1
     & (i2)%d(2)*s478*r4_78(i7)%b(2)
      r4_178(i2,i7)%a(2)=u3526_1(i2)%c(2)*s478*r4_78(i7)%b(1)+u3
     & 526_1(i2)%a(2)*r4_78(i7)%a(2)
      end do
      end do
  
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r4_178(i2,i7)%a,bb=r4_178(i2,i7)%b,a1=u3526_78(i7)%a,b1=u35
* 26_78(i7)%b,c1=u3526_78(i7)%c,d1=u3526_78(i7)%d,a2=r4_1(i2)%a,b2=r4_1(
* i2)%b,prq=s41,nsum=1
      r4_178(i2,i7)%a(1)=r4_178(i2,i7)%a(1)+u3526_78(i7)%a(1)*r4
     & _1(i2)%a(1)+u3526_78(i7)%c(1)*s41*r4_1(i2)%b(2)
      r4_178(i2,i7)%b(1)=r4_178(i2,i7)%b(1)+u3526_78(i7)%d(1)*s4
     & 1*r4_1(i2)%b(1)+u3526_78(i7)%b(1)*r4_1(i2)%a(2)
      r4_178(i2,i7)%b(2)=r4_178(i2,i7)%b(2)+u3526_78(i7)%b(2)*r4
     & _1(i2)%a(1)+u3526_78(i7)%d(2)*s41*r4_1(i2)%b(2)
      r4_178(i2,i7)%a(2)=r4_178(i2,i7)%a(2)+u3526_78(i7)%c(2)*s4
     & 1*r4_1(i2)%b(1)+u3526_78(i7)%a(2)*r4_1(i2)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r4_156(i2,i5)%a,bb=r4_156(i2,i5)%b,a1=u3728_1(i2)%a,b1=u372
* 8_1(i2)%b,c1=u3728_1(i2)%c,d1=u3728_1(i2)%d,a2=r4_56(i5)%a,b2=r4_56(i5
* )%b,prq=s456,nsum=0
      r4_156(i2,i5)%a(1)=u3728_1(i2)%a(1)*r4_56(i5)%a(1)+u3728_1
     & (i2)%c(1)*s456*r4_56(i5)%b(2)
      r4_156(i2,i5)%b(1)=u3728_1(i2)%d(1)*s456*r4_56(i5)%b(1)+u3
     & 728_1(i2)%b(1)*r4_56(i5)%a(2)
      r4_156(i2,i5)%b(2)=u3728_1(i2)%b(2)*r4_56(i5)%a(1)+u3728_1
     & (i2)%d(2)*s456*r4_56(i5)%b(2)
      r4_156(i2,i5)%a(2)=u3728_1(i2)%c(2)*s456*r4_56(i5)%b(1)+u3
     & 728_1(i2)%a(2)*r4_56(i5)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r4_156(i2,i5)%a,bb=r4_156(i2,i5)%b,a1=u3728_56(i5)%a,b1=u37
* 28_56(i5)%b,c1=u3728_56(i5)%c,d1=u3728_56(i5)%d,a2=r4_1(i2)%a,b2=r4_1(
* i2)%b,prq=s41,nsum=1
      r4_156(i2,i5)%a(1)=r4_156(i2,i5)%a(1)+u3728_56(i5)%a(1)*r4
     & _1(i2)%a(1)+u3728_56(i5)%c(1)*s41*r4_1(i2)%b(2)
      r4_156(i2,i5)%b(1)=r4_156(i2,i5)%b(1)+u3728_56(i5)%d(1)*s4
     & 1*r4_1(i2)%b(1)+u3728_56(i5)%b(1)*r4_1(i2)%a(2)
      r4_156(i2,i5)%b(2)=r4_156(i2,i5)%b(2)+u3728_56(i5)%b(2)*r4
     & _1(i2)%a(1)+u3728_56(i5)%d(2)*s41*r4_1(i2)%b(2)
      r4_156(i2,i5)%a(2)=r4_156(i2,i5)%a(2)+u3728_56(i5)%c(2)*s4
     & 1*r4_1(i2)%b(1)+u3728_56(i5)%a(2)*r4_1(i2)%a(2)
      end do
      end do
  
* IV                                                                    
* quqd -- p=p356,q=p478
      quqd=p356(0)*p478(0)-p356(1)*p478(1)-p356(2)*p478(2)-p356(
     & 3)*p478(3)
      ccr=1d0/(f478*s12)
      ccl=1d0/(f478*s12)
      do i1=1,2
      do i2=1,2
* T0 -- qu=p356,qd=p478,v=ctrip12(i1,i2)%e,a=u356_gg(i1,i2)%a,b=u356_gg(
* i1,i2)%b,c=u356_gg(i1,i2)%c,d=u356_gg(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p356(2)*p478(3)-p478(2)*p356(3
     & ))+p356k0*(ctrip12(i1,i2)%e(2)*p478(3)-p478(2)*ctrip12(i1
     & ,i2)%e(3))-p478k0*(ctrip12(i1,i2)%e(2)*p356(3)-p356(2)*ct
     & rip12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ctrip12(i1,i2)%e(3)*p356k0+p356(3)*ctrip12(i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      ceps_2=-ctrip12(i1,i2)%e(3)*p478k0+p478(3)*ctrip12(i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=ctrip12(i1,i2)%e(0)*p356(0)-ctrip12(i1,i2)%e(1)*p356(
     & 1)-ctrip12(i1,i2)%e(2)*p356(2)-ctrip12(i1,i2)%e(3)*p356(3
     & )
      cvqd=ctrip12(i1,i2)%e(0)*p478(0)-ctrip12(i1,i2)%e(1)*p478(
     & 1)-ctrip12(i1,i2)%e(2)*p478(2)-ctrip12(i1,i2)%e(3)*p478(3
     & )
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p356k0*cvqd+p478k0*cvqu
      cauxb=-ctrip12(i1,i2)%ek0*p478(2)+p478k0*ctrip12(i1,i2)%e(
     & 2)
      cauxc=+ctrip12(i1,i2)%ek0*p356(2)-p356k0*ctrip12(i1,i2)%e(
     & 2)
      u356_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      u356_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u356_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u356_gg(i1,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u356_gg(i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      u356_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u356_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      u356_gg(i1,i2)%d(2)=ccr*ctrip12(i1,i2)%ek0
      end do
      end do
  
* quqd -- p=p378,q=p456
      quqd=p378(0)*p456(0)-p378(1)*p456(1)-p378(2)*p456(2)-p378(
     & 3)*p456(3)
      ccr=1d0/(f456*s12)
      ccl=1d0/(f456*s12)
      do i1=1,2
      do i2=1,2
* T0 -- qu=p378,qd=p456,v=ctrip12(i1,i2)%e,a=u378_gg(i1,i2)%a,b=u378_gg(
* i1,i2)%b,c=u378_gg(i1,i2)%c,d=u378_gg(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p378(2)*p456(3)-p456(2)*p378(3
     & ))+p378k0*(ctrip12(i1,i2)%e(2)*p456(3)-p456(2)*ctrip12(i1
     & ,i2)%e(3))-p456k0*(ctrip12(i1,i2)%e(2)*p378(3)-p378(2)*ct
     & rip12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ctrip12(i1,i2)%e(3)*p378k0+p378(3)*ctrip12(i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      ceps_2=-ctrip12(i1,i2)%e(3)*p456k0+p456(3)*ctrip12(i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=ctrip12(i1,i2)%e(0)*p378(0)-ctrip12(i1,i2)%e(1)*p378(
     & 1)-ctrip12(i1,i2)%e(2)*p378(2)-ctrip12(i1,i2)%e(3)*p378(3
     & )
      cvqd=ctrip12(i1,i2)%e(0)*p456(0)-ctrip12(i1,i2)%e(1)*p456(
     & 1)-ctrip12(i1,i2)%e(2)*p456(2)-ctrip12(i1,i2)%e(3)*p456(3
     & )
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p378k0*cvqd+p456k0*cvqu
      cauxb=-ctrip12(i1,i2)%ek0*p456(2)+p456k0*ctrip12(i1,i2)%e(
     & 2)
      cauxc=+ctrip12(i1,i2)%ek0*p378(2)-p378k0*ctrip12(i1,i2)%e(
     & 2)
      u378_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      u378_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u378_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u378_gg(i1,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u378_gg(i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      u378_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u378_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      u378_gg(i1,i2)%d(2)=ccr*ctrip12(i1,i2)%ek0
      end do
      end do
  
      endif
  
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p512,q=p678
      quqd=p512(0)*p678(0)-p512(1)*p678(1)-p512(2)*p678(2)-p512(
     & 3)*p678(3)
      ccr=zcr(id5)/(f678)
      ccl=zcl(id5)/(f678)
      do i5=1,2
* T0 -- qu=p512,qd=p678,v=cz34(i5)%e,a=u512_34(i5)%a,b=u512_34(i5)%b,c=u
* 512_34(i5)%c,d=u512_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p512(2)*p678(3)-p678(2)*p512(3))+p51
     & 2k0*(cz34(i5)%e(2)*p678(3)-p678(2)*cz34(i5)%e(3))-p678k0*
     & (cz34(i5)%e(2)*p512(3)-p512(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p512k0+p512(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p678k0+p678(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p512(0)-cz34(i5)%e(1)*p512(1)-cz34(i5)%
     & e(2)*p512(2)-cz34(i5)%e(3)*p512(3)
      cvqd=cz34(i5)%e(0)*p678(0)-cz34(i5)%e(1)*p678(1)-cz34(i5)%
     & e(2)*p678(2)-cz34(i5)%e(3)*p678(3)
      cauxa=-cz34(i5)%ek0*quqd+p512k0*cvqd+p678k0*cvqu
      cauxb=-cz34(i5)%ek0*p678(2)+p678k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p512(2)-p512k0*cz34(i5)%e(2)
      u512_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u512_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u512_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u512_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u512_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u512_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u512_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u512_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f678)
      ccl=fcl(id5)/(f678)
      do i5=1,2
* T0 -- qu=p512,qd=p678,v=cf34(i5)%e,a=u512_34(i5)%a,b=u512_34(i5)%b,c=u
* 512_34(i5)%c,d=u512_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p512(2)*p678(3)-p678(2)*p512(3))+p51
     & 2k0*(cf34(i5)%e(2)*p678(3)-p678(2)*cf34(i5)%e(3))-p678k0*
     & (cf34(i5)%e(2)*p512(3)-p512(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p512k0+p512(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p678k0+p678(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p512(0)-cf34(i5)%e(1)*p512(1)-cf34(i5)%
     & e(2)*p512(2)-cf34(i5)%e(3)*p512(3)
      cvqd=cf34(i5)%e(0)*p678(0)-cf34(i5)%e(1)*p678(1)-cf34(i5)%
     & e(2)*p678(2)-cf34(i5)%e(3)*p678(3)
      cauxa=-cf34(i5)%ek0*quqd+p512k0*cvqd+p678k0*cvqu
      cauxb=-cf34(i5)%ek0*p678(2)+p678k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p512(2)-p512k0*cf34(i5)%e(2)
      u512_34(i5)%a(1)=u512_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u512_34(i5)%a(2)=u512_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u512_34(i5)%b(1)=u512_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u512_34(i5)%b(2)=u512_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u512_34(i5)%c(1)=u512_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u512_34(i5)%c(2)=u512_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u512_34(i5)%d(1)=u512_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u512_34(i5)%d(2)=u512_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TTR0 -- aa=r6_3478(i5,i7)%a,bb=r6_3478(i5,i7)%b,a1=u512_34(i5)%a,b1=u5
* 12_34(i5)%b,c1=u512_34(i5)%c,d1=u512_34(i5)%d,a2=r6_78(i7)%a,b2=r6_78(
* i7)%b,prq=s678,nsum=0
      r6_3478(i5,i7)%a(1)=u512_34(i5)%a(1)*r6_78(i7)%a(1)+u512_3
     & 4(i5)%c(1)*s678*r6_78(i7)%b(2)
      r6_3478(i5,i7)%b(1)=u512_34(i5)%d(1)*s678*r6_78(i7)%b(1)+u
     & 512_34(i5)%b(1)*r6_78(i7)%a(2)
      r6_3478(i5,i7)%b(2)=u512_34(i5)%b(2)*r6_78(i7)%a(1)+u512_3
     & 4(i5)%d(2)*s678*r6_78(i7)%b(2)
      r6_3478(i5,i7)%a(2)=u512_34(i5)%c(2)*s678*r6_78(i7)%b(1)+u
     & 512_34(i5)%a(2)*r6_78(i7)%a(2)
      end do
      end do
  
  
* quqd -- p=p512,q=p634
      quqd=p512(0)*p634(0)-p512(1)*p634(1)-p512(2)*p634(2)-p512(
     & 3)*p634(3)
      ccr=zcr(id5)/(f634)
      ccl=zcl(id5)/(f634)
      do i7=1,2
* T0 -- qu=p512,qd=p634,v=cz78(i7)%e,a=u512_78(i7)%a,b=u512_78(i7)%b,c=u
* 512_78(i7)%c,d=u512_78(i7)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i7)%ek0*(p512(2)*p634(3)-p634(2)*p512(3))+p51
     & 2k0*(cz78(i7)%e(2)*p634(3)-p634(2)*cz78(i7)%e(3))-p634k0*
     & (cz78(i7)%e(2)*p512(3)-p512(2)*cz78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i7)%e(3)*p512k0+p512(3)*cz78(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i7)%e(3)*p634k0+p634(3)*cz78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i7)%e(0)*p512(0)-cz78(i7)%e(1)*p512(1)-cz78(i7)%
     & e(2)*p512(2)-cz78(i7)%e(3)*p512(3)
      cvqd=cz78(i7)%e(0)*p634(0)-cz78(i7)%e(1)*p634(1)-cz78(i7)%
     & e(2)*p634(2)-cz78(i7)%e(3)*p634(3)
      cauxa=-cz78(i7)%ek0*quqd+p512k0*cvqd+p634k0*cvqu
      cauxb=-cz78(i7)%ek0*p634(2)+p634k0*cz78(i7)%e(2)
      cauxc=+cz78(i7)%ek0*p512(2)-p512k0*cz78(i7)%e(2)
      u512_78(i7)%a(1)=ccr*(cauxa+ceps_0)
      u512_78(i7)%a(2)=ccl*(cauxa-ceps_0)
      u512_78(i7)%b(1)=ccl*(cauxb-ceps_2)
      u512_78(i7)%b(2)=ccr*(-cauxb-ceps_2)
      u512_78(i7)%c(1)=ccr*(cauxc+ceps_1)
      u512_78(i7)%c(2)=ccl*(-cauxc+ceps_1)
      u512_78(i7)%d(1)=ccl*cz78(i7)%ek0
      u512_78(i7)%d(2)=ccr*cz78(i7)%ek0
      end do
  
      ccr=fcr(id5)/(f634)
      ccl=fcl(id5)/(f634)
      do i7=1,2
* T0 -- qu=p512,qd=p634,v=cf78(i7)%e,a=u512_78(i7)%a,b=u512_78(i7)%b,c=u
* 512_78(i7)%c,d=u512_78(i7)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i7)%ek0*(p512(2)*p634(3)-p634(2)*p512(3))+p51
     & 2k0*(cf78(i7)%e(2)*p634(3)-p634(2)*cf78(i7)%e(3))-p634k0*
     & (cf78(i7)%e(2)*p512(3)-p512(2)*cf78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i7)%e(3)*p512k0+p512(3)*cf78(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i7)%e(3)*p634k0+p634(3)*cf78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i7)%e(0)*p512(0)-cf78(i7)%e(1)*p512(1)-cf78(i7)%
     & e(2)*p512(2)-cf78(i7)%e(3)*p512(3)
      cvqd=cf78(i7)%e(0)*p634(0)-cf78(i7)%e(1)*p634(1)-cf78(i7)%
     & e(2)*p634(2)-cf78(i7)%e(3)*p634(3)
      cauxa=-cf78(i7)%ek0*quqd+p512k0*cvqd+p634k0*cvqu
      cauxb=-cf78(i7)%ek0*p634(2)+p634k0*cf78(i7)%e(2)
      cauxc=+cf78(i7)%ek0*p512(2)-p512k0*cf78(i7)%e(2)
      u512_78(i7)%a(1)=u512_78(i7)%a(1)+ccr*(cauxa+ceps_0)
      u512_78(i7)%a(2)=u512_78(i7)%a(2)+ccl*(cauxa-ceps_0)
      u512_78(i7)%b(1)=u512_78(i7)%b(1)+ccl*(cauxb-ceps_2)
      u512_78(i7)%b(2)=u512_78(i7)%b(2)+ccr*(-cauxb-ceps_2)
      u512_78(i7)%c(1)=u512_78(i7)%c(1)+ccr*(cauxc+ceps_1)
      u512_78(i7)%c(2)=u512_78(i7)%c(2)+ccl*(-cauxc+ceps_1)
      u512_78(i7)%d(1)=u512_78(i7)%d(1)+ccl*cf78(i7)%ek0
      u512_78(i7)%d(2)=u512_78(i7)%d(2)+ccr*cf78(i7)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TTR0 -- aa=r6_3478(i5,i7)%a,bb=r6_3478(i5,i7)%b,a1=u512_78(i7)%a,b1=u5
* 12_78(i7)%b,c1=u512_78(i7)%c,d1=u512_78(i7)%d,a2=r6_34(i5)%a,b2=r6_34(
* i5)%b,prq=s634,nsum=1
      r6_3478(i5,i7)%a(1)=r6_3478(i5,i7)%a(1)+u512_78(i7)%a(1)*r
     & 6_34(i5)%a(1)+u512_78(i7)%c(1)*s634*r6_34(i5)%b(2)
      r6_3478(i5,i7)%b(1)=r6_3478(i5,i7)%b(1)+u512_78(i7)%d(1)*s
     & 634*r6_34(i5)%b(1)+u512_78(i7)%b(1)*r6_34(i5)%a(2)
      r6_3478(i5,i7)%b(2)=r6_3478(i5,i7)%b(2)+u512_78(i7)%b(2)*r
     & 6_34(i5)%a(1)+u512_78(i7)%d(2)*s634*r6_34(i5)%b(2)
      r6_3478(i5,i7)%a(2)=r6_3478(i5,i7)%a(2)+u512_78(i7)%c(2)*s
     & 634*r6_34(i5)%b(1)+u512_78(i7)%a(2)*r6_34(i5)%a(2)
      end do
      end do
  
* II                                                                    
* quqd -- p=p534,q=p612
      quqd=p534(0)*p612(0)-p534(1)*p612(1)-p534(2)*p612(2)-p534(
     & 3)*p612(3)
      ccr=zcr(id5)/(f612)
      ccl=zcl(id5)/(f612)
      do i7=1,2
* T0 -- qu=p534,qd=p612,v=cz78(i7)%e,a=u534_78(i7)%a,b=u534_78(i7)%b,c=u
* 534_78(i7)%c,d=u534_78(i7)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz78(i7)%ek0*(p534(2)*p612(3)-p612(2)*p534(3))+p53
     & 4k0*(cz78(i7)%e(2)*p612(3)-p612(2)*cz78(i7)%e(3))-p612k0*
     & (cz78(i7)%e(2)*p534(3)-p534(2)*cz78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz78(i7)%e(3)*p534k0+p534(3)*cz78(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz78(i7)%e(3)*p612k0+p612(3)*cz78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz78(i7)%e(0)*p534(0)-cz78(i7)%e(1)*p534(1)-cz78(i7)%
     & e(2)*p534(2)-cz78(i7)%e(3)*p534(3)
      cvqd=cz78(i7)%e(0)*p612(0)-cz78(i7)%e(1)*p612(1)-cz78(i7)%
     & e(2)*p612(2)-cz78(i7)%e(3)*p612(3)
      cauxa=-cz78(i7)%ek0*quqd+p534k0*cvqd+p612k0*cvqu
      cauxb=-cz78(i7)%ek0*p612(2)+p612k0*cz78(i7)%e(2)
      cauxc=+cz78(i7)%ek0*p534(2)-p534k0*cz78(i7)%e(2)
      u534_78(i7)%a(1)=ccr*(cauxa+ceps_0)
      u534_78(i7)%a(2)=ccl*(cauxa-ceps_0)
      u534_78(i7)%b(1)=ccl*(cauxb-ceps_2)
      u534_78(i7)%b(2)=ccr*(-cauxb-ceps_2)
      u534_78(i7)%c(1)=ccr*(cauxc+ceps_1)
      u534_78(i7)%c(2)=ccl*(-cauxc+ceps_1)
      u534_78(i7)%d(1)=ccl*cz78(i7)%ek0
      u534_78(i7)%d(2)=ccr*cz78(i7)%ek0
      end do
  
      ccr=fcr(id5)/(f612)
      ccl=fcl(id5)/(f612)
      do i7=1,2
* T0 -- qu=p534,qd=p612,v=cf78(i7)%e,a=u534_78(i7)%a,b=u534_78(i7)%b,c=u
* 534_78(i7)%c,d=u534_78(i7)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf78(i7)%ek0*(p534(2)*p612(3)-p612(2)*p534(3))+p53
     & 4k0*(cf78(i7)%e(2)*p612(3)-p612(2)*cf78(i7)%e(3))-p612k0*
     & (cf78(i7)%e(2)*p534(3)-p534(2)*cf78(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf78(i7)%e(3)*p534k0+p534(3)*cf78(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf78(i7)%e(3)*p612k0+p612(3)*cf78(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf78(i7)%e(0)*p534(0)-cf78(i7)%e(1)*p534(1)-cf78(i7)%
     & e(2)*p534(2)-cf78(i7)%e(3)*p534(3)
      cvqd=cf78(i7)%e(0)*p612(0)-cf78(i7)%e(1)*p612(1)-cf78(i7)%
     & e(2)*p612(2)-cf78(i7)%e(3)*p612(3)
      cauxa=-cf78(i7)%ek0*quqd+p534k0*cvqd+p612k0*cvqu
      cauxb=-cf78(i7)%ek0*p612(2)+p612k0*cf78(i7)%e(2)
      cauxc=+cf78(i7)%ek0*p534(2)-p534k0*cf78(i7)%e(2)
      u534_78(i7)%a(1)=u534_78(i7)%a(1)+ccr*(cauxa+ceps_0)
      u534_78(i7)%a(2)=u534_78(i7)%a(2)+ccl*(cauxa-ceps_0)
      u534_78(i7)%b(1)=u534_78(i7)%b(1)+ccl*(cauxb-ceps_2)
      u534_78(i7)%b(2)=u534_78(i7)%b(2)+ccr*(-cauxb-ceps_2)
      u534_78(i7)%c(1)=u534_78(i7)%c(1)+ccr*(cauxc+ceps_1)
      u534_78(i7)%c(2)=u534_78(i7)%c(2)+ccl*(-cauxc+ceps_1)
      u534_78(i7)%d(1)=u534_78(i7)%d(1)+ccl*cf78(i7)%ek0
      u534_78(i7)%d(2)=u534_78(i7)%d(2)+ccr*cf78(i7)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TLT0 -- aa=l5_3478(i5,i7)%a,cc=l5_3478(i5,i7)%c,a1=l5_34(i5)%a,c1=l5_3
* 4(i5)%c,a2=u534_78(i7)%a,b2=u534_78(i7)%b,c2=u534_78(i7)%c,d2=u534_78(
* i7)%d,prq=s534,nsum=0
      l5_3478(i5,i7)%a(1)=l5_34(i5)%a(1)*u534_78(i7)%a(1)+l5_34(
     & i5)%c(1)*s534*u534_78(i7)%b(2)
      l5_3478(i5,i7)%c(1)=l5_34(i5)%a(1)*u534_78(i7)%c(1)+l5_34(
     & i5)%c(1)*s534*u534_78(i7)%d(2)
      l5_3478(i5,i7)%c(2)=l5_34(i5)%c(2)*s534*u534_78(i7)%d(1)+l
     & 5_34(i5)%a(2)*u534_78(i7)%c(2)
      l5_3478(i5,i7)%a(2)=l5_34(i5)%c(2)*s534*u534_78(i7)%b(1)+l
     & 5_34(i5)%a(2)*u534_78(i7)%a(2)
      end do
      end do
  
  
* quqd -- p=p578,q=p612
      quqd=p578(0)*p612(0)-p578(1)*p612(1)-p578(2)*p612(2)-p578(
     & 3)*p612(3)
      ccr=zcr(id5)/(f612)
      ccl=zcl(id5)/(f612)
      do i5=1,2
* T0 -- qu=p578,qd=p612,v=cz34(i5)%e,a=u578_34(i5)%a,b=u578_34(i5)%b,c=u
* 578_34(i5)%c,d=u578_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p578(2)*p612(3)-p612(2)*p578(3))+p57
     & 8k0*(cz34(i5)%e(2)*p612(3)-p612(2)*cz34(i5)%e(3))-p612k0*
     & (cz34(i5)%e(2)*p578(3)-p578(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p578k0+p578(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p612k0+p612(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p578(0)-cz34(i5)%e(1)*p578(1)-cz34(i5)%
     & e(2)*p578(2)-cz34(i5)%e(3)*p578(3)
      cvqd=cz34(i5)%e(0)*p612(0)-cz34(i5)%e(1)*p612(1)-cz34(i5)%
     & e(2)*p612(2)-cz34(i5)%e(3)*p612(3)
      cauxa=-cz34(i5)%ek0*quqd+p578k0*cvqd+p612k0*cvqu
      cauxb=-cz34(i5)%ek0*p612(2)+p612k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p578(2)-p578k0*cz34(i5)%e(2)
      u578_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u578_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u578_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u578_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u578_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u578_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u578_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u578_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id5)/(f612)
      ccl=fcl(id5)/(f612)
      do i5=1,2
* T0 -- qu=p578,qd=p612,v=cf34(i5)%e,a=u578_34(i5)%a,b=u578_34(i5)%b,c=u
* 578_34(i5)%c,d=u578_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p578(2)*p612(3)-p612(2)*p578(3))+p57
     & 8k0*(cf34(i5)%e(2)*p612(3)-p612(2)*cf34(i5)%e(3))-p612k0*
     & (cf34(i5)%e(2)*p578(3)-p578(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p578k0+p578(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p612k0+p612(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p578(0)-cf34(i5)%e(1)*p578(1)-cf34(i5)%
     & e(2)*p578(2)-cf34(i5)%e(3)*p578(3)
      cvqd=cf34(i5)%e(0)*p612(0)-cf34(i5)%e(1)*p612(1)-cf34(i5)%
     & e(2)*p612(2)-cf34(i5)%e(3)*p612(3)
      cauxa=-cf34(i5)%ek0*quqd+p578k0*cvqd+p612k0*cvqu
      cauxb=-cf34(i5)%ek0*p612(2)+p612k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p578(2)-p578k0*cf34(i5)%e(2)
      u578_34(i5)%a(1)=u578_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u578_34(i5)%a(2)=u578_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u578_34(i5)%b(1)=u578_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u578_34(i5)%b(2)=u578_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u578_34(i5)%c(1)=u578_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u578_34(i5)%c(2)=u578_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u578_34(i5)%d(1)=u578_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u578_34(i5)%d(2)=u578_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TLT0 -- aa=l5_3478(i5,i7)%a,cc=l5_3478(i5,i7)%c,a1=l5_78(i7)%a,c1=l5_7
* 8(i7)%c,a2=u578_34(i5)%a,b2=u578_34(i5)%b,c2=u578_34(i5)%c,d2=u578_34(
* i5)%d,prq=s578,nsum=1
      l5_3478(i5,i7)%a(1)=l5_3478(i5,i7)%a(1)+l5_78(i7)%a(1)*u57
     & 8_34(i5)%a(1)+l5_78(i7)%c(1)*s578*u578_34(i5)%b(2)
      l5_3478(i5,i7)%c(1)=l5_3478(i5,i7)%c(1)+l5_78(i7)%a(1)*u57
     & 8_34(i5)%c(1)+l5_78(i7)%c(1)*s578*u578_34(i5)%d(2)
      l5_3478(i5,i7)%c(2)=l5_3478(i5,i7)%c(2)+l5_78(i7)%c(2)*s57
     & 8*u578_34(i5)%d(1)+l5_78(i7)%a(2)*u578_34(i5)%c(2)
      l5_3478(i5,i7)%a(2)=l5_3478(i5,i7)%a(2)+l5_78(i7)%c(2)*s57
     & 8*u578_34(i5)%b(1)+l5_78(i7)%a(2)*u578_34(i5)%a(2)
      end do
      end do
  
* III                                                                   
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r6_278(i2,i7)%a,bb=r6_278(i2,i7)%b,a1=u5314_2(i2)%a,b1=u531
* 4_2(i2)%b,c1=u5314_2(i2)%c,d1=u5314_2(i2)%d,a2=r6_78(i7)%a,b2=r6_78(i7
* )%b,prq=s678,nsum=0
      r6_278(i2,i7)%a(1)=u5314_2(i2)%a(1)*r6_78(i7)%a(1)+u5314_2
     & (i2)%c(1)*s678*r6_78(i7)%b(2)
      r6_278(i2,i7)%b(1)=u5314_2(i2)%d(1)*s678*r6_78(i7)%b(1)+u5
     & 314_2(i2)%b(1)*r6_78(i7)%a(2)
      r6_278(i2,i7)%b(2)=u5314_2(i2)%b(2)*r6_78(i7)%a(1)+u5314_2
     & (i2)%d(2)*s678*r6_78(i7)%b(2)
      r6_278(i2,i7)%a(2)=u5314_2(i2)%c(2)*s678*r6_78(i7)%b(1)+u5
     & 314_2(i2)%a(2)*r6_78(i7)%a(2)
      end do
      end do
  
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r6_278(i2,i7)%a,bb=r6_278(i2,i7)%b,a1=u5314_78(i7)%a,b1=u53
* 14_78(i7)%b,c1=u5314_78(i7)%c,d1=u5314_78(i7)%d,a2=r6_2(i2)%a,b2=r6_2(
* i2)%b,prq=s62,nsum=1
      r6_278(i2,i7)%a(1)=r6_278(i2,i7)%a(1)+u5314_78(i7)%a(1)*r6
     & _2(i2)%a(1)+u5314_78(i7)%c(1)*s62*r6_2(i2)%b(2)
      r6_278(i2,i7)%b(1)=r6_278(i2,i7)%b(1)+u5314_78(i7)%d(1)*s6
     & 2*r6_2(i2)%b(1)+u5314_78(i7)%b(1)*r6_2(i2)%a(2)
      r6_278(i2,i7)%b(2)=r6_278(i2,i7)%b(2)+u5314_78(i7)%b(2)*r6
     & _2(i2)%a(1)+u5314_78(i7)%d(2)*s62*r6_2(i2)%b(2)
      r6_278(i2,i7)%a(2)=r6_278(i2,i7)%a(2)+u5314_78(i7)%c(2)*s6
     & 2*r6_2(i2)%b(1)+u5314_78(i7)%a(2)*r6_2(i2)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r6_234(i2,i5)%a,bb=r6_234(i2,i5)%b,a1=u5718_2(i2)%a,b1=u571
* 8_2(i2)%b,c1=u5718_2(i2)%c,d1=u5718_2(i2)%d,a2=r6_34(i5)%a,b2=r6_34(i5
* )%b,prq=s634,nsum=0
      r6_234(i2,i5)%a(1)=u5718_2(i2)%a(1)*r6_34(i5)%a(1)+u5718_2
     & (i2)%c(1)*s634*r6_34(i5)%b(2)
      r6_234(i2,i5)%b(1)=u5718_2(i2)%d(1)*s634*r6_34(i5)%b(1)+u5
     & 718_2(i2)%b(1)*r6_34(i5)%a(2)
      r6_234(i2,i5)%b(2)=u5718_2(i2)%b(2)*r6_34(i5)%a(1)+u5718_2
     & (i2)%d(2)*s634*r6_34(i5)%b(2)
      r6_234(i2,i5)%a(2)=u5718_2(i2)%c(2)*s634*r6_34(i5)%b(1)+u5
     & 718_2(i2)%a(2)*r6_34(i5)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r6_234(i2,i5)%a,bb=r6_234(i2,i5)%b,a1=u5718_34(i5)%a,b1=u57
* 18_34(i5)%b,c1=u5718_34(i5)%c,d1=u5718_34(i5)%d,a2=r6_2(i2)%a,b2=r6_2(
* i2)%b,prq=s62,nsum=1
      r6_234(i2,i5)%a(1)=r6_234(i2,i5)%a(1)+u5718_34(i5)%a(1)*r6
     & _2(i2)%a(1)+u5718_34(i5)%c(1)*s62*r6_2(i2)%b(2)
      r6_234(i2,i5)%b(1)=r6_234(i2,i5)%b(1)+u5718_34(i5)%d(1)*s6
     & 2*r6_2(i2)%b(1)+u5718_34(i5)%b(1)*r6_2(i2)%a(2)
      r6_234(i2,i5)%b(2)=r6_234(i2,i5)%b(2)+u5718_34(i5)%b(2)*r6
     & _2(i2)%a(1)+u5718_34(i5)%d(2)*s62*r6_2(i2)%b(2)
      r6_234(i2,i5)%a(2)=r6_234(i2,i5)%a(2)+u5718_34(i5)%c(2)*s6
     & 2*r6_2(i2)%b(1)+u5718_34(i5)%a(2)*r6_2(i2)%a(2)
      end do
      end do
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r6_178(i2,i7)%a,bb=r6_178(i2,i7)%b,a1=u5324_1(i2)%a,b1=u532
* 4_1(i2)%b,c1=u5324_1(i2)%c,d1=u5324_1(i2)%d,a2=r6_78(i7)%a,b2=r6_78(i7
* )%b,prq=s678,nsum=0
      r6_178(i2,i7)%a(1)=u5324_1(i2)%a(1)*r6_78(i7)%a(1)+u5324_1
     & (i2)%c(1)*s678*r6_78(i7)%b(2)
      r6_178(i2,i7)%b(1)=u5324_1(i2)%d(1)*s678*r6_78(i7)%b(1)+u5
     & 324_1(i2)%b(1)*r6_78(i7)%a(2)
      r6_178(i2,i7)%b(2)=u5324_1(i2)%b(2)*r6_78(i7)%a(1)+u5324_1
     & (i2)%d(2)*s678*r6_78(i7)%b(2)
      r6_178(i2,i7)%a(2)=u5324_1(i2)%c(2)*s678*r6_78(i7)%b(1)+u5
     & 324_1(i2)%a(2)*r6_78(i7)%a(2)
      end do
      end do
  
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r6_178(i2,i7)%a,bb=r6_178(i2,i7)%b,a1=u5324_78(i7)%a,b1=u53
* 24_78(i7)%b,c1=u5324_78(i7)%c,d1=u5324_78(i7)%d,a2=r6_1(i2)%a,b2=r6_1(
* i2)%b,prq=s61,nsum=1
      r6_178(i2,i7)%a(1)=r6_178(i2,i7)%a(1)+u5324_78(i7)%a(1)*r6
     & _1(i2)%a(1)+u5324_78(i7)%c(1)*s61*r6_1(i2)%b(2)
      r6_178(i2,i7)%b(1)=r6_178(i2,i7)%b(1)+u5324_78(i7)%d(1)*s6
     & 1*r6_1(i2)%b(1)+u5324_78(i7)%b(1)*r6_1(i2)%a(2)
      r6_178(i2,i7)%b(2)=r6_178(i2,i7)%b(2)+u5324_78(i7)%b(2)*r6
     & _1(i2)%a(1)+u5324_78(i7)%d(2)*s61*r6_1(i2)%b(2)
      r6_178(i2,i7)%a(2)=r6_178(i2,i7)%a(2)+u5324_78(i7)%c(2)*s6
     & 1*r6_1(i2)%b(1)+u5324_78(i7)%a(2)*r6_1(i2)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r6_134(i2,i5)%a,bb=r6_134(i2,i5)%b,a1=u5728_1(i2)%a,b1=u572
* 8_1(i2)%b,c1=u5728_1(i2)%c,d1=u5728_1(i2)%d,a2=r6_34(i5)%a,b2=r6_34(i5
* )%b,prq=s634,nsum=0
      r6_134(i2,i5)%a(1)=u5728_1(i2)%a(1)*r6_34(i5)%a(1)+u5728_1
     & (i2)%c(1)*s634*r6_34(i5)%b(2)
      r6_134(i2,i5)%b(1)=u5728_1(i2)%d(1)*s634*r6_34(i5)%b(1)+u5
     & 728_1(i2)%b(1)*r6_34(i5)%a(2)
      r6_134(i2,i5)%b(2)=u5728_1(i2)%b(2)*r6_34(i5)%a(1)+u5728_1
     & (i2)%d(2)*s634*r6_34(i5)%b(2)
      r6_134(i2,i5)%a(2)=u5728_1(i2)%c(2)*s634*r6_34(i5)%b(1)+u5
     & 728_1(i2)%a(2)*r6_34(i5)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r6_134(i2,i5)%a,bb=r6_134(i2,i5)%b,a1=u5728_34(i5)%a,b1=u57
* 28_34(i5)%b,c1=u5728_34(i5)%c,d1=u5728_34(i5)%d,a2=r6_1(i2)%a,b2=r6_1(
* i2)%b,prq=s61,nsum=1
      r6_134(i2,i5)%a(1)=r6_134(i2,i5)%a(1)+u5728_34(i5)%a(1)*r6
     & _1(i2)%a(1)+u5728_34(i5)%c(1)*s61*r6_1(i2)%b(2)
      r6_134(i2,i5)%b(1)=r6_134(i2,i5)%b(1)+u5728_34(i5)%d(1)*s6
     & 1*r6_1(i2)%b(1)+u5728_34(i5)%b(1)*r6_1(i2)%a(2)
      r6_134(i2,i5)%b(2)=r6_134(i2,i5)%b(2)+u5728_34(i5)%b(2)*r6
     & _1(i2)%a(1)+u5728_34(i5)%d(2)*s61*r6_1(i2)%b(2)
      r6_134(i2,i5)%a(2)=r6_134(i2,i5)%a(2)+u5728_34(i5)%c(2)*s6
     & 1*r6_1(i2)%b(1)+u5728_34(i5)%a(2)*r6_1(i2)%a(2)
      end do
      end do
  
* IV                                                                    
* quqd -- p=p534,q=p678
      quqd=p534(0)*p678(0)-p534(1)*p678(1)-p534(2)*p678(2)-p534(
     & 3)*p678(3)
      ccr=1d0/(f678*s12)
      ccl=1d0/(f678*s12)
      do i1=1,2
      do i2=1,2
* T0 -- qu=p534,qd=p678,v=ctrip12(i1,i2)%e,a=u534_gg(i1,i2)%a,b=u534_gg(
* i1,i2)%b,c=u534_gg(i1,i2)%c,d=u534_gg(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p534(2)*p678(3)-p678(2)*p534(3
     & ))+p534k0*(ctrip12(i1,i2)%e(2)*p678(3)-p678(2)*ctrip12(i1
     & ,i2)%e(3))-p678k0*(ctrip12(i1,i2)%e(2)*p534(3)-p534(2)*ct
     & rip12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ctrip12(i1,i2)%e(3)*p534k0+p534(3)*ctrip12(i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      ceps_2=-ctrip12(i1,i2)%e(3)*p678k0+p678(3)*ctrip12(i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=ctrip12(i1,i2)%e(0)*p534(0)-ctrip12(i1,i2)%e(1)*p534(
     & 1)-ctrip12(i1,i2)%e(2)*p534(2)-ctrip12(i1,i2)%e(3)*p534(3
     & )
      cvqd=ctrip12(i1,i2)%e(0)*p678(0)-ctrip12(i1,i2)%e(1)*p678(
     & 1)-ctrip12(i1,i2)%e(2)*p678(2)-ctrip12(i1,i2)%e(3)*p678(3
     & )
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p534k0*cvqd+p678k0*cvqu
      cauxb=-ctrip12(i1,i2)%ek0*p678(2)+p678k0*ctrip12(i1,i2)%e(
     & 2)
      cauxc=+ctrip12(i1,i2)%ek0*p534(2)-p534k0*ctrip12(i1,i2)%e(
     & 2)
      u534_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      u534_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u534_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u534_gg(i1,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u534_gg(i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      u534_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u534_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      u534_gg(i1,i2)%d(2)=ccr*ctrip12(i1,i2)%ek0
      end do
      end do
  
* quqd -- p=p578,q=p634
      quqd=p578(0)*p634(0)-p578(1)*p634(1)-p578(2)*p634(2)-p578(
     & 3)*p634(3)
      ccr=1d0/(f634*s12)
      ccl=1d0/(f634*s12)
      do i1=1,2
      do i2=1,2
* T0 -- qu=p578,qd=p634,v=ctrip12(i1,i2)%e,a=u578_gg(i1,i2)%a,b=u578_gg(
* i1,i2)%b,c=u578_gg(i1,i2)%c,d=u578_gg(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p578(2)*p634(3)-p634(2)*p578(3
     & ))+p578k0*(ctrip12(i1,i2)%e(2)*p634(3)-p634(2)*ctrip12(i1
     & ,i2)%e(3))-p634k0*(ctrip12(i1,i2)%e(2)*p578(3)-p578(2)*ct
     & rip12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ctrip12(i1,i2)%e(3)*p578k0+p578(3)*ctrip12(i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      ceps_2=-ctrip12(i1,i2)%e(3)*p634k0+p634(3)*ctrip12(i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=ctrip12(i1,i2)%e(0)*p578(0)-ctrip12(i1,i2)%e(1)*p578(
     & 1)-ctrip12(i1,i2)%e(2)*p578(2)-ctrip12(i1,i2)%e(3)*p578(3
     & )
      cvqd=ctrip12(i1,i2)%e(0)*p634(0)-ctrip12(i1,i2)%e(1)*p634(
     & 1)-ctrip12(i1,i2)%e(2)*p634(2)-ctrip12(i1,i2)%e(3)*p634(3
     & )
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p578k0*cvqd+p634k0*cvqu
      cauxb=-ctrip12(i1,i2)%ek0*p634(2)+p634k0*ctrip12(i1,i2)%e(
     & 2)
      cauxc=+ctrip12(i1,i2)%ek0*p578(2)-p578k0*ctrip12(i1,i2)%e(
     & 2)
      u578_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      u578_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u578_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u578_gg(i1,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u578_gg(i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      u578_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u578_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      u578_gg(i1,i2)%d(2)=ccr*ctrip12(i1,i2)%ek0
      end do
      end do
  
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p712,q=p856
      quqd=p712(0)*p856(0)-p712(1)*p856(1)-p712(2)*p856(2)-p712(
     & 3)*p856(3)
      ccr=zcr(id7)/(f856)
      ccl=zcl(id7)/(f856)
      do i5=1,2
* T0 -- qu=p712,qd=p856,v=cz34(i5)%e,a=u712_34(i5)%a,b=u712_34(i5)%b,c=u
* 712_34(i5)%c,d=u712_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p712(2)*p856(3)-p856(2)*p712(3))+p71
     & 2k0*(cz34(i5)%e(2)*p856(3)-p856(2)*cz34(i5)%e(3))-p856k0*
     & (cz34(i5)%e(2)*p712(3)-p712(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p712k0+p712(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p856k0+p856(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p712(0)-cz34(i5)%e(1)*p712(1)-cz34(i5)%
     & e(2)*p712(2)-cz34(i5)%e(3)*p712(3)
      cvqd=cz34(i5)%e(0)*p856(0)-cz34(i5)%e(1)*p856(1)-cz34(i5)%
     & e(2)*p856(2)-cz34(i5)%e(3)*p856(3)
      cauxa=-cz34(i5)%ek0*quqd+p712k0*cvqd+p856k0*cvqu
      cauxb=-cz34(i5)%ek0*p856(2)+p856k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p712(2)-p712k0*cz34(i5)%e(2)
      u712_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u712_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u712_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u712_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u712_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u712_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u712_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u712_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f856)
      ccl=fcl(id7)/(f856)
      do i5=1,2
* T0 -- qu=p712,qd=p856,v=cf34(i5)%e,a=u712_34(i5)%a,b=u712_34(i5)%b,c=u
* 712_34(i5)%c,d=u712_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p712(2)*p856(3)-p856(2)*p712(3))+p71
     & 2k0*(cf34(i5)%e(2)*p856(3)-p856(2)*cf34(i5)%e(3))-p856k0*
     & (cf34(i5)%e(2)*p712(3)-p712(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p712k0+p712(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p856k0+p856(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p712(0)-cf34(i5)%e(1)*p712(1)-cf34(i5)%
     & e(2)*p712(2)-cf34(i5)%e(3)*p712(3)
      cvqd=cf34(i5)%e(0)*p856(0)-cf34(i5)%e(1)*p856(1)-cf34(i5)%
     & e(2)*p856(2)-cf34(i5)%e(3)*p856(3)
      cauxa=-cf34(i5)%ek0*quqd+p712k0*cvqd+p856k0*cvqu
      cauxb=-cf34(i5)%ek0*p856(2)+p856k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p712(2)-p712k0*cf34(i5)%e(2)
      u712_34(i5)%a(1)=u712_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u712_34(i5)%a(2)=u712_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u712_34(i5)%b(1)=u712_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u712_34(i5)%b(2)=u712_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u712_34(i5)%c(1)=u712_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u712_34(i5)%c(2)=u712_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u712_34(i5)%d(1)=u712_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u712_34(i5)%d(2)=u712_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TTR0 -- aa=r8_3456(i5,i7)%a,bb=r8_3456(i5,i7)%b,a1=u712_34(i5)%a,b1=u7
* 12_34(i5)%b,c1=u712_34(i5)%c,d1=u712_34(i5)%d,a2=r8_56(i7)%a,b2=r8_56(
* i7)%b,prq=s856,nsum=0
      r8_3456(i5,i7)%a(1)=u712_34(i5)%a(1)*r8_56(i7)%a(1)+u712_3
     & 4(i5)%c(1)*s856*r8_56(i7)%b(2)
      r8_3456(i5,i7)%b(1)=u712_34(i5)%d(1)*s856*r8_56(i7)%b(1)+u
     & 712_34(i5)%b(1)*r8_56(i7)%a(2)
      r8_3456(i5,i7)%b(2)=u712_34(i5)%b(2)*r8_56(i7)%a(1)+u712_3
     & 4(i5)%d(2)*s856*r8_56(i7)%b(2)
      r8_3456(i5,i7)%a(2)=u712_34(i5)%c(2)*s856*r8_56(i7)%b(1)+u
     & 712_34(i5)%a(2)*r8_56(i7)%a(2)
      end do
      end do
  
  
* quqd -- p=p712,q=p834
      quqd=p712(0)*p834(0)-p712(1)*p834(1)-p712(2)*p834(2)-p712(
     & 3)*p834(3)
      ccr=zcr(id7)/(f834)
      ccl=zcl(id7)/(f834)
      do i7=1,2
* T0 -- qu=p712,qd=p834,v=cz56(i7)%e,a=u712_56(i7)%a,b=u712_56(i7)%b,c=u
* 712_56(i7)%c,d=u712_56(i7)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i7)%ek0*(p712(2)*p834(3)-p834(2)*p712(3))+p71
     & 2k0*(cz56(i7)%e(2)*p834(3)-p834(2)*cz56(i7)%e(3))-p834k0*
     & (cz56(i7)%e(2)*p712(3)-p712(2)*cz56(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i7)%e(3)*p712k0+p712(3)*cz56(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i7)%e(3)*p834k0+p834(3)*cz56(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i7)%e(0)*p712(0)-cz56(i7)%e(1)*p712(1)-cz56(i7)%
     & e(2)*p712(2)-cz56(i7)%e(3)*p712(3)
      cvqd=cz56(i7)%e(0)*p834(0)-cz56(i7)%e(1)*p834(1)-cz56(i7)%
     & e(2)*p834(2)-cz56(i7)%e(3)*p834(3)
      cauxa=-cz56(i7)%ek0*quqd+p712k0*cvqd+p834k0*cvqu
      cauxb=-cz56(i7)%ek0*p834(2)+p834k0*cz56(i7)%e(2)
      cauxc=+cz56(i7)%ek0*p712(2)-p712k0*cz56(i7)%e(2)
      u712_56(i7)%a(1)=ccr*(cauxa+ceps_0)
      u712_56(i7)%a(2)=ccl*(cauxa-ceps_0)
      u712_56(i7)%b(1)=ccl*(cauxb-ceps_2)
      u712_56(i7)%b(2)=ccr*(-cauxb-ceps_2)
      u712_56(i7)%c(1)=ccr*(cauxc+ceps_1)
      u712_56(i7)%c(2)=ccl*(-cauxc+ceps_1)
      u712_56(i7)%d(1)=ccl*cz56(i7)%ek0
      u712_56(i7)%d(2)=ccr*cz56(i7)%ek0
      end do
  
      ccr=fcr(id7)/(f834)
      ccl=fcl(id7)/(f834)
      do i7=1,2
* T0 -- qu=p712,qd=p834,v=cf56(i7)%e,a=u712_56(i7)%a,b=u712_56(i7)%b,c=u
* 712_56(i7)%c,d=u712_56(i7)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i7)%ek0*(p712(2)*p834(3)-p834(2)*p712(3))+p71
     & 2k0*(cf56(i7)%e(2)*p834(3)-p834(2)*cf56(i7)%e(3))-p834k0*
     & (cf56(i7)%e(2)*p712(3)-p712(2)*cf56(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i7)%e(3)*p712k0+p712(3)*cf56(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i7)%e(3)*p834k0+p834(3)*cf56(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i7)%e(0)*p712(0)-cf56(i7)%e(1)*p712(1)-cf56(i7)%
     & e(2)*p712(2)-cf56(i7)%e(3)*p712(3)
      cvqd=cf56(i7)%e(0)*p834(0)-cf56(i7)%e(1)*p834(1)-cf56(i7)%
     & e(2)*p834(2)-cf56(i7)%e(3)*p834(3)
      cauxa=-cf56(i7)%ek0*quqd+p712k0*cvqd+p834k0*cvqu
      cauxb=-cf56(i7)%ek0*p834(2)+p834k0*cf56(i7)%e(2)
      cauxc=+cf56(i7)%ek0*p712(2)-p712k0*cf56(i7)%e(2)
      u712_56(i7)%a(1)=u712_56(i7)%a(1)+ccr*(cauxa+ceps_0)
      u712_56(i7)%a(2)=u712_56(i7)%a(2)+ccl*(cauxa-ceps_0)
      u712_56(i7)%b(1)=u712_56(i7)%b(1)+ccl*(cauxb-ceps_2)
      u712_56(i7)%b(2)=u712_56(i7)%b(2)+ccr*(-cauxb-ceps_2)
      u712_56(i7)%c(1)=u712_56(i7)%c(1)+ccr*(cauxc+ceps_1)
      u712_56(i7)%c(2)=u712_56(i7)%c(2)+ccl*(-cauxc+ceps_1)
      u712_56(i7)%d(1)=u712_56(i7)%d(1)+ccl*cf56(i7)%ek0
      u712_56(i7)%d(2)=u712_56(i7)%d(2)+ccr*cf56(i7)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TTR0 -- aa=r8_3456(i5,i7)%a,bb=r8_3456(i5,i7)%b,a1=u712_56(i7)%a,b1=u7
* 12_56(i7)%b,c1=u712_56(i7)%c,d1=u712_56(i7)%d,a2=r8_34(i5)%a,b2=r8_34(
* i5)%b,prq=s834,nsum=1
      r8_3456(i5,i7)%a(1)=r8_3456(i5,i7)%a(1)+u712_56(i7)%a(1)*r
     & 8_34(i5)%a(1)+u712_56(i7)%c(1)*s834*r8_34(i5)%b(2)
      r8_3456(i5,i7)%b(1)=r8_3456(i5,i7)%b(1)+u712_56(i7)%d(1)*s
     & 834*r8_34(i5)%b(1)+u712_56(i7)%b(1)*r8_34(i5)%a(2)
      r8_3456(i5,i7)%b(2)=r8_3456(i5,i7)%b(2)+u712_56(i7)%b(2)*r
     & 8_34(i5)%a(1)+u712_56(i7)%d(2)*s834*r8_34(i5)%b(2)
      r8_3456(i5,i7)%a(2)=r8_3456(i5,i7)%a(2)+u712_56(i7)%c(2)*s
     & 834*r8_34(i5)%b(1)+u712_56(i7)%a(2)*r8_34(i5)%a(2)
      end do
      end do
  
* II                                                                    
* quqd -- p=p734,q=p812
      quqd=p734(0)*p812(0)-p734(1)*p812(1)-p734(2)*p812(2)-p734(
     & 3)*p812(3)
      ccr=zcr(id7)/(f812)
      ccl=zcl(id7)/(f812)
      do i7=1,2
* T0 -- qu=p734,qd=p812,v=cz56(i7)%e,a=u734_56(i7)%a,b=u734_56(i7)%b,c=u
* 734_56(i7)%c,d=u734_56(i7)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz56(i7)%ek0*(p734(2)*p812(3)-p812(2)*p734(3))+p73
     & 4k0*(cz56(i7)%e(2)*p812(3)-p812(2)*cz56(i7)%e(3))-p812k0*
     & (cz56(i7)%e(2)*p734(3)-p734(2)*cz56(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz56(i7)%e(3)*p734k0+p734(3)*cz56(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz56(i7)%e(3)*p812k0+p812(3)*cz56(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz56(i7)%e(0)*p734(0)-cz56(i7)%e(1)*p734(1)-cz56(i7)%
     & e(2)*p734(2)-cz56(i7)%e(3)*p734(3)
      cvqd=cz56(i7)%e(0)*p812(0)-cz56(i7)%e(1)*p812(1)-cz56(i7)%
     & e(2)*p812(2)-cz56(i7)%e(3)*p812(3)
      cauxa=-cz56(i7)%ek0*quqd+p734k0*cvqd+p812k0*cvqu
      cauxb=-cz56(i7)%ek0*p812(2)+p812k0*cz56(i7)%e(2)
      cauxc=+cz56(i7)%ek0*p734(2)-p734k0*cz56(i7)%e(2)
      u734_56(i7)%a(1)=ccr*(cauxa+ceps_0)
      u734_56(i7)%a(2)=ccl*(cauxa-ceps_0)
      u734_56(i7)%b(1)=ccl*(cauxb-ceps_2)
      u734_56(i7)%b(2)=ccr*(-cauxb-ceps_2)
      u734_56(i7)%c(1)=ccr*(cauxc+ceps_1)
      u734_56(i7)%c(2)=ccl*(-cauxc+ceps_1)
      u734_56(i7)%d(1)=ccl*cz56(i7)%ek0
      u734_56(i7)%d(2)=ccr*cz56(i7)%ek0
      end do
  
      ccr=fcr(id7)/(f812)
      ccl=fcl(id7)/(f812)
      do i7=1,2
* T0 -- qu=p734,qd=p812,v=cf56(i7)%e,a=u734_56(i7)%a,b=u734_56(i7)%b,c=u
* 734_56(i7)%c,d=u734_56(i7)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf56(i7)%ek0*(p734(2)*p812(3)-p812(2)*p734(3))+p73
     & 4k0*(cf56(i7)%e(2)*p812(3)-p812(2)*cf56(i7)%e(3))-p812k0*
     & (cf56(i7)%e(2)*p734(3)-p734(2)*cf56(i7)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf56(i7)%e(3)*p734k0+p734(3)*cf56(i7)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf56(i7)%e(3)*p812k0+p812(3)*cf56(i7)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf56(i7)%e(0)*p734(0)-cf56(i7)%e(1)*p734(1)-cf56(i7)%
     & e(2)*p734(2)-cf56(i7)%e(3)*p734(3)
      cvqd=cf56(i7)%e(0)*p812(0)-cf56(i7)%e(1)*p812(1)-cf56(i7)%
     & e(2)*p812(2)-cf56(i7)%e(3)*p812(3)
      cauxa=-cf56(i7)%ek0*quqd+p734k0*cvqd+p812k0*cvqu
      cauxb=-cf56(i7)%ek0*p812(2)+p812k0*cf56(i7)%e(2)
      cauxc=+cf56(i7)%ek0*p734(2)-p734k0*cf56(i7)%e(2)
      u734_56(i7)%a(1)=u734_56(i7)%a(1)+ccr*(cauxa+ceps_0)
      u734_56(i7)%a(2)=u734_56(i7)%a(2)+ccl*(cauxa-ceps_0)
      u734_56(i7)%b(1)=u734_56(i7)%b(1)+ccl*(cauxb-ceps_2)
      u734_56(i7)%b(2)=u734_56(i7)%b(2)+ccr*(-cauxb-ceps_2)
      u734_56(i7)%c(1)=u734_56(i7)%c(1)+ccr*(cauxc+ceps_1)
      u734_56(i7)%c(2)=u734_56(i7)%c(2)+ccl*(-cauxc+ceps_1)
      u734_56(i7)%d(1)=u734_56(i7)%d(1)+ccl*cf56(i7)%ek0
      u734_56(i7)%d(2)=u734_56(i7)%d(2)+ccr*cf56(i7)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TLT0 -- aa=l7_3456(i5,i7)%a,cc=l7_3456(i5,i7)%c,a1=l7_34(i5)%a,c1=l7_3
* 4(i5)%c,a2=u734_56(i7)%a,b2=u734_56(i7)%b,c2=u734_56(i7)%c,d2=u734_56(
* i7)%d,prq=s734,nsum=0
      l7_3456(i5,i7)%a(1)=l7_34(i5)%a(1)*u734_56(i7)%a(1)+l7_34(
     & i5)%c(1)*s734*u734_56(i7)%b(2)
      l7_3456(i5,i7)%c(1)=l7_34(i5)%a(1)*u734_56(i7)%c(1)+l7_34(
     & i5)%c(1)*s734*u734_56(i7)%d(2)
      l7_3456(i5,i7)%c(2)=l7_34(i5)%c(2)*s734*u734_56(i7)%d(1)+l
     & 7_34(i5)%a(2)*u734_56(i7)%c(2)
      l7_3456(i5,i7)%a(2)=l7_34(i5)%c(2)*s734*u734_56(i7)%b(1)+l
     & 7_34(i5)%a(2)*u734_56(i7)%a(2)
      end do
      end do
  
  
* quqd -- p=p756,q=p812
      quqd=p756(0)*p812(0)-p756(1)*p812(1)-p756(2)*p812(2)-p756(
     & 3)*p812(3)
      ccr=zcr(id7)/(f812)
      ccl=zcl(id7)/(f812)
      do i5=1,2
* T0 -- qu=p756,qd=p812,v=cz34(i5)%e,a=u756_34(i5)%a,b=u756_34(i5)%b,c=u
* 756_34(i5)%c,d=u756_34(i5)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i5)%ek0*(p756(2)*p812(3)-p812(2)*p756(3))+p75
     & 6k0*(cz34(i5)%e(2)*p812(3)-p812(2)*cz34(i5)%e(3))-p812k0*
     & (cz34(i5)%e(2)*p756(3)-p756(2)*cz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i5)%e(3)*p756k0+p756(3)*cz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i5)%e(3)*p812k0+p812(3)*cz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i5)%e(0)*p756(0)-cz34(i5)%e(1)*p756(1)-cz34(i5)%
     & e(2)*p756(2)-cz34(i5)%e(3)*p756(3)
      cvqd=cz34(i5)%e(0)*p812(0)-cz34(i5)%e(1)*p812(1)-cz34(i5)%
     & e(2)*p812(2)-cz34(i5)%e(3)*p812(3)
      cauxa=-cz34(i5)%ek0*quqd+p756k0*cvqd+p812k0*cvqu
      cauxb=-cz34(i5)%ek0*p812(2)+p812k0*cz34(i5)%e(2)
      cauxc=+cz34(i5)%ek0*p756(2)-p756k0*cz34(i5)%e(2)
      u756_34(i5)%a(1)=ccr*(cauxa+ceps_0)
      u756_34(i5)%a(2)=ccl*(cauxa-ceps_0)
      u756_34(i5)%b(1)=ccl*(cauxb-ceps_2)
      u756_34(i5)%b(2)=ccr*(-cauxb-ceps_2)
      u756_34(i5)%c(1)=ccr*(cauxc+ceps_1)
      u756_34(i5)%c(2)=ccl*(-cauxc+ceps_1)
      u756_34(i5)%d(1)=ccl*cz34(i5)%ek0
      u756_34(i5)%d(2)=ccr*cz34(i5)%ek0
      end do
  
      ccr=fcr(id7)/(f812)
      ccl=fcl(id7)/(f812)
      do i5=1,2
* T0 -- qu=p756,qd=p812,v=cf34(i5)%e,a=u756_34(i5)%a,b=u756_34(i5)%b,c=u
* 756_34(i5)%c,d=u756_34(i5)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i5)%ek0*(p756(2)*p812(3)-p812(2)*p756(3))+p75
     & 6k0*(cf34(i5)%e(2)*p812(3)-p812(2)*cf34(i5)%e(3))-p812k0*
     & (cf34(i5)%e(2)*p756(3)-p756(2)*cf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i5)%e(3)*p756k0+p756(3)*cf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i5)%e(3)*p812k0+p812(3)*cf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i5)%e(0)*p756(0)-cf34(i5)%e(1)*p756(1)-cf34(i5)%
     & e(2)*p756(2)-cf34(i5)%e(3)*p756(3)
      cvqd=cf34(i5)%e(0)*p812(0)-cf34(i5)%e(1)*p812(1)-cf34(i5)%
     & e(2)*p812(2)-cf34(i5)%e(3)*p812(3)
      cauxa=-cf34(i5)%ek0*quqd+p756k0*cvqd+p812k0*cvqu
      cauxb=-cf34(i5)%ek0*p812(2)+p812k0*cf34(i5)%e(2)
      cauxc=+cf34(i5)%ek0*p756(2)-p756k0*cf34(i5)%e(2)
      u756_34(i5)%a(1)=u756_34(i5)%a(1)+ccr*(cauxa+ceps_0)
      u756_34(i5)%a(2)=u756_34(i5)%a(2)+ccl*(cauxa-ceps_0)
      u756_34(i5)%b(1)=u756_34(i5)%b(1)+ccl*(cauxb-ceps_2)
      u756_34(i5)%b(2)=u756_34(i5)%b(2)+ccr*(-cauxb-ceps_2)
      u756_34(i5)%c(1)=u756_34(i5)%c(1)+ccr*(cauxc+ceps_1)
      u756_34(i5)%c(2)=u756_34(i5)%c(2)+ccl*(-cauxc+ceps_1)
      u756_34(i5)%d(1)=u756_34(i5)%d(1)+ccl*cf34(i5)%ek0
      u756_34(i5)%d(2)=u756_34(i5)%d(2)+ccr*cf34(i5)%ek0
      end do
  
      do i5=1,2
      do i7=1,2
* TLT0 -- aa=l7_3456(i5,i7)%a,cc=l7_3456(i5,i7)%c,a1=l7_56(i7)%a,c1=l7_5
* 6(i7)%c,a2=u756_34(i5)%a,b2=u756_34(i5)%b,c2=u756_34(i5)%c,d2=u756_34(
* i5)%d,prq=s756,nsum=1
      l7_3456(i5,i7)%a(1)=l7_3456(i5,i7)%a(1)+l7_56(i7)%a(1)*u75
     & 6_34(i5)%a(1)+l7_56(i7)%c(1)*s756*u756_34(i5)%b(2)
      l7_3456(i5,i7)%c(1)=l7_3456(i5,i7)%c(1)+l7_56(i7)%a(1)*u75
     & 6_34(i5)%c(1)+l7_56(i7)%c(1)*s756*u756_34(i5)%d(2)
      l7_3456(i5,i7)%c(2)=l7_3456(i5,i7)%c(2)+l7_56(i7)%c(2)*s75
     & 6*u756_34(i5)%d(1)+l7_56(i7)%a(2)*u756_34(i5)%c(2)
      l7_3456(i5,i7)%a(2)=l7_3456(i5,i7)%a(2)+l7_56(i7)%c(2)*s75
     & 6*u756_34(i5)%b(1)+l7_56(i7)%a(2)*u756_34(i5)%a(2)
      end do
      end do
  
* III                                                                   
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r8_256(i2,i7)%a,bb=r8_256(i2,i7)%b,a1=u7314_2(i2)%a,b1=u731
* 4_2(i2)%b,c1=u7314_2(i2)%c,d1=u7314_2(i2)%d,a2=r8_56(i7)%a,b2=r8_56(i7
* )%b,prq=s856,nsum=0
      r8_256(i2,i7)%a(1)=u7314_2(i2)%a(1)*r8_56(i7)%a(1)+u7314_2
     & (i2)%c(1)*s856*r8_56(i7)%b(2)
      r8_256(i2,i7)%b(1)=u7314_2(i2)%d(1)*s856*r8_56(i7)%b(1)+u7
     & 314_2(i2)%b(1)*r8_56(i7)%a(2)
      r8_256(i2,i7)%b(2)=u7314_2(i2)%b(2)*r8_56(i7)%a(1)+u7314_2
     & (i2)%d(2)*s856*r8_56(i7)%b(2)
      r8_256(i2,i7)%a(2)=u7314_2(i2)%c(2)*s856*r8_56(i7)%b(1)+u7
     & 314_2(i2)%a(2)*r8_56(i7)%a(2)
      end do
      end do
  
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r8_256(i2,i7)%a,bb=r8_256(i2,i7)%b,a1=u7314_56(i7)%a,b1=u73
* 14_56(i7)%b,c1=u7314_56(i7)%c,d1=u7314_56(i7)%d,a2=r8_2(i2)%a,b2=r8_2(
* i2)%b,prq=s82,nsum=1
      r8_256(i2,i7)%a(1)=r8_256(i2,i7)%a(1)+u7314_56(i7)%a(1)*r8
     & _2(i2)%a(1)+u7314_56(i7)%c(1)*s82*r8_2(i2)%b(2)
      r8_256(i2,i7)%b(1)=r8_256(i2,i7)%b(1)+u7314_56(i7)%d(1)*s8
     & 2*r8_2(i2)%b(1)+u7314_56(i7)%b(1)*r8_2(i2)%a(2)
      r8_256(i2,i7)%b(2)=r8_256(i2,i7)%b(2)+u7314_56(i7)%b(2)*r8
     & _2(i2)%a(1)+u7314_56(i7)%d(2)*s82*r8_2(i2)%b(2)
      r8_256(i2,i7)%a(2)=r8_256(i2,i7)%a(2)+u7314_56(i7)%c(2)*s8
     & 2*r8_2(i2)%b(1)+u7314_56(i7)%a(2)*r8_2(i2)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r8_234(i2,i5)%a,bb=r8_234(i2,i5)%b,a1=u7516_2(i2)%a,b1=u751
* 6_2(i2)%b,c1=u7516_2(i2)%c,d1=u7516_2(i2)%d,a2=r8_34(i5)%a,b2=r8_34(i5
* )%b,prq=s834,nsum=0
      r8_234(i2,i5)%a(1)=u7516_2(i2)%a(1)*r8_34(i5)%a(1)+u7516_2
     & (i2)%c(1)*s834*r8_34(i5)%b(2)
      r8_234(i2,i5)%b(1)=u7516_2(i2)%d(1)*s834*r8_34(i5)%b(1)+u7
     & 516_2(i2)%b(1)*r8_34(i5)%a(2)
      r8_234(i2,i5)%b(2)=u7516_2(i2)%b(2)*r8_34(i5)%a(1)+u7516_2
     & (i2)%d(2)*s834*r8_34(i5)%b(2)
      r8_234(i2,i5)%a(2)=u7516_2(i2)%c(2)*s834*r8_34(i5)%b(1)+u7
     & 516_2(i2)%a(2)*r8_34(i5)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r8_234(i2,i5)%a,bb=r8_234(i2,i5)%b,a1=u7516_34(i5)%a,b1=u75
* 16_34(i5)%b,c1=u7516_34(i5)%c,d1=u7516_34(i5)%d,a2=r8_2(i2)%a,b2=r8_2(
* i2)%b,prq=s82,nsum=1
      r8_234(i2,i5)%a(1)=r8_234(i2,i5)%a(1)+u7516_34(i5)%a(1)*r8
     & _2(i2)%a(1)+u7516_34(i5)%c(1)*s82*r8_2(i2)%b(2)
      r8_234(i2,i5)%b(1)=r8_234(i2,i5)%b(1)+u7516_34(i5)%d(1)*s8
     & 2*r8_2(i2)%b(1)+u7516_34(i5)%b(1)*r8_2(i2)%a(2)
      r8_234(i2,i5)%b(2)=r8_234(i2,i5)%b(2)+u7516_34(i5)%b(2)*r8
     & _2(i2)%a(1)+u7516_34(i5)%d(2)*s82*r8_2(i2)%b(2)
      r8_234(i2,i5)%a(2)=r8_234(i2,i5)%a(2)+u7516_34(i5)%c(2)*s8
     & 2*r8_2(i2)%b(1)+u7516_34(i5)%a(2)*r8_2(i2)%a(2)
      end do
      end do
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r8_156(i2,i7)%a,bb=r8_156(i2,i7)%b,a1=u7324_1(i2)%a,b1=u732
* 4_1(i2)%b,c1=u7324_1(i2)%c,d1=u7324_1(i2)%d,a2=r8_56(i7)%a,b2=r8_56(i7
* )%b,prq=s856,nsum=0
      r8_156(i2,i7)%a(1)=u7324_1(i2)%a(1)*r8_56(i7)%a(1)+u7324_1
     & (i2)%c(1)*s856*r8_56(i7)%b(2)
      r8_156(i2,i7)%b(1)=u7324_1(i2)%d(1)*s856*r8_56(i7)%b(1)+u7
     & 324_1(i2)%b(1)*r8_56(i7)%a(2)
      r8_156(i2,i7)%b(2)=u7324_1(i2)%b(2)*r8_56(i7)%a(1)+u7324_1
     & (i2)%d(2)*s856*r8_56(i7)%b(2)
      r8_156(i2,i7)%a(2)=u7324_1(i2)%c(2)*s856*r8_56(i7)%b(1)+u7
     & 324_1(i2)%a(2)*r8_56(i7)%a(2)
      end do
      end do
  
      do i2=1,2
      do i7=1,2
* TTR0 -- aa=r8_156(i2,i7)%a,bb=r8_156(i2,i7)%b,a1=u7324_56(i7)%a,b1=u73
* 24_56(i7)%b,c1=u7324_56(i7)%c,d1=u7324_56(i7)%d,a2=r8_1(i2)%a,b2=r8_1(
* i2)%b,prq=s81,nsum=1
      r8_156(i2,i7)%a(1)=r8_156(i2,i7)%a(1)+u7324_56(i7)%a(1)*r8
     & _1(i2)%a(1)+u7324_56(i7)%c(1)*s81*r8_1(i2)%b(2)
      r8_156(i2,i7)%b(1)=r8_156(i2,i7)%b(1)+u7324_56(i7)%d(1)*s8
     & 1*r8_1(i2)%b(1)+u7324_56(i7)%b(1)*r8_1(i2)%a(2)
      r8_156(i2,i7)%b(2)=r8_156(i2,i7)%b(2)+u7324_56(i7)%b(2)*r8
     & _1(i2)%a(1)+u7324_56(i7)%d(2)*s81*r8_1(i2)%b(2)
      r8_156(i2,i7)%a(2)=r8_156(i2,i7)%a(2)+u7324_56(i7)%c(2)*s8
     & 1*r8_1(i2)%b(1)+u7324_56(i7)%a(2)*r8_1(i2)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r8_134(i2,i5)%a,bb=r8_134(i2,i5)%b,a1=u7526_1(i2)%a,b1=u752
* 6_1(i2)%b,c1=u7526_1(i2)%c,d1=u7526_1(i2)%d,a2=r8_34(i5)%a,b2=r8_34(i5
* )%b,prq=s834,nsum=0
      r8_134(i2,i5)%a(1)=u7526_1(i2)%a(1)*r8_34(i5)%a(1)+u7526_1
     & (i2)%c(1)*s834*r8_34(i5)%b(2)
      r8_134(i2,i5)%b(1)=u7526_1(i2)%d(1)*s834*r8_34(i5)%b(1)+u7
     & 526_1(i2)%b(1)*r8_34(i5)%a(2)
      r8_134(i2,i5)%b(2)=u7526_1(i2)%b(2)*r8_34(i5)%a(1)+u7526_1
     & (i2)%d(2)*s834*r8_34(i5)%b(2)
      r8_134(i2,i5)%a(2)=u7526_1(i2)%c(2)*s834*r8_34(i5)%b(1)+u7
     & 526_1(i2)%a(2)*r8_34(i5)%a(2)
      end do
      end do
  
      do i2=1,2
      do i5=1,2
* TTR0 -- aa=r8_134(i2,i5)%a,bb=r8_134(i2,i5)%b,a1=u7526_34(i5)%a,b1=u75
* 26_34(i5)%b,c1=u7526_34(i5)%c,d1=u7526_34(i5)%d,a2=r8_1(i2)%a,b2=r8_1(
* i2)%b,prq=s81,nsum=1
      r8_134(i2,i5)%a(1)=r8_134(i2,i5)%a(1)+u7526_34(i5)%a(1)*r8
     & _1(i2)%a(1)+u7526_34(i5)%c(1)*s81*r8_1(i2)%b(2)
      r8_134(i2,i5)%b(1)=r8_134(i2,i5)%b(1)+u7526_34(i5)%d(1)*s8
     & 1*r8_1(i2)%b(1)+u7526_34(i5)%b(1)*r8_1(i2)%a(2)
      r8_134(i2,i5)%b(2)=r8_134(i2,i5)%b(2)+u7526_34(i5)%b(2)*r8
     & _1(i2)%a(1)+u7526_34(i5)%d(2)*s81*r8_1(i2)%b(2)
      r8_134(i2,i5)%a(2)=r8_134(i2,i5)%a(2)+u7526_34(i5)%c(2)*s8
     & 1*r8_1(i2)%b(1)+u7526_34(i5)%a(2)*r8_1(i2)%a(2)
      end do
      end do
  
* IV                                                                    
* quqd -- p=p734,q=p856
      quqd=p734(0)*p856(0)-p734(1)*p856(1)-p734(2)*p856(2)-p734(
     & 3)*p856(3)
      ccr=1d0/(f856*s12)
      ccl=1d0/(f856*s12)
      do i1=1,2
      do i2=1,2
* T0 -- qu=p734,qd=p856,v=ctrip12(i1,i2)%e,a=u734_gg(i1,i2)%a,b=u734_gg(
* i1,i2)%b,c=u734_gg(i1,i2)%c,d=u734_gg(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p734(2)*p856(3)-p856(2)*p734(3
     & ))+p734k0*(ctrip12(i1,i2)%e(2)*p856(3)-p856(2)*ctrip12(i1
     & ,i2)%e(3))-p856k0*(ctrip12(i1,i2)%e(2)*p734(3)-p734(2)*ct
     & rip12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ctrip12(i1,i2)%e(3)*p734k0+p734(3)*ctrip12(i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      ceps_2=-ctrip12(i1,i2)%e(3)*p856k0+p856(3)*ctrip12(i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=ctrip12(i1,i2)%e(0)*p734(0)-ctrip12(i1,i2)%e(1)*p734(
     & 1)-ctrip12(i1,i2)%e(2)*p734(2)-ctrip12(i1,i2)%e(3)*p734(3
     & )
      cvqd=ctrip12(i1,i2)%e(0)*p856(0)-ctrip12(i1,i2)%e(1)*p856(
     & 1)-ctrip12(i1,i2)%e(2)*p856(2)-ctrip12(i1,i2)%e(3)*p856(3
     & )
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p734k0*cvqd+p856k0*cvqu
      cauxb=-ctrip12(i1,i2)%ek0*p856(2)+p856k0*ctrip12(i1,i2)%e(
     & 2)
      cauxc=+ctrip12(i1,i2)%ek0*p734(2)-p734k0*ctrip12(i1,i2)%e(
     & 2)
      u734_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      u734_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u734_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u734_gg(i1,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u734_gg(i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      u734_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u734_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      u734_gg(i1,i2)%d(2)=ccr*ctrip12(i1,i2)%ek0
      end do
      end do
  
* quqd -- p=p756,q=p834
      quqd=p756(0)*p834(0)-p756(1)*p834(1)-p756(2)*p834(2)-p756(
     & 3)*p834(3)
      ccr=1d0/(f834*s12)
      ccl=1d0/(f834*s12)
      do i1=1,2
      do i2=1,2
* T0 -- qu=p756,qd=p834,v=ctrip12(i1,i2)%e,a=u756_gg(i1,i2)%a,b=u756_gg(
* i1,i2)%b,c=u756_gg(i1,i2)%c,d=u756_gg(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-ctrip12(i1,i2)%ek0*(p756(2)*p834(3)-p834(2)*p756(3
     & ))+p756k0*(ctrip12(i1,i2)%e(2)*p834(3)-p834(2)*ctrip12(i1
     & ,i2)%e(3))-p834k0*(ctrip12(i1,i2)%e(2)*p756(3)-p756(2)*ct
     & rip12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-ctrip12(i1,i2)%e(3)*p756k0+p756(3)*ctrip12(i1,i2)%
     & ek0
      ceps_1=ceps_1*cim
      ceps_2=-ctrip12(i1,i2)%e(3)*p834k0+p834(3)*ctrip12(i1,i2)%
     & ek0
      ceps_2=ceps_2*cim
      cvqu=ctrip12(i1,i2)%e(0)*p756(0)-ctrip12(i1,i2)%e(1)*p756(
     & 1)-ctrip12(i1,i2)%e(2)*p756(2)-ctrip12(i1,i2)%e(3)*p756(3
     & )
      cvqd=ctrip12(i1,i2)%e(0)*p834(0)-ctrip12(i1,i2)%e(1)*p834(
     & 1)-ctrip12(i1,i2)%e(2)*p834(2)-ctrip12(i1,i2)%e(3)*p834(3
     & )
      cauxa=-ctrip12(i1,i2)%ek0*quqd+p756k0*cvqd+p834k0*cvqu
      cauxb=-ctrip12(i1,i2)%ek0*p834(2)+p834k0*ctrip12(i1,i2)%e(
     & 2)
      cauxc=+ctrip12(i1,i2)%ek0*p756(2)-p756k0*ctrip12(i1,i2)%e(
     & 2)
      u756_gg(i1,i2)%a(1)=ccr*(cauxa+ceps_0)
      u756_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u756_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u756_gg(i1,i2)%b(2)=ccr*(-cauxb-ceps_2)
      u756_gg(i1,i2)%c(1)=ccr*(cauxc+ceps_1)
      u756_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u756_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      u756_gg(i1,i2)%d(2)=ccr*ctrip12(i1,i2)%ek0
      end do
      end do
  
      endif
  
  
      if (ilept(id3).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,7),a1=l3_12(i1,i2)%a,c1=l3_12(i1,i2)%c,
* a2=r4_5678(i5,i7)%a,b2=r4_5678(i5,i7)%b,prq=s312,bef=cres(i1,i2,&,i5,i
* 7,7)+,aft=
      cres(i1,i2,1,i5,i7,7)=cres(i1,i2,1,i5,i7,7)+(l3_12(i1,i2)%
     & a(1)*r4_5678(i5,i7)%a(1)+l3_12(i1,i2)%c(1)*s312*r4_5678(i
     & 5,i7)%b(2))
      cres(i1,i2,2,i5,i7,7)=cres(i1,i2,2,i5,i7,7)+(l3_12(i1,i2)%
     & c(2)*s312*r4_5678(i5,i7)%b(1)+l3_12(i1,i2)%a(2)*r4_5678(i
     & 5,i7)%a(2))
      end do
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,7),a1=l3_5678(i5,i7)%a,c1=l3_5678(i5,i7
* )%c,a2=r4_12(i1,i2)%a,b2=r4_12(i1,i2)%b,prq=s412,bef=cres(i1,i2,&,i5,i
* 7,7)+,aft=
      cres(i1,i2,1,i5,i7,7)=cres(i1,i2,1,i5,i7,7)+(l3_5678(i5,i7
     & )%a(1)*r4_12(i1,i2)%a(1)+l3_5678(i5,i7)%c(1)*s412*r4_12(i
     & 1,i2)%b(2))
      cres(i1,i2,2,i5,i7,7)=cres(i1,i2,2,i5,i7,7)+(l3_5678(i5,i7
     & )%c(2)*s412*r4_12(i1,i2)%b(1)+l3_5678(i5,i7)%a(2)*r4_12(i
     & 1,i2)%a(2))
      end do
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,7),a1=l3_156(i1,i5)%a,c1=l3_156(i1,i5)%
* c,a2=r4_278(i2,i7)%a,b2=r4_278(i2,i7)%b,prq=s3156,bef=cres(i1,i2,&,i5,
* i7,7)+,aft=
      cres(i1,i2,1,i5,i7,7)=cres(i1,i2,1,i5,i7,7)+(l3_156(i1,i5)
     & %a(1)*r4_278(i2,i7)%a(1)+l3_156(i1,i5)%c(1)*s3156*r4_278(
     & i2,i7)%b(2))
      cres(i1,i2,2,i5,i7,7)=cres(i1,i2,2,i5,i7,7)+(l3_156(i1,i5)
     & %c(2)*s3156*r4_278(i2,i7)%b(1)+l3_156(i1,i5)%a(2)*r4_278(
     & i2,i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,7),a1=l3_178(i1,i7)%a,c1=l3_178(i1,i7)%
* c,a2=r4_256(i2,i5)%a,b2=r4_256(i2,i5)%b,prq=s3178,bef=cres(i1,i2,&,i5,
* i7,7)+,aft=
      cres(i1,i2,1,i5,i7,7)=cres(i1,i2,1,i5,i7,7)+(l3_178(i1,i7)
     & %a(1)*r4_256(i2,i5)%a(1)+l3_178(i1,i7)%c(1)*s3178*r4_256(
     & i2,i5)%b(2))
      cres(i1,i2,2,i5,i7,7)=cres(i1,i2,2,i5,i7,7)+(l3_178(i1,i7)
     & %c(2)*s3178*r4_256(i2,i5)%b(1)+l3_178(i1,i7)%a(2)*r4_256(
     & i2,i5)%a(2))
      end do
      end do
      end do
      end do
  
* IV                                                                    
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i5,i1,i2)%a,cc=laux_iii(i5,i1,i2)%c,a1=l3_56(i5)%a
* ,c1=l3_56(i5)%c,a2=u356_gg(i1,i2)%a,b2=u356_gg(i1,i2)%b,c2=u356_gg(i1,
* i2)%c,d2=u356_gg(i1,i2)%d,prq=s356,nsum=0
      laux_iii(i5,i1,i2)%a(1)=l3_56(i5)%a(1)*u356_gg(i1,i2)%a(1)
     & +l3_56(i5)%c(1)*s356*u356_gg(i1,i2)%b(2)
      laux_iii(i5,i1,i2)%c(1)=l3_56(i5)%a(1)*u356_gg(i1,i2)%c(1)
     & +l3_56(i5)%c(1)*s356*u356_gg(i1,i2)%d(2)
      laux_iii(i5,i1,i2)%c(2)=l3_56(i5)%c(2)*s356*u356_gg(i1,i2)
     & %d(1)+l3_56(i5)%a(2)*u356_gg(i1,i2)%c(2)
      laux_iii(i5,i1,i2)%a(2)=l3_56(i5)%c(2)*s356*u356_gg(i1,i2)
     & %b(1)+l3_56(i5)%a(2)*u356_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,7),a1=laux_iii(i5,i1,i2)%a,c1=laux_iii(
* i5,i1,i2)%c,a2=r4_78(i7)%a,b2=r4_78(i7)%b,prq=s478,bef=cres(i1,i2,&,i5
* ,i7,7)+,aft=
      cres(i1,i2,1,i5,i7,7)=cres(i1,i2,1,i5,i7,7)+(laux_iii(i5,i
     & 1,i2)%a(1)*r4_78(i7)%a(1)+laux_iii(i5,i1,i2)%c(1)*s478*r4
     & _78(i7)%b(2))
      cres(i1,i2,2,i5,i7,7)=cres(i1,i2,2,i5,i7,7)+(laux_iii(i5,i
     & 1,i2)%c(2)*s478*r4_78(i7)%b(1)+laux_iii(i5,i1,i2)%a(2)*r4
     & _78(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i7,i1,i2)%a,cc=laux_iii(i7,i1,i2)%c,a1=l3_78(i7)%a
* ,c1=l3_78(i7)%c,a2=u378_gg(i1,i2)%a,b2=u378_gg(i1,i2)%b,c2=u378_gg(i1,
* i2)%c,d2=u378_gg(i1,i2)%d,prq=s378,nsum=0
      laux_iii(i7,i1,i2)%a(1)=l3_78(i7)%a(1)*u378_gg(i1,i2)%a(1)
     & +l3_78(i7)%c(1)*s378*u378_gg(i1,i2)%b(2)
      laux_iii(i7,i1,i2)%c(1)=l3_78(i7)%a(1)*u378_gg(i1,i2)%c(1)
     & +l3_78(i7)%c(1)*s378*u378_gg(i1,i2)%d(2)
      laux_iii(i7,i1,i2)%c(2)=l3_78(i7)%c(2)*s378*u378_gg(i1,i2)
     & %d(1)+l3_78(i7)%a(2)*u378_gg(i1,i2)%c(2)
      laux_iii(i7,i1,i2)%a(2)=l3_78(i7)%c(2)*s378*u378_gg(i1,i2)
     & %b(1)+l3_78(i7)%a(2)*u378_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,7),a1=laux_iii(i7,i1,i2)%a,c1=laux_iii(
* i7,i1,i2)%c,a2=r4_56(i5)%a,b2=r4_56(i5)%b,prq=s456,bef=cres(i1,i2,&,i5
* ,i7,7)+,aft=
      cres(i1,i2,1,i5,i7,7)=cres(i1,i2,1,i5,i7,7)+(laux_iii(i7,i
     & 1,i2)%a(1)*r4_56(i5)%a(1)+laux_iii(i7,i1,i2)%c(1)*s456*r4
     & _56(i5)%b(2))
      cres(i1,i2,2,i5,i7,7)=cres(i1,i2,2,i5,i7,7)+(laux_iii(i7,i
     & 1,i2)%c(2)*s456*r4_56(i5)%b(1)+laux_iii(i7,i1,i2)%a(2)*r4
     & _56(i5)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
  
      if (ilept(id3).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,8),a1=l3_21(i1,i2)%a,c1=l3_21(i1,i2)%c,
* a2=r4_5678(i5,i7)%a,b2=r4_5678(i5,i7)%b,prq=s312,bef=cres(i1,i2,&,i5,i
* 7,8)+,aft=
      cres(i1,i2,1,i5,i7,8)=cres(i1,i2,1,i5,i7,8)+(l3_21(i1,i2)%
     & a(1)*r4_5678(i5,i7)%a(1)+l3_21(i1,i2)%c(1)*s312*r4_5678(i
     & 5,i7)%b(2))
      cres(i1,i2,2,i5,i7,8)=cres(i1,i2,2,i5,i7,8)+(l3_21(i1,i2)%
     & c(2)*s312*r4_5678(i5,i7)%b(1)+l3_21(i1,i2)%a(2)*r4_5678(i
     & 5,i7)%a(2))
      end do
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,8),a1=l3_5678(i5,i7)%a,c1=l3_5678(i5,i7
* )%c,a2=r4_21(i1,i2)%a,b2=r4_21(i1,i2)%b,prq=s412,bef=cres(i1,i2,&,i5,i
* 7,8)+,aft=
      cres(i1,i2,1,i5,i7,8)=cres(i1,i2,1,i5,i7,8)+(l3_5678(i5,i7
     & )%a(1)*r4_21(i1,i2)%a(1)+l3_5678(i5,i7)%c(1)*s412*r4_21(i
     & 1,i2)%b(2))
      cres(i1,i2,2,i5,i7,8)=cres(i1,i2,2,i5,i7,8)+(l3_5678(i5,i7
     & )%c(2)*s412*r4_21(i1,i2)%b(1)+l3_5678(i5,i7)%a(2)*r4_21(i
     & 1,i2)%a(2))
      end do
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,8),a1=l3_256(i2,i5)%a,c1=l3_256(i2,i5)%
* c,a2=r4_178(i1,i7)%a,b2=r4_178(i1,i7)%b,prq=s3256,bef=cres(i1,i2,&,i5,
* i7,8)+,aft=
      cres(i1,i2,1,i5,i7,8)=cres(i1,i2,1,i5,i7,8)+(l3_256(i2,i5)
     & %a(1)*r4_178(i1,i7)%a(1)+l3_256(i2,i5)%c(1)*s3256*r4_178(
     & i1,i7)%b(2))
      cres(i1,i2,2,i5,i7,8)=cres(i1,i2,2,i5,i7,8)+(l3_256(i2,i5)
     & %c(2)*s3256*r4_178(i1,i7)%b(1)+l3_256(i2,i5)%a(2)*r4_178(
     & i1,i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,8),a1=l3_278(i2,i7)%a,c1=l3_278(i2,i7)%
* c,a2=r4_156(i1,i5)%a,b2=r4_156(i1,i5)%b,prq=s3278,bef=cres(i1,i2,&,i5,
* i7,8)+,aft=
      cres(i1,i2,1,i5,i7,8)=cres(i1,i2,1,i5,i7,8)+(l3_278(i2,i7)
     & %a(1)*r4_156(i1,i5)%a(1)+l3_278(i2,i7)%c(1)*s3278*r4_156(
     & i1,i5)%b(2))
      cres(i1,i2,2,i5,i7,8)=cres(i1,i2,2,i5,i7,8)+(l3_278(i2,i7)
     & %c(2)*s3278*r4_156(i1,i5)%b(1)+l3_278(i2,i7)%a(2)*r4_156(
     & i1,i5)%a(2))
      end do
      end do
      end do
      end do
  
* IV                                                                    
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i5,i1,i2)%a,cc=laux_iii(i5,i1,i2)%c,a1=l3_56(i5)%a
* ,c1=l3_56(i5)%c,a2=u356_gg(i1,i2)%a,b2=u356_gg(i1,i2)%b,c2=u356_gg(i1,
* i2)%c,d2=u356_gg(i1,i2)%d,prq=s356,nsum=0
      laux_iii(i5,i1,i2)%a(1)=l3_56(i5)%a(1)*u356_gg(i1,i2)%a(1)
     & +l3_56(i5)%c(1)*s356*u356_gg(i1,i2)%b(2)
      laux_iii(i5,i1,i2)%c(1)=l3_56(i5)%a(1)*u356_gg(i1,i2)%c(1)
     & +l3_56(i5)%c(1)*s356*u356_gg(i1,i2)%d(2)
      laux_iii(i5,i1,i2)%c(2)=l3_56(i5)%c(2)*s356*u356_gg(i1,i2)
     & %d(1)+l3_56(i5)%a(2)*u356_gg(i1,i2)%c(2)
      laux_iii(i5,i1,i2)%a(2)=l3_56(i5)%c(2)*s356*u356_gg(i1,i2)
     & %b(1)+l3_56(i5)%a(2)*u356_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,8),a1=laux_iii(i5,i1,i2)%a,c1=laux_iii(
* i5,i1,i2)%c,a2=r4_78(i7)%a,b2=r4_78(i7)%b,prq=s478,bef=cres(i1,i2,&,i5
* ,i7,8)-,aft=
      cres(i1,i2,1,i5,i7,8)=cres(i1,i2,1,i5,i7,8)-(laux_iii(i5,i
     & 1,i2)%a(1)*r4_78(i7)%a(1)+laux_iii(i5,i1,i2)%c(1)*s478*r4
     & _78(i7)%b(2))
      cres(i1,i2,2,i5,i7,8)=cres(i1,i2,2,i5,i7,8)-(laux_iii(i5,i
     & 1,i2)%c(2)*s478*r4_78(i7)%b(1)+laux_iii(i5,i1,i2)%a(2)*r4
     & _78(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i7,i1,i2)%a,cc=laux_iii(i7,i1,i2)%c,a1=l3_78(i7)%a
* ,c1=l3_78(i7)%c,a2=u378_gg(i1,i2)%a,b2=u378_gg(i1,i2)%b,c2=u378_gg(i1,
* i2)%c,d2=u378_gg(i1,i2)%d,prq=s378,nsum=0
      laux_iii(i7,i1,i2)%a(1)=l3_78(i7)%a(1)*u378_gg(i1,i2)%a(1)
     & +l3_78(i7)%c(1)*s378*u378_gg(i1,i2)%b(2)
      laux_iii(i7,i1,i2)%c(1)=l3_78(i7)%a(1)*u378_gg(i1,i2)%c(1)
     & +l3_78(i7)%c(1)*s378*u378_gg(i1,i2)%d(2)
      laux_iii(i7,i1,i2)%c(2)=l3_78(i7)%c(2)*s378*u378_gg(i1,i2)
     & %d(1)+l3_78(i7)%a(2)*u378_gg(i1,i2)%c(2)
      laux_iii(i7,i1,i2)%a(2)=l3_78(i7)%c(2)*s378*u378_gg(i1,i2)
     & %b(1)+l3_78(i7)%a(2)*u378_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,&,i5,i7,8),a1=laux_iii(i7,i1,i2)%a,c1=laux_iii(
* i7,i1,i2)%c,a2=r4_56(i5)%a,b2=r4_56(i5)%b,prq=s456,bef=cres(i1,i2,&,i5
* ,i7,8)-,aft=
      cres(i1,i2,1,i5,i7,8)=cres(i1,i2,1,i5,i7,8)-(laux_iii(i7,i
     & 1,i2)%a(1)*r4_56(i5)%a(1)+laux_iii(i7,i1,i2)%c(1)*s456*r4
     & _56(i5)%b(2))
      cres(i1,i2,2,i5,i7,8)=cres(i1,i2,2,i5,i7,8)-(laux_iii(i7,i
     & 1,i2)%c(2)*s456*r4_56(i5)%b(1)+laux_iii(i7,i1,i2)%a(2)*r4
     & _56(i5)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
  
      if (ilept(id5).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,9),a1=l5_12(i1,i2)%a,c1=l5_12(i1,i2)%c,
* a2=r6_3478(i3,i7)%a,b2=r6_3478(i3,i7)%b,prq=s512,bef=cres(i1,i2,i3,&,i
* 7,9)+,aft=
      cres(i1,i2,i3,1,i7,9)=cres(i1,i2,i3,1,i7,9)+(l5_12(i1,i2)%
     & a(1)*r6_3478(i3,i7)%a(1)+l5_12(i1,i2)%c(1)*s512*r6_3478(i
     & 3,i7)%b(2))
      cres(i1,i2,i3,2,i7,9)=cres(i1,i2,i3,2,i7,9)+(l5_12(i1,i2)%
     & c(2)*s512*r6_3478(i3,i7)%b(1)+l5_12(i1,i2)%a(2)*r6_3478(i
     & 3,i7)%a(2))
      end do
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,9),a1=l5_3478(i3,i7)%a,c1=l5_3478(i3,i7
* )%c,a2=r6_12(i1,i2)%a,b2=r6_12(i1,i2)%b,prq=s612,bef=cres(i1,i2,i3,&,i
* 7,9)+,aft=
      cres(i1,i2,i3,1,i7,9)=cres(i1,i2,i3,1,i7,9)+(l5_3478(i3,i7
     & )%a(1)*r6_12(i1,i2)%a(1)+l5_3478(i3,i7)%c(1)*s612*r6_12(i
     & 1,i2)%b(2))
      cres(i1,i2,i3,2,i7,9)=cres(i1,i2,i3,2,i7,9)+(l5_3478(i3,i7
     & )%c(2)*s612*r6_12(i1,i2)%b(1)+l5_3478(i3,i7)%a(2)*r6_12(i
     & 1,i2)%a(2))
      end do
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,9),a1=l5_134(i1,i3)%a,c1=l5_134(i1,i3)%
* c,a2=r6_278(i2,i7)%a,b2=r6_278(i2,i7)%b,prq=s5134,bef=cres(i1,i2,i3,&,
* i7,9)+,aft=
      cres(i1,i2,i3,1,i7,9)=cres(i1,i2,i3,1,i7,9)+(l5_134(i1,i3)
     & %a(1)*r6_278(i2,i7)%a(1)+l5_134(i1,i3)%c(1)*s5134*r6_278(
     & i2,i7)%b(2))
      cres(i1,i2,i3,2,i7,9)=cres(i1,i2,i3,2,i7,9)+(l5_134(i1,i3)
     & %c(2)*s5134*r6_278(i2,i7)%b(1)+l5_134(i1,i3)%a(2)*r6_278(
     & i2,i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,9),a1=l5_178(i1,i7)%a,c1=l5_178(i1,i7)%
* c,a2=r6_234(i2,i3)%a,b2=r6_234(i2,i3)%b,prq=s5178,bef=cres(i1,i2,i3,&,
* i7,9)+,aft=
      cres(i1,i2,i3,1,i7,9)=cres(i1,i2,i3,1,i7,9)+(l5_178(i1,i7)
     & %a(1)*r6_234(i2,i3)%a(1)+l5_178(i1,i7)%c(1)*s5178*r6_234(
     & i2,i3)%b(2))
      cres(i1,i2,i3,2,i7,9)=cres(i1,i2,i3,2,i7,9)+(l5_178(i1,i7)
     & %c(2)*s5178*r6_234(i2,i3)%b(1)+l5_178(i1,i7)%a(2)*r6_234(
     & i2,i3)%a(2))
      end do
      end do
      end do
      end do
  
* IV                                                                    
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i5,i1,i2)%a,cc=laux_iii(i5,i1,i2)%c,a1=l5_34(i5)%a
* ,c1=l5_34(i5)%c,a2=u534_gg(i1,i2)%a,b2=u534_gg(i1,i2)%b,c2=u534_gg(i1,
* i2)%c,d2=u534_gg(i1,i2)%d,prq=s534,nsum=0
      laux_iii(i5,i1,i2)%a(1)=l5_34(i5)%a(1)*u534_gg(i1,i2)%a(1)
     & +l5_34(i5)%c(1)*s534*u534_gg(i1,i2)%b(2)
      laux_iii(i5,i1,i2)%c(1)=l5_34(i5)%a(1)*u534_gg(i1,i2)%c(1)
     & +l5_34(i5)%c(1)*s534*u534_gg(i1,i2)%d(2)
      laux_iii(i5,i1,i2)%c(2)=l5_34(i5)%c(2)*s534*u534_gg(i1,i2)
     & %d(1)+l5_34(i5)%a(2)*u534_gg(i1,i2)%c(2)
      laux_iii(i5,i1,i2)%a(2)=l5_34(i5)%c(2)*s534*u534_gg(i1,i2)
     & %b(1)+l5_34(i5)%a(2)*u534_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,9),a1=laux_iii(i3,i1,i2)%a,c1=laux_iii(
* i3,i1,i2)%c,a2=r6_78(i7)%a,b2=r6_78(i7)%b,prq=s678,bef=cres(i1,i2,i3,&
* ,i7,9)+,aft=
      cres(i1,i2,i3,1,i7,9)=cres(i1,i2,i3,1,i7,9)+(laux_iii(i3,i
     & 1,i2)%a(1)*r6_78(i7)%a(1)+laux_iii(i3,i1,i2)%c(1)*s678*r6
     & _78(i7)%b(2))
      cres(i1,i2,i3,2,i7,9)=cres(i1,i2,i3,2,i7,9)+(laux_iii(i3,i
     & 1,i2)%c(2)*s678*r6_78(i7)%b(1)+laux_iii(i3,i1,i2)%a(2)*r6
     & _78(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i7,i1,i2)%a,cc=laux_iii(i7,i1,i2)%c,a1=l5_78(i7)%a
* ,c1=l5_78(i7)%c,a2=u578_gg(i1,i2)%a,b2=u578_gg(i1,i2)%b,c2=u578_gg(i1,
* i2)%c,d2=u578_gg(i1,i2)%d,prq=s578,nsum=0
      laux_iii(i7,i1,i2)%a(1)=l5_78(i7)%a(1)*u578_gg(i1,i2)%a(1)
     & +l5_78(i7)%c(1)*s578*u578_gg(i1,i2)%b(2)
      laux_iii(i7,i1,i2)%c(1)=l5_78(i7)%a(1)*u578_gg(i1,i2)%c(1)
     & +l5_78(i7)%c(1)*s578*u578_gg(i1,i2)%d(2)
      laux_iii(i7,i1,i2)%c(2)=l5_78(i7)%c(2)*s578*u578_gg(i1,i2)
     & %d(1)+l5_78(i7)%a(2)*u578_gg(i1,i2)%c(2)
      laux_iii(i7,i1,i2)%a(2)=l5_78(i7)%c(2)*s578*u578_gg(i1,i2)
     & %b(1)+l5_78(i7)%a(2)*u578_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,9),a1=laux_iii(i7,i1,i2)%a,c1=laux_iii(
* i7,i1,i2)%c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=cres(i1,i2,i3,&
* ,i7,9)+,aft=
      cres(i1,i2,i3,1,i7,9)=cres(i1,i2,i3,1,i7,9)+(laux_iii(i7,i
     & 1,i2)%a(1)*r6_34(i3)%a(1)+laux_iii(i7,i1,i2)%c(1)*s634*r6
     & _34(i3)%b(2))
      cres(i1,i2,i3,2,i7,9)=cres(i1,i2,i3,2,i7,9)+(laux_iii(i7,i
     & 1,i2)%c(2)*s634*r6_34(i3)%b(1)+laux_iii(i7,i1,i2)%a(2)*r6
     & _34(i3)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
  
      if (ilept(id5).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,10),a1=l5_21(i1,i2)%a,c1=l5_21(i1,i2)%c
* ,a2=r6_3478(i3,i7)%a,b2=r6_3478(i3,i7)%b,prq=s512,bef=cres(i1,i2,i3,&,
* i7,10)+,aft=
      cres(i1,i2,i3,1,i7,10)=cres(i1,i2,i3,1,i7,10)+(l5_21(i1,i2
     & )%a(1)*r6_3478(i3,i7)%a(1)+l5_21(i1,i2)%c(1)*s512*r6_3478
     & (i3,i7)%b(2))
      cres(i1,i2,i3,2,i7,10)=cres(i1,i2,i3,2,i7,10)+(l5_21(i1,i2
     & )%c(2)*s512*r6_3478(i3,i7)%b(1)+l5_21(i1,i2)%a(2)*r6_3478
     & (i3,i7)%a(2))
      end do
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,10),a1=l5_3478(i3,i7)%a,c1=l5_3478(i3,i
* 7)%c,a2=r6_21(i1,i2)%a,b2=r6_21(i1,i2)%b,prq=s612,bef=cres(i1,i2,i3,&,
* i7,10)+,aft=
      cres(i1,i2,i3,1,i7,10)=cres(i1,i2,i3,1,i7,10)+(l5_3478(i3,
     & i7)%a(1)*r6_21(i1,i2)%a(1)+l5_3478(i3,i7)%c(1)*s612*r6_21
     & (i1,i2)%b(2))
      cres(i1,i2,i3,2,i7,10)=cres(i1,i2,i3,2,i7,10)+(l5_3478(i3,
     & i7)%c(2)*s612*r6_21(i1,i2)%b(1)+l5_3478(i3,i7)%a(2)*r6_21
     & (i1,i2)%a(2))
      end do
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,10),a1=l5_234(i2,i3)%a,c1=l5_234(i2,i3)
* %c,a2=r6_178(i1,i7)%a,b2=r6_178(i1,i7)%b,prq=s5234,bef=cres(i1,i2,i3,&
* ,i7,10)+,aft=
      cres(i1,i2,i3,1,i7,10)=cres(i1,i2,i3,1,i7,10)+(l5_234(i2,i
     & 3)%a(1)*r6_178(i1,i7)%a(1)+l5_234(i2,i3)%c(1)*s5234*r6_17
     & 8(i1,i7)%b(2))
      cres(i1,i2,i3,2,i7,10)=cres(i1,i2,i3,2,i7,10)+(l5_234(i2,i
     & 3)%c(2)*s5234*r6_178(i1,i7)%b(1)+l5_234(i2,i3)%a(2)*r6_17
     & 8(i1,i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,10),a1=l5_278(i2,i7)%a,c1=l5_278(i2,i7)
* %c,a2=r6_134(i1,i3)%a,b2=r6_134(i1,i3)%b,prq=s5278,bef=cres(i1,i2,i3,&
* ,i7,10)+,aft=
      cres(i1,i2,i3,1,i7,10)=cres(i1,i2,i3,1,i7,10)+(l5_278(i2,i
     & 7)%a(1)*r6_134(i1,i3)%a(1)+l5_278(i2,i7)%c(1)*s5278*r6_13
     & 4(i1,i3)%b(2))
      cres(i1,i2,i3,2,i7,10)=cres(i1,i2,i3,2,i7,10)+(l5_278(i2,i
     & 7)%c(2)*s5278*r6_134(i1,i3)%b(1)+l5_278(i2,i7)%a(2)*r6_13
     & 4(i1,i3)%a(2))
      end do
      end do
      end do
      end do
  
* IV                                                                    
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i5,i1,i2)%a,cc=laux_iii(i5,i1,i2)%c,a1=l5_34(i5)%a
* ,c1=l5_34(i5)%c,a2=u534_gg(i1,i2)%a,b2=u534_gg(i1,i2)%b,c2=u534_gg(i1,
* i2)%c,d2=u534_gg(i1,i2)%d,prq=s534,nsum=0
      laux_iii(i5,i1,i2)%a(1)=l5_34(i5)%a(1)*u534_gg(i1,i2)%a(1)
     & +l5_34(i5)%c(1)*s534*u534_gg(i1,i2)%b(2)
      laux_iii(i5,i1,i2)%c(1)=l5_34(i5)%a(1)*u534_gg(i1,i2)%c(1)
     & +l5_34(i5)%c(1)*s534*u534_gg(i1,i2)%d(2)
      laux_iii(i5,i1,i2)%c(2)=l5_34(i5)%c(2)*s534*u534_gg(i1,i2)
     & %d(1)+l5_34(i5)%a(2)*u534_gg(i1,i2)%c(2)
      laux_iii(i5,i1,i2)%a(2)=l5_34(i5)%c(2)*s534*u534_gg(i1,i2)
     & %b(1)+l5_34(i5)%a(2)*u534_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,10),a1=laux_iii(i3,i1,i2)%a,c1=laux_iii
* (i3,i1,i2)%c,a2=r6_78(i7)%a,b2=r6_78(i7)%b,prq=s678,bef=cres(i1,i2,i3,
* &,i7,10)-,aft=
      cres(i1,i2,i3,1,i7,10)=cres(i1,i2,i3,1,i7,10)-(laux_iii(i3
     & ,i1,i2)%a(1)*r6_78(i7)%a(1)+laux_iii(i3,i1,i2)%c(1)*s678*
     & r6_78(i7)%b(2))
      cres(i1,i2,i3,2,i7,10)=cres(i1,i2,i3,2,i7,10)-(laux_iii(i3
     & ,i1,i2)%c(2)*s678*r6_78(i7)%b(1)+laux_iii(i3,i1,i2)%a(2)*
     & r6_78(i7)%a(2))
      end do
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i7,i1,i2)%a,cc=laux_iii(i7,i1,i2)%c,a1=l5_78(i7)%a
* ,c1=l5_78(i7)%c,a2=u578_gg(i1,i2)%a,b2=u578_gg(i1,i2)%b,c2=u578_gg(i1,
* i2)%c,d2=u578_gg(i1,i2)%d,prq=s578,nsum=0
      laux_iii(i7,i1,i2)%a(1)=l5_78(i7)%a(1)*u578_gg(i1,i2)%a(1)
     & +l5_78(i7)%c(1)*s578*u578_gg(i1,i2)%b(2)
      laux_iii(i7,i1,i2)%c(1)=l5_78(i7)%a(1)*u578_gg(i1,i2)%c(1)
     & +l5_78(i7)%c(1)*s578*u578_gg(i1,i2)%d(2)
      laux_iii(i7,i1,i2)%c(2)=l5_78(i7)%c(2)*s578*u578_gg(i1,i2)
     & %d(1)+l5_78(i7)%a(2)*u578_gg(i1,i2)%c(2)
      laux_iii(i7,i1,i2)%a(2)=l5_78(i7)%c(2)*s578*u578_gg(i1,i2)
     & %b(1)+l5_78(i7)%a(2)*u578_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i7=1,2
* TLTR0 -- aa=cres(i1,i2,i3,&,i7,10),a1=laux_iii(i7,i1,i2)%a,c1=laux_iii
* (i7,i1,i2)%c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=cres(i1,i2,i3,
* &,i7,10)-,aft=
      cres(i1,i2,i3,1,i7,10)=cres(i1,i2,i3,1,i7,10)-(laux_iii(i7
     & ,i1,i2)%a(1)*r6_34(i3)%a(1)+laux_iii(i7,i1,i2)%c(1)*s634*
     & r6_34(i3)%b(2))
      cres(i1,i2,i3,2,i7,10)=cres(i1,i2,i3,2,i7,10)-(laux_iii(i7
     & ,i1,i2)%c(2)*s634*r6_34(i3)%b(1)+laux_iii(i7,i1,i2)%a(2)*
     & r6_34(i3)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
  
      if (ilept(id7).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,11),a1=l7_12(i1,i2)%a,c1=l7_12(i1,i2)%c
* ,a2=r8_3456(i3,i5)%a,b2=r8_3456(i3,i5)%b,prq=s712,bef=cres(i1,i2,i3,i5
* ,&,11)+,aft=
      cres(i1,i2,i3,i5,1,11)=cres(i1,i2,i3,i5,1,11)+(l7_12(i1,i2
     & )%a(1)*r8_3456(i3,i5)%a(1)+l7_12(i1,i2)%c(1)*s712*r8_3456
     & (i3,i5)%b(2))
      cres(i1,i2,i3,i5,2,11)=cres(i1,i2,i3,i5,2,11)+(l7_12(i1,i2
     & )%c(2)*s712*r8_3456(i3,i5)%b(1)+l7_12(i1,i2)%a(2)*r8_3456
     & (i3,i5)%a(2))
      end do
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,11),a1=l7_3456(i3,i5)%a,c1=l7_3456(i3,i
* 5)%c,a2=r8_12(i1,i2)%a,b2=r8_12(i1,i2)%b,prq=s812,bef=cres(i1,i2,i3,i5
* ,&,11)+,aft=
      cres(i1,i2,i3,i5,1,11)=cres(i1,i2,i3,i5,1,11)+(l7_3456(i3,
     & i5)%a(1)*r8_12(i1,i2)%a(1)+l7_3456(i3,i5)%c(1)*s812*r8_12
     & (i1,i2)%b(2))
      cres(i1,i2,i3,i5,2,11)=cres(i1,i2,i3,i5,2,11)+(l7_3456(i3,
     & i5)%c(2)*s812*r8_12(i1,i2)%b(1)+l7_3456(i3,i5)%a(2)*r8_12
     & (i1,i2)%a(2))
      end do
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,11),a1=l7_134(i1,i3)%a,c1=l7_134(i1,i3)
* %c,a2=r8_256(i2,i5)%a,b2=r8_256(i2,i5)%b,prq=s7134,bef=cres(i1,i2,i3,i
* 5,&,11)+,aft=
      cres(i1,i2,i3,i5,1,11)=cres(i1,i2,i3,i5,1,11)+(l7_134(i1,i
     & 3)%a(1)*r8_256(i2,i5)%a(1)+l7_134(i1,i3)%c(1)*s7134*r8_25
     & 6(i2,i5)%b(2))
      cres(i1,i2,i3,i5,2,11)=cres(i1,i2,i3,i5,2,11)+(l7_134(i1,i
     & 3)%c(2)*s7134*r8_256(i2,i5)%b(1)+l7_134(i1,i3)%a(2)*r8_25
     & 6(i2,i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,11),a1=l7_156(i1,i5)%a,c1=l7_156(i1,i5)
* %c,a2=r8_234(i2,i3)%a,b2=r8_234(i2,i3)%b,prq=s7156,bef=cres(i1,i2,i3,i
* 5,&,11)+,aft=
      cres(i1,i2,i3,i5,1,11)=cres(i1,i2,i3,i5,1,11)+(l7_156(i1,i
     & 5)%a(1)*r8_234(i2,i3)%a(1)+l7_156(i1,i5)%c(1)*s7156*r8_23
     & 4(i2,i3)%b(2))
      cres(i1,i2,i3,i5,2,11)=cres(i1,i2,i3,i5,2,11)+(l7_156(i1,i
     & 5)%c(2)*s7156*r8_234(i2,i3)%b(1)+l7_156(i1,i5)%a(2)*r8_23
     & 4(i2,i3)%a(2))
      end do
      end do
      end do
      end do
  
* IV                                                                    
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i5,i1,i2)%a,cc=laux_iii(i5,i1,i2)%c,a1=l7_34(i5)%a
* ,c1=l7_34(i5)%c,a2=u734_gg(i1,i2)%a,b2=u734_gg(i1,i2)%b,c2=u734_gg(i1,
* i2)%c,d2=u734_gg(i1,i2)%d,prq=s734,nsum=0
      laux_iii(i5,i1,i2)%a(1)=l7_34(i5)%a(1)*u734_gg(i1,i2)%a(1)
     & +l7_34(i5)%c(1)*s734*u734_gg(i1,i2)%b(2)
      laux_iii(i5,i1,i2)%c(1)=l7_34(i5)%a(1)*u734_gg(i1,i2)%c(1)
     & +l7_34(i5)%c(1)*s734*u734_gg(i1,i2)%d(2)
      laux_iii(i5,i1,i2)%c(2)=l7_34(i5)%c(2)*s734*u734_gg(i1,i2)
     & %d(1)+l7_34(i5)%a(2)*u734_gg(i1,i2)%c(2)
      laux_iii(i5,i1,i2)%a(2)=l7_34(i5)%c(2)*s734*u734_gg(i1,i2)
     & %b(1)+l7_34(i5)%a(2)*u734_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,11),a1=laux_iii(i3,i1,i2)%a,c1=laux_iii
* (i3,i1,i2)%c,a2=r8_56(i5)%a,b2=r8_56(i5)%b,prq=s856,bef=cres(i1,i2,i3,
* i5,&,11)+,aft=
      cres(i1,i2,i3,i5,1,11)=cres(i1,i2,i3,i5,1,11)+(laux_iii(i3
     & ,i1,i2)%a(1)*r8_56(i5)%a(1)+laux_iii(i3,i1,i2)%c(1)*s856*
     & r8_56(i5)%b(2))
      cres(i1,i2,i3,i5,2,11)=cres(i1,i2,i3,i5,2,11)+(laux_iii(i3
     & ,i1,i2)%c(2)*s856*r8_56(i5)%b(1)+laux_iii(i3,i1,i2)%a(2)*
     & r8_56(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i7,i1,i2)%a,cc=laux_iii(i7,i1,i2)%c,a1=l7_56(i7)%a
* ,c1=l7_56(i7)%c,a2=u756_gg(i1,i2)%a,b2=u756_gg(i1,i2)%b,c2=u756_gg(i1,
* i2)%c,d2=u756_gg(i1,i2)%d,prq=s756,nsum=0
      laux_iii(i7,i1,i2)%a(1)=l7_56(i7)%a(1)*u756_gg(i1,i2)%a(1)
     & +l7_56(i7)%c(1)*s756*u756_gg(i1,i2)%b(2)
      laux_iii(i7,i1,i2)%c(1)=l7_56(i7)%a(1)*u756_gg(i1,i2)%c(1)
     & +l7_56(i7)%c(1)*s756*u756_gg(i1,i2)%d(2)
      laux_iii(i7,i1,i2)%c(2)=l7_56(i7)%c(2)*s756*u756_gg(i1,i2)
     & %d(1)+l7_56(i7)%a(2)*u756_gg(i1,i2)%c(2)
      laux_iii(i7,i1,i2)%a(2)=l7_56(i7)%c(2)*s756*u756_gg(i1,i2)
     & %b(1)+l7_56(i7)%a(2)*u756_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,11),a1=laux_iii(i5,i1,i2)%a,c1=laux_iii
* (i5,i1,i2)%c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=cres(i1,i2,i3,
* i5,&,11)+,aft=
      cres(i1,i2,i3,i5,1,11)=cres(i1,i2,i3,i5,1,11)+(laux_iii(i5
     & ,i1,i2)%a(1)*r8_34(i3)%a(1)+laux_iii(i5,i1,i2)%c(1)*s834*
     & r8_34(i3)%b(2))
      cres(i1,i2,i3,i5,2,11)=cres(i1,i2,i3,i5,2,11)+(laux_iii(i5
     & ,i1,i2)%c(2)*s834*r8_34(i3)%b(1)+laux_iii(i5,i1,i2)%a(2)*
     & r8_34(i3)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
  
      if (ilept(id7).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,12),a1=l7_21(i1,i2)%a,c1=l7_21(i1,i2)%c
* ,a2=r8_3456(i3,i5)%a,b2=r8_3456(i3,i5)%b,prq=s712,bef=cres(i1,i2,i3,i5
* ,&,12)+,aft=
      cres(i1,i2,i3,i5,1,12)=cres(i1,i2,i3,i5,1,12)+(l7_21(i1,i2
     & )%a(1)*r8_3456(i3,i5)%a(1)+l7_21(i1,i2)%c(1)*s712*r8_3456
     & (i3,i5)%b(2))
      cres(i1,i2,i3,i5,2,12)=cres(i1,i2,i3,i5,2,12)+(l7_21(i1,i2
     & )%c(2)*s712*r8_3456(i3,i5)%b(1)+l7_21(i1,i2)%a(2)*r8_3456
     & (i3,i5)%a(2))
      end do
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,12),a1=l7_3456(i3,i5)%a,c1=l7_3456(i3,i
* 5)%c,a2=r8_21(i1,i2)%a,b2=r8_21(i1,i2)%b,prq=s812,bef=cres(i1,i2,i3,i5
* ,&,12)+,aft=
      cres(i1,i2,i3,i5,1,12)=cres(i1,i2,i3,i5,1,12)+(l7_3456(i3,
     & i5)%a(1)*r8_21(i1,i2)%a(1)+l7_3456(i3,i5)%c(1)*s812*r8_21
     & (i1,i2)%b(2))
      cres(i1,i2,i3,i5,2,12)=cres(i1,i2,i3,i5,2,12)+(l7_3456(i3,
     & i5)%c(2)*s812*r8_21(i1,i2)%b(1)+l7_3456(i3,i5)%a(2)*r8_21
     & (i1,i2)%a(2))
      end do
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,12),a1=l7_234(i2,i3)%a,c1=l7_234(i2,i3)
* %c,a2=r8_156(i1,i5)%a,b2=r8_156(i1,i5)%b,prq=s7234,bef=cres(i1,i2,i3,i
* 5,&,12)+,aft=
      cres(i1,i2,i3,i5,1,12)=cres(i1,i2,i3,i5,1,12)+(l7_234(i2,i
     & 3)%a(1)*r8_156(i1,i5)%a(1)+l7_234(i2,i3)%c(1)*s7234*r8_15
     & 6(i1,i5)%b(2))
      cres(i1,i2,i3,i5,2,12)=cres(i1,i2,i3,i5,2,12)+(l7_234(i2,i
     & 3)%c(2)*s7234*r8_156(i1,i5)%b(1)+l7_234(i2,i3)%a(2)*r8_15
     & 6(i1,i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,12),a1=l7_256(i2,i5)%a,c1=l7_256(i2,i5)
* %c,a2=r8_134(i1,i3)%a,b2=r8_134(i1,i3)%b,prq=s7256,bef=cres(i1,i2,i3,i
* 5,&,12)+,aft=
      cres(i1,i2,i3,i5,1,12)=cres(i1,i2,i3,i5,1,12)+(l7_256(i2,i
     & 5)%a(1)*r8_134(i1,i3)%a(1)+l7_256(i2,i5)%c(1)*s7256*r8_13
     & 4(i1,i3)%b(2))
      cres(i1,i2,i3,i5,2,12)=cres(i1,i2,i3,i5,2,12)+(l7_256(i2,i
     & 5)%c(2)*s7256*r8_134(i1,i3)%b(1)+l7_256(i2,i5)%a(2)*r8_13
     & 4(i1,i3)%a(2))
      end do
      end do
      end do
      end do
  
* IV                                                                    
      do i5=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i5,i1,i2)%a,cc=laux_iii(i5,i1,i2)%c,a1=l7_34(i5)%a
* ,c1=l7_34(i5)%c,a2=u734_gg(i1,i2)%a,b2=u734_gg(i1,i2)%b,c2=u734_gg(i1,
* i2)%c,d2=u734_gg(i1,i2)%d,prq=s734,nsum=0
      laux_iii(i5,i1,i2)%a(1)=l7_34(i5)%a(1)*u734_gg(i1,i2)%a(1)
     & +l7_34(i5)%c(1)*s734*u734_gg(i1,i2)%b(2)
      laux_iii(i5,i1,i2)%c(1)=l7_34(i5)%a(1)*u734_gg(i1,i2)%c(1)
     & +l7_34(i5)%c(1)*s734*u734_gg(i1,i2)%d(2)
      laux_iii(i5,i1,i2)%c(2)=l7_34(i5)%c(2)*s734*u734_gg(i1,i2)
     & %d(1)+l7_34(i5)%a(2)*u734_gg(i1,i2)%c(2)
      laux_iii(i5,i1,i2)%a(2)=l7_34(i5)%c(2)*s734*u734_gg(i1,i2)
     & %b(1)+l7_34(i5)%a(2)*u734_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,12),a1=laux_iii(i3,i1,i2)%a,c1=laux_iii
* (i3,i1,i2)%c,a2=r8_56(i5)%a,b2=r8_56(i5)%b,prq=s856,bef=cres(i1,i2,i3,
* i5,&,12)-,aft=
      cres(i1,i2,i3,i5,1,12)=cres(i1,i2,i3,i5,1,12)-(laux_iii(i3
     & ,i1,i2)%a(1)*r8_56(i5)%a(1)+laux_iii(i3,i1,i2)%c(1)*s856*
     & r8_56(i5)%b(2))
      cres(i1,i2,i3,i5,2,12)=cres(i1,i2,i3,i5,2,12)-(laux_iii(i3
     & ,i1,i2)%c(2)*s856*r8_56(i5)%b(1)+laux_iii(i3,i1,i2)%a(2)*
     & r8_56(i5)%a(2))
      end do
      end do
      end do
      end do
  
      do i7=1,2
      do i1=1,2
      do i2=1,2
* TLT0 -- aa=laux_iii(i7,i1,i2)%a,cc=laux_iii(i7,i1,i2)%c,a1=l7_56(i7)%a
* ,c1=l7_56(i7)%c,a2=u756_gg(i1,i2)%a,b2=u756_gg(i1,i2)%b,c2=u756_gg(i1,
* i2)%c,d2=u756_gg(i1,i2)%d,prq=s756,nsum=0
      laux_iii(i7,i1,i2)%a(1)=l7_56(i7)%a(1)*u756_gg(i1,i2)%a(1)
     & +l7_56(i7)%c(1)*s756*u756_gg(i1,i2)%b(2)
      laux_iii(i7,i1,i2)%c(1)=l7_56(i7)%a(1)*u756_gg(i1,i2)%c(1)
     & +l7_56(i7)%c(1)*s756*u756_gg(i1,i2)%d(2)
      laux_iii(i7,i1,i2)%c(2)=l7_56(i7)%c(2)*s756*u756_gg(i1,i2)
     & %d(1)+l7_56(i7)%a(2)*u756_gg(i1,i2)%c(2)
      laux_iii(i7,i1,i2)%a(2)=l7_56(i7)%c(2)*s756*u756_gg(i1,i2)
     & %b(1)+l7_56(i7)%a(2)*u756_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i5=1,2
* TLTR0 -- aa=cres(i1,i2,i3,i5,&,12),a1=laux_iii(i5,i1,i2)%a,c1=laux_iii
* (i5,i1,i2)%c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=cres(i1,i2,i3,
* i5,&,12)-,aft=
      cres(i1,i2,i3,i5,1,12)=cres(i1,i2,i3,i5,1,12)-(laux_iii(i5
     & ,i1,i2)%a(1)*r8_34(i3)%a(1)+laux_iii(i5,i1,i2)%c(1)*s834*
     & r8_34(i3)%b(2))
      cres(i1,i2,i3,i5,2,12)=cres(i1,i2,i3,i5,2,12)-(laux_iii(i5
     & ,i1,i2)%c(2)*s834*r8_34(i3)%b(1)+laux_iii(i5,i1,i2)%a(2)*
     & r8_34(i3)%a(2))
      end do
      end do
      end do
      end do
  
      endif
  
  
  
  
      spk0=sqrt(abs(p3k0*p4k0*p5k0*p6k0*p7k0*p8k0))
  
* final result for color 1-6                                            
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,1)=cres(i1,i2,i3,i5,i7,1)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,1)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,2)=cres(i1,i2,i3,i5,i7,2)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,2)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,3)=cres(i1,i2,i3,i5,i7,3)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,3)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,4)=cres(i1,i2,i3,i5,i7,4)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,4)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,5)=cres(i1,i2,i3,i5,i7,5)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,5)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,6)=cres(i1,i2,i3,i5,i7,6)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,6)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
* final result for color 7-12                                           
  
      if (ilept(id3).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,7)=cres(i1,i2,i3,i5,i7,7)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,7)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id3).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,8)=cres(i1,i2,i3,i5,i7,8)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,8)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id5).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,9)=cres(i1,i2,i3,i5,i7,9)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,9)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id5).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,10)=cres(i1,i2,i3,i5,i7,10)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,10)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,11)=cres(i1,i2,i3,i5,i7,11)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,11)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
      if (ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,12)=cres(i1,i2,i3,i5,i7,12)/spk0
  
         enddo
         enddo
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
         do i5=1,2
         do i7=1,2
            cres(i1,i2,i3,i5,i7,12)=czero
  
         enddo
         enddo
         enddo
         enddo
         enddo
  
      endif
  
  
      return
      end
