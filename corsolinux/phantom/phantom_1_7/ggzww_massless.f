************************************************************************
*                                                                       
* This routine computes the basic amplitude corresponding to 2 gluons   
* plus 6 outgoing  fermions which can form a Z and 2 W-.                
* The 2 fermions which can form the Z are massive, the other massless.  
* All outgoing particles are considered different.                      
* The input particles are ordered in such a way that 1 and 2 are the glu
* From 3 to 8 odd are particles, even antiparticles.                    
* 34 corresponds to Z, 56 to W+, 78 to W-.                              
* p1, ...p8 are all outgoing momenta                                    
* id1....id8 give the identities of the outgoing particles              
* cres  is on output the complex helicity amplitude                     
* cres=cres(2,2,2,2,12). The first two indeces refer to the             
* helicities of the gluons, the third and fourth one to the chiralities 
* particles 3 and 4, while all other fermions have only chirality index 
* (corresponding to -)                                                  
* cres(i1,i2,i3,i4,1)  --> T[3,1,4]T[5,2,6][7,8]                        
* cres(i1,i2,i3,i4,2)  --> T[3,2,4]T[5,1,6][7,8]                        
* cres(i1,i2,i3,i4,3)  --> T[3,1,4]T[5,6][7,2,8]                        
* cres(i1,i2,i3,i4,4)  --> T[3,2,4]T[5,6][7,1,8]                        
* cres(i1,i2,i3,i4,5)  --> T[3,4]T[5,1,6][7,2,8]                        
* cres(i1,i2,i3,i4,6)  --> T[3,4]T[5,2,6][7,1,8]                        
* cres(i1,i2,i3,i4,7)  --> T[3,1,2,4]T[5,6][7,8]                        
* cres(i1,i2,i3,i4,8)  --> T[3,2,1,4]T[5,6][7,8]                        
* cres(i1,i2,i3,i4,9)  --> T[3,4]T[5,1,2,6][7,8]                        
* cres(i1,i2,i3,i4,10) --> T[3,4]T[5,2,1,6][7,8]                        
* cres(i1,i2,i3,i4,11) --> T[3,4]T[5,6][7,1,2,8]                        
* cres(i1,i2,i3,i4,12) --> T[3,4]T[5,6][7,2,1,8]                        
*                                                                       
************************************************************************
* Comments:  "W" and "Z"  stay for, respectivecly, a couple of quark    
* which can form a W  and a couple of quarks which can form a Z         
* Correpsondingly Wline and Zline stay for a fermion line which         
* ends with two quarks which can form a W or a Z                        
************************************************************************
      subroutine ggzww_massless(p1,p2,p3,p4,p5,p6,p7,p8,
     &              id1,id2,id3,id4,id5,id6,id7,id8,cres)
  
      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)
  
      include 'phact_data_types_massless.inc'
      include 'common.h'
  
      dimension cres(2,2,2,12), cres_aux(2,2,2,12)
  
* VARIABLE DEFINITIONS                                                  
*four momenta                                                           
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3)
     &         ,p5(0:3),p6(0:3),p7(0:3),p8(0:3)
* auxiliaries                                                           
      type(tl0) laux_imu(2,0:3),laux_iii(2,2,2)
  
* -> for forks                                                          
*    -simple fork                                                       
      dimension p56(0:3)
      type(polcom) cw56
      dimension p78(0:3)
      type(polcom) cw78
      dimension p34(0:3)
      type(polcom) cz34(2),cf34(2)
  
      dimension p12(0:3)
      type(polcom) ce1(2),ce2(2),ctrip12(2,2)
  
*    -one-gluon fork                                                    
      dimension p51(0:3),p61(0:3),p516(0:3)
      type(polcom) cw516(2)
      type(twl0) l5_1(2),lw561(0:3)
      type(twr0) r6_1(2),rw651(0:3)
      dimension p52(0:3),p62(0:3),p526(0:3)
      type(polcom) cw526(2)
      type(twl0) l5_2(2),lw562(0:3)
      type(twr0) r6_2(2),rw652(0:3)
      dimension p71(0:3),p81(0:3),p718(0:3)
      type(polcom) cw718(2)
      type(twl0) l7_1(2),lw781(0:3)
      type(twr0) r8_1(2),rw871(0:3)
      dimension p72(0:3),p82(0:3),p728(0:3)
      type(polcom) cw728(2)
      type(twl0) l7_2(2),lw782(0:3)
      type(twr0) r8_2(2),rw872(0:3)
  
      dimension p31(0:3),p41(0:3),p314(0:3)
      type(polcom) cz314(2,2),cf314(2,2)
      type(tl0) l3_1(2),lz341(0:3),lf341(0:3)
      type(tr0) r4_1(2),rz431(0:3),rf431(0:3)
      dimension p32(0:3),p42(0:3),p324(0:3)
      type(polcom) cz324(2,2),cf324(2,2)
      type(tl0) l3_2(2),lz342(0:3),lf342(0:3)
      type(tr0) r4_2(2),rz432(0:3),rf432(0:3)
  
*     -two-gluons fork                                                  
      dimension p512(0:3),p612(0:3),p5126(0:3)
      type(twl0) l5_gg(2,2),lw5612(0:3)
      type(twr0) r6_gg(2,2),rw6512(0:3)
      type(polcom) cw5126(2,2)
      type(tw0) u51_2(2),u61_2(2),uw51(0:3)
      type(twl0) l5_12(2,2)
      type(twr0) r6_12(2,2)
      type(polcom) cw5216(2,2)
      type(tw0) u52_1(2),u62_1(2),uw52(0:3)
      type(twl0) l5_21(2,2)
      type(twr0) r6_21(2,2)
      dimension p712(0:3),p812(0:3),p7128(0:3)
      type(twl0) l7_gg(2,2),lw7812(0:3)
      type(twr0) r8_gg(2,2),rw8712(0:3)
      type(polcom) cw7128(2,2)
      type(tw0) u71_2(2),u81_2(2),uw71(0:3)
      type(twl0) l7_12(2,2)
      type(twr0) r8_12(2,2)
      type(polcom) cw7218(2,2)
      type(tw0) u72_1(2),u82_1(2),uw72(0:3)
      type(twl0) l7_21(2,2)
      type(twr0) r8_21(2,2)
  
      dimension p312(0:3),p412(0:3),p3124(0:3)
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
  
* -> lines without gluon                                                
      dimension p356(0:3),p456(0:3)
      type(twl0) l3_56
      type(twr0) r4_56
  
      dimension p3156(0:3)
      type(twl0) l3_516(2),l3_5126(2,2)
      type(twr0) r4_516(2),r4_5126(2,2)
  
      dimension p3256(0:3)
      type(twl0) l3_526(2),l3_5216(2,2)
      type(twr0) r4_526(2),r4_5216(2,2)
      dimension p378(0:3),p478(0:3)
      type(twl0) l3_78
      type(twr0) r4_78
  
      dimension p3178(0:3)
      type(twl0) l3_718(2),l3_7128(2,2)
      type(twr0) r4_718(2),r4_7128(2,2)
  
      dimension p3278(0:3)
      type(twl0) l3_728(2),l3_7218(2,2)
      type(twr0) r4_728(2),r4_7218(2,2)
      dimension p534(0:3),p634(0:3)
      type(twl0) l5_34(2)
      type(twr0) r6_34(2)
  
      dimension p5134(0:3)
      type(twl0) l5_314(2,2),l5_3124(2,2,2)
      type(twr0) r6_314(2,2),r6_3124(2,2,2)
  
      dimension p5234(0:3)
      type(twl0) l5_324(2,2),l5_3214(2,2,2)
      type(twr0) r6_324(2,2),r6_3214(2,2,2)
      dimension p578(0:3),p678(0:3)
      type(twl0) l5_78
      type(twr0) r6_78
  
      dimension p5178(0:3)
      type(twl0) l5_718(2),l5_7128(2,2)
      type(twr0) r6_718(2),r6_7128(2,2)
  
      dimension p5278(0:3)
      type(twl0) l5_728(2),l5_7218(2,2)
      type(twr0) r6_728(2),r6_7218(2,2)
      dimension p734(0:3),p834(0:3)
      type(twl0) l7_34(2)
      type(twr0) r8_34(2)
  
      dimension p7134(0:3)
      type(twl0) l7_314(2,2),l7_3124(2,2,2)
      type(twr0) r8_314(2,2),r8_3124(2,2,2)
  
      dimension p7234(0:3)
      type(twl0) l7_324(2,2),l7_3214(2,2,2)
      type(twr0) r8_324(2,2),r8_3214(2,2,2)
      dimension p756(0:3),p856(0:3)
      type(twl0) l7_56
      type(twr0) r8_56
  
      dimension p7156(0:3)
      type(twl0) l7_516(2),l7_5126(2,2)
      type(twr0) r8_516(2),r8_5126(2,2)
  
      dimension p7256(0:3)
      type(twl0) l7_526(2),l7_5216(2,2)
      type(twr0) r8_526(2),r8_5216(2,2)
  
* -> lines with one gluon                                               
      type(twl0) l3_156(2),l3_1526(2,2),
     &  l3_56718(2)
      type(tw0) u31_56,u31_526(2),u356_1(2),
     &  u3516_2(2),u378_516(2),u3718_56
  
      type(twl0) l3_256(2),l3_2516(2,2),
     &  l3_56728(2)
      type(tw0) u32_56,u32_516(2),u356_2(2),
     &  u3526_1(2),u378_526(2),u3728_56
  
      type(twl0) l3_178(2),l3_1728(2,2),
     &  l3_78516(2)
      type(tw0) u31_78,u31_728(2),u378_1(2),
     &  u3718_2(2),u356_718(2),u3516_78
  
      type(twl0) l3_278(2),l3_2718(2,2),
     &  l3_78526(2)
      type(tw0) u32_78,u32_718(2),u378_2(2),
     &  u3728_1(2),u356_728(2),u3526_78
  
      type(twl0) l5_134(2,2),l5_1324(2,2,2),
     &  l5_34718(2,2)
      type(tw0) u51_34(2),u51_324(2,2),u534_1(2),
     &  u5314_2(2),u578_314(2,2),u5718_34(2)
  
      type(twl0) l5_234(2,2),l5_2314(2,2,2),
     &  l5_34728(2,2)
      type(tw0) u52_34(2),u52_314(2,2),u534_2(2),
     &  u5324_1(2),u578_324(2,2),u5728_34(2)
  
      type(twl0) l5_178(2),l5_1728(2,2),
     &  l5_78314(2,2)
      type(tw0) u51_78,u51_728(2),u578_1(2),
     &  u5718_2(2),u534_718(2),u5314_78
  
      type(twl0) l5_278(2),l5_2718(2,2),
     &  l5_78324(2,2)
      type(tw0) u52_78,u52_718(2),u578_2(2),
     &  u5728_1(2),u534_728(2),u5324_78
  
      type(twl0) l7_134(2,2),l7_1324(2,2,2),
     &  l7_34516(2,2)
      type(tw0) u71_34(2),u71_324(2,2),u734_1(2),
     &  u7314_2(2),u756_314(2,2),u7516_34(2)
  
      type(twl0) l7_234(2,2),l7_2314(2,2,2),
     &  l7_34526(2,2)
      type(tw0) u72_34(2),u72_314(2,2),u734_2(2),
     &  u7324_1(2),u756_324(2,2),u7526_34(2)
  
      type(twl0) l7_156(2),l7_1526(2,2),
     &  l7_56314(2,2)
      type(tw0) u71_56,u71_526(2),u756_1(2),
     &  u7516_2(2),u734_516(2),u7314_56
  
      type(twl0) l7_256(2),l7_2516(2,2),
     &  l7_56324(2,2)
      type(tw0) u72_56,u72_516(2),u756_2(2),
     &  u7526_1(2),u734_526(2),u7324_56
  
  
* -> lines with two gluons                                              
      type(tw0) u312_56,u312_78,u356_78,u378_56,
     &  u356_gg(2,2),u378_gg(2,2)
      type(twr0) r4_5678,r4_156(2),r4_256(2),
     &  r4_178(2),r4_278(2)
      type(twl0) l3_5678
      type(tw0) u512_34(2),u512_78,u534_78,u578_34(2),
     &  u534_gg(2,2),u578_gg(2,2)
      type(twr0) r6_3478(2),r6_134(2,2),r6_234(2,2),
     &  r6_178(2),r6_278(2)
      type(twl0) l5_3478(2)
      type(tw0) u712_34(2),u712_56,u734_56,u756_34(2),
     &  u734_gg(2,2),u756_gg(2,2)
      type(twr0) r8_3456(2),r8_134(2,2),r8_234(2,2),
     &  r8_156(2),r8_256(2)
      type(twl0) l7_3456(2)
  
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
  
* gluon polarizations                                                   
  
* pol -- p=p1,m=0,e1=ce1(1)%e,e2=ce1(2)%e,e3=0,nhz=0                    
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
* triple vertex gluon                                                   
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
      do mu=0,3
       ctrip12(i1,i2)%e(mu)=ctrip12(i1,i2)%e(mu)/s12
      enddo
      enddo
      enddo
  
      do i1=1,2
      do i2=1,2
* pk0 -- p=ctrip12(i1,i2)%e
      ctrip12(i1,i2)%ek0=ctrip12(i1,i2)%e(0)-ctrip12(i1,i2)%e(1)
      end do
      end do
  
* -> for forks                                                          
*    -simple fork ----------------------------------------------------  
  
* for triple vertex                                                     
      do mu=0,3
       p56(mu)=p5(mu)+p6(mu)
      enddo
  
* quqd -- p=p5,q=p6
      quqd=p5(0)*p6(0)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3)
      s56=2.d0*quqd
  
      ccl=wcl/(-s56+cmw2)
* TW10 -- qu=p5,qd=p6,v=0,a=cw56%e(0),cl=ccl,nsum=0
      eps_0=-p5(2)*p6(3)+p6(2)*p5(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p5k0*p6(0)+p6k0*p5(0)
      cw56%e(0)=ccl*(auxa-ceps_0)
* TW10 -- qu=p5,qd=p6,v=1,a=cw56%e(1),cl=ccl,nsum=0
      auxa=-quqd+p5k0*p6(1)+p6k0*p5(1)
      cw56%e(1)=ccl*(auxa-ceps_0)
* TW10 -- qu=p5,qd=p6,v=2,a=cw56%e(2),cl=ccl,nsum=0
      eps_0=-p5k0*p6(3)+p6k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(2)+p6k0*p5(2)
      cw56%e(2)=ccl*(auxa-ceps_0)
* TW10 -- qu=p5,qd=p6,v=3,a=cw56%e(3),cl=ccl,nsum=0
      eps_0=p5k0*p6(2)-p6k0*p5(2)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(3)+p6k0*p5(3)
      cw56%e(3)=ccl*(auxa-ceps_0)
  
* pk0 -- p=cw56%e
      cw56%ek0=cw56%e(0)-cw56%e(1)
  
  
* for triple vertex                                                     
      do mu=0,3
       p78(mu)=p7(mu)+p8(mu)
      enddo
  
* quqd -- p=p7,q=p8
      quqd=p7(0)*p8(0)-p7(1)*p8(1)-p7(2)*p8(2)-p7(3)*p8(3)
      s78=2.d0*quqd
  
      ccl=wcl/(-s78+cmw2)
* TW10 -- qu=p7,qd=p8,v=0,a=cw78%e(0),cl=ccl,nsum=0
      eps_0=-p7(2)*p8(3)+p8(2)*p7(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p7k0*p8(0)+p8k0*p7(0)
      cw78%e(0)=ccl*(auxa-ceps_0)
* TW10 -- qu=p7,qd=p8,v=1,a=cw78%e(1),cl=ccl,nsum=0
      auxa=-quqd+p7k0*p8(1)+p8k0*p7(1)
      cw78%e(1)=ccl*(auxa-ceps_0)
* TW10 -- qu=p7,qd=p8,v=2,a=cw78%e(2),cl=ccl,nsum=0
      eps_0=-p7k0*p8(3)+p8k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(2)+p8k0*p7(2)
      cw78%e(2)=ccl*(auxa-ceps_0)
* TW10 -- qu=p7,qd=p8,v=3,a=cw78%e(3),cl=ccl,nsum=0
      eps_0=p7k0*p8(2)-p8k0*p7(2)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(3)+p8k0*p7(3)
      cw78%e(3)=ccl*(auxa-ceps_0)
  
* pk0 -- p=cw78%e
      cw78%ek0=cw78%e(0)-cw78%e(1)
  
  
* for triple vertex                                                     
      do mu=0,3
       p34(mu)=p3(mu)+p4(mu)
      enddo
  
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
  
*    -one-gluon fork -------------------------------------------------  
  
  
      do mu=0,3
         p51(mu)=p5(mu)+p1(mu)
      enddo
  
* for triple vertex                                                     
      do mu=0,3
       p516(mu)=p51(mu)+p6(mu)
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
      ccl=1.d0/(f51)
      do i1=1,2
* TWL0 -- qu=p5,qd=p51,v=ce1(i1)%e,a=l5_1(i1)%a,c=l5_1(i1)%c,cl=ccl,nsum
* =0
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
      l5_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s61/2d0
      do i1=1,2
* TWR0 -- qu=p61,qd=p6,v=ce1(i1)%e,a=r6_1(i1)%a,b=r6_1(i1)%b,cl=1.d0,nsu
* m=0
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
      r6_1(i1)%a(2)=(cauxa-ceps_0)
      r6_1(i1)%b(1)=(cauxb-ceps_2)
      end do
  
      cden= (-s561+cmw2)
* quqd -- p=p51,q=p6
      quqd=p51(0)*p6(0)-p51(1)*p6(1)-p51(2)*p6(2)-p51(3)*p6(3)
      ccl=wcl/cden
* TWR0 -- qu=p51,qd=p6,v=0,a=rw651(0)%a,b=rw651(0)%b,cl=ccl,nsum=0
      eps_0=-p51(2)*p6(3)+p6(2)*p51(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p51k0*p6(0)+p6k0*p51(0)
      rw651(0)%a(2)=ccl*(auxa-ceps_0)
      rw651(0)%b(1)=-ccl*(p6(2)+ceps_2)
* TWR0 -- qu=p51,qd=p6,v=1,a=rw651(1)%a,b=rw651(1)%b,cl=ccl,nsum=0
      auxa=-quqd+p51k0*p6(1)+p6k0*p51(1)
      rw651(1)%a(2)=ccl*(auxa-ceps_0)
      rw651(1)%b(1)=-ccl*(p6(2)+ceps_2)
* TWR0 -- qu=p51,qd=p6,v=2,a=rw651(2)%a,b=rw651(2)%b,cl=ccl,nsum=0
      eps_0=-p51k0*p6(3)+p6k0*p51(3)
      ceps_0=eps_0*cim
      auxa=p51k0*p6(2)+p6k0*p51(2)
      rw651(2)%a(2)=ccl*(auxa-ceps_0)
      rw651(2)%b(1)=-ccl*p6k0
* TWR0 -- qu=p51,qd=p6,v=3,a=rw651(3)%a,b=rw651(3)%b,cl=ccl,nsum=0
      eps_0=p51k0*p6(2)-p6k0*p51(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p51k0*p6(3)+p6k0*p51(3)
      rw651(3)%a(2)=ccl*(auxa-ceps_0)
      rw651(3)%b(1)=-ccl*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw516(i1)%e(m),a1=l5_1(i1)%a,c1=l5_1(i1)%c,a2=rw651(m)%a
* ,b2=rw651(m)%b,prq=s51,bef=,aft=
      cw516(i1)%e(m)=(l5_1(i1)%c(2)*s51*rw651(m)%b(1)+l5_1(i1)%a
     & (2)*rw651(m)%a(2))
      end do
      end do
  
      cden=(-s561+cmw2)*f61
* quqd -- p=p5,q=p61
      quqd=p5(0)*p61(0)-p5(1)*p61(1)-p5(2)*p61(2)-p5(3)*p61(3)
      ccl=wcl/cden
* TWL0 -- qu=p5,qd=p61,v=0,a=lw561(0)%a,c=lw561(0)%c,cl=ccl,nsum=0
      eps_0=-p5(2)*p61(3)+p61(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p61(0)+p61k0*p5(0)
      lw561(0)%a(2)=ccl*(auxa-ceps_0)
      lw561(0)%c(2)=ccl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p61,v=1,a=lw561(1)%a,c=lw561(1)%c,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p61(1)+p61k0*p5(1)
      lw561(1)%a(2)=ccl*(auxa-ceps_0)
      lw561(1)%c(2)=ccl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p61,v=2,a=lw561(2)%a,c=lw561(2)%c,cl=ccl,nsum=0
      eps_0=-p5k0*p61(3)+p61k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p61(2)+p61k0*p5(2)
      lw561(2)%a(2)=ccl*(auxa-ceps_0)
      lw561(2)%c(2)=-ccl*p5k0
* TWL0 -- qu=p5,qd=p61,v=3,a=lw561(3)%a,c=lw561(3)%c,cl=ccl,nsum=0
      eps_0=p5k0*p61(2)-p61k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p61(3)+p61k0*p5(3)
      lw561(3)%a(2)=ccl*(auxa-ceps_0)
      lw561(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw516(i1)%e(m),a1=lw561(m)%a,c1=lw561(m)%c,a2=r6_1(i1)%a
* ,b2=r6_1(i1)%b,prq=s61,bef=cw516(i1)%e(m)+,aft=
      cw516(i1)%e(m)=cw516(i1)%e(m)+(lw561(m)%c(2)*s61*r6_1(i1)%
     & b(1)+lw561(m)%a(2)*r6_1(i1)%a(2))
      end do
      end do
  
      do i1=1,2
* pk0 -- p=cw516(i1)%e
      cw516(i1)%ek0=cw516(i1)%e(0)-cw516(i1)%e(1)
      end do
  
      endif  ! id5 = quark
  
      do mu=0,3
         p52(mu)=p5(mu)+p2(mu)
      enddo
  
* for triple vertex                                                     
      do mu=0,3
       p526(mu)=p52(mu)+p6(mu)
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
      ccl=1.d0/(f52)
      do i1=1,2
* TWL0 -- qu=p5,qd=p52,v=ce2(i1)%e,a=l5_2(i1)%a,c=l5_2(i1)%c,cl=ccl,nsum
* =0
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
      l5_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s62/2d0
      do i1=1,2
* TWR0 -- qu=p62,qd=p6,v=ce2(i1)%e,a=r6_2(i1)%a,b=r6_2(i1)%b,cl=1.d0,nsu
* m=0
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
      r6_2(i1)%a(2)=(cauxa-ceps_0)
      r6_2(i1)%b(1)=(cauxb-ceps_2)
      end do
  
      cden= (-s562+cmw2)
* quqd -- p=p52,q=p6
      quqd=p52(0)*p6(0)-p52(1)*p6(1)-p52(2)*p6(2)-p52(3)*p6(3)
      ccl=wcl/cden
* TWR0 -- qu=p52,qd=p6,v=0,a=rw652(0)%a,b=rw652(0)%b,cl=ccl,nsum=0
      eps_0=-p52(2)*p6(3)+p6(2)*p52(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p52k0*p6(0)+p6k0*p52(0)
      rw652(0)%a(2)=ccl*(auxa-ceps_0)
      rw652(0)%b(1)=-ccl*(p6(2)+ceps_2)
* TWR0 -- qu=p52,qd=p6,v=1,a=rw652(1)%a,b=rw652(1)%b,cl=ccl,nsum=0
      auxa=-quqd+p52k0*p6(1)+p6k0*p52(1)
      rw652(1)%a(2)=ccl*(auxa-ceps_0)
      rw652(1)%b(1)=-ccl*(p6(2)+ceps_2)
* TWR0 -- qu=p52,qd=p6,v=2,a=rw652(2)%a,b=rw652(2)%b,cl=ccl,nsum=0
      eps_0=-p52k0*p6(3)+p6k0*p52(3)
      ceps_0=eps_0*cim
      auxa=p52k0*p6(2)+p6k0*p52(2)
      rw652(2)%a(2)=ccl*(auxa-ceps_0)
      rw652(2)%b(1)=-ccl*p6k0
* TWR0 -- qu=p52,qd=p6,v=3,a=rw652(3)%a,b=rw652(3)%b,cl=ccl,nsum=0
      eps_0=p52k0*p6(2)-p6k0*p52(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p52k0*p6(3)+p6k0*p52(3)
      rw652(3)%a(2)=ccl*(auxa-ceps_0)
      rw652(3)%b(1)=-ccl*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw526(i1)%e(m),a1=l5_2(i1)%a,c1=l5_2(i1)%c,a2=rw652(m)%a
* ,b2=rw652(m)%b,prq=s52,bef=,aft=
      cw526(i1)%e(m)=(l5_2(i1)%c(2)*s52*rw652(m)%b(1)+l5_2(i1)%a
     & (2)*rw652(m)%a(2))
      end do
      end do
  
      cden=(-s562+cmw2)*f62
* quqd -- p=p5,q=p62
      quqd=p5(0)*p62(0)-p5(1)*p62(1)-p5(2)*p62(2)-p5(3)*p62(3)
      ccl=wcl/cden
* TWL0 -- qu=p5,qd=p62,v=0,a=lw562(0)%a,c=lw562(0)%c,cl=ccl,nsum=0
      eps_0=-p5(2)*p62(3)+p62(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p62(0)+p62k0*p5(0)
      lw562(0)%a(2)=ccl*(auxa-ceps_0)
      lw562(0)%c(2)=ccl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p62,v=1,a=lw562(1)%a,c=lw562(1)%c,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p62(1)+p62k0*p5(1)
      lw562(1)%a(2)=ccl*(auxa-ceps_0)
      lw562(1)%c(2)=ccl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p62,v=2,a=lw562(2)%a,c=lw562(2)%c,cl=ccl,nsum=0
      eps_0=-p5k0*p62(3)+p62k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p62(2)+p62k0*p5(2)
      lw562(2)%a(2)=ccl*(auxa-ceps_0)
      lw562(2)%c(2)=-ccl*p5k0
* TWL0 -- qu=p5,qd=p62,v=3,a=lw562(3)%a,c=lw562(3)%c,cl=ccl,nsum=0
      eps_0=p5k0*p62(2)-p62k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p62(3)+p62k0*p5(3)
      lw562(3)%a(2)=ccl*(auxa-ceps_0)
      lw562(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw526(i1)%e(m),a1=lw562(m)%a,c1=lw562(m)%c,a2=r6_2(i1)%a
* ,b2=r6_2(i1)%b,prq=s62,bef=cw526(i1)%e(m)+,aft=
      cw526(i1)%e(m)=cw526(i1)%e(m)+(lw562(m)%c(2)*s62*r6_2(i1)%
     & b(1)+lw562(m)%a(2)*r6_2(i1)%a(2))
      end do
      end do
  
      do i1=1,2
* pk0 -- p=cw526(i1)%e
      cw526(i1)%ek0=cw526(i1)%e(0)-cw526(i1)%e(1)
      end do
  
      endif  ! id5 = quark
  
  
      do mu=0,3
         p71(mu)=p7(mu)+p1(mu)
      enddo
  
* for triple vertex                                                     
      do mu=0,3
       p718(mu)=p71(mu)+p8(mu)
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
      ccl=1.d0/(f71)
      do i1=1,2
* TWL0 -- qu=p7,qd=p71,v=ce1(i1)%e,a=l7_1(i1)%a,c=l7_1(i1)%c,cl=ccl,nsum
* =0
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
      l7_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s81/2d0
      do i1=1,2
* TWR0 -- qu=p81,qd=p8,v=ce1(i1)%e,a=r8_1(i1)%a,b=r8_1(i1)%b,cl=1.d0,nsu
* m=0
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
      r8_1(i1)%a(2)=(cauxa-ceps_0)
      r8_1(i1)%b(1)=(cauxb-ceps_2)
      end do
  
      cden= (-s781+cmw2)
* quqd -- p=p71,q=p8
      quqd=p71(0)*p8(0)-p71(1)*p8(1)-p71(2)*p8(2)-p71(3)*p8(3)
      ccl=wcl/cden
* TWR0 -- qu=p71,qd=p8,v=0,a=rw871(0)%a,b=rw871(0)%b,cl=ccl,nsum=0
      eps_0=-p71(2)*p8(3)+p8(2)*p71(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p71k0*p8(0)+p8k0*p71(0)
      rw871(0)%a(2)=ccl*(auxa-ceps_0)
      rw871(0)%b(1)=-ccl*(p8(2)+ceps_2)
* TWR0 -- qu=p71,qd=p8,v=1,a=rw871(1)%a,b=rw871(1)%b,cl=ccl,nsum=0
      auxa=-quqd+p71k0*p8(1)+p8k0*p71(1)
      rw871(1)%a(2)=ccl*(auxa-ceps_0)
      rw871(1)%b(1)=-ccl*(p8(2)+ceps_2)
* TWR0 -- qu=p71,qd=p8,v=2,a=rw871(2)%a,b=rw871(2)%b,cl=ccl,nsum=0
      eps_0=-p71k0*p8(3)+p8k0*p71(3)
      ceps_0=eps_0*cim
      auxa=p71k0*p8(2)+p8k0*p71(2)
      rw871(2)%a(2)=ccl*(auxa-ceps_0)
      rw871(2)%b(1)=-ccl*p8k0
* TWR0 -- qu=p71,qd=p8,v=3,a=rw871(3)%a,b=rw871(3)%b,cl=ccl,nsum=0
      eps_0=p71k0*p8(2)-p8k0*p71(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p71k0*p8(3)+p8k0*p71(3)
      rw871(3)%a(2)=ccl*(auxa-ceps_0)
      rw871(3)%b(1)=-ccl*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw718(i1)%e(m),a1=l7_1(i1)%a,c1=l7_1(i1)%c,a2=rw871(m)%a
* ,b2=rw871(m)%b,prq=s71,bef=,aft=
      cw718(i1)%e(m)=(l7_1(i1)%c(2)*s71*rw871(m)%b(1)+l7_1(i1)%a
     & (2)*rw871(m)%a(2))
      end do
      end do
  
      cden=(-s781+cmw2)*f81
* quqd -- p=p7,q=p81
      quqd=p7(0)*p81(0)-p7(1)*p81(1)-p7(2)*p81(2)-p7(3)*p81(3)
      ccl=wcl/cden
* TWL0 -- qu=p7,qd=p81,v=0,a=lw781(0)%a,c=lw781(0)%c,cl=ccl,nsum=0
      eps_0=-p7(2)*p81(3)+p81(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p81(0)+p81k0*p7(0)
      lw781(0)%a(2)=ccl*(auxa-ceps_0)
      lw781(0)%c(2)=ccl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p81,v=1,a=lw781(1)%a,c=lw781(1)%c,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p81(1)+p81k0*p7(1)
      lw781(1)%a(2)=ccl*(auxa-ceps_0)
      lw781(1)%c(2)=ccl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p81,v=2,a=lw781(2)%a,c=lw781(2)%c,cl=ccl,nsum=0
      eps_0=-p7k0*p81(3)+p81k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p81(2)+p81k0*p7(2)
      lw781(2)%a(2)=ccl*(auxa-ceps_0)
      lw781(2)%c(2)=-ccl*p7k0
* TWL0 -- qu=p7,qd=p81,v=3,a=lw781(3)%a,c=lw781(3)%c,cl=ccl,nsum=0
      eps_0=p7k0*p81(2)-p81k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p81(3)+p81k0*p7(3)
      lw781(3)%a(2)=ccl*(auxa-ceps_0)
      lw781(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw718(i1)%e(m),a1=lw781(m)%a,c1=lw781(m)%c,a2=r8_1(i1)%a
* ,b2=r8_1(i1)%b,prq=s81,bef=cw718(i1)%e(m)+,aft=
      cw718(i1)%e(m)=cw718(i1)%e(m)+(lw781(m)%c(2)*s81*r8_1(i1)%
     & b(1)+lw781(m)%a(2)*r8_1(i1)%a(2))
      end do
      end do
  
      do i1=1,2
* pk0 -- p=cw718(i1)%e
      cw718(i1)%ek0=cw718(i1)%e(0)-cw718(i1)%e(1)
      end do
  
      endif  ! id7 = quark
  
      do mu=0,3
         p72(mu)=p7(mu)+p2(mu)
      enddo
  
* for triple vertex                                                     
      do mu=0,3
       p728(mu)=p72(mu)+p8(mu)
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
      ccl=1.d0/(f72)
      do i1=1,2
* TWL0 -- qu=p7,qd=p72,v=ce2(i1)%e,a=l7_2(i1)%a,c=l7_2(i1)%c,cl=ccl,nsum
* =0
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
      l7_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      quqd=-s82/2d0
      do i1=1,2
* TWR0 -- qu=p82,qd=p8,v=ce2(i1)%e,a=r8_2(i1)%a,b=r8_2(i1)%b,cl=1.d0,nsu
* m=0
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
      r8_2(i1)%a(2)=(cauxa-ceps_0)
      r8_2(i1)%b(1)=(cauxb-ceps_2)
      end do
  
      cden= (-s782+cmw2)
* quqd -- p=p72,q=p8
      quqd=p72(0)*p8(0)-p72(1)*p8(1)-p72(2)*p8(2)-p72(3)*p8(3)
      ccl=wcl/cden
* TWR0 -- qu=p72,qd=p8,v=0,a=rw872(0)%a,b=rw872(0)%b,cl=ccl,nsum=0
      eps_0=-p72(2)*p8(3)+p8(2)*p72(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p72k0*p8(0)+p8k0*p72(0)
      rw872(0)%a(2)=ccl*(auxa-ceps_0)
      rw872(0)%b(1)=-ccl*(p8(2)+ceps_2)
* TWR0 -- qu=p72,qd=p8,v=1,a=rw872(1)%a,b=rw872(1)%b,cl=ccl,nsum=0
      auxa=-quqd+p72k0*p8(1)+p8k0*p72(1)
      rw872(1)%a(2)=ccl*(auxa-ceps_0)
      rw872(1)%b(1)=-ccl*(p8(2)+ceps_2)
* TWR0 -- qu=p72,qd=p8,v=2,a=rw872(2)%a,b=rw872(2)%b,cl=ccl,nsum=0
      eps_0=-p72k0*p8(3)+p8k0*p72(3)
      ceps_0=eps_0*cim
      auxa=p72k0*p8(2)+p8k0*p72(2)
      rw872(2)%a(2)=ccl*(auxa-ceps_0)
      rw872(2)%b(1)=-ccl*p8k0
* TWR0 -- qu=p72,qd=p8,v=3,a=rw872(3)%a,b=rw872(3)%b,cl=ccl,nsum=0
      eps_0=p72k0*p8(2)-p8k0*p72(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p72k0*p8(3)+p8k0*p72(3)
      rw872(3)%a(2)=ccl*(auxa-ceps_0)
      rw872(3)%b(1)=-ccl*ceps_2
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw728(i1)%e(m),a1=l7_2(i1)%a,c1=l7_2(i1)%c,a2=rw872(m)%a
* ,b2=rw872(m)%b,prq=s72,bef=,aft=
      cw728(i1)%e(m)=(l7_2(i1)%c(2)*s72*rw872(m)%b(1)+l7_2(i1)%a
     & (2)*rw872(m)%a(2))
      end do
      end do
  
      cden=(-s782+cmw2)*f82
* quqd -- p=p7,q=p82
      quqd=p7(0)*p82(0)-p7(1)*p82(1)-p7(2)*p82(2)-p7(3)*p82(3)
      ccl=wcl/cden
* TWL0 -- qu=p7,qd=p82,v=0,a=lw782(0)%a,c=lw782(0)%c,cl=ccl,nsum=0
      eps_0=-p7(2)*p82(3)+p82(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p82(0)+p82k0*p7(0)
      lw782(0)%a(2)=ccl*(auxa-ceps_0)
      lw782(0)%c(2)=ccl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p82,v=1,a=lw782(1)%a,c=lw782(1)%c,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p82(1)+p82k0*p7(1)
      lw782(1)%a(2)=ccl*(auxa-ceps_0)
      lw782(1)%c(2)=ccl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p82,v=2,a=lw782(2)%a,c=lw782(2)%c,cl=ccl,nsum=0
      eps_0=-p7k0*p82(3)+p82k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p82(2)+p82k0*p7(2)
      lw782(2)%a(2)=ccl*(auxa-ceps_0)
      lw782(2)%c(2)=-ccl*p7k0
* TWL0 -- qu=p7,qd=p82,v=3,a=lw782(3)%a,c=lw782(3)%c,cl=ccl,nsum=0
      eps_0=p7k0*p82(2)-p82k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p82(3)+p82k0*p7(3)
      lw782(3)%a(2)=ccl*(auxa-ceps_0)
      lw782(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw728(i1)%e(m),a1=lw782(m)%a,c1=lw782(m)%c,a2=r8_2(i1)%a
* ,b2=r8_2(i1)%b,prq=s82,bef=cw728(i1)%e(m)+,aft=
      cw728(i1)%e(m)=cw728(i1)%e(m)+(lw782(m)%c(2)*s82*r8_2(i1)%
     & b(1)+lw782(m)%a(2)*r8_2(i1)%a(2))
      end do
      end do
  
      do i1=1,2
* pk0 -- p=cw728(i1)%e
      cw728(i1)%ek0=cw728(i1)%e(0)-cw728(i1)%e(1)
      end do
  
      endif  ! id7 = quark
  
  
      do mu=0,3
         p31(mu)=p3(mu)+p1(mu)
      enddo
  
* for triple vertex                                                     
      do mu=0,3
       p314(mu)=p31(mu)+p4(mu)
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
      enddo
  
* for triple vertex                                                     
      do mu=0,3
       p324(mu)=p32(mu)+p4(mu)
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
  
*     -two-gluons fork  ----------------------------------------------  
  
      do m=0,3
       p512(m)=p51(m) + p2(m)
      enddo
  
* for triple vertex                                                     
      do mu=0,3
       p5126(mu)=p512(mu)+p6(mu)
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
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p612,qd=p6,v=ctrip12(i1,i2)%e,a=r6_gg(i1,i2)%a,b=r6_gg(i1,i
* 2)%b,cl=1.d0,nsum=0
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
      r6_gg(i1,i2)%a(2)=(cauxa-ceps_0)
      r6_gg(i1,i2)%b(1)=(cauxb-ceps_2)
      end do
      end do
  
* quqd -- p=p5,q=p512
      quqd=p5(0)*p512(0)-p5(1)*p512(1)-p5(2)*p512(2)-p5(3)*p512(
     & 3)
      ccl=1.d0/(f512)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p512,v=ctrip12(i1,i2)%e,a=l5_gg(i1,i2)%a,c=l5_gg(i1,i
* 2)%c,cl=ccl,nsum=0
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
      l5_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      cden=(-s5612+cmw2)
* quqd -- p=p512,q=p6
      quqd=p512(0)*p6(0)-p512(1)*p6(1)-p512(2)*p6(2)-p512(3)*p6(
     & 3)
      ccl=wcl/cden
* TWR0 -- qu=p512,qd=p6,v=0,a=rw6512(0)%a,b=rw6512(0)%b,cl=ccl,nsum=0
      eps_0=-p512(2)*p6(3)+p6(2)*p512(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p512k0*p6(0)+p6k0*p512(0)
      rw6512(0)%a(2)=ccl*(auxa-ceps_0)
      rw6512(0)%b(1)=-ccl*(p6(2)+ceps_2)
* TWR0 -- qu=p512,qd=p6,v=1,a=rw6512(1)%a,b=rw6512(1)%b,cl=ccl,nsum=0
      auxa=-quqd+p512k0*p6(1)+p6k0*p512(1)
      rw6512(1)%a(2)=ccl*(auxa-ceps_0)
      rw6512(1)%b(1)=-ccl*(p6(2)+ceps_2)
* TWR0 -- qu=p512,qd=p6,v=2,a=rw6512(2)%a,b=rw6512(2)%b,cl=ccl,nsum=0
      eps_0=-p512k0*p6(3)+p6k0*p512(3)
      ceps_0=eps_0*cim
      auxa=p512k0*p6(2)+p6k0*p512(2)
      rw6512(2)%a(2)=ccl*(auxa-ceps_0)
      rw6512(2)%b(1)=-ccl*p6k0
* TWR0 -- qu=p512,qd=p6,v=3,a=rw6512(3)%a,b=rw6512(3)%b,cl=ccl,nsum=0
      eps_0=p512k0*p6(2)-p6k0*p512(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p512k0*p6(3)+p6k0*p512(3)
      rw6512(3)%a(2)=ccl*(auxa-ceps_0)
      rw6512(3)%b(1)=-ccl*ceps_2
  
      cden=(-s5612+cmw2)*f612
* quqd -- p=p5,q=p612
      quqd=p5(0)*p612(0)-p5(1)*p612(1)-p5(2)*p612(2)-p5(3)*p612(
     & 3)
      ccl=wcl/cden
* TWL0 -- qu=p5,qd=p612,v=0,a=lw5612(0)%a,c=lw5612(0)%c,cl=ccl,nsum=0
      eps_0=-p5(2)*p612(3)+p612(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p612(0)+p612k0*p5(0)
      lw5612(0)%a(2)=ccl*(auxa-ceps_0)
      lw5612(0)%c(2)=ccl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p612,v=1,a=lw5612(1)%a,c=lw5612(1)%c,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p612(1)+p612k0*p5(1)
      lw5612(1)%a(2)=ccl*(auxa-ceps_0)
      lw5612(1)%c(2)=ccl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p612,v=2,a=lw5612(2)%a,c=lw5612(2)%c,cl=ccl,nsum=0
      eps_0=-p5k0*p612(3)+p612k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p612(2)+p612k0*p5(2)
      lw5612(2)%a(2)=ccl*(auxa-ceps_0)
      lw5612(2)%c(2)=-ccl*p5k0
* TWL0 -- qu=p5,qd=p612,v=3,a=lw5612(3)%a,c=lw5612(3)%c,cl=ccl,nsum=0
      eps_0=p5k0*p612(2)-p612k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p612(3)+p612k0*p5(3)
      lw5612(3)%a(2)=ccl*(auxa-ceps_0)
      lw5612(3)%c(2)=ccl*ceps_1
  
  
* quqd -- p=p51,q=p512
      quqd=p51(0)*p512(0)-p51(1)*p512(1)-p51(2)*p512(2)-p51(3)*p
     & 512(3)
      ccl=1.d0/(f512)
      do i2=1,2
* TW0 -- qu=p51,qd=p512,v=ce2(i2)%e,a=u51_2(i2)%a,b=u51_2(i2)%b,c=u51_2(
* i2)%c,d=u51_2(i2)%d,cl=ccl,nsum=0
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
      u51_2(i2)%a(2)=ccl*(cauxa-ceps_0)
      u51_2(i2)%b(1)=ccl*(cauxb-ceps_2)
      u51_2(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u51_2(i2)%d(1)=ccl*ce2(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_12(i1,i2)%a,cc=l5_12(i1,i2)%c,a1=l5_1(i1)%a,c1=l5_1(i1
* )%c,a2=u51_2(i2)%a,b2=u51_2(i2)%b,c2=u51_2(i2)%c,d2=u51_2(i2)%d,prq=s5
* 1,nsum=0
      l5_12(i1,i2)%c(2)=l5_1(i1)%c(2)*s51*u51_2(i2)%d(1)+l5_1(i1
     & )%a(2)*u51_2(i2)%c(2)
      l5_12(i1,i2)%a(2)=l5_1(i1)%c(2)*s51*u51_2(i2)%b(1)+l5_1(i1
     & )%a(2)*u51_2(i2)%a(2)
      end do
      end do
  
* quqd -- p=p612,q=p62
      quqd=p612(0)*p62(0)-p612(1)*p62(1)-p612(2)*p62(2)-p612(3)*
     & p62(3)
      ccl=1.d0/(f62)
      do i1=1,2
* TW0 -- qu=p612,qd=p62,v=ce1(i1)%e,a=u62_1(i1)%a,b=u62_1(i1)%b,c=u62_1(
* i1)%c,d=u62_1(i1)%d,cl=ccl,nsum=0
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
      u62_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u62_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u62_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u62_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0_W -- aa=r6_12(i1,i2)%a,bb=r6_12(i1,i2)%b,a1=u62_1(i1)%a,b1=u62_1(
* i1)%b,c1=u62_1(i1)%c,d1=u62_1(i1)%d,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s6
* 2,nsum=0
      r6_12(i1,i2)%b(1)=u62_1(i1)%d(1)*s62*r6_2(i2)%b(1)+u62_1(i
     & 1)%b(1)*r6_2(i2)%a(2)
      r6_12(i1,i2)%a(2)=u62_1(i1)%c(2)*s62*r6_2(i2)%b(1)+u62_1(i
     & 1)%a(2)*r6_2(i2)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l5_12(i1,i2)%a(2)=l5_12(i1,i2)%a(2) +
     &   l5_gg(i1,i2)%a(2)
       l5_12(i1,i2)%c(2)=l5_12(i1,i2)%c(2) +
     &   l5_gg(i1,i2)%c(2)
  
       r6_12(i1,i2)%a(2)=r6_12(i1,i2)%a(2) +
     &   r6_gg(i1,i2)%a(2)
       r6_12(i1,i2)%b(1)=r6_12(i1,i2)%b(1) +
     &   r6_gg(i1,i2)%b(1)
  
      enddo
      enddo
  
  
      cden=(-s5612+cmw2)*f62
* quqd -- p=p51,q=p62
      quqd=p51(0)*p62(0)-p51(1)*p62(1)-p51(2)*p62(2)-p51(3)*p62(
     & 3)
      ccl=wcl/cden
* TW0 -- qu=p51,qd=p62,v=0,a=uw51(0)%a,b=uw51(0)%b,c=uw51(0)%c,d=uw51(0)
* %d,cl=ccl,nsum=0
      eps_0=-p51(2)*p62(3)+p62(2)*p51(3)
      ceps_0=eps_0*cim
      ceps_1=p51(3)*cim
      ceps_2=p62(3)*cim
      auxa=-quqd+p51k0*p62(0)+p62k0*p51(0)
      uw51(0)%a(2)=ccl*(auxa-ceps_0)
      uw51(0)%b(1)=-ccl*(p62(2)+ceps_2)
      uw51(0)%c(2)=ccl*(-p51(2)+ceps_1)
      uw51(0)%d(1)=ccl
* TW0 -- qu=p51,qd=p62,v=1,a=uw51(1)%a,b=uw51(1)%b,c=uw51(1)%c,d=uw51(1)
* %d,cl=ccl,nsum=0
      auxa=-quqd+p51k0*p62(1)+p62k0*p51(1)
      uw51(1)%a(2)=ccl*(auxa-ceps_0)
      uw51(1)%b(1)=-ccl*(p62(2)+ceps_2)
      uw51(1)%c(2)=ccl*(-p51(2)+ceps_1)
      uw51(1)%d(1)=ccl
* TW0 -- qu=p51,qd=p62,v=2,a=uw51(2)%a,b=uw51(2)%b,c=uw51(2)%c,d=uw51(2)
* %d,cl=ccl,nsum=0
      eps_0=-p51k0*p62(3)+p62k0*p51(3)
      ceps_0=eps_0*cim
      auxa=p51k0*p62(2)+p62k0*p51(2)
      uw51(2)%a(2)=ccl*(auxa-ceps_0)
      uw51(2)%b(1)=-ccl*p62k0
      uw51(2)%c(2)=-ccl*p51k0
* TW0 -- qu=p51,qd=p62,v=3,a=uw51(3)%a,b=uw51(3)%b,c=uw51(3)%c,d=uw51(3)
* %d,cl=ccl,nsum=0
      eps_0=p51k0*p62(2)-p62k0*p51(2)
      ceps_0=eps_0*cim
      ceps_1=p51k0*cim
      ceps_2=p62k0*cim
      auxa=+p51k0*p62(3)+p62k0*p51(3)
      uw51(3)%a(2)=ccl*(auxa-ceps_0)
      uw51(3)%b(1)=-ccl*ceps_2
      uw51(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0_W -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l5_1(i1)%a,c1=l
* 5_1(i1)%c,a2=uw51(mu)%a,b2=uw51(mu)%b,c2=uw51(mu)%c,d2=uw51(mu)%d,prq=
* s51,nsum=0
      laux_imu(i1,mu)%c(2)=l5_1(i1)%c(2)*s51*uw51(mu)%d(1)+l5_1(
     & i1)%a(2)*uw51(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l5_1(i1)%c(2)*s51*uw51(mu)%b(1)+l5_1(
     & i1)%a(2)*uw51(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw5126(i1,i2)%e(mu),a1=laux_imu(i1,mu)%a,c1=laux_imu(i1,
* mu)%c,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s62,bef=,aft=
      cw5126(i1,i2)%e(mu)=(laux_imu(i1,mu)%c(2)*s62*r6_2(i2)%b(1
     & )+laux_imu(i1,mu)%a(2)*r6_2(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw5126(i1,i2)%e(mu),a1=l5_12(i1,i2)%a,c1=l5_12(i1,i2)%c,
* a2=rw6512(mu)%a,b2=rw6512(mu)%b,prq=s512,bef=cw5126(i1,i2)%e(mu)+,aft=
      cw5126(i1,i2)%e(mu)=cw5126(i1,i2)%e(mu)+(l5_12(i1,i2)%c(2)
     & *s512*rw6512(mu)%b(1)+l5_12(i1,i2)%a(2)*rw6512(mu)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw5126(i1,i2)%e(mu),a1=lw5612(mu)%a,c1=lw5612(mu)%c,a2=r
* 6_12(i1,i2)%a,b2=r6_12(i1,i2)%b,prq=s612,bef=cw5126(i1,i2)%e(mu)+,aft=
      cw5126(i1,i2)%e(mu)=cw5126(i1,i2)%e(mu)+(lw5612(mu)%c(2)*s
     & 612*r6_12(i1,i2)%b(1)+lw5612(mu)%a(2)*r6_12(i1,i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* pk0 -- p=cw5126(i1,i2)%e
      cw5126(i1,i2)%ek0=cw5126(i1,i2)%e(0)-cw5126(i1,i2)%e(1)
      end do
      end do
  
* quqd -- p=p52,q=p512
      quqd=p52(0)*p512(0)-p52(1)*p512(1)-p52(2)*p512(2)-p52(3)*p
     & 512(3)
      ccl=1.d0/(f512)
      do i2=1,2
* TW0 -- qu=p52,qd=p512,v=ce1(i2)%e,a=u52_1(i2)%a,b=u52_1(i2)%b,c=u52_1(
* i2)%c,d=u52_1(i2)%d,cl=ccl,nsum=0
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
      u52_1(i2)%a(2)=ccl*(cauxa-ceps_0)
      u52_1(i2)%b(1)=ccl*(cauxb-ceps_2)
      u52_1(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u52_1(i2)%d(1)=ccl*ce1(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_21(i1,i2)%a,cc=l5_21(i1,i2)%c,a1=l5_2(i2)%a,c1=l5_2(i2
* )%c,a2=u52_1(i1)%a,b2=u52_1(i1)%b,c2=u52_1(i1)%c,d2=u52_1(i1)%d,prq=s5
* 2,nsum=0
      l5_21(i1,i2)%c(2)=l5_2(i2)%c(2)*s52*u52_1(i1)%d(1)+l5_2(i2
     & )%a(2)*u52_1(i1)%c(2)
      l5_21(i1,i2)%a(2)=l5_2(i2)%c(2)*s52*u52_1(i1)%b(1)+l5_2(i2
     & )%a(2)*u52_1(i1)%a(2)
      end do
      end do
  
* quqd -- p=p612,q=p61
      quqd=p612(0)*p61(0)-p612(1)*p61(1)-p612(2)*p61(2)-p612(3)*
     & p61(3)
      ccl=1.d0/(f61)
      do i1=1,2
* TW0 -- qu=p612,qd=p61,v=ce2(i1)%e,a=u61_2(i1)%a,b=u61_2(i1)%b,c=u61_2(
* i1)%c,d=u61_2(i1)%d,cl=ccl,nsum=0
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
      u61_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u61_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u61_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u61_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0_W -- aa=r6_21(i1,i2)%a,bb=r6_21(i1,i2)%b,a1=u61_2(i2)%a,b1=u61_2(
* i2)%b,c1=u61_2(i2)%c,d1=u61_2(i2)%d,a2=r6_1(i1)%a,b2=r6_1(i1)%b,prq=s6
* 1,nsum=0
      r6_21(i1,i2)%b(1)=u61_2(i2)%d(1)*s61*r6_1(i1)%b(1)+u61_2(i
     & 2)%b(1)*r6_1(i1)%a(2)
      r6_21(i1,i2)%a(2)=u61_2(i2)%c(2)*s61*r6_1(i1)%b(1)+u61_2(i
     & 2)%a(2)*r6_1(i1)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l5_21(i1,i2)%a(2)=l5_21(i1,i2)%a(2) -
     &   l5_gg(i1,i2)%a(2)
       l5_21(i1,i2)%c(2)=l5_21(i1,i2)%c(2) -
     &   l5_gg(i1,i2)%c(2)
  
       r6_21(i1,i2)%a(2)=r6_21(i1,i2)%a(2) -
     &   r6_gg(i1,i2)%a(2)
       r6_21(i1,i2)%b(1)=r6_21(i1,i2)%b(1) -
     &   r6_gg(i1,i2)%b(1)
  
      enddo
      enddo
  
  
      cden=(-s5612+cmw2)*f61
* quqd -- p=p52,q=p61
      quqd=p52(0)*p61(0)-p52(1)*p61(1)-p52(2)*p61(2)-p52(3)*p61(
     & 3)
      ccl=wcl/cden
* TW0 -- qu=p52,qd=p61,v=0,a=uw52(0)%a,b=uw52(0)%b,c=uw52(0)%c,d=uw52(0)
* %d,cl=ccl,nsum=0
      eps_0=-p52(2)*p61(3)+p61(2)*p52(3)
      ceps_0=eps_0*cim
      ceps_1=p52(3)*cim
      ceps_2=p61(3)*cim
      auxa=-quqd+p52k0*p61(0)+p61k0*p52(0)
      uw52(0)%a(2)=ccl*(auxa-ceps_0)
      uw52(0)%b(1)=-ccl*(p61(2)+ceps_2)
      uw52(0)%c(2)=ccl*(-p52(2)+ceps_1)
      uw52(0)%d(1)=ccl
* TW0 -- qu=p52,qd=p61,v=1,a=uw52(1)%a,b=uw52(1)%b,c=uw52(1)%c,d=uw52(1)
* %d,cl=ccl,nsum=0
      auxa=-quqd+p52k0*p61(1)+p61k0*p52(1)
      uw52(1)%a(2)=ccl*(auxa-ceps_0)
      uw52(1)%b(1)=-ccl*(p61(2)+ceps_2)
      uw52(1)%c(2)=ccl*(-p52(2)+ceps_1)
      uw52(1)%d(1)=ccl
* TW0 -- qu=p52,qd=p61,v=2,a=uw52(2)%a,b=uw52(2)%b,c=uw52(2)%c,d=uw52(2)
* %d,cl=ccl,nsum=0
      eps_0=-p52k0*p61(3)+p61k0*p52(3)
      ceps_0=eps_0*cim
      auxa=p52k0*p61(2)+p61k0*p52(2)
      uw52(2)%a(2)=ccl*(auxa-ceps_0)
      uw52(2)%b(1)=-ccl*p61k0
      uw52(2)%c(2)=-ccl*p52k0
* TW0 -- qu=p52,qd=p61,v=3,a=uw52(3)%a,b=uw52(3)%b,c=uw52(3)%c,d=uw52(3)
* %d,cl=ccl,nsum=0
      eps_0=p52k0*p61(2)-p61k0*p52(2)
      ceps_0=eps_0*cim
      ceps_1=p52k0*cim
      ceps_2=p61k0*cim
      auxa=+p52k0*p61(3)+p61k0*p52(3)
      uw52(3)%a(2)=ccl*(auxa-ceps_0)
      uw52(3)%b(1)=-ccl*ceps_2
      uw52(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0_W -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l5_2(i1)%a,c1=l
* 5_2(i1)%c,a2=uw52(mu)%a,b2=uw52(mu)%b,c2=uw52(mu)%c,d2=uw52(mu)%d,prq=
* s52,nsum=0
      laux_imu(i1,mu)%c(2)=l5_2(i1)%c(2)*s52*uw52(mu)%d(1)+l5_2(
     & i1)%a(2)*uw52(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l5_2(i1)%c(2)*s52*uw52(mu)%b(1)+l5_2(
     & i1)%a(2)*uw52(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw5216(i1,i2)%e(mu),a1=laux_imu(i2,mu)%a,c1=laux_imu(i2,
* mu)%c,a2=r6_1(i1)%a,b2=r6_1(i1)%b,prq=s61,bef=,aft=
      cw5216(i1,i2)%e(mu)=(laux_imu(i2,mu)%c(2)*s61*r6_1(i1)%b(1
     & )+laux_imu(i2,mu)%a(2)*r6_1(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw5216(i1,i2)%e(mu),a1=l5_21(i1,i2)%a,c1=l5_21(i1,i2)%c,
* a2=rw6512(mu)%a,b2=rw6512(mu)%b,prq=s512,bef=cw5216(i1,i2)%e(mu)+,aft=
      cw5216(i1,i2)%e(mu)=cw5216(i1,i2)%e(mu)+(l5_21(i1,i2)%c(2)
     & *s512*rw6512(mu)%b(1)+l5_21(i1,i2)%a(2)*rw6512(mu)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw5216(i1,i2)%e(mu),a1=lw5612(mu)%a,c1=lw5612(mu)%c,a2=r
* 6_21(i1,i2)%a,b2=r6_21(i1,i2)%b,prq=s612,bef=cw5216(i1,i2)%e(mu)+,aft=
      cw5216(i1,i2)%e(mu)=cw5216(i1,i2)%e(mu)+(lw5612(mu)%c(2)*s
     & 612*r6_21(i1,i2)%b(1)+lw5612(mu)%a(2)*r6_21(i1,i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* pk0 -- p=cw5216(i1,i2)%e
      cw5216(i1,i2)%ek0=cw5216(i1,i2)%e(0)-cw5216(i1,i2)%e(1)
      end do
      end do
  
       endif !id5 = quark
  
      do m=0,3
       p712(m)=p71(m) + p2(m)
      enddo
  
* for triple vertex                                                     
      do mu=0,3
       p7128(mu)=p712(mu)+p8(mu)
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
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p812,qd=p8,v=ctrip12(i1,i2)%e,a=r8_gg(i1,i2)%a,b=r8_gg(i1,i
* 2)%b,cl=1.d0,nsum=0
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
      r8_gg(i1,i2)%a(2)=(cauxa-ceps_0)
      r8_gg(i1,i2)%b(1)=(cauxb-ceps_2)
      end do
      end do
  
* quqd -- p=p7,q=p712
      quqd=p7(0)*p712(0)-p7(1)*p712(1)-p7(2)*p712(2)-p7(3)*p712(
     & 3)
      ccl=1.d0/(f712)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p712,v=ctrip12(i1,i2)%e,a=l7_gg(i1,i2)%a,c=l7_gg(i1,i
* 2)%c,cl=ccl,nsum=0
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
      l7_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      cden=(-s7812+cmw2)
* quqd -- p=p712,q=p8
      quqd=p712(0)*p8(0)-p712(1)*p8(1)-p712(2)*p8(2)-p712(3)*p8(
     & 3)
      ccl=wcl/cden
* TWR0 -- qu=p712,qd=p8,v=0,a=rw8712(0)%a,b=rw8712(0)%b,cl=ccl,nsum=0
      eps_0=-p712(2)*p8(3)+p8(2)*p712(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p712k0*p8(0)+p8k0*p712(0)
      rw8712(0)%a(2)=ccl*(auxa-ceps_0)
      rw8712(0)%b(1)=-ccl*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=1,a=rw8712(1)%a,b=rw8712(1)%b,cl=ccl,nsum=0
      auxa=-quqd+p712k0*p8(1)+p8k0*p712(1)
      rw8712(1)%a(2)=ccl*(auxa-ceps_0)
      rw8712(1)%b(1)=-ccl*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=2,a=rw8712(2)%a,b=rw8712(2)%b,cl=ccl,nsum=0
      eps_0=-p712k0*p8(3)+p8k0*p712(3)
      ceps_0=eps_0*cim
      auxa=p712k0*p8(2)+p8k0*p712(2)
      rw8712(2)%a(2)=ccl*(auxa-ceps_0)
      rw8712(2)%b(1)=-ccl*p8k0
* TWR0 -- qu=p712,qd=p8,v=3,a=rw8712(3)%a,b=rw8712(3)%b,cl=ccl,nsum=0
      eps_0=p712k0*p8(2)-p8k0*p712(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p712k0*p8(3)+p8k0*p712(3)
      rw8712(3)%a(2)=ccl*(auxa-ceps_0)
      rw8712(3)%b(1)=-ccl*ceps_2
  
      cden=(-s7812+cmw2)*f812
* quqd -- p=p7,q=p812
      quqd=p7(0)*p812(0)-p7(1)*p812(1)-p7(2)*p812(2)-p7(3)*p812(
     & 3)
      ccl=wcl/cden
* TWL0 -- qu=p7,qd=p812,v=0,a=lw7812(0)%a,c=lw7812(0)%c,cl=ccl,nsum=0
      eps_0=-p7(2)*p812(3)+p812(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p812(0)+p812k0*p7(0)
      lw7812(0)%a(2)=ccl*(auxa-ceps_0)
      lw7812(0)%c(2)=ccl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=1,a=lw7812(1)%a,c=lw7812(1)%c,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p812(1)+p812k0*p7(1)
      lw7812(1)%a(2)=ccl*(auxa-ceps_0)
      lw7812(1)%c(2)=ccl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=2,a=lw7812(2)%a,c=lw7812(2)%c,cl=ccl,nsum=0
      eps_0=-p7k0*p812(3)+p812k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p812(2)+p812k0*p7(2)
      lw7812(2)%a(2)=ccl*(auxa-ceps_0)
      lw7812(2)%c(2)=-ccl*p7k0
* TWL0 -- qu=p7,qd=p812,v=3,a=lw7812(3)%a,c=lw7812(3)%c,cl=ccl,nsum=0
      eps_0=p7k0*p812(2)-p812k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p812(3)+p812k0*p7(3)
      lw7812(3)%a(2)=ccl*(auxa-ceps_0)
      lw7812(3)%c(2)=ccl*ceps_1
  
  
* quqd -- p=p71,q=p712
      quqd=p71(0)*p712(0)-p71(1)*p712(1)-p71(2)*p712(2)-p71(3)*p
     & 712(3)
      ccl=1.d0/(f712)
      do i2=1,2
* TW0 -- qu=p71,qd=p712,v=ce2(i2)%e,a=u71_2(i2)%a,b=u71_2(i2)%b,c=u71_2(
* i2)%c,d=u71_2(i2)%d,cl=ccl,nsum=0
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
      u71_2(i2)%a(2)=ccl*(cauxa-ceps_0)
      u71_2(i2)%b(1)=ccl*(cauxb-ceps_2)
      u71_2(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u71_2(i2)%d(1)=ccl*ce2(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_12(i1,i2)%a,cc=l7_12(i1,i2)%c,a1=l7_1(i1)%a,c1=l7_1(i1
* )%c,a2=u71_2(i2)%a,b2=u71_2(i2)%b,c2=u71_2(i2)%c,d2=u71_2(i2)%d,prq=s7
* 1,nsum=0
      l7_12(i1,i2)%c(2)=l7_1(i1)%c(2)*s71*u71_2(i2)%d(1)+l7_1(i1
     & )%a(2)*u71_2(i2)%c(2)
      l7_12(i1,i2)%a(2)=l7_1(i1)%c(2)*s71*u71_2(i2)%b(1)+l7_1(i1
     & )%a(2)*u71_2(i2)%a(2)
      end do
      end do
  
* quqd -- p=p812,q=p82
      quqd=p812(0)*p82(0)-p812(1)*p82(1)-p812(2)*p82(2)-p812(3)*
     & p82(3)
      ccl=1.d0/(f82)
      do i1=1,2
* TW0 -- qu=p812,qd=p82,v=ce1(i1)%e,a=u82_1(i1)%a,b=u82_1(i1)%b,c=u82_1(
* i1)%c,d=u82_1(i1)%d,cl=ccl,nsum=0
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
      u82_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u82_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u82_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u82_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0_W -- aa=r8_12(i1,i2)%a,bb=r8_12(i1,i2)%b,a1=u82_1(i1)%a,b1=u82_1(
* i1)%b,c1=u82_1(i1)%c,d1=u82_1(i1)%d,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s8
* 2,nsum=0
      r8_12(i1,i2)%b(1)=u82_1(i1)%d(1)*s82*r8_2(i2)%b(1)+u82_1(i
     & 1)%b(1)*r8_2(i2)%a(2)
      r8_12(i1,i2)%a(2)=u82_1(i1)%c(2)*s82*r8_2(i2)%b(1)+u82_1(i
     & 1)%a(2)*r8_2(i2)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l7_12(i1,i2)%a(2)=l7_12(i1,i2)%a(2) +
     &   l7_gg(i1,i2)%a(2)
       l7_12(i1,i2)%c(2)=l7_12(i1,i2)%c(2) +
     &   l7_gg(i1,i2)%c(2)
  
       r8_12(i1,i2)%a(2)=r8_12(i1,i2)%a(2) +
     &   r8_gg(i1,i2)%a(2)
       r8_12(i1,i2)%b(1)=r8_12(i1,i2)%b(1) +
     &   r8_gg(i1,i2)%b(1)
  
      enddo
      enddo
  
  
      cden=(-s7812+cmw2)*f82
* quqd -- p=p71,q=p82
      quqd=p71(0)*p82(0)-p71(1)*p82(1)-p71(2)*p82(2)-p71(3)*p82(
     & 3)
      ccl=wcl/cden
* TW0 -- qu=p71,qd=p82,v=0,a=uw71(0)%a,b=uw71(0)%b,c=uw71(0)%c,d=uw71(0)
* %d,cl=ccl,nsum=0
      eps_0=-p71(2)*p82(3)+p82(2)*p71(3)
      ceps_0=eps_0*cim
      ceps_1=p71(3)*cim
      ceps_2=p82(3)*cim
      auxa=-quqd+p71k0*p82(0)+p82k0*p71(0)
      uw71(0)%a(2)=ccl*(auxa-ceps_0)
      uw71(0)%b(1)=-ccl*(p82(2)+ceps_2)
      uw71(0)%c(2)=ccl*(-p71(2)+ceps_1)
      uw71(0)%d(1)=ccl
* TW0 -- qu=p71,qd=p82,v=1,a=uw71(1)%a,b=uw71(1)%b,c=uw71(1)%c,d=uw71(1)
* %d,cl=ccl,nsum=0
      auxa=-quqd+p71k0*p82(1)+p82k0*p71(1)
      uw71(1)%a(2)=ccl*(auxa-ceps_0)
      uw71(1)%b(1)=-ccl*(p82(2)+ceps_2)
      uw71(1)%c(2)=ccl*(-p71(2)+ceps_1)
      uw71(1)%d(1)=ccl
* TW0 -- qu=p71,qd=p82,v=2,a=uw71(2)%a,b=uw71(2)%b,c=uw71(2)%c,d=uw71(2)
* %d,cl=ccl,nsum=0
      eps_0=-p71k0*p82(3)+p82k0*p71(3)
      ceps_0=eps_0*cim
      auxa=p71k0*p82(2)+p82k0*p71(2)
      uw71(2)%a(2)=ccl*(auxa-ceps_0)
      uw71(2)%b(1)=-ccl*p82k0
      uw71(2)%c(2)=-ccl*p71k0
* TW0 -- qu=p71,qd=p82,v=3,a=uw71(3)%a,b=uw71(3)%b,c=uw71(3)%c,d=uw71(3)
* %d,cl=ccl,nsum=0
      eps_0=p71k0*p82(2)-p82k0*p71(2)
      ceps_0=eps_0*cim
      ceps_1=p71k0*cim
      ceps_2=p82k0*cim
      auxa=+p71k0*p82(3)+p82k0*p71(3)
      uw71(3)%a(2)=ccl*(auxa-ceps_0)
      uw71(3)%b(1)=-ccl*ceps_2
      uw71(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0_W -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l7_1(i1)%a,c1=l
* 7_1(i1)%c,a2=uw71(mu)%a,b2=uw71(mu)%b,c2=uw71(mu)%c,d2=uw71(mu)%d,prq=
* s71,nsum=0
      laux_imu(i1,mu)%c(2)=l7_1(i1)%c(2)*s71*uw71(mu)%d(1)+l7_1(
     & i1)%a(2)*uw71(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l7_1(i1)%c(2)*s71*uw71(mu)%b(1)+l7_1(
     & i1)%a(2)*uw71(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw7128(i1,i2)%e(mu),a1=laux_imu(i1,mu)%a,c1=laux_imu(i1,
* mu)%c,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s82,bef=,aft=
      cw7128(i1,i2)%e(mu)=(laux_imu(i1,mu)%c(2)*s82*r8_2(i2)%b(1
     & )+laux_imu(i1,mu)%a(2)*r8_2(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw7128(i1,i2)%e(mu),a1=l7_12(i1,i2)%a,c1=l7_12(i1,i2)%c,
* a2=rw8712(mu)%a,b2=rw8712(mu)%b,prq=s712,bef=cw7128(i1,i2)%e(mu)+,aft=
      cw7128(i1,i2)%e(mu)=cw7128(i1,i2)%e(mu)+(l7_12(i1,i2)%c(2)
     & *s712*rw8712(mu)%b(1)+l7_12(i1,i2)%a(2)*rw8712(mu)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw7128(i1,i2)%e(mu),a1=lw7812(mu)%a,c1=lw7812(mu)%c,a2=r
* 8_12(i1,i2)%a,b2=r8_12(i1,i2)%b,prq=s812,bef=cw7128(i1,i2)%e(mu)+,aft=
      cw7128(i1,i2)%e(mu)=cw7128(i1,i2)%e(mu)+(lw7812(mu)%c(2)*s
     & 812*r8_12(i1,i2)%b(1)+lw7812(mu)%a(2)*r8_12(i1,i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* pk0 -- p=cw7128(i1,i2)%e
      cw7128(i1,i2)%ek0=cw7128(i1,i2)%e(0)-cw7128(i1,i2)%e(1)
      end do
      end do
  
* quqd -- p=p72,q=p712
      quqd=p72(0)*p712(0)-p72(1)*p712(1)-p72(2)*p712(2)-p72(3)*p
     & 712(3)
      ccl=1.d0/(f712)
      do i2=1,2
* TW0 -- qu=p72,qd=p712,v=ce1(i2)%e,a=u72_1(i2)%a,b=u72_1(i2)%b,c=u72_1(
* i2)%c,d=u72_1(i2)%d,cl=ccl,nsum=0
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
      u72_1(i2)%a(2)=ccl*(cauxa-ceps_0)
      u72_1(i2)%b(1)=ccl*(cauxb-ceps_2)
      u72_1(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u72_1(i2)%d(1)=ccl*ce1(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_21(i1,i2)%a,cc=l7_21(i1,i2)%c,a1=l7_2(i2)%a,c1=l7_2(i2
* )%c,a2=u72_1(i1)%a,b2=u72_1(i1)%b,c2=u72_1(i1)%c,d2=u72_1(i1)%d,prq=s7
* 2,nsum=0
      l7_21(i1,i2)%c(2)=l7_2(i2)%c(2)*s72*u72_1(i1)%d(1)+l7_2(i2
     & )%a(2)*u72_1(i1)%c(2)
      l7_21(i1,i2)%a(2)=l7_2(i2)%c(2)*s72*u72_1(i1)%b(1)+l7_2(i2
     & )%a(2)*u72_1(i1)%a(2)
      end do
      end do
  
* quqd -- p=p812,q=p81
      quqd=p812(0)*p81(0)-p812(1)*p81(1)-p812(2)*p81(2)-p812(3)*
     & p81(3)
      ccl=1.d0/(f81)
      do i1=1,2
* TW0 -- qu=p812,qd=p81,v=ce2(i1)%e,a=u81_2(i1)%a,b=u81_2(i1)%b,c=u81_2(
* i1)%c,d=u81_2(i1)%d,cl=ccl,nsum=0
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
      u81_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u81_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u81_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u81_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TTR0_W -- aa=r8_21(i1,i2)%a,bb=r8_21(i1,i2)%b,a1=u81_2(i2)%a,b1=u81_2(
* i2)%b,c1=u81_2(i2)%c,d1=u81_2(i2)%d,a2=r8_1(i1)%a,b2=r8_1(i1)%b,prq=s8
* 1,nsum=0
      r8_21(i1,i2)%b(1)=u81_2(i2)%d(1)*s81*r8_1(i1)%b(1)+u81_2(i
     & 2)%b(1)*r8_1(i1)%a(2)
      r8_21(i1,i2)%a(2)=u81_2(i2)%c(2)*s81*r8_1(i1)%b(1)+u81_2(i
     & 2)%a(2)*r8_1(i1)%a(2)
      end do
      end do
  
*add triple vertex                                                      
      do i1=1,2
      do i2=1,2
  
       l7_21(i1,i2)%a(2)=l7_21(i1,i2)%a(2) -
     &   l7_gg(i1,i2)%a(2)
       l7_21(i1,i2)%c(2)=l7_21(i1,i2)%c(2) -
     &   l7_gg(i1,i2)%c(2)
  
       r8_21(i1,i2)%a(2)=r8_21(i1,i2)%a(2) -
     &   r8_gg(i1,i2)%a(2)
       r8_21(i1,i2)%b(1)=r8_21(i1,i2)%b(1) -
     &   r8_gg(i1,i2)%b(1)
  
      enddo
      enddo
  
  
      cden=(-s7812+cmw2)*f81
* quqd -- p=p72,q=p81
      quqd=p72(0)*p81(0)-p72(1)*p81(1)-p72(2)*p81(2)-p72(3)*p81(
     & 3)
      ccl=wcl/cden
* TW0 -- qu=p72,qd=p81,v=0,a=uw72(0)%a,b=uw72(0)%b,c=uw72(0)%c,d=uw72(0)
* %d,cl=ccl,nsum=0
      eps_0=-p72(2)*p81(3)+p81(2)*p72(3)
      ceps_0=eps_0*cim
      ceps_1=p72(3)*cim
      ceps_2=p81(3)*cim
      auxa=-quqd+p72k0*p81(0)+p81k0*p72(0)
      uw72(0)%a(2)=ccl*(auxa-ceps_0)
      uw72(0)%b(1)=-ccl*(p81(2)+ceps_2)
      uw72(0)%c(2)=ccl*(-p72(2)+ceps_1)
      uw72(0)%d(1)=ccl
* TW0 -- qu=p72,qd=p81,v=1,a=uw72(1)%a,b=uw72(1)%b,c=uw72(1)%c,d=uw72(1)
* %d,cl=ccl,nsum=0
      auxa=-quqd+p72k0*p81(1)+p81k0*p72(1)
      uw72(1)%a(2)=ccl*(auxa-ceps_0)
      uw72(1)%b(1)=-ccl*(p81(2)+ceps_2)
      uw72(1)%c(2)=ccl*(-p72(2)+ceps_1)
      uw72(1)%d(1)=ccl
* TW0 -- qu=p72,qd=p81,v=2,a=uw72(2)%a,b=uw72(2)%b,c=uw72(2)%c,d=uw72(2)
* %d,cl=ccl,nsum=0
      eps_0=-p72k0*p81(3)+p81k0*p72(3)
      ceps_0=eps_0*cim
      auxa=p72k0*p81(2)+p81k0*p72(2)
      uw72(2)%a(2)=ccl*(auxa-ceps_0)
      uw72(2)%b(1)=-ccl*p81k0
      uw72(2)%c(2)=-ccl*p72k0
* TW0 -- qu=p72,qd=p81,v=3,a=uw72(3)%a,b=uw72(3)%b,c=uw72(3)%c,d=uw72(3)
* %d,cl=ccl,nsum=0
      eps_0=p72k0*p81(2)-p81k0*p72(2)
      ceps_0=eps_0*cim
      ceps_1=p72k0*cim
      ceps_2=p81k0*cim
      auxa=+p72k0*p81(3)+p81k0*p72(3)
      uw72(3)%a(2)=ccl*(auxa-ceps_0)
      uw72(3)%b(1)=-ccl*ceps_2
      uw72(3)%c(2)=ccl*ceps_1
  
      do i1=1,2
      do mu=0,3
* TLT0_W -- aa=laux_imu(i1,mu)%a,cc=laux_imu(i1,mu)%c,a1=l7_2(i1)%a,c1=l
* 7_2(i1)%c,a2=uw72(mu)%a,b2=uw72(mu)%b,c2=uw72(mu)%c,d2=uw72(mu)%d,prq=
* s72,nsum=0
      laux_imu(i1,mu)%c(2)=l7_2(i1)%c(2)*s72*uw72(mu)%d(1)+l7_2(
     & i1)%a(2)*uw72(mu)%c(2)
      laux_imu(i1,mu)%a(2)=l7_2(i1)%c(2)*s72*uw72(mu)%b(1)+l7_2(
     & i1)%a(2)*uw72(mu)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw7218(i1,i2)%e(mu),a1=laux_imu(i2,mu)%a,c1=laux_imu(i2,
* mu)%c,a2=r8_1(i1)%a,b2=r8_1(i1)%b,prq=s81,bef=,aft=
      cw7218(i1,i2)%e(mu)=(laux_imu(i2,mu)%c(2)*s81*r8_1(i1)%b(1
     & )+laux_imu(i2,mu)%a(2)*r8_1(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw7218(i1,i2)%e(mu),a1=l7_21(i1,i2)%a,c1=l7_21(i1,i2)%c,
* a2=rw8712(mu)%a,b2=rw8712(mu)%b,prq=s712,bef=cw7218(i1,i2)%e(mu)+,aft=
      cw7218(i1,i2)%e(mu)=cw7218(i1,i2)%e(mu)+(l7_21(i1,i2)%c(2)
     & *s712*rw8712(mu)%b(1)+l7_21(i1,i2)%a(2)*rw8712(mu)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw7218(i1,i2)%e(mu),a1=lw7812(mu)%a,c1=lw7812(mu)%c,a2=r
* 8_21(i1,i2)%a,b2=r8_21(i1,i2)%b,prq=s812,bef=cw7218(i1,i2)%e(mu)+,aft=
      cw7218(i1,i2)%e(mu)=cw7218(i1,i2)%e(mu)+(lw7812(mu)%c(2)*s
     & 812*r8_21(i1,i2)%b(1)+lw7812(mu)%a(2)*r8_21(i1,i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* pk0 -- p=cw7218(i1,i2)%e
      cw7218(i1,i2)%ek0=cw7218(i1,i2)%e(0)-cw7218(i1,i2)%e(1)
      end do
      end do
  
       endif !id5 = quark
  
  
      do m=0,3
       p312(m)=p31(m) + p2(m)
      enddo
* for triple vertex                                                     
      do mu=0,3
       p3124(mu)=p312(mu)+p4(mu)
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
      do i1=1,2
      do i2=1,2
* TR0 -- qu=p412,qd=p4,v=ctrip12(i1,i2)%e,a=r4_gg(i1,i2)%a,b=r4_gg(i1,i2
* )%b,cr=1.d0,cl=1.d0,nsum=0
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
      r4_gg(i1,i2)%a(1)=(cauxa+ceps_0)
      r4_gg(i1,i2)%a(2)=(cauxa-ceps_0)
      r4_gg(i1,i2)%b(1)=(cauxb-ceps_2)
      r4_gg(i1,i2)%b(2)=(-cauxb-ceps_2)
      end do
      end do
  
* quqd -- p=p3,q=p312
      quqd=p3(0)*p312(0)-p3(1)*p312(1)-p3(2)*p312(2)-p3(3)*p312(
     & 3)
      ccr=1.d0/(f312)
      ccl=1.d0/(f312)
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
  
  
* -> lines without gluon ---------------------------------------------  
  
  
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
  
  
       do m=0,3
        p3156(m)=p356(m) +p1(m)
       enddo
* pk0 -- p=p3156
      p3156k0=p3156(0)-p3156(1)
* p.q -- p.q=s3156,p=p3156,q=p3156,bef=,aft=
      s3156=(p3156(0)*p3156(0)-p3156(1)*p3156(1)-p3156(2)*p3156(
     & 2)-p3156(3)*p3156(3))
      f3156=s3156*p3156k0
  
       do m=0,3
        p3256(m)=p356(m) +p2(m)
       enddo
* pk0 -- p=p3256
      p3256k0=p3256(0)-p3256(1)
* p.q -- p.q=s3256,p=p3256,q=p3256,bef=,aft=
      s3256=(p3256(0)*p3256(0)-p3256(1)*p3256(1)-p3256(2)*p3256(
     & 2)-p3256(3)*p3256(3))
      f3256=s3256*p3256k0
  
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
  
  
       do m=0,3
        p3178(m)=p378(m) +p1(m)
       enddo
* pk0 -- p=p3178
      p3178k0=p3178(0)-p3178(1)
* p.q -- p.q=s3178,p=p3178,q=p3178,bef=,aft=
      s3178=(p3178(0)*p3178(0)-p3178(1)*p3178(1)-p3178(2)*p3178(
     & 2)-p3178(3)*p3178(3))
      f3178=s3178*p3178k0
  
       do m=0,3
        p3278(m)=p378(m) +p2(m)
       enddo
* pk0 -- p=p3278
      p3278k0=p3278(0)-p3278(1)
* p.q -- p.q=s3278,p=p3278,q=p3278,bef=,aft=
      s3278=(p3278(0)*p3278(0)-p3278(1)*p3278(1)-p3278(2)*p3278(
     & 2)-p3278(3)*p3278(3))
      f3278=s3278*p3278k0
  
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
  
  
       do m=0,3
        p5134(m)=p534(m) +p1(m)
       enddo
* pk0 -- p=p5134
      p5134k0=p5134(0)-p5134(1)
* p.q -- p.q=s5134,p=p5134,q=p5134,bef=,aft=
      s5134=(p5134(0)*p5134(0)-p5134(1)*p5134(1)-p5134(2)*p5134(
     & 2)-p5134(3)*p5134(3))
      f5134=s5134*p5134k0
  
       do m=0,3
        p5234(m)=p534(m) +p2(m)
       enddo
* pk0 -- p=p5234
      p5234k0=p5234(0)-p5234(1)
* p.q -- p.q=s5234,p=p5234,q=p5234,bef=,aft=
      s5234=(p5234(0)*p5234(0)-p5234(1)*p5234(1)-p5234(2)*p5234(
     & 2)-p5234(3)*p5234(3))
      f5234=s5234*p5234k0
  
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
  
  
       do m=0,3
        p5178(m)=p578(m) +p1(m)
       enddo
* pk0 -- p=p5178
      p5178k0=p5178(0)-p5178(1)
* p.q -- p.q=s5178,p=p5178,q=p5178,bef=,aft=
      s5178=(p5178(0)*p5178(0)-p5178(1)*p5178(1)-p5178(2)*p5178(
     & 2)-p5178(3)*p5178(3))
      f5178=s5178*p5178k0
  
       do m=0,3
        p5278(m)=p578(m) +p2(m)
       enddo
* pk0 -- p=p5278
      p5278k0=p5278(0)-p5278(1)
* p.q -- p.q=s5278,p=p5278,q=p5278,bef=,aft=
      s5278=(p5278(0)*p5278(0)-p5278(1)*p5278(1)-p5278(2)*p5278(
     & 2)-p5278(3)*p5278(3))
      f5278=s5278*p5278k0
  
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
  
  
       do m=0,3
        p7134(m)=p734(m) +p1(m)
       enddo
* pk0 -- p=p7134
      p7134k0=p7134(0)-p7134(1)
* p.q -- p.q=s7134,p=p7134,q=p7134,bef=,aft=
      s7134=(p7134(0)*p7134(0)-p7134(1)*p7134(1)-p7134(2)*p7134(
     & 2)-p7134(3)*p7134(3))
      f7134=s7134*p7134k0
  
       do m=0,3
        p7234(m)=p734(m) +p2(m)
       enddo
* pk0 -- p=p7234
      p7234k0=p7234(0)-p7234(1)
* p.q -- p.q=s7234,p=p7234,q=p7234,bef=,aft=
      s7234=(p7234(0)*p7234(0)-p7234(1)*p7234(1)-p7234(2)*p7234(
     & 2)-p7234(3)*p7234(3))
      f7234=s7234*p7234k0
  
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
  
  
       do m=0,3
        p7156(m)=p756(m) +p1(m)
       enddo
* pk0 -- p=p7156
      p7156k0=p7156(0)-p7156(1)
* p.q -- p.q=s7156,p=p7156,q=p7156,bef=,aft=
      s7156=(p7156(0)*p7156(0)-p7156(1)*p7156(1)-p7156(2)*p7156(
     & 2)-p7156(3)*p7156(3))
      f7156=s7156*p7156k0
  
       do m=0,3
        p7256(m)=p756(m) +p2(m)
       enddo
* pk0 -- p=p7256
      p7256k0=p7256(0)-p7256(1)
* p.q -- p.q=s7256,p=p7256,q=p7256,bef=,aft=
      s7256=(p7256(0)*p7256(0)-p7256(1)*p7256(1)-p7256(2)*p7256(
     & 2)-p7256(3)*p7256(3))
      f7256=s7256*p7256k0
  
* W                                                                     
  
      if (ilept(id5).ne.1.or.ilept(id7).ne.1) then
      quqd=(s534-s34)/2d0
      ccl=zcl(id5)/(f534)
      do i3=1,2
* TWL0 -- qu=p5,qd=p534,v=cz34(i3)%e,a=l5_34(i3)%a,c=l5_34(i3)%c,cl=ccl,
* nsum=0
      ceps_0=-cz34(i3)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*(c
     & z34(i3)%e(2)*p534(3)-p534(2)*cz34(i3)%e(3))-p534k0*(cz34(
     & i3)%e(2)*p5(3)-p5(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p5k0+p5(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i3)%e(0)*p5(0)-cz34(i3)%e(1)*p5(1)-cz34(i3)%e(2)
     & *p5(2)-cz34(i3)%e(3)*p5(3)
      cvqd=cz34(i3)%e(0)*p534(0)-cz34(i3)%e(1)*p534(1)-cz34(i3)%
     & e(2)*p534(2)-cz34(i3)%e(3)*p534(3)
      cauxa=-cz34(i3)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cz34(i3)%ek0*p5(2)-p5k0*cz34(i3)%e(2)
      l5_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      l5_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
        if (ineutri(id5).ne.1) then
      ccl=fcl(id5)/(f534)
      do i3=1,2
* TWL0 -- qu=p5,qd=p534,v=cf34(i3)%e,a=l5_34(i3)%a,c=l5_34(i3)%c,cl=ccl,
* nsum=1
      ceps_0=-cf34(i3)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*(c
     & f34(i3)%e(2)*p534(3)-p534(2)*cf34(i3)%e(3))-p534k0*(cf34(
     & i3)%e(2)*p5(3)-p5(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p5k0+p5(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf34(i3)%e(0)*p5(0)-cf34(i3)%e(1)*p5(1)-cf34(i3)%e(2)
     & *p5(2)-cf34(i3)%e(3)*p5(3)
      cvqd=cf34(i3)%e(0)*p534(0)-cf34(i3)%e(1)*p534(1)-cf34(i3)%
     & e(2)*p534(2)-cf34(i3)%e(3)*p534(3)
      cauxa=-cf34(i3)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cf34(i3)%ek0*p5(2)-p5k0*cf34(i3)%e(2)
      l5_34(i3)%a(2)=l5_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      l5_34(i3)%c(2)=l5_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
        endif
  
      endif
  
      if (ilept(id5).ne.1.or.ilept(id3).ne.1) then
      quqd=(s578-s78)/2d0
      ccl=wcl/(f578)
* TWL0 -- qu=p5,qd=p578,v=cw78%e,a=l5_78%a,c=l5_78%c,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p5(2)*p578(3)-p578(2)*p5(3))+p5k0*(cw78%
     & e(2)*p578(3)-p578(2)*cw78%e(3))-p578k0*(cw78%e(2)*p5(3)-p
     & 5(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p5k0+p5(3)*cw78%ek0
      ceps_1=ceps_1*cim
      cvqu=cw78%e(0)*p5(0)-cw78%e(1)*p5(1)-cw78%e(2)*p5(2)-cw78%
     & e(3)*p5(3)
      cvqd=cw78%e(0)*p578(0)-cw78%e(1)*p578(1)-cw78%e(2)*p578(2)
     & -cw78%e(3)*p578(3)
      cauxa=-cw78%ek0*quqd+p5k0*cvqd+p578k0*cvqu
      cauxc=+cw78%ek0*p5(2)-p5k0*cw78%e(2)
      l5_78%a(2)=ccl*(cauxa-ceps_0)
      l5_78%c(2)=ccl*(-cauxc+ceps_1)
      endif
  
      if (ilept(id6).ne.1.or.ilept(id3).ne.1) then
      quqd=(-s678+s78)/2d0
* TWR0 -- qu=p678,qd=p6,v=cw78%e,a=r6_78%a,b=r6_78%b,cl=wcl,nsum=0
      ceps_0=-cw78%ek0*(p678(2)*p6(3)-p6(2)*p678(3))+p678k0*(cw7
     & 8%e(2)*p6(3)-p6(2)*cw78%e(3))-p6k0*(cw78%e(2)*p678(3)-p67
     & 8(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw78%e(3)*p6k0+p6(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p678(0)-cw78%e(1)*p678(1)-cw78%e(2)*p678(2)
     & -cw78%e(3)*p678(3)
      cvqd=cw78%e(0)*p6(0)-cw78%e(1)*p6(1)-cw78%e(2)*p6(2)-cw78%
     & e(3)*p6(3)
      cauxa=-cw78%ek0*quqd+p678k0*cvqd+p6k0*cvqu
      cauxb=-cw78%ek0*p6(2)+p6k0*cw78%e(2)
      r6_78%a(2)=wcl*(cauxa-ceps_0)
      r6_78%b(1)=wcl*(cauxb-ceps_2)
      endif
  
      if (ilept(id6).ne.1.or.ilept(id7).ne.1) then
      quqd=(-s634+s34)/2d0
      do i3=1,2
* TWR0 -- qu=p634,qd=p6,v=cz34(i3)%e,a=r6_34(i3)%a,b=r6_34(i3)%b,cl=zcl(
* id6),nsum=0
      ceps_0=-cz34(i3)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k0*
     & (cz34(i3)%e(2)*p6(3)-p6(2)*cz34(i3)%e(3))-p6k0*(cz34(i3)%
     & e(2)*p634(3)-p634(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i3)%e(3)*p6k0+p6(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p634(0)-cz34(i3)%e(1)*p634(1)-cz34(i3)%
     & e(2)*p634(2)-cz34(i3)%e(3)*p634(3)
      cvqd=cz34(i3)%e(0)*p6(0)-cz34(i3)%e(1)*p6(1)-cz34(i3)%e(2)
     & *p6(2)-cz34(i3)%e(3)*p6(3)
      cauxa=-cz34(i3)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cz34(i3)%ek0*p6(2)+p6k0*cz34(i3)%e(2)
      r6_34(i3)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_34(i3)%b(1)=zcl(id6)*(cauxb-ceps_2)
      end do
  
        if (ineutri(id6).ne.1) then
      do i3=1,2
* TWR0 -- qu=p634,qd=p6,v=cf34(i3)%e,a=r6_34(i3)%a,b=r6_34(i3)%b,cl=fcl(
* id6),nsum=1
      ceps_0=-cf34(i3)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k0*
     & (cf34(i3)%e(2)*p6(3)-p6(2)*cf34(i3)%e(3))-p6k0*(cf34(i3)%
     & e(2)*p634(3)-p634(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf34(i3)%e(3)*p6k0+p6(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p634(0)-cf34(i3)%e(1)*p634(1)-cf34(i3)%
     & e(2)*p634(2)-cf34(i3)%e(3)*p634(3)
      cvqd=cf34(i3)%e(0)*p6(0)-cf34(i3)%e(1)*p6(1)-cf34(i3)%e(2)
     & *p6(2)-cf34(i3)%e(3)*p6(3)
      cauxa=-cf34(i3)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cf34(i3)%ek0*p6(2)+p6k0*cf34(i3)%e(2)
      r6_34(i3)%a(2)=r6_34(i3)%a(2)+fcl(id6)*(cauxa-ceps_0)
      r6_34(i3)%b(1)=r6_34(i3)%b(1)+fcl(id6)*(cauxb-ceps_2)
      end do
        endif
      endif
  
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p5,q=p678
      quqd=p5(0)*p678(0)-p5(1)*p678(1)-p5(2)*p678(2)-p5(3)*p678(
     & 3)
      ccl=zcl(id5)/(f678)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p678,v=cz3124(i3,i1,i2)%e,a=l5_3124(i3,i1,i2)%a,c=l5_
* 3124(i3,i1,i2)%c,cl=ccl,nsum=0
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
      l5_3124(i3,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_3124(i3,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
        if (ineutri(id5).ne.1) then
      ccl=fcl(id5)/(f678)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p678,v=cf3124(i3,i1,i2)%e,a=l5_3124(i3,i1,i2)%a,c=l5_
* 3124(i3,i1,i2)%c,cl=ccl,nsum=1
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
      l5_3124(i3,i1,i2)%a(2)=l5_3124(i3,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l5_3124(i3,i1,i2)%c(2)=l5_3124(i3,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
        endif
  
      if (ilept(id5).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p5,q=p5134
      quqd=p5(0)*p5134(0)-p5(1)*p5134(1)-p5(2)*p5134(2)-p5(3)*p5
     & 134(3)
      ccl=zcl(id5)/(f5134)
      do i3=1,2
      do i1=1,2
* TWL0 -- qu=p5,qd=p5134,v=cz314(i3,i1)%e,a=l5_314(i3,i1)%a,c=l5_314(i3,
* i1)%c,cl=ccl,nsum=0
      ceps_0=-cz314(i3,i1)%ek0*(p5(2)*p5134(3)-p5134(2)*p5(3))+p
     & 5k0*(cz314(i3,i1)%e(2)*p5134(3)-p5134(2)*cz314(i3,i1)%e(3
     & ))-p5134k0*(cz314(i3,i1)%e(2)*p5(3)-p5(2)*cz314(i3,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i3,i1)%e(3)*p5k0+p5(3)*cz314(i3,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz314(i3,i1)%e(0)*p5(0)-cz314(i3,i1)%e(1)*p5(1)-cz314
     & (i3,i1)%e(2)*p5(2)-cz314(i3,i1)%e(3)*p5(3)
      cvqd=cz314(i3,i1)%e(0)*p5134(0)-cz314(i3,i1)%e(1)*p5134(1)
     & -cz314(i3,i1)%e(2)*p5134(2)-cz314(i3,i1)%e(3)*p5134(3)
      cauxa=-cz314(i3,i1)%ek0*quqd+p5k0*cvqd+p5134k0*cvqu
      cauxc=+cz314(i3,i1)%ek0*p5(2)-p5k0*cz314(i3,i1)%e(2)
      l5_314(i3,i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_314(i3,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
       if (ineutri(id5).ne.1) then
      ccl=fcl(id5)/(f5134)
      do i3=1,2
      do i1=1,2
* TWL0 -- qu=p5,qd=p5134,v=cf314(i3,i1)%e,a=l5_314(i3,i1)%a,c=l5_314(i3,
* i1)%c,cl=ccl,nsum=1
      ceps_0=-cf314(i3,i1)%ek0*(p5(2)*p5134(3)-p5134(2)*p5(3))+p
     & 5k0*(cf314(i3,i1)%e(2)*p5134(3)-p5134(2)*cf314(i3,i1)%e(3
     & ))-p5134k0*(cf314(i3,i1)%e(2)*p5(3)-p5(2)*cf314(i3,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i3,i1)%e(3)*p5k0+p5(3)*cf314(i3,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf314(i3,i1)%e(0)*p5(0)-cf314(i3,i1)%e(1)*p5(1)-cf314
     & (i3,i1)%e(2)*p5(2)-cf314(i3,i1)%e(3)*p5(3)
      cvqd=cf314(i3,i1)%e(0)*p5134(0)-cf314(i3,i1)%e(1)*p5134(1)
     & -cf314(i3,i1)%e(2)*p5134(2)-cf314(i3,i1)%e(3)*p5134(3)
      cauxa=-cf314(i3,i1)%ek0*quqd+p5k0*cvqd+p5134k0*cvqu
      cauxc=+cf314(i3,i1)%ek0*p5(2)-p5k0*cf314(i3,i1)%e(2)
      l5_314(i3,i1)%a(2)=l5_314(i3,i1)%a(2)+ccl*(cauxa-ceps_0)
      l5_314(i3,i1)%c(2)=l5_314(i3,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
        endif
  
      endif
      endif
  
      if (ilept(id7).ne.1) then
* quqd -- p=p5,q=p634
      quqd=p5(0)*p634(0)-p5(1)*p634(1)-p5(2)*p634(2)-p5(3)*p634(
     & 3)
      ccl=wcl/(f634)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p634,v=cw7128(i1,i2)%e,a=l5_7128(i1,i2)%a,c=l5_7128(i
* 1,i2)%c,cl=ccl,nsum=0
      ceps_0=-cw7128(i1,i2)%ek0*(p5(2)*p634(3)-p634(2)*p5(3))+p5
     & k0*(cw7128(i1,i2)%e(2)*p634(3)-p634(2)*cw7128(i1,i2)%e(3)
     & )-p634k0*(cw7128(i1,i2)%e(2)*p5(3)-p5(2)*cw7128(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cw7128(i1,i2)%e(3)*p5k0+p5(3)*cw7128(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw7128(i1,i2)%e(0)*p5(0)-cw7128(i1,i2)%e(1)*p5(1)-cw7
     & 128(i1,i2)%e(2)*p5(2)-cw7128(i1,i2)%e(3)*p5(3)
      cvqd=cw7128(i1,i2)%e(0)*p634(0)-cw7128(i1,i2)%e(1)*p634(1)
     & -cw7128(i1,i2)%e(2)*p634(2)-cw7128(i1,i2)%e(3)*p634(3)
      cauxa=-cw7128(i1,i2)%ek0*quqd+p5k0*cvqd+p634k0*cvqu
      cauxc=+cw7128(i1,i2)%ek0*p5(2)-p5k0*cw7128(i1,i2)%e(2)
      l5_7128(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_7128(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      if (ilept(id5).ne.1.or.ilept(id3).ne.1) then
* quqd -- p=p5,q=p5178
      quqd=p5(0)*p5178(0)-p5(1)*p5178(1)-p5(2)*p5178(2)-p5(3)*p5
     & 178(3)
      ccl=wcl/(f5178)
      do i1=1,2
* TWL0 -- qu=p5,qd=p5178,v=cw718(i1)%e,a=l5_718(i1)%a,c=l5_718(i1)%c,cl=
* ccl,nsum=0
      ceps_0=-cw718(i1)%ek0*(p5(2)*p5178(3)-p5178(2)*p5(3))+p5k0
     & *(cw718(i1)%e(2)*p5178(3)-p5178(2)*cw718(i1)%e(3))-p5178k
     & 0*(cw718(i1)%e(2)*p5(3)-p5(2)*cw718(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw718(i1)%e(3)*p5k0+p5(3)*cw718(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw718(i1)%e(0)*p5(0)-cw718(i1)%e(1)*p5(1)-cw718(i1)%e
     & (2)*p5(2)-cw718(i1)%e(3)*p5(3)
      cvqd=cw718(i1)%e(0)*p5178(0)-cw718(i1)%e(1)*p5178(1)-cw718
     & (i1)%e(2)*p5178(2)-cw718(i1)%e(3)*p5178(3)
      cauxa=-cw718(i1)%ek0*quqd+p5k0*cvqd+p5178k0*cvqu
      cauxc=+cw718(i1)%ek0*p5(2)-p5k0*cw718(i1)%e(2)
      l5_718(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_718(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* quqd -- p=p534,q=p6
      quqd=p534(0)*p6(0)-p534(1)*p6(1)-p534(2)*p6(2)-p534(3)*p6(
     & 3)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p534,qd=p6,v=cw7128(i1,i2)%e,a=r6_7128(i1,i2)%a,b=r6_7128(i
* 1,i2)%b,cl=wcl,nsum=0
      ceps_0=-cw7128(i1,i2)%ek0*(p534(2)*p6(3)-p6(2)*p534(3))+p5
     & 34k0*(cw7128(i1,i2)%e(2)*p6(3)-p6(2)*cw7128(i1,i2)%e(3))-
     & p6k0*(cw7128(i1,i2)%e(2)*p534(3)-p534(2)*cw7128(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cw7128(i1,i2)%e(3)*p6k0+p6(3)*cw7128(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw7128(i1,i2)%e(0)*p534(0)-cw7128(i1,i2)%e(1)*p534(1)
     & -cw7128(i1,i2)%e(2)*p534(2)-cw7128(i1,i2)%e(3)*p534(3)
      cvqd=cw7128(i1,i2)%e(0)*p6(0)-cw7128(i1,i2)%e(1)*p6(1)-cw7
     & 128(i1,i2)%e(2)*p6(2)-cw7128(i1,i2)%e(3)*p6(3)
      cauxa=-cw7128(i1,i2)%ek0*quqd+p534k0*cvqd+p6k0*cvqu
      cauxb=-cw7128(i1,i2)%ek0*p6(2)+p6k0*cw7128(i1,i2)%e(2)
      r6_7128(i1,i2)%a(2)=wcl*(cauxa-ceps_0)
      r6_7128(i1,i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
      end do
  
      if (ilept(id6).ne.1.or.ilept(id3).ne.1) then
* quqd -- p=p5134,q=p6
      quqd=p5134(0)*p6(0)-p5134(1)*p6(1)-p5134(2)*p6(2)-p5134(3)
     & *p6(3)
      do i2=1,2
* TWR0 -- qu=p5134,qd=p6,v=cw728(i2)%e,a=r6_728(i2)%a,b=r6_728(i2)%b,cl=
* wcl,nsum=0
      ceps_0=-cw728(i2)%ek0*(p5134(2)*p6(3)-p6(2)*p5134(3))+p513
     & 4k0*(cw728(i2)%e(2)*p6(3)-p6(2)*cw728(i2)%e(3))-p6k0*(cw7
     & 28(i2)%e(2)*p5134(3)-p5134(2)*cw728(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw728(i2)%e(3)*p6k0+p6(3)*cw728(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw728(i2)%e(0)*p5134(0)-cw728(i2)%e(1)*p5134(1)-cw728
     & (i2)%e(2)*p5134(2)-cw728(i2)%e(3)*p5134(3)
      cvqd=cw728(i2)%e(0)*p6(0)-cw728(i2)%e(1)*p6(1)-cw728(i2)%e
     & (2)*p6(2)-cw728(i2)%e(3)*p6(3)
      cauxa=-cw728(i2)%ek0*quqd+p5134k0*cvqd+p6k0*cvqu
      cauxb=-cw728(i2)%ek0*p6(2)+p6k0*cw728(i2)%e(2)
      r6_728(i2)%a(2)=wcl*(cauxa-ceps_0)
      r6_728(i2)%b(1)=wcl*(cauxb-ceps_2)
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
* TWR0 -- qu=p578,qd=p6,v=cz3124(i3,i1,i2)%e,a=r6_3124(i3,i1,i2)%a,b=r6_
* 3124(i3,i1,i2)%b,cl=zcl(id6),nsum=0
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
      r6_3124(i3,i1,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_3124(i3,i1,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      end do
      end do
      end do
  
       if (ineutri(id6).ne.1) then
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p578,qd=p6,v=cf3124(i3,i1,i2)%e,a=r6_3124(i3,i1,i2)%a,b=r6_
* 3124(i3,i1,i2)%b,cl=fcl(id6),nsum=1
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
      r6_3124(i3,i1,i2)%a(2)=r6_3124(i3,i1,i2)%a(2)+fcl(id6)*(ca
     & uxa-ceps_0)
      r6_3124(i3,i1,i2)%b(1)=r6_3124(i3,i1,i2)%b(1)+fcl(id6)*(ca
     & uxb-ceps_2)
      end do
      end do
      end do
       endif
  
      if (ilept(id6).ne.1.or.ilept(id7).ne.1) then
* quqd -- p=p5178,q=p6
      quqd=p5178(0)*p6(0)-p5178(1)*p6(1)-p5178(2)*p6(2)-p5178(3)
     & *p6(3)
      do i3=1,2
      do i2=1,2
* TWR0 -- qu=p5178,qd=p6,v=cz324(i3,i2)%e,a=r6_324(i3,i2)%a,b=r6_324(i3,
* i2)%b,cl=zcl(id6),nsum=0
      ceps_0=-cz324(i3,i2)%ek0*(p5178(2)*p6(3)-p6(2)*p5178(3))+p
     & 5178k0*(cz324(i3,i2)%e(2)*p6(3)-p6(2)*cz324(i3,i2)%e(3))-
     & p6k0*(cz324(i3,i2)%e(2)*p5178(3)-p5178(2)*cz324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz324(i3,i2)%e(3)*p6k0+p6(3)*cz324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i3,i2)%e(0)*p5178(0)-cz324(i3,i2)%e(1)*p5178(1)
     & -cz324(i3,i2)%e(2)*p5178(2)-cz324(i3,i2)%e(3)*p5178(3)
      cvqd=cz324(i3,i2)%e(0)*p6(0)-cz324(i3,i2)%e(1)*p6(1)-cz324
     & (i3,i2)%e(2)*p6(2)-cz324(i3,i2)%e(3)*p6(3)
      cauxa=-cz324(i3,i2)%ek0*quqd+p5178k0*cvqd+p6k0*cvqu
      cauxb=-cz324(i3,i2)%ek0*p6(2)+p6k0*cz324(i3,i2)%e(2)
      r6_324(i3,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_324(i3,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      end do
      end do
  
       if (ineutri(id6).ne.1) then
      do i3=1,2
      do i2=1,2
* TWR0 -- qu=p5178,qd=p6,v=cf324(i3,i2)%e,a=r6_324(i3,i2)%a,b=r6_324(i3,
* i2)%b,cl=fcl(id6),nsum=1
      ceps_0=-cf324(i3,i2)%ek0*(p5178(2)*p6(3)-p6(2)*p5178(3))+p
     & 5178k0*(cf324(i3,i2)%e(2)*p6(3)-p6(2)*cf324(i3,i2)%e(3))-
     & p6k0*(cf324(i3,i2)%e(2)*p5178(3)-p5178(2)*cf324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf324(i3,i2)%e(3)*p6k0+p6(3)*cf324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i3,i2)%e(0)*p5178(0)-cf324(i3,i2)%e(1)*p5178(1)
     & -cf324(i3,i2)%e(2)*p5178(2)-cf324(i3,i2)%e(3)*p5178(3)
      cvqd=cf324(i3,i2)%e(0)*p6(0)-cf324(i3,i2)%e(1)*p6(1)-cf324
     & (i3,i2)%e(2)*p6(2)-cf324(i3,i2)%e(3)*p6(3)
      cauxa=-cf324(i3,i2)%ek0*quqd+p5178k0*cvqd+p6k0*cvqu
      cauxb=-cf324(i3,i2)%ek0*p6(2)+p6k0*cf324(i3,i2)%e(2)
      r6_324(i3,i2)%a(2)=r6_324(i3,i2)%a(2)+fcl(id6)*(cauxa-ceps
     & _0)
      r6_324(i3,i2)%b(1)=r6_324(i3,i2)%b(1)+fcl(id6)*(cauxb-ceps
     & _2)
      end do
      end do
       endif
  
      endif
      endif
  
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p5,q=p678
      quqd=p5(0)*p678(0)-p5(1)*p678(1)-p5(2)*p678(2)-p5(3)*p678(
     & 3)
      ccl=zcl(id5)/(f678)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p678,v=cz3214(i3,i1,i2)%e,a=l5_3214(i3,i1,i2)%a,c=l5_
* 3214(i3,i1,i2)%c,cl=ccl,nsum=0
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
      l5_3214(i3,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_3214(i3,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
        if (ineutri(id5).ne.1) then
      ccl=fcl(id5)/(f678)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p678,v=cf3214(i3,i1,i2)%e,a=l5_3214(i3,i1,i2)%a,c=l5_
* 3214(i3,i1,i2)%c,cl=ccl,nsum=1
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
      l5_3214(i3,i1,i2)%a(2)=l5_3214(i3,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l5_3214(i3,i1,i2)%c(2)=l5_3214(i3,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
        endif
  
      if (ilept(id5).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p5,q=p5234
      quqd=p5(0)*p5234(0)-p5(1)*p5234(1)-p5(2)*p5234(2)-p5(3)*p5
     & 234(3)
      ccl=zcl(id5)/(f5234)
      do i3=1,2
      do i1=1,2
* TWL0 -- qu=p5,qd=p5234,v=cz324(i3,i1)%e,a=l5_324(i3,i1)%a,c=l5_324(i3,
* i1)%c,cl=ccl,nsum=0
      ceps_0=-cz324(i3,i1)%ek0*(p5(2)*p5234(3)-p5234(2)*p5(3))+p
     & 5k0*(cz324(i3,i1)%e(2)*p5234(3)-p5234(2)*cz324(i3,i1)%e(3
     & ))-p5234k0*(cz324(i3,i1)%e(2)*p5(3)-p5(2)*cz324(i3,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i3,i1)%e(3)*p5k0+p5(3)*cz324(i3,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz324(i3,i1)%e(0)*p5(0)-cz324(i3,i1)%e(1)*p5(1)-cz324
     & (i3,i1)%e(2)*p5(2)-cz324(i3,i1)%e(3)*p5(3)
      cvqd=cz324(i3,i1)%e(0)*p5234(0)-cz324(i3,i1)%e(1)*p5234(1)
     & -cz324(i3,i1)%e(2)*p5234(2)-cz324(i3,i1)%e(3)*p5234(3)
      cauxa=-cz324(i3,i1)%ek0*quqd+p5k0*cvqd+p5234k0*cvqu
      cauxc=+cz324(i3,i1)%ek0*p5(2)-p5k0*cz324(i3,i1)%e(2)
      l5_324(i3,i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_324(i3,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
       if (ineutri(id5).ne.1) then
      ccl=fcl(id5)/(f5234)
      do i3=1,2
      do i1=1,2
* TWL0 -- qu=p5,qd=p5234,v=cf324(i3,i1)%e,a=l5_324(i3,i1)%a,c=l5_324(i3,
* i1)%c,cl=ccl,nsum=1
      ceps_0=-cf324(i3,i1)%ek0*(p5(2)*p5234(3)-p5234(2)*p5(3))+p
     & 5k0*(cf324(i3,i1)%e(2)*p5234(3)-p5234(2)*cf324(i3,i1)%e(3
     & ))-p5234k0*(cf324(i3,i1)%e(2)*p5(3)-p5(2)*cf324(i3,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i3,i1)%e(3)*p5k0+p5(3)*cf324(i3,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf324(i3,i1)%e(0)*p5(0)-cf324(i3,i1)%e(1)*p5(1)-cf324
     & (i3,i1)%e(2)*p5(2)-cf324(i3,i1)%e(3)*p5(3)
      cvqd=cf324(i3,i1)%e(0)*p5234(0)-cf324(i3,i1)%e(1)*p5234(1)
     & -cf324(i3,i1)%e(2)*p5234(2)-cf324(i3,i1)%e(3)*p5234(3)
      cauxa=-cf324(i3,i1)%ek0*quqd+p5k0*cvqd+p5234k0*cvqu
      cauxc=+cf324(i3,i1)%ek0*p5(2)-p5k0*cf324(i3,i1)%e(2)
      l5_324(i3,i1)%a(2)=l5_324(i3,i1)%a(2)+ccl*(cauxa-ceps_0)
      l5_324(i3,i1)%c(2)=l5_324(i3,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
        endif
  
      endif
      endif
  
      if (ilept(id7).ne.1) then
* quqd -- p=p5,q=p634
      quqd=p5(0)*p634(0)-p5(1)*p634(1)-p5(2)*p634(2)-p5(3)*p634(
     & 3)
      ccl=wcl/(f634)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p634,v=cw7218(i1,i2)%e,a=l5_7218(i1,i2)%a,c=l5_7218(i
* 1,i2)%c,cl=ccl,nsum=0
      ceps_0=-cw7218(i1,i2)%ek0*(p5(2)*p634(3)-p634(2)*p5(3))+p5
     & k0*(cw7218(i1,i2)%e(2)*p634(3)-p634(2)*cw7218(i1,i2)%e(3)
     & )-p634k0*(cw7218(i1,i2)%e(2)*p5(3)-p5(2)*cw7218(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cw7218(i1,i2)%e(3)*p5k0+p5(3)*cw7218(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw7218(i1,i2)%e(0)*p5(0)-cw7218(i1,i2)%e(1)*p5(1)-cw7
     & 218(i1,i2)%e(2)*p5(2)-cw7218(i1,i2)%e(3)*p5(3)
      cvqd=cw7218(i1,i2)%e(0)*p634(0)-cw7218(i1,i2)%e(1)*p634(1)
     & -cw7218(i1,i2)%e(2)*p634(2)-cw7218(i1,i2)%e(3)*p634(3)
      cauxa=-cw7218(i1,i2)%ek0*quqd+p5k0*cvqd+p634k0*cvqu
      cauxc=+cw7218(i1,i2)%ek0*p5(2)-p5k0*cw7218(i1,i2)%e(2)
      l5_7218(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_7218(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      if (ilept(id5).ne.1.or.ilept(id3).ne.1) then
* quqd -- p=p5,q=p5278
      quqd=p5(0)*p5278(0)-p5(1)*p5278(1)-p5(2)*p5278(2)-p5(3)*p5
     & 278(3)
      ccl=wcl/(f5278)
      do i1=1,2
* TWL0 -- qu=p5,qd=p5278,v=cw728(i1)%e,a=l5_728(i1)%a,c=l5_728(i1)%c,cl=
* ccl,nsum=0
      ceps_0=-cw728(i1)%ek0*(p5(2)*p5278(3)-p5278(2)*p5(3))+p5k0
     & *(cw728(i1)%e(2)*p5278(3)-p5278(2)*cw728(i1)%e(3))-p5278k
     & 0*(cw728(i1)%e(2)*p5(3)-p5(2)*cw728(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw728(i1)%e(3)*p5k0+p5(3)*cw728(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw728(i1)%e(0)*p5(0)-cw728(i1)%e(1)*p5(1)-cw728(i1)%e
     & (2)*p5(2)-cw728(i1)%e(3)*p5(3)
      cvqd=cw728(i1)%e(0)*p5278(0)-cw728(i1)%e(1)*p5278(1)-cw728
     & (i1)%e(2)*p5278(2)-cw728(i1)%e(3)*p5278(3)
      cauxa=-cw728(i1)%ek0*quqd+p5k0*cvqd+p5278k0*cvqu
      cauxc=+cw728(i1)%ek0*p5(2)-p5k0*cw728(i1)%e(2)
      l5_728(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_728(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* quqd -- p=p534,q=p6
      quqd=p534(0)*p6(0)-p534(1)*p6(1)-p534(2)*p6(2)-p534(3)*p6(
     & 3)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p534,qd=p6,v=cw7218(i1,i2)%e,a=r6_7218(i1,i2)%a,b=r6_7218(i
* 1,i2)%b,cl=wcl,nsum=0
      ceps_0=-cw7218(i1,i2)%ek0*(p534(2)*p6(3)-p6(2)*p534(3))+p5
     & 34k0*(cw7218(i1,i2)%e(2)*p6(3)-p6(2)*cw7218(i1,i2)%e(3))-
     & p6k0*(cw7218(i1,i2)%e(2)*p534(3)-p534(2)*cw7218(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cw7218(i1,i2)%e(3)*p6k0+p6(3)*cw7218(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw7218(i1,i2)%e(0)*p534(0)-cw7218(i1,i2)%e(1)*p534(1)
     & -cw7218(i1,i2)%e(2)*p534(2)-cw7218(i1,i2)%e(3)*p534(3)
      cvqd=cw7218(i1,i2)%e(0)*p6(0)-cw7218(i1,i2)%e(1)*p6(1)-cw7
     & 218(i1,i2)%e(2)*p6(2)-cw7218(i1,i2)%e(3)*p6(3)
      cauxa=-cw7218(i1,i2)%ek0*quqd+p534k0*cvqd+p6k0*cvqu
      cauxb=-cw7218(i1,i2)%ek0*p6(2)+p6k0*cw7218(i1,i2)%e(2)
      r6_7218(i1,i2)%a(2)=wcl*(cauxa-ceps_0)
      r6_7218(i1,i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
      end do
  
      if (ilept(id6).ne.1.or.ilept(id3).ne.1) then
* quqd -- p=p5234,q=p6
      quqd=p5234(0)*p6(0)-p5234(1)*p6(1)-p5234(2)*p6(2)-p5234(3)
     & *p6(3)
      do i2=1,2
* TWR0 -- qu=p5234,qd=p6,v=cw718(i2)%e,a=r6_718(i2)%a,b=r6_718(i2)%b,cl=
* wcl,nsum=0
      ceps_0=-cw718(i2)%ek0*(p5234(2)*p6(3)-p6(2)*p5234(3))+p523
     & 4k0*(cw718(i2)%e(2)*p6(3)-p6(2)*cw718(i2)%e(3))-p6k0*(cw7
     & 18(i2)%e(2)*p5234(3)-p5234(2)*cw718(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw718(i2)%e(3)*p6k0+p6(3)*cw718(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw718(i2)%e(0)*p5234(0)-cw718(i2)%e(1)*p5234(1)-cw718
     & (i2)%e(2)*p5234(2)-cw718(i2)%e(3)*p5234(3)
      cvqd=cw718(i2)%e(0)*p6(0)-cw718(i2)%e(1)*p6(1)-cw718(i2)%e
     & (2)*p6(2)-cw718(i2)%e(3)*p6(3)
      cauxa=-cw718(i2)%ek0*quqd+p5234k0*cvqd+p6k0*cvqu
      cauxb=-cw718(i2)%ek0*p6(2)+p6k0*cw718(i2)%e(2)
      r6_718(i2)%a(2)=wcl*(cauxa-ceps_0)
      r6_718(i2)%b(1)=wcl*(cauxb-ceps_2)
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
* TWR0 -- qu=p578,qd=p6,v=cz3214(i3,i1,i2)%e,a=r6_3214(i3,i1,i2)%a,b=r6_
* 3214(i3,i1,i2)%b,cl=zcl(id6),nsum=0
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
      r6_3214(i3,i1,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_3214(i3,i1,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      end do
      end do
      end do
  
       if (ineutri(id6).ne.1) then
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p578,qd=p6,v=cf3214(i3,i1,i2)%e,a=r6_3214(i3,i1,i2)%a,b=r6_
* 3214(i3,i1,i2)%b,cl=fcl(id6),nsum=1
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
      r6_3214(i3,i1,i2)%a(2)=r6_3214(i3,i1,i2)%a(2)+fcl(id6)*(ca
     & uxa-ceps_0)
      r6_3214(i3,i1,i2)%b(1)=r6_3214(i3,i1,i2)%b(1)+fcl(id6)*(ca
     & uxb-ceps_2)
      end do
      end do
      end do
       endif
  
      if (ilept(id6).ne.1.or.ilept(id7).ne.1) then
* quqd -- p=p5278,q=p6
      quqd=p5278(0)*p6(0)-p5278(1)*p6(1)-p5278(2)*p6(2)-p5278(3)
     & *p6(3)
      do i3=1,2
      do i2=1,2
* TWR0 -- qu=p5278,qd=p6,v=cz314(i3,i2)%e,a=r6_314(i3,i2)%a,b=r6_314(i3,
* i2)%b,cl=zcl(id6),nsum=0
      ceps_0=-cz314(i3,i2)%ek0*(p5278(2)*p6(3)-p6(2)*p5278(3))+p
     & 5278k0*(cz314(i3,i2)%e(2)*p6(3)-p6(2)*cz314(i3,i2)%e(3))-
     & p6k0*(cz314(i3,i2)%e(2)*p5278(3)-p5278(2)*cz314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz314(i3,i2)%e(3)*p6k0+p6(3)*cz314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i3,i2)%e(0)*p5278(0)-cz314(i3,i2)%e(1)*p5278(1)
     & -cz314(i3,i2)%e(2)*p5278(2)-cz314(i3,i2)%e(3)*p5278(3)
      cvqd=cz314(i3,i2)%e(0)*p6(0)-cz314(i3,i2)%e(1)*p6(1)-cz314
     & (i3,i2)%e(2)*p6(2)-cz314(i3,i2)%e(3)*p6(3)
      cauxa=-cz314(i3,i2)%ek0*quqd+p5278k0*cvqd+p6k0*cvqu
      cauxb=-cz314(i3,i2)%ek0*p6(2)+p6k0*cz314(i3,i2)%e(2)
      r6_314(i3,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      r6_314(i3,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      end do
      end do
  
       if (ineutri(id6).ne.1) then
      do i3=1,2
      do i2=1,2
* TWR0 -- qu=p5278,qd=p6,v=cf314(i3,i2)%e,a=r6_314(i3,i2)%a,b=r6_314(i3,
* i2)%b,cl=fcl(id6),nsum=1
      ceps_0=-cf314(i3,i2)%ek0*(p5278(2)*p6(3)-p6(2)*p5278(3))+p
     & 5278k0*(cf314(i3,i2)%e(2)*p6(3)-p6(2)*cf314(i3,i2)%e(3))-
     & p6k0*(cf314(i3,i2)%e(2)*p5278(3)-p5278(2)*cf314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf314(i3,i2)%e(3)*p6k0+p6(3)*cf314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i3,i2)%e(0)*p5278(0)-cf314(i3,i2)%e(1)*p5278(1)
     & -cf314(i3,i2)%e(2)*p5278(2)-cf314(i3,i2)%e(3)*p5278(3)
      cvqd=cf314(i3,i2)%e(0)*p6(0)-cf314(i3,i2)%e(1)*p6(1)-cf314
     & (i3,i2)%e(2)*p6(2)-cf314(i3,i2)%e(3)*p6(3)
      cauxa=-cf314(i3,i2)%ek0*quqd+p5278k0*cvqd+p6k0*cvqu
      cauxb=-cf314(i3,i2)%ek0*p6(2)+p6k0*cf314(i3,i2)%e(2)
      r6_314(i3,i2)%a(2)=r6_314(i3,i2)%a(2)+fcl(id6)*(cauxa-ceps
     & _0)
      r6_314(i3,i2)%b(1)=r6_314(i3,i2)%b(1)+fcl(id6)*(cauxb-ceps
     & _2)
      end do
      end do
       endif
  
      endif
      endif
  
  
      if (ilept(id7).ne.1.or.ilept(id5).ne.1) then
      quqd=(s734-s34)/2d0
      ccl=zcl(id7)/(f734)
      do i3=1,2
* TWL0 -- qu=p7,qd=p734,v=cz34(i3)%e,a=l7_34(i3)%a,c=l7_34(i3)%c,cl=ccl,
* nsum=0
      ceps_0=-cz34(i3)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*(c
     & z34(i3)%e(2)*p734(3)-p734(2)*cz34(i3)%e(3))-p734k0*(cz34(
     & i3)%e(2)*p7(3)-p7(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p7k0+p7(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i3)%e(0)*p7(0)-cz34(i3)%e(1)*p7(1)-cz34(i3)%e(2)
     & *p7(2)-cz34(i3)%e(3)*p7(3)
      cvqd=cz34(i3)%e(0)*p734(0)-cz34(i3)%e(1)*p734(1)-cz34(i3)%
     & e(2)*p734(2)-cz34(i3)%e(3)*p734(3)
      cauxa=-cz34(i3)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cz34(i3)%ek0*p7(2)-p7k0*cz34(i3)%e(2)
      l7_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      l7_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
        if (ineutri(id7).ne.1) then
      ccl=fcl(id7)/(f734)
      do i3=1,2
* TWL0 -- qu=p7,qd=p734,v=cf34(i3)%e,a=l7_34(i3)%a,c=l7_34(i3)%c,cl=ccl,
* nsum=1
      ceps_0=-cf34(i3)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*(c
     & f34(i3)%e(2)*p734(3)-p734(2)*cf34(i3)%e(3))-p734k0*(cf34(
     & i3)%e(2)*p7(3)-p7(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p7k0+p7(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf34(i3)%e(0)*p7(0)-cf34(i3)%e(1)*p7(1)-cf34(i3)%e(2)
     & *p7(2)-cf34(i3)%e(3)*p7(3)
      cvqd=cf34(i3)%e(0)*p734(0)-cf34(i3)%e(1)*p734(1)-cf34(i3)%
     & e(2)*p734(2)-cf34(i3)%e(3)*p734(3)
      cauxa=-cf34(i3)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cf34(i3)%ek0*p7(2)-p7k0*cf34(i3)%e(2)
      l7_34(i3)%a(2)=l7_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      l7_34(i3)%c(2)=l7_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
        endif
  
      endif
  
      if (ilept(id7).ne.1.or.ilept(id3).ne.1) then
      quqd=(s756-s56)/2d0
      ccl=wcl/(f756)
* TWL0 -- qu=p7,qd=p756,v=cw56%e,a=l7_56%a,c=l7_56%c,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p7(2)*p756(3)-p756(2)*p7(3))+p7k0*(cw56%
     & e(2)*p756(3)-p756(2)*cw56%e(3))-p756k0*(cw56%e(2)*p7(3)-p
     & 7(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p7k0+p7(3)*cw56%ek0
      ceps_1=ceps_1*cim
      cvqu=cw56%e(0)*p7(0)-cw56%e(1)*p7(1)-cw56%e(2)*p7(2)-cw56%
     & e(3)*p7(3)
      cvqd=cw56%e(0)*p756(0)-cw56%e(1)*p756(1)-cw56%e(2)*p756(2)
     & -cw56%e(3)*p756(3)
      cauxa=-cw56%ek0*quqd+p7k0*cvqd+p756k0*cvqu
      cauxc=+cw56%ek0*p7(2)-p7k0*cw56%e(2)
      l7_56%a(2)=ccl*(cauxa-ceps_0)
      l7_56%c(2)=ccl*(-cauxc+ceps_1)
      endif
  
      if (ilept(id8).ne.1.or.ilept(id3).ne.1) then
      quqd=(-s856+s56)/2d0
* TWR0 -- qu=p856,qd=p8,v=cw56%e,a=r8_56%a,b=r8_56%b,cl=wcl,nsum=0
      ceps_0=-cw56%ek0*(p856(2)*p8(3)-p8(2)*p856(3))+p856k0*(cw5
     & 6%e(2)*p8(3)-p8(2)*cw56%e(3))-p8k0*(cw56%e(2)*p856(3)-p85
     & 6(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw56%e(3)*p8k0+p8(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p856(0)-cw56%e(1)*p856(1)-cw56%e(2)*p856(2)
     & -cw56%e(3)*p856(3)
      cvqd=cw56%e(0)*p8(0)-cw56%e(1)*p8(1)-cw56%e(2)*p8(2)-cw56%
     & e(3)*p8(3)
      cauxa=-cw56%ek0*quqd+p856k0*cvqd+p8k0*cvqu
      cauxb=-cw56%ek0*p8(2)+p8k0*cw56%e(2)
      r8_56%a(2)=wcl*(cauxa-ceps_0)
      r8_56%b(1)=wcl*(cauxb-ceps_2)
      endif
  
      if (ilept(id8).ne.1.or.ilept(id5).ne.1) then
      quqd=(-s834+s34)/2d0
      do i3=1,2
* TWR0 -- qu=p834,qd=p8,v=cz34(i3)%e,a=r8_34(i3)%a,b=r8_34(i3)%b,cl=zcl(
* id8),nsum=0
      ceps_0=-cz34(i3)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834k0*
     & (cz34(i3)%e(2)*p8(3)-p8(2)*cz34(i3)%e(3))-p8k0*(cz34(i3)%
     & e(2)*p834(3)-p834(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i3)%e(3)*p8k0+p8(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p834(0)-cz34(i3)%e(1)*p834(1)-cz34(i3)%
     & e(2)*p834(2)-cz34(i3)%e(3)*p834(3)
      cvqd=cz34(i3)%e(0)*p8(0)-cz34(i3)%e(1)*p8(1)-cz34(i3)%e(2)
     & *p8(2)-cz34(i3)%e(3)*p8(3)
      cauxa=-cz34(i3)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cz34(i3)%ek0*p8(2)+p8k0*cz34(i3)%e(2)
      r8_34(i3)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_34(i3)%b(1)=zcl(id8)*(cauxb-ceps_2)
      end do
  
        if (ineutri(id8).ne.1) then
      do i3=1,2
* TWR0 -- qu=p834,qd=p8,v=cf34(i3)%e,a=r8_34(i3)%a,b=r8_34(i3)%b,cl=fcl(
* id8),nsum=1
      ceps_0=-cf34(i3)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834k0*
     & (cf34(i3)%e(2)*p8(3)-p8(2)*cf34(i3)%e(3))-p8k0*(cf34(i3)%
     & e(2)*p834(3)-p834(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf34(i3)%e(3)*p8k0+p8(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p834(0)-cf34(i3)%e(1)*p834(1)-cf34(i3)%
     & e(2)*p834(2)-cf34(i3)%e(3)*p834(3)
      cvqd=cf34(i3)%e(0)*p8(0)-cf34(i3)%e(1)*p8(1)-cf34(i3)%e(2)
     & *p8(2)-cf34(i3)%e(3)*p8(3)
      cauxa=-cf34(i3)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cf34(i3)%ek0*p8(2)+p8k0*cf34(i3)%e(2)
      r8_34(i3)%a(2)=r8_34(i3)%a(2)+fcl(id8)*(cauxa-ceps_0)
      r8_34(i3)%b(1)=r8_34(i3)%b(1)+fcl(id8)*(cauxb-ceps_2)
      end do
        endif
      endif
  
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p7,q=p856
      quqd=p7(0)*p856(0)-p7(1)*p856(1)-p7(2)*p856(2)-p7(3)*p856(
     & 3)
      ccl=zcl(id7)/(f856)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p856,v=cz3124(i3,i1,i2)%e,a=l7_3124(i3,i1,i2)%a,c=l7_
* 3124(i3,i1,i2)%c,cl=ccl,nsum=0
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
      l7_3124(i3,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_3124(i3,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
        if (ineutri(id7).ne.1) then
      ccl=fcl(id7)/(f856)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p856,v=cf3124(i3,i1,i2)%e,a=l7_3124(i3,i1,i2)%a,c=l7_
* 3124(i3,i1,i2)%c,cl=ccl,nsum=1
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
      l7_3124(i3,i1,i2)%a(2)=l7_3124(i3,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l7_3124(i3,i1,i2)%c(2)=l7_3124(i3,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
        endif
  
      if (ilept(id7).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p7,q=p7134
      quqd=p7(0)*p7134(0)-p7(1)*p7134(1)-p7(2)*p7134(2)-p7(3)*p7
     & 134(3)
      ccl=zcl(id7)/(f7134)
      do i3=1,2
      do i1=1,2
* TWL0 -- qu=p7,qd=p7134,v=cz314(i3,i1)%e,a=l7_314(i3,i1)%a,c=l7_314(i3,
* i1)%c,cl=ccl,nsum=0
      ceps_0=-cz314(i3,i1)%ek0*(p7(2)*p7134(3)-p7134(2)*p7(3))+p
     & 7k0*(cz314(i3,i1)%e(2)*p7134(3)-p7134(2)*cz314(i3,i1)%e(3
     & ))-p7134k0*(cz314(i3,i1)%e(2)*p7(3)-p7(2)*cz314(i3,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i3,i1)%e(3)*p7k0+p7(3)*cz314(i3,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz314(i3,i1)%e(0)*p7(0)-cz314(i3,i1)%e(1)*p7(1)-cz314
     & (i3,i1)%e(2)*p7(2)-cz314(i3,i1)%e(3)*p7(3)
      cvqd=cz314(i3,i1)%e(0)*p7134(0)-cz314(i3,i1)%e(1)*p7134(1)
     & -cz314(i3,i1)%e(2)*p7134(2)-cz314(i3,i1)%e(3)*p7134(3)
      cauxa=-cz314(i3,i1)%ek0*quqd+p7k0*cvqd+p7134k0*cvqu
      cauxc=+cz314(i3,i1)%ek0*p7(2)-p7k0*cz314(i3,i1)%e(2)
      l7_314(i3,i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_314(i3,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
       if (ineutri(id7).ne.1) then
      ccl=fcl(id7)/(f7134)
      do i3=1,2
      do i1=1,2
* TWL0 -- qu=p7,qd=p7134,v=cf314(i3,i1)%e,a=l7_314(i3,i1)%a,c=l7_314(i3,
* i1)%c,cl=ccl,nsum=1
      ceps_0=-cf314(i3,i1)%ek0*(p7(2)*p7134(3)-p7134(2)*p7(3))+p
     & 7k0*(cf314(i3,i1)%e(2)*p7134(3)-p7134(2)*cf314(i3,i1)%e(3
     & ))-p7134k0*(cf314(i3,i1)%e(2)*p7(3)-p7(2)*cf314(i3,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i3,i1)%e(3)*p7k0+p7(3)*cf314(i3,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf314(i3,i1)%e(0)*p7(0)-cf314(i3,i1)%e(1)*p7(1)-cf314
     & (i3,i1)%e(2)*p7(2)-cf314(i3,i1)%e(3)*p7(3)
      cvqd=cf314(i3,i1)%e(0)*p7134(0)-cf314(i3,i1)%e(1)*p7134(1)
     & -cf314(i3,i1)%e(2)*p7134(2)-cf314(i3,i1)%e(3)*p7134(3)
      cauxa=-cf314(i3,i1)%ek0*quqd+p7k0*cvqd+p7134k0*cvqu
      cauxc=+cf314(i3,i1)%ek0*p7(2)-p7k0*cf314(i3,i1)%e(2)
      l7_314(i3,i1)%a(2)=l7_314(i3,i1)%a(2)+ccl*(cauxa-ceps_0)
      l7_314(i3,i1)%c(2)=l7_314(i3,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
        endif
  
      endif
      endif
  
      if (ilept(id5).ne.1) then
* quqd -- p=p7,q=p834
      quqd=p7(0)*p834(0)-p7(1)*p834(1)-p7(2)*p834(2)-p7(3)*p834(
     & 3)
      ccl=wcl/(f834)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p834,v=cw5126(i1,i2)%e,a=l7_5126(i1,i2)%a,c=l7_5126(i
* 1,i2)%c,cl=ccl,nsum=0
      ceps_0=-cw5126(i1,i2)%ek0*(p7(2)*p834(3)-p834(2)*p7(3))+p7
     & k0*(cw5126(i1,i2)%e(2)*p834(3)-p834(2)*cw5126(i1,i2)%e(3)
     & )-p834k0*(cw5126(i1,i2)%e(2)*p7(3)-p7(2)*cw5126(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cw5126(i1,i2)%e(3)*p7k0+p7(3)*cw5126(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw5126(i1,i2)%e(0)*p7(0)-cw5126(i1,i2)%e(1)*p7(1)-cw5
     & 126(i1,i2)%e(2)*p7(2)-cw5126(i1,i2)%e(3)*p7(3)
      cvqd=cw5126(i1,i2)%e(0)*p834(0)-cw5126(i1,i2)%e(1)*p834(1)
     & -cw5126(i1,i2)%e(2)*p834(2)-cw5126(i1,i2)%e(3)*p834(3)
      cauxa=-cw5126(i1,i2)%ek0*quqd+p7k0*cvqd+p834k0*cvqu
      cauxc=+cw5126(i1,i2)%ek0*p7(2)-p7k0*cw5126(i1,i2)%e(2)
      l7_5126(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_5126(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      if (ilept(id7).ne.1.or.ilept(id3).ne.1) then
* quqd -- p=p7,q=p7156
      quqd=p7(0)*p7156(0)-p7(1)*p7156(1)-p7(2)*p7156(2)-p7(3)*p7
     & 156(3)
      ccl=wcl/(f7156)
      do i1=1,2
* TWL0 -- qu=p7,qd=p7156,v=cw516(i1)%e,a=l7_516(i1)%a,c=l7_516(i1)%c,cl=
* ccl,nsum=0
      ceps_0=-cw516(i1)%ek0*(p7(2)*p7156(3)-p7156(2)*p7(3))+p7k0
     & *(cw516(i1)%e(2)*p7156(3)-p7156(2)*cw516(i1)%e(3))-p7156k
     & 0*(cw516(i1)%e(2)*p7(3)-p7(2)*cw516(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw516(i1)%e(3)*p7k0+p7(3)*cw516(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw516(i1)%e(0)*p7(0)-cw516(i1)%e(1)*p7(1)-cw516(i1)%e
     & (2)*p7(2)-cw516(i1)%e(3)*p7(3)
      cvqd=cw516(i1)%e(0)*p7156(0)-cw516(i1)%e(1)*p7156(1)-cw516
     & (i1)%e(2)*p7156(2)-cw516(i1)%e(3)*p7156(3)
      cauxa=-cw516(i1)%ek0*quqd+p7k0*cvqd+p7156k0*cvqu
      cauxc=+cw516(i1)%ek0*p7(2)-p7k0*cw516(i1)%e(2)
      l7_516(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_516(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
* quqd -- p=p734,q=p8
      quqd=p734(0)*p8(0)-p734(1)*p8(1)-p734(2)*p8(2)-p734(3)*p8(
     & 3)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p734,qd=p8,v=cw5126(i1,i2)%e,a=r8_5126(i1,i2)%a,b=r8_5126(i
* 1,i2)%b,cl=wcl,nsum=0
      ceps_0=-cw5126(i1,i2)%ek0*(p734(2)*p8(3)-p8(2)*p734(3))+p7
     & 34k0*(cw5126(i1,i2)%e(2)*p8(3)-p8(2)*cw5126(i1,i2)%e(3))-
     & p8k0*(cw5126(i1,i2)%e(2)*p734(3)-p734(2)*cw5126(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cw5126(i1,i2)%e(3)*p8k0+p8(3)*cw5126(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw5126(i1,i2)%e(0)*p734(0)-cw5126(i1,i2)%e(1)*p734(1)
     & -cw5126(i1,i2)%e(2)*p734(2)-cw5126(i1,i2)%e(3)*p734(3)
      cvqd=cw5126(i1,i2)%e(0)*p8(0)-cw5126(i1,i2)%e(1)*p8(1)-cw5
     & 126(i1,i2)%e(2)*p8(2)-cw5126(i1,i2)%e(3)*p8(3)
      cauxa=-cw5126(i1,i2)%ek0*quqd+p734k0*cvqd+p8k0*cvqu
      cauxb=-cw5126(i1,i2)%ek0*p8(2)+p8k0*cw5126(i1,i2)%e(2)
      r8_5126(i1,i2)%a(2)=wcl*(cauxa-ceps_0)
      r8_5126(i1,i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
      end do
  
      if (ilept(id8).ne.1.or.ilept(id3).ne.1) then
* quqd -- p=p7134,q=p8
      quqd=p7134(0)*p8(0)-p7134(1)*p8(1)-p7134(2)*p8(2)-p7134(3)
     & *p8(3)
      do i2=1,2
* TWR0 -- qu=p7134,qd=p8,v=cw526(i2)%e,a=r8_526(i2)%a,b=r8_526(i2)%b,cl=
* wcl,nsum=0
      ceps_0=-cw526(i2)%ek0*(p7134(2)*p8(3)-p8(2)*p7134(3))+p713
     & 4k0*(cw526(i2)%e(2)*p8(3)-p8(2)*cw526(i2)%e(3))-p8k0*(cw5
     & 26(i2)%e(2)*p7134(3)-p7134(2)*cw526(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw526(i2)%e(3)*p8k0+p8(3)*cw526(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw526(i2)%e(0)*p7134(0)-cw526(i2)%e(1)*p7134(1)-cw526
     & (i2)%e(2)*p7134(2)-cw526(i2)%e(3)*p7134(3)
      cvqd=cw526(i2)%e(0)*p8(0)-cw526(i2)%e(1)*p8(1)-cw526(i2)%e
     & (2)*p8(2)-cw526(i2)%e(3)*p8(3)
      cauxa=-cw526(i2)%ek0*quqd+p7134k0*cvqd+p8k0*cvqu
      cauxb=-cw526(i2)%ek0*p8(2)+p8k0*cw526(i2)%e(2)
      r8_526(i2)%a(2)=wcl*(cauxa-ceps_0)
      r8_526(i2)%b(1)=wcl*(cauxb-ceps_2)
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
* TWR0 -- qu=p756,qd=p8,v=cz3124(i3,i1,i2)%e,a=r8_3124(i3,i1,i2)%a,b=r8_
* 3124(i3,i1,i2)%b,cl=zcl(id8),nsum=0
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
      r8_3124(i3,i1,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_3124(i3,i1,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      end do
      end do
      end do
  
       if (ineutri(id8).ne.1) then
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p756,qd=p8,v=cf3124(i3,i1,i2)%e,a=r8_3124(i3,i1,i2)%a,b=r8_
* 3124(i3,i1,i2)%b,cl=fcl(id8),nsum=1
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
      r8_3124(i3,i1,i2)%a(2)=r8_3124(i3,i1,i2)%a(2)+fcl(id8)*(ca
     & uxa-ceps_0)
      r8_3124(i3,i1,i2)%b(1)=r8_3124(i3,i1,i2)%b(1)+fcl(id8)*(ca
     & uxb-ceps_2)
      end do
      end do
      end do
       endif
  
      if (ilept(id8).ne.1.or.ilept(id5).ne.1) then
* quqd -- p=p7156,q=p8
      quqd=p7156(0)*p8(0)-p7156(1)*p8(1)-p7156(2)*p8(2)-p7156(3)
     & *p8(3)
      do i3=1,2
      do i2=1,2
* TWR0 -- qu=p7156,qd=p8,v=cz324(i3,i2)%e,a=r8_324(i3,i2)%a,b=r8_324(i3,
* i2)%b,cl=zcl(id8),nsum=0
      ceps_0=-cz324(i3,i2)%ek0*(p7156(2)*p8(3)-p8(2)*p7156(3))+p
     & 7156k0*(cz324(i3,i2)%e(2)*p8(3)-p8(2)*cz324(i3,i2)%e(3))-
     & p8k0*(cz324(i3,i2)%e(2)*p7156(3)-p7156(2)*cz324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz324(i3,i2)%e(3)*p8k0+p8(3)*cz324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i3,i2)%e(0)*p7156(0)-cz324(i3,i2)%e(1)*p7156(1)
     & -cz324(i3,i2)%e(2)*p7156(2)-cz324(i3,i2)%e(3)*p7156(3)
      cvqd=cz324(i3,i2)%e(0)*p8(0)-cz324(i3,i2)%e(1)*p8(1)-cz324
     & (i3,i2)%e(2)*p8(2)-cz324(i3,i2)%e(3)*p8(3)
      cauxa=-cz324(i3,i2)%ek0*quqd+p7156k0*cvqd+p8k0*cvqu
      cauxb=-cz324(i3,i2)%ek0*p8(2)+p8k0*cz324(i3,i2)%e(2)
      r8_324(i3,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_324(i3,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      end do
      end do
  
       if (ineutri(id8).ne.1) then
      do i3=1,2
      do i2=1,2
* TWR0 -- qu=p7156,qd=p8,v=cf324(i3,i2)%e,a=r8_324(i3,i2)%a,b=r8_324(i3,
* i2)%b,cl=fcl(id8),nsum=1
      ceps_0=-cf324(i3,i2)%ek0*(p7156(2)*p8(3)-p8(2)*p7156(3))+p
     & 7156k0*(cf324(i3,i2)%e(2)*p8(3)-p8(2)*cf324(i3,i2)%e(3))-
     & p8k0*(cf324(i3,i2)%e(2)*p7156(3)-p7156(2)*cf324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf324(i3,i2)%e(3)*p8k0+p8(3)*cf324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i3,i2)%e(0)*p7156(0)-cf324(i3,i2)%e(1)*p7156(1)
     & -cf324(i3,i2)%e(2)*p7156(2)-cf324(i3,i2)%e(3)*p7156(3)
      cvqd=cf324(i3,i2)%e(0)*p8(0)-cf324(i3,i2)%e(1)*p8(1)-cf324
     & (i3,i2)%e(2)*p8(2)-cf324(i3,i2)%e(3)*p8(3)
      cauxa=-cf324(i3,i2)%ek0*quqd+p7156k0*cvqd+p8k0*cvqu
      cauxb=-cf324(i3,i2)%ek0*p8(2)+p8k0*cf324(i3,i2)%e(2)
      r8_324(i3,i2)%a(2)=r8_324(i3,i2)%a(2)+fcl(id8)*(cauxa-ceps
     & _0)
      r8_324(i3,i2)%b(1)=r8_324(i3,i2)%b(1)+fcl(id8)*(cauxb-ceps
     & _2)
      end do
      end do
       endif
  
      endif
      endif
  
  
      if (ilept(id3).ne.1) then
  
* quqd -- p=p7,q=p856
      quqd=p7(0)*p856(0)-p7(1)*p856(1)-p7(2)*p856(2)-p7(3)*p856(
     & 3)
      ccl=zcl(id7)/(f856)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p856,v=cz3214(i3,i1,i2)%e,a=l7_3214(i3,i1,i2)%a,c=l7_
* 3214(i3,i1,i2)%c,cl=ccl,nsum=0
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
      l7_3214(i3,i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_3214(i3,i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      end do
  
        if (ineutri(id7).ne.1) then
      ccl=fcl(id7)/(f856)
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p856,v=cf3214(i3,i1,i2)%e,a=l7_3214(i3,i1,i2)%a,c=l7_
* 3214(i3,i1,i2)%c,cl=ccl,nsum=1
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
      l7_3214(i3,i1,i2)%a(2)=l7_3214(i3,i1,i2)%a(2)+ccl*(cauxa-c
     & eps_0)
      l7_3214(i3,i1,i2)%c(2)=l7_3214(i3,i1,i2)%c(2)+ccl*(-cauxc+
     & ceps_1)
      end do
      end do
      end do
        endif
  
      if (ilept(id7).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p7,q=p7234
      quqd=p7(0)*p7234(0)-p7(1)*p7234(1)-p7(2)*p7234(2)-p7(3)*p7
     & 234(3)
      ccl=zcl(id7)/(f7234)
      do i3=1,2
      do i1=1,2
* TWL0 -- qu=p7,qd=p7234,v=cz324(i3,i1)%e,a=l7_324(i3,i1)%a,c=l7_324(i3,
* i1)%c,cl=ccl,nsum=0
      ceps_0=-cz324(i3,i1)%ek0*(p7(2)*p7234(3)-p7234(2)*p7(3))+p
     & 7k0*(cz324(i3,i1)%e(2)*p7234(3)-p7234(2)*cz324(i3,i1)%e(3
     & ))-p7234k0*(cz324(i3,i1)%e(2)*p7(3)-p7(2)*cz324(i3,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i3,i1)%e(3)*p7k0+p7(3)*cz324(i3,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz324(i3,i1)%e(0)*p7(0)-cz324(i3,i1)%e(1)*p7(1)-cz324
     & (i3,i1)%e(2)*p7(2)-cz324(i3,i1)%e(3)*p7(3)
      cvqd=cz324(i3,i1)%e(0)*p7234(0)-cz324(i3,i1)%e(1)*p7234(1)
     & -cz324(i3,i1)%e(2)*p7234(2)-cz324(i3,i1)%e(3)*p7234(3)
      cauxa=-cz324(i3,i1)%ek0*quqd+p7k0*cvqd+p7234k0*cvqu
      cauxc=+cz324(i3,i1)%ek0*p7(2)-p7k0*cz324(i3,i1)%e(2)
      l7_324(i3,i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_324(i3,i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
       if (ineutri(id7).ne.1) then
      ccl=fcl(id7)/(f7234)
      do i3=1,2
      do i1=1,2
* TWL0 -- qu=p7,qd=p7234,v=cf324(i3,i1)%e,a=l7_324(i3,i1)%a,c=l7_324(i3,
* i1)%c,cl=ccl,nsum=1
      ceps_0=-cf324(i3,i1)%ek0*(p7(2)*p7234(3)-p7234(2)*p7(3))+p
     & 7k0*(cf324(i3,i1)%e(2)*p7234(3)-p7234(2)*cf324(i3,i1)%e(3
     & ))-p7234k0*(cf324(i3,i1)%e(2)*p7(3)-p7(2)*cf324(i3,i1)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i3,i1)%e(3)*p7k0+p7(3)*cf324(i3,i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf324(i3,i1)%e(0)*p7(0)-cf324(i3,i1)%e(1)*p7(1)-cf324
     & (i3,i1)%e(2)*p7(2)-cf324(i3,i1)%e(3)*p7(3)
      cvqd=cf324(i3,i1)%e(0)*p7234(0)-cf324(i3,i1)%e(1)*p7234(1)
     & -cf324(i3,i1)%e(2)*p7234(2)-cf324(i3,i1)%e(3)*p7234(3)
      cauxa=-cf324(i3,i1)%ek0*quqd+p7k0*cvqd+p7234k0*cvqu
      cauxc=+cf324(i3,i1)%ek0*p7(2)-p7k0*cf324(i3,i1)%e(2)
      l7_324(i3,i1)%a(2)=l7_324(i3,i1)%a(2)+ccl*(cauxa-ceps_0)
      l7_324(i3,i1)%c(2)=l7_324(i3,i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
        endif
  
      endif
      endif
  
      if (ilept(id5).ne.1) then
* quqd -- p=p7,q=p834
      quqd=p7(0)*p834(0)-p7(1)*p834(1)-p7(2)*p834(2)-p7(3)*p834(
     & 3)
      ccl=wcl/(f834)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p834,v=cw5216(i1,i2)%e,a=l7_5216(i1,i2)%a,c=l7_5216(i
* 1,i2)%c,cl=ccl,nsum=0
      ceps_0=-cw5216(i1,i2)%ek0*(p7(2)*p834(3)-p834(2)*p7(3))+p7
     & k0*(cw5216(i1,i2)%e(2)*p834(3)-p834(2)*cw5216(i1,i2)%e(3)
     & )-p834k0*(cw5216(i1,i2)%e(2)*p7(3)-p7(2)*cw5216(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cw5216(i1,i2)%e(3)*p7k0+p7(3)*cw5216(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw5216(i1,i2)%e(0)*p7(0)-cw5216(i1,i2)%e(1)*p7(1)-cw5
     & 216(i1,i2)%e(2)*p7(2)-cw5216(i1,i2)%e(3)*p7(3)
      cvqd=cw5216(i1,i2)%e(0)*p834(0)-cw5216(i1,i2)%e(1)*p834(1)
     & -cw5216(i1,i2)%e(2)*p834(2)-cw5216(i1,i2)%e(3)*p834(3)
      cauxa=-cw5216(i1,i2)%ek0*quqd+p7k0*cvqd+p834k0*cvqu
      cauxc=+cw5216(i1,i2)%ek0*p7(2)-p7k0*cw5216(i1,i2)%e(2)
      l7_5216(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_5216(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      if (ilept(id7).ne.1.or.ilept(id3).ne.1) then
* quqd -- p=p7,q=p7256
      quqd=p7(0)*p7256(0)-p7(1)*p7256(1)-p7(2)*p7256(2)-p7(3)*p7
     & 256(3)
      ccl=wcl/(f7256)
      do i1=1,2
* TWL0 -- qu=p7,qd=p7256,v=cw526(i1)%e,a=l7_526(i1)%a,c=l7_526(i1)%c,cl=
* ccl,nsum=0
      ceps_0=-cw526(i1)%ek0*(p7(2)*p7256(3)-p7256(2)*p7(3))+p7k0
     & *(cw526(i1)%e(2)*p7256(3)-p7256(2)*cw526(i1)%e(3))-p7256k
     & 0*(cw526(i1)%e(2)*p7(3)-p7(2)*cw526(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw526(i1)%e(3)*p7k0+p7(3)*cw526(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw526(i1)%e(0)*p7(0)-cw526(i1)%e(1)*p7(1)-cw526(i1)%e
     & (2)*p7(2)-cw526(i1)%e(3)*p7(3)
      cvqd=cw526(i1)%e(0)*p7256(0)-cw526(i1)%e(1)*p7256(1)-cw526
     & (i1)%e(2)*p7256(2)-cw526(i1)%e(3)*p7256(3)
      cauxa=-cw526(i1)%ek0*quqd+p7k0*cvqd+p7256k0*cvqu
      cauxc=+cw526(i1)%ek0*p7(2)-p7k0*cw526(i1)%e(2)
      l7_526(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_526(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
* quqd -- p=p734,q=p8
      quqd=p734(0)*p8(0)-p734(1)*p8(1)-p734(2)*p8(2)-p734(3)*p8(
     & 3)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p734,qd=p8,v=cw5216(i1,i2)%e,a=r8_5216(i1,i2)%a,b=r8_5216(i
* 1,i2)%b,cl=wcl,nsum=0
      ceps_0=-cw5216(i1,i2)%ek0*(p734(2)*p8(3)-p8(2)*p734(3))+p7
     & 34k0*(cw5216(i1,i2)%e(2)*p8(3)-p8(2)*cw5216(i1,i2)%e(3))-
     & p8k0*(cw5216(i1,i2)%e(2)*p734(3)-p734(2)*cw5216(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cw5216(i1,i2)%e(3)*p8k0+p8(3)*cw5216(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw5216(i1,i2)%e(0)*p734(0)-cw5216(i1,i2)%e(1)*p734(1)
     & -cw5216(i1,i2)%e(2)*p734(2)-cw5216(i1,i2)%e(3)*p734(3)
      cvqd=cw5216(i1,i2)%e(0)*p8(0)-cw5216(i1,i2)%e(1)*p8(1)-cw5
     & 216(i1,i2)%e(2)*p8(2)-cw5216(i1,i2)%e(3)*p8(3)
      cauxa=-cw5216(i1,i2)%ek0*quqd+p734k0*cvqd+p8k0*cvqu
      cauxb=-cw5216(i1,i2)%ek0*p8(2)+p8k0*cw5216(i1,i2)%e(2)
      r8_5216(i1,i2)%a(2)=wcl*(cauxa-ceps_0)
      r8_5216(i1,i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
      end do
  
      if (ilept(id8).ne.1.or.ilept(id3).ne.1) then
* quqd -- p=p7234,q=p8
      quqd=p7234(0)*p8(0)-p7234(1)*p8(1)-p7234(2)*p8(2)-p7234(3)
     & *p8(3)
      do i2=1,2
* TWR0 -- qu=p7234,qd=p8,v=cw516(i2)%e,a=r8_516(i2)%a,b=r8_516(i2)%b,cl=
* wcl,nsum=0
      ceps_0=-cw516(i2)%ek0*(p7234(2)*p8(3)-p8(2)*p7234(3))+p723
     & 4k0*(cw516(i2)%e(2)*p8(3)-p8(2)*cw516(i2)%e(3))-p8k0*(cw5
     & 16(i2)%e(2)*p7234(3)-p7234(2)*cw516(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw516(i2)%e(3)*p8k0+p8(3)*cw516(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw516(i2)%e(0)*p7234(0)-cw516(i2)%e(1)*p7234(1)-cw516
     & (i2)%e(2)*p7234(2)-cw516(i2)%e(3)*p7234(3)
      cvqd=cw516(i2)%e(0)*p8(0)-cw516(i2)%e(1)*p8(1)-cw516(i2)%e
     & (2)*p8(2)-cw516(i2)%e(3)*p8(3)
      cauxa=-cw516(i2)%ek0*quqd+p7234k0*cvqd+p8k0*cvqu
      cauxb=-cw516(i2)%ek0*p8(2)+p8k0*cw516(i2)%e(2)
      r8_516(i2)%a(2)=wcl*(cauxa-ceps_0)
      r8_516(i2)%b(1)=wcl*(cauxb-ceps_2)
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
* TWR0 -- qu=p756,qd=p8,v=cz3214(i3,i1,i2)%e,a=r8_3214(i3,i1,i2)%a,b=r8_
* 3214(i3,i1,i2)%b,cl=zcl(id8),nsum=0
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
      r8_3214(i3,i1,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_3214(i3,i1,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      end do
      end do
      end do
  
       if (ineutri(id8).ne.1) then
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p756,qd=p8,v=cf3214(i3,i1,i2)%e,a=r8_3214(i3,i1,i2)%a,b=r8_
* 3214(i3,i1,i2)%b,cl=fcl(id8),nsum=1
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
      r8_3214(i3,i1,i2)%a(2)=r8_3214(i3,i1,i2)%a(2)+fcl(id8)*(ca
     & uxa-ceps_0)
      r8_3214(i3,i1,i2)%b(1)=r8_3214(i3,i1,i2)%b(1)+fcl(id8)*(ca
     & uxb-ceps_2)
      end do
      end do
      end do
       endif
  
      if (ilept(id8).ne.1.or.ilept(id5).ne.1) then
* quqd -- p=p7256,q=p8
      quqd=p7256(0)*p8(0)-p7256(1)*p8(1)-p7256(2)*p8(2)-p7256(3)
     & *p8(3)
      do i3=1,2
      do i2=1,2
* TWR0 -- qu=p7256,qd=p8,v=cz314(i3,i2)%e,a=r8_314(i3,i2)%a,b=r8_314(i3,
* i2)%b,cl=zcl(id8),nsum=0
      ceps_0=-cz314(i3,i2)%ek0*(p7256(2)*p8(3)-p8(2)*p7256(3))+p
     & 7256k0*(cz314(i3,i2)%e(2)*p8(3)-p8(2)*cz314(i3,i2)%e(3))-
     & p8k0*(cz314(i3,i2)%e(2)*p7256(3)-p7256(2)*cz314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cz314(i3,i2)%e(3)*p8k0+p8(3)*cz314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i3,i2)%e(0)*p7256(0)-cz314(i3,i2)%e(1)*p7256(1)
     & -cz314(i3,i2)%e(2)*p7256(2)-cz314(i3,i2)%e(3)*p7256(3)
      cvqd=cz314(i3,i2)%e(0)*p8(0)-cz314(i3,i2)%e(1)*p8(1)-cz314
     & (i3,i2)%e(2)*p8(2)-cz314(i3,i2)%e(3)*p8(3)
      cauxa=-cz314(i3,i2)%ek0*quqd+p7256k0*cvqd+p8k0*cvqu
      cauxb=-cz314(i3,i2)%ek0*p8(2)+p8k0*cz314(i3,i2)%e(2)
      r8_314(i3,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      r8_314(i3,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      end do
      end do
  
       if (ineutri(id8).ne.1) then
      do i3=1,2
      do i2=1,2
* TWR0 -- qu=p7256,qd=p8,v=cf314(i3,i2)%e,a=r8_314(i3,i2)%a,b=r8_314(i3,
* i2)%b,cl=fcl(id8),nsum=1
      ceps_0=-cf314(i3,i2)%ek0*(p7256(2)*p8(3)-p8(2)*p7256(3))+p
     & 7256k0*(cf314(i3,i2)%e(2)*p8(3)-p8(2)*cf314(i3,i2)%e(3))-
     & p8k0*(cf314(i3,i2)%e(2)*p7256(3)-p7256(2)*cf314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cf314(i3,i2)%e(3)*p8k0+p8(3)*cf314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i3,i2)%e(0)*p7256(0)-cf314(i3,i2)%e(1)*p7256(1)
     & -cf314(i3,i2)%e(2)*p7256(2)-cf314(i3,i2)%e(3)*p7256(3)
      cvqd=cf314(i3,i2)%e(0)*p8(0)-cf314(i3,i2)%e(1)*p8(1)-cf314
     & (i3,i2)%e(2)*p8(2)-cf314(i3,i2)%e(3)*p8(3)
      cauxa=-cf314(i3,i2)%ek0*quqd+p7256k0*cvqd+p8k0*cvqu
      cauxb=-cf314(i3,i2)%ek0*p8(2)+p8k0*cf314(i3,i2)%e(2)
      r8_314(i3,i2)%a(2)=r8_314(i3,i2)%a(2)+fcl(id8)*(cauxa-ceps
     & _0)
      r8_314(i3,i2)%b(1)=r8_314(i3,i2)%b(1)+fcl(id8)*(cauxb-ceps
     & _2)
      end do
      end do
       endif
  
      endif
      endif
  
  
  
      if (iup(id3).eq.1) then
  
      if (ilept(id3).ne.1.or.ilept(id5).ne.1) then
      quqd=(s378-s78)/2d0
      ccl=wcl/(f378)
* TWL0 -- qu=p3,qd=p378,v=cw78%e,a=l3_78%a,c=l3_78%c,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p3(2)*p378(3)-p378(2)*p3(3))+p3k0*(cw78%
     & e(2)*p378(3)-p378(2)*cw78%e(3))-p378k0*(cw78%e(2)*p3(3)-p
     & 3(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p3k0+p3(3)*cw78%ek0
      ceps_1=ceps_1*cim
      cvqu=cw78%e(0)*p3(0)-cw78%e(1)*p3(1)-cw78%e(2)*p3(2)-cw78%
     & e(3)*p3(3)
      cvqd=cw78%e(0)*p378(0)-cw78%e(1)*p378(1)-cw78%e(2)*p378(2)
     & -cw78%e(3)*p378(3)
      cauxa=-cw78%ek0*quqd+p3k0*cvqd+p378k0*cvqu
      cauxc=+cw78%ek0*p3(2)-p3k0*cw78%e(2)
      l3_78%a(2)=ccl*(cauxa-ceps_0)
      l3_78%c(2)=ccl*(-cauxc+ceps_1)
      endif
  
      if (ilept(id4).ne.1.or.ilept(id7).ne.1) then
      quqd=(-s456+s56)/2d0
* TWR0 -- qu=p456,qd=p4,v=cw56%e,a=r4_56%a,b=r4_56%b,cl=wcl,nsum=0
      ceps_0=-cw56%ek0*(p456(2)*p4(3)-p4(2)*p456(3))+p456k0*(cw5
     & 6%e(2)*p4(3)-p4(2)*cw56%e(3))-p4k0*(cw56%e(2)*p456(3)-p45
     & 6(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw56%e(3)*p4k0+p4(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p456(0)-cw56%e(1)*p456(1)-cw56%e(2)*p456(2)
     & -cw56%e(3)*p456(3)
      cvqd=cw56%e(0)*p4(0)-cw56%e(1)*p4(1)-cw56%e(2)*p4(2)-cw56%
     & e(3)*p4(3)
      cauxa=-cw56%ek0*quqd+p456k0*cvqd+p4k0*cvqu
      cauxb=-cw56%ek0*p4(2)+p4k0*cw56%e(2)
      r4_56%a(2)=wcl*(cauxa-ceps_0)
      r4_56%b(1)=wcl*(cauxb-ceps_2)
      endif
  
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p3,q=p456
      quqd=p3(0)*p456(0)-p3(1)*p456(1)-p3(2)*p456(2)-p3(3)*p456(
     & 3)
      ccl=wcl/(f456)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p3,qd=p456,v=cw7128(i1,i2)%e,a=l3_7128(i1,i2)%a,c=l3_7128(i
* 1,i2)%c,cl=ccl,nsum=0
      ceps_0=-cw7128(i1,i2)%ek0*(p3(2)*p456(3)-p456(2)*p3(3))+p3
     & k0*(cw7128(i1,i2)%e(2)*p456(3)-p456(2)*cw7128(i1,i2)%e(3)
     & )-p456k0*(cw7128(i1,i2)%e(2)*p3(3)-p3(2)*cw7128(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cw7128(i1,i2)%e(3)*p3k0+p3(3)*cw7128(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw7128(i1,i2)%e(0)*p3(0)-cw7128(i1,i2)%e(1)*p3(1)-cw7
     & 128(i1,i2)%e(2)*p3(2)-cw7128(i1,i2)%e(3)*p3(3)
      cvqd=cw7128(i1,i2)%e(0)*p456(0)-cw7128(i1,i2)%e(1)*p456(1)
     & -cw7128(i1,i2)%e(2)*p456(2)-cw7128(i1,i2)%e(3)*p456(3)
      cauxa=-cw7128(i1,i2)%ek0*quqd+p3k0*cvqd+p456k0*cvqu
      cauxc=+cw7128(i1,i2)%ek0*p3(2)-p3k0*cw7128(i1,i2)%e(2)
      l3_7128(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l3_7128(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      if (ilept(id3).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p3,q=p3178
      quqd=p3(0)*p3178(0)-p3(1)*p3178(1)-p3(2)*p3178(2)-p3(3)*p3
     & 178(3)
      ccl=wcl/(f3178)
      do i1=1,2
* TWL0 -- qu=p3,qd=p3178,v=cw718(i1)%e,a=l3_718(i1)%a,c=l3_718(i1)%c,cl=
* ccl,nsum=0
      ceps_0=-cw718(i1)%ek0*(p3(2)*p3178(3)-p3178(2)*p3(3))+p3k0
     & *(cw718(i1)%e(2)*p3178(3)-p3178(2)*cw718(i1)%e(3))-p3178k
     & 0*(cw718(i1)%e(2)*p3(3)-p3(2)*cw718(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw718(i1)%e(3)*p3k0+p3(3)*cw718(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw718(i1)%e(0)*p3(0)-cw718(i1)%e(1)*p3(1)-cw718(i1)%e
     & (2)*p3(2)-cw718(i1)%e(3)*p3(3)
      cvqd=cw718(i1)%e(0)*p3178(0)-cw718(i1)%e(1)*p3178(1)-cw718
     & (i1)%e(2)*p3178(2)-cw718(i1)%e(3)*p3178(3)
      cauxa=-cw718(i1)%ek0*quqd+p3k0*cvqd+p3178k0*cvqu
      cauxc=+cw718(i1)%ek0*p3(2)-p3k0*cw718(i1)%e(2)
      l3_718(i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_718(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      endif
      endif
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p378,q=p4
      quqd=p378(0)*p4(0)-p378(1)*p4(1)-p378(2)*p4(2)-p378(3)*p4(
     & 3)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p378,qd=p4,v=cw5126(i1,i2)%e,a=r4_5126(i1,i2)%a,b=r4_5126(i
* 1,i2)%b,cl=wcl,nsum=0
      ceps_0=-cw5126(i1,i2)%ek0*(p378(2)*p4(3)-p4(2)*p378(3))+p3
     & 78k0*(cw5126(i1,i2)%e(2)*p4(3)-p4(2)*cw5126(i1,i2)%e(3))-
     & p4k0*(cw5126(i1,i2)%e(2)*p378(3)-p378(2)*cw5126(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cw5126(i1,i2)%e(3)*p4k0+p4(3)*cw5126(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw5126(i1,i2)%e(0)*p378(0)-cw5126(i1,i2)%e(1)*p378(1)
     & -cw5126(i1,i2)%e(2)*p378(2)-cw5126(i1,i2)%e(3)*p378(3)
      cvqd=cw5126(i1,i2)%e(0)*p4(0)-cw5126(i1,i2)%e(1)*p4(1)-cw5
     & 126(i1,i2)%e(2)*p4(2)-cw5126(i1,i2)%e(3)*p4(3)
      cauxa=-cw5126(i1,i2)%ek0*quqd+p378k0*cvqd+p4k0*cvqu
      cauxb=-cw5126(i1,i2)%ek0*p4(2)+p4k0*cw5126(i1,i2)%e(2)
      r4_5126(i1,i2)%a(2)=wcl*(cauxa-ceps_0)
      r4_5126(i1,i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
      end do
  
      if (ilept(id4).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p3178,q=p4
      quqd=p3178(0)*p4(0)-p3178(1)*p4(1)-p3178(2)*p4(2)-p3178(3)
     & *p4(3)
      do i2=1,2
* TWR0 -- qu=p3178,qd=p4,v=cw526(i2)%e,a=r4_526(i2)%a,b=r4_526(i2)%b,cl=
* wcl,nsum=0
      ceps_0=-cw526(i2)%ek0*(p3178(2)*p4(3)-p4(2)*p3178(3))+p317
     & 8k0*(cw526(i2)%e(2)*p4(3)-p4(2)*cw526(i2)%e(3))-p4k0*(cw5
     & 26(i2)%e(2)*p3178(3)-p3178(2)*cw526(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw526(i2)%e(3)*p4k0+p4(3)*cw526(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw526(i2)%e(0)*p3178(0)-cw526(i2)%e(1)*p3178(1)-cw526
     & (i2)%e(2)*p3178(2)-cw526(i2)%e(3)*p3178(3)
      cvqd=cw526(i2)%e(0)*p4(0)-cw526(i2)%e(1)*p4(1)-cw526(i2)%e
     & (2)*p4(2)-cw526(i2)%e(3)*p4(3)
      cauxa=-cw526(i2)%ek0*quqd+p3178k0*cvqd+p4k0*cvqu
      cauxb=-cw526(i2)%ek0*p4(2)+p4k0*cw526(i2)%e(2)
      r4_526(i2)%a(2)=wcl*(cauxa-ceps_0)
      r4_526(i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
  
      endif
      endif
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p3,q=p456
      quqd=p3(0)*p456(0)-p3(1)*p456(1)-p3(2)*p456(2)-p3(3)*p456(
     & 3)
      ccl=wcl/(f456)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p3,qd=p456,v=cw7218(i1,i2)%e,a=l3_7218(i1,i2)%a,c=l3_7218(i
* 1,i2)%c,cl=ccl,nsum=0
      ceps_0=-cw7218(i1,i2)%ek0*(p3(2)*p456(3)-p456(2)*p3(3))+p3
     & k0*(cw7218(i1,i2)%e(2)*p456(3)-p456(2)*cw7218(i1,i2)%e(3)
     & )-p456k0*(cw7218(i1,i2)%e(2)*p3(3)-p3(2)*cw7218(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cw7218(i1,i2)%e(3)*p3k0+p3(3)*cw7218(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw7218(i1,i2)%e(0)*p3(0)-cw7218(i1,i2)%e(1)*p3(1)-cw7
     & 218(i1,i2)%e(2)*p3(2)-cw7218(i1,i2)%e(3)*p3(3)
      cvqd=cw7218(i1,i2)%e(0)*p456(0)-cw7218(i1,i2)%e(1)*p456(1)
     & -cw7218(i1,i2)%e(2)*p456(2)-cw7218(i1,i2)%e(3)*p456(3)
      cauxa=-cw7218(i1,i2)%ek0*quqd+p3k0*cvqd+p456k0*cvqu
      cauxc=+cw7218(i1,i2)%ek0*p3(2)-p3k0*cw7218(i1,i2)%e(2)
      l3_7218(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l3_7218(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      if (ilept(id3).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p3,q=p3278
      quqd=p3(0)*p3278(0)-p3(1)*p3278(1)-p3(2)*p3278(2)-p3(3)*p3
     & 278(3)
      ccl=wcl/(f3278)
      do i1=1,2
* TWL0 -- qu=p3,qd=p3278,v=cw728(i1)%e,a=l3_728(i1)%a,c=l3_728(i1)%c,cl=
* ccl,nsum=0
      ceps_0=-cw728(i1)%ek0*(p3(2)*p3278(3)-p3278(2)*p3(3))+p3k0
     & *(cw728(i1)%e(2)*p3278(3)-p3278(2)*cw728(i1)%e(3))-p3278k
     & 0*(cw728(i1)%e(2)*p3(3)-p3(2)*cw728(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw728(i1)%e(3)*p3k0+p3(3)*cw728(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw728(i1)%e(0)*p3(0)-cw728(i1)%e(1)*p3(1)-cw728(i1)%e
     & (2)*p3(2)-cw728(i1)%e(3)*p3(3)
      cvqd=cw728(i1)%e(0)*p3278(0)-cw728(i1)%e(1)*p3278(1)-cw728
     & (i1)%e(2)*p3278(2)-cw728(i1)%e(3)*p3278(3)
      cauxa=-cw728(i1)%ek0*quqd+p3k0*cvqd+p3278k0*cvqu
      cauxc=+cw728(i1)%ek0*p3(2)-p3k0*cw728(i1)%e(2)
      l3_728(i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_728(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      endif
      endif
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p378,q=p4
      quqd=p378(0)*p4(0)-p378(1)*p4(1)-p378(2)*p4(2)-p378(3)*p4(
     & 3)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p378,qd=p4,v=cw5216(i1,i2)%e,a=r4_5216(i1,i2)%a,b=r4_5216(i
* 1,i2)%b,cl=wcl,nsum=0
      ceps_0=-cw5216(i1,i2)%ek0*(p378(2)*p4(3)-p4(2)*p378(3))+p3
     & 78k0*(cw5216(i1,i2)%e(2)*p4(3)-p4(2)*cw5216(i1,i2)%e(3))-
     & p4k0*(cw5216(i1,i2)%e(2)*p378(3)-p378(2)*cw5216(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cw5216(i1,i2)%e(3)*p4k0+p4(3)*cw5216(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw5216(i1,i2)%e(0)*p378(0)-cw5216(i1,i2)%e(1)*p378(1)
     & -cw5216(i1,i2)%e(2)*p378(2)-cw5216(i1,i2)%e(3)*p378(3)
      cvqd=cw5216(i1,i2)%e(0)*p4(0)-cw5216(i1,i2)%e(1)*p4(1)-cw5
     & 216(i1,i2)%e(2)*p4(2)-cw5216(i1,i2)%e(3)*p4(3)
      cauxa=-cw5216(i1,i2)%ek0*quqd+p378k0*cvqd+p4k0*cvqu
      cauxb=-cw5216(i1,i2)%ek0*p4(2)+p4k0*cw5216(i1,i2)%e(2)
      r4_5216(i1,i2)%a(2)=wcl*(cauxa-ceps_0)
      r4_5216(i1,i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
      end do
  
      if (ilept(id4).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p3278,q=p4
      quqd=p3278(0)*p4(0)-p3278(1)*p4(1)-p3278(2)*p4(2)-p3278(3)
     & *p4(3)
      do i2=1,2
* TWR0 -- qu=p3278,qd=p4,v=cw516(i2)%e,a=r4_516(i2)%a,b=r4_516(i2)%b,cl=
* wcl,nsum=0
      ceps_0=-cw516(i2)%ek0*(p3278(2)*p4(3)-p4(2)*p3278(3))+p327
     & 8k0*(cw516(i2)%e(2)*p4(3)-p4(2)*cw516(i2)%e(3))-p4k0*(cw5
     & 16(i2)%e(2)*p3278(3)-p3278(2)*cw516(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw516(i2)%e(3)*p4k0+p4(3)*cw516(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw516(i2)%e(0)*p3278(0)-cw516(i2)%e(1)*p3278(1)-cw516
     & (i2)%e(2)*p3278(2)-cw516(i2)%e(3)*p3278(3)
      cvqd=cw516(i2)%e(0)*p4(0)-cw516(i2)%e(1)*p4(1)-cw516(i2)%e
     & (2)*p4(2)-cw516(i2)%e(3)*p4(3)
      cauxa=-cw516(i2)%ek0*quqd+p3278k0*cvqd+p4k0*cvqu
      cauxb=-cw516(i2)%ek0*p4(2)+p4k0*cw516(i2)%e(2)
      r4_516(i2)%a(2)=wcl*(cauxa-ceps_0)
      r4_516(i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
  
      endif
      endif
  
      else
  
      if (ilept(id3).ne.1.or.ilept(id7).ne.1) then
      quqd=(s356-s56)/2d0
      ccl=wcl/(f356)
* TWL0 -- qu=p3,qd=p356,v=cw56%e,a=l3_56%a,c=l3_56%c,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p3(2)*p356(3)-p356(2)*p3(3))+p3k0*(cw56%
     & e(2)*p356(3)-p356(2)*cw56%e(3))-p356k0*(cw56%e(2)*p3(3)-p
     & 3(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p3k0+p3(3)*cw56%ek0
      ceps_1=ceps_1*cim
      cvqu=cw56%e(0)*p3(0)-cw56%e(1)*p3(1)-cw56%e(2)*p3(2)-cw56%
     & e(3)*p3(3)
      cvqd=cw56%e(0)*p356(0)-cw56%e(1)*p356(1)-cw56%e(2)*p356(2)
     & -cw56%e(3)*p356(3)
      cauxa=-cw56%ek0*quqd+p3k0*cvqd+p356k0*cvqu
      cauxc=+cw56%ek0*p3(2)-p3k0*cw56%e(2)
      l3_56%a(2)=ccl*(cauxa-ceps_0)
      l3_56%c(2)=ccl*(-cauxc+ceps_1)
  
      endif
  
      if (ilept(id4).ne.1.or.ilept(id5).ne.1) then
      quqd=(-s478+s78)/2d0
* TWR0 -- qu=p478,qd=p4,v=cw78%e,a=r4_78%a,b=r4_78%b,cl=wcl,nsum=0
      ceps_0=-cw78%ek0*(p478(2)*p4(3)-p4(2)*p478(3))+p478k0*(cw7
     & 8%e(2)*p4(3)-p4(2)*cw78%e(3))-p4k0*(cw78%e(2)*p478(3)-p47
     & 8(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw78%e(3)*p4k0+p4(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p478(0)-cw78%e(1)*p478(1)-cw78%e(2)*p478(2)
     & -cw78%e(3)*p478(3)
      cvqd=cw78%e(0)*p4(0)-cw78%e(1)*p4(1)-cw78%e(2)*p4(2)-cw78%
     & e(3)*p4(3)
      cauxa=-cw78%ek0*quqd+p478k0*cvqd+p4k0*cvqu
      cauxb=-cw78%ek0*p4(2)+p4k0*cw78%e(2)
      r4_78%a(2)=wcl*(cauxa-ceps_0)
      r4_78%b(1)=wcl*(cauxb-ceps_2)
  
      endif
  
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p3,q=p478
      quqd=p3(0)*p478(0)-p3(1)*p478(1)-p3(2)*p478(2)-p3(3)*p478(
     & 3)
      ccl=wcl/(f478)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p3,qd=p478,v=cw5126(i1,i2)%e,a=l3_5126(i1,i2)%a,c=l3_5126(i
* 1,i2)%c,cl=ccl,nsum=0
      ceps_0=-cw5126(i1,i2)%ek0*(p3(2)*p478(3)-p478(2)*p3(3))+p3
     & k0*(cw5126(i1,i2)%e(2)*p478(3)-p478(2)*cw5126(i1,i2)%e(3)
     & )-p478k0*(cw5126(i1,i2)%e(2)*p3(3)-p3(2)*cw5126(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cw5126(i1,i2)%e(3)*p3k0+p3(3)*cw5126(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw5126(i1,i2)%e(0)*p3(0)-cw5126(i1,i2)%e(1)*p3(1)-cw5
     & 126(i1,i2)%e(2)*p3(2)-cw5126(i1,i2)%e(3)*p3(3)
      cvqd=cw5126(i1,i2)%e(0)*p478(0)-cw5126(i1,i2)%e(1)*p478(1)
     & -cw5126(i1,i2)%e(2)*p478(2)-cw5126(i1,i2)%e(3)*p478(3)
      cauxa=-cw5126(i1,i2)%ek0*quqd+p3k0*cvqd+p478k0*cvqu
      cauxc=+cw5126(i1,i2)%ek0*p3(2)-p3k0*cw5126(i1,i2)%e(2)
      l3_5126(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l3_5126(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      if (ilept(id3).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p3,q=p3156
      quqd=p3(0)*p3156(0)-p3(1)*p3156(1)-p3(2)*p3156(2)-p3(3)*p3
     & 156(3)
      ccl=wcl/(f3156)
      do i1=1,2
* TWL0 -- qu=p3,qd=p3156,v=cw516(i1)%e,a=l3_516(i1)%a,c=l3_516(i1)%c,cl=
* ccl,nsum=0
      ceps_0=-cw516(i1)%ek0*(p3(2)*p3156(3)-p3156(2)*p3(3))+p3k0
     & *(cw516(i1)%e(2)*p3156(3)-p3156(2)*cw516(i1)%e(3))-p3156k
     & 0*(cw516(i1)%e(2)*p3(3)-p3(2)*cw516(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw516(i1)%e(3)*p3k0+p3(3)*cw516(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw516(i1)%e(0)*p3(0)-cw516(i1)%e(1)*p3(1)-cw516(i1)%e
     & (2)*p3(2)-cw516(i1)%e(3)*p3(3)
      cvqd=cw516(i1)%e(0)*p3156(0)-cw516(i1)%e(1)*p3156(1)-cw516
     & (i1)%e(2)*p3156(2)-cw516(i1)%e(3)*p3156(3)
      cauxa=-cw516(i1)%ek0*quqd+p3k0*cvqd+p3156k0*cvqu
      cauxc=+cw516(i1)%ek0*p3(2)-p3k0*cw516(i1)%e(2)
      l3_516(i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_516(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      endif
      endif
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p356,q=p4
      quqd=p356(0)*p4(0)-p356(1)*p4(1)-p356(2)*p4(2)-p356(3)*p4(
     & 3)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p356,qd=p4,v=cw7128(i1,i2)%e,a=r4_7128(i1,i2)%a,b=r4_7128(i
* 1,i2)%b,cl=wcl,nsum=0
      ceps_0=-cw7128(i1,i2)%ek0*(p356(2)*p4(3)-p4(2)*p356(3))+p3
     & 56k0*(cw7128(i1,i2)%e(2)*p4(3)-p4(2)*cw7128(i1,i2)%e(3))-
     & p4k0*(cw7128(i1,i2)%e(2)*p356(3)-p356(2)*cw7128(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cw7128(i1,i2)%e(3)*p4k0+p4(3)*cw7128(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw7128(i1,i2)%e(0)*p356(0)-cw7128(i1,i2)%e(1)*p356(1)
     & -cw7128(i1,i2)%e(2)*p356(2)-cw7128(i1,i2)%e(3)*p356(3)
      cvqd=cw7128(i1,i2)%e(0)*p4(0)-cw7128(i1,i2)%e(1)*p4(1)-cw7
     & 128(i1,i2)%e(2)*p4(2)-cw7128(i1,i2)%e(3)*p4(3)
      cauxa=-cw7128(i1,i2)%ek0*quqd+p356k0*cvqd+p4k0*cvqu
      cauxb=-cw7128(i1,i2)%ek0*p4(2)+p4k0*cw7128(i1,i2)%e(2)
      r4_7128(i1,i2)%a(2)=wcl*(cauxa-ceps_0)
      r4_7128(i1,i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
      end do
  
      if (ilept(id4).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p3156,q=p4
      quqd=p3156(0)*p4(0)-p3156(1)*p4(1)-p3156(2)*p4(2)-p3156(3)
     & *p4(3)
      do i2=1,2
* TWR0 -- qu=p3156,qd=p4,v=cw728(i2)%e,a=r4_728(i2)%a,b=r4_728(i2)%b,cl=
* wcl,nsum=0
      ceps_0=-cw728(i2)%ek0*(p3156(2)*p4(3)-p4(2)*p3156(3))+p315
     & 6k0*(cw728(i2)%e(2)*p4(3)-p4(2)*cw728(i2)%e(3))-p4k0*(cw7
     & 28(i2)%e(2)*p3156(3)-p3156(2)*cw728(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw728(i2)%e(3)*p4k0+p4(3)*cw728(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw728(i2)%e(0)*p3156(0)-cw728(i2)%e(1)*p3156(1)-cw728
     & (i2)%e(2)*p3156(2)-cw728(i2)%e(3)*p3156(3)
      cvqd=cw728(i2)%e(0)*p4(0)-cw728(i2)%e(1)*p4(1)-cw728(i2)%e
     & (2)*p4(2)-cw728(i2)%e(3)*p4(3)
      cauxa=-cw728(i2)%ek0*quqd+p3156k0*cvqd+p4k0*cvqu
      cauxb=-cw728(i2)%ek0*p4(2)+p4k0*cw728(i2)%e(2)
      r4_728(i2)%a(2)=wcl*(cauxa-ceps_0)
      r4_728(i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
  
      endif
      endif
  
      if (ilept(id5).ne.1) then
  
* quqd -- p=p3,q=p478
      quqd=p3(0)*p478(0)-p3(1)*p478(1)-p3(2)*p478(2)-p3(3)*p478(
     & 3)
      ccl=wcl/(f478)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p3,qd=p478,v=cw5216(i1,i2)%e,a=l3_5216(i1,i2)%a,c=l3_5216(i
* 1,i2)%c,cl=ccl,nsum=0
      ceps_0=-cw5216(i1,i2)%ek0*(p3(2)*p478(3)-p478(2)*p3(3))+p3
     & k0*(cw5216(i1,i2)%e(2)*p478(3)-p478(2)*cw5216(i1,i2)%e(3)
     & )-p478k0*(cw5216(i1,i2)%e(2)*p3(3)-p3(2)*cw5216(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cw5216(i1,i2)%e(3)*p3k0+p3(3)*cw5216(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw5216(i1,i2)%e(0)*p3(0)-cw5216(i1,i2)%e(1)*p3(1)-cw5
     & 216(i1,i2)%e(2)*p3(2)-cw5216(i1,i2)%e(3)*p3(3)
      cvqd=cw5216(i1,i2)%e(0)*p478(0)-cw5216(i1,i2)%e(1)*p478(1)
     & -cw5216(i1,i2)%e(2)*p478(2)-cw5216(i1,i2)%e(3)*p478(3)
      cauxa=-cw5216(i1,i2)%ek0*quqd+p3k0*cvqd+p478k0*cvqu
      cauxc=+cw5216(i1,i2)%ek0*p3(2)-p3k0*cw5216(i1,i2)%e(2)
      l3_5216(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l3_5216(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
  
      if (ilept(id3).ne.1.or.ilept(id7).ne.1) then
  
* quqd -- p=p3,q=p3256
      quqd=p3(0)*p3256(0)-p3(1)*p3256(1)-p3(2)*p3256(2)-p3(3)*p3
     & 256(3)
      ccl=wcl/(f3256)
      do i1=1,2
* TWL0 -- qu=p3,qd=p3256,v=cw526(i1)%e,a=l3_526(i1)%a,c=l3_526(i1)%c,cl=
* ccl,nsum=0
      ceps_0=-cw526(i1)%ek0*(p3(2)*p3256(3)-p3256(2)*p3(3))+p3k0
     & *(cw526(i1)%e(2)*p3256(3)-p3256(2)*cw526(i1)%e(3))-p3256k
     & 0*(cw526(i1)%e(2)*p3(3)-p3(2)*cw526(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw526(i1)%e(3)*p3k0+p3(3)*cw526(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cw526(i1)%e(0)*p3(0)-cw526(i1)%e(1)*p3(1)-cw526(i1)%e
     & (2)*p3(2)-cw526(i1)%e(3)*p3(3)
      cvqd=cw526(i1)%e(0)*p3256(0)-cw526(i1)%e(1)*p3256(1)-cw526
     & (i1)%e(2)*p3256(2)-cw526(i1)%e(3)*p3256(3)
      cauxa=-cw526(i1)%ek0*quqd+p3k0*cvqd+p3256k0*cvqu
      cauxc=+cw526(i1)%ek0*p3(2)-p3k0*cw526(i1)%e(2)
      l3_526(i1)%a(2)=ccl*(cauxa-ceps_0)
      l3_526(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
      endif
      endif
  
      if (ilept(id7).ne.1) then
  
* quqd -- p=p356,q=p4
      quqd=p356(0)*p4(0)-p356(1)*p4(1)-p356(2)*p4(2)-p356(3)*p4(
     & 3)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p356,qd=p4,v=cw7218(i1,i2)%e,a=r4_7218(i1,i2)%a,b=r4_7218(i
* 1,i2)%b,cl=wcl,nsum=0
      ceps_0=-cw7218(i1,i2)%ek0*(p356(2)*p4(3)-p4(2)*p356(3))+p3
     & 56k0*(cw7218(i1,i2)%e(2)*p4(3)-p4(2)*cw7218(i1,i2)%e(3))-
     & p4k0*(cw7218(i1,i2)%e(2)*p356(3)-p356(2)*cw7218(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_2=-cw7218(i1,i2)%e(3)*p4k0+p4(3)*cw7218(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw7218(i1,i2)%e(0)*p356(0)-cw7218(i1,i2)%e(1)*p356(1)
     & -cw7218(i1,i2)%e(2)*p356(2)-cw7218(i1,i2)%e(3)*p356(3)
      cvqd=cw7218(i1,i2)%e(0)*p4(0)-cw7218(i1,i2)%e(1)*p4(1)-cw7
     & 218(i1,i2)%e(2)*p4(2)-cw7218(i1,i2)%e(3)*p4(3)
      cauxa=-cw7218(i1,i2)%ek0*quqd+p356k0*cvqd+p4k0*cvqu
      cauxb=-cw7218(i1,i2)%ek0*p4(2)+p4k0*cw7218(i1,i2)%e(2)
      r4_7218(i1,i2)%a(2)=wcl*(cauxa-ceps_0)
      r4_7218(i1,i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
      end do
  
      if (ilept(id4).ne.1.or.ilept(id5).ne.1) then
  
* quqd -- p=p3256,q=p4
      quqd=p3256(0)*p4(0)-p3256(1)*p4(1)-p3256(2)*p4(2)-p3256(3)
     & *p4(3)
      do i2=1,2
* TWR0 -- qu=p3256,qd=p4,v=cw718(i2)%e,a=r4_718(i2)%a,b=r4_718(i2)%b,cl=
* wcl,nsum=0
      ceps_0=-cw718(i2)%ek0*(p3256(2)*p4(3)-p4(2)*p3256(3))+p325
     & 6k0*(cw718(i2)%e(2)*p4(3)-p4(2)*cw718(i2)%e(3))-p4k0*(cw7
     & 18(i2)%e(2)*p3256(3)-p3256(2)*cw718(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw718(i2)%e(3)*p4k0+p4(3)*cw718(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw718(i2)%e(0)*p3256(0)-cw718(i2)%e(1)*p3256(1)-cw718
     & (i2)%e(2)*p3256(2)-cw718(i2)%e(3)*p3256(3)
      cvqd=cw718(i2)%e(0)*p4(0)-cw718(i2)%e(1)*p4(1)-cw718(i2)%e
     & (2)*p4(2)-cw718(i2)%e(3)*p4(3)
      cauxa=-cw718(i2)%ek0*quqd+p3256k0*cvqd+p4k0*cvqu
      cauxb=-cw718(i2)%ek0*p4(2)+p4k0*cw718(i2)%e(2)
      r4_718(i2)%a(2)=wcl*(cauxa-ceps_0)
      r4_718(i2)%b(1)=wcl*(cauxb-ceps_2)
      end do
  
      endif
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
* W                                                                     
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,1),a1=l7_314(i3,i1)%a,c1=l7_314(i3,i1)%c,a
* 2=r8_526(i2)%a,b2=r8_526(i2)%b,prq=s7134,bef=,aft=
      cres(i1,i2,i3,1)=(l7_314(i3,i1)%c(2)*s7134*r8_526(i2)%b(1)
     & +l7_314(i3,i1)%a(2)*r8_526(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,1),a1=l7_526(i2)%a,c1=l7_526(i2)%c,a2=r8_3
* 14(i3,i1)%a,b2=r8_314(i3,i1)%b,prq=s7256,bef=cres(i1,i2,i3,1)+,aft=
      cres(i1,i2,i3,1)=cres(i1,i2,i3,1)+(l7_526(i2)%c(2)*s7256*r
     & 8_314(i3,i1)%b(1)+l7_526(i2)%a(2)*r8_314(i3,i1)%a(2))
      end do
      end do
      end do
      endif
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,2),a1=l7_324(i3,i2)%a,c1=l7_324(i3,i2)%c,a
* 2=r8_516(i1)%a,b2=r8_516(i1)%b,prq=s7234,bef=,aft=
      cres(i1,i2,i3,2)=(l7_324(i3,i2)%c(2)*s7234*r8_516(i1)%b(1)
     & +l7_324(i3,i2)%a(2)*r8_516(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,2),a1=l7_516(i1)%a,c1=l7_516(i1)%c,a2=r8_3
* 24(i3,i2)%a,b2=r8_324(i3,i2)%b,prq=s7156,bef=cres(i1,i2,i3,2)+,aft=
      cres(i1,i2,i3,2)=cres(i1,i2,i3,2)+(l7_516(i1)%c(2)*s7156*r
     & 8_324(i3,i2)%b(1)+l7_516(i1)%a(2)*r8_324(i3,i2)%a(2))
      end do
      end do
      end do
      endif
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,3),a1=l5_314(i3,i1)%a,c1=l5_314(i3,i1)%c,a
* 2=r6_728(i2)%a,b2=r6_728(i2)%b,prq=s5134,bef=,aft=
      cres(i1,i2,i3,3)=(l5_314(i3,i1)%c(2)*s5134*r6_728(i2)%b(1)
     & +l5_314(i3,i1)%a(2)*r6_728(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,3),a1=l5_728(i2)%a,c1=l5_728(i2)%c,a2=r6_3
* 14(i3,i1)%a,b2=r6_314(i3,i1)%b,prq=s5278,bef=cres(i1,i2,i3,3)+,aft=
      cres(i1,i2,i3,3)=cres(i1,i2,i3,3)+(l5_728(i2)%c(2)*s5278*r
     & 6_314(i3,i1)%b(1)+l5_728(i2)%a(2)*r6_314(i3,i1)%a(2))
      end do
      end do
      end do
      endif
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,4),a1=l5_324(i3,i2)%a,c1=l5_324(i3,i2)%c,a
* 2=r6_718(i1)%a,b2=r6_718(i1)%b,prq=s5234,bef=,aft=
      cres(i1,i2,i3,4)=(l5_324(i3,i2)%c(2)*s5234*r6_718(i1)%b(1)
     & +l5_324(i3,i2)%a(2)*r6_718(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,4),a1=l5_718(i1)%a,c1=l5_718(i1)%c,a2=r6_3
* 24(i3,i2)%a,b2=r6_324(i3,i2)%b,prq=s5178,bef=cres(i1,i2,i3,4)+,aft=
      cres(i1,i2,i3,4)=cres(i1,i2,i3,4)+(l5_718(i1)%c(2)*s5178*r
     & 6_324(i3,i2)%b(1)+l5_718(i1)%a(2)*r6_324(i3,i2)%a(2))
      end do
      end do
      end do
      endif
  
* Z                                                                     
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
        if (iup(id3).eq.1) then
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,5),a1=l3_728(i2)%a,c1=l3_728(i2)%c,a2=r4_51
* 6(i1)%a,b2=r4_516(i1)%b,prq=s3278,bef=,aft=
      cres(i1,i2,2,5)=(l3_728(i2)%c(2)*s3278*r4_516(i1)%b(1)+l3_
     & 728(i2)%a(2)*r4_516(i1)%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,5),a1=l3_516(i1)%a,c1=l3_516(i1)%c,a2=r4_72
* 8(i2)%a,b2=r4_728(i2)%b,prq=s3156,bef=,aft=
      cres(i1,i2,2,5)=(l3_516(i1)%c(2)*s3156*r4_728(i2)%b(1)+l3_
     & 516(i1)%a(2)*r4_728(i2)%a(2))
      end do
      end do
        endif
      endif
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
        if (iup(id3).eq.1) then
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,6),a1=l3_718(i1)%a,c1=l3_718(i1)%c,a2=r4_52
* 6(i2)%a,b2=r4_526(i2)%b,prq=s3178,bef=,aft=
      cres(i1,i2,2,6)=(l3_718(i1)%c(2)*s3178*r4_526(i2)%b(1)+l3_
     & 718(i1)%a(2)*r4_526(i2)%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,6),a1=l3_526(i2)%a,c1=l3_526(i2)%c,a2=r4_71
* 8(i1)%a,b2=r4_718(i1)%b,prq=s3256,bef=,aft=
      cres(i1,i2,2,6)=(l3_526(i2)%c(2)*s3256*r4_718(i1)%b(1)+l3_
     & 526(i2)%a(2)*r4_718(i1)%a(2))
      end do
      end do
        endif
      endif
  
      do i1=1,2
      do i2=1,2
      do icol=5,6
        cres(i1,i2,1,icol)=czero
      enddo
      enddo
      enddo
  
  
* line without gluon; col 7-12                                          
  
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,7),a1=l7_3124(i3,i1,i2)%a,c1=l7_3124(i3,i1
* ,i2)%c,a2=r8_56%a,b2=r8_56%b,prq=s856,bef=,aft=
      cres(i1,i2,i3,7)=(l7_3124(i3,i1,i2)%c(2)*s856*r8_56%b(1)+l
     & 7_3124(i3,i1,i2)%a(2)*r8_56%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,7),a1=l7_56%a,c1=l7_56%c,a2=r8_3124(i3,i1,
* i2)%a,b2=r8_3124(i3,i1,i2)%b,prq=s756,bef=cres(i1,i2,i3,7)+,aft=
      cres(i1,i2,i3,7)=cres(i1,i2,i3,7)+(l7_56%c(2)*s756*r8_3124
     & (i3,i1,i2)%b(1)+l7_56%a(2)*r8_3124(i3,i1,i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,7),a1=l5_3124(i3,i1,i2)%a,c1=l5_3124(i3,i1
* ,i2)%c,a2=r6_78%a,b2=r6_78%b,prq=s678,bef=cres(i1,i2,i3,7)+,aft=
      cres(i1,i2,i3,7)=cres(i1,i2,i3,7)+(l5_3124(i3,i1,i2)%c(2)*
     & s678*r6_78%b(1)+l5_3124(i3,i1,i2)%a(2)*r6_78%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,7),a1=l5_78%a,c1=l5_78%c,a2=r6_3124(i3,i1,
* i2)%a,b2=r6_3124(i3,i1,i2)%b,prq=s578,bef=cres(i1,i2,i3,7)+,aft=
      cres(i1,i2,i3,7)=cres(i1,i2,i3,7)+(l5_78%c(2)*s578*r6_3124
     & (i3,i1,i2)%b(1)+l5_78%a(2)*r6_3124(i3,i1,i2)%a(2))
      end do
      end do
      end do
      endif
  
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,8),a1=l7_3214(i3,i1,i2)%a,c1=l7_3214(i3,i1
* ,i2)%c,a2=r8_56%a,b2=r8_56%b,prq=s856,bef=,aft=
      cres(i1,i2,i3,8)=(l7_3214(i3,i1,i2)%c(2)*s856*r8_56%b(1)+l
     & 7_3214(i3,i1,i2)%a(2)*r8_56%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,8),a1=l7_56%a,c1=l7_56%c,a2=r8_3214(i3,i1,
* i2)%a,b2=r8_3214(i3,i1,i2)%b,prq=s756,bef=cres(i1,i2,i3,8)+,aft=
      cres(i1,i2,i3,8)=cres(i1,i2,i3,8)+(l7_56%c(2)*s756*r8_3214
     & (i3,i1,i2)%b(1)+l7_56%a(2)*r8_3214(i3,i1,i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,8),a1=l5_3214(i3,i1,i2)%a,c1=l5_3214(i3,i1
* ,i2)%c,a2=r6_78%a,b2=r6_78%b,prq=s678,bef=cres(i1,i2,i3,8)+,aft=
      cres(i1,i2,i3,8)=cres(i1,i2,i3,8)+(l5_3214(i3,i1,i2)%c(2)*
     & s678*r6_78%b(1)+l5_3214(i3,i1,i2)%a(2)*r6_78%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,8),a1=l5_78%a,c1=l5_78%c,a2=r6_3214(i3,i1,
* i2)%a,b2=r6_3214(i3,i1,i2)%b,prq=s578,bef=cres(i1,i2,i3,8)+,aft=
      cres(i1,i2,i3,8)=cres(i1,i2,i3,8)+(l5_78%c(2)*s578*r6_3214
     & (i3,i1,i2)%b(1)+l5_78%a(2)*r6_3214(i3,i1,i2)%a(2))
      end do
      end do
      end do
      endif
  
  
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,9),a1=l7_5126(i1,i2)%a,c1=l7_5126(i1,i2)%c
* ,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=,aft=
      cres(i1,i2,i3,9)=(l7_5126(i1,i2)%c(2)*s834*r8_34(i3)%b(1)+
     & l7_5126(i1,i2)%a(2)*r8_34(i3)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,9),a1=l7_34(i3)%a,c1=l7_34(i3)%c,a2=r8_512
* 6(i1,i2)%a,b2=r8_5126(i1,i2)%b,prq=s734,bef=cres(i1,i2,i3,9)+,aft=
      cres(i1,i2,i3,9)=cres(i1,i2,i3,9)+(l7_34(i3)%c(2)*s734*r8_
     & 5126(i1,i2)%b(1)+l7_34(i3)%a(2)*r8_5126(i1,i2)%a(2))
      end do
      end do
      end do
  
        if (iup(id3).eq.1) then
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,9),a1=l3_78%a,c1=l3_78%c,a2=r4_5126(i1,i2)%
* a,b2=r4_5126(i1,i2)%b,prq=s378,bef=cres(i1,i2,2,9)+,aft=
      cres(i1,i2,2,9)=cres(i1,i2,2,9)+(l3_78%c(2)*s378*r4_5126(i
     & 1,i2)%b(1)+l3_78%a(2)*r4_5126(i1,i2)%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,9),a1=l3_5126(i1,i2)%a,c1=l3_5126(i1,i2)%c,
* a2=r4_78%a,b2=r4_78%b,prq=s478,bef=cres(i1,i2,2,9)+,aft=
      cres(i1,i2,2,9)=cres(i1,i2,2,9)+(l3_5126(i1,i2)%c(2)*s478*
     & r4_78%b(1)+l3_5126(i1,i2)%a(2)*r4_78%a(2))
      end do
      end do
        endif
  
      endif
  
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,10),a1=l7_5216(i1,i2)%a,c1=l7_5216(i1,i2)%
* c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=,aft=
      cres(i1,i2,i3,10)=(l7_5216(i1,i2)%c(2)*s834*r8_34(i3)%b(1)
     & +l7_5216(i1,i2)%a(2)*r8_34(i3)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,10),a1=l7_34(i3)%a,c1=l7_34(i3)%c,a2=r8_52
* 16(i1,i2)%a,b2=r8_5216(i1,i2)%b,prq=s734,bef=cres(i1,i2,i3,10)+,aft=
      cres(i1,i2,i3,10)=cres(i1,i2,i3,10)+(l7_34(i3)%c(2)*s734*r
     & 8_5216(i1,i2)%b(1)+l7_34(i3)%a(2)*r8_5216(i1,i2)%a(2))
      end do
      end do
      end do
  
        if (iup(id3).eq.1) then
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,10),a1=l3_78%a,c1=l3_78%c,a2=r4_5216(i1,i2)
* %a,b2=r4_5216(i1,i2)%b,prq=s378,bef=cres(i1,i2,2,10)+,aft=
      cres(i1,i2,2,10)=cres(i1,i2,2,10)+(l3_78%c(2)*s378*r4_5216
     & (i1,i2)%b(1)+l3_78%a(2)*r4_5216(i1,i2)%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,10),a1=l3_5216(i1,i2)%a,c1=l3_5216(i1,i2)%c
* ,a2=r4_78%a,b2=r4_78%b,prq=s478,bef=cres(i1,i2,2,10)+,aft=
      cres(i1,i2,2,10)=cres(i1,i2,2,10)+(l3_5216(i1,i2)%c(2)*s47
     & 8*r4_78%b(1)+l3_5216(i1,i2)%a(2)*r4_78%a(2))
      end do
      end do
        endif
  
      endif
  
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,11),a1=l5_7128(i1,i2)%a,c1=l5_7128(i1,i2)%
* c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=,aft=
      cres(i1,i2,i3,11)=(l5_7128(i1,i2)%c(2)*s634*r6_34(i3)%b(1)
     & +l5_7128(i1,i2)%a(2)*r6_34(i3)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,11),a1=l5_34(i3)%a,c1=l5_34(i3)%c,a2=r6_71
* 28(i1,i2)%a,b2=r6_7128(i1,i2)%b,prq=s534,bef=cres(i1,i2,i3,11)+,aft=
      cres(i1,i2,i3,11)=cres(i1,i2,i3,11)+(l5_34(i3)%c(2)*s534*r
     & 6_7128(i1,i2)%b(1)+l5_34(i3)%a(2)*r6_7128(i1,i2)%a(2))
      end do
      end do
      end do
  
        if (.not.iup(id3).eq.1) then
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,11),a1=l3_56%a,c1=l3_56%c,a2=r4_7128(i1,i2)
* %a,b2=r4_7128(i1,i2)%b,prq=s356,bef=cres(i1,i2,2,11)+,aft=
      cres(i1,i2,2,11)=cres(i1,i2,2,11)+(l3_56%c(2)*s356*r4_7128
     & (i1,i2)%b(1)+l3_56%a(2)*r4_7128(i1,i2)%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,11),a1=l3_7128(i1,i2)%a,c1=l3_7128(i1,i2)%c
* ,a2=r4_56%a,b2=r4_56%b,prq=s456,bef=cres(i1,i2,2,11)+,aft=
      cres(i1,i2,2,11)=cres(i1,i2,2,11)+(l3_7128(i1,i2)%c(2)*s45
     & 6*r4_56%b(1)+l3_7128(i1,i2)%a(2)*r4_56%a(2))
      end do
      end do
        endif
  
      endif
  
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,12),a1=l5_7218(i1,i2)%a,c1=l5_7218(i1,i2)%
* c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=,aft=
      cres(i1,i2,i3,12)=(l5_7218(i1,i2)%c(2)*s634*r6_34(i3)%b(1)
     & +l5_7218(i1,i2)%a(2)*r6_34(i3)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,12),a1=l5_34(i3)%a,c1=l5_34(i3)%c,a2=r6_72
* 18(i1,i2)%a,b2=r6_7218(i1,i2)%b,prq=s534,bef=cres(i1,i2,i3,12)+,aft=
      cres(i1,i2,i3,12)=cres(i1,i2,i3,12)+(l5_34(i3)%c(2)*s534*r
     & 6_7218(i1,i2)%b(1)+l5_34(i3)%a(2)*r6_7218(i1,i2)%a(2))
      end do
      end do
      end do
  
        if (.not.iup(id3).eq.1) then
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,12),a1=l3_56%a,c1=l3_56%c,a2=r4_7218(i1,i2)
* %a,b2=r4_7218(i1,i2)%b,prq=s356,bef=cres(i1,i2,2,12)+,aft=
      cres(i1,i2,2,12)=cres(i1,i2,2,12)+(l3_56%c(2)*s356*r4_7218
     & (i1,i2)%b(1)+l3_56%a(2)*r4_7218(i1,i2)%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,12),a1=l3_7218(i1,i2)%a,c1=l3_7218(i1,i2)%c
* ,a2=r4_56%a,b2=r4_56%b,prq=s456,bef=cres(i1,i2,2,12)+,aft=
      cres(i1,i2,2,12)=cres(i1,i2,2,12)+(l3_7218(i1,i2)%c(2)*s45
     & 6*r4_56%b(1)+l3_7218(i1,i2)%a(2)*r4_56%a(2))
      end do
      end do
        endif
  
      endif
  
  
  
* -> lines with one gluon ------------------------------                
* W 3~gluon                                                             
  
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p51,q=p5134
      quqd=p51(0)*p5134(0)-p51(1)*p5134(1)-p51(2)*p5134(2)-p51(3
     & )*p5134(3)
      ccl=zcl(id5)/(f5134)
      do i3=1,2
* TW0 -- qu=p51,qd=p5134,v=cz34(i3)%e,a=u51_34(i3)%a,b=u51_34(i3)%b,c=u5
* 1_34(i3)%c,d=u51_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p51(2)*p5134(3)-p5134(2)*p51(3))+p51
     & k0*(cz34(i3)%e(2)*p5134(3)-p5134(2)*cz34(i3)%e(3))-p5134k
     & 0*(cz34(i3)%e(2)*p51(3)-p51(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p51k0+p51(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p5134k0+p5134(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p51(0)-cz34(i3)%e(1)*p51(1)-cz34(i3)%e(
     & 2)*p51(2)-cz34(i3)%e(3)*p51(3)
      cvqd=cz34(i3)%e(0)*p5134(0)-cz34(i3)%e(1)*p5134(1)-cz34(i3
     & )%e(2)*p5134(2)-cz34(i3)%e(3)*p5134(3)
      cauxa=-cz34(i3)%ek0*quqd+p51k0*cvqd+p5134k0*cvqu
      cauxb=-cz34(i3)%ek0*p5134(2)+p5134k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p51(2)-p51k0*cz34(i3)%e(2)
      u51_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u51_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u51_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u51_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id5)/(f5134)
      do i3=1,2
* TW0 -- qu=p51,qd=p5134,v=cf34(i3)%e,a=u51_34(i3)%a,b=u51_34(i3)%b,c=u5
* 1_34(i3)%c,d=u51_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p51(2)*p5134(3)-p5134(2)*p51(3))+p51
     & k0*(cf34(i3)%e(2)*p5134(3)-p5134(2)*cf34(i3)%e(3))-p5134k
     & 0*(cf34(i3)%e(2)*p51(3)-p51(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p51k0+p51(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p5134k0+p5134(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p51(0)-cf34(i3)%e(1)*p51(1)-cf34(i3)%e(
     & 2)*p51(2)-cf34(i3)%e(3)*p51(3)
      cvqd=cf34(i3)%e(0)*p5134(0)-cf34(i3)%e(1)*p5134(1)-cf34(i3
     & )%e(2)*p5134(2)-cf34(i3)%e(3)*p5134(3)
      cauxa=-cf34(i3)%ek0*quqd+p51k0*cvqd+p5134k0*cvqu
      cauxb=-cf34(i3)%ek0*p5134(2)+p5134k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p51(2)-p51k0*cf34(i3)%e(2)
      u51_34(i3)%a(2)=u51_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u51_34(i3)%b(1)=u51_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u51_34(i3)%c(2)=u51_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u51_34(i3)%d(1)=u51_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l5_134(i1,i3)%a,cc=l5_134(i1,i3)%c,a1=l5_1(i1)%a,c1=l5_1(
* i1)%c,a2=u51_34(i3)%a,b2=u51_34(i3)%b,c2=u51_34(i3)%c,d2=u51_34(i3)%d,
* prq=s51,nsum=0
      l5_134(i1,i3)%c(2)=l5_1(i1)%c(2)*s51*u51_34(i3)%d(1)+l5_1(
     & i1)%a(2)*u51_34(i3)%c(2)
      l5_134(i1,i3)%a(2)=l5_1(i1)%c(2)*s51*u51_34(i3)%b(1)+l5_1(
     & i1)%a(2)*u51_34(i3)%a(2)
      end do
      end do
  
* quqd -- p=p534,q=p5134
      quqd=p534(0)*p5134(0)-p534(1)*p5134(1)-p534(2)*p5134(2)-p5
     & 34(3)*p5134(3)
      ccl=1.d0/(f5134)
      do i1=1,2
* TW0 -- qu=p534,qd=p5134,v=ce1(i1)%e,a=u534_1(i1)%a,b=u534_1(i1)%b,c=u5
* 34_1(i1)%c,d=u534_1(i1)%d,cl=ccl,nsum=0
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
      u534_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u534_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u534_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u534_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l5_134(i1,i3)%a,cc=l5_134(i1,i3)%c,a1=l5_34(i3)%a,c1=l5_3
* 4(i3)%c,a2=u534_1(i1)%a,b2=u534_1(i1)%b,c2=u534_1(i1)%c,d2=u534_1(i1)%
* d,prq=s534,nsum=1
      l5_134(i1,i3)%c(2)=l5_134(i1,i3)%c(2)+l5_34(i3)%c(2)*s534*
     & u534_1(i1)%d(1)+l5_34(i3)%a(2)*u534_1(i1)%c(2)
      l5_134(i1,i3)%a(2)=l5_134(i1,i3)%a(2)+l5_34(i3)%c(2)*s534*
     & u534_1(i1)%b(1)+l5_34(i3)%a(2)*u534_1(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id7).ne.1) then
* quqd -- p=p51,q=p634
      quqd=p51(0)*p634(0)-p51(1)*p634(1)-p51(2)*p634(2)-p51(3)*p
     & 634(3)
      ccl=wcl/(f634)
      do i2=1,2
* TW0 -- qu=p51,qd=p634,v=cw728(i2)%e,a=u51_728(i2)%a,b=u51_728(i2)%b,c=
* u51_728(i2)%c,d=u51_728(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw728(i2)%ek0*(p51(2)*p634(3)-p634(2)*p51(3))+p51k
     & 0*(cw728(i2)%e(2)*p634(3)-p634(2)*cw728(i2)%e(3))-p634k0*
     & (cw728(i2)%e(2)*p51(3)-p51(2)*cw728(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw728(i2)%e(3)*p51k0+p51(3)*cw728(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw728(i2)%e(3)*p634k0+p634(3)*cw728(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw728(i2)%e(0)*p51(0)-cw728(i2)%e(1)*p51(1)-cw728(i2)
     & %e(2)*p51(2)-cw728(i2)%e(3)*p51(3)
      cvqd=cw728(i2)%e(0)*p634(0)-cw728(i2)%e(1)*p634(1)-cw728(i
     & 2)%e(2)*p634(2)-cw728(i2)%e(3)*p634(3)
      cauxa=-cw728(i2)%ek0*quqd+p51k0*cvqd+p634k0*cvqu
      cauxb=-cw728(i2)%ek0*p634(2)+p634k0*cw728(i2)%e(2)
      cauxc=+cw728(i2)%ek0*p51(2)-p51k0*cw728(i2)%e(2)
      u51_728(i2)%a(2)=ccl*(cauxa-ceps_0)
      u51_728(i2)%b(1)=ccl*(cauxb-ceps_2)
      u51_728(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u51_728(i2)%d(1)=ccl*cw728(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_1728(i1,i2)%a,cc=l5_1728(i1,i2)%c,a1=l5_1(i1)%a,c1=l5_
* 1(i1)%c,a2=u51_728(i2)%a,b2=u51_728(i2)%b,c2=u51_728(i2)%c,d2=u51_728(
* i2)%d,prq=s51,nsum=0
      l5_1728(i1,i2)%c(2)=l5_1(i1)%c(2)*s51*u51_728(i2)%d(1)+l5_
     & 1(i1)%a(2)*u51_728(i2)%c(2)
      l5_1728(i1,i2)%a(2)=l5_1(i1)%c(2)*s51*u51_728(i2)%b(1)+l5_
     & 1(i1)%a(2)*u51_728(i2)%a(2)
      end do
      end do
       endif
  
* To use also as u5278_1                                                
* quqd -- p=p5278,q=p634
      quqd=p5278(0)*p634(0)-p5278(1)*p634(1)-p5278(2)*p634(2)-p5
     & 278(3)*p634(3)
      ccl=1.d0/(f634)
      do i1=1,2
* TW0 -- qu=p5278,qd=p634,v=ce1(i1)%e,a=u5728_1(i1)%a,b=u5728_1(i1)%b,c=
* u5728_1(i1)%c,d=u5728_1(i1)%d,cl=ccl,nsum=0
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
      u5728_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u5728_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u5728_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u5728_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_1728(i1,i2)%a,cc=l5_1728(i1,i2)%c,a1=l5_728(i2)%a,c1=l
* 5_728(i2)%c,a2=u5728_1(i1)%a,b2=u5728_1(i1)%b,c2=u5728_1(i1)%c,d2=u572
* 8_1(i1)%d,prq=s5278,nsum=1
      l5_1728(i1,i2)%c(2)=l5_1728(i1,i2)%c(2)+l5_728(i2)%c(2)*s5
     & 278*u5728_1(i1)%d(1)+l5_728(i2)%a(2)*u5728_1(i1)%c(2)
      l5_1728(i1,i2)%a(2)=l5_1728(i1,i2)%a(2)+l5_728(i2)%c(2)*s5
     & 278*u5728_1(i1)%b(1)+l5_728(i2)%a(2)*u5728_1(i1)%a(2)
      end do
      end do
  
* III                                                                   
* quqd -- p=p534,q=p61
      quqd=p534(0)*p61(0)-p534(1)*p61(1)-p534(2)*p61(2)-p534(3)*
     & p61(3)
      ccl=wcl/(f61)
      do i2=1,2
* TW0 -- qu=p534,qd=p61,v=cw728(i2)%e,a=u534_728(i2)%a,b=u534_728(i2)%b,
* c=u534_728(i2)%c,d=u534_728(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw728(i2)%ek0*(p534(2)*p61(3)-p61(2)*p534(3))+p534
     & k0*(cw728(i2)%e(2)*p61(3)-p61(2)*cw728(i2)%e(3))-p61k0*(c
     & w728(i2)%e(2)*p534(3)-p534(2)*cw728(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw728(i2)%e(3)*p534k0+p534(3)*cw728(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw728(i2)%e(3)*p61k0+p61(3)*cw728(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw728(i2)%e(0)*p534(0)-cw728(i2)%e(1)*p534(1)-cw728(i
     & 2)%e(2)*p534(2)-cw728(i2)%e(3)*p534(3)
      cvqd=cw728(i2)%e(0)*p61(0)-cw728(i2)%e(1)*p61(1)-cw728(i2)
     & %e(2)*p61(2)-cw728(i2)%e(3)*p61(3)
      cauxa=-cw728(i2)%ek0*quqd+p534k0*cvqd+p61k0*cvqu
      cauxb=-cw728(i2)%ek0*p61(2)+p61k0*cw728(i2)%e(2)
      cauxc=+cw728(i2)%ek0*p534(2)-p534k0*cw728(i2)%e(2)
      u534_728(i2)%a(2)=ccl*(cauxa-ceps_0)
      u534_728(i2)%b(1)=ccl*(cauxb-ceps_2)
      u534_728(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u534_728(i2)%d(1)=ccl*cw728(i2)%ek0
      end do
  
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_34728(i3,i2)%a,cc=l5_34728(i3,i2)%c,a1=l5_34(i3)%a,c1=
* l5_34(i3)%c,a2=u534_728(i2)%a,b2=u534_728(i2)%b,c2=u534_728(i2)%c,d2=u
* 534_728(i2)%d,prq=s534,nsum=0
      l5_34728(i3,i2)%c(2)=l5_34(i3)%c(2)*s534*u534_728(i2)%d(1)
     & +l5_34(i3)%a(2)*u534_728(i2)%c(2)
      l5_34728(i3,i2)%a(2)=l5_34(i3)%c(2)*s534*u534_728(i2)%b(1)
     & +l5_34(i3)%a(2)*u534_728(i2)%a(2)
      end do
      end do
       endif
  
* To use also as u5278_34                                               
* quqd -- p=p5278,q=p61
      quqd=p5278(0)*p61(0)-p5278(1)*p61(1)-p5278(2)*p61(2)-p5278
     & (3)*p61(3)
      ccl=zcl(id6)/(f61)
      do i3=1,2
* TW0 -- qu=p5278,qd=p61,v=cz34(i3)%e,a=u5728_34(i3)%a,b=u5728_34(i3)%b,
* c=u5728_34(i3)%c,d=u5728_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p5278(2)*p61(3)-p61(2)*p5278(3))+p52
     & 78k0*(cz34(i3)%e(2)*p61(3)-p61(2)*cz34(i3)%e(3))-p61k0*(c
     & z34(i3)%e(2)*p5278(3)-p5278(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p5278k0+p5278(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p61k0+p61(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p5278(0)-cz34(i3)%e(1)*p5278(1)-cz34(i3
     & )%e(2)*p5278(2)-cz34(i3)%e(3)*p5278(3)
      cvqd=cz34(i3)%e(0)*p61(0)-cz34(i3)%e(1)*p61(1)-cz34(i3)%e(
     & 2)*p61(2)-cz34(i3)%e(3)*p61(3)
      cauxa=-cz34(i3)%ek0*quqd+p5278k0*cvqd+p61k0*cvqu
      cauxb=-cz34(i3)%ek0*p61(2)+p61k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p5278(2)-p5278k0*cz34(i3)%e(2)
      u5728_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u5728_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u5728_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u5728_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id6)/(f61)
      do i3=1,2
* TW0 -- qu=p5278,qd=p61,v=cf34(i3)%e,a=u5728_34(i3)%a,b=u5728_34(i3)%b,
* c=u5728_34(i3)%c,d=u5728_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p5278(2)*p61(3)-p61(2)*p5278(3))+p52
     & 78k0*(cf34(i3)%e(2)*p61(3)-p61(2)*cf34(i3)%e(3))-p61k0*(c
     & f34(i3)%e(2)*p5278(3)-p5278(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p5278k0+p5278(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p61k0+p61(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p5278(0)-cf34(i3)%e(1)*p5278(1)-cf34(i3
     & )%e(2)*p5278(2)-cf34(i3)%e(3)*p5278(3)
      cvqd=cf34(i3)%e(0)*p61(0)-cf34(i3)%e(1)*p61(1)-cf34(i3)%e(
     & 2)*p61(2)-cf34(i3)%e(3)*p61(3)
      cauxa=-cf34(i3)%ek0*quqd+p5278k0*cvqd+p61k0*cvqu
      cauxb=-cf34(i3)%ek0*p61(2)+p61k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p5278(2)-p5278k0*cf34(i3)%e(2)
      u5728_34(i3)%a(2)=u5728_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u5728_34(i3)%b(1)=u5728_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u5728_34(i3)%c(2)=u5728_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u5728_34(i3)%d(1)=u5728_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      if (ilept(id7).ne.1) then
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_34728(i3,i2)%a,cc=l5_34728(i3,i2)%c,a1=l5_728(i2)%a,c1
* =l5_728(i2)%c,a2=u5728_34(i3)%a,b2=u5728_34(i3)%b,c2=u5728_34(i3)%c,d2
* =u5728_34(i3)%d,prq=s5278,nsum=1
      l5_34728(i3,i2)%c(2)=l5_34728(i3,i2)%c(2)+l5_728(i2)%c(2)*
     & s5278*u5728_34(i3)%d(1)+l5_728(i2)%a(2)*u5728_34(i3)%c(2)
      l5_34728(i3,i2)%a(2)=l5_34728(i3,i2)%a(2)+l5_728(i2)%c(2)*
     & s5278*u5728_34(i3)%b(1)+l5_728(i2)%a(2)*u5728_34(i3)%a(2)
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p52,q=p5234
      quqd=p52(0)*p5234(0)-p52(1)*p5234(1)-p52(2)*p5234(2)-p52(3
     & )*p5234(3)
      ccl=zcl(id5)/(f5234)
      do i3=1,2
* TW0 -- qu=p52,qd=p5234,v=cz34(i3)%e,a=u52_34(i3)%a,b=u52_34(i3)%b,c=u5
* 2_34(i3)%c,d=u52_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p52(2)*p5234(3)-p5234(2)*p52(3))+p52
     & k0*(cz34(i3)%e(2)*p5234(3)-p5234(2)*cz34(i3)%e(3))-p5234k
     & 0*(cz34(i3)%e(2)*p52(3)-p52(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p52k0+p52(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p5234k0+p5234(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p52(0)-cz34(i3)%e(1)*p52(1)-cz34(i3)%e(
     & 2)*p52(2)-cz34(i3)%e(3)*p52(3)
      cvqd=cz34(i3)%e(0)*p5234(0)-cz34(i3)%e(1)*p5234(1)-cz34(i3
     & )%e(2)*p5234(2)-cz34(i3)%e(3)*p5234(3)
      cauxa=-cz34(i3)%ek0*quqd+p52k0*cvqd+p5234k0*cvqu
      cauxb=-cz34(i3)%ek0*p5234(2)+p5234k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p52(2)-p52k0*cz34(i3)%e(2)
      u52_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u52_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u52_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u52_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id5)/(f5234)
      do i3=1,2
* TW0 -- qu=p52,qd=p5234,v=cf34(i3)%e,a=u52_34(i3)%a,b=u52_34(i3)%b,c=u5
* 2_34(i3)%c,d=u52_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p52(2)*p5234(3)-p5234(2)*p52(3))+p52
     & k0*(cf34(i3)%e(2)*p5234(3)-p5234(2)*cf34(i3)%e(3))-p5234k
     & 0*(cf34(i3)%e(2)*p52(3)-p52(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p52k0+p52(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p5234k0+p5234(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p52(0)-cf34(i3)%e(1)*p52(1)-cf34(i3)%e(
     & 2)*p52(2)-cf34(i3)%e(3)*p52(3)
      cvqd=cf34(i3)%e(0)*p5234(0)-cf34(i3)%e(1)*p5234(1)-cf34(i3
     & )%e(2)*p5234(2)-cf34(i3)%e(3)*p5234(3)
      cauxa=-cf34(i3)%ek0*quqd+p52k0*cvqd+p5234k0*cvqu
      cauxb=-cf34(i3)%ek0*p5234(2)+p5234k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p52(2)-p52k0*cf34(i3)%e(2)
      u52_34(i3)%a(2)=u52_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u52_34(i3)%b(1)=u52_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u52_34(i3)%c(2)=u52_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u52_34(i3)%d(1)=u52_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l5_234(i1,i3)%a,cc=l5_234(i1,i3)%c,a1=l5_2(i1)%a,c1=l5_2(
* i1)%c,a2=u52_34(i3)%a,b2=u52_34(i3)%b,c2=u52_34(i3)%c,d2=u52_34(i3)%d,
* prq=s52,nsum=0
      l5_234(i1,i3)%c(2)=l5_2(i1)%c(2)*s52*u52_34(i3)%d(1)+l5_2(
     & i1)%a(2)*u52_34(i3)%c(2)
      l5_234(i1,i3)%a(2)=l5_2(i1)%c(2)*s52*u52_34(i3)%b(1)+l5_2(
     & i1)%a(2)*u52_34(i3)%a(2)
      end do
      end do
  
* quqd -- p=p534,q=p5234
      quqd=p534(0)*p5234(0)-p534(1)*p5234(1)-p534(2)*p5234(2)-p5
     & 34(3)*p5234(3)
      ccl=1.d0/(f5234)
      do i1=1,2
* TW0 -- qu=p534,qd=p5234,v=ce2(i1)%e,a=u534_2(i1)%a,b=u534_2(i1)%b,c=u5
* 34_2(i1)%c,d=u534_2(i1)%d,cl=ccl,nsum=0
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
      u534_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u534_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u534_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u534_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l5_234(i1,i3)%a,cc=l5_234(i1,i3)%c,a1=l5_34(i3)%a,c1=l5_3
* 4(i3)%c,a2=u534_2(i1)%a,b2=u534_2(i1)%b,c2=u534_2(i1)%c,d2=u534_2(i1)%
* d,prq=s534,nsum=1
      l5_234(i1,i3)%c(2)=l5_234(i1,i3)%c(2)+l5_34(i3)%c(2)*s534*
     & u534_2(i1)%d(1)+l5_34(i3)%a(2)*u534_2(i1)%c(2)
      l5_234(i1,i3)%a(2)=l5_234(i1,i3)%a(2)+l5_34(i3)%c(2)*s534*
     & u534_2(i1)%b(1)+l5_34(i3)%a(2)*u534_2(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id7).ne.1) then
* quqd -- p=p52,q=p634
      quqd=p52(0)*p634(0)-p52(1)*p634(1)-p52(2)*p634(2)-p52(3)*p
     & 634(3)
      ccl=wcl/(f634)
      do i2=1,2
* TW0 -- qu=p52,qd=p634,v=cw718(i2)%e,a=u52_718(i2)%a,b=u52_718(i2)%b,c=
* u52_718(i2)%c,d=u52_718(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw718(i2)%ek0*(p52(2)*p634(3)-p634(2)*p52(3))+p52k
     & 0*(cw718(i2)%e(2)*p634(3)-p634(2)*cw718(i2)%e(3))-p634k0*
     & (cw718(i2)%e(2)*p52(3)-p52(2)*cw718(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw718(i2)%e(3)*p52k0+p52(3)*cw718(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw718(i2)%e(3)*p634k0+p634(3)*cw718(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw718(i2)%e(0)*p52(0)-cw718(i2)%e(1)*p52(1)-cw718(i2)
     & %e(2)*p52(2)-cw718(i2)%e(3)*p52(3)
      cvqd=cw718(i2)%e(0)*p634(0)-cw718(i2)%e(1)*p634(1)-cw718(i
     & 2)%e(2)*p634(2)-cw718(i2)%e(3)*p634(3)
      cauxa=-cw718(i2)%ek0*quqd+p52k0*cvqd+p634k0*cvqu
      cauxb=-cw718(i2)%ek0*p634(2)+p634k0*cw718(i2)%e(2)
      cauxc=+cw718(i2)%ek0*p52(2)-p52k0*cw718(i2)%e(2)
      u52_718(i2)%a(2)=ccl*(cauxa-ceps_0)
      u52_718(i2)%b(1)=ccl*(cauxb-ceps_2)
      u52_718(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u52_718(i2)%d(1)=ccl*cw718(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_2718(i1,i2)%a,cc=l5_2718(i1,i2)%c,a1=l5_2(i1)%a,c1=l5_
* 2(i1)%c,a2=u52_718(i2)%a,b2=u52_718(i2)%b,c2=u52_718(i2)%c,d2=u52_718(
* i2)%d,prq=s52,nsum=0
      l5_2718(i1,i2)%c(2)=l5_2(i1)%c(2)*s52*u52_718(i2)%d(1)+l5_
     & 2(i1)%a(2)*u52_718(i2)%c(2)
      l5_2718(i1,i2)%a(2)=l5_2(i1)%c(2)*s52*u52_718(i2)%b(1)+l5_
     & 2(i1)%a(2)*u52_718(i2)%a(2)
      end do
      end do
       endif
  
* To use also as u5178_2                                                
* quqd -- p=p5178,q=p634
      quqd=p5178(0)*p634(0)-p5178(1)*p634(1)-p5178(2)*p634(2)-p5
     & 178(3)*p634(3)
      ccl=1.d0/(f634)
      do i1=1,2
* TW0 -- qu=p5178,qd=p634,v=ce2(i1)%e,a=u5718_2(i1)%a,b=u5718_2(i1)%b,c=
* u5718_2(i1)%c,d=u5718_2(i1)%d,cl=ccl,nsum=0
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
      u5718_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u5718_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u5718_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u5718_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_2718(i1,i2)%a,cc=l5_2718(i1,i2)%c,a1=l5_718(i2)%a,c1=l
* 5_718(i2)%c,a2=u5718_2(i1)%a,b2=u5718_2(i1)%b,c2=u5718_2(i1)%c,d2=u571
* 8_2(i1)%d,prq=s5178,nsum=1
      l5_2718(i1,i2)%c(2)=l5_2718(i1,i2)%c(2)+l5_718(i2)%c(2)*s5
     & 178*u5718_2(i1)%d(1)+l5_718(i2)%a(2)*u5718_2(i1)%c(2)
      l5_2718(i1,i2)%a(2)=l5_2718(i1,i2)%a(2)+l5_718(i2)%c(2)*s5
     & 178*u5718_2(i1)%b(1)+l5_718(i2)%a(2)*u5718_2(i1)%a(2)
      end do
      end do
  
* III                                                                   
* quqd -- p=p534,q=p62
      quqd=p534(0)*p62(0)-p534(1)*p62(1)-p534(2)*p62(2)-p534(3)*
     & p62(3)
      ccl=wcl/(f62)
      do i2=1,2
* TW0 -- qu=p534,qd=p62,v=cw718(i2)%e,a=u534_718(i2)%a,b=u534_718(i2)%b,
* c=u534_718(i2)%c,d=u534_718(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw718(i2)%ek0*(p534(2)*p62(3)-p62(2)*p534(3))+p534
     & k0*(cw718(i2)%e(2)*p62(3)-p62(2)*cw718(i2)%e(3))-p62k0*(c
     & w718(i2)%e(2)*p534(3)-p534(2)*cw718(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw718(i2)%e(3)*p534k0+p534(3)*cw718(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw718(i2)%e(3)*p62k0+p62(3)*cw718(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw718(i2)%e(0)*p534(0)-cw718(i2)%e(1)*p534(1)-cw718(i
     & 2)%e(2)*p534(2)-cw718(i2)%e(3)*p534(3)
      cvqd=cw718(i2)%e(0)*p62(0)-cw718(i2)%e(1)*p62(1)-cw718(i2)
     & %e(2)*p62(2)-cw718(i2)%e(3)*p62(3)
      cauxa=-cw718(i2)%ek0*quqd+p534k0*cvqd+p62k0*cvqu
      cauxb=-cw718(i2)%ek0*p62(2)+p62k0*cw718(i2)%e(2)
      cauxc=+cw718(i2)%ek0*p534(2)-p534k0*cw718(i2)%e(2)
      u534_718(i2)%a(2)=ccl*(cauxa-ceps_0)
      u534_718(i2)%b(1)=ccl*(cauxb-ceps_2)
      u534_718(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u534_718(i2)%d(1)=ccl*cw718(i2)%ek0
      end do
  
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_34718(i3,i2)%a,cc=l5_34718(i3,i2)%c,a1=l5_34(i3)%a,c1=
* l5_34(i3)%c,a2=u534_718(i2)%a,b2=u534_718(i2)%b,c2=u534_718(i2)%c,d2=u
* 534_718(i2)%d,prq=s534,nsum=0
      l5_34718(i3,i2)%c(2)=l5_34(i3)%c(2)*s534*u534_718(i2)%d(1)
     & +l5_34(i3)%a(2)*u534_718(i2)%c(2)
      l5_34718(i3,i2)%a(2)=l5_34(i3)%c(2)*s534*u534_718(i2)%b(1)
     & +l5_34(i3)%a(2)*u534_718(i2)%a(2)
      end do
      end do
       endif
  
* To use also as u5178_34                                               
* quqd -- p=p5178,q=p62
      quqd=p5178(0)*p62(0)-p5178(1)*p62(1)-p5178(2)*p62(2)-p5178
     & (3)*p62(3)
      ccl=zcl(id6)/(f62)
      do i3=1,2
* TW0 -- qu=p5178,qd=p62,v=cz34(i3)%e,a=u5718_34(i3)%a,b=u5718_34(i3)%b,
* c=u5718_34(i3)%c,d=u5718_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p5178(2)*p62(3)-p62(2)*p5178(3))+p51
     & 78k0*(cz34(i3)%e(2)*p62(3)-p62(2)*cz34(i3)%e(3))-p62k0*(c
     & z34(i3)%e(2)*p5178(3)-p5178(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p5178k0+p5178(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p62k0+p62(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p5178(0)-cz34(i3)%e(1)*p5178(1)-cz34(i3
     & )%e(2)*p5178(2)-cz34(i3)%e(3)*p5178(3)
      cvqd=cz34(i3)%e(0)*p62(0)-cz34(i3)%e(1)*p62(1)-cz34(i3)%e(
     & 2)*p62(2)-cz34(i3)%e(3)*p62(3)
      cauxa=-cz34(i3)%ek0*quqd+p5178k0*cvqd+p62k0*cvqu
      cauxb=-cz34(i3)%ek0*p62(2)+p62k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p5178(2)-p5178k0*cz34(i3)%e(2)
      u5718_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u5718_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u5718_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u5718_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id6)/(f62)
      do i3=1,2
* TW0 -- qu=p5178,qd=p62,v=cf34(i3)%e,a=u5718_34(i3)%a,b=u5718_34(i3)%b,
* c=u5718_34(i3)%c,d=u5718_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p5178(2)*p62(3)-p62(2)*p5178(3))+p51
     & 78k0*(cf34(i3)%e(2)*p62(3)-p62(2)*cf34(i3)%e(3))-p62k0*(c
     & f34(i3)%e(2)*p5178(3)-p5178(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p5178k0+p5178(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p62k0+p62(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p5178(0)-cf34(i3)%e(1)*p5178(1)-cf34(i3
     & )%e(2)*p5178(2)-cf34(i3)%e(3)*p5178(3)
      cvqd=cf34(i3)%e(0)*p62(0)-cf34(i3)%e(1)*p62(1)-cf34(i3)%e(
     & 2)*p62(2)-cf34(i3)%e(3)*p62(3)
      cauxa=-cf34(i3)%ek0*quqd+p5178k0*cvqd+p62k0*cvqu
      cauxb=-cf34(i3)%ek0*p62(2)+p62k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p5178(2)-p5178k0*cf34(i3)%e(2)
      u5718_34(i3)%a(2)=u5718_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u5718_34(i3)%b(1)=u5718_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u5718_34(i3)%c(2)=u5718_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u5718_34(i3)%d(1)=u5718_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      if (ilept(id7).ne.1) then
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_34718(i3,i2)%a,cc=l5_34718(i3,i2)%c,a1=l5_718(i2)%a,c1
* =l5_718(i2)%c,a2=u5718_34(i3)%a,b2=u5718_34(i3)%b,c2=u5718_34(i3)%c,d2
* =u5718_34(i3)%d,prq=s5178,nsum=1
      l5_34718(i3,i2)%c(2)=l5_34718(i3,i2)%c(2)+l5_718(i2)%c(2)*
     & s5178*u5718_34(i3)%d(1)+l5_718(i2)%a(2)*u5718_34(i3)%c(2)
      l5_34718(i3,i2)%a(2)=l5_34718(i3,i2)%a(2)+l5_718(i2)%c(2)*
     & s5178*u5718_34(i3)%b(1)+l5_718(i2)%a(2)*u5718_34(i3)%a(2)
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p71,q=p7134
      quqd=p71(0)*p7134(0)-p71(1)*p7134(1)-p71(2)*p7134(2)-p71(3
     & )*p7134(3)
      ccl=zcl(id7)/(f7134)
      do i3=1,2
* TW0 -- qu=p71,qd=p7134,v=cz34(i3)%e,a=u71_34(i3)%a,b=u71_34(i3)%b,c=u7
* 1_34(i3)%c,d=u71_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p71(2)*p7134(3)-p7134(2)*p71(3))+p71
     & k0*(cz34(i3)%e(2)*p7134(3)-p7134(2)*cz34(i3)%e(3))-p7134k
     & 0*(cz34(i3)%e(2)*p71(3)-p71(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p71k0+p71(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p7134k0+p7134(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p71(0)-cz34(i3)%e(1)*p71(1)-cz34(i3)%e(
     & 2)*p71(2)-cz34(i3)%e(3)*p71(3)
      cvqd=cz34(i3)%e(0)*p7134(0)-cz34(i3)%e(1)*p7134(1)-cz34(i3
     & )%e(2)*p7134(2)-cz34(i3)%e(3)*p7134(3)
      cauxa=-cz34(i3)%ek0*quqd+p71k0*cvqd+p7134k0*cvqu
      cauxb=-cz34(i3)%ek0*p7134(2)+p7134k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p71(2)-p71k0*cz34(i3)%e(2)
      u71_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u71_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u71_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u71_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id7)/(f7134)
      do i3=1,2
* TW0 -- qu=p71,qd=p7134,v=cf34(i3)%e,a=u71_34(i3)%a,b=u71_34(i3)%b,c=u7
* 1_34(i3)%c,d=u71_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p71(2)*p7134(3)-p7134(2)*p71(3))+p71
     & k0*(cf34(i3)%e(2)*p7134(3)-p7134(2)*cf34(i3)%e(3))-p7134k
     & 0*(cf34(i3)%e(2)*p71(3)-p71(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p71k0+p71(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p7134k0+p7134(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p71(0)-cf34(i3)%e(1)*p71(1)-cf34(i3)%e(
     & 2)*p71(2)-cf34(i3)%e(3)*p71(3)
      cvqd=cf34(i3)%e(0)*p7134(0)-cf34(i3)%e(1)*p7134(1)-cf34(i3
     & )%e(2)*p7134(2)-cf34(i3)%e(3)*p7134(3)
      cauxa=-cf34(i3)%ek0*quqd+p71k0*cvqd+p7134k0*cvqu
      cauxb=-cf34(i3)%ek0*p7134(2)+p7134k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p71(2)-p71k0*cf34(i3)%e(2)
      u71_34(i3)%a(2)=u71_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u71_34(i3)%b(1)=u71_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u71_34(i3)%c(2)=u71_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u71_34(i3)%d(1)=u71_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l7_134(i1,i3)%a,cc=l7_134(i1,i3)%c,a1=l7_1(i1)%a,c1=l7_1(
* i1)%c,a2=u71_34(i3)%a,b2=u71_34(i3)%b,c2=u71_34(i3)%c,d2=u71_34(i3)%d,
* prq=s71,nsum=0
      l7_134(i1,i3)%c(2)=l7_1(i1)%c(2)*s71*u71_34(i3)%d(1)+l7_1(
     & i1)%a(2)*u71_34(i3)%c(2)
      l7_134(i1,i3)%a(2)=l7_1(i1)%c(2)*s71*u71_34(i3)%b(1)+l7_1(
     & i1)%a(2)*u71_34(i3)%a(2)
      end do
      end do
  
* quqd -- p=p734,q=p7134
      quqd=p734(0)*p7134(0)-p734(1)*p7134(1)-p734(2)*p7134(2)-p7
     & 34(3)*p7134(3)
      ccl=1.d0/(f7134)
      do i1=1,2
* TW0 -- qu=p734,qd=p7134,v=ce1(i1)%e,a=u734_1(i1)%a,b=u734_1(i1)%b,c=u7
* 34_1(i1)%c,d=u734_1(i1)%d,cl=ccl,nsum=0
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
      u734_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u734_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u734_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u734_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l7_134(i1,i3)%a,cc=l7_134(i1,i3)%c,a1=l7_34(i3)%a,c1=l7_3
* 4(i3)%c,a2=u734_1(i1)%a,b2=u734_1(i1)%b,c2=u734_1(i1)%c,d2=u734_1(i1)%
* d,prq=s734,nsum=1
      l7_134(i1,i3)%c(2)=l7_134(i1,i3)%c(2)+l7_34(i3)%c(2)*s734*
     & u734_1(i1)%d(1)+l7_34(i3)%a(2)*u734_1(i1)%c(2)
      l7_134(i1,i3)%a(2)=l7_134(i1,i3)%a(2)+l7_34(i3)%c(2)*s734*
     & u734_1(i1)%b(1)+l7_34(i3)%a(2)*u734_1(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id5).ne.1) then
* quqd -- p=p71,q=p834
      quqd=p71(0)*p834(0)-p71(1)*p834(1)-p71(2)*p834(2)-p71(3)*p
     & 834(3)
      ccl=wcl/(f834)
      do i2=1,2
* TW0 -- qu=p71,qd=p834,v=cw526(i2)%e,a=u71_526(i2)%a,b=u71_526(i2)%b,c=
* u71_526(i2)%c,d=u71_526(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw526(i2)%ek0*(p71(2)*p834(3)-p834(2)*p71(3))+p71k
     & 0*(cw526(i2)%e(2)*p834(3)-p834(2)*cw526(i2)%e(3))-p834k0*
     & (cw526(i2)%e(2)*p71(3)-p71(2)*cw526(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw526(i2)%e(3)*p71k0+p71(3)*cw526(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw526(i2)%e(3)*p834k0+p834(3)*cw526(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw526(i2)%e(0)*p71(0)-cw526(i2)%e(1)*p71(1)-cw526(i2)
     & %e(2)*p71(2)-cw526(i2)%e(3)*p71(3)
      cvqd=cw526(i2)%e(0)*p834(0)-cw526(i2)%e(1)*p834(1)-cw526(i
     & 2)%e(2)*p834(2)-cw526(i2)%e(3)*p834(3)
      cauxa=-cw526(i2)%ek0*quqd+p71k0*cvqd+p834k0*cvqu
      cauxb=-cw526(i2)%ek0*p834(2)+p834k0*cw526(i2)%e(2)
      cauxc=+cw526(i2)%ek0*p71(2)-p71k0*cw526(i2)%e(2)
      u71_526(i2)%a(2)=ccl*(cauxa-ceps_0)
      u71_526(i2)%b(1)=ccl*(cauxb-ceps_2)
      u71_526(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u71_526(i2)%d(1)=ccl*cw526(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_1526(i1,i2)%a,cc=l7_1526(i1,i2)%c,a1=l7_1(i1)%a,c1=l7_
* 1(i1)%c,a2=u71_526(i2)%a,b2=u71_526(i2)%b,c2=u71_526(i2)%c,d2=u71_526(
* i2)%d,prq=s71,nsum=0
      l7_1526(i1,i2)%c(2)=l7_1(i1)%c(2)*s71*u71_526(i2)%d(1)+l7_
     & 1(i1)%a(2)*u71_526(i2)%c(2)
      l7_1526(i1,i2)%a(2)=l7_1(i1)%c(2)*s71*u71_526(i2)%b(1)+l7_
     & 1(i1)%a(2)*u71_526(i2)%a(2)
      end do
      end do
       endif
  
* To use also as u7256_1                                                
* quqd -- p=p7256,q=p834
      quqd=p7256(0)*p834(0)-p7256(1)*p834(1)-p7256(2)*p834(2)-p7
     & 256(3)*p834(3)
      ccl=1.d0/(f834)
      do i1=1,2
* TW0 -- qu=p7256,qd=p834,v=ce1(i1)%e,a=u7526_1(i1)%a,b=u7526_1(i1)%b,c=
* u7526_1(i1)%c,d=u7526_1(i1)%d,cl=ccl,nsum=0
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
      u7526_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u7526_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u7526_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u7526_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_1526(i1,i2)%a,cc=l7_1526(i1,i2)%c,a1=l7_526(i2)%a,c1=l
* 7_526(i2)%c,a2=u7526_1(i1)%a,b2=u7526_1(i1)%b,c2=u7526_1(i1)%c,d2=u752
* 6_1(i1)%d,prq=s7256,nsum=1
      l7_1526(i1,i2)%c(2)=l7_1526(i1,i2)%c(2)+l7_526(i2)%c(2)*s7
     & 256*u7526_1(i1)%d(1)+l7_526(i2)%a(2)*u7526_1(i1)%c(2)
      l7_1526(i1,i2)%a(2)=l7_1526(i1,i2)%a(2)+l7_526(i2)%c(2)*s7
     & 256*u7526_1(i1)%b(1)+l7_526(i2)%a(2)*u7526_1(i1)%a(2)
      end do
      end do
  
* III                                                                   
* quqd -- p=p734,q=p81
      quqd=p734(0)*p81(0)-p734(1)*p81(1)-p734(2)*p81(2)-p734(3)*
     & p81(3)
      ccl=wcl/(f81)
      do i2=1,2
* TW0 -- qu=p734,qd=p81,v=cw526(i2)%e,a=u734_526(i2)%a,b=u734_526(i2)%b,
* c=u734_526(i2)%c,d=u734_526(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw526(i2)%ek0*(p734(2)*p81(3)-p81(2)*p734(3))+p734
     & k0*(cw526(i2)%e(2)*p81(3)-p81(2)*cw526(i2)%e(3))-p81k0*(c
     & w526(i2)%e(2)*p734(3)-p734(2)*cw526(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw526(i2)%e(3)*p734k0+p734(3)*cw526(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw526(i2)%e(3)*p81k0+p81(3)*cw526(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw526(i2)%e(0)*p734(0)-cw526(i2)%e(1)*p734(1)-cw526(i
     & 2)%e(2)*p734(2)-cw526(i2)%e(3)*p734(3)
      cvqd=cw526(i2)%e(0)*p81(0)-cw526(i2)%e(1)*p81(1)-cw526(i2)
     & %e(2)*p81(2)-cw526(i2)%e(3)*p81(3)
      cauxa=-cw526(i2)%ek0*quqd+p734k0*cvqd+p81k0*cvqu
      cauxb=-cw526(i2)%ek0*p81(2)+p81k0*cw526(i2)%e(2)
      cauxc=+cw526(i2)%ek0*p734(2)-p734k0*cw526(i2)%e(2)
      u734_526(i2)%a(2)=ccl*(cauxa-ceps_0)
      u734_526(i2)%b(1)=ccl*(cauxb-ceps_2)
      u734_526(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u734_526(i2)%d(1)=ccl*cw526(i2)%ek0
      end do
  
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_34526(i3,i2)%a,cc=l7_34526(i3,i2)%c,a1=l7_34(i3)%a,c1=
* l7_34(i3)%c,a2=u734_526(i2)%a,b2=u734_526(i2)%b,c2=u734_526(i2)%c,d2=u
* 734_526(i2)%d,prq=s734,nsum=0
      l7_34526(i3,i2)%c(2)=l7_34(i3)%c(2)*s734*u734_526(i2)%d(1)
     & +l7_34(i3)%a(2)*u734_526(i2)%c(2)
      l7_34526(i3,i2)%a(2)=l7_34(i3)%c(2)*s734*u734_526(i2)%b(1)
     & +l7_34(i3)%a(2)*u734_526(i2)%a(2)
      end do
      end do
       endif
  
* To use also as u7256_34                                               
* quqd -- p=p7256,q=p81
      quqd=p7256(0)*p81(0)-p7256(1)*p81(1)-p7256(2)*p81(2)-p7256
     & (3)*p81(3)
      ccl=zcl(id8)/(f81)
      do i3=1,2
* TW0 -- qu=p7256,qd=p81,v=cz34(i3)%e,a=u7526_34(i3)%a,b=u7526_34(i3)%b,
* c=u7526_34(i3)%c,d=u7526_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p7256(2)*p81(3)-p81(2)*p7256(3))+p72
     & 56k0*(cz34(i3)%e(2)*p81(3)-p81(2)*cz34(i3)%e(3))-p81k0*(c
     & z34(i3)%e(2)*p7256(3)-p7256(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p7256k0+p7256(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p81k0+p81(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p7256(0)-cz34(i3)%e(1)*p7256(1)-cz34(i3
     & )%e(2)*p7256(2)-cz34(i3)%e(3)*p7256(3)
      cvqd=cz34(i3)%e(0)*p81(0)-cz34(i3)%e(1)*p81(1)-cz34(i3)%e(
     & 2)*p81(2)-cz34(i3)%e(3)*p81(3)
      cauxa=-cz34(i3)%ek0*quqd+p7256k0*cvqd+p81k0*cvqu
      cauxb=-cz34(i3)%ek0*p81(2)+p81k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p7256(2)-p7256k0*cz34(i3)%e(2)
      u7526_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u7526_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u7526_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u7526_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id8)/(f81)
      do i3=1,2
* TW0 -- qu=p7256,qd=p81,v=cf34(i3)%e,a=u7526_34(i3)%a,b=u7526_34(i3)%b,
* c=u7526_34(i3)%c,d=u7526_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p7256(2)*p81(3)-p81(2)*p7256(3))+p72
     & 56k0*(cf34(i3)%e(2)*p81(3)-p81(2)*cf34(i3)%e(3))-p81k0*(c
     & f34(i3)%e(2)*p7256(3)-p7256(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p7256k0+p7256(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p81k0+p81(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p7256(0)-cf34(i3)%e(1)*p7256(1)-cf34(i3
     & )%e(2)*p7256(2)-cf34(i3)%e(3)*p7256(3)
      cvqd=cf34(i3)%e(0)*p81(0)-cf34(i3)%e(1)*p81(1)-cf34(i3)%e(
     & 2)*p81(2)-cf34(i3)%e(3)*p81(3)
      cauxa=-cf34(i3)%ek0*quqd+p7256k0*cvqd+p81k0*cvqu
      cauxb=-cf34(i3)%ek0*p81(2)+p81k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p7256(2)-p7256k0*cf34(i3)%e(2)
      u7526_34(i3)%a(2)=u7526_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u7526_34(i3)%b(1)=u7526_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u7526_34(i3)%c(2)=u7526_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u7526_34(i3)%d(1)=u7526_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      if (ilept(id5).ne.1) then
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_34526(i3,i2)%a,cc=l7_34526(i3,i2)%c,a1=l7_526(i2)%a,c1
* =l7_526(i2)%c,a2=u7526_34(i3)%a,b2=u7526_34(i3)%b,c2=u7526_34(i3)%c,d2
* =u7526_34(i3)%d,prq=s7256,nsum=1
      l7_34526(i3,i2)%c(2)=l7_34526(i3,i2)%c(2)+l7_526(i2)%c(2)*
     & s7256*u7526_34(i3)%d(1)+l7_526(i2)%a(2)*u7526_34(i3)%c(2)
      l7_34526(i3,i2)%a(2)=l7_34526(i3,i2)%a(2)+l7_526(i2)%c(2)*
     & s7256*u7526_34(i3)%b(1)+l7_526(i2)%a(2)*u7526_34(i3)%a(2)
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p72,q=p7234
      quqd=p72(0)*p7234(0)-p72(1)*p7234(1)-p72(2)*p7234(2)-p72(3
     & )*p7234(3)
      ccl=zcl(id7)/(f7234)
      do i3=1,2
* TW0 -- qu=p72,qd=p7234,v=cz34(i3)%e,a=u72_34(i3)%a,b=u72_34(i3)%b,c=u7
* 2_34(i3)%c,d=u72_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p72(2)*p7234(3)-p7234(2)*p72(3))+p72
     & k0*(cz34(i3)%e(2)*p7234(3)-p7234(2)*cz34(i3)%e(3))-p7234k
     & 0*(cz34(i3)%e(2)*p72(3)-p72(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p72k0+p72(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p7234k0+p7234(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p72(0)-cz34(i3)%e(1)*p72(1)-cz34(i3)%e(
     & 2)*p72(2)-cz34(i3)%e(3)*p72(3)
      cvqd=cz34(i3)%e(0)*p7234(0)-cz34(i3)%e(1)*p7234(1)-cz34(i3
     & )%e(2)*p7234(2)-cz34(i3)%e(3)*p7234(3)
      cauxa=-cz34(i3)%ek0*quqd+p72k0*cvqd+p7234k0*cvqu
      cauxb=-cz34(i3)%ek0*p7234(2)+p7234k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p72(2)-p72k0*cz34(i3)%e(2)
      u72_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u72_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u72_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u72_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id7)/(f7234)
      do i3=1,2
* TW0 -- qu=p72,qd=p7234,v=cf34(i3)%e,a=u72_34(i3)%a,b=u72_34(i3)%b,c=u7
* 2_34(i3)%c,d=u72_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p72(2)*p7234(3)-p7234(2)*p72(3))+p72
     & k0*(cf34(i3)%e(2)*p7234(3)-p7234(2)*cf34(i3)%e(3))-p7234k
     & 0*(cf34(i3)%e(2)*p72(3)-p72(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p72k0+p72(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p7234k0+p7234(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p72(0)-cf34(i3)%e(1)*p72(1)-cf34(i3)%e(
     & 2)*p72(2)-cf34(i3)%e(3)*p72(3)
      cvqd=cf34(i3)%e(0)*p7234(0)-cf34(i3)%e(1)*p7234(1)-cf34(i3
     & )%e(2)*p7234(2)-cf34(i3)%e(3)*p7234(3)
      cauxa=-cf34(i3)%ek0*quqd+p72k0*cvqd+p7234k0*cvqu
      cauxb=-cf34(i3)%ek0*p7234(2)+p7234k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p72(2)-p72k0*cf34(i3)%e(2)
      u72_34(i3)%a(2)=u72_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u72_34(i3)%b(1)=u72_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u72_34(i3)%c(2)=u72_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u72_34(i3)%d(1)=u72_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l7_234(i1,i3)%a,cc=l7_234(i1,i3)%c,a1=l7_2(i1)%a,c1=l7_2(
* i1)%c,a2=u72_34(i3)%a,b2=u72_34(i3)%b,c2=u72_34(i3)%c,d2=u72_34(i3)%d,
* prq=s72,nsum=0
      l7_234(i1,i3)%c(2)=l7_2(i1)%c(2)*s72*u72_34(i3)%d(1)+l7_2(
     & i1)%a(2)*u72_34(i3)%c(2)
      l7_234(i1,i3)%a(2)=l7_2(i1)%c(2)*s72*u72_34(i3)%b(1)+l7_2(
     & i1)%a(2)*u72_34(i3)%a(2)
      end do
      end do
  
* quqd -- p=p734,q=p7234
      quqd=p734(0)*p7234(0)-p734(1)*p7234(1)-p734(2)*p7234(2)-p7
     & 34(3)*p7234(3)
      ccl=1.d0/(f7234)
      do i1=1,2
* TW0 -- qu=p734,qd=p7234,v=ce2(i1)%e,a=u734_2(i1)%a,b=u734_2(i1)%b,c=u7
* 34_2(i1)%c,d=u734_2(i1)%d,cl=ccl,nsum=0
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
      u734_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u734_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u734_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u734_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l7_234(i1,i3)%a,cc=l7_234(i1,i3)%c,a1=l7_34(i3)%a,c1=l7_3
* 4(i3)%c,a2=u734_2(i1)%a,b2=u734_2(i1)%b,c2=u734_2(i1)%c,d2=u734_2(i1)%
* d,prq=s734,nsum=1
      l7_234(i1,i3)%c(2)=l7_234(i1,i3)%c(2)+l7_34(i3)%c(2)*s734*
     & u734_2(i1)%d(1)+l7_34(i3)%a(2)*u734_2(i1)%c(2)
      l7_234(i1,i3)%a(2)=l7_234(i1,i3)%a(2)+l7_34(i3)%c(2)*s734*
     & u734_2(i1)%b(1)+l7_34(i3)%a(2)*u734_2(i1)%a(2)
      end do
      end do
  
* II                                                                    
      if (ilept(id5).ne.1) then
* quqd -- p=p72,q=p834
      quqd=p72(0)*p834(0)-p72(1)*p834(1)-p72(2)*p834(2)-p72(3)*p
     & 834(3)
      ccl=wcl/(f834)
      do i2=1,2
* TW0 -- qu=p72,qd=p834,v=cw516(i2)%e,a=u72_516(i2)%a,b=u72_516(i2)%b,c=
* u72_516(i2)%c,d=u72_516(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw516(i2)%ek0*(p72(2)*p834(3)-p834(2)*p72(3))+p72k
     & 0*(cw516(i2)%e(2)*p834(3)-p834(2)*cw516(i2)%e(3))-p834k0*
     & (cw516(i2)%e(2)*p72(3)-p72(2)*cw516(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw516(i2)%e(3)*p72k0+p72(3)*cw516(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw516(i2)%e(3)*p834k0+p834(3)*cw516(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw516(i2)%e(0)*p72(0)-cw516(i2)%e(1)*p72(1)-cw516(i2)
     & %e(2)*p72(2)-cw516(i2)%e(3)*p72(3)
      cvqd=cw516(i2)%e(0)*p834(0)-cw516(i2)%e(1)*p834(1)-cw516(i
     & 2)%e(2)*p834(2)-cw516(i2)%e(3)*p834(3)
      cauxa=-cw516(i2)%ek0*quqd+p72k0*cvqd+p834k0*cvqu
      cauxb=-cw516(i2)%ek0*p834(2)+p834k0*cw516(i2)%e(2)
      cauxc=+cw516(i2)%ek0*p72(2)-p72k0*cw516(i2)%e(2)
      u72_516(i2)%a(2)=ccl*(cauxa-ceps_0)
      u72_516(i2)%b(1)=ccl*(cauxb-ceps_2)
      u72_516(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u72_516(i2)%d(1)=ccl*cw516(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_2516(i1,i2)%a,cc=l7_2516(i1,i2)%c,a1=l7_2(i1)%a,c1=l7_
* 2(i1)%c,a2=u72_516(i2)%a,b2=u72_516(i2)%b,c2=u72_516(i2)%c,d2=u72_516(
* i2)%d,prq=s72,nsum=0
      l7_2516(i1,i2)%c(2)=l7_2(i1)%c(2)*s72*u72_516(i2)%d(1)+l7_
     & 2(i1)%a(2)*u72_516(i2)%c(2)
      l7_2516(i1,i2)%a(2)=l7_2(i1)%c(2)*s72*u72_516(i2)%b(1)+l7_
     & 2(i1)%a(2)*u72_516(i2)%a(2)
      end do
      end do
       endif
  
* To use also as u7156_2                                                
* quqd -- p=p7156,q=p834
      quqd=p7156(0)*p834(0)-p7156(1)*p834(1)-p7156(2)*p834(2)-p7
     & 156(3)*p834(3)
      ccl=1.d0/(f834)
      do i1=1,2
* TW0 -- qu=p7156,qd=p834,v=ce2(i1)%e,a=u7516_2(i1)%a,b=u7516_2(i1)%b,c=
* u7516_2(i1)%c,d=u7516_2(i1)%d,cl=ccl,nsum=0
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
      u7516_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u7516_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u7516_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u7516_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_2516(i1,i2)%a,cc=l7_2516(i1,i2)%c,a1=l7_516(i2)%a,c1=l
* 7_516(i2)%c,a2=u7516_2(i1)%a,b2=u7516_2(i1)%b,c2=u7516_2(i1)%c,d2=u751
* 6_2(i1)%d,prq=s7156,nsum=1
      l7_2516(i1,i2)%c(2)=l7_2516(i1,i2)%c(2)+l7_516(i2)%c(2)*s7
     & 156*u7516_2(i1)%d(1)+l7_516(i2)%a(2)*u7516_2(i1)%c(2)
      l7_2516(i1,i2)%a(2)=l7_2516(i1,i2)%a(2)+l7_516(i2)%c(2)*s7
     & 156*u7516_2(i1)%b(1)+l7_516(i2)%a(2)*u7516_2(i1)%a(2)
      end do
      end do
  
* III                                                                   
* quqd -- p=p734,q=p82
      quqd=p734(0)*p82(0)-p734(1)*p82(1)-p734(2)*p82(2)-p734(3)*
     & p82(3)
      ccl=wcl/(f82)
      do i2=1,2
* TW0 -- qu=p734,qd=p82,v=cw516(i2)%e,a=u734_516(i2)%a,b=u734_516(i2)%b,
* c=u734_516(i2)%c,d=u734_516(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw516(i2)%ek0*(p734(2)*p82(3)-p82(2)*p734(3))+p734
     & k0*(cw516(i2)%e(2)*p82(3)-p82(2)*cw516(i2)%e(3))-p82k0*(c
     & w516(i2)%e(2)*p734(3)-p734(2)*cw516(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw516(i2)%e(3)*p734k0+p734(3)*cw516(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw516(i2)%e(3)*p82k0+p82(3)*cw516(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw516(i2)%e(0)*p734(0)-cw516(i2)%e(1)*p734(1)-cw516(i
     & 2)%e(2)*p734(2)-cw516(i2)%e(3)*p734(3)
      cvqd=cw516(i2)%e(0)*p82(0)-cw516(i2)%e(1)*p82(1)-cw516(i2)
     & %e(2)*p82(2)-cw516(i2)%e(3)*p82(3)
      cauxa=-cw516(i2)%ek0*quqd+p734k0*cvqd+p82k0*cvqu
      cauxb=-cw516(i2)%ek0*p82(2)+p82k0*cw516(i2)%e(2)
      cauxc=+cw516(i2)%ek0*p734(2)-p734k0*cw516(i2)%e(2)
      u734_516(i2)%a(2)=ccl*(cauxa-ceps_0)
      u734_516(i2)%b(1)=ccl*(cauxb-ceps_2)
      u734_516(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u734_516(i2)%d(1)=ccl*cw516(i2)%ek0
      end do
  
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_34516(i3,i2)%a,cc=l7_34516(i3,i2)%c,a1=l7_34(i3)%a,c1=
* l7_34(i3)%c,a2=u734_516(i2)%a,b2=u734_516(i2)%b,c2=u734_516(i2)%c,d2=u
* 734_516(i2)%d,prq=s734,nsum=0
      l7_34516(i3,i2)%c(2)=l7_34(i3)%c(2)*s734*u734_516(i2)%d(1)
     & +l7_34(i3)%a(2)*u734_516(i2)%c(2)
      l7_34516(i3,i2)%a(2)=l7_34(i3)%c(2)*s734*u734_516(i2)%b(1)
     & +l7_34(i3)%a(2)*u734_516(i2)%a(2)
      end do
      end do
       endif
  
* To use also as u7156_34                                               
* quqd -- p=p7156,q=p82
      quqd=p7156(0)*p82(0)-p7156(1)*p82(1)-p7156(2)*p82(2)-p7156
     & (3)*p82(3)
      ccl=zcl(id8)/(f82)
      do i3=1,2
* TW0 -- qu=p7156,qd=p82,v=cz34(i3)%e,a=u7516_34(i3)%a,b=u7516_34(i3)%b,
* c=u7516_34(i3)%c,d=u7516_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p7156(2)*p82(3)-p82(2)*p7156(3))+p71
     & 56k0*(cz34(i3)%e(2)*p82(3)-p82(2)*cz34(i3)%e(3))-p82k0*(c
     & z34(i3)%e(2)*p7156(3)-p7156(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p7156k0+p7156(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p82k0+p82(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p7156(0)-cz34(i3)%e(1)*p7156(1)-cz34(i3
     & )%e(2)*p7156(2)-cz34(i3)%e(3)*p7156(3)
      cvqd=cz34(i3)%e(0)*p82(0)-cz34(i3)%e(1)*p82(1)-cz34(i3)%e(
     & 2)*p82(2)-cz34(i3)%e(3)*p82(3)
      cauxa=-cz34(i3)%ek0*quqd+p7156k0*cvqd+p82k0*cvqu
      cauxb=-cz34(i3)%ek0*p82(2)+p82k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p7156(2)-p7156k0*cz34(i3)%e(2)
      u7516_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u7516_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u7516_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u7516_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id8)/(f82)
      do i3=1,2
* TW0 -- qu=p7156,qd=p82,v=cf34(i3)%e,a=u7516_34(i3)%a,b=u7516_34(i3)%b,
* c=u7516_34(i3)%c,d=u7516_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p7156(2)*p82(3)-p82(2)*p7156(3))+p71
     & 56k0*(cf34(i3)%e(2)*p82(3)-p82(2)*cf34(i3)%e(3))-p82k0*(c
     & f34(i3)%e(2)*p7156(3)-p7156(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p7156k0+p7156(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p82k0+p82(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p7156(0)-cf34(i3)%e(1)*p7156(1)-cf34(i3
     & )%e(2)*p7156(2)-cf34(i3)%e(3)*p7156(3)
      cvqd=cf34(i3)%e(0)*p82(0)-cf34(i3)%e(1)*p82(1)-cf34(i3)%e(
     & 2)*p82(2)-cf34(i3)%e(3)*p82(3)
      cauxa=-cf34(i3)%ek0*quqd+p7156k0*cvqd+p82k0*cvqu
      cauxb=-cf34(i3)%ek0*p82(2)+p82k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p7156(2)-p7156k0*cf34(i3)%e(2)
      u7516_34(i3)%a(2)=u7516_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u7516_34(i3)%b(1)=u7516_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u7516_34(i3)%c(2)=u7516_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u7516_34(i3)%d(1)=u7516_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      if (ilept(id5).ne.1) then
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_34516(i3,i2)%a,cc=l7_34516(i3,i2)%c,a1=l7_516(i2)%a,c1
* =l7_516(i2)%c,a2=u7516_34(i3)%a,b2=u7516_34(i3)%b,c2=u7516_34(i3)%c,d2
* =u7516_34(i3)%d,prq=s7156,nsum=1
      l7_34516(i3,i2)%c(2)=l7_34516(i3,i2)%c(2)+l7_516(i2)%c(2)*
     & s7156*u7516_34(i3)%d(1)+l7_516(i2)%a(2)*u7516_34(i3)%c(2)
      l7_34516(i3,i2)%a(2)=l7_34516(i3,i2)%a(2)+l7_516(i2)%c(2)*
     & s7156*u7516_34(i3)%b(1)+l7_516(i2)%a(2)*u7516_34(i3)%a(2)
      end do
      end do
      endif
      endif
  
* W 7\~gluon                                                            
  
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p51,q=p5178
      quqd=p51(0)*p5178(0)-p51(1)*p5178(1)-p51(2)*p5178(2)-p51(3
     & )*p5178(3)
      ccl=wcl/(f5178)
* TW0 -- qu=p51,qd=p5178,v=cw78%e,a=u51_78%a,b=u51_78%b,c=u51_78%c,d=u51
* _78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p51(2)*p5178(3)-p5178(2)*p51(3))+p51k0*(
     & cw78%e(2)*p5178(3)-p5178(2)*cw78%e(3))-p5178k0*(cw78%e(2)
     & *p51(3)-p51(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p51k0+p51(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p5178k0+p5178(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p51(0)-cw78%e(1)*p51(1)-cw78%e(2)*p51(2)-cw
     & 78%e(3)*p51(3)
      cvqd=cw78%e(0)*p5178(0)-cw78%e(1)*p5178(1)-cw78%e(2)*p5178
     & (2)-cw78%e(3)*p5178(3)
      cauxa=-cw78%ek0*quqd+p51k0*cvqd+p5178k0*cvqu
      cauxb=-cw78%ek0*p5178(2)+p5178k0*cw78%e(2)
      cauxc=+cw78%ek0*p51(2)-p51k0*cw78%e(2)
      u51_78%a(2)=ccl*(cauxa-ceps_0)
      u51_78%b(1)=ccl*(cauxb-ceps_2)
      u51_78%c(2)=ccl*(-cauxc+ceps_1)
      u51_78%d(1)=ccl*cw78%ek0
  
      do i1=1,2
* TLT0_W -- aa=l5_178(i1)%a,cc=l5_178(i1)%c,a1=l5_1(i1)%a,c1=l5_1(i1)%c,
* a2=u51_78%a,b2=u51_78%b,c2=u51_78%c,d2=u51_78%d,prq=s51,nsum=0
      l5_178(i1)%c(2)=l5_1(i1)%c(2)*s51*u51_78%d(1)+l5_1(i1)%a(2
     & )*u51_78%c(2)
      l5_178(i1)%a(2)=l5_1(i1)%c(2)*s51*u51_78%b(1)+l5_1(i1)%a(2
     & )*u51_78%a(2)
      end do
  
* quqd -- p=p578,q=p5178
      quqd=p578(0)*p5178(0)-p578(1)*p5178(1)-p578(2)*p5178(2)-p5
     & 78(3)*p5178(3)
      ccl=1.d0/(f5178)
      do i1=1,2
* TW0 -- qu=p578,qd=p5178,v=ce1(i1)%e,a=u578_1(i1)%a,b=u578_1(i1)%b,c=u5
* 78_1(i1)%c,d=u578_1(i1)%d,cl=ccl,nsum=0
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
      u578_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u578_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u578_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u578_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      do i1=1,2
* TLT0_W -- aa=l5_178(i1)%a,cc=l5_178(i1)%c,a1=l5_78%a,c1=l5_78%c,a2=u57
* 8_1(i1)%a,b2=u578_1(i1)%b,c2=u578_1(i1)%c,d2=u578_1(i1)%d,prq=s578,nsu
* m=1
      l5_178(i1)%c(2)=l5_178(i1)%c(2)+l5_78%c(2)*s578*u578_1(i1)
     & %d(1)+l5_78%a(2)*u578_1(i1)%c(2)
      l5_178(i1)%a(2)=l5_178(i1)%a(2)+l5_78%c(2)*s578*u578_1(i1)
     & %b(1)+l5_78%a(2)*u578_1(i1)%a(2)
      end do
  
* II                                                                    
      if (ilept(id3).ne.1) then
* quqd -- p=p51,q=p678
      quqd=p51(0)*p678(0)-p51(1)*p678(1)-p51(2)*p678(2)-p51(3)*p
     & 678(3)
      ccl=zcl(id5)/(f678)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p51,qd=p678,v=cz324(i3,i2)%e,a=u51_324(i3,i2)%a,b=u51_324(i3
* ,i2)%b,c=u51_324(i3,i2)%c,d=u51_324(i3,i2)%d,cl=ccl,nsum=0
      ceps_0=-cz324(i3,i2)%ek0*(p51(2)*p678(3)-p678(2)*p51(3))+p
     & 51k0*(cz324(i3,i2)%e(2)*p678(3)-p678(2)*cz324(i3,i2)%e(3)
     & )-p678k0*(cz324(i3,i2)%e(2)*p51(3)-p51(2)*cz324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i3,i2)%e(3)*p51k0+p51(3)*cz324(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz324(i3,i2)%e(3)*p678k0+p678(3)*cz324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i3,i2)%e(0)*p51(0)-cz324(i3,i2)%e(1)*p51(1)-cz3
     & 24(i3,i2)%e(2)*p51(2)-cz324(i3,i2)%e(3)*p51(3)
      cvqd=cz324(i3,i2)%e(0)*p678(0)-cz324(i3,i2)%e(1)*p678(1)-c
     & z324(i3,i2)%e(2)*p678(2)-cz324(i3,i2)%e(3)*p678(3)
      cauxa=-cz324(i3,i2)%ek0*quqd+p51k0*cvqd+p678k0*cvqu
      cauxb=-cz324(i3,i2)%ek0*p678(2)+p678k0*cz324(i3,i2)%e(2)
      cauxc=+cz324(i3,i2)%ek0*p51(2)-p51k0*cz324(i3,i2)%e(2)
      u51_324(i3,i2)%a(2)=ccl*(cauxa-ceps_0)
      u51_324(i3,i2)%b(1)=ccl*(cauxb-ceps_2)
      u51_324(i3,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u51_324(i3,i2)%d(1)=ccl*cz324(i3,i2)%ek0
      end do
      end do
  
      ccl=fcl(id5)/(f678)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p51,qd=p678,v=cf324(i3,i2)%e,a=u51_324(i3,i2)%a,b=u51_324(i3
* ,i2)%b,c=u51_324(i3,i2)%c,d=u51_324(i3,i2)%d,cl=ccl,nsum=1
      ceps_0=-cf324(i3,i2)%ek0*(p51(2)*p678(3)-p678(2)*p51(3))+p
     & 51k0*(cf324(i3,i2)%e(2)*p678(3)-p678(2)*cf324(i3,i2)%e(3)
     & )-p678k0*(cf324(i3,i2)%e(2)*p51(3)-p51(2)*cf324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i3,i2)%e(3)*p51k0+p51(3)*cf324(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf324(i3,i2)%e(3)*p678k0+p678(3)*cf324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i3,i2)%e(0)*p51(0)-cf324(i3,i2)%e(1)*p51(1)-cf3
     & 24(i3,i2)%e(2)*p51(2)-cf324(i3,i2)%e(3)*p51(3)
      cvqd=cf324(i3,i2)%e(0)*p678(0)-cf324(i3,i2)%e(1)*p678(1)-c
     & f324(i3,i2)%e(2)*p678(2)-cf324(i3,i2)%e(3)*p678(3)
      cauxa=-cf324(i3,i2)%ek0*quqd+p51k0*cvqd+p678k0*cvqu
      cauxb=-cf324(i3,i2)%ek0*p678(2)+p678k0*cf324(i3,i2)%e(2)
      cauxc=+cf324(i3,i2)%ek0*p51(2)-p51k0*cf324(i3,i2)%e(2)
      u51_324(i3,i2)%a(2)=u51_324(i3,i2)%a(2)+ccl*(cauxa-ceps_0)
      u51_324(i3,i2)%b(1)=u51_324(i3,i2)%b(1)+ccl*(cauxb-ceps_2)
      u51_324(i3,i2)%c(2)=u51_324(i3,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u51_324(i3,i2)%d(1)=u51_324(i3,i2)%d(1)+ccl*cf324(i3,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_1324(i1,i3,i2)%a,cc=l5_1324(i1,i3,i2)%c,a1=l5_1(i1)%a,
* c1=l5_1(i1)%c,a2=u51_324(i3,i2)%a,b2=u51_324(i3,i2)%b,c2=u51_324(i3,i2
* )%c,d2=u51_324(i3,i2)%d,prq=s51,nsum=0
      l5_1324(i1,i3,i2)%c(2)=l5_1(i1)%c(2)*s51*u51_324(i3,i2)%d(
     & 1)+l5_1(i1)%a(2)*u51_324(i3,i2)%c(2)
      l5_1324(i1,i3,i2)%a(2)=l5_1(i1)%c(2)*s51*u51_324(i3,i2)%b(
     & 1)+l5_1(i1)%a(2)*u51_324(i3,i2)%a(2)
      end do
      end do
      end do
      endif
  
* To use also as u5234_1                                                
* quqd -- p=p5234,q=p678
      quqd=p5234(0)*p678(0)-p5234(1)*p678(1)-p5234(2)*p678(2)-p5
     & 234(3)*p678(3)
      ccl=1.d0/(f678)
      do i1=1,2
* TW0 -- qu=p5234,qd=p678,v=ce1(i1)%e,a=u5324_1(i1)%a,b=u5324_1(i1)%b,c=
* u5324_1(i1)%c,d=u5324_1(i1)%d,cl=ccl,nsum=0
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
      u5324_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u5324_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u5324_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u5324_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_1324(i1,i3,i2)%a,cc=l5_1324(i1,i3,i2)%c,a1=l5_324(i3,i
* 2)%a,c1=l5_324(i3,i2)%c,a2=u5324_1(i1)%a,b2=u5324_1(i1)%b,c2=u5324_1(i
* 1)%c,d2=u5324_1(i1)%d,prq=s5234,nsum=1
      l5_1324(i1,i3,i2)%c(2)=l5_1324(i1,i3,i2)%c(2)+l5_324(i3,i2
     & )%c(2)*s5234*u5324_1(i1)%d(1)+l5_324(i3,i2)%a(2)*u5324_1(
     & i1)%c(2)
      l5_1324(i1,i3,i2)%a(2)=l5_1324(i1,i3,i2)%a(2)+l5_324(i3,i2
     & )%c(2)*s5234*u5324_1(i1)%b(1)+l5_324(i3,i2)%a(2)*u5324_1(
     & i1)%a(2)
      end do
      end do
      end do
  
* III                                                                   
* quqd -- p=p578,q=p61
      quqd=p578(0)*p61(0)-p578(1)*p61(1)-p578(2)*p61(2)-p578(3)*
     & p61(3)
      ccl=zcl(id6)/(f61)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p578,qd=p61,v=cz324(i3,i2)%e,a=u578_324(i3,i2)%a,b=u578_324(
* i3,i2)%b,c=u578_324(i3,i2)%c,d=u578_324(i3,i2)%d,cl=ccl,nsum=0
      ceps_0=-cz324(i3,i2)%ek0*(p578(2)*p61(3)-p61(2)*p578(3))+p
     & 578k0*(cz324(i3,i2)%e(2)*p61(3)-p61(2)*cz324(i3,i2)%e(3))
     & -p61k0*(cz324(i3,i2)%e(2)*p578(3)-p578(2)*cz324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i3,i2)%e(3)*p578k0+p578(3)*cz324(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz324(i3,i2)%e(3)*p61k0+p61(3)*cz324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i3,i2)%e(0)*p578(0)-cz324(i3,i2)%e(1)*p578(1)-c
     & z324(i3,i2)%e(2)*p578(2)-cz324(i3,i2)%e(3)*p578(3)
      cvqd=cz324(i3,i2)%e(0)*p61(0)-cz324(i3,i2)%e(1)*p61(1)-cz3
     & 24(i3,i2)%e(2)*p61(2)-cz324(i3,i2)%e(3)*p61(3)
      cauxa=-cz324(i3,i2)%ek0*quqd+p578k0*cvqd+p61k0*cvqu
      cauxb=-cz324(i3,i2)%ek0*p61(2)+p61k0*cz324(i3,i2)%e(2)
      cauxc=+cz324(i3,i2)%ek0*p578(2)-p578k0*cz324(i3,i2)%e(2)
      u578_324(i3,i2)%a(2)=ccl*(cauxa-ceps_0)
      u578_324(i3,i2)%b(1)=ccl*(cauxb-ceps_2)
      u578_324(i3,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u578_324(i3,i2)%d(1)=ccl*cz324(i3,i2)%ek0
      end do
      end do
  
      ccl=fcl(id6)/(f61)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p578,qd=p61,v=cf324(i3,i2)%e,a=u578_324(i3,i2)%a,b=u578_324(
* i3,i2)%b,c=u578_324(i3,i2)%c,d=u578_324(i3,i2)%d,cl=ccl,nsum=1
      ceps_0=-cf324(i3,i2)%ek0*(p578(2)*p61(3)-p61(2)*p578(3))+p
     & 578k0*(cf324(i3,i2)%e(2)*p61(3)-p61(2)*cf324(i3,i2)%e(3))
     & -p61k0*(cf324(i3,i2)%e(2)*p578(3)-p578(2)*cf324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i3,i2)%e(3)*p578k0+p578(3)*cf324(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf324(i3,i2)%e(3)*p61k0+p61(3)*cf324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i3,i2)%e(0)*p578(0)-cf324(i3,i2)%e(1)*p578(1)-c
     & f324(i3,i2)%e(2)*p578(2)-cf324(i3,i2)%e(3)*p578(3)
      cvqd=cf324(i3,i2)%e(0)*p61(0)-cf324(i3,i2)%e(1)*p61(1)-cf3
     & 24(i3,i2)%e(2)*p61(2)-cf324(i3,i2)%e(3)*p61(3)
      cauxa=-cf324(i3,i2)%ek0*quqd+p578k0*cvqd+p61k0*cvqu
      cauxb=-cf324(i3,i2)%ek0*p61(2)+p61k0*cf324(i3,i2)%e(2)
      cauxc=+cf324(i3,i2)%ek0*p578(2)-p578k0*cf324(i3,i2)%e(2)
      u578_324(i3,i2)%a(2)=u578_324(i3,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u578_324(i3,i2)%b(1)=u578_324(i3,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u578_324(i3,i2)%c(2)=u578_324(i3,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u578_324(i3,i2)%d(1)=u578_324(i3,i2)%d(1)+ccl*cf324(i3,i2)
     & %ek0
      end do
      end do
  
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_78324(i3,i2)%a,cc=l5_78324(i3,i2)%c,a1=l5_78%a,c1=l5_7
* 8%c,a2=u578_324(i3,i2)%a,b2=u578_324(i3,i2)%b,c2=u578_324(i3,i2)%c,d2=
* u578_324(i3,i2)%d,prq=s578,nsum=0
      l5_78324(i3,i2)%c(2)=l5_78%c(2)*s578*u578_324(i3,i2)%d(1)+
     & l5_78%a(2)*u578_324(i3,i2)%c(2)
      l5_78324(i3,i2)%a(2)=l5_78%c(2)*s578*u578_324(i3,i2)%b(1)+
     & l5_78%a(2)*u578_324(i3,i2)%a(2)
      end do
      end do
       endif
  
* To use also as u5234_78                                               
* quqd -- p=p5234,q=p61
      quqd=p5234(0)*p61(0)-p5234(1)*p61(1)-p5234(2)*p61(2)-p5234
     & (3)*p61(3)
      ccl=wcl/(f61)
* TW0 -- qu=p5234,qd=p61,v=cw78%e,a=u5324_78%a,b=u5324_78%b,c=u5324_78%c
* ,d=u5324_78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p5234(2)*p61(3)-p61(2)*p5234(3))+p5234k0
     & *(cw78%e(2)*p61(3)-p61(2)*cw78%e(3))-p61k0*(cw78%e(2)*p52
     & 34(3)-p5234(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p5234k0+p5234(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p61k0+p61(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p5234(0)-cw78%e(1)*p5234(1)-cw78%e(2)*p5234
     & (2)-cw78%e(3)*p5234(3)
      cvqd=cw78%e(0)*p61(0)-cw78%e(1)*p61(1)-cw78%e(2)*p61(2)-cw
     & 78%e(3)*p61(3)
      cauxa=-cw78%ek0*quqd+p5234k0*cvqd+p61k0*cvqu
      cauxb=-cw78%ek0*p61(2)+p61k0*cw78%e(2)
      cauxc=+cw78%ek0*p5234(2)-p5234k0*cw78%e(2)
      u5324_78%a(2)=ccl*(cauxa-ceps_0)
      u5324_78%b(1)=ccl*(cauxb-ceps_2)
      u5324_78%c(2)=ccl*(-cauxc+ceps_1)
      u5324_78%d(1)=ccl*cw78%ek0
  
      if (ilept(id3).ne.1) then
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_78324(i3,i2)%a,cc=l5_78324(i3,i2)%c,a1=l5_324(i3,i2)%a
* ,c1=l5_324(i3,i2)%c,a2=u5324_78%a,b2=u5324_78%b,c2=u5324_78%c,d2=u5324
* _78%d,prq=s5234,nsum=1
      l5_78324(i3,i2)%c(2)=l5_78324(i3,i2)%c(2)+l5_324(i3,i2)%c(
     & 2)*s5234*u5324_78%d(1)+l5_324(i3,i2)%a(2)*u5324_78%c(2)
      l5_78324(i3,i2)%a(2)=l5_78324(i3,i2)%a(2)+l5_324(i3,i2)%c(
     & 2)*s5234*u5324_78%b(1)+l5_324(i3,i2)%a(2)*u5324_78%a(2)
      end do
      end do
      endif
      endif
  
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p52,q=p5278
      quqd=p52(0)*p5278(0)-p52(1)*p5278(1)-p52(2)*p5278(2)-p52(3
     & )*p5278(3)
      ccl=wcl/(f5278)
* TW0 -- qu=p52,qd=p5278,v=cw78%e,a=u52_78%a,b=u52_78%b,c=u52_78%c,d=u52
* _78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p52(2)*p5278(3)-p5278(2)*p52(3))+p52k0*(
     & cw78%e(2)*p5278(3)-p5278(2)*cw78%e(3))-p5278k0*(cw78%e(2)
     & *p52(3)-p52(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p52k0+p52(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p5278k0+p5278(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p52(0)-cw78%e(1)*p52(1)-cw78%e(2)*p52(2)-cw
     & 78%e(3)*p52(3)
      cvqd=cw78%e(0)*p5278(0)-cw78%e(1)*p5278(1)-cw78%e(2)*p5278
     & (2)-cw78%e(3)*p5278(3)
      cauxa=-cw78%ek0*quqd+p52k0*cvqd+p5278k0*cvqu
      cauxb=-cw78%ek0*p5278(2)+p5278k0*cw78%e(2)
      cauxc=+cw78%ek0*p52(2)-p52k0*cw78%e(2)
      u52_78%a(2)=ccl*(cauxa-ceps_0)
      u52_78%b(1)=ccl*(cauxb-ceps_2)
      u52_78%c(2)=ccl*(-cauxc+ceps_1)
      u52_78%d(1)=ccl*cw78%ek0
  
      do i1=1,2
* TLT0_W -- aa=l5_278(i1)%a,cc=l5_278(i1)%c,a1=l5_2(i1)%a,c1=l5_2(i1)%c,
* a2=u52_78%a,b2=u52_78%b,c2=u52_78%c,d2=u52_78%d,prq=s52,nsum=0
      l5_278(i1)%c(2)=l5_2(i1)%c(2)*s52*u52_78%d(1)+l5_2(i1)%a(2
     & )*u52_78%c(2)
      l5_278(i1)%a(2)=l5_2(i1)%c(2)*s52*u52_78%b(1)+l5_2(i1)%a(2
     & )*u52_78%a(2)
      end do
  
* quqd -- p=p578,q=p5278
      quqd=p578(0)*p5278(0)-p578(1)*p5278(1)-p578(2)*p5278(2)-p5
     & 78(3)*p5278(3)
      ccl=1.d0/(f5278)
      do i1=1,2
* TW0 -- qu=p578,qd=p5278,v=ce2(i1)%e,a=u578_2(i1)%a,b=u578_2(i1)%b,c=u5
* 78_2(i1)%c,d=u578_2(i1)%d,cl=ccl,nsum=0
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
      u578_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u578_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u578_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u578_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      do i1=1,2
* TLT0_W -- aa=l5_278(i1)%a,cc=l5_278(i1)%c,a1=l5_78%a,c1=l5_78%c,a2=u57
* 8_2(i1)%a,b2=u578_2(i1)%b,c2=u578_2(i1)%c,d2=u578_2(i1)%d,prq=s578,nsu
* m=1
      l5_278(i1)%c(2)=l5_278(i1)%c(2)+l5_78%c(2)*s578*u578_2(i1)
     & %d(1)+l5_78%a(2)*u578_2(i1)%c(2)
      l5_278(i1)%a(2)=l5_278(i1)%a(2)+l5_78%c(2)*s578*u578_2(i1)
     & %b(1)+l5_78%a(2)*u578_2(i1)%a(2)
      end do
  
* II                                                                    
      if (ilept(id3).ne.1) then
* quqd -- p=p52,q=p678
      quqd=p52(0)*p678(0)-p52(1)*p678(1)-p52(2)*p678(2)-p52(3)*p
     & 678(3)
      ccl=zcl(id5)/(f678)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p52,qd=p678,v=cz314(i3,i2)%e,a=u52_314(i3,i2)%a,b=u52_314(i3
* ,i2)%b,c=u52_314(i3,i2)%c,d=u52_314(i3,i2)%d,cl=ccl,nsum=0
      ceps_0=-cz314(i3,i2)%ek0*(p52(2)*p678(3)-p678(2)*p52(3))+p
     & 52k0*(cz314(i3,i2)%e(2)*p678(3)-p678(2)*cz314(i3,i2)%e(3)
     & )-p678k0*(cz314(i3,i2)%e(2)*p52(3)-p52(2)*cz314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i3,i2)%e(3)*p52k0+p52(3)*cz314(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz314(i3,i2)%e(3)*p678k0+p678(3)*cz314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i3,i2)%e(0)*p52(0)-cz314(i3,i2)%e(1)*p52(1)-cz3
     & 14(i3,i2)%e(2)*p52(2)-cz314(i3,i2)%e(3)*p52(3)
      cvqd=cz314(i3,i2)%e(0)*p678(0)-cz314(i3,i2)%e(1)*p678(1)-c
     & z314(i3,i2)%e(2)*p678(2)-cz314(i3,i2)%e(3)*p678(3)
      cauxa=-cz314(i3,i2)%ek0*quqd+p52k0*cvqd+p678k0*cvqu
      cauxb=-cz314(i3,i2)%ek0*p678(2)+p678k0*cz314(i3,i2)%e(2)
      cauxc=+cz314(i3,i2)%ek0*p52(2)-p52k0*cz314(i3,i2)%e(2)
      u52_314(i3,i2)%a(2)=ccl*(cauxa-ceps_0)
      u52_314(i3,i2)%b(1)=ccl*(cauxb-ceps_2)
      u52_314(i3,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u52_314(i3,i2)%d(1)=ccl*cz314(i3,i2)%ek0
      end do
      end do
  
      ccl=fcl(id5)/(f678)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p52,qd=p678,v=cf314(i3,i2)%e,a=u52_314(i3,i2)%a,b=u52_314(i3
* ,i2)%b,c=u52_314(i3,i2)%c,d=u52_314(i3,i2)%d,cl=ccl,nsum=1
      ceps_0=-cf314(i3,i2)%ek0*(p52(2)*p678(3)-p678(2)*p52(3))+p
     & 52k0*(cf314(i3,i2)%e(2)*p678(3)-p678(2)*cf314(i3,i2)%e(3)
     & )-p678k0*(cf314(i3,i2)%e(2)*p52(3)-p52(2)*cf314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i3,i2)%e(3)*p52k0+p52(3)*cf314(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf314(i3,i2)%e(3)*p678k0+p678(3)*cf314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i3,i2)%e(0)*p52(0)-cf314(i3,i2)%e(1)*p52(1)-cf3
     & 14(i3,i2)%e(2)*p52(2)-cf314(i3,i2)%e(3)*p52(3)
      cvqd=cf314(i3,i2)%e(0)*p678(0)-cf314(i3,i2)%e(1)*p678(1)-c
     & f314(i3,i2)%e(2)*p678(2)-cf314(i3,i2)%e(3)*p678(3)
      cauxa=-cf314(i3,i2)%ek0*quqd+p52k0*cvqd+p678k0*cvqu
      cauxb=-cf314(i3,i2)%ek0*p678(2)+p678k0*cf314(i3,i2)%e(2)
      cauxc=+cf314(i3,i2)%ek0*p52(2)-p52k0*cf314(i3,i2)%e(2)
      u52_314(i3,i2)%a(2)=u52_314(i3,i2)%a(2)+ccl*(cauxa-ceps_0)
      u52_314(i3,i2)%b(1)=u52_314(i3,i2)%b(1)+ccl*(cauxb-ceps_2)
      u52_314(i3,i2)%c(2)=u52_314(i3,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u52_314(i3,i2)%d(1)=u52_314(i3,i2)%d(1)+ccl*cf314(i3,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_2314(i1,i3,i2)%a,cc=l5_2314(i1,i3,i2)%c,a1=l5_2(i1)%a,
* c1=l5_2(i1)%c,a2=u52_314(i3,i2)%a,b2=u52_314(i3,i2)%b,c2=u52_314(i3,i2
* )%c,d2=u52_314(i3,i2)%d,prq=s52,nsum=0
      l5_2314(i1,i3,i2)%c(2)=l5_2(i1)%c(2)*s52*u52_314(i3,i2)%d(
     & 1)+l5_2(i1)%a(2)*u52_314(i3,i2)%c(2)
      l5_2314(i1,i3,i2)%a(2)=l5_2(i1)%c(2)*s52*u52_314(i3,i2)%b(
     & 1)+l5_2(i1)%a(2)*u52_314(i3,i2)%a(2)
      end do
      end do
      end do
      endif
  
* To use also as u5134_2                                                
* quqd -- p=p5134,q=p678
      quqd=p5134(0)*p678(0)-p5134(1)*p678(1)-p5134(2)*p678(2)-p5
     & 134(3)*p678(3)
      ccl=1.d0/(f678)
      do i1=1,2
* TW0 -- qu=p5134,qd=p678,v=ce2(i1)%e,a=u5314_2(i1)%a,b=u5314_2(i1)%b,c=
* u5314_2(i1)%c,d=u5314_2(i1)%d,cl=ccl,nsum=0
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
      u5314_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u5314_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u5314_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u5314_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_2314(i1,i3,i2)%a,cc=l5_2314(i1,i3,i2)%c,a1=l5_314(i3,i
* 2)%a,c1=l5_314(i3,i2)%c,a2=u5314_2(i1)%a,b2=u5314_2(i1)%b,c2=u5314_2(i
* 1)%c,d2=u5314_2(i1)%d,prq=s5134,nsum=1
      l5_2314(i1,i3,i2)%c(2)=l5_2314(i1,i3,i2)%c(2)+l5_314(i3,i2
     & )%c(2)*s5134*u5314_2(i1)%d(1)+l5_314(i3,i2)%a(2)*u5314_2(
     & i1)%c(2)
      l5_2314(i1,i3,i2)%a(2)=l5_2314(i1,i3,i2)%a(2)+l5_314(i3,i2
     & )%c(2)*s5134*u5314_2(i1)%b(1)+l5_314(i3,i2)%a(2)*u5314_2(
     & i1)%a(2)
      end do
      end do
      end do
  
* III                                                                   
* quqd -- p=p578,q=p62
      quqd=p578(0)*p62(0)-p578(1)*p62(1)-p578(2)*p62(2)-p578(3)*
     & p62(3)
      ccl=zcl(id6)/(f62)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p578,qd=p62,v=cz314(i3,i2)%e,a=u578_314(i3,i2)%a,b=u578_314(
* i3,i2)%b,c=u578_314(i3,i2)%c,d=u578_314(i3,i2)%d,cl=ccl,nsum=0
      ceps_0=-cz314(i3,i2)%ek0*(p578(2)*p62(3)-p62(2)*p578(3))+p
     & 578k0*(cz314(i3,i2)%e(2)*p62(3)-p62(2)*cz314(i3,i2)%e(3))
     & -p62k0*(cz314(i3,i2)%e(2)*p578(3)-p578(2)*cz314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i3,i2)%e(3)*p578k0+p578(3)*cz314(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz314(i3,i2)%e(3)*p62k0+p62(3)*cz314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i3,i2)%e(0)*p578(0)-cz314(i3,i2)%e(1)*p578(1)-c
     & z314(i3,i2)%e(2)*p578(2)-cz314(i3,i2)%e(3)*p578(3)
      cvqd=cz314(i3,i2)%e(0)*p62(0)-cz314(i3,i2)%e(1)*p62(1)-cz3
     & 14(i3,i2)%e(2)*p62(2)-cz314(i3,i2)%e(3)*p62(3)
      cauxa=-cz314(i3,i2)%ek0*quqd+p578k0*cvqd+p62k0*cvqu
      cauxb=-cz314(i3,i2)%ek0*p62(2)+p62k0*cz314(i3,i2)%e(2)
      cauxc=+cz314(i3,i2)%ek0*p578(2)-p578k0*cz314(i3,i2)%e(2)
      u578_314(i3,i2)%a(2)=ccl*(cauxa-ceps_0)
      u578_314(i3,i2)%b(1)=ccl*(cauxb-ceps_2)
      u578_314(i3,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u578_314(i3,i2)%d(1)=ccl*cz314(i3,i2)%ek0
      end do
      end do
  
      ccl=fcl(id6)/(f62)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p578,qd=p62,v=cf314(i3,i2)%e,a=u578_314(i3,i2)%a,b=u578_314(
* i3,i2)%b,c=u578_314(i3,i2)%c,d=u578_314(i3,i2)%d,cl=ccl,nsum=1
      ceps_0=-cf314(i3,i2)%ek0*(p578(2)*p62(3)-p62(2)*p578(3))+p
     & 578k0*(cf314(i3,i2)%e(2)*p62(3)-p62(2)*cf314(i3,i2)%e(3))
     & -p62k0*(cf314(i3,i2)%e(2)*p578(3)-p578(2)*cf314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i3,i2)%e(3)*p578k0+p578(3)*cf314(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf314(i3,i2)%e(3)*p62k0+p62(3)*cf314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i3,i2)%e(0)*p578(0)-cf314(i3,i2)%e(1)*p578(1)-c
     & f314(i3,i2)%e(2)*p578(2)-cf314(i3,i2)%e(3)*p578(3)
      cvqd=cf314(i3,i2)%e(0)*p62(0)-cf314(i3,i2)%e(1)*p62(1)-cf3
     & 14(i3,i2)%e(2)*p62(2)-cf314(i3,i2)%e(3)*p62(3)
      cauxa=-cf314(i3,i2)%ek0*quqd+p578k0*cvqd+p62k0*cvqu
      cauxb=-cf314(i3,i2)%ek0*p62(2)+p62k0*cf314(i3,i2)%e(2)
      cauxc=+cf314(i3,i2)%ek0*p578(2)-p578k0*cf314(i3,i2)%e(2)
      u578_314(i3,i2)%a(2)=u578_314(i3,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u578_314(i3,i2)%b(1)=u578_314(i3,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u578_314(i3,i2)%c(2)=u578_314(i3,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u578_314(i3,i2)%d(1)=u578_314(i3,i2)%d(1)+ccl*cf314(i3,i2)
     & %ek0
      end do
      end do
  
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_78314(i3,i2)%a,cc=l5_78314(i3,i2)%c,a1=l5_78%a,c1=l5_7
* 8%c,a2=u578_314(i3,i2)%a,b2=u578_314(i3,i2)%b,c2=u578_314(i3,i2)%c,d2=
* u578_314(i3,i2)%d,prq=s578,nsum=0
      l5_78314(i3,i2)%c(2)=l5_78%c(2)*s578*u578_314(i3,i2)%d(1)+
     & l5_78%a(2)*u578_314(i3,i2)%c(2)
      l5_78314(i3,i2)%a(2)=l5_78%c(2)*s578*u578_314(i3,i2)%b(1)+
     & l5_78%a(2)*u578_314(i3,i2)%a(2)
      end do
      end do
       endif
  
* To use also as u5134_78                                               
* quqd -- p=p5134,q=p62
      quqd=p5134(0)*p62(0)-p5134(1)*p62(1)-p5134(2)*p62(2)-p5134
     & (3)*p62(3)
      ccl=wcl/(f62)
* TW0 -- qu=p5134,qd=p62,v=cw78%e,a=u5314_78%a,b=u5314_78%b,c=u5314_78%c
* ,d=u5314_78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p5134(2)*p62(3)-p62(2)*p5134(3))+p5134k0
     & *(cw78%e(2)*p62(3)-p62(2)*cw78%e(3))-p62k0*(cw78%e(2)*p51
     & 34(3)-p5134(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p5134k0+p5134(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p62k0+p62(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p5134(0)-cw78%e(1)*p5134(1)-cw78%e(2)*p5134
     & (2)-cw78%e(3)*p5134(3)
      cvqd=cw78%e(0)*p62(0)-cw78%e(1)*p62(1)-cw78%e(2)*p62(2)-cw
     & 78%e(3)*p62(3)
      cauxa=-cw78%ek0*quqd+p5134k0*cvqd+p62k0*cvqu
      cauxb=-cw78%ek0*p62(2)+p62k0*cw78%e(2)
      cauxc=+cw78%ek0*p5134(2)-p5134k0*cw78%e(2)
      u5314_78%a(2)=ccl*(cauxa-ceps_0)
      u5314_78%b(1)=ccl*(cauxb-ceps_2)
      u5314_78%c(2)=ccl*(-cauxc+ceps_1)
      u5314_78%d(1)=ccl*cw78%ek0
  
      if (ilept(id3).ne.1) then
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l5_78314(i3,i2)%a,cc=l5_78314(i3,i2)%c,a1=l5_314(i3,i2)%a
* ,c1=l5_314(i3,i2)%c,a2=u5314_78%a,b2=u5314_78%b,c2=u5314_78%c,d2=u5314
* _78%d,prq=s5134,nsum=1
      l5_78314(i3,i2)%c(2)=l5_78314(i3,i2)%c(2)+l5_314(i3,i2)%c(
     & 2)*s5134*u5314_78%d(1)+l5_314(i3,i2)%a(2)*u5314_78%c(2)
      l5_78314(i3,i2)%a(2)=l5_78314(i3,i2)%a(2)+l5_314(i3,i2)%c(
     & 2)*s5134*u5314_78%b(1)+l5_314(i3,i2)%a(2)*u5314_78%a(2)
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p71,q=p7156
      quqd=p71(0)*p7156(0)-p71(1)*p7156(1)-p71(2)*p7156(2)-p71(3
     & )*p7156(3)
      ccl=wcl/(f7156)
* TW0 -- qu=p71,qd=p7156,v=cw56%e,a=u71_56%a,b=u71_56%b,c=u71_56%c,d=u71
* _56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p71(2)*p7156(3)-p7156(2)*p71(3))+p71k0*(
     & cw56%e(2)*p7156(3)-p7156(2)*cw56%e(3))-p7156k0*(cw56%e(2)
     & *p71(3)-p71(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p71k0+p71(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p7156k0+p7156(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p71(0)-cw56%e(1)*p71(1)-cw56%e(2)*p71(2)-cw
     & 56%e(3)*p71(3)
      cvqd=cw56%e(0)*p7156(0)-cw56%e(1)*p7156(1)-cw56%e(2)*p7156
     & (2)-cw56%e(3)*p7156(3)
      cauxa=-cw56%ek0*quqd+p71k0*cvqd+p7156k0*cvqu
      cauxb=-cw56%ek0*p7156(2)+p7156k0*cw56%e(2)
      cauxc=+cw56%ek0*p71(2)-p71k0*cw56%e(2)
      u71_56%a(2)=ccl*(cauxa-ceps_0)
      u71_56%b(1)=ccl*(cauxb-ceps_2)
      u71_56%c(2)=ccl*(-cauxc+ceps_1)
      u71_56%d(1)=ccl*cw56%ek0
  
      do i1=1,2
* TLT0_W -- aa=l7_156(i1)%a,cc=l7_156(i1)%c,a1=l7_1(i1)%a,c1=l7_1(i1)%c,
* a2=u71_56%a,b2=u71_56%b,c2=u71_56%c,d2=u71_56%d,prq=s71,nsum=0
      l7_156(i1)%c(2)=l7_1(i1)%c(2)*s71*u71_56%d(1)+l7_1(i1)%a(2
     & )*u71_56%c(2)
      l7_156(i1)%a(2)=l7_1(i1)%c(2)*s71*u71_56%b(1)+l7_1(i1)%a(2
     & )*u71_56%a(2)
      end do
  
* quqd -- p=p756,q=p7156
      quqd=p756(0)*p7156(0)-p756(1)*p7156(1)-p756(2)*p7156(2)-p7
     & 56(3)*p7156(3)
      ccl=1.d0/(f7156)
      do i1=1,2
* TW0 -- qu=p756,qd=p7156,v=ce1(i1)%e,a=u756_1(i1)%a,b=u756_1(i1)%b,c=u7
* 56_1(i1)%c,d=u756_1(i1)%d,cl=ccl,nsum=0
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
      u756_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u756_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u756_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u756_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      do i1=1,2
* TLT0_W -- aa=l7_156(i1)%a,cc=l7_156(i1)%c,a1=l7_56%a,c1=l7_56%c,a2=u75
* 6_1(i1)%a,b2=u756_1(i1)%b,c2=u756_1(i1)%c,d2=u756_1(i1)%d,prq=s756,nsu
* m=1
      l7_156(i1)%c(2)=l7_156(i1)%c(2)+l7_56%c(2)*s756*u756_1(i1)
     & %d(1)+l7_56%a(2)*u756_1(i1)%c(2)
      l7_156(i1)%a(2)=l7_156(i1)%a(2)+l7_56%c(2)*s756*u756_1(i1)
     & %b(1)+l7_56%a(2)*u756_1(i1)%a(2)
      end do
  
* II                                                                    
      if (ilept(id3).ne.1) then
* quqd -- p=p71,q=p856
      quqd=p71(0)*p856(0)-p71(1)*p856(1)-p71(2)*p856(2)-p71(3)*p
     & 856(3)
      ccl=zcl(id7)/(f856)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p71,qd=p856,v=cz324(i3,i2)%e,a=u71_324(i3,i2)%a,b=u71_324(i3
* ,i2)%b,c=u71_324(i3,i2)%c,d=u71_324(i3,i2)%d,cl=ccl,nsum=0
      ceps_0=-cz324(i3,i2)%ek0*(p71(2)*p856(3)-p856(2)*p71(3))+p
     & 71k0*(cz324(i3,i2)%e(2)*p856(3)-p856(2)*cz324(i3,i2)%e(3)
     & )-p856k0*(cz324(i3,i2)%e(2)*p71(3)-p71(2)*cz324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i3,i2)%e(3)*p71k0+p71(3)*cz324(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz324(i3,i2)%e(3)*p856k0+p856(3)*cz324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i3,i2)%e(0)*p71(0)-cz324(i3,i2)%e(1)*p71(1)-cz3
     & 24(i3,i2)%e(2)*p71(2)-cz324(i3,i2)%e(3)*p71(3)
      cvqd=cz324(i3,i2)%e(0)*p856(0)-cz324(i3,i2)%e(1)*p856(1)-c
     & z324(i3,i2)%e(2)*p856(2)-cz324(i3,i2)%e(3)*p856(3)
      cauxa=-cz324(i3,i2)%ek0*quqd+p71k0*cvqd+p856k0*cvqu
      cauxb=-cz324(i3,i2)%ek0*p856(2)+p856k0*cz324(i3,i2)%e(2)
      cauxc=+cz324(i3,i2)%ek0*p71(2)-p71k0*cz324(i3,i2)%e(2)
      u71_324(i3,i2)%a(2)=ccl*(cauxa-ceps_0)
      u71_324(i3,i2)%b(1)=ccl*(cauxb-ceps_2)
      u71_324(i3,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u71_324(i3,i2)%d(1)=ccl*cz324(i3,i2)%ek0
      end do
      end do
  
      ccl=fcl(id7)/(f856)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p71,qd=p856,v=cf324(i3,i2)%e,a=u71_324(i3,i2)%a,b=u71_324(i3
* ,i2)%b,c=u71_324(i3,i2)%c,d=u71_324(i3,i2)%d,cl=ccl,nsum=1
      ceps_0=-cf324(i3,i2)%ek0*(p71(2)*p856(3)-p856(2)*p71(3))+p
     & 71k0*(cf324(i3,i2)%e(2)*p856(3)-p856(2)*cf324(i3,i2)%e(3)
     & )-p856k0*(cf324(i3,i2)%e(2)*p71(3)-p71(2)*cf324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i3,i2)%e(3)*p71k0+p71(3)*cf324(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf324(i3,i2)%e(3)*p856k0+p856(3)*cf324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i3,i2)%e(0)*p71(0)-cf324(i3,i2)%e(1)*p71(1)-cf3
     & 24(i3,i2)%e(2)*p71(2)-cf324(i3,i2)%e(3)*p71(3)
      cvqd=cf324(i3,i2)%e(0)*p856(0)-cf324(i3,i2)%e(1)*p856(1)-c
     & f324(i3,i2)%e(2)*p856(2)-cf324(i3,i2)%e(3)*p856(3)
      cauxa=-cf324(i3,i2)%ek0*quqd+p71k0*cvqd+p856k0*cvqu
      cauxb=-cf324(i3,i2)%ek0*p856(2)+p856k0*cf324(i3,i2)%e(2)
      cauxc=+cf324(i3,i2)%ek0*p71(2)-p71k0*cf324(i3,i2)%e(2)
      u71_324(i3,i2)%a(2)=u71_324(i3,i2)%a(2)+ccl*(cauxa-ceps_0)
      u71_324(i3,i2)%b(1)=u71_324(i3,i2)%b(1)+ccl*(cauxb-ceps_2)
      u71_324(i3,i2)%c(2)=u71_324(i3,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u71_324(i3,i2)%d(1)=u71_324(i3,i2)%d(1)+ccl*cf324(i3,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_1324(i1,i3,i2)%a,cc=l7_1324(i1,i3,i2)%c,a1=l7_1(i1)%a,
* c1=l7_1(i1)%c,a2=u71_324(i3,i2)%a,b2=u71_324(i3,i2)%b,c2=u71_324(i3,i2
* )%c,d2=u71_324(i3,i2)%d,prq=s71,nsum=0
      l7_1324(i1,i3,i2)%c(2)=l7_1(i1)%c(2)*s71*u71_324(i3,i2)%d(
     & 1)+l7_1(i1)%a(2)*u71_324(i3,i2)%c(2)
      l7_1324(i1,i3,i2)%a(2)=l7_1(i1)%c(2)*s71*u71_324(i3,i2)%b(
     & 1)+l7_1(i1)%a(2)*u71_324(i3,i2)%a(2)
      end do
      end do
      end do
      endif
  
* To use also as u7234_1                                                
* quqd -- p=p7234,q=p856
      quqd=p7234(0)*p856(0)-p7234(1)*p856(1)-p7234(2)*p856(2)-p7
     & 234(3)*p856(3)
      ccl=1.d0/(f856)
      do i1=1,2
* TW0 -- qu=p7234,qd=p856,v=ce1(i1)%e,a=u7324_1(i1)%a,b=u7324_1(i1)%b,c=
* u7324_1(i1)%c,d=u7324_1(i1)%d,cl=ccl,nsum=0
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
      u7324_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u7324_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u7324_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u7324_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_1324(i1,i3,i2)%a,cc=l7_1324(i1,i3,i2)%c,a1=l7_324(i3,i
* 2)%a,c1=l7_324(i3,i2)%c,a2=u7324_1(i1)%a,b2=u7324_1(i1)%b,c2=u7324_1(i
* 1)%c,d2=u7324_1(i1)%d,prq=s7234,nsum=1
      l7_1324(i1,i3,i2)%c(2)=l7_1324(i1,i3,i2)%c(2)+l7_324(i3,i2
     & )%c(2)*s7234*u7324_1(i1)%d(1)+l7_324(i3,i2)%a(2)*u7324_1(
     & i1)%c(2)
      l7_1324(i1,i3,i2)%a(2)=l7_1324(i1,i3,i2)%a(2)+l7_324(i3,i2
     & )%c(2)*s7234*u7324_1(i1)%b(1)+l7_324(i3,i2)%a(2)*u7324_1(
     & i1)%a(2)
      end do
      end do
      end do
  
* III                                                                   
* quqd -- p=p756,q=p81
      quqd=p756(0)*p81(0)-p756(1)*p81(1)-p756(2)*p81(2)-p756(3)*
     & p81(3)
      ccl=zcl(id8)/(f81)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p756,qd=p81,v=cz324(i3,i2)%e,a=u756_324(i3,i2)%a,b=u756_324(
* i3,i2)%b,c=u756_324(i3,i2)%c,d=u756_324(i3,i2)%d,cl=ccl,nsum=0
      ceps_0=-cz324(i3,i2)%ek0*(p756(2)*p81(3)-p81(2)*p756(3))+p
     & 756k0*(cz324(i3,i2)%e(2)*p81(3)-p81(2)*cz324(i3,i2)%e(3))
     & -p81k0*(cz324(i3,i2)%e(2)*p756(3)-p756(2)*cz324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz324(i3,i2)%e(3)*p756k0+p756(3)*cz324(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz324(i3,i2)%e(3)*p81k0+p81(3)*cz324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz324(i3,i2)%e(0)*p756(0)-cz324(i3,i2)%e(1)*p756(1)-c
     & z324(i3,i2)%e(2)*p756(2)-cz324(i3,i2)%e(3)*p756(3)
      cvqd=cz324(i3,i2)%e(0)*p81(0)-cz324(i3,i2)%e(1)*p81(1)-cz3
     & 24(i3,i2)%e(2)*p81(2)-cz324(i3,i2)%e(3)*p81(3)
      cauxa=-cz324(i3,i2)%ek0*quqd+p756k0*cvqd+p81k0*cvqu
      cauxb=-cz324(i3,i2)%ek0*p81(2)+p81k0*cz324(i3,i2)%e(2)
      cauxc=+cz324(i3,i2)%ek0*p756(2)-p756k0*cz324(i3,i2)%e(2)
      u756_324(i3,i2)%a(2)=ccl*(cauxa-ceps_0)
      u756_324(i3,i2)%b(1)=ccl*(cauxb-ceps_2)
      u756_324(i3,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u756_324(i3,i2)%d(1)=ccl*cz324(i3,i2)%ek0
      end do
      end do
  
      ccl=fcl(id8)/(f81)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p756,qd=p81,v=cf324(i3,i2)%e,a=u756_324(i3,i2)%a,b=u756_324(
* i3,i2)%b,c=u756_324(i3,i2)%c,d=u756_324(i3,i2)%d,cl=ccl,nsum=1
      ceps_0=-cf324(i3,i2)%ek0*(p756(2)*p81(3)-p81(2)*p756(3))+p
     & 756k0*(cf324(i3,i2)%e(2)*p81(3)-p81(2)*cf324(i3,i2)%e(3))
     & -p81k0*(cf324(i3,i2)%e(2)*p756(3)-p756(2)*cf324(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf324(i3,i2)%e(3)*p756k0+p756(3)*cf324(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf324(i3,i2)%e(3)*p81k0+p81(3)*cf324(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf324(i3,i2)%e(0)*p756(0)-cf324(i3,i2)%e(1)*p756(1)-c
     & f324(i3,i2)%e(2)*p756(2)-cf324(i3,i2)%e(3)*p756(3)
      cvqd=cf324(i3,i2)%e(0)*p81(0)-cf324(i3,i2)%e(1)*p81(1)-cf3
     & 24(i3,i2)%e(2)*p81(2)-cf324(i3,i2)%e(3)*p81(3)
      cauxa=-cf324(i3,i2)%ek0*quqd+p756k0*cvqd+p81k0*cvqu
      cauxb=-cf324(i3,i2)%ek0*p81(2)+p81k0*cf324(i3,i2)%e(2)
      cauxc=+cf324(i3,i2)%ek0*p756(2)-p756k0*cf324(i3,i2)%e(2)
      u756_324(i3,i2)%a(2)=u756_324(i3,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u756_324(i3,i2)%b(1)=u756_324(i3,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u756_324(i3,i2)%c(2)=u756_324(i3,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u756_324(i3,i2)%d(1)=u756_324(i3,i2)%d(1)+ccl*cf324(i3,i2)
     & %ek0
      end do
      end do
  
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_56324(i3,i2)%a,cc=l7_56324(i3,i2)%c,a1=l7_56%a,c1=l7_5
* 6%c,a2=u756_324(i3,i2)%a,b2=u756_324(i3,i2)%b,c2=u756_324(i3,i2)%c,d2=
* u756_324(i3,i2)%d,prq=s756,nsum=0
      l7_56324(i3,i2)%c(2)=l7_56%c(2)*s756*u756_324(i3,i2)%d(1)+
     & l7_56%a(2)*u756_324(i3,i2)%c(2)
      l7_56324(i3,i2)%a(2)=l7_56%c(2)*s756*u756_324(i3,i2)%b(1)+
     & l7_56%a(2)*u756_324(i3,i2)%a(2)
      end do
      end do
       endif
  
* To use also as u7234_56                                               
* quqd -- p=p7234,q=p81
      quqd=p7234(0)*p81(0)-p7234(1)*p81(1)-p7234(2)*p81(2)-p7234
     & (3)*p81(3)
      ccl=wcl/(f81)
* TW0 -- qu=p7234,qd=p81,v=cw56%e,a=u7324_56%a,b=u7324_56%b,c=u7324_56%c
* ,d=u7324_56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p7234(2)*p81(3)-p81(2)*p7234(3))+p7234k0
     & *(cw56%e(2)*p81(3)-p81(2)*cw56%e(3))-p81k0*(cw56%e(2)*p72
     & 34(3)-p7234(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p7234k0+p7234(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p81k0+p81(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p7234(0)-cw56%e(1)*p7234(1)-cw56%e(2)*p7234
     & (2)-cw56%e(3)*p7234(3)
      cvqd=cw56%e(0)*p81(0)-cw56%e(1)*p81(1)-cw56%e(2)*p81(2)-cw
     & 56%e(3)*p81(3)
      cauxa=-cw56%ek0*quqd+p7234k0*cvqd+p81k0*cvqu
      cauxb=-cw56%ek0*p81(2)+p81k0*cw56%e(2)
      cauxc=+cw56%ek0*p7234(2)-p7234k0*cw56%e(2)
      u7324_56%a(2)=ccl*(cauxa-ceps_0)
      u7324_56%b(1)=ccl*(cauxb-ceps_2)
      u7324_56%c(2)=ccl*(-cauxc+ceps_1)
      u7324_56%d(1)=ccl*cw56%ek0
  
      if (ilept(id3).ne.1) then
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_56324(i3,i2)%a,cc=l7_56324(i3,i2)%c,a1=l7_324(i3,i2)%a
* ,c1=l7_324(i3,i2)%c,a2=u7324_56%a,b2=u7324_56%b,c2=u7324_56%c,d2=u7324
* _56%d,prq=s7234,nsum=1
      l7_56324(i3,i2)%c(2)=l7_56324(i3,i2)%c(2)+l7_324(i3,i2)%c(
     & 2)*s7234*u7324_56%d(1)+l7_324(i3,i2)%a(2)*u7324_56%c(2)
      l7_56324(i3,i2)%a(2)=l7_56324(i3,i2)%a(2)+l7_324(i3,i2)%c(
     & 2)*s7234*u7324_56%b(1)+l7_324(i3,i2)%a(2)*u7324_56%a(2)
      end do
      end do
      endif
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p72,q=p7256
      quqd=p72(0)*p7256(0)-p72(1)*p7256(1)-p72(2)*p7256(2)-p72(3
     & )*p7256(3)
      ccl=wcl/(f7256)
* TW0 -- qu=p72,qd=p7256,v=cw56%e,a=u72_56%a,b=u72_56%b,c=u72_56%c,d=u72
* _56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p72(2)*p7256(3)-p7256(2)*p72(3))+p72k0*(
     & cw56%e(2)*p7256(3)-p7256(2)*cw56%e(3))-p7256k0*(cw56%e(2)
     & *p72(3)-p72(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p72k0+p72(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p7256k0+p7256(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p72(0)-cw56%e(1)*p72(1)-cw56%e(2)*p72(2)-cw
     & 56%e(3)*p72(3)
      cvqd=cw56%e(0)*p7256(0)-cw56%e(1)*p7256(1)-cw56%e(2)*p7256
     & (2)-cw56%e(3)*p7256(3)
      cauxa=-cw56%ek0*quqd+p72k0*cvqd+p7256k0*cvqu
      cauxb=-cw56%ek0*p7256(2)+p7256k0*cw56%e(2)
      cauxc=+cw56%ek0*p72(2)-p72k0*cw56%e(2)
      u72_56%a(2)=ccl*(cauxa-ceps_0)
      u72_56%b(1)=ccl*(cauxb-ceps_2)
      u72_56%c(2)=ccl*(-cauxc+ceps_1)
      u72_56%d(1)=ccl*cw56%ek0
  
      do i1=1,2
* TLT0_W -- aa=l7_256(i1)%a,cc=l7_256(i1)%c,a1=l7_2(i1)%a,c1=l7_2(i1)%c,
* a2=u72_56%a,b2=u72_56%b,c2=u72_56%c,d2=u72_56%d,prq=s72,nsum=0
      l7_256(i1)%c(2)=l7_2(i1)%c(2)*s72*u72_56%d(1)+l7_2(i1)%a(2
     & )*u72_56%c(2)
      l7_256(i1)%a(2)=l7_2(i1)%c(2)*s72*u72_56%b(1)+l7_2(i1)%a(2
     & )*u72_56%a(2)
      end do
  
* quqd -- p=p756,q=p7256
      quqd=p756(0)*p7256(0)-p756(1)*p7256(1)-p756(2)*p7256(2)-p7
     & 56(3)*p7256(3)
      ccl=1.d0/(f7256)
      do i1=1,2
* TW0 -- qu=p756,qd=p7256,v=ce2(i1)%e,a=u756_2(i1)%a,b=u756_2(i1)%b,c=u7
* 56_2(i1)%c,d=u756_2(i1)%d,cl=ccl,nsum=0
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
      u756_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u756_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u756_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u756_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      do i1=1,2
* TLT0_W -- aa=l7_256(i1)%a,cc=l7_256(i1)%c,a1=l7_56%a,c1=l7_56%c,a2=u75
* 6_2(i1)%a,b2=u756_2(i1)%b,c2=u756_2(i1)%c,d2=u756_2(i1)%d,prq=s756,nsu
* m=1
      l7_256(i1)%c(2)=l7_256(i1)%c(2)+l7_56%c(2)*s756*u756_2(i1)
     & %d(1)+l7_56%a(2)*u756_2(i1)%c(2)
      l7_256(i1)%a(2)=l7_256(i1)%a(2)+l7_56%c(2)*s756*u756_2(i1)
     & %b(1)+l7_56%a(2)*u756_2(i1)%a(2)
      end do
  
* II                                                                    
      if (ilept(id3).ne.1) then
* quqd -- p=p72,q=p856
      quqd=p72(0)*p856(0)-p72(1)*p856(1)-p72(2)*p856(2)-p72(3)*p
     & 856(3)
      ccl=zcl(id7)/(f856)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p72,qd=p856,v=cz314(i3,i2)%e,a=u72_314(i3,i2)%a,b=u72_314(i3
* ,i2)%b,c=u72_314(i3,i2)%c,d=u72_314(i3,i2)%d,cl=ccl,nsum=0
      ceps_0=-cz314(i3,i2)%ek0*(p72(2)*p856(3)-p856(2)*p72(3))+p
     & 72k0*(cz314(i3,i2)%e(2)*p856(3)-p856(2)*cz314(i3,i2)%e(3)
     & )-p856k0*(cz314(i3,i2)%e(2)*p72(3)-p72(2)*cz314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i3,i2)%e(3)*p72k0+p72(3)*cz314(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz314(i3,i2)%e(3)*p856k0+p856(3)*cz314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i3,i2)%e(0)*p72(0)-cz314(i3,i2)%e(1)*p72(1)-cz3
     & 14(i3,i2)%e(2)*p72(2)-cz314(i3,i2)%e(3)*p72(3)
      cvqd=cz314(i3,i2)%e(0)*p856(0)-cz314(i3,i2)%e(1)*p856(1)-c
     & z314(i3,i2)%e(2)*p856(2)-cz314(i3,i2)%e(3)*p856(3)
      cauxa=-cz314(i3,i2)%ek0*quqd+p72k0*cvqd+p856k0*cvqu
      cauxb=-cz314(i3,i2)%ek0*p856(2)+p856k0*cz314(i3,i2)%e(2)
      cauxc=+cz314(i3,i2)%ek0*p72(2)-p72k0*cz314(i3,i2)%e(2)
      u72_314(i3,i2)%a(2)=ccl*(cauxa-ceps_0)
      u72_314(i3,i2)%b(1)=ccl*(cauxb-ceps_2)
      u72_314(i3,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u72_314(i3,i2)%d(1)=ccl*cz314(i3,i2)%ek0
      end do
      end do
  
      ccl=fcl(id7)/(f856)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p72,qd=p856,v=cf314(i3,i2)%e,a=u72_314(i3,i2)%a,b=u72_314(i3
* ,i2)%b,c=u72_314(i3,i2)%c,d=u72_314(i3,i2)%d,cl=ccl,nsum=1
      ceps_0=-cf314(i3,i2)%ek0*(p72(2)*p856(3)-p856(2)*p72(3))+p
     & 72k0*(cf314(i3,i2)%e(2)*p856(3)-p856(2)*cf314(i3,i2)%e(3)
     & )-p856k0*(cf314(i3,i2)%e(2)*p72(3)-p72(2)*cf314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i3,i2)%e(3)*p72k0+p72(3)*cf314(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf314(i3,i2)%e(3)*p856k0+p856(3)*cf314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i3,i2)%e(0)*p72(0)-cf314(i3,i2)%e(1)*p72(1)-cf3
     & 14(i3,i2)%e(2)*p72(2)-cf314(i3,i2)%e(3)*p72(3)
      cvqd=cf314(i3,i2)%e(0)*p856(0)-cf314(i3,i2)%e(1)*p856(1)-c
     & f314(i3,i2)%e(2)*p856(2)-cf314(i3,i2)%e(3)*p856(3)
      cauxa=-cf314(i3,i2)%ek0*quqd+p72k0*cvqd+p856k0*cvqu
      cauxb=-cf314(i3,i2)%ek0*p856(2)+p856k0*cf314(i3,i2)%e(2)
      cauxc=+cf314(i3,i2)%ek0*p72(2)-p72k0*cf314(i3,i2)%e(2)
      u72_314(i3,i2)%a(2)=u72_314(i3,i2)%a(2)+ccl*(cauxa-ceps_0)
      u72_314(i3,i2)%b(1)=u72_314(i3,i2)%b(1)+ccl*(cauxb-ceps_2)
      u72_314(i3,i2)%c(2)=u72_314(i3,i2)%c(2)+ccl*(-cauxc+ceps_1
     & )
      u72_314(i3,i2)%d(1)=u72_314(i3,i2)%d(1)+ccl*cf314(i3,i2)%e
     & k0
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_2314(i1,i3,i2)%a,cc=l7_2314(i1,i3,i2)%c,a1=l7_2(i1)%a,
* c1=l7_2(i1)%c,a2=u72_314(i3,i2)%a,b2=u72_314(i3,i2)%b,c2=u72_314(i3,i2
* )%c,d2=u72_314(i3,i2)%d,prq=s72,nsum=0
      l7_2314(i1,i3,i2)%c(2)=l7_2(i1)%c(2)*s72*u72_314(i3,i2)%d(
     & 1)+l7_2(i1)%a(2)*u72_314(i3,i2)%c(2)
      l7_2314(i1,i3,i2)%a(2)=l7_2(i1)%c(2)*s72*u72_314(i3,i2)%b(
     & 1)+l7_2(i1)%a(2)*u72_314(i3,i2)%a(2)
      end do
      end do
      end do
      endif
  
* To use also as u7134_2                                                
* quqd -- p=p7134,q=p856
      quqd=p7134(0)*p856(0)-p7134(1)*p856(1)-p7134(2)*p856(2)-p7
     & 134(3)*p856(3)
      ccl=1.d0/(f856)
      do i1=1,2
* TW0 -- qu=p7134,qd=p856,v=ce2(i1)%e,a=u7314_2(i1)%a,b=u7314_2(i1)%b,c=
* u7314_2(i1)%c,d=u7314_2(i1)%d,cl=ccl,nsum=0
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
      u7314_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u7314_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u7314_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u7314_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      if (ilept(id3).ne.1) then
      do i1=1,2
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_2314(i1,i3,i2)%a,cc=l7_2314(i1,i3,i2)%c,a1=l7_314(i3,i
* 2)%a,c1=l7_314(i3,i2)%c,a2=u7314_2(i1)%a,b2=u7314_2(i1)%b,c2=u7314_2(i
* 1)%c,d2=u7314_2(i1)%d,prq=s7134,nsum=1
      l7_2314(i1,i3,i2)%c(2)=l7_2314(i1,i3,i2)%c(2)+l7_314(i3,i2
     & )%c(2)*s7134*u7314_2(i1)%d(1)+l7_314(i3,i2)%a(2)*u7314_2(
     & i1)%c(2)
      l7_2314(i1,i3,i2)%a(2)=l7_2314(i1,i3,i2)%a(2)+l7_314(i3,i2
     & )%c(2)*s7134*u7314_2(i1)%b(1)+l7_314(i3,i2)%a(2)*u7314_2(
     & i1)%a(2)
      end do
      end do
      end do
  
* III                                                                   
* quqd -- p=p756,q=p82
      quqd=p756(0)*p82(0)-p756(1)*p82(1)-p756(2)*p82(2)-p756(3)*
     & p82(3)
      ccl=zcl(id8)/(f82)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p756,qd=p82,v=cz314(i3,i2)%e,a=u756_314(i3,i2)%a,b=u756_314(
* i3,i2)%b,c=u756_314(i3,i2)%c,d=u756_314(i3,i2)%d,cl=ccl,nsum=0
      ceps_0=-cz314(i3,i2)%ek0*(p756(2)*p82(3)-p82(2)*p756(3))+p
     & 756k0*(cz314(i3,i2)%e(2)*p82(3)-p82(2)*cz314(i3,i2)%e(3))
     & -p82k0*(cz314(i3,i2)%e(2)*p756(3)-p756(2)*cz314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz314(i3,i2)%e(3)*p756k0+p756(3)*cz314(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz314(i3,i2)%e(3)*p82k0+p82(3)*cz314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz314(i3,i2)%e(0)*p756(0)-cz314(i3,i2)%e(1)*p756(1)-c
     & z314(i3,i2)%e(2)*p756(2)-cz314(i3,i2)%e(3)*p756(3)
      cvqd=cz314(i3,i2)%e(0)*p82(0)-cz314(i3,i2)%e(1)*p82(1)-cz3
     & 14(i3,i2)%e(2)*p82(2)-cz314(i3,i2)%e(3)*p82(3)
      cauxa=-cz314(i3,i2)%ek0*quqd+p756k0*cvqd+p82k0*cvqu
      cauxb=-cz314(i3,i2)%ek0*p82(2)+p82k0*cz314(i3,i2)%e(2)
      cauxc=+cz314(i3,i2)%ek0*p756(2)-p756k0*cz314(i3,i2)%e(2)
      u756_314(i3,i2)%a(2)=ccl*(cauxa-ceps_0)
      u756_314(i3,i2)%b(1)=ccl*(cauxb-ceps_2)
      u756_314(i3,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u756_314(i3,i2)%d(1)=ccl*cz314(i3,i2)%ek0
      end do
      end do
  
      ccl=fcl(id8)/(f82)
      do i3=1,2
      do i2=1,2
* TW0 -- qu=p756,qd=p82,v=cf314(i3,i2)%e,a=u756_314(i3,i2)%a,b=u756_314(
* i3,i2)%b,c=u756_314(i3,i2)%c,d=u756_314(i3,i2)%d,cl=ccl,nsum=1
      ceps_0=-cf314(i3,i2)%ek0*(p756(2)*p82(3)-p82(2)*p756(3))+p
     & 756k0*(cf314(i3,i2)%e(2)*p82(3)-p82(2)*cf314(i3,i2)%e(3))
     & -p82k0*(cf314(i3,i2)%e(2)*p756(3)-p756(2)*cf314(i3,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf314(i3,i2)%e(3)*p756k0+p756(3)*cf314(i3,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf314(i3,i2)%e(3)*p82k0+p82(3)*cf314(i3,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf314(i3,i2)%e(0)*p756(0)-cf314(i3,i2)%e(1)*p756(1)-c
     & f314(i3,i2)%e(2)*p756(2)-cf314(i3,i2)%e(3)*p756(3)
      cvqd=cf314(i3,i2)%e(0)*p82(0)-cf314(i3,i2)%e(1)*p82(1)-cf3
     & 14(i3,i2)%e(2)*p82(2)-cf314(i3,i2)%e(3)*p82(3)
      cauxa=-cf314(i3,i2)%ek0*quqd+p756k0*cvqd+p82k0*cvqu
      cauxb=-cf314(i3,i2)%ek0*p82(2)+p82k0*cf314(i3,i2)%e(2)
      cauxc=+cf314(i3,i2)%ek0*p756(2)-p756k0*cf314(i3,i2)%e(2)
      u756_314(i3,i2)%a(2)=u756_314(i3,i2)%a(2)+ccl*(cauxa-ceps_
     & 0)
      u756_314(i3,i2)%b(1)=u756_314(i3,i2)%b(1)+ccl*(cauxb-ceps_
     & 2)
      u756_314(i3,i2)%c(2)=u756_314(i3,i2)%c(2)+ccl*(-cauxc+ceps
     & _1)
      u756_314(i3,i2)%d(1)=u756_314(i3,i2)%d(1)+ccl*cf314(i3,i2)
     & %ek0
      end do
      end do
  
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_56314(i3,i2)%a,cc=l7_56314(i3,i2)%c,a1=l7_56%a,c1=l7_5
* 6%c,a2=u756_314(i3,i2)%a,b2=u756_314(i3,i2)%b,c2=u756_314(i3,i2)%c,d2=
* u756_314(i3,i2)%d,prq=s756,nsum=0
      l7_56314(i3,i2)%c(2)=l7_56%c(2)*s756*u756_314(i3,i2)%d(1)+
     & l7_56%a(2)*u756_314(i3,i2)%c(2)
      l7_56314(i3,i2)%a(2)=l7_56%c(2)*s756*u756_314(i3,i2)%b(1)+
     & l7_56%a(2)*u756_314(i3,i2)%a(2)
      end do
      end do
       endif
  
* To use also as u7134_56                                               
* quqd -- p=p7134,q=p82
      quqd=p7134(0)*p82(0)-p7134(1)*p82(1)-p7134(2)*p82(2)-p7134
     & (3)*p82(3)
      ccl=wcl/(f82)
* TW0 -- qu=p7134,qd=p82,v=cw56%e,a=u7314_56%a,b=u7314_56%b,c=u7314_56%c
* ,d=u7314_56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p7134(2)*p82(3)-p82(2)*p7134(3))+p7134k0
     & *(cw56%e(2)*p82(3)-p82(2)*cw56%e(3))-p82k0*(cw56%e(2)*p71
     & 34(3)-p7134(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p7134k0+p7134(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p82k0+p82(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p7134(0)-cw56%e(1)*p7134(1)-cw56%e(2)*p7134
     & (2)-cw56%e(3)*p7134(3)
      cvqd=cw56%e(0)*p82(0)-cw56%e(1)*p82(1)-cw56%e(2)*p82(2)-cw
     & 56%e(3)*p82(3)
      cauxa=-cw56%ek0*quqd+p7134k0*cvqd+p82k0*cvqu
      cauxb=-cw56%ek0*p82(2)+p82k0*cw56%e(2)
      cauxc=+cw56%ek0*p7134(2)-p7134k0*cw56%e(2)
      u7314_56%a(2)=ccl*(cauxa-ceps_0)
      u7314_56%b(1)=ccl*(cauxb-ceps_2)
      u7314_56%c(2)=ccl*(-cauxc+ceps_1)
      u7314_56%d(1)=ccl*cw56%ek0
  
      if (ilept(id3).ne.1) then
      do i3=1,2
      do i2=1,2
* TLT0_W -- aa=l7_56314(i3,i2)%a,cc=l7_56314(i3,i2)%c,a1=l7_314(i3,i2)%a
* ,c1=l7_314(i3,i2)%c,a2=u7314_56%a,b2=u7314_56%b,c2=u7314_56%c,d2=u7314
* _56%d,prq=s7134,nsum=1
      l7_56314(i3,i2)%c(2)=l7_56314(i3,i2)%c(2)+l7_314(i3,i2)%c(
     & 2)*s7134*u7314_56%d(1)+l7_314(i3,i2)%a(2)*u7314_56%c(2)
      l7_56314(i3,i2)%a(2)=l7_56314(i3,i2)%a(2)+l7_314(i3,i2)%c(
     & 2)*s7134*u7314_56%b(1)+l7_314(i3,i2)%a(2)*u7314_56%a(2)
      end do
      end do
      endif
      endif
  
* Z                                                                     
  
      if (ilept(id3).ne.1) then
        if (.not.(iup(id3).eq.1)) then
* I                                                                     
* quqd -- p=p31,q=p3156
      quqd=p31(0)*p3156(0)-p31(1)*p3156(1)-p31(2)*p3156(2)-p31(3
     & )*p3156(3)
      ccl=wcl/(f3156)
* TW0 -- qu=p31,qd=p3156,v=cw56%e,a=u31_56%a,b=u31_56%b,c=u31_56%c,d=u31
* _56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p31(2)*p3156(3)-p3156(2)*p31(3))+p31k0*(
     & cw56%e(2)*p3156(3)-p3156(2)*cw56%e(3))-p3156k0*(cw56%e(2)
     & *p31(3)-p31(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p31k0+p31(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p3156k0+p3156(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p31(0)-cw56%e(1)*p31(1)-cw56%e(2)*p31(2)-cw
     & 56%e(3)*p31(3)
      cvqd=cw56%e(0)*p3156(0)-cw56%e(1)*p3156(1)-cw56%e(2)*p3156
     & (2)-cw56%e(3)*p3156(3)
      cauxa=-cw56%ek0*quqd+p31k0*cvqd+p3156k0*cvqu
      cauxb=-cw56%ek0*p3156(2)+p3156k0*cw56%e(2)
      cauxc=+cw56%ek0*p31(2)-p31k0*cw56%e(2)
      u31_56%a(2)=ccl*(cauxa-ceps_0)
      u31_56%b(1)=ccl*(cauxb-ceps_2)
      u31_56%c(2)=ccl*(-cauxc+ceps_1)
      u31_56%d(1)=ccl*cw56%ek0
  
      do i1=1,2
* TLT0_W -- aa=l3_156(i1)%a,cc=l3_156(i1)%c,a1=l3_1(i1)%a,c1=l3_1(i1)%c,
* a2=u31_56%a,b2=u31_56%b,c2=u31_56%c,d2=u31_56%d,prq=s31,nsum=0
      l3_156(i1)%c(2)=l3_1(i1)%c(2)*s31*u31_56%d(1)+l3_1(i1)%a(2
     & )*u31_56%c(2)
      l3_156(i1)%a(2)=l3_1(i1)%c(2)*s31*u31_56%b(1)+l3_1(i1)%a(2
     & )*u31_56%a(2)
      end do
  
* quqd -- p=p356,q=p3156
      quqd=p356(0)*p3156(0)-p356(1)*p3156(1)-p356(2)*p3156(2)-p3
     & 56(3)*p3156(3)
      ccl=1.d0/(f3156)
      do i1=1,2
* TW0 -- qu=p356,qd=p3156,v=ce1(i1)%e,a=u356_1(i1)%a,b=u356_1(i1)%b,c=u3
* 56_1(i1)%c,d=u356_1(i1)%d,cl=ccl,nsum=0
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
      u356_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u356_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u356_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u356_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      do i1=1,2
* TLT0_W -- aa=l3_156(i1)%a,cc=l3_156(i1)%c,a1=l3_56%a,c1=l3_56%c,a2=u35
* 6_1(i1)%a,b2=u356_1(i1)%b,c2=u356_1(i1)%c,d2=u356_1(i1)%d,prq=s356,nsu
* m=1
      l3_156(i1)%c(2)=l3_156(i1)%c(2)+l3_56%c(2)*s356*u356_1(i1)
     & %d(1)+l3_56%a(2)*u356_1(i1)%c(2)
      l3_156(i1)%a(2)=l3_156(i1)%a(2)+l3_56%c(2)*s356*u356_1(i1)
     & %b(1)+l3_56%a(2)*u356_1(i1)%a(2)
      end do
  
        else
  
* II                                                                    
      if (ilept(id7).ne.1) then
* quqd -- p=p31,q=p456
      quqd=p31(0)*p456(0)-p31(1)*p456(1)-p31(2)*p456(2)-p31(3)*p
     & 456(3)
      ccl=wcl/(f456)
      do i2=1,2
* TW0 -- qu=p31,qd=p456,v=cw728(i2)%e,a=u31_728(i2)%a,b=u31_728(i2)%b,c=
* u31_728(i2)%c,d=u31_728(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw728(i2)%ek0*(p31(2)*p456(3)-p456(2)*p31(3))+p31k
     & 0*(cw728(i2)%e(2)*p456(3)-p456(2)*cw728(i2)%e(3))-p456k0*
     & (cw728(i2)%e(2)*p31(3)-p31(2)*cw728(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw728(i2)%e(3)*p31k0+p31(3)*cw728(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw728(i2)%e(3)*p456k0+p456(3)*cw728(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw728(i2)%e(0)*p31(0)-cw728(i2)%e(1)*p31(1)-cw728(i2)
     & %e(2)*p31(2)-cw728(i2)%e(3)*p31(3)
      cvqd=cw728(i2)%e(0)*p456(0)-cw728(i2)%e(1)*p456(1)-cw728(i
     & 2)%e(2)*p456(2)-cw728(i2)%e(3)*p456(3)
      cauxa=-cw728(i2)%ek0*quqd+p31k0*cvqd+p456k0*cvqu
      cauxb=-cw728(i2)%ek0*p456(2)+p456k0*cw728(i2)%e(2)
      cauxc=+cw728(i2)%ek0*p31(2)-p31k0*cw728(i2)%e(2)
      u31_728(i2)%a(2)=ccl*(cauxa-ceps_0)
      u31_728(i2)%b(1)=ccl*(cauxb-ceps_2)
      u31_728(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u31_728(i2)%d(1)=ccl*cw728(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l3_1728(i1,i2)%a,cc=l3_1728(i1,i2)%c,a1=l3_1(i1)%a,c1=l3_
* 1(i1)%c,a2=u31_728(i2)%a,b2=u31_728(i2)%b,c2=u31_728(i2)%c,d2=u31_728(
* i2)%d,prq=s31,nsum=0
      l3_1728(i1,i2)%c(2)=l3_1(i1)%c(2)*s31*u31_728(i2)%d(1)+l3_
     & 1(i1)%a(2)*u31_728(i2)%c(2)
      l3_1728(i1,i2)%a(2)=l3_1(i1)%c(2)*s31*u31_728(i2)%b(1)+l3_
     & 1(i1)%a(2)*u31_728(i2)%a(2)
      end do
      end do
      endif
  
* To use also as u3278_1                                                
* quqd -- p=p3278,q=p456
      quqd=p3278(0)*p456(0)-p3278(1)*p456(1)-p3278(2)*p456(2)-p3
     & 278(3)*p456(3)
      ccl=1.d0/(f456)
      do i1=1,2
* TW0 -- qu=p3278,qd=p456,v=ce1(i1)%e,a=u3728_1(i1)%a,b=u3728_1(i1)%b,c=
* u3728_1(i1)%c,d=u3728_1(i1)%d,cl=ccl,nsum=0
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
      u3728_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u3728_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u3728_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u3728_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l3_1728(i1,i2)%a,cc=l3_1728(i1,i2)%c,a1=l3_728(i2)%a,c1=l
* 3_728(i2)%c,a2=u3728_1(i1)%a,b2=u3728_1(i1)%b,c2=u3728_1(i1)%c,d2=u372
* 8_1(i1)%d,prq=s3278,nsum=1
      l3_1728(i1,i2)%c(2)=l3_1728(i1,i2)%c(2)+l3_728(i2)%c(2)*s3
     & 278*u3728_1(i1)%d(1)+l3_728(i2)%a(2)*u3728_1(i1)%c(2)
      l3_1728(i1,i2)%a(2)=l3_1728(i1,i2)%a(2)+l3_728(i2)%c(2)*s3
     & 278*u3728_1(i1)%b(1)+l3_728(i2)%a(2)*u3728_1(i1)%a(2)
      end do
      end do
      endif
        endif !iup(id3)
  
* III                                                                   
       if (.not.(iup(id3).eq.1)) then
         if (ilept(id7).ne.1) then
* quqd -- p=p356,q=p41
      quqd=p356(0)*p41(0)-p356(1)*p41(1)-p356(2)*p41(2)-p356(3)*
     & p41(3)
      ccl=wcl/(f41)
      do i2=1,2
* TW0 -- qu=p356,qd=p41,v=cw728(i2)%e,a=u356_728(i2)%a,b=u356_728(i2)%b,
* c=u356_728(i2)%c,d=u356_728(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw728(i2)%ek0*(p356(2)*p41(3)-p41(2)*p356(3))+p356
     & k0*(cw728(i2)%e(2)*p41(3)-p41(2)*cw728(i2)%e(3))-p41k0*(c
     & w728(i2)%e(2)*p356(3)-p356(2)*cw728(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw728(i2)%e(3)*p356k0+p356(3)*cw728(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw728(i2)%e(3)*p41k0+p41(3)*cw728(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw728(i2)%e(0)*p356(0)-cw728(i2)%e(1)*p356(1)-cw728(i
     & 2)%e(2)*p356(2)-cw728(i2)%e(3)*p356(3)
      cvqd=cw728(i2)%e(0)*p41(0)-cw728(i2)%e(1)*p41(1)-cw728(i2)
     & %e(2)*p41(2)-cw728(i2)%e(3)*p41(3)
      cauxa=-cw728(i2)%ek0*quqd+p356k0*cvqd+p41k0*cvqu
      cauxb=-cw728(i2)%ek0*p41(2)+p41k0*cw728(i2)%e(2)
      cauxc=+cw728(i2)%ek0*p356(2)-p356k0*cw728(i2)%e(2)
      u356_728(i2)%a(2)=ccl*(cauxa-ceps_0)
      u356_728(i2)%b(1)=ccl*(cauxb-ceps_2)
      u356_728(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u356_728(i2)%d(1)=ccl*cw728(i2)%ek0
      end do
  
      do i2=1,2
* TLT0_W -- aa=l3_56728(i2)%a,cc=l3_56728(i2)%c,a1=l3_56%a,c1=l3_56%c,a2
* =u356_728(i2)%a,b2=u356_728(i2)%b,c2=u356_728(i2)%c,d2=u356_728(i2)%d,
* prq=s356,nsum=0
      l3_56728(i2)%c(2)=l3_56%c(2)*s356*u356_728(i2)%d(1)+l3_56%
     & a(2)*u356_728(i2)%c(2)
      l3_56728(i2)%a(2)=l3_56%c(2)*s356*u356_728(i2)%b(1)+l3_56%
     & a(2)*u356_728(i2)%a(2)
      end do
          endif
  
       else
* To use also as u3278_56                                               
* quqd -- p=p3278,q=p41
      quqd=p3278(0)*p41(0)-p3278(1)*p41(1)-p3278(2)*p41(2)-p3278
     & (3)*p41(3)
      ccl=wcl/(f41)
* TW0 -- qu=p3278,qd=p41,v=cw56%e,a=u3728_56%a,b=u3728_56%b,c=u3728_56%c
* ,d=u3728_56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p3278(2)*p41(3)-p41(2)*p3278(3))+p3278k0
     & *(cw56%e(2)*p41(3)-p41(2)*cw56%e(3))-p41k0*(cw56%e(2)*p32
     & 78(3)-p3278(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p3278k0+p3278(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p41k0+p41(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p3278(0)-cw56%e(1)*p3278(1)-cw56%e(2)*p3278
     & (2)-cw56%e(3)*p3278(3)
      cvqd=cw56%e(0)*p41(0)-cw56%e(1)*p41(1)-cw56%e(2)*p41(2)-cw
     & 56%e(3)*p41(3)
      cauxa=-cw56%ek0*quqd+p3278k0*cvqd+p41k0*cvqu
      cauxb=-cw56%ek0*p41(2)+p41k0*cw56%e(2)
      cauxc=+cw56%ek0*p3278(2)-p3278k0*cw56%e(2)
      u3728_56%a(2)=ccl*(cauxa-ceps_0)
      u3728_56%b(1)=ccl*(cauxb-ceps_2)
      u3728_56%c(2)=ccl*(-cauxc+ceps_1)
      u3728_56%d(1)=ccl*cw56%ek0
  
          if (ilept(id7).ne.1) then
      do i2=1,2
* TLT0_W -- aa=l3_56728(i2)%a,cc=l3_56728(i2)%c,a1=l3_728(i2)%a,c1=l3_72
* 8(i2)%c,a2=u3728_56%a,b2=u3728_56%b,c2=u3728_56%c,d2=u3728_56%d,prq=s3
* 278,nsum=0
      l3_56728(i2)%c(2)=l3_728(i2)%c(2)*s3278*u3728_56%d(1)+l3_7
     & 28(i2)%a(2)*u3728_56%c(2)
      l3_56728(i2)%a(2)=l3_728(i2)%c(2)*s3278*u3728_56%b(1)+l3_7
     & 28(i2)%a(2)*u3728_56%a(2)
      end do
          endif
  
       endif !iup(id3)
  
      endif
  
  
      if (ilept(id3).ne.1) then
        if (.not.(iup(id3).eq.1)) then
* I                                                                     
* quqd -- p=p32,q=p3256
      quqd=p32(0)*p3256(0)-p32(1)*p3256(1)-p32(2)*p3256(2)-p32(3
     & )*p3256(3)
      ccl=wcl/(f3256)
* TW0 -- qu=p32,qd=p3256,v=cw56%e,a=u32_56%a,b=u32_56%b,c=u32_56%c,d=u32
* _56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p32(2)*p3256(3)-p3256(2)*p32(3))+p32k0*(
     & cw56%e(2)*p3256(3)-p3256(2)*cw56%e(3))-p3256k0*(cw56%e(2)
     & *p32(3)-p32(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p32k0+p32(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p3256k0+p3256(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p32(0)-cw56%e(1)*p32(1)-cw56%e(2)*p32(2)-cw
     & 56%e(3)*p32(3)
      cvqd=cw56%e(0)*p3256(0)-cw56%e(1)*p3256(1)-cw56%e(2)*p3256
     & (2)-cw56%e(3)*p3256(3)
      cauxa=-cw56%ek0*quqd+p32k0*cvqd+p3256k0*cvqu
      cauxb=-cw56%ek0*p3256(2)+p3256k0*cw56%e(2)
      cauxc=+cw56%ek0*p32(2)-p32k0*cw56%e(2)
      u32_56%a(2)=ccl*(cauxa-ceps_0)
      u32_56%b(1)=ccl*(cauxb-ceps_2)
      u32_56%c(2)=ccl*(-cauxc+ceps_1)
      u32_56%d(1)=ccl*cw56%ek0
  
      do i1=1,2
* TLT0_W -- aa=l3_256(i1)%a,cc=l3_256(i1)%c,a1=l3_2(i1)%a,c1=l3_2(i1)%c,
* a2=u32_56%a,b2=u32_56%b,c2=u32_56%c,d2=u32_56%d,prq=s32,nsum=0
      l3_256(i1)%c(2)=l3_2(i1)%c(2)*s32*u32_56%d(1)+l3_2(i1)%a(2
     & )*u32_56%c(2)
      l3_256(i1)%a(2)=l3_2(i1)%c(2)*s32*u32_56%b(1)+l3_2(i1)%a(2
     & )*u32_56%a(2)
      end do
  
* quqd -- p=p356,q=p3256
      quqd=p356(0)*p3256(0)-p356(1)*p3256(1)-p356(2)*p3256(2)-p3
     & 56(3)*p3256(3)
      ccl=1.d0/(f3256)
      do i1=1,2
* TW0 -- qu=p356,qd=p3256,v=ce2(i1)%e,a=u356_2(i1)%a,b=u356_2(i1)%b,c=u3
* 56_2(i1)%c,d=u356_2(i1)%d,cl=ccl,nsum=0
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
      u356_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u356_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u356_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u356_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      do i1=1,2
* TLT0_W -- aa=l3_256(i1)%a,cc=l3_256(i1)%c,a1=l3_56%a,c1=l3_56%c,a2=u35
* 6_2(i1)%a,b2=u356_2(i1)%b,c2=u356_2(i1)%c,d2=u356_2(i1)%d,prq=s356,nsu
* m=1
      l3_256(i1)%c(2)=l3_256(i1)%c(2)+l3_56%c(2)*s356*u356_2(i1)
     & %d(1)+l3_56%a(2)*u356_2(i1)%c(2)
      l3_256(i1)%a(2)=l3_256(i1)%a(2)+l3_56%c(2)*s356*u356_2(i1)
     & %b(1)+l3_56%a(2)*u356_2(i1)%a(2)
      end do
  
        else
  
* II                                                                    
      if (ilept(id7).ne.1) then
* quqd -- p=p32,q=p456
      quqd=p32(0)*p456(0)-p32(1)*p456(1)-p32(2)*p456(2)-p32(3)*p
     & 456(3)
      ccl=wcl/(f456)
      do i2=1,2
* TW0 -- qu=p32,qd=p456,v=cw718(i2)%e,a=u32_718(i2)%a,b=u32_718(i2)%b,c=
* u32_718(i2)%c,d=u32_718(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw718(i2)%ek0*(p32(2)*p456(3)-p456(2)*p32(3))+p32k
     & 0*(cw718(i2)%e(2)*p456(3)-p456(2)*cw718(i2)%e(3))-p456k0*
     & (cw718(i2)%e(2)*p32(3)-p32(2)*cw718(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw718(i2)%e(3)*p32k0+p32(3)*cw718(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw718(i2)%e(3)*p456k0+p456(3)*cw718(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw718(i2)%e(0)*p32(0)-cw718(i2)%e(1)*p32(1)-cw718(i2)
     & %e(2)*p32(2)-cw718(i2)%e(3)*p32(3)
      cvqd=cw718(i2)%e(0)*p456(0)-cw718(i2)%e(1)*p456(1)-cw718(i
     & 2)%e(2)*p456(2)-cw718(i2)%e(3)*p456(3)
      cauxa=-cw718(i2)%ek0*quqd+p32k0*cvqd+p456k0*cvqu
      cauxb=-cw718(i2)%ek0*p456(2)+p456k0*cw718(i2)%e(2)
      cauxc=+cw718(i2)%ek0*p32(2)-p32k0*cw718(i2)%e(2)
      u32_718(i2)%a(2)=ccl*(cauxa-ceps_0)
      u32_718(i2)%b(1)=ccl*(cauxb-ceps_2)
      u32_718(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u32_718(i2)%d(1)=ccl*cw718(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l3_2718(i1,i2)%a,cc=l3_2718(i1,i2)%c,a1=l3_2(i1)%a,c1=l3_
* 2(i1)%c,a2=u32_718(i2)%a,b2=u32_718(i2)%b,c2=u32_718(i2)%c,d2=u32_718(
* i2)%d,prq=s32,nsum=0
      l3_2718(i1,i2)%c(2)=l3_2(i1)%c(2)*s32*u32_718(i2)%d(1)+l3_
     & 2(i1)%a(2)*u32_718(i2)%c(2)
      l3_2718(i1,i2)%a(2)=l3_2(i1)%c(2)*s32*u32_718(i2)%b(1)+l3_
     & 2(i1)%a(2)*u32_718(i2)%a(2)
      end do
      end do
      endif
  
* To use also as u3178_2                                                
* quqd -- p=p3178,q=p456
      quqd=p3178(0)*p456(0)-p3178(1)*p456(1)-p3178(2)*p456(2)-p3
     & 178(3)*p456(3)
      ccl=1.d0/(f456)
      do i1=1,2
* TW0 -- qu=p3178,qd=p456,v=ce2(i1)%e,a=u3718_2(i1)%a,b=u3718_2(i1)%b,c=
* u3718_2(i1)%c,d=u3718_2(i1)%d,cl=ccl,nsum=0
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
      u3718_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u3718_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u3718_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u3718_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      if (ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l3_2718(i1,i2)%a,cc=l3_2718(i1,i2)%c,a1=l3_718(i2)%a,c1=l
* 3_718(i2)%c,a2=u3718_2(i1)%a,b2=u3718_2(i1)%b,c2=u3718_2(i1)%c,d2=u371
* 8_2(i1)%d,prq=s3178,nsum=1
      l3_2718(i1,i2)%c(2)=l3_2718(i1,i2)%c(2)+l3_718(i2)%c(2)*s3
     & 178*u3718_2(i1)%d(1)+l3_718(i2)%a(2)*u3718_2(i1)%c(2)
      l3_2718(i1,i2)%a(2)=l3_2718(i1,i2)%a(2)+l3_718(i2)%c(2)*s3
     & 178*u3718_2(i1)%b(1)+l3_718(i2)%a(2)*u3718_2(i1)%a(2)
      end do
      end do
      endif
        endif !iup(id3)
  
* III                                                                   
       if (.not.(iup(id3).eq.1)) then
         if (ilept(id7).ne.1) then
* quqd -- p=p356,q=p42
      quqd=p356(0)*p42(0)-p356(1)*p42(1)-p356(2)*p42(2)-p356(3)*
     & p42(3)
      ccl=wcl/(f42)
      do i2=1,2
* TW0 -- qu=p356,qd=p42,v=cw718(i2)%e,a=u356_718(i2)%a,b=u356_718(i2)%b,
* c=u356_718(i2)%c,d=u356_718(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw718(i2)%ek0*(p356(2)*p42(3)-p42(2)*p356(3))+p356
     & k0*(cw718(i2)%e(2)*p42(3)-p42(2)*cw718(i2)%e(3))-p42k0*(c
     & w718(i2)%e(2)*p356(3)-p356(2)*cw718(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw718(i2)%e(3)*p356k0+p356(3)*cw718(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw718(i2)%e(3)*p42k0+p42(3)*cw718(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw718(i2)%e(0)*p356(0)-cw718(i2)%e(1)*p356(1)-cw718(i
     & 2)%e(2)*p356(2)-cw718(i2)%e(3)*p356(3)
      cvqd=cw718(i2)%e(0)*p42(0)-cw718(i2)%e(1)*p42(1)-cw718(i2)
     & %e(2)*p42(2)-cw718(i2)%e(3)*p42(3)
      cauxa=-cw718(i2)%ek0*quqd+p356k0*cvqd+p42k0*cvqu
      cauxb=-cw718(i2)%ek0*p42(2)+p42k0*cw718(i2)%e(2)
      cauxc=+cw718(i2)%ek0*p356(2)-p356k0*cw718(i2)%e(2)
      u356_718(i2)%a(2)=ccl*(cauxa-ceps_0)
      u356_718(i2)%b(1)=ccl*(cauxb-ceps_2)
      u356_718(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u356_718(i2)%d(1)=ccl*cw718(i2)%ek0
      end do
  
      do i2=1,2
* TLT0_W -- aa=l3_56718(i2)%a,cc=l3_56718(i2)%c,a1=l3_56%a,c1=l3_56%c,a2
* =u356_718(i2)%a,b2=u356_718(i2)%b,c2=u356_718(i2)%c,d2=u356_718(i2)%d,
* prq=s356,nsum=0
      l3_56718(i2)%c(2)=l3_56%c(2)*s356*u356_718(i2)%d(1)+l3_56%
     & a(2)*u356_718(i2)%c(2)
      l3_56718(i2)%a(2)=l3_56%c(2)*s356*u356_718(i2)%b(1)+l3_56%
     & a(2)*u356_718(i2)%a(2)
      end do
          endif
  
       else
* To use also as u3178_56                                               
* quqd -- p=p3178,q=p42
      quqd=p3178(0)*p42(0)-p3178(1)*p42(1)-p3178(2)*p42(2)-p3178
     & (3)*p42(3)
      ccl=wcl/(f42)
* TW0 -- qu=p3178,qd=p42,v=cw56%e,a=u3718_56%a,b=u3718_56%b,c=u3718_56%c
* ,d=u3718_56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p3178(2)*p42(3)-p42(2)*p3178(3))+p3178k0
     & *(cw56%e(2)*p42(3)-p42(2)*cw56%e(3))-p42k0*(cw56%e(2)*p31
     & 78(3)-p3178(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p3178k0+p3178(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p42k0+p42(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p3178(0)-cw56%e(1)*p3178(1)-cw56%e(2)*p3178
     & (2)-cw56%e(3)*p3178(3)
      cvqd=cw56%e(0)*p42(0)-cw56%e(1)*p42(1)-cw56%e(2)*p42(2)-cw
     & 56%e(3)*p42(3)
      cauxa=-cw56%ek0*quqd+p3178k0*cvqd+p42k0*cvqu
      cauxb=-cw56%ek0*p42(2)+p42k0*cw56%e(2)
      cauxc=+cw56%ek0*p3178(2)-p3178k0*cw56%e(2)
      u3718_56%a(2)=ccl*(cauxa-ceps_0)
      u3718_56%b(1)=ccl*(cauxb-ceps_2)
      u3718_56%c(2)=ccl*(-cauxc+ceps_1)
      u3718_56%d(1)=ccl*cw56%ek0
  
          if (ilept(id7).ne.1) then
      do i2=1,2
* TLT0_W -- aa=l3_56718(i2)%a,cc=l3_56718(i2)%c,a1=l3_718(i2)%a,c1=l3_71
* 8(i2)%c,a2=u3718_56%a,b2=u3718_56%b,c2=u3718_56%c,d2=u3718_56%d,prq=s3
* 178,nsum=0
      l3_56718(i2)%c(2)=l3_718(i2)%c(2)*s3178*u3718_56%d(1)+l3_7
     & 18(i2)%a(2)*u3718_56%c(2)
      l3_56718(i2)%a(2)=l3_718(i2)%c(2)*s3178*u3718_56%b(1)+l3_7
     & 18(i2)%a(2)*u3718_56%a(2)
      end do
          endif
  
       endif !iup(id3)
  
      endif
  
  
      if (ilept(id3).ne.1) then
        if ((iup(id3).eq.1)) then
* I                                                                     
* quqd -- p=p31,q=p3178
      quqd=p31(0)*p3178(0)-p31(1)*p3178(1)-p31(2)*p3178(2)-p31(3
     & )*p3178(3)
      ccl=wcl/(f3178)
* TW0 -- qu=p31,qd=p3178,v=cw78%e,a=u31_78%a,b=u31_78%b,c=u31_78%c,d=u31
* _78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p31(2)*p3178(3)-p3178(2)*p31(3))+p31k0*(
     & cw78%e(2)*p3178(3)-p3178(2)*cw78%e(3))-p3178k0*(cw78%e(2)
     & *p31(3)-p31(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p31k0+p31(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p3178k0+p3178(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p31(0)-cw78%e(1)*p31(1)-cw78%e(2)*p31(2)-cw
     & 78%e(3)*p31(3)
      cvqd=cw78%e(0)*p3178(0)-cw78%e(1)*p3178(1)-cw78%e(2)*p3178
     & (2)-cw78%e(3)*p3178(3)
      cauxa=-cw78%ek0*quqd+p31k0*cvqd+p3178k0*cvqu
      cauxb=-cw78%ek0*p3178(2)+p3178k0*cw78%e(2)
      cauxc=+cw78%ek0*p31(2)-p31k0*cw78%e(2)
      u31_78%a(2)=ccl*(cauxa-ceps_0)
      u31_78%b(1)=ccl*(cauxb-ceps_2)
      u31_78%c(2)=ccl*(-cauxc+ceps_1)
      u31_78%d(1)=ccl*cw78%ek0
  
      do i1=1,2
* TLT0_W -- aa=l3_178(i1)%a,cc=l3_178(i1)%c,a1=l3_1(i1)%a,c1=l3_1(i1)%c,
* a2=u31_78%a,b2=u31_78%b,c2=u31_78%c,d2=u31_78%d,prq=s31,nsum=0
      l3_178(i1)%c(2)=l3_1(i1)%c(2)*s31*u31_78%d(1)+l3_1(i1)%a(2
     & )*u31_78%c(2)
      l3_178(i1)%a(2)=l3_1(i1)%c(2)*s31*u31_78%b(1)+l3_1(i1)%a(2
     & )*u31_78%a(2)
      end do
  
* quqd -- p=p378,q=p3178
      quqd=p378(0)*p3178(0)-p378(1)*p3178(1)-p378(2)*p3178(2)-p3
     & 78(3)*p3178(3)
      ccl=1.d0/(f3178)
      do i1=1,2
* TW0 -- qu=p378,qd=p3178,v=ce1(i1)%e,a=u378_1(i1)%a,b=u378_1(i1)%b,c=u3
* 78_1(i1)%c,d=u378_1(i1)%d,cl=ccl,nsum=0
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
      u378_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u378_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u378_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u378_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      do i1=1,2
* TLT0_W -- aa=l3_178(i1)%a,cc=l3_178(i1)%c,a1=l3_78%a,c1=l3_78%c,a2=u37
* 8_1(i1)%a,b2=u378_1(i1)%b,c2=u378_1(i1)%c,d2=u378_1(i1)%d,prq=s378,nsu
* m=1
      l3_178(i1)%c(2)=l3_178(i1)%c(2)+l3_78%c(2)*s378*u378_1(i1)
     & %d(1)+l3_78%a(2)*u378_1(i1)%c(2)
      l3_178(i1)%a(2)=l3_178(i1)%a(2)+l3_78%c(2)*s378*u378_1(i1)
     & %b(1)+l3_78%a(2)*u378_1(i1)%a(2)
      end do
  
        else
  
* II                                                                    
      if (ilept(id5).ne.1) then
* quqd -- p=p31,q=p478
      quqd=p31(0)*p478(0)-p31(1)*p478(1)-p31(2)*p478(2)-p31(3)*p
     & 478(3)
      ccl=wcl/(f478)
      do i2=1,2
* TW0 -- qu=p31,qd=p478,v=cw526(i2)%e,a=u31_526(i2)%a,b=u31_526(i2)%b,c=
* u31_526(i2)%c,d=u31_526(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw526(i2)%ek0*(p31(2)*p478(3)-p478(2)*p31(3))+p31k
     & 0*(cw526(i2)%e(2)*p478(3)-p478(2)*cw526(i2)%e(3))-p478k0*
     & (cw526(i2)%e(2)*p31(3)-p31(2)*cw526(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw526(i2)%e(3)*p31k0+p31(3)*cw526(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw526(i2)%e(3)*p478k0+p478(3)*cw526(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw526(i2)%e(0)*p31(0)-cw526(i2)%e(1)*p31(1)-cw526(i2)
     & %e(2)*p31(2)-cw526(i2)%e(3)*p31(3)
      cvqd=cw526(i2)%e(0)*p478(0)-cw526(i2)%e(1)*p478(1)-cw526(i
     & 2)%e(2)*p478(2)-cw526(i2)%e(3)*p478(3)
      cauxa=-cw526(i2)%ek0*quqd+p31k0*cvqd+p478k0*cvqu
      cauxb=-cw526(i2)%ek0*p478(2)+p478k0*cw526(i2)%e(2)
      cauxc=+cw526(i2)%ek0*p31(2)-p31k0*cw526(i2)%e(2)
      u31_526(i2)%a(2)=ccl*(cauxa-ceps_0)
      u31_526(i2)%b(1)=ccl*(cauxb-ceps_2)
      u31_526(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u31_526(i2)%d(1)=ccl*cw526(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l3_1526(i1,i2)%a,cc=l3_1526(i1,i2)%c,a1=l3_1(i1)%a,c1=l3_
* 1(i1)%c,a2=u31_526(i2)%a,b2=u31_526(i2)%b,c2=u31_526(i2)%c,d2=u31_526(
* i2)%d,prq=s31,nsum=0
      l3_1526(i1,i2)%c(2)=l3_1(i1)%c(2)*s31*u31_526(i2)%d(1)+l3_
     & 1(i1)%a(2)*u31_526(i2)%c(2)
      l3_1526(i1,i2)%a(2)=l3_1(i1)%c(2)*s31*u31_526(i2)%b(1)+l3_
     & 1(i1)%a(2)*u31_526(i2)%a(2)
      end do
      end do
      endif
  
* To use also as u3256_1                                                
* quqd -- p=p3256,q=p478
      quqd=p3256(0)*p478(0)-p3256(1)*p478(1)-p3256(2)*p478(2)-p3
     & 256(3)*p478(3)
      ccl=1.d0/(f478)
      do i1=1,2
* TW0 -- qu=p3256,qd=p478,v=ce1(i1)%e,a=u3526_1(i1)%a,b=u3526_1(i1)%b,c=
* u3526_1(i1)%c,d=u3526_1(i1)%d,cl=ccl,nsum=0
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
      u3526_1(i1)%a(2)=ccl*(cauxa-ceps_0)
      u3526_1(i1)%b(1)=ccl*(cauxb-ceps_2)
      u3526_1(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u3526_1(i1)%d(1)=ccl*ce1(i1)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l3_1526(i1,i2)%a,cc=l3_1526(i1,i2)%c,a1=l3_526(i2)%a,c1=l
* 3_526(i2)%c,a2=u3526_1(i1)%a,b2=u3526_1(i1)%b,c2=u3526_1(i1)%c,d2=u352
* 6_1(i1)%d,prq=s3256,nsum=1
      l3_1526(i1,i2)%c(2)=l3_1526(i1,i2)%c(2)+l3_526(i2)%c(2)*s3
     & 256*u3526_1(i1)%d(1)+l3_526(i2)%a(2)*u3526_1(i1)%c(2)
      l3_1526(i1,i2)%a(2)=l3_1526(i1,i2)%a(2)+l3_526(i2)%c(2)*s3
     & 256*u3526_1(i1)%b(1)+l3_526(i2)%a(2)*u3526_1(i1)%a(2)
      end do
      end do
      endif
        endif !iup(id3)
  
* III                                                                   
       if ((iup(id3).eq.1)) then
         if (ilept(id5).ne.1) then
* quqd -- p=p378,q=p41
      quqd=p378(0)*p41(0)-p378(1)*p41(1)-p378(2)*p41(2)-p378(3)*
     & p41(3)
      ccl=wcl/(f41)
      do i2=1,2
* TW0 -- qu=p378,qd=p41,v=cw526(i2)%e,a=u378_526(i2)%a,b=u378_526(i2)%b,
* c=u378_526(i2)%c,d=u378_526(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw526(i2)%ek0*(p378(2)*p41(3)-p41(2)*p378(3))+p378
     & k0*(cw526(i2)%e(2)*p41(3)-p41(2)*cw526(i2)%e(3))-p41k0*(c
     & w526(i2)%e(2)*p378(3)-p378(2)*cw526(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw526(i2)%e(3)*p378k0+p378(3)*cw526(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw526(i2)%e(3)*p41k0+p41(3)*cw526(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw526(i2)%e(0)*p378(0)-cw526(i2)%e(1)*p378(1)-cw526(i
     & 2)%e(2)*p378(2)-cw526(i2)%e(3)*p378(3)
      cvqd=cw526(i2)%e(0)*p41(0)-cw526(i2)%e(1)*p41(1)-cw526(i2)
     & %e(2)*p41(2)-cw526(i2)%e(3)*p41(3)
      cauxa=-cw526(i2)%ek0*quqd+p378k0*cvqd+p41k0*cvqu
      cauxb=-cw526(i2)%ek0*p41(2)+p41k0*cw526(i2)%e(2)
      cauxc=+cw526(i2)%ek0*p378(2)-p378k0*cw526(i2)%e(2)
      u378_526(i2)%a(2)=ccl*(cauxa-ceps_0)
      u378_526(i2)%b(1)=ccl*(cauxb-ceps_2)
      u378_526(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u378_526(i2)%d(1)=ccl*cw526(i2)%ek0
      end do
  
      do i2=1,2
* TLT0_W -- aa=l3_78526(i2)%a,cc=l3_78526(i2)%c,a1=l3_78%a,c1=l3_78%c,a2
* =u378_526(i2)%a,b2=u378_526(i2)%b,c2=u378_526(i2)%c,d2=u378_526(i2)%d,
* prq=s378,nsum=0
      l3_78526(i2)%c(2)=l3_78%c(2)*s378*u378_526(i2)%d(1)+l3_78%
     & a(2)*u378_526(i2)%c(2)
      l3_78526(i2)%a(2)=l3_78%c(2)*s378*u378_526(i2)%b(1)+l3_78%
     & a(2)*u378_526(i2)%a(2)
      end do
          endif
  
       else
* To use also as u3256_78                                               
* quqd -- p=p3256,q=p41
      quqd=p3256(0)*p41(0)-p3256(1)*p41(1)-p3256(2)*p41(2)-p3256
     & (3)*p41(3)
      ccl=wcl/(f41)
* TW0 -- qu=p3256,qd=p41,v=cw78%e,a=u3526_78%a,b=u3526_78%b,c=u3526_78%c
* ,d=u3526_78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p3256(2)*p41(3)-p41(2)*p3256(3))+p3256k0
     & *(cw78%e(2)*p41(3)-p41(2)*cw78%e(3))-p41k0*(cw78%e(2)*p32
     & 56(3)-p3256(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p3256k0+p3256(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p41k0+p41(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p3256(0)-cw78%e(1)*p3256(1)-cw78%e(2)*p3256
     & (2)-cw78%e(3)*p3256(3)
      cvqd=cw78%e(0)*p41(0)-cw78%e(1)*p41(1)-cw78%e(2)*p41(2)-cw
     & 78%e(3)*p41(3)
      cauxa=-cw78%ek0*quqd+p3256k0*cvqd+p41k0*cvqu
      cauxb=-cw78%ek0*p41(2)+p41k0*cw78%e(2)
      cauxc=+cw78%ek0*p3256(2)-p3256k0*cw78%e(2)
      u3526_78%a(2)=ccl*(cauxa-ceps_0)
      u3526_78%b(1)=ccl*(cauxb-ceps_2)
      u3526_78%c(2)=ccl*(-cauxc+ceps_1)
      u3526_78%d(1)=ccl*cw78%ek0
  
          if (ilept(id5).ne.1) then
      do i2=1,2
* TLT0_W -- aa=l3_78526(i2)%a,cc=l3_78526(i2)%c,a1=l3_526(i2)%a,c1=l3_52
* 6(i2)%c,a2=u3526_78%a,b2=u3526_78%b,c2=u3526_78%c,d2=u3526_78%d,prq=s3
* 256,nsum=0
      l3_78526(i2)%c(2)=l3_526(i2)%c(2)*s3256*u3526_78%d(1)+l3_5
     & 26(i2)%a(2)*u3526_78%c(2)
      l3_78526(i2)%a(2)=l3_526(i2)%c(2)*s3256*u3526_78%b(1)+l3_5
     & 26(i2)%a(2)*u3526_78%a(2)
      end do
          endif
  
       endif !iup(id3)
  
      endif
  
  
      if (ilept(id3).ne.1) then
        if ((iup(id3).eq.1)) then
* I                                                                     
* quqd -- p=p32,q=p3278
      quqd=p32(0)*p3278(0)-p32(1)*p3278(1)-p32(2)*p3278(2)-p32(3
     & )*p3278(3)
      ccl=wcl/(f3278)
* TW0 -- qu=p32,qd=p3278,v=cw78%e,a=u32_78%a,b=u32_78%b,c=u32_78%c,d=u32
* _78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p32(2)*p3278(3)-p3278(2)*p32(3))+p32k0*(
     & cw78%e(2)*p3278(3)-p3278(2)*cw78%e(3))-p3278k0*(cw78%e(2)
     & *p32(3)-p32(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p32k0+p32(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p3278k0+p3278(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p32(0)-cw78%e(1)*p32(1)-cw78%e(2)*p32(2)-cw
     & 78%e(3)*p32(3)
      cvqd=cw78%e(0)*p3278(0)-cw78%e(1)*p3278(1)-cw78%e(2)*p3278
     & (2)-cw78%e(3)*p3278(3)
      cauxa=-cw78%ek0*quqd+p32k0*cvqd+p3278k0*cvqu
      cauxb=-cw78%ek0*p3278(2)+p3278k0*cw78%e(2)
      cauxc=+cw78%ek0*p32(2)-p32k0*cw78%e(2)
      u32_78%a(2)=ccl*(cauxa-ceps_0)
      u32_78%b(1)=ccl*(cauxb-ceps_2)
      u32_78%c(2)=ccl*(-cauxc+ceps_1)
      u32_78%d(1)=ccl*cw78%ek0
  
      do i1=1,2
* TLT0_W -- aa=l3_278(i1)%a,cc=l3_278(i1)%c,a1=l3_2(i1)%a,c1=l3_2(i1)%c,
* a2=u32_78%a,b2=u32_78%b,c2=u32_78%c,d2=u32_78%d,prq=s32,nsum=0
      l3_278(i1)%c(2)=l3_2(i1)%c(2)*s32*u32_78%d(1)+l3_2(i1)%a(2
     & )*u32_78%c(2)
      l3_278(i1)%a(2)=l3_2(i1)%c(2)*s32*u32_78%b(1)+l3_2(i1)%a(2
     & )*u32_78%a(2)
      end do
  
* quqd -- p=p378,q=p3278
      quqd=p378(0)*p3278(0)-p378(1)*p3278(1)-p378(2)*p3278(2)-p3
     & 78(3)*p3278(3)
      ccl=1.d0/(f3278)
      do i1=1,2
* TW0 -- qu=p378,qd=p3278,v=ce2(i1)%e,a=u378_2(i1)%a,b=u378_2(i1)%b,c=u3
* 78_2(i1)%c,d=u378_2(i1)%d,cl=ccl,nsum=0
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
      u378_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u378_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u378_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u378_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      do i1=1,2
* TLT0_W -- aa=l3_278(i1)%a,cc=l3_278(i1)%c,a1=l3_78%a,c1=l3_78%c,a2=u37
* 8_2(i1)%a,b2=u378_2(i1)%b,c2=u378_2(i1)%c,d2=u378_2(i1)%d,prq=s378,nsu
* m=1
      l3_278(i1)%c(2)=l3_278(i1)%c(2)+l3_78%c(2)*s378*u378_2(i1)
     & %d(1)+l3_78%a(2)*u378_2(i1)%c(2)
      l3_278(i1)%a(2)=l3_278(i1)%a(2)+l3_78%c(2)*s378*u378_2(i1)
     & %b(1)+l3_78%a(2)*u378_2(i1)%a(2)
      end do
  
        else
  
* II                                                                    
      if (ilept(id5).ne.1) then
* quqd -- p=p32,q=p478
      quqd=p32(0)*p478(0)-p32(1)*p478(1)-p32(2)*p478(2)-p32(3)*p
     & 478(3)
      ccl=wcl/(f478)
      do i2=1,2
* TW0 -- qu=p32,qd=p478,v=cw516(i2)%e,a=u32_516(i2)%a,b=u32_516(i2)%b,c=
* u32_516(i2)%c,d=u32_516(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw516(i2)%ek0*(p32(2)*p478(3)-p478(2)*p32(3))+p32k
     & 0*(cw516(i2)%e(2)*p478(3)-p478(2)*cw516(i2)%e(3))-p478k0*
     & (cw516(i2)%e(2)*p32(3)-p32(2)*cw516(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw516(i2)%e(3)*p32k0+p32(3)*cw516(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw516(i2)%e(3)*p478k0+p478(3)*cw516(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw516(i2)%e(0)*p32(0)-cw516(i2)%e(1)*p32(1)-cw516(i2)
     & %e(2)*p32(2)-cw516(i2)%e(3)*p32(3)
      cvqd=cw516(i2)%e(0)*p478(0)-cw516(i2)%e(1)*p478(1)-cw516(i
     & 2)%e(2)*p478(2)-cw516(i2)%e(3)*p478(3)
      cauxa=-cw516(i2)%ek0*quqd+p32k0*cvqd+p478k0*cvqu
      cauxb=-cw516(i2)%ek0*p478(2)+p478k0*cw516(i2)%e(2)
      cauxc=+cw516(i2)%ek0*p32(2)-p32k0*cw516(i2)%e(2)
      u32_516(i2)%a(2)=ccl*(cauxa-ceps_0)
      u32_516(i2)%b(1)=ccl*(cauxb-ceps_2)
      u32_516(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u32_516(i2)%d(1)=ccl*cw516(i2)%ek0
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l3_2516(i1,i2)%a,cc=l3_2516(i1,i2)%c,a1=l3_2(i1)%a,c1=l3_
* 2(i1)%c,a2=u32_516(i2)%a,b2=u32_516(i2)%b,c2=u32_516(i2)%c,d2=u32_516(
* i2)%d,prq=s32,nsum=0
      l3_2516(i1,i2)%c(2)=l3_2(i1)%c(2)*s32*u32_516(i2)%d(1)+l3_
     & 2(i1)%a(2)*u32_516(i2)%c(2)
      l3_2516(i1,i2)%a(2)=l3_2(i1)%c(2)*s32*u32_516(i2)%b(1)+l3_
     & 2(i1)%a(2)*u32_516(i2)%a(2)
      end do
      end do
      endif
  
* To use also as u3156_2                                                
* quqd -- p=p3156,q=p478
      quqd=p3156(0)*p478(0)-p3156(1)*p478(1)-p3156(2)*p478(2)-p3
     & 156(3)*p478(3)
      ccl=1.d0/(f478)
      do i1=1,2
* TW0 -- qu=p3156,qd=p478,v=ce2(i1)%e,a=u3516_2(i1)%a,b=u3516_2(i1)%b,c=
* u3516_2(i1)%c,d=u3516_2(i1)%d,cl=ccl,nsum=0
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
      u3516_2(i1)%a(2)=ccl*(cauxa-ceps_0)
      u3516_2(i1)%b(1)=ccl*(cauxb-ceps_2)
      u3516_2(i1)%c(2)=ccl*(-cauxc+ceps_1)
      u3516_2(i1)%d(1)=ccl*ce2(i1)%ek0
      end do
  
      if (ilept(id5).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l3_2516(i1,i2)%a,cc=l3_2516(i1,i2)%c,a1=l3_516(i2)%a,c1=l
* 3_516(i2)%c,a2=u3516_2(i1)%a,b2=u3516_2(i1)%b,c2=u3516_2(i1)%c,d2=u351
* 6_2(i1)%d,prq=s3156,nsum=1
      l3_2516(i1,i2)%c(2)=l3_2516(i1,i2)%c(2)+l3_516(i2)%c(2)*s3
     & 156*u3516_2(i1)%d(1)+l3_516(i2)%a(2)*u3516_2(i1)%c(2)
      l3_2516(i1,i2)%a(2)=l3_2516(i1,i2)%a(2)+l3_516(i2)%c(2)*s3
     & 156*u3516_2(i1)%b(1)+l3_516(i2)%a(2)*u3516_2(i1)%a(2)
      end do
      end do
      endif
        endif !iup(id3)
  
* III                                                                   
       if ((iup(id3).eq.1)) then
         if (ilept(id5).ne.1) then
* quqd -- p=p378,q=p42
      quqd=p378(0)*p42(0)-p378(1)*p42(1)-p378(2)*p42(2)-p378(3)*
     & p42(3)
      ccl=wcl/(f42)
      do i2=1,2
* TW0 -- qu=p378,qd=p42,v=cw516(i2)%e,a=u378_516(i2)%a,b=u378_516(i2)%b,
* c=u378_516(i2)%c,d=u378_516(i2)%d,cl=ccl,nsum=0
      ceps_0=-cw516(i2)%ek0*(p378(2)*p42(3)-p42(2)*p378(3))+p378
     & k0*(cw516(i2)%e(2)*p42(3)-p42(2)*cw516(i2)%e(3))-p42k0*(c
     & w516(i2)%e(2)*p378(3)-p378(2)*cw516(i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw516(i2)%e(3)*p378k0+p378(3)*cw516(i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw516(i2)%e(3)*p42k0+p42(3)*cw516(i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cw516(i2)%e(0)*p378(0)-cw516(i2)%e(1)*p378(1)-cw516(i
     & 2)%e(2)*p378(2)-cw516(i2)%e(3)*p378(3)
      cvqd=cw516(i2)%e(0)*p42(0)-cw516(i2)%e(1)*p42(1)-cw516(i2)
     & %e(2)*p42(2)-cw516(i2)%e(3)*p42(3)
      cauxa=-cw516(i2)%ek0*quqd+p378k0*cvqd+p42k0*cvqu
      cauxb=-cw516(i2)%ek0*p42(2)+p42k0*cw516(i2)%e(2)
      cauxc=+cw516(i2)%ek0*p378(2)-p378k0*cw516(i2)%e(2)
      u378_516(i2)%a(2)=ccl*(cauxa-ceps_0)
      u378_516(i2)%b(1)=ccl*(cauxb-ceps_2)
      u378_516(i2)%c(2)=ccl*(-cauxc+ceps_1)
      u378_516(i2)%d(1)=ccl*cw516(i2)%ek0
      end do
  
      do i2=1,2
* TLT0_W -- aa=l3_78516(i2)%a,cc=l3_78516(i2)%c,a1=l3_78%a,c1=l3_78%c,a2
* =u378_516(i2)%a,b2=u378_516(i2)%b,c2=u378_516(i2)%c,d2=u378_516(i2)%d,
* prq=s378,nsum=0
      l3_78516(i2)%c(2)=l3_78%c(2)*s378*u378_516(i2)%d(1)+l3_78%
     & a(2)*u378_516(i2)%c(2)
      l3_78516(i2)%a(2)=l3_78%c(2)*s378*u378_516(i2)%b(1)+l3_78%
     & a(2)*u378_516(i2)%a(2)
      end do
          endif
  
       else
* To use also as u3156_78                                               
* quqd -- p=p3156,q=p42
      quqd=p3156(0)*p42(0)-p3156(1)*p42(1)-p3156(2)*p42(2)-p3156
     & (3)*p42(3)
      ccl=wcl/(f42)
* TW0 -- qu=p3156,qd=p42,v=cw78%e,a=u3516_78%a,b=u3516_78%b,c=u3516_78%c
* ,d=u3516_78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p3156(2)*p42(3)-p42(2)*p3156(3))+p3156k0
     & *(cw78%e(2)*p42(3)-p42(2)*cw78%e(3))-p42k0*(cw78%e(2)*p31
     & 56(3)-p3156(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p3156k0+p3156(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p42k0+p42(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p3156(0)-cw78%e(1)*p3156(1)-cw78%e(2)*p3156
     & (2)-cw78%e(3)*p3156(3)
      cvqd=cw78%e(0)*p42(0)-cw78%e(1)*p42(1)-cw78%e(2)*p42(2)-cw
     & 78%e(3)*p42(3)
      cauxa=-cw78%ek0*quqd+p3156k0*cvqd+p42k0*cvqu
      cauxb=-cw78%ek0*p42(2)+p42k0*cw78%e(2)
      cauxc=+cw78%ek0*p3156(2)-p3156k0*cw78%e(2)
      u3516_78%a(2)=ccl*(cauxa-ceps_0)
      u3516_78%b(1)=ccl*(cauxb-ceps_2)
      u3516_78%c(2)=ccl*(-cauxc+ceps_1)
      u3516_78%d(1)=ccl*cw78%ek0
  
          if (ilept(id5).ne.1) then
      do i2=1,2
* TLT0_W -- aa=l3_78516(i2)%a,cc=l3_78516(i2)%c,a1=l3_516(i2)%a,c1=l3_51
* 6(i2)%c,a2=u3516_78%a,b2=u3516_78%b,c2=u3516_78%c,d2=u3516_78%d,prq=s3
* 156,nsum=0
      l3_78516(i2)%c(2)=l3_516(i2)%c(2)*s3156*u3516_78%d(1)+l3_5
     & 16(i2)%a(2)*u3516_78%c(2)
      l3_78516(i2)%a(2)=l3_516(i2)%c(2)*s3156*u3516_78%b(1)+l3_5
     & 16(i2)%a(2)*u3516_78%a(2)
      end do
          endif
  
       endif !iup(id3)
  
      endif
  
  
* W                                                                     
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,5),a1=l5_134(i1,i3)%a,c1=l5_134(i1,i3)%c,a
* 2=r6_728(i2)%a,b2=r6_728(i2)%b,prq=s5134,bef=cres(i1,i2,i3,5)+,aft=
      cres(i1,i2,i3,5)=cres(i1,i2,i3,5)+(l5_134(i1,i3)%c(2)*s513
     & 4*r6_728(i2)%b(1)+l5_134(i1,i3)%a(2)*r6_728(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,5),a1=l5_1728(i1,i2)%a,c1=l5_1728(i1,i2)%c
* ,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=cres(i1,i2,i3,5)+,aft=
      cres(i1,i2,i3,5)=cres(i1,i2,i3,5)+(l5_1728(i1,i2)%c(2)*s63
     & 4*r6_34(i3)%b(1)+l5_1728(i1,i2)%a(2)*r6_34(i3)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,5),a1=l5_34728(i3,i2)%a,c1=l5_34728(i3,i2)
* %c,a2=r6_1(i1)%a,b2=r6_1(i1)%b,prq=s61,bef=cres(i1,i2,i3,5)+,aft=
      cres(i1,i2,i3,5)=cres(i1,i2,i3,5)+(l5_34728(i3,i2)%c(2)*s6
     & 1*r6_1(i1)%b(1)+l5_34728(i3,i2)%a(2)*r6_1(i1)%a(2))
      end do
      end do
      end do
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,5),a1=l7_234(i2,i3)%a,c1=l7_234(i2,i3)%c,a
* 2=r8_516(i1)%a,b2=r8_516(i1)%b,prq=s7234,bef=cres(i1,i2,i3,5)+,aft=
      cres(i1,i2,i3,5)=cres(i1,i2,i3,5)+(l7_234(i2,i3)%c(2)*s723
     & 4*r8_516(i1)%b(1)+l7_234(i2,i3)%a(2)*r8_516(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,5),a1=l7_2516(i2,i1)%a,c1=l7_2516(i2,i1)%c
* ,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=cres(i1,i2,i3,5)+,aft=
      cres(i1,i2,i3,5)=cres(i1,i2,i3,5)+(l7_2516(i2,i1)%c(2)*s83
     & 4*r8_34(i3)%b(1)+l7_2516(i2,i1)%a(2)*r8_34(i3)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,5),a1=l7_34516(i3,i1)%a,c1=l7_34516(i3,i1)
* %c,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s82,bef=cres(i1,i2,i3,5)+,aft=
      cres(i1,i2,i3,5)=cres(i1,i2,i3,5)+(l7_34516(i3,i1)%c(2)*s8
     & 2*r8_2(i2)%b(1)+l7_34516(i3,i1)%a(2)*r8_2(i2)%a(2))
      end do
      end do
      end do
  
      endif
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,6),a1=l5_234(i2,i3)%a,c1=l5_234(i2,i3)%c,a
* 2=r6_718(i1)%a,b2=r6_718(i1)%b,prq=s5234,bef=cres(i1,i2,i3,6)+,aft=
      cres(i1,i2,i3,6)=cres(i1,i2,i3,6)+(l5_234(i2,i3)%c(2)*s523
     & 4*r6_718(i1)%b(1)+l5_234(i2,i3)%a(2)*r6_718(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,6),a1=l5_2718(i2,i1)%a,c1=l5_2718(i2,i1)%c
* ,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=cres(i1,i2,i3,6)+,aft=
      cres(i1,i2,i3,6)=cres(i1,i2,i3,6)+(l5_2718(i2,i1)%c(2)*s63
     & 4*r6_34(i3)%b(1)+l5_2718(i2,i1)%a(2)*r6_34(i3)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,6),a1=l5_34718(i3,i1)%a,c1=l5_34718(i3,i1)
* %c,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s62,bef=cres(i1,i2,i3,6)+,aft=
      cres(i1,i2,i3,6)=cres(i1,i2,i3,6)+(l5_34718(i3,i1)%c(2)*s6
     & 2*r6_2(i2)%b(1)+l5_34718(i3,i1)%a(2)*r6_2(i2)%a(2))
      end do
      end do
      end do
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,6),a1=l7_134(i1,i3)%a,c1=l7_134(i1,i3)%c,a
* 2=r8_526(i2)%a,b2=r8_526(i2)%b,prq=s7134,bef=cres(i1,i2,i3,6)+,aft=
      cres(i1,i2,i3,6)=cres(i1,i2,i3,6)+(l7_134(i1,i3)%c(2)*s713
     & 4*r8_526(i2)%b(1)+l7_134(i1,i3)%a(2)*r8_526(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,6),a1=l7_1526(i1,i2)%a,c1=l7_1526(i1,i2)%c
* ,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=cres(i1,i2,i3,6)+,aft=
      cres(i1,i2,i3,6)=cres(i1,i2,i3,6)+(l7_1526(i1,i2)%c(2)*s83
     & 4*r8_34(i3)%b(1)+l7_1526(i1,i2)%a(2)*r8_34(i3)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,6),a1=l7_34526(i3,i2)%a,c1=l7_34526(i3,i2)
* %c,a2=r8_1(i1)%a,b2=r8_1(i1)%b,prq=s81,bef=cres(i1,i2,i3,6)+,aft=
      cres(i1,i2,i3,6)=cres(i1,i2,i3,6)+(l7_34526(i3,i2)%c(2)*s8
     & 1*r8_1(i1)%b(1)+l7_34526(i3,i2)%a(2)*r8_1(i1)%a(2))
      end do
      end do
      end do
  
      endif
  
*Z                                                                      
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
       if ((iup(id3).eq.1)) then
* I                                                                     
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,1),a1=l3_178(i1)%a,c1=l3_178(i1)%c,a2=r4_52
* 6(i2)%a,b2=r4_526(i2)%b,prq=s3178,bef=cres(i1,i2,2,1)+,aft=
      cres(i1,i2,2,1)=cres(i1,i2,2,1)+(l3_178(i1)%c(2)*s3178*r4_
     & 526(i2)%b(1)+l3_178(i1)%a(2)*r4_526(i2)%a(2))
      end do
      end do
  
       else
* II                                                                    
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,1),a1=l3_1526(i1,i2)%a,c1=l3_1526(i1,i2)%c,
* a2=r4_78%a,b2=r4_78%b,prq=s478,bef=cres(i1,i2,2,1)+,aft=
      cres(i1,i2,2,1)=cres(i1,i2,2,1)+(l3_1526(i1,i2)%c(2)*s478*
     & r4_78%b(1)+l3_1526(i1,i2)%a(2)*r4_78%a(2))
      end do
      end do
       endif
  
* III                                                                   
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,1),a1=l3_78526(i2)%a,c1=l3_78526(i2)%c,a2=r
* 4_1(i1)%a,b2=r4_1(i1)%b,prq=s41,bef=cres(i1,i2,2,1)+,aft=
      cres(i1,i2,2,1)=cres(i1,i2,2,1)+(l3_78526(i2)%c(2)*s41*r4_
     & 1(i1)%b(1)+l3_78526(i2)%a(2)*r4_1(i1)%a(2))
      end do
      end do
  
  
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,1),a1=l5_278(i2)%a,c1=l5_278(i2)%c,a2=r6_3
* 14(i3,i1)%a,b2=r6_314(i3,i1)%b,prq=s5278,bef=cres(i1,i2,i3,1)+,aft=
      cres(i1,i2,i3,1)=cres(i1,i2,i3,1)+(l5_278(i2)%c(2)*s5278*r
     & 6_314(i3,i1)%b(1)+l5_278(i2)%a(2)*r6_314(i3,i1)%a(2))
      end do
      end do
      end do
  
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,1),a1=l5_2314(i2,i3,i1)%a,c1=l5_2314(i2,i3
* ,i1)%c,a2=r6_78%a,b2=r6_78%b,prq=s678,bef=cres(i1,i2,i3,1)+,aft=
      cres(i1,i2,i3,1)=cres(i1,i2,i3,1)+(l5_2314(i2,i3,i1)%c(2)*
     & s678*r6_78%b(1)+l5_2314(i2,i3,i1)%a(2)*r6_78%a(2))
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,1),a1=l5_78314(i3,i1)%a,c1=l5_78314(i3,i1)
* %c,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s62,bef=cres(i1,i2,i3,1)+,aft=
      cres(i1,i2,i3,1)=cres(i1,i2,i3,1)+(l5_78314(i3,i1)%c(2)*s6
     & 2*r6_2(i2)%b(1)+l5_78314(i3,i1)%a(2)*r6_2(i2)%a(2))
      end do
      end do
      end do
  
      endif
  
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
       if ((iup(id3).eq.1)) then
* I                                                                     
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,2),a1=l3_278(i2)%a,c1=l3_278(i2)%c,a2=r4_51
* 6(i1)%a,b2=r4_516(i1)%b,prq=s3278,bef=cres(i1,i2,2,2)+,aft=
      cres(i1,i2,2,2)=cres(i1,i2,2,2)+(l3_278(i2)%c(2)*s3278*r4_
     & 516(i1)%b(1)+l3_278(i2)%a(2)*r4_516(i1)%a(2))
      end do
      end do
  
       else
* II                                                                    
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,2),a1=l3_2516(i2,i1)%a,c1=l3_2516(i2,i1)%c,
* a2=r4_78%a,b2=r4_78%b,prq=s478,bef=cres(i1,i2,2,2)+,aft=
      cres(i1,i2,2,2)=cres(i1,i2,2,2)+(l3_2516(i2,i1)%c(2)*s478*
     & r4_78%b(1)+l3_2516(i2,i1)%a(2)*r4_78%a(2))
      end do
      end do
       endif
  
* III                                                                   
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,2),a1=l3_78516(i1)%a,c1=l3_78516(i1)%c,a2=r
* 4_2(i2)%a,b2=r4_2(i2)%b,prq=s42,bef=cres(i1,i2,2,2)+,aft=
      cres(i1,i2,2,2)=cres(i1,i2,2,2)+(l3_78516(i1)%c(2)*s42*r4_
     & 2(i2)%b(1)+l3_78516(i1)%a(2)*r4_2(i2)%a(2))
      end do
      end do
  
  
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,2),a1=l5_178(i1)%a,c1=l5_178(i1)%c,a2=r6_3
* 24(i3,i2)%a,b2=r6_324(i3,i2)%b,prq=s5178,bef=cres(i1,i2,i3,2)+,aft=
      cres(i1,i2,i3,2)=cres(i1,i2,i3,2)+(l5_178(i1)%c(2)*s5178*r
     & 6_324(i3,i2)%b(1)+l5_178(i1)%a(2)*r6_324(i3,i2)%a(2))
      end do
      end do
      end do
  
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,2),a1=l5_1324(i1,i3,i2)%a,c1=l5_1324(i1,i3
* ,i2)%c,a2=r6_78%a,b2=r6_78%b,prq=s678,bef=cres(i1,i2,i3,2)+,aft=
      cres(i1,i2,i3,2)=cres(i1,i2,i3,2)+(l5_1324(i1,i3,i2)%c(2)*
     & s678*r6_78%b(1)+l5_1324(i1,i3,i2)%a(2)*r6_78%a(2))
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,2),a1=l5_78324(i3,i2)%a,c1=l5_78324(i3,i2)
* %c,a2=r6_1(i1)%a,b2=r6_1(i1)%b,prq=s61,bef=cres(i1,i2,i3,2)+,aft=
      cres(i1,i2,i3,2)=cres(i1,i2,i3,2)+(l5_78324(i3,i2)%c(2)*s6
     & 1*r6_1(i1)%b(1)+l5_78324(i3,i2)%a(2)*r6_1(i1)%a(2))
      end do
      end do
      end do
  
      endif
  
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
       if (.not.(iup(id3).eq.1)) then
* I                                                                     
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,3),a1=l3_156(i1)%a,c1=l3_156(i1)%c,a2=r4_72
* 8(i2)%a,b2=r4_728(i2)%b,prq=s3156,bef=cres(i1,i2,2,3)+,aft=
      cres(i1,i2,2,3)=cres(i1,i2,2,3)+(l3_156(i1)%c(2)*s3156*r4_
     & 728(i2)%b(1)+l3_156(i1)%a(2)*r4_728(i2)%a(2))
      end do
      end do
  
       else
* II                                                                    
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,3),a1=l3_1728(i1,i2)%a,c1=l3_1728(i1,i2)%c,
* a2=r4_56%a,b2=r4_56%b,prq=s456,bef=cres(i1,i2,2,3)+,aft=
      cres(i1,i2,2,3)=cres(i1,i2,2,3)+(l3_1728(i1,i2)%c(2)*s456*
     & r4_56%b(1)+l3_1728(i1,i2)%a(2)*r4_56%a(2))
      end do
      end do
       endif
  
* III                                                                   
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,3),a1=l3_56728(i2)%a,c1=l3_56728(i2)%c,a2=r
* 4_1(i1)%a,b2=r4_1(i1)%b,prq=s41,bef=cres(i1,i2,2,3)+,aft=
      cres(i1,i2,2,3)=cres(i1,i2,2,3)+(l3_56728(i2)%c(2)*s41*r4_
     & 1(i1)%b(1)+l3_56728(i2)%a(2)*r4_1(i1)%a(2))
      end do
      end do
  
  
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,3),a1=l7_256(i2)%a,c1=l7_256(i2)%c,a2=r8_3
* 14(i3,i1)%a,b2=r8_314(i3,i1)%b,prq=s7256,bef=cres(i1,i2,i3,3)+,aft=
      cres(i1,i2,i3,3)=cres(i1,i2,i3,3)+(l7_256(i2)%c(2)*s7256*r
     & 8_314(i3,i1)%b(1)+l7_256(i2)%a(2)*r8_314(i3,i1)%a(2))
      end do
      end do
      end do
  
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,3),a1=l7_2314(i2,i3,i1)%a,c1=l7_2314(i2,i3
* ,i1)%c,a2=r8_56%a,b2=r8_56%b,prq=s856,bef=cres(i1,i2,i3,3)+,aft=
      cres(i1,i2,i3,3)=cres(i1,i2,i3,3)+(l7_2314(i2,i3,i1)%c(2)*
     & s856*r8_56%b(1)+l7_2314(i2,i3,i1)%a(2)*r8_56%a(2))
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,3),a1=l7_56314(i3,i1)%a,c1=l7_56314(i3,i1)
* %c,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s82,bef=cres(i1,i2,i3,3)+,aft=
      cres(i1,i2,i3,3)=cres(i1,i2,i3,3)+(l7_56314(i3,i1)%c(2)*s8
     & 2*r8_2(i2)%b(1)+l7_56314(i3,i1)%a(2)*r8_2(i2)%a(2))
      end do
      end do
      end do
  
      endif
  
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
       if (.not.(iup(id3).eq.1)) then
* I                                                                     
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,4),a1=l3_256(i2)%a,c1=l3_256(i2)%c,a2=r4_71
* 8(i1)%a,b2=r4_718(i1)%b,prq=s3256,bef=cres(i1,i2,2,4)+,aft=
      cres(i1,i2,2,4)=cres(i1,i2,2,4)+(l3_256(i2)%c(2)*s3256*r4_
     & 718(i1)%b(1)+l3_256(i2)%a(2)*r4_718(i1)%a(2))
      end do
      end do
  
       else
* II                                                                    
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,4),a1=l3_2718(i2,i1)%a,c1=l3_2718(i2,i1)%c,
* a2=r4_56%a,b2=r4_56%b,prq=s456,bef=cres(i1,i2,2,4)+,aft=
      cres(i1,i2,2,4)=cres(i1,i2,2,4)+(l3_2718(i2,i1)%c(2)*s456*
     & r4_56%b(1)+l3_2718(i2,i1)%a(2)*r4_56%a(2))
      end do
      end do
       endif
  
* III                                                                   
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,4),a1=l3_56718(i1)%a,c1=l3_56718(i1)%c,a2=r
* 4_2(i2)%a,b2=r4_2(i2)%b,prq=s42,bef=cres(i1,i2,2,4)+,aft=
      cres(i1,i2,2,4)=cres(i1,i2,2,4)+(l3_56718(i1)%c(2)*s42*r4_
     & 2(i2)%b(1)+l3_56718(i1)%a(2)*r4_2(i2)%a(2))
      end do
      end do
  
  
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,4),a1=l7_156(i1)%a,c1=l7_156(i1)%c,a2=r8_3
* 24(i3,i2)%a,b2=r8_324(i3,i2)%b,prq=s7156,bef=cres(i1,i2,i3,4)+,aft=
      cres(i1,i2,i3,4)=cres(i1,i2,i3,4)+(l7_156(i1)%c(2)*s7156*r
     & 8_324(i3,i2)%b(1)+l7_156(i1)%a(2)*r8_324(i3,i2)%a(2))
      end do
      end do
      end do
  
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,4),a1=l7_1324(i1,i3,i2)%a,c1=l7_1324(i1,i3
* ,i2)%c,a2=r8_56%a,b2=r8_56%b,prq=s856,bef=cres(i1,i2,i3,4)+,aft=
      cres(i1,i2,i3,4)=cres(i1,i2,i3,4)+(l7_1324(i1,i3,i2)%c(2)*
     & s856*r8_56%b(1)+l7_1324(i1,i3,i2)%a(2)*r8_56%a(2))
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,4),a1=l7_56324(i3,i2)%a,c1=l7_56324(i3,i2)
* %c,a2=r8_1(i1)%a,b2=r8_1(i1)%b,prq=s81,bef=cres(i1,i2,i3,4)+,aft=
      cres(i1,i2,i3,4)=cres(i1,i2,i3,4)+(l7_56324(i3,i2)%c(2)*s8
     & 1*r8_1(i1)%b(1)+l7_56324(i3,i2)%a(2)*r8_1(i1)%a(2))
      end do
      end do
      end do
  
      endif
  
  
* -> lines with two gluons                                              
  
*W                                                                      
      if (ilept(id5).ne.1) then
* I                                                                     
* quqd -- p=p512,q=p678
      quqd=p512(0)*p678(0)-p512(1)*p678(1)-p512(2)*p678(2)-p512(
     & 3)*p678(3)
      ccl=zcl(id5)/(f678)
      do i3=1,2
* TW0 -- qu=p512,qd=p678,v=cz34(i3)%e,a=u512_34(i3)%a,b=u512_34(i3)%b,c=
* u512_34(i3)%c,d=u512_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p512(2)*p678(3)-p678(2)*p512(3))+p51
     & 2k0*(cz34(i3)%e(2)*p678(3)-p678(2)*cz34(i3)%e(3))-p678k0*
     & (cz34(i3)%e(2)*p512(3)-p512(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p512k0+p512(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p678k0+p678(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p512(0)-cz34(i3)%e(1)*p512(1)-cz34(i3)%
     & e(2)*p512(2)-cz34(i3)%e(3)*p512(3)
      cvqd=cz34(i3)%e(0)*p678(0)-cz34(i3)%e(1)*p678(1)-cz34(i3)%
     & e(2)*p678(2)-cz34(i3)%e(3)*p678(3)
      cauxa=-cz34(i3)%ek0*quqd+p512k0*cvqd+p678k0*cvqu
      cauxb=-cz34(i3)%ek0*p678(2)+p678k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p512(2)-p512k0*cz34(i3)%e(2)
      u512_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u512_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u512_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u512_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id5)/(f678)
      do i3=1,2
* TW0 -- qu=p512,qd=p678,v=cf34(i3)%e,a=u512_34(i3)%a,b=u512_34(i3)%b,c=
* u512_34(i3)%c,d=u512_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p512(2)*p678(3)-p678(2)*p512(3))+p51
     & 2k0*(cf34(i3)%e(2)*p678(3)-p678(2)*cf34(i3)%e(3))-p678k0*
     & (cf34(i3)%e(2)*p512(3)-p512(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p512k0+p512(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p678k0+p678(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p512(0)-cf34(i3)%e(1)*p512(1)-cf34(i3)%
     & e(2)*p512(2)-cf34(i3)%e(3)*p512(3)
      cvqd=cf34(i3)%e(0)*p678(0)-cf34(i3)%e(1)*p678(1)-cf34(i3)%
     & e(2)*p678(2)-cf34(i3)%e(3)*p678(3)
      cauxa=-cf34(i3)%ek0*quqd+p512k0*cvqd+p678k0*cvqu
      cauxb=-cf34(i3)%ek0*p678(2)+p678k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p512(2)-p512k0*cf34(i3)%e(2)
      u512_34(i3)%a(2)=u512_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u512_34(i3)%b(1)=u512_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u512_34(i3)%c(2)=u512_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u512_34(i3)%d(1)=u512_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      do i3=1,2
* TTR0_W -- aa=r6_3478(i3)%a,bb=r6_3478(i3)%b,a1=u512_34(i3)%a,b1=u512_3
* 4(i3)%b,c1=u512_34(i3)%c,d1=u512_34(i3)%d,a2=r6_78%a,b2=r6_78%b,prq=s6
* 78,nsum=0
      r6_3478(i3)%b(1)=u512_34(i3)%d(1)*s678*r6_78%b(1)+u512_34(
     & i3)%b(1)*r6_78%a(2)
      r6_3478(i3)%a(2)=u512_34(i3)%c(2)*s678*r6_78%b(1)+u512_34(
     & i3)%a(2)*r6_78%a(2)
      end do
  
  
* quqd -- p=p512,q=p634
      quqd=p512(0)*p634(0)-p512(1)*p634(1)-p512(2)*p634(2)-p512(
     & 3)*p634(3)
      ccl=wcl/(f634)
* TW0 -- qu=p512,qd=p634,v=cw78%e,a=u512_78%a,b=u512_78%b,c=u512_78%c,d=
* u512_78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p512(2)*p634(3)-p634(2)*p512(3))+p512k0*
     & (cw78%e(2)*p634(3)-p634(2)*cw78%e(3))-p634k0*(cw78%e(2)*p
     & 512(3)-p512(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p512k0+p512(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p634k0+p634(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p512(0)-cw78%e(1)*p512(1)-cw78%e(2)*p512(2)
     & -cw78%e(3)*p512(3)
      cvqd=cw78%e(0)*p634(0)-cw78%e(1)*p634(1)-cw78%e(2)*p634(2)
     & -cw78%e(3)*p634(3)
      cauxa=-cw78%ek0*quqd+p512k0*cvqd+p634k0*cvqu
      cauxb=-cw78%ek0*p634(2)+p634k0*cw78%e(2)
      cauxc=+cw78%ek0*p512(2)-p512k0*cw78%e(2)
      u512_78%a(2)=ccl*(cauxa-ceps_0)
      u512_78%b(1)=ccl*(cauxb-ceps_2)
      u512_78%c(2)=ccl*(-cauxc+ceps_1)
      u512_78%d(1)=ccl*cw78%ek0
  
      do i3=1,2
* TTR0_W -- aa=r6_3478(i3)%a,bb=r6_3478(i3)%b,a1=u512_78%a,b1=u512_78%b,
* c1=u512_78%c,d1=u512_78%d,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,nsum=
* 1
      r6_3478(i3)%b(1)=r6_3478(i3)%b(1)+u512_78%d(1)*s634*r6_34(
     & i3)%b(1)+u512_78%b(1)*r6_34(i3)%a(2)
      r6_3478(i3)%a(2)=r6_3478(i3)%a(2)+u512_78%c(2)*s634*r6_34(
     & i3)%b(1)+u512_78%a(2)*r6_34(i3)%a(2)
      end do
  
* II                                                                    
* quqd -- p=p534,q=p612
      quqd=p534(0)*p612(0)-p534(1)*p612(1)-p534(2)*p612(2)-p534(
     & 3)*p612(3)
      ccl=wcl/(f612)
* TW0 -- qu=p534,qd=p612,v=cw78%e,a=u534_78%a,b=u534_78%b,c=u534_78%c,d=
* u534_78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p534(2)*p612(3)-p612(2)*p534(3))+p534k0*
     & (cw78%e(2)*p612(3)-p612(2)*cw78%e(3))-p612k0*(cw78%e(2)*p
     & 534(3)-p534(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p534k0+p534(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p612k0+p612(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p534(0)-cw78%e(1)*p534(1)-cw78%e(2)*p534(2)
     & -cw78%e(3)*p534(3)
      cvqd=cw78%e(0)*p612(0)-cw78%e(1)*p612(1)-cw78%e(2)*p612(2)
     & -cw78%e(3)*p612(3)
      cauxa=-cw78%ek0*quqd+p534k0*cvqd+p612k0*cvqu
      cauxb=-cw78%ek0*p612(2)+p612k0*cw78%e(2)
      cauxc=+cw78%ek0*p534(2)-p534k0*cw78%e(2)
      u534_78%a(2)=ccl*(cauxa-ceps_0)
      u534_78%b(1)=ccl*(cauxb-ceps_2)
      u534_78%c(2)=ccl*(-cauxc+ceps_1)
      u534_78%d(1)=ccl*cw78%ek0
  
      do i3=1,2
* TLT0_W -- aa=l5_3478(i3)%a,cc=l5_3478(i3)%c,a1=l5_34(i3)%a,c1=l5_34(i3
* )%c,a2=u534_78%a,b2=u534_78%b,c2=u534_78%c,d2=u534_78%d,prq=s534,nsum=
* 0
      l5_3478(i3)%c(2)=l5_34(i3)%c(2)*s534*u534_78%d(1)+l5_34(i3
     & )%a(2)*u534_78%c(2)
      l5_3478(i3)%a(2)=l5_34(i3)%c(2)*s534*u534_78%b(1)+l5_34(i3
     & )%a(2)*u534_78%a(2)
      end do
  
  
* quqd -- p=p578,q=p612
      quqd=p578(0)*p612(0)-p578(1)*p612(1)-p578(2)*p612(2)-p578(
     & 3)*p612(3)
      ccl=zcl(id6)/(f612)
      do i3=1,2
* TW0 -- qu=p578,qd=p612,v=cz34(i3)%e,a=u578_34(i3)%a,b=u578_34(i3)%b,c=
* u578_34(i3)%c,d=u578_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p578(2)*p612(3)-p612(2)*p578(3))+p57
     & 8k0*(cz34(i3)%e(2)*p612(3)-p612(2)*cz34(i3)%e(3))-p612k0*
     & (cz34(i3)%e(2)*p578(3)-p578(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p578k0+p578(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p612k0+p612(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p578(0)-cz34(i3)%e(1)*p578(1)-cz34(i3)%
     & e(2)*p578(2)-cz34(i3)%e(3)*p578(3)
      cvqd=cz34(i3)%e(0)*p612(0)-cz34(i3)%e(1)*p612(1)-cz34(i3)%
     & e(2)*p612(2)-cz34(i3)%e(3)*p612(3)
      cauxa=-cz34(i3)%ek0*quqd+p578k0*cvqd+p612k0*cvqu
      cauxb=-cz34(i3)%ek0*p612(2)+p612k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p578(2)-p578k0*cz34(i3)%e(2)
      u578_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u578_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u578_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u578_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id6)/(f612)
      do i3=1,2
* TW0 -- qu=p578,qd=p612,v=cf34(i3)%e,a=u578_34(i3)%a,b=u578_34(i3)%b,c=
* u578_34(i3)%c,d=u578_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p578(2)*p612(3)-p612(2)*p578(3))+p57
     & 8k0*(cf34(i3)%e(2)*p612(3)-p612(2)*cf34(i3)%e(3))-p612k0*
     & (cf34(i3)%e(2)*p578(3)-p578(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p578k0+p578(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p612k0+p612(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p578(0)-cf34(i3)%e(1)*p578(1)-cf34(i3)%
     & e(2)*p578(2)-cf34(i3)%e(3)*p578(3)
      cvqd=cf34(i3)%e(0)*p612(0)-cf34(i3)%e(1)*p612(1)-cf34(i3)%
     & e(2)*p612(2)-cf34(i3)%e(3)*p612(3)
      cauxa=-cf34(i3)%ek0*quqd+p578k0*cvqd+p612k0*cvqu
      cauxb=-cf34(i3)%ek0*p612(2)+p612k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p578(2)-p578k0*cf34(i3)%e(2)
      u578_34(i3)%a(2)=u578_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u578_34(i3)%b(1)=u578_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u578_34(i3)%c(2)=u578_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u578_34(i3)%d(1)=u578_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      do i3=1,2
* TLT0_W -- aa=l5_3478(i3)%a,cc=l5_3478(i3)%c,a1=l5_78%a,c1=l5_78%c,a2=u
* 578_34(i3)%a,b2=u578_34(i3)%b,c2=u578_34(i3)%c,d2=u578_34(i3)%d,prq=s5
* 78,nsum=1
      l5_3478(i3)%c(2)=l5_3478(i3)%c(2)+l5_78%c(2)*s578*u578_34(
     & i3)%d(1)+l5_78%a(2)*u578_34(i3)%c(2)
      l5_3478(i3)%a(2)=l5_3478(i3)%a(2)+l5_78%c(2)*s578*u578_34(
     & i3)%b(1)+l5_78%a(2)*u578_34(i3)%a(2)
      end do
  
* III                                                                   
      do i2=1,2
* TTR0_W -- aa=r6_278(i2)%a,bb=r6_278(i2)%b,a1=u5314_2(i2)%a,b1=u5314_2(
* i2)%b,c1=u5314_2(i2)%c,d1=u5314_2(i2)%d,a2=r6_78%a,b2=r6_78%b,prq=s678
* ,nsum=0
      r6_278(i2)%b(1)=u5314_2(i2)%d(1)*s678*r6_78%b(1)+u5314_2(i
     & 2)%b(1)*r6_78%a(2)
      r6_278(i2)%a(2)=u5314_2(i2)%c(2)*s678*r6_78%b(1)+u5314_2(i
     & 2)%a(2)*r6_78%a(2)
      end do
  
      do i2=1,2
* TTR0_W -- aa=r6_278(i2)%a,bb=r6_278(i2)%b,a1=u5314_78%a,b1=u5314_78%b,
* c1=u5314_78%c,d1=u5314_78%d,a2=r6_2(i2)%a,b2=r6_2(i2)%b,prq=s62,nsum=1
      r6_278(i2)%b(1)=r6_278(i2)%b(1)+u5314_78%d(1)*s62*r6_2(i2)
     & %b(1)+u5314_78%b(1)*r6_2(i2)%a(2)
      r6_278(i2)%a(2)=r6_278(i2)%a(2)+u5314_78%c(2)*s62*r6_2(i2)
     & %b(1)+u5314_78%a(2)*r6_2(i2)%a(2)
      end do
  
      do i2=1,2
      do i3=1,2
* TTR0_W -- aa=r6_234(i2,i3)%a,bb=r6_234(i2,i3)%b,a1=u5718_2(i2)%a,b1=u5
* 718_2(i2)%b,c1=u5718_2(i2)%c,d1=u5718_2(i2)%d,a2=r6_34(i3)%a,b2=r6_34(
* i3)%b,prq=s634,nsum=0
      r6_234(i2,i3)%b(1)=u5718_2(i2)%d(1)*s634*r6_34(i3)%b(1)+u5
     & 718_2(i2)%b(1)*r6_34(i3)%a(2)
      r6_234(i2,i3)%a(2)=u5718_2(i2)%c(2)*s634*r6_34(i3)%b(1)+u5
     & 718_2(i2)%a(2)*r6_34(i3)%a(2)
      end do
      end do
  
      do i2=1,2
      do i3=1,2
* TTR0_W -- aa=r6_234(i2,i3)%a,bb=r6_234(i2,i3)%b,a1=u5718_34(i3)%a,b1=u
* 5718_34(i3)%b,c1=u5718_34(i3)%c,d1=u5718_34(i3)%d,a2=r6_2(i2)%a,b2=r6_
* 2(i2)%b,prq=s62,nsum=1
      r6_234(i2,i3)%b(1)=r6_234(i2,i3)%b(1)+u5718_34(i3)%d(1)*s6
     & 2*r6_2(i2)%b(1)+u5718_34(i3)%b(1)*r6_2(i2)%a(2)
      r6_234(i2,i3)%a(2)=r6_234(i2,i3)%a(2)+u5718_34(i3)%c(2)*s6
     & 2*r6_2(i2)%b(1)+u5718_34(i3)%a(2)*r6_2(i2)%a(2)
      end do
      end do
      do i2=1,2
* TTR0_W -- aa=r6_178(i2)%a,bb=r6_178(i2)%b,a1=u5324_1(i2)%a,b1=u5324_1(
* i2)%b,c1=u5324_1(i2)%c,d1=u5324_1(i2)%d,a2=r6_78%a,b2=r6_78%b,prq=s678
* ,nsum=0
      r6_178(i2)%b(1)=u5324_1(i2)%d(1)*s678*r6_78%b(1)+u5324_1(i
     & 2)%b(1)*r6_78%a(2)
      r6_178(i2)%a(2)=u5324_1(i2)%c(2)*s678*r6_78%b(1)+u5324_1(i
     & 2)%a(2)*r6_78%a(2)
      end do
  
      do i2=1,2
* TTR0_W -- aa=r6_178(i2)%a,bb=r6_178(i2)%b,a1=u5324_78%a,b1=u5324_78%b,
* c1=u5324_78%c,d1=u5324_78%d,a2=r6_1(i2)%a,b2=r6_1(i2)%b,prq=s61,nsum=1
      r6_178(i2)%b(1)=r6_178(i2)%b(1)+u5324_78%d(1)*s61*r6_1(i2)
     & %b(1)+u5324_78%b(1)*r6_1(i2)%a(2)
      r6_178(i2)%a(2)=r6_178(i2)%a(2)+u5324_78%c(2)*s61*r6_1(i2)
     & %b(1)+u5324_78%a(2)*r6_1(i2)%a(2)
      end do
  
      do i2=1,2
      do i3=1,2
* TTR0_W -- aa=r6_134(i2,i3)%a,bb=r6_134(i2,i3)%b,a1=u5728_1(i2)%a,b1=u5
* 728_1(i2)%b,c1=u5728_1(i2)%c,d1=u5728_1(i2)%d,a2=r6_34(i3)%a,b2=r6_34(
* i3)%b,prq=s634,nsum=0
      r6_134(i2,i3)%b(1)=u5728_1(i2)%d(1)*s634*r6_34(i3)%b(1)+u5
     & 728_1(i2)%b(1)*r6_34(i3)%a(2)
      r6_134(i2,i3)%a(2)=u5728_1(i2)%c(2)*s634*r6_34(i3)%b(1)+u5
     & 728_1(i2)%a(2)*r6_34(i3)%a(2)
      end do
      end do
  
      do i2=1,2
      do i3=1,2
* TTR0_W -- aa=r6_134(i2,i3)%a,bb=r6_134(i2,i3)%b,a1=u5728_34(i3)%a,b1=u
* 5728_34(i3)%b,c1=u5728_34(i3)%c,d1=u5728_34(i3)%d,a2=r6_1(i2)%a,b2=r6_
* 1(i2)%b,prq=s61,nsum=1
      r6_134(i2,i3)%b(1)=r6_134(i2,i3)%b(1)+u5728_34(i3)%d(1)*s6
     & 1*r6_1(i2)%b(1)+u5728_34(i3)%b(1)*r6_1(i2)%a(2)
      r6_134(i2,i3)%a(2)=r6_134(i2,i3)%a(2)+u5728_34(i3)%c(2)*s6
     & 1*r6_1(i2)%b(1)+u5728_34(i3)%a(2)*r6_1(i2)%a(2)
      end do
      end do
  
* IV                                                                    
* quqd -- p=p534,q=p678
      quqd=p534(0)*p678(0)-p534(1)*p678(1)-p534(2)*p678(2)-p534(
     & 3)*p678(3)
      ccl=1d0/(f678)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p534,qd=p678,v=ctrip12(i1,i2)%e,a=u534_gg(i1,i2)%a,b=u534_gg
* (i1,i2)%b,c=u534_gg(i1,i2)%c,d=u534_gg(i1,i2)%d,cl=ccl,nsum=0
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
      u534_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u534_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u534_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u534_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      end do
      end do
  
* quqd -- p=p578,q=p634
      quqd=p578(0)*p634(0)-p578(1)*p634(1)-p578(2)*p634(2)-p578(
     & 3)*p634(3)
      ccl=1d0/(f634)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p578,qd=p634,v=ctrip12(i1,i2)%e,a=u578_gg(i1,i2)%a,b=u578_gg
* (i1,i2)%b,c=u578_gg(i1,i2)%c,d=u578_gg(i1,i2)%d,cl=ccl,nsum=0
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
      u578_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u578_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u578_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u578_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      end do
      end do
  
      endif
*W                                                                      
      if (ilept(id7).ne.1) then
* I                                                                     
* quqd -- p=p712,q=p856
      quqd=p712(0)*p856(0)-p712(1)*p856(1)-p712(2)*p856(2)-p712(
     & 3)*p856(3)
      ccl=zcl(id7)/(f856)
      do i3=1,2
* TW0 -- qu=p712,qd=p856,v=cz34(i3)%e,a=u712_34(i3)%a,b=u712_34(i3)%b,c=
* u712_34(i3)%c,d=u712_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p712(2)*p856(3)-p856(2)*p712(3))+p71
     & 2k0*(cz34(i3)%e(2)*p856(3)-p856(2)*cz34(i3)%e(3))-p856k0*
     & (cz34(i3)%e(2)*p712(3)-p712(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p712k0+p712(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p856k0+p856(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p712(0)-cz34(i3)%e(1)*p712(1)-cz34(i3)%
     & e(2)*p712(2)-cz34(i3)%e(3)*p712(3)
      cvqd=cz34(i3)%e(0)*p856(0)-cz34(i3)%e(1)*p856(1)-cz34(i3)%
     & e(2)*p856(2)-cz34(i3)%e(3)*p856(3)
      cauxa=-cz34(i3)%ek0*quqd+p712k0*cvqd+p856k0*cvqu
      cauxb=-cz34(i3)%ek0*p856(2)+p856k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p712(2)-p712k0*cz34(i3)%e(2)
      u712_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u712_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u712_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u712_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id7)/(f856)
      do i3=1,2
* TW0 -- qu=p712,qd=p856,v=cf34(i3)%e,a=u712_34(i3)%a,b=u712_34(i3)%b,c=
* u712_34(i3)%c,d=u712_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p712(2)*p856(3)-p856(2)*p712(3))+p71
     & 2k0*(cf34(i3)%e(2)*p856(3)-p856(2)*cf34(i3)%e(3))-p856k0*
     & (cf34(i3)%e(2)*p712(3)-p712(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p712k0+p712(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p856k0+p856(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p712(0)-cf34(i3)%e(1)*p712(1)-cf34(i3)%
     & e(2)*p712(2)-cf34(i3)%e(3)*p712(3)
      cvqd=cf34(i3)%e(0)*p856(0)-cf34(i3)%e(1)*p856(1)-cf34(i3)%
     & e(2)*p856(2)-cf34(i3)%e(3)*p856(3)
      cauxa=-cf34(i3)%ek0*quqd+p712k0*cvqd+p856k0*cvqu
      cauxb=-cf34(i3)%ek0*p856(2)+p856k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p712(2)-p712k0*cf34(i3)%e(2)
      u712_34(i3)%a(2)=u712_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u712_34(i3)%b(1)=u712_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u712_34(i3)%c(2)=u712_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u712_34(i3)%d(1)=u712_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      do i3=1,2
* TTR0_W -- aa=r8_3456(i3)%a,bb=r8_3456(i3)%b,a1=u712_34(i3)%a,b1=u712_3
* 4(i3)%b,c1=u712_34(i3)%c,d1=u712_34(i3)%d,a2=r8_56%a,b2=r8_56%b,prq=s8
* 56,nsum=0
      r8_3456(i3)%b(1)=u712_34(i3)%d(1)*s856*r8_56%b(1)+u712_34(
     & i3)%b(1)*r8_56%a(2)
      r8_3456(i3)%a(2)=u712_34(i3)%c(2)*s856*r8_56%b(1)+u712_34(
     & i3)%a(2)*r8_56%a(2)
      end do
  
  
* quqd -- p=p712,q=p834
      quqd=p712(0)*p834(0)-p712(1)*p834(1)-p712(2)*p834(2)-p712(
     & 3)*p834(3)
      ccl=wcl/(f834)
* TW0 -- qu=p712,qd=p834,v=cw56%e,a=u712_56%a,b=u712_56%b,c=u712_56%c,d=
* u712_56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p712(2)*p834(3)-p834(2)*p712(3))+p712k0*
     & (cw56%e(2)*p834(3)-p834(2)*cw56%e(3))-p834k0*(cw56%e(2)*p
     & 712(3)-p712(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p712k0+p712(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p834k0+p834(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p712(0)-cw56%e(1)*p712(1)-cw56%e(2)*p712(2)
     & -cw56%e(3)*p712(3)
      cvqd=cw56%e(0)*p834(0)-cw56%e(1)*p834(1)-cw56%e(2)*p834(2)
     & -cw56%e(3)*p834(3)
      cauxa=-cw56%ek0*quqd+p712k0*cvqd+p834k0*cvqu
      cauxb=-cw56%ek0*p834(2)+p834k0*cw56%e(2)
      cauxc=+cw56%ek0*p712(2)-p712k0*cw56%e(2)
      u712_56%a(2)=ccl*(cauxa-ceps_0)
      u712_56%b(1)=ccl*(cauxb-ceps_2)
      u712_56%c(2)=ccl*(-cauxc+ceps_1)
      u712_56%d(1)=ccl*cw56%ek0
  
      do i3=1,2
* TTR0_W -- aa=r8_3456(i3)%a,bb=r8_3456(i3)%b,a1=u712_56%a,b1=u712_56%b,
* c1=u712_56%c,d1=u712_56%d,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,nsum=
* 1
      r8_3456(i3)%b(1)=r8_3456(i3)%b(1)+u712_56%d(1)*s834*r8_34(
     & i3)%b(1)+u712_56%b(1)*r8_34(i3)%a(2)
      r8_3456(i3)%a(2)=r8_3456(i3)%a(2)+u712_56%c(2)*s834*r8_34(
     & i3)%b(1)+u712_56%a(2)*r8_34(i3)%a(2)
      end do
  
* II                                                                    
* quqd -- p=p734,q=p812
      quqd=p734(0)*p812(0)-p734(1)*p812(1)-p734(2)*p812(2)-p734(
     & 3)*p812(3)
      ccl=wcl/(f812)
* TW0 -- qu=p734,qd=p812,v=cw56%e,a=u734_56%a,b=u734_56%b,c=u734_56%c,d=
* u734_56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p734(2)*p812(3)-p812(2)*p734(3))+p734k0*
     & (cw56%e(2)*p812(3)-p812(2)*cw56%e(3))-p812k0*(cw56%e(2)*p
     & 734(3)-p734(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p734k0+p734(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p812k0+p812(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p734(0)-cw56%e(1)*p734(1)-cw56%e(2)*p734(2)
     & -cw56%e(3)*p734(3)
      cvqd=cw56%e(0)*p812(0)-cw56%e(1)*p812(1)-cw56%e(2)*p812(2)
     & -cw56%e(3)*p812(3)
      cauxa=-cw56%ek0*quqd+p734k0*cvqd+p812k0*cvqu
      cauxb=-cw56%ek0*p812(2)+p812k0*cw56%e(2)
      cauxc=+cw56%ek0*p734(2)-p734k0*cw56%e(2)
      u734_56%a(2)=ccl*(cauxa-ceps_0)
      u734_56%b(1)=ccl*(cauxb-ceps_2)
      u734_56%c(2)=ccl*(-cauxc+ceps_1)
      u734_56%d(1)=ccl*cw56%ek0
  
      do i3=1,2
* TLT0_W -- aa=l7_3456(i3)%a,cc=l7_3456(i3)%c,a1=l7_34(i3)%a,c1=l7_34(i3
* )%c,a2=u734_56%a,b2=u734_56%b,c2=u734_56%c,d2=u734_56%d,prq=s734,nsum=
* 0
      l7_3456(i3)%c(2)=l7_34(i3)%c(2)*s734*u734_56%d(1)+l7_34(i3
     & )%a(2)*u734_56%c(2)
      l7_3456(i3)%a(2)=l7_34(i3)%c(2)*s734*u734_56%b(1)+l7_34(i3
     & )%a(2)*u734_56%a(2)
      end do
  
  
* quqd -- p=p756,q=p812
      quqd=p756(0)*p812(0)-p756(1)*p812(1)-p756(2)*p812(2)-p756(
     & 3)*p812(3)
      ccl=zcl(id8)/(f812)
      do i3=1,2
* TW0 -- qu=p756,qd=p812,v=cz34(i3)%e,a=u756_34(i3)%a,b=u756_34(i3)%b,c=
* u756_34(i3)%c,d=u756_34(i3)%d,cl=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p756(2)*p812(3)-p812(2)*p756(3))+p75
     & 6k0*(cz34(i3)%e(2)*p812(3)-p812(2)*cz34(i3)%e(3))-p812k0*
     & (cz34(i3)%e(2)*p756(3)-p756(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p756k0+p756(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p812k0+p812(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p756(0)-cz34(i3)%e(1)*p756(1)-cz34(i3)%
     & e(2)*p756(2)-cz34(i3)%e(3)*p756(3)
      cvqd=cz34(i3)%e(0)*p812(0)-cz34(i3)%e(1)*p812(1)-cz34(i3)%
     & e(2)*p812(2)-cz34(i3)%e(3)*p812(3)
      cauxa=-cz34(i3)%ek0*quqd+p756k0*cvqd+p812k0*cvqu
      cauxb=-cz34(i3)%ek0*p812(2)+p812k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p756(2)-p756k0*cz34(i3)%e(2)
      u756_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      u756_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      u756_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      u756_34(i3)%d(1)=ccl*cz34(i3)%ek0
      end do
  
       if (ineutri(id3).ne.1) then
      ccl=fcl(id8)/(f812)
      do i3=1,2
* TW0 -- qu=p756,qd=p812,v=cf34(i3)%e,a=u756_34(i3)%a,b=u756_34(i3)%b,c=
* u756_34(i3)%c,d=u756_34(i3)%d,cl=ccl,nsum=1
      ceps_0=-cf34(i3)%ek0*(p756(2)*p812(3)-p812(2)*p756(3))+p75
     & 6k0*(cf34(i3)%e(2)*p812(3)-p812(2)*cf34(i3)%e(3))-p812k0*
     & (cf34(i3)%e(2)*p756(3)-p756(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p756k0+p756(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p812k0+p812(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p756(0)-cf34(i3)%e(1)*p756(1)-cf34(i3)%
     & e(2)*p756(2)-cf34(i3)%e(3)*p756(3)
      cvqd=cf34(i3)%e(0)*p812(0)-cf34(i3)%e(1)*p812(1)-cf34(i3)%
     & e(2)*p812(2)-cf34(i3)%e(3)*p812(3)
      cauxa=-cf34(i3)%ek0*quqd+p756k0*cvqd+p812k0*cvqu
      cauxb=-cf34(i3)%ek0*p812(2)+p812k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p756(2)-p756k0*cf34(i3)%e(2)
      u756_34(i3)%a(2)=u756_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      u756_34(i3)%b(1)=u756_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      u756_34(i3)%c(2)=u756_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      u756_34(i3)%d(1)=u756_34(i3)%d(1)+ccl*cf34(i3)%ek0
      end do
       endif
  
      do i3=1,2
* TLT0_W -- aa=l7_3456(i3)%a,cc=l7_3456(i3)%c,a1=l7_56%a,c1=l7_56%c,a2=u
* 756_34(i3)%a,b2=u756_34(i3)%b,c2=u756_34(i3)%c,d2=u756_34(i3)%d,prq=s7
* 56,nsum=1
      l7_3456(i3)%c(2)=l7_3456(i3)%c(2)+l7_56%c(2)*s756*u756_34(
     & i3)%d(1)+l7_56%a(2)*u756_34(i3)%c(2)
      l7_3456(i3)%a(2)=l7_3456(i3)%a(2)+l7_56%c(2)*s756*u756_34(
     & i3)%b(1)+l7_56%a(2)*u756_34(i3)%a(2)
      end do
  
* III                                                                   
      do i2=1,2
* TTR0_W -- aa=r8_256(i2)%a,bb=r8_256(i2)%b,a1=u7314_2(i2)%a,b1=u7314_2(
* i2)%b,c1=u7314_2(i2)%c,d1=u7314_2(i2)%d,a2=r8_56%a,b2=r8_56%b,prq=s856
* ,nsum=0
      r8_256(i2)%b(1)=u7314_2(i2)%d(1)*s856*r8_56%b(1)+u7314_2(i
     & 2)%b(1)*r8_56%a(2)
      r8_256(i2)%a(2)=u7314_2(i2)%c(2)*s856*r8_56%b(1)+u7314_2(i
     & 2)%a(2)*r8_56%a(2)
      end do
  
      do i2=1,2
* TTR0_W -- aa=r8_256(i2)%a,bb=r8_256(i2)%b,a1=u7314_56%a,b1=u7314_56%b,
* c1=u7314_56%c,d1=u7314_56%d,a2=r8_2(i2)%a,b2=r8_2(i2)%b,prq=s82,nsum=1
      r8_256(i2)%b(1)=r8_256(i2)%b(1)+u7314_56%d(1)*s82*r8_2(i2)
     & %b(1)+u7314_56%b(1)*r8_2(i2)%a(2)
      r8_256(i2)%a(2)=r8_256(i2)%a(2)+u7314_56%c(2)*s82*r8_2(i2)
     & %b(1)+u7314_56%a(2)*r8_2(i2)%a(2)
      end do
  
      do i2=1,2
      do i3=1,2
* TTR0_W -- aa=r8_234(i2,i3)%a,bb=r8_234(i2,i3)%b,a1=u7516_2(i2)%a,b1=u7
* 516_2(i2)%b,c1=u7516_2(i2)%c,d1=u7516_2(i2)%d,a2=r8_34(i3)%a,b2=r8_34(
* i3)%b,prq=s834,nsum=0
      r8_234(i2,i3)%b(1)=u7516_2(i2)%d(1)*s834*r8_34(i3)%b(1)+u7
     & 516_2(i2)%b(1)*r8_34(i3)%a(2)
      r8_234(i2,i3)%a(2)=u7516_2(i2)%c(2)*s834*r8_34(i3)%b(1)+u7
     & 516_2(i2)%a(2)*r8_34(i3)%a(2)
      end do
      end do
  
      do i2=1,2
      do i3=1,2
* TTR0_W -- aa=r8_234(i2,i3)%a,bb=r8_234(i2,i3)%b,a1=u7516_34(i3)%a,b1=u
* 7516_34(i3)%b,c1=u7516_34(i3)%c,d1=u7516_34(i3)%d,a2=r8_2(i2)%a,b2=r8_
* 2(i2)%b,prq=s82,nsum=1
      r8_234(i2,i3)%b(1)=r8_234(i2,i3)%b(1)+u7516_34(i3)%d(1)*s8
     & 2*r8_2(i2)%b(1)+u7516_34(i3)%b(1)*r8_2(i2)%a(2)
      r8_234(i2,i3)%a(2)=r8_234(i2,i3)%a(2)+u7516_34(i3)%c(2)*s8
     & 2*r8_2(i2)%b(1)+u7516_34(i3)%a(2)*r8_2(i2)%a(2)
      end do
      end do
      do i2=1,2
* TTR0_W -- aa=r8_156(i2)%a,bb=r8_156(i2)%b,a1=u7324_1(i2)%a,b1=u7324_1(
* i2)%b,c1=u7324_1(i2)%c,d1=u7324_1(i2)%d,a2=r8_56%a,b2=r8_56%b,prq=s856
* ,nsum=0
      r8_156(i2)%b(1)=u7324_1(i2)%d(1)*s856*r8_56%b(1)+u7324_1(i
     & 2)%b(1)*r8_56%a(2)
      r8_156(i2)%a(2)=u7324_1(i2)%c(2)*s856*r8_56%b(1)+u7324_1(i
     & 2)%a(2)*r8_56%a(2)
      end do
  
      do i2=1,2
* TTR0_W -- aa=r8_156(i2)%a,bb=r8_156(i2)%b,a1=u7324_56%a,b1=u7324_56%b,
* c1=u7324_56%c,d1=u7324_56%d,a2=r8_1(i2)%a,b2=r8_1(i2)%b,prq=s81,nsum=1
      r8_156(i2)%b(1)=r8_156(i2)%b(1)+u7324_56%d(1)*s81*r8_1(i2)
     & %b(1)+u7324_56%b(1)*r8_1(i2)%a(2)
      r8_156(i2)%a(2)=r8_156(i2)%a(2)+u7324_56%c(2)*s81*r8_1(i2)
     & %b(1)+u7324_56%a(2)*r8_1(i2)%a(2)
      end do
  
      do i2=1,2
      do i3=1,2
* TTR0_W -- aa=r8_134(i2,i3)%a,bb=r8_134(i2,i3)%b,a1=u7526_1(i2)%a,b1=u7
* 526_1(i2)%b,c1=u7526_1(i2)%c,d1=u7526_1(i2)%d,a2=r8_34(i3)%a,b2=r8_34(
* i3)%b,prq=s834,nsum=0
      r8_134(i2,i3)%b(1)=u7526_1(i2)%d(1)*s834*r8_34(i3)%b(1)+u7
     & 526_1(i2)%b(1)*r8_34(i3)%a(2)
      r8_134(i2,i3)%a(2)=u7526_1(i2)%c(2)*s834*r8_34(i3)%b(1)+u7
     & 526_1(i2)%a(2)*r8_34(i3)%a(2)
      end do
      end do
  
      do i2=1,2
      do i3=1,2
* TTR0_W -- aa=r8_134(i2,i3)%a,bb=r8_134(i2,i3)%b,a1=u7526_34(i3)%a,b1=u
* 7526_34(i3)%b,c1=u7526_34(i3)%c,d1=u7526_34(i3)%d,a2=r8_1(i2)%a,b2=r8_
* 1(i2)%b,prq=s81,nsum=1
      r8_134(i2,i3)%b(1)=r8_134(i2,i3)%b(1)+u7526_34(i3)%d(1)*s8
     & 1*r8_1(i2)%b(1)+u7526_34(i3)%b(1)*r8_1(i2)%a(2)
      r8_134(i2,i3)%a(2)=r8_134(i2,i3)%a(2)+u7526_34(i3)%c(2)*s8
     & 1*r8_1(i2)%b(1)+u7526_34(i3)%a(2)*r8_1(i2)%a(2)
      end do
      end do
  
* IV                                                                    
* quqd -- p=p734,q=p856
      quqd=p734(0)*p856(0)-p734(1)*p856(1)-p734(2)*p856(2)-p734(
     & 3)*p856(3)
      ccl=1d0/(f856)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p734,qd=p856,v=ctrip12(i1,i2)%e,a=u734_gg(i1,i2)%a,b=u734_gg
* (i1,i2)%b,c=u734_gg(i1,i2)%c,d=u734_gg(i1,i2)%d,cl=ccl,nsum=0
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
      u734_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u734_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u734_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u734_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      end do
      end do
  
* quqd -- p=p756,q=p834
      quqd=p756(0)*p834(0)-p756(1)*p834(1)-p756(2)*p834(2)-p756(
     & 3)*p834(3)
      ccl=1d0/(f834)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p756,qd=p834,v=ctrip12(i1,i2)%e,a=u756_gg(i1,i2)%a,b=u756_gg
* (i1,i2)%b,c=u756_gg(i1,i2)%c,d=u756_gg(i1,i2)%d,cl=ccl,nsum=0
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
      u756_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u756_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u756_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u756_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      end do
      end do
  
      endif
  
  
*Z                                                                      
      if (ilept(id3).ne.1) then
* I                                                                     
       if (iup(id3).ne.1) then
* quqd -- p=p312,q=p478
      quqd=p312(0)*p478(0)-p312(1)*p478(1)-p312(2)*p478(2)-p312(
     & 3)*p478(3)
      ccl=wcl/(f478)
* TW0 -- qu=p312,qd=p478,v=cw56%e,a=u312_56%a,b=u312_56%b,c=u312_56%c,d=
* u312_56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p312(2)*p478(3)-p478(2)*p312(3))+p312k0*
     & (cw56%e(2)*p478(3)-p478(2)*cw56%e(3))-p478k0*(cw56%e(2)*p
     & 312(3)-p312(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p312k0+p312(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p478k0+p478(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p312(0)-cw56%e(1)*p312(1)-cw56%e(2)*p312(2)
     & -cw56%e(3)*p312(3)
      cvqd=cw56%e(0)*p478(0)-cw56%e(1)*p478(1)-cw56%e(2)*p478(2)
     & -cw56%e(3)*p478(3)
      cauxa=-cw56%ek0*quqd+p312k0*cvqd+p478k0*cvqu
      cauxb=-cw56%ek0*p478(2)+p478k0*cw56%e(2)
      cauxc=+cw56%ek0*p312(2)-p312k0*cw56%e(2)
      u312_56%a(2)=ccl*(cauxa-ceps_0)
      u312_56%b(1)=ccl*(cauxb-ceps_2)
      u312_56%c(2)=ccl*(-cauxc+ceps_1)
      u312_56%d(1)=ccl*cw56%ek0
  
* TTR0_W -- aa=r4_5678%a,bb=r4_5678%b,a1=u312_56%a,b1=u312_56%b,c1=u312_
* 56%c,d1=u312_56%d,a2=r4_78%a,b2=r4_78%b,prq=s478,nsum=0
      r4_5678%b(1)=u312_56%d(1)*s478*r4_78%b(1)+u312_56%b(1)*r4_
     & 78%a(2)
      r4_5678%a(2)=u312_56%c(2)*s478*r4_78%b(1)+u312_56%a(2)*r4_
     & 78%a(2)
       else
* quqd -- p=p312,q=p456
      quqd=p312(0)*p456(0)-p312(1)*p456(1)-p312(2)*p456(2)-p312(
     & 3)*p456(3)
      ccl=wcl/(f456)
* TW0 -- qu=p312,qd=p456,v=cw78%e,a=u312_78%a,b=u312_78%b,c=u312_78%c,d=
* u312_78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p312(2)*p456(3)-p456(2)*p312(3))+p312k0*
     & (cw78%e(2)*p456(3)-p456(2)*cw78%e(3))-p456k0*(cw78%e(2)*p
     & 312(3)-p312(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p312k0+p312(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p456k0+p456(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p312(0)-cw78%e(1)*p312(1)-cw78%e(2)*p312(2)
     & -cw78%e(3)*p312(3)
      cvqd=cw78%e(0)*p456(0)-cw78%e(1)*p456(1)-cw78%e(2)*p456(2)
     & -cw78%e(3)*p456(3)
      cauxa=-cw78%ek0*quqd+p312k0*cvqd+p456k0*cvqu
      cauxb=-cw78%ek0*p456(2)+p456k0*cw78%e(2)
      cauxc=+cw78%ek0*p312(2)-p312k0*cw78%e(2)
      u312_78%a(2)=ccl*(cauxa-ceps_0)
      u312_78%b(1)=ccl*(cauxb-ceps_2)
      u312_78%c(2)=ccl*(-cauxc+ceps_1)
      u312_78%d(1)=ccl*cw78%ek0
  
* TTR0_W -- aa=r4_5678%a,bb=r4_5678%b,a1=u312_78%a,b1=u312_78%b,c1=u312_
* 78%c,d1=u312_78%d,a2=r4_56%a,b2=r4_56%b,prq=s456,nsum=0
      r4_5678%b(1)=u312_78%d(1)*s456*r4_56%b(1)+u312_78%b(1)*r4_
     & 56%a(2)
      r4_5678%a(2)=u312_78%c(2)*s456*r4_56%b(1)+u312_78%a(2)*r4_
     & 56%a(2)
       endif
  
* II                                                                    
       if (iup(id3).ne.1) then
* quqd -- p=p356,q=p412
      quqd=p356(0)*p412(0)-p356(1)*p412(1)-p356(2)*p412(2)-p356(
     & 3)*p412(3)
      ccl=wcl/(f412)
* TW0 -- qu=p356,qd=p412,v=cw78%e,a=u356_78%a,b=u356_78%b,c=u356_78%c,d=
* u356_78%d,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p356(2)*p412(3)-p412(2)*p356(3))+p356k0*
     & (cw78%e(2)*p412(3)-p412(2)*cw78%e(3))-p412k0*(cw78%e(2)*p
     & 356(3)-p356(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p356k0+p356(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p412k0+p412(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p356(0)-cw78%e(1)*p356(1)-cw78%e(2)*p356(2)
     & -cw78%e(3)*p356(3)
      cvqd=cw78%e(0)*p412(0)-cw78%e(1)*p412(1)-cw78%e(2)*p412(2)
     & -cw78%e(3)*p412(3)
      cauxa=-cw78%ek0*quqd+p356k0*cvqd+p412k0*cvqu
      cauxb=-cw78%ek0*p412(2)+p412k0*cw78%e(2)
      cauxc=+cw78%ek0*p356(2)-p356k0*cw78%e(2)
      u356_78%a(2)=ccl*(cauxa-ceps_0)
      u356_78%b(1)=ccl*(cauxb-ceps_2)
      u356_78%c(2)=ccl*(-cauxc+ceps_1)
      u356_78%d(1)=ccl*cw78%ek0
  
* TLT0_W -- aa=l3_5678%a,cc=l3_5678%c,a1=l3_56%a,c1=l3_56%c,a2=u356_78%a
* ,b2=u356_78%b,c2=u356_78%c,d2=u356_78%d,prq=s356,nsum=0
      l3_5678%c(2)=l3_56%c(2)*s356*u356_78%d(1)+l3_56%a(2)*u356_
     & 78%c(2)
      l3_5678%a(2)=l3_56%c(2)*s356*u356_78%b(1)+l3_56%a(2)*u356_
     & 78%a(2)
        else
* quqd -- p=p378,q=p412
      quqd=p378(0)*p412(0)-p378(1)*p412(1)-p378(2)*p412(2)-p378(
     & 3)*p412(3)
      ccl=wcl/(f412)
* TW0 -- qu=p378,qd=p412,v=cw56%e,a=u378_56%a,b=u378_56%b,c=u378_56%c,d=
* u378_56%d,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p378(2)*p412(3)-p412(2)*p378(3))+p378k0*
     & (cw56%e(2)*p412(3)-p412(2)*cw56%e(3))-p412k0*(cw56%e(2)*p
     & 378(3)-p378(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p378k0+p378(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p412k0+p412(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p378(0)-cw56%e(1)*p378(1)-cw56%e(2)*p378(2)
     & -cw56%e(3)*p378(3)
      cvqd=cw56%e(0)*p412(0)-cw56%e(1)*p412(1)-cw56%e(2)*p412(2)
     & -cw56%e(3)*p412(3)
      cauxa=-cw56%ek0*quqd+p378k0*cvqd+p412k0*cvqu
      cauxb=-cw56%ek0*p412(2)+p412k0*cw56%e(2)
      cauxc=+cw56%ek0*p378(2)-p378k0*cw56%e(2)
      u378_56%a(2)=ccl*(cauxa-ceps_0)
      u378_56%b(1)=ccl*(cauxb-ceps_2)
      u378_56%c(2)=ccl*(-cauxc+ceps_1)
      u378_56%d(1)=ccl*cw56%ek0
  
* TLT0_W -- aa=l3_5678%a,cc=l3_5678%c,a1=l3_78%a,c1=l3_78%c,a2=u378_56%a
* ,b2=u378_56%b,c2=u378_56%c,d2=u378_56%d,prq=s378,nsum=0
      l3_5678%c(2)=l3_78%c(2)*s378*u378_56%d(1)+l3_78%a(2)*u378_
     & 56%c(2)
      l3_5678%a(2)=l3_78%c(2)*s378*u378_56%b(1)+l3_78%a(2)*u378_
     & 56%a(2)
       endif
  
* III                                                                   
       if (iup(id3).ne.1) then
      do i2=1,2
* TTR0_W -- aa=r4_278(i2)%a,bb=r4_278(i2)%b,a1=u3516_2(i2)%a,b1=u3516_2(
* i2)%b,c1=u3516_2(i2)%c,d1=u3516_2(i2)%d,a2=r4_78%a,b2=r4_78%b,prq=s478
* ,nsum=0
      r4_278(i2)%b(1)=u3516_2(i2)%d(1)*s478*r4_78%b(1)+u3516_2(i
     & 2)%b(1)*r4_78%a(2)
      r4_278(i2)%a(2)=u3516_2(i2)%c(2)*s478*r4_78%b(1)+u3516_2(i
     & 2)%a(2)*r4_78%a(2)
      end do
  
      do i2=1,2
* TTR0_W -- aa=r4_278(i2)%a,bb=r4_278(i2)%b,a1=u3516_78%a,b1=u3516_78%b,
* c1=u3516_78%c,d1=u3516_78%d,a2=r4_2(i2)%a,b2=r4_2(i2)%b,prq=s42,nsum=1
      r4_278(i2)%b(1)=r4_278(i2)%b(1)+u3516_78%d(1)*s42*r4_2(i2)
     & %b(1)+u3516_78%b(1)*r4_2(i2)%a(2)
      r4_278(i2)%a(2)=r4_278(i2)%a(2)+u3516_78%c(2)*s42*r4_2(i2)
     & %b(1)+u3516_78%a(2)*r4_2(i2)%a(2)
      end do
       else
      do i2=1,2
* TTR0_W -- aa=r4_256(i2)%a,bb=r4_256(i2)%b,a1=u3718_2(i2)%a,b1=u3718_2(
* i2)%b,c1=u3718_2(i2)%c,d1=u3718_2(i2)%d,a2=r4_56%a,b2=r4_56%b,prq=s456
* ,nsum=0
      r4_256(i2)%b(1)=u3718_2(i2)%d(1)*s456*r4_56%b(1)+u3718_2(i
     & 2)%b(1)*r4_56%a(2)
      r4_256(i2)%a(2)=u3718_2(i2)%c(2)*s456*r4_56%b(1)+u3718_2(i
     & 2)%a(2)*r4_56%a(2)
      end do
  
      do i2=1,2
* TTR0_W -- aa=r4_256(i2)%a,bb=r4_256(i2)%b,a1=u3718_56%a,b1=u3718_56%b,
* c1=u3718_56%c,d1=u3718_56%d,a2=r4_2(i2)%a,b2=r4_2(i2)%b,prq=s42,nsum=1
      r4_256(i2)%b(1)=r4_256(i2)%b(1)+u3718_56%d(1)*s42*r4_2(i2)
     & %b(1)+u3718_56%b(1)*r4_2(i2)%a(2)
      r4_256(i2)%a(2)=r4_256(i2)%a(2)+u3718_56%c(2)*s42*r4_2(i2)
     & %b(1)+u3718_56%a(2)*r4_2(i2)%a(2)
      end do
        endif
       if (iup(id3).ne.1) then
      do i2=1,2
* TTR0_W -- aa=r4_178(i2)%a,bb=r4_178(i2)%b,a1=u3526_1(i2)%a,b1=u3526_1(
* i2)%b,c1=u3526_1(i2)%c,d1=u3526_1(i2)%d,a2=r4_78%a,b2=r4_78%b,prq=s478
* ,nsum=0
      r4_178(i2)%b(1)=u3526_1(i2)%d(1)*s478*r4_78%b(1)+u3526_1(i
     & 2)%b(1)*r4_78%a(2)
      r4_178(i2)%a(2)=u3526_1(i2)%c(2)*s478*r4_78%b(1)+u3526_1(i
     & 2)%a(2)*r4_78%a(2)
      end do
  
      do i2=1,2
* TTR0_W -- aa=r4_178(i2)%a,bb=r4_178(i2)%b,a1=u3526_78%a,b1=u3526_78%b,
* c1=u3526_78%c,d1=u3526_78%d,a2=r4_1(i2)%a,b2=r4_1(i2)%b,prq=s41,nsum=1
      r4_178(i2)%b(1)=r4_178(i2)%b(1)+u3526_78%d(1)*s41*r4_1(i2)
     & %b(1)+u3526_78%b(1)*r4_1(i2)%a(2)
      r4_178(i2)%a(2)=r4_178(i2)%a(2)+u3526_78%c(2)*s41*r4_1(i2)
     & %b(1)+u3526_78%a(2)*r4_1(i2)%a(2)
      end do
       else
      do i2=1,2
* TTR0_W -- aa=r4_156(i2)%a,bb=r4_156(i2)%b,a1=u3728_1(i2)%a,b1=u3728_1(
* i2)%b,c1=u3728_1(i2)%c,d1=u3728_1(i2)%d,a2=r4_56%a,b2=r4_56%b,prq=s456
* ,nsum=0
      r4_156(i2)%b(1)=u3728_1(i2)%d(1)*s456*r4_56%b(1)+u3728_1(i
     & 2)%b(1)*r4_56%a(2)
      r4_156(i2)%a(2)=u3728_1(i2)%c(2)*s456*r4_56%b(1)+u3728_1(i
     & 2)%a(2)*r4_56%a(2)
      end do
  
      do i2=1,2
* TTR0_W -- aa=r4_156(i2)%a,bb=r4_156(i2)%b,a1=u3728_56%a,b1=u3728_56%b,
* c1=u3728_56%c,d1=u3728_56%d,a2=r4_1(i2)%a,b2=r4_1(i2)%b,prq=s41,nsum=1
      r4_156(i2)%b(1)=r4_156(i2)%b(1)+u3728_56%d(1)*s41*r4_1(i2)
     & %b(1)+u3728_56%b(1)*r4_1(i2)%a(2)
      r4_156(i2)%a(2)=r4_156(i2)%a(2)+u3728_56%c(2)*s41*r4_1(i2)
     & %b(1)+u3728_56%a(2)*r4_1(i2)%a(2)
      end do
        endif
  
* IV                                                                    
        if (iup(id3).ne.1) then
* quqd -- p=p356,q=p478
      quqd=p356(0)*p478(0)-p356(1)*p478(1)-p356(2)*p478(2)-p356(
     & 3)*p478(3)
      ccl=1d0/(f478)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p356,qd=p478,v=ctrip12(i1,i2)%e,a=u356_gg(i1,i2)%a,b=u356_gg
* (i1,i2)%b,c=u356_gg(i1,i2)%c,d=u356_gg(i1,i2)%d,cl=ccl,nsum=0
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
      u356_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u356_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u356_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u356_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      end do
      end do
        else
* quqd -- p=p378,q=p456
      quqd=p378(0)*p456(0)-p378(1)*p456(1)-p378(2)*p456(2)-p378(
     & 3)*p456(3)
      ccl=1d0/(f456)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p378,qd=p456,v=ctrip12(i1,i2)%e,a=u378_gg(i1,i2)%a,b=u378_gg
* (i1,i2)%b,c=u378_gg(i1,i2)%c,d=u378_gg(i1,i2)%d,cl=ccl,nsum=0
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
      u378_gg(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      u378_gg(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      u378_gg(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      u378_gg(i1,i2)%d(1)=ccl*ctrip12(i1,i2)%ek0
      end do
      end do
        endif
  
      endif
  
* W line                                                                
  
      if (ilept(id5).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,9),a1=l5_12(i1,i2)%a,c1=l5_12(i1,i2)%c,a2=
* r6_3478(i3)%a,b2=r6_3478(i3)%b,prq=s512,bef=cres(i1,i2,i3,9)+,aft=
      cres(i1,i2,i3,9)=cres(i1,i2,i3,9)+(l5_12(i1,i2)%c(2)*s512*
     & r6_3478(i3)%b(1)+l5_12(i1,i2)%a(2)*r6_3478(i3)%a(2))
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,9),a1=l5_3478(i3)%a,c1=l5_3478(i3)%c,a2=r6
* _12(i1,i2)%a,b2=r6_12(i1,i2)%b,prq=s612,bef=cres(i1,i2,i3,9)+,aft=
      cres(i1,i2,i3,9)=cres(i1,i2,i3,9)+(l5_3478(i3)%c(2)*s612*r
     & 6_12(i1,i2)%b(1)+l5_3478(i3)%a(2)*r6_12(i1,i2)%a(2))
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,9),a1=l5_134(i1,i3)%a,c1=l5_134(i1,i3)%c,a
* 2=r6_278(i2)%a,b2=r6_278(i2)%b,prq=s5134,bef=cres(i1,i2,i3,9)+,aft=
      cres(i1,i2,i3,9)=cres(i1,i2,i3,9)+(l5_134(i1,i3)%c(2)*s513
     & 4*r6_278(i2)%b(1)+l5_134(i1,i3)%a(2)*r6_278(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,9),a1=l5_178(i1)%a,c1=l5_178(i1)%c,a2=r6_2
* 34(i2,i3)%a,b2=r6_234(i2,i3)%b,prq=s5178,bef=cres(i1,i2,i3,9)+,aft=
      cres(i1,i2,i3,9)=cres(i1,i2,i3,9)+(l5_178(i1)%c(2)*s5178*r
     & 6_234(i2,i3)%b(1)+l5_178(i1)%a(2)*r6_234(i2,i3)%a(2))
      end do
      end do
      end do
  
* IV                                                                    
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(i3,i1,i2)%a,cc=laux_iii(i3,i1,i2)%c,a1=l5_34(i3)
* %a,c1=l5_34(i3)%c,a2=u534_gg(i1,i2)%a,b2=u534_gg(i1,i2)%b,c2=u534_gg(i
* 1,i2)%c,d2=u534_gg(i1,i2)%d,prq=s534,nsum=0
      laux_iii(i3,i1,i2)%c(2)=l5_34(i3)%c(2)*s534*u534_gg(i1,i2)
     & %d(1)+l5_34(i3)%a(2)*u534_gg(i1,i2)%c(2)
      laux_iii(i3,i1,i2)%a(2)=l5_34(i3)%c(2)*s534*u534_gg(i1,i2)
     & %b(1)+l5_34(i3)%a(2)*u534_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,9),a1=laux_iii(i3,i1,i2)%a,c1=laux_iii(i3,
* i1,i2)%c,a2=r6_78%a,b2=r6_78%b,prq=s678,bef=cres(i1,i2,i3,9)+,aft=
      cres(i1,i2,i3,9)=cres(i1,i2,i3,9)+(laux_iii(i3,i1,i2)%c(2)
     & *s678*r6_78%b(1)+laux_iii(i3,i1,i2)%a(2)*r6_78%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(1,i1,i2)%a,cc=laux_iii(1,i1,i2)%c,a1=l5_78%a,c1=
* l5_78%c,a2=u578_gg(i1,i2)%a,b2=u578_gg(i1,i2)%b,c2=u578_gg(i1,i2)%c,d2
* =u578_gg(i1,i2)%d,prq=s578,nsum=0
      laux_iii(1,i1,i2)%c(2)=l5_78%c(2)*s578*u578_gg(i1,i2)%d(1)
     & +l5_78%a(2)*u578_gg(i1,i2)%c(2)
      laux_iii(1,i1,i2)%a(2)=l5_78%c(2)*s578*u578_gg(i1,i2)%b(1)
     & +l5_78%a(2)*u578_gg(i1,i2)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,9),a1=laux_iii(1,i1,i2)%a,c1=laux_iii(1,i1
* ,i2)%c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=cres(i1,i2,i3,9)+,af
* t=
      cres(i1,i2,i3,9)=cres(i1,i2,i3,9)+(laux_iii(1,i1,i2)%c(2)*
     & s634*r6_34(i3)%b(1)+laux_iii(1,i1,i2)%a(2)*r6_34(i3)%a(2)
     & )
      end do
      end do
      end do
  
      endif
  
      if (ilept(id5).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,10),a1=l5_21(i1,i2)%a,c1=l5_21(i1,i2)%c,a2
* =r6_3478(i3)%a,b2=r6_3478(i3)%b,prq=s512,bef=cres(i1,i2,i3,10)+,aft=
      cres(i1,i2,i3,10)=cres(i1,i2,i3,10)+(l5_21(i1,i2)%c(2)*s51
     & 2*r6_3478(i3)%b(1)+l5_21(i1,i2)%a(2)*r6_3478(i3)%a(2))
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,10),a1=l5_3478(i3)%a,c1=l5_3478(i3)%c,a2=r
* 6_21(i1,i2)%a,b2=r6_21(i1,i2)%b,prq=s612,bef=cres(i1,i2,i3,10)+,aft=
      cres(i1,i2,i3,10)=cres(i1,i2,i3,10)+(l5_3478(i3)%c(2)*s612
     & *r6_21(i1,i2)%b(1)+l5_3478(i3)%a(2)*r6_21(i1,i2)%a(2))
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,10),a1=l5_234(i2,i3)%a,c1=l5_234(i2,i3)%c,
* a2=r6_178(i1)%a,b2=r6_178(i1)%b,prq=s5234,bef=cres(i1,i2,i3,10)+,aft=
      cres(i1,i2,i3,10)=cres(i1,i2,i3,10)+(l5_234(i2,i3)%c(2)*s5
     & 234*r6_178(i1)%b(1)+l5_234(i2,i3)%a(2)*r6_178(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,10),a1=l5_278(i2)%a,c1=l5_278(i2)%c,a2=r6_
* 134(i1,i3)%a,b2=r6_134(i1,i3)%b,prq=s5278,bef=cres(i1,i2,i3,10)+,aft=
      cres(i1,i2,i3,10)=cres(i1,i2,i3,10)+(l5_278(i2)%c(2)*s5278
     & *r6_134(i1,i3)%b(1)+l5_278(i2)%a(2)*r6_134(i1,i3)%a(2))
      end do
      end do
      end do
  
* IV                                                                    
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(i3,i1,i2)%a,cc=laux_iii(i3,i1,i2)%c,a1=l5_34(i3)
* %a,c1=l5_34(i3)%c,a2=u534_gg(i1,i2)%a,b2=u534_gg(i1,i2)%b,c2=u534_gg(i
* 1,i2)%c,d2=u534_gg(i1,i2)%d,prq=s534,nsum=0
      laux_iii(i3,i1,i2)%c(2)=l5_34(i3)%c(2)*s534*u534_gg(i1,i2)
     & %d(1)+l5_34(i3)%a(2)*u534_gg(i1,i2)%c(2)
      laux_iii(i3,i1,i2)%a(2)=l5_34(i3)%c(2)*s534*u534_gg(i1,i2)
     & %b(1)+l5_34(i3)%a(2)*u534_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,10),a1=laux_iii(i3,i1,i2)%a,c1=laux_iii(i3
* ,i1,i2)%c,a2=r6_78%a,b2=r6_78%b,prq=s678,bef=cres(i1,i2,i3,10)-,aft=
      cres(i1,i2,i3,10)=cres(i1,i2,i3,10)-(laux_iii(i3,i1,i2)%c(
     & 2)*s678*r6_78%b(1)+laux_iii(i3,i1,i2)%a(2)*r6_78%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(1,i1,i2)%a,cc=laux_iii(1,i1,i2)%c,a1=l5_78%a,c1=
* l5_78%c,a2=u578_gg(i1,i2)%a,b2=u578_gg(i1,i2)%b,c2=u578_gg(i1,i2)%c,d2
* =u578_gg(i1,i2)%d,prq=s578,nsum=0
      laux_iii(1,i1,i2)%c(2)=l5_78%c(2)*s578*u578_gg(i1,i2)%d(1)
     & +l5_78%a(2)*u578_gg(i1,i2)%c(2)
      laux_iii(1,i1,i2)%a(2)=l5_78%c(2)*s578*u578_gg(i1,i2)%b(1)
     & +l5_78%a(2)*u578_gg(i1,i2)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,10),a1=laux_iii(1,i1,i2)%a,c1=laux_iii(1,i
* 1,i2)%c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=s634,bef=cres(i1,i2,i3,10)-,
* aft=
      cres(i1,i2,i3,10)=cres(i1,i2,i3,10)-(laux_iii(1,i1,i2)%c(2
     & )*s634*r6_34(i3)%b(1)+laux_iii(1,i1,i2)%a(2)*r6_34(i3)%a(
     & 2))
      end do
      end do
      end do
  
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,11),a1=l7_12(i1,i2)%a,c1=l7_12(i1,i2)%c,a2
* =r8_3456(i3)%a,b2=r8_3456(i3)%b,prq=s712,bef=cres(i1,i2,i3,11)+,aft=
      cres(i1,i2,i3,11)=cres(i1,i2,i3,11)+(l7_12(i1,i2)%c(2)*s71
     & 2*r8_3456(i3)%b(1)+l7_12(i1,i2)%a(2)*r8_3456(i3)%a(2))
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,11),a1=l7_3456(i3)%a,c1=l7_3456(i3)%c,a2=r
* 8_12(i1,i2)%a,b2=r8_12(i1,i2)%b,prq=s812,bef=cres(i1,i2,i3,11)+,aft=
      cres(i1,i2,i3,11)=cres(i1,i2,i3,11)+(l7_3456(i3)%c(2)*s812
     & *r8_12(i1,i2)%b(1)+l7_3456(i3)%a(2)*r8_12(i1,i2)%a(2))
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,11),a1=l7_134(i1,i3)%a,c1=l7_134(i1,i3)%c,
* a2=r8_256(i2)%a,b2=r8_256(i2)%b,prq=s7134,bef=cres(i1,i2,i3,11)+,aft=
      cres(i1,i2,i3,11)=cres(i1,i2,i3,11)+(l7_134(i1,i3)%c(2)*s7
     & 134*r8_256(i2)%b(1)+l7_134(i1,i3)%a(2)*r8_256(i2)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,11),a1=l7_156(i1)%a,c1=l7_156(i1)%c,a2=r8_
* 234(i2,i3)%a,b2=r8_234(i2,i3)%b,prq=s7156,bef=cres(i1,i2,i3,11)+,aft=
      cres(i1,i2,i3,11)=cres(i1,i2,i3,11)+(l7_156(i1)%c(2)*s7156
     & *r8_234(i2,i3)%b(1)+l7_156(i1)%a(2)*r8_234(i2,i3)%a(2))
      end do
      end do
      end do
  
* IV                                                                    
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(i3,i1,i2)%a,cc=laux_iii(i3,i1,i2)%c,a1=l7_34(i3)
* %a,c1=l7_34(i3)%c,a2=u734_gg(i1,i2)%a,b2=u734_gg(i1,i2)%b,c2=u734_gg(i
* 1,i2)%c,d2=u734_gg(i1,i2)%d,prq=s734,nsum=0
      laux_iii(i3,i1,i2)%c(2)=l7_34(i3)%c(2)*s734*u734_gg(i1,i2)
     & %d(1)+l7_34(i3)%a(2)*u734_gg(i1,i2)%c(2)
      laux_iii(i3,i1,i2)%a(2)=l7_34(i3)%c(2)*s734*u734_gg(i1,i2)
     & %b(1)+l7_34(i3)%a(2)*u734_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,11),a1=laux_iii(i3,i1,i2)%a,c1=laux_iii(i3
* ,i1,i2)%c,a2=r8_56%a,b2=r8_56%b,prq=s856,bef=cres(i1,i2,i3,11)+,aft=
      cres(i1,i2,i3,11)=cres(i1,i2,i3,11)+(laux_iii(i3,i1,i2)%c(
     & 2)*s856*r8_56%b(1)+laux_iii(i3,i1,i2)%a(2)*r8_56%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(1,i1,i2)%a,cc=laux_iii(1,i1,i2)%c,a1=l7_56%a,c1=
* l7_56%c,a2=u756_gg(i1,i2)%a,b2=u756_gg(i1,i2)%b,c2=u756_gg(i1,i2)%c,d2
* =u756_gg(i1,i2)%d,prq=s756,nsum=0
      laux_iii(1,i1,i2)%c(2)=l7_56%c(2)*s756*u756_gg(i1,i2)%d(1)
     & +l7_56%a(2)*u756_gg(i1,i2)%c(2)
      laux_iii(1,i1,i2)%a(2)=l7_56%c(2)*s756*u756_gg(i1,i2)%b(1)
     & +l7_56%a(2)*u756_gg(i1,i2)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,11),a1=laux_iii(1,i1,i2)%a,c1=laux_iii(1,i
* 1,i2)%c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=cres(i1,i2,i3,11)+,
* aft=
      cres(i1,i2,i3,11)=cres(i1,i2,i3,11)+(laux_iii(1,i1,i2)%c(2
     & )*s834*r8_34(i3)%b(1)+laux_iii(1,i1,i2)%a(2)*r8_34(i3)%a(
     & 2))
      end do
      end do
      end do
  
      endif
  
      if (ilept(id7).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,12),a1=l7_21(i1,i2)%a,c1=l7_21(i1,i2)%c,a2
* =r8_3456(i3)%a,b2=r8_3456(i3)%b,prq=s712,bef=cres(i1,i2,i3,12)+,aft=
      cres(i1,i2,i3,12)=cres(i1,i2,i3,12)+(l7_21(i1,i2)%c(2)*s71
     & 2*r8_3456(i3)%b(1)+l7_21(i1,i2)%a(2)*r8_3456(i3)%a(2))
      end do
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,12),a1=l7_3456(i3)%a,c1=l7_3456(i3)%c,a2=r
* 8_21(i1,i2)%a,b2=r8_21(i1,i2)%b,prq=s812,bef=cres(i1,i2,i3,12)+,aft=
      cres(i1,i2,i3,12)=cres(i1,i2,i3,12)+(l7_3456(i3)%c(2)*s812
     & *r8_21(i1,i2)%b(1)+l7_3456(i3)%a(2)*r8_21(i1,i2)%a(2))
      end do
      end do
      end do
  
* III                                                                   
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,12),a1=l7_234(i2,i3)%a,c1=l7_234(i2,i3)%c,
* a2=r8_156(i1)%a,b2=r8_156(i1)%b,prq=s7234,bef=cres(i1,i2,i3,12)+,aft=
      cres(i1,i2,i3,12)=cres(i1,i2,i3,12)+(l7_234(i2,i3)%c(2)*s7
     & 234*r8_156(i1)%b(1)+l7_234(i2,i3)%a(2)*r8_156(i1)%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,12),a1=l7_256(i2)%a,c1=l7_256(i2)%c,a2=r8_
* 134(i1,i3)%a,b2=r8_134(i1,i3)%b,prq=s7256,bef=cres(i1,i2,i3,12)+,aft=
      cres(i1,i2,i3,12)=cres(i1,i2,i3,12)+(l7_256(i2)%c(2)*s7256
     & *r8_134(i1,i3)%b(1)+l7_256(i2)%a(2)*r8_134(i1,i3)%a(2))
      end do
      end do
      end do
  
* IV                                                                    
      do i3=1,2
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(i3,i1,i2)%a,cc=laux_iii(i3,i1,i2)%c,a1=l7_34(i3)
* %a,c1=l7_34(i3)%c,a2=u734_gg(i1,i2)%a,b2=u734_gg(i1,i2)%b,c2=u734_gg(i
* 1,i2)%c,d2=u734_gg(i1,i2)%d,prq=s734,nsum=0
      laux_iii(i3,i1,i2)%c(2)=l7_34(i3)%c(2)*s734*u734_gg(i1,i2)
     & %d(1)+l7_34(i3)%a(2)*u734_gg(i1,i2)%c(2)
      laux_iii(i3,i1,i2)%a(2)=l7_34(i3)%c(2)*s734*u734_gg(i1,i2)
     & %b(1)+l7_34(i3)%a(2)*u734_gg(i1,i2)%a(2)
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,12),a1=laux_iii(i3,i1,i2)%a,c1=laux_iii(i3
* ,i1,i2)%c,a2=r8_56%a,b2=r8_56%b,prq=s856,bef=cres(i1,i2,i3,12)-,aft=
      cres(i1,i2,i3,12)=cres(i1,i2,i3,12)-(laux_iii(i3,i1,i2)%c(
     & 2)*s856*r8_56%b(1)+laux_iii(i3,i1,i2)%a(2)*r8_56%a(2))
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(1,i1,i2)%a,cc=laux_iii(1,i1,i2)%c,a1=l7_56%a,c1=
* l7_56%c,a2=u756_gg(i1,i2)%a,b2=u756_gg(i1,i2)%b,c2=u756_gg(i1,i2)%c,d2
* =u756_gg(i1,i2)%d,prq=s756,nsum=0
      laux_iii(1,i1,i2)%c(2)=l7_56%c(2)*s756*u756_gg(i1,i2)%d(1)
     & +l7_56%a(2)*u756_gg(i1,i2)%c(2)
      laux_iii(1,i1,i2)%a(2)=l7_56%c(2)*s756*u756_gg(i1,i2)%b(1)
     & +l7_56%a(2)*u756_gg(i1,i2)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
* TLTR0_W -- aa=cres(i1,i2,i3,12),a1=laux_iii(1,i1,i2)%a,c1=laux_iii(1,i
* 1,i2)%c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=s834,bef=cres(i1,i2,i3,12)-,
* aft=
      cres(i1,i2,i3,12)=cres(i1,i2,i3,12)-(laux_iii(1,i1,i2)%c(2
     & )*s834*r8_34(i3)%b(1)+laux_iii(1,i1,i2)%a(2)*r8_34(i3)%a(
     & 2))
      end do
      end do
      end do
  
      endif
  
* Z line                                                                
  
      if (ilept(id3).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,7),a1=l3_12(i1,i2)%a,c1=l3_12(i1,i2)%c,a2=r
* 4_5678%a,b2=r4_5678%b,prq=s312,bef=cres(i1,i2,2,7)+,aft=
      cres(i1,i2,2,7)=cres(i1,i2,2,7)+(l3_12(i1,i2)%c(2)*s312*r4
     & _5678%b(1)+l3_12(i1,i2)%a(2)*r4_5678%a(2))
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,7),a1=l3_5678%a,c1=l3_5678%c,a2=r4_12(i1,i2
* )%a,b2=r4_12(i1,i2)%b,prq=s412,bef=cres(i1,i2,2,7)+,aft=
      cres(i1,i2,2,7)=cres(i1,i2,2,7)+(l3_5678%c(2)*s412*r4_12(i
     & 1,i2)%b(1)+l3_5678%a(2)*r4_12(i1,i2)%a(2))
      end do
      end do
  
* III                                                                   
        if (iup(id3).ne.1) then
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,7),a1=l3_156(i1)%a,c1=l3_156(i1)%c,a2=r4_27
* 8(i2)%a,b2=r4_278(i2)%b,prq=s3156,bef=cres(i1,i2,2,7)+,aft=
      cres(i1,i2,2,7)=cres(i1,i2,2,7)+(l3_156(i1)%c(2)*s3156*r4_
     & 278(i2)%b(1)+l3_156(i1)%a(2)*r4_278(i2)%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,7),a1=l3_178(i1)%a,c1=l3_178(i1)%c,a2=r4_25
* 6(i2)%a,b2=r4_256(i2)%b,prq=s3178,bef=cres(i1,i2,2,7)+,aft=
      cres(i1,i2,2,7)=cres(i1,i2,2,7)+(l3_178(i1)%c(2)*s3178*r4_
     & 256(i2)%b(1)+l3_178(i1)%a(2)*r4_256(i2)%a(2))
      end do
      end do
        endif
* IV                                                                    
        if (iup(id3).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(2,i1,i2)%a,cc=laux_iii(2,i1,i2)%c,a1=l3_56%a,c1=
* l3_56%c,a2=u356_gg(i1,i2)%a,b2=u356_gg(i1,i2)%b,c2=u356_gg(i1,i2)%c,d2
* =u356_gg(i1,i2)%d,prq=s356,nsum=0
      laux_iii(2,i1,i2)%c(2)=l3_56%c(2)*s356*u356_gg(i1,i2)%d(1)
     & +l3_56%a(2)*u356_gg(i1,i2)%c(2)
      laux_iii(2,i1,i2)%a(2)=l3_56%c(2)*s356*u356_gg(i1,i2)%b(1)
     & +l3_56%a(2)*u356_gg(i1,i2)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,7),a1=laux_iii(2,i1,i2)%a,c1=laux_iii(2,i1,
* i2)%c,a2=r4_78%a,b2=r4_78%b,prq=s478,bef=cres(i1,i2,2,7)+,aft=
      cres(i1,i2,2,7)=cres(i1,i2,2,7)+(laux_iii(2,i1,i2)%c(2)*s4
     & 78*r4_78%b(1)+laux_iii(2,i1,i2)%a(2)*r4_78%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(1,i1,i2)%a,cc=laux_iii(1,i1,i2)%c,a1=l3_78%a,c1=
* l3_78%c,a2=u378_gg(i1,i2)%a,b2=u378_gg(i1,i2)%b,c2=u378_gg(i1,i2)%c,d2
* =u378_gg(i1,i2)%d,prq=s378,nsum=0
      laux_iii(1,i1,i2)%c(2)=l3_78%c(2)*s378*u378_gg(i1,i2)%d(1)
     & +l3_78%a(2)*u378_gg(i1,i2)%c(2)
      laux_iii(1,i1,i2)%a(2)=l3_78%c(2)*s378*u378_gg(i1,i2)%b(1)
     & +l3_78%a(2)*u378_gg(i1,i2)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,7),a1=laux_iii(1,i1,i2)%a,c1=laux_iii(1,i1,
* i2)%c,a2=r4_56%a,b2=r4_56%b,prq=s456,bef=cres(i1,i2,2,7)+,aft=
      cres(i1,i2,2,7)=cres(i1,i2,2,7)+(laux_iii(1,i1,i2)%c(2)*s4
     & 56*r4_56%b(1)+laux_iii(1,i1,i2)%a(2)*r4_56%a(2))
      end do
      end do
        endif
  
      endif
  
      if (ilept(id3).ne.1) then
* I                                                                     
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,8),a1=l3_21(i1,i2)%a,c1=l3_21(i1,i2)%c,a2=r
* 4_5678%a,b2=r4_5678%b,prq=s312,bef=cres(i1,i2,2,8)+,aft=
      cres(i1,i2,2,8)=cres(i1,i2,2,8)+(l3_21(i1,i2)%c(2)*s312*r4
     & _5678%b(1)+l3_21(i1,i2)%a(2)*r4_5678%a(2))
      end do
      end do
* II                                                                    
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,8),a1=l3_5678%a,c1=l3_5678%c,a2=r4_21(i1,i2
* )%a,b2=r4_21(i1,i2)%b,prq=s412,bef=cres(i1,i2,2,8)+,aft=
      cres(i1,i2,2,8)=cres(i1,i2,2,8)+(l3_5678%c(2)*s412*r4_21(i
     & 1,i2)%b(1)+l3_5678%a(2)*r4_21(i1,i2)%a(2))
      end do
      end do
  
* III                                                                   
        if (iup(id3).ne.1) then
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,8),a1=l3_256(i2)%a,c1=l3_256(i2)%c,a2=r4_17
* 8(i1)%a,b2=r4_178(i1)%b,prq=s3256,bef=cres(i1,i2,2,8)+,aft=
      cres(i1,i2,2,8)=cres(i1,i2,2,8)+(l3_256(i2)%c(2)*s3256*r4_
     & 178(i1)%b(1)+l3_256(i2)%a(2)*r4_178(i1)%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,8),a1=l3_278(i2)%a,c1=l3_278(i2)%c,a2=r4_15
* 6(i1)%a,b2=r4_156(i1)%b,prq=s3278,bef=cres(i1,i2,2,8)+,aft=
      cres(i1,i2,2,8)=cres(i1,i2,2,8)+(l3_278(i2)%c(2)*s3278*r4_
     & 156(i1)%b(1)+l3_278(i2)%a(2)*r4_156(i1)%a(2))
      end do
      end do
        endif
* IV                                                                    
        if (iup(id3).ne.1) then
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(2,i1,i2)%a,cc=laux_iii(2,i1,i2)%c,a1=l3_56%a,c1=
* l3_56%c,a2=u356_gg(i1,i2)%a,b2=u356_gg(i1,i2)%b,c2=u356_gg(i1,i2)%c,d2
* =u356_gg(i1,i2)%d,prq=s356,nsum=0
      laux_iii(2,i1,i2)%c(2)=l3_56%c(2)*s356*u356_gg(i1,i2)%d(1)
     & +l3_56%a(2)*u356_gg(i1,i2)%c(2)
      laux_iii(2,i1,i2)%a(2)=l3_56%c(2)*s356*u356_gg(i1,i2)%b(1)
     & +l3_56%a(2)*u356_gg(i1,i2)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,8),a1=laux_iii(2,i1,i2)%a,c1=laux_iii(2,i1,
* i2)%c,a2=r4_78%a,b2=r4_78%b,prq=s478,bef=cres(i1,i2,2,8)-,aft=
      cres(i1,i2,2,8)=cres(i1,i2,2,8)-(laux_iii(2,i1,i2)%c(2)*s4
     & 78*r4_78%b(1)+laux_iii(2,i1,i2)%a(2)*r4_78%a(2))
      end do
      end do
        else
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=laux_iii(1,i1,i2)%a,cc=laux_iii(1,i1,i2)%c,a1=l3_78%a,c1=
* l3_78%c,a2=u378_gg(i1,i2)%a,b2=u378_gg(i1,i2)%b,c2=u378_gg(i1,i2)%c,d2
* =u378_gg(i1,i2)%d,prq=s378,nsum=0
      laux_iii(1,i1,i2)%c(2)=l3_78%c(2)*s378*u378_gg(i1,i2)%d(1)
     & +l3_78%a(2)*u378_gg(i1,i2)%c(2)
      laux_iii(1,i1,i2)%a(2)=l3_78%c(2)*s378*u378_gg(i1,i2)%b(1)
     & +l3_78%a(2)*u378_gg(i1,i2)%a(2)
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TLTR0_W -- aa=cres(i1,i2,2,8),a1=laux_iii(1,i1,i2)%a,c1=laux_iii(1,i1,
* i2)%c,a2=r4_56%a,b2=r4_56%b,prq=s456,bef=cres(i1,i2,2,8)-,aft=
      cres(i1,i2,2,8)=cres(i1,i2,2,8)-(laux_iii(1,i1,i2)%c(2)*s4
     & 56*r4_56%b(1)+laux_iii(1,i1,i2)%a(2)*r4_56%a(2))
      end do
      end do
        endif
  
      endif
  
* Triple Vertices                                                       
  
  
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
*** triple vertex -- pfz(mu)=p314(mu),pwm(mu)=p78(mu),pwp(mu)=p526(mu),e
* fz=cz314(i3,i1),ewm=cw78,ewp=cw526(i2),res=cres_aux(i1?,i2?,i3?,1)
      do mu=0,3
      vfz(mu)=p78(mu)-p526(mu)
      vwm(mu)=p526(mu)-p314(mu)
      vwp(mu)=p314(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i1=1,2
* p.q -- p.q=cz314(i3,i1)%v,p=cz314(i3,i1)%e,q=vfz,bef=,aft=
      cz314(i3,i1)%v=(cz314(i3,i1)%e(0)*vfz(0)-cz314(i3,i1)%e(1)
     & *vfz(1)-cz314(i3,i1)%e(2)*vfz(2)-cz314(i3,i1)%e(3)*vfz(3)
     & )
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
      do i2=1,2
* p.q -- p.q=cw526(i2)%v,p=cw526(i2)%e,q=vwp,bef=,aft=
      cw526(i2)%v=(cw526(i2)%e(0)*vwp(0)-cw526(i2)%e(1)*vwp(1)-c
     & w526(i2)%e(2)*vwp(2)-cw526(i2)%e(3)*vwp(3))
      end do
* efz%ewm
      do i3=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cz314(i3,i1)%e,q=cw78%e,bef=,aft=
      caux=(cz314(i3,i1)%e(0)*cw78%e(0)-cz314(i3,i1)%e(1)*cw78%e
     & (1)-cz314(i3,i1)%e(2)*cw78%e(2)-cz314(i3,i1)%e(3)*cw78%e(
     & 3))
      do i2=1,2
      cres_aux(i1,i2,i3,1)=cw526(i2)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz314(i3,i1)%e,q=cw526(i2)%e,bef=,aft=
      caux=(cz314(i3,i1)%e(0)*cw526(i2)%e(0)-cz314(i3,i1)%e(1)*c
     & w526(i2)%e(1)-cz314(i3,i1)%e(2)*cw526(i2)%e(2)-cz314(i3,i
     & 1)%e(3)*cw526(i2)%e(3))
      cres_aux(i1,i2,i3,1)=cres_aux(i1,i2,i3,1)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i2=1,2
* p.q -- p.q=caux,p=cw78%e,q=cw526(i2)%e,bef=,aft=
      caux=(cw78%e(0)*cw526(i2)%e(0)-cw78%e(1)*cw526(i2)%e(1)-cw
     & 78%e(2)*cw526(i2)%e(2)-cw78%e(3)*cw526(i2)%e(3))
      do i1=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,1)=cres_aux(i1,i2,i3,1)+cz314(i3,i1)%v*c
     & aux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,1)=cres(i1,i2,i3,1)+
     &      rcotw*cres_aux(i1,i2,i3,1)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
*** triple vertex -- pfz(mu)=p314(mu),pwm(mu)=p78(mu),pwp(mu)=p526(mu),e
* fz=cf314(i3,i1),ewm=cw78,ewp=cw526(i2),res=cres_aux(i1?,i2?,i3?,1)
      do mu=0,3
      vfz(mu)=p78(mu)-p526(mu)
      vwm(mu)=p526(mu)-p314(mu)
      vwp(mu)=p314(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i1=1,2
* p.q -- p.q=cf314(i3,i1)%v,p=cf314(i3,i1)%e,q=vfz,bef=,aft=
      cf314(i3,i1)%v=(cf314(i3,i1)%e(0)*vfz(0)-cf314(i3,i1)%e(1)
     & *vfz(1)-cf314(i3,i1)%e(2)*vfz(2)-cf314(i3,i1)%e(3)*vfz(3)
     & )
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
      do i2=1,2
* p.q -- p.q=cw526(i2)%v,p=cw526(i2)%e,q=vwp,bef=,aft=
      cw526(i2)%v=(cw526(i2)%e(0)*vwp(0)-cw526(i2)%e(1)*vwp(1)-c
     & w526(i2)%e(2)*vwp(2)-cw526(i2)%e(3)*vwp(3))
      end do
* efz%ewm
      do i3=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cf314(i3,i1)%e,q=cw78%e,bef=,aft=
      caux=(cf314(i3,i1)%e(0)*cw78%e(0)-cf314(i3,i1)%e(1)*cw78%e
     & (1)-cf314(i3,i1)%e(2)*cw78%e(2)-cf314(i3,i1)%e(3)*cw78%e(
     & 3))
      do i2=1,2
      cres_aux(i1,i2,i3,1)=cw526(i2)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf314(i3,i1)%e,q=cw526(i2)%e,bef=,aft=
      caux=(cf314(i3,i1)%e(0)*cw526(i2)%e(0)-cf314(i3,i1)%e(1)*c
     & w526(i2)%e(1)-cf314(i3,i1)%e(2)*cw526(i2)%e(2)-cf314(i3,i
     & 1)%e(3)*cw526(i2)%e(3))
      cres_aux(i1,i2,i3,1)=cres_aux(i1,i2,i3,1)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i2=1,2
* p.q -- p.q=caux,p=cw78%e,q=cw526(i2)%e,bef=,aft=
      caux=(cw78%e(0)*cw526(i2)%e(0)-cw78%e(1)*cw526(i2)%e(1)-cw
     & 78%e(2)*cw526(i2)%e(2)-cw78%e(3)*cw526(i2)%e(3))
      do i1=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,1)=cres_aux(i1,i2,i3,1)+cf314(i3,i1)%v*c
     & aux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,1)=cres(i1,i2,i3,1)+
     &      cres_aux(i1,i2,i3,1)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
*** triple vertex -- pfz(mu)=p324(mu),pwm(mu)=p78(mu),pwp(mu)=p516(mu),e
* fz=cz324(i3,i2),ewm=cw78,ewp=cw516(i1),res=cres_aux(i1?,i2?,i3?,2)
      do mu=0,3
      vfz(mu)=p78(mu)-p516(mu)
      vwm(mu)=p516(mu)-p324(mu)
      vwp(mu)=p324(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i2=1,2
* p.q -- p.q=cz324(i3,i2)%v,p=cz324(i3,i2)%e,q=vfz,bef=,aft=
      cz324(i3,i2)%v=(cz324(i3,i2)%e(0)*vfz(0)-cz324(i3,i2)%e(1)
     & *vfz(1)-cz324(i3,i2)%e(2)*vfz(2)-cz324(i3,i2)%e(3)*vfz(3)
     & )
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
      do i1=1,2
* p.q -- p.q=cw516(i1)%v,p=cw516(i1)%e,q=vwp,bef=,aft=
      cw516(i1)%v=(cw516(i1)%e(0)*vwp(0)-cw516(i1)%e(1)*vwp(1)-c
     & w516(i1)%e(2)*vwp(2)-cw516(i1)%e(3)*vwp(3))
      end do
* efz%ewm
      do i3=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz324(i3,i2)%e,q=cw78%e,bef=,aft=
      caux=(cz324(i3,i2)%e(0)*cw78%e(0)-cz324(i3,i2)%e(1)*cw78%e
     & (1)-cz324(i3,i2)%e(2)*cw78%e(2)-cz324(i3,i2)%e(3)*cw78%e(
     & 3))
      do i1=1,2
      cres_aux(i1,i2,i3,2)=cw516(i1)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i2=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cz324(i3,i2)%e,q=cw516(i1)%e,bef=,aft=
      caux=(cz324(i3,i2)%e(0)*cw516(i1)%e(0)-cz324(i3,i2)%e(1)*c
     & w516(i1)%e(1)-cz324(i3,i2)%e(2)*cw516(i1)%e(2)-cz324(i3,i
     & 2)%e(3)*cw516(i1)%e(3))
      cres_aux(i1,i2,i3,2)=cres_aux(i1,i2,i3,2)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
* p.q -- p.q=caux,p=cw78%e,q=cw516(i1)%e,bef=,aft=
      caux=(cw78%e(0)*cw516(i1)%e(0)-cw78%e(1)*cw516(i1)%e(1)-cw
     & 78%e(2)*cw516(i1)%e(2)-cw78%e(3)*cw516(i1)%e(3))
      do i2=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,2)=cres_aux(i1,i2,i3,2)+cz324(i3,i2)%v*c
     & aux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,2)=cres(i1,i2,i3,2)+
     &      rcotw*cres_aux(i1,i2,i3,2)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
*** triple vertex -- pfz(mu)=p324(mu),pwm(mu)=p78(mu),pwp(mu)=p516(mu),e
* fz=cf324(i3,i2),ewm=cw78,ewp=cw516(i1),res=cres_aux(i1?,i2?,i3?,2)
      do mu=0,3
      vfz(mu)=p78(mu)-p516(mu)
      vwm(mu)=p516(mu)-p324(mu)
      vwp(mu)=p324(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i2=1,2
* p.q -- p.q=cf324(i3,i2)%v,p=cf324(i3,i2)%e,q=vfz,bef=,aft=
      cf324(i3,i2)%v=(cf324(i3,i2)%e(0)*vfz(0)-cf324(i3,i2)%e(1)
     & *vfz(1)-cf324(i3,i2)%e(2)*vfz(2)-cf324(i3,i2)%e(3)*vfz(3)
     & )
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
      do i1=1,2
* p.q -- p.q=cw516(i1)%v,p=cw516(i1)%e,q=vwp,bef=,aft=
      cw516(i1)%v=(cw516(i1)%e(0)*vwp(0)-cw516(i1)%e(1)*vwp(1)-c
     & w516(i1)%e(2)*vwp(2)-cw516(i1)%e(3)*vwp(3))
      end do
* efz%ewm
      do i3=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf324(i3,i2)%e,q=cw78%e,bef=,aft=
      caux=(cf324(i3,i2)%e(0)*cw78%e(0)-cf324(i3,i2)%e(1)*cw78%e
     & (1)-cf324(i3,i2)%e(2)*cw78%e(2)-cf324(i3,i2)%e(3)*cw78%e(
     & 3))
      do i1=1,2
      cres_aux(i1,i2,i3,2)=cw516(i1)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i2=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cf324(i3,i2)%e,q=cw516(i1)%e,bef=,aft=
      caux=(cf324(i3,i2)%e(0)*cw516(i1)%e(0)-cf324(i3,i2)%e(1)*c
     & w516(i1)%e(1)-cf324(i3,i2)%e(2)*cw516(i1)%e(2)-cf324(i3,i
     & 2)%e(3)*cw516(i1)%e(3))
      cres_aux(i1,i2,i3,2)=cres_aux(i1,i2,i3,2)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
* p.q -- p.q=caux,p=cw78%e,q=cw516(i1)%e,bef=,aft=
      caux=(cw78%e(0)*cw516(i1)%e(0)-cw78%e(1)*cw516(i1)%e(1)-cw
     & 78%e(2)*cw516(i1)%e(2)-cw78%e(3)*cw516(i1)%e(3))
      do i2=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,2)=cres_aux(i1,i2,i3,2)+cf324(i3,i2)%v*c
     & aux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,2)=cres(i1,i2,i3,2)+
     &      cres_aux(i1,i2,i3,2)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p314(mu),pwm(mu)=p728(mu),pwp(mu)=p56(mu),e
* fz=cz314(i3,i1),ewm=cw728(i2),ewp=cw56,res=cres_aux(i1?,i2?,i3?,3)
      do mu=0,3
      vfz(mu)=p728(mu)-p56(mu)
      vwm(mu)=p56(mu)-p314(mu)
      vwp(mu)=p314(mu)-p728(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i1=1,2
* p.q -- p.q=cz314(i3,i1)%v,p=cz314(i3,i1)%e,q=vfz,bef=,aft=
      cz314(i3,i1)%v=(cz314(i3,i1)%e(0)*vfz(0)-cz314(i3,i1)%e(1)
     & *vfz(1)-cz314(i3,i1)%e(2)*vfz(2)-cz314(i3,i1)%e(3)*vfz(3)
     & )
      end do
      end do
* vwm%ewm
      do i2=1,2
* p.q -- p.q=cw728(i2)%v,p=cw728(i2)%e,q=vwm,bef=,aft=
      cw728(i2)%v=(cw728(i2)%e(0)*vwm(0)-cw728(i2)%e(1)*vwm(1)-c
     & w728(i2)%e(2)*vwm(2)-cw728(i2)%e(3)*vwm(3))
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz314(i3,i1)%e,q=cw728(i2)%e,bef=,aft=
      caux=(cz314(i3,i1)%e(0)*cw728(i2)%e(0)-cz314(i3,i1)%e(1)*c
     & w728(i2)%e(1)-cz314(i3,i1)%e(2)*cw728(i2)%e(2)-cz314(i3,i
     & 1)%e(3)*cw728(i2)%e(3))
      cres_aux(i1,i2,i3,3)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cz314(i3,i1)%e,q=cw56%e,bef=,aft=
      caux=(cz314(i3,i1)%e(0)*cw56%e(0)-cz314(i3,i1)%e(1)*cw56%e
     & (1)-cz314(i3,i1)%e(2)*cw56%e(2)-cz314(i3,i1)%e(3)*cw56%e(
     & 3))
      do i2=1,2
      cres_aux(i1,i2,i3,3)=cres_aux(i1,i2,i3,3)+cw728(i2)%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i2=1,2
* p.q -- p.q=caux,p=cw728(i2)%e,q=cw56%e,bef=,aft=
      caux=(cw728(i2)%e(0)*cw56%e(0)-cw728(i2)%e(1)*cw56%e(1)-cw
     & 728(i2)%e(2)*cw56%e(2)-cw728(i2)%e(3)*cw56%e(3))
      do i1=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,3)=cres_aux(i1,i2,i3,3)+cz314(i3,i1)%v*c
     & aux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,3)=cres(i1,i2,i3,3)+
     &      rcotw*cres_aux(i1,i2,i3,3)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p314(mu),pwm(mu)=p728(mu),pwp(mu)=p56(mu),e
* fz=cf314(i3,i1),ewm=cw728(i2),ewp=cw56,res=cres_aux(i1?,i2?,i3?,3)
      do mu=0,3
      vfz(mu)=p728(mu)-p56(mu)
      vwm(mu)=p56(mu)-p314(mu)
      vwp(mu)=p314(mu)-p728(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i1=1,2
* p.q -- p.q=cf314(i3,i1)%v,p=cf314(i3,i1)%e,q=vfz,bef=,aft=
      cf314(i3,i1)%v=(cf314(i3,i1)%e(0)*vfz(0)-cf314(i3,i1)%e(1)
     & *vfz(1)-cf314(i3,i1)%e(2)*vfz(2)-cf314(i3,i1)%e(3)*vfz(3)
     & )
      end do
      end do
* vwm%ewm
      do i2=1,2
* p.q -- p.q=cw728(i2)%v,p=cw728(i2)%e,q=vwm,bef=,aft=
      cw728(i2)%v=(cw728(i2)%e(0)*vwm(0)-cw728(i2)%e(1)*vwm(1)-c
     & w728(i2)%e(2)*vwm(2)-cw728(i2)%e(3)*vwm(3))
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf314(i3,i1)%e,q=cw728(i2)%e,bef=,aft=
      caux=(cf314(i3,i1)%e(0)*cw728(i2)%e(0)-cf314(i3,i1)%e(1)*c
     & w728(i2)%e(1)-cf314(i3,i1)%e(2)*cw728(i2)%e(2)-cf314(i3,i
     & 1)%e(3)*cw728(i2)%e(3))
      cres_aux(i1,i2,i3,3)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cf314(i3,i1)%e,q=cw56%e,bef=,aft=
      caux=(cf314(i3,i1)%e(0)*cw56%e(0)-cf314(i3,i1)%e(1)*cw56%e
     & (1)-cf314(i3,i1)%e(2)*cw56%e(2)-cf314(i3,i1)%e(3)*cw56%e(
     & 3))
      do i2=1,2
      cres_aux(i1,i2,i3,3)=cres_aux(i1,i2,i3,3)+cw728(i2)%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i2=1,2
* p.q -- p.q=caux,p=cw728(i2)%e,q=cw56%e,bef=,aft=
      caux=(cw728(i2)%e(0)*cw56%e(0)-cw728(i2)%e(1)*cw56%e(1)-cw
     & 728(i2)%e(2)*cw56%e(2)-cw728(i2)%e(3)*cw56%e(3))
      do i1=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,3)=cres_aux(i1,i2,i3,3)+cf314(i3,i1)%v*c
     & aux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,3)=cres(i1,i2,i3,3)+
     &      cres_aux(i1,i2,i3,3)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p324(mu),pwm(mu)=p718(mu),pwp(mu)=p56(mu),e
* fz=cz324(i3,i2),ewm=cw718(i1),ewp=cw56,res=cres_aux(i1?,i2?,i3?,4)
      do mu=0,3
      vfz(mu)=p718(mu)-p56(mu)
      vwm(mu)=p56(mu)-p324(mu)
      vwp(mu)=p324(mu)-p718(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i2=1,2
* p.q -- p.q=cz324(i3,i2)%v,p=cz324(i3,i2)%e,q=vfz,bef=,aft=
      cz324(i3,i2)%v=(cz324(i3,i2)%e(0)*vfz(0)-cz324(i3,i2)%e(1)
     & *vfz(1)-cz324(i3,i2)%e(2)*vfz(2)-cz324(i3,i2)%e(3)*vfz(3)
     & )
      end do
      end do
* vwm%ewm
      do i1=1,2
* p.q -- p.q=cw718(i1)%v,p=cw718(i1)%e,q=vwm,bef=,aft=
      cw718(i1)%v=(cw718(i1)%e(0)*vwm(0)-cw718(i1)%e(1)*vwm(1)-c
     & w718(i1)%e(2)*vwm(2)-cw718(i1)%e(3)*vwm(3))
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i2=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cz324(i3,i2)%e,q=cw718(i1)%e,bef=,aft=
      caux=(cz324(i3,i2)%e(0)*cw718(i1)%e(0)-cz324(i3,i2)%e(1)*c
     & w718(i1)%e(1)-cz324(i3,i2)%e(2)*cw718(i1)%e(2)-cz324(i3,i
     & 2)%e(3)*cw718(i1)%e(3))
      cres_aux(i1,i2,i3,4)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz324(i3,i2)%e,q=cw56%e,bef=,aft=
      caux=(cz324(i3,i2)%e(0)*cw56%e(0)-cz324(i3,i2)%e(1)*cw56%e
     & (1)-cz324(i3,i2)%e(2)*cw56%e(2)-cz324(i3,i2)%e(3)*cw56%e(
     & 3))
      do i1=1,2
      cres_aux(i1,i2,i3,4)=cres_aux(i1,i2,i3,4)+cw718(i1)%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
* p.q -- p.q=caux,p=cw718(i1)%e,q=cw56%e,bef=,aft=
      caux=(cw718(i1)%e(0)*cw56%e(0)-cw718(i1)%e(1)*cw56%e(1)-cw
     & 718(i1)%e(2)*cw56%e(2)-cw718(i1)%e(3)*cw56%e(3))
      do i2=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,4)=cres_aux(i1,i2,i3,4)+cz324(i3,i2)%v*c
     & aux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,4)=cres(i1,i2,i3,4)+
     &      rcotw*cres_aux(i1,i2,i3,4)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p324(mu),pwm(mu)=p718(mu),pwp(mu)=p56(mu),e
* fz=cf324(i3,i2),ewm=cw718(i1),ewp=cw56,res=cres_aux(i1?,i2?,i3?,4)
      do mu=0,3
      vfz(mu)=p718(mu)-p56(mu)
      vwm(mu)=p56(mu)-p324(mu)
      vwp(mu)=p324(mu)-p718(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i2=1,2
* p.q -- p.q=cf324(i3,i2)%v,p=cf324(i3,i2)%e,q=vfz,bef=,aft=
      cf324(i3,i2)%v=(cf324(i3,i2)%e(0)*vfz(0)-cf324(i3,i2)%e(1)
     & *vfz(1)-cf324(i3,i2)%e(2)*vfz(2)-cf324(i3,i2)%e(3)*vfz(3)
     & )
      end do
      end do
* vwm%ewm
      do i1=1,2
* p.q -- p.q=cw718(i1)%v,p=cw718(i1)%e,q=vwm,bef=,aft=
      cw718(i1)%v=(cw718(i1)%e(0)*vwm(0)-cw718(i1)%e(1)*vwm(1)-c
     & w718(i1)%e(2)*vwm(2)-cw718(i1)%e(3)*vwm(3))
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i2=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cf324(i3,i2)%e,q=cw718(i1)%e,bef=,aft=
      caux=(cf324(i3,i2)%e(0)*cw718(i1)%e(0)-cf324(i3,i2)%e(1)*c
     & w718(i1)%e(1)-cf324(i3,i2)%e(2)*cw718(i1)%e(2)-cf324(i3,i
     & 2)%e(3)*cw718(i1)%e(3))
      cres_aux(i1,i2,i3,4)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf324(i3,i2)%e,q=cw56%e,bef=,aft=
      caux=(cf324(i3,i2)%e(0)*cw56%e(0)-cf324(i3,i2)%e(1)*cw56%e
     & (1)-cf324(i3,i2)%e(2)*cw56%e(2)-cf324(i3,i2)%e(3)*cw56%e(
     & 3))
      do i1=1,2
      cres_aux(i1,i2,i3,4)=cres_aux(i1,i2,i3,4)+cw718(i1)%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
* p.q -- p.q=caux,p=cw718(i1)%e,q=cw56%e,bef=,aft=
      caux=(cw718(i1)%e(0)*cw56%e(0)-cw718(i1)%e(1)*cw56%e(1)-cw
     & 718(i1)%e(2)*cw56%e(2)-cw718(i1)%e(3)*cw56%e(3))
      do i2=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,4)=cres_aux(i1,i2,i3,4)+cf324(i3,i2)%v*c
     & aux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,4)=cres(i1,i2,i3,4)+
     &      cres_aux(i1,i2,i3,4)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p728(mu),pwp(mu)=p516(mu),e
* fz=cz34(i3),ewm=cw728(i2),ewp=cw516(i1),res=cres_aux(i1?,i2?,i3?,5)
      do mu=0,3
      vfz(mu)=p728(mu)-p516(mu)
      vwm(mu)=p516(mu)-p34(mu)
      vwp(mu)=p34(mu)-p728(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cz34(i3)%v,p=cz34(i3)%e,q=vfz,bef=,aft=
      cz34(i3)%v=(cz34(i3)%e(0)*vfz(0)-cz34(i3)%e(1)*vfz(1)-cz34
     & (i3)%e(2)*vfz(2)-cz34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
      do i2=1,2
* p.q -- p.q=cw728(i2)%v,p=cw728(i2)%e,q=vwm,bef=,aft=
      cw728(i2)%v=(cw728(i2)%e(0)*vwm(0)-cw728(i2)%e(1)*vwm(1)-c
     & w728(i2)%e(2)*vwm(2)-cw728(i2)%e(3)*vwm(3))
      end do
* vwp%ewp
      do i1=1,2
* p.q -- p.q=cw516(i1)%v,p=cw516(i1)%e,q=vwp,bef=,aft=
      cw516(i1)%v=(cw516(i1)%e(0)*vwp(0)-cw516(i1)%e(1)*vwp(1)-c
     & w516(i1)%e(2)*vwp(2)-cw516(i1)%e(3)*vwp(3))
      end do
* efz%ewm
      do i3=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw728(i2)%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw728(i2)%e(0)-cz34(i3)%e(1)*cw728(i2)
     & %e(1)-cz34(i3)%e(2)*cw728(i2)%e(2)-cz34(i3)%e(3)*cw728(i2
     & )%e(3))
      do i1=1,2
      cres_aux(i1,i2,i3,5)=cw516(i1)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw516(i1)%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw516(i1)%e(0)-cz34(i3)%e(1)*cw516(i1)
     & %e(1)-cz34(i3)%e(2)*cw516(i1)%e(2)-cz34(i3)%e(3)*cw516(i1
     & )%e(3))
      do i2=1,2
      cres_aux(i1,i2,i3,5)=cres_aux(i1,i2,i3,5)+cw728(i2)%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i2=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cw728(i2)%e,q=cw516(i1)%e,bef=,aft=
      caux=(cw728(i2)%e(0)*cw516(i1)%e(0)-cw728(i2)%e(1)*cw516(i
     & 1)%e(1)-cw728(i2)%e(2)*cw516(i1)%e(2)-cw728(i2)%e(3)*cw51
     & 6(i1)%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,5)=cres_aux(i1,i2,i3,5)+cz34(i3)%v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,5)=cres(i1,i2,i3,5)+
     &      rcotw*cres_aux(i1,i2,i3,5)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p728(mu),pwp(mu)=p516(mu),e
* fz=cf34(i3),ewm=cw728(i2),ewp=cw516(i1),res=cres_aux(i1?,i2?,i3?,5)
      do mu=0,3
      vfz(mu)=p728(mu)-p516(mu)
      vwm(mu)=p516(mu)-p34(mu)
      vwp(mu)=p34(mu)-p728(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cf34(i3)%v,p=cf34(i3)%e,q=vfz,bef=,aft=
      cf34(i3)%v=(cf34(i3)%e(0)*vfz(0)-cf34(i3)%e(1)*vfz(1)-cf34
     & (i3)%e(2)*vfz(2)-cf34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
      do i2=1,2
* p.q -- p.q=cw728(i2)%v,p=cw728(i2)%e,q=vwm,bef=,aft=
      cw728(i2)%v=(cw728(i2)%e(0)*vwm(0)-cw728(i2)%e(1)*vwm(1)-c
     & w728(i2)%e(2)*vwm(2)-cw728(i2)%e(3)*vwm(3))
      end do
* vwp%ewp
      do i1=1,2
* p.q -- p.q=cw516(i1)%v,p=cw516(i1)%e,q=vwp,bef=,aft=
      cw516(i1)%v=(cw516(i1)%e(0)*vwp(0)-cw516(i1)%e(1)*vwp(1)-c
     & w516(i1)%e(2)*vwp(2)-cw516(i1)%e(3)*vwp(3))
      end do
* efz%ewm
      do i3=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw728(i2)%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw728(i2)%e(0)-cf34(i3)%e(1)*cw728(i2)
     & %e(1)-cf34(i3)%e(2)*cw728(i2)%e(2)-cf34(i3)%e(3)*cw728(i2
     & )%e(3))
      do i1=1,2
      cres_aux(i1,i2,i3,5)=cw516(i1)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw516(i1)%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw516(i1)%e(0)-cf34(i3)%e(1)*cw516(i1)
     & %e(1)-cf34(i3)%e(2)*cw516(i1)%e(2)-cf34(i3)%e(3)*cw516(i1
     & )%e(3))
      do i2=1,2
      cres_aux(i1,i2,i3,5)=cres_aux(i1,i2,i3,5)+cw728(i2)%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i2=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cw728(i2)%e,q=cw516(i1)%e,bef=,aft=
      caux=(cw728(i2)%e(0)*cw516(i1)%e(0)-cw728(i2)%e(1)*cw516(i
     & 1)%e(1)-cw728(i2)%e(2)*cw516(i1)%e(2)-cw728(i2)%e(3)*cw51
     & 6(i1)%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,5)=cres_aux(i1,i2,i3,5)+cf34(i3)%v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,5)=cres(i1,i2,i3,5)+
     &      cres_aux(i1,i2,i3,5)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p718(mu),pwp(mu)=p526(mu),e
* fz=cz34(i3),ewm=cw718(i1),ewp=cw526(i2),res=cres_aux(i1?,i2?,i3?,6)
      do mu=0,3
      vfz(mu)=p718(mu)-p526(mu)
      vwm(mu)=p526(mu)-p34(mu)
      vwp(mu)=p34(mu)-p718(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cz34(i3)%v,p=cz34(i3)%e,q=vfz,bef=,aft=
      cz34(i3)%v=(cz34(i3)%e(0)*vfz(0)-cz34(i3)%e(1)*vfz(1)-cz34
     & (i3)%e(2)*vfz(2)-cz34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
      do i1=1,2
* p.q -- p.q=cw718(i1)%v,p=cw718(i1)%e,q=vwm,bef=,aft=
      cw718(i1)%v=(cw718(i1)%e(0)*vwm(0)-cw718(i1)%e(1)*vwm(1)-c
     & w718(i1)%e(2)*vwm(2)-cw718(i1)%e(3)*vwm(3))
      end do
* vwp%ewp
      do i2=1,2
* p.q -- p.q=cw526(i2)%v,p=cw526(i2)%e,q=vwp,bef=,aft=
      cw526(i2)%v=(cw526(i2)%e(0)*vwp(0)-cw526(i2)%e(1)*vwp(1)-c
     & w526(i2)%e(2)*vwp(2)-cw526(i2)%e(3)*vwp(3))
      end do
* efz%ewm
      do i3=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw718(i1)%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw718(i1)%e(0)-cz34(i3)%e(1)*cw718(i1)
     & %e(1)-cz34(i3)%e(2)*cw718(i1)%e(2)-cz34(i3)%e(3)*cw718(i1
     & )%e(3))
      do i2=1,2
      cres_aux(i1,i2,i3,6)=cw526(i2)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw526(i2)%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw526(i2)%e(0)-cz34(i3)%e(1)*cw526(i2)
     & %e(1)-cz34(i3)%e(2)*cw526(i2)%e(2)-cz34(i3)%e(3)*cw526(i2
     & )%e(3))
      do i1=1,2
      cres_aux(i1,i2,i3,6)=cres_aux(i1,i2,i3,6)+cw718(i1)%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw718(i1)%e,q=cw526(i2)%e,bef=,aft=
      caux=(cw718(i1)%e(0)*cw526(i2)%e(0)-cw718(i1)%e(1)*cw526(i
     & 2)%e(1)-cw718(i1)%e(2)*cw526(i2)%e(2)-cw718(i1)%e(3)*cw52
     & 6(i2)%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,6)=cres_aux(i1,i2,i3,6)+cz34(i3)%v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,6)=cres(i1,i2,i3,6)+
     &      rcotw*cres_aux(i1,i2,i3,6)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p718(mu),pwp(mu)=p526(mu),e
* fz=cf34(i3),ewm=cw718(i1),ewp=cw526(i2),res=cres_aux(i1?,i2?,i3?,6)
      do mu=0,3
      vfz(mu)=p718(mu)-p526(mu)
      vwm(mu)=p526(mu)-p34(mu)
      vwp(mu)=p34(mu)-p718(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cf34(i3)%v,p=cf34(i3)%e,q=vfz,bef=,aft=
      cf34(i3)%v=(cf34(i3)%e(0)*vfz(0)-cf34(i3)%e(1)*vfz(1)-cf34
     & (i3)%e(2)*vfz(2)-cf34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
      do i1=1,2
* p.q -- p.q=cw718(i1)%v,p=cw718(i1)%e,q=vwm,bef=,aft=
      cw718(i1)%v=(cw718(i1)%e(0)*vwm(0)-cw718(i1)%e(1)*vwm(1)-c
     & w718(i1)%e(2)*vwm(2)-cw718(i1)%e(3)*vwm(3))
      end do
* vwp%ewp
      do i2=1,2
* p.q -- p.q=cw526(i2)%v,p=cw526(i2)%e,q=vwp,bef=,aft=
      cw526(i2)%v=(cw526(i2)%e(0)*vwp(0)-cw526(i2)%e(1)*vwp(1)-c
     & w526(i2)%e(2)*vwp(2)-cw526(i2)%e(3)*vwp(3))
      end do
* efz%ewm
      do i3=1,2
      do i1=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw718(i1)%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw718(i1)%e(0)-cf34(i3)%e(1)*cw718(i1)
     & %e(1)-cf34(i3)%e(2)*cw718(i1)%e(2)-cf34(i3)%e(3)*cw718(i1
     & )%e(3))
      do i2=1,2
      cres_aux(i1,i2,i3,6)=cw526(i2)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw526(i2)%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw526(i2)%e(0)-cf34(i3)%e(1)*cw526(i2)
     & %e(1)-cf34(i3)%e(2)*cw526(i2)%e(2)-cf34(i3)%e(3)*cw526(i2
     & )%e(3))
      do i1=1,2
      cres_aux(i1,i2,i3,6)=cres_aux(i1,i2,i3,6)+cw718(i1)%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw718(i1)%e,q=cw526(i2)%e,bef=,aft=
      caux=(cw718(i1)%e(0)*cw526(i2)%e(0)-cw718(i1)%e(1)*cw526(i
     & 2)%e(1)-cw718(i1)%e(2)*cw526(i2)%e(2)-cw718(i1)%e(3)*cw52
     & 6(i2)%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,6)=cres_aux(i1,i2,i3,6)+cf34(i3)%v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,6)=cres(i1,i2,i3,6)+
     &      cres_aux(i1,i2,i3,6)
       enddo
       enddo
       enddo
  
      endif
  
  
  
  
      if (ilept(id3).ne.1) then
*** triple vertex -- pfz(mu)=p3124(mu),pwm(mu)=p78(mu),pwp(mu)=p56(mu),e
* fz=cz3124(i3,i1,i2),ewm=cw78,ewp=cw56,res=cres_aux(i1?,i2?,i3?,7)
      do mu=0,3
      vfz(mu)=p78(mu)-p56(mu)
      vwm(mu)=p56(mu)-p3124(mu)
      vwp(mu)=p3124(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cz3124(i3,i1,i2)%v,p=cz3124(i3,i1,i2)%e,q=vfz,bef=,aft=
      cz3124(i3,i1,i2)%v=(cz3124(i3,i1,i2)%e(0)*vfz(0)-cz3124(i3
     & ,i1,i2)%e(1)*vfz(1)-cz3124(i3,i1,i2)%e(2)*vfz(2)-cz3124(i
     & 3,i1,i2)%e(3)*vfz(3))
      end do
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz3124(i3,i1,i2)%e,q=cw78%e,bef=,aft=
      caux=(cz3124(i3,i1,i2)%e(0)*cw78%e(0)-cz3124(i3,i1,i2)%e(1
     & )*cw78%e(1)-cz3124(i3,i1,i2)%e(2)*cw78%e(2)-cz3124(i3,i1,
     & i2)%e(3)*cw78%e(3))
      cres_aux(i1,i2,i3,7)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz3124(i3,i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cz3124(i3,i1,i2)%e(0)*cw56%e(0)-cz3124(i3,i1,i2)%e(1
     & )*cw56%e(1)-cz3124(i3,i1,i2)%e(2)*cw56%e(2)-cz3124(i3,i1,
     & i2)%e(3)*cw56%e(3))
      cres_aux(i1,i2,i3,7)=cres_aux(i1,i2,i3,7)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
* p.q -- p.q=caux,p=cw78%e,q=cw56%e,bef=,aft=
      caux=(cw78%e(0)*cw56%e(0)-cw78%e(1)*cw56%e(1)-cw78%e(2)*cw
     & 56%e(2)-cw78%e(3)*cw56%e(3))
      do i1=1,2
      do i2=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,7)=cres_aux(i1,i2,i3,7)+cz3124(i3,i1,i2)
     & %v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,7)=cres(i1,i2,i3,7)+
     &      rcotw*cres_aux(i1,i2,i3,7)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id3).ne.1) then
*** triple vertex -- pfz(mu)=p3124(mu),pwm(mu)=p78(mu),pwp(mu)=p56(mu),e
* fz=cf3124(i3,i1,i2),ewm=cw78,ewp=cw56,res=cres_aux(i1?,i2?,i3?,7)
      do mu=0,3
      vfz(mu)=p78(mu)-p56(mu)
      vwm(mu)=p56(mu)-p3124(mu)
      vwp(mu)=p3124(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cf3124(i3,i1,i2)%v,p=cf3124(i3,i1,i2)%e,q=vfz,bef=,aft=
      cf3124(i3,i1,i2)%v=(cf3124(i3,i1,i2)%e(0)*vfz(0)-cf3124(i3
     & ,i1,i2)%e(1)*vfz(1)-cf3124(i3,i1,i2)%e(2)*vfz(2)-cf3124(i
     & 3,i1,i2)%e(3)*vfz(3))
      end do
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf3124(i3,i1,i2)%e,q=cw78%e,bef=,aft=
      caux=(cf3124(i3,i1,i2)%e(0)*cw78%e(0)-cf3124(i3,i1,i2)%e(1
     & )*cw78%e(1)-cf3124(i3,i1,i2)%e(2)*cw78%e(2)-cf3124(i3,i1,
     & i2)%e(3)*cw78%e(3))
      cres_aux(i1,i2,i3,7)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf3124(i3,i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cf3124(i3,i1,i2)%e(0)*cw56%e(0)-cf3124(i3,i1,i2)%e(1
     & )*cw56%e(1)-cf3124(i3,i1,i2)%e(2)*cw56%e(2)-cf3124(i3,i1,
     & i2)%e(3)*cw56%e(3))
      cres_aux(i1,i2,i3,7)=cres_aux(i1,i2,i3,7)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
* p.q -- p.q=caux,p=cw78%e,q=cw56%e,bef=,aft=
      caux=(cw78%e(0)*cw56%e(0)-cw78%e(1)*cw56%e(1)-cw78%e(2)*cw
     & 56%e(2)-cw78%e(3)*cw56%e(3))
      do i1=1,2
      do i2=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,7)=cres_aux(i1,i2,i3,7)+cf3124(i3,i1,i2)
     & %v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,7)=cres(i1,i2,i3,7)+
     &      cres_aux(i1,i2,i3,7)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id3).ne.1) then
*** triple vertex -- pfz(mu)=p3124(mu),pwm(mu)=p78(mu),pwp(mu)=p56(mu),e
* fz=cz3214(i3,i1,i2),ewm=cw78,ewp=cw56,res=cres_aux(i1?,i2?,i3?,8)
      do mu=0,3
      vfz(mu)=p78(mu)-p56(mu)
      vwm(mu)=p56(mu)-p3124(mu)
      vwp(mu)=p3124(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cz3214(i3,i1,i2)%v,p=cz3214(i3,i1,i2)%e,q=vfz,bef=,aft=
      cz3214(i3,i1,i2)%v=(cz3214(i3,i1,i2)%e(0)*vfz(0)-cz3214(i3
     & ,i1,i2)%e(1)*vfz(1)-cz3214(i3,i1,i2)%e(2)*vfz(2)-cz3214(i
     & 3,i1,i2)%e(3)*vfz(3))
      end do
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz3214(i3,i1,i2)%e,q=cw78%e,bef=,aft=
      caux=(cz3214(i3,i1,i2)%e(0)*cw78%e(0)-cz3214(i3,i1,i2)%e(1
     & )*cw78%e(1)-cz3214(i3,i1,i2)%e(2)*cw78%e(2)-cz3214(i3,i1,
     & i2)%e(3)*cw78%e(3))
      cres_aux(i1,i2,i3,8)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz3214(i3,i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cz3214(i3,i1,i2)%e(0)*cw56%e(0)-cz3214(i3,i1,i2)%e(1
     & )*cw56%e(1)-cz3214(i3,i1,i2)%e(2)*cw56%e(2)-cz3214(i3,i1,
     & i2)%e(3)*cw56%e(3))
      cres_aux(i1,i2,i3,8)=cres_aux(i1,i2,i3,8)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
* p.q -- p.q=caux,p=cw78%e,q=cw56%e,bef=,aft=
      caux=(cw78%e(0)*cw56%e(0)-cw78%e(1)*cw56%e(1)-cw78%e(2)*cw
     & 56%e(2)-cw78%e(3)*cw56%e(3))
      do i1=1,2
      do i2=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,8)=cres_aux(i1,i2,i3,8)+cz3214(i3,i1,i2)
     & %v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,8)=cres(i1,i2,i3,8)+
     &      rcotw*cres_aux(i1,i2,i3,8)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id3).ne.1) then
*** triple vertex -- pfz(mu)=p3124(mu),pwm(mu)=p78(mu),pwp(mu)=p56(mu),e
* fz=cf3214(i3,i1,i2),ewm=cw78,ewp=cw56,res=cres_aux(i1?,i2?,i3?,8)
      do mu=0,3
      vfz(mu)=p78(mu)-p56(mu)
      vwm(mu)=p56(mu)-p3124(mu)
      vwp(mu)=p3124(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cf3214(i3,i1,i2)%v,p=cf3214(i3,i1,i2)%e,q=vfz,bef=,aft=
      cf3214(i3,i1,i2)%v=(cf3214(i3,i1,i2)%e(0)*vfz(0)-cf3214(i3
     & ,i1,i2)%e(1)*vfz(1)-cf3214(i3,i1,i2)%e(2)*vfz(2)-cf3214(i
     & 3,i1,i2)%e(3)*vfz(3))
      end do
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf3214(i3,i1,i2)%e,q=cw78%e,bef=,aft=
      caux=(cf3214(i3,i1,i2)%e(0)*cw78%e(0)-cf3214(i3,i1,i2)%e(1
     & )*cw78%e(1)-cf3214(i3,i1,i2)%e(2)*cw78%e(2)-cf3214(i3,i1,
     & i2)%e(3)*cw78%e(3))
      cres_aux(i1,i2,i3,8)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf3214(i3,i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cf3214(i3,i1,i2)%e(0)*cw56%e(0)-cf3214(i3,i1,i2)%e(1
     & )*cw56%e(1)-cf3214(i3,i1,i2)%e(2)*cw56%e(2)-cf3214(i3,i1,
     & i2)%e(3)*cw56%e(3))
      cres_aux(i1,i2,i3,8)=cres_aux(i1,i2,i3,8)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
* p.q -- p.q=caux,p=cw78%e,q=cw56%e,bef=,aft=
      caux=(cw78%e(0)*cw56%e(0)-cw78%e(1)*cw56%e(1)-cw78%e(2)*cw
     & 56%e(2)-cw78%e(3)*cw56%e(3))
      do i1=1,2
      do i2=1,2
      do i3=1,2
      cres_aux(i1,i2,i3,8)=cres_aux(i1,i2,i3,8)+cf3214(i3,i1,i2)
     & %v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,8)=cres(i1,i2,i3,8)+
     &      cres_aux(i1,i2,i3,8)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id5).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p78(mu),pwp(mu)=p5126(mu),e
* fz=cz34(i3),ewm=cw78,ewp=cw5126(i1,i2),res=cres_aux(i1?,i2?,i3?,9)
      do mu=0,3
      vfz(mu)=p78(mu)-p5126(mu)
      vwm(mu)=p5126(mu)-p34(mu)
      vwp(mu)=p34(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cz34(i3)%v,p=cz34(i3)%e,q=vfz,bef=,aft=
      cz34(i3)%v=(cz34(i3)%e(0)*vfz(0)-cz34(i3)%e(1)*vfz(1)-cz34
     & (i3)%e(2)*vfz(2)-cz34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cw5126(i1,i2)%v,p=cw5126(i1,i2)%e,q=vwp,bef=,aft=
      cw5126(i1,i2)%v=(cw5126(i1,i2)%e(0)*vwp(0)-cw5126(i1,i2)%e
     & (1)*vwp(1)-cw5126(i1,i2)%e(2)*vwp(2)-cw5126(i1,i2)%e(3)*v
     & wp(3))
      end do
      end do
* efz%ewm
      do i3=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw78%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw78%e(0)-cz34(i3)%e(1)*cw78%e(1)-cz34
     & (i3)%e(2)*cw78%e(2)-cz34(i3)%e(3)*cw78%e(3))
      do i1=1,2
      do i2=1,2
      cres_aux(i1,i2,i3,9)=cw5126(i1,i2)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw5126(i1,i2)%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw5126(i1,i2)%e(0)-cz34(i3)%e(1)*cw512
     & 6(i1,i2)%e(1)-cz34(i3)%e(2)*cw5126(i1,i2)%e(2)-cz34(i3)%e
     & (3)*cw5126(i1,i2)%e(3))
      cres_aux(i1,i2,i3,9)=cres_aux(i1,i2,i3,9)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw78%e,q=cw5126(i1,i2)%e,bef=,aft=
      caux=(cw78%e(0)*cw5126(i1,i2)%e(0)-cw78%e(1)*cw5126(i1,i2)
     & %e(1)-cw78%e(2)*cw5126(i1,i2)%e(2)-cw78%e(3)*cw5126(i1,i2
     & )%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,9)=cres_aux(i1,i2,i3,9)+cz34(i3)%v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,9)=cres(i1,i2,i3,9)+
     &      rcotw*cres_aux(i1,i2,i3,9)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id5).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p78(mu),pwp(mu)=p5126(mu),e
* fz=cf34(i3),ewm=cw78,ewp=cw5126(i1,i2),res=cres_aux(i1?,i2?,i3?,9)
      do mu=0,3
      vfz(mu)=p78(mu)-p5126(mu)
      vwm(mu)=p5126(mu)-p34(mu)
      vwp(mu)=p34(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cf34(i3)%v,p=cf34(i3)%e,q=vfz,bef=,aft=
      cf34(i3)%v=(cf34(i3)%e(0)*vfz(0)-cf34(i3)%e(1)*vfz(1)-cf34
     & (i3)%e(2)*vfz(2)-cf34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cw5126(i1,i2)%v,p=cw5126(i1,i2)%e,q=vwp,bef=,aft=
      cw5126(i1,i2)%v=(cw5126(i1,i2)%e(0)*vwp(0)-cw5126(i1,i2)%e
     & (1)*vwp(1)-cw5126(i1,i2)%e(2)*vwp(2)-cw5126(i1,i2)%e(3)*v
     & wp(3))
      end do
      end do
* efz%ewm
      do i3=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw78%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw78%e(0)-cf34(i3)%e(1)*cw78%e(1)-cf34
     & (i3)%e(2)*cw78%e(2)-cf34(i3)%e(3)*cw78%e(3))
      do i1=1,2
      do i2=1,2
      cres_aux(i1,i2,i3,9)=cw5126(i1,i2)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw5126(i1,i2)%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw5126(i1,i2)%e(0)-cf34(i3)%e(1)*cw512
     & 6(i1,i2)%e(1)-cf34(i3)%e(2)*cw5126(i1,i2)%e(2)-cf34(i3)%e
     & (3)*cw5126(i1,i2)%e(3))
      cres_aux(i1,i2,i3,9)=cres_aux(i1,i2,i3,9)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw78%e,q=cw5126(i1,i2)%e,bef=,aft=
      caux=(cw78%e(0)*cw5126(i1,i2)%e(0)-cw78%e(1)*cw5126(i1,i2)
     & %e(1)-cw78%e(2)*cw5126(i1,i2)%e(2)-cw78%e(3)*cw5126(i1,i2
     & )%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,9)=cres_aux(i1,i2,i3,9)+cf34(i3)%v*caux
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,9)=cres(i1,i2,i3,9)+
     &      cres_aux(i1,i2,i3,9)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id5).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p78(mu),pwp(mu)=p5126(mu),e
* fz=cz34(i3),ewm=cw78,ewp=cw5216(i1,i2),res=cres_aux(i1?,i2?,i3?,10)
      do mu=0,3
      vfz(mu)=p78(mu)-p5126(mu)
      vwm(mu)=p5126(mu)-p34(mu)
      vwp(mu)=p34(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cz34(i3)%v,p=cz34(i3)%e,q=vfz,bef=,aft=
      cz34(i3)%v=(cz34(i3)%e(0)*vfz(0)-cz34(i3)%e(1)*vfz(1)-cz34
     & (i3)%e(2)*vfz(2)-cz34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cw5216(i1,i2)%v,p=cw5216(i1,i2)%e,q=vwp,bef=,aft=
      cw5216(i1,i2)%v=(cw5216(i1,i2)%e(0)*vwp(0)-cw5216(i1,i2)%e
     & (1)*vwp(1)-cw5216(i1,i2)%e(2)*vwp(2)-cw5216(i1,i2)%e(3)*v
     & wp(3))
      end do
      end do
* efz%ewm
      do i3=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw78%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw78%e(0)-cz34(i3)%e(1)*cw78%e(1)-cz34
     & (i3)%e(2)*cw78%e(2)-cz34(i3)%e(3)*cw78%e(3))
      do i1=1,2
      do i2=1,2
      cres_aux(i1,i2,i3,10)=cw5216(i1,i2)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw5216(i1,i2)%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw5216(i1,i2)%e(0)-cz34(i3)%e(1)*cw521
     & 6(i1,i2)%e(1)-cz34(i3)%e(2)*cw5216(i1,i2)%e(2)-cz34(i3)%e
     & (3)*cw5216(i1,i2)%e(3))
      cres_aux(i1,i2,i3,10)=cres_aux(i1,i2,i3,10)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw78%e,q=cw5216(i1,i2)%e,bef=,aft=
      caux=(cw78%e(0)*cw5216(i1,i2)%e(0)-cw78%e(1)*cw5216(i1,i2)
     & %e(1)-cw78%e(2)*cw5216(i1,i2)%e(2)-cw78%e(3)*cw5216(i1,i2
     & )%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,10)=cres_aux(i1,i2,i3,10)+cz34(i3)%v*cau
     & x
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,10)=cres(i1,i2,i3,10)+
     &      rcotw*cres_aux(i1,i2,i3,10)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id5).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p78(mu),pwp(mu)=p5126(mu),e
* fz=cf34(i3),ewm=cw78,ewp=cw5216(i1,i2),res=cres_aux(i1?,i2?,i3?,10)
      do mu=0,3
      vfz(mu)=p78(mu)-p5126(mu)
      vwm(mu)=p5126(mu)-p34(mu)
      vwp(mu)=p34(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cf34(i3)%v,p=cf34(i3)%e,q=vfz,bef=,aft=
      cf34(i3)%v=(cf34(i3)%e(0)*vfz(0)-cf34(i3)%e(1)*vfz(1)-cf34
     & (i3)%e(2)*vfz(2)-cf34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cw5216(i1,i2)%v,p=cw5216(i1,i2)%e,q=vwp,bef=,aft=
      cw5216(i1,i2)%v=(cw5216(i1,i2)%e(0)*vwp(0)-cw5216(i1,i2)%e
     & (1)*vwp(1)-cw5216(i1,i2)%e(2)*vwp(2)-cw5216(i1,i2)%e(3)*v
     & wp(3))
      end do
      end do
* efz%ewm
      do i3=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw78%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw78%e(0)-cf34(i3)%e(1)*cw78%e(1)-cf34
     & (i3)%e(2)*cw78%e(2)-cf34(i3)%e(3)*cw78%e(3))
      do i1=1,2
      do i2=1,2
      cres_aux(i1,i2,i3,10)=cw5216(i1,i2)%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw5216(i1,i2)%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw5216(i1,i2)%e(0)-cf34(i3)%e(1)*cw521
     & 6(i1,i2)%e(1)-cf34(i3)%e(2)*cw5216(i1,i2)%e(2)-cf34(i3)%e
     & (3)*cw5216(i1,i2)%e(3))
      cres_aux(i1,i2,i3,10)=cres_aux(i1,i2,i3,10)+cw78%v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw78%e,q=cw5216(i1,i2)%e,bef=,aft=
      caux=(cw78%e(0)*cw5216(i1,i2)%e(0)-cw78%e(1)*cw5216(i1,i2)
     & %e(1)-cw78%e(2)*cw5216(i1,i2)%e(2)-cw78%e(3)*cw5216(i1,i2
     & )%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,10)=cres_aux(i1,i2,i3,10)+cf34(i3)%v*cau
     & x
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,10)=cres(i1,i2,i3,10)+
     &      cres_aux(i1,i2,i3,10)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p7128(mu),pwp(mu)=p56(mu),e
* fz=cz34(i3),ewm=cw7128(i1,i2),ewp=cw56,res=cres_aux(i1?,i2?,i3?,11)
      do mu=0,3
      vfz(mu)=p7128(mu)-p56(mu)
      vwm(mu)=p56(mu)-p34(mu)
      vwp(mu)=p34(mu)-p7128(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cz34(i3)%v,p=cz34(i3)%e,q=vfz,bef=,aft=
      cz34(i3)%v=(cz34(i3)%e(0)*vfz(0)-cz34(i3)%e(1)*vfz(1)-cz34
     & (i3)%e(2)*vfz(2)-cz34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cw7128(i1,i2)%v,p=cw7128(i1,i2)%e,q=vwm,bef=,aft=
      cw7128(i1,i2)%v=(cw7128(i1,i2)%e(0)*vwm(0)-cw7128(i1,i2)%e
     & (1)*vwm(1)-cw7128(i1,i2)%e(2)*vwm(2)-cw7128(i1,i2)%e(3)*v
     & wm(3))
      end do
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw7128(i1,i2)%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw7128(i1,i2)%e(0)-cz34(i3)%e(1)*cw712
     & 8(i1,i2)%e(1)-cz34(i3)%e(2)*cw7128(i1,i2)%e(2)-cz34(i3)%e
     & (3)*cw7128(i1,i2)%e(3))
      cres_aux(i1,i2,i3,11)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw56%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw56%e(0)-cz34(i3)%e(1)*cw56%e(1)-cz34
     & (i3)%e(2)*cw56%e(2)-cz34(i3)%e(3)*cw56%e(3))
      do i1=1,2
      do i2=1,2
      cres_aux(i1,i2,i3,11)=cres_aux(i1,i2,i3,11)+cw7128(i1,i2)%
     & v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw7128(i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cw7128(i1,i2)%e(0)*cw56%e(0)-cw7128(i1,i2)%e(1)*cw56
     & %e(1)-cw7128(i1,i2)%e(2)*cw56%e(2)-cw7128(i1,i2)%e(3)*cw5
     & 6%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,11)=cres_aux(i1,i2,i3,11)+cz34(i3)%v*cau
     & x
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,11)=cres(i1,i2,i3,11)+
     &      rcotw*cres_aux(i1,i2,i3,11)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p7128(mu),pwp(mu)=p56(mu),e
* fz=cf34(i3),ewm=cw7128(i1,i2),ewp=cw56,res=cres_aux(i1?,i2?,i3?,11)
      do mu=0,3
      vfz(mu)=p7128(mu)-p56(mu)
      vwm(mu)=p56(mu)-p34(mu)
      vwp(mu)=p34(mu)-p7128(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cf34(i3)%v,p=cf34(i3)%e,q=vfz,bef=,aft=
      cf34(i3)%v=(cf34(i3)%e(0)*vfz(0)-cf34(i3)%e(1)*vfz(1)-cf34
     & (i3)%e(2)*vfz(2)-cf34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cw7128(i1,i2)%v,p=cw7128(i1,i2)%e,q=vwm,bef=,aft=
      cw7128(i1,i2)%v=(cw7128(i1,i2)%e(0)*vwm(0)-cw7128(i1,i2)%e
     & (1)*vwm(1)-cw7128(i1,i2)%e(2)*vwm(2)-cw7128(i1,i2)%e(3)*v
     & wm(3))
      end do
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw7128(i1,i2)%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw7128(i1,i2)%e(0)-cf34(i3)%e(1)*cw712
     & 8(i1,i2)%e(1)-cf34(i3)%e(2)*cw7128(i1,i2)%e(2)-cf34(i3)%e
     & (3)*cw7128(i1,i2)%e(3))
      cres_aux(i1,i2,i3,11)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw56%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw56%e(0)-cf34(i3)%e(1)*cw56%e(1)-cf34
     & (i3)%e(2)*cw56%e(2)-cf34(i3)%e(3)*cw56%e(3))
      do i1=1,2
      do i2=1,2
      cres_aux(i1,i2,i3,11)=cres_aux(i1,i2,i3,11)+cw7128(i1,i2)%
     & v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw7128(i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cw7128(i1,i2)%e(0)*cw56%e(0)-cw7128(i1,i2)%e(1)*cw56
     & %e(1)-cw7128(i1,i2)%e(2)*cw56%e(2)-cw7128(i1,i2)%e(3)*cw5
     & 6%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,11)=cres_aux(i1,i2,i3,11)+cf34(i3)%v*cau
     & x
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,11)=cres(i1,i2,i3,11)+
     &      cres_aux(i1,i2,i3,11)
       enddo
       enddo
       enddo
  
      endif
  
  
  
      if (ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p7128(mu),pwp(mu)=p56(mu),e
* fz=cz34(i3),ewm=cw7218(i1,i2),ewp=cw56,res=cres_aux(i1?,i2?,i3?,12)
      do mu=0,3
      vfz(mu)=p7128(mu)-p56(mu)
      vwm(mu)=p56(mu)-p34(mu)
      vwp(mu)=p34(mu)-p7128(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cz34(i3)%v,p=cz34(i3)%e,q=vfz,bef=,aft=
      cz34(i3)%v=(cz34(i3)%e(0)*vfz(0)-cz34(i3)%e(1)*vfz(1)-cz34
     & (i3)%e(2)*vfz(2)-cz34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cw7218(i1,i2)%v,p=cw7218(i1,i2)%e,q=vwm,bef=,aft=
      cw7218(i1,i2)%v=(cw7218(i1,i2)%e(0)*vwm(0)-cw7218(i1,i2)%e
     & (1)*vwm(1)-cw7218(i1,i2)%e(2)*vwm(2)-cw7218(i1,i2)%e(3)*v
     & wm(3))
      end do
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw7218(i1,i2)%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw7218(i1,i2)%e(0)-cz34(i3)%e(1)*cw721
     & 8(i1,i2)%e(1)-cz34(i3)%e(2)*cw7218(i1,i2)%e(2)-cz34(i3)%e
     & (3)*cw7218(i1,i2)%e(3))
      cres_aux(i1,i2,i3,12)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
* p.q -- p.q=caux,p=cz34(i3)%e,q=cw56%e,bef=,aft=
      caux=(cz34(i3)%e(0)*cw56%e(0)-cz34(i3)%e(1)*cw56%e(1)-cz34
     & (i3)%e(2)*cw56%e(2)-cz34(i3)%e(3)*cw56%e(3))
      do i1=1,2
      do i2=1,2
      cres_aux(i1,i2,i3,12)=cres_aux(i1,i2,i3,12)+cw7218(i1,i2)%
     & v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw7218(i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cw7218(i1,i2)%e(0)*cw56%e(0)-cw7218(i1,i2)%e(1)*cw56
     & %e(1)-cw7218(i1,i2)%e(2)*cw56%e(2)-cw7218(i1,i2)%e(3)*cw5
     & 6%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,12)=cres_aux(i1,i2,i3,12)+cz34(i3)%v*cau
     & x
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,12)=cres(i1,i2,i3,12)+
     &      rcotw*cres_aux(i1,i2,i3,12)
       enddo
       enddo
       enddo
  
      endif
  
  
      if (ilept(id7).ne.1) then
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p7128(mu),pwp(mu)=p56(mu),e
* fz=cf34(i3),ewm=cw7218(i1,i2),ewp=cw56,res=cres_aux(i1?,i2?,i3?,12)
      do mu=0,3
      vfz(mu)=p7128(mu)-p56(mu)
      vwm(mu)=p56(mu)-p34(mu)
      vwp(mu)=p34(mu)-p7128(mu)
      end do !mu
* vfz%efz
      do i3=1,2
* p.q -- p.q=cf34(i3)%v,p=cf34(i3)%e,q=vfz,bef=,aft=
      cf34(i3)%v=(cf34(i3)%e(0)*vfz(0)-cf34(i3)%e(1)*vfz(1)-cf34
     & (i3)%e(2)*vfz(2)-cf34(i3)%e(3)*vfz(3))
      end do
* vwm%ewm
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cw7218(i1,i2)%v,p=cw7218(i1,i2)%e,q=vwm,bef=,aft=
      cw7218(i1,i2)%v=(cw7218(i1,i2)%e(0)*vwm(0)-cw7218(i1,i2)%e
     & (1)*vwm(1)-cw7218(i1,i2)%e(2)*vwm(2)-cw7218(i1,i2)%e(3)*v
     & wm(3))
      end do
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewm
      do i3=1,2
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw7218(i1,i2)%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw7218(i1,i2)%e(0)-cf34(i3)%e(1)*cw721
     & 8(i1,i2)%e(1)-cf34(i3)%e(2)*cw7218(i1,i2)%e(2)-cf34(i3)%e
     & (3)*cw7218(i1,i2)%e(3))
      cres_aux(i1,i2,i3,12)=cw56%v*caux
      end do
      end do
      end do
* efz%ewp
      do i3=1,2
* p.q -- p.q=caux,p=cf34(i3)%e,q=cw56%e,bef=,aft=
      caux=(cf34(i3)%e(0)*cw56%e(0)-cf34(i3)%e(1)*cw56%e(1)-cf34
     & (i3)%e(2)*cw56%e(2)-cf34(i3)%e(3)*cw56%e(3))
      do i1=1,2
      do i2=1,2
      cres_aux(i1,i2,i3,12)=cres_aux(i1,i2,i3,12)+cw7218(i1,i2)%
     & v*caux
      end do
      end do
      end do
* ewm%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cw7218(i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cw7218(i1,i2)%e(0)*cw56%e(0)-cw7218(i1,i2)%e(1)*cw56
     & %e(1)-cw7218(i1,i2)%e(2)*cw56%e(2)-cw7218(i1,i2)%e(3)*cw5
     & 6%e(3))
      do i3=1,2
      cres_aux(i1,i2,i3,12)=cres_aux(i1,i2,i3,12)+cf34(i3)%v*cau
     & x
      end do
      end do
      end do
  
       do i1=1,2
       do i2=1,2
       do i3=1,2
          cres(i1,i2,i3,12)=cres(i1,i2,i3,12)+
     &      cres_aux(i1,i2,i3,12)
       enddo
       enddo
       enddo
  
      endif
  
  
  
       spk0=sqrt(abs(p3k0*p4k0*p5k0*p6k0*p7k0*p8k0))
  
* final result for color 1-6                                            
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,1)=cres(i1,i2,i3,1)/spk0
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,1)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,2)=cres(i1,i2,i3,2)/spk0
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,2)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,3)=cres(i1,i2,i3,3)/spk0
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,3)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,4)=cres(i1,i2,i3,4)/spk0
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,4)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,5)=cres(i1,i2,i3,5)/spk0
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,5)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id5).ne.1.and.ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,6)=cres(i1,i2,i3,6)/spk0
         enddo
         enddo
         enddo
      else
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,6)=czero
         enddo
         enddo
         enddo
      endif
  
* final result for color 7-12                                           
  
      if (ilept(id3).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,7)=cres(i1,i2,i3,7)/spk0
         enddo
         enddo
         enddo
      else
      	 do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,7)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id3).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,8)=cres(i1,i2,i3,8)/spk0
         enddo
         enddo
         enddo
      else
      	 do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,8)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id5).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,9)=cres(i1,i2,i3,9)/spk0
         enddo
         enddo
         enddo
      else
      	 do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,9)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id5).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,10)=cres(i1,i2,i3,10)/spk0
         enddo
         enddo
         enddo
      else
      	 do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,10)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,11)=cres(i1,i2,i3,11)/spk0
         enddo
         enddo
         enddo
      else
      	 do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,11)=czero
         enddo
         enddo
         enddo
      endif
  
      if (ilept(id7).ne.1) then
         do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,12)=cres(i1,i2,i3,12)/spk0
         enddo
         enddo
         enddo
      else
      	 do i1=1,2
         do i2=1,2
         do i3=1,2
            cres(i1,i2,i3,12)=czero
         enddo
         enddo
         enddo
      endif
  
  
      return
      end
