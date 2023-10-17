************************************************************************
*                                                                       
*  This routine computes the basic amplitude corresponding to 8         
*  outgoing  fermions which can form 4Z                                 
*  All 8 fermions  are supposed  MASSLESS.                              
*  There is NO HIGGS coupling with fermions.                            
*  There are NO neutrinos                                               
*  The input particles are ordered from 1 to 8 in such a way that       
*  odd are particles, even are antiparticles.                           
*  p1, ...p8 are all outgoing momenta                                   
*  id1....id8 give the identities of the outgoing particles             
*                                                                       
* the routine gives the result "cres(2,2,2,2)" which is                 
* the amplitude for the 4 different possible helicities                 
*                                                                       
************************************************************************
  
      subroutine fourz_massless(p1,p2,p3,p4,p5,p6,p7,p8,
     &              id1,id2,id3,id4,id5,id6,id7,id8,cres)
  
      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)
  
* MOMENTA DEFINITION                                                    
  
*single momenta                                                         
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3),
     &  p7(0:3),p8(0:3)
  
*forks momenta                                                          
      dimension p12(0:3),p34(0:3),p56(0:3),p78(0:3)
  
*left insertion momenta                                                 
      dimension p134(0:3),p156(0:3),p178(0:3),
     &          p312(0:3),p356(0:3),p378(0:3),
     &          p512(0:3),p534(0:3),p578(0:3),
     &          p712(0:3),p734(0:3),p756(0:3)
  
*right insertion momenta                                                
      dimension p234(0:3),p256(0:3),p278(0:3),
     &          p412(0:3),p456(0:3),p478(0:3),
     &          p612(0:3),p634(0:3),p678(0:3),
     &          p812(0:3),p834(0:3),p856(0:3)
  
*middle propagator momenta                                              
      dimension p1234(0:3),p1256(0:3),p1278(0:3),p3456(0:3)
*zw                                                                     
      dimension iresz(8)
      dimension pz(0:3)
      dimension ceps1(0:3),ceps2(0:3)
      dimension cscalz(2),cscal1z(2),cscal2z(2)
*zwend                                                                  
  
  
* LINES                                                                 
  
*forks (vector)                                                         
      type pola
         double complex e(0:3),ek0
      end type
      type(pola) cfkz12(2),cfkz34(2),cfkz56(2),cfkz78(2),
     &            cfkf12(2),cfkf34(2),cfkf56(2),cfkf78(2)
  
*left                                                                   
      type l_line
         double complex a(2),c(2)
      end type
      type(l_line) l1_34(2),l1_56(2),l1_78(2),
     &              l3_12(2),l3_56(2),l3_78(2),
     &              l5_12(2),l5_34(2),l5_78(2),
     &              l7_12(2),l7_34(2),l7_56(2)
  
*right                                                                  
      type r_line
         double complex a(2),b(2)
      end type
      type(r_line) r2_34(2),r2_56(2),r2_78(2),
     &              r4_12(2),r4_56(2),r4_78(2),
     &              r6_12(2),r6_34(2),r6_78(2),
     &              r8_12(2),r8_34(2),r8_56(2)
  
*middle                                                                 
      type u_line
         double complex a(2),b(2),c(2),d(2)
      end type
      type(u_line) u134_56(2),u134_78(2),u156_34(2),
     &              u156_78(2),u178_34(2),u178_56(2),
     &              u312_56(2),u312_78(2),u356_12(2),
     &              u356_78(2),u378_12(2),u378_56(2),
     &              u512_34(2),u512_78(2),u534_12(2),
     &              u534_78(2),u578_12(2),u578_34(2),
     &              u712_34(2),u712_56(2),u734_12(2),
     &              u734_56(2),u756_12(2),u756_34(2)
  
*left*middle = left line  (2-forks)                                     
      type(l_line) l1_3456(2,2),l1_3478(2,2),l1_5634(2,2),
     &              l1_5678(2,2),l1_7834(2,2),l1_7856(2,2),
     &              l3_1256(2,2),l3_1278(2,2),l3_5612(2,2),
     &              l3_5678(2,2),l3_7812(2,2),l3_7856(2,2),
     &              l5_1234(2,2),l5_1278(2,2),l5_3412(2,2),
     &              l5_3478(2,2),l5_7812(2,2),l5_7834(2,2),
     &              l7_1234(2,2),l7_1256(2,2),l7_3412(2,2),
     &              l7_3456(2,2),l7_5612(2,2),l7_5634(2,2)
  
*left*midle*right = a (3-forks)                                         
      dimension c3fk12_345678(2,2,2,2),c3fk34_125678(2,2,2,2),
     &          c3fk56_123478(2,2,2,2),c3fk78_123456(2,2,2,2),
     &          c3fk_tot(2,2,2,2)
  
*left vector                                                            
      type(l_line) lz1_234(0:3),lz1_256(0:3),lz1_278(0:3),
     &              lz3_412(0:3),lz3_456(0:3),lz3_478(0:3),
     &              lz5_612(0:3),lz5_634(0:3),lz5_678(0:3),
     &              lz7_812(0:3),lz7_834(0:3),lz7_856(0:3),
     &              lf1_234(0:3),lf1_256(0:3),lf1_278(0:3),
     &              lf3_412(0:3),lf3_456(0:3),lf3_478(0:3),
     &              lf5_612(0:3),lf5_634(0:3),lf5_678(0:3),
     &              lf7_812(0:3),lf7_834(0:3),lf7_856(0:3)
  
*right vector                                                           
      type(r_line) rz2_134(0:3),rz2_156(0:3),rz2_178(0:3),
     &              rz4_312(0:3),rz4_356(0:3),rz4_378(0:3),
     &              rz6_512(0:3),rz6_534(0:3),rz6_578(0:3),
     &              rz8_712(0:3),rz8_734(0:3),rz8_756(0:3),
     &              rf2_134(0:3),rf2_156(0:3),rf2_178(0:3),
     &              rf4_312(0:3),rf4_356(0:3),rf4_378(0:3),
     &              rf6_512(0:3),rf6_534(0:3),rf6_578(0:3),
     &              rf8_712(0:3),rf8_734(0:3),rf8_756(0:3)
  
*left*right vector = a                                                  
      type(pola) cz1234(2,2),cz1256(2,2),cz1278(2,2),
     &           cz3456(2,2),cz3478(2,2),cz5678(2,2),
     &           cf1234(2,2),cf1256(2,2),cf1278(2,2),
     &           cf3456(2,2),cf3478(2,2),cf5678(2,2)
  
* boson connection                                                      
      dimension cbct1234(2,2,2,2),cbct1256(2,2,2,2),
     &      cbct1278(2,2,2,2), caux4(2,2,2,2),cbct_tot(2,2,2,2)
  
* boson fusion                                                          
      dimension ctrip7812(2,2),ctrip3456(2,2),ctrip7834(2,2),
     &          ctrip1256(2,2),ctrip7856(2,2),
     &          ctrip1234(2,2),cbf_tot(2,2,2,2)
  
* amplitude                                                             
      dimension cres(2,2,2,2)
      dimension cunit(2,2,2,2)
  
  
* GENERIC COMMON                                                        
      include 'common.h'
      include 'common_unitarization.h'
*zw                                                                     
      if (i_zz.gt.0) then
          iresz(1)=0
          if (p1(0).gt.0.d0) then
            if (id1.eq.idz(1).or.(i_zz.eq.2.and.id1.eq.idz(3)))
     &           iresz(1)=1
          endif
          iresz(3)=0
          if (p3(0).gt.0.d0) then
            if (id3.eq.idz(1).or.(i_zz.eq.2.and.id3.eq.idz(3)))
     &           iresz(3)=1
          endif
          iresz(5)=0
          if (p5(0).gt.0.d0) then
            if (id5.eq.idz(1).or.(i_zz.eq.2.and.id5.eq.idz(3)))
     &           iresz(5)=1
          endif
          iresz(7)=0
          if (p7(0).gt.0.d0) then
            if (id7.eq.idz(1).or.(i_zz.eq.2.and.id7.eq.idz(3)))
     &           iresz(7)=1
          endif
      endif
*zwend                                                                  
  
* MOMENTA INITIALIZATION                                                
  
*single momenta                                                         
* pk0 -- p=p1
      p1k0=p1(0)-p1(1)
* pk0 -- p=p2
      p2k0=p2(0)-p2(1)
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
  
*forks momenta                                                          
      do m=0,3
        p12(m)=p1(m)+p2(m)
        p34(m)=p3(m)+p4(m)
        p56(m)=p5(m)+p6(m)
        p78(m)=p7(m)+p8(m)
      enddo
  
*left insertion momenta                                                 
      do m=0,3
        p134(m)=p1(m)+p34(m)
      enddo
* pk0 -- p=p134
      p134k0=p134(0)-p134(1)
* p.q -- p.q=p134q,p=p134,q=p134,bef=,aft=
      p134q=(p134(0)*p134(0)-p134(1)*p134(1)-p134(2)*p134(2)-p13
     & 4(3)*p134(3))
      do m=0,3
        p156(m)=p1(m)+p56(m)
      enddo
* pk0 -- p=p156
      p156k0=p156(0)-p156(1)
* p.q -- p.q=p156q,p=p156,q=p156,bef=,aft=
      p156q=(p156(0)*p156(0)-p156(1)*p156(1)-p156(2)*p156(2)-p15
     & 6(3)*p156(3))
      do m=0,3
        p178(m)=p1(m)+p78(m)
      enddo
* pk0 -- p=p178
      p178k0=p178(0)-p178(1)
* p.q -- p.q=p178q,p=p178,q=p178,bef=,aft=
      p178q=(p178(0)*p178(0)-p178(1)*p178(1)-p178(2)*p178(2)-p17
     & 8(3)*p178(3))
      do m=0,3
        p312(m)=p3(m)+p12(m)
      enddo
* pk0 -- p=p312
      p312k0=p312(0)-p312(1)
* p.q -- p.q=p312q,p=p312,q=p312,bef=,aft=
      p312q=(p312(0)*p312(0)-p312(1)*p312(1)-p312(2)*p312(2)-p31
     & 2(3)*p312(3))
      do m=0,3
        p356(m)=p3(m)+p56(m)
      enddo
* pk0 -- p=p356
      p356k0=p356(0)-p356(1)
* p.q -- p.q=p356q,p=p356,q=p356,bef=,aft=
      p356q=(p356(0)*p356(0)-p356(1)*p356(1)-p356(2)*p356(2)-p35
     & 6(3)*p356(3))
      do m=0,3
        p378(m)=p3(m)+p78(m)
      enddo
* pk0 -- p=p378
      p378k0=p378(0)-p378(1)
* p.q -- p.q=p378q,p=p378,q=p378,bef=,aft=
      p378q=(p378(0)*p378(0)-p378(1)*p378(1)-p378(2)*p378(2)-p37
     & 8(3)*p378(3))
      do m=0,3
        p512(m)=p5(m)+p12(m)
      enddo
* pk0 -- p=p512
      p512k0=p512(0)-p512(1)
* p.q -- p.q=p512q,p=p512,q=p512,bef=,aft=
      p512q=(p512(0)*p512(0)-p512(1)*p512(1)-p512(2)*p512(2)-p51
     & 2(3)*p512(3))
      do m=0,3
        p534(m)=p5(m)+p34(m)
      enddo
* pk0 -- p=p534
      p534k0=p534(0)-p534(1)
* p.q -- p.q=p534q,p=p534,q=p534,bef=,aft=
      p534q=(p534(0)*p534(0)-p534(1)*p534(1)-p534(2)*p534(2)-p53
     & 4(3)*p534(3))
      do m=0,3
        p578(m)=p5(m)+p78(m)
      enddo
* pk0 -- p=p578
      p578k0=p578(0)-p578(1)
* p.q -- p.q=p578q,p=p578,q=p578,bef=,aft=
      p578q=(p578(0)*p578(0)-p578(1)*p578(1)-p578(2)*p578(2)-p57
     & 8(3)*p578(3))
      do m=0,3
        p712(m)=p7(m)+p12(m)
      enddo
* pk0 -- p=p712
      p712k0=p712(0)-p712(1)
* p.q -- p.q=p712q,p=p712,q=p712,bef=,aft=
      p712q=(p712(0)*p712(0)-p712(1)*p712(1)-p712(2)*p712(2)-p71
     & 2(3)*p712(3))
      do m=0,3
        p734(m)=p7(m)+p34(m)
      enddo
* pk0 -- p=p734
      p734k0=p734(0)-p734(1)
* p.q -- p.q=p734q,p=p734,q=p734,bef=,aft=
      p734q=(p734(0)*p734(0)-p734(1)*p734(1)-p734(2)*p734(2)-p73
     & 4(3)*p734(3))
      do m=0,3
        p756(m)=p7(m)+p56(m)
      enddo
* pk0 -- p=p756
      p756k0=p756(0)-p756(1)
* p.q -- p.q=p756q,p=p756,q=p756,bef=,aft=
      p756q=(p756(0)*p756(0)-p756(1)*p756(1)-p756(2)*p756(2)-p75
     & 6(3)*p756(3))
  
*right insertion momenta                                                
      do m=0,3
        p234(m)=-p2(m)-p34(m)
      enddo
* pk0 -- p=p234
      p234k0=p234(0)-p234(1)
* p.q -- p.q=p234q,p=p234,q=p234,bef=,aft=
      p234q=(p234(0)*p234(0)-p234(1)*p234(1)-p234(2)*p234(2)-p23
     & 4(3)*p234(3))
      do m=0,3
        p256(m)=-p2(m)-p56(m)
      enddo
* pk0 -- p=p256
      p256k0=p256(0)-p256(1)
* p.q -- p.q=p256q,p=p256,q=p256,bef=,aft=
      p256q=(p256(0)*p256(0)-p256(1)*p256(1)-p256(2)*p256(2)-p25
     & 6(3)*p256(3))
      do m=0,3
        p278(m)=-p2(m)-p78(m)
      enddo
* pk0 -- p=p278
      p278k0=p278(0)-p278(1)
* p.q -- p.q=p278q,p=p278,q=p278,bef=,aft=
      p278q=(p278(0)*p278(0)-p278(1)*p278(1)-p278(2)*p278(2)-p27
     & 8(3)*p278(3))
      do m=0,3
        p412(m)=-p4(m)-p12(m)
      enddo
* pk0 -- p=p412
      p412k0=p412(0)-p412(1)
* p.q -- p.q=p412q,p=p412,q=p412,bef=,aft=
      p412q=(p412(0)*p412(0)-p412(1)*p412(1)-p412(2)*p412(2)-p41
     & 2(3)*p412(3))
      do m=0,3
        p456(m)=-p4(m)-p56(m)
      enddo
* pk0 -- p=p456
      p456k0=p456(0)-p456(1)
* p.q -- p.q=p456q,p=p456,q=p456,bef=,aft=
      p456q=(p456(0)*p456(0)-p456(1)*p456(1)-p456(2)*p456(2)-p45
     & 6(3)*p456(3))
      do m=0,3
        p478(m)=-p4(m)-p78(m)
      enddo
* pk0 -- p=p478
      p478k0=p478(0)-p478(1)
* p.q -- p.q=p478q,p=p478,q=p478,bef=,aft=
      p478q=(p478(0)*p478(0)-p478(1)*p478(1)-p478(2)*p478(2)-p47
     & 8(3)*p478(3))
      do m=0,3
        p612(m)=-p6(m)-p12(m)
      enddo
* pk0 -- p=p612
      p612k0=p612(0)-p612(1)
* p.q -- p.q=p612q,p=p612,q=p612,bef=,aft=
      p612q=(p612(0)*p612(0)-p612(1)*p612(1)-p612(2)*p612(2)-p61
     & 2(3)*p612(3))
      do m=0,3
        p634(m)=-p6(m)-p34(m)
      enddo
* pk0 -- p=p634
      p634k0=p634(0)-p634(1)
* p.q -- p.q=p634q,p=p634,q=p634,bef=,aft=
      p634q=(p634(0)*p634(0)-p634(1)*p634(1)-p634(2)*p634(2)-p63
     & 4(3)*p634(3))
      do m=0,3
        p678(m)=-p6(m)-p78(m)
      enddo
* pk0 -- p=p678
      p678k0=p678(0)-p678(1)
* p.q -- p.q=p678q,p=p678,q=p678,bef=,aft=
      p678q=(p678(0)*p678(0)-p678(1)*p678(1)-p678(2)*p678(2)-p67
     & 8(3)*p678(3))
      do m=0,3
        p812(m)=-p8(m)-p12(m)
      enddo
* pk0 -- p=p812
      p812k0=p812(0)-p812(1)
* p.q -- p.q=p812q,p=p812,q=p812,bef=,aft=
      p812q=(p812(0)*p812(0)-p812(1)*p812(1)-p812(2)*p812(2)-p81
     & 2(3)*p812(3))
      do m=0,3
        p834(m)=-p8(m)-p34(m)
      enddo
* pk0 -- p=p834
      p834k0=p834(0)-p834(1)
* p.q -- p.q=p834q,p=p834,q=p834,bef=,aft=
      p834q=(p834(0)*p834(0)-p834(1)*p834(1)-p834(2)*p834(2)-p83
     & 4(3)*p834(3))
      do m=0,3
        p856(m)=-p8(m)-p56(m)
      enddo
* pk0 -- p=p856
      p856k0=p856(0)-p856(1)
* p.q -- p.q=p856q,p=p856,q=p856,bef=,aft=
      p856q=(p856(0)*p856(0)-p856(1)*p856(1)-p856(2)*p856(2)-p85
     & 6(3)*p856(3))
  
*boson connection propagator momenta                                    
      do m=0,3
        p1234(m)=p12(m)+p34(m)
      enddo
* p.q -- p.q=p1234q,p=p1234,q=p1234,bef=,aft=
      p1234q=(p1234(0)*p1234(0)-p1234(1)*p1234(1)-p1234(2)*p1234
     & (2)-p1234(3)*p1234(3))
      do m=0,3
        p1256(m)=p12(m)+p56(m)
      enddo
* p.q -- p.q=p1256q,p=p1256,q=p1256,bef=,aft=
      p1256q=(p1256(0)*p1256(0)-p1256(1)*p1256(1)-p1256(2)*p1256
     & (2)-p1256(3)*p1256(3))
      do m=0,3
        p1278(m)=p12(m)+p78(m)
      enddo
* p.q -- p.q=p1278q,p=p1278,q=p1278,bef=,aft=
      p1278q=(p1278(0)*p1278(0)-p1278(1)*p1278(1)-p1278(2)*p1278
     & (2)-p1278(3)*p1278(3))
      do m=0,3
        p3456(m)=p34(m)+p56(m)
      enddo
* p.q -- p.q=p3456q,p=p3456,q=p3456,bef=,aft=
      p3456q=(p3456(0)*p3456(0)-p3456(1)*p3456(1)-p3456(2)*p3456
     & (2)-p3456(3)*p3456(3))
  
  
*                                                                       
*  FORKS                                                                
*                                                                       
  
* quqd -- p=p1,q=p2
      quqd=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      s12=2.d0*quqd
  
*zw                                                                     
      cauxmz2=cmz2
      if ((i_zz.eq.2.or.(i_ww.eq.1.and.i_zz.eq.1)).
     &     and.i_osp.eq.2.and.iresz(1).eq.1) then
        cauxmz2=rmz2*cuno-cim*gamz*rmz
      endif
  
  
      ccr=zcr(id1)/(-s12+cauxmz2)
      ccl=zcl(id1)/(-s12+cauxmz2)
* T10 -- qu=p1,qd=p2,v=0,a=cfkz12(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p1(2)*p2(3)+p2(2)*p1(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p1k0*p2(0)+p2k0*p1(0)
      cfkz12(1)%e(0)=ccr*(auxa+ceps_0)
      cfkz12(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=1,a=cfkz12(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p1k0*p2(1)+p2k0*p1(1)
      cfkz12(1)%e(1)=ccr*(auxa+ceps_0)
      cfkz12(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=2,a=cfkz12(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p1k0*p2(3)+p2k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(2)+p2k0*p1(2)
      cfkz12(1)%e(2)=ccr*(auxa+ceps_0)
      cfkz12(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=3,a=cfkz12(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p1k0*p2(2)-p2k0*p1(2)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(3)+p2k0*p1(3)
      cfkz12(1)%e(3)=ccr*(auxa+ceps_0)
      cfkz12(2)%e(3)=ccl*(auxa-ceps_0)
  
  
  
      iceps1=0
      iceps2=0
      if (iresz(1).eq.1) then
*  compute the Z momentum                                               
        do mu=0,3
          pz(mu)=p1(mu)+p2(mu)
        enddo
        if (id1.eq.idz(1)) then
          jpol=ipolz(1)
        elseif (id1.eq.idz(3)) then
          jpol=ipolz(2)
        else
          print *, 'ERROR, MISMATCH IN POLARIZATIONS'
          stop
        endif
  
        if (jpol.ne.0.and.jpol.ne.4) then
          call pol(pz,jpol,ceps1)
          iceps1=1
        endif
        if (jpol.eq.4) then
* pol z transverse                                                      
          ip=2
          call pol(pz,ip,ceps1)
          ip=3
          call pol(pz,ip,ceps2)
          iceps2=1
        endif
  
        if (iceps1.eq.1) then
          do i=1,2
* p.q -- p.q=cscalz(i),p=ceps1,q=cfkz12(i)%e,bef=,aft=
      cscalz(i)=(ceps1(0)*cfkz12(i)%e(0)-ceps1(1)*cfkz12(i)%e(1)
     & -ceps1(2)*cfkz12(i)%e(2)-ceps1(3)*cfkz12(i)%e(3))
  
            do mu=0,3
              cfkz12(i)%e(mu)=cscalz(i)*conjg(ceps1(mu))
            enddo
          enddo
        endif
  
        if (iceps2.eq.1) then
          do i=1,2
* p.q -- p.q=cscal1z(i),p=ceps1,q=cfkz12(i)%e,bef=,aft=
      cscal1z(i)=(ceps1(0)*cfkz12(i)%e(0)-ceps1(1)*cfkz12(i)%e(1
     & )-ceps1(2)*cfkz12(i)%e(2)-ceps1(3)*cfkz12(i)%e(3))
  
* p.q -- p.q=cscal2z(i),p=ceps2,q=cfkz12(i)%e,bef=,aft=
      cscal2z(i)=(ceps2(0)*cfkz12(i)%e(0)-ceps2(1)*cfkz12(i)%e(1
     & )-ceps2(2)*cfkz12(i)%e(2)-ceps2(3)*cfkz12(i)%e(3))
  
            do mu=0,3
              cfkz12(i)%e(mu)=cscal1z(i)*conjg(ceps1(mu))
     &           +cscal2z(i)*conjg(ceps2(mu))
            enddo
          enddo
        endif
      endif
  
*zwend                                                                  
  
      do i1=1,2
* pk0 -- p=cfkz12(i1)%e
      cfkz12(i1)%ek0=cfkz12(i1)%e(0)-cfkz12(i1)%e(1)
      end do
  
  
      if (ineutri(id1).ne.1) then
*zw                                                                     
  
        if(iresz(1).eq.1) then
          do mu=0,3
            cfkf12(1)%e(mu)=czero
            cfkf12(2)%e(mu)=czero
          enddo
        else
  
      ccr=fcr(id1)/(-s12)
      ccl=fcl(id1)/(-s12)
* T10 -- qu=p1,qd=p2,v=0,a=cfkf12(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p1(2)*p2(3)+p2(2)*p1(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p1k0*p2(0)+p2k0*p1(0)
      cfkf12(1)%e(0)=ccr*(auxa+ceps_0)
      cfkf12(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=1,a=cfkf12(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p1k0*p2(1)+p2k0*p1(1)
      cfkf12(1)%e(1)=ccr*(auxa+ceps_0)
      cfkf12(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=2,a=cfkf12(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p1k0*p2(3)+p2k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(2)+p2k0*p1(2)
      cfkf12(1)%e(2)=ccr*(auxa+ceps_0)
      cfkf12(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=3,a=cfkf12(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p1k0*p2(2)-p2k0*p1(2)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(3)+p2k0*p1(3)
      cfkf12(1)%e(3)=ccr*(auxa+ceps_0)
      cfkf12(2)%e(3)=ccl*(auxa-ceps_0)
        endif
*zwend                                                                  
      do i1=1,2
* pk0 -- p=cfkf12(i1)%e
      cfkf12(i1)%ek0=cfkf12(i1)%e(0)-cfkf12(i1)%e(1)
      end do
      endif
  
* quqd -- p=p3,q=p4
      quqd=p3(0)*p4(0)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3)
      s34=2.d0*quqd
  
*zw                                                                     
      cauxmz2=cmz2
      if ((i_zz.eq.2.or.(i_ww.eq.1.and.i_zz.eq.1)).
     &     and.i_osp.eq.2.and.iresz(3).eq.1) then
        cauxmz2=rmz2*cuno-cim*gamz*rmz
      endif
  
  
      ccr=zcr(id3)/(-s34+cauxmz2)
      ccl=zcl(id3)/(-s34+cauxmz2)
* T10 -- qu=p3,qd=p4,v=0,a=cfkz34(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p4(3)+p4(2)*p3(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p3k0*p4(0)+p4k0*p3(0)
      cfkz34(1)%e(0)=ccr*(auxa+ceps_0)
      cfkz34(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=1,a=cfkz34(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p4(1)+p4k0*p3(1)
      cfkz34(1)%e(1)=ccr*(auxa+ceps_0)
      cfkz34(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=2,a=cfkz34(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p4(3)+p4k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(2)+p4k0*p3(2)
      cfkz34(1)%e(2)=ccr*(auxa+ceps_0)
      cfkz34(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=3,a=cfkz34(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p4(2)-p4k0*p3(2)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(3)+p4k0*p3(3)
      cfkz34(1)%e(3)=ccr*(auxa+ceps_0)
      cfkz34(2)%e(3)=ccl*(auxa-ceps_0)
  
  
  
      iceps1=0
      iceps2=0
      if (iresz(3).eq.1) then
*  compute the Z momentum                                               
        do mu=0,3
          pz(mu)=p3(mu)+p4(mu)
        enddo
        if (id3.eq.idz(1)) then
          jpol=ipolz(1)
        elseif (id3.eq.idz(3)) then
          jpol=ipolz(2)
        else
          print *, 'ERROR, MISMATCH IN POLARIZATIONS'
          stop
        endif
  
        if (jpol.ne.0.and.jpol.ne.4) then
          call pol(pz,jpol,ceps1)
          iceps1=1
        endif
        if (jpol.eq.4) then
* pol z transverse                                                      
          ip=2
          call pol(pz,ip,ceps1)
          ip=3
          call pol(pz,ip,ceps2)
          iceps2=1
        endif
  
        if (iceps1.eq.1) then
          do i=1,2
* p.q -- p.q=cscalz(i),p=ceps1,q=cfkz34(i)%e,bef=,aft=
      cscalz(i)=(ceps1(0)*cfkz34(i)%e(0)-ceps1(1)*cfkz34(i)%e(1)
     & -ceps1(2)*cfkz34(i)%e(2)-ceps1(3)*cfkz34(i)%e(3))
  
            do mu=0,3
              cfkz34(i)%e(mu)=cscalz(i)*conjg(ceps1(mu))
            enddo
          enddo
        endif
  
        if (iceps2.eq.1) then
          do i=1,2
* p.q -- p.q=cscal1z(i),p=ceps1,q=cfkz34(i)%e,bef=,aft=
      cscal1z(i)=(ceps1(0)*cfkz34(i)%e(0)-ceps1(1)*cfkz34(i)%e(1
     & )-ceps1(2)*cfkz34(i)%e(2)-ceps1(3)*cfkz34(i)%e(3))
  
* p.q -- p.q=cscal2z(i),p=ceps2,q=cfkz34(i)%e,bef=,aft=
      cscal2z(i)=(ceps2(0)*cfkz34(i)%e(0)-ceps2(1)*cfkz34(i)%e(1
     & )-ceps2(2)*cfkz34(i)%e(2)-ceps2(3)*cfkz34(i)%e(3))
  
            do mu=0,3
              cfkz34(i)%e(mu)=cscal1z(i)*conjg(ceps1(mu))
     &           +cscal2z(i)*conjg(ceps2(mu))
            enddo
          enddo
        endif
      endif
  
*zwend                                                                  
  
      do i1=1,2
* pk0 -- p=cfkz34(i1)%e
      cfkz34(i1)%ek0=cfkz34(i1)%e(0)-cfkz34(i1)%e(1)
      end do
  
  
      if (ineutri(id3).ne.1) then
*zw                                                                     
  
        if(iresz(3).eq.1) then
          do mu=0,3
            cfkf34(1)%e(mu)=czero
            cfkf34(2)%e(mu)=czero
          enddo
        else
  
      ccr=fcr(id3)/(-s34)
      ccl=fcl(id3)/(-s34)
* T10 -- qu=p3,qd=p4,v=0,a=cfkf34(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p4(3)+p4(2)*p3(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p3k0*p4(0)+p4k0*p3(0)
      cfkf34(1)%e(0)=ccr*(auxa+ceps_0)
      cfkf34(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=1,a=cfkf34(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p4(1)+p4k0*p3(1)
      cfkf34(1)%e(1)=ccr*(auxa+ceps_0)
      cfkf34(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=2,a=cfkf34(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p4(3)+p4k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(2)+p4k0*p3(2)
      cfkf34(1)%e(2)=ccr*(auxa+ceps_0)
      cfkf34(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p3,qd=p4,v=3,a=cfkf34(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p4(2)-p4k0*p3(2)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(3)+p4k0*p3(3)
      cfkf34(1)%e(3)=ccr*(auxa+ceps_0)
      cfkf34(2)%e(3)=ccl*(auxa-ceps_0)
        endif
*zwend                                                                  
      do i1=1,2
* pk0 -- p=cfkf34(i1)%e
      cfkf34(i1)%ek0=cfkf34(i1)%e(0)-cfkf34(i1)%e(1)
      end do
      endif
  
* quqd -- p=p5,q=p6
      quqd=p5(0)*p6(0)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3)
      s56=2.d0*quqd
  
*zw                                                                     
      cauxmz2=cmz2
      if ((i_zz.eq.2.or.(i_ww.eq.1.and.i_zz.eq.1)).
     &     and.i_osp.eq.2.and.iresz(5).eq.1) then
        cauxmz2=rmz2*cuno-cim*gamz*rmz
      endif
  
  
      ccr=zcr(id5)/(-s56+cauxmz2)
      ccl=zcl(id5)/(-s56+cauxmz2)
* T10 -- qu=p5,qd=p6,v=0,a=cfkz56(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p5(2)*p6(3)+p6(2)*p5(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p5k0*p6(0)+p6k0*p5(0)
      cfkz56(1)%e(0)=ccr*(auxa+ceps_0)
      cfkz56(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=1,a=cfkz56(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p6(1)+p6k0*p5(1)
      cfkz56(1)%e(1)=ccr*(auxa+ceps_0)
      cfkz56(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=2,a=cfkz56(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p5k0*p6(3)+p6k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(2)+p6k0*p5(2)
      cfkz56(1)%e(2)=ccr*(auxa+ceps_0)
      cfkz56(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=3,a=cfkz56(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p5k0*p6(2)-p6k0*p5(2)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(3)+p6k0*p5(3)
      cfkz56(1)%e(3)=ccr*(auxa+ceps_0)
      cfkz56(2)%e(3)=ccl*(auxa-ceps_0)
  
  
  
      iceps1=0
      iceps2=0
      if (iresz(5).eq.1) then
*  compute the Z momentum                                               
        do mu=0,3
          pz(mu)=p5(mu)+p6(mu)
        enddo
        if (id5.eq.idz(1)) then
          jpol=ipolz(1)
        elseif (id5.eq.idz(3)) then
          jpol=ipolz(2)
        else
          print *, 'ERROR, MISMATCH IN POLARIZATIONS'
          stop
        endif
  
        if (jpol.ne.0.and.jpol.ne.4) then
          call pol(pz,jpol,ceps1)
          iceps1=1
        endif
        if (jpol.eq.4) then
* pol z transverse                                                      
          ip=2
          call pol(pz,ip,ceps1)
          ip=3
          call pol(pz,ip,ceps2)
          iceps2=1
        endif
  
        if (iceps1.eq.1) then
          do i=1,2
* p.q -- p.q=cscalz(i),p=ceps1,q=cfkz56(i)%e,bef=,aft=
      cscalz(i)=(ceps1(0)*cfkz56(i)%e(0)-ceps1(1)*cfkz56(i)%e(1)
     & -ceps1(2)*cfkz56(i)%e(2)-ceps1(3)*cfkz56(i)%e(3))
  
            do mu=0,3
              cfkz56(i)%e(mu)=cscalz(i)*conjg(ceps1(mu))
            enddo
          enddo
        endif
  
        if (iceps2.eq.1) then
          do i=1,2
* p.q -- p.q=cscal1z(i),p=ceps1,q=cfkz56(i)%e,bef=,aft=
      cscal1z(i)=(ceps1(0)*cfkz56(i)%e(0)-ceps1(1)*cfkz56(i)%e(1
     & )-ceps1(2)*cfkz56(i)%e(2)-ceps1(3)*cfkz56(i)%e(3))
  
* p.q -- p.q=cscal2z(i),p=ceps2,q=cfkz56(i)%e,bef=,aft=
      cscal2z(i)=(ceps2(0)*cfkz56(i)%e(0)-ceps2(1)*cfkz56(i)%e(1
     & )-ceps2(2)*cfkz56(i)%e(2)-ceps2(3)*cfkz56(i)%e(3))
  
            do mu=0,3
              cfkz56(i)%e(mu)=cscal1z(i)*conjg(ceps1(mu))
     &           +cscal2z(i)*conjg(ceps2(mu))
            enddo
          enddo
        endif
      endif
  
*zwend                                                                  
  
      do i1=1,2
* pk0 -- p=cfkz56(i1)%e
      cfkz56(i1)%ek0=cfkz56(i1)%e(0)-cfkz56(i1)%e(1)
      end do
  
  
      if (ineutri(id5).ne.1) then
*zw                                                                     
  
        if(iresz(5).eq.1) then
          do mu=0,3
            cfkf56(1)%e(mu)=czero
            cfkf56(2)%e(mu)=czero
          enddo
        else
  
      ccr=fcr(id5)/(-s56)
      ccl=fcl(id5)/(-s56)
* T10 -- qu=p5,qd=p6,v=0,a=cfkf56(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p5(2)*p6(3)+p6(2)*p5(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p5k0*p6(0)+p6k0*p5(0)
      cfkf56(1)%e(0)=ccr*(auxa+ceps_0)
      cfkf56(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=1,a=cfkf56(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p5k0*p6(1)+p6k0*p5(1)
      cfkf56(1)%e(1)=ccr*(auxa+ceps_0)
      cfkf56(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=2,a=cfkf56(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p5k0*p6(3)+p6k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(2)+p6k0*p5(2)
      cfkf56(1)%e(2)=ccr*(auxa+ceps_0)
      cfkf56(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p5,qd=p6,v=3,a=cfkf56(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p5k0*p6(2)-p6k0*p5(2)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(3)+p6k0*p5(3)
      cfkf56(1)%e(3)=ccr*(auxa+ceps_0)
      cfkf56(2)%e(3)=ccl*(auxa-ceps_0)
        endif
*zwend                                                                  
      do i1=1,2
* pk0 -- p=cfkf56(i1)%e
      cfkf56(i1)%ek0=cfkf56(i1)%e(0)-cfkf56(i1)%e(1)
      end do
      endif
  
* quqd -- p=p7,q=p8
      quqd=p7(0)*p8(0)-p7(1)*p8(1)-p7(2)*p8(2)-p7(3)*p8(3)
      s78=2.d0*quqd
  
*zw                                                                     
      cauxmz2=cmz2
      if ((i_zz.eq.2.or.(i_ww.eq.1.and.i_zz.eq.1)).
     &     and.i_osp.eq.2.and.iresz(7).eq.1) then
        cauxmz2=rmz2*cuno-cim*gamz*rmz
      endif
  
  
      ccr=zcr(id7)/(-s78+cauxmz2)
      ccl=zcl(id7)/(-s78+cauxmz2)
* T10 -- qu=p7,qd=p8,v=0,a=cfkz78(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p7(2)*p8(3)+p8(2)*p7(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p7k0*p8(0)+p8k0*p7(0)
      cfkz78(1)%e(0)=ccr*(auxa+ceps_0)
      cfkz78(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=1,a=cfkz78(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p8(1)+p8k0*p7(1)
      cfkz78(1)%e(1)=ccr*(auxa+ceps_0)
      cfkz78(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=2,a=cfkz78(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p7k0*p8(3)+p8k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(2)+p8k0*p7(2)
      cfkz78(1)%e(2)=ccr*(auxa+ceps_0)
      cfkz78(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=3,a=cfkz78(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p7k0*p8(2)-p8k0*p7(2)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(3)+p8k0*p7(3)
      cfkz78(1)%e(3)=ccr*(auxa+ceps_0)
      cfkz78(2)%e(3)=ccl*(auxa-ceps_0)
  
  
  
      iceps1=0
      iceps2=0
      if (iresz(7).eq.1) then
*  compute the Z momentum                                               
        do mu=0,3
          pz(mu)=p7(mu)+p8(mu)
        enddo
        if (id7.eq.idz(1)) then
          jpol=ipolz(1)
        elseif (id7.eq.idz(3)) then
          jpol=ipolz(2)
        else
          print *, 'ERROR, MISMATCH IN POLARIZATIONS'
          stop
        endif
  
        if (jpol.ne.0.and.jpol.ne.4) then
          call pol(pz,jpol,ceps1)
          iceps1=1
        endif
        if (jpol.eq.4) then
* pol z transverse                                                      
          ip=2
          call pol(pz,ip,ceps1)
          ip=3
          call pol(pz,ip,ceps2)
          iceps2=1
        endif
  
        if (iceps1.eq.1) then
          do i=1,2
* p.q -- p.q=cscalz(i),p=ceps1,q=cfkz78(i)%e,bef=,aft=
      cscalz(i)=(ceps1(0)*cfkz78(i)%e(0)-ceps1(1)*cfkz78(i)%e(1)
     & -ceps1(2)*cfkz78(i)%e(2)-ceps1(3)*cfkz78(i)%e(3))
  
            do mu=0,3
              cfkz78(i)%e(mu)=cscalz(i)*conjg(ceps1(mu))
            enddo
          enddo
        endif
  
        if (iceps2.eq.1) then
          do i=1,2
* p.q -- p.q=cscal1z(i),p=ceps1,q=cfkz78(i)%e,bef=,aft=
      cscal1z(i)=(ceps1(0)*cfkz78(i)%e(0)-ceps1(1)*cfkz78(i)%e(1
     & )-ceps1(2)*cfkz78(i)%e(2)-ceps1(3)*cfkz78(i)%e(3))
  
* p.q -- p.q=cscal2z(i),p=ceps2,q=cfkz78(i)%e,bef=,aft=
      cscal2z(i)=(ceps2(0)*cfkz78(i)%e(0)-ceps2(1)*cfkz78(i)%e(1
     & )-ceps2(2)*cfkz78(i)%e(2)-ceps2(3)*cfkz78(i)%e(3))
  
            do mu=0,3
              cfkz78(i)%e(mu)=cscal1z(i)*conjg(ceps1(mu))
     &           +cscal2z(i)*conjg(ceps2(mu))
            enddo
          enddo
        endif
      endif
  
*zwend                                                                  
  
      do i1=1,2
* pk0 -- p=cfkz78(i1)%e
      cfkz78(i1)%ek0=cfkz78(i1)%e(0)-cfkz78(i1)%e(1)
      end do
  
  
      if (ineutri(id7).ne.1) then
*zw                                                                     
  
        if(iresz(7).eq.1) then
          do mu=0,3
            cfkf78(1)%e(mu)=czero
            cfkf78(2)%e(mu)=czero
          enddo
        else
  
      ccr=fcr(id7)/(-s78)
      ccl=fcl(id7)/(-s78)
* T10 -- qu=p7,qd=p8,v=0,a=cfkf78(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p7(2)*p8(3)+p8(2)*p7(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p7k0*p8(0)+p8k0*p7(0)
      cfkf78(1)%e(0)=ccr*(auxa+ceps_0)
      cfkf78(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=1,a=cfkf78(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p7k0*p8(1)+p8k0*p7(1)
      cfkf78(1)%e(1)=ccr*(auxa+ceps_0)
      cfkf78(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=2,a=cfkf78(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p7k0*p8(3)+p8k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(2)+p8k0*p7(2)
      cfkf78(1)%e(2)=ccr*(auxa+ceps_0)
      cfkf78(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p7,qd=p8,v=3,a=cfkf78(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p7k0*p8(2)-p8k0*p7(2)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(3)+p8k0*p7(3)
      cfkf78(1)%e(3)=ccr*(auxa+ceps_0)
      cfkf78(2)%e(3)=ccl*(auxa-ceps_0)
        endif
*zwend                                                                  
      do i1=1,2
* pk0 -- p=cfkf78(i1)%e
      cfkf78(i1)%ek0=cfkf78(i1)%e(0)-cfkf78(i1)%e(1)
      end do
      endif
  
  
*                                                                       
* LEFT INSERTIONS                                                       
*                                                                       
  
  
*zw                                                                     
  
  
        if(iresz(1).eq.1) then
          do i3=1,2
            l1_34(i3)%a(1)=czero
            l1_34(i3)%a(2)=czero
            l1_34(i3)%c(1)=czero
            l1_34(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p1,q=p134
      quqd=p1(0)*p134(0)-p1(1)*p134(1)-p1(2)*p134(2)-p1(3)*p134(
     & 3)
      ccr=zcr(id1)/(p134q*p134k0)
      ccl=zcl(id1)/(p134q*p134k0)
      do i3=1,2
* TL0 -- qu=p1,qd=p134,v=cfkz34(i3)%e,a=l1_34(i3)%a,c=l1_34(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz34(i3)%ek0*(p1(2)*p134(3)-p134(2)*p1(3))+p1k0*
     & (cfkz34(i3)%e(2)*p134(3)-p134(2)*cfkz34(i3)%e(3))-p134k0*
     & (cfkz34(i3)%e(2)*p1(3)-p1(2)*cfkz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz34(i3)%e(3)*p1k0+p1(3)*cfkz34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz34(i3)%e(0)*p1(0)-cfkz34(i3)%e(1)*p1(1)-cfkz34(i3
     & )%e(2)*p1(2)-cfkz34(i3)%e(3)*p1(3)
      cvqd=cfkz34(i3)%e(0)*p134(0)-cfkz34(i3)%e(1)*p134(1)-cfkz3
     & 4(i3)%e(2)*p134(2)-cfkz34(i3)%e(3)*p134(3)
      cauxa=-cfkz34(i3)%ek0*quqd+p1k0*cvqd+p134k0*cvqu
      cauxc=+cfkz34(i3)%ek0*p1(2)-p1k0*cfkz34(i3)%e(2)
      l1_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      l1_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      l1_34(i3)%c(1)=ccr*(cauxc+ceps_1)
      l1_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id1).ne.1.and.ineutri(id3).ne.1) then
      ccr=fcr(id1)/(p134q*p134k0)
      ccl=fcl(id1)/(p134q*p134k0)
      do i3=1,2
* TL0 -- qu=p1,qd=p134,v=cfkf34(i3)%e,a=l1_34(i3)%a,c=l1_34(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf34(i3)%ek0*(p1(2)*p134(3)-p134(2)*p1(3))+p1k0*
     & (cfkf34(i3)%e(2)*p134(3)-p134(2)*cfkf34(i3)%e(3))-p134k0*
     & (cfkf34(i3)%e(2)*p1(3)-p1(2)*cfkf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf34(i3)%e(3)*p1k0+p1(3)*cfkf34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf34(i3)%e(0)*p1(0)-cfkf34(i3)%e(1)*p1(1)-cfkf34(i3
     & )%e(2)*p1(2)-cfkf34(i3)%e(3)*p1(3)
      cvqd=cfkf34(i3)%e(0)*p134(0)-cfkf34(i3)%e(1)*p134(1)-cfkf3
     & 4(i3)%e(2)*p134(2)-cfkf34(i3)%e(3)*p134(3)
      cauxa=-cfkf34(i3)%ek0*quqd+p1k0*cvqd+p134k0*cvqu
      cauxc=+cfkf34(i3)%ek0*p1(2)-p1k0*cfkf34(i3)%e(2)
      l1_34(i3)%a(1)=l1_34(i3)%a(1)+ccr*(cauxa+ceps_0)
      l1_34(i3)%a(2)=l1_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      l1_34(i3)%c(1)=l1_34(i3)%c(1)+ccr*(cauxc+ceps_1)
      l1_34(i3)%c(2)=l1_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(1).eq.1) then
          do i3=1,2
            l1_56(i3)%a(1)=czero
            l1_56(i3)%a(2)=czero
            l1_56(i3)%c(1)=czero
            l1_56(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p1,q=p156
      quqd=p1(0)*p156(0)-p1(1)*p156(1)-p1(2)*p156(2)-p1(3)*p156(
     & 3)
      ccr=zcr(id1)/(p156q*p156k0)
      ccl=zcl(id1)/(p156q*p156k0)
      do i3=1,2
* TL0 -- qu=p1,qd=p156,v=cfkz56(i3)%e,a=l1_56(i3)%a,c=l1_56(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz56(i3)%ek0*(p1(2)*p156(3)-p156(2)*p1(3))+p1k0*
     & (cfkz56(i3)%e(2)*p156(3)-p156(2)*cfkz56(i3)%e(3))-p156k0*
     & (cfkz56(i3)%e(2)*p1(3)-p1(2)*cfkz56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz56(i3)%e(3)*p1k0+p1(3)*cfkz56(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz56(i3)%e(0)*p1(0)-cfkz56(i3)%e(1)*p1(1)-cfkz56(i3
     & )%e(2)*p1(2)-cfkz56(i3)%e(3)*p1(3)
      cvqd=cfkz56(i3)%e(0)*p156(0)-cfkz56(i3)%e(1)*p156(1)-cfkz5
     & 6(i3)%e(2)*p156(2)-cfkz56(i3)%e(3)*p156(3)
      cauxa=-cfkz56(i3)%ek0*quqd+p1k0*cvqd+p156k0*cvqu
      cauxc=+cfkz56(i3)%ek0*p1(2)-p1k0*cfkz56(i3)%e(2)
      l1_56(i3)%a(1)=ccr*(cauxa+ceps_0)
      l1_56(i3)%a(2)=ccl*(cauxa-ceps_0)
      l1_56(i3)%c(1)=ccr*(cauxc+ceps_1)
      l1_56(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id1).ne.1.and.ineutri(id5).ne.1) then
      ccr=fcr(id1)/(p156q*p156k0)
      ccl=fcl(id1)/(p156q*p156k0)
      do i3=1,2
* TL0 -- qu=p1,qd=p156,v=cfkf56(i3)%e,a=l1_56(i3)%a,c=l1_56(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf56(i3)%ek0*(p1(2)*p156(3)-p156(2)*p1(3))+p1k0*
     & (cfkf56(i3)%e(2)*p156(3)-p156(2)*cfkf56(i3)%e(3))-p156k0*
     & (cfkf56(i3)%e(2)*p1(3)-p1(2)*cfkf56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf56(i3)%e(3)*p1k0+p1(3)*cfkf56(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf56(i3)%e(0)*p1(0)-cfkf56(i3)%e(1)*p1(1)-cfkf56(i3
     & )%e(2)*p1(2)-cfkf56(i3)%e(3)*p1(3)
      cvqd=cfkf56(i3)%e(0)*p156(0)-cfkf56(i3)%e(1)*p156(1)-cfkf5
     & 6(i3)%e(2)*p156(2)-cfkf56(i3)%e(3)*p156(3)
      cauxa=-cfkf56(i3)%ek0*quqd+p1k0*cvqd+p156k0*cvqu
      cauxc=+cfkf56(i3)%ek0*p1(2)-p1k0*cfkf56(i3)%e(2)
      l1_56(i3)%a(1)=l1_56(i3)%a(1)+ccr*(cauxa+ceps_0)
      l1_56(i3)%a(2)=l1_56(i3)%a(2)+ccl*(cauxa-ceps_0)
      l1_56(i3)%c(1)=l1_56(i3)%c(1)+ccr*(cauxc+ceps_1)
      l1_56(i3)%c(2)=l1_56(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(1).eq.1) then
          do i3=1,2
            l1_78(i3)%a(1)=czero
            l1_78(i3)%a(2)=czero
            l1_78(i3)%c(1)=czero
            l1_78(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p1,q=p178
      quqd=p1(0)*p178(0)-p1(1)*p178(1)-p1(2)*p178(2)-p1(3)*p178(
     & 3)
      ccr=zcr(id1)/(p178q*p178k0)
      ccl=zcl(id1)/(p178q*p178k0)
      do i3=1,2
* TL0 -- qu=p1,qd=p178,v=cfkz78(i3)%e,a=l1_78(i3)%a,c=l1_78(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz78(i3)%ek0*(p1(2)*p178(3)-p178(2)*p1(3))+p1k0*
     & (cfkz78(i3)%e(2)*p178(3)-p178(2)*cfkz78(i3)%e(3))-p178k0*
     & (cfkz78(i3)%e(2)*p1(3)-p1(2)*cfkz78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz78(i3)%e(3)*p1k0+p1(3)*cfkz78(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz78(i3)%e(0)*p1(0)-cfkz78(i3)%e(1)*p1(1)-cfkz78(i3
     & )%e(2)*p1(2)-cfkz78(i3)%e(3)*p1(3)
      cvqd=cfkz78(i3)%e(0)*p178(0)-cfkz78(i3)%e(1)*p178(1)-cfkz7
     & 8(i3)%e(2)*p178(2)-cfkz78(i3)%e(3)*p178(3)
      cauxa=-cfkz78(i3)%ek0*quqd+p1k0*cvqd+p178k0*cvqu
      cauxc=+cfkz78(i3)%ek0*p1(2)-p1k0*cfkz78(i3)%e(2)
      l1_78(i3)%a(1)=ccr*(cauxa+ceps_0)
      l1_78(i3)%a(2)=ccl*(cauxa-ceps_0)
      l1_78(i3)%c(1)=ccr*(cauxc+ceps_1)
      l1_78(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id1).ne.1.and.ineutri(id7).ne.1) then
      ccr=fcr(id1)/(p178q*p178k0)
      ccl=fcl(id1)/(p178q*p178k0)
      do i3=1,2
* TL0 -- qu=p1,qd=p178,v=cfkf78(i3)%e,a=l1_78(i3)%a,c=l1_78(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf78(i3)%ek0*(p1(2)*p178(3)-p178(2)*p1(3))+p1k0*
     & (cfkf78(i3)%e(2)*p178(3)-p178(2)*cfkf78(i3)%e(3))-p178k0*
     & (cfkf78(i3)%e(2)*p1(3)-p1(2)*cfkf78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf78(i3)%e(3)*p1k0+p1(3)*cfkf78(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf78(i3)%e(0)*p1(0)-cfkf78(i3)%e(1)*p1(1)-cfkf78(i3
     & )%e(2)*p1(2)-cfkf78(i3)%e(3)*p1(3)
      cvqd=cfkf78(i3)%e(0)*p178(0)-cfkf78(i3)%e(1)*p178(1)-cfkf7
     & 8(i3)%e(2)*p178(2)-cfkf78(i3)%e(3)*p178(3)
      cauxa=-cfkf78(i3)%ek0*quqd+p1k0*cvqd+p178k0*cvqu
      cauxc=+cfkf78(i3)%ek0*p1(2)-p1k0*cfkf78(i3)%e(2)
      l1_78(i3)%a(1)=l1_78(i3)%a(1)+ccr*(cauxa+ceps_0)
      l1_78(i3)%a(2)=l1_78(i3)%a(2)+ccl*(cauxa-ceps_0)
      l1_78(i3)%c(1)=l1_78(i3)%c(1)+ccr*(cauxc+ceps_1)
      l1_78(i3)%c(2)=l1_78(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(3).eq.1) then
          do i3=1,2
            l3_12(i3)%a(1)=czero
            l3_12(i3)%a(2)=czero
            l3_12(i3)%c(1)=czero
            l3_12(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p3,q=p312
      quqd=p3(0)*p312(0)-p3(1)*p312(1)-p3(2)*p312(2)-p3(3)*p312(
     & 3)
      ccr=zcr(id3)/(p312q*p312k0)
      ccl=zcl(id3)/(p312q*p312k0)
      do i3=1,2
* TL0 -- qu=p3,qd=p312,v=cfkz12(i3)%e,a=l3_12(i3)%a,c=l3_12(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz12(i3)%ek0*(p3(2)*p312(3)-p312(2)*p3(3))+p3k0*
     & (cfkz12(i3)%e(2)*p312(3)-p312(2)*cfkz12(i3)%e(3))-p312k0*
     & (cfkz12(i3)%e(2)*p3(3)-p3(2)*cfkz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz12(i3)%e(3)*p3k0+p3(3)*cfkz12(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz12(i3)%e(0)*p3(0)-cfkz12(i3)%e(1)*p3(1)-cfkz12(i3
     & )%e(2)*p3(2)-cfkz12(i3)%e(3)*p3(3)
      cvqd=cfkz12(i3)%e(0)*p312(0)-cfkz12(i3)%e(1)*p312(1)-cfkz1
     & 2(i3)%e(2)*p312(2)-cfkz12(i3)%e(3)*p312(3)
      cauxa=-cfkz12(i3)%ek0*quqd+p3k0*cvqd+p312k0*cvqu
      cauxc=+cfkz12(i3)%ek0*p3(2)-p3k0*cfkz12(i3)%e(2)
      l3_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      l3_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      l3_12(i3)%c(1)=ccr*(cauxc+ceps_1)
      l3_12(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id3).ne.1.and.ineutri(id1).ne.1) then
      ccr=fcr(id3)/(p312q*p312k0)
      ccl=fcl(id3)/(p312q*p312k0)
      do i3=1,2
* TL0 -- qu=p3,qd=p312,v=cfkf12(i3)%e,a=l3_12(i3)%a,c=l3_12(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf12(i3)%ek0*(p3(2)*p312(3)-p312(2)*p3(3))+p3k0*
     & (cfkf12(i3)%e(2)*p312(3)-p312(2)*cfkf12(i3)%e(3))-p312k0*
     & (cfkf12(i3)%e(2)*p3(3)-p3(2)*cfkf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf12(i3)%e(3)*p3k0+p3(3)*cfkf12(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf12(i3)%e(0)*p3(0)-cfkf12(i3)%e(1)*p3(1)-cfkf12(i3
     & )%e(2)*p3(2)-cfkf12(i3)%e(3)*p3(3)
      cvqd=cfkf12(i3)%e(0)*p312(0)-cfkf12(i3)%e(1)*p312(1)-cfkf1
     & 2(i3)%e(2)*p312(2)-cfkf12(i3)%e(3)*p312(3)
      cauxa=-cfkf12(i3)%ek0*quqd+p3k0*cvqd+p312k0*cvqu
      cauxc=+cfkf12(i3)%ek0*p3(2)-p3k0*cfkf12(i3)%e(2)
      l3_12(i3)%a(1)=l3_12(i3)%a(1)+ccr*(cauxa+ceps_0)
      l3_12(i3)%a(2)=l3_12(i3)%a(2)+ccl*(cauxa-ceps_0)
      l3_12(i3)%c(1)=l3_12(i3)%c(1)+ccr*(cauxc+ceps_1)
      l3_12(i3)%c(2)=l3_12(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(3).eq.1) then
          do i3=1,2
            l3_56(i3)%a(1)=czero
            l3_56(i3)%a(2)=czero
            l3_56(i3)%c(1)=czero
            l3_56(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p3,q=p356
      quqd=p3(0)*p356(0)-p3(1)*p356(1)-p3(2)*p356(2)-p3(3)*p356(
     & 3)
      ccr=zcr(id3)/(p356q*p356k0)
      ccl=zcl(id3)/(p356q*p356k0)
      do i3=1,2
* TL0 -- qu=p3,qd=p356,v=cfkz56(i3)%e,a=l3_56(i3)%a,c=l3_56(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz56(i3)%ek0*(p3(2)*p356(3)-p356(2)*p3(3))+p3k0*
     & (cfkz56(i3)%e(2)*p356(3)-p356(2)*cfkz56(i3)%e(3))-p356k0*
     & (cfkz56(i3)%e(2)*p3(3)-p3(2)*cfkz56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz56(i3)%e(3)*p3k0+p3(3)*cfkz56(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz56(i3)%e(0)*p3(0)-cfkz56(i3)%e(1)*p3(1)-cfkz56(i3
     & )%e(2)*p3(2)-cfkz56(i3)%e(3)*p3(3)
      cvqd=cfkz56(i3)%e(0)*p356(0)-cfkz56(i3)%e(1)*p356(1)-cfkz5
     & 6(i3)%e(2)*p356(2)-cfkz56(i3)%e(3)*p356(3)
      cauxa=-cfkz56(i3)%ek0*quqd+p3k0*cvqd+p356k0*cvqu
      cauxc=+cfkz56(i3)%ek0*p3(2)-p3k0*cfkz56(i3)%e(2)
      l3_56(i3)%a(1)=ccr*(cauxa+ceps_0)
      l3_56(i3)%a(2)=ccl*(cauxa-ceps_0)
      l3_56(i3)%c(1)=ccr*(cauxc+ceps_1)
      l3_56(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id3).ne.1.and.ineutri(id5).ne.1) then
      ccr=fcr(id3)/(p356q*p356k0)
      ccl=fcl(id3)/(p356q*p356k0)
      do i3=1,2
* TL0 -- qu=p3,qd=p356,v=cfkf56(i3)%e,a=l3_56(i3)%a,c=l3_56(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf56(i3)%ek0*(p3(2)*p356(3)-p356(2)*p3(3))+p3k0*
     & (cfkf56(i3)%e(2)*p356(3)-p356(2)*cfkf56(i3)%e(3))-p356k0*
     & (cfkf56(i3)%e(2)*p3(3)-p3(2)*cfkf56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf56(i3)%e(3)*p3k0+p3(3)*cfkf56(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf56(i3)%e(0)*p3(0)-cfkf56(i3)%e(1)*p3(1)-cfkf56(i3
     & )%e(2)*p3(2)-cfkf56(i3)%e(3)*p3(3)
      cvqd=cfkf56(i3)%e(0)*p356(0)-cfkf56(i3)%e(1)*p356(1)-cfkf5
     & 6(i3)%e(2)*p356(2)-cfkf56(i3)%e(3)*p356(3)
      cauxa=-cfkf56(i3)%ek0*quqd+p3k0*cvqd+p356k0*cvqu
      cauxc=+cfkf56(i3)%ek0*p3(2)-p3k0*cfkf56(i3)%e(2)
      l3_56(i3)%a(1)=l3_56(i3)%a(1)+ccr*(cauxa+ceps_0)
      l3_56(i3)%a(2)=l3_56(i3)%a(2)+ccl*(cauxa-ceps_0)
      l3_56(i3)%c(1)=l3_56(i3)%c(1)+ccr*(cauxc+ceps_1)
      l3_56(i3)%c(2)=l3_56(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(3).eq.1) then
          do i3=1,2
            l3_78(i3)%a(1)=czero
            l3_78(i3)%a(2)=czero
            l3_78(i3)%c(1)=czero
            l3_78(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p3,q=p378
      quqd=p3(0)*p378(0)-p3(1)*p378(1)-p3(2)*p378(2)-p3(3)*p378(
     & 3)
      ccr=zcr(id3)/(p378q*p378k0)
      ccl=zcl(id3)/(p378q*p378k0)
      do i3=1,2
* TL0 -- qu=p3,qd=p378,v=cfkz78(i3)%e,a=l3_78(i3)%a,c=l3_78(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz78(i3)%ek0*(p3(2)*p378(3)-p378(2)*p3(3))+p3k0*
     & (cfkz78(i3)%e(2)*p378(3)-p378(2)*cfkz78(i3)%e(3))-p378k0*
     & (cfkz78(i3)%e(2)*p3(3)-p3(2)*cfkz78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz78(i3)%e(3)*p3k0+p3(3)*cfkz78(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz78(i3)%e(0)*p3(0)-cfkz78(i3)%e(1)*p3(1)-cfkz78(i3
     & )%e(2)*p3(2)-cfkz78(i3)%e(3)*p3(3)
      cvqd=cfkz78(i3)%e(0)*p378(0)-cfkz78(i3)%e(1)*p378(1)-cfkz7
     & 8(i3)%e(2)*p378(2)-cfkz78(i3)%e(3)*p378(3)
      cauxa=-cfkz78(i3)%ek0*quqd+p3k0*cvqd+p378k0*cvqu
      cauxc=+cfkz78(i3)%ek0*p3(2)-p3k0*cfkz78(i3)%e(2)
      l3_78(i3)%a(1)=ccr*(cauxa+ceps_0)
      l3_78(i3)%a(2)=ccl*(cauxa-ceps_0)
      l3_78(i3)%c(1)=ccr*(cauxc+ceps_1)
      l3_78(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id3).ne.1.and.ineutri(id7).ne.1) then
      ccr=fcr(id3)/(p378q*p378k0)
      ccl=fcl(id3)/(p378q*p378k0)
      do i3=1,2
* TL0 -- qu=p3,qd=p378,v=cfkf78(i3)%e,a=l3_78(i3)%a,c=l3_78(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf78(i3)%ek0*(p3(2)*p378(3)-p378(2)*p3(3))+p3k0*
     & (cfkf78(i3)%e(2)*p378(3)-p378(2)*cfkf78(i3)%e(3))-p378k0*
     & (cfkf78(i3)%e(2)*p3(3)-p3(2)*cfkf78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf78(i3)%e(3)*p3k0+p3(3)*cfkf78(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf78(i3)%e(0)*p3(0)-cfkf78(i3)%e(1)*p3(1)-cfkf78(i3
     & )%e(2)*p3(2)-cfkf78(i3)%e(3)*p3(3)
      cvqd=cfkf78(i3)%e(0)*p378(0)-cfkf78(i3)%e(1)*p378(1)-cfkf7
     & 8(i3)%e(2)*p378(2)-cfkf78(i3)%e(3)*p378(3)
      cauxa=-cfkf78(i3)%ek0*quqd+p3k0*cvqd+p378k0*cvqu
      cauxc=+cfkf78(i3)%ek0*p3(2)-p3k0*cfkf78(i3)%e(2)
      l3_78(i3)%a(1)=l3_78(i3)%a(1)+ccr*(cauxa+ceps_0)
      l3_78(i3)%a(2)=l3_78(i3)%a(2)+ccl*(cauxa-ceps_0)
      l3_78(i3)%c(1)=l3_78(i3)%c(1)+ccr*(cauxc+ceps_1)
      l3_78(i3)%c(2)=l3_78(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(5).eq.1) then
          do i3=1,2
            l5_12(i3)%a(1)=czero
            l5_12(i3)%a(2)=czero
            l5_12(i3)%c(1)=czero
            l5_12(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p5,q=p512
      quqd=p5(0)*p512(0)-p5(1)*p512(1)-p5(2)*p512(2)-p5(3)*p512(
     & 3)
      ccr=zcr(id5)/(p512q*p512k0)
      ccl=zcl(id5)/(p512q*p512k0)
      do i3=1,2
* TL0 -- qu=p5,qd=p512,v=cfkz12(i3)%e,a=l5_12(i3)%a,c=l5_12(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz12(i3)%ek0*(p5(2)*p512(3)-p512(2)*p5(3))+p5k0*
     & (cfkz12(i3)%e(2)*p512(3)-p512(2)*cfkz12(i3)%e(3))-p512k0*
     & (cfkz12(i3)%e(2)*p5(3)-p5(2)*cfkz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz12(i3)%e(3)*p5k0+p5(3)*cfkz12(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz12(i3)%e(0)*p5(0)-cfkz12(i3)%e(1)*p5(1)-cfkz12(i3
     & )%e(2)*p5(2)-cfkz12(i3)%e(3)*p5(3)
      cvqd=cfkz12(i3)%e(0)*p512(0)-cfkz12(i3)%e(1)*p512(1)-cfkz1
     & 2(i3)%e(2)*p512(2)-cfkz12(i3)%e(3)*p512(3)
      cauxa=-cfkz12(i3)%ek0*quqd+p5k0*cvqd+p512k0*cvqu
      cauxc=+cfkz12(i3)%ek0*p5(2)-p5k0*cfkz12(i3)%e(2)
      l5_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      l5_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      l5_12(i3)%c(1)=ccr*(cauxc+ceps_1)
      l5_12(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id5).ne.1.and.ineutri(id1).ne.1) then
      ccr=fcr(id5)/(p512q*p512k0)
      ccl=fcl(id5)/(p512q*p512k0)
      do i3=1,2
* TL0 -- qu=p5,qd=p512,v=cfkf12(i3)%e,a=l5_12(i3)%a,c=l5_12(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf12(i3)%ek0*(p5(2)*p512(3)-p512(2)*p5(3))+p5k0*
     & (cfkf12(i3)%e(2)*p512(3)-p512(2)*cfkf12(i3)%e(3))-p512k0*
     & (cfkf12(i3)%e(2)*p5(3)-p5(2)*cfkf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf12(i3)%e(3)*p5k0+p5(3)*cfkf12(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf12(i3)%e(0)*p5(0)-cfkf12(i3)%e(1)*p5(1)-cfkf12(i3
     & )%e(2)*p5(2)-cfkf12(i3)%e(3)*p5(3)
      cvqd=cfkf12(i3)%e(0)*p512(0)-cfkf12(i3)%e(1)*p512(1)-cfkf1
     & 2(i3)%e(2)*p512(2)-cfkf12(i3)%e(3)*p512(3)
      cauxa=-cfkf12(i3)%ek0*quqd+p5k0*cvqd+p512k0*cvqu
      cauxc=+cfkf12(i3)%ek0*p5(2)-p5k0*cfkf12(i3)%e(2)
      l5_12(i3)%a(1)=l5_12(i3)%a(1)+ccr*(cauxa+ceps_0)
      l5_12(i3)%a(2)=l5_12(i3)%a(2)+ccl*(cauxa-ceps_0)
      l5_12(i3)%c(1)=l5_12(i3)%c(1)+ccr*(cauxc+ceps_1)
      l5_12(i3)%c(2)=l5_12(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(5).eq.1) then
          do i3=1,2
            l5_34(i3)%a(1)=czero
            l5_34(i3)%a(2)=czero
            l5_34(i3)%c(1)=czero
            l5_34(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p5,q=p534
      quqd=p5(0)*p534(0)-p5(1)*p534(1)-p5(2)*p534(2)-p5(3)*p534(
     & 3)
      ccr=zcr(id5)/(p534q*p534k0)
      ccl=zcl(id5)/(p534q*p534k0)
      do i3=1,2
* TL0 -- qu=p5,qd=p534,v=cfkz34(i3)%e,a=l5_34(i3)%a,c=l5_34(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz34(i3)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*
     & (cfkz34(i3)%e(2)*p534(3)-p534(2)*cfkz34(i3)%e(3))-p534k0*
     & (cfkz34(i3)%e(2)*p5(3)-p5(2)*cfkz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz34(i3)%e(3)*p5k0+p5(3)*cfkz34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz34(i3)%e(0)*p5(0)-cfkz34(i3)%e(1)*p5(1)-cfkz34(i3
     & )%e(2)*p5(2)-cfkz34(i3)%e(3)*p5(3)
      cvqd=cfkz34(i3)%e(0)*p534(0)-cfkz34(i3)%e(1)*p534(1)-cfkz3
     & 4(i3)%e(2)*p534(2)-cfkz34(i3)%e(3)*p534(3)
      cauxa=-cfkz34(i3)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cfkz34(i3)%ek0*p5(2)-p5k0*cfkz34(i3)%e(2)
      l5_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      l5_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      l5_34(i3)%c(1)=ccr*(cauxc+ceps_1)
      l5_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id5).ne.1.and.ineutri(id3).ne.1) then
      ccr=fcr(id5)/(p534q*p534k0)
      ccl=fcl(id5)/(p534q*p534k0)
      do i3=1,2
* TL0 -- qu=p5,qd=p534,v=cfkf34(i3)%e,a=l5_34(i3)%a,c=l5_34(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf34(i3)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*
     & (cfkf34(i3)%e(2)*p534(3)-p534(2)*cfkf34(i3)%e(3))-p534k0*
     & (cfkf34(i3)%e(2)*p5(3)-p5(2)*cfkf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf34(i3)%e(3)*p5k0+p5(3)*cfkf34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf34(i3)%e(0)*p5(0)-cfkf34(i3)%e(1)*p5(1)-cfkf34(i3
     & )%e(2)*p5(2)-cfkf34(i3)%e(3)*p5(3)
      cvqd=cfkf34(i3)%e(0)*p534(0)-cfkf34(i3)%e(1)*p534(1)-cfkf3
     & 4(i3)%e(2)*p534(2)-cfkf34(i3)%e(3)*p534(3)
      cauxa=-cfkf34(i3)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cfkf34(i3)%ek0*p5(2)-p5k0*cfkf34(i3)%e(2)
      l5_34(i3)%a(1)=l5_34(i3)%a(1)+ccr*(cauxa+ceps_0)
      l5_34(i3)%a(2)=l5_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      l5_34(i3)%c(1)=l5_34(i3)%c(1)+ccr*(cauxc+ceps_1)
      l5_34(i3)%c(2)=l5_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(5).eq.1) then
          do i3=1,2
            l5_78(i3)%a(1)=czero
            l5_78(i3)%a(2)=czero
            l5_78(i3)%c(1)=czero
            l5_78(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p5,q=p578
      quqd=p5(0)*p578(0)-p5(1)*p578(1)-p5(2)*p578(2)-p5(3)*p578(
     & 3)
      ccr=zcr(id5)/(p578q*p578k0)
      ccl=zcl(id5)/(p578q*p578k0)
      do i3=1,2
* TL0 -- qu=p5,qd=p578,v=cfkz78(i3)%e,a=l5_78(i3)%a,c=l5_78(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz78(i3)%ek0*(p5(2)*p578(3)-p578(2)*p5(3))+p5k0*
     & (cfkz78(i3)%e(2)*p578(3)-p578(2)*cfkz78(i3)%e(3))-p578k0*
     & (cfkz78(i3)%e(2)*p5(3)-p5(2)*cfkz78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz78(i3)%e(3)*p5k0+p5(3)*cfkz78(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz78(i3)%e(0)*p5(0)-cfkz78(i3)%e(1)*p5(1)-cfkz78(i3
     & )%e(2)*p5(2)-cfkz78(i3)%e(3)*p5(3)
      cvqd=cfkz78(i3)%e(0)*p578(0)-cfkz78(i3)%e(1)*p578(1)-cfkz7
     & 8(i3)%e(2)*p578(2)-cfkz78(i3)%e(3)*p578(3)
      cauxa=-cfkz78(i3)%ek0*quqd+p5k0*cvqd+p578k0*cvqu
      cauxc=+cfkz78(i3)%ek0*p5(2)-p5k0*cfkz78(i3)%e(2)
      l5_78(i3)%a(1)=ccr*(cauxa+ceps_0)
      l5_78(i3)%a(2)=ccl*(cauxa-ceps_0)
      l5_78(i3)%c(1)=ccr*(cauxc+ceps_1)
      l5_78(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id5).ne.1.and.ineutri(id7).ne.1) then
      ccr=fcr(id5)/(p578q*p578k0)
      ccl=fcl(id5)/(p578q*p578k0)
      do i3=1,2
* TL0 -- qu=p5,qd=p578,v=cfkf78(i3)%e,a=l5_78(i3)%a,c=l5_78(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf78(i3)%ek0*(p5(2)*p578(3)-p578(2)*p5(3))+p5k0*
     & (cfkf78(i3)%e(2)*p578(3)-p578(2)*cfkf78(i3)%e(3))-p578k0*
     & (cfkf78(i3)%e(2)*p5(3)-p5(2)*cfkf78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf78(i3)%e(3)*p5k0+p5(3)*cfkf78(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf78(i3)%e(0)*p5(0)-cfkf78(i3)%e(1)*p5(1)-cfkf78(i3
     & )%e(2)*p5(2)-cfkf78(i3)%e(3)*p5(3)
      cvqd=cfkf78(i3)%e(0)*p578(0)-cfkf78(i3)%e(1)*p578(1)-cfkf7
     & 8(i3)%e(2)*p578(2)-cfkf78(i3)%e(3)*p578(3)
      cauxa=-cfkf78(i3)%ek0*quqd+p5k0*cvqd+p578k0*cvqu
      cauxc=+cfkf78(i3)%ek0*p5(2)-p5k0*cfkf78(i3)%e(2)
      l5_78(i3)%a(1)=l5_78(i3)%a(1)+ccr*(cauxa+ceps_0)
      l5_78(i3)%a(2)=l5_78(i3)%a(2)+ccl*(cauxa-ceps_0)
      l5_78(i3)%c(1)=l5_78(i3)%c(1)+ccr*(cauxc+ceps_1)
      l5_78(i3)%c(2)=l5_78(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(7).eq.1) then
          do i3=1,2
            l7_12(i3)%a(1)=czero
            l7_12(i3)%a(2)=czero
            l7_12(i3)%c(1)=czero
            l7_12(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p7,q=p712
      quqd=p7(0)*p712(0)-p7(1)*p712(1)-p7(2)*p712(2)-p7(3)*p712(
     & 3)
      ccr=zcr(id7)/(p712q*p712k0)
      ccl=zcl(id7)/(p712q*p712k0)
      do i3=1,2
* TL0 -- qu=p7,qd=p712,v=cfkz12(i3)%e,a=l7_12(i3)%a,c=l7_12(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz12(i3)%ek0*(p7(2)*p712(3)-p712(2)*p7(3))+p7k0*
     & (cfkz12(i3)%e(2)*p712(3)-p712(2)*cfkz12(i3)%e(3))-p712k0*
     & (cfkz12(i3)%e(2)*p7(3)-p7(2)*cfkz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz12(i3)%e(3)*p7k0+p7(3)*cfkz12(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz12(i3)%e(0)*p7(0)-cfkz12(i3)%e(1)*p7(1)-cfkz12(i3
     & )%e(2)*p7(2)-cfkz12(i3)%e(3)*p7(3)
      cvqd=cfkz12(i3)%e(0)*p712(0)-cfkz12(i3)%e(1)*p712(1)-cfkz1
     & 2(i3)%e(2)*p712(2)-cfkz12(i3)%e(3)*p712(3)
      cauxa=-cfkz12(i3)%ek0*quqd+p7k0*cvqd+p712k0*cvqu
      cauxc=+cfkz12(i3)%ek0*p7(2)-p7k0*cfkz12(i3)%e(2)
      l7_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      l7_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      l7_12(i3)%c(1)=ccr*(cauxc+ceps_1)
      l7_12(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id7).ne.1.and.ineutri(id1).ne.1) then
      ccr=fcr(id7)/(p712q*p712k0)
      ccl=fcl(id7)/(p712q*p712k0)
      do i3=1,2
* TL0 -- qu=p7,qd=p712,v=cfkf12(i3)%e,a=l7_12(i3)%a,c=l7_12(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf12(i3)%ek0*(p7(2)*p712(3)-p712(2)*p7(3))+p7k0*
     & (cfkf12(i3)%e(2)*p712(3)-p712(2)*cfkf12(i3)%e(3))-p712k0*
     & (cfkf12(i3)%e(2)*p7(3)-p7(2)*cfkf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf12(i3)%e(3)*p7k0+p7(3)*cfkf12(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf12(i3)%e(0)*p7(0)-cfkf12(i3)%e(1)*p7(1)-cfkf12(i3
     & )%e(2)*p7(2)-cfkf12(i3)%e(3)*p7(3)
      cvqd=cfkf12(i3)%e(0)*p712(0)-cfkf12(i3)%e(1)*p712(1)-cfkf1
     & 2(i3)%e(2)*p712(2)-cfkf12(i3)%e(3)*p712(3)
      cauxa=-cfkf12(i3)%ek0*quqd+p7k0*cvqd+p712k0*cvqu
      cauxc=+cfkf12(i3)%ek0*p7(2)-p7k0*cfkf12(i3)%e(2)
      l7_12(i3)%a(1)=l7_12(i3)%a(1)+ccr*(cauxa+ceps_0)
      l7_12(i3)%a(2)=l7_12(i3)%a(2)+ccl*(cauxa-ceps_0)
      l7_12(i3)%c(1)=l7_12(i3)%c(1)+ccr*(cauxc+ceps_1)
      l7_12(i3)%c(2)=l7_12(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(7).eq.1) then
          do i3=1,2
            l7_34(i3)%a(1)=czero
            l7_34(i3)%a(2)=czero
            l7_34(i3)%c(1)=czero
            l7_34(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p7,q=p734
      quqd=p7(0)*p734(0)-p7(1)*p734(1)-p7(2)*p734(2)-p7(3)*p734(
     & 3)
      ccr=zcr(id7)/(p734q*p734k0)
      ccl=zcl(id7)/(p734q*p734k0)
      do i3=1,2
* TL0 -- qu=p7,qd=p734,v=cfkz34(i3)%e,a=l7_34(i3)%a,c=l7_34(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz34(i3)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*
     & (cfkz34(i3)%e(2)*p734(3)-p734(2)*cfkz34(i3)%e(3))-p734k0*
     & (cfkz34(i3)%e(2)*p7(3)-p7(2)*cfkz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz34(i3)%e(3)*p7k0+p7(3)*cfkz34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz34(i3)%e(0)*p7(0)-cfkz34(i3)%e(1)*p7(1)-cfkz34(i3
     & )%e(2)*p7(2)-cfkz34(i3)%e(3)*p7(3)
      cvqd=cfkz34(i3)%e(0)*p734(0)-cfkz34(i3)%e(1)*p734(1)-cfkz3
     & 4(i3)%e(2)*p734(2)-cfkz34(i3)%e(3)*p734(3)
      cauxa=-cfkz34(i3)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cfkz34(i3)%ek0*p7(2)-p7k0*cfkz34(i3)%e(2)
      l7_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      l7_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      l7_34(i3)%c(1)=ccr*(cauxc+ceps_1)
      l7_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id7).ne.1.and.ineutri(id3).ne.1) then
      ccr=fcr(id7)/(p734q*p734k0)
      ccl=fcl(id7)/(p734q*p734k0)
      do i3=1,2
* TL0 -- qu=p7,qd=p734,v=cfkf34(i3)%e,a=l7_34(i3)%a,c=l7_34(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf34(i3)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*
     & (cfkf34(i3)%e(2)*p734(3)-p734(2)*cfkf34(i3)%e(3))-p734k0*
     & (cfkf34(i3)%e(2)*p7(3)-p7(2)*cfkf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf34(i3)%e(3)*p7k0+p7(3)*cfkf34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf34(i3)%e(0)*p7(0)-cfkf34(i3)%e(1)*p7(1)-cfkf34(i3
     & )%e(2)*p7(2)-cfkf34(i3)%e(3)*p7(3)
      cvqd=cfkf34(i3)%e(0)*p734(0)-cfkf34(i3)%e(1)*p734(1)-cfkf3
     & 4(i3)%e(2)*p734(2)-cfkf34(i3)%e(3)*p734(3)
      cauxa=-cfkf34(i3)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cfkf34(i3)%ek0*p7(2)-p7k0*cfkf34(i3)%e(2)
      l7_34(i3)%a(1)=l7_34(i3)%a(1)+ccr*(cauxa+ceps_0)
      l7_34(i3)%a(2)=l7_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      l7_34(i3)%c(1)=l7_34(i3)%c(1)+ccr*(cauxc+ceps_1)
      l7_34(i3)%c(2)=l7_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*zw                                                                     
  
  
        if(iresz(7).eq.1) then
          do i3=1,2
            l7_56(i3)%a(1)=czero
            l7_56(i3)%a(2)=czero
            l7_56(i3)%c(1)=czero
            l7_56(i3)%c(2)=czero
          enddo
        else
* quqd -- p=p7,q=p756
      quqd=p7(0)*p756(0)-p7(1)*p756(1)-p7(2)*p756(2)-p7(3)*p756(
     & 3)
      ccr=zcr(id7)/(p756q*p756k0)
      ccl=zcl(id7)/(p756q*p756k0)
      do i3=1,2
* TL0 -- qu=p7,qd=p756,v=cfkz56(i3)%e,a=l7_56(i3)%a,c=l7_56(i3)%c,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz56(i3)%ek0*(p7(2)*p756(3)-p756(2)*p7(3))+p7k0*
     & (cfkz56(i3)%e(2)*p756(3)-p756(2)*cfkz56(i3)%e(3))-p756k0*
     & (cfkz56(i3)%e(2)*p7(3)-p7(2)*cfkz56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz56(i3)%e(3)*p7k0+p7(3)*cfkz56(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkz56(i3)%e(0)*p7(0)-cfkz56(i3)%e(1)*p7(1)-cfkz56(i3
     & )%e(2)*p7(2)-cfkz56(i3)%e(3)*p7(3)
      cvqd=cfkz56(i3)%e(0)*p756(0)-cfkz56(i3)%e(1)*p756(1)-cfkz5
     & 6(i3)%e(2)*p756(2)-cfkz56(i3)%e(3)*p756(3)
      cauxa=-cfkz56(i3)%ek0*quqd+p7k0*cvqd+p756k0*cvqu
      cauxc=+cfkz56(i3)%ek0*p7(2)-p7k0*cfkz56(i3)%e(2)
      l7_56(i3)%a(1)=ccr*(cauxa+ceps_0)
      l7_56(i3)%a(2)=ccl*(cauxa-ceps_0)
      l7_56(i3)%c(1)=ccr*(cauxc+ceps_1)
      l7_56(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
          if (ineutri(id7).ne.1.and.ineutri(id5).ne.1) then
      ccr=fcr(id7)/(p756q*p756k0)
      ccl=fcl(id7)/(p756q*p756k0)
      do i3=1,2
* TL0 -- qu=p7,qd=p756,v=cfkf56(i3)%e,a=l7_56(i3)%a,c=l7_56(i3)%c,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf56(i3)%ek0*(p7(2)*p756(3)-p756(2)*p7(3))+p7k0*
     & (cfkf56(i3)%e(2)*p756(3)-p756(2)*cfkf56(i3)%e(3))-p756k0*
     & (cfkf56(i3)%e(2)*p7(3)-p7(2)*cfkf56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf56(i3)%e(3)*p7k0+p7(3)*cfkf56(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cfkf56(i3)%e(0)*p7(0)-cfkf56(i3)%e(1)*p7(1)-cfkf56(i3
     & )%e(2)*p7(2)-cfkf56(i3)%e(3)*p7(3)
      cvqd=cfkf56(i3)%e(0)*p756(0)-cfkf56(i3)%e(1)*p756(1)-cfkf5
     & 6(i3)%e(2)*p756(2)-cfkf56(i3)%e(3)*p756(3)
      cauxa=-cfkf56(i3)%ek0*quqd+p7k0*cvqd+p756k0*cvqu
      cauxc=+cfkf56(i3)%ek0*p7(2)-p7k0*cfkf56(i3)%e(2)
      l7_56(i3)%a(1)=l7_56(i3)%a(1)+ccr*(cauxa+ceps_0)
      l7_56(i3)%a(2)=l7_56(i3)%a(2)+ccl*(cauxa-ceps_0)
      l7_56(i3)%c(1)=l7_56(i3)%c(1)+ccr*(cauxc+ceps_1)
      l7_56(i3)%c(2)=l7_56(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
          endif
        endif
*zwend                                                                  
  
  
*                                                                       
* RIGHT INSERTIONS                                                      
*                                                                       
  
*right insertion momenta                                                
* quqd -- p=p234,q=p2
      quqd=p234(0)*p2(0)-p234(1)*p2(1)-p234(2)*p2(2)-p234(3)*p2(
     & 3)
      ccr=zcr(id2)/(p234q*p234k0)
      ccl=zcl(id2)/(p234q*p234k0)
      do i3=1,2
* TR0 -- qu=p234,qd=p2,v=cfkz34(i3)%e,a=r2_34(i3)%a,b=r2_34(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz34(i3)%ek0*(p234(2)*p2(3)-p2(2)*p234(3))+p234k
     & 0*(cfkz34(i3)%e(2)*p2(3)-p2(2)*cfkz34(i3)%e(3))-p2k0*(cfk
     & z34(i3)%e(2)*p234(3)-p234(2)*cfkz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz34(i3)%e(3)*p2k0+p2(3)*cfkz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz34(i3)%e(0)*p234(0)-cfkz34(i3)%e(1)*p234(1)-cfkz3
     & 4(i3)%e(2)*p234(2)-cfkz34(i3)%e(3)*p234(3)
      cvqd=cfkz34(i3)%e(0)*p2(0)-cfkz34(i3)%e(1)*p2(1)-cfkz34(i3
     & )%e(2)*p2(2)-cfkz34(i3)%e(3)*p2(3)
      cauxa=-cfkz34(i3)%ek0*quqd+p234k0*cvqd+p2k0*cvqu
      cauxb=-cfkz34(i3)%ek0*p2(2)+p2k0*cfkz34(i3)%e(2)
      r2_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      r2_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      r2_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      r2_34(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id2).ne.1.and.ineutri(id3).ne.1) then
      ccr=fcr(id2)/(p234q*p234k0)
      ccl=fcl(id2)/(p234q*p234k0)
      do i3=1,2
* TR0 -- qu=p234,qd=p2,v=cfkf34(i3)%e,a=r2_34(i3)%a,b=r2_34(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf34(i3)%ek0*(p234(2)*p2(3)-p2(2)*p234(3))+p234k
     & 0*(cfkf34(i3)%e(2)*p2(3)-p2(2)*cfkf34(i3)%e(3))-p2k0*(cfk
     & f34(i3)%e(2)*p234(3)-p234(2)*cfkf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf34(i3)%e(3)*p2k0+p2(3)*cfkf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf34(i3)%e(0)*p234(0)-cfkf34(i3)%e(1)*p234(1)-cfkf3
     & 4(i3)%e(2)*p234(2)-cfkf34(i3)%e(3)*p234(3)
      cvqd=cfkf34(i3)%e(0)*p2(0)-cfkf34(i3)%e(1)*p2(1)-cfkf34(i3
     & )%e(2)*p2(2)-cfkf34(i3)%e(3)*p2(3)
      cauxa=-cfkf34(i3)%ek0*quqd+p234k0*cvqd+p2k0*cvqu
      cauxb=-cfkf34(i3)%ek0*p2(2)+p2k0*cfkf34(i3)%e(2)
      r2_34(i3)%a(1)=r2_34(i3)%a(1)+ccr*(cauxa+ceps_0)
      r2_34(i3)%a(2)=r2_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      r2_34(i3)%b(1)=r2_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      r2_34(i3)%b(2)=r2_34(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p256,q=p2
      quqd=p256(0)*p2(0)-p256(1)*p2(1)-p256(2)*p2(2)-p256(3)*p2(
     & 3)
      ccr=zcr(id2)/(p256q*p256k0)
      ccl=zcl(id2)/(p256q*p256k0)
      do i3=1,2
* TR0 -- qu=p256,qd=p2,v=cfkz56(i3)%e,a=r2_56(i3)%a,b=r2_56(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz56(i3)%ek0*(p256(2)*p2(3)-p2(2)*p256(3))+p256k
     & 0*(cfkz56(i3)%e(2)*p2(3)-p2(2)*cfkz56(i3)%e(3))-p2k0*(cfk
     & z56(i3)%e(2)*p256(3)-p256(2)*cfkz56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz56(i3)%e(3)*p2k0+p2(3)*cfkz56(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz56(i3)%e(0)*p256(0)-cfkz56(i3)%e(1)*p256(1)-cfkz5
     & 6(i3)%e(2)*p256(2)-cfkz56(i3)%e(3)*p256(3)
      cvqd=cfkz56(i3)%e(0)*p2(0)-cfkz56(i3)%e(1)*p2(1)-cfkz56(i3
     & )%e(2)*p2(2)-cfkz56(i3)%e(3)*p2(3)
      cauxa=-cfkz56(i3)%ek0*quqd+p256k0*cvqd+p2k0*cvqu
      cauxb=-cfkz56(i3)%ek0*p2(2)+p2k0*cfkz56(i3)%e(2)
      r2_56(i3)%a(1)=ccr*(cauxa+ceps_0)
      r2_56(i3)%a(2)=ccl*(cauxa-ceps_0)
      r2_56(i3)%b(1)=ccl*(cauxb-ceps_2)
      r2_56(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id2).ne.1.and.ineutri(id5).ne.1) then
      ccr=fcr(id2)/(p256q*p256k0)
      ccl=fcl(id2)/(p256q*p256k0)
      do i3=1,2
* TR0 -- qu=p256,qd=p2,v=cfkf56(i3)%e,a=r2_56(i3)%a,b=r2_56(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf56(i3)%ek0*(p256(2)*p2(3)-p2(2)*p256(3))+p256k
     & 0*(cfkf56(i3)%e(2)*p2(3)-p2(2)*cfkf56(i3)%e(3))-p2k0*(cfk
     & f56(i3)%e(2)*p256(3)-p256(2)*cfkf56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf56(i3)%e(3)*p2k0+p2(3)*cfkf56(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf56(i3)%e(0)*p256(0)-cfkf56(i3)%e(1)*p256(1)-cfkf5
     & 6(i3)%e(2)*p256(2)-cfkf56(i3)%e(3)*p256(3)
      cvqd=cfkf56(i3)%e(0)*p2(0)-cfkf56(i3)%e(1)*p2(1)-cfkf56(i3
     & )%e(2)*p2(2)-cfkf56(i3)%e(3)*p2(3)
      cauxa=-cfkf56(i3)%ek0*quqd+p256k0*cvqd+p2k0*cvqu
      cauxb=-cfkf56(i3)%ek0*p2(2)+p2k0*cfkf56(i3)%e(2)
      r2_56(i3)%a(1)=r2_56(i3)%a(1)+ccr*(cauxa+ceps_0)
      r2_56(i3)%a(2)=r2_56(i3)%a(2)+ccl*(cauxa-ceps_0)
      r2_56(i3)%b(1)=r2_56(i3)%b(1)+ccl*(cauxb-ceps_2)
      r2_56(i3)%b(2)=r2_56(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p278,q=p2
      quqd=p278(0)*p2(0)-p278(1)*p2(1)-p278(2)*p2(2)-p278(3)*p2(
     & 3)
      ccr=zcr(id2)/(p278q*p278k0)
      ccl=zcl(id2)/(p278q*p278k0)
      do i3=1,2
* TR0 -- qu=p278,qd=p2,v=cfkz78(i3)%e,a=r2_78(i3)%a,b=r2_78(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz78(i3)%ek0*(p278(2)*p2(3)-p2(2)*p278(3))+p278k
     & 0*(cfkz78(i3)%e(2)*p2(3)-p2(2)*cfkz78(i3)%e(3))-p2k0*(cfk
     & z78(i3)%e(2)*p278(3)-p278(2)*cfkz78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz78(i3)%e(3)*p2k0+p2(3)*cfkz78(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz78(i3)%e(0)*p278(0)-cfkz78(i3)%e(1)*p278(1)-cfkz7
     & 8(i3)%e(2)*p278(2)-cfkz78(i3)%e(3)*p278(3)
      cvqd=cfkz78(i3)%e(0)*p2(0)-cfkz78(i3)%e(1)*p2(1)-cfkz78(i3
     & )%e(2)*p2(2)-cfkz78(i3)%e(3)*p2(3)
      cauxa=-cfkz78(i3)%ek0*quqd+p278k0*cvqd+p2k0*cvqu
      cauxb=-cfkz78(i3)%ek0*p2(2)+p2k0*cfkz78(i3)%e(2)
      r2_78(i3)%a(1)=ccr*(cauxa+ceps_0)
      r2_78(i3)%a(2)=ccl*(cauxa-ceps_0)
      r2_78(i3)%b(1)=ccl*(cauxb-ceps_2)
      r2_78(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id2).ne.1.and.ineutri(id7).ne.1) then
      ccr=fcr(id2)/(p278q*p278k0)
      ccl=fcl(id2)/(p278q*p278k0)
      do i3=1,2
* TR0 -- qu=p278,qd=p2,v=cfkf78(i3)%e,a=r2_78(i3)%a,b=r2_78(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf78(i3)%ek0*(p278(2)*p2(3)-p2(2)*p278(3))+p278k
     & 0*(cfkf78(i3)%e(2)*p2(3)-p2(2)*cfkf78(i3)%e(3))-p2k0*(cfk
     & f78(i3)%e(2)*p278(3)-p278(2)*cfkf78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf78(i3)%e(3)*p2k0+p2(3)*cfkf78(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf78(i3)%e(0)*p278(0)-cfkf78(i3)%e(1)*p278(1)-cfkf7
     & 8(i3)%e(2)*p278(2)-cfkf78(i3)%e(3)*p278(3)
      cvqd=cfkf78(i3)%e(0)*p2(0)-cfkf78(i3)%e(1)*p2(1)-cfkf78(i3
     & )%e(2)*p2(2)-cfkf78(i3)%e(3)*p2(3)
      cauxa=-cfkf78(i3)%ek0*quqd+p278k0*cvqd+p2k0*cvqu
      cauxb=-cfkf78(i3)%ek0*p2(2)+p2k0*cfkf78(i3)%e(2)
      r2_78(i3)%a(1)=r2_78(i3)%a(1)+ccr*(cauxa+ceps_0)
      r2_78(i3)%a(2)=r2_78(i3)%a(2)+ccl*(cauxa-ceps_0)
      r2_78(i3)%b(1)=r2_78(i3)%b(1)+ccl*(cauxb-ceps_2)
      r2_78(i3)%b(2)=r2_78(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p412,q=p4
      quqd=p412(0)*p4(0)-p412(1)*p4(1)-p412(2)*p4(2)-p412(3)*p4(
     & 3)
      ccr=zcr(id4)/(p412q*p412k0)
      ccl=zcl(id4)/(p412q*p412k0)
      do i3=1,2
* TR0 -- qu=p412,qd=p4,v=cfkz12(i3)%e,a=r4_12(i3)%a,b=r4_12(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz12(i3)%ek0*(p412(2)*p4(3)-p4(2)*p412(3))+p412k
     & 0*(cfkz12(i3)%e(2)*p4(3)-p4(2)*cfkz12(i3)%e(3))-p4k0*(cfk
     & z12(i3)%e(2)*p412(3)-p412(2)*cfkz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz12(i3)%e(3)*p4k0+p4(3)*cfkz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz12(i3)%e(0)*p412(0)-cfkz12(i3)%e(1)*p412(1)-cfkz1
     & 2(i3)%e(2)*p412(2)-cfkz12(i3)%e(3)*p412(3)
      cvqd=cfkz12(i3)%e(0)*p4(0)-cfkz12(i3)%e(1)*p4(1)-cfkz12(i3
     & )%e(2)*p4(2)-cfkz12(i3)%e(3)*p4(3)
      cauxa=-cfkz12(i3)%ek0*quqd+p412k0*cvqd+p4k0*cvqu
      cauxb=-cfkz12(i3)%ek0*p4(2)+p4k0*cfkz12(i3)%e(2)
      r4_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      r4_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      r4_12(i3)%b(1)=ccl*(cauxb-ceps_2)
      r4_12(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id4).ne.1.and.ineutri(id1).ne.1) then
      ccr=fcr(id4)/(p412q*p412k0)
      ccl=fcl(id4)/(p412q*p412k0)
      do i3=1,2
* TR0 -- qu=p412,qd=p4,v=cfkf12(i3)%e,a=r4_12(i3)%a,b=r4_12(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf12(i3)%ek0*(p412(2)*p4(3)-p4(2)*p412(3))+p412k
     & 0*(cfkf12(i3)%e(2)*p4(3)-p4(2)*cfkf12(i3)%e(3))-p4k0*(cfk
     & f12(i3)%e(2)*p412(3)-p412(2)*cfkf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf12(i3)%e(3)*p4k0+p4(3)*cfkf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf12(i3)%e(0)*p412(0)-cfkf12(i3)%e(1)*p412(1)-cfkf1
     & 2(i3)%e(2)*p412(2)-cfkf12(i3)%e(3)*p412(3)
      cvqd=cfkf12(i3)%e(0)*p4(0)-cfkf12(i3)%e(1)*p4(1)-cfkf12(i3
     & )%e(2)*p4(2)-cfkf12(i3)%e(3)*p4(3)
      cauxa=-cfkf12(i3)%ek0*quqd+p412k0*cvqd+p4k0*cvqu
      cauxb=-cfkf12(i3)%ek0*p4(2)+p4k0*cfkf12(i3)%e(2)
      r4_12(i3)%a(1)=r4_12(i3)%a(1)+ccr*(cauxa+ceps_0)
      r4_12(i3)%a(2)=r4_12(i3)%a(2)+ccl*(cauxa-ceps_0)
      r4_12(i3)%b(1)=r4_12(i3)%b(1)+ccl*(cauxb-ceps_2)
      r4_12(i3)%b(2)=r4_12(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p456,q=p4
      quqd=p456(0)*p4(0)-p456(1)*p4(1)-p456(2)*p4(2)-p456(3)*p4(
     & 3)
      ccr=zcr(id4)/(p456q*p456k0)
      ccl=zcl(id4)/(p456q*p456k0)
      do i3=1,2
* TR0 -- qu=p456,qd=p4,v=cfkz56(i3)%e,a=r4_56(i3)%a,b=r4_56(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz56(i3)%ek0*(p456(2)*p4(3)-p4(2)*p456(3))+p456k
     & 0*(cfkz56(i3)%e(2)*p4(3)-p4(2)*cfkz56(i3)%e(3))-p4k0*(cfk
     & z56(i3)%e(2)*p456(3)-p456(2)*cfkz56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz56(i3)%e(3)*p4k0+p4(3)*cfkz56(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz56(i3)%e(0)*p456(0)-cfkz56(i3)%e(1)*p456(1)-cfkz5
     & 6(i3)%e(2)*p456(2)-cfkz56(i3)%e(3)*p456(3)
      cvqd=cfkz56(i3)%e(0)*p4(0)-cfkz56(i3)%e(1)*p4(1)-cfkz56(i3
     & )%e(2)*p4(2)-cfkz56(i3)%e(3)*p4(3)
      cauxa=-cfkz56(i3)%ek0*quqd+p456k0*cvqd+p4k0*cvqu
      cauxb=-cfkz56(i3)%ek0*p4(2)+p4k0*cfkz56(i3)%e(2)
      r4_56(i3)%a(1)=ccr*(cauxa+ceps_0)
      r4_56(i3)%a(2)=ccl*(cauxa-ceps_0)
      r4_56(i3)%b(1)=ccl*(cauxb-ceps_2)
      r4_56(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id4).ne.1.and.ineutri(id5).ne.1) then
      ccr=fcr(id4)/(p456q*p456k0)
      ccl=fcl(id4)/(p456q*p456k0)
      do i3=1,2
* TR0 -- qu=p456,qd=p4,v=cfkf56(i3)%e,a=r4_56(i3)%a,b=r4_56(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf56(i3)%ek0*(p456(2)*p4(3)-p4(2)*p456(3))+p456k
     & 0*(cfkf56(i3)%e(2)*p4(3)-p4(2)*cfkf56(i3)%e(3))-p4k0*(cfk
     & f56(i3)%e(2)*p456(3)-p456(2)*cfkf56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf56(i3)%e(3)*p4k0+p4(3)*cfkf56(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf56(i3)%e(0)*p456(0)-cfkf56(i3)%e(1)*p456(1)-cfkf5
     & 6(i3)%e(2)*p456(2)-cfkf56(i3)%e(3)*p456(3)
      cvqd=cfkf56(i3)%e(0)*p4(0)-cfkf56(i3)%e(1)*p4(1)-cfkf56(i3
     & )%e(2)*p4(2)-cfkf56(i3)%e(3)*p4(3)
      cauxa=-cfkf56(i3)%ek0*quqd+p456k0*cvqd+p4k0*cvqu
      cauxb=-cfkf56(i3)%ek0*p4(2)+p4k0*cfkf56(i3)%e(2)
      r4_56(i3)%a(1)=r4_56(i3)%a(1)+ccr*(cauxa+ceps_0)
      r4_56(i3)%a(2)=r4_56(i3)%a(2)+ccl*(cauxa-ceps_0)
      r4_56(i3)%b(1)=r4_56(i3)%b(1)+ccl*(cauxb-ceps_2)
      r4_56(i3)%b(2)=r4_56(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p478,q=p4
      quqd=p478(0)*p4(0)-p478(1)*p4(1)-p478(2)*p4(2)-p478(3)*p4(
     & 3)
      ccr=zcr(id4)/(p478q*p478k0)
      ccl=zcl(id4)/(p478q*p478k0)
      do i3=1,2
* TR0 -- qu=p478,qd=p4,v=cfkz78(i3)%e,a=r4_78(i3)%a,b=r4_78(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz78(i3)%ek0*(p478(2)*p4(3)-p4(2)*p478(3))+p478k
     & 0*(cfkz78(i3)%e(2)*p4(3)-p4(2)*cfkz78(i3)%e(3))-p4k0*(cfk
     & z78(i3)%e(2)*p478(3)-p478(2)*cfkz78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz78(i3)%e(3)*p4k0+p4(3)*cfkz78(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz78(i3)%e(0)*p478(0)-cfkz78(i3)%e(1)*p478(1)-cfkz7
     & 8(i3)%e(2)*p478(2)-cfkz78(i3)%e(3)*p478(3)
      cvqd=cfkz78(i3)%e(0)*p4(0)-cfkz78(i3)%e(1)*p4(1)-cfkz78(i3
     & )%e(2)*p4(2)-cfkz78(i3)%e(3)*p4(3)
      cauxa=-cfkz78(i3)%ek0*quqd+p478k0*cvqd+p4k0*cvqu
      cauxb=-cfkz78(i3)%ek0*p4(2)+p4k0*cfkz78(i3)%e(2)
      r4_78(i3)%a(1)=ccr*(cauxa+ceps_0)
      r4_78(i3)%a(2)=ccl*(cauxa-ceps_0)
      r4_78(i3)%b(1)=ccl*(cauxb-ceps_2)
      r4_78(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id4).ne.1.and.ineutri(id7).ne.1) then
      ccr=fcr(id4)/(p478q*p478k0)
      ccl=fcl(id4)/(p478q*p478k0)
      do i3=1,2
* TR0 -- qu=p478,qd=p4,v=cfkf78(i3)%e,a=r4_78(i3)%a,b=r4_78(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf78(i3)%ek0*(p478(2)*p4(3)-p4(2)*p478(3))+p478k
     & 0*(cfkf78(i3)%e(2)*p4(3)-p4(2)*cfkf78(i3)%e(3))-p4k0*(cfk
     & f78(i3)%e(2)*p478(3)-p478(2)*cfkf78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf78(i3)%e(3)*p4k0+p4(3)*cfkf78(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf78(i3)%e(0)*p478(0)-cfkf78(i3)%e(1)*p478(1)-cfkf7
     & 8(i3)%e(2)*p478(2)-cfkf78(i3)%e(3)*p478(3)
      cvqd=cfkf78(i3)%e(0)*p4(0)-cfkf78(i3)%e(1)*p4(1)-cfkf78(i3
     & )%e(2)*p4(2)-cfkf78(i3)%e(3)*p4(3)
      cauxa=-cfkf78(i3)%ek0*quqd+p478k0*cvqd+p4k0*cvqu
      cauxb=-cfkf78(i3)%ek0*p4(2)+p4k0*cfkf78(i3)%e(2)
      r4_78(i3)%a(1)=r4_78(i3)%a(1)+ccr*(cauxa+ceps_0)
      r4_78(i3)%a(2)=r4_78(i3)%a(2)+ccl*(cauxa-ceps_0)
      r4_78(i3)%b(1)=r4_78(i3)%b(1)+ccl*(cauxb-ceps_2)
      r4_78(i3)%b(2)=r4_78(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p612,q=p6
      quqd=p612(0)*p6(0)-p612(1)*p6(1)-p612(2)*p6(2)-p612(3)*p6(
     & 3)
      ccr=zcr(id6)/(p612q*p612k0)
      ccl=zcl(id6)/(p612q*p612k0)
      do i3=1,2
* TR0 -- qu=p612,qd=p6,v=cfkz12(i3)%e,a=r6_12(i3)%a,b=r6_12(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz12(i3)%ek0*(p612(2)*p6(3)-p6(2)*p612(3))+p612k
     & 0*(cfkz12(i3)%e(2)*p6(3)-p6(2)*cfkz12(i3)%e(3))-p6k0*(cfk
     & z12(i3)%e(2)*p612(3)-p612(2)*cfkz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz12(i3)%e(3)*p6k0+p6(3)*cfkz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz12(i3)%e(0)*p612(0)-cfkz12(i3)%e(1)*p612(1)-cfkz1
     & 2(i3)%e(2)*p612(2)-cfkz12(i3)%e(3)*p612(3)
      cvqd=cfkz12(i3)%e(0)*p6(0)-cfkz12(i3)%e(1)*p6(1)-cfkz12(i3
     & )%e(2)*p6(2)-cfkz12(i3)%e(3)*p6(3)
      cauxa=-cfkz12(i3)%ek0*quqd+p612k0*cvqd+p6k0*cvqu
      cauxb=-cfkz12(i3)%ek0*p6(2)+p6k0*cfkz12(i3)%e(2)
      r6_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      r6_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      r6_12(i3)%b(1)=ccl*(cauxb-ceps_2)
      r6_12(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id6).ne.1.and.ineutri(id1).ne.1) then
      ccr=fcr(id6)/(p612q*p612k0)
      ccl=fcl(id6)/(p612q*p612k0)
      do i3=1,2
* TR0 -- qu=p612,qd=p6,v=cfkf12(i3)%e,a=r6_12(i3)%a,b=r6_12(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf12(i3)%ek0*(p612(2)*p6(3)-p6(2)*p612(3))+p612k
     & 0*(cfkf12(i3)%e(2)*p6(3)-p6(2)*cfkf12(i3)%e(3))-p6k0*(cfk
     & f12(i3)%e(2)*p612(3)-p612(2)*cfkf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf12(i3)%e(3)*p6k0+p6(3)*cfkf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf12(i3)%e(0)*p612(0)-cfkf12(i3)%e(1)*p612(1)-cfkf1
     & 2(i3)%e(2)*p612(2)-cfkf12(i3)%e(3)*p612(3)
      cvqd=cfkf12(i3)%e(0)*p6(0)-cfkf12(i3)%e(1)*p6(1)-cfkf12(i3
     & )%e(2)*p6(2)-cfkf12(i3)%e(3)*p6(3)
      cauxa=-cfkf12(i3)%ek0*quqd+p612k0*cvqd+p6k0*cvqu
      cauxb=-cfkf12(i3)%ek0*p6(2)+p6k0*cfkf12(i3)%e(2)
      r6_12(i3)%a(1)=r6_12(i3)%a(1)+ccr*(cauxa+ceps_0)
      r6_12(i3)%a(2)=r6_12(i3)%a(2)+ccl*(cauxa-ceps_0)
      r6_12(i3)%b(1)=r6_12(i3)%b(1)+ccl*(cauxb-ceps_2)
      r6_12(i3)%b(2)=r6_12(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p634,q=p6
      quqd=p634(0)*p6(0)-p634(1)*p6(1)-p634(2)*p6(2)-p634(3)*p6(
     & 3)
      ccr=zcr(id6)/(p634q*p634k0)
      ccl=zcl(id6)/(p634q*p634k0)
      do i3=1,2
* TR0 -- qu=p634,qd=p6,v=cfkz34(i3)%e,a=r6_34(i3)%a,b=r6_34(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz34(i3)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k
     & 0*(cfkz34(i3)%e(2)*p6(3)-p6(2)*cfkz34(i3)%e(3))-p6k0*(cfk
     & z34(i3)%e(2)*p634(3)-p634(2)*cfkz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz34(i3)%e(3)*p6k0+p6(3)*cfkz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz34(i3)%e(0)*p634(0)-cfkz34(i3)%e(1)*p634(1)-cfkz3
     & 4(i3)%e(2)*p634(2)-cfkz34(i3)%e(3)*p634(3)
      cvqd=cfkz34(i3)%e(0)*p6(0)-cfkz34(i3)%e(1)*p6(1)-cfkz34(i3
     & )%e(2)*p6(2)-cfkz34(i3)%e(3)*p6(3)
      cauxa=-cfkz34(i3)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cfkz34(i3)%ek0*p6(2)+p6k0*cfkz34(i3)%e(2)
      r6_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      r6_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      r6_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      r6_34(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id6).ne.1.and.ineutri(id3).ne.1) then
      ccr=fcr(id6)/(p634q*p634k0)
      ccl=fcl(id6)/(p634q*p634k0)
      do i3=1,2
* TR0 -- qu=p634,qd=p6,v=cfkf34(i3)%e,a=r6_34(i3)%a,b=r6_34(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf34(i3)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k
     & 0*(cfkf34(i3)%e(2)*p6(3)-p6(2)*cfkf34(i3)%e(3))-p6k0*(cfk
     & f34(i3)%e(2)*p634(3)-p634(2)*cfkf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf34(i3)%e(3)*p6k0+p6(3)*cfkf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf34(i3)%e(0)*p634(0)-cfkf34(i3)%e(1)*p634(1)-cfkf3
     & 4(i3)%e(2)*p634(2)-cfkf34(i3)%e(3)*p634(3)
      cvqd=cfkf34(i3)%e(0)*p6(0)-cfkf34(i3)%e(1)*p6(1)-cfkf34(i3
     & )%e(2)*p6(2)-cfkf34(i3)%e(3)*p6(3)
      cauxa=-cfkf34(i3)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cfkf34(i3)%ek0*p6(2)+p6k0*cfkf34(i3)%e(2)
      r6_34(i3)%a(1)=r6_34(i3)%a(1)+ccr*(cauxa+ceps_0)
      r6_34(i3)%a(2)=r6_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      r6_34(i3)%b(1)=r6_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      r6_34(i3)%b(2)=r6_34(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p678,q=p6
      quqd=p678(0)*p6(0)-p678(1)*p6(1)-p678(2)*p6(2)-p678(3)*p6(
     & 3)
      ccr=zcr(id6)/(p678q*p678k0)
      ccl=zcl(id6)/(p678q*p678k0)
      do i3=1,2
* TR0 -- qu=p678,qd=p6,v=cfkz78(i3)%e,a=r6_78(i3)%a,b=r6_78(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz78(i3)%ek0*(p678(2)*p6(3)-p6(2)*p678(3))+p678k
     & 0*(cfkz78(i3)%e(2)*p6(3)-p6(2)*cfkz78(i3)%e(3))-p6k0*(cfk
     & z78(i3)%e(2)*p678(3)-p678(2)*cfkz78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz78(i3)%e(3)*p6k0+p6(3)*cfkz78(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz78(i3)%e(0)*p678(0)-cfkz78(i3)%e(1)*p678(1)-cfkz7
     & 8(i3)%e(2)*p678(2)-cfkz78(i3)%e(3)*p678(3)
      cvqd=cfkz78(i3)%e(0)*p6(0)-cfkz78(i3)%e(1)*p6(1)-cfkz78(i3
     & )%e(2)*p6(2)-cfkz78(i3)%e(3)*p6(3)
      cauxa=-cfkz78(i3)%ek0*quqd+p678k0*cvqd+p6k0*cvqu
      cauxb=-cfkz78(i3)%ek0*p6(2)+p6k0*cfkz78(i3)%e(2)
      r6_78(i3)%a(1)=ccr*(cauxa+ceps_0)
      r6_78(i3)%a(2)=ccl*(cauxa-ceps_0)
      r6_78(i3)%b(1)=ccl*(cauxb-ceps_2)
      r6_78(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id6).ne.1.and.ineutri(id7).ne.1) then
      ccr=fcr(id6)/(p678q*p678k0)
      ccl=fcl(id6)/(p678q*p678k0)
      do i3=1,2
* TR0 -- qu=p678,qd=p6,v=cfkf78(i3)%e,a=r6_78(i3)%a,b=r6_78(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf78(i3)%ek0*(p678(2)*p6(3)-p6(2)*p678(3))+p678k
     & 0*(cfkf78(i3)%e(2)*p6(3)-p6(2)*cfkf78(i3)%e(3))-p6k0*(cfk
     & f78(i3)%e(2)*p678(3)-p678(2)*cfkf78(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf78(i3)%e(3)*p6k0+p6(3)*cfkf78(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf78(i3)%e(0)*p678(0)-cfkf78(i3)%e(1)*p678(1)-cfkf7
     & 8(i3)%e(2)*p678(2)-cfkf78(i3)%e(3)*p678(3)
      cvqd=cfkf78(i3)%e(0)*p6(0)-cfkf78(i3)%e(1)*p6(1)-cfkf78(i3
     & )%e(2)*p6(2)-cfkf78(i3)%e(3)*p6(3)
      cauxa=-cfkf78(i3)%ek0*quqd+p678k0*cvqd+p6k0*cvqu
      cauxb=-cfkf78(i3)%ek0*p6(2)+p6k0*cfkf78(i3)%e(2)
      r6_78(i3)%a(1)=r6_78(i3)%a(1)+ccr*(cauxa+ceps_0)
      r6_78(i3)%a(2)=r6_78(i3)%a(2)+ccl*(cauxa-ceps_0)
      r6_78(i3)%b(1)=r6_78(i3)%b(1)+ccl*(cauxb-ceps_2)
      r6_78(i3)%b(2)=r6_78(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p812,q=p8
      quqd=p812(0)*p8(0)-p812(1)*p8(1)-p812(2)*p8(2)-p812(3)*p8(
     & 3)
      ccr=zcr(id8)/(p812q*p812k0)
      ccl=zcl(id8)/(p812q*p812k0)
      do i3=1,2
* TR0 -- qu=p812,qd=p8,v=cfkz12(i3)%e,a=r8_12(i3)%a,b=r8_12(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz12(i3)%ek0*(p812(2)*p8(3)-p8(2)*p812(3))+p812k
     & 0*(cfkz12(i3)%e(2)*p8(3)-p8(2)*cfkz12(i3)%e(3))-p8k0*(cfk
     & z12(i3)%e(2)*p812(3)-p812(2)*cfkz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz12(i3)%e(3)*p8k0+p8(3)*cfkz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz12(i3)%e(0)*p812(0)-cfkz12(i3)%e(1)*p812(1)-cfkz1
     & 2(i3)%e(2)*p812(2)-cfkz12(i3)%e(3)*p812(3)
      cvqd=cfkz12(i3)%e(0)*p8(0)-cfkz12(i3)%e(1)*p8(1)-cfkz12(i3
     & )%e(2)*p8(2)-cfkz12(i3)%e(3)*p8(3)
      cauxa=-cfkz12(i3)%ek0*quqd+p812k0*cvqd+p8k0*cvqu
      cauxb=-cfkz12(i3)%ek0*p8(2)+p8k0*cfkz12(i3)%e(2)
      r8_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      r8_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      r8_12(i3)%b(1)=ccl*(cauxb-ceps_2)
      r8_12(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id8).ne.1.and.ineutri(id1).ne.1) then
      ccr=fcr(id8)/(p812q*p812k0)
      ccl=fcl(id8)/(p812q*p812k0)
      do i3=1,2
* TR0 -- qu=p812,qd=p8,v=cfkf12(i3)%e,a=r8_12(i3)%a,b=r8_12(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf12(i3)%ek0*(p812(2)*p8(3)-p8(2)*p812(3))+p812k
     & 0*(cfkf12(i3)%e(2)*p8(3)-p8(2)*cfkf12(i3)%e(3))-p8k0*(cfk
     & f12(i3)%e(2)*p812(3)-p812(2)*cfkf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf12(i3)%e(3)*p8k0+p8(3)*cfkf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf12(i3)%e(0)*p812(0)-cfkf12(i3)%e(1)*p812(1)-cfkf1
     & 2(i3)%e(2)*p812(2)-cfkf12(i3)%e(3)*p812(3)
      cvqd=cfkf12(i3)%e(0)*p8(0)-cfkf12(i3)%e(1)*p8(1)-cfkf12(i3
     & )%e(2)*p8(2)-cfkf12(i3)%e(3)*p8(3)
      cauxa=-cfkf12(i3)%ek0*quqd+p812k0*cvqd+p8k0*cvqu
      cauxb=-cfkf12(i3)%ek0*p8(2)+p8k0*cfkf12(i3)%e(2)
      r8_12(i3)%a(1)=r8_12(i3)%a(1)+ccr*(cauxa+ceps_0)
      r8_12(i3)%a(2)=r8_12(i3)%a(2)+ccl*(cauxa-ceps_0)
      r8_12(i3)%b(1)=r8_12(i3)%b(1)+ccl*(cauxb-ceps_2)
      r8_12(i3)%b(2)=r8_12(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p834,q=p8
      quqd=p834(0)*p8(0)-p834(1)*p8(1)-p834(2)*p8(2)-p834(3)*p8(
     & 3)
      ccr=zcr(id8)/(p834q*p834k0)
      ccl=zcl(id8)/(p834q*p834k0)
      do i3=1,2
* TR0 -- qu=p834,qd=p8,v=cfkz34(i3)%e,a=r8_34(i3)%a,b=r8_34(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz34(i3)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834k
     & 0*(cfkz34(i3)%e(2)*p8(3)-p8(2)*cfkz34(i3)%e(3))-p8k0*(cfk
     & z34(i3)%e(2)*p834(3)-p834(2)*cfkz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz34(i3)%e(3)*p8k0+p8(3)*cfkz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz34(i3)%e(0)*p834(0)-cfkz34(i3)%e(1)*p834(1)-cfkz3
     & 4(i3)%e(2)*p834(2)-cfkz34(i3)%e(3)*p834(3)
      cvqd=cfkz34(i3)%e(0)*p8(0)-cfkz34(i3)%e(1)*p8(1)-cfkz34(i3
     & )%e(2)*p8(2)-cfkz34(i3)%e(3)*p8(3)
      cauxa=-cfkz34(i3)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cfkz34(i3)%ek0*p8(2)+p8k0*cfkz34(i3)%e(2)
      r8_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      r8_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      r8_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      r8_34(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id8).ne.1.and.ineutri(id3).ne.1) then
      ccr=fcr(id8)/(p834q*p834k0)
      ccl=fcl(id8)/(p834q*p834k0)
      do i3=1,2
* TR0 -- qu=p834,qd=p8,v=cfkf34(i3)%e,a=r8_34(i3)%a,b=r8_34(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf34(i3)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834k
     & 0*(cfkf34(i3)%e(2)*p8(3)-p8(2)*cfkf34(i3)%e(3))-p8k0*(cfk
     & f34(i3)%e(2)*p834(3)-p834(2)*cfkf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf34(i3)%e(3)*p8k0+p8(3)*cfkf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf34(i3)%e(0)*p834(0)-cfkf34(i3)%e(1)*p834(1)-cfkf3
     & 4(i3)%e(2)*p834(2)-cfkf34(i3)%e(3)*p834(3)
      cvqd=cfkf34(i3)%e(0)*p8(0)-cfkf34(i3)%e(1)*p8(1)-cfkf34(i3
     & )%e(2)*p8(2)-cfkf34(i3)%e(3)*p8(3)
      cauxa=-cfkf34(i3)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cfkf34(i3)%ek0*p8(2)+p8k0*cfkf34(i3)%e(2)
      r8_34(i3)%a(1)=r8_34(i3)%a(1)+ccr*(cauxa+ceps_0)
      r8_34(i3)%a(2)=r8_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      r8_34(i3)%b(1)=r8_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      r8_34(i3)%b(2)=r8_34(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
* quqd -- p=p856,q=p8
      quqd=p856(0)*p8(0)-p856(1)*p8(1)-p856(2)*p8(2)-p856(3)*p8(
     & 3)
      ccr=zcr(id8)/(p856q*p856k0)
      ccl=zcl(id8)/(p856q*p856k0)
      do i3=1,2
* TR0 -- qu=p856,qd=p8,v=cfkz56(i3)%e,a=r8_56(i3)%a,b=r8_56(i3)%b,cr=ccr
* ,cl=ccl,nsum=0
      ceps_0=-cfkz56(i3)%ek0*(p856(2)*p8(3)-p8(2)*p856(3))+p856k
     & 0*(cfkz56(i3)%e(2)*p8(3)-p8(2)*cfkz56(i3)%e(3))-p8k0*(cfk
     & z56(i3)%e(2)*p856(3)-p856(2)*cfkz56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkz56(i3)%e(3)*p8k0+p8(3)*cfkz56(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz56(i3)%e(0)*p856(0)-cfkz56(i3)%e(1)*p856(1)-cfkz5
     & 6(i3)%e(2)*p856(2)-cfkz56(i3)%e(3)*p856(3)
      cvqd=cfkz56(i3)%e(0)*p8(0)-cfkz56(i3)%e(1)*p8(1)-cfkz56(i3
     & )%e(2)*p8(2)-cfkz56(i3)%e(3)*p8(3)
      cauxa=-cfkz56(i3)%ek0*quqd+p856k0*cvqd+p8k0*cvqu
      cauxb=-cfkz56(i3)%ek0*p8(2)+p8k0*cfkz56(i3)%e(2)
      r8_56(i3)%a(1)=ccr*(cauxa+ceps_0)
      r8_56(i3)%a(2)=ccl*(cauxa-ceps_0)
      r8_56(i3)%b(1)=ccl*(cauxb-ceps_2)
      r8_56(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
      if (ineutri(id8).ne.1.and.ineutri(id5).ne.1) then
      ccr=fcr(id8)/(p856q*p856k0)
      ccl=fcl(id8)/(p856q*p856k0)
      do i3=1,2
* TR0 -- qu=p856,qd=p8,v=cfkf56(i3)%e,a=r8_56(i3)%a,b=r8_56(i3)%b,cr=ccr
* ,cl=ccl,nsum=1
      ceps_0=-cfkf56(i3)%ek0*(p856(2)*p8(3)-p8(2)*p856(3))+p856k
     & 0*(cfkf56(i3)%e(2)*p8(3)-p8(2)*cfkf56(i3)%e(3))-p8k0*(cfk
     & f56(i3)%e(2)*p856(3)-p856(2)*cfkf56(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cfkf56(i3)%e(3)*p8k0+p8(3)*cfkf56(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf56(i3)%e(0)*p856(0)-cfkf56(i3)%e(1)*p856(1)-cfkf5
     & 6(i3)%e(2)*p856(2)-cfkf56(i3)%e(3)*p856(3)
      cvqd=cfkf56(i3)%e(0)*p8(0)-cfkf56(i3)%e(1)*p8(1)-cfkf56(i3
     & )%e(2)*p8(2)-cfkf56(i3)%e(3)*p8(3)
      cauxa=-cfkf56(i3)%ek0*quqd+p856k0*cvqd+p8k0*cvqu
      cauxb=-cfkf56(i3)%ek0*p8(2)+p8k0*cfkf56(i3)%e(2)
      r8_56(i3)%a(1)=r8_56(i3)%a(1)+ccr*(cauxa+ceps_0)
      r8_56(i3)%a(2)=r8_56(i3)%a(2)+ccl*(cauxa-ceps_0)
      r8_56(i3)%b(1)=r8_56(i3)%b(1)+ccl*(cauxb-ceps_2)
      r8_56(i3)%b(2)=r8_56(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      endif
  
  
*                                                                       
* MIDDLE INSERTIONS                                                     
*                                                                       
  
* quqd -- p=p134,q=p278
      quqd=p134(0)*p278(0)-p134(1)*p278(1)-p134(2)*p278(2)-p134(
     & 3)*p278(3)
      do i5=1,2
* T0 -- qu=p134,qd=p278,v=cfkz56(i5)%e,a=u134_56(i5)%a,b=u134_56(i5)%b,c
* =u134_56(i5)%c,d=u134_56(i5)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      ceps_0=-cfkz56(i5)%ek0*(p134(2)*p278(3)-p278(2)*p134(3))+p
     & 134k0*(cfkz56(i5)%e(2)*p278(3)-p278(2)*cfkz56(i5)%e(3))-p
     & 278k0*(cfkz56(i5)%e(2)*p134(3)-p134(2)*cfkz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz56(i5)%e(3)*p134k0+p134(3)*cfkz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz56(i5)%e(3)*p278k0+p278(3)*cfkz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz56(i5)%e(0)*p134(0)-cfkz56(i5)%e(1)*p134(1)-cfkz5
     & 6(i5)%e(2)*p134(2)-cfkz56(i5)%e(3)*p134(3)
      cvqd=cfkz56(i5)%e(0)*p278(0)-cfkz56(i5)%e(1)*p278(1)-cfkz5
     & 6(i5)%e(2)*p278(2)-cfkz56(i5)%e(3)*p278(3)
      cauxa=-cfkz56(i5)%ek0*quqd+p134k0*cvqd+p278k0*cvqu
      cauxb=-cfkz56(i5)%ek0*p278(2)+p278k0*cfkz56(i5)%e(2)
      cauxc=+cfkz56(i5)%ek0*p134(2)-p134k0*cfkz56(i5)%e(2)
      u134_56(i5)%a(1)=zcr(id1)*(cauxa+ceps_0)
      u134_56(i5)%a(2)=zcl(id1)*(cauxa-ceps_0)
      u134_56(i5)%b(1)=zcl(id1)*(cauxb-ceps_2)
      u134_56(i5)%b(2)=zcr(id1)*(-cauxb-ceps_2)
      u134_56(i5)%c(1)=zcr(id1)*(cauxc+ceps_1)
      u134_56(i5)%c(2)=zcl(id1)*(-cauxc+ceps_1)
      u134_56(i5)%d(1)=zcl(id1)*cfkz56(i5)%ek0
      u134_56(i5)%d(2)=zcr(id1)*cfkz56(i5)%ek0
      end do
  
      if (ineutri(id1).ne.1.and.ineutri(id5).ne.1) then
      do i5=1,2
* T0 -- qu=p134,qd=p278,v=cfkf56(i5)%e,a=u134_56(i5)%a,b=u134_56(i5)%b,c
* =u134_56(i5)%c,d=u134_56(i5)%d,cr=fcr(id1),cl=fcl(id1),nsum=1
      ceps_0=-cfkf56(i5)%ek0*(p134(2)*p278(3)-p278(2)*p134(3))+p
     & 134k0*(cfkf56(i5)%e(2)*p278(3)-p278(2)*cfkf56(i5)%e(3))-p
     & 278k0*(cfkf56(i5)%e(2)*p134(3)-p134(2)*cfkf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf56(i5)%e(3)*p134k0+p134(3)*cfkf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf56(i5)%e(3)*p278k0+p278(3)*cfkf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf56(i5)%e(0)*p134(0)-cfkf56(i5)%e(1)*p134(1)-cfkf5
     & 6(i5)%e(2)*p134(2)-cfkf56(i5)%e(3)*p134(3)
      cvqd=cfkf56(i5)%e(0)*p278(0)-cfkf56(i5)%e(1)*p278(1)-cfkf5
     & 6(i5)%e(2)*p278(2)-cfkf56(i5)%e(3)*p278(3)
      cauxa=-cfkf56(i5)%ek0*quqd+p134k0*cvqd+p278k0*cvqu
      cauxb=-cfkf56(i5)%ek0*p278(2)+p278k0*cfkf56(i5)%e(2)
      cauxc=+cfkf56(i5)%ek0*p134(2)-p134k0*cfkf56(i5)%e(2)
      u134_56(i5)%a(1)=u134_56(i5)%a(1)+fcr(id1)*(cauxa+ceps_0)
      u134_56(i5)%a(2)=u134_56(i5)%a(2)+fcl(id1)*(cauxa-ceps_0)
      u134_56(i5)%b(1)=u134_56(i5)%b(1)+fcl(id1)*(cauxb-ceps_2)
      u134_56(i5)%b(2)=u134_56(i5)%b(2)+fcr(id1)*(-cauxb-ceps_2)
      u134_56(i5)%c(1)=u134_56(i5)%c(1)+fcr(id1)*(cauxc+ceps_1)
      u134_56(i5)%c(2)=u134_56(i5)%c(2)+fcl(id1)*(-cauxc+ceps_1)
      u134_56(i5)%d(1)=u134_56(i5)%d(1)+fcl(id1)*cfkf56(i5)%ek0
      u134_56(i5)%d(2)=u134_56(i5)%d(2)+fcr(id1)*cfkf56(i5)%ek0
      end do
      endif
  
* quqd -- p=p134,q=p256
      quqd=p134(0)*p256(0)-p134(1)*p256(1)-p134(2)*p256(2)-p134(
     & 3)*p256(3)
      do i5=1,2
* T0 -- qu=p134,qd=p256,v=cfkz78(i5)%e,a=u134_78(i5)%a,b=u134_78(i5)%b,c
* =u134_78(i5)%c,d=u134_78(i5)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      ceps_0=-cfkz78(i5)%ek0*(p134(2)*p256(3)-p256(2)*p134(3))+p
     & 134k0*(cfkz78(i5)%e(2)*p256(3)-p256(2)*cfkz78(i5)%e(3))-p
     & 256k0*(cfkz78(i5)%e(2)*p134(3)-p134(2)*cfkz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz78(i5)%e(3)*p134k0+p134(3)*cfkz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz78(i5)%e(3)*p256k0+p256(3)*cfkz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz78(i5)%e(0)*p134(0)-cfkz78(i5)%e(1)*p134(1)-cfkz7
     & 8(i5)%e(2)*p134(2)-cfkz78(i5)%e(3)*p134(3)
      cvqd=cfkz78(i5)%e(0)*p256(0)-cfkz78(i5)%e(1)*p256(1)-cfkz7
     & 8(i5)%e(2)*p256(2)-cfkz78(i5)%e(3)*p256(3)
      cauxa=-cfkz78(i5)%ek0*quqd+p134k0*cvqd+p256k0*cvqu
      cauxb=-cfkz78(i5)%ek0*p256(2)+p256k0*cfkz78(i5)%e(2)
      cauxc=+cfkz78(i5)%ek0*p134(2)-p134k0*cfkz78(i5)%e(2)
      u134_78(i5)%a(1)=zcr(id1)*(cauxa+ceps_0)
      u134_78(i5)%a(2)=zcl(id1)*(cauxa-ceps_0)
      u134_78(i5)%b(1)=zcl(id1)*(cauxb-ceps_2)
      u134_78(i5)%b(2)=zcr(id1)*(-cauxb-ceps_2)
      u134_78(i5)%c(1)=zcr(id1)*(cauxc+ceps_1)
      u134_78(i5)%c(2)=zcl(id1)*(-cauxc+ceps_1)
      u134_78(i5)%d(1)=zcl(id1)*cfkz78(i5)%ek0
      u134_78(i5)%d(2)=zcr(id1)*cfkz78(i5)%ek0
      end do
  
      if (ineutri(id1).ne.1.and.ineutri(id7).ne.1) then
      do i5=1,2
* T0 -- qu=p134,qd=p256,v=cfkf78(i5)%e,a=u134_78(i5)%a,b=u134_78(i5)%b,c
* =u134_78(i5)%c,d=u134_78(i5)%d,cr=fcr(id1),cl=fcl(id1),nsum=1
      ceps_0=-cfkf78(i5)%ek0*(p134(2)*p256(3)-p256(2)*p134(3))+p
     & 134k0*(cfkf78(i5)%e(2)*p256(3)-p256(2)*cfkf78(i5)%e(3))-p
     & 256k0*(cfkf78(i5)%e(2)*p134(3)-p134(2)*cfkf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf78(i5)%e(3)*p134k0+p134(3)*cfkf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf78(i5)%e(3)*p256k0+p256(3)*cfkf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf78(i5)%e(0)*p134(0)-cfkf78(i5)%e(1)*p134(1)-cfkf7
     & 8(i5)%e(2)*p134(2)-cfkf78(i5)%e(3)*p134(3)
      cvqd=cfkf78(i5)%e(0)*p256(0)-cfkf78(i5)%e(1)*p256(1)-cfkf7
     & 8(i5)%e(2)*p256(2)-cfkf78(i5)%e(3)*p256(3)
      cauxa=-cfkf78(i5)%ek0*quqd+p134k0*cvqd+p256k0*cvqu
      cauxb=-cfkf78(i5)%ek0*p256(2)+p256k0*cfkf78(i5)%e(2)
      cauxc=+cfkf78(i5)%ek0*p134(2)-p134k0*cfkf78(i5)%e(2)
      u134_78(i5)%a(1)=u134_78(i5)%a(1)+fcr(id1)*(cauxa+ceps_0)
      u134_78(i5)%a(2)=u134_78(i5)%a(2)+fcl(id1)*(cauxa-ceps_0)
      u134_78(i5)%b(1)=u134_78(i5)%b(1)+fcl(id1)*(cauxb-ceps_2)
      u134_78(i5)%b(2)=u134_78(i5)%b(2)+fcr(id1)*(-cauxb-ceps_2)
      u134_78(i5)%c(1)=u134_78(i5)%c(1)+fcr(id1)*(cauxc+ceps_1)
      u134_78(i5)%c(2)=u134_78(i5)%c(2)+fcl(id1)*(-cauxc+ceps_1)
      u134_78(i5)%d(1)=u134_78(i5)%d(1)+fcl(id1)*cfkf78(i5)%ek0
      u134_78(i5)%d(2)=u134_78(i5)%d(2)+fcr(id1)*cfkf78(i5)%ek0
      end do
      endif
  
* quqd -- p=p156,q=p278
      quqd=p156(0)*p278(0)-p156(1)*p278(1)-p156(2)*p278(2)-p156(
     & 3)*p278(3)
      do i5=1,2
* T0 -- qu=p156,qd=p278,v=cfkz34(i5)%e,a=u156_34(i5)%a,b=u156_34(i5)%b,c
* =u156_34(i5)%c,d=u156_34(i5)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      ceps_0=-cfkz34(i5)%ek0*(p156(2)*p278(3)-p278(2)*p156(3))+p
     & 156k0*(cfkz34(i5)%e(2)*p278(3)-p278(2)*cfkz34(i5)%e(3))-p
     & 278k0*(cfkz34(i5)%e(2)*p156(3)-p156(2)*cfkz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz34(i5)%e(3)*p156k0+p156(3)*cfkz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz34(i5)%e(3)*p278k0+p278(3)*cfkz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz34(i5)%e(0)*p156(0)-cfkz34(i5)%e(1)*p156(1)-cfkz3
     & 4(i5)%e(2)*p156(2)-cfkz34(i5)%e(3)*p156(3)
      cvqd=cfkz34(i5)%e(0)*p278(0)-cfkz34(i5)%e(1)*p278(1)-cfkz3
     & 4(i5)%e(2)*p278(2)-cfkz34(i5)%e(3)*p278(3)
      cauxa=-cfkz34(i5)%ek0*quqd+p156k0*cvqd+p278k0*cvqu
      cauxb=-cfkz34(i5)%ek0*p278(2)+p278k0*cfkz34(i5)%e(2)
      cauxc=+cfkz34(i5)%ek0*p156(2)-p156k0*cfkz34(i5)%e(2)
      u156_34(i5)%a(1)=zcr(id1)*(cauxa+ceps_0)
      u156_34(i5)%a(2)=zcl(id1)*(cauxa-ceps_0)
      u156_34(i5)%b(1)=zcl(id1)*(cauxb-ceps_2)
      u156_34(i5)%b(2)=zcr(id1)*(-cauxb-ceps_2)
      u156_34(i5)%c(1)=zcr(id1)*(cauxc+ceps_1)
      u156_34(i5)%c(2)=zcl(id1)*(-cauxc+ceps_1)
      u156_34(i5)%d(1)=zcl(id1)*cfkz34(i5)%ek0
      u156_34(i5)%d(2)=zcr(id1)*cfkz34(i5)%ek0
      end do
  
      if (ineutri(id1).ne.1.and.ineutri(id3).ne.1) then
      do i5=1,2
* T0 -- qu=p156,qd=p278,v=cfkf34(i5)%e,a=u156_34(i5)%a,b=u156_34(i5)%b,c
* =u156_34(i5)%c,d=u156_34(i5)%d,cr=fcr(id1),cl=fcl(id1),nsum=1
      ceps_0=-cfkf34(i5)%ek0*(p156(2)*p278(3)-p278(2)*p156(3))+p
     & 156k0*(cfkf34(i5)%e(2)*p278(3)-p278(2)*cfkf34(i5)%e(3))-p
     & 278k0*(cfkf34(i5)%e(2)*p156(3)-p156(2)*cfkf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf34(i5)%e(3)*p156k0+p156(3)*cfkf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf34(i5)%e(3)*p278k0+p278(3)*cfkf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf34(i5)%e(0)*p156(0)-cfkf34(i5)%e(1)*p156(1)-cfkf3
     & 4(i5)%e(2)*p156(2)-cfkf34(i5)%e(3)*p156(3)
      cvqd=cfkf34(i5)%e(0)*p278(0)-cfkf34(i5)%e(1)*p278(1)-cfkf3
     & 4(i5)%e(2)*p278(2)-cfkf34(i5)%e(3)*p278(3)
      cauxa=-cfkf34(i5)%ek0*quqd+p156k0*cvqd+p278k0*cvqu
      cauxb=-cfkf34(i5)%ek0*p278(2)+p278k0*cfkf34(i5)%e(2)
      cauxc=+cfkf34(i5)%ek0*p156(2)-p156k0*cfkf34(i5)%e(2)
      u156_34(i5)%a(1)=u156_34(i5)%a(1)+fcr(id1)*(cauxa+ceps_0)
      u156_34(i5)%a(2)=u156_34(i5)%a(2)+fcl(id1)*(cauxa-ceps_0)
      u156_34(i5)%b(1)=u156_34(i5)%b(1)+fcl(id1)*(cauxb-ceps_2)
      u156_34(i5)%b(2)=u156_34(i5)%b(2)+fcr(id1)*(-cauxb-ceps_2)
      u156_34(i5)%c(1)=u156_34(i5)%c(1)+fcr(id1)*(cauxc+ceps_1)
      u156_34(i5)%c(2)=u156_34(i5)%c(2)+fcl(id1)*(-cauxc+ceps_1)
      u156_34(i5)%d(1)=u156_34(i5)%d(1)+fcl(id1)*cfkf34(i5)%ek0
      u156_34(i5)%d(2)=u156_34(i5)%d(2)+fcr(id1)*cfkf34(i5)%ek0
      end do
      endif
  
* quqd -- p=p156,q=p234
      quqd=p156(0)*p234(0)-p156(1)*p234(1)-p156(2)*p234(2)-p156(
     & 3)*p234(3)
      do i5=1,2
* T0 -- qu=p156,qd=p234,v=cfkz78(i5)%e,a=u156_78(i5)%a,b=u156_78(i5)%b,c
* =u156_78(i5)%c,d=u156_78(i5)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      ceps_0=-cfkz78(i5)%ek0*(p156(2)*p234(3)-p234(2)*p156(3))+p
     & 156k0*(cfkz78(i5)%e(2)*p234(3)-p234(2)*cfkz78(i5)%e(3))-p
     & 234k0*(cfkz78(i5)%e(2)*p156(3)-p156(2)*cfkz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz78(i5)%e(3)*p156k0+p156(3)*cfkz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz78(i5)%e(3)*p234k0+p234(3)*cfkz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz78(i5)%e(0)*p156(0)-cfkz78(i5)%e(1)*p156(1)-cfkz7
     & 8(i5)%e(2)*p156(2)-cfkz78(i5)%e(3)*p156(3)
      cvqd=cfkz78(i5)%e(0)*p234(0)-cfkz78(i5)%e(1)*p234(1)-cfkz7
     & 8(i5)%e(2)*p234(2)-cfkz78(i5)%e(3)*p234(3)
      cauxa=-cfkz78(i5)%ek0*quqd+p156k0*cvqd+p234k0*cvqu
      cauxb=-cfkz78(i5)%ek0*p234(2)+p234k0*cfkz78(i5)%e(2)
      cauxc=+cfkz78(i5)%ek0*p156(2)-p156k0*cfkz78(i5)%e(2)
      u156_78(i5)%a(1)=zcr(id1)*(cauxa+ceps_0)
      u156_78(i5)%a(2)=zcl(id1)*(cauxa-ceps_0)
      u156_78(i5)%b(1)=zcl(id1)*(cauxb-ceps_2)
      u156_78(i5)%b(2)=zcr(id1)*(-cauxb-ceps_2)
      u156_78(i5)%c(1)=zcr(id1)*(cauxc+ceps_1)
      u156_78(i5)%c(2)=zcl(id1)*(-cauxc+ceps_1)
      u156_78(i5)%d(1)=zcl(id1)*cfkz78(i5)%ek0
      u156_78(i5)%d(2)=zcr(id1)*cfkz78(i5)%ek0
      end do
  
      if (ineutri(id1).ne.1.and.ineutri(id7).ne.1) then
      do i5=1,2
* T0 -- qu=p156,qd=p234,v=cfkf78(i5)%e,a=u156_78(i5)%a,b=u156_78(i5)%b,c
* =u156_78(i5)%c,d=u156_78(i5)%d,cr=fcr(id1),cl=fcl(id1),nsum=1
      ceps_0=-cfkf78(i5)%ek0*(p156(2)*p234(3)-p234(2)*p156(3))+p
     & 156k0*(cfkf78(i5)%e(2)*p234(3)-p234(2)*cfkf78(i5)%e(3))-p
     & 234k0*(cfkf78(i5)%e(2)*p156(3)-p156(2)*cfkf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf78(i5)%e(3)*p156k0+p156(3)*cfkf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf78(i5)%e(3)*p234k0+p234(3)*cfkf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf78(i5)%e(0)*p156(0)-cfkf78(i5)%e(1)*p156(1)-cfkf7
     & 8(i5)%e(2)*p156(2)-cfkf78(i5)%e(3)*p156(3)
      cvqd=cfkf78(i5)%e(0)*p234(0)-cfkf78(i5)%e(1)*p234(1)-cfkf7
     & 8(i5)%e(2)*p234(2)-cfkf78(i5)%e(3)*p234(3)
      cauxa=-cfkf78(i5)%ek0*quqd+p156k0*cvqd+p234k0*cvqu
      cauxb=-cfkf78(i5)%ek0*p234(2)+p234k0*cfkf78(i5)%e(2)
      cauxc=+cfkf78(i5)%ek0*p156(2)-p156k0*cfkf78(i5)%e(2)
      u156_78(i5)%a(1)=u156_78(i5)%a(1)+fcr(id1)*(cauxa+ceps_0)
      u156_78(i5)%a(2)=u156_78(i5)%a(2)+fcl(id1)*(cauxa-ceps_0)
      u156_78(i5)%b(1)=u156_78(i5)%b(1)+fcl(id1)*(cauxb-ceps_2)
      u156_78(i5)%b(2)=u156_78(i5)%b(2)+fcr(id1)*(-cauxb-ceps_2)
      u156_78(i5)%c(1)=u156_78(i5)%c(1)+fcr(id1)*(cauxc+ceps_1)
      u156_78(i5)%c(2)=u156_78(i5)%c(2)+fcl(id1)*(-cauxc+ceps_1)
      u156_78(i5)%d(1)=u156_78(i5)%d(1)+fcl(id1)*cfkf78(i5)%ek0
      u156_78(i5)%d(2)=u156_78(i5)%d(2)+fcr(id1)*cfkf78(i5)%ek0
      end do
      endif
  
* quqd -- p=p178,q=p256
      quqd=p178(0)*p256(0)-p178(1)*p256(1)-p178(2)*p256(2)-p178(
     & 3)*p256(3)
      do i5=1,2
* T0 -- qu=p178,qd=p256,v=cfkz34(i5)%e,a=u178_34(i5)%a,b=u178_34(i5)%b,c
* =u178_34(i5)%c,d=u178_34(i5)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      ceps_0=-cfkz34(i5)%ek0*(p178(2)*p256(3)-p256(2)*p178(3))+p
     & 178k0*(cfkz34(i5)%e(2)*p256(3)-p256(2)*cfkz34(i5)%e(3))-p
     & 256k0*(cfkz34(i5)%e(2)*p178(3)-p178(2)*cfkz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz34(i5)%e(3)*p178k0+p178(3)*cfkz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz34(i5)%e(3)*p256k0+p256(3)*cfkz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz34(i5)%e(0)*p178(0)-cfkz34(i5)%e(1)*p178(1)-cfkz3
     & 4(i5)%e(2)*p178(2)-cfkz34(i5)%e(3)*p178(3)
      cvqd=cfkz34(i5)%e(0)*p256(0)-cfkz34(i5)%e(1)*p256(1)-cfkz3
     & 4(i5)%e(2)*p256(2)-cfkz34(i5)%e(3)*p256(3)
      cauxa=-cfkz34(i5)%ek0*quqd+p178k0*cvqd+p256k0*cvqu
      cauxb=-cfkz34(i5)%ek0*p256(2)+p256k0*cfkz34(i5)%e(2)
      cauxc=+cfkz34(i5)%ek0*p178(2)-p178k0*cfkz34(i5)%e(2)
      u178_34(i5)%a(1)=zcr(id1)*(cauxa+ceps_0)
      u178_34(i5)%a(2)=zcl(id1)*(cauxa-ceps_0)
      u178_34(i5)%b(1)=zcl(id1)*(cauxb-ceps_2)
      u178_34(i5)%b(2)=zcr(id1)*(-cauxb-ceps_2)
      u178_34(i5)%c(1)=zcr(id1)*(cauxc+ceps_1)
      u178_34(i5)%c(2)=zcl(id1)*(-cauxc+ceps_1)
      u178_34(i5)%d(1)=zcl(id1)*cfkz34(i5)%ek0
      u178_34(i5)%d(2)=zcr(id1)*cfkz34(i5)%ek0
      end do
  
      if (ineutri(id1).ne.1.and.ineutri(id3).ne.1) then
      do i5=1,2
* T0 -- qu=p178,qd=p256,v=cfkf34(i5)%e,a=u178_34(i5)%a,b=u178_34(i5)%b,c
* =u178_34(i5)%c,d=u178_34(i5)%d,cr=fcr(id1),cl=fcl(id1),nsum=1
      ceps_0=-cfkf34(i5)%ek0*(p178(2)*p256(3)-p256(2)*p178(3))+p
     & 178k0*(cfkf34(i5)%e(2)*p256(3)-p256(2)*cfkf34(i5)%e(3))-p
     & 256k0*(cfkf34(i5)%e(2)*p178(3)-p178(2)*cfkf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf34(i5)%e(3)*p178k0+p178(3)*cfkf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf34(i5)%e(3)*p256k0+p256(3)*cfkf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf34(i5)%e(0)*p178(0)-cfkf34(i5)%e(1)*p178(1)-cfkf3
     & 4(i5)%e(2)*p178(2)-cfkf34(i5)%e(3)*p178(3)
      cvqd=cfkf34(i5)%e(0)*p256(0)-cfkf34(i5)%e(1)*p256(1)-cfkf3
     & 4(i5)%e(2)*p256(2)-cfkf34(i5)%e(3)*p256(3)
      cauxa=-cfkf34(i5)%ek0*quqd+p178k0*cvqd+p256k0*cvqu
      cauxb=-cfkf34(i5)%ek0*p256(2)+p256k0*cfkf34(i5)%e(2)
      cauxc=+cfkf34(i5)%ek0*p178(2)-p178k0*cfkf34(i5)%e(2)
      u178_34(i5)%a(1)=u178_34(i5)%a(1)+fcr(id1)*(cauxa+ceps_0)
      u178_34(i5)%a(2)=u178_34(i5)%a(2)+fcl(id1)*(cauxa-ceps_0)
      u178_34(i5)%b(1)=u178_34(i5)%b(1)+fcl(id1)*(cauxb-ceps_2)
      u178_34(i5)%b(2)=u178_34(i5)%b(2)+fcr(id1)*(-cauxb-ceps_2)
      u178_34(i5)%c(1)=u178_34(i5)%c(1)+fcr(id1)*(cauxc+ceps_1)
      u178_34(i5)%c(2)=u178_34(i5)%c(2)+fcl(id1)*(-cauxc+ceps_1)
      u178_34(i5)%d(1)=u178_34(i5)%d(1)+fcl(id1)*cfkf34(i5)%ek0
      u178_34(i5)%d(2)=u178_34(i5)%d(2)+fcr(id1)*cfkf34(i5)%ek0
      end do
      endif
  
* quqd -- p=p178,q=p234
      quqd=p178(0)*p234(0)-p178(1)*p234(1)-p178(2)*p234(2)-p178(
     & 3)*p234(3)
      do i5=1,2
* T0 -- qu=p178,qd=p234,v=cfkz56(i5)%e,a=u178_56(i5)%a,b=u178_56(i5)%b,c
* =u178_56(i5)%c,d=u178_56(i5)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      ceps_0=-cfkz56(i5)%ek0*(p178(2)*p234(3)-p234(2)*p178(3))+p
     & 178k0*(cfkz56(i5)%e(2)*p234(3)-p234(2)*cfkz56(i5)%e(3))-p
     & 234k0*(cfkz56(i5)%e(2)*p178(3)-p178(2)*cfkz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz56(i5)%e(3)*p178k0+p178(3)*cfkz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz56(i5)%e(3)*p234k0+p234(3)*cfkz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz56(i5)%e(0)*p178(0)-cfkz56(i5)%e(1)*p178(1)-cfkz5
     & 6(i5)%e(2)*p178(2)-cfkz56(i5)%e(3)*p178(3)
      cvqd=cfkz56(i5)%e(0)*p234(0)-cfkz56(i5)%e(1)*p234(1)-cfkz5
     & 6(i5)%e(2)*p234(2)-cfkz56(i5)%e(3)*p234(3)
      cauxa=-cfkz56(i5)%ek0*quqd+p178k0*cvqd+p234k0*cvqu
      cauxb=-cfkz56(i5)%ek0*p234(2)+p234k0*cfkz56(i5)%e(2)
      cauxc=+cfkz56(i5)%ek0*p178(2)-p178k0*cfkz56(i5)%e(2)
      u178_56(i5)%a(1)=zcr(id1)*(cauxa+ceps_0)
      u178_56(i5)%a(2)=zcl(id1)*(cauxa-ceps_0)
      u178_56(i5)%b(1)=zcl(id1)*(cauxb-ceps_2)
      u178_56(i5)%b(2)=zcr(id1)*(-cauxb-ceps_2)
      u178_56(i5)%c(1)=zcr(id1)*(cauxc+ceps_1)
      u178_56(i5)%c(2)=zcl(id1)*(-cauxc+ceps_1)
      u178_56(i5)%d(1)=zcl(id1)*cfkz56(i5)%ek0
      u178_56(i5)%d(2)=zcr(id1)*cfkz56(i5)%ek0
      end do
  
      if (ineutri(id1).ne.1.and.ineutri(id5).ne.1) then
      do i5=1,2
* T0 -- qu=p178,qd=p234,v=cfkf56(i5)%e,a=u178_56(i5)%a,b=u178_56(i5)%b,c
* =u178_56(i5)%c,d=u178_56(i5)%d,cr=fcr(id1),cl=fcl(id1),nsum=1
      ceps_0=-cfkf56(i5)%ek0*(p178(2)*p234(3)-p234(2)*p178(3))+p
     & 178k0*(cfkf56(i5)%e(2)*p234(3)-p234(2)*cfkf56(i5)%e(3))-p
     & 234k0*(cfkf56(i5)%e(2)*p178(3)-p178(2)*cfkf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf56(i5)%e(3)*p178k0+p178(3)*cfkf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf56(i5)%e(3)*p234k0+p234(3)*cfkf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf56(i5)%e(0)*p178(0)-cfkf56(i5)%e(1)*p178(1)-cfkf5
     & 6(i5)%e(2)*p178(2)-cfkf56(i5)%e(3)*p178(3)
      cvqd=cfkf56(i5)%e(0)*p234(0)-cfkf56(i5)%e(1)*p234(1)-cfkf5
     & 6(i5)%e(2)*p234(2)-cfkf56(i5)%e(3)*p234(3)
      cauxa=-cfkf56(i5)%ek0*quqd+p178k0*cvqd+p234k0*cvqu
      cauxb=-cfkf56(i5)%ek0*p234(2)+p234k0*cfkf56(i5)%e(2)
      cauxc=+cfkf56(i5)%ek0*p178(2)-p178k0*cfkf56(i5)%e(2)
      u178_56(i5)%a(1)=u178_56(i5)%a(1)+fcr(id1)*(cauxa+ceps_0)
      u178_56(i5)%a(2)=u178_56(i5)%a(2)+fcl(id1)*(cauxa-ceps_0)
      u178_56(i5)%b(1)=u178_56(i5)%b(1)+fcl(id1)*(cauxb-ceps_2)
      u178_56(i5)%b(2)=u178_56(i5)%b(2)+fcr(id1)*(-cauxb-ceps_2)
      u178_56(i5)%c(1)=u178_56(i5)%c(1)+fcr(id1)*(cauxc+ceps_1)
      u178_56(i5)%c(2)=u178_56(i5)%c(2)+fcl(id1)*(-cauxc+ceps_1)
      u178_56(i5)%d(1)=u178_56(i5)%d(1)+fcl(id1)*cfkf56(i5)%ek0
      u178_56(i5)%d(2)=u178_56(i5)%d(2)+fcr(id1)*cfkf56(i5)%ek0
      end do
      endif
  
* quqd -- p=p312,q=p478
      quqd=p312(0)*p478(0)-p312(1)*p478(1)-p312(2)*p478(2)-p312(
     & 3)*p478(3)
      do i5=1,2
* T0 -- qu=p312,qd=p478,v=cfkz56(i5)%e,a=u312_56(i5)%a,b=u312_56(i5)%b,c
* =u312_56(i5)%c,d=u312_56(i5)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cfkz56(i5)%ek0*(p312(2)*p478(3)-p478(2)*p312(3))+p
     & 312k0*(cfkz56(i5)%e(2)*p478(3)-p478(2)*cfkz56(i5)%e(3))-p
     & 478k0*(cfkz56(i5)%e(2)*p312(3)-p312(2)*cfkz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz56(i5)%e(3)*p312k0+p312(3)*cfkz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz56(i5)%e(3)*p478k0+p478(3)*cfkz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz56(i5)%e(0)*p312(0)-cfkz56(i5)%e(1)*p312(1)-cfkz5
     & 6(i5)%e(2)*p312(2)-cfkz56(i5)%e(3)*p312(3)
      cvqd=cfkz56(i5)%e(0)*p478(0)-cfkz56(i5)%e(1)*p478(1)-cfkz5
     & 6(i5)%e(2)*p478(2)-cfkz56(i5)%e(3)*p478(3)
      cauxa=-cfkz56(i5)%ek0*quqd+p312k0*cvqd+p478k0*cvqu
      cauxb=-cfkz56(i5)%ek0*p478(2)+p478k0*cfkz56(i5)%e(2)
      cauxc=+cfkz56(i5)%ek0*p312(2)-p312k0*cfkz56(i5)%e(2)
      u312_56(i5)%a(1)=zcr(id3)*(cauxa+ceps_0)
      u312_56(i5)%a(2)=zcl(id3)*(cauxa-ceps_0)
      u312_56(i5)%b(1)=zcl(id3)*(cauxb-ceps_2)
      u312_56(i5)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      u312_56(i5)%c(1)=zcr(id3)*(cauxc+ceps_1)
      u312_56(i5)%c(2)=zcl(id3)*(-cauxc+ceps_1)
      u312_56(i5)%d(1)=zcl(id3)*cfkz56(i5)%ek0
      u312_56(i5)%d(2)=zcr(id3)*cfkz56(i5)%ek0
      end do
  
      if (ineutri(id3).ne.1.and.ineutri(id5).ne.1) then
      do i5=1,2
* T0 -- qu=p312,qd=p478,v=cfkf56(i5)%e,a=u312_56(i5)%a,b=u312_56(i5)%b,c
* =u312_56(i5)%c,d=u312_56(i5)%d,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cfkf56(i5)%ek0*(p312(2)*p478(3)-p478(2)*p312(3))+p
     & 312k0*(cfkf56(i5)%e(2)*p478(3)-p478(2)*cfkf56(i5)%e(3))-p
     & 478k0*(cfkf56(i5)%e(2)*p312(3)-p312(2)*cfkf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf56(i5)%e(3)*p312k0+p312(3)*cfkf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf56(i5)%e(3)*p478k0+p478(3)*cfkf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf56(i5)%e(0)*p312(0)-cfkf56(i5)%e(1)*p312(1)-cfkf5
     & 6(i5)%e(2)*p312(2)-cfkf56(i5)%e(3)*p312(3)
      cvqd=cfkf56(i5)%e(0)*p478(0)-cfkf56(i5)%e(1)*p478(1)-cfkf5
     & 6(i5)%e(2)*p478(2)-cfkf56(i5)%e(3)*p478(3)
      cauxa=-cfkf56(i5)%ek0*quqd+p312k0*cvqd+p478k0*cvqu
      cauxb=-cfkf56(i5)%ek0*p478(2)+p478k0*cfkf56(i5)%e(2)
      cauxc=+cfkf56(i5)%ek0*p312(2)-p312k0*cfkf56(i5)%e(2)
      u312_56(i5)%a(1)=u312_56(i5)%a(1)+fcr(id3)*(cauxa+ceps_0)
      u312_56(i5)%a(2)=u312_56(i5)%a(2)+fcl(id3)*(cauxa-ceps_0)
      u312_56(i5)%b(1)=u312_56(i5)%b(1)+fcl(id3)*(cauxb-ceps_2)
      u312_56(i5)%b(2)=u312_56(i5)%b(2)+fcr(id3)*(-cauxb-ceps_2)
      u312_56(i5)%c(1)=u312_56(i5)%c(1)+fcr(id3)*(cauxc+ceps_1)
      u312_56(i5)%c(2)=u312_56(i5)%c(2)+fcl(id3)*(-cauxc+ceps_1)
      u312_56(i5)%d(1)=u312_56(i5)%d(1)+fcl(id3)*cfkf56(i5)%ek0
      u312_56(i5)%d(2)=u312_56(i5)%d(2)+fcr(id3)*cfkf56(i5)%ek0
      end do
      endif
  
* quqd -- p=p312,q=p456
      quqd=p312(0)*p456(0)-p312(1)*p456(1)-p312(2)*p456(2)-p312(
     & 3)*p456(3)
      do i5=1,2
* T0 -- qu=p312,qd=p456,v=cfkz78(i5)%e,a=u312_78(i5)%a,b=u312_78(i5)%b,c
* =u312_78(i5)%c,d=u312_78(i5)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cfkz78(i5)%ek0*(p312(2)*p456(3)-p456(2)*p312(3))+p
     & 312k0*(cfkz78(i5)%e(2)*p456(3)-p456(2)*cfkz78(i5)%e(3))-p
     & 456k0*(cfkz78(i5)%e(2)*p312(3)-p312(2)*cfkz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz78(i5)%e(3)*p312k0+p312(3)*cfkz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz78(i5)%e(3)*p456k0+p456(3)*cfkz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz78(i5)%e(0)*p312(0)-cfkz78(i5)%e(1)*p312(1)-cfkz7
     & 8(i5)%e(2)*p312(2)-cfkz78(i5)%e(3)*p312(3)
      cvqd=cfkz78(i5)%e(0)*p456(0)-cfkz78(i5)%e(1)*p456(1)-cfkz7
     & 8(i5)%e(2)*p456(2)-cfkz78(i5)%e(3)*p456(3)
      cauxa=-cfkz78(i5)%ek0*quqd+p312k0*cvqd+p456k0*cvqu
      cauxb=-cfkz78(i5)%ek0*p456(2)+p456k0*cfkz78(i5)%e(2)
      cauxc=+cfkz78(i5)%ek0*p312(2)-p312k0*cfkz78(i5)%e(2)
      u312_78(i5)%a(1)=zcr(id3)*(cauxa+ceps_0)
      u312_78(i5)%a(2)=zcl(id3)*(cauxa-ceps_0)
      u312_78(i5)%b(1)=zcl(id3)*(cauxb-ceps_2)
      u312_78(i5)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      u312_78(i5)%c(1)=zcr(id3)*(cauxc+ceps_1)
      u312_78(i5)%c(2)=zcl(id3)*(-cauxc+ceps_1)
      u312_78(i5)%d(1)=zcl(id3)*cfkz78(i5)%ek0
      u312_78(i5)%d(2)=zcr(id3)*cfkz78(i5)%ek0
      end do
  
      if (ineutri(id3).ne.1.and.ineutri(id7).ne.1) then
      do i5=1,2
* T0 -- qu=p312,qd=p456,v=cfkf78(i5)%e,a=u312_78(i5)%a,b=u312_78(i5)%b,c
* =u312_78(i5)%c,d=u312_78(i5)%d,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cfkf78(i5)%ek0*(p312(2)*p456(3)-p456(2)*p312(3))+p
     & 312k0*(cfkf78(i5)%e(2)*p456(3)-p456(2)*cfkf78(i5)%e(3))-p
     & 456k0*(cfkf78(i5)%e(2)*p312(3)-p312(2)*cfkf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf78(i5)%e(3)*p312k0+p312(3)*cfkf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf78(i5)%e(3)*p456k0+p456(3)*cfkf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf78(i5)%e(0)*p312(0)-cfkf78(i5)%e(1)*p312(1)-cfkf7
     & 8(i5)%e(2)*p312(2)-cfkf78(i5)%e(3)*p312(3)
      cvqd=cfkf78(i5)%e(0)*p456(0)-cfkf78(i5)%e(1)*p456(1)-cfkf7
     & 8(i5)%e(2)*p456(2)-cfkf78(i5)%e(3)*p456(3)
      cauxa=-cfkf78(i5)%ek0*quqd+p312k0*cvqd+p456k0*cvqu
      cauxb=-cfkf78(i5)%ek0*p456(2)+p456k0*cfkf78(i5)%e(2)
      cauxc=+cfkf78(i5)%ek0*p312(2)-p312k0*cfkf78(i5)%e(2)
      u312_78(i5)%a(1)=u312_78(i5)%a(1)+fcr(id3)*(cauxa+ceps_0)
      u312_78(i5)%a(2)=u312_78(i5)%a(2)+fcl(id3)*(cauxa-ceps_0)
      u312_78(i5)%b(1)=u312_78(i5)%b(1)+fcl(id3)*(cauxb-ceps_2)
      u312_78(i5)%b(2)=u312_78(i5)%b(2)+fcr(id3)*(-cauxb-ceps_2)
      u312_78(i5)%c(1)=u312_78(i5)%c(1)+fcr(id3)*(cauxc+ceps_1)
      u312_78(i5)%c(2)=u312_78(i5)%c(2)+fcl(id3)*(-cauxc+ceps_1)
      u312_78(i5)%d(1)=u312_78(i5)%d(1)+fcl(id3)*cfkf78(i5)%ek0
      u312_78(i5)%d(2)=u312_78(i5)%d(2)+fcr(id3)*cfkf78(i5)%ek0
      end do
      endif
  
* quqd -- p=p356,q=p478
      quqd=p356(0)*p478(0)-p356(1)*p478(1)-p356(2)*p478(2)-p356(
     & 3)*p478(3)
      do i5=1,2
* T0 -- qu=p356,qd=p478,v=cfkz12(i5)%e,a=u356_12(i5)%a,b=u356_12(i5)%b,c
* =u356_12(i5)%c,d=u356_12(i5)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cfkz12(i5)%ek0*(p356(2)*p478(3)-p478(2)*p356(3))+p
     & 356k0*(cfkz12(i5)%e(2)*p478(3)-p478(2)*cfkz12(i5)%e(3))-p
     & 478k0*(cfkz12(i5)%e(2)*p356(3)-p356(2)*cfkz12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz12(i5)%e(3)*p356k0+p356(3)*cfkz12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz12(i5)%e(3)*p478k0+p478(3)*cfkz12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz12(i5)%e(0)*p356(0)-cfkz12(i5)%e(1)*p356(1)-cfkz1
     & 2(i5)%e(2)*p356(2)-cfkz12(i5)%e(3)*p356(3)
      cvqd=cfkz12(i5)%e(0)*p478(0)-cfkz12(i5)%e(1)*p478(1)-cfkz1
     & 2(i5)%e(2)*p478(2)-cfkz12(i5)%e(3)*p478(3)
      cauxa=-cfkz12(i5)%ek0*quqd+p356k0*cvqd+p478k0*cvqu
      cauxb=-cfkz12(i5)%ek0*p478(2)+p478k0*cfkz12(i5)%e(2)
      cauxc=+cfkz12(i5)%ek0*p356(2)-p356k0*cfkz12(i5)%e(2)
      u356_12(i5)%a(1)=zcr(id3)*(cauxa+ceps_0)
      u356_12(i5)%a(2)=zcl(id3)*(cauxa-ceps_0)
      u356_12(i5)%b(1)=zcl(id3)*(cauxb-ceps_2)
      u356_12(i5)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      u356_12(i5)%c(1)=zcr(id3)*(cauxc+ceps_1)
      u356_12(i5)%c(2)=zcl(id3)*(-cauxc+ceps_1)
      u356_12(i5)%d(1)=zcl(id3)*cfkz12(i5)%ek0
      u356_12(i5)%d(2)=zcr(id3)*cfkz12(i5)%ek0
      end do
  
      if (ineutri(id3).ne.1.and.ineutri(id1).ne.1) then
      do i5=1,2
* T0 -- qu=p356,qd=p478,v=cfkf12(i5)%e,a=u356_12(i5)%a,b=u356_12(i5)%b,c
* =u356_12(i5)%c,d=u356_12(i5)%d,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cfkf12(i5)%ek0*(p356(2)*p478(3)-p478(2)*p356(3))+p
     & 356k0*(cfkf12(i5)%e(2)*p478(3)-p478(2)*cfkf12(i5)%e(3))-p
     & 478k0*(cfkf12(i5)%e(2)*p356(3)-p356(2)*cfkf12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf12(i5)%e(3)*p356k0+p356(3)*cfkf12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf12(i5)%e(3)*p478k0+p478(3)*cfkf12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf12(i5)%e(0)*p356(0)-cfkf12(i5)%e(1)*p356(1)-cfkf1
     & 2(i5)%e(2)*p356(2)-cfkf12(i5)%e(3)*p356(3)
      cvqd=cfkf12(i5)%e(0)*p478(0)-cfkf12(i5)%e(1)*p478(1)-cfkf1
     & 2(i5)%e(2)*p478(2)-cfkf12(i5)%e(3)*p478(3)
      cauxa=-cfkf12(i5)%ek0*quqd+p356k0*cvqd+p478k0*cvqu
      cauxb=-cfkf12(i5)%ek0*p478(2)+p478k0*cfkf12(i5)%e(2)
      cauxc=+cfkf12(i5)%ek0*p356(2)-p356k0*cfkf12(i5)%e(2)
      u356_12(i5)%a(1)=u356_12(i5)%a(1)+fcr(id3)*(cauxa+ceps_0)
      u356_12(i5)%a(2)=u356_12(i5)%a(2)+fcl(id3)*(cauxa-ceps_0)
      u356_12(i5)%b(1)=u356_12(i5)%b(1)+fcl(id3)*(cauxb-ceps_2)
      u356_12(i5)%b(2)=u356_12(i5)%b(2)+fcr(id3)*(-cauxb-ceps_2)
      u356_12(i5)%c(1)=u356_12(i5)%c(1)+fcr(id3)*(cauxc+ceps_1)
      u356_12(i5)%c(2)=u356_12(i5)%c(2)+fcl(id3)*(-cauxc+ceps_1)
      u356_12(i5)%d(1)=u356_12(i5)%d(1)+fcl(id3)*cfkf12(i5)%ek0
      u356_12(i5)%d(2)=u356_12(i5)%d(2)+fcr(id3)*cfkf12(i5)%ek0
      end do
      endif
  
* quqd -- p=p356,q=p412
      quqd=p356(0)*p412(0)-p356(1)*p412(1)-p356(2)*p412(2)-p356(
     & 3)*p412(3)
      do i5=1,2
* T0 -- qu=p356,qd=p412,v=cfkz78(i5)%e,a=u356_78(i5)%a,b=u356_78(i5)%b,c
* =u356_78(i5)%c,d=u356_78(i5)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cfkz78(i5)%ek0*(p356(2)*p412(3)-p412(2)*p356(3))+p
     & 356k0*(cfkz78(i5)%e(2)*p412(3)-p412(2)*cfkz78(i5)%e(3))-p
     & 412k0*(cfkz78(i5)%e(2)*p356(3)-p356(2)*cfkz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz78(i5)%e(3)*p356k0+p356(3)*cfkz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz78(i5)%e(3)*p412k0+p412(3)*cfkz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz78(i5)%e(0)*p356(0)-cfkz78(i5)%e(1)*p356(1)-cfkz7
     & 8(i5)%e(2)*p356(2)-cfkz78(i5)%e(3)*p356(3)
      cvqd=cfkz78(i5)%e(0)*p412(0)-cfkz78(i5)%e(1)*p412(1)-cfkz7
     & 8(i5)%e(2)*p412(2)-cfkz78(i5)%e(3)*p412(3)
      cauxa=-cfkz78(i5)%ek0*quqd+p356k0*cvqd+p412k0*cvqu
      cauxb=-cfkz78(i5)%ek0*p412(2)+p412k0*cfkz78(i5)%e(2)
      cauxc=+cfkz78(i5)%ek0*p356(2)-p356k0*cfkz78(i5)%e(2)
      u356_78(i5)%a(1)=zcr(id3)*(cauxa+ceps_0)
      u356_78(i5)%a(2)=zcl(id3)*(cauxa-ceps_0)
      u356_78(i5)%b(1)=zcl(id3)*(cauxb-ceps_2)
      u356_78(i5)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      u356_78(i5)%c(1)=zcr(id3)*(cauxc+ceps_1)
      u356_78(i5)%c(2)=zcl(id3)*(-cauxc+ceps_1)
      u356_78(i5)%d(1)=zcl(id3)*cfkz78(i5)%ek0
      u356_78(i5)%d(2)=zcr(id3)*cfkz78(i5)%ek0
      end do
  
      if (ineutri(id3).ne.1.and.ineutri(id7).ne.1) then
      do i5=1,2
* T0 -- qu=p356,qd=p412,v=cfkf78(i5)%e,a=u356_78(i5)%a,b=u356_78(i5)%b,c
* =u356_78(i5)%c,d=u356_78(i5)%d,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cfkf78(i5)%ek0*(p356(2)*p412(3)-p412(2)*p356(3))+p
     & 356k0*(cfkf78(i5)%e(2)*p412(3)-p412(2)*cfkf78(i5)%e(3))-p
     & 412k0*(cfkf78(i5)%e(2)*p356(3)-p356(2)*cfkf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf78(i5)%e(3)*p356k0+p356(3)*cfkf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf78(i5)%e(3)*p412k0+p412(3)*cfkf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf78(i5)%e(0)*p356(0)-cfkf78(i5)%e(1)*p356(1)-cfkf7
     & 8(i5)%e(2)*p356(2)-cfkf78(i5)%e(3)*p356(3)
      cvqd=cfkf78(i5)%e(0)*p412(0)-cfkf78(i5)%e(1)*p412(1)-cfkf7
     & 8(i5)%e(2)*p412(2)-cfkf78(i5)%e(3)*p412(3)
      cauxa=-cfkf78(i5)%ek0*quqd+p356k0*cvqd+p412k0*cvqu
      cauxb=-cfkf78(i5)%ek0*p412(2)+p412k0*cfkf78(i5)%e(2)
      cauxc=+cfkf78(i5)%ek0*p356(2)-p356k0*cfkf78(i5)%e(2)
      u356_78(i5)%a(1)=u356_78(i5)%a(1)+fcr(id3)*(cauxa+ceps_0)
      u356_78(i5)%a(2)=u356_78(i5)%a(2)+fcl(id3)*(cauxa-ceps_0)
      u356_78(i5)%b(1)=u356_78(i5)%b(1)+fcl(id3)*(cauxb-ceps_2)
      u356_78(i5)%b(2)=u356_78(i5)%b(2)+fcr(id3)*(-cauxb-ceps_2)
      u356_78(i5)%c(1)=u356_78(i5)%c(1)+fcr(id3)*(cauxc+ceps_1)
      u356_78(i5)%c(2)=u356_78(i5)%c(2)+fcl(id3)*(-cauxc+ceps_1)
      u356_78(i5)%d(1)=u356_78(i5)%d(1)+fcl(id3)*cfkf78(i5)%ek0
      u356_78(i5)%d(2)=u356_78(i5)%d(2)+fcr(id3)*cfkf78(i5)%ek0
      end do
      endif
  
* quqd -- p=p378,q=p456
      quqd=p378(0)*p456(0)-p378(1)*p456(1)-p378(2)*p456(2)-p378(
     & 3)*p456(3)
      do i5=1,2
* T0 -- qu=p378,qd=p456,v=cfkz12(i5)%e,a=u378_12(i5)%a,b=u378_12(i5)%b,c
* =u378_12(i5)%c,d=u378_12(i5)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cfkz12(i5)%ek0*(p378(2)*p456(3)-p456(2)*p378(3))+p
     & 378k0*(cfkz12(i5)%e(2)*p456(3)-p456(2)*cfkz12(i5)%e(3))-p
     & 456k0*(cfkz12(i5)%e(2)*p378(3)-p378(2)*cfkz12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz12(i5)%e(3)*p378k0+p378(3)*cfkz12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz12(i5)%e(3)*p456k0+p456(3)*cfkz12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz12(i5)%e(0)*p378(0)-cfkz12(i5)%e(1)*p378(1)-cfkz1
     & 2(i5)%e(2)*p378(2)-cfkz12(i5)%e(3)*p378(3)
      cvqd=cfkz12(i5)%e(0)*p456(0)-cfkz12(i5)%e(1)*p456(1)-cfkz1
     & 2(i5)%e(2)*p456(2)-cfkz12(i5)%e(3)*p456(3)
      cauxa=-cfkz12(i5)%ek0*quqd+p378k0*cvqd+p456k0*cvqu
      cauxb=-cfkz12(i5)%ek0*p456(2)+p456k0*cfkz12(i5)%e(2)
      cauxc=+cfkz12(i5)%ek0*p378(2)-p378k0*cfkz12(i5)%e(2)
      u378_12(i5)%a(1)=zcr(id3)*(cauxa+ceps_0)
      u378_12(i5)%a(2)=zcl(id3)*(cauxa-ceps_0)
      u378_12(i5)%b(1)=zcl(id3)*(cauxb-ceps_2)
      u378_12(i5)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      u378_12(i5)%c(1)=zcr(id3)*(cauxc+ceps_1)
      u378_12(i5)%c(2)=zcl(id3)*(-cauxc+ceps_1)
      u378_12(i5)%d(1)=zcl(id3)*cfkz12(i5)%ek0
      u378_12(i5)%d(2)=zcr(id3)*cfkz12(i5)%ek0
      end do
  
      if (ineutri(id3).ne.1.and.ineutri(id1).ne.1) then
      do i5=1,2
* T0 -- qu=p378,qd=p456,v=cfkf12(i5)%e,a=u378_12(i5)%a,b=u378_12(i5)%b,c
* =u378_12(i5)%c,d=u378_12(i5)%d,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cfkf12(i5)%ek0*(p378(2)*p456(3)-p456(2)*p378(3))+p
     & 378k0*(cfkf12(i5)%e(2)*p456(3)-p456(2)*cfkf12(i5)%e(3))-p
     & 456k0*(cfkf12(i5)%e(2)*p378(3)-p378(2)*cfkf12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf12(i5)%e(3)*p378k0+p378(3)*cfkf12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf12(i5)%e(3)*p456k0+p456(3)*cfkf12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf12(i5)%e(0)*p378(0)-cfkf12(i5)%e(1)*p378(1)-cfkf1
     & 2(i5)%e(2)*p378(2)-cfkf12(i5)%e(3)*p378(3)
      cvqd=cfkf12(i5)%e(0)*p456(0)-cfkf12(i5)%e(1)*p456(1)-cfkf1
     & 2(i5)%e(2)*p456(2)-cfkf12(i5)%e(3)*p456(3)
      cauxa=-cfkf12(i5)%ek0*quqd+p378k0*cvqd+p456k0*cvqu
      cauxb=-cfkf12(i5)%ek0*p456(2)+p456k0*cfkf12(i5)%e(2)
      cauxc=+cfkf12(i5)%ek0*p378(2)-p378k0*cfkf12(i5)%e(2)
      u378_12(i5)%a(1)=u378_12(i5)%a(1)+fcr(id3)*(cauxa+ceps_0)
      u378_12(i5)%a(2)=u378_12(i5)%a(2)+fcl(id3)*(cauxa-ceps_0)
      u378_12(i5)%b(1)=u378_12(i5)%b(1)+fcl(id3)*(cauxb-ceps_2)
      u378_12(i5)%b(2)=u378_12(i5)%b(2)+fcr(id3)*(-cauxb-ceps_2)
      u378_12(i5)%c(1)=u378_12(i5)%c(1)+fcr(id3)*(cauxc+ceps_1)
      u378_12(i5)%c(2)=u378_12(i5)%c(2)+fcl(id3)*(-cauxc+ceps_1)
      u378_12(i5)%d(1)=u378_12(i5)%d(1)+fcl(id3)*cfkf12(i5)%ek0
      u378_12(i5)%d(2)=u378_12(i5)%d(2)+fcr(id3)*cfkf12(i5)%ek0
      end do
      endif
  
* quqd -- p=p378,q=p412
      quqd=p378(0)*p412(0)-p378(1)*p412(1)-p378(2)*p412(2)-p378(
     & 3)*p412(3)
      do i5=1,2
* T0 -- qu=p378,qd=p412,v=cfkz56(i5)%e,a=u378_56(i5)%a,b=u378_56(i5)%b,c
* =u378_56(i5)%c,d=u378_56(i5)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      ceps_0=-cfkz56(i5)%ek0*(p378(2)*p412(3)-p412(2)*p378(3))+p
     & 378k0*(cfkz56(i5)%e(2)*p412(3)-p412(2)*cfkz56(i5)%e(3))-p
     & 412k0*(cfkz56(i5)%e(2)*p378(3)-p378(2)*cfkz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz56(i5)%e(3)*p378k0+p378(3)*cfkz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz56(i5)%e(3)*p412k0+p412(3)*cfkz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz56(i5)%e(0)*p378(0)-cfkz56(i5)%e(1)*p378(1)-cfkz5
     & 6(i5)%e(2)*p378(2)-cfkz56(i5)%e(3)*p378(3)
      cvqd=cfkz56(i5)%e(0)*p412(0)-cfkz56(i5)%e(1)*p412(1)-cfkz5
     & 6(i5)%e(2)*p412(2)-cfkz56(i5)%e(3)*p412(3)
      cauxa=-cfkz56(i5)%ek0*quqd+p378k0*cvqd+p412k0*cvqu
      cauxb=-cfkz56(i5)%ek0*p412(2)+p412k0*cfkz56(i5)%e(2)
      cauxc=+cfkz56(i5)%ek0*p378(2)-p378k0*cfkz56(i5)%e(2)
      u378_56(i5)%a(1)=zcr(id3)*(cauxa+ceps_0)
      u378_56(i5)%a(2)=zcl(id3)*(cauxa-ceps_0)
      u378_56(i5)%b(1)=zcl(id3)*(cauxb-ceps_2)
      u378_56(i5)%b(2)=zcr(id3)*(-cauxb-ceps_2)
      u378_56(i5)%c(1)=zcr(id3)*(cauxc+ceps_1)
      u378_56(i5)%c(2)=zcl(id3)*(-cauxc+ceps_1)
      u378_56(i5)%d(1)=zcl(id3)*cfkz56(i5)%ek0
      u378_56(i5)%d(2)=zcr(id3)*cfkz56(i5)%ek0
      end do
  
      if (ineutri(id3).ne.1.and.ineutri(id5).ne.1) then
      do i5=1,2
* T0 -- qu=p378,qd=p412,v=cfkf56(i5)%e,a=u378_56(i5)%a,b=u378_56(i5)%b,c
* =u378_56(i5)%c,d=u378_56(i5)%d,cr=fcr(id3),cl=fcl(id3),nsum=1
      ceps_0=-cfkf56(i5)%ek0*(p378(2)*p412(3)-p412(2)*p378(3))+p
     & 378k0*(cfkf56(i5)%e(2)*p412(3)-p412(2)*cfkf56(i5)%e(3))-p
     & 412k0*(cfkf56(i5)%e(2)*p378(3)-p378(2)*cfkf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf56(i5)%e(3)*p378k0+p378(3)*cfkf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf56(i5)%e(3)*p412k0+p412(3)*cfkf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf56(i5)%e(0)*p378(0)-cfkf56(i5)%e(1)*p378(1)-cfkf5
     & 6(i5)%e(2)*p378(2)-cfkf56(i5)%e(3)*p378(3)
      cvqd=cfkf56(i5)%e(0)*p412(0)-cfkf56(i5)%e(1)*p412(1)-cfkf5
     & 6(i5)%e(2)*p412(2)-cfkf56(i5)%e(3)*p412(3)
      cauxa=-cfkf56(i5)%ek0*quqd+p378k0*cvqd+p412k0*cvqu
      cauxb=-cfkf56(i5)%ek0*p412(2)+p412k0*cfkf56(i5)%e(2)
      cauxc=+cfkf56(i5)%ek0*p378(2)-p378k0*cfkf56(i5)%e(2)
      u378_56(i5)%a(1)=u378_56(i5)%a(1)+fcr(id3)*(cauxa+ceps_0)
      u378_56(i5)%a(2)=u378_56(i5)%a(2)+fcl(id3)*(cauxa-ceps_0)
      u378_56(i5)%b(1)=u378_56(i5)%b(1)+fcl(id3)*(cauxb-ceps_2)
      u378_56(i5)%b(2)=u378_56(i5)%b(2)+fcr(id3)*(-cauxb-ceps_2)
      u378_56(i5)%c(1)=u378_56(i5)%c(1)+fcr(id3)*(cauxc+ceps_1)
      u378_56(i5)%c(2)=u378_56(i5)%c(2)+fcl(id3)*(-cauxc+ceps_1)
      u378_56(i5)%d(1)=u378_56(i5)%d(1)+fcl(id3)*cfkf56(i5)%ek0
      u378_56(i5)%d(2)=u378_56(i5)%d(2)+fcr(id3)*cfkf56(i5)%ek0
      end do
      endif
  
* quqd -- p=p512,q=p678
      quqd=p512(0)*p678(0)-p512(1)*p678(1)-p512(2)*p678(2)-p512(
     & 3)*p678(3)
      do i5=1,2
* T0 -- qu=p512,qd=p678,v=cfkz34(i5)%e,a=u512_34(i5)%a,b=u512_34(i5)%b,c
* =u512_34(i5)%c,d=u512_34(i5)%d,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cfkz34(i5)%ek0*(p512(2)*p678(3)-p678(2)*p512(3))+p
     & 512k0*(cfkz34(i5)%e(2)*p678(3)-p678(2)*cfkz34(i5)%e(3))-p
     & 678k0*(cfkz34(i5)%e(2)*p512(3)-p512(2)*cfkz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz34(i5)%e(3)*p512k0+p512(3)*cfkz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz34(i5)%e(3)*p678k0+p678(3)*cfkz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz34(i5)%e(0)*p512(0)-cfkz34(i5)%e(1)*p512(1)-cfkz3
     & 4(i5)%e(2)*p512(2)-cfkz34(i5)%e(3)*p512(3)
      cvqd=cfkz34(i5)%e(0)*p678(0)-cfkz34(i5)%e(1)*p678(1)-cfkz3
     & 4(i5)%e(2)*p678(2)-cfkz34(i5)%e(3)*p678(3)
      cauxa=-cfkz34(i5)%ek0*quqd+p512k0*cvqd+p678k0*cvqu
      cauxb=-cfkz34(i5)%ek0*p678(2)+p678k0*cfkz34(i5)%e(2)
      cauxc=+cfkz34(i5)%ek0*p512(2)-p512k0*cfkz34(i5)%e(2)
      u512_34(i5)%a(1)=zcr(id5)*(cauxa+ceps_0)
      u512_34(i5)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u512_34(i5)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u512_34(i5)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      u512_34(i5)%c(1)=zcr(id5)*(cauxc+ceps_1)
      u512_34(i5)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u512_34(i5)%d(1)=zcl(id5)*cfkz34(i5)%ek0
      u512_34(i5)%d(2)=zcr(id5)*cfkz34(i5)%ek0
      end do
  
      if (ineutri(id5).ne.1.and.ineutri(id3).ne.1) then
      do i5=1,2
* T0 -- qu=p512,qd=p678,v=cfkf34(i5)%e,a=u512_34(i5)%a,b=u512_34(i5)%b,c
* =u512_34(i5)%c,d=u512_34(i5)%d,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cfkf34(i5)%ek0*(p512(2)*p678(3)-p678(2)*p512(3))+p
     & 512k0*(cfkf34(i5)%e(2)*p678(3)-p678(2)*cfkf34(i5)%e(3))-p
     & 678k0*(cfkf34(i5)%e(2)*p512(3)-p512(2)*cfkf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf34(i5)%e(3)*p512k0+p512(3)*cfkf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf34(i5)%e(3)*p678k0+p678(3)*cfkf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf34(i5)%e(0)*p512(0)-cfkf34(i5)%e(1)*p512(1)-cfkf3
     & 4(i5)%e(2)*p512(2)-cfkf34(i5)%e(3)*p512(3)
      cvqd=cfkf34(i5)%e(0)*p678(0)-cfkf34(i5)%e(1)*p678(1)-cfkf3
     & 4(i5)%e(2)*p678(2)-cfkf34(i5)%e(3)*p678(3)
      cauxa=-cfkf34(i5)%ek0*quqd+p512k0*cvqd+p678k0*cvqu
      cauxb=-cfkf34(i5)%ek0*p678(2)+p678k0*cfkf34(i5)%e(2)
      cauxc=+cfkf34(i5)%ek0*p512(2)-p512k0*cfkf34(i5)%e(2)
      u512_34(i5)%a(1)=u512_34(i5)%a(1)+fcr(id5)*(cauxa+ceps_0)
      u512_34(i5)%a(2)=u512_34(i5)%a(2)+fcl(id5)*(cauxa-ceps_0)
      u512_34(i5)%b(1)=u512_34(i5)%b(1)+fcl(id5)*(cauxb-ceps_2)
      u512_34(i5)%b(2)=u512_34(i5)%b(2)+fcr(id5)*(-cauxb-ceps_2)
      u512_34(i5)%c(1)=u512_34(i5)%c(1)+fcr(id5)*(cauxc+ceps_1)
      u512_34(i5)%c(2)=u512_34(i5)%c(2)+fcl(id5)*(-cauxc+ceps_1)
      u512_34(i5)%d(1)=u512_34(i5)%d(1)+fcl(id5)*cfkf34(i5)%ek0
      u512_34(i5)%d(2)=u512_34(i5)%d(2)+fcr(id5)*cfkf34(i5)%ek0
      end do
      endif
  
* quqd -- p=p512,q=p634
      quqd=p512(0)*p634(0)-p512(1)*p634(1)-p512(2)*p634(2)-p512(
     & 3)*p634(3)
      do i5=1,2
* T0 -- qu=p512,qd=p634,v=cfkz78(i5)%e,a=u512_78(i5)%a,b=u512_78(i5)%b,c
* =u512_78(i5)%c,d=u512_78(i5)%d,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cfkz78(i5)%ek0*(p512(2)*p634(3)-p634(2)*p512(3))+p
     & 512k0*(cfkz78(i5)%e(2)*p634(3)-p634(2)*cfkz78(i5)%e(3))-p
     & 634k0*(cfkz78(i5)%e(2)*p512(3)-p512(2)*cfkz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz78(i5)%e(3)*p512k0+p512(3)*cfkz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz78(i5)%e(3)*p634k0+p634(3)*cfkz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz78(i5)%e(0)*p512(0)-cfkz78(i5)%e(1)*p512(1)-cfkz7
     & 8(i5)%e(2)*p512(2)-cfkz78(i5)%e(3)*p512(3)
      cvqd=cfkz78(i5)%e(0)*p634(0)-cfkz78(i5)%e(1)*p634(1)-cfkz7
     & 8(i5)%e(2)*p634(2)-cfkz78(i5)%e(3)*p634(3)
      cauxa=-cfkz78(i5)%ek0*quqd+p512k0*cvqd+p634k0*cvqu
      cauxb=-cfkz78(i5)%ek0*p634(2)+p634k0*cfkz78(i5)%e(2)
      cauxc=+cfkz78(i5)%ek0*p512(2)-p512k0*cfkz78(i5)%e(2)
      u512_78(i5)%a(1)=zcr(id5)*(cauxa+ceps_0)
      u512_78(i5)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u512_78(i5)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u512_78(i5)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      u512_78(i5)%c(1)=zcr(id5)*(cauxc+ceps_1)
      u512_78(i5)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u512_78(i5)%d(1)=zcl(id5)*cfkz78(i5)%ek0
      u512_78(i5)%d(2)=zcr(id5)*cfkz78(i5)%ek0
      end do
  
      if (ineutri(id5).ne.1.and.ineutri(id7).ne.1) then
      do i5=1,2
* T0 -- qu=p512,qd=p634,v=cfkf78(i5)%e,a=u512_78(i5)%a,b=u512_78(i5)%b,c
* =u512_78(i5)%c,d=u512_78(i5)%d,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cfkf78(i5)%ek0*(p512(2)*p634(3)-p634(2)*p512(3))+p
     & 512k0*(cfkf78(i5)%e(2)*p634(3)-p634(2)*cfkf78(i5)%e(3))-p
     & 634k0*(cfkf78(i5)%e(2)*p512(3)-p512(2)*cfkf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf78(i5)%e(3)*p512k0+p512(3)*cfkf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf78(i5)%e(3)*p634k0+p634(3)*cfkf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf78(i5)%e(0)*p512(0)-cfkf78(i5)%e(1)*p512(1)-cfkf7
     & 8(i5)%e(2)*p512(2)-cfkf78(i5)%e(3)*p512(3)
      cvqd=cfkf78(i5)%e(0)*p634(0)-cfkf78(i5)%e(1)*p634(1)-cfkf7
     & 8(i5)%e(2)*p634(2)-cfkf78(i5)%e(3)*p634(3)
      cauxa=-cfkf78(i5)%ek0*quqd+p512k0*cvqd+p634k0*cvqu
      cauxb=-cfkf78(i5)%ek0*p634(2)+p634k0*cfkf78(i5)%e(2)
      cauxc=+cfkf78(i5)%ek0*p512(2)-p512k0*cfkf78(i5)%e(2)
      u512_78(i5)%a(1)=u512_78(i5)%a(1)+fcr(id5)*(cauxa+ceps_0)
      u512_78(i5)%a(2)=u512_78(i5)%a(2)+fcl(id5)*(cauxa-ceps_0)
      u512_78(i5)%b(1)=u512_78(i5)%b(1)+fcl(id5)*(cauxb-ceps_2)
      u512_78(i5)%b(2)=u512_78(i5)%b(2)+fcr(id5)*(-cauxb-ceps_2)
      u512_78(i5)%c(1)=u512_78(i5)%c(1)+fcr(id5)*(cauxc+ceps_1)
      u512_78(i5)%c(2)=u512_78(i5)%c(2)+fcl(id5)*(-cauxc+ceps_1)
      u512_78(i5)%d(1)=u512_78(i5)%d(1)+fcl(id5)*cfkf78(i5)%ek0
      u512_78(i5)%d(2)=u512_78(i5)%d(2)+fcr(id5)*cfkf78(i5)%ek0
      end do
      endif
  
* quqd -- p=p534,q=p678
      quqd=p534(0)*p678(0)-p534(1)*p678(1)-p534(2)*p678(2)-p534(
     & 3)*p678(3)
      do i5=1,2
* T0 -- qu=p534,qd=p678,v=cfkz12(i5)%e,a=u534_12(i5)%a,b=u534_12(i5)%b,c
* =u534_12(i5)%c,d=u534_12(i5)%d,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cfkz12(i5)%ek0*(p534(2)*p678(3)-p678(2)*p534(3))+p
     & 534k0*(cfkz12(i5)%e(2)*p678(3)-p678(2)*cfkz12(i5)%e(3))-p
     & 678k0*(cfkz12(i5)%e(2)*p534(3)-p534(2)*cfkz12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz12(i5)%e(3)*p534k0+p534(3)*cfkz12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz12(i5)%e(3)*p678k0+p678(3)*cfkz12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz12(i5)%e(0)*p534(0)-cfkz12(i5)%e(1)*p534(1)-cfkz1
     & 2(i5)%e(2)*p534(2)-cfkz12(i5)%e(3)*p534(3)
      cvqd=cfkz12(i5)%e(0)*p678(0)-cfkz12(i5)%e(1)*p678(1)-cfkz1
     & 2(i5)%e(2)*p678(2)-cfkz12(i5)%e(3)*p678(3)
      cauxa=-cfkz12(i5)%ek0*quqd+p534k0*cvqd+p678k0*cvqu
      cauxb=-cfkz12(i5)%ek0*p678(2)+p678k0*cfkz12(i5)%e(2)
      cauxc=+cfkz12(i5)%ek0*p534(2)-p534k0*cfkz12(i5)%e(2)
      u534_12(i5)%a(1)=zcr(id5)*(cauxa+ceps_0)
      u534_12(i5)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u534_12(i5)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u534_12(i5)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      u534_12(i5)%c(1)=zcr(id5)*(cauxc+ceps_1)
      u534_12(i5)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u534_12(i5)%d(1)=zcl(id5)*cfkz12(i5)%ek0
      u534_12(i5)%d(2)=zcr(id5)*cfkz12(i5)%ek0
      end do
  
      if (ineutri(id5).ne.1.and.ineutri(id1).ne.1) then
      do i5=1,2
* T0 -- qu=p534,qd=p678,v=cfkf12(i5)%e,a=u534_12(i5)%a,b=u534_12(i5)%b,c
* =u534_12(i5)%c,d=u534_12(i5)%d,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cfkf12(i5)%ek0*(p534(2)*p678(3)-p678(2)*p534(3))+p
     & 534k0*(cfkf12(i5)%e(2)*p678(3)-p678(2)*cfkf12(i5)%e(3))-p
     & 678k0*(cfkf12(i5)%e(2)*p534(3)-p534(2)*cfkf12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf12(i5)%e(3)*p534k0+p534(3)*cfkf12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf12(i5)%e(3)*p678k0+p678(3)*cfkf12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf12(i5)%e(0)*p534(0)-cfkf12(i5)%e(1)*p534(1)-cfkf1
     & 2(i5)%e(2)*p534(2)-cfkf12(i5)%e(3)*p534(3)
      cvqd=cfkf12(i5)%e(0)*p678(0)-cfkf12(i5)%e(1)*p678(1)-cfkf1
     & 2(i5)%e(2)*p678(2)-cfkf12(i5)%e(3)*p678(3)
      cauxa=-cfkf12(i5)%ek0*quqd+p534k0*cvqd+p678k0*cvqu
      cauxb=-cfkf12(i5)%ek0*p678(2)+p678k0*cfkf12(i5)%e(2)
      cauxc=+cfkf12(i5)%ek0*p534(2)-p534k0*cfkf12(i5)%e(2)
      u534_12(i5)%a(1)=u534_12(i5)%a(1)+fcr(id5)*(cauxa+ceps_0)
      u534_12(i5)%a(2)=u534_12(i5)%a(2)+fcl(id5)*(cauxa-ceps_0)
      u534_12(i5)%b(1)=u534_12(i5)%b(1)+fcl(id5)*(cauxb-ceps_2)
      u534_12(i5)%b(2)=u534_12(i5)%b(2)+fcr(id5)*(-cauxb-ceps_2)
      u534_12(i5)%c(1)=u534_12(i5)%c(1)+fcr(id5)*(cauxc+ceps_1)
      u534_12(i5)%c(2)=u534_12(i5)%c(2)+fcl(id5)*(-cauxc+ceps_1)
      u534_12(i5)%d(1)=u534_12(i5)%d(1)+fcl(id5)*cfkf12(i5)%ek0
      u534_12(i5)%d(2)=u534_12(i5)%d(2)+fcr(id5)*cfkf12(i5)%ek0
      end do
      endif
  
* quqd -- p=p534,q=p612
      quqd=p534(0)*p612(0)-p534(1)*p612(1)-p534(2)*p612(2)-p534(
     & 3)*p612(3)
      do i5=1,2
* T0 -- qu=p534,qd=p612,v=cfkz78(i5)%e,a=u534_78(i5)%a,b=u534_78(i5)%b,c
* =u534_78(i5)%c,d=u534_78(i5)%d,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cfkz78(i5)%ek0*(p534(2)*p612(3)-p612(2)*p534(3))+p
     & 534k0*(cfkz78(i5)%e(2)*p612(3)-p612(2)*cfkz78(i5)%e(3))-p
     & 612k0*(cfkz78(i5)%e(2)*p534(3)-p534(2)*cfkz78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz78(i5)%e(3)*p534k0+p534(3)*cfkz78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz78(i5)%e(3)*p612k0+p612(3)*cfkz78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz78(i5)%e(0)*p534(0)-cfkz78(i5)%e(1)*p534(1)-cfkz7
     & 8(i5)%e(2)*p534(2)-cfkz78(i5)%e(3)*p534(3)
      cvqd=cfkz78(i5)%e(0)*p612(0)-cfkz78(i5)%e(1)*p612(1)-cfkz7
     & 8(i5)%e(2)*p612(2)-cfkz78(i5)%e(3)*p612(3)
      cauxa=-cfkz78(i5)%ek0*quqd+p534k0*cvqd+p612k0*cvqu
      cauxb=-cfkz78(i5)%ek0*p612(2)+p612k0*cfkz78(i5)%e(2)
      cauxc=+cfkz78(i5)%ek0*p534(2)-p534k0*cfkz78(i5)%e(2)
      u534_78(i5)%a(1)=zcr(id5)*(cauxa+ceps_0)
      u534_78(i5)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u534_78(i5)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u534_78(i5)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      u534_78(i5)%c(1)=zcr(id5)*(cauxc+ceps_1)
      u534_78(i5)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u534_78(i5)%d(1)=zcl(id5)*cfkz78(i5)%ek0
      u534_78(i5)%d(2)=zcr(id5)*cfkz78(i5)%ek0
      end do
  
      if (ineutri(id5).ne.1.and.ineutri(id7).ne.1) then
      do i5=1,2
* T0 -- qu=p534,qd=p612,v=cfkf78(i5)%e,a=u534_78(i5)%a,b=u534_78(i5)%b,c
* =u534_78(i5)%c,d=u534_78(i5)%d,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cfkf78(i5)%ek0*(p534(2)*p612(3)-p612(2)*p534(3))+p
     & 534k0*(cfkf78(i5)%e(2)*p612(3)-p612(2)*cfkf78(i5)%e(3))-p
     & 612k0*(cfkf78(i5)%e(2)*p534(3)-p534(2)*cfkf78(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf78(i5)%e(3)*p534k0+p534(3)*cfkf78(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf78(i5)%e(3)*p612k0+p612(3)*cfkf78(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf78(i5)%e(0)*p534(0)-cfkf78(i5)%e(1)*p534(1)-cfkf7
     & 8(i5)%e(2)*p534(2)-cfkf78(i5)%e(3)*p534(3)
      cvqd=cfkf78(i5)%e(0)*p612(0)-cfkf78(i5)%e(1)*p612(1)-cfkf7
     & 8(i5)%e(2)*p612(2)-cfkf78(i5)%e(3)*p612(3)
      cauxa=-cfkf78(i5)%ek0*quqd+p534k0*cvqd+p612k0*cvqu
      cauxb=-cfkf78(i5)%ek0*p612(2)+p612k0*cfkf78(i5)%e(2)
      cauxc=+cfkf78(i5)%ek0*p534(2)-p534k0*cfkf78(i5)%e(2)
      u534_78(i5)%a(1)=u534_78(i5)%a(1)+fcr(id5)*(cauxa+ceps_0)
      u534_78(i5)%a(2)=u534_78(i5)%a(2)+fcl(id5)*(cauxa-ceps_0)
      u534_78(i5)%b(1)=u534_78(i5)%b(1)+fcl(id5)*(cauxb-ceps_2)
      u534_78(i5)%b(2)=u534_78(i5)%b(2)+fcr(id5)*(-cauxb-ceps_2)
      u534_78(i5)%c(1)=u534_78(i5)%c(1)+fcr(id5)*(cauxc+ceps_1)
      u534_78(i5)%c(2)=u534_78(i5)%c(2)+fcl(id5)*(-cauxc+ceps_1)
      u534_78(i5)%d(1)=u534_78(i5)%d(1)+fcl(id5)*cfkf78(i5)%ek0
      u534_78(i5)%d(2)=u534_78(i5)%d(2)+fcr(id5)*cfkf78(i5)%ek0
      end do
      endif
  
* quqd -- p=p578,q=p634
      quqd=p578(0)*p634(0)-p578(1)*p634(1)-p578(2)*p634(2)-p578(
     & 3)*p634(3)
      do i5=1,2
* T0 -- qu=p578,qd=p634,v=cfkz12(i5)%e,a=u578_12(i5)%a,b=u578_12(i5)%b,c
* =u578_12(i5)%c,d=u578_12(i5)%d,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cfkz12(i5)%ek0*(p578(2)*p634(3)-p634(2)*p578(3))+p
     & 578k0*(cfkz12(i5)%e(2)*p634(3)-p634(2)*cfkz12(i5)%e(3))-p
     & 634k0*(cfkz12(i5)%e(2)*p578(3)-p578(2)*cfkz12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz12(i5)%e(3)*p578k0+p578(3)*cfkz12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz12(i5)%e(3)*p634k0+p634(3)*cfkz12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz12(i5)%e(0)*p578(0)-cfkz12(i5)%e(1)*p578(1)-cfkz1
     & 2(i5)%e(2)*p578(2)-cfkz12(i5)%e(3)*p578(3)
      cvqd=cfkz12(i5)%e(0)*p634(0)-cfkz12(i5)%e(1)*p634(1)-cfkz1
     & 2(i5)%e(2)*p634(2)-cfkz12(i5)%e(3)*p634(3)
      cauxa=-cfkz12(i5)%ek0*quqd+p578k0*cvqd+p634k0*cvqu
      cauxb=-cfkz12(i5)%ek0*p634(2)+p634k0*cfkz12(i5)%e(2)
      cauxc=+cfkz12(i5)%ek0*p578(2)-p578k0*cfkz12(i5)%e(2)
      u578_12(i5)%a(1)=zcr(id5)*(cauxa+ceps_0)
      u578_12(i5)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u578_12(i5)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u578_12(i5)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      u578_12(i5)%c(1)=zcr(id5)*(cauxc+ceps_1)
      u578_12(i5)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u578_12(i5)%d(1)=zcl(id5)*cfkz12(i5)%ek0
      u578_12(i5)%d(2)=zcr(id5)*cfkz12(i5)%ek0
      end do
  
      if (ineutri(id5).ne.1.and.ineutri(id1).ne.1) then
      do i5=1,2
* T0 -- qu=p578,qd=p634,v=cfkf12(i5)%e,a=u578_12(i5)%a,b=u578_12(i5)%b,c
* =u578_12(i5)%c,d=u578_12(i5)%d,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cfkf12(i5)%ek0*(p578(2)*p634(3)-p634(2)*p578(3))+p
     & 578k0*(cfkf12(i5)%e(2)*p634(3)-p634(2)*cfkf12(i5)%e(3))-p
     & 634k0*(cfkf12(i5)%e(2)*p578(3)-p578(2)*cfkf12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf12(i5)%e(3)*p578k0+p578(3)*cfkf12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf12(i5)%e(3)*p634k0+p634(3)*cfkf12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf12(i5)%e(0)*p578(0)-cfkf12(i5)%e(1)*p578(1)-cfkf1
     & 2(i5)%e(2)*p578(2)-cfkf12(i5)%e(3)*p578(3)
      cvqd=cfkf12(i5)%e(0)*p634(0)-cfkf12(i5)%e(1)*p634(1)-cfkf1
     & 2(i5)%e(2)*p634(2)-cfkf12(i5)%e(3)*p634(3)
      cauxa=-cfkf12(i5)%ek0*quqd+p578k0*cvqd+p634k0*cvqu
      cauxb=-cfkf12(i5)%ek0*p634(2)+p634k0*cfkf12(i5)%e(2)
      cauxc=+cfkf12(i5)%ek0*p578(2)-p578k0*cfkf12(i5)%e(2)
      u578_12(i5)%a(1)=u578_12(i5)%a(1)+fcr(id5)*(cauxa+ceps_0)
      u578_12(i5)%a(2)=u578_12(i5)%a(2)+fcl(id5)*(cauxa-ceps_0)
      u578_12(i5)%b(1)=u578_12(i5)%b(1)+fcl(id5)*(cauxb-ceps_2)
      u578_12(i5)%b(2)=u578_12(i5)%b(2)+fcr(id5)*(-cauxb-ceps_2)
      u578_12(i5)%c(1)=u578_12(i5)%c(1)+fcr(id5)*(cauxc+ceps_1)
      u578_12(i5)%c(2)=u578_12(i5)%c(2)+fcl(id5)*(-cauxc+ceps_1)
      u578_12(i5)%d(1)=u578_12(i5)%d(1)+fcl(id5)*cfkf12(i5)%ek0
      u578_12(i5)%d(2)=u578_12(i5)%d(2)+fcr(id5)*cfkf12(i5)%ek0
      end do
      endif
  
* quqd -- p=p578,q=p612
      quqd=p578(0)*p612(0)-p578(1)*p612(1)-p578(2)*p612(2)-p578(
     & 3)*p612(3)
      do i5=1,2
* T0 -- qu=p578,qd=p612,v=cfkz34(i5)%e,a=u578_34(i5)%a,b=u578_34(i5)%b,c
* =u578_34(i5)%c,d=u578_34(i5)%d,cr=zcr(id5),cl=zcl(id5),nsum=0
      ceps_0=-cfkz34(i5)%ek0*(p578(2)*p612(3)-p612(2)*p578(3))+p
     & 578k0*(cfkz34(i5)%e(2)*p612(3)-p612(2)*cfkz34(i5)%e(3))-p
     & 612k0*(cfkz34(i5)%e(2)*p578(3)-p578(2)*cfkz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz34(i5)%e(3)*p578k0+p578(3)*cfkz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz34(i5)%e(3)*p612k0+p612(3)*cfkz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz34(i5)%e(0)*p578(0)-cfkz34(i5)%e(1)*p578(1)-cfkz3
     & 4(i5)%e(2)*p578(2)-cfkz34(i5)%e(3)*p578(3)
      cvqd=cfkz34(i5)%e(0)*p612(0)-cfkz34(i5)%e(1)*p612(1)-cfkz3
     & 4(i5)%e(2)*p612(2)-cfkz34(i5)%e(3)*p612(3)
      cauxa=-cfkz34(i5)%ek0*quqd+p578k0*cvqd+p612k0*cvqu
      cauxb=-cfkz34(i5)%ek0*p612(2)+p612k0*cfkz34(i5)%e(2)
      cauxc=+cfkz34(i5)%ek0*p578(2)-p578k0*cfkz34(i5)%e(2)
      u578_34(i5)%a(1)=zcr(id5)*(cauxa+ceps_0)
      u578_34(i5)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u578_34(i5)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u578_34(i5)%b(2)=zcr(id5)*(-cauxb-ceps_2)
      u578_34(i5)%c(1)=zcr(id5)*(cauxc+ceps_1)
      u578_34(i5)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u578_34(i5)%d(1)=zcl(id5)*cfkz34(i5)%ek0
      u578_34(i5)%d(2)=zcr(id5)*cfkz34(i5)%ek0
      end do
  
      if (ineutri(id5).ne.1.and.ineutri(id3).ne.1) then
      do i5=1,2
* T0 -- qu=p578,qd=p612,v=cfkf34(i5)%e,a=u578_34(i5)%a,b=u578_34(i5)%b,c
* =u578_34(i5)%c,d=u578_34(i5)%d,cr=fcr(id5),cl=fcl(id5),nsum=1
      ceps_0=-cfkf34(i5)%ek0*(p578(2)*p612(3)-p612(2)*p578(3))+p
     & 578k0*(cfkf34(i5)%e(2)*p612(3)-p612(2)*cfkf34(i5)%e(3))-p
     & 612k0*(cfkf34(i5)%e(2)*p578(3)-p578(2)*cfkf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf34(i5)%e(3)*p578k0+p578(3)*cfkf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf34(i5)%e(3)*p612k0+p612(3)*cfkf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf34(i5)%e(0)*p578(0)-cfkf34(i5)%e(1)*p578(1)-cfkf3
     & 4(i5)%e(2)*p578(2)-cfkf34(i5)%e(3)*p578(3)
      cvqd=cfkf34(i5)%e(0)*p612(0)-cfkf34(i5)%e(1)*p612(1)-cfkf3
     & 4(i5)%e(2)*p612(2)-cfkf34(i5)%e(3)*p612(3)
      cauxa=-cfkf34(i5)%ek0*quqd+p578k0*cvqd+p612k0*cvqu
      cauxb=-cfkf34(i5)%ek0*p612(2)+p612k0*cfkf34(i5)%e(2)
      cauxc=+cfkf34(i5)%ek0*p578(2)-p578k0*cfkf34(i5)%e(2)
      u578_34(i5)%a(1)=u578_34(i5)%a(1)+fcr(id5)*(cauxa+ceps_0)
      u578_34(i5)%a(2)=u578_34(i5)%a(2)+fcl(id5)*(cauxa-ceps_0)
      u578_34(i5)%b(1)=u578_34(i5)%b(1)+fcl(id5)*(cauxb-ceps_2)
      u578_34(i5)%b(2)=u578_34(i5)%b(2)+fcr(id5)*(-cauxb-ceps_2)
      u578_34(i5)%c(1)=u578_34(i5)%c(1)+fcr(id5)*(cauxc+ceps_1)
      u578_34(i5)%c(2)=u578_34(i5)%c(2)+fcl(id5)*(-cauxc+ceps_1)
      u578_34(i5)%d(1)=u578_34(i5)%d(1)+fcl(id5)*cfkf34(i5)%ek0
      u578_34(i5)%d(2)=u578_34(i5)%d(2)+fcr(id5)*cfkf34(i5)%ek0
      end do
      endif
  
* quqd -- p=p712,q=p856
      quqd=p712(0)*p856(0)-p712(1)*p856(1)-p712(2)*p856(2)-p712(
     & 3)*p856(3)
      do i5=1,2
* T0 -- qu=p712,qd=p856,v=cfkz34(i5)%e,a=u712_34(i5)%a,b=u712_34(i5)%b,c
* =u712_34(i5)%c,d=u712_34(i5)%d,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cfkz34(i5)%ek0*(p712(2)*p856(3)-p856(2)*p712(3))+p
     & 712k0*(cfkz34(i5)%e(2)*p856(3)-p856(2)*cfkz34(i5)%e(3))-p
     & 856k0*(cfkz34(i5)%e(2)*p712(3)-p712(2)*cfkz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz34(i5)%e(3)*p712k0+p712(3)*cfkz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz34(i5)%e(3)*p856k0+p856(3)*cfkz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz34(i5)%e(0)*p712(0)-cfkz34(i5)%e(1)*p712(1)-cfkz3
     & 4(i5)%e(2)*p712(2)-cfkz34(i5)%e(3)*p712(3)
      cvqd=cfkz34(i5)%e(0)*p856(0)-cfkz34(i5)%e(1)*p856(1)-cfkz3
     & 4(i5)%e(2)*p856(2)-cfkz34(i5)%e(3)*p856(3)
      cauxa=-cfkz34(i5)%ek0*quqd+p712k0*cvqd+p856k0*cvqu
      cauxb=-cfkz34(i5)%ek0*p856(2)+p856k0*cfkz34(i5)%e(2)
      cauxc=+cfkz34(i5)%ek0*p712(2)-p712k0*cfkz34(i5)%e(2)
      u712_34(i5)%a(1)=zcr(id7)*(cauxa+ceps_0)
      u712_34(i5)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u712_34(i5)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u712_34(i5)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      u712_34(i5)%c(1)=zcr(id7)*(cauxc+ceps_1)
      u712_34(i5)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u712_34(i5)%d(1)=zcl(id7)*cfkz34(i5)%ek0
      u712_34(i5)%d(2)=zcr(id7)*cfkz34(i5)%ek0
      end do
  
      if (ineutri(id7).ne.1.and.ineutri(id3).ne.1) then
      do i5=1,2
* T0 -- qu=p712,qd=p856,v=cfkf34(i5)%e,a=u712_34(i5)%a,b=u712_34(i5)%b,c
* =u712_34(i5)%c,d=u712_34(i5)%d,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cfkf34(i5)%ek0*(p712(2)*p856(3)-p856(2)*p712(3))+p
     & 712k0*(cfkf34(i5)%e(2)*p856(3)-p856(2)*cfkf34(i5)%e(3))-p
     & 856k0*(cfkf34(i5)%e(2)*p712(3)-p712(2)*cfkf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf34(i5)%e(3)*p712k0+p712(3)*cfkf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf34(i5)%e(3)*p856k0+p856(3)*cfkf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf34(i5)%e(0)*p712(0)-cfkf34(i5)%e(1)*p712(1)-cfkf3
     & 4(i5)%e(2)*p712(2)-cfkf34(i5)%e(3)*p712(3)
      cvqd=cfkf34(i5)%e(0)*p856(0)-cfkf34(i5)%e(1)*p856(1)-cfkf3
     & 4(i5)%e(2)*p856(2)-cfkf34(i5)%e(3)*p856(3)
      cauxa=-cfkf34(i5)%ek0*quqd+p712k0*cvqd+p856k0*cvqu
      cauxb=-cfkf34(i5)%ek0*p856(2)+p856k0*cfkf34(i5)%e(2)
      cauxc=+cfkf34(i5)%ek0*p712(2)-p712k0*cfkf34(i5)%e(2)
      u712_34(i5)%a(1)=u712_34(i5)%a(1)+fcr(id7)*(cauxa+ceps_0)
      u712_34(i5)%a(2)=u712_34(i5)%a(2)+fcl(id7)*(cauxa-ceps_0)
      u712_34(i5)%b(1)=u712_34(i5)%b(1)+fcl(id7)*(cauxb-ceps_2)
      u712_34(i5)%b(2)=u712_34(i5)%b(2)+fcr(id7)*(-cauxb-ceps_2)
      u712_34(i5)%c(1)=u712_34(i5)%c(1)+fcr(id7)*(cauxc+ceps_1)
      u712_34(i5)%c(2)=u712_34(i5)%c(2)+fcl(id7)*(-cauxc+ceps_1)
      u712_34(i5)%d(1)=u712_34(i5)%d(1)+fcl(id7)*cfkf34(i5)%ek0
      u712_34(i5)%d(2)=u712_34(i5)%d(2)+fcr(id7)*cfkf34(i5)%ek0
      end do
      endif
  
* quqd -- p=p712,q=p834
      quqd=p712(0)*p834(0)-p712(1)*p834(1)-p712(2)*p834(2)-p712(
     & 3)*p834(3)
      do i5=1,2
* T0 -- qu=p712,qd=p834,v=cfkz56(i5)%e,a=u712_56(i5)%a,b=u712_56(i5)%b,c
* =u712_56(i5)%c,d=u712_56(i5)%d,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cfkz56(i5)%ek0*(p712(2)*p834(3)-p834(2)*p712(3))+p
     & 712k0*(cfkz56(i5)%e(2)*p834(3)-p834(2)*cfkz56(i5)%e(3))-p
     & 834k0*(cfkz56(i5)%e(2)*p712(3)-p712(2)*cfkz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz56(i5)%e(3)*p712k0+p712(3)*cfkz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz56(i5)%e(3)*p834k0+p834(3)*cfkz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz56(i5)%e(0)*p712(0)-cfkz56(i5)%e(1)*p712(1)-cfkz5
     & 6(i5)%e(2)*p712(2)-cfkz56(i5)%e(3)*p712(3)
      cvqd=cfkz56(i5)%e(0)*p834(0)-cfkz56(i5)%e(1)*p834(1)-cfkz5
     & 6(i5)%e(2)*p834(2)-cfkz56(i5)%e(3)*p834(3)
      cauxa=-cfkz56(i5)%ek0*quqd+p712k0*cvqd+p834k0*cvqu
      cauxb=-cfkz56(i5)%ek0*p834(2)+p834k0*cfkz56(i5)%e(2)
      cauxc=+cfkz56(i5)%ek0*p712(2)-p712k0*cfkz56(i5)%e(2)
      u712_56(i5)%a(1)=zcr(id7)*(cauxa+ceps_0)
      u712_56(i5)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u712_56(i5)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u712_56(i5)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      u712_56(i5)%c(1)=zcr(id7)*(cauxc+ceps_1)
      u712_56(i5)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u712_56(i5)%d(1)=zcl(id7)*cfkz56(i5)%ek0
      u712_56(i5)%d(2)=zcr(id7)*cfkz56(i5)%ek0
      end do
  
      if (ineutri(id7).ne.1.and.ineutri(id5).ne.1) then
      do i5=1,2
* T0 -- qu=p712,qd=p834,v=cfkf56(i5)%e,a=u712_56(i5)%a,b=u712_56(i5)%b,c
* =u712_56(i5)%c,d=u712_56(i5)%d,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cfkf56(i5)%ek0*(p712(2)*p834(3)-p834(2)*p712(3))+p
     & 712k0*(cfkf56(i5)%e(2)*p834(3)-p834(2)*cfkf56(i5)%e(3))-p
     & 834k0*(cfkf56(i5)%e(2)*p712(3)-p712(2)*cfkf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf56(i5)%e(3)*p712k0+p712(3)*cfkf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf56(i5)%e(3)*p834k0+p834(3)*cfkf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf56(i5)%e(0)*p712(0)-cfkf56(i5)%e(1)*p712(1)-cfkf5
     & 6(i5)%e(2)*p712(2)-cfkf56(i5)%e(3)*p712(3)
      cvqd=cfkf56(i5)%e(0)*p834(0)-cfkf56(i5)%e(1)*p834(1)-cfkf5
     & 6(i5)%e(2)*p834(2)-cfkf56(i5)%e(3)*p834(3)
      cauxa=-cfkf56(i5)%ek0*quqd+p712k0*cvqd+p834k0*cvqu
      cauxb=-cfkf56(i5)%ek0*p834(2)+p834k0*cfkf56(i5)%e(2)
      cauxc=+cfkf56(i5)%ek0*p712(2)-p712k0*cfkf56(i5)%e(2)
      u712_56(i5)%a(1)=u712_56(i5)%a(1)+fcr(id7)*(cauxa+ceps_0)
      u712_56(i5)%a(2)=u712_56(i5)%a(2)+fcl(id7)*(cauxa-ceps_0)
      u712_56(i5)%b(1)=u712_56(i5)%b(1)+fcl(id7)*(cauxb-ceps_2)
      u712_56(i5)%b(2)=u712_56(i5)%b(2)+fcr(id7)*(-cauxb-ceps_2)
      u712_56(i5)%c(1)=u712_56(i5)%c(1)+fcr(id7)*(cauxc+ceps_1)
      u712_56(i5)%c(2)=u712_56(i5)%c(2)+fcl(id7)*(-cauxc+ceps_1)
      u712_56(i5)%d(1)=u712_56(i5)%d(1)+fcl(id7)*cfkf56(i5)%ek0
      u712_56(i5)%d(2)=u712_56(i5)%d(2)+fcr(id7)*cfkf56(i5)%ek0
      end do
      endif
  
* quqd -- p=p734,q=p856
      quqd=p734(0)*p856(0)-p734(1)*p856(1)-p734(2)*p856(2)-p734(
     & 3)*p856(3)
      do i5=1,2
* T0 -- qu=p734,qd=p856,v=cfkz12(i5)%e,a=u734_12(i5)%a,b=u734_12(i5)%b,c
* =u734_12(i5)%c,d=u734_12(i5)%d,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cfkz12(i5)%ek0*(p734(2)*p856(3)-p856(2)*p734(3))+p
     & 734k0*(cfkz12(i5)%e(2)*p856(3)-p856(2)*cfkz12(i5)%e(3))-p
     & 856k0*(cfkz12(i5)%e(2)*p734(3)-p734(2)*cfkz12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz12(i5)%e(3)*p734k0+p734(3)*cfkz12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz12(i5)%e(3)*p856k0+p856(3)*cfkz12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz12(i5)%e(0)*p734(0)-cfkz12(i5)%e(1)*p734(1)-cfkz1
     & 2(i5)%e(2)*p734(2)-cfkz12(i5)%e(3)*p734(3)
      cvqd=cfkz12(i5)%e(0)*p856(0)-cfkz12(i5)%e(1)*p856(1)-cfkz1
     & 2(i5)%e(2)*p856(2)-cfkz12(i5)%e(3)*p856(3)
      cauxa=-cfkz12(i5)%ek0*quqd+p734k0*cvqd+p856k0*cvqu
      cauxb=-cfkz12(i5)%ek0*p856(2)+p856k0*cfkz12(i5)%e(2)
      cauxc=+cfkz12(i5)%ek0*p734(2)-p734k0*cfkz12(i5)%e(2)
      u734_12(i5)%a(1)=zcr(id7)*(cauxa+ceps_0)
      u734_12(i5)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u734_12(i5)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u734_12(i5)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      u734_12(i5)%c(1)=zcr(id7)*(cauxc+ceps_1)
      u734_12(i5)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u734_12(i5)%d(1)=zcl(id7)*cfkz12(i5)%ek0
      u734_12(i5)%d(2)=zcr(id7)*cfkz12(i5)%ek0
      end do
  
      if (ineutri(id7).ne.1.and.ineutri(id1).ne.1) then
      do i5=1,2
* T0 -- qu=p734,qd=p856,v=cfkf12(i5)%e,a=u734_12(i5)%a,b=u734_12(i5)%b,c
* =u734_12(i5)%c,d=u734_12(i5)%d,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cfkf12(i5)%ek0*(p734(2)*p856(3)-p856(2)*p734(3))+p
     & 734k0*(cfkf12(i5)%e(2)*p856(3)-p856(2)*cfkf12(i5)%e(3))-p
     & 856k0*(cfkf12(i5)%e(2)*p734(3)-p734(2)*cfkf12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf12(i5)%e(3)*p734k0+p734(3)*cfkf12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf12(i5)%e(3)*p856k0+p856(3)*cfkf12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf12(i5)%e(0)*p734(0)-cfkf12(i5)%e(1)*p734(1)-cfkf1
     & 2(i5)%e(2)*p734(2)-cfkf12(i5)%e(3)*p734(3)
      cvqd=cfkf12(i5)%e(0)*p856(0)-cfkf12(i5)%e(1)*p856(1)-cfkf1
     & 2(i5)%e(2)*p856(2)-cfkf12(i5)%e(3)*p856(3)
      cauxa=-cfkf12(i5)%ek0*quqd+p734k0*cvqd+p856k0*cvqu
      cauxb=-cfkf12(i5)%ek0*p856(2)+p856k0*cfkf12(i5)%e(2)
      cauxc=+cfkf12(i5)%ek0*p734(2)-p734k0*cfkf12(i5)%e(2)
      u734_12(i5)%a(1)=u734_12(i5)%a(1)+fcr(id7)*(cauxa+ceps_0)
      u734_12(i5)%a(2)=u734_12(i5)%a(2)+fcl(id7)*(cauxa-ceps_0)
      u734_12(i5)%b(1)=u734_12(i5)%b(1)+fcl(id7)*(cauxb-ceps_2)
      u734_12(i5)%b(2)=u734_12(i5)%b(2)+fcr(id7)*(-cauxb-ceps_2)
      u734_12(i5)%c(1)=u734_12(i5)%c(1)+fcr(id7)*(cauxc+ceps_1)
      u734_12(i5)%c(2)=u734_12(i5)%c(2)+fcl(id7)*(-cauxc+ceps_1)
      u734_12(i5)%d(1)=u734_12(i5)%d(1)+fcl(id7)*cfkf12(i5)%ek0
      u734_12(i5)%d(2)=u734_12(i5)%d(2)+fcr(id7)*cfkf12(i5)%ek0
      end do
      endif
  
* quqd -- p=p734,q=p812
      quqd=p734(0)*p812(0)-p734(1)*p812(1)-p734(2)*p812(2)-p734(
     & 3)*p812(3)
      do i5=1,2
* T0 -- qu=p734,qd=p812,v=cfkz56(i5)%e,a=u734_56(i5)%a,b=u734_56(i5)%b,c
* =u734_56(i5)%c,d=u734_56(i5)%d,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cfkz56(i5)%ek0*(p734(2)*p812(3)-p812(2)*p734(3))+p
     & 734k0*(cfkz56(i5)%e(2)*p812(3)-p812(2)*cfkz56(i5)%e(3))-p
     & 812k0*(cfkz56(i5)%e(2)*p734(3)-p734(2)*cfkz56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz56(i5)%e(3)*p734k0+p734(3)*cfkz56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz56(i5)%e(3)*p812k0+p812(3)*cfkz56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz56(i5)%e(0)*p734(0)-cfkz56(i5)%e(1)*p734(1)-cfkz5
     & 6(i5)%e(2)*p734(2)-cfkz56(i5)%e(3)*p734(3)
      cvqd=cfkz56(i5)%e(0)*p812(0)-cfkz56(i5)%e(1)*p812(1)-cfkz5
     & 6(i5)%e(2)*p812(2)-cfkz56(i5)%e(3)*p812(3)
      cauxa=-cfkz56(i5)%ek0*quqd+p734k0*cvqd+p812k0*cvqu
      cauxb=-cfkz56(i5)%ek0*p812(2)+p812k0*cfkz56(i5)%e(2)
      cauxc=+cfkz56(i5)%ek0*p734(2)-p734k0*cfkz56(i5)%e(2)
      u734_56(i5)%a(1)=zcr(id7)*(cauxa+ceps_0)
      u734_56(i5)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u734_56(i5)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u734_56(i5)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      u734_56(i5)%c(1)=zcr(id7)*(cauxc+ceps_1)
      u734_56(i5)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u734_56(i5)%d(1)=zcl(id7)*cfkz56(i5)%ek0
      u734_56(i5)%d(2)=zcr(id7)*cfkz56(i5)%ek0
      end do
  
      if (ineutri(id7).ne.1.and.ineutri(id5).ne.1) then
      do i5=1,2
* T0 -- qu=p734,qd=p812,v=cfkf56(i5)%e,a=u734_56(i5)%a,b=u734_56(i5)%b,c
* =u734_56(i5)%c,d=u734_56(i5)%d,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cfkf56(i5)%ek0*(p734(2)*p812(3)-p812(2)*p734(3))+p
     & 734k0*(cfkf56(i5)%e(2)*p812(3)-p812(2)*cfkf56(i5)%e(3))-p
     & 812k0*(cfkf56(i5)%e(2)*p734(3)-p734(2)*cfkf56(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf56(i5)%e(3)*p734k0+p734(3)*cfkf56(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf56(i5)%e(3)*p812k0+p812(3)*cfkf56(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf56(i5)%e(0)*p734(0)-cfkf56(i5)%e(1)*p734(1)-cfkf5
     & 6(i5)%e(2)*p734(2)-cfkf56(i5)%e(3)*p734(3)
      cvqd=cfkf56(i5)%e(0)*p812(0)-cfkf56(i5)%e(1)*p812(1)-cfkf5
     & 6(i5)%e(2)*p812(2)-cfkf56(i5)%e(3)*p812(3)
      cauxa=-cfkf56(i5)%ek0*quqd+p734k0*cvqd+p812k0*cvqu
      cauxb=-cfkf56(i5)%ek0*p812(2)+p812k0*cfkf56(i5)%e(2)
      cauxc=+cfkf56(i5)%ek0*p734(2)-p734k0*cfkf56(i5)%e(2)
      u734_56(i5)%a(1)=u734_56(i5)%a(1)+fcr(id7)*(cauxa+ceps_0)
      u734_56(i5)%a(2)=u734_56(i5)%a(2)+fcl(id7)*(cauxa-ceps_0)
      u734_56(i5)%b(1)=u734_56(i5)%b(1)+fcl(id7)*(cauxb-ceps_2)
      u734_56(i5)%b(2)=u734_56(i5)%b(2)+fcr(id7)*(-cauxb-ceps_2)
      u734_56(i5)%c(1)=u734_56(i5)%c(1)+fcr(id7)*(cauxc+ceps_1)
      u734_56(i5)%c(2)=u734_56(i5)%c(2)+fcl(id7)*(-cauxc+ceps_1)
      u734_56(i5)%d(1)=u734_56(i5)%d(1)+fcl(id7)*cfkf56(i5)%ek0
      u734_56(i5)%d(2)=u734_56(i5)%d(2)+fcr(id7)*cfkf56(i5)%ek0
      end do
      endif
  
* quqd -- p=p756,q=p834
      quqd=p756(0)*p834(0)-p756(1)*p834(1)-p756(2)*p834(2)-p756(
     & 3)*p834(3)
      do i5=1,2
* T0 -- qu=p756,qd=p834,v=cfkz12(i5)%e,a=u756_12(i5)%a,b=u756_12(i5)%b,c
* =u756_12(i5)%c,d=u756_12(i5)%d,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cfkz12(i5)%ek0*(p756(2)*p834(3)-p834(2)*p756(3))+p
     & 756k0*(cfkz12(i5)%e(2)*p834(3)-p834(2)*cfkz12(i5)%e(3))-p
     & 834k0*(cfkz12(i5)%e(2)*p756(3)-p756(2)*cfkz12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz12(i5)%e(3)*p756k0+p756(3)*cfkz12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz12(i5)%e(3)*p834k0+p834(3)*cfkz12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz12(i5)%e(0)*p756(0)-cfkz12(i5)%e(1)*p756(1)-cfkz1
     & 2(i5)%e(2)*p756(2)-cfkz12(i5)%e(3)*p756(3)
      cvqd=cfkz12(i5)%e(0)*p834(0)-cfkz12(i5)%e(1)*p834(1)-cfkz1
     & 2(i5)%e(2)*p834(2)-cfkz12(i5)%e(3)*p834(3)
      cauxa=-cfkz12(i5)%ek0*quqd+p756k0*cvqd+p834k0*cvqu
      cauxb=-cfkz12(i5)%ek0*p834(2)+p834k0*cfkz12(i5)%e(2)
      cauxc=+cfkz12(i5)%ek0*p756(2)-p756k0*cfkz12(i5)%e(2)
      u756_12(i5)%a(1)=zcr(id7)*(cauxa+ceps_0)
      u756_12(i5)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u756_12(i5)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u756_12(i5)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      u756_12(i5)%c(1)=zcr(id7)*(cauxc+ceps_1)
      u756_12(i5)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u756_12(i5)%d(1)=zcl(id7)*cfkz12(i5)%ek0
      u756_12(i5)%d(2)=zcr(id7)*cfkz12(i5)%ek0
      end do
  
      if (ineutri(id7).ne.1.and.ineutri(id1).ne.1) then
      do i5=1,2
* T0 -- qu=p756,qd=p834,v=cfkf12(i5)%e,a=u756_12(i5)%a,b=u756_12(i5)%b,c
* =u756_12(i5)%c,d=u756_12(i5)%d,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cfkf12(i5)%ek0*(p756(2)*p834(3)-p834(2)*p756(3))+p
     & 756k0*(cfkf12(i5)%e(2)*p834(3)-p834(2)*cfkf12(i5)%e(3))-p
     & 834k0*(cfkf12(i5)%e(2)*p756(3)-p756(2)*cfkf12(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf12(i5)%e(3)*p756k0+p756(3)*cfkf12(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf12(i5)%e(3)*p834k0+p834(3)*cfkf12(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf12(i5)%e(0)*p756(0)-cfkf12(i5)%e(1)*p756(1)-cfkf1
     & 2(i5)%e(2)*p756(2)-cfkf12(i5)%e(3)*p756(3)
      cvqd=cfkf12(i5)%e(0)*p834(0)-cfkf12(i5)%e(1)*p834(1)-cfkf1
     & 2(i5)%e(2)*p834(2)-cfkf12(i5)%e(3)*p834(3)
      cauxa=-cfkf12(i5)%ek0*quqd+p756k0*cvqd+p834k0*cvqu
      cauxb=-cfkf12(i5)%ek0*p834(2)+p834k0*cfkf12(i5)%e(2)
      cauxc=+cfkf12(i5)%ek0*p756(2)-p756k0*cfkf12(i5)%e(2)
      u756_12(i5)%a(1)=u756_12(i5)%a(1)+fcr(id7)*(cauxa+ceps_0)
      u756_12(i5)%a(2)=u756_12(i5)%a(2)+fcl(id7)*(cauxa-ceps_0)
      u756_12(i5)%b(1)=u756_12(i5)%b(1)+fcl(id7)*(cauxb-ceps_2)
      u756_12(i5)%b(2)=u756_12(i5)%b(2)+fcr(id7)*(-cauxb-ceps_2)
      u756_12(i5)%c(1)=u756_12(i5)%c(1)+fcr(id7)*(cauxc+ceps_1)
      u756_12(i5)%c(2)=u756_12(i5)%c(2)+fcl(id7)*(-cauxc+ceps_1)
      u756_12(i5)%d(1)=u756_12(i5)%d(1)+fcl(id7)*cfkf12(i5)%ek0
      u756_12(i5)%d(2)=u756_12(i5)%d(2)+fcr(id7)*cfkf12(i5)%ek0
      end do
      endif
  
* quqd -- p=p756,q=p812
      quqd=p756(0)*p812(0)-p756(1)*p812(1)-p756(2)*p812(2)-p756(
     & 3)*p812(3)
      do i5=1,2
* T0 -- qu=p756,qd=p812,v=cfkz34(i5)%e,a=u756_34(i5)%a,b=u756_34(i5)%b,c
* =u756_34(i5)%c,d=u756_34(i5)%d,cr=zcr(id7),cl=zcl(id7),nsum=0
      ceps_0=-cfkz34(i5)%ek0*(p756(2)*p812(3)-p812(2)*p756(3))+p
     & 756k0*(cfkz34(i5)%e(2)*p812(3)-p812(2)*cfkz34(i5)%e(3))-p
     & 812k0*(cfkz34(i5)%e(2)*p756(3)-p756(2)*cfkz34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkz34(i5)%e(3)*p756k0+p756(3)*cfkz34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkz34(i5)%e(3)*p812k0+p812(3)*cfkz34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkz34(i5)%e(0)*p756(0)-cfkz34(i5)%e(1)*p756(1)-cfkz3
     & 4(i5)%e(2)*p756(2)-cfkz34(i5)%e(3)*p756(3)
      cvqd=cfkz34(i5)%e(0)*p812(0)-cfkz34(i5)%e(1)*p812(1)-cfkz3
     & 4(i5)%e(2)*p812(2)-cfkz34(i5)%e(3)*p812(3)
      cauxa=-cfkz34(i5)%ek0*quqd+p756k0*cvqd+p812k0*cvqu
      cauxb=-cfkz34(i5)%ek0*p812(2)+p812k0*cfkz34(i5)%e(2)
      cauxc=+cfkz34(i5)%ek0*p756(2)-p756k0*cfkz34(i5)%e(2)
      u756_34(i5)%a(1)=zcr(id7)*(cauxa+ceps_0)
      u756_34(i5)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u756_34(i5)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u756_34(i5)%b(2)=zcr(id7)*(-cauxb-ceps_2)
      u756_34(i5)%c(1)=zcr(id7)*(cauxc+ceps_1)
      u756_34(i5)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u756_34(i5)%d(1)=zcl(id7)*cfkz34(i5)%ek0
      u756_34(i5)%d(2)=zcr(id7)*cfkz34(i5)%ek0
      end do
  
      if (ineutri(id7).ne.1.and.ineutri(id3).ne.1) then
      do i5=1,2
* T0 -- qu=p756,qd=p812,v=cfkf34(i5)%e,a=u756_34(i5)%a,b=u756_34(i5)%b,c
* =u756_34(i5)%c,d=u756_34(i5)%d,cr=fcr(id7),cl=fcl(id7),nsum=1
      ceps_0=-cfkf34(i5)%ek0*(p756(2)*p812(3)-p812(2)*p756(3))+p
     & 756k0*(cfkf34(i5)%e(2)*p812(3)-p812(2)*cfkf34(i5)%e(3))-p
     & 812k0*(cfkf34(i5)%e(2)*p756(3)-p756(2)*cfkf34(i5)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cfkf34(i5)%e(3)*p756k0+p756(3)*cfkf34(i5)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cfkf34(i5)%e(3)*p812k0+p812(3)*cfkf34(i5)%ek0
      ceps_2=ceps_2*cim
      cvqu=cfkf34(i5)%e(0)*p756(0)-cfkf34(i5)%e(1)*p756(1)-cfkf3
     & 4(i5)%e(2)*p756(2)-cfkf34(i5)%e(3)*p756(3)
      cvqd=cfkf34(i5)%e(0)*p812(0)-cfkf34(i5)%e(1)*p812(1)-cfkf3
     & 4(i5)%e(2)*p812(2)-cfkf34(i5)%e(3)*p812(3)
      cauxa=-cfkf34(i5)%ek0*quqd+p756k0*cvqd+p812k0*cvqu
      cauxb=-cfkf34(i5)%ek0*p812(2)+p812k0*cfkf34(i5)%e(2)
      cauxc=+cfkf34(i5)%ek0*p756(2)-p756k0*cfkf34(i5)%e(2)
      u756_34(i5)%a(1)=u756_34(i5)%a(1)+fcr(id7)*(cauxa+ceps_0)
      u756_34(i5)%a(2)=u756_34(i5)%a(2)+fcl(id7)*(cauxa-ceps_0)
      u756_34(i5)%b(1)=u756_34(i5)%b(1)+fcl(id7)*(cauxb-ceps_2)
      u756_34(i5)%b(2)=u756_34(i5)%b(2)+fcr(id7)*(-cauxb-ceps_2)
      u756_34(i5)%c(1)=u756_34(i5)%c(1)+fcr(id7)*(cauxc+ceps_1)
      u756_34(i5)%c(2)=u756_34(i5)%c(2)+fcl(id7)*(-cauxc+ceps_1)
      u756_34(i5)%d(1)=u756_34(i5)%d(1)+fcl(id7)*cfkf34(i5)%ek0
      u756_34(i5)%d(2)=u756_34(i5)%d(2)+fcr(id7)*cfkf34(i5)%ek0
      end do
      endif
  
  
*                                                                       
* LEFT*MIDDLE                                                           
*                                                                       
  
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l1_3456(i3,i5)%a,cc=l1_3456(i3,i5)%c,a1=l1_34(i3)%a,c1=l1_3
* 4(i3)%c,a2=u134_56(i5)%a,b2=u134_56(i5)%b,c2=u134_56(i5)%c,d2=u134_56(
* i5)%d,prq=p134q,nsum=0
      l1_3456(i3,i5)%a(1)=l1_34(i3)%a(1)*u134_56(i5)%a(1)+l1_34(
     & i3)%c(1)*p134q*u134_56(i5)%b(2)
      l1_3456(i3,i5)%c(1)=l1_34(i3)%a(1)*u134_56(i5)%c(1)+l1_34(
     & i3)%c(1)*p134q*u134_56(i5)%d(2)
      l1_3456(i3,i5)%c(2)=l1_34(i3)%c(2)*p134q*u134_56(i5)%d(1)+
     & l1_34(i3)%a(2)*u134_56(i5)%c(2)
      l1_3456(i3,i5)%a(2)=l1_34(i3)%c(2)*p134q*u134_56(i5)%b(1)+
     & l1_34(i3)%a(2)*u134_56(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l1_3478(i3,i5)%a,cc=l1_3478(i3,i5)%c,a1=l1_34(i3)%a,c1=l1_3
* 4(i3)%c,a2=u134_78(i5)%a,b2=u134_78(i5)%b,c2=u134_78(i5)%c,d2=u134_78(
* i5)%d,prq=p134q,nsum=0
      l1_3478(i3,i5)%a(1)=l1_34(i3)%a(1)*u134_78(i5)%a(1)+l1_34(
     & i3)%c(1)*p134q*u134_78(i5)%b(2)
      l1_3478(i3,i5)%c(1)=l1_34(i3)%a(1)*u134_78(i5)%c(1)+l1_34(
     & i3)%c(1)*p134q*u134_78(i5)%d(2)
      l1_3478(i3,i5)%c(2)=l1_34(i3)%c(2)*p134q*u134_78(i5)%d(1)+
     & l1_34(i3)%a(2)*u134_78(i5)%c(2)
      l1_3478(i3,i5)%a(2)=l1_34(i3)%c(2)*p134q*u134_78(i5)%b(1)+
     & l1_34(i3)%a(2)*u134_78(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l1_5634(i3,i5)%a,cc=l1_5634(i3,i5)%c,a1=l1_56(i3)%a,c1=l1_5
* 6(i3)%c,a2=u156_34(i5)%a,b2=u156_34(i5)%b,c2=u156_34(i5)%c,d2=u156_34(
* i5)%d,prq=p156q,nsum=0
      l1_5634(i3,i5)%a(1)=l1_56(i3)%a(1)*u156_34(i5)%a(1)+l1_56(
     & i3)%c(1)*p156q*u156_34(i5)%b(2)
      l1_5634(i3,i5)%c(1)=l1_56(i3)%a(1)*u156_34(i5)%c(1)+l1_56(
     & i3)%c(1)*p156q*u156_34(i5)%d(2)
      l1_5634(i3,i5)%c(2)=l1_56(i3)%c(2)*p156q*u156_34(i5)%d(1)+
     & l1_56(i3)%a(2)*u156_34(i5)%c(2)
      l1_5634(i3,i5)%a(2)=l1_56(i3)%c(2)*p156q*u156_34(i5)%b(1)+
     & l1_56(i3)%a(2)*u156_34(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l1_5678(i3,i5)%a,cc=l1_5678(i3,i5)%c,a1=l1_56(i3)%a,c1=l1_5
* 6(i3)%c,a2=u156_78(i5)%a,b2=u156_78(i5)%b,c2=u156_78(i5)%c,d2=u156_78(
* i5)%d,prq=p156q,nsum=0
      l1_5678(i3,i5)%a(1)=l1_56(i3)%a(1)*u156_78(i5)%a(1)+l1_56(
     & i3)%c(1)*p156q*u156_78(i5)%b(2)
      l1_5678(i3,i5)%c(1)=l1_56(i3)%a(1)*u156_78(i5)%c(1)+l1_56(
     & i3)%c(1)*p156q*u156_78(i5)%d(2)
      l1_5678(i3,i5)%c(2)=l1_56(i3)%c(2)*p156q*u156_78(i5)%d(1)+
     & l1_56(i3)%a(2)*u156_78(i5)%c(2)
      l1_5678(i3,i5)%a(2)=l1_56(i3)%c(2)*p156q*u156_78(i5)%b(1)+
     & l1_56(i3)%a(2)*u156_78(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l1_7834(i3,i5)%a,cc=l1_7834(i3,i5)%c,a1=l1_78(i3)%a,c1=l1_7
* 8(i3)%c,a2=u178_34(i5)%a,b2=u178_34(i5)%b,c2=u178_34(i5)%c,d2=u178_34(
* i5)%d,prq=p178q,nsum=0
      l1_7834(i3,i5)%a(1)=l1_78(i3)%a(1)*u178_34(i5)%a(1)+l1_78(
     & i3)%c(1)*p178q*u178_34(i5)%b(2)
      l1_7834(i3,i5)%c(1)=l1_78(i3)%a(1)*u178_34(i5)%c(1)+l1_78(
     & i3)%c(1)*p178q*u178_34(i5)%d(2)
      l1_7834(i3,i5)%c(2)=l1_78(i3)%c(2)*p178q*u178_34(i5)%d(1)+
     & l1_78(i3)%a(2)*u178_34(i5)%c(2)
      l1_7834(i3,i5)%a(2)=l1_78(i3)%c(2)*p178q*u178_34(i5)%b(1)+
     & l1_78(i3)%a(2)*u178_34(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l1_7856(i3,i5)%a,cc=l1_7856(i3,i5)%c,a1=l1_78(i3)%a,c1=l1_7
* 8(i3)%c,a2=u178_56(i5)%a,b2=u178_56(i5)%b,c2=u178_56(i5)%c,d2=u178_56(
* i5)%d,prq=p178q,nsum=0
      l1_7856(i3,i5)%a(1)=l1_78(i3)%a(1)*u178_56(i5)%a(1)+l1_78(
     & i3)%c(1)*p178q*u178_56(i5)%b(2)
      l1_7856(i3,i5)%c(1)=l1_78(i3)%a(1)*u178_56(i5)%c(1)+l1_78(
     & i3)%c(1)*p178q*u178_56(i5)%d(2)
      l1_7856(i3,i5)%c(2)=l1_78(i3)%c(2)*p178q*u178_56(i5)%d(1)+
     & l1_78(i3)%a(2)*u178_56(i5)%c(2)
      l1_7856(i3,i5)%a(2)=l1_78(i3)%c(2)*p178q*u178_56(i5)%b(1)+
     & l1_78(i3)%a(2)*u178_56(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l3_1256(i3,i5)%a,cc=l3_1256(i3,i5)%c,a1=l3_12(i3)%a,c1=l3_1
* 2(i3)%c,a2=u312_56(i5)%a,b2=u312_56(i5)%b,c2=u312_56(i5)%c,d2=u312_56(
* i5)%d,prq=p312q,nsum=0
      l3_1256(i3,i5)%a(1)=l3_12(i3)%a(1)*u312_56(i5)%a(1)+l3_12(
     & i3)%c(1)*p312q*u312_56(i5)%b(2)
      l3_1256(i3,i5)%c(1)=l3_12(i3)%a(1)*u312_56(i5)%c(1)+l3_12(
     & i3)%c(1)*p312q*u312_56(i5)%d(2)
      l3_1256(i3,i5)%c(2)=l3_12(i3)%c(2)*p312q*u312_56(i5)%d(1)+
     & l3_12(i3)%a(2)*u312_56(i5)%c(2)
      l3_1256(i3,i5)%a(2)=l3_12(i3)%c(2)*p312q*u312_56(i5)%b(1)+
     & l3_12(i3)%a(2)*u312_56(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l3_1278(i3,i5)%a,cc=l3_1278(i3,i5)%c,a1=l3_12(i3)%a,c1=l3_1
* 2(i3)%c,a2=u312_78(i5)%a,b2=u312_78(i5)%b,c2=u312_78(i5)%c,d2=u312_78(
* i5)%d,prq=p312q,nsum=0
      l3_1278(i3,i5)%a(1)=l3_12(i3)%a(1)*u312_78(i5)%a(1)+l3_12(
     & i3)%c(1)*p312q*u312_78(i5)%b(2)
      l3_1278(i3,i5)%c(1)=l3_12(i3)%a(1)*u312_78(i5)%c(1)+l3_12(
     & i3)%c(1)*p312q*u312_78(i5)%d(2)
      l3_1278(i3,i5)%c(2)=l3_12(i3)%c(2)*p312q*u312_78(i5)%d(1)+
     & l3_12(i3)%a(2)*u312_78(i5)%c(2)
      l3_1278(i3,i5)%a(2)=l3_12(i3)%c(2)*p312q*u312_78(i5)%b(1)+
     & l3_12(i3)%a(2)*u312_78(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l3_5612(i3,i5)%a,cc=l3_5612(i3,i5)%c,a1=l3_56(i3)%a,c1=l3_5
* 6(i3)%c,a2=u356_12(i5)%a,b2=u356_12(i5)%b,c2=u356_12(i5)%c,d2=u356_12(
* i5)%d,prq=p356q,nsum=0
      l3_5612(i3,i5)%a(1)=l3_56(i3)%a(1)*u356_12(i5)%a(1)+l3_56(
     & i3)%c(1)*p356q*u356_12(i5)%b(2)
      l3_5612(i3,i5)%c(1)=l3_56(i3)%a(1)*u356_12(i5)%c(1)+l3_56(
     & i3)%c(1)*p356q*u356_12(i5)%d(2)
      l3_5612(i3,i5)%c(2)=l3_56(i3)%c(2)*p356q*u356_12(i5)%d(1)+
     & l3_56(i3)%a(2)*u356_12(i5)%c(2)
      l3_5612(i3,i5)%a(2)=l3_56(i3)%c(2)*p356q*u356_12(i5)%b(1)+
     & l3_56(i3)%a(2)*u356_12(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l3_5678(i3,i5)%a,cc=l3_5678(i3,i5)%c,a1=l3_56(i3)%a,c1=l3_5
* 6(i3)%c,a2=u356_78(i5)%a,b2=u356_78(i5)%b,c2=u356_78(i5)%c,d2=u356_78(
* i5)%d,prq=p356q,nsum=0
      l3_5678(i3,i5)%a(1)=l3_56(i3)%a(1)*u356_78(i5)%a(1)+l3_56(
     & i3)%c(1)*p356q*u356_78(i5)%b(2)
      l3_5678(i3,i5)%c(1)=l3_56(i3)%a(1)*u356_78(i5)%c(1)+l3_56(
     & i3)%c(1)*p356q*u356_78(i5)%d(2)
      l3_5678(i3,i5)%c(2)=l3_56(i3)%c(2)*p356q*u356_78(i5)%d(1)+
     & l3_56(i3)%a(2)*u356_78(i5)%c(2)
      l3_5678(i3,i5)%a(2)=l3_56(i3)%c(2)*p356q*u356_78(i5)%b(1)+
     & l3_56(i3)%a(2)*u356_78(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l3_7812(i3,i5)%a,cc=l3_7812(i3,i5)%c,a1=l3_78(i3)%a,c1=l3_7
* 8(i3)%c,a2=u378_12(i5)%a,b2=u378_12(i5)%b,c2=u378_12(i5)%c,d2=u378_12(
* i5)%d,prq=p378q,nsum=0
      l3_7812(i3,i5)%a(1)=l3_78(i3)%a(1)*u378_12(i5)%a(1)+l3_78(
     & i3)%c(1)*p378q*u378_12(i5)%b(2)
      l3_7812(i3,i5)%c(1)=l3_78(i3)%a(1)*u378_12(i5)%c(1)+l3_78(
     & i3)%c(1)*p378q*u378_12(i5)%d(2)
      l3_7812(i3,i5)%c(2)=l3_78(i3)%c(2)*p378q*u378_12(i5)%d(1)+
     & l3_78(i3)%a(2)*u378_12(i5)%c(2)
      l3_7812(i3,i5)%a(2)=l3_78(i3)%c(2)*p378q*u378_12(i5)%b(1)+
     & l3_78(i3)%a(2)*u378_12(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l3_7856(i3,i5)%a,cc=l3_7856(i3,i5)%c,a1=l3_78(i3)%a,c1=l3_7
* 8(i3)%c,a2=u378_56(i5)%a,b2=u378_56(i5)%b,c2=u378_56(i5)%c,d2=u378_56(
* i5)%d,prq=p378q,nsum=0
      l3_7856(i3,i5)%a(1)=l3_78(i3)%a(1)*u378_56(i5)%a(1)+l3_78(
     & i3)%c(1)*p378q*u378_56(i5)%b(2)
      l3_7856(i3,i5)%c(1)=l3_78(i3)%a(1)*u378_56(i5)%c(1)+l3_78(
     & i3)%c(1)*p378q*u378_56(i5)%d(2)
      l3_7856(i3,i5)%c(2)=l3_78(i3)%c(2)*p378q*u378_56(i5)%d(1)+
     & l3_78(i3)%a(2)*u378_56(i5)%c(2)
      l3_7856(i3,i5)%a(2)=l3_78(i3)%c(2)*p378q*u378_56(i5)%b(1)+
     & l3_78(i3)%a(2)*u378_56(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l5_1234(i3,i5)%a,cc=l5_1234(i3,i5)%c,a1=l5_12(i3)%a,c1=l5_1
* 2(i3)%c,a2=u512_34(i5)%a,b2=u512_34(i5)%b,c2=u512_34(i5)%c,d2=u512_34(
* i5)%d,prq=p512q,nsum=0
      l5_1234(i3,i5)%a(1)=l5_12(i3)%a(1)*u512_34(i5)%a(1)+l5_12(
     & i3)%c(1)*p512q*u512_34(i5)%b(2)
      l5_1234(i3,i5)%c(1)=l5_12(i3)%a(1)*u512_34(i5)%c(1)+l5_12(
     & i3)%c(1)*p512q*u512_34(i5)%d(2)
      l5_1234(i3,i5)%c(2)=l5_12(i3)%c(2)*p512q*u512_34(i5)%d(1)+
     & l5_12(i3)%a(2)*u512_34(i5)%c(2)
      l5_1234(i3,i5)%a(2)=l5_12(i3)%c(2)*p512q*u512_34(i5)%b(1)+
     & l5_12(i3)%a(2)*u512_34(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l5_1278(i3,i5)%a,cc=l5_1278(i3,i5)%c,a1=l5_12(i3)%a,c1=l5_1
* 2(i3)%c,a2=u512_78(i5)%a,b2=u512_78(i5)%b,c2=u512_78(i5)%c,d2=u512_78(
* i5)%d,prq=p512q,nsum=0
      l5_1278(i3,i5)%a(1)=l5_12(i3)%a(1)*u512_78(i5)%a(1)+l5_12(
     & i3)%c(1)*p512q*u512_78(i5)%b(2)
      l5_1278(i3,i5)%c(1)=l5_12(i3)%a(1)*u512_78(i5)%c(1)+l5_12(
     & i3)%c(1)*p512q*u512_78(i5)%d(2)
      l5_1278(i3,i5)%c(2)=l5_12(i3)%c(2)*p512q*u512_78(i5)%d(1)+
     & l5_12(i3)%a(2)*u512_78(i5)%c(2)
      l5_1278(i3,i5)%a(2)=l5_12(i3)%c(2)*p512q*u512_78(i5)%b(1)+
     & l5_12(i3)%a(2)*u512_78(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l5_3412(i3,i5)%a,cc=l5_3412(i3,i5)%c,a1=l5_34(i3)%a,c1=l5_3
* 4(i3)%c,a2=u534_12(i5)%a,b2=u534_12(i5)%b,c2=u534_12(i5)%c,d2=u534_12(
* i5)%d,prq=p534q,nsum=0
      l5_3412(i3,i5)%a(1)=l5_34(i3)%a(1)*u534_12(i5)%a(1)+l5_34(
     & i3)%c(1)*p534q*u534_12(i5)%b(2)
      l5_3412(i3,i5)%c(1)=l5_34(i3)%a(1)*u534_12(i5)%c(1)+l5_34(
     & i3)%c(1)*p534q*u534_12(i5)%d(2)
      l5_3412(i3,i5)%c(2)=l5_34(i3)%c(2)*p534q*u534_12(i5)%d(1)+
     & l5_34(i3)%a(2)*u534_12(i5)%c(2)
      l5_3412(i3,i5)%a(2)=l5_34(i3)%c(2)*p534q*u534_12(i5)%b(1)+
     & l5_34(i3)%a(2)*u534_12(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l5_3478(i3,i5)%a,cc=l5_3478(i3,i5)%c,a1=l5_34(i3)%a,c1=l5_3
* 4(i3)%c,a2=u534_78(i5)%a,b2=u534_78(i5)%b,c2=u534_78(i5)%c,d2=u534_78(
* i5)%d,prq=p534q,nsum=0
      l5_3478(i3,i5)%a(1)=l5_34(i3)%a(1)*u534_78(i5)%a(1)+l5_34(
     & i3)%c(1)*p534q*u534_78(i5)%b(2)
      l5_3478(i3,i5)%c(1)=l5_34(i3)%a(1)*u534_78(i5)%c(1)+l5_34(
     & i3)%c(1)*p534q*u534_78(i5)%d(2)
      l5_3478(i3,i5)%c(2)=l5_34(i3)%c(2)*p534q*u534_78(i5)%d(1)+
     & l5_34(i3)%a(2)*u534_78(i5)%c(2)
      l5_3478(i3,i5)%a(2)=l5_34(i3)%c(2)*p534q*u534_78(i5)%b(1)+
     & l5_34(i3)%a(2)*u534_78(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l5_7812(i3,i5)%a,cc=l5_7812(i3,i5)%c,a1=l5_78(i3)%a,c1=l5_7
* 8(i3)%c,a2=u578_12(i5)%a,b2=u578_12(i5)%b,c2=u578_12(i5)%c,d2=u578_12(
* i5)%d,prq=p578q,nsum=0
      l5_7812(i3,i5)%a(1)=l5_78(i3)%a(1)*u578_12(i5)%a(1)+l5_78(
     & i3)%c(1)*p578q*u578_12(i5)%b(2)
      l5_7812(i3,i5)%c(1)=l5_78(i3)%a(1)*u578_12(i5)%c(1)+l5_78(
     & i3)%c(1)*p578q*u578_12(i5)%d(2)
      l5_7812(i3,i5)%c(2)=l5_78(i3)%c(2)*p578q*u578_12(i5)%d(1)+
     & l5_78(i3)%a(2)*u578_12(i5)%c(2)
      l5_7812(i3,i5)%a(2)=l5_78(i3)%c(2)*p578q*u578_12(i5)%b(1)+
     & l5_78(i3)%a(2)*u578_12(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l5_7834(i3,i5)%a,cc=l5_7834(i3,i5)%c,a1=l5_78(i3)%a,c1=l5_7
* 8(i3)%c,a2=u578_34(i5)%a,b2=u578_34(i5)%b,c2=u578_34(i5)%c,d2=u578_34(
* i5)%d,prq=p578q,nsum=0
      l5_7834(i3,i5)%a(1)=l5_78(i3)%a(1)*u578_34(i5)%a(1)+l5_78(
     & i3)%c(1)*p578q*u578_34(i5)%b(2)
      l5_7834(i3,i5)%c(1)=l5_78(i3)%a(1)*u578_34(i5)%c(1)+l5_78(
     & i3)%c(1)*p578q*u578_34(i5)%d(2)
      l5_7834(i3,i5)%c(2)=l5_78(i3)%c(2)*p578q*u578_34(i5)%d(1)+
     & l5_78(i3)%a(2)*u578_34(i5)%c(2)
      l5_7834(i3,i5)%a(2)=l5_78(i3)%c(2)*p578q*u578_34(i5)%b(1)+
     & l5_78(i3)%a(2)*u578_34(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l7_1234(i3,i5)%a,cc=l7_1234(i3,i5)%c,a1=l7_12(i3)%a,c1=l7_1
* 2(i3)%c,a2=u712_34(i5)%a,b2=u712_34(i5)%b,c2=u712_34(i5)%c,d2=u712_34(
* i5)%d,prq=p712q,nsum=0
      l7_1234(i3,i5)%a(1)=l7_12(i3)%a(1)*u712_34(i5)%a(1)+l7_12(
     & i3)%c(1)*p712q*u712_34(i5)%b(2)
      l7_1234(i3,i5)%c(1)=l7_12(i3)%a(1)*u712_34(i5)%c(1)+l7_12(
     & i3)%c(1)*p712q*u712_34(i5)%d(2)
      l7_1234(i3,i5)%c(2)=l7_12(i3)%c(2)*p712q*u712_34(i5)%d(1)+
     & l7_12(i3)%a(2)*u712_34(i5)%c(2)
      l7_1234(i3,i5)%a(2)=l7_12(i3)%c(2)*p712q*u712_34(i5)%b(1)+
     & l7_12(i3)%a(2)*u712_34(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l7_1256(i3,i5)%a,cc=l7_1256(i3,i5)%c,a1=l7_12(i3)%a,c1=l7_1
* 2(i3)%c,a2=u712_56(i5)%a,b2=u712_56(i5)%b,c2=u712_56(i5)%c,d2=u712_56(
* i5)%d,prq=p712q,nsum=0
      l7_1256(i3,i5)%a(1)=l7_12(i3)%a(1)*u712_56(i5)%a(1)+l7_12(
     & i3)%c(1)*p712q*u712_56(i5)%b(2)
      l7_1256(i3,i5)%c(1)=l7_12(i3)%a(1)*u712_56(i5)%c(1)+l7_12(
     & i3)%c(1)*p712q*u712_56(i5)%d(2)
      l7_1256(i3,i5)%c(2)=l7_12(i3)%c(2)*p712q*u712_56(i5)%d(1)+
     & l7_12(i3)%a(2)*u712_56(i5)%c(2)
      l7_1256(i3,i5)%a(2)=l7_12(i3)%c(2)*p712q*u712_56(i5)%b(1)+
     & l7_12(i3)%a(2)*u712_56(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l7_3412(i3,i5)%a,cc=l7_3412(i3,i5)%c,a1=l7_34(i3)%a,c1=l7_3
* 4(i3)%c,a2=u734_12(i5)%a,b2=u734_12(i5)%b,c2=u734_12(i5)%c,d2=u734_12(
* i5)%d,prq=p734q,nsum=0
      l7_3412(i3,i5)%a(1)=l7_34(i3)%a(1)*u734_12(i5)%a(1)+l7_34(
     & i3)%c(1)*p734q*u734_12(i5)%b(2)
      l7_3412(i3,i5)%c(1)=l7_34(i3)%a(1)*u734_12(i5)%c(1)+l7_34(
     & i3)%c(1)*p734q*u734_12(i5)%d(2)
      l7_3412(i3,i5)%c(2)=l7_34(i3)%c(2)*p734q*u734_12(i5)%d(1)+
     & l7_34(i3)%a(2)*u734_12(i5)%c(2)
      l7_3412(i3,i5)%a(2)=l7_34(i3)%c(2)*p734q*u734_12(i5)%b(1)+
     & l7_34(i3)%a(2)*u734_12(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l7_3456(i3,i5)%a,cc=l7_3456(i3,i5)%c,a1=l7_34(i3)%a,c1=l7_3
* 4(i3)%c,a2=u734_56(i5)%a,b2=u734_56(i5)%b,c2=u734_56(i5)%c,d2=u734_56(
* i5)%d,prq=p734q,nsum=0
      l7_3456(i3,i5)%a(1)=l7_34(i3)%a(1)*u734_56(i5)%a(1)+l7_34(
     & i3)%c(1)*p734q*u734_56(i5)%b(2)
      l7_3456(i3,i5)%c(1)=l7_34(i3)%a(1)*u734_56(i5)%c(1)+l7_34(
     & i3)%c(1)*p734q*u734_56(i5)%d(2)
      l7_3456(i3,i5)%c(2)=l7_34(i3)%c(2)*p734q*u734_56(i5)%d(1)+
     & l7_34(i3)%a(2)*u734_56(i5)%c(2)
      l7_3456(i3,i5)%a(2)=l7_34(i3)%c(2)*p734q*u734_56(i5)%b(1)+
     & l7_34(i3)%a(2)*u734_56(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l7_5612(i3,i5)%a,cc=l7_5612(i3,i5)%c,a1=l7_56(i3)%a,c1=l7_5
* 6(i3)%c,a2=u756_12(i5)%a,b2=u756_12(i5)%b,c2=u756_12(i5)%c,d2=u756_12(
* i5)%d,prq=p756q,nsum=0
      l7_5612(i3,i5)%a(1)=l7_56(i3)%a(1)*u756_12(i5)%a(1)+l7_56(
     & i3)%c(1)*p756q*u756_12(i5)%b(2)
      l7_5612(i3,i5)%c(1)=l7_56(i3)%a(1)*u756_12(i5)%c(1)+l7_56(
     & i3)%c(1)*p756q*u756_12(i5)%d(2)
      l7_5612(i3,i5)%c(2)=l7_56(i3)%c(2)*p756q*u756_12(i5)%d(1)+
     & l7_56(i3)%a(2)*u756_12(i5)%c(2)
      l7_5612(i3,i5)%a(2)=l7_56(i3)%c(2)*p756q*u756_12(i5)%b(1)+
     & l7_56(i3)%a(2)*u756_12(i5)%a(2)
      end do
      end do
      do i3=1,2
      do i5=1,2
* TLT0 -- aa=l7_5634(i3,i5)%a,cc=l7_5634(i3,i5)%c,a1=l7_56(i3)%a,c1=l7_5
* 6(i3)%c,a2=u756_34(i5)%a,b2=u756_34(i5)%b,c2=u756_34(i5)%c,d2=u756_34(
* i5)%d,prq=p756q,nsum=0
      l7_5634(i3,i5)%a(1)=l7_56(i3)%a(1)*u756_34(i5)%a(1)+l7_56(
     & i3)%c(1)*p756q*u756_34(i5)%b(2)
      l7_5634(i3,i5)%c(1)=l7_56(i3)%a(1)*u756_34(i5)%c(1)+l7_56(
     & i3)%c(1)*p756q*u756_34(i5)%d(2)
      l7_5634(i3,i5)%c(2)=l7_56(i3)%c(2)*p756q*u756_34(i5)%d(1)+
     & l7_56(i3)%a(2)*u756_34(i5)%c(2)
      l7_5634(i3,i5)%a(2)=l7_56(i3)%c(2)*p756q*u756_34(i5)%b(1)+
     & l7_56(i3)%a(2)*u756_34(i5)%a(2)
      end do
      end do
  
*                                                                       
*LEFT*MIDDLE*RIGHT (3 FORKS)                                            
*                                                                       
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk12_345678(&,i3,i5,i7),a1=l1_3456(i3,i5)%a,c1=l1_3456(i
* 3,i5)%c,a2=r2_78(i7)%a,b2=r2_78(i7)%b,prq=p278q,bef=,aft=
      c3fk12_345678(1,i3,i5,i7)=(l1_3456(i3,i5)%a(1)*r2_78(i7)%a
     & (1)+l1_3456(i3,i5)%c(1)*p278q*r2_78(i7)%b(2))
      c3fk12_345678(2,i3,i5,i7)=(l1_3456(i3,i5)%c(2)*p278q*r2_78
     & (i7)%b(1)+l1_3456(i3,i5)%a(2)*r2_78(i7)%a(2))
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk12_345678(&,i3,i5,i7),a1=l1_5634(i5,i3)%a,c1=l1_5634(i
* 5,i3)%c,a2=r2_78(i7)%a,b2=r2_78(i7)%b,prq=p278q,bef=,aft=+c3fk12_34567
* 8(&,i3,i5,i7)
      c3fk12_345678(1,i3,i5,i7)=(l1_5634(i5,i3)%a(1)*r2_78(i7)%a
     & (1)+l1_5634(i5,i3)%c(1)*p278q*r2_78(i7)%b(2))+c3fk12_3456
     & 78(1,i3,i5,i7)
      c3fk12_345678(2,i3,i5,i7)=(l1_5634(i5,i3)%c(2)*p278q*r2_78
     & (i7)%b(1)+l1_5634(i5,i3)%a(2)*r2_78(i7)%a(2))+c3fk12_3456
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk12_345678(&,i3,i5,i7),a1=l1_3478(i3,i7)%a,c1=l1_3478(i
* 3,i7)%c,a2=r2_56(i5)%a,b2=r2_56(i5)%b,prq=p256q,bef=,aft=+c3fk12_34567
* 8(&,i3,i5,i7)
      c3fk12_345678(1,i3,i5,i7)=(l1_3478(i3,i7)%a(1)*r2_56(i5)%a
     & (1)+l1_3478(i3,i7)%c(1)*p256q*r2_56(i5)%b(2))+c3fk12_3456
     & 78(1,i3,i5,i7)
      c3fk12_345678(2,i3,i5,i7)=(l1_3478(i3,i7)%c(2)*p256q*r2_56
     & (i5)%b(1)+l1_3478(i3,i7)%a(2)*r2_56(i5)%a(2))+c3fk12_3456
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk12_345678(&,i3,i5,i7),a1=l1_7834(i7,i3)%a,c1=l1_7834(i
* 7,i3)%c,a2=r2_56(i5)%a,b2=r2_56(i5)%b,prq=p256q,bef=,aft=+c3fk12_34567
* 8(&,i3,i5,i7)
      c3fk12_345678(1,i3,i5,i7)=(l1_7834(i7,i3)%a(1)*r2_56(i5)%a
     & (1)+l1_7834(i7,i3)%c(1)*p256q*r2_56(i5)%b(2))+c3fk12_3456
     & 78(1,i3,i5,i7)
      c3fk12_345678(2,i3,i5,i7)=(l1_7834(i7,i3)%c(2)*p256q*r2_56
     & (i5)%b(1)+l1_7834(i7,i3)%a(2)*r2_56(i5)%a(2))+c3fk12_3456
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk12_345678(&,i3,i5,i7),a1=l1_5678(i5,i7)%a,c1=l1_5678(i
* 5,i7)%c,a2=r2_34(i3)%a,b2=r2_34(i3)%b,prq=p234q,bef=,aft=+c3fk12_34567
* 8(&,i3,i5,i7)
      c3fk12_345678(1,i3,i5,i7)=(l1_5678(i5,i7)%a(1)*r2_34(i3)%a
     & (1)+l1_5678(i5,i7)%c(1)*p234q*r2_34(i3)%b(2))+c3fk12_3456
     & 78(1,i3,i5,i7)
      c3fk12_345678(2,i3,i5,i7)=(l1_5678(i5,i7)%c(2)*p234q*r2_34
     & (i3)%b(1)+l1_5678(i5,i7)%a(2)*r2_34(i3)%a(2))+c3fk12_3456
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk12_345678(&,i3,i5,i7),a1=l1_7856(i7,i5)%a,c1=l1_7856(i
* 7,i5)%c,a2=r2_34(i3)%a,b2=r2_34(i3)%b,prq=p234q,bef=,aft=+c3fk12_34567
* 8(&,i3,i5,i7)
      c3fk12_345678(1,i3,i5,i7)=(l1_7856(i7,i5)%a(1)*r2_34(i3)%a
     & (1)+l1_7856(i7,i5)%c(1)*p234q*r2_34(i3)%b(2))+c3fk12_3456
     & 78(1,i3,i5,i7)
      c3fk12_345678(2,i3,i5,i7)=(l1_7856(i7,i5)%c(2)*p234q*r2_34
     & (i3)%b(1)+l1_7856(i7,i5)%a(2)*r2_34(i3)%a(2))+c3fk12_3456
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk34_125678(&,i3,i5,i7),a1=l3_1256(i3,i5)%a,c1=l3_1256(i
* 3,i5)%c,a2=r4_78(i7)%a,b2=r4_78(i7)%b,prq=p478q,bef=,aft=
      c3fk34_125678(1,i3,i5,i7)=(l3_1256(i3,i5)%a(1)*r4_78(i7)%a
     & (1)+l3_1256(i3,i5)%c(1)*p478q*r4_78(i7)%b(2))
      c3fk34_125678(2,i3,i5,i7)=(l3_1256(i3,i5)%c(2)*p478q*r4_78
     & (i7)%b(1)+l3_1256(i3,i5)%a(2)*r4_78(i7)%a(2))
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk34_125678(&,i3,i5,i7),a1=l3_5612(i5,i3)%a,c1=l3_5612(i
* 5,i3)%c,a2=r4_78(i7)%a,b2=r4_78(i7)%b,prq=p478q,bef=,aft=+c3fk34_12567
* 8(&,i3,i5,i7)
      c3fk34_125678(1,i3,i5,i7)=(l3_5612(i5,i3)%a(1)*r4_78(i7)%a
     & (1)+l3_5612(i5,i3)%c(1)*p478q*r4_78(i7)%b(2))+c3fk34_1256
     & 78(1,i3,i5,i7)
      c3fk34_125678(2,i3,i5,i7)=(l3_5612(i5,i3)%c(2)*p478q*r4_78
     & (i7)%b(1)+l3_5612(i5,i3)%a(2)*r4_78(i7)%a(2))+c3fk34_1256
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk34_125678(&,i3,i5,i7),a1=l3_1278(i3,i7)%a,c1=l3_1278(i
* 3,i7)%c,a2=r4_56(i5)%a,b2=r4_56(i5)%b,prq=p456q,bef=,aft=+c3fk34_12567
* 8(&,i3,i5,i7)
      c3fk34_125678(1,i3,i5,i7)=(l3_1278(i3,i7)%a(1)*r4_56(i5)%a
     & (1)+l3_1278(i3,i7)%c(1)*p456q*r4_56(i5)%b(2))+c3fk34_1256
     & 78(1,i3,i5,i7)
      c3fk34_125678(2,i3,i5,i7)=(l3_1278(i3,i7)%c(2)*p456q*r4_56
     & (i5)%b(1)+l3_1278(i3,i7)%a(2)*r4_56(i5)%a(2))+c3fk34_1256
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk34_125678(&,i3,i5,i7),a1=l3_7812(i7,i3)%a,c1=l3_7812(i
* 7,i3)%c,a2=r4_56(i5)%a,b2=r4_56(i5)%b,prq=p456q,bef=,aft=+c3fk34_12567
* 8(&,i3,i5,i7)
      c3fk34_125678(1,i3,i5,i7)=(l3_7812(i7,i3)%a(1)*r4_56(i5)%a
     & (1)+l3_7812(i7,i3)%c(1)*p456q*r4_56(i5)%b(2))+c3fk34_1256
     & 78(1,i3,i5,i7)
      c3fk34_125678(2,i3,i5,i7)=(l3_7812(i7,i3)%c(2)*p456q*r4_56
     & (i5)%b(1)+l3_7812(i7,i3)%a(2)*r4_56(i5)%a(2))+c3fk34_1256
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk34_125678(&,i3,i5,i7),a1=l3_5678(i5,i7)%a,c1=l3_5678(i
* 5,i7)%c,a2=r4_12(i3)%a,b2=r4_12(i3)%b,prq=p412q,bef=,aft=+c3fk34_12567
* 8(&,i3,i5,i7)
      c3fk34_125678(1,i3,i5,i7)=(l3_5678(i5,i7)%a(1)*r4_12(i3)%a
     & (1)+l3_5678(i5,i7)%c(1)*p412q*r4_12(i3)%b(2))+c3fk34_1256
     & 78(1,i3,i5,i7)
      c3fk34_125678(2,i3,i5,i7)=(l3_5678(i5,i7)%c(2)*p412q*r4_12
     & (i3)%b(1)+l3_5678(i5,i7)%a(2)*r4_12(i3)%a(2))+c3fk34_1256
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk34_125678(&,i3,i5,i7),a1=l3_7856(i7,i5)%a,c1=l3_7856(i
* 7,i5)%c,a2=r4_12(i3)%a,b2=r4_12(i3)%b,prq=p412q,bef=,aft=+c3fk34_12567
* 8(&,i3,i5,i7)
      c3fk34_125678(1,i3,i5,i7)=(l3_7856(i7,i5)%a(1)*r4_12(i3)%a
     & (1)+l3_7856(i7,i5)%c(1)*p412q*r4_12(i3)%b(2))+c3fk34_1256
     & 78(1,i3,i5,i7)
      c3fk34_125678(2,i3,i5,i7)=(l3_7856(i7,i5)%c(2)*p412q*r4_12
     & (i3)%b(1)+l3_7856(i7,i5)%a(2)*r4_12(i3)%a(2))+c3fk34_1256
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk56_123478(&,i3,i5,i7),a1=l5_1234(i3,i5)%a,c1=l5_1234(i
* 3,i5)%c,a2=r6_78(i7)%a,b2=r6_78(i7)%b,prq=p678q,bef=,aft=
      c3fk56_123478(1,i3,i5,i7)=(l5_1234(i3,i5)%a(1)*r6_78(i7)%a
     & (1)+l5_1234(i3,i5)%c(1)*p678q*r6_78(i7)%b(2))
      c3fk56_123478(2,i3,i5,i7)=(l5_1234(i3,i5)%c(2)*p678q*r6_78
     & (i7)%b(1)+l5_1234(i3,i5)%a(2)*r6_78(i7)%a(2))
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk56_123478(&,i3,i5,i7),a1=l5_3412(i5,i3)%a,c1=l5_3412(i
* 5,i3)%c,a2=r6_78(i7)%a,b2=r6_78(i7)%b,prq=p678q,bef=,aft=+c3fk56_12347
* 8(&,i3,i5,i7)
      c3fk56_123478(1,i3,i5,i7)=(l5_3412(i5,i3)%a(1)*r6_78(i7)%a
     & (1)+l5_3412(i5,i3)%c(1)*p678q*r6_78(i7)%b(2))+c3fk56_1234
     & 78(1,i3,i5,i7)
      c3fk56_123478(2,i3,i5,i7)=(l5_3412(i5,i3)%c(2)*p678q*r6_78
     & (i7)%b(1)+l5_3412(i5,i3)%a(2)*r6_78(i7)%a(2))+c3fk56_1234
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk56_123478(&,i3,i5,i7),a1=l5_1278(i3,i7)%a,c1=l5_1278(i
* 3,i7)%c,a2=r6_34(i5)%a,b2=r6_34(i5)%b,prq=p634q,bef=,aft=+c3fk56_12347
* 8(&,i3,i5,i7)
      c3fk56_123478(1,i3,i5,i7)=(l5_1278(i3,i7)%a(1)*r6_34(i5)%a
     & (1)+l5_1278(i3,i7)%c(1)*p634q*r6_34(i5)%b(2))+c3fk56_1234
     & 78(1,i3,i5,i7)
      c3fk56_123478(2,i3,i5,i7)=(l5_1278(i3,i7)%c(2)*p634q*r6_34
     & (i5)%b(1)+l5_1278(i3,i7)%a(2)*r6_34(i5)%a(2))+c3fk56_1234
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk56_123478(&,i3,i5,i7),a1=l5_7812(i7,i3)%a,c1=l5_7812(i
* 7,i3)%c,a2=r6_34(i5)%a,b2=r6_34(i5)%b,prq=p634q,bef=,aft=+c3fk56_12347
* 8(&,i3,i5,i7)
      c3fk56_123478(1,i3,i5,i7)=(l5_7812(i7,i3)%a(1)*r6_34(i5)%a
     & (1)+l5_7812(i7,i3)%c(1)*p634q*r6_34(i5)%b(2))+c3fk56_1234
     & 78(1,i3,i5,i7)
      c3fk56_123478(2,i3,i5,i7)=(l5_7812(i7,i3)%c(2)*p634q*r6_34
     & (i5)%b(1)+l5_7812(i7,i3)%a(2)*r6_34(i5)%a(2))+c3fk56_1234
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk56_123478(&,i3,i5,i7),a1=l5_3478(i5,i7)%a,c1=l5_3478(i
* 5,i7)%c,a2=r6_12(i3)%a,b2=r6_12(i3)%b,prq=p612q,bef=,aft=+c3fk56_12347
* 8(&,i3,i5,i7)
      c3fk56_123478(1,i3,i5,i7)=(l5_3478(i5,i7)%a(1)*r6_12(i3)%a
     & (1)+l5_3478(i5,i7)%c(1)*p612q*r6_12(i3)%b(2))+c3fk56_1234
     & 78(1,i3,i5,i7)
      c3fk56_123478(2,i3,i5,i7)=(l5_3478(i5,i7)%c(2)*p612q*r6_12
     & (i3)%b(1)+l5_3478(i5,i7)%a(2)*r6_12(i3)%a(2))+c3fk56_1234
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk56_123478(&,i3,i5,i7),a1=l5_7834(i7,i5)%a,c1=l5_7834(i
* 7,i5)%c,a2=r6_12(i3)%a,b2=r6_12(i3)%b,prq=p612q,bef=,aft=+c3fk56_12347
* 8(&,i3,i5,i7)
      c3fk56_123478(1,i3,i5,i7)=(l5_7834(i7,i5)%a(1)*r6_12(i3)%a
     & (1)+l5_7834(i7,i5)%c(1)*p612q*r6_12(i3)%b(2))+c3fk56_1234
     & 78(1,i3,i5,i7)
      c3fk56_123478(2,i3,i5,i7)=(l5_7834(i7,i5)%c(2)*p612q*r6_12
     & (i3)%b(1)+l5_7834(i7,i5)%a(2)*r6_12(i3)%a(2))+c3fk56_1234
     & 78(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk78_123456(&,i3,i5,i7),a1=l7_1234(i3,i5)%a,c1=l7_1234(i
* 3,i5)%c,a2=r8_56(i7)%a,b2=r8_56(i7)%b,prq=p856q,bef=,aft=
      c3fk78_123456(1,i3,i5,i7)=(l7_1234(i3,i5)%a(1)*r8_56(i7)%a
     & (1)+l7_1234(i3,i5)%c(1)*p856q*r8_56(i7)%b(2))
      c3fk78_123456(2,i3,i5,i7)=(l7_1234(i3,i5)%c(2)*p856q*r8_56
     & (i7)%b(1)+l7_1234(i3,i5)%a(2)*r8_56(i7)%a(2))
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk78_123456(&,i3,i5,i7),a1=l7_3412(i5,i3)%a,c1=l7_3412(i
* 5,i3)%c,a2=r8_56(i7)%a,b2=r8_56(i7)%b,prq=p856q,bef=,aft=+c3fk78_12345
* 6(&,i3,i5,i7)
      c3fk78_123456(1,i3,i5,i7)=(l7_3412(i5,i3)%a(1)*r8_56(i7)%a
     & (1)+l7_3412(i5,i3)%c(1)*p856q*r8_56(i7)%b(2))+c3fk78_1234
     & 56(1,i3,i5,i7)
      c3fk78_123456(2,i3,i5,i7)=(l7_3412(i5,i3)%c(2)*p856q*r8_56
     & (i7)%b(1)+l7_3412(i5,i3)%a(2)*r8_56(i7)%a(2))+c3fk78_1234
     & 56(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk78_123456(&,i3,i5,i7),a1=l7_1256(i3,i7)%a,c1=l7_1256(i
* 3,i7)%c,a2=r8_34(i5)%a,b2=r8_34(i5)%b,prq=p834q,bef=,aft=+c3fk78_12345
* 6(&,i3,i5,i7)
      c3fk78_123456(1,i3,i5,i7)=(l7_1256(i3,i7)%a(1)*r8_34(i5)%a
     & (1)+l7_1256(i3,i7)%c(1)*p834q*r8_34(i5)%b(2))+c3fk78_1234
     & 56(1,i3,i5,i7)
      c3fk78_123456(2,i3,i5,i7)=(l7_1256(i3,i7)%c(2)*p834q*r8_34
     & (i5)%b(1)+l7_1256(i3,i7)%a(2)*r8_34(i5)%a(2))+c3fk78_1234
     & 56(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk78_123456(&,i3,i5,i7),a1=l7_5612(i7,i3)%a,c1=l7_5612(i
* 7,i3)%c,a2=r8_34(i5)%a,b2=r8_34(i5)%b,prq=p834q,bef=,aft=+c3fk78_12345
* 6(&,i3,i5,i7)
      c3fk78_123456(1,i3,i5,i7)=(l7_5612(i7,i3)%a(1)*r8_34(i5)%a
     & (1)+l7_5612(i7,i3)%c(1)*p834q*r8_34(i5)%b(2))+c3fk78_1234
     & 56(1,i3,i5,i7)
      c3fk78_123456(2,i3,i5,i7)=(l7_5612(i7,i3)%c(2)*p834q*r8_34
     & (i5)%b(1)+l7_5612(i7,i3)%a(2)*r8_34(i5)%a(2))+c3fk78_1234
     & 56(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk78_123456(&,i3,i5,i7),a1=l7_3456(i5,i7)%a,c1=l7_3456(i
* 5,i7)%c,a2=r8_12(i3)%a,b2=r8_12(i3)%b,prq=p812q,bef=,aft=+c3fk78_12345
* 6(&,i3,i5,i7)
      c3fk78_123456(1,i3,i5,i7)=(l7_3456(i5,i7)%a(1)*r8_12(i3)%a
     & (1)+l7_3456(i5,i7)%c(1)*p812q*r8_12(i3)%b(2))+c3fk78_1234
     & 56(1,i3,i5,i7)
      c3fk78_123456(2,i3,i5,i7)=(l7_3456(i5,i7)%c(2)*p812q*r8_12
     & (i3)%b(1)+l7_3456(i5,i7)%a(2)*r8_12(i3)%a(2))+c3fk78_1234
     & 56(2,i3,i5,i7)
      end do
      end do
      end do
  
      do i3=1,2
      do i5=1,2
      do i7=1,2
* TLTR0 -- aa=c3fk78_123456(&,i3,i5,i7),a1=l7_5634(i7,i5)%a,c1=l7_5634(i
* 7,i5)%c,a2=r8_12(i3)%a,b2=r8_12(i3)%b,prq=p812q,bef=,aft=+c3fk78_12345
* 6(&,i3,i5,i7)
      c3fk78_123456(1,i3,i5,i7)=(l7_5634(i7,i5)%a(1)*r8_12(i3)%a
     & (1)+l7_5634(i7,i5)%c(1)*p812q*r8_12(i3)%b(2))+c3fk78_1234
     & 56(1,i3,i5,i7)
      c3fk78_123456(2,i3,i5,i7)=(l7_5634(i7,i5)%c(2)*p812q*r8_12
     & (i3)%b(1)+l7_5634(i7,i5)%a(2)*r8_12(i3)%a(2))+c3fk78_1234
     & 56(2,i3,i5,i7)
      end do
      end do
      end do
  
  
  
* 3forks total amplitude                                                
  
*     do i1=1,2                                                         
*     do i3=1,2                                                         
*     do i5=1,2                                                         
*     do i7=1,2                                                         
*       c3fk_tot(i1,i3,i5,i7)= c3fk12_345678(i1,i3,i5,i7)+              
*    &   c3fk34_125678(i3,i1,i5,i7)+c3fk56_123478(i5,i1,i3,i7)+         
*    &   c3fk78_123456(i7,i1,i3,i5)                                     
*     enddo                                                             
*     enddo                                                             
*     enddo                                                             
*     enddo                                                             
  
*                                                                       
* LEFT VECTOR                                                           
*                                                                       
  
  
*zw                                                                     
  
  
      if(iresz(1).eq.1) then
          do m=0,3
          do i=1,2
            lz1_234(m)%a(i)=czero
            lz1_234(m)%c(i)=czero
            lf1_234(m)%a(i)=czero
            lf1_234(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p1,q=p234
      quqd=p1(0)*p234(0)-p1(1)*p234(1)-p1(2)*p234(2)-p1(3)*p234(
     & 3)
* TL0 -- qu=p1,qd=p234,v=0,a=lz1_234(0)%a,c=lz1_234(0)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      eps_0=-p1(2)*p234(3)+p234(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p234(0)+p234k0*p1(0)
      lz1_234(0)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_234(0)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(0)%c(1)=zcr(id1)*(p1(2)+ceps_1)
      lz1_234(0)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p234,v=1,a=lz1_234(1)%a,c=lz1_234(1)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      auxa=-quqd+p1k0*p234(1)+p234k0*p1(1)
      lz1_234(1)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_234(1)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(1)%c(1)=zcr(id1)*(p1(2)+ceps_1)
      lz1_234(1)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p234,v=2,a=lz1_234(2)%a,c=lz1_234(2)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      eps_0=-p1k0*p234(3)+p234k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p234(2)+p234k0*p1(2)
      lz1_234(2)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_234(2)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(2)%c(1)=zcr(id1)*p1k0
      lz1_234(2)%c(2)=-zcl(id1)*p1k0
* TL0 -- qu=p1,qd=p234,v=3,a=lz1_234(3)%a,c=lz1_234(3)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      eps_0=p1k0*p234(2)-p234k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p234(3)+p234k0*p1(3)
      lz1_234(3)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_234(3)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(3)%c(1)=zcr(id1)*ceps_1
      lz1_234(3)%c(2)=zcl(id1)*ceps_1
  
        if (ineutri(id1).ne.1) then
* TL0 -- qu=p1,qd=p234,v=0,a=lf1_234(0)%a,c=lf1_234(0)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      eps_0=-p1(2)*p234(3)+p234(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p234(0)+p234k0*p1(0)
      lf1_234(0)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_234(0)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(0)%c(1)=fcr(id1)*(p1(2)+ceps_1)
      lf1_234(0)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p234,v=1,a=lf1_234(1)%a,c=lf1_234(1)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      auxa=-quqd+p1k0*p234(1)+p234k0*p1(1)
      lf1_234(1)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_234(1)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(1)%c(1)=fcr(id1)*(p1(2)+ceps_1)
      lf1_234(1)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p234,v=2,a=lf1_234(2)%a,c=lf1_234(2)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      eps_0=-p1k0*p234(3)+p234k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p234(2)+p234k0*p1(2)
      lf1_234(2)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_234(2)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(2)%c(1)=fcr(id1)*p1k0
      lf1_234(2)%c(2)=-fcl(id1)*p1k0
* TL0 -- qu=p1,qd=p234,v=3,a=lf1_234(3)%a,c=lf1_234(3)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      eps_0=p1k0*p234(2)-p234k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p234(3)+p234k0*p1(3)
      lf1_234(3)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_234(3)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(3)%c(1)=fcr(id1)*ceps_1
      lf1_234(3)%c(2)=fcl(id1)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(1).eq.1) then
          do m=0,3
          do i=1,2
            lz1_256(m)%a(i)=czero
            lz1_256(m)%c(i)=czero
            lf1_256(m)%a(i)=czero
            lf1_256(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p1,q=p256
      quqd=p1(0)*p256(0)-p1(1)*p256(1)-p1(2)*p256(2)-p1(3)*p256(
     & 3)
* TL0 -- qu=p1,qd=p256,v=0,a=lz1_256(0)%a,c=lz1_256(0)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      eps_0=-p1(2)*p256(3)+p256(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p256(0)+p256k0*p1(0)
      lz1_256(0)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_256(0)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_256(0)%c(1)=zcr(id1)*(p1(2)+ceps_1)
      lz1_256(0)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p256,v=1,a=lz1_256(1)%a,c=lz1_256(1)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      auxa=-quqd+p1k0*p256(1)+p256k0*p1(1)
      lz1_256(1)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_256(1)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_256(1)%c(1)=zcr(id1)*(p1(2)+ceps_1)
      lz1_256(1)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p256,v=2,a=lz1_256(2)%a,c=lz1_256(2)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      eps_0=-p1k0*p256(3)+p256k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p256(2)+p256k0*p1(2)
      lz1_256(2)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_256(2)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_256(2)%c(1)=zcr(id1)*p1k0
      lz1_256(2)%c(2)=-zcl(id1)*p1k0
* TL0 -- qu=p1,qd=p256,v=3,a=lz1_256(3)%a,c=lz1_256(3)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      eps_0=p1k0*p256(2)-p256k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p256(3)+p256k0*p1(3)
      lz1_256(3)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_256(3)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_256(3)%c(1)=zcr(id1)*ceps_1
      lz1_256(3)%c(2)=zcl(id1)*ceps_1
  
        if (ineutri(id1).ne.1) then
* TL0 -- qu=p1,qd=p256,v=0,a=lf1_256(0)%a,c=lf1_256(0)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      eps_0=-p1(2)*p256(3)+p256(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p256(0)+p256k0*p1(0)
      lf1_256(0)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_256(0)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_256(0)%c(1)=fcr(id1)*(p1(2)+ceps_1)
      lf1_256(0)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p256,v=1,a=lf1_256(1)%a,c=lf1_256(1)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      auxa=-quqd+p1k0*p256(1)+p256k0*p1(1)
      lf1_256(1)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_256(1)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_256(1)%c(1)=fcr(id1)*(p1(2)+ceps_1)
      lf1_256(1)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p256,v=2,a=lf1_256(2)%a,c=lf1_256(2)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      eps_0=-p1k0*p256(3)+p256k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p256(2)+p256k0*p1(2)
      lf1_256(2)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_256(2)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_256(2)%c(1)=fcr(id1)*p1k0
      lf1_256(2)%c(2)=-fcl(id1)*p1k0
* TL0 -- qu=p1,qd=p256,v=3,a=lf1_256(3)%a,c=lf1_256(3)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      eps_0=p1k0*p256(2)-p256k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p256(3)+p256k0*p1(3)
      lf1_256(3)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_256(3)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_256(3)%c(1)=fcr(id1)*ceps_1
      lf1_256(3)%c(2)=fcl(id1)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(1).eq.1) then
          do m=0,3
          do i=1,2
            lz1_278(m)%a(i)=czero
            lz1_278(m)%c(i)=czero
            lf1_278(m)%a(i)=czero
            lf1_278(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p1,q=p278
      quqd=p1(0)*p278(0)-p1(1)*p278(1)-p1(2)*p278(2)-p1(3)*p278(
     & 3)
* TL0 -- qu=p1,qd=p278,v=0,a=lz1_278(0)%a,c=lz1_278(0)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      eps_0=-p1(2)*p278(3)+p278(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p278(0)+p278k0*p1(0)
      lz1_278(0)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_278(0)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_278(0)%c(1)=zcr(id1)*(p1(2)+ceps_1)
      lz1_278(0)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p278,v=1,a=lz1_278(1)%a,c=lz1_278(1)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      auxa=-quqd+p1k0*p278(1)+p278k0*p1(1)
      lz1_278(1)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_278(1)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_278(1)%c(1)=zcr(id1)*(p1(2)+ceps_1)
      lz1_278(1)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p278,v=2,a=lz1_278(2)%a,c=lz1_278(2)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      eps_0=-p1k0*p278(3)+p278k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p278(2)+p278k0*p1(2)
      lz1_278(2)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_278(2)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_278(2)%c(1)=zcr(id1)*p1k0
      lz1_278(2)%c(2)=-zcl(id1)*p1k0
* TL0 -- qu=p1,qd=p278,v=3,a=lz1_278(3)%a,c=lz1_278(3)%c,cr=zcr(id1),cl=
* zcl(id1),nsum=0
      eps_0=p1k0*p278(2)-p278k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p278(3)+p278k0*p1(3)
      lz1_278(3)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_278(3)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_278(3)%c(1)=zcr(id1)*ceps_1
      lz1_278(3)%c(2)=zcl(id1)*ceps_1
  
        if (ineutri(id1).ne.1) then
* TL0 -- qu=p1,qd=p278,v=0,a=lf1_278(0)%a,c=lf1_278(0)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      eps_0=-p1(2)*p278(3)+p278(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p278(0)+p278k0*p1(0)
      lf1_278(0)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_278(0)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_278(0)%c(1)=fcr(id1)*(p1(2)+ceps_1)
      lf1_278(0)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p278,v=1,a=lf1_278(1)%a,c=lf1_278(1)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      auxa=-quqd+p1k0*p278(1)+p278k0*p1(1)
      lf1_278(1)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_278(1)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_278(1)%c(1)=fcr(id1)*(p1(2)+ceps_1)
      lf1_278(1)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p278,v=2,a=lf1_278(2)%a,c=lf1_278(2)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      eps_0=-p1k0*p278(3)+p278k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p278(2)+p278k0*p1(2)
      lf1_278(2)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_278(2)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_278(2)%c(1)=fcr(id1)*p1k0
      lf1_278(2)%c(2)=-fcl(id1)*p1k0
* TL0 -- qu=p1,qd=p278,v=3,a=lf1_278(3)%a,c=lf1_278(3)%c,cr=fcr(id1),cl=
* fcl(id1),nsum=0
      eps_0=p1k0*p278(2)-p278k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p278(3)+p278k0*p1(3)
      lf1_278(3)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_278(3)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_278(3)%c(1)=fcr(id1)*ceps_1
      lf1_278(3)%c(2)=fcl(id1)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(3).eq.1) then
          do m=0,3
          do i=1,2
            lz3_412(m)%a(i)=czero
            lz3_412(m)%c(i)=czero
            lf3_412(m)%a(i)=czero
            lf3_412(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p3,q=p412
      quqd=p3(0)*p412(0)-p3(1)*p412(1)-p3(2)*p412(2)-p3(3)*p412(
     & 3)
* TL0 -- qu=p3,qd=p412,v=0,a=lz3_412(0)%a,c=lz3_412(0)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lz3_412(0)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_412(0)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(0)%c(1)=zcr(id3)*(p3(2)+ceps_1)
      lz3_412(0)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=1,a=lz3_412(1)%a,c=lz3_412(1)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lz3_412(1)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_412(1)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(1)%c(1)=zcr(id3)*(p3(2)+ceps_1)
      lz3_412(1)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=2,a=lz3_412(2)%a,c=lz3_412(2)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lz3_412(2)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_412(2)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(2)%c(1)=zcr(id3)*p3k0
      lz3_412(2)%c(2)=-zcl(id3)*p3k0
* TL0 -- qu=p3,qd=p412,v=3,a=lz3_412(3)%a,c=lz3_412(3)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p412(3)+p412k0*p3(3)
      lz3_412(3)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_412(3)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(3)%c(1)=zcr(id3)*ceps_1
      lz3_412(3)%c(2)=zcl(id3)*ceps_1
  
        if (ineutri(id3).ne.1) then
* TL0 -- qu=p3,qd=p412,v=0,a=lf3_412(0)%a,c=lf3_412(0)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lf3_412(0)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_412(0)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(0)%c(1)=fcr(id3)*(p3(2)+ceps_1)
      lf3_412(0)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=1,a=lf3_412(1)%a,c=lf3_412(1)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lf3_412(1)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_412(1)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(1)%c(1)=fcr(id3)*(p3(2)+ceps_1)
      lf3_412(1)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=2,a=lf3_412(2)%a,c=lf3_412(2)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lf3_412(2)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_412(2)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(2)%c(1)=fcr(id3)*p3k0
      lf3_412(2)%c(2)=-fcl(id3)*p3k0
* TL0 -- qu=p3,qd=p412,v=3,a=lf3_412(3)%a,c=lf3_412(3)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p412(3)+p412k0*p3(3)
      lf3_412(3)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_412(3)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(3)%c(1)=fcr(id3)*ceps_1
      lf3_412(3)%c(2)=fcl(id3)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(3).eq.1) then
          do m=0,3
          do i=1,2
            lz3_456(m)%a(i)=czero
            lz3_456(m)%c(i)=czero
            lf3_456(m)%a(i)=czero
            lf3_456(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p3,q=p456
      quqd=p3(0)*p456(0)-p3(1)*p456(1)-p3(2)*p456(2)-p3(3)*p456(
     & 3)
* TL0 -- qu=p3,qd=p456,v=0,a=lz3_456(0)%a,c=lz3_456(0)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      eps_0=-p3(2)*p456(3)+p456(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p456(0)+p456k0*p3(0)
      lz3_456(0)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_456(0)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_456(0)%c(1)=zcr(id3)*(p3(2)+ceps_1)
      lz3_456(0)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p456,v=1,a=lz3_456(1)%a,c=lz3_456(1)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      auxa=-quqd+p3k0*p456(1)+p456k0*p3(1)
      lz3_456(1)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_456(1)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_456(1)%c(1)=zcr(id3)*(p3(2)+ceps_1)
      lz3_456(1)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p456,v=2,a=lz3_456(2)%a,c=lz3_456(2)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      eps_0=-p3k0*p456(3)+p456k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p456(2)+p456k0*p3(2)
      lz3_456(2)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_456(2)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_456(2)%c(1)=zcr(id3)*p3k0
      lz3_456(2)%c(2)=-zcl(id3)*p3k0
* TL0 -- qu=p3,qd=p456,v=3,a=lz3_456(3)%a,c=lz3_456(3)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      eps_0=p3k0*p456(2)-p456k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p456(3)+p456k0*p3(3)
      lz3_456(3)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_456(3)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_456(3)%c(1)=zcr(id3)*ceps_1
      lz3_456(3)%c(2)=zcl(id3)*ceps_1
  
        if (ineutri(id3).ne.1) then
* TL0 -- qu=p3,qd=p456,v=0,a=lf3_456(0)%a,c=lf3_456(0)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      eps_0=-p3(2)*p456(3)+p456(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p456(0)+p456k0*p3(0)
      lf3_456(0)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_456(0)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_456(0)%c(1)=fcr(id3)*(p3(2)+ceps_1)
      lf3_456(0)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p456,v=1,a=lf3_456(1)%a,c=lf3_456(1)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      auxa=-quqd+p3k0*p456(1)+p456k0*p3(1)
      lf3_456(1)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_456(1)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_456(1)%c(1)=fcr(id3)*(p3(2)+ceps_1)
      lf3_456(1)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p456,v=2,a=lf3_456(2)%a,c=lf3_456(2)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      eps_0=-p3k0*p456(3)+p456k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p456(2)+p456k0*p3(2)
      lf3_456(2)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_456(2)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_456(2)%c(1)=fcr(id3)*p3k0
      lf3_456(2)%c(2)=-fcl(id3)*p3k0
* TL0 -- qu=p3,qd=p456,v=3,a=lf3_456(3)%a,c=lf3_456(3)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      eps_0=p3k0*p456(2)-p456k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p456(3)+p456k0*p3(3)
      lf3_456(3)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_456(3)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_456(3)%c(1)=fcr(id3)*ceps_1
      lf3_456(3)%c(2)=fcl(id3)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(3).eq.1) then
          do m=0,3
          do i=1,2
            lz3_478(m)%a(i)=czero
            lz3_478(m)%c(i)=czero
            lf3_478(m)%a(i)=czero
            lf3_478(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p3,q=p478
      quqd=p3(0)*p478(0)-p3(1)*p478(1)-p3(2)*p478(2)-p3(3)*p478(
     & 3)
* TL0 -- qu=p3,qd=p478,v=0,a=lz3_478(0)%a,c=lz3_478(0)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      eps_0=-p3(2)*p478(3)+p478(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p478(0)+p478k0*p3(0)
      lz3_478(0)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_478(0)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_478(0)%c(1)=zcr(id3)*(p3(2)+ceps_1)
      lz3_478(0)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p478,v=1,a=lz3_478(1)%a,c=lz3_478(1)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      auxa=-quqd+p3k0*p478(1)+p478k0*p3(1)
      lz3_478(1)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_478(1)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_478(1)%c(1)=zcr(id3)*(p3(2)+ceps_1)
      lz3_478(1)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p478,v=2,a=lz3_478(2)%a,c=lz3_478(2)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      eps_0=-p3k0*p478(3)+p478k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p478(2)+p478k0*p3(2)
      lz3_478(2)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_478(2)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_478(2)%c(1)=zcr(id3)*p3k0
      lz3_478(2)%c(2)=-zcl(id3)*p3k0
* TL0 -- qu=p3,qd=p478,v=3,a=lz3_478(3)%a,c=lz3_478(3)%c,cr=zcr(id3),cl=
* zcl(id3),nsum=0
      eps_0=p3k0*p478(2)-p478k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p478(3)+p478k0*p3(3)
      lz3_478(3)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_478(3)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_478(3)%c(1)=zcr(id3)*ceps_1
      lz3_478(3)%c(2)=zcl(id3)*ceps_1
  
        if (ineutri(id3).ne.1) then
* TL0 -- qu=p3,qd=p478,v=0,a=lf3_478(0)%a,c=lf3_478(0)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      eps_0=-p3(2)*p478(3)+p478(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p478(0)+p478k0*p3(0)
      lf3_478(0)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_478(0)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_478(0)%c(1)=fcr(id3)*(p3(2)+ceps_1)
      lf3_478(0)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p478,v=1,a=lf3_478(1)%a,c=lf3_478(1)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      auxa=-quqd+p3k0*p478(1)+p478k0*p3(1)
      lf3_478(1)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_478(1)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_478(1)%c(1)=fcr(id3)*(p3(2)+ceps_1)
      lf3_478(1)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p478,v=2,a=lf3_478(2)%a,c=lf3_478(2)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      eps_0=-p3k0*p478(3)+p478k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p478(2)+p478k0*p3(2)
      lf3_478(2)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_478(2)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_478(2)%c(1)=fcr(id3)*p3k0
      lf3_478(2)%c(2)=-fcl(id3)*p3k0
* TL0 -- qu=p3,qd=p478,v=3,a=lf3_478(3)%a,c=lf3_478(3)%c,cr=fcr(id3),cl=
* fcl(id3),nsum=0
      eps_0=p3k0*p478(2)-p478k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p478(3)+p478k0*p3(3)
      lf3_478(3)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_478(3)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_478(3)%c(1)=fcr(id3)*ceps_1
      lf3_478(3)%c(2)=fcl(id3)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(5).eq.1) then
          do m=0,3
          do i=1,2
            lz5_612(m)%a(i)=czero
            lz5_612(m)%c(i)=czero
            lf5_612(m)%a(i)=czero
            lf5_612(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p5,q=p612
      quqd=p5(0)*p612(0)-p5(1)*p612(1)-p5(2)*p612(2)-p5(3)*p612(
     & 3)
* TL0 -- qu=p5,qd=p612,v=0,a=lz5_612(0)%a,c=lz5_612(0)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      eps_0=-p5(2)*p612(3)+p612(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p612(0)+p612k0*p5(0)
      lz5_612(0)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_612(0)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_612(0)%c(1)=zcr(id5)*(p5(2)+ceps_1)
      lz5_612(0)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p612,v=1,a=lz5_612(1)%a,c=lz5_612(1)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      auxa=-quqd+p5k0*p612(1)+p612k0*p5(1)
      lz5_612(1)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_612(1)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_612(1)%c(1)=zcr(id5)*(p5(2)+ceps_1)
      lz5_612(1)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p612,v=2,a=lz5_612(2)%a,c=lz5_612(2)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      eps_0=-p5k0*p612(3)+p612k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p612(2)+p612k0*p5(2)
      lz5_612(2)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_612(2)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_612(2)%c(1)=zcr(id5)*p5k0
      lz5_612(2)%c(2)=-zcl(id5)*p5k0
* TL0 -- qu=p5,qd=p612,v=3,a=lz5_612(3)%a,c=lz5_612(3)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      eps_0=p5k0*p612(2)-p612k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p612(3)+p612k0*p5(3)
      lz5_612(3)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_612(3)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_612(3)%c(1)=zcr(id5)*ceps_1
      lz5_612(3)%c(2)=zcl(id5)*ceps_1
  
        if (ineutri(id5).ne.1) then
* TL0 -- qu=p5,qd=p612,v=0,a=lf5_612(0)%a,c=lf5_612(0)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      eps_0=-p5(2)*p612(3)+p612(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p612(0)+p612k0*p5(0)
      lf5_612(0)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_612(0)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_612(0)%c(1)=fcr(id5)*(p5(2)+ceps_1)
      lf5_612(0)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p612,v=1,a=lf5_612(1)%a,c=lf5_612(1)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      auxa=-quqd+p5k0*p612(1)+p612k0*p5(1)
      lf5_612(1)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_612(1)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_612(1)%c(1)=fcr(id5)*(p5(2)+ceps_1)
      lf5_612(1)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p612,v=2,a=lf5_612(2)%a,c=lf5_612(2)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      eps_0=-p5k0*p612(3)+p612k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p612(2)+p612k0*p5(2)
      lf5_612(2)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_612(2)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_612(2)%c(1)=fcr(id5)*p5k0
      lf5_612(2)%c(2)=-fcl(id5)*p5k0
* TL0 -- qu=p5,qd=p612,v=3,a=lf5_612(3)%a,c=lf5_612(3)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      eps_0=p5k0*p612(2)-p612k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p612(3)+p612k0*p5(3)
      lf5_612(3)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_612(3)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_612(3)%c(1)=fcr(id5)*ceps_1
      lf5_612(3)%c(2)=fcl(id5)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(5).eq.1) then
          do m=0,3
          do i=1,2
            lz5_634(m)%a(i)=czero
            lz5_634(m)%c(i)=czero
            lf5_634(m)%a(i)=czero
            lf5_634(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p5,q=p634
      quqd=p5(0)*p634(0)-p5(1)*p634(1)-p5(2)*p634(2)-p5(3)*p634(
     & 3)
* TL0 -- qu=p5,qd=p634,v=0,a=lz5_634(0)%a,c=lz5_634(0)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      eps_0=-p5(2)*p634(3)+p634(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p634(0)+p634k0*p5(0)
      lz5_634(0)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_634(0)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_634(0)%c(1)=zcr(id5)*(p5(2)+ceps_1)
      lz5_634(0)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p634,v=1,a=lz5_634(1)%a,c=lz5_634(1)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      auxa=-quqd+p5k0*p634(1)+p634k0*p5(1)
      lz5_634(1)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_634(1)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_634(1)%c(1)=zcr(id5)*(p5(2)+ceps_1)
      lz5_634(1)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p634,v=2,a=lz5_634(2)%a,c=lz5_634(2)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      eps_0=-p5k0*p634(3)+p634k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p634(2)+p634k0*p5(2)
      lz5_634(2)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_634(2)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_634(2)%c(1)=zcr(id5)*p5k0
      lz5_634(2)%c(2)=-zcl(id5)*p5k0
* TL0 -- qu=p5,qd=p634,v=3,a=lz5_634(3)%a,c=lz5_634(3)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      eps_0=p5k0*p634(2)-p634k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p634(3)+p634k0*p5(3)
      lz5_634(3)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_634(3)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_634(3)%c(1)=zcr(id5)*ceps_1
      lz5_634(3)%c(2)=zcl(id5)*ceps_1
  
        if (ineutri(id5).ne.1) then
* TL0 -- qu=p5,qd=p634,v=0,a=lf5_634(0)%a,c=lf5_634(0)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      eps_0=-p5(2)*p634(3)+p634(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p634(0)+p634k0*p5(0)
      lf5_634(0)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_634(0)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_634(0)%c(1)=fcr(id5)*(p5(2)+ceps_1)
      lf5_634(0)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p634,v=1,a=lf5_634(1)%a,c=lf5_634(1)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      auxa=-quqd+p5k0*p634(1)+p634k0*p5(1)
      lf5_634(1)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_634(1)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_634(1)%c(1)=fcr(id5)*(p5(2)+ceps_1)
      lf5_634(1)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p634,v=2,a=lf5_634(2)%a,c=lf5_634(2)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      eps_0=-p5k0*p634(3)+p634k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p634(2)+p634k0*p5(2)
      lf5_634(2)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_634(2)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_634(2)%c(1)=fcr(id5)*p5k0
      lf5_634(2)%c(2)=-fcl(id5)*p5k0
* TL0 -- qu=p5,qd=p634,v=3,a=lf5_634(3)%a,c=lf5_634(3)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      eps_0=p5k0*p634(2)-p634k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p634(3)+p634k0*p5(3)
      lf5_634(3)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_634(3)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_634(3)%c(1)=fcr(id5)*ceps_1
      lf5_634(3)%c(2)=fcl(id5)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(5).eq.1) then
          do m=0,3
          do i=1,2
            lz5_678(m)%a(i)=czero
            lz5_678(m)%c(i)=czero
            lf5_678(m)%a(i)=czero
            lf5_678(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p5,q=p678
      quqd=p5(0)*p678(0)-p5(1)*p678(1)-p5(2)*p678(2)-p5(3)*p678(
     & 3)
* TL0 -- qu=p5,qd=p678,v=0,a=lz5_678(0)%a,c=lz5_678(0)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      eps_0=-p5(2)*p678(3)+p678(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p678(0)+p678k0*p5(0)
      lz5_678(0)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_678(0)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_678(0)%c(1)=zcr(id5)*(p5(2)+ceps_1)
      lz5_678(0)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p678,v=1,a=lz5_678(1)%a,c=lz5_678(1)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      auxa=-quqd+p5k0*p678(1)+p678k0*p5(1)
      lz5_678(1)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_678(1)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_678(1)%c(1)=zcr(id5)*(p5(2)+ceps_1)
      lz5_678(1)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p678,v=2,a=lz5_678(2)%a,c=lz5_678(2)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      eps_0=-p5k0*p678(3)+p678k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p678(2)+p678k0*p5(2)
      lz5_678(2)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_678(2)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_678(2)%c(1)=zcr(id5)*p5k0
      lz5_678(2)%c(2)=-zcl(id5)*p5k0
* TL0 -- qu=p5,qd=p678,v=3,a=lz5_678(3)%a,c=lz5_678(3)%c,cr=zcr(id5),cl=
* zcl(id5),nsum=0
      eps_0=p5k0*p678(2)-p678k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p678(3)+p678k0*p5(3)
      lz5_678(3)%a(1)=zcr(id5)*(auxa+ceps_0)
      lz5_678(3)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_678(3)%c(1)=zcr(id5)*ceps_1
      lz5_678(3)%c(2)=zcl(id5)*ceps_1
  
        if (ineutri(id5).ne.1) then
* TL0 -- qu=p5,qd=p678,v=0,a=lf5_678(0)%a,c=lf5_678(0)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      eps_0=-p5(2)*p678(3)+p678(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p678(0)+p678k0*p5(0)
      lf5_678(0)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_678(0)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_678(0)%c(1)=fcr(id5)*(p5(2)+ceps_1)
      lf5_678(0)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p678,v=1,a=lf5_678(1)%a,c=lf5_678(1)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      auxa=-quqd+p5k0*p678(1)+p678k0*p5(1)
      lf5_678(1)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_678(1)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_678(1)%c(1)=fcr(id5)*(p5(2)+ceps_1)
      lf5_678(1)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TL0 -- qu=p5,qd=p678,v=2,a=lf5_678(2)%a,c=lf5_678(2)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      eps_0=-p5k0*p678(3)+p678k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p678(2)+p678k0*p5(2)
      lf5_678(2)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_678(2)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_678(2)%c(1)=fcr(id5)*p5k0
      lf5_678(2)%c(2)=-fcl(id5)*p5k0
* TL0 -- qu=p5,qd=p678,v=3,a=lf5_678(3)%a,c=lf5_678(3)%c,cr=fcr(id5),cl=
* fcl(id5),nsum=0
      eps_0=p5k0*p678(2)-p678k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p678(3)+p678k0*p5(3)
      lf5_678(3)%a(1)=fcr(id5)*(auxa+ceps_0)
      lf5_678(3)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_678(3)%c(1)=fcr(id5)*ceps_1
      lf5_678(3)%c(2)=fcl(id5)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(7).eq.1) then
          do m=0,3
          do i=1,2
            lz7_812(m)%a(i)=czero
            lz7_812(m)%c(i)=czero
            lf7_812(m)%a(i)=czero
            lf7_812(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p7,q=p812
      quqd=p7(0)*p812(0)-p7(1)*p812(1)-p7(2)*p812(2)-p7(3)*p812(
     & 3)
* TL0 -- qu=p7,qd=p812,v=0,a=lz7_812(0)%a,c=lz7_812(0)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      eps_0=-p7(2)*p812(3)+p812(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p812(0)+p812k0*p7(0)
      lz7_812(0)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_812(0)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_812(0)%c(1)=zcr(id7)*(p7(2)+ceps_1)
      lz7_812(0)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p812,v=1,a=lz7_812(1)%a,c=lz7_812(1)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      auxa=-quqd+p7k0*p812(1)+p812k0*p7(1)
      lz7_812(1)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_812(1)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_812(1)%c(1)=zcr(id7)*(p7(2)+ceps_1)
      lz7_812(1)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p812,v=2,a=lz7_812(2)%a,c=lz7_812(2)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      eps_0=-p7k0*p812(3)+p812k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p812(2)+p812k0*p7(2)
      lz7_812(2)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_812(2)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_812(2)%c(1)=zcr(id7)*p7k0
      lz7_812(2)%c(2)=-zcl(id7)*p7k0
* TL0 -- qu=p7,qd=p812,v=3,a=lz7_812(3)%a,c=lz7_812(3)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      eps_0=p7k0*p812(2)-p812k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p812(3)+p812k0*p7(3)
      lz7_812(3)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_812(3)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_812(3)%c(1)=zcr(id7)*ceps_1
      lz7_812(3)%c(2)=zcl(id7)*ceps_1
  
        if (ineutri(id7).ne.1) then
* TL0 -- qu=p7,qd=p812,v=0,a=lf7_812(0)%a,c=lf7_812(0)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      eps_0=-p7(2)*p812(3)+p812(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p812(0)+p812k0*p7(0)
      lf7_812(0)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_812(0)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_812(0)%c(1)=fcr(id7)*(p7(2)+ceps_1)
      lf7_812(0)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p812,v=1,a=lf7_812(1)%a,c=lf7_812(1)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      auxa=-quqd+p7k0*p812(1)+p812k0*p7(1)
      lf7_812(1)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_812(1)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_812(1)%c(1)=fcr(id7)*(p7(2)+ceps_1)
      lf7_812(1)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p812,v=2,a=lf7_812(2)%a,c=lf7_812(2)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      eps_0=-p7k0*p812(3)+p812k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p812(2)+p812k0*p7(2)
      lf7_812(2)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_812(2)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_812(2)%c(1)=fcr(id7)*p7k0
      lf7_812(2)%c(2)=-fcl(id7)*p7k0
* TL0 -- qu=p7,qd=p812,v=3,a=lf7_812(3)%a,c=lf7_812(3)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      eps_0=p7k0*p812(2)-p812k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p812(3)+p812k0*p7(3)
      lf7_812(3)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_812(3)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_812(3)%c(1)=fcr(id7)*ceps_1
      lf7_812(3)%c(2)=fcl(id7)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(7).eq.1) then
          do m=0,3
          do i=1,2
            lz7_834(m)%a(i)=czero
            lz7_834(m)%c(i)=czero
            lf7_834(m)%a(i)=czero
            lf7_834(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p7,q=p834
      quqd=p7(0)*p834(0)-p7(1)*p834(1)-p7(2)*p834(2)-p7(3)*p834(
     & 3)
* TL0 -- qu=p7,qd=p834,v=0,a=lz7_834(0)%a,c=lz7_834(0)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      eps_0=-p7(2)*p834(3)+p834(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p834(0)+p834k0*p7(0)
      lz7_834(0)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_834(0)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_834(0)%c(1)=zcr(id7)*(p7(2)+ceps_1)
      lz7_834(0)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p834,v=1,a=lz7_834(1)%a,c=lz7_834(1)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      auxa=-quqd+p7k0*p834(1)+p834k0*p7(1)
      lz7_834(1)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_834(1)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_834(1)%c(1)=zcr(id7)*(p7(2)+ceps_1)
      lz7_834(1)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p834,v=2,a=lz7_834(2)%a,c=lz7_834(2)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      eps_0=-p7k0*p834(3)+p834k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p834(2)+p834k0*p7(2)
      lz7_834(2)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_834(2)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_834(2)%c(1)=zcr(id7)*p7k0
      lz7_834(2)%c(2)=-zcl(id7)*p7k0
* TL0 -- qu=p7,qd=p834,v=3,a=lz7_834(3)%a,c=lz7_834(3)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      eps_0=p7k0*p834(2)-p834k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p834(3)+p834k0*p7(3)
      lz7_834(3)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_834(3)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_834(3)%c(1)=zcr(id7)*ceps_1
      lz7_834(3)%c(2)=zcl(id7)*ceps_1
  
        if (ineutri(id7).ne.1) then
* TL0 -- qu=p7,qd=p834,v=0,a=lf7_834(0)%a,c=lf7_834(0)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      eps_0=-p7(2)*p834(3)+p834(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p834(0)+p834k0*p7(0)
      lf7_834(0)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_834(0)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_834(0)%c(1)=fcr(id7)*(p7(2)+ceps_1)
      lf7_834(0)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p834,v=1,a=lf7_834(1)%a,c=lf7_834(1)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      auxa=-quqd+p7k0*p834(1)+p834k0*p7(1)
      lf7_834(1)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_834(1)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_834(1)%c(1)=fcr(id7)*(p7(2)+ceps_1)
      lf7_834(1)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p834,v=2,a=lf7_834(2)%a,c=lf7_834(2)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      eps_0=-p7k0*p834(3)+p834k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p834(2)+p834k0*p7(2)
      lf7_834(2)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_834(2)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_834(2)%c(1)=fcr(id7)*p7k0
      lf7_834(2)%c(2)=-fcl(id7)*p7k0
* TL0 -- qu=p7,qd=p834,v=3,a=lf7_834(3)%a,c=lf7_834(3)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      eps_0=p7k0*p834(2)-p834k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p834(3)+p834k0*p7(3)
      lf7_834(3)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_834(3)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_834(3)%c(1)=fcr(id7)*ceps_1
      lf7_834(3)%c(2)=fcl(id7)*ceps_1
        endif
      endif
*zwend                                                                  
  
*zw                                                                     
  
  
      if(iresz(7).eq.1) then
          do m=0,3
          do i=1,2
            lz7_856(m)%a(i)=czero
            lz7_856(m)%c(i)=czero
            lf7_856(m)%a(i)=czero
            lf7_856(m)%c(i)=czero
          enddo
          enddo
      else
* quqd -- p=p7,q=p856
      quqd=p7(0)*p856(0)-p7(1)*p856(1)-p7(2)*p856(2)-p7(3)*p856(
     & 3)
* TL0 -- qu=p7,qd=p856,v=0,a=lz7_856(0)%a,c=lz7_856(0)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      eps_0=-p7(2)*p856(3)+p856(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p856(0)+p856k0*p7(0)
      lz7_856(0)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_856(0)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_856(0)%c(1)=zcr(id7)*(p7(2)+ceps_1)
      lz7_856(0)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p856,v=1,a=lz7_856(1)%a,c=lz7_856(1)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      auxa=-quqd+p7k0*p856(1)+p856k0*p7(1)
      lz7_856(1)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_856(1)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_856(1)%c(1)=zcr(id7)*(p7(2)+ceps_1)
      lz7_856(1)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p856,v=2,a=lz7_856(2)%a,c=lz7_856(2)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      eps_0=-p7k0*p856(3)+p856k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p856(2)+p856k0*p7(2)
      lz7_856(2)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_856(2)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_856(2)%c(1)=zcr(id7)*p7k0
      lz7_856(2)%c(2)=-zcl(id7)*p7k0
* TL0 -- qu=p7,qd=p856,v=3,a=lz7_856(3)%a,c=lz7_856(3)%c,cr=zcr(id7),cl=
* zcl(id7),nsum=0
      eps_0=p7k0*p856(2)-p856k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p856(3)+p856k0*p7(3)
      lz7_856(3)%a(1)=zcr(id7)*(auxa+ceps_0)
      lz7_856(3)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_856(3)%c(1)=zcr(id7)*ceps_1
      lz7_856(3)%c(2)=zcl(id7)*ceps_1
  
        if (ineutri(id7).ne.1) then
* TL0 -- qu=p7,qd=p856,v=0,a=lf7_856(0)%a,c=lf7_856(0)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      eps_0=-p7(2)*p856(3)+p856(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p856(0)+p856k0*p7(0)
      lf7_856(0)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_856(0)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_856(0)%c(1)=fcr(id7)*(p7(2)+ceps_1)
      lf7_856(0)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p856,v=1,a=lf7_856(1)%a,c=lf7_856(1)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      auxa=-quqd+p7k0*p856(1)+p856k0*p7(1)
      lf7_856(1)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_856(1)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_856(1)%c(1)=fcr(id7)*(p7(2)+ceps_1)
      lf7_856(1)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TL0 -- qu=p7,qd=p856,v=2,a=lf7_856(2)%a,c=lf7_856(2)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      eps_0=-p7k0*p856(3)+p856k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p856(2)+p856k0*p7(2)
      lf7_856(2)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_856(2)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_856(2)%c(1)=fcr(id7)*p7k0
      lf7_856(2)%c(2)=-fcl(id7)*p7k0
* TL0 -- qu=p7,qd=p856,v=3,a=lf7_856(3)%a,c=lf7_856(3)%c,cr=fcr(id7),cl=
* fcl(id7),nsum=0
      eps_0=p7k0*p856(2)-p856k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p856(3)+p856k0*p7(3)
      lf7_856(3)%a(1)=fcr(id7)*(auxa+ceps_0)
      lf7_856(3)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_856(3)%c(1)=fcr(id7)*ceps_1
      lf7_856(3)%c(2)=fcl(id7)*ceps_1
        endif
      endif
*zwend                                                                  
  
*                                                                       
* RIGHT VECTOR                                                          
*                                                                       
  
* quqd -- p=p134,q=p2
      quqd=p134(0)*p2(0)-p134(1)*p2(1)-p134(2)*p2(2)-p134(3)*p2(
     & 3)
* TR0 -- qu=p134,qd=p2,v=0,a=rz2_134(0)%a,b=rz2_134(0)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      eps_0=-p134(2)*p2(3)+p2(2)*p134(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p134k0*p2(0)+p2k0*p134(0)
      rz2_134(0)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_134(0)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(0)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_134(0)%b(2)=zcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p134,qd=p2,v=1,a=rz2_134(1)%a,b=rz2_134(1)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      auxa=-quqd+p134k0*p2(1)+p2k0*p134(1)
      rz2_134(1)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_134(1)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(1)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_134(1)%b(2)=zcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p134,qd=p2,v=2,a=rz2_134(2)%a,b=rz2_134(2)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      eps_0=-p134k0*p2(3)+p2k0*p134(3)
      ceps_0=eps_0*cim
      auxa=p134k0*p2(2)+p2k0*p134(2)
      rz2_134(2)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_134(2)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(2)%b(1)=-zcl(id2)*p2k0
      rz2_134(2)%b(2)=zcr(id2)*p2k0
* TR0 -- qu=p134,qd=p2,v=3,a=rz2_134(3)%a,b=rz2_134(3)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      eps_0=p134k0*p2(2)-p2k0*p134(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p134k0*p2(3)+p2k0*p134(3)
      rz2_134(3)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_134(3)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(3)%b(1)=-zcl(id2)*ceps_2
      rz2_134(3)%b(2)=-zcr(id2)*ceps_2
  
      if (ineutri(id1).ne.1) then
* TR0 -- qu=p134,qd=p2,v=0,a=rf2_134(0)%a,b=rf2_134(0)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      eps_0=-p134(2)*p2(3)+p2(2)*p134(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p134k0*p2(0)+p2k0*p134(0)
      rf2_134(0)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_134(0)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(0)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_134(0)%b(2)=fcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p134,qd=p2,v=1,a=rf2_134(1)%a,b=rf2_134(1)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      auxa=-quqd+p134k0*p2(1)+p2k0*p134(1)
      rf2_134(1)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_134(1)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(1)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_134(1)%b(2)=fcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p134,qd=p2,v=2,a=rf2_134(2)%a,b=rf2_134(2)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      eps_0=-p134k0*p2(3)+p2k0*p134(3)
      ceps_0=eps_0*cim
      auxa=p134k0*p2(2)+p2k0*p134(2)
      rf2_134(2)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_134(2)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(2)%b(1)=-fcl(id2)*p2k0
      rf2_134(2)%b(2)=fcr(id2)*p2k0
* TR0 -- qu=p134,qd=p2,v=3,a=rf2_134(3)%a,b=rf2_134(3)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      eps_0=p134k0*p2(2)-p2k0*p134(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p134k0*p2(3)+p2k0*p134(3)
      rf2_134(3)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_134(3)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(3)%b(1)=-fcl(id2)*ceps_2
      rf2_134(3)%b(2)=-fcr(id2)*ceps_2
      endif
* quqd -- p=p156,q=p2
      quqd=p156(0)*p2(0)-p156(1)*p2(1)-p156(2)*p2(2)-p156(3)*p2(
     & 3)
* TR0 -- qu=p156,qd=p2,v=0,a=rz2_156(0)%a,b=rz2_156(0)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      eps_0=-p156(2)*p2(3)+p2(2)*p156(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p156k0*p2(0)+p2k0*p156(0)
      rz2_156(0)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_156(0)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_156(0)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_156(0)%b(2)=zcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p156,qd=p2,v=1,a=rz2_156(1)%a,b=rz2_156(1)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      auxa=-quqd+p156k0*p2(1)+p2k0*p156(1)
      rz2_156(1)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_156(1)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_156(1)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_156(1)%b(2)=zcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p156,qd=p2,v=2,a=rz2_156(2)%a,b=rz2_156(2)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      eps_0=-p156k0*p2(3)+p2k0*p156(3)
      ceps_0=eps_0*cim
      auxa=p156k0*p2(2)+p2k0*p156(2)
      rz2_156(2)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_156(2)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_156(2)%b(1)=-zcl(id2)*p2k0
      rz2_156(2)%b(2)=zcr(id2)*p2k0
* TR0 -- qu=p156,qd=p2,v=3,a=rz2_156(3)%a,b=rz2_156(3)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      eps_0=p156k0*p2(2)-p2k0*p156(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p156k0*p2(3)+p2k0*p156(3)
      rz2_156(3)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_156(3)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_156(3)%b(1)=-zcl(id2)*ceps_2
      rz2_156(3)%b(2)=-zcr(id2)*ceps_2
  
      if (ineutri(id1).ne.1) then
* TR0 -- qu=p156,qd=p2,v=0,a=rf2_156(0)%a,b=rf2_156(0)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      eps_0=-p156(2)*p2(3)+p2(2)*p156(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p156k0*p2(0)+p2k0*p156(0)
      rf2_156(0)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_156(0)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_156(0)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_156(0)%b(2)=fcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p156,qd=p2,v=1,a=rf2_156(1)%a,b=rf2_156(1)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      auxa=-quqd+p156k0*p2(1)+p2k0*p156(1)
      rf2_156(1)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_156(1)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_156(1)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_156(1)%b(2)=fcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p156,qd=p2,v=2,a=rf2_156(2)%a,b=rf2_156(2)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      eps_0=-p156k0*p2(3)+p2k0*p156(3)
      ceps_0=eps_0*cim
      auxa=p156k0*p2(2)+p2k0*p156(2)
      rf2_156(2)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_156(2)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_156(2)%b(1)=-fcl(id2)*p2k0
      rf2_156(2)%b(2)=fcr(id2)*p2k0
* TR0 -- qu=p156,qd=p2,v=3,a=rf2_156(3)%a,b=rf2_156(3)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      eps_0=p156k0*p2(2)-p2k0*p156(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p156k0*p2(3)+p2k0*p156(3)
      rf2_156(3)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_156(3)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_156(3)%b(1)=-fcl(id2)*ceps_2
      rf2_156(3)%b(2)=-fcr(id2)*ceps_2
      endif
* quqd -- p=p178,q=p2
      quqd=p178(0)*p2(0)-p178(1)*p2(1)-p178(2)*p2(2)-p178(3)*p2(
     & 3)
* TR0 -- qu=p178,qd=p2,v=0,a=rz2_178(0)%a,b=rz2_178(0)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      eps_0=-p178(2)*p2(3)+p2(2)*p178(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p178k0*p2(0)+p2k0*p178(0)
      rz2_178(0)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_178(0)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_178(0)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_178(0)%b(2)=zcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p178,qd=p2,v=1,a=rz2_178(1)%a,b=rz2_178(1)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      auxa=-quqd+p178k0*p2(1)+p2k0*p178(1)
      rz2_178(1)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_178(1)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_178(1)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_178(1)%b(2)=zcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p178,qd=p2,v=2,a=rz2_178(2)%a,b=rz2_178(2)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      eps_0=-p178k0*p2(3)+p2k0*p178(3)
      ceps_0=eps_0*cim
      auxa=p178k0*p2(2)+p2k0*p178(2)
      rz2_178(2)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_178(2)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_178(2)%b(1)=-zcl(id2)*p2k0
      rz2_178(2)%b(2)=zcr(id2)*p2k0
* TR0 -- qu=p178,qd=p2,v=3,a=rz2_178(3)%a,b=rz2_178(3)%b,cr=zcr(id2),cl=
* zcl(id2),nsum=0
      eps_0=p178k0*p2(2)-p2k0*p178(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p178k0*p2(3)+p2k0*p178(3)
      rz2_178(3)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_178(3)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_178(3)%b(1)=-zcl(id2)*ceps_2
      rz2_178(3)%b(2)=-zcr(id2)*ceps_2
  
      if (ineutri(id1).ne.1) then
* TR0 -- qu=p178,qd=p2,v=0,a=rf2_178(0)%a,b=rf2_178(0)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      eps_0=-p178(2)*p2(3)+p2(2)*p178(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p178k0*p2(0)+p2k0*p178(0)
      rf2_178(0)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_178(0)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_178(0)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_178(0)%b(2)=fcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p178,qd=p2,v=1,a=rf2_178(1)%a,b=rf2_178(1)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      auxa=-quqd+p178k0*p2(1)+p2k0*p178(1)
      rf2_178(1)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_178(1)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_178(1)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_178(1)%b(2)=fcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p178,qd=p2,v=2,a=rf2_178(2)%a,b=rf2_178(2)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      eps_0=-p178k0*p2(3)+p2k0*p178(3)
      ceps_0=eps_0*cim
      auxa=p178k0*p2(2)+p2k0*p178(2)
      rf2_178(2)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_178(2)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_178(2)%b(1)=-fcl(id2)*p2k0
      rf2_178(2)%b(2)=fcr(id2)*p2k0
* TR0 -- qu=p178,qd=p2,v=3,a=rf2_178(3)%a,b=rf2_178(3)%b,cr=fcr(id2),cl=
* fcl(id2),nsum=0
      eps_0=p178k0*p2(2)-p2k0*p178(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p178k0*p2(3)+p2k0*p178(3)
      rf2_178(3)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_178(3)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_178(3)%b(1)=-fcl(id2)*ceps_2
      rf2_178(3)%b(2)=-fcr(id2)*ceps_2
      endif
* quqd -- p=p312,q=p4
      quqd=p312(0)*p4(0)-p312(1)*p4(1)-p312(2)*p4(2)-p312(3)*p4(
     & 3)
* TR0 -- qu=p312,qd=p4,v=0,a=rz4_312(0)%a,b=rz4_312(0)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rz4_312(0)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_312(0)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(0)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_312(0)%b(2)=zcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=1,a=rz4_312(1)%a,b=rz4_312(1)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rz4_312(1)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_312(1)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(1)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_312(1)%b(2)=zcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=2,a=rz4_312(2)%a,b=rz4_312(2)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rz4_312(2)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_312(2)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(2)%b(1)=-zcl(id4)*p4k0
      rz4_312(2)%b(2)=zcr(id4)*p4k0
* TR0 -- qu=p312,qd=p4,v=3,a=rz4_312(3)%a,b=rz4_312(3)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p312k0*p4(3)+p4k0*p312(3)
      rz4_312(3)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_312(3)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(3)%b(1)=-zcl(id4)*ceps_2
      rz4_312(3)%b(2)=-zcr(id4)*ceps_2
  
      if (ineutri(id3).ne.1) then
* TR0 -- qu=p312,qd=p4,v=0,a=rf4_312(0)%a,b=rf4_312(0)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rf4_312(0)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_312(0)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(0)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_312(0)%b(2)=fcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=1,a=rf4_312(1)%a,b=rf4_312(1)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rf4_312(1)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_312(1)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(1)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_312(1)%b(2)=fcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=2,a=rf4_312(2)%a,b=rf4_312(2)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rf4_312(2)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_312(2)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(2)%b(1)=-fcl(id4)*p4k0
      rf4_312(2)%b(2)=fcr(id4)*p4k0
* TR0 -- qu=p312,qd=p4,v=3,a=rf4_312(3)%a,b=rf4_312(3)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p312k0*p4(3)+p4k0*p312(3)
      rf4_312(3)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_312(3)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(3)%b(1)=-fcl(id4)*ceps_2
      rf4_312(3)%b(2)=-fcr(id4)*ceps_2
      endif
* quqd -- p=p356,q=p4
      quqd=p356(0)*p4(0)-p356(1)*p4(1)-p356(2)*p4(2)-p356(3)*p4(
     & 3)
* TR0 -- qu=p356,qd=p4,v=0,a=rz4_356(0)%a,b=rz4_356(0)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      eps_0=-p356(2)*p4(3)+p4(2)*p356(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p356k0*p4(0)+p4k0*p356(0)
      rz4_356(0)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_356(0)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_356(0)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_356(0)%b(2)=zcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p356,qd=p4,v=1,a=rz4_356(1)%a,b=rz4_356(1)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      auxa=-quqd+p356k0*p4(1)+p4k0*p356(1)
      rz4_356(1)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_356(1)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_356(1)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_356(1)%b(2)=zcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p356,qd=p4,v=2,a=rz4_356(2)%a,b=rz4_356(2)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      eps_0=-p356k0*p4(3)+p4k0*p356(3)
      ceps_0=eps_0*cim
      auxa=p356k0*p4(2)+p4k0*p356(2)
      rz4_356(2)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_356(2)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_356(2)%b(1)=-zcl(id4)*p4k0
      rz4_356(2)%b(2)=zcr(id4)*p4k0
* TR0 -- qu=p356,qd=p4,v=3,a=rz4_356(3)%a,b=rz4_356(3)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      eps_0=p356k0*p4(2)-p4k0*p356(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p356k0*p4(3)+p4k0*p356(3)
      rz4_356(3)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_356(3)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_356(3)%b(1)=-zcl(id4)*ceps_2
      rz4_356(3)%b(2)=-zcr(id4)*ceps_2
  
      if (ineutri(id3).ne.1) then
* TR0 -- qu=p356,qd=p4,v=0,a=rf4_356(0)%a,b=rf4_356(0)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      eps_0=-p356(2)*p4(3)+p4(2)*p356(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p356k0*p4(0)+p4k0*p356(0)
      rf4_356(0)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_356(0)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_356(0)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_356(0)%b(2)=fcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p356,qd=p4,v=1,a=rf4_356(1)%a,b=rf4_356(1)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      auxa=-quqd+p356k0*p4(1)+p4k0*p356(1)
      rf4_356(1)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_356(1)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_356(1)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_356(1)%b(2)=fcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p356,qd=p4,v=2,a=rf4_356(2)%a,b=rf4_356(2)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      eps_0=-p356k0*p4(3)+p4k0*p356(3)
      ceps_0=eps_0*cim
      auxa=p356k0*p4(2)+p4k0*p356(2)
      rf4_356(2)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_356(2)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_356(2)%b(1)=-fcl(id4)*p4k0
      rf4_356(2)%b(2)=fcr(id4)*p4k0
* TR0 -- qu=p356,qd=p4,v=3,a=rf4_356(3)%a,b=rf4_356(3)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      eps_0=p356k0*p4(2)-p4k0*p356(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p356k0*p4(3)+p4k0*p356(3)
      rf4_356(3)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_356(3)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_356(3)%b(1)=-fcl(id4)*ceps_2
      rf4_356(3)%b(2)=-fcr(id4)*ceps_2
      endif
* quqd -- p=p378,q=p4
      quqd=p378(0)*p4(0)-p378(1)*p4(1)-p378(2)*p4(2)-p378(3)*p4(
     & 3)
* TR0 -- qu=p378,qd=p4,v=0,a=rz4_378(0)%a,b=rz4_378(0)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      eps_0=-p378(2)*p4(3)+p4(2)*p378(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p378k0*p4(0)+p4k0*p378(0)
      rz4_378(0)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_378(0)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_378(0)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_378(0)%b(2)=zcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p378,qd=p4,v=1,a=rz4_378(1)%a,b=rz4_378(1)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      auxa=-quqd+p378k0*p4(1)+p4k0*p378(1)
      rz4_378(1)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_378(1)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_378(1)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_378(1)%b(2)=zcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p378,qd=p4,v=2,a=rz4_378(2)%a,b=rz4_378(2)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      eps_0=-p378k0*p4(3)+p4k0*p378(3)
      ceps_0=eps_0*cim
      auxa=p378k0*p4(2)+p4k0*p378(2)
      rz4_378(2)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_378(2)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_378(2)%b(1)=-zcl(id4)*p4k0
      rz4_378(2)%b(2)=zcr(id4)*p4k0
* TR0 -- qu=p378,qd=p4,v=3,a=rz4_378(3)%a,b=rz4_378(3)%b,cr=zcr(id4),cl=
* zcl(id4),nsum=0
      eps_0=p378k0*p4(2)-p4k0*p378(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p378k0*p4(3)+p4k0*p378(3)
      rz4_378(3)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_378(3)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_378(3)%b(1)=-zcl(id4)*ceps_2
      rz4_378(3)%b(2)=-zcr(id4)*ceps_2
  
      if (ineutri(id3).ne.1) then
* TR0 -- qu=p378,qd=p4,v=0,a=rf4_378(0)%a,b=rf4_378(0)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      eps_0=-p378(2)*p4(3)+p4(2)*p378(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p378k0*p4(0)+p4k0*p378(0)
      rf4_378(0)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_378(0)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_378(0)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_378(0)%b(2)=fcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p378,qd=p4,v=1,a=rf4_378(1)%a,b=rf4_378(1)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      auxa=-quqd+p378k0*p4(1)+p4k0*p378(1)
      rf4_378(1)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_378(1)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_378(1)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_378(1)%b(2)=fcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p378,qd=p4,v=2,a=rf4_378(2)%a,b=rf4_378(2)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      eps_0=-p378k0*p4(3)+p4k0*p378(3)
      ceps_0=eps_0*cim
      auxa=p378k0*p4(2)+p4k0*p378(2)
      rf4_378(2)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_378(2)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_378(2)%b(1)=-fcl(id4)*p4k0
      rf4_378(2)%b(2)=fcr(id4)*p4k0
* TR0 -- qu=p378,qd=p4,v=3,a=rf4_378(3)%a,b=rf4_378(3)%b,cr=fcr(id4),cl=
* fcl(id4),nsum=0
      eps_0=p378k0*p4(2)-p4k0*p378(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p378k0*p4(3)+p4k0*p378(3)
      rf4_378(3)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_378(3)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_378(3)%b(1)=-fcl(id4)*ceps_2
      rf4_378(3)%b(2)=-fcr(id4)*ceps_2
      endif
* quqd -- p=p512,q=p6
      quqd=p512(0)*p6(0)-p512(1)*p6(1)-p512(2)*p6(2)-p512(3)*p6(
     & 3)
* TR0 -- qu=p512,qd=p6,v=0,a=rz6_512(0)%a,b=rz6_512(0)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      eps_0=-p512(2)*p6(3)+p6(2)*p512(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p512k0*p6(0)+p6k0*p512(0)
      rz6_512(0)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_512(0)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_512(0)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
      rz6_512(0)%b(2)=zcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p512,qd=p6,v=1,a=rz6_512(1)%a,b=rz6_512(1)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      auxa=-quqd+p512k0*p6(1)+p6k0*p512(1)
      rz6_512(1)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_512(1)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_512(1)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
      rz6_512(1)%b(2)=zcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p512,qd=p6,v=2,a=rz6_512(2)%a,b=rz6_512(2)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      eps_0=-p512k0*p6(3)+p6k0*p512(3)
      ceps_0=eps_0*cim
      auxa=p512k0*p6(2)+p6k0*p512(2)
      rz6_512(2)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_512(2)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_512(2)%b(1)=-zcl(id6)*p6k0
      rz6_512(2)%b(2)=zcr(id6)*p6k0
* TR0 -- qu=p512,qd=p6,v=3,a=rz6_512(3)%a,b=rz6_512(3)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      eps_0=p512k0*p6(2)-p6k0*p512(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p512k0*p6(3)+p6k0*p512(3)
      rz6_512(3)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_512(3)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_512(3)%b(1)=-zcl(id6)*ceps_2
      rz6_512(3)%b(2)=-zcr(id6)*ceps_2
  
      if (ineutri(id5).ne.1) then
* TR0 -- qu=p512,qd=p6,v=0,a=rf6_512(0)%a,b=rf6_512(0)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      eps_0=-p512(2)*p6(3)+p6(2)*p512(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p512k0*p6(0)+p6k0*p512(0)
      rf6_512(0)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_512(0)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_512(0)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
      rf6_512(0)%b(2)=fcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p512,qd=p6,v=1,a=rf6_512(1)%a,b=rf6_512(1)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      auxa=-quqd+p512k0*p6(1)+p6k0*p512(1)
      rf6_512(1)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_512(1)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_512(1)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
      rf6_512(1)%b(2)=fcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p512,qd=p6,v=2,a=rf6_512(2)%a,b=rf6_512(2)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      eps_0=-p512k0*p6(3)+p6k0*p512(3)
      ceps_0=eps_0*cim
      auxa=p512k0*p6(2)+p6k0*p512(2)
      rf6_512(2)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_512(2)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_512(2)%b(1)=-fcl(id6)*p6k0
      rf6_512(2)%b(2)=fcr(id6)*p6k0
* TR0 -- qu=p512,qd=p6,v=3,a=rf6_512(3)%a,b=rf6_512(3)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      eps_0=p512k0*p6(2)-p6k0*p512(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p512k0*p6(3)+p6k0*p512(3)
      rf6_512(3)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_512(3)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_512(3)%b(1)=-fcl(id6)*ceps_2
      rf6_512(3)%b(2)=-fcr(id6)*ceps_2
      endif
* quqd -- p=p534,q=p6
      quqd=p534(0)*p6(0)-p534(1)*p6(1)-p534(2)*p6(2)-p534(3)*p6(
     & 3)
* TR0 -- qu=p534,qd=p6,v=0,a=rz6_534(0)%a,b=rz6_534(0)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      eps_0=-p534(2)*p6(3)+p6(2)*p534(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p534k0*p6(0)+p6k0*p534(0)
      rz6_534(0)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_534(0)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_534(0)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
      rz6_534(0)%b(2)=zcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p534,qd=p6,v=1,a=rz6_534(1)%a,b=rz6_534(1)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      auxa=-quqd+p534k0*p6(1)+p6k0*p534(1)
      rz6_534(1)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_534(1)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_534(1)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
      rz6_534(1)%b(2)=zcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p534,qd=p6,v=2,a=rz6_534(2)%a,b=rz6_534(2)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      eps_0=-p534k0*p6(3)+p6k0*p534(3)
      ceps_0=eps_0*cim
      auxa=p534k0*p6(2)+p6k0*p534(2)
      rz6_534(2)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_534(2)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_534(2)%b(1)=-zcl(id6)*p6k0
      rz6_534(2)%b(2)=zcr(id6)*p6k0
* TR0 -- qu=p534,qd=p6,v=3,a=rz6_534(3)%a,b=rz6_534(3)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      eps_0=p534k0*p6(2)-p6k0*p534(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p534k0*p6(3)+p6k0*p534(3)
      rz6_534(3)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_534(3)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_534(3)%b(1)=-zcl(id6)*ceps_2
      rz6_534(3)%b(2)=-zcr(id6)*ceps_2
  
      if (ineutri(id5).ne.1) then
* TR0 -- qu=p534,qd=p6,v=0,a=rf6_534(0)%a,b=rf6_534(0)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      eps_0=-p534(2)*p6(3)+p6(2)*p534(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p534k0*p6(0)+p6k0*p534(0)
      rf6_534(0)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_534(0)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_534(0)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
      rf6_534(0)%b(2)=fcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p534,qd=p6,v=1,a=rf6_534(1)%a,b=rf6_534(1)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      auxa=-quqd+p534k0*p6(1)+p6k0*p534(1)
      rf6_534(1)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_534(1)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_534(1)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
      rf6_534(1)%b(2)=fcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p534,qd=p6,v=2,a=rf6_534(2)%a,b=rf6_534(2)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      eps_0=-p534k0*p6(3)+p6k0*p534(3)
      ceps_0=eps_0*cim
      auxa=p534k0*p6(2)+p6k0*p534(2)
      rf6_534(2)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_534(2)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_534(2)%b(1)=-fcl(id6)*p6k0
      rf6_534(2)%b(2)=fcr(id6)*p6k0
* TR0 -- qu=p534,qd=p6,v=3,a=rf6_534(3)%a,b=rf6_534(3)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      eps_0=p534k0*p6(2)-p6k0*p534(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p534k0*p6(3)+p6k0*p534(3)
      rf6_534(3)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_534(3)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_534(3)%b(1)=-fcl(id6)*ceps_2
      rf6_534(3)%b(2)=-fcr(id6)*ceps_2
      endif
* quqd -- p=p578,q=p6
      quqd=p578(0)*p6(0)-p578(1)*p6(1)-p578(2)*p6(2)-p578(3)*p6(
     & 3)
* TR0 -- qu=p578,qd=p6,v=0,a=rz6_578(0)%a,b=rz6_578(0)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      eps_0=-p578(2)*p6(3)+p6(2)*p578(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p578k0*p6(0)+p6k0*p578(0)
      rz6_578(0)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_578(0)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_578(0)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
      rz6_578(0)%b(2)=zcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p578,qd=p6,v=1,a=rz6_578(1)%a,b=rz6_578(1)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      auxa=-quqd+p578k0*p6(1)+p6k0*p578(1)
      rz6_578(1)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_578(1)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_578(1)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
      rz6_578(1)%b(2)=zcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p578,qd=p6,v=2,a=rz6_578(2)%a,b=rz6_578(2)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      eps_0=-p578k0*p6(3)+p6k0*p578(3)
      ceps_0=eps_0*cim
      auxa=p578k0*p6(2)+p6k0*p578(2)
      rz6_578(2)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_578(2)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_578(2)%b(1)=-zcl(id6)*p6k0
      rz6_578(2)%b(2)=zcr(id6)*p6k0
* TR0 -- qu=p578,qd=p6,v=3,a=rz6_578(3)%a,b=rz6_578(3)%b,cr=zcr(id6),cl=
* zcl(id6),nsum=0
      eps_0=p578k0*p6(2)-p6k0*p578(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p578k0*p6(3)+p6k0*p578(3)
      rz6_578(3)%a(1)=zcr(id6)*(auxa+ceps_0)
      rz6_578(3)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_578(3)%b(1)=-zcl(id6)*ceps_2
      rz6_578(3)%b(2)=-zcr(id6)*ceps_2
  
      if (ineutri(id5).ne.1) then
* TR0 -- qu=p578,qd=p6,v=0,a=rf6_578(0)%a,b=rf6_578(0)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      eps_0=-p578(2)*p6(3)+p6(2)*p578(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p578k0*p6(0)+p6k0*p578(0)
      rf6_578(0)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_578(0)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_578(0)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
      rf6_578(0)%b(2)=fcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p578,qd=p6,v=1,a=rf6_578(1)%a,b=rf6_578(1)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      auxa=-quqd+p578k0*p6(1)+p6k0*p578(1)
      rf6_578(1)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_578(1)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_578(1)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
      rf6_578(1)%b(2)=fcr(id6)*(p6(2)-ceps_2)
* TR0 -- qu=p578,qd=p6,v=2,a=rf6_578(2)%a,b=rf6_578(2)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      eps_0=-p578k0*p6(3)+p6k0*p578(3)
      ceps_0=eps_0*cim
      auxa=p578k0*p6(2)+p6k0*p578(2)
      rf6_578(2)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_578(2)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_578(2)%b(1)=-fcl(id6)*p6k0
      rf6_578(2)%b(2)=fcr(id6)*p6k0
* TR0 -- qu=p578,qd=p6,v=3,a=rf6_578(3)%a,b=rf6_578(3)%b,cr=fcr(id6),cl=
* fcl(id6),nsum=0
      eps_0=p578k0*p6(2)-p6k0*p578(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p578k0*p6(3)+p6k0*p578(3)
      rf6_578(3)%a(1)=fcr(id6)*(auxa+ceps_0)
      rf6_578(3)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_578(3)%b(1)=-fcl(id6)*ceps_2
      rf6_578(3)%b(2)=-fcr(id6)*ceps_2
      endif
* quqd -- p=p712,q=p8
      quqd=p712(0)*p8(0)-p712(1)*p8(1)-p712(2)*p8(2)-p712(3)*p8(
     & 3)
* TR0 -- qu=p712,qd=p8,v=0,a=rz8_712(0)%a,b=rz8_712(0)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      eps_0=-p712(2)*p8(3)+p8(2)*p712(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p712k0*p8(0)+p8k0*p712(0)
      rz8_712(0)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_712(0)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_712(0)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
      rz8_712(0)%b(2)=zcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p712,qd=p8,v=1,a=rz8_712(1)%a,b=rz8_712(1)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      auxa=-quqd+p712k0*p8(1)+p8k0*p712(1)
      rz8_712(1)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_712(1)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_712(1)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
      rz8_712(1)%b(2)=zcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p712,qd=p8,v=2,a=rz8_712(2)%a,b=rz8_712(2)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      eps_0=-p712k0*p8(3)+p8k0*p712(3)
      ceps_0=eps_0*cim
      auxa=p712k0*p8(2)+p8k0*p712(2)
      rz8_712(2)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_712(2)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_712(2)%b(1)=-zcl(id8)*p8k0
      rz8_712(2)%b(2)=zcr(id8)*p8k0
* TR0 -- qu=p712,qd=p8,v=3,a=rz8_712(3)%a,b=rz8_712(3)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      eps_0=p712k0*p8(2)-p8k0*p712(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p712k0*p8(3)+p8k0*p712(3)
      rz8_712(3)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_712(3)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_712(3)%b(1)=-zcl(id8)*ceps_2
      rz8_712(3)%b(2)=-zcr(id8)*ceps_2
  
      if (ineutri(id7).ne.1) then
* TR0 -- qu=p712,qd=p8,v=0,a=rf8_712(0)%a,b=rf8_712(0)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      eps_0=-p712(2)*p8(3)+p8(2)*p712(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p712k0*p8(0)+p8k0*p712(0)
      rf8_712(0)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_712(0)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_712(0)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
      rf8_712(0)%b(2)=fcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p712,qd=p8,v=1,a=rf8_712(1)%a,b=rf8_712(1)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      auxa=-quqd+p712k0*p8(1)+p8k0*p712(1)
      rf8_712(1)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_712(1)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_712(1)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
      rf8_712(1)%b(2)=fcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p712,qd=p8,v=2,a=rf8_712(2)%a,b=rf8_712(2)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      eps_0=-p712k0*p8(3)+p8k0*p712(3)
      ceps_0=eps_0*cim
      auxa=p712k0*p8(2)+p8k0*p712(2)
      rf8_712(2)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_712(2)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_712(2)%b(1)=-fcl(id8)*p8k0
      rf8_712(2)%b(2)=fcr(id8)*p8k0
* TR0 -- qu=p712,qd=p8,v=3,a=rf8_712(3)%a,b=rf8_712(3)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      eps_0=p712k0*p8(2)-p8k0*p712(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p712k0*p8(3)+p8k0*p712(3)
      rf8_712(3)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_712(3)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_712(3)%b(1)=-fcl(id8)*ceps_2
      rf8_712(3)%b(2)=-fcr(id8)*ceps_2
      endif
* quqd -- p=p734,q=p8
      quqd=p734(0)*p8(0)-p734(1)*p8(1)-p734(2)*p8(2)-p734(3)*p8(
     & 3)
* TR0 -- qu=p734,qd=p8,v=0,a=rz8_734(0)%a,b=rz8_734(0)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      eps_0=-p734(2)*p8(3)+p8(2)*p734(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p734k0*p8(0)+p8k0*p734(0)
      rz8_734(0)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_734(0)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_734(0)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
      rz8_734(0)%b(2)=zcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p734,qd=p8,v=1,a=rz8_734(1)%a,b=rz8_734(1)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      auxa=-quqd+p734k0*p8(1)+p8k0*p734(1)
      rz8_734(1)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_734(1)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_734(1)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
      rz8_734(1)%b(2)=zcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p734,qd=p8,v=2,a=rz8_734(2)%a,b=rz8_734(2)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      eps_0=-p734k0*p8(3)+p8k0*p734(3)
      ceps_0=eps_0*cim
      auxa=p734k0*p8(2)+p8k0*p734(2)
      rz8_734(2)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_734(2)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_734(2)%b(1)=-zcl(id8)*p8k0
      rz8_734(2)%b(2)=zcr(id8)*p8k0
* TR0 -- qu=p734,qd=p8,v=3,a=rz8_734(3)%a,b=rz8_734(3)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      eps_0=p734k0*p8(2)-p8k0*p734(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p734k0*p8(3)+p8k0*p734(3)
      rz8_734(3)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_734(3)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_734(3)%b(1)=-zcl(id8)*ceps_2
      rz8_734(3)%b(2)=-zcr(id8)*ceps_2
  
      if (ineutri(id7).ne.1) then
* TR0 -- qu=p734,qd=p8,v=0,a=rf8_734(0)%a,b=rf8_734(0)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      eps_0=-p734(2)*p8(3)+p8(2)*p734(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p734k0*p8(0)+p8k0*p734(0)
      rf8_734(0)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_734(0)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_734(0)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
      rf8_734(0)%b(2)=fcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p734,qd=p8,v=1,a=rf8_734(1)%a,b=rf8_734(1)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      auxa=-quqd+p734k0*p8(1)+p8k0*p734(1)
      rf8_734(1)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_734(1)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_734(1)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
      rf8_734(1)%b(2)=fcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p734,qd=p8,v=2,a=rf8_734(2)%a,b=rf8_734(2)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      eps_0=-p734k0*p8(3)+p8k0*p734(3)
      ceps_0=eps_0*cim
      auxa=p734k0*p8(2)+p8k0*p734(2)
      rf8_734(2)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_734(2)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_734(2)%b(1)=-fcl(id8)*p8k0
      rf8_734(2)%b(2)=fcr(id8)*p8k0
* TR0 -- qu=p734,qd=p8,v=3,a=rf8_734(3)%a,b=rf8_734(3)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      eps_0=p734k0*p8(2)-p8k0*p734(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p734k0*p8(3)+p8k0*p734(3)
      rf8_734(3)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_734(3)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_734(3)%b(1)=-fcl(id8)*ceps_2
      rf8_734(3)%b(2)=-fcr(id8)*ceps_2
      endif
* quqd -- p=p756,q=p8
      quqd=p756(0)*p8(0)-p756(1)*p8(1)-p756(2)*p8(2)-p756(3)*p8(
     & 3)
* TR0 -- qu=p756,qd=p8,v=0,a=rz8_756(0)%a,b=rz8_756(0)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      eps_0=-p756(2)*p8(3)+p8(2)*p756(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p756k0*p8(0)+p8k0*p756(0)
      rz8_756(0)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_756(0)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_756(0)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
      rz8_756(0)%b(2)=zcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p756,qd=p8,v=1,a=rz8_756(1)%a,b=rz8_756(1)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      auxa=-quqd+p756k0*p8(1)+p8k0*p756(1)
      rz8_756(1)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_756(1)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_756(1)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
      rz8_756(1)%b(2)=zcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p756,qd=p8,v=2,a=rz8_756(2)%a,b=rz8_756(2)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      eps_0=-p756k0*p8(3)+p8k0*p756(3)
      ceps_0=eps_0*cim
      auxa=p756k0*p8(2)+p8k0*p756(2)
      rz8_756(2)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_756(2)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_756(2)%b(1)=-zcl(id8)*p8k0
      rz8_756(2)%b(2)=zcr(id8)*p8k0
* TR0 -- qu=p756,qd=p8,v=3,a=rz8_756(3)%a,b=rz8_756(3)%b,cr=zcr(id8),cl=
* zcl(id8),nsum=0
      eps_0=p756k0*p8(2)-p8k0*p756(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p756k0*p8(3)+p8k0*p756(3)
      rz8_756(3)%a(1)=zcr(id8)*(auxa+ceps_0)
      rz8_756(3)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_756(3)%b(1)=-zcl(id8)*ceps_2
      rz8_756(3)%b(2)=-zcr(id8)*ceps_2
  
      if (ineutri(id7).ne.1) then
* TR0 -- qu=p756,qd=p8,v=0,a=rf8_756(0)%a,b=rf8_756(0)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      eps_0=-p756(2)*p8(3)+p8(2)*p756(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p756k0*p8(0)+p8k0*p756(0)
      rf8_756(0)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_756(0)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_756(0)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
      rf8_756(0)%b(2)=fcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p756,qd=p8,v=1,a=rf8_756(1)%a,b=rf8_756(1)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      auxa=-quqd+p756k0*p8(1)+p8k0*p756(1)
      rf8_756(1)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_756(1)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_756(1)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
      rf8_756(1)%b(2)=fcr(id8)*(p8(2)-ceps_2)
* TR0 -- qu=p756,qd=p8,v=2,a=rf8_756(2)%a,b=rf8_756(2)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      eps_0=-p756k0*p8(3)+p8k0*p756(3)
      ceps_0=eps_0*cim
      auxa=p756k0*p8(2)+p8k0*p756(2)
      rf8_756(2)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_756(2)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_756(2)%b(1)=-fcl(id8)*p8k0
      rf8_756(2)%b(2)=fcr(id8)*p8k0
* TR0 -- qu=p756,qd=p8,v=3,a=rf8_756(3)%a,b=rf8_756(3)%b,cr=fcr(id8),cl=
* fcl(id8),nsum=0
      eps_0=p756k0*p8(2)-p8k0*p756(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p756k0*p8(3)+p8k0*p756(3)
      rf8_756(3)%a(1)=fcr(id8)*(auxa+ceps_0)
      rf8_756(3)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_756(3)%b(1)=-fcl(id8)*ceps_2
      rf8_756(3)%b(2)=-fcr(id8)*ceps_2
      endif
  
*                                                                       
* LEFT*RIGHT VECTOR                                                     
*                                                                       
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz1234(&,i3)%e(mu),a1=lz1_234(mu)%a,c1=lz1_234(mu)%c,a2=r2
* _34(i3)%a,b2=r2_34(i3)%b,prq=p234q,bef=,aft=
      cz1234(1,i3)%e(mu)=(lz1_234(mu)%a(1)*r2_34(i3)%a(1)+lz1_23
     & 4(mu)%c(1)*p234q*r2_34(i3)%b(2))
      cz1234(2,i3)%e(mu)=(lz1_234(mu)%c(2)*p234q*r2_34(i3)%b(1)+
     & lz1_234(mu)%a(2)*r2_34(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz1234(&,i3)%e(mu),a1=l1_34(i3)%a,c1=l1_34(i3)%c,a2=rz2_13
* 4(mu)%a,b2=rz2_134(mu)%b,prq=p134q,bef=cz1234(&,i3)%e(mu)+,aft=
      cz1234(1,i3)%e(mu)=cz1234(1,i3)%e(mu)+(l1_34(i3)%a(1)*rz2_
     & 134(mu)%a(1)+l1_34(i3)%c(1)*p134q*rz2_134(mu)%b(2))
      cz1234(2,i3)%e(mu)=cz1234(2,i3)%e(mu)+(l1_34(i3)%c(2)*p134
     & q*rz2_134(mu)%b(1)+l1_34(i3)%a(2)*rz2_134(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz1234(i1,&)%e(mu),a1=lz3_412(mu)%a,c1=lz3_412(mu)%c,a2=r4
* _12(i1)%a,b2=r4_12(i1)%b,prq=p412q,bef=cz1234(i1,&)%e(mu)+,aft=
      cz1234(i1,1)%e(mu)=cz1234(i1,1)%e(mu)+(lz3_412(mu)%a(1)*r4
     & _12(i1)%a(1)+lz3_412(mu)%c(1)*p412q*r4_12(i1)%b(2))
      cz1234(i1,2)%e(mu)=cz1234(i1,2)%e(mu)+(lz3_412(mu)%c(2)*p4
     & 12q*r4_12(i1)%b(1)+lz3_412(mu)%a(2)*r4_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz1234(i1,&)%e(mu),a1=l3_12(i1)%a,c1=l3_12(i1)%c,a2=rz4_31
* 2(mu)%a,b2=rz4_312(mu)%b,prq=p312q,bef=cz1234(i1,&)%e(mu)+,aft=
      cz1234(i1,1)%e(mu)=cz1234(i1,1)%e(mu)+(l3_12(i1)%a(1)*rz4_
     & 312(mu)%a(1)+l3_12(i1)%c(1)*p312q*rz4_312(mu)%b(2))
      cz1234(i1,2)%e(mu)=cz1234(i1,2)%e(mu)+(l3_12(i1)%c(2)*p312
     & q*rz4_312(mu)%b(1)+l3_12(i1)%a(2)*rz4_312(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cz1234(i1,i3)%e
      cz1234(i1,i3)%ek0=cz1234(i1,i3)%e(0)-cz1234(i1,i3)%e(1)
      end do
      end do
  
      if (ineutri(id1).ne.1) then
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf1234(&,i3)%e(mu),a1=lf1_234(mu)%a,c1=lf1_234(mu)%c,a2=r2
* _34(i3)%a,b2=r2_34(i3)%b,prq=p234q,bef=,aft=
      cf1234(1,i3)%e(mu)=(lf1_234(mu)%a(1)*r2_34(i3)%a(1)+lf1_23
     & 4(mu)%c(1)*p234q*r2_34(i3)%b(2))
      cf1234(2,i3)%e(mu)=(lf1_234(mu)%c(2)*p234q*r2_34(i3)%b(1)+
     & lf1_234(mu)%a(2)*r2_34(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf1234(&,i3)%e(mu),a1=l1_34(i3)%a,c1=l1_34(i3)%c,a2=rf2_13
* 4(mu)%a,b2=rf2_134(mu)%b,prq=p134q,bef=cf1234(&,i3)%e(mu)+,aft=
      cf1234(1,i3)%e(mu)=cf1234(1,i3)%e(mu)+(l1_34(i3)%a(1)*rf2_
     & 134(mu)%a(1)+l1_34(i3)%c(1)*p134q*rf2_134(mu)%b(2))
      cf1234(2,i3)%e(mu)=cf1234(2,i3)%e(mu)+(l1_34(i3)%c(2)*p134
     & q*rf2_134(mu)%b(1)+l1_34(i3)%a(2)*rf2_134(mu)%a(2))
      end do
      end do
  
      else
  
       do i1=1,2
       do i3=1,2
       do mu=0,3
        cf1234(i1,i3)%e(mu)=czero
       enddo
       enddo
       enddo
  
      endif
  
      if (ineutri(id3).ne.1) then
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf1234(i1,&)%e(mu),a1=lf3_412(mu)%a,c1=lf3_412(mu)%c,a2=r4
* _12(i1)%a,b2=r4_12(i1)%b,prq=p412q,bef=cf1234(i1,&)%e(mu)+,aft=
      cf1234(i1,1)%e(mu)=cf1234(i1,1)%e(mu)+(lf3_412(mu)%a(1)*r4
     & _12(i1)%a(1)+lf3_412(mu)%c(1)*p412q*r4_12(i1)%b(2))
      cf1234(i1,2)%e(mu)=cf1234(i1,2)%e(mu)+(lf3_412(mu)%c(2)*p4
     & 12q*r4_12(i1)%b(1)+lf3_412(mu)%a(2)*r4_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf1234(i1,&)%e(mu),a1=l3_12(i1)%a,c1=l3_12(i1)%c,a2=rf4_31
* 2(mu)%a,b2=rf4_312(mu)%b,prq=p312q,bef=cf1234(i1,&)%e(mu)+,aft=
      cf1234(i1,1)%e(mu)=cf1234(i1,1)%e(mu)+(l3_12(i1)%a(1)*rf4_
     & 312(mu)%a(1)+l3_12(i1)%c(1)*p312q*rf4_312(mu)%b(2))
      cf1234(i1,2)%e(mu)=cf1234(i1,2)%e(mu)+(l3_12(i1)%c(2)*p312
     & q*rf4_312(mu)%b(1)+l3_12(i1)%a(2)*rf4_312(mu)%a(2))
      end do
      end do
      endif
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cf1234(i1,i3)%e
      cf1234(i1,i3)%ek0=cf1234(i1,i3)%e(0)-cf1234(i1,i3)%e(1)
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz1256(&,i3)%e(mu),a1=lz1_256(mu)%a,c1=lz1_256(mu)%c,a2=r2
* _56(i3)%a,b2=r2_56(i3)%b,prq=p256q,bef=,aft=
      cz1256(1,i3)%e(mu)=(lz1_256(mu)%a(1)*r2_56(i3)%a(1)+lz1_25
     & 6(mu)%c(1)*p256q*r2_56(i3)%b(2))
      cz1256(2,i3)%e(mu)=(lz1_256(mu)%c(2)*p256q*r2_56(i3)%b(1)+
     & lz1_256(mu)%a(2)*r2_56(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz1256(&,i3)%e(mu),a1=l1_56(i3)%a,c1=l1_56(i3)%c,a2=rz2_15
* 6(mu)%a,b2=rz2_156(mu)%b,prq=p156q,bef=cz1256(&,i3)%e(mu)+,aft=
      cz1256(1,i3)%e(mu)=cz1256(1,i3)%e(mu)+(l1_56(i3)%a(1)*rz2_
     & 156(mu)%a(1)+l1_56(i3)%c(1)*p156q*rz2_156(mu)%b(2))
      cz1256(2,i3)%e(mu)=cz1256(2,i3)%e(mu)+(l1_56(i3)%c(2)*p156
     & q*rz2_156(mu)%b(1)+l1_56(i3)%a(2)*rz2_156(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz1256(i1,&)%e(mu),a1=lz5_612(mu)%a,c1=lz5_612(mu)%c,a2=r6
* _12(i1)%a,b2=r6_12(i1)%b,prq=p612q,bef=cz1256(i1,&)%e(mu)+,aft=
      cz1256(i1,1)%e(mu)=cz1256(i1,1)%e(mu)+(lz5_612(mu)%a(1)*r6
     & _12(i1)%a(1)+lz5_612(mu)%c(1)*p612q*r6_12(i1)%b(2))
      cz1256(i1,2)%e(mu)=cz1256(i1,2)%e(mu)+(lz5_612(mu)%c(2)*p6
     & 12q*r6_12(i1)%b(1)+lz5_612(mu)%a(2)*r6_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz1256(i1,&)%e(mu),a1=l5_12(i1)%a,c1=l5_12(i1)%c,a2=rz6_51
* 2(mu)%a,b2=rz6_512(mu)%b,prq=p512q,bef=cz1256(i1,&)%e(mu)+,aft=
      cz1256(i1,1)%e(mu)=cz1256(i1,1)%e(mu)+(l5_12(i1)%a(1)*rz6_
     & 512(mu)%a(1)+l5_12(i1)%c(1)*p512q*rz6_512(mu)%b(2))
      cz1256(i1,2)%e(mu)=cz1256(i1,2)%e(mu)+(l5_12(i1)%c(2)*p512
     & q*rz6_512(mu)%b(1)+l5_12(i1)%a(2)*rz6_512(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cz1256(i1,i3)%e
      cz1256(i1,i3)%ek0=cz1256(i1,i3)%e(0)-cz1256(i1,i3)%e(1)
      end do
      end do
  
      if (ineutri(id1).ne.1) then
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf1256(&,i3)%e(mu),a1=lf1_256(mu)%a,c1=lf1_256(mu)%c,a2=r2
* _56(i3)%a,b2=r2_56(i3)%b,prq=p256q,bef=,aft=
      cf1256(1,i3)%e(mu)=(lf1_256(mu)%a(1)*r2_56(i3)%a(1)+lf1_25
     & 6(mu)%c(1)*p256q*r2_56(i3)%b(2))
      cf1256(2,i3)%e(mu)=(lf1_256(mu)%c(2)*p256q*r2_56(i3)%b(1)+
     & lf1_256(mu)%a(2)*r2_56(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf1256(&,i3)%e(mu),a1=l1_56(i3)%a,c1=l1_56(i3)%c,a2=rf2_15
* 6(mu)%a,b2=rf2_156(mu)%b,prq=p156q,bef=cf1256(&,i3)%e(mu)+,aft=
      cf1256(1,i3)%e(mu)=cf1256(1,i3)%e(mu)+(l1_56(i3)%a(1)*rf2_
     & 156(mu)%a(1)+l1_56(i3)%c(1)*p156q*rf2_156(mu)%b(2))
      cf1256(2,i3)%e(mu)=cf1256(2,i3)%e(mu)+(l1_56(i3)%c(2)*p156
     & q*rf2_156(mu)%b(1)+l1_56(i3)%a(2)*rf2_156(mu)%a(2))
      end do
      end do
  
      else
  
       do i1=1,2
       do i3=1,2
       do mu=0,3
        cf1256(i1,i3)%e(mu)=czero
       enddo
       enddo
       enddo
  
      endif
  
      if (ineutri(id5).ne.1) then
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf1256(i1,&)%e(mu),a1=lf5_612(mu)%a,c1=lf5_612(mu)%c,a2=r6
* _12(i1)%a,b2=r6_12(i1)%b,prq=p612q,bef=cf1256(i1,&)%e(mu)+,aft=
      cf1256(i1,1)%e(mu)=cf1256(i1,1)%e(mu)+(lf5_612(mu)%a(1)*r6
     & _12(i1)%a(1)+lf5_612(mu)%c(1)*p612q*r6_12(i1)%b(2))
      cf1256(i1,2)%e(mu)=cf1256(i1,2)%e(mu)+(lf5_612(mu)%c(2)*p6
     & 12q*r6_12(i1)%b(1)+lf5_612(mu)%a(2)*r6_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf1256(i1,&)%e(mu),a1=l5_12(i1)%a,c1=l5_12(i1)%c,a2=rf6_51
* 2(mu)%a,b2=rf6_512(mu)%b,prq=p512q,bef=cf1256(i1,&)%e(mu)+,aft=
      cf1256(i1,1)%e(mu)=cf1256(i1,1)%e(mu)+(l5_12(i1)%a(1)*rf6_
     & 512(mu)%a(1)+l5_12(i1)%c(1)*p512q*rf6_512(mu)%b(2))
      cf1256(i1,2)%e(mu)=cf1256(i1,2)%e(mu)+(l5_12(i1)%c(2)*p512
     & q*rf6_512(mu)%b(1)+l5_12(i1)%a(2)*rf6_512(mu)%a(2))
      end do
      end do
      endif
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cf1256(i1,i3)%e
      cf1256(i1,i3)%ek0=cf1256(i1,i3)%e(0)-cf1256(i1,i3)%e(1)
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz1278(&,i3)%e(mu),a1=lz1_278(mu)%a,c1=lz1_278(mu)%c,a2=r2
* _78(i3)%a,b2=r2_78(i3)%b,prq=p278q,bef=,aft=
      cz1278(1,i3)%e(mu)=(lz1_278(mu)%a(1)*r2_78(i3)%a(1)+lz1_27
     & 8(mu)%c(1)*p278q*r2_78(i3)%b(2))
      cz1278(2,i3)%e(mu)=(lz1_278(mu)%c(2)*p278q*r2_78(i3)%b(1)+
     & lz1_278(mu)%a(2)*r2_78(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz1278(&,i3)%e(mu),a1=l1_78(i3)%a,c1=l1_78(i3)%c,a2=rz2_17
* 8(mu)%a,b2=rz2_178(mu)%b,prq=p178q,bef=cz1278(&,i3)%e(mu)+,aft=
      cz1278(1,i3)%e(mu)=cz1278(1,i3)%e(mu)+(l1_78(i3)%a(1)*rz2_
     & 178(mu)%a(1)+l1_78(i3)%c(1)*p178q*rz2_178(mu)%b(2))
      cz1278(2,i3)%e(mu)=cz1278(2,i3)%e(mu)+(l1_78(i3)%c(2)*p178
     & q*rz2_178(mu)%b(1)+l1_78(i3)%a(2)*rz2_178(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz1278(i1,&)%e(mu),a1=lz7_812(mu)%a,c1=lz7_812(mu)%c,a2=r8
* _12(i1)%a,b2=r8_12(i1)%b,prq=p812q,bef=cz1278(i1,&)%e(mu)+,aft=
      cz1278(i1,1)%e(mu)=cz1278(i1,1)%e(mu)+(lz7_812(mu)%a(1)*r8
     & _12(i1)%a(1)+lz7_812(mu)%c(1)*p812q*r8_12(i1)%b(2))
      cz1278(i1,2)%e(mu)=cz1278(i1,2)%e(mu)+(lz7_812(mu)%c(2)*p8
     & 12q*r8_12(i1)%b(1)+lz7_812(mu)%a(2)*r8_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz1278(i1,&)%e(mu),a1=l7_12(i1)%a,c1=l7_12(i1)%c,a2=rz8_71
* 2(mu)%a,b2=rz8_712(mu)%b,prq=p712q,bef=cz1278(i1,&)%e(mu)+,aft=
      cz1278(i1,1)%e(mu)=cz1278(i1,1)%e(mu)+(l7_12(i1)%a(1)*rz8_
     & 712(mu)%a(1)+l7_12(i1)%c(1)*p712q*rz8_712(mu)%b(2))
      cz1278(i1,2)%e(mu)=cz1278(i1,2)%e(mu)+(l7_12(i1)%c(2)*p712
     & q*rz8_712(mu)%b(1)+l7_12(i1)%a(2)*rz8_712(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cz1278(i1,i3)%e
      cz1278(i1,i3)%ek0=cz1278(i1,i3)%e(0)-cz1278(i1,i3)%e(1)
      end do
      end do
  
      if (ineutri(id1).ne.1) then
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf1278(&,i3)%e(mu),a1=lf1_278(mu)%a,c1=lf1_278(mu)%c,a2=r2
* _78(i3)%a,b2=r2_78(i3)%b,prq=p278q,bef=,aft=
      cf1278(1,i3)%e(mu)=(lf1_278(mu)%a(1)*r2_78(i3)%a(1)+lf1_27
     & 8(mu)%c(1)*p278q*r2_78(i3)%b(2))
      cf1278(2,i3)%e(mu)=(lf1_278(mu)%c(2)*p278q*r2_78(i3)%b(1)+
     & lf1_278(mu)%a(2)*r2_78(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf1278(&,i3)%e(mu),a1=l1_78(i3)%a,c1=l1_78(i3)%c,a2=rf2_17
* 8(mu)%a,b2=rf2_178(mu)%b,prq=p178q,bef=cf1278(&,i3)%e(mu)+,aft=
      cf1278(1,i3)%e(mu)=cf1278(1,i3)%e(mu)+(l1_78(i3)%a(1)*rf2_
     & 178(mu)%a(1)+l1_78(i3)%c(1)*p178q*rf2_178(mu)%b(2))
      cf1278(2,i3)%e(mu)=cf1278(2,i3)%e(mu)+(l1_78(i3)%c(2)*p178
     & q*rf2_178(mu)%b(1)+l1_78(i3)%a(2)*rf2_178(mu)%a(2))
      end do
      end do
  
      else
  
       do i1=1,2
       do i3=1,2
       do mu=0,3
        cf1278(i1,i3)%e(mu)=czero
       enddo
       enddo
       enddo
  
      endif
  
      if (ineutri(id7).ne.1) then
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf1278(i1,&)%e(mu),a1=lf7_812(mu)%a,c1=lf7_812(mu)%c,a2=r8
* _12(i1)%a,b2=r8_12(i1)%b,prq=p812q,bef=cf1278(i1,&)%e(mu)+,aft=
      cf1278(i1,1)%e(mu)=cf1278(i1,1)%e(mu)+(lf7_812(mu)%a(1)*r8
     & _12(i1)%a(1)+lf7_812(mu)%c(1)*p812q*r8_12(i1)%b(2))
      cf1278(i1,2)%e(mu)=cf1278(i1,2)%e(mu)+(lf7_812(mu)%c(2)*p8
     & 12q*r8_12(i1)%b(1)+lf7_812(mu)%a(2)*r8_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf1278(i1,&)%e(mu),a1=l7_12(i1)%a,c1=l7_12(i1)%c,a2=rf8_71
* 2(mu)%a,b2=rf8_712(mu)%b,prq=p712q,bef=cf1278(i1,&)%e(mu)+,aft=
      cf1278(i1,1)%e(mu)=cf1278(i1,1)%e(mu)+(l7_12(i1)%a(1)*rf8_
     & 712(mu)%a(1)+l7_12(i1)%c(1)*p712q*rf8_712(mu)%b(2))
      cf1278(i1,2)%e(mu)=cf1278(i1,2)%e(mu)+(l7_12(i1)%c(2)*p712
     & q*rf8_712(mu)%b(1)+l7_12(i1)%a(2)*rf8_712(mu)%a(2))
      end do
      end do
      endif
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cf1278(i1,i3)%e
      cf1278(i1,i3)%ek0=cf1278(i1,i3)%e(0)-cf1278(i1,i3)%e(1)
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz3456(&,i3)%e(mu),a1=lz3_456(mu)%a,c1=lz3_456(mu)%c,a2=r4
* _56(i3)%a,b2=r4_56(i3)%b,prq=p456q,bef=,aft=
      cz3456(1,i3)%e(mu)=(lz3_456(mu)%a(1)*r4_56(i3)%a(1)+lz3_45
     & 6(mu)%c(1)*p456q*r4_56(i3)%b(2))
      cz3456(2,i3)%e(mu)=(lz3_456(mu)%c(2)*p456q*r4_56(i3)%b(1)+
     & lz3_456(mu)%a(2)*r4_56(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz3456(&,i3)%e(mu),a1=l3_56(i3)%a,c1=l3_56(i3)%c,a2=rz4_35
* 6(mu)%a,b2=rz4_356(mu)%b,prq=p356q,bef=cz3456(&,i3)%e(mu)+,aft=
      cz3456(1,i3)%e(mu)=cz3456(1,i3)%e(mu)+(l3_56(i3)%a(1)*rz4_
     & 356(mu)%a(1)+l3_56(i3)%c(1)*p356q*rz4_356(mu)%b(2))
      cz3456(2,i3)%e(mu)=cz3456(2,i3)%e(mu)+(l3_56(i3)%c(2)*p356
     & q*rz4_356(mu)%b(1)+l3_56(i3)%a(2)*rz4_356(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz3456(i1,&)%e(mu),a1=lz5_634(mu)%a,c1=lz5_634(mu)%c,a2=r6
* _34(i1)%a,b2=r6_34(i1)%b,prq=p634q,bef=cz3456(i1,&)%e(mu)+,aft=
      cz3456(i1,1)%e(mu)=cz3456(i1,1)%e(mu)+(lz5_634(mu)%a(1)*r6
     & _34(i1)%a(1)+lz5_634(mu)%c(1)*p634q*r6_34(i1)%b(2))
      cz3456(i1,2)%e(mu)=cz3456(i1,2)%e(mu)+(lz5_634(mu)%c(2)*p6
     & 34q*r6_34(i1)%b(1)+lz5_634(mu)%a(2)*r6_34(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz3456(i1,&)%e(mu),a1=l5_34(i1)%a,c1=l5_34(i1)%c,a2=rz6_53
* 4(mu)%a,b2=rz6_534(mu)%b,prq=p534q,bef=cz3456(i1,&)%e(mu)+,aft=
      cz3456(i1,1)%e(mu)=cz3456(i1,1)%e(mu)+(l5_34(i1)%a(1)*rz6_
     & 534(mu)%a(1)+l5_34(i1)%c(1)*p534q*rz6_534(mu)%b(2))
      cz3456(i1,2)%e(mu)=cz3456(i1,2)%e(mu)+(l5_34(i1)%c(2)*p534
     & q*rz6_534(mu)%b(1)+l5_34(i1)%a(2)*rz6_534(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cz3456(i1,i3)%e
      cz3456(i1,i3)%ek0=cz3456(i1,i3)%e(0)-cz3456(i1,i3)%e(1)
      end do
      end do
  
      if (ineutri(id3).ne.1) then
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf3456(&,i3)%e(mu),a1=lf3_456(mu)%a,c1=lf3_456(mu)%c,a2=r4
* _56(i3)%a,b2=r4_56(i3)%b,prq=p456q,bef=,aft=
      cf3456(1,i3)%e(mu)=(lf3_456(mu)%a(1)*r4_56(i3)%a(1)+lf3_45
     & 6(mu)%c(1)*p456q*r4_56(i3)%b(2))
      cf3456(2,i3)%e(mu)=(lf3_456(mu)%c(2)*p456q*r4_56(i3)%b(1)+
     & lf3_456(mu)%a(2)*r4_56(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf3456(&,i3)%e(mu),a1=l3_56(i3)%a,c1=l3_56(i3)%c,a2=rf4_35
* 6(mu)%a,b2=rf4_356(mu)%b,prq=p356q,bef=cf3456(&,i3)%e(mu)+,aft=
      cf3456(1,i3)%e(mu)=cf3456(1,i3)%e(mu)+(l3_56(i3)%a(1)*rf4_
     & 356(mu)%a(1)+l3_56(i3)%c(1)*p356q*rf4_356(mu)%b(2))
      cf3456(2,i3)%e(mu)=cf3456(2,i3)%e(mu)+(l3_56(i3)%c(2)*p356
     & q*rf4_356(mu)%b(1)+l3_56(i3)%a(2)*rf4_356(mu)%a(2))
      end do
      end do
  
      else
  
       do i1=1,2
       do i3=1,2
       do mu=0,3
        cf3456(i1,i3)%e(mu)=czero
       enddo
       enddo
       enddo
  
      endif
  
      if (ineutri(id5).ne.1) then
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf3456(i1,&)%e(mu),a1=lf5_634(mu)%a,c1=lf5_634(mu)%c,a2=r6
* _34(i1)%a,b2=r6_34(i1)%b,prq=p634q,bef=cf3456(i1,&)%e(mu)+,aft=
      cf3456(i1,1)%e(mu)=cf3456(i1,1)%e(mu)+(lf5_634(mu)%a(1)*r6
     & _34(i1)%a(1)+lf5_634(mu)%c(1)*p634q*r6_34(i1)%b(2))
      cf3456(i1,2)%e(mu)=cf3456(i1,2)%e(mu)+(lf5_634(mu)%c(2)*p6
     & 34q*r6_34(i1)%b(1)+lf5_634(mu)%a(2)*r6_34(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf3456(i1,&)%e(mu),a1=l5_34(i1)%a,c1=l5_34(i1)%c,a2=rf6_53
* 4(mu)%a,b2=rf6_534(mu)%b,prq=p534q,bef=cf3456(i1,&)%e(mu)+,aft=
      cf3456(i1,1)%e(mu)=cf3456(i1,1)%e(mu)+(l5_34(i1)%a(1)*rf6_
     & 534(mu)%a(1)+l5_34(i1)%c(1)*p534q*rf6_534(mu)%b(2))
      cf3456(i1,2)%e(mu)=cf3456(i1,2)%e(mu)+(l5_34(i1)%c(2)*p534
     & q*rf6_534(mu)%b(1)+l5_34(i1)%a(2)*rf6_534(mu)%a(2))
      end do
      end do
      endif
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cf3456(i1,i3)%e
      cf3456(i1,i3)%ek0=cf3456(i1,i3)%e(0)-cf3456(i1,i3)%e(1)
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz3478(&,i3)%e(mu),a1=lz3_478(mu)%a,c1=lz3_478(mu)%c,a2=r4
* _78(i3)%a,b2=r4_78(i3)%b,prq=p478q,bef=,aft=
      cz3478(1,i3)%e(mu)=(lz3_478(mu)%a(1)*r4_78(i3)%a(1)+lz3_47
     & 8(mu)%c(1)*p478q*r4_78(i3)%b(2))
      cz3478(2,i3)%e(mu)=(lz3_478(mu)%c(2)*p478q*r4_78(i3)%b(1)+
     & lz3_478(mu)%a(2)*r4_78(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz3478(&,i3)%e(mu),a1=l3_78(i3)%a,c1=l3_78(i3)%c,a2=rz4_37
* 8(mu)%a,b2=rz4_378(mu)%b,prq=p378q,bef=cz3478(&,i3)%e(mu)+,aft=
      cz3478(1,i3)%e(mu)=cz3478(1,i3)%e(mu)+(l3_78(i3)%a(1)*rz4_
     & 378(mu)%a(1)+l3_78(i3)%c(1)*p378q*rz4_378(mu)%b(2))
      cz3478(2,i3)%e(mu)=cz3478(2,i3)%e(mu)+(l3_78(i3)%c(2)*p378
     & q*rz4_378(mu)%b(1)+l3_78(i3)%a(2)*rz4_378(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz3478(i1,&)%e(mu),a1=lz7_834(mu)%a,c1=lz7_834(mu)%c,a2=r8
* _34(i1)%a,b2=r8_34(i1)%b,prq=p834q,bef=cz3478(i1,&)%e(mu)+,aft=
      cz3478(i1,1)%e(mu)=cz3478(i1,1)%e(mu)+(lz7_834(mu)%a(1)*r8
     & _34(i1)%a(1)+lz7_834(mu)%c(1)*p834q*r8_34(i1)%b(2))
      cz3478(i1,2)%e(mu)=cz3478(i1,2)%e(mu)+(lz7_834(mu)%c(2)*p8
     & 34q*r8_34(i1)%b(1)+lz7_834(mu)%a(2)*r8_34(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz3478(i1,&)%e(mu),a1=l7_34(i1)%a,c1=l7_34(i1)%c,a2=rz8_73
* 4(mu)%a,b2=rz8_734(mu)%b,prq=p734q,bef=cz3478(i1,&)%e(mu)+,aft=
      cz3478(i1,1)%e(mu)=cz3478(i1,1)%e(mu)+(l7_34(i1)%a(1)*rz8_
     & 734(mu)%a(1)+l7_34(i1)%c(1)*p734q*rz8_734(mu)%b(2))
      cz3478(i1,2)%e(mu)=cz3478(i1,2)%e(mu)+(l7_34(i1)%c(2)*p734
     & q*rz8_734(mu)%b(1)+l7_34(i1)%a(2)*rz8_734(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cz3478(i1,i3)%e
      cz3478(i1,i3)%ek0=cz3478(i1,i3)%e(0)-cz3478(i1,i3)%e(1)
      end do
      end do
  
      if (ineutri(id3).ne.1) then
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf3478(&,i3)%e(mu),a1=lf3_478(mu)%a,c1=lf3_478(mu)%c,a2=r4
* _78(i3)%a,b2=r4_78(i3)%b,prq=p478q,bef=,aft=
      cf3478(1,i3)%e(mu)=(lf3_478(mu)%a(1)*r4_78(i3)%a(1)+lf3_47
     & 8(mu)%c(1)*p478q*r4_78(i3)%b(2))
      cf3478(2,i3)%e(mu)=(lf3_478(mu)%c(2)*p478q*r4_78(i3)%b(1)+
     & lf3_478(mu)%a(2)*r4_78(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf3478(&,i3)%e(mu),a1=l3_78(i3)%a,c1=l3_78(i3)%c,a2=rf4_37
* 8(mu)%a,b2=rf4_378(mu)%b,prq=p378q,bef=cf3478(&,i3)%e(mu)+,aft=
      cf3478(1,i3)%e(mu)=cf3478(1,i3)%e(mu)+(l3_78(i3)%a(1)*rf4_
     & 378(mu)%a(1)+l3_78(i3)%c(1)*p378q*rf4_378(mu)%b(2))
      cf3478(2,i3)%e(mu)=cf3478(2,i3)%e(mu)+(l3_78(i3)%c(2)*p378
     & q*rf4_378(mu)%b(1)+l3_78(i3)%a(2)*rf4_378(mu)%a(2))
      end do
      end do
  
      else
  
       do i1=1,2
       do i3=1,2
       do mu=0,3
        cf3478(i1,i3)%e(mu)=czero
       enddo
       enddo
       enddo
  
      endif
  
      if (ineutri(id7).ne.1) then
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf3478(i1,&)%e(mu),a1=lf7_834(mu)%a,c1=lf7_834(mu)%c,a2=r8
* _34(i1)%a,b2=r8_34(i1)%b,prq=p834q,bef=cf3478(i1,&)%e(mu)+,aft=
      cf3478(i1,1)%e(mu)=cf3478(i1,1)%e(mu)+(lf7_834(mu)%a(1)*r8
     & _34(i1)%a(1)+lf7_834(mu)%c(1)*p834q*r8_34(i1)%b(2))
      cf3478(i1,2)%e(mu)=cf3478(i1,2)%e(mu)+(lf7_834(mu)%c(2)*p8
     & 34q*r8_34(i1)%b(1)+lf7_834(mu)%a(2)*r8_34(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf3478(i1,&)%e(mu),a1=l7_34(i1)%a,c1=l7_34(i1)%c,a2=rf8_73
* 4(mu)%a,b2=rf8_734(mu)%b,prq=p734q,bef=cf3478(i1,&)%e(mu)+,aft=
      cf3478(i1,1)%e(mu)=cf3478(i1,1)%e(mu)+(l7_34(i1)%a(1)*rf8_
     & 734(mu)%a(1)+l7_34(i1)%c(1)*p734q*rf8_734(mu)%b(2))
      cf3478(i1,2)%e(mu)=cf3478(i1,2)%e(mu)+(l7_34(i1)%c(2)*p734
     & q*rf8_734(mu)%b(1)+l7_34(i1)%a(2)*rf8_734(mu)%a(2))
      end do
      end do
      endif
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cf3478(i1,i3)%e
      cf3478(i1,i3)%ek0=cf3478(i1,i3)%e(0)-cf3478(i1,i3)%e(1)
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz5678(&,i3)%e(mu),a1=lz5_678(mu)%a,c1=lz5_678(mu)%c,a2=r6
* _78(i3)%a,b2=r6_78(i3)%b,prq=p678q,bef=,aft=
      cz5678(1,i3)%e(mu)=(lz5_678(mu)%a(1)*r6_78(i3)%a(1)+lz5_67
     & 8(mu)%c(1)*p678q*r6_78(i3)%b(2))
      cz5678(2,i3)%e(mu)=(lz5_678(mu)%c(2)*p678q*r6_78(i3)%b(1)+
     & lz5_678(mu)%a(2)*r6_78(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cz5678(&,i3)%e(mu),a1=l5_78(i3)%a,c1=l5_78(i3)%c,a2=rz6_57
* 8(mu)%a,b2=rz6_578(mu)%b,prq=p578q,bef=cz5678(&,i3)%e(mu)+,aft=
      cz5678(1,i3)%e(mu)=cz5678(1,i3)%e(mu)+(l5_78(i3)%a(1)*rz6_
     & 578(mu)%a(1)+l5_78(i3)%c(1)*p578q*rz6_578(mu)%b(2))
      cz5678(2,i3)%e(mu)=cz5678(2,i3)%e(mu)+(l5_78(i3)%c(2)*p578
     & q*rz6_578(mu)%b(1)+l5_78(i3)%a(2)*rz6_578(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz5678(i1,&)%e(mu),a1=lz7_856(mu)%a,c1=lz7_856(mu)%c,a2=r8
* _56(i1)%a,b2=r8_56(i1)%b,prq=p856q,bef=cz5678(i1,&)%e(mu)+,aft=
      cz5678(i1,1)%e(mu)=cz5678(i1,1)%e(mu)+(lz7_856(mu)%a(1)*r8
     & _56(i1)%a(1)+lz7_856(mu)%c(1)*p856q*r8_56(i1)%b(2))
      cz5678(i1,2)%e(mu)=cz5678(i1,2)%e(mu)+(lz7_856(mu)%c(2)*p8
     & 56q*r8_56(i1)%b(1)+lz7_856(mu)%a(2)*r8_56(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cz5678(i1,&)%e(mu),a1=l7_56(i1)%a,c1=l7_56(i1)%c,a2=rz8_75
* 6(mu)%a,b2=rz8_756(mu)%b,prq=p756q,bef=cz5678(i1,&)%e(mu)+,aft=
      cz5678(i1,1)%e(mu)=cz5678(i1,1)%e(mu)+(l7_56(i1)%a(1)*rz8_
     & 756(mu)%a(1)+l7_56(i1)%c(1)*p756q*rz8_756(mu)%b(2))
      cz5678(i1,2)%e(mu)=cz5678(i1,2)%e(mu)+(l7_56(i1)%c(2)*p756
     & q*rz8_756(mu)%b(1)+l7_56(i1)%a(2)*rz8_756(mu)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cz5678(i1,i3)%e
      cz5678(i1,i3)%ek0=cz5678(i1,i3)%e(0)-cz5678(i1,i3)%e(1)
      end do
      end do
  
      if (ineutri(id5).ne.1) then
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf5678(&,i3)%e(mu),a1=lf5_678(mu)%a,c1=lf5_678(mu)%c,a2=r6
* _78(i3)%a,b2=r6_78(i3)%b,prq=p678q,bef=,aft=
      cf5678(1,i3)%e(mu)=(lf5_678(mu)%a(1)*r6_78(i3)%a(1)+lf5_67
     & 8(mu)%c(1)*p678q*r6_78(i3)%b(2))
      cf5678(2,i3)%e(mu)=(lf5_678(mu)%c(2)*p678q*r6_78(i3)%b(1)+
     & lf5_678(mu)%a(2)*r6_78(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf5678(&,i3)%e(mu),a1=l5_78(i3)%a,c1=l5_78(i3)%c,a2=rf6_57
* 8(mu)%a,b2=rf6_578(mu)%b,prq=p578q,bef=cf5678(&,i3)%e(mu)+,aft=
      cf5678(1,i3)%e(mu)=cf5678(1,i3)%e(mu)+(l5_78(i3)%a(1)*rf6_
     & 578(mu)%a(1)+l5_78(i3)%c(1)*p578q*rf6_578(mu)%b(2))
      cf5678(2,i3)%e(mu)=cf5678(2,i3)%e(mu)+(l5_78(i3)%c(2)*p578
     & q*rf6_578(mu)%b(1)+l5_78(i3)%a(2)*rf6_578(mu)%a(2))
      end do
      end do
  
      else
  
       do i1=1,2
       do i3=1,2
       do mu=0,3
        cf5678(i1,i3)%e(mu)=czero
       enddo
       enddo
       enddo
  
      endif
  
      if (ineutri(id7).ne.1) then
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf5678(i1,&)%e(mu),a1=lf7_856(mu)%a,c1=lf7_856(mu)%c,a2=r8
* _56(i1)%a,b2=r8_56(i1)%b,prq=p856q,bef=cf5678(i1,&)%e(mu)+,aft=
      cf5678(i1,1)%e(mu)=cf5678(i1,1)%e(mu)+(lf7_856(mu)%a(1)*r8
     & _56(i1)%a(1)+lf7_856(mu)%c(1)*p856q*r8_56(i1)%b(2))
      cf5678(i1,2)%e(mu)=cf5678(i1,2)%e(mu)+(lf7_856(mu)%c(2)*p8
     & 56q*r8_56(i1)%b(1)+lf7_856(mu)%a(2)*r8_56(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cf5678(i1,&)%e(mu),a1=l7_56(i1)%a,c1=l7_56(i1)%c,a2=rf8_75
* 6(mu)%a,b2=rf8_756(mu)%b,prq=p756q,bef=cf5678(i1,&)%e(mu)+,aft=
      cf5678(i1,1)%e(mu)=cf5678(i1,1)%e(mu)+(l7_56(i1)%a(1)*rf8_
     & 756(mu)%a(1)+l7_56(i1)%c(1)*p756q*rf8_756(mu)%b(2))
      cf5678(i1,2)%e(mu)=cf5678(i1,2)%e(mu)+(l7_56(i1)%c(2)*p756
     & q*rf8_756(mu)%b(1)+l7_56(i1)%a(2)*rf8_756(mu)%a(2))
      end do
      end do
      endif
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cf5678(i1,i3)%e
      cf5678(i1,i3)%ek0=cf5678(i1,i3)%e(0)-cf5678(i1,i3)%e(1)
      end do
      end do
  
  
  
*                                                                       
* BOSON CONNECTION                                                      
*                                                                       
  
*1234                                                                   
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
* p.q -- p.q=caux4(i1,i3,i5,i7),p=cz1234(i1,i3)%e,q=cz5678(i5,i7)%e,bef=
* ,aft=
      caux4(i1,i3,i5,i7)=(cz1234(i1,i3)%e(0)*cz5678(i5,i7)%e(0)-
     & cz1234(i1,i3)%e(1)*cz5678(i5,i7)%e(1)-cz1234(i1,i3)%e(2)*
     & cz5678(i5,i7)%e(2)-cz1234(i1,i3)%e(3)*cz5678(i5,i7)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
        cbct1234(i1,i3,i5,i7) = caux4(i1,i3,i5,i7)/(-p1234q+cmz2)
      enddo
      enddo
      enddo
      enddo
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
* p.q -- p.q=caux4(i1,i3,i5,i7),p=cf1234(i1,i3)%e,q=cf5678(i5,i7)%e,bef=
* ,aft=
      caux4(i1,i3,i5,i7)=(cf1234(i1,i3)%e(0)*cf5678(i5,i7)%e(0)-
     & cf1234(i1,i3)%e(1)*cf5678(i5,i7)%e(1)-cf1234(i1,i3)%e(2)*
     & cf5678(i5,i7)%e(2)-cf1234(i1,i3)%e(3)*cf5678(i5,i7)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
        cbct1234(i1,i3,i5,i7) = cbct1234(i1,i3,i5,i7)+
     &    caux4(i1,i3,i5,i7)/(-p1234q)
      enddo
      enddo
      enddo
      enddo
  
*1256                                                                   
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
* p.q -- p.q=caux4(i1,i3,i5,i7),p=cz1256(i1,i5)%e,q=cz3478(i3,i7)%e,bef=
* ,aft=
      caux4(i1,i3,i5,i7)=(cz1256(i1,i5)%e(0)*cz3478(i3,i7)%e(0)-
     & cz1256(i1,i5)%e(1)*cz3478(i3,i7)%e(1)-cz1256(i1,i5)%e(2)*
     & cz3478(i3,i7)%e(2)-cz1256(i1,i5)%e(3)*cz3478(i3,i7)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
        cbct1256(i1,i3,i5,i7) = caux4(i1,i3,i5,i7)/(-p1256q+cmz2)
      enddo
      enddo
      enddo
      enddo
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
* p.q -- p.q=caux4(i1,i3,i5,i7),p=cf1256(i1,i5)%e,q=cf3478(i3,i7)%e,bef=
* ,aft=
      caux4(i1,i3,i5,i7)=(cf1256(i1,i5)%e(0)*cf3478(i3,i7)%e(0)-
     & cf1256(i1,i5)%e(1)*cf3478(i3,i7)%e(1)-cf1256(i1,i5)%e(2)*
     & cf3478(i3,i7)%e(2)-cf1256(i1,i5)%e(3)*cf3478(i3,i7)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
        cbct1256(i1,i3,i5,i7) = cbct1256(i1,i3,i5,i7)+
     &    caux4(i1,i3,i5,i7)/(-p1256q)
      enddo
      enddo
      enddo
      enddo
  
*1278                                                                   
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
* p.q -- p.q=caux4(i1,i3,i5,i7),p=cz1278(i1,i7)%e,q=cz3456(i3,i5)%e,bef=
* ,aft=
      caux4(i1,i3,i5,i7)=(cz1278(i1,i7)%e(0)*cz3456(i3,i5)%e(0)-
     & cz1278(i1,i7)%e(1)*cz3456(i3,i5)%e(1)-cz1278(i1,i7)%e(2)*
     & cz3456(i3,i5)%e(2)-cz1278(i1,i7)%e(3)*cz3456(i3,i5)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
        cbct1278(i1,i3,i5,i7) = caux4(i1,i3,i5,i7)/(-p1278q+cmz2)
      enddo
      enddo
      enddo
      enddo
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
* p.q -- p.q=caux4(i1,i3,i5,i7),p=cf1278(i1,i7)%e,q=cf3456(i3,i5)%e,bef=
* ,aft=
      caux4(i1,i3,i5,i7)=(cf1278(i1,i7)%e(0)*cf3456(i3,i5)%e(0)-
     & cf1278(i1,i7)%e(1)*cf3456(i3,i5)%e(1)-cf1278(i1,i7)%e(2)*
     & cf3456(i3,i5)%e(2)-cf1278(i1,i7)%e(3)*cf3456(i3,i5)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
        cbct1278(i1,i3,i5,i7) = cbct1278(i1,i3,i5,i7)+
     &    caux4(i1,i3,i5,i7)/(-p1278q)
      enddo
      enddo
      enddo
      enddo
  
  
*                                                                       
* BOSON FUSION                                                          
*                                                                       
       if (rmh.ge.0.d0) then
  
      do i1=1,2
      do i7=1,2
* p.q -- p.q=ctrip7812(i1,i7),p=cfkz78(i7)%e,q=cfkz12(i1)%e,bef=,aft=
      ctrip7812(i1,i7)=(cfkz78(i7)%e(0)*cfkz12(i1)%e(0)-cfkz78(i
     & 7)%e(1)*cfkz12(i1)%e(1)-cfkz78(i7)%e(2)*cfkz12(i1)%e(2)-c
     & fkz78(i7)%e(3)*cfkz12(i1)%e(3))
      end do
      end do
      do i3=1,2
      do i5=1,2
* p.q -- p.q=ctrip3456(i3,i5),p=cfkz34(i3)%e,q=cfkz56(i5)%e,bef=,aft=
      ctrip3456(i3,i5)=(cfkz34(i3)%e(0)*cfkz56(i5)%e(0)-cfkz34(i
     & 3)%e(1)*cfkz56(i5)%e(1)-cfkz34(i3)%e(2)*cfkz56(i5)%e(2)-c
     & fkz34(i3)%e(3)*cfkz56(i5)%e(3))
      end do
      end do
  
      do i3=1,2
      do i7=1,2
* p.q -- p.q=ctrip7834(i3,i7),p=cfkz78(i7)%e,q=cfkz34(i3)%e,bef=,aft=
      ctrip7834(i3,i7)=(cfkz78(i7)%e(0)*cfkz34(i3)%e(0)-cfkz78(i
     & 7)%e(1)*cfkz34(i3)%e(1)-cfkz78(i7)%e(2)*cfkz34(i3)%e(2)-c
     & fkz78(i7)%e(3)*cfkz34(i3)%e(3))
      end do
      end do
      do i1=1,2
      do i5=1,2
* p.q -- p.q=ctrip1256(i1,i5),p=cfkz12(i1)%e,q=cfkz56(i5)%e,bef=,aft=
      ctrip1256(i1,i5)=(cfkz12(i1)%e(0)*cfkz56(i5)%e(0)-cfkz12(i
     & 1)%e(1)*cfkz56(i5)%e(1)-cfkz12(i1)%e(2)*cfkz56(i5)%e(2)-c
     & fkz12(i1)%e(3)*cfkz56(i5)%e(3))
      end do
      end do
  
      do i5=1,2
      do i7=1,2
* p.q -- p.q=ctrip7856(i5,i7),p=cfkz78(i7)%e,q=cfkz56(i5)%e,bef=,aft=
      ctrip7856(i5,i7)=(cfkz78(i7)%e(0)*cfkz56(i5)%e(0)-cfkz78(i
     & 7)%e(1)*cfkz56(i5)%e(1)-cfkz78(i7)%e(2)*cfkz56(i5)%e(2)-c
     & fkz78(i7)%e(3)*cfkz56(i5)%e(3))
      end do
      end do
      do i1=1,2
      do i3=1,2
* p.q -- p.q=ctrip1234(i1,i3),p=cfkz12(i1)%e,q=cfkz34(i3)%e,bef=,aft=
      ctrip1234(i1,i3)=(cfkz12(i1)%e(0)*cfkz34(i3)%e(0)-cfkz12(i
     & 1)%e(1)*cfkz34(i3)%e(1)-cfkz12(i1)%e(2)*cfkz34(i3)%e(2)-c
     & fkz12(i1)%e(3)*cfkz34(i3)%e(3))
      end do
      end do
  
*      do i1=1,2                                                        
*      do i3=1,2                                                        
*      do i5=1,2                                                        
*      do i7=1,2                                                        
  
*        cbf_tot(i1,i3,i5,i7)=(rmz2/s2w/rc2w)*                          
*     &    (ctrip7812(i1,i7)*ctrip3456(i3,i5)/(p1278q-cmh2)+            
*     &     ctrip7834(i3,i7)*ctrip1256(i1,i5)/(p1256q-cmh2)+            
*     &     ctrip7856(i5,i7)*ctrip1234(i1,i3)/(p1234q-cmh2))            
*                                                                       
*      enddo                                                            
*      enddo                                                            
*      enddo                                                            
*      enddo                                                            
  
      endif
  
*                                                                       
* RESULT                                                                
*                                                                       
      spk0=sqrt(p1k0*p2k0*p3k0*p4k0*p5k0*p6k0*p7k0*p8k0)
  
*     if (rmh.ge.0.d0) then                                             
*                                                                       
*     do i1=1,2                                                         
*     do i3=1,2                                                         
*     do i5=1,2                                                         
*     do i7=1,2                                                         
  
*      cres(i1,i3,i5,i7)=(c3fk_tot(i1,i3,i5,i7)+                        
*    &   cbct_tot(i1,i3,i5,i7)+cbf_tot(i1,i3,i5,i7))/spk0               
*                                                                       
*     enddo                                                             
*     enddo                                                             
*     enddo                                                             
*     enddo                                                             
*                                                                       
*     else                                                              
*                                                                       
*     do i1=1,2                                                         
*     do i3=1,2                                                         
*     do i5=1,2                                                         
*     do i7=1,2                                                         
  
*      cres(i1,i3,i5,i7)=(c3fk_tot(i1,i3,i5,i7)+                        
*    &   cbct_tot(i1,i3,i5,i7))/spk0                                    
*                                                                       
*     enddo                                                             
*     enddo                                                             
*     enddo                                                             
*     enddo                                                             
*                                                                       
*      endif                                                            
  
      if (rmh.ge.0.d0) then
  
* sig                                                                   
        if (i_signal.eq.0) then
          do i1=1,2
          do i3=1,2
          do i5=1,2
          do i7=1,2
  
            cres(i1,i3,i5,i7)=( c3fk12_345678(i1,i3,i5,i7)+
     &           c3fk34_125678(i3,i1,i5,i7)+c3fk56_123478(i5,i1,i3,i7)+
     &           c3fk78_123456(i7,i1,i3,i5)+
     &           cbct1234(i1,i3,i5,i7)+cbct1256(i1,i3,i5,i7)+
     &           cbct1278(i1,i3,i5,i7)+
  
     &           ctrip7812(i1,i7)*ctrip3456(i3,i5)*
     &            (rhzz**2/(p1278q-cmh2)+rhhzz**2/(p1278q-cmhh2))+
     &           ctrip7834(i3,i7)*ctrip1256(i1,i5)*
     &            (rhzz**2/(p1256q-cmh2)+rhhzz**2/(p1256q-cmhh2))+
     &           ctrip7856(i5,i7)*ctrip1234(i1,i3)*
     &            (rhzz**2/(p1234q-cmh2)+rhhzz**2/(p1234q-cmhh2)) )/spk0
          enddo
          enddo
          enddo
          enddo
        else if (i_signal.eq.1)then
          if (p1(0).lt.0.d0.and.p2(0).lt.0.d0.or.
     &        p3(0).lt.0.d0.and.p4(0).lt.0.d0.or.
     &        p5(0).lt.0.d0.and.p6(0).lt.0.d0.or.
     &        p7(0).lt.0.d0.and.p8(0).lt.0.d0)then
            do i1=1,2
            do i3=1,2
            do i5=1,2
            do i7=1,2
              cres(i1,i3,i5,i7)=czero
            enddo
            enddo
            enddo
            enddo
          else
            if (p1234q.gt.0.d0) then
              do i1=1,2
              do i3=1,2
              do i5=1,2
              do i7=1,2
                cres(i1,i3,i5,i7)= ctrip7856(i5,i7)*ctrip1234(i1,i3)*
     &            (rhzz**2/(p1234q-cmh2)+rhhzz**2/(p1234q-cmhh2))/spk0
              enddo
              enddo
              enddo
              enddo
            else if (p1256q.gt.0.d0) then
              do i1=1,2
              do i3=1,2
              do i5=1,2
              do i7=1,2
                cres(i1,i3,i5,i7)= ctrip7834(i3,i7)*ctrip1256(i1,i5)*
     &            (rhzz**2/(p1256q-cmh2)+rhhzz**2/(p1256q-cmhh2))/spk0
              enddo
              enddo
              enddo
              enddo
            else if (p1278q.gt.0.d0) then
              do i1=1,2
              do i3=1,2
              do i5=1,2
              do i7=1,2
                cres(i1,i3,i5,i7)= ctrip7812(i1,i7)*ctrip3456(i3,i5)*
     &            (rhzz**2/(p1278q-cmh2)+rhhzz**2/(p1278q-cmhh2))/spk0
              enddo
              enddo
              enddo
              enddo
            else
*              print*, 'ERROR: check if there are higgs in s channel'   
*              stop                                                     
              do i1=1,2
                do i3=1,2
                  do i5=1,2
                    do i7=1,2
                      cres(i1,i3,i5,i7)=czero
                    enddo
                  enddo
                enddo
              enddo
            endif
          endif
* six                                                                   
        else if (i_signal.eq.2)then
          if (p1(0).lt.0.d0.and.p2(0).lt.0.d0.or.
     &        p3(0).lt.0.d0.and.p4(0).lt.0.d0.or.
     &        p5(0).lt.0.d0.and.p6(0).lt.0.d0.or.
     &        p7(0).lt.0.d0.and.p8(0).lt.0.d0)then
            do i1=1,2
            do i3=1,2
            do i5=1,2
            do i7=1,2
              cres(i1,i3,i5,i7)=czero
            enddo
            enddo
            enddo
            enddo
          else
              do i1=1,2
              do i3=1,2
              do i5=1,2
              do i7=1,2
                cres(i1,i3,i5,i7)=
     &              ( ctrip7812(i1,i7)*ctrip3456(i3,i5)*
     &                (rhzz**2/(p1278q-cmh2)+rhhzz**2/(p1278q-cmhh2))+
     &                ctrip7834(i3,i7)*ctrip1256(i1,i5)*
     &                (rhzz**2/(p1256q-cmh2)+rhhzz**2/(p1256q-cmhh2))+
     &                ctrip7856(i5,i7)*ctrip1234(i1,i3)*
     &                (rhzz**2/(p1234q-cmh2)+rhhzz**2/(p1234q-cmhh2))
     &                )/spk0
              enddo
              enddo
              enddo
              enddo
          endif
  
        else if (i_signal.eq.3)then
          do i1=1,2
            do i3=1,2
              do i5=1,2
                do i7=1,2
                cres(i1,i3,i5,i7)=
     &              ( ctrip7812(i1,i7)*ctrip3456(i3,i5)*
     &                (rhzz**2/(p1278q-cmh2)+rhhzz**2/(p1278q-cmhh2))+
     &                ctrip7834(i3,i7)*ctrip1256(i1,i5)*
     &                (rhzz**2/(p1256q-cmh2)+rhhzz**2/(p1256q-cmhh2))+
     &                ctrip7856(i5,i7)*ctrip1234(i1,i3)*
     &                (rhzz**2/(p1234q-cmh2)+rhhzz**2/(p1234q-cmhh2))
     &                )/spk0
                enddo
              enddo
            enddo
          enddo
  
* sixend                                                                
        else
            print*, 'ERROR: check i_signal value'
            stop
        endif
  
* sigend                                                                
  
      else
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
  
       cres(i1,i3,i5,i7)=( c3fk12_345678(i1,i3,i5,i7)+
     &   c3fk34_125678(i3,i1,i5,i7)+c3fk56_123478(i5,i1,i3,i7)+
     &   c3fk78_123456(i7,i1,i3,i5)+
     &   cbct1234(i1,i3,i5,i7)+cbct1256(i1,i3,i5,i7)+
     &   cbct1278(i1,i3,i5,i7) )/spk0
  
      enddo
      enddo
      enddo
      enddo
  
  
  
      endif
  
*diogo 10/2009 unitarization part                                       
  
*                                                                       
* ZZZZ Quartic effective coupling and unitarization                     
*                                                                       
  
      if (i_unitarize.eq.1) then
  
       call unit_4z_massless(p12,s12,cfkz12,p34,s34,cfkz34,
     &    p56,s56,cfkz56,p78,s78,cfkz78,cunit)
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
  
       cres(i1,i3,i5,i7)=cres(i1,i3,i5,i7) + cunit(i1,i3,i5,i7)/spk0
  
      enddo
      enddo
      enddo
      enddo
  
      endif
*end diogo                                                              
  
      return
      end
  
