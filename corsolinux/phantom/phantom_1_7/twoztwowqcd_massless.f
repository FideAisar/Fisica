************************************************************************
*                                                                       
* This routine computes the basic amplitude corresponding to 8          
*  outgoing  fermions which can form 2 Z and 2 W-.                      
*  All fermions are massless.                                           
*  All outgoing particles are considered different.                     
*  The input particles are ordered from 1 to 8 in such a way that       
*  odd are particles, even antiparticles. 12 corresponds to Z,          
*  34 to Z, 56 to W+, 78 to W-                                          
*  p1, ...p8 are all outgoing momenta                                   
*  id1....id8 give the identities of the outgoing particles             
*  cres  is on output the complex helicity amplitude                    
*  cres=cres(2,2)                                                       
*   the  two indeces refer to the chiralities of the two Z lines,       
*   while the W ones have only                                          
*   chirality indeces  = 2   (corresponding to -)                       
*                                                                       
************************************************************************
*                                                                       
*cres(1):   1| 3| 5| 7|   (pure electroweak contribution)               
*             |  |  |  |                                                
*            2| 4| 6| 8|                                                
*                                                                       
*                                                                       
* cres(2):   5| 1| 3| 7|   (QCD contribution)                           
*             |  |~~|  |                                                
*            6| 2| 4| 8|                                                
*                                                                       
*                                                                       
* cres(3):   3| 1| 5| 7|   (QCD contribution)                           
*             |  |~~|  |                                                
*            4| 2| 6| 8|                                                
*                                                                       
*                                                                       
* cres(4):   1| 3| 5| 7|   (QCD contribution)                           
*             |  |~~|  |                                                
*            2| 4| 6| 8|                                                
*                                                                       
*                                                                       
* cres(5):   3| 1| 7| 5|   (QCD contribution)                           
*             |  |~~|  |                                                
*            4| 2| 8| 6|                                                
*                                                                       
*                                                                       
* cres(6):   1| 3| 7| 5|   (QCD contribution)                           
*             |  |~~|  |                                                
*            2| 4| 8| 6|                                                
*                                                                       
*                                                                       
* Notice that the 7th configuration                                     
*                                                                       
*            1| 5| 7| 3|   (QCD contribution)                           
*             |  |~~|  |                                                
*            2| 6| 8| 4|                                                
*                                                                       
* Comments:  "W" and "Z"  stay for, respectivecly, a couple of fermion  
*   which can form a W  and a couple of fermion which can form a Z      
*   Correpsondingly Wline and Zline stay for a fermion line which       
*    ends with two quarks which can form a W or a Z                     
*                                                                       
  
  
      subroutine twoztwowqcd_massless(p1,p2,p3,p4,p5,p6,p7,p8,
     &              id1,id2,id3,id4,id5,id6,id7,id8,cres)
  
      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)
  
* FOUR MOMENTA                                                          
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),
     &          p5(0:3),p6(0:3),p7(0:3),p8(0:3)
* forks                                                                 
      dimension p12(0:3),p34(0:3),p56(0:3),p78(0:3)
* left                                                                  
      dimension p134(0:3),p156(0:3),p178(0:3),
     &          p312(0:3),p356(0:3),p378(0:3),
     &          p512(0:3),p534(0:3),p578(0:3),
     &          p712(0:3),p734(0:3),p756(0:3)
* right                                                                 
      dimension p234(0:3),p256(0:3),p278(0:3),
     &          p412(0:3),p456(0:3),p478(0:3),
     &          p612(0:3),p634(0:3),p678(0:3),
     &          p812(0:3),p834(0:3),p856(0:3)
* four2four                                                             
      dimension p1234(0:3),p1256(0:3),p1278(0:3),
     &   p3456(0:3),p3478(0:3),p5678(0:3)
  
* triple vertices auxiliary                                             
      dimension vfz(0:3), vwm(0:3), vwp(0:3)
  
* LINES                                                                 
*forks (vector)                                                         
      type pol
         double complex e(0:3),ek0,v
      end type
      type(pol) cz12(2),cf12(2),cz34(2),cf34(2),cw56,cw78
* left                                                                  
      type l_line
         double complex a(2),c(2)
      end type
      type(l_line) l1_34(2),l1_56   ,l1_78,
     &              l3_12(2),l3_56   ,l3_78,
     &              l5_12(2),l5_34(2),l5_78,
     &              l7_12(2),l7_34(2),l7_56,
**QCD                                                                   
     &              lg1_34(2),lg3_12(2),lg5_12(2),lg5_34(2),
     &              lg7_12(2),lg7_34(2)
*right                                                                  
      type r_line
         double complex a(2),b(2)
      end type
      type(r_line) r2_34(2),r2_56   ,r2_78,
     &              r4_12(2),r4_56   ,r4_78,
     &              r6_12(2),r6_34(2),r6_78,
     &              r8_12(2),r8_34(2),r8_56,
**QCD                                                                   
     &              rg2_34(2),rg4_12(2),rg6_12(2),rg6_34(2),
     &              rg8_12(2),rg8_34(2)
*middle                                                                 
      type u_line
         double complex a(2),b(2),c(2),d(2)
      end type
      type(u_line) u134_56   ,u134_78   ,u156_34(2),
     &              u156_78   ,u178_34(2),u178_56   ,
     &              u312_56   ,u312_78   ,u356_12(2),
     &              u356_78   ,u378_12(2),u378_56   ,
     &              u512_34(2),u512_78   ,u534_12(2),
     &              u534_78   ,u578_12(2),u578_34(2),
     &              u712_34(2),u712_56   ,u734_12(2),
     &              u734_56   ,u756_12(2),u756_34(2),
**QCD                                                                   
     &              ug156_34(2),ug178_34(2),ug356_12(2),
     &              ug378_12(2),ug512_34(2),ug534_12(2),
     &              ug578_12(2),ug578_34(2),ug712_34(2),
     &              ug734_12(2),ug756_12(2),ug756_34(2)
  
*left*middle = left line  (2-forks)                                     
      type(l_line) l1_3456(2)  ,l1_3478(2)  ,l1_5634(2)  ,
     &              l1_5678     ,l1_7834(2)  ,l1_7856     ,
     &              l3_1256(2)  ,l3_1278(2)  ,l3_5612(2)  ,
     &              l3_5678     ,l3_7812(2)  ,l3_7856     ,
     &              l5_1234(2,2),l5_1278(2)  ,l5_3412(2,2),
     &              l5_3478(2)  ,l5_7812(2)  ,l5_7834(2)  ,
     &              l7_1234(2,2),l7_1256(2)  ,l7_3412(2,2),
     &              l7_3456(2)  ,l7_5612(2)  ,l7_5634(2),
**QCD                                                                   
     &              lg1_3456(2)  ,lg1_3478(2)  ,lg3_1256(2)  ,
     &              lg3_1278(2)  ,lg5_1234(2,2),lg5_1278(2)  ,
     &              lg5_3412(2,2),lg5_3478(2)  ,lg7_1234(2,2),
     &              lg7_1256(2)  ,lg7_3412(2,2),lg7_3456(2)
  
*left*midle*right = a (3-forks)                                         
      dimension c3fk12_345678(2)  ,c3fk34_125678(2)  ,
     &          c3fk56_123478(2,2),c3fk78_123456(2,2),c3fk_tot(2,2),
**QCD                                                                   
     &          cg3fk12_345678(2),cg3fk34_125678(2),cg3fk1234(2,2),
     &          cg3fk56_123478(2,2),cg3fk56_341278(2,2),
     &          cg3fk78_123456(2,2),cg3fk78_341256(2,2)
  
* left vector                                                           
      type(l_line) lz1_34(0:3),lf1_34(0:3),lw1_56(0:3),lw1_78(0:3),
     &              lz3_12(0:3),lf3_12(0:3),lw3_56(0:3),lw3_78(0:3),
     &              lw5_12(0:3),lw5_34(0:3),lz5_78(0:3),lf5_78(0:3),
     &              lw7_12(0:3),lw7_34(0:3),lz7_56(0:3),lf7_56(0:3)
* right vector                                                          
      type(r_line) rz2_34(0:3),rf2_34(0:3),rw2_56(0:3),rw2_78(0:3),
     &              rz4_12(0:3),rf4_12(0:3),rw4_56(0:3),rw4_78(0:3),
     &              rw6_12(0:3),rw6_34(0:3),rz6_78(0:3),rf6_78(0:3),
     &              rw8_12(0:3),rw8_34(0:3),rz8_56(0:3),rf8_56(0:3)
* left*right vector                                                     
      type(pol) cz1234(2,2),cf1234(2,2),cw1256(2),cw1278(2),
     &           cz5678,     cf5678,     cw3478(2),cw3456(2),
**QCD                                                                   
     &           cgz1234(2,2),cgf1234(2,2),cgw5612(2),cgw7812(2),
     &           cgw7834(2),cgw5634(2),
     &           cg1234(2,2),cg3412(2,2),cg5678,cg7856
  
      dimension ch1234(2,2)   ! ch5678
  
* triple vertex auxiliary vector                                        
      dimension cwaux(2,0:3),c4aux(2,2),c4aux_2(2,2),
     &          c2aux1(2),c2aux2(2),cvaux(0:3)
  
* four2four BOSON CONNECTION                                            
      dimension cbct1234(2,2),cbct1256(2,2),cbct1278(2,2),cquart(2,2)
**QCD                                                                   
      dimension cgbct1234(2,2),cgbct1256(2,2),cgbct1278(2,2),
     &          cgbct3456(2,2),cgbct3478(2,2),
* gluon connection                                                      
     &          cgct1256(2,2),cgct1278(2,2),
     &          cgct3456(2,2),cgct3478(2,2)
  
  
      dimension cres(2,2,6)
  
      include 'common.h'
  
* pk0 of single momenta                                                 
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
  
*                                                                       
* FORKS, with propagators                                               
*                                                                       
  
*forks momenta                                                          
      do m=0,3
        p12(m)=p1(m)+p2(m)
        p34(m)=p3(m)+p4(m)
        p56(m)=p5(m)+p6(m)
        p78(m)=p7(m)+p8(m)
      enddo
* Z forks                                                               
* quqd -- p=p1,q=p2
      quqd=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      s12=2.d0*quqd
  
      ccr=zcr(id1)/(-s12+cmz2)
      ccl=zcl(id1)/(-s12+cmz2)
* T10 -- qu=p1,qd=p2,v=0,a=cz12(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p1(2)*p2(3)+p2(2)*p1(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p1k0*p2(0)+p2k0*p1(0)
      cz12(1)%e(0)=ccr*(auxa+ceps_0)
      cz12(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=1,a=cz12(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p1k0*p2(1)+p2k0*p1(1)
      cz12(1)%e(1)=ccr*(auxa+ceps_0)
      cz12(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=2,a=cz12(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p1k0*p2(3)+p2k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(2)+p2k0*p1(2)
      cz12(1)%e(2)=ccr*(auxa+ceps_0)
      cz12(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=3,a=cz12(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p1k0*p2(2)-p2k0*p1(2)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(3)+p2k0*p1(3)
      cz12(1)%e(3)=ccr*(auxa+ceps_0)
      cz12(2)%e(3)=ccl*(auxa-ceps_0)
  
      do i1=1,2
* pk0 -- p=cz12(i1)%e
      cz12(i1)%ek0=cz12(i1)%e(0)-cz12(i1)%e(1)
      end do
  
      ccr=fcr(id1)/(-s12)
      ccl=fcl(id1)/(-s12)
* T10 -- qu=p1,qd=p2,v=0,a=cf12(&)%e(0),cr=ccr,cl=ccl,nsum=0
      eps_0=-p1(2)*p2(3)+p2(2)*p1(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p1k0*p2(0)+p2k0*p1(0)
      cf12(1)%e(0)=ccr*(auxa+ceps_0)
      cf12(2)%e(0)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=1,a=cf12(&)%e(1),cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p1k0*p2(1)+p2k0*p1(1)
      cf12(1)%e(1)=ccr*(auxa+ceps_0)
      cf12(2)%e(1)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=2,a=cf12(&)%e(2),cr=ccr,cl=ccl,nsum=0
      eps_0=-p1k0*p2(3)+p2k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(2)+p2k0*p1(2)
      cf12(1)%e(2)=ccr*(auxa+ceps_0)
      cf12(2)%e(2)=ccl*(auxa-ceps_0)
* T10 -- qu=p1,qd=p2,v=3,a=cf12(&)%e(3),cr=ccr,cl=ccl,nsum=0
      eps_0=p1k0*p2(2)-p2k0*p1(2)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(3)+p2k0*p1(3)
      cf12(1)%e(3)=ccr*(auxa+ceps_0)
      cf12(2)%e(3)=ccl*(auxa-ceps_0)
  
      do i1=1,2
* pk0 -- p=cf12(i1)%e
      cf12(i1)%ek0=cf12(i1)%e(0)-cf12(i1)%e(1)
      end do
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
  
      do i1=1,2
* pk0 -- p=cz34(i1)%e
      cz34(i1)%ek0=cz34(i1)%e(0)-cz34(i1)%e(1)
      end do
  
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
  
      do i1=1,2
* pk0 -- p=cf34(i1)%e
      cf34(i1)%ek0=cf34(i1)%e(0)-cf34(i1)%e(1)
      end do
  
* W forks                                                               
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
  
*                                                                       
* LEFT INSERTIONS, with fermion propagator and pk0                      
*                                                                       
  
* main line Z (12 or 34)                                                
  
*    W insertion                                                        
  
      if (iup(id1).eq.1) then !W = W- (78)
       do m=0,3
        p178(m) = p1(m)+p78(m)
       enddo
* pk0 -- p=p178
      p178k0=p178(0)-p178(1)
* p.q -- p.q=p178q,p=p178,q=p178,bef=,aft=
      p178q=(p178(0)*p178(0)-p178(1)*p178(1)-p178(2)*p178(2)-p17
     & 8(3)*p178(3))
  
* quqd -- p=p1,q=p178
      quqd=p1(0)*p178(0)-p1(1)*p178(1)-p1(2)*p178(2)-p1(3)*p178(
     & 3)
      ccl=wcl/(p178q*p178k0)
* TWL0 -- qu=p1,qd=p178,v=cw78%e,a=l1_78%a,c=l1_78%c,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p1(2)*p178(3)-p178(2)*p1(3))+p1k0*(cw78%
     & e(2)*p178(3)-p178(2)*cw78%e(3))-p178k0*(cw78%e(2)*p1(3)-p
     & 1(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p1k0+p1(3)*cw78%ek0
      ceps_1=ceps_1*cim
      cvqu=cw78%e(0)*p1(0)-cw78%e(1)*p1(1)-cw78%e(2)*p1(2)-cw78%
     & e(3)*p1(3)
      cvqd=cw78%e(0)*p178(0)-cw78%e(1)*p178(1)-cw78%e(2)*p178(2)
     & -cw78%e(3)*p178(3)
      cauxa=-cw78%ek0*quqd+p1k0*cvqd+p178k0*cvqu
      cauxc=+cw78%ek0*p1(2)-p1k0*cw78%e(2)
      l1_78%a(2)=ccl*(cauxa-ceps_0)
      l1_78%c(2)=ccl*(-cauxc+ceps_1)
  
      else !W = W+ (56)
       do m=0,3
        p156(m) = p1(m)+p56(m)
       enddo
* pk0 -- p=p156
      p156k0=p156(0)-p156(1)
* p.q -- p.q=p156q,p=p156,q=p156,bef=,aft=
      p156q=(p156(0)*p156(0)-p156(1)*p156(1)-p156(2)*p156(2)-p15
     & 6(3)*p156(3))
  
* quqd -- p=p1,q=p156
      quqd=p1(0)*p156(0)-p1(1)*p156(1)-p1(2)*p156(2)-p1(3)*p156(
     & 3)
      ccl=wcl/(p156q*p156k0)
* TWL0 -- qu=p1,qd=p156,v=cw56%e,a=l1_56%a,c=l1_56%c,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p1(2)*p156(3)-p156(2)*p1(3))+p1k0*(cw56%
     & e(2)*p156(3)-p156(2)*cw56%e(3))-p156k0*(cw56%e(2)*p1(3)-p
     & 1(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p1k0+p1(3)*cw56%ek0
      ceps_1=ceps_1*cim
      cvqu=cw56%e(0)*p1(0)-cw56%e(1)*p1(1)-cw56%e(2)*p1(2)-cw56%
     & e(3)*p1(3)
      cvqd=cw56%e(0)*p156(0)-cw56%e(1)*p156(1)-cw56%e(2)*p156(2)
     & -cw56%e(3)*p156(3)
      cauxa=-cw56%ek0*quqd+p1k0*cvqd+p156k0*cvqu
      cauxc=+cw56%ek0*p1(2)-p1k0*cw56%e(2)
      l1_56%a(2)=ccl*(cauxa-ceps_0)
      l1_56%c(2)=ccl*(-cauxc+ceps_1)
  
      endif
  
      if (iup(id3).eq.1) then !W = W- (78)
       do m=0,3
        p378(m) = p3(m)+p78(m)
       enddo
* pk0 -- p=p378
      p378k0=p378(0)-p378(1)
* p.q -- p.q=p378q,p=p378,q=p378,bef=,aft=
      p378q=(p378(0)*p378(0)-p378(1)*p378(1)-p378(2)*p378(2)-p37
     & 8(3)*p378(3))
  
* quqd -- p=p3,q=p378
      quqd=p3(0)*p378(0)-p3(1)*p378(1)-p3(2)*p378(2)-p3(3)*p378(
     & 3)
      ccl=wcl/(p378q*p378k0)
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
  
      else !W = W+ (56)
       do m=0,3
        p356(m) = p3(m)+p56(m)
       enddo
* pk0 -- p=p356
      p356k0=p356(0)-p356(1)
* p.q -- p.q=p356q,p=p356,q=p356,bef=,aft=
      p356q=(p356(0)*p356(0)-p356(1)*p356(1)-p356(2)*p356(2)-p35
     & 6(3)*p356(3))
  
* quqd -- p=p3,q=p356
      quqd=p3(0)*p356(0)-p3(1)*p356(1)-p3(2)*p356(2)-p3(3)*p356(
     & 3)
      ccl=wcl/(p356q*p356k0)
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
  
*     Z insertion                                                       
       do m=0,3
        p134(m) = p1(m)+p34(m)
       enddo
* pk0 -- p=p134
      p134k0=p134(0)-p134(1)
* p.q -- p.q=p134q,p=p134,q=p134,bef=,aft=
      p134q=(p134(0)*p134(0)-p134(1)*p134(1)-p134(2)*p134(2)-p13
     & 4(3)*p134(3))
  
  
      if (ineutri(id1).ne.1.and.ineutri(id3).ne.1) then
  
* quqd -- p=p1,q=p134
      quqd=p1(0)*p134(0)-p1(1)*p134(1)-p1(2)*p134(2)-p1(3)*p134(
     & 3)
      ccr=fcr(id1)/(p134q*p134k0)
      ccl=fcl(id1)/(p134q*p134k0)
      do i3=1,2
* TL0 -- qu=p1,qd=p134,v=cf34(i3)%e,a=l1_34(i3)%a,c=l1_34(i3)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cf34(i3)%ek0*(p1(2)*p134(3)-p134(2)*p1(3))+p1k0*(c
     & f34(i3)%e(2)*p134(3)-p134(2)*cf34(i3)%e(3))-p134k0*(cf34(
     & i3)%e(2)*p1(3)-p1(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p1k0+p1(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf34(i3)%e(0)*p1(0)-cf34(i3)%e(1)*p1(1)-cf34(i3)%e(2)
     & *p1(2)-cf34(i3)%e(3)*p1(3)
      cvqd=cf34(i3)%e(0)*p134(0)-cf34(i3)%e(1)*p134(1)-cf34(i3)%
     & e(2)*p134(2)-cf34(i3)%e(3)*p134(3)
      cauxa=-cf34(i3)%ek0*quqd+p1k0*cvqd+p134k0*cvqu
      cauxc=+cf34(i3)%ek0*p1(2)-p1k0*cf34(i3)%e(2)
      l1_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      l1_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      l1_34(i3)%c(1)=ccr*(cauxc+ceps_1)
      l1_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
**QCD                                                                   
       if (ilept(id1).ne.1.and.ilept(id3).ne.1) then
         do iut=1,2
            lg1_34(iut)%a(1) = l1_34(iut)%a(1)
     &           /((fcl(id1))*(fcl(id3)))
            lg1_34(iut)%a(2) = l1_34(iut)%a(2)
     &           /((fcl(id1))*(fcl(id3)))
            lg1_34(iut)%c(1) = l1_34(iut)%c(1)
     &           /((fcl(id1))*(fcl(id3)))
            lg1_34(iut)%c(2) = l1_34(iut)%c(2)
     &           /((fcl(id1))*(fcl(id3)))
         enddo
       endif
  
      ccr=zcr(id1)/(p134q*p134k0)
      ccl=zcl(id1)/(p134q*p134k0)
      do i3=1,2
* TL0 -- qu=p1,qd=p134,v=cz34(i3)%e,a=l1_34(i3)%a,c=l1_34(i3)%c,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cz34(i3)%ek0*(p1(2)*p134(3)-p134(2)*p1(3))+p1k0*(c
     & z34(i3)%e(2)*p134(3)-p134(2)*cz34(i3)%e(3))-p134k0*(cz34(
     & i3)%e(2)*p1(3)-p1(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p1k0+p1(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i3)%e(0)*p1(0)-cz34(i3)%e(1)*p1(1)-cz34(i3)%e(2)
     & *p1(2)-cz34(i3)%e(3)*p1(3)
      cvqd=cz34(i3)%e(0)*p134(0)-cz34(i3)%e(1)*p134(1)-cz34(i3)%
     & e(2)*p134(2)-cz34(i3)%e(3)*p134(3)
      cauxa=-cz34(i3)%ek0*quqd+p1k0*cvqd+p134k0*cvqu
      cauxc=+cz34(i3)%ek0*p1(2)-p1k0*cz34(i3)%e(2)
      l1_34(i3)%a(1)=l1_34(i3)%a(1)+ccr*(cauxa+ceps_0)
      l1_34(i3)%a(2)=l1_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      l1_34(i3)%c(1)=l1_34(i3)%c(1)+ccr*(cauxc+ceps_1)
      l1_34(i3)%c(2)=l1_34(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      else
* quqd -- p=p1,q=p134
      quqd=p1(0)*p134(0)-p1(1)*p134(1)-p1(2)*p134(2)-p1(3)*p134(
     & 3)
      ccr=zcr(id1)/(p134q*p134k0)
      ccl=zcl(id1)/(p134q*p134k0)
      do i3=1,2
* TL0 -- qu=p1,qd=p134,v=cz34(i3)%e,a=l1_34(i3)%a,c=l1_34(i3)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p1(2)*p134(3)-p134(2)*p1(3))+p1k0*(c
     & z34(i3)%e(2)*p134(3)-p134(2)*cz34(i3)%e(3))-p134k0*(cz34(
     & i3)%e(2)*p1(3)-p1(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p1k0+p1(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i3)%e(0)*p1(0)-cz34(i3)%e(1)*p1(1)-cz34(i3)%e(2)
     & *p1(2)-cz34(i3)%e(3)*p1(3)
      cvqd=cz34(i3)%e(0)*p134(0)-cz34(i3)%e(1)*p134(1)-cz34(i3)%
     & e(2)*p134(2)-cz34(i3)%e(3)*p134(3)
      cauxa=-cz34(i3)%ek0*quqd+p1k0*cvqd+p134k0*cvqu
      cauxc=+cz34(i3)%ek0*p1(2)-p1k0*cz34(i3)%e(2)
      l1_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      l1_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      l1_34(i3)%c(1)=ccr*(cauxc+ceps_1)
      l1_34(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
  
       do m=0,3
        p312(m) = p3(m)+p12(m)
       enddo
* pk0 -- p=p312
      p312k0=p312(0)-p312(1)
* p.q -- p.q=p312q,p=p312,q=p312,bef=,aft=
      p312q=(p312(0)*p312(0)-p312(1)*p312(1)-p312(2)*p312(2)-p31
     & 2(3)*p312(3))
  
  
      if (ineutri(id1).ne.1.and.ineutri(id3).ne.1) then
  
* quqd -- p=p3,q=p312
      quqd=p3(0)*p312(0)-p3(1)*p312(1)-p3(2)*p312(2)-p3(3)*p312(
     & 3)
      ccr=fcr(id3)/(p312q*p312k0)
      ccl=fcl(id3)/(p312q*p312k0)
      do i3=1,2
* TL0 -- qu=p3,qd=p312,v=cf12(i3)%e,a=l3_12(i3)%a,c=l3_12(i3)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cf12(i3)%ek0*(p3(2)*p312(3)-p312(2)*p3(3))+p3k0*(c
     & f12(i3)%e(2)*p312(3)-p312(2)*cf12(i3)%e(3))-p312k0*(cf12(
     & i3)%e(2)*p3(3)-p3(2)*cf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i3)%e(3)*p3k0+p3(3)*cf12(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf12(i3)%e(0)*p3(0)-cf12(i3)%e(1)*p3(1)-cf12(i3)%e(2)
     & *p3(2)-cf12(i3)%e(3)*p3(3)
      cvqd=cf12(i3)%e(0)*p312(0)-cf12(i3)%e(1)*p312(1)-cf12(i3)%
     & e(2)*p312(2)-cf12(i3)%e(3)*p312(3)
      cauxa=-cf12(i3)%ek0*quqd+p3k0*cvqd+p312k0*cvqu
      cauxc=+cf12(i3)%ek0*p3(2)-p3k0*cf12(i3)%e(2)
      l3_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      l3_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      l3_12(i3)%c(1)=ccr*(cauxc+ceps_1)
      l3_12(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
**QCD                                                                   
       if (ilept(id1).ne.1.and.ilept(id3).ne.1) then
         do iut=1,2
            lg3_12(iut)%a(1) = l3_12(iut)%a(1)
     &           /((fcl(id1))*(fcl(id3)))
            lg3_12(iut)%a(2) = l3_12(iut)%a(2)
     &           /((fcl(id1))*(fcl(id3)))
            lg3_12(iut)%c(1) = l3_12(iut)%c(1)
     &           /((fcl(id1))*(fcl(id3)))
            lg3_12(iut)%c(2) = l3_12(iut)%c(2)
     &           /((fcl(id1))*(fcl(id3)))
         enddo
       endif
  
      ccr=zcr(id3)/(p312q*p312k0)
      ccl=zcl(id3)/(p312q*p312k0)
      do i3=1,2
* TL0 -- qu=p3,qd=p312,v=cz12(i3)%e,a=l3_12(i3)%a,c=l3_12(i3)%c,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cz12(i3)%ek0*(p3(2)*p312(3)-p312(2)*p3(3))+p3k0*(c
     & z12(i3)%e(2)*p312(3)-p312(2)*cz12(i3)%e(3))-p312k0*(cz12(
     & i3)%e(2)*p3(3)-p3(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p3k0+p3(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz12(i3)%e(0)*p3(0)-cz12(i3)%e(1)*p3(1)-cz12(i3)%e(2)
     & *p3(2)-cz12(i3)%e(3)*p3(3)
      cvqd=cz12(i3)%e(0)*p312(0)-cz12(i3)%e(1)*p312(1)-cz12(i3)%
     & e(2)*p312(2)-cz12(i3)%e(3)*p312(3)
      cauxa=-cz12(i3)%ek0*quqd+p3k0*cvqd+p312k0*cvqu
      cauxc=+cz12(i3)%ek0*p3(2)-p3k0*cz12(i3)%e(2)
      l3_12(i3)%a(1)=l3_12(i3)%a(1)+ccr*(cauxa+ceps_0)
      l3_12(i3)%a(2)=l3_12(i3)%a(2)+ccl*(cauxa-ceps_0)
      l3_12(i3)%c(1)=l3_12(i3)%c(1)+ccr*(cauxc+ceps_1)
      l3_12(i3)%c(2)=l3_12(i3)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      else
* quqd -- p=p3,q=p312
      quqd=p3(0)*p312(0)-p3(1)*p312(1)-p3(2)*p312(2)-p3(3)*p312(
     & 3)
      ccr=zcr(id3)/(p312q*p312k0)
      ccl=zcl(id3)/(p312q*p312k0)
      do i3=1,2
* TL0 -- qu=p3,qd=p312,v=cz12(i3)%e,a=l3_12(i3)%a,c=l3_12(i3)%c,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz12(i3)%ek0*(p3(2)*p312(3)-p312(2)*p3(3))+p3k0*(c
     & z12(i3)%e(2)*p312(3)-p312(2)*cz12(i3)%e(3))-p312k0*(cz12(
     & i3)%e(2)*p3(3)-p3(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p3k0+p3(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz12(i3)%e(0)*p3(0)-cz12(i3)%e(1)*p3(1)-cz12(i3)%e(2)
     & *p3(2)-cz12(i3)%e(3)*p3(3)
      cvqd=cz12(i3)%e(0)*p312(0)-cz12(i3)%e(1)*p312(1)-cz12(i3)%
     & e(2)*p312(2)-cz12(i3)%e(3)*p312(3)
      cauxa=-cz12(i3)%ek0*quqd+p3k0*cvqd+p312k0*cvqu
      cauxc=+cz12(i3)%ek0*p3(2)-p3k0*cz12(i3)%e(2)
      l3_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      l3_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      l3_12(i3)%c(1)=ccr*(cauxc+ceps_1)
      l3_12(i3)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
  
  
* main line W (56 or 78)                                                
  
*     Z insertion                                                       
       do m=0,3
        p512(m) = p5(m)+p12(m)
       enddo
* pk0 -- p=p512
      p512k0=p512(0)-p512(1)
* p.q -- p.q=p512q,p=p512,q=p512,bef=,aft=
      p512q=(p512(0)*p512(0)-p512(1)*p512(1)-p512(2)*p512(2)-p51
     & 2(3)*p512(3))
  
      if (ineutri(id5).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p5,q=p512
      quqd=p5(0)*p512(0)-p5(1)*p512(1)-p5(2)*p512(2)-p5(3)*p512(
     & 3)
      ccl=fcl(id5)/(p512q*p512k0)
      do i1=1,2
* TWL0 -- qu=p5,qd=p512,v=cf12(i1)%e,a=l5_12(i1)%a,c=l5_12(i1)%c,cl=ccl,
* nsum=0
      ceps_0=-cf12(i1)%ek0*(p5(2)*p512(3)-p512(2)*p5(3))+p5k0*(c
     & f12(i1)%e(2)*p512(3)-p512(2)*cf12(i1)%e(3))-p512k0*(cf12(
     & i1)%e(2)*p5(3)-p5(2)*cf12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1)%e(3)*p5k0+p5(3)*cf12(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf12(i1)%e(0)*p5(0)-cf12(i1)%e(1)*p5(1)-cf12(i1)%e(2)
     & *p5(2)-cf12(i1)%e(3)*p5(3)
      cvqd=cf12(i1)%e(0)*p512(0)-cf12(i1)%e(1)*p512(1)-cf12(i1)%
     & e(2)*p512(2)-cf12(i1)%e(3)*p512(3)
      cauxa=-cf12(i1)%ek0*quqd+p5k0*cvqd+p512k0*cvqu
      cauxc=+cf12(i1)%ek0*p5(2)-p5k0*cf12(i1)%e(2)
      l5_12(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_12(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
**QCD                                                                   
       if (ilept(id5).ne.1.and.ilept(id1).ne.1) then
         do iut=1,2
            lg5_12(iut)%a(1) = l5_12(iut)%a(1)
     &           /((fcl(id5))*(fcl(id1)))
            lg5_12(iut)%a(2) = l5_12(iut)%a(2)
     &           /((fcl(id5))*(fcl(id1)))
            lg5_12(iut)%c(1) = l5_12(iut)%c(1)
     &           /((fcl(id5))*(fcl(id1)))
            lg5_12(iut)%c(2) = l5_12(iut)%c(2)
     &           /((fcl(id5))*(fcl(id1)))
         enddo
       endif
  
      ccl=zcl(id5)/(p512q*p512k0)
      do i1=1,2
* TWL0 -- qu=p5,qd=p512,v=cz12(i1)%e,a=l5_12(i1)%a,c=l5_12(i1)%c,cl=ccl,
* nsum=1
      ceps_0=-cz12(i1)%ek0*(p5(2)*p512(3)-p512(2)*p5(3))+p5k0*(c
     & z12(i1)%e(2)*p512(3)-p512(2)*cz12(i1)%e(3))-p512k0*(cz12(
     & i1)%e(2)*p5(3)-p5(2)*cz12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1)%e(3)*p5k0+p5(3)*cz12(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz12(i1)%e(0)*p5(0)-cz12(i1)%e(1)*p5(1)-cz12(i1)%e(2)
     & *p5(2)-cz12(i1)%e(3)*p5(3)
      cvqd=cz12(i1)%e(0)*p512(0)-cz12(i1)%e(1)*p512(1)-cz12(i1)%
     & e(2)*p512(2)-cz12(i1)%e(3)*p512(3)
      cauxa=-cz12(i1)%ek0*quqd+p5k0*cvqd+p512k0*cvqu
      cauxc=+cz12(i1)%ek0*p5(2)-p5k0*cz12(i1)%e(2)
      l5_12(i1)%a(2)=l5_12(i1)%a(2)+ccl*(cauxa-ceps_0)
      l5_12(i1)%c(2)=l5_12(i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      else
* quqd -- p=p5,q=p512
      quqd=p5(0)*p512(0)-p5(1)*p512(1)-p5(2)*p512(2)-p5(3)*p512(
     & 3)
      ccl=zcl(id5)/(p512q*p512k0)
      do i1=1,2
* TWL0 -- qu=p5,qd=p512,v=cz12(i1)%e,a=l5_12(i1)%a,c=l5_12(i1)%c,cl=ccl,
* nsum=0
      ceps_0=-cz12(i1)%ek0*(p5(2)*p512(3)-p512(2)*p5(3))+p5k0*(c
     & z12(i1)%e(2)*p512(3)-p512(2)*cz12(i1)%e(3))-p512k0*(cz12(
     & i1)%e(2)*p5(3)-p5(2)*cz12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1)%e(3)*p5k0+p5(3)*cz12(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz12(i1)%e(0)*p5(0)-cz12(i1)%e(1)*p5(1)-cz12(i1)%e(2)
     & *p5(2)-cz12(i1)%e(3)*p5(3)
      cvqd=cz12(i1)%e(0)*p512(0)-cz12(i1)%e(1)*p512(1)-cz12(i1)%
     & e(2)*p512(2)-cz12(i1)%e(3)*p512(3)
      cauxa=-cz12(i1)%ek0*quqd+p5k0*cvqd+p512k0*cvqu
      cauxc=+cz12(i1)%ek0*p5(2)-p5k0*cz12(i1)%e(2)
      l5_12(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_12(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
  
       do m=0,3
        p534(m) = p5(m)+p34(m)
       enddo
* pk0 -- p=p534
      p534k0=p534(0)-p534(1)
* p.q -- p.q=p534q,p=p534,q=p534,bef=,aft=
      p534q=(p534(0)*p534(0)-p534(1)*p534(1)-p534(2)*p534(2)-p53
     & 4(3)*p534(3))
  
      if (ineutri(id5).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p5,q=p534
      quqd=p5(0)*p534(0)-p5(1)*p534(1)-p5(2)*p534(2)-p5(3)*p534(
     & 3)
      ccl=fcl(id5)/(p534q*p534k0)
      do i1=1,2
* TWL0 -- qu=p5,qd=p534,v=cf34(i1)%e,a=l5_34(i1)%a,c=l5_34(i1)%c,cl=ccl,
* nsum=0
      ceps_0=-cf34(i1)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*(c
     & f34(i1)%e(2)*p534(3)-p534(2)*cf34(i1)%e(3))-p534k0*(cf34(
     & i1)%e(2)*p5(3)-p5(2)*cf34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1)%e(3)*p5k0+p5(3)*cf34(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf34(i1)%e(0)*p5(0)-cf34(i1)%e(1)*p5(1)-cf34(i1)%e(2)
     & *p5(2)-cf34(i1)%e(3)*p5(3)
      cvqd=cf34(i1)%e(0)*p534(0)-cf34(i1)%e(1)*p534(1)-cf34(i1)%
     & e(2)*p534(2)-cf34(i1)%e(3)*p534(3)
      cauxa=-cf34(i1)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cf34(i1)%ek0*p5(2)-p5k0*cf34(i1)%e(2)
      l5_34(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_34(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
**QCD                                                                   
       if (ilept(id5).ne.1.and.ilept(id3).ne.1) then
         do iut=1,2
            lg5_34(iut)%a(1) = l5_34(iut)%a(1)
     &           /((fcl(id5))*(fcl(id3)))
            lg5_34(iut)%a(2) = l5_34(iut)%a(2)
     &           /((fcl(id5))*(fcl(id3)))
            lg5_34(iut)%c(1) = l5_34(iut)%c(1)
     &           /((fcl(id5))*(fcl(id3)))
            lg5_34(iut)%c(2) = l5_34(iut)%c(2)
     &           /((fcl(id5))*(fcl(id3)))
         enddo
       endif
  
      ccl=zcl(id5)/(p534q*p534k0)
      do i1=1,2
* TWL0 -- qu=p5,qd=p534,v=cz34(i1)%e,a=l5_34(i1)%a,c=l5_34(i1)%c,cl=ccl,
* nsum=1
      ceps_0=-cz34(i1)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*(c
     & z34(i1)%e(2)*p534(3)-p534(2)*cz34(i1)%e(3))-p534k0*(cz34(
     & i1)%e(2)*p5(3)-p5(2)*cz34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1)%e(3)*p5k0+p5(3)*cz34(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i1)%e(0)*p5(0)-cz34(i1)%e(1)*p5(1)-cz34(i1)%e(2)
     & *p5(2)-cz34(i1)%e(3)*p5(3)
      cvqd=cz34(i1)%e(0)*p534(0)-cz34(i1)%e(1)*p534(1)-cz34(i1)%
     & e(2)*p534(2)-cz34(i1)%e(3)*p534(3)
      cauxa=-cz34(i1)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cz34(i1)%ek0*p5(2)-p5k0*cz34(i1)%e(2)
      l5_34(i1)%a(2)=l5_34(i1)%a(2)+ccl*(cauxa-ceps_0)
      l5_34(i1)%c(2)=l5_34(i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      else
* quqd -- p=p5,q=p534
      quqd=p5(0)*p534(0)-p5(1)*p534(1)-p5(2)*p534(2)-p5(3)*p534(
     & 3)
      ccl=zcl(id5)/(p534q*p534k0)
      do i1=1,2
* TWL0 -- qu=p5,qd=p534,v=cz34(i1)%e,a=l5_34(i1)%a,c=l5_34(i1)%c,cl=ccl,
* nsum=0
      ceps_0=-cz34(i1)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*(c
     & z34(i1)%e(2)*p534(3)-p534(2)*cz34(i1)%e(3))-p534k0*(cz34(
     & i1)%e(2)*p5(3)-p5(2)*cz34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1)%e(3)*p5k0+p5(3)*cz34(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i1)%e(0)*p5(0)-cz34(i1)%e(1)*p5(1)-cz34(i1)%e(2)
     & *p5(2)-cz34(i1)%e(3)*p5(3)
      cvqd=cz34(i1)%e(0)*p534(0)-cz34(i1)%e(1)*p534(1)-cz34(i1)%
     & e(2)*p534(2)-cz34(i1)%e(3)*p534(3)
      cauxa=-cz34(i1)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cz34(i1)%ek0*p5(2)-p5k0*cz34(i1)%e(2)
      l5_34(i1)%a(2)=ccl*(cauxa-ceps_0)
      l5_34(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
  
       do m=0,3
        p712(m) = p7(m)+p12(m)
       enddo
* pk0 -- p=p712
      p712k0=p712(0)-p712(1)
* p.q -- p.q=p712q,p=p712,q=p712,bef=,aft=
      p712q=(p712(0)*p712(0)-p712(1)*p712(1)-p712(2)*p712(2)-p71
     & 2(3)*p712(3))
  
      if (ineutri(id7).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p7,q=p712
      quqd=p7(0)*p712(0)-p7(1)*p712(1)-p7(2)*p712(2)-p7(3)*p712(
     & 3)
      ccl=fcl(id7)/(p712q*p712k0)
      do i1=1,2
* TWL0 -- qu=p7,qd=p712,v=cf12(i1)%e,a=l7_12(i1)%a,c=l7_12(i1)%c,cl=ccl,
* nsum=0
      ceps_0=-cf12(i1)%ek0*(p7(2)*p712(3)-p712(2)*p7(3))+p7k0*(c
     & f12(i1)%e(2)*p712(3)-p712(2)*cf12(i1)%e(3))-p712k0*(cf12(
     & i1)%e(2)*p7(3)-p7(2)*cf12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1)%e(3)*p7k0+p7(3)*cf12(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf12(i1)%e(0)*p7(0)-cf12(i1)%e(1)*p7(1)-cf12(i1)%e(2)
     & *p7(2)-cf12(i1)%e(3)*p7(3)
      cvqd=cf12(i1)%e(0)*p712(0)-cf12(i1)%e(1)*p712(1)-cf12(i1)%
     & e(2)*p712(2)-cf12(i1)%e(3)*p712(3)
      cauxa=-cf12(i1)%ek0*quqd+p7k0*cvqd+p712k0*cvqu
      cauxc=+cf12(i1)%ek0*p7(2)-p7k0*cf12(i1)%e(2)
      l7_12(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_12(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
**QCD                                                                   
       if (ilept(id7).ne.1.and.ilept(id1).ne.1) then
         do iut=1,2
            lg7_12(iut)%a(1) = l7_12(iut)%a(1)
     &           /((fcl(id7))*(fcl(id1)))
            lg7_12(iut)%a(2) = l7_12(iut)%a(2)
     &           /((fcl(id7))*(fcl(id1)))
            lg7_12(iut)%c(1) = l7_12(iut)%c(1)
     &           /((fcl(id7))*(fcl(id1)))
            lg7_12(iut)%c(2) = l7_12(iut)%c(2)
     &           /((fcl(id7))*(fcl(id1)))
         enddo
       endif
  
      ccl=zcl(id7)/(p712q*p712k0)
      do i1=1,2
* TWL0 -- qu=p7,qd=p712,v=cz12(i1)%e,a=l7_12(i1)%a,c=l7_12(i1)%c,cl=ccl,
* nsum=1
      ceps_0=-cz12(i1)%ek0*(p7(2)*p712(3)-p712(2)*p7(3))+p7k0*(c
     & z12(i1)%e(2)*p712(3)-p712(2)*cz12(i1)%e(3))-p712k0*(cz12(
     & i1)%e(2)*p7(3)-p7(2)*cz12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1)%e(3)*p7k0+p7(3)*cz12(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz12(i1)%e(0)*p7(0)-cz12(i1)%e(1)*p7(1)-cz12(i1)%e(2)
     & *p7(2)-cz12(i1)%e(3)*p7(3)
      cvqd=cz12(i1)%e(0)*p712(0)-cz12(i1)%e(1)*p712(1)-cz12(i1)%
     & e(2)*p712(2)-cz12(i1)%e(3)*p712(3)
      cauxa=-cz12(i1)%ek0*quqd+p7k0*cvqd+p712k0*cvqu
      cauxc=+cz12(i1)%ek0*p7(2)-p7k0*cz12(i1)%e(2)
      l7_12(i1)%a(2)=l7_12(i1)%a(2)+ccl*(cauxa-ceps_0)
      l7_12(i1)%c(2)=l7_12(i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      else
* quqd -- p=p7,q=p712
      quqd=p7(0)*p712(0)-p7(1)*p712(1)-p7(2)*p712(2)-p7(3)*p712(
     & 3)
      ccl=zcl(id7)/(p712q*p712k0)
      do i1=1,2
* TWL0 -- qu=p7,qd=p712,v=cz12(i1)%e,a=l7_12(i1)%a,c=l7_12(i1)%c,cl=ccl,
* nsum=0
      ceps_0=-cz12(i1)%ek0*(p7(2)*p712(3)-p712(2)*p7(3))+p7k0*(c
     & z12(i1)%e(2)*p712(3)-p712(2)*cz12(i1)%e(3))-p712k0*(cz12(
     & i1)%e(2)*p7(3)-p7(2)*cz12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1)%e(3)*p7k0+p7(3)*cz12(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz12(i1)%e(0)*p7(0)-cz12(i1)%e(1)*p7(1)-cz12(i1)%e(2)
     & *p7(2)-cz12(i1)%e(3)*p7(3)
      cvqd=cz12(i1)%e(0)*p712(0)-cz12(i1)%e(1)*p712(1)-cz12(i1)%
     & e(2)*p712(2)-cz12(i1)%e(3)*p712(3)
      cauxa=-cz12(i1)%ek0*quqd+p7k0*cvqd+p712k0*cvqu
      cauxc=+cz12(i1)%ek0*p7(2)-p7k0*cz12(i1)%e(2)
      l7_12(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_12(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
  
       do m=0,3
        p734(m) = p7(m)+p34(m)
       enddo
* pk0 -- p=p734
      p734k0=p734(0)-p734(1)
* p.q -- p.q=p734q,p=p734,q=p734,bef=,aft=
      p734q=(p734(0)*p734(0)-p734(1)*p734(1)-p734(2)*p734(2)-p73
     & 4(3)*p734(3))
  
      if (ineutri(id7).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p7,q=p734
      quqd=p7(0)*p734(0)-p7(1)*p734(1)-p7(2)*p734(2)-p7(3)*p734(
     & 3)
      ccl=fcl(id7)/(p734q*p734k0)
      do i1=1,2
* TWL0 -- qu=p7,qd=p734,v=cf34(i1)%e,a=l7_34(i1)%a,c=l7_34(i1)%c,cl=ccl,
* nsum=0
      ceps_0=-cf34(i1)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*(c
     & f34(i1)%e(2)*p734(3)-p734(2)*cf34(i1)%e(3))-p734k0*(cf34(
     & i1)%e(2)*p7(3)-p7(2)*cf34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1)%e(3)*p7k0+p7(3)*cf34(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf34(i1)%e(0)*p7(0)-cf34(i1)%e(1)*p7(1)-cf34(i1)%e(2)
     & *p7(2)-cf34(i1)%e(3)*p7(3)
      cvqd=cf34(i1)%e(0)*p734(0)-cf34(i1)%e(1)*p734(1)-cf34(i1)%
     & e(2)*p734(2)-cf34(i1)%e(3)*p734(3)
      cauxa=-cf34(i1)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cf34(i1)%ek0*p7(2)-p7k0*cf34(i1)%e(2)
      l7_34(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_34(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
  
**QCD                                                                   
       if (ilept(id7).ne.1.and.ilept(id3).ne.1) then
         do iut=1,2
            lg7_34(iut)%a(1) = l7_34(iut)%a(1)
     &           /((fcl(id7))*(fcl(id3)))
            lg7_34(iut)%a(2) = l7_34(iut)%a(2)
     &           /((fcl(id7))*(fcl(id3)))
            lg7_34(iut)%c(1) = l7_34(iut)%c(1)
     &           /((fcl(id7))*(fcl(id3)))
            lg7_34(iut)%c(2) = l7_34(iut)%c(2)
     &           /((fcl(id7))*(fcl(id3)))
         enddo
       endif
  
      ccl=zcl(id7)/(p734q*p734k0)
      do i1=1,2
* TWL0 -- qu=p7,qd=p734,v=cz34(i1)%e,a=l7_34(i1)%a,c=l7_34(i1)%c,cl=ccl,
* nsum=1
      ceps_0=-cz34(i1)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*(c
     & z34(i1)%e(2)*p734(3)-p734(2)*cz34(i1)%e(3))-p734k0*(cz34(
     & i1)%e(2)*p7(3)-p7(2)*cz34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1)%e(3)*p7k0+p7(3)*cz34(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i1)%e(0)*p7(0)-cz34(i1)%e(1)*p7(1)-cz34(i1)%e(2)
     & *p7(2)-cz34(i1)%e(3)*p7(3)
      cvqd=cz34(i1)%e(0)*p734(0)-cz34(i1)%e(1)*p734(1)-cz34(i1)%
     & e(2)*p734(2)-cz34(i1)%e(3)*p734(3)
      cauxa=-cz34(i1)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cz34(i1)%ek0*p7(2)-p7k0*cz34(i1)%e(2)
      l7_34(i1)%a(2)=l7_34(i1)%a(2)+ccl*(cauxa-ceps_0)
      l7_34(i1)%c(2)=l7_34(i1)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      else
* quqd -- p=p7,q=p734
      quqd=p7(0)*p734(0)-p7(1)*p734(1)-p7(2)*p734(2)-p7(3)*p734(
     & 3)
      ccl=zcl(id7)/(p734q*p734k0)
      do i1=1,2
* TWL0 -- qu=p7,qd=p734,v=cz34(i1)%e,a=l7_34(i1)%a,c=l7_34(i1)%c,cl=ccl,
* nsum=0
      ceps_0=-cz34(i1)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*(c
     & z34(i1)%e(2)*p734(3)-p734(2)*cz34(i1)%e(3))-p734k0*(cz34(
     & i1)%e(2)*p7(3)-p7(2)*cz34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1)%e(3)*p7k0+p7(3)*cz34(i1)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i1)%e(0)*p7(0)-cz34(i1)%e(1)*p7(1)-cz34(i1)%e(2)
     & *p7(2)-cz34(i1)%e(3)*p7(3)
      cvqd=cz34(i1)%e(0)*p734(0)-cz34(i1)%e(1)*p734(1)-cz34(i1)%
     & e(2)*p734(2)-cz34(i1)%e(3)*p734(3)
      cauxa=-cz34(i1)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cz34(i1)%ek0*p7(2)-p7k0*cz34(i1)%e(2)
      l7_34(i1)%a(2)=ccl*(cauxa-ceps_0)
      l7_34(i1)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      endif
  
  
*     W insertion                                                       
       do m=0,3
        p578(m) = p5(m)+p78(m)
       enddo
* pk0 -- p=p578
      p578k0=p578(0)-p578(1)
* p.q -- p.q=p578q,p=p578,q=p578,bef=,aft=
      p578q=(p578(0)*p578(0)-p578(1)*p578(1)-p578(2)*p578(2)-p57
     & 8(3)*p578(3))
  
* quqd -- p=p5,q=p578
      quqd=p5(0)*p578(0)-p5(1)*p578(1)-p5(2)*p578(2)-p5(3)*p578(
     & 3)
      ccl=wcl/(p578q*p578k0)
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
  
       do m=0,3
        p756(m) = p7(m)+p56(m)
       enddo
* pk0 -- p=p756
      p756k0=p756(0)-p756(1)
* p.q -- p.q=p756q,p=p756,q=p756,bef=,aft=
      p756q=(p756(0)*p756(0)-p756(1)*p756(1)-p756(2)*p756(2)-p75
     & 6(3)*p756(3))
  
* quqd -- p=p7,q=p756
      quqd=p7(0)*p756(0)-p7(1)*p756(1)-p7(2)*p756(2)-p7(3)*p756(
     & 3)
      ccl=wcl/(p756q*p756k0)
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
  
  
*                                                                       
* RIGHT INSERTIONS, with fermion propagator and pk0                     
*                                                                       
  
* main line Z (12 or 34)                                                
  
*     W insertion                                                       
  
      if (iup(id2).eq.1) then !W = W+ (56)
       do m=0,3
        p256(m) = -p2(m)-p56(m)
       enddo
* pk0 -- p=p256
      p256k0=p256(0)-p256(1)
* p.q -- p.q=p256q,p=p256,q=p256,bef=,aft=
      p256q=(p256(0)*p256(0)-p256(1)*p256(1)-p256(2)*p256(2)-p25
     & 6(3)*p256(3))
  
* quqd -- p=p256,q=p2
      quqd=p256(0)*p2(0)-p256(1)*p2(1)-p256(2)*p2(2)-p256(3)*p2(
     & 3)
      ccl=wcl/(p256q*p256k0)
* TWR0 -- qu=p256,qd=p2,v=cw56%e,a=r2_56%a,b=r2_56%b,cl=ccl,nsum=0
      ceps_0=-cw56%ek0*(p256(2)*p2(3)-p2(2)*p256(3))+p256k0*(cw5
     & 6%e(2)*p2(3)-p2(2)*cw56%e(3))-p2k0*(cw56%e(2)*p256(3)-p25
     & 6(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw56%e(3)*p2k0+p2(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p256(0)-cw56%e(1)*p256(1)-cw56%e(2)*p256(2)
     & -cw56%e(3)*p256(3)
      cvqd=cw56%e(0)*p2(0)-cw56%e(1)*p2(1)-cw56%e(2)*p2(2)-cw56%
     & e(3)*p2(3)
      cauxa=-cw56%ek0*quqd+p256k0*cvqd+p2k0*cvqu
      cauxb=-cw56%ek0*p2(2)+p2k0*cw56%e(2)
      r2_56%a(2)=ccl*(cauxa-ceps_0)
      r2_56%b(1)=ccl*(cauxb-ceps_2)
  
      else !W = W- (78)
       do m=0,3
        p278(m) = -p2(m)-p78(m)
       enddo
* pk0 -- p=p278
      p278k0=p278(0)-p278(1)
* p.q -- p.q=p278q,p=p278,q=p278,bef=,aft=
      p278q=(p278(0)*p278(0)-p278(1)*p278(1)-p278(2)*p278(2)-p27
     & 8(3)*p278(3))
  
* quqd -- p=p278,q=p2
      quqd=p278(0)*p2(0)-p278(1)*p2(1)-p278(2)*p2(2)-p278(3)*p2(
     & 3)
      ccl=wcl/(p278q*p278k0)
* TWR0 -- qu=p278,qd=p2,v=cw78%e,a=r2_78%a,b=r2_78%b,cl=ccl,nsum=0
      ceps_0=-cw78%ek0*(p278(2)*p2(3)-p2(2)*p278(3))+p278k0*(cw7
     & 8%e(2)*p2(3)-p2(2)*cw78%e(3))-p2k0*(cw78%e(2)*p278(3)-p27
     & 8(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw78%e(3)*p2k0+p2(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p278(0)-cw78%e(1)*p278(1)-cw78%e(2)*p278(2)
     & -cw78%e(3)*p278(3)
      cvqd=cw78%e(0)*p2(0)-cw78%e(1)*p2(1)-cw78%e(2)*p2(2)-cw78%
     & e(3)*p2(3)
      cauxa=-cw78%ek0*quqd+p278k0*cvqd+p2k0*cvqu
      cauxb=-cw78%ek0*p2(2)+p2k0*cw78%e(2)
      r2_78%a(2)=ccl*(cauxa-ceps_0)
      r2_78%b(1)=ccl*(cauxb-ceps_2)
  
      endif
  
      if (iup(id4).eq.1) then !W = W+ (56)
       do m=0,3
        p456(m) = -p4(m)-p56(m)
       enddo
* pk0 -- p=p456
      p456k0=p456(0)-p456(1)
* p.q -- p.q=p456q,p=p456,q=p456,bef=,aft=
      p456q=(p456(0)*p456(0)-p456(1)*p456(1)-p456(2)*p456(2)-p45
     & 6(3)*p456(3))
  
* quqd -- p=p456,q=p4
      quqd=p456(0)*p4(0)-p456(1)*p4(1)-p456(2)*p4(2)-p456(3)*p4(
     & 3)
      ccl=wcl/(p456q*p456k0)
* TWR0 -- qu=p456,qd=p4,v=cw56%e,a=r4_56%a,b=r4_56%b,cl=ccl,nsum=0
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
      r4_56%a(2)=ccl*(cauxa-ceps_0)
      r4_56%b(1)=ccl*(cauxb-ceps_2)
  
      else !W = W- (78)
       do m=0,3
        p478(m) = -p4(m)-p78(m)
       enddo
* pk0 -- p=p478
      p478k0=p478(0)-p478(1)
* p.q -- p.q=p478q,p=p478,q=p478,bef=,aft=
      p478q=(p478(0)*p478(0)-p478(1)*p478(1)-p478(2)*p478(2)-p47
     & 8(3)*p478(3))
  
* quqd -- p=p478,q=p4
      quqd=p478(0)*p4(0)-p478(1)*p4(1)-p478(2)*p4(2)-p478(3)*p4(
     & 3)
      ccl=wcl/(p478q*p478k0)
* TWR0 -- qu=p478,qd=p4,v=cw78%e,a=r4_78%a,b=r4_78%b,cl=ccl,nsum=0
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
      r4_78%a(2)=ccl*(cauxa-ceps_0)
      r4_78%b(1)=ccl*(cauxb-ceps_2)
  
      endif
  
*     Z insertion                                                       
       do m=0,3
        p234(m) = -p2(m)-p34(m)
       enddo
* pk0 -- p=p234
      p234k0=p234(0)-p234(1)
* p.q -- p.q=p234q,p=p234,q=p234,bef=,aft=
      p234q=(p234(0)*p234(0)-p234(1)*p234(1)-p234(2)*p234(2)-p23
     & 4(3)*p234(3))
  
      if (ineutri(id2).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p234,q=p2
      quqd=p234(0)*p2(0)-p234(1)*p2(1)-p234(2)*p2(2)-p234(3)*p2(
     & 3)
      ccr=fcr(id2)/(p234q*p234k0)
      ccl=fcl(id2)/(p234q*p234k0)
      do i3=1,2
* TR0 -- qu=p234,qd=p2,v=cf34(i3)%e,a=r2_34(i3)%a,b=r2_34(i3)%b,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cf34(i3)%ek0*(p234(2)*p2(3)-p2(2)*p234(3))+p234k0*
     & (cf34(i3)%e(2)*p2(3)-p2(2)*cf34(i3)%e(3))-p2k0*(cf34(i3)%
     & e(2)*p234(3)-p234(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf34(i3)%e(3)*p2k0+p2(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p234(0)-cf34(i3)%e(1)*p234(1)-cf34(i3)%
     & e(2)*p234(2)-cf34(i3)%e(3)*p234(3)
      cvqd=cf34(i3)%e(0)*p2(0)-cf34(i3)%e(1)*p2(1)-cf34(i3)%e(2)
     & *p2(2)-cf34(i3)%e(3)*p2(3)
      cauxa=-cf34(i3)%ek0*quqd+p234k0*cvqd+p2k0*cvqu
      cauxb=-cf34(i3)%ek0*p2(2)+p2k0*cf34(i3)%e(2)
      r2_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      r2_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      r2_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      r2_34(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
**QCD                                                                   
       if (ilept(id2).ne.1.and.ilept(id3).ne.1) then
         do iut=1,2
            rg2_34(iut)%a(1) = r2_34(iut)%a(1)
     &           /((fcl(id2))*(fcl(id3)))
            rg2_34(iut)%a(2) = r2_34(iut)%a(2)
     &           /((fcl(id2))*(fcl(id3)))
            rg2_34(iut)%b(1) = r2_34(iut)%b(1)
     &           /((fcl(id2))*(fcl(id3)))
            rg2_34(iut)%b(2) = r2_34(iut)%b(2)
     &           /((fcl(id2))*(fcl(id3)))
         enddo
       endif
      ccr=zcr(id2)/(p234q*p234k0)
      ccl=zcl(id2)/(p234q*p234k0)
      do i3=1,2
* TR0 -- qu=p234,qd=p2,v=cz34(i3)%e,a=r2_34(i3)%a,b=r2_34(i3)%b,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cz34(i3)%ek0*(p234(2)*p2(3)-p2(2)*p234(3))+p234k0*
     & (cz34(i3)%e(2)*p2(3)-p2(2)*cz34(i3)%e(3))-p2k0*(cz34(i3)%
     & e(2)*p234(3)-p234(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i3)%e(3)*p2k0+p2(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p234(0)-cz34(i3)%e(1)*p234(1)-cz34(i3)%
     & e(2)*p234(2)-cz34(i3)%e(3)*p234(3)
      cvqd=cz34(i3)%e(0)*p2(0)-cz34(i3)%e(1)*p2(1)-cz34(i3)%e(2)
     & *p2(2)-cz34(i3)%e(3)*p2(3)
      cauxa=-cz34(i3)%ek0*quqd+p234k0*cvqd+p2k0*cvqu
      cauxb=-cz34(i3)%ek0*p2(2)+p2k0*cz34(i3)%e(2)
      r2_34(i3)%a(1)=r2_34(i3)%a(1)+ccr*(cauxa+ceps_0)
      r2_34(i3)%a(2)=r2_34(i3)%a(2)+ccl*(cauxa-ceps_0)
      r2_34(i3)%b(1)=r2_34(i3)%b(1)+ccl*(cauxb-ceps_2)
      r2_34(i3)%b(2)=r2_34(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      else
* quqd -- p=p234,q=p2
      quqd=p234(0)*p2(0)-p234(1)*p2(1)-p234(2)*p2(2)-p234(3)*p2(
     & 3)
      ccr=zcr(id2)/(p234q*p234k0)
      ccl=zcl(id2)/(p234q*p234k0)
      do i3=1,2
* TR0 -- qu=p234,qd=p2,v=cz34(i3)%e,a=r2_34(i3)%a,b=r2_34(i3)%b,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz34(i3)%ek0*(p234(2)*p2(3)-p2(2)*p234(3))+p234k0*
     & (cz34(i3)%e(2)*p2(3)-p2(2)*cz34(i3)%e(3))-p2k0*(cz34(i3)%
     & e(2)*p234(3)-p234(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i3)%e(3)*p2k0+p2(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p234(0)-cz34(i3)%e(1)*p234(1)-cz34(i3)%
     & e(2)*p234(2)-cz34(i3)%e(3)*p234(3)
      cvqd=cz34(i3)%e(0)*p2(0)-cz34(i3)%e(1)*p2(1)-cz34(i3)%e(2)
     & *p2(2)-cz34(i3)%e(3)*p2(3)
      cauxa=-cz34(i3)%ek0*quqd+p234k0*cvqd+p2k0*cvqu
      cauxb=-cz34(i3)%ek0*p2(2)+p2k0*cz34(i3)%e(2)
      r2_34(i3)%a(1)=ccr*(cauxa+ceps_0)
      r2_34(i3)%a(2)=ccl*(cauxa-ceps_0)
      r2_34(i3)%b(1)=ccl*(cauxb-ceps_2)
      r2_34(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
      endif
  
       do m=0,3
        p412(m) = -p4(m)-p12(m)
       enddo
* pk0 -- p=p412
      p412k0=p412(0)-p412(1)
* p.q -- p.q=p412q,p=p412,q=p412,bef=,aft=
      p412q=(p412(0)*p412(0)-p412(1)*p412(1)-p412(2)*p412(2)-p41
     & 2(3)*p412(3))
  
      if (ineutri(id4).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p412,q=p4
      quqd=p412(0)*p4(0)-p412(1)*p4(1)-p412(2)*p4(2)-p412(3)*p4(
     & 3)
      ccr=fcr(id4)/(p412q*p412k0)
      ccl=fcl(id4)/(p412q*p412k0)
      do i3=1,2
* TR0 -- qu=p412,qd=p4,v=cf12(i3)%e,a=r4_12(i3)%a,b=r4_12(i3)%b,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cf12(i3)%ek0*(p412(2)*p4(3)-p4(2)*p412(3))+p412k0*
     & (cf12(i3)%e(2)*p4(3)-p4(2)*cf12(i3)%e(3))-p4k0*(cf12(i3)%
     & e(2)*p412(3)-p412(2)*cf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf12(i3)%e(3)*p4k0+p4(3)*cf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i3)%e(0)*p412(0)-cf12(i3)%e(1)*p412(1)-cf12(i3)%
     & e(2)*p412(2)-cf12(i3)%e(3)*p412(3)
      cvqd=cf12(i3)%e(0)*p4(0)-cf12(i3)%e(1)*p4(1)-cf12(i3)%e(2)
     & *p4(2)-cf12(i3)%e(3)*p4(3)
      cauxa=-cf12(i3)%ek0*quqd+p412k0*cvqd+p4k0*cvqu
      cauxb=-cf12(i3)%ek0*p4(2)+p4k0*cf12(i3)%e(2)
      r4_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      r4_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      r4_12(i3)%b(1)=ccl*(cauxb-ceps_2)
      r4_12(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
  
**QCD                                                                   
       if (ilept(id2).ne.1.and.ilept(id3).ne.1) then
         do iut=1,2
            rg4_12(iut)%a(1) = r4_12(iut)%a(1)
     &           /((fcl(id2))*(fcl(id3)))
            rg4_12(iut)%a(2) = r4_12(iut)%a(2)
     &           /((fcl(id2))*(fcl(id3)))
            rg4_12(iut)%b(1) = r4_12(iut)%b(1)
     &           /((fcl(id2))*(fcl(id3)))
            rg4_12(iut)%b(2) = r4_12(iut)%b(2)
     &           /((fcl(id2))*(fcl(id3)))
         enddo
       endif
      ccr=zcr(id4)/(p412q*p412k0)
      ccl=zcl(id4)/(p412q*p412k0)
      do i3=1,2
* TR0 -- qu=p412,qd=p4,v=cz12(i3)%e,a=r4_12(i3)%a,b=r4_12(i3)%b,cr=ccr,c
* l=ccl,nsum=1
      ceps_0=-cz12(i3)%ek0*(p412(2)*p4(3)-p4(2)*p412(3))+p412k0*
     & (cz12(i3)%e(2)*p4(3)-p4(2)*cz12(i3)%e(3))-p4k0*(cz12(i3)%
     & e(2)*p412(3)-p412(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz12(i3)%e(3)*p4k0+p4(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p412(0)-cz12(i3)%e(1)*p412(1)-cz12(i3)%
     & e(2)*p412(2)-cz12(i3)%e(3)*p412(3)
      cvqd=cz12(i3)%e(0)*p4(0)-cz12(i3)%e(1)*p4(1)-cz12(i3)%e(2)
     & *p4(2)-cz12(i3)%e(3)*p4(3)
      cauxa=-cz12(i3)%ek0*quqd+p412k0*cvqd+p4k0*cvqu
      cauxb=-cz12(i3)%ek0*p4(2)+p4k0*cz12(i3)%e(2)
      r4_12(i3)%a(1)=r4_12(i3)%a(1)+ccr*(cauxa+ceps_0)
      r4_12(i3)%a(2)=r4_12(i3)%a(2)+ccl*(cauxa-ceps_0)
      r4_12(i3)%b(1)=r4_12(i3)%b(1)+ccl*(cauxb-ceps_2)
      r4_12(i3)%b(2)=r4_12(i3)%b(2)+ccr*(-cauxb-ceps_2)
      end do
      else
* quqd -- p=p412,q=p4
      quqd=p412(0)*p4(0)-p412(1)*p4(1)-p412(2)*p4(2)-p412(3)*p4(
     & 3)
      ccr=zcr(id4)/(p412q*p412k0)
      ccl=zcl(id4)/(p412q*p412k0)
      do i3=1,2
* TR0 -- qu=p412,qd=p4,v=cz12(i3)%e,a=r4_12(i3)%a,b=r4_12(i3)%b,cr=ccr,c
* l=ccl,nsum=0
      ceps_0=-cz12(i3)%ek0*(p412(2)*p4(3)-p4(2)*p412(3))+p412k0*
     & (cz12(i3)%e(2)*p4(3)-p4(2)*cz12(i3)%e(3))-p4k0*(cz12(i3)%
     & e(2)*p412(3)-p412(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz12(i3)%e(3)*p4k0+p4(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p412(0)-cz12(i3)%e(1)*p412(1)-cz12(i3)%
     & e(2)*p412(2)-cz12(i3)%e(3)*p412(3)
      cvqd=cz12(i3)%e(0)*p4(0)-cz12(i3)%e(1)*p4(1)-cz12(i3)%e(2)
     & *p4(2)-cz12(i3)%e(3)*p4(3)
      cauxa=-cz12(i3)%ek0*quqd+p412k0*cvqd+p4k0*cvqu
      cauxb=-cz12(i3)%ek0*p4(2)+p4k0*cz12(i3)%e(2)
      r4_12(i3)%a(1)=ccr*(cauxa+ceps_0)
      r4_12(i3)%a(2)=ccl*(cauxa-ceps_0)
      r4_12(i3)%b(1)=ccl*(cauxb-ceps_2)
      r4_12(i3)%b(2)=ccr*(-cauxb-ceps_2)
      end do
      endif
  
  
* main line W (56 or 78)                                                
  
*     Z insertion                                                       
       do m=0,3
        p612(m) = -p6(m)-p12(m)
       enddo
* pk0 -- p=p612
      p612k0=p612(0)-p612(1)
* p.q -- p.q=p612q,p=p612,q=p612,bef=,aft=
      p612q=(p612(0)*p612(0)-p612(1)*p612(1)-p612(2)*p612(2)-p61
     & 2(3)*p612(3))
  
      if (ineutri(id6).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p612,q=p6
      quqd=p612(0)*p6(0)-p612(1)*p6(1)-p612(2)*p6(2)-p612(3)*p6(
     & 3)
      ccl=fcl(id6)/(p612q*p612k0)
      do i1=1,2
* TWR0 -- qu=p612,qd=p6,v=cf12(i1)%e,a=r6_12(i1)%a,b=r6_12(i1)%b,cl=ccl,
* nsum=0
      ceps_0=-cf12(i1)%ek0*(p612(2)*p6(3)-p6(2)*p612(3))+p612k0*
     & (cf12(i1)%e(2)*p6(3)-p6(2)*cf12(i1)%e(3))-p6k0*(cf12(i1)%
     & e(2)*p612(3)-p612(2)*cf12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf12(i1)%e(3)*p6k0+p6(3)*cf12(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1)%e(0)*p612(0)-cf12(i1)%e(1)*p612(1)-cf12(i1)%
     & e(2)*p612(2)-cf12(i1)%e(3)*p612(3)
      cvqd=cf12(i1)%e(0)*p6(0)-cf12(i1)%e(1)*p6(1)-cf12(i1)%e(2)
     & *p6(2)-cf12(i1)%e(3)*p6(3)
      cauxa=-cf12(i1)%ek0*quqd+p612k0*cvqd+p6k0*cvqu
      cauxb=-cf12(i1)%ek0*p6(2)+p6k0*cf12(i1)%e(2)
      r6_12(i1)%a(2)=ccl*(cauxa-ceps_0)
      r6_12(i1)%b(1)=ccl*(cauxb-ceps_2)
      end do
  
**QCD                                                                   
       if (ilept(id6).ne.1.and.ilept(id1).ne.1) then
         do iut=1,2
            rg6_12(iut)%a(1) = r6_12(iut)%a(1)
     &           /((fcl(id6))*(fcl(id1)))
            rg6_12(iut)%a(2) = r6_12(iut)%a(2)
     &           /((fcl(id6))*(fcl(id1)))
            rg6_12(iut)%b(1) = r6_12(iut)%b(1)
     &           /((fcl(id6))*(fcl(id1)))
            rg6_12(iut)%b(2) = r6_12(iut)%b(2)
     &           /((fcl(id6))*(fcl(id1)))
         enddo
       endif
      ccl=zcl(id6)/(p612q*p612k0)
      do i1=1,2
* TWR0 -- qu=p612,qd=p6,v=cz12(i1)%e,a=r6_12(i1)%a,b=r6_12(i1)%b,cl=ccl,
* nsum=1
      ceps_0=-cz12(i1)%ek0*(p612(2)*p6(3)-p6(2)*p612(3))+p612k0*
     & (cz12(i1)%e(2)*p6(3)-p6(2)*cz12(i1)%e(3))-p6k0*(cz12(i1)%
     & e(2)*p612(3)-p612(2)*cz12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz12(i1)%e(3)*p6k0+p6(3)*cz12(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1)%e(0)*p612(0)-cz12(i1)%e(1)*p612(1)-cz12(i1)%
     & e(2)*p612(2)-cz12(i1)%e(3)*p612(3)
      cvqd=cz12(i1)%e(0)*p6(0)-cz12(i1)%e(1)*p6(1)-cz12(i1)%e(2)
     & *p6(2)-cz12(i1)%e(3)*p6(3)
      cauxa=-cz12(i1)%ek0*quqd+p612k0*cvqd+p6k0*cvqu
      cauxb=-cz12(i1)%ek0*p6(2)+p6k0*cz12(i1)%e(2)
      r6_12(i1)%a(2)=r6_12(i1)%a(2)+ccl*(cauxa-ceps_0)
      r6_12(i1)%b(1)=r6_12(i1)%b(1)+ccl*(cauxb-ceps_2)
      end do
      else
* quqd -- p=p612,q=p6
      quqd=p612(0)*p6(0)-p612(1)*p6(1)-p612(2)*p6(2)-p612(3)*p6(
     & 3)
      ccl=zcl(id6)/(p612q*p612k0)
      do i1=1,2
* TWR0 -- qu=p612,qd=p6,v=cz12(i1)%e,a=r6_12(i1)%a,b=r6_12(i1)%b,cl=ccl,
* nsum=0
      ceps_0=-cz12(i1)%ek0*(p612(2)*p6(3)-p6(2)*p612(3))+p612k0*
     & (cz12(i1)%e(2)*p6(3)-p6(2)*cz12(i1)%e(3))-p6k0*(cz12(i1)%
     & e(2)*p612(3)-p612(2)*cz12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz12(i1)%e(3)*p6k0+p6(3)*cz12(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1)%e(0)*p612(0)-cz12(i1)%e(1)*p612(1)-cz12(i1)%
     & e(2)*p612(2)-cz12(i1)%e(3)*p612(3)
      cvqd=cz12(i1)%e(0)*p6(0)-cz12(i1)%e(1)*p6(1)-cz12(i1)%e(2)
     & *p6(2)-cz12(i1)%e(3)*p6(3)
      cauxa=-cz12(i1)%ek0*quqd+p612k0*cvqd+p6k0*cvqu
      cauxb=-cz12(i1)%ek0*p6(2)+p6k0*cz12(i1)%e(2)
      r6_12(i1)%a(2)=ccl*(cauxa-ceps_0)
      r6_12(i1)%b(1)=ccl*(cauxb-ceps_2)
      end do
      endif
  
       do m=0,3
        p634(m) = -p6(m)-p34(m)
       enddo
* pk0 -- p=p634
      p634k0=p634(0)-p634(1)
* p.q -- p.q=p634q,p=p634,q=p634,bef=,aft=
      p634q=(p634(0)*p634(0)-p634(1)*p634(1)-p634(2)*p634(2)-p63
     & 4(3)*p634(3))
  
      if (ineutri(id6).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p634,q=p6
      quqd=p634(0)*p6(0)-p634(1)*p6(1)-p634(2)*p6(2)-p634(3)*p6(
     & 3)
      ccl=fcl(id6)/(p634q*p634k0)
      do i1=1,2
* TWR0 -- qu=p634,qd=p6,v=cf34(i1)%e,a=r6_34(i1)%a,b=r6_34(i1)%b,cl=ccl,
* nsum=0
      ceps_0=-cf34(i1)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k0*
     & (cf34(i1)%e(2)*p6(3)-p6(2)*cf34(i1)%e(3))-p6k0*(cf34(i1)%
     & e(2)*p634(3)-p634(2)*cf34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf34(i1)%e(3)*p6k0+p6(3)*cf34(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1)%e(0)*p634(0)-cf34(i1)%e(1)*p634(1)-cf34(i1)%
     & e(2)*p634(2)-cf34(i1)%e(3)*p634(3)
      cvqd=cf34(i1)%e(0)*p6(0)-cf34(i1)%e(1)*p6(1)-cf34(i1)%e(2)
     & *p6(2)-cf34(i1)%e(3)*p6(3)
      cauxa=-cf34(i1)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cf34(i1)%ek0*p6(2)+p6k0*cf34(i1)%e(2)
      r6_34(i1)%a(2)=ccl*(cauxa-ceps_0)
      r6_34(i1)%b(1)=ccl*(cauxb-ceps_2)
      end do
  
**QCD                                                                   
       if (ilept(id6).ne.1.and.ilept(id3).ne.1) then
         do iut=1,2
            rg6_34(iut)%a(1) = r6_34(iut)%a(1)
     &           /((fcl(id6))*(fcl(id3)))
            rg6_34(iut)%a(2) = r6_34(iut)%a(2)
     &           /((fcl(id6))*(fcl(id3)))
            rg6_34(iut)%b(1) = r6_34(iut)%b(1)
     &           /((fcl(id6))*(fcl(id3)))
            rg6_34(iut)%b(2) = r6_34(iut)%b(2)
     &           /((fcl(id6))*(fcl(id3)))
         enddo
       endif
      ccl=zcl(id6)/(p634q*p634k0)
      do i1=1,2
* TWR0 -- qu=p634,qd=p6,v=cz34(i1)%e,a=r6_34(i1)%a,b=r6_34(i1)%b,cl=ccl,
* nsum=1
      ceps_0=-cz34(i1)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k0*
     & (cz34(i1)%e(2)*p6(3)-p6(2)*cz34(i1)%e(3))-p6k0*(cz34(i1)%
     & e(2)*p634(3)-p634(2)*cz34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i1)%e(3)*p6k0+p6(3)*cz34(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1)%e(0)*p634(0)-cz34(i1)%e(1)*p634(1)-cz34(i1)%
     & e(2)*p634(2)-cz34(i1)%e(3)*p634(3)
      cvqd=cz34(i1)%e(0)*p6(0)-cz34(i1)%e(1)*p6(1)-cz34(i1)%e(2)
     & *p6(2)-cz34(i1)%e(3)*p6(3)
      cauxa=-cz34(i1)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cz34(i1)%ek0*p6(2)+p6k0*cz34(i1)%e(2)
      r6_34(i1)%a(2)=r6_34(i1)%a(2)+ccl*(cauxa-ceps_0)
      r6_34(i1)%b(1)=r6_34(i1)%b(1)+ccl*(cauxb-ceps_2)
      end do
      else
* quqd -- p=p634,q=p6
      quqd=p634(0)*p6(0)-p634(1)*p6(1)-p634(2)*p6(2)-p634(3)*p6(
     & 3)
      ccl=zcl(id6)/(p634q*p634k0)
      do i1=1,2
* TWR0 -- qu=p634,qd=p6,v=cz34(i1)%e,a=r6_34(i1)%a,b=r6_34(i1)%b,cl=ccl,
* nsum=0
      ceps_0=-cz34(i1)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k0*
     & (cz34(i1)%e(2)*p6(3)-p6(2)*cz34(i1)%e(3))-p6k0*(cz34(i1)%
     & e(2)*p634(3)-p634(2)*cz34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i1)%e(3)*p6k0+p6(3)*cz34(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1)%e(0)*p634(0)-cz34(i1)%e(1)*p634(1)-cz34(i1)%
     & e(2)*p634(2)-cz34(i1)%e(3)*p634(3)
      cvqd=cz34(i1)%e(0)*p6(0)-cz34(i1)%e(1)*p6(1)-cz34(i1)%e(2)
     & *p6(2)-cz34(i1)%e(3)*p6(3)
      cauxa=-cz34(i1)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cz34(i1)%ek0*p6(2)+p6k0*cz34(i1)%e(2)
      r6_34(i1)%a(2)=ccl*(cauxa-ceps_0)
      r6_34(i1)%b(1)=ccl*(cauxb-ceps_2)
      end do
      endif
  
       do m=0,3
        p812(m) = -p8(m)-p12(m)
       enddo
* pk0 -- p=p812
      p812k0=p812(0)-p812(1)
* p.q -- p.q=p812q,p=p812,q=p812,bef=,aft=
      p812q=(p812(0)*p812(0)-p812(1)*p812(1)-p812(2)*p812(2)-p81
     & 2(3)*p812(3))
  
      if (ineutri(id8).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p812,q=p8
      quqd=p812(0)*p8(0)-p812(1)*p8(1)-p812(2)*p8(2)-p812(3)*p8(
     & 3)
      ccl=fcl(id8)/(p812q*p812k0)
      do i1=1,2
* TWR0 -- qu=p812,qd=p8,v=cf12(i1)%e,a=r8_12(i1)%a,b=r8_12(i1)%b,cl=ccl,
* nsum=0
      ceps_0=-cf12(i1)%ek0*(p812(2)*p8(3)-p8(2)*p812(3))+p812k0*
     & (cf12(i1)%e(2)*p8(3)-p8(2)*cf12(i1)%e(3))-p8k0*(cf12(i1)%
     & e(2)*p812(3)-p812(2)*cf12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf12(i1)%e(3)*p8k0+p8(3)*cf12(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1)%e(0)*p812(0)-cf12(i1)%e(1)*p812(1)-cf12(i1)%
     & e(2)*p812(2)-cf12(i1)%e(3)*p812(3)
      cvqd=cf12(i1)%e(0)*p8(0)-cf12(i1)%e(1)*p8(1)-cf12(i1)%e(2)
     & *p8(2)-cf12(i1)%e(3)*p8(3)
      cauxa=-cf12(i1)%ek0*quqd+p812k0*cvqd+p8k0*cvqu
      cauxb=-cf12(i1)%ek0*p8(2)+p8k0*cf12(i1)%e(2)
      r8_12(i1)%a(2)=ccl*(cauxa-ceps_0)
      r8_12(i1)%b(1)=ccl*(cauxb-ceps_2)
      end do
  
**QCD                                                                   
       if (ilept(id8).ne.1.and.ilept(id1).ne.1) then
         do iut=1,2
            rg8_12(iut)%a(1) = r8_12(iut)%a(1)
     &           /((fcl(id8))*(fcl(id1)))
            rg8_12(iut)%a(2) = r8_12(iut)%a(2)
     &           /((fcl(id8))*(fcl(id1)))
            rg8_12(iut)%b(1) = r8_12(iut)%b(1)
     &           /((fcl(id8))*(fcl(id1)))
            rg8_12(iut)%b(2) = r8_12(iut)%b(2)
     &           /((fcl(id8))*(fcl(id1)))
         enddo
       endif
      ccl=zcl(id8)/(p812q*p812k0)
      do i1=1,2
* TWR0 -- qu=p812,qd=p8,v=cz12(i1)%e,a=r8_12(i1)%a,b=r8_12(i1)%b,cl=ccl,
* nsum=1
      ceps_0=-cz12(i1)%ek0*(p812(2)*p8(3)-p8(2)*p812(3))+p812k0*
     & (cz12(i1)%e(2)*p8(3)-p8(2)*cz12(i1)%e(3))-p8k0*(cz12(i1)%
     & e(2)*p812(3)-p812(2)*cz12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz12(i1)%e(3)*p8k0+p8(3)*cz12(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1)%e(0)*p812(0)-cz12(i1)%e(1)*p812(1)-cz12(i1)%
     & e(2)*p812(2)-cz12(i1)%e(3)*p812(3)
      cvqd=cz12(i1)%e(0)*p8(0)-cz12(i1)%e(1)*p8(1)-cz12(i1)%e(2)
     & *p8(2)-cz12(i1)%e(3)*p8(3)
      cauxa=-cz12(i1)%ek0*quqd+p812k0*cvqd+p8k0*cvqu
      cauxb=-cz12(i1)%ek0*p8(2)+p8k0*cz12(i1)%e(2)
      r8_12(i1)%a(2)=r8_12(i1)%a(2)+ccl*(cauxa-ceps_0)
      r8_12(i1)%b(1)=r8_12(i1)%b(1)+ccl*(cauxb-ceps_2)
      end do
      else
* quqd -- p=p812,q=p8
      quqd=p812(0)*p8(0)-p812(1)*p8(1)-p812(2)*p8(2)-p812(3)*p8(
     & 3)
      ccl=zcl(id8)/(p812q*p812k0)
      do i1=1,2
* TWR0 -- qu=p812,qd=p8,v=cz12(i1)%e,a=r8_12(i1)%a,b=r8_12(i1)%b,cl=ccl,
* nsum=0
      ceps_0=-cz12(i1)%ek0*(p812(2)*p8(3)-p8(2)*p812(3))+p812k0*
     & (cz12(i1)%e(2)*p8(3)-p8(2)*cz12(i1)%e(3))-p8k0*(cz12(i1)%
     & e(2)*p812(3)-p812(2)*cz12(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz12(i1)%e(3)*p8k0+p8(3)*cz12(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1)%e(0)*p812(0)-cz12(i1)%e(1)*p812(1)-cz12(i1)%
     & e(2)*p812(2)-cz12(i1)%e(3)*p812(3)
      cvqd=cz12(i1)%e(0)*p8(0)-cz12(i1)%e(1)*p8(1)-cz12(i1)%e(2)
     & *p8(2)-cz12(i1)%e(3)*p8(3)
      cauxa=-cz12(i1)%ek0*quqd+p812k0*cvqd+p8k0*cvqu
      cauxb=-cz12(i1)%ek0*p8(2)+p8k0*cz12(i1)%e(2)
      r8_12(i1)%a(2)=ccl*(cauxa-ceps_0)
      r8_12(i1)%b(1)=ccl*(cauxb-ceps_2)
      end do
      endif
  
       do m=0,3
        p834(m) = -p8(m)-p34(m)
       enddo
* pk0 -- p=p834
      p834k0=p834(0)-p834(1)
* p.q -- p.q=p834q,p=p834,q=p834,bef=,aft=
      p834q=(p834(0)*p834(0)-p834(1)*p834(1)-p834(2)*p834(2)-p83
     & 4(3)*p834(3))
  
      if (ineutri(id8).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p834,q=p8
      quqd=p834(0)*p8(0)-p834(1)*p8(1)-p834(2)*p8(2)-p834(3)*p8(
     & 3)
      ccl=fcl(id8)/(p834q*p834k0)
      do i1=1,2
* TWR0 -- qu=p834,qd=p8,v=cf34(i1)%e,a=r8_34(i1)%a,b=r8_34(i1)%b,cl=ccl,
* nsum=0
      ceps_0=-cf34(i1)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834k0*
     & (cf34(i1)%e(2)*p8(3)-p8(2)*cf34(i1)%e(3))-p8k0*(cf34(i1)%
     & e(2)*p834(3)-p834(2)*cf34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf34(i1)%e(3)*p8k0+p8(3)*cf34(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1)%e(0)*p834(0)-cf34(i1)%e(1)*p834(1)-cf34(i1)%
     & e(2)*p834(2)-cf34(i1)%e(3)*p834(3)
      cvqd=cf34(i1)%e(0)*p8(0)-cf34(i1)%e(1)*p8(1)-cf34(i1)%e(2)
     & *p8(2)-cf34(i1)%e(3)*p8(3)
      cauxa=-cf34(i1)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cf34(i1)%ek0*p8(2)+p8k0*cf34(i1)%e(2)
      r8_34(i1)%a(2)=ccl*(cauxa-ceps_0)
      r8_34(i1)%b(1)=ccl*(cauxb-ceps_2)
      end do
  
**QCD                                                                   
       if (ilept(id8).ne.1.and.ilept(id3).ne.1) then
         do iut=1,2
            rg8_34(iut)%a(1) = r8_34(iut)%a(1)
     &           /((fcl(id8))*(fcl(id3)))
            rg8_34(iut)%a(2) = r8_34(iut)%a(2)
     &           /((fcl(id8))*(fcl(id3)))
            rg8_34(iut)%b(1) = r8_34(iut)%b(1)
     &           /((fcl(id8))*(fcl(id3)))
            rg8_34(iut)%b(2) = r8_34(iut)%b(2)
     &           /((fcl(id8))*(fcl(id3)))
         enddo
       endif
      ccl=zcl(id8)/(p834q*p834k0)
      do i1=1,2
* TWR0 -- qu=p834,qd=p8,v=cz34(i1)%e,a=r8_34(i1)%a,b=r8_34(i1)%b,cl=ccl,
* nsum=1
      ceps_0=-cz34(i1)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834k0*
     & (cz34(i1)%e(2)*p8(3)-p8(2)*cz34(i1)%e(3))-p8k0*(cz34(i1)%
     & e(2)*p834(3)-p834(2)*cz34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i1)%e(3)*p8k0+p8(3)*cz34(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1)%e(0)*p834(0)-cz34(i1)%e(1)*p834(1)-cz34(i1)%
     & e(2)*p834(2)-cz34(i1)%e(3)*p834(3)
      cvqd=cz34(i1)%e(0)*p8(0)-cz34(i1)%e(1)*p8(1)-cz34(i1)%e(2)
     & *p8(2)-cz34(i1)%e(3)*p8(3)
      cauxa=-cz34(i1)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cz34(i1)%ek0*p8(2)+p8k0*cz34(i1)%e(2)
      r8_34(i1)%a(2)=r8_34(i1)%a(2)+ccl*(cauxa-ceps_0)
      r8_34(i1)%b(1)=r8_34(i1)%b(1)+ccl*(cauxb-ceps_2)
      end do
      else
* quqd -- p=p834,q=p8
      quqd=p834(0)*p8(0)-p834(1)*p8(1)-p834(2)*p8(2)-p834(3)*p8(
     & 3)
      ccl=zcl(id8)/(p834q*p834k0)
      do i1=1,2
* TWR0 -- qu=p834,qd=p8,v=cz34(i1)%e,a=r8_34(i1)%a,b=r8_34(i1)%b,cl=ccl,
* nsum=0
      ceps_0=-cz34(i1)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834k0*
     & (cz34(i1)%e(2)*p8(3)-p8(2)*cz34(i1)%e(3))-p8k0*(cz34(i1)%
     & e(2)*p834(3)-p834(2)*cz34(i1)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i1)%e(3)*p8k0+p8(3)*cz34(i1)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1)%e(0)*p834(0)-cz34(i1)%e(1)*p834(1)-cz34(i1)%
     & e(2)*p834(2)-cz34(i1)%e(3)*p834(3)
      cvqd=cz34(i1)%e(0)*p8(0)-cz34(i1)%e(1)*p8(1)-cz34(i1)%e(2)
     & *p8(2)-cz34(i1)%e(3)*p8(3)
      cauxa=-cz34(i1)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cz34(i1)%ek0*p8(2)+p8k0*cz34(i1)%e(2)
      r8_34(i1)%a(2)=ccl*(cauxa-ceps_0)
      r8_34(i1)%b(1)=ccl*(cauxb-ceps_2)
      end do
      endif
  
  
*     W insertion                                                       
       do m=0,3
        p678(m) = -p6(m)-p78(m)
       enddo
* pk0 -- p=p678
      p678k0=p678(0)-p678(1)
* p.q -- p.q=p678q,p=p678,q=p678,bef=,aft=
      p678q=(p678(0)*p678(0)-p678(1)*p678(1)-p678(2)*p678(2)-p67
     & 8(3)*p678(3))
  
* quqd -- p=p678,q=p6
      quqd=p678(0)*p6(0)-p678(1)*p6(1)-p678(2)*p6(2)-p678(3)*p6(
     & 3)
      ccl=wcl/(p678q*p678k0)
* TWR0 -- qu=p678,qd=p6,v=cw78%e,a=r6_78%a,b=r6_78%b,cl=ccl,nsum=0
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
      r6_78%a(2)=ccl*(cauxa-ceps_0)
      r6_78%b(1)=ccl*(cauxb-ceps_2)
  
       do m=0,3
        p856(m) = -p8(m)-p56(m)
       enddo
* pk0 -- p=p856
      p856k0=p856(0)-p856(1)
* p.q -- p.q=p856q,p=p856,q=p856,bef=,aft=
      p856q=(p856(0)*p856(0)-p856(1)*p856(1)-p856(2)*p856(2)-p85
     & 6(3)*p856(3))
  
* quqd -- p=p856,q=p8
      quqd=p856(0)*p8(0)-p856(1)*p8(1)-p856(2)*p8(2)-p856(3)*p8(
     & 3)
      ccl=wcl/(p856q*p856k0)
* TWR0 -- qu=p856,qd=p8,v=cw56%e,a=r8_56%a,b=r8_56%b,cl=ccl,nsum=0
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
      r8_56%a(2)=ccl*(cauxa-ceps_0)
      r8_56%b(1)=ccl*(cauxb-ceps_2)
  
  
  
*                                                                       
* MIDDLE                                                                
*                                                                       
  
* main line Z (12 or 34)                                                
  
*  Z insertion                                                          
      if (iup(id1).eq.1) then
        if (ineutri(id1-1).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p178,q=p256
      quqd=p178(0)*p256(0)-p178(1)*p256(1)-p178(2)*p256(2)-p178(
     & 3)*p256(3)
      do i3=1,2
* TW0 -- qu=p178,qd=p256,v=cf34(i3)%e,a=u178_34(i3)%a,b=u178_34(i3)%b,c=
* u178_34(i3)%c,d=u178_34(i3)%d,cl=fcl(id1-1),nsum=0
      ceps_0=-cf34(i3)%ek0*(p178(2)*p256(3)-p256(2)*p178(3))+p17
     & 8k0*(cf34(i3)%e(2)*p256(3)-p256(2)*cf34(i3)%e(3))-p256k0*
     & (cf34(i3)%e(2)*p178(3)-p178(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p178k0+p178(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p256k0+p256(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p178(0)-cf34(i3)%e(1)*p178(1)-cf34(i3)%
     & e(2)*p178(2)-cf34(i3)%e(3)*p178(3)
      cvqd=cf34(i3)%e(0)*p256(0)-cf34(i3)%e(1)*p256(1)-cf34(i3)%
     & e(2)*p256(2)-cf34(i3)%e(3)*p256(3)
      cauxa=-cf34(i3)%ek0*quqd+p178k0*cvqd+p256k0*cvqu
      cauxb=-cf34(i3)%ek0*p256(2)+p256k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p178(2)-p178k0*cf34(i3)%e(2)
      u178_34(i3)%a(2)=fcl(id1-1)*(cauxa-ceps_0)
      u178_34(i3)%b(1)=fcl(id1-1)*(cauxb-ceps_2)
      u178_34(i3)%c(2)=fcl(id1-1)*(-cauxc+ceps_1)
      u178_34(i3)%d(1)=fcl(id1-1)*cf34(i3)%ek0
      end do
  
**QCD                                                                   
          if (ilept(id1).ne.1.and.ilept(id3).ne.1) then
           do iut=1,2
             ug178_34(iut)%a(2) = u178_34(iut)%a(2)
     &           /((fcl(id1-1))*(fcl(id3)))
             ug178_34(iut)%b(1) = u178_34(iut)%b(1)
     &           /((fcl(id1-1))*(fcl(id3)))
             ug178_34(iut)%c(2) = u178_34(iut)%c(2)
     &           /((fcl(id1-1))*(fcl(id3)))
             ug178_34(iut)%d(1) = u178_34(iut)%d(1)
     &           /((fcl(id1-1))*(fcl(id3)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p178,qd=p256,v=cz34(i3)%e,a=u178_34(i3)%a,b=u178_34(i3)%b,c=
* u178_34(i3)%c,d=u178_34(i3)%d,cl=zcl(id1-1),nsum=1
      ceps_0=-cz34(i3)%ek0*(p178(2)*p256(3)-p256(2)*p178(3))+p17
     & 8k0*(cz34(i3)%e(2)*p256(3)-p256(2)*cz34(i3)%e(3))-p256k0*
     & (cz34(i3)%e(2)*p178(3)-p178(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p178k0+p178(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p256k0+p256(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p178(0)-cz34(i3)%e(1)*p178(1)-cz34(i3)%
     & e(2)*p178(2)-cz34(i3)%e(3)*p178(3)
      cvqd=cz34(i3)%e(0)*p256(0)-cz34(i3)%e(1)*p256(1)-cz34(i3)%
     & e(2)*p256(2)-cz34(i3)%e(3)*p256(3)
      cauxa=-cz34(i3)%ek0*quqd+p178k0*cvqd+p256k0*cvqu
      cauxb=-cz34(i3)%ek0*p256(2)+p256k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p178(2)-p178k0*cz34(i3)%e(2)
      u178_34(i3)%a(2)=u178_34(i3)%a(2)+zcl(id1-1)*(cauxa-ceps_0
     & )
      u178_34(i3)%b(1)=u178_34(i3)%b(1)+zcl(id1-1)*(cauxb-ceps_2
     & )
      u178_34(i3)%c(2)=u178_34(i3)%c(2)+zcl(id1-1)*(-cauxc+ceps_
     & 1)
      u178_34(i3)%d(1)=u178_34(i3)%d(1)+zcl(id1-1)*cz34(i3)%ek0
      end do
        else ! if it's a neutrino
* quqd -- p=p178,q=p256
      quqd=p178(0)*p256(0)-p178(1)*p256(1)-p178(2)*p256(2)-p178(
     & 3)*p256(3)
      do i3=1,2
* TW0 -- qu=p178,qd=p256,v=cz34(i3)%e,a=u178_34(i3)%a,b=u178_34(i3)%b,c=
* u178_34(i3)%c,d=u178_34(i3)%d,cl=zcl(id1-1),nsum=0
      ceps_0=-cz34(i3)%ek0*(p178(2)*p256(3)-p256(2)*p178(3))+p17
     & 8k0*(cz34(i3)%e(2)*p256(3)-p256(2)*cz34(i3)%e(3))-p256k0*
     & (cz34(i3)%e(2)*p178(3)-p178(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p178k0+p178(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p256k0+p256(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p178(0)-cz34(i3)%e(1)*p178(1)-cz34(i3)%
     & e(2)*p178(2)-cz34(i3)%e(3)*p178(3)
      cvqd=cz34(i3)%e(0)*p256(0)-cz34(i3)%e(1)*p256(1)-cz34(i3)%
     & e(2)*p256(2)-cz34(i3)%e(3)*p256(3)
      cauxa=-cz34(i3)%ek0*quqd+p178k0*cvqd+p256k0*cvqu
      cauxb=-cz34(i3)%ek0*p256(2)+p256k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p178(2)-p178k0*cz34(i3)%e(2)
      u178_34(i3)%a(2)=zcl(id1-1)*(cauxa-ceps_0)
      u178_34(i3)%b(1)=zcl(id1-1)*(cauxb-ceps_2)
      u178_34(i3)%c(2)=zcl(id1-1)*(-cauxc+ceps_1)
      u178_34(i3)%d(1)=zcl(id1-1)*cz34(i3)%ek0
      end do
        endif
      else  ! if it's down type
        if (ineutri(id1+1).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p156,q=p278
      quqd=p156(0)*p278(0)-p156(1)*p278(1)-p156(2)*p278(2)-p156(
     & 3)*p278(3)
      do i3=1,2
* TW0 -- qu=p156,qd=p278,v=cf34(i3)%e,a=u156_34(i3)%a,b=u156_34(i3)%b,c=
* u156_34(i3)%c,d=u156_34(i3)%d,cl=fcl(id1+1),nsum=0
      ceps_0=-cf34(i3)%ek0*(p156(2)*p278(3)-p278(2)*p156(3))+p15
     & 6k0*(cf34(i3)%e(2)*p278(3)-p278(2)*cf34(i3)%e(3))-p278k0*
     & (cf34(i3)%e(2)*p156(3)-p156(2)*cf34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i3)%e(3)*p156k0+p156(3)*cf34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i3)%e(3)*p278k0+p278(3)*cf34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i3)%e(0)*p156(0)-cf34(i3)%e(1)*p156(1)-cf34(i3)%
     & e(2)*p156(2)-cf34(i3)%e(3)*p156(3)
      cvqd=cf34(i3)%e(0)*p278(0)-cf34(i3)%e(1)*p278(1)-cf34(i3)%
     & e(2)*p278(2)-cf34(i3)%e(3)*p278(3)
      cauxa=-cf34(i3)%ek0*quqd+p156k0*cvqd+p278k0*cvqu
      cauxb=-cf34(i3)%ek0*p278(2)+p278k0*cf34(i3)%e(2)
      cauxc=+cf34(i3)%ek0*p156(2)-p156k0*cf34(i3)%e(2)
      u156_34(i3)%a(2)=fcl(id1+1)*(cauxa-ceps_0)
      u156_34(i3)%b(1)=fcl(id1+1)*(cauxb-ceps_2)
      u156_34(i3)%c(2)=fcl(id1+1)*(-cauxc+ceps_1)
      u156_34(i3)%d(1)=fcl(id1+1)*cf34(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id1).ne.1.and.ilept(id3).ne.1) then
           do iut=1,2
             ug156_34(iut)%a(2) = u156_34(iut)%a(2)
     &           /((fcl(id1+1))*(fcl(id3)))
             ug156_34(iut)%b(1) = u156_34(iut)%b(1)
     &           /((fcl(id1+1))*(fcl(id3)))
             ug156_34(iut)%c(2) = u156_34(iut)%c(2)
     &           /((fcl(id1+1))*(fcl(id3)))
             ug156_34(iut)%d(1) = u156_34(iut)%d(1)
     &           /((fcl(id1+1))*(fcl(id3)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p156,qd=p278,v=cz34(i3)%e,a=u156_34(i3)%a,b=u156_34(i3)%b,c=
* u156_34(i3)%c,d=u156_34(i3)%d,cl=zcl(id1+1),nsum=1
      ceps_0=-cz34(i3)%ek0*(p156(2)*p278(3)-p278(2)*p156(3))+p15
     & 6k0*(cz34(i3)%e(2)*p278(3)-p278(2)*cz34(i3)%e(3))-p278k0*
     & (cz34(i3)%e(2)*p156(3)-p156(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p156k0+p156(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p278k0+p278(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p156(0)-cz34(i3)%e(1)*p156(1)-cz34(i3)%
     & e(2)*p156(2)-cz34(i3)%e(3)*p156(3)
      cvqd=cz34(i3)%e(0)*p278(0)-cz34(i3)%e(1)*p278(1)-cz34(i3)%
     & e(2)*p278(2)-cz34(i3)%e(3)*p278(3)
      cauxa=-cz34(i3)%ek0*quqd+p156k0*cvqd+p278k0*cvqu
      cauxb=-cz34(i3)%ek0*p278(2)+p278k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p156(2)-p156k0*cz34(i3)%e(2)
      u156_34(i3)%a(2)=u156_34(i3)%a(2)+zcl(id1+1)*(cauxa-ceps_0
     & )
      u156_34(i3)%b(1)=u156_34(i3)%b(1)+zcl(id1+1)*(cauxb-ceps_2
     & )
      u156_34(i3)%c(2)=u156_34(i3)%c(2)+zcl(id1+1)*(-cauxc+ceps_
     & 1)
      u156_34(i3)%d(1)=u156_34(i3)%d(1)+zcl(id1+1)*cz34(i3)%ek0
      end do
        else ! if it's a neutrino
* quqd -- p=p156,q=p278
      quqd=p156(0)*p278(0)-p156(1)*p278(1)-p156(2)*p278(2)-p156(
     & 3)*p278(3)
      do i3=1,2
* TW0 -- qu=p156,qd=p278,v=cz34(i3)%e,a=u156_34(i3)%a,b=u156_34(i3)%b,c=
* u156_34(i3)%c,d=u156_34(i3)%d,cl=zcl(id1+1),nsum=0
      ceps_0=-cz34(i3)%ek0*(p156(2)*p278(3)-p278(2)*p156(3))+p15
     & 6k0*(cz34(i3)%e(2)*p278(3)-p278(2)*cz34(i3)%e(3))-p278k0*
     & (cz34(i3)%e(2)*p156(3)-p156(2)*cz34(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i3)%e(3)*p156k0+p156(3)*cz34(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i3)%e(3)*p278k0+p278(3)*cz34(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i3)%e(0)*p156(0)-cz34(i3)%e(1)*p156(1)-cz34(i3)%
     & e(2)*p156(2)-cz34(i3)%e(3)*p156(3)
      cvqd=cz34(i3)%e(0)*p278(0)-cz34(i3)%e(1)*p278(1)-cz34(i3)%
     & e(2)*p278(2)-cz34(i3)%e(3)*p278(3)
      cauxa=-cz34(i3)%ek0*quqd+p156k0*cvqd+p278k0*cvqu
      cauxb=-cz34(i3)%ek0*p278(2)+p278k0*cz34(i3)%e(2)
      cauxc=+cz34(i3)%ek0*p156(2)-p156k0*cz34(i3)%e(2)
      u156_34(i3)%a(2)=zcl(id1+1)*(cauxa-ceps_0)
      u156_34(i3)%b(1)=zcl(id1+1)*(cauxb-ceps_2)
      u156_34(i3)%c(2)=zcl(id1+1)*(-cauxc+ceps_1)
      u156_34(i3)%d(1)=zcl(id1+1)*cz34(i3)%ek0
      end do
        endif
      endif
  
*  W insertion                                                          
      if (iup(id1).eq.1) then
* quqd -- p=p178,q=p234
      quqd=p178(0)*p234(0)-p178(1)*p234(1)-p178(2)*p234(2)-p178(
     & 3)*p234(3)
* TW0 -- qu=p178,qd=p234,v=cw56%e,a=u178_56%a,b=u178_56%b,c=u178_56%c,d=
* u178_56%d,cl=wcl,nsum=0
      ceps_0=-cw56%ek0*(p178(2)*p234(3)-p234(2)*p178(3))+p178k0*
     & (cw56%e(2)*p234(3)-p234(2)*cw56%e(3))-p234k0*(cw56%e(2)*p
     & 178(3)-p178(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p178k0+p178(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p234k0+p234(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p178(0)-cw56%e(1)*p178(1)-cw56%e(2)*p178(2)
     & -cw56%e(3)*p178(3)
      cvqd=cw56%e(0)*p234(0)-cw56%e(1)*p234(1)-cw56%e(2)*p234(2)
     & -cw56%e(3)*p234(3)
      cauxa=-cw56%ek0*quqd+p178k0*cvqd+p234k0*cvqu
      cauxb=-cw56%ek0*p234(2)+p234k0*cw56%e(2)
      cauxc=+cw56%ek0*p178(2)-p178k0*cw56%e(2)
      u178_56%a(2)=wcl*(cauxa-ceps_0)
      u178_56%b(1)=wcl*(cauxb-ceps_2)
      u178_56%c(2)=wcl*(-cauxc+ceps_1)
      u178_56%d(1)=wcl*cw56%ek0
  
* quqd -- p=p134,q=p256
      quqd=p134(0)*p256(0)-p134(1)*p256(1)-p134(2)*p256(2)-p134(
     & 3)*p256(3)
* TW0 -- qu=p134,qd=p256,v=cw78%e,a=u134_78%a,b=u134_78%b,c=u134_78%c,d=
* u134_78%d,cl=wcl,nsum=0
      ceps_0=-cw78%ek0*(p134(2)*p256(3)-p256(2)*p134(3))+p134k0*
     & (cw78%e(2)*p256(3)-p256(2)*cw78%e(3))-p256k0*(cw78%e(2)*p
     & 134(3)-p134(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p134k0+p134(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p256k0+p256(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p134(0)-cw78%e(1)*p134(1)-cw78%e(2)*p134(2)
     & -cw78%e(3)*p134(3)
      cvqd=cw78%e(0)*p256(0)-cw78%e(1)*p256(1)-cw78%e(2)*p256(2)
     & -cw78%e(3)*p256(3)
      cauxa=-cw78%ek0*quqd+p134k0*cvqd+p256k0*cvqu
      cauxb=-cw78%ek0*p256(2)+p256k0*cw78%e(2)
      cauxc=+cw78%ek0*p134(2)-p134k0*cw78%e(2)
      u134_78%a(2)=wcl*(cauxa-ceps_0)
      u134_78%b(1)=wcl*(cauxb-ceps_2)
      u134_78%c(2)=wcl*(-cauxc+ceps_1)
      u134_78%d(1)=wcl*cw78%ek0
      else
* quqd -- p=p156,q=p234
      quqd=p156(0)*p234(0)-p156(1)*p234(1)-p156(2)*p234(2)-p156(
     & 3)*p234(3)
* TW0 -- qu=p156,qd=p234,v=cw78%e,a=u156_78%a,b=u156_78%b,c=u156_78%c,d=
* u156_78%d,cl=wcl,nsum=0
      ceps_0=-cw78%ek0*(p156(2)*p234(3)-p234(2)*p156(3))+p156k0*
     & (cw78%e(2)*p234(3)-p234(2)*cw78%e(3))-p234k0*(cw78%e(2)*p
     & 156(3)-p156(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p156k0+p156(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p234k0+p234(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p156(0)-cw78%e(1)*p156(1)-cw78%e(2)*p156(2)
     & -cw78%e(3)*p156(3)
      cvqd=cw78%e(0)*p234(0)-cw78%e(1)*p234(1)-cw78%e(2)*p234(2)
     & -cw78%e(3)*p234(3)
      cauxa=-cw78%ek0*quqd+p156k0*cvqd+p234k0*cvqu
      cauxb=-cw78%ek0*p234(2)+p234k0*cw78%e(2)
      cauxc=+cw78%ek0*p156(2)-p156k0*cw78%e(2)
      u156_78%a(2)=wcl*(cauxa-ceps_0)
      u156_78%b(1)=wcl*(cauxb-ceps_2)
      u156_78%c(2)=wcl*(-cauxc+ceps_1)
      u156_78%d(1)=wcl*cw78%ek0
  
* quqd -- p=p134,q=p278
      quqd=p134(0)*p278(0)-p134(1)*p278(1)-p134(2)*p278(2)-p134(
     & 3)*p278(3)
* TW0 -- qu=p134,qd=p278,v=cw56%e,a=u134_56%a,b=u134_56%b,c=u134_56%c,d=
* u134_56%d,cl=wcl,nsum=0
      ceps_0=-cw56%ek0*(p134(2)*p278(3)-p278(2)*p134(3))+p134k0*
     & (cw56%e(2)*p278(3)-p278(2)*cw56%e(3))-p278k0*(cw56%e(2)*p
     & 134(3)-p134(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p134k0+p134(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p278k0+p278(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p134(0)-cw56%e(1)*p134(1)-cw56%e(2)*p134(2)
     & -cw56%e(3)*p134(3)
      cvqd=cw56%e(0)*p278(0)-cw56%e(1)*p278(1)-cw56%e(2)*p278(2)
     & -cw56%e(3)*p278(3)
      cauxa=-cw56%ek0*quqd+p134k0*cvqd+p278k0*cvqu
      cauxb=-cw56%ek0*p278(2)+p278k0*cw56%e(2)
      cauxc=+cw56%ek0*p134(2)-p134k0*cw56%e(2)
      u134_56%a(2)=wcl*(cauxa-ceps_0)
      u134_56%b(1)=wcl*(cauxb-ceps_2)
      u134_56%c(2)=wcl*(-cauxc+ceps_1)
      u134_56%d(1)=wcl*cw56%ek0
      endif
  
  
*  Z insertion                                                          
      if (iup(id3).eq.1) then
        if (ineutri(id3-1).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p378,q=p456
      quqd=p378(0)*p456(0)-p378(1)*p456(1)-p378(2)*p456(2)-p378(
     & 3)*p456(3)
      do i3=1,2
* TW0 -- qu=p378,qd=p456,v=cf12(i3)%e,a=u378_12(i3)%a,b=u378_12(i3)%b,c=
* u378_12(i3)%c,d=u378_12(i3)%d,cl=fcl(id3-1),nsum=0
      ceps_0=-cf12(i3)%ek0*(p378(2)*p456(3)-p456(2)*p378(3))+p37
     & 8k0*(cf12(i3)%e(2)*p456(3)-p456(2)*cf12(i3)%e(3))-p456k0*
     & (cf12(i3)%e(2)*p378(3)-p378(2)*cf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i3)%e(3)*p378k0+p378(3)*cf12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i3)%e(3)*p456k0+p456(3)*cf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i3)%e(0)*p378(0)-cf12(i3)%e(1)*p378(1)-cf12(i3)%
     & e(2)*p378(2)-cf12(i3)%e(3)*p378(3)
      cvqd=cf12(i3)%e(0)*p456(0)-cf12(i3)%e(1)*p456(1)-cf12(i3)%
     & e(2)*p456(2)-cf12(i3)%e(3)*p456(3)
      cauxa=-cf12(i3)%ek0*quqd+p378k0*cvqd+p456k0*cvqu
      cauxb=-cf12(i3)%ek0*p456(2)+p456k0*cf12(i3)%e(2)
      cauxc=+cf12(i3)%ek0*p378(2)-p378k0*cf12(i3)%e(2)
      u378_12(i3)%a(2)=fcl(id3-1)*(cauxa-ceps_0)
      u378_12(i3)%b(1)=fcl(id3-1)*(cauxb-ceps_2)
      u378_12(i3)%c(2)=fcl(id3-1)*(-cauxc+ceps_1)
      u378_12(i3)%d(1)=fcl(id3-1)*cf12(i3)%ek0
      end do
  
**QCD                                                                   
          if (ilept(id3).ne.1.and.ilept(id1).ne.1) then
           do iut=1,2
             ug378_12(iut)%a(2) = u378_12(iut)%a(2)
     &           /((fcl(id3-1))*(fcl(id1)))
             ug378_12(iut)%b(1) = u378_12(iut)%b(1)
     &           /((fcl(id3-1))*(fcl(id1)))
             ug378_12(iut)%c(2) = u378_12(iut)%c(2)
     &           /((fcl(id3-1))*(fcl(id1)))
             ug378_12(iut)%d(1) = u378_12(iut)%d(1)
     &           /((fcl(id3-1))*(fcl(id1)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p378,qd=p456,v=cz12(i3)%e,a=u378_12(i3)%a,b=u378_12(i3)%b,c=
* u378_12(i3)%c,d=u378_12(i3)%d,cl=zcl(id3-1),nsum=1
      ceps_0=-cz12(i3)%ek0*(p378(2)*p456(3)-p456(2)*p378(3))+p37
     & 8k0*(cz12(i3)%e(2)*p456(3)-p456(2)*cz12(i3)%e(3))-p456k0*
     & (cz12(i3)%e(2)*p378(3)-p378(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p378k0+p378(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p456k0+p456(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p378(0)-cz12(i3)%e(1)*p378(1)-cz12(i3)%
     & e(2)*p378(2)-cz12(i3)%e(3)*p378(3)
      cvqd=cz12(i3)%e(0)*p456(0)-cz12(i3)%e(1)*p456(1)-cz12(i3)%
     & e(2)*p456(2)-cz12(i3)%e(3)*p456(3)
      cauxa=-cz12(i3)%ek0*quqd+p378k0*cvqd+p456k0*cvqu
      cauxb=-cz12(i3)%ek0*p456(2)+p456k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p378(2)-p378k0*cz12(i3)%e(2)
      u378_12(i3)%a(2)=u378_12(i3)%a(2)+zcl(id3-1)*(cauxa-ceps_0
     & )
      u378_12(i3)%b(1)=u378_12(i3)%b(1)+zcl(id3-1)*(cauxb-ceps_2
     & )
      u378_12(i3)%c(2)=u378_12(i3)%c(2)+zcl(id3-1)*(-cauxc+ceps_
     & 1)
      u378_12(i3)%d(1)=u378_12(i3)%d(1)+zcl(id3-1)*cz12(i3)%ek0
      end do
        else ! if it's a neutrino
* quqd -- p=p378,q=p456
      quqd=p378(0)*p456(0)-p378(1)*p456(1)-p378(2)*p456(2)-p378(
     & 3)*p456(3)
      do i3=1,2
* TW0 -- qu=p378,qd=p456,v=cz12(i3)%e,a=u378_12(i3)%a,b=u378_12(i3)%b,c=
* u378_12(i3)%c,d=u378_12(i3)%d,cl=zcl(id3-1),nsum=0
      ceps_0=-cz12(i3)%ek0*(p378(2)*p456(3)-p456(2)*p378(3))+p37
     & 8k0*(cz12(i3)%e(2)*p456(3)-p456(2)*cz12(i3)%e(3))-p456k0*
     & (cz12(i3)%e(2)*p378(3)-p378(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p378k0+p378(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p456k0+p456(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p378(0)-cz12(i3)%e(1)*p378(1)-cz12(i3)%
     & e(2)*p378(2)-cz12(i3)%e(3)*p378(3)
      cvqd=cz12(i3)%e(0)*p456(0)-cz12(i3)%e(1)*p456(1)-cz12(i3)%
     & e(2)*p456(2)-cz12(i3)%e(3)*p456(3)
      cauxa=-cz12(i3)%ek0*quqd+p378k0*cvqd+p456k0*cvqu
      cauxb=-cz12(i3)%ek0*p456(2)+p456k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p378(2)-p378k0*cz12(i3)%e(2)
      u378_12(i3)%a(2)=zcl(id3-1)*(cauxa-ceps_0)
      u378_12(i3)%b(1)=zcl(id3-1)*(cauxb-ceps_2)
      u378_12(i3)%c(2)=zcl(id3-1)*(-cauxc+ceps_1)
      u378_12(i3)%d(1)=zcl(id3-1)*cz12(i3)%ek0
      end do
        endif
      else  ! if it's down type
        if (ineutri(id3+1).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p356,q=p478
      quqd=p356(0)*p478(0)-p356(1)*p478(1)-p356(2)*p478(2)-p356(
     & 3)*p478(3)
      do i3=1,2
* TW0 -- qu=p356,qd=p478,v=cf12(i3)%e,a=u356_12(i3)%a,b=u356_12(i3)%b,c=
* u356_12(i3)%c,d=u356_12(i3)%d,cl=fcl(id3+1),nsum=0
      ceps_0=-cf12(i3)%ek0*(p356(2)*p478(3)-p478(2)*p356(3))+p35
     & 6k0*(cf12(i3)%e(2)*p478(3)-p478(2)*cf12(i3)%e(3))-p478k0*
     & (cf12(i3)%e(2)*p356(3)-p356(2)*cf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i3)%e(3)*p356k0+p356(3)*cf12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i3)%e(3)*p478k0+p478(3)*cf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i3)%e(0)*p356(0)-cf12(i3)%e(1)*p356(1)-cf12(i3)%
     & e(2)*p356(2)-cf12(i3)%e(3)*p356(3)
      cvqd=cf12(i3)%e(0)*p478(0)-cf12(i3)%e(1)*p478(1)-cf12(i3)%
     & e(2)*p478(2)-cf12(i3)%e(3)*p478(3)
      cauxa=-cf12(i3)%ek0*quqd+p356k0*cvqd+p478k0*cvqu
      cauxb=-cf12(i3)%ek0*p478(2)+p478k0*cf12(i3)%e(2)
      cauxc=+cf12(i3)%ek0*p356(2)-p356k0*cf12(i3)%e(2)
      u356_12(i3)%a(2)=fcl(id3+1)*(cauxa-ceps_0)
      u356_12(i3)%b(1)=fcl(id3+1)*(cauxb-ceps_2)
      u356_12(i3)%c(2)=fcl(id3+1)*(-cauxc+ceps_1)
      u356_12(i3)%d(1)=fcl(id3+1)*cf12(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id3).ne.1.and.ilept(id1).ne.1) then
           do iut=1,2
             ug356_12(iut)%a(2) = u356_12(iut)%a(2)
     &           /((fcl(id3+1))*(fcl(id1)))
             ug356_12(iut)%b(1) = u356_12(iut)%b(1)
     &           /((fcl(id3+1))*(fcl(id1)))
             ug356_12(iut)%c(2) = u356_12(iut)%c(2)
     &           /((fcl(id3+1))*(fcl(id1)))
             ug356_12(iut)%d(1) = u356_12(iut)%d(1)
     &           /((fcl(id3+1))*(fcl(id1)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p356,qd=p478,v=cz12(i3)%e,a=u356_12(i3)%a,b=u356_12(i3)%b,c=
* u356_12(i3)%c,d=u356_12(i3)%d,cl=zcl(id3+1),nsum=1
      ceps_0=-cz12(i3)%ek0*(p356(2)*p478(3)-p478(2)*p356(3))+p35
     & 6k0*(cz12(i3)%e(2)*p478(3)-p478(2)*cz12(i3)%e(3))-p478k0*
     & (cz12(i3)%e(2)*p356(3)-p356(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p356k0+p356(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p478k0+p478(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p356(0)-cz12(i3)%e(1)*p356(1)-cz12(i3)%
     & e(2)*p356(2)-cz12(i3)%e(3)*p356(3)
      cvqd=cz12(i3)%e(0)*p478(0)-cz12(i3)%e(1)*p478(1)-cz12(i3)%
     & e(2)*p478(2)-cz12(i3)%e(3)*p478(3)
      cauxa=-cz12(i3)%ek0*quqd+p356k0*cvqd+p478k0*cvqu
      cauxb=-cz12(i3)%ek0*p478(2)+p478k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p356(2)-p356k0*cz12(i3)%e(2)
      u356_12(i3)%a(2)=u356_12(i3)%a(2)+zcl(id3+1)*(cauxa-ceps_0
     & )
      u356_12(i3)%b(1)=u356_12(i3)%b(1)+zcl(id3+1)*(cauxb-ceps_2
     & )
      u356_12(i3)%c(2)=u356_12(i3)%c(2)+zcl(id3+1)*(-cauxc+ceps_
     & 1)
      u356_12(i3)%d(1)=u356_12(i3)%d(1)+zcl(id3+1)*cz12(i3)%ek0
      end do
        else ! if it's a neutrino
* quqd -- p=p356,q=p478
      quqd=p356(0)*p478(0)-p356(1)*p478(1)-p356(2)*p478(2)-p356(
     & 3)*p478(3)
      do i3=1,2
* TW0 -- qu=p356,qd=p478,v=cz12(i3)%e,a=u356_12(i3)%a,b=u356_12(i3)%b,c=
* u356_12(i3)%c,d=u356_12(i3)%d,cl=zcl(id3+1),nsum=0
      ceps_0=-cz12(i3)%ek0*(p356(2)*p478(3)-p478(2)*p356(3))+p35
     & 6k0*(cz12(i3)%e(2)*p478(3)-p478(2)*cz12(i3)%e(3))-p478k0*
     & (cz12(i3)%e(2)*p356(3)-p356(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p356k0+p356(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p478k0+p478(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p356(0)-cz12(i3)%e(1)*p356(1)-cz12(i3)%
     & e(2)*p356(2)-cz12(i3)%e(3)*p356(3)
      cvqd=cz12(i3)%e(0)*p478(0)-cz12(i3)%e(1)*p478(1)-cz12(i3)%
     & e(2)*p478(2)-cz12(i3)%e(3)*p478(3)
      cauxa=-cz12(i3)%ek0*quqd+p356k0*cvqd+p478k0*cvqu
      cauxb=-cz12(i3)%ek0*p478(2)+p478k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p356(2)-p356k0*cz12(i3)%e(2)
      u356_12(i3)%a(2)=zcl(id3+1)*(cauxa-ceps_0)
      u356_12(i3)%b(1)=zcl(id3+1)*(cauxb-ceps_2)
      u356_12(i3)%c(2)=zcl(id3+1)*(-cauxc+ceps_1)
      u356_12(i3)%d(1)=zcl(id3+1)*cz12(i3)%ek0
      end do
        endif
      endif
  
*  W insertion                                                          
      if (iup(id3).eq.1) then
* quqd -- p=p378,q=p412
      quqd=p378(0)*p412(0)-p378(1)*p412(1)-p378(2)*p412(2)-p378(
     & 3)*p412(3)
* TW0 -- qu=p378,qd=p412,v=cw56%e,a=u378_56%a,b=u378_56%b,c=u378_56%c,d=
* u378_56%d,cl=wcl,nsum=0
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
      u378_56%a(2)=wcl*(cauxa-ceps_0)
      u378_56%b(1)=wcl*(cauxb-ceps_2)
      u378_56%c(2)=wcl*(-cauxc+ceps_1)
      u378_56%d(1)=wcl*cw56%ek0
  
* quqd -- p=p312,q=p456
      quqd=p312(0)*p456(0)-p312(1)*p456(1)-p312(2)*p456(2)-p312(
     & 3)*p456(3)
* TW0 -- qu=p312,qd=p456,v=cw78%e,a=u312_78%a,b=u312_78%b,c=u312_78%c,d=
* u312_78%d,cl=wcl,nsum=0
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
      u312_78%a(2)=wcl*(cauxa-ceps_0)
      u312_78%b(1)=wcl*(cauxb-ceps_2)
      u312_78%c(2)=wcl*(-cauxc+ceps_1)
      u312_78%d(1)=wcl*cw78%ek0
      else
* quqd -- p=p356,q=p412
      quqd=p356(0)*p412(0)-p356(1)*p412(1)-p356(2)*p412(2)-p356(
     & 3)*p412(3)
* TW0 -- qu=p356,qd=p412,v=cw78%e,a=u356_78%a,b=u356_78%b,c=u356_78%c,d=
* u356_78%d,cl=wcl,nsum=0
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
      u356_78%a(2)=wcl*(cauxa-ceps_0)
      u356_78%b(1)=wcl*(cauxb-ceps_2)
      u356_78%c(2)=wcl*(-cauxc+ceps_1)
      u356_78%d(1)=wcl*cw78%ek0
  
* quqd -- p=p312,q=p478
      quqd=p312(0)*p478(0)-p312(1)*p478(1)-p312(2)*p478(2)-p312(
     & 3)*p478(3)
* TW0 -- qu=p312,qd=p478,v=cw56%e,a=u312_56%a,b=u312_56%b,c=u312_56%c,d=
* u312_56%d,cl=wcl,nsum=0
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
      u312_56%a(2)=wcl*(cauxa-ceps_0)
      u312_56%b(1)=wcl*(cauxb-ceps_2)
      u312_56%c(2)=wcl*(-cauxc+ceps_1)
      u312_56%d(1)=wcl*cw56%ek0
      endif
  
  
* main line W (56 or 78)                                                
*  Z insertion                                                          
      if (ineutri(id6).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p578,q=p612
      quqd=p578(0)*p612(0)-p578(1)*p612(1)-p578(2)*p612(2)-p578(
     & 3)*p612(3)
      do i3=1,2
* TW0 -- qu=p578,qd=p612,v=cf34(i3)%e,a=u578_34(i3)%a,b=u578_34(i3)%b,c=
* u578_34(i3)%c,d=u578_34(i3)%d,cl=fcl(id6),nsum=0
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
      u578_34(i3)%a(2)=fcl(id6)*(cauxa-ceps_0)
      u578_34(i3)%b(1)=fcl(id6)*(cauxb-ceps_2)
      u578_34(i3)%c(2)=fcl(id6)*(-cauxc+ceps_1)
      u578_34(i3)%d(1)=fcl(id6)*cf34(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id6).ne.1.and.ilept(id3).ne.1) then
           do iut=1,2
             ug578_34(iut)%a(2) = u578_34(iut)%a(2)
     &           /((fcl(id6))*(fcl(id3)))
             ug578_34(iut)%b(1) = u578_34(iut)%b(1)
     &           /((fcl(id6))*(fcl(id3)))
             ug578_34(iut)%c(2) = u578_34(iut)%c(2)
     &           /((fcl(id6))*(fcl(id3)))
             ug578_34(iut)%d(1) = u578_34(iut)%d(1)
     &           /((fcl(id6))*(fcl(id3)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p578,qd=p612,v=cz34(i3)%e,a=u578_34(i3)%a,b=u578_34(i3)%b,c=
* u578_34(i3)%c,d=u578_34(i3)%d,cl=zcl(id6),nsum=1
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
      u578_34(i3)%a(2)=u578_34(i3)%a(2)+zcl(id6)*(cauxa-ceps_0)
      u578_34(i3)%b(1)=u578_34(i3)%b(1)+zcl(id6)*(cauxb-ceps_2)
      u578_34(i3)%c(2)=u578_34(i3)%c(2)+zcl(id6)*(-cauxc+ceps_1)
      u578_34(i3)%d(1)=u578_34(i3)%d(1)+zcl(id6)*cz34(i3)%ek0
      end do
      else
* quqd -- p=p578,q=p612
      quqd=p578(0)*p612(0)-p578(1)*p612(1)-p578(2)*p612(2)-p578(
     & 3)*p612(3)
      do i3=1,2
* TW0 -- qu=p578,qd=p612,v=cz34(i3)%e,a=u578_34(i3)%a,b=u578_34(i3)%b,c=
* u578_34(i3)%c,d=u578_34(i3)%d,cl=zcl(id6),nsum=0
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
      u578_34(i3)%a(2)=zcl(id6)*(cauxa-ceps_0)
      u578_34(i3)%b(1)=zcl(id6)*(cauxb-ceps_2)
      u578_34(i3)%c(2)=zcl(id6)*(-cauxc+ceps_1)
      u578_34(i3)%d(1)=zcl(id6)*cz34(i3)%ek0
      end do
      endif
  
      if (ineutri(id5).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p512,q=p678
      quqd=p512(0)*p678(0)-p512(1)*p678(1)-p512(2)*p678(2)-p512(
     & 3)*p678(3)
      do i3=1,2
* TW0 -- qu=p512,qd=p678,v=cf34(i3)%e,a=u512_34(i3)%a,b=u512_34(i3)%b,c=
* u512_34(i3)%c,d=u512_34(i3)%d,cl=fcl(id5),nsum=0
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
      u512_34(i3)%a(2)=fcl(id5)*(cauxa-ceps_0)
      u512_34(i3)%b(1)=fcl(id5)*(cauxb-ceps_2)
      u512_34(i3)%c(2)=fcl(id5)*(-cauxc+ceps_1)
      u512_34(i3)%d(1)=fcl(id5)*cf34(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id5).ne.1.and.ilept(id3).ne.1) then
           do iut=1,2
             ug512_34(iut)%a(2) = u512_34(iut)%a(2)
     &           /((fcl(id5))*(fcl(id3)))
             ug512_34(iut)%b(1) = u512_34(iut)%b(1)
     &           /((fcl(id5))*(fcl(id3)))
             ug512_34(iut)%c(2) = u512_34(iut)%c(2)
     &           /((fcl(id5))*(fcl(id3)))
             ug512_34(iut)%d(1) = u512_34(iut)%d(1)
     &           /((fcl(id5))*(fcl(id3)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p512,qd=p678,v=cz34(i3)%e,a=u512_34(i3)%a,b=u512_34(i3)%b,c=
* u512_34(i3)%c,d=u512_34(i3)%d,cl=zcl(id5),nsum=1
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
      u512_34(i3)%a(2)=u512_34(i3)%a(2)+zcl(id5)*(cauxa-ceps_0)
      u512_34(i3)%b(1)=u512_34(i3)%b(1)+zcl(id5)*(cauxb-ceps_2)
      u512_34(i3)%c(2)=u512_34(i3)%c(2)+zcl(id5)*(-cauxc+ceps_1)
      u512_34(i3)%d(1)=u512_34(i3)%d(1)+zcl(id5)*cz34(i3)%ek0
      end do
      else
* quqd -- p=p512,q=p678
      quqd=p512(0)*p678(0)-p512(1)*p678(1)-p512(2)*p678(2)-p512(
     & 3)*p678(3)
      do i3=1,2
* TW0 -- qu=p512,qd=p678,v=cz34(i3)%e,a=u512_34(i3)%a,b=u512_34(i3)%b,c=
* u512_34(i3)%c,d=u512_34(i3)%d,cl=zcl(id5),nsum=0
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
      u512_34(i3)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u512_34(i3)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u512_34(i3)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u512_34(i3)%d(1)=zcl(id5)*cz34(i3)%ek0
      end do
      endif
  
*  W insertion                                                          
* quqd -- p=p512,q=p634
      quqd=p512(0)*p634(0)-p512(1)*p634(1)-p512(2)*p634(2)-p512(
     & 3)*p634(3)
* TW0 -- qu=p512,qd=p634,v=cw78%e,a=u512_78%a,b=u512_78%b,c=u512_78%c,d=
* u512_78%d,cl=wcl,nsum=0
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
      u512_78%a(2)=wcl*(cauxa-ceps_0)
      u512_78%b(1)=wcl*(cauxb-ceps_2)
      u512_78%c(2)=wcl*(-cauxc+ceps_1)
      u512_78%d(1)=wcl*cw78%ek0
  
*  Z insertion                                                          
      if (ineutri(id6).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p578,q=p634
      quqd=p578(0)*p634(0)-p578(1)*p634(1)-p578(2)*p634(2)-p578(
     & 3)*p634(3)
      do i3=1,2
* TW0 -- qu=p578,qd=p634,v=cf12(i3)%e,a=u578_12(i3)%a,b=u578_12(i3)%b,c=
* u578_12(i3)%c,d=u578_12(i3)%d,cl=fcl(id6),nsum=0
      ceps_0=-cf12(i3)%ek0*(p578(2)*p634(3)-p634(2)*p578(3))+p57
     & 8k0*(cf12(i3)%e(2)*p634(3)-p634(2)*cf12(i3)%e(3))-p634k0*
     & (cf12(i3)%e(2)*p578(3)-p578(2)*cf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i3)%e(3)*p578k0+p578(3)*cf12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i3)%e(3)*p634k0+p634(3)*cf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i3)%e(0)*p578(0)-cf12(i3)%e(1)*p578(1)-cf12(i3)%
     & e(2)*p578(2)-cf12(i3)%e(3)*p578(3)
      cvqd=cf12(i3)%e(0)*p634(0)-cf12(i3)%e(1)*p634(1)-cf12(i3)%
     & e(2)*p634(2)-cf12(i3)%e(3)*p634(3)
      cauxa=-cf12(i3)%ek0*quqd+p578k0*cvqd+p634k0*cvqu
      cauxb=-cf12(i3)%ek0*p634(2)+p634k0*cf12(i3)%e(2)
      cauxc=+cf12(i3)%ek0*p578(2)-p578k0*cf12(i3)%e(2)
      u578_12(i3)%a(2)=fcl(id6)*(cauxa-ceps_0)
      u578_12(i3)%b(1)=fcl(id6)*(cauxb-ceps_2)
      u578_12(i3)%c(2)=fcl(id6)*(-cauxc+ceps_1)
      u578_12(i3)%d(1)=fcl(id6)*cf12(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id6).ne.1.and.ilept(id1).ne.1) then
           do iut=1,2
             ug578_12(iut)%a(2) = u578_12(iut)%a(2)
     &           /((fcl(id6))*(fcl(id1)))
             ug578_12(iut)%b(1) = u578_12(iut)%b(1)
     &           /((fcl(id6))*(fcl(id1)))
             ug578_12(iut)%c(2) = u578_12(iut)%c(2)
     &           /((fcl(id6))*(fcl(id1)))
             ug578_12(iut)%d(1) = u578_12(iut)%d(1)
     &           /((fcl(id6))*(fcl(id1)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p578,qd=p634,v=cz12(i3)%e,a=u578_12(i3)%a,b=u578_12(i3)%b,c=
* u578_12(i3)%c,d=u578_12(i3)%d,cl=zcl(id6),nsum=1
      ceps_0=-cz12(i3)%ek0*(p578(2)*p634(3)-p634(2)*p578(3))+p57
     & 8k0*(cz12(i3)%e(2)*p634(3)-p634(2)*cz12(i3)%e(3))-p634k0*
     & (cz12(i3)%e(2)*p578(3)-p578(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p578k0+p578(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p634k0+p634(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p578(0)-cz12(i3)%e(1)*p578(1)-cz12(i3)%
     & e(2)*p578(2)-cz12(i3)%e(3)*p578(3)
      cvqd=cz12(i3)%e(0)*p634(0)-cz12(i3)%e(1)*p634(1)-cz12(i3)%
     & e(2)*p634(2)-cz12(i3)%e(3)*p634(3)
      cauxa=-cz12(i3)%ek0*quqd+p578k0*cvqd+p634k0*cvqu
      cauxb=-cz12(i3)%ek0*p634(2)+p634k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p578(2)-p578k0*cz12(i3)%e(2)
      u578_12(i3)%a(2)=u578_12(i3)%a(2)+zcl(id6)*(cauxa-ceps_0)
      u578_12(i3)%b(1)=u578_12(i3)%b(1)+zcl(id6)*(cauxb-ceps_2)
      u578_12(i3)%c(2)=u578_12(i3)%c(2)+zcl(id6)*(-cauxc+ceps_1)
      u578_12(i3)%d(1)=u578_12(i3)%d(1)+zcl(id6)*cz12(i3)%ek0
      end do
      else
* quqd -- p=p578,q=p634
      quqd=p578(0)*p634(0)-p578(1)*p634(1)-p578(2)*p634(2)-p578(
     & 3)*p634(3)
      do i3=1,2
* TW0 -- qu=p578,qd=p634,v=cz12(i3)%e,a=u578_12(i3)%a,b=u578_12(i3)%b,c=
* u578_12(i3)%c,d=u578_12(i3)%d,cl=zcl(id6),nsum=0
      ceps_0=-cz12(i3)%ek0*(p578(2)*p634(3)-p634(2)*p578(3))+p57
     & 8k0*(cz12(i3)%e(2)*p634(3)-p634(2)*cz12(i3)%e(3))-p634k0*
     & (cz12(i3)%e(2)*p578(3)-p578(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p578k0+p578(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p634k0+p634(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p578(0)-cz12(i3)%e(1)*p578(1)-cz12(i3)%
     & e(2)*p578(2)-cz12(i3)%e(3)*p578(3)
      cvqd=cz12(i3)%e(0)*p634(0)-cz12(i3)%e(1)*p634(1)-cz12(i3)%
     & e(2)*p634(2)-cz12(i3)%e(3)*p634(3)
      cauxa=-cz12(i3)%ek0*quqd+p578k0*cvqd+p634k0*cvqu
      cauxb=-cz12(i3)%ek0*p634(2)+p634k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p578(2)-p578k0*cz12(i3)%e(2)
      u578_12(i3)%a(2)=zcl(id6)*(cauxa-ceps_0)
      u578_12(i3)%b(1)=zcl(id6)*(cauxb-ceps_2)
      u578_12(i3)%c(2)=zcl(id6)*(-cauxc+ceps_1)
      u578_12(i3)%d(1)=zcl(id6)*cz12(i3)%ek0
      end do
      endif
  
      if (ineutri(id5).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p534,q=p678
      quqd=p534(0)*p678(0)-p534(1)*p678(1)-p534(2)*p678(2)-p534(
     & 3)*p678(3)
      do i3=1,2
* TW0 -- qu=p534,qd=p678,v=cf12(i3)%e,a=u534_12(i3)%a,b=u534_12(i3)%b,c=
* u534_12(i3)%c,d=u534_12(i3)%d,cl=fcl(id5),nsum=0
      ceps_0=-cf12(i3)%ek0*(p534(2)*p678(3)-p678(2)*p534(3))+p53
     & 4k0*(cf12(i3)%e(2)*p678(3)-p678(2)*cf12(i3)%e(3))-p678k0*
     & (cf12(i3)%e(2)*p534(3)-p534(2)*cf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i3)%e(3)*p534k0+p534(3)*cf12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i3)%e(3)*p678k0+p678(3)*cf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i3)%e(0)*p534(0)-cf12(i3)%e(1)*p534(1)-cf12(i3)%
     & e(2)*p534(2)-cf12(i3)%e(3)*p534(3)
      cvqd=cf12(i3)%e(0)*p678(0)-cf12(i3)%e(1)*p678(1)-cf12(i3)%
     & e(2)*p678(2)-cf12(i3)%e(3)*p678(3)
      cauxa=-cf12(i3)%ek0*quqd+p534k0*cvqd+p678k0*cvqu
      cauxb=-cf12(i3)%ek0*p678(2)+p678k0*cf12(i3)%e(2)
      cauxc=+cf12(i3)%ek0*p534(2)-p534k0*cf12(i3)%e(2)
      u534_12(i3)%a(2)=fcl(id5)*(cauxa-ceps_0)
      u534_12(i3)%b(1)=fcl(id5)*(cauxb-ceps_2)
      u534_12(i3)%c(2)=fcl(id5)*(-cauxc+ceps_1)
      u534_12(i3)%d(1)=fcl(id5)*cf12(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id5).ne.1.and.ilept(id1).ne.1) then
           do iut=1,2
             ug534_12(iut)%a(2) = u534_12(iut)%a(2)
     &           /((fcl(id5))*(fcl(id1)))
             ug534_12(iut)%b(1) = u534_12(iut)%b(1)
     &           /((fcl(id5))*(fcl(id1)))
             ug534_12(iut)%c(2) = u534_12(iut)%c(2)
     &           /((fcl(id5))*(fcl(id1)))
             ug534_12(iut)%d(1) = u534_12(iut)%d(1)
     &           /((fcl(id5))*(fcl(id1)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p534,qd=p678,v=cz12(i3)%e,a=u534_12(i3)%a,b=u534_12(i3)%b,c=
* u534_12(i3)%c,d=u534_12(i3)%d,cl=zcl(id5),nsum=1
      ceps_0=-cz12(i3)%ek0*(p534(2)*p678(3)-p678(2)*p534(3))+p53
     & 4k0*(cz12(i3)%e(2)*p678(3)-p678(2)*cz12(i3)%e(3))-p678k0*
     & (cz12(i3)%e(2)*p534(3)-p534(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p534k0+p534(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p678k0+p678(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p534(0)-cz12(i3)%e(1)*p534(1)-cz12(i3)%
     & e(2)*p534(2)-cz12(i3)%e(3)*p534(3)
      cvqd=cz12(i3)%e(0)*p678(0)-cz12(i3)%e(1)*p678(1)-cz12(i3)%
     & e(2)*p678(2)-cz12(i3)%e(3)*p678(3)
      cauxa=-cz12(i3)%ek0*quqd+p534k0*cvqd+p678k0*cvqu
      cauxb=-cz12(i3)%ek0*p678(2)+p678k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p534(2)-p534k0*cz12(i3)%e(2)
      u534_12(i3)%a(2)=u534_12(i3)%a(2)+zcl(id5)*(cauxa-ceps_0)
      u534_12(i3)%b(1)=u534_12(i3)%b(1)+zcl(id5)*(cauxb-ceps_2)
      u534_12(i3)%c(2)=u534_12(i3)%c(2)+zcl(id5)*(-cauxc+ceps_1)
      u534_12(i3)%d(1)=u534_12(i3)%d(1)+zcl(id5)*cz12(i3)%ek0
      end do
      else
* quqd -- p=p534,q=p678
      quqd=p534(0)*p678(0)-p534(1)*p678(1)-p534(2)*p678(2)-p534(
     & 3)*p678(3)
      do i3=1,2
* TW0 -- qu=p534,qd=p678,v=cz12(i3)%e,a=u534_12(i3)%a,b=u534_12(i3)%b,c=
* u534_12(i3)%c,d=u534_12(i3)%d,cl=zcl(id5),nsum=0
      ceps_0=-cz12(i3)%ek0*(p534(2)*p678(3)-p678(2)*p534(3))+p53
     & 4k0*(cz12(i3)%e(2)*p678(3)-p678(2)*cz12(i3)%e(3))-p678k0*
     & (cz12(i3)%e(2)*p534(3)-p534(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p534k0+p534(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p678k0+p678(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p534(0)-cz12(i3)%e(1)*p534(1)-cz12(i3)%
     & e(2)*p534(2)-cz12(i3)%e(3)*p534(3)
      cvqd=cz12(i3)%e(0)*p678(0)-cz12(i3)%e(1)*p678(1)-cz12(i3)%
     & e(2)*p678(2)-cz12(i3)%e(3)*p678(3)
      cauxa=-cz12(i3)%ek0*quqd+p534k0*cvqd+p678k0*cvqu
      cauxb=-cz12(i3)%ek0*p678(2)+p678k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p534(2)-p534k0*cz12(i3)%e(2)
      u534_12(i3)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u534_12(i3)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u534_12(i3)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u534_12(i3)%d(1)=zcl(id5)*cz12(i3)%ek0
      end do
      endif
  
*  W insertion                                                          
* quqd -- p=p534,q=p612
      quqd=p534(0)*p612(0)-p534(1)*p612(1)-p534(2)*p612(2)-p534(
     & 3)*p612(3)
* TW0 -- qu=p534,qd=p612,v=cw78%e,a=u534_78%a,b=u534_78%b,c=u534_78%c,d=
* u534_78%d,cl=wcl,nsum=0
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
      u534_78%a(2)=wcl*(cauxa-ceps_0)
      u534_78%b(1)=wcl*(cauxb-ceps_2)
      u534_78%c(2)=wcl*(-cauxc+ceps_1)
      u534_78%d(1)=wcl*cw78%ek0
  
*  Z insertion                                                          
      if (ineutri(id8).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p756,q=p812
      quqd=p756(0)*p812(0)-p756(1)*p812(1)-p756(2)*p812(2)-p756(
     & 3)*p812(3)
      do i3=1,2
* TW0 -- qu=p756,qd=p812,v=cf34(i3)%e,a=u756_34(i3)%a,b=u756_34(i3)%b,c=
* u756_34(i3)%c,d=u756_34(i3)%d,cl=fcl(id8),nsum=0
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
      u756_34(i3)%a(2)=fcl(id8)*(cauxa-ceps_0)
      u756_34(i3)%b(1)=fcl(id8)*(cauxb-ceps_2)
      u756_34(i3)%c(2)=fcl(id8)*(-cauxc+ceps_1)
      u756_34(i3)%d(1)=fcl(id8)*cf34(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id8).ne.1.and.ilept(id3).ne.1) then
           do iut=1,2
             ug756_34(iut)%a(2) = u756_34(iut)%a(2)
     &           /((fcl(id8))*(fcl(id3)))
             ug756_34(iut)%b(1) = u756_34(iut)%b(1)
     &           /((fcl(id8))*(fcl(id3)))
             ug756_34(iut)%c(2) = u756_34(iut)%c(2)
     &           /((fcl(id8))*(fcl(id3)))
             ug756_34(iut)%d(1) = u756_34(iut)%d(1)
     &           /((fcl(id8))*(fcl(id3)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p756,qd=p812,v=cz34(i3)%e,a=u756_34(i3)%a,b=u756_34(i3)%b,c=
* u756_34(i3)%c,d=u756_34(i3)%d,cl=zcl(id8),nsum=1
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
      u756_34(i3)%a(2)=u756_34(i3)%a(2)+zcl(id8)*(cauxa-ceps_0)
      u756_34(i3)%b(1)=u756_34(i3)%b(1)+zcl(id8)*(cauxb-ceps_2)
      u756_34(i3)%c(2)=u756_34(i3)%c(2)+zcl(id8)*(-cauxc+ceps_1)
      u756_34(i3)%d(1)=u756_34(i3)%d(1)+zcl(id8)*cz34(i3)%ek0
      end do
      else
* quqd -- p=p756,q=p812
      quqd=p756(0)*p812(0)-p756(1)*p812(1)-p756(2)*p812(2)-p756(
     & 3)*p812(3)
      do i3=1,2
* TW0 -- qu=p756,qd=p812,v=cz34(i3)%e,a=u756_34(i3)%a,b=u756_34(i3)%b,c=
* u756_34(i3)%c,d=u756_34(i3)%d,cl=zcl(id8),nsum=0
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
      u756_34(i3)%a(2)=zcl(id8)*(cauxa-ceps_0)
      u756_34(i3)%b(1)=zcl(id8)*(cauxb-ceps_2)
      u756_34(i3)%c(2)=zcl(id8)*(-cauxc+ceps_1)
      u756_34(i3)%d(1)=zcl(id8)*cz34(i3)%ek0
      end do
      endif
  
      if (ineutri(id7).ne.1.and.ineutri(id3).ne.1) then
* quqd -- p=p712,q=p856
      quqd=p712(0)*p856(0)-p712(1)*p856(1)-p712(2)*p856(2)-p712(
     & 3)*p856(3)
      do i3=1,2
* TW0 -- qu=p712,qd=p856,v=cf34(i3)%e,a=u712_34(i3)%a,b=u712_34(i3)%b,c=
* u712_34(i3)%c,d=u712_34(i3)%d,cl=fcl(id7),nsum=0
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
      u712_34(i3)%a(2)=fcl(id7)*(cauxa-ceps_0)
      u712_34(i3)%b(1)=fcl(id7)*(cauxb-ceps_2)
      u712_34(i3)%c(2)=fcl(id7)*(-cauxc+ceps_1)
      u712_34(i3)%d(1)=fcl(id7)*cf34(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id7).ne.1.and.ilept(id3).ne.1) then
           do iut=1,2
             ug712_34(iut)%a(2) = u712_34(iut)%a(2)
     &           /((fcl(id7))*(fcl(id3)))
             ug712_34(iut)%b(1) = u712_34(iut)%b(1)
     &           /((fcl(id7))*(fcl(id3)))
             ug712_34(iut)%c(2) = u712_34(iut)%c(2)
     &           /((fcl(id7))*(fcl(id3)))
             ug712_34(iut)%d(1) = u712_34(iut)%d(1)
     &           /((fcl(id7))*(fcl(id3)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p712,qd=p856,v=cz34(i3)%e,a=u712_34(i3)%a,b=u712_34(i3)%b,c=
* u712_34(i3)%c,d=u712_34(i3)%d,cl=zcl(id7),nsum=1
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
      u712_34(i3)%a(2)=u712_34(i3)%a(2)+zcl(id7)*(cauxa-ceps_0)
      u712_34(i3)%b(1)=u712_34(i3)%b(1)+zcl(id7)*(cauxb-ceps_2)
      u712_34(i3)%c(2)=u712_34(i3)%c(2)+zcl(id7)*(-cauxc+ceps_1)
      u712_34(i3)%d(1)=u712_34(i3)%d(1)+zcl(id7)*cz34(i3)%ek0
      end do
      else
* quqd -- p=p712,q=p856
      quqd=p712(0)*p856(0)-p712(1)*p856(1)-p712(2)*p856(2)-p712(
     & 3)*p856(3)
      do i3=1,2
* TW0 -- qu=p712,qd=p856,v=cz34(i3)%e,a=u712_34(i3)%a,b=u712_34(i3)%b,c=
* u712_34(i3)%c,d=u712_34(i3)%d,cl=zcl(id7),nsum=0
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
      u712_34(i3)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u712_34(i3)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u712_34(i3)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u712_34(i3)%d(1)=zcl(id7)*cz34(i3)%ek0
      end do
      endif
  
*  W insertion                                                          
* quqd -- p=p712,q=p834
      quqd=p712(0)*p834(0)-p712(1)*p834(1)-p712(2)*p834(2)-p712(
     & 3)*p834(3)
* TW0 -- qu=p712,qd=p834,v=cw56%e,a=u712_56%a,b=u712_56%b,c=u712_56%c,d=
* u712_56%d,cl=wcl,nsum=0
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
      u712_56%a(2)=wcl*(cauxa-ceps_0)
      u712_56%b(1)=wcl*(cauxb-ceps_2)
      u712_56%c(2)=wcl*(-cauxc+ceps_1)
      u712_56%d(1)=wcl*cw56%ek0
  
*  Z insertion                                                          
      if (ineutri(id8).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p756,q=p834
      quqd=p756(0)*p834(0)-p756(1)*p834(1)-p756(2)*p834(2)-p756(
     & 3)*p834(3)
      do i3=1,2
* TW0 -- qu=p756,qd=p834,v=cf12(i3)%e,a=u756_12(i3)%a,b=u756_12(i3)%b,c=
* u756_12(i3)%c,d=u756_12(i3)%d,cl=fcl(id8),nsum=0
      ceps_0=-cf12(i3)%ek0*(p756(2)*p834(3)-p834(2)*p756(3))+p75
     & 6k0*(cf12(i3)%e(2)*p834(3)-p834(2)*cf12(i3)%e(3))-p834k0*
     & (cf12(i3)%e(2)*p756(3)-p756(2)*cf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i3)%e(3)*p756k0+p756(3)*cf12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i3)%e(3)*p834k0+p834(3)*cf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i3)%e(0)*p756(0)-cf12(i3)%e(1)*p756(1)-cf12(i3)%
     & e(2)*p756(2)-cf12(i3)%e(3)*p756(3)
      cvqd=cf12(i3)%e(0)*p834(0)-cf12(i3)%e(1)*p834(1)-cf12(i3)%
     & e(2)*p834(2)-cf12(i3)%e(3)*p834(3)
      cauxa=-cf12(i3)%ek0*quqd+p756k0*cvqd+p834k0*cvqu
      cauxb=-cf12(i3)%ek0*p834(2)+p834k0*cf12(i3)%e(2)
      cauxc=+cf12(i3)%ek0*p756(2)-p756k0*cf12(i3)%e(2)
      u756_12(i3)%a(2)=fcl(id8)*(cauxa-ceps_0)
      u756_12(i3)%b(1)=fcl(id8)*(cauxb-ceps_2)
      u756_12(i3)%c(2)=fcl(id8)*(-cauxc+ceps_1)
      u756_12(i3)%d(1)=fcl(id8)*cf12(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id8).ne.1.and.ilept(id1).ne.1) then
           do iut=1,2
             ug756_12(iut)%a(2) = u756_12(iut)%a(2)
     &           /((fcl(id8))*(fcl(id1)))
             ug756_12(iut)%b(1) = u756_12(iut)%b(1)
     &           /((fcl(id8))*(fcl(id1)))
             ug756_12(iut)%c(2) = u756_12(iut)%c(2)
     &           /((fcl(id8))*(fcl(id1)))
             ug756_12(iut)%d(1) = u756_12(iut)%d(1)
     &           /((fcl(id8))*(fcl(id1)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p756,qd=p834,v=cz12(i3)%e,a=u756_12(i3)%a,b=u756_12(i3)%b,c=
* u756_12(i3)%c,d=u756_12(i3)%d,cl=zcl(id8),nsum=1
      ceps_0=-cz12(i3)%ek0*(p756(2)*p834(3)-p834(2)*p756(3))+p75
     & 6k0*(cz12(i3)%e(2)*p834(3)-p834(2)*cz12(i3)%e(3))-p834k0*
     & (cz12(i3)%e(2)*p756(3)-p756(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p756k0+p756(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p834k0+p834(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p756(0)-cz12(i3)%e(1)*p756(1)-cz12(i3)%
     & e(2)*p756(2)-cz12(i3)%e(3)*p756(3)
      cvqd=cz12(i3)%e(0)*p834(0)-cz12(i3)%e(1)*p834(1)-cz12(i3)%
     & e(2)*p834(2)-cz12(i3)%e(3)*p834(3)
      cauxa=-cz12(i3)%ek0*quqd+p756k0*cvqd+p834k0*cvqu
      cauxb=-cz12(i3)%ek0*p834(2)+p834k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p756(2)-p756k0*cz12(i3)%e(2)
      u756_12(i3)%a(2)=u756_12(i3)%a(2)+zcl(id8)*(cauxa-ceps_0)
      u756_12(i3)%b(1)=u756_12(i3)%b(1)+zcl(id8)*(cauxb-ceps_2)
      u756_12(i3)%c(2)=u756_12(i3)%c(2)+zcl(id8)*(-cauxc+ceps_1)
      u756_12(i3)%d(1)=u756_12(i3)%d(1)+zcl(id8)*cz12(i3)%ek0
      end do
      else
* quqd -- p=p756,q=p834
      quqd=p756(0)*p834(0)-p756(1)*p834(1)-p756(2)*p834(2)-p756(
     & 3)*p834(3)
      do i3=1,2
* TW0 -- qu=p756,qd=p834,v=cz12(i3)%e,a=u756_12(i3)%a,b=u756_12(i3)%b,c=
* u756_12(i3)%c,d=u756_12(i3)%d,cl=zcl(id8),nsum=0
      ceps_0=-cz12(i3)%ek0*(p756(2)*p834(3)-p834(2)*p756(3))+p75
     & 6k0*(cz12(i3)%e(2)*p834(3)-p834(2)*cz12(i3)%e(3))-p834k0*
     & (cz12(i3)%e(2)*p756(3)-p756(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p756k0+p756(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p834k0+p834(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p756(0)-cz12(i3)%e(1)*p756(1)-cz12(i3)%
     & e(2)*p756(2)-cz12(i3)%e(3)*p756(3)
      cvqd=cz12(i3)%e(0)*p834(0)-cz12(i3)%e(1)*p834(1)-cz12(i3)%
     & e(2)*p834(2)-cz12(i3)%e(3)*p834(3)
      cauxa=-cz12(i3)%ek0*quqd+p756k0*cvqd+p834k0*cvqu
      cauxb=-cz12(i3)%ek0*p834(2)+p834k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p756(2)-p756k0*cz12(i3)%e(2)
      u756_12(i3)%a(2)=zcl(id8)*(cauxa-ceps_0)
      u756_12(i3)%b(1)=zcl(id8)*(cauxb-ceps_2)
      u756_12(i3)%c(2)=zcl(id8)*(-cauxc+ceps_1)
      u756_12(i3)%d(1)=zcl(id8)*cz12(i3)%ek0
      end do
      endif
  
      if (ineutri(id7).ne.1.and.ineutri(id1).ne.1) then
* quqd -- p=p734,q=p856
      quqd=p734(0)*p856(0)-p734(1)*p856(1)-p734(2)*p856(2)-p734(
     & 3)*p856(3)
      do i3=1,2
* TW0 -- qu=p734,qd=p856,v=cf12(i3)%e,a=u734_12(i3)%a,b=u734_12(i3)%b,c=
* u734_12(i3)%c,d=u734_12(i3)%d,cl=fcl(id7),nsum=0
      ceps_0=-cf12(i3)%ek0*(p734(2)*p856(3)-p856(2)*p734(3))+p73
     & 4k0*(cf12(i3)%e(2)*p856(3)-p856(2)*cf12(i3)%e(3))-p856k0*
     & (cf12(i3)%e(2)*p734(3)-p734(2)*cf12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i3)%e(3)*p734k0+p734(3)*cf12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i3)%e(3)*p856k0+p856(3)*cf12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i3)%e(0)*p734(0)-cf12(i3)%e(1)*p734(1)-cf12(i3)%
     & e(2)*p734(2)-cf12(i3)%e(3)*p734(3)
      cvqd=cf12(i3)%e(0)*p856(0)-cf12(i3)%e(1)*p856(1)-cf12(i3)%
     & e(2)*p856(2)-cf12(i3)%e(3)*p856(3)
      cauxa=-cf12(i3)%ek0*quqd+p734k0*cvqd+p856k0*cvqu
      cauxb=-cf12(i3)%ek0*p856(2)+p856k0*cf12(i3)%e(2)
      cauxc=+cf12(i3)%ek0*p734(2)-p734k0*cf12(i3)%e(2)
      u734_12(i3)%a(2)=fcl(id7)*(cauxa-ceps_0)
      u734_12(i3)%b(1)=fcl(id7)*(cauxb-ceps_2)
      u734_12(i3)%c(2)=fcl(id7)*(-cauxc+ceps_1)
      u734_12(i3)%d(1)=fcl(id7)*cf12(i3)%ek0
      end do
**QCD                                                                   
          if (ilept(id7).ne.1.and.ilept(id1).ne.1) then
           do iut=1,2
             ug734_12(iut)%a(2) = u734_12(iut)%a(2)
     &           /((fcl(id7))*(fcl(id1)))
             ug734_12(iut)%b(1) = u734_12(iut)%b(1)
     &           /((fcl(id7))*(fcl(id1)))
             ug734_12(iut)%c(2) = u734_12(iut)%c(2)
     &           /((fcl(id7))*(fcl(id1)))
             ug734_12(iut)%d(1) = u734_12(iut)%d(1)
     &           /((fcl(id7))*(fcl(id1)))
           enddo
          endif
      do i3=1,2
* TW0 -- qu=p734,qd=p856,v=cz12(i3)%e,a=u734_12(i3)%a,b=u734_12(i3)%b,c=
* u734_12(i3)%c,d=u734_12(i3)%d,cl=zcl(id7),nsum=1
      ceps_0=-cz12(i3)%ek0*(p734(2)*p856(3)-p856(2)*p734(3))+p73
     & 4k0*(cz12(i3)%e(2)*p856(3)-p856(2)*cz12(i3)%e(3))-p856k0*
     & (cz12(i3)%e(2)*p734(3)-p734(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p734k0+p734(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p856k0+p856(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p734(0)-cz12(i3)%e(1)*p734(1)-cz12(i3)%
     & e(2)*p734(2)-cz12(i3)%e(3)*p734(3)
      cvqd=cz12(i3)%e(0)*p856(0)-cz12(i3)%e(1)*p856(1)-cz12(i3)%
     & e(2)*p856(2)-cz12(i3)%e(3)*p856(3)
      cauxa=-cz12(i3)%ek0*quqd+p734k0*cvqd+p856k0*cvqu
      cauxb=-cz12(i3)%ek0*p856(2)+p856k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p734(2)-p734k0*cz12(i3)%e(2)
      u734_12(i3)%a(2)=u734_12(i3)%a(2)+zcl(id7)*(cauxa-ceps_0)
      u734_12(i3)%b(1)=u734_12(i3)%b(1)+zcl(id7)*(cauxb-ceps_2)
      u734_12(i3)%c(2)=u734_12(i3)%c(2)+zcl(id7)*(-cauxc+ceps_1)
      u734_12(i3)%d(1)=u734_12(i3)%d(1)+zcl(id7)*cz12(i3)%ek0
      end do
      else
* quqd -- p=p734,q=p856
      quqd=p734(0)*p856(0)-p734(1)*p856(1)-p734(2)*p856(2)-p734(
     & 3)*p856(3)
      do i3=1,2
* TW0 -- qu=p734,qd=p856,v=cz12(i3)%e,a=u734_12(i3)%a,b=u734_12(i3)%b,c=
* u734_12(i3)%c,d=u734_12(i3)%d,cl=zcl(id7),nsum=0
      ceps_0=-cz12(i3)%ek0*(p734(2)*p856(3)-p856(2)*p734(3))+p73
     & 4k0*(cz12(i3)%e(2)*p856(3)-p856(2)*cz12(i3)%e(3))-p856k0*
     & (cz12(i3)%e(2)*p734(3)-p734(2)*cz12(i3)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i3)%e(3)*p734k0+p734(3)*cz12(i3)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i3)%e(3)*p856k0+p856(3)*cz12(i3)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i3)%e(0)*p734(0)-cz12(i3)%e(1)*p734(1)-cz12(i3)%
     & e(2)*p734(2)-cz12(i3)%e(3)*p734(3)
      cvqd=cz12(i3)%e(0)*p856(0)-cz12(i3)%e(1)*p856(1)-cz12(i3)%
     & e(2)*p856(2)-cz12(i3)%e(3)*p856(3)
      cauxa=-cz12(i3)%ek0*quqd+p734k0*cvqd+p856k0*cvqu
      cauxb=-cz12(i3)%ek0*p856(2)+p856k0*cz12(i3)%e(2)
      cauxc=+cz12(i3)%ek0*p734(2)-p734k0*cz12(i3)%e(2)
      u734_12(i3)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u734_12(i3)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u734_12(i3)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u734_12(i3)%d(1)=zcl(id7)*cz12(i3)%ek0
      end do
      endif
  
*  W insertion                                                          
* quqd -- p=p734,q=p812
      quqd=p734(0)*p812(0)-p734(1)*p812(1)-p734(2)*p812(2)-p734(
     & 3)*p812(3)
* TW0 -- qu=p734,qd=p812,v=cw56%e,a=u734_56%a,b=u734_56%b,c=u734_56%c,d=
* u734_56%d,cl=wcl,nsum=0
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
      u734_56%a(2)=wcl*(cauxa-ceps_0)
      u734_56%b(1)=wcl*(cauxb-ceps_2)
      u734_56%c(2)=wcl*(-cauxc+ceps_1)
      u734_56%d(1)=wcl*cw56%ek0
  
  
*                                                                       
* LEFT*MIDDLE                                                           
*                                                                       
  
**EW**                                                                  
  
* main Z line                                                           
      if (iup(id1).eq.1) then
* TLT0_W -- aa=l1_7856%a,cc=l1_7856%c,a1=l1_78%a,c1=l1_78%c,a2=u178_56%a
* ,b2=u178_56%b,c2=u178_56%c,d2=u178_56%d,prq=p178q,nsum=0
      l1_7856%c(2)=l1_78%c(2)*p178q*u178_56%d(1)+l1_78%a(2)*u178
     & _56%c(2)
      l1_7856%a(2)=l1_78%c(2)*p178q*u178_56%b(1)+l1_78%a(2)*u178
     & _56%a(2)
  
      do i3=1,2
* TLT0_W -- aa=l1_3478(i3)%a,cc=l1_3478(i3)%c,a1=l1_78%a,c1=l1_78%c,a2=u
* 178_34(i3)%a,b2=u178_34(i3)%b,c2=u178_34(i3)%c,d2=u178_34(i3)%d,prq=p1
* 78q,nsum=0
      l1_3478(i3)%c(2)=l1_78%c(2)*p178q*u178_34(i3)%d(1)+l1_78%a
     & (2)*u178_34(i3)%c(2)
      l1_3478(i3)%a(2)=l1_78%c(2)*p178q*u178_34(i3)%b(1)+l1_78%a
     & (2)*u178_34(i3)%a(2)
      end do
  
      do i3=1,2
* TLT0_W -- aa=l1_3478(i3)%a,cc=l1_3478(i3)%c,a1=l1_34(i3)%a,c1=l1_34(i3
* )%c,a2=u134_78%a,b2=u134_78%b,c2=u134_78%c,d2=u134_78%d,prq=p134q,nsum
* =1
      l1_3478(i3)%c(2)=l1_3478(i3)%c(2)+l1_34(i3)%c(2)*p134q*u13
     & 4_78%d(1)+l1_34(i3)%a(2)*u134_78%c(2)
      l1_3478(i3)%a(2)=l1_3478(i3)%a(2)+l1_34(i3)%c(2)*p134q*u13
     & 4_78%b(1)+l1_34(i3)%a(2)*u134_78%a(2)
      end do
      else
* TLT0_W -- aa=l1_5678%a,cc=l1_5678%c,a1=l1_56%a,c1=l1_56%c,a2=u156_78%a
* ,b2=u156_78%b,c2=u156_78%c,d2=u156_78%d,prq=p156q,nsum=0
      l1_5678%c(2)=l1_56%c(2)*p156q*u156_78%d(1)+l1_56%a(2)*u156
     & _78%c(2)
      l1_5678%a(2)=l1_56%c(2)*p156q*u156_78%b(1)+l1_56%a(2)*u156
     & _78%a(2)
  
      do i3=1,2
* TLT0_W -- aa=l1_3456(i3)%a,cc=l1_3456(i3)%c,a1=l1_56%a,c1=l1_56%c,a2=u
* 156_34(i3)%a,b2=u156_34(i3)%b,c2=u156_34(i3)%c,d2=u156_34(i3)%d,prq=p1
* 56q,nsum=0
      l1_3456(i3)%c(2)=l1_56%c(2)*p156q*u156_34(i3)%d(1)+l1_56%a
     & (2)*u156_34(i3)%c(2)
      l1_3456(i3)%a(2)=l1_56%c(2)*p156q*u156_34(i3)%b(1)+l1_56%a
     & (2)*u156_34(i3)%a(2)
      end do
  
      do i3=1,2
* TLT0_W -- aa=l1_3456(i3)%a,cc=l1_3456(i3)%c,a1=l1_34(i3)%a,c1=l1_34(i3
* )%c,a2=u134_56%a,b2=u134_56%b,c2=u134_56%c,d2=u134_56%d,prq=p134q,nsum
* =1
      l1_3456(i3)%c(2)=l1_3456(i3)%c(2)+l1_34(i3)%c(2)*p134q*u13
     & 4_56%d(1)+l1_34(i3)%a(2)*u134_56%c(2)
      l1_3456(i3)%a(2)=l1_3456(i3)%a(2)+l1_34(i3)%c(2)*p134q*u13
     & 4_56%b(1)+l1_34(i3)%a(2)*u134_56%a(2)
      end do
      endif
  
      if (iup(id3).eq.1) then
* TLT0_W -- aa=l3_7856%a,cc=l3_7856%c,a1=l3_78%a,c1=l3_78%c,a2=u378_56%a
* ,b2=u378_56%b,c2=u378_56%c,d2=u378_56%d,prq=p378q,nsum=0
      l3_7856%c(2)=l3_78%c(2)*p378q*u378_56%d(1)+l3_78%a(2)*u378
     & _56%c(2)
      l3_7856%a(2)=l3_78%c(2)*p378q*u378_56%b(1)+l3_78%a(2)*u378
     & _56%a(2)
  
      do i3=1,2
* TLT0_W -- aa=l3_1278(i3)%a,cc=l3_1278(i3)%c,a1=l3_78%a,c1=l3_78%c,a2=u
* 378_12(i3)%a,b2=u378_12(i3)%b,c2=u378_12(i3)%c,d2=u378_12(i3)%d,prq=p3
* 78q,nsum=0
      l3_1278(i3)%c(2)=l3_78%c(2)*p378q*u378_12(i3)%d(1)+l3_78%a
     & (2)*u378_12(i3)%c(2)
      l3_1278(i3)%a(2)=l3_78%c(2)*p378q*u378_12(i3)%b(1)+l3_78%a
     & (2)*u378_12(i3)%a(2)
      end do
  
      do i3=1,2
* TLT0_W -- aa=l3_1278(i3)%a,cc=l3_1278(i3)%c,a1=l3_12(i3)%a,c1=l3_12(i3
* )%c,a2=u312_78%a,b2=u312_78%b,c2=u312_78%c,d2=u312_78%d,prq=p312q,nsum
* =1
      l3_1278(i3)%c(2)=l3_1278(i3)%c(2)+l3_12(i3)%c(2)*p312q*u31
     & 2_78%d(1)+l3_12(i3)%a(2)*u312_78%c(2)
      l3_1278(i3)%a(2)=l3_1278(i3)%a(2)+l3_12(i3)%c(2)*p312q*u31
     & 2_78%b(1)+l3_12(i3)%a(2)*u312_78%a(2)
      end do
      else
* TLT0_W -- aa=l3_5678%a,cc=l3_5678%c,a1=l3_56%a,c1=l3_56%c,a2=u356_78%a
* ,b2=u356_78%b,c2=u356_78%c,d2=u356_78%d,prq=p356q,nsum=0
      l3_5678%c(2)=l3_56%c(2)*p356q*u356_78%d(1)+l3_56%a(2)*u356
     & _78%c(2)
      l3_5678%a(2)=l3_56%c(2)*p356q*u356_78%b(1)+l3_56%a(2)*u356
     & _78%a(2)
  
      do i3=1,2
* TLT0_W -- aa=l3_1256(i3)%a,cc=l3_1256(i3)%c,a1=l3_56%a,c1=l3_56%c,a2=u
* 356_12(i3)%a,b2=u356_12(i3)%b,c2=u356_12(i3)%c,d2=u356_12(i3)%d,prq=p3
* 56q,nsum=0
      l3_1256(i3)%c(2)=l3_56%c(2)*p356q*u356_12(i3)%d(1)+l3_56%a
     & (2)*u356_12(i3)%c(2)
      l3_1256(i3)%a(2)=l3_56%c(2)*p356q*u356_12(i3)%b(1)+l3_56%a
     & (2)*u356_12(i3)%a(2)
      end do
  
      do i3=1,2
* TLT0_W -- aa=l3_1256(i3)%a,cc=l3_1256(i3)%c,a1=l3_12(i3)%a,c1=l3_12(i3
* )%c,a2=u312_56%a,b2=u312_56%b,c2=u312_56%c,d2=u312_56%d,prq=p312q,nsum
* =1
      l3_1256(i3)%c(2)=l3_1256(i3)%c(2)+l3_12(i3)%c(2)*p312q*u31
     & 2_56%d(1)+l3_12(i3)%a(2)*u312_56%c(2)
      l3_1256(i3)%a(2)=l3_1256(i3)%a(2)+l3_12(i3)%c(2)*p312q*u31
     & 2_56%b(1)+l3_12(i3)%a(2)*u312_56%a(2)
      end do
      endif
  
  
* main W                                                                
  
      do i1=1,2
* TLT0_W -- aa=l5_1278(i1)%a,cc=l5_1278(i1)%c,a1=l5_78%a,c1=l5_78%c,a2=u
* 578_12(i1)%a,b2=u578_12(i1)%b,c2=u578_12(i1)%c,d2=u578_12(i1)%d,prq=p5
* 78q,nsum=0
      l5_1278(i1)%c(2)=l5_78%c(2)*p578q*u578_12(i1)%d(1)+l5_78%a
     & (2)*u578_12(i1)%c(2)
      l5_1278(i1)%a(2)=l5_78%c(2)*p578q*u578_12(i1)%b(1)+l5_78%a
     & (2)*u578_12(i1)%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=l5_1278(i1)%a,cc=l5_1278(i1)%c,a1=l5_12(i1)%a,c1=l5_12(i1
* )%c,a2=u512_78%a,b2=u512_78%b,c2=u512_78%c,d2=u512_78%d,prq=p512q,nsum
* =1
      l5_1278(i1)%c(2)=l5_1278(i1)%c(2)+l5_12(i1)%c(2)*p512q*u51
     & 2_78%d(1)+l5_12(i1)%a(2)*u512_78%c(2)
      l5_1278(i1)%a(2)=l5_1278(i1)%a(2)+l5_12(i1)%c(2)*p512q*u51
     & 2_78%b(1)+l5_12(i1)%a(2)*u512_78%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=l5_3478(i1)%a,cc=l5_3478(i1)%c,a1=l5_78%a,c1=l5_78%c,a2=u
* 578_34(i1)%a,b2=u578_34(i1)%b,c2=u578_34(i1)%c,d2=u578_34(i1)%d,prq=p5
* 78q,nsum=0
      l5_3478(i1)%c(2)=l5_78%c(2)*p578q*u578_34(i1)%d(1)+l5_78%a
     & (2)*u578_34(i1)%c(2)
      l5_3478(i1)%a(2)=l5_78%c(2)*p578q*u578_34(i1)%b(1)+l5_78%a
     & (2)*u578_34(i1)%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=l5_3478(i1)%a,cc=l5_3478(i1)%c,a1=l5_34(i1)%a,c1=l5_34(i1
* )%c,a2=u534_78%a,b2=u534_78%b,c2=u534_78%c,d2=u534_78%d,prq=p534q,nsum
* =1
      l5_3478(i1)%c(2)=l5_3478(i1)%c(2)+l5_34(i1)%c(2)*p534q*u53
     & 4_78%d(1)+l5_34(i1)%a(2)*u534_78%c(2)
      l5_3478(i1)%a(2)=l5_3478(i1)%a(2)+l5_34(i1)%c(2)*p534q*u53
     & 4_78%b(1)+l5_34(i1)%a(2)*u534_78%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=l7_1256(i1)%a,cc=l7_1256(i1)%c,a1=l7_56%a,c1=l7_56%c,a2=u
* 756_12(i1)%a,b2=u756_12(i1)%b,c2=u756_12(i1)%c,d2=u756_12(i1)%d,prq=p7
* 56q,nsum=0
      l7_1256(i1)%c(2)=l7_56%c(2)*p756q*u756_12(i1)%d(1)+l7_56%a
     & (2)*u756_12(i1)%c(2)
      l7_1256(i1)%a(2)=l7_56%c(2)*p756q*u756_12(i1)%b(1)+l7_56%a
     & (2)*u756_12(i1)%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=l7_1256(i1)%a,cc=l7_1256(i1)%c,a1=l7_12(i1)%a,c1=l7_12(i1
* )%c,a2=u712_56%a,b2=u712_56%b,c2=u712_56%c,d2=u712_56%d,prq=p712q,nsum
* =1
      l7_1256(i1)%c(2)=l7_1256(i1)%c(2)+l7_12(i1)%c(2)*p712q*u71
     & 2_56%d(1)+l7_12(i1)%a(2)*u712_56%c(2)
      l7_1256(i1)%a(2)=l7_1256(i1)%a(2)+l7_12(i1)%c(2)*p712q*u71
     & 2_56%b(1)+l7_12(i1)%a(2)*u712_56%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=l7_3456(i1)%a,cc=l7_3456(i1)%c,a1=l7_56%a,c1=l7_56%c,a2=u
* 756_34(i1)%a,b2=u756_34(i1)%b,c2=u756_34(i1)%c,d2=u756_34(i1)%d,prq=p7
* 56q,nsum=0
      l7_3456(i1)%c(2)=l7_56%c(2)*p756q*u756_34(i1)%d(1)+l7_56%a
     & (2)*u756_34(i1)%c(2)
      l7_3456(i1)%a(2)=l7_56%c(2)*p756q*u756_34(i1)%b(1)+l7_56%a
     & (2)*u756_34(i1)%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=l7_3456(i1)%a,cc=l7_3456(i1)%c,a1=l7_34(i1)%a,c1=l7_34(i1
* )%c,a2=u734_56%a,b2=u734_56%b,c2=u734_56%c,d2=u734_56%d,prq=p734q,nsum
* =1
      l7_3456(i1)%c(2)=l7_3456(i1)%c(2)+l7_34(i1)%c(2)*p734q*u73
     & 4_56%d(1)+l7_34(i1)%a(2)*u734_56%c(2)
      l7_3456(i1)%a(2)=l7_3456(i1)%a(2)+l7_34(i1)%c(2)*p734q*u73
     & 4_56%b(1)+l7_34(i1)%a(2)*u734_56%a(2)
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l5_1234(i1,i3)%a,cc=l5_1234(i1,i3)%c,a1=l5_12(i1)%a,c1=l5
* _12(i1)%c,a2=u512_34(i3)%a,b2=u512_34(i3)%b,c2=u512_34(i3)%c,d2=u512_3
* 4(i3)%d,prq=p512q,nsum=0
      l5_1234(i1,i3)%c(2)=l5_12(i1)%c(2)*p512q*u512_34(i3)%d(1)+
     & l5_12(i1)%a(2)*u512_34(i3)%c(2)
      l5_1234(i1,i3)%a(2)=l5_12(i1)%c(2)*p512q*u512_34(i3)%b(1)+
     & l5_12(i1)%a(2)*u512_34(i3)%a(2)
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l5_1234(i1,i3)%a,cc=l5_1234(i1,i3)%c,a1=l5_34(i3)%a,c1=l5
* _34(i3)%c,a2=u534_12(i1)%a,b2=u534_12(i1)%b,c2=u534_12(i1)%c,d2=u534_1
* 2(i1)%d,prq=p534q,nsum=1
      l5_1234(i1,i3)%c(2)=l5_1234(i1,i3)%c(2)+l5_34(i3)%c(2)*p53
     & 4q*u534_12(i1)%d(1)+l5_34(i3)%a(2)*u534_12(i1)%c(2)
      l5_1234(i1,i3)%a(2)=l5_1234(i1,i3)%a(2)+l5_34(i3)%c(2)*p53
     & 4q*u534_12(i1)%b(1)+l5_34(i3)%a(2)*u534_12(i1)%a(2)
      end do
      end do
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l7_1234(i1,i3)%a,cc=l7_1234(i1,i3)%c,a1=l7_12(i1)%a,c1=l7
* _12(i1)%c,a2=u712_34(i3)%a,b2=u712_34(i3)%b,c2=u712_34(i3)%c,d2=u712_3
* 4(i3)%d,prq=p712q,nsum=0
      l7_1234(i1,i3)%c(2)=l7_12(i1)%c(2)*p712q*u712_34(i3)%d(1)+
     & l7_12(i1)%a(2)*u712_34(i3)%c(2)
      l7_1234(i1,i3)%a(2)=l7_12(i1)%c(2)*p712q*u712_34(i3)%b(1)+
     & l7_12(i1)%a(2)*u712_34(i3)%a(2)
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=l7_1234(i1,i3)%a,cc=l7_1234(i1,i3)%c,a1=l7_34(i3)%a,c1=l7
* _34(i3)%c,a2=u734_12(i1)%a,b2=u734_12(i1)%b,c2=u734_12(i1)%c,d2=u734_1
* 2(i1)%d,prq=p734q,nsum=1
      l7_1234(i1,i3)%c(2)=l7_1234(i1,i3)%c(2)+l7_34(i3)%c(2)*p73
     & 4q*u734_12(i1)%d(1)+l7_34(i3)%a(2)*u734_12(i1)%c(2)
      l7_1234(i1,i3)%a(2)=l7_1234(i1,i3)%a(2)+l7_34(i3)%c(2)*p73
     & 4q*u734_12(i1)%b(1)+l7_34(i3)%a(2)*u734_12(i1)%a(2)
      end do
      end do
  
**QCD**                                                                 
* main Z line                                                           
  
      if (ilept(id1).ne.1.and.ilept(id3).ne.1) then
        if (iup(id1).eq.1) then
      do i3=1,2
* TLT0_W -- aa=lg1_3478(i3)%a,cc=lg1_3478(i3)%c,a1=l1_78%a,c1=l1_78%c,a2
* =ug178_34(i3)%a,b2=ug178_34(i3)%b,c2=ug178_34(i3)%c,d2=ug178_34(i3)%d,
* prq=p178q,nsum=0
      lg1_3478(i3)%c(2)=l1_78%c(2)*p178q*ug178_34(i3)%d(1)+l1_78
     & %a(2)*ug178_34(i3)%c(2)
      lg1_3478(i3)%a(2)=l1_78%c(2)*p178q*ug178_34(i3)%b(1)+l1_78
     & %a(2)*ug178_34(i3)%a(2)
      end do
  
      do i3=1,2
* TLT0_W -- aa=lg1_3478(i3)%a,cc=lg1_3478(i3)%c,a1=lg1_34(i3)%a,c1=lg1_3
* 4(i3)%c,a2=u134_78%a,b2=u134_78%b,c2=u134_78%c,d2=u134_78%d,prq=p134q,
* nsum=1
      lg1_3478(i3)%c(2)=lg1_3478(i3)%c(2)+lg1_34(i3)%c(2)*p134q*
     & u134_78%d(1)+lg1_34(i3)%a(2)*u134_78%c(2)
      lg1_3478(i3)%a(2)=lg1_3478(i3)%a(2)+lg1_34(i3)%c(2)*p134q*
     & u134_78%b(1)+lg1_34(i3)%a(2)*u134_78%a(2)
      end do
        else
      do i3=1,2
* TLT0_W -- aa=lg1_3456(i3)%a,cc=lg1_3456(i3)%c,a1=l1_56%a,c1=l1_56%c,a2
* =ug156_34(i3)%a,b2=ug156_34(i3)%b,c2=ug156_34(i3)%c,d2=ug156_34(i3)%d,
* prq=p156q,nsum=0
      lg1_3456(i3)%c(2)=l1_56%c(2)*p156q*ug156_34(i3)%d(1)+l1_56
     & %a(2)*ug156_34(i3)%c(2)
      lg1_3456(i3)%a(2)=l1_56%c(2)*p156q*ug156_34(i3)%b(1)+l1_56
     & %a(2)*ug156_34(i3)%a(2)
      end do
  
      do i3=1,2
* TLT0_W -- aa=lg1_3456(i3)%a,cc=lg1_3456(i3)%c,a1=lg1_34(i3)%a,c1=lg1_3
* 4(i3)%c,a2=u134_56%a,b2=u134_56%b,c2=u134_56%c,d2=u134_56%d,prq=p134q,
* nsum=1
      lg1_3456(i3)%c(2)=lg1_3456(i3)%c(2)+lg1_34(i3)%c(2)*p134q*
     & u134_56%d(1)+lg1_34(i3)%a(2)*u134_56%c(2)
      lg1_3456(i3)%a(2)=lg1_3456(i3)%a(2)+lg1_34(i3)%c(2)*p134q*
     & u134_56%b(1)+lg1_34(i3)%a(2)*u134_56%a(2)
      end do
        endif
      endif
  
      if (ilept(id3).ne.1.and.ilept(id1).ne.1) then
        if (iup(id3).eq.1) then
      do i3=1,2
* TLT0_W -- aa=lg3_1278(i3)%a,cc=lg3_1278(i3)%c,a1=l3_78%a,c1=l3_78%c,a2
* =ug378_12(i3)%a,b2=ug378_12(i3)%b,c2=ug378_12(i3)%c,d2=ug378_12(i3)%d,
* prq=p378q,nsum=0
      lg3_1278(i3)%c(2)=l3_78%c(2)*p378q*ug378_12(i3)%d(1)+l3_78
     & %a(2)*ug378_12(i3)%c(2)
      lg3_1278(i3)%a(2)=l3_78%c(2)*p378q*ug378_12(i3)%b(1)+l3_78
     & %a(2)*ug378_12(i3)%a(2)
      end do
  
      do i3=1,2
* TLT0_W -- aa=lg3_1278(i3)%a,cc=lg3_1278(i3)%c,a1=lg3_12(i3)%a,c1=lg3_1
* 2(i3)%c,a2=u312_78%a,b2=u312_78%b,c2=u312_78%c,d2=u312_78%d,prq=p312q,
* nsum=1
      lg3_1278(i3)%c(2)=lg3_1278(i3)%c(2)+lg3_12(i3)%c(2)*p312q*
     & u312_78%d(1)+lg3_12(i3)%a(2)*u312_78%c(2)
      lg3_1278(i3)%a(2)=lg3_1278(i3)%a(2)+lg3_12(i3)%c(2)*p312q*
     & u312_78%b(1)+lg3_12(i3)%a(2)*u312_78%a(2)
      end do
        else
      do i3=1,2
* TLT0_W -- aa=lg3_1256(i3)%a,cc=lg3_1256(i3)%c,a1=l3_56%a,c1=l3_56%c,a2
* =ug356_12(i3)%a,b2=ug356_12(i3)%b,c2=ug356_12(i3)%c,d2=ug356_12(i3)%d,
* prq=p356q,nsum=0
      lg3_1256(i3)%c(2)=l3_56%c(2)*p356q*ug356_12(i3)%d(1)+l3_56
     & %a(2)*ug356_12(i3)%c(2)
      lg3_1256(i3)%a(2)=l3_56%c(2)*p356q*ug356_12(i3)%b(1)+l3_56
     & %a(2)*ug356_12(i3)%a(2)
      end do
  
      do i3=1,2
* TLT0_W -- aa=lg3_1256(i3)%a,cc=lg3_1256(i3)%c,a1=lg3_12(i3)%a,c1=lg3_1
* 2(i3)%c,a2=u312_56%a,b2=u312_56%b,c2=u312_56%c,d2=u312_56%d,prq=p312q,
* nsum=1
      lg3_1256(i3)%c(2)=lg3_1256(i3)%c(2)+lg3_12(i3)%c(2)*p312q*
     & u312_56%d(1)+lg3_12(i3)%a(2)*u312_56%c(2)
      lg3_1256(i3)%a(2)=lg3_1256(i3)%a(2)+lg3_12(i3)%c(2)*p312q*
     & u312_56%b(1)+lg3_12(i3)%a(2)*u312_56%a(2)
      end do
        endif
      endif
  
  
**QCD main W                                                            
  
      do i1=1,2
* TLT0_W -- aa=lg5_1278(i1)%a,cc=lg5_1278(i1)%c,a1=l5_78%a,c1=l5_78%c,a2
* =ug578_12(i1)%a,b2=ug578_12(i1)%b,c2=ug578_12(i1)%c,d2=ug578_12(i1)%d,
* prq=p578q,nsum=0
      lg5_1278(i1)%c(2)=l5_78%c(2)*p578q*ug578_12(i1)%d(1)+l5_78
     & %a(2)*ug578_12(i1)%c(2)
      lg5_1278(i1)%a(2)=l5_78%c(2)*p578q*ug578_12(i1)%b(1)+l5_78
     & %a(2)*ug578_12(i1)%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=lg5_1278(i1)%a,cc=lg5_1278(i1)%c,a1=lg5_12(i1)%a,c1=lg5_1
* 2(i1)%c,a2=u512_78%a,b2=u512_78%b,c2=u512_78%c,d2=u512_78%d,prq=p512q,
* nsum=1
      lg5_1278(i1)%c(2)=lg5_1278(i1)%c(2)+lg5_12(i1)%c(2)*p512q*
     & u512_78%d(1)+lg5_12(i1)%a(2)*u512_78%c(2)
      lg5_1278(i1)%a(2)=lg5_1278(i1)%a(2)+lg5_12(i1)%c(2)*p512q*
     & u512_78%b(1)+lg5_12(i1)%a(2)*u512_78%a(2)
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=lg5_1234(i1,i3)%a,cc=lg5_1234(i1,i3)%c,a1=lg5_12(i1)%a,c1
* =lg5_12(i1)%c,a2=u512_34(i3)%a,b2=u512_34(i3)%b,c2=u512_34(i3)%c,d2=u5
* 12_34(i3)%d,prq=p512q,nsum=0
      lg5_1234(i1,i3)%c(2)=lg5_12(i1)%c(2)*p512q*u512_34(i3)%d(1
     & )+lg5_12(i1)%a(2)*u512_34(i3)%c(2)
      lg5_1234(i1,i3)%a(2)=lg5_12(i1)%c(2)*p512q*u512_34(i3)%b(1
     & )+lg5_12(i1)%a(2)*u512_34(i3)%a(2)
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=lg5_1234(i1,i3)%a,cc=lg5_1234(i1,i3)%c,a1=l5_34(i3)%a,c1=
* l5_34(i3)%c,a2=ug534_12(i1)%a,b2=ug534_12(i1)%b,c2=ug534_12(i1)%c,d2=u
* g534_12(i1)%d,prq=p534q,nsum=1
      lg5_1234(i1,i3)%c(2)=lg5_1234(i1,i3)%c(2)+l5_34(i3)%c(2)*p
     & 534q*ug534_12(i1)%d(1)+l5_34(i3)%a(2)*ug534_12(i1)%c(2)
      lg5_1234(i1,i3)%a(2)=lg5_1234(i1,i3)%a(2)+l5_34(i3)%c(2)*p
     & 534q*ug534_12(i1)%b(1)+l5_34(i3)%a(2)*ug534_12(i1)%a(2)
      end do
      end do
  
      do i1=1,2
* TLT0_W -- aa=lg5_3478(i1)%a,cc=lg5_3478(i1)%c,a1=l5_78%a,c1=l5_78%c,a2
* =ug578_34(i1)%a,b2=ug578_34(i1)%b,c2=ug578_34(i1)%c,d2=ug578_34(i1)%d,
* prq=p578q,nsum=0
      lg5_3478(i1)%c(2)=l5_78%c(2)*p578q*ug578_34(i1)%d(1)+l5_78
     & %a(2)*ug578_34(i1)%c(2)
      lg5_3478(i1)%a(2)=l5_78%c(2)*p578q*ug578_34(i1)%b(1)+l5_78
     & %a(2)*ug578_34(i1)%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=lg5_3478(i1)%a,cc=lg5_3478(i1)%c,a1=lg5_34(i1)%a,c1=lg5_3
* 4(i1)%c,a2=u534_78%a,b2=u534_78%b,c2=u534_78%c,d2=u534_78%d,prq=p534q,
* nsum=1
      lg5_3478(i1)%c(2)=lg5_3478(i1)%c(2)+lg5_34(i1)%c(2)*p534q*
     & u534_78%d(1)+lg5_34(i1)%a(2)*u534_78%c(2)
      lg5_3478(i1)%a(2)=lg5_3478(i1)%a(2)+lg5_34(i1)%c(2)*p534q*
     & u534_78%b(1)+lg5_34(i1)%a(2)*u534_78%a(2)
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=lg5_3412(i1,i3)%a,cc=lg5_3412(i1,i3)%c,a1=lg5_34(i1)%a,c1
* =lg5_34(i1)%c,a2=u534_12(i3)%a,b2=u534_12(i3)%b,c2=u534_12(i3)%c,d2=u5
* 34_12(i3)%d,prq=p534q,nsum=0
      lg5_3412(i1,i3)%c(2)=lg5_34(i1)%c(2)*p534q*u534_12(i3)%d(1
     & )+lg5_34(i1)%a(2)*u534_12(i3)%c(2)
      lg5_3412(i1,i3)%a(2)=lg5_34(i1)%c(2)*p534q*u534_12(i3)%b(1
     & )+lg5_34(i1)%a(2)*u534_12(i3)%a(2)
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=lg5_3412(i1,i3)%a,cc=lg5_3412(i1,i3)%c,a1=l5_12(i3)%a,c1=
* l5_12(i3)%c,a2=ug512_34(i1)%a,b2=ug512_34(i1)%b,c2=ug512_34(i1)%c,d2=u
* g512_34(i1)%d,prq=p512q,nsum=1
      lg5_3412(i1,i3)%c(2)=lg5_3412(i1,i3)%c(2)+l5_12(i3)%c(2)*p
     & 512q*ug512_34(i1)%d(1)+l5_12(i3)%a(2)*ug512_34(i1)%c(2)
      lg5_3412(i1,i3)%a(2)=lg5_3412(i1,i3)%a(2)+l5_12(i3)%c(2)*p
     & 512q*ug512_34(i1)%b(1)+l5_12(i3)%a(2)*ug512_34(i1)%a(2)
      end do
      end do
  
      do i1=1,2
* TLT0_W -- aa=lg7_1256(i1)%a,cc=lg7_1256(i1)%c,a1=l7_56%a,c1=l7_56%c,a2
* =ug756_12(i1)%a,b2=ug756_12(i1)%b,c2=ug756_12(i1)%c,d2=ug756_12(i1)%d,
* prq=p756q,nsum=0
      lg7_1256(i1)%c(2)=l7_56%c(2)*p756q*ug756_12(i1)%d(1)+l7_56
     & %a(2)*ug756_12(i1)%c(2)
      lg7_1256(i1)%a(2)=l7_56%c(2)*p756q*ug756_12(i1)%b(1)+l7_56
     & %a(2)*ug756_12(i1)%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=lg7_1256(i1)%a,cc=lg7_1256(i1)%c,a1=lg7_12(i1)%a,c1=lg7_1
* 2(i1)%c,a2=u712_56%a,b2=u712_56%b,c2=u712_56%c,d2=u712_56%d,prq=p712q,
* nsum=1
      lg7_1256(i1)%c(2)=lg7_1256(i1)%c(2)+lg7_12(i1)%c(2)*p712q*
     & u712_56%d(1)+lg7_12(i1)%a(2)*u712_56%c(2)
      lg7_1256(i1)%a(2)=lg7_1256(i1)%a(2)+lg7_12(i1)%c(2)*p712q*
     & u712_56%b(1)+lg7_12(i1)%a(2)*u712_56%a(2)
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=lg7_1234(i1,i3)%a,cc=lg7_1234(i1,i3)%c,a1=lg7_12(i1)%a,c1
* =lg7_12(i1)%c,a2=u712_34(i3)%a,b2=u712_34(i3)%b,c2=u712_34(i3)%c,d2=u7
* 12_34(i3)%d,prq=p712q,nsum=0
      lg7_1234(i1,i3)%c(2)=lg7_12(i1)%c(2)*p712q*u712_34(i3)%d(1
     & )+lg7_12(i1)%a(2)*u712_34(i3)%c(2)
      lg7_1234(i1,i3)%a(2)=lg7_12(i1)%c(2)*p712q*u712_34(i3)%b(1
     & )+lg7_12(i1)%a(2)*u712_34(i3)%a(2)
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=lg7_1234(i1,i3)%a,cc=lg7_1234(i1,i3)%c,a1=l7_34(i3)%a,c1=
* l7_34(i3)%c,a2=ug734_12(i1)%a,b2=ug734_12(i1)%b,c2=ug734_12(i1)%c,d2=u
* g734_12(i1)%d,prq=p734q,nsum=1
      lg7_1234(i1,i3)%c(2)=lg7_1234(i1,i3)%c(2)+l7_34(i3)%c(2)*p
     & 734q*ug734_12(i1)%d(1)+l7_34(i3)%a(2)*ug734_12(i1)%c(2)
      lg7_1234(i1,i3)%a(2)=lg7_1234(i1,i3)%a(2)+l7_34(i3)%c(2)*p
     & 734q*ug734_12(i1)%b(1)+l7_34(i3)%a(2)*ug734_12(i1)%a(2)
      end do
      end do
  
      do i1=1,2
* TLT0_W -- aa=lg7_3456(i1)%a,cc=lg7_3456(i1)%c,a1=l7_56%a,c1=l7_56%c,a2
* =ug756_34(i1)%a,b2=ug756_34(i1)%b,c2=ug756_34(i1)%c,d2=ug756_34(i1)%d,
* prq=p756q,nsum=0
      lg7_3456(i1)%c(2)=l7_56%c(2)*p756q*ug756_34(i1)%d(1)+l7_56
     & %a(2)*ug756_34(i1)%c(2)
      lg7_3456(i1)%a(2)=l7_56%c(2)*p756q*ug756_34(i1)%b(1)+l7_56
     & %a(2)*ug756_34(i1)%a(2)
      end do
  
      do i1=1,2
* TLT0_W -- aa=lg7_3456(i1)%a,cc=lg7_3456(i1)%c,a1=lg7_34(i1)%a,c1=lg7_3
* 4(i1)%c,a2=u734_56%a,b2=u734_56%b,c2=u734_56%c,d2=u734_56%d,prq=p734q,
* nsum=1
      lg7_3456(i1)%c(2)=lg7_3456(i1)%c(2)+lg7_34(i1)%c(2)*p734q*
     & u734_56%d(1)+lg7_34(i1)%a(2)*u734_56%c(2)
      lg7_3456(i1)%a(2)=lg7_3456(i1)%a(2)+lg7_34(i1)%c(2)*p734q*
     & u734_56%b(1)+lg7_34(i1)%a(2)*u734_56%a(2)
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=lg7_3412(i1,i3)%a,cc=lg7_3412(i1,i3)%c,a1=lg7_34(i1)%a,c1
* =lg7_34(i1)%c,a2=u734_12(i3)%a,b2=u734_12(i3)%b,c2=u734_12(i3)%c,d2=u7
* 34_12(i3)%d,prq=p734q,nsum=0
      lg7_3412(i1,i3)%c(2)=lg7_34(i1)%c(2)*p734q*u734_12(i3)%d(1
     & )+lg7_34(i1)%a(2)*u734_12(i3)%c(2)
      lg7_3412(i1,i3)%a(2)=lg7_34(i1)%c(2)*p734q*u734_12(i3)%b(1
     & )+lg7_34(i1)%a(2)*u734_12(i3)%a(2)
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLT0_W -- aa=lg7_3412(i1,i3)%a,cc=lg7_3412(i1,i3)%c,a1=l7_12(i3)%a,c1=
* l7_12(i3)%c,a2=ug712_34(i1)%a,b2=ug712_34(i1)%b,c2=ug712_34(i1)%c,d2=u
* g712_34(i1)%d,prq=p712q,nsum=1
      lg7_3412(i1,i3)%c(2)=lg7_3412(i1,i3)%c(2)+l7_12(i3)%c(2)*p
     & 712q*ug712_34(i1)%d(1)+l7_12(i3)%a(2)*ug712_34(i1)%c(2)
      lg7_3412(i1,i3)%a(2)=lg7_3412(i1,i3)%a(2)+l7_12(i3)%c(2)*p
     & 712q*ug712_34(i1)%b(1)+l7_12(i3)%a(2)*ug712_34(i1)%a(2)
      end do
      end do
  
*                                                                       
* LEFT*MIDDLE*RIGHT (3forks)                                            
*                                                                       
  
**EW**                                                                  
* main line Z                                                           
      if (iup(id1).eq.1) then
      do i3=1,2
* TLTR0_W -- aa=c3fk12_345678(i3),a1=l1_7856%a,c1=l1_7856%c,a2=r2_34(i3)
* %a,b2=r2_34(i3)%b,prq=p234q,bef=,aft=
      c3fk12_345678(i3)=(l1_7856%c(2)*p234q*r2_34(i3)%b(1)+l1_78
     & 56%a(2)*r2_34(i3)%a(2))
      end do
  
      do i3=1,2
* TLTR0_W -- aa=c3fk12_345678(i3),a1=l1_3478(i3)%a,c1=l1_3478(i3)%c,a2=r
* 2_56%a,b2=r2_56%b,prq=p256q,bef=c3fk12_345678(i3)+,aft=
      c3fk12_345678(i3)=c3fk12_345678(i3)+(l1_3478(i3)%c(2)*p256
     & q*r2_56%b(1)+l1_3478(i3)%a(2)*r2_56%a(2))
      end do
      else
      do i3=1,2
* TLTR0_W -- aa=c3fk12_345678(i3),a1=l1_5678%a,c1=l1_5678%c,a2=r2_34(i3)
* %a,b2=r2_34(i3)%b,prq=p234q,bef=,aft=
      c3fk12_345678(i3)=(l1_5678%c(2)*p234q*r2_34(i3)%b(1)+l1_56
     & 78%a(2)*r2_34(i3)%a(2))
      end do
  
      do i3=1,2
* TLTR0_W -- aa=c3fk12_345678(i3),a1=l1_3456(i3)%a,c1=l1_3456(i3)%c,a2=r
* 2_78%a,b2=r2_78%b,prq=p278q,bef=c3fk12_345678(i3)+,aft=
      c3fk12_345678(i3)=c3fk12_345678(i3)+(l1_3456(i3)%c(2)*p278
     & q*r2_78%b(1)+l1_3456(i3)%a(2)*r2_78%a(2))
      end do
      endif
      if (iup(id3).eq.1) then
      do i3=1,2
* TLTR0_W -- aa=c3fk34_125678(i3),a1=l3_7856%a,c1=l3_7856%c,a2=r4_12(i3)
* %a,b2=r4_12(i3)%b,prq=p412q,bef=,aft=
      c3fk34_125678(i3)=(l3_7856%c(2)*p412q*r4_12(i3)%b(1)+l3_78
     & 56%a(2)*r4_12(i3)%a(2))
      end do
  
      do i3=1,2
* TLTR0_W -- aa=c3fk34_125678(i3),a1=l3_1278(i3)%a,c1=l3_1278(i3)%c,a2=r
* 4_56%a,b2=r4_56%b,prq=p456q,bef=c3fk34_125678(i3)+,aft=
      c3fk34_125678(i3)=c3fk34_125678(i3)+(l3_1278(i3)%c(2)*p456
     & q*r4_56%b(1)+l3_1278(i3)%a(2)*r4_56%a(2))
      end do
      else
      do i3=1,2
* TLTR0_W -- aa=c3fk34_125678(i3),a1=l3_5678%a,c1=l3_5678%c,a2=r4_12(i3)
* %a,b2=r4_12(i3)%b,prq=p412q,bef=,aft=
      c3fk34_125678(i3)=(l3_5678%c(2)*p412q*r4_12(i3)%b(1)+l3_56
     & 78%a(2)*r4_12(i3)%a(2))
      end do
  
      do i3=1,2
* TLTR0_W -- aa=c3fk34_125678(i3),a1=l3_1256(i3)%a,c1=l3_1256(i3)%c,a2=r
* 4_78%a,b2=r4_78%b,prq=p478q,bef=c3fk34_125678(i3)+,aft=
      c3fk34_125678(i3)=c3fk34_125678(i3)+(l3_1256(i3)%c(2)*p478
     & q*r4_78%b(1)+l3_1256(i3)%a(2)*r4_78%a(2))
      end do
      endif
  
* main line W                                                           
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=c3fk56_123478(i1,i3),a1=l5_1278(i1)%a,c1=l5_1278(i1)%c,a
* 2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=p634q,bef=,aft=
      c3fk56_123478(i1,i3)=(l5_1278(i1)%c(2)*p634q*r6_34(i3)%b(1
     & )+l5_1278(i1)%a(2)*r6_34(i3)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=c3fk56_123478(i1,i3),a1=l5_3478(i3)%a,c1=l5_3478(i3)%c,a
* 2=r6_12(i1)%a,b2=r6_12(i1)%b,prq=p612q,bef=c3fk56_123478(i1,i3)+,aft=
      c3fk56_123478(i1,i3)=c3fk56_123478(i1,i3)+(l5_3478(i3)%c(2
     & )*p612q*r6_12(i1)%b(1)+l5_3478(i3)%a(2)*r6_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=c3fk56_123478(i1,i3),a1=l5_1234(i1,i3)%a,c1=l5_1234(i1,i
* 3)%c,a2=r6_78%a,b2=r6_78%b,prq=p678q,bef=c3fk56_123478(i1,i3)+,aft=
      c3fk56_123478(i1,i3)=c3fk56_123478(i1,i3)+(l5_1234(i1,i3)%
     & c(2)*p678q*r6_78%b(1)+l5_1234(i1,i3)%a(2)*r6_78%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=c3fk78_123456(i1,i3),a1=l7_1256(i1)%a,c1=l7_1256(i1)%c,a
* 2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=p834q,bef=,aft=
      c3fk78_123456(i1,i3)=(l7_1256(i1)%c(2)*p834q*r8_34(i3)%b(1
     & )+l7_1256(i1)%a(2)*r8_34(i3)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=c3fk78_123456(i1,i3),a1=l7_3456(i3)%a,c1=l7_3456(i3)%c,a
* 2=r8_12(i1)%a,b2=r8_12(i1)%b,prq=p812q,bef=c3fk78_123456(i1,i3)+,aft=
      c3fk78_123456(i1,i3)=c3fk78_123456(i1,i3)+(l7_3456(i3)%c(2
     & )*p812q*r8_12(i1)%b(1)+l7_3456(i3)%a(2)*r8_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=c3fk78_123456(i1,i3),a1=l7_1234(i1,i3)%a,c1=l7_1234(i1,i
* 3)%c,a2=r8_56%a,b2=r8_56%b,prq=p856q,bef=c3fk78_123456(i1,i3)+,aft=
      c3fk78_123456(i1,i3)=c3fk78_123456(i1,i3)+(l7_1234(i1,i3)%
     & c(2)*p856q*r8_56%b(1)+l7_1234(i1,i3)%a(2)*r8_56%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
        c3fk_tot(i1,i3)= c3fk56_123478(i1,i3)+c3fk78_123456(i1,i3)
      enddo
      enddo
  
      do i3=1,2
       c3fk_tot(2,i3) = c3fk_tot(2,i3)+c3fk12_345678(i3)
      enddo
      do i1=1,2
       c3fk_tot(i1,2) = c3fk_tot(i1,2)+c3fk34_125678(i1)
      enddo
  
**QCD**                                                                 
* main line Z                                                           
      if (iup(id1).eq.1) then
      do i3=1,2
* TLTR0_W -- aa=cg3fk12_345678(i3),a1=l1_7856%a,c1=l1_7856%c,a2=rg2_34(i
* 3)%a,b2=rg2_34(i3)%b,prq=p234q,bef=,aft=
      cg3fk12_345678(i3)=(l1_7856%c(2)*p234q*rg2_34(i3)%b(1)+l1_
     & 7856%a(2)*rg2_34(i3)%a(2))
      end do
  
      do i3=1,2
* TLTR0_W -- aa=cg3fk12_345678(i3),a1=lg1_3478(i3)%a,c1=lg1_3478(i3)%c,a
* 2=r2_56%a,b2=r2_56%b,prq=p256q,bef=cg3fk12_345678(i3)+,aft=
      cg3fk12_345678(i3)=cg3fk12_345678(i3)+(lg1_3478(i3)%c(2)*p
     & 256q*r2_56%b(1)+lg1_3478(i3)%a(2)*r2_56%a(2))
      end do
      else
      do i3=1,2
* TLTR0_W -- aa=cg3fk12_345678(i3),a1=l1_5678%a,c1=l1_5678%c,a2=rg2_34(i
* 3)%a,b2=rg2_34(i3)%b,prq=p234q,bef=,aft=
      cg3fk12_345678(i3)=(l1_5678%c(2)*p234q*rg2_34(i3)%b(1)+l1_
     & 5678%a(2)*rg2_34(i3)%a(2))
      end do
  
      do i3=1,2
* TLTR0_W -- aa=cg3fk12_345678(i3),a1=lg1_3456(i3)%a,c1=lg1_3456(i3)%c,a
* 2=r2_78%a,b2=r2_78%b,prq=p278q,bef=cg3fk12_345678(i3)+,aft=
      cg3fk12_345678(i3)=cg3fk12_345678(i3)+(lg1_3456(i3)%c(2)*p
     & 278q*r2_78%b(1)+lg1_3456(i3)%a(2)*r2_78%a(2))
      end do
      endif
      if (iup(id3).eq.1) then
      do i3=1,2
* TLTR0_W -- aa=cg3fk34_125678(i3),a1=l3_7856%a,c1=l3_7856%c,a2=rg4_12(i
* 3)%a,b2=rg4_12(i3)%b,prq=p412q,bef=,aft=
      cg3fk34_125678(i3)=(l3_7856%c(2)*p412q*rg4_12(i3)%b(1)+l3_
     & 7856%a(2)*rg4_12(i3)%a(2))
      end do
  
      do i3=1,2
* TLTR0_W -- aa=cg3fk34_125678(i3),a1=lg3_1278(i3)%a,c1=lg3_1278(i3)%c,a
* 2=r4_56%a,b2=r4_56%b,prq=p456q,bef=cg3fk34_125678(i3)+,aft=
      cg3fk34_125678(i3)=cg3fk34_125678(i3)+(lg3_1278(i3)%c(2)*p
     & 456q*r4_56%b(1)+lg3_1278(i3)%a(2)*r4_56%a(2))
      end do
      else
      do i3=1,2
* TLTR0_W -- aa=cg3fk34_125678(i3),a1=l3_5678%a,c1=l3_5678%c,a2=rg4_12(i
* 3)%a,b2=rg4_12(i3)%b,prq=p412q,bef=,aft=
      cg3fk34_125678(i3)=(l3_5678%c(2)*p412q*rg4_12(i3)%b(1)+l3_
     & 5678%a(2)*rg4_12(i3)%a(2))
      end do
  
      do i3=1,2
* TLTR0_W -- aa=cg3fk34_125678(i3),a1=lg3_1256(i3)%a,c1=lg3_1256(i3)%c,a
* 2=r4_78%a,b2=r4_78%b,prq=p478q,bef=cg3fk34_125678(i3)+,aft=
      cg3fk34_125678(i3)=cg3fk34_125678(i3)+(lg3_1256(i3)%c(2)*p
     & 478q*r4_78%b(1)+lg3_1256(i3)%a(2)*r4_78%a(2))
      end do
      endif
  
        cg3fk1234(2,2) = cg3fk12_345678(2)+cg3fk34_125678(2)
        cg3fk1234(1,2) = cg3fk34_125678(1)
        cg3fk1234(2,1) = cg3fk12_345678(1)
        cg3fk1234(1,1) = czero
  
**QCD main line W                                                       
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk56_123478(i1,i3),a1=lg5_1278(i1)%a,c1=lg5_1278(i1)%
* c,a2=r6_34(i3)%a,b2=r6_34(i3)%b,prq=p634q,bef=,aft=
      cg3fk56_123478(i1,i3)=(lg5_1278(i1)%c(2)*p634q*r6_34(i3)%b
     & (1)+lg5_1278(i1)%a(2)*r6_34(i3)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk56_123478(i1,i3),a1=l5_3478(i3)%a,c1=l5_3478(i3)%c,
* a2=rg6_12(i1)%a,b2=rg6_12(i1)%b,prq=p612q,bef=cg3fk56_123478(i1,i3)+,a
* ft=
      cg3fk56_123478(i1,i3)=cg3fk56_123478(i1,i3)+(l5_3478(i3)%c
     & (2)*p612q*rg6_12(i1)%b(1)+l5_3478(i3)%a(2)*rg6_12(i1)%a(2
     & ))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk56_123478(i1,i3),a1=lg5_1234(i1,i3)%a,c1=lg5_1234(i
* 1,i3)%c,a2=r6_78%a,b2=r6_78%b,prq=p678q,bef=cg3fk56_123478(i1,i3)+,aft
* =
      cg3fk56_123478(i1,i3)=cg3fk56_123478(i1,i3)+(lg5_1234(i1,i
     & 3)%c(2)*p678q*r6_78%b(1)+lg5_1234(i1,i3)%a(2)*r6_78%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk56_341278(i1,i3),a1=lg5_3478(i3)%a,c1=lg5_3478(i3)%
* c,a2=r6_12(i1)%a,b2=r6_12(i1)%b,prq=p612q,bef=,aft=
      cg3fk56_341278(i1,i3)=(lg5_3478(i3)%c(2)*p612q*r6_12(i1)%b
     & (1)+lg5_3478(i3)%a(2)*r6_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk56_341278(i1,i3),a1=l5_1278(i1)%a,c1=l5_1278(i1)%c,
* a2=rg6_34(i3)%a,b2=rg6_34(i3)%b,prq=p634q,bef=cg3fk56_341278(i1,i3)+,a
* ft=
      cg3fk56_341278(i1,i3)=cg3fk56_341278(i1,i3)+(l5_1278(i1)%c
     & (2)*p634q*rg6_34(i3)%b(1)+l5_1278(i1)%a(2)*rg6_34(i3)%a(2
     & ))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk56_341278(i1,i3),a1=lg5_3412(i3,i1)%a,c1=lg5_3412(i
* 3,i1)%c,a2=r6_78%a,b2=r6_78%b,prq=p678q,bef=cg3fk56_341278(i1,i3)+,aft
* =
      cg3fk56_341278(i1,i3)=cg3fk56_341278(i1,i3)+(lg5_3412(i3,i
     & 1)%c(2)*p678q*r6_78%b(1)+lg5_3412(i3,i1)%a(2)*r6_78%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk78_123456(i1,i3),a1=lg7_1256(i1)%a,c1=lg7_1256(i1)%
* c,a2=r8_34(i3)%a,b2=r8_34(i3)%b,prq=p834q,bef=,aft=
      cg3fk78_123456(i1,i3)=(lg7_1256(i1)%c(2)*p834q*r8_34(i3)%b
     & (1)+lg7_1256(i1)%a(2)*r8_34(i3)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk78_123456(i1,i3),a1=l7_3456(i3)%a,c1=l7_3456(i3)%c,
* a2=rg8_12(i1)%a,b2=rg8_12(i1)%b,prq=p812q,bef=cg3fk78_123456(i1,i3)+,a
* ft=
      cg3fk78_123456(i1,i3)=cg3fk78_123456(i1,i3)+(l7_3456(i3)%c
     & (2)*p812q*rg8_12(i1)%b(1)+l7_3456(i3)%a(2)*rg8_12(i1)%a(2
     & ))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk78_123456(i1,i3),a1=lg7_1234(i1,i3)%a,c1=lg7_1234(i
* 1,i3)%c,a2=r8_56%a,b2=r8_56%b,prq=p856q,bef=cg3fk78_123456(i1,i3)+,aft
* =
      cg3fk78_123456(i1,i3)=cg3fk78_123456(i1,i3)+(lg7_1234(i1,i
     & 3)%c(2)*p856q*r8_56%b(1)+lg7_1234(i1,i3)%a(2)*r8_56%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk78_341256(i1,i3),a1=lg7_3456(i3)%a,c1=lg7_3456(i3)%
* c,a2=r8_12(i1)%a,b2=r8_12(i1)%b,prq=p812q,bef=,aft=
      cg3fk78_341256(i1,i3)=(lg7_3456(i3)%c(2)*p812q*r8_12(i1)%b
     & (1)+lg7_3456(i3)%a(2)*r8_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk78_341256(i1,i3),a1=l7_1256(i1)%a,c1=l7_1256(i1)%c,
* a2=rg8_34(i3)%a,b2=rg8_34(i3)%b,prq=p834q,bef=cg3fk78_341256(i1,i3)+,a
* ft=
      cg3fk78_341256(i1,i3)=cg3fk78_341256(i1,i3)+(l7_1256(i1)%c
     & (2)*p834q*rg8_34(i3)%b(1)+l7_1256(i1)%a(2)*rg8_34(i3)%a(2
     & ))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* TLTR0_W -- aa=cg3fk78_341256(i1,i3),a1=lg7_3412(i3,i1)%a,c1=lg7_3412(i
* 3,i1)%c,a2=r8_56%a,b2=r8_56%b,prq=p856q,bef=cg3fk78_341256(i1,i3)+,aft
* =
      cg3fk78_341256(i1,i3)=cg3fk78_341256(i1,i3)+(lg7_3412(i3,i
     & 1)%c(2)*p856q*r8_56%b(1)+lg7_3412(i3,i1)%a(2)*r8_56%a(2))
      end do
      end do
  
  
  
*                                                                       
*  BOSON CONNECTION                                                     
*  FOUR-> FOUR ---------------                                          
*                                                                       
  
      do m=0,3
       p1234(m)=-p12(m)-p34(m)
      enddo
* p.q -- p.q=p1234q,p=p1234,q=p1234,bef=,aft=
      p1234q=(p1234(0)*p1234(0)-p1234(1)*p1234(1)-p1234(2)*p1234
     & (2)-p1234(3)*p1234(3))
  
      do m=0,3
       p1256(m)=-p12(m)-p56(m)
      enddo
* p.q -- p.q=p1256q,p=p1256,q=p1256,bef=,aft=
      p1256q=(p1256(0)*p1256(0)-p1256(1)*p1256(1)-p1256(2)*p1256
     & (2)-p1256(3)*p1256(3))
  
      do m=0,3
       p1278(m)=-p12(m)-p78(m)
      enddo
* p.q -- p.q=p1278q,p=p1278,q=p1278,bef=,aft=
      p1278q=(p1278(0)*p1278(0)-p1278(1)*p1278(1)-p1278(2)*p1278
     & (2)-p1278(3)*p1278(3))
  
      do m=0,3
       p3456(m)=-p34(m)-p56(m)
      enddo
* p.q -- p.q=p3456q,p=p3456,q=p3456,bef=,aft=
      p3456q=(p3456(0)*p3456(0)-p3456(1)*p3456(1)-p3456(2)*p3456
     & (2)-p3456(3)*p3456(3))
  
      do m=0,3
       p3478(m)=-p34(m)-p78(m)
      enddo
* p.q -- p.q=p3478q,p=p3478,q=p3478,bef=,aft=
      p3478q=(p3478(0)*p3478(0)-p3478(1)*p3478(1)-p3478(2)*p3478
     & (2)-p3478(3)*p3478(3))
  
      do m=0,3
       p5678(m)=-p56(m)-p78(m)
      enddo
* p.q -- p.q=p5678q,p=p5678,q=p5678,bef=,aft=
      p5678q=(p5678(0)*p5678(0)-p5678(1)*p5678(1)-p5678(2)*p5678
     & (2)-p5678(3)*p5678(3))
  
  
*                                                                       
* LEFT VECTOR                                                           
*                                                                       
  
*     Z boson                                                           
* quqd -- p=p1,q=p234
      quqd=p1(0)*p234(0)-p1(1)*p234(1)-p1(2)*p234(2)-p1(3)*p234(
     & 3)
* TL0 -- qu=p1,qd=p234,v=0,a=lz1_34(0)%a,c=lz1_34(0)%c,cr=zcr(id1),cl=zc
* l(id1),nsum=0
      eps_0=-p1(2)*p234(3)+p234(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p234(0)+p234k0*p1(0)
      lz1_34(0)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_34(0)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_34(0)%c(1)=zcr(id1)*(p1(2)+ceps_1)
      lz1_34(0)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p234,v=1,a=lz1_34(1)%a,c=lz1_34(1)%c,cr=zcr(id1),cl=zc
* l(id1),nsum=0
      auxa=-quqd+p1k0*p234(1)+p234k0*p1(1)
      lz1_34(1)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_34(1)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_34(1)%c(1)=zcr(id1)*(p1(2)+ceps_1)
      lz1_34(1)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p234,v=2,a=lz1_34(2)%a,c=lz1_34(2)%c,cr=zcr(id1),cl=zc
* l(id1),nsum=0
      eps_0=-p1k0*p234(3)+p234k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p234(2)+p234k0*p1(2)
      lz1_34(2)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_34(2)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_34(2)%c(1)=zcr(id1)*p1k0
      lz1_34(2)%c(2)=-zcl(id1)*p1k0
* TL0 -- qu=p1,qd=p234,v=3,a=lz1_34(3)%a,c=lz1_34(3)%c,cr=zcr(id1),cl=zc
* l(id1),nsum=0
      eps_0=p1k0*p234(2)-p234k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p234(3)+p234k0*p1(3)
      lz1_34(3)%a(1)=zcr(id1)*(auxa+ceps_0)
      lz1_34(3)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_34(3)%c(1)=zcr(id1)*ceps_1
      lz1_34(3)%c(2)=zcl(id1)*ceps_1
  
      if (ineutri(id1).ne.1) then
* TL0 -- qu=p1,qd=p234,v=0,a=lf1_34(0)%a,c=lf1_34(0)%c,cr=fcr(id1),cl=fc
* l(id1),nsum=0
      eps_0=-p1(2)*p234(3)+p234(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p234(0)+p234k0*p1(0)
      lf1_34(0)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_34(0)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_34(0)%c(1)=fcr(id1)*(p1(2)+ceps_1)
      lf1_34(0)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p234,v=1,a=lf1_34(1)%a,c=lf1_34(1)%c,cr=fcr(id1),cl=fc
* l(id1),nsum=0
      auxa=-quqd+p1k0*p234(1)+p234k0*p1(1)
      lf1_34(1)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_34(1)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_34(1)%c(1)=fcr(id1)*(p1(2)+ceps_1)
      lf1_34(1)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TL0 -- qu=p1,qd=p234,v=2,a=lf1_34(2)%a,c=lf1_34(2)%c,cr=fcr(id1),cl=fc
* l(id1),nsum=0
      eps_0=-p1k0*p234(3)+p234k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p234(2)+p234k0*p1(2)
      lf1_34(2)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_34(2)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_34(2)%c(1)=fcr(id1)*p1k0
      lf1_34(2)%c(2)=-fcl(id1)*p1k0
* TL0 -- qu=p1,qd=p234,v=3,a=lf1_34(3)%a,c=lf1_34(3)%c,cr=fcr(id1),cl=fc
* l(id1),nsum=0
      eps_0=p1k0*p234(2)-p234k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p234(3)+p234k0*p1(3)
      lf1_34(3)%a(1)=fcr(id1)*(auxa+ceps_0)
      lf1_34(3)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_34(3)%c(1)=fcr(id1)*ceps_1
      lf1_34(3)%c(2)=fcl(id1)*ceps_1
      endif
* quqd -- p=p3,q=p412
      quqd=p3(0)*p412(0)-p3(1)*p412(1)-p3(2)*p412(2)-p3(3)*p412(
     & 3)
* TL0 -- qu=p3,qd=p412,v=0,a=lz3_12(0)%a,c=lz3_12(0)%c,cr=zcr(id3),cl=zc
* l(id3),nsum=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lz3_12(0)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_12(0)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_12(0)%c(1)=zcr(id3)*(p3(2)+ceps_1)
      lz3_12(0)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=1,a=lz3_12(1)%a,c=lz3_12(1)%c,cr=zcr(id3),cl=zc
* l(id3),nsum=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lz3_12(1)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_12(1)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_12(1)%c(1)=zcr(id3)*(p3(2)+ceps_1)
      lz3_12(1)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=2,a=lz3_12(2)%a,c=lz3_12(2)%c,cr=zcr(id3),cl=zc
* l(id3),nsum=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lz3_12(2)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_12(2)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_12(2)%c(1)=zcr(id3)*p3k0
      lz3_12(2)%c(2)=-zcl(id3)*p3k0
* TL0 -- qu=p3,qd=p412,v=3,a=lz3_12(3)%a,c=lz3_12(3)%c,cr=zcr(id3),cl=zc
* l(id3),nsum=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p412(3)+p412k0*p3(3)
      lz3_12(3)%a(1)=zcr(id3)*(auxa+ceps_0)
      lz3_12(3)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_12(3)%c(1)=zcr(id3)*ceps_1
      lz3_12(3)%c(2)=zcl(id3)*ceps_1
  
      if (ineutri(id3).ne.1) then
* TL0 -- qu=p3,qd=p412,v=0,a=lf3_12(0)%a,c=lf3_12(0)%c,cr=fcr(id3),cl=fc
* l(id3),nsum=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lf3_12(0)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_12(0)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_12(0)%c(1)=fcr(id3)*(p3(2)+ceps_1)
      lf3_12(0)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=1,a=lf3_12(1)%a,c=lf3_12(1)%c,cr=fcr(id3),cl=fc
* l(id3),nsum=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lf3_12(1)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_12(1)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_12(1)%c(1)=fcr(id3)*(p3(2)+ceps_1)
      lf3_12(1)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TL0 -- qu=p3,qd=p412,v=2,a=lf3_12(2)%a,c=lf3_12(2)%c,cr=fcr(id3),cl=fc
* l(id3),nsum=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lf3_12(2)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_12(2)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_12(2)%c(1)=fcr(id3)*p3k0
      lf3_12(2)%c(2)=-fcl(id3)*p3k0
* TL0 -- qu=p3,qd=p412,v=3,a=lf3_12(3)%a,c=lf3_12(3)%c,cr=fcr(id3),cl=fc
* l(id3),nsum=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p412(3)+p412k0*p3(3)
      lf3_12(3)%a(1)=fcr(id3)*(auxa+ceps_0)
      lf3_12(3)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_12(3)%c(1)=fcr(id3)*ceps_1
      lf3_12(3)%c(2)=fcl(id3)*ceps_1
      endif
  
* quqd -- p=p5,q=p678
      quqd=p5(0)*p678(0)-p5(1)*p678(1)-p5(2)*p678(2)-p5(3)*p678(
     & 3)
* TWL0 -- qu=p5,qd=p678,v=0,a=lz5_78(0)%a,c=lz5_78(0)%c,cl=zcl(id5),nsum
* =0
      eps_0=-p5(2)*p678(3)+p678(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p678(0)+p678k0*p5(0)
      lz5_78(0)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_78(0)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p678,v=1,a=lz5_78(1)%a,c=lz5_78(1)%c,cl=zcl(id5),nsum
* =0
      auxa=-quqd+p5k0*p678(1)+p678k0*p5(1)
      lz5_78(1)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_78(1)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p678,v=2,a=lz5_78(2)%a,c=lz5_78(2)%c,cl=zcl(id5),nsum
* =0
      eps_0=-p5k0*p678(3)+p678k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p678(2)+p678k0*p5(2)
      lz5_78(2)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_78(2)%c(2)=-zcl(id5)*p5k0
* TWL0 -- qu=p5,qd=p678,v=3,a=lz5_78(3)%a,c=lz5_78(3)%c,cl=zcl(id5),nsum
* =0
      eps_0=p5k0*p678(2)-p678k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p678(3)+p678k0*p5(3)
      lz5_78(3)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_78(3)%c(2)=zcl(id5)*ceps_1
  
      if (ineutri(id5).ne.1) then
* TWL0 -- qu=p5,qd=p678,v=0,a=lf5_78(0)%a,c=lf5_78(0)%c,cl=fcl(id5),nsum
* =0
      eps_0=-p5(2)*p678(3)+p678(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p678(0)+p678k0*p5(0)
      lf5_78(0)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_78(0)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p678,v=1,a=lf5_78(1)%a,c=lf5_78(1)%c,cl=fcl(id5),nsum
* =0
      auxa=-quqd+p5k0*p678(1)+p678k0*p5(1)
      lf5_78(1)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_78(1)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p678,v=2,a=lf5_78(2)%a,c=lf5_78(2)%c,cl=fcl(id5),nsum
* =0
      eps_0=-p5k0*p678(3)+p678k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p678(2)+p678k0*p5(2)
      lf5_78(2)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_78(2)%c(2)=-fcl(id5)*p5k0
* TWL0 -- qu=p5,qd=p678,v=3,a=lf5_78(3)%a,c=lf5_78(3)%c,cl=fcl(id5),nsum
* =0
      eps_0=p5k0*p678(2)-p678k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p678(3)+p678k0*p5(3)
      lf5_78(3)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_78(3)%c(2)=fcl(id5)*ceps_1
      endif
* quqd -- p=p7,q=p856
      quqd=p7(0)*p856(0)-p7(1)*p856(1)-p7(2)*p856(2)-p7(3)*p856(
     & 3)
* TWL0 -- qu=p7,qd=p856,v=0,a=lz7_56(0)%a,c=lz7_56(0)%c,cl=zcl(id7),nsum
* =0
      eps_0=-p7(2)*p856(3)+p856(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p856(0)+p856k0*p7(0)
      lz7_56(0)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_56(0)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p856,v=1,a=lz7_56(1)%a,c=lz7_56(1)%c,cl=zcl(id7),nsum
* =0
      auxa=-quqd+p7k0*p856(1)+p856k0*p7(1)
      lz7_56(1)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_56(1)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p856,v=2,a=lz7_56(2)%a,c=lz7_56(2)%c,cl=zcl(id7),nsum
* =0
      eps_0=-p7k0*p856(3)+p856k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p856(2)+p856k0*p7(2)
      lz7_56(2)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_56(2)%c(2)=-zcl(id7)*p7k0
* TWL0 -- qu=p7,qd=p856,v=3,a=lz7_56(3)%a,c=lz7_56(3)%c,cl=zcl(id7),nsum
* =0
      eps_0=p7k0*p856(2)-p856k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p856(3)+p856k0*p7(3)
      lz7_56(3)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_56(3)%c(2)=zcl(id7)*ceps_1
  
      if (ineutri(id7).ne.1) then
* TWL0 -- qu=p7,qd=p856,v=0,a=lf7_56(0)%a,c=lf7_56(0)%c,cl=fcl(id7),nsum
* =0
      eps_0=-p7(2)*p856(3)+p856(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p856(0)+p856k0*p7(0)
      lf7_56(0)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_56(0)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p856,v=1,a=lf7_56(1)%a,c=lf7_56(1)%c,cl=fcl(id7),nsum
* =0
      auxa=-quqd+p7k0*p856(1)+p856k0*p7(1)
      lf7_56(1)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_56(1)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p856,v=2,a=lf7_56(2)%a,c=lf7_56(2)%c,cl=fcl(id7),nsum
* =0
      eps_0=-p7k0*p856(3)+p856k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p856(2)+p856k0*p7(2)
      lf7_56(2)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_56(2)%c(2)=-fcl(id7)*p7k0
* TWL0 -- qu=p7,qd=p856,v=3,a=lf7_56(3)%a,c=lf7_56(3)%c,cl=fcl(id7),nsum
* =0
      eps_0=p7k0*p856(2)-p856k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p856(3)+p856k0*p7(3)
      lf7_56(3)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_56(3)%c(2)=fcl(id7)*ceps_1
      endif
  
*     W boson                                                           
* quqd -- p=p5,q=p612
      quqd=p5(0)*p612(0)-p5(1)*p612(1)-p5(2)*p612(2)-p5(3)*p612(
     & 3)
* TWL0 -- qu=p5,qd=p612,v=0,a=lw5_12(0)%a,c=lw5_12(0)%c,cl=wcl,nsum=0
      eps_0=-p5(2)*p612(3)+p612(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p612(0)+p612k0*p5(0)
      lw5_12(0)%a(2)=wcl*(auxa-ceps_0)
      lw5_12(0)%c(2)=wcl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p612,v=1,a=lw5_12(1)%a,c=lw5_12(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p5k0*p612(1)+p612k0*p5(1)
      lw5_12(1)%a(2)=wcl*(auxa-ceps_0)
      lw5_12(1)%c(2)=wcl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p612,v=2,a=lw5_12(2)%a,c=lw5_12(2)%c,cl=wcl,nsum=0
      eps_0=-p5k0*p612(3)+p612k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p612(2)+p612k0*p5(2)
      lw5_12(2)%a(2)=wcl*(auxa-ceps_0)
      lw5_12(2)%c(2)=-wcl*p5k0
* TWL0 -- qu=p5,qd=p612,v=3,a=lw5_12(3)%a,c=lw5_12(3)%c,cl=wcl,nsum=0
      eps_0=p5k0*p612(2)-p612k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p612(3)+p612k0*p5(3)
      lw5_12(3)%a(2)=wcl*(auxa-ceps_0)
      lw5_12(3)%c(2)=wcl*ceps_1
* quqd -- p=p5,q=p634
      quqd=p5(0)*p634(0)-p5(1)*p634(1)-p5(2)*p634(2)-p5(3)*p634(
     & 3)
* TWL0 -- qu=p5,qd=p634,v=0,a=lw5_34(0)%a,c=lw5_34(0)%c,cl=wcl,nsum=0
      eps_0=-p5(2)*p634(3)+p634(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p634(0)+p634k0*p5(0)
      lw5_34(0)%a(2)=wcl*(auxa-ceps_0)
      lw5_34(0)%c(2)=wcl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p634,v=1,a=lw5_34(1)%a,c=lw5_34(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p5k0*p634(1)+p634k0*p5(1)
      lw5_34(1)%a(2)=wcl*(auxa-ceps_0)
      lw5_34(1)%c(2)=wcl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p634,v=2,a=lw5_34(2)%a,c=lw5_34(2)%c,cl=wcl,nsum=0
      eps_0=-p5k0*p634(3)+p634k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p634(2)+p634k0*p5(2)
      lw5_34(2)%a(2)=wcl*(auxa-ceps_0)
      lw5_34(2)%c(2)=-wcl*p5k0
* TWL0 -- qu=p5,qd=p634,v=3,a=lw5_34(3)%a,c=lw5_34(3)%c,cl=wcl,nsum=0
      eps_0=p5k0*p634(2)-p634k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p634(3)+p634k0*p5(3)
      lw5_34(3)%a(2)=wcl*(auxa-ceps_0)
      lw5_34(3)%c(2)=wcl*ceps_1
* quqd -- p=p7,q=p812
      quqd=p7(0)*p812(0)-p7(1)*p812(1)-p7(2)*p812(2)-p7(3)*p812(
     & 3)
* TWL0 -- qu=p7,qd=p812,v=0,a=lw7_12(0)%a,c=lw7_12(0)%c,cl=wcl,nsum=0
      eps_0=-p7(2)*p812(3)+p812(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p812(0)+p812k0*p7(0)
      lw7_12(0)%a(2)=wcl*(auxa-ceps_0)
      lw7_12(0)%c(2)=wcl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=1,a=lw7_12(1)%a,c=lw7_12(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p7k0*p812(1)+p812k0*p7(1)
      lw7_12(1)%a(2)=wcl*(auxa-ceps_0)
      lw7_12(1)%c(2)=wcl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=2,a=lw7_12(2)%a,c=lw7_12(2)%c,cl=wcl,nsum=0
      eps_0=-p7k0*p812(3)+p812k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p812(2)+p812k0*p7(2)
      lw7_12(2)%a(2)=wcl*(auxa-ceps_0)
      lw7_12(2)%c(2)=-wcl*p7k0
* TWL0 -- qu=p7,qd=p812,v=3,a=lw7_12(3)%a,c=lw7_12(3)%c,cl=wcl,nsum=0
      eps_0=p7k0*p812(2)-p812k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p812(3)+p812k0*p7(3)
      lw7_12(3)%a(2)=wcl*(auxa-ceps_0)
      lw7_12(3)%c(2)=wcl*ceps_1
* quqd -- p=p7,q=p834
      quqd=p7(0)*p834(0)-p7(1)*p834(1)-p7(2)*p834(2)-p7(3)*p834(
     & 3)
* TWL0 -- qu=p7,qd=p834,v=0,a=lw7_34(0)%a,c=lw7_34(0)%c,cl=wcl,nsum=0
      eps_0=-p7(2)*p834(3)+p834(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p834(0)+p834k0*p7(0)
      lw7_34(0)%a(2)=wcl*(auxa-ceps_0)
      lw7_34(0)%c(2)=wcl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p834,v=1,a=lw7_34(1)%a,c=lw7_34(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p7k0*p834(1)+p834k0*p7(1)
      lw7_34(1)%a(2)=wcl*(auxa-ceps_0)
      lw7_34(1)%c(2)=wcl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p834,v=2,a=lw7_34(2)%a,c=lw7_34(2)%c,cl=wcl,nsum=0
      eps_0=-p7k0*p834(3)+p834k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p834(2)+p834k0*p7(2)
      lw7_34(2)%a(2)=wcl*(auxa-ceps_0)
      lw7_34(2)%c(2)=-wcl*p7k0
* TWL0 -- qu=p7,qd=p834,v=3,a=lw7_34(3)%a,c=lw7_34(3)%c,cl=wcl,nsum=0
      eps_0=p7k0*p834(2)-p834k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p834(3)+p834k0*p7(3)
      lw7_34(3)%a(2)=wcl*(auxa-ceps_0)
      lw7_34(3)%c(2)=wcl*ceps_1
  
      if (iup(id1).eq.1) then
* quqd -- p=p1,q=p256
      quqd=p1(0)*p256(0)-p1(1)*p256(1)-p1(2)*p256(2)-p1(3)*p256(
     & 3)
* TWL0 -- qu=p1,qd=p256,v=0,a=lw1_56(0)%a,c=lw1_56(0)%c,cl=wcl,nsum=0
      eps_0=-p1(2)*p256(3)+p256(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p256(0)+p256k0*p1(0)
      lw1_56(0)%a(2)=wcl*(auxa-ceps_0)
      lw1_56(0)%c(2)=wcl*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p256,v=1,a=lw1_56(1)%a,c=lw1_56(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p1k0*p256(1)+p256k0*p1(1)
      lw1_56(1)%a(2)=wcl*(auxa-ceps_0)
      lw1_56(1)%c(2)=wcl*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p256,v=2,a=lw1_56(2)%a,c=lw1_56(2)%c,cl=wcl,nsum=0
      eps_0=-p1k0*p256(3)+p256k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p256(2)+p256k0*p1(2)
      lw1_56(2)%a(2)=wcl*(auxa-ceps_0)
      lw1_56(2)%c(2)=-wcl*p1k0
* TWL0 -- qu=p1,qd=p256,v=3,a=lw1_56(3)%a,c=lw1_56(3)%c,cl=wcl,nsum=0
      eps_0=p1k0*p256(2)-p256k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p256(3)+p256k0*p1(3)
      lw1_56(3)%a(2)=wcl*(auxa-ceps_0)
      lw1_56(3)%c(2)=wcl*ceps_1
      else
* quqd -- p=p1,q=p278
      quqd=p1(0)*p278(0)-p1(1)*p278(1)-p1(2)*p278(2)-p1(3)*p278(
     & 3)
* TWL0 -- qu=p1,qd=p278,v=0,a=lw1_78(0)%a,c=lw1_78(0)%c,cl=wcl,nsum=0
      eps_0=-p1(2)*p278(3)+p278(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p278(0)+p278k0*p1(0)
      lw1_78(0)%a(2)=wcl*(auxa-ceps_0)
      lw1_78(0)%c(2)=wcl*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p278,v=1,a=lw1_78(1)%a,c=lw1_78(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p1k0*p278(1)+p278k0*p1(1)
      lw1_78(1)%a(2)=wcl*(auxa-ceps_0)
      lw1_78(1)%c(2)=wcl*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p278,v=2,a=lw1_78(2)%a,c=lw1_78(2)%c,cl=wcl,nsum=0
      eps_0=-p1k0*p278(3)+p278k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p278(2)+p278k0*p1(2)
      lw1_78(2)%a(2)=wcl*(auxa-ceps_0)
      lw1_78(2)%c(2)=-wcl*p1k0
* TWL0 -- qu=p1,qd=p278,v=3,a=lw1_78(3)%a,c=lw1_78(3)%c,cl=wcl,nsum=0
      eps_0=p1k0*p278(2)-p278k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p278(3)+p278k0*p1(3)
      lw1_78(3)%a(2)=wcl*(auxa-ceps_0)
      lw1_78(3)%c(2)=wcl*ceps_1
      endif
      if (iup(id3).eq.1) then
* quqd -- p=p3,q=p456
      quqd=p3(0)*p456(0)-p3(1)*p456(1)-p3(2)*p456(2)-p3(3)*p456(
     & 3)
* TWL0 -- qu=p3,qd=p456,v=0,a=lw3_56(0)%a,c=lw3_56(0)%c,cl=wcl,nsum=0
      eps_0=-p3(2)*p456(3)+p456(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p456(0)+p456k0*p3(0)
      lw3_56(0)%a(2)=wcl*(auxa-ceps_0)
      lw3_56(0)%c(2)=wcl*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p456,v=1,a=lw3_56(1)%a,c=lw3_56(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p3k0*p456(1)+p456k0*p3(1)
      lw3_56(1)%a(2)=wcl*(auxa-ceps_0)
      lw3_56(1)%c(2)=wcl*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p456,v=2,a=lw3_56(2)%a,c=lw3_56(2)%c,cl=wcl,nsum=0
      eps_0=-p3k0*p456(3)+p456k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p456(2)+p456k0*p3(2)
      lw3_56(2)%a(2)=wcl*(auxa-ceps_0)
      lw3_56(2)%c(2)=-wcl*p3k0
* TWL0 -- qu=p3,qd=p456,v=3,a=lw3_56(3)%a,c=lw3_56(3)%c,cl=wcl,nsum=0
      eps_0=p3k0*p456(2)-p456k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p456(3)+p456k0*p3(3)
      lw3_56(3)%a(2)=wcl*(auxa-ceps_0)
      lw3_56(3)%c(2)=wcl*ceps_1
      else
* quqd -- p=p3,q=p478
      quqd=p3(0)*p478(0)-p3(1)*p478(1)-p3(2)*p478(2)-p3(3)*p478(
     & 3)
* TWL0 -- qu=p3,qd=p478,v=0,a=lw3_78(0)%a,c=lw3_78(0)%c,cl=wcl,nsum=0
      eps_0=-p3(2)*p478(3)+p478(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p478(0)+p478k0*p3(0)
      lw3_78(0)%a(2)=wcl*(auxa-ceps_0)
      lw3_78(0)%c(2)=wcl*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p478,v=1,a=lw3_78(1)%a,c=lw3_78(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p3k0*p478(1)+p478k0*p3(1)
      lw3_78(1)%a(2)=wcl*(auxa-ceps_0)
      lw3_78(1)%c(2)=wcl*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p478,v=2,a=lw3_78(2)%a,c=lw3_78(2)%c,cl=wcl,nsum=0
      eps_0=-p3k0*p478(3)+p478k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p478(2)+p478k0*p3(2)
      lw3_78(2)%a(2)=wcl*(auxa-ceps_0)
      lw3_78(2)%c(2)=-wcl*p3k0
* TWL0 -- qu=p3,qd=p478,v=3,a=lw3_78(3)%a,c=lw3_78(3)%c,cl=wcl,nsum=0
      eps_0=p3k0*p478(2)-p478k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p478(3)+p478k0*p3(3)
      lw3_78(3)%a(2)=wcl*(auxa-ceps_0)
      lw3_78(3)%c(2)=wcl*ceps_1
      endif
  
*                                                                       
* RIGHT VECTOR                                                          
*                                                                       
  
*     Z boson                                                           
* quqd -- p=p134,q=p2
      quqd=p134(0)*p2(0)-p134(1)*p2(1)-p134(2)*p2(2)-p134(3)*p2(
     & 3)
* TR0 -- qu=p134,qd=p2,v=0,a=rz2_34(0)%a,b=rz2_34(0)%b,cr=zcr(id2),cl=zc
* l(id2),nsum=0
      eps_0=-p134(2)*p2(3)+p2(2)*p134(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p134k0*p2(0)+p2k0*p134(0)
      rz2_34(0)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_34(0)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_34(0)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_34(0)%b(2)=zcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p134,qd=p2,v=1,a=rz2_34(1)%a,b=rz2_34(1)%b,cr=zcr(id2),cl=zc
* l(id2),nsum=0
      auxa=-quqd+p134k0*p2(1)+p2k0*p134(1)
      rz2_34(1)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_34(1)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_34(1)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_34(1)%b(2)=zcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p134,qd=p2,v=2,a=rz2_34(2)%a,b=rz2_34(2)%b,cr=zcr(id2),cl=zc
* l(id2),nsum=0
      eps_0=-p134k0*p2(3)+p2k0*p134(3)
      ceps_0=eps_0*cim
      auxa=p134k0*p2(2)+p2k0*p134(2)
      rz2_34(2)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_34(2)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_34(2)%b(1)=-zcl(id2)*p2k0
      rz2_34(2)%b(2)=zcr(id2)*p2k0
* TR0 -- qu=p134,qd=p2,v=3,a=rz2_34(3)%a,b=rz2_34(3)%b,cr=zcr(id2),cl=zc
* l(id2),nsum=0
      eps_0=p134k0*p2(2)-p2k0*p134(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p134k0*p2(3)+p2k0*p134(3)
      rz2_34(3)%a(1)=zcr(id2)*(auxa+ceps_0)
      rz2_34(3)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_34(3)%b(1)=-zcl(id2)*ceps_2
      rz2_34(3)%b(2)=-zcr(id2)*ceps_2
  
      if (ineutri(id2).ne.1) then
* TR0 -- qu=p134,qd=p2,v=0,a=rf2_34(0)%a,b=rf2_34(0)%b,cr=fcr(id2),cl=fc
* l(id2),nsum=0
      eps_0=-p134(2)*p2(3)+p2(2)*p134(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p134k0*p2(0)+p2k0*p134(0)
      rf2_34(0)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_34(0)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_34(0)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_34(0)%b(2)=fcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p134,qd=p2,v=1,a=rf2_34(1)%a,b=rf2_34(1)%b,cr=fcr(id2),cl=fc
* l(id2),nsum=0
      auxa=-quqd+p134k0*p2(1)+p2k0*p134(1)
      rf2_34(1)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_34(1)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_34(1)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_34(1)%b(2)=fcr(id2)*(p2(2)-ceps_2)
* TR0 -- qu=p134,qd=p2,v=2,a=rf2_34(2)%a,b=rf2_34(2)%b,cr=fcr(id2),cl=fc
* l(id2),nsum=0
      eps_0=-p134k0*p2(3)+p2k0*p134(3)
      ceps_0=eps_0*cim
      auxa=p134k0*p2(2)+p2k0*p134(2)
      rf2_34(2)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_34(2)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_34(2)%b(1)=-fcl(id2)*p2k0
      rf2_34(2)%b(2)=fcr(id2)*p2k0
* TR0 -- qu=p134,qd=p2,v=3,a=rf2_34(3)%a,b=rf2_34(3)%b,cr=fcr(id2),cl=fc
* l(id2),nsum=0
      eps_0=p134k0*p2(2)-p2k0*p134(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p134k0*p2(3)+p2k0*p134(3)
      rf2_34(3)%a(1)=fcr(id2)*(auxa+ceps_0)
      rf2_34(3)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_34(3)%b(1)=-fcl(id2)*ceps_2
      rf2_34(3)%b(2)=-fcr(id2)*ceps_2
      endif
* quqd -- p=p312,q=p4
      quqd=p312(0)*p4(0)-p312(1)*p4(1)-p312(2)*p4(2)-p312(3)*p4(
     & 3)
* TR0 -- qu=p312,qd=p4,v=0,a=rz4_12(0)%a,b=rz4_12(0)%b,cr=zcr(id4),cl=zc
* l(id4),nsum=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rz4_12(0)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_12(0)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_12(0)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_12(0)%b(2)=zcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=1,a=rz4_12(1)%a,b=rz4_12(1)%b,cr=zcr(id4),cl=zc
* l(id4),nsum=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rz4_12(1)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_12(1)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_12(1)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_12(1)%b(2)=zcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=2,a=rz4_12(2)%a,b=rz4_12(2)%b,cr=zcr(id4),cl=zc
* l(id4),nsum=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rz4_12(2)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_12(2)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_12(2)%b(1)=-zcl(id4)*p4k0
      rz4_12(2)%b(2)=zcr(id4)*p4k0
* TR0 -- qu=p312,qd=p4,v=3,a=rz4_12(3)%a,b=rz4_12(3)%b,cr=zcr(id4),cl=zc
* l(id4),nsum=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p312k0*p4(3)+p4k0*p312(3)
      rz4_12(3)%a(1)=zcr(id4)*(auxa+ceps_0)
      rz4_12(3)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_12(3)%b(1)=-zcl(id4)*ceps_2
      rz4_12(3)%b(2)=-zcr(id4)*ceps_2
  
      if (ineutri(id4).ne.1) then
* TR0 -- qu=p312,qd=p4,v=0,a=rf4_12(0)%a,b=rf4_12(0)%b,cr=fcr(id4),cl=fc
* l(id4),nsum=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rf4_12(0)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_12(0)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_12(0)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_12(0)%b(2)=fcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=1,a=rf4_12(1)%a,b=rf4_12(1)%b,cr=fcr(id4),cl=fc
* l(id4),nsum=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rf4_12(1)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_12(1)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_12(1)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_12(1)%b(2)=fcr(id4)*(p4(2)-ceps_2)
* TR0 -- qu=p312,qd=p4,v=2,a=rf4_12(2)%a,b=rf4_12(2)%b,cr=fcr(id4),cl=fc
* l(id4),nsum=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rf4_12(2)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_12(2)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_12(2)%b(1)=-fcl(id4)*p4k0
      rf4_12(2)%b(2)=fcr(id4)*p4k0
* TR0 -- qu=p312,qd=p4,v=3,a=rf4_12(3)%a,b=rf4_12(3)%b,cr=fcr(id4),cl=fc
* l(id4),nsum=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p312k0*p4(3)+p4k0*p312(3)
      rf4_12(3)%a(1)=fcr(id4)*(auxa+ceps_0)
      rf4_12(3)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_12(3)%b(1)=-fcl(id4)*ceps_2
      rf4_12(3)%b(2)=-fcr(id4)*ceps_2
      endif
  
* quqd -- p=p578,q=p6
      quqd=p578(0)*p6(0)-p578(1)*p6(1)-p578(2)*p6(2)-p578(3)*p6(
     & 3)
* TWR0 -- qu=p578,qd=p6,v=0,a=rz6_78(0)%a,b=rz6_78(0)%b,cl=zcl(id6),nsum
* =0
      eps_0=-p578(2)*p6(3)+p6(2)*p578(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p578k0*p6(0)+p6k0*p578(0)
      rz6_78(0)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_78(0)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p578,qd=p6,v=1,a=rz6_78(1)%a,b=rz6_78(1)%b,cl=zcl(id6),nsum
* =0
      auxa=-quqd+p578k0*p6(1)+p6k0*p578(1)
      rz6_78(1)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_78(1)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p578,qd=p6,v=2,a=rz6_78(2)%a,b=rz6_78(2)%b,cl=zcl(id6),nsum
* =0
      eps_0=-p578k0*p6(3)+p6k0*p578(3)
      ceps_0=eps_0*cim
      auxa=p578k0*p6(2)+p6k0*p578(2)
      rz6_78(2)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_78(2)%b(1)=-zcl(id6)*p6k0
* TWR0 -- qu=p578,qd=p6,v=3,a=rz6_78(3)%a,b=rz6_78(3)%b,cl=zcl(id6),nsum
* =0
      eps_0=p578k0*p6(2)-p6k0*p578(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p578k0*p6(3)+p6k0*p578(3)
      rz6_78(3)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_78(3)%b(1)=-zcl(id6)*ceps_2
  
      if (ineutri(id6).ne.1) then
* TWR0 -- qu=p578,qd=p6,v=0,a=rf6_78(0)%a,b=rf6_78(0)%b,cl=fcl(id6),nsum
* =0
      eps_0=-p578(2)*p6(3)+p6(2)*p578(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p578k0*p6(0)+p6k0*p578(0)
      rf6_78(0)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_78(0)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p578,qd=p6,v=1,a=rf6_78(1)%a,b=rf6_78(1)%b,cl=fcl(id6),nsum
* =0
      auxa=-quqd+p578k0*p6(1)+p6k0*p578(1)
      rf6_78(1)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_78(1)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p578,qd=p6,v=2,a=rf6_78(2)%a,b=rf6_78(2)%b,cl=fcl(id6),nsum
* =0
      eps_0=-p578k0*p6(3)+p6k0*p578(3)
      ceps_0=eps_0*cim
      auxa=p578k0*p6(2)+p6k0*p578(2)
      rf6_78(2)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_78(2)%b(1)=-fcl(id6)*p6k0
* TWR0 -- qu=p578,qd=p6,v=3,a=rf6_78(3)%a,b=rf6_78(3)%b,cl=fcl(id6),nsum
* =0
      eps_0=p578k0*p6(2)-p6k0*p578(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p578k0*p6(3)+p6k0*p578(3)
      rf6_78(3)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_78(3)%b(1)=-fcl(id6)*ceps_2
      endif
* quqd -- p=p756,q=p8
      quqd=p756(0)*p8(0)-p756(1)*p8(1)-p756(2)*p8(2)-p756(3)*p8(
     & 3)
* TWR0 -- qu=p756,qd=p8,v=0,a=rz8_56(0)%a,b=rz8_56(0)%b,cl=zcl(id8),nsum
* =0
      eps_0=-p756(2)*p8(3)+p8(2)*p756(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p756k0*p8(0)+p8k0*p756(0)
      rz8_56(0)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_56(0)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p756,qd=p8,v=1,a=rz8_56(1)%a,b=rz8_56(1)%b,cl=zcl(id8),nsum
* =0
      auxa=-quqd+p756k0*p8(1)+p8k0*p756(1)
      rz8_56(1)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_56(1)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p756,qd=p8,v=2,a=rz8_56(2)%a,b=rz8_56(2)%b,cl=zcl(id8),nsum
* =0
      eps_0=-p756k0*p8(3)+p8k0*p756(3)
      ceps_0=eps_0*cim
      auxa=p756k0*p8(2)+p8k0*p756(2)
      rz8_56(2)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_56(2)%b(1)=-zcl(id8)*p8k0
* TWR0 -- qu=p756,qd=p8,v=3,a=rz8_56(3)%a,b=rz8_56(3)%b,cl=zcl(id8),nsum
* =0
      eps_0=p756k0*p8(2)-p8k0*p756(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p756k0*p8(3)+p8k0*p756(3)
      rz8_56(3)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_56(3)%b(1)=-zcl(id8)*ceps_2
  
      if (ineutri(id8).ne.1) then
* TWR0 -- qu=p756,qd=p8,v=0,a=rf8_56(0)%a,b=rf8_56(0)%b,cl=fcl(id8),nsum
* =0
      eps_0=-p756(2)*p8(3)+p8(2)*p756(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p756k0*p8(0)+p8k0*p756(0)
      rf8_56(0)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_56(0)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p756,qd=p8,v=1,a=rf8_56(1)%a,b=rf8_56(1)%b,cl=fcl(id8),nsum
* =0
      auxa=-quqd+p756k0*p8(1)+p8k0*p756(1)
      rf8_56(1)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_56(1)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p756,qd=p8,v=2,a=rf8_56(2)%a,b=rf8_56(2)%b,cl=fcl(id8),nsum
* =0
      eps_0=-p756k0*p8(3)+p8k0*p756(3)
      ceps_0=eps_0*cim
      auxa=p756k0*p8(2)+p8k0*p756(2)
      rf8_56(2)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_56(2)%b(1)=-fcl(id8)*p8k0
* TWR0 -- qu=p756,qd=p8,v=3,a=rf8_56(3)%a,b=rf8_56(3)%b,cl=fcl(id8),nsum
* =0
      eps_0=p756k0*p8(2)-p8k0*p756(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p756k0*p8(3)+p8k0*p756(3)
      rf8_56(3)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_56(3)%b(1)=-fcl(id8)*ceps_2
      endif
  
  
*     W boson                                                           
* quqd -- p=p512,q=p6
      quqd=p512(0)*p6(0)-p512(1)*p6(1)-p512(2)*p6(2)-p512(3)*p6(
     & 3)
* TWR0 -- qu=p512,qd=p6,v=0,a=rw6_12(0)%a,b=rw6_12(0)%b,cl=wcl,nsum=0
      eps_0=-p512(2)*p6(3)+p6(2)*p512(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p512k0*p6(0)+p6k0*p512(0)
      rw6_12(0)%a(2)=wcl*(auxa-ceps_0)
      rw6_12(0)%b(1)=-wcl*(p6(2)+ceps_2)
* TWR0 -- qu=p512,qd=p6,v=1,a=rw6_12(1)%a,b=rw6_12(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p512k0*p6(1)+p6k0*p512(1)
      rw6_12(1)%a(2)=wcl*(auxa-ceps_0)
      rw6_12(1)%b(1)=-wcl*(p6(2)+ceps_2)
* TWR0 -- qu=p512,qd=p6,v=2,a=rw6_12(2)%a,b=rw6_12(2)%b,cl=wcl,nsum=0
      eps_0=-p512k0*p6(3)+p6k0*p512(3)
      ceps_0=eps_0*cim
      auxa=p512k0*p6(2)+p6k0*p512(2)
      rw6_12(2)%a(2)=wcl*(auxa-ceps_0)
      rw6_12(2)%b(1)=-wcl*p6k0
* TWR0 -- qu=p512,qd=p6,v=3,a=rw6_12(3)%a,b=rw6_12(3)%b,cl=wcl,nsum=0
      eps_0=p512k0*p6(2)-p6k0*p512(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p512k0*p6(3)+p6k0*p512(3)
      rw6_12(3)%a(2)=wcl*(auxa-ceps_0)
      rw6_12(3)%b(1)=-wcl*ceps_2
* quqd -- p=p534,q=p6
      quqd=p534(0)*p6(0)-p534(1)*p6(1)-p534(2)*p6(2)-p534(3)*p6(
     & 3)
* TWR0 -- qu=p534,qd=p6,v=0,a=rw6_34(0)%a,b=rw6_34(0)%b,cl=wcl,nsum=0
      eps_0=-p534(2)*p6(3)+p6(2)*p534(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p534k0*p6(0)+p6k0*p534(0)
      rw6_34(0)%a(2)=wcl*(auxa-ceps_0)
      rw6_34(0)%b(1)=-wcl*(p6(2)+ceps_2)
* TWR0 -- qu=p534,qd=p6,v=1,a=rw6_34(1)%a,b=rw6_34(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p534k0*p6(1)+p6k0*p534(1)
      rw6_34(1)%a(2)=wcl*(auxa-ceps_0)
      rw6_34(1)%b(1)=-wcl*(p6(2)+ceps_2)
* TWR0 -- qu=p534,qd=p6,v=2,a=rw6_34(2)%a,b=rw6_34(2)%b,cl=wcl,nsum=0
      eps_0=-p534k0*p6(3)+p6k0*p534(3)
      ceps_0=eps_0*cim
      auxa=p534k0*p6(2)+p6k0*p534(2)
      rw6_34(2)%a(2)=wcl*(auxa-ceps_0)
      rw6_34(2)%b(1)=-wcl*p6k0
* TWR0 -- qu=p534,qd=p6,v=3,a=rw6_34(3)%a,b=rw6_34(3)%b,cl=wcl,nsum=0
      eps_0=p534k0*p6(2)-p6k0*p534(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p534k0*p6(3)+p6k0*p534(3)
      rw6_34(3)%a(2)=wcl*(auxa-ceps_0)
      rw6_34(3)%b(1)=-wcl*ceps_2
* quqd -- p=p712,q=p8
      quqd=p712(0)*p8(0)-p712(1)*p8(1)-p712(2)*p8(2)-p712(3)*p8(
     & 3)
* TWR0 -- qu=p712,qd=p8,v=0,a=rw8_12(0)%a,b=rw8_12(0)%b,cl=wcl,nsum=0
      eps_0=-p712(2)*p8(3)+p8(2)*p712(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p712k0*p8(0)+p8k0*p712(0)
      rw8_12(0)%a(2)=wcl*(auxa-ceps_0)
      rw8_12(0)%b(1)=-wcl*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=1,a=rw8_12(1)%a,b=rw8_12(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p712k0*p8(1)+p8k0*p712(1)
      rw8_12(1)%a(2)=wcl*(auxa-ceps_0)
      rw8_12(1)%b(1)=-wcl*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=2,a=rw8_12(2)%a,b=rw8_12(2)%b,cl=wcl,nsum=0
      eps_0=-p712k0*p8(3)+p8k0*p712(3)
      ceps_0=eps_0*cim
      auxa=p712k0*p8(2)+p8k0*p712(2)
      rw8_12(2)%a(2)=wcl*(auxa-ceps_0)
      rw8_12(2)%b(1)=-wcl*p8k0
* TWR0 -- qu=p712,qd=p8,v=3,a=rw8_12(3)%a,b=rw8_12(3)%b,cl=wcl,nsum=0
      eps_0=p712k0*p8(2)-p8k0*p712(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p712k0*p8(3)+p8k0*p712(3)
      rw8_12(3)%a(2)=wcl*(auxa-ceps_0)
      rw8_12(3)%b(1)=-wcl*ceps_2
* quqd -- p=p734,q=p8
      quqd=p734(0)*p8(0)-p734(1)*p8(1)-p734(2)*p8(2)-p734(3)*p8(
     & 3)
* TWR0 -- qu=p734,qd=p8,v=0,a=rw8_34(0)%a,b=rw8_34(0)%b,cl=wcl,nsum=0
      eps_0=-p734(2)*p8(3)+p8(2)*p734(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p734k0*p8(0)+p8k0*p734(0)
      rw8_34(0)%a(2)=wcl*(auxa-ceps_0)
      rw8_34(0)%b(1)=-wcl*(p8(2)+ceps_2)
* TWR0 -- qu=p734,qd=p8,v=1,a=rw8_34(1)%a,b=rw8_34(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p734k0*p8(1)+p8k0*p734(1)
      rw8_34(1)%a(2)=wcl*(auxa-ceps_0)
      rw8_34(1)%b(1)=-wcl*(p8(2)+ceps_2)
* TWR0 -- qu=p734,qd=p8,v=2,a=rw8_34(2)%a,b=rw8_34(2)%b,cl=wcl,nsum=0
      eps_0=-p734k0*p8(3)+p8k0*p734(3)
      ceps_0=eps_0*cim
      auxa=p734k0*p8(2)+p8k0*p734(2)
      rw8_34(2)%a(2)=wcl*(auxa-ceps_0)
      rw8_34(2)%b(1)=-wcl*p8k0
* TWR0 -- qu=p734,qd=p8,v=3,a=rw8_34(3)%a,b=rw8_34(3)%b,cl=wcl,nsum=0
      eps_0=p734k0*p8(2)-p8k0*p734(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p734k0*p8(3)+p8k0*p734(3)
      rw8_34(3)%a(2)=wcl*(auxa-ceps_0)
      rw8_34(3)%b(1)=-wcl*ceps_2
  
      if (iup(id1).eq.1) then
* quqd -- p=p178,q=p2
      quqd=p178(0)*p2(0)-p178(1)*p2(1)-p178(2)*p2(2)-p178(3)*p2(
     & 3)
* TWR0 -- qu=p178,qd=p2,v=0,a=rw2_78(0)%a,b=rw2_78(0)%b,cl=wcl,nsum=0
      eps_0=-p178(2)*p2(3)+p2(2)*p178(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p178k0*p2(0)+p2k0*p178(0)
      rw2_78(0)%a(2)=wcl*(auxa-ceps_0)
      rw2_78(0)%b(1)=-wcl*(p2(2)+ceps_2)
* TWR0 -- qu=p178,qd=p2,v=1,a=rw2_78(1)%a,b=rw2_78(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p178k0*p2(1)+p2k0*p178(1)
      rw2_78(1)%a(2)=wcl*(auxa-ceps_0)
      rw2_78(1)%b(1)=-wcl*(p2(2)+ceps_2)
* TWR0 -- qu=p178,qd=p2,v=2,a=rw2_78(2)%a,b=rw2_78(2)%b,cl=wcl,nsum=0
      eps_0=-p178k0*p2(3)+p2k0*p178(3)
      ceps_0=eps_0*cim
      auxa=p178k0*p2(2)+p2k0*p178(2)
      rw2_78(2)%a(2)=wcl*(auxa-ceps_0)
      rw2_78(2)%b(1)=-wcl*p2k0
* TWR0 -- qu=p178,qd=p2,v=3,a=rw2_78(3)%a,b=rw2_78(3)%b,cl=wcl,nsum=0
      eps_0=p178k0*p2(2)-p2k0*p178(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p178k0*p2(3)+p2k0*p178(3)
      rw2_78(3)%a(2)=wcl*(auxa-ceps_0)
      rw2_78(3)%b(1)=-wcl*ceps_2
      else
* quqd -- p=p156,q=p2
      quqd=p156(0)*p2(0)-p156(1)*p2(1)-p156(2)*p2(2)-p156(3)*p2(
     & 3)
* TWR0 -- qu=p156,qd=p2,v=0,a=rw2_56(0)%a,b=rw2_56(0)%b,cl=wcl,nsum=0
      eps_0=-p156(2)*p2(3)+p2(2)*p156(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p156k0*p2(0)+p2k0*p156(0)
      rw2_56(0)%a(2)=wcl*(auxa-ceps_0)
      rw2_56(0)%b(1)=-wcl*(p2(2)+ceps_2)
* TWR0 -- qu=p156,qd=p2,v=1,a=rw2_56(1)%a,b=rw2_56(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p156k0*p2(1)+p2k0*p156(1)
      rw2_56(1)%a(2)=wcl*(auxa-ceps_0)
      rw2_56(1)%b(1)=-wcl*(p2(2)+ceps_2)
* TWR0 -- qu=p156,qd=p2,v=2,a=rw2_56(2)%a,b=rw2_56(2)%b,cl=wcl,nsum=0
      eps_0=-p156k0*p2(3)+p2k0*p156(3)
      ceps_0=eps_0*cim
      auxa=p156k0*p2(2)+p2k0*p156(2)
      rw2_56(2)%a(2)=wcl*(auxa-ceps_0)
      rw2_56(2)%b(1)=-wcl*p2k0
* TWR0 -- qu=p156,qd=p2,v=3,a=rw2_56(3)%a,b=rw2_56(3)%b,cl=wcl,nsum=0
      eps_0=p156k0*p2(2)-p2k0*p156(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p156k0*p2(3)+p2k0*p156(3)
      rw2_56(3)%a(2)=wcl*(auxa-ceps_0)
      rw2_56(3)%b(1)=-wcl*ceps_2
      endif
      if (iup(id3).eq.1) then
* quqd -- p=p378,q=p4
      quqd=p378(0)*p4(0)-p378(1)*p4(1)-p378(2)*p4(2)-p378(3)*p4(
     & 3)
* TWR0 -- qu=p378,qd=p4,v=0,a=rw4_78(0)%a,b=rw4_78(0)%b,cl=wcl,nsum=0
      eps_0=-p378(2)*p4(3)+p4(2)*p378(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p378k0*p4(0)+p4k0*p378(0)
      rw4_78(0)%a(2)=wcl*(auxa-ceps_0)
      rw4_78(0)%b(1)=-wcl*(p4(2)+ceps_2)
* TWR0 -- qu=p378,qd=p4,v=1,a=rw4_78(1)%a,b=rw4_78(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p378k0*p4(1)+p4k0*p378(1)
      rw4_78(1)%a(2)=wcl*(auxa-ceps_0)
      rw4_78(1)%b(1)=-wcl*(p4(2)+ceps_2)
* TWR0 -- qu=p378,qd=p4,v=2,a=rw4_78(2)%a,b=rw4_78(2)%b,cl=wcl,nsum=0
      eps_0=-p378k0*p4(3)+p4k0*p378(3)
      ceps_0=eps_0*cim
      auxa=p378k0*p4(2)+p4k0*p378(2)
      rw4_78(2)%a(2)=wcl*(auxa-ceps_0)
      rw4_78(2)%b(1)=-wcl*p4k0
* TWR0 -- qu=p378,qd=p4,v=3,a=rw4_78(3)%a,b=rw4_78(3)%b,cl=wcl,nsum=0
      eps_0=p378k0*p4(2)-p4k0*p378(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p378k0*p4(3)+p4k0*p378(3)
      rw4_78(3)%a(2)=wcl*(auxa-ceps_0)
      rw4_78(3)%b(1)=-wcl*ceps_2
      else
* quqd -- p=p356,q=p4
      quqd=p356(0)*p4(0)-p356(1)*p4(1)-p356(2)*p4(2)-p356(3)*p4(
     & 3)
* TWR0 -- qu=p356,qd=p4,v=0,a=rw4_56(0)%a,b=rw4_56(0)%b,cl=wcl,nsum=0
      eps_0=-p356(2)*p4(3)+p4(2)*p356(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p356k0*p4(0)+p4k0*p356(0)
      rw4_56(0)%a(2)=wcl*(auxa-ceps_0)
      rw4_56(0)%b(1)=-wcl*(p4(2)+ceps_2)
* TWR0 -- qu=p356,qd=p4,v=1,a=rw4_56(1)%a,b=rw4_56(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p356k0*p4(1)+p4k0*p356(1)
      rw4_56(1)%a(2)=wcl*(auxa-ceps_0)
      rw4_56(1)%b(1)=-wcl*(p4(2)+ceps_2)
* TWR0 -- qu=p356,qd=p4,v=2,a=rw4_56(2)%a,b=rw4_56(2)%b,cl=wcl,nsum=0
      eps_0=-p356k0*p4(3)+p4k0*p356(3)
      ceps_0=eps_0*cim
      auxa=p356k0*p4(2)+p4k0*p356(2)
      rw4_56(2)%a(2)=wcl*(auxa-ceps_0)
      rw4_56(2)%b(1)=-wcl*p4k0
* TWR0 -- qu=p356,qd=p4,v=3,a=rw4_56(3)%a,b=rw4_56(3)%b,cl=wcl,nsum=0
      eps_0=p356k0*p4(2)-p4k0*p356(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p356k0*p4(3)+p4k0*p356(3)
      rw4_56(3)%a(2)=wcl*(auxa-ceps_0)
      rw4_56(3)%b(1)=-wcl*ceps_2
      endif
  
*                                                                       
* LEFT*RIGHT VECTOR                                                     
*                                                                       
  
**EW**                                                                  
  
* 1234 Z,A,g                                                            
  
      do i3=1,2
      do m=0,3
* TLTR0 -- aa=cz1234(&,i3)%e(m),a1=lz1_34(m)%a,c1=lz1_34(m)%c,a2=r2_34(i
* 3)%a,b2=r2_34(i3)%b,prq=p234q,bef=,aft=
      cz1234(1,i3)%e(m)=(lz1_34(m)%a(1)*r2_34(i3)%a(1)+lz1_34(m)
     & %c(1)*p234q*r2_34(i3)%b(2))
      cz1234(2,i3)%e(m)=(lz1_34(m)%c(2)*p234q*r2_34(i3)%b(1)+lz1
     & _34(m)%a(2)*r2_34(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do m=0,3
* TLTR0 -- aa=cz1234(&,i3)%e(m),a1=l1_34(i3)%a,c1=l1_34(i3)%c,a2=rz2_34(
* m)%a,b2=rz2_34(m)%b,prq=p134q,bef=cz1234(&,i3)%e(m)+,aft=
      cz1234(1,i3)%e(m)=cz1234(1,i3)%e(m)+(l1_34(i3)%a(1)*rz2_34
     & (m)%a(1)+l1_34(i3)%c(1)*p134q*rz2_34(m)%b(2))
      cz1234(2,i3)%e(m)=cz1234(2,i3)%e(m)+(l1_34(i3)%c(2)*p134q*
     & rz2_34(m)%b(1)+l1_34(i3)%a(2)*rz2_34(m)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz1234(i1,&)%e(m),a1=lz3_12(m)%a,c1=lz3_12(m)%c,a2=r4_12(i
* 1)%a,b2=r4_12(i1)%b,prq=p412q,bef=cz1234(i1,&)%e(m)+,aft=
      cz1234(i1,1)%e(m)=cz1234(i1,1)%e(m)+(lz3_12(m)%a(1)*r4_12(
     & i1)%a(1)+lz3_12(m)%c(1)*p412q*r4_12(i1)%b(2))
      cz1234(i1,2)%e(m)=cz1234(i1,2)%e(m)+(lz3_12(m)%c(2)*p412q*
     & r4_12(i1)%b(1)+lz3_12(m)%a(2)*r4_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cz1234(i1,&)%e(m),a1=l3_12(i1)%a,c1=l3_12(i1)%c,a2=rz4_12(
* m)%a,b2=rz4_12(m)%b,prq=p312q,bef=cz1234(i1,&)%e(m)+,aft=
      cz1234(i1,1)%e(m)=cz1234(i1,1)%e(m)+(l3_12(i1)%a(1)*rz4_12
     & (m)%a(1)+l3_12(i1)%c(1)*p312q*rz4_12(m)%b(2))
      cz1234(i1,2)%e(m)=cz1234(i1,2)%e(m)+(l3_12(i1)%c(2)*p312q*
     & rz4_12(m)%b(1)+l3_12(i1)%a(2)*rz4_12(m)%a(2))
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
* TLTR0 -- aa=cf1234(&,i3)%e(mu),a1=lf1_34(mu)%a,c1=lf1_34(mu)%c,a2=r2_3
* 4(i3)%a,b2=r2_34(i3)%b,prq=p234q,bef=,aft=
      cf1234(1,i3)%e(mu)=(lf1_34(mu)%a(1)*r2_34(i3)%a(1)+lf1_34(
     & mu)%c(1)*p234q*r2_34(i3)%b(2))
      cf1234(2,i3)%e(mu)=(lf1_34(mu)%c(2)*p234q*r2_34(i3)%b(1)+l
     & f1_34(mu)%a(2)*r2_34(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do mu=0,3
* TLTR0 -- aa=cf1234(&,i3)%e(mu),a1=l1_34(i3)%a,c1=l1_34(i3)%c,a2=rf2_34
* (mu)%a,b2=rf2_34(mu)%b,prq=p134q,bef=cf1234(&,i3)%e(mu)+,aft=
      cf1234(1,i3)%e(mu)=cf1234(1,i3)%e(mu)+(l1_34(i3)%a(1)*rf2_
     & 34(mu)%a(1)+l1_34(i3)%c(1)*p134q*rf2_34(mu)%b(2))
      cf1234(2,i3)%e(mu)=cf1234(2,i3)%e(mu)+(l1_34(i3)%c(2)*p134
     & q*rf2_34(mu)%b(1)+l1_34(i3)%a(2)*rf2_34(mu)%a(2))
      end do
      end do
  
**QCD gluon connection                                                  
       if (ilept(id1).ne.1) then
        do i1=1,2
        do i3=1,2
        do mu=0,3
         cg1234(i1,i3)%e(mu)=cf1234(i1,i3)%e(mu)/fcl(id1)
        enddo
        enddo
        enddo
       endif
  
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
* TLTR0 -- aa=cg3412(&,i1)%e(mu),a1=lf3_12(mu)%a,c1=lf3_12(mu)%c,a2=r4_1
* 2(i1)%a,b2=r4_12(i1)%b,prq=p412q,bef=,aft=/fcl(id3)
      cg3412(1,i1)%e(mu)=(lf3_12(mu)%a(1)*r4_12(i1)%a(1)+lf3_12(
     & mu)%c(1)*p412q*r4_12(i1)%b(2))/fcl(id3)
      cg3412(2,i1)%e(mu)=(lf3_12(mu)%c(2)*p412q*r4_12(i1)%b(1)+l
     & f3_12(mu)%a(2)*r4_12(i1)%a(2))/fcl(id3)
      end do
      end do
  
      do i1=1,2
      do mu=0,3
* TLTR0 -- aa=cg3412(&,i1)%e(mu),a1=l3_12(i1)%a,c1=l3_12(i1)%c,a2=rf4_12
* (mu)%a,b2=rf4_12(mu)%b,prq=p312q,bef=cg3412(&,i1)%e(mu)+,aft=/fcl(id3)
      cg3412(1,i1)%e(mu)=cg3412(1,i1)%e(mu)+(l3_12(i1)%a(1)*rf4_
     & 12(mu)%a(1)+l3_12(i1)%c(1)*p312q*rf4_12(mu)%b(2))/fcl(id3
     & )
      cg3412(2,i1)%e(mu)=cg3412(2,i1)%e(mu)+(l3_12(i1)%c(2)*p312
     & q*rf4_12(mu)%b(1)+l3_12(i1)%a(2)*rf4_12(mu)%a(2))/fcl(id3
     & )
      end do
      end do
  
        do i1=1,2
        do i3=1,2
        do mu=0,3
         cf1234(i1,i3)%e(mu)=cf1234(i1,i3)%e(mu)+
     &      cg3412(i3,i1)%e(mu)*fcl(id3)
        enddo
        enddo
        enddo
  
      endif
  
      do i1=1,2
      do i3=1,2
* pk0 -- p=cf1234(i1,i3)%e
      cf1234(i1,i3)%ek0=cf1234(i1,i3)%e(0)-cf1234(i1,i3)%e(1)
      end do
      end do
  
*   higgs                                                               
      if (rmh.ge.0.d0) then
  
       do i1=1,2
       do i3=1,2
  
* p.q -- p.q=caux,p=cz12(i1)%e,q=cz34(i3)%e,bef=,aft=
      caux=(cz12(i1)%e(0)*cz34(i3)%e(0)-cz12(i1)%e(1)*cz34(i3)%e
     & (1)-cz12(i1)%e(2)*cz34(i3)%e(2)-cz12(i1)%e(3)*cz34(i3)%e(
     & 3))
  
        ch1234(i1,i3)=rmz/sw/rcw*caux
  
       enddo
       enddo
      endif
  
* 5678 Z                                                                
  
* triple vertex                                                         
*** triple vertex -- pfz(mu)=p5678(mu),pwm(mu)=p78(mu),pwp(mu)=p56(mu),e
* fz=#,ewm=cw78,ewp=cw56,res=cvaux(mu#)
      do mu=0,3
      vfz(mu)=p78(mu)-p56(mu)
      vwm(mu)=p56(mu)-p5678(mu)
      vwp(mu)=p5678(mu)-p78(mu)
      end do !mu
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* ewm%ewp
* p.q -- p.q=caux,p=cw78%e,q=cw56%e,bef=,aft=
      caux=(cw78%e(0)*cw56%e(0)-cw78%e(1)*cw56%e(1)-cw78%e(2)*cw
     & 56%e(2)-cw78%e(3)*cw56%e(3))
      do mu=0,3
      cvaux(mu)=vfz(mu)*caux+cw78%v*cw56%e(mu)+cw56%v*cw78%e(mu)
      end do
  
      do m=0,3
       cz5678%e(m)= rcotw*cvaux(m)
      enddo
  
      do m=0,3
* TLTR0_W -- aa=cz5678%e(m),a1=lz5_78(m)%a,c1=lz5_78(m)%c,a2=r6_78%a,b2=
* r6_78%b,prq=p678q,bef=cz5678%e(m)+,aft=
      cz5678%e(m)=cz5678%e(m)+(lz5_78(m)%c(2)*p678q*r6_78%b(1)+l
     & z5_78(m)%a(2)*r6_78%a(2))
      end do
  
      do m=0,3
* TLTR0_W -- aa=cz5678%e(m),a1=lz7_56(m)%a,c1=lz7_56(m)%c,a2=r8_56%a,b2=
* r8_56%b,prq=p856q,bef=cz5678%e(m)+,aft=
      cz5678%e(m)=cz5678%e(m)+(lz7_56(m)%c(2)*p856q*r8_56%b(1)+l
     & z7_56(m)%a(2)*r8_56%a(2))
      end do
  
      do m=0,3
* TLTR0_W -- aa=cz5678%e(m),a1=l5_78%a,c1=l5_78%c,a2=rz6_78(m)%a,b2=rz6_
* 78(m)%b,prq=p578q,bef=cz5678%e(m)+,aft=
      cz5678%e(m)=cz5678%e(m)+(l5_78%c(2)*p578q*rz6_78(m)%b(1)+l
     & 5_78%a(2)*rz6_78(m)%a(2))
      end do
  
      do m=0,3
* TLTR0_W -- aa=cz5678%e(m),a1=l7_56%a,c1=l7_56%c,a2=rz8_56(m)%a,b2=rz8_
* 56(m)%b,prq=p756q,bef=cz5678%e(m)+,aft=
      cz5678%e(m)=cz5678%e(m)+(l7_56%c(2)*p756q*rz8_56(m)%b(1)+l
     & 7_56%a(2)*rz8_56(m)%a(2))
      end do
  
* foton                                                                 
      do m=0,3
       cf5678%e(m)= cvaux(m)
      enddo
  
**QCD gluon connection                                                  
  
      if (ilept(id5).ne.1) then
      do m=0,3
* TLTR0_W -- aa=cg5678%e(m),a1=lf5_78(m)%a,c1=lf5_78(m)%c,a2=r6_78%a,b2=
* r6_78%b,prq=p678q,bef=,aft=/fcl(id5)
      cg5678%e(m)=(lf5_78(m)%c(2)*p678q*r6_78%b(1)+lf5_78(m)%a(2
     & )*r6_78%a(2))/fcl(id5)
      end do
  
        do m=0,3
         cf5678%e(m)= cf5678%e(m)+cg5678%e(m)*fcl(id5)
        enddo
  
      do m=0,3
* TLTR0_W -- aa=cvaux(m),a1=l5_78%a,c1=l5_78%c,a2=rf6_78(m)%a,b2=rf6_78(
* m)%b,prq=p578q,bef=,aft=
      cvaux(m)=(l5_78%c(2)*p578q*rf6_78(m)%b(1)+l5_78%a(2)*rf6_7
     & 8(m)%a(2))
      end do
  
        do m=0,3
         cf5678%e(m)= cf5678%e(m)+cvaux(m)
         cg5678%e(m)= cg5678%e(m)+cvaux(m)/fcl(id6)
        enddo
  
      else
        if (ineutri(id5).ne.1) then
      do m=0,3
* TLTR0_W -- aa=cf5678%e(m),a1=lf5_78(m)%a,c1=lf5_78(m)%c,a2=r6_78%a,b2=
* r6_78%b,prq=p678q,bef=cf5678%e(m)+,aft=
      cf5678%e(m)=cf5678%e(m)+(lf5_78(m)%c(2)*p678q*r6_78%b(1)+l
     & f5_78(m)%a(2)*r6_78%a(2))
      end do
        endif
  
        if (ineutri(id6).ne.1) then
      do m=0,3
* TLTR0_W -- aa=cf5678%e(m),a1=l5_78%a,c1=l5_78%c,a2=rf6_78(m)%a,b2=rf6_
* 78(m)%b,prq=p578q,bef=cf5678%e(m)+,aft=
      cf5678%e(m)=cf5678%e(m)+(l5_78%c(2)*p578q*rf6_78(m)%b(1)+l
     & 5_78%a(2)*rf6_78(m)%a(2))
      end do
        endif
  
      endif
  
      if (ilept(id7).ne.1) then
  
      do m=0,3
* TLTR0_W -- aa=cg7856%e(m),a1=lf7_56(m)%a,c1=lf7_56(m)%c,a2=r8_56%a,b2=
* r8_56%b,prq=p856q,bef=,aft=/fcl(id7)
      cg7856%e(m)=(lf7_56(m)%c(2)*p856q*r8_56%b(1)+lf7_56(m)%a(2
     & )*r8_56%a(2))/fcl(id7)
      end do
  
        do m=0,3
         cf5678%e(m)= cf5678%e(m)+cg7856%e(m)*fcl(id7)
        enddo
  
      do m=0,3
* TLTR0_W -- aa=cvaux(m),a1=l7_56%a,c1=l7_56%c,a2=rf8_56(m)%a,b2=rf8_56(
* m)%b,prq=p756q,bef=,aft=
      cvaux(m)=(l7_56%c(2)*p756q*rf8_56(m)%b(1)+l7_56%a(2)*rf8_5
     & 6(m)%a(2))
      end do
  
        do m=0,3
         cf5678%e(m)= cf5678%e(m)+cvaux(m)
         cg7856%e(m)= cg7856%e(m)+cvaux(m)/fcl(id8)
        enddo
  
      else
        if (ineutri(id7).ne.1) then
      do m=0,3
* TLTR0_W -- aa=cf5678%e(m),a1=lf7_56(m)%a,c1=lf7_56(m)%c,a2=r8_56%a,b2=
* r8_56%b,prq=p856q,bef=cf5678%e(m)+,aft=
      cf5678%e(m)=cf5678%e(m)+(lf7_56(m)%c(2)*p856q*r8_56%b(1)+l
     & f7_56(m)%a(2)*r8_56%a(2))
      end do
        endif
  
        if (ineutri(id8).ne.1) then
      do m=0,3
* TLTR0_W -- aa=cf5678%e(m),a1=l7_56%a,c1=l7_56%c,a2=rf8_56(m)%a,b2=rf8_
* 56(m)%b,prq=p756q,bef=cf5678%e(m)+,aft=
      cf5678%e(m)=cf5678%e(m)+(l7_56%c(2)*p756q*rf8_56(m)%b(1)+l
     & 7_56%a(2)*rf8_56(m)%a(2))
      end do
        endif
  
      endif
  
*   higgs                                                               
      if (rmh.ge.0.d0) then
  
* p.q -- p.q=ch5678,p=cw56%e,q=cw78%e,bef=rmw/sw*,aft=
      ch5678=rmw/sw*(cw56%e(0)*cw78%e(0)-cw56%e(1)*cw78%e(1)-cw5
     & 6%e(2)*cw78%e(2)-cw56%e(3)*cw78%e(3))
  
      endif
  
* (1256,3456) W+                                                        
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw1256(i1)%e(m),a1=lw5_12(m)%a,c1=lw5_12(m)%c,a2=r6_12(i
* 1)%a,b2=r6_12(i1)%b,prq=p612q,bef=,aft=
      cw1256(i1)%e(m)=(lw5_12(m)%c(2)*p612q*r6_12(i1)%b(1)+lw5_1
     & 2(m)%a(2)*r6_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw1256(i1)%e(m),a1=l5_12(i1)%a,c1=l5_12(i1)%c,a2=rw6_12(
* m)%a,b2=rw6_12(m)%b,prq=p512q,bef=cw1256(i1)%e(m)+,aft=
      cw1256(i1)%e(m)=cw1256(i1)%e(m)+(l5_12(i1)%c(2)*p512q*rw6_
     & 12(m)%b(1)+l5_12(i1)%a(2)*rw6_12(m)%a(2))
      end do
      end do
  
      if (iup(id1).eq.1) then
      do m=0,3
* TLTR0_W -- aa=cw1256(2)%e(m),a1=lw1_56(m)%a,c1=lw1_56(m)%c,a2=r2_56%a,
* b2=r2_56%b,prq=p256q,bef=cw1256(2)%e(m)+,aft=
      cw1256(2)%e(m)=cw1256(2)%e(m)+(lw1_56(m)%c(2)*p256q*r2_56%
     & b(1)+lw1_56(m)%a(2)*r2_56%a(2))
      end do
      else
      do m=0,3
* TLTR0_W -- aa=cw1256(2)%e(m),a1=l1_56%a,c1=l1_56%c,a2=rw2_56(m)%a,b2=r
* w2_56(m)%b,prq=p156q,bef=cw1256(2)%e(m)+,aft=
      cw1256(2)%e(m)=cw1256(2)%e(m)+(l1_56%c(2)*p156q*rw2_56(m)%
     & b(1)+l1_56%a(2)*rw2_56(m)%a(2))
      end do
      endif
  
*** triple vertex -- pfz(mu)=p12(mu),pwm(mu)=p1256(mu),pwp(mu)=p56(mu),e
* fz=cz12(i1),ewm=#,ewp=cw56,res=cwaux(i1?,mu#)
      do mu=0,3
      vfz(mu)=p1256(mu)-p56(mu)
      vwm(mu)=p56(mu)-p12(mu)
      vwp(mu)=p12(mu)-p1256(mu)
      end do !mu
* vfz%efz
      do i1=1,2
* p.q -- p.q=cz12(i1)%v,p=cz12(i1)%e,q=vfz,bef=,aft=
      cz12(i1)%v=(cz12(i1)%e(0)*vfz(0)-cz12(i1)%e(1)*vfz(1)-cz12
     & (i1)%e(2)*vfz(2)-cz12(i1)%e(3)*vfz(3))
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewp
      do i1=1,2
* p.q -- p.q=caux,p=cz12(i1)%e,q=cw56%e,bef=,aft=
      caux=(cz12(i1)%e(0)*cw56%e(0)-cz12(i1)%e(1)*cw56%e(1)-cz12
     & (i1)%e(2)*cw56%e(2)-cz12(i1)%e(3)*cw56%e(3))
      do mu=0,3
      cwaux(i1,mu)=cz12(i1)%v*cw56%e(mu)+vwm(mu)*caux+cw56%v*cz1
     & 2(i1)%e(mu)
      end do
      end do
  
      do i1=1,2
      do m=0,3
        cw1256(i1)%e(m) = cw1256(i1)%e(m)+rcotw*cwaux(i1,m)
      enddo
      enddo
  
*** triple vertex -- pfz(mu)=p12(mu),pwm(mu)=p1256(mu),pwp(mu)=p56(mu),e
* fz=cf12(i1),ewm=#,ewp=cw56,res=cwaux(i1?,mu#)
      do mu=0,3
      vfz(mu)=p1256(mu)-p56(mu)
      vwm(mu)=p56(mu)-p12(mu)
      vwp(mu)=p12(mu)-p1256(mu)
      end do !mu
* vfz%efz
      do i1=1,2
* p.q -- p.q=cf12(i1)%v,p=cf12(i1)%e,q=vfz,bef=,aft=
      cf12(i1)%v=(cf12(i1)%e(0)*vfz(0)-cf12(i1)%e(1)*vfz(1)-cf12
     & (i1)%e(2)*vfz(2)-cf12(i1)%e(3)*vfz(3))
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewp
      do i1=1,2
* p.q -- p.q=caux,p=cf12(i1)%e,q=cw56%e,bef=,aft=
      caux=(cf12(i1)%e(0)*cw56%e(0)-cf12(i1)%e(1)*cw56%e(1)-cf12
     & (i1)%e(2)*cw56%e(2)-cf12(i1)%e(3)*cw56%e(3))
      do mu=0,3
      cwaux(i1,mu)=cf12(i1)%v*cw56%e(mu)+vwm(mu)*caux+cw56%v*cf1
     & 2(i1)%e(mu)
      end do
      end do
  
      do i1=1,2
      do m=0,3
        cw1256(i1)%e(m) = cw1256(i1)%e(m)+cwaux(i1,m)
      enddo
      enddo
  
  
      do i1=1,2
* pk0 -- p=cw1256(i1)%e
      cw1256(i1)%ek0=cw1256(i1)%e(0)-cw1256(i1)%e(1)
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw3456(i1)%e(m),a1=lw5_34(m)%a,c1=lw5_34(m)%c,a2=r6_34(i
* 1)%a,b2=r6_34(i1)%b,prq=p634q,bef=,aft=
      cw3456(i1)%e(m)=(lw5_34(m)%c(2)*p634q*r6_34(i1)%b(1)+lw5_3
     & 4(m)%a(2)*r6_34(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw3456(i1)%e(m),a1=l5_34(i1)%a,c1=l5_34(i1)%c,a2=rw6_34(
* m)%a,b2=rw6_34(m)%b,prq=p534q,bef=cw3456(i1)%e(m)+,aft=
      cw3456(i1)%e(m)=cw3456(i1)%e(m)+(l5_34(i1)%c(2)*p534q*rw6_
     & 34(m)%b(1)+l5_34(i1)%a(2)*rw6_34(m)%a(2))
      end do
      end do
  
      if (iup(id3).eq.1) then
      do m=0,3
* TLTR0_W -- aa=cw3456(2)%e(m),a1=lw3_56(m)%a,c1=lw3_56(m)%c,a2=r4_56%a,
* b2=r4_56%b,prq=p456q,bef=cw3456(2)%e(m)+,aft=
      cw3456(2)%e(m)=cw3456(2)%e(m)+(lw3_56(m)%c(2)*p456q*r4_56%
     & b(1)+lw3_56(m)%a(2)*r4_56%a(2))
      end do
      else
      do m=0,3
* TLTR0_W -- aa=cw3456(2)%e(m),a1=l3_56%a,c1=l3_56%c,a2=rw4_56(m)%a,b2=r
* w4_56(m)%b,prq=p356q,bef=cw3456(2)%e(m)+,aft=
      cw3456(2)%e(m)=cw3456(2)%e(m)+(l3_56%c(2)*p356q*rw4_56(m)%
     & b(1)+l3_56%a(2)*rw4_56(m)%a(2))
      end do
      endif
  
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p3456(mu),pwp(mu)=p56(mu),e
* fz=cz34(i1),ewm=#,ewp=cw56,res=cwaux(i1?,mu#)
      do mu=0,3
      vfz(mu)=p3456(mu)-p56(mu)
      vwm(mu)=p56(mu)-p34(mu)
      vwp(mu)=p34(mu)-p3456(mu)
      end do !mu
* vfz%efz
      do i1=1,2
* p.q -- p.q=cz34(i1)%v,p=cz34(i1)%e,q=vfz,bef=,aft=
      cz34(i1)%v=(cz34(i1)%e(0)*vfz(0)-cz34(i1)%e(1)*vfz(1)-cz34
     & (i1)%e(2)*vfz(2)-cz34(i1)%e(3)*vfz(3))
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewp
      do i1=1,2
* p.q -- p.q=caux,p=cz34(i1)%e,q=cw56%e,bef=,aft=
      caux=(cz34(i1)%e(0)*cw56%e(0)-cz34(i1)%e(1)*cw56%e(1)-cz34
     & (i1)%e(2)*cw56%e(2)-cz34(i1)%e(3)*cw56%e(3))
      do mu=0,3
      cwaux(i1,mu)=cz34(i1)%v*cw56%e(mu)+vwm(mu)*caux+cw56%v*cz3
     & 4(i1)%e(mu)
      end do
      end do
  
      do i1=1,2
      do m=0,3
        cw3456(i1)%e(m) = cw3456(i1)%e(m)+rcotw*cwaux(i1,m)
      enddo
      enddo
  
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p3456(mu),pwp(mu)=p56(mu),e
* fz=cf34(i1),ewm=#,ewp=cw56,res=cwaux(i1?,mu#)
      do mu=0,3
      vfz(mu)=p3456(mu)-p56(mu)
      vwm(mu)=p56(mu)-p34(mu)
      vwp(mu)=p34(mu)-p3456(mu)
      end do !mu
* vfz%efz
      do i1=1,2
* p.q -- p.q=cf34(i1)%v,p=cf34(i1)%e,q=vfz,bef=,aft=
      cf34(i1)%v=(cf34(i1)%e(0)*vfz(0)-cf34(i1)%e(1)*vfz(1)-cf34
     & (i1)%e(2)*vfz(2)-cf34(i1)%e(3)*vfz(3))
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewp
      do i1=1,2
* p.q -- p.q=caux,p=cf34(i1)%e,q=cw56%e,bef=,aft=
      caux=(cf34(i1)%e(0)*cw56%e(0)-cf34(i1)%e(1)*cw56%e(1)-cf34
     & (i1)%e(2)*cw56%e(2)-cf34(i1)%e(3)*cw56%e(3))
      do mu=0,3
      cwaux(i1,mu)=cf34(i1)%v*cw56%e(mu)+vwm(mu)*caux+cw56%v*cf3
     & 4(i1)%e(mu)
      end do
      end do
  
      do i1=1,2
      do m=0,3
        cw3456(i1)%e(m) = cw3456(i1)%e(m)+cwaux(i1,m)
      enddo
      enddo
  
  
      do i1=1,2
* pk0 -- p=cw3456(i1)%e
      cw3456(i1)%ek0=cw3456(i1)%e(0)-cw3456(i1)%e(1)
      end do
  
  
* (1278,3478) W-                                                        
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw1278(i1)%e(m),a1=lw7_12(m)%a,c1=lw7_12(m)%c,a2=r8_12(i
* 1)%a,b2=r8_12(i1)%b,prq=p812q,bef=,aft=
      cw1278(i1)%e(m)=(lw7_12(m)%c(2)*p812q*r8_12(i1)%b(1)+lw7_1
     & 2(m)%a(2)*r8_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw1278(i1)%e(m),a1=l7_12(i1)%a,c1=l7_12(i1)%c,a2=rw8_12(
* m)%a,b2=rw8_12(m)%b,prq=p712q,bef=cw1278(i1)%e(m)+,aft=
      cw1278(i1)%e(m)=cw1278(i1)%e(m)+(l7_12(i1)%c(2)*p712q*rw8_
     & 12(m)%b(1)+l7_12(i1)%a(2)*rw8_12(m)%a(2))
      end do
      end do
  
      if (iup(id1).eq.1) then
      do m=0,3
* TLTR0_W -- aa=cw1278(2)%e(m),a1=l1_78%a,c1=l1_78%c,a2=rw2_78(m)%a,b2=r
* w2_78(m)%b,prq=p178q,bef=cw1278(2)%e(m)+,aft=
      cw1278(2)%e(m)=cw1278(2)%e(m)+(l1_78%c(2)*p178q*rw2_78(m)%
     & b(1)+l1_78%a(2)*rw2_78(m)%a(2))
      end do
      else
      do m=0,3
* TLTR0_W -- aa=cw1278(2)%e(m),a1=lw1_78(m)%a,c1=lw1_78(m)%c,a2=r2_78%a,
* b2=r2_78%b,prq=p278q,bef=cw1278(2)%e(m)+,aft=
      cw1278(2)%e(m)=cw1278(2)%e(m)+(lw1_78(m)%c(2)*p278q*r2_78%
     & b(1)+lw1_78(m)%a(2)*r2_78%a(2))
      end do
      endif
  
*** triple vertex -- pfz(mu)=p12(mu),pwm(mu)=p78(mu),pwp(mu)=p1278(mu),e
* fz=cz12(i1),ewm=cw78,ewp=#,res=cwaux(i1?,mu#)
      do mu=0,3
      vfz(mu)=p78(mu)-p1278(mu)
      vwm(mu)=p1278(mu)-p12(mu)
      vwp(mu)=p12(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i1=1,2
* p.q -- p.q=cz12(i1)%v,p=cz12(i1)%e,q=vfz,bef=,aft=
      cz12(i1)%v=(cz12(i1)%e(0)*vfz(0)-cz12(i1)%e(1)*vfz(1)-cz12
     & (i1)%e(2)*vfz(2)-cz12(i1)%e(3)*vfz(3))
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* efz%ewm
      do i1=1,2
* p.q -- p.q=caux,p=cz12(i1)%e,q=cw78%e,bef=,aft=
      caux=(cz12(i1)%e(0)*cw78%e(0)-cz12(i1)%e(1)*cw78%e(1)-cz12
     & (i1)%e(2)*cw78%e(2)-cz12(i1)%e(3)*cw78%e(3))
      do mu=0,3
      cwaux(i1,mu)=cz12(i1)%v*cw78%e(mu)+cw78%v*cz12(i1)%e(mu)+v
     & wp(mu)*caux
      end do
      end do
  
      do i1=1,2
      do m=0,3
        cw1278(i1)%e(m) = cw1278(i1)%e(m)+rcotw*cwaux(i1,m)
      enddo
      enddo
*** triple vertex -- pfz(mu)=p12(mu),pwm(mu)=p78(mu),pwp(mu)=p1278(mu),e
* fz=cf12(i1),ewm=cw78,ewp=#,res=cwaux(i1?,mu#)
      do mu=0,3
      vfz(mu)=p78(mu)-p1278(mu)
      vwm(mu)=p1278(mu)-p12(mu)
      vwp(mu)=p12(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i1=1,2
* p.q -- p.q=cf12(i1)%v,p=cf12(i1)%e,q=vfz,bef=,aft=
      cf12(i1)%v=(cf12(i1)%e(0)*vfz(0)-cf12(i1)%e(1)*vfz(1)-cf12
     & (i1)%e(2)*vfz(2)-cf12(i1)%e(3)*vfz(3))
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* efz%ewm
      do i1=1,2
* p.q -- p.q=caux,p=cf12(i1)%e,q=cw78%e,bef=,aft=
      caux=(cf12(i1)%e(0)*cw78%e(0)-cf12(i1)%e(1)*cw78%e(1)-cf12
     & (i1)%e(2)*cw78%e(2)-cf12(i1)%e(3)*cw78%e(3))
      do mu=0,3
      cwaux(i1,mu)=cf12(i1)%v*cw78%e(mu)+cw78%v*cf12(i1)%e(mu)+v
     & wp(mu)*caux
      end do
      end do
  
      do i1=1,2
      do m=0,3
        cw1278(i1)%e(m) = cw1278(i1)%e(m)+cwaux(i1,m)
      enddo
      enddo
  
      do i1=1,2
* pk0 -- p=cw1278(i1)%e
      cw1278(i1)%ek0=cw1278(i1)%e(0)-cw1278(i1)%e(1)
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw3478(i1)%e(m),a1=lw7_34(m)%a,c1=lw7_34(m)%c,a2=r8_34(i
* 1)%a,b2=r8_34(i1)%b,prq=p834q,bef=,aft=
      cw3478(i1)%e(m)=(lw7_34(m)%c(2)*p834q*r8_34(i1)%b(1)+lw7_3
     & 4(m)%a(2)*r8_34(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cw3478(i1)%e(m),a1=l7_34(i1)%a,c1=l7_34(i1)%c,a2=rw8_34(
* m)%a,b2=rw8_34(m)%b,prq=p734q,bef=cw3478(i1)%e(m)+,aft=
      cw3478(i1)%e(m)=cw3478(i1)%e(m)+(l7_34(i1)%c(2)*p734q*rw8_
     & 34(m)%b(1)+l7_34(i1)%a(2)*rw8_34(m)%a(2))
      end do
      end do
  
      if (iup(id3).eq.1) then
      do m=0,3
* TLTR0_W -- aa=cw3478(2)%e(m),a1=l3_78%a,c1=l3_78%c,a2=rw4_78(m)%a,b2=r
* w4_78(m)%b,prq=p378q,bef=cw3478(2)%e(m)+,aft=
      cw3478(2)%e(m)=cw3478(2)%e(m)+(l3_78%c(2)*p378q*rw4_78(m)%
     & b(1)+l3_78%a(2)*rw4_78(m)%a(2))
      end do
      else
      do m=0,3
* TLTR0_W -- aa=cw3478(2)%e(m),a1=lw3_78(m)%a,c1=lw3_78(m)%c,a2=r4_78%a,
* b2=r4_78%b,prq=p478q,bef=cw3478(2)%e(m)+,aft=
      cw3478(2)%e(m)=cw3478(2)%e(m)+(lw3_78(m)%c(2)*p478q*r4_78%
     & b(1)+lw3_78(m)%a(2)*r4_78%a(2))
      end do
      endif
  
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p78(mu),pwp(mu)=p3478(mu),e
* fz=cz34(i1),ewm=cw78,ewp=#,res=cwaux(i1?,mu#)
      do mu=0,3
      vfz(mu)=p78(mu)-p3478(mu)
      vwm(mu)=p3478(mu)-p34(mu)
      vwp(mu)=p34(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i1=1,2
* p.q -- p.q=cz34(i1)%v,p=cz34(i1)%e,q=vfz,bef=,aft=
      cz34(i1)%v=(cz34(i1)%e(0)*vfz(0)-cz34(i1)%e(1)*vfz(1)-cz34
     & (i1)%e(2)*vfz(2)-cz34(i1)%e(3)*vfz(3))
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* efz%ewm
      do i1=1,2
* p.q -- p.q=caux,p=cz34(i1)%e,q=cw78%e,bef=,aft=
      caux=(cz34(i1)%e(0)*cw78%e(0)-cz34(i1)%e(1)*cw78%e(1)-cz34
     & (i1)%e(2)*cw78%e(2)-cz34(i1)%e(3)*cw78%e(3))
      do mu=0,3
      cwaux(i1,mu)=cz34(i1)%v*cw78%e(mu)+cw78%v*cz34(i1)%e(mu)+v
     & wp(mu)*caux
      end do
      end do
  
      do i1=1,2
      do m=0,3
        cw3478(i1)%e(m) = cw3478(i1)%e(m)+rcotw*cwaux(i1,m)
      enddo
      enddo
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p78(mu),pwp(mu)=p3478(mu),e
* fz=cf34(i1),ewm=cw78,ewp=#,res=cwaux(i1?,mu#)
      do mu=0,3
      vfz(mu)=p78(mu)-p3478(mu)
      vwm(mu)=p3478(mu)-p34(mu)
      vwp(mu)=p34(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i1=1,2
* p.q -- p.q=cf34(i1)%v,p=cf34(i1)%e,q=vfz,bef=,aft=
      cf34(i1)%v=(cf34(i1)%e(0)*vfz(0)-cf34(i1)%e(1)*vfz(1)-cf34
     & (i1)%e(2)*vfz(2)-cf34(i1)%e(3)*vfz(3))
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* efz%ewm
      do i1=1,2
* p.q -- p.q=caux,p=cf34(i1)%e,q=cw78%e,bef=,aft=
      caux=(cf34(i1)%e(0)*cw78%e(0)-cf34(i1)%e(1)*cw78%e(1)-cf34
     & (i1)%e(2)*cw78%e(2)-cf34(i1)%e(3)*cw78%e(3))
      do mu=0,3
      cwaux(i1,mu)=cf34(i1)%v*cw78%e(mu)+cw78%v*cf34(i1)%e(mu)+v
     & wp(mu)*caux
      end do
      end do
  
      do i1=1,2
      do m=0,3
        cw3478(i1)%e(m) = cw3478(i1)%e(m)+cwaux(i1,m)
      enddo
      enddo
  
      do i1=1,2
* pk0 -- p=cw3478(i1)%e
      cw3478(i1)%ek0=cw3478(i1)%e(0)-cw3478(i1)%e(1)
      end do
  
  
**QCD Z,f connection                                                    
  
      if (ilept(id1).ne.1.and.ilept(id3).ne.1) then
      do i3=1,2
      do m=0,3
* TLTR0 -- aa=cgz1234(&,i3)%e(m),a1=lz1_34(m)%a,c1=lz1_34(m)%c,a2=rg2_34
* (i3)%a,b2=rg2_34(i3)%b,prq=p234q,bef=,aft=
      cgz1234(1,i3)%e(m)=(lz1_34(m)%a(1)*rg2_34(i3)%a(1)+lz1_34(
     & m)%c(1)*p234q*rg2_34(i3)%b(2))
      cgz1234(2,i3)%e(m)=(lz1_34(m)%c(2)*p234q*rg2_34(i3)%b(1)+l
     & z1_34(m)%a(2)*rg2_34(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do m=0,3
* TLTR0 -- aa=cgz1234(&,i3)%e(m),a1=lg1_34(i3)%a,c1=lg1_34(i3)%c,a2=rz2_
* 34(m)%a,b2=rz2_34(m)%b,prq=p134q,bef=cgz1234(&,i3)%e(m)+,aft=
      cgz1234(1,i3)%e(m)=cgz1234(1,i3)%e(m)+(lg1_34(i3)%a(1)*rz2
     & _34(m)%a(1)+lg1_34(i3)%c(1)*p134q*rz2_34(m)%b(2))
      cgz1234(2,i3)%e(m)=cgz1234(2,i3)%e(m)+(lg1_34(i3)%c(2)*p13
     & 4q*rz2_34(m)%b(1)+lg1_34(i3)%a(2)*rz2_34(m)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cgz1234(i1,&)%e(m),a1=lz3_12(m)%a,c1=lz3_12(m)%c,a2=rg4_12
* (i1)%a,b2=rg4_12(i1)%b,prq=p412q,bef=cgz1234(i1,&)%e(m)+,aft=
      cgz1234(i1,1)%e(m)=cgz1234(i1,1)%e(m)+(lz3_12(m)%a(1)*rg4_
     & 12(i1)%a(1)+lz3_12(m)%c(1)*p412q*rg4_12(i1)%b(2))
      cgz1234(i1,2)%e(m)=cgz1234(i1,2)%e(m)+(lz3_12(m)%c(2)*p412
     & q*rg4_12(i1)%b(1)+lz3_12(m)%a(2)*rg4_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cgz1234(i1,&)%e(m),a1=lg3_12(i1)%a,c1=lg3_12(i1)%c,a2=rz4_
* 12(m)%a,b2=rz4_12(m)%b,prq=p312q,bef=cgz1234(i1,&)%e(m)+,aft=
      cgz1234(i1,1)%e(m)=cgz1234(i1,1)%e(m)+(lg3_12(i1)%a(1)*rz4
     & _12(m)%a(1)+lg3_12(i1)%c(1)*p312q*rz4_12(m)%b(2))
      cgz1234(i1,2)%e(m)=cgz1234(i1,2)%e(m)+(lg3_12(i1)%c(2)*p31
     & 2q*rz4_12(m)%b(1)+lg3_12(i1)%a(2)*rz4_12(m)%a(2))
      end do
      end do
      do i3=1,2
      do m=0,3
* TLTR0 -- aa=cgf1234(&,i3)%e(m),a1=lf1_34(m)%a,c1=lf1_34(m)%c,a2=rg2_34
* (i3)%a,b2=rg2_34(i3)%b,prq=p234q,bef=,aft=
      cgf1234(1,i3)%e(m)=(lf1_34(m)%a(1)*rg2_34(i3)%a(1)+lf1_34(
     & m)%c(1)*p234q*rg2_34(i3)%b(2))
      cgf1234(2,i3)%e(m)=(lf1_34(m)%c(2)*p234q*rg2_34(i3)%b(1)+l
     & f1_34(m)%a(2)*rg2_34(i3)%a(2))
      end do
      end do
  
      do i3=1,2
      do m=0,3
* TLTR0 -- aa=cgf1234(&,i3)%e(m),a1=lg1_34(i3)%a,c1=lg1_34(i3)%c,a2=rf2_
* 34(m)%a,b2=rf2_34(m)%b,prq=p134q,bef=cgf1234(&,i3)%e(m)+,aft=
      cgf1234(1,i3)%e(m)=cgf1234(1,i3)%e(m)+(lg1_34(i3)%a(1)*rf2
     & _34(m)%a(1)+lg1_34(i3)%c(1)*p134q*rf2_34(m)%b(2))
      cgf1234(2,i3)%e(m)=cgf1234(2,i3)%e(m)+(lg1_34(i3)%c(2)*p13
     & 4q*rf2_34(m)%b(1)+lg1_34(i3)%a(2)*rf2_34(m)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cgf1234(i1,&)%e(m),a1=lf3_12(m)%a,c1=lf3_12(m)%c,a2=rg4_12
* (i1)%a,b2=rg4_12(i1)%b,prq=p412q,bef=cgf1234(i1,&)%e(m)+,aft=
      cgf1234(i1,1)%e(m)=cgf1234(i1,1)%e(m)+(lf3_12(m)%a(1)*rg4_
     & 12(i1)%a(1)+lf3_12(m)%c(1)*p412q*rg4_12(i1)%b(2))
      cgf1234(i1,2)%e(m)=cgf1234(i1,2)%e(m)+(lf3_12(m)%c(2)*p412
     & q*rg4_12(i1)%b(1)+lf3_12(m)%a(2)*rg4_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0 -- aa=cgf1234(i1,&)%e(m),a1=lg3_12(i1)%a,c1=lg3_12(i1)%c,a2=rf4_
* 12(m)%a,b2=rf4_12(m)%b,prq=p312q,bef=cgf1234(i1,&)%e(m)+,aft=
      cgf1234(i1,1)%e(m)=cgf1234(i1,1)%e(m)+(lg3_12(i1)%a(1)*rf4
     & _12(m)%a(1)+lg3_12(i1)%c(1)*p312q*rf4_12(m)%b(2))
      cgf1234(i1,2)%e(m)=cgf1234(i1,2)%e(m)+(lg3_12(i1)%c(2)*p31
     & 2q*rf4_12(m)%b(1)+lg3_12(i1)%a(2)*rf4_12(m)%a(2))
      end do
      end do
  
      endif
  
**QCD W connections                                                     
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cgw5612(i1)%e(m),a1=lw5_12(m)%a,c1=lw5_12(m)%c,a2=rg6_12
* (i1)%a,b2=rg6_12(i1)%b,prq=p612q,bef=,aft=
      cgw5612(i1)%e(m)=(lw5_12(m)%c(2)*p612q*rg6_12(i1)%b(1)+lw5
     & _12(m)%a(2)*rg6_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cgw5612(i1)%e(m),a1=lg5_12(i1)%a,c1=lg5_12(i1)%c,a2=rw6_
* 12(m)%a,b2=rw6_12(m)%b,prq=p512q,bef=cgw5612(i1)%e(m)+,aft=
      cgw5612(i1)%e(m)=cgw5612(i1)%e(m)+(lg5_12(i1)%c(2)*p512q*r
     & w6_12(m)%b(1)+lg5_12(i1)%a(2)*rw6_12(m)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cgw7812(i1)%e(m),a1=lw7_12(m)%a,c1=lw7_12(m)%c,a2=rg8_12
* (i1)%a,b2=rg8_12(i1)%b,prq=p812q,bef=,aft=
      cgw7812(i1)%e(m)=(lw7_12(m)%c(2)*p812q*rg8_12(i1)%b(1)+lw7
     & _12(m)%a(2)*rg8_12(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cgw7812(i1)%e(m),a1=lg7_12(i1)%a,c1=lg7_12(i1)%c,a2=rw8_
* 12(m)%a,b2=rw8_12(m)%b,prq=p712q,bef=cgw7812(i1)%e(m)+,aft=
      cgw7812(i1)%e(m)=cgw7812(i1)%e(m)+(lg7_12(i1)%c(2)*p712q*r
     & w8_12(m)%b(1)+lg7_12(i1)%a(2)*rw8_12(m)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cgw5634(i1)%e(m),a1=lw5_34(m)%a,c1=lw5_34(m)%c,a2=rg6_34
* (i1)%a,b2=rg6_34(i1)%b,prq=p634q,bef=,aft=
      cgw5634(i1)%e(m)=(lw5_34(m)%c(2)*p634q*rg6_34(i1)%b(1)+lw5
     & _34(m)%a(2)*rg6_34(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cgw5634(i1)%e(m),a1=lg5_34(i1)%a,c1=lg5_34(i1)%c,a2=rw6_
* 34(m)%a,b2=rw6_34(m)%b,prq=p534q,bef=cgw5634(i1)%e(m)+,aft=
      cgw5634(i1)%e(m)=cgw5634(i1)%e(m)+(lg5_34(i1)%c(2)*p534q*r
     & w6_34(m)%b(1)+lg5_34(i1)%a(2)*rw6_34(m)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cgw7834(i1)%e(m),a1=lw7_34(m)%a,c1=lw7_34(m)%c,a2=rg8_34
* (i1)%a,b2=rg8_34(i1)%b,prq=p834q,bef=,aft=
      cgw7834(i1)%e(m)=(lw7_34(m)%c(2)*p834q*rg8_34(i1)%b(1)+lw7
     & _34(m)%a(2)*rg8_34(i1)%a(2))
      end do
      end do
  
      do i1=1,2
      do m=0,3
* TLTR0_W -- aa=cgw7834(i1)%e(m),a1=lg7_34(i1)%a,c1=lg7_34(i1)%c,a2=rw8_
* 34(m)%a,b2=rw8_34(m)%b,prq=p734q,bef=cgw7834(i1)%e(m)+,aft=
      cgw7834(i1)%e(m)=cgw7834(i1)%e(m)+(lg7_34(i1)%c(2)*p734q*r
     & w8_34(m)%b(1)+lg7_34(i1)%a(2)*rw8_34(m)%a(2))
      end do
      end do
  
  
*                                                                       
*  BOSON CONNECTION                                                     
*  four->four                                                           
*                                                                       
*      dimension c1234(2,2),c1256(2,2),c1278(2,2)                       
  
* -->   1234 -> 5678                                                    
  
**EW                                                                    
*   Z,f                                                                 
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux(i1,i3),p=cz1234(i1,i3)%e,q=cz5678%e,bef=,aft=
      c4aux(i1,i3)=(cz1234(i1,i3)%e(0)*cz5678%e(0)-cz1234(i1,i3)
     & %e(1)*cz5678%e(1)-cz1234(i1,i3)%e(2)*cz5678%e(2)-cz1234(i
     & 1,i3)%e(3)*cz5678%e(3))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux_2(i1,i3),p=cf1234(i1,i3)%e,q=cf5678%e,bef=,aft=
      c4aux_2(i1,i3)=(cf1234(i1,i3)%e(0)*cf5678%e(0)-cf1234(i1,i
     & 3)%e(1)*cf5678%e(1)-cf1234(i1,i3)%e(2)*cf5678%e(2)-cf1234
     & (i1,i3)%e(3)*cf5678%e(3))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
  
       cbct1234(i1,i3)= c4aux(i1,i3)/(-p1234q+cmz2)
     &    + c4aux_2(i1,i3)/(-p1234q)
      enddo
      enddo
  
*   higgs                                                               
  
      if (rmh.ge.0.d0) then
  
       do i1=1,2
       do i3=1,2
        cbct1234(i1,i3)= cbct1234(i1,i3)+
     &    ch1234(i1,i3)*ch5678/(p1234q-cmh2)
       enddo
       enddo
  
      endif
  
  
* --> 1256 ->347\8\                                                     
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux(i1,i3),p=cw1256(i1)%e,q=cw3478(i3)%e,bef=,aft=
      c4aux(i1,i3)=(cw1256(i1)%e(0)*cw3478(i3)%e(0)-cw1256(i1)%e
     & (1)*cw3478(i3)%e(1)-cw1256(i1)%e(2)*cw3478(i3)%e(2)-cw125
     & 6(i1)%e(3)*cw3478(i3)%e(3))
      end do
      end do
  
      do i1=1,2
* p.q -- p.q=c2aux1(i1),p=cw1256(i1)%e,q=p1256,bef=,aft=
      c2aux1(i1)=(cw1256(i1)%e(0)*p1256(0)-cw1256(i1)%e(1)*p1256
     & (1)-cw1256(i1)%e(2)*p1256(2)-cw1256(i1)%e(3)*p1256(3))
      end do
  
      do i3=1,2
* p.q -- p.q=c2aux2(i3),p=cw3478(i3)%e,q=p1256,bef=,aft=
      c2aux2(i3)=(cw3478(i3)%e(0)*p1256(0)-cw3478(i3)%e(1)*p1256
     & (1)-cw3478(i3)%e(2)*p1256(2)-cw3478(i3)%e(3)*p1256(3))
      end do
  
      do i1=1,2
      do i3=1,2
       cbct1256(i1,i3)= (c4aux(i1,i3)-c2aux1(i1)*c2aux2(i3)/cmw2)
     &     /(-p1256q+cmw2)
      enddo
      enddo
  
  
* --> 1278 ->347\8\                                                     
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux(i1,i3),p=cw1278(i1)%e,q=cw3456(i3)%e,bef=,aft=
      c4aux(i1,i3)=(cw1278(i1)%e(0)*cw3456(i3)%e(0)-cw1278(i1)%e
     & (1)*cw3456(i3)%e(1)-cw1278(i1)%e(2)*cw3456(i3)%e(2)-cw127
     & 8(i1)%e(3)*cw3456(i3)%e(3))
      end do
      end do
  
      do i1=1,2
* p.q -- p.q=c2aux1(i1),p=cw1278(i1)%e,q=p1278,bef=,aft=
      c2aux1(i1)=(cw1278(i1)%e(0)*p1278(0)-cw1278(i1)%e(1)*p1278
     & (1)-cw1278(i1)%e(2)*p1278(2)-cw1278(i1)%e(3)*p1278(3))
      end do
  
      do i3=1,2
* p.q -- p.q=c2aux2(i3),p=cw3456(i3)%e,q=p1278,bef=,aft=
      c2aux2(i3)=(cw3456(i3)%e(0)*p1278(0)-cw3456(i3)%e(1)*p1278
     & (1)-cw3456(i3)%e(2)*p1278(2)-cw3456(i3)%e(3)*p1278(3))
      end do
  
      do i1=1,2
      do i3=1,2
       cbct1278(i1,i3)= (c4aux(i1,i3)-c2aux1(i1)*c2aux2(i3)/cmw2)
     &     /(-p1278q+cmw2)
      enddo
      enddo
  
  
**QCD                                                                   
  
* -->   1234 -> 5678                                                    
  
*   Z,f 1~3                                                             
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux(i1,i3),p=cgz1234(i1,i3)%e,q=cz5678%e,bef=,aft=
      c4aux(i1,i3)=(cgz1234(i1,i3)%e(0)*cz5678%e(0)-cgz1234(i1,i
     & 3)%e(1)*cz5678%e(1)-cgz1234(i1,i3)%e(2)*cz5678%e(2)-cgz12
     & 34(i1,i3)%e(3)*cz5678%e(3))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux_2(i1,i3),p=cgf1234(i1,i3)%e,q=cf5678%e,bef=,aft=
      c4aux_2(i1,i3)=(cgf1234(i1,i3)%e(0)*cf5678%e(0)-cgf1234(i1
     & ,i3)%e(1)*cf5678%e(1)-cgf1234(i1,i3)%e(2)*cf5678%e(2)-cgf
     & 1234(i1,i3)%e(3)*cf5678%e(3))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
       cgbct1234(i1,i3)= c4aux(i1,i3)/(-p1234q+cmz2)
     &    + c4aux_2(i1,i3)/(-p1234q)
      enddo
      enddo
  
* gluon connection 34\ 12\~56\ 78\                                      
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=cgct1256(i1,i3),p=cg1234(i1,i3)%e,q=cg5678%e,bef=,aft=/(-p1
* 234q)
      cgct1256(i1,i3)=(cg1234(i1,i3)%e(0)*cg5678%e(0)-cg1234(i1,
     & i3)%e(1)*cg5678%e(1)-cg1234(i1,i3)%e(2)*cg5678%e(2)-cg123
     & 4(i1,i3)%e(3)*cg5678%e(3))/(-p1234q)
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=cgct1278(i1,i3),p=cg1234(i1,i3)%e,q=cg7856%e,bef=,aft=/(-p1
* 234q)
      cgct1278(i1,i3)=(cg1234(i1,i3)%e(0)*cg7856%e(0)-cg1234(i1,
     & i3)%e(1)*cg7856%e(1)-cg1234(i1,i3)%e(2)*cg7856%e(2)-cg123
     & 4(i1,i3)%e(3)*cg7856%e(3))/(-p1234q)
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=cgct3456(i1,i3),p=cg3412(i3,i1)%e,q=cg5678%e,bef=,aft=/(-p1
* 234q)
      cgct3456(i1,i3)=(cg3412(i3,i1)%e(0)*cg5678%e(0)-cg3412(i3,
     & i1)%e(1)*cg5678%e(1)-cg3412(i3,i1)%e(2)*cg5678%e(2)-cg341
     & 2(i3,i1)%e(3)*cg5678%e(3))/(-p1234q)
      end do
      end do
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=cgct3478(i1,i3),p=cg3412(i3,i1)%e,q=cg7856%e,bef=,aft=/(-p1
* 234q)
      cgct3478(i1,i3)=(cg3412(i3,i1)%e(0)*cg7856%e(0)-cg3412(i3,
     & i1)%e(1)*cg7856%e(1)-cg3412(i3,i1)%e(2)*cg7856%e(2)-cg341
     & 2(i3,i1)%e(3)*cg7856%e(3))/(-p1234q)
      end do
      end do
  
  
* W connection  12\~56\ ->34\7\8\                                       
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux(i1,i3),p=cgw5612(i1)%e,q=cw3478(i3)%e,bef=,aft=
      c4aux(i1,i3)=(cgw5612(i1)%e(0)*cw3478(i3)%e(0)-cgw5612(i1)
     & %e(1)*cw3478(i3)%e(1)-cgw5612(i1)%e(2)*cw3478(i3)%e(2)-cg
     & w5612(i1)%e(3)*cw3478(i3)%e(3))
      end do
      end do
  
      do i1=1,2
* p.q -- p.q=c2aux1(i1),p=cgw5612(i1)%e,q=p1256,bef=,aft=
      c2aux1(i1)=(cgw5612(i1)%e(0)*p1256(0)-cgw5612(i1)%e(1)*p12
     & 56(1)-cgw5612(i1)%e(2)*p1256(2)-cgw5612(i1)%e(3)*p1256(3)
     & )
      end do
  
      do i3=1,2
* p.q -- p.q=c2aux2(i3),p=cw3478(i3)%e,q=p1256,bef=,aft=
      c2aux2(i3)=(cw3478(i3)%e(0)*p1256(0)-cw3478(i3)%e(1)*p1256
     & (1)-cw3478(i3)%e(2)*p1256(2)-cw3478(i3)%e(3)*p1256(3))
      end do
  
      do i1=1,2
      do i3=1,2
       cgbct1256(i1,i3)=(c4aux(i1,i3)-c2aux1(i1)*c2aux2(i3)/cmw2)
     &     /(-p1256q+cmw2)
      enddo
      enddo
  
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux(i1,i3),p=cgw7812(i1)%e,q=cw3456(i3)%e,bef=,aft=
      c4aux(i1,i3)=(cgw7812(i1)%e(0)*cw3456(i3)%e(0)-cgw7812(i1)
     & %e(1)*cw3456(i3)%e(1)-cgw7812(i1)%e(2)*cw3456(i3)%e(2)-cg
     & w7812(i1)%e(3)*cw3456(i3)%e(3))
      end do
      end do
  
      do i1=1,2
* p.q -- p.q=c2aux1(i1),p=cgw7812(i1)%e,q=p1278,bef=,aft=
      c2aux1(i1)=(cgw7812(i1)%e(0)*p1278(0)-cgw7812(i1)%e(1)*p12
     & 78(1)-cgw7812(i1)%e(2)*p1278(2)-cgw7812(i1)%e(3)*p1278(3)
     & )
      end do
  
      do i3=1,2
* p.q -- p.q=c2aux2(i3),p=cw3456(i3)%e,q=p1278,bef=,aft=
      c2aux2(i3)=(cw3456(i3)%e(0)*p1278(0)-cw3456(i3)%e(1)*p1278
     & (1)-cw3456(i3)%e(2)*p1278(2)-cw3456(i3)%e(3)*p1278(3))
      end do
  
      do i1=1,2
      do i3=1,2
       cgbct1278(i1,i3)=(c4aux(i1,i3)-c2aux1(i1)*c2aux2(i3)/cmw2)
     &     /(-p1278q+cmw2)
      enddo
      enddo
  
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux(i1,i3),p=cgw5634(i3)%e,q=cw1278(i1)%e,bef=,aft=
      c4aux(i1,i3)=(cgw5634(i3)%e(0)*cw1278(i1)%e(0)-cgw5634(i3)
     & %e(1)*cw1278(i1)%e(1)-cgw5634(i3)%e(2)*cw1278(i1)%e(2)-cg
     & w5634(i3)%e(3)*cw1278(i1)%e(3))
      end do
      end do
  
      do i1=1,2
* p.q -- p.q=c2aux1(i1),p=cgw5634(i1)%e,q=p3456,bef=,aft=
      c2aux1(i1)=(cgw5634(i1)%e(0)*p3456(0)-cgw5634(i1)%e(1)*p34
     & 56(1)-cgw5634(i1)%e(2)*p3456(2)-cgw5634(i1)%e(3)*p3456(3)
     & )
      end do
  
      do i3=1,2
* p.q -- p.q=c2aux2(i3),p=cw1278(i3)%e,q=p3456,bef=,aft=
      c2aux2(i3)=(cw1278(i3)%e(0)*p3456(0)-cw1278(i3)%e(1)*p3456
     & (1)-cw1278(i3)%e(2)*p3456(2)-cw1278(i3)%e(3)*p3456(3))
      end do
  
      do i1=1,2
      do i3=1,2
       cgbct3456(i1,i3)=(c4aux(i1,i3)-c2aux1(i3)*c2aux2(i1)/cmw2)
     &     /(-p3456q+cmw2)
      enddo
      enddo
  
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=c4aux(i1,i3),p=cgw7834(i3)%e,q=cw1256(i1)%e,bef=,aft=
      c4aux(i1,i3)=(cgw7834(i3)%e(0)*cw1256(i1)%e(0)-cgw7834(i3)
     & %e(1)*cw1256(i1)%e(1)-cgw7834(i3)%e(2)*cw1256(i1)%e(2)-cg
     & w7834(i3)%e(3)*cw1256(i1)%e(3))
      end do
      end do
  
      do i1=1,2
* p.q -- p.q=c2aux1(i1),p=cgw7834(i1)%e,q=p3478,bef=,aft=
      c2aux1(i1)=(cgw7834(i1)%e(0)*p3478(0)-cgw7834(i1)%e(1)*p34
     & 78(1)-cgw7834(i1)%e(2)*p3478(2)-cgw7834(i1)%e(3)*p3478(3)
     & )
      end do
  
      do i3=1,2
* p.q -- p.q=c2aux2(i3),p=cw1256(i3)%e,q=p3478,bef=,aft=
      c2aux2(i3)=(cw1256(i3)%e(0)*p3478(0)-cw1256(i3)%e(1)*p3478
     & (1)-cw1256(i3)%e(2)*p3478(2)-cw1256(i3)%e(3)*p3478(3))
      end do
  
      do i1=1,2
      do i3=1,2
       cgbct3478(i1,i3)=(c4aux(i1,i3)-c2aux1(i3)*c2aux2(i1)/cmw2)
     &     /(-p3478q+cmw2)
      enddo
      enddo
  
  
  
*                                                                       
* QUARTIC COUPLINGS                                                     
*                                                                       
  
      do i1=1,2
      do i3=1,2
        cquart(i1,i3)=czero
      enddo
      enddo
  
  
* p.q -- p.q=cwwaux,p=cw56%e,q=cw78%e,bef=,aft=
      cwwaux=(cw56%e(0)*cw78%e(0)-cw56%e(1)*cw78%e(1)-cw56%e(2)*
     & cw78%e(2)-cw56%e(3)*cw78%e(3))
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=czzaux,p=cz12(i1)%e,q=cz34(i3)%e,bef=,aft=
      czzaux=(cz12(i1)%e(0)*cz34(i3)%e(0)-cz12(i1)%e(1)*cz34(i3)
     & %e(1)-cz12(i1)%e(2)*cz34(i3)%e(2)-cz12(i1)%e(3)*cz34(i3)%
     & e(3))
  
* p.q -- p.q=cz1wpaux,p=cz12(i1)%e,q=cw56%e,bef=,aft=
      cz1wpaux=(cz12(i1)%e(0)*cw56%e(0)-cz12(i1)%e(1)*cw56%e(1)-
     & cz12(i1)%e(2)*cw56%e(2)-cz12(i1)%e(3)*cw56%e(3))
  
* p.q -- p.q=cz1wmaux,p=cz12(i1)%e,q=cw78%e,bef=,aft=
      cz1wmaux=(cz12(i1)%e(0)*cw78%e(0)-cz12(i1)%e(1)*cw78%e(1)-
     & cz12(i1)%e(2)*cw78%e(2)-cz12(i1)%e(3)*cw78%e(3))
  
* p.q -- p.q=cz2wpaux,p=cz34(i3)%e,q=cw56%e,bef=,aft=
      cz2wpaux=(cz34(i3)%e(0)*cw56%e(0)-cz34(i3)%e(1)*cw56%e(1)-
     & cz34(i3)%e(2)*cw56%e(2)-cz34(i3)%e(3)*cw56%e(3))
  
* p.q -- p.q=cz2wmaux,p=cz34(i3)%e,q=cw78%e,bef=,aft=
      cz2wmaux=(cz34(i3)%e(0)*cw78%e(0)-cz34(i3)%e(1)*cw78%e(1)-
     & cz34(i3)%e(2)*cw78%e(2)-cz34(i3)%e(3)*cw78%e(3))
  
      cquart(i1,i3)=cquart(i1,i3)+rcotw**2*(2*czzaux*cwwaux -
     &   cz1wpaux*cz2wmaux - cz1wmaux*cz2wpaux)
  
      enddo
      enddo
  
* p.q -- p.q=cwwaux,p=cw56%e,q=cw78%e,bef=,aft=
      cwwaux=(cw56%e(0)*cw78%e(0)-cw56%e(1)*cw78%e(1)-cw56%e(2)*
     & cw78%e(2)-cw56%e(3)*cw78%e(3))
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=czzaux,p=cz12(i1)%e,q=cf34(i3)%e,bef=,aft=
      czzaux=(cz12(i1)%e(0)*cf34(i3)%e(0)-cz12(i1)%e(1)*cf34(i3)
     & %e(1)-cz12(i1)%e(2)*cf34(i3)%e(2)-cz12(i1)%e(3)*cf34(i3)%
     & e(3))
  
* p.q -- p.q=cz1wpaux,p=cz12(i1)%e,q=cw56%e,bef=,aft=
      cz1wpaux=(cz12(i1)%e(0)*cw56%e(0)-cz12(i1)%e(1)*cw56%e(1)-
     & cz12(i1)%e(2)*cw56%e(2)-cz12(i1)%e(3)*cw56%e(3))
  
* p.q -- p.q=cz1wmaux,p=cz12(i1)%e,q=cw78%e,bef=,aft=
      cz1wmaux=(cz12(i1)%e(0)*cw78%e(0)-cz12(i1)%e(1)*cw78%e(1)-
     & cz12(i1)%e(2)*cw78%e(2)-cz12(i1)%e(3)*cw78%e(3))
  
* p.q -- p.q=cz2wpaux,p=cf34(i3)%e,q=cw56%e,bef=,aft=
      cz2wpaux=(cf34(i3)%e(0)*cw56%e(0)-cf34(i3)%e(1)*cw56%e(1)-
     & cf34(i3)%e(2)*cw56%e(2)-cf34(i3)%e(3)*cw56%e(3))
  
* p.q -- p.q=cz2wmaux,p=cf34(i3)%e,q=cw78%e,bef=,aft=
      cz2wmaux=(cf34(i3)%e(0)*cw78%e(0)-cf34(i3)%e(1)*cw78%e(1)-
     & cf34(i3)%e(2)*cw78%e(2)-cf34(i3)%e(3)*cw78%e(3))
  
      cquart(i1,i3)=cquart(i1,i3)+rcotw*(2*czzaux*cwwaux -
     &   cz1wpaux*cz2wmaux - cz1wmaux*cz2wpaux)
  
      enddo
      enddo
  
* p.q -- p.q=cwwaux,p=cw56%e,q=cw78%e,bef=,aft=
      cwwaux=(cw56%e(0)*cw78%e(0)-cw56%e(1)*cw78%e(1)-cw56%e(2)*
     & cw78%e(2)-cw56%e(3)*cw78%e(3))
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=czzaux,p=cf12(i1)%e,q=cz34(i3)%e,bef=,aft=
      czzaux=(cf12(i1)%e(0)*cz34(i3)%e(0)-cf12(i1)%e(1)*cz34(i3)
     & %e(1)-cf12(i1)%e(2)*cz34(i3)%e(2)-cf12(i1)%e(3)*cz34(i3)%
     & e(3))
  
* p.q -- p.q=cz1wpaux,p=cf12(i1)%e,q=cw56%e,bef=,aft=
      cz1wpaux=(cf12(i1)%e(0)*cw56%e(0)-cf12(i1)%e(1)*cw56%e(1)-
     & cf12(i1)%e(2)*cw56%e(2)-cf12(i1)%e(3)*cw56%e(3))
  
* p.q -- p.q=cz1wmaux,p=cf12(i1)%e,q=cw78%e,bef=,aft=
      cz1wmaux=(cf12(i1)%e(0)*cw78%e(0)-cf12(i1)%e(1)*cw78%e(1)-
     & cf12(i1)%e(2)*cw78%e(2)-cf12(i1)%e(3)*cw78%e(3))
  
* p.q -- p.q=cz2wpaux,p=cz34(i3)%e,q=cw56%e,bef=,aft=
      cz2wpaux=(cz34(i3)%e(0)*cw56%e(0)-cz34(i3)%e(1)*cw56%e(1)-
     & cz34(i3)%e(2)*cw56%e(2)-cz34(i3)%e(3)*cw56%e(3))
  
* p.q -- p.q=cz2wmaux,p=cz34(i3)%e,q=cw78%e,bef=,aft=
      cz2wmaux=(cz34(i3)%e(0)*cw78%e(0)-cz34(i3)%e(1)*cw78%e(1)-
     & cz34(i3)%e(2)*cw78%e(2)-cz34(i3)%e(3)*cw78%e(3))
  
      cquart(i1,i3)=cquart(i1,i3)+rcotw*(2*czzaux*cwwaux -
     &   cz1wpaux*cz2wmaux - cz1wmaux*cz2wpaux)
  
      enddo
      enddo
  
* p.q -- p.q=cwwaux,p=cw56%e,q=cw78%e,bef=,aft=
      cwwaux=(cw56%e(0)*cw78%e(0)-cw56%e(1)*cw78%e(1)-cw56%e(2)*
     & cw78%e(2)-cw56%e(3)*cw78%e(3))
  
      do i1=1,2
      do i3=1,2
* p.q -- p.q=czzaux,p=cf12(i1)%e,q=cf34(i3)%e,bef=,aft=
      czzaux=(cf12(i1)%e(0)*cf34(i3)%e(0)-cf12(i1)%e(1)*cf34(i3)
     & %e(1)-cf12(i1)%e(2)*cf34(i3)%e(2)-cf12(i1)%e(3)*cf34(i3)%
     & e(3))
  
* p.q -- p.q=cz1wpaux,p=cf12(i1)%e,q=cw56%e,bef=,aft=
      cz1wpaux=(cf12(i1)%e(0)*cw56%e(0)-cf12(i1)%e(1)*cw56%e(1)-
     & cf12(i1)%e(2)*cw56%e(2)-cf12(i1)%e(3)*cw56%e(3))
  
* p.q -- p.q=cz1wmaux,p=cf12(i1)%e,q=cw78%e,bef=,aft=
      cz1wmaux=(cf12(i1)%e(0)*cw78%e(0)-cf12(i1)%e(1)*cw78%e(1)-
     & cf12(i1)%e(2)*cw78%e(2)-cf12(i1)%e(3)*cw78%e(3))
  
* p.q -- p.q=cz2wpaux,p=cf34(i3)%e,q=cw56%e,bef=,aft=
      cz2wpaux=(cf34(i3)%e(0)*cw56%e(0)-cf34(i3)%e(1)*cw56%e(1)-
     & cf34(i3)%e(2)*cw56%e(2)-cf34(i3)%e(3)*cw56%e(3))
  
* p.q -- p.q=cz2wmaux,p=cf34(i3)%e,q=cw78%e,bef=,aft=
      cz2wmaux=(cf34(i3)%e(0)*cw78%e(0)-cf34(i3)%e(1)*cw78%e(1)-
     & cf34(i3)%e(2)*cw78%e(2)-cf34(i3)%e(3)*cw78%e(3))
  
      cquart(i1,i3)=cquart(i1,i3)+(2*czzaux*cwwaux -
     &   cz1wpaux*cz2wmaux - cz1wmaux*cz2wpaux)
  
      enddo
      enddo
  
  
*                                                                       
* RESULT                                                                
*                                                                       
      spk0=sqrt(p1k0*p2k0*p3k0*p4k0*p5k0*p6k0*p7k0*p8k0)
*** electroweak conf 1                                                  
      do i1=1,2
      do i3=1,2
        cres(i1,i3,1)= (c3fk_tot(i1,i3)+cbct1234(i1,i3)+
     &        cbct1256(i1,i3)+cbct1278(i1,i3)+cquart(i1,i3))/spk0
  
      enddo
      enddo
  
**QCD                                                                   
* conf 2 12~34                                                          
      if (ilept(id1).ne.1.and.ilept(id3).ne.1) then
        do i1=1,2
        do i3=1,2
          cres(i1,i3,2)=(cg3fk1234(i1,i3)+cgbct1234(i1,i3))/spk0
        enddo
        enddo
      else
        do i1=1,2
        do i3=1,2
          cres(i1,i3,2)=czero
        enddo
        enddo
      endif
  
*conf 3                                                                 
      if (ilept(id1).ne.1.and.ilept(id5).ne.1) then
        do i1=1,2
        do i3=1,2
          cres(i1,i3,3)=(cg3fk56_123478(i1,i3)+cgbct1256(i1,i3)+
     &      cgct1256(i1,i3))/spk0
  
        enddo
        enddo
      else
        do i1=1,2
        do i3=1,2
          cres(i1,i3,3)=czero
  
        enddo
        enddo
      endif
*conf 4                                                                 
      if (ilept(id3).ne.1.and.ilept(id5).ne.1) then
        do i1=1,2
        do i3=1,2
          cres(i1,i3,4)=(cg3fk56_341278(i1,i3)+cgbct3456(i1,i3)+
     &      cgct3456(i1,i3))/spk0
  
        enddo
        enddo
      else
        do i1=1,2
        do i3=1,2
          cres(i1,i3,4)=czero
  
        enddo
        enddo
      endif
*conf 5                                                                 
      if (ilept(id1).ne.1.and.ilept(id7).ne.1) then
        do i1=1,2
        do i3=1,2
          cres(i1,i3,5)=(cg3fk78_123456(i1,i3)+cgbct1278(i1,i3)+
     &      cgct1278(i1,i3))/spk0
  
        enddo
        enddo
      else
        do i1=1,2
        do i3=1,2
          cres(i1,i3,5)=czero
  
        enddo
        enddo
      endif
*conf 6                                                                 
      if (ilept(id3).ne.1.and.ilept(id7).ne.1) then
        do i1=1,2
        do i3=1,2
          cres(i1,i3,6)=(cg3fk78_341256(i1,i3)+cgbct3478(i1,i3)+
     &      cgct3478(i1,i3))/spk0
  
        enddo
        enddo
      else
        do i1=1,2
        do i3=1,2
          cres(i1,i3,6)=czero
  
        enddo
        enddo
      endif
  
      return
      end
  
