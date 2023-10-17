************************************************************************
*                                                                       
* fourwqcd.f                                                            
*                                                                       
* Last update: Jul 26, 2006                                             
*                                                                       
* This routine computes the basic amplitude corresponding to 8          
* outgoing  massless fermions which can form 2 W+ and 2 W-.             
* All outgoing particles are considered different and no one has a      
* coupling with higgs (no b's or t's)                                   
* The input particles are ordered from 1 to 8 in such a way that        
* odd are particles, even antiparticles. 12 corresponds to W+,          
* 34 to W-, 56 to W+, 78 to W-                                          
* p1, ...p8 are all outgoing momenta                                    
* id1....id8 give the identities of the outgoing particles              
* cres  is on output the complex helicity amplitude                     
* For this amplitude all chirality indeces  = 2   (corresponding to -)  
*                                                                       
*                                                                       
* Added QCD LO contribution                                             
*                                                                       
* e.g.:                   |     |                                       
*                     >-W-|     |-W-<                                   
*                         |~~g~~|                                       
*                         |     |                                       
*                                                                       
* (QCD lines are labelled by **QCD)                                     
*                                                                       
* This subroutine returns the 7 contributions to the total amplitude    
* corresponding to 7 different colour flow configurations.              
* cres(1:7) contains the 7 colour configurations of a generic 4w        
* process. Several configurations may be zero, depending on the number  
* of leptons lines:                                                     
* 0 lepton line -> no configuration is identically 0                    
* 1 lepton line ->  3 configurations equal 0                            
* 2 lepton line ->  5 configurations equal 0                            
* 3 lepton line ->  6 configurations equal 0 (no QCD contributions)     
* 4 lepton line ->  6 configurations equal 0 (no QCD contributions)     
*                                                                       
* The colour flow configurations are the following (expressed in terms  
* of colour lines):                                                     
*                                                                       
* cres(1):   1| 3| 5| 7|   (pure electroweak contribution)              
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
* cres(7):   1| 5| 7| 3|   (QCD contribution)                           
*             |  |~~|  |                                                
*            2| 6| 8| 4|                                                
*                                                                       
*                                                                       
* As to the QCD contributions, it is clear that for each configuration  
* two different quark lines are connected by the gluon. This will be the
* criterion for separating the different colour configurations.         
* Following this convention, it is straightforward to state which       
* components of cres(1:7) are identically zero if one or more lepton    
* lines are present.                                                    
*                                                                       
* e.g.:    c s~ d u~ c' s'~ mu- vm~                                     
*          1 2  3 4  5  6    7   8    => cres(5)=cres(6)=cres(7)=0      
*                                                                       
* More generally, if one has n lepton lines and cares to put them in    
* the last positions, then the first 7-n components of cres(1:7) will be
* nonzero                                                               
*********************************************************************** 
  
      subroutine fourwqcd(p1,p2,p3,p4,p5,p6,p7,p8,
     &              id1,id2,id3,id4,id5,id6,id7,id8,cres)
  
  
      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)
  
  
      dimension cres(7)
  
  
*four momenta                                                           
  
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3),p7(0:3)
     &     ,p8(0:3)
  
      dimension p12(0:3),p34(0:3),p56(0:3),p78(0:3)
      dimension p134(0:3),p178(0:3),p312(0:3),p356(0:3),p534(0:3)
     &     ,p578(0:3),p712(0:3),p756(0:3)
      dimension p234(0:3),p278(0:3),p412(0:3),p456(0:3),p634(0:3)
     &     ,p678(0:3),p812(0:3),p856(0:3)
      dimension p1234(0:3),p1278(0:3)
  
* Z and gamma "decaying" to two w's                                     
      dimension cz1234(0:3),cz1278(0:3),cz5634(0:3),cz5678(0:3)
     &     ,cf1234(0:3),cf1278(0:3),cf5634(0:3),cf5678(0:3)
**QCD: g "decaying" to two w's                                          
     &     ,cg1234(0:3),cg1278(0:3),cg5634(0:3),cg5678(0:3)
     &     ,cg3412(0:3),cg7812(0:3),cg3456(0:3),cg7856(0:3)
  
*temporary auxiliary array                                              
      dimension cftemp(0:3)
  
* triple vertices                                                       
      dimension vfz(0:3), vwm(0:3), vwp(0:3)
     &     ,ctrip1234(0:3),ctrip1278(0:3),ctrip5634(0:3),ctrip5678(0:3)
  
* forks                                                                 
  
      type polcom
        double complex e(0:3),ek0,v
      end type
      type(polcom) cw12,cw34,cw56,cw78
  
* left t                                                                
  
      type l
         double complex a(2:2),c(2:2)
      end type
      type(l) l1_34,l1_78,l3_12,l3_56,l5_34,l5_78,l7_12,l7_56
     &     ,l1_3456,l1_7856,l3_1278,l3_5678,l5_3412,l5_7812,l7_1234
     &     ,l7_5634
* with z, gamma etc the second indices represent down four momentum     
*  and not incoming as for li_kl  et similar                            
     &     ,lz1_234(0:3),lz1_278(0:3),lz3_412(0:3),lz3_456(0:3),
     &     lz5_634(0:3),lz5_678(0:3),lz7_812(0:3),lz7_856(0:3),
     &     lf1_234(0:3),lf1_278(0:3),lf3_412(0:3),lf3_456(0:3),
     &     lf5_634(0:3),lf5_678(0:3),lf7_812(0:3),lf7_856(0:3)
  
* u t                                                                   
* only to be used with other w insertions, hence tw0 is enough          
      type u
         double complex a(2:2),b(1:1),c(2:2),d(1:1)
      end type
      type(u) u134_56,u178_56,u312_78,u356_78,u534_12,u578_12,u712_34
     &     ,u756_34
  
* right t                                                               
  
      type r
         double complex a(2:2),b(1:1)
      end type
      type(r) r2_34,r2_78,r4_12,r4_56,r6_34,r6_78,r8_12,r8_56
* with z, gamma etc the second indices represent down four momentum     
*  and not incoming as for ri_kl  et similar                            
     &     ,rz2_134(0:3),rz2_178(0:3),rz4_312(0:3),rz4_356(0:3),
     &     rz6_534(0:3),rz6_578(0:3),rz8_712(0:3),rz8_756(0:3),
     &     rf2_134(0:3),rf2_178(0:3),rf4_312(0:3),rf4_356(0:3),
     &     rf6_534(0:3),rf6_578(0:3),rf8_712(0:3),rf8_756(0:3)
  
  
*common blocks                                                          
  
* generic common                                                        
      include 'common.h'
  
  
* control structure                                                     
  
      if(ilept(id1).ne.ilept(id2)) then
        print*,'**FOURW ERROR: cannot consider a fermion line with a
     .quark and a lepton !'
        stop
      endif
      if(ilept(id3).ne.ilept(id4)) then
        print*,'**FOURW ERROR: cannot consider a fermion line with a
     .quark and a lepton !'
        stop
      endif
      if(ilept(id5).ne.ilept(id6)) then
        print*,'**FOURW ERROR: cannot consider a fermion line with a
     .quark and a lepton !'
        stop
      endif
      if(ilept(id7).ne.ilept(id8)) then
        print*,'**FOURW ERROR: cannot consider a fermion line with a
     .quark and a lepton !'
        stop
      endif
  
* 4momenta and their sums                                               
  
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
  
      do m=0,3
        p12(m)=p1(m)+p2(m)
      enddo
      do m=0,3
        p34(m)=p3(m)+p4(m)
      enddo
      do m=0,3
        p56(m)=p5(m)+p6(m)
      enddo
      do m=0,3
        p78(m)=p7(m)+p8(m)
      enddo
  
* compute all forks ( W's "decaying" to two fermions)                   
*                                                                       
  
* quqd -- p=p1,q=p2
      quqd=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      s12=2.d0*quqd
      cdw=-wcl/(s12-cmw2)
* TW10 -- qu=p1,qd=p2,v=0,a=cw12%e(0),cl=cdw,nsum=0
      eps_0=-p1(2)*p2(3)+p2(2)*p1(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p1k0*p2(0)+p2k0*p1(0)
      cw12%e(0)=cdw*(auxa-ceps_0)
* TW10 -- qu=p1,qd=p2,v=1,a=cw12%e(1),cl=cdw,nsum=0
      auxa=-quqd+p1k0*p2(1)+p2k0*p1(1)
      cw12%e(1)=cdw*(auxa-ceps_0)
* TW10 -- qu=p1,qd=p2,v=2,a=cw12%e(2),cl=cdw,nsum=0
      eps_0=-p1k0*p2(3)+p2k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(2)+p2k0*p1(2)
      cw12%e(2)=cdw*(auxa-ceps_0)
* TW10 -- qu=p1,qd=p2,v=3,a=cw12%e(3),cl=cdw,nsum=0
      eps_0=p1k0*p2(2)-p2k0*p1(2)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(3)+p2k0*p1(3)
      cw12%e(3)=cdw*(auxa-ceps_0)
* pk0 -- p=cw12%e
      cw12%ek0=cw12%e(0)-cw12%e(1)
* quqd -- p=p3,q=p4
      quqd=p3(0)*p4(0)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3)
      s34=2.d0*quqd
      cdw=-wcl/(s34-cmw2)
* TW10 -- qu=p3,qd=p4,v=0,a=cw34%e(0),cl=cdw,nsum=0
      eps_0=-p3(2)*p4(3)+p4(2)*p3(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p3k0*p4(0)+p4k0*p3(0)
      cw34%e(0)=cdw*(auxa-ceps_0)
* TW10 -- qu=p3,qd=p4,v=1,a=cw34%e(1),cl=cdw,nsum=0
      auxa=-quqd+p3k0*p4(1)+p4k0*p3(1)
      cw34%e(1)=cdw*(auxa-ceps_0)
* TW10 -- qu=p3,qd=p4,v=2,a=cw34%e(2),cl=cdw,nsum=0
      eps_0=-p3k0*p4(3)+p4k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(2)+p4k0*p3(2)
      cw34%e(2)=cdw*(auxa-ceps_0)
* TW10 -- qu=p3,qd=p4,v=3,a=cw34%e(3),cl=cdw,nsum=0
      eps_0=p3k0*p4(2)-p4k0*p3(2)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(3)+p4k0*p3(3)
      cw34%e(3)=cdw*(auxa-ceps_0)
* pk0 -- p=cw34%e
      cw34%ek0=cw34%e(0)-cw34%e(1)
* quqd -- p=p5,q=p6
      quqd=p5(0)*p6(0)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3)
      s56=2.d0*quqd
      cdw=-wcl/(s56-cmw2)
* TW10 -- qu=p5,qd=p6,v=0,a=cw56%e(0),cl=cdw,nsum=0
      eps_0=-p5(2)*p6(3)+p6(2)*p5(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p5k0*p6(0)+p6k0*p5(0)
      cw56%e(0)=cdw*(auxa-ceps_0)
* TW10 -- qu=p5,qd=p6,v=1,a=cw56%e(1),cl=cdw,nsum=0
      auxa=-quqd+p5k0*p6(1)+p6k0*p5(1)
      cw56%e(1)=cdw*(auxa-ceps_0)
* TW10 -- qu=p5,qd=p6,v=2,a=cw56%e(2),cl=cdw,nsum=0
      eps_0=-p5k0*p6(3)+p6k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(2)+p6k0*p5(2)
      cw56%e(2)=cdw*(auxa-ceps_0)
* TW10 -- qu=p5,qd=p6,v=3,a=cw56%e(3),cl=cdw,nsum=0
      eps_0=p5k0*p6(2)-p6k0*p5(2)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(3)+p6k0*p5(3)
      cw56%e(3)=cdw*(auxa-ceps_0)
* pk0 -- p=cw56%e
      cw56%ek0=cw56%e(0)-cw56%e(1)
* quqd -- p=p7,q=p8
      quqd=p7(0)*p8(0)-p7(1)*p8(1)-p7(2)*p8(2)-p7(3)*p8(3)
      s78=2.d0*quqd
      cdw=-wcl/(s78-cmw2)
* TW10 -- qu=p7,qd=p8,v=0,a=cw78%e(0),cl=cdw,nsum=0
      eps_0=-p7(2)*p8(3)+p8(2)*p7(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p7k0*p8(0)+p8k0*p7(0)
      cw78%e(0)=cdw*(auxa-ceps_0)
* TW10 -- qu=p7,qd=p8,v=1,a=cw78%e(1),cl=cdw,nsum=0
      auxa=-quqd+p7k0*p8(1)+p8k0*p7(1)
      cw78%e(1)=cdw*(auxa-ceps_0)
* TW10 -- qu=p7,qd=p8,v=2,a=cw78%e(2),cl=cdw,nsum=0
      eps_0=-p7k0*p8(3)+p8k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(2)+p8k0*p7(2)
      cw78%e(2)=cdw*(auxa-ceps_0)
* TW10 -- qu=p7,qd=p8,v=3,a=cw78%e(3),cl=cdw,nsum=0
      eps_0=p7k0*p8(2)-p8k0*p7(2)
      ceps_0=eps_0*cim
      auxa=p7k0*p8(3)+p8k0*p7(3)
      cw78%e(3)=cdw*(auxa-ceps_0)
* pk0 -- p=cw78%e
      cw78%ek0=cw78%e(0)-cw78%e(1)
  
* compute all single insertions of the type li_ (i=1,3,5,7)             
*                                                                       
*        i __                                                           
*            |_W__/                                                     
*            |    \                                                     
*                                                                       
*together with its propagator and pk0                                   
  
      do m=0,3
        p134(m)=p1(m)+p34(m)
      enddo
* pk0 -- p=p134
      p134k0=p134(0)-p134(1)
* p.q -- p.q=p134q,p=p134,q=p134,bef=,aft=
      p134q=(p134(0)*p134(0)-p134(1)*p134(1)-p134(2)*p134(2)-p13
     & 4(3)*p134(3))
* quqd -- p=p1,q=p134
      quqd=p1(0)*p134(0)-p1(1)*p134(1)-p1(2)*p134(2)-p1(3)*p134(
     & 3)
      ccl=wcl/(p134q*p134k0)
* TWL0 -- qu=p1,qd=p134,v=cw34%e,a=l1_34%a,c=l1_34%c,cl=ccl,nsum=0
      ceps_0=-cw34%ek0*(p1(2)*p134(3)-p134(2)*p1(3))+p1k0*(cw34%
     & e(2)*p134(3)-p134(2)*cw34%e(3))-p134k0*(cw34%e(2)*p1(3)-p
     & 1(2)*cw34%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw34%e(3)*p1k0+p1(3)*cw34%ek0
      ceps_1=ceps_1*cim
      cvqu=cw34%e(0)*p1(0)-cw34%e(1)*p1(1)-cw34%e(2)*p1(2)-cw34%
     & e(3)*p1(3)
      cvqd=cw34%e(0)*p134(0)-cw34%e(1)*p134(1)-cw34%e(2)*p134(2)
     & -cw34%e(3)*p134(3)
      cauxa=-cw34%ek0*quqd+p1k0*cvqd+p134k0*cvqu
      cauxc=+cw34%ek0*p1(2)-p1k0*cw34%e(2)
      l1_34%a(2)=ccl*(cauxa-ceps_0)
      l1_34%c(2)=ccl*(-cauxc+ceps_1)
      do m=0,3
        p178(m)=p1(m)+p78(m)
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
      do m=0,3
        p312(m)=p3(m)+p12(m)
      enddo
* pk0 -- p=p312
      p312k0=p312(0)-p312(1)
* p.q -- p.q=p312q,p=p312,q=p312,bef=,aft=
      p312q=(p312(0)*p312(0)-p312(1)*p312(1)-p312(2)*p312(2)-p31
     & 2(3)*p312(3))
* quqd -- p=p3,q=p312
      quqd=p3(0)*p312(0)-p3(1)*p312(1)-p3(2)*p312(2)-p3(3)*p312(
     & 3)
      ccl=wcl/(p312q*p312k0)
* TWL0 -- qu=p3,qd=p312,v=cw12%e,a=l3_12%a,c=l3_12%c,cl=ccl,nsum=0
      ceps_0=-cw12%ek0*(p3(2)*p312(3)-p312(2)*p3(3))+p3k0*(cw12%
     & e(2)*p312(3)-p312(2)*cw12%e(3))-p312k0*(cw12%e(2)*p3(3)-p
     & 3(2)*cw12%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw12%e(3)*p3k0+p3(3)*cw12%ek0
      ceps_1=ceps_1*cim
      cvqu=cw12%e(0)*p3(0)-cw12%e(1)*p3(1)-cw12%e(2)*p3(2)-cw12%
     & e(3)*p3(3)
      cvqd=cw12%e(0)*p312(0)-cw12%e(1)*p312(1)-cw12%e(2)*p312(2)
     & -cw12%e(3)*p312(3)
      cauxa=-cw12%ek0*quqd+p3k0*cvqd+p312k0*cvqu
      cauxc=+cw12%ek0*p3(2)-p3k0*cw12%e(2)
      l3_12%a(2)=ccl*(cauxa-ceps_0)
      l3_12%c(2)=ccl*(-cauxc+ceps_1)
      do m=0,3
        p356(m)=p3(m)+p56(m)
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
      do m=0,3
        p534(m)=p5(m)+p34(m)
      enddo
* pk0 -- p=p534
      p534k0=p534(0)-p534(1)
* p.q -- p.q=p534q,p=p534,q=p534,bef=,aft=
      p534q=(p534(0)*p534(0)-p534(1)*p534(1)-p534(2)*p534(2)-p53
     & 4(3)*p534(3))
* quqd -- p=p5,q=p534
      quqd=p5(0)*p534(0)-p5(1)*p534(1)-p5(2)*p534(2)-p5(3)*p534(
     & 3)
      ccl=wcl/(p534q*p534k0)
* TWL0 -- qu=p5,qd=p534,v=cw34%e,a=l5_34%a,c=l5_34%c,cl=ccl,nsum=0
      ceps_0=-cw34%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0*(cw34%
     & e(2)*p534(3)-p534(2)*cw34%e(3))-p534k0*(cw34%e(2)*p5(3)-p
     & 5(2)*cw34%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw34%e(3)*p5k0+p5(3)*cw34%ek0
      ceps_1=ceps_1*cim
      cvqu=cw34%e(0)*p5(0)-cw34%e(1)*p5(1)-cw34%e(2)*p5(2)-cw34%
     & e(3)*p5(3)
      cvqd=cw34%e(0)*p534(0)-cw34%e(1)*p534(1)-cw34%e(2)*p534(2)
     & -cw34%e(3)*p534(3)
      cauxa=-cw34%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cw34%ek0*p5(2)-p5k0*cw34%e(2)
      l5_34%a(2)=ccl*(cauxa-ceps_0)
      l5_34%c(2)=ccl*(-cauxc+ceps_1)
      do m=0,3
        p578(m)=p5(m)+p78(m)
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
        p712(m)=p7(m)+p12(m)
      enddo
* pk0 -- p=p712
      p712k0=p712(0)-p712(1)
* p.q -- p.q=p712q,p=p712,q=p712,bef=,aft=
      p712q=(p712(0)*p712(0)-p712(1)*p712(1)-p712(2)*p712(2)-p71
     & 2(3)*p712(3))
* quqd -- p=p7,q=p712
      quqd=p7(0)*p712(0)-p7(1)*p712(1)-p7(2)*p712(2)-p7(3)*p712(
     & 3)
      ccl=wcl/(p712q*p712k0)
* TWL0 -- qu=p7,qd=p712,v=cw12%e,a=l7_12%a,c=l7_12%c,cl=ccl,nsum=0
      ceps_0=-cw12%ek0*(p7(2)*p712(3)-p712(2)*p7(3))+p7k0*(cw12%
     & e(2)*p712(3)-p712(2)*cw12%e(3))-p712k0*(cw12%e(2)*p7(3)-p
     & 7(2)*cw12%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw12%e(3)*p7k0+p7(3)*cw12%ek0
      ceps_1=ceps_1*cim
      cvqu=cw12%e(0)*p7(0)-cw12%e(1)*p7(1)-cw12%e(2)*p7(2)-cw12%
     & e(3)*p7(3)
      cvqd=cw12%e(0)*p712(0)-cw12%e(1)*p712(1)-cw12%e(2)*p712(2)
     & -cw12%e(3)*p712(3)
      cauxa=-cw12%ek0*quqd+p7k0*cvqd+p712k0*cvqu
      cauxc=+cw12%ek0*p7(2)-p7k0*cw12%e(2)
      l7_12%a(2)=ccl*(cauxa-ceps_0)
      l7_12%c(2)=ccl*(-cauxc+ceps_1)
      do m=0,3
        p756(m)=p7(m)+p56(m)
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
  
* compute all single insertions of the type ri_ (i=2,4,6,8)             
*                                                                       
*            |                                                          
*            |_W__/                                                     
*        i __|    \                                                     
*                                                                       
*together with its propagator and pk0                                   
  
      do m=0,3
        p234(m)=-p2(m)-p34(m)
      enddo
* pk0 -- p=p234
      p234k0=p234(0)-p234(1)
* p.q -- p.q=p234q,p=p234,q=p234,bef=,aft=
      p234q=(p234(0)*p234(0)-p234(1)*p234(1)-p234(2)*p234(2)-p23
     & 4(3)*p234(3))
* quqd -- p=p234,q=p2
      quqd=p234(0)*p2(0)-p234(1)*p2(1)-p234(2)*p2(2)-p234(3)*p2(
     & 3)
      ccl=wcl/(p234q*p234k0)
* TWR0 -- qu=p234,qd=p2,v=cw34%e,a=r2_34%a,b=r2_34%b,cl=ccl,nsum=0
      ceps_0=-cw34%ek0*(p234(2)*p2(3)-p2(2)*p234(3))+p234k0*(cw3
     & 4%e(2)*p2(3)-p2(2)*cw34%e(3))-p2k0*(cw34%e(2)*p234(3)-p23
     & 4(2)*cw34%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw34%e(3)*p2k0+p2(3)*cw34%ek0
      ceps_2=ceps_2*cim
      cvqu=cw34%e(0)*p234(0)-cw34%e(1)*p234(1)-cw34%e(2)*p234(2)
     & -cw34%e(3)*p234(3)
      cvqd=cw34%e(0)*p2(0)-cw34%e(1)*p2(1)-cw34%e(2)*p2(2)-cw34%
     & e(3)*p2(3)
      cauxa=-cw34%ek0*quqd+p234k0*cvqd+p2k0*cvqu
      cauxb=-cw34%ek0*p2(2)+p2k0*cw34%e(2)
      r2_34%a(2)=ccl*(cauxa-ceps_0)
      r2_34%b(1)=ccl*(cauxb-ceps_2)
      do m=0,3
        p278(m)=-p2(m)-p78(m)
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
      do m=0,3
        p412(m)=-p4(m)-p12(m)
      enddo
* pk0 -- p=p412
      p412k0=p412(0)-p412(1)
* p.q -- p.q=p412q,p=p412,q=p412,bef=,aft=
      p412q=(p412(0)*p412(0)-p412(1)*p412(1)-p412(2)*p412(2)-p41
     & 2(3)*p412(3))
* quqd -- p=p412,q=p4
      quqd=p412(0)*p4(0)-p412(1)*p4(1)-p412(2)*p4(2)-p412(3)*p4(
     & 3)
      ccl=wcl/(p412q*p412k0)
* TWR0 -- qu=p412,qd=p4,v=cw12%e,a=r4_12%a,b=r4_12%b,cl=ccl,nsum=0
      ceps_0=-cw12%ek0*(p412(2)*p4(3)-p4(2)*p412(3))+p412k0*(cw1
     & 2%e(2)*p4(3)-p4(2)*cw12%e(3))-p4k0*(cw12%e(2)*p412(3)-p41
     & 2(2)*cw12%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw12%e(3)*p4k0+p4(3)*cw12%ek0
      ceps_2=ceps_2*cim
      cvqu=cw12%e(0)*p412(0)-cw12%e(1)*p412(1)-cw12%e(2)*p412(2)
     & -cw12%e(3)*p412(3)
      cvqd=cw12%e(0)*p4(0)-cw12%e(1)*p4(1)-cw12%e(2)*p4(2)-cw12%
     & e(3)*p4(3)
      cauxa=-cw12%ek0*quqd+p412k0*cvqd+p4k0*cvqu
      cauxb=-cw12%ek0*p4(2)+p4k0*cw12%e(2)
      r4_12%a(2)=ccl*(cauxa-ceps_0)
      r4_12%b(1)=ccl*(cauxb-ceps_2)
      do m=0,3
        p456(m)=-p4(m)-p56(m)
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
      do m=0,3
        p634(m)=-p6(m)-p34(m)
      enddo
* pk0 -- p=p634
      p634k0=p634(0)-p634(1)
* p.q -- p.q=p634q,p=p634,q=p634,bef=,aft=
      p634q=(p634(0)*p634(0)-p634(1)*p634(1)-p634(2)*p634(2)-p63
     & 4(3)*p634(3))
* quqd -- p=p634,q=p6
      quqd=p634(0)*p6(0)-p634(1)*p6(1)-p634(2)*p6(2)-p634(3)*p6(
     & 3)
      ccl=wcl/(p634q*p634k0)
* TWR0 -- qu=p634,qd=p6,v=cw34%e,a=r6_34%a,b=r6_34%b,cl=ccl,nsum=0
      ceps_0=-cw34%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634k0*(cw3
     & 4%e(2)*p6(3)-p6(2)*cw34%e(3))-p6k0*(cw34%e(2)*p634(3)-p63
     & 4(2)*cw34%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw34%e(3)*p6k0+p6(3)*cw34%ek0
      ceps_2=ceps_2*cim
      cvqu=cw34%e(0)*p634(0)-cw34%e(1)*p634(1)-cw34%e(2)*p634(2)
     & -cw34%e(3)*p634(3)
      cvqd=cw34%e(0)*p6(0)-cw34%e(1)*p6(1)-cw34%e(2)*p6(2)-cw34%
     & e(3)*p6(3)
      cauxa=-cw34%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cw34%ek0*p6(2)+p6k0*cw34%e(2)
      r6_34%a(2)=ccl*(cauxa-ceps_0)
      r6_34%b(1)=ccl*(cauxb-ceps_2)
      do m=0,3
        p678(m)=-p6(m)-p78(m)
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
        p812(m)=-p8(m)-p12(m)
      enddo
* pk0 -- p=p812
      p812k0=p812(0)-p812(1)
* p.q -- p.q=p812q,p=p812,q=p812,bef=,aft=
      p812q=(p812(0)*p812(0)-p812(1)*p812(1)-p812(2)*p812(2)-p81
     & 2(3)*p812(3))
* quqd -- p=p812,q=p8
      quqd=p812(0)*p8(0)-p812(1)*p8(1)-p812(2)*p8(2)-p812(3)*p8(
     & 3)
      ccl=wcl/(p812q*p812k0)
* TWR0 -- qu=p812,qd=p8,v=cw12%e,a=r8_12%a,b=r8_12%b,cl=ccl,nsum=0
      ceps_0=-cw12%ek0*(p812(2)*p8(3)-p8(2)*p812(3))+p812k0*(cw1
     & 2%e(2)*p8(3)-p8(2)*cw12%e(3))-p8k0*(cw12%e(2)*p812(3)-p81
     & 2(2)*cw12%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cw12%e(3)*p8k0+p8(3)*cw12%ek0
      ceps_2=ceps_2*cim
      cvqu=cw12%e(0)*p812(0)-cw12%e(1)*p812(1)-cw12%e(2)*p812(2)
     & -cw12%e(3)*p812(3)
      cvqd=cw12%e(0)*p8(0)-cw12%e(1)*p8(1)-cw12%e(2)*p8(2)-cw12%
     & e(3)*p8(3)
      cauxa=-cw12%ek0*quqd+p812k0*cvqd+p8k0*cvqu
      cauxb=-cw12%ek0*p8(2)+p8k0*cw12%e(2)
      r8_12%a(2)=ccl*(cauxa-ceps_0)
      r8_12%b(1)=ccl*(cauxb-ceps_2)
      do m=0,3
        p856(m)=-p8(m)-p56(m)
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
  
* compute all single insertions in the middle of a triple insertion     
*   line:  ex. u134_56                                                  
*                                                                       
*            |                                                          
*            |_W__/                                                     
*            |    \                                                     
*            |                                                          
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
* quqd -- p=p534,q=p678
      quqd=p534(0)*p678(0)-p534(1)*p678(1)-p534(2)*p678(2)-p534(
     & 3)*p678(3)
* TW0 -- qu=p534,qd=p678,v=cw12%e,a=u534_12%a,b=u534_12%b,c=u534_12%c,d=
* u534_12%d,cl=wcl,nsum=0
      ceps_0=-cw12%ek0*(p534(2)*p678(3)-p678(2)*p534(3))+p534k0*
     & (cw12%e(2)*p678(3)-p678(2)*cw12%e(3))-p678k0*(cw12%e(2)*p
     & 534(3)-p534(2)*cw12%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw12%e(3)*p534k0+p534(3)*cw12%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw12%e(3)*p678k0+p678(3)*cw12%ek0
      ceps_2=ceps_2*cim
      cvqu=cw12%e(0)*p534(0)-cw12%e(1)*p534(1)-cw12%e(2)*p534(2)
     & -cw12%e(3)*p534(3)
      cvqd=cw12%e(0)*p678(0)-cw12%e(1)*p678(1)-cw12%e(2)*p678(2)
     & -cw12%e(3)*p678(3)
      cauxa=-cw12%ek0*quqd+p534k0*cvqd+p678k0*cvqu
      cauxb=-cw12%ek0*p678(2)+p678k0*cw12%e(2)
      cauxc=+cw12%ek0*p534(2)-p534k0*cw12%e(2)
      u534_12%a(2)=wcl*(cauxa-ceps_0)
      u534_12%b(1)=wcl*(cauxb-ceps_2)
      u534_12%c(2)=wcl*(-cauxc+ceps_1)
      u534_12%d(1)=wcl*cw12%ek0
* quqd -- p=p578,q=p634
      quqd=p578(0)*p634(0)-p578(1)*p634(1)-p578(2)*p634(2)-p578(
     & 3)*p634(3)
* TW0 -- qu=p578,qd=p634,v=cw12%e,a=u578_12%a,b=u578_12%b,c=u578_12%c,d=
* u578_12%d,cl=wcl,nsum=0
      ceps_0=-cw12%ek0*(p578(2)*p634(3)-p634(2)*p578(3))+p578k0*
     & (cw12%e(2)*p634(3)-p634(2)*cw12%e(3))-p634k0*(cw12%e(2)*p
     & 578(3)-p578(2)*cw12%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw12%e(3)*p578k0+p578(3)*cw12%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw12%e(3)*p634k0+p634(3)*cw12%ek0
      ceps_2=ceps_2*cim
      cvqu=cw12%e(0)*p578(0)-cw12%e(1)*p578(1)-cw12%e(2)*p578(2)
     & -cw12%e(3)*p578(3)
      cvqd=cw12%e(0)*p634(0)-cw12%e(1)*p634(1)-cw12%e(2)*p634(2)
     & -cw12%e(3)*p634(3)
      cauxa=-cw12%ek0*quqd+p578k0*cvqd+p634k0*cvqu
      cauxb=-cw12%ek0*p634(2)+p634k0*cw12%e(2)
      cauxc=+cw12%ek0*p578(2)-p578k0*cw12%e(2)
      u578_12%a(2)=wcl*(cauxa-ceps_0)
      u578_12%b(1)=wcl*(cauxb-ceps_2)
      u578_12%c(2)=wcl*(-cauxc+ceps_1)
      u578_12%d(1)=wcl*cw12%ek0
* quqd -- p=p712,q=p856
      quqd=p712(0)*p856(0)-p712(1)*p856(1)-p712(2)*p856(2)-p712(
     & 3)*p856(3)
* TW0 -- qu=p712,qd=p856,v=cw34%e,a=u712_34%a,b=u712_34%b,c=u712_34%c,d=
* u712_34%d,cl=wcl,nsum=0
      ceps_0=-cw34%ek0*(p712(2)*p856(3)-p856(2)*p712(3))+p712k0*
     & (cw34%e(2)*p856(3)-p856(2)*cw34%e(3))-p856k0*(cw34%e(2)*p
     & 712(3)-p712(2)*cw34%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw34%e(3)*p712k0+p712(3)*cw34%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw34%e(3)*p856k0+p856(3)*cw34%ek0
      ceps_2=ceps_2*cim
      cvqu=cw34%e(0)*p712(0)-cw34%e(1)*p712(1)-cw34%e(2)*p712(2)
     & -cw34%e(3)*p712(3)
      cvqd=cw34%e(0)*p856(0)-cw34%e(1)*p856(1)-cw34%e(2)*p856(2)
     & -cw34%e(3)*p856(3)
      cauxa=-cw34%ek0*quqd+p712k0*cvqd+p856k0*cvqu
      cauxb=-cw34%ek0*p856(2)+p856k0*cw34%e(2)
      cauxc=+cw34%ek0*p712(2)-p712k0*cw34%e(2)
      u712_34%a(2)=wcl*(cauxa-ceps_0)
      u712_34%b(1)=wcl*(cauxb-ceps_2)
      u712_34%c(2)=wcl*(-cauxc+ceps_1)
      u712_34%d(1)=wcl*cw34%ek0
* quqd -- p=p756,q=p812
      quqd=p756(0)*p812(0)-p756(1)*p812(1)-p756(2)*p812(2)-p756(
     & 3)*p812(3)
* TW0 -- qu=p756,qd=p812,v=cw34%e,a=u756_34%a,b=u756_34%b,c=u756_34%c,d=
* u756_34%d,cl=wcl,nsum=0
      ceps_0=-cw34%ek0*(p756(2)*p812(3)-p812(2)*p756(3))+p756k0*
     & (cw34%e(2)*p812(3)-p812(2)*cw34%e(3))-p812k0*(cw34%e(2)*p
     & 756(3)-p756(2)*cw34%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw34%e(3)*p756k0+p756(3)*cw34%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw34%e(3)*p812k0+p812(3)*cw34%ek0
      ceps_2=ceps_2*cim
      cvqu=cw34%e(0)*p756(0)-cw34%e(1)*p756(1)-cw34%e(2)*p756(2)
     & -cw34%e(3)*p756(3)
      cvqd=cw34%e(0)*p812(0)-cw34%e(1)*p812(1)-cw34%e(2)*p812(2)
     & -cw34%e(3)*p812(3)
      cauxa=-cw34%ek0*quqd+p756k0*cvqd+p812k0*cvqu
      cauxb=-cw34%ek0*p812(2)+p812k0*cw34%e(2)
      cauxc=+cw34%ek0*p756(2)-p756k0*cw34%e(2)
      u756_34%a(2)=wcl*(cauxa-ceps_0)
      u756_34%b(1)=wcl*(cauxb-ceps_2)
      u756_34%c(2)=wcl*(-cauxc+ceps_1)
      u756_34%d(1)=wcl*cw34%ek0
  
* compute all double insertions  of the type li_jklm                    
*                                                                       
*        i __                                                           
*            |_W__/                                                     
*            |    \                                                     
*            |                                                          
*            |_W__/                                                     
*            |    \                                                     
*            |                                                          
*                                                                       
* ex. l1_3456                                                           
  
* TLT0_W -- aa=l1_3456%a,cc=l1_3456%c,a1=l1_34%a,c1=l1_34%c,a2=u134_56%a
* ,b2=u134_56%b,c2=u134_56%c,d2=u134_56%d,prq=p134q,nsum=0
      l1_3456%c(2)=l1_34%c(2)*p134q*u134_56%d(1)+l1_34%a(2)*u134
     & _56%c(2)
      l1_3456%a(2)=l1_34%c(2)*p134q*u134_56%b(1)+l1_34%a(2)*u134
     & _56%a(2)
* TLT0_W -- aa=l1_7856%a,cc=l1_7856%c,a1=l1_78%a,c1=l1_78%c,a2=u178_56%a
* ,b2=u178_56%b,c2=u178_56%c,d2=u178_56%d,prq=p178q,nsum=0
      l1_7856%c(2)=l1_78%c(2)*p178q*u178_56%d(1)+l1_78%a(2)*u178
     & _56%c(2)
      l1_7856%a(2)=l1_78%c(2)*p178q*u178_56%b(1)+l1_78%a(2)*u178
     & _56%a(2)
* TLT0_W -- aa=l3_1278%a,cc=l3_1278%c,a1=l3_12%a,c1=l3_12%c,a2=u312_78%a
* ,b2=u312_78%b,c2=u312_78%c,d2=u312_78%d,prq=p312q,nsum=0
      l3_1278%c(2)=l3_12%c(2)*p312q*u312_78%d(1)+l3_12%a(2)*u312
     & _78%c(2)
      l3_1278%a(2)=l3_12%c(2)*p312q*u312_78%b(1)+l3_12%a(2)*u312
     & _78%a(2)
* TLT0_W -- aa=l3_5678%a,cc=l3_5678%c,a1=l3_56%a,c1=l3_56%c,a2=u356_78%a
* ,b2=u356_78%b,c2=u356_78%c,d2=u356_78%d,prq=p356q,nsum=0
      l3_5678%c(2)=l3_56%c(2)*p356q*u356_78%d(1)+l3_56%a(2)*u356
     & _78%c(2)
      l3_5678%a(2)=l3_56%c(2)*p356q*u356_78%b(1)+l3_56%a(2)*u356
     & _78%a(2)
* TLT0_W -- aa=l5_3412%a,cc=l5_3412%c,a1=l5_34%a,c1=l5_34%c,a2=u534_12%a
* ,b2=u534_12%b,c2=u534_12%c,d2=u534_12%d,prq=p534q,nsum=0
      l5_3412%c(2)=l5_34%c(2)*p534q*u534_12%d(1)+l5_34%a(2)*u534
     & _12%c(2)
      l5_3412%a(2)=l5_34%c(2)*p534q*u534_12%b(1)+l5_34%a(2)*u534
     & _12%a(2)
* TLT0_W -- aa=l5_7812%a,cc=l5_7812%c,a1=l5_78%a,c1=l5_78%c,a2=u578_12%a
* ,b2=u578_12%b,c2=u578_12%c,d2=u578_12%d,prq=p578q,nsum=0
      l5_7812%c(2)=l5_78%c(2)*p578q*u578_12%d(1)+l5_78%a(2)*u578
     & _12%c(2)
      l5_7812%a(2)=l5_78%c(2)*p578q*u578_12%b(1)+l5_78%a(2)*u578
     & _12%a(2)
* TLT0_W -- aa=l7_1234%a,cc=l7_1234%c,a1=l7_12%a,c1=l7_12%c,a2=u712_34%a
* ,b2=u712_34%b,c2=u712_34%c,d2=u712_34%d,prq=p712q,nsum=0
      l7_1234%c(2)=l7_12%c(2)*p712q*u712_34%d(1)+l7_12%a(2)*u712
     & _34%c(2)
      l7_1234%a(2)=l7_12%c(2)*p712q*u712_34%b(1)+l7_12%a(2)*u712
     & _34%a(2)
* TLT0_W -- aa=l7_5634%a,cc=l7_5634%c,a1=l7_56%a,c1=l7_56%c,a2=u756_34%a
* ,b2=u756_34%b,c2=u756_34%c,d2=u756_34%d,prq=p756q,nsum=0
      l7_5634%c(2)=l7_56%c(2)*p756q*u756_34%d(1)+l7_56%a(2)*u756
     & _34%c(2)
      l7_5634%a(2)=l7_56%c(2)*p756q*u756_34%b(1)+l7_56%a(2)*u756
     & _34%a(2)
  
* compute all single insertions of the type lzi_ and lfi_ (i=1,3,5,7)   
*        i __              i __                                         
*            |_Z__(mu)         |__gamma__(mu)                           
*            |                 |                                        
*                                                                       
* As they are then connected with w insertions, they are calculated     
* as w insertions with the appropriate coupling cl only.                
  
  
* quqd -- p=p1,q=p234
      quqd=p1(0)*p234(0)-p1(1)*p234(1)-p1(2)*p234(2)-p1(3)*p234(
     & 3)
* TWL0 -- qu=p1,qd=p234,v=0,a=lz1_234(0)%a,c=lz1_234(0)%c,cl=zcl(id1),ns
* um=0
      eps_0=-p1(2)*p234(3)+p234(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p234(0)+p234k0*p1(0)
      lz1_234(0)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(0)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p234,v=1,a=lz1_234(1)%a,c=lz1_234(1)%c,cl=zcl(id1),ns
* um=0
      auxa=-quqd+p1k0*p234(1)+p234k0*p1(1)
      lz1_234(1)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(1)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p234,v=2,a=lz1_234(2)%a,c=lz1_234(2)%c,cl=zcl(id1),ns
* um=0
      eps_0=-p1k0*p234(3)+p234k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p234(2)+p234k0*p1(2)
      lz1_234(2)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(2)%c(2)=-zcl(id1)*p1k0
* TWL0 -- qu=p1,qd=p234,v=3,a=lz1_234(3)%a,c=lz1_234(3)%c,cl=zcl(id1),ns
* um=0
      eps_0=p1k0*p234(2)-p234k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p234(3)+p234k0*p1(3)
      lz1_234(3)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(3)%c(2)=zcl(id1)*ceps_1
  
  
  
      if(ineutri(id1).ne.1) then
* TWL0 -- qu=p1,qd=p234,v=0,a=lf1_234(0)%a,c=lf1_234(0)%c,cl=fcl(id1),ns
* um=0
      eps_0=-p1(2)*p234(3)+p234(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p234(0)+p234k0*p1(0)
      lf1_234(0)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(0)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p234,v=1,a=lf1_234(1)%a,c=lf1_234(1)%c,cl=fcl(id1),ns
* um=0
      auxa=-quqd+p1k0*p234(1)+p234k0*p1(1)
      lf1_234(1)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(1)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p234,v=2,a=lf1_234(2)%a,c=lf1_234(2)%c,cl=fcl(id1),ns
* um=0
      eps_0=-p1k0*p234(3)+p234k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p234(2)+p234k0*p1(2)
      lf1_234(2)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(2)%c(2)=-fcl(id1)*p1k0
* TWL0 -- qu=p1,qd=p234,v=3,a=lf1_234(3)%a,c=lf1_234(3)%c,cl=fcl(id1),ns
* um=0
      eps_0=p1k0*p234(2)-p234k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p234(3)+p234k0*p1(3)
      lf1_234(3)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(3)%c(2)=fcl(id1)*ceps_1
  
      else
        do mu=0,3
          lf1_234(mu)%a(2) = czero
          lf1_234(mu)%c(2) = czero
        enddo
      endif
  
  
* quqd -- p=p1,q=p278
      quqd=p1(0)*p278(0)-p1(1)*p278(1)-p1(2)*p278(2)-p1(3)*p278(
     & 3)
* TWL0 -- qu=p1,qd=p278,v=0,a=lz1_278(0)%a,c=lz1_278(0)%c,cl=zcl(id1),ns
* um=0
      eps_0=-p1(2)*p278(3)+p278(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p278(0)+p278k0*p1(0)
      lz1_278(0)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_278(0)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p278,v=1,a=lz1_278(1)%a,c=lz1_278(1)%c,cl=zcl(id1),ns
* um=0
      auxa=-quqd+p1k0*p278(1)+p278k0*p1(1)
      lz1_278(1)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_278(1)%c(2)=zcl(id1)*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p278,v=2,a=lz1_278(2)%a,c=lz1_278(2)%c,cl=zcl(id1),ns
* um=0
      eps_0=-p1k0*p278(3)+p278k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p278(2)+p278k0*p1(2)
      lz1_278(2)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_278(2)%c(2)=-zcl(id1)*p1k0
* TWL0 -- qu=p1,qd=p278,v=3,a=lz1_278(3)%a,c=lz1_278(3)%c,cl=zcl(id1),ns
* um=0
      eps_0=p1k0*p278(2)-p278k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p278(3)+p278k0*p1(3)
      lz1_278(3)%a(2)=zcl(id1)*(auxa-ceps_0)
      lz1_278(3)%c(2)=zcl(id1)*ceps_1
  
  
  
      if(ineutri(id1).ne.1) then
* TWL0 -- qu=p1,qd=p278,v=0,a=lf1_278(0)%a,c=lf1_278(0)%c,cl=fcl(id1),ns
* um=0
      eps_0=-p1(2)*p278(3)+p278(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      auxa=-quqd+p1k0*p278(0)+p278k0*p1(0)
      lf1_278(0)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_278(0)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p278,v=1,a=lf1_278(1)%a,c=lf1_278(1)%c,cl=fcl(id1),ns
* um=0
      auxa=-quqd+p1k0*p278(1)+p278k0*p1(1)
      lf1_278(1)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_278(1)%c(2)=fcl(id1)*(-p1(2)+ceps_1)
* TWL0 -- qu=p1,qd=p278,v=2,a=lf1_278(2)%a,c=lf1_278(2)%c,cl=fcl(id1),ns
* um=0
      eps_0=-p1k0*p278(3)+p278k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p278(2)+p278k0*p1(2)
      lf1_278(2)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_278(2)%c(2)=-fcl(id1)*p1k0
* TWL0 -- qu=p1,qd=p278,v=3,a=lf1_278(3)%a,c=lf1_278(3)%c,cl=fcl(id1),ns
* um=0
      eps_0=p1k0*p278(2)-p278k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      auxa=p1k0*p278(3)+p278k0*p1(3)
      lf1_278(3)%a(2)=fcl(id1)*(auxa-ceps_0)
      lf1_278(3)%c(2)=fcl(id1)*ceps_1
  
      else
        do mu=0,3
          lf1_278(mu)%a(2) = czero
          lf1_278(mu)%c(2) = czero
        enddo
      endif
  
  
* quqd -- p=p3,q=p412
      quqd=p3(0)*p412(0)-p3(1)*p412(1)-p3(2)*p412(2)-p3(3)*p412(
     & 3)
* TWL0 -- qu=p3,qd=p412,v=0,a=lz3_412(0)%a,c=lz3_412(0)%c,cl=zcl(id3),ns
* um=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lz3_412(0)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(0)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p412,v=1,a=lz3_412(1)%a,c=lz3_412(1)%c,cl=zcl(id3),ns
* um=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lz3_412(1)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(1)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p412,v=2,a=lz3_412(2)%a,c=lz3_412(2)%c,cl=zcl(id3),ns
* um=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lz3_412(2)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(2)%c(2)=-zcl(id3)*p3k0
* TWL0 -- qu=p3,qd=p412,v=3,a=lz3_412(3)%a,c=lz3_412(3)%c,cl=zcl(id3),ns
* um=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p412(3)+p412k0*p3(3)
      lz3_412(3)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(3)%c(2)=zcl(id3)*ceps_1
  
  
  
      if(ineutri(id3).ne.1) then
* TWL0 -- qu=p3,qd=p412,v=0,a=lf3_412(0)%a,c=lf3_412(0)%c,cl=fcl(id3),ns
* um=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lf3_412(0)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(0)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p412,v=1,a=lf3_412(1)%a,c=lf3_412(1)%c,cl=fcl(id3),ns
* um=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lf3_412(1)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(1)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p412,v=2,a=lf3_412(2)%a,c=lf3_412(2)%c,cl=fcl(id3),ns
* um=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lf3_412(2)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(2)%c(2)=-fcl(id3)*p3k0
* TWL0 -- qu=p3,qd=p412,v=3,a=lf3_412(3)%a,c=lf3_412(3)%c,cl=fcl(id3),ns
* um=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p412(3)+p412k0*p3(3)
      lf3_412(3)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(3)%c(2)=fcl(id3)*ceps_1
  
      else
        do mu=0,3
          lf3_412(mu)%a(2) = czero
          lf3_412(mu)%c(2) = czero
        enddo
      endif
  
  
* quqd -- p=p3,q=p456
      quqd=p3(0)*p456(0)-p3(1)*p456(1)-p3(2)*p456(2)-p3(3)*p456(
     & 3)
* TWL0 -- qu=p3,qd=p456,v=0,a=lz3_456(0)%a,c=lz3_456(0)%c,cl=zcl(id3),ns
* um=0
      eps_0=-p3(2)*p456(3)+p456(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p456(0)+p456k0*p3(0)
      lz3_456(0)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_456(0)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p456,v=1,a=lz3_456(1)%a,c=lz3_456(1)%c,cl=zcl(id3),ns
* um=0
      auxa=-quqd+p3k0*p456(1)+p456k0*p3(1)
      lz3_456(1)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_456(1)%c(2)=zcl(id3)*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p456,v=2,a=lz3_456(2)%a,c=lz3_456(2)%c,cl=zcl(id3),ns
* um=0
      eps_0=-p3k0*p456(3)+p456k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p456(2)+p456k0*p3(2)
      lz3_456(2)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_456(2)%c(2)=-zcl(id3)*p3k0
* TWL0 -- qu=p3,qd=p456,v=3,a=lz3_456(3)%a,c=lz3_456(3)%c,cl=zcl(id3),ns
* um=0
      eps_0=p3k0*p456(2)-p456k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p456(3)+p456k0*p3(3)
      lz3_456(3)%a(2)=zcl(id3)*(auxa-ceps_0)
      lz3_456(3)%c(2)=zcl(id3)*ceps_1
  
  
  
      if(ineutri(id3).ne.1) then
* TWL0 -- qu=p3,qd=p456,v=0,a=lf3_456(0)%a,c=lf3_456(0)%c,cl=fcl(id3),ns
* um=0
      eps_0=-p3(2)*p456(3)+p456(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      auxa=-quqd+p3k0*p456(0)+p456k0*p3(0)
      lf3_456(0)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_456(0)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p456,v=1,a=lf3_456(1)%a,c=lf3_456(1)%c,cl=fcl(id3),ns
* um=0
      auxa=-quqd+p3k0*p456(1)+p456k0*p3(1)
      lf3_456(1)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_456(1)%c(2)=fcl(id3)*(-p3(2)+ceps_1)
* TWL0 -- qu=p3,qd=p456,v=2,a=lf3_456(2)%a,c=lf3_456(2)%c,cl=fcl(id3),ns
* um=0
      eps_0=-p3k0*p456(3)+p456k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p456(2)+p456k0*p3(2)
      lf3_456(2)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_456(2)%c(2)=-fcl(id3)*p3k0
* TWL0 -- qu=p3,qd=p456,v=3,a=lf3_456(3)%a,c=lf3_456(3)%c,cl=fcl(id3),ns
* um=0
      eps_0=p3k0*p456(2)-p456k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      auxa=p3k0*p456(3)+p456k0*p3(3)
      lf3_456(3)%a(2)=fcl(id3)*(auxa-ceps_0)
      lf3_456(3)%c(2)=fcl(id3)*ceps_1
  
      else
        do mu=0,3
          lf3_456(mu)%a(2) = czero
          lf3_456(mu)%c(2) = czero
        enddo
      endif
  
  
* quqd -- p=p5,q=p634
      quqd=p5(0)*p634(0)-p5(1)*p634(1)-p5(2)*p634(2)-p5(3)*p634(
     & 3)
* TWL0 -- qu=p5,qd=p634,v=0,a=lz5_634(0)%a,c=lz5_634(0)%c,cl=zcl(id5),ns
* um=0
      eps_0=-p5(2)*p634(3)+p634(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p634(0)+p634k0*p5(0)
      lz5_634(0)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_634(0)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p634,v=1,a=lz5_634(1)%a,c=lz5_634(1)%c,cl=zcl(id5),ns
* um=0
      auxa=-quqd+p5k0*p634(1)+p634k0*p5(1)
      lz5_634(1)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_634(1)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p634,v=2,a=lz5_634(2)%a,c=lz5_634(2)%c,cl=zcl(id5),ns
* um=0
      eps_0=-p5k0*p634(3)+p634k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p634(2)+p634k0*p5(2)
      lz5_634(2)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_634(2)%c(2)=-zcl(id5)*p5k0
* TWL0 -- qu=p5,qd=p634,v=3,a=lz5_634(3)%a,c=lz5_634(3)%c,cl=zcl(id5),ns
* um=0
      eps_0=p5k0*p634(2)-p634k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p634(3)+p634k0*p5(3)
      lz5_634(3)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_634(3)%c(2)=zcl(id5)*ceps_1
  
  
  
      if(ineutri(id5).ne.1) then
* TWL0 -- qu=p5,qd=p634,v=0,a=lf5_634(0)%a,c=lf5_634(0)%c,cl=fcl(id5),ns
* um=0
      eps_0=-p5(2)*p634(3)+p634(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p634(0)+p634k0*p5(0)
      lf5_634(0)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_634(0)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p634,v=1,a=lf5_634(1)%a,c=lf5_634(1)%c,cl=fcl(id5),ns
* um=0
      auxa=-quqd+p5k0*p634(1)+p634k0*p5(1)
      lf5_634(1)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_634(1)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p634,v=2,a=lf5_634(2)%a,c=lf5_634(2)%c,cl=fcl(id5),ns
* um=0
      eps_0=-p5k0*p634(3)+p634k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p634(2)+p634k0*p5(2)
      lf5_634(2)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_634(2)%c(2)=-fcl(id5)*p5k0
* TWL0 -- qu=p5,qd=p634,v=3,a=lf5_634(3)%a,c=lf5_634(3)%c,cl=fcl(id5),ns
* um=0
      eps_0=p5k0*p634(2)-p634k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p634(3)+p634k0*p5(3)
      lf5_634(3)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_634(3)%c(2)=fcl(id5)*ceps_1
  
      else
        do mu=0,3
          lf5_634(mu)%a(2) = czero
          lf5_634(mu)%c(2) = czero
        enddo
      endif
  
  
* quqd -- p=p5,q=p678
      quqd=p5(0)*p678(0)-p5(1)*p678(1)-p5(2)*p678(2)-p5(3)*p678(
     & 3)
* TWL0 -- qu=p5,qd=p678,v=0,a=lz5_678(0)%a,c=lz5_678(0)%c,cl=zcl(id5),ns
* um=0
      eps_0=-p5(2)*p678(3)+p678(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p678(0)+p678k0*p5(0)
      lz5_678(0)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_678(0)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p678,v=1,a=lz5_678(1)%a,c=lz5_678(1)%c,cl=zcl(id5),ns
* um=0
      auxa=-quqd+p5k0*p678(1)+p678k0*p5(1)
      lz5_678(1)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_678(1)%c(2)=zcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p678,v=2,a=lz5_678(2)%a,c=lz5_678(2)%c,cl=zcl(id5),ns
* um=0
      eps_0=-p5k0*p678(3)+p678k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p678(2)+p678k0*p5(2)
      lz5_678(2)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_678(2)%c(2)=-zcl(id5)*p5k0
* TWL0 -- qu=p5,qd=p678,v=3,a=lz5_678(3)%a,c=lz5_678(3)%c,cl=zcl(id5),ns
* um=0
      eps_0=p5k0*p678(2)-p678k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p678(3)+p678k0*p5(3)
      lz5_678(3)%a(2)=zcl(id5)*(auxa-ceps_0)
      lz5_678(3)%c(2)=zcl(id5)*ceps_1
  
  
  
      if(ineutri(id5).ne.1) then
* TWL0 -- qu=p5,qd=p678,v=0,a=lf5_678(0)%a,c=lf5_678(0)%c,cl=fcl(id5),ns
* um=0
      eps_0=-p5(2)*p678(3)+p678(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p678(0)+p678k0*p5(0)
      lf5_678(0)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_678(0)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p678,v=1,a=lf5_678(1)%a,c=lf5_678(1)%c,cl=fcl(id5),ns
* um=0
      auxa=-quqd+p5k0*p678(1)+p678k0*p5(1)
      lf5_678(1)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_678(1)%c(2)=fcl(id5)*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p678,v=2,a=lf5_678(2)%a,c=lf5_678(2)%c,cl=fcl(id5),ns
* um=0
      eps_0=-p5k0*p678(3)+p678k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p678(2)+p678k0*p5(2)
      lf5_678(2)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_678(2)%c(2)=-fcl(id5)*p5k0
* TWL0 -- qu=p5,qd=p678,v=3,a=lf5_678(3)%a,c=lf5_678(3)%c,cl=fcl(id5),ns
* um=0
      eps_0=p5k0*p678(2)-p678k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p678(3)+p678k0*p5(3)
      lf5_678(3)%a(2)=fcl(id5)*(auxa-ceps_0)
      lf5_678(3)%c(2)=fcl(id5)*ceps_1
  
      else
        do mu=0,3
          lf5_678(mu)%a(2) = czero
          lf5_678(mu)%c(2) = czero
        enddo
      endif
  
  
* quqd -- p=p7,q=p812
      quqd=p7(0)*p812(0)-p7(1)*p812(1)-p7(2)*p812(2)-p7(3)*p812(
     & 3)
* TWL0 -- qu=p7,qd=p812,v=0,a=lz7_812(0)%a,c=lz7_812(0)%c,cl=zcl(id7),ns
* um=0
      eps_0=-p7(2)*p812(3)+p812(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p812(0)+p812k0*p7(0)
      lz7_812(0)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_812(0)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=1,a=lz7_812(1)%a,c=lz7_812(1)%c,cl=zcl(id7),ns
* um=0
      auxa=-quqd+p7k0*p812(1)+p812k0*p7(1)
      lz7_812(1)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_812(1)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=2,a=lz7_812(2)%a,c=lz7_812(2)%c,cl=zcl(id7),ns
* um=0
      eps_0=-p7k0*p812(3)+p812k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p812(2)+p812k0*p7(2)
      lz7_812(2)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_812(2)%c(2)=-zcl(id7)*p7k0
* TWL0 -- qu=p7,qd=p812,v=3,a=lz7_812(3)%a,c=lz7_812(3)%c,cl=zcl(id7),ns
* um=0
      eps_0=p7k0*p812(2)-p812k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p812(3)+p812k0*p7(3)
      lz7_812(3)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_812(3)%c(2)=zcl(id7)*ceps_1
  
  
  
      if(ineutri(id7).ne.1) then
* TWL0 -- qu=p7,qd=p812,v=0,a=lf7_812(0)%a,c=lf7_812(0)%c,cl=fcl(id7),ns
* um=0
      eps_0=-p7(2)*p812(3)+p812(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p812(0)+p812k0*p7(0)
      lf7_812(0)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_812(0)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=1,a=lf7_812(1)%a,c=lf7_812(1)%c,cl=fcl(id7),ns
* um=0
      auxa=-quqd+p7k0*p812(1)+p812k0*p7(1)
      lf7_812(1)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_812(1)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=2,a=lf7_812(2)%a,c=lf7_812(2)%c,cl=fcl(id7),ns
* um=0
      eps_0=-p7k0*p812(3)+p812k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p812(2)+p812k0*p7(2)
      lf7_812(2)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_812(2)%c(2)=-fcl(id7)*p7k0
* TWL0 -- qu=p7,qd=p812,v=3,a=lf7_812(3)%a,c=lf7_812(3)%c,cl=fcl(id7),ns
* um=0
      eps_0=p7k0*p812(2)-p812k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p812(3)+p812k0*p7(3)
      lf7_812(3)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_812(3)%c(2)=fcl(id7)*ceps_1
  
      else
        do mu=0,3
          lf7_812(mu)%a(2) = czero
          lf7_812(mu)%c(2) = czero
        enddo
      endif
  
  
* quqd -- p=p7,q=p856
      quqd=p7(0)*p856(0)-p7(1)*p856(1)-p7(2)*p856(2)-p7(3)*p856(
     & 3)
* TWL0 -- qu=p7,qd=p856,v=0,a=lz7_856(0)%a,c=lz7_856(0)%c,cl=zcl(id7),ns
* um=0
      eps_0=-p7(2)*p856(3)+p856(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p856(0)+p856k0*p7(0)
      lz7_856(0)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_856(0)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p856,v=1,a=lz7_856(1)%a,c=lz7_856(1)%c,cl=zcl(id7),ns
* um=0
      auxa=-quqd+p7k0*p856(1)+p856k0*p7(1)
      lz7_856(1)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_856(1)%c(2)=zcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p856,v=2,a=lz7_856(2)%a,c=lz7_856(2)%c,cl=zcl(id7),ns
* um=0
      eps_0=-p7k0*p856(3)+p856k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p856(2)+p856k0*p7(2)
      lz7_856(2)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_856(2)%c(2)=-zcl(id7)*p7k0
* TWL0 -- qu=p7,qd=p856,v=3,a=lz7_856(3)%a,c=lz7_856(3)%c,cl=zcl(id7),ns
* um=0
      eps_0=p7k0*p856(2)-p856k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p856(3)+p856k0*p7(3)
      lz7_856(3)%a(2)=zcl(id7)*(auxa-ceps_0)
      lz7_856(3)%c(2)=zcl(id7)*ceps_1
  
  
  
      if(ineutri(id7).ne.1) then
* TWL0 -- qu=p7,qd=p856,v=0,a=lf7_856(0)%a,c=lf7_856(0)%c,cl=fcl(id7),ns
* um=0
      eps_0=-p7(2)*p856(3)+p856(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p856(0)+p856k0*p7(0)
      lf7_856(0)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_856(0)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p856,v=1,a=lf7_856(1)%a,c=lf7_856(1)%c,cl=fcl(id7),ns
* um=0
      auxa=-quqd+p7k0*p856(1)+p856k0*p7(1)
      lf7_856(1)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_856(1)%c(2)=fcl(id7)*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p856,v=2,a=lf7_856(2)%a,c=lf7_856(2)%c,cl=fcl(id7),ns
* um=0
      eps_0=-p7k0*p856(3)+p856k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p856(2)+p856k0*p7(2)
      lf7_856(2)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_856(2)%c(2)=-fcl(id7)*p7k0
* TWL0 -- qu=p7,qd=p856,v=3,a=lf7_856(3)%a,c=lf7_856(3)%c,cl=fcl(id7),ns
* um=0
      eps_0=p7k0*p856(2)-p856k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p856(3)+p856k0*p7(3)
      lf7_856(3)%a(2)=fcl(id7)*(auxa-ceps_0)
      lf7_856(3)%c(2)=fcl(id7)*ceps_1
  
      else
        do mu=0,3
          lf7_856(mu)%a(2) = czero
          lf7_856(mu)%c(2) = czero
        enddo
      endif
  
  
  
* compute all single insertions of the type rzi_ and rfi_ (i=2,4,6,8)   
*                                                                       
*            |_Z__(mu)         |__gamma__(mu)                           
*        i __|             i __|                                        
*                                                                       
* As they are then connected with w insertions, they are calculated     
* as w insertions with the appropriate coupling cl only                 
  
  
* quqd -- p=p134,q=p2
      quqd=p134(0)*p2(0)-p134(1)*p2(1)-p134(2)*p2(2)-p134(3)*p2(
     & 3)
* TWR0 -- qu=p134,qd=p2,v=0,a=rz2_134(0)%a,b=rz2_134(0)%b,cl=zcl(id2),ns
* um=0
      eps_0=-p134(2)*p2(3)+p2(2)*p134(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p134k0*p2(0)+p2k0*p134(0)
      rz2_134(0)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(0)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
* TWR0 -- qu=p134,qd=p2,v=1,a=rz2_134(1)%a,b=rz2_134(1)%b,cl=zcl(id2),ns
* um=0
      auxa=-quqd+p134k0*p2(1)+p2k0*p134(1)
      rz2_134(1)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(1)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
* TWR0 -- qu=p134,qd=p2,v=2,a=rz2_134(2)%a,b=rz2_134(2)%b,cl=zcl(id2),ns
* um=0
      eps_0=-p134k0*p2(3)+p2k0*p134(3)
      ceps_0=eps_0*cim
      auxa=p134k0*p2(2)+p2k0*p134(2)
      rz2_134(2)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(2)%b(1)=-zcl(id2)*p2k0
* TWR0 -- qu=p134,qd=p2,v=3,a=rz2_134(3)%a,b=rz2_134(3)%b,cl=zcl(id2),ns
* um=0
      eps_0=p134k0*p2(2)-p2k0*p134(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p134k0*p2(3)+p2k0*p134(3)
      rz2_134(3)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(3)%b(1)=-zcl(id2)*ceps_2
  
  
  
      if(ineutri(id2).ne.1) then
* TWR0 -- qu=p134,qd=p2,v=0,a=rf2_134(0)%a,b=rf2_134(0)%b,cl=fcl(id2),ns
* um=0
      eps_0=-p134(2)*p2(3)+p2(2)*p134(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p134k0*p2(0)+p2k0*p134(0)
      rf2_134(0)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(0)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
* TWR0 -- qu=p134,qd=p2,v=1,a=rf2_134(1)%a,b=rf2_134(1)%b,cl=fcl(id2),ns
* um=0
      auxa=-quqd+p134k0*p2(1)+p2k0*p134(1)
      rf2_134(1)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(1)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
* TWR0 -- qu=p134,qd=p2,v=2,a=rf2_134(2)%a,b=rf2_134(2)%b,cl=fcl(id2),ns
* um=0
      eps_0=-p134k0*p2(3)+p2k0*p134(3)
      ceps_0=eps_0*cim
      auxa=p134k0*p2(2)+p2k0*p134(2)
      rf2_134(2)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(2)%b(1)=-fcl(id2)*p2k0
* TWR0 -- qu=p134,qd=p2,v=3,a=rf2_134(3)%a,b=rf2_134(3)%b,cl=fcl(id2),ns
* um=0
      eps_0=p134k0*p2(2)-p2k0*p134(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p134k0*p2(3)+p2k0*p134(3)
      rf2_134(3)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(3)%b(1)=-fcl(id2)*ceps_2
  
      else
        do mu=0,3
          rf2_134(mu)%a(2) = czero
          rf2_134(mu)%b(1) = czero
        enddo
      endif
  
  
* quqd -- p=p178,q=p2
      quqd=p178(0)*p2(0)-p178(1)*p2(1)-p178(2)*p2(2)-p178(3)*p2(
     & 3)
* TWR0 -- qu=p178,qd=p2,v=0,a=rz2_178(0)%a,b=rz2_178(0)%b,cl=zcl(id2),ns
* um=0
      eps_0=-p178(2)*p2(3)+p2(2)*p178(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p178k0*p2(0)+p2k0*p178(0)
      rz2_178(0)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_178(0)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
* TWR0 -- qu=p178,qd=p2,v=1,a=rz2_178(1)%a,b=rz2_178(1)%b,cl=zcl(id2),ns
* um=0
      auxa=-quqd+p178k0*p2(1)+p2k0*p178(1)
      rz2_178(1)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_178(1)%b(1)=-zcl(id2)*(p2(2)+ceps_2)
* TWR0 -- qu=p178,qd=p2,v=2,a=rz2_178(2)%a,b=rz2_178(2)%b,cl=zcl(id2),ns
* um=0
      eps_0=-p178k0*p2(3)+p2k0*p178(3)
      ceps_0=eps_0*cim
      auxa=p178k0*p2(2)+p2k0*p178(2)
      rz2_178(2)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_178(2)%b(1)=-zcl(id2)*p2k0
* TWR0 -- qu=p178,qd=p2,v=3,a=rz2_178(3)%a,b=rz2_178(3)%b,cl=zcl(id2),ns
* um=0
      eps_0=p178k0*p2(2)-p2k0*p178(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p178k0*p2(3)+p2k0*p178(3)
      rz2_178(3)%a(2)=zcl(id2)*(auxa-ceps_0)
      rz2_178(3)%b(1)=-zcl(id2)*ceps_2
  
  
  
      if(ineutri(id2).ne.1) then
* TWR0 -- qu=p178,qd=p2,v=0,a=rf2_178(0)%a,b=rf2_178(0)%b,cl=fcl(id2),ns
* um=0
      eps_0=-p178(2)*p2(3)+p2(2)*p178(3)
      ceps_0=eps_0*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p178k0*p2(0)+p2k0*p178(0)
      rf2_178(0)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_178(0)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
* TWR0 -- qu=p178,qd=p2,v=1,a=rf2_178(1)%a,b=rf2_178(1)%b,cl=fcl(id2),ns
* um=0
      auxa=-quqd+p178k0*p2(1)+p2k0*p178(1)
      rf2_178(1)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_178(1)%b(1)=-fcl(id2)*(p2(2)+ceps_2)
* TWR0 -- qu=p178,qd=p2,v=2,a=rf2_178(2)%a,b=rf2_178(2)%b,cl=fcl(id2),ns
* um=0
      eps_0=-p178k0*p2(3)+p2k0*p178(3)
      ceps_0=eps_0*cim
      auxa=p178k0*p2(2)+p2k0*p178(2)
      rf2_178(2)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_178(2)%b(1)=-fcl(id2)*p2k0
* TWR0 -- qu=p178,qd=p2,v=3,a=rf2_178(3)%a,b=rf2_178(3)%b,cl=fcl(id2),ns
* um=0
      eps_0=p178k0*p2(2)-p2k0*p178(2)
      ceps_0=eps_0*cim
      ceps_2=p2k0*cim
      auxa=p178k0*p2(3)+p2k0*p178(3)
      rf2_178(3)%a(2)=fcl(id2)*(auxa-ceps_0)
      rf2_178(3)%b(1)=-fcl(id2)*ceps_2
  
      else
        do mu=0,3
          rf2_178(mu)%a(2) = czero
          rf2_178(mu)%b(1) = czero
        enddo
      endif
  
  
* quqd -- p=p312,q=p4
      quqd=p312(0)*p4(0)-p312(1)*p4(1)-p312(2)*p4(2)-p312(3)*p4(
     & 3)
* TWR0 -- qu=p312,qd=p4,v=0,a=rz4_312(0)%a,b=rz4_312(0)%b,cl=zcl(id4),ns
* um=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rz4_312(0)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(0)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
* TWR0 -- qu=p312,qd=p4,v=1,a=rz4_312(1)%a,b=rz4_312(1)%b,cl=zcl(id4),ns
* um=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rz4_312(1)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(1)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
* TWR0 -- qu=p312,qd=p4,v=2,a=rz4_312(2)%a,b=rz4_312(2)%b,cl=zcl(id4),ns
* um=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rz4_312(2)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(2)%b(1)=-zcl(id4)*p4k0
* TWR0 -- qu=p312,qd=p4,v=3,a=rz4_312(3)%a,b=rz4_312(3)%b,cl=zcl(id4),ns
* um=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p312k0*p4(3)+p4k0*p312(3)
      rz4_312(3)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(3)%b(1)=-zcl(id4)*ceps_2
  
  
  
      if(ineutri(id4).ne.1) then
* TWR0 -- qu=p312,qd=p4,v=0,a=rf4_312(0)%a,b=rf4_312(0)%b,cl=fcl(id4),ns
* um=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rf4_312(0)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(0)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
* TWR0 -- qu=p312,qd=p4,v=1,a=rf4_312(1)%a,b=rf4_312(1)%b,cl=fcl(id4),ns
* um=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rf4_312(1)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(1)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
* TWR0 -- qu=p312,qd=p4,v=2,a=rf4_312(2)%a,b=rf4_312(2)%b,cl=fcl(id4),ns
* um=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rf4_312(2)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(2)%b(1)=-fcl(id4)*p4k0
* TWR0 -- qu=p312,qd=p4,v=3,a=rf4_312(3)%a,b=rf4_312(3)%b,cl=fcl(id4),ns
* um=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p312k0*p4(3)+p4k0*p312(3)
      rf4_312(3)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(3)%b(1)=-fcl(id4)*ceps_2
  
      else
        do mu=0,3
          rf4_312(mu)%a(2) = czero
          rf4_312(mu)%b(1) = czero
        enddo
      endif
  
  
* quqd -- p=p356,q=p4
      quqd=p356(0)*p4(0)-p356(1)*p4(1)-p356(2)*p4(2)-p356(3)*p4(
     & 3)
* TWR0 -- qu=p356,qd=p4,v=0,a=rz4_356(0)%a,b=rz4_356(0)%b,cl=zcl(id4),ns
* um=0
      eps_0=-p356(2)*p4(3)+p4(2)*p356(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p356k0*p4(0)+p4k0*p356(0)
      rz4_356(0)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_356(0)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
* TWR0 -- qu=p356,qd=p4,v=1,a=rz4_356(1)%a,b=rz4_356(1)%b,cl=zcl(id4),ns
* um=0
      auxa=-quqd+p356k0*p4(1)+p4k0*p356(1)
      rz4_356(1)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_356(1)%b(1)=-zcl(id4)*(p4(2)+ceps_2)
* TWR0 -- qu=p356,qd=p4,v=2,a=rz4_356(2)%a,b=rz4_356(2)%b,cl=zcl(id4),ns
* um=0
      eps_0=-p356k0*p4(3)+p4k0*p356(3)
      ceps_0=eps_0*cim
      auxa=p356k0*p4(2)+p4k0*p356(2)
      rz4_356(2)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_356(2)%b(1)=-zcl(id4)*p4k0
* TWR0 -- qu=p356,qd=p4,v=3,a=rz4_356(3)%a,b=rz4_356(3)%b,cl=zcl(id4),ns
* um=0
      eps_0=p356k0*p4(2)-p4k0*p356(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p356k0*p4(3)+p4k0*p356(3)
      rz4_356(3)%a(2)=zcl(id4)*(auxa-ceps_0)
      rz4_356(3)%b(1)=-zcl(id4)*ceps_2
  
  
  
      if(ineutri(id4).ne.1) then
* TWR0 -- qu=p356,qd=p4,v=0,a=rf4_356(0)%a,b=rf4_356(0)%b,cl=fcl(id4),ns
* um=0
      eps_0=-p356(2)*p4(3)+p4(2)*p356(3)
      ceps_0=eps_0*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p356k0*p4(0)+p4k0*p356(0)
      rf4_356(0)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_356(0)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
* TWR0 -- qu=p356,qd=p4,v=1,a=rf4_356(1)%a,b=rf4_356(1)%b,cl=fcl(id4),ns
* um=0
      auxa=-quqd+p356k0*p4(1)+p4k0*p356(1)
      rf4_356(1)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_356(1)%b(1)=-fcl(id4)*(p4(2)+ceps_2)
* TWR0 -- qu=p356,qd=p4,v=2,a=rf4_356(2)%a,b=rf4_356(2)%b,cl=fcl(id4),ns
* um=0
      eps_0=-p356k0*p4(3)+p4k0*p356(3)
      ceps_0=eps_0*cim
      auxa=p356k0*p4(2)+p4k0*p356(2)
      rf4_356(2)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_356(2)%b(1)=-fcl(id4)*p4k0
* TWR0 -- qu=p356,qd=p4,v=3,a=rf4_356(3)%a,b=rf4_356(3)%b,cl=fcl(id4),ns
* um=0
      eps_0=p356k0*p4(2)-p4k0*p356(2)
      ceps_0=eps_0*cim
      ceps_2=p4k0*cim
      auxa=p356k0*p4(3)+p4k0*p356(3)
      rf4_356(3)%a(2)=fcl(id4)*(auxa-ceps_0)
      rf4_356(3)%b(1)=-fcl(id4)*ceps_2
  
      else
        do mu=0,3
          rf4_356(mu)%a(2) = czero
          rf4_356(mu)%b(1) = czero
        enddo
      endif
  
  
* quqd -- p=p534,q=p6
      quqd=p534(0)*p6(0)-p534(1)*p6(1)-p534(2)*p6(2)-p534(3)*p6(
     & 3)
* TWR0 -- qu=p534,qd=p6,v=0,a=rz6_534(0)%a,b=rz6_534(0)%b,cl=zcl(id6),ns
* um=0
      eps_0=-p534(2)*p6(3)+p6(2)*p534(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p534k0*p6(0)+p6k0*p534(0)
      rz6_534(0)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_534(0)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p534,qd=p6,v=1,a=rz6_534(1)%a,b=rz6_534(1)%b,cl=zcl(id6),ns
* um=0
      auxa=-quqd+p534k0*p6(1)+p6k0*p534(1)
      rz6_534(1)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_534(1)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p534,qd=p6,v=2,a=rz6_534(2)%a,b=rz6_534(2)%b,cl=zcl(id6),ns
* um=0
      eps_0=-p534k0*p6(3)+p6k0*p534(3)
      ceps_0=eps_0*cim
      auxa=p534k0*p6(2)+p6k0*p534(2)
      rz6_534(2)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_534(2)%b(1)=-zcl(id6)*p6k0
* TWR0 -- qu=p534,qd=p6,v=3,a=rz6_534(3)%a,b=rz6_534(3)%b,cl=zcl(id6),ns
* um=0
      eps_0=p534k0*p6(2)-p6k0*p534(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p534k0*p6(3)+p6k0*p534(3)
      rz6_534(3)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_534(3)%b(1)=-zcl(id6)*ceps_2
  
  
  
      if(ineutri(id6).ne.1) then
* TWR0 -- qu=p534,qd=p6,v=0,a=rf6_534(0)%a,b=rf6_534(0)%b,cl=fcl(id6),ns
* um=0
      eps_0=-p534(2)*p6(3)+p6(2)*p534(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p534k0*p6(0)+p6k0*p534(0)
      rf6_534(0)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_534(0)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p534,qd=p6,v=1,a=rf6_534(1)%a,b=rf6_534(1)%b,cl=fcl(id6),ns
* um=0
      auxa=-quqd+p534k0*p6(1)+p6k0*p534(1)
      rf6_534(1)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_534(1)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p534,qd=p6,v=2,a=rf6_534(2)%a,b=rf6_534(2)%b,cl=fcl(id6),ns
* um=0
      eps_0=-p534k0*p6(3)+p6k0*p534(3)
      ceps_0=eps_0*cim
      auxa=p534k0*p6(2)+p6k0*p534(2)
      rf6_534(2)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_534(2)%b(1)=-fcl(id6)*p6k0
* TWR0 -- qu=p534,qd=p6,v=3,a=rf6_534(3)%a,b=rf6_534(3)%b,cl=fcl(id6),ns
* um=0
      eps_0=p534k0*p6(2)-p6k0*p534(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p534k0*p6(3)+p6k0*p534(3)
      rf6_534(3)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_534(3)%b(1)=-fcl(id6)*ceps_2
  
      else
        do mu=0,3
          rf6_534(mu)%a(2) = czero
          rf6_534(mu)%b(1) = czero
        enddo
      endif
  
  
* quqd -- p=p578,q=p6
      quqd=p578(0)*p6(0)-p578(1)*p6(1)-p578(2)*p6(2)-p578(3)*p6(
     & 3)
* TWR0 -- qu=p578,qd=p6,v=0,a=rz6_578(0)%a,b=rz6_578(0)%b,cl=zcl(id6),ns
* um=0
      eps_0=-p578(2)*p6(3)+p6(2)*p578(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p578k0*p6(0)+p6k0*p578(0)
      rz6_578(0)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_578(0)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p578,qd=p6,v=1,a=rz6_578(1)%a,b=rz6_578(1)%b,cl=zcl(id6),ns
* um=0
      auxa=-quqd+p578k0*p6(1)+p6k0*p578(1)
      rz6_578(1)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_578(1)%b(1)=-zcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p578,qd=p6,v=2,a=rz6_578(2)%a,b=rz6_578(2)%b,cl=zcl(id6),ns
* um=0
      eps_0=-p578k0*p6(3)+p6k0*p578(3)
      ceps_0=eps_0*cim
      auxa=p578k0*p6(2)+p6k0*p578(2)
      rz6_578(2)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_578(2)%b(1)=-zcl(id6)*p6k0
* TWR0 -- qu=p578,qd=p6,v=3,a=rz6_578(3)%a,b=rz6_578(3)%b,cl=zcl(id6),ns
* um=0
      eps_0=p578k0*p6(2)-p6k0*p578(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p578k0*p6(3)+p6k0*p578(3)
      rz6_578(3)%a(2)=zcl(id6)*(auxa-ceps_0)
      rz6_578(3)%b(1)=-zcl(id6)*ceps_2
  
  
  
      if(ineutri(id6).ne.1) then
* TWR0 -- qu=p578,qd=p6,v=0,a=rf6_578(0)%a,b=rf6_578(0)%b,cl=fcl(id6),ns
* um=0
      eps_0=-p578(2)*p6(3)+p6(2)*p578(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p578k0*p6(0)+p6k0*p578(0)
      rf6_578(0)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_578(0)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p578,qd=p6,v=1,a=rf6_578(1)%a,b=rf6_578(1)%b,cl=fcl(id6),ns
* um=0
      auxa=-quqd+p578k0*p6(1)+p6k0*p578(1)
      rf6_578(1)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_578(1)%b(1)=-fcl(id6)*(p6(2)+ceps_2)
* TWR0 -- qu=p578,qd=p6,v=2,a=rf6_578(2)%a,b=rf6_578(2)%b,cl=fcl(id6),ns
* um=0
      eps_0=-p578k0*p6(3)+p6k0*p578(3)
      ceps_0=eps_0*cim
      auxa=p578k0*p6(2)+p6k0*p578(2)
      rf6_578(2)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_578(2)%b(1)=-fcl(id6)*p6k0
* TWR0 -- qu=p578,qd=p6,v=3,a=rf6_578(3)%a,b=rf6_578(3)%b,cl=fcl(id6),ns
* um=0
      eps_0=p578k0*p6(2)-p6k0*p578(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p578k0*p6(3)+p6k0*p578(3)
      rf6_578(3)%a(2)=fcl(id6)*(auxa-ceps_0)
      rf6_578(3)%b(1)=-fcl(id6)*ceps_2
  
      else
        do mu=0,3
          rf6_578(mu)%a(2) = czero
          rf6_578(mu)%b(1) = czero
        enddo
      endif
  
  
* quqd -- p=p712,q=p8
      quqd=p712(0)*p8(0)-p712(1)*p8(1)-p712(2)*p8(2)-p712(3)*p8(
     & 3)
* TWR0 -- qu=p712,qd=p8,v=0,a=rz8_712(0)%a,b=rz8_712(0)%b,cl=zcl(id8),ns
* um=0
      eps_0=-p712(2)*p8(3)+p8(2)*p712(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p712k0*p8(0)+p8k0*p712(0)
      rz8_712(0)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_712(0)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=1,a=rz8_712(1)%a,b=rz8_712(1)%b,cl=zcl(id8),ns
* um=0
      auxa=-quqd+p712k0*p8(1)+p8k0*p712(1)
      rz8_712(1)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_712(1)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=2,a=rz8_712(2)%a,b=rz8_712(2)%b,cl=zcl(id8),ns
* um=0
      eps_0=-p712k0*p8(3)+p8k0*p712(3)
      ceps_0=eps_0*cim
      auxa=p712k0*p8(2)+p8k0*p712(2)
      rz8_712(2)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_712(2)%b(1)=-zcl(id8)*p8k0
* TWR0 -- qu=p712,qd=p8,v=3,a=rz8_712(3)%a,b=rz8_712(3)%b,cl=zcl(id8),ns
* um=0
      eps_0=p712k0*p8(2)-p8k0*p712(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p712k0*p8(3)+p8k0*p712(3)
      rz8_712(3)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_712(3)%b(1)=-zcl(id8)*ceps_2
  
  
  
      if(ineutri(id8).ne.1) then
* TWR0 -- qu=p712,qd=p8,v=0,a=rf8_712(0)%a,b=rf8_712(0)%b,cl=fcl(id8),ns
* um=0
      eps_0=-p712(2)*p8(3)+p8(2)*p712(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p712k0*p8(0)+p8k0*p712(0)
      rf8_712(0)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_712(0)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=1,a=rf8_712(1)%a,b=rf8_712(1)%b,cl=fcl(id8),ns
* um=0
      auxa=-quqd+p712k0*p8(1)+p8k0*p712(1)
      rf8_712(1)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_712(1)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=2,a=rf8_712(2)%a,b=rf8_712(2)%b,cl=fcl(id8),ns
* um=0
      eps_0=-p712k0*p8(3)+p8k0*p712(3)
      ceps_0=eps_0*cim
      auxa=p712k0*p8(2)+p8k0*p712(2)
      rf8_712(2)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_712(2)%b(1)=-fcl(id8)*p8k0
* TWR0 -- qu=p712,qd=p8,v=3,a=rf8_712(3)%a,b=rf8_712(3)%b,cl=fcl(id8),ns
* um=0
      eps_0=p712k0*p8(2)-p8k0*p712(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p712k0*p8(3)+p8k0*p712(3)
      rf8_712(3)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_712(3)%b(1)=-fcl(id8)*ceps_2
  
      else
        do mu=0,3
          rf8_712(mu)%a(2) = czero
          rf8_712(mu)%b(1) = czero
        enddo
      endif
  
  
* quqd -- p=p756,q=p8
      quqd=p756(0)*p8(0)-p756(1)*p8(1)-p756(2)*p8(2)-p756(3)*p8(
     & 3)
* TWR0 -- qu=p756,qd=p8,v=0,a=rz8_756(0)%a,b=rz8_756(0)%b,cl=zcl(id8),ns
* um=0
      eps_0=-p756(2)*p8(3)+p8(2)*p756(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p756k0*p8(0)+p8k0*p756(0)
      rz8_756(0)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_756(0)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p756,qd=p8,v=1,a=rz8_756(1)%a,b=rz8_756(1)%b,cl=zcl(id8),ns
* um=0
      auxa=-quqd+p756k0*p8(1)+p8k0*p756(1)
      rz8_756(1)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_756(1)%b(1)=-zcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p756,qd=p8,v=2,a=rz8_756(2)%a,b=rz8_756(2)%b,cl=zcl(id8),ns
* um=0
      eps_0=-p756k0*p8(3)+p8k0*p756(3)
      ceps_0=eps_0*cim
      auxa=p756k0*p8(2)+p8k0*p756(2)
      rz8_756(2)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_756(2)%b(1)=-zcl(id8)*p8k0
* TWR0 -- qu=p756,qd=p8,v=3,a=rz8_756(3)%a,b=rz8_756(3)%b,cl=zcl(id8),ns
* um=0
      eps_0=p756k0*p8(2)-p8k0*p756(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p756k0*p8(3)+p8k0*p756(3)
      rz8_756(3)%a(2)=zcl(id8)*(auxa-ceps_0)
      rz8_756(3)%b(1)=-zcl(id8)*ceps_2
  
  
  
      if(ineutri(id8).ne.1) then
* TWR0 -- qu=p756,qd=p8,v=0,a=rf8_756(0)%a,b=rf8_756(0)%b,cl=fcl(id8),ns
* um=0
      eps_0=-p756(2)*p8(3)+p8(2)*p756(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p756k0*p8(0)+p8k0*p756(0)
      rf8_756(0)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_756(0)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p756,qd=p8,v=1,a=rf8_756(1)%a,b=rf8_756(1)%b,cl=fcl(id8),ns
* um=0
      auxa=-quqd+p756k0*p8(1)+p8k0*p756(1)
      rf8_756(1)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_756(1)%b(1)=-fcl(id8)*(p8(2)+ceps_2)
* TWR0 -- qu=p756,qd=p8,v=2,a=rf8_756(2)%a,b=rf8_756(2)%b,cl=fcl(id8),ns
* um=0
      eps_0=-p756k0*p8(3)+p8k0*p756(3)
      ceps_0=eps_0*cim
      auxa=p756k0*p8(2)+p8k0*p756(2)
      rf8_756(2)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_756(2)%b(1)=-fcl(id8)*p8k0
* TWR0 -- qu=p756,qd=p8,v=3,a=rf8_756(3)%a,b=rf8_756(3)%b,cl=fcl(id8),ns
* um=0
      eps_0=p756k0*p8(2)-p8k0*p756(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p756k0*p8(3)+p8k0*p756(3)
      rf8_756(3)%a(2)=fcl(id8)*(auxa-ceps_0)
      rf8_756(3)%b(1)=-fcl(id8)*ceps_2
  
      else
        do mu=0,3
          rf8_756(mu)%a(2) = czero
          rf8_756(mu)%b(1) = czero
        enddo
      endif
  
  
  
* compute all subdiagrams with a Z or a gamma "decaying" to 4 fermions: 
*                  czijkl(mu)  and cfijkl(mu) (ijkl integers)           
*                                                                       
*               _____ i                       _____ i                   
*              |  |__ j                      |  |__ j                   
*    (mu) --Z--|  |__ k        (mu) --gamma--|  |__ k                   
*              |__|__ l                      |__|__ l                   
*                                                                       
*                                                                       
  
*                                                                       
*    First all diagrams with a Z of the type                            
*                                                                       
*                              |_W__/                                   
*                    (mu) _Z___|    \                                   
*                              |                                        
  
*                 first four without sum                                
  
      do mu=0,3
* TLTR0_W -- aa=cz1234(mu),a1=l1_34%a,c1=l1_34%c,a2=rz2_134(mu)%a,b2=rz2
* _134(mu)%b,prq=p134q,bef=,aft=
      cz1234(mu)=(l1_34%c(2)*p134q*rz2_134(mu)%b(1)+l1_34%a(2)*r
     & z2_134(mu)%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz1278(mu),a1=l1_78%a,c1=l1_78%c,a2=rz2_178(mu)%a,b2=rz2
* _178(mu)%b,prq=p178q,bef=,aft=
      cz1278(mu)=(l1_78%c(2)*p178q*rz2_178(mu)%b(1)+l1_78%a(2)*r
     & z2_178(mu)%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz5634(mu),a1=l3_56%a,c1=l3_56%c,a2=rz4_356(mu)%a,b2=rz4
* _356(mu)%b,prq=p356q,bef=,aft=
      cz5634(mu)=(l3_56%c(2)*p356q*rz4_356(mu)%b(1)+l3_56%a(2)*r
     & z4_356(mu)%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz5678(mu),a1=l5_78%a,c1=l5_78%c,a2=rz6_578(mu)%a,b2=rz6
* _578(mu)%b,prq=p578q,bef=,aft=
      cz5678(mu)=(l5_78%c(2)*p578q*rz6_578(mu)%b(1)+l5_78%a(2)*r
     & z6_578(mu)%a(2))
      end do
  
*                 second four  with  sum                                
  
      do mu=0,3
* TLTR0_W -- aa=cz1234(mu),a1=l3_12%a,c1=l3_12%c,a2=rz4_312(mu)%a,b2=rz4
* _312(mu)%b,prq=p312q,bef=cz1234(mu)+,aft=
      cz1234(mu)=cz1234(mu)+(l3_12%c(2)*p312q*rz4_312(mu)%b(1)+l
     & 3_12%a(2)*rz4_312(mu)%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz5634(mu),a1=l5_34%a,c1=l5_34%c,a2=rz6_534(mu)%a,b2=rz6
* _534(mu)%b,prq=p534q,bef=cz5634(mu)+,aft=
      cz5634(mu)=cz5634(mu)+(l5_34%c(2)*p534q*rz6_534(mu)%b(1)+l
     & 5_34%a(2)*rz6_534(mu)%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz1278(mu),a1=l7_12%a,c1=l7_12%c,a2=rz8_712(mu)%a,b2=rz8
* _712(mu)%b,prq=p712q,bef=cz1278(mu)+,aft=
      cz1278(mu)=cz1278(mu)+(l7_12%c(2)*p712q*rz8_712(mu)%b(1)+l
     & 7_12%a(2)*rz8_712(mu)%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz5678(mu),a1=l7_56%a,c1=l7_56%c,a2=rz8_756(mu)%a,b2=rz8
* _756(mu)%b,prq=p756q,bef=cz5678(mu)+,aft=
      cz5678(mu)=cz5678(mu)+(l7_56%c(2)*p756q*rz8_756(mu)%b(1)+l
     & 7_56%a(2)*rz8_756(mu)%a(2))
      end do
  
  
*    Then diagrams with a gamma of the type                             
*                                                                       
*                          |_W__/                                       
*             (mu) _gamma__|    \                                       
*                          |                                            
*                                                                       
**QCD                                                                   
*    and their QCD counterparts                                         
*                                                                       
*                          |_W__/                                       
*             (mu)  ~~~g~~~|    \                                       
*                          |                                            
  
*                 first four                                            
  
      do mu=0,3
* TLTR0_W -- aa=cf1234(mu),a1=l1_34%a,c1=l1_34%c,a2=rf2_134(mu)%a,b2=rf2
* _134(mu)%b,prq=p134q,bef=,aft=
      cf1234(mu)=(l1_34%c(2)*p134q*rf2_134(mu)%b(1)+l1_34%a(2)*r
     & f2_134(mu)%a(2))
      end do
      if(ilept(id2).ne.1) then
        do m=0,3
          cg1234(m)=cf1234(m)/fcl(id2)
        enddo
      endif
      do mu=0,3
* TLTR0_W -- aa=cf1278(mu),a1=l1_78%a,c1=l1_78%c,a2=rf2_178(mu)%a,b2=rf2
* _178(mu)%b,prq=p178q,bef=,aft=
      cf1278(mu)=(l1_78%c(2)*p178q*rf2_178(mu)%b(1)+l1_78%a(2)*r
     & f2_178(mu)%a(2))
      end do
      if(ilept(id2).ne.1) then
        do m=0,3
          cg1278(m)=cf1278(m)/fcl(id2)
        enddo
      endif
      do mu=0,3
* TLTR0_W -- aa=cf5634(mu),a1=l3_56%a,c1=l3_56%c,a2=rf4_356(mu)%a,b2=rf4
* _356(mu)%b,prq=p356q,bef=,aft=
      cf5634(mu)=(l3_56%c(2)*p356q*rf4_356(mu)%b(1)+l3_56%a(2)*r
     & f4_356(mu)%a(2))
      end do
      if(ilept(id4).ne.1) then
        do m=0,3
          cg3456(m)=cf5634(m)/fcl(id4)
        enddo
      endif
      do mu=0,3
* TLTR0_W -- aa=cf5678(mu),a1=l5_78%a,c1=l5_78%c,a2=rf6_578(mu)%a,b2=rf6
* _578(mu)%b,prq=p578q,bef=,aft=
      cf5678(mu)=(l5_78%c(2)*p578q*rf6_578(mu)%b(1)+l5_78%a(2)*r
     & f6_578(mu)%a(2))
      end do
      if(ilept(id6).ne.1) then
        do m=0,3
          cg5678(m)=cf5678(m)/fcl(id6)
        enddo
      endif
  
*                 second four                                           
  
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=l3_12%a,c1=l3_12%c,a2=rf4_312(mu)%a,b2=rf4
* _312(mu)%b,prq=p312q,bef=,aft=
      cftemp(mu)=(l3_12%c(2)*p312q*rf4_312(mu)%b(1)+l3_12%a(2)*r
     & f4_312(mu)%a(2))
      end do
      do m=0,3
        cf1234(m)=cf1234(m)+cftemp(m)
        if(ilept(id4).ne.1) then
          cg3412(m)=cftemp(m)/fcl(id4)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=l5_34%a,c1=l5_34%c,a2=rf6_534(mu)%a,b2=rf6
* _534(mu)%b,prq=p534q,bef=,aft=
      cftemp(mu)=(l5_34%c(2)*p534q*rf6_534(mu)%b(1)+l5_34%a(2)*r
     & f6_534(mu)%a(2))
      end do
      do m=0,3
        cf5634(m)=cf5634(m)+cftemp(m)
        if(ilept(id6).ne.1) then
          cg5634(m)=cftemp(m)/fcl(id6)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=l7_12%a,c1=l7_12%c,a2=rf8_712(mu)%a,b2=rf8
* _712(mu)%b,prq=p712q,bef=,aft=
      cftemp(mu)=(l7_12%c(2)*p712q*rf8_712(mu)%b(1)+l7_12%a(2)*r
     & f8_712(mu)%a(2))
      end do
      do m=0,3
        cf1278(m)=cf1278(m)+cftemp(m)
        if(ilept(id8).ne.1) then
          cg7812(m)=cftemp(m)/fcl(id8)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=l7_56%a,c1=l7_56%c,a2=rf8_756(mu)%a,b2=rf8
* _756(mu)%b,prq=p756q,bef=,aft=
      cftemp(mu)=(l7_56%c(2)*p756q*rf8_756(mu)%b(1)+l7_56%a(2)*r
     & f8_756(mu)%a(2))
      end do
      do m=0,3
        cf5678(m)=cf5678(m)+cftemp(m)
        if(ilept(id8).ne.1) then
          cg7856(m)=cftemp(m)/fcl(id8)
        endif
      enddo
  
  
*    Then those with a Z of the type                                    
*                                                                       
*                  (mu) _Z___|                                          
*                            |_W__/                                     
*                            |    \                                     
*                                                                       
  
*                        all 8 with sum                                 
  
      do mu=0,3
* TLTR0_W -- aa=cz1234(mu),a1=lz1_234(mu)%a,c1=lz1_234(mu)%c,a2=r2_34%a,
* b2=r2_34%b,prq=p234q,bef=cz1234(mu)+,aft=
      cz1234(mu)=cz1234(mu)+(lz1_234(mu)%c(2)*p234q*r2_34%b(1)+l
     & z1_234(mu)%a(2)*r2_34%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz1278(mu),a1=lz1_278(mu)%a,c1=lz1_278(mu)%c,a2=r2_78%a,
* b2=r2_78%b,prq=p278q,bef=cz1278(mu)+,aft=
      cz1278(mu)=cz1278(mu)+(lz1_278(mu)%c(2)*p278q*r2_78%b(1)+l
     & z1_278(mu)%a(2)*r2_78%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz1234(mu),a1=lz3_412(mu)%a,c1=lz3_412(mu)%c,a2=r4_12%a,
* b2=r4_12%b,prq=p412q,bef=cz1234(mu)+,aft=
      cz1234(mu)=cz1234(mu)+(lz3_412(mu)%c(2)*p412q*r4_12%b(1)+l
     & z3_412(mu)%a(2)*r4_12%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz5634(mu),a1=lz3_456(mu)%a,c1=lz3_456(mu)%c,a2=r4_56%a,
* b2=r4_56%b,prq=p456q,bef=cz5634(mu)+,aft=
      cz5634(mu)=cz5634(mu)+(lz3_456(mu)%c(2)*p456q*r4_56%b(1)+l
     & z3_456(mu)%a(2)*r4_56%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz5634(mu),a1=lz5_634(mu)%a,c1=lz5_634(mu)%c,a2=r6_34%a,
* b2=r6_34%b,prq=p634q,bef=cz5634(mu)+,aft=
      cz5634(mu)=cz5634(mu)+(lz5_634(mu)%c(2)*p634q*r6_34%b(1)+l
     & z5_634(mu)%a(2)*r6_34%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz5678(mu),a1=lz5_678(mu)%a,c1=lz5_678(mu)%c,a2=r6_78%a,
* b2=r6_78%b,prq=p678q,bef=cz5678(mu)+,aft=
      cz5678(mu)=cz5678(mu)+(lz5_678(mu)%c(2)*p678q*r6_78%b(1)+l
     & z5_678(mu)%a(2)*r6_78%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz1278(mu),a1=lz7_812(mu)%a,c1=lz7_812(mu)%c,a2=r8_12%a,
* b2=r8_12%b,prq=p812q,bef=cz1278(mu)+,aft=
      cz1278(mu)=cz1278(mu)+(lz7_812(mu)%c(2)*p812q*r8_12%b(1)+l
     & z7_812(mu)%a(2)*r8_12%a(2))
      end do
      do mu=0,3
* TLTR0_W -- aa=cz5678(mu),a1=lz7_856(mu)%a,c1=lz7_856(mu)%c,a2=r8_56%a,
* b2=r8_56%b,prq=p856q,bef=cz5678(mu)+,aft=
      cz5678(mu)=cz5678(mu)+(lz7_856(mu)%c(2)*p856q*r8_56%b(1)+l
     & z7_856(mu)%a(2)*r8_56%a(2))
      end do
  
  
*    Finally those with a gamma of the type                             
*                                                                       
*                (mu) _gamma__|                                         
*                             |_W__/                                    
*                             |    \                                    
*                                                                       
*                                                                       
**QCD                                                                   
*    and their QCD counterparts                                         
*                                                                       
*                   ~~~g~~~|                                            
*             (mu)         |_W__/                                       
*                          |    \                                       
  
*                        all 8 with sum                                 
  
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=lf1_234(mu)%a,c1=lf1_234(mu)%c,a2=r2_34%a,
* b2=r2_34%b,prq=p234q,bef=,aft=
      cftemp(mu)=(lf1_234(mu)%c(2)*p234q*r2_34%b(1)+lf1_234(mu)%
     & a(2)*r2_34%a(2))
      end do
      do m=0,3
        cf1234(m)=cf1234(m)+cftemp(m)
        if(ilept(id1).ne.1) then
          cg1234(m)=cg1234(m)+cftemp(m)/fcl(id1)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=lf1_278(mu)%a,c1=lf1_278(mu)%c,a2=r2_78%a,
* b2=r2_78%b,prq=p278q,bef=,aft=
      cftemp(mu)=(lf1_278(mu)%c(2)*p278q*r2_78%b(1)+lf1_278(mu)%
     & a(2)*r2_78%a(2))
      end do
      do m=0,3
        cf1278(m)=cf1278(m)+cftemp(m)
        if(ilept(id1).ne.1) then
          cg1278(m)=cg1278(m)+cftemp(m)/fcl(id1)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=lf3_412(mu)%a,c1=lf3_412(mu)%c,a2=r4_12%a,
* b2=r4_12%b,prq=p412q,bef=,aft=
      cftemp(mu)=(lf3_412(mu)%c(2)*p412q*r4_12%b(1)+lf3_412(mu)%
     & a(2)*r4_12%a(2))
      end do
      do m=0,3
        cf1234(m)=cf1234(m)+cftemp(m)
        if(ilept(id3).ne.1) then
          cg3412(m)=cg3412(m)+cftemp(m)/fcl(id3)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=lf3_456(mu)%a,c1=lf3_456(mu)%c,a2=r4_56%a,
* b2=r4_56%b,prq=p456q,bef=,aft=
      cftemp(mu)=(lf3_456(mu)%c(2)*p456q*r4_56%b(1)+lf3_456(mu)%
     & a(2)*r4_56%a(2))
      end do
      do m=0,3
        cf5634(m)=cf5634(m)+cftemp(m)
        if(ilept(id3).ne.1) then
          cg3456(m)=cg3456(m)+cftemp(m)/fcl(id3)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=lf5_634(mu)%a,c1=lf5_634(mu)%c,a2=r6_34%a,
* b2=r6_34%b,prq=p634q,bef=,aft=
      cftemp(mu)=(lf5_634(mu)%c(2)*p634q*r6_34%b(1)+lf5_634(mu)%
     & a(2)*r6_34%a(2))
      end do
      do m=0,3
        cf5634(m)=cf5634(m)+cftemp(m)
        if(ilept(id5).ne.1) then
          cg5634(m)=cg5634(m)+cftemp(m)/fcl(id5)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=lf5_678(mu)%a,c1=lf5_678(mu)%c,a2=r6_78%a,
* b2=r6_78%b,prq=p678q,bef=,aft=
      cftemp(mu)=(lf5_678(mu)%c(2)*p678q*r6_78%b(1)+lf5_678(mu)%
     & a(2)*r6_78%a(2))
      end do
      do m=0,3
        cf5678(m)=cf5678(m)+cftemp(m)
        if(ilept(id5).ne.1) then
          cg5678(m)=cg5678(m)+cftemp(m)/fcl(id5)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=lf7_812(mu)%a,c1=lf7_812(mu)%c,a2=r8_12%a,
* b2=r8_12%b,prq=p812q,bef=,aft=
      cftemp(mu)=(lf7_812(mu)%c(2)*p812q*r8_12%b(1)+lf7_812(mu)%
     & a(2)*r8_12%a(2))
      end do
      do m=0,3
        cf1278(m)=cf1278(m)+cftemp(m)
        if(ilept(id7).ne.1) then
          cg7812(m)=cg7812(m)+cftemp(m)/fcl(id7)
        endif
      enddo
      do mu=0,3
* TLTR0_W -- aa=cftemp(mu),a1=lf7_856(mu)%a,c1=lf7_856(mu)%c,a2=r8_56%a,
* b2=r8_56%b,prq=p856q,bef=,aft=
      cftemp(mu)=(lf7_856(mu)%c(2)*p856q*r8_56%b(1)+lf7_856(mu)%
     & a(2)*r8_56%a(2))
      end do
      do m=0,3
        cf5678(m)=cf5678(m)+cftemp(m)
        if(ilept(id7).ne.1) then
          cg7856(m)=cg7856(m)+cftemp(m)/fcl(id7)
        endif
      enddo
  
*    Then those with triple vertex                                      
*                                                                       
*                    /_                            /_                   
*                  W/                            W/                     
*         (mu)_Z___/          and    (mu)_gamma__/                      
*                  \                             \                      
*                  W\_                           W\_                    
*                    \                             \                    
*                                                                       
  
      do mu=0,3
        p1234(mu)=p12(mu)+p34(mu)
      enddo
      do mu=0,3
        p1278(mu)=p12(mu)+p78(mu)
      enddo
  
*** triple vertex -- pfz(mu)=(-p1234(mu)),pwm(mu)=p34(mu),pwp(mu)=p12(mu
* ),efz=#,ewm=cw34,ewp=cw12,res=ctrip1234(mu#)
      do mu=0,3
      vfz(mu)=p34(mu)-p12(mu)
      vwm(mu)=p12(mu)-(-p1234(mu))
      vwp(mu)=(-p1234(mu))-p34(mu)
      end do !mu
* vwm%ewm
* p.q -- p.q=cw34%v,p=cw34%e,q=vwm,bef=,aft=
      cw34%v=(cw34%e(0)*vwm(0)-cw34%e(1)*vwm(1)-cw34%e(2)*vwm(2)
     & -cw34%e(3)*vwm(3))
* vwp%ewp
* p.q -- p.q=cw12%v,p=cw12%e,q=vwp,bef=,aft=
      cw12%v=(cw12%e(0)*vwp(0)-cw12%e(1)*vwp(1)-cw12%e(2)*vwp(2)
     & -cw12%e(3)*vwp(3))
* ewm%ewp
* p.q -- p.q=caux,p=cw34%e,q=cw12%e,bef=,aft=
      caux=(cw34%e(0)*cw12%e(0)-cw34%e(1)*cw12%e(1)-cw34%e(2)*cw
     & 12%e(2)-cw34%e(3)*cw12%e(3))
      do mu=0,3
      ctrip1234(mu)=vfz(mu)*caux+cw34%v*cw12%e(mu)+cw12%v*cw34%e
     & (mu)
      end do
*** triple vertex -- pfz(mu)=(-p1278(mu)),pwm(mu)=p78(mu),pwp(mu)=p12(mu
* ),efz=#,ewm=cw78,ewp=cw12,res=ctrip1278(mu#)
      do mu=0,3
      vfz(mu)=p78(mu)-p12(mu)
      vwm(mu)=p12(mu)-(-p1278(mu))
      vwp(mu)=(-p1278(mu))-p78(mu)
      end do !mu
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* vwp%ewp
* p.q -- p.q=cw12%v,p=cw12%e,q=vwp,bef=,aft=
      cw12%v=(cw12%e(0)*vwp(0)-cw12%e(1)*vwp(1)-cw12%e(2)*vwp(2)
     & -cw12%e(3)*vwp(3))
* ewm%ewp
* p.q -- p.q=caux,p=cw78%e,q=cw12%e,bef=,aft=
      caux=(cw78%e(0)*cw12%e(0)-cw78%e(1)*cw12%e(1)-cw78%e(2)*cw
     & 12%e(2)-cw78%e(3)*cw12%e(3))
      do mu=0,3
      ctrip1278(mu)=vfz(mu)*caux+cw78%v*cw12%e(mu)+cw12%v*cw78%e
     & (mu)
      end do
*** triple vertex -- pfz(mu)=p1278(mu),pwm(mu)=p34(mu),pwp(mu)=p56(mu),e
* fz=#,ewm=cw34,ewp=cw56,res=ctrip5634(mu#)
      do mu=0,3
      vfz(mu)=p34(mu)-p56(mu)
      vwm(mu)=p56(mu)-p1278(mu)
      vwp(mu)=p1278(mu)-p34(mu)
      end do !mu
* vwm%ewm
* p.q -- p.q=cw34%v,p=cw34%e,q=vwm,bef=,aft=
      cw34%v=(cw34%e(0)*vwm(0)-cw34%e(1)*vwm(1)-cw34%e(2)*vwm(2)
     & -cw34%e(3)*vwm(3))
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* ewm%ewp
* p.q -- p.q=caux,p=cw34%e,q=cw56%e,bef=,aft=
      caux=(cw34%e(0)*cw56%e(0)-cw34%e(1)*cw56%e(1)-cw34%e(2)*cw
     & 56%e(2)-cw34%e(3)*cw56%e(3))
      do mu=0,3
      ctrip5634(mu)=vfz(mu)*caux+cw34%v*cw56%e(mu)+cw56%v*cw34%e
     & (mu)
      end do
*** triple vertex -- pfz(mu)=p1234(mu),pwm(mu)=p78(mu),pwp(mu)=p56(mu),e
* fz=#,ewm=cw78,ewp=cw56,res=ctrip5678(mu#)
      do mu=0,3
      vfz(mu)=p78(mu)-p56(mu)
      vwm(mu)=p56(mu)-p1234(mu)
      vwp(mu)=p1234(mu)-p78(mu)
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
      ctrip5678(mu)=vfz(mu)*caux+cw78%v*cw56%e(mu)+cw56%v*cw78%e
     & (mu)
      end do
  
      do mu=0,3
      cz1234(mu)=rcotw*ctrip1234(mu)+cz1234(mu)
      enddo
      do mu=0,3
      cf1234(mu)= ctrip1234(mu)+cf1234(mu)
      enddo
      do mu=0,3
      cz1278(mu)=rcotw*ctrip1278(mu)+cz1278(mu)
      enddo
      do mu=0,3
      cf1278(mu)= ctrip1278(mu)+cf1278(mu)
      enddo
      do mu=0,3
      cz5634(mu)=rcotw*ctrip5634(mu)+cz5634(mu)
      enddo
      do mu=0,3
      cf5634(mu)= ctrip5634(mu)+cf5634(mu)
      enddo
      do mu=0,3
      cz5678(mu)=rcotw*ctrip5678(mu)+cz5678(mu)
      enddo
      do mu=0,3
      cf5678(mu)= ctrip5678(mu)+cf5678(mu)
      enddo
  
* compute all subdiagrams with a higgs decaying to 4 fermions           
*                                                                       
*               /_                                                      
*            w+/                                                        
*       ___h__/                                                         
*             \                                                         
*            w-\_                                                       
*               \                                                       
*                                                                       
  
* We factorize out from the amplitude the modulus of the electric       
* charge from every coupling, hence instead of g  we have               
* coupling g/|e|=1/sw  The complete coupling is g*rmw, i.e. rmw/sw      
  
*rmh < 0 in our convention means no higgs coupling                      
      if (rmh.ge.0.d0) then
  
* p.q -- p.q=ch1234,p=cw12%e,q=cw34%e,bef=,aft=*rmw/sw
      ch1234=(cw12%e(0)*cw34%e(0)-cw12%e(1)*cw34%e(1)-cw12%e(2)*
     & cw34%e(2)-cw12%e(3)*cw34%e(3))*rmw/sw
* p.q -- p.q=ch1278,p=cw12%e,q=cw78%e,bef=,aft=*rmw/sw
      ch1278=(cw12%e(0)*cw78%e(0)-cw12%e(1)*cw78%e(1)-cw12%e(2)*
     & cw78%e(2)-cw12%e(3)*cw78%e(3))*rmw/sw
* p.q -- p.q=ch5634,p=cw56%e,q=cw34%e,bef=,aft=*rmw/sw
      ch5634=(cw56%e(0)*cw34%e(0)-cw56%e(1)*cw34%e(1)-cw56%e(2)*
     & cw34%e(2)-cw56%e(3)*cw34%e(3))*rmw/sw
* p.q -- p.q=ch5678,p=cw56%e,q=cw78%e,bef=,aft=*rmw/sw
      ch5678=(cw56%e(0)*cw78%e(0)-cw56%e(1)*cw78%e(1)-cw56%e(2)*
     & cw78%e(2)-cw56%e(3)*cw78%e(3))*rmw/sw
  
      endif
  
*  compute 8 diagrams                                                   
*          __                                                           
*            |_W__/                                                     
*            |    \                                                     
*            |                                                          
*            |_W__/                                                     
*            |    \                                                     
*            |                                                          
*            |_W__/                                                     
*          __|    \                                                     
*                                                                       
*                                                                       
  
* TLTR0_W -- aa=cresUD_3w,a1=l1_3456%a,c1=l1_3456%c,a2=r2_78%a,b2=r2_78%
* b,prq=p278q,bef=,aft=
      cresUD_3w=(l1_3456%c(2)*p278q*r2_78%b(1)+l1_3456%a(2)*r2_7
     & 8%a(2))
* TLTR0_W -- aa=cresUD_3w,a1=l1_7856%a,c1=l1_7856%c,a2=r2_34%a,b2=r2_34%
* b,prq=p234q,bef=cresUD_3w+,aft=
      cresUD_3w=cresUD_3w+(l1_7856%c(2)*p234q*r2_34%b(1)+l1_7856
     & %a(2)*r2_34%a(2))
* TLTR0_W -- aa=cresUD_3w,a1=l3_1278%a,c1=l3_1278%c,a2=r4_56%a,b2=r4_56%
* b,prq=p456q,bef=cresUD_3w+,aft=
      cresUD_3w=cresUD_3w+(l3_1278%c(2)*p456q*r4_56%b(1)+l3_1278
     & %a(2)*r4_56%a(2))
* TLTR0_W -- aa=cresUD_3w,a1=l3_5678%a,c1=l3_5678%c,a2=r4_12%a,b2=r4_12%
* b,prq=p412q,bef=cresUD_3w+,aft=
      cresUD_3w=cresUD_3w+(l3_5678%c(2)*p412q*r4_12%b(1)+l3_5678
     & %a(2)*r4_12%a(2))
* TLTR0_W -- aa=cresUD_3w,a1=l5_3412%a,c1=l5_3412%c,a2=r6_78%a,b2=r6_78%
* b,prq=p678q,bef=cresUD_3w+,aft=
      cresUD_3w=cresUD_3w+(l5_3412%c(2)*p678q*r6_78%b(1)+l5_3412
     & %a(2)*r6_78%a(2))
* TLTR0_W -- aa=cresUD_3w,a1=l5_7812%a,c1=l5_7812%c,a2=r6_34%a,b2=r6_34%
* b,prq=p634q,bef=cresUD_3w+,aft=
      cresUD_3w=cresUD_3w+(l5_7812%c(2)*p634q*r6_34%b(1)+l5_7812
     & %a(2)*r6_34%a(2))
* TLTR0_W -- aa=cresUD_3w,a1=l7_1234%a,c1=l7_1234%c,a2=r8_56%a,b2=r8_56%
* b,prq=p856q,bef=cresUD_3w+,aft=
      cresUD_3w=cresUD_3w+(l7_1234%c(2)*p856q*r8_56%b(1)+l7_1234
     & %a(2)*r8_56%a(2))
* TLTR0_W -- aa=cresUD_3w,a1=l7_5634%a,c1=l7_5634%c,a2=r8_12%a,b2=r8_12%
* b,prq=p812q,bef=cresUD_3w+,aft=
      cresUD_3w=cresUD_3w+(l7_5634%c(2)*p812q*r8_12%b(1)+l7_5634
     & %a(2)*r8_12%a(2))
  
* compute the 50 (some may be zero with neutrinos) with a photon        
*  propagator connecting 2 w's on both sides                            
  
* p.q -- p.q=p1234q,p=p1234,q=p1234,bef=,aft=
      p1234q=(p1234(0)*p1234(0)-p1234(1)*p1234(1)-p1234(2)*p1234
     & (2)-p1234(3)*p1234(3))
* p.q -- p.q=p1278q,p=p1278,q=p1278,bef=,aft=
      p1278q=(p1278(0)*p1278(0)-p1278(1)*p1278(1)-p1278(2)*p1278
     & (2)-p1278(3)*p1278(3))
* p.q -- p.q=cres2w_f_2w,p=cf1234,q=cf5678,bef=-,aft=/p1234q
      cres2w_f_2w=-(cf1234(0)*cf5678(0)-cf1234(1)*cf5678(1)-cf12
     & 34(2)*cf5678(2)-cf1234(3)*cf5678(3))/p1234q
* p.q -- p.q=cres2w_f_2w,p=cf1278,q=cf5634,bef=cres2w_f_2w-,aft=/p1278q
      cres2w_f_2w=cres2w_f_2w-(cf1278(0)*cf5634(0)-cf1278(1)*cf5
     & 634(1)-cf1278(2)*cf5634(2)-cf1278(3)*cf5634(3))/p1278q
  
  
* compute the 50 diagrams with a z propagator connecting 2 w's on       
*  both sides                                                           
  
* p.q -- p.q=cp1234dotcz1234,p=p1234,q=cz1234,bef=,aft=
      cp1234dotcz1234=(p1234(0)*cz1234(0)-p1234(1)*cz1234(1)-p12
     & 34(2)*cz1234(2)-p1234(3)*cz1234(3))
* p.q -- p.q=cp1234dotcz5678,p=p1234,q=cz5678,bef=,aft=
      cp1234dotcz5678=(p1234(0)*cz5678(0)-p1234(1)*cz5678(1)-p12
     & 34(2)*cz5678(2)-p1234(3)*cz5678(3))
* p.q -- p.q=cp1278dotcz1278,p=p1278,q=cz1278,bef=,aft=
      cp1278dotcz1278=(p1278(0)*cz1278(0)-p1278(1)*cz1278(1)-p12
     & 78(2)*cz1278(2)-p1278(3)*cz1278(3))
* p.q -- p.q=cp1278dotcz5634,p=p1278,q=cz5634,bef=,aft=
      cp1278dotcz5634=(p1278(0)*cz5634(0)-p1278(1)*cz5634(1)-p12
     & 78(2)*cz5634(2)-p1278(3)*cz5634(3))
  
* p.q -- p.q=cres2w_z_2w,p=cz1234,q=cz5678,bef=-,aft=/(p1234q-cmz2)
      cres2w_z_2w=-(cz1234(0)*cz5678(0)-cz1234(1)*cz5678(1)-cz12
     & 34(2)*cz5678(2)-cz1234(3)*cz5678(3))/(p1234q-cmz2)
* p.q -- p.q=cres2w_z_2w,p=cz1278,q=cz5634,bef=cres2w_z_2w-,aft=/(p1278q
* -cmz2)
      cres2w_z_2w=cres2w_z_2w-(cz1278(0)*cz5634(0)-cz1278(1)*cz5
     & 634(1)-cz1278(2)*cz5634(2)-cz1278(3)*cz5634(3))/(p1278q-c
     & mz2)
  
      cres2w_z_2w=cres2w_z_2w+
     &   (cp1234dotcz1234*cp1234dotcz5678)/cmz2/(p1234q-cmz2)
     &  +(cp1278dotcz1278*cp1278dotcz5634)/cmz2/(p1278q-cmz2)
  
  
  
* compute the 2 diagrams with the w+w-_h_w+w- coupling                  
  
* rmh < 0 in our convention means no higgs coupling                     
      if (rmh.ge.0.d0) then
  
      cres2w_h_2w=(ch1234*ch5678)/(p1234q-cmh2)+
     &            (ch1278*ch5634)/(p1278q-cmh2)
  
      endif
  
* compute the quadruple coupling diagram  W+W-W+W-                      
* p.q -- p.q=c12dot34,p=cw12%e,q=cw34%e,bef=,aft=
      c12dot34=(cw12%e(0)*cw34%e(0)-cw12%e(1)*cw34%e(1)-cw12%e(2
     & )*cw34%e(2)-cw12%e(3)*cw34%e(3))
* p.q -- p.q=c12dot56,p=cw12%e,q=cw56%e,bef=,aft=
      c12dot56=(cw12%e(0)*cw56%e(0)-cw12%e(1)*cw56%e(1)-cw12%e(2
     & )*cw56%e(2)-cw12%e(3)*cw56%e(3))
* p.q -- p.q=c12dot78,p=cw12%e,q=cw78%e,bef=,aft=
      c12dot78=(cw12%e(0)*cw78%e(0)-cw12%e(1)*cw78%e(1)-cw12%e(2
     & )*cw78%e(2)-cw12%e(3)*cw78%e(3))
* p.q -- p.q=c34dot56,p=cw34%e,q=cw56%e,bef=,aft=
      c34dot56=(cw34%e(0)*cw56%e(0)-cw34%e(1)*cw56%e(1)-cw34%e(2
     & )*cw56%e(2)-cw34%e(3)*cw56%e(3))
* p.q -- p.q=c34dot78,p=cw34%e,q=cw78%e,bef=,aft=
      c34dot78=(cw34%e(0)*cw78%e(0)-cw34%e(1)*cw78%e(1)-cw34%e(2
     & )*cw78%e(2)-cw34%e(3)*cw78%e(3))
* p.q -- p.q=c56dot78,p=cw56%e,q=cw78%e,bef=,aft=
      c56dot78=(cw56%e(0)*cw78%e(0)-cw56%e(1)*cw78%e(1)-cw56%e(2
     & )*cw78%e(2)-cw56%e(3)*cw78%e(3))
  
  
      cres4wquad=2.d0*c12dot56*c34dot78-c12dot34*c56dot78-
     &     c12dot78*c34dot56
* The - sign to the coupling comes from the fact that we are neglecting 
* cim in every propagator and coupling. The diagram with quartic        
* coupling should have cim**9 and the others cim**11, hence we          
* give a - sign to it with respect to the others.                       
* We factorize out from the amplitude the modulus of the electric       
* charge from every coupling, hence instead of g**2  we have            
* coupling g**2/|e|**2=1/s2w                                            
  
      cres4wquad=-cres4wquad/s2w
  
  
      do m=1,7
        cres(m)=czero
      enddo
  
  
* sum the ew results                                                    
* cres(1):   1| 3| 5| 7|                                                
*             |  |  |  |   (pure electroweak contribution)              
*            2| 4| 6| 8|                                                
  
      cres(1) = cresUD_3w+cres2w_f_2w+cres2w_z_2w+cres2w_h_2w+cres4wquad
  
  
**QCD                                                                   
* compute the diagrams with a gluon propagator connecting 2 w's on both 
* sides:                                                                
*                                                                       
* e.g.:                    |     |                                      
*                      >-W-|     |-W-<                                  
*                          |~~g~~|                                      
*                          |     |                                      
*                                                                       
  
**QCD                                                                   
* compute separately the different colour flow configurations           
  
* cres(2):   5| 1| 3| 7|   (QCD contribution)                           
*             |  |~~|  |                                                
*            6| 2| 4| 8|   cg1278.cg3456                                
  
      if(ilept(id1).ne.1 .and. ilept(id3).ne.1)then
* p.q -- p.q=cres(2),p=cg1278,q=cg3456,bef=-,aft=/p1278q
      cres(2)=-(cg1278(0)*cg3456(0)-cg1278(1)*cg3456(1)-cg1278(2
     & )*cg3456(2)-cg1278(3)*cg3456(3))/p1278q
      endif
  
  
* cres(3):   3| 1| 5| 7|   (QCD contribution)                           
*             |  |~~|  |                                                
*            4| 2| 6| 8|   cg1234.cg5678 + cg1278.cg5634                
  
      if(ilept(id1).ne.1 .and. ilept(id5).ne.1)then
* p.q -- p.q=cres(3),p=cg1234,q=cg5678,bef=-,aft=/p1234q
      cres(3)=-(cg1234(0)*cg5678(0)-cg1234(1)*cg5678(1)-cg1234(2
     & )*cg5678(2)-cg1234(3)*cg5678(3))/p1234q
  
* p.q -- p.q=cres(3),p=cg1278,q=cg5634,bef=cres(3)-,aft=/p1278q
      cres(3)=cres(3)-(cg1278(0)*cg5634(0)-cg1278(1)*cg5634(1)-c
     & g1278(2)*cg5634(2)-cg1278(3)*cg5634(3))/p1278q
      endif
  
  
* cres(4):   1| 3| 5| 7|   (QCD contribution)                           
*             |  |~~|  |                                                
*            2| 4| 6| 8|   cg3412.cg5678                                
  
      if(ilept(id3).ne.1 .and. ilept(id5).ne.1)then
* p.q -- p.q=cres(4),p=cg3412,q=cg5678,bef=-,aft=/p1234q
      cres(4)=-(cg3412(0)*cg5678(0)-cg3412(1)*cg5678(1)-cg3412(2
     & )*cg5678(2)-cg3412(3)*cg5678(3))/p1234q
      endif
  
  
* cres(5):   3| 1| 7| 5|   (QCD contribution)                           
*             |  |~~|  |                                                
*            4| 2| 8| 6|   cg1234.cg7856                                
  
      if(ilept(id1).ne.1 .and. ilept(id7).ne.1)then
* p.q -- p.q=cres(5),p=cg1234,q=cg7856,bef=-,aft=/p1234q
      cres(5)=-(cg1234(0)*cg7856(0)-cg1234(1)*cg7856(1)-cg1234(2
     & )*cg7856(2)-cg1234(3)*cg7856(3))/p1234q
      endif
  
  
* cres(6):   1| 3| 7| 5|   (QCD contribution)                           
*             |  |~~|  |                                                
*            2| 4| 8| 6|   cg3412.cg7856 + cg3456.cg7812                
  
      if(ilept(id3).ne.1 .and. ilept(id7).ne.1)then
* p.q -- p.q=cres(6),p=cg3412,q=cg7856,bef=-,aft=/p1234q
      cres(6)=-(cg3412(0)*cg7856(0)-cg3412(1)*cg7856(1)-cg3412(2
     & )*cg7856(2)-cg3412(3)*cg7856(3))/p1234q
  
* p.q -- p.q=cres(6),p=cg3456,q=cg7812,bef=cres(6)-,aft=/p1278q
      cres(6)=cres(6)-(cg3456(0)*cg7812(0)-cg3456(1)*cg7812(1)-c
     & g3456(2)*cg7812(2)-cg3456(3)*cg7812(3))/p1278q
      endif
  
  
* cres(7):   1| 5| 7| 3|   (QCD contribution)                           
*             |  |~~|  |                                                
*            2| 6| 8| 4|   cg5634.cg7812                                
  
      if(ilept(id5).ne.1 .and. ilept(id7).ne.1)then
* p.q -- p.q=cres(7),p=cg5634,q=cg7812,bef=-,aft=/p1278q
      cres(7)=-(cg5634(0)*cg7812(0)-cg5634(1)*cg7812(1)-cg5634(2
     & )*cg7812(2)-cg5634(3)*cg7812(3))/p1278q
      endif
  
  
* the resulting amplitudes are now divided   by sqrt(abs(p.k0))         
* for every external particle                                           
  
      do i=1,7
        cres(i) = cres(i)/sqrt(p1k0*p2k0*p3k0*p4k0*p5k0*p6k0*p7k0*p8k0)
      enddo
  
  
* Sum and average over spin and colour factors                          
* has to be performed externally, when the |amp|**2 is computed         
  
  
      return
      end
  
  
  
  
