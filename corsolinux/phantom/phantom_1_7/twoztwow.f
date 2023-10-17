************************************************************************
*                                                                       
* This routine computes the basic amplitude corresponding to 8          
*  outgoing  fermions which can form 2 Z and 2 W-.                      
*  The 4 fermion which can for two Z are massive, the other massless.   
*  All outgoing particles are considered different.                     
*  The input particles are ordered from 1 to 8 in such a way that       
*  odd are particles, even antiparticles. 12 corresponds to Z,          
*  34 to Z, 56 to W+, 78 to W-                                          
*  p1, ...p8 are all outgoing momenta                                   
*  id1....id8 give the identities of the outgoing particles             
*  cres  is on output the complex helicity amplitude                    
*  cres=cres(2,2,2,2)                                                   
*   the  four indeces refer to the chiralities of the first four        
*   particles, while the others four have only                          
*   chirality indeces  = 2   (corresponding to -)                       
*                                                                       
************************************************************************
*  Comments:  "W" and "Z"  stay for, respectivecly, a couple of quark   
*   which can form a W  and a couple of quarks which can form a Z       
*   Correpsondingly Wline and Zline stay for a fermion line which       
*    ends with two quarks which can form a W or a Z                     
*                                                                       
  
      subroutine twoztwow(p1,p2,p3,p4,p5,p6,p7,p8,
     &              id1,id2,id3,id4,id5,id6,id7,id8,cres)
  
      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)
  
*four momenta                                                           
  
      dimension p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3),p7(0:3)
     &     ,p8(0:3)
  
      dimension p12(0:3),p34(0:3),p56(0:3),p78(0:3)
      dimension p578(0:3),p756(0:3)
     &     ,p512(0:3),p534(0:3),p712(0:3),p734(0:3)
     &     ,p134(0:3),p156(0:3),p178(0:3),p312(0:3),p356(0:3),p378(0:3)
      dimension p678(0:3),p856(0:3)
     &     ,p612(0:3),p634(0:3),p812(0:3),p834(0:3)
     &     ,p234(0:3),p256(0:3),p278(0:3),p412(0:3),p456(0:3),p478(0:3)
      dimension p1234(0:3),p1256(0:3),p1278(0:3)
  
*results                                                                
      dimension cres(2,2,2,2),cunit(2,2,2,2)
      dimension cresw_3fork(2,2,2,2),cresz_3fork(2,2,2,2)
     &     ,creswz_w_wz(2,2,2,2),creszz_z_ww(2,2,2,2)
     &     ,creszz_f_ww(2,2,2,2),creszz_h_ww(2,2,2,2)
     &     ,creszz_hh_ww(2,2,2,2)
     &     ,cresquad(2,2,2,2)
*scalar products for w propagators                                      
      dimension cp1256dotcw1256(2,2),cp1256dotcw3478(2,2)
     &     ,cp1278dotcw1278(2,2),cp1278dotcw3456(2,2)
* triple vertices auxiliary                                             
      dimension vfz(0:3), vwm(0:3), vwp(0:3)
  
* Z and gamma "decaying" to two w's                                     
      dimension cz5678(0:3),cf5678(0:3),ctrip5678(0:3)
  
* W "decaying" to one w and one z                                       
      dimension cw1256(2,2,0:3),cw1278(2,2,0:3),cw3456(2,2,0:3)
     &     ,cw3478(2,2,0:3),ctrip1256(2,2,0:3),ctrip1278(2,2,0:3)
     &     ,ctrip3456(2,2,0:3),ctrip3478(2,2,0:3)
* Z and gamma "decaying" to two z's                                     
      dimension cz1234(2,2,2,2,0:3),cf1234(2,2,2,2,0:3)
  
* higgs "decaying" to two z's                                           
      dimension ch1234(2,2,2,2)
* higgs "decaying" to two w's                                           
      double complex ch5678
* auxiliary                                                             
      dimension cdum(2,2,2,2)
* heavy higgs "decaying" to two z's                                     
      dimension chh1234(2,2,2,2)
* heavy higgs "decaying" to two w's                                     
      double complex chh5678
  
*zw                                                                     
      dimension iresw(8),iresz(8)
      dimension pw(0:3),pz(0:3)
      dimension ceps1(0:3),ceps2(0:3)
      dimension cscalz(2,2),cscal1z(2,2),cscal2z(2,2)
*zwend                                                                  
  
  
* auxiliary structure used in mline                                     
      type aux
        double complex a(2,2),b(2,2),c(2,2),d(2,2)
      end type
      type(aux) clineth,clinet_mu(0:3),clinetww_mu(0:3),
     &     clinetzz(2,2,0:3), clinetfz(2,2,0:3),clineth12(2,2)
     &     ,clineth34(2,2),clinet12_3fork(2,2),clinet34_3fork(2,2)
  
* forks                                                                 
      dimension ch12(2,2),ch34(2,2)
  
      type polcom
        double complex e(0:3),ek0,v
      end type
      type(polcom) cw56,cw78
     &     ,cz12(2,2),cz34(2,2),cf12(2,2),cf34(2,2)
  
* left t Wline                                                          
* only to be used with  w insertions, hence tw0 is enough               
      type l_wline
         double complex a(2:2),c(2:2)
      end type
      type(l_wline) l5_78,l7_56
     &     ,l5_12(2,2),l5_34(2,2),l7_12(2,2),l7_34(2,2)
     &     ,l5_1234(2,2,2,2),l5_1278(2,2),l5_3478(2,2)
     &     ,l7_1234(2,2,2,2),l7_1256(2,2),l7_3456(2,2)
* with z, gamma etc the second indices represent down four momentum     
*  and not incoming as for li_kl  et similar                            
     &     ,lw5_612(0:3),lw5_634(0:3),lw7_812(0:3),lw7_834(0:3)
     &     ,lz5_678(0:3),lz7_856(0:3),lf5_678(0:3),lf7_856(0:3)
  
* left t Zline                                                          
      type l_zline
         double complex a(2,2),b(2,2),c(2,2),d(2,2)
      end type
      type(l_zline) l1_34(2,2),l1_56,l1_78,l3_12(2,2),l3_56,l3_78
     &     ,l1_3456(2,2),l1_3478(2,2),l1_5678
     &     ,l3_1256(2,2),l3_1278(2,2),l3_5678
* with z, gamma etc the second indices represent down four momentum     
*  and not incoming as for li_kl  et similar                            
     &     ,lw1_256(0:3),lw1_278(0:3),lz1_234(0:3),lf1_234(0:3)
     &     ,lw3_456(0:3),lw3_478(0:3),lz3_412(0:3),lf3_412(0:3)
     &     ,lh1_234,lh3_412
* u t Wline                                                             
* only to be used with w insertions, hence tw0 is enough                
      type u_wline
         double complex a(2:2),b(1:1),c(2:2),d(1:1)
      end type
      type(u_wline) u512_34(2,2),u512_78,u534_12(2,2),u534_78
     &     ,u578_12(2,2),u578_34(2,2)
     &     ,u712_34(2,2),u712_56,u734_12(2,2),u734_56
     &     ,u756_12(2,2),u756_34(2,2)
* u t Zline                                                             
      type u_zline
         double complex a(2,2),b(2,2),c(2,2),d(2,2)
      end type
      type(u_zline) u134_56,u134_78,u156_34(2,2),u156_78
     &     ,u178_34(2,2),u178_56
     &     ,u312_56,u312_78,u356_12(2,2),u356_78
     &     ,u378_12(2,2),u378_56
  
* right t Wline                                                         
* only to be used with w insertions, hence tw0 is enough                
      type r_wline
         double complex a(2:2),b(1:1)
      end type
      type(r_wline) r6_12(2,2),r6_34(2,2),r6_78,r8_12(2,2),r8_34(2,2)
     &     ,r8_56
* with z, gamma etc the second indices represent down four momentum     
*  and not incoming as for ri_kl  et similar                            
     &     ,rw6_512(0:3),rw6_534(0:3),rw8_712(0:3),rw8_734(0:3)
     &     ,rz6_578(0:3),rz8_756(0:3),rf6_578(0:3),rf8_756(0:3)
  
* right t Zline                                                         
      type r_zline
         double complex a(2,2),b(2,2),c(2,2),d(2,2)
      end type
      type(r_zline) r2_34(2,2),r2_56,r2_78,r4_12(2,2),r4_56,r4_78
* with z, gamma etc the second indices represent down four momentum     
*  and not incoming as for ri_kl  et similar                            
     &     ,rw2_156(0:3),rw2_178(0:3),rz2_134(0:3),rf2_134(0:3)
     &     ,rw4_356(0:3),rw4_378(0:3),rz4_312(0:3),rf4_312(0:3)
     &     ,rh2_134,rh4_312
  
*common blocks                                                          
* generic common                                                        
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
      endif
  
      if (i_ww.gt.0) then
          iresw(5)=0
          if (p5(0).gt.0.d0) then
            if (id5.eq.idw(1).or.(i_ww.eq.2.and.id5.eq.idw(3)))
     &           iresw(5)=1
          endif
          iresw(7)=0
          if (p7(0).gt.0.d0) then
            if (id7.eq.idw(1).or.(i_ww.eq.2.and.id7.eq.idw(3)))
     &           iresw(7)=1
          endif
      endif
*zwend                                                                  
  
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
  
      do mu=0,3
        p1256(mu)=p12(mu)+p56(mu)
      enddo
* p.q -- p.q=p1256q,p=p1256,q=p1256,bef=,aft=
      p1256q=(p1256(0)*p1256(0)-p1256(1)*p1256(1)-p1256(2)*p1256
     & (2)-p1256(3)*p1256(3))
      do mu=0,3
        p1278(mu)=p12(mu)+p78(mu)
      enddo
* p.q -- p.q=p1278q,p=p1278,q=p1278,bef=,aft=
      p1278q=(p1278(0)*p1278(0)-p1278(1)*p1278(1)-p1278(2)*p1278
     & (2)-p1278(3)*p1278(3))
  
  
* compute all forks ( W's, Z's, gamma's, higgs's                        
*   "decaying" to two fermions)                                         
*                                                                       
  
*      W's                                                              
* quqd -- p=p5,q=p6
      quqd=p5(0)*p6(0)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3)
      s56=2.d0*quqd
*zw                                                                     
      cauxmw2=cmw2
      if (i_ww.eq.2.and.i_osp.eq.2.and.iresw(5).eq.1) then
        cauxmw2=rmw2*cuno-cim*gamw*rmw
      endif
      cdw=-wcl/(s56-cauxmw2)
*zwend                                                                  
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
  
*zw                                                                     
      iceps1=0
      iceps2=0
      if (iresw(5).eq.1) then
*  compute the W momentum                                               
        do mu=0,3
          pw(mu)=p5(mu)+p6(mu)
        enddo
        if (id5.eq.idw(1)) then
          jpol=ipolw(1)
        elseif (id5.eq.idw(3)) then
          jpol=ipolw(2)
        else
          print *, 'ERROR, MISMATCH IN POLARIZATIONS'
          stop
        endif
  
        if (jpol.ne.0.and.jpol.ne.4) then
          call pol(pw,jpol,ceps1)
          iceps1=1
        endif
        if (jpol.eq.4) then
* pol w transverse                                                      
          ip=2
          call pol(pw,ip,ceps1)
          ip=3
          call pol(pw,ip,ceps2)
          iceps2=1
        endif
  
        if (iceps1.eq.1) then
* p.q -- p.q=cscal,p=ceps1,q=cw56%e,bef=,aft=
      cscal=(ceps1(0)*cw56%e(0)-ceps1(1)*cw56%e(1)-ceps1(2)*cw56
     & %e(2)-ceps1(3)*cw56%e(3))
  
          do mu=0,3
            cw56%e(mu)=cscal*conjg(ceps1(mu))
          enddo
        endif
  
        if (iceps2.eq.1) then
* p.q -- p.q=cscal1,p=ceps1,q=cw56%e,bef=,aft=
      cscal1=(ceps1(0)*cw56%e(0)-ceps1(1)*cw56%e(1)-ceps1(2)*cw5
     & 6%e(2)-ceps1(3)*cw56%e(3))
  
* p.q -- p.q=cscal2,p=ceps2,q=cw56%e,bef=,aft=
      cscal2=(ceps2(0)*cw56%e(0)-ceps2(1)*cw56%e(1)-ceps2(2)*cw5
     & 6%e(2)-ceps2(3)*cw56%e(3))
  
          do mu=0,3
            cw56%e(mu)=cscal1*conjg(ceps1(mu))
     &           +cscal2*conjg(ceps2(mu))
          enddo
        endif
      endif
  
*zwend                                                                  
  
* pk0 -- p=cw56%e
      cw56%ek0=cw56%e(0)-cw56%e(1)
* quqd -- p=p7,q=p8
      quqd=p7(0)*p8(0)-p7(1)*p8(1)-p7(2)*p8(2)-p7(3)*p8(3)
      s78=2.d0*quqd
*zw                                                                     
      cauxmw2=cmw2
      if (i_ww.eq.2.and.i_osp.eq.2.and.iresw(7).eq.1) then
        cauxmw2=rmw2*cuno-cim*gamw*rmw
      endif
      cdw=-wcl/(s78-cauxmw2)
*zwend                                                                  
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
  
*zw                                                                     
      iceps1=0
      iceps2=0
      if (iresw(7).eq.1) then
*  compute the W momentum                                               
        do mu=0,3
          pw(mu)=p7(mu)+p8(mu)
        enddo
        if (id7.eq.idw(1)) then
          jpol=ipolw(1)
        elseif (id7.eq.idw(3)) then
          jpol=ipolw(2)
        else
          print *, 'ERROR, MISMATCH IN POLARIZATIONS'
          stop
        endif
  
        if (jpol.ne.0.and.jpol.ne.4) then
          call pol(pw,jpol,ceps1)
          iceps1=1
        endif
        if (jpol.eq.4) then
* pol w transverse                                                      
          ip=2
          call pol(pw,ip,ceps1)
          ip=3
          call pol(pw,ip,ceps2)
          iceps2=1
        endif
  
        if (iceps1.eq.1) then
* p.q -- p.q=cscal,p=ceps1,q=cw78%e,bef=,aft=
      cscal=(ceps1(0)*cw78%e(0)-ceps1(1)*cw78%e(1)-ceps1(2)*cw78
     & %e(2)-ceps1(3)*cw78%e(3))
  
          do mu=0,3
            cw78%e(mu)=cscal*conjg(ceps1(mu))
          enddo
        endif
  
        if (iceps2.eq.1) then
* p.q -- p.q=cscal1,p=ceps1,q=cw78%e,bef=,aft=
      cscal1=(ceps1(0)*cw78%e(0)-ceps1(1)*cw78%e(1)-ceps1(2)*cw7
     & 8%e(2)-ceps1(3)*cw78%e(3))
  
* p.q -- p.q=cscal2,p=ceps2,q=cw78%e,bef=,aft=
      cscal2=(ceps2(0)*cw78%e(0)-ceps2(1)*cw78%e(1)-ceps2(2)*cw7
     & 8%e(2)-ceps2(3)*cw78%e(3))
  
          do mu=0,3
            cw78%e(mu)=cscal1*conjg(ceps1(mu))
     &           +cscal2*conjg(ceps2(mu))
          enddo
        endif
      endif
  
*zwend                                                                  
  
* pk0 -- p=cw78%e
      cw78%ek0=cw78%e(0)-cw78%e(1)
  
*      Z's and gammas  and higgs                                        
* quqd -- p=p1,q=p2
      quqd=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      s12=2.d0*quqd+rmass2(id1)+rmass2(id2)
  
*zw                                                                     
  
  
*zw1                                                                    
  
  
      cauxmz2=cmz2
      if ((i_zz.eq.2.or.(i_ww.eq.1.and.i_zz.eq.1)).
     &     and.i_osp.eq.2.and.iresz(1).eq.1) then
        cauxmz2=rmz2*cuno-cim*gamz*rmz
      endif
  
  
      ccr=zcr(id1)/(-s12+cmauxz2)
      ccl=zcl(id1)/(-s12+cmauxz2)
* T -- qu=p1,qd=p2,v=0,a=clinet_mu(0)%a,b=clinet_mu(0)%b,c=clinet_mu(0)%
* c,d=clinet_mu(0)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p1(2)*p2(3)+p2(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p1k0*p2(0)+p2k0*p1(0)
      clinet_mu(0)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(0)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(0)%b(1,2)=-ccl*(p2(2)+ceps_2)
      clinet_mu(0)%b(2,1)=ccr*(p2(2)-ceps_2)
      clinet_mu(0)%c(1,2)=ccr*(p1(2)+ceps_1)
      clinet_mu(0)%c(2,1)=ccl*(-p1(2)+ceps_1)
      clinet_mu(0)%d(1,1)=ccl
      clinet_mu(0)%d(2,2)=ccr
* T -- qu=p1,qd=p2,v=1,a=clinet_mu(1)%a,b=clinet_mu(1)%b,c=clinet_mu(1)%
* c,d=clinet_mu(1)%d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p1k0*p2(1)+p2k0*p1(1)
      clinet_mu(1)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(1)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(1)%b(1,2)=-ccl*(p2(2)+ceps_2)
      clinet_mu(1)%b(2,1)=ccr*(p2(2)-ceps_2)
      clinet_mu(1)%c(1,2)=ccr*(p1(2)+ceps_1)
      clinet_mu(1)%c(2,1)=ccl*(-p1(2)+ceps_1)
      clinet_mu(1)%d(1,1)=ccl
      clinet_mu(1)%d(2,2)=ccr
* T -- qu=p1,qd=p2,v=2,a=clinet_mu(2)%a,b=clinet_mu(2)%b,c=clinet_mu(2)%
* c,d=clinet_mu(2)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p1k0*p2(3)+p2k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(2)+p2k0*p1(2)
      clinet_mu(2)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(2)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(2)%b(1,2)=-ccl*p2k0
      clinet_mu(2)%b(2,1)=ccr*p2k0
      clinet_mu(2)%c(1,2)=ccr*p1k0
      clinet_mu(2)%c(2,1)=-ccl*p1k0
* T -- qu=p1,qd=p2,v=3,a=clinet_mu(3)%a,b=clinet_mu(3)%b,c=clinet_mu(3)%
* c,d=clinet_mu(3)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=p1k0*p2(2)-p2k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      ceps_2=p2k0*cim
      auxa=+p1k0*p2(3)+p2k0*p1(3)
      clinet_mu(3)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(3)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(3)%b(1,2)=-ccl*ceps_2
      clinet_mu(3)%b(2,1)=-ccr*ceps_2
      clinet_mu(3)%c(1,2)=ccr*ceps_1
      clinet_mu(3)%c(2,1)=ccl*ceps_1
  
*zw1end                                                                 
  
      rmassl=rmass(id1)
      rmassr=-rmassl
  
      do mu=0,3
* mline -- res=cz12(&1,&2)%e(mu),abcd=clinet_mu(mu)%,m1=rmassl,m2=rmassr
* ,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cz12(iut,jut)%e(mu)=clinet_mu(mu)%a(iut,jut)+rmassl*clinet
     & _mu(mu)%b(iut,jut)+rmassr*clinet_mu(mu)%c(iut,jut)+rmassl
     & *rmassr*clinet_mu(mu)%d(iut,jut)
      enddo
      enddo
      end do
  
      if (imass(id1).eq.1) then
        do i1=1,2
          do i2=1,2
* p.q -- p.q=cauxdot,p=p12,q=cz12(i1,i2)%e,bef=,aft=
      cauxdot=(p12(0)*cz12(i1,i2)%e(0)-p12(1)*cz12(i1,i2)%e(1)-p
     & 12(2)*cz12(i1,i2)%e(2)-p12(3)*cz12(i1,i2)%e(3))
            do mu=0,3
             cz12(i1,i2)%e(mu)=cz12(i1,i2)%e(mu)-p12(mu)*cauxdot
     &             /cmz2
            enddo
          enddo
        enddo
      endif
*zw1                                                                    
  
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
            doj=1,2
* p.q -- p.q=cscalz(i,j),p=ceps1,q=cz12(i,j)%e,bef=,aft=
      cscalz(i,j)=(ceps1(0)*cz12(i,j)%e(0)-ceps1(1)*cz12(i,j)%e(
     & 1)-ceps1(2)*cz12(i,j)%e(2)-ceps1(3)*cz12(i,j)%e(3))
  
              do mu=0,3
                cz12(i,j)%e(mu)=cscalz(i,j)*conjg(ceps1(mu))
              enddo
            enddo
          enddo
        endif
  
        if (iceps2.eq.1) then
          do i=1,2
            doj=1,2
* p.q -- p.q=cscal1z(i,j),p=ceps1,q=cz12(i,j)%e,bef=,aft=
      cscal1z(i,j)=(ceps1(0)*cz12(i,j)%e(0)-ceps1(1)*cz12(i,j)%e
     & (1)-ceps1(2)*cz12(i,j)%e(2)-ceps1(3)*cz12(i,j)%e(3))
  
* p.q -- p.q=cscal2z(i,j),p=ceps2,q=cz12(i,j)%e,bef=,aft=
      cscal2z(i,j)=(ceps2(0)*cz12(i,j)%e(0)-ceps2(1)*cz12(i,j)%e
     & (1)-ceps2(2)*cz12(i,j)%e(2)-ceps2(3)*cz12(i,j)%e(3))
  
              do mu=0,3
              cz12(i,j)%e(mu)=cscal1z(i,j)*conjg(ceps1(mu))
     &           +cscal2z(i,j)*conjg(ceps2(mu))
              enddo
            enddo
          enddo
        endif
      endif
  
*zw1end                                                                 
  
      do i1=1,2
      do i2=1,2
* pk0 -- p=cz12(i1,i2)%e
      cz12(i1,i2)%ek0=cz12(i1,i2)%e(0)-cz12(i1,i2)%e(1)
      end do
      end do
  
*zw1                                                                    
  
      if(iresz(1).eq.1) then
        do i=1,2
          do j=1,2
            do mu=0,3
              cf12(i,j)%e(mu)=czero
            enddo
          enddo
        enddo
      else
  
      ccr=fcr(id1)/(-s12)
      ccl=fcl(id1)/(-s12)
* T -- qu=p1,qd=p2,v=0,a=clinet_mu(0)%a,b=clinet_mu(0)%b,c=clinet_mu(0)%
* c,d=clinet_mu(0)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p1(2)*p2(3)+p2(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p1k0*p2(0)+p2k0*p1(0)
      clinet_mu(0)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(0)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(0)%b(1,2)=-ccl*(p2(2)+ceps_2)
      clinet_mu(0)%b(2,1)=ccr*(p2(2)-ceps_2)
      clinet_mu(0)%c(1,2)=ccr*(p1(2)+ceps_1)
      clinet_mu(0)%c(2,1)=ccl*(-p1(2)+ceps_1)
      clinet_mu(0)%d(1,1)=ccl
      clinet_mu(0)%d(2,2)=ccr
* T -- qu=p1,qd=p2,v=1,a=clinet_mu(1)%a,b=clinet_mu(1)%b,c=clinet_mu(1)%
* c,d=clinet_mu(1)%d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p1k0*p2(1)+p2k0*p1(1)
      clinet_mu(1)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(1)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(1)%b(1,2)=-ccl*(p2(2)+ceps_2)
      clinet_mu(1)%b(2,1)=ccr*(p2(2)-ceps_2)
      clinet_mu(1)%c(1,2)=ccr*(p1(2)+ceps_1)
      clinet_mu(1)%c(2,1)=ccl*(-p1(2)+ceps_1)
      clinet_mu(1)%d(1,1)=ccl
      clinet_mu(1)%d(2,2)=ccr
* T -- qu=p1,qd=p2,v=2,a=clinet_mu(2)%a,b=clinet_mu(2)%b,c=clinet_mu(2)%
* c,d=clinet_mu(2)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p1k0*p2(3)+p2k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p2(2)+p2k0*p1(2)
      clinet_mu(2)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(2)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(2)%b(1,2)=-ccl*p2k0
      clinet_mu(2)%b(2,1)=ccr*p2k0
      clinet_mu(2)%c(1,2)=ccr*p1k0
      clinet_mu(2)%c(2,1)=-ccl*p1k0
* T -- qu=p1,qd=p2,v=3,a=clinet_mu(3)%a,b=clinet_mu(3)%b,c=clinet_mu(3)%
* c,d=clinet_mu(3)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=p1k0*p2(2)-p2k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      ceps_2=p2k0*cim
      auxa=+p1k0*p2(3)+p2k0*p1(3)
      clinet_mu(3)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(3)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(3)%b(1,2)=-ccl*ceps_2
      clinet_mu(3)%b(2,1)=-ccr*ceps_2
      clinet_mu(3)%c(1,2)=ccr*ceps_1
      clinet_mu(3)%c(2,1)=ccl*ceps_1
      rmassl=rmass(id1)
      rmassr=-rmassl
      do mu=0,3
* mline -- res=cf12(&1,&2)%e(mu),abcd=clinet_mu(mu)%,m1=rmassl,m2=rmassr
* ,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cf12(iut,jut)%e(mu)=clinet_mu(mu)%a(iut,jut)+rmassl*clinet
     & _mu(mu)%b(iut,jut)+rmassr*clinet_mu(mu)%c(iut,jut)+rmassl
     & *rmassr*clinet_mu(mu)%d(iut,jut)
      enddo
      enddo
      end do
  
      do i1=1,2
      do i2=1,2
* pk0 -- p=cf12(i1,i2)%e
      cf12(i1,i2)%ek0=cf12(i1,i2)%e(0)-cf12(i1,i2)%e(1)
      end do
      end do
  
      endif
  
*zw1end                                                                 
  
*zwend                                                                  
  
*      higgs                                                            
* rmh < 0 in our convention means no higgs coupling                     
  
      if (id1.eq.5.and.rmh.ge.0.d0) then
  
  
*zw                                                                     
  
        if(iresz(1).eq.1) then
          do i=1,2
            do j=1,2
                ch12(i,j)=czero
            enddo
          enddo
        else
  
* TH -- qu=p1,qd=p2,a=clineth%a,b=clineth%b,c=clineth%c,coupl=1.d0
      auxa=-p1k0*p2(2)+p2k0*p1(2)
      cauxa=auxa-cim*(p2(3)*p1k0-p1(3)*p2k0)
      clineth%a(1,2)=cauxa
      clineth%a(2,1)=-conjg(cauxa)
      clineth%b(1,1)=p2k0
      clineth%b(2,2)=clineth%b(1,1)
      clineth%c(1,1)=p1k0
      clineth%c(2,2)=clineth%c(1,1)
* mline -- res=ch12(&1,&2),abcd=clineth%,m1=rmassl,m2=rmassr,den=0,nsum=
* 0
      do iut=1,2
      do jut=1,2
      ch12(iut,jut)=clineth%a(iut,jut)+rmassl*clineth%b(iut,jut)
     & +rmassr*clineth%c(iut,jut)+rmassl*rmassr*clineth%d(iut,ju
     & t)
      enddo
      enddo
        endif
*zwend                                                                  
  
      else
        do i1=1,2
          do i2=1,2
            ch12(i1,i2)=czero
          enddo
        enddo
      endif
* quqd -- p=p3,q=p4
      quqd=p3(0)*p4(0)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3)
      s34=2.d0*quqd+rmass2(id3)+rmass2(id4)
  
*zw                                                                     
  
  
*zw1                                                                    
  
  
      cauxmz2=cmz2
      if ((i_zz.eq.2.or.(i_ww.eq.1.and.i_zz.eq.1)).
     &     and.i_osp.eq.2.and.iresz(3).eq.1) then
        cauxmz2=rmz2*cuno-cim*gamz*rmz
      endif
  
  
      ccr=zcr(id3)/(-s34+cmauxz2)
      ccl=zcl(id3)/(-s34+cmauxz2)
* T -- qu=p3,qd=p4,v=0,a=clinet_mu(0)%a,b=clinet_mu(0)%b,c=clinet_mu(0)%
* c,d=clinet_mu(0)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p4(3)+p4(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p3k0*p4(0)+p4k0*p3(0)
      clinet_mu(0)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(0)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(0)%b(1,2)=-ccl*(p4(2)+ceps_2)
      clinet_mu(0)%b(2,1)=ccr*(p4(2)-ceps_2)
      clinet_mu(0)%c(1,2)=ccr*(p3(2)+ceps_1)
      clinet_mu(0)%c(2,1)=ccl*(-p3(2)+ceps_1)
      clinet_mu(0)%d(1,1)=ccl
      clinet_mu(0)%d(2,2)=ccr
* T -- qu=p3,qd=p4,v=1,a=clinet_mu(1)%a,b=clinet_mu(1)%b,c=clinet_mu(1)%
* c,d=clinet_mu(1)%d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p4(1)+p4k0*p3(1)
      clinet_mu(1)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(1)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(1)%b(1,2)=-ccl*(p4(2)+ceps_2)
      clinet_mu(1)%b(2,1)=ccr*(p4(2)-ceps_2)
      clinet_mu(1)%c(1,2)=ccr*(p3(2)+ceps_1)
      clinet_mu(1)%c(2,1)=ccl*(-p3(2)+ceps_1)
      clinet_mu(1)%d(1,1)=ccl
      clinet_mu(1)%d(2,2)=ccr
* T -- qu=p3,qd=p4,v=2,a=clinet_mu(2)%a,b=clinet_mu(2)%b,c=clinet_mu(2)%
* c,d=clinet_mu(2)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p4(3)+p4k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(2)+p4k0*p3(2)
      clinet_mu(2)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(2)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(2)%b(1,2)=-ccl*p4k0
      clinet_mu(2)%b(2,1)=ccr*p4k0
      clinet_mu(2)%c(1,2)=ccr*p3k0
      clinet_mu(2)%c(2,1)=-ccl*p3k0
* T -- qu=p3,qd=p4,v=3,a=clinet_mu(3)%a,b=clinet_mu(3)%b,c=clinet_mu(3)%
* c,d=clinet_mu(3)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p4(2)-p4k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      ceps_2=p4k0*cim
      auxa=+p3k0*p4(3)+p4k0*p3(3)
      clinet_mu(3)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(3)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(3)%b(1,2)=-ccl*ceps_2
      clinet_mu(3)%b(2,1)=-ccr*ceps_2
      clinet_mu(3)%c(1,2)=ccr*ceps_1
      clinet_mu(3)%c(2,1)=ccl*ceps_1
  
*zw1end                                                                 
  
      rmassl=rmass(id3)
      rmassr=-rmassl
  
      do mu=0,3
* mline -- res=cz34(&1,&2)%e(mu),abcd=clinet_mu(mu)%,m1=rmassl,m2=rmassr
* ,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cz34(iut,jut)%e(mu)=clinet_mu(mu)%a(iut,jut)+rmassl*clinet
     & _mu(mu)%b(iut,jut)+rmassr*clinet_mu(mu)%c(iut,jut)+rmassl
     & *rmassr*clinet_mu(mu)%d(iut,jut)
      enddo
      enddo
      end do
  
      if (imass(id3).eq.1) then
        do i1=1,2
          do i2=1,2
* p.q -- p.q=cauxdot,p=p34,q=cz34(i1,i2)%e,bef=,aft=
      cauxdot=(p34(0)*cz34(i1,i2)%e(0)-p34(1)*cz34(i1,i2)%e(1)-p
     & 34(2)*cz34(i1,i2)%e(2)-p34(3)*cz34(i1,i2)%e(3))
            do mu=0,3
             cz34(i1,i2)%e(mu)=cz34(i1,i2)%e(mu)-p34(mu)*cauxdot
     &             /cmz2
            enddo
          enddo
        enddo
      endif
*zw1                                                                    
  
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
            doj=1,2
* p.q -- p.q=cscalz(i,j),p=ceps1,q=cz34(i,j)%e,bef=,aft=
      cscalz(i,j)=(ceps1(0)*cz34(i,j)%e(0)-ceps1(1)*cz34(i,j)%e(
     & 1)-ceps1(2)*cz34(i,j)%e(2)-ceps1(3)*cz34(i,j)%e(3))
  
              do mu=0,3
                cz34(i,j)%e(mu)=cscalz(i,j)*conjg(ceps1(mu))
              enddo
            enddo
          enddo
        endif
  
        if (iceps2.eq.1) then
          do i=1,2
            doj=1,2
* p.q -- p.q=cscal1z(i,j),p=ceps1,q=cz34(i,j)%e,bef=,aft=
      cscal1z(i,j)=(ceps1(0)*cz34(i,j)%e(0)-ceps1(1)*cz34(i,j)%e
     & (1)-ceps1(2)*cz34(i,j)%e(2)-ceps1(3)*cz34(i,j)%e(3))
  
* p.q -- p.q=cscal2z(i,j),p=ceps2,q=cz34(i,j)%e,bef=,aft=
      cscal2z(i,j)=(ceps2(0)*cz34(i,j)%e(0)-ceps2(1)*cz34(i,j)%e
     & (1)-ceps2(2)*cz34(i,j)%e(2)-ceps2(3)*cz34(i,j)%e(3))
  
              do mu=0,3
              cz34(i,j)%e(mu)=cscal1z(i,j)*conjg(ceps1(mu))
     &           +cscal2z(i,j)*conjg(ceps2(mu))
              enddo
            enddo
          enddo
        endif
      endif
  
*zw1end                                                                 
  
      do i1=1,2
      do i2=1,2
* pk0 -- p=cz34(i1,i2)%e
      cz34(i1,i2)%ek0=cz34(i1,i2)%e(0)-cz34(i1,i2)%e(1)
      end do
      end do
  
*zw1                                                                    
  
      if(iresz(3).eq.1) then
        do i=1,2
          do j=1,2
            do mu=0,3
              cf34(i,j)%e(mu)=czero
            enddo
          enddo
        enddo
      else
  
      ccr=fcr(id3)/(-s34)
      ccl=fcl(id3)/(-s34)
* T -- qu=p3,qd=p4,v=0,a=clinet_mu(0)%a,b=clinet_mu(0)%b,c=clinet_mu(0)%
* c,d=clinet_mu(0)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3(2)*p4(3)+p4(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p3k0*p4(0)+p4k0*p3(0)
      clinet_mu(0)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(0)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(0)%b(1,2)=-ccl*(p4(2)+ceps_2)
      clinet_mu(0)%b(2,1)=ccr*(p4(2)-ceps_2)
      clinet_mu(0)%c(1,2)=ccr*(p3(2)+ceps_1)
      clinet_mu(0)%c(2,1)=ccl*(-p3(2)+ceps_1)
      clinet_mu(0)%d(1,1)=ccl
      clinet_mu(0)%d(2,2)=ccr
* T -- qu=p3,qd=p4,v=1,a=clinet_mu(1)%a,b=clinet_mu(1)%b,c=clinet_mu(1)%
* c,d=clinet_mu(1)%d,cr=ccr,cl=ccl,nsum=0
      auxa=-quqd+p3k0*p4(1)+p4k0*p3(1)
      clinet_mu(1)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(1)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(1)%b(1,2)=-ccl*(p4(2)+ceps_2)
      clinet_mu(1)%b(2,1)=ccr*(p4(2)-ceps_2)
      clinet_mu(1)%c(1,2)=ccr*(p3(2)+ceps_1)
      clinet_mu(1)%c(2,1)=ccl*(-p3(2)+ceps_1)
      clinet_mu(1)%d(1,1)=ccl
      clinet_mu(1)%d(2,2)=ccr
* T -- qu=p3,qd=p4,v=2,a=clinet_mu(2)%a,b=clinet_mu(2)%b,c=clinet_mu(2)%
* c,d=clinet_mu(2)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=-p3k0*p4(3)+p4k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(2)+p4k0*p3(2)
      clinet_mu(2)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(2)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(2)%b(1,2)=-ccl*p4k0
      clinet_mu(2)%b(2,1)=ccr*p4k0
      clinet_mu(2)%c(1,2)=ccr*p3k0
      clinet_mu(2)%c(2,1)=-ccl*p3k0
* T -- qu=p3,qd=p4,v=3,a=clinet_mu(3)%a,b=clinet_mu(3)%b,c=clinet_mu(3)%
* c,d=clinet_mu(3)%d,cr=ccr,cl=ccl,nsum=0
      eps_0=p3k0*p4(2)-p4k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      ceps_2=p4k0*cim
      auxa=+p3k0*p4(3)+p4k0*p3(3)
      clinet_mu(3)%a(1,1)=ccr*(auxa+ceps_0)
      clinet_mu(3)%a(2,2)=ccl*(auxa-ceps_0)
      clinet_mu(3)%b(1,2)=-ccl*ceps_2
      clinet_mu(3)%b(2,1)=-ccr*ceps_2
      clinet_mu(3)%c(1,2)=ccr*ceps_1
      clinet_mu(3)%c(2,1)=ccl*ceps_1
      rmassl=rmass(id3)
      rmassr=-rmassl
      do mu=0,3
* mline -- res=cf34(&1,&2)%e(mu),abcd=clinet_mu(mu)%,m1=rmassl,m2=rmassr
* ,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cf34(iut,jut)%e(mu)=clinet_mu(mu)%a(iut,jut)+rmassl*clinet
     & _mu(mu)%b(iut,jut)+rmassr*clinet_mu(mu)%c(iut,jut)+rmassl
     & *rmassr*clinet_mu(mu)%d(iut,jut)
      enddo
      enddo
      end do
  
      do i1=1,2
      do i2=1,2
* pk0 -- p=cf34(i1,i2)%e
      cf34(i1,i2)%ek0=cf34(i1,i2)%e(0)-cf34(i1,i2)%e(1)
      end do
      end do
  
      endif
  
*zw1end                                                                 
  
*zwend                                                                  
  
*      higgs                                                            
* rmh < 0 in our convention means no higgs coupling                     
  
      if (id3.eq.5.and.rmh.ge.0.d0) then
  
  
*zw                                                                     
  
        if(iresz(3).eq.1) then
          do i=1,2
            do j=1,2
                ch34(i,j)=czero
            enddo
          enddo
        else
  
* TH -- qu=p3,qd=p4,a=clineth%a,b=clineth%b,c=clineth%c,coupl=1.d0
      auxa=-p3k0*p4(2)+p4k0*p3(2)
      cauxa=auxa-cim*(p4(3)*p3k0-p3(3)*p4k0)
      clineth%a(1,2)=cauxa
      clineth%a(2,1)=-conjg(cauxa)
      clineth%b(1,1)=p4k0
      clineth%b(2,2)=clineth%b(1,1)
      clineth%c(1,1)=p3k0
      clineth%c(2,2)=clineth%c(1,1)
* mline -- res=ch34(&1,&2),abcd=clineth%,m1=rmassl,m2=rmassr,den=0,nsum=
* 0
      do iut=1,2
      do jut=1,2
      ch34(iut,jut)=clineth%a(iut,jut)+rmassl*clineth%b(iut,jut)
     & +rmassr*clineth%c(iut,jut)+rmassl*rmassr*clineth%d(iut,ju
     & t)
      enddo
      enddo
        endif
*zwend                                                                  
  
      else
        do i3=1,2
          do i4=1,2
            ch34(i3,i4)=czero
          enddo
        enddo
      endif
  
  
  
* compute all single insertions of the type li_ (i=5,7)  for a Wline    
* and a W insertion                                                     
*        i __                                                           
*            |_W__/                                                     
*            |    \                                                     
*                                                                       
*together with its propagator and pk0                                   
  
      do m=0,3
        p578(m)=p5(m)+p78(m)
      enddo
* pk0 -- p=p578
      p578k0=p578(0)-p578(1)
* p.q -- p.q=p578q,p=p578,q=p578,bef=,aft=
      p578q=(p578(0)*p578(0)-p578(1)*p578(1)-p578(2)*p578(2)-p57
     & 8(3)*p578(3))
  
*zw                                                                     
  
      if(iresw(5).eq.1) then
       l5_78%a(2)=czero
       l5_78%c(2)=czero
      else
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
      endif
*zwend                                                                  
      do m=0,3
        p756(m)=p7(m)+p56(m)
      enddo
* pk0 -- p=p756
      p756k0=p756(0)-p756(1)
* p.q -- p.q=p756q,p=p756,q=p756,bef=,aft=
      p756q=(p756(0)*p756(0)-p756(1)*p756(1)-p756(2)*p756(2)-p75
     & 6(3)*p756(3))
  
*zw                                                                     
  
      if(iresw(7).eq.1) then
       l7_56%a(2)=czero
       l7_56%c(2)=czero
      else
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
      endif
*zwend                                                                  
  
* compute all single insertions of the type li_ (i=5,7)  for a Wline    
* and a Z and gamma insertion                                           
*        i __                                                           
*            |_Z,f_/                                                    
*            |     \                                                    
*                                                                       
*together with its propagator and pk0                                   
  
      do m=0,3
        p512(m)=p5(m)+p12(m)
      enddo
* pk0 -- p=p512
      p512k0=p512(0)-p512(1)
* p.q -- p.q=p512q,p=p512,q=p512,bef=,aft=
      p512q=(p512(0)*p512(0)-p512(1)*p512(1)-p512(2)*p512(2)-p51
     & 2(3)*p512(3))
*zw                                                                     
      if(iresw(5).eq.1) then
        do i1=1,2
          do i2=1,2
            l5_12(i1,i2)%a(2)=czero
            l5_12(i1,i2)%c(2)=czero
          enddo
        enddo
      else
*zwend                                                                  
* quqd -- p=p5,q=p512
      quqd=p5(0)*p512(0)-p5(1)*p512(1)-p5(2)*p512(2)-p5(3)*p512(
     & 3)
      ccl=zcl(id5)/(p512q*p512k0)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p512,v=cz12(i1,i2)%e,a=l5_12(i1,i2)%a,c=l5_12(i1,i2)%
* c,cl=ccl,nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p5(2)*p512(3)-p512(2)*p5(3))+p5k0
     & *(cz12(i1,i2)%e(2)*p512(3)-p512(2)*cz12(i1,i2)%e(3))-p512
     & k0*(cz12(i1,i2)%e(2)*p5(3)-p5(2)*cz12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p5k0+p5(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz12(i1,i2)%e(0)*p5(0)-cz12(i1,i2)%e(1)*p5(1)-cz12(i1
     & ,i2)%e(2)*p5(2)-cz12(i1,i2)%e(3)*p5(3)
      cvqd=cz12(i1,i2)%e(0)*p512(0)-cz12(i1,i2)%e(1)*p512(1)-cz1
     & 2(i1,i2)%e(2)*p512(2)-cz12(i1,i2)%e(3)*p512(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p5k0*cvqd+p512k0*cvqu
      cauxc=+cz12(i1,i2)%ek0*p5(2)-p5k0*cz12(i1,i2)%e(2)
      l5_12(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_12(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      ccl=fcl(id5)/(p512q*p512k0)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p512,v=cf12(i1,i2)%e,a=l5_12(i1,i2)%a,c=l5_12(i1,i2)%
* c,cl=ccl,nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p5(2)*p512(3)-p512(2)*p5(3))+p5k0
     & *(cf12(i1,i2)%e(2)*p512(3)-p512(2)*cf12(i1,i2)%e(3))-p512
     & k0*(cf12(i1,i2)%e(2)*p5(3)-p5(2)*cf12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p5k0+p5(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf12(i1,i2)%e(0)*p5(0)-cf12(i1,i2)%e(1)*p5(1)-cf12(i1
     & ,i2)%e(2)*p5(2)-cf12(i1,i2)%e(3)*p5(3)
      cvqd=cf12(i1,i2)%e(0)*p512(0)-cf12(i1,i2)%e(1)*p512(1)-cf1
     & 2(i1,i2)%e(2)*p512(2)-cf12(i1,i2)%e(3)*p512(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p5k0*cvqd+p512k0*cvqu
      cauxc=+cf12(i1,i2)%ek0*p5(2)-p5k0*cf12(i1,i2)%e(2)
      l5_12(i1,i2)%a(2)=l5_12(i1,i2)%a(2)+ccl*(cauxa-ceps_0)
      l5_12(i1,i2)%c(2)=l5_12(i1,i2)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
*zw                                                                     
      endif ! if(iresw(5).eq.1)
*zwend                                                                  
      do m=0,3
        p534(m)=p5(m)+p34(m)
      enddo
* pk0 -- p=p534
      p534k0=p534(0)-p534(1)
* p.q -- p.q=p534q,p=p534,q=p534,bef=,aft=
      p534q=(p534(0)*p534(0)-p534(1)*p534(1)-p534(2)*p534(2)-p53
     & 4(3)*p534(3))
*zw                                                                     
      if(iresw(5).eq.1) then
        do i1=1,2
          do i2=1,2
            l5_34(i1,i2)%a(2)=czero
            l5_34(i1,i2)%c(2)=czero
          enddo
        enddo
      else
*zwend                                                                  
* quqd -- p=p5,q=p534
      quqd=p5(0)*p534(0)-p5(1)*p534(1)-p5(2)*p534(2)-p5(3)*p534(
     & 3)
      ccl=zcl(id5)/(p534q*p534k0)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p534,v=cz34(i1,i2)%e,a=l5_34(i1,i2)%a,c=l5_34(i1,i2)%
* c,cl=ccl,nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0
     & *(cz34(i1,i2)%e(2)*p534(3)-p534(2)*cz34(i1,i2)%e(3))-p534
     & k0*(cz34(i1,i2)%e(2)*p5(3)-p5(2)*cz34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p5k0+p5(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i1,i2)%e(0)*p5(0)-cz34(i1,i2)%e(1)*p5(1)-cz34(i1
     & ,i2)%e(2)*p5(2)-cz34(i1,i2)%e(3)*p5(3)
      cvqd=cz34(i1,i2)%e(0)*p534(0)-cz34(i1,i2)%e(1)*p534(1)-cz3
     & 4(i1,i2)%e(2)*p534(2)-cz34(i1,i2)%e(3)*p534(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cz34(i1,i2)%ek0*p5(2)-p5k0*cz34(i1,i2)%e(2)
      l5_34(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l5_34(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      ccl=fcl(id5)/(p534q*p534k0)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p5,qd=p534,v=cf34(i1,i2)%e,a=l5_34(i1,i2)%a,c=l5_34(i1,i2)%
* c,cl=ccl,nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p5(2)*p534(3)-p534(2)*p5(3))+p5k0
     & *(cf34(i1,i2)%e(2)*p534(3)-p534(2)*cf34(i1,i2)%e(3))-p534
     & k0*(cf34(i1,i2)%e(2)*p5(3)-p5(2)*cf34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p5k0+p5(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf34(i1,i2)%e(0)*p5(0)-cf34(i1,i2)%e(1)*p5(1)-cf34(i1
     & ,i2)%e(2)*p5(2)-cf34(i1,i2)%e(3)*p5(3)
      cvqd=cf34(i1,i2)%e(0)*p534(0)-cf34(i1,i2)%e(1)*p534(1)-cf3
     & 4(i1,i2)%e(2)*p534(2)-cf34(i1,i2)%e(3)*p534(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p5k0*cvqd+p534k0*cvqu
      cauxc=+cf34(i1,i2)%ek0*p5(2)-p5k0*cf34(i1,i2)%e(2)
      l5_34(i1,i2)%a(2)=l5_34(i1,i2)%a(2)+ccl*(cauxa-ceps_0)
      l5_34(i1,i2)%c(2)=l5_34(i1,i2)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
*zw                                                                     
      endif ! if(iresw(5).eq.1)
*zwend                                                                  
      do m=0,3
        p712(m)=p7(m)+p12(m)
      enddo
* pk0 -- p=p712
      p712k0=p712(0)-p712(1)
* p.q -- p.q=p712q,p=p712,q=p712,bef=,aft=
      p712q=(p712(0)*p712(0)-p712(1)*p712(1)-p712(2)*p712(2)-p71
     & 2(3)*p712(3))
*zw                                                                     
      if(iresw(7).eq.1) then
        do i1=1,2
          do i2=1,2
            l7_12(i1,i2)%a(2)=czero
            l7_12(i1,i2)%c(2)=czero
          enddo
        enddo
      else
*zwend                                                                  
* quqd -- p=p7,q=p712
      quqd=p7(0)*p712(0)-p7(1)*p712(1)-p7(2)*p712(2)-p7(3)*p712(
     & 3)
      ccl=zcl(id7)/(p712q*p712k0)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p712,v=cz12(i1,i2)%e,a=l7_12(i1,i2)%a,c=l7_12(i1,i2)%
* c,cl=ccl,nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p7(2)*p712(3)-p712(2)*p7(3))+p7k0
     & *(cz12(i1,i2)%e(2)*p712(3)-p712(2)*cz12(i1,i2)%e(3))-p712
     & k0*(cz12(i1,i2)%e(2)*p7(3)-p7(2)*cz12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p7k0+p7(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz12(i1,i2)%e(0)*p7(0)-cz12(i1,i2)%e(1)*p7(1)-cz12(i1
     & ,i2)%e(2)*p7(2)-cz12(i1,i2)%e(3)*p7(3)
      cvqd=cz12(i1,i2)%e(0)*p712(0)-cz12(i1,i2)%e(1)*p712(1)-cz1
     & 2(i1,i2)%e(2)*p712(2)-cz12(i1,i2)%e(3)*p712(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p7k0*cvqd+p712k0*cvqu
      cauxc=+cz12(i1,i2)%ek0*p7(2)-p7k0*cz12(i1,i2)%e(2)
      l7_12(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_12(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      ccl=fcl(id7)/(p712q*p712k0)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p712,v=cf12(i1,i2)%e,a=l7_12(i1,i2)%a,c=l7_12(i1,i2)%
* c,cl=ccl,nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p7(2)*p712(3)-p712(2)*p7(3))+p7k0
     & *(cf12(i1,i2)%e(2)*p712(3)-p712(2)*cf12(i1,i2)%e(3))-p712
     & k0*(cf12(i1,i2)%e(2)*p7(3)-p7(2)*cf12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p7k0+p7(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf12(i1,i2)%e(0)*p7(0)-cf12(i1,i2)%e(1)*p7(1)-cf12(i1
     & ,i2)%e(2)*p7(2)-cf12(i1,i2)%e(3)*p7(3)
      cvqd=cf12(i1,i2)%e(0)*p712(0)-cf12(i1,i2)%e(1)*p712(1)-cf1
     & 2(i1,i2)%e(2)*p712(2)-cf12(i1,i2)%e(3)*p712(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p7k0*cvqd+p712k0*cvqu
      cauxc=+cf12(i1,i2)%ek0*p7(2)-p7k0*cf12(i1,i2)%e(2)
      l7_12(i1,i2)%a(2)=l7_12(i1,i2)%a(2)+ccl*(cauxa-ceps_0)
      l7_12(i1,i2)%c(2)=l7_12(i1,i2)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
*zw                                                                     
      endif ! if(iresw(7).eq.1)
*zwend                                                                  
      do m=0,3
        p734(m)=p7(m)+p34(m)
      enddo
* pk0 -- p=p734
      p734k0=p734(0)-p734(1)
* p.q -- p.q=p734q,p=p734,q=p734,bef=,aft=
      p734q=(p734(0)*p734(0)-p734(1)*p734(1)-p734(2)*p734(2)-p73
     & 4(3)*p734(3))
*zw                                                                     
      if(iresw(7).eq.1) then
        do i1=1,2
          do i2=1,2
            l7_34(i1,i2)%a(2)=czero
            l7_34(i1,i2)%c(2)=czero
          enddo
        enddo
      else
*zwend                                                                  
* quqd -- p=p7,q=p734
      quqd=p7(0)*p734(0)-p7(1)*p734(1)-p7(2)*p734(2)-p7(3)*p734(
     & 3)
      ccl=zcl(id7)/(p734q*p734k0)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p734,v=cz34(i1,i2)%e,a=l7_34(i1,i2)%a,c=l7_34(i1,i2)%
* c,cl=ccl,nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0
     & *(cz34(i1,i2)%e(2)*p734(3)-p734(2)*cz34(i1,i2)%e(3))-p734
     & k0*(cz34(i1,i2)%e(2)*p7(3)-p7(2)*cz34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p7k0+p7(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cz34(i1,i2)%e(0)*p7(0)-cz34(i1,i2)%e(1)*p7(1)-cz34(i1
     & ,i2)%e(2)*p7(2)-cz34(i1,i2)%e(3)*p7(3)
      cvqd=cz34(i1,i2)%e(0)*p734(0)-cz34(i1,i2)%e(1)*p734(1)-cz3
     & 4(i1,i2)%e(2)*p734(2)-cz34(i1,i2)%e(3)*p734(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cz34(i1,i2)%ek0*p7(2)-p7k0*cz34(i1,i2)%e(2)
      l7_34(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      l7_34(i1,i2)%c(2)=ccl*(-cauxc+ceps_1)
      end do
      end do
      ccl=fcl(id7)/(p734q*p734k0)
      do i1=1,2
      do i2=1,2
* TWL0 -- qu=p7,qd=p734,v=cf34(i1,i2)%e,a=l7_34(i1,i2)%a,c=l7_34(i1,i2)%
* c,cl=ccl,nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0
     & *(cf34(i1,i2)%e(2)*p734(3)-p734(2)*cf34(i1,i2)%e(3))-p734
     & k0*(cf34(i1,i2)%e(2)*p7(3)-p7(2)*cf34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p7k0+p7(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      cvqu=cf34(i1,i2)%e(0)*p7(0)-cf34(i1,i2)%e(1)*p7(1)-cf34(i1
     & ,i2)%e(2)*p7(2)-cf34(i1,i2)%e(3)*p7(3)
      cvqd=cf34(i1,i2)%e(0)*p734(0)-cf34(i1,i2)%e(1)*p734(1)-cf3
     & 4(i1,i2)%e(2)*p734(2)-cf34(i1,i2)%e(3)*p734(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxc=+cf34(i1,i2)%ek0*p7(2)-p7k0*cf34(i1,i2)%e(2)
      l7_34(i1,i2)%a(2)=l7_34(i1,i2)%a(2)+ccl*(cauxa-ceps_0)
      l7_34(i1,i2)%c(2)=l7_34(i1,i2)%c(2)+ccl*(-cauxc+ceps_1)
      end do
      end do
*zw                                                                     
      endif ! if(iresw(7).eq.1)
*zwend                                                                  
  
  
* compute all single insertions of the type li_ (i=1,3)  for a Zline    
* and a w  insertion                                                    
*        i __                                                           
*            |_W__/                                                     
*            |    \                                                     
*                                                                       
*together with its propagator and pk0                                   
  
  
      if (iup(id1).eq.1) then  ! W is W- : 78
        do m=0,3
          p178(m)=p1(m)+p78(m)
        enddo
* pk0 -- p=p178
      p178k0=p178(0)-p178(1)
* p.q -- p.q=p178q,p=p178,q=p178,bef=,aft=
      p178q=(p178(0)*p178(0)-p178(1)*p178(1)-p178(2)*p178(2)-p17
     & 8(3)*p178(3))
  
*zw                                                                     
        if(iresz(1).eq.1) then
          l1_78%a(2,2)=czero
          l1_78%b(1,2)=czero
          l1_78%c(2,1)=czero
          l1_78%d(1,1)=czero
        else
* quqd -- p=p1,q=p178
      quqd=p1(0)*p178(0)-p1(1)*p178(1)-p1(2)*p178(2)-p1(3)*p178(
     & 3)
      ccl=wcl/((p178q-rmass2(id1-1))*p178k0)
* TW -- qu=p1,qd=p178,v=cw78%e,a=l1_78%a,b=l1_78%b,c=l1_78%c,d=l1_78%d,c
* l=ccl,nsum=0
      ceps_0=-cw78%ek0*(p1(2)*p178(3)-p178(2)*p1(3))+p1k0*(cw78%
     & e(2)*p178(3)-p178(2)*cw78%e(3))-p178k0*(cw78%e(2)*p1(3)-p
     & 1(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p1k0+p1(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p178k0+p178(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p1(0)-cw78%e(1)*p1(1)-cw78%e(2)*p1(2)-cw78%
     & e(3)*p1(3)
      cvqd=cw78%e(0)*p178(0)-cw78%e(1)*p178(1)-cw78%e(2)*p178(2)
     & -cw78%e(3)*p178(3)
      cauxa=-cw78%ek0*quqd+p1k0*cvqd+p178k0*cvqu
      cauxb=-cw78%ek0*p178(2)+p178k0*cw78%e(2)
      cauxc=+cw78%ek0*p1(2)-p1k0*cw78%e(2)
      l1_78%a(2,2)=ccl*(cauxa-ceps_0)
      l1_78%b(1,2)=ccl*(cauxb-ceps_2)
      l1_78%c(2,1)=ccl*(-cauxc+ceps_1)
      l1_78%d(1,1)=ccl*cw78%ek0
        endif
*zwend                                                                  
      else  ! W is W+ : 56
        do m=0,3
          p156(m)=p1(m)+p56(m)
        enddo
* pk0 -- p=p156
      p156k0=p156(0)-p156(1)
* p.q -- p.q=p156q,p=p156,q=p156,bef=,aft=
      p156q=(p156(0)*p156(0)-p156(1)*p156(1)-p156(2)*p156(2)-p15
     & 6(3)*p156(3))
  
*zw                                                                     
        if(iresz(1).eq.1) then
          l1_56%a(2,2)=czero
          l1_56%b(1,2)=czero
          l1_56%c(2,1)=czero
          l1_56%d(1,1)=czero
        else
* quqd -- p=p1,q=p156
      quqd=p1(0)*p156(0)-p1(1)*p156(1)-p1(2)*p156(2)-p1(3)*p156(
     & 3)
      ccl=wcl/((p156q-cmass2(id1+1))*p156k0)
* TW -- qu=p1,qd=p156,v=cw56%e,a=l1_56%a,b=l1_56%b,c=l1_56%c,d=l1_56%d,c
* l=ccl,nsum=0
      ceps_0=-cw56%ek0*(p1(2)*p156(3)-p156(2)*p1(3))+p1k0*(cw56%
     & e(2)*p156(3)-p156(2)*cw56%e(3))-p156k0*(cw56%e(2)*p1(3)-p
     & 1(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p1k0+p1(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p156k0+p156(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p1(0)-cw56%e(1)*p1(1)-cw56%e(2)*p1(2)-cw56%
     & e(3)*p1(3)
      cvqd=cw56%e(0)*p156(0)-cw56%e(1)*p156(1)-cw56%e(2)*p156(2)
     & -cw56%e(3)*p156(3)
      cauxa=-cw56%ek0*quqd+p1k0*cvqd+p156k0*cvqu
      cauxb=-cw56%ek0*p156(2)+p156k0*cw56%e(2)
      cauxc=+cw56%ek0*p1(2)-p1k0*cw56%e(2)
      l1_56%a(2,2)=ccl*(cauxa-ceps_0)
      l1_56%b(1,2)=ccl*(cauxb-ceps_2)
      l1_56%c(2,1)=ccl*(-cauxc+ceps_1)
      l1_56%d(1,1)=ccl*cw56%ek0
        endif
*zwend                                                                  
      endif
  
      if (iup(id3).eq.1) then  ! W is W- : 78
        do m=0,3
          p378(m)=p3(m)+p78(m)
        enddo
* pk0 -- p=p378
      p378k0=p378(0)-p378(1)
* p.q -- p.q=p378q,p=p378,q=p378,bef=,aft=
      p378q=(p378(0)*p378(0)-p378(1)*p378(1)-p378(2)*p378(2)-p37
     & 8(3)*p378(3))
  
*zw                                                                     
        if(iresz(3).eq.1) then
          l3_78%a(2,2)=czero
          l3_78%b(1,2)=czero
          l3_78%c(2,1)=czero
          l3_78%d(1,1)=czero
        else
* quqd -- p=p3,q=p378
      quqd=p3(0)*p378(0)-p3(1)*p378(1)-p3(2)*p378(2)-p3(3)*p378(
     & 3)
      ccl=wcl/((p378q-rmass2(id3-1))*p378k0)
* TW -- qu=p3,qd=p378,v=cw78%e,a=l3_78%a,b=l3_78%b,c=l3_78%c,d=l3_78%d,c
* l=ccl,nsum=0
      ceps_0=-cw78%ek0*(p3(2)*p378(3)-p378(2)*p3(3))+p3k0*(cw78%
     & e(2)*p378(3)-p378(2)*cw78%e(3))-p378k0*(cw78%e(2)*p3(3)-p
     & 3(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p3k0+p3(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p378k0+p378(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p3(0)-cw78%e(1)*p3(1)-cw78%e(2)*p3(2)-cw78%
     & e(3)*p3(3)
      cvqd=cw78%e(0)*p378(0)-cw78%e(1)*p378(1)-cw78%e(2)*p378(2)
     & -cw78%e(3)*p378(3)
      cauxa=-cw78%ek0*quqd+p3k0*cvqd+p378k0*cvqu
      cauxb=-cw78%ek0*p378(2)+p378k0*cw78%e(2)
      cauxc=+cw78%ek0*p3(2)-p3k0*cw78%e(2)
      l3_78%a(2,2)=ccl*(cauxa-ceps_0)
      l3_78%b(1,2)=ccl*(cauxb-ceps_2)
      l3_78%c(2,1)=ccl*(-cauxc+ceps_1)
      l3_78%d(1,1)=ccl*cw78%ek0
        endif
*zwend                                                                  
      else  ! W is W+ : 56
        do m=0,3
          p356(m)=p3(m)+p56(m)
        enddo
* pk0 -- p=p356
      p356k0=p356(0)-p356(1)
* p.q -- p.q=p356q,p=p356,q=p356,bef=,aft=
      p356q=(p356(0)*p356(0)-p356(1)*p356(1)-p356(2)*p356(2)-p35
     & 6(3)*p356(3))
  
*zw                                                                     
        if(iresz(3).eq.1) then
          l3_56%a(2,2)=czero
          l3_56%b(1,2)=czero
          l3_56%c(2,1)=czero
          l3_56%d(1,1)=czero
        else
* quqd -- p=p3,q=p356
      quqd=p3(0)*p356(0)-p3(1)*p356(1)-p3(2)*p356(2)-p3(3)*p356(
     & 3)
      ccl=wcl/((p356q-cmass2(id3+1))*p356k0)
* TW -- qu=p3,qd=p356,v=cw56%e,a=l3_56%a,b=l3_56%b,c=l3_56%c,d=l3_56%d,c
* l=ccl,nsum=0
      ceps_0=-cw56%ek0*(p3(2)*p356(3)-p356(2)*p3(3))+p3k0*(cw56%
     & e(2)*p356(3)-p356(2)*cw56%e(3))-p356k0*(cw56%e(2)*p3(3)-p
     & 3(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p3k0+p3(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p356k0+p356(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p3(0)-cw56%e(1)*p3(1)-cw56%e(2)*p3(2)-cw56%
     & e(3)*p3(3)
      cvqd=cw56%e(0)*p356(0)-cw56%e(1)*p356(1)-cw56%e(2)*p356(2)
     & -cw56%e(3)*p356(3)
      cauxa=-cw56%ek0*quqd+p3k0*cvqd+p356k0*cvqu
      cauxb=-cw56%ek0*p356(2)+p356k0*cw56%e(2)
      cauxc=+cw56%ek0*p3(2)-p3k0*cw56%e(2)
      l3_56%a(2,2)=ccl*(cauxa-ceps_0)
      l3_56%b(1,2)=ccl*(cauxb-ceps_2)
      l3_56%c(2,1)=ccl*(-cauxc+ceps_1)
      l3_56%d(1,1)=ccl*cw56%ek0
        endif
*zwend                                                                  
      endif
  
* compute all single insertions of the type li_ (i=1,3)  for a Zline    
* and a Z  gamma and eventually higgs insertion                         
*        i __                                                           
*            |_Z,f,h_/                                                  
*            |       \                                                  
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
  
*zw                                                                     
  
      if(iresz(1).eq.1) then
        do i1=1,2
          do i2=1,2
            l1_34(i1,i2)%a(i3,i4)=czero
            l1_34(i1,i2)%b(i3,i4)=czero
            l1_34(i1,i2)%c(i3,i4)=czero
            l1_34(i1,i2)%d(i3,i4)=czero
          enddo
        enddo
      else
*zwend                                                                  
  
* quqd -- p=p1,q=p134
      quqd=p1(0)*p134(0)-p1(1)*p134(1)-p1(2)*p134(2)-p1(3)*p134(
     & 3)
      ccr=zcr(id1)/((p134q-rmass2(id1))*p134k0)
      ccl=zcl(id1)/((p134q-rmass2(id1))*p134k0)
      do i1=1,2
      do i2=1,2
* T -- qu=p1,qd=p134,v=cz34(i1,i2)%e,a=l1_34(i1,i2)%a,b=l1_34(i1,i2)%b,c
* =l1_34(i1,i2)%c,d=l1_34(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p1(2)*p134(3)-p134(2)*p1(3))+p1k0
     & *(cz34(i1,i2)%e(2)*p134(3)-p134(2)*cz34(i1,i2)%e(3))-p134
     & k0*(cz34(i1,i2)%e(2)*p1(3)-p1(2)*cz34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p1k0+p1(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i1,i2)%e(3)*p134k0+p134(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p1(0)-cz34(i1,i2)%e(1)*p1(1)-cz34(i1
     & ,i2)%e(2)*p1(2)-cz34(i1,i2)%e(3)*p1(3)
      cvqd=cz34(i1,i2)%e(0)*p134(0)-cz34(i1,i2)%e(1)*p134(1)-cz3
     & 4(i1,i2)%e(2)*p134(2)-cz34(i1,i2)%e(3)*p134(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p1k0*cvqd+p134k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p134(2)+p134k0*cz34(i1,i2)%e(2)
      cauxc=+cz34(i1,i2)%ek0*p1(2)-p1k0*cz34(i1,i2)%e(2)
      l1_34(i1,i2)%a(1,1)=ccr*(cauxa+ceps_0)
      l1_34(i1,i2)%a(2,2)=ccl*(cauxa-ceps_0)
      l1_34(i1,i2)%b(1,2)=ccl*(cauxb-ceps_2)
      l1_34(i1,i2)%b(2,1)=ccr*(-cauxb-ceps_2)
      l1_34(i1,i2)%c(1,2)=ccr*(cauxc+ceps_1)
      l1_34(i1,i2)%c(2,1)=ccl*(-cauxc+ceps_1)
      l1_34(i1,i2)%d(1,1)=ccl*cz34(i1,i2)%ek0
      l1_34(i1,i2)%d(2,2)=ccr*cz34(i1,i2)%ek0
      end do
      end do
      ccr=fcr(id1)/((p134q-rmass2(id1))*p134k0)
      ccl=fcl(id1)/((p134q-rmass2(id1))*p134k0)
      do i1=1,2
      do i2=1,2
* T -- qu=p1,qd=p134,v=cf34(i1,i2)%e,a=l1_34(i1,i2)%a,b=l1_34(i1,i2)%b,c
* =l1_34(i1,i2)%c,d=l1_34(i1,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p1(2)*p134(3)-p134(2)*p1(3))+p1k0
     & *(cf34(i1,i2)%e(2)*p134(3)-p134(2)*cf34(i1,i2)%e(3))-p134
     & k0*(cf34(i1,i2)%e(2)*p1(3)-p1(2)*cf34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p1k0+p1(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i1,i2)%e(3)*p134k0+p134(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p1(0)-cf34(i1,i2)%e(1)*p1(1)-cf34(i1
     & ,i2)%e(2)*p1(2)-cf34(i1,i2)%e(3)*p1(3)
      cvqd=cf34(i1,i2)%e(0)*p134(0)-cf34(i1,i2)%e(1)*p134(1)-cf3
     & 4(i1,i2)%e(2)*p134(2)-cf34(i1,i2)%e(3)*p134(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p1k0*cvqd+p134k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p134(2)+p134k0*cf34(i1,i2)%e(2)
      cauxc=+cf34(i1,i2)%ek0*p1(2)-p1k0*cf34(i1,i2)%e(2)
      l1_34(i1,i2)%a(1,1)=l1_34(i1,i2)%a(1,1)+ccr*(cauxa+ceps_0)
      l1_34(i1,i2)%a(2,2)=l1_34(i1,i2)%a(2,2)+ccl*(cauxa-ceps_0)
      l1_34(i1,i2)%b(1,2)=l1_34(i1,i2)%b(1,2)+ccl*(cauxb-ceps_2)
      l1_34(i1,i2)%b(2,1)=l1_34(i1,i2)%b(2,1)+ccr*(-cauxb-ceps_2
     & )
      l1_34(i1,i2)%c(1,2)=l1_34(i1,i2)%c(1,2)+ccr*(cauxc+ceps_1)
      l1_34(i1,i2)%c(2,1)=l1_34(i1,i2)%c(2,1)+ccl*(-cauxc+ceps_1
     & )
      l1_34(i1,i2)%d(1,1)=l1_34(i1,i2)%d(1,1)+ccl*cf34(i1,i2)%ek
     & 0
      l1_34(i1,i2)%d(2,2)=l1_34(i1,i2)%d(2,2)+ccr*cf34(i1,i2)%ek
     & 0
      end do
      end do
  
* rmh < 0 in our convention means no higgs coupling                     
      if (id1.eq.5.and.id3.eq.5.and.rmh.ge.0.d0) then
  
* TH -- qu=p1,qd=p134,a=clineth%a,b=clineth%b,c=clineth%c,coupl=1.
      auxa=-p1k0*p134(2)+p134k0*p1(2)
      cauxa=auxa-cim*(p134(3)*p1k0-p1(3)*p134k0)
      clineth%a(1,2)=cauxa
      clineth%a(2,1)=-conjg(cauxa)
      clineth%b(1,1)=p134k0
      clineth%b(2,2)=clineth%b(1,1)
      clineth%c(1,1)=p1k0
      clineth%c(2,2)=clineth%c(1,1)
  
      cffactor=rhbb**2/(s34-cmh2)+rhhbb**2/(s34-cmhh2)
      cfactor= cffactor/((p134q-rmass2(id1))*p134k0)
  
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  l1_34(i1,i2)%a(i3,i4)=ch34(i1,i2)*clineth%a(i3,i4)
     &                 *cfactor
                else
                  l1_34(i1,i2)%b(i3,i4)=ch34(i1,i2)*clineth%b(i3,i4)
     &                 *cfactor
                  l1_34(i1,i2)%c(i3,i4)=ch34(i1,i2)*clineth%c(i3,i4)
     &                 *cfactor
                endif
              enddo
            enddo
          enddo
        enddo
      else
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  l1_34(i1,i2)%a(i3,i4)=czero
                else
                  l1_34(i1,i2)%b(i3,i4)=czero
                  l1_34(i1,i2)%c(i3,i4)=czero
                endif
              enddo
            enddo
          enddo
        enddo
      endif
*zw                                                                     
      endif   !  if(iresz(1).eq.1)
*zwend                                                                  
  
      do m=0,3
        p312(m)=p3(m)+p12(m)
      enddo
* pk0 -- p=p312
      p312k0=p312(0)-p312(1)
* p.q -- p.q=p312q,p=p312,q=p312,bef=,aft=
      p312q=(p312(0)*p312(0)-p312(1)*p312(1)-p312(2)*p312(2)-p31
     & 2(3)*p312(3))
  
*zw                                                                     
  
      if(iresz(3).eq.1) then
        do i1=1,2
          do i2=1,2
            l3_12(i1,i2)%a(i3,i4)=czero
            l3_12(i1,i2)%b(i3,i4)=czero
            l3_12(i1,i2)%c(i3,i4)=czero
            l3_12(i1,i2)%d(i3,i4)=czero
          enddo
        enddo
      else
*zwend                                                                  
  
* quqd -- p=p3,q=p312
      quqd=p3(0)*p312(0)-p3(1)*p312(1)-p3(2)*p312(2)-p3(3)*p312(
     & 3)
      ccr=zcr(id3)/((p312q-rmass2(id3))*p312k0)
      ccl=zcl(id3)/((p312q-rmass2(id3))*p312k0)
      do i1=1,2
      do i2=1,2
* T -- qu=p3,qd=p312,v=cz12(i1,i2)%e,a=l3_12(i1,i2)%a,b=l3_12(i1,i2)%b,c
* =l3_12(i1,i2)%c,d=l3_12(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p3(2)*p312(3)-p312(2)*p3(3))+p3k0
     & *(cz12(i1,i2)%e(2)*p312(3)-p312(2)*cz12(i1,i2)%e(3))-p312
     & k0*(cz12(i1,i2)%e(2)*p3(3)-p3(2)*cz12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p3k0+p3(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i1,i2)%e(3)*p312k0+p312(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p3(0)-cz12(i1,i2)%e(1)*p3(1)-cz12(i1
     & ,i2)%e(2)*p3(2)-cz12(i1,i2)%e(3)*p3(3)
      cvqd=cz12(i1,i2)%e(0)*p312(0)-cz12(i1,i2)%e(1)*p312(1)-cz1
     & 2(i1,i2)%e(2)*p312(2)-cz12(i1,i2)%e(3)*p312(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p3k0*cvqd+p312k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p312(2)+p312k0*cz12(i1,i2)%e(2)
      cauxc=+cz12(i1,i2)%ek0*p3(2)-p3k0*cz12(i1,i2)%e(2)
      l3_12(i1,i2)%a(1,1)=ccr*(cauxa+ceps_0)
      l3_12(i1,i2)%a(2,2)=ccl*(cauxa-ceps_0)
      l3_12(i1,i2)%b(1,2)=ccl*(cauxb-ceps_2)
      l3_12(i1,i2)%b(2,1)=ccr*(-cauxb-ceps_2)
      l3_12(i1,i2)%c(1,2)=ccr*(cauxc+ceps_1)
      l3_12(i1,i2)%c(2,1)=ccl*(-cauxc+ceps_1)
      l3_12(i1,i2)%d(1,1)=ccl*cz12(i1,i2)%ek0
      l3_12(i1,i2)%d(2,2)=ccr*cz12(i1,i2)%ek0
      end do
      end do
      ccr=fcr(id3)/((p312q-rmass2(id3))*p312k0)
      ccl=fcl(id3)/((p312q-rmass2(id3))*p312k0)
      do i1=1,2
      do i2=1,2
* T -- qu=p3,qd=p312,v=cf12(i1,i2)%e,a=l3_12(i1,i2)%a,b=l3_12(i1,i2)%b,c
* =l3_12(i1,i2)%c,d=l3_12(i1,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p3(2)*p312(3)-p312(2)*p3(3))+p3k0
     & *(cf12(i1,i2)%e(2)*p312(3)-p312(2)*cf12(i1,i2)%e(3))-p312
     & k0*(cf12(i1,i2)%e(2)*p3(3)-p3(2)*cf12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p3k0+p3(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i1,i2)%e(3)*p312k0+p312(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p3(0)-cf12(i1,i2)%e(1)*p3(1)-cf12(i1
     & ,i2)%e(2)*p3(2)-cf12(i1,i2)%e(3)*p3(3)
      cvqd=cf12(i1,i2)%e(0)*p312(0)-cf12(i1,i2)%e(1)*p312(1)-cf1
     & 2(i1,i2)%e(2)*p312(2)-cf12(i1,i2)%e(3)*p312(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p3k0*cvqd+p312k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p312(2)+p312k0*cf12(i1,i2)%e(2)
      cauxc=+cf12(i1,i2)%ek0*p3(2)-p3k0*cf12(i1,i2)%e(2)
      l3_12(i1,i2)%a(1,1)=l3_12(i1,i2)%a(1,1)+ccr*(cauxa+ceps_0)
      l3_12(i1,i2)%a(2,2)=l3_12(i1,i2)%a(2,2)+ccl*(cauxa-ceps_0)
      l3_12(i1,i2)%b(1,2)=l3_12(i1,i2)%b(1,2)+ccl*(cauxb-ceps_2)
      l3_12(i1,i2)%b(2,1)=l3_12(i1,i2)%b(2,1)+ccr*(-cauxb-ceps_2
     & )
      l3_12(i1,i2)%c(1,2)=l3_12(i1,i2)%c(1,2)+ccr*(cauxc+ceps_1)
      l3_12(i1,i2)%c(2,1)=l3_12(i1,i2)%c(2,1)+ccl*(-cauxc+ceps_1
     & )
      l3_12(i1,i2)%d(1,1)=l3_12(i1,i2)%d(1,1)+ccl*cf12(i1,i2)%ek
     & 0
      l3_12(i1,i2)%d(2,2)=l3_12(i1,i2)%d(2,2)+ccr*cf12(i1,i2)%ek
     & 0
      end do
      end do
  
* rmh < 0 in our convention means no higgs coupling                     
      if (id1.eq.5.and.id3.eq.5.and.rmh.ge.0.d0) then
  
* TH -- qu=p3,qd=p312,a=clineth%a,b=clineth%b,c=clineth%c,coupl=1.
      auxa=-p3k0*p312(2)+p312k0*p3(2)
      cauxa=auxa-cim*(p312(3)*p3k0-p3(3)*p312k0)
      clineth%a(1,2)=cauxa
      clineth%a(2,1)=-conjg(cauxa)
      clineth%b(1,1)=p312k0
      clineth%b(2,2)=clineth%b(1,1)
      clineth%c(1,1)=p3k0
      clineth%c(2,2)=clineth%c(1,1)
  
      cffactor=rhbb**2/(s12-cmh2)+rhhbb**2/(s12-cmhh2)
      cfactor= cffactor/((p312q-rmass2(id3))*p312k0)
  
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  l3_12(i1,i2)%a(i3,i4)=ch12(i1,i2)*clineth%a(i3,i4)
     &                 *cfactor
                else
                  l3_12(i1,i2)%b(i3,i4)=ch12(i1,i2)*clineth%b(i3,i4)
     &                 *cfactor
                  l3_12(i1,i2)%c(i3,i4)=ch12(i1,i2)*clineth%c(i3,i4)
     &                 *cfactor
                endif
              enddo
            enddo
          enddo
        enddo
      else
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  l3_12(i1,i2)%a(i3,i4)=czero
                else
                  l3_12(i1,i2)%b(i3,i4)=czero
                  l3_12(i1,i2)%c(i3,i4)=czero
                endif
              enddo
            enddo
          enddo
        enddo
      endif
*zw                                                                     
      endif   !  if(iresz(3).eq.1)
*zwend                                                                  
  
  
  
* compute all single insertions of the type ri_ (i=6,8)  for a Wline    
* and a W insertion                                                     
*                                                                       
*            |_W__/                                                     
*        i __|    \                                                     
*                                                                       
*together with its propagator and pk0                                   
  
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
  
* compute all single insertions of the type ri_ (i=6,8)  for a Wline    
* and a Z and gamma insertion                                           
*                                                                       
*            |_Z,f_/                                                    
*         i__|     \                                                    
*                                                                       
*together with its propagator and pk0                                   
  
      do m=0,3
        p612(m)=-p6(m)-p12(m)
      enddo
* pk0 -- p=p612
      p612k0=p612(0)-p612(1)
* p.q -- p.q=p612q,p=p612,q=p612,bef=,aft=
      p612q=(p612(0)*p612(0)-p612(1)*p612(1)-p612(2)*p612(2)-p61
     & 2(3)*p612(3))
* quqd -- p=p612,q=p6
      quqd=p612(0)*p6(0)-p612(1)*p6(1)-p612(2)*p6(2)-p612(3)*p6(
     & 3)
      ccl=zcl(id6)/(p612q*p612k0)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p612,qd=p6,v=cz12(i1,i2)%e,a=r6_12(i1,i2)%a,b=r6_12(i1,i2)%
* b,cl=ccl,nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p612(2)*p6(3)-p6(2)*p612(3))+p612
     & k0*(cz12(i1,i2)%e(2)*p6(3)-p6(2)*cz12(i1,i2)%e(3))-p6k0*(
     & cz12(i1,i2)%e(2)*p612(3)-p612(2)*cz12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz12(i1,i2)%e(3)*p6k0+p6(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p612(0)-cz12(i1,i2)%e(1)*p612(1)-cz1
     & 2(i1,i2)%e(2)*p612(2)-cz12(i1,i2)%e(3)*p612(3)
      cvqd=cz12(i1,i2)%e(0)*p6(0)-cz12(i1,i2)%e(1)*p6(1)-cz12(i1
     & ,i2)%e(2)*p6(2)-cz12(i1,i2)%e(3)*p6(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p612k0*cvqd+p6k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p6(2)+p6k0*cz12(i1,i2)%e(2)
      r6_12(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      r6_12(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      end do
      end do
      ccl=fcl(id6)/(p612q*p612k0)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p612,qd=p6,v=cf12(i1,i2)%e,a=r6_12(i1,i2)%a,b=r6_12(i1,i2)%
* b,cl=ccl,nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p612(2)*p6(3)-p6(2)*p612(3))+p612
     & k0*(cf12(i1,i2)%e(2)*p6(3)-p6(2)*cf12(i1,i2)%e(3))-p6k0*(
     & cf12(i1,i2)%e(2)*p612(3)-p612(2)*cf12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf12(i1,i2)%e(3)*p6k0+p6(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p612(0)-cf12(i1,i2)%e(1)*p612(1)-cf1
     & 2(i1,i2)%e(2)*p612(2)-cf12(i1,i2)%e(3)*p612(3)
      cvqd=cf12(i1,i2)%e(0)*p6(0)-cf12(i1,i2)%e(1)*p6(1)-cf12(i1
     & ,i2)%e(2)*p6(2)-cf12(i1,i2)%e(3)*p6(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p612k0*cvqd+p6k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p6(2)+p6k0*cf12(i1,i2)%e(2)
      r6_12(i1,i2)%a(2)=r6_12(i1,i2)%a(2)+ccl*(cauxa-ceps_0)
      r6_12(i1,i2)%b(1)=r6_12(i1,i2)%b(1)+ccl*(cauxb-ceps_2)
      end do
      end do
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
      ccl=zcl(id6)/(p634q*p634k0)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p634,qd=p6,v=cz34(i1,i2)%e,a=r6_34(i1,i2)%a,b=r6_34(i1,i2)%
* b,cl=ccl,nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634
     & k0*(cz34(i1,i2)%e(2)*p6(3)-p6(2)*cz34(i1,i2)%e(3))-p6k0*(
     & cz34(i1,i2)%e(2)*p634(3)-p634(2)*cz34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i1,i2)%e(3)*p6k0+p6(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p634(0)-cz34(i1,i2)%e(1)*p634(1)-cz3
     & 4(i1,i2)%e(2)*p634(2)-cz34(i1,i2)%e(3)*p634(3)
      cvqd=cz34(i1,i2)%e(0)*p6(0)-cz34(i1,i2)%e(1)*p6(1)-cz34(i1
     & ,i2)%e(2)*p6(2)-cz34(i1,i2)%e(3)*p6(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p6(2)+p6k0*cz34(i1,i2)%e(2)
      r6_34(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      r6_34(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      end do
      end do
      ccl=fcl(id6)/(p634q*p634k0)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p634,qd=p6,v=cf34(i1,i2)%e,a=r6_34(i1,i2)%a,b=r6_34(i1,i2)%
* b,cl=ccl,nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p634(2)*p6(3)-p6(2)*p634(3))+p634
     & k0*(cf34(i1,i2)%e(2)*p6(3)-p6(2)*cf34(i1,i2)%e(3))-p6k0*(
     & cf34(i1,i2)%e(2)*p634(3)-p634(2)*cf34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf34(i1,i2)%e(3)*p6k0+p6(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p634(0)-cf34(i1,i2)%e(1)*p634(1)-cf3
     & 4(i1,i2)%e(2)*p634(2)-cf34(i1,i2)%e(3)*p634(3)
      cvqd=cf34(i1,i2)%e(0)*p6(0)-cf34(i1,i2)%e(1)*p6(1)-cf34(i1
     & ,i2)%e(2)*p6(2)-cf34(i1,i2)%e(3)*p6(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p634k0*cvqd+p6k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p6(2)+p6k0*cf34(i1,i2)%e(2)
      r6_34(i1,i2)%a(2)=r6_34(i1,i2)%a(2)+ccl*(cauxa-ceps_0)
      r6_34(i1,i2)%b(1)=r6_34(i1,i2)%b(1)+ccl*(cauxb-ceps_2)
      end do
      end do
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
      ccl=zcl(id8)/(p812q*p812k0)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p812,qd=p8,v=cz12(i1,i2)%e,a=r8_12(i1,i2)%a,b=r8_12(i1,i2)%
* b,cl=ccl,nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p812(2)*p8(3)-p8(2)*p812(3))+p812
     & k0*(cz12(i1,i2)%e(2)*p8(3)-p8(2)*cz12(i1,i2)%e(3))-p8k0*(
     & cz12(i1,i2)%e(2)*p812(3)-p812(2)*cz12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz12(i1,i2)%e(3)*p8k0+p8(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p812(0)-cz12(i1,i2)%e(1)*p812(1)-cz1
     & 2(i1,i2)%e(2)*p812(2)-cz12(i1,i2)%e(3)*p812(3)
      cvqd=cz12(i1,i2)%e(0)*p8(0)-cz12(i1,i2)%e(1)*p8(1)-cz12(i1
     & ,i2)%e(2)*p8(2)-cz12(i1,i2)%e(3)*p8(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p812k0*cvqd+p8k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p8(2)+p8k0*cz12(i1,i2)%e(2)
      r8_12(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      r8_12(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      end do
      end do
      ccl=fcl(id8)/(p812q*p812k0)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p812,qd=p8,v=cf12(i1,i2)%e,a=r8_12(i1,i2)%a,b=r8_12(i1,i2)%
* b,cl=ccl,nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p812(2)*p8(3)-p8(2)*p812(3))+p812
     & k0*(cf12(i1,i2)%e(2)*p8(3)-p8(2)*cf12(i1,i2)%e(3))-p8k0*(
     & cf12(i1,i2)%e(2)*p812(3)-p812(2)*cf12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf12(i1,i2)%e(3)*p8k0+p8(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p812(0)-cf12(i1,i2)%e(1)*p812(1)-cf1
     & 2(i1,i2)%e(2)*p812(2)-cf12(i1,i2)%e(3)*p812(3)
      cvqd=cf12(i1,i2)%e(0)*p8(0)-cf12(i1,i2)%e(1)*p8(1)-cf12(i1
     & ,i2)%e(2)*p8(2)-cf12(i1,i2)%e(3)*p8(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p812k0*cvqd+p8k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p8(2)+p8k0*cf12(i1,i2)%e(2)
      r8_12(i1,i2)%a(2)=r8_12(i1,i2)%a(2)+ccl*(cauxa-ceps_0)
      r8_12(i1,i2)%b(1)=r8_12(i1,i2)%b(1)+ccl*(cauxb-ceps_2)
      end do
      end do
      do m=0,3
        p834(m)=-p8(m)-p34(m)
      enddo
* pk0 -- p=p834
      p834k0=p834(0)-p834(1)
* p.q -- p.q=p834q,p=p834,q=p834,bef=,aft=
      p834q=(p834(0)*p834(0)-p834(1)*p834(1)-p834(2)*p834(2)-p83
     & 4(3)*p834(3))
* quqd -- p=p834,q=p8
      quqd=p834(0)*p8(0)-p834(1)*p8(1)-p834(2)*p8(2)-p834(3)*p8(
     & 3)
      ccl=zcl(id8)/(p834q*p834k0)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p834,qd=p8,v=cz34(i1,i2)%e,a=r8_34(i1,i2)%a,b=r8_34(i1,i2)%
* b,cl=ccl,nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834
     & k0*(cz34(i1,i2)%e(2)*p8(3)-p8(2)*cz34(i1,i2)%e(3))-p8k0*(
     & cz34(i1,i2)%e(2)*p834(3)-p834(2)*cz34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cz34(i1,i2)%e(3)*p8k0+p8(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p834(0)-cz34(i1,i2)%e(1)*p834(1)-cz3
     & 4(i1,i2)%e(2)*p834(2)-cz34(i1,i2)%e(3)*p834(3)
      cvqd=cz34(i1,i2)%e(0)*p8(0)-cz34(i1,i2)%e(1)*p8(1)-cz34(i1
     & ,i2)%e(2)*p8(2)-cz34(i1,i2)%e(3)*p8(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p8(2)+p8k0*cz34(i1,i2)%e(2)
      r8_34(i1,i2)%a(2)=ccl*(cauxa-ceps_0)
      r8_34(i1,i2)%b(1)=ccl*(cauxb-ceps_2)
      end do
      end do
      ccl=fcl(id8)/(p834q*p834k0)
      do i1=1,2
      do i2=1,2
* TWR0 -- qu=p834,qd=p8,v=cf34(i1,i2)%e,a=r8_34(i1,i2)%a,b=r8_34(i1,i2)%
* b,cl=ccl,nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p834(2)*p8(3)-p8(2)*p834(3))+p834
     & k0*(cf34(i1,i2)%e(2)*p8(3)-p8(2)*cf34(i1,i2)%e(3))-p8k0*(
     & cf34(i1,i2)%e(2)*p834(3)-p834(2)*cf34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_2=-cf34(i1,i2)%e(3)*p8k0+p8(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p834(0)-cf34(i1,i2)%e(1)*p834(1)-cf3
     & 4(i1,i2)%e(2)*p834(2)-cf34(i1,i2)%e(3)*p834(3)
      cvqd=cf34(i1,i2)%e(0)*p8(0)-cf34(i1,i2)%e(1)*p8(1)-cf34(i1
     & ,i2)%e(2)*p8(2)-cf34(i1,i2)%e(3)*p8(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p834k0*cvqd+p8k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p8(2)+p8k0*cf34(i1,i2)%e(2)
      r8_34(i1,i2)%a(2)=r8_34(i1,i2)%a(2)+ccl*(cauxa-ceps_0)
      r8_34(i1,i2)%b(1)=r8_34(i1,i2)%b(1)+ccl*(cauxb-ceps_2)
      end do
      end do
  
  
* compute all single insertions of the type ri_ (i=2,4)  for a Zline    
* and a W  insertion                                                    
*                                                                       
*            |_W__/                                                     
*         i__|    \                                                     
*                                                                       
*together with its propagator and pk0                                   
  
  
      if (iup(id2).eq.1) then ! W is W+ : 56
        do m=0,3
          p256(m)=-p2(m)-p56(m)
        enddo
* pk0 -- p=p256
      p256k0=p256(0)-p256(1)
* p.q -- p.q=p256q,p=p256,q=p256,bef=,aft=
      p256q=(p256(0)*p256(0)-p256(1)*p256(1)-p256(2)*p256(2)-p25
     & 6(3)*p256(3))
* id2 has a negative value                                              
* quqd -- p=p256,q=p2
      quqd=p256(0)*p2(0)-p256(1)*p2(1)-p256(2)*p2(2)-p256(3)*p2(
     & 3)
      ccl=wcl/((p256q-rmass2(id2+1))*p256k0)
* TW -- qu=p256,qd=p2,v=cw56%e,a=r2_56%a,b=r2_56%b,c=r2_56%c,d=r2_56%d,c
* l=ccl,nsum=0
      ceps_0=-cw56%ek0*(p256(2)*p2(3)-p2(2)*p256(3))+p256k0*(cw5
     & 6%e(2)*p2(3)-p2(2)*cw56%e(3))-p2k0*(cw56%e(2)*p256(3)-p25
     & 6(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p256k0+p256(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p2k0+p2(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p256(0)-cw56%e(1)*p256(1)-cw56%e(2)*p256(2)
     & -cw56%e(3)*p256(3)
      cvqd=cw56%e(0)*p2(0)-cw56%e(1)*p2(1)-cw56%e(2)*p2(2)-cw56%
     & e(3)*p2(3)
      cauxa=-cw56%ek0*quqd+p256k0*cvqd+p2k0*cvqu
      cauxb=-cw56%ek0*p2(2)+p2k0*cw56%e(2)
      cauxc=+cw56%ek0*p256(2)-p256k0*cw56%e(2)
      r2_56%a(2,2)=ccl*(cauxa-ceps_0)
      r2_56%b(1,2)=ccl*(cauxb-ceps_2)
      r2_56%c(2,1)=ccl*(-cauxc+ceps_1)
      r2_56%d(1,1)=ccl*cw56%ek0
      else   ! W is W- : 78
        do m=0,3
          p278(m)=-p2(m)-p78(m)
        enddo
* pk0 -- p=p278
      p278k0=p278(0)-p278(1)
* p.q -- p.q=p278q,p=p278,q=p278,bef=,aft=
      p278q=(p278(0)*p278(0)-p278(1)*p278(1)-p278(2)*p278(2)-p27
     & 8(3)*p278(3))
* id2 has a negative value                                              
* quqd -- p=p278,q=p2
      quqd=p278(0)*p2(0)-p278(1)*p2(1)-p278(2)*p2(2)-p278(3)*p2(
     & 3)
      ccl=wcl/((p278q-cmass2(id2-1))*p278k0)
* TW -- qu=p278,qd=p2,v=cw78%e,a=r2_78%a,b=r2_78%b,c=r2_78%c,d=r2_78%d,c
* l=ccl,nsum=0
      ceps_0=-cw78%ek0*(p278(2)*p2(3)-p2(2)*p278(3))+p278k0*(cw7
     & 8%e(2)*p2(3)-p2(2)*cw78%e(3))-p2k0*(cw78%e(2)*p278(3)-p27
     & 8(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p278k0+p278(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p2k0+p2(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p278(0)-cw78%e(1)*p278(1)-cw78%e(2)*p278(2)
     & -cw78%e(3)*p278(3)
      cvqd=cw78%e(0)*p2(0)-cw78%e(1)*p2(1)-cw78%e(2)*p2(2)-cw78%
     & e(3)*p2(3)
      cauxa=-cw78%ek0*quqd+p278k0*cvqd+p2k0*cvqu
      cauxb=-cw78%ek0*p2(2)+p2k0*cw78%e(2)
      cauxc=+cw78%ek0*p278(2)-p278k0*cw78%e(2)
      r2_78%a(2,2)=ccl*(cauxa-ceps_0)
      r2_78%b(1,2)=ccl*(cauxb-ceps_2)
      r2_78%c(2,1)=ccl*(-cauxc+ceps_1)
      r2_78%d(1,1)=ccl*cw78%ek0
  
      endif
  
      if (iup(id4).eq.1) then ! W is W+ : 56
        do m=0,3
          p456(m)=-p4(m)-p56(m)
        enddo
* pk0 -- p=p456
      p456k0=p456(0)-p456(1)
* p.q -- p.q=p456q,p=p456,q=p456,bef=,aft=
      p456q=(p456(0)*p456(0)-p456(1)*p456(1)-p456(2)*p456(2)-p45
     & 6(3)*p456(3))
* id4 has a negative value                                              
* quqd -- p=p456,q=p4
      quqd=p456(0)*p4(0)-p456(1)*p4(1)-p456(2)*p4(2)-p456(3)*p4(
     & 3)
      ccl=wcl/((p456q-rmass2(id4+1))*p456k0)
* TW -- qu=p456,qd=p4,v=cw56%e,a=r4_56%a,b=r4_56%b,c=r4_56%c,d=r4_56%d,c
* l=ccl,nsum=0
      ceps_0=-cw56%ek0*(p456(2)*p4(3)-p4(2)*p456(3))+p456k0*(cw5
     & 6%e(2)*p4(3)-p4(2)*cw56%e(3))-p4k0*(cw56%e(2)*p456(3)-p45
     & 6(2)*cw56%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56%e(3)*p456k0+p456(3)*cw56%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56%e(3)*p4k0+p4(3)*cw56%ek0
      ceps_2=ceps_2*cim
      cvqu=cw56%e(0)*p456(0)-cw56%e(1)*p456(1)-cw56%e(2)*p456(2)
     & -cw56%e(3)*p456(3)
      cvqd=cw56%e(0)*p4(0)-cw56%e(1)*p4(1)-cw56%e(2)*p4(2)-cw56%
     & e(3)*p4(3)
      cauxa=-cw56%ek0*quqd+p456k0*cvqd+p4k0*cvqu
      cauxb=-cw56%ek0*p4(2)+p4k0*cw56%e(2)
      cauxc=+cw56%ek0*p456(2)-p456k0*cw56%e(2)
      r4_56%a(2,2)=ccl*(cauxa-ceps_0)
      r4_56%b(1,2)=ccl*(cauxb-ceps_2)
      r4_56%c(2,1)=ccl*(-cauxc+ceps_1)
      r4_56%d(1,1)=ccl*cw56%ek0
      else   ! W is W- : 78
        do m=0,3
          p478(m)=-p4(m)-p78(m)
        enddo
* pk0 -- p=p478
      p478k0=p478(0)-p478(1)
* p.q -- p.q=p478q,p=p478,q=p478,bef=,aft=
      p478q=(p478(0)*p478(0)-p478(1)*p478(1)-p478(2)*p478(2)-p47
     & 8(3)*p478(3))
* id4 has a negative value                                              
* quqd -- p=p478,q=p4
      quqd=p478(0)*p4(0)-p478(1)*p4(1)-p478(2)*p4(2)-p478(3)*p4(
     & 3)
      ccl=wcl/((p478q-cmass2(id4-1))*p478k0)
* TW -- qu=p478,qd=p4,v=cw78%e,a=r4_78%a,b=r4_78%b,c=r4_78%c,d=r4_78%d,c
* l=ccl,nsum=0
      ceps_0=-cw78%ek0*(p478(2)*p4(3)-p4(2)*p478(3))+p478k0*(cw7
     & 8%e(2)*p4(3)-p4(2)*cw78%e(3))-p4k0*(cw78%e(2)*p478(3)-p47
     & 8(2)*cw78%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw78%e(3)*p478k0+p478(3)*cw78%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw78%e(3)*p4k0+p4(3)*cw78%ek0
      ceps_2=ceps_2*cim
      cvqu=cw78%e(0)*p478(0)-cw78%e(1)*p478(1)-cw78%e(2)*p478(2)
     & -cw78%e(3)*p478(3)
      cvqd=cw78%e(0)*p4(0)-cw78%e(1)*p4(1)-cw78%e(2)*p4(2)-cw78%
     & e(3)*p4(3)
      cauxa=-cw78%ek0*quqd+p478k0*cvqd+p4k0*cvqu
      cauxb=-cw78%ek0*p4(2)+p4k0*cw78%e(2)
      cauxc=+cw78%ek0*p478(2)-p478k0*cw78%e(2)
      r4_78%a(2,2)=ccl*(cauxa-ceps_0)
      r4_78%b(1,2)=ccl*(cauxb-ceps_2)
      r4_78%c(2,1)=ccl*(-cauxc+ceps_1)
      r4_78%d(1,1)=ccl*cw78%ek0
  
      endif
  
* compute all single insertions of the type ri_ (i=2,4)  for a Zline    
* and a Z  gamma and eventually higgs insertion                         
*                                                                       
*            |_Z,f,h_/                                                  
*         i__|       \                                                  
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
      ccr=zcr(id2)/((p234q-rmass2(id2))*p234k0)
      ccl=zcl(id2)/((p234q-rmass2(id2))*p234k0)
      do i1=1,2
      do i2=1,2
* T -- qu=p234,qd=p2,v=cz34(i1,i2)%e,a=r2_34(i1,i2)%a,b=r2_34(i1,i2)%b,c
* =r2_34(i1,i2)%c,d=r2_34(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p234(2)*p2(3)-p2(2)*p234(3))+p234
     & k0*(cz34(i1,i2)%e(2)*p2(3)-p2(2)*cz34(i1,i2)%e(3))-p2k0*(
     & cz34(i1,i2)%e(2)*p234(3)-p234(2)*cz34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p234k0+p234(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i1,i2)%e(3)*p2k0+p2(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p234(0)-cz34(i1,i2)%e(1)*p234(1)-cz3
     & 4(i1,i2)%e(2)*p234(2)-cz34(i1,i2)%e(3)*p234(3)
      cvqd=cz34(i1,i2)%e(0)*p2(0)-cz34(i1,i2)%e(1)*p2(1)-cz34(i1
     & ,i2)%e(2)*p2(2)-cz34(i1,i2)%e(3)*p2(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p234k0*cvqd+p2k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p2(2)+p2k0*cz34(i1,i2)%e(2)
      cauxc=+cz34(i1,i2)%ek0*p234(2)-p234k0*cz34(i1,i2)%e(2)
      r2_34(i1,i2)%a(1,1)=ccr*(cauxa+ceps_0)
      r2_34(i1,i2)%a(2,2)=ccl*(cauxa-ceps_0)
      r2_34(i1,i2)%b(1,2)=ccl*(cauxb-ceps_2)
      r2_34(i1,i2)%b(2,1)=ccr*(-cauxb-ceps_2)
      r2_34(i1,i2)%c(1,2)=ccr*(cauxc+ceps_1)
      r2_34(i1,i2)%c(2,1)=ccl*(-cauxc+ceps_1)
      r2_34(i1,i2)%d(1,1)=ccl*cz34(i1,i2)%ek0
      r2_34(i1,i2)%d(2,2)=ccr*cz34(i1,i2)%ek0
      end do
      end do
      ccr=fcr(id2)/((p234q-rmass2(id2))*p234k0)
      ccl=fcl(id2)/((p234q-rmass2(id2))*p234k0)
      do i1=1,2
      do i2=1,2
* T -- qu=p234,qd=p2,v=cf34(i1,i2)%e,a=r2_34(i1,i2)%a,b=r2_34(i1,i2)%b,c
* =r2_34(i1,i2)%c,d=r2_34(i1,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p234(2)*p2(3)-p2(2)*p234(3))+p234
     & k0*(cf34(i1,i2)%e(2)*p2(3)-p2(2)*cf34(i1,i2)%e(3))-p2k0*(
     & cf34(i1,i2)%e(2)*p234(3)-p234(2)*cf34(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p234k0+p234(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i1,i2)%e(3)*p2k0+p2(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p234(0)-cf34(i1,i2)%e(1)*p234(1)-cf3
     & 4(i1,i2)%e(2)*p234(2)-cf34(i1,i2)%e(3)*p234(3)
      cvqd=cf34(i1,i2)%e(0)*p2(0)-cf34(i1,i2)%e(1)*p2(1)-cf34(i1
     & ,i2)%e(2)*p2(2)-cf34(i1,i2)%e(3)*p2(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p234k0*cvqd+p2k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p2(2)+p2k0*cf34(i1,i2)%e(2)
      cauxc=+cf34(i1,i2)%ek0*p234(2)-p234k0*cf34(i1,i2)%e(2)
      r2_34(i1,i2)%a(1,1)=r2_34(i1,i2)%a(1,1)+ccr*(cauxa+ceps_0)
      r2_34(i1,i2)%a(2,2)=r2_34(i1,i2)%a(2,2)+ccl*(cauxa-ceps_0)
      r2_34(i1,i2)%b(1,2)=r2_34(i1,i2)%b(1,2)+ccl*(cauxb-ceps_2)
      r2_34(i1,i2)%b(2,1)=r2_34(i1,i2)%b(2,1)+ccr*(-cauxb-ceps_2
     & )
      r2_34(i1,i2)%c(1,2)=r2_34(i1,i2)%c(1,2)+ccr*(cauxc+ceps_1)
      r2_34(i1,i2)%c(2,1)=r2_34(i1,i2)%c(2,1)+ccl*(-cauxc+ceps_1
     & )
      r2_34(i1,i2)%d(1,1)=r2_34(i1,i2)%d(1,1)+ccl*cf34(i1,i2)%ek
     & 0
      r2_34(i1,i2)%d(2,2)=r2_34(i1,i2)%d(2,2)+ccr*cf34(i1,i2)%ek
     & 0
      end do
      end do
  
* rmh < 0 in our convention means no higgs coupling                     
      if (id1.eq.5.and.id3.eq.5.and.rmh.ge.0.d0) then
  
* TH -- qu=p234,qd=p2,a=clineth%a,b=clineth%b,c=clineth%c,coupl=1.
      auxa=-p234k0*p2(2)+p2k0*p234(2)
      cauxa=auxa-cim*(p2(3)*p234k0-p234(3)*p2k0)
      clineth%a(1,2)=cauxa
      clineth%a(2,1)=-conjg(cauxa)
      clineth%b(1,1)=p2k0
      clineth%b(2,2)=clineth%b(1,1)
      clineth%c(1,1)=p234k0
      clineth%c(2,2)=clineth%c(1,1)
  
      cffactor=rhbb**2/(s34-cmh2)+rhhbb**2/(s34-cmhh2)
      cfactor= cffactor/((p234q-rmass2(id2))*p234k0)
  
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  r2_34(i1,i2)%a(i3,i4)=ch34(i1,i2)*clineth%a(i3,i4)
     &                 *cfactor
                else
                  r2_34(i1,i2)%b(i3,i4)=ch34(i1,i2)*clineth%b(i3,i4)
     &                 *cfactor
                  r2_34(i1,i2)%c(i3,i4)=ch34(i1,i2)*clineth%c(i3,i4)
     &                 *cfactor
                endif
              enddo
            enddo
          enddo
        enddo
      else
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  r2_34(i1,i2)%a(i3,i4)=czero
                else
                  r2_34(i1,i2)%b(i3,i4)=czero
                  r2_34(i1,i2)%c(i3,i4)=czero
                endif
              enddo
            enddo
          enddo
        enddo
      endif
  
  
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
      ccr=zcr(id4)/((p412q-rmass2(id4))*p412k0)
      ccl=zcl(id4)/((p412q-rmass2(id4))*p412k0)
      do i1=1,2
      do i2=1,2
* T -- qu=p412,qd=p4,v=cz12(i1,i2)%e,a=r4_12(i1,i2)%a,b=r4_12(i1,i2)%b,c
* =r4_12(i1,i2)%c,d=r4_12(i1,i2)%d,cr=ccr,cl=ccl,nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p412(2)*p4(3)-p4(2)*p412(3))+p412
     & k0*(cz12(i1,i2)%e(2)*p4(3)-p4(2)*cz12(i1,i2)%e(3))-p4k0*(
     & cz12(i1,i2)%e(2)*p412(3)-p412(2)*cz12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p412k0+p412(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i1,i2)%e(3)*p4k0+p4(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p412(0)-cz12(i1,i2)%e(1)*p412(1)-cz1
     & 2(i1,i2)%e(2)*p412(2)-cz12(i1,i2)%e(3)*p412(3)
      cvqd=cz12(i1,i2)%e(0)*p4(0)-cz12(i1,i2)%e(1)*p4(1)-cz12(i1
     & ,i2)%e(2)*p4(2)-cz12(i1,i2)%e(3)*p4(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p412k0*cvqd+p4k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p4(2)+p4k0*cz12(i1,i2)%e(2)
      cauxc=+cz12(i1,i2)%ek0*p412(2)-p412k0*cz12(i1,i2)%e(2)
      r4_12(i1,i2)%a(1,1)=ccr*(cauxa+ceps_0)
      r4_12(i1,i2)%a(2,2)=ccl*(cauxa-ceps_0)
      r4_12(i1,i2)%b(1,2)=ccl*(cauxb-ceps_2)
      r4_12(i1,i2)%b(2,1)=ccr*(-cauxb-ceps_2)
      r4_12(i1,i2)%c(1,2)=ccr*(cauxc+ceps_1)
      r4_12(i1,i2)%c(2,1)=ccl*(-cauxc+ceps_1)
      r4_12(i1,i2)%d(1,1)=ccl*cz12(i1,i2)%ek0
      r4_12(i1,i2)%d(2,2)=ccr*cz12(i1,i2)%ek0
      end do
      end do
      ccr=fcr(id4)/((p412q-rmass2(id4))*p412k0)
      ccl=fcl(id4)/((p412q-rmass2(id4))*p412k0)
      do i1=1,2
      do i2=1,2
* T -- qu=p412,qd=p4,v=cf12(i1,i2)%e,a=r4_12(i1,i2)%a,b=r4_12(i1,i2)%b,c
* =r4_12(i1,i2)%c,d=r4_12(i1,i2)%d,cr=ccr,cl=ccl,nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p412(2)*p4(3)-p4(2)*p412(3))+p412
     & k0*(cf12(i1,i2)%e(2)*p4(3)-p4(2)*cf12(i1,i2)%e(3))-p4k0*(
     & cf12(i1,i2)%e(2)*p412(3)-p412(2)*cf12(i1,i2)%e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p412k0+p412(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i1,i2)%e(3)*p4k0+p4(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p412(0)-cf12(i1,i2)%e(1)*p412(1)-cf1
     & 2(i1,i2)%e(2)*p412(2)-cf12(i1,i2)%e(3)*p412(3)
      cvqd=cf12(i1,i2)%e(0)*p4(0)-cf12(i1,i2)%e(1)*p4(1)-cf12(i1
     & ,i2)%e(2)*p4(2)-cf12(i1,i2)%e(3)*p4(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p412k0*cvqd+p4k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p4(2)+p4k0*cf12(i1,i2)%e(2)
      cauxc=+cf12(i1,i2)%ek0*p412(2)-p412k0*cf12(i1,i2)%e(2)
      r4_12(i1,i2)%a(1,1)=r4_12(i1,i2)%a(1,1)+ccr*(cauxa+ceps_0)
      r4_12(i1,i2)%a(2,2)=r4_12(i1,i2)%a(2,2)+ccl*(cauxa-ceps_0)
      r4_12(i1,i2)%b(1,2)=r4_12(i1,i2)%b(1,2)+ccl*(cauxb-ceps_2)
      r4_12(i1,i2)%b(2,1)=r4_12(i1,i2)%b(2,1)+ccr*(-cauxb-ceps_2
     & )
      r4_12(i1,i2)%c(1,2)=r4_12(i1,i2)%c(1,2)+ccr*(cauxc+ceps_1)
      r4_12(i1,i2)%c(2,1)=r4_12(i1,i2)%c(2,1)+ccl*(-cauxc+ceps_1
     & )
      r4_12(i1,i2)%d(1,1)=r4_12(i1,i2)%d(1,1)+ccl*cf12(i1,i2)%ek
     & 0
      r4_12(i1,i2)%d(2,2)=r4_12(i1,i2)%d(2,2)+ccr*cf12(i1,i2)%ek
     & 0
      end do
      end do
  
* rmh < 0 in our convention means no higgs coupling                     
      if (id1.eq.5.and.id3.eq.5.and.rmh.ge.0.d0) then
  
* TH -- qu=p412,qd=p4,a=clineth%a,b=clineth%b,c=clineth%c,coupl=1.
      auxa=-p412k0*p4(2)+p4k0*p412(2)
      cauxa=auxa-cim*(p4(3)*p412k0-p412(3)*p4k0)
      clineth%a(1,2)=cauxa
      clineth%a(2,1)=-conjg(cauxa)
      clineth%b(1,1)=p4k0
      clineth%b(2,2)=clineth%b(1,1)
      clineth%c(1,1)=p412k0
      clineth%c(2,2)=clineth%c(1,1)
  
      cffactor=rhbb**2/(s12-cmh2)+rhhbb**2/(s12-cmhh2)
      cfactor= cffactor/((p412q-rmass2(id4))*p412k0)
  
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  r4_12(i1,i2)%a(i3,i4)=ch12(i1,i2)*clineth%a(i3,i4)
     &                 *cfactor
                else
                  r4_12(i1,i2)%b(i3,i4)=ch12(i1,i2)*clineth%b(i3,i4)
     &                 *cfactor
                  r4_12(i1,i2)%c(i3,i4)=ch12(i1,i2)*clineth%c(i3,i4)
     &                 *cfactor
                endif
              enddo
            enddo
          enddo
        enddo
      else
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  r4_12(i1,i2)%a(i3,i4)=czero
                else
                  r4_12(i1,i2)%b(i3,i4)=czero
                  r4_12(i1,i2)%c(i3,i4)=czero
                endif
              enddo
            enddo
          enddo
        enddo
      endif
  
  
  
  
  
* compute all single W insertions in the middle of a triple insertion   
*  Wline                                                                
*                                                                       
*            |                                                          
*            |_W__/                                                     
*            |    \                                                     
*            |                                                          
  
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
  
  
  
* compute all single Z or gamma insertions in the middle of a           
*  triple insertion Wline                                               
*                                                                       
*            |                                                          
*            |_Z,f__/                                                   
*            |      \                                                   
*            |                                                          
  
  
* quqd -- p=p512,q=p678
      quqd=p512(0)*p678(0)-p512(1)*p678(1)-p512(2)*p678(2)-p512(
     & 3)*p678(3)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p512,qd=p678,v=cz34(i1,i2)%e,a=u512_34(i1,i2)%a,b=u512_34(i1
* ,i2)%b,c=u512_34(i1,i2)%c,d=u512_34(i1,i2)%d,cl=zcl(id5),nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p512(2)*p678(3)-p678(2)*p512(3))+
     & p512k0*(cz34(i1,i2)%e(2)*p678(3)-p678(2)*cz34(i1,i2)%e(3)
     & )-p678k0*(cz34(i1,i2)%e(2)*p512(3)-p512(2)*cz34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p512k0+p512(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i1,i2)%e(3)*p678k0+p678(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p512(0)-cz34(i1,i2)%e(1)*p512(1)-cz3
     & 4(i1,i2)%e(2)*p512(2)-cz34(i1,i2)%e(3)*p512(3)
      cvqd=cz34(i1,i2)%e(0)*p678(0)-cz34(i1,i2)%e(1)*p678(1)-cz3
     & 4(i1,i2)%e(2)*p678(2)-cz34(i1,i2)%e(3)*p678(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p512k0*cvqd+p678k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p678(2)+p678k0*cz34(i1,i2)%e(2)
      cauxc=+cz34(i1,i2)%ek0*p512(2)-p512k0*cz34(i1,i2)%e(2)
      u512_34(i1,i2)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u512_34(i1,i2)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u512_34(i1,i2)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u512_34(i1,i2)%d(1)=zcl(id5)*cz34(i1,i2)%ek0
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p512,qd=p678,v=cf34(i1,i2)%e,a=u512_34(i1,i2)%a,b=u512_34(i1
* ,i2)%b,c=u512_34(i1,i2)%c,d=u512_34(i1,i2)%d,cl=fcl(id5),nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p512(2)*p678(3)-p678(2)*p512(3))+
     & p512k0*(cf34(i1,i2)%e(2)*p678(3)-p678(2)*cf34(i1,i2)%e(3)
     & )-p678k0*(cf34(i1,i2)%e(2)*p512(3)-p512(2)*cf34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p512k0+p512(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i1,i2)%e(3)*p678k0+p678(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p512(0)-cf34(i1,i2)%e(1)*p512(1)-cf3
     & 4(i1,i2)%e(2)*p512(2)-cf34(i1,i2)%e(3)*p512(3)
      cvqd=cf34(i1,i2)%e(0)*p678(0)-cf34(i1,i2)%e(1)*p678(1)-cf3
     & 4(i1,i2)%e(2)*p678(2)-cf34(i1,i2)%e(3)*p678(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p512k0*cvqd+p678k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p678(2)+p678k0*cf34(i1,i2)%e(2)
      cauxc=+cf34(i1,i2)%ek0*p512(2)-p512k0*cf34(i1,i2)%e(2)
      u512_34(i1,i2)%a(2)=u512_34(i1,i2)%a(2)+fcl(id5)*(cauxa-ce
     & ps_0)
      u512_34(i1,i2)%b(1)=u512_34(i1,i2)%b(1)+fcl(id5)*(cauxb-ce
     & ps_2)
      u512_34(i1,i2)%c(2)=u512_34(i1,i2)%c(2)+fcl(id5)*(-cauxc+c
     & eps_1)
      u512_34(i1,i2)%d(1)=u512_34(i1,i2)%d(1)+fcl(id5)*cf34(i1,i
     & 2)%ek0
      end do
      end do
  
* quqd -- p=p534,q=p678
      quqd=p534(0)*p678(0)-p534(1)*p678(1)-p534(2)*p678(2)-p534(
     & 3)*p678(3)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p534,qd=p678,v=cz12(i1,i2)%e,a=u534_12(i1,i2)%a,b=u534_12(i1
* ,i2)%b,c=u534_12(i1,i2)%c,d=u534_12(i1,i2)%d,cl=zcl(id5),nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p534(2)*p678(3)-p678(2)*p534(3))+
     & p534k0*(cz12(i1,i2)%e(2)*p678(3)-p678(2)*cz12(i1,i2)%e(3)
     & )-p678k0*(cz12(i1,i2)%e(2)*p534(3)-p534(2)*cz12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p534k0+p534(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i1,i2)%e(3)*p678k0+p678(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p534(0)-cz12(i1,i2)%e(1)*p534(1)-cz1
     & 2(i1,i2)%e(2)*p534(2)-cz12(i1,i2)%e(3)*p534(3)
      cvqd=cz12(i1,i2)%e(0)*p678(0)-cz12(i1,i2)%e(1)*p678(1)-cz1
     & 2(i1,i2)%e(2)*p678(2)-cz12(i1,i2)%e(3)*p678(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p534k0*cvqd+p678k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p678(2)+p678k0*cz12(i1,i2)%e(2)
      cauxc=+cz12(i1,i2)%ek0*p534(2)-p534k0*cz12(i1,i2)%e(2)
      u534_12(i1,i2)%a(2)=zcl(id5)*(cauxa-ceps_0)
      u534_12(i1,i2)%b(1)=zcl(id5)*(cauxb-ceps_2)
      u534_12(i1,i2)%c(2)=zcl(id5)*(-cauxc+ceps_1)
      u534_12(i1,i2)%d(1)=zcl(id5)*cz12(i1,i2)%ek0
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p534,qd=p678,v=cf12(i1,i2)%e,a=u534_12(i1,i2)%a,b=u534_12(i1
* ,i2)%b,c=u534_12(i1,i2)%c,d=u534_12(i1,i2)%d,cl=fcl(id5),nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p534(2)*p678(3)-p678(2)*p534(3))+
     & p534k0*(cf12(i1,i2)%e(2)*p678(3)-p678(2)*cf12(i1,i2)%e(3)
     & )-p678k0*(cf12(i1,i2)%e(2)*p534(3)-p534(2)*cf12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p534k0+p534(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i1,i2)%e(3)*p678k0+p678(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p534(0)-cf12(i1,i2)%e(1)*p534(1)-cf1
     & 2(i1,i2)%e(2)*p534(2)-cf12(i1,i2)%e(3)*p534(3)
      cvqd=cf12(i1,i2)%e(0)*p678(0)-cf12(i1,i2)%e(1)*p678(1)-cf1
     & 2(i1,i2)%e(2)*p678(2)-cf12(i1,i2)%e(3)*p678(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p534k0*cvqd+p678k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p678(2)+p678k0*cf12(i1,i2)%e(2)
      cauxc=+cf12(i1,i2)%ek0*p534(2)-p534k0*cf12(i1,i2)%e(2)
      u534_12(i1,i2)%a(2)=u534_12(i1,i2)%a(2)+fcl(id5)*(cauxa-ce
     & ps_0)
      u534_12(i1,i2)%b(1)=u534_12(i1,i2)%b(1)+fcl(id5)*(cauxb-ce
     & ps_2)
      u534_12(i1,i2)%c(2)=u534_12(i1,i2)%c(2)+fcl(id5)*(-cauxc+c
     & eps_1)
      u534_12(i1,i2)%d(1)=u534_12(i1,i2)%d(1)+fcl(id5)*cf12(i1,i
     & 2)%ek0
      end do
      end do
  
* quqd -- p=p578,q=p634
      quqd=p578(0)*p634(0)-p578(1)*p634(1)-p578(2)*p634(2)-p578(
     & 3)*p634(3)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p578,qd=p634,v=cz12(i1,i2)%e,a=u578_12(i1,i2)%a,b=u578_12(i1
* ,i2)%b,c=u578_12(i1,i2)%c,d=u578_12(i1,i2)%d,cl=zcl(id6),nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p578(2)*p634(3)-p634(2)*p578(3))+
     & p578k0*(cz12(i1,i2)%e(2)*p634(3)-p634(2)*cz12(i1,i2)%e(3)
     & )-p634k0*(cz12(i1,i2)%e(2)*p578(3)-p578(2)*cz12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p578k0+p578(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i1,i2)%e(3)*p634k0+p634(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p578(0)-cz12(i1,i2)%e(1)*p578(1)-cz1
     & 2(i1,i2)%e(2)*p578(2)-cz12(i1,i2)%e(3)*p578(3)
      cvqd=cz12(i1,i2)%e(0)*p634(0)-cz12(i1,i2)%e(1)*p634(1)-cz1
     & 2(i1,i2)%e(2)*p634(2)-cz12(i1,i2)%e(3)*p634(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p578k0*cvqd+p634k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p634(2)+p634k0*cz12(i1,i2)%e(2)
      cauxc=+cz12(i1,i2)%ek0*p578(2)-p578k0*cz12(i1,i2)%e(2)
      u578_12(i1,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      u578_12(i1,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      u578_12(i1,i2)%c(2)=zcl(id6)*(-cauxc+ceps_1)
      u578_12(i1,i2)%d(1)=zcl(id6)*cz12(i1,i2)%ek0
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p578,qd=p634,v=cf12(i1,i2)%e,a=u578_12(i1,i2)%a,b=u578_12(i1
* ,i2)%b,c=u578_12(i1,i2)%c,d=u578_12(i1,i2)%d,cl=fcl(id6),nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p578(2)*p634(3)-p634(2)*p578(3))+
     & p578k0*(cf12(i1,i2)%e(2)*p634(3)-p634(2)*cf12(i1,i2)%e(3)
     & )-p634k0*(cf12(i1,i2)%e(2)*p578(3)-p578(2)*cf12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p578k0+p578(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i1,i2)%e(3)*p634k0+p634(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p578(0)-cf12(i1,i2)%e(1)*p578(1)-cf1
     & 2(i1,i2)%e(2)*p578(2)-cf12(i1,i2)%e(3)*p578(3)
      cvqd=cf12(i1,i2)%e(0)*p634(0)-cf12(i1,i2)%e(1)*p634(1)-cf1
     & 2(i1,i2)%e(2)*p634(2)-cf12(i1,i2)%e(3)*p634(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p578k0*cvqd+p634k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p634(2)+p634k0*cf12(i1,i2)%e(2)
      cauxc=+cf12(i1,i2)%ek0*p578(2)-p578k0*cf12(i1,i2)%e(2)
      u578_12(i1,i2)%a(2)=u578_12(i1,i2)%a(2)+fcl(id6)*(cauxa-ce
     & ps_0)
      u578_12(i1,i2)%b(1)=u578_12(i1,i2)%b(1)+fcl(id6)*(cauxb-ce
     & ps_2)
      u578_12(i1,i2)%c(2)=u578_12(i1,i2)%c(2)+fcl(id6)*(-cauxc+c
     & eps_1)
      u578_12(i1,i2)%d(1)=u578_12(i1,i2)%d(1)+fcl(id6)*cf12(i1,i
     & 2)%ek0
      end do
      end do
  
* quqd -- p=p578,q=p612
      quqd=p578(0)*p612(0)-p578(1)*p612(1)-p578(2)*p612(2)-p578(
     & 3)*p612(3)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p578,qd=p612,v=cz34(i1,i2)%e,a=u578_34(i1,i2)%a,b=u578_34(i1
* ,i2)%b,c=u578_34(i1,i2)%c,d=u578_34(i1,i2)%d,cl=zcl(id6),nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p578(2)*p612(3)-p612(2)*p578(3))+
     & p578k0*(cz34(i1,i2)%e(2)*p612(3)-p612(2)*cz34(i1,i2)%e(3)
     & )-p612k0*(cz34(i1,i2)%e(2)*p578(3)-p578(2)*cz34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p578k0+p578(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i1,i2)%e(3)*p612k0+p612(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p578(0)-cz34(i1,i2)%e(1)*p578(1)-cz3
     & 4(i1,i2)%e(2)*p578(2)-cz34(i1,i2)%e(3)*p578(3)
      cvqd=cz34(i1,i2)%e(0)*p612(0)-cz34(i1,i2)%e(1)*p612(1)-cz3
     & 4(i1,i2)%e(2)*p612(2)-cz34(i1,i2)%e(3)*p612(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p578k0*cvqd+p612k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p612(2)+p612k0*cz34(i1,i2)%e(2)
      cauxc=+cz34(i1,i2)%ek0*p578(2)-p578k0*cz34(i1,i2)%e(2)
      u578_34(i1,i2)%a(2)=zcl(id6)*(cauxa-ceps_0)
      u578_34(i1,i2)%b(1)=zcl(id6)*(cauxb-ceps_2)
      u578_34(i1,i2)%c(2)=zcl(id6)*(-cauxc+ceps_1)
      u578_34(i1,i2)%d(1)=zcl(id6)*cz34(i1,i2)%ek0
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p578,qd=p612,v=cf34(i1,i2)%e,a=u578_34(i1,i2)%a,b=u578_34(i1
* ,i2)%b,c=u578_34(i1,i2)%c,d=u578_34(i1,i2)%d,cl=fcl(id6),nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p578(2)*p612(3)-p612(2)*p578(3))+
     & p578k0*(cf34(i1,i2)%e(2)*p612(3)-p612(2)*cf34(i1,i2)%e(3)
     & )-p612k0*(cf34(i1,i2)%e(2)*p578(3)-p578(2)*cf34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p578k0+p578(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i1,i2)%e(3)*p612k0+p612(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p578(0)-cf34(i1,i2)%e(1)*p578(1)-cf3
     & 4(i1,i2)%e(2)*p578(2)-cf34(i1,i2)%e(3)*p578(3)
      cvqd=cf34(i1,i2)%e(0)*p612(0)-cf34(i1,i2)%e(1)*p612(1)-cf3
     & 4(i1,i2)%e(2)*p612(2)-cf34(i1,i2)%e(3)*p612(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p578k0*cvqd+p612k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p612(2)+p612k0*cf34(i1,i2)%e(2)
      cauxc=+cf34(i1,i2)%ek0*p578(2)-p578k0*cf34(i1,i2)%e(2)
      u578_34(i1,i2)%a(2)=u578_34(i1,i2)%a(2)+fcl(id6)*(cauxa-ce
     & ps_0)
      u578_34(i1,i2)%b(1)=u578_34(i1,i2)%b(1)+fcl(id6)*(cauxb-ce
     & ps_2)
      u578_34(i1,i2)%c(2)=u578_34(i1,i2)%c(2)+fcl(id6)*(-cauxc+c
     & eps_1)
      u578_34(i1,i2)%d(1)=u578_34(i1,i2)%d(1)+fcl(id6)*cf34(i1,i
     & 2)%ek0
      end do
      end do
  
* quqd -- p=p712,q=p856
      quqd=p712(0)*p856(0)-p712(1)*p856(1)-p712(2)*p856(2)-p712(
     & 3)*p856(3)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p712,qd=p856,v=cz34(i1,i2)%e,a=u712_34(i1,i2)%a,b=u712_34(i1
* ,i2)%b,c=u712_34(i1,i2)%c,d=u712_34(i1,i2)%d,cl=zcl(id7),nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p712(2)*p856(3)-p856(2)*p712(3))+
     & p712k0*(cz34(i1,i2)%e(2)*p856(3)-p856(2)*cz34(i1,i2)%e(3)
     & )-p856k0*(cz34(i1,i2)%e(2)*p712(3)-p712(2)*cz34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p712k0+p712(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i1,i2)%e(3)*p856k0+p856(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p712(0)-cz34(i1,i2)%e(1)*p712(1)-cz3
     & 4(i1,i2)%e(2)*p712(2)-cz34(i1,i2)%e(3)*p712(3)
      cvqd=cz34(i1,i2)%e(0)*p856(0)-cz34(i1,i2)%e(1)*p856(1)-cz3
     & 4(i1,i2)%e(2)*p856(2)-cz34(i1,i2)%e(3)*p856(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p712k0*cvqd+p856k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p856(2)+p856k0*cz34(i1,i2)%e(2)
      cauxc=+cz34(i1,i2)%ek0*p712(2)-p712k0*cz34(i1,i2)%e(2)
      u712_34(i1,i2)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u712_34(i1,i2)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u712_34(i1,i2)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u712_34(i1,i2)%d(1)=zcl(id7)*cz34(i1,i2)%ek0
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p712,qd=p856,v=cf34(i1,i2)%e,a=u712_34(i1,i2)%a,b=u712_34(i1
* ,i2)%b,c=u712_34(i1,i2)%c,d=u712_34(i1,i2)%d,cl=fcl(id7),nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p712(2)*p856(3)-p856(2)*p712(3))+
     & p712k0*(cf34(i1,i2)%e(2)*p856(3)-p856(2)*cf34(i1,i2)%e(3)
     & )-p856k0*(cf34(i1,i2)%e(2)*p712(3)-p712(2)*cf34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p712k0+p712(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i1,i2)%e(3)*p856k0+p856(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p712(0)-cf34(i1,i2)%e(1)*p712(1)-cf3
     & 4(i1,i2)%e(2)*p712(2)-cf34(i1,i2)%e(3)*p712(3)
      cvqd=cf34(i1,i2)%e(0)*p856(0)-cf34(i1,i2)%e(1)*p856(1)-cf3
     & 4(i1,i2)%e(2)*p856(2)-cf34(i1,i2)%e(3)*p856(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p712k0*cvqd+p856k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p856(2)+p856k0*cf34(i1,i2)%e(2)
      cauxc=+cf34(i1,i2)%ek0*p712(2)-p712k0*cf34(i1,i2)%e(2)
      u712_34(i1,i2)%a(2)=u712_34(i1,i2)%a(2)+fcl(id7)*(cauxa-ce
     & ps_0)
      u712_34(i1,i2)%b(1)=u712_34(i1,i2)%b(1)+fcl(id7)*(cauxb-ce
     & ps_2)
      u712_34(i1,i2)%c(2)=u712_34(i1,i2)%c(2)+fcl(id7)*(-cauxc+c
     & eps_1)
      u712_34(i1,i2)%d(1)=u712_34(i1,i2)%d(1)+fcl(id7)*cf34(i1,i
     & 2)%ek0
      end do
      end do
  
* quqd -- p=p734,q=p856
      quqd=p734(0)*p856(0)-p734(1)*p856(1)-p734(2)*p856(2)-p734(
     & 3)*p856(3)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p734,qd=p856,v=cz12(i1,i2)%e,a=u734_12(i1,i2)%a,b=u734_12(i1
* ,i2)%b,c=u734_12(i1,i2)%c,d=u734_12(i1,i2)%d,cl=zcl(id7),nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p734(2)*p856(3)-p856(2)*p734(3))+
     & p734k0*(cz12(i1,i2)%e(2)*p856(3)-p856(2)*cz12(i1,i2)%e(3)
     & )-p856k0*(cz12(i1,i2)%e(2)*p734(3)-p734(2)*cz12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p734k0+p734(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i1,i2)%e(3)*p856k0+p856(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p734(0)-cz12(i1,i2)%e(1)*p734(1)-cz1
     & 2(i1,i2)%e(2)*p734(2)-cz12(i1,i2)%e(3)*p734(3)
      cvqd=cz12(i1,i2)%e(0)*p856(0)-cz12(i1,i2)%e(1)*p856(1)-cz1
     & 2(i1,i2)%e(2)*p856(2)-cz12(i1,i2)%e(3)*p856(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p734k0*cvqd+p856k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p856(2)+p856k0*cz12(i1,i2)%e(2)
      cauxc=+cz12(i1,i2)%ek0*p734(2)-p734k0*cz12(i1,i2)%e(2)
      u734_12(i1,i2)%a(2)=zcl(id7)*(cauxa-ceps_0)
      u734_12(i1,i2)%b(1)=zcl(id7)*(cauxb-ceps_2)
      u734_12(i1,i2)%c(2)=zcl(id7)*(-cauxc+ceps_1)
      u734_12(i1,i2)%d(1)=zcl(id7)*cz12(i1,i2)%ek0
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p734,qd=p856,v=cf12(i1,i2)%e,a=u734_12(i1,i2)%a,b=u734_12(i1
* ,i2)%b,c=u734_12(i1,i2)%c,d=u734_12(i1,i2)%d,cl=fcl(id7),nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p734(2)*p856(3)-p856(2)*p734(3))+
     & p734k0*(cf12(i1,i2)%e(2)*p856(3)-p856(2)*cf12(i1,i2)%e(3)
     & )-p856k0*(cf12(i1,i2)%e(2)*p734(3)-p734(2)*cf12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p734k0+p734(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i1,i2)%e(3)*p856k0+p856(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p734(0)-cf12(i1,i2)%e(1)*p734(1)-cf1
     & 2(i1,i2)%e(2)*p734(2)-cf12(i1,i2)%e(3)*p734(3)
      cvqd=cf12(i1,i2)%e(0)*p856(0)-cf12(i1,i2)%e(1)*p856(1)-cf1
     & 2(i1,i2)%e(2)*p856(2)-cf12(i1,i2)%e(3)*p856(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p734k0*cvqd+p856k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p856(2)+p856k0*cf12(i1,i2)%e(2)
      cauxc=+cf12(i1,i2)%ek0*p734(2)-p734k0*cf12(i1,i2)%e(2)
      u734_12(i1,i2)%a(2)=u734_12(i1,i2)%a(2)+fcl(id7)*(cauxa-ce
     & ps_0)
      u734_12(i1,i2)%b(1)=u734_12(i1,i2)%b(1)+fcl(id7)*(cauxb-ce
     & ps_2)
      u734_12(i1,i2)%c(2)=u734_12(i1,i2)%c(2)+fcl(id7)*(-cauxc+c
     & eps_1)
      u734_12(i1,i2)%d(1)=u734_12(i1,i2)%d(1)+fcl(id7)*cf12(i1,i
     & 2)%ek0
      end do
      end do
  
* quqd -- p=p756,q=p834
      quqd=p756(0)*p834(0)-p756(1)*p834(1)-p756(2)*p834(2)-p756(
     & 3)*p834(3)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p756,qd=p834,v=cz12(i1,i2)%e,a=u756_12(i1,i2)%a,b=u756_12(i1
* ,i2)%b,c=u756_12(i1,i2)%c,d=u756_12(i1,i2)%d,cl=zcl(id8),nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p756(2)*p834(3)-p834(2)*p756(3))+
     & p756k0*(cz12(i1,i2)%e(2)*p834(3)-p834(2)*cz12(i1,i2)%e(3)
     & )-p834k0*(cz12(i1,i2)%e(2)*p756(3)-p756(2)*cz12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p756k0+p756(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i1,i2)%e(3)*p834k0+p834(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p756(0)-cz12(i1,i2)%e(1)*p756(1)-cz1
     & 2(i1,i2)%e(2)*p756(2)-cz12(i1,i2)%e(3)*p756(3)
      cvqd=cz12(i1,i2)%e(0)*p834(0)-cz12(i1,i2)%e(1)*p834(1)-cz1
     & 2(i1,i2)%e(2)*p834(2)-cz12(i1,i2)%e(3)*p834(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p756k0*cvqd+p834k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p834(2)+p834k0*cz12(i1,i2)%e(2)
      cauxc=+cz12(i1,i2)%ek0*p756(2)-p756k0*cz12(i1,i2)%e(2)
      u756_12(i1,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      u756_12(i1,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      u756_12(i1,i2)%c(2)=zcl(id8)*(-cauxc+ceps_1)
      u756_12(i1,i2)%d(1)=zcl(id8)*cz12(i1,i2)%ek0
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p756,qd=p834,v=cf12(i1,i2)%e,a=u756_12(i1,i2)%a,b=u756_12(i1
* ,i2)%b,c=u756_12(i1,i2)%c,d=u756_12(i1,i2)%d,cl=fcl(id8),nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p756(2)*p834(3)-p834(2)*p756(3))+
     & p756k0*(cf12(i1,i2)%e(2)*p834(3)-p834(2)*cf12(i1,i2)%e(3)
     & )-p834k0*(cf12(i1,i2)%e(2)*p756(3)-p756(2)*cf12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p756k0+p756(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i1,i2)%e(3)*p834k0+p834(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p756(0)-cf12(i1,i2)%e(1)*p756(1)-cf1
     & 2(i1,i2)%e(2)*p756(2)-cf12(i1,i2)%e(3)*p756(3)
      cvqd=cf12(i1,i2)%e(0)*p834(0)-cf12(i1,i2)%e(1)*p834(1)-cf1
     & 2(i1,i2)%e(2)*p834(2)-cf12(i1,i2)%e(3)*p834(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p756k0*cvqd+p834k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p834(2)+p834k0*cf12(i1,i2)%e(2)
      cauxc=+cf12(i1,i2)%ek0*p756(2)-p756k0*cf12(i1,i2)%e(2)
      u756_12(i1,i2)%a(2)=u756_12(i1,i2)%a(2)+fcl(id8)*(cauxa-ce
     & ps_0)
      u756_12(i1,i2)%b(1)=u756_12(i1,i2)%b(1)+fcl(id8)*(cauxb-ce
     & ps_2)
      u756_12(i1,i2)%c(2)=u756_12(i1,i2)%c(2)+fcl(id8)*(-cauxc+c
     & eps_1)
      u756_12(i1,i2)%d(1)=u756_12(i1,i2)%d(1)+fcl(id8)*cf12(i1,i
     & 2)%ek0
      end do
      end do
  
* quqd -- p=p756,q=p812
      quqd=p756(0)*p812(0)-p756(1)*p812(1)-p756(2)*p812(2)-p756(
     & 3)*p812(3)
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p756,qd=p812,v=cz34(i1,i2)%e,a=u756_34(i1,i2)%a,b=u756_34(i1
* ,i2)%b,c=u756_34(i1,i2)%c,d=u756_34(i1,i2)%d,cl=zcl(id8),nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p756(2)*p812(3)-p812(2)*p756(3))+
     & p756k0*(cz34(i1,i2)%e(2)*p812(3)-p812(2)*cz34(i1,i2)%e(3)
     & )-p812k0*(cz34(i1,i2)%e(2)*p756(3)-p756(2)*cz34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p756k0+p756(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i1,i2)%e(3)*p812k0+p812(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p756(0)-cz34(i1,i2)%e(1)*p756(1)-cz3
     & 4(i1,i2)%e(2)*p756(2)-cz34(i1,i2)%e(3)*p756(3)
      cvqd=cz34(i1,i2)%e(0)*p812(0)-cz34(i1,i2)%e(1)*p812(1)-cz3
     & 4(i1,i2)%e(2)*p812(2)-cz34(i1,i2)%e(3)*p812(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p756k0*cvqd+p812k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p812(2)+p812k0*cz34(i1,i2)%e(2)
      cauxc=+cz34(i1,i2)%ek0*p756(2)-p756k0*cz34(i1,i2)%e(2)
      u756_34(i1,i2)%a(2)=zcl(id8)*(cauxa-ceps_0)
      u756_34(i1,i2)%b(1)=zcl(id8)*(cauxb-ceps_2)
      u756_34(i1,i2)%c(2)=zcl(id8)*(-cauxc+ceps_1)
      u756_34(i1,i2)%d(1)=zcl(id8)*cz34(i1,i2)%ek0
      end do
      end do
  
      do i1=1,2
      do i2=1,2
* TW0 -- qu=p756,qd=p812,v=cf34(i1,i2)%e,a=u756_34(i1,i2)%a,b=u756_34(i1
* ,i2)%b,c=u756_34(i1,i2)%c,d=u756_34(i1,i2)%d,cl=fcl(id8),nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p756(2)*p812(3)-p812(2)*p756(3))+
     & p756k0*(cf34(i1,i2)%e(2)*p812(3)-p812(2)*cf34(i1,i2)%e(3)
     & )-p812k0*(cf34(i1,i2)%e(2)*p756(3)-p756(2)*cf34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p756k0+p756(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i1,i2)%e(3)*p812k0+p812(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p756(0)-cf34(i1,i2)%e(1)*p756(1)-cf3
     & 4(i1,i2)%e(2)*p756(2)-cf34(i1,i2)%e(3)*p756(3)
      cvqd=cf34(i1,i2)%e(0)*p812(0)-cf34(i1,i2)%e(1)*p812(1)-cf3
     & 4(i1,i2)%e(2)*p812(2)-cf34(i1,i2)%e(3)*p812(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p756k0*cvqd+p812k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p812(2)+p812k0*cf34(i1,i2)%e(2)
      cauxc=+cf34(i1,i2)%ek0*p756(2)-p756k0*cf34(i1,i2)%e(2)
      u756_34(i1,i2)%a(2)=u756_34(i1,i2)%a(2)+fcl(id8)*(cauxa-ce
     & ps_0)
      u756_34(i1,i2)%b(1)=u756_34(i1,i2)%b(1)+fcl(id8)*(cauxb-ce
     & ps_2)
      u756_34(i1,i2)%c(2)=u756_34(i1,i2)%c(2)+fcl(id8)*(-cauxc+c
     & eps_1)
      u756_34(i1,i2)%d(1)=u756_34(i1,i2)%d(1)+fcl(id8)*cf34(i1,i
     & 2)%ek0
      end do
      end do
  
  
* compute all single W insertions in the middle of a triple insertion   
*  Zline                                                                
*                                                                       
*            |                                                          
*            |_W__/                                                     
*            |    \                                                     
*            |                                                          
  
      if (iup(id1).eq.1) then ! first W- and then W+
* quqd -- p=p134,q=p256
      quqd=p134(0)*p256(0)-p134(1)*p256(1)-p134(2)*p256(2)-p134(
     & 3)*p256(3)
* TW -- qu=p134,qd=p256,v=cw78%e,a=u134_78%a,b=u134_78%b,c=u134_78%c,d=u
* 134_78%d,cl=wcl,nsum=0
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
      u134_78%a(2,2)=wcl*(cauxa-ceps_0)
      u134_78%b(1,2)=wcl*(cauxb-ceps_2)
      u134_78%c(2,1)=wcl*(-cauxc+ceps_1)
      u134_78%d(1,1)=wcl*cw78%ek0
* quqd -- p=p178,q=p234
      quqd=p178(0)*p234(0)-p178(1)*p234(1)-p178(2)*p234(2)-p178(
     & 3)*p234(3)
* TW -- qu=p178,qd=p234,v=cw56%e,a=u178_56%a,b=u178_56%b,c=u178_56%c,d=u
* 178_56%d,cl=wcl,nsum=0
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
      u178_56%a(2,2)=wcl*(cauxa-ceps_0)
      u178_56%b(1,2)=wcl*(cauxb-ceps_2)
      u178_56%c(2,1)=wcl*(-cauxc+ceps_1)
      u178_56%d(1,1)=wcl*cw56%ek0
      else   !first W+ and then W-
* quqd -- p=p134,q=p278
      quqd=p134(0)*p278(0)-p134(1)*p278(1)-p134(2)*p278(2)-p134(
     & 3)*p278(3)
* TW -- qu=p134,qd=p278,v=cw56%e,a=u134_56%a,b=u134_56%b,c=u134_56%c,d=u
* 134_56%d,cl=wcl,nsum=0
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
      u134_56%a(2,2)=wcl*(cauxa-ceps_0)
      u134_56%b(1,2)=wcl*(cauxb-ceps_2)
      u134_56%c(2,1)=wcl*(-cauxc+ceps_1)
      u134_56%d(1,1)=wcl*cw56%ek0
* quqd -- p=p156,q=p234
      quqd=p156(0)*p234(0)-p156(1)*p234(1)-p156(2)*p234(2)-p156(
     & 3)*p234(3)
* TW -- qu=p156,qd=p234,v=cw78%e,a=u156_78%a,b=u156_78%b,c=u156_78%c,d=u
* 156_78%d,cl=wcl,nsum=0
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
      u156_78%a(2,2)=wcl*(cauxa-ceps_0)
      u156_78%b(1,2)=wcl*(cauxb-ceps_2)
      u156_78%c(2,1)=wcl*(-cauxc+ceps_1)
      u156_78%d(1,1)=wcl*cw78%ek0
      endif
      if (iup(id3).eq.1) then ! first W- and then W+
* quqd -- p=p312,q=p456
      quqd=p312(0)*p456(0)-p312(1)*p456(1)-p312(2)*p456(2)-p312(
     & 3)*p456(3)
* TW -- qu=p312,qd=p456,v=cw78%e,a=u312_78%a,b=u312_78%b,c=u312_78%c,d=u
* 312_78%d,cl=wcl,nsum=0
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
      u312_78%a(2,2)=wcl*(cauxa-ceps_0)
      u312_78%b(1,2)=wcl*(cauxb-ceps_2)
      u312_78%c(2,1)=wcl*(-cauxc+ceps_1)
      u312_78%d(1,1)=wcl*cw78%ek0
* quqd -- p=p378,q=p412
      quqd=p378(0)*p412(0)-p378(1)*p412(1)-p378(2)*p412(2)-p378(
     & 3)*p412(3)
* TW -- qu=p378,qd=p412,v=cw56%e,a=u378_56%a,b=u378_56%b,c=u378_56%c,d=u
* 378_56%d,cl=wcl,nsum=0
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
      u378_56%a(2,2)=wcl*(cauxa-ceps_0)
      u378_56%b(1,2)=wcl*(cauxb-ceps_2)
      u378_56%c(2,1)=wcl*(-cauxc+ceps_1)
      u378_56%d(1,1)=wcl*cw56%ek0
      else   !first W+ and then W-
* quqd -- p=p312,q=p478
      quqd=p312(0)*p478(0)-p312(1)*p478(1)-p312(2)*p478(2)-p312(
     & 3)*p478(3)
* TW -- qu=p312,qd=p478,v=cw56%e,a=u312_56%a,b=u312_56%b,c=u312_56%c,d=u
* 312_56%d,cl=wcl,nsum=0
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
      u312_56%a(2,2)=wcl*(cauxa-ceps_0)
      u312_56%b(1,2)=wcl*(cauxb-ceps_2)
      u312_56%c(2,1)=wcl*(-cauxc+ceps_1)
      u312_56%d(1,1)=wcl*cw56%ek0
* quqd -- p=p356,q=p412
      quqd=p356(0)*p412(0)-p356(1)*p412(1)-p356(2)*p412(2)-p356(
     & 3)*p412(3)
* TW -- qu=p356,qd=p412,v=cw78%e,a=u356_78%a,b=u356_78%b,c=u356_78%c,d=u
* 356_78%d,cl=wcl,nsum=0
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
      u356_78%a(2,2)=wcl*(cauxa-ceps_0)
      u356_78%b(1,2)=wcl*(cauxb-ceps_2)
      u356_78%c(2,1)=wcl*(-cauxc+ceps_1)
      u356_78%d(1,1)=wcl*cw78%ek0
      endif
  
* compute all single Z or gamma or higgs insertions in the middle of a  
*  triple insertion Zline                                               
*                                                                       
*            |                                                          
*            |_Z,f,h__/                                                 
*            |        \                                                 
*            |                                                          
  
      if (iup(id1).eq.1) then  ! first W- and then W+
  
* quqd -- p=p178,q=p256
      quqd=p178(0)*p256(0)-p178(1)*p256(1)-p178(2)*p256(2)-p178(
     & 3)*p256(3)
      do i1=1,2
      do i2=1,2
* T -- qu=p178,qd=p256,v=cz34(i1,i2)%e,a=u178_34(i1,i2)%a,b=u178_34(i1,i
* 2)%b,c=u178_34(i1,i2)%c,d=u178_34(i1,i2)%d,cr=zcr(id1-1),cl=zcl(id1-1)
* ,nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p178(2)*p256(3)-p256(2)*p178(3))+
     & p178k0*(cz34(i1,i2)%e(2)*p256(3)-p256(2)*cz34(i1,i2)%e(3)
     & )-p256k0*(cz34(i1,i2)%e(2)*p178(3)-p178(2)*cz34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p178k0+p178(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i1,i2)%e(3)*p256k0+p256(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p178(0)-cz34(i1,i2)%e(1)*p178(1)-cz3
     & 4(i1,i2)%e(2)*p178(2)-cz34(i1,i2)%e(3)*p178(3)
      cvqd=cz34(i1,i2)%e(0)*p256(0)-cz34(i1,i2)%e(1)*p256(1)-cz3
     & 4(i1,i2)%e(2)*p256(2)-cz34(i1,i2)%e(3)*p256(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p178k0*cvqd+p256k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p256(2)+p256k0*cz34(i1,i2)%e(2)
      cauxc=+cz34(i1,i2)%ek0*p178(2)-p178k0*cz34(i1,i2)%e(2)
      u178_34(i1,i2)%a(1,1)=zcr(id1-1)*(cauxa+ceps_0)
      u178_34(i1,i2)%a(2,2)=zcl(id1-1)*(cauxa-ceps_0)
      u178_34(i1,i2)%b(1,2)=zcl(id1-1)*(cauxb-ceps_2)
      u178_34(i1,i2)%b(2,1)=zcr(id1-1)*(-cauxb-ceps_2)
      u178_34(i1,i2)%c(1,2)=zcr(id1-1)*(cauxc+ceps_1)
      u178_34(i1,i2)%c(2,1)=zcl(id1-1)*(-cauxc+ceps_1)
      u178_34(i1,i2)%d(1,1)=zcl(id1-1)*cz34(i1,i2)%ek0
      u178_34(i1,i2)%d(2,2)=zcr(id1-1)*cz34(i1,i2)%ek0
      end do
      end do
      do i1=1,2
      do i2=1,2
* T -- qu=p178,qd=p256,v=cf34(i1,i2)%e,a=u178_34(i1,i2)%a,b=u178_34(i1,i
* 2)%b,c=u178_34(i1,i2)%c,d=u178_34(i1,i2)%d,cr=fcr(id1-1),cl=fcl(id1-1)
* ,nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p178(2)*p256(3)-p256(2)*p178(3))+
     & p178k0*(cf34(i1,i2)%e(2)*p256(3)-p256(2)*cf34(i1,i2)%e(3)
     & )-p256k0*(cf34(i1,i2)%e(2)*p178(3)-p178(2)*cf34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p178k0+p178(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i1,i2)%e(3)*p256k0+p256(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p178(0)-cf34(i1,i2)%e(1)*p178(1)-cf3
     & 4(i1,i2)%e(2)*p178(2)-cf34(i1,i2)%e(3)*p178(3)
      cvqd=cf34(i1,i2)%e(0)*p256(0)-cf34(i1,i2)%e(1)*p256(1)-cf3
     & 4(i1,i2)%e(2)*p256(2)-cf34(i1,i2)%e(3)*p256(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p178k0*cvqd+p256k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p256(2)+p256k0*cf34(i1,i2)%e(2)
      cauxc=+cf34(i1,i2)%ek0*p178(2)-p178k0*cf34(i1,i2)%e(2)
      u178_34(i1,i2)%a(1,1)=u178_34(i1,i2)%a(1,1)+fcr(id1-1)*(ca
     & uxa+ceps_0)
      u178_34(i1,i2)%a(2,2)=u178_34(i1,i2)%a(2,2)+fcl(id1-1)*(ca
     & uxa-ceps_0)
      u178_34(i1,i2)%b(1,2)=u178_34(i1,i2)%b(1,2)+fcl(id1-1)*(ca
     & uxb-ceps_2)
      u178_34(i1,i2)%b(2,1)=u178_34(i1,i2)%b(2,1)+fcr(id1-1)*(-c
     & auxb-ceps_2)
      u178_34(i1,i2)%c(1,2)=u178_34(i1,i2)%c(1,2)+fcr(id1-1)*(ca
     & uxc+ceps_1)
      u178_34(i1,i2)%c(2,1)=u178_34(i1,i2)%c(2,1)+fcl(id1-1)*(-c
     & auxc+ceps_1)
      u178_34(i1,i2)%d(1,1)=u178_34(i1,i2)%d(1,1)+fcl(id1-1)*cf3
     & 4(i1,i2)%ek0
      u178_34(i1,i2)%d(2,2)=u178_34(i1,i2)%d(2,2)+fcr(id1-1)*cf3
     & 4(i1,i2)%ek0
      end do
      end do
  
  
      else   !first W+ and then W-
  
* quqd -- p=p156,q=p278
      quqd=p156(0)*p278(0)-p156(1)*p278(1)-p156(2)*p278(2)-p156(
     & 3)*p278(3)
      do i1=1,2
      do i2=1,2
* T -- qu=p156,qd=p278,v=cz34(i1,i2)%e,a=u156_34(i1,i2)%a,b=u156_34(i1,i
* 2)%b,c=u156_34(i1,i2)%c,d=u156_34(i1,i2)%d,cr=zcr(id1+1),cl=zcl(id1+1)
* ,nsum=0
      ceps_0=-cz34(i1,i2)%ek0*(p156(2)*p278(3)-p278(2)*p156(3))+
     & p156k0*(cz34(i1,i2)%e(2)*p278(3)-p278(2)*cz34(i1,i2)%e(3)
     & )-p278k0*(cz34(i1,i2)%e(2)*p156(3)-p156(2)*cz34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz34(i1,i2)%e(3)*p156k0+p156(3)*cz34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz34(i1,i2)%e(3)*p278k0+p278(3)*cz34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz34(i1,i2)%e(0)*p156(0)-cz34(i1,i2)%e(1)*p156(1)-cz3
     & 4(i1,i2)%e(2)*p156(2)-cz34(i1,i2)%e(3)*p156(3)
      cvqd=cz34(i1,i2)%e(0)*p278(0)-cz34(i1,i2)%e(1)*p278(1)-cz3
     & 4(i1,i2)%e(2)*p278(2)-cz34(i1,i2)%e(3)*p278(3)
      cauxa=-cz34(i1,i2)%ek0*quqd+p156k0*cvqd+p278k0*cvqu
      cauxb=-cz34(i1,i2)%ek0*p278(2)+p278k0*cz34(i1,i2)%e(2)
      cauxc=+cz34(i1,i2)%ek0*p156(2)-p156k0*cz34(i1,i2)%e(2)
      u156_34(i1,i2)%a(1,1)=zcr(id1+1)*(cauxa+ceps_0)
      u156_34(i1,i2)%a(2,2)=zcl(id1+1)*(cauxa-ceps_0)
      u156_34(i1,i2)%b(1,2)=zcl(id1+1)*(cauxb-ceps_2)
      u156_34(i1,i2)%b(2,1)=zcr(id1+1)*(-cauxb-ceps_2)
      u156_34(i1,i2)%c(1,2)=zcr(id1+1)*(cauxc+ceps_1)
      u156_34(i1,i2)%c(2,1)=zcl(id1+1)*(-cauxc+ceps_1)
      u156_34(i1,i2)%d(1,1)=zcl(id1+1)*cz34(i1,i2)%ek0
      u156_34(i1,i2)%d(2,2)=zcr(id1+1)*cz34(i1,i2)%ek0
      end do
      end do
      do i1=1,2
      do i2=1,2
* T -- qu=p156,qd=p278,v=cf34(i1,i2)%e,a=u156_34(i1,i2)%a,b=u156_34(i1,i
* 2)%b,c=u156_34(i1,i2)%c,d=u156_34(i1,i2)%d,cr=fcr(id1+1),cl=fcl(id1+1)
* ,nsum=1
      ceps_0=-cf34(i1,i2)%ek0*(p156(2)*p278(3)-p278(2)*p156(3))+
     & p156k0*(cf34(i1,i2)%e(2)*p278(3)-p278(2)*cf34(i1,i2)%e(3)
     & )-p278k0*(cf34(i1,i2)%e(2)*p156(3)-p156(2)*cf34(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf34(i1,i2)%e(3)*p156k0+p156(3)*cf34(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf34(i1,i2)%e(3)*p278k0+p278(3)*cf34(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf34(i1,i2)%e(0)*p156(0)-cf34(i1,i2)%e(1)*p156(1)-cf3
     & 4(i1,i2)%e(2)*p156(2)-cf34(i1,i2)%e(3)*p156(3)
      cvqd=cf34(i1,i2)%e(0)*p278(0)-cf34(i1,i2)%e(1)*p278(1)-cf3
     & 4(i1,i2)%e(2)*p278(2)-cf34(i1,i2)%e(3)*p278(3)
      cauxa=-cf34(i1,i2)%ek0*quqd+p156k0*cvqd+p278k0*cvqu
      cauxb=-cf34(i1,i2)%ek0*p278(2)+p278k0*cf34(i1,i2)%e(2)
      cauxc=+cf34(i1,i2)%ek0*p156(2)-p156k0*cf34(i1,i2)%e(2)
      u156_34(i1,i2)%a(1,1)=u156_34(i1,i2)%a(1,1)+fcr(id1+1)*(ca
     & uxa+ceps_0)
      u156_34(i1,i2)%a(2,2)=u156_34(i1,i2)%a(2,2)+fcl(id1+1)*(ca
     & uxa-ceps_0)
      u156_34(i1,i2)%b(1,2)=u156_34(i1,i2)%b(1,2)+fcl(id1+1)*(ca
     & uxb-ceps_2)
      u156_34(i1,i2)%b(2,1)=u156_34(i1,i2)%b(2,1)+fcr(id1+1)*(-c
     & auxb-ceps_2)
      u156_34(i1,i2)%c(1,2)=u156_34(i1,i2)%c(1,2)+fcr(id1+1)*(ca
     & uxc+ceps_1)
      u156_34(i1,i2)%c(2,1)=u156_34(i1,i2)%c(2,1)+fcl(id1+1)*(-c
     & auxc+ceps_1)
      u156_34(i1,i2)%d(1,1)=u156_34(i1,i2)%d(1,1)+fcl(id1+1)*cf3
     & 4(i1,i2)%ek0
      u156_34(i1,i2)%d(2,2)=u156_34(i1,i2)%d(2,2)+fcr(id1+1)*cf3
     & 4(i1,i2)%ek0
      end do
      end do
  
* rmh < 0 in our convention means no higgs coupling                     
  
      if (id1.eq.5.and.id3.eq.5.and.rmh.ge.0.d0) then
  
* TH -- qu=p156,qd=p278,a=clineth%a,b=clineth%b,c=clineth%c,coupl=1.d0
      auxa=-p156k0*p278(2)+p278k0*p156(2)
      cauxa=auxa-cim*(p278(3)*p156k0-p156(3)*p278k0)
      clineth%a(1,2)=cauxa
      clineth%a(2,1)=-conjg(cauxa)
      clineth%b(1,1)=p278k0
      clineth%b(2,2)=clineth%b(1,1)
      clineth%c(1,1)=p156k0
      clineth%c(2,2)=clineth%c(1,1)
  
      cfactor= rhbb*rhtt/(s34-cmh2)+rhhbb*rhhtt/(s34-cmhh2)
         do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  u156_34(i1,i2)%a(i3,i4)=ch34(i1,i2)
     &                 *clineth%a(i3,i4)
     &                 *cfactor
                else
                  u156_34(i1,i2)%b(i3,i4)=ch34(i1,i2)
     &                 *clineth%b(i3,i4)
     &                 *cfactor
                  u156_34(i1,i2)%c(i3,i4)=ch34(i1,i2)
     &                 *clineth%c(i3,i4)
     &                 *cfactor
                endif
              enddo
            enddo
          enddo
        enddo
      else
         do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  u156_34(i1,i2)%a(i3,i4)=czero
                else
                  u156_34(i1,i2)%b(i3,i4)=czero
                  u156_34(i1,i2)%c(i3,i4)=czero
                endif
              enddo
            enddo
          enddo
        enddo
      endif
  
  
  
      endif
      if (iup(id3).eq.1) then  ! first W- and then W+
  
* quqd -- p=p378,q=p456
      quqd=p378(0)*p456(0)-p378(1)*p456(1)-p378(2)*p456(2)-p378(
     & 3)*p456(3)
      do i1=1,2
      do i2=1,2
* T -- qu=p378,qd=p456,v=cz12(i1,i2)%e,a=u378_12(i1,i2)%a,b=u378_12(i1,i
* 2)%b,c=u378_12(i1,i2)%c,d=u378_12(i1,i2)%d,cr=zcr(id3-1),cl=zcl(id3-1)
* ,nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p378(2)*p456(3)-p456(2)*p378(3))+
     & p378k0*(cz12(i1,i2)%e(2)*p456(3)-p456(2)*cz12(i1,i2)%e(3)
     & )-p456k0*(cz12(i1,i2)%e(2)*p378(3)-p378(2)*cz12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p378k0+p378(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i1,i2)%e(3)*p456k0+p456(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p378(0)-cz12(i1,i2)%e(1)*p378(1)-cz1
     & 2(i1,i2)%e(2)*p378(2)-cz12(i1,i2)%e(3)*p378(3)
      cvqd=cz12(i1,i2)%e(0)*p456(0)-cz12(i1,i2)%e(1)*p456(1)-cz1
     & 2(i1,i2)%e(2)*p456(2)-cz12(i1,i2)%e(3)*p456(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p378k0*cvqd+p456k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p456(2)+p456k0*cz12(i1,i2)%e(2)
      cauxc=+cz12(i1,i2)%ek0*p378(2)-p378k0*cz12(i1,i2)%e(2)
      u378_12(i1,i2)%a(1,1)=zcr(id3-1)*(cauxa+ceps_0)
      u378_12(i1,i2)%a(2,2)=zcl(id3-1)*(cauxa-ceps_0)
      u378_12(i1,i2)%b(1,2)=zcl(id3-1)*(cauxb-ceps_2)
      u378_12(i1,i2)%b(2,1)=zcr(id3-1)*(-cauxb-ceps_2)
      u378_12(i1,i2)%c(1,2)=zcr(id3-1)*(cauxc+ceps_1)
      u378_12(i1,i2)%c(2,1)=zcl(id3-1)*(-cauxc+ceps_1)
      u378_12(i1,i2)%d(1,1)=zcl(id3-1)*cz12(i1,i2)%ek0
      u378_12(i1,i2)%d(2,2)=zcr(id3-1)*cz12(i1,i2)%ek0
      end do
      end do
      do i1=1,2
      do i2=1,2
* T -- qu=p378,qd=p456,v=cf12(i1,i2)%e,a=u378_12(i1,i2)%a,b=u378_12(i1,i
* 2)%b,c=u378_12(i1,i2)%c,d=u378_12(i1,i2)%d,cr=fcr(id3-1),cl=fcl(id3-1)
* ,nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p378(2)*p456(3)-p456(2)*p378(3))+
     & p378k0*(cf12(i1,i2)%e(2)*p456(3)-p456(2)*cf12(i1,i2)%e(3)
     & )-p456k0*(cf12(i1,i2)%e(2)*p378(3)-p378(2)*cf12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p378k0+p378(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i1,i2)%e(3)*p456k0+p456(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p378(0)-cf12(i1,i2)%e(1)*p378(1)-cf1
     & 2(i1,i2)%e(2)*p378(2)-cf12(i1,i2)%e(3)*p378(3)
      cvqd=cf12(i1,i2)%e(0)*p456(0)-cf12(i1,i2)%e(1)*p456(1)-cf1
     & 2(i1,i2)%e(2)*p456(2)-cf12(i1,i2)%e(3)*p456(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p378k0*cvqd+p456k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p456(2)+p456k0*cf12(i1,i2)%e(2)
      cauxc=+cf12(i1,i2)%ek0*p378(2)-p378k0*cf12(i1,i2)%e(2)
      u378_12(i1,i2)%a(1,1)=u378_12(i1,i2)%a(1,1)+fcr(id3-1)*(ca
     & uxa+ceps_0)
      u378_12(i1,i2)%a(2,2)=u378_12(i1,i2)%a(2,2)+fcl(id3-1)*(ca
     & uxa-ceps_0)
      u378_12(i1,i2)%b(1,2)=u378_12(i1,i2)%b(1,2)+fcl(id3-1)*(ca
     & uxb-ceps_2)
      u378_12(i1,i2)%b(2,1)=u378_12(i1,i2)%b(2,1)+fcr(id3-1)*(-c
     & auxb-ceps_2)
      u378_12(i1,i2)%c(1,2)=u378_12(i1,i2)%c(1,2)+fcr(id3-1)*(ca
     & uxc+ceps_1)
      u378_12(i1,i2)%c(2,1)=u378_12(i1,i2)%c(2,1)+fcl(id3-1)*(-c
     & auxc+ceps_1)
      u378_12(i1,i2)%d(1,1)=u378_12(i1,i2)%d(1,1)+fcl(id3-1)*cf1
     & 2(i1,i2)%ek0
      u378_12(i1,i2)%d(2,2)=u378_12(i1,i2)%d(2,2)+fcr(id3-1)*cf1
     & 2(i1,i2)%ek0
      end do
      end do
  
  
      else   !first W+ and then W-
  
* quqd -- p=p356,q=p478
      quqd=p356(0)*p478(0)-p356(1)*p478(1)-p356(2)*p478(2)-p356(
     & 3)*p478(3)
      do i1=1,2
      do i2=1,2
* T -- qu=p356,qd=p478,v=cz12(i1,i2)%e,a=u356_12(i1,i2)%a,b=u356_12(i1,i
* 2)%b,c=u356_12(i1,i2)%c,d=u356_12(i1,i2)%d,cr=zcr(id3+1),cl=zcl(id3+1)
* ,nsum=0
      ceps_0=-cz12(i1,i2)%ek0*(p356(2)*p478(3)-p478(2)*p356(3))+
     & p356k0*(cz12(i1,i2)%e(2)*p478(3)-p478(2)*cz12(i1,i2)%e(3)
     & )-p478k0*(cz12(i1,i2)%e(2)*p356(3)-p356(2)*cz12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cz12(i1,i2)%e(3)*p356k0+p356(3)*cz12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cz12(i1,i2)%e(3)*p478k0+p478(3)*cz12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cz12(i1,i2)%e(0)*p356(0)-cz12(i1,i2)%e(1)*p356(1)-cz1
     & 2(i1,i2)%e(2)*p356(2)-cz12(i1,i2)%e(3)*p356(3)
      cvqd=cz12(i1,i2)%e(0)*p478(0)-cz12(i1,i2)%e(1)*p478(1)-cz1
     & 2(i1,i2)%e(2)*p478(2)-cz12(i1,i2)%e(3)*p478(3)
      cauxa=-cz12(i1,i2)%ek0*quqd+p356k0*cvqd+p478k0*cvqu
      cauxb=-cz12(i1,i2)%ek0*p478(2)+p478k0*cz12(i1,i2)%e(2)
      cauxc=+cz12(i1,i2)%ek0*p356(2)-p356k0*cz12(i1,i2)%e(2)
      u356_12(i1,i2)%a(1,1)=zcr(id3+1)*(cauxa+ceps_0)
      u356_12(i1,i2)%a(2,2)=zcl(id3+1)*(cauxa-ceps_0)
      u356_12(i1,i2)%b(1,2)=zcl(id3+1)*(cauxb-ceps_2)
      u356_12(i1,i2)%b(2,1)=zcr(id3+1)*(-cauxb-ceps_2)
      u356_12(i1,i2)%c(1,2)=zcr(id3+1)*(cauxc+ceps_1)
      u356_12(i1,i2)%c(2,1)=zcl(id3+1)*(-cauxc+ceps_1)
      u356_12(i1,i2)%d(1,1)=zcl(id3+1)*cz12(i1,i2)%ek0
      u356_12(i1,i2)%d(2,2)=zcr(id3+1)*cz12(i1,i2)%ek0
      end do
      end do
      do i1=1,2
      do i2=1,2
* T -- qu=p356,qd=p478,v=cf12(i1,i2)%e,a=u356_12(i1,i2)%a,b=u356_12(i1,i
* 2)%b,c=u356_12(i1,i2)%c,d=u356_12(i1,i2)%d,cr=fcr(id3+1),cl=fcl(id3+1)
* ,nsum=1
      ceps_0=-cf12(i1,i2)%ek0*(p356(2)*p478(3)-p478(2)*p356(3))+
     & p356k0*(cf12(i1,i2)%e(2)*p478(3)-p478(2)*cf12(i1,i2)%e(3)
     & )-p478k0*(cf12(i1,i2)%e(2)*p356(3)-p356(2)*cf12(i1,i2)%e(
     & 3))
      ceps_0=ceps_0*cim
      ceps_1=-cf12(i1,i2)%e(3)*p356k0+p356(3)*cf12(i1,i2)%ek0
      ceps_1=ceps_1*cim
      ceps_2=-cf12(i1,i2)%e(3)*p478k0+p478(3)*cf12(i1,i2)%ek0
      ceps_2=ceps_2*cim
      cvqu=cf12(i1,i2)%e(0)*p356(0)-cf12(i1,i2)%e(1)*p356(1)-cf1
     & 2(i1,i2)%e(2)*p356(2)-cf12(i1,i2)%e(3)*p356(3)
      cvqd=cf12(i1,i2)%e(0)*p478(0)-cf12(i1,i2)%e(1)*p478(1)-cf1
     & 2(i1,i2)%e(2)*p478(2)-cf12(i1,i2)%e(3)*p478(3)
      cauxa=-cf12(i1,i2)%ek0*quqd+p356k0*cvqd+p478k0*cvqu
      cauxb=-cf12(i1,i2)%ek0*p478(2)+p478k0*cf12(i1,i2)%e(2)
      cauxc=+cf12(i1,i2)%ek0*p356(2)-p356k0*cf12(i1,i2)%e(2)
      u356_12(i1,i2)%a(1,1)=u356_12(i1,i2)%a(1,1)+fcr(id3+1)*(ca
     & uxa+ceps_0)
      u356_12(i1,i2)%a(2,2)=u356_12(i1,i2)%a(2,2)+fcl(id3+1)*(ca
     & uxa-ceps_0)
      u356_12(i1,i2)%b(1,2)=u356_12(i1,i2)%b(1,2)+fcl(id3+1)*(ca
     & uxb-ceps_2)
      u356_12(i1,i2)%b(2,1)=u356_12(i1,i2)%b(2,1)+fcr(id3+1)*(-c
     & auxb-ceps_2)
      u356_12(i1,i2)%c(1,2)=u356_12(i1,i2)%c(1,2)+fcr(id3+1)*(ca
     & uxc+ceps_1)
      u356_12(i1,i2)%c(2,1)=u356_12(i1,i2)%c(2,1)+fcl(id3+1)*(-c
     & auxc+ceps_1)
      u356_12(i1,i2)%d(1,1)=u356_12(i1,i2)%d(1,1)+fcl(id3+1)*cf1
     & 2(i1,i2)%ek0
      u356_12(i1,i2)%d(2,2)=u356_12(i1,i2)%d(2,2)+fcr(id3+1)*cf1
     & 2(i1,i2)%ek0
      end do
      end do
  
* rmh < 0 in our convention means no higgs coupling                     
  
      if (id1.eq.5.and.id3.eq.5.and.rmh.ge.0.d0) then
  
* TH -- qu=p356,qd=p478,a=clineth%a,b=clineth%b,c=clineth%c,coupl=1.d0
      auxa=-p356k0*p478(2)+p478k0*p356(2)
      cauxa=auxa-cim*(p478(3)*p356k0-p356(3)*p478k0)
      clineth%a(1,2)=cauxa
      clineth%a(2,1)=-conjg(cauxa)
      clineth%b(1,1)=p478k0
      clineth%b(2,2)=clineth%b(1,1)
      clineth%c(1,1)=p356k0
      clineth%c(2,2)=clineth%c(1,1)
  
      cfactor= rhbb*rhtt/(s12-cmh2)+rhhbb*rhhtt/(s12-cmhh2)
         do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  u356_12(i1,i2)%a(i3,i4)=ch12(i1,i2)
     &                 *clineth%a(i3,i4)
     &                 *cfactor
                else
                  u356_12(i1,i2)%b(i3,i4)=ch12(i1,i2)
     &                 *clineth%b(i3,i4)
     &                 *cfactor
                  u356_12(i1,i2)%c(i3,i4)=ch12(i1,i2)
     &                 *clineth%c(i3,i4)
     &                 *cfactor
                endif
              enddo
            enddo
          enddo
        enddo
      else
         do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                if (i3.ne.i4) then
                  u356_12(i1,i2)%a(i3,i4)=czero
                else
                  u356_12(i1,i2)%b(i3,i4)=czero
                  u356_12(i1,i2)%c(i3,i4)=czero
                endif
              enddo
            enddo
          enddo
        enddo
      endif
  
  
  
      endif
  
  
* compute all double insertions  of the type li_jklm                    
*     for a Wline                                                       
*                                                                       
*        i __                                                           
*            |___/                                                      
*            |   \                                                      
*            |                                                          
*            |___/                                                      
*            |   \                                                      
*            |                                                          
*                                                                       
  
*     with two z/gamma insertions                                       
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLT0_W -- aa=l5_1234(i1,i2,i3,i4)%a,cc=l5_1234(i1,i2,i3,i4)%c,a1=l5_12
* (i1,i2)%a,c1=l5_12(i1,i2)%c,a2=u512_34(i3,i4)%a,b2=u512_34(i3,i4)%b,c2
* =u512_34(i3,i4)%c,d2=u512_34(i3,i4)%d,prq=p512q,nsum=0
      l5_1234(i1,i2,i3,i4)%c(2)=l5_12(i1,i2)%c(2)*p512q*u512_34(
     & i3,i4)%d(1)+l5_12(i1,i2)%a(2)*u512_34(i3,i4)%c(2)
      l5_1234(i1,i2,i3,i4)%a(2)=l5_12(i1,i2)%c(2)*p512q*u512_34(
     & i3,i4)%b(1)+l5_12(i1,i2)%a(2)*u512_34(i3,i4)%a(2)
      end do
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLT0_W -- aa=l5_1234(i1,i2,i3,i4)%a,cc=l5_1234(i1,i2,i3,i4)%c,a1=l5_34
* (i3,i4)%a,c1=l5_34(i3,i4)%c,a2=u534_12(i1,i2)%a,b2=u534_12(i1,i2)%b,c2
* =u534_12(i1,i2)%c,d2=u534_12(i1,i2)%d,prq=p534q,nsum=1
      l5_1234(i1,i2,i3,i4)%c(2)=l5_1234(i1,i2,i3,i4)%c(2)+l5_34(
     & i3,i4)%c(2)*p534q*u534_12(i1,i2)%d(1)+l5_34(i3,i4)%a(2)*u
     & 534_12(i1,i2)%c(2)
      l5_1234(i1,i2,i3,i4)%a(2)=l5_1234(i1,i2,i3,i4)%a(2)+l5_34(
     & i3,i4)%c(2)*p534q*u534_12(i1,i2)%b(1)+l5_34(i3,i4)%a(2)*u
     & 534_12(i1,i2)%a(2)
      end do
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLT0_W -- aa=l7_1234(i1,i2,i3,i4)%a,cc=l7_1234(i1,i2,i3,i4)%c,a1=l7_12
* (i1,i2)%a,c1=l7_12(i1,i2)%c,a2=u712_34(i3,i4)%a,b2=u712_34(i3,i4)%b,c2
* =u712_34(i3,i4)%c,d2=u712_34(i3,i4)%d,prq=p712q,nsum=0
      l7_1234(i1,i2,i3,i4)%c(2)=l7_12(i1,i2)%c(2)*p712q*u712_34(
     & i3,i4)%d(1)+l7_12(i1,i2)%a(2)*u712_34(i3,i4)%c(2)
      l7_1234(i1,i2,i3,i4)%a(2)=l7_12(i1,i2)%c(2)*p712q*u712_34(
     & i3,i4)%b(1)+l7_12(i1,i2)%a(2)*u712_34(i3,i4)%a(2)
      end do
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLT0_W -- aa=l7_1234(i1,i2,i3,i4)%a,cc=l7_1234(i1,i2,i3,i4)%c,a1=l7_34
* (i3,i4)%a,c1=l7_34(i3,i4)%c,a2=u734_12(i1,i2)%a,b2=u734_12(i1,i2)%b,c2
* =u734_12(i1,i2)%c,d2=u734_12(i1,i2)%d,prq=p734q,nsum=1
      l7_1234(i1,i2,i3,i4)%c(2)=l7_1234(i1,i2,i3,i4)%c(2)+l7_34(
     & i3,i4)%c(2)*p734q*u734_12(i1,i2)%d(1)+l7_34(i3,i4)%a(2)*u
     & 734_12(i1,i2)%c(2)
      l7_1234(i1,i2,i3,i4)%a(2)=l7_1234(i1,i2,i3,i4)%a(2)+l7_34(
     & i3,i4)%c(2)*p734q*u734_12(i1,i2)%b(1)+l7_34(i3,i4)%a(2)*u
     & 734_12(i1,i2)%a(2)
      end do
      end do
      end do
      end do
  
*     with one z/gamma and one w insertion                              
  
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_1278(i1,i2)%a,cc=l5_1278(i1,i2)%c,a1=l5_12(i1,i2)%a,c1
* =l5_12(i1,i2)%c,a2=u512_78%a,b2=u512_78%b,c2=u512_78%c,d2=u512_78%d,pr
* q=p512q,nsum=0
      l5_1278(i1,i2)%c(2)=l5_12(i1,i2)%c(2)*p512q*u512_78%d(1)+l
     & 5_12(i1,i2)%a(2)*u512_78%c(2)
      l5_1278(i1,i2)%a(2)=l5_12(i1,i2)%c(2)*p512q*u512_78%b(1)+l
     & 5_12(i1,i2)%a(2)*u512_78%a(2)
      end do
      end do
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_1278(i1,i2)%a,cc=l5_1278(i1,i2)%c,a1=l5_78%a,c1=l5_78%
* c,a2=u578_12(i1,i2)%a,b2=u578_12(i1,i2)%b,c2=u578_12(i1,i2)%c,d2=u578_
* 12(i1,i2)%d,prq=p578q,nsum=1
      l5_1278(i1,i2)%c(2)=l5_1278(i1,i2)%c(2)+l5_78%c(2)*p578q*u
     & 578_12(i1,i2)%d(1)+l5_78%a(2)*u578_12(i1,i2)%c(2)
      l5_1278(i1,i2)%a(2)=l5_1278(i1,i2)%a(2)+l5_78%c(2)*p578q*u
     & 578_12(i1,i2)%b(1)+l5_78%a(2)*u578_12(i1,i2)%a(2)
      end do
      end do
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_3478(i1,i2)%a,cc=l5_3478(i1,i2)%c,a1=l5_34(i1,i2)%a,c1
* =l5_34(i1,i2)%c,a2=u534_78%a,b2=u534_78%b,c2=u534_78%c,d2=u534_78%d,pr
* q=p534q,nsum=0
      l5_3478(i1,i2)%c(2)=l5_34(i1,i2)%c(2)*p534q*u534_78%d(1)+l
     & 5_34(i1,i2)%a(2)*u534_78%c(2)
      l5_3478(i1,i2)%a(2)=l5_34(i1,i2)%c(2)*p534q*u534_78%b(1)+l
     & 5_34(i1,i2)%a(2)*u534_78%a(2)
      end do
      end do
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l5_3478(i1,i2)%a,cc=l5_3478(i1,i2)%c,a1=l5_78%a,c1=l5_78%
* c,a2=u578_34(i1,i2)%a,b2=u578_34(i1,i2)%b,c2=u578_34(i1,i2)%c,d2=u578_
* 34(i1,i2)%d,prq=p578q,nsum=1
      l5_3478(i1,i2)%c(2)=l5_3478(i1,i2)%c(2)+l5_78%c(2)*p578q*u
     & 578_34(i1,i2)%d(1)+l5_78%a(2)*u578_34(i1,i2)%c(2)
      l5_3478(i1,i2)%a(2)=l5_3478(i1,i2)%a(2)+l5_78%c(2)*p578q*u
     & 578_34(i1,i2)%b(1)+l5_78%a(2)*u578_34(i1,i2)%a(2)
      end do
      end do
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_1256(i1,i2)%a,cc=l7_1256(i1,i2)%c,a1=l7_12(i1,i2)%a,c1
* =l7_12(i1,i2)%c,a2=u712_56%a,b2=u712_56%b,c2=u712_56%c,d2=u712_56%d,pr
* q=p712q,nsum=0
      l7_1256(i1,i2)%c(2)=l7_12(i1,i2)%c(2)*p712q*u712_56%d(1)+l
     & 7_12(i1,i2)%a(2)*u712_56%c(2)
      l7_1256(i1,i2)%a(2)=l7_12(i1,i2)%c(2)*p712q*u712_56%b(1)+l
     & 7_12(i1,i2)%a(2)*u712_56%a(2)
      end do
      end do
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_1256(i1,i2)%a,cc=l7_1256(i1,i2)%c,a1=l7_56%a,c1=l7_56%
* c,a2=u756_12(i1,i2)%a,b2=u756_12(i1,i2)%b,c2=u756_12(i1,i2)%c,d2=u756_
* 12(i1,i2)%d,prq=p756q,nsum=1
      l7_1256(i1,i2)%c(2)=l7_1256(i1,i2)%c(2)+l7_56%c(2)*p756q*u
     & 756_12(i1,i2)%d(1)+l7_56%a(2)*u756_12(i1,i2)%c(2)
      l7_1256(i1,i2)%a(2)=l7_1256(i1,i2)%a(2)+l7_56%c(2)*p756q*u
     & 756_12(i1,i2)%b(1)+l7_56%a(2)*u756_12(i1,i2)%a(2)
      end do
      end do
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_3456(i1,i2)%a,cc=l7_3456(i1,i2)%c,a1=l7_34(i1,i2)%a,c1
* =l7_34(i1,i2)%c,a2=u734_56%a,b2=u734_56%b,c2=u734_56%c,d2=u734_56%d,pr
* q=p734q,nsum=0
      l7_3456(i1,i2)%c(2)=l7_34(i1,i2)%c(2)*p734q*u734_56%d(1)+l
     & 7_34(i1,i2)%a(2)*u734_56%c(2)
      l7_3456(i1,i2)%a(2)=l7_34(i1,i2)%c(2)*p734q*u734_56%b(1)+l
     & 7_34(i1,i2)%a(2)*u734_56%a(2)
      end do
      end do
      do i1=1,2
      do i2=1,2
* TLT0_W -- aa=l7_3456(i1,i2)%a,cc=l7_3456(i1,i2)%c,a1=l7_56%a,c1=l7_56%
* c,a2=u756_34(i1,i2)%a,b2=u756_34(i1,i2)%b,c2=u756_34(i1,i2)%c,d2=u756_
* 34(i1,i2)%d,prq=p756q,nsum=1
      l7_3456(i1,i2)%c(2)=l7_3456(i1,i2)%c(2)+l7_56%c(2)*p756q*u
     & 756_34(i1,i2)%d(1)+l7_56%a(2)*u756_34(i1,i2)%c(2)
      l7_3456(i1,i2)%a(2)=l7_3456(i1,i2)%a(2)+l7_56%c(2)*p756q*u
     & 756_34(i1,i2)%b(1)+l7_56%a(2)*u756_34(i1,i2)%a(2)
      end do
      end do
  
* compute all double insertions  of the type li_jklm                    
*     for a Zline                                                       
*                                                                       
*        i __                                                           
*            |___/                                                      
*            |   \                                                      
*            |                                                          
*            |___/                                                      
*            |   \                                                      
*            |                                                          
*                                                                       
  
*     with two w insertions                                             
  
      if(iup(id1).eq.1) then
* TWTW -- aa=l1_5678%a,bb=l1_5678%b,cc=l1_5678%c,dd=l1_5678%d,a1=l1_78%a
* ,b1=l1_78%b,c1=l1_78%c,d1=l1_78%d,a2=u178_56%a,b2=u178_56%b,c2=u178_56
* %c,d2=u178_56%d,prq=p178q,nsum=0
      l1_5678%d(1,1)=l1_78%d(1,1)*p178q*u178_56%d(1,1)+l1_78%b(1
     & ,2)*u178_56%c(2,1)
      l1_5678%b(1,2)=l1_78%d(1,1)*p178q*u178_56%b(1,2)+l1_78%b(1
     & ,2)*u178_56%a(2,2)
      l1_5678%c(2,1)=l1_78%c(2,1)*p178q*u178_56%d(1,1)+l1_78%a(2
     & ,2)*u178_56%c(2,1)
      l1_5678%a(2,2)=l1_78%c(2,1)*p178q*u178_56%b(1,2)+l1_78%a(2
     & ,2)*u178_56%a(2,2)
      else
* TWTW -- aa=l1_5678%a,bb=l1_5678%b,cc=l1_5678%c,dd=l1_5678%d,a1=l1_56%a
* ,b1=l1_56%b,c1=l1_56%c,d1=l1_56%d,a2=u156_78%a,b2=u156_78%b,c2=u156_78
* %c,d2=u156_78%d,prq=p156q,nsum=0
      l1_5678%d(1,1)=l1_56%d(1,1)*p156q*u156_78%d(1,1)+l1_56%b(1
     & ,2)*u156_78%c(2,1)
      l1_5678%b(1,2)=l1_56%d(1,1)*p156q*u156_78%b(1,2)+l1_56%b(1
     & ,2)*u156_78%a(2,2)
      l1_5678%c(2,1)=l1_56%c(2,1)*p156q*u156_78%d(1,1)+l1_56%a(2
     & ,2)*u156_78%c(2,1)
      l1_5678%a(2,2)=l1_56%c(2,1)*p156q*u156_78%b(1,2)+l1_56%a(2
     & ,2)*u156_78%a(2,2)
      endif
      if(iup(id3).eq.1) then
* TWTW -- aa=l3_5678%a,bb=l3_5678%b,cc=l3_5678%c,dd=l3_5678%d,a1=l3_78%a
* ,b1=l3_78%b,c1=l3_78%c,d1=l3_78%d,a2=u378_56%a,b2=u378_56%b,c2=u378_56
* %c,d2=u378_56%d,prq=p378q,nsum=0
      l3_5678%d(1,1)=l3_78%d(1,1)*p378q*u378_56%d(1,1)+l3_78%b(1
     & ,2)*u378_56%c(2,1)
      l3_5678%b(1,2)=l3_78%d(1,1)*p378q*u378_56%b(1,2)+l3_78%b(1
     & ,2)*u378_56%a(2,2)
      l3_5678%c(2,1)=l3_78%c(2,1)*p378q*u378_56%d(1,1)+l3_78%a(2
     & ,2)*u378_56%c(2,1)
      l3_5678%a(2,2)=l3_78%c(2,1)*p378q*u378_56%b(1,2)+l3_78%a(2
     & ,2)*u378_56%a(2,2)
      else
* TWTW -- aa=l3_5678%a,bb=l3_5678%b,cc=l3_5678%c,dd=l3_5678%d,a1=l3_56%a
* ,b1=l3_56%b,c1=l3_56%c,d1=l3_56%d,a2=u356_78%a,b2=u356_78%b,c2=u356_78
* %c,d2=u356_78%d,prq=p356q,nsum=0
      l3_5678%d(1,1)=l3_56%d(1,1)*p356q*u356_78%d(1,1)+l3_56%b(1
     & ,2)*u356_78%c(2,1)
      l3_5678%b(1,2)=l3_56%d(1,1)*p356q*u356_78%b(1,2)+l3_56%b(1
     & ,2)*u356_78%a(2,2)
      l3_5678%c(2,1)=l3_56%c(2,1)*p356q*u356_78%d(1,1)+l3_56%a(2
     & ,2)*u356_78%c(2,1)
      l3_5678%a(2,2)=l3_56%c(2,1)*p356q*u356_78%b(1,2)+l3_56%a(2
     & ,2)*u356_78%a(2,2)
      endif
  
*     with one w and one z,gamma,h insertion                            
  
  
      if(iup(id1).eq.1) then
        rmassexc=rmass(id1-1)
      else
        rmassexc=rmass(id1+1)
      endif
  
      if (id1.eq.5.and.id3.eq.5.and.rmh.ge.0.d0) then
* rmh < 0 in our convention means no higgs coupling                     
* use twts because t is the sum of gamma zeta and higgs                 
      do i3=1,2
      do i4=1,2
* TWTs -- aa=l1_3456(i3,i4)%a,bb=l1_3456(i3,i4)%b,cc=l1_3456(i3,i4)%c,dd
* =l1_3456(i3,i4)%d,a1=l1_56%a,b1=l1_56%b,c1=l1_56%c,d1=l1_56%d,a2=u156_
* 34(i3,i4)%a,b2=u156_34(i3,i4)%b,c2=u156_34(i3,i4)%c,d2=u156_34(i3,i4)%
* d,prq=p156q,m=rmassexc,nsum=0
      do jut=1,2
      cx1=u156_34(i3,i4)%a(1,jut)+rmassexc*u156_34(i3,i4)%b(1,ju
     & t)
      cx2=u156_34(i3,i4)%a(2,jut)+rmassexc*u156_34(i3,i4)%b(2,ju
     & t)
      cy1=p156q*u156_34(i3,i4)%b(1,jut)+rmassexc*u156_34(i3,i4)%
     & a(1,jut)
      cy2=p156q*u156_34(i3,i4)%b(2,jut)+rmassexc*u156_34(i3,i4)%
     & a(2,jut)
      cw1=u156_34(i3,i4)%c(1,jut)+rmassexc*u156_34(i3,i4)%d(1,ju
     & t)
      cw2=u156_34(i3,i4)%c(2,jut)+rmassexc*u156_34(i3,i4)%d(2,ju
     & t)
      cz1=p156q*u156_34(i3,i4)%d(1,jut)+rmassexc*u156_34(i3,i4)%
     & c(1,jut)
      cz2=p156q*u156_34(i3,i4)%d(2,jut)+rmassexc*u156_34(i3,i4)%
     & c(2,jut)
      l1_3456(i3,i4)%b(1,jut)=l1_56%d(1,1)*cy1+l1_56%b(1,2)*cx2
      l1_3456(i3,i4)%d(1,jut)=l1_56%d(1,1)*cz1+l1_56%b(1,2)*cw2
      l1_3456(i3,i4)%a(2,jut)=l1_56%c(2,1)*cy1+l1_56%a(2,2)*cx2
      l1_3456(i3,i4)%c(2,jut)=l1_56%c(2,1)*cz1+l1_56%a(2,2)*cw2
      end do
      end do
      end do
  
* set = 0 components of  l1_3456(i3,i4) which have not been             
*  just computed                                                        
  
      do i3=1,2
        do i4=1,2
          do iut=1,2
            l1_3456(i3,i4)%a(1,iut)=czero
            l1_3456(i3,i4)%b(2,iut)=czero
            l1_3456(i3,i4)%c(1,iut)=czero
            l1_3456(i3,i4)%d(2,iut)=czero
          enddo
        enddo
      enddo
  
      do i3=1,2
      do i4=1,2
* TsTW -- aa=l1_3456(i3,i4)%a,bb=l1_3456(i3,i4)%b,cc=l1_3456(i3,i4)%c,dd
* =l1_3456(i3,i4)%d,a1=l1_34(i3,i4)%a,b1=l1_34(i3,i4)%b,c1=l1_34(i3,i4)%
* c,d1=l1_34(i3,i4)%d,a2=u134_56%a,b2=u134_56%b,c2=u134_56%c,d2=u134_56%
* d,prq=p134q,m=rmass(id1),nsum=1
      do iut=1,2
      cw1=l1_34(i3,i4)%c(iut,1)*p134q+l1_34(i3,i4)%a(iut,1)*rmas
     & s(id1)
      cw2=l1_34(i3,i4)%a(iut,2)+l1_34(i3,i4)%c(iut,2)*rmass(id1)
      cz1=l1_34(i3,i4)%d(iut,1)*p134q+l1_34(i3,i4)%b(iut,1)*rmas
     & s(id1)
      cz2=l1_34(i3,i4)%b(iut,2)+l1_34(i3,i4)%d(iut,2)*rmass(id1)
      l1_3456(i3,i4)%c(iut,1)=l1_3456(i3,i4)%c(iut,1)+cw1*u134_5
     & 6%d(1,1)+cw2*u134_56%c(2,1)
      l1_3456(i3,i4)%d(iut,1)=l1_3456(i3,i4)%d(iut,1)+cz1*u134_5
     & 6%d(1,1)+cz2*u134_56%c(2,1)
      l1_3456(i3,i4)%a(iut,2)=l1_3456(i3,i4)%a(iut,2)+cw1*u134_5
     & 6%b(1,2)+cw2*u134_56%a(2,2)
      l1_3456(i3,i4)%b(iut,2)=l1_3456(i3,i4)%b(iut,2)+cz1*u134_5
     & 6%b(1,2)+cz2*u134_56%a(2,2)
      end do
      end do
      end do
  
      else
  
        if(iup(id1).eq.1) then
  
      do i3=1,2
      do i4=1,2
* TWT -- aa=l1_3478(i3,i4)%a,bb=l1_3478(i3,i4)%b,cc=l1_3478(i3,i4)%c,dd=
* l1_3478(i3,i4)%d,a1=l1_78%a,b1=l1_78%b,c1=l1_78%c,d1=l1_78%d,a2=u178_3
* 4(i3,i4)%a,b2=u178_34(i3,i4)%b,c2=u178_34(i3,i4)%c,d2=u178_34(i3,i4)%d
* ,prq=p178q,m=rmassexc,nsum=0
      l1_3478(i3,i4)%b(1,1)=rmassexc*(l1_78%d(1,1)*u178_34(i3,i4
     & )%a(1,1)+l1_78%b(1,2)*u178_34(i3,i4)%b(2,1))
      l1_3478(i3,i4)%d(1,1)=l1_78%d(1,1)*p178q*u178_34(i3,i4)%d(
     & 1,1)+l1_78%b(1,2)*u178_34(i3,i4)%c(2,1)
      l1_3478(i3,i4)%b(1,2)=l1_78%d(1,1)*p178q*u178_34(i3,i4)%b(
     & 1,2)+l1_78%b(1,2)*u178_34(i3,i4)%a(2,2)
      l1_3478(i3,i4)%d(1,2)=rmassexc*(l1_78%d(1,1)*u178_34(i3,i4
     & )%c(1,2)+l1_78%b(1,2)*u178_34(i3,i4)%d(2,2))
      l1_3478(i3,i4)%a(2,1)=rmassexc*(l1_78%c(2,1)*u178_34(i3,i4
     & )%a(1,1)+l1_78%a(2,2)*u178_34(i3,i4)%b(2,1))
      l1_3478(i3,i4)%c(2,1)=l1_78%c(2,1)*p178q*u178_34(i3,i4)%d(
     & 1,1)+l1_78%a(2,2)*u178_34(i3,i4)%c(2,1)
      l1_3478(i3,i4)%a(2,2)=l1_78%c(2,1)*p178q*u178_34(i3,i4)%b(
     & 1,2)+l1_78%a(2,2)*u178_34(i3,i4)%a(2,2)
      l1_3478(i3,i4)%c(2,2)=rmassexc*(l1_78%c(2,1)*u178_34(i3,i4
     & )%c(1,2)+l1_78%a(2,2)*u178_34(i3,i4)%d(2,2))
      end do
      end do
  
* set = 0 components of  l1_3478(i3,i4) which have not been             
*  just computed                                                        
  
      do i3=1,2
        do i4=1,2
          do iut=1,2
            l1_3478(i3,i4)%a(1,iut)=czero
            l1_3478(i3,i4)%b(2,iut)=czero
            l1_3478(i3,i4)%c(1,iut)=czero
            l1_3478(i3,i4)%d(2,iut)=czero
          enddo
        enddo
      enddo
  
  
      do i3=1,2
      do i4=1,2
* TTW -- aa=l1_3478(i3,i4)%a,bb=l1_3478(i3,i4)%b,cc=l1_3478(i3,i4)%c,dd=
* l1_3478(i3,i4)%d,a1=l1_34(i3,i4)%a,b1=l1_34(i3,i4)%b,c1=l1_34(i3,i4)%c
* ,d1=l1_34(i3,i4)%d,a2=u134_78%a,b2=u134_78%b,c2=u134_78%c,d2=u134_78%d
* ,prq=p134q,m=rmass(id1),nsum=1
      l1_3478(i3,i4)%c(1,1)=l1_3478(i3,i4)%c(1,1)+rmass(id1)*(l1
     & _34(i3,i4)%a(1,1)*u134_78%d(1,1)+l1_34(i3,i4)%c(1,2)*u134
     & _78%c(2,1))
      l1_3478(i3,i4)%d(1,1)=l1_3478(i3,i4)%d(1,1)+l1_34(i3,i4)%d
     & (1,1)*p134q*u134_78%d(1,1)+l1_34(i3,i4)%b(1,2)*u134_78%c(
     & 2,1)
      l1_3478(i3,i4)%a(1,2)=l1_3478(i3,i4)%a(1,2)+rmass(id1)*(l1
     & _34(i3,i4)%a(1,1)*u134_78%b(1,2)+l1_34(i3,i4)%c(1,2)*u134
     & _78%a(2,2))
      l1_3478(i3,i4)%b(1,2)=l1_3478(i3,i4)%b(1,2)+l1_34(i3,i4)%d
     & (1,1)*p134q*u134_78%b(1,2)+l1_34(i3,i4)%b(1,2)*u134_78%a(
     & 2,2)
      l1_3478(i3,i4)%c(2,1)=l1_3478(i3,i4)%c(2,1)+l1_34(i3,i4)%c
     & (2,1)*p134q*u134_78%d(1,1)+l1_34(i3,i4)%a(2,2)*u134_78%c(
     & 2,1)
      l1_3478(i3,i4)%d(2,1)=l1_3478(i3,i4)%d(2,1)+rmass(id1)*(l1
     & _34(i3,i4)%b(2,1)*u134_78%d(1,1)+l1_34(i3,i4)%d(2,2)*u134
     & _78%c(2,1))
      l1_3478(i3,i4)%a(2,2)=l1_3478(i3,i4)%a(2,2)+l1_34(i3,i4)%c
     & (2,1)*p134q*u134_78%b(1,2)+l1_34(i3,i4)%a(2,2)*u134_78%a(
     & 2,2)
      l1_3478(i3,i4)%b(2,2)=l1_3478(i3,i4)%b(2,2)+rmass(id1)*(l1
     & _34(i3,i4)%b(2,1)*u134_78%b(1,2)+l1_34(i3,i4)%d(2,2)*u134
     & _78%a(2,2))
      end do
      end do
  
        else
  
      do i3=1,2
      do i4=1,2
* TWT -- aa=l1_3456(i3,i4)%a,bb=l1_3456(i3,i4)%b,cc=l1_3456(i3,i4)%c,dd=
* l1_3456(i3,i4)%d,a1=l1_56%a,b1=l1_56%b,c1=l1_56%c,d1=l1_56%d,a2=u156_3
* 4(i3,i4)%a,b2=u156_34(i3,i4)%b,c2=u156_34(i3,i4)%c,d2=u156_34(i3,i4)%d
* ,prq=p156q,m=rmassexc,nsum=0
      l1_3456(i3,i4)%b(1,1)=rmassexc*(l1_56%d(1,1)*u156_34(i3,i4
     & )%a(1,1)+l1_56%b(1,2)*u156_34(i3,i4)%b(2,1))
      l1_3456(i3,i4)%d(1,1)=l1_56%d(1,1)*p156q*u156_34(i3,i4)%d(
     & 1,1)+l1_56%b(1,2)*u156_34(i3,i4)%c(2,1)
      l1_3456(i3,i4)%b(1,2)=l1_56%d(1,1)*p156q*u156_34(i3,i4)%b(
     & 1,2)+l1_56%b(1,2)*u156_34(i3,i4)%a(2,2)
      l1_3456(i3,i4)%d(1,2)=rmassexc*(l1_56%d(1,1)*u156_34(i3,i4
     & )%c(1,2)+l1_56%b(1,2)*u156_34(i3,i4)%d(2,2))
      l1_3456(i3,i4)%a(2,1)=rmassexc*(l1_56%c(2,1)*u156_34(i3,i4
     & )%a(1,1)+l1_56%a(2,2)*u156_34(i3,i4)%b(2,1))
      l1_3456(i3,i4)%c(2,1)=l1_56%c(2,1)*p156q*u156_34(i3,i4)%d(
     & 1,1)+l1_56%a(2,2)*u156_34(i3,i4)%c(2,1)
      l1_3456(i3,i4)%a(2,2)=l1_56%c(2,1)*p156q*u156_34(i3,i4)%b(
     & 1,2)+l1_56%a(2,2)*u156_34(i3,i4)%a(2,2)
      l1_3456(i3,i4)%c(2,2)=rmassexc*(l1_56%c(2,1)*u156_34(i3,i4
     & )%c(1,2)+l1_56%a(2,2)*u156_34(i3,i4)%d(2,2))
      end do
      end do
  
* set = 0 components of  l1_3456(i3,i4) which have not been             
*  just computed                                                        
  
      do i3=1,2
        do i4=1,2
          do iut=1,2
            l1_3456(i3,i4)%a(1,iut)=czero
            l1_3456(i3,i4)%b(2,iut)=czero
            l1_3456(i3,i4)%c(1,iut)=czero
            l1_3456(i3,i4)%d(2,iut)=czero
          enddo
        enddo
      enddo
  
  
      do i3=1,2
      do i4=1,2
* TTW -- aa=l1_3456(i3,i4)%a,bb=l1_3456(i3,i4)%b,cc=l1_3456(i3,i4)%c,dd=
* l1_3456(i3,i4)%d,a1=l1_34(i3,i4)%a,b1=l1_34(i3,i4)%b,c1=l1_34(i3,i4)%c
* ,d1=l1_34(i3,i4)%d,a2=u134_56%a,b2=u134_56%b,c2=u134_56%c,d2=u134_56%d
* ,prq=p134q,m=rmass(id1),nsum=1
      l1_3456(i3,i4)%c(1,1)=l1_3456(i3,i4)%c(1,1)+rmass(id1)*(l1
     & _34(i3,i4)%a(1,1)*u134_56%d(1,1)+l1_34(i3,i4)%c(1,2)*u134
     & _56%c(2,1))
      l1_3456(i3,i4)%d(1,1)=l1_3456(i3,i4)%d(1,1)+l1_34(i3,i4)%d
     & (1,1)*p134q*u134_56%d(1,1)+l1_34(i3,i4)%b(1,2)*u134_56%c(
     & 2,1)
      l1_3456(i3,i4)%a(1,2)=l1_3456(i3,i4)%a(1,2)+rmass(id1)*(l1
     & _34(i3,i4)%a(1,1)*u134_56%b(1,2)+l1_34(i3,i4)%c(1,2)*u134
     & _56%a(2,2))
      l1_3456(i3,i4)%b(1,2)=l1_3456(i3,i4)%b(1,2)+l1_34(i3,i4)%d
     & (1,1)*p134q*u134_56%b(1,2)+l1_34(i3,i4)%b(1,2)*u134_56%a(
     & 2,2)
      l1_3456(i3,i4)%c(2,1)=l1_3456(i3,i4)%c(2,1)+l1_34(i3,i4)%c
     & (2,1)*p134q*u134_56%d(1,1)+l1_34(i3,i4)%a(2,2)*u134_56%c(
     & 2,1)
      l1_3456(i3,i4)%d(2,1)=l1_3456(i3,i4)%d(2,1)+rmass(id1)*(l1
     & _34(i3,i4)%b(2,1)*u134_56%d(1,1)+l1_34(i3,i4)%d(2,2)*u134
     & _56%c(2,1))
      l1_3456(i3,i4)%a(2,2)=l1_3456(i3,i4)%a(2,2)+l1_34(i3,i4)%c
     & (2,1)*p134q*u134_56%b(1,2)+l1_34(i3,i4)%a(2,2)*u134_56%a(
     & 2,2)
      l1_3456(i3,i4)%b(2,2)=l1_3456(i3,i4)%b(2,2)+rmass(id1)*(l1
     & _34(i3,i4)%b(2,1)*u134_56%b(1,2)+l1_34(i3,i4)%d(2,2)*u134
     & _56%a(2,2))
      end do
      end do
  
        endif
      endif
  
  
      if(iup(id3).eq.1) then
        rmassexc=rmass(id3-1)
      else
        rmassexc=rmass(id3+1)
      endif
  
      if (id1.eq.5.and.id3.eq.5.and.rmh.ge.0.d0) then
* rmh < 0 in our convention means no higgs coupling                     
* use twts because t is the sum of gamma zeta and higgs                 
      do i1=1,2
      do i2=1,2
* TWTs -- aa=l3_1256(i1,i2)%a,bb=l3_1256(i1,i2)%b,cc=l3_1256(i1,i2)%c,dd
* =l3_1256(i1,i2)%d,a1=l3_56%a,b1=l3_56%b,c1=l3_56%c,d1=l3_56%d,a2=u356_
* 12(i1,i2)%a,b2=u356_12(i1,i2)%b,c2=u356_12(i1,i2)%c,d2=u356_12(i1,i2)%
* d,prq=p356q,m=rmassexc,nsum=0
      do jut=1,2
      cx1=u356_12(i1,i2)%a(1,jut)+rmassexc*u356_12(i1,i2)%b(1,ju
     & t)
      cx2=u356_12(i1,i2)%a(2,jut)+rmassexc*u356_12(i1,i2)%b(2,ju
     & t)
      cy1=p356q*u356_12(i1,i2)%b(1,jut)+rmassexc*u356_12(i1,i2)%
     & a(1,jut)
      cy2=p356q*u356_12(i1,i2)%b(2,jut)+rmassexc*u356_12(i1,i2)%
     & a(2,jut)
      cw1=u356_12(i1,i2)%c(1,jut)+rmassexc*u356_12(i1,i2)%d(1,ju
     & t)
      cw2=u356_12(i1,i2)%c(2,jut)+rmassexc*u356_12(i1,i2)%d(2,ju
     & t)
      cz1=p356q*u356_12(i1,i2)%d(1,jut)+rmassexc*u356_12(i1,i2)%
     & c(1,jut)
      cz2=p356q*u356_12(i1,i2)%d(2,jut)+rmassexc*u356_12(i1,i2)%
     & c(2,jut)
      l3_1256(i1,i2)%b(1,jut)=l3_56%d(1,1)*cy1+l3_56%b(1,2)*cx2
      l3_1256(i1,i2)%d(1,jut)=l3_56%d(1,1)*cz1+l3_56%b(1,2)*cw2
      l3_1256(i1,i2)%a(2,jut)=l3_56%c(2,1)*cy1+l3_56%a(2,2)*cx2
      l3_1256(i1,i2)%c(2,jut)=l3_56%c(2,1)*cz1+l3_56%a(2,2)*cw2
      end do
      end do
      end do
  
* set = 0 components of  l3_1256(i1,i2) which have not been             
*  just computed                                                        
  
      do i1=1,2
        do i2=1,2
          do iut=1,2
            l3_1256(i1,i2)%a(1,iut)=czero
            l3_1256(i1,i2)%b(2,iut)=czero
            l3_1256(i1,i2)%c(1,iut)=czero
            l3_1256(i1,i2)%d(2,iut)=czero
          enddo
        enddo
      enddo
  
      do i1=1,2
      do i2=1,2
* TsTW -- aa=l3_1256(i1,i2)%a,bb=l3_1256(i1,i2)%b,cc=l3_1256(i1,i2)%c,dd
* =l3_1256(i1,i2)%d,a1=l3_12(i1,i2)%a,b1=l3_12(i1,i2)%b,c1=l3_12(i1,i2)%
* c,d1=l3_12(i1,i2)%d,a2=u312_56%a,b2=u312_56%b,c2=u312_56%c,d2=u312_56%
* d,prq=p312q,m=rmass(id3),nsum=1
      do iut=1,2
      cw1=l3_12(i1,i2)%c(iut,1)*p312q+l3_12(i1,i2)%a(iut,1)*rmas
     & s(id3)
      cw2=l3_12(i1,i2)%a(iut,2)+l3_12(i1,i2)%c(iut,2)*rmass(id3)
      cz1=l3_12(i1,i2)%d(iut,1)*p312q+l3_12(i1,i2)%b(iut,1)*rmas
     & s(id3)
      cz2=l3_12(i1,i2)%b(iut,2)+l3_12(i1,i2)%d(iut,2)*rmass(id3)
      l3_1256(i1,i2)%c(iut,1)=l3_1256(i1,i2)%c(iut,1)+cw1*u312_5
     & 6%d(1,1)+cw2*u312_56%c(2,1)
      l3_1256(i1,i2)%d(iut,1)=l3_1256(i1,i2)%d(iut,1)+cz1*u312_5
     & 6%d(1,1)+cz2*u312_56%c(2,1)
      l3_1256(i1,i2)%a(iut,2)=l3_1256(i1,i2)%a(iut,2)+cw1*u312_5
     & 6%b(1,2)+cw2*u312_56%a(2,2)
      l3_1256(i1,i2)%b(iut,2)=l3_1256(i1,i2)%b(iut,2)+cz1*u312_5
     & 6%b(1,2)+cz2*u312_56%a(2,2)
      end do
      end do
      end do
  
      else
  
        if(iup(id3).eq.1) then
  
      do i1=1,2
      do i2=1,2
* TWT -- aa=l3_1278(i1,i2)%a,bb=l3_1278(i1,i2)%b,cc=l3_1278(i1,i2)%c,dd=
* l3_1278(i1,i2)%d,a1=l3_78%a,b1=l3_78%b,c1=l3_78%c,d1=l3_78%d,a2=u378_1
* 2(i1,i2)%a,b2=u378_12(i1,i2)%b,c2=u378_12(i1,i2)%c,d2=u378_12(i1,i2)%d
* ,prq=p378q,m=rmassexc,nsum=0
      l3_1278(i1,i2)%b(1,1)=rmassexc*(l3_78%d(1,1)*u378_12(i1,i2
     & )%a(1,1)+l3_78%b(1,2)*u378_12(i1,i2)%b(2,1))
      l3_1278(i1,i2)%d(1,1)=l3_78%d(1,1)*p378q*u378_12(i1,i2)%d(
     & 1,1)+l3_78%b(1,2)*u378_12(i1,i2)%c(2,1)
      l3_1278(i1,i2)%b(1,2)=l3_78%d(1,1)*p378q*u378_12(i1,i2)%b(
     & 1,2)+l3_78%b(1,2)*u378_12(i1,i2)%a(2,2)
      l3_1278(i1,i2)%d(1,2)=rmassexc*(l3_78%d(1,1)*u378_12(i1,i2
     & )%c(1,2)+l3_78%b(1,2)*u378_12(i1,i2)%d(2,2))
      l3_1278(i1,i2)%a(2,1)=rmassexc*(l3_78%c(2,1)*u378_12(i1,i2
     & )%a(1,1)+l3_78%a(2,2)*u378_12(i1,i2)%b(2,1))
      l3_1278(i1,i2)%c(2,1)=l3_78%c(2,1)*p378q*u378_12(i1,i2)%d(
     & 1,1)+l3_78%a(2,2)*u378_12(i1,i2)%c(2,1)
      l3_1278(i1,i2)%a(2,2)=l3_78%c(2,1)*p378q*u378_12(i1,i2)%b(
     & 1,2)+l3_78%a(2,2)*u378_12(i1,i2)%a(2,2)
      l3_1278(i1,i2)%c(2,2)=rmassexc*(l3_78%c(2,1)*u378_12(i1,i2
     & )%c(1,2)+l3_78%a(2,2)*u378_12(i1,i2)%d(2,2))
      end do
      end do
  
* set = 0 components of  l3_1278(i1,i2) which have not been             
*  just computed                                                        
  
      do i1=1,2
        do i2=1,2
          do iut=1,2
            l3_1278(i1,i2)%a(1,iut)=czero
            l3_1278(i1,i2)%b(2,iut)=czero
            l3_1278(i1,i2)%c(1,iut)=czero
            l3_1278(i1,i2)%d(2,iut)=czero
          enddo
        enddo
      enddo
  
  
      do i1=1,2
      do i2=1,2
* TTW -- aa=l3_1278(i1,i2)%a,bb=l3_1278(i1,i2)%b,cc=l3_1278(i1,i2)%c,dd=
* l3_1278(i1,i2)%d,a1=l3_12(i1,i2)%a,b1=l3_12(i1,i2)%b,c1=l3_12(i1,i2)%c
* ,d1=l3_12(i1,i2)%d,a2=u312_78%a,b2=u312_78%b,c2=u312_78%c,d2=u312_78%d
* ,prq=p312q,m=rmass(id3),nsum=1
      l3_1278(i1,i2)%c(1,1)=l3_1278(i1,i2)%c(1,1)+rmass(id3)*(l3
     & _12(i1,i2)%a(1,1)*u312_78%d(1,1)+l3_12(i1,i2)%c(1,2)*u312
     & _78%c(2,1))
      l3_1278(i1,i2)%d(1,1)=l3_1278(i1,i2)%d(1,1)+l3_12(i1,i2)%d
     & (1,1)*p312q*u312_78%d(1,1)+l3_12(i1,i2)%b(1,2)*u312_78%c(
     & 2,1)
      l3_1278(i1,i2)%a(1,2)=l3_1278(i1,i2)%a(1,2)+rmass(id3)*(l3
     & _12(i1,i2)%a(1,1)*u312_78%b(1,2)+l3_12(i1,i2)%c(1,2)*u312
     & _78%a(2,2))
      l3_1278(i1,i2)%b(1,2)=l3_1278(i1,i2)%b(1,2)+l3_12(i1,i2)%d
     & (1,1)*p312q*u312_78%b(1,2)+l3_12(i1,i2)%b(1,2)*u312_78%a(
     & 2,2)
      l3_1278(i1,i2)%c(2,1)=l3_1278(i1,i2)%c(2,1)+l3_12(i1,i2)%c
     & (2,1)*p312q*u312_78%d(1,1)+l3_12(i1,i2)%a(2,2)*u312_78%c(
     & 2,1)
      l3_1278(i1,i2)%d(2,1)=l3_1278(i1,i2)%d(2,1)+rmass(id3)*(l3
     & _12(i1,i2)%b(2,1)*u312_78%d(1,1)+l3_12(i1,i2)%d(2,2)*u312
     & _78%c(2,1))
      l3_1278(i1,i2)%a(2,2)=l3_1278(i1,i2)%a(2,2)+l3_12(i1,i2)%c
     & (2,1)*p312q*u312_78%b(1,2)+l3_12(i1,i2)%a(2,2)*u312_78%a(
     & 2,2)
      l3_1278(i1,i2)%b(2,2)=l3_1278(i1,i2)%b(2,2)+rmass(id3)*(l3
     & _12(i1,i2)%b(2,1)*u312_78%b(1,2)+l3_12(i1,i2)%d(2,2)*u312
     & _78%a(2,2))
      end do
      end do
  
        else
  
      do i1=1,2
      do i2=1,2
* TWT -- aa=l3_1256(i1,i2)%a,bb=l3_1256(i1,i2)%b,cc=l3_1256(i1,i2)%c,dd=
* l3_1256(i1,i2)%d,a1=l3_56%a,b1=l3_56%b,c1=l3_56%c,d1=l3_56%d,a2=u356_1
* 2(i1,i2)%a,b2=u356_12(i1,i2)%b,c2=u356_12(i1,i2)%c,d2=u356_12(i1,i2)%d
* ,prq=p356q,m=rmassexc,nsum=0
      l3_1256(i1,i2)%b(1,1)=rmassexc*(l3_56%d(1,1)*u356_12(i1,i2
     & )%a(1,1)+l3_56%b(1,2)*u356_12(i1,i2)%b(2,1))
      l3_1256(i1,i2)%d(1,1)=l3_56%d(1,1)*p356q*u356_12(i1,i2)%d(
     & 1,1)+l3_56%b(1,2)*u356_12(i1,i2)%c(2,1)
      l3_1256(i1,i2)%b(1,2)=l3_56%d(1,1)*p356q*u356_12(i1,i2)%b(
     & 1,2)+l3_56%b(1,2)*u356_12(i1,i2)%a(2,2)
      l3_1256(i1,i2)%d(1,2)=rmassexc*(l3_56%d(1,1)*u356_12(i1,i2
     & )%c(1,2)+l3_56%b(1,2)*u356_12(i1,i2)%d(2,2))
      l3_1256(i1,i2)%a(2,1)=rmassexc*(l3_56%c(2,1)*u356_12(i1,i2
     & )%a(1,1)+l3_56%a(2,2)*u356_12(i1,i2)%b(2,1))
      l3_1256(i1,i2)%c(2,1)=l3_56%c(2,1)*p356q*u356_12(i1,i2)%d(
     & 1,1)+l3_56%a(2,2)*u356_12(i1,i2)%c(2,1)
      l3_1256(i1,i2)%a(2,2)=l3_56%c(2,1)*p356q*u356_12(i1,i2)%b(
     & 1,2)+l3_56%a(2,2)*u356_12(i1,i2)%a(2,2)
      l3_1256(i1,i2)%c(2,2)=rmassexc*(l3_56%c(2,1)*u356_12(i1,i2
     & )%c(1,2)+l3_56%a(2,2)*u356_12(i1,i2)%d(2,2))
      end do
      end do
  
* set = 0 components of  l3_1256(i1,i2) which have not been             
*  just computed                                                        
  
      do i1=1,2
        do i2=1,2
          do iut=1,2
            l3_1256(i1,i2)%a(1,iut)=czero
            l3_1256(i1,i2)%b(2,iut)=czero
            l3_1256(i1,i2)%c(1,iut)=czero
            l3_1256(i1,i2)%d(2,iut)=czero
          enddo
        enddo
      enddo
  
  
      do i1=1,2
      do i2=1,2
* TTW -- aa=l3_1256(i1,i2)%a,bb=l3_1256(i1,i2)%b,cc=l3_1256(i1,i2)%c,dd=
* l3_1256(i1,i2)%d,a1=l3_12(i1,i2)%a,b1=l3_12(i1,i2)%b,c1=l3_12(i1,i2)%c
* ,d1=l3_12(i1,i2)%d,a2=u312_56%a,b2=u312_56%b,c2=u312_56%c,d2=u312_56%d
* ,prq=p312q,m=rmass(id3),nsum=1
      l3_1256(i1,i2)%c(1,1)=l3_1256(i1,i2)%c(1,1)+rmass(id3)*(l3
     & _12(i1,i2)%a(1,1)*u312_56%d(1,1)+l3_12(i1,i2)%c(1,2)*u312
     & _56%c(2,1))
      l3_1256(i1,i2)%d(1,1)=l3_1256(i1,i2)%d(1,1)+l3_12(i1,i2)%d
     & (1,1)*p312q*u312_56%d(1,1)+l3_12(i1,i2)%b(1,2)*u312_56%c(
     & 2,1)
      l3_1256(i1,i2)%a(1,2)=l3_1256(i1,i2)%a(1,2)+rmass(id3)*(l3
     & _12(i1,i2)%a(1,1)*u312_56%b(1,2)+l3_12(i1,i2)%c(1,2)*u312
     & _56%a(2,2))
      l3_1256(i1,i2)%b(1,2)=l3_1256(i1,i2)%b(1,2)+l3_12(i1,i2)%d
     & (1,1)*p312q*u312_56%b(1,2)+l3_12(i1,i2)%b(1,2)*u312_56%a(
     & 2,2)
      l3_1256(i1,i2)%c(2,1)=l3_1256(i1,i2)%c(2,1)+l3_12(i1,i2)%c
     & (2,1)*p312q*u312_56%d(1,1)+l3_12(i1,i2)%a(2,2)*u312_56%c(
     & 2,1)
      l3_1256(i1,i2)%d(2,1)=l3_1256(i1,i2)%d(2,1)+rmass(id3)*(l3
     & _12(i1,i2)%b(2,1)*u312_56%d(1,1)+l3_12(i1,i2)%d(2,2)*u312
     & _56%c(2,1))
      l3_1256(i1,i2)%a(2,2)=l3_1256(i1,i2)%a(2,2)+l3_12(i1,i2)%c
     & (2,1)*p312q*u312_56%b(1,2)+l3_12(i1,i2)%a(2,2)*u312_56%a(
     & 2,2)
      l3_1256(i1,i2)%b(2,2)=l3_1256(i1,i2)%b(2,2)+rmass(id3)*(l3
     & _12(i1,i2)%b(2,1)*u312_56%b(1,2)+l3_12(i1,i2)%d(2,2)*u312
     & _56%a(2,2))
      end do
      end do
  
        endif
      endif
  
  
* compute all single insertions of the type lwi_ (i=5,7)                
*    to a Wline                                                         
*                                                                       
*        i __                                                           
*            |_W__(mu)                                                  
*            |                                                          
  
*zw                                                                     
  
      if(iresw(5).eq.1) then
        do mu=0,3
          lw5_612(mu)%a(2)=czero
          lw5_612(mu)%c(2)=czero
        enddo
      else
* quqd -- p=p5,q=p612
      quqd=p5(0)*p612(0)-p5(1)*p612(1)-p5(2)*p612(2)-p5(3)*p612(
     & 3)
* TWL0 -- qu=p5,qd=p612,v=0,a=lw5_612(0)%a,c=lw5_612(0)%c,cl=wcl,nsum=0
      eps_0=-p5(2)*p612(3)+p612(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p612(0)+p612k0*p5(0)
      lw5_612(0)%a(2)=wcl*(auxa-ceps_0)
      lw5_612(0)%c(2)=wcl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p612,v=1,a=lw5_612(1)%a,c=lw5_612(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p5k0*p612(1)+p612k0*p5(1)
      lw5_612(1)%a(2)=wcl*(auxa-ceps_0)
      lw5_612(1)%c(2)=wcl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p612,v=2,a=lw5_612(2)%a,c=lw5_612(2)%c,cl=wcl,nsum=0
      eps_0=-p5k0*p612(3)+p612k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p612(2)+p612k0*p5(2)
      lw5_612(2)%a(2)=wcl*(auxa-ceps_0)
      lw5_612(2)%c(2)=-wcl*p5k0
* TWL0 -- qu=p5,qd=p612,v=3,a=lw5_612(3)%a,c=lw5_612(3)%c,cl=wcl,nsum=0
      eps_0=p5k0*p612(2)-p612k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p612(3)+p612k0*p5(3)
      lw5_612(3)%a(2)=wcl*(auxa-ceps_0)
      lw5_612(3)%c(2)=wcl*ceps_1
      endif
*zwend                                                                  
*zw                                                                     
  
      if(iresw(5).eq.1) then
        do mu=0,3
          lw5_634(mu)%a(2)=czero
          lw5_634(mu)%c(2)=czero
        enddo
      else
* quqd -- p=p5,q=p634
      quqd=p5(0)*p634(0)-p5(1)*p634(1)-p5(2)*p634(2)-p5(3)*p634(
     & 3)
* TWL0 -- qu=p5,qd=p634,v=0,a=lw5_634(0)%a,c=lw5_634(0)%c,cl=wcl,nsum=0
      eps_0=-p5(2)*p634(3)+p634(2)*p5(3)
      ceps_0=eps_0*cim
      ceps_1=p5(3)*cim
      auxa=-quqd+p5k0*p634(0)+p634k0*p5(0)
      lw5_634(0)%a(2)=wcl*(auxa-ceps_0)
      lw5_634(0)%c(2)=wcl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p634,v=1,a=lw5_634(1)%a,c=lw5_634(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p5k0*p634(1)+p634k0*p5(1)
      lw5_634(1)%a(2)=wcl*(auxa-ceps_0)
      lw5_634(1)%c(2)=wcl*(-p5(2)+ceps_1)
* TWL0 -- qu=p5,qd=p634,v=2,a=lw5_634(2)%a,c=lw5_634(2)%c,cl=wcl,nsum=0
      eps_0=-p5k0*p634(3)+p634k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p634(2)+p634k0*p5(2)
      lw5_634(2)%a(2)=wcl*(auxa-ceps_0)
      lw5_634(2)%c(2)=-wcl*p5k0
* TWL0 -- qu=p5,qd=p634,v=3,a=lw5_634(3)%a,c=lw5_634(3)%c,cl=wcl,nsum=0
      eps_0=p5k0*p634(2)-p634k0*p5(2)
      ceps_0=eps_0*cim
      ceps_1=p5k0*cim
      auxa=p5k0*p634(3)+p634k0*p5(3)
      lw5_634(3)%a(2)=wcl*(auxa-ceps_0)
      lw5_634(3)%c(2)=wcl*ceps_1
      endif
*zwend                                                                  
*zw                                                                     
  
      if(iresw(7).eq.1) then
        do mu=0,3
          lw7_812(mu)%a(2)=czero
          lw7_812(mu)%c(2)=czero
        enddo
      else
* quqd -- p=p7,q=p812
      quqd=p7(0)*p812(0)-p7(1)*p812(1)-p7(2)*p812(2)-p7(3)*p812(
     & 3)
* TWL0 -- qu=p7,qd=p812,v=0,a=lw7_812(0)%a,c=lw7_812(0)%c,cl=wcl,nsum=0
      eps_0=-p7(2)*p812(3)+p812(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p812(0)+p812k0*p7(0)
      lw7_812(0)%a(2)=wcl*(auxa-ceps_0)
      lw7_812(0)%c(2)=wcl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=1,a=lw7_812(1)%a,c=lw7_812(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p7k0*p812(1)+p812k0*p7(1)
      lw7_812(1)%a(2)=wcl*(auxa-ceps_0)
      lw7_812(1)%c(2)=wcl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p812,v=2,a=lw7_812(2)%a,c=lw7_812(2)%c,cl=wcl,nsum=0
      eps_0=-p7k0*p812(3)+p812k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p812(2)+p812k0*p7(2)
      lw7_812(2)%a(2)=wcl*(auxa-ceps_0)
      lw7_812(2)%c(2)=-wcl*p7k0
* TWL0 -- qu=p7,qd=p812,v=3,a=lw7_812(3)%a,c=lw7_812(3)%c,cl=wcl,nsum=0
      eps_0=p7k0*p812(2)-p812k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p812(3)+p812k0*p7(3)
      lw7_812(3)%a(2)=wcl*(auxa-ceps_0)
      lw7_812(3)%c(2)=wcl*ceps_1
      endif
*zwend                                                                  
*zw                                                                     
  
      if(iresw(7).eq.1) then
        do mu=0,3
          lw7_834(mu)%a(2)=czero
          lw7_834(mu)%c(2)=czero
        enddo
      else
* quqd -- p=p7,q=p834
      quqd=p7(0)*p834(0)-p7(1)*p834(1)-p7(2)*p834(2)-p7(3)*p834(
     & 3)
* TWL0 -- qu=p7,qd=p834,v=0,a=lw7_834(0)%a,c=lw7_834(0)%c,cl=wcl,nsum=0
      eps_0=-p7(2)*p834(3)+p834(2)*p7(3)
      ceps_0=eps_0*cim
      ceps_1=p7(3)*cim
      auxa=-quqd+p7k0*p834(0)+p834k0*p7(0)
      lw7_834(0)%a(2)=wcl*(auxa-ceps_0)
      lw7_834(0)%c(2)=wcl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p834,v=1,a=lw7_834(1)%a,c=lw7_834(1)%c,cl=wcl,nsum=0
      auxa=-quqd+p7k0*p834(1)+p834k0*p7(1)
      lw7_834(1)%a(2)=wcl*(auxa-ceps_0)
      lw7_834(1)%c(2)=wcl*(-p7(2)+ceps_1)
* TWL0 -- qu=p7,qd=p834,v=2,a=lw7_834(2)%a,c=lw7_834(2)%c,cl=wcl,nsum=0
      eps_0=-p7k0*p834(3)+p834k0*p7(3)
      ceps_0=eps_0*cim
      auxa=p7k0*p834(2)+p834k0*p7(2)
      lw7_834(2)%a(2)=wcl*(auxa-ceps_0)
      lw7_834(2)%c(2)=-wcl*p7k0
* TWL0 -- qu=p7,qd=p834,v=3,a=lw7_834(3)%a,c=lw7_834(3)%c,cl=wcl,nsum=0
      eps_0=p7k0*p834(2)-p834k0*p7(2)
      ceps_0=eps_0*cim
      ceps_1=p7k0*cim
      auxa=p7k0*p834(3)+p834k0*p7(3)
      lw7_834(3)%a(2)=wcl*(auxa-ceps_0)
      lw7_834(3)%c(2)=wcl*ceps_1
      endif
*zwend                                                                  
  
  
* compute all single insertions of the type rwi_ (i=6,8)                
*    to a Wline                                                         
*                                                                       
*            |_W__(mu)                                                  
*        i __|                                                          
*                                                                       
  
* quqd -- p=p512,q=p6
      quqd=p512(0)*p6(0)-p512(1)*p6(1)-p512(2)*p6(2)-p512(3)*p6(
     & 3)
* TWR0 -- qu=p512,qd=p6,v=0,a=rw6_512(0)%a,b=rw6_512(0)%b,cl=wcl,nsum=0
      eps_0=-p512(2)*p6(3)+p6(2)*p512(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p512k0*p6(0)+p6k0*p512(0)
      rw6_512(0)%a(2)=wcl*(auxa-ceps_0)
      rw6_512(0)%b(1)=-wcl*(p6(2)+ceps_2)
* TWR0 -- qu=p512,qd=p6,v=1,a=rw6_512(1)%a,b=rw6_512(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p512k0*p6(1)+p6k0*p512(1)
      rw6_512(1)%a(2)=wcl*(auxa-ceps_0)
      rw6_512(1)%b(1)=-wcl*(p6(2)+ceps_2)
* TWR0 -- qu=p512,qd=p6,v=2,a=rw6_512(2)%a,b=rw6_512(2)%b,cl=wcl,nsum=0
      eps_0=-p512k0*p6(3)+p6k0*p512(3)
      ceps_0=eps_0*cim
      auxa=p512k0*p6(2)+p6k0*p512(2)
      rw6_512(2)%a(2)=wcl*(auxa-ceps_0)
      rw6_512(2)%b(1)=-wcl*p6k0
* TWR0 -- qu=p512,qd=p6,v=3,a=rw6_512(3)%a,b=rw6_512(3)%b,cl=wcl,nsum=0
      eps_0=p512k0*p6(2)-p6k0*p512(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p512k0*p6(3)+p6k0*p512(3)
      rw6_512(3)%a(2)=wcl*(auxa-ceps_0)
      rw6_512(3)%b(1)=-wcl*ceps_2
* quqd -- p=p534,q=p6
      quqd=p534(0)*p6(0)-p534(1)*p6(1)-p534(2)*p6(2)-p534(3)*p6(
     & 3)
* TWR0 -- qu=p534,qd=p6,v=0,a=rw6_534(0)%a,b=rw6_534(0)%b,cl=wcl,nsum=0
      eps_0=-p534(2)*p6(3)+p6(2)*p534(3)
      ceps_0=eps_0*cim
      ceps_2=p6(3)*cim
      auxa=-quqd+p534k0*p6(0)+p6k0*p534(0)
      rw6_534(0)%a(2)=wcl*(auxa-ceps_0)
      rw6_534(0)%b(1)=-wcl*(p6(2)+ceps_2)
* TWR0 -- qu=p534,qd=p6,v=1,a=rw6_534(1)%a,b=rw6_534(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p534k0*p6(1)+p6k0*p534(1)
      rw6_534(1)%a(2)=wcl*(auxa-ceps_0)
      rw6_534(1)%b(1)=-wcl*(p6(2)+ceps_2)
* TWR0 -- qu=p534,qd=p6,v=2,a=rw6_534(2)%a,b=rw6_534(2)%b,cl=wcl,nsum=0
      eps_0=-p534k0*p6(3)+p6k0*p534(3)
      ceps_0=eps_0*cim
      auxa=p534k0*p6(2)+p6k0*p534(2)
      rw6_534(2)%a(2)=wcl*(auxa-ceps_0)
      rw6_534(2)%b(1)=-wcl*p6k0
* TWR0 -- qu=p534,qd=p6,v=3,a=rw6_534(3)%a,b=rw6_534(3)%b,cl=wcl,nsum=0
      eps_0=p534k0*p6(2)-p6k0*p534(2)
      ceps_0=eps_0*cim
      ceps_2=p6k0*cim
      auxa=p534k0*p6(3)+p6k0*p534(3)
      rw6_534(3)%a(2)=wcl*(auxa-ceps_0)
      rw6_534(3)%b(1)=-wcl*ceps_2
* quqd -- p=p712,q=p8
      quqd=p712(0)*p8(0)-p712(1)*p8(1)-p712(2)*p8(2)-p712(3)*p8(
     & 3)
* TWR0 -- qu=p712,qd=p8,v=0,a=rw8_712(0)%a,b=rw8_712(0)%b,cl=wcl,nsum=0
      eps_0=-p712(2)*p8(3)+p8(2)*p712(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p712k0*p8(0)+p8k0*p712(0)
      rw8_712(0)%a(2)=wcl*(auxa-ceps_0)
      rw8_712(0)%b(1)=-wcl*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=1,a=rw8_712(1)%a,b=rw8_712(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p712k0*p8(1)+p8k0*p712(1)
      rw8_712(1)%a(2)=wcl*(auxa-ceps_0)
      rw8_712(1)%b(1)=-wcl*(p8(2)+ceps_2)
* TWR0 -- qu=p712,qd=p8,v=2,a=rw8_712(2)%a,b=rw8_712(2)%b,cl=wcl,nsum=0
      eps_0=-p712k0*p8(3)+p8k0*p712(3)
      ceps_0=eps_0*cim
      auxa=p712k0*p8(2)+p8k0*p712(2)
      rw8_712(2)%a(2)=wcl*(auxa-ceps_0)
      rw8_712(2)%b(1)=-wcl*p8k0
* TWR0 -- qu=p712,qd=p8,v=3,a=rw8_712(3)%a,b=rw8_712(3)%b,cl=wcl,nsum=0
      eps_0=p712k0*p8(2)-p8k0*p712(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p712k0*p8(3)+p8k0*p712(3)
      rw8_712(3)%a(2)=wcl*(auxa-ceps_0)
      rw8_712(3)%b(1)=-wcl*ceps_2
* quqd -- p=p734,q=p8
      quqd=p734(0)*p8(0)-p734(1)*p8(1)-p734(2)*p8(2)-p734(3)*p8(
     & 3)
* TWR0 -- qu=p734,qd=p8,v=0,a=rw8_734(0)%a,b=rw8_734(0)%b,cl=wcl,nsum=0
      eps_0=-p734(2)*p8(3)+p8(2)*p734(3)
      ceps_0=eps_0*cim
      ceps_2=p8(3)*cim
      auxa=-quqd+p734k0*p8(0)+p8k0*p734(0)
      rw8_734(0)%a(2)=wcl*(auxa-ceps_0)
      rw8_734(0)%b(1)=-wcl*(p8(2)+ceps_2)
* TWR0 -- qu=p734,qd=p8,v=1,a=rw8_734(1)%a,b=rw8_734(1)%b,cl=wcl,nsum=0
      auxa=-quqd+p734k0*p8(1)+p8k0*p734(1)
      rw8_734(1)%a(2)=wcl*(auxa-ceps_0)
      rw8_734(1)%b(1)=-wcl*(p8(2)+ceps_2)
* TWR0 -- qu=p734,qd=p8,v=2,a=rw8_734(2)%a,b=rw8_734(2)%b,cl=wcl,nsum=0
      eps_0=-p734k0*p8(3)+p8k0*p734(3)
      ceps_0=eps_0*cim
      auxa=p734k0*p8(2)+p8k0*p734(2)
      rw8_734(2)%a(2)=wcl*(auxa-ceps_0)
      rw8_734(2)%b(1)=-wcl*p8k0
* TWR0 -- qu=p734,qd=p8,v=3,a=rw8_734(3)%a,b=rw8_734(3)%b,cl=wcl,nsum=0
      eps_0=p734k0*p8(2)-p8k0*p734(2)
      ceps_0=eps_0*cim
      ceps_2=p8k0*cim
      auxa=p734k0*p8(3)+p8k0*p734(3)
      rw8_734(3)%a(2)=wcl*(auxa-ceps_0)
      rw8_734(3)%b(1)=-wcl*ceps_2
  
  
* compute all single insertions of the type lwi_ (i=1,3)                
*    to a Zline                                                         
*                                                                       
*        i __                                                           
*            |_W__(mu)                                                  
*            |                                                          
*                                                                       
*CAREFUL: also the second attachment is a W                             
* so (as 56=+ and 78=-) if i is up one can have only lwi_j56,           
*    if i is down only lwi_j78                                          
  
*zw                                                                     
  
      if(iresz(1).eq.1) then
        if (iup(id1).eq.1) then
          do mu=0,3
            lw1_256(mu)%a(2,2)=czero
            lw1_256(mu)%b(1,2)=czero
            lw1_256(mu)%c(2,1)=czero
          enddo
        else
          do mu=0,3
            lw1_278(mu)%a(2,2)=czero
            lw1_278(mu)%b(1,2)=czero
            lw1_278(mu)%c(2,1)=czero
          enddo
        endif
      else
  
        if (iup(id1).eq.1) then
* quqd -- p=p1,q=p256
      quqd=p1(0)*p256(0)-p1(1)*p256(1)-p1(2)*p256(2)-p1(3)*p256(
     & 3)
* TW -- qu=p1,qd=p256,v=0,a=lw1_256(0)%a,b=lw1_256(0)%b,c=lw1_256(0)%c,d
* =lw1_256(0)%d,cl=wcl,nsum=0
      eps_0=-p1(2)*p256(3)+p256(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      ceps_2=p256(3)*cim
      auxa=-quqd+p1k0*p256(0)+p256k0*p1(0)
      lw1_256(0)%a(2,2)=wcl*(auxa-ceps_0)
      lw1_256(0)%b(1,2)=-wcl*(p256(2)+ceps_2)
      lw1_256(0)%c(2,1)=wcl*(-p1(2)+ceps_1)
      lw1_256(0)%d(1,1)=wcl
* TW -- qu=p1,qd=p256,v=1,a=lw1_256(1)%a,b=lw1_256(1)%b,c=lw1_256(1)%c,d
* =lw1_256(1)%d,cl=wcl,nsum=0
      auxa=-quqd+p1k0*p256(1)+p256k0*p1(1)
      lw1_256(1)%a(2,2)=wcl*(auxa-ceps_0)
      lw1_256(1)%b(1,2)=-wcl*(p256(2)+ceps_2)
      lw1_256(1)%c(2,1)=wcl*(-p1(2)+ceps_1)
      lw1_256(1)%d(1,1)=wcl
* TW -- qu=p1,qd=p256,v=2,a=lw1_256(2)%a,b=lw1_256(2)%b,c=lw1_256(2)%c,d
* =lw1_256(2)%d,cl=wcl,nsum=0
      eps_0=-p1k0*p256(3)+p256k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p256(2)+p256k0*p1(2)
      lw1_256(2)%a(2,2)=wcl*(auxa-ceps_0)
      lw1_256(2)%b(1,2)=-wcl*p256k0
      lw1_256(2)%c(2,1)=-wcl*p1k0
* TW -- qu=p1,qd=p256,v=3,a=lw1_256(3)%a,b=lw1_256(3)%b,c=lw1_256(3)%c,d
* =lw1_256(3)%d,cl=wcl,nsum=0
      eps_0=p1k0*p256(2)-p256k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      ceps_2=p256k0*cim
      auxa=+p1k0*p256(3)+p256k0*p1(3)
      lw1_256(3)%a(2,2)=wcl*(auxa-ceps_0)
      lw1_256(3)%b(1,2)=-wcl*ceps_2
      lw1_256(3)%c(2,1)=wcl*ceps_1
        else
* quqd -- p=p1,q=p278
      quqd=p1(0)*p278(0)-p1(1)*p278(1)-p1(2)*p278(2)-p1(3)*p278(
     & 3)
* TW -- qu=p1,qd=p278,v=0,a=lw1_278(0)%a,b=lw1_278(0)%b,c=lw1_278(0)%c,d
* =lw1_278(0)%d,cl=wcl,nsum=0
      eps_0=-p1(2)*p278(3)+p278(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      ceps_2=p278(3)*cim
      auxa=-quqd+p1k0*p278(0)+p278k0*p1(0)
      lw1_278(0)%a(2,2)=wcl*(auxa-ceps_0)
      lw1_278(0)%b(1,2)=-wcl*(p278(2)+ceps_2)
      lw1_278(0)%c(2,1)=wcl*(-p1(2)+ceps_1)
      lw1_278(0)%d(1,1)=wcl
* TW -- qu=p1,qd=p278,v=1,a=lw1_278(1)%a,b=lw1_278(1)%b,c=lw1_278(1)%c,d
* =lw1_278(1)%d,cl=wcl,nsum=0
      auxa=-quqd+p1k0*p278(1)+p278k0*p1(1)
      lw1_278(1)%a(2,2)=wcl*(auxa-ceps_0)
      lw1_278(1)%b(1,2)=-wcl*(p278(2)+ceps_2)
      lw1_278(1)%c(2,1)=wcl*(-p1(2)+ceps_1)
      lw1_278(1)%d(1,1)=wcl
* TW -- qu=p1,qd=p278,v=2,a=lw1_278(2)%a,b=lw1_278(2)%b,c=lw1_278(2)%c,d
* =lw1_278(2)%d,cl=wcl,nsum=0
      eps_0=-p1k0*p278(3)+p278k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p278(2)+p278k0*p1(2)
      lw1_278(2)%a(2,2)=wcl*(auxa-ceps_0)
      lw1_278(2)%b(1,2)=-wcl*p278k0
      lw1_278(2)%c(2,1)=-wcl*p1k0
* TW -- qu=p1,qd=p278,v=3,a=lw1_278(3)%a,b=lw1_278(3)%b,c=lw1_278(3)%c,d
* =lw1_278(3)%d,cl=wcl,nsum=0
      eps_0=p1k0*p278(2)-p278k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      ceps_2=p278k0*cim
      auxa=+p1k0*p278(3)+p278k0*p1(3)
      lw1_278(3)%a(2,2)=wcl*(auxa-ceps_0)
      lw1_278(3)%b(1,2)=-wcl*ceps_2
      lw1_278(3)%c(2,1)=wcl*ceps_1
        endif
      endif
*zwend                                                                  
*zw                                                                     
  
      if(iresz(3).eq.1) then
        if (iup(id3).eq.1) then
          do mu=0,3
            lw3_456(mu)%a(2,2)=czero
            lw3_456(mu)%b(1,2)=czero
            lw3_456(mu)%c(2,1)=czero
          enddo
        else
          do mu=0,3
            lw3_478(mu)%a(2,2)=czero
            lw3_478(mu)%b(1,2)=czero
            lw3_478(mu)%c(2,1)=czero
          enddo
        endif
      else
  
        if (iup(id3).eq.1) then
* quqd -- p=p3,q=p456
      quqd=p3(0)*p456(0)-p3(1)*p456(1)-p3(2)*p456(2)-p3(3)*p456(
     & 3)
* TW -- qu=p3,qd=p456,v=0,a=lw3_456(0)%a,b=lw3_456(0)%b,c=lw3_456(0)%c,d
* =lw3_456(0)%d,cl=wcl,nsum=0
      eps_0=-p3(2)*p456(3)+p456(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      ceps_2=p456(3)*cim
      auxa=-quqd+p3k0*p456(0)+p456k0*p3(0)
      lw3_456(0)%a(2,2)=wcl*(auxa-ceps_0)
      lw3_456(0)%b(1,2)=-wcl*(p456(2)+ceps_2)
      lw3_456(0)%c(2,1)=wcl*(-p3(2)+ceps_1)
      lw3_456(0)%d(1,1)=wcl
* TW -- qu=p3,qd=p456,v=1,a=lw3_456(1)%a,b=lw3_456(1)%b,c=lw3_456(1)%c,d
* =lw3_456(1)%d,cl=wcl,nsum=0
      auxa=-quqd+p3k0*p456(1)+p456k0*p3(1)
      lw3_456(1)%a(2,2)=wcl*(auxa-ceps_0)
      lw3_456(1)%b(1,2)=-wcl*(p456(2)+ceps_2)
      lw3_456(1)%c(2,1)=wcl*(-p3(2)+ceps_1)
      lw3_456(1)%d(1,1)=wcl
* TW -- qu=p3,qd=p456,v=2,a=lw3_456(2)%a,b=lw3_456(2)%b,c=lw3_456(2)%c,d
* =lw3_456(2)%d,cl=wcl,nsum=0
      eps_0=-p3k0*p456(3)+p456k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p456(2)+p456k0*p3(2)
      lw3_456(2)%a(2,2)=wcl*(auxa-ceps_0)
      lw3_456(2)%b(1,2)=-wcl*p456k0
      lw3_456(2)%c(2,1)=-wcl*p3k0
* TW -- qu=p3,qd=p456,v=3,a=lw3_456(3)%a,b=lw3_456(3)%b,c=lw3_456(3)%c,d
* =lw3_456(3)%d,cl=wcl,nsum=0
      eps_0=p3k0*p456(2)-p456k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      ceps_2=p456k0*cim
      auxa=+p3k0*p456(3)+p456k0*p3(3)
      lw3_456(3)%a(2,2)=wcl*(auxa-ceps_0)
      lw3_456(3)%b(1,2)=-wcl*ceps_2
      lw3_456(3)%c(2,1)=wcl*ceps_1
        else
* quqd -- p=p3,q=p478
      quqd=p3(0)*p478(0)-p3(1)*p478(1)-p3(2)*p478(2)-p3(3)*p478(
     & 3)
* TW -- qu=p3,qd=p478,v=0,a=lw3_478(0)%a,b=lw3_478(0)%b,c=lw3_478(0)%c,d
* =lw3_478(0)%d,cl=wcl,nsum=0
      eps_0=-p3(2)*p478(3)+p478(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      ceps_2=p478(3)*cim
      auxa=-quqd+p3k0*p478(0)+p478k0*p3(0)
      lw3_478(0)%a(2,2)=wcl*(auxa-ceps_0)
      lw3_478(0)%b(1,2)=-wcl*(p478(2)+ceps_2)
      lw3_478(0)%c(2,1)=wcl*(-p3(2)+ceps_1)
      lw3_478(0)%d(1,1)=wcl
* TW -- qu=p3,qd=p478,v=1,a=lw3_478(1)%a,b=lw3_478(1)%b,c=lw3_478(1)%c,d
* =lw3_478(1)%d,cl=wcl,nsum=0
      auxa=-quqd+p3k0*p478(1)+p478k0*p3(1)
      lw3_478(1)%a(2,2)=wcl*(auxa-ceps_0)
      lw3_478(1)%b(1,2)=-wcl*(p478(2)+ceps_2)
      lw3_478(1)%c(2,1)=wcl*(-p3(2)+ceps_1)
      lw3_478(1)%d(1,1)=wcl
* TW -- qu=p3,qd=p478,v=2,a=lw3_478(2)%a,b=lw3_478(2)%b,c=lw3_478(2)%c,d
* =lw3_478(2)%d,cl=wcl,nsum=0
      eps_0=-p3k0*p478(3)+p478k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p478(2)+p478k0*p3(2)
      lw3_478(2)%a(2,2)=wcl*(auxa-ceps_0)
      lw3_478(2)%b(1,2)=-wcl*p478k0
      lw3_478(2)%c(2,1)=-wcl*p3k0
* TW -- qu=p3,qd=p478,v=3,a=lw3_478(3)%a,b=lw3_478(3)%b,c=lw3_478(3)%c,d
* =lw3_478(3)%d,cl=wcl,nsum=0
      eps_0=p3k0*p478(2)-p478k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      ceps_2=p478k0*cim
      auxa=+p3k0*p478(3)+p478k0*p3(3)
      lw3_478(3)%a(2,2)=wcl*(auxa-ceps_0)
      lw3_478(3)%b(1,2)=-wcl*ceps_2
      lw3_478(3)%c(2,1)=wcl*ceps_1
        endif
      endif
*zwend                                                                  
  
* compute all single insertions of the type rwi_ (i=2,4)                
*    to a Zline                                                         
*                                                                       
*            |_W__(mu)                                                  
*        i __|                                                          
*                                                                       
*CAREFUL: also the upper attachment is a W                              
* so (as 56=+ and 78=-) if i is up one can have only rwi_j78,           
*    if i is down only rwi_j56                                          
  
      if (iup(id2).eq.1) then
* quqd -- p=p178,q=p2
      quqd=p178(0)*p2(0)-p178(1)*p2(1)-p178(2)*p2(2)-p178(3)*p2(
     & 3)
* TW -- qu=p178,qd=p2,v=0,a=rw2_178(0)%a,b=rw2_178(0)%b,c=rw2_178(0)%c,d
* =rw2_178(0)%d,cl=wcl,nsum=0
      eps_0=-p178(2)*p2(3)+p2(2)*p178(3)
      ceps_0=eps_0*cim
      ceps_1=p178(3)*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p178k0*p2(0)+p2k0*p178(0)
      rw2_178(0)%a(2,2)=wcl*(auxa-ceps_0)
      rw2_178(0)%b(1,2)=-wcl*(p2(2)+ceps_2)
      rw2_178(0)%c(2,1)=wcl*(-p178(2)+ceps_1)
      rw2_178(0)%d(1,1)=wcl
* TW -- qu=p178,qd=p2,v=1,a=rw2_178(1)%a,b=rw2_178(1)%b,c=rw2_178(1)%c,d
* =rw2_178(1)%d,cl=wcl,nsum=0
      auxa=-quqd+p178k0*p2(1)+p2k0*p178(1)
      rw2_178(1)%a(2,2)=wcl*(auxa-ceps_0)
      rw2_178(1)%b(1,2)=-wcl*(p2(2)+ceps_2)
      rw2_178(1)%c(2,1)=wcl*(-p178(2)+ceps_1)
      rw2_178(1)%d(1,1)=wcl
* TW -- qu=p178,qd=p2,v=2,a=rw2_178(2)%a,b=rw2_178(2)%b,c=rw2_178(2)%c,d
* =rw2_178(2)%d,cl=wcl,nsum=0
      eps_0=-p178k0*p2(3)+p2k0*p178(3)
      ceps_0=eps_0*cim
      auxa=p178k0*p2(2)+p2k0*p178(2)
      rw2_178(2)%a(2,2)=wcl*(auxa-ceps_0)
      rw2_178(2)%b(1,2)=-wcl*p2k0
      rw2_178(2)%c(2,1)=-wcl*p178k0
* TW -- qu=p178,qd=p2,v=3,a=rw2_178(3)%a,b=rw2_178(3)%b,c=rw2_178(3)%c,d
* =rw2_178(3)%d,cl=wcl,nsum=0
      eps_0=p178k0*p2(2)-p2k0*p178(2)
      ceps_0=eps_0*cim
      ceps_1=p178k0*cim
      ceps_2=p2k0*cim
      auxa=+p178k0*p2(3)+p2k0*p178(3)
      rw2_178(3)%a(2,2)=wcl*(auxa-ceps_0)
      rw2_178(3)%b(1,2)=-wcl*ceps_2
      rw2_178(3)%c(2,1)=wcl*ceps_1
      else
* quqd -- p=p156,q=p2
      quqd=p156(0)*p2(0)-p156(1)*p2(1)-p156(2)*p2(2)-p156(3)*p2(
     & 3)
* TW -- qu=p156,qd=p2,v=0,a=rw2_156(0)%a,b=rw2_156(0)%b,c=rw2_156(0)%c,d
* =rw2_156(0)%d,cl=wcl,nsum=0
      eps_0=-p156(2)*p2(3)+p2(2)*p156(3)
      ceps_0=eps_0*cim
      ceps_1=p156(3)*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p156k0*p2(0)+p2k0*p156(0)
      rw2_156(0)%a(2,2)=wcl*(auxa-ceps_0)
      rw2_156(0)%b(1,2)=-wcl*(p2(2)+ceps_2)
      rw2_156(0)%c(2,1)=wcl*(-p156(2)+ceps_1)
      rw2_156(0)%d(1,1)=wcl
* TW -- qu=p156,qd=p2,v=1,a=rw2_156(1)%a,b=rw2_156(1)%b,c=rw2_156(1)%c,d
* =rw2_156(1)%d,cl=wcl,nsum=0
      auxa=-quqd+p156k0*p2(1)+p2k0*p156(1)
      rw2_156(1)%a(2,2)=wcl*(auxa-ceps_0)
      rw2_156(1)%b(1,2)=-wcl*(p2(2)+ceps_2)
      rw2_156(1)%c(2,1)=wcl*(-p156(2)+ceps_1)
      rw2_156(1)%d(1,1)=wcl
* TW -- qu=p156,qd=p2,v=2,a=rw2_156(2)%a,b=rw2_156(2)%b,c=rw2_156(2)%c,d
* =rw2_156(2)%d,cl=wcl,nsum=0
      eps_0=-p156k0*p2(3)+p2k0*p156(3)
      ceps_0=eps_0*cim
      auxa=p156k0*p2(2)+p2k0*p156(2)
      rw2_156(2)%a(2,2)=wcl*(auxa-ceps_0)
      rw2_156(2)%b(1,2)=-wcl*p2k0
      rw2_156(2)%c(2,1)=-wcl*p156k0
* TW -- qu=p156,qd=p2,v=3,a=rw2_156(3)%a,b=rw2_156(3)%b,c=rw2_156(3)%c,d
* =rw2_156(3)%d,cl=wcl,nsum=0
      eps_0=p156k0*p2(2)-p2k0*p156(2)
      ceps_0=eps_0*cim
      ceps_1=p156k0*cim
      ceps_2=p2k0*cim
      auxa=+p156k0*p2(3)+p2k0*p156(3)
      rw2_156(3)%a(2,2)=wcl*(auxa-ceps_0)
      rw2_156(3)%b(1,2)=-wcl*ceps_2
      rw2_156(3)%c(2,1)=wcl*ceps_1
      endif
      if (iup(id4).eq.1) then
* quqd -- p=p378,q=p4
      quqd=p378(0)*p4(0)-p378(1)*p4(1)-p378(2)*p4(2)-p378(3)*p4(
     & 3)
* TW -- qu=p378,qd=p4,v=0,a=rw4_378(0)%a,b=rw4_378(0)%b,c=rw4_378(0)%c,d
* =rw4_378(0)%d,cl=wcl,nsum=0
      eps_0=-p378(2)*p4(3)+p4(2)*p378(3)
      ceps_0=eps_0*cim
      ceps_1=p378(3)*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p378k0*p4(0)+p4k0*p378(0)
      rw4_378(0)%a(2,2)=wcl*(auxa-ceps_0)
      rw4_378(0)%b(1,2)=-wcl*(p4(2)+ceps_2)
      rw4_378(0)%c(2,1)=wcl*(-p378(2)+ceps_1)
      rw4_378(0)%d(1,1)=wcl
* TW -- qu=p378,qd=p4,v=1,a=rw4_378(1)%a,b=rw4_378(1)%b,c=rw4_378(1)%c,d
* =rw4_378(1)%d,cl=wcl,nsum=0
      auxa=-quqd+p378k0*p4(1)+p4k0*p378(1)
      rw4_378(1)%a(2,2)=wcl*(auxa-ceps_0)
      rw4_378(1)%b(1,2)=-wcl*(p4(2)+ceps_2)
      rw4_378(1)%c(2,1)=wcl*(-p378(2)+ceps_1)
      rw4_378(1)%d(1,1)=wcl
* TW -- qu=p378,qd=p4,v=2,a=rw4_378(2)%a,b=rw4_378(2)%b,c=rw4_378(2)%c,d
* =rw4_378(2)%d,cl=wcl,nsum=0
      eps_0=-p378k0*p4(3)+p4k0*p378(3)
      ceps_0=eps_0*cim
      auxa=p378k0*p4(2)+p4k0*p378(2)
      rw4_378(2)%a(2,2)=wcl*(auxa-ceps_0)
      rw4_378(2)%b(1,2)=-wcl*p4k0
      rw4_378(2)%c(2,1)=-wcl*p378k0
* TW -- qu=p378,qd=p4,v=3,a=rw4_378(3)%a,b=rw4_378(3)%b,c=rw4_378(3)%c,d
* =rw4_378(3)%d,cl=wcl,nsum=0
      eps_0=p378k0*p4(2)-p4k0*p378(2)
      ceps_0=eps_0*cim
      ceps_1=p378k0*cim
      ceps_2=p4k0*cim
      auxa=+p378k0*p4(3)+p4k0*p378(3)
      rw4_378(3)%a(2,2)=wcl*(auxa-ceps_0)
      rw4_378(3)%b(1,2)=-wcl*ceps_2
      rw4_378(3)%c(2,1)=wcl*ceps_1
      else
* quqd -- p=p356,q=p4
      quqd=p356(0)*p4(0)-p356(1)*p4(1)-p356(2)*p4(2)-p356(3)*p4(
     & 3)
* TW -- qu=p356,qd=p4,v=0,a=rw4_356(0)%a,b=rw4_356(0)%b,c=rw4_356(0)%c,d
* =rw4_356(0)%d,cl=wcl,nsum=0
      eps_0=-p356(2)*p4(3)+p4(2)*p356(3)
      ceps_0=eps_0*cim
      ceps_1=p356(3)*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p356k0*p4(0)+p4k0*p356(0)
      rw4_356(0)%a(2,2)=wcl*(auxa-ceps_0)
      rw4_356(0)%b(1,2)=-wcl*(p4(2)+ceps_2)
      rw4_356(0)%c(2,1)=wcl*(-p356(2)+ceps_1)
      rw4_356(0)%d(1,1)=wcl
* TW -- qu=p356,qd=p4,v=1,a=rw4_356(1)%a,b=rw4_356(1)%b,c=rw4_356(1)%c,d
* =rw4_356(1)%d,cl=wcl,nsum=0
      auxa=-quqd+p356k0*p4(1)+p4k0*p356(1)
      rw4_356(1)%a(2,2)=wcl*(auxa-ceps_0)
      rw4_356(1)%b(1,2)=-wcl*(p4(2)+ceps_2)
      rw4_356(1)%c(2,1)=wcl*(-p356(2)+ceps_1)
      rw4_356(1)%d(1,1)=wcl
* TW -- qu=p356,qd=p4,v=2,a=rw4_356(2)%a,b=rw4_356(2)%b,c=rw4_356(2)%c,d
* =rw4_356(2)%d,cl=wcl,nsum=0
      eps_0=-p356k0*p4(3)+p4k0*p356(3)
      ceps_0=eps_0*cim
      auxa=p356k0*p4(2)+p4k0*p356(2)
      rw4_356(2)%a(2,2)=wcl*(auxa-ceps_0)
      rw4_356(2)%b(1,2)=-wcl*p4k0
      rw4_356(2)%c(2,1)=-wcl*p356k0
* TW -- qu=p356,qd=p4,v=3,a=rw4_356(3)%a,b=rw4_356(3)%b,c=rw4_356(3)%c,d
* =rw4_356(3)%d,cl=wcl,nsum=0
      eps_0=p356k0*p4(2)-p4k0*p356(2)
      ceps_0=eps_0*cim
      ceps_1=p356k0*cim
      ceps_2=p4k0*cim
      auxa=+p356k0*p4(3)+p4k0*p356(3)
      rw4_356(3)%a(2,2)=wcl*(auxa-ceps_0)
      rw4_356(3)%b(1,2)=-wcl*ceps_2
      rw4_356(3)%c(2,1)=wcl*ceps_1
      endif
  
* compute all single insertions of the type lzi_ and lfi_ (i=5,7)       
*    to a Wline                                                         
*                                                                       
*        i __              i __                                         
*            |_Z__(mu)         |__f__(mu)                               
*            |                 |                                        
*                                                                       
* As they are then connected with w insertions, they are calculated     
* as w insertions with the appropriate coupling cl only                 
  
*zw                                                                     
      if(iresw(5).eq.1) then
        do mu=0,3
          lz5_678(mu)%a(2)=czero
          lz5_678(mu)%c(2)=czero
          lf5_678(mu)%a(2)=czero
          lf5_678(mu)%c(2)=czero
        enddo
      else
*zwend                                                                  
  
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
*zw                                                                     
      endif
*zwend                                                                  
*zw                                                                     
      if(iresw(7).eq.1) then
        do mu=0,3
          lz7_856(mu)%a(2)=czero
          lz7_856(mu)%c(2)=czero
          lf7_856(mu)%a(2)=czero
          lf7_856(mu)%c(2)=czero
        enddo
      else
*zwend                                                                  
  
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
*zw                                                                     
      endif
*zwend                                                                  
  
* compute all single insertions of the type rzi_ and rfi_ (i=6,8)       
*    to a Wline                                                         
*                                                                       
*            |_Z__(mu)         |__f__(mu)                               
*        i __|             i __|                                        
*                                                                       
* As they are then connected with w insertions, they are calculated     
* as w insertions with the appropriate coupling cl only                 
  
  
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
  
  
* compute all single insertions of the type lzi_ and lfi_ (i=1,3)       
*    to a Zline                                                         
*                                                                       
*        i __              i __                                         
*            |_Z__(mu)         |__f__(mu)                               
*            |                 |                                        
  
  
*zw                                                                     
      if(iresz(1).eq.1) then
        do m=0,3
          lz1_234(m)%a(1,1)=czero
          lz1_234(m)%a(2,2)=czero
          lz1_234(m)%b(1,2)=czero
          lz1_234(m)%b(2,1)=czero
          lz1_234(m)%c(1,2)=czero
          lz1_234(m)%c(2,1)=czero
          lf1_234(m)%a(1,1)=czero
          lf1_234(m)%a(2,2)=czero
          lf1_234(m)%b(1,2)=czero
          lf1_234(m)%b(2,1)=czero
          lf1_234(m)%c(1,2)=czero
          lf1_234(m)%c(2,1)=czero
        enddo
      else
*zwend                                                                  
  
* quqd -- p=p1,q=p234
      quqd=p1(0)*p234(0)-p1(1)*p234(1)-p1(2)*p234(2)-p1(3)*p234(
     & 3)
* T -- qu=p1,qd=p234,v=0,a=lz1_234(0)%a,b=lz1_234(0)%b,c=lz1_234(0)%c,d=
* lz1_234(0)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      eps_0=-p1(2)*p234(3)+p234(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      ceps_2=p234(3)*cim
      auxa=-quqd+p1k0*p234(0)+p234k0*p1(0)
      lz1_234(0)%a(1,1)=zcr(id1)*(auxa+ceps_0)
      lz1_234(0)%a(2,2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(0)%b(1,2)=-zcl(id1)*(p234(2)+ceps_2)
      lz1_234(0)%b(2,1)=zcr(id1)*(p234(2)-ceps_2)
      lz1_234(0)%c(1,2)=zcr(id1)*(p1(2)+ceps_1)
      lz1_234(0)%c(2,1)=zcl(id1)*(-p1(2)+ceps_1)
      lz1_234(0)%d(1,1)=zcl(id1)
      lz1_234(0)%d(2,2)=zcr(id1)
* T -- qu=p1,qd=p234,v=1,a=lz1_234(1)%a,b=lz1_234(1)%b,c=lz1_234(1)%c,d=
* lz1_234(1)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      auxa=-quqd+p1k0*p234(1)+p234k0*p1(1)
      lz1_234(1)%a(1,1)=zcr(id1)*(auxa+ceps_0)
      lz1_234(1)%a(2,2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(1)%b(1,2)=-zcl(id1)*(p234(2)+ceps_2)
      lz1_234(1)%b(2,1)=zcr(id1)*(p234(2)-ceps_2)
      lz1_234(1)%c(1,2)=zcr(id1)*(p1(2)+ceps_1)
      lz1_234(1)%c(2,1)=zcl(id1)*(-p1(2)+ceps_1)
      lz1_234(1)%d(1,1)=zcl(id1)
      lz1_234(1)%d(2,2)=zcr(id1)
* T -- qu=p1,qd=p234,v=2,a=lz1_234(2)%a,b=lz1_234(2)%b,c=lz1_234(2)%c,d=
* lz1_234(2)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      eps_0=-p1k0*p234(3)+p234k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p234(2)+p234k0*p1(2)
      lz1_234(2)%a(1,1)=zcr(id1)*(auxa+ceps_0)
      lz1_234(2)%a(2,2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(2)%b(1,2)=-zcl(id1)*p234k0
      lz1_234(2)%b(2,1)=zcr(id1)*p234k0
      lz1_234(2)%c(1,2)=zcr(id1)*p1k0
      lz1_234(2)%c(2,1)=-zcl(id1)*p1k0
* T -- qu=p1,qd=p234,v=3,a=lz1_234(3)%a,b=lz1_234(3)%b,c=lz1_234(3)%c,d=
* lz1_234(3)%d,cr=zcr(id1),cl=zcl(id1),nsum=0
      eps_0=p1k0*p234(2)-p234k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      ceps_2=p234k0*cim
      auxa=+p1k0*p234(3)+p234k0*p1(3)
      lz1_234(3)%a(1,1)=zcr(id1)*(auxa+ceps_0)
      lz1_234(3)%a(2,2)=zcl(id1)*(auxa-ceps_0)
      lz1_234(3)%b(1,2)=-zcl(id1)*ceps_2
      lz1_234(3)%b(2,1)=-zcr(id1)*ceps_2
      lz1_234(3)%c(1,2)=zcr(id1)*ceps_1
      lz1_234(3)%c(2,1)=zcl(id1)*ceps_1
  
* T -- qu=p1,qd=p234,v=0,a=lf1_234(0)%a,b=lf1_234(0)%b,c=lf1_234(0)%c,d=
* lf1_234(0)%d,cr=fcr(id1),cl=fcl(id1),nsum=0
      eps_0=-p1(2)*p234(3)+p234(2)*p1(3)
      ceps_0=eps_0*cim
      ceps_1=p1(3)*cim
      ceps_2=p234(3)*cim
      auxa=-quqd+p1k0*p234(0)+p234k0*p1(0)
      lf1_234(0)%a(1,1)=fcr(id1)*(auxa+ceps_0)
      lf1_234(0)%a(2,2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(0)%b(1,2)=-fcl(id1)*(p234(2)+ceps_2)
      lf1_234(0)%b(2,1)=fcr(id1)*(p234(2)-ceps_2)
      lf1_234(0)%c(1,2)=fcr(id1)*(p1(2)+ceps_1)
      lf1_234(0)%c(2,1)=fcl(id1)*(-p1(2)+ceps_1)
      lf1_234(0)%d(1,1)=fcl(id1)
      lf1_234(0)%d(2,2)=fcr(id1)
* T -- qu=p1,qd=p234,v=1,a=lf1_234(1)%a,b=lf1_234(1)%b,c=lf1_234(1)%c,d=
* lf1_234(1)%d,cr=fcr(id1),cl=fcl(id1),nsum=0
      auxa=-quqd+p1k0*p234(1)+p234k0*p1(1)
      lf1_234(1)%a(1,1)=fcr(id1)*(auxa+ceps_0)
      lf1_234(1)%a(2,2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(1)%b(1,2)=-fcl(id1)*(p234(2)+ceps_2)
      lf1_234(1)%b(2,1)=fcr(id1)*(p234(2)-ceps_2)
      lf1_234(1)%c(1,2)=fcr(id1)*(p1(2)+ceps_1)
      lf1_234(1)%c(2,1)=fcl(id1)*(-p1(2)+ceps_1)
      lf1_234(1)%d(1,1)=fcl(id1)
      lf1_234(1)%d(2,2)=fcr(id1)
* T -- qu=p1,qd=p234,v=2,a=lf1_234(2)%a,b=lf1_234(2)%b,c=lf1_234(2)%c,d=
* lf1_234(2)%d,cr=fcr(id1),cl=fcl(id1),nsum=0
      eps_0=-p1k0*p234(3)+p234k0*p1(3)
      ceps_0=eps_0*cim
      auxa=p1k0*p234(2)+p234k0*p1(2)
      lf1_234(2)%a(1,1)=fcr(id1)*(auxa+ceps_0)
      lf1_234(2)%a(2,2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(2)%b(1,2)=-fcl(id1)*p234k0
      lf1_234(2)%b(2,1)=fcr(id1)*p234k0
      lf1_234(2)%c(1,2)=fcr(id1)*p1k0
      lf1_234(2)%c(2,1)=-fcl(id1)*p1k0
* T -- qu=p1,qd=p234,v=3,a=lf1_234(3)%a,b=lf1_234(3)%b,c=lf1_234(3)%c,d=
* lf1_234(3)%d,cr=fcr(id1),cl=fcl(id1),nsum=0
      eps_0=p1k0*p234(2)-p234k0*p1(2)
      ceps_0=eps_0*cim
      ceps_1=p1k0*cim
      ceps_2=p234k0*cim
      auxa=+p1k0*p234(3)+p234k0*p1(3)
      lf1_234(3)%a(1,1)=fcr(id1)*(auxa+ceps_0)
      lf1_234(3)%a(2,2)=fcl(id1)*(auxa-ceps_0)
      lf1_234(3)%b(1,2)=-fcl(id1)*ceps_2
      lf1_234(3)%b(2,1)=-fcr(id1)*ceps_2
      lf1_234(3)%c(1,2)=fcr(id1)*ceps_1
      lf1_234(3)%c(2,1)=fcl(id1)*ceps_1
*zw                                                                     
      endif   !  if(iresz(1).eq.1)
*zwend                                                                  
  
  
*zw                                                                     
      if(iresz(3).eq.1) then
        do m=0,3
          lz3_412(m)%a(1,1)=czero
          lz3_412(m)%a(2,2)=czero
          lz3_412(m)%b(1,2)=czero
          lz3_412(m)%b(2,1)=czero
          lz3_412(m)%c(1,2)=czero
          lz3_412(m)%c(2,1)=czero
          lf3_412(m)%a(1,1)=czero
          lf3_412(m)%a(2,2)=czero
          lf3_412(m)%b(1,2)=czero
          lf3_412(m)%b(2,1)=czero
          lf3_412(m)%c(1,2)=czero
          lf3_412(m)%c(2,1)=czero
        enddo
      else
*zwend                                                                  
  
* quqd -- p=p3,q=p412
      quqd=p3(0)*p412(0)-p3(1)*p412(1)-p3(2)*p412(2)-p3(3)*p412(
     & 3)
* T -- qu=p3,qd=p412,v=0,a=lz3_412(0)%a,b=lz3_412(0)%b,c=lz3_412(0)%c,d=
* lz3_412(0)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      ceps_2=p412(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lz3_412(0)%a(1,1)=zcr(id3)*(auxa+ceps_0)
      lz3_412(0)%a(2,2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(0)%b(1,2)=-zcl(id3)*(p412(2)+ceps_2)
      lz3_412(0)%b(2,1)=zcr(id3)*(p412(2)-ceps_2)
      lz3_412(0)%c(1,2)=zcr(id3)*(p3(2)+ceps_1)
      lz3_412(0)%c(2,1)=zcl(id3)*(-p3(2)+ceps_1)
      lz3_412(0)%d(1,1)=zcl(id3)
      lz3_412(0)%d(2,2)=zcr(id3)
* T -- qu=p3,qd=p412,v=1,a=lz3_412(1)%a,b=lz3_412(1)%b,c=lz3_412(1)%c,d=
* lz3_412(1)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lz3_412(1)%a(1,1)=zcr(id3)*(auxa+ceps_0)
      lz3_412(1)%a(2,2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(1)%b(1,2)=-zcl(id3)*(p412(2)+ceps_2)
      lz3_412(1)%b(2,1)=zcr(id3)*(p412(2)-ceps_2)
      lz3_412(1)%c(1,2)=zcr(id3)*(p3(2)+ceps_1)
      lz3_412(1)%c(2,1)=zcl(id3)*(-p3(2)+ceps_1)
      lz3_412(1)%d(1,1)=zcl(id3)
      lz3_412(1)%d(2,2)=zcr(id3)
* T -- qu=p3,qd=p412,v=2,a=lz3_412(2)%a,b=lz3_412(2)%b,c=lz3_412(2)%c,d=
* lz3_412(2)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lz3_412(2)%a(1,1)=zcr(id3)*(auxa+ceps_0)
      lz3_412(2)%a(2,2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(2)%b(1,2)=-zcl(id3)*p412k0
      lz3_412(2)%b(2,1)=zcr(id3)*p412k0
      lz3_412(2)%c(1,2)=zcr(id3)*p3k0
      lz3_412(2)%c(2,1)=-zcl(id3)*p3k0
* T -- qu=p3,qd=p412,v=3,a=lz3_412(3)%a,b=lz3_412(3)%b,c=lz3_412(3)%c,d=
* lz3_412(3)%d,cr=zcr(id3),cl=zcl(id3),nsum=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      ceps_2=p412k0*cim
      auxa=+p3k0*p412(3)+p412k0*p3(3)
      lz3_412(3)%a(1,1)=zcr(id3)*(auxa+ceps_0)
      lz3_412(3)%a(2,2)=zcl(id3)*(auxa-ceps_0)
      lz3_412(3)%b(1,2)=-zcl(id3)*ceps_2
      lz3_412(3)%b(2,1)=-zcr(id3)*ceps_2
      lz3_412(3)%c(1,2)=zcr(id3)*ceps_1
      lz3_412(3)%c(2,1)=zcl(id3)*ceps_1
  
* T -- qu=p3,qd=p412,v=0,a=lf3_412(0)%a,b=lf3_412(0)%b,c=lf3_412(0)%c,d=
* lf3_412(0)%d,cr=fcr(id3),cl=fcl(id3),nsum=0
      eps_0=-p3(2)*p412(3)+p412(2)*p3(3)
      ceps_0=eps_0*cim
      ceps_1=p3(3)*cim
      ceps_2=p412(3)*cim
      auxa=-quqd+p3k0*p412(0)+p412k0*p3(0)
      lf3_412(0)%a(1,1)=fcr(id3)*(auxa+ceps_0)
      lf3_412(0)%a(2,2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(0)%b(1,2)=-fcl(id3)*(p412(2)+ceps_2)
      lf3_412(0)%b(2,1)=fcr(id3)*(p412(2)-ceps_2)
      lf3_412(0)%c(1,2)=fcr(id3)*(p3(2)+ceps_1)
      lf3_412(0)%c(2,1)=fcl(id3)*(-p3(2)+ceps_1)
      lf3_412(0)%d(1,1)=fcl(id3)
      lf3_412(0)%d(2,2)=fcr(id3)
* T -- qu=p3,qd=p412,v=1,a=lf3_412(1)%a,b=lf3_412(1)%b,c=lf3_412(1)%c,d=
* lf3_412(1)%d,cr=fcr(id3),cl=fcl(id3),nsum=0
      auxa=-quqd+p3k0*p412(1)+p412k0*p3(1)
      lf3_412(1)%a(1,1)=fcr(id3)*(auxa+ceps_0)
      lf3_412(1)%a(2,2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(1)%b(1,2)=-fcl(id3)*(p412(2)+ceps_2)
      lf3_412(1)%b(2,1)=fcr(id3)*(p412(2)-ceps_2)
      lf3_412(1)%c(1,2)=fcr(id3)*(p3(2)+ceps_1)
      lf3_412(1)%c(2,1)=fcl(id3)*(-p3(2)+ceps_1)
      lf3_412(1)%d(1,1)=fcl(id3)
      lf3_412(1)%d(2,2)=fcr(id3)
* T -- qu=p3,qd=p412,v=2,a=lf3_412(2)%a,b=lf3_412(2)%b,c=lf3_412(2)%c,d=
* lf3_412(2)%d,cr=fcr(id3),cl=fcl(id3),nsum=0
      eps_0=-p3k0*p412(3)+p412k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p412(2)+p412k0*p3(2)
      lf3_412(2)%a(1,1)=fcr(id3)*(auxa+ceps_0)
      lf3_412(2)%a(2,2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(2)%b(1,2)=-fcl(id3)*p412k0
      lf3_412(2)%b(2,1)=fcr(id3)*p412k0
      lf3_412(2)%c(1,2)=fcr(id3)*p3k0
      lf3_412(2)%c(2,1)=-fcl(id3)*p3k0
* T -- qu=p3,qd=p412,v=3,a=lf3_412(3)%a,b=lf3_412(3)%b,c=lf3_412(3)%c,d=
* lf3_412(3)%d,cr=fcr(id3),cl=fcl(id3),nsum=0
      eps_0=p3k0*p412(2)-p412k0*p3(2)
      ceps_0=eps_0*cim
      ceps_1=p3k0*cim
      ceps_2=p412k0*cim
      auxa=+p3k0*p412(3)+p412k0*p3(3)
      lf3_412(3)%a(1,1)=fcr(id3)*(auxa+ceps_0)
      lf3_412(3)%a(2,2)=fcl(id3)*(auxa-ceps_0)
      lf3_412(3)%b(1,2)=-fcl(id3)*ceps_2
      lf3_412(3)%b(2,1)=-fcr(id3)*ceps_2
      lf3_412(3)%c(1,2)=fcr(id3)*ceps_1
      lf3_412(3)%c(2,1)=fcl(id3)*ceps_1
*zw                                                                     
      endif   !  if(iresz(3).eq.1)
*zwend                                                                  
  
  
  
  
* compute all single insertions of the type rzi_ and rfi_ (i=2,4)       
*    to a Zline                                                         
*                                                                       
*            |_Z__(mu)         |__f__(mu)                               
*        i __|             i __|                                        
*                                                                       
  
  
* quqd -- p=p134,q=p2
      quqd=p134(0)*p2(0)-p134(1)*p2(1)-p134(2)*p2(2)-p134(3)*p2(
     & 3)
* T -- qu=p134,qd=p2,v=0,a=rz2_134(0)%a,b=rz2_134(0)%b,c=rz2_134(0)%c,d=
* rz2_134(0)%d,cr=zcr(id2),cl=zcl(id2),nsum=0
      eps_0=-p134(2)*p2(3)+p2(2)*p134(3)
      ceps_0=eps_0*cim
      ceps_1=p134(3)*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p134k0*p2(0)+p2k0*p134(0)
      rz2_134(0)%a(1,1)=zcr(id2)*(auxa+ceps_0)
      rz2_134(0)%a(2,2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(0)%b(1,2)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_134(0)%b(2,1)=zcr(id2)*(p2(2)-ceps_2)
      rz2_134(0)%c(1,2)=zcr(id2)*(p134(2)+ceps_1)
      rz2_134(0)%c(2,1)=zcl(id2)*(-p134(2)+ceps_1)
      rz2_134(0)%d(1,1)=zcl(id2)
      rz2_134(0)%d(2,2)=zcr(id2)
* T -- qu=p134,qd=p2,v=1,a=rz2_134(1)%a,b=rz2_134(1)%b,c=rz2_134(1)%c,d=
* rz2_134(1)%d,cr=zcr(id2),cl=zcl(id2),nsum=0
      auxa=-quqd+p134k0*p2(1)+p2k0*p134(1)
      rz2_134(1)%a(1,1)=zcr(id2)*(auxa+ceps_0)
      rz2_134(1)%a(2,2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(1)%b(1,2)=-zcl(id2)*(p2(2)+ceps_2)
      rz2_134(1)%b(2,1)=zcr(id2)*(p2(2)-ceps_2)
      rz2_134(1)%c(1,2)=zcr(id2)*(p134(2)+ceps_1)
      rz2_134(1)%c(2,1)=zcl(id2)*(-p134(2)+ceps_1)
      rz2_134(1)%d(1,1)=zcl(id2)
      rz2_134(1)%d(2,2)=zcr(id2)
* T -- qu=p134,qd=p2,v=2,a=rz2_134(2)%a,b=rz2_134(2)%b,c=rz2_134(2)%c,d=
* rz2_134(2)%d,cr=zcr(id2),cl=zcl(id2),nsum=0
      eps_0=-p134k0*p2(3)+p2k0*p134(3)
      ceps_0=eps_0*cim
      auxa=p134k0*p2(2)+p2k0*p134(2)
      rz2_134(2)%a(1,1)=zcr(id2)*(auxa+ceps_0)
      rz2_134(2)%a(2,2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(2)%b(1,2)=-zcl(id2)*p2k0
      rz2_134(2)%b(2,1)=zcr(id2)*p2k0
      rz2_134(2)%c(1,2)=zcr(id2)*p134k0
      rz2_134(2)%c(2,1)=-zcl(id2)*p134k0
* T -- qu=p134,qd=p2,v=3,a=rz2_134(3)%a,b=rz2_134(3)%b,c=rz2_134(3)%c,d=
* rz2_134(3)%d,cr=zcr(id2),cl=zcl(id2),nsum=0
      eps_0=p134k0*p2(2)-p2k0*p134(2)
      ceps_0=eps_0*cim
      ceps_1=p134k0*cim
      ceps_2=p2k0*cim
      auxa=+p134k0*p2(3)+p2k0*p134(3)
      rz2_134(3)%a(1,1)=zcr(id2)*(auxa+ceps_0)
      rz2_134(3)%a(2,2)=zcl(id2)*(auxa-ceps_0)
      rz2_134(3)%b(1,2)=-zcl(id2)*ceps_2
      rz2_134(3)%b(2,1)=-zcr(id2)*ceps_2
      rz2_134(3)%c(1,2)=zcr(id2)*ceps_1
      rz2_134(3)%c(2,1)=zcl(id2)*ceps_1
  
* T -- qu=p134,qd=p2,v=0,a=rf2_134(0)%a,b=rf2_134(0)%b,c=rf2_134(0)%c,d=
* rf2_134(0)%d,cr=fcr(id2),cl=fcl(id2),nsum=0
      eps_0=-p134(2)*p2(3)+p2(2)*p134(3)
      ceps_0=eps_0*cim
      ceps_1=p134(3)*cim
      ceps_2=p2(3)*cim
      auxa=-quqd+p134k0*p2(0)+p2k0*p134(0)
      rf2_134(0)%a(1,1)=fcr(id2)*(auxa+ceps_0)
      rf2_134(0)%a(2,2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(0)%b(1,2)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_134(0)%b(2,1)=fcr(id2)*(p2(2)-ceps_2)
      rf2_134(0)%c(1,2)=fcr(id2)*(p134(2)+ceps_1)
      rf2_134(0)%c(2,1)=fcl(id2)*(-p134(2)+ceps_1)
      rf2_134(0)%d(1,1)=fcl(id2)
      rf2_134(0)%d(2,2)=fcr(id2)
* T -- qu=p134,qd=p2,v=1,a=rf2_134(1)%a,b=rf2_134(1)%b,c=rf2_134(1)%c,d=
* rf2_134(1)%d,cr=fcr(id2),cl=fcl(id2),nsum=0
      auxa=-quqd+p134k0*p2(1)+p2k0*p134(1)
      rf2_134(1)%a(1,1)=fcr(id2)*(auxa+ceps_0)
      rf2_134(1)%a(2,2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(1)%b(1,2)=-fcl(id2)*(p2(2)+ceps_2)
      rf2_134(1)%b(2,1)=fcr(id2)*(p2(2)-ceps_2)
      rf2_134(1)%c(1,2)=fcr(id2)*(p134(2)+ceps_1)
      rf2_134(1)%c(2,1)=fcl(id2)*(-p134(2)+ceps_1)
      rf2_134(1)%d(1,1)=fcl(id2)
      rf2_134(1)%d(2,2)=fcr(id2)
* T -- qu=p134,qd=p2,v=2,a=rf2_134(2)%a,b=rf2_134(2)%b,c=rf2_134(2)%c,d=
* rf2_134(2)%d,cr=fcr(id2),cl=fcl(id2),nsum=0
      eps_0=-p134k0*p2(3)+p2k0*p134(3)
      ceps_0=eps_0*cim
      auxa=p134k0*p2(2)+p2k0*p134(2)
      rf2_134(2)%a(1,1)=fcr(id2)*(auxa+ceps_0)
      rf2_134(2)%a(2,2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(2)%b(1,2)=-fcl(id2)*p2k0
      rf2_134(2)%b(2,1)=fcr(id2)*p2k0
      rf2_134(2)%c(1,2)=fcr(id2)*p134k0
      rf2_134(2)%c(2,1)=-fcl(id2)*p134k0
* T -- qu=p134,qd=p2,v=3,a=rf2_134(3)%a,b=rf2_134(3)%b,c=rf2_134(3)%c,d=
* rf2_134(3)%d,cr=fcr(id2),cl=fcl(id2),nsum=0
      eps_0=p134k0*p2(2)-p2k0*p134(2)
      ceps_0=eps_0*cim
      ceps_1=p134k0*cim
      ceps_2=p2k0*cim
      auxa=+p134k0*p2(3)+p2k0*p134(3)
      rf2_134(3)%a(1,1)=fcr(id2)*(auxa+ceps_0)
      rf2_134(3)%a(2,2)=fcl(id2)*(auxa-ceps_0)
      rf2_134(3)%b(1,2)=-fcl(id2)*ceps_2
      rf2_134(3)%b(2,1)=-fcr(id2)*ceps_2
      rf2_134(3)%c(1,2)=fcr(id2)*ceps_1
      rf2_134(3)%c(2,1)=fcl(id2)*ceps_1
  
* quqd -- p=p312,q=p4
      quqd=p312(0)*p4(0)-p312(1)*p4(1)-p312(2)*p4(2)-p312(3)*p4(
     & 3)
* T -- qu=p312,qd=p4,v=0,a=rz4_312(0)%a,b=rz4_312(0)%b,c=rz4_312(0)%c,d=
* rz4_312(0)%d,cr=zcr(id4),cl=zcl(id4),nsum=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_1=p312(3)*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rz4_312(0)%a(1,1)=zcr(id4)*(auxa+ceps_0)
      rz4_312(0)%a(2,2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(0)%b(1,2)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_312(0)%b(2,1)=zcr(id4)*(p4(2)-ceps_2)
      rz4_312(0)%c(1,2)=zcr(id4)*(p312(2)+ceps_1)
      rz4_312(0)%c(2,1)=zcl(id4)*(-p312(2)+ceps_1)
      rz4_312(0)%d(1,1)=zcl(id4)
      rz4_312(0)%d(2,2)=zcr(id4)
* T -- qu=p312,qd=p4,v=1,a=rz4_312(1)%a,b=rz4_312(1)%b,c=rz4_312(1)%c,d=
* rz4_312(1)%d,cr=zcr(id4),cl=zcl(id4),nsum=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rz4_312(1)%a(1,1)=zcr(id4)*(auxa+ceps_0)
      rz4_312(1)%a(2,2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(1)%b(1,2)=-zcl(id4)*(p4(2)+ceps_2)
      rz4_312(1)%b(2,1)=zcr(id4)*(p4(2)-ceps_2)
      rz4_312(1)%c(1,2)=zcr(id4)*(p312(2)+ceps_1)
      rz4_312(1)%c(2,1)=zcl(id4)*(-p312(2)+ceps_1)
      rz4_312(1)%d(1,1)=zcl(id4)
      rz4_312(1)%d(2,2)=zcr(id4)
* T -- qu=p312,qd=p4,v=2,a=rz4_312(2)%a,b=rz4_312(2)%b,c=rz4_312(2)%c,d=
* rz4_312(2)%d,cr=zcr(id4),cl=zcl(id4),nsum=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rz4_312(2)%a(1,1)=zcr(id4)*(auxa+ceps_0)
      rz4_312(2)%a(2,2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(2)%b(1,2)=-zcl(id4)*p4k0
      rz4_312(2)%b(2,1)=zcr(id4)*p4k0
      rz4_312(2)%c(1,2)=zcr(id4)*p312k0
      rz4_312(2)%c(2,1)=-zcl(id4)*p312k0
* T -- qu=p312,qd=p4,v=3,a=rz4_312(3)%a,b=rz4_312(3)%b,c=rz4_312(3)%c,d=
* rz4_312(3)%d,cr=zcr(id4),cl=zcl(id4),nsum=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_1=p312k0*cim
      ceps_2=p4k0*cim
      auxa=+p312k0*p4(3)+p4k0*p312(3)
      rz4_312(3)%a(1,1)=zcr(id4)*(auxa+ceps_0)
      rz4_312(3)%a(2,2)=zcl(id4)*(auxa-ceps_0)
      rz4_312(3)%b(1,2)=-zcl(id4)*ceps_2
      rz4_312(3)%b(2,1)=-zcr(id4)*ceps_2
      rz4_312(3)%c(1,2)=zcr(id4)*ceps_1
      rz4_312(3)%c(2,1)=zcl(id4)*ceps_1
  
* T -- qu=p312,qd=p4,v=0,a=rf4_312(0)%a,b=rf4_312(0)%b,c=rf4_312(0)%c,d=
* rf4_312(0)%d,cr=fcr(id4),cl=fcl(id4),nsum=0
      eps_0=-p312(2)*p4(3)+p4(2)*p312(3)
      ceps_0=eps_0*cim
      ceps_1=p312(3)*cim
      ceps_2=p4(3)*cim
      auxa=-quqd+p312k0*p4(0)+p4k0*p312(0)
      rf4_312(0)%a(1,1)=fcr(id4)*(auxa+ceps_0)
      rf4_312(0)%a(2,2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(0)%b(1,2)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_312(0)%b(2,1)=fcr(id4)*(p4(2)-ceps_2)
      rf4_312(0)%c(1,2)=fcr(id4)*(p312(2)+ceps_1)
      rf4_312(0)%c(2,1)=fcl(id4)*(-p312(2)+ceps_1)
      rf4_312(0)%d(1,1)=fcl(id4)
      rf4_312(0)%d(2,2)=fcr(id4)
* T -- qu=p312,qd=p4,v=1,a=rf4_312(1)%a,b=rf4_312(1)%b,c=rf4_312(1)%c,d=
* rf4_312(1)%d,cr=fcr(id4),cl=fcl(id4),nsum=0
      auxa=-quqd+p312k0*p4(1)+p4k0*p312(1)
      rf4_312(1)%a(1,1)=fcr(id4)*(auxa+ceps_0)
      rf4_312(1)%a(2,2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(1)%b(1,2)=-fcl(id4)*(p4(2)+ceps_2)
      rf4_312(1)%b(2,1)=fcr(id4)*(p4(2)-ceps_2)
      rf4_312(1)%c(1,2)=fcr(id4)*(p312(2)+ceps_1)
      rf4_312(1)%c(2,1)=fcl(id4)*(-p312(2)+ceps_1)
      rf4_312(1)%d(1,1)=fcl(id4)
      rf4_312(1)%d(2,2)=fcr(id4)
* T -- qu=p312,qd=p4,v=2,a=rf4_312(2)%a,b=rf4_312(2)%b,c=rf4_312(2)%c,d=
* rf4_312(2)%d,cr=fcr(id4),cl=fcl(id4),nsum=0
      eps_0=-p312k0*p4(3)+p4k0*p312(3)
      ceps_0=eps_0*cim
      auxa=p312k0*p4(2)+p4k0*p312(2)
      rf4_312(2)%a(1,1)=fcr(id4)*(auxa+ceps_0)
      rf4_312(2)%a(2,2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(2)%b(1,2)=-fcl(id4)*p4k0
      rf4_312(2)%b(2,1)=fcr(id4)*p4k0
      rf4_312(2)%c(1,2)=fcr(id4)*p312k0
      rf4_312(2)%c(2,1)=-fcl(id4)*p312k0
* T -- qu=p312,qd=p4,v=3,a=rf4_312(3)%a,b=rf4_312(3)%b,c=rf4_312(3)%c,d=
* rf4_312(3)%d,cr=fcr(id4),cl=fcl(id4),nsum=0
      eps_0=p312k0*p4(2)-p4k0*p312(2)
      ceps_0=eps_0*cim
      ceps_1=p312k0*cim
      ceps_2=p4k0*cim
      auxa=+p312k0*p4(3)+p4k0*p312(3)
      rf4_312(3)%a(1,1)=fcr(id4)*(auxa+ceps_0)
      rf4_312(3)%a(2,2)=fcl(id4)*(auxa-ceps_0)
      rf4_312(3)%b(1,2)=-fcl(id4)*ceps_2
      rf4_312(3)%b(2,1)=-fcr(id4)*ceps_2
      rf4_312(3)%c(1,2)=fcr(id4)*ceps_1
      rf4_312(3)%c(2,1)=fcl(id4)*ceps_1
  
  
* compute all single insertions of the type lhi_ (i=1,3)                
*    to a Zline                                                         
*                                                                       
*        i __                                                           
*            |_h__                                                      
*            |                                                          
* rmh < 0 in our convention means no higgs coupling                     
      if (id1.eq.5.and.rmh.ge.0.d0) then
*zw                                                                     
      if(iresz(1).eq.1) then
        lh1_234%a(1,2)=czero
        lh1_234%a(2,1)=czero
        lh1_234%b(1,1)=czero
        lh1_234%b(2,2)=czero
        lh1_234%c(1,1)=czero
        lh1_234%c(2,2)=czero
      else
*zwend                                                                  
* TH -- qu=p1,qd=p234,a=lh1_234%a,b=lh1_234%b,c=lh1_234%c,coupl=1.d0
      auxa=-p1k0*p234(2)+p234k0*p1(2)
      cauxa=auxa-cim*(p234(3)*p1k0-p1(3)*p234k0)
      lh1_234%a(1,2)=cauxa
      lh1_234%a(2,1)=-conjg(cauxa)
      lh1_234%b(1,1)=p234k0
      lh1_234%b(2,2)=lh1_234%b(1,1)
      lh1_234%c(1,1)=p1k0
      lh1_234%c(2,2)=lh1_234%c(1,1)
*zw                                                                     
      endif   !  if(iresz(1).eq.1)
*zwend                                                                  
  
      endif
* rmh < 0 in our convention means no higgs coupling                     
      if (id3.eq.5.and.rmh.ge.0.d0) then
*zw                                                                     
      if(iresz(3).eq.1) then
        lh3_412%a(1,2)=czero
        lh3_412%a(2,1)=czero
        lh3_412%b(1,1)=czero
        lh3_412%b(2,2)=czero
        lh3_412%c(1,1)=czero
        lh3_412%c(2,2)=czero
      else
*zwend                                                                  
* TH -- qu=p3,qd=p412,a=lh3_412%a,b=lh3_412%b,c=lh3_412%c,coupl=1.d0
      auxa=-p3k0*p412(2)+p412k0*p3(2)
      cauxa=auxa-cim*(p412(3)*p3k0-p3(3)*p412k0)
      lh3_412%a(1,2)=cauxa
      lh3_412%a(2,1)=-conjg(cauxa)
      lh3_412%b(1,1)=p412k0
      lh3_412%b(2,2)=lh3_412%b(1,1)
      lh3_412%c(1,1)=p3k0
      lh3_412%c(2,2)=lh3_412%c(1,1)
*zw                                                                     
      endif   !  if(iresz(3).eq.1)
*zwend                                                                  
  
      endif
* compute all single insertions of the type rhi_ (i=2,4)                
*    to a Zline                                                         
*                                                                       
*            |_h__                                                      
*        i __|                                                          
*                                                                       
* rmh < 0 in our convention means no higgs coupling                     
      if (id2.eq.-5.and.rmh.ge.0.d0) then
* TH -- qu=p134,qd=p2,a=rh2_134%a,b=rh2_134%b,c=rh2_134%c,coupl=1.d0
      auxa=-p134k0*p2(2)+p2k0*p134(2)
      cauxa=auxa-cim*(p2(3)*p134k0-p134(3)*p2k0)
      rh2_134%a(1,2)=cauxa
      rh2_134%a(2,1)=-conjg(cauxa)
      rh2_134%b(1,1)=p2k0
      rh2_134%b(2,2)=rh2_134%b(1,1)
      rh2_134%c(1,1)=p134k0
      rh2_134%c(2,2)=rh2_134%c(1,1)
      endif
* rmh < 0 in our convention means no higgs coupling                     
      if (id4.eq.-5.and.rmh.ge.0.d0) then
* TH -- qu=p312,qd=p4,a=rh4_312%a,b=rh4_312%b,c=rh4_312%c,coupl=1.d0
      auxa=-p312k0*p4(2)+p4k0*p312(2)
      cauxa=auxa-cim*(p4(3)*p312k0-p312(3)*p4k0)
      rh4_312%a(1,2)=cauxa
      rh4_312%a(2,1)=-conjg(cauxa)
      rh4_312%b(1,1)=p4k0
      rh4_312%b(2,2)=rh4_312%b(1,1)
      rh4_312%c(1,1)=p312k0
      rh4_312%c(2,2)=rh4_312%c(1,1)
      endif
  
  
  
* compute all subdiagrams with a Z or a gamma "decaying" to 4 fermions  
*  which correspond to two w's                                          
*                  czijkl(mu)  and cfijkl(mu) (ijkl particle numbers)   
*                                                                       
*               _____ i                       _____ i                   
*              |  |__ j                      |  |__ j                   
*    (mu) --Z--|  |__ k            (mu) --f--|  |__ k                   
*              |__|__ l                      |__|__ l                   
*                                                                       
  
  
      do m=0,3
        cz5678(m)=czero
      enddo
  
*                                                                       
*    First all diagrams of the type                                     
*                                                                       
*                 |_W__/                           |_W__/               
*       (mu) _Z___|    \       and        (mu) _f__|    \               
*                 |                                |                    
  
      do m=0,3
* TLTR0_W -- aa=cz5678(m),a1=l5_78%a,c1=l5_78%c,a2=rz6_578(m)%a,b2=rz6_5
* 78(m)%b,prq=p578q,bef=cz5678(m)+,aft=
      cz5678(m)=cz5678(m)+(l5_78%c(2)*p578q*rz6_578(m)%b(1)+l5_7
     & 8%a(2)*rz6_578(m)%a(2))
      end do
  
*    Then those of the type                                             
*                                                                       
*       (mu) _Z___|                      (mu) _f__|                     
*                 |_W__/       and                |_W__/                
*                 |    \                          |    \                
*                                                                       
  
      do m=0,3
* TLTR0_W -- aa=cz5678(m),a1=lz5_678(m)%a,c1=lz5_678(m)%c,a2=r6_78%a,b2=
* r6_78%b,prq=p678q,bef=cz5678(m)+,aft=
      cz5678(m)=cz5678(m)+(lz5_678(m)%c(2)*p678q*r6_78%b(1)+lz5_
     & 678(m)%a(2)*r6_78%a(2))
      end do
*                                                                       
*    First all diagrams of the type                                     
*                                                                       
*                 |_W__/                           |_W__/               
*       (mu) _Z___|    \       and        (mu) _f__|    \               
*                 |                                |                    
  
      do m=0,3
* TLTR0_W -- aa=cz5678(m),a1=l7_56%a,c1=l7_56%c,a2=rz8_756(m)%a,b2=rz8_7
* 56(m)%b,prq=p756q,bef=cz5678(m)+,aft=
      cz5678(m)=cz5678(m)+(l7_56%c(2)*p756q*rz8_756(m)%b(1)+l7_5
     & 6%a(2)*rz8_756(m)%a(2))
      end do
  
*    Then those of the type                                             
*                                                                       
*       (mu) _Z___|                      (mu) _f__|                     
*                 |_W__/       and                |_W__/                
*                 |    \                          |    \                
*                                                                       
  
      do m=0,3
* TLTR0_W -- aa=cz5678(m),a1=lz7_856(m)%a,c1=lz7_856(m)%c,a2=r8_56%a,b2=
* r8_56%b,prq=p856q,bef=cz5678(m)+,aft=
      cz5678(m)=cz5678(m)+(lz7_856(m)%c(2)*p856q*r8_56%b(1)+lz7_
     & 856(m)%a(2)*r8_56%a(2))
      end do
  
*    Then those with triple vertex                                      
*                                                                       
*                    /_                         /_                      
*                  W/                         W/                        
*         (mu)_Z___/          and     (mu)_f__/                         
*                  \                          \                         
*                  W\_                        W\_                       
*                    \                          \                       
*                                                                       
      do mu=0,3
        p1234(mu)=-p56(mu)-p78(mu)
      enddo
* p.q -- p.q=p1234q,p=p1234,q=p1234,bef=,aft=
      p1234q=(p1234(0)*p1234(0)-p1234(1)*p1234(1)-p1234(2)*p1234
     & (2)-p1234(3)*p1234(3))
  
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
        cz5678(mu)=rcotw*ctrip5678(mu)+cz5678(mu)
      enddo
  
  
      do m=0,3
        cf5678(m)=czero
      enddo
  
*                                                                       
*    First all diagrams of the type                                     
*                                                                       
*                 |_W__/                           |_W__/               
*       (mu) _Z___|    \       and        (mu) _f__|    \               
*                 |                                |                    
  
      do m=0,3
* TLTR0_W -- aa=cf5678(m),a1=l5_78%a,c1=l5_78%c,a2=rf6_578(m)%a,b2=rf6_5
* 78(m)%b,prq=p578q,bef=cf5678(m)+,aft=
      cf5678(m)=cf5678(m)+(l5_78%c(2)*p578q*rf6_578(m)%b(1)+l5_7
     & 8%a(2)*rf6_578(m)%a(2))
      end do
  
*    Then those of the type                                             
*                                                                       
*       (mu) _Z___|                      (mu) _f__|                     
*                 |_W__/       and                |_W__/                
*                 |    \                          |    \                
*                                                                       
  
      do m=0,3
* TLTR0_W -- aa=cf5678(m),a1=lf5_678(m)%a,c1=lf5_678(m)%c,a2=r6_78%a,b2=
* r6_78%b,prq=p678q,bef=cf5678(m)+,aft=
      cf5678(m)=cf5678(m)+(lf5_678(m)%c(2)*p678q*r6_78%b(1)+lf5_
     & 678(m)%a(2)*r6_78%a(2))
      end do
*                                                                       
*    First all diagrams of the type                                     
*                                                                       
*                 |_W__/                           |_W__/               
*       (mu) _Z___|    \       and        (mu) _f__|    \               
*                 |                                |                    
  
      do m=0,3
* TLTR0_W -- aa=cf5678(m),a1=l7_56%a,c1=l7_56%c,a2=rf8_756(m)%a,b2=rf8_7
* 56(m)%b,prq=p756q,bef=cf5678(m)+,aft=
      cf5678(m)=cf5678(m)+(l7_56%c(2)*p756q*rf8_756(m)%b(1)+l7_5
     & 6%a(2)*rf8_756(m)%a(2))
      end do
  
*    Then those of the type                                             
*                                                                       
*       (mu) _Z___|                      (mu) _f__|                     
*                 |_W__/       and                |_W__/                
*                 |    \                          |    \                
*                                                                       
  
      do m=0,3
* TLTR0_W -- aa=cf5678(m),a1=lf7_856(m)%a,c1=lf7_856(m)%c,a2=r8_56%a,b2=
* r8_56%b,prq=p856q,bef=cf5678(m)+,aft=
      cf5678(m)=cf5678(m)+(lf7_856(m)%c(2)*p856q*r8_56%b(1)+lf7_
     & 856(m)%a(2)*r8_56%a(2))
      end do
  
*    Then those with triple vertex                                      
*                                                                       
*                    /_                         /_                      
*                  W/                         W/                        
*         (mu)_Z___/          and     (mu)_f__/                         
*                  \                          \                         
*                  W\_                        W\_                       
*                    \                          \                       
*                                                                       
  
  
      do mu=0,3
        cf5678(mu)=ctrip5678(mu)+cf5678(mu)
      enddo
  
  
* compute  subdiagram with a higgs decaying to 4 fermions (two w's)     
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
* coupling g/|e|=1/sw                                                   
  
* rmh < 0 in our convention means no higgs coupling                     
      if(rmh.ge.0.d0) then
* p.q -- p.q=ch5678,p=cw56%e,q=cw78%e,bef=,aft=
      ch5678=(cw56%e(0)*cw78%e(0)-cw56%e(1)*cw78%e(1)-cw56%e(2)*
     & cw78%e(2)-cw56%e(3)*cw78%e(3))
      chh5678=ch5678*rhhww
      ch5678=ch5678*rhww
      endif
  
* compute all subdiagrams with a W "decaying" to 4 fermions             
*  which correspond to one Z and one W                                  
*                  cwijkl(mu)                                           
*                                                                       
*                      _____ i                                          
*                     |  |__ j                                          
*           (mu) --W--|  |__ k                                          
*                     |__|__ l                                          
*                                                                       
  
  
  
*                                                                       
*    First all diagrams of the type                                     
*                                                                       
*                 |_W__/                                                
*       (mu) _W___|    \                                                
*                 |                                                     
  
*                 first    without sum                                  
  
  
      if (iup(id1).eq.0) then
* 56 in the upper part                                                  
      do mu=0,3
* TWTW -- aa=clinetww_mu(mu)%a,bb=clinetww_mu(mu)%b,cc=clinetww_mu(mu)%c
* ,dd=clinetww_mu(mu)%d,a1=l1_56%a,b1=l1_56%b,c1=l1_56%c,d1=l1_56%d,a2=r
* w2_156(mu)%a,b2=rw2_156(mu)%b,c2=rw2_156(mu)%c,d2=rw2_156(mu)%d,prq=p1
* 56q,nsum=0
      clinetww_mu(mu)%d(1,1)=l1_56%d(1,1)*p156q*rw2_156(mu)%d(1,
     & 1)+l1_56%b(1,2)*rw2_156(mu)%c(2,1)
      clinetww_mu(mu)%b(1,2)=l1_56%d(1,1)*p156q*rw2_156(mu)%b(1,
     & 2)+l1_56%b(1,2)*rw2_156(mu)%a(2,2)
      clinetww_mu(mu)%c(2,1)=l1_56%c(2,1)*p156q*rw2_156(mu)%d(1,
     & 1)+l1_56%a(2,2)*rw2_156(mu)%c(2,1)
      clinetww_mu(mu)%a(2,2)=l1_56%c(2,1)*p156q*rw2_156(mu)%b(1,
     & 2)+l1_56%a(2,2)*rw2_156(mu)%a(2,2)
      end do
      else
* 56 in the lower part                                                  
      do mu=0,3
* TWTW -- aa=clinetww_mu(mu)%a,bb=clinetww_mu(mu)%b,cc=clinetww_mu(mu)%c
* ,dd=clinetww_mu(mu)%d,a1=lw1_256(mu)%a,b1=lw1_256(mu)%b,c1=lw1_256(mu)
* %c,d1=lw1_256(mu)%d,a2=r2_56%a,b2=r2_56%b,c2=r2_56%c,d2=r2_56%d,prq=p2
* 56q,nsum=0
      clinetww_mu(mu)%d(1,1)=lw1_256(mu)%d(1,1)*p256q*r2_56%d(1,
     & 1)+lw1_256(mu)%b(1,2)*r2_56%c(2,1)
      clinetww_mu(mu)%b(1,2)=lw1_256(mu)%d(1,1)*p256q*r2_56%b(1,
     & 2)+lw1_256(mu)%b(1,2)*r2_56%a(2,2)
      clinetww_mu(mu)%c(2,1)=lw1_256(mu)%c(2,1)*p256q*r2_56%d(1,
     & 1)+lw1_256(mu)%a(2,2)*r2_56%c(2,1)
      clinetww_mu(mu)%a(2,2)=lw1_256(mu)%c(2,1)*p256q*r2_56%b(1,
     & 2)+lw1_256(mu)%a(2,2)*r2_56%a(2,2)
      end do
      endif
      rmassl1=rmass(id1)
      rmassr2=-rmass(id2)
      do mu=0,3
* mline -- res=cw1256(&1,&2,mu),abcd=clinetww_mu(mu)%,m1=rmassl1,m2=rmas
* sr2,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cw1256(iut,jut,mu)=clinetww_mu(mu)%a(iut,jut)+rmassl1*clin
     & etww_mu(mu)%b(iut,jut)+rmassr2*clinetww_mu(mu)%c(iut,jut)
     & +rmassl1*rmassr2*clinetww_mu(mu)%d(iut,jut)
      enddo
      enddo
      end do
      if (iup(id3).eq.0) then
* 56 in the upper part                                                  
      do mu=0,3
* TWTW -- aa=clinetww_mu(mu)%a,bb=clinetww_mu(mu)%b,cc=clinetww_mu(mu)%c
* ,dd=clinetww_mu(mu)%d,a1=l3_56%a,b1=l3_56%b,c1=l3_56%c,d1=l3_56%d,a2=r
* w4_356(mu)%a,b2=rw4_356(mu)%b,c2=rw4_356(mu)%c,d2=rw4_356(mu)%d,prq=p3
* 56q,nsum=0
      clinetww_mu(mu)%d(1,1)=l3_56%d(1,1)*p356q*rw4_356(mu)%d(1,
     & 1)+l3_56%b(1,2)*rw4_356(mu)%c(2,1)
      clinetww_mu(mu)%b(1,2)=l3_56%d(1,1)*p356q*rw4_356(mu)%b(1,
     & 2)+l3_56%b(1,2)*rw4_356(mu)%a(2,2)
      clinetww_mu(mu)%c(2,1)=l3_56%c(2,1)*p356q*rw4_356(mu)%d(1,
     & 1)+l3_56%a(2,2)*rw4_356(mu)%c(2,1)
      clinetww_mu(mu)%a(2,2)=l3_56%c(2,1)*p356q*rw4_356(mu)%b(1,
     & 2)+l3_56%a(2,2)*rw4_356(mu)%a(2,2)
      end do
      else
* 56 in the lower part                                                  
      do mu=0,3
* TWTW -- aa=clinetww_mu(mu)%a,bb=clinetww_mu(mu)%b,cc=clinetww_mu(mu)%c
* ,dd=clinetww_mu(mu)%d,a1=lw3_456(mu)%a,b1=lw3_456(mu)%b,c1=lw3_456(mu)
* %c,d1=lw3_456(mu)%d,a2=r4_56%a,b2=r4_56%b,c2=r4_56%c,d2=r4_56%d,prq=p4
* 56q,nsum=0
      clinetww_mu(mu)%d(1,1)=lw3_456(mu)%d(1,1)*p456q*r4_56%d(1,
     & 1)+lw3_456(mu)%b(1,2)*r4_56%c(2,1)
      clinetww_mu(mu)%b(1,2)=lw3_456(mu)%d(1,1)*p456q*r4_56%b(1,
     & 2)+lw3_456(mu)%b(1,2)*r4_56%a(2,2)
      clinetww_mu(mu)%c(2,1)=lw3_456(mu)%c(2,1)*p456q*r4_56%d(1,
     & 1)+lw3_456(mu)%a(2,2)*r4_56%c(2,1)
      clinetww_mu(mu)%a(2,2)=lw3_456(mu)%c(2,1)*p456q*r4_56%b(1,
     & 2)+lw3_456(mu)%a(2,2)*r4_56%a(2,2)
      end do
      endif
      rmassl3=rmass(id3)
      rmassr4=-rmass(id4)
      do mu=0,3
* mline -- res=cw3456(&1,&2,mu),abcd=clinetww_mu(mu)%,m1=rmassl3,m2=rmas
* sr4,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cw3456(iut,jut,mu)=clinetww_mu(mu)%a(iut,jut)+rmassl3*clin
     & etww_mu(mu)%b(iut,jut)+rmassr4*clinetww_mu(mu)%c(iut,jut)
     & +rmassl3*rmassr4*clinetww_mu(mu)%d(iut,jut)
      enddo
      enddo
      end do
      if (iup(id1).eq.1) then
* 78 in the upper part                                                  
      do mu=0,3
* TWTW -- aa=clinetww_mu(mu)%a,bb=clinetww_mu(mu)%b,cc=clinetww_mu(mu)%c
* ,dd=clinetww_mu(mu)%d,a1=l1_78%a,b1=l1_78%b,c1=l1_78%c,d1=l1_78%d,a2=r
* w2_178(mu)%a,b2=rw2_178(mu)%b,c2=rw2_178(mu)%c,d2=rw2_178(mu)%d,prq=p1
* 78q,nsum=0
      clinetww_mu(mu)%d(1,1)=l1_78%d(1,1)*p178q*rw2_178(mu)%d(1,
     & 1)+l1_78%b(1,2)*rw2_178(mu)%c(2,1)
      clinetww_mu(mu)%b(1,2)=l1_78%d(1,1)*p178q*rw2_178(mu)%b(1,
     & 2)+l1_78%b(1,2)*rw2_178(mu)%a(2,2)
      clinetww_mu(mu)%c(2,1)=l1_78%c(2,1)*p178q*rw2_178(mu)%d(1,
     & 1)+l1_78%a(2,2)*rw2_178(mu)%c(2,1)
      clinetww_mu(mu)%a(2,2)=l1_78%c(2,1)*p178q*rw2_178(mu)%b(1,
     & 2)+l1_78%a(2,2)*rw2_178(mu)%a(2,2)
      end do
      else
* 78 in the lower part                                                  
      do mu=0,3
* TWTW -- aa=clinetww_mu(mu)%a,bb=clinetww_mu(mu)%b,cc=clinetww_mu(mu)%c
* ,dd=clinetww_mu(mu)%d,a1=lw1_278(mu)%a,b1=lw1_278(mu)%b,c1=lw1_278(mu)
* %c,d1=lw1_278(mu)%d,a2=r2_78%a,b2=r2_78%b,c2=r2_78%c,d2=r2_78%d,prq=p2
* 78q,nsum=0
      clinetww_mu(mu)%d(1,1)=lw1_278(mu)%d(1,1)*p278q*r2_78%d(1,
     & 1)+lw1_278(mu)%b(1,2)*r2_78%c(2,1)
      clinetww_mu(mu)%b(1,2)=lw1_278(mu)%d(1,1)*p278q*r2_78%b(1,
     & 2)+lw1_278(mu)%b(1,2)*r2_78%a(2,2)
      clinetww_mu(mu)%c(2,1)=lw1_278(mu)%c(2,1)*p278q*r2_78%d(1,
     & 1)+lw1_278(mu)%a(2,2)*r2_78%c(2,1)
      clinetww_mu(mu)%a(2,2)=lw1_278(mu)%c(2,1)*p278q*r2_78%b(1,
     & 2)+lw1_278(mu)%a(2,2)*r2_78%a(2,2)
      end do
      endif
      do mu=0,3
* mline -- res=cw1278(&1,&2,mu),abcd=clinetww_mu(mu)%,m1=rmassl1,m2=rmas
* sr2,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cw1278(iut,jut,mu)=clinetww_mu(mu)%a(iut,jut)+rmassl1*clin
     & etww_mu(mu)%b(iut,jut)+rmassr2*clinetww_mu(mu)%c(iut,jut)
     & +rmassl1*rmassr2*clinetww_mu(mu)%d(iut,jut)
      enddo
      enddo
      end do
      if (iup(id3).eq.1) then
* 78 in the upper part                                                  
      do mu=0,3
* TWTW -- aa=clinetww_mu(mu)%a,bb=clinetww_mu(mu)%b,cc=clinetww_mu(mu)%c
* ,dd=clinetww_mu(mu)%d,a1=l3_78%a,b1=l3_78%b,c1=l3_78%c,d1=l3_78%d,a2=r
* w4_378(mu)%a,b2=rw4_378(mu)%b,c2=rw4_378(mu)%c,d2=rw4_378(mu)%d,prq=p3
* 78q,nsum=0
      clinetww_mu(mu)%d(1,1)=l3_78%d(1,1)*p378q*rw4_378(mu)%d(1,
     & 1)+l3_78%b(1,2)*rw4_378(mu)%c(2,1)
      clinetww_mu(mu)%b(1,2)=l3_78%d(1,1)*p378q*rw4_378(mu)%b(1,
     & 2)+l3_78%b(1,2)*rw4_378(mu)%a(2,2)
      clinetww_mu(mu)%c(2,1)=l3_78%c(2,1)*p378q*rw4_378(mu)%d(1,
     & 1)+l3_78%a(2,2)*rw4_378(mu)%c(2,1)
      clinetww_mu(mu)%a(2,2)=l3_78%c(2,1)*p378q*rw4_378(mu)%b(1,
     & 2)+l3_78%a(2,2)*rw4_378(mu)%a(2,2)
      end do
      else
* 78 in the lower part                                                  
      do mu=0,3
* TWTW -- aa=clinetww_mu(mu)%a,bb=clinetww_mu(mu)%b,cc=clinetww_mu(mu)%c
* ,dd=clinetww_mu(mu)%d,a1=lw3_478(mu)%a,b1=lw3_478(mu)%b,c1=lw3_478(mu)
* %c,d1=lw3_478(mu)%d,a2=r4_78%a,b2=r4_78%b,c2=r4_78%c,d2=r4_78%d,prq=p4
* 78q,nsum=0
      clinetww_mu(mu)%d(1,1)=lw3_478(mu)%d(1,1)*p478q*r4_78%d(1,
     & 1)+lw3_478(mu)%b(1,2)*r4_78%c(2,1)
      clinetww_mu(mu)%b(1,2)=lw3_478(mu)%d(1,1)*p478q*r4_78%b(1,
     & 2)+lw3_478(mu)%b(1,2)*r4_78%a(2,2)
      clinetww_mu(mu)%c(2,1)=lw3_478(mu)%c(2,1)*p478q*r4_78%d(1,
     & 1)+lw3_478(mu)%a(2,2)*r4_78%c(2,1)
      clinetww_mu(mu)%a(2,2)=lw3_478(mu)%c(2,1)*p478q*r4_78%b(1,
     & 2)+lw3_478(mu)%a(2,2)*r4_78%a(2,2)
      end do
      endif
      do mu=0,3
* mline -- res=cw3478(&1,&2,mu),abcd=clinetww_mu(mu)%,m1=rmassl3,m2=rmas
* sr4,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cw3478(iut,jut,mu)=clinetww_mu(mu)%a(iut,jut)+rmassl3*clin
     & etww_mu(mu)%b(iut,jut)+rmassr4*clinetww_mu(mu)%c(iut,jut)
     & +rmassl3*rmassr4*clinetww_mu(mu)%d(iut,jut)
      enddo
      enddo
      end do
  
*    Then those of the type                                             
*                                                                       
*       (mu) _W___|                                                     
*                 |_Z,f__/                                              
*                 |      \                                              
*                                                                       
  
*                        all  with sum                                  
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw1256(i1,i2,mu),a1=lw5_612(mu)%a,c1=lw5_612(mu)%c,a2=r6
* _12(i1,i2)%a,b2=r6_12(i1,i2)%b,prq=p612q,bef=cw1256(i1,i2,mu)+,aft=
      cw1256(i1,i2,mu)=cw1256(i1,i2,mu)+(lw5_612(mu)%c(2)*p612q*
     & r6_12(i1,i2)%b(1)+lw5_612(mu)%a(2)*r6_12(i1,i2)%a(2))
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw1278(i1,i2,mu),a1=lw7_812(mu)%a,c1=lw7_812(mu)%c,a2=r8
* _12(i1,i2)%a,b2=r8_12(i1,i2)%b,prq=p812q,bef=cw1278(i1,i2,mu)+,aft=
      cw1278(i1,i2,mu)=cw1278(i1,i2,mu)+(lw7_812(mu)%c(2)*p812q*
     & r8_12(i1,i2)%b(1)+lw7_812(mu)%a(2)*r8_12(i1,i2)%a(2))
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw3456(i1,i2,mu),a1=lw5_634(mu)%a,c1=lw5_634(mu)%c,a2=r6
* _34(i1,i2)%a,b2=r6_34(i1,i2)%b,prq=p634q,bef=cw3456(i1,i2,mu)+,aft=
      cw3456(i1,i2,mu)=cw3456(i1,i2,mu)+(lw5_634(mu)%c(2)*p634q*
     & r6_34(i1,i2)%b(1)+lw5_634(mu)%a(2)*r6_34(i1,i2)%a(2))
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw3478(i1,i2,mu),a1=lw7_834(mu)%a,c1=lw7_834(mu)%c,a2=r8
* _34(i1,i2)%a,b2=r8_34(i1,i2)%b,prq=p834q,bef=cw3478(i1,i2,mu)+,aft=
      cw3478(i1,i2,mu)=cw3478(i1,i2,mu)+(lw7_834(mu)%c(2)*p834q*
     & r8_34(i1,i2)%b(1)+lw7_834(mu)%a(2)*r8_34(i1,i2)%a(2))
      end do
      end do
      end do
  
*     And of the type                                                   
*                                                                       
*                 |_Z,f__/                                              
*       (mu) _W___|      \                                              
*                 |                                                     
  
*                        all  with sum                                  
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw1256(i1,i2,mu),a1=l5_12(i1,i2)%a,c1=l5_12(i1,i2)%c,a2=
* rw6_512(mu)%a,b2=rw6_512(mu)%b,prq=p512q,bef=cw1256(i1,i2,mu)+,aft=
      cw1256(i1,i2,mu)=cw1256(i1,i2,mu)+(l5_12(i1,i2)%c(2)*p512q
     & *rw6_512(mu)%b(1)+l5_12(i1,i2)%a(2)*rw6_512(mu)%a(2))
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw1278(i1,i2,mu),a1=l7_12(i1,i2)%a,c1=l7_12(i1,i2)%c,a2=
* rw8_712(mu)%a,b2=rw8_712(mu)%b,prq=p712q,bef=cw1278(i1,i2,mu)+,aft=
      cw1278(i1,i2,mu)=cw1278(i1,i2,mu)+(l7_12(i1,i2)%c(2)*p712q
     & *rw8_712(mu)%b(1)+l7_12(i1,i2)%a(2)*rw8_712(mu)%a(2))
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw3456(i1,i2,mu),a1=l5_34(i1,i2)%a,c1=l5_34(i1,i2)%c,a2=
* rw6_534(mu)%a,b2=rw6_534(mu)%b,prq=p534q,bef=cw3456(i1,i2,mu)+,aft=
      cw3456(i1,i2,mu)=cw3456(i1,i2,mu)+(l5_34(i1,i2)%c(2)*p534q
     & *rw6_534(mu)%b(1)+l5_34(i1,i2)%a(2)*rw6_534(mu)%a(2))
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TLTR0_W -- aa=cw3478(i1,i2,mu),a1=l7_34(i1,i2)%a,c1=l7_34(i1,i2)%c,a2=
* rw8_734(mu)%a,b2=rw8_734(mu)%b,prq=p734q,bef=cw3478(i1,i2,mu)+,aft=
      cw3478(i1,i2,mu)=cw3478(i1,i2,mu)+(l7_34(i1,i2)%c(2)*p734q
     & *rw8_734(mu)%b(1)+l7_34(i1,i2)%a(2)*rw8_734(mu)%a(2))
      end do
      end do
      end do
  
  
  
*    Then those with triple vertex                                      
*                                                                       
*                    /_                                                 
*                  W/                                                   
*         (mu)_W___/                                                    
*                  \                                                    
*                Z,f\_                                                  
*                    \                                                  
*                                                                       
*  (mu)_W__  negative charge                                            
*** triple vertex -- pfz(mu)=p12(mu),pwm(mu)=(-p1256(mu)),pwp(mu)=p56(mu
* ),efz=cz12(i1,i2),ewm=#,ewp=cw56,res=ctrip1256(i1?,i2?,m#)
      do mu=0,3
      vfz(mu)=(-p1256(mu))-p56(mu)
      vwm(mu)=p56(mu)-p12(mu)
      vwp(mu)=p12(mu)-(-p1256(mu))
      end do !mu
* vfz%efz
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cz12(i1,i2)%v,p=cz12(i1,i2)%e,q=vfz,bef=,aft=
      cz12(i1,i2)%v=(cz12(i1,i2)%e(0)*vfz(0)-cz12(i1,i2)%e(1)*vf
     & z(1)-cz12(i1,i2)%e(2)*vfz(2)-cz12(i1,i2)%e(3)*vfz(3))
      end do
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz12(i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cz12(i1,i2)%e(0)*cw56%e(0)-cz12(i1,i2)%e(1)*cw56%e(1
     & )-cz12(i1,i2)%e(2)*cw56%e(2)-cz12(i1,i2)%e(3)*cw56%e(3))
      do m=0,3
      ctrip1256(i1,i2,m)=cz12(i1,i2)%v*cw56%e(m)+vwm(m)*caux+cw5
     & 6%v*cz12(i1,i2)%e(m)
      end do
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do mu =0,3
            cw1256(i1,i2,mu)= cw1256(i1,i2,mu)
     &           +rcotw*ctrip1256(i1,i2,mu)
          enddo
        enddo
      enddo
*** triple vertex -- pfz(mu)=p12(mu),pwm(mu)=(-p1256(mu)),pwp(mu)=p56(mu
* ),efz=cf12(i1,i2),ewm=#,ewp=cw56,res=ctrip1256(i1?,i2?,m#)
      do mu=0,3
      vfz(mu)=(-p1256(mu))-p56(mu)
      vwm(mu)=p56(mu)-p12(mu)
      vwp(mu)=p12(mu)-(-p1256(mu))
      end do !mu
* vfz%efz
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cf12(i1,i2)%v,p=cf12(i1,i2)%e,q=vfz,bef=,aft=
      cf12(i1,i2)%v=(cf12(i1,i2)%e(0)*vfz(0)-cf12(i1,i2)%e(1)*vf
     & z(1)-cf12(i1,i2)%e(2)*vfz(2)-cf12(i1,i2)%e(3)*vfz(3))
      end do
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf12(i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cf12(i1,i2)%e(0)*cw56%e(0)-cf12(i1,i2)%e(1)*cw56%e(1
     & )-cf12(i1,i2)%e(2)*cw56%e(2)-cf12(i1,i2)%e(3)*cw56%e(3))
      do m=0,3
      ctrip1256(i1,i2,m)=cf12(i1,i2)%v*cw56%e(m)+vwm(m)*caux+cw5
     & 6%v*cf12(i1,i2)%e(m)
      end do
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do mu =0,3
            cw1256(i1,i2,mu)= cw1256(i1,i2,mu)
     &           +ctrip1256(i1,i2,mu)
          enddo
        enddo
      enddo
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=(p1278(mu)),pwp(mu)=p56(mu)
* ,efz=cz34(i1,i2),ewm=#,ewp=cw56,res=ctrip3456(i1?,i2?,m#)
      do mu=0,3
      vfz(mu)=(p1278(mu))-p56(mu)
      vwm(mu)=p56(mu)-p34(mu)
      vwp(mu)=p34(mu)-(p1278(mu))
      end do !mu
* vfz%efz
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cz34(i1,i2)%v,p=cz34(i1,i2)%e,q=vfz,bef=,aft=
      cz34(i1,i2)%v=(cz34(i1,i2)%e(0)*vfz(0)-cz34(i1,i2)%e(1)*vf
     & z(1)-cz34(i1,i2)%e(2)*vfz(2)-cz34(i1,i2)%e(3)*vfz(3))
      end do
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz34(i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cz34(i1,i2)%e(0)*cw56%e(0)-cz34(i1,i2)%e(1)*cw56%e(1
     & )-cz34(i1,i2)%e(2)*cw56%e(2)-cz34(i1,i2)%e(3)*cw56%e(3))
      do m=0,3
      ctrip3456(i1,i2,m)=cz34(i1,i2)%v*cw56%e(m)+vwm(m)*caux+cw5
     & 6%v*cz34(i1,i2)%e(m)
      end do
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do mu =0,3
            cw3456(i1,i2,mu)= cw3456(i1,i2,mu)
     &           +rcotw*ctrip3456(i1,i2,mu)
          enddo
        enddo
      enddo
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=(p1278(mu)),pwp(mu)=p56(mu)
* ,efz=cf34(i1,i2),ewm=#,ewp=cw56,res=ctrip3456(i1?,i2?,m#)
      do mu=0,3
      vfz(mu)=(p1278(mu))-p56(mu)
      vwm(mu)=p56(mu)-p34(mu)
      vwp(mu)=p34(mu)-(p1278(mu))
      end do !mu
* vfz%efz
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cf34(i1,i2)%v,p=cf34(i1,i2)%e,q=vfz,bef=,aft=
      cf34(i1,i2)%v=(cf34(i1,i2)%e(0)*vfz(0)-cf34(i1,i2)%e(1)*vf
     & z(1)-cf34(i1,i2)%e(2)*vfz(2)-cf34(i1,i2)%e(3)*vfz(3))
      end do
      end do
* vwp%ewp
* p.q -- p.q=cw56%v,p=cw56%e,q=vwp,bef=,aft=
      cw56%v=(cw56%e(0)*vwp(0)-cw56%e(1)*vwp(1)-cw56%e(2)*vwp(2)
     & -cw56%e(3)*vwp(3))
* efz%ewp
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf34(i1,i2)%e,q=cw56%e,bef=,aft=
      caux=(cf34(i1,i2)%e(0)*cw56%e(0)-cf34(i1,i2)%e(1)*cw56%e(1
     & )-cf34(i1,i2)%e(2)*cw56%e(2)-cf34(i1,i2)%e(3)*cw56%e(3))
      do m=0,3
      ctrip3456(i1,i2,m)=cf34(i1,i2)%v*cw56%e(m)+vwm(m)*caux+cw5
     & 6%v*cf34(i1,i2)%e(m)
      end do
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do mu =0,3
            cw3456(i1,i2,mu)= cw3456(i1,i2,mu)
     &           +ctrip3456(i1,i2,mu)
          enddo
        enddo
      enddo
  
*  (mu)_W__  positive charge                                            
*** triple vertex -- pfz(mu)=p12(mu),pwm(mu)=p78(mu),pwp(mu)=(-p1278(mu)
* ),efz=cz12(i1,i2),ewm=cw78,ewp=#,res=ctrip1278(i1?,i2?,m#)
      do mu=0,3
      vfz(mu)=p78(mu)-(-p1278(mu))
      vwm(mu)=(-p1278(mu))-p12(mu)
      vwp(mu)=p12(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cz12(i1,i2)%v,p=cz12(i1,i2)%e,q=vfz,bef=,aft=
      cz12(i1,i2)%v=(cz12(i1,i2)%e(0)*vfz(0)-cz12(i1,i2)%e(1)*vf
     & z(1)-cz12(i1,i2)%e(2)*vfz(2)-cz12(i1,i2)%e(3)*vfz(3))
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* efz%ewm
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz12(i1,i2)%e,q=cw78%e,bef=,aft=
      caux=(cz12(i1,i2)%e(0)*cw78%e(0)-cz12(i1,i2)%e(1)*cw78%e(1
     & )-cz12(i1,i2)%e(2)*cw78%e(2)-cz12(i1,i2)%e(3)*cw78%e(3))
      do m=0,3
      ctrip1278(i1,i2,m)=cz12(i1,i2)%v*cw78%e(m)+cw78%v*cz12(i1,
     & i2)%e(m)+vwp(m)*caux
      end do
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do mu =0,3
            cw1278(i1,i2,mu)= cw1278(i1,i2,mu)
     &           +rcotw*ctrip1278(i1,i2,mu)
          enddo
        enddo
      enddo
*** triple vertex -- pfz(mu)=p12(mu),pwm(mu)=p78(mu),pwp(mu)=(-p1278(mu)
* ),efz=cf12(i1,i2),ewm=cw78,ewp=#,res=ctrip1278(i1?,i2?,m#)
      do mu=0,3
      vfz(mu)=p78(mu)-(-p1278(mu))
      vwm(mu)=(-p1278(mu))-p12(mu)
      vwp(mu)=p12(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cf12(i1,i2)%v,p=cf12(i1,i2)%e,q=vfz,bef=,aft=
      cf12(i1,i2)%v=(cf12(i1,i2)%e(0)*vfz(0)-cf12(i1,i2)%e(1)*vf
     & z(1)-cf12(i1,i2)%e(2)*vfz(2)-cf12(i1,i2)%e(3)*vfz(3))
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* efz%ewm
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf12(i1,i2)%e,q=cw78%e,bef=,aft=
      caux=(cf12(i1,i2)%e(0)*cw78%e(0)-cf12(i1,i2)%e(1)*cw78%e(1
     & )-cf12(i1,i2)%e(2)*cw78%e(2)-cf12(i1,i2)%e(3)*cw78%e(3))
      do m=0,3
      ctrip1278(i1,i2,m)=cf12(i1,i2)%v*cw78%e(m)+cw78%v*cf12(i1,
     & i2)%e(m)+vwp(m)*caux
      end do
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do mu =0,3
            cw1278(i1,i2,mu)= cw1278(i1,i2,mu)
     &           +ctrip1278(i1,i2,mu)
          enddo
        enddo
      enddo
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p78(mu),pwp(mu)=(p1256(mu))
* ,efz=cz34(i1,i2),ewm=cw78,ewp=#,res=ctrip3478(i1?,i2?,m#)
      do mu=0,3
      vfz(mu)=p78(mu)-(p1256(mu))
      vwm(mu)=(p1256(mu))-p34(mu)
      vwp(mu)=p34(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cz34(i1,i2)%v,p=cz34(i1,i2)%e,q=vfz,bef=,aft=
      cz34(i1,i2)%v=(cz34(i1,i2)%e(0)*vfz(0)-cz34(i1,i2)%e(1)*vf
     & z(1)-cz34(i1,i2)%e(2)*vfz(2)-cz34(i1,i2)%e(3)*vfz(3))
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* efz%ewm
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cz34(i1,i2)%e,q=cw78%e,bef=,aft=
      caux=(cz34(i1,i2)%e(0)*cw78%e(0)-cz34(i1,i2)%e(1)*cw78%e(1
     & )-cz34(i1,i2)%e(2)*cw78%e(2)-cz34(i1,i2)%e(3)*cw78%e(3))
      do m=0,3
      ctrip3478(i1,i2,m)=cz34(i1,i2)%v*cw78%e(m)+cw78%v*cz34(i1,
     & i2)%e(m)+vwp(m)*caux
      end do
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do mu =0,3
            cw3478(i1,i2,mu)= cw3478(i1,i2,mu)
     &           +rcotw*ctrip3478(i1,i2,mu)
          enddo
        enddo
      enddo
*** triple vertex -- pfz(mu)=p34(mu),pwm(mu)=p78(mu),pwp(mu)=(p1256(mu))
* ,efz=cf34(i1,i2),ewm=cw78,ewp=#,res=ctrip3478(i1?,i2?,m#)
      do mu=0,3
      vfz(mu)=p78(mu)-(p1256(mu))
      vwm(mu)=(p1256(mu))-p34(mu)
      vwp(mu)=p34(mu)-p78(mu)
      end do !mu
* vfz%efz
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cf34(i1,i2)%v,p=cf34(i1,i2)%e,q=vfz,bef=,aft=
      cf34(i1,i2)%v=(cf34(i1,i2)%e(0)*vfz(0)-cf34(i1,i2)%e(1)*vf
     & z(1)-cf34(i1,i2)%e(2)*vfz(2)-cf34(i1,i2)%e(3)*vfz(3))
      end do
      end do
* vwm%ewm
* p.q -- p.q=cw78%v,p=cw78%e,q=vwm,bef=,aft=
      cw78%v=(cw78%e(0)*vwm(0)-cw78%e(1)*vwm(1)-cw78%e(2)*vwm(2)
     & -cw78%e(3)*vwm(3))
* efz%ewm
      do i1=1,2
      do i2=1,2
* p.q -- p.q=caux,p=cf34(i1,i2)%e,q=cw78%e,bef=,aft=
      caux=(cf34(i1,i2)%e(0)*cw78%e(0)-cf34(i1,i2)%e(1)*cw78%e(1
     & )-cf34(i1,i2)%e(2)*cw78%e(2)-cf34(i1,i2)%e(3)*cw78%e(3))
      do m=0,3
      ctrip3478(i1,i2,m)=cf34(i1,i2)%v*cw78%e(m)+cw78%v*cf34(i1,
     & i2)%e(m)+vwp(m)*caux
      end do
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do mu =0,3
            cw3478(i1,i2,mu)= cw3478(i1,i2,mu)
     &           +ctrip3478(i1,i2,mu)
          enddo
        enddo
      enddo
  
*    And finally the triple vertex                                      
*                                                                       
*                    /_                                                 
*                  W/                                                   
*         (mu)_W___/                                                    
*                  \                                                    
*                  h\_                                                  
*                    \                                                  
*                                                                       
  
* rmh < 0 in our convention means no higgs coupling                     
      if(id1.eq.5.and.rmh.ge.0.d0) then
        do i1=1,2
          do i2=1,2
            do mu =0,3
              cw1256(i1,i2,mu)= cw1256(i1,i2,mu)
     &             +ch12(i1,i2)*cw56%e(mu)*
     &       (rhww*rhbb/(s12-cmh2)+rhhww*rhhbb/(s12-cmhh2))
  
            enddo
          enddo
        enddo
      endif
  
* rmh < 0 in our convention means no higgs coupling                     
      if(id3.eq.5.and.rmh.ge.0.d0) then
        do i1=1,2
          do i2=1,2
            do mu =0,3
              cw3456(i1,i2,mu)= cw3456(i1,i2,mu)
     &             +ch34(i1,i2)*cw56%e(mu)*
     &       (rhww*rhbb/(s34-cmh2)+rhhww*rhhbb/(s34-cmhh2))
  
            enddo
          enddo
        enddo
      endif
  
* rmh < 0 in our convention means no higgs coupling                     
      if(id1.eq.5.and.rmh.ge.0.d0) then
        do i1=1,2
          do i2=1,2
            do mu =0,3
              cw1278(i1,i2,mu)= cw1278(i1,i2,mu)
     &             +ch12(i1,i2)*cw78%e(mu)*
     &       (rhww*rhbb/(s12-cmh2)+rhhww*rhhbb/(s12-cmhh2))
  
            enddo
          enddo
        enddo
      endif
  
* rmh < 0 in our convention means no higgs coupling                     
      if(id3.eq.5.and.rmh.ge.0.d0) then
        do i1=1,2
          do i2=1,2
            do mu =0,3
              cw3478(i1,i2,mu)= cw3478(i1,i2,mu)
     &             +ch34(i1,i2)*cw78%e(mu)*
     &       (rhww*rhbb/(s34-cmh2)+rhhww*rhhbb/(s34-cmhh2))
  
            enddo
          enddo
        enddo
      endif
  
* compute all subdiagrams with a Z or a gamma "decaying" to 4 fermions  
*  which correspond to two z's                                          
*                  czijkl(mu)  and cfijkl(mu) (ijkl particle numbers)   
*                                                                       
*               _____ i                     _____ i                     
*              |  |__ j                    |  |__ j                     
*    (mu) --Z--|  |__ k          (mu) --f--|  |__ k                     
*              |__|__ l                    |__|__ l                     
*                                                                       
*                                                                       
  
  
*    First all diagrams of the type                                     
*                                                                       
*                 |_Z,f,h__/                        |_Z,f,h__/          
*       (mu) _Z___|        \    and    (mu) _gamma__|        \          
*                 |                                 |                   
  
* Use tst because Z,f,h correspond to ts                                
  
  
      do i3=1,2
      do i4=1,2
      do mu=0,3
* TsT -- aa=clinetzz(i3,i4,mu)%a,bb=clinetzz(i3,i4,mu)%b,cc=clinetzz(i3,
* i4,mu)%c,dd=clinetzz(i3,i4,mu)%d,a1=l1_34(i3,i4)%a,b1=l1_34(i3,i4)%b,c
* 1=l1_34(i3,i4)%c,d1=l1_34(i3,i4)%d,a2=rz2_134(mu)%a,b2=rz2_134(mu)%b,c
* 2=rz2_134(mu)%c,d2=rz2_134(mu)%d,prq=p134q,m=rmass(id1),nsum=0
      do iut=1,2
      cx1=l1_34(i3,i4)%a(iut,1)+l1_34(i3,i4)%c(iut,1)*rmass(id1)
      cx2=l1_34(i3,i4)%c(iut,2)*p134q+l1_34(i3,i4)%a(iut,2)*rmas
     & s(id1)
      cy1=l1_34(i3,i4)%b(iut,1)+l1_34(i3,i4)%d(iut,1)*rmass(id1)
      cy2=l1_34(i3,i4)%d(iut,2)*p134q+l1_34(i3,i4)%b(iut,2)*rmas
     & s(id1)
      cw1=l1_34(i3,i4)%c(iut,1)*p134q+l1_34(i3,i4)%a(iut,1)*rmas
     & s(id1)
      cw2=l1_34(i3,i4)%a(iut,2)+l1_34(i3,i4)%c(iut,2)*rmass(id1)
      cz1=l1_34(i3,i4)%d(iut,1)*p134q+l1_34(i3,i4)%b(iut,1)*rmas
     & s(id1)
      cz2=l1_34(i3,i4)%b(iut,2)+l1_34(i3,i4)%d(iut,2)*rmass(id1)
      clinetzz(i3,i4,mu)%a(iut,1)=cx1*rz2_134(mu)%a(1,1)+cx2*rz2
     & _134(mu)%b(2,1)
      clinetzz(i3,i4,mu)%b(iut,1)=cy1*rz2_134(mu)%a(1,1)+cy2*rz2
     & _134(mu)%b(2,1)
      clinetzz(i3,i4,mu)%c(iut,1)=cw1*rz2_134(mu)%d(1,1)+cw2*rz2
     & _134(mu)%c(2,1)
      clinetzz(i3,i4,mu)%d(iut,1)=cz1*rz2_134(mu)%d(1,1)+cz2*rz2
     & _134(mu)%c(2,1)
      clinetzz(i3,i4,mu)%a(iut,2)=cw1*rz2_134(mu)%b(1,2)+cw2*rz2
     & _134(mu)%a(2,2)
      clinetzz(i3,i4,mu)%b(iut,2)=cz1*rz2_134(mu)%b(1,2)+cz2*rz2
     & _134(mu)%a(2,2)
      clinetzz(i3,i4,mu)%c(iut,2)=cx1*rz2_134(mu)%c(1,2)+cx2*rz2
     & _134(mu)%d(2,2)
      clinetzz(i3,i4,mu)%d(iut,2)=cy1*rz2_134(mu)%c(1,2)+cy2*rz2
     & _134(mu)%d(2,2)
      end do
      end do
      end do
      end do
  
*    Then those of the type                                             
*                                                                       
*       (mu) _Z___|                      (mu) _f__|                     
*                 |_Z,f,h__/     and              |_Z,f,h__/            
*                 |        \                      |        \            
*                                                                       
* Use tts because Z,f,h correspond to ts                                
  
  
      do i3=1,2
      do i4=1,2
      do mu=0,3
* TTs -- aa=clinetzz(i3,i4,mu)%a,bb=clinetzz(i3,i4,mu)%b,cc=clinetzz(i3,
* i4,mu)%c,dd=clinetzz(i3,i4,mu)%d,a1=lz1_234(mu)%a,b1=lz1_234(mu)%b,c1=
* lz1_234(mu)%c,d1=lz1_234(mu)%d,a2=r2_34(i3,i4)%a,b2=r2_34(i3,i4)%b,c2=
* r2_34(i3,i4)%c,d2=r2_34(i3,i4)%d,prq=p234q,m=rmass(id1),nsum=1
      do jut=1,2
      cx1=r2_34(i3,i4)%a(1,jut)+rmass(id1)*r2_34(i3,i4)%b(1,jut)
      cx2=r2_34(i3,i4)%a(2,jut)+rmass(id1)*r2_34(i3,i4)%b(2,jut)
      cy1=p234q*r2_34(i3,i4)%b(1,jut)+rmass(id1)*r2_34(i3,i4)%a(
     & 1,jut)
      cy2=p234q*r2_34(i3,i4)%b(2,jut)+rmass(id1)*r2_34(i3,i4)%a(
     & 2,jut)
      cw1=r2_34(i3,i4)%c(1,jut)+rmass(id1)*r2_34(i3,i4)%d(1,jut)
      cw2=r2_34(i3,i4)%c(2,jut)+rmass(id1)*r2_34(i3,i4)%d(2,jut)
      cz1=p234q*r2_34(i3,i4)%d(1,jut)+rmass(id1)*r2_34(i3,i4)%c(
     & 1,jut)
      cz2=p234q*r2_34(i3,i4)%d(2,jut)+rmass(id1)*r2_34(i3,i4)%c(
     & 2,jut)
      clinetzz(i3,i4,mu)%a(1,jut)=clinetzz(i3,i4,mu)%a(1,jut)+lz
     & 1_234(mu)%a(1,1)*cx1+lz1_234(mu)%c(1,2)*cy2
      clinetzz(i3,i4,mu)%b(1,jut)=clinetzz(i3,i4,mu)%b(1,jut)+lz
     & 1_234(mu)%d(1,1)*cy1+lz1_234(mu)%b(1,2)*cx2
      clinetzz(i3,i4,mu)%c(1,jut)=clinetzz(i3,i4,mu)%c(1,jut)+lz
     & 1_234(mu)%a(1,1)*cw1+lz1_234(mu)%c(1,2)*cz2
      clinetzz(i3,i4,mu)%d(1,jut)=clinetzz(i3,i4,mu)%d(1,jut)+lz
     & 1_234(mu)%d(1,1)*cz1+lz1_234(mu)%b(1,2)*cw2
      clinetzz(i3,i4,mu)%a(2,jut)=clinetzz(i3,i4,mu)%a(2,jut)+lz
     & 1_234(mu)%c(2,1)*cy1+lz1_234(mu)%a(2,2)*cx2
      clinetzz(i3,i4,mu)%b(2,jut)=clinetzz(i3,i4,mu)%b(2,jut)+lz
     & 1_234(mu)%b(2,1)*cx1+lz1_234(mu)%d(2,2)*cy2
      clinetzz(i3,i4,mu)%c(2,jut)=clinetzz(i3,i4,mu)%c(2,jut)+lz
     & 1_234(mu)%c(2,1)*cz1+lz1_234(mu)%a(2,2)*cw2
      clinetzz(i3,i4,mu)%d(2,jut)=clinetzz(i3,i4,mu)%d(2,jut)+lz
     & 1_234(mu)%b(2,1)*cw1+lz1_234(mu)%d(2,2)*cz2
      end do
      end do
      end do
      end do
  
      rmassl=rmass(id1)
      rmassr=-rmass(id2)
  
      do i3=1,2
      do i4=1,2
      do mu=0,3
* mline -- res=cz1234(&1,&2,i3,i4,mu),abcd=clinetzz(i3,i4,mu)%,m1=rmassl
* ,m2=rmassr,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cz1234(iut,jut,i3,i4,mu)=clinetzz(i3,i4,mu)%a(iut,jut)+rma
     & ssl*clinetzz(i3,i4,mu)%b(iut,jut)+rmassr*clinetzz(i3,i4,m
     & u)%c(iut,jut)+rmassl*rmassr*clinetzz(i3,i4,mu)%d(iut,jut)
      enddo
      enddo
      end do
      end do
      end do
  
  
*    First all diagrams of the type                                     
*                                                                       
*                 |_Z,f,h__/                        |_Z,f,h__/          
*       (mu) _Z___|        \    and    (mu) _gamma__|        \          
*                 |                                 |                   
  
* Use tst because Z,f,h correspond to ts                                
  
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TsT -- aa=clinetzz(i1,i2,mu)%a,bb=clinetzz(i1,i2,mu)%b,cc=clinetzz(i1,
* i2,mu)%c,dd=clinetzz(i1,i2,mu)%d,a1=l3_12(i1,i2)%a,b1=l3_12(i1,i2)%b,c
* 1=l3_12(i1,i2)%c,d1=l3_12(i1,i2)%d,a2=rz4_312(mu)%a,b2=rz4_312(mu)%b,c
* 2=rz4_312(mu)%c,d2=rz4_312(mu)%d,prq=p312q,m=rmass(id3),nsum=0
      do iut=1,2
      cx1=l3_12(i1,i2)%a(iut,1)+l3_12(i1,i2)%c(iut,1)*rmass(id3)
      cx2=l3_12(i1,i2)%c(iut,2)*p312q+l3_12(i1,i2)%a(iut,2)*rmas
     & s(id3)
      cy1=l3_12(i1,i2)%b(iut,1)+l3_12(i1,i2)%d(iut,1)*rmass(id3)
      cy2=l3_12(i1,i2)%d(iut,2)*p312q+l3_12(i1,i2)%b(iut,2)*rmas
     & s(id3)
      cw1=l3_12(i1,i2)%c(iut,1)*p312q+l3_12(i1,i2)%a(iut,1)*rmas
     & s(id3)
      cw2=l3_12(i1,i2)%a(iut,2)+l3_12(i1,i2)%c(iut,2)*rmass(id3)
      cz1=l3_12(i1,i2)%d(iut,1)*p312q+l3_12(i1,i2)%b(iut,1)*rmas
     & s(id3)
      cz2=l3_12(i1,i2)%b(iut,2)+l3_12(i1,i2)%d(iut,2)*rmass(id3)
      clinetzz(i1,i2,mu)%a(iut,1)=cx1*rz4_312(mu)%a(1,1)+cx2*rz4
     & _312(mu)%b(2,1)
      clinetzz(i1,i2,mu)%b(iut,1)=cy1*rz4_312(mu)%a(1,1)+cy2*rz4
     & _312(mu)%b(2,1)
      clinetzz(i1,i2,mu)%c(iut,1)=cw1*rz4_312(mu)%d(1,1)+cw2*rz4
     & _312(mu)%c(2,1)
      clinetzz(i1,i2,mu)%d(iut,1)=cz1*rz4_312(mu)%d(1,1)+cz2*rz4
     & _312(mu)%c(2,1)
      clinetzz(i1,i2,mu)%a(iut,2)=cw1*rz4_312(mu)%b(1,2)+cw2*rz4
     & _312(mu)%a(2,2)
      clinetzz(i1,i2,mu)%b(iut,2)=cz1*rz4_312(mu)%b(1,2)+cz2*rz4
     & _312(mu)%a(2,2)
      clinetzz(i1,i2,mu)%c(iut,2)=cx1*rz4_312(mu)%c(1,2)+cx2*rz4
     & _312(mu)%d(2,2)
      clinetzz(i1,i2,mu)%d(iut,2)=cy1*rz4_312(mu)%c(1,2)+cy2*rz4
     & _312(mu)%d(2,2)
      end do
      end do
      end do
      end do
  
*    Then those of the type                                             
*                                                                       
*       (mu) _Z___|                      (mu) _f__|                     
*                 |_Z,f,h__/     and              |_Z,f,h__/            
*                 |        \                      |        \            
*                                                                       
* Use tts because Z,f,h correspond to ts                                
  
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TTs -- aa=clinetzz(i1,i2,mu)%a,bb=clinetzz(i1,i2,mu)%b,cc=clinetzz(i1,
* i2,mu)%c,dd=clinetzz(i1,i2,mu)%d,a1=lz3_412(mu)%a,b1=lz3_412(mu)%b,c1=
* lz3_412(mu)%c,d1=lz3_412(mu)%d,a2=r4_12(i1,i2)%a,b2=r4_12(i1,i2)%b,c2=
* r4_12(i1,i2)%c,d2=r4_12(i1,i2)%d,prq=p412q,m=rmass(id3),nsum=1
      do jut=1,2
      cx1=r4_12(i1,i2)%a(1,jut)+rmass(id3)*r4_12(i1,i2)%b(1,jut)
      cx2=r4_12(i1,i2)%a(2,jut)+rmass(id3)*r4_12(i1,i2)%b(2,jut)
      cy1=p412q*r4_12(i1,i2)%b(1,jut)+rmass(id3)*r4_12(i1,i2)%a(
     & 1,jut)
      cy2=p412q*r4_12(i1,i2)%b(2,jut)+rmass(id3)*r4_12(i1,i2)%a(
     & 2,jut)
      cw1=r4_12(i1,i2)%c(1,jut)+rmass(id3)*r4_12(i1,i2)%d(1,jut)
      cw2=r4_12(i1,i2)%c(2,jut)+rmass(id3)*r4_12(i1,i2)%d(2,jut)
      cz1=p412q*r4_12(i1,i2)%d(1,jut)+rmass(id3)*r4_12(i1,i2)%c(
     & 1,jut)
      cz2=p412q*r4_12(i1,i2)%d(2,jut)+rmass(id3)*r4_12(i1,i2)%c(
     & 2,jut)
      clinetzz(i1,i2,mu)%a(1,jut)=clinetzz(i1,i2,mu)%a(1,jut)+lz
     & 3_412(mu)%a(1,1)*cx1+lz3_412(mu)%c(1,2)*cy2
      clinetzz(i1,i2,mu)%b(1,jut)=clinetzz(i1,i2,mu)%b(1,jut)+lz
     & 3_412(mu)%d(1,1)*cy1+lz3_412(mu)%b(1,2)*cx2
      clinetzz(i1,i2,mu)%c(1,jut)=clinetzz(i1,i2,mu)%c(1,jut)+lz
     & 3_412(mu)%a(1,1)*cw1+lz3_412(mu)%c(1,2)*cz2
      clinetzz(i1,i2,mu)%d(1,jut)=clinetzz(i1,i2,mu)%d(1,jut)+lz
     & 3_412(mu)%d(1,1)*cz1+lz3_412(mu)%b(1,2)*cw2
      clinetzz(i1,i2,mu)%a(2,jut)=clinetzz(i1,i2,mu)%a(2,jut)+lz
     & 3_412(mu)%c(2,1)*cy1+lz3_412(mu)%a(2,2)*cx2
      clinetzz(i1,i2,mu)%b(2,jut)=clinetzz(i1,i2,mu)%b(2,jut)+lz
     & 3_412(mu)%b(2,1)*cx1+lz3_412(mu)%d(2,2)*cy2
      clinetzz(i1,i2,mu)%c(2,jut)=clinetzz(i1,i2,mu)%c(2,jut)+lz
     & 3_412(mu)%c(2,1)*cz1+lz3_412(mu)%a(2,2)*cw2
      clinetzz(i1,i2,mu)%d(2,jut)=clinetzz(i1,i2,mu)%d(2,jut)+lz
     & 3_412(mu)%b(2,1)*cw1+lz3_412(mu)%d(2,2)*cz2
      end do
      end do
      end do
      end do
  
      rmassl=rmass(id3)
      rmassr=-rmass(id4)
  
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* mline -- res=cz1234(i1,i2,&1,&2,mu),abcd=clinetzz(i1,i2,mu)%,m1=rmassl
* ,m2=rmassr,den=0,nsum=1
      do iut=1,2
      do jut=1,2
      cz1234(i1,i2,iut,jut,mu)=cz1234(i1,i2,iut,jut,mu)+clinetzz
     & (i1,i2,mu)%a(iut,jut)+rmassl*clinetzz(i1,i2,mu)%b(iut,jut
     & )+rmassr*clinetzz(i1,i2,mu)%c(iut,jut)+rmassl*rmassr*clin
     & etzz(i1,i2,mu)%d(iut,jut)
      enddo
      enddo
      end do
      end do
      end do
  
*    First all diagrams of the type                                     
*                                                                       
*                 |_Z,f,h__/                        |_Z,f,h__/          
*       (mu) _Z___|        \    and    (mu) _gamma__|        \          
*                 |                                 |                   
  
* Use tst because Z,f,h correspond to ts                                
  
  
      do i3=1,2
      do i4=1,2
      do mu=0,3
* TsT -- aa=clinetfz(i3,i4,mu)%a,bb=clinetfz(i3,i4,mu)%b,cc=clinetfz(i3,
* i4,mu)%c,dd=clinetfz(i3,i4,mu)%d,a1=l1_34(i3,i4)%a,b1=l1_34(i3,i4)%b,c
* 1=l1_34(i3,i4)%c,d1=l1_34(i3,i4)%d,a2=rf2_134(mu)%a,b2=rf2_134(mu)%b,c
* 2=rf2_134(mu)%c,d2=rf2_134(mu)%d,prq=p134q,m=rmass(id1),nsum=0
      do iut=1,2
      cx1=l1_34(i3,i4)%a(iut,1)+l1_34(i3,i4)%c(iut,1)*rmass(id1)
      cx2=l1_34(i3,i4)%c(iut,2)*p134q+l1_34(i3,i4)%a(iut,2)*rmas
     & s(id1)
      cy1=l1_34(i3,i4)%b(iut,1)+l1_34(i3,i4)%d(iut,1)*rmass(id1)
      cy2=l1_34(i3,i4)%d(iut,2)*p134q+l1_34(i3,i4)%b(iut,2)*rmas
     & s(id1)
      cw1=l1_34(i3,i4)%c(iut,1)*p134q+l1_34(i3,i4)%a(iut,1)*rmas
     & s(id1)
      cw2=l1_34(i3,i4)%a(iut,2)+l1_34(i3,i4)%c(iut,2)*rmass(id1)
      cz1=l1_34(i3,i4)%d(iut,1)*p134q+l1_34(i3,i4)%b(iut,1)*rmas
     & s(id1)
      cz2=l1_34(i3,i4)%b(iut,2)+l1_34(i3,i4)%d(iut,2)*rmass(id1)
      clinetfz(i3,i4,mu)%a(iut,1)=cx1*rf2_134(mu)%a(1,1)+cx2*rf2
     & _134(mu)%b(2,1)
      clinetfz(i3,i4,mu)%b(iut,1)=cy1*rf2_134(mu)%a(1,1)+cy2*rf2
     & _134(mu)%b(2,1)
      clinetfz(i3,i4,mu)%c(iut,1)=cw1*rf2_134(mu)%d(1,1)+cw2*rf2
     & _134(mu)%c(2,1)
      clinetfz(i3,i4,mu)%d(iut,1)=cz1*rf2_134(mu)%d(1,1)+cz2*rf2
     & _134(mu)%c(2,1)
      clinetfz(i3,i4,mu)%a(iut,2)=cw1*rf2_134(mu)%b(1,2)+cw2*rf2
     & _134(mu)%a(2,2)
      clinetfz(i3,i4,mu)%b(iut,2)=cz1*rf2_134(mu)%b(1,2)+cz2*rf2
     & _134(mu)%a(2,2)
      clinetfz(i3,i4,mu)%c(iut,2)=cx1*rf2_134(mu)%c(1,2)+cx2*rf2
     & _134(mu)%d(2,2)
      clinetfz(i3,i4,mu)%d(iut,2)=cy1*rf2_134(mu)%c(1,2)+cy2*rf2
     & _134(mu)%d(2,2)
      end do
      end do
      end do
      end do
  
*    Then those of the type                                             
*                                                                       
*       (mu) _Z___|                      (mu) _f__|                     
*                 |_Z,f,h__/     and              |_Z,f,h__/            
*                 |        \                      |        \            
*                                                                       
* Use tts because Z,f,h correspond to ts                                
  
  
      do i3=1,2
      do i4=1,2
      do mu=0,3
* TTs -- aa=clinetfz(i3,i4,mu)%a,bb=clinetfz(i3,i4,mu)%b,cc=clinetfz(i3,
* i4,mu)%c,dd=clinetfz(i3,i4,mu)%d,a1=lf1_234(mu)%a,b1=lf1_234(mu)%b,c1=
* lf1_234(mu)%c,d1=lf1_234(mu)%d,a2=r2_34(i3,i4)%a,b2=r2_34(i3,i4)%b,c2=
* r2_34(i3,i4)%c,d2=r2_34(i3,i4)%d,prq=p234q,m=rmass(id1),nsum=1
      do jut=1,2
      cx1=r2_34(i3,i4)%a(1,jut)+rmass(id1)*r2_34(i3,i4)%b(1,jut)
      cx2=r2_34(i3,i4)%a(2,jut)+rmass(id1)*r2_34(i3,i4)%b(2,jut)
      cy1=p234q*r2_34(i3,i4)%b(1,jut)+rmass(id1)*r2_34(i3,i4)%a(
     & 1,jut)
      cy2=p234q*r2_34(i3,i4)%b(2,jut)+rmass(id1)*r2_34(i3,i4)%a(
     & 2,jut)
      cw1=r2_34(i3,i4)%c(1,jut)+rmass(id1)*r2_34(i3,i4)%d(1,jut)
      cw2=r2_34(i3,i4)%c(2,jut)+rmass(id1)*r2_34(i3,i4)%d(2,jut)
      cz1=p234q*r2_34(i3,i4)%d(1,jut)+rmass(id1)*r2_34(i3,i4)%c(
     & 1,jut)
      cz2=p234q*r2_34(i3,i4)%d(2,jut)+rmass(id1)*r2_34(i3,i4)%c(
     & 2,jut)
      clinetfz(i3,i4,mu)%a(1,jut)=clinetfz(i3,i4,mu)%a(1,jut)+lf
     & 1_234(mu)%a(1,1)*cx1+lf1_234(mu)%c(1,2)*cy2
      clinetfz(i3,i4,mu)%b(1,jut)=clinetfz(i3,i4,mu)%b(1,jut)+lf
     & 1_234(mu)%d(1,1)*cy1+lf1_234(mu)%b(1,2)*cx2
      clinetfz(i3,i4,mu)%c(1,jut)=clinetfz(i3,i4,mu)%c(1,jut)+lf
     & 1_234(mu)%a(1,1)*cw1+lf1_234(mu)%c(1,2)*cz2
      clinetfz(i3,i4,mu)%d(1,jut)=clinetfz(i3,i4,mu)%d(1,jut)+lf
     & 1_234(mu)%d(1,1)*cz1+lf1_234(mu)%b(1,2)*cw2
      clinetfz(i3,i4,mu)%a(2,jut)=clinetfz(i3,i4,mu)%a(2,jut)+lf
     & 1_234(mu)%c(2,1)*cy1+lf1_234(mu)%a(2,2)*cx2
      clinetfz(i3,i4,mu)%b(2,jut)=clinetfz(i3,i4,mu)%b(2,jut)+lf
     & 1_234(mu)%b(2,1)*cx1+lf1_234(mu)%d(2,2)*cy2
      clinetfz(i3,i4,mu)%c(2,jut)=clinetfz(i3,i4,mu)%c(2,jut)+lf
     & 1_234(mu)%c(2,1)*cz1+lf1_234(mu)%a(2,2)*cw2
      clinetfz(i3,i4,mu)%d(2,jut)=clinetfz(i3,i4,mu)%d(2,jut)+lf
     & 1_234(mu)%b(2,1)*cw1+lf1_234(mu)%d(2,2)*cz2
      end do
      end do
      end do
      end do
  
      rmassl=rmass(id1)
      rmassr=-rmass(id2)
  
      do i3=1,2
      do i4=1,2
      do mu=0,3
* mline -- res=cf1234(&1,&2,i3,i4,mu),abcd=clinetfz(i3,i4,mu)%,m1=rmassl
* ,m2=rmassr,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      cf1234(iut,jut,i3,i4,mu)=clinetfz(i3,i4,mu)%a(iut,jut)+rma
     & ssl*clinetfz(i3,i4,mu)%b(iut,jut)+rmassr*clinetfz(i3,i4,m
     & u)%c(iut,jut)+rmassl*rmassr*clinetfz(i3,i4,mu)%d(iut,jut)
      enddo
      enddo
      end do
      end do
      end do
  
  
*    First all diagrams of the type                                     
*                                                                       
*                 |_Z,f,h__/                        |_Z,f,h__/          
*       (mu) _Z___|        \    and    (mu) _gamma__|        \          
*                 |                                 |                   
  
* Use tst because Z,f,h correspond to ts                                
  
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TsT -- aa=clinetfz(i1,i2,mu)%a,bb=clinetfz(i1,i2,mu)%b,cc=clinetfz(i1,
* i2,mu)%c,dd=clinetfz(i1,i2,mu)%d,a1=l3_12(i1,i2)%a,b1=l3_12(i1,i2)%b,c
* 1=l3_12(i1,i2)%c,d1=l3_12(i1,i2)%d,a2=rf4_312(mu)%a,b2=rf4_312(mu)%b,c
* 2=rf4_312(mu)%c,d2=rf4_312(mu)%d,prq=p312q,m=rmass(id3),nsum=0
      do iut=1,2
      cx1=l3_12(i1,i2)%a(iut,1)+l3_12(i1,i2)%c(iut,1)*rmass(id3)
      cx2=l3_12(i1,i2)%c(iut,2)*p312q+l3_12(i1,i2)%a(iut,2)*rmas
     & s(id3)
      cy1=l3_12(i1,i2)%b(iut,1)+l3_12(i1,i2)%d(iut,1)*rmass(id3)
      cy2=l3_12(i1,i2)%d(iut,2)*p312q+l3_12(i1,i2)%b(iut,2)*rmas
     & s(id3)
      cw1=l3_12(i1,i2)%c(iut,1)*p312q+l3_12(i1,i2)%a(iut,1)*rmas
     & s(id3)
      cw2=l3_12(i1,i2)%a(iut,2)+l3_12(i1,i2)%c(iut,2)*rmass(id3)
      cz1=l3_12(i1,i2)%d(iut,1)*p312q+l3_12(i1,i2)%b(iut,1)*rmas
     & s(id3)
      cz2=l3_12(i1,i2)%b(iut,2)+l3_12(i1,i2)%d(iut,2)*rmass(id3)
      clinetfz(i1,i2,mu)%a(iut,1)=cx1*rf4_312(mu)%a(1,1)+cx2*rf4
     & _312(mu)%b(2,1)
      clinetfz(i1,i2,mu)%b(iut,1)=cy1*rf4_312(mu)%a(1,1)+cy2*rf4
     & _312(mu)%b(2,1)
      clinetfz(i1,i2,mu)%c(iut,1)=cw1*rf4_312(mu)%d(1,1)+cw2*rf4
     & _312(mu)%c(2,1)
      clinetfz(i1,i2,mu)%d(iut,1)=cz1*rf4_312(mu)%d(1,1)+cz2*rf4
     & _312(mu)%c(2,1)
      clinetfz(i1,i2,mu)%a(iut,2)=cw1*rf4_312(mu)%b(1,2)+cw2*rf4
     & _312(mu)%a(2,2)
      clinetfz(i1,i2,mu)%b(iut,2)=cz1*rf4_312(mu)%b(1,2)+cz2*rf4
     & _312(mu)%a(2,2)
      clinetfz(i1,i2,mu)%c(iut,2)=cx1*rf4_312(mu)%c(1,2)+cx2*rf4
     & _312(mu)%d(2,2)
      clinetfz(i1,i2,mu)%d(iut,2)=cy1*rf4_312(mu)%c(1,2)+cy2*rf4
     & _312(mu)%d(2,2)
      end do
      end do
      end do
      end do
  
*    Then those of the type                                             
*                                                                       
*       (mu) _Z___|                      (mu) _f__|                     
*                 |_Z,f,h__/     and              |_Z,f,h__/            
*                 |        \                      |        \            
*                                                                       
* Use tts because Z,f,h correspond to ts                                
  
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* TTs -- aa=clinetfz(i1,i2,mu)%a,bb=clinetfz(i1,i2,mu)%b,cc=clinetfz(i1,
* i2,mu)%c,dd=clinetfz(i1,i2,mu)%d,a1=lf3_412(mu)%a,b1=lf3_412(mu)%b,c1=
* lf3_412(mu)%c,d1=lf3_412(mu)%d,a2=r4_12(i1,i2)%a,b2=r4_12(i1,i2)%b,c2=
* r4_12(i1,i2)%c,d2=r4_12(i1,i2)%d,prq=p412q,m=rmass(id3),nsum=1
      do jut=1,2
      cx1=r4_12(i1,i2)%a(1,jut)+rmass(id3)*r4_12(i1,i2)%b(1,jut)
      cx2=r4_12(i1,i2)%a(2,jut)+rmass(id3)*r4_12(i1,i2)%b(2,jut)
      cy1=p412q*r4_12(i1,i2)%b(1,jut)+rmass(id3)*r4_12(i1,i2)%a(
     & 1,jut)
      cy2=p412q*r4_12(i1,i2)%b(2,jut)+rmass(id3)*r4_12(i1,i2)%a(
     & 2,jut)
      cw1=r4_12(i1,i2)%c(1,jut)+rmass(id3)*r4_12(i1,i2)%d(1,jut)
      cw2=r4_12(i1,i2)%c(2,jut)+rmass(id3)*r4_12(i1,i2)%d(2,jut)
      cz1=p412q*r4_12(i1,i2)%d(1,jut)+rmass(id3)*r4_12(i1,i2)%c(
     & 1,jut)
      cz2=p412q*r4_12(i1,i2)%d(2,jut)+rmass(id3)*r4_12(i1,i2)%c(
     & 2,jut)
      clinetfz(i1,i2,mu)%a(1,jut)=clinetfz(i1,i2,mu)%a(1,jut)+lf
     & 3_412(mu)%a(1,1)*cx1+lf3_412(mu)%c(1,2)*cy2
      clinetfz(i1,i2,mu)%b(1,jut)=clinetfz(i1,i2,mu)%b(1,jut)+lf
     & 3_412(mu)%d(1,1)*cy1+lf3_412(mu)%b(1,2)*cx2
      clinetfz(i1,i2,mu)%c(1,jut)=clinetfz(i1,i2,mu)%c(1,jut)+lf
     & 3_412(mu)%a(1,1)*cw1+lf3_412(mu)%c(1,2)*cz2
      clinetfz(i1,i2,mu)%d(1,jut)=clinetfz(i1,i2,mu)%d(1,jut)+lf
     & 3_412(mu)%d(1,1)*cz1+lf3_412(mu)%b(1,2)*cw2
      clinetfz(i1,i2,mu)%a(2,jut)=clinetfz(i1,i2,mu)%a(2,jut)+lf
     & 3_412(mu)%c(2,1)*cy1+lf3_412(mu)%a(2,2)*cx2
      clinetfz(i1,i2,mu)%b(2,jut)=clinetfz(i1,i2,mu)%b(2,jut)+lf
     & 3_412(mu)%b(2,1)*cx1+lf3_412(mu)%d(2,2)*cy2
      clinetfz(i1,i2,mu)%c(2,jut)=clinetfz(i1,i2,mu)%c(2,jut)+lf
     & 3_412(mu)%c(2,1)*cz1+lf3_412(mu)%a(2,2)*cw2
      clinetfz(i1,i2,mu)%d(2,jut)=clinetfz(i1,i2,mu)%d(2,jut)+lf
     & 3_412(mu)%b(2,1)*cw1+lf3_412(mu)%d(2,2)*cz2
      end do
      end do
      end do
      end do
  
      rmassl=rmass(id3)
      rmassr=-rmass(id4)
  
  
      do i1=1,2
      do i2=1,2
      do mu=0,3
* mline -- res=cf1234(i1,i2,&1,&2,mu),abcd=clinetfz(i1,i2,mu)%,m1=rmassl
* ,m2=rmassr,den=0,nsum=1
      do iut=1,2
      do jut=1,2
      cf1234(i1,i2,iut,jut,mu)=cf1234(i1,i2,iut,jut,mu)+clinetfz
     & (i1,i2,mu)%a(iut,jut)+rmassl*clinetfz(i1,i2,mu)%b(iut,jut
     & )+rmassr*clinetfz(i1,i2,mu)%c(iut,jut)+rmassl*rmassr*clin
     & etfz(i1,i2,mu)%d(iut,jut)
      enddo
      enddo
      end do
      end do
      end do
  
*    Then those with triple vertex                                      
*                                                                       
*                    /_                                                 
*                  Z/                                                   
*         (mu)_Z___/                                                    
*                  \                                                    
*                  h\_                                                  
*                    \                                                  
*                                                                       
  
  
* rmh < 0 in our convention means no higgs coupling                     
      if(id1.eq.5.and.rmh.ge.0.d0) then
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                do mu =0,3
                  cz1234(i1,i2,i3,i4,mu)= cz1234(i1,i2,i3,i4,mu)
     &             +ch12(i1,i2)*cz34(i3,i4)%e(mu)*
     &        (rhzz*rhbb/(s12-cmh2)+rhhzz*rhhbb/(s12-cmhh2))
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
  
* rmh < 0 in our convention means no higgs coupling                     
      if(id3.eq.5.and.rmh.ge.0.d0) then
        do i1=1,2
          do i2=1,2
            do i3=1,2
              do i4=1,2
                do mu =0,3
                  cz1234(i1,i2,i3,i4,mu)= cz1234(i1,i2,i3,i4,mu)
     &             +ch34(i3,i4)*cz12(i1,i2)%e(mu)*
     &        (rhzz*rhbb/(s34-cmh2)+rhhzz*rhhbb/(s34-cmhh2))
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
  
  
* compute all subdiagrams with a higgs "decaying" to 4 fermions         
*  which correspond to two z's                                          
*                                                                       
*                   _____ i                                             
*                  |  |__ j                                             
*             --h--|  |__ k                                             
*                  |__|__ l                                             
*                                                                       
*                  ch1234(i1,i2,i3,i4) con t aux clineth(2,2)           
  
*                                                                       
*    First all diagrams of the type                                     
*                                                                       
*                   |_Z,f,h__/                                          
*              _h___|        \                                          
*                   |                                                   
  
* rmh < 0 in our convention means no higgs coupling                     
      if (id1.eq.5.and.rmh.ge.0.d0) then
      do i3=1,2
      do i4=1,2
* TsTSC -- aa=clineth12(i3,i4)%a,bb=clineth12(i3,i4)%b,cc=clineth12(i3,i
* 4)%c,dd=clineth12(i3,i4)%d,a1=l1_34(i3,i4)%a,b1=l1_34(i3,i4)%b,c1=l1_3
* 4(i3,i4)%c,d1=l1_34(i3,i4)%d,a2=rh2_134%a,b2=rh2_134%b,c2=rh2_134%c,pr
* q=p134q,m=rmass(id1),nsum=0
      do iut=1,2
      cx1=l1_34(i3,i4)%a(iut,1)+l1_34(i3,i4)%c(iut,1)*rmass(id1)
      cx2=l1_34(i3,i4)%c(iut,2)*p134q+l1_34(i3,i4)%a(iut,2)*rmas
     & s(id1)
      cy1=l1_34(i3,i4)%b(iut,1)+l1_34(i3,i4)%d(iut,1)*rmass(id1)
      cy2=l1_34(i3,i4)%d(iut,2)*p134q+l1_34(i3,i4)%b(iut,2)*rmas
     & s(id1)
      cw1=l1_34(i3,i4)%c(iut,1)*p134q+l1_34(i3,i4)%a(iut,1)*rmas
     & s(id1)
      cw2=l1_34(i3,i4)%a(iut,2)+l1_34(i3,i4)%c(iut,2)*rmass(id1)
      cz1=l1_34(i3,i4)%d(iut,1)*p134q+l1_34(i3,i4)%b(iut,1)*rmas
     & s(id1)
      cz2=l1_34(i3,i4)%b(iut,2)+l1_34(i3,i4)%d(iut,2)*rmass(id1)
      clineth12(i3,i4)%a(iut,2)=cx1*rh2_134%a(1,2)+cx2*rh2_134%b
     & (2,2)
      clineth12(i3,i4)%b(iut,2)=cy1*rh2_134%a(1,2)+cy2*rh2_134%b
     & (2,2)
      clineth12(i3,i4)%c(iut,2)=cw2*rh2_134%c(2,2)
      clineth12(i3,i4)%d(iut,2)=cz2*rh2_134%c(2,2)
      clineth12(i3,i4)%a(iut,1)=cw1*rh2_134%b(1,1)+cw2*rh2_134%a
     & (2,1)
      clineth12(i3,i4)%b(iut,1)=cz1*rh2_134%b(1,1)+cz2*rh2_134%a
     & (2,1)
      clineth12(i3,i4)%c(iut,1)=cx1*rh2_134%c(1,1)
      clineth12(i3,i4)%d(iut,1)=cy1*rh2_134%c(1,1)
      end do
      end do
      end do
  
  
*    Then those of the type                                             
*                                                                       
*              _h___|                                                   
*                   |_Z,f,h__/                                          
*                   |        \                                          
*                                                                       
  
      do i3=1,2
      do i4=1,2
* TSCTs -- aa=clineth12(i3,i4)%a,bb=clineth12(i3,i4)%b,cc=clineth12(i3,i
* 4)%c,dd=clineth12(i3,i4)%d,a1=lh1_234%a,b1=lh1_234%b,c1=lh1_234%c,a2=r
* 2_34(i3,i4)%a,b2=r2_34(i3,i4)%b,c2=r2_34(i3,i4)%c,d2=r2_34(i3,i4)%d,pr
* q=p234q,m=rmass(id1),nsum=1
      do jut=1,2
      cx1=r2_34(i3,i4)%a(1,jut)+rmass(id1)*r2_34(i3,i4)%b(1,jut)
      cx2=r2_34(i3,i4)%a(2,jut)+rmass(id1)*r2_34(i3,i4)%b(2,jut)
      cy1=p234q*r2_34(i3,i4)%b(1,jut)+rmass(id1)*r2_34(i3,i4)%a(
     & 1,jut)
      cy2=p234q*r2_34(i3,i4)%b(2,jut)+rmass(id1)*r2_34(i3,i4)%a(
     & 2,jut)
      cw1=r2_34(i3,i4)%c(1,jut)+rmass(id1)*r2_34(i3,i4)%d(1,jut)
      cw2=r2_34(i3,i4)%c(2,jut)+rmass(id1)*r2_34(i3,i4)%d(2,jut)
      cz1=p234q*r2_34(i3,i4)%d(1,jut)+rmass(id1)*r2_34(i3,i4)%c(
     & 1,jut)
      cz2=p234q*r2_34(i3,i4)%d(2,jut)+rmass(id1)*r2_34(i3,i4)%c(
     & 2,jut)
      clineth12(i3,i4)%a(2,jut)=clineth12(i3,i4)%a(2,jut)+lh1_23
     & 4%a(2,1)*cx1+lh1_234%c(2,2)*cy2
      clineth12(i3,i4)%b(2,jut)=clineth12(i3,i4)%b(2,jut)+lh1_23
     & 4%b(2,2)*cx2
      clineth12(i3,i4)%c(2,jut)=clineth12(i3,i4)%c(2,jut)+lh1_23
     & 4%a(2,1)*cw1+lh1_234%c(2,2)*cz2
      clineth12(i3,i4)%d(2,jut)=clineth12(i3,i4)%d(2,jut)+lh1_23
     & 4%b(2,2)*cw2
      clineth12(i3,i4)%a(1,jut)=clineth12(i3,i4)%a(1,jut)+lh1_23
     & 4%c(1,1)*cy1+lh1_234%a(1,2)*cx2
      clineth12(i3,i4)%b(1,jut)=clineth12(i3,i4)%b(1,jut)+lh1_23
     & 4%b(1,1)*cx1
      clineth12(i3,i4)%c(1,jut)=clineth12(i3,i4)%c(1,jut)+lh1_23
     & 4%c(1,1)*cz1+lh1_234%a(1,2)*cw2
      clineth12(i3,i4)%d(1,jut)=clineth12(i3,i4)%d(1,jut)+lh1_23
     & 4%b(1,1)*cw1
      end do
      end do
      end do
      endif
  
*                                                                       
*    First all diagrams of the type                                     
*                                                                       
*                   |_Z,f,h__/                                          
*              _h___|        \                                          
*                   |                                                   
  
* rmh < 0 in our convention means no higgs coupling                     
      if (id3.eq.5.and.rmh.ge.0.d0) then
      do i1=1,2
      do i2=1,2
* TsTSC -- aa=clineth34(i1,i2)%a,bb=clineth34(i1,i2)%b,cc=clineth34(i1,i
* 2)%c,dd=clineth34(i1,i2)%d,a1=l3_12(i1,i2)%a,b1=l3_12(i1,i2)%b,c1=l3_1
* 2(i1,i2)%c,d1=l3_12(i1,i2)%d,a2=rh4_312%a,b2=rh4_312%b,c2=rh4_312%c,pr
* q=p312q,m=rmass(id3),nsum=0
      do iut=1,2
      cx1=l3_12(i1,i2)%a(iut,1)+l3_12(i1,i2)%c(iut,1)*rmass(id3)
      cx2=l3_12(i1,i2)%c(iut,2)*p312q+l3_12(i1,i2)%a(iut,2)*rmas
     & s(id3)
      cy1=l3_12(i1,i2)%b(iut,1)+l3_12(i1,i2)%d(iut,1)*rmass(id3)
      cy2=l3_12(i1,i2)%d(iut,2)*p312q+l3_12(i1,i2)%b(iut,2)*rmas
     & s(id3)
      cw1=l3_12(i1,i2)%c(iut,1)*p312q+l3_12(i1,i2)%a(iut,1)*rmas
     & s(id3)
      cw2=l3_12(i1,i2)%a(iut,2)+l3_12(i1,i2)%c(iut,2)*rmass(id3)
      cz1=l3_12(i1,i2)%d(iut,1)*p312q+l3_12(i1,i2)%b(iut,1)*rmas
     & s(id3)
      cz2=l3_12(i1,i2)%b(iut,2)+l3_12(i1,i2)%d(iut,2)*rmass(id3)
      clineth34(i1,i2)%a(iut,2)=cx1*rh4_312%a(1,2)+cx2*rh4_312%b
     & (2,2)
      clineth34(i1,i2)%b(iut,2)=cy1*rh4_312%a(1,2)+cy2*rh4_312%b
     & (2,2)
      clineth34(i1,i2)%c(iut,2)=cw2*rh4_312%c(2,2)
      clineth34(i1,i2)%d(iut,2)=cz2*rh4_312%c(2,2)
      clineth34(i1,i2)%a(iut,1)=cw1*rh4_312%b(1,1)+cw2*rh4_312%a
     & (2,1)
      clineth34(i1,i2)%b(iut,1)=cz1*rh4_312%b(1,1)+cz2*rh4_312%a
     & (2,1)
      clineth34(i1,i2)%c(iut,1)=cx1*rh4_312%c(1,1)
      clineth34(i1,i2)%d(iut,1)=cy1*rh4_312%c(1,1)
      end do
      end do
      end do
  
  
*    Then those of the type                                             
*                                                                       
*              _h___|                                                   
*                   |_Z,f,h__/                                          
*                   |        \                                          
*                                                                       
  
      do i1=1,2
      do i2=1,2
* TSCTs -- aa=clineth34(i1,i2)%a,bb=clineth34(i1,i2)%b,cc=clineth34(i1,i
* 2)%c,dd=clineth34(i1,i2)%d,a1=lh3_412%a,b1=lh3_412%b,c1=lh3_412%c,a2=r
* 4_12(i1,i2)%a,b2=r4_12(i1,i2)%b,c2=r4_12(i1,i2)%c,d2=r4_12(i1,i2)%d,pr
* q=p412q,m=rmass(id3),nsum=1
      do jut=1,2
      cx1=r4_12(i1,i2)%a(1,jut)+rmass(id3)*r4_12(i1,i2)%b(1,jut)
      cx2=r4_12(i1,i2)%a(2,jut)+rmass(id3)*r4_12(i1,i2)%b(2,jut)
      cy1=p412q*r4_12(i1,i2)%b(1,jut)+rmass(id3)*r4_12(i1,i2)%a(
     & 1,jut)
      cy2=p412q*r4_12(i1,i2)%b(2,jut)+rmass(id3)*r4_12(i1,i2)%a(
     & 2,jut)
      cw1=r4_12(i1,i2)%c(1,jut)+rmass(id3)*r4_12(i1,i2)%d(1,jut)
      cw2=r4_12(i1,i2)%c(2,jut)+rmass(id3)*r4_12(i1,i2)%d(2,jut)
      cz1=p412q*r4_12(i1,i2)%d(1,jut)+rmass(id3)*r4_12(i1,i2)%c(
     & 1,jut)
      cz2=p412q*r4_12(i1,i2)%d(2,jut)+rmass(id3)*r4_12(i1,i2)%c(
     & 2,jut)
      clineth34(i1,i2)%a(2,jut)=clineth34(i1,i2)%a(2,jut)+lh3_41
     & 2%a(2,1)*cx1+lh3_412%c(2,2)*cy2
      clineth34(i1,i2)%b(2,jut)=clineth34(i1,i2)%b(2,jut)+lh3_41
     & 2%b(2,2)*cx2
      clineth34(i1,i2)%c(2,jut)=clineth34(i1,i2)%c(2,jut)+lh3_41
     & 2%a(2,1)*cw1+lh3_412%c(2,2)*cz2
      clineth34(i1,i2)%d(2,jut)=clineth34(i1,i2)%d(2,jut)+lh3_41
     & 2%b(2,2)*cw2
      clineth34(i1,i2)%a(1,jut)=clineth34(i1,i2)%a(1,jut)+lh3_41
     & 2%c(1,1)*cy1+lh3_412%a(1,2)*cx2
      clineth34(i1,i2)%b(1,jut)=clineth34(i1,i2)%b(1,jut)+lh3_41
     & 2%b(1,1)*cx1
      clineth34(i1,i2)%c(1,jut)=clineth34(i1,i2)%c(1,jut)+lh3_41
     & 2%c(1,1)*cz1+lh3_412%a(1,2)*cw2
      clineth34(i1,i2)%d(1,jut)=clineth34(i1,i2)%d(1,jut)+lh3_41
     & 2%b(1,1)*cw1
      end do
      end do
      end do
      endif
  
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              ch1234(i1,i2,i3,i4)=czero
              chh1234(i1,i2,i3,i4)=czero
            enddo
          enddo
        enddo
      enddo
  
  
* rmh < 0 in our convention means no higgs coupling                     
      if (id1.eq.5.and.rmh.ge.0.d0)then
      rmassl=rmass(id1)
      rmassr=-rmass(id2)
      do i3=1,2
      do i4=1,2
* mline -- res=ch1234(&1,&2,i3,i4),abcd=clineth12(i3,i4)%,m1=rmassl,m2=r
* massr,den=0,nsum=0
      do iut=1,2
      do jut=1,2
      ch1234(iut,jut,i3,i4)=clineth12(i3,i4)%a(iut,jut)+rmassl*c
     & lineth12(i3,i4)%b(iut,jut)+rmassr*clineth12(i3,i4)%c(iut,
     & jut)+rmassl*rmassr*clineth12(i3,i4)%d(iut,jut)
      enddo
      enddo
      end do
      end do
      endif
  
* rmh < 0 in our convention means no higgs coupling                     
      if (id3.eq.5.and.rmh.ge.0.d0)then
      rmassl=rmass(id3)
      rmassr=-rmass(id4)
      do i1=1,2
      do i2=1,2
* mline -- res=ch1234(i1,i2,&1,&2),abcd=clineth34(i1,i2)%,m1=rmassl,m2=r
* massr,den=0,nsum=1
      do iut=1,2
      do jut=1,2
      ch1234(i1,i2,iut,jut)=ch1234(i1,i2,iut,jut)+clineth34(i1,i
     & 2)%a(iut,jut)+rmassl*clineth34(i1,i2)%b(iut,jut)+rmassr*c
     & lineth34(i1,i2)%c(iut,jut)+rmassl*rmassr*clineth34(i1,i2)
     & %d(iut,jut)
      enddo
      enddo
      end do
      end do
      endif
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              chh1234(i1,i2,i3,i4)=ch1234(i1,i2,i3,i4)*rhhbb
              ch1234(i1,i2,i3,i4)=ch1234(i1,i2,i3,i4)*rhbb
            enddo
          enddo
        enddo
      enddo
  
  
* rmh < 0 in our convention means no higgs coupling                     
      if (rmh.ge.0.d0)then
*    Then the triple vertex                                             
*                                                                       
*                    /_                                                 
*                  Z/                                                   
*             _h___/                                                    
*                  \                                                    
*                  Z\_                                                  
*                    \                                                  
*                                                                       
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* p.q -- p.q=cdum(i1,i2,i3,i4),p=cz12(i1,i2)%e,q=cz34(i3,i4)%e,bef=,aft=
      cdum(i1,i2,i3,i4)=(cz12(i1,i2)%e(0)*cz34(i3,i4)%e(0)-cz12(
     & i1,i2)%e(1)*cz34(i3,i4)%e(1)-cz12(i1,i2)%e(2)*cz34(i3,i4)
     & %e(2)-cz12(i1,i2)%e(3)*cz34(i3,i4)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              ch1234(i1,i2,i3,i4)=ch1234(i1,i2,i3,i4)+
     &             rhzz*cdum(i1,i2,i3,i4)
              chh1234(i1,i2,i3,i4)=chh1234(i1,i2,i3,i4)+
     &             rhhzz*cdum(i1,i2,i3,i4)
            enddo
          enddo
        enddo
      enddo
  
  
*     And finally  the triple vertex                                    
*                                                                       
*                    /_                                                 
*                  h/                                                   
*             _h___/                                                    
*                  \                                                    
*                  h\_                                                  
*                    \                                                  
*                                                                       
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              cdum(i1,i2,i3,i4)=ch12(i1,i2)*ch34(i3,i4)
  
              ch1234(i1,i2,i3,i4)=ch1234(i1,i2,i3,i4)
     &             +cdum(i1,i2,i3,i4)*
     &             (r3h*rhbb/(s12-cmh2)*rhbb/(s34-cmh2)+
     &             r2h1hh*rhbb/(s12-cmh2)*rhhbb/(s34-cmhh2)+
     &             r2h1hh*rhhbb/(s12-cmhh2)*rhbb/(s34-cmh2)+
     &             r1h2hh*rhhbb/(s12-cmhh2)*rhhbb/(s34-cmhh2))
  
              chh1234(i1,i2,i3,i4)=chh1234(i1,i2,i3,i4)
     &             +cdum(i1,i2,i3,i4)*
     &             (r3hh*rhhbb/(s12-cmhh2)*rhhbb/(s34-cmhh2)+
     &             r1h2hh*rhhbb/(s12-cmhh2)*rhbb/(s34-cmh2)+
     &             r1h2hh*rhbb/(s12-cmh2)*rhhbb/(s34-cmhh2)+
     &             r2h1hh*rhbb/(s12-cmh2)*rhbb/(s34-cmh2))
  
  
            enddo
          enddo
        enddo
      enddo
  
      endif
  
*  compute  diagrams  with 3 forks attached to a Wline                  
*               __                                                      
*                 |_Z,f_/                                               
*                 |     \                                               
*                 |                                                     
*                 |__W__/                                               
*                 |     \                                               
*                 |                                                     
*                 |_Z,f_/                                               
*               __|     \                                               
*                                                                       
*                                                                       
*     l5_1234(2,2,2,2),l5_1278(2,2),l5_3478(2,2)                        
*     ,l7_1234(2,2,2,2),l7_1256(2,2),l7_3456(2,2)                       
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              cresw_3fork(i1,i2,i3,i4)=czero
            enddo
          enddo
        enddo
      enddo
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLTR0_W -- aa=cresw_3fork(i1,i2,i3,i4),a1=l5_1278(i1,i2)%a,c1=l5_1278(
* i1,i2)%c,a2=r6_34(i3,i4)%a,b2=r6_34(i3,i4)%b,prq=p634q,bef=cresw_3fork
* (i1,i2,i3,i4)+,aft=
      cresw_3fork(i1,i2,i3,i4)=cresw_3fork(i1,i2,i3,i4)+(l5_1278
     & (i1,i2)%c(2)*p634q*r6_34(i3,i4)%b(1)+l5_1278(i1,i2)%a(2)*
     & r6_34(i3,i4)%a(2))
      end do
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLTR0_W -- aa=cresw_3fork(i1,i2,i3,i4),a1=l5_3478(i3,i4)%a,c1=l5_3478(
* i3,i4)%c,a2=r6_12(i1,i2)%a,b2=r6_12(i1,i2)%b,prq=p612q,bef=cresw_3fork
* (i1,i2,i3,i4)+,aft=
      cresw_3fork(i1,i2,i3,i4)=cresw_3fork(i1,i2,i3,i4)+(l5_3478
     & (i3,i4)%c(2)*p612q*r6_12(i1,i2)%b(1)+l5_3478(i3,i4)%a(2)*
     & r6_12(i1,i2)%a(2))
      end do
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLTR0_W -- aa=cresw_3fork(i1,i2,i3,i4),a1=l5_1234(i1,i2,i3,i4)%a,c1=l5
* _1234(i1,i2,i3,i4)%c,a2=r6_78%a,b2=r6_78%b,prq=p678q,bef=cresw_3fork(i
* 1,i2,i3,i4)+,aft=
      cresw_3fork(i1,i2,i3,i4)=cresw_3fork(i1,i2,i3,i4)+(l5_1234
     & (i1,i2,i3,i4)%c(2)*p678q*r6_78%b(1)+l5_1234(i1,i2,i3,i4)%
     & a(2)*r6_78%a(2))
      end do
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLTR0_W -- aa=cresw_3fork(i1,i2,i3,i4),a1=l7_1256(i1,i2)%a,c1=l7_1256(
* i1,i2)%c,a2=r8_34(i3,i4)%a,b2=r8_34(i3,i4)%b,prq=p834q,bef=cresw_3fork
* (i1,i2,i3,i4)+,aft=
      cresw_3fork(i1,i2,i3,i4)=cresw_3fork(i1,i2,i3,i4)+(l7_1256
     & (i1,i2)%c(2)*p834q*r8_34(i3,i4)%b(1)+l7_1256(i1,i2)%a(2)*
     & r8_34(i3,i4)%a(2))
      end do
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLTR0_W -- aa=cresw_3fork(i1,i2,i3,i4),a1=l7_3456(i3,i4)%a,c1=l7_3456(
* i3,i4)%c,a2=r8_12(i1,i2)%a,b2=r8_12(i1,i2)%b,prq=p812q,bef=cresw_3fork
* (i1,i2,i3,i4)+,aft=
      cresw_3fork(i1,i2,i3,i4)=cresw_3fork(i1,i2,i3,i4)+(l7_3456
     & (i3,i4)%c(2)*p812q*r8_12(i1,i2)%b(1)+l7_3456(i3,i4)%a(2)*
     & r8_12(i1,i2)%a(2))
      end do
      end do
      end do
      end do
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* TLTR0_W -- aa=cresw_3fork(i1,i2,i3,i4),a1=l7_1234(i1,i2,i3,i4)%a,c1=l7
* _1234(i1,i2,i3,i4)%c,a2=r8_56%a,b2=r8_56%b,prq=p856q,bef=cresw_3fork(i
* 1,i2,i3,i4)+,aft=
      cresw_3fork(i1,i2,i3,i4)=cresw_3fork(i1,i2,i3,i4)+(l7_1234
     & (i1,i2,i3,i4)%c(2)*p856q*r8_56%b(1)+l7_1234(i1,i2,i3,i4)%
     & a(2)*r8_56%a(2))
      end do
      end do
      end do
      end do
  
*  compute  diagrams  with 3 forks attached to a Zline                  
*               __                                                      
*                 |_Z,f,h_/                                             
*                 |       \                                             
*                 |                                                     
*                 |__W__/                                               
*                 |     \                                               
*                 |                                                     
*                 |__W__/                                               
*               __|     \                                               
*                                                                       
*                                                                       
*          l1_3456(2,2),l1_3478(2,2),l1_5678                            
*          ,l3_1256(2,2),l3_1278(2,2),l3_5678                           
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              cresz_3fork(i1,i2,i3,i4)=czero
            enddo
          enddo
        enddo
      enddo
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              clinet12_3fork(i1,i2)%a(i3,i4)=czero
              clinet12_3fork(i1,i2)%b(i3,i4)=czero
              clinet12_3fork(i1,i2)%c(i3,i4)=czero
              clinet12_3fork(i1,i2)%d(i3,i4)=czero
            enddo
          enddo
        enddo
      enddo
      do i3=1,2
      do i4=1,2
* TsTs -- aa=clinet12_3fork(i3,i4)%a,bb=clinet12_3fork(i3,i4)%b,cc=cline
* t12_3fork(i3,i4)%c,dd=clinet12_3fork(i3,i4)%d,a1=l1_5678%a,b1=l1_5678%
* b,c1=l1_5678%c,d1=l1_5678%d,a2=r2_34(i3,i4)%a,b2=r2_34(i3,i4)%b,c2=r2_
* 34(i3,i4)%c,d2=r2_34(i3,i4)%d,prq=p234q,m=rmass(id2),nsum=1
      do iut=1,2
      do jut=1,2
      cx1=r2_34(i3,i4)%a(1,jut)+rmass(id2)*r2_34(i3,i4)%b(1,jut)
      cx2=r2_34(i3,i4)%a(2,jut)+rmass(id2)*r2_34(i3,i4)%b(2,jut)
      cy1=p234q*r2_34(i3,i4)%b(1,jut)+rmass(id2)*r2_34(i3,i4)%a(
     & 1,jut)
      cy2=p234q*r2_34(i3,i4)%b(2,jut)+rmass(id2)*r2_34(i3,i4)%a(
     & 2,jut)
      clinet12_3fork(i3,i4)%a(iut,jut)=clinet12_3fork(i3,i4)%a(i
     & ut,jut)+l1_5678%a(iut,1)*cx1+l1_5678%c(iut,1)*cy1+l1_5678
     & %a(iut,2)*cx2+l1_5678%c(iut,2)*cy2
      clinet12_3fork(i3,i4)%b(iut,jut)=clinet12_3fork(i3,i4)%b(i
     & ut,jut)+l1_5678%b(iut,1)*cx1+l1_5678%d(iut,1)*cy1+l1_5678
     & %b(iut,2)*cx2+l1_5678%d(iut,2)*cy2
      cw1=r2_34(i3,i4)%c(1,jut)+rmass(id2)*r2_34(i3,i4)%d(1,jut)
      cw2=r2_34(i3,i4)%c(2,jut)+rmass(id2)*r2_34(i3,i4)%d(2,jut)
      cz1=p234q*r2_34(i3,i4)%d(1,jut)+rmass(id2)*r2_34(i3,i4)%c(
     & 1,jut)
      cz2=p234q*r2_34(i3,i4)%d(2,jut)+rmass(id2)*r2_34(i3,i4)%c(
     & 2,jut)
      clinet12_3fork(i3,i4)%c(iut,jut)=clinet12_3fork(i3,i4)%c(i
     & ut,jut)+l1_5678%a(iut,1)*cw1+l1_5678%c(iut,1)*cz1+l1_5678
     & %a(iut,2)*cw2+l1_5678%c(iut,2)*cz2
      clinet12_3fork(i3,i4)%d(iut,jut)=clinet12_3fork(i3,i4)%d(i
     & ut,jut)+l1_5678%b(iut,1)*cw1+l1_5678%d(iut,1)*cz1+l1_5678
     & %b(iut,2)*cw2+l1_5678%d(iut,2)*cz2
      end do
      end do
      end do
      end do
      if (iup(id1).eq.1) then
      do i3=1,2
      do i4=1,2
* TsTW -- aa=clinet12_3fork(i3,i4)%a,bb=clinet12_3fork(i3,i4)%b,cc=cline
* t12_3fork(i3,i4)%c,dd=clinet12_3fork(i3,i4)%d,a1=l1_3478(i3,i4)%a,b1=l
* 1_3478(i3,i4)%b,c1=l1_3478(i3,i4)%c,d1=l1_3478(i3,i4)%d,a2=r2_56%a,b2=
* r2_56%b,c2=r2_56%c,d2=r2_56%d,prq=p256q,m=rmass(id1-1),nsum=1
      do iut=1,2
      cw1=l1_3478(i3,i4)%c(iut,1)*p256q+l1_3478(i3,i4)%a(iut,1)*
     & rmass(id1-1)
      cw2=l1_3478(i3,i4)%a(iut,2)+l1_3478(i3,i4)%c(iut,2)*rmass(
     & id1-1)
      cz1=l1_3478(i3,i4)%d(iut,1)*p256q+l1_3478(i3,i4)%b(iut,1)*
     & rmass(id1-1)
      cz2=l1_3478(i3,i4)%b(iut,2)+l1_3478(i3,i4)%d(iut,2)*rmass(
     & id1-1)
      clinet12_3fork(i3,i4)%c(iut,1)=clinet12_3fork(i3,i4)%c(iut
     & ,1)+cw1*r2_56%d(1,1)+cw2*r2_56%c(2,1)
      clinet12_3fork(i3,i4)%d(iut,1)=clinet12_3fork(i3,i4)%d(iut
     & ,1)+cz1*r2_56%d(1,1)+cz2*r2_56%c(2,1)
      clinet12_3fork(i3,i4)%a(iut,2)=clinet12_3fork(i3,i4)%a(iut
     & ,2)+cw1*r2_56%b(1,2)+cw2*r2_56%a(2,2)
      clinet12_3fork(i3,i4)%b(iut,2)=clinet12_3fork(i3,i4)%b(iut
     & ,2)+cz1*r2_56%b(1,2)+cz2*r2_56%a(2,2)
      end do
      end do
      end do
      else
      do i3=1,2
      do i4=1,2
* TsTW -- aa=clinet12_3fork(i3,i4)%a,bb=clinet12_3fork(i3,i4)%b,cc=cline
* t12_3fork(i3,i4)%c,dd=clinet12_3fork(i3,i4)%d,a1=l1_3456(i3,i4)%a,b1=l
* 1_3456(i3,i4)%b,c1=l1_3456(i3,i4)%c,d1=l1_3456(i3,i4)%d,a2=r2_78%a,b2=
* r2_78%b,c2=r2_78%c,d2=r2_78%d,prq=p278q,m=rmass(id1+1),nsum=1
      do iut=1,2
      cw1=l1_3456(i3,i4)%c(iut,1)*p278q+l1_3456(i3,i4)%a(iut,1)*
     & rmass(id1+1)
      cw2=l1_3456(i3,i4)%a(iut,2)+l1_3456(i3,i4)%c(iut,2)*rmass(
     & id1+1)
      cz1=l1_3456(i3,i4)%d(iut,1)*p278q+l1_3456(i3,i4)%b(iut,1)*
     & rmass(id1+1)
      cz2=l1_3456(i3,i4)%b(iut,2)+l1_3456(i3,i4)%d(iut,2)*rmass(
     & id1+1)
      clinet12_3fork(i3,i4)%c(iut,1)=clinet12_3fork(i3,i4)%c(iut
     & ,1)+cw1*r2_78%d(1,1)+cw2*r2_78%c(2,1)
      clinet12_3fork(i3,i4)%d(iut,1)=clinet12_3fork(i3,i4)%d(iut
     & ,1)+cz1*r2_78%d(1,1)+cz2*r2_78%c(2,1)
      clinet12_3fork(i3,i4)%a(iut,2)=clinet12_3fork(i3,i4)%a(iut
     & ,2)+cw1*r2_78%b(1,2)+cw2*r2_78%a(2,2)
      clinet12_3fork(i3,i4)%b(iut,2)=clinet12_3fork(i3,i4)%b(iut
     & ,2)+cz1*r2_78%b(1,2)+cz2*r2_78%a(2,2)
      end do
      end do
      end do
      endif
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              clinet34_3fork(i1,i2)%a(i3,i4)=czero
              clinet34_3fork(i1,i2)%b(i3,i4)=czero
              clinet34_3fork(i1,i2)%c(i3,i4)=czero
              clinet34_3fork(i1,i2)%d(i3,i4)=czero
            enddo
          enddo
        enddo
      enddo
      do i1=1,2
      do i2=1,2
* TsTs -- aa=clinet34_3fork(i1,i2)%a,bb=clinet34_3fork(i1,i2)%b,cc=cline
* t34_3fork(i1,i2)%c,dd=clinet34_3fork(i1,i2)%d,a1=l3_5678%a,b1=l3_5678%
* b,c1=l3_5678%c,d1=l3_5678%d,a2=r4_12(i1,i2)%a,b2=r4_12(i1,i2)%b,c2=r4_
* 12(i1,i2)%c,d2=r4_12(i1,i2)%d,prq=p412q,m=rmass(id4),nsum=1
      do iut=1,2
      do jut=1,2
      cx1=r4_12(i1,i2)%a(1,jut)+rmass(id4)*r4_12(i1,i2)%b(1,jut)
      cx2=r4_12(i1,i2)%a(2,jut)+rmass(id4)*r4_12(i1,i2)%b(2,jut)
      cy1=p412q*r4_12(i1,i2)%b(1,jut)+rmass(id4)*r4_12(i1,i2)%a(
     & 1,jut)
      cy2=p412q*r4_12(i1,i2)%b(2,jut)+rmass(id4)*r4_12(i1,i2)%a(
     & 2,jut)
      clinet34_3fork(i1,i2)%a(iut,jut)=clinet34_3fork(i1,i2)%a(i
     & ut,jut)+l3_5678%a(iut,1)*cx1+l3_5678%c(iut,1)*cy1+l3_5678
     & %a(iut,2)*cx2+l3_5678%c(iut,2)*cy2
      clinet34_3fork(i1,i2)%b(iut,jut)=clinet34_3fork(i1,i2)%b(i
     & ut,jut)+l3_5678%b(iut,1)*cx1+l3_5678%d(iut,1)*cy1+l3_5678
     & %b(iut,2)*cx2+l3_5678%d(iut,2)*cy2
      cw1=r4_12(i1,i2)%c(1,jut)+rmass(id4)*r4_12(i1,i2)%d(1,jut)
      cw2=r4_12(i1,i2)%c(2,jut)+rmass(id4)*r4_12(i1,i2)%d(2,jut)
      cz1=p412q*r4_12(i1,i2)%d(1,jut)+rmass(id4)*r4_12(i1,i2)%c(
     & 1,jut)
      cz2=p412q*r4_12(i1,i2)%d(2,jut)+rmass(id4)*r4_12(i1,i2)%c(
     & 2,jut)
      clinet34_3fork(i1,i2)%c(iut,jut)=clinet34_3fork(i1,i2)%c(i
     & ut,jut)+l3_5678%a(iut,1)*cw1+l3_5678%c(iut,1)*cz1+l3_5678
     & %a(iut,2)*cw2+l3_5678%c(iut,2)*cz2
      clinet34_3fork(i1,i2)%d(iut,jut)=clinet34_3fork(i1,i2)%d(i
     & ut,jut)+l3_5678%b(iut,1)*cw1+l3_5678%d(iut,1)*cz1+l3_5678
     & %b(iut,2)*cw2+l3_5678%d(iut,2)*cz2
      end do
      end do
      end do
      end do
      if (iup(id3).eq.1) then
      do i1=1,2
      do i2=1,2
* TsTW -- aa=clinet34_3fork(i1,i2)%a,bb=clinet34_3fork(i1,i2)%b,cc=cline
* t34_3fork(i1,i2)%c,dd=clinet34_3fork(i1,i2)%d,a1=l3_1278(i1,i2)%a,b1=l
* 3_1278(i1,i2)%b,c1=l3_1278(i1,i2)%c,d1=l3_1278(i1,i2)%d,a2=r4_56%a,b2=
* r4_56%b,c2=r4_56%c,d2=r4_56%d,prq=p456q,m=rmass(id3-1),nsum=1
      do iut=1,2
      cw1=l3_1278(i1,i2)%c(iut,1)*p456q+l3_1278(i1,i2)%a(iut,1)*
     & rmass(id3-1)
      cw2=l3_1278(i1,i2)%a(iut,2)+l3_1278(i1,i2)%c(iut,2)*rmass(
     & id3-1)
      cz1=l3_1278(i1,i2)%d(iut,1)*p456q+l3_1278(i1,i2)%b(iut,1)*
     & rmass(id3-1)
      cz2=l3_1278(i1,i2)%b(iut,2)+l3_1278(i1,i2)%d(iut,2)*rmass(
     & id3-1)
      clinet34_3fork(i1,i2)%c(iut,1)=clinet34_3fork(i1,i2)%c(iut
     & ,1)+cw1*r4_56%d(1,1)+cw2*r4_56%c(2,1)
      clinet34_3fork(i1,i2)%d(iut,1)=clinet34_3fork(i1,i2)%d(iut
     & ,1)+cz1*r4_56%d(1,1)+cz2*r4_56%c(2,1)
      clinet34_3fork(i1,i2)%a(iut,2)=clinet34_3fork(i1,i2)%a(iut
     & ,2)+cw1*r4_56%b(1,2)+cw2*r4_56%a(2,2)
      clinet34_3fork(i1,i2)%b(iut,2)=clinet34_3fork(i1,i2)%b(iut
     & ,2)+cz1*r4_56%b(1,2)+cz2*r4_56%a(2,2)
      end do
      end do
      end do
      else
      do i1=1,2
      do i2=1,2
* TsTW -- aa=clinet34_3fork(i1,i2)%a,bb=clinet34_3fork(i1,i2)%b,cc=cline
* t34_3fork(i1,i2)%c,dd=clinet34_3fork(i1,i2)%d,a1=l3_1256(i1,i2)%a,b1=l
* 3_1256(i1,i2)%b,c1=l3_1256(i1,i2)%c,d1=l3_1256(i1,i2)%d,a2=r4_78%a,b2=
* r4_78%b,c2=r4_78%c,d2=r4_78%d,prq=p478q,m=rmass(id3+1),nsum=1
      do iut=1,2
      cw1=l3_1256(i1,i2)%c(iut,1)*p478q+l3_1256(i1,i2)%a(iut,1)*
     & rmass(id3+1)
      cw2=l3_1256(i1,i2)%a(iut,2)+l3_1256(i1,i2)%c(iut,2)*rmass(
     & id3+1)
      cz1=l3_1256(i1,i2)%d(iut,1)*p478q+l3_1256(i1,i2)%b(iut,1)*
     & rmass(id3+1)
      cz2=l3_1256(i1,i2)%b(iut,2)+l3_1256(i1,i2)%d(iut,2)*rmass(
     & id3+1)
      clinet34_3fork(i1,i2)%c(iut,1)=clinet34_3fork(i1,i2)%c(iut
     & ,1)+cw1*r4_78%d(1,1)+cw2*r4_78%c(2,1)
      clinet34_3fork(i1,i2)%d(iut,1)=clinet34_3fork(i1,i2)%d(iut
     & ,1)+cz1*r4_78%d(1,1)+cz2*r4_78%c(2,1)
      clinet34_3fork(i1,i2)%a(iut,2)=clinet34_3fork(i1,i2)%a(iut
     & ,2)+cw1*r4_78%b(1,2)+cw2*r4_78%a(2,2)
      clinet34_3fork(i1,i2)%b(iut,2)=clinet34_3fork(i1,i2)%b(iut
     & ,2)+cz1*r4_78%b(1,2)+cz2*r4_78%a(2,2)
      end do
      end do
      end do
      endif
  
      rmassl=rmass(id1)
      rmassr=-rmass(id2)
  
      do i3=1,2
      do i4=1,2
* mline -- res=cresz_3fork(&1,&2,i3,i4),abcd=clinet12_3fork(i3,i4)%,m1=r
* massl,m2=rmassr,den=0,nsum=1
      do iut=1,2
      do jut=1,2
      cresz_3fork(iut,jut,i3,i4)=cresz_3fork(iut,jut,i3,i4)+clin
     & et12_3fork(i3,i4)%a(iut,jut)+rmassl*clinet12_3fork(i3,i4)
     & %b(iut,jut)+rmassr*clinet12_3fork(i3,i4)%c(iut,jut)+rmass
     & l*rmassr*clinet12_3fork(i3,i4)%d(iut,jut)
      enddo
      enddo
      end do
      end do
  
      rmassl=rmass(id3)
      rmassr=-rmass(id4)
  
      do i1=1,2
      do i2=1,2
* mline -- res=cresz_3fork(i1,i2,&1,&2),abcd=clinet34_3fork(i1,i2)%,m1=r
* massl,m2=rmassr,den=0,nsum=1
      do iut=1,2
      do jut=1,2
      cresz_3fork(i1,i2,iut,jut)=cresz_3fork(i1,i2,iut,jut)+clin
     & et34_3fork(i1,i2)%a(iut,jut)+rmassl*clinet34_3fork(i1,i2)
     & %b(iut,jut)+rmassr*clinet34_3fork(i1,i2)%c(iut,jut)+rmass
     & l*rmassr*clinet34_3fork(i1,i2)%d(iut,jut)
      enddo
      enddo
      end do
      end do
  
  
*  compute diagrams with one W propagator connecting one "W" and        
*  one "Z" on both sides                                                
  
*  start computing cp1256dotcw1256(2,2),cp1256dotcw3478(2,2)            
*     cp1278dotcw1278(2,2),cp1278dotcw3456(2,2)                         
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cp1256dotcw1256(i1,i2),p=p1256(#),q=cw1256(i1,i2,#),bef=,af
* t=
      cp1256dotcw1256(i1,i2)=(p1256(0)*cw1256(i1,i2,0)-p1256(1)*
     & cw1256(i1,i2,1)-p1256(2)*cw1256(i1,i2,2)-p1256(3)*cw1256(
     & i1,i2,3))
      end do
      end do
      do i3=1,2
      do i4=1,2
* p.q -- p.q=cp1256dotcw3478(i3,i4),p=p1256(#),q=cw3478(i3,i4,#),bef=,af
* t=
      cp1256dotcw3478(i3,i4)=(p1256(0)*cw3478(i3,i4,0)-p1256(1)*
     & cw3478(i3,i4,1)-p1256(2)*cw3478(i3,i4,2)-p1256(3)*cw3478(
     & i3,i4,3))
      end do
      end do
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cp1278dotcw1278(i1,i2),p=p1278(#),q=cw1278(i1,i2,#),bef=,af
* t=
      cp1278dotcw1278(i1,i2)=(p1278(0)*cw1278(i1,i2,0)-p1278(1)*
     & cw1278(i1,i2,1)-p1278(2)*cw1278(i1,i2,2)-p1278(3)*cw1278(
     & i1,i2,3))
      end do
      end do
      do i3=1,2
      do i4=1,2
* p.q -- p.q=cp1278dotcw3456(i3,i4),p=p1278(#),q=cw3456(i3,i4,#),bef=,af
* t=
      cp1278dotcw3456(i3,i4)=(p1278(0)*cw3456(i3,i4,0)-p1278(1)*
     & cw3456(i3,i4,1)-p1278(2)*cw3456(i3,i4,2)-p1278(3)*cw3456(
     & i3,i4,3))
      end do
      end do
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              creswz_w_wz(i1,i2,i3,i4)=czero
* p.q -- p.q=creswz_w_wz(i1,i2,i3,i4),p=cw1256(i1,i2,#),q=cw3478(i3,i4,#
* ),bef=-,aft=/(p1256q-cmw2)
      creswz_w_wz(i1,i2,i3,i4)=-(cw1256(i1,i2,0)*cw3478(i3,i4,0)
     & -cw1256(i1,i2,1)*cw3478(i3,i4,1)-cw1256(i1,i2,2)*cw3478(i
     & 3,i4,2)-cw1256(i1,i2,3)*cw3478(i3,i4,3))/(p1256q-cmw2)
* p.q -- p.q=creswz_w_wz(i1,i2,i3,i4),p=cw1278(i1,i2,#),q=cw3456(i3,i4,#
* ),bef=creswz_w_wz(i1,i2,i3,i4)-,aft=/(p1278q-cmw2)
      creswz_w_wz(i1,i2,i3,i4)=creswz_w_wz(i1,i2,i3,i4)-(cw1278(
     & i1,i2,0)*cw3456(i3,i4,0)-cw1278(i1,i2,1)*cw3456(i3,i4,1)-
     & cw1278(i1,i2,2)*cw3456(i3,i4,2)-cw1278(i1,i2,3)*cw3456(i3
     & ,i4,3))/(p1278q-cmw2)
  
              creswz_w_wz(i1,i2,i3,i4)=creswz_w_wz(i1,i2,i3,i4)+
     &          (cp1256dotcw1256(i1,i2)*cp1256dotcw3478(i3,i4))
     &           /cmw2/(p1256q-cmw2)
     &          +(cp1278dotcw1278(i1,i2)*cp1278dotcw3456(i3,i4))
     &           /cmw2/(p1278q-cmw2)
            enddo
          enddo
        enddo
      enddo
  
*  compute diagrams with one Z propagator connecting two "W" on one side
*    and two  "Z" on the other                                          
  
* creszz_z_ww(2,2,2,2)                                                  
  
* p.q -- p.q=cp1234dotcz5678,p=p1234,q=cz5678,bef=,aft=
      cp1234dotcz5678=(p1234(0)*cz5678(0)-p1234(1)*cz5678(1)-p12
     & 34(2)*cz5678(2)-p1234(3)*cz5678(3))
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              creszz_z_ww(i1,i2,i3,i4)=czero
* p.q -- p.q=cp1234dotcz1234,p=p1234(#),q=cz1234(i1,i2,i3,i4,#),bef=,aft
* =
      cp1234dotcz1234=(p1234(0)*cz1234(i1,i2,i3,i4,0)-p1234(1)*c
     & z1234(i1,i2,i3,i4,1)-p1234(2)*cz1234(i1,i2,i3,i4,2)-p1234
     & (3)*cz1234(i1,i2,i3,i4,3))
  
* p.q -- p.q=creszz_z_ww(i1,i2,i3,i4),p=cz1234(i1,i2,i3,i4,#),q=cz5678(#
* ),bef=-,aft=/(p1234q-cmz2)
      creszz_z_ww(i1,i2,i3,i4)=-(cz1234(i1,i2,i3,i4,0)*cz5678(0)
     & -cz1234(i1,i2,i3,i4,1)*cz5678(1)-cz1234(i1,i2,i3,i4,2)*cz
     & 5678(2)-cz1234(i1,i2,i3,i4,3)*cz5678(3))/(p1234q-cmz2)
  
              creszz_z_ww(i1,i2,i3,i4)=creszz_z_ww(i1,i2,i3,i4)+
     &          (cp1234dotcz1234*cp1234dotcz5678)
     &           /cmz2/(p1234q-cmz2)
            enddo
          enddo
        enddo
      enddo
  
  
*  compute diagrams with one gamma propagator connecting two "W"  on    
*   one side and two  "Z" on the other                                  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* p.q -- p.q=creszz_f_ww(i1,i2,i3,i4),p=cf1234(i1,i2,i3,i4,#),q=cf5678(#
* ),bef=-,aft=/p1234q
      creszz_f_ww(i1,i2,i3,i4)=-(cf1234(i1,i2,i3,i4,0)*cf5678(0)
     & -cf1234(i1,i2,i3,i4,1)*cf5678(1)-cf1234(i1,i2,i3,i4,2)*cf
     & 5678(2)-cf1234(i1,i2,i3,i4,3)*cf5678(3))/p1234q
      end do
      end do
      end do
      end do
  
* rmh < 0 in our convention means no higgs coupling                     
      if (rmh.ge.0.d0) then
  
*  compute diagrams with one higgs propagator connecting two "W"  on    
*   one side and two  "Z" on the other                                  
  
* sig                                                                   
  
*      do i1=1,2                                                        
*        do i2=1,2                                                      
*          do i3=1,2                                                    
*            do i4=1,2                                                  
*              creszz_h_ww(i1,i2,i3,i4)=(ch1234(i1,i2,i3,i4)*ch5678)    
*     &             /(p1234q-cmh2)                                      
*            enddo                                                      
*          enddo                                                        
*        enddo                                                          
*      enddo                                                            
        if (i_signal.eq.0) then
          do i1=1,2
            do i2=1,2
              do i3=1,2
                do i4=1,2
                  creszz_h_ww(i1,i2,i3,i4)=(ch1234(i1,i2,i3,i4)*ch5678)
     &                 /(p1234q-cmh2)
               creszz_hh_ww(i1,i2,i3,i4)=(chh1234(i1,i2,i3,i4)*chh5678)
     &                 /(p1234q-cmhh2)
                enddo
              enddo
            enddo
          enddo
        elseif (i_signal.eq.1)then
          if (p1(0).lt.0.d0.and.p2(0).lt.0.d0.or.
     &        p3(0).lt.0.d0.and.p4(0).lt.0.d0.or.
     &        p5(0).lt.0.d0.and.p6(0).lt.0.d0.or.
     &        p7(0).lt.0.d0.and.p8(0).lt.0.d0)then
            do i1=1,2
              do i2=1,2
                do i3=1,2
                  do i4=1,2
                    creszz_h_ww(i1,i2,i3,i4)=czero
                    creszz_hh_ww(i1,i2,i3,i4)=czero
                  enddo
                enddo
              enddo
            enddo
          else
            if (p1234q.gt.0.d0) then
              do i1=1,2
                do i2=1,2
                  do i3=1,2
                    do i4=1,2
                  creszz_h_ww(i1,i2,i3,i4)=(ch1234(i1,i2,i3,i4)*ch5678)
     &                 /(p1234q-cmh2)
               creszz_hh_ww(i1,i2,i3,i4)=(chh1234(i1,i2,i3,i4)*chh5678)
     &                 /(p1234q-cmhh2)
                    enddo
                  enddo
                enddo
              enddo
            else
              do i1=1,2
                do i2=1,2
                  do i3=1,2
                    do i4=1,2
                  creszz_h_ww(i1,i2,i3,i4)=czero
                  creszz_hh_ww(i1,i2,i3,i4)=czero
                    enddo
                  enddo
                enddo
              enddo
            endif
          endif
        else
            print*, 'ERROR: check i_signal value'
            stop
        endif
  
  
* sigend                                                                
      endif
  
* compute the five quadruple coupling diagrams :                        
  
* The - sign to the coupling comes from the fact that we are neglecting 
* cim in every propagator and coupling. The diagram with quartic        
* coupling should have cim**9 and the others cim**11, hence we          
* give a - sign to it with respect to the others.                       
* We factorize out from the amplitude the modulus of the electric       
* charge from every coupling, hence instead of g**2  we have            
* coupling g**2/|e|**2=1/s2w                                            
  
* cresquad(2,2,2,2)                                                     
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              cresquad(i1,i2,i3,i4)=czero
  
*  gamma gamma W+ W-                                                    
* p.q -- p.q=c12fdot34f,p=cf12(i1,i2)%e,q=cf34(i3,i4)%e,bef=,aft=
      c12fdot34f=(cf12(i1,i2)%e(0)*cf34(i3,i4)%e(0)-cf12(i1,i2)%
     & e(1)*cf34(i3,i4)%e(1)-cf12(i1,i2)%e(2)*cf34(i3,i4)%e(2)-c
     & f12(i1,i2)%e(3)*cf34(i3,i4)%e(3))
* p.q -- p.q=c12fdot34z,p=cf12(i1,i2)%e,q=cz34(i3,i4)%e,bef=,aft=
      c12fdot34z=(cf12(i1,i2)%e(0)*cz34(i3,i4)%e(0)-cf12(i1,i2)%
     & e(1)*cz34(i3,i4)%e(1)-cf12(i1,i2)%e(2)*cz34(i3,i4)%e(2)-c
     & f12(i1,i2)%e(3)*cz34(i3,i4)%e(3))
* p.q -- p.q=c12zdot34f,p=cz12(i1,i2)%e,q=cf34(i3,i4)%e,bef=,aft=
      c12zdot34f=(cz12(i1,i2)%e(0)*cf34(i3,i4)%e(0)-cz12(i1,i2)%
     & e(1)*cf34(i3,i4)%e(1)-cz12(i1,i2)%e(2)*cf34(i3,i4)%e(2)-c
     & z12(i1,i2)%e(3)*cf34(i3,i4)%e(3))
* p.q -- p.q=c12zdot34z,p=cz12(i1,i2)%e,q=cz34(i3,i4)%e,bef=,aft=
      c12zdot34z=(cz12(i1,i2)%e(0)*cz34(i3,i4)%e(0)-cz12(i1,i2)%
     & e(1)*cz34(i3,i4)%e(1)-cz12(i1,i2)%e(2)*cz34(i3,i4)%e(2)-c
     & z12(i1,i2)%e(3)*cz34(i3,i4)%e(3))
* p.q -- p.q=c12fdot56,p=cf12(i1,i2)%e,q=cw56%e,bef=,aft=
      c12fdot56=(cf12(i1,i2)%e(0)*cw56%e(0)-cf12(i1,i2)%e(1)*cw5
     & 6%e(1)-cf12(i1,i2)%e(2)*cw56%e(2)-cf12(i1,i2)%e(3)*cw56%e
     & (3))
* p.q -- p.q=c12zdot56,p=cz12(i1,i2)%e,q=cw56%e,bef=,aft=
      c12zdot56=(cz12(i1,i2)%e(0)*cw56%e(0)-cz12(i1,i2)%e(1)*cw5
     & 6%e(1)-cz12(i1,i2)%e(2)*cw56%e(2)-cz12(i1,i2)%e(3)*cw56%e
     & (3))
* p.q -- p.q=c12fdot78,p=cf12(i1,i2)%e,q=cw78%e,bef=,aft=
      c12fdot78=(cf12(i1,i2)%e(0)*cw78%e(0)-cf12(i1,i2)%e(1)*cw7
     & 8%e(1)-cf12(i1,i2)%e(2)*cw78%e(2)-cf12(i1,i2)%e(3)*cw78%e
     & (3))
* p.q -- p.q=c12zdot78,p=cz12(i1,i2)%e,q=cw78%e,bef=,aft=
      c12zdot78=(cz12(i1,i2)%e(0)*cw78%e(0)-cz12(i1,i2)%e(1)*cw7
     & 8%e(1)-cz12(i1,i2)%e(2)*cw78%e(2)-cz12(i1,i2)%e(3)*cw78%e
     & (3))
* p.q -- p.q=c34fdot56,p=cf34(i3,i4)%e,q=cw56%e,bef=,aft=
      c34fdot56=(cf34(i3,i4)%e(0)*cw56%e(0)-cf34(i3,i4)%e(1)*cw5
     & 6%e(1)-cf34(i3,i4)%e(2)*cw56%e(2)-cf34(i3,i4)%e(3)*cw56%e
     & (3))
* p.q -- p.q=c34zdot56,p=cz34(i3,i4)%e,q=cw56%e,bef=,aft=
      c34zdot56=(cz34(i3,i4)%e(0)*cw56%e(0)-cz34(i3,i4)%e(1)*cw5
     & 6%e(1)-cz34(i3,i4)%e(2)*cw56%e(2)-cz34(i3,i4)%e(3)*cw56%e
     & (3))
* p.q -- p.q=c34fdot78,p=cf34(i3,i4)%e,q=cw78%e,bef=,aft=
      c34fdot78=(cf34(i3,i4)%e(0)*cw78%e(0)-cf34(i3,i4)%e(1)*cw7
     & 8%e(1)-cf34(i3,i4)%e(2)*cw78%e(2)-cf34(i3,i4)%e(3)*cw78%e
     & (3))
* p.q -- p.q=c34zdot78,p=cz34(i3,i4)%e,q=cw78%e,bef=,aft=
      c34zdot78=(cz34(i3,i4)%e(0)*cw78%e(0)-cz34(i3,i4)%e(1)*cw7
     & 8%e(1)-cz34(i3,i4)%e(2)*cw78%e(2)-cz34(i3,i4)%e(3)*cw78%e
     & (3))
* p.q -- p.q=c56dot78,p=cw56%e,q=cw78%e,bef=,aft=
      c56dot78=(cw56%e(0)*cw78%e(0)-cw56%e(1)*cw78%e(1)-cw56%e(2
     & )*cw78%e(2)-cw56%e(3)*cw78%e(3))
  
*  coupling 1 after factorizing e**2                                    
              cresquad(i1,i2,i3,i4)=cresquad(i1,i2,i3,i4)+
     &             (2.d0*c12fdot34f*c56dot78-c12fdot56*c34fdot78-
     &             c12fdot78*c34fdot56)
  
*    Z   gamma W+ W-                                                    
  
*  coupling rcotw after factorizing e**2                                
              cresquad(i1,i2,i3,i4)=cresquad(i1,i2,i3,i4)+
     &             (2.d0*c12zdot34f*c56dot78-c12zdot56*c34fdot78-
     &             c12zdot78*c34fdot56)*rcotw
  
*  gamma   Z   W+ W-                                                    
  
*  coupling rcotw after factorizing e**2                                
              cresquad(i1,i2,i3,i4)=cresquad(i1,i2,i3,i4)+
     &             (2.d0*c12fdot34z*c56dot78-c12fdot56*c34zdot78-
     &             c12fdot78*c34zdot56)*rcotw
  
*    Z     Z   W+ W-                                                    
  
*  coupling rcot2w after factorizing e**2                               
              cresquad(i1,i2,i3,i4)=cresquad(i1,i2,i3,i4)+
     &             (2.d0*c12zdot34z*c56dot78-c12zdot56*c34zdot78-
     &             c12zdot78*c34zdot56)*rcot2w
  
*    h     h   W+ W-                                                    
*  coupling -1/(2*s2w) after factorizing e**2                           
* rmh < 0 in our convention means no higgs coupling                     
              if (id1.eq.5.and.id3.eq.5.and.rmh.ge.0.d0) then
                cresquad(i1,i2,i3,i4)=cresquad(i1,i2,i3,i4)
     &               +c56dot78*ch12(i1,i2)*ch34(i3,i4)
     &  *(r2h2w*rhbb**2/(s12-cmh2)/(s34-cmh2)+
     &  r2hh2w*rhhbb**2/(s12-cmhh2)/(s34-cmhh2)+
     &  r1h1hh2w*rhbb*rhhbb/(s12-cmhh2)/(s34-cmh2)+
     &  r1h1hh2w*rhbb*rhhbb/(s12-cmh2)/(s34-cmhh2))
  
              endif
            enddo
          enddo
        enddo
      enddo
  
* sum the results                                                       
  
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
  
              cres(i1,i2,i3,i4)=cresw_3fork(i1,i2,i3,i4)
     &             +cresz_3fork(i1,i2,i3,i4)+creswz_w_wz(i1,i2,i3,i4)
     &             +creszz_z_ww(i1,i2,i3,i4)+creszz_f_ww(i1,i2,i3,i4)
     &             +cresquad(i1,i2,i3,i4)
* rmh < 0 in our convention means no higgs coupling                     
              if (rmh.ge.0.d0) then
* sig                                                                   
                if (i_signal.eq.0) then
                  cres(i1,i2,i3,i4)=cres(i1,i2,i3,i4)
     &                 +creszz_h_ww(i1,i2,i3,i4)
     &                 +creszz_hh_ww(i1,i2,i3,i4)
                elseif (i_signal.eq.1) then
                  cres(i1,i2,i3,i4)=creszz_h_ww(i1,i2,i3,i4)
     &                 +creszz_hh_ww(i1,i2,i3,i4)
                endif
              endif
* sigend                                                                
            enddo
          enddo
        enddo
      enddo
  
  
*                                                                       
* Unitarization                                                         
*                                                                       
  
       if (i_unitarize.eq.1) then
  
        call unit_2z2w(p12,s12,cz12,p34,s34,cz34,
     &    p56,s56,cw56,p78,s78,cw78,cunit)
  
  
        do i1=1,2
        do i2=1,2
        do i3=1,2
        do i4=1,2
         cres(i1,i2,i3,i4)=cres(i1,i2,i3,i4)+cunit(i1,i2,i3,i4)
        enddo
        enddo
        enddo
        enddo
  
       endif
  
  
  
* the resulting amplitudes are now divided   by sqrt(abs(p.k0))         
* for every external particle                                           
  
      sqpk0=sqrt(p1k0*p2k0*p3k0*p4k0*p5k0*p6k0*p7k0*p8k0)
  
      do i1=1,2
        do i2=1,2
          do i3=1,2
            do i4=1,2
              cres(i1,i2,i3,i4)=cres(i1,i2,i3,i4)/sqpk0
            enddo
          enddo
        enddo
      enddo
  
* Sum and average over spin and colour factors                          
* has to be performed externally, when the |amp|**2 is computed         
  
  
      return
      end
  
