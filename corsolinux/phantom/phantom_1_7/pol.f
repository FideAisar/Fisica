********************************************************************
* This routine computes the complex polarization vectors ceps(0:3)
* longitudial, right and left for a W+ or a W- of four momentum p.
*zw
* old * The flag ip =  0 corresponds to W longitudinal
* old *          ip =  1 corresponds to W right
* old *          ip = -1 corresponds to W left 

* The flag ip =  1 corresponds to W longitudinal
*          ip =  3 corresponds to W right
*          ip =  2 corresponds to W left 
******************************************************************** 

      subroutine pol(p,ip,ceps)
  
      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)
  
      dimension p(0:3)
      dimension ceps(0:3)

      DATA czero/(0.d0,0.d0)/, cuno/(1.d0,0.d0)/, cim/(0.d0,1.d0)/      

      rpt2=p(1)**2+p(2)**2
      rp2=rpt2+p(3)**2
      rpt=sqrt(rpt2)
      rp=sqrt(rp2)
*zw
*      if (ip.eq.0) then
      if (ip.eq.1) then
*zwend
        rm2=p(0)**2-rp2
        if (rm2.lt.0.d0) then
          print*,'ERROR in subroutine pol: rm2.lt.0'
        endif
        rm=sqrt(rm2)
        fact=p(0)/rm/rp
        ceps(0)=rp/rm
        ceps(1)=fact*p(1)
        ceps(2)=fact*p(2)
        ceps(3)=fact*p(3)
      else
        rcosthe=p(3)/rp
        sinthe=sqrt(1.d0-rcosthe**2)
        rcosfi=p(1)/rpt
        sinfi=sqrt(1.d0-rcosfi**2)
        sinfi=sign(sinfi,p(2))
        
        sq2=sqrt(2.d0)

        ceps(0)=czero

*zw
*        if (ip.eq.1) then
        if (ip.eq.3) then
*zwend
          ceps(1)=(-rcosthe*rcosfi+cim*sinfi)/sq2
          ceps(2)=(-rcosthe*sinfi-cim*rcosfi)/sq2
          ceps(3)=(+sinthe)/sq2
*zw
*        elseif (ip.eq.-1) then
        elseif (ip.eq.2) then
*zwend
          ceps(1)=(+rcosthe*rcosfi+cim*sinfi)/sq2
          ceps(2)=(+rcosthe*sinfi-cim*rcosfi)/sq2
          ceps(3)=(-sinthe)/sq2
        else
          print *, 'ERROR ip not valid'
          stop
        endif
      endif

      return
      end
