********************************************************************
* This routine computes the complex polarization vectors ceps(0:3)
* longitudial, right and left for a W+ or a W- of four momentum p.
* The flag ip =  0 corresponds to W longitudinal
*          ip =  1 corresponds to W+ right and W- left
*          ip = -1 corresponds to W+ left and W- right
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

      if (ip.eq.0) then
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

        ceps(0)=czero

        if (ip.eq.1) then
          ceps(1)=p(1)*p(3)/rp/rpt
          ceps(2)=p(2)*p(3)/rp/rpt
          ceps(3)=-rpt/rp
        elseif (ip.eq.-1) then
          ceps(1)=-p(2)/rpt
          ceps(2)=p(1)/rpt
          ceps(3)=czero
        else
          print *, 'ERROR ip not valid'
          stop
        endif
      endif

      return
      end
