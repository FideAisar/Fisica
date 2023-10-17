
      subroutine unit(skin,cdelta_aij)
      
      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT COMPLEX*16 (c)
      
      include 'common.h'      
      include 'common_unitarization.h'
      
      dimension cdelta_aij(0:2,0:3)
      dimension fij(0:2,0:3),aij0(0:2,0:3),gij(0:2,0:3),
     &  ssj(0:3),ppj(0:1),ddj(0:1),resom(0:2,0:3),rgam(0:2,0:3)
      
      skin2=skin*skin
      g2=elcharge2/s2w
      
       do i=0,2
       do j=0,3
	aij0(i,j)=0d0       
        fij(i,j)=0d0
	gij(i,j)=0d0
	cdelta_aij(i,j)=czero
       enddo
       enddo      

* LO part
       aij0(0,0)=2d0*skin/vev2
       aij0(1,1)=skin/vev2/3d0
       aij0(2,0)=-skin/vev2
       
* NLO part  ! We can use the perturbative relation Imt_NLO = |t|^2 -->
*   Ima_NLO=1/32pi*|a|^2 to compute the imaginary part
*  below there is only the real part
       if (inlo.eq.1) then
       
       rmu2=rmu*rmu
       
       fij(0,0)=fij(0,0)+ (skin2/vev4)*(8d0/3d0*(7d0*alpha4+11d0*alpha5)
     &    +1d0/(16d0*pi*pi)*(25d0/9d0*log(rmu2/skin) + 11d0/54d0))
     
       fij(0,2)=fij(0,2)+ (skin2/vev4)*(8d0/15d0*(2d0*alpha4+alpha5)
     &    +1d0/(16d0*pi*pi)*(1d0/9d0*log(rmu2/skin) - 7d0/135d0))
       
       fij(1,1)=fij(1,1)+ (skin2/vev4)*(4d0/3d0*(alpha4-2d0*alpha5)
     &    +1d0/(16d0*pi*pi)*(-1d0/54d0))
     
       fij(1,3)=fij(1,3)+ (skin2/vev4)*(
     &    +1d0/(16d0*pi*pi)*(7d0/1080d0))

       fij(2,0)=fij(2,0)+ (skin2/vev4)*(16d0/3d0*(2d0*alpha4+alpha5)
     &    +1d0/(16d0*pi*pi)*(10d0/9d0*log(rmu2/skin) + 25d0/108d0))
       
       fij(2,2)=fij(2,2)+ (skin2/vev4)*(4d0/15d0*(alpha4+2d0*alpha5)
     &    +1d0/(16d0*pi*pi)*(2d0/45d0*log(rmu2/skin) - 247d0/5400d0))
       
       endif       

* resonant part

      if (ireson.eq.1) then

      if (isigma.eq.1) then
      
       rg_sigma2=rg_sigma*rg_sigma
       rm_sigma2=rm_sigma*rm_sigma
       rm_sigma4=rm_sigma2*rm_sigma2
       
       resom(0,0)=rm_sigma
       rgam(0,0)=gam_sigma
      
       ssj(0)=rm_sigma2-skin/2d0+rm_sigma4/skin*
     &   log(rm_sigma2/(skin+rm_sigma2))
       ssj(1)=2d0*rm_sigma4/skin+skin/6d0+rm_sigma4/skin2*
     &  (2d0*rm_sigma2+skin)*log(rm_sigma2/(skin+rm_sigma2))
       ssj(2)=rm_sigma4/skin2*(6d0*rm_sigma2+3d0*skin)+
     &   rm_sigma4/(skin2*skin)*(6d0*rm_sigma4+6d0*rm_sigma2*skin+skin2)
     &   *log(rm_sigma2/(skin+rm_sigma2))
       ssj(3)=rm_sigma4/(3d0*skin2*skin)*
     &   (60d0*rm_sigma4+60d0*rm_sigma2*skin+11d0*skin2)+
     &   rm_sigma4/(skin2*skin2)*(2d0*rm_sigma2+skin)*
     &  (10d0*rm_sigma4+10d0*rm_sigma2*skin+skin2)*
     &  log(rm_sigma2/(skin+rm_sigma2))
              
*       fij(0,0)=fij(0,0)-2d0*g2/vev2*ssj(0)     
       fij(0,0)=fij(0,0)-2d0*rg_sigma2/vev2*ssj(0) 
         !suppose there is a misprinting at p.16
       fij(0,2)=fij(0,2)-2d0*rg_sigma2/vev2*ssj(2)
       fij(1,1)=fij(1,1)-2d0*rg_sigma2/vev2*ssj(1)
       fij(1,3)=fij(1,3)-2d0*rg_sigma2/vev2*ssj(3)
       fij(2,0)=fij(2,0)-2d0*rg_sigma2/vev2*ssj(0)
       fij(2,2)=fij(2,2)-2d0*rg_sigma2/vev2*ssj(2)
       
       gij(0,0)=gij(0,0)-3d0*rg_sigma2/vev2*skin2
       
*       print*, 'f20=',fij(2,0)
*       print*, 'f22=',fij(2,2)
*       print*, 'gij=',gij
      endif
      
      if (irho.eq.1) then
      
       rg_rho2=rg_rho*rg_rho
       rm_rho2=rm_rho*rm_rho
       rm_rho4=rm_rho2*rm_rho2      

       resom(1,1)=rm_rho
       rgam(1,1)=gam_rho
                    
       ssj(0)=rm_rho2-skin/2d0+rm_rho4/skin*
     &   log(rm_rho2/(skin+rm_rho2))
       ssj(1)=2d0*rm_rho4/skin+skin/6d0+rm_rho4/skin2*
     &  (2d0*rm_rho2+skin)*log(rm_rho2/(skin+rm_rho2))
       ssj(2)=rm_rho4/skin2*(6d0*rm_rho2+3d0*skin)+
     &   rm_rho4/(skin2*skin)*(6d0*rm_rho4+6d0*rm_rho2*skin+skin2)
     &   *log(rm_rho2/(skin+rm_rho2))
       ssj(3)=rm_rho4/(3d0*skin2*skin)*
     &   (60d0*rm_rho4+60d0*rm_rho2*skin+11d0*skin2)+
     &   rm_rho4/(skin2*skin2)*(2d0*rm_rho2+skin)*
     &  (10d0*rm_rho4+10d0*rm_rho2*skin+skin2)*
     &  log(rm_rho2/(skin+rm_rho2))
     
       ppj(0)=1d0+(2d0*skin+rm_rho2)/skin*log(rm_rho2/(skin+rm_rho2))
       ppj(1)=(rm_rho2+2d0*skin)/skin2*(2d0*skin+(2d0*rm_rho2+skin)*
     &  log(rm_rho2/(skin+rm_rho2)))
     
       fij(0,0)=fij(0,0)-4d0*rg_rho2*ppj(0)-3d0*rg_rho2*skin/rm_rho2
       fij(0,2)=fij(0,2)-4d0*rg_rho2*(2d0*skin+rm_rho2)/rm_rho4*ssj(2)
       fij(1,1)=fij(1,1)-rg_rho2*skin/rm_rho2-2d0*rg_rho2*ppj(1)
       fij(1,3)=fij(1,3)-2d0*rg_rho2*(2d0*skin+rm_rho2)/rm_rho4*ssj(3)
       fij(2,0)=fij(2,0)+2d0*rg_rho2*ppj(0)+3d0*rg_rho2*skin/rm_rho2
       fij(2,2)=fij(2,2)+2d0*rg_rho2*(2d0*skin+rm_rho2)/rm_rho4*ssj(2)       
       
       gij(1,1)=gij(1,1)-2d0/3d0*rg_rho2*skin
       
      endif

      if (iphi.eq.1) then
      
       rg_phi2=rg_phi*rg_phi
       rm_phi2=rm_phi*rm_phi
       rm_phi4=rm_phi2*rm_phi2      

       resom(2,0)=rm_phi
       rgam(2,0)=gam_phi
                    
       ssj(0)=rm_phi2-skin/2d0+rm_phi4/skin*
     &   log(rm_phi2/(skin+rm_phi2))
       ssj(1)=2d0*rm_phi4/skin+skin/6d0+rm_phi4/skin2*
     &  (2d0*rm_phi2+skin)*log(rm_phi2/(skin+rm_phi2))
       ssj(2)=rm_phi4/skin2*(6d0*rm_phi2+3d0*skin)+
     &   rm_phi4/(skin2*skin)*(6d0*rm_phi4+6d0*rm_phi2*skin+skin2)
     &   *log(rm_phi2/(skin+rm_phi2))
       ssj(3)=rm_phi4/(3d0*skin2*skin)*
     &   (60d0*rm_phi4+60d0*rm_phi2*skin+11d0*skin2)+
     &   rm_phi4/(skin2*skin2)*(2d0*rm_phi2+skin)*
     &  (10d0*rm_phi4+10d0*rm_phi2*skin+skin2)*
     &  log(rm_phi2/(skin+rm_phi2))
              
       fij(0,0)=fij(0,0)-g2f       
       
       fij(0,2)=fij(0,2)-5d0/3d0*rg_phi2/vev2*ssj(2)
       fij(1,1)=fij(1,1)+5d0/6d0*rg_phi2/vev2*ssj(1)
       fij(1,3)=fij(1,3)+5d0/6d0*rg_phi2/vev2*ssj(3)
       fij(2,0)=fij(2,0)-1d0/6d0*rg_phi2/vev2*ssj(0)
       fij(2,2)=fij(2,2)-1d0/6d0*rg_phi2/vev2*ssj(2)
       
       gij(2,0)=gij(2,0)-1d0/2d0*rg_phi2/vev2*skin2
       
      endif

      if (iff.eq.1) then
      
       rg_ff2=rg_ff*rg_ff
       rm_ff2=rm_ff*rm_ff
       rm_ff4=rm_ff2*rm_ff2      

       resom(0,2)=rm_ff
       rgam(0,2)=gam_ff
                    
       ssj(0)=rm_ff2-skin/2d0+rm_ff4/skin*
     &   log(rm_ff2/(skin+rm_ff2))
       ssj(1)=2d0*rm_ff4/skin+skin/6d0+rm_ff4/skin2*
     &  (2d0*rm_ff2+skin)*log(rm_ff2/(skin+rm_ff2))
       ssj(2)=rm_ff4/skin2*(6d0*rm_ff2+3d0*skin)+
     &   rm_ff4/(skin2*skin)*(6d0*rm_ff4+6d0*rm_ff2*skin+skin2)
     &   *log(rm_ff2/(skin+rm_ff2))
       ssj(3)=rm_ff4/(3d0*skin2*skin)*
     &   (60d0*rm_ff4+60d0*rm_ff2*skin+11d0*skin2)+
     &   rm_ff4/(skin2*skin2)*(2d0*rm_ff2+skin)*
     &  (10d0*rm_ff4+10d0*rm_ff2*skin+skin2)*
     &  log(rm_ff2/(skin+rm_ff2))
     
       ddj(0)=1d0/2d0*(2d0*rm_ff2+11d0*skin)+
     &   1d0/skin*(rm_ff4+6d0*rm_ff2*skin+6d0*skin2)*
     &   log(rm_ff2/(skin+rm_ff2))
       ddj(1)=1d0/6d0/skin2*
     &   (skin*(12d0*rm_ff4+72d0*rm_ff2*skin+73d0*skin2)+
     & 6d0*(2d0*rm_ff2+skin)*(rm_ff4+6d0*rm_ff2*skin+6d0*skin2)*
     &   log(rm_ff2/(skin+rm_ff2)))
              

       fij(0,0)=fij(0,0)-rg_ff2/3d0/vev2*ddj(0)-
     &   11d0/36d0*rg_ff2/vev2*skin2/rm_ff2
       fij(0,2)=fij(0,2)-
     &   rg_ff2/3d0/vev2*(1d0+6d0*skin/rm_ff2+6d0*skin2/rm_ff4)*ssj(2)-
     &   1d0/180d0*rg_ff2/vev2*skin2/rm_ff2
       fij(1,1)=fij(1,1)-rg_ff2/3d0/vev2*ddj(1)+
     &   1d0/36d0*rg_ff2/vev2*skin2/rm_ff2
       fij(1,3)=fij(1,3)-
     &   rg_ff2/3d0/vev2*(1d0+6d0*skin/rm_ff2+6d0*skin2/rm_ff4)*ssj(3)
       fij(2,0)=fij(2,0)-rg_ff2/3d0/vev2*ddj(0)-
     &   1d0/18d0*rg_ff2/vev2*skin2/rm_ff2
       fij(2,2)=fij(2,2)-
     &  rg_ff2/3d0/vev2*(1d0+6d0*skin/rm_ff2+6d0*skin2/rm_ff4)*ssj(2)-
     &  1d0/180d0*rg_ff2/vev2*skin2/rm_ff2
       
       gij(0,2)=gij(0,2)-rg_ff2/10d0/vev2*skin2
       
      endif            
      
      
      if (ia.eq.1) then
      
       rg_a2=rg_a*rg_a
       rm_a2=rm_a*rm_a
       rm_a4=rm_a2*rm_a2      

       resom(2,2)=rm_a
       rgam(2,2)=gam_a
                    
       ssj(0)=rm_a2-skin/2d0+rm_a4/skin*
     &   log(rm_a2/(skin+rm_a2))
       ssj(1)=2d0*rm_a4/skin+skin/6d0+rm_a4/skin2*
     &  (2d0*rm_a2+skin)*log(rm_a2/(skin+rm_a2))
       ssj(2)=rm_a4/skin2*(6d0*rm_a2+3d0*skin)+
     &   rm_a4/(skin2*skin)*(6d0*rm_a4+6d0*rm_a2*skin+skin2)
     &   *log(rm_a2/(skin+rm_a2))
       ssj(3)=rm_a4/(3d0*skin2*skin)*
     &   (60d0*rm_a4+60d0*rm_a2*skin+11d0*skin2)+
     &   rm_a4/(skin2*skin2)*(2d0*rm_a2+skin)*
     &  (10d0*rm_a4+10d0*rm_a2*skin+skin2)*
     &  log(rm_a2/(skin+rm_a2))
     
       ddj(0)=1d0/2d0*(2d0*rm_a2+11d0*skin)+
     &   1d0/skin*(rm_a4+6d0*rm_a2*skin+6d0*skin2)*
     &   log(rm_a2/(skin+rm_a2))
       ddj(1)=1d0/6d0/skin2*
     &   (skin*(12d0*rm_a4+72d0*rm_a2*skin+73d0*skin2)+
     & 6d0*(2d0*rm_a2+skin)*(rm_a4+6d0*rm_a2*skin+6d0*skin2)*
     &   log(rm_a2/(skin+rm_a2)))
              

       fij(0,0)=fij(0,0)-5d0/6d0*rg_a2/3d0/vev2*ddj(0)-
     &   5d0/108d0*rg_a2/vev2*skin2/rm_a2
       fij(0,2)=fij(0,2)-5d0/6d0*rg_a2/3d0/vev2
     &   *(1d0+6d0*skin/rm_a2+6d0*skin2/rm_a4)*ssj(2)-
     &   1d0/216d0*rg_a2/vev2*skin2/rm_a2
       fij(1,1)=fij(1,1)-5d0/12d0*rg_a2/3d0/vev2*ddj(1)
     &   -5d0/432d0*rg_a2/vev2*skin2/rm_a2
       fij(1,3)=fij(1,3)-5d0/12d0*
     &   rg_a2/3d0/vev2*(1d0+6d0*skin/rm_a2+6d0*skin2/rm_a4)*ssj(3)
       fij(2,0)=fij(2,0)-1d0/12d0*rg_a2/3d0/vev2*ddj(0)-
     &   5d0/108d0*rg_a2/vev2*skin2/rm_a2
       fij(2,2)=fij(2,2)-1d0/12d0*
     &  rg_a2/3d0/vev2*(1d0+6d0*skin/rm_a2+6d0*skin2/rm_a4)*ssj(2)-
     &  1d0/2160d0*rg_a2/vev2*skin2/rm_a2
       
       gij(2,2)=gij(2,2)-rg_a2/60d0/vev2*skin2
       
      endif       
      
      endif
      
* unitarization
C     itype= 1 kmatrix
C            2 ipade
C            3 ind
C            4 ilargen
      
* no unitarization
* PS: if a resonance is requested, it must be introduced its width.
      if (iunittype.eq.0) then       
      
*      print*, 'no unitarization'
*       gam_sigma = gamh  ! a la Higgs
*       gam_sigma = 6.d0*rg_sigma/(64.d0*pi)*rm_sigma2*rm_sigma/vev2

       if (ireson.eq.1) then   
       
        do i=0,2
        do j=0,3
         if (gij(i,j).eq.0d0) then
          cdelta_aij(i,j)=fij(i,j)
	 else
          cresom2=resom(i,j)**2-cim*rgam(i,j)*resom(i,j)
	  cdelta_aij(i,j)=fij(i,j)+gij(i,j)/(skin-cresom2)
         endif
        enddo
        enddo
	
       else       
       
        do i=0,2
        do j=0,3
         cdelta_aij(i,j)=fij(i,j)

        enddo
        enddo
       
       endif
       
*       print*, 'cdelta_aij=',cdelta_aij

      
* K-matrix      
      elseif (iunittype.eq.1) then !kmatrix
      
*      print*, 'k-matrix'
      
       do i=0,2
       do j=0,3

        if (gij(i,j).eq.0d0) then 
         cdelta_aij(i,j)=32d0*pi*cim*(1+cim/(32d0*pi)*aij0(i,j)-
     &    1d0/(1d0-cim/(32d0*pi)*(aij0(i,j)+fij(i,j))))
        else
         cdelta_aij(i,j)=32d0*pi*cim*(1+cim/(32d0*pi)*aij0(i,j)+
     &     (skin-resom(i,j)**2)/(cim*gij(i,j)/(32d0*pi)-
     &     (skin-resom(i,j)**2)*
     &      (1d0-cim/(32d0*pi)*(aij0(i,j)+fij(i,j)))))
	endif

       enddo
       enddo

*	 print*, 'cdelta_a20=',-cdelta_aij(2,0)/elcharge2
*	 print*, 'cdelta_a22=',-cdelta_aij(2,2)/elcharge2
	 

* IAM unitarization.
* it must be included the imaginary part of the amplitude 
* Imt^(4) = |t^(2)|^2 --> ImA^(1)=1/(32pi)|A^(0)|^2
      elseif (iunittype.eq.2) then !pade approximant
      
*      print*, 'IAM '
      
      
        if (inlo.eq.1.and.ireson.eq.0) then
	  
	 cdelta_aij(0,0)=aij0(0,0)*(fij(0,0)+cim/32d0/pi*aij0(0,0)**2)
     &    /(aij0(0,0)-(fij(0,0)+cim/32d0/pi*aij0(0,0)**2))
	 
	 cdelta_aij(1,1)=aij0(1,1)*(fij(1,1)+cim/32d0/pi*aij0(1,1)**2)
     &    /(aij0(1,1)-(fij(1,1)+cim/32d0/pi*aij0(1,1)**2))

*          print*, sqrt(skin),aij0(1,1)-(fij(1,1)+cim*aij0(1,1)**2)
	  
	 cdelta_aij(2,0)=aij0(2,0)*(fij(2,0)+cim/32d0/pi*aij0(2,0)**2)
     &    /(aij0(2,0)-(fij(2,0)+cim/32d0/pi*aij0(2,0)**2))

	 
	else
	 print*, 'error: The Pade unitarization must
     & have only NLO flag turned on.'
         stop
	endif
      
      elseif (iunittype.eq.3) then !N/D
      
        if (inlo.eq.1.and.ireson.eq.0) then

* definition is g=1/pi*log(-s/M^2)
* gbar = 1/pi*log(s/M^2)

       gbn=1d0/32d0/pi/pi*log(skin/rmnd**2)
       
        xij = aij0(0,0)+fij(0,0)+gbn*aij0(0,0)**2
	cdij = 1d0+gbn*xij-cim/32d0/pi*xij
	cdelta_aij(0,0)=xij/cdij-aij0(0,0)

        xij = aij0(1,1)+fij(1,1)+gbn*aij0(1,1)**2
	cdij = 1d0+gbn*xij-cim/32d0/pi*xij
	cdelta_aij(1,1)=xij/cdij-aij0(1,1)

        xij = aij0(2,0)+fij(2,0)+gbn*aij0(2,0)**2
	cdij = 1d0+gbn*xij-cim/32d0/pi*xij
	cdelta_aij(2,0)=xij/cdij-aij0(2,0)
       
       	else
	 print*, 'error: The N/D unitarization must
     & have only NLO flag turned on.'
         stop
	endif



      elseif (iunittype.eq.4) then !Large N limit
       print*, 'Large N unitarization not implemented yet'
       stop     
      else
       print*, 'wrong unitarization type'
       stop
      endif
      
*            print*, 'cdeltaA_00=',cdelta_aij(0,0)
*            print*, 'cdeltaA_02=',cdelta_aij(0,2)
*            print*, 'cdeltaA_11=',cdelta_aij(1,1)
*            print*, 'cdeltaA_20=',cdelta_aij(2,0)
*            print*, 'cdeltaA_22=',cdelta_aij(2,2)

      
            
      
      
      end
