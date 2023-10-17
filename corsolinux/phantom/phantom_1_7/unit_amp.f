  
* The - sign to the coupling comes from the fact that we are neglecting 
* cim in every propagator and coupling. The diagram with quartic        
* coupling should have cim**9 and the others cim**11, hence we          
* give a - sign to it with respect to the others.                       
* We factorize out from the amplitude the modulus of the electric       
* charge from every coupling, hence we do at the end                    
* cunit = -cunit/elcharge2                                              
  
      subroutine unit_4w(pwp1,swp1,cewp1,pwm1,swm1,cewm1,
     &   pwp2,swp2,cewp2,pwm2,swm2,cewm2,cunit)
  
      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT COMPLEX*16 (c)
  
      include 'common.h'
      include 'common_unitarization.h'
*four momenta                                                           
      DIMENSION pwp1(0:3),pwm1(0:3),pwp2(0:3),pwm2(0:3),paux(0:3)
  
      type polcom
       double complex e(0:3),ek0,v
      end type
      type(polcom) cewp1,cewm1,cewp2,cewm2
  
      dimension cdelta_aij(0:2,0:3)
      
*      print*, '4w'
  
* check if we have ww scattering                                        
      nout=0
      IF(swp1.GT.0.d0)THEN
        nout=nout+1
        iout1=1
      ELSE
        iout1=-1
      ENDIF
      IF(swm1.GT.0.d0)THEN
        nout=nout+1
        iout2=1
      ELSE
        iout2=-1
      ENDIF
      IF(swp2.GT.0.d0)THEN
        nout=nout+1
        iout3=1
      ELSE
        iout3=-1
      ENDIF
      IF(swm2.GT.0.d0)THEN
        nout=nout+1
        iout4=1
      ELSE
        iout4=-1
      ENDIF
  
      IF(nout.NE.2)THEN
        cunit=czero
        RETURN
      ENDIF
  
* Lorentz variables (in terms of w+w+ -> w-w-)                          
  
* p.q -- p.q=cdot1,p=cewp1.e,q=cewp2.e,bef=,aft=
      cdot1=(cewp1%e(0)*cewp2%e(0)-cewp1%e(1)*cewp2%e(1)-cewp1%e
     & (2)*cewp2%e(2)-cewp1%e(3)*cewp2%e(3))
* p.q -- p.q=cdot2,p=cewm1.e,q=cewm2.e,bef=,aft=
      cdot2=(cewm1%e(0)*cewm2%e(0)-cewm1%e(1)*cewm2%e(1)-cewm1%e
     & (2)*cewm2%e(2)-cewm1%e(3)*cewm2%e(3))
  
      cp1p2lor = 4d0*rmw2*rmw2*cdot1*cdot2/vev4
  
* p.q -- p.q=cdot1,p=cewp1.e,q=cewm1.e,bef=,aft=
      cdot1=(cewp1%e(0)*cewm1%e(0)-cewp1%e(1)*cewm1%e(1)-cewp1%e
     & (2)*cewm1%e(2)-cewp1%e(3)*cewm1%e(3))
* p.q -- p.q=cdot2,p=cewp2.e,q=cewm2.e,bef=,aft=
      cdot2=(cewp2%e(0)*cewm2%e(0)-cewp2%e(1)*cewm2%e(1)-cewp2%e
     & (2)*cewm2%e(2)-cewp2%e(3)*cewm2%e(3))
  
      cp1m1lor = 4d0*rmw2*rmw2*cdot1*cdot2/vev4
  
* p.q -- p.q=cdot1,p=cewp1.e,q=cewm2.e,bef=,aft=
      cdot1=(cewp1%e(0)*cewm2%e(0)-cewp1%e(1)*cewm2%e(1)-cewp1%e
     & (2)*cewm2%e(2)-cewp1%e(3)*cewm2%e(3))
* p.q -- p.q=cdot2,p=cewp2.e,q=cewm1.e,bef=,aft=
      cdot2=(cewp2%e(0)*cewm1%e(0)-cewp2%e(1)*cewm1%e(1)-cewp2%e
     & (2)*cewm1%e(2)-cewp2%e(3)*cewm1%e(3))
  
      cp1m2lor = 4d0*rmw2*rmw2*cdot1*cdot2/vev4
  
  
      IF(iout1*iout3.EQ.1)THEN   ! W+W+ or W-W-
      
      
*      print*, 'right place w+w+ <-> w-w-'
      
* kinematical s                                                         
       do m=0,3
        paux(m)=pwp1(m)+pwp2(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then
*      if (skin.lt.10000d0) then      
       cunit=czero
       return      
      endif
        
*      print*, 'skin=', skin                                            
*      print*, 'skin2=', skin2                                          
  
      call unit(skin,cdelta_aij)
  
      cunit= 8d0*(vev4/(8d0*skin2)*
     & (cdelta_aij(2,0)-10d0*cdelta_aij(2,2)))*cp1p2lor+
     & 4d0*(15d0*vev4/(4d0*skin2)*cdelta_aij(2,2))*
     & (cp1m1lor+cp1m2lor)
  
      else if (iout1*iout2.eq.1) then !  w+1w-1 <-> w+2w-2
      
*      print*, 'wrong place w+1w-1 <-> w+2w-2'
* kinematical s                                                         
       do m=0,3
        paux(m)=pwp1(m)+pwm1(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
       cunit=czero
       return      
      endif
        
      call unit(skin,cdelta_aij)

* inversion t,u fixed 2/8/10
      cunit= 4d0*(vev4/(24d0*skin2)*
     & (2d0*cdelta_aij(0,0)+cdelta_aij(2,0))-
     & 5d0*vev4/(12d0*skin2)*(2d0*cdelta_aij(0,2)+cdelta_aij(2,2)))*
     & cp1m1lor +
     & 4d0*(vev4/(8d0*skin2)*
     & (10d0*cdelta_aij(0,2)-3d0*cdelta_aij(1,1)+5d0*cdelta_aij(2,2)))
     & *cp1m2lor +
     & 8d0*(vev4/(16d0*skin2)*
     & (10d0*cdelta_aij(0,2)+3d0*cdelta_aij(1,1)+5d0*cdelta_aij(2,2)))
     & *cp1p2lor
  
      else  !  w+1w-2 <-> w+2w-1
      
*       print*, 'wrong place w+1w-2 <-> w+2w-1'
	    
* kinematical s                                                         
       do m=0,3
        paux(m)=pwp1(m)+pwm2(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin
  
* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
       cunit=czero
       return      
      endif
      
        
      call unit(skin,cdelta_aij)
  
      cunit= 4d0*(vev4/(24d0*skin2)*
     & (2d0*cdelta_aij(0,0)+cdelta_aij(2,0))-
     & 5d0*vev4/(12d0*skin2)*(2d0*cdelta_aij(0,2)+cdelta_aij(2,2)))*
     & cp1m2lor +
     & 4d0*(vev4/(8d0*skin2)*
     & (10d0*cdelta_aij(0,2)-3d0*cdelta_aij(1,1)+5d0*cdelta_aij(2,2)))
     & *cp1m1lor +
     & 8d0*(vev4/(16d0*skin2)*
     & (10d0*cdelta_aij(0,2)+3d0*cdelta_aij(1,1)+5d0*cdelta_aij(2,2)))
     & *cp1p2lor
  
  
      endif
  
       cunit = -cunit/elcharge2
       
       
*       stop
       
      end
  
      subroutine unit_4z(pz1,sz1,cez1,pz2,sz2,cez2,
     &   pz3,sz3,cez3,pz4,sz4,cez4,cunit)
  
      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT COMPLEX*16 (c)
  
      include 'common.h'
      include 'common_unitarization.h'
*four momenta                                                           
      DIMENSION pz1(0:3),pz2(0:3),pz3(0:3),pz4(0:3),paux(0:3)
  
      type pol
         double complex e(0:3),ek0,v
      end type
      type(pol) cez1(2,2),cez2(2,2),cez3(2,2),cez4(2,2)
  
      dimension cdelta_aij(0:2,0:3)
  
      dimension cdot1(2,2,2,2),cdot2(2,2,2,2)
  
      type three_forks_amp
        double complex pol(2,2)
       end type
       type(three_forks_amp) cz1z2lor(2,2,2,2,2,2),
     &  cz1z3lor(2,2,2,2,2,2),cz1z4lor(2,2,2,2,2,2),cunit(2,2,2,2,2,2)
  
* check if we have zz scattering                                        
      nout=0
      IF(sz1.GT.0.d0)THEN
        nout=nout+1
        iout1=1
      ELSE
        iout1=-1
      ENDIF
      IF(sz2.GT.0.d0)THEN
        nout=nout+1
        iout2=1
      ELSE
        iout2=-1
      ENDIF
      IF(sz3.GT.0.d0)THEN
        nout=nout+1
        iout3=1
      ELSE
        iout3=-1
      ENDIF
      IF(sz4.GT.0.d0)THEN
        nout=nout+1
        iout4=1
      ELSE
        iout4=-1
      ENDIF
  
      IF(nout.NE.2)THEN
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
      do i7=1,2
      do i8=1,2
        cunit(i1,i2,i3,i4,i5,i6)%pol(i7,i8)=czero
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
  
  
        RETURN
      ENDIF
  
* Lorentz variables (in terms of z1z2 -> z3z4)                          
  
* s                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* p.q -- p.q=cdot1(i1,i2,i3,i4),p=cez1(i1,i2).e,q=cez2(i3,i4).e,bef=,aft
* =
      cdot1(i1,i2,i3,i4)=(cez1(i1,i2)%e(0)*cez2(i3,i4)%e(0)-cez1
     & (i1,i2)%e(1)*cez2(i3,i4)%e(1)-cez1(i1,i2)%e(2)*cez2(i3,i4
     & )%e(2)-cez1(i1,i2)%e(3)*cez2(i3,i4)%e(3))
      end do
      end do
      end do
      end do
      do i5=1,2
      do i6=1,2
      do i7=1,2
      do i8=1,2
* p.q -- p.q=cdot2(i5,i6,i7,i8),p=cez3(i5,i6).e,q=cez4(i7,i8).e,bef=,aft
* =
      cdot2(i5,i6,i7,i8)=(cez3(i5,i6)%e(0)*cez4(i7,i8)%e(0)-cez3
     & (i5,i6)%e(1)*cez4(i7,i8)%e(1)-cez3(i5,i6)%e(2)*cez4(i7,i8
     & )%e(2)-cez3(i5,i6)%e(3)*cez4(i7,i8)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
      do i7=1,2
      do i8=1,2
       cz1z2lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8) =
     &  4d0*rmz2*rmz2*cdot1(i1,i2,i3,i4)*cdot2(i5,i6,i7,i8)/vev4
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
  
* t                                                                     
      do i1=1,2
      do i2=1,2
      do i5=1,2
      do i6=1,2
* p.q -- p.q=cdot1(i1,i2,i5,i6),p=cez1(i1,i2).e,q=cez3(i5,i6).e,bef=,aft
* =
      cdot1(i1,i2,i5,i6)=(cez1(i1,i2)%e(0)*cez3(i5,i6)%e(0)-cez1
     & (i1,i2)%e(1)*cez3(i5,i6)%e(1)-cez1(i1,i2)%e(2)*cez3(i5,i6
     & )%e(2)-cez1(i1,i2)%e(3)*cez3(i5,i6)%e(3))
      end do
      end do
      end do
      end do
      do i3=1,2
      do i4=1,2
      do i7=1,2
      do i8=1,2
* p.q -- p.q=cdot2(i3,i4,i7,i8),p=cez2(i3,i4).e,q=cez4(i7,i8).e,bef=,aft
* =
      cdot2(i3,i4,i7,i8)=(cez2(i3,i4)%e(0)*cez4(i7,i8)%e(0)-cez2
     & (i3,i4)%e(1)*cez4(i7,i8)%e(1)-cez2(i3,i4)%e(2)*cez4(i7,i8
     & )%e(2)-cez2(i3,i4)%e(3)*cez4(i7,i8)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
      do i7=1,2
      do i8=1,2
        cz1z3lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8) =
     &  4d0*rmz2*rmz2*cdot1(i1,i2,i5,i6)*cdot2(i3,i4,i7,i8)/vev4
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
  
* u                                                                     
      do i1=1,2
      do i2=1,2
      do i7=1,2
      do i8=1,2
* p.q -- p.q=cdot1(i1,i2,i7,i8),p=cez1(i1,i2).e,q=cez4(i7,i8).e,bef=,aft
* =
      cdot1(i1,i2,i7,i8)=(cez1(i1,i2)%e(0)*cez4(i7,i8)%e(0)-cez1
     & (i1,i2)%e(1)*cez4(i7,i8)%e(1)-cez1(i1,i2)%e(2)*cez4(i7,i8
     & )%e(2)-cez1(i1,i2)%e(3)*cez4(i7,i8)%e(3))
      end do
      end do
      end do
      end do
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
* p.q -- p.q=cdot2(i3,i4,i5,i6),p=cez2(i3,i4).e,q=cez3(i5,i6).e,bef=,aft
* =
      cdot2(i3,i4,i5,i6)=(cez2(i3,i4)%e(0)*cez3(i5,i6)%e(0)-cez2
     & (i3,i4)%e(1)*cez3(i5,i6)%e(1)-cez2(i3,i4)%e(2)*cez3(i5,i6
     & )%e(2)-cez2(i3,i4)%e(3)*cez3(i5,i6)%e(3))
      end do
      end do
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
      do i7=1,2
      do i8=1,2
        cz1z4lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8) =
     &  4d0*rmz2*rmz2*cdot1(i1,i2,i7,i8)*cdot2(i3,i4,i5,i6)/vev4
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
  
      if (iout1*iout2.eq.1) then ! z1z2 <->z3z4
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pz2(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
         do i1=1,2
	 do i2=1,2
	 do i3=1,2
	 do i4=1,2
	 do i5=1,2
	 do i6=1,2
	 do i7=1,2
	 do i8=1,2
           cunit(i1,i2,i3,i4,i5,i6)%pol(i7,i8)=czero
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
       return      
      endif
      
        
      call unit(skin,cdelta_aij)
  
      cfact1=8d0*(vev4/(24d0*skin2)*
     & (cdelta_aij(0,0)+2d0*cdelta_aij(2,0))-
     &   5d0*vev4/(12d0*skin2)*
     & (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      cfact2=8d0*(5d0*vev4/(8d0*skin2)*
     &  (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
      do i7=1,2
      do i8=1,2
  
      cunit(i1,i2,i3,i4,i5,i6)%pol(i7,i8)=
     &  cfact1*cz1z2lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8)
     & +cfact2*(cz1z3lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8)
     & +cz1z4lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8))
  
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
  
  
      elseif (iout1*iout3.eq.1) then ! z1z3 <->z2z4
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pz3(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
         do i1=1,2
	 do i2=1,2
	 do i3=1,2
	 do i4=1,2
	 do i5=1,2
	 do i6=1,2
	 do i7=1,2
	 do i8=1,2
           cunit(i1,i2,i3,i4,i5,i6)%pol(i7,i8)=czero
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=8d0*(vev4/(24d0*skin2)*
     & (cdelta_aij(0,0)+2d0*cdelta_aij(2,0))-
     &   5d0*vev4/(12d0*skin2)*
     & (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      cfact2=8d0*(5d0*vev4/(8d0*skin2)*
     &  (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
      do i7=1,2
      do i8=1,2
  
      cunit(i1,i2,i3,i4,i5,i6)%pol(i7,i8)=
     &  cfact1*cz1z3lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8)
     & +cfact2*(cz1z2lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8)
     & +cz1z4lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8))
  
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
  
      elseif (iout1*iout4.eq.1) then ! z1z4 <->z2z3
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pz4(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
         do i1=1,2
	 do i2=1,2
	 do i3=1,2
	 do i4=1,2
	 do i5=1,2
	 do i6=1,2
	 do i7=1,2
	 do i8=1,2
           cunit(i1,i2,i3,i4,i5,i6)%pol(i7,i8)=czero
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
	 enddo
       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=8d0*(vev4/(24d0*skin2)*
     & (cdelta_aij(0,0)+2d0*cdelta_aij(2,0))-
     &   5d0*vev4/(12d0*skin2)*
     & (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      cfact2=8d0*(5d0*vev4/(8d0*skin2)*
     &  (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
      do i7=1,2
      do i8=1,2
  
      cunit(i1,i2,i3,i4,i5,i6)%pol(i7,i8)=
     & cfact1*cz1z4lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8)
     & +cfact2*(cz1z2lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8)+
     & cz1z3lor(i1,i2,i3,i4,i5,i6)%pol(i7,i8))
  
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
  
      endif
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      do i5=1,2
      do i6=1,2
      do i7=1,2
      do i8=1,2
  
      cunit(i1,i2,i3,i4,i5,i6)%pol(i7,i8)=
     &  -cunit(i1,i2,i3,i4,i5,i6)%pol(i7,i8)/elcharge2
  
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
  
  
      end
  
      subroutine unit_4z_massless(pz1,sz1,cez1,pz2,sz2,cez2,
     &   pz3,sz3,cez3,pz4,sz4,cez4,cunit)
  
      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT COMPLEX*16 (c)
  
      include 'common.h'
      include 'common_unitarization.h'
*four momenta                                                           
      DIMENSION pz1(0:3),pz2(0:3),pz3(0:3),pz4(0:3),paux(0:3)
  
      type pol
         double complex e(0:3),ek0
      end type
      type(pol) cez1(2),cez2(2),cez3(2),cez4(2)
  
      dimension cdelta_aij(0:2,0:3)
  
      dimension cdot1(2,2),cdot2(2,2)
      dimension cz1z2lor(2,2,2,2),cz1z3lor(2,2,2,2),cz1z4lor(2,2,2,2),
     &  cunit(2,2,2,2)
  
* check if we have zz scattering                                        
      nout=0
      IF(sz1.GT.0.d0)THEN
        nout=nout+1
        iout1=1
      ELSE
        iout1=-1
      ENDIF
      IF(sz2.GT.0.d0)THEN
        nout=nout+1
        iout2=1
      ELSE
        iout2=-1
      ENDIF
      IF(sz3.GT.0.d0)THEN
        nout=nout+1
        iout3=1
      ELSE
        iout3=-1
      ENDIF
      IF(sz4.GT.0.d0)THEN
        nout=nout+1
        iout4=1
      ELSE
        iout4=-1
      ENDIF
  
      IF(nout.NE.2)THEN
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
        cunit(i1,i3,i5,i7)=czero
      enddo
      enddo
      enddo
      enddo
  
  
        RETURN
      ENDIF
  
* Lorentz variables (in terms of z1z2 -> z3z4)                          
  
* s                                                                     
      do i1=1,2
      do i3=1,2
* p.q -- p.q=cdot1(i1,i3),p=cez1(i1).e,q=cez2(i3).e,bef=,aft=
      cdot1(i1,i3)=(cez1(i1)%e(0)*cez2(i3)%e(0)-cez1(i1)%e(1)*ce
     & z2(i3)%e(1)-cez1(i1)%e(2)*cez2(i3)%e(2)-cez1(i1)%e(3)*cez
     & 2(i3)%e(3))
      end do
      end do
      do i5=1,2
      do i7=1,2
* p.q -- p.q=cdot2(i5,i7),p=cez3(i5).e,q=cez4(i7).e,bef=,aft=
      cdot2(i5,i7)=(cez3(i5)%e(0)*cez4(i7)%e(0)-cez3(i5)%e(1)*ce
     & z4(i7)%e(1)-cez3(i5)%e(2)*cez4(i7)%e(2)-cez3(i5)%e(3)*cez
     & 4(i7)%e(3))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
       cz1z2lor(i1,i3,i5,i7) =
     &  4d0*rmz2*rmz2*cdot1(i1,i3)*cdot2(i5,i7)/vev4
      enddo
      enddo
      enddo
      enddo
  
* t                                                                     
      do i1=1,2
      do i5=1,2
* p.q -- p.q=cdot1(i1,i5),p=cez1(i1)%e,q=cez3(i5)%e,bef=,aft=
      cdot1(i1,i5)=(cez1(i1)%e(0)*cez3(i5)%e(0)-cez1(i1)%e(1)*ce
     & z3(i5)%e(1)-cez1(i1)%e(2)*cez3(i5)%e(2)-cez1(i1)%e(3)*cez
     & 3(i5)%e(3))
      end do
      end do
      do i3=1,2
      do i7=1,2
* p.q -- p.q=cdot2(i3,i7),p=cez2(i3)%e,q=cez4(i7)%e,bef=,aft=
      cdot2(i3,i7)=(cez2(i3)%e(0)*cez4(i7)%e(0)-cez2(i3)%e(1)*ce
     & z4(i7)%e(1)-cez2(i3)%e(2)*cez4(i7)%e(2)-cez2(i3)%e(3)*cez
     & 4(i7)%e(3))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
       cz1z3lor(i1,i3,i5,i7) =
     &  4d0*rmz2*rmz2*cdot1(i1,i5)*cdot2(i3,i7)/vev4
      enddo
      enddo
      enddo
      enddo
  
* u                                                                     
      do i1=1,2
      do i7=1,2
* p.q -- p.q=cdot1(i1,i7),p=cez1(i1)%e,q=cez4(i7)%e,bef=,aft=
      cdot1(i1,i7)=(cez1(i1)%e(0)*cez4(i7)%e(0)-cez1(i1)%e(1)*ce
     & z4(i7)%e(1)-cez1(i1)%e(2)*cez4(i7)%e(2)-cez1(i1)%e(3)*cez
     & 4(i7)%e(3))
      end do
      end do
      do i3=1,2
      do i5=1,2
* p.q -- p.q=cdot2(i3,i5),p=cez2(i3)%e,q=cez3(i5)%e,bef=,aft=
      cdot2(i3,i5)=(cez2(i3)%e(0)*cez3(i5)%e(0)-cez2(i3)%e(1)*ce
     & z3(i5)%e(1)-cez2(i3)%e(2)*cez3(i5)%e(2)-cez2(i3)%e(3)*cez
     & 3(i5)%e(3))
      end do
      end do
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
       cz1z4lor(i1,i3,i5,i7) =
     &  4d0*rmz2*rmz2*cdot1(i1,i7)*cdot2(i3,i5)/vev4
      enddo
      enddo
      enddo
      enddo
  
      if (iout1*iout2.eq.1) then ! z1z2 <->z3z4
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pz2(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
	do i1=1,2
	do i3=1,2
	do i5=1,2
	do i7=1,2
          cunit(i1,i3,i5,i7)=czero
	enddo
	enddo
	enddo
	enddo

       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=8d0*(vev4/(24d0*skin2)*
     & (cdelta_aij(0,0)+2d0*cdelta_aij(2,0))-
     &   5d0*vev4/(12d0*skin2)*
     & (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      cfact2=8d0*(5d0*vev4/(8d0*skin2)*
     &  (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
  
      cunit(i1,i3,i5,i7)=cfact1*cz1z2lor(i1,i3,i5,i7)
     & +cfact2*(cz1z3lor(i1,i3,i5,i7)+cz1z4lor(i1,i3,i5,i7))
  
      enddo
      enddo
      enddo
      enddo
  
      elseif (iout1*iout3.eq.1) then ! z1z3 <->z2z4
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pz3(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
	do i1=1,2
	do i3=1,2
	do i5=1,2
	do i7=1,2
          cunit(i1,i3,i5,i7)=czero
	enddo
	enddo
	enddo
	enddo

       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=8d0*(vev4/(24d0*skin2)*
     & (cdelta_aij(0,0)+2d0*cdelta_aij(2,0))-
     &   5d0*vev4/(12d0*skin2)*
     & (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      cfact2=8d0*(5d0*vev4/(8d0*skin2)*
     &  (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
  
      cunit(i1,i3,i5,i7)=cfact1*cz1z3lor(i1,i3,i5,i7)
     & +cfact2*(cz1z2lor(i1,i3,i5,i7)+cz1z4lor(i1,i3,i5,i7))
  
      enddo
      enddo
      enddo
      enddo
  
      elseif (iout1*iout4.eq.1) then ! z1z4 <->z2z3
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pz4(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
	do i1=1,2
	do i3=1,2
	do i5=1,2
	do i7=1,2
          cunit(i1,i3,i5,i7)=czero
	enddo
	enddo
	enddo
	enddo

       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=8d0*(vev4/(24d0*skin2)*
     & (cdelta_aij(0,0)+2d0*cdelta_aij(2,0))-
     &   5d0*vev4/(12d0*skin2)*
     & (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      cfact2=8d0*(5d0*vev4/(8d0*skin2)*
     &  (cdelta_aij(0,2)+2d0*cdelta_aij(2,2)))
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
  
      cunit(i1,i3,i5,i7)=cfact1*cz1z4lor(i1,i3,i5,i7)
     & +cfact2*(cz1z2lor(i1,i3,i5,i7)+cz1z3lor(i1,i3,i5,i7))
  
      enddo
      enddo
      enddo
      enddo
  
      endif
  
      do i1=1,2
      do i3=1,2
      do i5=1,2
      do i7=1,2
  
      cunit(i1,i3,i5,i7)=-cunit(i1,i3,i5,i7)/elcharge2
  
      enddo
      enddo
      enddo
      enddo
  
  
      end
  
      subroutine unit_2z2w(pz1,sz1,cz1,pz2,sz2,cz2,
     &   pwp,swp,cwp,pwm,swm,cwm,cunit)
  
      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT COMPLEX*16 (c)
  
      include 'common.h'
      include 'common_unitarization.h'
*four momenta                                                           
      DIMENSION pz1(0:3),pz2(0:3),pwp(0:3),pwm(0:3),paux(0:3)
  
      type pol
         double complex e(0:3),ek0,v
      end type
      type(pol) cz1(2,2),cz2(2,2),cwp,cwm
  
      dimension cdelta_aij(0:2,0:3)
  
      dimension cdotzz(2,2,2,2),cdotwz1(2,2),cdotwz2(2,2)
      dimension cz1z2lor(2,2,2,2),cz1wplor(2,2,2,2),
     &  cz1wmlor(2,2,2,2),cunit(2,2,2,2)
  
* check if we have vv scattering                                        
      nout=0
      IF(sz1.GT.0.d0)THEN
        nout=nout+1
        iout1=1
      ELSE
        iout1=-1
      ENDIF
      IF(sz2.GT.0.d0)THEN
        nout=nout+1
        iout2=1
      ELSE
        iout2=-1
      ENDIF
      IF(swp.GT.0.d0)THEN
        nout=nout+1
        iout3=1
      ELSE
        iout3=-1
      ENDIF
      IF(swm.GT.0.d0)THEN
        nout=nout+1
        iout4=1
      ELSE
        iout4=-1
      ENDIF
  
      IF(nout.NE.2)THEN
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
        cunit(i1,i2,i3,i4)=czero
      enddo
      enddo
      enddo
      enddo
  
  
        RETURN
      ENDIF
  
* Lorentz variables (in terms of z1z2 -> w+w-)                          
  
* s                                                                     
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
* p.q -- p.q=cdotzz(i1,i2,i3,i4),p=cz1(i1,i2)%e,q=cz2(i3,i4)%e,bef=,aft=
      cdotzz(i1,i2,i3,i4)=(cz1(i1,i2)%e(0)*cz2(i3,i4)%e(0)-cz1(i
     & 1,i2)%e(1)*cz2(i3,i4)%e(1)-cz1(i1,i2)%e(2)*cz2(i3,i4)%e(2
     & )-cz1(i1,i2)%e(3)*cz2(i3,i4)%e(3))
      end do
      end do
      end do
      end do
* p.q -- p.q=cdotww,p=cwp%e,q=cwm%e,bef=,aft=
      cdotww=(cwp%e(0)*cwm%e(0)-cwp%e(1)*cwm%e(1)-cwp%e(2)*cwm%e
     & (2)-cwp%e(3)*cwm%e(3))
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
        cz1z2lor(i1,i2,i3,i4) =
     &  4d0*rmw2*rmz2*cdotzz(i1,i2,i3,i4)*cdotww/vev4
      enddo
      enddo
      enddo
      enddo
  
* t                                                                     
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cdotwz1(i1,i2),p=cz1(i1,i2)%e,q=cwp%e,bef=,aft=
      cdotwz1(i1,i2)=(cz1(i1,i2)%e(0)*cwp%e(0)-cz1(i1,i2)%e(1)*c
     & wp%e(1)-cz1(i1,i2)%e(2)*cwp%e(2)-cz1(i1,i2)%e(3)*cwp%e(3)
     & )
      end do
      end do
      do i3=1,2
      do i4=1,2
* p.q -- p.q=cdotwz2(i3,i4),p=cz2(i3,i4)%e,q=cwm%e,bef=,aft=
      cdotwz2(i3,i4)=(cz2(i3,i4)%e(0)*cwm%e(0)-cz2(i3,i4)%e(1)*c
     & wm%e(1)-cz2(i3,i4)%e(2)*cwm%e(2)-cz2(i3,i4)%e(3)*cwm%e(3)
     & )
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
       cz1wplor(i1,i2,i3,i4) =
     &  4d0*rmw2*rmz2*cdotwz1(i1,i2)*cdotwz2(i3,i4)/vev4
      enddo
      enddo
      enddo
      enddo
  
* u                                                                     
      do i1=1,2
      do i2=1,2
* p.q -- p.q=cdotwz1(i1,i2),p=cz1(i1,i2)%e,q=cwm%e,bef=,aft=
      cdotwz1(i1,i2)=(cz1(i1,i2)%e(0)*cwm%e(0)-cz1(i1,i2)%e(1)*c
     & wm%e(1)-cz1(i1,i2)%e(2)*cwm%e(2)-cz1(i1,i2)%e(3)*cwm%e(3)
     & )
      end do
      end do
      do i3=1,2
      do i4=1,2
* p.q -- p.q=cdotwz2(i3,i4),p=cz2(i3,i4)%e,q=cwp%e,bef=,aft=
      cdotwz2(i3,i4)=(cz2(i3,i4)%e(0)*cwp%e(0)-cz2(i3,i4)%e(1)*c
     & wp%e(1)-cz2(i3,i4)%e(2)*cwp%e(2)-cz2(i3,i4)%e(3)*cwp%e(3)
     & )
      end do
      end do
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
        cz1wmlor(i1,i2,i3,i4) =
     &  4d0*rmw2*rmz2*cdotwz1(i1,i2)*cdotwz2(i3,i4)/vev4
      enddo
      enddo
      enddo
      enddo
  
      if (iout1*iout2.eq.1) then ! z1z2 <->w+w-
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pz2(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
	do i1=1,2
	do i2=1,2
	do i3=1,2
	do i4=1,2
          cunit(i1,i2,i3,i4)=czero
	enddo
	enddo
	enddo
	enddo

       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=8d0*(vev4/(24d0*skin2)*
     & (cdelta_aij(0,0)-cdelta_aij(2,0))-
     &   5d0*vev4/(12d0*skin2)*
     & (cdelta_aij(0,2)-cdelta_aij(2,2)))
  
      cfact2=4d0*(5d0*vev4/(4d0*skin2)*
     &  (cdelta_aij(0,2)-cdelta_aij(2,2)))
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
  
      cunit(i1,i2,i3,i4)=cfact1*cz1z2lor(i1,i2,i3,i4)
     & +cfact2*(cz1wplor(i1,i2,i3,i4)+cz1wmlor(i1,i2,i3,i4))
  
      enddo
      enddo
      enddo
      enddo
  
      elseif (iout1*iout3.eq.1) then ! z2w+ ->z1w+ or z1w- ->z2w-
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pwp(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
	do i1=1,2
	do i2=1,2
	do i3=1,2
	do i4=1,2
          cunit(i1,i2,i3,i4)=czero
	enddo
	enddo
	enddo
	enddo

       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=4d0*(vev4/(8d0*skin2)*cdelta_aij(2,0)
     &   -5d0*vev4/(4d0*skin2)*cdelta_aij(2,2))
  
      cfact2=8d0*(-3d0*vev4/(16d0*skin2)*cdelta_aij(1,1)
     &  +15d0*vev4/(16d0*skin2)*cdelta_aij(2,2))
  
      cfact3=4d0*(3d0*vev4/(8d0*skin2)*cdelta_aij(1,1)
     &   +15d0*vev4/(8d0*skin2)*cdelta_aij(2,2))
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
  
      cunit(i1,i2,i3,i4)=cfact1*cz1wplor(i1,i2,i3,i4)
     & +cfact2*cz1z2lor(i1,i2,i3,i4)+cfact3*cz1wmlor(i1,i2,i3,i4)
  
      enddo
      enddo
      enddo
      enddo
  
      elseif (iout1*iout4.eq.1) then ! z2w- ->z1w- or z1w+ ->z2w+
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pwm(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
	do i1=1,2
	do i2=1,2
	do i3=1,2
	do i4=1,2
          cunit(i1,i2,i3,i4)=czero
	enddo
	enddo
	enddo
	enddo

       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=4d0*(vev4/(8d0*skin2)*cdelta_aij(2,0)
     &   -5d0*vev4/(4d0*skin2)*cdelta_aij(2,2))
  
      cfact2=8d0*(-3d0*vev4/(16d0*skin2)*cdelta_aij(1,1)
     &  +15d0*vev4/(16d0*skin2)*cdelta_aij(2,2))
  
      cfact3=4d0*(3d0*vev4/(8d0*skin2)*cdelta_aij(1,1)
     &   +15d0*vev4/(8d0*skin2)*cdelta_aij(2,2))
  
  
      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
  
      cunit(i1,i2,i3,i4)=cfact1*cz1wmlor(i1,i2,i3,i4)
     & +cfact2*cz1z2lor(i1,i2,i3,i4)+cfact3*cz1wplor(i1,i2,i3,i4)
  
      enddo
      enddo
      enddo
      enddo
  
      endif
  
* FORTRAN 90 can operate with arrays                                    
      cunit=-cunit/elcharge2
  
      end
  
  
      subroutine unit_2z2w_massless(pz1,sz1,cz1,pz2,sz2,cz2,
     &   pwp,swp,cwp,pwm,swm,cwm,cunit)
  
      IMPLICIT REAL*8 (a-b,d-h,o-z)
      IMPLICIT COMPLEX*16 (c)
  
      include 'common.h'
      include 'common_unitarization.h'
*four momenta                                                           
      DIMENSION pz1(0:3),pz2(0:3),pwp(0:3),pwm(0:3),paux(0:3)
  
      type pol
         double complex e(0:3),ek0,v
      end type
      type(pol) cz1(2),cz2(2),cwp,cwm
  
      dimension cdelta_aij(0:2,0:3)
  
      dimension cdotzz(2,2),cdotwz1(2),cdotwz2(2)
      dimension cz1z2lor(2,2),cz1wplor(2,2),cz1wmlor(2,2),cunit(2,2)
  
* check if we have vv scattering                                        
      nout=0
      IF(sz1.GT.0.d0)THEN
        nout=nout+1
        iout1=1
      ELSE
        iout1=-1
      ENDIF
      IF(sz2.GT.0.d0)THEN
        nout=nout+1
        iout2=1
      ELSE
        iout2=-1
      ENDIF
      IF(swp.GT.0.d0)THEN
        nout=nout+1
        iout3=1
      ELSE
        iout3=-1
      ENDIF
      IF(swm.GT.0.d0)THEN
        nout=nout+1
        iout4=1
      ELSE
        iout4=-1
      ENDIF
  
      IF(nout.NE.2)THEN
      do i1=1,2
      do i3=1,2
        cunit(i1,i3)=czero
      enddo
      enddo
  
  
        RETURN
      ENDIF
  
* Lorentz variables (in terms of z1z2 -> w+w-)                          
  
* s                                                                     
      do i1=1,2
      do i3=1,2
* p.q -- p.q=cdotzz(i1,i3),p=cz1(i1)%e,q=cz2(i3)%e,bef=,aft=
      cdotzz(i1,i3)=(cz1(i1)%e(0)*cz2(i3)%e(0)-cz1(i1)%e(1)*cz2(
     & i3)%e(1)-cz1(i1)%e(2)*cz2(i3)%e(2)-cz1(i1)%e(3)*cz2(i3)%e
     & (3))
      end do
      end do
* p.q -- p.q=cdotww,p=cwp%e,q=cwm%e,bef=,aft=
      cdotww=(cwp%e(0)*cwm%e(0)-cwp%e(1)*cwm%e(1)-cwp%e(2)*cwm%e
     & (2)-cwp%e(3)*cwm%e(3))
  
      do i1=1,2
      do i3=1,2
       cz1z2lor(i1,i3) =
     &  4d0*rmw2*rmz2*cdotzz(i1,i3)*cdotww/vev4
      enddo
      enddo
  
* t                                                                     
      do i1=1,2
* p.q -- p.q=cdotwz1(i1),p=cz1(i1)%e,q=cwp%e,bef=,aft=
      cdotwz1(i1)=(cz1(i1)%e(0)*cwp%e(0)-cz1(i1)%e(1)*cwp%e(1)-c
     & z1(i1)%e(2)*cwp%e(2)-cz1(i1)%e(3)*cwp%e(3))
      end do
      do i3=1,2
* p%q -- p%q=cdotwz2(i3),p=cz2(i3)%e,q=cwm%e,bef=,aft=
      cdotwz2(i3)=(cz2(i3)%e(0)*cwm%e(0)-cz2(i3)%e(1)*cwm%e(1)-c
     & z2(i3)%e(2)*cwm%e(2)-cz2(i3)%e(3)*cwm%e(3))
      end do
  
      do i1=1,2
      do i3=1,2
       cz1wplor(i1,i3) =
     &  4d0*rmw2*rmz2*cdotwz1(i1)*cdotwz2(i3)/vev4
      enddo
      enddo
  
* u                                                                     
      do i1=1,2
* p.q -- p.q=cdotwz1(i1),p=cz1(i1)%e,q=cwm%e,bef=,aft=
      cdotwz1(i1)=(cz1(i1)%e(0)*cwm%e(0)-cz1(i1)%e(1)*cwm%e(1)-c
     & z1(i1)%e(2)*cwm%e(2)-cz1(i1)%e(3)*cwm%e(3))
      end do
      do i3=1,2
* p.q -- p.q=cdotwz2(i3),p=cz2(i3)%e,q=cwp%e,bef=,aft=
      cdotwz2(i3)=(cz2(i3)%e(0)*cwp%e(0)-cz2(i3)%e(1)*cwp%e(1)-c
     & z2(i3)%e(2)*cwp%e(2)-cz2(i3)%e(3)*cwp%e(3))
      end do
  
      do i1=1,2
      do i3=1,2
       cz1wmlor(i1,i3) =
     &  4d0*rmw2*rmz2*cdotwz1(i1)*cdotwz2(i3)/vev4
      enddo
      enddo
  
      if (iout1*iout2.eq.1) then ! z1z2 <->w+w-
*       print*, 'zz<->ww'

* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pz2(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
	do i1=1,2
	do i3=1,2
          cunit(i1,i3)=czero
	enddo
	enddo

       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=8d0*(vev4/(24d0*skin2)*
     & (cdelta_aij(0,0)-cdelta_aij(2,0))-
     &   5d0*vev4/(12d0*skin2)*
     & (cdelta_aij(0,2)-cdelta_aij(2,2)))
  
      cfact2=4d0*(5d0*vev4/(4d0*skin2)*
     &  (cdelta_aij(0,2)-cdelta_aij(2,2)))
     
      do i1=1,2
      do i3=1,2
  
      cunit(i1,i3)=cfact1*cz1z2lor(i1,i3)
     & +cfact2*(cz1wplor(i1,i3)+cz1wmlor(i1,i3))
  
      enddo
      enddo
  
      elseif (iout1*iout3.eq.1) then ! z2w+ ->z1w+ or z1w- ->z2w-
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pwp(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
	do i1=1,2
	do i3=1,2
          cunit(i1,i3)=czero
	enddo
	enddo

       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=4d0*(vev4/(8d0*skin2)*cdelta_aij(2,0)
     &   -5d0*vev4/(4d0*skin2)*cdelta_aij(2,2))
  
      cfact2=8d0*(-3d0*vev4/(16d0*skin2)*cdelta_aij(1,1)
     &  +15d0*vev4/(16d0*skin2)*cdelta_aij(2,2))
  
      cfact3=4d0*(3d0*vev4/(8d0*skin2)*cdelta_aij(1,1)
     &   +15d0*vev4/(8d0*skin2)*cdelta_aij(2,2))
  
  
      do i1=1,2
      do i3=1,2
  
      cunit(i1,i3)=cfact1*cz1wplor(i1,i3)
     & +cfact2*cz1z2lor(i1,i3)+cfact3*cz1wmlor(i1,i3)
  
      enddo
      enddo
  
      elseif (iout1*iout4.eq.1) then ! z2w- ->z1w- or z1w+ ->z2w+
* kinematical s                                                         
       do m=0,3
        paux(m)=pz1(m)+pwm(m)
       enddo
* p.q -- p.q=skin,p=paux,q=paux,bef=,aft=
      skin=(paux(0)*paux(0)-paux(1)*paux(1)-paux(2)*paux(2)-paux
     & (3)*paux(3))
  
      skin2 =skin*skin

* minimum "s" for applying corrections to the amplitudes
      if (skin.lt.100d0) then      
	do i1=1,2
	do i3=1,2
          cunit(i1,i3)=czero
	enddo
	enddo

       return      
      endif
        
      call unit(skin,cdelta_aij)
  
      cfact1=4d0*(vev4/(8d0*skin2)*cdelta_aij(2,0)
     &   -5d0*vev4/(4d0*skin2)*cdelta_aij(2,2))
  
      cfact2=8d0*(-3d0*vev4/(16d0*skin2)*cdelta_aij(1,1)
     &  +15d0*vev4/(16d0*skin2)*cdelta_aij(2,2))
  
      cfact3=4d0*(3d0*vev4/(8d0*skin2)*cdelta_aij(1,1)
     &   +15d0*vev4/(8d0*skin2)*cdelta_aij(2,2))
  
  
      do i1=1,2
      do i3=1,2
  
      cunit(i1,i3)=cfact1*cz1wmlor(i1,i3)
     & +cfact2*cz1z2lor(i1,i3)+cfact3*cz1wplor(i1,i3)
  
      enddo
      enddo
  
      endif
      
      do i1=1,2
      do i3=1,2
       cunit(i1,i3)=-cunit(i1,i3)/elcharge2
      enddo
      enddo
      
      end
