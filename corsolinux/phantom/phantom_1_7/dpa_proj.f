      real*8 function vecmod(v)
      implicit real*8 (a-h,o-z)
c  given Lorentz vector v returns its 3-momentum absolute value 
      dimension v(0:3)
      vecmod = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      end
      
***********************************************************************
*                      SUBROUTINE doubleosp                           *
*                                                                     *
*     Purpose:  conserves the total VV' four-momentum,                *
*     the W/Z bosons direction in the VV' CoM reference frame         *
*     and the l/v/q direction in the correspondent W(Z)               *
*     rest frame.                                                     *
*     This routine requires one of the following final states (+jj)   *  
*                                                                     *
*     Warning: this osp works only for M4l/M2l2q > MV+MV'.            *
*                                                                     *
***********************************************************************

      subroutine doubleosp (plab,offactor)     
      implicit real*8 (a-h,o-z)     

      INCLUDE 'common.h'
      INCLUDE 'common_subproc.h'

      dimension plab(0:3,8)
      real*8 offactor
      real*8 rmr1,rmr2, wwinv, beta, fx, fy, gamr1, gamr2
      real*8 fact1,fact2,fact3
      integer i1id,i2id,i3id,i4id      
      dimension  e(0:3), ve(0:3), u(0:3), vu(0:3)
      dimension  We(0:3), Wu(0:3), ww(0:3)
      dimension  Wec(0:3), Wuc(0:3)
c projected leptons in lab frame
      dimension  e_on(0:3), ve_on(0:3), u_on(0:3), vu_on(0:3)
c leptons in WW CoM frame      
      dimension  e_w(0:3), ve_w(0:3), u_w(0:3), vu_w(0:3)
c leptons in WW CoM frame : projected  
      dimension  e_w_on(0:3), ve_w_on(0:3), u_w_on(0:3), vu_w_on(0:3)
c bosons in WW CoM frame : projected  
      dimension  We_w_on(0:3), Wu_w_on(0:3)
c bosons in lab frame : projected  
      dimension  We_l_on(0:3), Wu_l_on(0:3)
c leptons versors in WW CoM frame
      dimension  vew(1:3), vuw(1:3)
c We boson versor in WW CoM frame
      dimension  vW(1:3)
c auxiliary 3-vectors      
      dimension e3(1:3),We3(1:3),u3(1:3),Wu3(1:3)
     
c WW/WZ/ZZ      
      i1id = idosp(1)
      i2id = idosp(2)
      i3id = idosp(3)
      i4id = idosp(4)

      if((i1id+i2id).eq.0.and.(i1id*i2id).lt.0) then
         rmr1 = rmz
         gamr1 = gamz
      elseif(abs(i1id+i2id).eq.1.and.(i1id*i2id).lt.0)then
         rmr1 = rmw
         gamr1 = gamw
      else
         print*,'ERROR:  first two indices dont reconstruct a W/Z'
         stop
      endif   
      if((i3id+i4id).eq.0.and.(i3id*i4id).lt.0) then
         rmr2 = rmz
         gamr2 = gamz
      elseif(abs(i3id+i4id).eq.1.and.(i3id*i4id).lt.0)then
         rmr2 = rmw
         gamr2 = gamw
      else
         print*,'ERROR: second two indices dont reconstruct a W/Z'
         stop
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     performing the on shell projection
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO i=1,8
         IF ((idp(i)).EQ.i1id) THEN
            DO mu=0,3
               e(mu)=plab(mu,i) 
            ENDDO
         ELSEIF ((idp(i)).EQ.i2id) THEN 
            DO mu=0,3
               ve(mu)=plab(mu,i) 
            ENDDO	    
         ELSEIF ((idp(i)).EQ.i3id) THEN 
            DO mu=0,3
               u(mu)=plab(mu,i) 
            ENDDO	    
         ELSEIF ((idp(i)).EQ.i4id) THEN 
            DO mu=0,3
               vu(mu)=plab(mu,i) 
            ENDDO 
         ENDIF
      ENDDO
      
      We = e + ve 
      Wu = u + vu    
      ww = e + ve + u + vu 
      
c proietto i bosoni on-shell
      wwinv= sqrt(ww(0)*ww(0)-ww(1)*ww(1)-ww(2)*ww(2)-ww(3)*ww(3))
c      if(irr.eq.2.and.wwinv.le.175.d0)then
c         print*,"ERROR in OSP WZ: choose M4l > 175 GeV (M4l > Mw + Mz)"
c         stop
c      elseif(irr.eq.1.and.wwinv.le.185.d0)then
c         print*,"ERROR in OSP ZZ: choose M4l > 185 GeV (M4l > 2 Mz)"
c         stop
c      elseif(irr.eq.0.and.wwinv.le.165.d0)then
c         print*,"ERROR in OSP WW: choose M4l > 165 GeV (M4l > 2 Mw)"
c         stop
c      endif

    
      fact1 = 0.5d0*wwinv
      fact2 = 0.5d0*(rmr2**2-rmr1**2)/wwinv
      fact3 = 0.5d0*(rmr2**2+rmr1**2)
      beta= sqrt(fact1**2 + fact2**2 - fact3)
      call boostinv (We, ww, Wec)
      call boostinv (Wu, ww, Wuc)
      vW(1) = Wec(1)/vecmod(Wec)
      vW(2) = Wec(2)/vecmod(Wec)
      vW(3) = Wec(3)/vecmod(Wec)                 
      We_w_on(0) =fact1-fact2; We_w_on(1) = beta*vW(1)
      We_w_on(2) = beta*vW(2); We_w_on(3) = beta*vW(3)
      Wu_w_on(0) =fact1+fact2; Wu_w_on(1) = -beta*vW(1)
      Wu_w_on(2) = -beta*vW(2); Wu_w_on(3) = -beta*vW(3)
      call boost (We_w_on, ww, We_l_on)
      call boost (Wu_w_on, ww, Wu_l_on)
c ora i bosoni sono proiettati on-shell           
c proietto i leptoni
      Winv=We(0)*We(0)-We(1)*We(1)-We(2)*We(2)-We(3)*We(3)
      if(Winv.le.0.d0)then
         print*,"ERROR: W invariant mass < 0"
         stop
      endif
      call boostinv (e, We, e_w)
      Winv=Wu(0)*Wu(0)-Wu(1)*Wu(1)-Wu(2)*Wu(2)-Wu(3)*Wu(3)
      if(Winv.le.0.d0)then
         print*,"ERROR: Z invariant mass < 0"
         stop
      endif
      call boostinv (u, Wu, u_w)
c porto i leptoni nel sistema in cui i rispettivi W sono lungo l'asse z
      e3(1) = e_w(1); e3(2) = e_w(2); e3(3) = e_w(3)
      u3(1) = u_w(1); u3(2) = u_w(2); u3(3) = u_w(3)
      We3(1) = We(1); We3(2) = We(2); We3(3) = We(3)
      Wu3(1) = Wu(1); Wu3(2) = Wu(2); Wu3(3) = Wu(3)
      call rotinv(e3,We3,e3)
      call rotinv(u3,Wu3,u3)
      e_w(1) = e3(1) ; e_w(2) = e3(2) ;  e_w(3) = e3(3)
      u_w(1) = u3(1) ; u_w(2) = u3(2) ;  u_w(3) = u3(3)
      
c++++++++++++++++++++++++++++++++++++
c normalize to rmr*0.5    
       vew(1) = e_w(1)/vecmod(e_w)
       vew(2) = e_w(2)/vecmod(e_w)
       vew(3) = e_w(3)/vecmod(e_w)
       vuw(1) = u_w(1)/vecmod(u_w)
       vuw(2) = u_w(2)/vecmod(u_w)
       vuw(3) = u_w(3)/vecmod(u_w)                               
       fx = rmr1*0.5d0
       fy = rmr2*0.5d0
       e_w_on(0) = fx; e_w_on(1) = fx*vew(1)
       e_w_on(2) = fx*vew(2); e_w_on(3) = fx*vew(3)     
       u_w_on(0) = fy; u_w_on(1) = fy*vuw(1)
       u_w_on(2) = fy*vuw(2); u_w_on(3) = fy*vuw(3)        
       e_w = e_w_on
       u_w = u_w_on
       
c++++++++++++++++++++++++++++++++++++
c ruoto in modo che l'asse zeta sia nella direzione dei rispettivi bosoni on-shell
       e3(1) = e_w(1); e3(2) = e_w(2); e3(3) = e_w(3)
       u3(1) = u_w(1); u3(2) = u_w(2); u3(3) = u_w(3)
       We3(1) = We_l_on(1); We3(2) = We_l_on(2); We3(3) = We_l_on(3)
       Wu3(1) = Wu_l_on(1); Wu3(2) = Wu_l_on(2); Wu3(3) = Wu_l_on(3)
       call rot(e3,We3,e3)
       call rot(u3,Wu3,u3)
       e_w(1) = e3(1) ; e_w(2) = e3(2) ;  e_w(3) = e3(3)
       u_w(1) = u3(1) ; u_w(2) = u3(2) ;  u_w(3) = u3(3)
c boost nel sistema dei W on shell nel lab
      call boost (e_w, We_l_on, e_on)
      call boost (u_w, Wu_l_on, u_on)

      ve_on = We_l_on - e_on
      vu_on = Wu_l_on - u_on      
      
c ora i leptoni sono proiettati                    

         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
      DO i=1,8
         IF ((idp(i)).EQ.i1id) THEN
            DO mu=0,3
               plab(mu,i)=e_on(mu)
            ENDDO
         ELSEIF ((idp(i)).EQ.i2id) THEN
            DO mu=0,3
               plab(mu,i)=ve_on(mu) 
            ENDDO	    
         ELSEIF ((idp(i)).EQ.i3id) THEN
            DO mu=0,3
               plab(mu,i)=u_on(mu) 
            ENDDO	    
         ELSEIF ((idp(i)).EQ.i4id) THEN  
            DO mu=0,3
               plab(mu,i)=vu_on(mu) 
            ENDDO 
         ENDIF
      ENDDO
      offWe = (e(0)+ve(0))**2.-(e(1)+ve(1))**2.-
     &     (e(2)+ve(2))**2.-(e(3)+ve(3))**2.
      offWu = (u(0)+vu(0))**2.-(u(1)+vu(1))**2.-
     &     (u(2)+vu(2))**2.-(u(3)+vu(3))**2.
     
      offactor = ((gamr1*rmr1)**2*(gamr2*rmr2)**2)/
     &    (((offWe-rmr1**2)**2.+(gamr1*rmr1)**2.)*
     &     ((offWu-rmr2**2)**2.+(gamr2*rmr2)**2.))
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      return
      end
