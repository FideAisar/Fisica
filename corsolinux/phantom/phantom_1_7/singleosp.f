***********************************************************************
*                        SUBROUTINE singleosp                         *
*                                                                     *
* Purpose: It performs the on-shell projection of W/Z bosons          *
*          in processes p p > W/Z(->ll/lv/qq) + X	 	      *
*	   This projection conserves 	        	              *
*		- the X system four-momentum,			      *
*		- the W/Z three-momentum in the lab. ref. frame,      *	
*		- the l/v/q direction in the W/Z  rest frame.         *
*	   This projection modifies the initial parton four-momenta,  *
*          in order to conserve the total four-momentum.              *
*                                                                     * 
***********************************************************************
  
      subroutine singleosp (plab,offactor)    
      implicit real*8 (a-h,o-z)
      
      INCLUDE 'common.h'
      INCLUDE 'common_subproc.h'
      
      dimension plab(0:3,8)
      real*8 wmass,DeltaW,rmr,gamr
      real*8 offactor,offW
      integer i,j
      integer i1id,i2id
c particles in lab frame
      dimension rl(0:3),rv(0:3),W(0:3),q1(0:3),q2(0:3),rI(0:3,2)
      dimension rlo(0:3),rvo(0:3),wo(0:3),q1o(0:3),q2o(0:3),rIo(0:3,2)    
c temporary particles fourvector for calculations      
      dimension ec(0:3)
c lepton versor in W CoM frame, W boson versor in overall CoM frame 
      dimension Vl(1:3)
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      i1id = idosp(1)
      i2id = idosp(2)
c      print*,'singleosp.f: i1id = ',i1id
c      print*,'singleosp.f: i2id = ',i2id
c     Z single osp
      if((i1id+i2id).eq.0.and.(i1id*i2id).lt.0) then
         wmass = rmz
         rmr = rmz
         gamr = gamz
c     W single osp
      elseif(abs(i1id+i2id).eq.1.and.(i1id*i2id).lt.0)then
         wmass = rmw
         rmr = rmw
         gamr = gamw
      else
         print*,'ERROR:  the two indices dont reconstruct a W/Z'
         stop   
      endif 
      
c      print*,'singleosp.f: wmass = ',wmass
c      print*,'singleosp.f: rmr = ',rmr
c      print*,'singleosp.f: gamr = ',gamr
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	 j=1
         DO i=1,8          
           IF ((idp(i)).eq.i1id) THEN
            DO mu=0,3
             rl(mu)=plab(mu,i)
            ENDDO
           ELSEIF ((idp(i)).EQ.i2id) THEN
            DO mu=0,3
	     rv(mu)=plab(mu,i) 
	    ENDDO
           ELSEIF (idp_inout(i).EQ. -1) THEN         
            DO mu=0,3
             rI(mu,j)=-plab(mu,i)
            ENDDO 
            j=j+1                              
           ENDIF          
         ENDDO    
         
c         print*,'singleosp.f: MOMENTA BEFORE PROJECTION'
c         print*,'singleosp.f: DECAY PRODUCTS TO PROJECT'
c         print*,'rl(mu)=',rl(0),rl(1),rl(2),rl(3)          
c	 print*,'rv(mu)=',rv(0),rv(1),rv(2),rv(3) 
c	 print*,'singleosp.f: INITIAL MOMENTA'
c	 print*,'rI(mu,1)=',rI(0,1),rI(1,1),rI(2,1),rI(3,1) 
c	 print*,'rI(mu,2)=',rI(0,2),rI(1,2),rI(2,2),rI(3,2)
	
	
	
         do mu=0,3
            q1(mu) = rI(mu,1)
            q2(mu) = rI(mu,2)
         enddo	
         
         W = rl + rv         
         call boostinv(rl,W,ec)            
         Vl(1) = ec(1)/((ec(1)*ec(1)+ec(2)*ec(2)+ec(3)*ec(3))**(0.5))
         Vl(2) = ec(2)/((ec(1)*ec(1)+ec(2)*ec(2)+ec(3)*ec(3))**(0.5))
         Vl(3) = ec(3)/((ec(1)*ec(1)+ec(2)*ec(2)+ec(3)*ec(3))**(0.5))
         ec(0)=wmass*0.5; ec(1)=wmass*0.5*Vl(1)
         ec(2)=wmass*0.5*Vl(2); ec(3)=wmass*0.5*Vl(3)     
         wo(0)=(wmass*wmass+W(1)*W(1)+W(2)*W(2)+W(3)*W(3))**0.5
         do i=1,3
            wo(i)=W(i)
         enddo
         call boost(ec,wo,rlo)
         rvo = wo - rlo
         
         DeltaW = wo(0) - W(0)
         
         IF (q1(3) .LT. 0.) THEN
            q1o(0)= q1(0)+0.5*DeltaW
            q1o(1)= 0.
            q1o(2)= 0.
            q1o(3)= -q1o(0)
            q2o(0)= q2(0)+0.5*DeltaW
            q2o(1)= 0.
            q2o(2)= 0.
            q2o(3)= q2o(0)
         ELSE
            q1o(0)= q1(0)+0.5*DeltaW
            q1o(1)= 0.
            q1o(2)= 0.
            q1o(3)= q1o(0)
            q2o(0)= q2(0)+0.5*DeltaW
            q2o(1)= 0.
            q2o(2)= 0.
            q2o(3)= -q2o(0)      
         ENDIF
         
         do mu=0,3
            rIo(mu,1) = q1o(mu)
            rIo(mu,2) = q2o(mu)
         enddo    
         
         
c         print*,'singleosp.f: MOMENTA AFTER PROJECTION'
c         print*,'singleosp.f: DECAY PRODUCTS PROJECTED'
c         print*,'rl(mu)=',rlo(0),rlo(1),rlo(2),rlo(3)          
c	 print*,'rv(mu)=',rvo(0),rvo(1),rvo(2),rvo(3) 
c	 print*,'singleosp.f: INITIAL MOMENTA'
c	 print*,'rI(mu,1)=',rIo(0,1),rIo(1,1),rIo(2,1),rIo(3,1) 
c	 print*,'rI(mu,2)=',rIo(0,2),rIo(1,2),rIo(2,2),rIo(3,2)
         
                
	 j=1
         DO i=1,8          
           IF ((idp(i)).EQ. i1id) THEN
            DO mu=0,3
             plab(mu,i)=rlo(mu)
            ENDDO
           ELSEIF ((idp(i)).EQ. i2id) THEN
            DO mu=0,3
	     plab(mu,i)=rvo(mu) 
	    ENDDO
           ELSEIF (idp_inout(i).EQ. -1) THEN           
            DO mu=0,3
             plab(mu,i)=-rIo(mu,j)
            ENDDO 
            j=j+1                              
           ENDIF          
         ENDDO                
	 offW = (rl(0)+rv(0))**2.-(rl(1)+rv(1))**2.-
     &	(rl(2)+rv(2))**2.-(rl(3)+rv(3))**2.
         offactor = ((gamr*rmr)**2.)/
     &  ((offW-rmr*rmr)**2.+(gamr*rmr)**2.)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
c         print*,'singleosp.f: offactor = ',offactor            
      return
      end
