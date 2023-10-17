      REAL*8 v1(0:2),v2(0:2),v3(0:2),v4(0:2),vout(0:2),r1,r2,a(0:2,0:2)

      v1(0)=5.d0
      v1(1)=1.d0
      v1(2)=2.d0

      v2(0)=8.d0
      v2(1)=-2.d0
      v2(2)=-4.d0

      a(0,0)= 0.d0
      a(0,1)= 1.d0
      a(0,2)= 2.d0

      a(1,0)= 10.d0
      a(1,1)= 11.d0
      a(1,2)= 12.d0

      v3= a(0,:)
      write(*,*)'v3: ',v3

      v4= a(1,:)
      write(*,*)'v4: ',v4

      write(*,*)'a(1,:): ',a(1,:)

c      r1 = v3*a(1,:) doesn't work
c v1*v2 returns a vector, multiplying v1 and v2 component by component

      vout = v3*v4
      WRITE(*,*)'v3.v4: ',vout
      r1 = SUM(vout)
      WRITE(*,*)'r1: ',r1

c      r1=(v1(0)*v1(0)-v1(1)*v1(1)-v1(2)*v1(2)-v1(3)*v1(3))**0.5
c      r2=sqrt(v1(0)*v1(0)-v1(1)*v1(1)-v1(2)*v1(2)-v1(3)*v1(3))

c      WRITE(*,*)r1,r2

      end
