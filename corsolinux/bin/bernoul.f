************************************************************************
*                        SUBROUTINE BERNOULLI                          *
*                                                                      *
* Purpose:                                                             *
* It computes Bernoulli numbers for dilog (Spence) function evaluation.*
*                                                                      *
* A call to this subroutine (call Bernoulli(n) must be made before     *
* using the function rli2(x,n).  n must be the sambe in both case and  *
* not greater than 20. It is the order at which the series is truncated*
*                                                                      *
*                                                                      *
************************************************************************

      Subroutine Bernoulli(n)
      implicit real*8(a-h,o-z)
      dimension Be(0:20),B(0:20)
      parameter (Zeta2=1.64493406684823d0)
      common/coeff/B

* This routine computes the Bernoulli numbers for rli2
* B(n)= Bernoulli_numbers(2*n)

      Be(0)=1.d0
      Be(1)=1.d0/6.d0
      Be(2)=-1.d0/30.d0
      Be(3)=1.d0/42.d0
      Be(4)=Be(2)
      Be(5)=5.d0/66.d0
      Be(6)=-691.d0/2730.d0
      Be(7)=7.d0*Be(1)
      Be(8)=-3617.d0/510.d0
      Be(9)=43867.d0/798.d0
      Be(10)=-174611.d0/330.d0
      Be(11)=854513.d0/138.d0
      Be(12)=-236364091.d0/2730.d0
      Be(13)=8553103.d0/6.d0
      Be(14)=-23749461029.d0/870.d0
      Be(15)=8615841276005.d0/14322.d0
      Be(16)=-7709321041217.d0/510.d0
      Be(17)=2577687858367.d0/6.d0
      Be(18)=-26315271553053477373.d0/1919190.d0
      Be(19)=2929993913841559.d0/6.d0
      Be(20)=-261082718496449122051.d0/13530.d0

* Divido per (2n+1)!
      fn1=1.d0
      fn2=1.d0
      do i=1,n
        fn1=fn1*(2.d0*i+1)
        fn2=fn2*2.d0*i
        B(i)=Be(i)/fn1
        B(i)=B(i)/fn2
      enddo
        B(0)=1.d0
      return
      end

      function rli2(x,n)
*n is the order to which it is computed the series. Must be at most 20 
* with the above Bernoulli numbers, and not greater that the one used
* for the call to Bernoulli subroutine.   
      implicit real*8(a-h,o-z)
      dimension B(0:20)
      parameter (Zeta2=1.64493406684823d0)
      common/coeff/B

      if(x.gt.1.d0) then
        print*,'Variable out of range'
        stop
      endif

      if(abs(x).le.0.5d0) then
        t=-log(1.d0-x)
        rli2= -t/4.d0
        f=1.d0
        do i=0,n
          rli2=rli2+B(i)*f
          f=f*t*t
        enddo
        rli2=rli2*t
      elseif(x.gt.0.5d0) then
        q=1.d0-x
        t=-log(1.d0-q)
        rli2= -t/4.d0
        f=1.d0
        do i=0,n
          rli2=rli2+B(i)*f
          f=f*t*t
        enddo
        rli2=-rli2*t+Zeta2
        if (x.ne.1.d0) then
          rli2=rli2-log(x)*log(q)
        endif
      elseif(x.ge.-1.d0) then
        q=-x/(1.d0-x)
        t=-log(1.d0-q)
        rli2= -t/4.d0
        f=1.d0
        do i=0,n
          rli2=rli2+B(i)*f
          f=f*t*t
        enddo
        rli2=-rli2*t-log(1-x)**2/2.d0
      elseif(x.ge.-2.d0) then
        q=1.d0/(1.d0-x)
        t=-log(1.d0-q)
        rli2= -t/4.d0
        f=1.d0
        do i=0,n
          rli2=rli2+B(i)*f
          f=f*t*t
        enddo
        rli2=rli2*t-log(-x)*log(1-x)+log(1-x)**2/2.d0-Zeta2
      else
        q=1.d0/x
        t=-log(1.d0-q)
        rli2= -t/4.d0
        f=1.d0
        do i=0,n
          rli2=rli2+B(i)*f
          f=f*t*t
        enddo
        rli2=-rli2*t-log(-x)**2/2.d0-Zeta2
      endif
      return
      end

