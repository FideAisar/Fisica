***********************************************************************
*                        SUBROUTINE COUPLING                          *
*                                                                     *
*                                                                     *
*  Purpose:  It contains DATA (for masses, W and Z widths and         *
*                                                     alfa_couplings) *
*            It computes Higgs and top widths, and couplings between  *
*            particles.                                               *
*                                                                     *
*  Call to Subroutine Bernoulli                                       *
*                                                                     *
***********************************************************************


      PROGRAM H2hh_Width

      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)

      DATA rmw/80.40d0/, rmz/91.187d0/, rmt/175.d0/, rmb/4.8d0/,
     &     rmc/0.75d0/,rmtau/1.78d0/,
     &     gamw/0.2042774d+01/, gamz/0.25007d+01/
      DATA gf/1.16639d-05/, alfainv_me/137.0359895d0/,
     &  alfas_w/0.1225d0/, alfas_z/0.123d0/,
     &  alfas_t/0.100d0/, alfas_h/0.100d0/

      DATA pi/3.141592653589793238462643d0/

      DATA czero/(0.d0,0.d0)/, cuno/(1.d0,0.d0)/, cim/(0.d0,1.d0)/


      WRITE(*,*)'>>> Enter rmh in GeV'
      READ(*,*) rmh

      WRITE(*,*)'>>> Enter rmhh in GeV'
      READ(*,*) rmhh

      WRITE(*,*) '>>> Enter ghhfactor = sin(alfa)  in Sinlet Model'
      READ(*,*) ghhfactor

      WRITE(*,*) '>>> Enter tgbeta'
      READ(*,*) tgbeta

      ghfactor = sqrt(1.d0-ghhfactor**2)
      print*,'sinalf=ghhfactor=',ghhfactor,
     &           ' cosalf=ghfactor=',ghfactor
      
      alfa=1.d0/alfainv_me
      esquared=alfa*4.d0*pi
      gamH2hh=esquared*rmhh**3/(128.d0*pi*rmw**2*(1.d0-rmw**2/rmz**2))*
     &      sqrt(1.d0-4.d0*rmh**2/rmhh**2)*
     &      (1.d0+2.d0*rmh**2/rmhh**2)**2*
     &      ghfactor**2*ghhfactor**2*(ghfactor+ghhfactor*tgbeta)**2


      print*,'rmh  =',rmh,'  rmhh =',rmhh
      print*,'ghhfactor = sin(alfa) =',ghhfactor
      print*,'tgbeta=',tgbeta

      print*,'gamH2hh=',gamH2hh,' GeV'

      END



