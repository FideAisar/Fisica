c
c     IP    = PDG CODE FOR THE PARTICLE TO DECAY
c     IMODE = DECAY MODE
c     ND    = NUMBER OF DECAY PRODUCTS
c
      INTEGER       IP,IMODE,ND	
      COMMON/DEC_ID/IP,IMODE,ND
c
c     masses and couplings for the decay routines
c
      REAL*8          M1,M2,M3,M4,M5,MV,GV
      COMMON/m_decays/M1,M2,M3,M4,M5,MV,GV
      DOUBLE COMPLEX       GX,GXX(2),GXX1(2)  
      COMMON/WEAK_COUPL/   GX,GXX   ,GXX1
c
c    MC width, analytic width and best branching ratio
c     	
      real*8          MC_width,calc_width,bratio
      common/to_width/MC_width,calc_width,bratio
