C     iunittype= 0 no unitarization
C            1 kmatrix
C            2 pade
C            3 nd
C            4 largen     
      COMMON/unittype/i_unitarize,iunittype
      COMMON/unitnlo/inlo
      COMMON/unitreson/ireson,isigma,irho,iphi,iff,ia
      COMMON/unitresonpar/rg_sigma,rm_sigma,gam_sigma,
     &  rg_rho,rm_rho,gam_rho,rg_phi,rm_phi,gam_phi,
     &  rg_ff,rm_ff,gam_ff,rg_a,rm_a,gam_a
      COMMON/uniconst/vev,vev2,vev4,rmu,alpha4,alpha5
      COMMON/unitother/rmnd

