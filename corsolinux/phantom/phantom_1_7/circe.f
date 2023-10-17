c circe.f -- canonical beam spectra for linear collider physics
c   Copyright (C) 1996-1999 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
c
c   Circe is free software; you can redistribute it and/or modify it
c   under the terms of the GNU General Public License as published by
c   the Free Software Foundation; either version 2, or (at your option)
c   any later version.
c
c   Circe is distributed in the hope that it will be useful, but
c   WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c   GNU General Public License for more details.
c
c   You should have received a copy of the GNU General Public License
c   along with this program; if not, write to the Free Software
c   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
      double precision function circe (x1, x2, p1, p2)
      implicit none
      double precision x1, x2
      integer p1, p2
      double precision circee, circeg, circgg
      integer ELECTR, POSITR, PHOTON
      parameter (ELECTR =  11)
      parameter (POSITR = -11)
      parameter (PHOTON =  22)
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      circe = -1.0
      if (abs(p1) .eq. ELECTR) then
         if (abs(p2) .eq. ELECTR) then
            circe = circee (x1, x2)
         elseif (p2 .eq. PHOTON) then
            circe = circeg (x1, x2)
         endif
      elseif (p1 .eq. PHOTON) then
         if (abs(p2) .eq. ELECTR) then
            circe = circeg (x2, x1)
         elseif (p2 .eq. PHOTON) then
            circe = circgg (x1, x2)
         endif
      endif
      end
      subroutine circes (xx1m, xx2m, xroots, xacc, xver, xrev, xchat)
      implicit none
      double precision xx1m, xx2m, xroots
      double precision beta
      integer xacc, xver, xrev, xchat
      integer SBAND, TESLA, XBAND
      parameter (SBAND  =  1, TESLA  =  2, XBAND  =  3)
      integer JLCNLC
      parameter (JLCNLC =  3)
      integer SBNDEE, TESLEE, XBNDEE
      parameter (SBNDEE =  4, TESLEE =  5, XBNDEE =  6)
      integer NACC
      parameter (NACC = 6)
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      character*60 msgbuf
      character*6 accnam(NACC)
      integer ver34
      integer GEV090, GEV170, GEV350, GEV500, GEV800, TEV1, TEV16
      parameter (GEV090 = -1, GEV170 = 0, GEV350 = 1, GEV500 = 2,
     $           GEV800 =  3, TEV1   = 4, TEV16  = 5)
      integer A1NEGY, A1NREV
      parameter (A1NEGY = 5, A1NREV = 5)
      integer i
      real xa1lum(A1NEGY,NACC,0:A1NREV)
      real xa1(0:7,A1NEGY,NACC,0:A1NREV)
      integer A3NEGY, A3NREV
      parameter (A3NEGY = 5, A3NREV = 5)
      real xa3lum(A3NEGY,NACC,0:A3NREV)
      real xa3(0:7,A3NEGY,NACC,0:A3NREV)
      integer A5NEGY, A5NREV
      parameter (A5NEGY = 5, A5NREV = 1)
      real xa5lum(A5NEGY,NACC,0:A5NREV)
      real xa5(0:7,A5NEGY,NACC,0:A5NREV)
      integer A6NEGY, A6NREV
      parameter (A6NEGY = 2, A6NREV = 1)
      real xa6lum(GEV090:A6NEGY,NACC,0:A6NREV)
      real xa6(0:7,GEV090:A6NEGY,NACC,0:A6NREV)
      double precision eloval, ehival
      integer A7NEGY, A7NREV
      parameter (A7NEGY = TEV1, A7NREV = 1)
      real xa7lum(GEV090:A7NEGY,NACC,0:A7NREV)
      real xa7(0:7,GEV090:A7NEGY,NACC,0:A7NREV)
      data accnam(SBAND)  /'SBAND'/
      data accnam(TESLA)  /'TESLA'/
      data accnam(JLCNLC) /'JLCNLC'/
      data accnam(SBNDEE) /'SBNDEE'/
      data accnam(TESLEE) /'TESLEE'/
      data accnam(XBNDEE) /'XBNDEE'/
      data xa1lum(GEV500,SBAND,1) /  5.212299E+01 /
      data (xa1(i,GEV500,SBAND,1),i=0,7) /
     $    .39192E+00,   .66026E+00,   .11828E+02,  -.62543E+00, 
     $    .52292E+00,  -.69245E+00,   .14983E+02,   .65421E+00 /
      data xa1lum(GEV500,TESLA,1) /  6.066178E+01 /
      data (xa1(i,GEV500,TESLA,1),i=0,7) /
     $    .30196E+00,   .12249E+01,   .21423E+02,  -.57848E+00, 
     $    .68766E+00,  -.69788E+00,   .23121E+02,   .78399E+00 /
      data xa1lum(GEV500,XBAND,1) /  5.884699E+01 /
      data (xa1(i,GEV500,XBAND,1),i=0,7) /
     $    .48594E+00,   .52435E+00,   .83585E+01,  -.61347E+00, 
     $    .30703E+00,  -.68804E+00,   .84109E+01,   .44312E+00 /
      data xa1lum(TEV1,SBAND,1)   /  1.534650E+02 /
      data (xa1(i,TEV1,SBAND,1),i=0,7) /
     $    .24399E+00,   .87464E+00,   .66751E+01,  -.56808E+00, 
     $    .59295E+00,  -.68921E+00,   .94232E+01,   .83351E+00 /
      data xa1lum(TEV1,TESLA,1)   /  1.253381E+03 /
      data (xa1(i,TEV1,TESLA,1),i=0,7) /
     $    .39843E+00,   .70097E+00,   .11602E+02,  -.61061E+00, 
     $    .40737E+00,  -.69319E+00,   .14800E+02,   .51382E+00 /
      data xa1lum(TEV1,XBAND,1)   /  1.901783E+02 /
      data (xa1(i,TEV1,XBAND,1),i=0,7) /
     $    .32211E+00,   .61798E+00,   .28298E+01,  -.54644E+00, 
     $    .45674E+00,  -.67301E+00,   .41703E+01,   .74536E+00 /
      data (xa1lum(GEV350,i,1),i=1,NACC) / NACC*-1d0 /
      data (xa1lum(GEV800,i,1),i=1,NACC) / NACC*-1d0 /
      data (xa1lum(GEV500,i,1),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV1,i,1),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV16,i,1),i=1,NACC) / 6*-1d0 /
      data xa1lum(GEV500,SBAND,2) /   .31057E+02 /
      data (xa1(i,GEV500,SBAND,2),i=0,7) /
     $    .38504E+00,   .79723E+00,   .14191E+02,  -.60456E+00, 
     $    .53411E+00,  -.68873E+00,   .15105E+02,   .65151E+00 /
      data xa1lum(TEV1,SBAND,2) /   .24297E+03 /
      data (xa1(i,TEV1,SBAND,2),i=0,7) /
     $    .24374E+00,   .89466E+00,   .70242E+01,  -.56754E+00, 
     $    .60910E+00,  -.68682E+00,   .96083E+01,   .83985E+00 /
      data xa1lum(GEV350,TESLA,2) /   .73369E+02 /
      data (xa1(i,GEV350,TESLA,2),i=0,7) /
     $    .36083E+00,   .12819E+01,   .37880E+02,  -.59492E+00, 
     $    .69109E+00,  -.69379E+00,   .40061E+02,   .65036E+00 /
      data xa1lum(GEV500,TESLA,2) /   .10493E+03 /
      data (xa1(i,GEV500,TESLA,2),i=0,7) /
     $    .29569E+00,   .11854E+01,   .21282E+02,  -.58553E+00, 
     $    .71341E+00,  -.69279E+00,   .24061E+02,   .77709E+00 /
      data xa1lum(GEV800,TESLA,2) /   .28010E+03 /
      data (xa1(i,GEV800,TESLA,2),i=0,7) /
     $    .22745E+00,   .11265E+01,   .10483E+02,  -.55711E+00, 
     $    .69579E+00,  -.69068E+00,   .13093E+02,   .89605E+00 /
      data xa1lum(TEV1,TESLA,2) /   .10992E+03 /
      data (xa1(i,TEV1,TESLA,2),i=0,7) /
     $    .40969E+00,   .66105E+00,   .11972E+02,  -.62041E+00, 
     $    .40463E+00,  -.69354E+00,   .14669E+02,   .51281E+00 /
      data xa1lum(GEV500,XBAND,2) /   .35689E+02 /
      data (xa1(i,GEV500,XBAND,2),i=0,7) /
     $    .48960E+00,   .46815E+00,   .75249E+01,  -.62769E+00, 
     $    .30341E+00,  -.68754E+00,   .85545E+01,   .43453E+00 /
      data xa1lum(TEV1,XBAND,2) /   .11724E+03 /
      data (xa1(i,TEV1,XBAND,2),i=0,7) /
     $    .31939E+00,   .62415E+00,   .30763E+01,  -.55314E+00, 
     $    .45634E+00,  -.67089E+00,   .41529E+01,   .73807E+00 /
      data xa1lum(GEV350,SBAND,2) / -1d0 /
      data xa1lum(GEV350,XBAND,2) / -1d0 /
      data xa1lum(GEV800,SBAND,2) / -1d0 /
      data xa1lum(GEV800,XBAND,2) / -1d0 /
      data (xa1lum(GEV350,i,2),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV500,i,2),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV800,i,2),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV1,i,2),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV16,i,2),i=1,NACC) / 6*-1d0 /
      data xa1lum(GEV500,SBAND, 3) /   .31469E+02 /
      data (xa1(i,GEV500,SBAND, 3),i=0,7) /
     $    .38299E+00,   .72035E+00,   .12618E+02,  -.61611E+00, 
     $    .51971E+00,  -.68960E+00,   .15066E+02,   .63784E+00 /
      data xa1lum(TEV1,  SBAND, 3) /   .24566E+03 /
      data (xa1(i,TEV1,  SBAND, 3),i=0,7) /
     $    .24013E+00,   .95763E+00,   .69085E+01,  -.55151E+00, 
     $    .59497E+00,  -.68622E+00,   .94494E+01,   .82158E+00 /
      data xa1lum(GEV350,TESLA, 3) /   .74700E+02 /
      data (xa1(i,GEV350,TESLA, 3),i=0,7) /
     $    .34689E+00,   .12484E+01,   .33720E+02,  -.59523E+00, 
     $    .66266E+00,  -.69524E+00,   .38488E+02,   .63775E+00 /
      data xa1lum(GEV500,TESLA, 3) /   .10608E+03 /
      data (xa1(i,GEV500,TESLA, 3),i=0,7) /
     $    .28282E+00,   .11700E+01,   .19258E+02,  -.58390E+00, 
     $    .68777E+00,  -.69402E+00,   .23638E+02,   .75929E+00 /
      data xa1lum(GEV800,TESLA, 3) /   .28911E+03 /
      data (xa1(i,GEV800,TESLA, 3),i=0,7) /
     $    .21018E+00,   .12039E+01,   .96763E+01,  -.54024E+00, 
     $    .67220E+00,  -.69083E+00,   .12733E+02,   .87355E+00 /
      data xa1lum(TEV1,  TESLA, 3) /   .10936E+03 /
      data (xa1(i,TEV1,  TESLA, 3),i=0,7) /
     $    .41040E+00,   .68099E+00,   .11610E+02,  -.61237E+00, 
     $    .40155E+00,  -.69073E+00,   .14698E+02,   .49989E+00 /
      data xa1lum(GEV500,XBAND, 3) /   .36145E+02 /
      data (xa1(i,GEV500,XBAND, 3),i=0,7) /
     $    .51285E+00,   .45812E+00,   .75135E+01,  -.62247E+00, 
     $    .30444E+00,  -.68530E+00,   .85519E+01,   .43062E+00 /
      data xa1lum(TEV1,  XBAND, 3) /   .11799E+03 /
      data (xa1(i,TEV1,  XBAND, 3),i=0,7) /
     $    .31241E+00,   .61241E+00,   .29938E+01,  -.55848E+00, 
     $    .44801E+00,  -.67116E+00,   .41119E+01,   .72753E+00 /
      data xa1lum(GEV350,SBAND,3) / -1d0 /
      data xa1lum(GEV350,XBAND,3) / -1d0 /
      data xa1lum(GEV800,SBAND,3) / -1d0 /
      data xa1lum(GEV800,XBAND,3) / -1d0 /
      data (xa1lum(GEV350,i,3),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV500,i,3),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV800,i,3),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV1,i,3),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV16,i,3),i=1,NACC) / 6*-1d0 /
      data xa1lum(GEV500,SBAND, 4) /   .31528E+02 /
      data (xa1(i,GEV500,SBAND, 4),i=0,7) /
     $    .38169E+00,   .73949E+00,   .12543E+02,  -.61112E+00, 
     $    .51256E+00,  -.69009E+00,   .14892E+02,   .63314E+00 /
      data xa1lum(TEV1,  SBAND, 4) /   .24613E+03 /
      data (xa1(i,TEV1,  SBAND, 4),i=0,7) /
     $    .24256E+00,   .94117E+00,   .66775E+01,  -.55160E+00, 
     $    .57484E+00,  -.68891E+00,   .92271E+01,   .81162E+00 /
      data xa1lum(GEV350,TESLA, 4) /   .74549E+02 /
      data (xa1(i,GEV350,TESLA, 4),i=0,7) /
     $    .34120E+00,   .12230E+01,   .32932E+02,  -.59850E+00, 
     $    .65947E+00,  -.69574E+00,   .38116E+02,   .63879E+00 /
      data xa1lum(GEV500,TESLA, 4) /   .10668E+03 /
      data (xa1(i,GEV500,TESLA, 4),i=0,7) /
     $    .28082E+00,   .11074E+01,   .18399E+02,  -.59118E+00, 
     $    .68880E+00,  -.69375E+00,   .23463E+02,   .76073E+00 /
      data xa1lum(GEV800,TESLA, 4) /   .29006E+03 /
      data (xa1(i,GEV800,TESLA, 4),i=0,7) /
     $    .21272E+00,   .11443E+01,   .92564E+01,  -.54657E+00, 
     $    .66799E+00,  -.69137E+00,   .12498E+02,   .87571E+00 /
      data xa1lum(TEV1,  TESLA, 4) /   .11009E+03 /
      data (xa1(i,TEV1,  TESLA, 4),i=0,7) /
     $    .41058E+00,   .64745E+00,   .11271E+02,  -.61996E+00, 
     $    .39801E+00,  -.69150E+00,   .14560E+02,   .49924E+00 /
      data xa1lum(GEV500,XBAND, 4) /   .36179E+02 /
      data (xa1(i,GEV500,XBAND, 4),i=0,7) /
     $    .51155E+00,   .43313E+00,   .70446E+01,  -.63003E+00, 
     $    .29449E+00,  -.68747E+00,   .83489E+01,   .42458E+00 /
      data xa1lum(TEV1,  XBAND, 4) /   .11748E+03 /
      data (xa1(i,TEV1,  XBAND, 4),i=0,7) /
     $    .32917E+00,   .54322E+00,   .28493E+01,  -.57959E+00, 
     $    .39266E+00,  -.68217E+00,   .38475E+01,   .68478E+00 /
      data xa1lum(GEV350,SBAND,4) / -1d0 /
      data xa1lum(GEV350,XBAND,4) / -1d0 /
      data xa1lum(GEV800,SBAND,4) / -1d0 /
      data xa1lum(GEV800,XBAND,4) / -1d0 /
      data (xa1lum(GEV350,i,4),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV500,i,4),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV800,i,4),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV1,i,4),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV16,i,4),i=1,NACC) / 6*-1d0 /
      data xa1lum(GEV350,SBAND, 5) /  0.21897E+02 /
      data (xa1(i,GEV350,SBAND, 5),i=0,7) /
     $   0.57183E+00,  0.53877E+00,  0.19422E+02, -0.63064E+00, 
     $   0.49112E+00, -0.69109E+00,  0.24331E+02,  0.52718E+00 /
      data xa1lum(GEV500,SBAND, 5) /  0.31383E+02 /
      data (xa1(i,GEV500,SBAND, 5),i=0,7) /
     $   0.51882E+00,  0.49915E+00,  0.11153E+02, -0.63017E+00, 
     $   0.50217E+00, -0.69113E+00,  0.14935E+02,  0.62373E+00 /
      data xa1lum(GEV800,SBAND, 5) /  0.95091E+02 /
      data (xa1(i,GEV800,SBAND, 5),i=0,7) /
     $   0.47137E+00,  0.46150E+00,  0.56562E+01, -0.61758E+00, 
     $   0.46863E+00, -0.68897E+00,  0.85876E+01,  0.67577E+00 /
      data xa1lum(TEV1,SBAND, 5) /  0.11900E+03 /
      data (xa1(i,TEV1,SBAND, 5),i=0,7) /
     $   0.43956E+00,  0.45471E+00,  0.42170E+01, -0.61180E+00, 
     $   0.48711E+00, -0.68696E+00,  0.67145E+01,  0.74551E+00 /
      data xa1lum(TEV16,SBAND, 5) /  0.11900E+03 /
      data (xa1(i,TEV16,SBAND, 5),i=0,7) /
     $   0.43956E+00,  0.45471E+00,  0.42170E+01, -0.61180E+00, 
     $   0.48711E+00, -0.68696E+00,  0.67145E+01,  0.74551E+00 /
      data xa1lum(GEV350,TESLA, 5) /  0.97452E+02 /
      data (xa1(i,GEV350,TESLA, 5),i=0,7) /
     $   0.39071E+00,  0.84996E+00,  0.17614E+02, -0.60609E+00, 
     $   0.73920E+00, -0.69490E+00,  0.28940E+02,  0.77286E+00 /
      data xa1lum(GEV500,TESLA, 5) /  0.10625E+03 /
      data (xa1(i,GEV500,TESLA, 5),i=0,7) /
     $   0.42770E+00,  0.71457E+00,  0.15284E+02, -0.61664E+00, 
     $   0.68166E+00, -0.69208E+00,  0.24165E+02,  0.73806E+00 /
      data xa1lum(GEV800,TESLA, 5) /  0.17086E+03 /
      data (xa1(i,GEV800,TESLA, 5),i=0,7) /
     $   0.36025E+00,  0.69118E+00,  0.76221E+01, -0.59440E+00, 
     $   0.71269E+00, -0.69077E+00,  0.13117E+02,  0.91780E+00 /
      data xa1lum(TEV1,TESLA, 5) /  0.21433E+03 /
      data (xa1(i,TEV1,TESLA, 5),i=0,7) /
     $   0.33145E+00,  0.67075E+00,  0.55438E+01, -0.58468E+00, 
     $   0.72503E+00, -0.69084E+00,  0.99992E+01,  0.10112E+01 /
      data xa1lum(TEV16,TESLA, 5) /  0.34086E+03 /
      data (xa1(i,TEV16,TESLA, 5),i=0,7) /
     $   0.49058E+00,  0.42609E+00,  0.50550E+01, -0.61867E+00, 
     $   0.39225E+00, -0.68916E+00,  0.75514E+01,  0.58754E+00 /
      data xa1lum(GEV350,XBAND, 5) /  0.31901E+02 /
      data (xa1(i,GEV350,XBAND, 5),i=0,7) /
     $   0.65349E+00,  0.31752E+00,  0.94342E+01, -0.64291E+00, 
     $   0.30364E+00, -0.68989E+00,  0.11446E+02,  0.40486E+00 /
      data xa1lum(GEV500,XBAND, 5) /  0.36386E+02 /
      data (xa1(i,GEV500,XBAND, 5),i=0,7) /
     $   0.65132E+00,  0.28728E+00,  0.69853E+01, -0.64440E+00, 
     $   0.28736E+00, -0.68758E+00,  0.83227E+01,  0.41492E+00 /
      data xa1lum(GEV800,XBAND, 5) /  0.10854E+03 /
      data (xa1(i,GEV800,XBAND, 5),i=0,7) /
     $   0.49478E+00,  0.36221E+00,  0.30116E+01, -0.61548E+00, 
     $   0.39890E+00, -0.68418E+00,  0.45183E+01,  0.67243E+00 /
      data xa1lum(TEV1,XBAND, 5) /  0.11899E+03 /
      data (xa1(i,TEV1,XBAND, 5),i=0,7) /
     $   0.49992E+00,  0.34299E+00,  0.26184E+01, -0.61584E+00, 
     $   0.38450E+00, -0.68342E+00,  0.38589E+01,  0.67408E+00 /
      data xa1lum(TEV16,XBAND, 5) /  0.13675E+03 /
      data (xa1(i,TEV16,XBAND, 5),i=0,7) /
     $   0.50580E+00,  0.30760E+00,  0.18339E+01, -0.61421E+00, 
     $   0.35233E+00, -0.68315E+00,  0.26708E+01,  0.67918E+00 /
      data xa1lum(GEV500,SBNDEE, 0) /   .92914E+01 /
      data (xa1(i,GEV500,SBNDEE, 0),i=0,7) /
     $    .34866E+00,   .78710E+00,   .10304E+02,  -.59464E+00, 
     $    .40234E+00,  -.69741E+00,   .20645E+02,   .47274E+00 /
      data xa1lum(TEV1,  SBNDEE, 0) /   .45586E+02 /
      data (xa1(i,TEV1,  SBNDEE, 0),i=0,7) /
     $    .21084E+00,   .99168E+00,   .54407E+01,  -.52851E+00, 
     $    .47493E+00,  -.69595E+00,   .12480E+02,   .64027E+00 /
      data xa1lum(GEV350,TESLEE, 0) /   .15175E+02 /
      data (xa1(i,GEV350,TESLEE, 0),i=0,7) /
     $    .33093E+00,   .11137E+01,   .25275E+02,  -.59942E+00, 
     $    .49623E+00,  -.70403E+00,   .60188E+02,   .44637E+00 /
      data xa1lum(GEV500,TESLEE, 0) /   .21622E+02 /
      data (xa1(i,GEV500,TESLEE, 0),i=0,7) /
     $    .27175E+00,   .10697E+01,   .14858E+02,  -.58418E+00, 
     $    .50824E+00,  -.70387E+00,   .36129E+02,   .53002E+00 /
      data xa1lum(GEV800,TESLEE, 0) /   .43979E+02 /
      data (xa1(i,GEV800,TESLEE, 0),i=0,7) /
     $    .22994E+00,   .10129E+01,   .81905E+01,  -.55751E+00, 
     $    .46551E+00,  -.70461E+00,   .19394E+02,   .58387E+00 /
      data xa1lum(TEV1,  TESLEE, 0) /   .25465E+02 /
      data (xa1(i,TEV1,  TESLEE, 0),i=0,7) /
     $    .37294E+00,   .67522E+00,   .87504E+01,  -.60576E+00, 
     $    .35095E+00,  -.69821E+00,   .18526E+02,   .42784E+00 /
      data xa1lum(GEV500,XBNDEE, 0) /   .13970E+02 /
      data (xa1(i,GEV500,XBNDEE, 0),i=0,7) /
     $    .47296E+00,   .46800E+00,   .58897E+01,  -.61689E+00, 
     $    .27181E+00,  -.68923E+00,   .10087E+02,   .37462E+00 /
      data xa1lum(TEV1,  XBNDEE, 0) /   .41056E+02 /
      data (xa1(i,TEV1,  XBNDEE, 0),i=0,7) /
     $    .27965E+00,   .74816E+00,   .27415E+01,  -.50491E+00, 
     $    .38320E+00,  -.67945E+00,   .47506E+01,   .62218E+00 /
      data xa1lum(GEV350,SBNDEE,0) / -1d0 /
      data xa1lum(GEV350,XBNDEE,0) / -1d0 /
      data xa1lum(GEV800,SBNDEE,0) / -1d0 /
      data xa1lum(GEV800,XBNDEE,0) / -1d0 /
      data xa1lum(GEV500,SBAND, 0) /   .31528E+02 /
      data (xa1(i,GEV500,SBAND, 0),i=0,7) /
     $    .38169E+00,   .73949E+00,   .12543E+02,  -.61112E+00, 
     $    .51256E+00,  -.69009E+00,   .14892E+02,   .63314E+00 /
      data xa1lum(TEV1,  SBAND, 0) /   .24613E+03 /
      data (xa1(i,TEV1,  SBAND, 0),i=0,7) /
     $    .24256E+00,   .94117E+00,   .66775E+01,  -.55160E+00, 
     $    .57484E+00,  -.68891E+00,   .92271E+01,   .81162E+00 /
      data xa1lum(GEV350,TESLA, 0) /   .74549E+02 /
      data (xa1(i,GEV350,TESLA, 0),i=0,7) /
     $    .34120E+00,   .12230E+01,   .32932E+02,  -.59850E+00, 
     $    .65947E+00,  -.69574E+00,   .38116E+02,   .63879E+00 /
      data xa1lum(GEV500,TESLA, 0) /   .10668E+03 /
      data (xa1(i,GEV500,TESLA, 0),i=0,7) /
     $    .28082E+00,   .11074E+01,   .18399E+02,  -.59118E+00, 
     $    .68880E+00,  -.69375E+00,   .23463E+02,   .76073E+00 /
      data xa1lum(GEV800,TESLA, 0) /   .29006E+03 /
      data (xa1(i,GEV800,TESLA, 0),i=0,7) /
     $    .21272E+00,   .11443E+01,   .92564E+01,  -.54657E+00, 
     $    .66799E+00,  -.69137E+00,   .12498E+02,   .87571E+00 /
      data xa1lum(TEV1,  TESLA, 0) /   .11009E+03 /
      data (xa1(i,TEV1,  TESLA, 0),i=0,7) /
     $    .41058E+00,   .64745E+00,   .11271E+02,  -.61996E+00, 
     $    .39801E+00,  -.69150E+00,   .14560E+02,   .49924E+00 /
      data xa1lum(GEV500,XBAND, 0) /   .36179E+02 /
      data (xa1(i,GEV500,XBAND, 0),i=0,7) /
     $    .51155E+00,   .43313E+00,   .70446E+01,  -.63003E+00, 
     $    .29449E+00,  -.68747E+00,   .83489E+01,   .42458E+00 /
      data xa1lum(TEV1,  XBAND, 0) /   .11748E+03 /
      data (xa1(i,TEV1,  XBAND, 0),i=0,7) /
     $    .32917E+00,   .54322E+00,   .28493E+01,  -.57959E+00, 
     $    .39266E+00,  -.68217E+00,   .38475E+01,   .68478E+00 /
      data xa1lum(GEV350,SBAND,0) / -1d0 /
      data xa1lum(GEV350,XBAND,0) / -1d0 /
      data xa1lum(GEV800,SBAND,0) / -1d0 /
      data xa1lum(GEV800,XBAND,0) / -1d0 /
      data xa3lum(GEV800,TESLA, 3) /   .17196E+03 /
      data (xa3(i,GEV800,TESLA, 3),i=0,7) /
     $    .21633E+00,   .11333E+01,   .95928E+01,  -.55095E+00, 
     $    .73044E+00,  -.69101E+00,   .12868E+02,   .94737E+00 /
      data xa3lum(GEV800,TESLA, 4) /   .16408E+03 /
      data (xa3(i,GEV800,TESLA, 4),i=0,7) /
     $    .41828E+00,   .72418E+00,   .14137E+02,  -.61189E+00, 
     $    .36697E+00,  -.69205E+00,   .17713E+02,   .43583E+00 /
      data xa3lum(GEV350,TESLA, 5) /  0.66447E+02 /
      data (xa3(i,GEV350,TESLA, 5),i=0,7) /
     $   0.69418E+00,  0.50553E+00,  0.48430E+02, -0.63911E+00, 
     $   0.34074E+00, -0.69533E+00,  0.55502E+02,  0.29397E+00 /
      data xa3lum(GEV500,TESLA, 5) /  0.95241E+02 /
      data (xa3(i,GEV500,TESLA, 5),i=0,7) /
     $   0.64882E+00,  0.45462E+00,  0.27103E+02, -0.64535E+00, 
     $   0.35101E+00, -0.69467E+00,  0.33658E+02,  0.35024E+00 /
      data xa3lum(GEV800,TESLA, 5) /  0.16974E+03 /
      data (xa3(i,GEV800,TESLA, 5),i=0,7) /
     $   0.58706E+00,  0.43771E+00,  0.13422E+02, -0.63804E+00, 
     $   0.35541E+00, -0.69467E+00,  0.17528E+02,  0.43051E+00 /
      data xa3lum(TEV1,TESLA, 5) /  0.21222E+03 /
      data (xa3(i,TEV1,TESLA, 5),i=0,7) /
     $   0.55525E+00,  0.42577E+00,  0.96341E+01, -0.63587E+00, 
     $   0.36448E+00, -0.69365E+00,  0.13161E+02,  0.47715E+00 /
      data xa3lum(TEV16,TESLA, 5) /  0.34086E+03 /
      data (xa3(i,TEV16,TESLA, 5),i=0,7) /
     $   0.49058E+00,  0.42609E+00,  0.50550E+01, -0.61867E+00, 
     $   0.39225E+00, -0.68916E+00,  0.75514E+01,  0.58754E+00 /
      data xa3lum(GEV350,TESLA, 0) /  0.66447E+02 /
      data (xa3(i,GEV350,TESLA, 0),i=0,7) /
     $   0.69418E+00,  0.50553E+00,  0.48430E+02, -0.63911E+00, 
     $   0.34074E+00, -0.69533E+00,  0.55502E+02,  0.29397E+00 /
      data xa3lum(GEV500,TESLA, 0) /  0.95241E+02 /
      data (xa3(i,GEV500,TESLA, 0),i=0,7) /
     $   0.64882E+00,  0.45462E+00,  0.27103E+02, -0.64535E+00, 
     $   0.35101E+00, -0.69467E+00,  0.33658E+02,  0.35024E+00 /
      data xa3lum(GEV800,TESLA, 0) /  0.16974E+03 /
      data (xa3(i,GEV800,TESLA, 0),i=0,7) /
     $   0.58706E+00,  0.43771E+00,  0.13422E+02, -0.63804E+00, 
     $   0.35541E+00, -0.69467E+00,  0.17528E+02,  0.43051E+00 /
      data xa3lum(TEV1,TESLA, 0) /  0.21222E+03 /
      data (xa3(i,TEV1,TESLA, 0),i=0,7) /
     $   0.55525E+00,  0.42577E+00,  0.96341E+01, -0.63587E+00, 
     $   0.36448E+00, -0.69365E+00,  0.13161E+02,  0.47715E+00 /
      data xa3lum(TEV16,TESLA, 0) /  0.34086E+03 /
      data (xa3(i,TEV16,TESLA, 0),i=0,7) /
     $   0.49058E+00,  0.42609E+00,  0.50550E+01, -0.61867E+00, 
     $   0.39225E+00, -0.68916E+00,  0.75514E+01,  0.58754E+00 /
      data xa5lum(GEV350,TESLA, 1) /  -1.0 /
      data xa5lum(GEV500,TESLA, 1) /  0.33980E+03 /
      data (xa5(i,GEV500,TESLA, 1),i=0,7) /
     $   0.49808E+00,  0.54613E+00,  0.12287E+02, -0.62756E+00, 
     $   0.42817E+00, -0.69120E+00,  0.17067E+02,  0.51143E+00 /
      data xa5lum(GEV800,TESLA, 1) /  0.35936E+03 /
      data (xa5(i,GEV800,TESLA, 1),i=0,7) /
     $   0.58751E+00,  0.43128E+00,  0.13324E+02, -0.64006E+00, 
     $   0.30682E+00, -0.69235E+00,  0.16815E+02,  0.37078E+00 /
      data xa5lum(TEV1,  TESLA, 1) /  -1.0 /
      data xa5lum(TEV16, TESLA, 1) /  -1.0 /
      data xa5lum(GEV350,TESLA, 0) /  -1.0 /
      data xa5lum(GEV500,TESLA, 0) /  0.33980E+03 /
      data (xa5(i,GEV500,TESLA, 0),i=0,7) /
     $   0.49808E+00,  0.54613E+00,  0.12287E+02, -0.62756E+00, 
     $   0.42817E+00, -0.69120E+00,  0.17067E+02,  0.51143E+00 /
      data xa5lum(GEV800,TESLA, 0) /  0.35936E+03 /
      data (xa5(i,GEV800,TESLA, 0),i=0,7) /
     $   0.58751E+00,  0.43128E+00,  0.13324E+02, -0.64006E+00, 
     $   0.30682E+00, -0.69235E+00,  0.16815E+02,  0.37078E+00 /
      data xa5lum(TEV1,  TESLA, 0) /  -1.0 /
      data xa5lum(TEV16, TESLA, 0) /  -1.0 /
      data xa6lum(GEV090,TESLA, 1) /  0.62408E+02 /
      data (xa6(i,GEV090,TESLA, 1),i=0,7) /
     $   0.72637E+00,  0.75534E+00,  0.18180E+03, -0.63426E+00,
     $   0.36829E+00, -0.69653E+00,  0.18908E+03,  0.22157E+00 /
      data xa6lum(GEV170,TESLA, 1) /  0.11532E+02 /
      data (xa6(i,GEV170,TESLA, 1),i=0,7) /
     $   0.65232E+00,  0.67249E+00,  0.66862E+02, -0.63315E+00,
     $   0.38470E+00, -0.69477E+00,  0.75120E+02,  0.30162E+00 /
      data xa6lum(GEV350,TESLA, 1) /  0.24641E+03 /
      data (xa6(i,GEV350,TESLA, 1),i=0,7) /
     $   0.54610E+00,  0.59105E+00,  0.20297E+02, -0.62747E+00,
     $   0.41588E+00, -0.69188E+00,  0.26345E+02,  0.43818E+00 /
      data xa6lum(GEV500,TESLA, 1) /  0.30340E+03 /
      data (xa6(i,GEV500,TESLA, 1),i=0,7) /
     $   0.52744E+00,  0.52573E+00,  0.13895E+02, -0.63145E+00,
     $   0.40824E+00, -0.69150E+00,  0.18645E+02,  0.47585E+00 /
      data xa6lum(GEV090,TESLA, 0) /  0.62408E+02 /
      data (xa6(i,GEV090,TESLA, 0),i=0,7) /
     $   0.72637E+00,  0.75534E+00,  0.18180E+03, -0.63426E+00,
     $   0.36829E+00, -0.69653E+00,  0.18908E+03,  0.22157E+00 /
      data xa6lum(GEV170,TESLA, 0) /  0.11532E+02 /
      data (xa6(i,GEV170,TESLA, 0),i=0,7) /
     $   0.65232E+00,  0.67249E+00,  0.66862E+02, -0.63315E+00,
     $   0.38470E+00, -0.69477E+00,  0.75120E+02,  0.30162E+00 /
      data xa6lum(GEV350,TESLA, 0) /  0.24641E+03 /
      data (xa6(i,GEV350,TESLA, 0),i=0,7) /
     $   0.54610E+00,  0.59105E+00,  0.20297E+02, -0.62747E+00,
     $   0.41588E+00, -0.69188E+00,  0.26345E+02,  0.43818E+00 /
      data xa6lum(GEV500,TESLA, 0) /  0.30340E+03 /
      data (xa6(i,GEV500,TESLA, 0),i=0,7) /
     $   0.52744E+00,  0.52573E+00,  0.13895E+02, -0.63145E+00,
     $   0.40824E+00, -0.69150E+00,  0.18645E+02,  0.47585E+00 /
      data xa7lum(GEV090,TESLA,1) /  0.62408E+02 /
      data (xa7(i,GEV090,TESLA,1),i=0,7) /
     $   0.72637E+00,  0.75534E+00,  0.18180E+03, -0.63426E+00,
     $   0.36829E+00, -0.69653E+00,  0.18908E+03,  0.22157E+00 /
      data xa7lum(GEV170,TESLA,1) /  0.11532E+02 /
      data (xa7(i,GEV170,TESLA,1),i=0,7) /
     $   0.65232E+00,  0.67249E+00,  0.66862E+02, -0.63315E+00,
     $   0.38470E+00, -0.69477E+00,  0.75120E+02,  0.30162E+00 /
      data xa7lum(GEV350,TESLA,1) /  0.24641E+03 /
      data (xa7(i,GEV350,TESLA,1),i=0,7) /
     $   0.54610E+00,  0.59105E+00,  0.20297E+02, -0.62747E+00,
     $   0.41588E+00, -0.69188E+00,  0.26345E+02,  0.43818E+00 /
      data xa7lum(GEV500,TESLA,1) /  0.34704E+03 /
      data (xa7(i,GEV500,TESLA,1),i=0,7) /
     $   0.51288E+00,  0.49025E+00,  0.99716E+01, -0.62850E+00, 
     $   0.41048E+00, -0.69065E+00,  0.13922E+02,  0.51902E+00 /
      data xa7lum(GEV800,TESLA,1) /  0.57719E+03 /
      data (xa7(i,GEV800,TESLA,1),i=0,7) /
     $   0.52490E+00,  0.42573E+00,  0.69069E+01, -0.62649E+00, 
     $   0.32380E+00, -0.68958E+00,  0.93819E+01,  0.45671E+00 /
      data xa7lum(TEV1,  TESLA,1) /  -1.0 /
      data xa7lum(GEV090,JLCNLC,1) /  -1.0 /
      data xa7lum(GEV170,JLCNLC,1) /  -1.0 /
      data xa7lum(GEV350,JLCNLC,1) /  -1.0 /
      data xa7lum(GEV500,JLCNLC,1) /  0.63039E+02 /
      data (xa7(i,GEV500,JLCNLC,1),i=0,7) /
     $   0.58967E+00,  0.34035E+00,  0.63631E+01, -0.63683E+00, 
     $   0.33383E+00, -0.68803E+00,  0.81005E+01,  0.48702E+00 /
      data xa7lum(TEV1,JLCNLC,1) /  0.12812E+03 /
      data (xa7(i,TEV1,JLCNLC,1),i=0,7) /
     $   0.50222E+00,  0.33773E+00,  0.25681E+01, -0.61711E+00, 
     $   0.36826E+00, -0.68335E+00,  0.36746E+01,  0.65393E+00 /
      data xa7lum(GEV090,TESLA,0) /  0.62408E+02 /
      data (xa7(i,GEV090,TESLA,0),i=0,7) /
     $   0.72637E+00,  0.75534E+00,  0.18180E+03, -0.63426E+00,
     $   0.36829E+00, -0.69653E+00,  0.18908E+03,  0.22157E+00 /
      data xa7lum(GEV170,TESLA,0) /  0.11532E+02 /
      data (xa7(i,GEV170,TESLA,0),i=0,7) /
     $   0.65232E+00,  0.67249E+00,  0.66862E+02, -0.63315E+00,
     $   0.38470E+00, -0.69477E+00,  0.75120E+02,  0.30162E+00 /
      data xa7lum(GEV350,TESLA,0) /  0.24641E+03 /
      data (xa7(i,GEV350,TESLA,0),i=0,7) /
     $   0.54610E+00,  0.59105E+00,  0.20297E+02, -0.62747E+00,
     $   0.41588E+00, -0.69188E+00,  0.26345E+02,  0.43818E+00 /
      data xa7lum(GEV500,TESLA,0) /  0.34704E+03 /
      data (xa7(i,GEV500,TESLA,0),i=0,7) /
     $   0.51288E+00,  0.49025E+00,  0.99716E+01, -0.62850E+00, 
     $   0.41048E+00, -0.69065E+00,  0.13922E+02,  0.51902E+00 /
      data xa7lum(GEV800,TESLA,0) /  0.57719E+03 /
      data (xa7(i,GEV800,TESLA,0),i=0,7) /
     $   0.52490E+00,  0.42573E+00,  0.69069E+01, -0.62649E+00, 
     $   0.32380E+00, -0.68958E+00,  0.93819E+01,  0.45671E+00 /
      data xa7lum(TEV1,  TESLA,0) /  -1.0 /
      data xa7lum(GEV090,JLCNLC,0) /  -1.0 /
      data xa7lum(GEV170,JLCNLC,0) /  -1.0 /
      data xa7lum(GEV350,JLCNLC,0) /  -1.0 /
      data xa7lum(GEV500,JLCNLC,0) /  0.63039E+02 /
      data (xa7(i,GEV500,JLCNLC,0),i=0,7) /
     $   0.58967E+00,  0.34035E+00,  0.63631E+01, -0.63683E+00, 
     $   0.33383E+00, -0.68803E+00,  0.81005E+01,  0.48702E+00 /
      data xa7lum(TEV1,JLCNLC,0) /  0.12812E+03 /
      data (xa7(i,TEV1,JLCNLC,0),i=0,7) /
     $   0.50222E+00,  0.33773E+00,  0.25681E+01, -0.61711E+00, 
     $   0.36826E+00, -0.68335E+00,  0.36746E+01,  0.65393E+00 /
      if (magic .ne. 1904 06 16) then
         magic = 1904 06 16
      x1m = 0d0
      x2m = 0d0
      roots = 500D0
      acc = TESLA
      ver = 0
      rev = 0
      chat = 1
      if (xchat .ne. 0) then
         call circem ('MESSAGE', 'starting up ...')
         call circem ('MESSAGE',
     $        '$Id: circe.f,v 1.3 2004/09/20 13:53:03 kilian Exp $')
      endif
      endif
      if ((xchat .ge. 0) .and. (xchat .ne. chat)) then
         chat = xchat
         if (chat .ge. 1) then
            write (msgbuf, 1000) 'chat', chat
 1000       format ('updating `', A, ''' to ', I2)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1100) 'chat', chat
 1100       format ('keeping `', A, ''' at ', I2)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xx1m .ge. 0d0) .and. (xx1m .ne. x1m)) then
         x1m = xx1m
         if (chat .ge. 1) then
            write (msgbuf, 1001) 'x1min', x1m
 1001       format ('updating `', A, ''' to ', E12.4)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1101) 'x1min', x1m
 1101       format ('keeping `', A, ''' at ', E12.4)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xx2m .ge. 0d0) .and. (xx2m .ne. x2m)) then
         x2m = xx2m
         if (chat .ge. 1) then
            write (msgbuf, 1001) 'x2min', x2m
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1101) 'x2min', x2m
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xroots .ge. 0d0) .and.(xroots .ne. roots)) then
         roots = xroots
         if (chat .ge. 1) then
            write (msgbuf, 1002) 'roots', roots
 1002       format ('updating `', A, ''' to ', F6.1)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1102) 'roots', roots
 1102       format ('keeping `', A, ''' at ', F6.1)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xacc .ge. 0) .and.(xacc .ne. acc)) then
         if ((xacc .ge. 1) .and. (xacc .le. NACC)) then
            acc = xacc
            if (chat .ge. 1) then
               write (msgbuf, 1003) 'acc', accnam(acc)
 1003          format ('updating `', A, ''' to ', A)
               call circem ('MESSAGE', msgbuf)
            endif
         else
            write (msgbuf, 1203) xacc
 1203       format ('invalid `acc'': ', I8)
            call circem ('ERROR', msgbuf)
            write (msgbuf, 1103) 'acc', accnam(acc)
 1103       format ('keeping `', A, ''' at ', A)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1003) 'acc', accnam(acc)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((acc .eq. SBNDEE) .or. (acc .eq. TESLEE)
     $    .or. (acc .eq. XBNDEE)) then
      call circem ('WARNING', '***********************************')
      call circem ('WARNING', '* The accelerator parameters have *')
      call circem ('WARNING', '* not been endorsed for use in    *')
      call circem ('WARNING', '* an e-e- collider yet!!!         *')
      call circem ('WARNING', '***********************************')
      endif
      if (xver .ge. 0) then
         ver = xver
         if (chat .ge. 1) then
            write (msgbuf, 1000) 'ver', ver
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1100) 'ver', ver
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xrev .ge. 0) .and.(xrev .ne. rev)) then
         rev = xrev
         if (chat .ge. 1) then
            write (msgbuf, 1004) 'rev', rev
 1004       format ('updating `', A, ''' to ', I8)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1104) 'rev', rev
 1104       format ('keeping `', A, ''' at ', I8)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      ver34 = 0
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      if (rev .eq. 0) then
         r = 0
c       call circem ('WARNING', '*************************************')
c       call circem ('WARNING', '* This release is not official yet, *')
c       call circem ('WARNING', '* do not use it in publications!    *')
c       call circem ('WARNING', '*************************************')
      elseif (rev .ge. 1997 04 17) then
         r = 5
      elseif (rev .ge. 1996 09 02) then
         r = 4
      elseif (rev .ge. 1996 07 29) then
         r = 3
      elseif (rev .ge. 1996 07 11) then
         r = 2
      elseif (rev .ge. 1996 04 01) then
         r = 1
      elseif (rev .lt. 1996 04 01) then
         call circem ('ERROR',
     $        'no revision of version 1 before 96/04/01 available')
         call circem ('MESSAGE', 'falling back to default')
         r = 1
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2000) rev, r
 2000    format ('mapping date ', I8, ' to revision index ', I2)
         call circem ('MESSAGE', msgbuf)
      endif
      if (roots .eq. 350d0) then
         e = GEV350
      elseif ((roots .ge. 340d0) .and. (roots .le. 370d0)) then
         write (msgbuf, 2001) roots, 350d0
         call circem ('MESSAGE', msgbuf)
         e = GEV350
      elseif (roots .eq. 500d0) then
         e = GEV500
      elseif ((roots .ge. 480d0) .and. (roots .le. 520d0)) then
         write (msgbuf, 2001) roots, 500d0
         call circem ('MESSAGE', msgbuf)
         e = GEV500
      elseif (roots .eq. 800d0) then
         e = GEV800
      elseif ((roots .ge. 750d0) .and. (roots .le. 850d0)) then
         write (msgbuf, 2001) roots, 800d0
         call circem ('MESSAGE', msgbuf)
         e = GEV800
      elseif (roots .eq. 1000d0) then
         e = TEV1
      elseif ((roots .ge. 900d0) .and. (roots .le. 1100d0)) then
         write (msgbuf, 2001) roots, 1000d0
         call circem ('MESSAGE', msgbuf)
         e = TEV1
      elseif (roots .eq. 1600d0) then
         e = TEV16
      elseif ((roots .ge. 1500d0) .and. (roots .le. 1700d0)) then
         write (msgbuf, 2001) roots, 1600d0
         call circem ('MESSAGE', msgbuf)
         e = TEV16
      else
         call circem ('ERROR',
     $        'only ROOTS = 350, 500, 800, 1000 and 1600GeV available')
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (xa1lum(e,acc,r) .lt. 0d0) then
         write (msgbuf, 2002) roots, accnam(acc), r
         call circem ('ERROR', msgbuf)
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (chat .ge. 2) then
         if (e .ge. GEV090) then
            write (msgbuf, 2003) roots, e
            call circem ('MESSAGE', msgbuf)
         else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
            write (msgbuf, 2013) roots, elo, ehi
            call circem ('MESSAGE', msgbuf)
         end if
      endif
      lumi = xa1lum (e,acc,r)
      do 10 i = 0, 7
         a1(i) = xa1(i,e,acc,r)
 10   continue
      elseif ((ver .eq. 3) .or. (ver .eq. 4)) then
         ver34 = ver
         ver = 1
      if (rev .eq. 0) then
         r = 0
      call circem ('WARNING', '*************************************')
      call circem ('WARNING', '* This release is not official yet, *')
      call circem ('WARNING', '* do not use it in publications!    *')
      call circem ('WARNING', '*************************************')
      elseif (rev .ge. 1997 04 17) then
         r = 5
         if (ver34 .eq. 3) then
            call circem ('WARNING', 'version 3 retired after 97/04/17')
            call circem ('MESSAGE', 'falling back to version 4')
         endif
      elseif (rev .ge. 1996 10 22) then
         r = ver34
         if ((roots .ne. 800d0) .or. (acc .ne. TESLA)) then
            call circem ('ERROR', 'versions 3 and 4 before 97/04/17')
            call circem ('ERROR', 'apply to TESLA at 800 GeV only')
            call circem ('MESSAGE', 'falling back to TESLA at 800GeV')
            acc = TESLA
            e = GEV800
         endif
      elseif (rev .lt. 1996 10 22) then
         call circem ('ERROR',
     $     'no revision of versions 3 and 4 available before 96/10/22')
         call circem ('MESSAGE', 'falling back to default')
         r = 5
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2000) rev, r
         call circem ('MESSAGE', msgbuf)
      endif
      if (roots .eq. 350d0) then
         e = GEV350
      elseif ((roots .ge. 340d0) .and. (roots .le. 370d0)) then
         write (msgbuf, 2001) roots, 350d0
         call circem ('MESSAGE', msgbuf)
         e = GEV350
      elseif (roots .eq. 500d0) then
         e = GEV500
      elseif ((roots .ge. 480d0) .and. (roots .le. 520d0)) then
         write (msgbuf, 2001) roots, 500d0
         call circem ('MESSAGE', msgbuf)
         e = GEV500
      elseif (roots .eq. 800d0) then
         e = GEV800
      elseif ((roots .ge. 750d0) .and. (roots .le. 850d0)) then
         write (msgbuf, 2001) roots, 800d0
         call circem ('MESSAGE', msgbuf)
         e = GEV800
      elseif (roots .eq. 1000d0) then
         e = TEV1
      elseif ((roots .ge. 900d0) .and. (roots .le. 1100d0)) then
         write (msgbuf, 2001) roots, 1000d0
         call circem ('MESSAGE', msgbuf)
         e = TEV1
      elseif (roots .eq. 1600d0) then
         e = TEV16
      elseif ((roots .ge. 1500d0) .and. (roots .le. 1700d0)) then
         write (msgbuf, 2001) roots, 1600d0
         call circem ('MESSAGE', msgbuf)
         e = TEV16
      else
         call circem ('ERROR',
     $        'only ROOTS = 350, 500, 800, 1000 and 1600GeV available')
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (xa3lum(e,acc,r) .lt. 0d0) then
         write (msgbuf, 2002) roots, accnam(acc), r
         call circem ('ERROR', msgbuf)
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (chat .ge. 2) then
         if (e .ge. GEV090) then
            write (msgbuf, 2003) roots, e
            call circem ('MESSAGE', msgbuf)
         else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
            write (msgbuf, 2013) roots, elo, ehi
            call circem ('MESSAGE', msgbuf)
         end if
      endif
      lumi = xa3lum (e,acc,r)
      do 20 i = 0, 7
         a1(i) = xa3(i,e,acc,r)
 20   continue
      elseif (ver .eq. 5) then
         ver = 1
      if (rev .eq. 0) then
         r = 0
      call circem ('WARNING', '*************************************')
      call circem ('WARNING', '* This release is not official yet, *')
      call circem ('WARNING', '* do not use it in publications!    *')
      call circem ('WARNING', '*************************************')
      elseif (rev .ge. 1998 05 05) then
         r = 1
      elseif (rev .lt. 1998 05 05) then
         call circem ('ERROR',
     $     'no revision of version 5 available before 98/05/05')
         call circem ('MESSAGE', 'falling back to default')
         r = 1
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2000) rev, r
         call circem ('MESSAGE', msgbuf)
      endif
      if (acc .ne. TESLA) then
         call circem ('ERROR', 'versions 5 applies to TESLA only')
         acc = TESLA
      end if
      if (roots .eq. 350d0) then
         e = GEV350
      elseif ((roots .ge. 340d0) .and. (roots .le. 370d0)) then
         write (msgbuf, 2001) roots, 350d0
         call circem ('MESSAGE', msgbuf)
         e = GEV350
      elseif (roots .eq. 500d0) then
         e = GEV500
      elseif ((roots .ge. 480d0) .and. (roots .le. 520d0)) then
         write (msgbuf, 2001) roots, 500d0
         call circem ('MESSAGE', msgbuf)
         e = GEV500
      elseif (roots .eq. 800d0) then
         e = GEV800
      elseif ((roots .ge. 750d0) .and. (roots .le. 850d0)) then
         write (msgbuf, 2001) roots, 800d0
         call circem ('MESSAGE', msgbuf)
         e = GEV800
      elseif (roots .eq. 1000d0) then
         e = TEV1
      elseif ((roots .ge. 900d0) .and. (roots .le. 1100d0)) then
         write (msgbuf, 2001) roots, 1000d0
         call circem ('MESSAGE', msgbuf)
         e = TEV1
      elseif (roots .eq. 1600d0) then
         e = TEV16
      elseif ((roots .ge. 1500d0) .and. (roots .le. 1700d0)) then
         write (msgbuf, 2001) roots, 1600d0
         call circem ('MESSAGE', msgbuf)
         e = TEV16
      else
         call circem ('ERROR',
     $        'only ROOTS = 350, 500, 800, 1000 and 1600GeV available')
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (xa5lum(e,acc,r) .lt. 0d0) then
         write (msgbuf, 2002) roots, accnam(acc), r
         call circem ('ERROR', msgbuf)
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (chat .ge. 2) then
         if (e .ge. GEV090) then
            write (msgbuf, 2003) roots, e
            call circem ('MESSAGE', msgbuf)
         else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
            write (msgbuf, 2013) roots, elo, ehi
            call circem ('MESSAGE', msgbuf)
         end if
      endif
      lumi = xa5lum (e,acc,r)
      do 30 i = 0, 7
         a1(i) = xa5(i,e,acc,r)
 30   continue
      elseif (ver .eq. 6) then
         ver = 1
      if (rev .eq. 0) then
         r = 0
      call circem ('WARNING', '*************************************')
      call circem ('WARNING', '* This release is not official yet, *')
      call circem ('WARNING', '* do not use it in publications!    *')
      call circem ('WARNING', '*************************************')
      elseif (rev .ge. 1999 04 15) then
         r = 1
      elseif (rev .lt. 1999 04 15) then
         call circem ('ERROR',
     $     'no revision of version 6 available before 1999/04/15')
         call circem ('MESSAGE', 'falling back to default')
         r = 1
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2000) rev, r
         call circem ('MESSAGE', msgbuf)
      endif
      if (acc .ne. TESLA) then
         call circem ('ERROR', 'versions 6 applies to TESLA only')
         acc = TESLA
      end if
      if (roots .eq.  90d0) then
         e = GEV090
      elseif ((roots .ge. 85d0) .and. (roots .le. 95d0)) then
         write (msgbuf, 2001) roots, 90d0
         call circem ('MESSAGE', msgbuf)
         e = GEV090
      elseif (roots .eq. 170d0) then
         e = GEV170
      elseif ((roots .ge. 160d0) .and. (roots .le. 180d0)) then
         write (msgbuf, 2001) roots, 170d0
         call circem ('MESSAGE', msgbuf)
         e = GEV170
      elseif (roots .eq. 350d0) then
         e = GEV350
      elseif ((roots .ge. 340d0) .and. (roots .le. 370d0)) then
         write (msgbuf, 2001) roots, 350d0
         call circem ('MESSAGE', msgbuf)
         e = GEV350
      elseif (roots .eq. 500d0) then
         e = GEV500
      elseif ((roots .ge. 480d0) .and. (roots .le. 520d0)) then
         write (msgbuf, 2001) roots, 500d0
         call circem ('MESSAGE', msgbuf)
         e = GEV500
      else
         call circem ('ERROR',
     $        'only ROOTS = 90, 170, 350, and 500GeV available')
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (xa6lum(e,acc,r) .lt. 0d0) then
         write (msgbuf, 2002) roots, accnam(acc), r
         call circem ('ERROR', msgbuf)
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (chat .ge. 2) then
         if (e .ge. GEV090) then
            write (msgbuf, 2003) roots, e
            call circem ('MESSAGE', msgbuf)
         else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
            write (msgbuf, 2013) roots, elo, ehi
            call circem ('MESSAGE', msgbuf)
         end if
      endif
      lumi = xa6lum (e,acc,r)
      do 40 i = 0, 7
         a1(i) = xa6(i,e,acc,r)
 40   continue
      elseif (ver .eq. 7) then
         ver = 1
      if (rev .eq. 0) then
         r = 0
      call circem ('WARNING', '*************************************')
      call circem ('WARNING', '* This release is not official yet, *')
      call circem ('WARNING', '* do not use it in publications!    *')
      call circem ('WARNING', '*************************************')
      elseif (rev .ge. 2000 04 26) then
         r = 1
      elseif (rev .lt. 2000 04 26) then
         call circem ('ERROR',
     $     'no revision of version 7 available before 2000/04/26')
         call circem ('MESSAGE', 'falling back to default')
         r = 1
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2000) rev, r
         call circem ('MESSAGE', msgbuf)
      endif
      if (acc .ne. TESLA .and. acc .ne. JLCNLC) then
         call circem ('ERROR',
     $                'version 7 applies to TESLA and JLCNLC only')
         call circem ('ERROR', 'falling back to TESLA')
         acc = TESLA
      end if
      e = GEV090 - 1
      elo = e
      ehi = e
      if (acc .eq. TESLA) then
         if (roots .lt.  90d0) then
            write (msgbuf, 2004) roots, 90d0
            call circem ('MESSAGE', msgbuf)
            e = GEV090
         elseif (roots .eq. 090d0) then
            e = GEV090
         elseif (roots .lt. 170d0) then
            write (msgbuf, 2005) roots, 170d0
            call circem ('MESSAGE', msgbuf)
            e = GEV170
         elseif (roots .eq. 170d0) then
            e = GEV170
         elseif (roots .lt. 350d0) then
            write (msgbuf, 2006) roots, 170d0, 350d0
            call circem ('MESSAGE', msgbuf)
            elo = GEV170
            ehi = GEV350
            eloval = 170d0
            ehival = 350d0
         elseif (roots .eq. 350d0) then
            e = GEV350
         elseif (roots .lt. 500d0) then
            write (msgbuf, 2006) roots, 350d0, 500d0
            call circem ('MESSAGE', msgbuf)
            elo = GEV350
            ehi = GEV500
            eloval = 350d0
            ehival = 500d0
         elseif (roots .eq. 500d0) then
            e = GEV500
         elseif (roots .lt. 800d0) then
            write (msgbuf, 2006) roots, 500d0, 800d0
            call circem ('MESSAGE', msgbuf)
            elo = GEV500
            ehi = GEV800
            eloval = 500d0
            ehival = 800d0
         elseif (roots .eq. 800d0) then
            e = GEV800
         else
            write (msgbuf, 2005) roots, 800d0
            call circem ('MESSAGE', msgbuf)
            e = GEV800
         endif
      elseif (acc .eq. XBAND) then
         if (roots .lt.  500d0) then
            write (msgbuf, 2004) roots, 500d0
            call circem ('MESSAGE', msgbuf)
            e = GEV500
         elseif (roots .eq. 500d0) then
            e = GEV500
         elseif (roots .lt. 1000d0) then
            write (msgbuf, 2006) roots, 500d0, 1000d0
            call circem ('MESSAGE', msgbuf)
            elo = GEV500
            ehi = TEV1
            eloval =  500d0
            ehival = 1000d0
         elseif (roots .eq. 1000d0) then
            e = TEV1
         else
            write (msgbuf, 2005) roots, 1000d0
            call circem ('MESSAGE', msgbuf)
            e = TEV1
         endif
      endif
      if (chat .ge. 2) then
         if (e .ge. GEV090) then
            write (msgbuf, 2003) roots, e
            call circem ('MESSAGE', msgbuf)
         else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
            write (msgbuf, 2013) roots, elo, ehi
            call circem ('MESSAGE', msgbuf)
         end if
      endif
      if (e .ge. GEV090) then
         lumi = xa7lum(e,acc,r)
         do 50 i = 0, 7
            a1(i) = xa7(i,e,acc,r)
 50      continue
      elseif (elo .ge. GEV090 .and. ehi .ge. GEV090) then
         lumi = ((roots-eloval)*xa7lum(ehi,acc,r)
     $        + (ehival-roots)*xa7lum(elo,acc,r)) / (ehival - eloval)
         do 51 i = 1, 6
            a1(i) = ((roots-eloval)*xa7(i,ehi,acc,r)
     $           + (ehival-roots)*xa7(i,elo,acc,r)) / (ehival - eloval)
 51      continue
         a1(0) = 1d0 - a1(1) * beta(a1(2)+1d0,a1(3)+1d0)
         a1(7) = a1(4) * beta(a1(5)+1d0,a1(6)+1d0)
      endif
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      elect0 = a1(0) + a1(1) * KIREPS**(a1(3)+1) / (a1(3)+1)
      elect0 = elect0 / KIREPS
      gamma0 = a1(4) * KIREPS**(a1(5)+1) / (a1(5)+1)
      gamma0 = gamma0 / KIREPS
 2001 format ('treating energy ', F6.1, 'GeV as ',  F6.1, 'GeV')
 2002 format ('energy ', F6.1, ' not available for ', A6,
     $        ' in revison ', I2)
 2003 format ('mapping energy ', F6.1, ' to energy index ', I2)
 2013 format ('mapping energy ', F6.1, ' to energy indices ',
     $        I2, ' and ', I2)
 2004 format ('energy ', F6.1, 'GeV too low, using spectrum for ',
     *                   F6.1, 'GeV')
 2005 format ('energy ', F6.1, 'GeV too high, using spectrum for ',
     *                   F6.1, 'GeV')
 2006 format ('energy ', F6.1, 'GeV interpolated between ',
     *                   F6.1, ' and ', F6.1, 'GeV')
      end
      subroutine circel (l)
      implicit none
      double precision l
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      l = lumi
      end
      double precision function circee (x1, x2)
      implicit none
      double precision x1, x2
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision d1, d2
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      circee = -1.0
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      if (x1 .eq. 1d0) then
         d1 = a1(0)
      elseif (x1 .lt. 1d0 .and. x1 .gt. 0d0) then
         d1 = a1(1) * x1**a1(2) * (1d0 - x1)**a1(3)
      elseif (x1 .eq. -1d0) then
         d1 = 1d0 - a1(0)
      else
         d1 = 0d0
      endif
      if (x2 .eq. 1d0) then
         d2 = a1(0)
      elseif (x2 .lt. 1d0 .and. x2 .gt. 0d0) then
         d2 = a1(1) * x2**a1(2) * (1d0 - x2)**a1(3)
      elseif (x2 .eq. -1d0) then
         d2 = 1d0 - a1(0)
      else
         d2 = 0d0
      endif
      circee = d1 * d2
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
      double precision function circeg (x1, x2)
      implicit none
      double precision x1, x2
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision d1, d2
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      circeg = -1.0
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      if (x1 .eq. 1d0) then
         d1 = a1(0)
      elseif (x1 .lt. 1d0 .and. x1 .gt. 0d0) then
         d1 = a1(1) * x1**a1(2) * (1d0 - x1)**a1(3)
      elseif (x1 .eq. -1d0) then
         d1 = 1d0 - a1(0)
      else
         d1 = 0d0
      endif
      if (x2 .lt. 1d0 .and. x2 .gt. 0d0) then
         d2 = a1(4) * x2**a1(5) * (1d0 - x2)**a1(6)
      elseif (x2 .eq. -1d0) then
         d2 = a1(7)
      else
         d2 = 0d0
      endif
      circeg = d1 * d2
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
      double precision function circgg (x1, x2)
      implicit none
      double precision x1, x2
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision d1, d2
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      circgg = -1.0
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      if (x1 .lt. 1d0 .and. x1 .gt. 0d0) then
         d1 = a1(4) * x1**a1(5) * (1d0 - x1)**a1(6)
      elseif (x1 .eq. -1d0) then
         d1 = a1(7)
      else
         d1 = 0d0
      endif
      if (x2 .lt. 1d0 .and. x2 .gt. 0d0) then
         d2 = a1(4) * x2**a1(5) * (1d0 - x2)**a1(6)
      elseif (x2 .eq. -1d0) then
         d2 = a1(7)
      else
         d2 = 0d0
      endif
      circgg = d1 * d2
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
      double precision function beta (a, b)
      implicit none
      double precision a, b
      double precision dlogam
      beta = exp (dlogam(a) + dlogam(b) - dlogam(a+b))
      end
CERNLIB C304
      DOUBLE PRECISION FUNCTION DLOGAM(X)
      IMPLICIT NONE
      DOUBLE PRECISION P1(7),Q1(7),P2(7),Q2(7),P3(7),Q3(7),C(5),XL(5)
      DOUBLE PRECISION X,Y,ZERO,ONE,TWO,HALF,AP,AQ
      INTEGER I
      DATA ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/, HALF /0.5D0/
      DATA XL /0.0D0,0.5D0,1.5D0,4.0D0,12.0D0/
      DATA P1
     1/+3.84287 36567 460D+0, +5.27068 93753 010D+1,
     2 +5.55840 45723 515D+1, -2.15135 13573 726D+2,
     3 -2.45872 61722 292D+2, -5.75008 93603 041D+1,
     4 -2.33590 98949 513D+0/
      DATA Q1
     1/+1.00000 00000 000D+0, +3.37330 47907 071D+1,
     2 +1.93877 84034 377D+2, +3.08829 54973 424D+2,
     3 +1.50068 39064 891D+2, +2.01068 51344 334D+1,
     4 +4.57174 20282 503D-1/
      DATA P2
     1/+4.87402 01396 839D+0, +2.48845 25168 574D+2,
     2 +2.17973 66058 896D+3, +3.79751 24011 525D+3,
     3 -1.97780 70769 842D+3, -3.69298 34005 591D+3,
     4 -5.60177 73537 804D+2/
      DATA Q2
     1/+1.00000 00000 000D+0, +9.50999 17418 209D+1,
     2 +1.56120 45277 929D+3, +7.23400 87928 948D+3,
     3 +1.04595 76594 059D+4, +4.16994 15153 200D+3,
     4 +2.76785 83623 804D+2/
      DATA P3
     1/-6.88062 40094 594D+3, -4.30699 69819 571D+5,
     2 -4.75045 94653 440D+6, -2.94234 45930 322D+6,
     3 +3.63218 04931 543D+7, -3.35677 82814 546D+6,
     4 -2.48043 69488 286D+7/
      DATA Q3
     1/+1.00000 00000 000D+0, -1.42168 29839 651D+3,
     2 -1.55528 90280 854D+5, -3.41525 17108 011D+6,
     3 -2.09696 23255 804D+7, -3.45441 75093 344D+7,
     4 -9.16055 82863 713D+6/
      DATA C
     1/ 1.12249 21356 561D-1,  7.95916 92961 204D-2,
     1 -1.70877 94611 020D-3,  9.18938 53320 467D-1,
     2  1.34699 05627 879D+0/
      IF(X .LE. XL(1)) THEN
       print *, 'ERROR: DLOGAM non positive argument: ', X
       DLOGAM=ZERO
      ENDIF
      IF(X .LE. XL(2)) THEN
       Y=X+ONE
       AP=P1(1)
       AQ=Q1(1)
       DO 2 I = 2,7
          AP=P1(I)+Y*AP
    2     AQ=Q1(I)+Y*AQ
       Y=-LOG(X)+X*AP/AQ
      ELSEIF(X .LE. XL(3)) THEN
       AP=P1(1)
       AQ=Q1(1)
       DO 3 I = 2,7
          AP=P1(I)+X*AP
    3     AQ=Q1(I)+X*AQ
       Y=(X-ONE)*AP/AQ
      ELSEIF(X .LE. XL(4)) THEN
       AP=P2(1)
       AQ=Q2(1)
       DO 4 I = 2,7
          AP=P2(I)+X*AP
    4     AQ=Q2(I)+X*AQ
       Y=(X-TWO)*AP/AQ
      ELSEIF(X .LE. XL(5)) THEN
       AP=P3(1)
       AQ=Q3(1)
       DO 5 I = 2,7
          AP=P3(I)+X*AP
    5     AQ=Q3(I)+X*AQ
       Y=AP/AQ
      ELSE
       Y=ONE/X**2
       Y=(X-HALF)*LOG(X)-X+C(4)+(C(1)+Y*(C(2)+Y*C(3)))/
     1                                        ((C(5)+Y)*X)
      ENDIF
      DLOGAM=Y
      END
      double precision function kirke (x1, x2, p1, p2)
      implicit none
      double precision x1, x2
      integer p1, p2
      double precision kirkee, kirkeg, kirkgg
      integer ELECTR, POSITR, PHOTON
      parameter (ELECTR =  11)
      parameter (POSITR = -11)
      parameter (PHOTON =  22)
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      kirke = -1.0
      if (abs(p1) .eq. ELECTR) then
         if (abs(p2) .eq. ELECTR) then
            kirke = kirkee (x1, x2)
         elseif (p2 .eq. PHOTON) then
            kirke = kirkeg (x1, x2)
         endif
      elseif (p1 .eq. PHOTON) then
         if (abs(p2) .eq. ELECTR) then
            kirke = kirkeg (x2, x1)
         elseif (p2 .eq. PHOTON) then
            kirke = kirkgg (x1, x2)
         endif
      endif
      end
      double precision function kirkee (x1, x2)
      implicit none
      double precision x1, x2
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision d1, d2
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      kirkee = -1.0
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      if (x1 .gt. 1d0) then
         d1 = 0d0
      elseif (x1 .ge. (1d0 - KIREPS)) then
         d1 = elect0
      elseif (x1 .ge. 0d0) then
         d1 = a1(1) * x1**a1(2) * (1d0 - x1)**a1(3)
      else
         d1 = 0d0
      endif
      if (x2 .gt. 1d0) then
         d2 = 0d0
      elseif (x2 .ge. (1d0 - KIREPS)) then
         d2 = elect0
      elseif (x2 .ge. 0d0) then
         d2 = a1(1) * x2**a1(2) * (1d0 - x2)**a1(3)
      else
         d2 = 0d0
      endif
      kirkee = d1 * d2         
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
      double precision function kirkeg (x1, x2)
      implicit none
      double precision x1, x2
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision d1, d2
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      kirkeg = -1.0
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      if (x1 .gt. 1d0) then
         d1 = 0d0
      elseif (x1 .ge. (1d0 - KIREPS)) then
         d1 = elect0
      elseif (x1 .ge. 0d0) then
         d1 = a1(1) * x1**a1(2) * (1d0 - x1)**a1(3)
      else
         d1 = 0d0
      endif
      if (x2 .gt. 1d0) then
         d2 = 0d0
      elseif (x2 .gt. KIREPS) then
         d2 = a1(4) * x2**a1(5) * (1d0 - x2)**a1(6)
      elseif (x2 .ge. 0d0) then
         d2 = gamma0
      else
         d2 = 0d0
      endif
      kirkeg = d1 * d2         
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
      double precision function kirkgg (x1, x2)
      implicit none
      double precision x1, x2
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision d1, d2
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      kirkgg = -1.0
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      if (x1 .gt. 1d0) then
         d1 = 0d0
      elseif (x1 .gt. KIREPS) then
         d1 = a1(4) * x1**a1(5) * (1d0 - x1)**a1(6)
      elseif (x1 .ge. 0d0) then
         d1 = gamma0
      else
         d1 = 0d0
      endif
      if (x2 .gt. 1d0) then
         d2 = 0d0
      elseif (x2 .gt. KIREPS) then
         d2 = a1(4) * x2**a1(5) * (1d0 - x2)**a1(6)
      elseif (x2 .ge. 0d0) then
         d2 = gamma0
      else
         d2 = 0d0
      endif
      kirkgg = d1 * d2         
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
      subroutine girce (x1, x2, p1, p2, rng)
      implicit none
      double precision x1, x2
      integer p1, p2
      external rng
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision u, w, circgg
      integer ELECTR, POSITR, PHOTON
      parameter (ELECTR =  11)
      parameter (POSITR = -11)
      parameter (PHOTON =  22)
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
 99   continue
      w = 1d0 / (1d0 + circgg (-1d0, -1d0))
      call rng (u)
      if (u*u .le. w) then
         p1 = POSITR
      else
         p1 = PHOTON
      endif
      call rng (u)
      if (u*u .le. w) then
         p2 = ELECTR
      else
         p2 = PHOTON
      endif
      if (abs(p1) .eq. ELECTR) then
         if (abs(p2) .eq. ELECTR) then
            call gircee (x1, x2, rng)
         elseif (p2 .eq. PHOTON) then
            call girceg (x1, x2, rng)
         endif
      elseif (p1 .eq. PHOTON) then
         if (abs(p2) .eq. ELECTR) then
            call girceg (x2, x1, rng)
         elseif (p2 .eq. PHOTON) then
            call gircgg (x1, x2, rng)
         endif
      endif
      if ((x1 .lt. x1m) .or. (x2 .lt. x2m)) goto 99
      end
      subroutine gircee (x1, x2, rng)
      implicit none
      double precision x1, x2
      external rng
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision u, girceb
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      call rng (u)
      if (u .le. a1(0)) then
         x1 = 1d0
      else
         x1 = 1d0 - girceb (0d0, 1d0-x1m, a1(3)+1d0, a1(2)+1d0, rng)
      endif
      call rng (u)
      if (u .le. a1(0)) then
         x2 = 1d0
      else
         x2 = 1d0 - girceb (0d0, 1d0-x2m, a1(3)+1d0, a1(2)+1d0, rng)
      endif
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
      subroutine girceg (x1, x2, rng)
      implicit none
      double precision x1, x2
      external rng
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision u, girceb
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      call rng (u)
      if (u .le. a1(0)) then
         x1 = 1d0
      else
         x1 = 1d0 - girceb (0d0, 1d0-x1m, a1(3)+1d0, a1(2)+1d0, rng)
      endif
      x2 = girceb (x2m, 1d0, a1(5)+1d0, a1(6)+1d0, rng)
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
      subroutine gircgg (x1, x2, rng)
      implicit none
      double precision x1, x2
      external rng
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      double precision girceb
      if (magic .ne. MAGIC0) then
         call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
      endif
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      x1 = girceb (x1m, 1d0, a1(5)+1d0, a1(6)+1d0, rng)
      x2 = girceb (x2m, 1d0, a1(5)+1d0, a1(6)+1d0, rng)
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 7) then
         call circem ('PANIC', 'versions >6 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
      end
      double precision function girceb (xmin, xmax, a, b, rng)
      implicit none
      double precision xmin, xmax, a, b
      external rng
      double precision t, p, u, umin, umax, x, w
      if ((a .gt. 1d0) .or. (b .lt. 1d0)) then
         girceb = -1d0
         call circem ('ERROR', 'beta-distribution expects a<=1<=b')
         return
      endif
      t = (1d0 - a) / (b + 1d0 - a)
      p = b*t / (b*t + a * (1d0 - t)**b)
      if (xmin .le. 0d0) then
         umin = 0d0
      elseif (xmin .lt. t) then
         umin = p * (xmin/t)**a
      elseif (xmin .eq. t) then
         umin = p
      elseif (xmin .lt. 1d0) then
         umin = 1d0 - (1d0 - p) * ((1d0 - xmin)/(1d0 - t))**b
      else
         umin = 1d0
      endif
      if (xmax .ge. 1d0) then
         umax = 1d0
      elseif (xmax .gt. t) then
         umax = 1d0 - (1d0 - p) * ((1d0 - xmax)/(1d0 - t))**b
      elseif (xmax .eq. t) then
         umax = p
      elseif (xmax .gt. 0d0) then
         umax = p * (xmax/t)**a
      else
         umax = 0d0
      endif
      if (umax .lt. umin) then
         girceb = -1d0
         return
      endif
 10   continue
      call rng (u)
      u = umin + (umax - umin) * u
      if (u .le. p) then
         x = t * (u/p)**(1d0/a)
         w = (1d0 - x)**(b-1d0)
      else
         x = 1d0 - (1d0 - t) * ((1d0 - u)/(1d0 - p))**(1d0/b)
         w = (x/t)**(a-1d0)
      endif
         call rng (u)
      if (w .le. u) goto 10
      girceb = x
      end
      subroutine circem (errlvl, errmsg)
      implicit none
      character*(*) errlvl, errmsg
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision KIREPS
      parameter (KIREPS = 1D-6)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      double precision elect0, gamma0
      common /circom/ elect0, gamma0
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r, ehi, elo
      common /circom/ e, r, ehi, elo
      save /circom/
      integer errcnt
      save errcnt
      data errcnt /0/
      if (errlvl .eq. 'MESSAGE') then
         print *, 'circe:message: ', errmsg
      elseif (errlvl .eq. 'WARNING') then
         if (errcnt .lt. 100) then
            errcnt = errcnt + 1
            print *, 'circe:warning: ', errmsg
         elseif (errcnt .eq. 100) then
            errcnt = errcnt + 1
            print *, 'circe:message: more than 100 messages'
            print *, 'circe:message: turning warnings off'
         endif
      elseif (errlvl .eq. 'ERROR') then
         if (errcnt .lt. 200) then
            errcnt = errcnt + 1
            print *, 'circe:error:   ', errmsg
         elseif (errcnt .eq. 200) then
            errcnt = errcnt + 1
            print *, 'circe:message: more than 200 messages'
            print *, 'circe:message: turning error messages off'
         endif
      elseif (errlvl .eq. 'PANIC') then
         if (errcnt .lt. 300) then
            errcnt = errcnt + 1
            print *, 'circe:panic:   ', errmsg
         elseif (errcnt .eq. 300) then
            errcnt = errcnt + 1
            print *, 'circe:message: more than 300 messages'
            print *, 'circe:message: turning panic messages off'
         endif
      else
         print *, 'circe:panic:    invalid error code ', errlvl
      endif
      end
