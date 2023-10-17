v p1,p2,p3,p4,e1,e2,e3,e4,q13,q14;
s s,t,u,MW,MZ,MH,d13W,d14W,d12H,cw;
i i,j,k;

* All momenta incoming ZZWW --> 0

* factored out i_*g^2*cw^2; MW/cw=MZ


*        W+ (1) ----------------- Z (3)
*                       +
*                       +   W
*                       +
*                       +
*        W- (2) ----------------- Z (4)


L AZZWW1 = ( -2*p1.e3*e1(i) + 2*p3.e1*e3(i) + e1.e3*(p1(i)-p3(i)) )
          *( d_(i,j) )
          *( 2*p2.e4*e2(j) - 2*p4.e2*e4(j) + e2.e4*(p4(j)-p2(j)) );


 L AZZWW2 = -4*p1.e3*p2.e4*e1.e2 +4*p1.e3*p4.e2*e1.e4
         +4*p3.e1*p2.e4*e2.e3 -4*p3.e1*p4.e2*e3.e4
         -2*p1.e3*e2.e4*(p4.e1-p2.e1) +2*p3.e1*e2.e4*(p4.e3-p2.e3)
         +2*p2.e4*e1.e3*(p1.e2-p3.e2) -2*p4.e2*e1.e3*(p1.e4-p3.e4)
         + e1.e3*e2.e4*(p1.p4-p1.p2-p3.p4+p2.p3);

L DIFF = AZZWW1 - AZZWW2;

id p1.p1 = MW^2;
id p2.p2 = MW^2;
id p3.p3 = MZ^2;
id p4.p4 = MZ^2;

id p1.e1 = 0;
id p2.e2 = 0;
id p3.e3 = 0;
id p4.e4 = 0;

print;
.end
