v p1,p2,p3,p4,e1,e2,e3,e4,q13,q14;
v r1,r2,r3,r4;
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
          *( d_(i,j)-(p1(i)+p3(i))*(-p2(j)-p4(j))/MW^2 )
          *( 2*p2.e4*e2(j) - 2*p4.e2*e4(j) + e2.e4*(p4(j)-p2(j)) )/d13W;


*        W+ (1) ----------------- Z (4)
*                       +
*                       +  W
*                       +
*                       +
*        W- (2) ----------------- Z (3)


L AZZWW2 = ( -2*p1.e4*e1(i) + 2*p4.e1*e4(i) + e1.e4*(p1(i)-p4(i)) )
          *( d_(i,j)-(p1(i)+p4(i))*(-p2(j)-p3(j))/MW^2 )
          *( 2*p2.e3*e2(j) - 2*p3.e2*e3(j) + e2.e3*(p3(j)-p2(j)) )/d14W;


*                 Z (3)
*                       +
*                       +
*                       +
*        W+ (1) ----------------- Z (4)
*                       +
*                       +
*                       +
*                W- (2) 


L AZZWW3 = - ( 2*e1.e2*e3.e4 - e1.e3*e2.e4 - e1.e4*e2.e3);


*        W+ (1) ----------------- W- (2)
*                       =
*                       =   H
*                       =
*                       =
*         Z (3) ----------------- Z (4)


* L AZZWW4 = - MW*MZ/cw^3*(e1.e2*e3.e4)/d12H;
L AZZWW4 = - MZ^2/cw^2*(e1.e2*e3.e4)/d12H;

* Verify cancellation of terms O(s^2): A1 does not include terms O(s^4)
* L A1 = (AZZWW1 + AZZWW2 + AZZWW3)*d14W*d13W;

* Verify cancellation of terms O(s): A2 does not include terms O(s^4)
L A2 = (AZZWW1 + AZZWW2 + AZZWW3 + AZZWW4)*d14W*d13W*d12H;

id d14W = MZ^2 + 2*p1.p4; 
id d13W = MZ^2 + 2*p1.p3;
id d12H = 2*MW^2 + 2*p1.p2 - MH^2;

id p1.p1 = MW^2;
id p2.p2 = MW^2;
id p3.p3 = MZ^2;
id p4.p4 = MZ^2;

id p1.e1 = 0;
id p2.e2 = 0;
id p3.e3 = 0;
id p4.e4 = 0;

print;
.sort

* ri --> 0 as E --> infinity

id e1 = p1/MW + r1;
id e2 = p2/MW + r2;
id e3 = p3/MZ + r3;
id e4 = p4/MZ + r4;

* High Energy test

id p1.p2 = (s - 2*MW^2)/2;
id p3.p4 = (s - 2*MZ^2)/2;
id p1.p3 = (t - MW^2 - MZ^2)/2;
id p2.p4 = (t - MW^2 - MZ^2)/2;
id p1.p4 = (u - MW^2 - MZ^2)/2;
id p2.p3 = (u - MW^2 - MZ^2)/2;

id u = -s -t + 2*MW^2 + 2*MZ^2;

id p1.r1 = -MW;
id p2.r2 = -MW;
id p3.r3 = -MZ;
id p4.r4 = -MZ;

id p4 = -p1 -p2 -p3;

id p1.r1 = -MW;
id p2.r2 = -MW;
id p3.r3 = -MZ;


Bracket s,t;

print;
.end





**************************************************

* High Energy test: E^4 terms

L A = AZZWW1 + AZZWW2 + AZZWW3;

id MZ = 0;

print;

id MW =0;

id e1 = p1;
id e2 = p2;
id e3 = p3;
id e4 = p4;
id p1.p3/d13W = 1/2;
id p1.p4/d14W = 1/2;

id p1.p2=s;
id p3.p4=s;
id p1.p3=t;
id p2.p4=t;
id p1.p4=u;
id p2.p3=u;
id u=-s-t;

print;

**************************************************