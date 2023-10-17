v p1,p2,p3,p4,e1,e2,e3,e4;
v r1,r2,r3,r4;
s s,t,u,MW,MZ,MH,d12Z,d14Z,d12g,d14g,d12H,d14H,cw,sw;
i i,j,k,l,m,n;

* All momenta incoming WWWW --> 0

* The k(mu)*k(nu)/MZ^2 terms in the Z ropagator give zero diagram by diagram 

* factored out i_*g^2; MW/cw=MZ

*        W+ (1) ----------------- W- (2)
*                       +
*                       +   gamma,Z
*                       +
*                       +
*        W+ (3) ----------------- W- (4)


L AWWWW1Z = cw^2*e1(i)*e2(j)*
            ( d_(i,j)*(p2(k)-p1(k)) + d_(j,k)*(-p1(i)-2*p2(i)) + d_(k,i)*(2*p1(j)+p2(j)))
           *( d_(k,l)-(p1(k)+p2(k))*(-p3(l)-p4(l))/MZ^2 )
           *( d_(m,n)*(p4(l)-p3(l)) + d_(n,l)*(-2*p4(m)-p3(m)) + d_(l,m)*(2*p3(n)+p4(n)))
           *e3(m)*e4(n)/d12Z;

L AWWWW1g = sw^2*e1(i)*e2(j)*
            ( d_(i,j)*(p2(k)-p1(k)) + d_(j,k)*(-p1(i)-2*p2(i)) + d_(k,i)*(2*p1(j)+p2(j)))
           *( d_(k,l) )
           *( d_(m,n)*(p4(l)-p3(l)) + d_(n,l)*(-2*p4(m)-p3(m)) + d_(l,m)*(2*p3(n)+p4(n)))
           *e3(m)*e4(n)/d12g;

*        W+ (1) ----------------- W- (4)
*                       +
*                       +   gamma,Z
*                       +
*                       +
*        W+ (3) ----------------- W- (2)


L AWWWW2Z = cw^2*e1(i)*e4(j)*
            ( d_(i,j)*(p4(k)-p1(k)) + d_(j,k)*(-p1(i)-2*p4(i)) + d_(k,i)*(2*p1(j)+p4(j)))
           *( d_(k,l)-(p1(k)+p4(k))*(-p3(l)-p2(l))/MZ^2 )
           *( d_(m,n)*(p2(l)-p3(l)) + d_(n,l)*(-2*p2(m)-p3(m)) + d_(l,m)*(2*p3(n)+p2(n)))
           *e3(m)*e2(n)/d14Z;

L AWWWW2g = sw^2*e1(i)*e4(j)*
            ( d_(i,j)*(p4(k)-p1(k)) + d_(j,k)*(-p1(i)-2*p4(i)) + d_(k,i)*(2*p1(j)+p4(j)))
           *( d_(k,l) )
           *( d_(m,n)*(p2(l)-p3(l)) + d_(n,l)*(-2*p2(m)-p3(m)) + d_(l,m)*(2*p3(n)+p2(n)))
           *e3(m)*e2(n)/d14g;


*                W+ (3)
*                       +
*                       +
*                       +
*        W+ (1) ----------------- W- (2)
*                       +
*                       +
*                       +
*                W- (4) 

L AWWWW3 =  2*e1.e3*e2.e4 - e1.e2*e3.e4 - e1.e4*e2.e3;


*        W+ (1) ----------------- W- (2)
*                       =
*                       =   H
*                       =
*                       =
*        W+ (3) ----------------- W- (4)


L AWWWW4 = - MW^2/cw^2*(e1.e2*e3.e4)/d12H;


*        W+ (1) ----------------- W- (4)
*                       =
*                       =   H
*                       =
*                       =
*        W+ (3) ----------------- W- (2)


L AWWWW5 = - MW^2/cw^2*(e1.e4*e2.e3)/d14H;


* Verify cancellation of terms O(s^2): A1 does not include terms O(s^6)
* L A = AWWWW1g + AWWWW1Z + AWWWW2g + AWWWW2Z + AWWWW3;
 L A1 = (AWWWW1g + AWWWW1Z + AWWWW2g + AWWWW2Z + AWWWW3)*d12Z*d12g*d14Z*d14g;
 L A2 = ( AWWWW4 + AWWWW5)*d12H*d14H;

L A3 = (AWWWW1g + AWWWW1Z + AWWWW2g + AWWWW2Z + AWWWW3 + AWWWW4 + AWWWW5)*d12Z*d12g*d14Z*d14g*d12H*d14H;
id d12Z = 2*MW^2 + 2*p1.p2 - MZ^2; 
id d12g = 2*MW^2 + 2*p1.p2; 
id d14Z = 2*MW^2 + 2*p1.p4 - MZ^2; 
id d14g = 2*MW^2 + 2*p1.p4; 

id d12H = 2*MW^2 + 2*p1.p2 - MH^2; 
id d14H = 2*MW^2 + 2*p1.p4 - MH^2; 

id p1.p1 = MW^2;
id p2.p2 = MW^2;
id p3.p3 = MW^2;
id p4.p4 = MW^2;

id p1.e1 = 0;
id p2.e2 = 0;
id p3.e3 = 0;
id p4.e4 = 0;

print;
.sort

* ri --> 0 as E --> infinity

id e1 = p1/MW + r1;
id e2 = p2/MW + r2;
id e3 = p3/MW + r3;
id e4 = p4/MW + r4;

id p4 = -p1 -p2 -p3;

id p1.p1 = MW^2;
id p2.p2 = MW^2;
id p3.p3 = MW^2;


* High Energy test

id p1.p2 = (s - 2*MW^2)/2;
* id p3.p4 = (s - 2*MW^2)/2;
id p1.p3 = (t - 2*MW^2)/2;
* id p2.p4 = (t - 2*MW^2)/2;
* id p1.p4 = (u - 2*MW^2)/2;
id p2.p3 = (u - 2*MW^2)/2;

id u = -s -t + 4*MW^2;

id p1.r1 = -MW;
id p2.r2 = -MW;
id p3.r3 = -MW;
* id p4.r4 = -MW;

id MZ = MW/cw;

id sw^2 = 1 - cw^2;

id p1.r4 = - p2.r4 - p3.r4 + MW;


Bracket s,t;

print;
.end






**************************************************

* High Energy test

id MW = 0;

id MZ =0;

id e1 = p1;
id e2 = p2;
id e3 = p3;
id e4 = p4;
id p1.p2/d12Z = 1/2;
id p1.p2/d12g = 1/2;
id p1.p4/d14Z = 1/2;
id p1.p4/d14g = 1/2;

id p1.p2=s;
id p3.p4=s;
id p1.p3=t;
id p2.p4=t;
id p1.p4=u;
id p2.p3=u;
id u=-s-t;

id cw^2 = 1 - sw^2;

print;

**************************************************
