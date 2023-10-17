v p1,p2,p3,p4;
s s,t,u;

Local A = -3*p1.p4*p2.p4 + p1.p2*(3/2*p2.p3+3/2*p2.p4-2*p3.p4)
    -3/2*p1.p4*p2.p4 + 3/2*p2.p4*p3.p4
    -3/2*p2.p3*p1.p3 + 3/2*p2.p3*p3.p4
    -2*p1.p2*p3.p4 + p1.p3*p2.p4 + p1.p4*p2.p3;

id p1.p2=s;
id p3.p4=s;
id p1.p3=t;
id p2.p4=t;
id p1.p4=u;
id p2.p3=u;
id u=-s-t;
print;
.end