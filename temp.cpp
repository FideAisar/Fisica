#include <iostream>
using namespace std;



int g=1;
int* f(int x, int &y, int *z) {        // funzione che returna un array dinamico
    if(x>=0) y++;
        x=x+4;
    for(int i1=1;i1<=3;i1++){
        z[i1]+=(y+i1+g);
        g++;
    }

return z+1;
}

int main() {
    int a=1,b=2,c[4]={0,1,2,3};
    int *d;
    d=f(a,b,c);
    c[0]=*d;
    cout << a << " " << b << " " << c[0] << " " << c[1] << " " << c[2] << " "
<< c[3] << " ";

    return 0;
}