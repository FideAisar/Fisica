#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
#include<cmath>

using namespace std;

int fac(int x){
    int i,f=1;
    while(i<=x){
        f *= i;
        i++;
    }
    return f;
}

float func(float x, float e){
    int x1=0, x2=0, n=0;

    while(abs(x2-x1) < e){
        x1 += (pow(-1,n)*pow(x,2*n+1))/fac(2*n+1);
        x2 += (pow(-1,n+1)*pow(x,2*(n+1)+1))/fac(2*(n+1)+1);
        n++;
    }
    return x2;
}


int main(){
	
    float e, x;

    cout << "caloclo sin di: ";
    cin >> x;

    cout << "Inserrire livello di errore e: ";
    cin >> e;

    cout << func(x,e) << endl;
	
	return 0;
}