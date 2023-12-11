#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
#include<cmath>

using namespace std;

int pow(int x,int n){
    int x1 = x;
    for(int i=2;i<=n;i++){
        x = x1*x;
    }
    return x;
}

int approx(int n){
    int x=0,a=0;
    while(a < 1000){
        a = pow(x,n);
        if(a<1000){
            x++; 
        }   
    }
    return x;
}

int main(){

    int p;
    cout << "inserire n (x^n = 1000): ";
    cin >> p;

    cout << approx(p) << endl;
    return 0;
}