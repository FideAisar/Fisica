#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
#include<cmath>

using namespace std;


int ones_counter(int v[],int n){
    int c=0;
    for(int i=0;i<n;i++){
        if(v[i] == 1){
            c++;
        }
        else{
            continue;
        }
    }
    return c;
}

void nuovo_vettore(int n,int v[],int v1[]){   // prende la grandezza e il vettore precedente
    int j=0;
    for (int i = 0; i < n; i++, j++) {
        if(v[i] == 1){
            v1[j] = v[i];
            v1[j+1] = 0;
            j++;
        } else {
            v1[j] = v[i];
        }
    }
}

void stampa_vettore(int v[], int n){
    for(int i=0;i<n;i++){
        cout << v[i] << " ";
    }
}


int main(){

    int n=6;
    int a[n] = {0, 1, 5, 1, 0, -3};
    
    int n1=ones_counter(a,n) + n;
    int a1[n1];

    nuovo_vettore(n1,a,a1);
    stampa_vettore(a1,n1);

    return 0;
}