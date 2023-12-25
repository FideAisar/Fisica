#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
#include<cmath>

using namespace std;

const int r = 3;
const int c = 2;

void stampa_vettore(int v[]){
    for(int i=0;i<r;i++){
        cout << v[i] << " ";
    }
    cout << endl;
}

void stampa_matrice(int m[][c]){
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            cout << " " << m[i][j] << " ";
        }
        cout << endl;
    }
}

void matrice1(int m[][c], int m1[][c], int v[], int r, int c){
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            m1[i][j] = m[v[i]-1][j];
        }
    }
}

int main(){

    int v[r] = {2,1,3};
    int m[r][c] = {{1,2},{2,4},{3,3}};
    int m1[r][c];
    matrice1(m,m1,v,r,c);    
    

    cout << "Matrice originale: " << endl;
    stampa_matrice(m);

    cout << "La matrice deve essere permutata nel modo descritt dal vettore " << endl;
    stampa_vettore(v);

    cout << "Matrice permutata: " << endl;
    stampa_matrice(m1);

	return 0;
}