#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
#include<cmath>

using namespace std;

// alloca memoria per un array di n int
// restituisce l'indirizzo dell'array
int* alloca_vettore(int n){
  int *v;
  v=new int[n];
  return v;
}

// libera la memoria occupata da un array di int
void dealloca_vettore(int *v){
  delete[] v;
}

// legge dalla tastiera un array di n int
// l'array è creata nella funzione con allocazione dinamica
// restituisce l'indirizzo dell'array
int* legge_vettore(int n){
  int *v;
  v=alloca_vettore(n); // alloca memoria per l'array
  for(int i=0;i<n;i++){
    cout << "elemento numero " << i+1 << ":";
    cin >> v[i];
  }
  return v; // restituisce l'indirizzo dell'array
}

// stampa un array di k int a video
void stampa_vettore(int *s, int k){
  for(int j=0;j<k;j++)
    cout << s[j] << " ";
  cout << endl;
}

// restituisce un int casuale fra m e M
int numero_casuale(int m, int M){
  int r;
  r=rand()%(M-m+1)+m;
  return r;
}

int* vettore_casuale(int k, int mm, int MM){
  int *r;

  r = alloca_vettore(k);

  for(int i=0;i<=k;i++){
    r[i] = numero_casuale(mm,MM);	  
  }
  return r;
}

int* somma_vettori(int *v1, int *v2, int k){
  int *ris;

  ris = alloca_vettore(k);

  for(int i=0;i<=k;i++){
    ris[i] = v1[i] + v2[i];
  }

  return ris;
}

int* divverenza_vettori(int *v1, int *v2, int k){
  int *ris;

  ris = alloca_vettore(k);

  for(int i=0;i<=k;i++){
    ris[i] = v1[i] - v2[i];
  }

  return ris;
}

int* prod_con_scalare(int *v1, int l, int k){
  int *ris;

  ris = alloca_vettore(k);

  for(int i=0;i<=k;i++){
    ris[i] = v1[i] * l;
  }

  return ris;
}

double* norma_vettore(int *v1, int k){
  double *ris;
  double m = 0;

  for(int i=0;i<=k;i++){
     m = m + v1[i]*v1[i];
  }

  return sqrt(ris);
}













int main()
{
  int *v1,*v2,s;

  v2=vettore_casuale(20, 5, 10);
  stampa_vettore(v2,20);
  delete[] v2;

  return 0;
}

