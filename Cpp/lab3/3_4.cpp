#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>

using namespace std;

int* alloca_vettore(int n){
  int *v;
  v=new int[n];
  return v;
}

// libera la memoria occupata da un array di int
void dealloca_vettore(int *v){
  delete[] v;
}

int* legge_vettore(int n){
  int *v;
  v=alloca_vettore(n); // alloca memoria per l'array
  for(int i=0;i<n;i++){
    cout << "elemento numero " << i+1 << ":";
    cin >> v[i];
  }
  return v; // restituisce l'indirizzo dell'array
}

void stampa_vettore(int *s, int k){
  for(int j=0;j<k;j++)
    cout << s[j] << " ";
  cout << endl;
}

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

int index_max(int *v, int inf, int sup){
  int m=v[inf], index;

  for(int i=inf;i<=sup;i++){
    if(m<=v[i]){
      m = v[i];
      index = i;
    }
    else{
      m = m;
    }
  }
  return index;
}

int index_min(int *v, int inf, int sup){
  int m=v[inf], index;

  for(int i=inf;i<=sup;i++){
    if(m>=v[i]){
      m = v[i];
      index = i;
    }
    else{
      m = m;
    }
  }
  return index;
}

int* ord(int *v,int n){
  int b,c;
  for(int i=1;i<=n;i++){
    b = v[index_max(v,0,n-i)];
    v[index_max(v,0,n-i)] = v[n-i];
    v[n-i] = b;

    
  }
  return v;
}



int main()
{
  int *v, *r, n;
  n = 100;
  v = vettore_casuale(n,0,n);
  stampa_vettore(v,10);

  r = ord(v,10);
  stampa_vettore(r,10);

  double start = clock();
  selection_sort(n,v);
  double end = clock();
  double seconds = (end - start) / CLOCKS_PER_SEC;
  cout << seconds;
  return 0;
}

