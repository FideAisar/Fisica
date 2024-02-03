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

int index_max_el(int *v, int inf, int sup){
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

int index_min_el(int *v, int inf, int sup){
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


int main()
{
    int *v;	
    
    v = legge_vettore(10);
   
    cout << "Vettore v: " << endl;
    stampa_vettore(v,10);

    cout << "Massimo e minimo " << endl;	
    cout << index_max_el(v,3,5) << endl;
    cout << index_min_el(v,3,5) << endl;
    return 0;
}

