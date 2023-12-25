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

int* ord_per(int *v, int inf, int sup){
  int p,i,j,b,c;

  p = v[inf];
  i = inf + 1;
  j = sup;

  while(i<=j){
    if(v[i]<=p){
      i++;
    }
    else if(v[i]>p && v[j]>p){
      j--;
    }
    else if(v[i]>p && v[j]<=p){
      b = v[j];
      v[j] = v[i];
      v[i] = b;
      i++;
      j--;
    }
  }
  c = v[inf];
  v[inf] = v[j];
  v[j] = c;
  return j;
}



int main()
{
   int *v,*r;
   v = vettore_causale(10,0,10);
   r = ord_per(v,2,3);
   stampa_vettore(r,10);
   return 0;
}

