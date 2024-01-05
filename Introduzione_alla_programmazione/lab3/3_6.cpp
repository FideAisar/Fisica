#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>

using namespace std;

void swap(int *v,int a,int b){
	int c = v[a];
	v[a] = v[b];
	v[b] = c;
}
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
int* vettore_casuale(int k, int mm, int MM){
  	int *r;

  	r = alloca_vettore(k);

  	for(int i=0;i<=k;i++){
    		r[i] = numero_casuale(mm,MM);	  
  	}
  	return r;
}
int index_of_the_max(int *v, int inf, int sup){
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
int* ord(int *v,int n){
  	int b,c;
  	for(int i=1;i<=n;i++){
		swap(v,index_of_the_max(v,0,n-i),n-i);   
  	}
  	return v;
}
int index_perno(int *v, int inf, int sup){
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
			swap(v,i,j);
      			i++;
      			j--;
    		}
  	}
   	swap(v,inf,j);
  	return j;
}

void quick_sort(int *v, int inf, int sup){
  	int max;
	
	if(inf >= sup){
		return;
	}
	max = index_perno(v,inf,sup);
	quick_sort(v,inf,max);
	quick_sort(v,max+1,sup);
}

// MAIN // MAIN // MAIN // MAIN // MAIN // MAIN // MAIN // MAIN



int main()
{	
	int n = 6;
   	int *v;
	
	v = legge_vettore(n);
	stampa_vettore(v,n);


	quick_sort(v,0,n-1);
	stampa_vettore(v,n);
	delete[] v;
   	return 0;
}

