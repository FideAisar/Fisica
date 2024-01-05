#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<ctime>
using namespace std;


void stampa_vettore(int *v, int n){
	for(int i=0;i<n;i++){
		cout << v[i] << " ";
	}
}

bool func(int *v, int n){
	if(n==1){
		return true;
	}
	for(int i=0;i<n;i++){
		if(v[i]<=v[i+1]){
			continue;
		}
		else{
			return false;
			break;
		}
	}
	return true;
}


int main(){
	
	int n, *v;
	cout << "scegli grandezza vettore: ";
	cin >> n;
	v=new int[n];

	for(int i=0;i<n;i++){
			cout << "elemento" << i+1 << ": ";
			cin >> v[i];
	}

	if(func(v,n)){
		cout << "Il vettore è crescente: "; 
		stampa_vettore(v,n);
	}
	else{
		cout << "Il vettore non è crescente: "; 
		stampa_vettore(v,n);
	}
	
	delete[] v;
	return 0;
}
