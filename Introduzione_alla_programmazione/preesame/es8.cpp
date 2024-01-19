#include<iostream>

using namespace std;


void f(int* v,int n, int d){
	int* w = new int[n];

	for(int i=0;i<d;i++){
		w[i] = v[n-d+i];
	}
	for(int k=d;k<n;k++){
		w[k] = v[k-d];
	}
	for(int j=0;j<n;j++){
		v[j] = w[j];
	}
}


int main(){

	int n = 5;
	int *v=new int[n];

	for(int i=0;i<n;i++){
		cout << i <<": ";
		cin >> v[i];
	}

	f(v,n,2);

	for(int i=0;i<n;i++){
		cout << v[i] << " " ;
	}
	return 0;
}
