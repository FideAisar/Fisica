#include<iostream>

using namespace std;

void f(int *v, int n){
	int i=0,j=0,sum=0;
	int *w = new int[n];

	for(i;i<n;i++){
		for(j;j<=i;j++){
			sum += v[j];
		}
		w[i] = sum;
		sum = 0;
		j = 0;
	}
	for(int k=0;k<n;k++){
		v[k] = w[k];
	}
}


int main(){
	int n=5;
	int *v=new int[n];

	cout << "inserire valori vettore" << endl;
	for(int j=0;j<n;j++){
		cout << j << ": ";
		cin >> v[j];
	}

	f(v,n);

	for(int i=0;i<n;i++){
		cout << v[i] << " " << endl;
	}

	return 0;
}
