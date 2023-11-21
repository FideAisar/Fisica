#include <iostream>
using namespace std;

int fib(int n){
	int v[n],i;
	v[0] = 1;
	v[1] = 1;
	for(i=2;i<=n;i++){
		v[i] = v[i-1] + v [i-2];
	}
	return v[n-1];
}

int main(){
	int n;
	cout << "Inserire numero intero n: " << endl;
	cin >> n;
	cout << "n-esimo numero di Fibonacci: " << fib(n) << endl;
	return 0;
}
