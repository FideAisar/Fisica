#include <iostream>
using namespace std;

int count = 0;

int fib(int n){
	count++;
	if(n<=1){
		return n;
	}
	else{
		return fib(n-1) + fib(n-2);
	}
}


int main(){
	int n;
	cout << "Inserire numero intero n: " << endl;
	cin >> n;
	cout << fib(n) << endl;
	cout << "Chiamate: " << count << endl;
	return 0;
}
