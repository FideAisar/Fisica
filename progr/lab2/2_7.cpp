#include <iostream>
using namespace std;

int pow(int a, int n){
	if(n==0){
		return 1;
	}
	return a * pow(a,n-1);
}

int pow1(int a, int n){
	if(n%2 == 0){
		return pow1(a,n/2)*pow1(a,n/2);
	}
	else {
		return a*pow1(a,n/2)*pow1(a,n/2);
	}
	
}

int main(){	
	int a,n;
	cout << "Inserire base potenza 'a' e espnente 'n': " << endl;
	cin >> a >> n;
	cout << pow(a,n) << endl;
	cout << pow1(a,n) << endl;
	return 0;
}	
