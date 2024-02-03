#include<iostream>
#include<cmath>
using namespace std;


double somma(int n){
	if(n==0)
		return 0;
	return sin(n) + somma(n-1);
}

int main(){
	cout << somma(5) << endl;

	return 0;
}
