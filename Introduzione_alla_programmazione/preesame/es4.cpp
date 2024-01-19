#include<iostream>

using namespace std;

double fracIter(int n){
	double sum = 0.;
	int i=1;
	while(i<=n){
		sum += 1./i;
		i++;
	}
	return sum;
}
double fracRec(int n){
	if(n==0){
		return 0;
	}
	return 1./n + fracRec(n-1);
}


int main(){
	cout << fracIter(3) << endl;
	cout << fracRec(3) << endl;
	return 0;
}
