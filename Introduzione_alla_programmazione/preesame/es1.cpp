#include<iostream>

using namespace std;

double func(int n, double x){
	double a=1;
	int i=0;

	if(n==0){
		return 1;
	}

        return x*func(n-1,x);

        //while(i <= n){
	//	x *= a;
	//	i++;
	//}
	//return x;
}

int main(){
	int n ;
	double x;

	cout << "X alla n" << endl;
	cin >> x;
	cin >> n;
        cout << func(n,x) << endl;

	return 0;
}
