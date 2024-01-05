#include <iostream>
using namespace std;

int main(){
	double i,n,max,min;
	cout << "n: " << endl;
	cin >> n;
	max = n;
	min = n;
	for(i=0;i<=9;i++){
		cout << "n: " << endl;
		cin >> n;
		if(n >= max)
			max = n;
		cout << "max attuale: " << max << endl;
		if(n <= min)
			min = n;
		cout << "min attuale; " << min << endl;
	}
	return 0;
}	
