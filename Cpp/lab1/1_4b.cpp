#include <iostream>
using namespace std;

int main(){
	int r=1,i=1,righe;
	cout << "numero righe: " << endl;
	cin >> righe;
	
	r = righe;
	while(i <= righe){
		while(righe <= r){
			cout << "*";
			r--;
		}
		i++;
		r = r + i;
		cout << "" << endl;
	}
	return 0;
}
























