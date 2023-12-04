#include <iostream>
using namespace std;

int main(){
	int k,i,l;
	k = 1;
	cout << "lato quadrato: ";
	cin >> l;
	
	while(k<=l){
		i = 1;
		while(i<=l){
			cout << "*";	
			i++;
		}
		cout << endl;
		k++;
	}
	return 0;
}
