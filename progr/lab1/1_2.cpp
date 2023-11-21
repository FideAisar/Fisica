#include <iostream>
using namespace std;

int main(){
	int a,b;
	cout << "Verifico che b sia multiplo di a" << endl;
	cout << "a: ";
	cin >> a;
	cout << "b: ";
	cin >> b;

	if(a%b == 0)
		cout << "a è multiplo di b" << endl;
	else{
		cout << "a non è multiplo di b" << endl;
	}
	
	return 0;
}	
