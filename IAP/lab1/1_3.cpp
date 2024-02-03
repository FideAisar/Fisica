#include <iostream>
using namespace std;

int main(){
	int i,n;
	i = 1;
	cout << "n:";
	cin >> n;
	while(i <= n){
		cout << "*";
		i = i+1;
	}
	cout << "" << endl;
	return 0;
}
