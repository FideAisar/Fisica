#include <iostream>
#include <cmath>
using namespace std;

int main(){
	int n;
	double p;

	cout << "n: ";
	cin >> n;

	for(int j=0;j<=n;j++){
		p = p + pow(-1,j)*4./((2.*j) + 1);
	}

	cout << p << endl;

	return 0;
}
