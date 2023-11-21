#include <iostream>
#include <cmath>

using namespace std;

int main(){
	int j=1;
	double p,p1=0,e;
	cout << "e: ";
	cin >> e;

	p = pow(-1,0)*4/((2*(0)) + 1);

	while(abs(p-p1) >= e){
		p = p + pow(-1,j)*4/((2*(j)) + 1);
		p1 =  p1 + pow(-1,(j-1))*4/((2*(j-1)) + 1);
		j++;
	}
	
	cout << p << endl;

	return 0;
}
