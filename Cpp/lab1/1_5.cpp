#include <iostream>
using namespace std;

int main(){
	int n,c,i;
	
	cout << "n: ";
	cin >> n;
	c = 1;
	i = 0;
	while(i<n){
		c = c*n;
		n--;
		i++;
	}
	cout << c << endl;
	return 0;
}
