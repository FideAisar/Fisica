#include<iostream>

using namespace std;

int cifre(int x){
	int sum=0;
	
	while(x>0){
		sum += x % 10;
		x /= 10;
	}
	return sum;
}

int main(){
	cout << cifre(338) << endl;

	return 0;}
