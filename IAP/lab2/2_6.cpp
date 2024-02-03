#include <iostream>
using namespace std;

int fun(int &x, int &y){
	int b = x;
	x = y;
	y = b;
	return x,y; 
}

int main(){
	int x,y;
	cout << "Inserire due numeri  interi: " << endl;
	cin >> x >> y;
	x,y = fun(x,y);
	cout << x << " " << y << endl;
	return 0;
}
