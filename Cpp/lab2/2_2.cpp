#include <iostream>
#include <cmath>
using namespace std;

bool dist(float xp,float yp,float xc,float yc,float r){
	if(sqrt((xp - xc)*(xp - xc) + (yp - yc)*(yp - yc)) < r){
		return true;}
	else{
		return false;}
}

int main(){
	float xp,yp,xc,yc,r;
	cout << "coordinate puno: " << endl;
	cin >> xp >> yp;
	cout << "coordinate centro cerchio: " << endl;
	cin >> xc >> yc;
	cout << "raggio: " << endl;
	cin >> r;
	
	if(dist(xp,yp,xc,yc,r) == true){ 
		cout << "Nel cerchio" << endl;}
	else{
		cout << "non nel cerchio" << endl;}
	
	return 0;
}
