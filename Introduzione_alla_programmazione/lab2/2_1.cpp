#include <iostream>
#include <cmath>
using namespace std;

float pos(float t, float a, float v0, float p0){
	float p = a*(t*t)/2 + v0*t + p0;
	return p;
}

int main(){
	int n;
	float dt,t=0,a,v0,p0;
	cout << "Dammi a, v0, p0: " << endl;
	cin >> a >> v0 >> p0;
	cout << "Ogni quanti secondi vuoi stampare la posizione? " << endl;
	cin >> dt;
	cout << "Quanti tempi vuoi stampare? " << endl;
	cin >> n;
	
	while(t<=n*dt){
		cout << "posizione = " << pos(t,a,v0,p0) << t << endl;
		t = t + dt;
	}
	return 0;
}
