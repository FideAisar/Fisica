#include<iostream>

using namespace std;

int f(int x, int y){
	int max,min,sum=0;
	if(x>y){
		max = x;
		min = y;
	}
	else{
		max = y;
		min = x;
	}
	for(int i=min;i<=max;i++){
		sum += i;
	}
	return sum;
}


int main(){
	cout << f(1,2) << endl;
	return 0;
}
