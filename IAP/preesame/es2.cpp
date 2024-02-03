#include<iostream>

using namespace std;

void table(int x){
	int i=1;
	int j=0;
	while(i<=x){
		while(j<=x){
			cout << i*j << " ";
			j++;
		}
		cout << endl;
		i++;
		j=1;
	}

}


int main(){

	table(10);

	return 0;
}
