#include <iostream>
using namespace std;

void fun(int riga, int j=1,int m=0){
	int i;
	i = riga;
	while(j<=riga){
		while(riga <= i){
			m++;
			cout << m << " ";
			i--;
		}
		cout << endl;
		j++;
		i += j;
	}
}

int main(){
	cout << fun(5);
	return 0;
 }
