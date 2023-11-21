#include <iostream>
using namespace std;

int main(){
	int i=1,j=1,n;
	cout << "l: ";
	cin >> n;
	

	while(i <= n){
		while(j <= n){		 // finchè j <= n stampa un "*"
			cout << "*";
			j++;             // esce dal ciclo quando ha stampato n righe
		}
		i++;			 // incremento di 1 per cambiare riga
		j = i;			 // rientro nel while annidato (partendo da j = 1+1)
		cout << "" << endl;	 
	}

	return 0;
}
