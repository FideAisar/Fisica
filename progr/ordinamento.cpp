#include <iostream>
using namespace std;

						// ==FUNZIONI== //
void read_array(int v[], int n){		// read array
	for(int i=0;i<n;i++)
		cin >> v[i];
}

void print_array(const int v[], int n){		// print array
	for(int i=0;i<n;i++){
		cout << v[i] << " ";
	}
}

void swap(int &a, int &b){			// swap	
	int t=a;
	a=b;
	b=t;
}

int search_min(int v[], int first, int last){	// search minimum
	int pos=first;
	for(int k=first;k<=last;k++)
		if(v[k]<v[pos]) pos = k;
	return pos;
}

	int p;
	for(int i=0;i<=n-2;i++){
		p = search_min(v,i,n-1);
		swap(v[i],v[p]);
	}
}

						// ==MAIN== //
int main(){
	int x[10];
	read_array(x,10);
	selection(x,10);

	print_array(x,10);
	return 0;
}
