#include <iostream>
using namespace std;



int main() {

  int a;
  cout << "Insert number: ";
  cin >> a;
  int b = a;

  do{
    
    cout << "value of a: " << b << endl;
    b++;
  
  }  while( b < b*2 );
  
  return 0;
}

