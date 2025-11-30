#include <iostream>
using namespace std;


// POTENZA
double pow(double base, int pow) {

  double num = 1.;
  for(int i=0; i<pow; i++) {
    num = num * base;
  }

  return num;
}



// LOGARITMO NATURALE 
double ln(double y) {
  double x = y - 1;
  return x - (pow(x,2)/2) + (pow(x,3)/3) - (pow(x,4)/4) + (pow(x,5)/5); 
}




int main ()
  { 
    
    double log;
    cout << "Logaritmo naturale di: ";
    cin >> log;
    cout << ln(log) << endl;


    return 0 ;
  }  
