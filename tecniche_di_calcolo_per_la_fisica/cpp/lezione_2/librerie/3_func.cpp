#include <iostream>
#include "mylibrary.h"
using namespace std;


int main ()
  { 
    
    double log;
    cout << "Logaritmo naturale di: ";
    cin >> log;
    cout << ln(log,1000000) << endl;


    return 0 ;
  }  
