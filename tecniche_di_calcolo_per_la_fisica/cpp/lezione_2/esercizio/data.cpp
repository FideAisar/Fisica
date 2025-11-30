#include <iostream>
#include "dataanalysis.h"
using namespace std;

// prova commmento
int main ()
  {

  int dim = 4;
  double cen_val[dim] = {1,2,3,4};
  double unc[dim] = {0.3,0.4,0.5,0.1};
  
  cout << wavg(cen_val,unc,dim) << endl;

  return 0 ;
  }  
