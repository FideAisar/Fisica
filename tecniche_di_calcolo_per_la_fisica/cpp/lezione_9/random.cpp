#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"        
#include "TF2.h"        
#include "TLegend.h"   
#include "TRandom3.h"
using namespace std;



int main() {

  TRandom3 randGen(42);  // 42 is the seed value
  double r = randGen.Rndm();  // Uniform random number in [0, 1)
  double r10 = randGen.Rndm() * 10.0;  // Random number in [0, 10)
  double r_gauss = randGen.Gaus(0, 1);  // Gaussian with mean=0 and std=1
  double r_exp   = randGen.Exp(2);  // Exponential with mean=2
  double r_poiss = randGen.Poisson(0.5);  // Poisson with mean=0.5  

  // global pointer to random number generator
  // gRandom set to point to a TRandom3 variable
  gRandom = new TRandom3(42);

  // Define a custom distribution (exponential distribution)
  TF1* f1 = new TF1("f1", "exp(-x)", 0, 10);  // Exponential decay from 0 to 10

  // Set normalization to make the total probability equal to 1 (not really needed)
  f1->SetNormalized(1);

  // Generate a random number following this distribution
  double randomValue = f1->GetRandom();

  cout << "Random value from custom distribution: " << randomValue << endl;

  delete f1;

  cout << "Random number: " << r << endl;

  return 0;
}
