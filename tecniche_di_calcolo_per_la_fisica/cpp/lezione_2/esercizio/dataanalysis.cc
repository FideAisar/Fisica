#include "dataanalysis.h"

double wavg(double* cent_val, double* unc, int dim) {
  double num_media_pesata = 0.;
  double den_media_pesata = 0.;
  
  for(int i = 0; i < dim; i++) {
    double val = cent_val[i];  
    double peso = 1.0 / unc[i] * unc[i];  
    
    num_media_pesata += val * peso;  
    den_media_pesata += peso;
  }

  return num_media_pesata / den_media_pesata;
}
