#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"        
#include "TF2.h"        
#include "TLegend.h"   
#include "TGraphErrors.h"
#include "TFitResult.h"
using namespace std;



int main() {
  TApplication app("app",0,0);


  double t[20] = { 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5 } ;
  double V[20] = { 5.09, 4.16, 3.53, 3.13, 2.72, 2.09, 1.82, 1.84, 1.34, 0.99, 1.03, 0.74, 0.64, 
            0.60, 0.54, 0.44, 0.55, 0.30, 0.32, 0.18 } ;
  
  // grafici
  double eV[20] = {0.1};
  TGraphErrors *g = new TGraphErrors(20, t, V, nullptr, eV);
  g->SetTitle("Capacitor Discharge; time [s]; voltage [V]");
  g->SetMarkerStyle(20);
  
  // i valori del fit sono salvati in r
  TF1 *fexp = new TF1("fexp", "[0] * exp(-x/[1])", 0, 10);
  fexp->SetParameters(4.0, 1.0);  // initial guesses
  TFitResultPtr r = g->Fit(fexp, "S");

  // frame
  TH1I *frame = new TH1I("frame","Capacitor Discharge; time [s]; voltage [V]",1,-0.5,10.5);
  frame->GetYaxis()->SetRangeUser(0,5.9);

  // canvas
  TCanvas *can = new TCanvas("can","can",1000,1000);
  can->cd();

  // draw
  frame->Draw();
  g->Draw("P");

  // grafico
  cout << "V0  = " << r->Parameter(0) << " ± "   << r->ParError(0) << endl;
  cout << "tau = " << r->Parameter(1) << " ± "   << r->ParError(1) << endl;
  cout << "chi2/ndf = " << r->Chi2() << "/" << r->Ndf() << "   Prob = " << r->Prob() << endl;


  app.Run();
  return 0;
}
