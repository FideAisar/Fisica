#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"        
#include "TF2.h"        
#include "TLegend.h"   
#include "TRandom3.h"
#include <TFile.h>
using namespace std;



int main() {
  
  // numeri random e distribuzioni
  TRandom3 randGen(42);

  TF1* f1 = new TF1("f1", "exp(-x)", 0, 10);
  f1->SetNormalized(1);

  TF1* f2 = new TF1("f2", "exp(-0.2*(x-4)^2)", -5, 15);
  f1->SetNormalized(1);



  TApplication app("app",0,0);
  TCanvas *c = new TCanvas("","",1600,800);
  c->Divide(3,1);


  // istogramma riempito con random gaus
  c->cd(1);
  TH1D* hist = new TH1D("hist", "Gaus Histogram", 100, -5, 5);
  for (int i = 0; i < 10000; ++i) {
    hist->Fill(randGen.Gaus(0, 1));  
  }
  hist->Draw();

  // istogramma riempito con exp 1
  c->cd(2);
  TH1D* hist_exp1 = new TH1D("hist", "Exp Histogram", 100, 0, 10);
  for (int i = 0; i < 10000; ++i) {
    hist_exp1->Fill(f1->GetRandom());  
  }
  hist_exp1->Draw();

  // istogramma riempito con exp 2
  c->cd(3);
  TH1D* hist_exp2 = new TH1D("hist", "Exp Histogram", 100, 0, 10);
  for (int i = 0; i < 10000; ++i) {
    hist_exp2->Fill(f2->GetRandom());  
  }
  hist_exp2->Draw();


  app.Run();

  delete hist;
  delete hist_exp1;
  delete hist_exp2;

  return 0;
}
