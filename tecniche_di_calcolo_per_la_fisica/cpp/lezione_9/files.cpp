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

  TApplication app("app",0,0);
  TCanvas *c = new TCanvas("","",1600,800);
  c->Divide(1,1);


  // istogramma riempito con random gaus
  c->cd(1);
  TH1D *hist = new TH1D("hist", "gaus histogram", 100, -5, 5);
  hist->FillRandom("gaus", 5000);
  hist->Draw();
  c->Update();
  
  // file
  TFile *file = new TFile("output.root", "RECREATE");
  hist->Write();
  file->Close();


  app.Run();
  

  delete c;

  return 0;
}
