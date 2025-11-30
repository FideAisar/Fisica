
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"        
#include "TF2.h"        
#include "TLegend.h"   
using namespace std;


void signalextraction(){

  // setup histogram
  TH1D *hist = new TH1D("hist","Invariant mass distribution",100,1.3,5);
  
  // read from file
  double readvalue;
  ifstream infile("jpsimass.txt");
  while(infile >> readvalue){
      hist->Fill(readvalue);
  }
  
  // fit
  TF1 *modelfit = new TF1("modelfit","gaus(0)+gaus(3)+pol1(6)",1,5);
  modelfit->SetNpx(5000);
  modelfit->SetParameters(6000,3.1,0.14,1500,3.6,0.1,48200,-1621);
  modelfit->SetParLimits(4,3.4,3.7);

  for(int i=0; i<8; i++){
    hist->Fit(modelfit,"0");
  }

  // canvas
  TCanvas *can = new TCanvas("can","can",600,400);
  can->SetMargin(0.15,0.05,0.1,0.1);
  can->cd();
  hist->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}}");
  hist->GetYaxis()->SetTitle("counts");

  // draw
  hist->Draw("E");
  modelfit->Draw("same");
}


//void interpolation()




int main() {
  TApplication app("app",0,0);



  signalextraction();



  app.Run();
  return 0;
}
