#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"        
#include "TF2.h"        
#include "TLegend.h"   
using namespace std;


void 3d_and_2d_visual(){
  
  // funzione
  TF2 *func2d1 = new TF2("func2d1","[0] * exp(-(x*x + y*y/[2])/[1])",-5,5,-5,5);
  func2d1->SetParameters(1.0, 2.0,0.3);
  func2d1->SetNpx(10000);

  // canvas
  TCanvas *can2 = new TCanvas("can2","can2",2000,1000);
  can2->Divide(2,1);
  can2->cd(1);
  func2d1->Draw("SURF2");
  can2->cd(2);
  func2d1->Draw("COLZ");
}

int main() {
  
  TApplication app("app",0,0);
  
  3d_and_2d_visual();

  app.Run();

  return 0;
}
