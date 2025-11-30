#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace std;

int main() {

  // apro il file
  ifstream file;
  file.open("InterfPattern.txt",ios::in);
  
  // cerco il minimo e il massimo
  double value = 0;
  double min = 0;
  double max = 0;
  while(true){
    file >> value;                      
    if(value < min) min = value;
    if(value > max) max = value;
    if(file.eof() == true) break;    
  }
  cout << "Minimo:  " << min << endl;
  cout << "Massimo: " << max << endl;
  
  // riavvolgo il file
  file.clear();
  file.seekg(0, ios::beg);


  // creo l'istogramma e lo riempio
  TH1D *hist = new TH1D("hist","Titolo",600,min,max);
  double new_value = 0;
  while(true){
    file >> new_value;
    if(file.eof() == true) break;
    hist->Fill(new_value);    
  }

  double max_1 = 3.202e-3;
  double max_2 = 6.394e-3;
  double delta_max = max_2 - max_1;
  double d = 1.;
  double a = 3e-6;
  double lambda = ( delta_max * a ) / d;

  cout << "lambda:  " << lambda << endl;


  // =================== ROOT =================== //
  TApplication app("app",0,0);

  // creazione canvas
  TCanvas *c = new TCanvas("","",800,800);
  c->cd();

  // stampo l'istogramma
  hist->Draw("HIST");
  
  app.Run();
  // =================== ROOT =================== //
  
  
  

  return 0;
  }
