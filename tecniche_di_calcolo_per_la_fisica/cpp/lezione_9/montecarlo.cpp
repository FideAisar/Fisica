#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"        
#include "TF2.h"        
#include "TLegend.h"   
#include "TRandom3.h"
#include <TFile.h>
using namespace std;


double xcube(double x){return x*x*x;}
double xsquare(double x){return x*x;}


double MonteCarloIntegration(double (*f)(double), int npoints, double a, double b) {

  TCanvas *c = new TCanvas("c","Monte Carlo Integration",1600,800);
  c->Divide(2,1);

  // calocolo il massimo della funzione nell'intervallo
  vector<double> domain;
  vector<double> function;
  
  int f_points = 1000;
  double domain_len = b - a;
  double dx = domain_len / f_points;

  for(int i=0; i<1000; i++) {
    double x = a + i*dx;
    domain.push_back(x);
    function.push_back(f(x));
  }


  // grafico la funzone
  c->cd(1);
  TGraph *f_plot = new TGraph(domain.size(), domain.data(), function.data());
  f_plot->SetLineColor(kBlack);
  f_plot->SetLineWidth(3);
  f_plot->SetTitle("f(x)");
  f_plot->Draw("AL");


  double fmax=function[0];
  double fmin=function[0];
  for(int i=0; i<1000; i++) {
    if(function[i] > fmax){
      fmax = function[i];
    }
    if(function[i] < fmin) {
      fmin = function[i];
    }
  }
  

  // colpi casuali nel rettangolo costruito  
  vector<double> x_hits, y_hits;
  vector<double> x_misses, y_misses;

  TRandom3 rand(42);
  int hits = 0;

  for(int i=0; i<npoints; i++) {
    double x = a + (b-a)*rand.Rndm();
    double y = fmin + (fmax - fmin)*rand.Rndm();

    double fx = f(x);
    
    if(fx >= 0) {
      if(y <= fx && y >= 0) {
        hits++;
        x_hits.push_back(x);
        y_hits.push_back(y);
      } else {
        x_misses.push_back(x);  
        y_misses.push_back(y);
      }
    } 
    else {
      if(y >= fx && y <= 0) {
        hits--;
        x_hits.push_back(x);
        y_hits.push_back(y);
      } else {
        x_misses.push_back(x);  
        y_misses.push_back(y);
      }
    }  
  }

  c->cd(2);
  TGraph *hits_plot = new TGraph(hits, x_hits.data(), y_hits.data());
  TGraph *misses_plot = new TGraph(npoints - hits, x_misses.data(), y_misses.data());
  
  hits_plot->SetMarkerColor(kRed);
  hits_plot->SetMarkerStyle(8);
  hits_plot->SetMarkerSize(0.7);
  
  misses_plot->SetMarkerColor(kBlue);
  misses_plot->SetMarkerStyle(8);
  misses_plot->SetMarkerSize(0.7);

  misses_plot->Draw("P SAME");
  hits_plot->Draw("P SAME");

  // conclusione
  double areatot = (b-a)*(fmax-fmin);
  double integral = areatot * hits / npoints;
  
  
  return integral;
}


int main() {
  
  TApplication app("app",0,0);
  
  double a = -1;
  double b = 1;
  int npoints = 10000;

  double integral = MonteCarloIntegration(xsquare, npoints, a, b);
  cout << "Integrale di x^2 tra " << a << " e " << b << ":  " << integral << endl;
  
  app.Run();

  return 0;
}
