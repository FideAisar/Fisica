#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"        
#include "TLegend.h"   
#include "TGraph.h"
#include "TStyle.h"

using namespace std;

// ================ FUNZIONI ================ //

// f(x) = x^2
double x_square(double x) {
    return x*x;
}
vector<double> getDomainVector(double min, double max, int steps) {

  double domain_len = max - min;
  double dx = domain_len / steps;
  
  vector<double> x;
  for(int i=min; i<steps; i++) {
    x.push_back(min + i*dx);
  }

  return x;
}
vector<double> getFuncVector(vector<double> dom, double (*f)(double)) {
  
  vector<double> func;
  for(int i=0; i<dom.size(); i++) {
    func.push_back(f(dom[i]));
  }

  return func;
}




// ================ FOURIER ================ //

vector<double> fourier(double (*f)(double), double min, double max, int N, int steps) {

    double L = (max - min) / 2.0;  // L è la metà della lunghezza del dominio
    double domain_len = max - min;
    double dx = domain_len / steps;


    // calcolo  a0
    double a0 = 0;
    for(int i = 0; i <= steps; i++) {
        double x = min + i * dx;
        a0 += f(x) * dx;
    }
    a0 = a0 / (2 * L);


    // calocolo an e bn
    vector<double> an(N, 0.0);
    vector<double> bn(N, 0.0);
    
    for(int n = 1; n <= N; n++) {
        for(int i = 0; i <= steps; i++) {
            double x = min + i * dx;
            double func_val = f(x);
            an[n-1] += func_val * cos(n * M_PI * x / L) * dx;
            bn[n-1] += func_val * sin(n * M_PI * x / L) * dx;
        }
        an[n-1] /= L;
        bn[n-1] /= L;
    }


    // calcolo fourier
    vector<double> fourier_series;
    for(int i = 0; i <= steps; i++) {
        double x = min + i * dx;
        double fourier_val = a0;
        
        for(int n = 1; n <= N; n++) {
            fourier_val += an[n-1] * cos(n * M_PI * x / L) + bn[n-1] * sin(n * M_PI * x / L);
        }
        fourier_series.push_back(fourier_val);
    }


    return fourier_series;
}
void printFourier() {
  TApplication app("app", 0, 0);
  double min = -1.0;
  double max = 1.0;
  int N = 2;                // grado della serie 
  int steps = 1000;         // setp degli integrali
  double dx = (max - min) / steps;

  vector<double> domain = getDomainVector(min, max, steps);
  vector<double> funzione = getFuncVector(domain, x_square);
  vector<double> fourier_serie_2 = fourier(x_square, min, max, 2, steps);
  vector<double> fourier_serie_4 = fourier(x_square, min, max, 4, steps);
  vector<double> fourier_serie_8 = fourier(x_square, min, max, 8, steps);
  vector<double> fourier_serie_12 = fourier(x_square, min, max, 12, steps);
  
  
  // ================ ROOT ================ //
  TCanvas *c = new TCanvas("","",800,700);
  c->Divide(2,2); // Divide il canvas in 2x2 quadranti

  // Primo quadrante: funzione originale + troncamento n=2
  c->cd(1);
  TGraph *g1 = new TGraph(domain.size(), domain.data(), funzione.data());
  TGraph *g2 = new TGraph(domain.size(), domain.data(), fourier_serie_2.data());
  g1->SetLineColor(1);
  g1->SetLineWidth(3);
  g1->SetTitle("Troncamento n=2");
  g2->SetLineColor(2);
  g2->SetLineWidth(2);
  g1->Draw("AC");
  g2->Draw("C");

  // Secondo quadrante: funzione originale + troncamento n=4
  c->cd(2);
  TGraph *g1_2 = new TGraph(domain.size(), domain.data(), funzione.data());
  TGraph *g4 = new TGraph(domain.size(), domain.data(), fourier_serie_4.data());
  g1_2->SetLineColor(1);
  g1_2->SetLineWidth(3);
  g1_2->SetTitle("Troncamento n=4");
  g4->SetLineColor(3);
  g4->SetLineWidth(2);
  g1_2->Draw("AC");
  g4->Draw("C");

  // Terzo quadrante: funzione originale + troncamento n=8
  c->cd(3);
  TGraph *g1_3 = new TGraph(domain.size(), domain.data(), funzione.data());
  TGraph *g8 = new TGraph(domain.size(), domain.data(), fourier_serie_8.data());
  g1_3->SetLineColor(1);
  g1_3->SetLineWidth(3);
  g1_3->SetTitle("Troncamento n=8");
  g8->SetLineColor(4);
  g8->SetLineWidth(2);
  g1_3->Draw("AC");
  g8->Draw("C");

  // Quarto quadrante: funzione originale + troncamento n=12
  c->cd(4);
  TGraph *g1_4 = new TGraph(domain.size(), domain.data(), funzione.data());
  TGraph *g12 = new TGraph(domain.size(), domain.data(), fourier_serie_12.data());
  g1_4->SetLineColor(1);
  g1_4->SetLineWidth(3);
  g1_4->SetTitle("Troncamento n=12");
  g12->SetLineColor(5);
  g12->SetLineWidth(2);
  g1_4->Draw("AC");
  g12->Draw("C");

  app.Run();
}
void printDerivativeAndIntegral() {
  TApplication app("app", 0, 0);

  app.Run();
}



int main() {
  


  return 0;
}
