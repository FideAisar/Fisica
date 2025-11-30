#include <fstream>
#include "TH1D.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace std;

int main(){
  
  // apro il file per scriverci dentro, poi lo chiudo
  ofstream out_file("numbers.txt") ;
  for(int i=0; i<10; i++) out_file << i*3 << endl;
  out_file.close();
  

  double value = 0;
  
  
  // inizializzo un istogramma di double
  // - TH1D è la classe histogram di double
  // - Creo gli oggetti come puntatori perché ROOT 
  //   utilizza un sistema di gestione della memoria 
  //   che funziona meglio con oggetti allocati dinamicamente sull'heap
  //
  //                                                      bins  min max
  TH1D *hist = new TH1D("hist","Histogram;x-value;counts",50,   0,  50);
  

  // apro il file per leggerlo
  ifstream in_file;
  in_file.open("numbers.txt",ios::in);
  
  // itero sulle righe finché non raggiunto end of file e riempio
  // l'istogramma con i valori che trovo nelle righe.
  // >> operator reads a number (before a space, a tab or a to-the-line)
  while(true){
    file >> value;
    if(file.eof() == true) break;
    hist->Fill(value);
  }
  
  // If you want to activate in the same way the interactive mode, 
  // you need to use an object of the class TApplication.
  // Between the declaration of the object and the call of the 
  // method Run() you need to insert all the code that you want 
  // to use interactively
  TApplication app("app",0,0);
  TCanvas *c = new TCanvas("","",800,800);

  // l'operatore -> è usato per accedere a membri (funzioni o variabili) 
  // di un oggetto tramite un puntatore.
  //
  // Se hai un OGGETTO diretto, usi il punto (.)
  // TCanvas can;           // oggetto
  // can.cd();             // usa il punto

  // Se hai un PUNTATORE a un oggetto, usi la freccia (->)
  // TCanvas *can = new TCanvas();  // puntatore
  // can->cd();                     // usa la freccia
  //
  // Seleziona questa canvas come corrente
  c->cd();
  hist->Draw("lego");
  app.Run();
  
  return 0;

}
