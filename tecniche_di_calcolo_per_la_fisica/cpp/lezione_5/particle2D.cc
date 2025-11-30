#include "particle2D.h"
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGraph2D.h"  
#include "TPolyLine3D.h"  
#include "TLegend.h"  // Aggiungi questo
using namespace std;

// Costruttore
particle2D::particle2D( vector<double> init_position, 
                        vector<double> init_velocity,
                        vector<double> init_acceleration,
                        double init_mass):
  position(init_position), 
  velocity(init_velocity), 
  acceleration(init_acceleration), 
  mass(init_mass)
{
}

// Operatori
particle2D& particle2D::operator=(const particle2D& orig) {
  if (this != &orig) {
    for(int i=0; i<3; i++) {
      position.at(i) = orig.position.at(i);
      velocity.at(i) = orig.velocity.at(i);
      acceleration.at(i) = orig.acceleration.at(i);
    } 
  }
  return *this;
}

particle2D particle2D::operator+(const particle2D& merger) const {
    double total_mass = mass + merger.mass;
    vector<double> new_position(3);
    vector<double> new_velocity(3);
    vector<double> new_acceleration = {0., 0., 0.}; 

    for(int i=0; i<3; i++) {
      double new_x = (position.at(i) * mass + merger.position.at(i) * merger.mass) / total_mass;
      double new_v = (velocity.at(i) * mass + merger.velocity.at(i) * merger.mass) / total_mass;
      new_position[i] = new_x;
      new_velocity[i] = new_v;
    }
    return particle2D(new_position, new_velocity, new_acceleration, total_mass);
}
// Setters
void particle2D::setPosition(const vector<double>& new_pos) {
    position = new_pos;
}

void particle2D::setVelocity(const vector<double>& new_vel) {
    velocity = new_vel;
}

void particle2D::setGravity(float g) {
  force[2] += g;
}

void particle2D::setElastic(double k, double x0, double y0, double z0) {
    force[0] += -k*(position[0] - x0);
    force[1] += -k*(position[1] - y0);
    force[2] += -k*(position[2] - z0);
}

void particle2D::setCustomForce(const vector<double>& new_force) {
  for(int i=0; i<3; i++){
    force[i] += new_force[i];
  }
}

void particle2D::setLorentz(const vector<double>& B, const vector<double>& E) {
  vector<double> electric_force = {E[0], E[1], E[2]};
  
  vector<double> magnetic_force = {
      (velocity[1] * B[2] - velocity[2] * B[1]),
      (velocity[2] * B[0] - velocity[0] * B[2]),
      (velocity[0] * B[1] - velocity[1] * B[0])
  };
  
  force[0] += electric_force[0] + magnetic_force[0];
  force[1] += electric_force[1] + magnetic_force[1];
  force[2] += electric_force[2] + magnetic_force[2];
}

// Getters
void particle2D::getPosition(vector<double>& pos) const { 
    pos = position; 
}

void particle2D::getVelocity(vector<double>& vel) const { 
    vel = velocity; 
}

// Getters per dati stored - RIMUOVI 'const' dalla fine
const vector<double>& particle2D::getStoredX() const { 
    return stored_position_x; 
}

const vector<double>& particle2D::getStoredY() const{ 
    return stored_position_y; 
}

const vector<double>& particle2D::getStoredZ() const { 
    return stored_position_z; 
}

// Update
void particle2D::update(double dt) {
  for(int i=0; i<3; i++) {
    acceleration[i] = force[i] / mass;
    velocity[i] += acceleration[i] * dt;
    position[i] += velocity[i] * dt;
  }
}

// Reset
void particle2D::resetForces() {
  force = {0.,0.,0.};
}

// State
void particle2D::storeState() {
  stored_position_x.push_back(position[0]);
  stored_position_y.push_back(position[1]);
  stored_position_z.push_back(position[2]);
  stored_velocity_x.push_back(velocity[0]);
  stored_velocity_y.push_back(velocity[1]);
  stored_velocity_z.push_back(velocity[2]);
}

void particle2D::printState() {
    cout << "Pos  " << position[0] << "      " << position[1] << "      " << position[2] << endl;
}

double particle2D::distanceTo(const particle2D & other) const {
  double dx = position[0] - other.position[0];
  double dy = position[1] - other.position[1];
  double dz = position[2] - other.position[2];
  return sqrt(dx*dx + dy*dy + dz*dz);
}



// Elastic Scattering - Rendi statico
double particle2D::randomAngle() {
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0.0, 360.0);
  return dis(gen);
}

void particle2D::elasticScattering(particle2D& p1, particle2D& p2) {
  double theta = randomAngle() * M_PI / 180.0;
  double phi = randomAngle() * M_PI / 180.0;
  vector<double> versor = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};

  vector<double> p1_velocity;
  vector<double> p2_velocity;
  p1.getVelocity(p1_velocity);
  p2.getVelocity(p2_velocity);

  vector<double> cm_velocity;
  particle2D cm = p1 + p2;
  cm.getVelocity(cm_velocity);

  vector<double> p1_cm_velocity(3);
  vector<double> p2_cm_velocity(3);
  vector<double> p1_new_velocity(3);
  vector<double> p2_new_velocity(3);

  for(int i=0; i<3; i++) {
    p1_cm_velocity[i] = p1_velocity[i] - cm_velocity[i];
    p2_cm_velocity[i] = p2_velocity[i] - cm_velocity[i];
  }
  
  double p1_speed_cm = sqrt(p1_cm_velocity[0]*p1_cm_velocity[0] + 
                           p1_cm_velocity[1]*p1_cm_velocity[1] + 
                           p1_cm_velocity[2]*p1_cm_velocity[2]);
  double p2_speed_cm = sqrt(p2_cm_velocity[0]*p2_cm_velocity[0] + 
                           p2_cm_velocity[1]*p2_cm_velocity[1] + 
                           p2_cm_velocity[2]*p2_cm_velocity[2]);
  
  for(int i=0; i<3; i++) {
    p1_new_velocity[i] = p1_speed_cm * versor[i] + cm_velocity[i];
    p2_new_velocity[i] = -p2_speed_cm * versor[i] + cm_velocity[i];
  }
  
  p1.setVelocity(p1_new_velocity);
  p2.setVelocity(p2_new_velocity);
}





// Simulazione di una particella
void particle2D::simulateOne(double total_time, double dt, double g, double k, 
                          double x0, double y0, double z0,
                          const vector<double>& B, const vector<double>& E,
                          const vector<double>& custom_force) {
  cout << " == INIZIO SIMULAZIONE == " << endl;
  cout << " ======================== " << endl;
  cout << " Tempo totale: " << total_time << "s" << endl;
  cout << " dt:           " << dt << "s" << endl;
  cout << " ======================== " << endl;
  printState();

  int total_steps = total_time / dt;

  for(int step=0; step < total_steps; step++) {
    resetForces();
    setGravity(g);
    setElastic(k,x0,y0,z0);
    setLorentz(B, E);  
    setCustomForce(custom_force);

    update(dt);
    storeState();
    
    if(step % 20 == 0 || step == total_steps - 1) {
      printState();
    }
  }
}

// Simulazione due particelle
void particle2D::simulateTwo(particle2D& p1, particle2D& p2, 
                            double total_time, double dt, double g, double k,
                            double x0, double y0, double z0,
                            const vector<double>& B, const vector<double>& E,
                            const vector<double>& custom_force,
                            double scattering_distance) {
  cout << " == INIZIO SIMULAZIONE DUE PARTICELLE == " << endl;
  cout << " ======================================= " << endl;
  cout << " Tempo totale: " << total_time << "s" << endl;
  cout << " dt:           " << dt << "s" << endl;
  cout << " Distanza scattering: " << scattering_distance << "m" << endl;
  cout << " ======================================= " << endl;
  
  cout << "Stato iniziale:" << endl;
  cout << "P1: "; p1.printState();
  cout << "P2: "; p2.printState();
  cout << "Distanza iniziale: " << p1.distanceTo(p2) << "m" << endl;
  cout << "------------------------" << endl;

  int total_steps = total_time / dt;
  int scattering_count = 0;
  bool is_scattering_possible = 1;

  p1.storeState();
  p2.storeState();

  for(int step=0; step < total_steps; step++) {
    double current_time = step * dt;
    
    p1.resetForces();
    p2.resetForces();
    
    p1.setGravity(g);
    p2.setGravity(g);
    p1.setElastic(k,x0,y0,z0);
    p2.setElastic(k,x0,y0,z0);
    p1.setLorentz(B, E);  
    p2.setLorentz(B, E);
    p1.setCustomForce(custom_force);
    p2.setCustomForce(custom_force);
   
    if(is_scattering_possible && p1.distanceTo(p2) < scattering_distance) {
        scattering_count++;
        cout << "SCATTERING #" << scattering_count << " al tempo " << current_time << "s" << endl;
        particle2D::elasticScattering(p1, p2);  // Chiama come statico
        cout << "Nuove velocità dopo scattering:" << endl;
        cout << "P1: "; p1.printState();
        cout << "P2: "; p2.printState();
        cout << "------------------------" << endl;
        is_scattering_possible = 0;
    }
    
    p1.update(dt);
    p2.update(dt);
    
    p1.storeState();
    p2.storeState();

    if(step % 1000 == 0 || step == total_steps - 1) {
      cout << "Tempo: " << current_time << "s" << endl;
      cout << "P1: "; p1.printState();
      cout << "P2: "; p2.printState();
      cout << "Distanza: " << p1.distanceTo(p2) << "m" << endl;
      cout << "Scattering totali: " << scattering_count << endl;
      cout << "------------------------" << endl;
    }
  }
  
  cout << " == FINE SIMULAZIONE == " << endl;
  cout << "Scattering totali: " << scattering_count << endl;
}




// ROOT //////////////////////////////////////////////////////////////////////////////////////////

void particle2D::plotZvsTime(double total_time, double dt) {
    // Crea l'applicazione ROOT
    TApplication app("app", 0, 0);
    
    // Crea un canvas per il grafico
    TCanvas *canvas = new TCanvas("canvas", "Posizione Z vs Tempo", 800, 600);
    
    int data_points = stored_position_z.size();
    vector<double> time(data_points);
    vector<double> pos_z(data_points);
    
    // Usa i dati reali dalla simulazione
    for(int i = 0; i < data_points; i++) {
        time[i] = i * dt;
        pos_z[i] = stored_position_z[i];
    }
    
    // Crea il grafico
    TGraph *graph = new TGraph(data_points, &time[0], &pos_z[0]);
    graph->SetTitle("Posizione Z vs Tempo;Tempo (s);Posizione Z (m)");
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.5);
    graph->SetMarkerColor(kBlue);
    
    // Disegna il grafico
    graph->Draw("APL");
    
    // Styling aggiuntivo
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();
    
    // Aggiorna il canvas
    canvas->Update();
    canvas->Modified();
    
    // Salva il grafico
    canvas->SaveAs("z_vs_tempo.png");
    
    cout << "Grafico salvato come 'z_vs_tempo.png'" << endl;
    cout << "Premi Ctrl+C per uscire..." << endl;
    
    // Esegui l'applicazione ROOT
    app.Run();
}


void particle2D::plot3DPosition(double total_time, double dt) {
    // Crea l'applicazione ROOT
    TApplication app("app", 0, 0);
    
    // Crea un canvas per il grafico 3D
    TCanvas *canvas = new TCanvas("canvas", "Traiettoria 3D della Particella", 1000, 800);
    
    // Prepara i dati per il grafico 3D
    int data_points = stored_position_x.size();
    
    cout << "Creazione grafico 3D con " << data_points << " punti..." << endl;
    
    // Crea un TGraph2D per la traiettoria 3D
    TGraph2D *graph3d = new TGraph2D();
    graph3d->SetTitle("Traiettoria 3D della Particella;Posizione X (m);Posizione Y (m);Posizione Z (m)");
    
    // Riempie il grafico con i dati della simulazione
    for(int i = 0; i < data_points; i++) {
        graph3d->SetPoint(i, stored_position_x[i], stored_position_y[i], stored_position_z[i]);
    }
    
    // Configura l'aspetto del grafico
    graph3d->SetLineColor(kBlue);
    graph3d->SetLineWidth(2);
    graph3d->SetMarkerStyle(20);
    graph3d->SetMarkerSize(0.3);
    graph3d->SetMarkerColor(kRed);
    
    // Disegna il grafico 3D
    graph3d->Draw("LINE P");
    
    // Styling semplificato (evita GetHistogram() che causa problemi)
    graph3d->GetXaxis()->SetTitleOffset(1.5);
    graph3d->GetYaxis()->SetTitleOffset(1.5);
    graph3d->GetZaxis()->SetTitleOffset(1.5);
    
    graph3d->GetXaxis()->CenterTitle();
    graph3d->GetYaxis()->CenterTitle();
    graph3d->GetZaxis()->CenterTitle();
    
    // Aggiorna il canvas
    canvas->Update();
    canvas->Modified();
    
    // Salva il grafico
    canvas->SaveAs("traiettoria_3d.png");
    
    cout << "Grafico 3D salvato come 'traiettoria_3d.png'" << endl;
    cout << "Punti tracciati: " << data_points << endl;
    cout << "Premi Ctrl+C per uscire..." << endl;
    
    // Esegui l'applicazione ROOT
    app.Run();
}


void particle2D::plotTwoParticles3D(particle2D& p1, particle2D& p2, double total_time, double dt) {
    // Crea l'applicazione ROOT
    int argc = 0;
    char** argv = nullptr;
    TApplication app("app", &argc, argv);
    
    // Prepara i dati per le due particelle
    const vector<double>& p1_x = p1.getStoredX();
    const vector<double>& p1_y = p1.getStoredY();
    const vector<double>& p1_z = p1.getStoredZ();
    
    const vector<double>& p2_x = p2.getStoredX();
    const vector<double>& p2_y = p2.getStoredY();
    const vector<double>& p2_z = p2.getStoredZ();
    
    int data_points = min(p1_x.size(), p2_x.size());
    
    cout << "Creazione grafico 3D con " << data_points << " punti per particella..." << endl;
    
    // Verifica che ci siano dati da plottare
    if(data_points == 0) {
        cout << "ERRORE: Nessun dato da plottare!" << endl;
        cout << "p1_x size: " << p1_x.size() << ", p2_x size: " << p2_x.size() << endl;
        return;
    }
    
    // Crea TGraph2D per la prima particella (rossa)
    TGraph2D *graph3d_p1 = new TGraph2D();
    graph3d_p1->SetTitle("Traiettorie 3D delle Due Particelle;Posizione X (m);Posizione Y (m);Posizione Z (m)");
    
    // Crea TGraph2D per la seconda particella (blu)
    TGraph2D *graph3d_p2 = new TGraph2D();
    
    // Riempie i grafici con i dati delle simulazioni
    for(int i = 0; i < data_points; i++) {
        graph3d_p1->SetPoint(i, p1_x[i], p1_y[i], p1_z[i]);
        graph3d_p2->SetPoint(i, p2_x[i], p2_y[i], p2_z[i]);
    }
    
    // Configura l'aspetto della prima particella (rossa)
    graph3d_p1->SetLineColor(kRed);
    graph3d_p1->SetLineWidth(2);
    graph3d_p1->SetMarkerStyle(20);
    graph3d_p1->SetMarkerSize(0.5);
    graph3d_p1->SetMarkerColor(kRed);
    
    // Configura l'aspetto della seconda particella (blu)
    graph3d_p2->SetLineColor(kBlue);
    graph3d_p2->SetLineWidth(2);
    graph3d_p2->SetMarkerStyle(20);
    graph3d_p2->SetMarkerSize(0.5);
    graph3d_p2->SetMarkerColor(kBlue);
    
    // Crea un canvas per il grafico 3D
    TCanvas *canvas = new TCanvas("canvas", "Traiettorie 3D delle Due Particelle", 1200, 800);
    
    // Disegna entrambi i grafici
    graph3d_p1->Draw("LINE P");
    graph3d_p2->Draw("LINE P SAME");
    
    // Aggiungi legenda
    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(graph3d_p1, "Particella 1", "lp");
    legend->AddEntry(graph3d_p2, "Particella 2", "lp");
    legend->Draw();
    
    // Styling degli assi
    graph3d_p1->GetXaxis()->SetTitleOffset(1.5);
    graph3d_p1->GetYaxis()->SetTitleOffset(1.5);
    graph3d_p1->GetZaxis()->SetTitleOffset(1.5);
    
    graph3d_p1->GetXaxis()->CenterTitle();
    graph3d_p1->GetYaxis()->CenterTitle();
    graph3d_p1->GetZaxis()->CenterTitle();
    
    // Aggiorna il canvas
    canvas->Update();
    canvas->Modified();
    
    // Salva il grafico
    canvas->SaveAs("traiettorie_3d_due_particelle.png");
    
    cout << "Grafico 3D salvato come 'traiettorie_3d_due_particelle.png'" << endl;
    cout << "Punti tracciati: " << data_points << " per particella" << endl;
    cout << "Premi Ctrl+C per chiudere la finestra..." << endl;
    
    // Esegui l'applicazione ROOT (questo tiene aperto il programma)
    app.Run();
}

// ROOT //////////////////////////////////////////////////////////////////////////////////////////
