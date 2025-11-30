#ifndef PARTICLE2D_H
#define PARTICLE2D_H

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
#include "TH2.h"
#include "TLegend.h"
using namespace std;

class particle2D {
  
  private:
    // particella
    vector<double> position;
    vector<double> velocity;
    vector<double> acceleration;
    double mass;

    // campo esterno
    vector<double> force = {0.,0.,0.};
    
    vector<double> stored_position_x;
    vector<double> stored_position_y;
    vector<double> stored_position_z;
    vector<double> stored_velocity_x;
    vector<double> stored_velocity_y;
    vector<double> stored_velocity_z;

  public:
    particle2D(vector<double> init_position = {0.,0.,0.}, 
               vector<double> init_velocity = {0.,0.,0.},
               vector<double> init_acceleration = {0.,0.,0.},
               double mass = 1.);

    
    particle2D& operator=(const particle2D& orig);
    particle2D operator+(const particle2D& merger) const;

    void setPosition(const vector<double>& new_position);
    void setVelocity(const vector<double>& new_velocity);
    void setGravity(float g=9.81);
    void setElastic(double k, double x0, double y0, double z0);
    void setCustomForce(const vector<double>& new_force);
    void setLorentz(const vector<double>& B, const vector<double>& E = {0,0,0});

    void getPosition(vector<double>& pos) const;
    void getVelocity(vector<double>& vel) const;
    
    const vector<double>& getStoredX() const;
    const vector<double>& getStoredY() const;
    const vector<double>& getStoredZ() const;
    void update(double dt);
    void resetForces();

    void storeState();
    void printState();
    double distanceTo(const particle2D & other) const;

    void simulateOne(double total_time, double dt, double g, double k, 
                     double x0, double y0, double z0, 
                     const vector<double>& B, const vector<double>& E = {0,0,0},   
                     const vector<double>& custom_force = {0,0,0});  
    
    static void simulateTwo(particle2D& p1, particle2D& p2, 
                           double total_time, double dt, double g, double k,
                           double x0, double y0, double z0,
                           const vector<double>& B, const vector<double>& E,
                           const vector<double>& custom_force,
                           double scattering_distance = 1.0);
                                                                  
    static double randomAngle();
    static void elasticScattering(particle2D& p1, particle2D& p2);

    void plotZvsTime(double total_time, double dt);
    void plot3DPosition(double total_time, double dt);
    static void plotTwoParticles3D(particle2D& p1, particle2D& p2, double total_time, double dt);
    
    // Aggiungi funzione di debug
    void printStoredData();
};

#endif
