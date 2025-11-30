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
using namespace std;

int main() {
  // particella 1
  vector<double> p1_init_position = {10.,0.,0.};
  vector<double> p1_init_velocity = {-1.,0.,0.};
  vector<double> p1_init_acceleration= {0.,0.,0.};
  double p1_mass = 2;

  particle2D p1(p1_init_position, p1_init_velocity, p1_init_acceleration, p1_mass);

  vector<double> p2_init_position = {0.,0.,0.};
  vector<double> p2_init_velocity = {1.,0.,0.};
  vector<double> p2_init_acceleration= {0.,0.,0.};
  double p2_mass = 10;  

  particle2D p2(p2_init_position, p2_init_velocity, p2_init_acceleration, p2_mass);  // CORREGGI: passa p2_init

  // simulation parameters
  double total_time = 10;
  double dt = 0.02;
  double g = 0.;
  double k = 0.;
  double x0 = 0.;
  double y0 = 0.;
  double z0 = 0.;
  
  vector<double> B = {0., 0., 0.};
  vector<double> E = {0., 0., 0.};
  vector<double> custom_force = {0, 0., 0.};
  double scattering_distance = 0.1;

  particle2D::simulateTwo(p1,p2,total_time,dt,g,k,x0,y0,z0,B,E,custom_force,scattering_distance);
  particle2D::plotTwoParticles3D(p1, p2, total_time, dt);

  return 0;
}
