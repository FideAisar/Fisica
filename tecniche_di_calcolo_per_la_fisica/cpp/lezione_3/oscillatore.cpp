#include <iostream>
#include <cmath>
using namespace std;

//---------------------------------------------
// Function to compute acceleration (by value)
double acceleration(double k, double m, double x) {
  return -k * x / m;
}

// Function to update velocity (pass by reference)
void updateVelocity(double& v, double a, double dt) {
  v += a * dt;
  return;
}

// Function to update position (pass by pointer)
void updatePosition(double* x, double v, double dt) {
  *x += v * dt;
  return;
}

//---------------------------------------------
int main() {
  
  // --- Stack variables ---
  double k = 10.0;    // spring constant (N/m)
  double m = 2.0;     // mass (kg)
  double x0 = 0.1;    // initial displacement (m)
  double v0 = 0.0;    // initial velocity (m/s)
  double dt = 0.01;   // time step (s)
  int steps = 1000;   // number of time steps

  // --- Heap allocation for data storage ---
  double* positions = new double[steps];
  double* velocities = new double[steps];
  double* times = new double[steps];

  // initialize first values
  positions[0] = x0;
  velocities[0] = v0;
  times[0] = 0.0;

  // --- Simulation loop ---
  double x = x0;  // local variables (stack)
  double v = v0;

  for (int i = 1; i < steps; ++i) {
      double a = acceleration(k, m, x);   // compute acceleration

      updateVelocity(v, a, dt);           // pass by reference
      updatePosition(&x, v, dt);          // pass by pointer

      positions[i] = x;
      velocities[i] = v;
      times[i] = i * dt;

      cout << positions[i] << " " ;
  }

  
  // --- Cleanup (free heap memory) ---
  delete[] positions;
  delete[] velocities;
  delete[] times;

  positions = nullptr;
  velocities = nullptr;
  times = nullptr;
  
  return 0;
}
