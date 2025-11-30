#include <iostream>
#include "particle2D.h"
using namespace std;

int main() {
    float dt = 0.001;
    float g = 9.81;
    float k = 9.;
    float x0 = 0.;
    float y0 = 0.;
              

    //            x  y  vx  vy  ax  ay  m  
    particle2D p1(0, 0, 0,  2,  0,  0,  13);
    
    cout << "PARTICELLA 1:" << endl;
    p1.printState();  // Senza parametro

    cout << endl;
    cout << "NUOVA SIMULAZIONE:" << endl;
    
    // 20sec per 1000 steps
    p1.simulate(20, 1000);

    return 0;
}
