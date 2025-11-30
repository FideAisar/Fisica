#include <iostream>
#include "particle2D.h"
using namespace std;

// ============ COSTRUTTORI ============ //

particle2D::particle2D(float x0, float y0, float vx0, float vy0, float ax0, float ay0, float m, int steps) {
    m_position = new float[2];
    m_velocity = new float[2];
    m_acceleration = new float[2];
    m_mass = new float;
    
    m_force= new float[2];
    
    m_steps = steps;
    m_positions = new double[m_steps];
    m_velocities = new double[m_steps];


    m_position[0] = x0;
    m_position[1] = y0;
    m_velocity[0] = vx0;
    m_velocity[1] = vy0;
    m_acceleration[0] = ax0;
    m_acceleration[1] = ay0;
    *m_mass = m;

    m_force[0] = 0;
    m_force[1] = 0;

}

// Copy constructor
particle2D::particle2D(const particle2D& other) {
    m_position = new float[2];
    m_velocity = new float[2];
    m_acceleration = new float[2];
    m_mass = new float;
    
    m_position[0] = other.m_position[0];
    m_position[1] = other.m_position[1];
    m_velocity[0] = other.m_velocity[0];
    m_velocity[1] = other.m_velocity[1];
    m_acceleration[0] = other.m_acceleration[0];
    m_acceleration[1] = other.m_acceleration[1];
    *m_mass = *other.m_mass;

    m_force[0] = other.m_force[0];
    m_force[1] = other.m_force[1];
}

particle2D::~particle2D() {
    delete[] m_position;
    delete[] m_velocity;
    delete[] m_acceleration;
    delete m_mass;
    delete[] m_force;
    delete[] m_positions;
    delete[] m_velocities;
}

// ============ OPERATORI ============ //

particle2D& particle2D::operator=(const particle2D& orig) {
    if (this != &orig) {
        m_position[0] = orig.m_position[0];
        m_position[1] = orig.m_position[1];
        m_velocity[0] = orig.m_velocity[0];
        m_velocity[1] = orig.m_velocity[1];
        m_acceleration[0] = orig.m_acceleration[0];
        m_acceleration[1] = orig.m_acceleration[1];
        *m_mass = *orig.m_mass;

        m_force[0] = orig.m_force[0];
        m_force[1] = orig.m_force[1];
    }
    return *this;
}

particle2D particle2D::operator+(const particle2D& merger) const {
    float total_mass = *m_mass + *merger.m_mass;
    
    // Calcolo centro di massa per posizione
    float new_x = (m_position[0] * (*m_mass) + merger.m_position[0] * (*merger.m_mass)) / total_mass;
    float new_y = (m_position[1] * (*m_mass) + merger.m_position[1] * (*merger.m_mass)) / total_mass;
    
    // Calcolo velocità risultante
    float new_vx = (*m_mass * m_velocity[0] + *merger.m_mass * merger.m_velocity[0]) / total_mass;
    float new_vy = (*m_mass * m_velocity[1] + *merger.m_mass * merger.m_velocity[1]) / total_mass;
    
    // Calcolo accelerazione risultante
    float new_ax = (*m_mass * m_acceleration[0] + *merger.m_mass * merger.m_acceleration[0]) / total_mass;
    float new_ay = (*m_mass * m_acceleration[1] + *merger.m_mass * merger.m_acceleration[1]) / total_mass;
    
    return particle2D(new_x, new_y, new_vx, new_vy, new_ax, new_ay, total_mass);
}

// ============ METODI ============ //

void particle2D::printState(int step) const {
    cout << m_position[0] << "      " << m_position[1] << endl;
}



// BASICS SETS
void particle2D::setPosition(float x, float y) {
    m_position[0] = x;
    m_position[1] = y;
}
void particle2D::setVelocity(float vx, float vy) {
    m_velocity[0] = vx;
    m_velocity[1] = vy;
}
void particle2D::setAcceleration(float ax, float ay) {
    m_acceleration[0] = ax;
    m_acceleration[1] = ay;
}
void particle2D::setMass(float mass) {
    *m_mass = mass;
}



// BASICS GETS
void particle2D::getPosition(float& x, float& y) const {
    x = m_position[0];
    y = m_position[1];
}
void particle2D::getVelocity(float& vx, float& vy) const {
    vx = m_velocity[0];
    vy = m_velocity[1];
}
void particle2D::getAcceleration(float& ax, float& ay) const {
    ax = m_acceleration[0];
    ay = m_acceleration[1];
}
float particle2D::getMass() const {
    return *m_mass;
}



// FORCES
void particle2D::setGravity(float g) {
    m_force[0] += 0;
    m_force[1] += - *m_mass * g;
    updateAcc();
}
void particle2D::setElastic(float k, float x0, float y0) {
    m_force[0] += -k*(m_position[0] - x0);
    m_force[1] += -k*(m_position[1] - y0);
    updateAcc();
}
void particle2D::setCustomForce(float fx, float fy) {
    m_force[0] += fx;
    m_force[1] += fy;
    updateAcc();
}




// UPDATE
void particle2D::updateAcc() {
    m_acceleration[0] = m_force[0] / (*m_mass);
    m_acceleration[1] = m_force[1] / (*m_mass);
}
void particle2D::update(float dt) {
    m_velocity[0] += m_acceleration[0] * dt;
    m_velocity[1] += m_acceleration[1] * dt;
    
    m_position[0] += m_velocity[0] * dt;
    m_position[1] += m_velocity[1] * dt;
}
 


// SIMULA MOVIMENTO
void particle2D::resetForces() {
    m_force[0] = 0;
    m_force[1] = 0;
    updateAcc();
}

void particle2D::storeState(int step) {
    if(step < m_steps) {
        m_positions[step * 2] = m_position[0];     
        m_positions[step * 2 + 1] = m_position[1]; 
        m_velocities[step * 2] = m_velocity[0];    
        m_velocities[step * 2 + 1] = m_velocity[1]; 
    }
}

void particle2D::simulate(float total_time, int num_steps) {
    float dt = total_time / num_steps;
    float g = 9.81;

    float k1 = 5;
    float x01 = 0;
    float y01 = 0;

    float k2 = 5;
    float x02 = -10;
    float y02 = 0;
    
    cout << "╔════════════════════════════════════════╗\n";
    cout << "║        INIZIO SIMULAZIONE              ║\n";
    cout << "╠════════════════════════════════════════╣\n";
    cout << "║ Tempo totale: " << total_time << "s                 ║\n";
    cout << "║ Numero di steps: " << num_steps << "               ║\n";
    cout << "║ Delta t: " << dt << "s                    ║\n";
    cout << "╚════════════════════════════════════════╝\n";
    
    cout << "\nSTATO INIZIALE:\n";
    printState();
    
    for(int step = 0; step < num_steps; step++) {
        resetForces();
        
        setGravity(g);
        setElastic(k1, x01, y01);
        setElastic(k2, x02, y02);
        
        update(dt);
        
        // Stampa ogni 20 steps o all'ultimo step
        if(step % 20== 0 || step == num_steps - 1) {
            printState(step);
        }
    }
    
    cout << "╔════════════════════════════════════════╗\n";
    cout << "║         FINE SIMULAZIONE               ║\n";
    cout << "╚════════════════════════════════════════╝\n";
}
