#include "ForceField.h"
#include <cmath>


// ============ COSTRUTTORI ============ //

ForceField::ForceField(): 
    m_g(0), 
    m_k(0), 
    m_center_x(0), 
    m_center_y(0), 
    m_is_gravity(true) 
{
}

ForceField::ForceField(float g_acc) : 
    m_g(g_acc),
    m_k(0),
    m_center_x(0), 
    m_center_y(0), 
    m_is_gravity(true) 
{
}

ForceField::ForceField(float k_cost, float center_x, float center_y) : 
    m_g(0),
    m_k(k_cost),
    m_center_x(center_x),
    m_center_y(center_y), 
    m_is_gravity(false) {
}



// ============ METODI ============ //

void ForceField::applyToParticle(particle2D& particle) const {
    if (m_is_gravity) {
        // a = -g verso il basso
        particle.setAcceleration(0, -m_g);
    } else {
        // a = -k/m * (x - x0)
        float x, y;
        particle.getPosition(x, y);
        float mass = particle.getMass();
        
        float ax = -m_spring_constant / mass * (x - m_center_x);
        float ay = -m_spring_constant / mass * (y - m_center_y);
        
        particle.setAcceleration(ax, ay);
    }
}
