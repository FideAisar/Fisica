#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include "particle2D.h"

class ForceField {
private:
    float m_g;     
    float m_k;      
    float m_center_x, m_center_y;  
    bool m_is_gravity;           
public:
    ForceField(); // Nulla
    ForceField(float g_acceleration); // Gravità
    ForceField(float center_x, float center_y, float spring_constant); // Armonico
    
    // Applica accelerazione a una particella
    void applyToParticle(particle2D& particle) const;
};

#endif


