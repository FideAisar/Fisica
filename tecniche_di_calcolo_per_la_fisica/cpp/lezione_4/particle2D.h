#ifndef PARTICLE2D_H
#define PARTICLE2D_H

class particle2D {
private:
    float* m_position;    // [x, y]
    float* m_velocity;    // [vx, vy]  
    float* m_acceleration; // [ax, ay]
    float* m_mass;        // puntatore a massa
    
    float* m_force;

    int m_steps;
    double* m_positions;
    double* m_velocities;

    
public:
    // ============ COSTRUTTORI ============ //
    particle2D(float x0, float y0, float vx0, float vy0, float ax0, float ay0, float m, int steps = 2000);
    particle2D(const particle2D& other); // Copy constructor
    ~particle2D();

    // ============ OPERATORI ============ //
    particle2D& operator=(const particle2D& orig);
    particle2D operator+(const particle2D& merger) const;

    // ============ METODI ============ //
    void printState(int step = -1) const;  // CAMBIA QUESTA LINEA
    void setPosition(float x, float y);
    void setVelocity(float vx, float vy);
    void setAcceleration(float ax, float ay);
    void setMass(float mass);

    void setGravity(float g);
    void setElastic(float k, float x0, float y0);
    void setCustomForce(float fx, float fy);

    void getPosition(float& x, float& y) const;
    void getVelocity(float& vx, float& vy) const;
    void getAcceleration(float& ax, float& ay) const;
    float getMass() const;

    void updateAcc();
    void update(float dt);

    void simulate(float total_time, int num_steps);
    void resetForces(); 
    void storeState(int step); 
    void printTrajectory() const; 
};

#endif
