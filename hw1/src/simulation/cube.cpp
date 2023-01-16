#include "cube.h"

#include "Eigen/Dense"

#include "../util/helper.h"
namespace simulation {
constexpr float g_cdK = 2500.0f;
constexpr float g_cdD = 50.0f;

Cube::Cube()
    : particleNumPerEdge(10),
      cubeLength(2.0),
      initialPosition(Eigen::Vector3f(0.0, 0.0, 0.0)),
      springCoefStruct(g_cdK),
      springCoefShear(g_cdK),
      springCoefBending(g_cdK),
      damperCoefStruct(g_cdD),
      damperCoefShear(g_cdD),
      damperCoefBending(g_cdD) {
    particleNumPerFace = particleNumPerEdge * particleNumPerEdge;
    initializeParticle();
    initializeSpring();
}

Cube::Cube(const Eigen::Vector3f &a_kInitPos, const float cubeLength, const int numAtEdge, const float dSpringCoef,
           const float dDamperCoef)
    : particleNumPerEdge(numAtEdge),
      cubeLength(cubeLength),
      initialPosition(a_kInitPos),
      springCoefStruct(dSpringCoef),
      springCoefShear(dSpringCoef),
      springCoefBending(dSpringCoef),
      damperCoefStruct(dDamperCoef),
      damperCoefShear(dDamperCoef),
      damperCoefBending(dDamperCoef) {
    particleNumPerFace = numAtEdge * numAtEdge;
    initializeParticle();
    initializeSpring();
}

int Cube::getParticleNum() const { return static_cast<int>(particles.size()); }

int Cube::getSpringNum() const { return static_cast<int>(springs.size()); }

int Cube::getNumAtEdge() const { return particleNumPerEdge; }

unsigned int Cube::getPointMap(const int a_ciSide, const int a_ciI, const int a_ciJ) {
    int r = -1;

    switch (a_ciSide) {
        case 1:  // [a_ciI][a_ciJ][0] bottom face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ;
            break;
        case 6:  // [a_ciI][a_ciJ][9] top face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ + particleNumPerEdge - 1;
            break;
        case 2:  // [a_ciI][0][a_ciJ] front face
            r = particleNumPerFace * a_ciI + a_ciJ;
            break;
        case 5:  // [a_ciI][9][a_ciJ] back face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * (particleNumPerEdge - 1) + a_ciJ;
            break;
        case 3:  // [0][a_ciI][a_ciJ] left face
            r = particleNumPerEdge * a_ciI + a_ciJ;
            break;
        case 4:  // [9][a_ciI][a_ciJ] ra_ciIght face
            r = particleNumPerFace * (particleNumPerEdge - 1) + particleNumPerEdge * a_ciI + a_ciJ;
            break;
    }

    return r;
}

Particle &Cube::getParticle(int particleIdx) { return particles[particleIdx]; }

std::vector<Particle> *Cube::getParticlePointer() { return &particles; }

Spring &Cube::getSpring(int springIdx) { return springs[springIdx]; }

void Cube::setSpringCoef(const float springCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        springCoefStruct = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        springCoefShear = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        springCoefBending = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::BENDING);
    }
}

void Cube::setDamperCoef(const float damperCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        damperCoefStruct = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        damperCoefShear = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        damperCoefBending = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::BENDING);
    }
}

void Cube::resetCube(const Eigen::Vector3f &offset, const float &rotate) {
    float dTheta = util::radians(rotate);  //  change angle from degree to
                                           //  radian

    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        int i = uiI / particleNumPerFace;
        int j = (uiI / particleNumPerEdge) % particleNumPerEdge;
        int k = uiI % particleNumPerEdge;
        float offset_x = (float)((i - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
        float offset_y = (float)((j - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
        float offset_z = (float)((k - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));

        Eigen::Vector3f RotateVec(offset_x, offset_y,
                                  offset_z);  //  vector from center of cube to the particle

        Eigen::AngleAxis<float> rotation(dTheta, Eigen::Vector3f(1.0f, 0.0f, 1.0f).normalized());

        RotateVec = rotation * RotateVec;

        particles[uiI].setPosition(initialPosition + offset + RotateVec);
        particles[uiI].setForce(Eigen::Vector3f::Zero());
        particles[uiI].setVelocity(Eigen::Vector3f::Zero());
    }
}

void Cube::addForceField(const Eigen::Vector3f &force) {
    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        particles[uiI].setAcceleration(force);
    }
}

void Cube::computeInternalForce() {
    // TODO
    for (int i = 0; i < springs.size(); i++) {

        int A_id = springs[i].getSpringStartID();
        int B_id = springs[i].getSpringEndID();
        Eigen::Vector3f A_position = particles[A_id].getPosition();
        Eigen::Vector3f B_position = particles[B_id].getPosition();
        Eigen::Vector3f A_velocity = particles[A_id].getVelocity();
        Eigen::Vector3f B_velocity = particles[B_id].getVelocity();

        Eigen::Vector3f spring_force =
            computeSpringForce(A_position, B_position, springs[i].getSpringCoef(), springs[i].getSpringRestLength());
        Eigen::Vector3f damper_force =
            computeDamperForce(A_position, B_position, A_velocity, B_velocity, springs[i].getDamperCoef());

        Eigen::Vector3f net_force = spring_force + damper_force;

        particles[A_id].addForce(net_force);
        particles[B_id].addForce(-net_force);


    }
}

Eigen::Vector3f Cube::computeSpringForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const float springCoef, const float restLength) {
    // TODO
    Eigen::Vector3f distance_vector = positionA - positionB;
    float distance = sqrt(distance_vector[0] * distance_vector[0] + distance_vector[1] * distance_vector[1] +
                             distance_vector[2] * distance_vector[2]);
    
    float deltaLength = distance - restLength;
    Eigen::Vector3f spring_force = -springCoef * deltaLength * distance_vector / distance;
    return spring_force;
}

Eigen::Vector3f Cube::computeDamperForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const Eigen::Vector3f &velocityA, const Eigen::Vector3f &velocityB,
                                         const float damperCoef) {
    // TODO

    Eigen::Vector3f distance_vector = positionA - positionB;
    float distance = sqrt(distance_vector[0] * distance_vector[0] + distance_vector[1] * distance_vector[1] +
                          distance_vector[2] * distance_vector[2]);
    Eigen::Vector3f relative_velocity = velocityA - velocityB;
    Eigen::Vector3f damper_force =
        -damperCoef * relative_velocity.dot(distance_vector) / distance * distance_vector / distance;
    return damper_force;
}

void Cube::initializeParticle() {
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
            //for (int k = 0; k < particleNumPerEdge; k++) {
                Particle Particle;
                float offset_x = (float)((i - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
                float offset_y = (float)((j - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
                float offset_z = (float)((k - particleNumPerEdge / 2) * cubeLength / (particleNumPerEdge - 1));
                Particle.setPosition(Eigen::Vector3f(initialPosition(0) + offset_x, initialPosition(1) + offset_y,
                                                     initialPosition(2) + offset_z));
                particles.push_back(Particle);

            }
        }
    }
}

void Cube::initializeSpring() {
    /*
    int iParticleID = 0;
    int iNeighborID = 0;
    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
    float length = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

    springs.push_back(
        Spring(iParticleID, iNeighborID, length, springCoefStruct, damperCoefStruct, Spring::SpringType::STRUCT));

    enum {
        X_id = 1,
        Y_id,  //  = 2
        Z_id   //  = 3
    };
    */


    // TODO

    // struct
    // z-direction
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge - 1; k++) {

                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = iParticleID + 1;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefStruct, damperCoefStruct,
                                         Spring::SpringType::STRUCT));


            }
        }
    }

    // y-direction
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int k = 0; k < particleNumPerEdge; k++) {
            for (int j = 0; j < particleNumPerEdge - 1; j++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = iParticleID + particleNumPerEdge;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefStruct, damperCoefStruct,
                                         Spring::SpringType::STRUCT));
            }
        }
    }

    // x-direction
    for (int j = 0; j < particleNumPerEdge; j++) {
        for (int k = 0; k < particleNumPerEdge; k++) {
            for (int i = 0; i < particleNumPerEdge - 1; i++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = iParticleID + particleNumPerFace;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefStruct, damperCoefStruct,
                                         Spring::SpringType::STRUCT));
            }
        }
    }



    // bend
    // z-direction
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge - 2; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = iParticleID + 2;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefBending, damperCoefBending,
                                         Spring::SpringType::BENDING));
            }
        }
    }
    
    // y-direction
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int k = 0; k < particleNumPerEdge; k++) {
            for (int j = 0; j < particleNumPerEdge - 2; j++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = iParticleID + 2 * particleNumPerEdge;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefBending, damperCoefBending,
                                         Spring::SpringType::BENDING));
            }
        }
    }

    // x-direction
    for (int j = 0; j < particleNumPerEdge; j++) {
        for (int k = 0; k < particleNumPerEdge; k++) {
            for (int i = 0; i < particleNumPerEdge - 2; i++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID = iParticleID + 2 * particleNumPerFace;
                Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
                Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
                Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
                float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefBending, damperCoefBending,
                                         Spring::SpringType::BENDING));
            }
        }
    }


    // shear
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {

                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID;
                Eigen::Vector3f SpringStartPos;
                Eigen::Vector3f SpringEndPos;
                Eigen::Vector3f Length;
                float absLength;

                if (k != 0) {
                    if (i != particleNumPerEdge - 1) {
                        // 1
                        iNeighborID = iParticleID + particleNumPerFace - 1;
                        SpringStartPos = particles[iParticleID].getPosition();
                        SpringEndPos = particles[iNeighborID].getPosition();
                        Length = SpringStartPos - SpringEndPos;
                        absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                        springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                                                 damperCoefShear, Spring::SpringType::SHEAR));
                        

                        // 2
                        if (j != 0) {

                            iNeighborID = iParticleID + particleNumPerFace - particleNumPerEdge - 1;
                            SpringStartPos = particles[iParticleID].getPosition();
                            SpringEndPos = particles[iNeighborID].getPosition();
                            Length = SpringStartPos - SpringEndPos;
                            absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                            springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                                                     damperCoefShear, Spring::SpringType::SHEAR));

                        
                        }

                        // 3
                        if (j != particleNumPerEdge - 1) {
                            iNeighborID = iParticleID + particleNumPerFace + particleNumPerEdge - 1;
                            SpringStartPos = particles[iParticleID].getPosition();
                            SpringEndPos = particles[iNeighborID].getPosition();
                            Length = SpringStartPos - SpringEndPos;
                            absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                            springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                                                     damperCoefShear, Spring::SpringType::SHEAR));
                        }

                    }

                    //4
                    if (j != particleNumPerEdge - 1) {
                        iNeighborID = iParticleID + particleNumPerEdge - 1;
                        SpringStartPos = particles[iParticleID].getPosition();
                        SpringEndPos = particles[iNeighborID].getPosition();
                        Length = SpringStartPos - SpringEndPos;
                        absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                        springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                                                 damperCoefShear, Spring::SpringType::SHEAR));
                    }
                }

                // 5
                if (i != particleNumPerEdge - 1 && j != particleNumPerEdge - 1) {
                    iNeighborID = iParticleID + particleNumPerFace + particleNumPerEdge;
                    SpringStartPos = particles[iParticleID].getPosition();
                    SpringEndPos = particles[iNeighborID].getPosition();
                    Length = SpringStartPos - SpringEndPos;
                    absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                }

                // 6
                if (i != particleNumPerEdge - 1 && j != 0) {
                    iNeighborID = iParticleID + particleNumPerFace - particleNumPerEdge;
                    SpringStartPos = particles[iParticleID].getPosition();
                    SpringEndPos = particles[iNeighborID].getPosition();
                    Length = SpringStartPos - SpringEndPos;
                    absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear, damperCoefShear,
                                             Spring::SpringType::SHEAR));
                }

                
                if (k != particleNumPerEdge - 1) {
                    if (i != particleNumPerEdge - 1) {

                        // 7
                        iNeighborID = iParticleID + particleNumPerFace + 1;
                        SpringStartPos = particles[iParticleID].getPosition();
                        SpringEndPos = particles[iNeighborID].getPosition();
                        Length = SpringStartPos - SpringEndPos;
                        absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                        springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                                                 damperCoefShear, Spring::SpringType::SHEAR));
                        

                        // 8
                        if (j != 0) {

                            iNeighborID = iParticleID + particleNumPerFace - particleNumPerEdge + 1;
                            SpringStartPos = particles[iParticleID].getPosition();
                            SpringEndPos = particles[iNeighborID].getPosition();
                            Length = SpringStartPos - SpringEndPos;
                            absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                            springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                                                     damperCoefShear, Spring::SpringType::SHEAR));

                        
                        }

                        // 9
                        if (j != particleNumPerEdge - 1) {
                            iNeighborID = iParticleID + particleNumPerFace + particleNumPerEdge + 1;
                            SpringStartPos = particles[iParticleID].getPosition();
                            SpringEndPos = particles[iNeighborID].getPosition();
                            Length = SpringStartPos - SpringEndPos;
                            absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                            springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                                                     damperCoefShear, Spring::SpringType::SHEAR));
                        }
                    }

                    // 10
                    if (j != particleNumPerEdge - 1) {
                        iNeighborID = iParticleID + particleNumPerEdge + 1;
                        SpringStartPos = particles[iParticleID].getPosition();
                        SpringEndPos = particles[iNeighborID].getPosition();
                        Length = SpringStartPos - SpringEndPos;
                        absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);

                        springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoefShear,
                                                 damperCoefShear, Spring::SpringType::SHEAR));
                    }
                }
                
            }
        }
    }

}

void Cube::updateSpringCoef(const float a_cdSpringCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setSpringCoef(a_cdSpringCoef);
        }
    }
}

void Cube::updateDamperCoef(const float a_cdDamperCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setDamperCoef(a_cdDamperCoef);
        }
    }
}
}  //  namespace simulation
