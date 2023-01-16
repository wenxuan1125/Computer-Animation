#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO
    // float deltaTime = particleSystem.deltaTime;
    for (int i = 0; i < particleSystem.getCubeCount(); i++) {
        Cube* cube = particleSystem.getCubePointer(i);

        for (int j = 0; j < cube->getParticleNum(); j++) {
            Particle& particle = cube->getParticle(j);

            // x=x0+vt
            particle.addPosition(particleSystem.deltaTime * particle.getVelocity());

            Eigen::Vector3f acceleration = particle.getAcceleration();
            // v=v0+at
            particle.addVelocity(particleSystem.deltaTime * acceleration);

            // clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO

    for (int i = 0; i < particleSystem.getCubeCount(); i++) {
        Cube* cube = particleSystem.getCubePointer(i);

        std::vector<Eigen::Vector3f> current_position, current_velocity;

        for (int j = 0; j < cube->getParticleNum(); j++) {
            Particle& particle = cube->getParticle(j);

            // retain the current position and velocity
            current_position.push_back(particle.getPosition());
            current_velocity.push_back(particle.getVelocity());

            // clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }

        // calculate the next point force
        particleSystem.computeCubeForce(*cube);

        // next point information
        for (int j = 0; j < cube->getParticleNum(); j++) {
            Particle& particle = cube->getParticle(j);

            particle.setPosition(current_position[j] + particleSystem.deltaTime * particle.getVelocity());
            particle.setVelocity(current_velocity[j] + particleSystem.deltaTime * particle.getAcceleration());

            // clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO
    // For midpoint euler, the deltaTime passed in is correct.
    // But this deltaTime is for a full step.
    // So you may need to adjust it before computing, but don't forget to restore original value.

    for (int i = 0; i < particleSystem.getCubeCount(); i++) {
        Cube* cube = particleSystem.getCubePointer(i);

        std::vector<Eigen::Vector3f> current_position, current_velocity;

        for (int j = 0; j < cube->getParticleNum(); j++) {
            Particle& particle = cube->getParticle(j);

            // retain the current position and velocity
            current_position.push_back(particle.getPosition());
            current_velocity.push_back(particle.getVelocity());

            // change the position to mid point position
            // change velocity to mid point velocity
            particle.addPosition(0.5 * particleSystem.deltaTime * particle.getVelocity());
            particle.addVelocity(0.5 * particleSystem.deltaTime * particle.getAcceleration());

            // clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }

        // calculate the mid point force
        particleSystem.computeCubeForce(*cube);

        // mid point information
        for (int j = 0; j < cube->getParticleNum(); j++) {
            Particle& particle = cube->getParticle(j);

            particle.setPosition(current_position[j] + particleSystem.deltaTime * particle.getVelocity());
            particle.setVelocity(current_velocity[j] + particleSystem.deltaTime * particle.getAcceleration());

            // clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }
}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };
    // TODO
    // StateStep struct is just a hint, you can use whatever you want.

    for (int i = 0; i < particleSystem.getCubeCount(); i++) {
        Cube* cube = particleSystem.getCubePointer(i);

        std::vector<Eigen::Vector3f> current_position, current_velocity;
        std::vector<StateStep> k1, k2, k3;

        for (int j = 0; j < cube->getParticleNum(); j++) {
            Particle& particle = cube->getParticle(j);

            // retain the current position and velocity
            current_position.push_back(particle.getPosition());
            current_velocity.push_back(particle.getVelocity());

            StateStep temp;
            temp.deltaPos = particleSystem.deltaTime * particle.getVelocity();
            temp.deltaVel = particleSystem.deltaTime * particle.getAcceleration();
            k1.push_back(temp);

            // change the position and the velocity to k2 state
            particle.addPosition(0.5 * particleSystem.deltaTime * particle.getVelocity());
            particle.addVelocity(0.5 * particleSystem.deltaTime * particle.getAcceleration());

            // clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }

        // calculate the k2 force
        particleSystem.computeCubeForce(*cube);

        for (int j = 0; j < cube->getParticleNum(); j++) {
            Particle& particle = cube->getParticle(j);

            StateStep temp;
            temp.deltaPos = particleSystem.deltaTime * particle.getVelocity();
            temp.deltaVel = particleSystem.deltaTime * particle.getAcceleration();
            k2.push_back(temp);

            // change the position and the velocity to k3 state
            particle.setPosition(current_position[j] + 0.5 * temp.deltaPos);
            particle.setVelocity(current_velocity[j] + 0.5 * temp.deltaVel);

            // clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }

        // calculate the k3 force
        particleSystem.computeCubeForce(*cube);

        for (int j = 0; j < cube->getParticleNum(); j++) {
            Particle& particle = cube->getParticle(j);

            StateStep temp;
            temp.deltaPos = particleSystem.deltaTime * particle.getVelocity();
            temp.deltaVel = particleSystem.deltaTime * particle.getAcceleration();
            k3.push_back(temp);

            // change the position and the velocity to k4 state
            particle.setPosition(current_position[j] + temp.deltaPos);
            particle.setVelocity(current_velocity[j] + temp.deltaVel);

            // clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }

        // calculate the k4 force
        particleSystem.computeCubeForce(*cube);

        for (int j = 0; j < cube->getParticleNum(); j++) {
            Particle& particle = cube->getParticle(j);

            StateStep k4, total;
            k4.deltaPos = particleSystem.deltaTime * particle.getVelocity();
            k4.deltaVel = particleSystem.deltaTime * particle.getAcceleration();

            total.deltaPos = (k1[j].deltaPos + 2 * k2[j].deltaPos + 2 * k3[j].deltaPos + k4.deltaPos) / 6;
            total.deltaVel = (k1[j].deltaVel + 2 * k2[j].deltaVel + 2 * k3[j].deltaVel + k4.deltaVel) / 6;

            // change the position and the velocity to the weighted result
            particle.setPosition(current_position[j] + total.deltaPos);
            particle.setVelocity(current_velocity[j] + total.deltaVel);

            // clear force
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }
}
}  // namespace simulation
