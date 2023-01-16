#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Sphere:
            return std::make_unique<SphereTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        case simulation::TerrainType::TiltedPlane:
            return std::make_unique<TiltedPlaneTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    // TODO
    for (int i = 0; i < cube.getParticleNum(); i++) {
        Particle& particle = cube.getParticle(i);

        Eigen::Vector3f relative_position = particle.getPosition() - position;

        // if the particle and the wall is close enough and if the particle is moving forward to the wall
        // collision happens
        if (normal.dot(relative_position) < eEPSILON && normal.dot(particle.getVelocity()) < 0) {

            Eigen::Vector3f normal_velocity = normal.dot(particle.getVelocity()) * normal;
            Eigen::Vector3f tangent_velocity = particle.getVelocity() - normal_velocity;

            particle.setVelocity(-coefResist * normal_velocity + tangent_velocity);

            // if a force pushes the particle into the wall 
            // exert a contact force (resist + friction) to resist it
            if (normal.dot(particle.getForce()) < 0) {

                Eigen::Vector3f resist = -normal.dot(particle.getForce()) * normal;
                Eigen::Vector3f friction =
                    coefFriction * normal.dot(particle.getForce()) * tangent_velocity / tangent_velocity.norm();
                particle.addForce(resist + friction);

            }
        }
        
    }
    

}

// SphereTerrain //

SphereTerrain::SphereTerrain() { modelMatrix = util::translate(position) * util::scale(radius, radius, radius); }

TerrainType SphereTerrain::getType() { return TerrainType::Sphere; }

void SphereTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    // TODO

    for (int i = 0; i < cube.getParticleNum(); i++) {
        Particle& particle = cube.getParticle(i);

        Eigen::Vector3f relative_position = particle.getPosition() - position;
        Eigen::Vector3f unit = relative_position / relative_position.norm();

        // if the particle and the sphere is close enough and if the particle is moving forward to the sphere
        // collision happens
        if (relative_position.norm() - radius < eEPSILON && relative_position.dot(particle.getVelocity()) < 0) {
            Eigen::Vector3f normal_velocity = unit.dot(particle.getVelocity()) * unit;
            Eigen::Vector3f tangent_velocity = particle.getVelocity() - normal_velocity;

            particle.setVelocity(coefResist * normal_velocity * (particle.getMass() - mass) /
                                     (particle.getMass() + mass) +
                                 tangent_velocity);

            // if a force pushes the particle into the sphere
            // exert a contact force (resist + friction) to resist it
            if (unit.dot(particle.getForce()) < 0) {
                Eigen::Vector3f resist = -unit.dot(particle.getForce()) * unit;
                Eigen::Vector3f friction =
                    coefFriction * unit.dot(particle.getForce()) * tangent_velocity / tangent_velocity.norm();
                particle.addForce(resist + friction);
            }
        }
    }
}

// BowlTerrain //

BowlTerrain::BowlTerrain() { modelMatrix = util::translate(position) * util::scale(radius, radius, radius); }

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    // TODO

    for (int i = 0; i < cube.getParticleNum(); i++) {
        Particle& particle = cube.getParticle(i);

        Eigen::Vector3f relative_position = particle.getPosition() - position;
        Eigen::Vector3f unit = relative_position / relative_position.norm();

        // if the particle and the bowl is close enough and if the particle is moving forward to the bowl
        // collision happens
        if (radius - relative_position.norm() < eEPSILON && relative_position.dot(particle.getVelocity()) > 0) {
            Eigen::Vector3f normal_velocity = unit.dot(particle.getVelocity()) * unit;
            Eigen::Vector3f tangent_velocity = particle.getVelocity() - normal_velocity;


            particle.setVelocity(coefResist * normal_velocity * (particle.getMass() - mass) /
                                     (particle.getMass() + mass) +
                                 tangent_velocity);

            // if a force pushes the particle into the bowl
            // exert a contact force (resist + friction) to resist it
            if (unit.dot(particle.getForce()) > 0) {
                Eigen::Vector3f resist = - unit.dot(particle.getForce()) * unit;
                Eigen::Vector3f friction =
                    -coefFriction * unit.dot(particle.getForce()) * tangent_velocity / tangent_velocity.norm();
                particle.addForce(resist + friction);
            }
        }
    }
}

// TiltedPlaneTerrain //

TiltedPlaneTerrain::TiltedPlaneTerrain() { modelMatrix = util::rotateDegree(0, 0, -45) * util::scale(60, 1, 60); }

TerrainType TiltedPlaneTerrain::getType() { return TerrainType::TiltedPlane; }

void TiltedPlaneTerrain::handleCollision(const float delta_T, Cube& cube) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.3f;
    // TODO

    Eigen::Vector3f unit_normal = normal / normal.norm();

    for (int i = 0; i < cube.getParticleNum(); i++) {
        Particle& particle = cube.getParticle(i);

        Eigen::Vector3f relative_position = particle.getPosition() - position;

        // if the particle and the wall is close enough and if the particle is moving forward to the wall
        // collision happens
        if (unit_normal.dot(relative_position) < eEPSILON && unit_normal.dot(particle.getVelocity()) < 0) {
            Eigen::Vector3f normal_velocity = unit_normal.dot(particle.getVelocity()) * unit_normal;
            Eigen::Vector3f tangent_velocity = particle.getVelocity() - normal_velocity;

            particle.setVelocity(-coefResist * normal_velocity + tangent_velocity);

            // if a force pushes the particle into the wall
            // exert a contact force (resist + friction) to resist it
            if (unit_normal.dot(particle.getForce()) < 0) {
                Eigen::Vector3f resist = -unit_normal.dot(particle.getForce()) * unit_normal;
                Eigen::Vector3f friction =
                    coefFriction * unit_normal.dot(particle.getForce()) * tangent_velocity / tangent_velocity.norm();
                particle.addForce(resist + friction);
            }
        }
    }
}
}  // namespace simulation
