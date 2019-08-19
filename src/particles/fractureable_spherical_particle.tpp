//
// Created by erolsson on 2019-01-19.
//

#include "fractureable_spherical_particle.h"

#include <random>
#include <sstream>

#include "../materials/stone_material.h"

template<typename ForceModel>
DEM::FractureableSphericalParticle<ForceModel>::FractureableSphericalParticle(double radius, const Vec3 &position,
                                                                              const Vec3 &velocity,
                                                                              MaterialBase *material, unsigned id)
        :SphericalParticleBase<ForceModel>(radius, position, velocity, material, id) {

}


template<typename ForceModel>
void DEM::FractureableSphericalParticle<ForceModel>::fracture_particle(const DEM::Vec3& position, double force,
        std::size_t id_of_impacter, const Vec3& normal) {
    auto crack_iter = has_crack_at_position(position);
    if (crack_iter != cracks_.end()) {
        if (force > crack_iter->get_force()) {
            crack_iter->set_force(force);
        }
    }
    else {
        std::cout << "Particle " << get_id() << " Fractures at position " << position
                  << " with force " << force << " from object " << id_of_impacter << "\n";
        # pragma omp critical
        {
            cracks_.emplace_back(position, force, id_of_impacter, normal);
        }
    }

}

template<typename ForceModel>
std::vector<DEM::ParticleCrack>::iterator
DEM::FractureableSphericalParticle<ForceModel>::has_crack_at_position(const DEM::Vec3& position) {
    auto mat = dynamic_cast<const DEM::StoneMaterial*>(material_);
    double min_crack_distance = mat->min_crack_distance;
    for(auto iter = cracks_.begin();  iter != cracks_.end(); ++iter) {
        const auto& pos = iter->get_position();
        if ( (pos - position).length() < min_crack_distance) {
            return iter;
        }
    }
    return cracks_.end();
}



