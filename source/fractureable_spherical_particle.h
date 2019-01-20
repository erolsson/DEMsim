//
// Created by erolsson on 2019-01-19.
//

#ifndef DEMSIM_FRACTUREABLE_SPHERICAL_PARTICLE_H
#define DEMSIM_FRACTUREABLE_SPHERICAL_PARTICLE_H

#include "spherical_particle.h"

namespace DEM {

    template<typename ForceModel>
    class FractureableSphericalParticle : public SphericalParticle<ForceModel> {

    };
}

#include "fractureable_spherical_particle.tpp"

#endif //DEMSIM_FRACTUREABLE_SPHERICAL_PARTICLE_H
