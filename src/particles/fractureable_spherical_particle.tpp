//
// Created by erolsson on 2019-01-19.
//

#include "fractureable_spherical_particle.h"

template<typename ForceModel>
DEM::FractureableSphericalParticle<ForceModel>::FractureableSphericalParticle(double radius, const Vec3 &position,
                                                                              const Vec3 &velocity,
                                                                              MaterialBase *material, unsigned id)
        :SphericalParticleBase<ForceModel>(radius, position, velocity, material, id) {

}