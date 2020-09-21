//
// Created by erolsson on 2018-09-07.
//

#include "spherical_particle.h"

#include "../utilities/file_reading_functions.h"
#include "../utilities/vec3.h"

template<typename ForceModel>
DEM::SphericalParticle<ForceModel>::SphericalParticle(double radius, const DEM::Vec3& position,
                                                      const DEM::Vec3& velocity, const MaterialBase* material,
                                                      std::size_t object_id, std::size_t collision_id):
   SphericalParticleBase<ForceModel>(radius, position, velocity, material, object_id, collision_id)

{
    // Empty constructor
}

template<typename ForceModel>
DEM::SphericalParticle<ForceModel>::SphericalParticle(const DEM::ParameterMap& parameters, DEM::MaterialBase *material):
    SphericalParticleBase<ForceModel>(parameters, material)
{
    // Empty constructor
}





