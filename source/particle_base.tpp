//
// Created by erolsson on 2018-09-02.
//

#include "particle_base.h"

template<typename ForceModel>
DEM::ParticleBase<ForceModel>::ParticleBase(double mass, const Vec3& pos, const Vec3& velocity, MaterialBase* m,
                                       unsigned id) :
        id_(id),
        mass_(mass),
        material_(m),
        position_(pos),
        velocity_(velocity)
{
    // Empty constructor
}