//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_ENGINE_H
#define DEMSIM_ENGINE_H

#include <cstddef>
#include <iostream>
#include <map>
#include <vector>

#include "vec3.h"
#include "contact_matrix.h"
#include "material_base.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Engine {
    public:
        Engine();

        template<typename MaterialType>
        MaterialType* create_material();

        ParticleType* create_particle(double radius, const Vec3& position, const Vec3& velocity,
                                      MaterialBase* material);


    private:
        std::size_t number_of_objects_;
        std::vector<MaterialBase*> materials_;
        std::vector<ParticleType*> particles_;

    };

    template<typename ForceModel, typename ParticleType>
    Engine<ForceModel, ParticleType>::Engine() :
        number_of_objects_(0),
        particles_()
    {
        // Empty constructor
    }

    template<typename ForceModel, typename ParticleType>
    template<typename MaterialType>
    MaterialType* Engine<ForceModel, ParticleType>::create_material()
    {
        auto* mat = new MaterialType();
        materials_.push_back(mat);
        return mat;
    }

    template<typename ForceModel, typename ParticleType>
    ParticleType* Engine<ForceModel, ParticleType>::create_particle(double radius, const Vec3& position,
                                                                    const Vec3& velocity, MaterialBase* material)
    {
        auto p = new ParticleType(radius, position, velocity, material, number_of_objects_);
        particles_.push_back(p);
        ++number_of_objects_;
        return p;
    }




}

#endif //DEMSIM_ENGINE_H
