//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_ENGINE_H
#define DEMSIM_ENGINE_H

#include <cstddef>
#include <iostream>
#include <map>
#include <vector>

#include "contact.h"
#include "contact_matrix.h"
#include "cylinder.h"
#include "material_base.h"
#include "point_surface.h"
#include "surface.h"
#include "vec3.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Engine {
    public:
        using ParticlePointer = ParticleType*;
        using PointSurfacePointer = PointSurface<ForceModel, ParticleType>*;
        using CylinderPointer = Cylinder<ForceModel, ParticleType>*;
        Engine();

        //Object creation functions
        template<typename MaterialType>
        MaterialType* create_material(double density);
        ParticlePointer create_particle(double radius, const Vec3& position, const Vec3& velocity,
                                      MaterialBase* material);
        PointSurfacePointer create_point_surface(const std::vector<Vec3>& points, bool infinite);
        CylinderPointer create_cylinder(double radius, const Vec3& axis, const Vec3& base_point, double length,
                                        bool inward=true, bool infinite=false);

    private:
        using ContactType = Contact<ForceModel, ParticleType>;
        using SurfaceType = Surface<ForceModel, ParticleType>;

        std::size_t number_of_objects_;
        std::vector<MaterialBase*> materials_;
        std::vector<ParticleType*> particles_;
        std::vector<SurfaceType*> surfaces_;
        ContactMatrix<ContactType> contacts_;

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
    MaterialType* Engine<ForceModel, ParticleType>::create_material(double density)
    {
        auto m = new MaterialType(materials_.size(), density);
        materials_.push_back(m);
        return m;
    }

    template<typename ForceModel, typename ParticleType>
    typename Engine<ForceModel, ParticleType>::ParticlePointer
    Engine<ForceModel, ParticleType>::create_particle(double radius, const Vec3& position,
                                                      const Vec3& velocity, MaterialBase* material)
    {
        auto p = new ParticleType(radius, position, velocity, material, number_of_objects_);
        particles_.push_back(p);
        ++number_of_objects_;
        return p;
    }

    template<typename ForceModel, typename ParticleType>
    typename Engine<ForceModel, ParticleType>::PointSurfacePointer
    Engine<ForceModel, ParticleType>::create_point_surface(const std::vector<Vec3>& points, bool infinite)
    {
        auto ps = new PointSurface<ForceModel, ParticleType>(number_of_objects_, points, infinite);
        surfaces_.push_back(ps);
        ++number_of_objects_;
        return ps;
    }

    template<typename ForceModel, typename ParticleType>
    typename Engine<ForceModel, ParticleType>::CylinderPointer
    Engine<ForceModel, ParticleType>::create_cylinder(double radius, const Vec3& axis, const Vec3& base_point,
                                                      double length, bool inward, bool infinite)
    {
        auto c = new Cylinder<ForceModel, ParticleType>(number_of_objects_, radius, axis, base_point, length,
                                                          inward, infinite);
        surfaces_.push_back(c);
        ++number_of_objects_;
        return c;
    }


}

#endif //DEMSIM_ENGINE_H
