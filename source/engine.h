//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_ENGINE_H
#define DEMSIM_ENGINE_H

#include <cstddef>
#include <iostream>
#include <map>
#include <vector>

#include "collision_detector.h"
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

        std::size_t number_of_objects_{ 0 };
        std::vector<MaterialBase*> materials_{};
        std::vector<ParticleType*> particles_{};
        std::vector<SurfaceType*> surfaces_{};
        ContactMatrix<ContactType> contacts_{};

        CollisionDetector<ForceModel, ParticleType> collision_detector_;
    };
}

#include "engine.tpp"

#endif //DEMSIM_ENGINE_H
