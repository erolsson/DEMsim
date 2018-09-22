//
// Created by erolsson on 2018-08-24.
//

#ifndef DEMSIM_CYLINDER_H
#define DEMSIM_CYLINDER_H

#include <array>
#include <sstream>

#include "surface_base.h"
#include "vec3.h"

namespace DEM {
    // ToDo: update this class to allow for cylinders with arbitrary orientation, currently cylinders aligned in the
    //       z-axis are only working properly in the collision handling
    template<typename ForceModel, typename ParticleType>
    class Cylinder : public Surface<ForceModel, ParticleType> {
    public:
        Cylinder(std::size_t id, double radius, const Vec3& axis, const Vec3& base_point, double length,
                bool inward=true, bool infinite=false);
        ~Cylinder() override = default;

        Vec3 get_normal(const Vec3& position) const override;

        double distance_to_point(const Vec3& point) const override;
        Vec3 vector_to_point(const Vec3& point) const override;
        Vec3 displacement_this_inc(const Vec3& position) const override;

        void move(const Vec3& distance, const Vec3& velocity) override;
        void rotate(const Vec3& position, const Vec3& rotation_vector) override;

        std::string output_data() const override;

        void expand(double radius_increase);
        double get_radius() const { return radius_; }

    private:
        using Surface<ForceModel, ParticleType>::id_;
        using Surface<ForceModel, ParticleType>::displacement_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_point_;
        using Surface<ForceModel, ParticleType>::velocity_;
        using Surface<ForceModel, ParticleType>::bbox_values_;

        double radius_;
        Vec3 axis_;
        Vec3 point_;
        double length_;
        bool inward_;
        bool infinite_;
        // To allow fast functions for a common case, currently only z-aligned is supported in the collision detector!!!
        bool z_aligned_;

        void update_bounding_box();
    };
}

#include "cylinder.tpp"

#endif //DEMSIM_CYLINDER_H
