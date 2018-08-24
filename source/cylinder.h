//
// Created by erolsson on 2018-08-24.
//

#ifndef DEMSIM_CYLINDER_H
#define DEMSIM_CYLINDER_H

#include "surface.h"
#include "vec3.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Cylinder : public Surface<ForceModel, ParticleType> {
    public:
        Cylinder(unsigned id, double radius, Vec3 axis, Vec3 center_point, bool invard=false);
        ~Cylinder() override = default;
        Vec3 get_normal(const Vec3& position) const override;
        double distance_to_point(const Vec3& point) const override;
        Vec3 vector_to_point(const Vec3& point) const override;
        Vec3 displacement_this_inc(const Vec3& position) const override;
        void move(const Vec3& distance, const Vec3& velocity) override;
        void rotate(const Vec3& position, const Vec3& rotation_vector) override;
        std::string output_data() const override;
        std::pair<Vec3, Vec3> bounding_box_values() const override;

    private:
        using Surface<ForceModel, ParticleType>::id_;
        using Surface<ForceModel, ParticleType>::displacement_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_point_;
        using Surface<ForceModel, ParticleType>::velocity_;

        double radius_;
        Vec3 axis_;
        Vec3 point_;
        bool invard_;
    };


    template<typename ForceModel, typename ParticleType>
    Cylinder<ForceModel, ParticleType>::Cylinder(unsigned id, double radius, Vec3 axis,
                                                 Vec3 center_point, bool invard) :
        Surface<ForceModel, ParticleType>::Surface(id),
        radius_(radius),
        axis_(axis),
        point_(center_point),
        invard_(invard)
    {
        // Empty constructor
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Cylinder<ForceModel, ParticleType>::get_normal(const Vec3& position) const
    {
        return Vec3();
    }

    template<typename ForceModel, typename ParticleType>
    double Cylinder<ForceModel, ParticleType>::distance_to_point(const Vec3& point) const
    {
        return 0;
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Cylinder<ForceModel, ParticleType>::vector_to_point(const Vec3& point) const
    {
        return Vec3();
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Cylinder<ForceModel, ParticleType>::displacement_this_inc(const Vec3& position) const
    {
        return Vec3();
    }

    template<typename ForceModel, typename ParticleType>
    void Cylinder<ForceModel, ParticleType>::move(const Vec3& distance, const Vec3& velocity)
    {

    }

    template<typename ForceModel, typename ParticleType>
    void Cylinder<ForceModel, ParticleType>::rotate(const Vec3& position, const Vec3& rotation_vector)
    {

    }

    template<typename ForceModel, typename ParticleType>
    std::string Cylinder<ForceModel, ParticleType>::output_data() const
    {
        return std::__cxx11::string();
    }

    template<typename ForceModel, typename ParticleType>
    std::pair<Vec3, Vec3> Cylinder<ForceModel, ParticleType>::bounding_box_values() const
    {
        return std::pair<Vec3, Vec3>();
    }



}

#endif //DEMSIM_CYLINDER_H
