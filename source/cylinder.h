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
        Cylinder(unsigned id, double radius, Vec3 axis, Vec3 center_point, double length, bool inward=false);
        ~Cylinder() override = default;
        Vec3 get_normal(const Vec3& position) const override;
        double distance_to_point(const Vec3& point) const override;
        Vec3 vector_to_point(const Vec3& point) const override;
        Vec3 displacement_this_inc(const Vec3& position) const override;
        void move(const Vec3& distance, const Vec3& velocity) override;
        void rotate(const Vec3& position, const Vec3& rotation_vector) override;
        std::string output_data() const override;
        std::pair<Vec3, Vec3> bounding_box_values() const override;

        void expand(double radius_increase);
    private:
        using Surface<ForceModel, ParticleType>::id_;
        using Surface<ForceModel, ParticleType>::displacement_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_point_;
        using Surface<ForceModel, ParticleType>::velocity_;

        double radius_;
        Vec3 axis_;
        Vec3 point_;
        double length_;
        bool inward_;
    };


    template<typename ForceModel, typename ParticleType>
    Cylinder<ForceModel, ParticleType>::Cylinder(unsigned id, double radius, Vec3 axis,
                                                 Vec3 center_point, double length, bool inward) :
        Surface<ForceModel, ParticleType>::Surface(id),
        radius_(radius),
        axis_(axis),
        point_(center_point),
        length_(length),
        inward_(inward)
    {
        // Empty constructor
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Cylinder<ForceModel, ParticleType>::get_normal(const Vec3& position) const
    {
        Vec3 n = (position - point_) - dot_product(axis_, position)*axis_;
        if (inward_)
            return -n.normalize();
        return n.normalize();
    }

    template<typename ForceModel, typename ParticleType>
    double Cylinder<ForceModel, ParticleType>::distance_to_point(const Vec3& point) const
    {
        Vec3 n = get_normal(point);
        return radius_ + dot_product((point - point_), n);
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Cylinder<ForceModel, ParticleType>::vector_to_point(const Vec3& point) const
    {
        Vec3 n = get_normal(point);
        return (radius_ + dot_product((point - point_), n))*n;
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Cylinder<ForceModel, ParticleType>::displacement_this_inc(const Vec3& position) const
    {
        return displacement_this_inc_;
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

    template<typename ForceModel, typename ParticleType>
    void Cylinder<ForceModel, ParticleType>::expand(double radius_increase)
    {
        radius_ += radius_increase;
    }


}

#endif //DEMSIM_CYLINDER_H
