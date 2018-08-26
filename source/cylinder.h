//
// Created by erolsson on 2018-08-24.
//

#ifndef DEMSIM_CYLINDER_H
#define DEMSIM_CYLINDER_H

#include <sstream>

#include "surface.h"
#include "vec3.h"

namespace DEM {
    // ToDo: update this class to allow for cylinders with arbitrary orientation, currently cylinders aligned in the
    //       z-axis are only working properly in the collision handling
    template<typename ForceModel, typename ParticleType>
    class Cylinder : public Surface<ForceModel, ParticleType> {
    public:
        Cylinder(unsigned id, double radius, Vec3 axis, Vec3 base_point, double length,
                bool inward=true, bool infinite=false);
        ~Cylinder() override = default;

        Vec3 get_normal(const Vec3& position) const override;

        double distance_to_point(const Vec3& point) const override;
        Vec3 vector_to_point(const Vec3& point) const override;
        Vec3 displacement_this_inc(const Vec3& position) const override;

        void move(const Vec3& distance, const Vec3& velocity) override;
        void rotate(const Vec3& position, const Vec3& rotation_vector) override;

        std::string output_data() const override;
        std::pair<Vec3, Vec3> bounding_box_values() const override;

        void expand(double radius_increase) { radius_ += radius_increase; }
        double get_radius() const { return radius_; }

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
        bool infinite_;
        // To allow fast functions for a common case, currently only z-aligned is supported in the collision detector!!!
        bool z_aligned_;
    };


    template<typename ForceModel, typename ParticleType>
    Cylinder<ForceModel, ParticleType>::Cylinder(unsigned id, double radius, Vec3 axis,
                                                 Vec3 base_point, double length, bool inward, bool infinite) :
        Surface<ForceModel, ParticleType>::Surface(id),
        radius_(radius),
        axis_(axis.normalize()),
        point_(base_point),
        length_(length),
        inward_(inward),
        infinite_(infinite),
        z_aligned_(axis_ == Vec3(0, 0, 1))
    {
        // Empty constructor
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Cylinder<ForceModel, ParticleType>::get_normal(const Vec3& position) const
    {
        Vec3 n = (position-point_) - dot_product(axis_, position)*axis_;
        if (inward_)
            n*= -1;
        if (n.is_zero())
            return Vec3(1, 0, 0);  //Special case, we are on the central axis, any unit vector can be used as normal
        return n.normalize();
    }

    template<typename ForceModel, typename ParticleType>
    double Cylinder<ForceModel, ParticleType>::distance_to_point(const Vec3& point) const
    {
        Vec3 n = get_normal(point);
        if (!infinite_) {
            double position_on_axis = dot_product(point-point_, axis_);
            if (position_on_axis<0 || position_on_axis>length_) {
                Vec3 point_on_surface = point_+axis_*position_on_axis;
                return (point-point_on_surface).length();
            }
        }
        return (radius_ + dot_product((point-point_), n));
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Cylinder<ForceModel, ParticleType>::vector_to_point(const Vec3& point) const
    {
        Vec3 n = get_normal(point);
        if (!infinite_) {
            double position_on_axis = dot_product(point - point_, axis_);
            if (position_on_axis < 0 || position_on_axis > length_) {
                Vec3 point_on_surface = point_+axis_*position_on_axis;
                return point - point_on_surface;
            }
        }
        return (radius_ + dot_product((point-point_), n))*n;
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Cylinder<ForceModel, ParticleType>::displacement_this_inc(const Vec3& position) const
    {
        if (rotation_this_inc_.is_zero())
            return displacement_this_inc_;
        return cross_product(rotation_this_inc_, position - rotation_point_) + displacement_this_inc_;
    }

    template<typename ForceModel, typename ParticleType>
    void Cylinder<ForceModel, ParticleType>::move(const Vec3& distance, const Vec3& velocity)
    {
        point_ += distance;
        velocity_ = velocity;
        displacement_this_inc_ = distance;
    }

    template<typename ForceModel, typename ParticleType>
    void Cylinder<ForceModel, ParticleType>::rotate(const Vec3& position, const Vec3& rotation_vector)
    {
        point_ += cross_product(rotation_vector, point_ - position);
        Vec3 p1 = cross_product(rotation_vector, point_ - position + axis_);
        axis_ = (p1 + axis_).normal();

        rotation_this_inc_ = rotation_vector;
        rotation_point_ = position;
        z_aligned_ = rotation_vector == Vec3(0, 0, 0);
    }

    template<typename ForceModel, typename ParticleType>
    std::string Cylinder<ForceModel, ParticleType>::output_data() const
    {
        std::ostringstream stream;
        stream << id_ << ", " << radius_;
        return stream.str();
    }

    template<typename ForceModel, typename ParticleType>
    std::pair<Vec3, Vec3> Cylinder<ForceModel, ParticleType>::bounding_box_values() const
    {
        //Cylinders not aligned with the z-axis not implemented but support exists
        double x_max = point_.x;
        double x_min = point_.x;
        double y_min = point_.y;
        double y_max = point_.y;
        double z_min = -1e99;
        double z_max = 1e99;
        if (inward_){
            x_min -= radius_/sqrt(2);
            x_max += radius_/sqrt(2);
            y_min -= radius_/sqrt(2);
            y_max += radius_/sqrt(2);
        } else {
            x_min -= radius_;
            x_max += radius_;
            y_min -= radius_;
            y_max += radius_;
        }
        if ( !infinite_ ) {
            z_min = point_.z;
            z_max = point_.z + length_;
        }
        return std::make_pair(Vec3(x_min, y_min, z_min), Vec3(x_max, y_max, z_max));
    }

}

#endif //DEMSIM_CYLINDER_H
