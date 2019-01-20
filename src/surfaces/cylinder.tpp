//
// Created by erolsson on 2018-09-02.
//

#include "cylinder.h"

#include <sstream>

template<typename ForceModel, typename ParticleType>
DEM::Cylinder<ForceModel, ParticleType>::Cylinder(std::size_t id, double radius, const Vec3& axis,
                                             const Vec3& base_point, double length, bool inward, bool infinite) :
        Surface<ForceModel, ParticleType>::Surface(id),
        radius_(radius),
        axis_(axis.normal()),
        point_(base_point),
        length_(length),
        inward_(inward),
        infinite_(infinite),
        z_aligned_(axis_ == Vec3(0, 0, 1))
{
    update_bounding_box();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Cylinder<ForceModel, ParticleType>::get_normal(const Vec3& position) const
{
    Vec3 n = (position-point_) - dot_product(axis_, position-point_)*axis_;
    if (inward_)
        n*= -1;
    if (n.is_zero()) {
        return Vec3(1, 0, 0);  //Special case, we are on the central axis, any unit vector can be used as normal
    }
    return n.normalize();
}

template<typename ForceModel, typename ParticleType>
double DEM::Cylinder<ForceModel, ParticleType>::distance_to_point(const Vec3& point) const
{
    Vec3 n = get_normal(point);
    if (!infinite_) {
        double position_on_axis = dot_product(point-point_, axis_);
        if (position_on_axis < 0 || position_on_axis > length_) {
            Vec3 point_on_surface = point_+axis_*position_on_axis;
            return (point-point_on_surface).length();
        }
    }
    return (radius_ + dot_product((point-point_), n));
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Cylinder<ForceModel, ParticleType>::vector_to_point(const Vec3& point) const
{
    Vec3 n = get_normal(point);
    if (!infinite_) {
        double position_on_axis = dot_product(point - point_, axis_);
        if (position_on_axis < 0 || position_on_axis > length_) {
            double d = position_on_axis < 0 ? 0. : length_;     // Should we be on the upper or lower plate
            Vec3 point_on_surface = point_ + axis_*d - get_normal(point)*radius_;
            return point - point_on_surface;
        }
    }

    return (radius_ + dot_product((point-point_), n))*n;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Cylinder<ForceModel, ParticleType>::displacement_this_inc(const Vec3& position) const
{
    if (rotation_this_inc_.is_zero())
        return displacement_this_inc_;
    return cross_product(rotation_this_inc_, position - rotation_point_) + displacement_this_inc_;
}

template<typename ForceModel, typename ParticleType>
void DEM::Cylinder<ForceModel, ParticleType>::move(const Vec3& distance, const Vec3& velocity)
{
    point_ += distance;
    velocity_ = velocity;
    displacement_this_inc_ = distance;
    update_bounding_box();
}

template<typename ForceModel, typename ParticleType>
void DEM::Cylinder<ForceModel, ParticleType>::rotate(const Vec3& position, const Vec3& rotation_vector)
{
    point_ += cross_product(rotation_vector, point_ - position);
    Vec3 p1 = cross_product(rotation_vector, point_ - position + axis_);
    axis_ = (p1 + axis_).normal();

    rotation_this_inc_ = rotation_vector;
    rotation_point_ = position;
    z_aligned_ = rotation_vector == Vec3(0, 0, 0);
    update_bounding_box();
}

template<typename ForceModel, typename ParticleType>
std::string DEM::Cylinder<ForceModel, ParticleType>::get_output_string() const
{
    std::ostringstream stream;
    stream << "ID=" << id_ << ", TYPE=Cylinder, " << radius_  << ", " << axis_.x() << ", " << axis_.y() << ", "
           << axis_.z() << ", " << point_.x() << ", " << point_.y() << ", " << point_.z() << ", " << length_;
    return stream.str();
}

template<typename ForceModel, typename ParticleType>
void DEM::Cylinder<ForceModel, ParticleType>::update_bounding_box()
{
    if (inward_) {
        bbox_values_[0] = point_.x() - radius_/sqrt(2);
        bbox_values_[1] = point_.x() + radius_/sqrt(2);
        bbox_values_[2] = point_.y() - radius_/sqrt(2);
        bbox_values_[3] = point_.y() + radius_/sqrt(2);
    }

    else {
        bbox_values_[0] = point_.x() - radius_;
        bbox_values_[1] = point_.x() + radius_;
        bbox_values_[2] = point_.y() - radius_;
        bbox_values_[3] = point_.y() + radius_;
    }

    if (!infinite_) {
        bbox_values_[4] = point_.z();
        bbox_values_[5] = point_.z() + length_;
    }

    else {
        bbox_values_[4] = -1e99;
        bbox_values_[5] = 1e99;
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Cylinder<ForceModel, ParticleType>::expand(double radius_increase)
{
    radius_ += radius_increase;
    update_bounding_box();
}
