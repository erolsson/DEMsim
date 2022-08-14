//
// Created by erolsson on 2018-09-02.
//

#include "cylinder.h"

#include <sstream>
#include <string>

#include "../utilities/printing_functions.h"

template<typename ForceModel, typename ParticleType>
DEM::Cylinder<ForceModel, ParticleType>::Cylinder(std::size_t id, double radius, const Vec3& axis,
                                             const Vec3& base_point, double length, const std::string& name,
                                             bool inward, bool infinite, bool closed_ends, std::size_t collision_id) :
        Surface<ForceModel, ParticleType>::Surface(id, collision_id,  name),
        radius_(radius),
        axis_(axis.normal()),
        point_(base_point),
        length_(length),
        inward_(inward),
        infinite_(infinite),
        closed_ends_(closed_ends),
        z_aligned_(axis_ == Vec3(0, 0, 1))
{
    update_bounding_box();
}

template<typename ForceModel, typename ParticleType>
DEM::Cylinder<ForceModel, ParticleType>::Cylinder(const DEM::ParameterMap& parameters) :
        Surface<ForceModel, ParticleType>::Surface(parameters),
        radius_(parameters.get_parameter<double>("radius")),
        axis_(parameters.get_vec3("axis")),
        point_(parameters.get_vec3("point")),
        length_(parameters.get_parameter<double>("length")),
        inward_(parameters.get_parameter<bool>("inward")),
        infinite_(parameters.get_parameter<bool>("infinite")),
        closed_ends_(parameters.get_parameter<bool>("closed_ends")),
        z_aligned_(parameters.get_parameter<bool>("z_aligned"))
{
    update_bounding_box();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Cylinder<ForceModel, ParticleType>::get_normal(const Vec3& position) const
{
    Vec3 n = (position - point_) - dot_product(axis_, position - point_)*axis_;
    if (n.is_zero()) {
        // Special case, we are on the central axis, any unit vector can be used as normal
        return Vec3(1, 0, 0);
    }
    if (inward_)
        n *= -1;

    // Only cylinders aligned with the z-axis are currently implemented
    if (closed_ends_ && !inward_) {
        double position_on_axis = dot_product(position-point_, axis_);
        if (position_on_axis < 0 || position_on_axis > length_) {
            double r2 = pow(position.x() - point_.x(), 2) + pow(position.y() - point_.y(), 2);
            if (r2 <= radius_*radius_) {
                int sgn = position_on_axis > length_ ? -1 : 1;
                return Vec3(0., 0., 1.)*sgn;
            }
            else {
                double d = position_on_axis <= point_.z() ? 0 : length_;
                Vec3 point_on_surface = point_ + axis_*d +
                        radius_*Vec3(position.x() - point_.x(), position.y() - point_.y(), 0.)/sqrt(r2);
                n = position - point_on_surface;
            }
        }

    }
    n.normalize();
    return n;
}

template<typename ForceModel, typename ParticleType>
double DEM::Cylinder<ForceModel, ParticleType>::distance_to_point(const Vec3& point) const
{
    if (!infinite_) {
        return vector_to_point(point).length();
    }
    Vec3 n = get_normal(point);
    int sgn = inward_? 1 : -1;
    return (sgn*radius_ + dot_product((point-point_), n));
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Cylinder<ForceModel, ParticleType>::vector_to_point(const Vec3& point) const
{
    if (!infinite_) {
        double position_on_axis = dot_product(point - point_, axis_);
        if (position_on_axis < 0 || position_on_axis > length_) {
            double d =  position_on_axis < point_.z() ? 0 : length_;
            Vec3 point_on_surface;
            if (!closed_ends_){
                point_on_surface = point_ + axis_*d - get_normal(point)*radius_;
            }
            else {
                double r2 = pow(point.x() - point_.x(), 2) + pow(point.y() - point_.y(), 2);
                if (r2 <= radius_*radius_) {
                    return {0., 0., position_on_axis};
                }
                point_on_surface = point_ + axis_*d +
                        radius_*Vec3(point.x() - point_.x(), point.y() - point_.y(), 0.)/sqrt(r2);
            }
            return point - point_on_surface;
        }
    }
    Vec3 n = get_normal(point);
    int sgn = inward_? 1 : -1;
    return (sgn*radius_ + dot_product((point-point_), n))*n;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Cylinder<ForceModel, ParticleType>::get_displacement_this_increment(const Vec3& position) const
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
    stream << "ID=" << get_id() << ", TYPE=Cylinder, " << radius_  << ", " << axis_.x() << ", " << axis_.y() << ", "
           << axis_.z() << ", " << point_.x() << ", " << point_.y() << ", " << point_.z() << ", " << length_;
    return stream.str();
}

template<typename ForceModel, typename ParticleType>
std::string DEM::Cylinder<ForceModel, ParticleType>::restart_data() const {
    using DEM::named_print;
    std::ostringstream ss;
    ss << DEM::Surface<ForceModel, ParticleType>::restart_data() << ", "
       << named_print(radius_, "radius") << ", " << named_print(axis_, "axis") << ", "
       << named_print(point_, "point") << ", " << named_print(length_, "length") << ", "
       << named_print(inward_, "inward") << ", " << named_print(infinite_, "infinite") << ","
       << named_print(closed_ends_, "closed_ends") << ", " << named_print(z_aligned_, "z_aligned");
    return ss.str();
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
