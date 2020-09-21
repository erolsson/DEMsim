//
// Created by erolsson on 2018-09-07.
//

#include "point_surface.h"

#include <sstream>
#include "string"

#include "../utilities/file_reading_functions.h"
#include "../utilities/printing_functions.h"
#include "../utilities/vec3.h"


template<typename ForceModel, typename ParticleType>
DEM::PointSurface<ForceModel, ParticleType>::PointSurface(std::size_t id, std::vector<DEM::Vec3> points, bool infinite,
                                                          const std::string& name, bool adhesive,
                                                          std::size_t collision_id) :
        Surface<ForceModel, ParticleType>::Surface(id, collision_id, name, adhesive),
        points_(std::move(points)),
        normal_(calculate_normal()),
        infinite_(infinite)
{
    update_bounding_box();
}


template<typename ForceModel, typename ParticleType>
DEM::PointSurface<ForceModel, ParticleType>::PointSurface(const DEM::ParameterMap& parameters) :
    Surface<ForceModel, ParticleType>::Surface(parameters),
    points_(),
    infinite_(parameters.get_parameter<bool>("infinite"))
{
    auto no_points = parameters.get_parameter<std::size_t>("no_points");
    for (unsigned i = 0; i != no_points; ++i) {
        std::ostringstream point;
        point << "p" << i;
        points_.push_back(parameters.get_vec3(point.str()));
    }
    normal_ = calculate_normal();
    update_bounding_box();
}


template<typename ForceModel, typename ParticleType>
double DEM::PointSurface<ForceModel, ParticleType>::distance_to_point(const DEM::Vec3& point) const
{
    return vector_to_point(point).length();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::PointSurface<ForceModel, ParticleType>::vector_to_point(const Vec3& point) const
{
    DEM::Vec3 v = dot_product(point-points_[0], normal_)*normal_;
    if(infinite_)
        return v;

    DEM::Vec3 plane_vector = point - v;  // Vector in the same plane as the surface going from origo to the point

    double min_distance = 1E99;
    DEM::Vec3 min_vector = Vec3(0, 0, 0);
    bool inside = true;

    for (unsigned i = 0; i != points_.size(); ++i){
        DEM::Vec3 vec{};
        DEM::Vec3 ps = plane_vector - points_[i];    // Vector from a corner to the point in the plane of the surfacve
        DEM::Vec3 edge{};
        if (i != points_.size() - 1)
            edge = points_[i+1] - points_[i];
        else
            edge = points_[points_.size() - 1] - points_[0];

        double l = dot_product(ps, edge.normal());

        if ( l < 0 )
            vec = ps;
        else if ( l > edge.length() )
            vec = ps-edge;
        else
            vec = ps-l*edge.normal();
        l = vec.length();
        inside = inside && dot_product(cross_product(vec, edge), normal_) < 0;
        if ( l < min_distance ){
            min_distance = l;
            min_vector = vec+v;
        }
    }

    if ( inside )
        return v;
    return min_vector;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::PointSurface<ForceModel, ParticleType>::get_displacement_this_increment(const Vec3& position) const
{
    if (rotation_this_inc_.is_zero())
        return displacement_this_inc_;
    return displacement_this_inc_ + cross_product(rotation_this_inc_, position - rotation_point_);
}

template<typename ForceModel, typename ParticleType>
void DEM::PointSurface<ForceModel, ParticleType>::move(const Vec3& distance, const Vec3& velocity)
{
    for (auto& p: points_){
        p += distance;
    }
    displacement_this_inc_ = distance;
    velocity_ = velocity;
    normal_ = calculate_normal();
    update_bounding_box();
}

template<typename ForceModel, typename ParticleType>
void DEM::PointSurface<ForceModel, ParticleType>::rotate(const Vec3& point, const Vec3& rotation_vector)
{
    for(auto& p: points_) {
        p += cross_product(rotation_vector, p-point);
    }
    rotation_this_inc_ = rotation_vector;
    rotation_point_ = point;
    normal_ = calculate_normal();
    update_bounding_box();
}

template<typename ForceModel, typename ParticleType>
std::string DEM::PointSurface<ForceModel, ParticleType>::get_output_string() const
{
    std::ostringstream stream;
    stream << "ID=" << get_id() << ", TYPE=PointSurface, " << points_.size();
    for(auto& p: points_) {
        stream << ", " << p.x() << ", " << p.y() << ", " << p.z();
    }
    return stream.str();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::PointSurface<ForceModel, ParticleType>::calculate_normal() const
{
    DEM::Vec3 normal = cross_product(points_[1]-points_[0], points_[2]-points_[1]);
    return normal.normalize();
}

template<typename ForceModel, typename ParticleType>
void DEM::PointSurface<ForceModel, ParticleType>::update_bounding_box()
{
    auto x_cmp = [](const DEM::Vec3& v1, const DEM::Vec3& v2) -> bool {
        return v1.x() < v2.x();
    };

    auto y_cmp = [](const DEM::Vec3& v1, const DEM::Vec3& v2) -> bool {
        return v1.y() < v2.y();
    };

    auto z_cmp = [](const DEM::Vec3& v1, const DEM::Vec3& v2) -> bool {
        return v1.z() < v2.z();
    };

    if (infinite_) {
        bbox_values_ = {-1e99, 1e99, -1e99, 1e99, -1e99, 1e99};
        if (normal_[0] != 0) {
            bbox_values_[0] = std::min_element(points_.begin(), points_.end(), x_cmp)->x();
            bbox_values_[1] = std::max_element(points_.begin(), points_.end(), x_cmp)->x();
        }
        if (normal_[1] != 0) {
            bbox_values_[2] = std::min_element(points_.begin(), points_.end(), y_cmp)->y();
            bbox_values_[3] = std::max_element(points_.begin(), points_.end(), y_cmp)->y();
        }
        if (normal_[2] != 0) {
            bbox_values_[4] = std::min_element(points_.begin(), points_.end(), z_cmp)->z();
            bbox_values_[5] = std::max_element(points_.begin(), points_.end(), z_cmp)->z();
        }
    }
    else {
        bbox_values_[0] = std::min_element(points_.begin(), points_.end(), x_cmp)->x();
        bbox_values_[1] = std::max_element(points_.begin(), points_.end(), x_cmp)->x();

        bbox_values_[2] = std::min_element(points_.begin(), points_.end(), y_cmp)->y();
        bbox_values_[3] = std::max_element(points_.begin(), points_.end(), y_cmp)->y();

        bbox_values_[4] = std::min_element(points_.begin(), points_.end(), z_cmp)->z();
        bbox_values_[5] = std::max_element(points_.begin(), points_.end(), z_cmp)->z();
    }
}

template<typename ForceModel, typename ParticleType>
std::string DEM::PointSurface<ForceModel, ParticleType>::restart_data() const {
    using DEM::named_print;
    std::ostringstream ss;
    ss  << DEM::Surface<ForceModel, ParticleType>::restart_data() << ", "
        << named_print(points_.size(), "no_points") << ", " ;
    for (unsigned i = 0; i != points_.size(); ++i) {
        std::stringstream pname;
        pname << "p" << i;
        ss << named_print(points_[i], pname.str()) << ", ";
    }
    ss << named_print(infinite_, "infinite");
    return ss.str();
}
