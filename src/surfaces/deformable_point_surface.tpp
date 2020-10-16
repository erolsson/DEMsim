//
// Created by erolsson on 14/09/2020.
//

#include "deformable_point_surface.h"

#include "../utilities/vec3.h"

template<typename ForceModel, typename ParticleType>
DEM::DeformablePointSurface<ForceModel, ParticleType>::DeformablePointSurface(std::size_t id,
                                                                              const std::vector<Vec3>& points,
                                                                              bool infinite, const std::string& name,
                                                                              bool adhesive, std::size_t collision_id):
    DEM::PointSurface<ForceModel, ParticleType>(id, points, infinite, name, adhesive, collision_id),
    nodal_displacements_(std::vector<Vec3>(points.size()))
{
    if (points.size() != 4) {
        throw std::invalid_argument("Only 4-node surfaces are currently supported");
    }
}

template<typename ForceModel, typename ParticleType>
DEM::DeformablePointSurface<ForceModel, ParticleType>::DeformablePointSurface(const ParameterMap& parameters):
    PointSurface<ForceModel, ParticleType>(parameters) ,
    strain_x_(parameters.get_parameter("strain_x")),
    strain_y_(parameters.get_parameter("strain_y"))
{
    auto no_points = parameters.get_parameter<std::size_t>("no_points");
    for (unsigned i = 0; i != no_points; ++i) {
        std::ostringstream dp;
        dp << "d" << i;
        nodal_displacements_.push_back(parameters.get_vec3(dp.str()));
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::DeformablePointSurface<ForceModel, ParticleType>::deform(const std::vector<Vec3>& nodal_displacements)
{
    for (unsigned i = 0; i != nodal_displacements.size(); ++i) {
        points_[i] += nodal_displacements[i];
        nodal_displacements_[i] = nodal_displacements[i];
    }
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3
DEM::DeformablePointSurface<ForceModel, ParticleType>::get_displacement_this_increment(const DEM::Vec3& position) const
{
    auto x_axis = (points_[1] - points_[0]).normalize();
    auto y_axis = -cross_product(x_axis, get_normal());

    // Calculating the mid point
    auto mid_point = Vec3 {};
    for (const auto& p: points_) {
        mid_point += p;
    }
    mid_point /= points_.size();

    double x = dot_product(position - mid_point, x_axis);
    double y = dot_product(position - mid_point, y_axis);

    // Calculating the nodal positions
    std::vector<double> nodal_coords_x;
    std::vector<double> nodal_coords_y;
    for (const auto& point: points_) {
        nodal_coords_x.push_back(dot_product(point - mid_point, x_axis));
        nodal_coords_y.push_back(dot_product(point - mid_point, y_axis));
    }

    // Assuming a rectangular surface, could be changed here later
    double xi = 2*x/(nodal_coords_x[2] - nodal_coords_x[0]);
    double eta = 2*y/(nodal_coords_y[2] - nodal_coords_y[0]);

    std::array<double, 4> N = { 1./4*(1-xi)*(1-eta), 1./4*(1+xi)*(1-eta),
                                1./4*(1+xi)*(1+eta), 1./4*(1-xi)*(1+eta)};

    // Calculating the in_plane nodal displacements and from that the in-plane displacement at point
    double ux = 0;
    double uy = 0;
    for (unsigned i = 0; i!= 4; ++i) {
        ux += dot_product(nodal_displacements_[i], x_axis)*N[i];
        uy += dot_product(nodal_displacements_[i], y_axis)*N[i];
    }
    Vec3 u = ux*x_axis + uy*y_axis;
    return DEM::PointSurface<ForceModel, ParticleType>::get_displacement_this_increment(position) + u;
}

template<typename ForceModel, typename ParticleType>
std::string DEM::DeformablePointSurface<ForceModel, ParticleType>::get_output_string() const
{
    std::ostringstream stream;
    stream << DEM::PointSurface<ForceModel, ParticleType>::get_output_string();
    for (const auto& dp: nodal_displacements_) {
        stream << ", " << dp.x() << ", " << dp.y() << ", " << dp.z();
    }
    return stream.str();
}

template<typename ForceModel, typename ParticleType>
std::string DEM::DeformablePointSurface<ForceModel, ParticleType>::restart_data() const
{
    using DEM::named_print;
    std::ostringstream ss;
    ss  << DEM::PointSurface<ForceModel, ParticleType>::restart_data() << ", ";
    for (unsigned i = 0; i != nodal_displacements_.size(); ++i) {
        std::stringstream pname;
        pname << "d" << i;
        ss << named_print(nodal_displacements_[i], pname.str()) << ", ";
    }
    ss << named_print(strain_x_, "strain_x") << ", "
       << named_print(strain_y_, "strain_y");
    return ss.str();
}

template<typename ForceModel, typename ParticleType>
void DEM::DeformablePointSurface<ForceModel, ParticleType>::deform(std::chrono::duration<double> time_increment) {
    auto x_axis = (points_[1] - points_[0]).normalize();
    auto y_axis = -cross_product(x_axis, get_normal());

    double x = points_[2].x() - points_[0].x();
    double y = points_[2].y() - points_[0].y();

    // The nodal displacements if each nodes moves -u/2 and u/2 the total deformation is u
    double ux = x*strain_x_/2*time_increment.count();
    double uy = y*strain_y_/2*time_increment.count();
    nodal_displacements_ = {-ux*x_axis - uy*y_axis,  ux*x_axis - uy*y_axis,
                             ux*x_axis + uy*y_axis, -ux*x_axis + uy*y_axis};
        
    for (unsigned i = 0; i != nodal_displacements_.size(); ++i) {
        points_[i] += nodal_displacements_[i];
    }
    
    normal_ = calculate_normal();
    update_bounding_box();
}


