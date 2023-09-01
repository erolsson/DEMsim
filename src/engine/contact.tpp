//
// Created by erolsson on 2018-09-02.
//

#include "contact.h"

#include <fstream>
#include <sstream>

#include "../utilities/printing_functions.h"
#include "../utilities/vec3.h"

template<typename ForceModel, typename ParticleType>
DEM::Contact<ForceModel, ParticleType>::Contact(ParticleType* particle1, ParticleType* particle2,
                                                std::chrono::duration<double> increment) :
    p1_(particle1),
    p2_(particle2),
    surface_(nullptr),
    r2_(particle1->get_radius() + particle2->get_radius()),
    position_divider_(2),
    force_model_(particle1, particle2, increment),
    distance_function(&Contact::calculate_distance_vector_particle),
    tangential_function(&Contact::calculate_tangential_vector_particle),
    rotation_function(&Contact::calculate_rotation_vector_particle)
{
    normal_ = calculate_distance_vector().normal();
}
template<typename ForceModel, typename ParticleType>
DEM::Contact<ForceModel, ParticleType>::Contact(ParticleType* particle1, SurfaceType* surface,
                                                std::chrono::duration<double> increment) :
    p1_(particle1),
    p2_(nullptr),
    surface_(surface),
    r2_(particle1->get_radius()),
    position_divider_(1),
    force_model_(particle1, surface, increment),
    distance_function(&Contact::calculate_distance_vector_surface),
    tangential_function(&Contact::calculate_tangential_vector_surface),
    rotation_function(&Contact::calculate_rotation_vector_surface)
{
    normal_ = calculate_distance_vector().normal();
}

template<typename ForceModel, typename ParticleType>
DEM::Contact<ForceModel, ParticleType>::Contact(ParticleType *particle1, ParticleType *particle2,
                                                std::chrono::duration<double> increment,
                                                const DEM::ParameterMap& parameters) :
    p1_(particle1),
    p2_(particle2),
    surface_(nullptr),
    r2_(particle1->get_radius() + particle2->get_radius()),
    position_divider_(2),
    force_model_(particle1, particle2, increment, parameters),
    distance_function(&Contact::calculate_distance_vector_particle),
    tangential_function(&Contact::calculate_tangential_vector_particle),
    rotation_function(&Contact::calculate_rotation_vector_particle)
{
    normal_ = calculate_distance_vector().normal();
}

template<typename ForceModel, typename ParticleType>
DEM::Contact<ForceModel, ParticleType>::Contact(ParticleType *particle1, SurfaceType *surface,
                                                std::chrono::duration<double> increment,
                                                const DEM::ParameterMap& parameters) :
    p1_(particle1),
    p2_(nullptr),
    surface_(surface),
    r2_(particle1->get_radius()),
    position_divider_(1),
    force_model_(particle1, surface, increment, parameters),
    distance_function(&Contact::calculate_distance_vector_surface),
    tangential_function(&Contact::calculate_tangential_vector_surface),
    rotation_function(&Contact::calculate_rotation_vector_surface)
{
    normal_ = calculate_distance_vector().normal();
}

template<typename ForceModel, typename ParticleType>
void DEM::Contact<ForceModel, ParticleType>::update()
{
    auto distance_vector = calculate_distance_vector();
    auto h = r2_ - distance_vector.length();
    normal_ = distance_vector.normal();
    Vec3 dt(0, 0, 0);
    Vec3 w(0, 0, 0);
    
    dt = calculate_tangential_displacement_this_inc();
    w = calculate_rotation_this_inc();
    
    force_model_.update(h, dt, w, get_normal());
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::get_torque(int direction) const
{
    Vec3 point;
    if (direction == 1) {
        point = p1_->get_position();
    }
    else if (direction == -1) {
        point = p2_->get_position();
    }
    return cross_product(get_position() - point, get_tangential_force())
        + force_model_.get_rolling_resistance_torque();
}


template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::calculate_distance_vector_particle() const
{
    return p1_->get_position() - p2_->get_position();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::calculate_distance_vector_surface() const
{
    return surface_->vector_to_point(p1_->get_position());
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::calculate_tangential_vector_particle() const
{
    Vec3 r1 = -get_normal()*(p1_->get_radius() - get_overlap()/2);
    Vec3 u1 = p1_->get_displacement_this_increment() + cross_product(p1_->get_rotation_this_increment(), r1);

    Vec3 r2 = get_normal()*(p2_->get_radius() - get_overlap()/2);
    Vec3 u2 = p2_->get_displacement_this_increment() + cross_product(p2_->get_rotation_this_increment(), r2);
    return (u1 - u2) - dot_product((u1 - u2), get_normal())*get_normal();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::calculate_tangential_vector_surface() const
{
    Vec3 r1 = -get_normal()*(p1_->get_radius() - get_overlap()/2);
    Vec3 u1 = p1_->get_displacement_this_increment() + cross_product(p1_->get_rotation_this_increment(), r1);

    Vec3 u2 = surface_->get_displacement_this_increment(p1_->get_position() + r1);

    return (u1 - u2) - dot_product(u1 - u2, get_normal())*get_normal();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::get_position() const
{
    return p1_->get_position() - get_normal()*(p1_->get_radius() - get_overlap()/position_divider_);
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::calculate_rotation_vector_particle() const
{
    return p1_->get_rotation_this_increment() - p2_->get_rotation_this_increment();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::calculate_rotation_vector_surface() const
{
    return p1_->get_rotation_this_increment();
}

template<typename ForceModel, typename ParticleType>
std::pair<std::size_t, size_t> DEM::Contact<ForceModel, ParticleType>::get_id_pair() const
{
    std::size_t id2;
    if (surface_ != nullptr) {
        id2 = surface_->get_id();
    }
    else {
        id2 = p2_->get_id();
    }
    return std::make_pair(p1_->get_id(), id2);
}


template<typename ForceModel, typename ParticleType>
std::string DEM::Contact<ForceModel, ParticleType>::get_output_string() const
{
    std::stringstream ss;
    ss << p1_->get_id() << ", ";
    if (p2_ != nullptr) {
        ss << p2_->get_id();
    }
    else {
        ss << surface_->get_id();
    }
    ss << ", " << get_normal().x() << ", " << get_normal().y() << ", " << get_normal().z();
    ss << ", " << get_overlap() << ", " << force_model_.get_output_string();
    return ss.str();
}

template<typename ForceModel, typename ParticleType>
std::string DEM::Contact<ForceModel, ParticleType>::restart_data() const
{
    std::stringstream ss;
    using DEM::named_print;
    if (p2_ != nullptr) {
        ss << named_print("particle-particle", "type") << ", "
           << named_print(p1_->get_id(), "object1") << ", " << named_print(p2_->get_id(), "object2") << ", ";
    }
    else {
        ss << named_print("particle-surface", "type") << ", "
           << named_print(p1_->get_id(), "object1") << ", " << named_print(surface_->get_id(), "object2") << ", ";
    }
    ss << force_model_.restart_data();
    return ss.str();
}

template<typename ForceModel, typename ParticleType>
Eigen::Matrix<double, 3, 3> DEM::Contact<ForceModel, ParticleType>::get_force_fabric_tensor() const
{
    Eigen::Matrix<double, 3, 3> tensor = Eigen::Matrix<double, 3, 3>::Zero();
    double d = calculate_distance_vector().length();
    Vec3 fn = get_normal_force();
    Vec3 ft = get_tangential_force();
    for (unsigned i = 0; i != 3; ++i) {
        for (unsigned j = 0; j != 3; ++j) {
            tensor(i, j) += d*fn[i]*get_normal()[j] + d*get_normal()[i]*ft[j];
        }
    }
    return tensor;
}

template<typename ForceModel, typename ParticleType>
void DEM::Contact<ForceModel, ParticleType>::assign_new_contact_particles(ParticleType* p1, ParticleType* p2)
{
    p1_ = p1;
    p2_ = p2;
}

template<typename ForceModel, typename ParticleType>
void DEM::Contact<ForceModel, ParticleType>::assign_new_contact_particles(ParticleType* p1, SurfaceType* s){
    p1_ = p1;
    surface_ = s;
}





