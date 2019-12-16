//
// Created by erolsson on 2018-09-02.
//

#include "contact.h"

#include <fstream>
#include <sstream>

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
void DEM::Contact<ForceModel, ParticleType>::update()
{
    auto distance_vector = calculate_distance_vector();
    auto h = r2_ - distance_vector.length();
    normal_ = distance_vector.normalize();
    Vec3 dt(0, 0, 0);
    Vec3 w(0, 0, 0);
    if (h > 0) {
        dt = calculate_tangential_displacement_this_inc();
        w = calculate_rotation_this_inc();
    }
    force_model_.update(h, dt, w, normal_);
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::get_torque(const DEM::Vec3 &point) const {
    return cross_product(get_position() - point, get_tangential_force()) + force_model_.get_rolling_resistance_torque();
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
    Vec3 r1 = -normal_*(p1_->get_radius() - get_overlap()/2);
    Vec3 u1 = p1_->get_displacement_this_increment() + cross_product(p1_->get_rotation_this_increment(), r1);

    Vec3 r2 = normal_*(p2_->get_radius() - get_overlap()/2);
    Vec3 u2 = p2_->get_displacement_this_increment() + cross_product(p2_->get_rotation_this_increment(), r2);
    return (u1 - u2) - dot_product((u1 - u2), normal_)*normal_;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::calculate_tangential_vector_surface() const
{
    Vec3 r1 = -normal_*(p1_->get_radius() - get_overlap()/2);
    Vec3 u1 = p1_->get_displacement_this_increment() + cross_product(p1_->get_rotation_this_increment(), r1);

    return u1 - dot_product(u1, normal_)*normal_;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::get_position() const
{
    return p1_->get_position() - normal_*(p1_->get_radius() - get_overlap()/position_divider_);
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::calculate_rotation_vector_particle() const {
    return p1_->get_rotation_this_increment() - p2_->get_rotation_this_increment();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Contact<ForceModel, ParticleType>::calculate_rotation_vector_surface() const {
    return p1_->get_rotation_this_increment();
}

template<typename ForceModel, typename ParticleType>
std::string DEM::Contact<ForceModel, ParticleType>::get_output_string() const {
    std::stringstream ss;
    ss << p1_->get_id() << ", ";
    if (p2_ != nullptr) {
        ss << p2_->get_id();
    }
    else {
        ss << surface_->get_id();
    }
    ss << ", " << get_overlap() << ", " << force_model_.get_output_string();
    return ss.str();
}

