//
// Created by erolsson on 2018-09-02.
//

#include "contact.h"

template<typename ForceModel, typename ParticleType>
DEM::Contact<ForceModel, ParticleType>::Contact(ParticleType* p1, ParticleType* p2, double increment) :
        p1_(p1),
        p2_(p2),
        surface_(nullptr),
        r2_(p1->get_radius() + p2->get_radius()),
        position_divider_(2),
        force_model_(p1, p2, increment),
        distance_function(&Contact::calculate_distance_vector_particle),
        tangential_function(&Contact::calculate_tangential_vector_particle)
{
    normal_ = calculate_distance_vector().normal();
}


template<typename ForceModel, typename ParticleType>
DEM::Contact<ForceModel, ParticleType>::Contact(ParticleType* p, SurfaceType* s,  double increment) :
        p1_(p),
        p2_(nullptr),
        surface_(s),
        r2_(p->get_radius()),
        position_divider_(1),
        force_model_(p, s, increment),
        distance_function(&Contact::calculate_distance_vector_surface),
        tangential_function(&Contact::calculate_tangential_vector_surface)

{
    normal_ = calculate_distance_vector().normal();
}


template<typename ForceModel, typename ParticleType>
DEM::void Contact<ForceModel, ParticleType>::update()
{
    auto distance_vector = calculate_distance_vector();
    double h = r2_ - distance_vector.length();
    normal_ = distance_vector.normalize();
    Vec3 dt = calculate_tangential_displacement_this_inc();
    force_model_.update(h, dt);
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 Contact<ForceModel, ParticleType>::calculate_distance_vector_particle() const
{
    return p1_->get_position() - p2_->get_position();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 Contact<ForceModel, ParticleType>::calculate_distance_vector_surface() const
{
    return surface_->vector_to_point(p1_->get_position());
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 Contact<ForceModel, ParticleType>::calculate_tangential_vector_particle() const
{
    Vec3 r1 = -normal_*(p1_->get_radius() - get_overlap()/2);
    Vec3 u1 = p1_->get_displacement_this_increment() + cross_product(p1_->get_rotation_this_increment(), r1);

    Vec3 r2 = normal_*(p2_->get_radius() - get_overlap()/2);
    Vec3 u2 = p2_->get_displacement_this_increment() + cross_product(p2_->get_rotation_this_increment(), r2);
    return (u1 - u2) - dot_product((u1 - u2), normal_)*normal_;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 Contact<ForceModel, ParticleType>::calculate_tangential_vector_surface() const
{
    Vec3 r1 = -normal_*(p1_->get_radius() - get_overlap()/2);
    Vec3 u1 = p1_->get_displacement_this_increment() + cross_product(p1_->get_rotation_this_increment(), r1);

    return u1 - dot_product(u1, normal_)*normal_;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 Contact<ForceModel, ParticleType>::position() const
{
    return p1_->get_position - normal_*(p1_->get_radius() - get_overlap()/position_divider_);
}
