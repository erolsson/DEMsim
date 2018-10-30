//
// Created by erolsson on 2018-09-02.
//

#include "surface_base.h"

template<typename ForceModel, typename ParticleType>
DEM::Surface<ForceModel, ParticleType>::Surface(std::size_t id) :
        id_(id)
{
    //Empty constructor
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Surface<ForceModel, ParticleType>::get_tangential_displacement_this_inc(const Vec3& point) const
{
    return displacement_this_inc(point)-
            dot_product(displacement_this_inc(point), get_normal(point))*get_normal(point);
}

template<typename ForceModel, typename ParticleType>
void DEM::Surface<ForceModel, ParticleType>::rest()
{
    velocity_ *= 0;
    rotation_this_inc_ *= 0;
    rotation_point_ *= 0;
}

template<typename ForceModel, typename ParticleType>
std::vector<ParticleType*> DEM::Surface<ForceModel, ParticleType>::get_contacting_particles() const
{
    std::vector<ParticleType*> contacting_particles;
    for (auto const& c : contacts_) {
        //I'm not proud over the cast but it works, at least complies
        //    if(c->active())
        contacting_particles.push_back(const_cast<ParticleType*>(c->get_particles().first));
    }
    return contacting_particles;
}

template<typename ForceModel, typename ParticleType>
double DEM::Surface<ForceModel, ParticleType>::get_normal_force() const
{
    double normal_force = 0.0;
    for (auto const& c : contacts_) {
        normal_force += dot_product(c->get_normal_force(), get_normal(c->get_position()));
    }
    return normal_force;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Surface<ForceModel, ParticleType>::get_tangential_force() const
{
    Vec3 tangential_force = Vec3(0, 0, 0);
    for (auto const& c : contacts_) {
        tangential_force += c->get_tangential_force();
    }
    return -1*tangential_force;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Surface<ForceModel, ParticleType>::get_total_force() const
{
    Vec3 force = Vec3(0, 0, 0);
    for (auto const& c : contacts_) {
        force += c->get_normal_force()+c->get_tangential_force();
    }
    return force;
}

template<typename ForceModel, typename ParticleType>
void DEM::Surface<ForceModel, ParticleType>::add_contact(ContactPointerType contact, size_t index_of_other_object)
{
    contacts_.insert(index_of_other_object, contact);
}

template<typename ForceModel, typename ParticleType>
void DEM::Surface<ForceModel, ParticleType>::remove_contact(size_t index_of_other_object)
{
    contacts_.erase(index_of_other_object);
}

template<typename ForceModel, typename ParticleType>
void DEM::Surface<ForceModel, ParticleType>::set_force_amplitude(DEM::Surface<ForceModel, ParticleType>::ForceAmpPtr
    amplitude, char direction)
{
    if( direction == 'x') {
        force_control_amplitudes_[0] == amplitude;
    }
    if( direction == 'y') {
        force_control_amplitudes_[1] == amplitude;
    }
    if( direction == 'z') {
        force_control_amplitudes_[2] == amplitude;
    }
}
