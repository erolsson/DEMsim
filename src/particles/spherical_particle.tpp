//
// Created by erolsson on 2018-09-07.
//

#include "spherical_particle.h"

#include <sstream>

#include "../utilities/vec3.h"

template<typename ForceModel>
DEM::SphericalParticle<ForceModel>::SphericalParticle(double radius, const DEM::Vec3& position,
                                                      const DEM::Vec3& velocity, MaterialBase* material, unsigned id):
        ParticleBase<ForceModel>(4.*pi*radius*radius*radius/3*material->density, position, velocity, material, id),
        radius_(radius),
        inertia_(2*mass_*radius_*radius_/5),
        contacts_(ContactVector<std::pair<ContactPointerType, int> >())
{
    // Empty constructor
}


template<typename ForceModel>
void DEM::SphericalParticle<ForceModel>::move(const DEM::Vec3& new_disp_this_inc)
{
    position_ += new_disp_this_inc;
    displacement_this_inc_ = new_disp_this_inc;
}

template<typename ForceModel>
void DEM::SphericalParticle<ForceModel>::rotate(const DEM::Vec3& new_rot_this_inc)
{
    rot_this_inc_ = new_rot_this_inc;
    rot_ += new_rot_this_inc;
}

template<typename ForceModel>
double DEM::SphericalParticle<ForceModel>::kinetic_energy() const
{
    return translational_energy() + rotational_energy();
}

template<typename ForceModel>
double DEM::SphericalParticle<ForceModel>::translational_energy() const
{
    return mass_*dot_product(velocity_, velocity_)/2;
}

template<typename ForceModel>
double DEM::SphericalParticle<ForceModel>::rotational_energy() const
{
    return inertia_*dot_product(ang_vel_, ang_vel_)/2;
}

template<typename ForceModel>
void DEM::SphericalParticle<ForceModel>::reset_contacts() {
    number_of_contacts_ = 0;
    f_.set_zero();
    torque_.set_zero();
}

template<typename ForceModel>
void DEM::SphericalParticle<ForceModel>::sum_contact_forces()
{
    reset_contacts();
    auto contacts = contacts_.get_objects();
    for (const auto c_data : contacts){
        const auto c = c_data.first;
        const auto dir = c_data.second;
        f_ += (c->get_normal_force() + c->get_tangential_force())*dir;
        torque_ += (c->get_torque(position_))*dir;
    }
}

template<typename ForceModel>
std::size_t DEM::SphericalParticle<ForceModel>::number_of_contacts() const
{
    return 0;
}

template<typename ForceModel>
void DEM::SphericalParticle<ForceModel>::add_contact(DEM::SphericalParticle<ForceModel>::ContactPointerType contact,
                                                     std::size_t index, int direction)
{
    contacts_.insert(index, std::make_pair(contact, direction));
}

template<typename ForceModel>
void DEM::SphericalParticle<ForceModel>::remove_contact(std::size_t index)
{
    contacts_.erase(index);
}

template<typename ForceModel>
std::string DEM::SphericalParticle<ForceModel>::get_output_string() const
{
    std::stringstream ss;
    ss << id_ << ", " << position_.x() << ", " << position_.y() << ", " << position_.z() << ", ";
    ss << rot_.x() << ", " << rot_.y() << ", " << rot_.z() << ", " << radius_<< ", ";
    ss << kinetic_energy() << ", " << material_->id;
    return ss.str();
}




