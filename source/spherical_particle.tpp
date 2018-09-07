//
// Created by erolsson on 2018-09-07.
//

#include "spherical_particle.h"

#include <sstream>

#include "vec3.h"

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
    return mass_*dot_product(velocity_, velocity_)/2 + inertia_*dot_product(ang_vel_, ang_vel_)/2;
}

template<typename ForceModel>
void DEM::SphericalParticle<ForceModel>::reset_contacts() {
    number_of_contacts_ = 0;
    fn_ = 0*fn_;
    ft_ = 0*ft_;
    torque_ = 0*torque_;
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
    ss << id_ << ", " << position_.x << ", " << position_.y << ", " << position_.z << ", ";
    ss << rot_.x << ", " << rot_.y << ", " << rot_.z << ", " << radius_<< ", ";
    ss << kinetic_energy() << ", " << material_->id;
    return ss.str();
}
