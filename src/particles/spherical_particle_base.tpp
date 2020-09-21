//
// Created by erolsson on 16/05/19.
//

#include "spherical_particle_base.h"

#include <sstream>

#include "../utilities/file_reading_functions.h"
#include "../utilities/printing_functions.h"
#include "../utilities/vec3.h"

template<typename ForceModel>
DEM::SphericalParticleBase<ForceModel>::SphericalParticleBase(double radius, const DEM::Vec3& position,
                                                              const DEM::Vec3& velocity, const MaterialBase* material,
                                                              std::size_t object_id, std::size_t collision_id):
        ParticleBase<ForceModel>(4.*3.1415*radius*radius*radius/3*material->density, position, velocity, material,
                                 object_id),
        collision_id_(collision_id),
        radius_(radius),
        inertia_(2*mass_*radius_*radius_/5)
{
    // Empty constructor
}

template<typename ForceModel>
DEM::SphericalParticleBase<ForceModel>::SphericalParticleBase(const DEM::ParameterMap& parameters,
                                                              DEM::MaterialBase *material):
    ParticleBase<ForceModel>(parameters, material),
    collision_id_(parameters.get_parameter<std::size_t>("collision_id")),
    radius_(parameters.get_parameter<double>("r")),
    inertia_(parameters.get_parameter<double>("I"))
{
    // Empty constructor
}


template<typename ForceModel>
void DEM::SphericalParticleBase<ForceModel>::move(const DEM::Vec3& new_disp_this_inc)
{
    position_ += new_disp_this_inc;
    displacement_this_inc_ = new_disp_this_inc;
}

template<typename ForceModel>
void DEM::SphericalParticleBase<ForceModel>::rotate(const DEM::Vec3& new_rot_this_inc)
{
    rot_this_inc_ = new_rot_this_inc;
    rot_ += new_rot_this_inc;
}

template<typename ForceModel>
double DEM::SphericalParticleBase<ForceModel>::kinetic_energy() const
{
    return translational_energy() + rotational_energy();
}

template<typename ForceModel>
double DEM::SphericalParticleBase<ForceModel>::translational_energy() const
{
    return mass_*dot_product(velocity_, velocity_)/2;
}

template<typename ForceModel>
double DEM::SphericalParticleBase<ForceModel>::rotational_energy() const
{
    return inertia_*dot_product(ang_vel_, ang_vel_)/2;
}

template<typename ForceModel>
void DEM::SphericalParticleBase<ForceModel>::reset_contacts() {
    number_of_contacts_ = 0;
    f_.set_zero();
    torque_.set_zero();
}

template<typename ForceModel>
template<typename T>
void DEM::SphericalParticleBase<ForceModel>::sum_contact_forces(const T& contacts)
{
    reset_contacts();
    for (const auto c_data : contacts.get_objects()){
        const auto c = c_data.first;
        const auto dir = c_data.second;
        f_ += (c->get_normal_force() + c->get_tangential_force())*dir;
        torque_ += (c->get_torque(dir))*dir;
    }
}

template<typename ForceModel>
std::size_t DEM::SphericalParticleBase<ForceModel>::number_of_contacts() const
{
    return 0;
}

template<typename ForceModel>
std::string DEM::SphericalParticleBase<ForceModel>::get_output_string() const
{
    std::ostringstream ss;
    ss << id_ << ", " << position_.x() << ", " << position_.y() << ", " << position_.z() << ", ";
    ss << rot_.x() << ", " << rot_.y() << ", " << rot_.z() << ", " << radius_<< ", ";
    ss << kinetic_energy() << ", " << material_->id << ", " << f_.x() << ", " << f_.y() << ", " << f_.z();
    return ss.str();
}

template<typename ForceModel>
template<typename ContactType, typename VectorType>
void DEM::SphericalParticleBase<ForceModel>::add_contact(ContactType contact, std::size_t index, int direction,
                                                         VectorType& contacts) {
    contacts.insert(index, std::make_pair(contact, direction));
}

template<typename ForceModel>
template<typename VectorType>
void DEM::SphericalParticleBase<ForceModel>::remove_contact(std::size_t index, VectorType& contacts) {
    contacts.erase(index);
}

template<typename ForceModel>
std::string DEM::SphericalParticleBase<ForceModel>::restart_data() const {
    using DEM::named_print;
    std::ostringstream ss;
    ss << ParticleBase<ForceModel>::restart_data() << ", "
       << named_print(collision_id_, "collision_id") << ", "
       << named_print(radius_, "r") << ", "
       << named_print(inertia_, "I");
    return ss.str();
}

