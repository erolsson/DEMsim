//
// Created by erolsson on 2018-09-02.
//

#include "surface_base.h"

#include <sstream>
#include <string>
#include <utility>

#include "../utilities/file_reading_functions.h"
#include "../utilities/printing_functions.h"

template<typename ForceModel, typename ParticleType>
DEM::Surface<ForceModel, ParticleType>::Surface(std::size_t id, std::size_t collision_id, std::string  name,
                                                bool adhesive) :
        object_id_(id), collision_id_(collision_id), name_(std::move(name)), adhesive_(adhesive)
{
    //Empty constructor
}

template<typename ForceModel, typename ParticleType>
DEM::Surface<ForceModel, ParticleType>::Surface(const DEM::ParameterMap& parameters) :
    velocity_{ parameters.get_vec3("v") },
    acceleration_{ parameters.get_vec3("a")},
    displacement_this_inc_{parameters.get_vec3("disp_this_inc") },
    rotation_this_inc_{ parameters.get_vec3("rot_this_inc")},
    rotation_point_{ parameters.get_vec3("rot_point")},
    object_id_(parameters.get_parameter<std::size_t>("id")),
    collision_id_(parameters.get_parameter<std::size_t>("collision_id")),
    name_(parameters.get_parameter<std::string>("name")),
    mass_(parameters.get_parameter<double>("mass")),
    adhesive_(false),
    force_control_amplitudes_ {nullptr, nullptr, nullptr}
{
    if (parameters.exist("adhesive")) {
        adhesive_ = parameters.get_parameter<bool>("adhesive");
        std::string directions = "xyz";
        for (std::size_t i = 0; i != force_control_amplitudes_.size(); ++i) {
            if (parameters.get_parameter<bool>("force_control_" + directions.substr(i, i+1))) {
                std::cout << "Warning: force control specified in direction " << directions[i]
                          << "for surface " << object_id_ << "\n";
                std::cout << "The force control has to be specified manually in the restart file\n";
            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Surface<ForceModel, ParticleType>::get_tangential_displacement_this_inc(const Vec3& point) const
{
    return get_displacement_this_increment(point) -
           dot_product(get_displacement_this_increment(point), get_normal(point))*get_normal(point);
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
    for (auto const& c : contacts_.get_objects()) {
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
    for (auto const& c : contacts_.get_objects()) {
        Vec3 F = c->get_normal_force();
        if (!F.is_zero()) {
            normal_force += dot_product(c->get_normal_force(), get_normal(c->get_particles().first->get_position()));
        }
    }
    return normal_force;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Surface<ForceModel, ParticleType>::get_tangential_force() const
{
    Vec3 tangential_force = Vec3(0, 0, 0);
    for (auto const& c : contacts_.get_objects()) {
        tangential_force += c->get_tangential_force();
    }
    return -1*tangential_force;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::Surface<ForceModel, ParticleType>::get_total_force() const
{
    Vec3 force = Vec3(0, 0, 0);
    for (auto const& c : contacts_.get_objects()) {
        force += c->get_normal_force() + c->get_tangential_force();
    }
    return -force;
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
        force_control_amplitudes_[0] = amplitude;
    }
    if( direction == 'y') {
        force_control_amplitudes_[1] = amplitude;
    }
    if( direction == 'z') {
        force_control_amplitudes_[2] = amplitude;
    }
}

template<typename ForceModel, typename ParticleType>
std::string DEM::Surface<ForceModel, ParticleType>::restart_data() const {
    using DEM::named_print;
    std::ostringstream ss;
    ss << named_print(object_id_, "id") << ", "
       << named_print(collision_id_, "collision_id") << ", "
       << named_print(type(), "type") << ", "
       << named_print(name_, "name") << ", "
       << named_print(mass_, "mass") << ", "
       << named_print(velocity_, "v") << ", "
       << named_print(acceleration_, "a") << ", "
       << named_print(displacement_this_inc_, "disp_this_inc") << ", "
       << named_print(rotation_this_inc_, "rot_this_inc") << ", "
       << named_print(rotation_point_, "rot_point");
    std::string directions = "xyz";
    for (unsigned i = 0; i != force_control_amplitudes_.size(); ++i) {
        ss << ", force_control_" << directions[i] << "=";
        if (force_control_amplitudes_[i] == nullptr) {
            ss << 0;
        }
        else {
            ss << 1;
        }
    }
    return ss.str();
}
