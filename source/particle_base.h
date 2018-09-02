//
// Created by erolsson on 2018-07-29.
//

#ifndef DEMSIM_PARTICLE_BASE_NEW_H
#define DEMSIM_PARTICLE_BASE_NEW_H

#include <memory>
#include <string>
#include <vector>

#include "material_base.h"
#include "vec3.h"

namespace DEM {

    template<typename ForceModel>
    class ParticleBase {
    public:
        ParticleBase(const ParticleBase&) = delete;
        ParticleBase& operator=(const ParticleBase&) = delete;
        ParticleBase(double, const Vec3&, const Vec3&, MaterialBase*, unsigned );

        unsigned get_id() const { return id_; }
        const MaterialBase* get_material() const { return material_; }

        const Vec3& get_normal_force() const { return fn_; }
        const Vec3& get_tangential_force() const { return ft_; }
        const Vec3& get_torque() const { return torque_; }

        double get_mass() const { return mass_; }


        const Vec3& get_position() const { return position_; }
        const Vec3& get_velocity_() const { return velocity_; }
        void set_velocity(const Vec3& new_velocity) { velocity_ = new_velocity; }

        const Vec3& get_rotation() const { return rot_; }
        const Vec3& get_angular_velocity() const { return ang_vel_; }
        void set_angular_velocity(Vec3 new_ang_vel) { ang_vel_ = new_ang_vel; }

        const Vec3& get_displacement_this_increment() const { return displacement_this_inc_; }
        const Vec3& get_rotation_this_increment() const { return rot_this_inc_; }

    protected:
        const unsigned id_;
        double mass_;
        const MaterialBase* material_;

        Vec3 position_;
        Vec3 velocity_;
        Vec3 rot_{ Vec3(0., 0., 0.) };
        Vec3 ang_vel_{ Vec3(0., 0., 0.) };

        // Forces
        Vec3 fn_{ Vec3(0., 0., 0.) };
        Vec3 ft_{ Vec3(0., 0., 0.) };
        Vec3 torque_{ Vec3(0., 0., 0.) };

        // Needed for sticking friction model
        Vec3 displacement_this_inc_{ Vec3(0., 0., 0.)};
        Vec3 rot_this_inc_{ Vec3(0., 0., 0.)};

        std::vector<Vec3> contact_forces_ { std::vector<Vec3>() };
        unsigned number_of_contacts_ { 0 };
    };
}

#include "particle_base.tpp"

#endif //DEMSIM_PARTICLE_BASE_NEW_H
