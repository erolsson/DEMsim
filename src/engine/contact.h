//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_CONTACT_H
#define DEMSIM_CONTACT_H

#include <chrono>
#include <memory>

#include "Eigen/Dense"

#include "../surfaces/surface_base.h"


namespace DEM {
    class ParameterMap;
    template<typename ForceModel, typename ParticleType> class Surface;

    template<typename ForceModel, typename ParticleType>
    class Contact {
        using ContactPointerType = Contact<ForceModel, ParticleType>*;
        using SurfaceType = Surface<ForceModel, ParticleType>;

    public:
        //Constructors
        Contact(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double> increment);
        Contact(ParticleType* particle1, SurfaceType* surface,  std::chrono::duration<double> increment);

        Contact(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double> increment,
                const ParameterMap& parameters);
        Contact(ParticleType* particle1, SurfaceType* surface,  std::chrono::duration<double> increment,
                const ParameterMap& parameters);

        void update();
        [[nodiscard]] Vec3 get_position() const;

        [[nodiscard]] Vec3 get_normal_force() const { return force_model_.get_normal_force()*get_normal(); }
        [[nodiscard]] Vec3 get_tangential_force() const { return force_model_.get_tangential_force(); }
        [[nodiscard]] Vec3 get_torque(int direction) const;   // 1 for p1 -1 for p2

        std::pair<ParticleType*, ParticleType*> get_particles() const  { return std::make_pair(p1_, p2_); }
        const SurfaceType* get_surface() const { return surface_; }

        [[nodiscard]] bool active() const { return force_model_.active();  }
        [[nodiscard]] double get_overlap() const { return force_model_.get_overlap(); }
        [[nodiscard]] double get_contact_radius() const { return force_model_.get_contact_area();}
        [[nodiscard]] const Vec3& get_normal() const { return normal_; }

        [[nodiscard]] std::pair<std::size_t, size_t> get_id_pair() const;

        void set_increment(std::chrono::duration<double> increment) {force_model_.set_increment(increment); }
        [[nodiscard]] std::string get_output_string() const;
        [[nodiscard]] std::string restart_data() const;
        Eigen::Matrix<double, 3, 3> get_force_fabric_tensor() const;

    private:
        ParticleType* const p1_;
        ParticleType* const p2_;
        SurfaceType*  const surface_;

        double r2_;                       // r2 - distance between particles or particle plane is overlap
        char position_divider_;
        Vec3 normal_;                     // Contact plane normal, same direction as normal_force_

        ForceModel force_model_;

        Vec3 (Contact<ForceModel, ParticleType>::*distance_function)() const;
        Vec3 (Contact<ForceModel, ParticleType>::*tangential_function)() const;
        Vec3 (Contact<ForceModel, ParticleType>::*rotation_function)() const;

        [[nodiscard]] Vec3 calculate_distance_vector() const { return (this->*distance_function)(); }
        [[nodiscard]] Vec3 calculate_tangential_displacement_this_inc() const { return (this->*tangential_function)(); };
        [[nodiscard]] Vec3 calculate_rotation_this_inc() const { return (this->*rotation_function)(); }

        [[nodiscard]] Vec3 calculate_distance_vector_particle() const;
        [[nodiscard]] Vec3 calculate_distance_vector_surface() const;

        [[nodiscard]] Vec3 calculate_tangential_vector_particle() const;
        [[nodiscard]] Vec3 calculate_tangential_vector_surface() const;

        [[nodiscard]] Vec3 calculate_rotation_vector_particle() const;
        [[nodiscard]] Vec3 calculate_rotation_vector_surface() const;
    };
}

#include "contact.tpp"



#endif //DEMSIM_CONTACT_H
