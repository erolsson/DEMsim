//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_CONTACT_H
#define DEMSIM_CONTACT_H

#include <chrono>
#include <memory>

#include "surface_base.h"


namespace DEM{
    template<typename ForceModel, typename ParticleType> class Surface;

    template<typename ForceModel, typename ParticleType>
    class Contact {
        using ContactPointerType = Contact<ForceModel, ParticleType>*;
        using SurfaceType = Surface<ForceModel, ParticleType>;

    public:
        //Constructors
        Contact(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double> increment);
        Contact(ParticleType* particle1, SurfaceType* surface,  std::chrono::duration<double> increment);

        void update();
        Vec3 position() const;

        Vec3 get_normal_force() const { return force_model_.get_normal_force()*get_normal(); }
        Vec3 get_tangential_force() const { return force_model_.get_tangential_force(); }
        Vec3 get_torque(const Vec3& point) const { return cross_product(position() - point, get_tangential_force()); }

        std::pair<ParticleType*, ParticleType*> get_particles() const  { return std::make_pair(p1_, p2_); }
        const SurfaceType* get_surface() const { return surface_; }

        bool active() const { return force_model_.active();  }
        double get_overlap() const { return force_model_.get_overlap(); }
        double get_contact_radius() const { return force_model_.get_contact_area();}
        const Vec3& get_normal() const { return normal_; }

        void set_increment(std::chrono::duration<double> increment) {force_model_.set_increment(increment); }

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

        Vec3 calculate_distance_vector() const { return (this->*distance_function)(); }
        Vec3 calculate_tangential_displacement_this_inc() const { return (this->*tangential_function)(); };

        Vec3 calculate_distance_vector_particle() const;
        Vec3 calculate_distance_vector_surface() const;

        Vec3 calculate_tangential_vector_particle() const;
        Vec3 calculate_tangential_vector_surface() const;

    };
}

#include "contact.tpp"


#endif //DEMSIM_CONTACT_H
