//
// Created by erolsson on 16/05/19.
//

#ifndef DEMSIM_SPHERICAL_PARTICLE_BASE_H
#define DEMSIM_SPHERICAL_PARTICLE_BASE_H

#include "particle_base.h"

namespace DEM {
    template<typename ForceModel>
    class SphericalParticleBase : public ParticleBase<ForceModel> {
        using ParticleBase<ForceModel>::id_;

        using ParticleBase<ForceModel>::position_;
        using ParticleBase<ForceModel>::displacement_this_inc_;
        using ParticleBase<ForceModel>::rot_this_inc_;
        using ParticleBase<ForceModel>::rot_;
        using ParticleBase<ForceModel>::mass_;
        using ParticleBase<ForceModel>::velocity_;
        using ParticleBase<ForceModel>::ang_vel_;

        using ParticleBase<ForceModel>::number_of_contacts_;
        using ParticleBase<ForceModel>::f_;
        using ParticleBase<ForceModel>::torque_;
    protected:
        using ParticleBase<ForceModel>::material_;
    public:

        // No assignment of particles and no plain copies

        SphericalParticleBase(double radius, const Vec3& position, const Vec3& velocity, MaterialBase* material,
                          unsigned id);

        SphericalParticleBase(const SphericalParticleBase&) = delete;
        SphericalParticleBase& operator=(const SphericalParticleBase&) = delete;

        [[nodiscard]] double get_radius() const { return radius_; }
        [[nodiscard]] double get_inertia() const { return inertia_; }

        void move(const Vec3& new_disp_this_inc);
        void rotate(const Vec3& new_rot_this_inc);

        [[nodiscard]] double kinetic_energy() const;
        [[nodiscard]] double translational_energy() const;
        [[nodiscard]] double rotational_energy() const;

        [[nodiscard]] std::size_t number_of_contacts() const;

        void reset_contacts();

        [[nodiscard]] virtual std::string get_output_string() const;

    protected:
        template<typename T>
        void sum_contact_forces(const T& contacts);

        template<typename ContactType, typename VectorType>
        void add_contact(ContactType contact, std::size_t index, int direction, VectorType& contacts);

        template<typename VectorType>
        void remove_contact(std::size_t index, VectorType& contacts);

    private:
        double radius_;   // Not const due to particle swelling
        double inertia_;
        /*
          Vector of all contacts, first is a pointer to the contact
          second is a multiplier (1, -1) to get_parameter the correct direction.
          Might be better done with another type than int
        */
    };
}

#include "spherical_particle_base.tpp"

#endif //DEMSIM_SPHERICAL_PARTICLE_BASE_H
