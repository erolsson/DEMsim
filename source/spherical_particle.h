//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_SPHERICAL_PARTICLE_H
#define DEMSIM_SPHERICAL_PARTICLE_H

#include <sstream>
#include <vector>

#include "contact.h"
#include "contact_vector.h"
#include "material_base.h"
#include "particle_base.h"
#include "vec3.h"

namespace DEM {
    // Please don't change this one
    constexpr double pi = 3.1415;
    template<typename ForceModel, typename SphericalParticle> class Contact;

    template<typename ForceModel>
    class SphericalParticle : public ParticleBase<ForceModel> {
        using ContactType = Contact<ForceModel, SphericalParticle<ForceModel>>;
        using ContactPointerType = typename ContactMatrix<ContactType>::PointerType;

        using ParticleBase<ForceModel>::id_;
        using ParticleBase<ForceModel>::material_;
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

    public:

        // No assignment of particles and no plain copies

        SphericalParticle(double radius, const Vec3& position, const Vec3& velocity, MaterialBase* material,
                          unsigned id);

        SphericalParticle(const SphericalParticle&) = delete;
        SphericalParticle& operator=(const SphericalParticle&) = delete;

        double get_radius() const { return radius_; }
        double get_inertia() const { return inertia_; }

        void move(const Vec3& new_disp_this_inc);
        void rotate(const Vec3& new_rot_this_inc);

        double kinetic_energy() const;

        void sum_contact_forces();
        std::size_t number_of_contacts() const;

        void add_contact(ContactPointerType contact, std::size_t index, int direction);
        void remove_contact(std::size_t index);
        void reset_contacts();

        std::string get_output_string() const;
    private:
        double radius_;   // Not const due to particle swelling
        double inertia_;
        /*
          Vector of all contacts, first is a pointer to the contact
          second is a multiplier (1, -1) to get the correct direction.
          Might be better done with another type than int
        */
        ContactVector<std::pair<ContactPointerType, int> >  contacts_;
    };




}

#include "spherical_particle.tpp"

#endif //DEMSIM_SPHERICAL_PARTICLE_H
