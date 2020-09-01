//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_SPHERICAL_PARTICLE_H
#define DEMSIM_SPHERICAL_PARTICLE_H

#include <sstream>
#include <vector>

#include "../engine/contact.h"
#include "../utilities/contact_vector.h"
#include "../materials/material_base.h"
#include "particle_base.h"
#include "../utilities/vec3.h"
#include "spherical_particle_base.h"

namespace DEM {
    class ParameterMap;
    // Please don't change this one
    constexpr double pi = 3.1415;
    template<typename ForceModel, typename SphericalParticleBase> class Contact;

    template<typename ForceModel>
    class SphericalParticle : public SphericalParticleBase<ForceModel> {
        using ContactType = Contact<ForceModel, SphericalParticle<ForceModel>>;
        using ContactPointerType = typename ContactMatrix<ContactType>::PointerType;

    public:
        // No assignment of particles and no plain copies

        SphericalParticle(double radius, const Vec3& position, const Vec3& velocity, const MaterialBase* material,
                          unsigned id);
        SphericalParticle(const ParameterMap& parameters, MaterialBase* material);
        virtual ~SphericalParticle() = default;
        void sum_contact_forces() {
            SphericalParticleBase<ForceModel>::sum_contact_forces(contacts_);
    }

        void add_contact(ContactPointerType contact, std::size_t index, int direction) {
            SphericalParticleBase<ForceModel>::add_contact(contact, index, direction, contacts_);
        }

        void remove_contact(std::size_t index) {
            SphericalParticleBase<ForceModel>::remove_contact(index, contacts_);
        }

    private:
        ContactVector<std::pair<ContactPointerType, int> >  contacts_;
    };

}

#include "spherical_particle.tpp"

#endif //DEMSIM_SPHERICAL_PARTICLE_H
