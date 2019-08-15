//
// Created by erolsson on 2019-01-19.
//

#ifndef DEMSIM_FRACTUREABLE_SPHERICAL_PARTICLE_H
#define DEMSIM_FRACTUREABLE_SPHERICAL_PARTICLE_H

#include "spherical_particle_base.h"

#include "../engine/contact.h"

namespace DEM {

    template<typename ForceModel>
    class FractureableSphericalParticle : public SphericalParticleBase<ForceModel> {
        using ContactType = Contact<ForceModel, FractureableSphericalParticle<ForceModel>>;
        using ContactPointerType = typename ContactMatrix<ContactType>::PointerType;
    public:
        FractureableSphericalParticle(double radius, const Vec3& position, const Vec3& velocity, MaterialBase* material,
                                      unsigned id);
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

#include "fractureable_spherical_particle.tpp"

#endif //DEMSIM_FRACTUREABLE_SPHERICAL_PARTICLE_H
