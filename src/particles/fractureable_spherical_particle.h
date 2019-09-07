//
// Created by erolsson on 2019-01-19.
//

#ifndef DEMSIM_FRACTUREABLE_SPHERICAL_PARTICLE_H
#define DEMSIM_FRACTUREABLE_SPHERICAL_PARTICLE_H

#include "spherical_particle_base.h"

#include "../engine/contact.h"

#include <random>

namespace DEM {
    class ParticleCrack {
    public:
        ParticleCrack(const Vec3& position, double force, std::size_t id_of_impacter, const Vec3& normal) :
        position_(position), force_(force), id_impacter_(id_of_impacter), normal_(normal){}
        [[nodiscard]] std::string get_output_string() const {
            std::stringstream ss;
            ss << position_.x() << ", " << position_.y() << ", " << position_.z() << ", " << force_;
            return ss.str();
        }
        [[nodiscard]] const Vec3& get_position() const { return position_; }
        [[nodiscard]] double get_force() const { return force_; }
        [[nodiscard]] std::size_t get_impacter_id() const { return  id_impacter_; }
        [[nodiscard]] const Vec3& get_normal() const { return normal_; }
        void set_force(double f) { force_ = f; }

    private:
        Vec3 position_;
        double force_;
        std::size_t id_impacter_;
        Vec3 normal_;
    };

    template<typename ForceModel>
    class FractureableSphericalParticle : public SphericalParticleBase<ForceModel> {
        using ContactType = Contact<ForceModel, FractureableSphericalParticle<ForceModel>>;
        using ContactPointerType = typename ContactMatrix<ContactType>::PointerType;
    public:
        FractureableSphericalParticle(double radius, const Vec3& position, const Vec3& velocity, MaterialBase* material,
                                      unsigned id);
        using SphericalParticleBase<ForceModel>::get_id;
        using SphericalParticleBase<ForceModel>::get_position;
        void sum_contact_forces() {
            SphericalParticleBase<ForceModel>::sum_contact_forces(contacts_);
        }

        void add_contact(ContactPointerType contact, std::size_t index, int direction) {
            SphericalParticleBase<ForceModel>::add_contact(contact, index, direction, contacts_);
        }

        void remove_contact(std::size_t index) {
            SphericalParticleBase<ForceModel>::remove_contact(index, contacts_);
        }

        void fracture_particle(const Vec3& position, double force, std::size_t id_of_impacter, const Vec3& normal);
        [[nodiscard]] const std::vector<ParticleCrack>& get_particle_cracks() const { return cracks_; }

    private:
        ContactVector<std::pair<ContactPointerType, int> >  contacts_;
        std::vector<ParticleCrack> cracks_;
        using SphericalParticleBase<ForceModel>::material_;
        [[nodiscard]] std::vector<ParticleCrack>::iterator has_crack_at_position(const Vec3& position);
    };


}

#include "fractureable_spherical_particle.tpp"

#endif //DEMSIM_FRACTUREABLE_SPHERICAL_PARTICLE_H
