//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_SPHERICAL_PARTICLE_H
#define DEMSIM_SPHERICAL_PARTICLE_H

#include <vector>
#include <memory>
#include <sstream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

#include "material_base.h"
#include "vec3.h"
#include "contact_vector.h"
#include "contact.h"
#include "particle_base.h"

namespace DEM {
    // Please don't change this one
    constexpr double pi = 3.1415;

    template<typename ForceModel, typename SphericalParticle> class Contact;

    template<typename ForceModel>
    class SphericalParticle : public ParticleBase<ForceModel> {
        using ContactPointerType = std::shared_ptr<Contact<ForceModel, SphericalParticle<ForceModel> >>;

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
        using ParticleBase<ForceModel>::fn_;
        using ParticleBase<ForceModel>::ft_;
        using ParticleBase<ForceModel>::torque_;
    public:

        // No assignment of particles and no plain copies
        SphericalParticle(const SphericalParticle&) = delete;
        SphericalParticle& operator=(const SphericalParticle&) = delete;
        SphericalParticle(double, const Vec3&, const Vec3&, MaterialBase*, unsigned );
        virtual ~SphericalParticle() = default;

        double get_radius() const { return radius_; }
        double get_inertia() const { return inertia_; }

        void move(const Vec3& new_disp_this_inc)
        {
            position_ += new_disp_this_inc;
            displacement_this_inc_ = new_disp_this_inc;
        }

        void rotate(const Vec3& new_rot)
        {
            rot_this_inc_ = new_rot;
            rot_ += new_rot;
        }

        double kinetic_energy() const
        {
            return mass_*dot_product(velocity_, velocity_)/2 + inertia_*dot_product(ang_vel_, ang_vel_)/2;
        }

        void reset_contacts() {
            number_of_contacts_ = 0;
            fn_ = 0*fn_;
            ft_ = 0*ft_;
            torque_ = 0*torque_;
        }

        void sum_contact_forces();
        std::vector<SphericalParticle*> get_neighbours() const;
        std::size_t number_of_contacts() const;

        void add_contact(ContactPointerType, std::size_t, int);
        void remove_contact(std::size_t);

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

    template<typename ForceModel>
    SphericalParticle<ForceModel>::SphericalParticle(double radius, const Vec3& pos, const Vec3& v, MaterialBase* mat,
                                                     unsigned id):
            ParticleBase<ForceModel>(4.*pi*radius*radius*radius/3*mat->density, pos, v, mat, id),
            radius_(radius),
            inertia_(2*mass_*radius_*radius_/5),
            contacts_(ContactVector<std::pair<ContactPointerType, int> >())
    {
        // Empty constructor
    }

    template<typename ForceModel>
    void
    SphericalParticle<ForceModel>::add_contact(SphericalParticle::ContactPointerType c, std::size_t pos,
                                               int direction)
    {
        contacts_.insert(pos, std::make_pair(c, direction));
    }

    template<typename ForceModel>
    void SphericalParticle<ForceModel>::remove_contact(std::size_t pos)
    {
        contacts_.erase(pos);
    }

    template<typename ForceModel>
    std::string SphericalParticle<ForceModel>::get_output_string() const
    {
        std::stringstream ss;
        ss << id_ << ", " << position_.x << ", " << position_.y << ", " << position_.z << ", ";
        ss << rot_.x << ", " << rot_.y << ", " << rot_.z << ", " << radius_<< ", ";
        ss << kinetic_energy() << ", " << material_->id;
        return ss.str();
    }


}


#endif //DEMSIM_SPHERICAL_PARTICLE_H
