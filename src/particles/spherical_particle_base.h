//
// Created by erolsson on 16/05/19.
//

#ifndef DEMSIM_SPHERICAL_PARTICLE_BASE_H
#define DEMSIM_SPHERICAL_PARTICLE_BASE_H

#include "particle_base.h"

namespace DEM {
    class ParameterMap;
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

        SphericalParticleBase(double radius, const Vec3& position, const Vec3& velocity, const MaterialBase* material,
                              std::size_t object_id, std::size_t collision_id);
        SphericalParticleBase(const ParameterMap& parameters, MaterialBase* material);
        virtual ~SphericalParticleBase() = default;
        SphericalParticleBase& operator=(const SphericalParticleBase&) = delete;

        [[nodiscard]] std::size_t get_collision_id() const { return collision_id_; }
        void set_collision_id(std::size_t collision_id) { collision_id_ = collision_id;}
        [[nodiscard]] double get_radius() const { return radius_; }
        [[nodiscard]] double get_inertia() const { return inertia_; }

        void move(const Vec3& new_disp_this_inc);
        void set_position(const Vec3& new_position) {position_ = new_position; }
        void rotate(const Vec3& new_rot_this_inc);

        [[nodiscard]] double kinetic_energy() const;
        [[nodiscard]] double translational_energy() const;
        [[nodiscard]] double rotational_energy() const;

        [[nodiscard]] std::size_t number_of_contacts() const;

        void reset_contacts();

        [[nodiscard]] virtual std::string get_output_string() const;
        [[nodiscard]] virtual std::string restart_data() const;

    protected:
        template<typename T>
        void sum_contact_forces(const T& contacts);

        template<typename ContactType, typename VectorType>
        void add_contact(ContactType contact, std::size_t index, int direction, VectorType& contacts);

        template<typename VectorType>
        void remove_contact(std::size_t index, VectorType& contacts);

    private:
        std::size_t collision_id_;
        double radius_;   // Not const due to particle swelling
        double inertia_;

    };
}

#include "spherical_particle_base.tpp"

#endif //DEMSIM_SPHERICAL_PARTICLE_BASE_H
