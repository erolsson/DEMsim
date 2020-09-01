//
// Created by erolsson on 24/08/2020.
//

#ifndef DEMSIM_PERIODIC_BC_HANDLER_H
#define DEMSIM_PERIODIC_BC_HANDLER_H

#include <array>
#include <map>
#include <vector>

namespace DEM {
    struct Interval {
        double max = 0;
        double min = 0;
    };
    class Vec3;


    template<typename ForceModel, typename ParticleType>
    class CollisionDetector;

    template<typename ForceModel, typename ParticleType>
    class Contact;

    template<typename ForceModel, typename ParticleType>
    class Engine;

    template<typename ContactType>
    class ContactMatrix;


    template<typename ForceModel, typename ParticleType>
    class PeriodicBCHandler {

    public:
        using CollisionDetectorType = CollisionDetector<ForceModel, ParticleType>;
        using ContactType = Contact<ForceModel, ParticleType>;
        using EngineType = Engine<ForceModel, ParticleType>;
        PeriodicBCHandler(EngineType& engine,
                          std::vector<ParticleType*>& simulation_particles,
                          CollisionDetectorType& collision_detector,
                          ContactMatrix<ContactType>& contacts);
        void fulfill_periodic_bc();
        void add_periodic_bc(char axis, double boundary_min, double boundary_max);
        void set_periodic_bc_strain_rate(char axis, double strain_rate);

        void create_periodic_bc_contacts();
        void destroy_periodic_bc_contacts();

    private:

        constexpr static auto particle_id_sort = [](const auto& p1, const auto& rhs) {return p1->get_id() < rhs; };
        double stretch_ = 1e-2;
        std::map<std::size_t, std::array<ParticleType*, 7>> mirror_particles_;
        std::array<Interval, 3> boundaries_ {};
        std::array<double, 3> velocities_ {0., 0., 0.};
        std::array<bool, 3> active_directions_ = {false, false, false};
        EngineType& engine_;
        std::vector<ParticleType*>& simulation_particles_;
        CollisionDetector<ForceModel, ParticleType>& collision_detector_;
        ContactMatrix<ContactType>& contacts_;
        std::size_t no_active_directions_ = 0;

        void move_periodic_boundaries();
        void move_mirror_particles(std::size_t particle_idx);
        void create_corner_particles(std::size_t particle_idx);
        void create_mirror_particles(std::size_t particle_idx);
        void remove_mirror_particles(std::size_t particle_idx);
        void respect_boundaries(std::size_t particle_idx);

        ParticleType* get_mirror_particle(ParticleType* simulation_particle, std::size_t direction);
        ParticleType* create_mirror_particle(const ParticleType* simulation_particle, std::size_t direction,
                                             const Vec3& position);

        void remove_mirror_particle(ParticleType* simulation_particle, std::size_t direction);


        void position_corner_particle(ParticleType* simulation_particle,
                                      std::size_t direction,
                                      const std::array<bool, 3>& mirror_directions);

        bool has_mirror_particle(const ParticleType* particle) const;

    };
}

#include "periodic_bc_handler.tpp"

#endif //DEMSIM_PERIODIC_BC_HANDLER_H
