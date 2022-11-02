//
// Created by erolsson on 24/08/2020.
//

#ifndef DEMSIM_PERIODIC_BC_HANDLER_H
#define DEMSIM_PERIODIC_BC_HANDLER_H

#include <array>
#include <map>
#include <set>
#include <vector>

namespace DEM {
    struct Interval {
        double max = 0;
        double min = 0;
    };
    class Vec3;
    class ParameterMap;

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

        PeriodicBCHandler(EngineType& engine,
                          std::vector<ParticleType*>& simulation_particles,
                          CollisionDetectorType& collision_detector,
                          ContactMatrix<ContactType>& contacts,
                          const std::vector<ParameterMap>& restart_data);
        void fulfill_periodic_bc();
        void add_periodic_bc(char axis, double boundary_min, double boundary_max);
        void set_periodic_bc_strain_rate(char axis, double strain_rate);
        void set_periodic_bc_velocity(char axis, double velocity);
        void set_boundary_stretch(double stretch) { stretch_ = stretch; }
        void create_periodic_bc_contacts();
        void destroy_periodic_bc_contacts();
        std::string print_periodic_bc() const;
        std::vector<std::string> mirror_particles_output() const;
        std::vector<std::string> restart_data() const;
        [[nodiscard]] std::array<Interval, 3> get_periodic_boundaries() const { return boundaries_; }

    private:

        double stretch_ = 1e-9;
        std::map<std::size_t, std::array<ParticleType*, 7>> mirror_particles_;
        std::array<Interval, 3> boundaries_ {};
        std::array<double, 3> velocities_ {0., 0., 0.};
        std::array<double, 3> strain_rates_ {0., 0., 0.};
        std::array<bool, 3> active_directions_ {false, false, false};
        EngineType& engine_;
        std::vector<ParticleType*>& simulation_particles_;
        CollisionDetector<ForceModel, ParticleType>& collision_detector_;
        ContactMatrix<ContactType>& contacts_;
        std::size_t no_active_directions_ = 0;
        std::set<std::size_t> jump_particles_;

        void move_periodic_boundaries();
        void move_mirror_particles(ParticleType* simulation_particle);
        void create_corner_particles(ParticleType* simulation_particle);
        void create_mirror_particles(ParticleType* simulation_particle);
        void remove_mirror_particles(ParticleType* particle);
        void respect_boundaries(ParticleType* simulation_particle);

        ParticleType* get_mirror_particle(ParticleType* simulation_particle, std::size_t direction);
        ParticleType* get_simulation_particle(std::size_t particle_id);
        ParticleType* create_mirror_particle(const ParticleType* simulation_particle, std::size_t direction,
                                             const Vec3& position);

        ParticleType* create_mirror_particle(const ParticleType* simulation_particle, std::size_t direction,
                                             const Vec3& position, std::size_t collision_id);

        void remove_mirror_particle(ParticleType* simulation_particle, std::size_t direction);


        void position_corner_particle(ParticleType* simulation_particle,
                                      std::size_t direction,
                                      const std::array<bool, 3>& mirror_directions);

        bool has_mirror_particle(const ParticleType* particle) const;
        bool is_mirror_particle(const ParticleType* particle) const;

        void handle_mirror_particles_add_contact(Contact<ForceModel, ParticleType>* contact,
                                                 std::size_t id1, std::size_t id2, ParticleType* particle,
                                                 int contact_direction);

        static std::size_t direction_idx(char axis);
    };
}

#include "periodic_bc_handler.tpp"

#endif //DEMSIM_PERIODIC_BC_HANDLER_H
