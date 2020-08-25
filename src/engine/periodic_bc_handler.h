//
// Created by erolsson on 24/08/2020.
//

#ifndef DEMSIM_PERIODIC_BC_HANDLER_H
#define DEMSIM_PERIODIC_BC_HANDLER_H

#include <array>
#include <vector>

namespace DEM {
    struct Interval {
        double max = 0;
        double min = 0;
    };

    template<typename ForceModel, typename ParticleType>
    class CollisionDetector;
    template<typename ForceModel, typename ParticleType>
    class PeriodicBCHandler {

    public:
        using CollisionDetectorType = CollisionDetector<ForceModel, ParticleType>;
        PeriodicBCHandler(std::vector<ParticleType*>& simulation_particles,
                          CollisionDetectorType & collision_detector);
        void add_particle(ParticleType* p);
        void fulfill_periodic_bc();
        void add_periodic_bc(char axis, double boundary_min, double boundary_max);
    private:
        double stretch_ = 0;
        std::vector<std::array<ParticleType*, 7>> mirror_particles_ {};
        std::array<Interval, 3> boundaries_ {};
        std::array<bool, 3> active_directions_ = {false, false, false};
        std::vector<ParticleType*>& simulation_particles_;
        CollisionDetector<ForceModel, ParticleType>& collision_detector;

        void respect_boundaries(ParticleType* particle);
        void handle_mirror_particles(std::size_t particle_idx);

        ParticleType* get_mirror_particle(std::size_t idx, std::size_t direction);
        void remove_mirror_particle(std::size_t idx, std::size_t direction);

        void handle_corner_particles(std::size_t particle_idx);
        void position_corner_particle(std::size_t particle_idx,
                                      std::size_t direction,
                                      const std::array<bool, 3>& mirror_directions);
    };
}

#include "periodic_bc_handler.tpp"

#endif //DEMSIM_PERIODIC_BC_HANDLER_H
