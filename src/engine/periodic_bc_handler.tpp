//
// Created by erolsson on 24/08/2020.
//

#include "periodic_bc_handler.h"

#include "../utilities/vec3.h"

using namespace DEM;
template<typename ForceModel, typename ParticleType>
PeriodicBCHandler<ForceModel, ParticleType>::PeriodicBCHandler(std::vector<ParticleType*>& simulation_particles,
                                                                CollisionDetectorType& collision_detector) :
    simulation_particles_(simulation_particles),
    collision_detector(collision_detector)
{
    for(std::size_t i = 0; i != simulation_particles_.size(); ++i) {
        mirror_particles_.push_back({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr});
    }
}


template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::add_particle(ParticleType* p) {
    mirror_particles_.push_back({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr});
}


template<typename ForceModel, typename ParticleType>
void DEM::PeriodicBCHandler<ForceModel, ParticleType>::fulfill_periodic_bc()
{
    for (std::size_t i = 0; i != simulation_particles_.size(); ++i) {
        respect_boundaries(simulation_particles_[i]);
        handle_mirror_particles(i);
        handle_corner_particles(i);
    }
}


template<typename ForceModel, typename ParticleType>
void DEM::PeriodicBCHandler<ForceModel, ParticleType>::respect_boundaries(ParticleType* particle) {
    for (unsigned direction = 0; direction != active_directions_.size(); ++direction) {
        if (active_directions_[direction]) {
            double d1 = particle->get_position()[direction] - boundaries_[direction].min;
            double d2 = boundaries_[direction].max - particle->get_position()[direction];
            if (d1 < 0 || d2 < 0) {
                double d = -(boundaries_[direction].max - boundaries_[direction].min)*d1/std::abs(d1);
                Vec3 distance_to_move = Vec3(0., 0., 0.);
                distance_to_move[direction] = d;
                particle->move(distance_to_move);
            }
        }
    }
}


template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::add_periodic_bc(char axis, double boundary_min, double boundary_max) {
    std::string directions = "xyz";
    auto axis_idx = directions.find(axis);
    if (axis_idx == std::string::npos) {
        throw std::invalid_argument("axis argument must be x, y or z");
    }
    active_directions_[axis_idx] = true;
    boundaries_[axis_idx].min = boundary_min;
    boundaries_[axis_idx].max = boundary_max;
}


template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::handle_mirror_particles(std::size_t particle_idx) {
    for (unsigned direction = 0; direction != active_directions_.size(); ++direction) {
        if (active_directions_[direction]) {
            auto particle = simulation_particles_[particle_idx];
            double d1 = particle->get_position()[direction] - particle->get_radius() - stretch_
                        - boundaries_[direction].min;
            double d2 = boundaries_[direction].max - particle->get_radius() - stretch_
                        - particle->get_position()[direction];
            if (d1 < 0 || d2 < 0) {
                auto mirror_particle = get_mirror_particle(particle_idx, direction);
                Vec3 distance_to_move = Vec3(0, 0, 0);
                distance_to_move[direction] = -(boundaries_[direction].max
                                                - boundaries_[direction].min)*d1/(std::abs(d1));
                mirror_particle->set_position(particle->get_position() + distance_to_move);
            }
            else if (d1 > 0 && d2 > 0) {
                remove_mirror_particle(particle_idx, direction);
            }
        }

    }
}

template<typename ForceModel, typename ParticleType>
ParticleType* PeriodicBCHandler<ForceModel, ParticleType>::get_mirror_particle(std::size_t idx, std::size_t direction) {
    if (mirror_particles_[idx][direction] == nullptr) {
        mirror_particles_[idx][direction] = new ParticleType(*simulation_particles_[idx]);
    }
    return mirror_particles_[idx][direction];
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::remove_mirror_particle(std::size_t idx, std::size_t direction) {
    if (mirror_particles_[idx][direction] != nullptr) {
        delete mirror_particles_[idx][direction];
        mirror_particles_[idx][direction] = nullptr;
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::handle_corner_particles(std::size_t particle_idx) {
    position_corner_particle(particle_idx, 3, {true, true, false});
    position_corner_particle(particle_idx, 4, {true, false, true});
    position_corner_particle(particle_idx, 5, {false, true, true});
    position_corner_particle(particle_idx, 6, {true, true, true});
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::position_corner_particle(std::size_t particle_idx,
                                                                           std::size_t direction,
                                                                           const std::array<bool, 3>& mirror_directions) {
    bool create = true;
    for (unsigned axis = 0; axis != mirror_directions.size(); ++axis) {
        if (create && mirror_directions[axis]) {
            create = create && (mirror_particles_[particle_idx][axis] != nullptr);
        }
    }
    if (create) {
        auto mirror_p = get_mirror_particle(particle_idx, direction);
        Vec3 mirror_pos = Vec3(0, 0, 0);
        for (unsigned axis = 0; axis != mirror_directions.size(); ++axis) {
            mirror_pos[axis] = !mirror_directions[axis]*simulation_particles_[particle_idx]->get_position()[axis];
            auto axis_mirior_p = mirror_particles_[particle_idx][axis];
            if (axis_mirior_p != nullptr) {
                mirror_pos[axis] += mirror_directions[axis]*axis_mirior_p->get_position()[axis];
            }
            mirror_p->set_position(mirror_pos);
        }
    }
    else {
        remove_mirror_particle(particle_idx, direction);
    }
}


