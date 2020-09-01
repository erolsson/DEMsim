//
// Created by erolsson on 24/08/2020.
//

#include "periodic_bc_handler.h"

#include "collision_detection/collision_detector.h"
#include "engine.h"
#include "../utilities/contact_matrix.h"
#include "../utilities/vec3.h"

using namespace DEM;
template<typename ForceModel, typename ParticleType>
PeriodicBCHandler<ForceModel, ParticleType>::PeriodicBCHandler(EngineType& engine,
                                                               std::vector<ParticleType*>& simulation_particles,
                                                               CollisionDetectorType& collision_detector,
                                                               ContactMatrix<ContactType>& contacts) :
    engine_(engine),
    simulation_particles_(simulation_particles),
    collision_detector_(collision_detector),
    contacts_(contacts)
{

}


template<typename ForceModel, typename ParticleType>
void DEM::PeriodicBCHandler<ForceModel, ParticleType>::fulfill_periodic_bc()
{
    move_periodic_boundaries();
    for (std::size_t i = 0; i != simulation_particles_.size(); ++i) {
        move_mirror_particles(i);
        create_mirror_particles(i);
        create_corner_particles(i);
        respect_boundaries(i);
        if (mirror_particles_.count(simulation_particles_[i]->get_id()) > 0) {
            for (std::size_t j = 0; j != 7; ++j){
                const auto& mp = mirror_particles_[i][j];
                if (mp != nullptr) {
                    std::cout << "Mirror particle " << mp->get_id() << " at position " << mp->get_position()
                              << " with address " << mp << "\n";
                }
            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::move_periodic_boundaries() {
    for (unsigned i = 0; i != 3; ++i) {
        if (active_directions_[i]) {
            boundaries_[i].max += velocities_[i]*engine_.get_time_increment().count();
            boundaries_[i].min -= velocities_[i]*engine_.get_time_increment().count();
        }
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::move_mirror_particles(std::size_t particle_idx) {
    auto it = mirror_particles_.find(particle_idx);
    if (it != mirror_particles_.end()) {
        auto sim_particle = *std::lower_bound(simulation_particles_.begin(), simulation_particles_.end(), particle_idx,
                                              particle_id_sort);
        for (const auto& mp: it->second) {
            if (mp != nullptr) {
                mp->move(sim_particle->get_displacement_this_increment());
            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::create_mirror_particles(std::size_t particle_idx) {
    auto particle = simulation_particles_[particle_idx];
    for (unsigned direction = 0; direction != active_directions_.size(); ++direction) {
        if (active_directions_[direction] && (mirror_particles_.count(particle_idx) == 0
            || mirror_particles_[particle_idx][direction] == nullptr)) {
            auto d1 = particle->get_position()[direction] - particle->get_radius() - stretch_
                      - boundaries_[direction].min;
            auto d2 = boundaries_[direction].max - particle->get_radius() - stretch_
                       - particle->get_position()[direction];
            if (d1 < 0 || d2 < 0) {
                Vec3 position = particle->get_position();
                position[direction] -= (boundaries_[direction].max
                                        - boundaries_[direction].min)*d1/(std::abs(d1));

                auto p = create_mirror_particle(particle, direction, position);
                std::cout << "Creating mirror particle at position: " << position
                          << " with address " << p << std::endl;

            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::remove_mirror_particles(std::size_t particle_idx) {
    auto particle = simulation_particles_[particle_idx];
    if (mirror_particles_.count(particle->get_id()) > 0) {
        bool inside = true;
        for (std::size_t direction = 0; direction != 3; ++direction) {
            auto d1 = particle->get_position()[direction] - particle->get_radius() - stretch_
                      - boundaries_[direction].min;
            auto d2 = boundaries_[direction].max - particle->get_radius() - stretch_
                      - particle->get_position()[direction];
            inside = inside && (!active_directions_[direction] || (d1 > 0 && d2 > 0));
        }
        if (inside) {
            for (std::size_t idx = 0; idx != 7; ++idx) {
                remove_mirror_particle(particle, idx);
                std::cout << "Mirror particle on idx " << idx << " removed\n";
            }
            mirror_particles_.erase(particle->get_id());
        }
    }
}


template<typename ForceModel, typename ParticleType>
void DEM::PeriodicBCHandler<ForceModel, ParticleType>::respect_boundaries(std::size_t particle_idx) {
    std::array<double, 3> d1;
    std::array<double, 3> d2;
    std::array<bool, 3> directions = {false, false, false};
    ParticleType* particle = nullptr;
    for (std::size_t direction = 0; direction != active_directions_.size(); ++direction) {
        if (active_directions_[direction]) {
            particle = simulation_particles_[particle_idx];
            d1[direction] = particle->get_position()[direction] - boundaries_[direction].min;
            d2[direction] = boundaries_[direction].max - particle->get_position()[direction];
            if (d1[direction] < 0 || d2[direction] < 0) {
                directions[direction] = true;
            }
        }
    }

    std::size_t overlapping_directions = std::count(directions.begin(), directions.end(), true);
    if (overlapping_directions > 0) {
        std::size_t mirror_idx;
        Vec3 distance_to_move = Vec3(0., 0., 0.);
        std::cout << "overlapping_directions: " << overlapping_directions << "\n";
        if (overlapping_directions == 1) {
            mirror_idx = std::find(directions.begin(), directions.end(), true) - directions.begin();
        }
        else if (directions[0] == true && directions[1] == true) {
            mirror_idx = 3;
        }
        else if (directions[0] == true && directions[2] == true) {
            mirror_idx = 4;
        }
        else if (directions[1] == true && directions[2] == true) {
            mirror_idx = 5;
        }
        else {
            mirror_idx = 6;
        }

        for (std::size_t direction = 0; direction != active_directions_.size(); ++direction) {
            if (directions[direction]) {
                distance_to_move[direction] = -(boundaries_[direction].max
                                              - boundaries_[direction].min)*d1[direction]/(std::abs(d1[direction]));
            }
        }
        auto mirror_particle = get_mirror_particle(particle, mirror_idx);
        mirror_particles_[particle->get_id()][mirror_idx] = particle;
        simulation_particles_[particle_idx] = mirror_particle;
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
    ++no_active_directions_;
    boundaries_[axis_idx].min = boundary_min;
    boundaries_[axis_idx].max = boundary_max;
}


template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::set_periodic_bc_strain_rate(char axis, double strain_rate) {
    std::string directions = "xyz";
    auto axis_idx = directions.find(axis);
    if (axis_idx == std::string::npos) {
        throw std::invalid_argument("axis argument must be x, y or z");
    }
    velocities_[axis_idx] = strain_rate*(boundaries_[axis_idx].max - boundaries_[axis_idx].min);
}


template<typename ForceModel, typename ParticleType>
ParticleType* PeriodicBCHandler<ForceModel, ParticleType>::get_mirror_particle(ParticleType* simulation_particle,
                                                                               std::size_t direction) {
    return mirror_particles_[simulation_particle->get_id()][direction];
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::remove_mirror_particle(ParticleType* simulation_particle,
                                                                         std::size_t direction) {
    auto p = mirror_particles_[simulation_particle->get_id()][direction];
    if (p != nullptr) {
        collision_detector_.remove_particle(p);
        delete p;
        mirror_particles_[simulation_particle->get_id()][direction] = nullptr;
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::create_corner_particles(std::size_t particle_idx) {
    auto particle = simulation_particles_[particle_idx];
    position_corner_particle(particle, 3, {true, true, false});
    position_corner_particle(particle, 4, {true, false, true});
    position_corner_particle(particle, 5, {false, true, true});
    position_corner_particle(particle, 6, {true, true, true});
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::position_corner_particle(ParticleType* simulation_particle,
                                                                           std::size_t direction,
                                                                           const std::array<bool, 3>& mirror_directions) {
    auto p_id = simulation_particle->get_id();
    bool create = mirror_particles_.count(p_id) > 0 && mirror_particles_[p_id][direction] == nullptr;
    for (unsigned axis = 0; axis != mirror_directions.size(); ++axis) {
        if (create && mirror_directions[axis] ) {
            create = create && (mirror_particles_.count(p_id) > 0) &&  (mirror_particles_[p_id][axis] != nullptr);
        }
    }

    if (create) {
        Vec3 position = Vec3(0, 0, 0);
        for (unsigned axis = 0; axis != mirror_directions.size(); ++axis) {
            position[axis] = !mirror_directions[axis]*simulation_particles_[p_id]->get_position()[axis];
            auto axis_mirior_p = mirror_particles_[p_id][axis];
            if (axis_mirior_p != nullptr) {
                position[axis] += mirror_directions[axis]*axis_mirior_p->get_position()[axis];
            }
        }
        auto p = create_mirror_particle(simulation_particle, direction, position);
        std::cout << "Creating corner particle at: " << position << " with address " << p << "\n";
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::create_periodic_bc_contacts() {
    for (const auto& collision_pair: collision_detector_.contacts_to_create()){
        auto p1 = collision_pair.particle1;
        auto p2 = collision_pair.particle2;

        bool mirror_p1 = has_mirror_particle(p1);
        bool mirror_p2 = has_mirror_particle(p2);

        if (mirror_p1 || mirror_p2) {
            auto id1 = collision_pair.get_id_pair().first;
            auto id2 = collision_pair.get_id_pair().second;
            std::cout << "Creating contact between " << id1 << " and " << id2 << "\n";
            Contact<ForceModel, ParticleType>* c = nullptr;
            auto sim_particle_1 = *std::lower_bound(simulation_particles_.begin(), simulation_particles_.end(), id1,
                                                    particle_id_sort);
            if ( p2 != nullptr) {
                c = contacts_.create_item_inplace(id1, id2, p1, p2, engine_.get_time_increment());
                auto sim_particle_2 = *std::lower_bound(simulation_particles_.begin(), simulation_particles_.end(), id2,
                                                        particle_id_sort);
                sim_particle_2->add_contact(c, id1, -1);
            }
            else {
                auto s = collision_pair.surface;
                c = contacts_.create_item_inplace(id1, id2, p1, s, engine_.get_time_increment());
                s->add_contact(c, id1);
            }
            sim_particle_1->add_contact(c, id2, 1);
            if (mirror_p1) {
                for (const auto& mp: mirror_particles_[id1]) {
                    if (mp != nullptr) {
                        mp->add_contact(c, id2, 1);
                    }
                }
            }

            if (mirror_p2) {
                for (const auto& mp: mirror_particles_[id2]) {
                    if (mp != nullptr) {
                        mp->add_contact(c, id1, -1);
                    }
                }
            }
        }
    }

}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::destroy_periodic_bc_contacts() {
    std::cout << "Removing " << collision_detector_.contacts_to_destroy().size() << " contacts\n";
    for (const auto& collision_pair: collision_detector_.contacts_to_destroy()){
        auto p1 = collision_pair.particle1;
        auto p2 = collision_pair.particle2;

        bool mirror_p1 = has_mirror_particle(p1);
        bool mirror_p2 = has_mirror_particle(p2);

        if (mirror_p1 || mirror_p2) {

            auto id1 = collision_pair.get_id_pair().first;
            auto id2 = collision_pair.get_id_pair().second;
            std::cout << "Removing: " << id1 << "  " << id2 << "\n";
            auto sim_particle_1 = *std::lower_bound(simulation_particles_.begin(), simulation_particles_.end(), id1,
                                                    particle_id_sort);
            contacts_.erase(id1, id2);
            if (p2 != nullptr) {
                auto sim_particle_2 = *std::lower_bound(simulation_particles_.begin(), simulation_particles_.end(), id2,
                                                        particle_id_sort);
                sim_particle_2->remove_contact(id1);
            }
            else {
                auto s = collision_pair.surface;
                s->remove_contact(id1);
            }
            sim_particle_1->remove_contact(id2);
            if (mirror_p1) {
                for (const auto& mp: mirror_particles_[id1]) {
                    if (mp != nullptr) {
                        mp->remove_contact(id2);
                    }
                }
            }

            if (mirror_p2) {
                for (const auto& mp: mirror_particles_[id2]) {
                    if (mp != nullptr) {
                        mp->remove_contact(id1);
                    }
                }
            }

        }

    }
    for (std::size_t idx = 0; idx != simulation_particles_.size(); ++idx) {
        remove_mirror_particles(idx);
    }
}

template<typename ForceModel, typename ParticleType>
bool PeriodicBCHandler<ForceModel, ParticleType>::has_mirror_particle(const ParticleType* particle) const {
    return mirror_particles_.find(particle->get_id()) != mirror_particles_.end();
}

template<typename ForceModel, typename ParticleType>
ParticleType* PeriodicBCHandler<ForceModel, ParticleType>::create_mirror_particle(
        const ParticleType* simulation_particle, std::size_t direction, const Vec3& position) {
    auto p = new ParticleType(simulation_particle->get_radius(), position, simulation_particle->get_velocity(),
                              simulation_particle->get_material(), simulation_particle->get_id());
    collision_detector_.add_particle(p);
    mirror_particles_[simulation_particle->get_id()][direction] = p;
    return p;
}
