//
// Created by erolsson on 24/08/2020.
//

#include "periodic_bc_handler.h"

#include <array>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "collision_detection/collision_detector.h"
#include "engine.h"
#include "../utilities/contact_matrix.h"
#include "../utilities/file_reading_functions.h"
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
PeriodicBCHandler<ForceModel, ParticleType>::PeriodicBCHandler(PeriodicBCHandler::EngineType& engine,
                                                               std::vector<ParticleType*>& simulation_particles,
                                                               PeriodicBCHandler::CollisionDetectorType& collision_detector,
                                                               ContactMatrix<ContactType>& contacts,
                                                               const std::vector<ParameterMap>& restart_data) :
        engine_(engine),
        simulation_particles_(simulation_particles),
        collision_detector_(collision_detector),
        contacts_(contacts)
{
    for (const auto& restart_line: restart_data) {
        if (restart_line.get_parameter<std::string>("data") == "stretch") {
            stretch_ = restart_line.get_parameter("stretch");
        }

        std::string axes = "xyz";
        if (restart_line.get_parameter<std::string>("data") == "boundary") {
            for (std::size_t dir = 0; dir != 3; ++dir) {
                boundaries_[dir].min = restart_line.get_parameter<>(axes.substr(dir, 1) + "min");
                boundaries_[dir].max = restart_line.get_parameter<>(axes.substr(dir, 1) + "max");
                if (boundaries_[dir].max - boundaries_[dir].min > 0) {
                    active_directions_[dir] = true;
                    ++no_active_directions_;
                }
            }
        }
        if (restart_line.get_parameter<std::string>("data") == "velocity") {
            for (std::size_t dir = 0; dir != 3; ++dir) {
                velocities_[dir] = restart_line.get_parameter<>("v" + axes.substr(dir, 1));
            }
        }

        if (restart_line.get_parameter<std::string>("data") == "strain_rate") {
            for (std::size_t dir = 0; dir != 3; ++dir) {
                strain_rates_[dir] = restart_line.get_parameter<>("e" + axes.substr(dir, 1));
            }
        }
        if (restart_line.get_parameter<std::string>("data") == "mirror_particle") {
            auto id = restart_line.get_parameter<std::size_t>("id");
            auto axis = restart_line.get_parameter<std::size_t>("axis");
            auto position = restart_line.get_vec3("pos");
            auto collision_id = restart_line.get_parameter<std::size_t>("collision_id");
            auto sim_particle = get_simulation_particle(id);

            create_mirror_particle(sim_particle, axis, position, collision_id);
        }

        if (restart_line.get_parameter<std::string>("data") == "mirror_contact") {
            auto contact_type = restart_line.get_parameter<std::string>("pair");
            std::regex contact_type_re ("(.+)-(.+)");
            std::smatch sm_type;
            std::regex_match(contact_type, sm_type, contact_type_re);

            std::regex object_re ("(\\D+)(\\d)?");
            std::array<ParticleType*, 2> contact_particles = {nullptr, nullptr};
            std::array<std::size_t, 2> object_id;
            Surface<ForceModel, ParticleType>* surface = nullptr;
            for (unsigned i = 0; i != 2; ++ i) {
                std::smatch sm_obj;
                object_id[i] = restart_line.get_parameter<std::size_t>("object" + std::to_string(i+1));
                std::string obj_string = sm_type[i+1];
                std::regex_match(obj_string, sm_obj, object_re);
                std::string contact_object = sm_obj[1];
                if (contact_object == "mirror") {
                    std::size_t mirror_axis = std::stoi(sm_obj[2]);
                    contact_particles[i] = mirror_particles_[object_id[i]][mirror_axis];
                }
                else if (contact_object == "particle") {
                    contact_particles[i] = get_simulation_particle(object_id[i]);
                }
                else if (contact_object == "surface") {
                    surface = engine_.get_surface(object_id[i]);
                }
            }
            // Remove the corresponding contact from the contact matrix and replace it with the mirror particle
            // contact
            auto contact_to_remove = contacts_.get(object_id[0], object_id[1]);
            auto contact_to_remove_particles = contact_to_remove->get_particles();
            auto contact_to_remove_surface = contact_to_remove->get_surface();
            contact_to_remove_particles.first->remove_contact(object_id[1]);

            if (contact_to_remove_surface == nullptr) {
                contact_to_remove_particles.second->remove_contact(object_id[0]);
            }
            else {
                surface->remove_contact(object_id[0]);
            }

            auto remove_contact_from_mirror_particles = [this](auto particle, auto id2) mutable {
                if (has_mirror_particle(particle)) {
                    for (auto& mp: mirror_particles_[particle->get_id()]) {
                        if (mp != nullptr) {
                            mp->remove_contact(id2);
                        }
                    }
                }
            };

            remove_contact_from_mirror_particles(contact_to_remove_particles.first, object_id[1]);
            remove_contact_from_mirror_particles(contact_to_remove_particles.second, object_id[0]);
            contacts_.erase(object_id[0], object_id[1]);

            typename ContactMatrix<Contact<ForceModel, ParticleType>>::PointerType c = nullptr;
            if (surface == nullptr) {
                c = contacts_.create_item_inplace(object_id[0], object_id[1], contact_particles[0],
                                                  contact_particles[1], engine_.get_time_increment(), restart_line);
                contact_particles[1]->add_contact(c, object_id[0], -1);
                if (is_mirror_particle(contact_particles[1])) {
                    get_simulation_particle(object_id[1])->add_contact(c, object_id[0], -1);
                }
            }
            else {
                c = contacts_.create_item_inplace(object_id[0], object_id[1], contact_particles[0], surface,
                                                  engine_.get_time_increment(), restart_line);
                surface->add_contact(c, object_id[0]);
            }
            contact_particles[0]->add_contact(c, object_id[1], 1);
            if (is_mirror_particle(contact_particles[0])) {
                get_simulation_particle(object_id[0])->add_contact(c, object_id[1], 1);
            }
        }
    }

    // Finally iterate over all contacts to connect all mirror particles
    for (const auto& contact: contacts_.get_objects()) {
        auto contact_particles = contact->get_particles();
        auto id_pair = contact->get_id_pair();
        handle_mirror_particles_add_contact(contact, id_pair.first, id_pair.second,
                                            contact_particles.first, 1);
        if (contact_particles.second != nullptr) {
            handle_mirror_particles_add_contact(contact, id_pair.second, id_pair.first,
                                                contact_particles.second, -1);
        }
    }
}



template<typename ForceModel, typename ParticleType>
void DEM::PeriodicBCHandler<ForceModel, ParticleType>::fulfill_periodic_bc()
{
    jump_particles_.clear();
    move_periodic_boundaries();
    for (auto& p: simulation_particles_) {
        move_mirror_particles(p);
        respect_boundaries(p);
        create_mirror_particles(p);
        create_corner_particles(p);
        // remove_mirror_particles(p);
        handle_jump_contacts();
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::move_periodic_boundaries() {
    for (unsigned i = 0; i != 3; ++i) {
        if (active_directions_[i]) {
            double v = velocities_[i];
            if (strain_rates_[i] != 0.) {
                v = strain_rates_[i]*(boundaries_[i].max - boundaries_[i].min)/2;
                velocities_[i] = v;
            }
            boundaries_[i].max += v*engine_.get_time_increment().count();
            boundaries_[i].min -= v*engine_.get_time_increment().count();
        }
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::move_mirror_particles(ParticleType* simulation_particle) {
    auto calc_bc_velocity = [vel=this->velocities_, dt=engine_.get_time_increment().count()](const auto& axis,
                                                                                             const auto& d) {
        Vec3 boundary_vel = Vec3(0, 0, 0);
        if (d[axis] < 0) {
            boundary_vel[axis] = -2*vel[axis]*dt;
        }
        else {
            boundary_vel[axis] = 2*vel[axis]*dt;
        }
        return boundary_vel;
    };

    auto it = mirror_particles_.find(simulation_particle->get_id());
    if (it != mirror_particles_.end()) {
        for (std::size_t i = 0; i != 7; ++i) {
            auto& mp = (*it).second[i];
            if (mp != nullptr) {
                Vec3 d = mp->get_position() - simulation_particle->get_position();
                Vec3 bc_velocity = Vec3(0, 0, 0);
                if(i < 3) {
                    bc_velocity = calc_bc_velocity(i, d);
                }
                else if (i == 3) {
                    bc_velocity = calc_bc_velocity(0, d) + calc_bc_velocity(1, d);
                }
                else if (i == 4) {
                    bc_velocity = calc_bc_velocity(0, d) + calc_bc_velocity(2, d);
                }
                else if (i == 5) {
                    bc_velocity = calc_bc_velocity(1, d) + calc_bc_velocity(2, d);
                }
                else if (i == 6) {
                    bc_velocity = calc_bc_velocity(0, d) + calc_bc_velocity(1, d) + calc_bc_velocity(2, d);
                }

                mp->move(bc_velocity);
                mp->move(simulation_particle->get_displacement_this_increment());
                mp->set_velocity(simulation_particle->get_velocity());

                mp->rotate(simulation_particle->get_rotation_this_increment());
                mp->set_angular_velocity(simulation_particle->get_angular_velocity());

                // Moves the particle to the other side of the periodic BC if it has moved to far away
                for (unsigned j = 0; j != 3; ++j) {
                    if (abs(mp->get_position()[j]) > 2*boundaries_[j].max && active_directions_[j]) {
                        Vec3 new_pos = mp->get_position();
                        new_pos[j] -= 4*boundaries_[j].max*mp->get_position()[j]/abs(mp->get_position()[j]);
                        mp->set_position(new_pos);
                        jump_particles_.push_back(mp);
                    }
                }
            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::create_mirror_particles(ParticleType* simulation_particle) {
    for (unsigned direction = 0; direction != active_directions_.size(); ++direction) {
        if (active_directions_[direction] && (mirror_particles_.count(simulation_particle->get_id()) == 0
                                              || mirror_particles_[simulation_particle->get_id()][direction] == nullptr)) {
            const auto d1 = simulation_particle->get_position()[direction] - simulation_particle->get_radius()
                            - stretch_ - boundaries_[direction].min;
            const auto d2 = boundaries_[direction].max - simulation_particle->get_radius() - stretch_
                            - simulation_particle->get_position()[direction];
            if (d1 < 0 || d2 < 0) {
                Vec3 position = simulation_particle->get_position();
                position[direction] -= (boundaries_[direction].max
                                        - boundaries_[direction].min)*d1/(std::abs(d1));

                create_mirror_particle(simulation_particle, direction, position);
            }
        }
    }

}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::remove_mirror_particles(ParticleType* particle) {
    if (mirror_particles_.count(particle->get_id()) > 0) {
        for (std::size_t direction = 0; direction != 3; ++direction) {
            const auto d1 = particle->get_position()[direction] - particle->get_radius() - stretch_
                            - boundaries_[direction].min;
            const auto d2 = boundaries_[direction].max - (particle->get_radius() + stretch_
                                                          + particle->get_position()[direction]);
            if (d1 > 0 && d2 > 0) {
                switch (direction) {
                    case 0:
                        remove_mirror_particle(particle, 0);
                        remove_mirror_particle(particle, 3);
                        remove_mirror_particle(particle, 4);
                        remove_mirror_particle(particle, 6);
                        break;
                    case 1:
                        remove_mirror_particle(particle, 1);
                        remove_mirror_particle(particle, 3);
                        remove_mirror_particle(particle, 5);
                        remove_mirror_particle(particle, 6);
                        break;
                    case 2:
                        remove_mirror_particle(particle, 2);
                        remove_mirror_particle(particle, 4);
                        remove_mirror_particle(particle, 5);
                        remove_mirror_particle(particle, 6);
                        break;
                }

            }
        }

    }
}


template<typename ForceModel, typename ParticleType>
void DEM::PeriodicBCHandler<ForceModel, ParticleType>::respect_boundaries(ParticleType* simulation_particle) {
    std::array<double, 3> d1;
    std::array<double, 3> d2;
    std::array<bool, 3> directions = {false, false, false};
    std::size_t overlapping_directions = 0;
    for (std::size_t direction = 0; direction != active_directions_.size(); ++direction) {
        if (active_directions_[direction]) {
            d1[direction] = simulation_particle->get_position()[direction] - boundaries_[direction].min;
            d2[direction] = boundaries_[direction].max - simulation_particle->get_position()[direction];
            if (d1[direction] < 0 || d2[direction] < 0) {
                directions[direction] = true;
                ++overlapping_directions;
            }
        }
    }

    if (overlapping_directions > 0) {
        std::size_t mirror_idx;
        Vec3 distance_to_move = Vec3(0., 0., 0.);
        if (overlapping_directions == 1) {
            mirror_idx = std::find(directions.begin(), directions.end(), true) - directions.begin();
        }
        else if (directions[0] && directions[1]) {
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
        auto mirror_particle = get_mirror_particle(simulation_particle, mirror_idx);
        mirror_particle->sum_contact_forces();
        simulation_particle->reset_contact_forces();
        auto p_id = simulation_particle->get_id();
        mirror_particles_[p_id][mirror_idx] = simulation_particle;
        auto particle_idx = (std::find(simulation_particles_.begin(), simulation_particles_.end(), simulation_particle)
                             - simulation_particles_.begin());
        simulation_particles_[particle_idx] = mirror_particle;

        //  Check so that the new mirror particle does not contain any mirror - mirror contacts

        switch (mirror_idx) {
            case 0:
                std::swap(mirror_particles_[p_id][1], mirror_particles_[p_id][3]);
                std::swap(mirror_particles_[p_id][2], mirror_particles_[p_id][4]);
                std::swap(mirror_particles_[p_id][5], mirror_particles_[p_id][6]);
                break;
            case 1:
                std::swap(mirror_particles_[p_id][0], mirror_particles_[p_id][3]);
                std::swap(mirror_particles_[p_id][2], mirror_particles_[p_id][5]);
                std::swap(mirror_particles_[p_id][4], mirror_particles_[p_id][6]);
                break;
            case 2:
                std::swap(mirror_particles_[p_id][0], mirror_particles_[p_id][4]);
                std::swap(mirror_particles_[p_id][1], mirror_particles_[p_id][5]);
                std::swap(mirror_particles_[p_id][3], mirror_particles_[p_id][6]);
                break;
            case 3:
                std::swap(mirror_particles_[p_id][0], mirror_particles_[p_id][1]);
                std::swap(mirror_particles_[p_id][4], mirror_particles_[p_id][5]);
                std::swap(mirror_particles_[p_id][2], mirror_particles_[p_id][6]);
                break;
            case 4:
                std::swap(mirror_particles_[p_id][0], mirror_particles_[p_id][2]);
                std::swap(mirror_particles_[p_id][5], mirror_particles_[p_id][3]);
                std::swap(mirror_particles_[p_id][1], mirror_particles_[p_id][6]);
                break;
            case 5:
                std::swap(mirror_particles_[p_id][0], mirror_particles_[p_id][6]);
                std::swap(mirror_particles_[p_id][1], mirror_particles_[p_id][2]);
                std::swap(mirror_particles_[p_id][3], mirror_particles_[p_id][4]);
                break;
            case 6:
                std::swap(mirror_particles_[p_id][0], mirror_particles_[p_id][5]);
                std::swap(mirror_particles_[p_id][1], mirror_particles_[p_id][4]);
                std::swap(mirror_particles_[p_id][2], mirror_particles_[p_id][3]);
                break;
        }
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::handle_jump_contacts() {
    for (auto& jump_particle: jump_particles_){
        for (auto& c: jump_particle->get_contacts().get_objects()) {
            auto contact_pair = c.first->get_particles();
            auto old_distance = (contact_pair.first->get_position() - contact_pair.second->get_position()).length();
            // The particle p1 is arbitrary chosen as the one belonging to the simulation box
            auto p1 = get_simulation_particle(contact_pair.first->get_id());

            auto p2 = get_simulation_particle(contact_pair.second->get_id());
            for (unsigned i = 0; i != 7; ++i) {
                auto p = get_mirror_particle(p2, i);
                if (p != nullptr) {
                    auto new_distance = (p1->get_position() - p->get_position()).length();
                    if (new_distance - old_distance < 1e-12) {
                        p2 = p;
                    }
                }
            }
            c.first->assign_new_contact_particles(p1, p2);
        }
    }
}


template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::add_periodic_bc(char axis, double boundary_min, double boundary_max) {
    auto axis_idx = direction_idx(axis);
    if (!active_directions_[axis_idx]) {
        active_directions_[axis_idx] = true;
        ++no_active_directions_;
    }
    boundaries_[axis_idx].min = boundary_min;
    boundaries_[axis_idx].max = boundary_max;
}


template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::set_periodic_bc_velocity(char axis, double velocity) {
    auto axis_idx = direction_idx(axis);
    velocities_[axis_idx] = velocity;
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::set_periodic_bc_strain_rate(char axis, double strain_rate) {
    auto axis_idx = direction_idx(axis);
    velocities_[axis_idx] = 0;
    strain_rates_[axis_idx] = strain_rate;
}


template<typename ForceModel, typename ParticleType>
ParticleType* PeriodicBCHandler<ForceModel, ParticleType>::get_mirror_particle(ParticleType* simulation_particle,
                                                                               std::size_t direction) {
    return mirror_particles_[simulation_particle->get_id()][direction];
}

template<typename ForceModel, typename ParticleType>
ParticleType* PeriodicBCHandler<ForceModel, ParticleType>::get_simulation_particle(std::size_t particle_id) {
    return *std::lower_bound(simulation_particles_.begin(), simulation_particles_.end(), particle_id,
                             [](const auto& p1, const auto& rhs) {return p1->get_id() < rhs; });
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::remove_mirror_particle(ParticleType* simulation_particle,
                                                                         std::size_t direction) {
    auto p = mirror_particles_[simulation_particle->get_id()][direction];
    if (p != nullptr) {
        for (auto& c: p->get_contacts().get_objects()) {
            if (c.first->get_particles().first == p ||
                (c.first->get_particles().second != nullptr && c.first->get_particles().second == p)) {
                std::cout << "Problem with mp\n";
                return;
            }
        }
        collision_detector_.remove_particle(p);
        delete p;
        mirror_particles_[simulation_particle->get_id()][direction] = nullptr;
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::create_corner_particles(ParticleType* simulation_particle) {
    position_corner_particle(simulation_particle, 3, {true, true, false});
    position_corner_particle(simulation_particle, 4, {true, false, true});
    position_corner_particle(simulation_particle, 5, {false, true, true});
    position_corner_particle(simulation_particle, 6, {true, true, true});
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::position_corner_particle(ParticleType* simulation_particle,
                                                                           std::size_t direction,
                                                                           const std::array<bool, 3>& mirror_directions) {
    const auto p_id = simulation_particle->get_id();
    bool create = mirror_particles_.count(p_id) > 0 && mirror_particles_[p_id][direction] == nullptr;
    for (unsigned axis = 0; axis != mirror_directions.size(); ++axis) {
        if (create && mirror_directions[axis] ) {
            create = create && (mirror_particles_.count(p_id) > 0) &&  (mirror_particles_[p_id][axis] != nullptr);
        }
    }

    if (create) {
        Vec3 position = Vec3(0, 0, 0);
        for (unsigned axis = 0; axis != mirror_directions.size(); ++axis) {
            position[axis] = !mirror_directions[axis]*simulation_particle->get_position()[axis];
            auto axis_mirior_p = mirror_particles_[p_id][axis];
            if (axis_mirior_p != nullptr) {
                position[axis] += mirror_directions[axis]*axis_mirior_p->get_position()[axis];
            }
        }
        create_mirror_particle(simulation_particle, direction, position);
    }
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::create_periodic_bc_contacts() {
    auto& contacts_to_create = collision_detector_.contacts_to_create();
    for (auto& c_data : contacts_to_create) {
        typename ContactMatrix<Contact<ForceModel, ParticleType>>::PointerType c = nullptr;
        std::size_t id1 = c_data.get_id_pair().first;
        std::size_t id2 = c_data.get_id_pair().second;
        auto p1 = c_data.particle1;
        auto p2 = c_data.particle2;
        if (!(is_mirror_particle(p1) && is_mirror_particle(p2))) {
            auto s = c_data.surface;
            if (s == nullptr) {
                if (is_mirror_particle(p1) && !is_mirror_particle(p2)) {
                    std::swap(p1, p2);
                    std::swap(id1, id2);
                }
                c = contacts_.create_item_inplace(id1, id2, p1, p2, engine_.get_time_increment());
            }
            else {
                c = contacts_.create_item_inplace(id1, id2, p1, s, engine_.get_time_increment());
            }
            if (c != nullptr) {
                p1->add_contact(c, id2, 1.);
                if (p2 != nullptr) {
                    p2->add_contact(c, id1, -1);
                }
                else {
                    s->add_contact(c, id1);
                }
                handle_mirror_particles_add_contact(c, id1, id2, p1, 1);
                handle_mirror_particles_add_contact(c, id2, id1, p2, -1);
            }
        }
    }
    contacts_to_create.clear();

}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::destroy_periodic_bc_contacts() {
    auto& contacts_to_destroy = collision_detector_.contacts_to_destroy();
    for (const auto& c_data : contacts_to_destroy) {
        const auto id1 = c_data.get_id_pair().first;
        const auto id2 = c_data.get_id_pair().second;
        auto p1 = c_data.particle1;
        auto p2 = c_data.particle2;
        auto s = c_data.surface;
        auto jump_contact = false;
        for (auto& jp: jump_particles_) {
            if (jp->get_id() == id1 || jp->get_id() == id2) {
                jump_contact = true;
            }
        }
        if (!jump_contact) {
            if (contacts_.erase(id1, id2)) {
                p1->remove_contact(id2);
                if (s == nullptr) {
                    p2->remove_contact(id1);
                }
                else {
                    s->remove_contact(id1);
                }

                auto handle_mirror_particles = [this](std::size_t id1, std::size_t id2, const ParticleType* p) mutable {
                    if (has_mirror_particle(p)) {
                        for (auto& mp: mirror_particles_[id1]) {
                            if (mp != nullptr) {
                                mp->remove_contact(id2);
                            }
                        }
                    }
                    if (is_mirror_particle(p)) {
                        auto sim_particle_1 = get_simulation_particle(id1);
                        sim_particle_1->remove_contact(id2);
                    }
                };
                handle_mirror_particles(id1, id2, p1);
                handle_mirror_particles(id2, id1, p2);
            }
        }
    }
    contacts_to_destroy.clear();
}


template<typename ForceModel, typename ParticleType>
std::string PeriodicBCHandler<ForceModel, ParticleType>::print_periodic_bc() const {
    std::ostringstream ss;
    ss << boundaries_[0].min << ", " << boundaries_[0].max << ", "
       << boundaries_[1].min << ", " << boundaries_[1].max << ", "
       << boundaries_[2].min << ", " << boundaries_[2].max;
    return ss.str();
}

template<typename ForceModel, typename ParticleType>
bool PeriodicBCHandler<ForceModel, ParticleType>::has_mirror_particle(const ParticleType* particle) const {
    if (particle == nullptr) {
        return false;
    }
    return mirror_particles_.find(particle->get_id()) != mirror_particles_.end();
}

template<typename ForceModel, typename ParticleType>
bool PeriodicBCHandler<ForceModel, ParticleType>::is_mirror_particle(const ParticleType* particle) const {
    if (particle == nullptr) {
        return false;
    }
    auto it = std::lower_bound(simulation_particles_.begin(),  simulation_particles_.end(), particle->get_id(),
                               [](const auto& p1, const auto& rhs) {return p1->get_id() < rhs; });
    return *it != particle;
}

template<typename ForceModel, typename ParticleType>
ParticleType* PeriodicBCHandler<ForceModel, ParticleType>::create_mirror_particle(
        const ParticleType* simulation_particle, std::size_t direction, const Vec3& position, std::size_t collision_id) {
    auto p = new ParticleType(*simulation_particle);
    p->set_position(position);
    p->reset_contact_forces();
    p->set_collision_id(collision_id);
    collision_detector_.add_particle(p);
    mirror_particles_[simulation_particle->get_id()][direction] = p;
    return p;
}

template<typename ForceModel, typename ParticleType>
ParticleType* PeriodicBCHandler<ForceModel, ParticleType>::create_mirror_particle(
        const ParticleType* simulation_particle, std::size_t direction, const Vec3& position) {
    auto p = create_mirror_particle(simulation_particle, direction, position, engine_.get_collision_object_count());
    engine_.increment_collision_counter();
    return p;
}

template<typename ForceModel, typename ParticleType>
std::vector<std::string> PeriodicBCHandler<ForceModel, ParticleType>::mirror_particles_output() const {
    std::vector<std::string> output_data;
    for (const auto& [idx, mp_array]: mirror_particles_) {
        for (const auto mp: mp_array) {
            if (mp != nullptr) {
                output_data.push_back(mp->get_output_string());
            }
        }
    }
    return output_data;
}

template<typename ForceModel, typename ParticleType>
std::vector<std::string> PeriodicBCHandler<ForceModel, ParticleType>::restart_data() const {
    std::vector<std::string> restart_data;
    // Writing boundaries and velocities
    std::ostringstream boundary_info;
    std::ostringstream velocity_info;
    std::ostringstream strain_rate_info;
    boundary_info << "data=boundary";
    velocity_info << "data=velocity";
    strain_rate_info << "data=strain_rate";
    std::string axes = "xyz";
    for (std::size_t dir = 0; dir != boundaries_.size(); ++dir) {
        boundary_info << ", " << named_print(boundaries_[dir].min, axes.substr(dir, 1) + "min")
                      << ", " << named_print(boundaries_[dir].max, axes.substr(dir, 1) + "max");

        strain_rate_info << ", " << named_print(strain_rates_[dir], "e" + axes.substr(dir, 1));
        if (strain_rates_[dir] != 0.) {
            velocity_info << ", " << named_print(0, "v" + axes.substr(dir, 1));
        }
        else {
            velocity_info << ", " << named_print(velocities_[dir], "v" + axes.substr(dir, 1));
        }
    }
    restart_data.push_back(boundary_info.str());
    restart_data.push_back(velocity_info.str());
    restart_data.push_back(strain_rate_info.str());
    restart_data.push_back("data=stretch," + named_print(stretch_, "stretch"));
    for (const auto&[id, mp_array]: mirror_particles_) {
        for (unsigned axis = 0; axis != 7; ++axis) {
            auto mp = mp_array[axis];
            if (mp != nullptr) {
                std::ostringstream ss;
                ss << "data=mirror_particle, axis=" << axis << ", " << mp->restart_data();
                restart_data.push_back(ss.str());
            }
        }
    }

    // Writing all contacts that involve a mirror particle
    for (const auto& contact: contacts_.get_objects_sorted()) {
        const auto& [p1, p2] = contact->get_particles();
        if (is_mirror_particle(p1) || is_mirror_particle(p2)) {
            std::ostringstream ss;
            ss << "data=mirror_contact, pair=";
            auto write_particle_info = [this, &ss](const auto& p) mutable {
                if (is_mirror_particle(p)) {
                    const auto& mp_vec = mirror_particles_.at(p->get_id());
                    std::size_t direction = std::find(mp_vec.begin(), mp_vec.end(), p) - mp_vec.begin();
                    return std::string("mirror") + std::to_string(direction);
                }
                else {
                    if (p == nullptr) {
                        return std::string("surface");
                    }
                    return std::string("particle");
                }
            };
            ss << write_particle_info(p1) << "-" << write_particle_info(p2) << ", ";
            ss << contact->restart_data();
            restart_data.push_back(ss.str());
        }
    }
    return restart_data;
}

template<typename ForceModel, typename ParticleType>
std::size_t PeriodicBCHandler<ForceModel, ParticleType>::direction_idx(char axis) {
    const std::string directions = "xyz";
    const auto axis_idx = directions.find(axis);
    if (axis_idx == std::string::npos) {
        throw std::invalid_argument("axis argument must be x, y or z");
    }
    return axis_idx;
}

template<typename ForceModel, typename ParticleType>
void PeriodicBCHandler<ForceModel, ParticleType>::handle_mirror_particles_add_contact(
        Contact<ForceModel, ParticleType>* contact, std::size_t id1, std::size_t id2, ParticleType* particle,
        int contact_direction) {
    if (has_mirror_particle(particle)) {
        for (auto& mp: mirror_particles_[id1]) {
            if (mp != nullptr) {
                mp->add_contact(contact, id2, contact_direction);
            }
        }
    }
    if (is_mirror_particle(particle)) {
        auto sim_particle_1 = get_simulation_particle(id1);
        sim_particle_1->add_contact(contact, id2, contact_direction);
    }

}