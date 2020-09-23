//
// Created by erolsson on 21/09/2020.
//

#include "../simulations.h"

#include "../../contact_models/viscoelastic.h"
#include "../../engine/engine.h"

void DEM::battery_rve_compaction(const std::string& settings_file_name)
{
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);
    auto bt = parameters.get_parameter<double>("bt");
    auto compaction_density = parameters.get_parameter<double>("compaction_density");
    auto compaction_time = parameters.get_parameter<double>("compaction_time");
    auto unloading_time = parameters.get_parameter<double>("unloading_time");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");

    EngineType simulator = EngineType(output_directory + "/filling_state.res");
    auto filling_output = simulator.get_output("filling_output");
    simulator.remove_output(filling_output);

    auto compaction_output = simulator.create_output(output_directory + "/compaction", 0.001s,
                                                     "compaction_output");
    compaction_output->print_fabric_force_tensor = true;
    compaction_output->print_periodic_bc = true;
    compaction_output->print_mirror_particles = true;
    compaction_output->print_particles = true;
    compaction_output->print_contacts = true;
    compaction_output->print_surface_positions = true;
    compaction_output->print_kinetic_energy = true;
    compaction_output->print_surface_forces = true;

    auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");

    // Moving the top surface to the top position of the particles + bt
    auto bounding_box = simulator.get_bounding_box();
    auto z_max = bounding_box[5] + 1.01*bt;
    top_surface->move(Vec3(0, 0, z_max - top_surface->get_points()[0].z()), Vec3(0, 0, 0));

    // Calculating the volume of the particles to calculate the final surface position
    auto simulation_particles = simulator.get_particles();
    double particle_volume = 0;
    for (const auto& p: simulation_particles) {
        auto r = p->get_radius();
        particle_volume += 4*pi*r*r*r/3;
    }

    // Calculate the height of the final box
    auto box_size = simulator.get_periodic_boundaries()[0].max - simulator.get_periodic_boundaries()[0].min;
    auto compaction_height = particle_volume/compaction_density/box_size/box_size;
    auto compaction_velocity = (z_max - compaction_height)/compaction_time;
    top_surface->set_velocity(Vec3(0, 0, -compaction_velocity));
    // simulator.set_rotation(false);
    auto compaction = EngineType::RunForTime(simulator, std::chrono::duration<double>(compaction_time));
    simulator.run(compaction);
    auto unloading = EngineType::RunForTime(simulator, std::chrono::duration<double>(unloading_time));
    top_surface->set_velocity(Vec3(0, 0, compaction_velocity));
    simulator.run(unloading);
    auto max_velocity = EngineType::ParticleVelocityLess(simulator, 0.01, 1us);
    simulator.run(max_velocity);
    simulator.write_restart_file(output_directory + "/compaction_state.res");
}