//
// Created by erolsson on 21/08/2020.
//

#include "../simulations.h"

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../particles/spherical_particle.h"
#include "../../materials/electrode_material.h"
#include "../../utilities/vec3.h"

void DEM::electrode_cylinder_compaction(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto filling_file_name = parameters.get_parameter<std::string>("filling_file_name");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto compaction_velocity = parameters.get_parameter<double>("compaction_velocity");
    auto unload_velocity = parameters.get_parameter<double>("unload_velocity");
    auto compaction_time = std::chrono::duration<double>(parameters.get_parameter<double>("compaction_time"));
    auto unloading_time = std::chrono::duration<double>(parameters.get_parameter<double>("unload_time"));
    auto simulator = EngineType(filling_file_name);
    auto bbox = simulator.get_bounding_box();
    auto z_max = bbox[5] + 3e-3;
    auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_surface");
    top_surface->move(Vec3(0, 0, z_max -  top_surface->get_points()[0].z()),
                       Vec3(0, 0, 0));
    top_surface->set_velocity(Vec3(0, 0, -compaction_velocity));

    auto filling_output = simulator.get_output("output");
    simulator.remove_output(filling_output);
    auto compaction_output = simulator.create_output(output_directory + "/compaction", 0.001s);
    compaction_output->print_particles = true;
    compaction_output->print_surface_positions = true;
    compaction_output->print_kinetic_energy = true;
    compaction_output->print_contacts = true;
    compaction_output->print_surface_forces = true;

    auto compaction_time_runner = EngineType::RunForTime(simulator, compaction_time);
    simulator.run(compaction_time_runner);

    auto unloading_time_runner = EngineType::RunForTime(simulator, unloading_time);
    top_surface->set_velocity(Vec3(0, 0, unload_velocity));
    simulator.run(unloading_time_runner);
    simulator.write_restart_file(output_directory + "/compact_restart_file.res");
}