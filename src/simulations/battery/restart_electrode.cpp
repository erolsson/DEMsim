//
// Created by erolsson on 16/10/2020.
//

//
// Created by erolsson on 21/09/2020.
//

#include "../simulations.h"
#include "../../contact_models/viscoelastic.h"
#include "../../engine/engine.h"

void DEM::restart_electrode(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto output_directory = parameters.get_parameter<std::string>("output_dir");

    EngineType simulator = EngineType("compact_restart_file.res");
    simulator.write_restart_file("new_restart.res");
    auto output = simulator.get_output("output_0");
    simulator.remove_output(output);

    auto compaction_output = simulator.create_output(output_directory, 0.001s,
                                                     "compaction_output");
    compaction_output->print_fabric_force_tensor = true;
    compaction_output->print_periodic_bc = true;
    compaction_output->print_mirror_particles = true;
    compaction_output->print_particles = true;
    compaction_output->print_contacts = true;
    compaction_output->print_surface_positions = true;
    compaction_output->print_kinetic_energy = true;
    compaction_output->print_surface_forces = true;

    auto compaction = EngineType::RunForTime(simulator, std::chrono::duration<double>(10.));
    simulator.run(compaction);
}