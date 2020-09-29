//
// Created by elahe on 2020-09-11.
//
#include "simulations.h"

#include "../engine/engine.h"
#include "../contact_models/viscoelastic.h"
#include "../particles/spherical_particle.h"
#include "../materials/electrode_material.h"
#include "../utilities/vec3.h"

void DEM::Cathode_mechanical_simulations(const std::string &settings_file_name) {
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);
    auto restart_file_name = parameters.get_parameter<std::string>("restart_file_name");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto simulator = EngineType(restart_file_name);
    //auto mat = simulator.create_material<ElectrodeMaterial>(4800);
    //mat->active_particle_height = parameters.get_parameter<double>("active_particle_height");


    //More compaction
    EngineType::RunForTime run_for_time(simulator, 1s);
    //auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("point_surface_4501");
    //double surface_velocity = 0.01;

    auto Cathode_output = simulator.get_output("filling_output");
    simulator.remove_output(Cathode_output);
    auto compaction_output = simulator.create_output(output_directory + "/compaction", 0.001s);
    compaction_output->print_particles = true;
    compaction_output->print_surface_positions = true;
    compaction_output->print_kinetic_energy = true;
    compaction_output->print_contacts = true;
    compaction_output->print_surface_forces = true;
    compaction_output->print_fabric_force_tensor =true;

    //Move the side_lid
    //simulator.add_periodic_boundary_condition('x', -box_side/2, box_side/2);
    //simulator.add_periodic_boundary_condition('y', -box_side/2, box_side/2);
    std::cout<<"Biginning of simulations"<< std::endl;;
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.001);
    simulator.run(run_for_time);



    //unload extra compaction

    std::cout<<"beginning of unloading"<< std::endl;
    //top_surface->set_velocity(Vec3(0, 0, surface_velocity));
    //EngineType::SurfaceNormalForceLess zero_force(top_surface, 0.);
    //simulator.run(run_for_time);

    simulator.set_periodic_boundary_condition_strain_rate('x',-0.001);
    simulator.run(run_for_time);
    simulator.write_restart_file(output_directory + "/Cathode_final_compacted");
}


