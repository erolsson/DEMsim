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
    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
    auto Cathode_output = simulator.get_output("output_0");
    simulator.remove_output(Cathode_output);
    auto compaction_output = simulator.create_output(output_directory + "/fixar_porosity", 0.005s);
    compaction_output->print_particles = true;
    compaction_output->print_surface_positions = true;
    compaction_output->print_kinetic_energy = true;
    compaction_output->print_contacts = true;
    compaction_output->print_surface_forces = true;
    compaction_output->print_fabric_force_tensor =true;
    compaction_output->print_periodic_bc = true;
    compaction_output->print_mirror_particles= true;
    //mat-> adhesive = true;
    auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("deformable_point_surface_0");
    //double surface_velocity = 0.01;
    top_surface->set_velocity(Vec3(0, 0,0));
    EngineType::RunForTime run_for_time_unload_compact(simulator,1s);
    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.01s);
    simulator.set_mass_scale_factor(10.0);
    mat-> adhesive = true;
    simulator.run(run_for_time_unload_compact);
    //simulator.write_restart_file(output_directory + "/stability_check.res");

    double surface_velocity = 0.01;


    EngineType::RunForTime run_for_time_relax(simulator,15s);
    //simulator.set_rotation(false);
   // mat-> adhesive = true;
    top_surface->set_velocity(Vec3(0, 0, surface_velocity));
    simulator.run(run_for_time_relax);
    //simulator.write_restart_file(output_directory + "/new_porosity.res");

    std::cout<<"Biginning of simulation 4"<< std::endl;
    EngineType::RunForTime run_for_time_compact_4(simulator,0.003s);

    simulator.set_periodic_boundary_condition_strain_rate('x',-1.0);
    deformable_surface -> set_in_plane_strain_rates(-1.0, 0.);
    //simulator.set_mass_scale_factor(10.0);
    mat-> adhesive = true;
    simulator.run(run_for_time_compact_4);

    simulator.write_restart_file(output_directory + "/tryck_4.res");

    //unload extra compaction

    std::cout<<"beginning of unloading 4"<< std::endl;
    simulator.set_periodic_boundary_condition_strain_rate('x',0.);
    deformable_surface -> set_in_plane_strain_rates(0., 0.);
    EngineType::RunForTime run_for_time_relax_4(simulator,50s);
    //simulator.set_mass_scale_factor(1.0);
    mat-> adhesive = true;
    simulator.run(run_for_time_relax_4);
    simulator.write_restart_file(output_directory + "/relaxation_4.res");



}


