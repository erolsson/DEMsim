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
    auto Cathode_output = simulator.get_output("output_0");
    simulator.remove_output(Cathode_output);
    auto compaction_output = simulator.create_output(output_directory + "/compaction", 0.001s);
    compaction_output->print_particles = true;
    compaction_output->print_surface_positions = true;
    compaction_output->print_kinetic_energy = true;
    compaction_output->print_contacts = true;
    compaction_output->print_surface_forces = true;
    compaction_output->print_fabric_force_tensor =true;
    compaction_output->print_periodic_bc = true;


    auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    double surface_velocity = 0.01;
    top_surface->set_velocity(Vec3(0, 0,0));
    EngineType::RunForTime run_for_time_unload_compact(simulator,5s);
    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.01s);
    simulator.set_mass_scale_factor(10.0);
    simulator.run(run_for_time_unload_compact);
    simulator.write_restart_file(output_directory + "/unload_restart_file.res");






    //Beginning of tryck

    std::cout<<"Biginning of simulation 1"<< std::endl;
    EngineType::RunForTime run_for_time_compact_1(simulator,1.75s);


    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    auto bottom_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("deformable_point_surface_0");

    bottom_surface -> set_in_plane_strain_rates(-0.01, 0.);


    simulator.run(run_for_time_compact_1);

    simulator.write_restart_file(output_directory + "/tryck_1.res");

    //unload extra compaction

    std::cout<<"beginning of unloading 1"<< std::endl;
    //top_surface->set_velocity(Vec3(0, 0, surface_velocity));
    //EngineType::SurfaceNormalForceLess zero_force(top_surface, 0.);
    //simulator.run(run_for_time);

    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    bottom_surface -> set_in_plane_strain_rates(0.01, 0.);
    EngineType::RunForTime run_for_time_relax_1(simulator,1.75s);

    simulator.run(run_for_time_relax_1);
    simulator.write_restart_file(output_directory + "/relaxation_1.res");




    std::cout<<"Biginning of simulation 2"<< std::endl;
    EngineType::RunForTime run_for_time_compact_2(simulator,1.9s);


    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    bottom_surface -> set_in_plane_strain_rates(-0.01, 0.);


    simulator.run(run_for_time_compact_2);

    simulator.write_restart_file(output_directory + "/tryck_2.res");

    //unload extra compaction

    std::cout<<"beginning of unloading 2"<< std::endl;
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    bottom_surface -> set_in_plane_strain_rates(0.01, 0.);
    EngineType::RunForTime run_for_time_relax_2(simulator,1.9s);

    simulator.run(run_for_time_relax_2);
    simulator.write_restart_file(output_directory + "/relaxation_2.res");



    std::cout<<"Biginning of simulation 3"<< std::endl;
    EngineType::RunForTime run_for_time_compact_3(simulator,2.125s);


    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    bottom_surface -> set_in_plane_strain_rates(-0.01, 0.);


    simulator.run(run_for_time_compact_3);

    simulator.write_restart_file(output_directory + "/tryck_3.res");

    //unload extra compaction

    std::cout<<"beginning of unloading 3"<< std::endl;
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    bottom_surface -> set_in_plane_strain_rates(0.01, 0.);
    EngineType::RunForTime run_for_time_relax_3(simulator,2.125s);

    simulator.run(run_for_time_relax_3);
    simulator.write_restart_file(output_directory + "/relaxation_3.res");





    std::cout<<"Biginning of simulation 4"<< std::endl;
    EngineType::RunForTime run_for_time_compact_4(simulator,2.43s);


    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    bottom_surface -> set_in_plane_strain_rates(-0.01, 0.);


    simulator.run(run_for_time_compact_4);

    simulator.write_restart_file(output_directory + "/tryck_4.res");

    //unload extra compaction

    std::cout<<"beginning of unloading 4"<< std::endl;
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    bottom_surface -> set_in_plane_strain_rates(0.01, 0.);
    EngineType::RunForTime run_for_time_relax_4(simulator,2.43s);

    simulator.run(run_for_time_relax_4);
    simulator.write_restart_file(output_directory + "/relaxation_4.res");






    std::cout<<"Biginning of simulation 5"<< std::endl;
    EngineType::RunForTime run_for_time_compact_5(simulator,2.85s);


    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    bottom_surface -> set_in_plane_strain_rates(-0.01, 0.);


    simulator.run(run_for_time_compact_5);

    simulator.write_restart_file(output_directory + "/tryck_5.res");

    //unload extra compaction

    std::cout<<"beginning of unloading 5"<< std::endl;
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    bottom_surface -> set_in_plane_strain_rates(0.01, 0.);
    EngineType::RunForTime run_for_time_relax_5(simulator,2.85s);

    simulator.run(run_for_time_relax_5);
    simulator.write_restart_file(output_directory + "/relaxation_5.res");


}


