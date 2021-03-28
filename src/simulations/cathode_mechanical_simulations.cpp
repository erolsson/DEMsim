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
    auto compaction_output = simulator.create_output(output_directory + "/unload_restart_file", 0.005s);
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
    EngineType::RunForTime run_for_time_pause(simulator,5.0s);
    simulator.run(run_for_time_pause);

    EngineType::RunForTime run_for_time_compact_10(simulator,0.125s);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_10);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_10(simulator,0.125s);
    simulator.run(run_for_time_relax_10);

    EngineType::RunForTime run_for_time_compact_9(simulator,0.25s);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_9);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_9(simulator,0.125s);
    simulator.run(run_for_time_relax_9);

    EngineType::RunForTime run_for_time_compact_8(simulator,0.375s);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_8);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_8(simulator,0.25s);
    simulator.run(run_for_time_relax_8);

    EngineType::RunForTime run_for_time_compact_7(simulator,0.5s);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_7);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_7(simulator,0.375s);
    simulator.run(run_for_time_relax_7);

    EngineType::RunForTime run_for_time_compact_6(simulator,0.625s);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_6);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_6(simulator,0.5s);
    simulator.run(run_for_time_relax_6);

    EngineType::RunForTime run_for_time_compact_5(simulator,0.6s);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_5);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_5(simulator,0.55s);
    simulator.run(run_for_time_relax_5);

    EngineType::RunForTime run_for_time_compact_4(simulator,0.68s);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_4);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_4(simulator,0.615s);
    simulator.run(run_for_time_relax_4);

    EngineType::RunForTime run_for_time_compact_3(simulator,0.525s);
    simulator.set_mass_scale_factor(1.0);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_3);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_3(simulator,0.615s);
    simulator.run(run_for_time_relax_3);

    EngineType::RunForTime run_for_time_compact_2(simulator,0.785s);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_2);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_2(simulator,0.7s);
    simulator.run(run_for_time_relax_2);

    EngineType::RunForTime run_for_time_compact_1(simulator,0.95s);
    simulator.set_periodic_boundary_condition_strain_rate('x',0.01);
    deformable_surface -> set_in_plane_strain_rates(0.01, 0.);
    simulator.run(run_for_time_compact_1);
    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface -> set_in_plane_strain_rates(-0.01, 0.);
    EngineType::RunForTime run_for_time_relax_1(simulator,0.8250s);
    simulator.run(run_for_time_relax_1);





}


