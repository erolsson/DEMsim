//
// Created by erolsson on 05/11/2020.
//

#include "../simulations.h"

#include "../../contact_models/porous_electrode_contact.h"
#include "../../materials/porous_electrode_material.h"
#include "../../particles/spherical_particle.h"
#include "../../engine/engine.h"
#include "../../utilities/filling_functions.h"

void DEM::porous_electrode_rve(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = PorousElectrodeContact;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto radius = parameters.get_parameter<double>("R");
    auto N = parameters.get_parameter<unsigned >("N");
    auto output_directory = parameters.get_parameter<std::string>("output_directory");
    auto initial_particle_density = parameters.get_parameter<double>("initial_particle_density");
    auto electrode_particle_density = parameters.get_parameter<double>("electrode_particle_density");
    auto compaction_time = parameters.get_parameter<double>("compaction_time");

    EngineType simulator(1us);
    auto mat = simulator.create_material<PorousElectrodeMaterial>(4800);

    mat->E_binder = parameters.get_parameter<double>("E_binder");
    mat->v_binder = parameters.get_parameter<double>("v_binder");
    mat->E_particle = parameters.get_parameter<double>("E_particle");
    mat->v_particle = parameters.get_parameter<double>("v_particle");
    mat->binder_thickness = parameters.get_parameter<double>("binder_thickness");
    mat->binder_radius_fraction =parameters.get_parameter<double>("binder_radius_fraction");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->tau_i =parameters.get_vector<double>("tau_i");
    mat->fraction_binder_contacts = parameters.get_parameter<double>("fraction_binder_contacts");

    auto particle_radii = std::vector<double>(N, radius);
    auto particle_volume = 0.;
    for (const auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }
    auto initial_box_volume = particle_volume/initial_particle_density;
    auto box_side = pow(initial_box_volume, 1./3);
    auto particle_positions = random_fill_box_periodic(-box_side/2, box_side/2, -box_side/2, box_side/2,
                                                       -box_side/2, box_side/2,  particle_radii,
                                                       mat->binder_thickness, "xyz");

    for (std::size_t i = 0; i!= N; ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(), mat);
    }
    simulator.add_periodic_boundary_condition('x', -box_side/2, box_side/2);
    simulator.add_periodic_boundary_condition('y', -box_side/2, box_side/2);
    simulator.add_periodic_boundary_condition('z', -box_side/2, box_side/2);
    auto final_box_volume = particle_volume/electrode_particle_density;
    auto final_box_side = pow(final_box_volume, 1./3);

    auto boundary_velocity = (final_box_side - box_side)/2/compaction_time;
    simulator.set_periodic_boundary_condition_velocity('x', boundary_velocity);
    simulator.set_periodic_boundary_condition_velocity('y', boundary_velocity);
    simulator.set_periodic_boundary_condition_velocity('z', boundary_velocity);

    auto compaction_output = simulator.create_output(output_directory + "/compaction/", 0.01s,
                                                     "compaction_output");

    compaction_output->print_particles = true;
    compaction_output->print_kinetic_energy = true;
    compaction_output->print_contacts = true;
    compaction_output->print_periodic_bc = true;
    compaction_output->print_mirror_particles = true;
    compaction_output->print_fabric_force_tensor = true;

    simulator.setup(1.01*mat->binder_thickness);
    EngineType::RunForTime run_for_time(simulator, std::chrono::duration<double>(compaction_time));
    simulator.run(run_for_time);
    simulator.write_restart_file(output_directory + "/compaction/restart.res");
}
