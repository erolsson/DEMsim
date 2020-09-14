//
// Created by erolsson on 01/09/2020.
//

#include "../simulations.h"

#include "../../contact_models/storakers_mesarovic_johnson.h"
#include "../../engine/engine.h"
#include "../../materials/elastic_ideal_plastic_material.h"
#include <fenv.h>

void DEM::periodic_bc_simulation(const std::string& settings_file_name) {
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    using namespace DEM;
    using ForceModel = StorakersMesarovicJohnson;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto N = parameters.get_parameter<double>("N");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");
    auto filling_density = parameters.get_parameter<double>("filling_density");
    auto aspect_ratio_at_dense = parameters.get_parameter<double>("aspect_ratio_at_dense");

    EngineType simulator(1us);
    auto material = simulator.create_material<ElasticIdealPlasticMaterial>(4800);
    material->E = parameters.get_parameter<double>("E");;
    material->sY = parameters.get_parameter<double>("sY");;
    material->nu = parameters.get_parameter<double>("nu");
    material->mu = parameters.get_parameter<double>("mu");
    material->mu_wall = parameters.get_parameter<double>("mu_wall");
    material->kT = parameters.get_parameter<double>("kT");

    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
    std::sort(particle_radii.rbegin(), particle_radii.rend());

    double particle_volume = 0.;
    for(auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }
    std::cout << "Volume of simulated particles is " << particle_volume << "\n";

    auto box_side = pow(particle_volume, 1./3)/aspect_ratio_at_dense;
    auto box_height = aspect_ratio_at_dense*box_side/filling_density;
    auto particle_positions = random_fill_box(-box_side/2, box_side/2,
                                              -box_side/2, box_side/2,
                                              -box_height/2, box_height/2,
                                              particle_radii);
    for (std::size_t i=0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), material);
    }
    auto output1 = simulator.create_output(output_directory, 0.001s, "output1");
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_contacts = true;
    output1->print_periodic_bc = true;
    output1->print_mirror_particles = true;
    output1->print_fabric_force_tensor = true;

    simulator.add_periodic_boundary_condition('x', -box_side/2, box_side/2);
    simulator.add_periodic_boundary_condition('y', -box_side/2, box_side/2);
    simulator.add_periodic_boundary_condition('z', -box_height/2, box_height/2);

    //simulator.set_periodic_boundary_condition_strain_rate('x',-1.);
    // simulator.set_periodic_boundary_condition_strain_rate('y',-1.);
    simulator.set_periodic_boundary_condition_strain_rate('z',-2/3.);

    simulator.set_mass_scale_factor(1.);

    simulator.set_rotation(false);
    simulator.setup();
    std::chrono::duration<double> restart_time = 0.8s;
    EngineType::RunForTime run_for_time(simulator, restart_time);
    simulator.run(run_for_time);
    simulator.write_restart_file("restart.res");
    simulator.remove_output(output1);
    output1 = simulator.create_output(output_directory + "/original/", 0.001s, "output2");
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_contacts = true;
    output1->print_periodic_bc = true;
    output1->print_mirror_particles = true;
    output1->print_fabric_force_tensor = true;

    EngineType::RunForTime run_for_time2(simulator, 1.s - restart_time);
    simulator.run(run_for_time2);

    EngineType restart_engine("restart.res");
    restart_engine.write_restart_file("restart2.res");
    auto output = restart_engine.get_output("output1");
    restart_engine.remove_output(output);
    output1 = restart_engine.create_output(output_directory + "/restart/", 0.001s, "output2");
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_contacts = true;
    output1->print_periodic_bc = true;
    output1->print_mirror_particles = true;
    output1->print_fabric_force_tensor = true;

    EngineType::RunForTime run_for_time3(restart_engine, 1.s - restart_time);
    restart_engine.run(run_for_time3);
}

