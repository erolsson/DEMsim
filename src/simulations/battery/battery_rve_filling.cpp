//
// Created by erolsson on 16/09/2020.
//

#include "../simulations.h"

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../particles/spherical_particle.h"
#include "../../materials/electrode_material.h"

void DEM::battery_rve_filling(const std::string& settings_file_name)
{
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto N = parameters.get_parameter<double>("N");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");

    EngineType simulator(1us);
    auto mat = simulator.create_material<ElectrodeMaterial>(4800);

    mat->E = parameters.get_parameter<double>("E");
    mat->kT = parameters.get_parameter<double>("kT");
    mat->Ep = parameters.get_parameter<double>("Ep");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->fb = parameters.get_parameter<double>("fb");
    mat->nup = parameters.get_parameter<double>("nup");
    mat->mu = parameters.get_parameter<double>("mu");
    mat->mu_wall = parameters.get_parameter<double>("mu_wall");
    mat->tau_i = parameters.get_vector<double>("tau_i");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->bt = parameters.get_parameter<double>("bt");
    mat->mu_binder = parameters.get_parameter<double>("mu_binder");
    auto particle_density_at_filling = parameters.get_parameter<double>("filling_density");
    auto particle_density_at_cube = parameters.get_parameter<double>("particle_density_at_cube");

    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
    std::sort(particle_radii.rbegin(), particle_radii.rend());

    double particle_volume = 0.;
    for(const auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }
    std::cout << "Volume of simulated particles is " << particle_volume << "\n";

    auto box_side = pow(particle_volume/particle_density_at_cube, 1./3);
    auto box_height = particle_density_at_cube*box_side/particle_density_at_filling;

    auto p1 = Vec3(-box_side/2, -box_side/2, 0);
    auto p2 = Vec3( box_side/2, -box_side/2, 0);
    auto p3 = Vec3( box_side/2,  box_side/2, 0);
    auto p4 = Vec3(-box_side/2,  box_side/2, 0);

    auto p5 = Vec3(-box_side/2, -box_side/2, box_height);
    auto p6 = Vec3( box_side/2, -box_side/2, box_height);
    auto p7 = Vec3( box_side/2,  box_side/2, box_height);
    auto p8 = Vec3(-box_side/2,  box_side/2, box_height);

    std::vector<Vec3> bottom_points{p1, p2, p3, p4};
    std::vector<Vec3> top_points{p8, p7, p6, p5};

    auto particle_positions = random_fill_box(-box_side/2, box_side/2, -box_side/2, box_side/2,
                                              0, box_height, particle_radii, mat->bt);
    auto bottom_surface = simulator.create_deformable_point_surface(bottom_points, "bottom_plate");
    auto top_surface = simulator.create_point_surface(top_points, true, "top_plate", false);
    std::cout << "Normal of bottom surface: " << bottom_surface->get_normal() << "\n";
    std::cout << "Normal of bottom surface: " << top_surface->get_normal() << "\n";

    for (std::size_t i = 0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), mat);
    }
    auto filling_output = simulator.create_output(output_directory + "/filling/", 0.001s,
                                                  "filling_output");

    filling_output->print_particles = true;
    filling_output->print_kinetic_energy = true;
    filling_output->print_surface_positions = true;
    filling_output->print_surface_forces = true;
    filling_output->print_contacts = true;
    filling_output->print_periodic_bc = true;
    filling_output->print_mirror_particles = true;

    simulator.add_periodic_boundary_condition('x', -box_side/2, box_side/2);
    simulator.add_periodic_boundary_condition('y', -box_side/2, box_side/2);

    simulator.set_gravity(Vec3(0, 0, -9.82));
    simulator.setup(1.01*mat->bt);
    EngineType::RunForTime run_for_time(simulator, 0.1s);
    simulator.run(run_for_time);

    auto max_velocity = EngineType::ParticleVelocityLess(simulator, 0.01, 1us);
    simulator.run(max_velocity);
    simulator.write_restart_file(output_directory + "/filling_state.res");
}