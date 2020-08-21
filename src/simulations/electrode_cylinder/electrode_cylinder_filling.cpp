//
// Created by erolsson on 14/08/2020.
//

#include "../simulations.h"

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../particles/spherical_particle.h"
#include "../../materials/electrode_material.h"

void DEM::electrode_cylinder_filling(const std::string& settings_file_name) {
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
    auto filling_density = parameters.get_parameter<double>("filling_density");

    simulator.setup(mat->bt);
    // auto target_relative_density = parameters.get_parameter<double>("density");
    std::chrono::duration<double> compaction_time {parameters.get_parameter<double>("compaction_time")};
    auto unloading_velocity = parameters.get_parameter<double>("unloading_velocity");
    std::chrono::duration<double> unloading_time {parameters.get_parameter<double>("unloading_time")};

    // Read particle radii from file
    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
    std::sort(particle_radii.rbegin(), particle_radii.rend());
    auto aspect_ratio_at_dense = parameters.get_parameter<double>("aspect_ratio_at_dense");
    double particle_volume = 0.;
    for(auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }
    std::cout << "Volume of simulated particles is " << particle_volume << "\n";
    auto cylinder_radius = pow(4*particle_volume/pi/aspect_ratio_at_dense, 1./3)/2;
    auto cylinder_height = 2*cylinder_radius*aspect_ratio_at_dense/filling_density;
    simulator.create_cylinder(cylinder_radius, Vec3(0, 0, 1), Vec3(0, 0, 0), cylinder_height, true, true);
    std::cout << "The simulated inward_cylinder has a radius of " << cylinder_radius << " and a height of "
              << cylinder_height << "\n";
    auto particle_positions = random_fill_cylinder(0, cylinder_height, cylinder_radius, particle_radii, mat->bt);

    for (std::size_t i=0; i != N; ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), mat);
    }

    // Creating The bottom plate surface
    Vec3 p1(-cylinder_radius, -cylinder_radius, 0.);
    Vec3 p2(-cylinder_radius,  cylinder_radius, 0.);
    Vec3 p3(cylinder_radius,   cylinder_radius, 0.);
    Vec3 p4(cylinder_radius,  -cylinder_radius, 0.);

    // Creating The top plate surface
    Vec3 p5(-cylinder_radius, -cylinder_radius, cylinder_height);
    Vec3 p6(-cylinder_radius,  cylinder_radius, cylinder_height);
    Vec3 p7(cylinder_radius,   cylinder_radius, cylinder_height);
    Vec3 p8(cylinder_radius,  -cylinder_radius, cylinder_height);

    std::vector<Vec3> bottom_points{p4, p3, p2, p1};
    std::vector<Vec3> top_points{p5, p6, p7, p8};

    auto bottom_surface = simulator.create_point_surface(bottom_points, true);
    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << std::endl;

    auto top_surface = simulator.create_point_surface(top_points, true, "top_surface");
    std::cout << "Normal of top surface is " << top_surface->get_normal() << std::endl;

    auto output1 = simulator.create_output(output_directory, 0.001s, "output");
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_contacts = true;

    simulator.set_gravity(Vec3(0, 0, -9.820));
    simulator.set_mass_scale_factor(10.);
    simulator.setup(2*mat->bt);
    simulator.write_restart_file(output_directory + "/starting_state.res");
    EngineType::RunForTime run_for_time(simulator, 0.1s);
    simulator.run(run_for_time);

    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.01s);
    simulator.run(max_velocity);
    simulator.write_restart_file(output_directory + "/resting_state.res");
}