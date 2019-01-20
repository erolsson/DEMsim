//
// Created by erolsson on 2019-01-05.
//

#include <atomic>
#include <string>

#include "../engine/engine.h"
#include "../materials/elastic_ideal_plastic_material.h"
#include "../utilities/file_reading_functions.h"
#include "../utilities/filling_functions.h"
#include "simulations.h"
#include "../contact_models/storakers_mesarovic_johnson.h"

void DEM::closed_die_compaction(const std::string& settings_file_name){
    using namespace DEM;
    using ForceModel = StorakersMesarovicJohnson;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;

    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);

    auto N = parameters.get<std::size_t>("N");
    auto output_directory = parameters.get<std::string>("output_dir");
    auto particle_file = parameters.get<std::string>("radius_file");
    auto filling_density = parameters.get<double>("filling_density");

    EngineType simulator(1us);

    auto material = simulator.create_material<ElasticIdealPlasticMaterial>(2630.);
    material->sY = parameters.get<double>("sY");
    material->E = parameters.get<double>("E");
    material->nu = parameters.get<double>("nu");

    material->mu = parameters.get<double>("mu");
    material->mu_wall = parameters.get<double>("mu_wall");
    material->kT = parameters.get<double>("kT");

    auto target_relative_density = parameters.get<double>("density");
    std::chrono::duration<double> compaction_time {parameters.get<double>("compaction_time")};
    auto unloading_velocity = parameters.get<double>("unloading_velocity");
    std::chrono::duration<double> unloading_time {parameters.get<double>("unloading_time")};


    // Read particle radii from file
    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
    std::sort(particle_radii.rbegin(), particle_radii.rend());
    auto aspect_ratio_at_dense = parameters.get<double>("aspect_ratio_at_dense");
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
    auto particle_positions = random_fill_cylinder(0, cylinder_height, cylinder_radius, particle_radii);

    for (std::size_t i=0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), material);
    }
    // simulator.create_particle(0.005, Vec3(0, 0.0, 0.005), Vec3(0,0,0), material);

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

    auto top_surface = simulator.create_point_surface(top_points, true);
    std::cout << "Normal of top surface is " << top_surface->get_normal() << std::endl;

    auto output1 = simulator.create_output(output_directory, 0.001s);
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;

    simulator.set_gravity(Vec3(0, 0, -9.820));
    simulator.set_mass_scale_factor(10.);
    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, 0.1s);
    /*
    simulator.run(run_for_time);
    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.01s);
    simulator.run(max_velocity);

    // Move the lid to the uppermost particle
    auto bbox = simulator.get_bounding_box();
    double h = bbox[5];
    top_surface->move(-Vec3(0, 0, cylinder_height - h), Vec3(0, 0, 0));

    // Compress the compact
    double h_target = (particle_volume/target_relative_density)/(pi*cylinder_radius*cylinder_radius);
    double surface_velocity = (h_target - h)/(compaction_time.count());
    top_surface->set_velocity(Vec3(0, 0, surface_velocity));
    run_for_time.reset(compaction_time);
    simulator.run(run_for_time);

    // Unload the compact
    top_surface->set_velocity(Vec3(0, 0, unloading_velocity));
    run_for_time.reset(unloading_time);
    simulator.run(run_for_time);
     */
}
