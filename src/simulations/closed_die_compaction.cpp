//
// Created by erolsson on 2019-01-05.
//

#include <string>

#include "../engine/engine.h"
#include "../materials/elastic_ideal_plastic_material.h"
#include "../surfaces/point_surface.h"
#include "../utilities/file_reading_functions.h"
#include "../utilities/filling_functions.h"
#include "simulations.h"
#include "../contact_models/storakers_mesarovic_johnson.h"

void DEM::closed_die_compaction(const std::string& settings_file_name){
    using namespace DEM;
    using ForceModel = StorakersMesarovicJohnson;
    using ParticleType = SphericalParticle<ForceModel>;

    // Defines the engine as an engine with StorakersMesarovicJohnson as force model and SphericalParticle as
    // particle type
    using EngineType = Engine<ForceModel, ParticleType>;

    // Makes it possible to write 1s for one second
    using namespace std::chrono_literals;

    // Reads the file defining different parameters which is passed as the second argument to the program
    // The typ of simulation is the first argument
    SimulationParameters parameters(settings_file_name);

    // Reads the parameter N from the simulation file which will be the number of particles
    auto N = parameters.get_parameter<std::size_t>("N");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");

    // The relative density of the "particle gas" at start of the simulation
    auto gas_density = parameters.get_parameter<double>("gas_density");

    // Creates a DEM engine with a time step of one microsecond
    EngineType simulator(1us);

    // Reads a vector with different relative densities where unloading will be made
    auto density_levels = parameters.get_vector<double>( "density_levels" );

    // Creates the material, Elastic - ideal plastic where the elastic part upon loading is assumed to be negligible
    // in the Storakers-Mesarovic-Johnson model
    auto material = simulator.create_material<ElasticIdealPlasticMaterial>(2630.);
    material->sY = parameters.get_parameter<double>("sY");
    material->E = parameters.get_parameter<double>("E");
    material->nu = parameters.get_parameter<double>("nu");

    material->mu = parameters.get_parameter<double>("mu");
    material->mu_wall = parameters.get_parameter<double>("mu_wall");
    material->kT = parameters.get_parameter<double>("kT");

    std::chrono::duration<double> compaction_time {parameters.get_parameter<double>("compaction_time")};
    auto unloading_velocity = parameters.get_parameter<double>("unloading_velocity");
    std::chrono::duration<double> unloading_time {parameters.get_parameter<double>("unloading_time")};

    // Read particle radii from file
    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
    std::sort(particle_radii.rbegin(), particle_radii.rend());

    // The ratio between cylinder diameter and height at zero porosity
    auto aspect_ratio_at_dense = parameters.get_parameter<double>("aspect_ratio_at_dense");

    // ================================================================================================================
    // Creating the particles
    // ================================================================================================================

    //  Calculating the volume of all particles
    double particle_volume = 0.;
    for(auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }

    std::cout << "Volume of simulated particles is " << particle_volume << "\n";

    // Creates the cylinder so that the volume corresponds to an initial density of gas_density
    auto cylinder_radius = pow(4*particle_volume/pi/aspect_ratio_at_dense, 1./3)/2;
    auto cylinder_height = 2*cylinder_radius*aspect_ratio_at_dense/gas_density;
    simulator.create_cylinder(cylinder_radius, Vec3(0, 0, 1), Vec3(0, 0, 0),
                              cylinder_height, true, true);
    std::cout << "The simulated inward_cylinder has a radius of " << cylinder_radius << " and a height of "
              << cylinder_height << "\n";

    // This creates the random particle gas with particles at particle_position
    auto particle_positions = random_fill_cylinder(0, cylinder_height, cylinder_radius, particle_radii);

    // Creating the particles at particle_position
    for (std::size_t i=0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), material);
    }

    // ================================================================================================================
    // Creating the upper and bottom plate by defining two quadratic plates
    // ================================================================================================================

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

    // ================================================================================================================
    // Defining output, gravity and mass scaling
    // ================================================================================================================

    // output every 1 millisecond
    auto output1 = simulator.create_output(output_directory, 0.001s, "output");
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;

    simulator.set_gravity(Vec3(0, 0, -9.820));
    simulator.set_mass_scale_factor(10.);

    // ================================================================================================================
    // Running the simulation
    // ================================================================================================================

    simulator.setup();  // Needed to initialize everything

    // Running everything for 0.1 second to get things moving
    EngineType::RunForTime run_for_time(simulator, 0.1s);
    simulator.run(run_for_time);

    // Run until the fastest particle have a velocity of 0.1 and this is checked every 10 millisecond
    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.01s);
    simulator.run(max_velocity);

    // Move the lid to the uppermost particle
    auto bbox = simulator.get_bounding_box();
    double h = bbox[5];     // The sixth component (zero indexing) is the z_max
    top_surface->move(-Vec3(0, 0, cylinder_height - h), Vec3(0, 0, 0));

    // This is done in several steps
    // 1. Compacting to the density level
    // 2. writing a restart file
    // 3. creating a new simulator (Engine) from that restart and unload that
    // 4. Start over, second iteration in the loop to compact to a new level

    for (unsigned level = 0; level != density_levels.size(); ++level) {
        double h_target = (particle_volume/density_levels[level])/(pi*cylinder_radius*cylinder_radius);
        double surface_velocity = (h_target - top_surface->get_points()[0].z())/(compaction_time.count());
        top_surface->set_velocity(Vec3(0, 0, surface_velocity));
        run_for_time.reset(compaction_time);
        simulator.run(run_for_time);

        std::stringstream restart_file_name;
        restart_file_name << output_directory << "/restart_D=" << density_levels[level] << ".res";
        simulator.write_restart_file(restart_file_name.str());

        auto unloading_simulator = Engine<ForceModel, ParticleType>(restart_file_name.str());
        // Unload the compact
        auto unload_top_surface = unloading_simulator.get_surface<EngineType::PointSurfacePointer>("top_surface");
        unload_top_surface->set_velocity(Vec3(0, 0, unloading_velocity));
        auto compaction_output = unloading_simulator.get_output("output");
        unloading_simulator.remove_output(compaction_output);
        std::stringstream unloading_output_name;
        unloading_output_name << output_directory << "/unload_D=" << density_levels[level];
        auto output2 = unloading_simulator.create_output(unloading_output_name.str(), 0.001s);
        output2->print_particles = true;
        output2->print_kinetic_energy = true;
        output2->print_surface_positions = true;
        output2->print_surface_forces = true;
        restart_file_name << 1;
        unloading_simulator.write_restart_file(restart_file_name.str());
        EngineType::RunForTime time_to_unload(unloading_simulator, unloading_time);
        unloading_simulator.run(time_to_unload);
    }
}
