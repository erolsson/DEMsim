//
// Created by erolsson on 12/08/2019.
//

#include <sstream>

#include "simulations.h"

#include "../engine/engine.h"

#include "../contact_models/stone_material_contact.h"
#include "../materials/stone_material.h"

#include "../utilities/file_reading_functions.h"
#include "../utilities/filling_functions.h"

void DEM::proctor_test(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = StoneMaterialContact;
    using ParticleType = FractureableSphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;

    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto output_directory = parameters.get<std::string>("output_dir");
    auto particle_file = parameters.get<std::string>("radius_file");
    auto gas_density = parameters.get<double>("gas_density");
    auto filling_density = parameters.get<double>("filling_density");

    EngineType simulator(1us);

    auto mat = simulator.create_material<StoneMaterial>(2370.);
    mat->E = parameters.get<double>("E");
    mat->nu = parameters.get<double>("nu");
    mat->unloading_exponent = parameters.get<double>("unloading_exponent");
    mat->mu = parameters.get<double>("mu");
    mat->mu_wall = parameters.get<double>("mu_wall");
    mat->weibull_fracture_stress = parameters.get<double>("weibull_fracture_stress");
    mat->weibull_exponent = parameters.get<double>("weibull_exponent");

    auto particle_radii = read_vector_from_file<double>(particle_file);

    double particle_volume = 0.;
    for(auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }

    auto mean_particle_volume = particle_volume/particle_radii.size();

    double main_cylinder_height = 0.11643;
    double main_cylinder_radius = 0.0254*2;

    simulator.create_cylinder(main_cylinder_radius, Vec3(0, 0, 1), Vec3(0, 0,0),
            main_cylinder_height, true, true, false);

    auto cylinder_volume = main_cylinder_height*main_cylinder_radius*main_cylinder_radius*3.1415;
    int particles_per_layer = int(cylinder_volume/5*filling_density/mean_particle_volume);
    std::cout << "Particles per layer: " << particles_per_layer << "\n";
    double hammer_radius = 0.025;
    double hammer_height = 0.08;
    Vec3 hammer_resting_position = Vec3(0., 0., -1.01*hammer_height);
    auto hammer = simulator.create_cylinder(hammer_radius, Vec3(0, 0, 1),
            hammer_resting_position, hammer_height,
            false, false, true);

    // Creating The bottom plate surface
    Vec3 p1(-main_cylinder_radius, -main_cylinder_radius, 0.);
    Vec3 p2(-main_cylinder_radius,  main_cylinder_radius, 0.);
    Vec3 p3(main_cylinder_radius,   main_cylinder_radius, 0.);
    Vec3 p4(main_cylinder_radius,  -main_cylinder_radius, 0.);

    std::vector<Vec3> bottom_points{p4, p3, p2, p1};
    auto bottom_surface = simulator.create_point_surface(bottom_points, true);
    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << std::endl;

    auto particle_output = simulator.create_output(output_directory + "/animation", 0.01s);
    particle_output->print_particles = true;
    particle_output->print_surface_positions = true;

    simulator.set_gravity(Vec3(0, 0, -9.820));
    simulator.set_mass_scale_factor(1.);
    std::vector<ParticleType*> particles;
    simulator.set_rotation(false);


    for (unsigned layer = 0; layer != 5; ++layer) {
        std::vector<double> layer_particles;
        layer_particles.assign(particle_radii.begin() + layer*particles_per_layer,
                               particle_radii.begin() + (layer + 1)*particles_per_layer);
        std::sort(layer_particles.begin(), layer_particles.end());
        double layer_volume = 0.;
        for (auto r: layer_particles) {
            layer_volume += 4*3.1415*r*r*r/3;
        }
        std::cout << "Stone volume of layer: " << layer_volume << "\n";
        auto particle_bead_height = layer_volume/gas_density/
                                    main_cylinder_radius/main_cylinder_radius/3.1415;
        std::cout << "starting height of layer: " << particle_bead_height << "\n";
        auto particle_positions = random_fill_cylinder(main_cylinder_height,
                                                       particle_bead_height + main_cylinder_height,
                                                       main_cylinder_radius, layer_particles);

        for (std::size_t i = 0; i != particle_positions.size(); ++i) {
            ParticleType *p = simulator.create_particle(layer_particles[i], particle_positions[i],
                                                        Vec3(0, 0, 0), mat);
            particles.push_back(p);
        }

        std::stringstream directory_name;
        directory_name << output_directory << "/layer_" << layer << "/filling";
        auto output_filling = simulator.create_output(directory_name.str(), 0.001s);
        output_filling->print_kinetic_energy = true;
        output_filling->print_surface_forces = true;
        output_filling->print_particle_cracks = true;

        simulator.setup();
        EngineType::RunForTime run_for_time01(simulator, 0.1s);
        simulator.run(run_for_time01);

        EngineType::ParticleVelocityLess max_velocity(simulator, 0.1, 0.01s);
        simulator.run(max_velocity);
        // Moving the hammer to the impact position
        simulator.remove_output(output_filling);
        for (unsigned stroke = 0; stroke != 25; ++stroke) {
            double hammer_position_phi = 2*3.1415/25*stroke;
            double hammer_position_x = hammer_radius*cos(hammer_position_phi);
            double hammer_position_y = hammer_radius*sin(hammer_position_phi);

            // Finding the highest material point under the hammer position
            double z_max = 0;
            for (const auto& p: particles) {
                const auto& pos = p->get_position();
                double x0 = pos.x() - hammer_position_x;
                double y0 = pos.y() - hammer_position_y;
                if (pow(x0, 2) + pow(y0, 2) > hammer_radius*hammer_radius) {
                    if (pos.z() + p->get_radius() > z_max) {
                        z_max = pos.z() + p->get_radius();
                    }
                }
            }
            std::cout << "z_max at hammer position " << hammer_position_x << ", " << hammer_position_y << " is "
                      << z_max << "\n";
            Vec3 new_hammer_position =
                    Vec3(hammer_position_x, hammer_position_y, z_max + 0.457) - hammer->get_position();
            hammer->move(new_hammer_position, Vec3(0, 0, 0));

            hammer->set_mass(4.54);
            simulator.set_force_control_on_surface(hammer, 'z');
            std::chrono::duration<double> falling_time(1.1*sqrt(2*0.457/9.82));
            std::cout << "Hammer falling time is: " << falling_time.count() << std::endl;

            directory_name.str(std::string());
            directory_name << output_directory << "/layer_" << layer << "/stroke_" << stroke;
            auto output_stroke = simulator.create_output(directory_name.str(), 0.001s);
            output_stroke->print_kinetic_energy = true;
            output_stroke->print_surface_forces = true;
            output_stroke->print_particle_cracks = true;

            EngineType::RunForTime run_for_time_falling(simulator, falling_time);
            EngineType::SurfaceVelocityLessThan max_velocity_surface(0.01, hammer);

            EngineType::CombinedConditions run_condition({&max_velocity_surface,
                                                          &run_for_time_falling});
            simulator.run(run_condition);

            hammer->move(Vec3(0, 0, -1.01*hammer_height) - hammer->get_position(), Vec3(0, 0, 0));
            simulator.remove_force_control_on_surface(hammer, 'z');

            simulator.run(max_velocity);
            /*
            EngineType::ObjectVelocityLess max_velocity_object(simulator, 0.1, 0.01s);
            EngineType::SurfaceNormalForceWithinInterval surface_force(simulator, hammer, 10., 100., 0.1s);
            EngineType::CombinedConditions run_condition({&max_velocity_object,
                                                          &surface_force,
                                                          &run_for_time_falling});
            simulator.run(run_condition);
             */
            simulator.remove_output(output_stroke);
        }

    }
}