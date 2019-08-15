//
// Created by erolsson on 12/08/2019.
//

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

    auto output1 = simulator.create_output(output_directory, 0.001s);
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;

    simulator.set_gravity(Vec3(0, 0, -9.820));
    simulator.set_mass_scale_factor(1.);
    std::vector<ParticleType*> particles;
    for (unsigned layer = 0; layer != 1; ++layer) {
        std::vector<double > layer_particles;
        layer_particles.assign(particle_radii.begin() + layer*particles_per_layer,
                particle_radii.begin() + (layer+1)*particles_per_layer);
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
                particle_bead_height+main_cylinder_height, main_cylinder_radius, layer_particles);
        
        for (std::size_t i = 0; i != particle_positions.size(); ++i) {
            ParticleType* p = simulator.create_particle(layer_particles[i], particle_positions[i],
                   Vec3(0, 0, 0), mat);
            particles.push_back(p);
        }
        
        // ParticleType* p = simulator.create_particle(0.00625, Vec3(0, 0, 0.1),
        //                                            Vec3(0, 0, 0), mat);
        // particles.push_back(p);
        simulator.setup();
        EngineType::RunForTime run_for_time01(simulator, 5s);
        simulator.run(run_for_time01);

        EngineType::ParticleVelocityLess max_velocity(simulator, 0.1, 0.01s);
        simulator.run(max_velocity);

        // Moving the hammer to the impact position
        double hammer_position_phi = 0.;
        double hammer_position_x = hammer_radius*cos(hammer_position_phi);
        double hammer_position_y = hammer_radius*sin(hammer_position_phi);

        // Finding the highest material point under the hammer position
        double z_max = 0;
        for (const auto& p: particles) {
            const auto& pos = p->get_position();
            double x0 = pos.x()-hammer_position_x;
            double y0 = pos.y()-hammer_position_y;
            if (pow(x0, 2) + pow(y0, 2) > hammer_radius*hammer_radius) {
                if (pos.z() + p->get_radius() > z_max) {
                    z_max = pos.z() + p->get_radius();
                }
            }
        }
        std::cout << "z_max at hammer position " << hammer_position_x << ", " << hammer_position_y << " is "
                  << z_max << "\n";
        Vec3 new_hammer_position = Vec3(hammer_position_x, hammer_position_y, z_max + 0.457) - hammer->get_position();
        hammer->move(new_hammer_position, Vec3(0,0 ,0));

        hammer->set_mass(4.54);
        simulator.set_force_control_on_surface(hammer, 'z');
        simulator.run(run_for_time01);
        simulator.run(max_velocity);
    }
}