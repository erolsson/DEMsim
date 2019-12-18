//
// Created by elahe on 2019-12-03.
//
#include "simulations.h"
#include <vector>


#include "../engine/engine.h"

#include "../contact_models/viscoelastic.h"
#include "../materials/ViscoelasticMaterial.h"

#include "../utilities/file_reading_functions.h"

void DEM::electrode_box(const std::string &settings_file_name) {
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto N = parameters.get_parameter<std::size_t>("N");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");

    //auto aspect_ratio_after_filling = parameters.get_parameter<double>("aspect_ratio_after_filling");
    EngineType simulator(1us);
    auto mat = simulator.create_material<ViscoelasticMaterial>(4.8);
    mat->E = parameters.get_parameter<double>("E");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->mu = parameters.get_parameter<double>("mu");//VAD ?
    mat->mu_wall = parameters.get_parameter<double>("mu_wall");
    mat->tau_i=parameters.get_vector<double>( "tau_i" );
    mat->alpha_i=parameters.get_vector<double>( "alpha_i" );
    mat->bt = parameters.get_parameter<double>("bt");




    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
    std::sort(particle_radii.rbegin(), particle_radii.rend());

    double particle_volume = 0.;
    for(auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }
    std::cout << "Volume of simulated particles is " << particle_volume << "\n";

    auto box_width = pow(3*particle_volume/4*pi, 1./3)*2;
    auto box_height =box_width/0.2;
    std::cout << "The simulated box has a width of " << box_width << " and a height of "
              << box_height << "\n";


    auto particle_positions = random_fill_box(0.0, box_height, box_width, particle_radii, mat->bt);
    for (std::size_t i=0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), mat);
    }
    double scale=2;

    // Creating The bottom plate surface
    Vec3 p4(-scale, -scale, 0.);
    Vec3 p2(-scale,  scale, 0.);
    Vec3 p3(scale,   -scale , 0.);
    Vec3 p1(scale,  scale, 0.);

    // Creating The top plate surface
    Vec3 p5(-scale, -scale, box_height);
    Vec3 p6(-scale,  scale, box_height);
    Vec3 p7(scale,   -scale , box_height);
    Vec3 p8(scale,  scale, box_height);


    std::vector<Vec3> bottom_points{p4, p3, p2, p1};
    std::vector<Vec3> top_points{p5, p6, p7, p8};
    std::vector<Vec3> side_points_1{p8, p6, p1, p2};
    std::vector<Vec3> side_points_2{ p6, p5,p2, p4};
    std::vector<Vec3> side_points_3{ p5, p7,p4, p3};
    std::vector<Vec3> side_points_4{p7, p8, p3, p1};


    auto bottom_surface = simulator.create_point_surface(bottom_points, true);
    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << std::endl;

    auto top_surface = simulator.create_point_surface(top_points, true);
    std::cout << "Normal of top surface is " << top_surface->get_normal() << std::endl;

    auto side_surface_1 = simulator.create_point_surface(side_points_1, true);
    std::cout << "Normal of first side surface is " << side_surface_1->get_normal() << std::endl;

    auto side_surface_2 = simulator.create_point_surface(side_points_2, true);
    std::cout << "Normal of second side surface is " << side_surface_2->get_normal() << std::endl;

    auto side_surface_3 = simulator.create_point_surface(side_points_3, true);
    std::cout << "Normal of third side surface is " << side_surface_3->get_normal() << std::endl;

    auto side_surface_4 = simulator.create_point_surface(side_points_4, true);
    std::cout << "Normal of fourth side surface is " << side_surface_4->get_normal() << std::endl;

    auto output1 = simulator.create_output(output_directory, 0.01s);
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_contacts = true;

    simulator.set_gravity(Vec3(0, 0, -9.820));
    simulator.set_mass_scale_factor(10.0e5);
    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, 0.1s);

    simulator.run(run_for_time);
    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.01s);
    simulator.run(max_velocity);
    auto bbox = simulator.get_bounding_box();
    double h = bbox[5];
    std::cout<< h<< "h surface"<< std::endl;
    bottom_surface->move(-Vec3(0, 0, box_height - h), Vec3(0, 0, 0));
}
