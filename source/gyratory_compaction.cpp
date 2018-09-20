//
// Created by erolsson on 2018-09-18.
//

#include "simulations.h"

#include "engine.h"
#include "file_reading_functions.h"
#include "linear_contact_material.h"
#include "linear_stick_slip_model.h"
#include "vec3.h"


void DEM::gyratory_compaction(const std::string& settings_file_name){
    using namespace DEM;
    using ForceModel = LinearStickSlipModel;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;

    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);

    auto N = parameters.get<double>("N");
    auto output_directory = parameters.get<std::string>("output_dir");
    auto particle_file = parameters.get<std::string>("radius_file");

    EngineType simulator(1us);

    auto m = simulator.create_material<LinearContactMaterial>(1000.);
    m->k = 10;

    // Read particle radii from file
    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
    auto aspect_ratio_at_dense = 1.;
    double particle_volume = 0.;
    for(auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }

    std::cout << "Volume of simulated particles is " << particle_volume << "\n";
    auto cylinder_radius = pow(4*particle_volume/pi/aspect_ratio_at_dense, 1./3)/2;
    auto particle_1 = simulator.create_particle(1., Vec3(-1.1,0,0), Vec3(1.,0,0), m);
    simulator.create_particle(2., Vec3(1,0,0), Vec3(0,0,0), m);

    // Creating a surface
    Vec3 p1(-1, -1, 2);
    Vec3 p2(-1, 1, 2);
    Vec3 p3(1, 1, 2);
    Vec3 p4(1, -1, 2);
    std::vector<Vec3> points{p1, p2, p3, p4};

    auto surf = simulator.create_point_surface(points, true);
    auto cylinder = simulator.create_cylinder(cylinder_radius, Vec3(0, 0, 1), Vec3(1, 0, 0), 2, true, true);
    // Testing the cylinder class

    std::cout << "cylinder normal " << cylinder->get_normal(Vec3(0.5, 0.5, 0)) << std::endl;
    std::cout << "vector to cylinder " << cylinder->vector_to_point(Vec3(3, 0, 3)) << std::endl;

    std::cout << "Surface normal " << surf->get_normal() << std::endl;

    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, 0.5s);
    simulator.run(run_for_time);

    std::cout << simulator.get_time().count() << std::endl;
    std::cout << particle_1->get_position() << "\n";
}