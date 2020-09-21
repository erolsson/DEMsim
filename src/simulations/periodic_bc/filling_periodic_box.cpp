//
// Created by erolsson on 01/09/2020.
//

#include "../simulations.h"

#include "../../contact_models/storakers_mesarovic_johnson.h"
#include "../../engine/engine.h"
#include "../../materials/elastic_ideal_plastic_material.h"
#include "../../utilities/vec3.h"
#include <fenv.h>
void DEM::filling_periodic_box(const std::string& settings_file_name) {
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

    EngineType simulator(1us);
    auto material = simulator.create_material<ElasticIdealPlasticMaterial>(2630);
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

    auto box_volume = particle_volume/filling_density;
    auto box_side = pow(box_volume, 1./3);

    auto p1 = Vec3(-box_side, -box_side, 0);
    auto p2 = Vec3( box_side, -box_side, 0);
    auto p3 = Vec3( box_side,  box_side, 0);
    auto p4 = Vec3(-box_side,  box_side, 0);
    std::vector<Vec3> bottom_points{p1, p2, p3, p4};
    auto bottom_surface = simulator.create_point_surface(bottom_points, true , true);
    auto particle_positions = random_fill_box(-box_side/2, box_side/2,
                                              -box_side/2, box_side/2,
                                              0, box_side,
                                              particle_radii);
    for (std::size_t i=0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), material);
    }

    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << std::endl;
    auto output1 = simulator.create_output(output_directory, 0.001s);
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_contacts = true;
    output1->print_periodic_bc = true;
    output1->print_mirror_particles = true;

    simulator.add_periodic_boundary_condition('x', -box_side/2, box_side/2);
    simulator.add_periodic_boundary_condition('y', -box_side/2, box_side/2);

    simulator.set_mass_scale_factor(1);
    simulator.set_gravity(Vec3(0, 0, -9.82));
    // simulator.set_rotation(false);
    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, 0.5s);
    simulator.run(run_for_time);
}

