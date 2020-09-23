//
// Created by erolsson on 15/09/2020.
//

#include "../simulations.h"
#include "../../engine/engine.h"

#include "../../contact_models/viscoelastic.h"

void DEM::deformable_surface_tester(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto simulator = EngineType(1us);
    auto mat = simulator.create_material<ElectrodeMaterial>(4800);

    auto output_directory = parameters.get_parameter<std::string>("output_dir");
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

    std::vector<Vec3> points = {Vec3(-0.05, -0.02, 0), Vec3(0.05, -0.02, 0),
                                Vec3(0.05, 0.02, 0),   Vec3(-0.05, 0.02, 0)};

    simulator.set_gravity(Vec3(0, 0, -9.82));
    auto deformable_surface = simulator.create_deformable_point_surface(points);

    std::cout << "normal: " << deformable_surface->get_normal() << "\n";
    simulator.create_particle(0.01, Vec3(0.03, 0, 0.015), Vec3(0, 0, 0), mat);
    simulator.create_particle(0.01, Vec3(0, 0, 0.015), Vec3(0, 0, 0), mat);
    EngineType::RunForTime run_for_time(simulator, 1s);
    auto output1 = simulator.create_output(output_directory, 0.001s, "output");
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_contacts = true;
    output1->add_particle_to_follow(1);
    simulator.setup(0.005);
    simulator.run(run_for_time);

    deformable_surface->set_in_plane_strain_rates(1., 0.);
    run_for_time.reset(1.s);
    simulator.run(run_for_time);
}
