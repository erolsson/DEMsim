//
// Created by erolsson on 2019-04-13.
//

#include <chrono>

#include "simulations.h"

#include "../engine/engine.h"
#include "../particles/fractureable_spherical_particle.h"
#include "../engine/contact.h"
#include "../contact_models/stone_material_contact.h"
#include "../materials/stone_material.h"
#include "../utilities/file_reading_functions.h"


void DEM::contact_tester(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = StoneMaterialContact;
    using ParticleType = FractureableSphericalParticle<ForceModel>;
    using namespace std::chrono_literals;
    namespace fs = std::filesystem;

    using Material = StoneMaterial;
    Material mat {0, 2370.};

    SimulationParameters parameters{settings_file_name};
    auto radius = parameters.get_parameter<double>("R");
    auto increments = parameters.get_parameter<unsigned>("N");
    auto h1 = parameters.get_parameter<double>("h1");
    auto dh = parameters.get_parameter<double>("dh");
    auto dt = parameters.get_parameter<double>("dt");
    auto filename= parameters.get_parameter<std::string>("output_file");
    mat.E = parameters.get_parameter<double>("E");
    mat.nu = parameters.get_parameter<double>("nu");
    mat.unloading_exponent = parameters.get_parameter<double>("unloading_exponent");
    mat.mu = parameters.get_parameter<double>("mu");

    auto p1 = FractureableSphericalParticle<ForceModel>(radius, Vec3{-radius, 0, 0}, Vec3{}, &mat, 0);
    auto p2 = FractureableSphericalParticle<ForceModel>(radius, Vec3{radius, 0, 0}, Vec3{}, &mat, 1);

    auto c = Contact<ForceModel, ParticleType>(&p2, &p1, 1s);

    p1.move(Vec3{h1, 0, 0});
    p2.move(Vec3{-h1, 0, 0});
    c.update();
    c.update();

    fs::path path_to_output_file {filename};
    fs::create_directories(path_to_output_file.parent_path());
    std::ofstream output_file;
    output_file.open(filename);

    for(unsigned i = 0; i != increments; ++i) {
        p1.move(Vec3{0, dt/2, 0});
        p2.move(Vec3{0, -dt/2, 0});
        c.update();
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p1.get_position().y() - p2.get_position().y() << ", "
                    << c.get_tangential_force().y() << std::endl;
    }
    dt *= 0;
    for(unsigned i = 0; i != 2*increments; ++i) {
        p1.move(Vec3{-dh/2, -dt/2, 0});
        p2.move(Vec3{dh/2, dt/2, 0});
        c.update();
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p1.get_position().y() - p2.get_position().y() << ", "
                    << c.get_tangential_force().y() << std::endl;
    }

    for(unsigned i = 0; i != 2*increments; ++i) {
        p1.move(Vec3{dh/2, dt/2, 0});
        p2.move(Vec3{-dh/2, -dt/2, 0});
        c.update();
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p1.get_position().y() - p2.get_position().y() << ", "
                    << c.get_tangential_force().y() << std::endl;
    }
}
