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
    namespace fs = std::experimental::filesystem;

    using Material = StoneMaterial;
    Material mat {0, 2370.};

    SimulationParameters parameters{settings_file_name};
    auto radius = parameters.get<double>("R");
    auto increments = parameters.get<unsigned>("N");
    auto h0 = parameters.get<double>("h0");
    auto h1 = parameters.get<double>("h1");
    auto h2 = parameters.get<double>("h2");
    auto filename=parameters.get<std::string>("output_file");
    mat.E = parameters.get<double>("E");
    mat.nu = parameters.get<double>("nu");
    mat.unloading_exponent = parameters.get<double>("unloading_exponent");

    auto p1 = FractureableSphericalParticle<ForceModel>(radius, Vec3{-radius, 0, 0}, Vec3{}, &mat, 0);
    auto p2 = FractureableSphericalParticle<ForceModel>(radius, Vec3{radius, 0, 0}, Vec3{}, &mat, 1);

    auto c = Contact<ForceModel, ParticleType>(&p2, &p1, 1s);

    p1.move(Vec3{h0, 0, 0});
    p2.move(Vec3{-h0, 0, 0});
    c.update();

    fs::path path_to_output_file {filename};
    fs::create_directories(path_to_output_file.parent_path());
    std::ofstream output_file;
    output_file.open(filename);

    for(unsigned i = 0; i != increments; ++i) {
        p1.move(Vec3{h1/increments, 0, 0});
        p2.move(Vec3{-h1/increments, 0, 0});
        c.update();
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << std::endl;
    }

    for(unsigned i = 0; i != increments; ++i) {
        p1.move(Vec3{-h2/increments, 0, 0});
        p2.move(Vec3{h2/increments, 0, 0});
        c.update();
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << std::endl;
    }

    for(unsigned j=0; j !=2 ; ++j) {
        for(unsigned i = 0; i != increments; ++i) {
            p1.move(Vec3{h2/increments, 0, 0});
            p2.move(Vec3{-h2/increments, 0, 0});
            c.update();
            output_file << c.get_overlap() << ", " << c.get_normal_force().x() << std::endl;
        }

        for(unsigned i = 0; i != increments; ++i) {
            p1.move(Vec3{-h2/increments, 0, 0});
            p2.move(Vec3{h2/increments, 0, 0});
            c.update();
            output_file << c.get_overlap() << ", " << c.get_normal_force().x() << std::endl;
        }
    }
}
