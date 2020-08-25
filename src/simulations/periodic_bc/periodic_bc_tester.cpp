//
// Created by erolsson on 25/08/2020.
//

#include <iostream>

#include "../simulations.h"
#include "../../engine/periodic_bc_handler.h"
#include "../../engine/collision_detection/collision_detector.h"
#include "../../utilities/contact_matrix.h"

#include "../../engine/contact.h"
#include "../../particles/spherical_particle.h"
#include "../../surfaces/surface_base.h"
#include "../../materials/elastic_ideal_plastic_material.h"
#include "../../contact_models/storakers_mesarovic_johnson.h"


void DEM::periodic_bc_tester(const std::string& settings_file_name) {
    using ForceModel = StorakersMesarovicJohnson;
    using ParticleType = SphericalParticle<ForceModel>;
    std::vector<ParticleType*> particles;
    auto mat = new ElasticIdealPlasticMaterial(0, 7200.);
    auto p = new ParticleType(0.1, Vec3(0, 0, 0), Vec3(0, 0, 0), mat, 0);
    particles.push_back(p);
    std::vector<Surface<ForceModel, ParticleType>*> surfaces;
    ContactMatrix<Contact<ForceModel, ParticleType>> contacts;
    auto collision_detector = CollisionDetector<ForceModel, ParticleType>(particles, surfaces, contacts);

    auto periodic_bc_handler = PeriodicBCHandler<ForceModel, ParticleType>(particles, collision_detector);
    periodic_bc_handler.add_periodic_bc('x', -1., 1.);
    periodic_bc_handler.add_periodic_bc('z', -1., 1.);


    for (unsigned i = 0; i != 100; ++i) {
        p->move(Vec3(-0.05, 0., 0.));
        p->move(Vec3(0.0, 0., 0.05));
        periodic_bc_handler.fulfill_periodic_bc();
        std::cout << p->get_position() << "\n";
    }

}