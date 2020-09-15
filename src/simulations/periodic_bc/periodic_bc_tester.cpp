//
// Created by erolsson on 25/08/2020.
//

#include <chrono>
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


void DEM::periodic_bc_tester(const std::string&) {
    using ForceModel = StorakersMesarovicJohnson;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;

    auto simulator = EngineType(1us);
    std::vector<ParticleType*> particles;

    auto mat = simulator.create_material<ElasticIdealPlasticMaterial>(7200.);
    mat->E = 200E3;
    mat->sY = 200;
    auto p = simulator.create_particle(0.1, Vec3(0, 0, 0), Vec3(0, 0, 0), mat);
    auto p2 = simulator.create_particle(0.1, Vec3(0.89, 0, 0.89), Vec3(0, 0, 0), mat);

    particles.push_back(p);
    particles.push_back(p2);
    std::vector<Surface<ForceModel, ParticleType>*> surfaces;
    ContactMatrix<Contact<ForceModel, ParticleType>> contacts;
    contacts.resize(particles.size());
    auto collision_detector = CollisionDetector<ForceModel, ParticleType>(particles, surfaces);

    auto periodic_bc_handler = PeriodicBCHandler<ForceModel, ParticleType>(simulator, particles, collision_detector,
                                                                           contacts);
    periodic_bc_handler.add_periodic_bc('x', -1., 1.);
    periodic_bc_handler.add_periodic_bc('z', -1., 1.);
    collision_detector.setup(0.01);
    collision_detector.do_check();

    for (unsigned i = 0; i != 100; ++i) {
        auto particle = particles[0];
        particle->move(Vec3(-0.05, 0., 0.05));
        std::cout << "Iteration " << i << "\n";
        std::cout << "Particle position before periodic BC: " << particle->get_position()
                  << " address " << particle <<  "\n";
        periodic_bc_handler.fulfill_periodic_bc();
        collision_detector.do_check();
        std::cout << "New collisions: " << collision_detector.contacts_to_create().size() << std::endl;
        std::cout << "Collisions to remove: " << collision_detector.contacts_to_destroy().size() << std::endl;
        periodic_bc_handler.create_periodic_bc_contacts();
        periodic_bc_handler.destroy_periodic_bc_contacts();
        for (const auto& c: contacts.get_objects()) {
            c->update();
            std::cout << "p1 position: " << c->get_particles().first->get_position() << "  "
                      << "p2 position: " << c->get_particles().second->get_position() << "\n";
            std::cout << "contact overlap: " << c->get_overlap() << "\n";
        }
        for (const auto& part: particles) {
            part->sum_contact_forces();
            std::cout << "Particle " << part->get_id() << " position after periodic BC: " << part->get_position()
                      << " Force " << part->get_force() <<  "\n";
        }
        std::cout << "Number of contacts: " << contacts.get_objects().size() << std::endl;



    }
    std::cout << "Done\n";

}