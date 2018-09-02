#include <iostream>
#include <vector>

#include "collision_detector.h"
#include "contact.h"
#include "contact_matrix.h"
#include "cylinder.h"
#include "engine.h"
#include "linear_stick_slip_model.h"
#include "linear_contact_material.h"
#include "point_surface.h"
#include "spherical_particle.h"
#include "vec3.h"

int main(int, char**)
{
    using namespace DEM;
    using ForceModel = LinearStickSlipModel;
    using ParticleType = SphericalParticle<ForceModel>;

    DEM::Engine<ForceModel, ParticleType> simulator;
    auto m = simulator.create_material<LinearContactMaterial>(1000.);

    auto particle1 = simulator.create_particle(1., Vec3(0,0,0), Vec3(0,0,0), m);
    auto particle2 = simulator.create_particle(2., Vec3(0,0,0), Vec3(0,0,0), m);
    m->k = 10;


    // Creating a surface
    Vec3 p1(-1, -1, 2);
    Vec3 p2(-1, 1, 2);
    Vec3 p3(1, 1, 2);
    Vec3 p4(1, -1, 2);
    std::vector<Vec3> points{p1, p2, p3, p4};

    auto surf = simulator.create_point_surface(points, true);
    auto cylinder = simulator.create_cylinder(1., Vec3(0, 0, 1), Vec3(1, 0, 0), 2, true);
    // Testing the cylinder class

    std::cout << "cylinder normal " << cylinder->get_normal(Vec3(0.5, 0.5, 0)) << std::endl;
    std::cout << "vector to cylinder " << cylinder->vector_to_point(Vec3(3, 0, 3)) << std::endl;

    std::cout << "Surface normal " << surf->get_normal() << std::endl;
    /*
    particles.push_back(&particle1);
    particles.push_back(&particle2);
    surfaces.push_back(&surf);
    surfaces.push_back(&cylinder);

    CollisionDetector<ForceModel, ParticleType> collision_detector(particles, surfaces, matrix);
    collision_detector.setup();

    for (unsigned i = 0; i!= 100; ++i) {
        collision_detector.do_check();

        auto contacts_to_create = collision_detector.contacts_to_create();
        auto contacts_to_destroy = collision_detector.contacts_to_destroy();

        for (const auto& contact : contacts_to_create) {

            std::cout << "Creating contact: " << contact.get_id_pair().first << ", "
                      << contact.get_id_pair().second << " Particle 0 position "
                      << particle1.get_position() << std::endl;
            if (contact.surface == nullptr) {
                matrix.create_item_inplace(contact.get_id_pair().first, contact.get_id_pair().second,
                        contact.particle1, contact.particle2, 0.);
            }
            else {
                matrix.create_item_inplace(contact.get_id_pair().first, contact.get_id_pair().second,
                        contact.particle1, contact.surface, 0.);
            }

        }

        for (const auto& contact : contacts_to_destroy) {
            std::cout << "Destroying contact: " << contact.get_id_pair().first << ", "
                      << contact.get_id_pair().second << " Particle 0 position "
                      << particle1.get_position() << std::endl;
            matrix.erase(contact.get_id_pair().first, contact.get_id_pair().second);
        }

        particle1.move(Vec3(0.05, 0., 0.0));
    }
    */
    return 0;
}