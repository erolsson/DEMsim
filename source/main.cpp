#include <iostream>
#include <vector>

#include "linear_contact_material.h"
#include "contact_matrix.h"
#include "vec3.h"
#include "contact.h"
#include "spherical_particle.h"
#include "linear_stick_slip_model.h"
#include "point_surface.h"
#include "collision_detector.h"

int main(int, char**)
{
    using namespace DEM;
    using ForceModel = LinearStickSlipModel;
    using ParticleType = SphericalParticle<ForceModel>;
    using ContactType = std::shared_ptr<Contact<ForceModel, ParticleType>>;
    using PointSurfaceType = PointSurface<ForceModel, ParticleType>;

    DEM::ContactMatrix<ContactType> matrix = DEM::ContactMatrix<ContactType>(3);
    DEM::LinearContactMaterial m = DEM::LinearContactMaterial(0, 1000);
    m.density = 1;
    m.k = 10;

    ParticleType p1(1, Vec3(-0.9, 0, 0), Vec3(0, 0, 0), &m, 0);
    ParticleType p2(1, Vec3(0.9, 0, 0), Vec3(0, 0, 0), &m, 1);

    std::vector<ParticleType*> particles;
    std::vector<PointSurfaceType*> surfaces;

    particles.push_back(&p1);
    particles.push_back(&p2);

    CollisionDetector<ForceModel, ParticleType> collision_detector(particles, surfaces, matrix);
    collision_detector.setup();
    collision_detector.do_check();

    auto contacts_to_create = collision_detector.contacts_to_create();
    std::cout << "Contacts to create" << std::endl;
    for (const auto& contact : contacts_to_create) {
        std::cout << "\t" << contact.first->get_id() << ", " << contact.second->get_id() << std::endl;
    }

    return 0;
}