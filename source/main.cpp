#include <iostream>

#include "linear_contact_material.h"
#include "contact_matrix.h"
#include "vec3.h"
#include "contact.h"
#include "spherical_particle.h"
#include "linear_stick_slip_model.h"

int main(int argc, char** argv)
{
    using namespace DEM;
    using ParticleType = SphericalParticle<LinearStickSlipModel>;
    using ContactType = Contact<LinearStickSlipModel, ParticleType>;

    DEM::ContactMatrix<ContactType> matrix = DEM::ContactMatrix<ContactType>(100);
    DEM::LinearContactMaterial m = DEM::LinearContactMaterial(0, 1000);
    m.density = 1;
    m.k = 10;

    ParticleType p1(1, Vec3(-0.9, 0, 0), Vec3(0, 0, 0), &m, 0);
    ParticleType p2(1, Vec3(0.9, 0, 0), Vec3(0, 0, 0), &m, 1);

    DEM::Contact<DEM::LinearStickSlipModel, ParticleType> c(&p1, &p2, 0.);
    matrix.insert(0, 1, c);

    c.update();
    std::cout << c.get_normal_force() << " " << c.get_overlap() << std::endl;
    std::cout << p1.get_position() << std::endl;

    return 0;
}