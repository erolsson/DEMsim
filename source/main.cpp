#include <iostream>
#include <vector>

#include "engine.h"
#include "linear_stick_slip_model.h"
#include "linear_contact_material.h"
#include "spherical_particle.h"
#include "vec3.h"

int main(int, char**)
{
    using namespace DEM;
    using ForceModel = LinearStickSlipModel;
    using ParticleType = SphericalParticle<ForceModel>;

    DEM::Engine<ForceModel, ParticleType> simulator;
    auto settings = simulator.get_settings();
    settings->increment = 1e-6;

    auto m = simulator.create_material<LinearContactMaterial>(1000.);
    m->k = 10;

    simulator.create_particle(1., Vec3(0,0,0), Vec3(0,0,0), m);
    simulator.create_particle(2., Vec3(0,0,0), Vec3(0,0,0), m);

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

    simulator.setup();

    return 0;
}