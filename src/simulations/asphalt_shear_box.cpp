//
// Created by erolsson on 05/08/22.
//

# include "simulations.h"
# include "../contact_models/hertz_with_bonds.h"
# include "../materials/elastic_bonded_material.h"
# include "../engine/engine.h"
# include "../utilities/filling_functions.h"
# include "../utilities/amplitude.h"

void DEM::asphalt_shear_box(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = HertzWithBonds;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);
    auto output_directory = parameters.get_parameter<std::string>("output_dir");

    EngineType simulator(1us);
    auto material = simulator.create_material<ElasticBondedMaterial>(7800);
    material->E = 200e9;
    material->nu = 0.3;
    material->mu = parameters.get_parameter<double>("mu");
    material->sY = parameters.get_parameter<double>("yield_stress");
    material->mu_wall = parameters.get_parameter<double>("mu_wall");
    material->kT = parameters.get_parameter<double>("kT");
    auto pressure = parameters.get_parameter<double>("pressure");

    double cylinder_diameter = 0.1;
    double cylinder_height = 0.035;
    double cylinder_volume = pi*cylinder_diameter*cylinder_diameter/4*cylinder_height;
    double packing_fraction = 0.625;
    double gas_density = 0.2;
    double gas_volume = cylinder_volume*packing_fraction/gas_density;
    double gas_height = 4*gas_volume/pi/cylinder_diameter/cylinder_diameter;

    double compaction_force = pressure*pi*cylinder_diameter*cylinder_diameter/4;
    std::cout << "compaction force: " << compaction_force << "\n";
    auto radius_1 = parameters.get_parameter<double>("diameter_1")/2;
    auto radius_2 = parameters.get_parameter<double>("diameter_2")/2;
    double v1 = 4*pi*radius_1*radius_1*radius_1/3;
    double v2 = 4*pi*radius_2*radius_2*radius_2/3;

    unsigned n1 = static_cast<int>(cylinder_volume/2*packing_fraction/v1);
    unsigned n2 = static_cast<int>(cylinder_volume/2*packing_fraction/v2);

    std::cout << n2 << ", " << n1 << std::endl;
    std::cout << "Gas Height: " << gas_height << "\n";

    auto p1 = Vec3(-cylinder_diameter/2, -cylinder_diameter/2, 0);
    auto p2 = Vec3( cylinder_diameter/2, -cylinder_diameter/2, 0);
    auto p3 = Vec3( cylinder_diameter/2,  cylinder_diameter/2, 0);
    auto p4 = Vec3(-cylinder_diameter/2,  cylinder_diameter/2, 0);

    auto p5 = Vec3(-cylinder_diameter/2, -cylinder_diameter/2, gas_height);
    auto p6 = Vec3( cylinder_diameter/2, -cylinder_diameter/2, gas_height);
    auto p7 = Vec3( cylinder_diameter/2,  cylinder_diameter/2, gas_height);
    auto p8 = Vec3(-cylinder_diameter/2,  cylinder_diameter/2, gas_height);

    std::vector<Vec3> bottom_points{p1, p2, p3, p4};
    std::vector<Vec3> top_points{p8, p7, p6, p5};
    auto bottom_surface = simulator.create_point_surface(bottom_points, true , false);
    auto top_surface = simulator.create_point_surface(top_points, true , false);
    std::cout << "Normal of bottom surface: " << bottom_surface->get_normal() << "\n";
    std::cout << "Normal of top surface: " << top_surface->get_normal() << "\n";

    auto bottom_cylinder = simulator.create_cylinder(cylinder_diameter/2, Vec3(0, 0, 1),
                                                Vec3(0, 0, 0), cylinder_height/2);
    auto top_cylinder = simulator.create_cylinder(cylinder_diameter/2, Vec3(0, 0, 1),
                              Vec3(0, 0, cylinder_height/2), gas_height - cylinder_height/2);

    std::vector<double> radii_1 = std::vector<double>(n1, radius_1);
    std::vector<double> radii_2 = std::vector<double>(n2, radius_2);
    auto positions_1 = random_fill_cylinder(0, gas_height, cylinder_diameter/2, radii_1);
    std::vector<ParticleType*> particles_1 {};
    for (std::size_t i = 0; i != radii_1.size(); ++i) {
        particles_1.push_back(simulator.create_particle(radii_1[i], positions_1[i], Vec3(0,0,0), material));
    }

    auto output1 = simulator.create_output(output_directory, 0.001s);
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_contacts = true;
    simulator.set_gravity(Vec3(0, 0, -9.82));
    simulator.setup(1.01);  // Needed to initialize everything

    // Running everything for 0.1 second to get things moving
    EngineType::RunForTime run_for_time(simulator, 0.1s);
    simulator.run(run_for_time);

    // Run until the fastest particle have a velocity of 0.1 and this is checked every 10 millisecond
    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.01s);
    simulator.run(max_velocity);

    auto bbox = simulator.get_bounding_box();
    auto positions_2 = random_fill_cylinder(bbox[5], gas_height, cylinder_diameter/2, radii_2);
    std::vector<ParticleType*> particles_2 {};
    for (std::size_t i = 0; i != radii_2.size(); ++i) {
        particles_2.push_back(simulator.create_particle(radii_2[i], positions_2[i], Vec3(0,0,0), material));
    }

    simulator.setup();

    // Running everything for 0.1 second to get things moving
    run_for_time.reset(0.1s);
    simulator.run(run_for_time);

    // Run until the fastest particle have a velocity of 0.1 and this is checked every 10 millisecond
    simulator.run(max_velocity);

    bbox = simulator.get_bounding_box();
    top_surface->move(Vec3(0, 0, -(gas_height - bbox[5])), Vec3());
    top_surface->set_mass(1e-3);

    auto t0 = simulator.get_time();
    auto amp_func = [compaction_force, &simulator, t0]() {
        auto time = simulator.get_time() -t0;
        if(time < 2s) {
            return -compaction_force*time/2s;
        }
        return -compaction_force;
    };

    auto amp = std::make_shared<DEM::Amplitude>(amp_func);
    top_surface->set_force_amplitude(amp, 'z');
    run_for_time.reset(2s);
    auto output2 = simulator.create_output(output_directory + "/compaction", 0.001s);
    output2->print_kinetic_energy = true;
    output2->print_surface_positions = true;
    output2->print_surface_forces = true;

    simulator.run(run_for_time);

    // The shear test
    auto top_particle = std::max_element(particles_1.begin(), particles_1.end(),
                                         [](const auto& p1, const auto& p2) {
       return p1->get_position().z() < p2->get_position().z();
    });
    double z_max = (*top_particle)->get_position().z();
    bottom_cylinder->set_length(z_max);
    top_cylinder->set_point(Vec3(0, 0, z_max));
    top_cylinder->set_velocity(Vec3(0.025/60, 0, 0));
    run_for_time.reset(24s);

    auto output3 = simulator.create_output(output_directory + "/shear_test", 0.001s);
    output3->print_kinetic_energy = true;
    output3->print_surface_positions = true;
    output3->print_surface_forces = true;

    simulator.run(run_for_time);
}