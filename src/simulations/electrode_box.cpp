//
// Created by elahe on 2019-12-03.
//
#include "simulations.h"
#include <vector>


#include "../engine/engine.h"

#include "../contact_models/viscoelastic.h"
#include "../materials/ViscoelasticMaterial.h"
#include "../utilities/file_reading_functions.h"

void DEM::electrode_box(const std::string &settings_file_name) {
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto N = parameters.get_parameter<unsigned>("N");
    auto delta = parameters.get_parameter<double>("delta");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");


    EngineType simulator(1us);
    auto mat = simulator.create_material<ViscoelasticMaterial>(4800);

    mat->E = parameters.get_parameter<double>("E");
    mat->kT=parameters.get_parameter<double>("kT");
    mat->contact = parameters.get_parameter<double>("contact");
    mat->Ep= parameters.get_parameter <double> ("Ep");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->fb = parameters.get_parameter<double>("fb");
    mat->nup = parameters.get_parameter<double>("nup");
    mat->mu = parameters.get_parameter<double>("mu");//VAD ?
    mat->mu_wall = parameters.get_parameter<double>("mu_wall");
    mat->tau_i=parameters.get_vector<double>( "tau_i" );
    mat->alpha_i=parameters.get_vector<double>( "alpha_i" );
    mat->bindervolume = parameters.get_parameter<double>("bindervolume");
    mat->active_particle_height=parameters.get_parameter<double>("active_particle_height");
    mat->bt = parameters.get_parameter<double>("bt");

    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
    std::sort(particle_radii.rbegin(), particle_radii.rend());

    std::cout << "Number of particles" <<N<< "\n";
    std::cout << "delta" <<delta<< "\n";
    double just_particle_volume=0.;
    double particle_surface_area= 0.;
    for(auto& r: particle_radii) {
        particle_surface_area += 4.*pi*(r)*(r);
        just_particle_volume += 4.*pi*(r)*(r)*(r)/3.;
    }

    std::cout << "Volume of simulated particles is " <<just_particle_volume<< "\n";
    double box_width = pow(3*(just_particle_volume)/4*pi, 1./3)/0.9770;
    double box_height =box_width*2.0;
    std::cout << "The simulated box has a width of " << box_width << " and a height of "
              << box_height << "\n";



     auto particle_positions = random_fill_box(0.0, box_height, box_width, particle_radii, mat->bt);

    for (std::size_t i=0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), mat);
    }

    // Creating The bottom plate surface
    Vec3 p4(0, 0, 0.);
    Vec3 p2(0,  box_width, 0.);
    Vec3 p3(box_width,   0 , 0.);
    Vec3 p1(box_width,  box_width, 0.);

    // Creating The top plate surface
    Vec3 p5(0, 0, box_height);
    Vec3 p6(0, box_width, box_height);
    Vec3 p7(box_width,   0 , box_height);
    Vec3 p8(box_width,  box_width, box_height);


    std::vector<Vec3> bottom_points{p4, p3, p2, p1};
    std::vector<Vec3> top_points{p5, p6, p7, p8};
    std::vector<Vec3> side_points_1{p8, p6, p1, p2};
    std::vector<Vec3> side_points_2{ p6, p5,p2, p4};
    std::vector<Vec3> side_points_3{ p5, p7,p4, p3};
    std::vector<Vec3> side_points_4{p7, p8, p3, p1};


    auto bottom_surface = simulator.create_point_surface(bottom_points, true , true);
    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << std::endl;


    auto top_surface = simulator.create_point_surface(top_points, true , false);
    std::cout << "Normal of top surface is " << top_surface->get_normal() << std::endl;

    auto side_surface_1 = simulator.create_point_surface(side_points_1, true, false);
    std::cout << "Normal of first side surface is " << side_surface_1->get_normal() << std::endl;

    auto side_surface_2 = simulator.create_point_surface(side_points_2, true, false);
    std::cout << "Normal of second side surface is " << side_surface_2->get_normal() << std::endl;

    auto side_surface_3 = simulator.create_point_surface(side_points_3, true , false);
    std::cout << "Normal of third side surface is " << side_surface_3->get_normal() << std::endl;

    auto side_surface_4 = simulator.create_point_surface(side_points_4, true, false);
    std::cout << "Normal of fourth side surface is " << side_surface_4->get_normal() << std::endl;


    auto output1 = simulator.create_output(output_directory, 0.01s);
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_contacts = true;

    simulator.set_gravity(Vec3(0, 0, -9.820));
    simulator.set_mass_scale_factor(1.0);
    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, 0.1s);

    simulator.run(run_for_time);
    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.01s);
    simulator.run(max_velocity);

    // Move the lid to the uppermost particle
    std::cout<<"beginning of compaction"<< std::endl;
    auto bbox = simulator.get_bounding_box();
    double h = bbox[5];
    std::cout<<"h"<< h << std::endl;
    top_surface->move(-Vec3(0, 0, box_height - h), Vec3(0, 0, 0));
    std::cout<<"h"<< h<< std::endl;
    double surface_velocity = 0.05;
    top_surface->set_velocity(Vec3(0, 0, 0.-surface_velocity));
    std::chrono::duration<double> compaction_time {((h - mat->active_particle_height) / surface_velocity)};
    run_for_time.reset(compaction_time);
    simulator.run(run_for_time);


    std::cout<<"beginning of unloading"<< std::endl;
    top_surface->set_velocity(Vec3(0, 0, surface_velocity*1000));
    simulator.run(max_velocity);


    std::cout<<"Calculation Porosity"<< std::endl;
    bbox = simulator.get_bounding_box();
    h = bbox[5];
    std::cout<<"h:"<< h << std::endl;
    std::cout<<"box_width:"<< box_width << std::endl;
    std::cout << "Volume of simulated particles is " <<just_particle_volume << "\n";
    double Prorosity= (1-((just_particle_volume)/ (box_width*box_width*h+0.075*(box_width*box_width*h))))*100;
    std::cout<<"Prosity is:"<< Prorosity <<std::endl;
    std::cout<<"h is:"<< h <<std::endl;


    //std::cout<<"Move the lid to the uppermost particle "<< std::endl;
    //std::vector<Vec3> points_=top_surface->get_points();
    //std::cout<<"surface height5:"<< points_[1].z() <<std::endl;
    //top_surface->move(-Vec3(0,0,points_[1].z()-h),Vec3(0,0,0));
    //std::vector<Vec3> points_=top_surface->get_points();
    std::cout<<"Moving the side surface to get force-deformation "<< std::endl;
    double side_surface_velocity=0.0005;
    side_surface_2->set_velocity(Vec3(side_surface_velocity-0. , 0, 0.));
    //double delta_=points_[1].z()*1.6/100.0;
    std::chrono::duration<double> side_surface_time {((0.0240) / surface_velocity)};


    run_for_time.reset(side_surface_time);
    simulator.run(run_for_time);
    std::vector<Vec3> points_side_=side_surface_2->get_points();
    std::cout<<"side surface:"<< points_side_[1].x() <<std::endl;

    //std::cout<<"top surface"<< points_[1].z()<<std::endl;
    //std::cout<<"side surface"
    std::cout<<"Relaxation "<< std::endl;
    side_surface_2->set_velocity(Vec3(0. , 0, 0.));
    run_for_time.reset(side_surface_time*1000);
    simulator.run(run_for_time);
    std::cout<<"side surface:"<< points_side_[1].x() <<std::endl;

}

