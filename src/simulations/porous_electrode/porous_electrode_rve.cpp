//
// Created by erolsson on 05/11/2020.
//

#include "../simulations.h"

#include "../../contact_models/porous_electrode_contact.h"
#include "../../materials/porous_electrode_material.h"
#include "../../particles/spherical_particle.h"
#include "../../engine/engine.h"
#include "../../utilities/filling_functions.h"

void DEM::porous_electrode_rve(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = PorousElectrodeContact;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto radius = parameters.get_parameter<double>("R");
    auto N = parameters.get_parameter<unsigned >("N");
    auto output_directory = parameters.get_parameter<std::string>("output_directory");
    auto initial_particle_density = parameters.get_parameter<double>("initial_particle_density");
    auto electrode_particle_density = parameters.get_parameter<double>("electrode_particle_density");
    auto compaction_time = parameters.get_parameter<double>("compaction_time");

    EngineType simulator(1us);
    auto mat = simulator.create_material<PorousElectrodeMaterial>(4800);

    mat->E_binder = parameters.get_parameter<double>("E_binder");
    mat->v_binder = parameters.get_parameter<double>("v_binder");
    mat->E_particle = parameters.get_parameter<double>("E_particle");
    mat->v_particle = parameters.get_parameter<double>("v_particle");
    mat->binder_thickness = parameters.get_parameter<double>("binder_thickness");
    mat->binder_radius_fraction =parameters.get_parameter<double>("binder_radius_fraction");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->tau_i =parameters.get_vector<double>("tau_i");

    auto particle_radii = std::vector<double>(N, radius);
    auto particle_volume = std::accumulate(particle_radii.begin(), particle_radii.end(), 0.,
                                           [](const auto& r, double vol) {
                                              return vol + 4*pi*r*r*r/3; });
    auto box_volume = particle_volume/initial_particle_density;
    auto box_side = pow(box_volume, 1./3);

    auto particle_positions = random_fill_box(-box_side/2, box_side/2, -box_side/2, box_side/2,
                                              -box_side/2, box_side/2,  particle_radii, mat->binder_thickness);


}
