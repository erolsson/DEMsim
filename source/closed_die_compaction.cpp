//
// Created by erolsson on 2019-01-05.
//

#include <string>

#include "engine.h"
#include "elastic_ideal_plastic_material.h"
#include "file_reading_functions.h"
#include "simulations.h"
#include "storakers_mesarovic_johnson.h"

void DEM::closed_die_compaction(const std::string& settings_file_name){
    using namespace DEM;
    using ForceModel = StorakersMesarovicJohnson;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;

    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);

    auto N = parameters.get<std::size_t>("N");
    auto output_directory = parameters.get<std::string>("output_dir");
    auto particle_file = parameters.get<std::string>("radius_file");
    auto filling_density = parameters.get<double>("filling_density");

    EngineType simulator(1us);

    auto material = simulator.create_material<ElasticIdealPlasticMaterial>(2630.);
    material->sY = parameters.get<double>("sY");
    material->E = parameters.get<double>("E");
    material->nu = parameters.get<double>("nu");

    material->mu = parameters.get<double>("mu");
    material->mu_wall = parameters.get<double>("mu_wall");
    material->kT = parameters.get<double>("kT");
}
