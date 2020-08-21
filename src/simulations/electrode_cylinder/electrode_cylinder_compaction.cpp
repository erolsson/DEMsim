//
// Created by erolsson on 21/08/2020.
//

#include "../simulations.h"

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../particles/spherical_particle.h"
#include "../../materials/electrode_material.h"

void DEM::electrode_cylinder_filling(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto filling_file_name = parameters.get_parameter<std::string>("filling_file_name");

    auto simulator = EngineType(filling_file_name);
}