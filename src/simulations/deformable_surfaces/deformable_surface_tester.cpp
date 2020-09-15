//
// Created by erolsson on 15/09/2020.
//

#include "../simulations.h"
#include "../../engine/engine.h"

#include "../../contact_models/viscoelastic.h"

void DEM::deformable_surface_tester(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);

    auto simulator = EngineType(1us);
    simulator.create_def
}
