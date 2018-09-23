//
// Created by erolsson on 2018-09-08.
//

#include "output.h"

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

template<typename ForceModel, typename ParticleType>
DEM::Output<ForceModel, ParticleType>::Output(std::string directory, std::chrono::duration<double> interval,
                                              const Engine<ForceModel, ParticleType>& engine) :
    particles_(engine.particles_), surfaces_(engine.surfaces_), contacts_(engine.contacts_),
    directory_(directory), current_time_(engine.get_time()), time_until_output_(interval), interval_(interval)
{
    fs::create_directories(directory);
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::run_output(const std::chrono::duration<double>& increment)
{
    current_time_ += increment;
    time_until_output_ -= increment;
    if (time_until_output_ < increment) {
        time_until_output_ = interval_;
        for (const auto& func_pair: output_functions_) {
            if (func_pair.first) {    // Output-function is activated
                (this->*(func_pair.second))();
            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_particles() const
{
    std::cout << "writing particle data" << std::endl;
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_kinetic_energy() const
{
    std::cout << "writing kinetic energy" << std::endl;
}


