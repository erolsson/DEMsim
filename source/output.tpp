//
// Created by erolsson on 2018-09-08.
//

#include "output.h"

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

template<typename ForceModel, typename ParticleType>
DEM::Output<ForceModel, ParticleType>::Output(std::chrono::duration<double> interval, std::string directory,
                                              const Engine<ForceModel, ParticleType>& engine) :
    particles_(engine.particles_), surfaces_(engine.surfaces_), contacts_(engine.contacts_),
    directory_(directory), current_time_(engine.get_time()), time_until_output_(interval), interval_(interval)
{
    fs::create_directories(directory);
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::run_output(const std::chrono::duration<double>& increment)
{

}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_particles() const
{

}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_kinetic_energy() const
{

}


