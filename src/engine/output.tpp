//
// Created by erolsson on 2018-09-08.
//

#include "output.h"

#include <chrono>
#include <experimental/filesystem>
#include <fstream>
#include <sstream>

#include "../particles/fractureable_spherical_particle.h"

namespace fs = std::experimental::filesystem;

template<typename ForceModel, typename ParticleType>
DEM::Output<ForceModel, ParticleType>::Output(std::string directory, std::chrono::duration<double> interval,
                                              const Engine<ForceModel, ParticleType>& engine) :
    particles_(engine.particles_), surfaces_(engine.surfaces_), contacts_(engine.contacts_),
    directory_(directory), current_time_(engine.get_time()), time_until_output_(interval), interval_(interval)
{
    if (!fs::exists(directory_)) {
        fs::create_directories(directory_);
    }
    // Remove all output files in the directory
    else {
        for (auto& iter : fs::directory_iterator(directory_)) {
            if (iter.path().extension() == ".dou") {
                fs::remove(iter.path());
            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::run_output(const std::chrono::duration<double>& increment)
{
    current_time_ += increment;
    time_until_output_ -= increment;
    using namespace std::chrono_literals;
    if (time_until_output_ - increment < -1ns) {
        time_until_output_ = interval_;

        // Looping over all output functions and checking if they are enabled, if so call it
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
    std::ostringstream filename_stream;
    filename_stream << directory_ << "/" << "particles_" << current_time_.count() << ".dou";
    std::ofstream output_file;
    output_file.open(filename_stream.str());
    for (const auto& p: particles_) {
        output_file << p->get_output_string() << "\n";
    }
    output_file.close();
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_kinetic_energy() const
{
    // Sum energies of all particles
    auto transl_e = 0.;
    auto rot_e = 0.;
    for (const auto& p : particles_){
        transl_e += p->translational_energy();
        rot_e += p->rotational_energy();
    }

    std::string filename = directory_ + "/kinetic_energy.dou";
    std::ofstream output_file;
    output_file.open(filename, std::fstream::app);
    output_file << transl_e << ", " << rot_e << ", " << transl_e + rot_e << ", " << current_time_.count() << "\n";
    output_file.close();
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_surface_positions() const
{
    std::string filename = directory_ + "/surface_positions.dou";
    std::ofstream output_file;
    output_file.open(filename, std::fstream::app);
    for (auto& surface : surfaces_) {
        output_file << surface->get_output_string() << ", ";
    }
    output_file << current_time_.count() << "\n";
    output_file.close();

}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_surface_forces() const
{
    std::string filename = directory_ + "/surface_forces.dou";
    std::ofstream output_file;
    output_file.open(filename, std::fstream::app);
    for (auto& surface : surfaces_) {
        auto F = surface->get_total_force();
        output_file << "ID=" << surface->get_id() << ", " << surface->get_normal_force() << ", "
                    << F.x() << ", " << F.y() << ", " << F.z() << ", ";
    }
    output_file << current_time_.count() << "\n";
    output_file.close();
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_particle_cracks() const {
    std::string filename = directory_ + "/particle_cracks.dou";
    std::ofstream output_file;
    output_file.open(filename, std::fstream::app);
    for (const auto& particle: particles_) {
        const auto p = dynamic_cast<const FractureableSphericalParticle<ForceModel>*>(particle);
        const auto& cracks = p->get_particle_cracks();
        for (const auto& crack : cracks) {
            output_file << "ID=" << p->get_id() << ", ID_IMPACTER=" << crack.get_impacter_id() << ", "
                        << crack.get_position().x() << ", " << crack.get_position().y()
                        << ", " << crack.get_position().z() << ", " << crack.get_force() << ", "
                        << crack.get_normal().x() << ", " << crack.get_normal().y() << ", "
                        << crack.get_normal().z() << ", " << current_time_.count() << "\n";
        }
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_contacts() const {
    std::ostringstream filename_stream;
    filename_stream << directory_ << "/" << "contacts_" << current_time_.count() << ".dou";
    std::ofstream output_file;
    output_file.open(filename_stream.str());
    for (const auto& c: contacts_.get_objects()) {
        if(c->get_overlap() > 0) {
            output_file << c->get_output_string() << "\n";
        }
    }
    output_file.close();
}


