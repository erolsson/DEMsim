//
// Created by erolsson on 2018-09-08.
//

#include "output.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "engine.h"
#include "../particles/fractureable_spherical_particle.h"
#include "../utilities/printing_functions.h"

namespace fs = std::filesystem;

template<typename ForceModel, typename ParticleType>
DEM::Output<ForceModel, ParticleType>::Output(std::string directory, std::chrono::duration<double> interval,
                                              const Engine<ForceModel, ParticleType>& engine,
                                              const std::string& name,
                                              bool remove_old_files) :
    name_(name),
    particles_(engine.particles_),
    surfaces_(engine.surfaces_),
    contacts_(engine.contacts_),
    engine_(engine),
    directory_(directory),
    current_time_(engine.get_time()),
    time_until_output_(interval),
    interval_(interval)
{
        std::cout << "Creating output: " << directory << std::endl;
    if (!fs::exists(directory_)) {
        fs::create_directories(directory_);
    }
    // Remove all output files in the directory
    else if (remove_old_files){
        for (auto& iter : fs::directory_iterator(directory_)) {
            if (iter.path().extension() == ".dou") {
                fs::remove(iter.path());
            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
DEM::Output<ForceModel, ParticleType>::Output(const DEM::ParameterMap& parameters,
                                              const Engine<ForceModel, ParticleType>& engine):
        print_particles(parameters.get_parameter<bool>("print_particles")),
        print_kinetic_energy(parameters.get_parameter<bool>("print_kinetic_energy")),
        print_surface_positions(parameters.get_parameter<bool>("print_surface_positions")),
        print_surface_forces(parameters.get_parameter<bool>("print_surface_forces")),
        print_particle_cracks(parameters.get_parameter<bool>("print_particle_cracks")),
        print_contacts(parameters.get_parameter<bool>("print_contacts")),
        print_bounding_box(parameters.get_parameter<bool>("print_bounding_box")),
        name_(parameters.get_parameter<std::string>("name")),
        particles_(engine.particles_),
        surfaces_(engine.surfaces_),
        contacts_(engine.contacts_),
        engine_(engine),
        directory_(parameters.get_parameter<std::string>("directory")),
        current_time_(engine.get_time()),
        time_until_output_(parameters.get_parameter<double>("time_until_output")),
        interval_(parameters.get_parameter<double>("interval"))
{

}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::run_output() {
    // Looping over all output functions and checking if they are enabled, if so call it
    for (const auto& func_pair: output_functions_) {
        if (func_pair.first) {    // Output-function is activated
            (this->*(func_pair.second))();
        }
    }
}


template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::run_output(const std::chrono::duration<double>& increment)
{
    using namespace std::chrono_literals;
    if (time_until_output_ - increment < -1ns || current_time_ == 0s) {
        time_until_output_ = interval_;
        run_output();
    }
    current_time_ += increment;
    time_until_output_ -= increment;
}

template<typename ForceModel, typename ParticleType>
std::string DEM::Output<ForceModel, ParticleType>::restart_data() const {
    using DEM::named_print;
    std::ostringstream ss;
    ss << named_print(directory_, "directory") << ", "
       << named_print(name_, "name") << ", " << ", "
       << named_print(current_time_.count(), "current_time") << ", "
       << named_print(time_until_output_.count(), "time_until_output") << ", "
       << named_print(interval_.count(), "interval") << ", "
       << named_print(print_particles, "print_particles") << ", "
       << named_print(print_kinetic_energy, "print_kinetic_energy") << ", "
       << named_print(print_surface_positions, "print_surface_positions") << ", "
       << named_print(print_surface_forces, "print_surface_forces") << ", "
       << named_print(print_particle_cracks, "print_particle_cracks") << ", "
       << named_print(print_contacts, "print_contacts") << ", "
       << named_print(print_bounding_box, "print_bounding_box");
    return ss.str();
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
        output_file << c->get_output_string() << "\n";
    }
    output_file.close();
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_bounding_box() const {
    std::string filename = directory_ + "/bounding_box.dou";
    std::ofstream output_file;
    output_file.open(filename, std::fstream::app);
    auto bbox = engine_.get_bounding_box();
    output_file << current_time_.count()
                << ", " << bbox[0] << ", " << bbox[1] << ", " << bbox[2]
                << ", " << bbox[3] << ", " << bbox[4] << ", " << bbox[5] << "\n";
    output_file.close();
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::set_new_directory(const std::string& directory) {
    directory_ = directory;
    if (!fs::exists(directory_)) {
        fs::create_directories(directory_);
    }
}
