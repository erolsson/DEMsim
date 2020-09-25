//
// Created by erolsson on 2018-09-08.
//

#include "output.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <utility>

#include "Eigen/Dense"

#include "engine.h"
#include "../particles/fractureable_spherical_particle.h"
#include "../utilities/printing_functions.h"

namespace fs = std::filesystem;

template<typename ForceModel, typename ParticleType>
DEM::Output<ForceModel, ParticleType>::Output(const std::string& directory, std::chrono::duration<double> interval,
                                              const Engine<ForceModel, ParticleType>& engine,
                                              std::string  name,
                                              bool remove_old_files) :
    name_(std::move(name)),
    particles_(engine.particles_),
    surfaces_(engine.surfaces_),
    contacts_(engine.contacts_),
    engine_(engine),
    current_time_(engine.get_time()),
    time_until_output_(interval),
    interval_(interval)
{
    std::string directory_name = std::regex_replace(directory, std::regex("//"), "/");
    directory_ = (directory_name);
    std::cout << "Creating output: " << directory << std::endl;
    if (!fs::exists(directory_)) {
        fs::create_directories(directory_);
    }
    // Remove all output files in the directory
    else if (remove_old_files) {
        for (auto& iter : fs::directory_iterator(directory_)) {
            if (iter.path().extension() == ".dou") {
                fs::remove(iter.path());
            }
        }
        fs::remove_all(directory_ / fs::path("particles"));
        fs::remove_all(directory_ / fs::path("contacts/"));
        fs::remove_all(directory_ / fs::path("mirror_particles/"));
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
        print_periodic_bc(parameters.get_parameter<bool>("print_periodic_bc")),
        print_mirror_particles(parameters.get_parameter<bool>("print_mirror_particles")),
        print_fabric_force_tensor(parameters.get_parameter<bool>("print_fabric_force_tensor")),
        name_(parameters.get_parameter<std::string>("name")),
        particles_(engine.particles_),
        surfaces_(engine.surfaces_),
        contacts_(engine.contacts_),
        engine_(engine),
        directory_(fs::path(parameters.get_parameter<std::string>("directory"))),
        current_time_(engine.get_time()),
        time_until_output_(parameters.get_parameter<double>("time_until_output")),
        interval_(parameters.get_parameter<double>("interval"))
{
    auto no_particles_to_print = parameters.get_parameter<std::size_t>("no_particles_to_print");
    for (std::size_t i = 0; i != no_particles_to_print; ++i) {
        add_particle_to_follow(parameters.get_parameter<std::size_t>("p_to_print_" + std::to_string(i)));
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::run_output() {
    // Looping over all output functions and checking if they are enabled, if so call it
    for (const auto& func_pair: output_functions_) {
        if (func_pair.first) {    // Output-function is activated
            (this->*(func_pair.second))();
        }
    }
    write_particles_to_follow();
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
       << named_print(print_bounding_box, "print_bounding_box") << ", "
       << named_print(print_periodic_bc, "print_periodic_bc") << ", "
       << named_print(print_mirror_particles, "print_mirror_particles") << ", "
       << named_print(print_fabric_force_tensor, "print_fabric_force_tensor") << ", "
       << named_print(particles_to_print_.size(), "no_particles_to_print");
    for (std::size_t i = 0; i != particles_to_print_.size(); ++i) {
        ss << ", " << named_print(particles_to_print_[i], "p_to_print_" + std::to_string(i));
    }
    return ss.str();
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::add_particle_to_follow(std::size_t particle_id) {
    auto p_it = std::find_if(particles_.begin(), particles_.end(),
                          [particle_id](const auto& p) {return p->get_id() == particle_id; });
    if (p_it != particles_.end()) {
        particles_to_print_.push_back(*p_it);
    }
    else {
        throw std::invalid_argument("particle " + std::to_string(particle_id) + " does not exist");
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_particles() const
{
    std::ostringstream filename_stream;
    fs::path particle_directory = directory_ / fs::path("particles");
    if (!fs::exists(particle_directory)) {
        fs::create_directories(particle_directory);
    }
    filename_stream << "particles_" << current_time_.count() << ".dou";
    std::ofstream output_file;
    output_file.open(particle_directory / fs::path(filename_stream.str()));
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

    fs::path filename = directory_ / fs::path("kinetic_energy.dou");
    std::ofstream output_file;
    output_file.open(filename, std::fstream::app);
    output_file << transl_e << ", " << rot_e << ", " << transl_e + rot_e << ", " << current_time_.count() << "\n";
    output_file.close();
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_surface_positions() const
{
    fs::path filename = directory_ / fs::path("surface_positions.dou");
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
    fs::path filename = directory_ / "surface_forces.dou";
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
    fs::path filename = directory_ / "/particle_cracks.dou";
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
    fs::path contacts_directory = directory_ / "contacts";
    if (!fs::exists(contacts_directory)) {
        fs::create_directories(contacts_directory);
    }
    filename_stream << "contacts_" << current_time_.count() << ".dou";
    std::ofstream output_file;
    output_file.open(contacts_directory / filename_stream.str());
    for (const auto& c: contacts_.get_objects_sorted()) {
        output_file << c->get_output_string() << "\n";
    }
    output_file.close();
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_bounding_box() const {
    fs::path filename = directory_ / "bounding_box.dou";
    std::ofstream output_file;
    output_file.open(filename, std::fstream::app);
    auto bbox = engine_.get_bounding_box();
    output_file << current_time_.count()
                << ", " << bbox[0] << ", " << bbox[1] << ", " << bbox[2]
                << ", " << bbox[3] << ", " << bbox[4] << ", " << bbox[5] << "\n";
    output_file.close();
}


template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_periodic_bc() const {
    if (engine_.periodic_bc_handler_ != nullptr) {
        fs::path filename = directory_ / "periodic_bc.dou";
        std::ofstream output_file;
        output_file.open(filename, std::fstream::app);
        output_file << current_time_.count() << ", " << engine_.periodic_bc_handler_->print_periodic_bc() << "\n";
        output_file.close();
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_mirror_particles() const {
    if (engine_.periodic_bc_handler_ != nullptr) {
        auto mirror_particles = engine_.periodic_bc_handler_->mirror_particles_output();
        fs::path mirror_particles_directory = directory_ / "mirror_particles/";
        if (!fs::exists(mirror_particles_directory)) {
            fs::create_directories(mirror_particles_directory);
        }
        std::ostringstream filename_stream;
        filename_stream << "mirror_particles_" << current_time_.count() << ".dou";

        std::ofstream output_file;
        output_file.open(mirror_particles_directory / filename_stream.str());
        for (const auto& mp: mirror_particles) {
            output_file << mp << "\n";
        }
        output_file.close();
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_fabric_force_tensor() const {
    Eigen::Matrix<double, 3, 3> force_tensor = Eigen::Matrix<double, 3, 3>::Zero();
    for (const auto& c: contacts_.get_objects()) {
        force_tensor += c->get_force_fabric_tensor();
    }
    fs::path filename = directory_ / "force_fabric_tensor.dou";
    std::ofstream output_file;
    output_file.open(filename, std::fstream::app);
    output_file << current_time_.count();
    for (unsigned i = 0; i != 3; ++i) {
        for (unsigned j = 0; j != 3; ++j) {
           output_file << ", " << force_tensor(i, j);
        }
    }
    output_file << "\n";
    output_file.close();

}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::set_new_directory(const std::string& directory) {
    directory_ = directory;
    if (!fs::exists(directory_)) {
        fs::create_directories(directory_);
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Output<ForceModel, ParticleType>::write_particles_to_follow() const {
    for (const auto p: particles_to_print_) {
        fs::path filename = directory_ / fs::path(std::string("particle_" + std::to_string(p->get_id()) + ".dou"));
        std::ofstream output_file;
        output_file.open(filename, std::fstream::app);
        output_file << current_time_.count() << ", " <<  p->get_output_string() << "\n";
        output_file.close();
    }
}



