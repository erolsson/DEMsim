//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_OUTPUT_H
#define DEMSIM_OUTPUT_H

#include <chrono>
#include <map>
#include <string>
#include <vector>

#include "contact.h"
#include "surface_base.h"
#include "spherical_particle.h"

namespace DEM {

    template<typename ForceModel, typename ParticleType>
    class Engine;

    template<typename ForceModel, typename ParticleType>
    class Output {
        using SurfaceType = Surface<ForceModel, ParticleType>;
        using ContactPointerType = std::shared_ptr<Contact<ForceModel, ParticleType>>;

    public:
        bool print_particles;
        bool print_kinetic_energy;
        std::vector<ParticleType*> particles_to_track;

        Output(std::chrono::duration<double> interval, std::string directory,
               const Engine<ForceModel, ParticleType>& engine);

        void run_output(const std::chrono::duration<double>& increment); //Prints info
        void activate(const std::string& name);
        void deactivate(const std::string& name);

    private:
        using OutputFunPtr = void (Output<ForceModel, ParticleType>::*)();
        using FuncMap = std::map<std::string, std::pair<bool, OutputFunPtr>>;

        const std::vector<ParticleType*>& particles_;
        const std::vector<SurfaceType*>& surfaces_;
        const ContactMatrix<ContactPointerType>& contacts_;
        std::string directory_;

        FuncMap output_functions_ {{"particles",      {false, &Output::write_particles}},
                                   {"kinetic_energy", {false, &Output::write_kinetic_energy}}};

        std::chrono::duration<double> current_time_;
        std::chrono::duration<double> time_until_output_;
        std::chrono::duration<double> interval_;

        void write_particles() const;
        void write_kinetic_energy() const;

    };


}

#include "output.tpp"



#endif //DEMSIM_OUTPUT_H
