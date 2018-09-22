//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_OUTPUT_H
#define DEMSIM_OUTPUT_H

#include <chrono>
#include <string>
#include <vector>

#include "contact.h"
#include "contact_matrix.h"
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
        Output(std::chrono::duration<double> interval, std::string directory,
               const Engine<ForceModel, ParticleType>& engine);

        void run_output(const std::chrono::duration<double>& increment); //Prints info
        bool print_particles = false;
        bool print_kinetic_energy = false;

    private:
        using OutputFunPtr = void (Output<ForceModel, ParticleType>::*)();
        using FuncVec = std::vector<std::pair<const bool&, OutputFunPtr>>;

        const std::vector<ParticleType*>& particles_;
        const std::vector<SurfaceType*>& surfaces_;
        const ContactMatrix<ContactPointerType>& contacts_;
        std::string directory_;

        FuncVec output_functions_ {{Output::print_particles, &Output::write_particles},
                                   {Output::print_kinetic_energy, &Output::write_kinetic_energy}};

        std::chrono::duration<double> current_time_;
        std::chrono::duration<double> time_until_output_;
        std::chrono::duration<double> interval_;

        void write_particles() const;
        void write_kinetic_energy() const;

    };
}

#include "output.tpp"

#endif //DEMSIM_OUTPUT_H
