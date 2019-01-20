//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_OUTPUT_H
#define DEMSIM_OUTPUT_H

#include <chrono>
#include <string>
#include <vector>

#include "contact.h"
#include "../utilities/contact_matrix.h"
#include "../surfaces/surface_base.h"
#include "../particles/spherical_particle.h"

namespace DEM {

    template<typename ForceModel, typename ParticleType>
    class Engine;

    template<typename ForceModel, typename ParticleType>
    class Output {
        using SurfaceType = Surface<ForceModel, ParticleType>;
        using ContactType = Contact<ForceModel, ParticleType>;;

    public:
        Output(std::string directory, std::chrono::duration<double> interval,
               const Engine<ForceModel, ParticleType>& engine);

        void run_output(const std::chrono::duration<double>& increment); //Prints info
        bool print_particles = false;
        bool print_kinetic_energy = false;
        bool print_surface_positions = false;
        bool print_surface_forces = false;

    private:
        using OutputFunPtr = void (Output<ForceModel, ParticleType>::*)() const;
        using FuncVec = std::vector<std::pair<bool&, OutputFunPtr>>;

        const std::vector<ParticleType*>& particles_;
        const std::vector<SurfaceType*>& surfaces_;
        const ContactMatrix<ContactType>& contacts_;
        std::string directory_;

        FuncVec output_functions_ {{Output::print_particles,         &Output::write_particles},
                                   {Output::print_kinetic_energy,    &Output::write_kinetic_energy},
                                   {Output::print_surface_positions, &Output::write_surface_positions},
                                   {Output::print_surface_forces,    &Output::write_surface_forces}};

        std::chrono::duration<double> current_time_;
        std::chrono::duration<double> time_until_output_;
        std::chrono::duration<double> interval_;

        void write_particles() const;
        void write_kinetic_energy() const;
        void write_surface_positions() const;
        void write_surface_forces() const;
    };
}

#include "output.tpp"

#endif //DEMSIM_OUTPUT_H
