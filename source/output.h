//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_OUTPUT_H
#define DEMSIM_OUTPUT_H


#include <vector>
#include <string>
#include <memory>
#include <fstream>

#include "settings.h"
#include "contact_matrix.h"
#include "spherical_particle.h"
#include "contact.h"
#include "surface.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Engine;

    template<typename ForceModel, typename ParticleType>
    class Output {
        using SurfaceType = Surface<ForceModel, ParticleType>;
        using ContactPointerType = std::shared_ptr<Contact<ForceModel, ParticleType>>;

    public:
        double interval; //Frequency for the output;
        std::string name;
        bool active;
        bool activated;
        bool print_particles;
        bool print_surfaces;
        bool print_surface_forces;
        bool print_friction_forces;
        bool print_kinetic_energy;
        bool print_average_contacts;
        bool print_velocity_info;
        bool print_wall_positions;
        bool print_contacts;
        bool print_tangential_displacement;
        bool print_track_particles;
        bool print_fractured_particles;
        std::vector<ParticleType*> particles_to_track;

        Output(Engine<ForceModel, ParticleType>, double);

        void run_output(); //Prints info

    private:
        std::vector<ParticleType*>& particles_;
        std::vector<SurfaceType*>& surfaces_;
        ContactMatrix<ContactPointerType>& contacts_;
        std::vector<bool> fractured_particles_;
        const Engine<ForceModel, ParticleType>& engine_;
        Settings*& settings_;
        unsigned counter_;
        double start_time_;
        double end_time_;

        void write_particles() const;
        //void Surfaces() const;  //broken!
        void write_surface_forces() const;
        void write_friction_forces() const;
        void write_kinetic_energy() const;
        void write_average_contacts() const;
        void write_velocity_info() const;
        void write_tangential_displacement() const;
        void write_wall_positions() const;
        void write_track_particles() const;
        void write_contacts() const;
        void write_fractured_particles();
    };


}



#endif //DEMSIM_OUTPUT_H
