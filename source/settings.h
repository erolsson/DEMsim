/*
  This file is part of the DEM program DEMsim. The code is free to run and
  develop after getting permission from the original developers Erik Olsson
  and Per-Lennart Larsson, KTH Solid Mechanics, Stockholm, Sweden

  erolsson@kth.se
  pelle@kth.se
 */

#ifndef SETTINGS_H
#define SETTINGS_H

namespace DEM {
    class Settings {
    public:
        double increment;
        double end_time;
        bool mass_scaling;
        bool velocity_limitation;
        bool air_resistance;
        double max_velocity;
        double damping_time;
        double mass_scale_time;
        double mass_scale_factor;
        double friction;
        double rotation;
        double bounding_box_stretch_ratio;
        double viscosity;
        bool particle_fracture;
        double friction_time;

        Settings() :
                increment(0),
                end_time(0),
                mass_scaling(false),
                velocity_limitation(false),
                air_resistance(false),
                max_velocity(1E9),
                damping_time(end_time),
                mass_scale_time(end_time),
                mass_scale_factor(1),
                friction(true),
                rotation(true),
                bounding_box_stretch_ratio(1),
                viscosity(0.0),
                particle_fracture(false),
                friction_time(0) {} //Empty constructor
    };
}


#endif
