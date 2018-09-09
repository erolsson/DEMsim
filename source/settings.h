//
// Created by erolsson on 2018-09-09.
//

#ifndef DEMSIM_SETTINGS_H
#define DEMSIM_SETTINGS_H

#include <chrono>

#include "vec3.h"

namespace DEM {
    struct Settings {
        std::chrono::duration<double> increment { 0. };
        Vec3 gravity { Vec3(0,0,0) };
        double mass_scale_factor { 1. };
    };
}

#endif //DEMSIM_SETTINGS_H
