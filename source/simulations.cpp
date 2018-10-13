//
// Created by erolsson on 2018-09-18.
//

#include "simulations.h"

std::set<std::string> DEM::valid_simulations() {
    // Add all valid simulation routines here
    return std::set<std::string>{"gyratory_compaction", "cylinder_tester"};
}

