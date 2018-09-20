//
// Created by erolsson on 2018-09-18.
//

#include "simulations.h"

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

std::set<std::string> DEM::valid_simulations() {
    // Add all valid simulation routines here
    return std::set<std::string>{"gyratory_compaction"};
}

