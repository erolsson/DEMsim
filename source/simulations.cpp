//
// Created by erolsson on 2018-09-18.
//

#include "simulations.h"

#include <set>
#include <string>
std::set<std::string> DEM::valid_simulations() {
    return std::set<std::string>{"gyratory_compaction"};
}
