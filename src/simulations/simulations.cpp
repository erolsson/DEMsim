//
// Created by erolsson on 2018-09-18.
//

#include <map>

#include "simulations.h"

std::map<std::string, DEM::SimulationFunctionPtr> DEM::valid_simulations() {
    // Add all valid simulation routines here
    return std::map<std::string, SimulationFunctionPtr> {
            {"gyratory_compaction",   DEM::gyratory_compaction},
            {"closed_die_compaction", DEM::closed_die_compaction}};
}

