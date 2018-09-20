//
// Created by erolsson on 2018-09-18.
//

#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include <set>
#include <string>

namespace DEM {
    void gyratory_compaction(const std::string& settings_file_name);
    std::set<std::string> valid_simulations();
}

#endif //DEMSIM_SIMULATIONS_H
