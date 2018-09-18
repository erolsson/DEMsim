//
// Created by erolsson on 2018-09-18.
//

#ifndef DEMSIM_SIMULATIONS_H
#define DEMSIM_SIMULATIONS_H

#include <fstream>
#include <set>
#include <string>

namespace DEM {
    void gyratory_compaction(const std::ifstream& settings_file);
    std::set<std::string> valid_simulations();
}



#endif //DEMSIM_SIMULATIONS_H
