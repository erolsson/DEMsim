//
// Created by erolsson on 2018-09-18.
//

#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include <atomic>
#include <map>
#include <string>

namespace DEM {
    using SimulationFunctionPtr = void (*)(const std::string&);

    void gyratory_compaction(const std::string& settings_file_name);
    void closed_die_compaction(const std::string& settings_file_name);
    void contact_tester(const std::string& settings_file_name);
    void cyclic_triaxial(const std::string& settings_file_name);
    void proctor_test(const std::string& settings_file_name);
    void stone_compaction(const std::string& settings_file_name);
    void electrode_box(const std::string& settings_file_name);

    std::map<std::string, SimulationFunctionPtr> valid_simulations();


}

#endif //DEMSIM_SIMULATIONS_H
