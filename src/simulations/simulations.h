//
// Created by erolsson on 2018-09-18.
//

#ifndef SIMULATIONS_H
#define SIMULATIONS_H

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
    //void electrode_cylinder_filling(const std::string& settings_file_name);
    void electrode_cylinder_compaction(const std::string& settings_file_name);
    void periodic_bc_tester(const std::string& settings_file_name);
    void periodic_bc_simulation(const std::string& settings_file_name);
    void filling_periodic_box(const std::string& settings_file_name);
    //void deformable_surface_tester(const std::string& settings_file_name);
    //void battery_rve_filling(const std::string& settings_file_name);
    void battery_rve_compaction(const std::string &settings_file_name);
    void battery_rve_loading(const std::string &settings_file_name);
    void Cathode_mechanical_simulations(const std::string& settings_file_name);
    void restart_electrode(const std::string& settings_file_name);
    void porous_electrode_rve(const std::string& settings_file_name);

    std::map<std::string, SimulationFunctionPtr> valid_simulations();


}

#endif //DEMSIM_SIMULATIONS_H
