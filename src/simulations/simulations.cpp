//
// Created by erolsson on 2018-09-18.
//

#include <map>

#include "simulations.h"

std::map<std::string, DEM::SimulationFunctionPtr> DEM::valid_simulations() {
    // Add all valid simulation routines here
    return std::map<std::string, SimulationFunctionPtr>{
            {"gyratory_compaction",            DEM::gyratory_compaction},
            {"closed_die_compaction",          DEM::closed_die_compaction},
            {"contact_tester",                 DEM::contact_tester},
            {"cyclic_triaxial",                DEM::cyclic_triaxial},
            {"proctor_test",                   DEM::proctor_test},
            {"stone_compaction",               DEM::stone_compaction},
            {"electrode_box",                  DEM::electrode_box},
            {"electrode_cylinder_filling",     DEM::electrode_cylinder_filling},
            {"electrode_cylinder_compaction",  DEM::electrode_cylinder_compaction},
            {"periodic_bc_tester",             DEM::periodic_bc_tester}
            // {"binder_contact_calculator",   DEM::binder_contact_calculator}
    };
}

