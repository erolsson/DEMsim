//
// Created by erolsson on 02/11/2020.
//

#ifndef DEMSIM_POROUS_ELECTRODE_MATERIAL_H
#define DEMSIM_POROUS_ELECTRODE_MATERIAL_H

#include "material_base.h"

namespace DEM {
    class PorousElectrodeMaterial : public MaterialBase {
    public:
        PorousElectrodeMaterial(std::size_t id, double density) :
            MaterialBase(id, density) {}
        double E_particle = 0;
        double E_binder = 0;

        double v_particle = 0;
        double v_binder = 0;

        double binder_thickness = 0;
        double binder_radius_fraction;

        double fraction_binder_contacts = 0;

        std::vector<double> tau_i {};
        std::vector<double> alpha_i {};
    };
}

#endif //DEMSIM_POROUS_ELECTRODE_MATERIAL_H
