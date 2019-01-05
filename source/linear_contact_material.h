//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_MATERIAL_H
#define DEMSIM_MATERIAL_H

#include "material_base.h"

namespace DEM {
    class LinearContactMaterial : public MaterialBase {
    public:
        LinearContactMaterial(unsigned id_number, double density) : MaterialBase(id_number, density) {}
        ~LinearContactMaterial() override = default;
        double k{ 0. };
        double kT{ 0. };
        double mu{ 0. };
        double mu_wall{ 0. };
    };
}

#endif //DEMSIM_MATERIAL_H

