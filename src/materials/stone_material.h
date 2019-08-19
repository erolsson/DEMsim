//
// Created by erolsson on 2019-04-12.
//

#ifndef DEMSIM_STONE_MATERIAL_H
#define DEMSIM_STONE_MATERIAL_H

#include "material_base.h"

namespace DEM {
    class StoneMaterial : public MaterialBase {
    public:
        StoneMaterial(unsigned id_number, double density) : MaterialBase(id_number, density) {}
        ~StoneMaterial() override = default;
        double E { 0. };
        double nu { 0. };

        // For implementing cyclic dissipation according to the results in Celma Cervera et al
        double unloading_exponent {1.5};
        double mu { 0. };
        double mu_wall { 0. };

        double hs { 20e-6 };      // Thickness of the surface layer
        double Fs { 100. };        // Force at indentation depth hs
        double min_crack_distance { 3e-3 };   // Minimum distance between cracks to treat them as separate cracks

        // Parameters for Weibull model for fracture strength
        double weibull_fracture_stress { 1e99 };
        double weibull_exponent {1. };
        double weibull_ref_volume {pow(0.00625, 3)};
    };
}

#endif //DEMSIM_STONE_MATERIAL_H
