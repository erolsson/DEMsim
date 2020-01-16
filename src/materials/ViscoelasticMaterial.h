//
// Created by elahe on 2019-11-14.
//
#include "material_base.h"

#include <vector>

#ifndef DEMSIM_VISCOELASTICMATERIAL_H
#define DEMSIM_VISCOELASTICMATERIAL_H

namespace DEM {
    class ViscoelasticMaterial : public MaterialBase {
    public:
        ViscoelasticMaterial(unsigned id_number, double density) : MaterialBase(id_number, density) {}

        ~ViscoelasticMaterial() override = default;

        double E;
        double mu ;
        double nu;
        double nup;
        double Ep;
        double bt;
        std::vector<double> tau_i;
        std::vector<double> alpha_i;
        [[nodiscard]] unsigned M() const { return tau_i.size(); }
        double kT;
        double mu_wall { 0. };
        double unloading_exponent;

    };
}
#endif