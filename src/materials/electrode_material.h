//
// Created by elahe on 2019-11-14.
//
#include "material_base.h"

#include <vector>

#ifndef DEMSIM_VISCOELASTICMATERIAL_H
#define DEMSIM_VISCOELASTICMATERIAL_H

namespace DEM {
    class ParameterMap;
    class ElectrodeMaterial : public MaterialBase {
    public:
        ElectrodeMaterial(unsigned id_number, double density) : MaterialBase(id_number, density) {}
        explicit ElectrodeMaterial(const ParameterMap& parameters);
        ~ElectrodeMaterial() override = default;
        [[nodiscard]] std::string restart_data() const override;

        double E;    // E0 for the binder material
        double nu;   // nu for the binder material

        double Ep;   // Young's modulus for the particles
        double nup;  // Poisson's ratio for the binder material

        double fb;   // Fraction binder contacts
        double bt;   // Binder thickness

        // double active_particle_height

        double yield_displacement_coeff { 8.59e-3 };

        // double N;


        // double contact;
        std::vector<double> tau_i;
        std::vector<double> alpha_i;
        double kT;
        double mu{ 0. };
        double mu_wall{ 0. };
        double mu_binder { 0. };
        [[nodiscard]] unsigned M() const { return tau_i.size(); }

        double active_particle_height;
    };
}
#endif