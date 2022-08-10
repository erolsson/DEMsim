//
// Created by erolsson on 05/08/22.
//

#ifndef DEMSIM_ELASTIC_BONDED_MATERIAL_H
#define DEMSIM_ELASTIC_BONDED_MATERIAL_H

#include "material_base.h"

namespace DEM {
    class ParameterMap;
    class ElasticBondedMaterial : public MaterialBase {
    public:
        ElasticBondedMaterial(unsigned id_number, double density) : MaterialBase(id_number, density) {}
        explicit ElasticBondedMaterial(const ParameterMap& parameters);
        ~ElasticBondedMaterial() override = default;
        double E { 0. };
        double nu { 0. };
        double sY { 1e20 };
        double kT{ 0. };
        double mu { 0. };
        double mu_wall { 0. };
        double k_bond { 0. };
        double c_bond { 0. };
        double bond_radius_fraction { 0.};
        double fracture_stress { 0. };
        bool bonded  { false };
        [[nodiscard]] std::string restart_data() const override;
    };
}


#endif //DEMSIM_ELASTIC_BONDED_MATERIAL_H
