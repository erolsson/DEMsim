//
// Created by erolsson on 2019-01-05.
//

#ifndef DEMSIM_ELASTIC_IDEAL_PLASTIC_MATERIAL_H
#define DEMSIM_ELASTIC_IDEAL_PLASTIC_MATERIAL_H

#include "material_base.h"

namespace DEM {
    class ParameterMap;
    class ElasticIdealPlasticMaterial : public MaterialBase {
    public:
        ElasticIdealPlasticMaterial(unsigned id_number, double density) : MaterialBase(id_number, density) {}
        explicit ElasticIdealPlasticMaterial(const ParameterMap& parameters);
        ~ElasticIdealPlasticMaterial() override = default;
        double E { 0. };
        double nu { 0. };
        double sY { 0. };
        double kT{ 0. };
        double mu { 0. };
        double mu_wall { 0. };
        [[nodiscard]] std::string restart_data() const override;
    };
}

#endif //DEMSIM_ELASTIC_IDEAL_PLASTIC_MATERIAL_H
