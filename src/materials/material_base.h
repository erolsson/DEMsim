//
// Created by erolsson on 2018-07-30.
//

#ifndef DEMSIM_MATERIAL_BASE_H
#define DEMSIM_MATERIAL_BASE_H

#include <string>
#include <sstream>
#include <memory>
#include "../utilities/printing_functions.h"

namespace DEM {
    class ParameterMap;
    class MaterialBase {
    public:
        MaterialBase(unsigned id_number, double dens) : id(id_number), density(dens) {}
        MaterialBase(const ParameterMap& parameters);
        virtual ~MaterialBase() = default;
        unsigned id ;
        double density;
        [[nodiscard]] virtual std::string restart_data() const;
    };


}

#endif //DEMSIM_MATERIAL_BASE_H
