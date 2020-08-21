//
// Created by erolsson on 30/07/2020.
//

#include "material_base.h"
#include "../utilities/file_reading_functions.h"

std::string DEM::MaterialBase::restart_data() const {
    std::ostringstream ss;
    ss << named_print(id, "id") << ", " << named_print(density, "density");
    return ss.str();
}

DEM::MaterialBase::MaterialBase(const DEM::ParameterMap& parameters) :
    id(parameters.get_parameter<std::size_t>("id")),
    density(parameters.get_parameter<double>("density"))
{
    // Empty constructor
}

