//
// Created by erolsson on 30/07/2020.
//

#include "elastic_ideal_plastic_material.h"

#include <sstream>
#include <string>

#include "../utilities/file_reading_functions.h"

DEM::ElasticIdealPlasticMaterial::ElasticIdealPlasticMaterial(const ParameterMap& parameters)
    : MaterialBase(parameters),
      E(parameters.get_parameter<double>("E")),
      nu(parameters.get_parameter<double>("nu")),
      sY(parameters.get_parameter<double>("sY")),
      kT(parameters.get_parameter<double>("kT")),
      mu(parameters.get_parameter<double>("mu")),
      mu_wall(parameters.get_parameter<double>("mu_wall")){
    //  Empty constructor body
}

std::string DEM::ElasticIdealPlasticMaterial::restart_data() const {
    std::ostringstream ss;
    ss << named_print("elastic_ideal_plastic_material", "type") << ", "
       << MaterialBase::restart_data() << ", "
       << named_print(E, "E") << ", "
       << named_print(nu, "nu") << ", "
       << named_print(sY, "sY") << ", "
       << named_print(kT, "kT") << ", "
       << named_print(mu, "mu") << ", "
       << named_print(mu_wall, "mu_wall");
    return ss.str();
}

