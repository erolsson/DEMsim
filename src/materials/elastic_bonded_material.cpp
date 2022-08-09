//
// Created by erolsson on 05/08/22.
//

#include "elastic_bonded_material.h"

#include "../utilities/file_reading_functions.h"
#include "../utilities/printing_functions.h"

DEM::ElasticBondedMaterial::ElasticBondedMaterial(const DEM::ParameterMap& parameters) :
    DEM::MaterialBase(parameters),
    E(parameters.get_parameter<double>("E")),
    nu(parameters.get_parameter<double>("nu")),
    kT(parameters.get_parameter<double>("kT")),
    mu(parameters.get_parameter<double>("mu")),
    mu_wall(parameters.get_parameter<double>("mu_wall")),
    k_bond(parameters.get_parameter<double>("k_bond")),
    c_bond(parameters.get_parameter<double>("c_bond")),
    bond_radius_fraction(parameters.get_parameter<double>("bond_radius_fraction")),
    fracture_stress(parameters.get_parameter<double>("fracture_stress"))
    {}

std::string DEM::ElasticBondedMaterial::restart_data() const {
    std::ostringstream ss;
    ss << named_print("elastic_bonded_material", "type") << ", "
       << MaterialBase::restart_data() << ", "
       << named_print(E, "E") << ", "
       << named_print(nu, "nu") << ", "
       << named_print(kT, "kT") << ", "
       << named_print(mu, "mu") << ", "
       << named_print(mu_wall, "mu_wall") << ", "
       << named_print(k_bond, "k_bond") << ", "
       << named_print(c_bond, "c_bond") << ", "
       << named_print(bond_radius_fraction, "bond_radius_fraction") << ", "
       << named_print(fracture_stress, "fracture_stress");
    return ss.str();
}