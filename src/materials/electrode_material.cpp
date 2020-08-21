//
// Created by erolsson on 12/08/2020.
//

#include "electrode_material.h"

#include "../utilities/file_reading_functions.h"
#include "../utilities/printing_functions.h"

DEM::ElectrodeMaterial::ElectrodeMaterial(const ParameterMap& parameters) :
    DEM::MaterialBase(parameters),
    E(parameters.get_parameter<double>("E")),
    nu(parameters.get_parameter<double>("nu")),
    Ep(parameters.get_parameter<double>("Ep")),
    nup(parameters.get_parameter<double>("nup")),
    fb(parameters.get_parameter<double>("fb")),
    bt(parameters.get_parameter<double>("bt")),
    yield_stress(parameters.get_parameter<double>("yield_stress")),
    tau_i(),
    alpha_i(),
    kT(parameters.get_parameter<double>("kT")),
    mu(parameters.get_parameter<double>("mu")),
    mu_wall(parameters.get_parameter<double>("mu_wall"))
{
    auto M = parameters.get_parameter<std::size_t>("M");
    for (std::size_t i = 0; i != M; ++i) {
        tau_i.push_back(parameters.get_parameter<double>("tau_" + std::to_string(i)));
        alpha_i.push_back(parameters.get_parameter<double>("alpha_" + std::to_string(i)));
    }
}

std::string DEM::ElectrodeMaterial::restart_data() const {
    std::ostringstream ss;
    ss << named_print("electrode_material", "type") << ", "
       << MaterialBase::restart_data() << ", "
       << named_print(E, "E") << ", "
       << named_print(nu, "nu") << ", "
       << named_print(fb, "fb") << ", "
       << named_print(nup, "nup") << ", "
       << named_print(Ep, "Ep") << ", "
       << named_print(yield_stress, "yield_stress") << ", "
       << named_print(bt, "bt") << ", "
       << named_print(kT, "kT") << ", "
       << named_print(mu, "mu") << ", "
       << named_print(mu_wall, "mu_wall") << ", "
       << named_print(M(), "M");

    for (std::size_t i = 0; i != M(); ++i) {
        ss << ", " << named_print(tau_i[i], "tau_"+ std::to_string(i)) << ", "
           << named_print(alpha_i[i], "alpha_"+ std::to_string(i));
    }
    return ss.str();
}

