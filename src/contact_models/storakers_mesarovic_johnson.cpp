//
// Created by erolsson on 2019-01-05.
//

#include "storakers_mesarovic_johnson.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "../materials/elastic_ideal_plastic_material.h"
#include "../utilities/printing_functions.h"

DEM::StorakersMesarovicJohnson::StorakersMesarovicJohnson(DEM::StorakersMesarovicJohnson::ParticleType* particle1,
                                                          DEM::StorakersMesarovicJohnson::ParticleType* particle2,
                                                          std::chrono::duration<double>)
{
    auto mat1 = dynamic_cast<const ElasticIdealPlasticMaterial*>(particle1->get_material());
    auto mat2 = dynamic_cast<const ElasticIdealPlasticMaterial*>(particle2->get_material());
    R0_ = 1./(1./particle1->get_radius() + 1./particle2->get_radius());

    double E1 = mat1->E;
    double E2 = mat2->E;
    double v1 = mat1->nu;
    double v2 = mat2->nu;

    double E0 = 1./((1-v1*v1)/E1 + (1-v2*v2)/E2);
    double s0 = std::min(mat1->sY, mat2->sY);

    ku_ = 6*s0/E0;

    k_ = 6.*pi*1.43*R0_*s0;
    kT_ = (mat1->kT + mat2->kT)/2*R0_;
    mu_ = (mat1->mu + mat2->mu)/2;
}

DEM::StorakersMesarovicJohnson::StorakersMesarovicJohnson(DEM::StorakersMesarovicJohnson::ParticleType* particle1,
                                                          DEM::StorakersMesarovicJohnson::SurfaceType*,
                                                          std::chrono::duration<double>)
{
    auto mat1 = dynamic_cast<const ElasticIdealPlasticMaterial*>(particle1->get_material());
    R0_ = particle1->get_radius();

    double E1 = mat1->E;
    double v1 = mat1->nu;
    double E0 = E1/(1-v1*v1);
    double s0 = mat1->sY;

    ku_ = 6*s0/E0;

    k_ = 6.*pi*1.43*R0_*s0;
    kT_ = mat1->kT*R0_;
    mu_ = mat1->mu_wall;
}

DEM::StorakersMesarovicJohnson::StorakersMesarovicJohnson(DEM::StorakersMesarovicJohnson::ParticleType*,
                                                          DEM::StorakersMesarovicJohnson::ParticleType*,
                                                          std::chrono::duration<double>,
                                                          const DEM::ParameterMap& parameters) :
    h_(parameters.get_parameter<double>("h")),
    h_max_(parameters.get_parameter<double>("h_max")),
    a_(parameters.get_parameter<double>("a")),
    a_max_(parameters.get_parameter<double>("a_max")),
    k_(parameters.get_parameter<double>("k")),
    ku_(parameters.get_parameter<double>("ku")),
    hu_max_(parameters.get_parameter<double>("hu_max")),
    R0_(parameters.get_parameter<double>("R0")),
    kT_(parameters.get_parameter<double>("kT")),
    mu_(parameters.get_parameter<double>("mu")),
    F_(parameters.get_parameter<double>("F")),
    FT_(parameters.get_vec3("FT")),
    uT_(parameters.get_vec3("uT"))
{
    // Empty constructor
}

DEM::StorakersMesarovicJohnson::StorakersMesarovicJohnson(DEM::StorakersMesarovicJohnson::ParticleType*,
                                                          DEM::StorakersMesarovicJohnson::SurfaceType*,
                                                          std::chrono::duration<double>,
                                                          const DEM::ParameterMap& parameters):
    h_(parameters.get_parameter<double>("h")),
    h_max_(parameters.get_parameter<double>("h_max")),
    a_(parameters.get_parameter<double>("a")),
    a_max_(parameters.get_parameter<double>("a_max")),
    k_(parameters.get_parameter<double>("k")),
    ku_(parameters.get_parameter<double>("ku")),
    hu_max_(parameters.get_parameter<double>("hu_max")),
    R0_(parameters.get_parameter<double>("R0")),
    kT_(parameters.get_parameter<double>("kT")),
    mu_(parameters.get_parameter<double>("mu")),
    F_(parameters.get_parameter<double>("F")),
    FT_(parameters.get_vec3("FT")),
    uT_(parameters.get_vec3("uT"))
{
    // Empty constructor
}


void DEM::StorakersMesarovicJohnson::update(double dh, const DEM::Vec3& dt, const Vec3&, const DEM::Vec3& normal)
{
    update_normal_force(dh);
    update_tangential_force(dt, normal);

}

void DEM::StorakersMesarovicJohnson::update_normal_force(double h)
{
    h_ = h;
    if (h_ > h_max_) {  // Plastic loading
        F_ = k_*h_;
        a_ = sqrt(1.43*R0_*h_);
        a_max_ = a_;
        hu_max_ = ku_*a_max_;
        h_max_ = h_;
    }
    else { // Elastic unloading
        double hu = h_max_ - h_;
        if (hu < hu_max_) {
           double x = sqrt(1 - hu*hu/hu_max_/hu_max_);   // x = a / a_max_
           F_ = k_*h_max_*2/pi*(asin(x) - x*sqrt(1-x));
        } else {
            F_ = 0;
            a_ = 0;
        }
    }
}

void DEM::StorakersMesarovicJohnson::update_tangential_force(const DEM::Vec3& dt, const DEM::Vec3& normal)
{
    if (F_ != 0.) {
        // Projecting uT on the new contact plane by removing the component in the contact normal direction
        uT_ -= dot_product(uT_, normal)*normal;
        uT_ += dt;
        if (kT_*uT_.length() > mu_*F_) { // Slip
            uT_ = mu_*F_/kT_*uT_.normal();
        }
        FT_ = -kT_*uT_ *0.9;
        if(!dt.is_zero())
            FT_ -= mu_*F_*dt.normal()*0.1;
    }
    else {
        FT_.set_zero();
        uT_.set_zero();
    }
}

std::string DEM::StorakersMesarovicJohnson::get_output_string() const {
    std::stringstream ss;
    ss << h_max_ << ", " << F_;
    return ss.str();
}

std::string DEM::StorakersMesarovicJohnson::restart_data() const {
    std::ostringstream ss;
    ss << named_print(h_, "h") << ", " << named_print(h_max_, "h_max") << ", "
       << named_print(a_, "a") << ", " << named_print(a_max_, "a_max") << ", "
       << named_print(k_, "k") << ", " << named_print(ku_, "ku") << ", "
       << named_print(hu_max_, "hu_max") << ", " << named_print(R0_, "R0") << ", "
       << named_print(kT_, "kT") << ", " << named_print(mu_, "mu") << ", "
       << named_print(F_, "F") << ", " << named_print(FT_, "FT") << ", "
       << named_print(uT_, "uT");
    return ss.str();
}

