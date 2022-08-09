//
// Created by erolsson on 05/08/22.
//
#include "hertz_with_bonds.h"

# include "../materials/elastic_bonded_material.h"

DEM::HertzWithBonds::HertzWithBonds(SphericalParticle<HertzWithBonds>* p1, SphericalParticle<HertzWithBonds>* p2,
                                    std::chrono::duration<double>) :
        R0_(1./(1./p1->get_radius() + 1./p2->get_radius()))
{
    auto mat1 = dynamic_cast<const ElasticBondedMaterial*>(p1->get_material());
    auto mat2 = dynamic_cast<const ElasticBondedMaterial*>(p2->get_material());

    double E0 = 1/((1 - mat1->nu*mat1->nu)/mat1->E + (1 - mat2->nu*mat2->nu)/mat2->E);
    kHertz_ = 4./3*E0*sqrt(R0_);
    k_bond_ = (mat1->k_bond + mat2->k_bond)/2;
    c_bond_ = (mat1->c_bond + mat2->c_bond)/2;
    material = mat1;
    kT_ = (mat1->kT + mat2->kT)/2*R0_;
    mu_ = (mat1->mu + mat2->mu)/2;;
}

DEM::HertzWithBonds::HertzWithBonds(SphericalParticle<HertzWithBonds>* p1,
                                    Surface<HertzWithBonds, SphericalParticle<HertzWithBonds>>*,
                                    std::chrono::duration<double>) :
        R0_(p1->get_radius())

{
    auto mat1 = dynamic_cast<const ElasticBondedMaterial*>(p1->get_material());
    double E0 = mat1->E/(1-mat1->nu*mat1->nu);
    kHertz_ = 4./3*E0*sqrt(R0_);
    k_bond_ = mat1->k_bond;
    c_bond_ = mat1->c_bond;
    material = mat1;
    kT_ = 2*mat1->kT*R0_;
    mu_ = mat1->mu_wall;
}

void DEM::HertzWithBonds::update(double h, const Vec3& dt, const Vec3&, const Vec3& normal) {
    double dh = h - h_;
    h_ = h;
    if (h_ > 0) {
        F_ = kHertz_*pow(h_, 1.5);
        if (F_ < 0) {
            F_ = 0;
        }
        // Projecting uT on the new contact plane by removing the component in the contact normal direction
        uT_ -= dot_product(uT_, normal)*normal;
        uT_ += dt;
        if (kT_*uT_.length() > mu_*F_) { // Slip
            uT_ = mu_*F_/kT_*uT_.normal();
        }
        FT_ = -kT_*uT_;
    }
    else {
        h_ = 0;
        F_ = 0;
        FT_.set_zero();
        uT_.set_zero();
    }
}

std::string DEM::HertzWithBonds::get_output_string() const {
    std::stringstream ss;
    ss << F_;
    return ss.str();
}

bool DEM::HertzWithBonds::bonded() const {
    return material->bonded;
}
