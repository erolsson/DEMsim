//
// Created by erolsson on 05/08/22.
//
#include "hertz_with_bonds.h"

# include "../materials/elastic_bonded_material.h"

DEM::HertzWithBonds::HertzWithBonds(SphericalParticle<HertzWithBonds>* p1, SphericalParticle<HertzWithBonds>* p2,
                                    std::chrono::duration<double> increment) :
        R0_(1./(1./p1->get_radius() + 1./p2->get_radius())),
        increment_(increment.count())

{
    auto mat1 = dynamic_cast<const ElasticBondedMaterial*>(p1->get_material());
    auto mat2 = dynamic_cast<const ElasticBondedMaterial*>(p2->get_material());

    double E0 = 1/((1 - mat1->nu*mat1->nu)/mat1->E + (1 - mat2->nu*mat2->nu)/mat2->E);
    kHertz_ = 4./3*E0*sqrt(R0_);
    double bond_radius = (mat1->bond_radius_fraction + mat2->bond_radius_fraction)*(p1->get_radius() + p2->get_radius())/4;
    double h1 = p1->get_radius() - sqrt(p1->get_radius()*p1->get_radius() - bond_radius*bond_radius);
    double h2 = p2->get_radius() - sqrt(p2->get_radius()*p2->get_radius() - bond_radius*bond_radius);
    double bond_height = h1 + h2;
    bond_area_ = bond_radius*bond_radius*pi;
    k_bond_ = (mat1->k_bond + mat2->k_bond)/2*bond_area_/bond_height;
    c_bond_ = (mat1->c_bond + mat2->c_bond)/2*bond_area_/bond_height;
    material1 = mat1;
    material2 = mat2;
    kT_ = (mat1->kT + mat2->kT)/2*R0_;
    mu_ = (mat1->mu + mat2->mu)/2;
    ky_ = 6*pi*1.43*mat1->sY*R0_;
    hy_ = pow(3*pi*1.43*mat1->sY/E0,2)*R0_;
    fracture_stress_ = (mat1->fracture_stress + mat2->fracture_stress);
}

DEM::HertzWithBonds::HertzWithBonds(SphericalParticle<HertzWithBonds>* p1,
                                    Surface<HertzWithBonds, SphericalParticle<HertzWithBonds>>*,
                                    std::chrono::duration<double> increment) :
        k_bond_(0),
        c_bond_(0),
        R0_(p1->get_radius()),
        increment_(increment.count()),
        bond_area_(0.)

{
    auto mat1 = dynamic_cast<const ElasticBondedMaterial*>(p1->get_material());
    double E0 = mat1->E/(1-mat1->nu*mat1->nu);
    kHertz_ = 4./3*E0*sqrt(R0_);
    material1 = mat1;
    material2 = mat1;
    kT_ = 2*mat1->kT*R0_;
    mu_ = mat1->mu_wall;
    ky_ = 6*pi*1.43*mat1->sY*R0_;
    hy_ = pow(3*pi*1.43*mat1->sY/E0, 2)*R0_;
    double bond_radius = mat1->bond_radius_fraction*R0_;
    bond_area_ = bond_radius*bond_radius*pi;
    double bond_height = p1->get_radius() - sqrt(p1->get_radius()*p1->get_radius() - bond_radius*bond_radius);
    k_bond_ = mat1->k_bond*bond_area_/bond_height;
    c_bond_ = mat1->c_bond*bond_area_/bond_height;
    fracture_stress_ = mat1->fracture_stress;
}

void DEM::HertzWithBonds::update(double h, const Vec3& dt, const Vec3&, const Vec3& normal) {
    double dh = h - h_;
    h_ = h;
    if (h_ > 0) {
        if (!fractured_) {
           bonded_ = true;
        }
        if (h < hy_) {
            F_ = kHertz_*pow(h_, 1.5);
        }
        else {
            double Fy = kHertz_*pow(hy_, 1.5);
            F_ = ky_*(h - hy_) + Fy;
        }
        if (F_ < 0) {
            F_ = 0;
        }
        // Projecting uT on the new contact plane by removing the component in the contact normal direction
        uT_ -= dot_product(uT_, normal)*normal;
        uT_ += dt;
        if (kT_*uT_.length() > mu_*F_) { // Slip
            FT_ = -mu_*F_*uT_.normal();
            uT_ = mu_*F_*uT_.normal()/kT_;
        }
        else {
            FT_ = -kT_*uT_;
        }
    }

    else {
        F_ = 0;
        FT_.set_zero();
        uT_.set_zero();
    }
    if (bonded() && k_bond_ > 0) {
        F_bond_ += k_bond_*dh - k_bond_/c_bond_*F_bond_*increment_;
        FT_bond_ -= (k_bond_*dt - k_bond_/c_bond_*FT_bond_*increment_);
        double effective_contact_stress = (-0.5*F_bond_
                + sqrt(F_bond_*F_bond_/4 + FT_bond_.length()*FT_bond_.length()))/bond_area_;
        if (effective_contact_stress > fracture_stress_) {
            fractured_ = true;
        }
        F_ += F_bond_;
        FT_ += FT_bond_;
    }
}

std::string DEM::HertzWithBonds::get_output_string() const {
    std::stringstream ss;
    ss << F_ << ", " << FT_.x() << ", " << FT_.y() << ", " << FT_.z() << ", " << FT_.length() << ", "
       << F_bond_ << ", " << FT_bond_.x() << ", " << FT_bond_.y() << ", " << FT_bond_.z() << ", " << FT_bond_.length();
    return ss.str();
}

bool DEM::HertzWithBonds::bonded() const {
    return material1->bonded && material2->bonded && !fractured_ && bonded_ && k_bond_ > 0;
}

