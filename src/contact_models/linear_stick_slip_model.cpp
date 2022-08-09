//
// Created by erolsson on 2018-07-30.
//

#include <iostream>

#include "linear_stick_slip_model.h"
#include "../materials/linear_contact_material.h"


DEM::LinearStickSlipModel::LinearStickSlipModel(DEM::LinearStickSlipModel::ParticleType* p1,
        DEM::LinearStickSlipModel::ParticleType* p2, std::chrono::duration<double>):
        R0_(1./(1./p1->get_radius() + 1./p2->get_radius()))
{
    auto mat1 = dynamic_cast<const LinearContactMaterial*>(p1->get_material());
    auto mat2 = dynamic_cast<const LinearContactMaterial*>(p2->get_material());
    k_ = (mat1->k + mat2->k)/2*R0_;
    kT_ = (mat1->kT + mat2->kT)/2*R0_;
    mu_ = (mat1->mu + mat2->mu)/2;
}

DEM::LinearStickSlipModel::LinearStickSlipModel(DEM::LinearStickSlipModel::ParticleType* p1,
        DEM::LinearStickSlipModel::SurfaceType*, std::chrono::duration<double>):
        R0_(p1->get_radius())
{
    auto mat1 = dynamic_cast<const LinearContactMaterial*>(p1->get_material());
    k_ = mat1->k*R0_;
    kT_ = mat1->kT*R0_;
    mu_ = mat1->mu_wall;
}


void DEM::LinearStickSlipModel::update(double h, const Vec3& dt, const Vec3&, const Vec3& normal)
{
    h_ = h;

    if (h_ > 0) {
        F_ = k_*h_;

        // Projecting uT on the new contact plane by removing the component in the contact normal direction
        uT_ -= dot_product(uT_, normal)*normal;
        uT_ += dt;
        if (kT_*uT_.length() > mu_*F_) { // Slip
            uT_ = mu_*F_/kT_*uT_.normal();
        }
        FT_ = -kT_*uT_;
    }
    else {
        F_ = 0;
        FT_.set_zero();
        uT_.set_zero();
    }
}

std::string DEM::LinearStickSlipModel::get_output_string() const {
    std::stringstream ss;
    ss << F_;
    return ss.str();
}
