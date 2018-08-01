//
// Created by erolsson on 2018-07-30.
//

#include "linear_stick_slip_model.h"
#include "linear_contact_material.h"
#include "spherical_particle.h"


DEM::LinearStickSlipModel::LinearStickSlipModel(DEM::LinearStickSlipModel::ParticleType* p1,
        DEM::LinearStickSlipModel::ParticleType* p2, double)
{
    auto mat1 = dynamic_cast<const LinearContactMaterial*>(p1->get_material());
    auto mat2 = dynamic_cast<const LinearContactMaterial*>(p2->get_material());
    double R0 = 1./(1./p1->get_radius() + 1./p2->get_radius());
    k_ = (mat1->k + mat2->k)/2*R0;
    kT_ = (mat1->kT + mat2->kT)/2*R0;
    mu_ = (mat1->mu + mat2->mu)/2;
}

DEM::LinearStickSlipModel::LinearStickSlipModel(DEM::LinearStickSlipModel::ParticleType* p1,
        DEM::LinearStickSlipModel::SurfaceType*, double)
{
    auto mat1 = dynamic_cast<const LinearContactMaterial*>(p1->get_material());
    double R0 = p1->get_radius();
    k_ = mat1->k*R0;
    kT_ = mat1->kT*R0;
    mu_ = mat1->mu_wall;
}


void DEM::LinearStickSlipModel::update(double h, const Vec3& dt)
{
    h_ = h;
    if (h > 0) {
        F_ = k_*h_;

        FT_ -= kT_*dt;
        uT_ += dt;
        if (FT_.length() > mu_*F_){ // Slip
            uT_ = mu_*F_/kT_*uT_.normal();
            FT_ = -kT_*uT_;
        }
    }
    else {
        F_ = 0;
        FT_.set_zero();
        uT_.set_zero();
    }
}
