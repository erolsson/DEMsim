//
// Created by erolsson on 2019-01-20.
//

#include "stone_material_contact.h"

#include "../materials/stone_material.h"

DEM::StoneMaterialContact::StoneMaterialContact(DEM::StoneMaterialContact::ParticleType* particle1,
                                                DEM::StoneMaterialContact::ParticleType* particle2,
                                                std::chrono::duration<double>)
{
    auto mat1 = dynamic_cast<const StoneMaterial*>(particle1->get_material());
    auto mat2 = dynamic_cast<const StoneMaterial*>(particle2->get_material());

    double E1 = mat1->E;
    double E2 = mat2->E;
    double v1 = mat1->nu;
    double v2 = mat2->nu;
    double E0 = 1/((1 - v1*v1)/E1 + (1 - v2*v2)/E2);

    double R1 = particle1->get_radius();
    double R2 = particle2->get_radius();
    R0_ = 1/(1/R1 + 1/R2);

    kp_ = E0*4./3*sqrt(R0_);           //Loading stiffness
    ke_ = E0*4./3*sqrt(R0_)/0.9;
    kl_ = E0*4./3*sqrt(R0_)/0.9;
    ku_ = E0*4./3*sqrt(R0_)/0.9;           //Unloading stiffness

    unloading_exp_ = (mat1->unloading_exponent + mat2->unloading_exponent)/2;
}

DEM::StoneMaterialContact::StoneMaterialContact(DEM::StoneMaterialContact::ParticleType* particle1,
                                                DEM::StoneMaterialContact::SurfaceType*,
                                                std::chrono::duration<double>)
{
    auto mat1 = dynamic_cast<const StoneMaterial*>(particle1->get_material());

    double E1 = mat1->E;
    double v1 = mat1->nu;
    double E0 = 1/((1 - v1*v1)/E1);

    double R1 = particle1->get_radius();
    R0_ = R1;

    kp_ = E0*4./3*sqrt(R0_);           //Loading stiffness
    ke_ = E0*4./3*sqrt(R0_)/0.9;
    kl_ = E0*4./3*sqrt(R0_)/0.9;
    ku_ = E0*4./3*sqrt(R0_)/0.9;           //Unloading stiffness

    unloading_exp_ = mat1->unloading_exponent;
}

void DEM::StoneMaterialContact::update(double h, const DEM::Vec3& dt, const DEM::Vec3& normal)
{
    update_normal_force(h);
    update_tangential_force(dt, normal);
}

void DEM::StoneMaterialContact::update_normal_force(double h)
{
    if (h > 0) {
        a_ = sqrt(h*R0_);
        if (h >= hmax_) {
            F_ = kp_*pow(h, 1.5);
            hmax_ = h;
            hp_ = (hmax_ - pow(F_/ke_, 1./1.5)); // Plastic indentation depth
            ku_ = F_/pow(h-hp_, unloading_exp_);
        }
        else if (h > hp_) {
            if (h > h_) {
                F_ = kl_*pow(h-hl_, 1.5);
                ku_ = F_/pow(h-hp_, unloading_exp_);
            }
            else {
                F_ = ku_*pow(h-hp_, unloading_exp_);
                double A = pow(F_/kp_, 1./1.5)/hmax_;
                hl_ = (h - A*hmax_)/(1-A);
                kl_ = kp_*pow(hmax_/(hmax_-hl_), 1.5);
            }
        }
        else {
            F_ = 0.;
            a_ = 0;
            hl_ = hp_;
            kl_ = ke_;
        }
    }
    else {
        F_ = 0.;
        a_ = 0;
        hl_ = hp_;
        kl_ = ke_;
    }
    h_ = h;
}

void DEM::StoneMaterialContact::update_tangential_force(const DEM::Vec3& dt, const DEM::Vec3& normal)
{

}

