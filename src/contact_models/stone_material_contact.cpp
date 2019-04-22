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

    mu_ = (mat1->mu + mat2->mu)/2;    
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
    
    mu_ = mat1->mu_wall;
}

void DEM::StoneMaterialContact::update(double dh, const DEM::Vec3& dt, const DEM::Vec3& normal)
{
    h_ += dh;
    double old_F = F_;
    if (h_ > 0) {
        a_ = sqrt(h_*R0_);
        if (h_ >= hmax_) {
            F_ = kp_*pow(h_, 1.5);
            hmax_ = h_;
            hp_ = (hmax_ - pow(F_/ke_, 1./1.5)); // Plastic indentation depth
            ku_ = F_/pow(h_-hp_, unloading_exp_);
        }
        else if (h_ > hp_) {
            if (h_ > h_) {
                F_ = kl_*pow(h_-hl_, 1.5);
                ku_ = F_/pow(h_-hp_, unloading_exp_);
            }
            else {
                F_ = ku_*pow(h_-hp_, unloading_exp_);
                double A = pow(F_/kp_, 1./1.5)/hmax_;
                hl_ = (h_ - A*hmax_)/(1-A);
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

    // Friction model according to simplified mindlin in
    // "An investigation of the comparative behaviour of alternative contact force models during elastic collisions"

    if (F_ > 0) {
        // Project previous contact force on the contact plane
        FT_ -= dot_product(FT_, normal)*normal;

    }
}

