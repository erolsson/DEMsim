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
    old_mu_ = (mat1->mu + mat2->mu)/2;
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
    old_mu_ = mat1->mu_wall;
}

void DEM::StoneMaterialContact::update(double dh, const DEM::Vec3& dt, const DEM::Vec3& normal)
{
    double new_FN = update_normal_force(dh);
    double dFN = new_FN - FN_;
    FN_ = new_FN;
    update_tangential_force(dt, normal, dFN);
}

double DEM::StoneMaterialContact::update_normal_force(double dh) {
    h_ += dh;
    if (h_ > 0) {
        a_ = sqrt(h_*R0_);
        if (h_ >= hmax_) {
            hmax_ = h_;
            hp_ = (hmax_ - pow(FN_/ke_, 1./1.5)); // Plastic indentation depth
            ku_ = FN_/pow(h_-hp_, unloading_exp_);
            return kp_*pow(h_, 1.5);
        }
        else if (h_ > hp_) {
            if (h_ > h_) {
                ku_ = FN_/pow(h_-hp_, unloading_exp_);
                return kl_*pow(h_-hl_, 1.5);
            }
            else {
                double A = pow(FN_/kp_, 1./1.5)/hmax_;
                hl_ = (h_ - A*hmax_)/(1-A);
                kl_ = kp_*pow(hmax_/(hmax_-hl_), 1.5);
                return ku_*pow(h_-hp_, unloading_exp_);
            }
        }
        else {
            a_ = 0;
            hl_ = hp_;
            kl_ = ke_;
            return 0;
        }
    }
    else {
        a_ = 0;
        hl_ = hp_;
        kl_ = ke_;
        return 0;
    }
}

void DEM::StoneMaterialContact::update_tangential_force(const Vec3& dt, const Vec3& normal, double d_mu_FN) {
    /*
     * Friction model according to simplified Mindlin in
     * "An investigation of the comparative behaviour of alternative contact force models during elastic collisions"
     * All equation numbers refer to that paper
    */

    if (FN_ > 0) {  // Only care about friction if we have a normal force
        // Project previous contact force on the contact plane
        FT_ -= dot_product(FT_, normal)*normal;
        // We have a change in loading direction, from "loading" to unloading
        if (dot_product(dt, old_dT_) < 0) {
            FT0[(turning_point_ + 1)/2] = FT_;
            turning_point_ *= -1;
        }

        // Update the turning points with changes in normal force
        for (auto& F0: FT0) {
            if (!F0.is_zero()) {
                // Equation (23) for the 2D case F** is negative and a - sign is used
                // This is accounted by the normal instead of F0
                F0 += mu_*dFN*F0.normal();
            }
        }



    }
    else  {
        // Reset all state variables
        FT_.set_zero();
        for (auto& F0: FT0) {
            F0.set_zero();
        }
    }
}


