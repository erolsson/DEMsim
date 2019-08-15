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
    ke_ = E0*4./3*sqrt(R0_)/0.95;
    kl_ = E0*4./3*sqrt(R0_)/0.95;

    hs_ = (mat1->hs + mat2->hs)/2;
    Fs_ = (mat1->Fs + mat2->Fs)/2/sqrt(0.00625/R0_);
    h1_ = hs_ - pow(Fs_/kp_, 2./3);

    unloading_exp_ = (mat1->unloading_exponent + mat2->unloading_exponent)/2;
    ku_ = E0*4./3*pow(R0_, 2-unloading_exp_)/0.95;
    mu_ = (mat1->mu + mat2->mu)/2;
    old_mu_ = (mat1->mu + mat2->mu)/2;

    double G1 = E1/2/(1+v1);
    double G2 = E2/2/(1+v2);
    kT_ = 8/((2-v1)/G1 + (2-v2)/G2);
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
    ke_ = E0*4./3*sqrt(R0_)/0.95;
    kl_ = E0*4./3*sqrt(R0_)/0.95;

    hs_ = mat1->hs;
    Fs_ = mat1->Fs/sqrt(0.00625/R0_);
    h1_ = hs_ - pow(Fs_/kp_, 2./3);

    unloading_exp_ = mat1->unloading_exponent;
    ku_ = E0*4./3*pow(R0_, 2-unloading_exp_)/0.95;           //Unloading stiffness

    mu_ = mat1->mu_wall;
    old_mu_ = mat1->mu_wall;

    double G1 = E1/2/(1+v1);
    kT_ = 8/((2-v1)/G1);
}

void DEM::StoneMaterialContact::update(double h, const DEM::Vec3& dt, const Vec3& rot, const DEM::Vec3& normal)
{
    double new_FN = update_normal_force(h - h_);
    double d_mu_FN = new_FN*mu_ - FN_*old_mu_;
    FN_ = new_FN;
    update_tangential_force(dt, normal, d_mu_FN);
    update_tangential_resistance(rot);

}

double DEM::StoneMaterialContact::update_normal_force(double dh) {
    h_ += dh;
    double F = 0;
    if (h_ > 0) {
        a_ = sqrt(h_*R0_);
        if (h_ > hmax_) {
            hmax_ = h_;
            if (h_ < hs_) {
                F = Fs_/hs_*h_;
            }
            else {
                F = kp_*pow(h_ - h1_, 1.5);
            }
            hp_ = (h_ - pow(F/ke_, 1./1.5)); // Plastic indentation depth
            ku_ = F/pow(h_-hp_, unloading_exp_);
        }
        else if (h_ > 0) {
            if (dh > 0) {
                if (h_ > hl_){
                    F = kl_*pow(h_-hl_, 1.5);
                    ku_ = F/pow(h_-hp_, unloading_exp_);
                }
            }
            else if (h_ > hp_) {
                double Fmax;
                if (hmax_ > hs_) {
                    Fmax = kp_*pow(hmax_ - h1_, 1.5);
                }
                else {
                    Fmax = Fs_/hs_*hmax_;
                }
                F = ku_*pow(h_-hp_, unloading_exp_);
                double A = pow(F/Fmax, 1./1.5);
                hl_ = (h_ - A*hmax_)/(1-A);
                kl_ = Fmax/pow(hmax_-hl_, 1.5);
            }
            else {
                kl_ = ke_;
            }
        }
        else {
            a_ = 0;
            hl_ = hp_;
            kl_ = ke_;
        }
    }
    else if (h_ < -R0_){
        hmax_ = 0;
        hp_ = 0;
        F = 0;
        hl_ = 0;
        kl_ = ke_;
    }
    else {
        hmax_ = 0;
        a_ = 0;
        hl_ = hp_;
        kl_ = ke_;
        F = 0;
    }
    return F;
}

void DEM::StoneMaterialContact::update_tangential_force(const Vec3& dt, const Vec3& normal, double d_mu_FN) {
    /*
     * Friction model according to simplified Mindlin in
     * "An investigation of the comparative behaviour of alternative contact force models during elastic collisions"
     * All equation numbers refer to that paper
    */

    if (mu_*FN_ > 0) {  // Only care about friction if we have a normal force
        if (!dt.is_zero()) {
            // Project previous contact force on the contact plane
            FT_ -= dot_product(FT_, normal)*normal;
            // We have a change in loading direction, from "loading" to unloading
            if (dot_product(dt, old_dT_) < 0) {
                turning_point_ *= -1;
                FT0[(turning_point_ + 1)/2] = FT_;
            }

            // Update the turning points with changes in normal force
            for (auto &F0: FT0) {
                if (!F0.is_zero()) {
                    // Equation (23) for the 2D case F** is negative and a - sign is used
                    // This is accounted by the normal instead of F0
                    F0 += d_mu_FN*F0.normal();
                }
            }
            int sgn = 1;
            double q = 0;
            if (d_mu_FN <= 0 || load_var_tangential_.length() > load_var_normal_) {
                if (dot_product(dt, old_dT_) >= 0. && FT0[0].is_zero()) {     // Loading of the contact
                    q = 1 - (FT_.length() + d_mu_FN)/(mu_*FN_);
                    // std::cout << "loading\n";
                }
                else if (turning_point_ < 0) {                    // Unloading
                    q = 1 - ((- FT0[0] + FT_).length() + 2*d_mu_FN)/(2*mu_*FN_);
                    // std::cout << "unloading " << (FT0[0] - FT_).length() << ", " << 2*d_mu_FN << ", "
                    //          << FT0[0] << "," << FT_ << "\n";
                    sgn = -1;
                }
                else {                                                        // Reloading
                    q = 1 - ((- FT_ + FT0[1]).length() + 2*d_mu_FN)/(2*mu_*FN_);
                    // std::cout << "reloading\n";
                }
                if (d_mu_FN < 0) {
                    load_var_tangential_.set_zero();
                    load_var_normal_ = 0;
                }
            }
            else {
                q = 1;
                load_var_tangential_ += kT_*a_*dt;
                load_var_normal_ += d_mu_FN;
            }
            if (q > 0) {
                q = pow(q, 1./3);
            }
            else {
                q = 0;
            }

            double kt = kT_*a_*q + sgn*d_mu_FN*(1 - q)/dt.length();
            FT_ -= kt*dt;
            old_dT_ = dt;
            if (FT_.length() > mu_*FN_)
                FT_ = mu_*FN_*FT_.normalize();
            // std::cout << FT_ << ", " << dt << ", " << h_ << ", "  << d_mu_FN*(1 - q) << ", " << kt
            //          << ", "  << q << ", " << FN_*0.3 - FT_.length() <<  "\n";
        }
    }
    else  {
        // Reset all state variables
        FT_.set_zero();
        for (auto& F0: FT0) {
            F0.set_zero();
        }
        turning_point_ = 1;
        old_dT_.set_zero();
        load_var_tangential_.set_zero();
        load_var_normal_ = 0;
    }
}

void DEM::StoneMaterialContact::update_tangential_resistance(const DEM::Vec3 &rot) {
    rot.length();
}


