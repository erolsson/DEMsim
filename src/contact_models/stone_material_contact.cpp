//
// Created by erolsson on 2019-01-20.
//

#include "stone_material_contact.h"

DEM::StoneMaterialContact::StoneMaterialContact(DEM::StoneMaterialContact::ParticleType* particle1,
                                                DEM::StoneMaterialContact::ParticleType* particle2,
                                                std::chrono::duration<double>)
{

}

DEM::StoneMaterialContact::StoneMaterialContact(DEM::StoneMaterialContact::ParticleType* particle1,
                                                DEM::StoneMaterialContact::SurfaceType* surface,
                                                std::chrono::duration<double>)
{

}

void DEM::StoneMaterialContact::update(double h, const DEM::Vec3& dt, const DEM::Vec3& normal)
{
    update_normal_force(h);
    update_tangential_force(dt, normal);
}

void DEM::StoneMaterialContact::update_normal_force(double h)
{
    if (h > 0) {
        if (h >= hmax_) {
            F_ = kp_*pow(h, 1.5);
            hmax_ = h;
            hp_ = (hmax_ - pow(F_/ke_, 1./1.5)); // Plastic indentation depth
            ku_ = F_/pow(h-hp_, unloading_exp_);
        }
        a_ = sqrt(h*R0_);
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

