//
// Created by elahe on 2019-11-14.
//

#include "viscoelastic.h"
#include <cmath>
#include <vector>

#include "../materials/ViscoelasticMaterial.h"

#include <iostream>

DEM::Viscoelastic::Viscoelastic (DEM::Viscoelastic::ParticleType *particle1,DEM::Viscoelastic::ParticleType* particle2,
                                std::chrono::duration<double> dt)  {
    //extracting from material
    auto mat1 = dynamic_cast<const ViscoelasticMaterial *>(particle1->get_material());
    auto mat2 = dynamic_cast<const ViscoelasticMaterial *>(particle2->get_material());
    R0_ = 1. / (1. / particle1->get_radius() + 1. / particle2->get_radius());

    double E1 = mat1->E;
    double E2 = mat2->E;
    double v1 = mat1->nu;
    double v2 = mat2->nu;
    M=mat1->M();
    tau_i=mat1->tau_i;
    alpha_i=mat1->alpha_i;
    kT_=mat1->kT;
    bt_= mat1->bt;
    dt_ = dt.count();  // time increment
    tsi0_ = 1. / (((1 - v1 * v1) / E1) + ((1 - v2 * v2) / E2));
    k_=4.*tsi0_*sqrt(R0_)/3; //initial contact stiffness
    for (unsigned i=0; i!=M; ++i)
    {
        di_.push_back(0);
        ddti_.emplace_back(0., 0., 0.);
        dti_.emplace_back(0., 0., 0.);
        ddi_.push_back(0);
        ai.push_back(1-exp((-dt_/tau_i[i])));
        bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-ai[i]));
    }

//Normal force
}
DEM::Viscoelastic::Viscoelastic(DEM::Viscoelastic::ParticleType *particle1, DEM::Viscoelastic::SurfaceType *,
                                std::chrono::duration<double>dt){
    auto mat1 = dynamic_cast<const ViscoelasticMaterial *>(particle1->get_material());

    R0_ = 1. / (1. / particle1->get_radius());
    double E1 = mat1->E;

    double v1 = mat1->nu;

    M=mat1->M();
    tau_i=mat1->tau_i;
    alpha_i=mat1->alpha_i;
    kT_=mat1->kT;
    mu_=mat1->mu;
    bt_= mat1->bt;
    dt_ = dt.count();  // time increment
    tsi0_ = 1. / ((1 - v1 * v1) / E1);

    k_=4.*tsi0_*sqrt(R0_)/3; //initial contact stiffness
    for (unsigned i=0; i!=M; ++i)
    {
        di_.push_back(0);
        ddi_.push_back(0);
        ddti_.emplace_back(0., 0., 0.);
        dti_.emplace_back(0., 0., 0.);
        ai.push_back(1-exp((-dt_/tau_i[i])));
        bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-ai[i]));
    }
}
void DEM::Viscoelastic::update(double h, const DEM::Vec3& dt,const Vec3& , const DEM::Vec3& normal){
    update_normal_force(h);
    update_tangential_force(dt, normal);

}
unsigned DEM::Viscoelastic::M;

double DEM::Viscoelastic::update_normal_force(double h)
{
    double dh=h-h_;
    //std::cout << "h:" << h << std::endl;
    if (h> -bt_ && h_ > -bt_ ){
        if(h>0){
            area_ = pi*R0_*(h+bt_ );
        } else{
            area_=0;
        }

        dF_=3./2*sqrt(h_+bt_ )*dh;
        //std::cout << "dF_:" << dF_ << std::endl;
        auto hn32 = pow(h_+bt_, 3./2);
        auto h_32diff = pow(h+bt_ , 3./ 2)-hn32;
        for (unsigned i=0 ; i != M; ++i)
        {
            ddi_[i]= bi[i]*h_32diff + ai[i]*(hn32-di_[i]);
            dF_ -= alpha_i[i]*ddi_[i];
            di_[i]+=ddi_[i];
        }

         F_visc+=k_*dF_;
         //std::cout << F_visc << std::endl;
    }
    else{
        F_visc=0;
    }
    h_+= dh;
    //std::cout << "h_:" << h_ << std::endl;
    //std::cout << "dh" << dh << std::endl;

    return -F_visc;
}

void DEM::Viscoelastic::update_tangential_force(const DEM::Vec3 &dt, const DEM::Vec3 &normal) {
    uT_ -= dot_product(uT_, normal) * normal;
    uT_ += dt;
   //std::cout<<uT_<<"uT"<<std::endl;
    if (-F_visc> 0.0  && !uT_.is_zero()) {
        // Projecting uT on the new contact plane by removing the component in the contact normal direction

        dFT_=-2*area_*tsi0_*dt/bt_;
        //std::cout<<dFT_<<"dFT"<<std::endl;

        for (unsigned i=0; i!=M; ++i)
        {
            ddti_[i] = bi[i]*dt + ai[i]*(uT_-dti_[i]);
            dFT_ -= alpha_i[i]*ddti_[i];
            dti_[i]+=ddti_[i];
        }
        FT_-=dFT_;

        if (FT_.length() > -0.05*F_visc) { // contact aborted
            FT_.set_zero();
            uT_.set_zero();
            for (unsigned i=0; i!=M; ++i)
            {
                dti_[i].set_zero();
            }
        }


    } else {
        FT_.set_zero();
        uT_.set_zero();
    }
    if (std::isnan(FT_.length())) {
        //std::cout<<FT_<<"FT"<<std::endl;
        std::abort();
    }
    //std::cout<<FT_<<"FT"<<std::endl;

}

std::string DEM::Viscoelastic::get_output_string() const {
    std::stringstream ss;
    ss  << F_visc;
    //ss  << FT_;
    //ss  << uT_;
    return ss.str();
}
