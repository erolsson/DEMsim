//
// Created by elahe on 2019-11-14.
//

#include "viscoelastic.h"
#include <cmath>
#include "../materials/ViscoelasticMaterial.h"

#include <iostream>

DEM::Viscoelastic::Viscoelastic (DEM::Viscoelastic::ParticleType *particle1,DEM::Viscoelastic::ParticleType* particle2,
                                std::chrono::duration<double> dt)  {
    //extracting from material
    auto mat1 = dynamic_cast<const ViscoelasticMaterial *>(particle1->get_material());
    auto mat2 = dynamic_cast<const ViscoelasticMaterial *>(particle2->get_material());
    R0_ = 1. / (1. / particle1->get_radius() + 1. / particle2->get_radius());
    std::cout<<R0_<<"R";
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
        ddi_.push_back(0);
        ai.push_back(1-exp((-dt_/tau_i[i])));
        std::cout << ai[i] << std::endl;
        bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-ai[i]));
        std::cout << bi[i]<< std::endl;
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
        ai.push_back(1-exp((-dt_/tau_i[i])));
        bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-ai[i]));
    }
}
void DEM::Viscoelastic::update(double h, const DEM::Vec3& dt,const Vec3& , const DEM::Vec3& normal){
    update_normal_force(h);
    update_tangential_force(dt, normal);

}

double DEM::Viscoelastic::update_normal_force(double h)
{
    double dh=h-h_;

    if (h> bt_ && h_ > bt_){

        dF_=3./2*sqrt(h_ + bt_)*dh;

        std::cout << dF_ << std::endl;
        auto hn32 = pow(h_ + bt_, 3./2);
        auto h_32diff = pow(h + bt_, 3./ 2)-hn32;

        for (unsigned i=0 ; i != M; ++i)
        {
            ddi_ [i]= bi[i]*h_32diff + ai[i]*(hn32-di_[i]);
            dF_ -= alpha_i[i]*ddi_[i];

            di_[i]+=ddi_[i];
        }

        F_visc+=k_*dF_;


        //std::cout << F_visc << std::endl;
    }
    h_+= dh;
    return F_visc;
}
unsigned DEM::Viscoelastic::M;
//tangential

void DEM::Viscoelastic::update_tangential_force(const DEM::Vec3& dt, const DEM::Vec3& normal) {
    if (F_visc > 0.0) {
        // Projecting uT on the new contact plane by removing the component in the contact normal direction
        uT_ -= dot_product(uT_, normal)*normal;
        uT_ += dt;
        if (kT_*uT_.length() > mu_*F_visc) { // Slip
            uT_ = mu_*F_visc/kT_*uT_.normal();
        }
        FT_ = -kT_*uT_ *0.9;
        if(!dt.is_zero())
            FT_ -= mu_*F_visc*dt.normal()*0.1;
    }
    else {
        FT_.set_zero();
        uT_.set_zero();
    }


}

std::string DEM::Viscoelastic::get_output_string() const {
    std::stringstream ss;
    ss  << F_visc;
    return ss.str();
}