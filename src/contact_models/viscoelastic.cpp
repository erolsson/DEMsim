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
    double random = rand() % 10 + 1;


    //std::cout << "random:" << random << std::endl;
    if (random < mat1->contact){
        procent_=true;
    }
    double E1 = mat1->E;
    double E2 = mat2->E;
    double v1 = mat1->nu;
    double v2 = mat2->nu;
    double vp1 = mat1->nup;
    double vp2 = mat2->nup;
    double Ep2 = mat2->Ep;
    //id2_= particle1->get_id();
    //std::cout << "id:" << id2_ << std::endl;


    double Ep1 = mat1->Ep;
    //binder_radii_= (particle1->get_radius()*particle2->get_radius())/(particle1->get_radius()+particle2->get_radius());
    mu_ = (mat1->mu + mat2->mu)/2;
    adhesive_=true;
    M=mat1->M();
    tau_i=mat1->tau_i;
    alpha_i=mat1->alpha_i;
    kT_=mat1->kT;
    bt_= mat1->bt;
    dt_ = dt.count();  // time increment
    tsi0_ = 1. / (((1 - v1 * v1) / E1) + ((1 - v2 * v2) / E2));
    tsi0particle_=1./(((1-vp1*vp1)/Ep1)+((1-vp2*vp2)/Ep2));
    k_=4.*tsi0_*sqrt(R0_)/3; //initial contact stiffness
    kparticle_=4*tsi0particle_*sqrt(R0_)/3;
    yield_h_ = pow(std::min(mat1->yield_stress, mat2->yield_stress)*R0_*R0_/kparticle_, 2./3);
    //std::cout << "procent:" << procent_ << std::endl;

    double G1 = E1/2/(1+v1);
    double G2 = E2/2/(1+v2);
    kT_ = 8/((2-v1)/G1 + (2-v2)/G2)*0.001*R0_;

    for (unsigned i=0; i!=M; ++i)
    {
        di_.push_back(0);
        ddti_.emplace_back(0., 0., 0.);
        dti_.emplace_back(0., 0., 0.);
        ddi_.push_back(0);
        ai.push_back(1-exp((-dt_/tau_i[i])));
        //std::cout << "ai:" << ai[i] << std::endl;
        bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-ai[i]));
    }


}
DEM::Viscoelastic::Viscoelastic(DEM::Viscoelastic::ParticleType *particle1, DEM::Viscoelastic::SurfaceType * surface,
                                std::chrono::duration<double>dt){
    auto mat1 = dynamic_cast<const ViscoelasticMaterial *>(particle1->get_material());

    R0_ = 1. / (1. / particle1->get_radius());
    double E1 = mat1->E;
    double v1 = mat1->nu;
    double vp1=mat1->nup;
    double Ep1 = mat1->Ep;
    mu_ = mat1->mu_wall;


    adhesive_ = surface->adhesive();

    M=mat1->M();
    tau_i=mat1->tau_i;
    alpha_i=mat1->alpha_i;
    kT_=mat1->kT;
    bt_= mat1->bt;
    double random = rand() % 10 + 1;
    //std::cout << "bt model:" <<bt_ << std::endl;
    if (random < mat1->contact){
        procent_=true;
    }


    dt_ = dt.count();  // time increment
    tsi0_ = 1. / ((1 - v1 * v1) / E1);
    //std::cout << "tsi0:" << tsi0_ << std::endl;
    tsi0particle_=1./((1-vp1*vp1)/Ep1);
    k_=4.*tsi0_*sqrt(R0_)/3; //initial contact stiffness
    //std::cout << "k:" << k_ << std::endl;
    kparticle_=4*tsi0particle_*sqrt(R0_)/3;
    yield_h_ = pow(mat1->yield_stress*R0_*R0_/kparticle_, 2./3);
    id2_= surface->get_id();

    double G1 = E1/2/(1+v1);
    kT_ = 8/((2-v1)/G1)*0.001*R0_;

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
    F_ = update_normal_force(h);
    update_tangential_force(dt, normal);

}
unsigned DEM::Viscoelastic::M;
//Normal force
double DEM::Viscoelastic::update_normal_force(double h)
{
    double dh=h-h_;
    h_ = h - dh;
    if (h > hmax_) {
        hmax_ = h;
    }
    if(procent_ ) {
        if (h > -bt_ && h_ > -bt_) {
            dF_ = 3. / 2 * sqrt(h + bt_) * dh;
            area_ = pi * R0_ * (h + bt_);
            auto hn32 = pow(h_ + bt_, 3. / 2);
            auto h_32diff = pow(h + bt_, 3. / 2) - hn32;
            for (unsigned i = 0; i != M; ++i) {
                ddi_[i] = bi[i] * h_32diff + ai[i] * (hn32 - di_[i]);
                dF_ -= alpha_i[i] * ddi_[i];
                di_[i] += ddi_[i];
            }
            F_visc += k_ * dF_;
        }
        else {
            F_particle = 0;
            F_visc = 0;
            dF_ = 0;
            area_ = 0;
            for (unsigned i = 0; i != M; ++i) {
                ddi_[i] = 0;
                di_[i] = 0;
            }
            return 0;
        }
    }

    if (h > 0.0 && h_ > 0.0) {
        // Particles in contact
        // If unloading or displacement smaller than yield displacement, use Hertz
        if (h < hmax_ || h_ + dh < yield_h_) {
            F_particle += 1.5*kparticle_*sqrt(h_)*dh;
        }
            // Plastic contact, use a linear relationship with the stiffness obtained at the yield point
        else {
            F_particle += 1.5*kparticle_*sqrt(yield_h_)*dh;
        }
    }
    else {
        F_particle = 0;
    }

    h_ += dh;
    if (adhesive_) {
        return F_visc + std::max(F_particle, 0.);
    }
    return std::max(F_visc, 0.) + std::max(F_particle, 0.);
}



void DEM::Viscoelastic::update_tangential_force(const DEM::Vec3 &dt, const DEM::Vec3 &normal) {
    uT_ -= dot_product(uT_, normal) * normal;
    FT_part_ -= dot_product(FT_part_, normal) * normal;
    FT_visc_ -= dot_product(FT_visc_, normal) * normal;
    uT_ += dt;
   //std::cout<<uT_<<"uT"<<std::endl;
    if (F_visc != 0.0  && !dt.is_zero()) {
        dFT_ = 2*area_*tsi0_*dt/bt_;
        //std::cout<<dFT_<<"dFT"<<std::endl;

        for (unsigned i=0; i!=M; ++i)
        {
            ddti_[i] = bi[i]*dt + ai[i]*(uT_-dti_[i]);
            dFT_ -= alpha_i[i]*ddti_[i];
            dti_[i]+=ddti_[i];
        }
        FT_visc_ += dFT_;

        if (FT_visc_.length() > F_visc) { // contact aborted
            FT_visc_.set_zero();
            uT_.set_zero();
            for (unsigned i=0; i!=M; ++i)
            {
                dti_[i].set_zero();
            }
        }
    }
     else if (F_particle > 0.0) {
        FT_part_ += kT_*dt;
        //kT_: G*8*sqrt(R*angel of impact) G=60GPa
        if (FT_part_.length() > mu_*F_) { // Slip
            FT_part_ = mu_*F_*FT_part_.normal();
        }
    }
    else {
        FT_.set_zero();
        uT_.set_zero();
        FT_part_.set_zero();
        FT_visc_.set_zero();
        for (unsigned i=0; i!=M; ++i)
        {
            dti_[i].set_zero();
            ddti_[i].set_zero();
        }
    }
    FT_ = - FT_visc_ - FT_part_;

}

std::string DEM::Viscoelastic::get_output_string() const {
    std::stringstream ss;
    ss  << F_visc+F_particle;
    return ss.str();
}

void DEM::Viscoelastic::set_increment(std::chrono::duration<double> dt) {
    dt_ = dt.count();
    ai = {};
    bi = {};
    for (unsigned i=0; i!=M; ++i) {
        ai.push_back(1-exp((-dt_/tau_i[i])));
        bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-ai[i]));
    }
}
