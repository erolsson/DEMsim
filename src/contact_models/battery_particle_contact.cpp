//
// Created by erolsson on 15/04/2020.
//

#include "battery_particle_contact.h"

#include <cmath>
#include <iostream>

#include "../materials/ViscoelasticMaterial.h"


DEM::BatteryParticleContact::BatteryParticleContact (DEM::BatteryParticleContact::ParticleType* particle1,
                                                     DEM::BatteryParticleContact::ParticleType* particle2,
                                                     std::chrono::duration<double> dt)  {
    //extracting from material
    auto mat1 = dynamic_cast<const ViscoelasticMaterial *>(particle1->get_material());
    auto mat2 = dynamic_cast<const ViscoelasticMaterial *>(particle2->get_material());


    R0_ = 1. / (1. / particle1->get_radius() + 1. / particle2->get_radius());
    double random = rand() % 10 + 1;


    //std::cout << "random:" << random << std::endl;
    procent_= true;

    double E1 = mat1->E;
    double E2 = mat2->E;
    double v1 = mat1->nu;
    double v2 = mat2->nu;
    double vp1 = mat1->nup;
    double vp2 = mat2->nup;
    double Ep2 = mat2->Ep;
    double Ep1 = mat1->Ep;
    //binder_radii_= (particle1->get_radius()*particle2->get_radius())/(particle1->get_radius()+particle2->get_radius());
    mu_ = (mat1->mu + mat2->mu)/2;
    adhesive_=true;
    M=mat1->M();
    tau_i=mat1->tau_i;
    alpha_i=mat1->alpha_i;
    kT_=mat1->kT;
    bt_= mat1->bt;

    //bindervolume_= mat1->bindervolume;
    dt_ = dt.count();  // time increment
    tsi0_ = 1. / (((1 - v1 * v1) / E1) + ((1 - v2 * v2) / E2));
    tsi0particle_=1./(((1-vp1*vp1)/Ep1)+((1-vp2*vp2)/Ep2));
    k_=4.*tsi0_*sqrt(R0_)/3; //initial contact stiffness
    kparticle_=4*tsi0particle_*sqrt(R0_)/3;


    yield_h_ = pow(std::min(mat1->yield_stress, mat2->yield_stress)*R0_*R0_/kparticle_, 2./3);
    std::cout << "procent:" << procent_ << std::endl;

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
DEM::BatteryParticleContact::BatteryParticleContact(DEM::BatteryParticleContact::ParticleType* particle1,
                                                    DEM::BatteryParticleContact::SurfaceType* surface,
                                                    std::chrono::duration<double>dt){
    auto mat1 = dynamic_cast<const ViscoelasticMaterial* >(particle1->get_material());

    R0_ = 1. / (1. / particle1->get_radius());
    double E1 = mat1->E;
    double v1 = mat1->nu;
    double vp1=mat1->nup;
    double Ep1 = mat1->Ep;
    mu_ = mat1->mu_wall;

    //bindervolume_= mat1->bindervolume;
    //binder_radii_= (particle1->get_radius()*(1.0/2));

    adhesive_ = surface->adhesive();

    M=mat1->M();
    tau_i=mat1->tau_i;
    alpha_i=mat1->alpha_i;
    kT_=mat1->kT;
    bt_= mat1->bt;
    double random = rand() % 10 + 1;
    //std::cout << "random:" << random << std::endl;
    procent_=true;


    dt_ = dt.count();  // time increment
    tsi0_ = 1. / ((1 - v1 * v1) / E1);
    //std::cout << "tsi0:" << tsi0_ << std::endl;
    tsi0particle_=1./((1-vp1*vp1)/Ep1);
    k_=4.*tsi0_*sqrt(R0_)/3; //initial contact stiffness
    //std::cout << "k:" << k_ << std::endl;
    kparticle_=4*tsi0particle_*sqrt(R0_)/3;
    id2_= surface->get_id();

    yield_h_ = pow(mat1->yield_stress*R0_*R0_/kparticle_, 2./3);

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
void DEM::BatteryParticleContact::update(double h, const DEM::Vec3& dt,const Vec3& , const DEM::Vec3& normal){
    F_ = update_normal_force(h);
    update_tangential_force(dt, normal);
}

unsigned DEM::BatteryParticleContact::M;

//Normal force
double DEM::BatteryParticleContact::update_normal_force(double h)
{
    double dh=h-h_;
    h_ = h - dh;
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

            if (!adhesive_ && F_visc < 0) {
                return 0;
            }
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
        }
        h_ += dh;
    }
    else if (h > 0){
        area_ = pi * R0_ * (h);
    }

    if (h > 0.0 && h_ > 0) {
        // Particles in contact
        // If unloading or displacement smaller than yield displacement, use Hertz
        if (dh < 0 || h_ + dh < yield_h_) {
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
    return F_visc + F_particle;
}



void DEM::BatteryParticleContact::update_tangential_force(const DEM::Vec3 &dt, const DEM::Vec3 &normal) {
    // Projecting uT on the new contact plane by removing the component in the contact normal direction
    uT_ -= dot_product(uT_, normal) * normal;
    uT_ += dt;
    if (F_visc != 0.0  && !uT_.is_zero()) {
        dFT_=-2*area_*tsi0_*dt/bt_;
        for (unsigned i=0; i!=M; ++i)
        {
            ddti_[i] = bi[i]*dt + ai[i]*(uT_-dti_[i]);
            dFT_ -= alpha_i[i]*ddti_[i];
            dti_[i]+=ddti_[i];
        }
        FT_+=dFT_;

        if (FT_.length() > 0.05*F_visc) { // contact aborted
            FT_.set_zero();
            uT_.set_zero();
            for (unsigned i=0; i!=M; ++i)
            {
                dti_[i].set_zero();
            }
        }
    }
    else if (F_particle > 0.0 ) {

        uT_ -= dot_product(uT_, normal)*normal;
        uT_ += dt;
        if (kT_*uT_.length() > mu_*F_) { // Slip
            uT_ = mu_*F_/kT_*uT_.normal();
        }
        FT_ = -kT_*uT_;
    }
    else {
        FT_.set_zero();
        uT_.set_zero();
        for (unsigned i=0; i!=M; ++i)
        {
            dti_[i].set_zero();
            ddti_[i].set_zero();
        }
    }

}

std::string DEM::BatteryParticleContact::get_output_string() const {
    std::stringstream ss;
    ss  << F_visc+F_particle;
    return ss.str();
}

void DEM::BatteryParticleContact::set_increment(std::chrono::duration<double> dt) {
    dt_ = dt.count();
    ai = {};
    bi = {};
    for (unsigned i=0; i!=M; ++i) {
        ai.push_back(1-exp((-dt_/tau_i[i])));
        bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-ai[i]));
    }
}
