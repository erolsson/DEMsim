//
// Created by erolsson on 23/10/2020.
//

#include "porous_electrode_contact.h"

#include <random>

#include "../materials/porous_electrode_material.h"

DEM::PorousElectrodeContact::PorousElectrodeContact(DEM::PorousElectrodeContact::ParticleType* particle1,
                                                    DEM::PorousElectrodeContact::ParticleType* particle2,
                                                    std::chrono::duration<double> dt) :
    dt_(dt.count()){
    auto mat1 = dynamic_cast<const DEM::PorousElectrodeMaterial*>(particle1->get_material());
    auto mat2 = dynamic_cast<const DEM::PorousElectrodeMaterial*>(particle2->get_material());
    double Ep1 = mat1->E_particle;
    double vp1 = mat1->v_particle;

    double Ep2 = mat2->E_particle;
    double vp2 = mat2->v_particle;

    double Eb1 = mat1->E_binder;
    double vb1 = mat1->v_binder;
    bt_ = mat1->binder_thickness;
    R0_ = 1/(1/particle1->get_radius() + 1/particle2->get_radius());
    br_ = mat1->binder_radius_fraction*2*R0_;
    double A = DEM::pi*br_*br_;
    double E0 = 1/((1-vp1*vp1)/Ep1 + (1-vp2*vp2)/Ep2);
    kparticle_ = 2*E0*sqrt(R0_);
    kbinder_ = (1-vb1)/(1+vb1)/(1-2*vb1)*Eb1/bt_*A;

    binder_contact_ = create_binder_contact(mat1);
    M = mat1->alpha_i.size();
    alpha_i = mat1->alpha_i;

    for (unsigned i=0; i!=M; ++i)
    {
        di_.push_back(0);
        ddi_.push_back(0);
        ai.push_back(1-exp((-dt_/mat1->tau_i[i])));
        bi.push_back(mat1->tau_i[i]/dt_*((dt_/mat1->tau_i[i])-mat1->alpha_i[i]));
    }
}

DEM::PorousElectrodeContact::PorousElectrodeContact(DEM::PorousElectrodeContact::ParticleType* particle1,
                                                    DEM::PorousElectrodeContact::SurfaceType* surface,
                                                    std::chrono::duration<double> dt) {

}

DEM::PorousElectrodeContact::PorousElectrodeContact(DEM::PorousElectrodeContact::ParticleType*,
                                                    DEM::PorousElectrodeContact::ParticleType*,
                                                    std::chrono::duration<double>,
                                                    const DEM::ParameterMap& parameters) {

}

DEM::PorousElectrodeContact::PorousElectrodeContact(DEM::PorousElectrodeContact::ParticleType*,
                                                    DEM::PorousElectrodeContact::SurfaceType*,
                                                    std::chrono::duration<double>,
                                                    const DEM::ParameterMap& parameters) {

}

void
DEM::PorousElectrodeContact::update(double h, const DEM::Vec3& dt, const DEM::Vec3& drot, const DEM::Vec3& normal) {
    double dh = h - h_;
    h_ = h;
    if (binder_contact_) {
        if ((h > -bt_) || activated_) {
            activated_ = true;
            double dF = dh;
            for (unsigned i = 0; i != M; ++i) {
                ddi_[i] = bi[i]*dh + ai[i]*(h_ + bt_ -di_[i]);
                dF -= alpha_i[i]*ddi_[i];
                di_[i] += ddi_[i];
            }
            Fvisc_ += kbinder_*dF;
        }
    }
    if (h_ > 0) {
        Fparticle_ += kparticle_*sqrt(h_)*dh;
    }
    F_ = std::max(Fparticle_, 0.) + Fvisc_;
}

bool DEM::PorousElectrodeContact::create_binder_contact(const DEM::PorousElectrodeMaterial* mat) {
    std::random_device random_device;
    std::default_random_engine rand_engine { random_device() };
    std::uniform_real_distribution<double> distribution{0., 1.};
    double random_value = distribution(rand_engine);
    if (random_value < mat->fraction_binder_contacts){
        return true;
    }
    return false;
}

unsigned DEM::PorousElectrodeContact::M;
