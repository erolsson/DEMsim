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
    knp_ = 2*E0*sqrt(R0_);
    knb_ = (1 - vb1)/(1 + vb1)/(1 - 2*vb1)*Eb1/bt_*A;
    ktb_ = Eb1/bt_*A/2/(1+vb1);
    double G = Ep1/2/(1+vp1);
    ktp_ = 8*G*0.001*R0_/(2-vp1);
    mup_ = mat1->mu_particle;

    binder_contact_ = create_binder_contact(mat1);
    M = mat1->alpha_i.size();
    alpha_i = mat1->alpha_i;

    for (unsigned i=0; i!=M; ++i)
    {
        di_.push_back(0);
        ddi_.push_back(0);
        dti_.emplace_back(0, 0, 0);
        ddti_.emplace_back(0, 0, 0);
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

DEM::Vec3 DEM::PorousElectrodeContact::get_rolling_resistance_torque() const {
    return DEM::Vec3();
}

std::string DEM::PorousElectrodeContact::get_output_string() const {
    std::stringstream ss;
    ss << F_ << ", " << FNb_ << ", " << FNp_ << ", "
        << FT_.x() << ", "  << FT_.y() << ", " << FT_.z() << ", "
        << binder_contact_ << ", " << bt_;
    return ss.str();
}

std::string DEM::PorousElectrodeContact::restart_data() const {
    std::ostringstream ss;
    ss << named_print(h_, "h") << ", "
       << named_print(F_, "F") << ", "
       << named_print(FNb_, "Fvisc") << ", "
       << named_print(FNp_, "Fparticle") << ", "
       << named_print(FT_, "FT") << ", "
       << named_print(a_, "a") << ", "
       << named_print(binder_contact_, "binder_contact") << ", "
       << named_print(activated_, "activated") << ", "
       << named_print(knp_, "kparticle") << ", "
       << named_print(knb_, "kbinder") << ", "
       << named_print(R0_, "R0") << ", "
       << named_print(bt_, "bt") << ", "
       << named_print(br_, "br") << ", "
       << named_print(dt_, "dt") << ", "
       << named_print(M, "M");
    for (unsigned i=0; i != M; ++i) {
        ss <<  ", "
           << named_print(alpha_i[i], "alpha_" + std::to_string(i)) << ", "
           << named_print(ai[i], "a_" + std::to_string(i)) << ", "
           << named_print(bi[i], "b_" + std::to_string(i)) << ", "
           << named_print(di_[i], "d_" + std::to_string(i)) << ", "
           << named_print(ddi_[i], "dd_" + std::to_string(i));
    }
    return ss.str();
}

void DEM::PorousElectrodeContact::update_normal_force(double h) {
    double dh = h - h_;
    h_ = h;
    if (binder_contact_) {
        if ((h > -bt_) || activated_) {
            activated_ = true;
            double dF = dh;
            for (unsigned i = 0; i != M; ++i) {
                ddi_[i] = bi[i]*dh + ai[i]*(h_ + bt_ - di_[i]);
                dF -= alpha_i[i]*ddi_[i];
                di_[i] += ddi_[i];
            }
            FNb_ += knb_*dF;
        }
    }
    if (h_ > 0) {
        FNp_ += knp_*sqrt(h_)*dh;
    }
    else {
        FNp_ = 0;
    }
    F_ = std::max(FNp_, 0.) + FNb_;
}

void DEM::PorousElectrodeContact::update_tangential_force(const DEM::Vec3& dt, const Vec3& normal) {
    if (activated_) {
        FTb_ -= dot_product(FTb_, normal)*normal;
        uT_ -= dot_product(uT_, normal)*normal;
        uT_ += dt;

        Vec3 dFT_ = dt;
        for (unsigned i = 0; i != M; ++i) {
            ddti_[i] = bi[i]*dt + ai[i]*(uT_ - dti_[i]);
            dFT_ -= alpha_i[i]*ddti_[i];
            dti_[i] += ddti_[i];
        }
        FTb_ += ktb_*dFT_;
    }
    FT_ = -FTb_;

    if (FNp_ > 0) {
        FTp_ += ktp_*dt;
        if (FTp_.length() > mup_*FNp_) {
            FTp_ = mup_*FNp_*FTp_.normal();
        }
    }
    else {
        FTp_.set_zero();
    }
    FT_ -= FTp_;
}

unsigned DEM::PorousElectrodeContact::M;
