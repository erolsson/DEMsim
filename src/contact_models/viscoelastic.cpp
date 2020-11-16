//
// Created by elahe on 2019-11-14.
//

#include "viscoelastic.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "../materials/electrode_material.h"



DEM::Viscoelastic::Viscoelastic(DEM::Viscoelastic::ParticleType *particle1,DEM::Viscoelastic::ParticleType* particle2,
                                std::chrono::duration<double> dt)
{
    //extracting from material
    auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
    auto mat2 = dynamic_cast<const ElectrodeMaterial *>(particle2->get_material());

    id_1 = particle1->get_id();
    id_2 = particle2->get_id();

    R0_ = 1. / (1. / particle1->get_radius() + 1. / particle2->get_radius());
    Rb_ = 1. / (1. / (particle1->get_radius() + mat1->bt/2) + 1. / (particle2->get_radius() + mat2->bt/2));


    double E1 = mat1->E;
    //double E2 = mat2->E;
    double v1 = mat1->nu;
    double v2 = mat2->nu;
    double vp1 = mat1->nup;
    double vp2 = mat2->nup;
    double Ep2 = mat2->Ep;
    double Ep1 = mat1->Ep;
    bt_ = mat1->bt; //Thickness of the binder disk

    kT_B_=E1*0.3*0.0016/(2*(1+v1)*bt_);
    //std::cout << "KT_B " <<kT_B_;

    stiff_b_=((1-v1)*E1)/(1+v1)/(1-2*v1);
    kB_=(0.3*0.0016*stiff_b_)/(bt_);
    //std::cout << "stiffness " << stiff_b_;

    //double tsi0 = 1. / (((1 - v1 * v1) / E1) + ((1 - v2 * v2) / E2));
    double tsi0particle = 1./(((1-vp1*vp1)/Ep1)+((1-vp2*vp2)/Ep2));

    binder_contact_ = create_binder_contact(mat1);
    mu_particle_ = (mat1->mu + mat2->mu)/2;
    mu_binder_ = std::min(mat1->mu_binder, mat2->mu_binder);
    adhesive_ = true;
    M = mat1->M();
    tau_i = mat1->tau_i;
    alpha_i = mat1->alpha_i;

    dt_ = dt.count();  // time increment

    //k_=4.*tsi0*sqrt(R0_ + bt_)/3; //initial contact stiffness
    //std::cout << "K_" << k_;

    kparticle_=4*tsi0particle*sqrt(R0_)/3;
    yield_h_ = 2*mat1->yield_displacement_coeff*R0_;

    double G1p = Ep1/2/(1+vp1);
    double G2p = Ep2/2/(1+vp2);
    kT_part_ = 8/((2 - v1)/G1p + (2 - v2)/G2p)*0.001*R0_;

    for (unsigned i=0; i!=M; ++i)
    {
        di_.push_back(0);
        ddti_.emplace_back(0., 0., 0.);
        dti_.emplace_back(0., 0., 0.);
        ddi_.push_back(0);
        ai.push_back(1-exp((-dt_/tau_i[i])));
        bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-ai[i]));
    }
}

DEM::Viscoelastic::Viscoelastic(DEM::Viscoelastic::ParticleType *particle1, DEM::Viscoelastic::SurfaceType * surface,
                                std::chrono::duration<double>dt){
    auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
    id_1 = particle1->get_id();
    id_2 = surface->get_id();
    R0_ = particle1->get_radius();
    Rb_ = particle1->get_radius() + mat1->bt/2;

    double E1 = mat1->E;
    double v1 = mat1->nu;
    double vp1=mat1->nup;
    double Ep1 = mat1->Ep;
    //double G1 = E1/2/(1+v1);
    bt_= mat1->bt;
    stiff_b_=((1-v1)*E1)/(1+v1)/(1-2*v1);
    kT_B_=E1*0.3*0.0016/(2*(1+v1)*bt_);
    kB_=(0.3*0.0016*stiff_b_)/(bt_);

    //std::cout << "KT_B " <<kT_B_;
    //std::cout << "stiffness " << stiff_b_;

    //double tsi0 = 1. / ((1 - v1 * v1) / E1);
    double tsi0particle = 1./((1-vp1*vp1)/Ep1);
    mu_particle_ = mat1->mu_wall;
    mu_binder_ = mat1->mu_binder;
    adhesive_ = surface->adhesive();

    M = mat1->M();
    tau_i=mat1->tau_i;
    alpha_i=mat1->alpha_i;


    binder_contact_ = create_binder_contact(mat1);
    //std::cout<< "binder_contact " <<binder_contact_;

    dt_ = dt.count();  // time increment
    //k_=4.*tsi0*sqrt(R0_ + bt_)/3; //initial contact stiffness
    kparticle_=4*tsi0particle*sqrt(R0_)/3;
    yield_h_ = 2*mat1->yield_displacement_coeff*R0_;

    double G1p = Ep1/2/(1+vp1);
    kT_part_ = 8/((2 - v1)/G1p)*0.001*R0_;

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

DEM::Viscoelastic::Viscoelastic(DEM::Viscoelastic::ParticleType* p1, DEM::Viscoelastic::ParticleType* p2,
                                std::chrono::duration<double>, const DEM::ParameterMap& parameters) :
        id_1(p1->get_id()),
        id_2(p2->get_id()),
        dt_(parameters.get_parameter<double>("dt")),  // Time increment
        kT_part_(parameters.get_parameter<double>("kT_part")),
        bt_(parameters.get_parameter<double>("bt")),
        h_(parameters.get_parameter<double>("h")),
        hmax_(parameters.get_parameter<double>("hmax")),
        yield_h_(parameters.get_parameter<double>("yield_h")),
        //k_(parameters.get_parameter<double>("k")),
        kB_(parameters.get_parameter<double>("kB")),
        kT_B_(parameters.get_parameter<double>("kT_B")),
        kparticle_(parameters.get_parameter<double>("kparticle")),
        R0_(parameters.get_parameter<double>("R0")),
        Rb_(parameters.get_parameter<double>("Rb")),
        F_(parameters.get_parameter<double>("F")),
        mu_particle_(parameters.get_parameter<double>("mu_particle")),
        mu_binder_(parameters.get_parameter<double>("mu_binder")),
        dF_(parameters.get_parameter<double>("dF")),
        F_visc(parameters.get_parameter<double>("F_visc")),
        F_particle(parameters.get_parameter<double>("F_particle")),
        dFT_(parameters.get_vec3("dFT")),
        FT_(parameters.get_vec3("FT")),
        FT_visc_(parameters.get_vec3("FT_visc")),
        FT_part_(parameters.get_vec3("FT_part")),
        uT_(parameters.get_vec3("uT")),
        rot_(parameters.get_vec3("rot")),
        adhesive_(parameters.get_parameter<bool>("adhesive")),
        binder_contact_(parameters.get_parameter<bool>("binder_contact")),
        fractured_(parameters.get_parameter<bool>("fractured"))

{
    M = parameters.get_parameter<unsigned>("M");
    for (unsigned i=0; i != M; ++i) {
        tau_i.push_back(parameters.get_parameter<double>("tau_" + std::to_string(i)));
        alpha_i.push_back(parameters.get_parameter<double>("alpha_" + std::to_string(i)));
        ai.push_back(parameters.get_parameter<double>("a_" + std::to_string(i)));
        bi.push_back(parameters.get_parameter<double>("b_" + std::to_string(i)));
        di_.push_back(parameters.get_parameter<double>("d_" + std::to_string(i)));
        ddi_.push_back(parameters.get_parameter<double>("dd_" + std::to_string(i)));
        dti_.push_back(parameters.get_vec3("dt_" + std::to_string(i)));
        ddti_.push_back(parameters.get_vec3("ddt_" + std::to_string(i)));
    }
}

DEM::Viscoelastic::Viscoelastic(DEM::Viscoelastic::ParticleType* p, DEM::Viscoelastic::SurfaceType* s,
                                std::chrono::duration<double>, const DEM::ParameterMap& parameters) :
        id_1(p->get_id()),
        id_2(s->get_id()),
        dt_(parameters.get_parameter<double>("dt")),  // Time increment
        kT_part_(parameters.get_parameter<double>("kT_part")),
        bt_(parameters.get_parameter<double>("bt")),
        h_(parameters.get_parameter<double>("h")),
        hmax_(parameters.get_parameter<double>("hmax")),
        yield_h_(parameters.get_parameter<double>("yield_h")),
        //k_(parameters.get_parameter<double>("k")),
        kB_(parameters.get_parameter<double>("kB")),
        kT_B_(parameters.get_parameter<double>("kT_B_")),
        kparticle_(parameters.get_parameter<double>("kparticle")),
        R0_(parameters.get_parameter<double>("R0")),
        Rb_(parameters.get_parameter<double>("Rb")),
        F_(parameters.get_parameter<double>("F")),
        mu_particle_(parameters.get_parameter<double>("mu_particle")),
        mu_binder_(parameters.get_parameter<double>("mu_binder")),
        dF_(parameters.get_parameter<double>("dF")),
        F_visc(parameters.get_parameter<double>("F_visc")),
        F_particle(parameters.get_parameter<double>("F_particle")),
        dFT_(parameters.get_vec3("dFT")),
        FT_(parameters.get_vec3("FT")),
        FT_visc_(parameters.get_vec3("FT_visc")),
        FT_part_(parameters.get_vec3("FT_part")),
        uT_(parameters.get_vec3("uT")),
        rot_(parameters.get_vec3("rot")),
        adhesive_(parameters.get_parameter<bool>("adhesive")),
        binder_contact_(parameters.get_parameter<bool>("binder_contact")),
        fractured_(parameters.get_parameter<bool>("fractured"))
{
    M = parameters.get_parameter<std::size_t>("M");
    for (unsigned i=0; i != M; ++i) {
        tau_i.push_back(parameters.get_parameter<double>("tau_" + std::to_string(i)));
        alpha_i.push_back(parameters.get_parameter<double>("alpha_" + std::to_string(i)));
        ai.push_back(parameters.get_parameter<double>("a_" + std::to_string(i)));
        bi.push_back(parameters.get_parameter<double>("b_" + std::to_string(i)));
        di_.push_back(parameters.get_parameter<double>("d_" + std::to_string(i)));
        ddi_.push_back(parameters.get_parameter<double>("dd_" + std::to_string(i)));
        dti_.push_back(parameters.get_vec3("dt_" + std::to_string(i)));
        ddti_.push_back(parameters.get_vec3("ddt_" + std::to_string(i)));
    }
}

void DEM::Viscoelastic::update(double h, const DEM::Vec3& dt, const Vec3& drot, const DEM::Vec3& normal){
    rot_ += drot;
    F_ = update_normal_force(h);
    update_tangential_force(dt, normal);
}

unsigned DEM::Viscoelastic::M;


//Normal force
double DEM::Viscoelastic::update_normal_force(double h)
{
    double dh = h - h_;
    h_ = h - dh;
    if (h > hmax_) {
        hmax_ = h;
        fractured_ = false;
    }
    if (binder_contact_) {

        if (h > -bt_ && h_ > -bt_) {
            dF_ = dh;

            for (unsigned i = 0; i != M; ++i) {
                ddi_[i] = bi[i] * dh + ai[i] * (h_+bt_ - di_[i]);
                dF_ -= alpha_i[i] * ddi_[i];
                di_[i] += ddi_[i];
            }

            F_visc += kB_*dF_;

        }
        else {
            F_visc = 0;
            dF_ = 0;
            rot_.set_zero();
            for (unsigned i = 0; i != M; ++i) {
                ddi_[i] = 0;
                di_[i] = 0;
            }
        }
    }

    if (h > 0.0 && h_ > 0.0) {
        // Particles in contact
        // If h > yield_h and loading i. e that dh > 0 plastic loading, use stiffness at yield point

        // Otherwise use elastic formulation based on hertz

        if (h > yield_h_ && h >= hmax_) {
            F_particle += 1.5*kparticle_*sqrt(yield_h_)*dh;
        }
        else {
            F_particle += 1.5*kparticle_*sqrt(h_)*dh;
        }
    }
    else {
        F_particle = 0;
    }

    h_ += dh;
    // std::cout << "Adhesive: " << adhesive_ << ", Fvisc: " << F_visc << ", fractured: " << fractured_ <<  "\n";
    //std::cout << ": " <<
    if (adhesive_ && !fractured_) {
        return F_visc + std::max(F_particle, 0.);
    }
    return std::max(F_visc, 0.) + std::max(F_particle, 0.);
}

void DEM::Viscoelastic::update_tangential_force(const DEM::Vec3 &dt, const DEM::Vec3 &normal) {

    uT_ -= dot_product(uT_, normal)*normal;
    FT_part_ -= dot_product(FT_part_, normal)*normal;
    if (!FT_visc_.is_zero()) {
        FT_visc_ -= dot_product(FT_visc_, normal)*normal;
        for (unsigned i = 0; i != M; ++i) {
            dti_[i] -= dot_product(dti_[i], normal)*normal;
        }
    }
    uT_ += dt;

    if (F_visc != 0.0) {
        dFT_ = dt;
        for (unsigned i = 0; i != M; ++i) {
            ddti_[i] = bi[i]*dt + ai[i]*(uT_-dti_[i]);
            dFT_ -= alpha_i[i]*ddti_[i];
            dti_[i] += ddti_[i];
        }
        FT_visc_ += kT_B_*dFT_;
        /*
        if (FT_visc_.length() > mu_binder_*abs(F_visc)) { // contact aborted
            fractured_ = true;
        }
         */
    }
    else {
        rot_.set_zero();
        FT_visc_.set_zero();
        for (unsigned i = 0; i != M; ++i) {
            dti_[i].set_zero();
        }
    }

    if (F_particle > 0.0) {
        FT_part_ += kT_part_*dt;
        if (FT_part_.length() > mu_particle_*F_particle && !uT_.is_zero() && !dt.is_zero()) { // Slip
            FT_part_ = mu_particle_*F_particle*(0.5*uT_.normal() + 0.5*dt.normal());
        }
    }
    else {
        uT_.set_zero();
        FT_part_.set_zero();
    }
    FT_.set_zero();
    if (fractured_ && !uT_.is_zero() && !dt.is_zero()){
        FT_visc_ = (0.5*uT_.normal() + 0.5*dt.normal())*mu_binder_*abs(F_visc);
        // std::cout << fractured_ << "  " << FT_visc_ << "  " << F_visc << "\n";
        for (unsigned i = 0; i != M; ++i) {
            dti_[i].set_zero();
        }
    }
    FT_ -= FT_visc_;
    FT_ -= FT_part_;
}

std::string DEM::Viscoelastic::get_output_string() const {
    std::stringstream ss;
    ss  << F_ << ", " << F_visc << ", " << F_particle << ", " << hmax_ << ", "
        << FT_.x() << ", "  << FT_.y() << ", " << FT_.z() << ", "
        << FT_visc_.x() << ", "  << FT_visc_.y() << ", " << FT_visc_.z() << ", "
        << FT_part_.x() << ", "  << FT_part_.y() << ", " << FT_part_.z() << ", "
        << binder_contact_ << ", " << bt_;
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

std::string DEM::Viscoelastic::restart_data() const {
    std::ostringstream ss;
    ss << named_print(dt_, "dt") << ", "
       << named_print(bt_, "bt") << ", "
       << named_print(kT_part_, "kT_part") << ", "
       << named_print(h_, "h") << ", "
       << named_print(hmax_, "hmax") << ", "
       << named_print(yield_h_, "yield_h") << ", "
       //<< named_print(k_, "k") << ", "
       << named_print(kparticle_, "kparticle") << ", "
       << named_print(R0_, "R0") << ", "
       << named_print(Rb_, "Rb") << ", "
       << named_print(F_, "F") << ", "
       << named_print(mu_particle_, "mu_particle") << ", "
       << named_print(mu_binder_, "mu_binder") << ", "
       << named_print(dF_, "dF") << ", "
       << named_print(F_visc, "F_visc") << ", "
       << named_print(F_particle, "F_particle") << ", "
       << named_print(dFT_, "dFT") << ", "
       << named_print(FT_, "FT") << ", "
       << named_print(FT_visc_, "FT_visc") << ", "
       << named_print(FT_part_, "FT_part") << ", "
       << named_print(uT_, "uT") << ", "
       << named_print(rot_, "rot") << ", "
       << named_print(adhesive_, "adhesive") << ", "
       << named_print(binder_contact_, "binder_contact") << ", "
       << named_print(fractured_, "fractured") << ", "
       << named_print(kB_, "kB") << ", "
       << named_print(kT_B_, "kT_B_") << ", "
       << named_print(M, "M");

    for (unsigned i=0; i != M; ++i) {
        ss <<  ", "
           << named_print(tau_i[i], "tau_" + std::to_string(i)) << ", "
           << named_print(alpha_i[i], "alpha_" + std::to_string(i)) << ", "
           << named_print(ai[i], "a_" + std::to_string(i)) << ", "
           << named_print(bi[i], "b_" + std::to_string(i)) << ", "
           << named_print(di_[i], "d_" + std::to_string(i)) << ", "
           << named_print(ddi_[i], "dd_" + std::to_string(i)) << ", "
           << named_print(dti_[i], "dt_" + std::to_string(i)) << ", "
           << named_print(ddti_[i], "ddt_" + std::to_string(i));
    }
    return ss.str();
}

bool DEM::Viscoelastic::create_binder_contact(const ElectrodeMaterial* mat) {
    std::random_device random_device;
    std::default_random_engine rand_engine { random_device() };
    std::uniform_real_distribution<double> distribution{0., 1.};
    double random_value = distribution(rand_engine);
    if (random_value < mat->fb){
        return true;
    }
    return false;
}

DEM::Vec3 DEM::Viscoelastic::get_rolling_resistance_torque() const {
    // if (binder_contact_) {
    //    return -Rb_*Rb_*0.01*kB_*rot_;
    // }
    // else {
    return DEM::Vec3(0, 0, 0);
    // }
}
