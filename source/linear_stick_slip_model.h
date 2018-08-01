//
// Created by erolsson on 2018-07-30.
//

#ifndef DEMSIM_LINEAR_STICK_SLIP_MODEL_H
#define DEMSIM_LINEAR_STICK_SLIP_MODEL_H

#include "vec3.h"
#include "spherical_particle.h"
#include "surface.h"
namespace DEM {
    class LinearStickSlipModel {
        using ParticleType = SphericalParticle<LinearStickSlipModel>;
        using SurfaceType = Surface<LinearStickSlipModel, ParticleType>;

    public:
        LinearStickSlipModel(ParticleType*, ParticleType*, double);
        LinearStickSlipModel(ParticleType*, SurfaceType*, double);

        void update(double h, const Vec3& dt);

        double get_overlap() const { return h_; }
        double get_normal_force() const { return F_; }
        const Vec3& get_tangential_force() const { return FT_; }
        double get_contact_area() const {return sqrt(h_*R0_); }
        double active() const {return h_>0; }

    private:
        double h_{ 0. };
        double k_;
        double R0_;

        double kT_;
        double mu_;
        double F_{ 0 };

        Vec3 FT_{ Vec3(0., 0., 0.) };
        Vec3 uT_{ Vec3(0., 0., 0.) };
    };
}

#endif //DEMSIM_LINEAR_STICK_SLIP_MODEL_H
