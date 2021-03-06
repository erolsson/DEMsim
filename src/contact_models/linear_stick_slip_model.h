//
// Created by erolsson on 2018-07-30.
//

#ifndef DEMSIM_LINEAR_STICK_SLIP_MODEL_H
#define DEMSIM_LINEAR_STICK_SLIP_MODEL_H

#include <chrono>

#include "../utilities/vec3.h"
#include "../particles/spherical_particle.h"
#include "../surfaces/surface_base.h"

namespace DEM {
    class LinearStickSlipModel {
        using ParticleType = SphericalParticle<LinearStickSlipModel>;
        using SurfaceType = Surface<LinearStickSlipModel, ParticleType>;

    public:
        LinearStickSlipModel(ParticleType*, ParticleType*, std::chrono::duration<double>);
        LinearStickSlipModel(ParticleType*, SurfaceType*, std::chrono::duration<double>);

        void update(double dh, const Vec3& dt, const Vec3&, const Vec3& normal);

        [[nodiscard]] double get_overlap() const { return h_; }
        [[nodiscard]] double get_normal_force() const { return F_; }
        [[nodiscard]] const Vec3& get_tangential_force() const { return FT_; }
        [[nodiscard]] Vec3 get_rolling_resistance_torque() const { return Vec3{};};
        [[nodiscard]] double get_contact_area() const {return sqrt(h_*R0_); }
        [[nodiscard]] bool active() const {return h_>0; }

        static void set_increment(std::chrono::duration<double>) {}
        [[nodiscard]] std::string get_output_string() const;
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
