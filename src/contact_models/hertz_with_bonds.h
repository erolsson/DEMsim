//
// Created by erolsson on 05/08/22.
//

#ifndef DEMSIM_HERTZ_WITH_BONDS_H
#define DEMSIM_HERTZ_WITH_BONDS_H
#include <chrono>

#include "../utilities/vec3.h"
#include "../particles/spherical_particle.h"
#include "../surfaces/surface_base.h"

namespace DEM {
    class ElasticBondedMaterial;
    class HertzWithBonds {
        using ParticleType = SphericalParticle<HertzWithBonds>;
        using SurfaceType = Surface<HertzWithBonds, ParticleType>;

    public:
        HertzWithBonds(ParticleType*, ParticleType*, std::chrono::duration<double>);
        HertzWithBonds(ParticleType*, SurfaceType*, std::chrono::duration<double>);

        void update(double d, const Vec3& dt, const Vec3&, const Vec3& normal);

        [[nodiscard]] double get_overlap() const { return h_; }
        [[nodiscard]] double get_normal_force() const { return F_; }
        [[nodiscard]] const Vec3& get_tangential_force() const { return FT_; }
        [[nodiscard]] Vec3 get_rolling_resistance_torque() const { return Vec3{};};
        [[nodiscard]] double get_contact_area() const {return sqrt(h_*R0_); }
        [[nodiscard]] bool active() const {return h_ > 0; }

        static void set_increment(std::chrono::duration<double>) {}
        [[nodiscard]] std::string get_output_string() const;
    private:
        double h_{ 0. };
        double kHertz_;
        double h_max = 0;
        double k_bond_ = 0;
        double kt_bond_ = 0;
        double c_bond_ = 1e99;
        double R0_;
        double ky_;
        double hy_;

        double kT_;
        double mu_;
        double F_{ 0 };
        double F_bond_ { 0 };
        Vec3 FT_bond_ {0, 0, 0};
        bool fractured_ {false };
        bool bonded_ {false};
        double increment_;
        double bond_area_;
        double fracture_stress_ = 0;

        const ElasticBondedMaterial* material1;
        const ElasticBondedMaterial* material2;

        Vec3 FT_{ Vec3(0., 0., 0.) };
        Vec3 uT_{ Vec3(0., 0., 0.) };

        bool bonded() const;
    };
}

#endif //DEMSIM_HERTZ_WITH_BONDS_H
