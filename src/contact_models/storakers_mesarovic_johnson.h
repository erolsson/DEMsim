//
// Created by erolsson on 2019-01-05.
//

#ifndef DEMSIM_STORAKERSMESAROVICJOHNSON_H
#define DEMSIM_STORAKERSMESAROVICJOHNSON_H

#include <chrono>

#include "../particles/spherical_particle.h"
#include "../surfaces/surface_base.h"
#include "../utilities/vec3.h"


namespace DEM {
    class StorakersMesarovicJohnson {
        using ParticleType = SphericalParticle<StorakersMesarovicJohnson>;
        using SurfaceType = Surface<StorakersMesarovicJohnson, ParticleType>;

    public:
        StorakersMesarovicJohnson(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double>);
        StorakersMesarovicJohnson(ParticleType* particle1, SurfaceType* surface, std::chrono::duration<double>);

        void update(double h, const Vec3& dt, const Vec3& normal);

        double get_overlap() const { return h_; }
        double get_normal_force() const { return F_; }
        const Vec3& get_tangential_force() const { return FT_; }
        double get_contact_area() const {return sqrt(a_); }
        bool active() const {return F_ != 0; }

        void set_increment(std::chrono::duration<double>) {}

    private:
        double h_{ 0. };
        double h_max_ { 0. };
        double a_ { 0. };
        double a_max_ { 0. };

        double k_;
        double ku_;
        double hu_max_ { 0. };

        double R0_;
        double kT_;
        double mu_;
        double F_{ 0 };

        Vec3 FT_{ Vec3(0., 0., 0.) };
        Vec3 uT_{ Vec3(0., 0., 0.) };

        void update_normal_force(double h);
        void update_tangential_force(const Vec3& dt, const Vec3& normal);
    };
}
#endif //DEMSIM_STORAKERSMESAROVICJOHNSON_H
