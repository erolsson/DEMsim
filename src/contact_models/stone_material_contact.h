//
// Created by erolsson on 2019-01-20.
//

#ifndef DEMSIM_STONE_MATERIAL_CONTACT_H
#define DEMSIM_STONE_MATERIAL_CONTACT_H

#include "../particles/fractureable_spherical_particle.h"
#include "../utilities/vec3.h"

namespace DEM {
    class StoneMaterialContact {
        using ParticleType = FractureableSphericalParticle<StoneMaterialContact>;
        using SurfaceType = Surface<StoneMaterialContact, ParticleType>;

    public:
        StoneMaterialContact(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double>);
        StoneMaterialContact(ParticleType* particle1, SurfaceType* surface, std::chrono::duration<double>);
        void update(double h, const Vec3& dt, const Vec3& normal);
        double get_overlap() const { return h_; }
        double get_normal_force() const { return F_; }
        const Vec3& get_tangential_force() const { return FT_; }
        double get_contact_area() const { return sqrt(a_); }
        bool active() const { return F_!=0; }
        void set_increment(std::chrono::duration<double>) { }

    private:
        void update_normal_force(double h);
        void update_tangential_force(const Vec3& dt, const Vec3& normal);

        double h_ { 0. };
        double hmax_ { 0. };
        double hp_ { 0. };
        double hl_ {0. };

        double F_ { 0. };
        double a_ { 0, };
        double R0_;

        Vec3 FT_ {Vec3(0., 0., 0.)};

        // Material parameters for the stone contact model
        double kp_;           //Loading stiffness
        double ke_;
        double kl_;
        double ku_;           //Unloading stiffness

        double unloading_exp_;
    };
}

#endif //DEMSIM_STONE_MATERIAL_CONTACT_H
