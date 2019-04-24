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
        void update(double dh, const Vec3& dt, const Vec3& normal);
        double get_overlap() const { return h_; }
        double get_normal_force() const { return FN_; }
        const Vec3& get_tangential_force() const { return FT_; }
        double get_contact_area() const { return sqrt(a_); }
        bool active() const { return FN_!=0; }
        void set_increment(std::chrono::duration<double>) { }

    private:
        double update_normal_force(double dh);
        void update_tangential_force(const Vec3& dt, const Vec3& normal, double dFN);

        double FN_ { 0. };
        double a_ { 0, };
        double R0_;

        // state parameters in the normal direction
        double h_ { 0. };
        double hmax_ { 0. };
        double hp_ { 0. };
        double hl_ {0. };

        // state parameters in the tangential direction
        Vec3 FT_ {Vec3{}};
        Vec3 uT_ {Vec3{}};
        Vec3 old_dT_ {Vec3{}};    // Tangential displacement from previous increment

        char turning_point_ { -1 };  //First sign shift -1, second 1
        std::array<Vec3, 2> FT0 {Vec3{}, Vec3{}};

        // Material parameters for the stone contact model in normal direction
        double kp_;           //Loading stiffness
        double ke_;
        double kl_;
        double ku_;           //Unloading stiffness
        double unloading_exp_;

        // Material parameters for the stone contact model in the tangential direction
        double mu_;
        double old_mu_;
        double kT_;
    };
}

#endif //DEMSIM_STONE_MATERIAL_CONTACT_H
