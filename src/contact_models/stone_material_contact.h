//
// Created by erolsson on 2019-01-20.
//

#ifndef DEMSIM_STONE_MATERIAL_CONTACT_H
#define DEMSIM_STONE_MATERIAL_CONTACT_H

#include <random>
#include "../particles/fractureable_spherical_particle.h"
#include "../utilities/vec3.h"

namespace DEM {
    class StoneMaterialContact {
        using ParticleType = FractureableSphericalParticle<StoneMaterialContact>;
        using SurfaceType = Surface<StoneMaterialContact, ParticleType>;

    public:
        StoneMaterialContact(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double>);
        StoneMaterialContact(ParticleType* particle1, SurfaceType* surface, std::chrono::duration<double>);
        void update(double dh, const Vec3& dt, const Vec3& rot, const Vec3& normal);
        [[nodiscard]] double get_overlap() const { return h_; }
        [[nodiscard]] double get_normal_force() const { return FN_; }
        [[nodiscard]] const Vec3& get_tangential_force() const { return FT_; }
        [[nodiscard]] const Vec3& get_rolling_resistance_torque() const { return Tr_; };
        [[nodiscard]] double get_contact_area() const { return sqrt(a_); }
        [[nodiscard]] bool active() const { return FN_!=0; }

        static void set_increment(std::chrono::duration<double>) {}
        [[nodiscard]] std::string get_output_string() const;

    private:
        double update_normal_force(double dh);
        void update_tangential_force(const Vec3& dt, const Vec3& normal, double dFN);
        void update_tangential_resistance(const Vec3& rot);

        void assign_fracture_strengths();
        void check_fracture(const Vec3& normal);

        ParticleType* p1_;
        ParticleType* p2_;

        double FN_ { 0. };
        double a_ { 0, };
        double R0_;

        double Fs_;
        double hs_;

        double hlinear_;

        // state parameters in the normal direction
        double h_ { 0. };
        double hmax_ { 0. };
        double Fmax_ {0. };
        double hp_ { 0. };
        double hl_ {0. };

        double h1_;

        // state parameters in the tangential direction
        Vec3 FT_ {Vec3{}};
        Vec3 old_dT_ {Vec3{}};                      // Tangential displacement from previous increment
        Vec3 load_var_tangential_ {Vec3{}};         // Accumulation variable for solving lhs in Eq. (?)
        double load_var_normal_ {0.};

        int turning_point_ { 1 };                 //First sign shift -1, second 1
        std::array<Vec3, 2> FT0 {Vec3{}, Vec3{}};   //Tangential loads at the turning points

        // Material parameters for the stone contact model in normal direction
        double kp_;                                 //Loading stiffness
        double ke_;
        double kl_;
        double ku_;                                 //Unloading stiffness
        double unloading_exp_;

        // Material parameters for the stone contact model in the tangential direction
        double mu_;
        double old_mu_;
        double kT_;

        // Parameters for the rolling resistance model
        double mu_r_;
        double kr_;
        Vec3 Tr_ {Vec3{}};

        double fracture_strength_p1_;
        double fracture_strength_p2_;
        std::size_t id2_;
    };
}

#endif //DEMSIM_STONE_MATERIAL_CONTACT_H
