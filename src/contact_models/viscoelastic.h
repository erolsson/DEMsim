//
// Created by elahe on 2019-11-14.
//

#ifndef DEMSIM_VISCOELASTIC_H
#define DEMSIM_VISCOELASTIC_H

#include <chrono>
#include <vector>

#include "../particles/spherical_particle.h"
#include "../surfaces/surface_base.h"
#include "../surfaces/point_surface.h"
#include "../utilities/vec3.h"


namespace DEM {
    class Viscoelastic {
        using ParticleType = SphericalParticle<Viscoelastic>;
        using SurfaceType = Surface<Viscoelastic, ParticleType>;

    public:

        Viscoelastic(ParticleType* particle1, ParticleType* particle2,  std::chrono::duration<double> dt);
        Viscoelastic(ParticleType* particle1, SurfaceType* surface, std::chrono::duration<double>dt );

        void update(double h, const Vec3& dt, const Vec3& rot, const Vec3& normal);

        [[nodiscard]] double get_overlap() const { return h_; }
        [[nodiscard]] double get_normal_force() const {
            //std::cout << F_visc + F_particle << std::endl;
            return F_; }
        [[nodiscard]] const Vec3& get_tangential_force() const { return FT_; }
        [[nodiscard]] double get_contact_area() const {return area_; }
        [[nodiscard]] Vec3 get_rolling_resistance_torque() const { return Vec3{};};
        [[nodiscard]] bool active() const {return F_ != 0; }
        [[nodiscard]] std::string get_output_string() const;
        void set_increment(std::chrono::duration<double>);

    private:
        double dt_;   // Time increment
        double kT_;
        double bt_;
        double h_ {0. };
        double hmax_ { 0. };
        double area_ { 0. };
        static unsigned M;
        double yield_h_ { 1e99 };
        double k_;
        double kparticle_;
        //double binder_radii_;
        //double bindervolume_;
        double R0_;
        double F_{ 0 };
        double mu_;
        std::size_t id2_{};

        std::vector<double> tau_i {};
        std::vector<double> alpha_i {};
        std::vector<double> ai {};
        std::vector<double> bi {};

        std::vector<double> di_ {};
        std::vector<double > ddi_ {};
        std::vector<DEM::Vec3> ddti_ {};
        std::vector<DEM::Vec3> dti_ {};
        std::vector<double> x {};
        double dF_{0.};
        double F_visc{0.};
        double F_particle{0.};


        double tsi0_;
        double tsi0particle_;
        Vec3 dFT_{Vec3(0., 0., 0.)};
        Vec3 FT_{Vec3(0., 0., 0.)};
        Vec3 uT_{ Vec3(0., 0., 0.) };
        bool adhesive_;
        bool procent_;

        double update_normal_force(double h);
        void update_tangential_force(const Vec3& dt, const Vec3& normal);
    };
}




#endif //DEMSIM_VISCOELASTIC_H