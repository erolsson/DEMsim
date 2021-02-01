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
    class ElectrodeMaterial;
    class ParameterMap;
    class Viscoelastic {
        using ParticleType = SphericalParticle<Viscoelastic>;
        using SurfaceType = Surface<Viscoelastic, ParticleType>;

    public:
        Viscoelastic(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double> dt);
        Viscoelastic(ParticleType* particle1, SurfaceType* surface, std::chrono::duration<double>dt );

        Viscoelastic(ParticleType*, ParticleType*, std::chrono::duration<double>, const ParameterMap& parameters);
        Viscoelastic(ParticleType*, SurfaceType*, std::chrono::duration<double>, const ParameterMap& parameters);

        void update(double h, const Vec3& dt, const Vec3& drot, const Vec3& normal);

        [[nodiscard]] double get_overlap() const { return h_; }
        [[nodiscard]] double get_normal_force() const { return F_; }
        [[nodiscard]] const Vec3& get_tangential_force() const { return FT_; }
        [[nodiscard]] Vec3 get_rolling_resistance_torque() const;
        [[nodiscard]] bool active() const {return F_ != 0; }
        [[nodiscard]] std::string get_output_string() const;
        [[nodiscard]] std::string restart_data() const;
        void set_increment(std::chrono::duration<double>);

    private:
        const static ElectrodeMaterial* material;
        double kT_part_;
        double kB_;
        double kT_B_;
        double kT_;
        double kparticle_;
        double R0_;
        // double Rb_;
        double bt_;
        double h_ = 0.;
        double yield_h_ = 1e99 ;
        double hmax_ = -1e99 ;
        double mu_particle_;

        bool bonded_ = false;
        bool adhesive_ = true;
        bool binder_contact_ ;
        bool fractured_ = true;

        static unsigned M;
        double dt_;   // Time increment
        std::vector<double> tau_i {};
        std::vector<double> alpha_i {};
        std::vector<double> ai {};
        std::vector<double> bi {};

        std::vector<double> di_ {};
        std::vector<double > ddi_ {};
        std::vector<DEM::Vec3> ddti_ {};
        std::vector<DEM::Vec3> dti_ {};

        double F_ = 0.;
        double dF_ = 0.;
        double F_visc = 0.;
        double F_particle = 0.;
        Vec3 dFT_{Vec3(0., 0., 0.)};
        Vec3 FT_{Vec3(0., 0., 0.)};
        Vec3 FT_visc_ {Vec3(0., 0., 0.)};
        Vec3 FT_part_ {Vec3(0., 0., 0.)};
        Vec3 uT_{ Vec3(0., 0., 0.) };

        Vec3 rot_ {Vec3(0., 0., 0.)};


        double update_normal_force(double h);
        void update_tangential_force(const Vec3& dt, const Vec3& normal);
        static bool create_binder_contact(const ElectrodeMaterial* mat);
        bool adhesive() const;
    };
}




#endif //DEMSIM_VISCOELASTIC_H