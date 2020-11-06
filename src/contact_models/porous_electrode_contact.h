//
// Created by erolsson on 23/10/2020.
//

#ifndef DEMSIM_POROUS_ELECTRODE_CONTACT_H
#define DEMSIM_POROUS_ELECTRODE_CONTACT_H

#include <chrono>
#include <vector>

#include "../particles/spherical_particle.h"
#include "../surfaces/surface_base.h"
#include "../surfaces/point_surface.h"
#include "../utilities/vec3.h"

namespace DEM {
    class PorousElectrodeMaterial;

    class PorousElectrodeContact {
        using ParticleType = SphericalParticle<PorousElectrodeContact>;
        using SurfaceType = Surface<PorousElectrodeContact, ParticleType>;
    public:
        PorousElectrodeContact(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double> dt);
        PorousElectrodeContact(ParticleType* particle1, SurfaceType* surface, std::chrono::duration<double>dt );

        PorousElectrodeContact(ParticleType*, ParticleType*, std::chrono::duration<double>,
                const ParameterMap& parameters);
        PorousElectrodeContact(ParticleType*, SurfaceType*, std::chrono::duration<double>,
                const ParameterMap& parameters);

        void update(double h, const Vec3& dt, const Vec3& drot, const Vec3& normal);

        [[nodiscard]] double get_overlap() const { return h_; }
        [[nodiscard]] double get_normal_force() const { return F_; }
        [[nodiscard]] const Vec3& get_tangential_force() const { return FT_; }
        [[nodiscard]] double get_contact_area() const {return pi*a_*a_; }
        [[nodiscard]] Vec3 get_rolling_resistance_torque() const;
        [[nodiscard]] bool active() const {return F_ != 0; }
        [[nodiscard]] std::string get_output_string() const;
        [[nodiscard]] std::string restart_data() const;
        void set_increment(std::chrono::duration<double>);

    private:
        double h_ = -1e99;
        double F_ = 0;
        double Fvisc_ = 0;
        double Fparticle_ = 0;
        Vec3 FT_ = Vec3(0, 0, 0);
        double a_ = 0;

        bool binder_contact_;
        bool activated_ = false;

        double kparticle_;
        double kbinder_;
        double R0_;
        double bt_;         // Thickness of the binder disk
        double br_;         // Radius of the binder disk

        static unsigned M;
        double dt_;
        std::vector<double> tau_i {};
        std::vector<double> alpha_i {};
        std::vector<double> ai {};
        std::vector<double> bi {};

        std::vector<double> di_ {};
        std::vector<double > ddi_ {};

        static bool create_binder_contact(const PorousElectrodeMaterial* mat);
    };

}

#endif //DEMSIM_POROUS_ELECTRODE_CONTACT_H
