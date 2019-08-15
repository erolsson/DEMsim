//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_ENGINE_H
#define DEMSIM_ENGINE_H

#include <atomic>
#include <chrono>
#include <iostream>
#include <memory>
#include <vector>

#include "../utilities/amplitude.h"
#include "collision_detection/collision_detector.h"
#include "contact.h"
#include "../utilities/contact_matrix.h"
#include "../surfaces/cylinder.h"
#include "../materials/material_base.h"
#include "output.h"
#include "../surfaces/point_surface.h"
#include "../surfaces/surface_base.h"
#include "../utilities/vec3.h"

namespace DEM {

    template<typename ForceModel, typename ParticleType>
    class Engine {
    public:
        using ParticlePointer = ParticleType*;
        using PointSurfacePointer = PointSurface<ForceModel, ParticleType>*;
        using CylinderPointer = Cylinder<ForceModel, ParticleType>*;
        using OutputPointerType = std::shared_ptr<Output<ForceModel, ParticleType>>;

        explicit Engine(std::chrono::duration<double> dt);

        void setup();
        template<typename Condition>
        void run(Condition& condition);

        //Object creation functions
        template<typename MaterialType>
        MaterialType* create_material(double density);

        ParticlePointer create_particle(double radius, const Vec3& position, const Vec3& velocity,
                                      MaterialBase* material);

        PointSurfacePointer create_point_surface(const std::vector<Vec3>& points, bool infinite);

        CylinderPointer create_cylinder(double radius, const Vec3& axis, const Vec3& base_point, double length,
                                        bool inward=true, bool infinite=false, bool closed_ends=false);

        OutputPointerType create_output(std::string directory, std::chrono::duration<double> interval);
        void remove_output(const OutputPointerType& output_to_remove);

        std::shared_ptr<Amplitude> set_force_control_on_surface(Surface<ForceModel, ParticleType>* surface,
                char direction,  bool global_time=false);
        void remove_force_control_on_surface(Surface<ForceModel, ParticleType>* surface, char direction);

        std::pair<double, std::size_t> set_viscocity_parameters(double viscosity, size_t order=1);
        void remove_viscosity_parameters(std::pair<double, std::size_t> parameter_pair);

        // Getters
        std::chrono::duration<double> get_time() const { return time_; }
        double get_kinetic_energy() const;
        std::pair<size_t, double> max_particle_velocity() const;
        std::array<double, 6> get_bounding_box() const;

        // Setters
        void set_gravity(const Vec3& g) { gravity_ = g; }

        void set_mass_scale_factor(double factor) { mass_scale_factor_ = factor; }

        // Functors for running a simulation until a condition is fulfilled
        class RunForTime {
        public:
            RunForTime(const Engine& e, std::chrono::duration<double> time) :
                engine_{e}, start_time_{engine_.get_time()}, time_to_run_{time} {}

            void reset(std::chrono::duration<double> new_run_time) {
                start_time_ = engine_.get_time();
                time_to_run_ = new_run_time;
            }

            bool operator()() const
            {
                using namespace std::chrono_literals;
                return std::chrono::abs(time_to_run_  - (engine_.get_time() - start_time_ )) > 0.1ns;
            }

        private:
            const Engine& engine_;
            std::chrono::duration<double> start_time_ ;
            std::chrono::duration<double> time_to_run_;
        };

        class KineticEnergyLess {
        public:
            KineticEnergyLess(const Engine& e, double kinetic_energy, std::chrono::duration<double> update_time) :
                    engine_{e}, kinetic_energy_{kinetic_energy}, start_time_(engine_.get_time()),
                    update_time_{update_time} {}

            void set_new_value(double new_kinetic_energy) { kinetic_energy_ = new_kinetic_energy; }
            bool operator()() {
                if (engine_.get_time() - start_time_ < update_time_) {
                    return true;
                }
                start_time_ += update_time_;
                return engine_.get_kinetic_energy() > kinetic_energy_;
            }

        private:
            const Engine& engine_;
            double kinetic_energy_;

            std::chrono::duration<double> start_time_ ;
            std::chrono::duration<double> update_time_;
        };

        class ParticleVelocityLess {
        public:
            ParticleVelocityLess(const Engine& e, double max_velocity, std::chrono::duration<double> update_time) :
                    engine_{e}, max_velocity_{max_velocity}, start_time_(engine_.get_time()),
                    update_time_{update_time} {}

            void set_new_value(double new_max_velocity) { max_velocity_ = new_max_velocity; }
            bool operator()() {
                if (engine_.get_time() - start_time_ < update_time_) {
                    return true;
                }
                start_time_ += update_time_;
                return engine_.max_particle_velocity().second > max_velocity_;
            }

        private:
            const Engine& engine_;
            double max_velocity_;

            std::chrono::duration<double> start_time_ ;
            std::chrono::duration<double> update_time_;
        };

    private:
        using ContactType = Contact<ForceModel, ParticleType>;
        using SurfaceType = Surface<ForceModel, ParticleType>;

        std::size_t number_of_objects_{ 0 };
        std::chrono::duration<double> time_ { 0. };

        std::vector<MaterialBase*> materials_{};
        std::vector<ParticleType*> particles_{};
        std::vector<SurfaceType*> surfaces_{};
        ContactMatrix<ContactType> contacts_{};
        std::vector<OutputPointerType> outputs_{};

        CollisionDetector<ForceModel, ParticleType> collision_detector_;

        // Settings type of private data
        Vec3 gravity_ {Vec3{0,0,0}};
        std::vector<std::pair<double, std::size_t>> viscocity_parameters_;
        std::chrono::duration<double> increment_;
        double mass_scale_factor_ { 1. };


        void do_step();

        // Helper functions
        void move_particles();
        void move_surfaces();
        void create_contacts();
        void destroy_contacts();
        void update_contacts();
        void sum_contact_forces();
        void run_output();

        friend class Output<ForceModel, ParticleType>;
    };

    // Functors for different running conditions

}

#include "engine.tpp"

#endif //DEMSIM_ENGINE_H
