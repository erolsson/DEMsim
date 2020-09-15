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
#include "../surfaces/point_surface.h"
#include "collision_detection/collision_detector.h"
#include "contact.h"
#include "../utilities/contact_matrix.h"
#include "../surfaces/cylinder.h"

#include "../materials/material_base.h"
#include "output.h"
#include "../surfaces/surface_base.h"
#include "../utilities/vec3.h"

namespace DEM {
    class Amplitude;
    class ParameterMap;

    template<typename ForceModel, typename ParticleType>
    class PeriodicBCHandler;
    template<typename ForceModel, typename ParticleType>
    class Engine {
    public:
        using ParticlePointer = ParticleType*;
        using PointSurfacePointer = PointSurface<ForceModel, ParticleType>*;
        using DeformablePointSurfacePointer = DeformablePointSurface<ForceModel, ParticleType>*;
        using CylinderPointer = Cylinder<ForceModel, ParticleType>*;

        using OutputPointerType = std::shared_ptr<Output<ForceModel, ParticleType>>;
        using SurfaceType = Surface<ForceModel, ParticleType>;
        using PeriodicBCHandlerType = PeriodicBCHandler<ForceModel, ParticleType>;
        explicit Engine(std::chrono::duration<double> dt);
        explicit Engine(const std::string& restart_file_name);

        void setup();
        void setup(double bounding_box_stretch);
        template<typename Condition>
        void run(Condition& condition);

        //Object creation functions
        template<typename MaterialType>
        MaterialType* create_material(double density);

        ParticlePointer create_particle(double radius, const Vec3& position, const Vec3& velocity,
                                      MaterialBase* material);

        PointSurfacePointer create_point_surface(const std::vector<Vec3>& points, bool infinite, bool adhesive=true);
        PointSurfacePointer create_point_surface(const std::vector<Vec3>& points, bool infinite,
                                                 const std::string& name, bool adhesive=true);

        PointSurfacePointer create_point_surface(const std::vector<Vec3>& points, bool infinite,
                                                 const char* name, bool adhesive=true);

        CylinderPointer create_cylinder(double radius, const Vec3& axis, const Vec3& base_point, double length,
                                        bool inward=true, bool infinite=false, bool closed_ends=false);

        CylinderPointer create_cylinder(double radius, const Vec3& axis, const Vec3& base_point, double length,
                                        const std::string& name, bool inward=true, bool infinite=false,
                                        bool closed_ends=false);

        OutputPointerType create_output(std::string directory, std::chrono::duration<double> interval);
        OutputPointerType create_output(std::string directory, std::chrono::duration<double> interval,
                                        const std::string& name);
        void remove_output(const OutputPointerType& output_to_remove);

        std::shared_ptr<Amplitude> set_force_control_on_surface(Surface<ForceModel, ParticleType>* surface,
                char direction,  bool global_time=false);
        void remove_force_control_on_surface(Surface<ForceModel, ParticleType>* surface, char direction);

        [[maybe_unused]] std::pair<double, std::size_t> set_viscocity_parameters(double viscosity, size_t order=1);
        [[maybe_unused]] void remove_viscosity_parameters(std::pair<double, std::size_t> parameter_pair);

        [[maybe_unused]] void add_periodic_boundary_condition(char axis, double boundary_min, double boundary_max);
        [[maybe_unused]] void set_periodic_boundary_condition_strain_rate(char axis, double strain_rate);

        // Getters
        [[nodiscard]] std::chrono::duration<double> get_time() const { return time_; }
        [[nodiscard]] std::chrono::duration<double> get_time_increment() const { return increment_; }
        [[nodiscard]] std::size_t get_collision_object_count() const { return collision_id_counter_; }
        void increment_collision_counter() { ++collision_id_counter_; }

        [[nodiscard]] double get_kinetic_energy() const;
        [[nodiscard]] std::pair<size_t, double> max_particle_velocity() const;
        [[nodiscard]] std::pair<size_t, double> max_surface_velocity() const;
        [[nodiscard]] std::array<double, 6> get_bounding_box() const;

        [[maybe_unused]] [[nodiscard]] MaterialBase* get_material(std::size_t id) const;
        [[maybe_unused]] [[nodiscard]] SurfaceType* get_surface(std::size_t id) const;
        template<typename SurfaceT>
        [[maybe_unused]] [[nodiscard]] SurfaceT get_surface(const std::string& surface_name) const;
        [[maybe_unused]] [[nodiscard]] OutputPointerType get_output(const std::string& output_name) const;
        [[maybe_unused]] [[nodiscard]] ParticlePointer get_particle(std::size_t id) const;

        // Setters
        void set_gravity(const Vec3& g) { gravity_ = g; }
        void set_mass_scale_factor(double factor) { mass_scale_factor_ = factor; }
        void set_rotation(bool particle_rotation) {rotation_ = particle_rotation; }

        [[maybe_unused]] void set_time_increment(std::chrono::duration<double> dt);

        void write_restart_file(const std::string& filename) const;

        class RunFunctorBase {
        public:
            virtual ~RunFunctorBase() = default;
            virtual bool operator()() = 0;
        };

        // Functors for running a simulation until a condition is fulfilled
        class RunForTime : public RunFunctorBase {
        public:
            RunForTime(const Engine& e, std::chrono::duration<double> time) :
                engine_{e}, start_time_{engine_.get_time()}, time_to_run_{time} {}
            ~RunForTime() = default;
            void reset(std::chrono::duration<double> new_run_time) {
                start_time_ = engine_.get_time();
                time_to_run_ = new_run_time;
            }

            bool operator()() override
            {
                using namespace std::chrono_literals;
                return time_to_run_  - (engine_.get_time() - start_time_ ) > 0.1ns;
            }

        private:
            const Engine& engine_;
            std::chrono::duration<double> start_time_ ;
            std::chrono::duration<double> time_to_run_;
        };

        class KineticEnergyLess : public RunFunctorBase {
        public:
            KineticEnergyLess(const Engine& e, double kinetic_energy, std::chrono::duration<double> update_time) :
                    engine_{e}, kinetic_energy_{kinetic_energy}, start_time_(engine_.get_time()),
                    update_time_{update_time} {}
            ~KineticEnergyLess() = default;
            void set_new_value(double new_kinetic_energy) { kinetic_energy_ = new_kinetic_energy; }
            bool operator()() override {
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

        class ParticleVelocityLess : public RunFunctorBase {
        public:
            ParticleVelocityLess(const Engine& e, double max_velocity, std::chrono::duration<double> update_time) :
                    engine_{e}, max_velocity_{max_velocity}, start_time_(engine_.get_time()),
                    update_time_{update_time} {}
            ~ParticleVelocityLess() = default;
            void set_new_value(double new_max_velocity) { max_velocity_ = new_max_velocity; }
            bool operator()() override {
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

        class ObjectVelocityLess : public RunFunctorBase {
        public:
            ObjectVelocityLess(const Engine& e, double max_velocity, std::chrono::duration<double> update_time) :
                    engine_{e}, max_velocity_{max_velocity}, start_time_(engine_.get_time()),
                    update_time_{update_time} {}
            ~ObjectVelocityLess() = default;
            void set_new_value(double new_max_velocity) { max_velocity_ = new_max_velocity; }
            bool operator()() override {
                if (engine_.get_time() - start_time_ < update_time_) {
                    return true;
                }
                start_time_ += update_time_;
                double max_velocity = std::max(engine_.max_particle_velocity().second,
                        engine_.max_surface_velocity().second);
                return max_velocity > max_velocity_;
            }

        private:
            const Engine& engine_;
            double max_velocity_;

            std::chrono::duration<double> start_time_ ;
            std::chrono::duration<double> update_time_;
        };

        class SurfaceNormalForceGreater : public RunFunctorBase {
        public:
            SurfaceNormalForceGreater(Engine::SurfaceType * surface, double force):
            surface_(surface), force_(force){}
            bool operator()() override  {
                return abs(surface_->get_normal_force()) < force_;
            }

        private:
            Engine::SurfaceType* surface_;
            double force_;
        };

        class SurfaceNormalForceLess : public RunFunctorBase {
        public:
            SurfaceNormalForceLess(Engine::SurfaceType * surface, double force):
                    surface_(surface), force_(force){}
            bool operator()() override  {
                return abs(surface_->get_normal_force()) > force_;
            }

        private:
            Engine::SurfaceType* surface_;
            double force_;
        };

        class CombinedConditions {
        public:
            CombinedConditions(std::initializer_list<RunFunctorBase*> conditions) :
            conditions_(conditions){}

            bool operator()()  {
                bool cond_value = false;
                for(const auto& cond: conditions_){
                    cond_value = cond_value || cond->operator()();
                }
                return cond_value;
            }
        private:
            std::vector<RunFunctorBase*> conditions_;
        };

        class SurfaceNormalForceWithinInterval : public RunFunctorBase {
        public:
            SurfaceNormalForceWithinInterval(const Engine& engine, const SurfaceType* surface, double fmin, double fmax,
                    std::chrono::duration<double>time_interval) :
                    engine_(engine), surface_(surface), fmin_(fmin), fmax_(fmax), time_interval_(time_interval) ,
                    start_time_(engine.get_time()){}
            ~SurfaceNormalForceWithinInterval() = default;
            bool operator()() override {
                using namespace std::chrono_literals;
                // ToDo Fix that it is just force in z that is studied!!!
                double fn = abs(surface_->get_total_force().z());
                if (fn > fmin_ && fn < fmax_) {
                    time_count_ = engine_.get_time() - start_time_;
                }
                else {
                    time_count_ = 0s;
                    start_time_ = engine_.get_time();
                }
                return time_count_ < time_interval_;
            }

        private:
            const Engine& engine_;
            const SurfaceType * surface_;
            double fmin_;
            double fmax_;
            std::chrono::duration<double> time_interval_;
            std::chrono::duration<double> start_time_;
            std::chrono::duration<double> time_count_ { 0. };
        };

        class SurfaceVelocityLessThan : public RunFunctorBase {
        public:
            SurfaceVelocityLessThan(double max_velocity, SurfaceType* surface) :
            velocity_(max_velocity), surface_(surface) {}
            ~SurfaceVelocityLessThan() = default;
            bool operator()() override {
                return surface_->get_velocity().length() > velocity_;
            }
        private:
            double velocity_;
            SurfaceType* surface_;
        };

    private:
        using ContactType = Contact<ForceModel, ParticleType>;


        std::size_t object_id_counter_ { 0 };
        std::size_t collision_id_counter_ { 0 };
        std::chrono::duration<double> time_ { 0. };

        std::vector<MaterialBase*> materials_{};
        std::vector<ParticleType*> particles_{};
        std::vector<SurfaceType*> surfaces_{};
        ContactMatrix<ContactType> contacts_{};
        std::vector<OutputPointerType> outputs_{};

        CollisionDetector<ForceModel, ParticleType> collision_detector_;
        std::unique_ptr<PeriodicBCHandlerType> periodic_bc_handler_ = nullptr;


        // Settings type of private data
        Vec3 gravity_ {Vec3{0,0,0}};
        std::vector<std::pair<double, std::size_t>> viscocity_parameters_ {};
        std::chrono::duration<double> increment_{};
        double bounding_box_stretch_ { 1e-4};
        double mass_scale_factor_ { 1. };
        bool rotation_ = true;

        void do_step();

        void make_material_from_restart_data(const ParameterMap& parameters);
        void make_output_from_restart_data(const ParameterMap& parameters);
        void make_surface_from_restart_data(const ParameterMap& parameters);
        void make_particle_from_restart_data(const ParameterMap& parameters);
        void make_contact_from_restart_data(const ParameterMap& parameters);

        // Helper functions
        void move_particles();
        void move_surfaces();
        void create_contacts();
        void destroy_contacts();
        void update_contacts();
        void sum_contact_forces();
        void run_output();

        friend class Output<ForceModel, ParticleType>;
        friend class PeriodicBCHandler<ForceModel, ParticleType>;



    };

}

#include "engine.tpp"

#endif //DEMSIM_ENGINE_H
