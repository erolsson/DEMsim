//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_ENGINE_H
#define DEMSIM_ENGINE_H

#include <iostream>
#include <vector>

#include "collision_detector.h"
#include "contact.h"
#include "contact_matrix.h"
#include "cylinder.h"
#include "material_base.h"
#include "point_surface.h"
#include "settings.h"
#include "surface.h"
#include "vec3.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Engine {
    public:
        using ParticlePointer = ParticleType*;
        using PointSurfacePointer = PointSurface<ForceModel, ParticleType>*;
        using CylinderPointer = Cylinder<ForceModel, ParticleType>*;

        Engine();

        void setup();
        template<typename Condition>
        void run(const Condition& condition);

        //Object creation functions
        template<typename MaterialType>
        MaterialType* create_material(double density);

        ParticlePointer create_particle(double radius, const Vec3& position, const Vec3& velocity,
                                      MaterialBase* material);

        PointSurfacePointer create_point_surface(const std::vector<Vec3>& points, bool infinite);

        CylinderPointer create_cylinder(double radius, const Vec3& axis, const Vec3& base_point, double length,
                                        bool inward=true, bool infinite=false);
        // Getters
        Settings* get_settings() { return  &settings_; }
        double get_time() const { return time_; }

        // Functors for running a simulation until a condition is fulfilled
        class RunForTime {
        public:
            RunForTime(const Engine& e, double time) :
                engine_{e}, start_time_{engine_.get_time()}, time_to_run_{time} {}

            void reset(double new_run_time) {
                start_time_ = engine_.get_time();
                time_to_run_ = new_run_time;
            }
            bool operator()() const { return (engine_.get_time() - start_time_) < (time_to_run_); }
        private:

            const Engine& engine_;
            double start_time_ ;
            double time_to_run_;
        };


    private:
        using ContactType = Contact<ForceModel, ParticleType>;
        using SurfaceType = Surface<ForceModel, ParticleType>;

        std::size_t number_of_objects_{ 0 };
        double time_ { 0. };
        Settings settings_ {};

        std::vector<MaterialBase*> materials_{};
        std::vector<ParticleType*> particles_{};
        std::vector<SurfaceType*> surfaces_{};
        ContactMatrix<ContactType> contacts_{};

        CollisionDetector<ForceModel, ParticleType> collision_detector_;

        void do_step();

        // Helper functions
        void move_particles();
        void update_contacts();
    };

    // Functors for different running conditions

}

#include "engine.tpp"

#endif //DEMSIM_ENGINE_H
