//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_ENGINE_H
#define DEMSIM_ENGINE_H

#include <vector>
#include <map>
#include <cstddef>
#include <iostream>
#include <memory>

#include "vec3.h"
#include "output.h"
#include "contact_matrix.h"
#include "material_base.h"

using std::size_t;

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Engine {
        friend class Output<ForceModel, ParticleType>;
        using ContactPointerType = std::shared_ptr<ForceModel*>;
        using SurfaceType = Surface<ContactPointerType, ParticleType>;

        std::vector<ParticleType*> particles;
        ContactMatrix<ContactPointerType> contacts;
        std::vector<MaterialBase*> materials;
        std::vector<SurfaceType*> surfaces;
        std::vector<Edge> edges;

        std::vector<Output<ForceModel, ParticleType>> outputs;
        CollisionDetector collision_detector;
        double time;
        Settings* settings;

        Vec3 gravity;
        unsigned damped_particles;
        unsigned material_counter;
        unsigned number_of_objects;

        void do_step();
        void detect_collision();
        void move_surfaces();
        void move_particles();
        void run_output();
        void clear_model();
        void check_fractured_particles();
        void updateContacts();

    public:
        Engine();
        Engine(const Engine&);
        ~Engine();

        void setup();

        template <typename T>
        void run(T);
        void run(double, bool, bool tensileTest = false);
        void run(unsigned);

        ParticleType* createParticle(Vec3, Vec3, double, MaterialBase*);
        Rectangle* createRectangle(Vec3, Vec3,Vec3, Vec3, bool inf=false);
        Edge* createEdge(Vec3, Vec3,Vec3, Vec3);
        Cylinder* createCylinder(double, Vec3, Vec3);
        CylinderPlate* createPlateCylinder(double, Vec3, Vec3, Vec3);
        Edge* createEdge(Vec3, Vec3);
        std::size_t numberOfParticles() const;
        std::size_t nOfObjects() const;

        Material* createMaterial();

        void add_surface_movement(Rectangle*, Vec3, double, double);
        void add_surface_movement(Rectangle*, double); //used for perpendicular move
        void add_surface_movement(Rectangle*, Vec3); //used for perpendicular move
        void add_surface_movement(Rectangle*, double, double);
        void add_surface_movement(Rectangle*, double, double, double);
        void add_surface_movement(Rectangle*, double, double, double, bool);
        void addSurfaceRotation(Rectangle*, const Vec3&, const Vec3&);
        void addControlableSurface(SurfaceControlData*);
        void addExpandableCylinder(Cylinder*, double);
        void addExpandableCylinder(CylinderPlate*, double);
        void addMoveableParticles(const std::vector<SphericalParticle*>&, Vec3);
        void clearMoveSchemes();
        void shakeWalls(char, double, double, double, double);
        void shakeWalls(std::vector<Rectangle*>,
                        char, double, double, double, double);

        bool createContact(SphericalParticle*, SphericalParticle*);
        bool createContact(SphericalParticle*, Surface*);
        bool createContact(SphericalParticle*, Edge*);
        bool destroyContact(SphericalParticle*, SphericalParticle*);
        bool destroyContact(SphericalParticle*, Surface*);
        bool destroyContact(SphericalParticle*, Edge*);
        SphericalParticle* getParticle(size_t) const;
        Surface* getSurface(size_t) const;
        double getTime() const;
        double totalKineticEnergy() const;
        double avgNumberOfContacts() const;
        double getStickRatio() const;
        double getTotalWornVolume() const;
        double getMinFrictionCoefficient() const;
        double relativeDensity() const;
        void removeWall(Surface*);
        void setGravity(Vec3);
        double velLimitedParticles() const;
        Settings* getSettings();
        Output* createOutput(std::string);
        Output* createOutput(std::string, std::vector<SphericalParticle*>&);
        Output* createOutput(std::string, double);
        Output* createOutput(std::string, double, double);
        void removeOutput(Output*);
        unsigned noOfObjects() const;
        void restart(std::string,unsigned);
    };

    template<typename T>
    void Engine::run(T pred) {
        unsigned counter = 0;
        while(pred()) {
            doStep();                   //Do one step
            time+=settings->increment;  //increase time
            ++counter;
            if(counter % 1000 == 0) {
                std::cout << "Fraction of Sticking Contacts = "
                          << getStickRatio() << std::endl;
                std::cout << "Total worn volume = "
                          << getTotalWornVolume() << std::endl;
                std::cout << "The minimum friction coefficient is "
                          << getMinFrictionCoefficient() << std::endl;
                std::cout << "Simulation Time is " << getTime() << std::endl;

            }
        }
    }

}

#endif //DEMSIM_ENGINE_H
