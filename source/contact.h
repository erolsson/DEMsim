//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_CONTACT_H
#define DEMSIM_CONTACT_H

#include <memory>

#include "surface.h"
#include "edge.h"


namespace DEM{
    class Edge;
    template<typename ForceModel, typename ParticleType> class Surface;

    template<typename ForceModel, typename ParticleType>
    class Contact {
        using ContactPointerType = std::shared_ptr<Contact<ForceModel, ParticleType> >;
        using SurfaceType = Surface<ContactPointerType, ParticleType>;

    public:
        //Constructors
        Contact(ParticleType*, ParticleType*, double);
        Contact(ParticleType*, SurfaceType*,  double);
        Contact(ParticleType*, Edge*,     double);
        ~Contact() = default;                        //No dynamic allocated memory

        void update();
        Vec3 position() const;

        Vec3 get_normal_force() const { return force_model_.get_normal_force()*get_normal(); }
        Vec3 get_tangential_force() const { return force_model_.get_tangential_force(); }
        Vec3 get_torque(const Vec3& point) const { return cross_product(position() - point, get_tangential_force()); }

        std::pair<ParticleType*, ParticleType*> get_particles() const  { return std::make_pair(p1_, p2_); }
        const SurfaceType* get_surface() const { return surface_; }
        const Edge* get_edge() const { return edge_; }

        bool active() const { return force_model_.acitve();  }
        double get_overlap() const { return force_model_.get_overlap(); }
        double get_contact_radius() const { return force_model_.get_contact_area();}
        const Vec3& get_normal() const { return normal_; }

    private:
        ParticleType* const p1_;
        ParticleType* const p2_;
        SurfaceType*  const surface_;
        Edge*         const edge_;

        double r2_;                       // r2 - distance between particles or particle plane is overlap
        Vec3 normal_;                     //Contact plane normal, same direction as normal_force_

        ForceModel force_model_;
        bool const particle_contact_;

        Vec3 calculate_distance_vector() const;
        Vec3 calculate_tangential_displacement_this_inc() const;
    };



    template<typename ForceModel, typename ParticleType>
    Contact<ForceModel, ParticleType>::Contact(ParticleType* p1, ParticleType* p2, double increment) :
            p1_(p1),
            p2_(p2),
            surface_(nullptr),
            edge_(nullptr),
            r2_(p1->get_radius() + p2->get_radius()),
            force_model_(p1, p2, increment),
            particle_contact_(true)
    {
        normal_ = calculate_distance_vector().normal();
    }


    template<typename ForceModel, typename ParticleType>
    Contact<ForceModel, ParticleType>::Contact(ParticleType* p, SurfaceType* s,  double increment) :
            p1_(p),
            p2_(nullptr),
            surface_(s),
            edge_(nullptr),
            r2_(p->get_radius()),
            force_model_(p, s, increment),
            particle_contact_(false)

    {
        normal_ = calculate_distance_vector().normal();
    }


    template<typename ForceModel, typename ParticleType>
    Contact<ForceModel, ParticleType>::Contact(ParticleType* p, Edge* e, double increment) :
            p1_(p),
            p2_(nullptr),
            surface_(nullptr),
            edge_(e),
            r2_(p->get_radius()),
            force_model_(p, e, increment),
            particle_contact_(false)

    {
        normal_ = calculate_distance_vector().normal();
    }


    template<typename ForceModel, typename ParticleType>
    void Contact<ForceModel, ParticleType>::update()
    {
        auto distance_vector = calculate_distance_vector();
        double h = r2_ - distance_vector.length();
        normal_ = distance_vector.normalize();
        Vec3 dt = calculate_tangential_displacement_this_inc();
        force_model_.update(h, dt);
    }


    template<typename ForceModel, typename ParticleType>
    Vec3 Contact<ForceModel, ParticleType>::position() const
    {
        return Vec3();
    }


    template<typename ForceModel, typename ParticleType>
    Vec3 Contact<ForceModel, ParticleType>::calculate_distance_vector() const
    {
        if (particle_contact_)
            return p1_->get_position() - p2_->get_position();
        else if ( surface_ !=nullptr)
            return surface_->vector_to_point(p1_->get_position());
        else
            return edge_->vector_to_point(p1_->get_position());
    }


    template<typename ForceModel, typename ParticleType>
    Vec3 Contact<ForceModel, ParticleType>::calculate_tangential_displacement_this_inc() const
    {
        Vec3 r1 = -normal_*(p1_->get_radius() - get_overlap()/2);
        Vec3 u1 = p1_->get_displacement_this_increment() + cross_product(p1_->get_rotation_this_increment(), r1);
        Vec3 u2 = Vec3(0., 0. , 0.);
        if (particle_contact_) {
            Vec3 r2 = normal_*(p2_->get_radius() - get_overlap()/2);
            u2 = p2_->get_displacement_this_increment() + cross_product(p2_->get_rotation_this_increment(), r2);
        }
        else if (surface_ !=nullptr) {
            u2 = surface_->displacement_this_inc(position());
        }
        return (u1 - u2) - dot_product((u1 - u2), normal_)*normal_;
    }

}




#endif //DEMSIM_CONTACT_H
