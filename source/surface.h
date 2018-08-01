//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_SURFACE_H
#define DEMSIM_SURFACE_H

#include <vector>
#include <map>
#include <memory>

#include "vec3.h"
#include "contact_vector.h"



namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Contact;

    class Amplitude;

    template<typename ForceModel, typename ParticleType>
    class Surface {
        using ContactPointerType = std::shared_ptr<Contact<ForceModel, ParticleType>>;

    public:
        bool no_stiffness{false};
        bool no_fracture{false};
        bool tensile{true};
        double adhesive_multiplier{1.0};
        double mass{0.};
        double max_velocity{0.};
        bool force_control{false};

        Amplitude* fx{nullptr};
        Amplitude* fy{nullptr};
        Amplitude* fz{nullptr};

        explicit Surface(unsigned);

        virtual ~Surface() = default;

        unsigned get_id() const { return id_; }

        virtual Vec3 get_normal(const Vec3&) const = 0;

        // virtual double get_curvature_radius() const = 0;  //  ToDo Implement later
        virtual double distance_to_point(const Vec3&) const = 0;

        virtual Vec3 vector_to_point(const Vec3&) const = 0;

        virtual Vec3 displacement_this_inc(const Vec3&) const = 0;

        virtual void move(const Vec3&, const Vec3&) = 0;

        virtual void rotate(const Vec3&, const Vec3&) = 0;

        virtual std::string output_data() const = 0;

        const Vec3& get_velocity() const { return velocity_; }

        const Vec3& get_acceleration() const { return acceleration_; }

        void set_velocity(const Vec3& v) { velocity_ = v; }

        void set_acceleration(const Vec3& a) { acceleration_ = a; }

        Vec3 get_tangential_displacement_this_inc(const Vec3& p) const
        {
            return displacement_this_inc(p)-dot_product(displacement_this_inc(p), get_normal(p))*get_normal(p);
        }

        void rest()
        {
            velocity_ *= 0;
            rotation_this_inc_ *= 0;
            rotation_point_ *= 0;
        }

        std::vector<ParticleType*> get_contacting_particles() const;

        double get_normal_force() const;

        Vec3 get_tangential_force() const;

        Vec3 get_total_force() const;

        void add_contact(ContactPointerType, std::size_t);

        void remove_contact(std::size_t);

    protected:
        Vec3 velocity_{Vec3(0, 0, 0)};
        Vec3 acceleration_{Vec3(0, 0, 0)};
        Vec3 displacement_this_inc_{Vec3(0, 0, 0)};

        // Might be updated to std::vector<std::pair<Vec3, Vec3> >
        // to allow for multiple rotations
        Vec3 rotation_this_inc_{Vec3(0, 0, 0)};
        Vec3 rotation_point_{Vec3(0, 0, 0)};

    private:
        unsigned id_;
        ContactVector<ContactPointerType> contacts_{ContactVector<ContactPointerType>()};
    };

    template<typename ForceModel, typename ParticleType>
    Surface<ForceModel, ParticleType>::Surface(unsigned id)
            :
            id_(id)
    {
        //Empty constructor
    }

    template<typename ForceModel, typename ParticleType>
    std::vector<ParticleType*> Surface<ForceModel, ParticleType>::get_contacting_particles() const
    {
        std::vector<ParticleType*> contacting_particles;
        for (auto const& c : contacts_) {
            //I'm not proud over the cast but it works, at least complies
            //    if(c->active())
            contacting_particles.push_back(const_cast<ParticleType*>(c->get_particles().first));
        }
        return contacting_particles;
    }

    template<typename ForceModel, typename ParticleType>
    double Surface<ForceModel, ParticleType>::get_normal_force() const
    {
        double normal_force = 0.0;
        for (auto const& c : contacts_) {
            normal_force += dot_product(c->get_normal_force(), get_normal(c->get_position()));
        }
        return normal_force;
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Surface<ForceModel, ParticleType>::get_tangential_force() const
    {
        Vec3 tangential_force = Vec3(0, 0, 0);
        for (auto const& c : contacts_) {
            tangential_force += c->get_tangential_force();
        }
        return -1*tangential_force;
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 Surface<ForceModel, ParticleType>::get_total_force() const
    {
        Vec3 force = Vec3(0, 0, 0);
        for (auto const& c : contacts_) {
            force += c->get_normal_force()+c->get_tangential_force();
        }
        return force;
    }

    template<typename ForceModel, typename ParticleType>
    void Surface<ForceModel, ParticleType>::add_contact(ContactPointerType contact, size_t pos)
    {
        contacts_.insert(pos, contact);
    }

    template<typename ForceModel, typename ParticleType>
    void Surface<ForceModel, ParticleType>::remove_contact(size_t pos)
    {
        contacts_.erase(pos);
    }

}
#endif //DEMSIM_SURFACE_H
