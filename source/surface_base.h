//
// Created by erolsson on 2018-07-27.
//

#ifndef DEMSIM_SURFACE_H
#define DEMSIM_SURFACE_H

#include <array>
#include <map>
#include <memory>
#include <vector>

#include "amplitude.h"
#include "contact_vector.h"
#include "contact_matrix.h"
#include "vec3.h"


namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Contact;

    template<typename ForceModel, typename ParticleType>
    class Engine;

    template<typename ForceModel, typename ParticleType>
    class Surface {
        using ContactType = Contact<ForceModel, ParticleType>;
        using ContactPointerType = typename ContactMatrix<ContactType>::PointerType;
        using ForceAmpPtr = std::shared_ptr<Amplitude<Engine<ForceModel, ParticleType>>>;

    public:
        double mass{0.};
        double max_velocity{0.};
        bool force_control{ false };

        explicit Surface(std::size_t id);
        virtual ~Surface() = default;
        std::size_t get_id() const { return id_; }
        virtual Vec3 get_normal(const Vec3& position) const = 0;
        // virtual double get_curvature_radius() const = 0;  //  ToDo Implement later
        virtual double distance_to_point(const Vec3& point) const = 0;
        virtual Vec3 vector_to_point(const Vec3& point) const = 0;
        virtual Vec3 displacement_this_inc(const Vec3& position) const = 0;

        virtual void move(const Vec3& distance, const Vec3& velocity) = 0;
        virtual void rotate(const Vec3& position, const Vec3& rotation_vector) = 0;

        virtual std::string get_output_string() const = 0;
        virtual const std::array<double, 6>& get_bounding_box_values() const { return bbox_values_; };
        const Vec3& get_velocity() const { return velocity_; }
        const Vec3& get_acceleration() const { return acceleration_; }

        void set_velocity(const Vec3& v) { velocity_ = v; }
        void set_acceleration(const Vec3& a) { acceleration_ = a; }

        Vec3 get_tangential_displacement_this_inc(const Vec3& point) const;
        void rest();

        std::vector<ParticleType*> get_contacting_particles() const;
        double get_normal_force() const;
        Vec3 get_tangential_force() const;
        Vec3 get_total_force() const;
        void add_contact(ContactPointerType contact, std::size_t index_of_other_object);
        void remove_contact(std::size_t index_of_other_object);

        const std::array<ForceAmpPtr, 3>& get_applied_forces() const { return force_control_amplitudes_; }
        void set_force_amplitude(ForceAmpPtr amplitude, char direction);
    protected:
        Vec3 velocity_{ Vec3(0, 0, 0) };
        Vec3 acceleration_{ Vec3(0, 0, 0) };
        Vec3 displacement_this_inc_{ Vec3(0, 0, 0) };

        // Might be updated to std::vector<std::pair<Vec3, Vec3> >
        // to allow for multiple rotations
        Vec3 rotation_this_inc_{ Vec3(0, 0, 0) };
        Vec3 rotation_point_{ Vec3(0, 0, 0) };

        std::size_t id_;

        std::array<double, 6> bbox_values_{0, 0, 0, 0, 0, 0};

        std::array<ForceAmpPtr, 3> force_control_amplitudes_ = {nullptr, nullptr, nullptr};

        virtual void update_bounding_box() = 0;

    private:
        ContactVector<ContactPointerType> contacts_{ContactVector<ContactPointerType>()};
    };
}

#include "surface_base.tpp"
#endif //DEMSIM_SURFACE_H
