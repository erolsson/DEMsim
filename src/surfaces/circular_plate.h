//
// Created by erolsson on 09/07/19.
//

#ifndef DEMSIM_CIRCULAR_PLATE_H
#define DEMSIM_CIRCULAR_PLATE_H

#include "surface_base.h"
#include "../utilities/vec3.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class CircularPlate : Surface<ForceModel, ParticleType> {
    public:
        CircularPlate(std::size_t id, double radius, const Vec3& normal, const Vec3& mid_point);

        ~CircularPlate() override = default;

        const Vec3& get_normal() { return normal_; }

        [[nodiscard]] double distance_to_point(const Vec3& point) const override;

        [[nodiscard]] Vec3 vector_to_point(const Vec3& point) const override;

        [[nodiscard]] Vec3 displacement_this_inc(const Vec3& position) const override;

        void move(const Vec3& distance, const Vec3& velocity) override;

        void rotate(const Vec3& position, const Vec3& rotation_vector) override;

        [[nodiscard]] std::string get_output_string() const override;

        [[nodiscard]] double get_radius() const { return radius_; }

    private:
        using Surface<ForceModel, ParticleType>::id_;
        using Surface<ForceModel, ParticleType>::displacement_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_point_;
        using Surface<ForceModel, ParticleType>::velocity_;
        using Surface<ForceModel, ParticleType>::bbox_values_;

        double radius_;
        Vec3 normal_;
        Vec3 mid_point_;

        void update_bounding_box();
    };
}

#include "circular_plate.tpp"


#endif //DEMSIM_CIRCULAR_PLATE_H
