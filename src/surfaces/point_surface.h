//
// Created by erolsson on 2018-08-12.
//

#ifndef DEMSIM_POINT_SURFACE_H
#define DEMSIM_POINT_SURFACE_H

#include <algorithm>
#include <array>
#include <sstream>
#include <vector>
#include <utility>
#include "../utilities/vec3.h"
#include "surface_base.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class PointSurface : public Surface<ForceModel, ParticleType> {
    public:
        PointSurface(std::size_t id, std::vector<Vec3> points, bool infinite, bool adhesive=true);
        ~PointSurface() override = default;

        [[nodiscard]] Vec3 get_normal(const Vec3&) const override { return normal_; }
        [[nodiscard]] Vec3 get_normal() const { return normal_; }     //Convenience function for avoiding to pass an arbitrary point

        [[nodiscard]] double distance_to_point(const Vec3& point) const override;
        [[nodiscard]] Vec3 vector_to_point(const Vec3& point) const override;
        [[nodiscard]] Vec3 displacement_this_inc(const Vec3& position) const override;

        void move(const Vec3& distance, const Vec3& velocity) override;
        void rotate(const Vec3& point, const Vec3& rotation_vector) override;

        [[nodiscard]] std::string get_output_string() const override;
        [[nodiscard]] bool adhesive() const { return adhesive_; }

    private:
        std::vector<Vec3> points_;
        bool infinite_;
        bool adhesive_ {false};
        Vec3 normal_;

        using Surface<ForceModel, ParticleType>::id_;
        using Surface<ForceModel, ParticleType>::displacement_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_point_;
        using Surface<ForceModel, ParticleType>::velocity_;
        using Surface<ForceModel, ParticleType>::bbox_values_;

        [[nodiscard]] Vec3 calculate_normal() const;
        void update_bounding_box() override;
    };

}

#include "point_surface.tpp"

#endif //DEMSIM_POINT_SURFACE_H
