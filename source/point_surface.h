//
// Created by erolsson on 2018-08-12.
//

#ifndef DEMSIM_POINT_SURFACE_H
#define DEMSIM_POINT_SURFACE_H

#include <algorithm>
#include <sstream>
#include <vector>
#include <utility>

#include "vec3.h"
#include "surface.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class PointSurface : public Surface<ForceModel, ParticleType> {
    public:
        PointSurface(unsigned id, std::vector<Vec3> points, bool infinite);
        ~PointSurface() override = default;

        Vec3 get_normal(const Vec3&) const override { return normal_; }
        Vec3 get_normal() const { return normal_; }     //Convenience function for avoiding to pass an arbitrary point

        double distance_to_point(const Vec3& point) const override;
        Vec3 vector_to_point(const Vec3& point) const override;
        Vec3 displacement_this_inc(const Vec3& position) const override;

        void move(const Vec3& distance, const Vec3& velocity) override;
        void rotate(const Vec3& point, const Vec3& rotation_vector) override;

        std::string output_data() const override;

        const std::vector<double>& bounding_box_values() const override;

    private:
        std::vector<Vec3> points_;
        bool infinite_;
        Vec3 normal_;

        using Surface<ForceModel, ParticleType>::id_;
        using Surface<ForceModel, ParticleType>::displacement_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_this_inc_;
        using Surface<ForceModel, ParticleType>::rotation_point_;
        using Surface<ForceModel, ParticleType>::velocity_;

        std::vector<double> bbox_values_;

        Vec3 calculate_normal() const;
        void update_bounding_box();

    };

    template<typename ForceModel, typename ParticleType>
    PointSurface<ForceModel, ParticleType>::PointSurface(unsigned id, std::vector<Vec3> points, bool infinite) :
        Surface<ForceModel, ParticleType>::Surface(id),
        points_(std::move(points)),
        infinite_(infinite),
        normal_(calculate_normal()),
        bbox_values_(6, 0)
    {
            update_bounding_box();
    }


    template<typename ForceModel, typename ParticleType>
    double PointSurface<ForceModel, ParticleType>::distance_to_point(const Vec3& point) const
    {
        return vector_to_point(point).length();
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 PointSurface<ForceModel, ParticleType>::vector_to_point(const Vec3& point) const
    {
        Vec3 v = dot_product(point-points_[0], normal_)*normal_;
        if(infinite_)
            return v;

        Vec3 plane_vector = point - v;  // Vector in the same plane as the surface going from origo to the point

        double min_distance = 1E99;
        Vec3 min_vector = Vec3(0, 0, 0);
        bool inside = true;

        for (unsigned i = 0; points_.size(); ++i){
            Vec3 vec{};
            Vec3 ps = plane_vector - points_[i];    // Vector from a corner to the point in the plane of the surfacve
            Vec3 edge{};
            if (i != points_.size() - 1)
                edge = points_[i+1] - points_[i];
            else
                edge = points_[points_.size() - 1] - points_[0];

            double l = dot_product(ps, edge.normal());

            if(l<0)
                vec = ps;
            else if(l>edge.length())
                vec = ps-edge;
            else
                vec = ps-l*edge.normal();
            l = vec.length();
            inside = inside && dot_product(cross_product(vec, edge), normal_) < 0;
            if(l < min_distance){
                min_distance = l;
                min_vector = vec+v;
            }
        }

        if(inside)
            return v;
        return min_vector;
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 PointSurface<ForceModel, ParticleType>::displacement_this_inc(const Vec3& position) const
    {
        if (rotation_this_inc_.is_zero())
            return displacement_this_inc_;
        return displacement_this_inc_ + cross_product(rotation_this_inc_, position - rotation_point_);
    }

    template<typename ForceModel, typename ParticleType>
    void PointSurface<ForceModel, ParticleType>::move(const Vec3& distance, const Vec3& velocity)
    {
        for (auto& p: points_){
            p += distance;
        }
        displacement_this_inc_ = distance;
        velocity_ = velocity;
        normal_ = calculate_normal();
        update_bounding_box();
    }

    template<typename ForceModel, typename ParticleType>
    void PointSurface<ForceModel, ParticleType>::rotate(const Vec3& point, const Vec3& rotation_vector)
    {
        for(auto& p: points_) {
            p += cross_product(rotation_vector, p-point);
        }
        rotation_this_inc_ = rotation_vector;
        rotation_point_ = point;
        normal_ = calculate_normal();
        update_bounding_box();
    }

    template<typename ForceModel, typename ParticleType>
    std::string PointSurface<ForceModel, ParticleType>::output_data() const
    {
        std::ostringstream stream;
        stream << id_;
        for(auto& p: points_) {
            stream << "," << p.x << "," << p.y << "," << p.z;
        }
        return stream.str();
    }

    template<typename ForceModel, typename ParticleType>
    Vec3 PointSurface<ForceModel, ParticleType>::calculate_normal() const
    {
        Vec3 normal = cross_product(points_[1]-points_[0], points_[2]-points_[1]);
        return normal.normalize();
    }

    template<typename ForceModel, typename ParticleType>
    const std::vector<double>& PointSurface<ForceModel, ParticleType>::bounding_box_values() const
    {
        return bbox_values_;
    }

    template<typename ForceModel, typename ParticleType>
    void PointSurface<ForceModel, ParticleType>::update_bounding_box()
    {
        auto x_cmp = [](const Vec3& v1, const Vec3& v2) -> bool {
            return v1.x < v2.x;
        };

        auto y_cmp = [](const Vec3& v1, const Vec3& v2) -> bool {
            return v1.y < v2.y;
        };

        auto z_cmp = [](const Vec3& v1, const Vec3& v2) -> bool {
            return v1.z < v2.z;
        };

        bbox_values_[0] = std::min_element(points_.begin(), points_.end(), x_cmp)->x;
        bbox_values_[1] = std::max_element(points_.begin(), points_.end(), x_cmp)->x;

        bbox_values_[2] = std::min_element(points_.begin(), points_.end(), y_cmp)->y;
        bbox_values_[3] = std::max_element(points_.begin(), points_.end(), y_cmp)->y;

        bbox_values_[4] = std::min_element(points_.begin(), points_.end(), z_cmp)->z;
        bbox_values_[5] = std::max_element(points_.begin(), points_.end(), z_cmp)->z;
    }
}


#endif //DEMSIM_POINT_SURFACE_H
