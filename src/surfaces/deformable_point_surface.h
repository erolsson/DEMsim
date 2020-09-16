//
// Created by erolsson on 14/09/2020.
//

#ifndef DEMSIM_DEFORMABLE_POINT_SURFACE_H
#define DEMSIM_DEFORMABLE_POINT_SURFACE_H

#include "point_surface.h"
#include "../utilities/vec3.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class DeformablePointSurface : public PointSurface<ForceModel, ParticleType> {
    public:
        DeformablePointSurface(std::size_t id, const std::vector<Vec3>& points, bool infinite, const std::string& name,
                     bool adhesive=true, std::size_t collision_id=0);
        explicit DeformablePointSurface(const ParameterMap& parameters);
        ~DeformablePointSurface() override = default;
        void deform(const std::vector<Vec3>& nodal_displacements);
        [[nodiscard]] Vec3 get_displacement_this_increment(const Vec3& position) const override;

        [[nodiscard]] std::string get_output_string() const override;
        [[nodiscard]] std::string restart_data() const override;

        using PointSurface<ForceModel, ParticleType>::get_normal;
        void set_in_plane_strain_rates(double ex, double ey)
        {
            strain_x_ = ex;
            strain_y_ = ey;
        }

        void deform(std::chrono::duration<double> time_increment);

    private:
        std::vector<Vec3> nodal_displacements_;
        double strain_x_ { 0. };
        double strain_y_ { 0. };
        using PointSurface<ForceModel, ParticleType>::points_;
        using PointSurface<ForceModel, ParticleType>::displacement_this_inc_;
        using PointSurface<ForceModel, ParticleType>::velocity_;
        using PointSurface<ForceModel, ParticleType>::normal_;
        using PointSurface<ForceModel, ParticleType>::calculate_normal;
        using PointSurface<ForceModel, ParticleType>::update_bounding_box;
    };
}

#include "deformable_point_surface.tpp"


#endif //DEMSIM_DEFORMABLE_POINT_SURFACE_H
