//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_BOUNDING_BOX_H
#define DEMSIM_BOUNDING_BOX_H

#include <cmath>

#include "bounding_box_projection.h"
#include "../../surfaces/cylinder.h"
#include "../../surfaces/point_surface.h"
#include "../../surfaces/surface_base.h"

namespace DEM {
    template <typename ForceModel, typename ParticleType>
    class BoundingBox {
    using BProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;
    using PointSurfaceType = PointSurface<ForceModel, ParticleType>;
    using SurfaceType = Surface<ForceModel, ParticleType>;
    using CylinderType = Cylinder<ForceModel, ParticleType>;

    public:
        BoundingBox(ParticleType* particle, std::size_t index, double stretch);
        BoundingBox(SurfaceType* surface,  std::size_t index, double stretch);
        BoundingBox(CylinderType* cylinder, std::size_t index, double stretch, bool); //Special bounding box for inward cylinders

        // Copy constructor and assignment operator needed due to pointers between different projection vectors
        // which becomes invalid when different bounding boxes are re-allocated due to vector-over-capacity
        // These constructors repairs the pointers
        BoundingBox(const BoundingBox& rhs) = delete;
        BoundingBox& operator=(const BoundingBox& rhs) = delete;

        BoundingBox(BoundingBox&& rhs) noexcept;
        BoundingBox& operator=(BoundingBox&& rhs) noexcept;

        [[nodiscard]] std::size_t get_collision_id() const;
        [[nodiscard]] std::size_t get_object_id() const;
        void update();
        void set_stretch(double stretch) {stretch_ = stretch; }

        ParticleType* get_particle() const { return particle_;}
        SurfaceType* get_surface() const { return surface_;}

        std::vector<BProjectionType> bounding_box_projections;

    private:
        ParticleType* particle_;
        SurfaceType* surface_;
        double stretch_;

        void (BoundingBox<ForceModel, ParticleType>::*update_function)();
        void particle_update();
        void surface_update();
    };

    template<typename ForceModel, typename ParticleType>
    bool operator==(const BoundingBox<ForceModel, ParticleType> b1, const BoundingBox<ForceModel, ParticleType> b2) {
        return b1.get_collision_id() == b2.get_collision_id();
    }

    template<typename ForceModel, typename ParticleType>
    bool operator!=(const BoundingBox<ForceModel, ParticleType> b1, const BoundingBox<ForceModel, ParticleType> b2) {
        return !(b1 == b2);
    }
}

#include "bounding_box.tpp"
#endif //DEMSIM_BOUNDING_BOX_H
