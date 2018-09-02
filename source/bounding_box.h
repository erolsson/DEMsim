//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_BOUNDING_BOX_H
#define DEMSIM_BOUNDING_BOX_H

#include <cmath>

#include "bounding_box_projection.h"
#include "cylinder.h"
#include "point_surface.h"
#include "surface.h"

namespace DEM {
    template <typename ForceModel, typename ParticleType>
    class BoundingBox {
    using BProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;
    using PointSurfaceType = PointSurface<ForceModel, ParticleType>;
    using SurfaceType = Surface<ForceModel, ParticleType>;
    using CylinderType = Cylinder<ForceModel, ParticleType>;
    public:
        BoundingBox(ParticleType* particle, std::size_t index);
        BoundingBox(SurfaceType* surface,  std::size_t index);
        BoundingBox(CylinderType* cylinder,  std::size_t index, bool);  //Special bounding box for inward cylinders

        // Copy constructor and assignment operator needed due to pointers between different projection vectors
        // which becomes invalid when different boundingboxes are re-allocated due to vector-over-capacity
        // These constructors repairs the pointers
        BoundingBox(const BoundingBox& rhs);
        BoundingBox& operator=(const BoundingBox& rhs);


        std::size_t get_id() const;
        void update();
        void set_stretch(double stretch) {stretch_ = stretch; }

        ParticleType* get_particle() const { return particle_;}
        SurfaceType* get_surface() const { return surface_;}

        BProjectionType bx;
        BProjectionType ex;
        BProjectionType by;
        BProjectionType ey;
        BProjectionType bz;
        BProjectionType ez;

    private:
        ParticleType* particle_;
        SurfaceType* surface_;
        double stretch_{ 0. };

        // Function pointer to the update function, set in the construction of the bounding box
        // In doing so, if statements is avoided at each update
        void (BoundingBox<ForceModel, ParticleType>::*update_function)();
        void particle_update();
        void surface_update();
        void setup_projections();
    };
}

#include "bounding_box.h"
#endif //DEMSIM_BOUNDING_BOX_H
