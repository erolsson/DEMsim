//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_BOUNDING_BOX_PROJECTION_H
#define DEMSIM_BOUNDING_BOX_PROJECTION_H

#include <array>

#include "surface.h"

namespace DEM {
    template<typename ForceModel, typename particleType> class BoundingBox;

    template<typename ForceModel, typename ParticleType>
    class BoundingBoxProjection {
        using SurfaceType = Surface<ForceModel, ParticleType>;

    public:
        double value;

    private:
        char position_char_;
        char axis_;
        int sign;

        BoundingBox<ForceModel, ParticleType>* bbox;
        std::array<std::size_t*, 4> other_indices;
        std::size_t index;
    };
}

#endif //DEMSIM_BOUNDING_BOX_PROJECTION_H
