//
// Created by erolsson on 2018-08-12.
//

#ifndef DEMSIM_POINT_SURFACE_H
#define DEMSIM_POINT_SURFACE_H

#include <vector>

#include "vec3.h"
#include "surface.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class PointSurface : public Surface<ForceModel, ParticleType> {
    public:
        PointSurface(const std::vector<Vec3>&, bool);

    private:
        std::vector<Vec3> points_;
        bool infinte_;

    };
}

#endif //DEMSIM_POINT_SURFACE_H
