//
// Created by erolsson on 2018-07-29.
//

#include "edge.h"
#include "vec3.h"

Vec3 DEM::Edge::vector_to_point(const Vec3& p) const
{
    Vec3 v = p2_ - p1_;
    Vec3 pp = p - p1_;
    double l = dot_product(pp, v.normal());

    if (l < 0)
        return pp;
    else if (l>v.length())
        return p2_ - p;
    return pp-l*v.normal();
}
