//
// Created by erolsson on 2018-07-29.
//

#ifndef DEMSIM_EDGE_H
#define DEMSIM_EDGE_H

#include "vec3.h"

namespace DEM {
    class Edge {
        Vec3 p1_;
        Vec3 p2_;
        unsigned id_;

    public:
        Edge(unsigned id, const Vec3& p1, const Vec3& p2) : p1_(p1), p2_(p2), id_(id) {}
        double distance_to_point(const Vec3& p) const { return vector_to_point(p).length();}
        Vec3 vector_to_point(const Vec3& p) const;
        unsigned get_id() const{ return id_;}
    };
}

#endif //DEMSIM_EDGE_H
