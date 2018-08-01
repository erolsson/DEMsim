//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_COLLISION_DETECTOR_H
#define DEMSIM_COLLISION_DETECTOR_H

#include <vector>

namespace DEM {
    class CollisionDetector {
        std::vector<BBProjection*> xproj;
        std::vector<BBProjection*> yproj;
        std::vector<BBProjection*> zproj;
        std::vector<BoundingBox*> boxes;
        std::vector<Particle*>& particles;
        std::vector<Surface*>& surfaces;
        std::vector<Edge*>& edges;
        Engine& engine;
        ContactMatrix<BBPair> activeCollisions;
        //BoolMatrix activeCollisions;
        void generateBoundingBoxes();
        void updateBoundingBoxes(double k);

        double overlap(BBProjection*, BBProjection*) const;
        double overlap(Particle*, Particle*) const;

        void checkPossibleCollisionPairs(std::vector<cType>&);
        void checkBoundingBoxVector(std::vector<BBProjection*>&);
        bool overlapping(BBProjection*, BBProjection*) const;
        double distanceParticleWall(const Particle*, const Surface*) const;
        bool insideTriangle(const Vec3&, const Surface*) const;
        void destroyContact(BBProjection*, BBProjection*);
        void createContact(BBProjection*, BBProjection*);
        inline bool checkOtherAxes(BBProjection&, BBProjection&) const;
        void checkSurfaces();
        void checkEdges();
    };
}
#endif //DEMSIM_COLLISION_DETECTOR_H
