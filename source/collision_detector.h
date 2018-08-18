//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_COLLISION_DETECTOR_H
#define DEMSIM_COLLISION_DETECTOR_H

#include <vector>
#include <omp.h>

#include "bounding_box.h"
#include "bounding_box_projection.h"
#include "contact_matrix.h"
#include "point_surface.h"

namespace DEM {
    template <typename ForceModel, typename ParticleType>
    class CollisionDetector {
        using BoundingBoxType = BoundingBox<ForceModel, ParticleType>;
        using BoundingBoxProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;

        using ContactPointerType = std::shared_ptr<Contact<ForceModel, ParticleType> >;

    public:
        using CollisionPair = std::pair<BoundingBoxType*, BoundingBoxType*>;
        CollisionDetector(const std::vector<ParticleType*>& particles,
                const std::vector<PointSurface<ForceModel, ParticleType>*>& point_surfaces);

        void setup();
        std::vector<CollisionPair> do_check();  //Not const due to re-ordering of the proj vectors

    private:
        std::vector<BoundingBox<ForceModel, ParticleType> > bounding_boxes_;
        std::vector<BoundingBoxProjection<ForceModel, ParticleType>* > xproj_;
        std::vector<BoundingBoxProjection<ForceModel, ParticleType>* > yproj_;
        std::vector<BoundingBoxProjection<ForceModel, ParticleType>* > zproj_;

        const std::vector<ParticleType*>& particles_;
        const std::vector<PointSurface<ForceModel, ParticleType>*>& point_surfaces_;

        ContactMatrix<CollisionPair> active_collisions_;
        //BoolMatrix activeCollisions;
        void update_bounding_boxes();


        // void check_possible_collision_pairs(std::vector<CollisionPair>&);
        // void check_bounding_box_vector(std::vector<BoundingBoxProjectionType*>&);
        // bool overlapping(BBProjection*, BBProjection*) const;
        // void destroyContact(BBProjection*, BBProjection*);
        // void createContact(BBProjection*, BBProjection*);
        // inline bool checkOtherAxes(BBProjection&, BBProjection&) const;
    };

    template<typename ForceModel, typename ParticleType>
    std::vector<typename CollisionDetector<ForceModel, ParticleType>::CollisionPair>
            CollisionDetector<ForceModel, ParticleType>::do_check()
    {
        update_bounding_boxes();
        return std::vector<CollisionDetector::CollisionPair>();
    }

    template<typename ForceModel, typename ParticleType>
    CollisionDetector<ForceModel, ParticleType>::CollisionDetector(const std::vector<ParticleType*>& particles,
            const std::vector<PointSurface<ForceModel, ParticleType>*>& point_surfaces) :
            particles_(particles), point_surfaces_(point_surfaces)
    {
        //Empty constructor
    }

    template<typename ForceModel, typename ParticleType>
    void CollisionDetector<ForceModel, ParticleType>::setup()
    {
        for(std::size_t i = 0; i != particles_.size(); ++i){
            bounding_boxes_.emplace_back(particles_[i], i);
        }

        for(std::size_t i = 0; i != point_surfaces_.size(); ++i){
            bounding_boxes_.emplace_back(point_surfaces_[i], i);
        }


    }

    template<typename ForceModel, typename ParticleType>
    void CollisionDetector<ForceModel, ParticleType>::update_bounding_boxes()
    {
        #pragma omp parallel for
        for(std::size_t i = 0; i < bounding_boxes_.size(); ++i){
            bounding_boxes_[i].update();
        }
    }




}

#endif //DEMSIM_COLLISION_DETECTOR_H
