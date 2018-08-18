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
        using CollisionPair = std::pair<BoundingBoxType*, BoundingBoxType*>;
        using ContactPointerType = std::shared_ptr<Contact<ForceModel, ParticleType> >;

    public:
        CollisionDetector(const std::vector<ParticleType*>&,
                const std::vector<PointSurface<ForceModel, ParticleType>*>&);

        void setup();
        std::vector<CollisionPair> do_check();  //Not const due to re-ordering of the proj vectors

    private:
        std::vector<BoundingBox<ForceModel, ParticleType> > bounding_boxes_;
        std::vector<BoundingBoxProjection<ForceModel, ParticleType>* > xproj_;
        std::vector<BoundingBoxProjection<ForceModel, ParticleType>* > yproj_;
        std::vector<BoundingBoxProjection<ForceModel, ParticleType>* > zproj_;

        std::vector<ParticleType*>& particles_;
        std::vector<PointSurface<ForceModel, ParticleType>*>& point_surfaces_;

        ContactMatrix<CollisionPair> active_collisions_;
        //BoolMatrix activeCollisions;
        void generate_boundingBoxes();
        void update_bounding_boxes();


        void check_possible_collision_pairs(std::vector<CollisionPair>&);
        void check_bounding_box_vector(std::vector<BoundingBoxProjectionType*>&);
        bool overlapping(BBProjection*, BBProjection*) const;
        void destroyContact(BBProjection*, BBProjection*);
        void createContact(BBProjection*, BBProjection*);
        inline bool checkOtherAxes(BBProjection&, BBProjection&) const;
    };

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
        for(const auto& iter: particles_){
            bounding_boxes_.emplace_back(*iter);
        }

        for(const auto& iter: point_surfaces_){
            bounding_boxes_.emplace_back(*iter);
        }


    }

    template<typename ForceModel, typename ParticleType>
    void CollisionDetector<ForceModel, ParticleType>::update_bounding_boxes()
    {
        #pragma omp parallel for
        for(auto iter=bounding_boxes_.begin(); iter!=bounding_boxes_.end(); ++iter){
            iter->update();
        }
    }


}

#endif //DEMSIM_COLLISION_DETECTOR_H
