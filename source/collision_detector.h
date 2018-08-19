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
    template <typename ForceModel, typename ParticleType> class Contact;
    template <typename ForceModel, typename ParticleType>
    class CollisionDetector {
    public:

        using BoundingBoxType = BoundingBox<ForceModel, ParticleType>;
        using BoundingBoxProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;
        using CollisionPair = std::pair<BoundingBoxType*, BoundingBoxType*>;
        CollisionDetector(const std::vector<ParticleType*>& particles,
                const std::vector<PointSurface<ForceModel, ParticleType>*>& point_surfaces);

        void setup();
        void do_check();  //Not const due to re-ordering of the proj vectors
        std::vector<CollisionPair> contacts_to_create() const;
        std::vector<CollisionPair> contacts_to_destroy() const;

    private:
        using ContactPointerType = std::shared_ptr<Contact<ForceModel, ParticleType> >;
        std::vector<BoundingBox<ForceModel, ParticleType> > bounding_boxes_;
        std::vector<BoundingBoxProjectionType*> xproj_;
        std::vector<BoundingBoxProjectionType*> yproj_;
        std::vector<BoundingBoxProjectionType*> zproj_;

        const std::vector<ParticleType*>& particles_;
        const std::vector<PointSurface<ForceModel, ParticleType>*>& point_surfaces_;

        ContactMatrix<ContactPointerType>& contacts_;
        //BoolMatrix activeCollisions;
        void update_bounding_boxes();


        // void check_possible_collision_pairs(std::vector<CollisionPair>&);
        void check_bounding_box_vector(std::vector<BoundingBoxProjectionType*>& vector);
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
        check_bounding_box_vector(xproj_);
        check_bounding_box_vector(yproj_);
        check_bounding_box_vector(zproj_);
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

        for(auto& bounding_box: bounding_boxes_){
            xproj_.push_back(&bounding_box.bx);
            xproj_.push_back(&bounding_box.ex);

            yproj_.push_back(&bounding_box.by);
            yproj_.push_back(&bounding_box.ey);

            zproj_.push_back(&bounding_box.bz);
            zproj_.push_back(&bounding_box.ez);
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

    template<typename ForceModel, typename ParticleType>
    void CollisionDetector<ForceModel, ParticleType>::check_bounding_box_vector(
            std::vector<CollisionDetector::BoundingBoxProjectionType*>& vector)
    {
        for(unsigned i = 0; i != vector.size(); ++i){
            unsigned j = i;
            while(j != 0 && (vector[j-1]->val() > vector[j]->val())){
                BoundingBoxProjectionType* BBm = vector[j];
                BoundingBoxProjectionType* BBn = vector[j-1];
                //depending on de beginnings and endings of the swapping
                // remove or add contact
                if(BBm->position_char_() == 'e' && BBn->position_char_() == 'b')
                    destroyContact(BBm, BBn);
                if(BBm->position_char_() == 'b' && BBn->position_char_() == 'e'){
                    if(checkOtherAxes(*BBm, *BBn))
                        active_collisions_.insert(BBm->ID(), BBn->ID(), BBPair(BBm,BBn));
                }

                BoundingBoxProjectionType* temp = vector[j];
                vector[j] = vector[j-1];
                vector[j-1] = temp;
                --j;
                BBn->increase_index();
                BBm->decrease_index();
            }

        }
    }


}

#endif //DEMSIM_COLLISION_DETECTOR_H
