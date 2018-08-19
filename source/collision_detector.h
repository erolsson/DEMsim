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
#include "contact_vector.h"
#include "point_surface.h"

namespace DEM {
    template <typename ForceModel, typename ParticleType> class Contact;
    template <typename ForceModel, typename ParticleType>
    class CollisionDetector {
        using ContactPointerType = std::shared_ptr<Contact<ForceModel, ParticleType> >;
    public:
        using BoundingBoxType = BoundingBox<ForceModel, ParticleType>;
        using BoundingBoxProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;
        using CollisionPair = std::pair<const BoundingBoxType*, const BoundingBoxType*>;
        CollisionDetector(const std::vector<ParticleType*>& particles,
                          const std::vector<PointSurface<ForceModel, ParticleType>*>& point_surfaces,
                          const ContactMatrix<ContactPointerType>& contacts);

        void setup();
        void do_check();  //Not const due to re-ordering of the proj vectors
        std::vector<CollisionPair> contacts_to_create() const { return contacts_to_create_.get_objects(); }
        std::vector<CollisionPair> contacts_to_destroy() const { return contacts_to_destroy_;}

    private:
        std::vector<BoundingBox<ForceModel, ParticleType> > bounding_boxes_{};
        std::vector<BoundingBoxProjectionType*> xproj_{};
        std::vector<BoundingBoxProjectionType*> yproj_{};
        std::vector<BoundingBoxProjectionType*> zproj_{};

        const std::vector<ParticleType*>& particles_;
        const std::vector<PointSurface<ForceModel, ParticleType>*>& point_surfaces_;
        const ContactMatrix<ContactPointerType>& contacts_;

        ContactVector<CollisionPair, std::pair<std::size_t, std::size_t>> contacts_to_create_{};
        std::vector<CollisionPair> contacts_to_destroy_ {};

        //BoolMatrix activeCollisions;
        void update_bounding_boxes();

        void check_bounding_box_vector(std::vector<BoundingBoxProjectionType*>& vector);
        bool check_other_axes(const BoundingBoxProjectionType* b1, const BoundingBoxProjectionType* b2) const;
    };

    template<typename ForceModel, typename ParticleType>
    void CollisionDetector<ForceModel, ParticleType>::do_check()
    {
        contacts_to_create_.clear();
        contacts_to_destroy_.clear();
        update_bounding_boxes();
        check_bounding_box_vector(xproj_);
        check_bounding_box_vector(yproj_);
        check_bounding_box_vector(zproj_);
    }

    template<typename ForceModel, typename ParticleType>
    CollisionDetector<ForceModel, ParticleType>::CollisionDetector(const std::vector<ParticleType*>& particles,
            const std::vector<PointSurface<ForceModel, ParticleType>*>& point_surfaces,
            const ContactMatrix<ContactPointerType>& contacts) :
            particles_(particles), point_surfaces_(point_surfaces), contacts_(contacts)
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
        for (unsigned i = 0; i != vector.size(); ++i){
             unsigned j = i;
             while (j != 0 && vector[j-1]->get_value() > vector[j]->get_value()){
                BoundingBoxProjectionType* BBm = vector[j];
                BoundingBoxProjectionType* BBn = vector[j-1];

                //depending on de beginnings and endings of the swapping
                // remove or add contact
                char c1 = BBm->get_position_char();
                char c2 = BBn->get_position_char();

                if (c1 == 'e' && c2 == 'b') {
                    auto id_pair = std::make_pair(BBm->get_id(), BBn->get_id());
                    if (!contacts_to_create_.erase(id_pair)) {
                        if (contacts_.exist(id_pair.first, id_pair.second)){
                            contacts_to_destroy_.push_back(std::make_pair(BBm->get_bounding_box(),
                                    BBn->get_bounding_box()));
                        }
                    }
                }
                else if (c1 == 'b' && c2 == 'e'){
                    if(check_other_axes(BBm, BBn)) {
                        contacts_to_create_.insert(std::make_pair(BBm->get_id(), BBn->get_id()),
                                                   std::make_pair(BBm->get_bounding_box(), BBn->get_bounding_box()));
                    }
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

    template<typename ForceModel, typename ParticleType>
    bool
    CollisionDetector<ForceModel, ParticleType>::check_other_axes(const CollisionDetector::BoundingBoxProjectionType* b1,
            const CollisionDetector::BoundingBoxProjectionType* b2) const
    {
        auto idx1 = b1->get_indices_on_other_axes();
        auto idx2 = b2->get_indices_on_other_axes();
        //checking the first of the axes
        if ( (*idx1[0] < *idx2[0] && *idx2[0] < *idx1[1]) || (*idx2[0] < *idx1[0] && *idx1[0] < *idx2[1]) ) {
            if ( (*idx1[2] < *idx2[2] && *idx2[2] < *idx1[3]) || (*idx2[2] < *idx1[2] && *idx1[2] < *idx2[3]) ) {
                return true;
            }
            return false;
        }
        return false;
    }


}

#endif //DEMSIM_COLLISION_DETECTOR_H
