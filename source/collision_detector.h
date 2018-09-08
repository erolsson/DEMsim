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
#include "cylinder.h"
#include "point_surface.h"

namespace DEM {
    template <typename ForceModel, typename ParticleType> class Contact;
    template <typename ForceModel, typename ParticleType>
    class CollisionDetector {
        using ContactType = Contact<ForceModel, ParticleType>;
    public:

        using BoundingBoxType = BoundingBox<ForceModel, ParticleType>;
        // using CollisionPair = std::pair<const BoundingBoxType*, const BoundingBoxType*>;
        using BoundingBoxProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;
        using SurfaceType = Surface<ForceModel, ParticleType>;
        using CylinderType = Cylinder<ForceModel, ParticleType>;

        class CollisionPair {
        public:
            ParticleType* particle1;
            ParticleType* particle2;
            SurfaceType* surface;

            CollisionPair(ParticleType* p1, ParticleType* p2) :
                particle1(p1), particle2(p2), surface(nullptr)
            {
                // Empty constructor
            }

            CollisionPair(ParticleType* p1, SurfaceType* surf) :
                    particle1(p1), particle2(nullptr), surface(surf)
            {
                // Empty constructor
            }

            std::pair<std::size_t, size_t> get_id_pair() const {
                if (surface == nullptr)
                    return std::make_pair(particle1->get_id(), particle2->get_id());
                else
                    return std::make_pair(particle1->get_id(), surface->get_id());
            }
        };

        CollisionDetector(const std::vector<ParticleType*>& particles,
                          const std::vector<SurfaceType*>& surfaces,
                          const ContactMatrix<ContactType>& contacts);

        void setup();
        void do_check();  //Not const due to re-ordering of the proj vectors
        std::vector<CollisionPair> contacts_to_create() const { return contacts_to_create_.get_objects(); }
        std::vector<CollisionPair> contacts_to_destroy() const { return contacts_to_destroy_;}

    private:
        std::vector<BoundingBox<ForceModel, ParticleType> > bounding_boxes_{};

        std::vector<BoundingBoxProjectionType*> xproj_{};
        std::vector<BoundingBoxProjectionType*> yproj_{};
        std::vector<BoundingBoxProjectionType*> zproj_{};

        std::size_t n_ = 0;

        const std::vector<ParticleType*>& particles_;
        const std::vector<SurfaceType*>& surfaces_;
        const ContactMatrix<ContactType>& contacts_;

        ContactVector<CollisionPair, std::pair<std::size_t, std::size_t>> contacts_to_create_{};
        std::vector<CollisionPair> contacts_to_destroy_ {};

        //BoolMatrix activeCollisions;
        void update_bounding_boxes();

        void check_bounding_box_vector(std::vector<BoundingBoxProjectionType*>& vector);
        bool check_other_axes(const BoundingBoxProjectionType* b1, const BoundingBoxProjectionType* b2) const;

    };

}

#include "collision_detector.tpp"

#endif //DEMSIM_COLLISION_DETECTOR_H
