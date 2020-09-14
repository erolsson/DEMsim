//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_COLLISION_DETECTOR_H
#define DEMSIM_COLLISION_DETECTOR_H

#include <vector>
#include <omp.h>

#include "bounding_box.h"
#include "bounding_box_projection.h"
#include "../../utilities/contact_matrix.h"
#include "../../utilities/contact_vector.h"
#include "../../surfaces/cylinder.h"
#include "../../surfaces/point_surface.h"

namespace DEM {
    class ParameterMap;

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
                          const std::vector<SurfaceType*>& surfaces);

        void setup(double stretch);
        void restart(std::vector<ParameterMap>& restart_parameters);
        void add_particle(ParticleType* particle);
        // void add_periodic_bc(std::pair<PointSurface<ForceModel, ParticleType>,
        //         PointSurface<ForceModel, ParticleType>>& boundaries, char direction);
        void remove_particle(ParticleType* particle);
        void do_check();  //Not const due to re-ordering of the proj vectors
        std::vector<CollisionPair>& contacts_to_create() { return contacts_to_create_.get_objects(); }
        std::vector<CollisionPair>& contacts_to_destroy() { return contacts_to_destroy_;}

        std::vector<std::string> restart_data() const;

    private:
        std::vector<BoundingBox<ForceModel, ParticleType> > bounding_boxes_{};

        std::vector<BoundingBoxProjectionType*> xproj_{};
        std::vector<BoundingBoxProjectionType*> yproj_{};
        std::vector<BoundingBoxProjectionType*> zproj_{};


        std::size_t n_ = 0;
        std::size_t collision_id_counter_ = 0;
        double bounding_box_stretch_ = 1e-6;

        const std::vector<ParticleType*>& particles_;
        const std::vector<SurfaceType*>& surfaces_;
        ContactMatrix<bool> current_contacts_;

        ContactVector<CollisionPair, std::pair<std::size_t, std::size_t>> contacts_to_create_{};
        std::vector<CollisionPair> contacts_to_destroy_ {};

        //BoolMatrix activeCollisions;
        void update_bounding_boxes();

        void check_bounding_box_vector(std::vector<BoundingBoxProjectionType*>& vector, char axis);
        bool check_other_axes(const BoundingBoxProjectionType* b1, const BoundingBoxProjectionType* b2,
                char axis) const;

        bool cylinder_overlap(const BoundingBoxProjectionType* b1,
                const BoundingBoxProjectionType* b2) const;

        void create_contact_pair(const BoundingBoxProjectionType* b1, const BoundingBoxProjectionType* b2);
        void destroy_contact_pair(const BoundingBoxProjectionType* b1, const BoundingBoxProjectionType* b2);

        void add_bounding_box_projections(BoundingBox<ForceModel, ParticleType>& bounding_box);
    };

    inline std::ostream& operator<<(std::ostream& os, const std::pair<std::size_t, std::size_t>& id_pair) {
        os << id_pair.first << ", " << id_pair.second;
        return os;
    }

}

#include "collision_detector.tpp"

#endif //DEMSIM_COLLISION_DETECTOR_H
