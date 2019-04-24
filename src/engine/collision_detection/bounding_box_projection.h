//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_BOUNDING_BOX_PROJECTION_H
#define DEMSIM_BOUNDING_BOX_PROJECTION_H

#include <array>

#include "../../surfaces/surface_base.h"

namespace DEM {
    template<typename ForceModel, typename particleType> class BoundingBox;

    template<typename ForceModel, typename ParticleType>
    class BoundingBoxProjection {
        friend class BoundingBox<ForceModel, ParticleType>;
        using SurfaceType = Surface<ForceModel, ParticleType>;

    public:
        BoundingBoxProjection(BoundingBox<ForceModel, ParticleType>* bbox, std::size_t idx, char position,
                bool inward_cylinder=false);
        void increase_index() { ++index_; }
        void decrease_index() { --index_; }

        double get_value() const { return value_; }
        char get_position_char() const { return position_char_; };
        // char get_axis() const { return axis_; }
        std::array<std::size_t, 4> get_indices_on_other_axes(char axis) const;
        const BoundingBox<ForceModel, ParticleType>* get_bounding_box() const { return bbox_;}
        std::size_t get_id() const {return bbox_->get_id(); }
        std::size_t get_index() const { return index_; };

        bool inward_cylinder() const { return  inward_cylinder_; }

        ParticleType* get_particle() const { return bbox_->get_particle(); }
        SurfaceType* get_surface() const { return bbox_->get_surface(); }

    private:
        double value_ = 0;
        const char position_char_;
        const BoundingBox<ForceModel, ParticleType>* bbox_;
        std::size_t index_;
        bool inward_cylinder_;
    };

}

#include "bounding_box_projection.tpp"
#endif //DEMSIM_BOUNDING_BOX_PROJECTION_H

