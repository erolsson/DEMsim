//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_BOUNDING_BOX_PROJECTION_H
#define DEMSIM_BOUNDING_BOX_PROJECTION_H

#include <array>

#include "surface.h"

namespace DEM {
    template<typename ForceModel, typename particleType> class BoundingBox;

    template<typename ForceModel, typename ParticleType>
    class BoundingBoxProjection {
        friend class BoundingBox<ForceModel, ParticleType>;
        using SurfaceType = Surface<ForceModel, ParticleType>;

    public:
        BoundingBoxProjection(BoundingBox<ForceModel, ParticleType>* bbox, std::size_t idx, char position,
                char axis);
        void setup();
        void increase_index() { ++index_; }
        void decrease_index() { --index_; }

        double get_value() const { return value; }
        char get_position_char() const { return position_char_; };
        std::array<const std::size_t*, 4> get_indices_on_other_axes() const { return other_indices_; }
        const BoundingBox<ForceModel, ParticleType>* get_bounding_box() const { return bbox_;}
        std::size_t get_id() const {return bbox_->get_id(); }
        std::size_t get_index() const { return  index_; }

    private:
        double value = 0;
        const char position_char_;
        const char axis_;
        const BoundingBox<ForceModel, ParticleType>* bbox_;
        std::array<const std::size_t*, 4> other_indices_ = { nullptr, nullptr, nullptr, nullptr };
        std::size_t index_;
    };

    template<typename ForceModel, typename ParticleType>
    BoundingBoxProjection<ForceModel, ParticleType>::BoundingBoxProjection(BoundingBox<ForceModel, ParticleType>* bbox,
                                                                           std::size_t idx, char position, char axis) :
            position_char_(position), axis_(axis), bbox_(bbox), index_(idx)
    {
        // Empty constructor
    }

    template<typename ForceModel, typename ParticleType>
    void BoundingBoxProjection<ForceModel, ParticleType>::setup()
    {
        if (axis_ == 'x') {
            other_indices_[0] = &(bbox_->by.index_);
            other_indices_[1] = &(bbox_->ey.index_);
            other_indices_[2] = &(bbox_->bz.index_);
            other_indices_[3] = &(bbox_->ez.index_);
        }
        else if (axis_ == 'y') {
            other_indices_[0] = &(bbox_->bx.index_);
            other_indices_[1] = &(bbox_->ex.index_);
            other_indices_[2] = &(bbox_->bz.index_);
            other_indices_[3] = &(bbox_->ez.index_);
        }
        else if (axis_ == 'z') {
            other_indices_[0] = &(bbox_->bx.index_);
            other_indices_[1] = &(bbox_->ex.index_);
            other_indices_[2] = &(bbox_->by.index_);
            other_indices_[3] = &(bbox_->ey.index_);
        }
    }
}
#endif //DEMSIM_BOUNDING_BOX_PROJECTION_H

