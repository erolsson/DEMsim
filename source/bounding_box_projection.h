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
        using SurfaceType = Surface<ForceModel, ParticleType>;

    public:
        BoundingBoxProjection(BoundingBox<ForceModel, ParticleType>*, std::size_t, char);
        double value = 0;
        void setup();

    private:
        char position_char_;
        BoundingBox<ForceModel, ParticleType>* bbox_;
        std::array<std::size_t*, 4> other_indices_ = { nullptr, nullptr, nullptr, nullptr };
        std::size_t index_;
    };

    template<typename ForceModel, typename ParticleType>
    BoundingBoxProjection<ForceModel, ParticleType>::BoundingBoxProjection(BoundingBox<ForceModel, ParticleType>* bbox,
                                                                           std::size_t idx, char position) :
            bbox_(bbox), index_(idx), position_char_(position)
    {
        // Empty constructor
    }

    template<typename ForceModel, typename ParticleType>
    void BoundingBoxProjection<ForceModel, ParticleType>::setup()
    {
        if (position_char_ == 'x') {
            other_indices_[0] = &bbox_->by.index_;
            other_indices_[1] = &bbox_->ey.index_;
            other_indices_[2] = &bbox_->bz.index_;
            other_indices_[3] = &bbox_->ez.index_;
        }
        else if (position_char_ == 'y') {
            other_indices_[0] = &bbox_->bx.index_;
            other_indices_[1] = &bbox_->ex.index_;
            other_indices_[2] = &bbox_->bz.index_;
            other_indices_[3] = &bbox_->ez.index_;
        }
        else if (position_char_ == 'z') {
            other_indices_[0] = &bbox_->bx.index_;
            other_indices_[1] = &bbox_->ex.index_;
            other_indices_[2] = &bbox_->by.index_;
            other_indices_[3] = &bbox_->ey.index_;
        }
    }

    template<typename ForceModel, typename ParticleType>
    bool operator<(const BoundingBoxProjection<ForceModel, ParticleType>& b1,
                   const BoundingBoxProjection<ForceModel, ParticleType>& b2)
    {
        return b1.value < b2.value;
    }
}
#endif //DEMSIM_BOUNDING_BOX_PROJECTION_H

