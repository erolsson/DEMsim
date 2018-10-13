//
// Created by erolsson on 2018-09-02.
//

#include "bounding_box_projection.h"

template<typename ForceModel, typename ParticleType>
DEM::BoundingBoxProjection<ForceModel, ParticleType>::BoundingBoxProjection(BoundingBox<ForceModel, ParticleType>* bbox,
                                                                       std::size_t idx, char position, char axis,
                                                                       bool inward_cylinder) :
        position_char_(position), axis_(axis), bbox_(bbox), index_(idx), inward_cylinder_(inward_cylinder)
{
    // Empty constructor
}

template<typename ForceModel, typename ParticleType>
void DEM::BoundingBoxProjection<ForceModel, ParticleType>::setup()
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
