//
// Created by erolsson on 2018-09-02.
//

#include "bounding_box_projection.h"
#include "bounding_box.h"

template<typename ForceModel, typename ParticleType>
DEM::BoundingBoxProjection<ForceModel, ParticleType>::BoundingBoxProjection(BoundingBox<ForceModel, ParticleType>* bbox,
                                                                       std::size_t idx, char position,
                                                                       bool inward_cylinder) :
        position_char_(position), bbox_(bbox), index_(idx), inward_cylinder_(inward_cylinder)
{
    // Empty constructor
}

template<typename ForceModel, typename ParticleType>
std::array<std::size_t, 4>
DEM::BoundingBoxProjection<ForceModel, ParticleType>::get_indices_on_other_axes(char axis) const {
    if (axis == 'x') {
        return std::array<std::size_t, 4>{ bbox_->bounding_box_projections[2].index_,
                                           bbox_->bounding_box_projections[3].index_,
                                           bbox_->bounding_box_projections[4].index_,
                                           bbox_->bounding_box_projections[5].index_};
    }
    else if (axis == 'y'){
        return std::array<std::size_t, 4>{ bbox_->bounding_box_projections[0].index_,
                                           bbox_->bounding_box_projections[1].index_,
                                           bbox_->bounding_box_projections[4].index_,
                                           bbox_->bounding_box_projections[5].index_};
    }

    else {
        return std::array<std::size_t, 4>{ bbox_->bounding_box_projections[0].index_,
                                           bbox_->bounding_box_projections[1].index_,
                                           bbox_->bounding_box_projections[2].index_,
                                           bbox_->bounding_box_projections[3].index_};
    }
}
