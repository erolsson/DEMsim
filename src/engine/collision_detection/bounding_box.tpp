//
// Created by erolsson on 2018-09-02.
//

#include "bounding_box.h"

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(ParticleType* particle, std::size_t index) :
        bounding_box_projections{BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e'),
                                 BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e'),
                                 BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e')},
        particle_(particle),
        surface_(nullptr),
        update_function(&BoundingBox<ForceModel, ParticleType>::particle_update)
{

}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox::SurfaceType* surface, std::size_t index) :
        bounding_box_projections{BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e'),
                                 BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e'),
                                 BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e')},
        particle_(nullptr),
        surface_(surface),
        update_function(&BoundingBox<ForceModel, ParticleType>::surface_update)
{

}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox::CylinderType* cylinder, std::size_t index, bool) :
        bounding_box_projections{BoundingBoxProjection(this, 2*index,   'b', true),
                                 BoundingBoxProjection(this, 2*index+1, 'e', true),
                                 BoundingBoxProjection(this, 2*index,   'b', true),
                                 BoundingBoxProjection(this, 2*index+1, 'e', true),
                                 BoundingBoxProjection(this, 2*index,   'b', true),
                                 BoundingBoxProjection(this, 2*index+1, 'e', true)},
        particle_(nullptr),
        surface_(cylinder),
        update_function(&BoundingBox<ForceModel, ParticleType>::surface_update)
{

}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(const BoundingBox& rhs) :
        bounding_box_projections{BoundingBoxProjection(this, rhs.bounding_box_projections[0].get_index(), 'b',
                                                       rhs.bounding_box_projections[0].inward_cylinder()),
                                 BoundingBoxProjection(this, rhs.bounding_box_projections[1].get_index(), 'e',
                                                       rhs.bounding_box_projections[1].inward_cylinder()),
                                 BoundingBoxProjection(this, rhs.bounding_box_projections[2].get_index(), 'b',
                                                       rhs.bounding_box_projections[2].inward_cylinder()),
                                 BoundingBoxProjection(this, rhs.bounding_box_projections[3].get_index(), 'e',
                                                       rhs.bounding_box_projections[3].inward_cylinder()),
                                 BoundingBoxProjection(this, rhs.bounding_box_projections[4].get_index(), 'b',
                                                       rhs.bounding_box_projections[4].inward_cylinder()),
                                 BoundingBoxProjection(this, rhs.bounding_box_projections[5].get_index(), 'e',
                                                       rhs.bounding_box_projections[5].inward_cylinder())},
        particle_(rhs.particle_),
        surface_(rhs.surface_),
        update_function(rhs.update_function)
{

}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>&
DEM::BoundingBox<ForceModel, ParticleType>::operator=(const BoundingBox& rhs)
{
    if (*this != rhs) {  // Avoiding x=x
        for (unsigned i = 0; i != 6; ++i) {
            bounding_box_projections[i] = BoundingBoxProjection(this,
                                                                rhs.bounding_box_projections[i].get_index(),
                                                                rhs.bounding_box_projections[i].get_position_char(),
                                                                rhs.bounding_box_projections[i].inward_cylinder());
        }
    }
    return *this;
}

//==================================================================================================================
//                                        Public member functions
//==================================================================================================================


template<typename ForceModel, typename ParticleType>
std::size_t DEM::BoundingBox<ForceModel, ParticleType>::get_id() const
{
    if (surface_ == nullptr) {
        return particle_->get_id();
    }
    return surface_->get_id();
}

template<typename ForceModel, typename ParticleType>
void DEM::BoundingBox<ForceModel, ParticleType>::update()
{
    (this->*update_function)();
}

//==================================================================================================================
//                                        Private member functions
//==================================================================================================================


template<typename ForceModel, typename ParticleType>
void DEM::BoundingBox<ForceModel, ParticleType>::particle_update()
{
    Vec3 particle_position = particle_->get_position();
    double R = particle_->get_radius();
    int sign = -1;
    for(unsigned i = 0; i != bounding_box_projections.size(); ++i) {
        BProjectionType& b = bounding_box_projections[i];
        b.value_ = particle_position[i/2] + sign*(R + stretch_);
        sign *= -1;
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::BoundingBox<ForceModel, ParticleType>::surface_update()
{
    auto bbox = surface_->get_bounding_box_values();
    int sign = -1;
    for(unsigned i = 0; i != bounding_box_projections.size(); ++i) {
        BProjectionType& b = bounding_box_projections[i];
        b.value_ = bbox[i] + sign*stretch_;
        sign *= -1;
    }
}
