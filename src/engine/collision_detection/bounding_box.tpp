//
// Created by erolsson on 2018-09-02.
//

#include "bounding_box.h"

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(ParticleType* particle, std::size_t index, double stretch) :
        bounding_box_projections{BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e'),
                                 BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e'),
                                 BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e')},
        particle_(particle),
        surface_(nullptr),
        stretch_(stretch),
        update_function(&BoundingBox<ForceModel, ParticleType>::particle_update)
{

}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox::SurfaceType* surface, std::size_t index,
                                                        double stretch) :
        bounding_box_projections{BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e'),
                                 BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e'),
                                 BoundingBoxProjection(this, 2*index,   'b'),
                                 BoundingBoxProjection(this, 2*index+1, 'e')},
        particle_(nullptr),
        surface_(surface),
        stretch_(stretch),
        update_function(&BoundingBox<ForceModel, ParticleType>::surface_update)
{

}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox::CylinderType* cylinder, std::size_t index,
                                                        double stretch, bool) :
        bounding_box_projections{BoundingBoxProjection(this, 2*index,   'b', true),
                                 BoundingBoxProjection(this, 2*index+1, 'e', true),
                                 BoundingBoxProjection(this, 2*index,   'b', true),
                                 BoundingBoxProjection(this, 2*index+1, 'e', true),
                                 BoundingBoxProjection(this, 2*index,   'b', true),
                                 BoundingBoxProjection(this, 2*index+1, 'e', true)},
        particle_(nullptr),
        surface_(cylinder),
        stretch_(stretch),
        update_function(&BoundingBox<ForceModel, ParticleType>::surface_update)
{

}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox&& rhs) noexcept :
    bounding_box_projections(std::move(rhs.bounding_box_projections)),
    particle_(rhs.particle_),
    surface_(rhs.surface_),
    stretch_(rhs.stretch_),
    update_function(rhs.update_function)
{
    for (auto& bbox : bounding_box_projections) {
        bbox.set_bounding_box_pointer(this);
    }
}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>&
        DEM::BoundingBox<ForceModel, ParticleType>::operator=(BoundingBox&& rhs) noexcept
{
    if (this != &rhs) {  // Avoiding x=x
        bounding_box_projections = std::move(rhs.bounding_box_projections);
        particle_ = rhs.particle_;
        surface_ = rhs.surface_;
        stretch_ = rhs.stretch_;
        update_function = rhs.update_function;

        for (auto& bbox : bounding_box_projections) {
            bbox.set_bounding_box_pointer(this);
        }
    }
    return *this;

}

//==================================================================================================================
//                                        Public member functions
//==================================================================================================================


template<typename ForceModel, typename ParticleType>
std::size_t DEM::BoundingBox<ForceModel, ParticleType>::get_collision_id() const
{
    if (particle_ != nullptr) {
        return particle_->get_collision_id();
    }
    return surface_->get_collision_id();
}

template<typename ForceModel, typename ParticleType>
std::size_t DEM::BoundingBox<ForceModel, ParticleType>::get_object_id() const {
    if (particle_ != nullptr) {
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




