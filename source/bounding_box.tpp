//
// Created by erolsson on 2018-09-02.
//

#include "bounding_box.h"

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(ParticleType* particle, std::size_t index) :
        bx(this, 2*index,   'b', 'x'),
        ex(this, 2*index+1, 'e', 'x'),
        by(this, 2*index,   'b', 'y'),
        ey(this, 2*index+1, 'e', 'y'),
        bz(this, 2*index,   'b', 'z'),
        ez(this, 2*index+1, 'e', 'z'),
        particle_(particle),
        surface_(nullptr),
        update_function(&BoundingBox<ForceModel, ParticleType>::particle_update)
{
    setup_projections();
}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox::SurfaceType* surface, std::size_t index) :
        bx(this, 2*index,   'b', 'x'),
        ex(this, 2*index+1, 'e', 'x'),
        by(this, 2*index,   'b', 'y'),
        ey(this, 2*index+1, 'e', 'y'),
        bz(this, 2*index,   'b', 'z'),
        ez(this, 2*index+1, 'e', 'z'),
        particle_(nullptr),
        surface_(surface),
        update_function(&BoundingBox<ForceModel, ParticleType>::surface_update)
{
    setup_projections();
}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox::CylinderType* cylinder, std::size_t index, bool) :
        bx(this, 2*index,   'e', 'x'),  //Note the reversed ordering of b e for x and y axis, for inward cylinders
        ex(this, 2*index+1, 'b', 'x'),
        by(this, 2*index,   'e', 'y'),
        ey(this, 2*index+1, 'b', 'y'),
        bz(this, 2*index,   'b', 'z'),
        ez(this, 2*index+1, 'e', 'z'),
        particle_(nullptr),
        surface_(cylinder),
        update_function(&BoundingBox<ForceModel, ParticleType>::surface_update)
{
    setup_projections();
}

//==================================================================================================================
//                                  Copy constructor and assignment operator
//==================================================================================================================

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>::BoundingBox(const BoundingBox& rhs) :
        bx(this, rhs.bx.get_index(), 'b', 'x'),
        ex(this, rhs.ex.get_index(), 'e', 'x'),
        by(this, rhs.by.get_index(), 'b', 'y'),
        ey(this, rhs.ey.get_index(), 'e', 'y'),
        bz(this, rhs.bz.get_index(), 'b', 'z'),
        ez(this, rhs.ez.get_index(), 'e', 'z'),
        particle_(rhs.particle_),
        surface_(rhs.surface_),
        update_function(rhs.update_function)
{
    setup_projections();
}

template<typename ForceModel, typename ParticleType>
DEM::BoundingBox<ForceModel, ParticleType>&
        DEM::BoundingBox<ForceModel, ParticleType>::operator=(const BoundingBox& rhs)
{
    if (*this != rhs) {  // Avoiding x=x
        bx(this, rhs.bx.get_index(), 'b', 'x');
        ex(this, rhs.ex.get_index(), 'e', 'x'),
        by(this, rhs.by.get_index(), 'b', 'y'),
        ey(this, rhs.ey.get_index(), 'e', 'y'),
        bz(this, rhs.bz.get_index(), 'b', 'z'),
        ez(this, rhs.ez.get_index(), 'e', 'z'),
        particle_(rhs.particle_);
        surface_(rhs.surface_);
        update_function(rhs.update_function);
        setup_projections();
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
    bx.value = particle_position.x - R - stretch_;
    ex.value = particle_position.x + R + stretch_;

    by.value = particle_position.y - R - stretch_;
    ey.value = particle_position.y + R + stretch_;

    bz.value = particle_position.z - R - stretch_;
    ez.value = particle_position.z + R + stretch_;
}

template<typename ForceModel, typename ParticleType>
void DEM::BoundingBox<ForceModel, ParticleType>::surface_update()
{
    const auto& bbox = surface_->get_bounding_box_values();
    bx.value = bbox[0] - stretch_;
    ex.value = bbox[1] +  stretch_;

    by.value = bbox[2] - stretch_;
    ey.value = bbox[3] +  stretch_;

    bz.value = bbox[4] - stretch_;
    ez.value = bbox[5] +  stretch_;
}

template<typename ForceModel, typename ParticleType>
void DEM::BoundingBox<ForceModel, ParticleType>::setup_projections()
{
    bx.setup();
    by.setup();
    bz.setup();
    ex.setup();
    ey.setup();
    ez.setup();
}
