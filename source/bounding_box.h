//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_BOUNDING_BOX_H
#define DEMSIM_BOUNDING_BOX_H

#include <algorithm>

#include "bounding_box_projection.h"
#include "cylinder.h"
#include "point_surface.h"
#include "surface.h"

namespace DEM {
    template <typename ForceModel, typename ParticleType>
    class BoundingBox {
    using BProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;
    using PointSurfaceType = PointSurface<ForceModel, ParticleType>;
    using SurfaceType = Surface<ForceModel, ParticleType>;
    using CylinderType = Cylinder<ForceModel, ParticleType>;
    public:
        explicit BoundingBox(ParticleType* particle, std::size_t index);
        explicit BoundingBox(SurfaceType* surface,  std::size_t index);
        explicit BoundingBox(CylinderType* surface,  std::size_t index);

        // Copy constructor and assignment operator needed due to pointers between different projection vectors
        // which becomes invalid when different boundingboxes are re-allocated due to vector-over-capacity
        // These constructors repairs the pointers
        BoundingBox(const BoundingBox& rhs);
        BoundingBox& operator=(const BoundingBox& rhs);


        std::size_t get_id() const;
        void update();
        void set_stretch(double stretch) {stretch_ = stretch; }

        ParticleType* get_particle() const { return particle_;}
        SurfaceType* get_surface() const { return surface_;}

        BProjectionType bx;
        BProjectionType ex;
        BProjectionType by;
        BProjectionType ey;
        BProjectionType bz;
        BProjectionType ez;

    private:
        ParticleType* particle_;
        SurfaceType* surface_;
        double stretch_{ 0. };

        // Function pointer to the update function, set in the construction of the bounding box
        // In doing so, if statements is avoided at each update
        void (BoundingBox<ForceModel, ParticleType>::*update_function)();
        void particle_update();
        void surface_update();
    };

    template<typename ForceModel, typename ParticleType>
    BoundingBox<ForceModel, ParticleType>::BoundingBox(ParticleType* particle, std::size_t index) :
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
        bx.setup();
        by.setup();
        bz.setup();
        ex.setup();
        ey.setup();
        ez.setup();
    }

    template<typename ForceModel, typename ParticleType>
    BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox::SurfaceType* surface, std::size_t index) :
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
        bx.setup();
        by.setup();
        bz.setup();
        ex.setup();
        ey.setup();
        ez.setup();
    }

    template<typename ForceModel, typename ParticleType>
    BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox::CylinderType* surface, std::size_t index) :
        particle_(nullptr),
        surface_(surface),
        update_function(&BoundingBox<ForceModel, ParticleType>::surface_update)    {

        std::cout << "Cylinder bounding box constructed " << std::endl;
    }

    //==================================================================================================================
    //                                  Copy constructor and assignment operator
    //==================================================================================================================

    template<typename ForceModel, typename ParticleType>
    BoundingBox<ForceModel, ParticleType>::BoundingBox(const BoundingBox& rhs) :
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
        bx.setup();
        by.setup();
        bz.setup();
        ex.setup();
        ey.setup();
        ez.setup();
    }

    template<typename ForceModel, typename ParticleType>
    BoundingBox<ForceModel, ParticleType>& BoundingBox<ForceModel, ParticleType>::operator=(const BoundingBox& rhs)
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
            bx.setup();
            by.setup();
            bz.setup();
            ex.setup();
            ey.setup();
            ez.setup();
        }
        return *this;
    }

    //==================================================================================================================
    //                                        Public member functions
    //==================================================================================================================


    template<typename ForceModel, typename ParticleType>
    std::size_t BoundingBox<ForceModel, ParticleType>::get_id() const
    {
        if (surface_ == nullptr) {
            return particle_->get_id();
        }
        return surface_->get_id();
    }

    template<typename ForceModel, typename ParticleType>
    void BoundingBox<ForceModel, ParticleType>::update()
    {
        (this->*update_function)();
    }

    //==================================================================================================================
    //                                        Private member functions
    //==================================================================================================================


    template<typename ForceModel, typename ParticleType>
    void BoundingBox<ForceModel, ParticleType>::particle_update()
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
    void BoundingBox<ForceModel, ParticleType>::surface_update()
    {
        auto bbox = surface_->bounding_box_values();
        bx.value = bbox[0] - stretch_;
        ex.value = bbox[1] +  stretch_;

        by.value = bbox[2] - stretch_;
        ey.value = bbox[3] +  stretch_;

        bz.value = bbox[4] - stretch_;
        ez.value = bbox[5] +  stretch_;
    }


}

#endif //DEMSIM_BOUNDING_BOX_H
