//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_BOUNDING_BOX_H
#define DEMSIM_BOUNDING_BOX_H

#include <algorithm>

#include "bounding_box_projection.h"
#include "point_surface.h"

namespace DEM {
    template <typename ForceModel, typename ParticleType>
    class BoundingBox {
    using BProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;
    using SurfaceType = PointSurface<ForceModel, ParticleType>;
    public:
        explicit BoundingBox(ParticleType* particle, std::size_t index);
        explicit BoundingBox(SurfaceType* surface,  std::size_t index);

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
    std::size_t BoundingBox<ForceModel, ParticleType>::get_id() const
    {
        if (surface_ == nullptr) {
            return particle_->get_id();
        }
        return surface_->get_id();
    }


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
        bx.value = bbox.first.x - stretch_;
        ex.value = bbox.second.x +  stretch_;

        by.value = bbox.first.y - stretch_;
        ey.value = bbox.second.y +  stretch_;

        bz.value = bbox.first.z - stretch_;
        ez.value = bbox.second.z +  stretch_;
    }

    template<typename ForceModel, typename ParticleType>
    void BoundingBox<ForceModel, ParticleType>::update()
    {
        (this->*update_function)();
    }


}

#endif //DEMSIM_BOUNDING_BOX_H
