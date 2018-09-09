//
// Created by erolsson on 2018-09-02.
//

#include "engine.h"

//=====================================================================================================================
//                        *** *** *** *** Constructors *** *** *** ***
//=====================================================================================================================

template<typename ForceModel, typename ParticleType>
DEM::Engine<ForceModel, ParticleType>::Engine() :
        collision_detector_(particles_, surfaces_, contacts_)
{
    // Empty constructor
}

//=====================================================================================================================
//                        *** *** *** *** Setup and run *** *** *** ***
//=====================================================================================================================

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::setup()
{
    contacts_.resize(number_of_objects_);
    collision_detector_.setup();
}


template<typename ForceModel, typename ParticleType>
template<typename Condition>
void DEM::Engine<ForceModel, ParticleType>::run(const Condition& condition)
{
    while (condition()) {
        time_ += settings_.increment;
        do_step();
    }
}

//=====================================================================================================================
//                        *** *** *** *** Object creation functions *** *** *** ***
//=====================================================================================================================

template<typename ForceModel, typename ParticleType>
template<typename MaterialType>
MaterialType* DEM::Engine<ForceModel, ParticleType>::create_material(double density)
{
    auto m = new MaterialType(materials_.size(), density);
    materials_.push_back(m);
    return m;
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::ParticlePointer
DEM::Engine<ForceModel, ParticleType>::create_particle(double radius, const Vec3& position,
                                                  const Vec3& velocity, MaterialBase* material)
{
    auto p = new ParticleType(radius, position, velocity, material, number_of_objects_);
    particles_.push_back(p);
    ++number_of_objects_;
    return p;
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::PointSurfacePointer
DEM::Engine<ForceModel, ParticleType>::create_point_surface(const std::vector<Vec3>& points, bool infinite)
{
    auto ps = new PointSurface<ForceModel, ParticleType>(number_of_objects_, points, infinite);
    surfaces_.push_back(ps);
    ++number_of_objects_;
    return ps;
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::CylinderPointer
DEM::Engine<ForceModel, ParticleType>::create_cylinder(double radius, const Vec3& axis, const Vec3& base_point,
                                                  double length, bool inward, bool infinite)
{
    auto c = new Cylinder<ForceModel, ParticleType>(number_of_objects_, radius, axis, base_point, length,
                                                    inward, infinite);
    surfaces_.push_back(c);
    ++number_of_objects_;
    return c;
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::do_step()
{
    move_particles();
    collision_detector_.do_check();
    update_contacts();
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::move_particles()
{

}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::update_contacts()
{

}

