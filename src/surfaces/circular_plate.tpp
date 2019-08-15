//
// Created by erolsson on 09/07/19.
//

#include "circular_plate.h"

using namespace DEM;

template<typename ForceModel, typename ParticleType>
DEM::CircularPlate<ForceModel, ParticleType>::CircularPlate(std::size_t id, double radius, const DEM::Vec3& normal,
                                                            const DEM::Vec3& mid_point) :
        Surface<ForceModel, ParticleType>::Surface(id),
        radius_(radius),
        normal_(normal),
        mid_point_(mid_point) {

}

template<typename ForceModel, typename ParticleType>
double DEM::CircularPlate<ForceModel, ParticleType>::distance_to_point(const DEM::Vec3& point) const {
    return 0;
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::CircularPlate<ForceModel, ParticleType>::vector_to_point(const DEM::Vec3& point) const {
    return DEM::Vec3();
}

template<typename ForceModel, typename ParticleType>
DEM::Vec3 DEM::CircularPlate<ForceModel, ParticleType>::displacement_this_inc(const DEM::Vec3& position) const {
    return DEM::Vec3();
}

template<typename ForceModel, typename ParticleType>
void DEM::CircularPlate<ForceModel, ParticleType>::move(const DEM::Vec3& distance, const DEM::Vec3& velocity) {
    velocity_ = velocity;
    mid_point_ += distance;
    update_bounding_box();
}

template<typename ForceModel, typename ParticleType>
void DEM::CircularPlate<ForceModel, ParticleType>::rotate(const DEM::Vec3& position, const DEM::Vec3& rotation_vector) {

}

template<typename ForceModel, typename ParticleType>
std::string DEM::CircularPlate<ForceModel, ParticleType>::get_output_string() const {
    return std::__cxx11::string();
}

template<typename ForceModel, typename ParticleType>
void CircularPlate<ForceModel, ParticleType>::update_bounding_box() {

}

#include "circular_plate.h"
