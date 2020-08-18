//
// Created by erolsson on 2018-09-02.
//

#include "particle_base.h"

#include <sstream>

#include "../utilities/file_reading_functions.h"
#include "../utilities/printing_functions.h"

template<typename ForceModel>
DEM::ParticleBase<ForceModel>::ParticleBase(double mass, const Vec3& pos, const Vec3& velocity, MaterialBase* m,
                                       unsigned id) :
        id_(id),
        mass_(mass),
        material_(m),
        position_(pos),
        velocity_(velocity)
{
    // Empty constructor
}

template<typename ForceModel>
DEM::ParticleBase<ForceModel>::ParticleBase(const DEM::ParameterMap& parameters, DEM::MaterialBase *material) :
    id_(parameters.get_parameter<std::size_t>("id")),
    mass_(parameters.get_parameter<double>("mass")),
    material_(material),
    position_(parameters.get_vec3("pos")),
    velocity_(parameters.get_vec3("vel")),
    rot_(parameters.get_vec3("rot")),
    ang_vel_(parameters.get_vec3("ang_vel")),
    displacement_this_inc_(parameters.get_vec3("disp_this_inc")),
    rot_this_inc_(parameters.get_vec3("rot_this_inc"))
{

}

template<typename ForceModel>
std::string DEM::ParticleBase<ForceModel>::restart_data() const {
    using DEM::named_print;
    std::ostringstream ss;
    ss << named_print(id_, "id") << ", "
       << named_print(mass_, "mass") << ", "
       << named_print(material_->id, "material") << ", "
       << named_print(position_, "pos") << ", "
       << named_print(velocity_, "vel") << ", "
       << named_print(rot_, "rot") << ", "
       << named_print(ang_vel_, "ang_vel") << ", "
       << named_print(displacement_this_inc_, "disp_this_inc") << ", "
       << named_print(rot_this_inc_, "rot_this_inc");
    return ss.str();
}
