//
// Created by erolsson on 2018-09-02.
//

#include "engine.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <memory>
#include <omp.h>
#include <regex>
#include <vector>

#include "../materials/electrode_material.h"
#include "../materials/elastic_ideal_plastic_material.h"
#include "../materials/stone_material.h"
#include "../materials/linear_contact_material.h"

#include "../utilities/amplitude.h"
#include "../utilities/contact_matrix.h"
#include "../utilities/filling_functions.h"
#include "../utilities/printing_functions.h"
#include "collision_detection/collision_detector.h"
#include "output.h"
#include "periodic_bc_handler.h"
#include "../simulations/simulations.h"
#include "../utilities/file_reading_functions.h"

//=====================================================================================================================
//                        *** *** *** *** Constructors *** *** *** ***
//=====================================================================================================================

template<typename ForceModel, typename ParticleType>
DEM::Engine<ForceModel, ParticleType>::Engine(std::chrono::duration<double> dt) :
        collision_detector_(particles_, surfaces_),
        increment_{dt}
{
    // Empty constructor
}

template<typename ForceModel, typename ParticleType>
DEM::Engine<ForceModel, ParticleType>::Engine(const std::string& restart_file_name) :
        collision_detector_(particles_, surfaces_) {
    std::ifstream restart_file(restart_file_name);
    if (!restart_file.good()) {
        std::stringstream error_ss;
        error_ss << "The restart file" << restart_file_name << " is not found\n";
        throw std::invalid_argument(error_ss.str());
    }

    std::string data_string;
    std::map<std::string,      std::vector<DEM::ParameterMap>> keyword_data {
        {"*material",           std::vector<DEM::ParameterMap> {} },
        {"*output",             std::vector<DEM::ParameterMap>{} },
        {"*surface",            std::vector<DEM::ParameterMap>{} },
        {"*particle",           std::vector<DEM::ParameterMap>{} },
        {"*contact",            std::vector<DEM::ParameterMap>{} },
        {"*collision_detector", std::vector<DEM::ParameterMap>{} },
        {"*periodicbc",         std::vector<DEM::ParameterMap>{} },
    };

    DEM::ParameterMap engine_parameters;
    std::regex keyword_re ("(\\*.+?):(.+)");  // *Word: data
    std::smatch sm;
    std::size_t line_count = 0;
    while (getline(restart_file, data_string)) {
        data_string.erase(remove_if(data_string.begin(), data_string.end(), isspace), data_string.end());
        ++line_count;
        if (std::regex_match(data_string, sm, keyword_re)) {
            std::string keyword = sm[1];
            //transforming keyword to lower case
            std::transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);
            std::string data = sm[2];
            data.erase(remove_if(data.begin(), data.end(), isspace), data.end());
            if (keyword_data.find(keyword) == keyword_data.end()) {
                std::stringstream error_ss;
                error_ss << "Unsupported keyword " << keyword << " on line " << line_count << "\n";
                throw std::invalid_argument(error_ss.str());
            }
            DEM::ParameterMap par;
            par.add_csv_data_string(sm[2]);
            keyword_data[keyword].push_back(par);
        }
        else {
            // We have a parameter for engine
            std::string data;
            std::istringstream data_stream {data_string};
            while (getline(data_stream, data, ',')) {
                engine_parameters.add_data(data);
            }
        }
    }

    // Adding the parameters for engine
    increment_ = std::chrono::duration<double>(engine_parameters.get_parameter<double>("dt"));
    time_ = std::chrono::duration<double>(engine_parameters.get_parameter<double>("time"));
    mass_scale_factor_ = engine_parameters.get_parameter<double>("mass_scale_factor");
    bounding_box_stretch_ = engine_parameters.get_parameter<double>("bounding_box_stretch");
    object_id_counter_ = engine_parameters.get_parameter<double>("number_of_objects");
    collision_id_counter_= engine_parameters.get_parameter<double>("number_of_collision_objects");
    rotation_= engine_parameters.get_parameter<bool>("rotation");
    gravity_ = {engine_parameters.get_parameter<double>("gravity_x"),
                engine_parameters.get_parameter<double>("gravity_y"),
                engine_parameters.get_parameter<double>("gravity_z")};
    contacts_.resize(object_id_counter_);

    for (const auto& material_data: keyword_data["*material"]) {
        make_material_from_restart_data(material_data);
    }

    for (const auto& surface_data: keyword_data["*surface"]) {
        make_surface_from_restart_data(surface_data);
    }

    for (const auto& particle_data: keyword_data["*particle"]) {
        make_particle_from_restart_data(particle_data);
    }

    for (const auto& contact_data: keyword_data["*contact"]) {
        make_contact_from_restart_data(contact_data);
    }

    if (keyword_data["*periodicbc"].size() > 0 ){
        periodic_bc_handler_ = std::make_unique<PeriodicBCHandlerType>(*this, particles_, collision_detector_,contacts_,
                                                                       keyword_data["*periodicbc"]);
    }

    for (const auto& output_data: keyword_data["*output"]) {
        make_output_from_restart_data(output_data);
    }
    /*
    for (auto& c: contacts_.get_objects()) {
        c->update();
    }
    */
    for (auto& p: particles_) {
        p->sum_contact_forces();
    }
    collision_detector_.setup(bounding_box_stretch_);
    collision_detector_.restart(keyword_data["*collision_detector"]);
    collision_detector_.do_check();
}

//=====================================================================================================================
//                        *** *** *** *** Setup and run *** *** *** ***
//=====================================================================================================================

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::setup()
{
    auto max_radius = [](const ParticlePointer& p1, const ParticlePointer& p2) -> bool {
        return p1->get_radius() < p2->get_radius();
    };

    auto pr_max = *std::max_element(particles_.begin(), particles_.end(), max_radius);
    setup(pr_max->get_radius()*1e-2);

}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::setup(double bounding_box_stretch) {
    std::cout << "Number of objects is " << object_id_counter_ << "\n";
    bounding_box_stretch_ = bounding_box_stretch;
    contacts_.resize(object_id_counter_);
    collision_detector_.setup(bounding_box_stretch);
    if (periodic_bc_handler_ != nullptr) {
        periodic_bc_handler_->set_boundary_stretch(bounding_box_stretch);
    }
}

template<typename ForceModel, typename ParticleType>
template<typename Condition>
void DEM::Engine<ForceModel, ParticleType>::run(Condition& condition)
{
    using namespace std::chrono_literals;
    std::chrono::duration<double> logging_interval = 0.01s;
    std::chrono::duration<double> time_to_log = 0.01s;
    // Run all outputs in the beginning of the simulation
    for (auto& o : outputs_) {
        o->run_output();
    }
    while (condition()) {
        time_ += increment_;
        do_step();
        time_to_log -= increment_;
        if (time_to_log <= increment_) {
            time_to_log = logging_interval;
            std::cout << "Simulation time is " << get_time().count() << std::endl;
            auto velocity_pair = max_particle_velocity();
            std::cout << "Fastest particle is " <<  velocity_pair.first << " with a speed of "
                      <<  velocity_pair.second << std::endl;
        }
    }
    // Running all outputs in the end of the simulation
    for (auto& o : outputs_) {
        o->run_output();
    }
    std::cout << "Simulation finalized at " << get_time().count() << std::endl;
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
    auto p = new ParticleType(radius, position, velocity, material, object_id_counter_, collision_id_counter_);
    particles_.push_back(p);
    ++object_id_counter_;
    ++collision_id_counter_;
    return p;
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::PointSurfacePointer
DEM::Engine<ForceModel, ParticleType>::create_point_surface(const std::vector<Vec3>& points, bool infinite,
                                                            const std::string& name, bool adhesive)
{
    auto ps = new PointSurface<ForceModel, ParticleType>(object_id_counter_, points, infinite, name, adhesive,
                                                         collision_id_counter_);
    surfaces_.push_back(ps);
    ++object_id_counter_;
    ++collision_id_counter_;
    return ps;
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::PointSurfacePointer
DEM::Engine<ForceModel, ParticleType>::create_point_surface(const std::vector<Vec3>& points, bool infinite,
                                                            const char* name, bool adhesive)
{
    return create_point_surface(points, infinite, std::string(name), adhesive);
}



template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::PointSurfacePointer
DEM::Engine<ForceModel, ParticleType>::create_point_surface(const std::vector<Vec3>& points, bool infinite,
                                                            bool adhesive)
{
    std::stringstream name_ss;
    name_ss << "point_surface_" << object_id_counter_;
    return create_point_surface(points, infinite, name_ss.str().c_str(), adhesive);
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::DeformablePointSurfacePointer
DEM::Engine<ForceModel, ParticleType>::create_deformable_point_surface(const std::vector<Vec3>& points, bool adhesive) {
    std::stringstream name_ss;
    name_ss << "deformable_point_surface_" << object_id_counter_;
    return create_deformable_point_surface(points, name_ss.str().c_str(), adhesive);
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::DeformablePointSurfacePointer
DEM::Engine<ForceModel, ParticleType>::create_deformable_point_surface(const std::vector<Vec3>& points,
                                                                       const char* name, bool adhesive) {

    auto dps = new DeformablePointSurface<ForceModel, ParticleType>(object_id_counter_, points, true, name, adhesive,
                                                                    collision_id_counter_);
    surfaces_.push_back(dps);
    ++object_id_counter_;
    ++collision_id_counter_;
    return dps;
}


template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::CylinderPointer
DEM::Engine<ForceModel, ParticleType>::create_cylinder(double radius, const Vec3& axis, const Vec3& base_point,
                                                  double length, const std::string& name, bool inward, bool infinite,
                                                  bool closed_ends)
{
    auto c = new Cylinder<ForceModel, ParticleType>(object_id_counter_, radius, axis, base_point, length,
                                                    name, inward, infinite, closed_ends, collision_id_counter_);
    surfaces_.push_back(c);
    ++object_id_counter_;
    ++collision_id_counter_;
    return c;
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::CylinderPointer
DEM::Engine<ForceModel, ParticleType>::create_cylinder(double radius, const Vec3& axis, const Vec3& base_point,
                                                       double length, bool inward, bool infinite,
                                                       bool closed_ends)
{
    std::stringstream name_ss;
    name_ss << "cylinder_" << object_id_counter_;
    return create_cylinder(radius, axis, base_point, length, name_ss.str(), inward, infinite, closed_ends);
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::OutputPointerType
DEM::Engine<ForceModel, ParticleType>::create_output(std::string directory, std::chrono::duration<double> interval){
    std::ostringstream ss;
    ss << "output_" << outputs_.size();
    return create_output(directory, interval, ss.str());
}


template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::OutputPointerType
DEM::Engine<ForceModel, ParticleType>::create_output(std::string directory, std::chrono::duration<double> interval,
                                                     const std::string& name)
{
    using OutputType = DEM::Output<ForceModel, ParticleType>;
    auto output_ptr = std::make_shared<OutputType>(directory, interval, *this, name);
    outputs_.push_back(output_ptr);
    return output_ptr;
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::
        remove_output(const DEM::Engine<ForceModel, ParticleType>::OutputPointerType& output_to_remove)
{
    outputs_.erase(std::remove(outputs_.begin(), outputs_.end(), output_to_remove), outputs_.end());
}


template<typename ForceModel, typename ParticleType>
std::shared_ptr<DEM::Amplitude>
DEM::Engine<ForceModel, ParticleType>::set_force_control_on_surface(DEM::Surface<ForceModel, ParticleType>* surface,
                                                                    char direction, bool global_time)
{
    if (direction == 'x' || direction == 'y' || direction == 'z') {
        auto amp = std::make_shared<DEM::Amplitude>([](){ return 0.; });
        surface->set_force_amplitude(amp, direction);
        return amp;
    }
    else {
        throw std::invalid_argument("Axis must be x, y or z");
    }
}

template<typename ForceModel, typename ParticleType>
void
DEM::Engine<ForceModel, ParticleType>::remove_force_control_on_surface(DEM::Surface<ForceModel, ParticleType>* surface,
                                                                       char direction)
{
    if (direction == 'x' || direction == 'y' || direction == 'z') {
        surface->set_force_amplitude(nullptr, direction);
    }
    else {
        throw std::invalid_argument("Axis must be x, y or z");
    }
}

template<typename ForceModel, typename ParticleType>
[[maybe_unused]] std::pair<double, std::size_t>
DEM::Engine<ForceModel, ParticleType>::set_viscocity_parameters(double viscosity,
                                                                                               size_t order)
{
    std::pair<double, std::size_t> parameters {viscosity, order};
    viscocity_parameters_.push_back(parameters);
    return parameters;
}

template<typename ForceModel, typename ParticleType>
[[maybe_unused]] void
DEM::Engine<ForceModel, ParticleType>::remove_viscosity_parameters(std::pair<double, std::size_t> parameter_pair)
{
    viscocity_parameters_.erase(std::remove(viscocity_parameters_.begin(), viscocity_parameters_.end(),
                                                 parameter_pair), viscocity_parameters_.end());
}

template<typename ForceModel, typename ParticleType>
[[maybe_unused]] void DEM::Engine<ForceModel, ParticleType>::add_periodic_boundary_condition(char axis,
                                                                                             double boundary_min,
                                                                                             double boundary_max) {
    if (periodic_bc_handler_ == nullptr) {
        periodic_bc_handler_ = std::make_unique<PeriodicBCHandlerType>(*this, particles_, collision_detector_,
                                                                       contacts_);
    }
    periodic_bc_handler_->add_periodic_bc(axis, boundary_min, boundary_max);
}

template<typename ForceModel, typename ParticleType>
[[maybe_unused]] void DEM::Engine<ForceModel, ParticleType>::set_periodic_boundary_condition_strain_rate(
        char axis, double strain_rate) {
    periodic_bc_handler_->set_periodic_bc_strain_rate(axis, strain_rate);
}

template<typename ForceModel, typename ParticleType>
[[maybe_unused]] void DEM::Engine<ForceModel, ParticleType>::set_periodic_boundary_condition_velocity(char axis,
                                                                                                      double velocity) {
    periodic_bc_handler_->set_periodic_bc_velocity(axis, velocity);
}

//=====================================================================================================================
//                        *** *** *** *** Get simulation data *** *** *** ***
//=====================================================================================================================


template<typename ForceModel, typename ParticleType>
double DEM::Engine<ForceModel, ParticleType>::get_kinetic_energy() const
{
    double energy = 0.;
    for(const auto& p: particles_){
        energy += p->kinetic_energy();
    }
    return energy;
}

template<typename ForceModel, typename ParticleType>
std::pair<size_t, double> DEM::Engine<ForceModel, ParticleType>::max_particle_velocity() const
{
    auto particle_velocity = [](const ParticlePointer& p1, const ParticlePointer& p2) -> bool
    { return p1->get_velocity().length() < p2->get_velocity().length(); };

    auto p = *std::max_element(particles_.begin(), particles_.end(), particle_velocity);
    return std::make_pair(p->get_id(), p->get_velocity().length());
}

template<typename ForceModel, typename ParticleType>
std::pair<size_t, double> DEM::Engine<ForceModel, ParticleType>::max_surface_velocity() const {
    auto surface_velocity = [](const SurfaceType* s1, const SurfaceType* s2) -> bool
    { return s1->get_velocity().length() < s2->get_velocity().length(); };
    auto s = *std::max_element(surfaces_.begin(), surfaces_.end(), surface_velocity);
    return std::make_pair(s->get_id(), s->get_velocity().length());
}

template<typename ForceModel, typename ParticleType>
std::array<double, 6> DEM::Engine<ForceModel, ParticleType>::get_bounding_box() const
{
    auto x_min = [](const ParticlePointer& p1, const ParticlePointer& p2) -> bool {
        return p1->get_position().x() - p1->get_radius() < p2->get_position().x() - p2->get_radius();
    };

    auto x_max = [](const ParticlePointer& p1, const ParticlePointer& p2) -> bool {
        return p1->get_position().x() + p1->get_radius() < p2->get_position().x() + p2->get_radius();
    };

    auto y_min = [](const ParticlePointer& p1, const ParticlePointer& p2) -> bool {
        return p1->get_position().y() - p1->get_radius() < p2->get_position().y() - p2->get_radius();
    };

    auto y_max = [](const ParticlePointer& p1, const ParticlePointer& p2) -> bool {
        return p1->get_position().y() + p1->get_radius() < p2->get_position().y() + p2->get_radius();
    };

    auto z_min = [](const ParticlePointer& p1, const ParticlePointer& p2) -> bool {
        return p1->get_position().z() - p1->get_radius() < p2->get_position().z() - p2->get_radius();
    };

    auto z_max = [](const ParticlePointer& p1, const ParticlePointer& p2) -> bool {
        return p1->get_position().z() + p1->get_radius() < p2->get_position().z() + p2->get_radius();
    };

    auto px_min = *std::min_element(particles_.begin(), particles_.end(), x_min);
    auto px_max = *std::max_element(particles_.begin(), particles_.end(), x_max);
    auto py_min = *std::min_element(particles_.begin(), particles_.end(), y_min);
    auto py_max = *std::max_element(particles_.begin(), particles_.end(), y_max);
    auto pz_min = *std::min_element(particles_.begin(), particles_.end(), z_min);
    auto pz_max = *std::max_element(particles_.begin(), particles_.end(), z_max);

    return std::array<double, 6> {px_min->get_position().x() - px_min->get_radius(),
                                  px_max->get_position().x() + px_max->get_radius(),
                                  py_min->get_position().y() - py_min->get_radius(),
                                  py_max->get_position().y() + py_max->get_radius(),
                                  pz_min->get_position().z() - pz_min->get_radius(),
                                  pz_max->get_position().z() + pz_max->get_radius()};
}

template<typename ForceModel, typename ParticleType>
DEM::MaterialBase* DEM::Engine<ForceModel, ParticleType>::get_material(std::size_t id) const {
    auto it = std::lower_bound(materials_.begin(), materials_.end(), id,
                               [](const auto& m1, const auto& rhs) {return m1->id < rhs; } );
    if ((*it)->id == id) {
        return *it;
    }
    else {
        std::ostringstream error_ss;
        error_ss << "Material " << id << " does not exist";
        throw std::invalid_argument(error_ss.str());
    }
}

template<typename ForceModel, typename ParticleType> typename DEM::Engine<ForceModel, ParticleType>::SurfaceType*
DEM::Engine<ForceModel, ParticleType>::get_surface(std::size_t id) const {
    auto it = std::lower_bound(surfaces_.begin(), surfaces_.end(), id,
                               [](const auto& s1, const auto& rhs) {return s1->get_id() < rhs; } );
    if ((*it)->get_id() == id) {
        return *it;
    }
    else {
        std::ostringstream error_ss;
        error_ss << "Surface " << id << " does not exist";
        throw std::invalid_argument(error_ss.str());
    }
}

template<typename ForceModel, typename ParticleType>
template<typename SurfaceT>
SurfaceT DEM::Engine<ForceModel, ParticleType>::get_surface(const std::string& surface_name) const {
    auto it = std::find_if(surfaces_.begin(), surfaces_.end(),
                        [&surface_name](const auto& s1) { return s1->get_name() == surface_name; } );
    if (it == surfaces_.end()) {
        std::ostringstream error_ss;
        error_ss << "Surface " << surface_name << " does not exist";
        throw std::invalid_argument(error_ss.str());
    }

    SurfaceT surface = dynamic_cast<SurfaceT>(*it);
    if (surface) {
        return surface;
    }
    std::ostringstream error_ss;
    error_ss << "Surface " << surface_name << " has wrong type";
    throw std::invalid_argument(error_ss.str());
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::OutputPointerType
DEM::Engine<ForceModel, ParticleType>::get_output(const std::string& output_name) const {
    auto it = std::find_if(outputs_.begin(), outputs_.end(),
                               [&output_name](const auto& o1) {return o1->get_name() == output_name; } );
    if (it != outputs_.end()) {
        return *it;
    }
    else {
        std::ostringstream error_ss;
        error_ss << "Output " << output_name << " does not exist";
        throw std::invalid_argument(error_ss.str());
    }
}

template<typename ForceModel, typename ParticleType>
typename DEM::Engine<ForceModel, ParticleType>::ParticlePointer
DEM::Engine<ForceModel, ParticleType>::get_particle(std::size_t id) const {
    auto it = std::lower_bound(particles_.begin(), particles_.end(), id,
                               [](const auto& p1, const auto& rhs) {return p1->get_id() < rhs; } );
    if ((*it)->get_id() == id) {
        return *it;
    }
    else {
        std::ostringstream error_ss;
        error_ss << "Particle " << id << " does not exist";
        throw std::invalid_argument(error_ss.str());
    }
}


//=====================================================================================================================
//                        *** *** *** *** Private functions *** *** *** ***
//=====================================================================================================================


template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::do_step()
{
    // std::cout << "new step at time " << get_time().count() << "\n";
    move_particles();
    move_surfaces();
    collision_detector_.do_check();
    destroy_contacts();
    create_contacts();
    update_contacts();
    sum_contact_forces();
    run_output();
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::make_material_from_restart_data(const DEM::ParameterMap& parameters) {
    DEM::MaterialBase* mat;
    auto material_type = parameters.get_parameter<std::string>("type");
    if (material_type == "elastic_ideal_plastic_material") {
        mat = new DEM::ElasticIdealPlasticMaterial(parameters);
    }
    else if (material_type == "electrode_material") {
        mat = new DEM::ElectrodeMaterial(parameters);
    }
    else {
        std::ostringstream error_ss;
        error_ss << "Material " << material_type << " is currently not supported and must be added in the function "
                                                    "make_material_from_restart_data in engine.tpp";
        throw std::invalid_argument(error_ss.str());
    }
    materials_.push_back(mat);
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::make_surface_from_restart_data(const DEM::ParameterMap& parameters) {
    DEM::Surface<ForceModel, ParticleType>* surface;
    auto surface_type = parameters.get_parameter<std::string>("type");
    if (surface_type == "PointSurface") {
        surface = new DEM::PointSurface<ForceModel, ParticleType>(parameters);
    }
    else if (surface_type == "Cylinder") {
        surface = new DEM::Cylinder<ForceModel, ParticleType>(parameters);
    }
    else if (surface_type == "DeformablePointSurface") {
        surface = new DEM::DeformablePointSurface<ForceModel, ParticleType>(parameters);
    }
    else {
        std::ostringstream error_ss;
        error_ss << "Surface " << surface_type << " is currently not supported and must be added in the function "
                                                    "make_surface_from_restart_data in engine.tpp";
        throw std::invalid_argument(error_ss.str());
    }
    surfaces_.push_back(surface);
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::make_output_from_restart_data(const DEM::ParameterMap& parameters) {
    OutputPointerType output = std::make_shared<DEM::Output<ForceModel, ParticleType>>(parameters, *this);
    outputs_.push_back(output);
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::make_particle_from_restart_data(const DEM::ParameterMap& parameters) {
    DEM::MaterialBase* mat = get_material(parameters.get_parameter<std::size_t>("material"));
    auto p = new ParticleType(parameters, mat);
    particles_.push_back(p);
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::make_contact_from_restart_data(const DEM::ParameterMap& parameters) {
    auto contact_type = parameters.get_parameter<std::string>("type");
    typename ContactMatrix<Contact<ForceModel, ParticleType>>::PointerType c = nullptr;

    auto id1 = parameters.get_parameter<std::size_t>("object1");
    auto id2 = parameters.get_parameter<std::size_t>("object2");
    auto p1 = get_particle(id1);

    if (contact_type == "particle-particle") {
        auto p2 = get_particle(id2);
        c = contacts_.create_item_inplace(id1, id2, p1, p2, increment_, parameters);
        p2->add_contact(c, id1, -1.);
    }
    else if (contact_type == "particle-surface") {
        auto s = get_surface(id2);
        c = contacts_.create_item_inplace(id1, id2, p1, s, increment_, parameters);
        s->add_contact(c, id1);
    }
    else {
        std::ostringstream error_ss;
        error_ss << "Contact type " << contact_type << " is currently not supported and must be added in the function "
                                                       "make_surface_from_restart_data in engine.tpp";
        throw std::invalid_argument(error_ss.str());
    }
    p1->add_contact(c, id2, 1.);
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::create_contacts()
{
    if (periodic_bc_handler_ != nullptr) {
        periodic_bc_handler_->create_periodic_bc_contacts();
    }
    const auto& contacts_to_create = collision_detector_.contacts_to_create();
    // std::cout << contacts_to_create.size() << " contacts to create after periodic bc\n";
    for (const auto& c_data : contacts_to_create) {
        typename ContactMatrix<Contact<ForceModel, ParticleType>>::PointerType c = nullptr;
        auto id1 = c_data.get_id_pair().first;
        auto id2 = c_data.get_id_pair().second;
        auto p1 = c_data.particle1;
        auto p2 = c_data.particle2;
        auto s = c_data.surface;
        if (s == nullptr) {
            c = contacts_.create_item_inplace(id1, id2, p1, p2, increment_);
        }
        else {
            c = contacts_.create_item_inplace(id1, id2, p1, s, increment_);
        }
        if (c != nullptr) {
            p1->add_contact(c, id2, 1.);
            if (p2 != nullptr) {
                p2->add_contact(c, id1, -1);
            }
            else {
                s->add_contact(c, id1);
            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::destroy_contacts()
{
    if (periodic_bc_handler_ != nullptr) {
        periodic_bc_handler_->destroy_periodic_bc_contacts();
    }
    const auto& contacts_to_destroy = collision_detector_.contacts_to_destroy();
    for (const auto& c_data : contacts_to_destroy) {
        auto id1 = c_data.get_id_pair().first;
        auto id2 = c_data.get_id_pair().second;
        auto p1 = c_data.particle1;
        auto p2 = c_data.particle2;
        auto s = c_data.surface;
        if ( contacts_.erase(id1, id2)) {
            p1->remove_contact(id2);
            if (s == nullptr) {
                p2->remove_contact(id1);
            }
            else {
                s->remove_contact(id1);
            }
        }
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::sum_contact_forces()
{
    #pragma omp parallel for default(none)
    for (unsigned i = 0; i < particles_.size(); ++i) {
        particles_[i]->sum_contact_forces();
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::move_particles()
{
    Vec3 F;
    Vec3 M;
    Vec3 new_a;
    Vec3 new_v;
    Vec3 new_ang_a;
    Vec3 new_ang_v;
    Vec3 new_disp;
    Vec3 new_rot;
    Vec3 v;
    double dt = increment_.count();
    #pragma omp parallel for private(F, M, new_a, new_v, new_ang_a,  new_ang_v, new_disp, new_rot, v)
    for (unsigned i = 0; i < particles_.size(); ++i) {
        ParticleType* p = particles_[i];
        F = p->get_force();
        v = p->get_velocity();

        if (!viscocity_parameters_.empty() && !v.is_zero()) {
            for (const auto& v_par: viscocity_parameters_) {
                F -= pow(v.length(), v_par.second)*v.normal()*v_par.first*p->get_radius()*p->get_radius();
            }
        }

        double m = p->get_mass()*mass_scale_factor_;

        new_a = F/m + gravity_;
        new_v = v + new_a*dt;
        new_disp = new_v*dt;
        p->set_velocity(new_v);
        p->set_acceleration(new_a);
        p->move(new_disp);

        if (rotation_) {
            M = p->get_torque();
            double I = p->get_inertia()*mass_scale_factor_;
            new_ang_a = M/I;
            new_ang_v = p->get_angular_velocity() + new_ang_a*dt;
            new_rot = new_ang_v*dt;
            p->set_angular_velocity(new_ang_v);
            p->set_angular_acceleration(new_ang_a);
            p->rotate(new_rot);
        }
    }

    if (periodic_bc_handler_ != nullptr) {
        periodic_bc_handler_->fulfill_periodic_bc();
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::move_surfaces()
{
    double dt = increment_.count();
    for (auto& surface : surfaces_) {
        auto surface_forces = surface->get_applied_forces();
        Vec3 velocity = surface->get_velocity();
        Vec3 distance;
        for(unsigned axis = 0; axis != 3; ++axis) {
            auto force_amp = surface_forces[axis];
            if (force_amp != nullptr) {
                double f = force_amp->value() + surface->get_total_force()[axis];
                double a = f/(surface->get_mass()*mass_scale_factor_) + gravity_[axis];
                velocity[axis] +=  a*dt;
                distance[axis] = velocity[axis]*dt;
            }
            else {
                distance[axis] = velocity[axis]*dt;
            }
        }
        surface->move(distance, velocity);
        auto deformable_surface = dynamic_cast<DeformablePointSurface<ForceModel, ParticleType>*>(surface);
        if (deformable_surface != nullptr) {
            deformable_surface->deform(increment_);
        }
    }
}


template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::update_contacts()
{
    auto& contact_vector = contacts_.get_objects();
    #pragma omp parallel for
    for (unsigned i = 0; i < contact_vector.size(); ++i) {
        auto c = contact_vector[i];
        c->update();
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::run_output()
{
    for (auto& o : outputs_) {
        o->run_output(increment_);
    }
}

template<typename ForceModel, typename ParticleType>
[[maybe_unused]] void DEM::Engine<ForceModel, ParticleType>::set_time_increment(std::chrono::duration<double> dt) {
    increment_ = dt;
    auto& contact_vector = contacts_.get_objects();
    for (unsigned i = 0; i < contact_vector.size(); ++i) {
        auto c = contact_vector[i];
        c->set_increment(dt);
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::write_restart_file(const std::string& filename) const {
    std::ofstream restart_file;
    restart_file.open(filename);
    restart_file << "dt=" << increment_.count() << "\n";
    restart_file << "time=" << time_.count() << "\n";
    restart_file << "mass_scale_factor=" << mass_scale_factor_ << std::endl;
    restart_file << "bounding_box_stretch=" << bounding_box_stretch_ << std::endl;
    restart_file << "number_of_objects=" << object_id_counter_ << std::endl;
    restart_file << "number_of_collision_objects=" << collision_id_counter_ << std::endl;
    restart_file << "rotation=" << rotation_ << std::endl;
    restart_file << DEM::named_print(gravity_, "gravity") << "\n";

    for (const auto& m: materials_) {
        restart_file << "*material: " << m->restart_data() << "\n";
    }

    for (const auto& o: outputs_) {
        restart_file << "*output: " << o->restart_data() << "\n";
    }

    for (const auto& s: surfaces_) {
        restart_file << "*surface: " << s->restart_data() << "\n";
    }

    for (const auto& p: particles_) {
        restart_file << "*particle: " << p->restart_data() << "\n";
    }

    for (const auto& c: contacts_.get_objects_sorted()) {
        restart_file << "*contact: " << c->restart_data() << "\n";
    }

    for (const auto& collision_restart_data: collision_detector_.restart_data()) {
        restart_file << "*collision_detector: " << collision_restart_data << "\n";
    }

    if (periodic_bc_handler_ != nullptr) {
        for (const auto& periodic_restart_line: periodic_bc_handler_->restart_data()) {
            restart_file << "*periodic bc: " << periodic_restart_line << "\n";
        }
    }

    restart_file.close();
}

template<typename ForceModel, typename ParticleType>
decltype(auto) DEM::Engine<ForceModel, ParticleType>::get_periodic_boundaries() const {
    if (periodic_bc_handler_ != nullptr) {
        return periodic_bc_handler_->get_periodic_boundaries();
    }
    else {
        throw std::logic_error("No periodic boundary conditions active");
    }
}
