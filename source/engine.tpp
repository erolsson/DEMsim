//
// Created by erolsson on 2018-09-02.
//

#include "engine.h"

#include <algorithm>
#include <chrono>
#include <memory>
#include <vector>

#include "contact_matrix.h"
#include "collision_detector.h"
#include "output.h"

//=====================================================================================================================
//                        *** *** *** *** Constructors *** *** *** ***
//=====================================================================================================================

template<typename ForceModel, typename ParticleType>
DEM::Engine<ForceModel, ParticleType>::Engine(std::chrono::duration<double> dt) :
        collision_detector_(particles_, surfaces_, contacts_), increment_{dt}
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
    using namespace std::chrono_literals;
    std::chrono::duration<double> logging_interval = 0.01s;
    std::chrono::duration<double> time_to_log = 0.01s;
    while (condition()) {
        time_ += increment_;
        do_step();
        time_to_log -= increment_;
        if (time_to_log <= increment_){
            time_to_log = logging_interval;
            std::cout << "Simulation time is " << get_time().count() << "\n";
        }
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
typename DEM::Engine<ForceModel, ParticleType>::OutputPointerType
DEM::Engine<ForceModel, ParticleType>::create_output(std::string directory, std::chrono::duration<double> interval)
{
    using OutputType = DEM::Output<ForceModel, ParticleType>;
    auto output_ptr = std::make_shared<OutputType>(directory, interval, *this);
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
typename DEM::Engine<ForceModel, ParticleType>::AmplitudePtrType
DEM::Engine<ForceModel, ParticleType>::set_force_control_on_surface(DEM::Surface<ForceModel, ParticleType>* surface,
                                                                    char direction)
{
    if (direction == 'x' || direction == 'y' || direction == 'z') {
        auto amp = std::make_shared<Amplitude<Engine>>
        (*this);
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


//=====================================================================================================================
//                        *** *** *** *** Private functions *** *** *** ***
//=====================================================================================================================


template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::do_step()
{
    sum_contact_forces();
    move_particles();
    collision_detector_.do_check();
    destroy_contacts();
    create_contacts();
    update_contacts();
    run_output();
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::create_contacts()
{
    const auto& contacts_to_create = collision_detector_.contacts_to_create();
    for (const auto& c_data : contacts_to_create) {
        typename ContactMatrix<Contact<ForceModel, ParticleType>>::PointerType c = nullptr;
        auto id1 = c_data.get_id_pair().first;
        auto id2 = c_data.get_id_pair().second;
        auto p1 = c_data.particle1;
        auto p2 = c_data.particle2;
        auto s = c_data.surface;
        if (s == nullptr) {
            c = contacts_.create_item_inplace(id1, id2, p1, p2, increment_);
            p2->add_contact(c, id1, -1.);
        }
        else {
            c = contacts_.create_item_inplace(id1, id2, p1, s, increment_);
            s->add_contact(c, id1);
        }
        p1->add_contact(c, id2, 1.);
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::destroy_contacts()
{
    const auto& contacts_to_destroy = collision_detector_.contacts_to_destroy();
    for (const auto& c_data : contacts_to_destroy) {
        auto id1 = c_data.get_id_pair().first;
        auto id2 = c_data.get_id_pair().second;
        auto p1 = c_data.particle1;
        auto p2 = c_data.particle2;
        auto s = c_data.surface;
        p1->remove_contact(id2);
        if (s == nullptr) {
            p2->remove_contact(id1);
        }
        else {
            s->remove_contact(id1);
        }
        contacts_.erase(id1, id2);
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::sum_contact_forces()
{
    // #pragma omp parallel for
    for (unsigned i =0; i < particles_.size(); ++i) {
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
    double dt = increment_.count();
    //#pragma omp parallel for private(F, M, new_a, new_v, new_ang_a,  new_ang_v, new_disp, new_rot)
    for (unsigned i=0; i < particles_.size(); ++i) {
        ParticleType* p = particles_[i];
        F = p->get_force();
        M = p->get_torque();

        double m = p->get_mass()*mass_scale_factor_;
        double I = p->get_inertia()*mass_scale_factor_;

        new_a = F/m + gravity_;
        new_ang_a = M/I;

        new_v = p->get_velocity()+ new_a*dt;
        new_ang_v = p->get_angular_velocity() + new_ang_a*dt;

        new_disp = new_v*dt;
        new_rot = new_ang_v*dt;

        p->set_velocity(new_v);
        p->set_angular_velocity(new_ang_v);

        p->move(new_disp);
        p->rotate(new_rot);
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::Engine<ForceModel, ParticleType>::update_contacts()
{
    auto& contact_vector = contacts_.get_objects();
    //#pragma omp parallel for
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
double DEM::Engine<ForceModel, ParticleType>::get_kinetic_energy() const
{
    double energy = 0.;
    for(const auto& p: particles_){
        energy += p->kinetic_energy();
    }
    return energy;
}
