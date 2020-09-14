//
// Created by erolsson on 2018-09-02.
//

#include "collision_detector.h"

#include <vector>

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::do_check()
{
    contacts_to_create_.clear();
    contacts_to_destroy_.clear();
    update_bounding_boxes();
    check_bounding_box_vector(xproj_, 'x');
    check_bounding_box_vector(yproj_, 'y');
    check_bounding_box_vector(zproj_, 'z');
}

template<typename ForceModel, typename ParticleType>
DEM::CollisionDetector<ForceModel, ParticleType>::CollisionDetector(const std::vector<ParticleType*>& particles,
                                                                    const std::vector<SurfaceType*>& surfaces) :
        particles_(particles),
        surfaces_(surfaces),
        current_contacts_()
{
    current_contacts_.resize(particles.size() + surfaces.size());
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::setup(double stretch)
{
    bounding_boxes_.clear();
    xproj_.clear();
    yproj_.clear();
    zproj_.clear();
    // bounding_boxes_.reserve(particles_.size() + surfaces_.size());
    bounding_box_stretch_ = stretch;
    for(const auto& p: particles_){
        bounding_boxes_.emplace_back(p, collision_id_counter_, stretch);
        ++collision_id_counter_;
    }

    for(const auto& s: surfaces_){
        auto cylinder_ptr = dynamic_cast<CylinderType*>(s);
        if (cylinder_ptr != nullptr) {
            // Inspect the inward_cylinder bounding box to figure out if it is inward or not
            const auto& cyl_bounding_box = cylinder_ptr->get_bounding_box_values();
            if (std::abs(cyl_bounding_box[1]-cyl_bounding_box[0]) < 2*cylinder_ptr->get_radius()) {
                bounding_boxes_.emplace_back(cylinder_ptr, collision_id_counter_, stretch, true);
            }
            else {
                bounding_boxes_.emplace_back(s, collision_id_counter_, stretch);
            }

        }
        else {
            bounding_boxes_.emplace_back(s, collision_id_counter_, stretch);
        }
        ++collision_id_counter_;
    }

    for(auto& bounding_box: bounding_boxes_){
        add_bounding_box_projections(bounding_box);
    }
    n_ = xproj_.size();
    current_contacts_.resize(collision_id_counter_);
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::restart(std::vector<ParameterMap>& restart_parameters) {
    for (const auto& parameter_line: restart_parameters) {
        auto data = parameter_line.get_parameter<std::string>("data");
        if (data == "active_collisions") {
            auto id1 = parameter_line.get_parameter<std::size_t>("object1");
            auto id2 = parameter_line.get_parameter<std::size_t>("object2");
            current_contacts_.create_item_inplace(id1, id2, true);
        }
    }
    for (auto& bbox: bounding_boxes_) {
        bbox.update();
    }

    auto bounding_box_compare = [](const auto& b1, const auto& b2) {return b1->get_value() < b2->get_value(); };
    std::sort(xproj_.begin(), xproj_.end(), bounding_box_compare);
    std::sort(yproj_.begin(), yproj_.end(), bounding_box_compare);
    std::sort(zproj_.begin(), zproj_.end(), bounding_box_compare);
    auto set_bbproj_indices = [this](auto& vec) mutable {
        for (std::size_t idx = 0; idx != vec.size(); ++idx) {
            vec[idx]->set_index(idx);
        }
    };
    set_bbproj_indices(xproj_);
    set_bbproj_indices(yproj_);
    set_bbproj_indices(zproj_);
}


template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::update_bounding_boxes()
{
    #pragma omp parallel for
    for(std::size_t i = 0; i < bounding_boxes_.size(); ++i){
        bounding_boxes_[i].update();
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::check_bounding_box_vector(
        std::vector<CollisionDetector::BoundingBoxProjectionType*>& vector, char axis)
{
    for (unsigned i = 0; i != n_; ++i) {
         unsigned j = i;
         while (j != 0 && vector[j-1]->get_value() > vector[j]->get_value()){

            BoundingBoxProjectionType* bbm = vector[j];
            BoundingBoxProjectionType* bbn = vector[j-1];

            // depending on de beginnings and endings of the swapping
            // remove or add contact
            char c1 = bbm->get_position_char();
            char c2 = bbn->get_position_char();

            if (c1 == 'e' && c2 == 'b') {
                if (! (bbm->inward_cylinder() || bbn->inward_cylinder()) || !cylinder_overlap(bbm, bbn)) {
                    auto id_pair = std::make_pair(bbm->get_collision_id(), bbn->get_collision_id());
                    if (!contacts_to_create_.erase(id_pair)) {
                        destroy_contact_pair(bbm, bbn);
                    }
                }
            }
            else if (c1 == 'b' && c2 == 'e') {
                if (((bbm->inward_cylinder() || bbn->inward_cylinder()) && cylinder_overlap(bbm, bbn)) ||
                check_other_axes(bbm, bbn, axis)) {
                    create_contact_pair(bbm, bbn);
                }
            }

            std::swap(vector[j], vector[j-1]);
            bbn->increase_index();
            bbm->decrease_index();
            --j;
         }
    }
}


template<typename ForceModel, typename ParticleType>
bool
DEM::CollisionDetector<ForceModel, ParticleType>::check_other_axes(
        const CollisionDetector::BoundingBoxProjectionType* b1,
        const CollisionDetector::BoundingBoxProjectionType* b2,
        char axis) const
{
    auto idx1 = b1->get_indices_on_other_axes(axis);
    auto idx2 = b2->get_indices_on_other_axes(axis);

    // checking the first of the axes
    for (unsigned i=0; i!=2; ++i) {
        if (!( (idx1[2*i]<idx2[2*i] && idx2[2*i]<idx1[2*i+1]) || (idx2[2*i]<idx1[2*i] && idx1[2*i]<idx2[2*i+1]) )) {
            return false;
        }
    }
    return true;
}

template<typename ForceModel, typename ParticleType>
bool DEM::CollisionDetector<ForceModel, ParticleType>::cylinder_overlap(
        const CollisionDetector::BoundingBoxProjectionType* b1,
        const CollisionDetector::BoundingBoxProjectionType* b2) const {
    const BoundingBoxProjectionType* cyl;
    const BoundingBoxProjectionType* other;
    if (b1->inward_cylinder()){
        cyl = b1;
        other = b2;
    }
    else {
        other = b1;
        cyl = b2;
    }
    auto cyl_bbox_proj = cyl->get_bounding_box()->bounding_box_projections;
    auto other_bbox_proj = other->get_bounding_box()->bounding_box_projections;

    for(unsigned i=0; i != 3; ++i) {
        if(!( cyl_bbox_proj[2*i].get_value() < other_bbox_proj[2*i].get_value() &&
            cyl_bbox_proj[2*i+1].get_value() > other_bbox_proj[2*i+1].get_value() )){
            return true;
        }
    }
    return false;
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::create_contact_pair(const BoundingBoxProjectionType* b1, 
        const BoundingBoxProjectionType* b2)
{
    if (b1->get_particle() != nullptr && b2->get_particle() != nullptr)
        contacts_to_create_.insert(std::make_pair(b1->get_collision_id(), b2->get_collision_id()),
                                   CollisionPair(b1->get_particle(), b2->get_particle()));

    else if (b1->get_particle() != nullptr && b2->get_surface() != nullptr)
        contacts_to_create_.insert(std::make_pair(b1->get_collision_id(), b2->get_collision_id()),
                                   CollisionPair(b1->get_particle(), b2->get_surface()));
    else if (b2->get_particle() != nullptr && b1->get_surface() != nullptr)
        contacts_to_create_.insert(std::make_pair(b1->get_collision_id(), b2->get_collision_id()),
                                   CollisionPair(b2->get_particle(), b1->get_surface()));
    current_contacts_.create_item_inplace(b1->get_collision_id(), b2->get_collision_id(), true);
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::destroy_contact_pair(const BoundingBoxProjectionType* b1,
        const BoundingBoxProjectionType* b2)
{
    if (current_contacts_.erase(b1->get_collision_id(), b2->get_collision_id())) {
        // There is actually a contact to destroy
        if (b1->get_particle() != nullptr && b2->get_particle() != nullptr)
            contacts_to_destroy_.push_back(CollisionPair(b1->get_particle(), b2->get_particle()));
        else if (b1->get_particle() != nullptr && b2->get_surface() != nullptr)
            contacts_to_destroy_.push_back(CollisionPair(b1->get_particle(), b2->get_surface()));
        else if (b2->get_particle() != nullptr && b1->get_surface() != nullptr)
            contacts_to_destroy_.push_back(CollisionPair(b2->get_particle(), b1->get_surface()));
    }
 }

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::add_particle(ParticleType* particle) {
    bounding_boxes_.emplace_back(particle, collision_id_counter_, bounding_box_stretch_);
    add_bounding_box_projections(bounding_boxes_.back());
    n_ += 2;
    ++collision_id_counter_;
    current_contacts_.resize(particle->get_collision_id()+1);
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::remove_particle(ParticleType* particle) {
    auto bbox_to_remove = std::find_if(bounding_boxes_.begin(), bounding_boxes_.end(),
                                       [particle](auto& bbox) {return bbox.get_particle() == particle; });
    auto bbroj_pred = [particle](const auto& bbproj){return bbproj->get_bounding_box()->get_particle() == particle;};
    xproj_.erase(std::remove_if(xproj_.begin(), xproj_.end(), bbroj_pred), xproj_.end());
    yproj_.erase(std::remove_if(yproj_.begin(), yproj_.end(), bbroj_pred), yproj_.end());
    zproj_.erase(std::remove_if(zproj_.begin(), zproj_.end(), bbroj_pred), zproj_.end());
    bounding_boxes_.erase(bbox_to_remove);
    n_ -= 2;
}



template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::add_bounding_box_projections(
        DEM::BoundingBox<ForceModel, ParticleType>& bounding_box) {
    xproj_.push_back(&bounding_box.bounding_box_projections[0]);
    xproj_.push_back(&bounding_box.bounding_box_projections[1]);

    yproj_.push_back(&bounding_box.bounding_box_projections[2]);
    yproj_.push_back(&bounding_box.bounding_box_projections[3]);

    zproj_.push_back(&bounding_box.bounding_box_projections[4]);
    zproj_.push_back(&bounding_box.bounding_box_projections[5]);
}

template<typename ForceModel, typename ParticleType>
std::vector<std::string> DEM::CollisionDetector<ForceModel, ParticleType>::restart_data() const {
    std::vector<std::string> restart_strings;
    restart_strings.push_back("data=stretch, stretch=" + std::to_string(bounding_box_stretch_));
    auto contact_indices = current_contacts_.get_matrix_indices();
    std::sort(contact_indices.begin(), contact_indices.end(),
              [](const auto& pair1, const auto& pair2) {
                  if (pair1.first == pair2.first) {
                      return pair1.second < pair2.second;
                  }
                  else {
                      return pair1.first < pair2.first;
                  }
    });
    for (const auto& id_pair: contact_indices) {
        restart_strings.push_back("data=active_collisions, object1=" + std::to_string(id_pair.first) + ", object2="
                                  + std::to_string(id_pair.second));
    }
    return restart_strings;
}
