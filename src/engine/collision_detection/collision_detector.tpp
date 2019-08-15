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
                                                                    const std::vector<SurfaceType*>& surfaces,
                                                                    const ContactMatrix<ContactType>& contacts) :
        particles_(particles),
        surfaces_(surfaces),
        contacts_(contacts)
{
    //Empty constructor
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::setup()
{
    bounding_boxes_.clear();
    xproj_.clear();
    yproj_.clear();
    zproj_.clear();
    std::size_t counter = 0;
    // bounding_boxes_.reserve(particles_.size() + surfaces_.size());
    for(const auto& p: particles_){
        bounding_boxes_.emplace_back(p, counter);
        ++counter;
    }

    for(const auto& s: surfaces_){
        auto cylinder_ptr = dynamic_cast<CylinderType*>(s);
        if (cylinder_ptr != nullptr) {
            // Inspect the inward_cylinder bounding box to figure out if it is inward or not
            const auto& cyl_bounding_box = cylinder_ptr->get_bounding_box_values();
            if (std::abs(cyl_bounding_box[1]-cyl_bounding_box[0]) < 2*cylinder_ptr->get_radius()) {
                bounding_boxes_.emplace_back(cylinder_ptr, counter, true);
            }
            else {
                bounding_boxes_.emplace_back(s, counter);
            }

        }
        else {
            bounding_boxes_.emplace_back(s, counter);
        }
        ++counter;
    }

    for(auto& bounding_box: bounding_boxes_){
        xproj_.push_back(&bounding_box.bounding_box_projections[0]);
        xproj_.push_back(&bounding_box.bounding_box_projections[1]);

        yproj_.push_back(&bounding_box.bounding_box_projections[2]);
        yproj_.push_back(&bounding_box.bounding_box_projections[3]);

        zproj_.push_back(&bounding_box.bounding_box_projections[4]);
        zproj_.push_back(&bounding_box.bounding_box_projections[5]);
    }
    n_ = xproj_.size();
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::update_bounding_boxes()
{
    // #pragma omp parallel for
    for(std::size_t i = 0; i < bounding_boxes_.size(); ++i){
        bounding_boxes_[i].update();
    }
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::check_bounding_box_vector(
        std::vector<CollisionDetector::BoundingBoxProjectionType*>& vector, char axis)
{
    //std::cout << "\n";
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
                    auto id_pair = std::make_pair(bbm->get_id(), bbn->get_id());
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
    if (b1->get_particle() != nullptr && b2->get_particle() !=nullptr)
        contacts_to_create_.insert(std::make_pair(b1->get_id(), b2->get_id()),
                                   CollisionPair(b1->get_particle(), b2->get_particle()));
    else if (b1->get_particle() != nullptr && b2->get_surface() != nullptr)
        contacts_to_create_.insert(std::make_pair(b1->get_id(), b2->get_id()),
                                   CollisionPair(b1->get_particle(), b2->get_surface()));
    else if (b2->get_particle() != nullptr && b1->get_surface() != nullptr)
        contacts_to_create_.insert(std::make_pair(b2->get_id(), b1->get_id()),
                                   CollisionPair(b2->get_particle(), b1->get_surface()));
}

template<typename ForceModel, typename ParticleType>
void DEM::CollisionDetector<ForceModel, ParticleType>::destroy_contact_pair(const BoundingBoxProjectionType* b1,
        const BoundingBoxProjectionType* b2)
{
    if (contacts_.exist(b1->get_id(), b2->get_id())) {
       // There is actually a contact to destroy
        if (b1->get_particle() != nullptr && b2->get_particle() !=nullptr)
            contacts_to_destroy_.push_back(CollisionPair(b1->get_particle() ,b2->get_particle()));
        else if (b1->get_particle() != nullptr && b2->get_surface() != nullptr)
            contacts_to_destroy_.push_back(CollisionPair(b1->get_particle() ,b2->get_surface()));
        else if (b2->get_particle() != nullptr && b1->get_surface() != nullptr)
            contacts_to_destroy_.push_back(CollisionPair(b2->get_particle() ,b1->get_surface()));
    }
 }



