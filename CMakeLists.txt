cmake_minimum_required(VERSION 3.9)
project(DEMsim)

set(CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(DEMsim
        source/contact_matrix.h
        source/engine.h
        source/main.cpp
        source/linear_contact_material.h
        source/output.h
        source/spherical_particle.h
        source/settings.h
        source/vec3.h
        source/contact.h
        source/surface.h
        source/contact_vector.h
        source/edge.h
        source/edge.cpp source/particle_base.h source/linear_stick_slip_model.cpp source/linear_stick_slip_model.h source/material_base.h source/material_base.h source/collision_detector.h source/bounding_box.h source/bounding_box_projection.cpp source/bounding_box_projection.h)

target_link_libraries (DEMsim Eigen3::Eigen)