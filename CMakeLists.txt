cmake_minimum_required(VERSION 3.9)
project(DEMsim)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g -msse4 -Ofast -fopenmp")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(DEMsim
        source/contact_matrix.h
        source/engine.h
        source/main.cpp
        source/linear_contact_material.h
        source/output.h
        source/spherical_particle.h
        source/vec3.h
        source/contact.h
        source/surface_base.h
        source/contact_vector.h
        source/particle_base.h
        source/linear_stick_slip_model.cpp
        source/linear_stick_slip_model.h
        source/material_base.h
        source/collision_detector.h
        source/bounding_box.h
        source/bounding_box_projection.h
        source/point_surface.h
        source/cylinder.h
        source/engine.tpp
        source/collision_detector.tpp
        source/bounding_box.tpp
        source/bounding_box_projection.tpp
        source/contact.tpp
        source/surface_base.tpp
        source/contact_matrix.tpp
        source/contact_vector.tpp
        source/cylinder.tpp
        source/particle_base.tpp
        source/point_surface.tpp
        source/spherical_particle.tpp
        source/output.tpp
        source/settings.h
        source/filling_functions.h
        source/filling_functions.cpp
        source/gyratory_compaction.cpp
        source/simulations.h
        source/simulations.cpp
        source/file_reading_functions.h
        source/file_reading_functions.cpp
        source/amplitude.h
        source/amplitude.cpp
        source/closed_die_compaction.cpp
        source/storakers_mesarovic_johnson.cpp
        source/storakers_mesarovic_johnson.h
        source/elastic_ideal_plastic_material.h source/fractureable_spherical_particle.tpp source/fractureable_spherical_particle.h)

target_link_libraries (DEMsim Eigen3::Eigen stdc++fs)