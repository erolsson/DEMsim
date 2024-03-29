cmake_minimum_required(VERSION 3.9)
project(DEMsim)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -msse4 -Ofast -g -fopenmp")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(DEMsim
               src/utilities/contact_matrix.h
               src/engine/engine.h
               src/main.cpp
               src/materials/linear_contact_material.h
               src/engine/output.h
               src/particles/spherical_particle.h
               src/utilities/vec3.h
               src/engine/contact.h
               src/surfaces/surface_base.h
               src/utilities/contact_vector.h
               src/particles/particle_base.h
               src/contact_models/linear_stick_slip_model.cpp
               src/contact_models/linear_stick_slip_model.h
               src/materials/material_base.h
               src/materials/material_base.cpp
               src/engine/collision_detection/collision_detector.h
               src/engine/collision_detection/bounding_box.h
               src/engine/collision_detection/bounding_box_projection.h
               src/surfaces/point_surface.h
               src/surfaces/cylinder.h
               src/engine/engine.tpp
               src/engine/collision_detection/collision_detector.tpp
               src/engine/collision_detection/bounding_box.tpp
               src/engine/collision_detection/bounding_box_projection.tpp
               src/engine/contact.tpp
               src/surfaces/surface_base.tpp
               src/utilities/contact_matrix.tpp
               src/utilities/contact_vector.tpp
               src/surfaces/cylinder.tpp
               src/particles/particle_base.tpp
               src/surfaces/point_surface.tpp
               src/particles/spherical_particle.tpp
               src/engine/output.tpp
               src/utilities/filling_functions.h
               src/utilities/filling_functions.cpp
               src/simulations/gyratory_compaction.cpp
               src/simulations/simulations.h
               src/simulations/simulations.cpp
               src/utilities/file_reading_functions.h
               src/utilities/file_reading_functions.cpp
               src/utilities/amplitude.h
               src/utilities/amplitude.cpp
               src/simulations/closed_die_compaction.cpp
               src/contact_models/storakers_mesarovic_johnson.cpp
               src/contact_models/storakers_mesarovic_johnson.h
               src/materials/elastic_ideal_plastic_material.h
               src/materials/elastic_ideal_plastic_material.cpp
               src/particles/fractureable_spherical_particle.tpp
               src/particles/fractureable_spherical_particle.h
               src/contact_models/stone_material_contact.cpp
               src/contact_models/stone_material_contact.h
               src/materials/stone_material.h
               src/simulations/contact_tester.cpp
               src/simulations/cyclic_triaxial.cpp
               src/particles/spherical_particle_base.tpp
               src/particles/spherical_particle_base.h
               src/simulations/proctor_test.cpp
               src/simulations/stone_compaction.cpp
               src/contact_models/viscoelastic.cpp
               src/contact_models/viscoelastic.h
               src/materials/electrode_material.h
               src/materials/electrode_material.cpp
               src/utilities/printing_functions.h
               src/simulations/electrode_compaction.cpp
               src/simulations/electrode_cylinder/electrode_cylinder_compaction.cpp
               src/simulations/periodic_bc/periodic_bc_tester.cpp
               src/engine/periodic_bc_handler.h
               src/engine/periodic_bc_handler.tpp
               src/simulations/periodic_bc/periodic_bc_simulation.cpp
               src/simulations/periodic_bc/filling_periodic_box.cpp
               src/surfaces/deformable_point_surface.h
               src/surfaces/deformable_point_surface.tpp
               src/simulations/cathode_mechanical_simulations.cpp
               src/simulations/battery/battery_rve_compaction.cpp
               src/simulations/battery/battery_rve_loading.cpp
               src/simulations/battery/restart_electrode.cpp
               src/contact_models/porous_electrode_contact.cpp
               src/contact_models/porous_electrode_contact.h
               src/simulations/porous_electrode/porous_electrode_rve.cpp
               src/simulations/electrode_mechanical_test.cpp
               src/materials/elastic_bonded_material.h
               src/materials/elastic_bonded_material.cpp
               src/contact_models/hertz_with_bonds.h
               src/contact_models/hertz_with_bonds.cpp
               src/simulations/asphalt_shear_box.cpp
               src/simulations/asphalt_shear_box_bonded.cpp
               )

target_link_libraries (DEMsim stdc++fs "${OpenMP_CXX_FLAGS}" Eigen3::Eigen)