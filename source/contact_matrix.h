/*
  This file is part of the DEM program DEMsim. The code is free to run and
  develop after getting permission from the original developers Erik Olsson
  and Per-Lennart Larsson, KTH Solid Mechanics, Stockholm, Sweden

  erolsson@kth.se
  pelle@kth.se
 */

#ifndef DEMSIM_CONTACT_MATRIX_H
#define DEMSIM_CONTACT_MATRIX_H

#include <map>
#include <vector>

/*

Class for storage of contacts in a symmetric matrix with the ID's
of the objects as index

The data structure uses a vector of maps with makes insertion, and erase
in log(size(map)) which would be approx C*log (20)

When proper used, no map should have a size larger than ~20

*/

namespace DEM {
    template<typename T>
    class ContactMatrix {
    public:
        // Pointer class to be able to do in-place storage of objects instead of pointers
        using PointerType = T*;

        ContactMatrix() = default;
        explicit ContactMatrix(std::size_t);
        void resize(size_t new_size);
        std::vector<T*>& get_objects() {return data_;};
        const std::vector<T*>& get_objects() const {return data_;}
        bool erase(std::size_t idx1, std::size_t idx2);
        bool exist(std::size_t idx1, std::size_t idx2) const;
        template<typename ...Args>
        PointerType create_item_inplace(std::size_t idx1, std::size_t idx2, Args&&... args);
        PointerType get(size_t idx1, size_t idx2);

    private:
        std::vector<std::map<std::size_t, std::size_t> > data_indices_;
        std::vector<std::pair<std::size_t, std::size_t> > matrix_indices_;
        std::vector<T*> data_;
    };
}

#include "contact_matrix.tpp"

#endif // DEMSIM_CONTACT_MATRIX_H
