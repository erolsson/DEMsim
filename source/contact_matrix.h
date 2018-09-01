/*
  This file is part of the DEM program DEMsim. The code is free to run and
  develop after getting permission from the original developers Erik Olsson
  and Per-Lennart Larsson, KTH Solid Mechanics, Stockholm, Sweden

  erolsson@kth.se
  pelle@kth.se
 */

#ifndef DEMSIM_CONTACT_MATRIX_H
#define DEMSIM_CONTACT_MATRIX_H

#include <vector>
#include <map>


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
        explicit ContactMatrix(size_t);
        std::vector<T>& get_objects() {return data_;};
        const std::vector<T>& get_objects() const {return data_;}
        void insert(size_t idx1, size_t idx2, const T&);
        void erase(size_t idx1, size_t idx2);
        bool exist(size_t idx1, size_t idx2) const;
        PointerType get(size_t idx1, size_t idx2) const;

    private:
        std::vector<std::map<std::size_t, std::size_t> > data_indices_;
        std::vector<std::pair<std::size_t, std::size_t> > matrix_indices_;
        std::vector<T> data_;
    };
}


template<typename T>
DEM::ContactMatrix<T>::ContactMatrix(size_t N):
        data_indices_(std::vector<std::map<std::size_t, std::size_t> >(N, std::map<std::size_t, std::size_t>())),
        matrix_indices_(std::vector<std::pair<std::size_t, std::size_t> >()),
        data_(std::vector<T>())
{
    // Empty constructor
}


template<typename T>
void DEM::ContactMatrix<T>::insert(std::size_t idx1, std::size_t idx2, const T& obj)
{
    if(!exist(idx1, idx2) && !exist(idx2, idx1)){
        data_indices_[idx1].insert(std::pair<std::size_t, std::size_t>(idx2, data_.size()));
        matrix_indices_.push_back(std::pair<std::size_t, std::size_t>(idx1, idx2));
        data_.push_back(obj);
    }
}


template<typename T>
void DEM::ContactMatrix<T>::erase(size_t idx1, size_t idx2)
{
    //Finding the index
    auto pos = data_indices_[idx1].find(idx2);
    if(pos != data_indices_[idx1].end()) {
        data_indices_[idx1].erase(pos);
    }
    else {
        pos = data_indices_[idx2].find(idx1);
        if(pos!=data_indices_[idx2].end())
            data_indices_[idx2].erase(pos);
        else
            return;
    }

    // Tidying up the data vector by swapping elements and placing the removed element last and then pop it
    size_t i = pos->second;
    if(i != data_.size() - 1) {
        data_[i] = data_[data_.size() - 1];
        matrix_indices_[i] = matrix_indices_[matrix_indices_.size() -1];
        std::pair<std::size_t, std::size_t> ind = matrix_indices_[i];
        data_indices_[ind.first][ind.second] = i;
    }
    data_.pop_back();
    matrix_indices_.pop_back();
}


template<typename T>
bool DEM::ContactMatrix<T>::exist(size_t idx1, size_t idx2) const
{
     return data_indices_[idx1].find(idx2) != data_indices_[idx1].end() ||
         data_indices_[idx2].find(idx1) != data_indices_[idx2].end();
}


#endif // DEMSIM_CONTACT_MATRIX_H
