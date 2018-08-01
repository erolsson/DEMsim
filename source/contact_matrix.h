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
        ContactMatrix() = default;
        explicit ContactMatrix(size_t);
        std::vector<T>& get_objects() {return data_;};
        const std::vector<T>& get_objects() const {return data_;}
        void insert(size_t, size_t, const T&);
        void erase(size_t, size_t);
        bool exist(size_t, size_t) const;
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
void DEM::ContactMatrix<T>::insert(std::size_t x, std::size_t y, const T& obj)
{
    if(!exist(x, y) && !exist(y, x)){
        data_indices_[x].insert(std::pair<std::size_t, std::size_t>(y, data_.size()));
        matrix_indices_.push_back(std::pair<std::size_t, std::size_t>(x, y));
        data_.push_back(obj);
    }
}


template<typename T>
void DEM::ContactMatrix<T>::erase(size_t x, size_t y)
{
    //Finding the index
    auto pos = data_indices_[x].find(y);
    if(pos != data_indices_[x].end()) {
        data_indices_[x].erase(pos);
    }
    else {
        pos = data_indices_[y].find(x);
        if(pos!=data_indices_[y].end())
            data_indices_[y].erase(pos);
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
bool DEM::ContactMatrix<T>::exist(size_t x, size_t y) const
{
     return data_indices_[x].find(y) != data_indices_[x].end() || data_indices_[y].find(x) != data_indices_[y].end();
}


#endif // DEMSIM_CONTACT_MATRIX_H
