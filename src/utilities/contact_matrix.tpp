//
// Created by erolsson on 2018-09-02.
//

#include "contact_matrix.h"

#include <map>
#include <vector>

template<typename T>
DEM::ContactMatrix<T>::ContactMatrix(size_t N):
        data_indices_(std::vector<std::map<std::size_t, std::size_t> >(N, std::map<std::size_t, std::size_t>())),
        matrix_indices_(std::vector<std::pair<std::size_t, std::size_t> >()),
        data_(std::vector<T*>())
{
    // Empty constructor
}

template<typename T>
void DEM::ContactMatrix<T>::resize(size_t new_size)
{
    if (new_size > data_indices_.size()) {
        data_indices_.resize(new_size);
    }
}

template<typename T>
bool DEM::ContactMatrix<T>::erase(size_t idx1, size_t idx2)
{
    //Finding the index
    auto pos = data_indices_[idx1].find(idx2);
    size_t i = 0;
    if (pos != data_indices_[idx1].end()) {
        i = pos->second;
        data_indices_[idx1].erase(pos);
    }
    else {
        pos = data_indices_[idx2].find(idx1);
        if(pos!=data_indices_[idx2].end()) {
            i = pos->second;
            data_indices_[idx2].erase(pos);
        }
        else
            return false;
    }

    // freeing the memory allocated by the pointer at idx1 and idx2
    delete data_[i];

    // Tidying up the data vector by swapping elements and placing the removed element last and then pop it
    if (i != data_.size() - 1) {
        data_[i] = data_[data_.size() - 1];
        matrix_indices_[i] = matrix_indices_[matrix_indices_.size() -1];
        std::pair<std::size_t, std::size_t> ind = matrix_indices_[i];
        data_indices_[ind.first][ind.second] = i;
    }

    data_.pop_back();
    matrix_indices_.pop_back();
    return true;
}


template<typename T>
bool DEM::ContactMatrix<T>::exist(size_t idx1, size_t idx2) const
{
    return data_indices_[idx1].find(idx2) != data_indices_[idx1].end() ||
            data_indices_[idx2].find(idx1) != data_indices_[idx2].end();
}


template<typename T>
template<typename... Args>
typename DEM::ContactMatrix<T>::PointerType
DEM::ContactMatrix<T>::create_item_inplace(std::size_t idx1, std::size_t idx2, Args&& ... args)
{
    if (!exist(idx1, idx2) && !exist(idx2, idx1)){
        data_indices_[idx1].insert(std::pair<std::size_t, std::size_t>(idx2, data_.size()));
        matrix_indices_.push_back(std::pair<std::size_t, std::size_t>(idx1, idx2));
        T* obj = new T(std::forward<Args>(args)...);
        data_.push_back(obj);
        return obj;
    }
    else {
        return nullptr;
    }
}


template<typename T>
typename DEM::ContactMatrix<T>::PointerType DEM::ContactMatrix<T>::get(size_t idx1, size_t idx2)
{
    if (data_indices_[idx1].find(idx2) != data_indices_[idx1].end()) {
        return  data_[data_indices_[idx1][idx2]];
    }
    else if (data_indices_[idx2].find(idx1) != data_indices_[idx2].end()) {
        return  data_[data_indices_[idx2][idx1]];
    }
    else {
        return nullptr;
    }
}

