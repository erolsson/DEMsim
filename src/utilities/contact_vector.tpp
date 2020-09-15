//
// Created by erolsson on 2018-09-02.
//

#include "contact_vector.h"

template<typename T, typename KeyType>
void DEM::ContactVector<T, KeyType>::insert(KeyType key, const T& obj)
{
    if (!exist(key)) {
        data_indices_.insert(std::pair<KeyType, std::size_t>(key, data_.size()));
        vector_indices_.push_back(key);
        data_.push_back(obj);
    }
}

template<typename T, typename KeyType>
bool DEM::ContactVector<T, KeyType>::erase(KeyType key)
{
    if (exist(key) ){
        auto pos = data_indices_.find(key);
        std::size_t i = pos->second;   // Index where the removed element was, swapping with the last element
        if (i != data_.size() - 1) {
            vector_indices_[i] = vector_indices_[vector_indices_.size()-1];
            data_[i] = data_[data_.size()-1];
            data_indices_[vector_indices_[i]] = i;
        }
        data_indices_.erase(pos);
        vector_indices_.pop_back();
        data_.pop_back();
        return true;
    }
    return false;
}

template<typename T, typename KeyType>
bool DEM::ContactVector<T, KeyType>::exist(KeyType key) const
{
    return data_indices_.find(key) != data_indices_.end();
}

template<typename T, typename KeyType>
void DEM::ContactVector<T, KeyType>::clear()
{
    data_indices_.clear();
    vector_indices_.clear();
    data_.clear();
}


