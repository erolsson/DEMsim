//
// Created by erolsson on 2018-07-28.
//

#ifndef DEMSIM_CONTACT_VECTOR_H
#define DEMSIM_CONTACT_VECTOR_H

#include <map>
#include <vector>

namespace DEM{
    template<typename T, typename KeyType=std::size_t>
    class ContactVector {
    public:
        ContactVector() = default;
        std::vector<T>& get_objects() {return data_;}
        const std::vector<T>& get_objects() const {return data_;}
        void insert(KeyType key, const T& obj);
        bool erase(KeyType key);
        bool exist(KeyType key) const;

    private:
        std::map<std::size_t, std::size_t> data_indices_;
        std::vector<std::size_t> vector_indices_;
        std::vector<T> data_;
    };


    template<typename T, typename KeyType>
    void ContactVector<T, KeyType>::insert(KeyType key, const T& obj)
    {
        if (!exist(key)) {
            data_indices_.insert(std::pair<std::size_t, std::size_t>(key, data_.size()));
            vector_indices_.push_back(key);
            data_.push_back(obj);
        }
    }

    template<typename T, typename KeyType>
    bool ContactVector<T, KeyType>::erase(KeyType key)
    {
        if (exist(key) ){
            auto pos = data_indices_.find(key);
            data_indices_.erase(pos);
            std::size_t i = pos->second;   // Index where the removed element was, swapping with the last element
            if (i != data_.size() - 1) {
                vector_indices_[i] = vector_indices_[vector_indices_.size()-1];
                data_[i] = data_[data_.size()-1];
                data_indices_[vector_indices_[i]] = i;
            }
            vector_indices_.pop_back();
            data_.pop_back();
            return true;
        }
        return false;
    }

    template<typename T, typename KeyType>
    bool ContactVector<T, KeyType>::exist(KeyType key) const
    {
        return data_indices_.find(key) != data_indices_.end();
    }
}

#endif //DEMSIM_CONTACT_VECTOR_H
