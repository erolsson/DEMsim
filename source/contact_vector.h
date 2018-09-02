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
        void clear();

    private:
        std::map<KeyType, std::size_t> data_indices_;
        std::vector<KeyType> vector_indices_;
        std::vector<T> data_;
    };
}

#include "contact_vector.tpp"

#endif //DEMSIM_CONTACT_VECTOR_H
