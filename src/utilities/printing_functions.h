//
// Created by erolsson on 25/07/2020.
//

#ifndef DEMSIM_PRINTING_FUNCTIONS_H
#define DEMSIM_PRINTING_FUNCTIONS_H

#include "vec3.h"

namespace DEM {
    [[nodiscard]] inline std::string named_print(const Vec3& vec, const std::string& name)
    {
        std::stringstream ss;
        ss << name << "_x=" << vec.x() << ", " << name << "_y=" <<  vec.y() << ", " << name << "_z=" << vec.z();
        return ss.str();
    }

    template <typename T>
    [[nodiscard]] inline std::string named_print(T val, const std::string& name){
        std::stringstream ss;
        ss << name << "=" << val;
        return ss.str();
    }
}

#endif //DEMSIM_PRINTING_FUNCTIONS_H
