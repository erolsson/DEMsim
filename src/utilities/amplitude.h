//
// Created by erolsson on 2018-10-27.
//

#ifndef DEMSIM_AMPLITUDE_H
#define DEMSIM_AMPLITUDE_H

#include <chrono>
#include <memory>
#include <utility>
#include <vector>
#include <functional>
#include "vec3.h"

namespace DEM {
    class Amplitude {
    public:
        explicit Amplitude(std::function<double()> time_function, double factor=1.) :
        func_(std::move(time_function)), factor_(factor) {}
        [[nodiscard]] double value() const { return  func_()*factor_; }
    private:
        std::function<double()> func_;
        double factor_;
    };
}

#endif //DEMSIM_AMPLITUDE_H
