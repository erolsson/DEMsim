//
// Created by erolsson on 2018-10-27.
//

#ifndef DEMSIM_AMPLITUDE_H
#define DEMSIM_AMPLITUDE_H

#include <chrono>
#include <memory>
#include <vector>
#include <functional>
#include "vec3.h"

namespace DEM {
    class Amplitude {
    public:
        Amplitude(std::function<std::chrono::duration<double>()> time_function, bool global_time);
        [[nodiscard]] double value() const;
        void constant_term(double value) { constant_ = value; }
        void linear_term(double dfdt) { dfdt_ = dfdt; }
        void add_sine_term(double amp, double frequency, double phase_shift = 0);
        void multiplying_factor(double multiplying_factor) {factor_ = multiplying_factor;};

    private:
        double constant_ = 0;
        double dfdt_ = 0;
        std::chrono::duration<double> t0_;
        struct Sine_term {
            double amplitude;
            double frequency;
            double phase;
        };
        std::vector<Sine_term> sine_terms_ {};
        double factor_ = 1.;
        std::function<std::chrono::duration<double>()> time_func_;
    };
}

#endif //DEMSIM_AMPLITUDE_H
