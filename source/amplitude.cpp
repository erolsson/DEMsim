//
// Created by erolsson on 2018-10-28.
//

#include "amplitude.h"

#include <cmath>
#include <iostream>

DEM::Amplitude::Amplitude(std::function<std::chrono::duration<double>()> time_function, bool global_time) :
        time_func_(std::move(time_function))
{
    using namespace std::chrono_literals;
    if (global_time) {
        t0_ = 0s;
    } else {
        t0_ = time_func_();
    }

}

double DEM::Amplitude::value() const
{
    double f = constant_;
    double t = time_func_().count() - - t0_.count();
    f += t*dfdt_;
    for (const auto& term : sine_terms_) {
        f += term.amplitude*(sin(2*3.1415*term.frequency*t-term.phase));
    }
    return f*factor_;
}

void DEM::Amplitude::add_sine_term(double amp, double frequency, double phase_shift)
{
    Sine_term term {amp ,frequency, phase_shift};
    sine_terms_.push_back(term);
}