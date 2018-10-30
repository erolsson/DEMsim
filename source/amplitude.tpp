//
// Created by erolsson on 2018-10-28.
//

#include "amplitude.h"

#include <cmath>

template<class EngineType>
DEM::Amplitude<EngineType>::Amplitude(const EngineType& engine, bool global_time=false) :
    engine_(engine),
{
    using namespace std::chrono_literals;
    if (global_time) {
        t0_ = 0s;
    } else {
        t0_ = engine_.get_time();
    }

}

template<class EngineType>
double DEM::Amplitude<EngineType>::value() const
{
    double f = constant;
    double t = engine_.get_time().count();
    f += (t - t0)*dfdt_;
    for (const auto& term : sine_terms_) {
        f += term.amplitude_*(sin(2*3.1415*term.frequency_-term.phase));
    }
    return f*factor_;
}

template<class EngineType>
void DEM::Amplitude<EngineType>::add_sine_term(double amp, double frequency, double phase_shift)
{
    Sine_term term {amp ,frequency, phase_shift};
    sine_terms_.push_back(term);
}