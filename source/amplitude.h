//
// Created by erolsson on 2018-10-27.
//

#ifndef DEMSIM_AMPLITUDE_H
#define DEMSIM_AMPLITUDE_H

#include <chrono>
#include <memory>
#include <vector>
#include <functional>

namespace DEM {
    template<class EngineType>
    class Amplitude {
    public:
        Amplitude(const EngineType& engine, bool global_time);
        double value() const;
        void constant_term(double value) { constant = value; }
        void linear_term(double dfdt) { dfdt_ = dfdt; }
        void add_sine_term(double amp, double frequency, double phase_shift = 0);
        void add_multiplying_factor(double multiplying_factor) {factor_ = multiplying_factor;};

    private:
        double constant = 0;
        double dfdt_ = 0;
        std::chrono::duration<double> t0_;
        struct Sine_term {
            double amplitude ;
            double frequency ;
            double phase_ ;
        };
        std::vector<Sine_term> sine_terms_ {};
        double factor_;


        const EngineType& engine_;
    };




}

#endif //DEMSIM_AMPLITUDE_H
