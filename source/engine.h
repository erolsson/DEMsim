//
// Created by erolsson on 2018-07-26.
//

#ifndef DEMSIM_ENGINE_H
#define DEMSIM_ENGINE_H

#include <cstddef>
#include <iostream>
#include <map>
#include <vector>

#include "vec3.h"
#include "contact_matrix.h"
#include "material_base.h"

namespace DEM {
    template<typename ForceModel, typename ParticleType>
    class Engine {
    public:
        Engine();

        template <typename Predicate>
        void run(Predicate predicate);

    private:

    };


}

#endif //DEMSIM_ENGINE_H
