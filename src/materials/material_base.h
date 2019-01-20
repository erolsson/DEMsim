//
// Created by erolsson on 2018-07-30.
//

#ifndef DEMSIM_MATERIAL_BASE_H
#define DEMSIM_MATERIAL_BASE_H
namespace DEM {
    class MaterialBase {
    public:
        MaterialBase(unsigned id_number, double dens) : id(id_number), density(dens) {}
        virtual ~MaterialBase() = default;
        unsigned id ;
        double density;
    };
}

#endif //DEMSIM_MATERIAL_BASE_H
