//
// Created by erolsson on 2018-09-18.
//

#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <sstream>

namespace DEM {
    void gyratory_compaction(const std::string& settings_file_name);
    std::set<std::string> valid_simulations();

    class SimulationParameters {
    public:
        explicit SimulationParameters(const std::string& settings_file_name);

        template<typename DataType>
        DataType get(const std::string& name) const;
    private:
        std::map<std::string, std::string> data_;
    };


    template<typename DataType>
    DataType SimulationParameters::get(const std::string& name) const
    {
        auto data_position = data_.find(name);
        std::stringstream ss;
        ss << (*data_position).second;
        DataType data;
        ss >> data;
        return data;
    }

}

#endif //DEMSIM_SIMULATIONS_H
