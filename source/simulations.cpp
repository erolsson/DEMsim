//
// Created by erolsson on 2018-09-18.
//

#include "simulations.h"

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

std::set<std::string> DEM::valid_simulations() {
    return std::set<std::string>{"gyratory_compaction"};
}

DEM::SimulationParameters::SimulationParameters(const std::string& settings_file_name)
{
    std::ifstream settings_file(settings_file_name);
    std::string data_string;
    auto line(1);
    while(getline(settings_file, data_string)){
        auto del_position = data_string.find('=');
        if (del_position > data_string.size()-1)
            throw std::invalid_argument("Settings parameters has to be on the form x=data");
        std::string key = data_string.substr(0, del_position);
        std::string data = data_string.substr(del_position+1, data_string.size() - 1 - del_position);

        if (data_.count(key) != 0) {
            std::stringstream error_ss;
            error_ss << "Parameter " << key << " on line " << line << " is already defined";
            throw std::invalid_argument(error_ss.str());
        }

        data_.insert(std::pair<std::string, std::string>(key, data));

        ++line;
    }


}