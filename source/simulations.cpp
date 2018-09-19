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
    // Add all valid simulation routines here
    return std::set<std::string>{"gyratory_compaction"};
}

DEM::SimulationParameters::SimulationParameters(const std::string& settings_file_name) :
    filename_(settings_file_name)
{
    std::ifstream settings_file(settings_file_name);
    std::string data_string;

    // Read all lines in the file
    auto line_count(1);
    while(getline(settings_file, data_string)){
        auto del_position = data_string.find('=');

        // If format is not identifier=value, throw exception
        if (del_position > data_string.size()-1) {
            throw std::invalid_argument("Settings parameters has to be on the form x=data");
        }

        // Split string in identifier, value pair
        std::string key = data_string.substr(0, del_position);
        std::string data = data_string.substr(del_position+1, data_string.size() - 1 - del_position);

        // If the same identifier occurs twice, throw exception
        if (data_.count(key) != 0) {
            std::stringstream error_ss;
            error_ss << "Parameter " << key << " on line " << line_count << " is already defined";
            throw std::invalid_argument(error_ss.str());
        }

        data_.insert(std::pair<std::string, std::string>(key, data));

        ++line_count;
    }
}