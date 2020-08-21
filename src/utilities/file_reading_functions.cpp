//
// Created by erolsson on 2018-09-20.
//

#include "file_reading_functions.h"

#include <algorithm>
#include <exception>
#include <fstream>
#include <map>
#include <sstream>
#include <string>


DEM::SimulationParameters::SimulationParameters(const std::string& settings_file_name) :
        filename_(settings_file_name), data_("=")
{
    std::ifstream settings_file(settings_file_name);
    std::string data_string;

    // Read all lines in the file
    while (getline(settings_file, data_string)) {
        data_.add_data(data_string);
    }
}

DEM::ParameterMap::ParameterMap(const std::string& delimiter) : delimiter_(delimiter) {}

void DEM::ParameterMap::add_data(std::string data_string) {
    // removing whitespace
    data_string.erase(remove_if(data_string.begin(), data_string.end(), isspace), data_string.end());
    auto del_position = data_string.find(delimiter_);
    // If format is not identifier=value, throw exception
    if (del_position > data_string.size()-1) {
        std::stringstream error_ss;
        error_ss << "string " << data_string << " is not on the form name=parameter";
        throw std::invalid_argument(error_ss.str());
    }

    // Split string in identifier, value pair
    std::string key = data_string.substr(0, del_position);
    std::string data = data_string.substr(del_position+1, data_string.size() - 1 - del_position);

    // If the same identifier occurs twice, throw exception
    if (data_.count(key) != 0) {
        std::stringstream error_ss;
        error_ss << "Parameter " << key << " is already defined";
        throw std::invalid_argument(error_ss.str());
    }

    data_.insert(std::pair<std::string, std::string>(key, data));
}


void DEM::ParameterMap::add_csv_data_string(std::string csv_string) {
    csv_string.erase(remove_if(csv_string.begin(), csv_string.end(), isspace), csv_string.end());
    std::string data;
    std::istringstream csv_stream(csv_string);
    while (getline(csv_stream, data, ',')) {
        add_data(data);
    }
}