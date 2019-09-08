//
// Created by erolsson on 2018-09-20.
//

#ifndef DEMSIM_FILE_READING_FUNCTIONS_H
#define DEMSIM_FILE_READING_FUNCTIONS_H

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

namespace DEM {
    class SimulationParameters {
    public:
        explicit SimulationParameters(const std::string& settings_file_name);

        template<typename DataType>
        DataType get_parameter(const std::string& name) const;

        template<typename DataType>
        std::vector<DataType> get_vector(const std::string& name);
    private:
        std::string filename_;
        std::map<std::string, std::string> data_;
    };


    template<typename DataType>
    DataType SimulationParameters::get_parameter(const std::string& name) const
    {
        auto data_position = data_.find(name);
        if (data_position == data_.end()) {
            std::stringstream error_ss;
            error_ss << "Parameter " << name << " is not defined in file " << filename_;
            throw std::invalid_argument(error_ss.str());
        }
        // Read data back and forth to convert it to the requested data type
        std::stringstream ss;
        ss << data_position->second;
        DataType data;
        ss >> data;
        return data;
    }

    template<typename DataType>
    std::vector<DataType> SimulationParameters::get_vector(const std::string& name) {
        auto data_position = data_.find(name);
        if (data_position == data_.end()) {
            std::stringstream error_ss;
            error_ss << "Parameter " << name << " is not defined in file " << filename_;
            throw std::invalid_argument(error_ss.str());
        }
        std::string value;
        std::vector<DataType> data;
        std::stringstream ss(data_position->second);

        while (std::getline(ss, value, ',')){
            std::stringstream data_converter;
            data_converter << value;
            DataType val;
            data_converter >> val;
            data.push_back(val);        }
        return data;
    }

    template<typename DataType>
    std::vector<DataType> read_vector_from_file(const std::string& filename) {
        std::ifstream data_file(filename);
        std::string data_string;
        std::vector<DataType> data;

        while (getline(data_file, data_string)) {
            std::stringstream ss;
            ss << data_string;
            DataType data_point;
            ss >> data_point;
            data.push_back(data_point);
        }
        return data;
    }

}

#endif //DEMSIM_FILE_READING_FUNCTIONS_H
