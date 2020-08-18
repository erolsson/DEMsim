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

#include "vec3.h"

namespace DEM {
    class ParameterMap {
    public:
        explicit ParameterMap(const std::string& delimiter="=");
        void add_data(std::string data_string);
        void add_csv_data_string (std::string csv_string);

        template<typename DataType>
        DataType get_parameter(const std::string& name) const;

        template<typename DataType>
        std::vector<DataType> get_vector(const std::string& name) const;
        [[nodiscard]] Vec3 get_vec3(const std::string& name) const {
            return Vec3(get_parameter<double>(name + "_x"),
                        get_parameter<double>(name + "_y"),
                        get_parameter<double>(name + "_z"));
        }

        [[nodiscard]] bool exist(const std::string& name) const {
            return data_.find(name) != data_.end();
        }

    private:
        std::string delimiter_;
        std::map<std::string, std::string> data_ {};
    };

    class SimulationParameters {
    public:
        explicit SimulationParameters(const std::string& settings_file_name);

        template<typename DataType>
        inline DataType get_parameter(const std::string& name) const {
            return data_.get_parameter<DataType>(name);
        }

        template<typename DataType>
        inline std::vector<DataType> get_vector(const std::string& name) const {
            return data_.get_vector<DataType>(name);
        }
    private:
        std::string filename_;
        ParameterMap data_;
    };

    template<typename DataType>
    DataType ParameterMap::get_parameter(const std::string& name) const
    {
        auto data_position = data_.find(name);
        if (data_position == data_.end()) {
            std::ostringstream error_ss;
            error_ss << "Parameter " << name << " is not defined";
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
    std::vector<DataType> ParameterMap::get_vector(const std::string& name) const {
        auto data_position = data_.find(name);
        if (data_position == data_.end()) {
            std::stringstream error_ss;
            error_ss << "Parameter " << name << " is not defined";
            throw std::invalid_argument(error_ss.str());
        }
        std::string value;
        std::vector<DataType> data;
        std::istringstream ss(data_position->second);

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
