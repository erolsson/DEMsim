//
// Created by erolsson on 2018-09-20.
//

#ifndef DEMSIM_FILE_READING_FUNCTIONS_H
#define DEMSIM_FILE_READING_FUNCTIONS_H

#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

namespace DEM {
    class SimulationParameters {
    public:
        explicit SimulationParameters(const std::string& settings_file_name);

        template<typename DataType>
        DataType get(const std::string& name) const;
    private:
        std::string filename_;
        std::map<std::string, std::string> data_;
    };


    template<typename DataType>
    DataType SimulationParameters::get(const std::string& name) const
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
    std::vector<DataType> read_vector_from_file(const std::string& filename) {
        std::ifstream settings_file(filename);
        std::string data_string;
        std::vector<DataType> data;

        while(getline(settings_file, data_string)){
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
