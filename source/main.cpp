#include <fstream>
#include <iostream>
#include <vector>

#include "simulations.h"

int main(int argc, char** argv)
{
    if (argc < 2 || DEM::valid_simulations().count(argv[1]) == 0) {
        std::cerr << "Please provide a valid program name as first argument" << '\n';
        return 0;
    }

    std::vector<std::string> arguments(argv+1, argv+argc);
    std::string program_name = arguments[0];

    std::string settings_file_name = arguments[1];

    if (argc < 3) {
        std::cerr << "Please provide a path to a simulation settings file as second argument" << '\n';
        return 0;
    }

    std::ifstream settings_file(arguments[1]);
    if (!settings_file.good()) {
        std::cerr << "Please provide a valid path to a simulation settings file as "
                     "second argument having extension .sim" << '\n';
    }

    if (program_name == "gyratory_compaction") {
        DEM::gyratory_compaction(settings_file_name);
    }

    return 0;
}