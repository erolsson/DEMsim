#include <fstream>
#include <iostream>
#include <vector>

#include <fenv.h>

#include "simulations/simulations.h"

int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);
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

    // Find the simulation program in the map of available simulations and execute it.
    std::cout << "Running simulation program " << program_name << "\n";
    DEM::valid_simulations()[program_name](settings_file_name);


    return 0;
}