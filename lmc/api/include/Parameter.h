#ifndef LMC_LMC_API_INCLUDE_PARAMETER_H_
#define LMC_LMC_API_INCLUDE_PARAMETER_H_
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include "VectorMatrix.hpp"
namespace api {
struct Parameter {
  public:
    Parameter(int argc, char *argv[]);
    explicit Parameter(const std::string &param_filename);

    // parse arguments
    void ParseArgs(int argc, char *argv[]);

    // read parameter file
    void ReadParam(const std::string &param_filename);

    std::string parameters_filename{};
    std::string method{};

    std::string config_filename_{};
    std::string map_filename_{};
    std::string json_coefficients_filename_{};
    std::string time_temperature_filename_{};
    std::string log_type_{};
    std::string config_type_{};
    unsigned long long int log_dump_steps_{};
    unsigned long long int config_dump_steps_{};
    unsigned long long int maximum_steps_{};
    unsigned long long int thermodynamic_averaging_steps_{};
    double temperature_{};
    double initial_temperature_{};
    double decrement_temperature_{};
    std::vector<std::string> element_set_{};
    unsigned long long int restart_steps_{};
    double restart_energy_{};
    double restart_time_{};
    bool rate_corrector_{};
    bool early_stop_{};
    bool solute_disp_{};
    Vector_t vacancy_trajectory_{};
    unsigned long long int initial_steps_{};
    unsigned long long int increment_steps_{};
    size_t smallest_cluster_criteria_{};
    size_t solvent_bond_criteria_{};
    double escape_temperature_{};

    size_t factor_{};
    std::string solvent_element_{};
    std::vector<std::string> solute_element_set_{};
    std::vector<size_t> solute_number_set_{};
    unsigned long long int early_stop_steps_{};
};
} // api
#endif //LMC_LMC_API_INCLUDE_PARAMETER_H_
