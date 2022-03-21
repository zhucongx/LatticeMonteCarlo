#ifndef LKMC_LKMC_API_INCLUDE_PARAMETER_H_
#define LKMC_LKMC_API_INCLUDE_PARAMETER_H_
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
namespace api {
struct Parameter {
  public:
    Parameter();
    Parameter(int argc, char *argv[]);
    explicit Parameter(const std::string &param_filename);

    // parse arguments
    void ParseArgs(int argc, char *argv[]);

    // read parameter file
    void ReadParam(const std::string &param_filename);

    // show parameters
    void PrintParameters() const;

    std::string parameters_filename;
    std::string config_filename_;
    std::string json_parameters_filename_;
    unsigned long long int log_dump_steps_;
    unsigned long long int config_dump_steps_;
    unsigned long long int maximum_number_;
    double temperature_;
    std::vector<std::string> element_symbols_;
    unsigned long long int restart_steps_;
    double restart_energy_;
    double restart_time_;
};
} // namespace api
#endif //LKMC_LKMC_API_INCLUDE_PARAMETER_H_
