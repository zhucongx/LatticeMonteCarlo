#ifndef LKMC_LKMC_API_INCLUDE_UTILITY_H_
#define LKMC_LKMC_API_INCLUDE_UTILITY_H_
#include <cstring>
#include <numeric>
#include <string>
#include <vector>
#include <sys/stat.h>
namespace api {
inline std::vector<std::string> split(const std::string &s, const char *delim) {
  std::vector<std::string> res;
  char *dup = strdup(s.c_str());
  char *token = std::strtok(dup, delim);
  while (token != nullptr) {
    res.emplace_back(token);
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = std::strtok(nullptr, delim);
  }
  free(dup);
  return res;
}
} // namespace api

#endif //LKMC_LKMC_API_INCLUDE_UTILITY_H_
