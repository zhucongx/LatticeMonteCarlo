#include "Iterator.h"
namespace ansys {

Iterator::Iterator(unsigned long long int initial_number,
                         unsigned long long int increment_number,
                         unsigned long long int finial_number)
    : initial_number_(initial_number),
      increment_number_(increment_number),
      finial_number_(finial_number) {}
Iterator::~Iterator() = default;
} // namespace ansys